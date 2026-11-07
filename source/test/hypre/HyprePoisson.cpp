/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example Hypre Poisson solver
 *
 ************************************************************************/
#include "HyprePoisson.h"

#if defined(HAVE_HYPRE)

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

extern "C" {
void SAMRAI_F77_FUNC(setexactandrhs2d, SETEXACTANDRHS2D) (const int& ifirst0,
   const int& ilast0,
   const int& ifirst1,
   const int& ilast1,
   double* exact,
   double* rhs,
   const double* dx,
   const double* xlower);
void SAMRAI_F77_FUNC(setexactandrhs3d, SETEXACTANDRHS3D) (const int& ifirst0,
   const int& ilast0,
   const int& ifirst1,
   const int& ilast1,
   const int& ifirst2,
   const int& ilast2,
   double* exact,
   double* rhs,
   const double* dx,
   const double* xlower);
}

namespace SAMRAI {

/*
 *************************************************************************
 * Constructor creates a unique context for the object and register
 * all its internal variables with the variable database.
 *************************************************************************
 */
HyprePoisson::HyprePoisson(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<solv::CellPoissonHypreSolver>& hypre_solver,
   std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs):
   d_object_name(object_name),
   d_dim(dim),
   d_poisson_hypre(hypre_solver),
   d_bc_coefs(bc_coefs)
{

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   /*
    * Get a unique context for this object.
    */
   d_context = vdb->getContext(d_object_name + ":Context");

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    */
   std::shared_ptr<pdat::CellVariable<double> > comp_soln(
      new pdat::CellVariable<double>(
         d_dim,
         object_name + ":computed solution",
         1));
   d_comp_soln_id =
      vdb->registerVariableAndContext(
         comp_soln,
         d_context,
         hier::IntVector(d_dim, 1) /* ghost cell width is 1 for stencil widths */);
   std::shared_ptr<pdat::CellVariable<double> > exact_solution(
      new pdat::CellVariable<double>(
         d_dim,
         object_name + ":exact solution"));
   d_exact_id =
      vdb->registerVariableAndContext(
         exact_solution,
         d_context,
         hier::IntVector(d_dim, 1) /* ghost cell width is 1 in case needed */);
   std::shared_ptr<pdat::CellVariable<double> > rhs_variable(
      new pdat::CellVariable<double>(
         d_dim,
         object_name
         + ":linear system right hand side"));
   d_rhs_id =
      vdb->registerVariableAndContext(
         rhs_variable,
         d_context,
         hier::IntVector(d_dim, 0) /* ghost cell width is 0 */);
}

/*
 *************************************************************************
 * Destructor does nothing interesting
 *************************************************************************
 */
HyprePoisson::~HyprePoisson()
{
}

/*
 *************************************************************************
 * Initialize data on a level.
 *
 * Allocate the solution, exact solution and rhs memory.
 * Fill the rhs and exact solution.
 *************************************************************************
 */
void HyprePoisson::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(init_data_time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);

   std::shared_ptr<hier::PatchHierarchy> patch_hierarchy = hierarchy;
   std::shared_ptr<geom::CartesianGridGeometry> grid_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
         patch_hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geom);

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   /*
    * If required, allocate all patch data on the level.
    */
   if (allocate_data) {
      level->allocatePatchData(d_comp_soln_id);
      level->allocatePatchData(d_rhs_id);
      level->allocatePatchData(d_exact_id);
   }

   /*
    * Initialize data in all patches in the level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {

      const std::shared_ptr<hier::Patch>& patch = *pi;
      if (!patch) {
         TBOX_ERROR(d_object_name
            << ": Cannot find patch.  Null patch pointer.");
      }
      hier::Box pbox = patch->getBox();
      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));

      std::shared_ptr<pdat::CellData<double> > exact_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_exact_id)));
      std::shared_ptr<pdat::CellData<double> > rhs_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_rhs_id)));
      TBOX_ASSERT(patch_geom);
      TBOX_ASSERT(exact_data);
      TBOX_ASSERT(rhs_data);

      /*
       * Set source function and exact solution.
       */
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(setexactandrhs2d, SETEXACTANDRHS2D) (
            pbox.lower()[0],
            pbox.upper()[0],
            pbox.lower()[1],
            pbox.upper()[1],
            exact_data->getPointer(),
            rhs_data->getPointer(),
            grid_geom->getDx(),
            patch_geom->getXLower());
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(setexactandrhs3d, SETEXACTANDRHS3D) (
            pbox.lower()[0],
            pbox.upper()[0],
            pbox.lower()[1],
            pbox.upper()[1],
            pbox.lower()[2],
            pbox.upper()[2],
            exact_data->getPointer(),
            rhs_data->getPointer(),
            grid_geom->getDx(),
            patch_geom->getXLower());
      }

   }    // End patch loop.
}

/*
 *************************************************************************
 * Reset the hierarchy-dependent internal information.
 *************************************************************************
 */
void HyprePoisson::resetHierarchyConfiguration(
   const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   int coarsest_level,
   int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   d_hierarchy = new_hierarchy;
}

/*
 *************************************************************************
 * Solve the Poisson problem.
 *************************************************************************
 */
bool HyprePoisson::solvePoisson()
{

   if (!d_hierarchy) {
      TBOX_ERROR("Cannot solve using an uninitialized object.\n");
   }

   const int level_number = 0;

   /*
    * Fill in the initial guess and Dirichlet boundary condition data.
    * For this example, we want u=0 on all boundaries.
    * The easiest way to do this is to just write 0 everywhere,
    * simultaneous setting the boundary values and initial guess.
    */
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(
                                                level_number));
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;
      std::shared_ptr<pdat::CellData<double> > data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_comp_soln_id)));
      TBOX_ASSERT(data);
      data->fill(0.0);
   }
   // d_poisson_hypre->setBoundaries( "Dirichlet" );
   d_poisson_hypre->setPhysicalBcCoefObject(d_bc_coefs.get());

   /*
    * Set up HYPRE solver object.
    * The problem specification is set using the
    * CellPoissonSpecifications object then passed to the solver
    * for setting the coefficients.
    */
   d_poisson_hypre->initializeSolverState(d_hierarchy,
      level_number);
   solv::PoissonSpecifications sps("Hypre Poisson solver");
   sps.setCZero();
   sps.setDConstant(1.0);
   d_poisson_hypre->setMatrixCoefficients(sps);

   /*
    * Solve the system.
    */
   tbox::plog << "solving..." << std::endl;
   int solver_ret;
   solver_ret = d_poisson_hypre->solveSystem(d_comp_soln_id,
         d_rhs_id, false, true);
   /*
    * Present data on the solve.
    */
   tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
              << "      iterations: " << d_poisson_hypre->getNumberOfIterations()
              << "\n"
              << "      residual: " << d_poisson_hypre->getRelativeResidualNorm()
              << "\n"
              << std::flush;

   /*
    * Deallocate state.
    */
   d_poisson_hypre->deallocateSolverState();

   /*
    * Return whether solver converged.
    */
   return solver_ret ? true : false;
}

#ifdef HAVE_HDF5
/*
 *************************************************************************
 * Set up VisIt to plot internal data from this class.
 * Tell the plotter about the refinement ratios.  Register variables
 * appropriate for plotting.
 *************************************************************************
 */
int HyprePoisson::registerVariablesWithPlotter(
   appu::VisItDataWriter& visit_writer) const {

   /*
    * This must be done once.
    */
   if (!d_hierarchy) {
      TBOX_ERROR(
         d_object_name << ": No hierarchy in\n"
                       << " HyprePoisson::registerVariablesWithPlotter\n"
                       << "The hierarchy must be built before calling\n"
                       << "this function.\n");
   }
   /*
    * Register variables with plotter.
    */
   visit_writer.registerPlotQuantity("Computed solution",
      "SCALAR",
      d_comp_soln_id);
   visit_writer.registerDerivedPlotQuantity("Error",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);
   visit_writer.registerPlotQuantity("Exact solution",
      "SCALAR",
      d_exact_id);
   visit_writer.registerPlotQuantity("Poisson source",
      "SCALAR",
      d_rhs_id);
   visit_writer.registerDerivedPlotQuantity("Patch level number",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);

   return 0;
}
#endif

/*
 *************************************************************************
 * Write derived data to the given stream.
 *************************************************************************
 */
bool HyprePoisson::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(region);
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   pdat::CellData<double>::iterator icell(pdat::CellGeometry::begin(patch.getBox()));
   pdat::CellData<double>::iterator icellend(pdat::CellGeometry::end(patch.getBox()));

   if (variable_name == "Error") {
      std::shared_ptr<pdat::CellData<double> > current_solution_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_comp_soln_id)));
      std::shared_ptr<pdat::CellData<double> > exact_solution_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_exact_id)));
      TBOX_ASSERT(current_solution_);
      TBOX_ASSERT(exact_solution_);
      pdat::CellData<double>& current_solution = *current_solution_;
      pdat::CellData<double>& exact_solution = *exact_solution_;
      for ( ; icell != icellend; ++icell) {
         double diff = (current_solution(*icell) - exact_solution(*icell));
         *buffer = diff;
         buffer += 1;
      }
   } else if (variable_name == "Patch level number") {
      double pln = patch.getPatchLevelNumber();
      for ( ; icell != icellend; ++icell) {
         *buffer = pln;
         buffer += 1;
      }
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        << "HyprePoissonX::writeDerivedDataToStream");

   }
   return true;
}

}

#endif
