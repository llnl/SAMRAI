/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Poisson solver
 *
 ************************************************************************/
#include "FACPoisson.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
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
FACPoisson::FACPoisson(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<solv::CellPoissonFACSolver>& fac_solver,
   const std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs):
   d_object_name(object_name),
   d_dim(dim),
   d_poisson_fac_solver(fac_solver),
   d_bc_coefs(bc_coefs)
{

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();

   /*
    * Get a unique context for variables owned by this object.
    */
   d_context = vdb->getContext(d_object_name + ":Context");

   tbox::ResourceAllocator allocator =
      tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator();

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    */

   std::shared_ptr<pdat::CellVariable<double> > comp_soln(
      new pdat::CellVariable<double>(
         dim,
         object_name + ":computed solution",
         allocator));
   d_comp_soln_id =
      vdb->registerVariableAndContext(
         comp_soln,
         d_context,
         hier::IntVector(dim, 1) /* ghost cell width is 1 for stencil widths */);

   std::shared_ptr<pdat::CellVariable<double> > exact_solution(
      new pdat::CellVariable<double>(
         dim,
         object_name + ":exact solution",
         allocator));
   d_exact_id =
      vdb->registerVariableAndContext(
         exact_solution,
         d_context,
         hier::IntVector(dim, 1) /* ghost cell width is 1 in case needed */);

   std::shared_ptr<pdat::CellVariable<double> > rhs_variable(
      new pdat::CellVariable<double>(
         dim,
         object_name
         + ":linear system right hand side",
         allocator));
   d_rhs_id =
      vdb->registerVariableAndContext(
         rhs_variable,
         d_context,
         hier::IntVector(dim, 0) /* ghost cell width is 0 */);

   /*
    * Specify an implementation of solv::RobinBcCoefStrategy for the solver to use.
    * We use the implementation solv::LocationIndexRobinBcCoefs, but other
    * implementations are possible, including user-implemented.
    */
   d_poisson_fac_solver->setBcObject(d_bc_coefs.get());
}

/*
 *************************************************************************
 * Destructor does nothing interesting
 *************************************************************************
 */
FACPoisson::~FACPoisson()
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
void FACPoisson::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
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

   std::shared_ptr<hier::PatchHierarchy> hierarchy = patch_hierarchy;

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

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
            patch_geom->getDx(),
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
            patch_geom->getDx(),
            patch_geom->getXLower());
      }

   }    // End patch loop.
}

/*
 *************************************************************************
 * Reset the hierarchy-dependent internal information.
 *************************************************************************
 */
void FACPoisson::resetHierarchyConfiguration(
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
 * Set up the initial guess and problem parameters
 * and solve the Poisson problem.  We explicitly initialize and
 * deallocate the solver state in this example.
 *************************************************************************
 */
int FACPoisson::solvePoisson()
{

   if (!d_hierarchy) {
      TBOX_ERROR(d_object_name
         << "Cannot solve using an uninitialized object.\n");
   }

   int ln;
   /*
    * Fill in the initial guess.
    */
   for (ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         std::shared_ptr<pdat::CellData<double> > data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_comp_soln_id)));
         TBOX_ASSERT(data);
         data->fill(0.0);
      }
   }

   /*
    * Set the parameters for the Poisson equation.
    * See classes solv::CellPoissonFACSolver or
    * solv::PoissonSpecifications.
    * (D is the diffusion coefficient.
    * C is the source term which is not used in this example.)
    */
   d_poisson_fac_solver->setDConstant(1.0);
   d_poisson_fac_solver->setCConstant(0.0);

   d_poisson_fac_solver->initializeSolverState(
      d_comp_soln_id,
      d_rhs_id,
      d_hierarchy,
      0,
      d_hierarchy->getFinestLevelNumber());

   tbox::plog << "solving..." << std::endl;
   int solver_ret;
   solver_ret = d_poisson_fac_solver->solveSystem(d_comp_soln_id, d_rhs_id);
   /*
    * Present data on the solve.
    */
   double avg_factor, final_factor;
   d_poisson_fac_solver->getConvergenceFactors(avg_factor, final_factor);
   tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
              << "      iterations: "
              << d_poisson_fac_solver->getNumberOfIterations() << "\n"
              << "      residual: " << d_poisson_fac_solver->getResidualNorm()
              << "\n"
              << "      average convergence: " << avg_factor << "\n"
              << "      final convergence: " << final_factor << "\n"
              << std::flush;

   d_poisson_fac_solver->deallocateSolverState();

   return 0;
}

#ifdef HAVE_HDF5
/*
 *************************************************************************
 * Set up external plotter to plot internal data from this class.
 * Register variables appropriate for plotting.
 *************************************************************************
 */
int FACPoisson::setupPlotter(
   appu::VisItDataWriter& plotter) const {
   if (!d_hierarchy) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                               << " FACPoisson::setupPlotter\n"
                               << "The hierarchy must be set before calling\n"
                               << "this function.\n");
   }
   plotter.registerPlotQuantity("Computed solution",
      "SCALAR",
      d_comp_soln_id);
   plotter.registerDerivedPlotQuantity("Error",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);
   plotter.registerPlotQuantity("Exact solution",
      "SCALAR",
      d_exact_id);
   plotter.registerPlotQuantity("Poisson source",
      "SCALAR",
      d_rhs_id);

   return 0;
}
#endif

/*
 *************************************************************************
 * Write derived data to the given stream.
 *************************************************************************
 */
bool FACPoisson::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   pdat::CellData<double>::iterator icell(pdat::CellGeometry::begin(region));
   pdat::CellData<double>::iterator icellend(pdat::CellGeometry::end(region));

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
         buffer = buffer + 1;
      }
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        << "FACPoissonX::packDerivedDataIntoDoubleBuffer");

   }
   // Return true if this patch has derived data on it.
   // False otherwise.
   return true;
}

}
