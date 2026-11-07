/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class containing numerical routines for modified Bratu problem
 *
 ************************************************************************/

#include "ModifiedBratuProblem.h"
#include "ModifiedBratuFort.h"

#if defined(HAVE_PETSC) && defined(HAVE_SUNDIALS) && defined(HAVE_HYPRE)

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif

using namespace SAMRAI;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <float.h>
#include <vector>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"

#define MODIFIED_BRATU_PROBLEM (1)

// Define for number of ghost cells on solution quantity
#define NUM_GHOSTS_U (1)

std::shared_ptr<tbox::Timer> ModifiedBratuProblem::s_copy_timer;
std::shared_ptr<tbox::Timer> ModifiedBratuProblem::s_pc_timer;

/*
 *************************************************************************
 *
 * Constructor and destructor for Modified Bratu Problem class.
 *
 * The problem we solve is:
 *
 *     du/dt = div( D(x,t)*grad(u) ) + lambda * exp(u) + f(u,x,t).
 *
 * The unknown u is a cell-centered variable.  The time discretization
 * uses the backward Euler strategy.  We use a standard 7-point stencil
 * (based on the finite volume meethod) for the div(D*grad(u)) term.
 *
 * The variables used here to manage data for the discrete problem are:
 *
 *    solution ........... unknown quantity "u"
 *    source_term ........ source term "f"
 *    exponential_term ... product "lambda * exp(u)"
 *    diffusion_coef ..... diffusion coefficient "D"
 *
 * Other quantities used in the solution process:
 *
 *    weight.............. weights for solution vector entries on grid
 *    flux................ side-centered fluxes "D * grad(u)"
 *    coarse_fine_flux.... fluxes at coarse-fine interfaces
 *    jacobian_a.......... Jacobian entries for FAC solver
 *    jacobian_b.......... Jacobian entries for FAC solver
 *
 * The constructor creates these variables to represent the solution
 * and other quantities on the patch hierarchy.  The solution quantity
 * is managed using three variable contexts, "CURRENT", "NEW", and
 * "SCRATCH".  CURRENT and NEW represent the current and new values of
 * the solution, respectively.  They have no ghost cells.  The SCRATCH
 * context is where most of the computations involving the spatial
 * discretization stencil occur.  This storage has a ghost cell width
 * of one.  All other quantities are managed using only the SCRATCH
 * context.  However, they all have zero ghost cell widths.
 *
 *************************************************************************
 */

ModifiedBratuProblem::ModifiedBratuProblem(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<solv::CellPoissonFACSolver> fac_solver,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
   std::shared_ptr<appu::VisItDataWriter> visit_writer):
   RefinePatchStrategy(),
   CoarsenPatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_grid_geometry(grid_geometry),
   d_lambda(tbox::MathUtilities<double>::getSignalingNaN()),
   d_input_dt(tbox::MathUtilities<double>::getSignalingNaN()),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_solution(new pdat::CellVariable<double>(
                 dim, object_name + "solution", d_allocator)),
   d_source_term(new pdat::CellVariable<double>(
                    dim, object_name + "source_term", d_allocator)),
   d_exponential_term(new pdat::CellVariable<double>(
                         dim, object_name + "exponential_term", d_allocator)),
   d_diffusion_coef(new pdat::SideVariable<double>(
                       dim, object_name + "diffusion_coef",
                       hier::IntVector::getOne(dim), d_allocator)),
   d_flux(new pdat::SideVariable<double>(
             dim, object_name + "flux",
             hier::IntVector::getOne(dim), d_allocator)),
   d_coarse_fine_flux(new pdat::OutersideVariable<double>(
                         dim, object_name + "coarse_fine_flux", d_allocator)),
   d_jacobian_a(new pdat::CellVariable<double>(
                   dim, object_name + ":jacobian_a", d_allocator)),
   d_jacobian_b(new pdat::FaceVariable<double>(
                   dim, object_name + ":jacobian_b", d_allocator)),
   d_jacobian_a_id(-1),
   d_jacobian_b_id(-1),
   d_precond_a(new pdat::CellVariable<double>(
                  dim, object_name + "precond_a", d_allocator)),
   d_precond_b(new pdat::FaceVariable<double>(
                  dim, object_name + "precond_b", d_allocator)),
   d_precond_a_id(-1),
   d_precond_b_id(-1),
   d_nghosts(hier::IntVector(dim, NUM_GHOSTS_U)),
   d_weight(new pdat::CellVariable<double>(
               dim, object_name + "weight", d_allocator, 1)),
   d_fill_new_level(),
   d_soln_fill(),
   d_flux_coarsen(dim),
   d_soln_coarsen(dim),
   d_scratch_soln_coarsen(dim),
   d_current_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_new_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_current_dt(tbox::MathUtilities<double>::getSignalingNaN()),
   d_FAC_solver(fac_solver),
   d_max_precond_its(tbox::MathUtilities<int>::getMax()),
   d_precond_tol(tbox::MathUtilities<double>::getSignalingNaN())
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(grid_geometry);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   getFromInput(input_db, false);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_current = variable_db->getContext("CURRENT");
   d_new = variable_db->getContext("NEW");
   d_scratch = variable_db->getContext("SCRATCH");

   /*
    * Create variables used in the discrete problem.  Register variables
    * and contexts with variable database to retrieve their descriptor
    * indices for data management.
    */

   int soln_id = variable_db->registerVariableAndContext(d_solution,
         d_current,
         hier::IntVector(d_dim, 0));
   if (visit_writer) {
      visit_writer->registerPlotQuantity("U", "SCALAR", soln_id);
   }
   int soln_new_id = variable_db->registerVariableAndContext(d_solution,
         d_new,
         hier::IntVector(d_dim, 0));
   d_soln_scratch_id = variable_db->registerVariableAndContext(d_solution,
         d_scratch,
         d_nghosts);

   /*
    * *hier::Variable to weight solution vector entries on a composite grid.
    */

   d_weight_id = variable_db->registerVariableAndContext(d_weight,
         d_scratch,
         hier::IntVector(d_dim, 0));

   /*
    * Other variables used in discrete problem.
    */

   int source_id = variable_db->registerVariableAndContext(d_source_term,
         d_scratch,
         hier::IntVector(d_dim, 0));

   int exp_id = variable_db->registerVariableAndContext(d_exponential_term,
         d_scratch,
         hier::IntVector(d_dim, 0));

   int diffcoef_id = variable_db->registerVariableAndContext(d_diffusion_coef,
         d_scratch,
         hier::IntVector(d_dim, 0));

   d_flux_id = variable_db->registerVariableAndContext(d_flux,
         d_scratch,
         hier::IntVector(d_dim, 0));

   d_coarse_fine_flux_id =
      variable_db->registerVariableAndContext(d_coarse_fine_flux,
         d_scratch,
         hier::IntVector(d_dim, 0));

   /*
    * Variables for A(x)z=r preconditioning.
    */
   d_precond_a_id = variable_db->registerVariableAndContext(d_precond_a,
         d_scratch,
         hier::IntVector(d_dim, 0));
   d_precond_b_id = variable_db->registerVariableAndContext(d_precond_b,
         d_scratch,
         hier::IntVector(d_dim, 0));

   /*
    * Variables for A(x)*v operations.
    */
   d_jacobian_a_id = variable_db->registerVariableAndContext(d_jacobian_a,
         d_scratch,
         hier::IntVector(d_dim, 0));
   d_jacobian_b_id = variable_db->registerVariableAndContext(d_jacobian_b,
         d_scratch,
         hier::IntVector(d_dim, 0));

   /*
    * ComponentSelectors are used to allocate and deallocate storage for
    * variables on the patch hierarchy in groups.  This reduces the
    * complexity of tracking each individual data index.   Storage
    * for the time evolving solution is managed collectively and
    * auxiliary quantities used in nonlinear function evaluation.  These
    * latter quantities could be placed in a separate hier::ComponentSelector
    * and allocated only when needed, but we chose not to add that
    * for this simple problem.  Data needed to correctly compute
    * fluxes at coarse/fine interfaces is also managed collectively
    * and is only allocated when needed.  Finally data to form the
    * Jacobian is also grouped together.
    */

   d_new_patch_problem_data.setFlag(soln_id);

   d_problem_data.setFlag(soln_new_id);
   d_problem_data.setFlag(d_soln_scratch_id);

   d_problem_data.setFlag(d_weight_id);

   d_problem_data.setFlag(source_id);
   d_problem_data.setFlag(exp_id);
   d_problem_data.setFlag(diffcoef_id);

   d_jacobian_data.setFlag(d_jacobian_a_id);
   d_jacobian_data.setFlag(d_jacobian_b_id);

   d_precond_data.setFlag(d_precond_a_id);
   d_precond_data.setFlag(d_precond_b_id);

   /*
    * Refine and coarsen algorithms used to solve the problem.
    *
    *    d_fill_new_level :
    *
    *       When new levels are created in the hierarchy (e.g., after
    *       regridding), the "CURRENT" solution state is initialized in
    *       new fine regions by using linear interpolation from coarser
    *       levels.
    *
    *    d_flux_coarsen:
    *
    *       Used to coarsen fluxes at coarse-fine interfaces during
    *       nonlinear function evaluation.
    */

   d_soln_refine_op =
      grid_geometry->lookupRefineOperator(d_solution,
         "LINEAR_REFINE");

   d_fill_new_level.registerRefine(soln_id,
      soln_id,
      d_soln_scratch_id,
      d_soln_refine_op);

   d_soln_coarsen_op =
      grid_geometry->lookupCoarsenOperator(d_solution,
         "CONSERVATIVE_COARSEN");

   std::shared_ptr<CoarsenOperator> flux_coarsen_op(
      grid_geometry->lookupCoarsenOperator(d_coarse_fine_flux,
         "CONSERVATIVE_COARSEN"));

   /*
    * Register the transfer algorithms so we can create
    * schedules when we have info on the hierarchy (in
    * resetHierarchyConfiguration).
    */

   d_scratch_soln_coarsen.registerCoarsen(d_soln_scratch_id,
      d_soln_scratch_id,
      d_soln_coarsen_op);

   d_soln_coarsen.registerCoarsen(soln_id,
      soln_id,
      d_soln_coarsen_op);

   d_flux_coarsen.registerCoarsen(d_flux_id,
      d_coarse_fine_flux_id,
      flux_coarsen_op);

   d_soln_fill.registerRefine(d_soln_scratch_id,
      soln_new_id,
      d_soln_scratch_id,
      d_soln_refine_op);

   s_copy_timer = tbox::TimerManager::getManager()->
      getTimer("apps::usrFcns::evaluateBratuFunction");
   s_pc_timer = tbox::TimerManager::getManager()
      ->getTimer("apps::usrFcns::applyBratuPreconditioner");

}

/*
 *************************************************************************
 *
 * The destructor does nothing interesting.
 *
 *************************************************************************
 */

ModifiedBratuProblem::~ModifiedBratuProblem()
{
   s_copy_timer.reset();
   s_pc_timer.reset();
}

/*
 *************************************************************************
 *
 * The vector used in the nonlinear iteration to advance the data
 * includes only the "NEW" solution values.  In other words, we
 * know the "CURRENT" solution and we advance to the "NEW" solution.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::setupSolutionVector(
   const std::shared_ptr<solv::SAMRAIVectorReal<double> >& solution)
{
   TBOX_ASSERT(solution);

   d_solution_vector = solution;

   solution->addComponent(d_solution,
      hier::VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(d_solution, d_new),
      d_weight_id);
}

/*
 * Vector weights are used to calculate inner products and norms.
 * In cells not covered by finer cells, the weights are set to the
 * the cell volume.  In those covered by finer cells, the weights
 * are set to 0.  This latter assignment ensures that data in
 * covered cells contribute to these sums only once.
 */

void ModifiedBratuProblem::setVectorWeights(
   std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {

      /*
       * On every level, first assign cell volume to vector weight.
       */

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geometry);
         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (d_dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
         }
         if (d_dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
         }
         std::shared_ptr<pdat::CellData<double> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_weight_id)));
         TBOX_ASSERT(w);
         w->fillAll(cell_vol);
      }

      /*
       * On all but the finest level, assign 0 to vector
       * weight to cells covered by finer cells.
       */

      if (amr_level < hierarchy->getFinestLevelNumber()) {

         /*
          * First get the boxes that describe index space of the next finer
          * level and coarsen them to describe corresponding index space
          * at this level.
          */

         std::shared_ptr<hier::PatchLevel> next_finer_level(
            hierarchy->getPatchLevel(amr_level + 1));
         hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
         hier::IntVector coarsen_ratio = next_finer_level->getRatioToLevelZero();
         coarsen_ratio /= level->getRatioToLevelZero();
         coarsened_boxes.coarsen(coarsen_ratio);

         /*
          * Then set vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for (hier::PatchLevel::iterator p(level->begin());
              p != level->end(); ++p) {

            const std::shared_ptr<hier::Patch>& patch = *p;
            for (hier::BoxContainer::const_iterator i = coarsened_boxes.begin();
                 i != coarsened_boxes.end(); ++i) {

               const hier::Box& coarse_box = *i;
               hier::Box intersection = coarse_box * patch->getBox();
               if (!intersection.empty()) {
                  std::shared_ptr<pdat::CellData<double> > w(
                     SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(d_weight_id)));
                  TBOX_ASSERT(w);
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
         }  // loop over patches in level
      }  // all levels except finest
   }  // loop over levels
}

/*
 *************************************************************************
 *
 * Set the initial guess for the nonlinear iteration.  If we are at the
 * first step on the current hierarchy configuration, we initialize the
 * the new solution with the current soltuion values.  Note that we do
 * the same if we are not at the first step on the current hierarchy
 * configuration.  In the future, we may employ a more sophisticated
 * strategy involving some sort of extrapolation.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::setInitialGuess(
   const bool first_step,
   const double current_time,
   const double current_dt,
   const double old_dt)
{
   NULL_USE(old_dt);

   d_current_time = current_time;
   d_current_dt = current_dt;
   d_new_time = d_current_time + d_current_dt;

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      hierarchy->getPatchLevel(ln)->
      allocatePatchData(d_jacobian_data);
      hierarchy->getPatchLevel(ln)->
      allocatePatchData(d_precond_data);
   }

   if (first_step) {

      for (int amr_level = 0;
           amr_level < hierarchy->getNumberOfLevels();
           ++amr_level) {

         std::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(amr_level));
         for (hier::PatchLevel::iterator p(patch_level->begin());
              p != patch_level->end(); ++p) {
            const std::shared_ptr<hier::Patch>& patch = *p;

            std::shared_ptr<pdat::CellData<double> > y_cur(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(d_solution, d_current)));
            std::shared_ptr<pdat::CellData<double> > y_new(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(d_solution, d_new)));
            TBOX_ASSERT(y_cur);
            TBOX_ASSERT(y_new);
            y_new->copy(*y_cur);

            y_new->setTime(d_new_time);

            patch->getPatchData(d_solution, d_scratch)->setTime(d_new_time);

            std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
               SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                  patch->getPatchGeometry()));
            TBOX_ASSERT(patch_geometry);
            const double* dx = patch_geometry->getDx();
            const double* xlo = patch_geometry->getXLower();
            const double* xhi = patch_geometry->getXUpper();
            const hier::Index ifirst = patch->getBox().lower();
            const hier::Index ilast = patch->getBox().upper();

            std::shared_ptr<pdat::SideData<double> > diffusion(
               SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                  patch->getPatchData(d_diffusion_coef, d_scratch)));
            TBOX_ASSERT(diffusion);

            if (d_dim == tbox::Dimension(1)) {
               FORT_EVALDIFFUSION1D(ifirst(0), ilast(0),
                  dx, xlo, xhi,
                  diffusion->getPointer(0));
            } else if (d_dim == tbox::Dimension(2)) {
               FORT_EVALDIFFUSION2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  dx, xlo, xhi,
                  diffusion->getPointer(0),
                  diffusion->getPointer(1));
            } else if (d_dim == tbox::Dimension(3)) {
               FORT_EVALDIFFUSION3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  dx, xlo, xhi,
                  diffusion->getPointer(0),
                  diffusion->getPointer(1),
                  diffusion->getPointer(2));
            }
         }
      }
   } else {

      for (int amr_level = 0;
           amr_level < hierarchy->getNumberOfLevels();
           ++amr_level) {

         std::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(amr_level));
         for (hier::PatchLevel::iterator p(patch_level->begin());
              p != patch_level->end(); ++p) {
            const std::shared_ptr<hier::Patch>& patch = *p;

            std::shared_ptr<pdat::CellData<double> > y_cur(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(d_solution, d_current)));
            std::shared_ptr<pdat::CellData<double> > y_new(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(d_solution, d_new)));
            TBOX_ASSERT(y_cur);
            TBOX_ASSERT(y_new);
            y_new->copy(*y_cur);

            y_new->setTime(d_new_time);

            patch->getPatchData(d_solution, d_scratch)->setTime(d_new_time);

            std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
               SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                  patch->getPatchGeometry()));
            TBOX_ASSERT(patch_geometry);
            const double* dx = patch_geometry->getDx();
            const double* xlo = patch_geometry->getXLower();
            const double* xhi = patch_geometry->getXUpper();
            const hier::Index ifirst = patch->getBox().lower();
            const hier::Index ilast = patch->getBox().upper();

            std::shared_ptr<pdat::SideData<double> > diffusion(
               SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                  patch->getPatchData(d_diffusion_coef, d_scratch)));
            TBOX_ASSERT(diffusion);
            if (d_dim == tbox::Dimension(1)) {
               FORT_EVALDIFFUSION1D(ifirst(0), ilast(0),
                  dx, xlo, xhi,
                  diffusion->getPointer(0));
            } else if (d_dim == tbox::Dimension(2)) {
               FORT_EVALDIFFUSION2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  dx, xlo, xhi,
                  diffusion->getPointer(0),
                  diffusion->getPointer(1));
            } else if (d_dim == tbox::Dimension(3)) {
               FORT_EVALDIFFUSION3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  dx, xlo, xhi,
                  diffusion->getPointer(0),
                  diffusion->getPointer(1),
                  diffusion->getPointer(2));
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Return the time increment used for the first solution advance step.
 *
 *************************************************************************
 */

double ModifiedBratuProblem::getInitialDt()
{
   return d_input_dt;
}

/*
 *************************************************************************
 *
 * Compute and return the next time increment through which to advance
 * the solution.  Note that good_solution is the value returned by the
 * preceeding call to checkNewSolution().  The integer solver_retcode
 * is the return code generated by the nonlinear solver.
 *
 * In the future, we will have a more intelligent strategy.  When
 * good_solution is true, we would like to take the largest timestep
 * possible.   When good_solution is false, we will cut the timestep
 * depending on the solver return code.  In this case, if the code
 * indicates that that nonlinear iteration converged, we will cut the
 * timestep by some amount since the temporal error is probably too
 * large.  For a first order method like BE, determining a suitable
 * timestep can be based on an estimate of the second time derivative
 * of the solution.  See somments before checkNewSolution() below.
 * If the nonlinear iteration did not converge, we will have to try
 * another approach.
 *
 * At this time, we always return a constant time step read from input.
 *
 *************************************************************************
 */

double ModifiedBratuProblem::getNextDt(
   const bool good_solution,
   const int solver_retcode)
{
   NULL_USE(good_solution);
   NULL_USE(solver_retcode);

   return d_input_dt;
}

/*
 *************************************************************************
 *
 * Check the computed solution and return a boolean value of true if it
 * is acceptable; otherwise return false.  The integer solver_retcode
 * is the return code generated by the nonlinear solver.  This value
 * must be interpreted in a manner consistant with the solver in use.
 *
 * Ordinarily we would estimate the temporal error in the solution
 * before accepting it.  For backward Euler, this is straightforward.
 * For example, if we had both (d_solution, d_current) and (d_solution,
 * d_old), we could evaluate dy/dt at these two timesteps by evaluating
 * the rhs at these values.  Then the second derivative in the solution
 * can be estimated by using finite differences and compared to some
 * prescribed tolerance.  Other heuristics may be possible.
 *
 * For now, we accept every solution we compute.
 *
 *************************************************************************
 */

bool ModifiedBratuProblem::checkNewSolution(
   const int solver_retcode)
{
   NULL_USE(solver_retcode);

   double new_time = d_current_time + d_current_dt;
   double maxerror = 0.0;
   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());
   for (int amr_level = 0;
        amr_level < hierarchy->getNumberOfLevels();
        ++amr_level) {
      std::shared_ptr<hier::PatchLevel> patch_level(
         hierarchy->getPatchLevel(amr_level));
      double levelerror = 0.0;
      double levell2error = 0.0;
      for (hier::PatchLevel::iterator p(patch_level->begin());
           p != patch_level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;
         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();
         const double* xlo = patch_geom->getXLower();
         const double* xhi = patch_geom->getXUpper();
         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_new)));
         std::shared_ptr<pdat::CellData<double> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_weight_id)));
         TBOX_ASSERT(u);
         TBOX_ASSERT(w);

         if (d_dim == tbox::Dimension(1)) {
            FORT_ERROR1D(ifirst(0), ilast(0),
               u->getPointer(), w->getPointer(),
               d_lambda,
               xlo, xhi, dx,
               new_time,
               levelerror, levell2error);
         } else if (d_dim == tbox::Dimension(2)) {
            FORT_ERROR2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               u->getPointer(), w->getPointer(),
               d_lambda,
               xlo, xhi, dx,
               new_time,
               levelerror, levell2error);
         } else if (d_dim == tbox::Dimension(3)) {
            FORT_ERROR3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               u->getPointer(), w->getPointer(),
               d_lambda,
               xlo, xhi, dx,
               new_time,
               levelerror, levell2error);
         }
      }
      maxerror = (levelerror > maxerror) ? levelerror : maxerror;
      tbox::plog << " At " << new_time
                 << " on level " << amr_level
                 << " max err is " << levelerror
                 << " wtd l2 err is " << levell2error
                 << std::endl;
   }

   tbox::plog << " At " << new_time << " err is " << maxerror << std::endl;

   return true;
}

/*
 *************************************************************************
 *
 * Update soltution quantities after computing an acceptable time
 * advanced solution.   The new_time value is the new solution time.
 *
 * Since we have no dependent variables in this problem and our memory
 * management is extremely simple, we just copy new into current.
 * Then, we synchronize the times for each set of solution vlaues.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::updateSolution(
   const double new_time)
{
   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   d_new_time = d_current_time = new_time;

   for (int amr_level = 0;
        amr_level < hierarchy->getNumberOfLevels();
        ++amr_level) {

      std::shared_ptr<hier::PatchLevel> patch_level(
         hierarchy->getPatchLevel(amr_level));
      for (hier::PatchLevel::iterator p(patch_level->begin());
           p != patch_level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         std::shared_ptr<pdat::CellData<double> > y_cur(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_current)));
         std::shared_ptr<pdat::CellData<double> > y_new(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_new)));
         TBOX_ASSERT(y_cur);
         TBOX_ASSERT(y_new);
         y_cur->copy(*y_new);

         y_cur->setTime(new_time);

         patch->getPatchData(d_solution, d_scratch)->setTime(new_time);

      }

      hierarchy->getPatchLevel(amr_level)->
      deallocatePatchData(d_jacobian_data);
      hierarchy->getPatchLevel(amr_level)->
      deallocatePatchData(d_precond_data);

   }

}

/*
 *************************************************************************
 *
 * Allocate and initialize data for a new level in the patch hierarchy.
 *
 * At initial time only, we initialize the solution to zero.  At all
 * other times the solution is either copied from an old level or
 * interpolated from coarser levels.  The cell weights and diffusion
 * coefficients are always set here since they live as long as a
 * patch and are not time dependent.
 *
 *************************************************************************
 */
void ModifiedBratuProblem::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(allocate_data);

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(!old_level || level_number == old_level->getLevelNumber());
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   level->allocatePatchData(d_new_patch_problem_data, time);

   if ((level_number > 0) || old_level) {

      std::shared_ptr<RefineSchedule> sched(
         d_fill_new_level.createSchedule(
            level,
            old_level,
            level_number - 1,
            hierarchy,
            0));
      sched->fillData(time);

      if (old_level) {
         old_level->deallocatePatchData(d_new_patch_problem_data);
         old_level->deallocatePatchData(d_problem_data);
      }
   }

   level->allocatePatchData(d_problem_data);

   /*
    * Initialize data on patch interiors.  At initial time only, we
    * initialize solution to zero.  At all other times the solution
    * is either copied from an old level, or interpolated from some
    * coarser level using the "d_fill_new_level" algorithm above.
    */

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      if (initial_time) {

         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_current)));
         TBOX_ASSERT(u);
         u->fillAll(0.0);

      }

      std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));
      TBOX_ASSERT(patch_geometry);
      const double* dx = patch_geometry->getDx();
      const double* xlo = patch_geometry->getXLower();
      const double* xhi = patch_geometry->getXUpper();
      const hier::Index ifirst = patch->getBox().lower();
      const hier::Index ilast = patch->getBox().upper();

      std::shared_ptr<pdat::SideData<double> > diffusion(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch->getPatchData(d_diffusion_coef, d_scratch)));
      TBOX_ASSERT(diffusion);
      if (d_dim == tbox::Dimension(1)) {
         FORT_EVALDIFFUSION1D(ifirst(0), ilast(0),
            dx, xlo, xhi,
            diffusion->getPointer(0));
      } else if (d_dim == tbox::Dimension(2)) {
         FORT_EVALDIFFUSION2D(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            dx, xlo, xhi,
            diffusion->getPointer(0),
            diffusion->getPointer(1));
      } else if (d_dim == tbox::Dimension(3)) {
         FORT_EVALDIFFUSION3D(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            ifirst(2), ilast(2),
            dx, xlo, xhi,
            diffusion->getPointer(0),
            diffusion->getPointer(1),
            diffusion->getPointer(2));
      }
   }
}

/*
 *************************************************************************
 *
 * After hierarchy levels whose numbers lie in the given range have
 * been either added or changed in the hierarchy, this routine is
 * called.  Typically, it is used to reset any data members that
 * remain constant until the levels change again via adaptive gridding
 * (e.g., communication schedules), thus allowing some efficiencies to
 * be achieved.  Here, we do nothing, since the function evalutation
 * and Jacobian-vector product routines have indeterminate data ids.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::resetHierarchyConfiguration(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level <= finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int ln0 = 0; ln0 <= finest_level; ++ln0) {
      TBOX_ASSERT(hierarchy->getPatchLevel(ln0));
   }
#endif

   d_flux_coarsen_schedule.resize(hierarchy->getNumberOfLevels());
   d_soln_fill_schedule.resize(hierarchy->getNumberOfLevels());
   d_soln_coarsen_schedule.resize(hierarchy->getNumberOfLevels());
   d_scratch_soln_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   int ln;
   /*
    * Rebuild schedules affected by the hierarchy change.
    * In the following loops, ln is the level number of
    * the destination level affected by the change.
    */
   {
      const int ln_beg = coarsest_level - (coarsest_level > 0);
      const int ln_end = finest_level;
      for (ln = ln_beg; ln < ln_end; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         std::shared_ptr<hier::PatchLevel> finer_level(
            hierarchy->getPatchLevel(ln + 1));

         d_flux_coarsen_schedule[ln] =
            d_flux_coarsen.createSchedule(level,
               finer_level);
      }
   }
   {
      const int ln_beg = coarsest_level;
      const int ln_end = finest_level;
      for (ln = ln_beg; ln <= ln_end; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         d_soln_fill_schedule[ln] = d_soln_fill.createSchedule(level,
               ln - 1,
               hierarchy,
               this);
      }
   }
   {
      const int ln_beg = coarsest_level - (coarsest_level > 0);
      const int ln_end = finest_level;
      for (ln = ln_beg; ln < ln_end; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         std::shared_ptr<hier::PatchLevel> finer_level(
            hierarchy->getPatchLevel(ln + 1));

         d_soln_coarsen_schedule[ln] =
            d_soln_coarsen.createSchedule(level,
               finer_level);
         d_scratch_soln_coarsen_schedule[ln] =
            d_scratch_soln_coarsen.createSchedule(level,
               finer_level);
      }
   }
}

/*
 *************************************************************************
 *
 * KINSOL interface for user-supplied routines.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::evaluateNonlinearFunction(
   solv::SundialsAbstractVector* soln,
   solv::SundialsAbstractVector* fval)
{
   TBOX_ASSERT(soln != 0);
   TBOX_ASSERT(fval != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > x(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(soln));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > f(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(fval));

   evaluateBratuFunction(x, f);
}

int ModifiedBratuProblem::precondSetup(
   solv::SundialsAbstractVector* soln,
   solv::SundialsAbstractVector* soln_scale,
   solv::SundialsAbstractVector* fval,
   solv::SundialsAbstractVector* fval_scale,
   int& num_feval)
{
   NULL_USE(soln_scale);
   NULL_USE(fval);
   NULL_USE(fval_scale);

   TBOX_ASSERT(soln != 0);

   num_feval += 0;

   std::shared_ptr<solv::SAMRAIVectorReal<double> > x(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(soln));

   setupBratuPreconditioner(x);

   return 0;
}

int ModifiedBratuProblem::precondSolve(
   solv::SundialsAbstractVector* soln,
   solv::SundialsAbstractVector* soln_scale,
   solv::SundialsAbstractVector* fval,
   solv::SundialsAbstractVector* fval_scale,
   solv::SundialsAbstractVector* rhs,
   int& num_feval)
{
   NULL_USE(soln);
   NULL_USE(soln_scale);
   NULL_USE(fval);
   NULL_USE(fval_scale);

   TBOX_ASSERT(rhs != 0);

   num_feval += 0;

   std::shared_ptr<solv::SAMRAIVectorReal<double> > r(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(rhs));

   return applyBratuPreconditioner(r, r);
}

int
ModifiedBratuProblem::jacobianTimesVector(
   solv::SundialsAbstractVector* vector,
   solv::SundialsAbstractVector* product,
   const bool soln_changed,
   solv::SundialsAbstractVector* soln)
{
   TBOX_ASSERT(vector != 0);
   TBOX_ASSERT(product != 0);
   TBOX_ASSERT(soln != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > v(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(vector));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > Jv(
      solv::Sundials_SAMRAIVector::getSAMRAIVector(product));

   if (soln_changed) {
      std::shared_ptr<solv::SAMRAIVectorReal<double> > ucur(
         solv::Sundials_SAMRAIVector::getSAMRAIVector(soln));
      evaluateBratuJacobian(ucur);
   }
   return jacobianTimesVector(v, Jv);
}

/*
 *************************************************************************
 *
 * PETSc/SNES interface for user-supplied routines.
 *
 *************************************************************************
 */

int ModifiedBratuProblem::evaluateNonlinearFunction(
   Vec xcur,
   Vec fcur)
{
   TBOX_ASSERT(xcur != 0);
   TBOX_ASSERT(fcur != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > x(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(xcur));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > f(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(fcur));

   evaluateBratuFunction(x, f);

   return 0;
}

int ModifiedBratuProblem::evaluateJacobian(
   Vec x)
{
   TBOX_ASSERT(x != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > xvec(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(x));

   evaluateBratuJacobian(xvec);

   return 0;
}

int ModifiedBratuProblem::jacobianTimesVector(
   Vec xin,
   Vec xout)
{
   TBOX_ASSERT(xin != 0);
   TBOX_ASSERT(xout != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > xinvec(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(xin));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > xoutvec(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(xout));

   jacobianTimesVector(xinvec, xoutvec);

   return 0;
}

int ModifiedBratuProblem::setupPreconditioner(
   Vec x)
{
   std::shared_ptr<solv::SAMRAIVectorReal<double> > uvec(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(x));
   setupBratuPreconditioner(uvec);

   return 0;
}

int ModifiedBratuProblem::applyPreconditioner(
   Vec r,
   Vec z)
{
   TBOX_ASSERT(r != 0);
   TBOX_ASSERT(z != 0);

   std::shared_ptr<solv::SAMRAIVectorReal<double> > rhs(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(r));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > soln(
      solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(z));

   return applyBratuPreconditioner(rhs, soln);
}

/*
 *************************************************************************
 *
 * Generic user-supplied functions.  Here the interfaces are expressed
 * in terms of generic SAMRAIVectors.  These methods are invoked from
 * within a nonlinear solver package (through an appropriate interface).
 * Since we have no control over the vectors that appear as actual
 * arguments, each of these methods must create communication algorithms
 * and schedules to fill ghost cells before stencil operations can be
 * applied.
 *
 *************************************************************************
 */

/*
 *************************************************************************
 *
 * Evaluate nonlinear residual at the vector "x" and place the result in
 * the vector "f".
 *
 *************************************************************************
 */

void ModifiedBratuProblem::evaluateBratuFunction(
   std::shared_ptr<solv::SAMRAIVectorReal<double> > x,
   std::shared_ptr<solv::SAMRAIVectorReal<double> > f)
{
   TBOX_ASSERT(x);
   TBOX_ASSERT(f);

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   /*
    * Create a coarsen algorithm to coarsen data on fine patch interiors.
    * Then, cycle through the levels and coarsen the data.
    */
   CoarsenAlgorithm eval_average(d_dim);
   eval_average.registerCoarsen(x->getComponentDescriptorIndex(0),
      x->getComponentDescriptorIndex(0),
      d_soln_coarsen_op);
   int amr_level = 0;
   for (amr_level = hierarchy->getFinestLevelNumber() - 1;
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      std::shared_ptr<hier::PatchLevel> finer_level(
         hierarchy->getPatchLevel(amr_level + 1));
      eval_average.resetSchedule(d_soln_coarsen_schedule[amr_level]);
      d_soln_coarsen_schedule[amr_level]->coarsenData();
      d_soln_coarsen.resetSchedule(d_soln_coarsen_schedule[amr_level]);
   }

   /*
    * Create a refine algorithm to fill ghost cells of solution variable.
    */

   s_copy_timer->start();
   RefineAlgorithm eval_fill;
   eval_fill.registerRefine(d_soln_scratch_id,
      x->getComponentDescriptorIndex(0),
      d_soln_scratch_id,
      d_soln_refine_op);
   s_copy_timer->stop();

   for (amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));

      s_copy_timer->start();
      eval_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      d_soln_fill_schedule[amr_level]->fillData(d_new_time);
      d_soln_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      s_copy_timer->stop();

      level->allocatePatchData(d_flux_id);
      if (amr_level > 0) {
         level->allocatePatchData(d_coarse_fine_flux_id);
      }

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();
         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_diffusion_coef, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > flux(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_flux_id)));
         std::shared_ptr<pdat::OutersideData<double> > coarse_fine_flux(
            SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
               patch->getPatchData(d_coarse_fine_flux_id)));
         TBOX_ASSERT(u);
         TBOX_ASSERT(diffusion);
         TBOX_ASSERT(flux);
         if (d_dim == tbox::Dimension(1)) {
            FORT_EVALFACEFLUXES1D(ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               u->getPointer(),
               dx,
               flux->getPointer(0));
         } else if (d_dim == tbox::Dimension(2)) {
            FORT_EVALFACEFLUXES2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               u->getPointer(),
               dx,
               flux->getPointer(0),
               flux->getPointer(1));
         } else if (d_dim == tbox::Dimension(3)) {
            FORT_EVALFACEFLUXES3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               diffusion->getPointer(2),
               u->getPointer(),
               dx,
               flux->getPointer(0),
               flux->getPointer(1),
               flux->getPointer(2));
         }
         const std::vector<hier::BoundaryBox>& bdry_faces =
            patch_geom->getCodimensionBoundaries(1);
#if 0
         if (d_dim == tbox::Dimension(1)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getNodeBoundaries();
         } else if (d_dim == tbox::Dimension(2)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getEdgeBoundaries();
         } else if (d_dim == tbox::Dimension(3)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getFaceBoundaries();
         }
#endif
         for (int i = 0; i < static_cast<int>(bdry_faces.size()); ++i) {

            hier::Box bbox = bdry_faces[i].getBox();
            const hier::Index ibeg = bbox.lower();
            const hier::Index iend = bbox.upper();
            int face = bdry_faces[i].getLocationIndex();

            if (d_dim == tbox::Dimension(1)) {
               FORT_EWBCFLUXFIX1D(ifirst(0), ilast(0),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            } else if (d_dim == tbox::Dimension(2)) {
               FORT_EWBCFLUXFIX2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_NSBCFLUXFIX2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(1),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            } else if (d_dim == tbox::Dimension(3)) {
               FORT_EWBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_NSBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(1),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_TBBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  u->getPointer(),
                  flux->getPointer(2),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            }

         } // end loop over boundary faces

         // correctPatchFlux( level, patch, u );

         if (amr_level > 0) {

            TBOX_ASSERT(coarse_fine_flux);

            for (int side = 0; side < 2; ++side) {
               if (d_dim == tbox::Dimension(1)) {
                  FORT_EWFLUXCOPY1D(ifirst(0), ilast(0),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
               } else if (d_dim == tbox::Dimension(2)) {
                  FORT_EWFLUXCOPY2D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
                  FORT_NSFLUXCOPY2D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     flux->getPointer(1),
                     coarse_fine_flux->getPointer(1, side),
                     side);
               } else if (d_dim == tbox::Dimension(3)) {
                  FORT_EWFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
                  FORT_NSFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(1),
                     coarse_fine_flux->getPointer(1, side),
                     side);
                  FORT_TBFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(2),
                     coarse_fine_flux->getPointer(2, side),
                     side);
               }

            }
         }
      } // end first pass over patches in level

      /*
       * Now that fluxes have been computed in each patch on the level,
       *  those that lie at coarse/fine interfaces must be replaced by
       *  the sum of fluxes on coincident faces from the next finer level.
       */

      if (amr_level < hierarchy->getFinestLevelNumber()) {

         std::shared_ptr<hier::PatchLevel> finer_level(
            hierarchy->getPatchLevel(amr_level + 1));

         d_flux_coarsen_schedule[amr_level]->coarsenData();

         finer_level->deallocatePatchData(d_coarse_fine_flux_id);

      } // end face flux fixup

      /*
       * Now that the right fluxes reside at coarse/fine interfaces,
       *  complete function evaluation by differencing.
       */

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);

         const double* dx = patch_geom->getDx();
         const double* xlo = patch_geom->getXLower();
         const double* xhi = patch_geom->getXUpper();
         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > u_cur(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_current)));
         std::shared_ptr<pdat::CellData<double> > source(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_source_term, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > exponential(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_exponential_term, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > flux(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_flux_id)));
         std::shared_ptr<pdat::CellData<double> > fcur(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               f->getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > xdat(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               x->getComponentPatchData(0, *patch)));
         TBOX_ASSERT(u);
         TBOX_ASSERT(u_cur);
         TBOX_ASSERT(source);
         TBOX_ASSERT(exponential);
         TBOX_ASSERT(flux);
         TBOX_ASSERT(fcur);
         TBOX_ASSERT(xdat);
         if (d_dim == tbox::Dimension(1)) {
            FORT_EVALEXPONENTIAL1D(ifirst(0), ilast(0),
               xdat->getPointer(),
               d_lambda,
               exponential->getPointer());
            FORT_EVALSOURCE1D(ifirst(0), ilast(0),
               d_lambda,
               xlo, xhi, dx, d_new_time,
               source->getPointer());
            FORT_EVALBRATU1D(ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               flux->getPointer(0),
               source->getPointer(),
               exponential->getPointer(),
               u_cur->getPointer(),
               u->getPointer(),
               dx, d_current_dt,
               fcur->getPointer());
         } else if (d_dim == tbox::Dimension(2)) {
            FORT_EVALEXPONENTIAL2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               xdat->getPointer(),
               d_lambda,
               exponential->getPointer());
            FORT_EVALSOURCE2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               d_lambda,
               xlo, xhi, dx, d_new_time,
               source->getPointer());
            FORT_EVALBRATU2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               flux->getPointer(0),
               flux->getPointer(1),
               source->getPointer(),
               exponential->getPointer(),
               u_cur->getPointer(),
               u->getPointer(),
               dx, d_current_dt,
               fcur->getPointer());
         } else if (d_dim == tbox::Dimension(3)) {
            FORT_EVALEXPONENTIAL3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               xdat->getPointer(),
               d_lambda,
               exponential->getPointer());
            FORT_EVALSOURCE3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               d_lambda,
               xlo, xhi, dx, d_new_time,
               source->getPointer());
            FORT_EVALBRATU3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               flux->getPointer(0),
               flux->getPointer(1),
               flux->getPointer(2),
               source->getPointer(),
               exponential->getPointer(),
               u_cur->getPointer(),
               u->getPointer(),
               dx, d_current_dt,
               fcur->getPointer());
         }
      } // end second pass over patches in level

      level->deallocatePatchData(d_flux_id);

   } // end loop over levels

   /*
    * Create a coarsen algorithm to coarsen data on fine patch interiors.
    * Then, cycle through the levels and coarsen the data.
    */

   CoarsenAlgorithm f_average(d_dim);
   f_average.registerCoarsen(f->getComponentDescriptorIndex(0),
      f->getComponentDescriptorIndex(0),
      d_soln_coarsen_op);
   for (amr_level = hierarchy->getFinestLevelNumber() - 1;
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      std::shared_ptr<hier::PatchLevel> finer_level(
         hierarchy->getPatchLevel(amr_level + 1));
      /*
       * We take advantage of the knowlege that f has the same
       * structure as the solution and reset the solution coarsening
       * schedule instead of creating a new schedule for f.
       */
      f_average.resetSchedule(d_soln_coarsen_schedule[amr_level]);
      d_soln_coarsen_schedule[amr_level]->coarsenData();
      d_soln_coarsen.resetSchedule(d_soln_coarsen_schedule[amr_level]);
   }
}

/*
 *************************************************************************
 *
 * Evaluate Jacobian-vector product on input vector "v" and place the
 * result in the vector "Jv".
 *
 *************************************************************************
 */

int
ModifiedBratuProblem::jacobianTimesVector(
   std::shared_ptr<solv::SAMRAIVectorReal<double> > v,
   std::shared_ptr<solv::SAMRAIVectorReal<double> > Jv)
{
   TBOX_ASSERT(v);
   TBOX_ASSERT(Jv);

   /*
    * Create a coarsen algorithm to average down fine-to-coarse.
    */

   CoarsenAlgorithm jacv_average(d_dim);
   jacv_average.registerCoarsen(v->getComponentDescriptorIndex(0),
      v->getComponentDescriptorIndex(0),
      d_soln_coarsen_op);

   /*
    * Now cycle through the levels, finest to coarsest, and average
    * down the data.
    */

   std::shared_ptr<hier::PatchHierarchy> hierarchy(Jv->getPatchHierarchy());

   for (int amr_level = hierarchy->getFinestLevelNumber() - 1;
        amr_level >= 0;
        --amr_level) {

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      std::shared_ptr<hier::PatchLevel> finer_level(
         hierarchy->getPatchLevel(amr_level + 1));
      jacv_average.resetSchedule(d_soln_coarsen_schedule[amr_level]);
      d_soln_coarsen_schedule[amr_level]->coarsenData();
      d_soln_coarsen.resetSchedule(d_soln_coarsen_schedule[amr_level]);

   }

   /*
    * Create a refine algorithm needed to fill ghost cells of vector being
    * multiplied.
    */
   s_copy_timer->start();
   RefineAlgorithm jacv_fill;
   jacv_fill.registerRefine(d_soln_scratch_id,
      v->getComponentDescriptorIndex(0),
      d_soln_scratch_id,
      d_soln_refine_op);
   s_copy_timer->stop();

   /*
    * Perform Jacobian-vector product by looping over levels in the hierarchy.
    * The RefineSchedule created by the above RefineAlgorithm will copy
    * the data in the input vector into the interior of the specified scratch
    * space as well as filling the scratch space ghost cells.  Within each
    * level storage locations for each patch are extracted and passed to
    * computation routines.
    *
    * The product is implemented in two passes.  In the first pass, fluxes
    * are computed.  On all levels except the finest, the fluxes that
    * reside at coarse-fine interfaces are replaced by an average of the
    * fluxes from the next finer level.  This ensures proper flux matching
    * at the interfaces.
    *
    * In the second pass, the fluxes that reside at the interfaces are
    * differenced to complete the Jacobian-vector product.  The exponential
    * term, which is a cell-centered quantity, is also added in.
    */

   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));

      /*
       * Allocate scratch space on this level used to store the fluxes.  Also,
       * on all but the coarsest level, allocate space for copies of these
       * fluxes to be averaged down to the next coarser level.
       */

      level->allocatePatchData(d_flux_id);

      if (amr_level > 0) {
         level->allocatePatchData(d_coarse_fine_flux_id);
      }

      /*
       * Fill ghost cell locations on this level.
       */

      s_copy_timer->start();
      jacv_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      d_soln_fill_schedule[amr_level]->fillData(d_new_time);
      d_soln_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      s_copy_timer->stop();

      /*
       * Now sweep through patches on a level, evaluating fluxes on faces
       * of cells.
       */

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();

         std::shared_ptr<pdat::CellData<double> > vdat(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_soln_scratch_id)));
         std::shared_ptr<pdat::SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_diffusion_coef, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > flux(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_flux_id)));
         std::shared_ptr<pdat::OutersideData<double> > coarse_fine_flux(
            SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
               patch->getPatchData(d_coarse_fine_flux_id)));
         TBOX_ASSERT(vdat);
         TBOX_ASSERT(diffusion);
         TBOX_ASSERT(flux);

         if (d_dim == tbox::Dimension(1)) {
            FORT_EVALFACEFLUXES1D(ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               vdat->getPointer(),
               dx,
               flux->getPointer(0));
         } else if (d_dim == tbox::Dimension(2)) {
            FORT_EVALFACEFLUXES2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               vdat->getPointer(),
               dx,
               flux->getPointer(0),
               flux->getPointer(1));
         } else if (d_dim == tbox::Dimension(3)) {
            FORT_EVALFACEFLUXES3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               diffusion->getPointer(2),
               vdat->getPointer(),
               dx,
               flux->getPointer(0),
               flux->getPointer(1),
               flux->getPointer(2));
         }

         /*
          * We employ a convention in which boundary conditions are located on
          * the physical boundaries, but are stored at the centers of ghost
          * cells.  As a result the grid spacing is smaller along the physical
          * boundaries, and the fluxes at these locations have to be modified.
          */

         const std::vector<hier::BoundaryBox>& bdry_faces =
            patch_geom->getCodimensionBoundaries(1);
#if 0
         if (d_dim == tbox::Dimension(1)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getNodeBoundaries();
         } else if (d_dim == tbox::Dimension(2)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getEdgeBoundaries();
         } else if (d_dim == tbox::Dimension(3)) {
            const std::vector<hier::BoundaryBox>& bdry_faces =
               patch_geom->getFaceBoundaries();
         }
#endif
         for (int i = 0; i < static_cast<int>(bdry_faces.size()); ++i) {

            hier::Box bbox = bdry_faces[i].getBox();
            const hier::Index ibeg = bbox.lower();
            const hier::Index iend = bbox.upper();
            int face = bdry_faces[i].getLocationIndex();

            if (d_dim == tbox::Dimension(1)) {
               FORT_EWBCFLUXFIX1D(ifirst(0), ilast(0),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            } else if (d_dim == tbox::Dimension(2)) {
               FORT_EWBCFLUXFIX2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_NSBCFLUXFIX2D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(1),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            } else if (d_dim == tbox::Dimension(3)) {
               FORT_EWBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(0),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_NSBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(1),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
               FORT_TBBCFLUXFIX3D(ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  NUM_GHOSTS_U,
                  dx,
                  vdat->getPointer(),
                  flux->getPointer(2),
                  &bbox.lower()[0], &bbox.upper()[0],
                  face);
            }
         } // end loop over boundary faces

         //         correctPatchFlux( level, patch, vdat );

         /*
          * On all but the coarsest level, we save the face fluxes on the outer edge
          * of the grid level so that they can be averaged down to the next
          * coarser level.
          */

         if (amr_level > 0) {
            TBOX_ASSERT(coarse_fine_flux);
            for (int side = 0; side < 2; ++side) {
               if (d_dim == tbox::Dimension(1)) {
                  FORT_EWFLUXCOPY1D(ifirst(0), ilast(0),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
               } else if (d_dim == tbox::Dimension(2)) {
                  FORT_EWFLUXCOPY2D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
                  FORT_NSFLUXCOPY2D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     flux->getPointer(1),
                     coarse_fine_flux->getPointer(1, side),
                     side);
               } else if (d_dim == tbox::Dimension(3)) {
                  FORT_EWFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(0),
                     coarse_fine_flux->getPointer(0, side),
                     side);
                  FORT_NSFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(1),
                     coarse_fine_flux->getPointer(1, side),
                     side);
                  FORT_TBFLUXCOPY3D(ifirst(0), ilast(0),
                     ifirst(1), ilast(1),
                     ifirst(2), ilast(2),
                     flux->getPointer(2),
                     coarse_fine_flux->getPointer(2, side),
                     side);
               }
            }
         }

      } // end first pass over patches in level

      /*
       * Now that fluxes have been computed in each patch on the level,
       *  those that lie at coarse/fine interfaces must be replaced by
       *  the sum of fluxes on coincident faces from the next finer level.
       */

      if (amr_level < hierarchy->getFinestLevelNumber()) {

         std::shared_ptr<hier::PatchLevel> finer_level(
            hierarchy->getPatchLevel(amr_level + 1));

         d_flux_coarsen_schedule[amr_level]->coarsenData();

         finer_level->deallocatePatchData(d_coarse_fine_flux_id);

      } // end face flux fixup

      /*
       * Now that the right fluxes reside at coarse/fine interfaces,
       *  complete function evaluation by differencing.
       */

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);

         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         const double* dx = patch_geom->getDx();

         std::shared_ptr<pdat::CellData<double> > jac_a(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_jacobian_a_id)));
         std::shared_ptr<pdat::CellData<double> > Jvdat(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               Jv->getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > vdat(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_soln_scratch_id)));
         std::shared_ptr<pdat::SideData<double> > flux(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_flux_id)));
         TBOX_ASSERT(jac_a);
         TBOX_ASSERT(Jvdat);
         TBOX_ASSERT(vdat);
         TBOX_ASSERT(flux);

         TBOX_ASSERT(vdat->getGhostCellWidth() ==
            hier::IntVector(d_dim, NUM_GHOSTS_U));
         TBOX_ASSERT(Jvdat->getGhostCellWidth() == hier::IntVector(d_dim, 0));
         TBOX_ASSERT(jac_a->getGhostCellWidth() == hier::IntVector(d_dim, 0));

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(compjv2d, COMPJV2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               jac_a->getPointer(),
               flux->getPointer(0),
               flux->getPointer(1),
               vdat->getPointer(),
               dx, d_current_dt,
               Jvdat->getPointer());
#if 0
            FORT_BRATUJV2D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               flux->getPointer(0),
               flux->getPointer(1),
               exponential->getPointer(),
               vdat->getPointer(),
               dx, d_current_dt,
               Jvdat->getPointer());
#endif
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(compjv3d, COMPJV3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               jac_a->getPointer(),
               flux->getPointer(0),
               flux->getPointer(1),
               flux->getPointer(2),
               vdat->getPointer(),
               dx, d_current_dt,
               Jvdat->getPointer());
#if 0
            FORT_BRATUJV3D(ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               flux->getPointer(0),
               flux->getPointer(1),
               flux->getPointer(2),
               exponential->getPointer(),
               vdat->getPointer(),
               dx, d_current_dt,
               Jvdat->getPointer());
#endif
         }

      }   // end second pass over patches

      level->deallocatePatchData(d_flux_id);

   } // end pass over levels

   return 0;

}

/*
 *************************************************************************
 *
 * Set up FAC hierarchy preconditioner for Jacobian system.  Here we
 * use the FAC hierarchy solver in SAMRAI which automatically sets
 * up the composite grid system and uses hypre as a solver on each
 * level.  Fortunately, the FAC hierarchy solver suits the problem
 * we are solving, as the Jacobian matrix arises by discretizing:
 *
 *    div( D(x,t)*grad(u) ) + 1/dt - lambda * exp(u) - df/du
 *
 * evaluated at some iterate of the solution u.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::setupBratuPreconditioner(
   std::shared_ptr<solv::SAMRAIVectorReal<double> > x)
{
   TBOX_ASSERT(x);

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   RefineAlgorithm eval_fill;
   eval_fill.registerRefine(d_soln_scratch_id,
      x->getComponentDescriptorIndex(0),
      d_soln_scratch_id,
      d_soln_refine_op);

   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));

      eval_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      d_soln_fill_schedule[amr_level]->fillData(d_new_time);
      d_soln_fill.resetSchedule(d_soln_fill_schedule[amr_level]);

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);

         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         const double* dx = patch_geom->getDx();
         const double* xlo = patch_geom->getXLower();

         double cell_vol = dx[0];
         if (d_dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
         }
         if (d_dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
         }
         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > exponential(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_exponential_term, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > source(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_source_term, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_diffusion_coef, d_scratch)));

         std::shared_ptr<pdat::CellData<double> > a(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_precond_a_id)));
         std::shared_ptr<pdat::FaceData<double> > b(
            SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
               patch->getPatchData(d_precond_b_id)));

         TBOX_ASSERT(exponential);
         TBOX_ASSERT(source);
         TBOX_ASSERT(u);
         TBOX_ASSERT(diffusion);
         TBOX_ASSERT(a);
         TBOX_ASSERT(b);

         /*
          * Compute exponential = lambda * exp(u)
          */
         if (d_dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(compexpu1d, COMPEXPU1D) (ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv1d, COMPSRCDERV1D) (ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         } else if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(compexpu2d, COMPEXPU2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv2d, COMPSRCDERV2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(compexpu3d, COMPEXPU3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv3d, COMPSRCDERV3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         }
         /*
          * Compute a = cell_vol * ( 1 - d_current_dt*(exponential + source))
          */
         if (d_dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(compfacdiag1d, COMPFACDIAG1D) (ifirst(0), ilast(0),
               d_current_dt,
               cell_vol,
               exponential->getPointer(),
               source->getPointer(),
               a->getPointer());
         } else if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(compfacdiag2d, COMPFACDIAG2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               d_current_dt,
               cell_vol,
               exponential->getPointer(),
               source->getPointer(),
               a->getPointer());
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(compfacdiag3d, COMPFACDIAG3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               d_current_dt,
               cell_vol,
               exponential->getPointer(),
               source->getPointer(),
               a->getPointer());
         }

         /*
          * Copy side-centered diffusion coefficients to face-centered
          * and multiply by cell_vol*d_current_dt.
          */
         if (d_dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(compfacoffdiag1d, COMPFACOFFDIAG1D) (ifirst(0), ilast(0),
               d_current_dt,
               cell_vol,
               diffusion->getPointer(0),
               b->getPointer(0));
         } else if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(compfacoffdiag2d, COMPFACOFFDIAG2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               d_current_dt,
               cell_vol,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               b->getPointer(0),
               b->getPointer(1));
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(compfacoffdiag3d, COMPFACOFFDIAG3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               d_current_dt,
               cell_vol,
               diffusion->getPointer(0),
               diffusion->getPointer(1),
               diffusion->getPointer(2),
               b->getPointer(0),
               b->getPointer(1),
               b->getPointer(2));
         }

      }

   }

   // This will have to generalized if we add in other boundary conditions

   d_FAC_solver->setBoundaries("Dirichlet");

   if (d_precond_b_id == -1) {
      d_FAC_solver->setDConstant(1.0);
   } else {
      d_FAC_solver->setDConstant(1.0);
   }
   if (d_precond_a_id == -1) {
      d_FAC_solver->setCConstant(0.0);
   } else {
      d_FAC_solver->setCPatchDataId(d_precond_a_id);
   }

}

/*
 *************************************************************************
 *
 * Apply preconditioner where right-hand-side is "r" and "z" is the
 * solution.   This routine assumes that the preconditioner setup call
 * has already been invoked.  Return 0 if preconditioner fails;
 * return 1 otherwise.
 *
 *************************************************************************
 */

int ModifiedBratuProblem::applyBratuPreconditioner(
   std::shared_ptr<solv::SAMRAIVectorReal<double> > r,
   std::shared_ptr<solv::SAMRAIVectorReal<double> > z)
{
   TBOX_ASSERT(r);
   TBOX_ASSERT(z);

   int ret_val = 0;

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   int r_indx = r->getComponentDescriptorIndex(0);
   int z_indx = z->getComponentDescriptorIndex(0);

   /*
    * Create a coarsen algorithm to coarsen rhs data on fine patch interiors.
    * Then, cycle through the levels and coarsen the data.
    */
   CoarsenAlgorithm pc_rhs_average(d_dim);
   pc_rhs_average.registerCoarsen(r_indx, r_indx, d_soln_coarsen_op);

   for (int amr_level = hierarchy->getFinestLevelNumber() - 1;
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      std::shared_ptr<hier::PatchLevel> finer_level(
         hierarchy->getPatchLevel(amr_level + 1));
      pc_rhs_average.resetSchedule(d_soln_coarsen_schedule[amr_level]);
      d_soln_coarsen_schedule[amr_level]->coarsenData();
      d_soln_coarsen.resetSchedule(d_soln_coarsen_schedule[amr_level]);
   }

   math::HierarchyCellDataOpsReal<double>
   math_ops(hierarchy,
            0,
            hierarchy->getFinestLevelNumber());
   math_ops.setToScalar(d_soln_scratch_id,
      0.0);

   s_pc_timer->start();

   d_FAC_solver->solveSystem(d_soln_scratch_id,
      r_indx,
      hierarchy,
      0,
      hierarchy->getFinestLevelNumber());

   s_pc_timer->stop();

   /*
    * Create a coarsen algorithm to coarsen soln data on fine patch interiors.
    * Then, cycle through the levels and coarsen the data.
    */
   CoarsenAlgorithm pc_sol_average(d_dim);
   pc_sol_average.registerCoarsen(d_soln_scratch_id,
      d_soln_scratch_id,
      d_soln_coarsen_op);

   for (int amr_level = hierarchy->getFinestLevelNumber() - 1;
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));
      std::shared_ptr<hier::PatchLevel> finer_level(
         hierarchy->getPatchLevel(amr_level + 1));
      pc_sol_average.resetSchedule(d_scratch_soln_coarsen_schedule[amr_level]);
      d_scratch_soln_coarsen_schedule[amr_level]->coarsenData();
      d_scratch_soln_coarsen.resetSchedule(d_scratch_soln_coarsen_schedule[
            amr_level]);
   }

   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         std::shared_ptr<pdat::CellData<double> > src_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_soln_scratch_id)));
         std::shared_ptr<pdat::CellData<double> > dst_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(z_indx)));
         TBOX_ASSERT(src_data);
         TBOX_ASSERT(dst_data);

         dst_data->copy(*src_data);
      }
   }

   return ret_val;
}

/*
 *************************************************************************
 *
 * Evaluate the Jacobian matrix A(x) for a given solution x.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::evaluateBratuJacobian(
   std::shared_ptr<solv::SAMRAIVectorReal<double> > x)
{
   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      d_solution_vector->getPatchHierarchy());

   RefineAlgorithm eval_fill;
   eval_fill.registerRefine(d_soln_scratch_id,
      x->getComponentDescriptorIndex(0),
      d_soln_scratch_id,
      d_soln_refine_op);

   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(
                                                   amr_level));

      eval_fill.resetSchedule(d_soln_fill_schedule[amr_level]);
      d_soln_fill_schedule[amr_level]->fillData(d_new_time);
      d_soln_fill.resetSchedule(d_soln_fill_schedule[amr_level]);

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);

         const hier::Index ifirst = patch->getBox().lower();
         const hier::Index ilast = patch->getBox().upper();

         const double* dx = patch_geom->getDx();
         const double* xlo = patch_geom->getXLower();

         double cell_vol = dx[0];
         if (d_dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
         }
         if (d_dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
         }

         std::shared_ptr<pdat::CellData<double> > u(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_solution, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > exponential(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_exponential_term, d_scratch)));
         std::shared_ptr<pdat::CellData<double> > source(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_source_term, d_scratch)));
         std::shared_ptr<pdat::SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(d_diffusion_coef, d_scratch)));

         std::shared_ptr<pdat::CellData<double> > a(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_jacobian_a_id)));
         std::shared_ptr<pdat::FaceData<double> > b(
            SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
               patch->getPatchData(d_jacobian_b_id)));
         TBOX_ASSERT(u);
         TBOX_ASSERT(exponential);
         TBOX_ASSERT(source);
         TBOX_ASSERT(diffusion);
         TBOX_ASSERT(a);
         TBOX_ASSERT(b);

         if (d_dim == tbox::Dimension(1)) {
            /*
             * Compute exponential = lambda * exp(u)
             */
            SAMRAI_F77_FUNC(compexpu1d, COMPEXPU1D) (ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv1d, COMPSRCDERV1D) (ifirst(0), ilast(0),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         } else if (d_dim == tbox::Dimension(2)) {
            /*
             * Compute exponential = lambda * exp(u)
             */
            SAMRAI_F77_FUNC(compexpu2d, COMPEXPU2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv2d, COMPSRCDERV2D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         } else if (d_dim == tbox::Dimension(3)) {
            /*
             * Compute exponential = lambda * exp(u)
             */
            SAMRAI_F77_FUNC(compexpu3d, COMPEXPU3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               d_lambda,
               u->getPointer(),
               exponential->getPointer());

            /*
             * Compute source = df/du
             */
            SAMRAI_F77_FUNC(compsrcderv3d, COMPSRCDERV3D) (ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               NUM_GHOSTS_U,
               xlo, dx,
               d_new_time,
               d_lambda,
               u->getPointer(),
               source->getPointer());
         }
         /*
          * Compute a = lambda*exp(u) + df/du
          */
         math::PatchCellDataOpsReal<double> cell_basic_ops;
         cell_basic_ops.add(a,
            exponential,
            source,
            a->getBox());

      }
   }
}

/*
 *************************************************************************
 *
 * Set physical boundary conditions for patch at a given time.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(time);
   NULL_USE(ghost_width_to_fill);

   /*
    * Grab data to operate on.
    */

   std::shared_ptr<pdat::CellData<double> > u(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_soln_scratch_id)));
   TBOX_ASSERT(u);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   /*
    * Determine boxes that touch the physical boundary.
    */

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const std::vector<hier::BoundaryBox>& boundary =
      patch_geom->getCodimensionBoundaries(1);
#if 0
   if (d_dim == tbox::Dimension(1)) {
      const std::vector<hier::BoundaryBox>& boundary =
         patch_geom->getNodeBoundaries();
   } else if (d_dim == tbox::Dimension(2)) {
      const std::vector<hier::BoundaryBox>& boundary =
         patch_geom->getEdgeBoundaries();
   } else if (d_dim == tbox::Dimension(3)) {
      const std::vector<hier::BoundaryBox>& boundary =
         patch_geom->getFaceBoundaries();
   }
#endif

   /*
    * Walk the list of boxes that describe the boundary, and apply
    * boundary conditions to each segment.
    */

   int face;
   for (int i = 0; i < static_cast<int>(boundary.size()); ++i) {
      hier::Box bbox = hier::Box(boundary[i].getBox());
      face = boundary[i].getLocationIndex();
      if (d_dim == tbox::Dimension(1)) {
         FORT_SETBC1D(ifirst(0), ilast(0),
            NUM_GHOSTS_U,
            u->getPointer(),
            &bbox.lower()[0], &bbox.upper()[0],
            face);
      } else if (d_dim == tbox::Dimension(2)) {
         FORT_SETBC2D(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            NUM_GHOSTS_U,
            u->getPointer(),
            &bbox.lower()[0], &bbox.upper()[0],
            face);
      } else if (d_dim == tbox::Dimension(3)) {
         FORT_SETBC3D(ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            ifirst(2), ilast(2),
            NUM_GHOSTS_U,
            u->getPointer(),
            &bbox.lower()[0], &bbox.upper()[0],
            face);
      }
   }
}

/*
 *************************************************************************
 *
 * Read data from input database.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   if (input_db->keyExists("lambda")) {
      d_lambda = input_db->getDouble("lambda");
   } else {
      TBOX_ERROR(d_object_name << " -- Key data `lambda'"
                               << " missing in input.");
   }

   if (input_db->keyExists("timestep")) {
      d_input_dt = input_db->getDouble("timestep");
      if (d_input_dt < 0.0) {
         TBOX_ERROR(d_object_name << " Input error: timestep < 0.0");
      }
   } else {
      if (!is_from_restart) {
         TBOX_ERROR(d_object_name << " -- Key data `timestep'"
                                  << " missing in input.");
      }
   }

   if (input_db->keyExists("max_precond_its")) {
      d_max_precond_its = input_db->getInteger("max_precond_its");
      if (d_max_precond_its < 0) {
         TBOX_ERROR(d_object_name << " Input error: max_precond_its < 0");
      }
   } else {
      if (!is_from_restart) {
         TBOX_ERROR(d_object_name << " -- Key data `max_precond_its'"
                                  << " missing in input.");
      }
   }

   if (input_db->keyExists("precond_tol")) {
      d_precond_tol = input_db->getDouble("precond_tol");
      if (d_precond_tol <= 0.0) {
         TBOX_ERROR(d_object_name << " Input error: precond_tol <= 0.0");
      }
   } else {
      if (!is_from_restart) {
         TBOX_ERROR(d_object_name << " -- Key data `precond_tol'"
                                  << " missing in input.");
      }
   }

}

/*
 *************************************************************************
 *
 * Put class version number and data members in restart database.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("MODIFIED_BRATU_PROBLEM", MODIFIED_BRATU_PROBLEM);

   restart_db->putDouble("d_lambda", d_lambda);
   restart_db->putDouble("d_input_dt", d_input_dt);
   restart_db->putInteger("d_max_precond_its", d_max_precond_its);
   restart_db->putDouble("d_precond_tol", d_precond_tol);

}

/*
 *************************************************************************
 *
 * Print all class data members to given output stream.
 *
 *************************************************************************
 */

void ModifiedBratuProblem::printClassData(
   std::ostream& os) const
{
   os << "\nModifiedBratuProblem::printClassData..." << std::endl;
   os << "ModifiedBratuProblem: this = " << (ModifiedBratuProblem *)this
      << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_grid_geometry = "
      << d_grid_geometry.get() << std::endl;

   os << "d_soln_scratch_id =   " << d_soln_scratch_id << std::endl;
   os << "d_flux_id =   " << d_flux_id << std::endl;
   os << "d_coarse_fine_flux_id =   " << d_coarse_fine_flux_id << std::endl;
   os << "d_jacobian_a_id =   " << d_jacobian_a_id << std::endl;
   os << "d_jacobian_b_id =   " << d_jacobian_b_id << std::endl;
   os << "d_weight_id =   " << d_weight_id << std::endl;
   os << "d_nghosts =   " << d_nghosts << std::endl;

   os << "d_lambda =   " << d_lambda << std::endl;
   os << "d_input_dt = " << d_input_dt << std::endl;

   os << "d_current_time =   " << d_current_time << std::endl;
   os << "d_new_time = " << d_new_time << std::endl;
   os << "d_current_dt =   " << d_current_dt << std::endl;

   os << "d_max_precond_its =   " << d_max_precond_its << std::endl;
   os << "d_precond_tol = " << d_precond_tol << std::endl;

}

void ModifiedBratuProblem::getLevelEdges(
   hier::BoxContainer& boxes,
   std::shared_ptr<hier::Patch> patch,
   std::shared_ptr<hier::PatchLevel> level,
   const tbox::Dimension::dir_t dim,
   const int face)
{

   /*
    * Shift the box associated with the patch on the side indicated by
    * (dim,face).
    */

   boxes.clear();
   hier::Box box = patch->getBox();
   hier::Box boundary = box;
   if (face == 0) {
      boundary.setLower(dim, box.lower(dim) - 1);
      boundary.setUpper(dim, box.lower(dim) - 1);
   } else {
      boundary.setLower(dim, box.upper(dim) + 1);
      boundary.setUpper(dim, box.upper(dim) + 1);
   }

   /*
    * Initialize the list of boxes that constitute the edge of the
    * level, and delete the portion of the box that lies interior to
    * the level.
    */

   boxes.pushBack(boundary);
   boxes.removeIntersections(level->getBoxes());

   /*
    * Finally delete the part of the level edge that meets the
    * physical boundary.
    */

   std::shared_ptr<geom::CartesianPatchGeometry> geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch->getPatchGeometry()));
   TBOX_ASSERT(geometry);
   const std::vector<hier::BoundaryBox>& boundary_boxes =
      geometry->getCodimensionBoundaries(1);
   for (int i = 0; i < static_cast<int>(boundary_boxes.size()); ++i) {
      boxes.removeIntersections(boundary_boxes[i].getBox());
   }

   /*
    * Put the result into canonical order.
    */

   boxes.simplify();
}

void ModifiedBratuProblem::correctLevelFlux(
   std::shared_ptr<hier::PatchLevel> level)
{
   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;
      const hier::Box box = patch->getBox();
      std::shared_ptr<pdat::SideData<double> > flux_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch->getPatchData(d_flux_id)));
      TBOX_ASSERT(flux_data);

      /*
       * For each direction, for each side:  compute the index space that
       * describes the part of the patch that lies on the outer edge of
       * the refinement level.  Then, scale the fluxes along that edge
       * to account for the different grid spacing when CONSTANT_REFINE
       * is used to fill the ghost cells.
       */

      for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
         for (int s = 0; s <= 1; ++s) {
            hier::Index delta(d_dim, 0);
            delta(d) = ((s == 0) ? 1 : -1);
            hier::Index twodelta(d_dim, 0);
            twodelta(d) = ((s == 0) ? 2 : -2);
            hier::BoxContainer level_edges;
            getLevelEdges(level_edges, patch, level, d, s);
            for (hier::BoxContainer::iterator l = level_edges.begin();
                 l != level_edges.end(); ++l) {
               pdat::CellIterator icend(pdat::CellGeometry::end(*l));
               for (pdat::CellIterator ic(pdat::CellGeometry::begin(*l));
                    ic != icend; ++ic) {
                  pdat::SideIndex iside(*ic + delta, d, s);
                  (*flux_data)(iside) = 2.0 * (*flux_data)(iside) / 3.0;
               }   // cell loop
            }   // box loop
         }   // side loop
      }   // direction loop
   }   // patch loop
}

void ModifiedBratuProblem::correctPatchFlux(
   std::shared_ptr<hier::PatchLevel> level,
   std::shared_ptr<hier::Patch> patch,
   std::shared_ptr<pdat::CellData<double> > u)
{
   const hier::Box box = patch->getBox();
   std::shared_ptr<pdat::SideData<double> > flux_data(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
         patch->getPatchData(d_flux_id)));
   const std::shared_ptr<geom::CartesianPatchGeometry> geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch->getPatchGeometry()));
   TBOX_ASSERT(flux_data);
   TBOX_ASSERT(geometry);
   const double* dx = geometry->getDx();

   for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
      for (int s = 0; s <= 1; ++s) {
         hier::Index delta1(d_dim, 0);
         delta1(d) = ((s == 0) ? 1 : -1);
         hier::Index delta2(d_dim, 0);
         delta2(d) = ((s == 0) ? 2 : -2);
         double factor = ((s == 0) ? 1.0 / dx[d] : -1.0 / dx[d]);
         hier::BoxContainer level_edges;
         getLevelEdges(level_edges, patch, level, d, s);
         for (hier::BoxContainer::iterator l = level_edges.begin();
              l != level_edges.end(); ++l) {
            pdat::CellIterator icend(pdat::CellGeometry::end(*l));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(*l));
                 ic != icend; ++ic) {
               pdat::SideIndex iside(*ic + delta1, d, s);
               (*flux_data)(iside) = factor * (-8.0 * (*u)(*ic) / 15.0
                                               + (*u)(*ic + delta1) / 3.0
                                               + (*u)(*ic + delta2) / 5.0);
            } // cell loop
         } // box loop
      } // side loop
   } // direction loop
}

#endif
