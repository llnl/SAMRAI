/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic method-of-lines time integration algorithm
 *
 ************************************************************************/
#include "SAMRAI/algs/MethodOfLinesIntegrator.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstdlib>
#include <fstream>

namespace SAMRAI {
namespace algs {

const int MethodOfLinesIntegrator::ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION = 2;

/*
 *************************************************************************
 *
 * The constructor and destructor for MethodOfLinesIntegrator.
 *
 *************************************************************************
 */

MethodOfLinesIntegrator::MethodOfLinesIntegrator(
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db,
   MethodOfLinesPatchStrategy* patch_strategy):
   d_object_name(object_name),
   d_order(3),
   d_patch_strategy(patch_strategy),
   d_current(hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
   d_scratch(hier::VariableDatabase::getDatabase()->getContext("SCRATCH"))
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(patch_strategy != 0);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
      this);

   /*
    * hier::Variable contexts used in algorithm.
    */
   d_patch_strategy->setInteriorContext(d_current);
   d_patch_strategy->setInteriorWithGhostsContext(d_scratch);

   /*
    * Set default to third-order SSP Runge-Kutta method.
    */
   d_alpha_1.resize(d_order);
   d_alpha_1[0] = 1.0;
   d_alpha_1[1] = 0.75;
   d_alpha_1[2] = 1.0 / 3.0;
   d_alpha_2.resize(d_order);
   d_alpha_2[0] = 0.0;
   d_alpha_2[1] = 0.25;
   d_alpha_2[2] = 2.0 / 3.0;
   d_beta.resize(d_order);
   d_beta[0] = 1.0;
   d_beta[1] = 0.25;
   d_beta[2] = 2.0 / 3.0;

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

}

/*
 *************************************************************************
 *
 * Destructor tells tbox::RestartManager to remove this object from the
 * list of restart items.
 *
 *************************************************************************
 */

MethodOfLinesIntegrator::~MethodOfLinesIntegrator()
{
   tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/*
 *************************************************************************
 *
 *
 * Initialize integrator by:
 *
 *
 *
 *   (1) Setting the number of time data levels based on needs of
 *
 *       the gridding algorithm
 *
 *   (2) Invoking variable registration in patch strategy.
 *
 *
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::initializeIntegrator(
   const std::shared_ptr<mesh::GriddingAlgorithm>& gridding_alg)
{
   NULL_USE(gridding_alg);
   TBOX_ASSERT(gridding_alg);

   /*
    * We may eventually need support for three (or more) time
    * levels, information which may be accessed from the
    * gridding algorithm.  Since we don't yet support this,
    * this method simply registers variables with the integrator.
    */

   /*
    * Call variable registration in patch strategy.
    */
   d_patch_strategy->registerModelVariables(this);
}
/*
 *************************************************************************
 *
 * Calculate the stable time increment by taking the minimum over
 * all patches on all levels in the hierarchy.
 *
 *************************************************************************
 */

double
MethodOfLinesIntegrator::getTimestep(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const double time) const
{
   TBOX_ASSERT(hierarchy);

   double dt = tbox::MathUtilities<double>::getMax();
   const int nlevels = hierarchy->getNumberOfLevels();

   for (int l = 0; l < nlevels; ++l) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(l);

      TBOX_ASSERT(level);

      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {

         const std::shared_ptr<hier::Patch>& patch = *p;

         const double dt_patch =
            d_patch_strategy->computeStableDtOnPatch(*patch, time);
         if (dt_patch < dt) {
            dt = dt_patch;
         }
      }
   }

   const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&dt, 1, MPI_MIN);
   }

   return dt;
}

/*
 *************************************************************************
 *
 * Advance the solution through the specified time increment using the
 * general RK algorithm.  Each of the following steps is performed over
 * all hierarchy levels.
 *
 * (1) Copy solution values from current context to scratch context.
 *
 * (2) RK multistep loop for d(U)/dt = F(U):
 *
 *    do i = 1, order
 *       U_i = U_n + alpha_i * dt/(order) * F(U_i)
 *    end do
 *
 * (3) Copy last update of scratch solution to current context.
 *
 * Note that each update is performed by the concrete patch strategy
 * in which the numerical routines are defined.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::advanceHierarchy(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const double time,
   const double dt)
{
   TBOX_ASSERT(hierarchy);

   /*
    * Stamp data on all levels to current simulation time.
    */
   const int nlevels = hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < nlevels; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      TBOX_ASSERT(level);

      level->setTime(time, d_current_data);
      level->setTime(time, d_scratch_data);
      level->setTime(time, d_rhs_data);

      /*
       * Allocate memory for U_scratch and rhs data
       */
      level->allocatePatchData(d_scratch_data, time);
      level->allocatePatchData(d_rhs_data, time);

      copyCurrentToScratch(level);
   }

   /*
    * Loop through Runge-Kutta steps
    */
   for (int rkstep = 0; rkstep < d_order; ++rkstep) {


      /*
       * Fill ghost cells for all levels
       */
      for (int ln = 0; ln < nlevels; ++ln) {

         d_bdry_sched_advance[ln]->fillData(time);

      }

      /*
       * Loop through levels in the patch hierarchy and advance data on
       * each level by a single RK step.
       */
      for (int ln = 0; ln < nlevels; ++ln) {

         /*
          * Loop through patches in current level and "singleStep" on each
          * patch.
          */
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));

         TBOX_ASSERT(level);

         for (hier::PatchLevel::iterator p(level->begin());
              p != level->end(); ++p) {

            const std::shared_ptr<hier::Patch>& patch = *p;
            d_patch_strategy->singleStep(*patch,
               dt,
               d_alpha_1[rkstep],
               d_alpha_2[rkstep],
               d_beta[rkstep]);

         }  // patch loop
      }  // levels loop

      /*
       * Coarsen data from finest to coarsest
       */
      for (int ln = nlevels-1; ln > 0; --ln) {
         d_coarsen_schedule[ln]->coarsenData();
      }

   }  // rksteps loop

   for (int ln = 0; ln < nlevels; ++ln) {
      copyScratchToCurrent(hierarchy->getPatchLevel(ln));

      /*
       * update timestamp to time after advance
       */
      hierarchy->getPatchLevel(ln)->setTime(time + dt, d_current_data);
   }

   /*
    * dallocate U_scratch and rhs data
    */
   for (int ln = 0; ln < nlevels; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      level->deallocatePatchData(d_scratch_data);
      level->deallocatePatchData(d_rhs_data);
   }

}

/*
 *************************************************************************
 *
 * Register the variables with the method of lines solution algorithm
 * according to specified algorithm role (i.e., MOL_VAR_TYPE).
 *
 *       du/dt = F(u)
 *           u  - defined as "solution"
 *         F(u) - defined as "right-hand-side"
 *
 * Assignment of descriptor indices to variable lists, component
 * selectors, and communication  algorithms takes place here.
 *
 * The different cases are:
 *
 * SOLN:
 *            One time level of data is maintained for the current
 *            solution and a "scratch" copy is used during the update
 *            process.
 *
 *            Two factories are needed: SCRATCH, CURRENT.
 *
 *            SCRATCH index is added to d_scratch_data.
 *            CURRENT index is added to d_current_data.
 *
 * RHS:
 *            Only one time level of data is stored and no scratch space
 *            is used.  Data may be set and manipulated at will in user
 *            routines.  Data (including ghost values) is never touched
 *            outside of user routines.
 *
 *            One factory needed: CURRENT.
 *
 *            CURRENT index is added to d_current_data.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::registerVariable(
   const std::shared_ptr<hier::Variable>& variable,
   const hier::IntVector& ghosts,
   const MOL_VAR_TYPE m_v_type,
   const std::shared_ptr<hier::BaseGridGeometry>& transfer_geom,
   const std::string& coarsen_name,
   const std::string& refine_name)
{
   TBOX_ASSERT(variable);
   TBOX_ASSERT(transfer_geom);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*variable, ghosts);

   tbox::Dimension dim(ghosts.getDim());

   if (!d_bdry_fill_advance) {
      /*
       * One-time set-up for communication algorithms.
       * We wait until this point to do this because we need a dimension.
       */
      d_bdry_fill_advance.reset(new xfer::RefineAlgorithm());
      d_fill_after_regrid.reset(new xfer::RefineAlgorithm());
      d_fill_before_tagging.reset(new xfer::RefineAlgorithm());
      d_coarsen_algorithm.reset(new xfer::CoarsenAlgorithm(dim));
   }

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   switch (m_v_type) {

      case SOLN: {
         /*
          * Associate the current and scratch contexts with the variable in the
          * database and get the patch data identifiers.  The flag arrays will
          * be used to manage the allocation and deallocation of current and
          * scratch data.
          */
         d_soln_variables.push_back(variable);

         const hier::IntVector no_ghosts(dim, 0);

         const int current = variable_db->registerVariableAndContext(variable,
               d_current,
               no_ghosts);

         const int scratch = variable_db->registerVariableAndContext(variable,
               d_scratch,
               ghosts);

         d_current_data.setFlag(current);

         d_scratch_data.setFlag(scratch);

         /*
          * Register variable and context needed for restart.
          */
         hier::PatchDataRestartManager::getManager()->
         registerPatchDataForRestart(current);

         /*
          * Ask the geometry for the appropriate refinement operator and
          * register that operator and the variables with the communication
          * algorithms.  Two different communication algorithms are required by
          * the RK method.  The Fillghosts algorithm is called during the
          * normal Runge-Kutta time stepping and fills the ghost cells of the
          * scratch variables.  The regrid algorithm is called after regrid and
          * fills the current data on the new level.
          */

         std::shared_ptr<hier::RefineOperator> refine_operator(
            transfer_geom->lookupRefineOperator(variable, refine_name));

         //  Fill ghosts for a variable using always the "scratch" context
         d_bdry_fill_advance->registerRefine(
            scratch,    // destination
            scratch,    // source
            scratch,    // temporary work space
            refine_operator);

         //  After regrid, use "current" context to communicate information
         //  to updated patches.  Use "scratch" as the temporary storage.
         d_fill_after_regrid->registerRefine(
            current,    // destination
            current,    // source
            scratch,    // temporary work space
            refine_operator);

         //  Before tagging error cells, copy data in current context to
         //  scratch context.  Note that this operation is not a simple
         //  copy - it also requires filling of ghost cells.  This is why
         //  it is designated as a refine operation.
         d_fill_before_tagging->registerRefine(
            scratch,    // destination
            current,    // source
            scratch,    // temporary work space
            refine_operator);

         std::shared_ptr<hier::CoarsenOperator> coarsen_operator(
            transfer_geom->lookupCoarsenOperator(variable, coarsen_name));

         //  Coarsen solution between levels during RK process so that
         //  coarser levels see the fine solution during integration.
         d_coarsen_algorithm->registerCoarsen(scratch,
            scratch,
            coarsen_operator);

         break;
      }

      case RHS: {
         /*
          * Associate the current context with the RHS variable in the
          * database.  The d_rhs_data component selector will be used to
          * allocate and de-allocate rhs data.
          * NOTE:  The d_rhs_data component selector was added 3/23/00 to
          * facilitate allocation and de-allocation of rhs data for restarts.
          */
         d_rhs_variables.push_back(variable);

         const int current = variable_db->registerVariableAndContext(variable,
               d_current,
               ghosts);

         d_rhs_data.setFlag(current);

         break;
      }

      default: {

         TBOX_ERROR(d_object_name << ":  "
                                  << "unknown MOL_VAR_TYPE = " << m_v_type
                                  << std::endl);

      }

   }

}

/*
 *************************************************************************
 *
 * Allocate data for new level in hierarchy and initialize that data.
 * If the new level replaces a pre-existing level in the hierarchy,
 * data is copied from that level to the new level on their intersection.
 * Other data on the new level is set by interpolating from coarser
 * levels in the hierarchy.  Then, user-defined initialization routines
 * are called.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::initializeLevelData(
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
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(level_number >= 0);
#ifdef DEBUG_CHECK_ASSERTIONS
   if (old_level) {
      TBOX_ASSERT(level_number == old_level->getLevelNumber());
      TBOX_ASSERT_OBJDIM_EQUALITY2(*hierarchy, *old_level);
   }
#endif

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Allocate storage needed to initialize level and fill data from
    * coarser levels in AMR hierarchy.
    */
   level->allocatePatchData(d_current_data, time);
   level->allocatePatchData(d_scratch_data, time);

   if ((level_number > 0) || old_level) {
      d_fill_after_regrid->createSchedule(
         level,
         old_level,
         level_number - 1,
         hierarchy,
         d_patch_strategy)->fillData(time);
   }

   level->deallocatePatchData(d_scratch_data);

   /*
    * Initialize current data for new level.
    */
   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      d_patch_strategy->initializeDataOnPatch(*patch,
         time,
         initial_time);
   }
}

/*
 *************************************************************************
 *
 * Re-generate communication schedule after changes to the specified
 * range of levels in the hierarchy.
 *
 *************************************************************************
 */
void
MethodOfLinesIntegrator::resetHierarchyConfiguration(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(finest_level);

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level <= finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));

   int finest_hiera_level = hierarchy->getFinestLevelNumber();

   //  If we have added or removed a level, resize the schedule arrays
   d_bdry_sched_advance.resize(finest_hiera_level + 1);
   d_coarsen_schedule.resize(finest_hiera_level + 1);

   //  Build coarsen and refine communication schedules.
   for (int ln = coarsest_level; ln <= finest_hiera_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      TBOX_ASSERT(level);

      d_bdry_sched_advance[ln] =
         d_bdry_fill_advance->createSchedule(
            level,
            ln - 1,
            hierarchy,
            d_patch_strategy);

      // coarsen schedule only for levels > 0
      if (ln > 0) {
         std::shared_ptr<hier::PatchLevel> coarser_level(
            hierarchy->getPatchLevel(ln - 1));
         d_coarsen_schedule[ln] =
            d_coarsen_algorithm->createSchedule(
               coarser_level,
               level,
               0);
      }

   }
}

/*
 *************************************************************************
 *
 * Fill ghost cells for patches on level and call application-specific
 * cell tagging routines.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int ln,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(hierarchy->getPatchLevel(ln));

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));

   level->allocatePatchData(d_scratch_data, time);

   /*
    * Transfer information from the "current" context to the "scratch"
    * context, on the current level. Note that ghosts will be filled
    * in this process.  We create and apply the schedule at the same
    * time because this routine is only called during
    * a regrid step, and the changing grid system means the schedule will
    * change since the last time it was called.
    */
   d_fill_before_tagging->createSchedule(level,
      level,
      ln - 1,
      hierarchy,
      d_patch_strategy)->fillData(time);

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      d_patch_strategy->tagGradientDetectorCells(*patch,
         time,
         initial_time,
         tag_index,
         uses_richardson_extrapolation_too);
   }

   level->deallocatePatchData(d_scratch_data);

}

/*
 *************************************************************************
 *
 * Writes the class version number, order, and
 * alpha array to the restart database.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION",
      ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION);

   restart_db->putDoubleVector("alpha_1", d_alpha_1);
   restart_db->putDoubleVector("alpha_2", d_alpha_2);
   restart_db->putDoubleVector("beta", d_beta);
}

/*
 *************************************************************************
 *
 * Reads in paramemters from the database overriding any values
 * read in from the restart database. Also checks to make sure that
 * number of alpha values specified equals order of Runga-Kutta scheme.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (input_db) {

      bool read_on_restart =
         input_db->getBoolWithDefault("read_on_restart", false);
      if (!is_from_restart || read_on_restart) {

         if (input_db->keyExists("alpha_1")) {
            size_t array_size = input_db->getArraySize("alpha_1");
            if (array_size > 3) {
               TBOX_ERROR("MethodOfLinesIntegrator::getFromInput() error...\n"
                  << "number of alpha_1 entries must be <=3." << std::endl);
            }
            d_alpha_1 = input_db->getDoubleVector("alpha_1");
         }

         if (input_db->keyExists("alpha_2")) {
            size_t array_size = input_db->getArraySize("alpha_2");
            if (array_size > 3) {
               TBOX_ERROR("MethodOfLinesIntegrator::getFromInput() error...\n"
                  << "number of alpha_2 entries must be <=3." << std::endl);
            }
            d_alpha_2 = input_db->getDoubleVector("alpha_2");
         }

         if (input_db->keyExists("beta")) {
            size_t array_size = input_db->getArraySize("beta");
            if (array_size > 3) {
               TBOX_ERROR("MethodOfLinesIntegrator::getFromInput() error...\n"
                  << "number of beta entries must be <=3." << std::endl);
            }
            d_beta = input_db->getDoubleVector("beta");
         }

         if (d_alpha_1.size() != d_alpha_2.size() ||
             d_alpha_2.size() != d_beta.size()) {
            TBOX_ERROR(
               d_object_name << ":  "
                             << "The number of alpha_1, alpha_2, and beta "
                             << "values specified in input is not consistent");
         }

         d_order = static_cast<int>(d_alpha_1.size());
      }
   }
}

/*
 *************************************************************************
 *
 * Checks that class and restart file version numbers are equal.  If so,
 * reads in d_alpha_1, d_alpha_2, and d_beta from the database.  Also,
 * dooes a consistency check to make sure that the number of alpha values
 * specified equals the order of the Runga-Kutta scheme.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file." << std::endl);
   }
   std::shared_ptr<tbox::Database> restart_db(
      root_db->getDatabase(d_object_name));

   int ver = restart_db->getInteger("ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION");
   if (ver != ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version."
                       << std::endl);
   }

   d_alpha_1 = restart_db->getDoubleVector("alpha_1");
   d_alpha_2 = restart_db->getDoubleVector("alpha_2");
   d_beta = restart_db->getDoubleVector("beta");

   if (d_alpha_1.size() != d_alpha_2.size() ||
       d_alpha_2.size() != d_beta.size()) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "The number of alpha_1, alpha_2, and beta values "
                       << "specified in restart is not consistent"
                       << std::endl);
   }

   d_order = static_cast<int>(d_alpha_1.size());

}

/*
 *************************************************************************
 *
 * Copy all solution data from current context to scratch context.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::copyCurrentToScratch(
   const std::shared_ptr<hier::PatchLevel>& level) const
{
   TBOX_ASSERT(level);

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      std::list<std::shared_ptr<hier::Variable> >::const_iterator soln_var =
         d_soln_variables.begin();
      while (soln_var != d_soln_variables.end()) {

         std::shared_ptr<hier::PatchData> src_data(
            patch->getPatchData(*soln_var, d_current));

         std::shared_ptr<hier::PatchData> dst_data(
            patch->getPatchData(*soln_var, d_scratch));

         dst_data->copy(*src_data);
         ++soln_var;

      }

   }

}

/*
 *************************************************************************
 *
 * Copy all solution data from scratch context to current context.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::copyScratchToCurrent(
   const std::shared_ptr<hier::PatchLevel>& level) const
{
   TBOX_ASSERT(level);

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      std::list<std::shared_ptr<hier::Variable> >::const_iterator soln_var =
         d_soln_variables.begin();
      while (soln_var != d_soln_variables.end()) {

         std::shared_ptr<hier::PatchData> src_data(
            patch->getPatchData(*soln_var, d_scratch));

         std::shared_ptr<hier::PatchData> dst_data(
            patch->getPatchData(*soln_var, d_current));

         dst_data->copy(*src_data);
         ++soln_var;

      }

   }

}

/*
 *************************************************************************
 *
 * Print all class data members for MethodOfLinesIntegrator object.
 *
 *************************************************************************
 */

void
MethodOfLinesIntegrator::printClassData(
   std::ostream& os) const
{
   os << "\nMethodOfLinesIntegrator::printClassData..." << std::endl;
   os << "\nMethodOfLinesIntegrator: this = "
      << (MethodOfLinesIntegrator *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_order = " << d_order << std::endl;

   for (int j = 0; j < d_order; ++j) {
      os << "d_alpha_1[" << j << "] = " << d_alpha_1[j] << std::endl;
      os << "d_alpha_2[" << j << "] = " << d_alpha_2[j] << std::endl;
      os << "d_beta[" << j << "] = " << d_beta[j] << std::endl;
   }

   os << "d_patch_strategy = "
      << (MethodOfLinesPatchStrategy *)d_patch_strategy << std::endl;
}

}
}
