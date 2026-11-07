/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Time integration manager for AMR with local time stepping.
 *
 ************************************************************************/
#include "SAMRAI/algs/TimeRefinementIntegrator.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

#include <cstdlib>
#include <fstream>

#include <cmath>

// #define DEBUG_TIMES

namespace SAMRAI {
namespace algs {

const int TimeRefinementIntegrator::ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION = 2;

tbox::StartupShutdownManager::Handler
TimeRefinementIntegrator::s_initialize_handler(
   TimeRefinementIntegrator::initializeCallback,
   0,
   0,
   TimeRefinementIntegrator::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 * tbox::Timer objects for performance measurement.
 */
std::shared_ptr<tbox::Timer> TimeRefinementIntegrator::t_initialize_hier;
std::shared_ptr<tbox::Timer> TimeRefinementIntegrator::t_advance_hier;
std::shared_ptr<tbox::Timer> TimeRefinementIntegrator::t_advance_level;

/*
 *************************************************************************
 *
 * The constructor for TimeRefinementIntegrator checks for valid
 * input data, initializes time stepping data to undefined values,
 * and forces certain parameters in the level strategy, regridding
 * algorithm to be consistent with this object's data members.  This
 * constructor also invokes the variable registration process in the
 * level strategy object.
 *
 *************************************************************************
 */

TimeRefinementIntegrator::TimeRefinementIntegrator(
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::shared_ptr<TimeRefinementLevelStrategy>& level_integrator,
   const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_algorithm):
   d_object_name(object_name),
   d_patch_hierarchy(hierarchy),
   d_refine_level_integrator(level_integrator),
   d_gridding_algorithm(gridding_algorithm),
   d_start_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_end_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_grow_dt(1.0),
   d_integrator_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_just_regridded(false),
   d_level_0_advanced(false),
   d_hierarchy_advanced(false),
   d_connector_width_requestor(),
   d_barrier_and_time(false)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_integrator);
   TBOX_ASSERT(gridding_algorithm);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
      this);

   d_use_refined_timestepping = level_integrator->usingRefinedTimestepping();

   const int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
   d_regrid_interval.resize(max_levels);

   d_level_old_old_time.resize(max_levels);
   d_level_old_time.resize(max_levels);
   d_level_sim_time.resize(max_levels);
   d_dt_max_level.resize(max_levels);
   d_dt_actual_level.resize(max_levels);
   d_step_level.resize(max_levels);
   d_max_steps_level.resize(max_levels);

   int level_number;

   for (level_number = 0; level_number < max_levels; ++level_number) {
      d_regrid_interval[level_number] = tbox::MathUtilities<int>::getMax();
      d_level_old_old_time[level_number] =
         tbox::MathUtilities<double>::getSignalingNaN();
      d_level_old_time[level_number] =
         tbox::MathUtilities<double>::getSignalingNaN();
      d_level_sim_time[level_number] =
         tbox::MathUtilities<double>::getSignalingNaN();
      d_dt_max_level[level_number] =
         tbox::MathUtilities<double>::getSignalingNaN();
      d_dt_actual_level[level_number] =
         tbox::MathUtilities<double>::getSignalingNaN();
      d_step_level[level_number] = tbox::MathUtilities<int>::getMax();
      d_max_steps_level[level_number] = tbox::MathUtilities<int>::getMax();
   }

   /*
    * Set regrid interval data based on ratios between levels or input value.
    */

   if (d_use_refined_timestepping) {
      if (max_levels > 1) {

         for (level_number = 1; level_number < max_levels; ++level_number) {
            const hier::IntVector ratio(d_patch_hierarchy->
                                        getRatioToCoarserLevel(level_number));

            if (((ratio.max() % 2) == 0) || (ratio.max() == 1)) {
               d_regrid_interval[level_number] = 2;
            } else {
               if ((ratio.max() % 3) == 0) {
                  d_regrid_interval[level_number] = 3;
               } else {
                  TBOX_ERROR(
                     d_object_name << ":  "
                                   << " integrator cannot set regrid interval"
                                   << " based on ratios between levels."
                                   << std::endl);
               }
            }
         }

         d_regrid_interval[0] = d_regrid_interval[1];

      } else {
         for (level_number = 0; level_number < max_levels; ++level_number) {
            d_regrid_interval[level_number] = 1;
         }
      }
   } else {
      /* Set the regrid interval to the default.  This can be overridden by any
       * input database in getFromInput.
       */
      setRegridInterval(1);
   }

   /*
    * Initialize this time refinement integration object with data read
    * from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   d_connector_width_requestor.setTagBuffer(d_tag_buffer);
   hierarchy->registerConnectorWidthRequestor(
      d_connector_width_requestor);

   /*
    * Initialize remaining integrator data members.
    */

   if (!is_from_restart) {
      d_integrator_time = d_start_time;
      d_step_level[0] = 0;
      d_last_finest_level = 0;
   }

   tbox::plog << "TimeRefinementIntegrator constructor setting regrid intervals:";
   for (size_t i = 0; i < d_regrid_interval.size(); ++i) {
      tbox::plog << "  [" << i << "]=" << d_regrid_interval[i];
   }
   tbox::plog << "\n";
}

/*
 *************************************************************************
 *
 * Destructor tells tbox::RestartManager to remove this object from the
 * list of restart items.
 *
 *************************************************************************
 */

TimeRefinementIntegrator::~TimeRefinementIntegrator()
{
   tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/*
 *************************************************************************
 *
 * Create patch level and patches for coarsest level in AMR hierarchy
 * at initial simulation time.  Then, set time stepping data and set
 * initial patch data on the coarsest level.  The actual level data
 * initialization process is invoked by the recursive private member
 * function initialize().  This function will create and initialize
 * successively finer hierarchy levels until either the maximum number
 * allowable level is reached or no further refinement is needed.
 *
 *************************************************************************
 */

double
TimeRefinementIntegrator::initializeHierarchy()
{
   if (d_barrier_and_time) {
      t_initialize_hier->barrierAndStart();
   }

   d_refine_level_integrator->initializeLevelIntegrator(d_gridding_algorithm);

   if (tbox::RestartManager::getManager()->isFromRestart()) {

      d_patch_hierarchy->initializeHierarchy();

      d_gridding_algorithm->
      getTagAndInitializeStrategy()->
      resetHierarchyConfiguration(d_patch_hierarchy,
         0,
         d_patch_hierarchy->getFinestLevelNumber());

   } else {

      d_gridding_algorithm->makeCoarsestLevel(
         d_start_time);

      if (d_use_refined_timestepping) {
         initializeRefinedTimesteppingLevelData(0);
      } else {
         initializeSynchronizedTimesteppingLevelData(0);
      }

      /*
       * After data on each level is initialized at simulation start time,
       * coarser levels are synchronized with finer levels that didn't exist
       * when the coarser level initial data was set.  This synchronization
       * process is defined by the integration algorithm.
       */
      if (d_patch_hierarchy->getFinestLevelNumber() > 0) {
         // "true" argument: const bool initial_time = true;
         d_refine_level_integrator->synchronizeNewLevels(d_patch_hierarchy,
            0,
            d_patch_hierarchy->getFinestLevelNumber(),
            d_start_time,
            true);
      }
   }

   d_level_0_advanced = true;
   d_hierarchy_advanced = true;

   if (d_barrier_and_time) {
      t_initialize_hier->stop();
   }

   return d_dt_max_level[0];

}

/*
 *************************************************************************
 *
 * Advance all levels in hierarchy through specified time interval dt
 * by invoking the recursive timestepping process at the coarsest
 * hierarchy level (level 0).  If synchronized timestepping is used,
 * then loop through each level advance and then synchronize all levels.
 * The return value is the proper time increment for a subsequent
 * advance of level 0.
 *
 *************************************************************************
 */

double
TimeRefinementIntegrator::advanceHierarchy(
   const double dt,
   const bool rebalance_coarsest)
{
   TBOX_ASSERT(dt >= 0.0);

   if (d_barrier_and_time) {
      t_advance_hier->barrierAndStart();
   }

   d_level_sim_time[0] = d_integrator_time;

   SAMRAI_CALI_MARK_BEGIN("makeCoarsestLevel");
   if (rebalance_coarsest) {
      d_gridding_algorithm->makeCoarsestLevel(
         d_level_sim_time[0]);
   }
   SAMRAI_CALI_MARK_END("makeCoarsestLevel");

   double dt_new;

   if (!d_use_refined_timestepping) {

      SAMRAI_CALI_MARK_BEGIN("adavanceForSync");
      dt_new = advanceForSynchronizedTimestepping(d_level_sim_time[0] + dt);
      SAMRAI_CALI_MARK_END("adavanceForSync");

      if (d_integrator_time + dt_new > d_end_time) {
         dt_new = d_end_time - d_integrator_time;
      }

   } else {
      SAMRAI_CALI_MARK_BEGIN("adavanceRecursively");
      advanceRecursivelyForRefinedTimestepping(0, d_level_sim_time[0] + dt);
      SAMRAI_CALI_MARK_END("adavanceRecursively");

      d_integrator_time += dt;
      dt_new = tbox::MathUtilities<double>::Min(d_dt_actual_level[0],
            d_end_time - d_integrator_time);
   }

   if (d_barrier_and_time) {
      t_advance_hier->stop();
   }

   return dt_new;
}

/*
 *************************************************************************
 *
 * Initialize data on given level.  Then, invoke regridding and
 * recursively initialize finer levels, if necessary.  The coarsest
 * level (i.e., level 0) is created by initializeHierarchy().  Each
 * finer level is constructed by calls this function.  The gridding
 * algorithm data member build finer each level and its patches.
 *
 * This function assumes that the patch level and patches exist before
 * it is called.  Upon leaving this routine, initial data is set
 * according to the methods in the level integration routines, which
 * are generally invoked from the gridding algorithm class.  Also,
 * basic time increment data is set for each level.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::initializeRefinedTimesteppingLevelData(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   const std::shared_ptr<hier::PatchLevel> patch_level(
      d_patch_hierarchy->getPatchLevel(level_number));

   /*
    * Initialize step count and start time for current level.
    * Then, query level strategy for proper time step.
    */

   if (level_number > 0) {
      d_step_level[level_number] = 0;
      d_max_steps_level[level_number] = 1;
   }
   d_level_sim_time[level_number] = d_start_time;

   bool initial_time = true;
   double dt_level;

   dt_level = d_refine_level_integrator->getLevelDt(patch_level,
         d_level_sim_time[level_number],
         initial_time);

   /*
    * If a coarser level exists, we require that the time increment on
    * the current level be consistent with the needs of the numerical
    * integration algorithm (i.e., the level strategy).  Futhermore,
    * if the current level can be refined and the level is not the coarsest
    * in the hierarchy, we require that the time increment be adjusted so
    * that the number of steps taken on the level before synchronization
    * with the next coarser level is an integer multiple of the grid
    * refinement ratio between the two levels.
    */

   if (level_number > 0) {
      dt_level = tbox::MathUtilities<double>::Min(dt_level,
            d_refine_level_integrator->getMaxFinerLevelDt(
               level_number,
               d_dt_actual_level[level_number - 1],
               d_patch_hierarchy->getRatioToCoarserLevel(level_number)));
   }

   d_dt_max_level[level_number] = d_dt_actual_level[level_number] = dt_level;

   /*
    * Create a finer level if appropriate.
    */

   if (d_patch_hierarchy->levelCanBeRefined(level_number)) {

      int tag_buffer;
      if (d_gridding_algorithm->getTagAndInitializeStrategy()->
          usesTimeIntegration(d_step_level[0], d_integrator_time)) {
         tag_buffer = d_regrid_interval[level_number];
      } else {
         tag_buffer = d_tag_buffer[level_number];
      }

      double regrid_start_time =
         d_level_sim_time[level_number] - d_dt_actual_level[level_number];
      // "true" argument: const bool initial_time = true;
      d_gridding_algorithm->makeFinerLevel(
         tag_buffer,
         true,
         d_step_level[0],
         d_level_sim_time[level_number],
         regrid_start_time);

      /*
       * If new finer level is made, data on its patches is initialized.
       * Also, if new level can be refined and time integration is used
       * during regridding process, data on the current level is advanced
       * through a single time increment to provide boundary data for the
       * regridding process on finer levels.
       */

      if (d_patch_hierarchy->finerLevelExists(level_number)) {

         /*
          * If time integration is used in the refinement process and
          * the newly created finer level can be further refined,
          * data is advanced on the current level through one time
          * increment to provide time interpolated boundary data for
          * regridding on the finer level.
          */
         if (d_patch_hierarchy->levelCanBeRefined(level_number + 1) &&
             d_gridding_algorithm->getTagAndInitializeStrategy()->
             usesTimeIntegration(d_step_level[0], d_integrator_time)) {

            if (d_barrier_and_time) {
               t_advance_level->barrierAndStart();
            }
            // "false" argument: bool last_step = false;
            // "true" argument: bool regrid_advance = true;
            d_refine_level_integrator->
            advanceLevel(patch_level,
               d_patch_hierarchy,
               d_level_sim_time[level_number],
               d_level_sim_time[level_number]
               + d_dt_actual_level[level_number],
               firstLevelStep(level_number),
               false,
               true);
            if (d_barrier_and_time) {
               t_advance_level->stop();
            }
         }

         /*
          * RECURSIVE invocation of initialization on next finer level.
          */

         initializeRefinedTimesteppingLevelData(level_number + 1);

         /*
          * If current level patch data was advanced above, all new
          * data is discarded so that advance routines all start with the
          * same patch data on each level, and because AMR synchronization
          * may require a different initial time increment size.
          */

         if (d_patch_hierarchy->levelCanBeRefined(level_number + 1) &&
             d_gridding_algorithm->getTagAndInitializeStrategy()->
             usesTimeIntegration(d_step_level[0], d_integrator_time)) {
            d_refine_level_integrator->
            resetDataToPreadvanceState(patch_level);
         }
      }
   }
}

void
TimeRefinementIntegrator::initializeSynchronizedTimesteppingLevelData(
   const int level_number)
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   const std::shared_ptr<hier::PatchLevel> patch_level(
      d_patch_hierarchy->getPatchLevel(level_number));

   /*
    * Initialize step count and start time for current level.
    * Then, query level strategy for proper time step.
    */

   if (level_number > 0) {
      d_step_level[level_number] = 0;
      d_max_steps_level[level_number] = 1;
   }
   d_level_sim_time[level_number] = d_start_time;
   d_level_old_time[level_number] = d_start_time;

   bool initial_time = true;
   double dt_level;

   dt_level = d_refine_level_integrator->getLevelDt(patch_level,
         d_level_sim_time[level_number],
         initial_time);

   if (level_number == 0) {
      d_dt = dt_level;
   } else {
      d_dt = tbox::MathUtilities<double>::Min(d_dt, dt_level);
   }
   d_dt_max_level[level_number] = d_dt_actual_level[level_number] = d_dt;

   for (int i = 0; i < level_number; ++i) {
      if (d_dt_max_level[i] > d_dt_max_level[level_number]) {
         d_dt_max_level[i] = d_dt_actual_level[i] =
               d_dt_max_level[level_number];
      }
   }

   /*
    * Create a finer level if appropriate.
    */

   if (d_patch_hierarchy->levelCanBeRefined(level_number)) {

      int tag_buffer;
      if (d_gridding_algorithm->getTagAndInitializeStrategy()->
          usesTimeIntegration(d_step_level[0], d_integrator_time)) {
         tag_buffer = d_regrid_interval[level_number];
      } else {
         tag_buffer = d_tag_buffer[level_number];
      }

      double regrid_start_time =
         d_level_sim_time[level_number] - d_dt_actual_level[level_number];
      // "true" argument: const bool initial_time = true;
      d_gridding_algorithm->makeFinerLevel(
         tag_buffer,
         true,
         d_step_level[0],
         d_level_sim_time[level_number],
         regrid_start_time);

      /*
       * If new finer level is made, data on its patches is initialized.
       * Also, if new level can be refined and time integration is used
       * during regridding process, data on the current level is advanced
       * through a single time increment to provide boundary data for the
       * regridding process on finer levels.
       */

      if (d_patch_hierarchy->finerLevelExists(level_number)) {

         /*
          * If time integration is used in the refinement process and
          * the newly created finer level can be further refined,
          * data is advanced on the current level through one time
          * increment to provide time interpolated boundary data for
          * regridding on the finer level.
          */
         if (d_patch_hierarchy->levelCanBeRefined(level_number + 1) &&
             d_gridding_algorithm->getTagAndInitializeStrategy()->
             usesTimeIntegration(d_step_level[0], d_integrator_time)) {

            if (d_barrier_and_time) {
               t_advance_level->barrierAndStart();
            }
            // "false" argument: bool last_step = false;
            // "true" argument: bool regrid_advance = true;
            d_refine_level_integrator->
            advanceLevel(patch_level,
               d_patch_hierarchy,
               d_level_sim_time[level_number],
               d_level_sim_time[level_number]
               + d_dt_actual_level[level_number],
               firstLevelStep(level_number),
               false,
               true);

            if (d_barrier_and_time) {
               t_advance_level->stop();
            }
         }

         /*
          * RECURSIVE invocation of initialization on next finer level.
          */

         initializeSynchronizedTimesteppingLevelData(level_number + 1);

         /*
          * If current level patch data was advanced above, all new
          * data is discarded so that advance routines all start with the
          * same patch data on each level, and because AMR synchronization
          * may require a different initial time increment size.
          */
         if (d_patch_hierarchy->levelCanBeRefined(level_number + 1) &&
             d_gridding_algorithm->getTagAndInitializeStrategy()->
             usesTimeIntegration(d_step_level[0], d_integrator_time)) {
            d_refine_level_integrator->
            resetDataToPreadvanceState(patch_level);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Advance data on level to specified time.  Then, advance each finer
 * hierarchy level to the same end time using recursive function calls.
 * It is assumed that when this function is called, only the data needed
 * for initialization exists on the level.  Also, the solutions on
 * specified level and all finer levels are synchronized at the the
 * current simulation time on level level_number-1.  The integration
 * process implemented in this function is outlined as follows:
 *
 *    1) Initialize timestep count, simulation time, time step increment
 *       for level level_number.
 *    2) Adjust time first time increment if necessary and estimate
 *       number of timesteps required to advance the level to end time.
 *    3) Iterate over sequence of timesteps for level level_number:
 *       a) Set end time for current step, record current time for
 *          synchronization.
 *       b) Advance solution on level using given time increment.
 *          Note that the level strategy performs the advance.
 *       c) Increment step count information.
 *       d) If finer level exists, advance finer level to end time of
 *          most recent time step using recursive function call.
 *       e) Synchronize data between levels as necessary.  Note that
 *          new solution on level level_number as well as each finer
 *          level corresponds to the same time.  Also, synchronization
 *          involves several levels, typically, including all levels
 *          finer than level_number, level level_number itself, and the
 *          possibly level level_number-1.  Which levels synchronize
 *          here depends on when regridding occurs.  Note that the
 *          level strategy performs the data synchronization.
 *       f) Increment the simulation time and determine the size of the
 *          next timestep for the level and estimate the number of
 *          timesteps needed to advance the level to end time.
 *       g) If appropriate, regrid all finer levels.
 *          Note that the gridding algorithm regrids the mesh.
 *       h) Synchronize all levels involved in the regridding process,
 *          if it occurred.  The sync is performed by level strategy.
 *       i) If regridding did not occur, then reset current and new
 *          solution data on the level.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::advanceRecursivelyForRefinedTimestepping(
   const int level_number,
   const double end_time)
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(end_time >= d_integrator_time);

   const std::shared_ptr<hier::PatchLevel> patch_level(
      d_patch_hierarchy->getPatchLevel(level_number));

   /*
    * Initialize step count, start time for current level.
    * Determine time remaining on level and, if necessary, set
    * the max time increment for the level.
    */

   if (level_number > 0) {
      d_step_level[level_number] = 0;
      d_max_steps_level[level_number] = 1;
   } else {
      d_level_0_advanced = false;
      d_hierarchy_advanced = false;
   }
   double time_remaining = 0.0;

   if (level_number > 0) {
      d_level_sim_time[level_number] = d_level_sim_time[level_number - 1];
      time_remaining = end_time - d_level_sim_time[level_number];

      /*
       * If this level was created during last regrid, compute max time
       * increment for level.
       */

      if (d_last_finest_level < level_number) {
         d_dt_max_level[level_number] =
            tbox::MathUtilities<double>::Min(
               d_dt_actual_level[level_number - 1],
               d_refine_level_integrator->getMaxFinerLevelDt(
                  level_number,
                  d_dt_actual_level[level_number - 1],
                  d_patch_hierarchy->getRatioToCoarserLevel(level_number)));

         if (d_patch_hierarchy->levelCanBeRefined(level_number)) {
            d_dt_max_level[level_number] =
               tbox::MathUtilities<double>::Min(d_dt_max_level[level_number],
                  time_remaining
                  / double(d_regrid_interval[level_number]));
         }
         d_last_finest_level = level_number;
      }

   } else {
      time_remaining = end_time - d_level_sim_time[level_number];
   }

   /*
    * Determine time increment for first step in integration sequence,
    * and estimate the number of steps before synchronization with the
    * next coarser level occurs.
    */

   bool sync_after_step =
      findNextDtAndStepsRemaining(level_number,
         time_remaining,
         d_dt_max_level[level_number]);

   /*
    * Loop over a dynamically determined sequence of timesteps on the
    * current level (level_number).  Note that if level is coarsest in
    * AMR hierarchy, we must determine whether there will only be a
    * single advance step before the function returns.
    */

   double new_level_time, dt_new;
   while (!lastLevelStep(level_number)) {

      /*
       * Determine new time after this step, set old time for synchronization.
       */

      new_level_time = (sync_after_step ? end_time
                        : d_level_sim_time[level_number]
                        + d_dt_actual_level[level_number]);
      d_level_old_old_time[level_number] = d_level_old_time[level_number];
      d_level_old_time[level_number] = d_level_sim_time[level_number];
      d_just_regridded = false;

#ifdef DEBUG_TIMES
      tbox::plog << "\nAdvancing level number = " << level_number << std::endl;
      tbox::plog << "step number = " << d_step_level[level_number] << std::endl;
      tbox::plog << "max steps = " << d_max_steps_level[level_number]
                 << std::endl;
      tbox::plog << "current time = " << d_level_sim_time[level_number]
                 << std::endl;
      tbox::plog << "dt used = " << d_dt_actual_level[level_number]
                 << std::endl;
      tbox::plog << "new level time = " << new_level_time << std::endl;
      tbox::plog << "dt max = " << d_dt_max_level[level_number] << std::endl;
      tbox::plog << "end time = " << end_time << std::endl;
      tbox::plog << "sync_after_step = " << sync_after_step << std::endl;
#endif

      /*
       * Advance level from current simulation time to new_level_time
       * using a single time advance step.  Note that the level strategy
       * returns the next time increment for the level.  Also, we keep both
       * new and previous data on level so that time interpolation can be
       * used to set boundary conditions for finer levels and for proper
       * data synchronization once all finer levels have been advanced.
       */

      if (d_barrier_and_time) {
         t_advance_level->barrierAndStart();
      }
      // "sync_after_step" is same as "last_step" in level strategy.
      dt_new = d_refine_level_integrator->advanceLevel(patch_level,
            d_patch_hierarchy,
            d_level_sim_time[level_number],
            new_level_time,
            firstLevelStep(level_number),
            sync_after_step);

      if (d_barrier_and_time) {
         t_advance_level->stop();
      }

      /*
       * Update step count information.  Then, advance all finer levels.
       */

      if (level_number == 0) {
         d_level_0_advanced = true;
      } else {
         ++d_step_level[level_number];
      }

      if (d_patch_hierarchy->finerLevelExists(level_number)) {
         advanceRecursivelyForRefinedTimestepping(level_number + 1,
            new_level_time);
      }

      if (level_number == 0) {
         ++d_step_level[level_number];
         d_hierarchy_advanced = true;
      }

      /*
       * Synchronize data between levels in the hierarchy as necessary.
       * Note that this process synchronizes data between this level,
       * several finer levels, and the next coarser level, potentially.
       */

      int coarsest_sync_level = -1;
      int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();

      if (atRegridPoint(level_number)) {

         if (!lastLevelStep(level_number)
             || !coarserLevelRegridsToo(level_number)) {

            coarsest_sync_level = (((level_number > 0)
                                    && lastLevelStep(level_number))
                                   ? level_number - 1 : level_number);

            if (coarsest_sync_level < finest_level_number) {
#ifdef DEBUG_TIMES
               tbox::plog << "\nSynchronizing levels " << coarsest_sync_level
                          << " to "
                          << finest_level_number << std::endl;
#endif
               d_refine_level_integrator->
               standardLevelSynchronization(d_patch_hierarchy,
                  coarsest_sync_level,
                  finest_level_number,
                  new_level_time,
                  d_level_old_time);
            }

         }

      } else {

         if (level_number < finest_level_number) {
            if ((!lastLevelStep(level_number)
                 || (level_number == 0)) && !d_just_regridded) {
#ifdef DEBUG_TIMES
               tbox::plog << "\nSynchronizing levels " << level_number
                          << " to "
                          << finest_level_number << std::endl;
#endif
               d_refine_level_integrator->
               standardLevelSynchronization(d_patch_hierarchy,
                  level_number,
                  finest_level_number,
                  new_level_time,
                  d_level_old_time);

               d_refine_level_integrator->
               resetTimeDependentData(d_patch_hierarchy->
                  getPatchLevel(finest_level_number),
                  new_level_time,
                  d_patch_hierarchy->levelCanBeRefined(finest_level_number));
            }
         }

      }

      /*
       * Update level simulation time and time remaining until
       * synchronization with next coarser level.  Then, adjust
       * time increment and step sequence for current level.
       */

      d_level_sim_time[level_number] = new_level_time;
      time_remaining = end_time - new_level_time;

      sync_after_step = findNextDtAndStepsRemaining(level_number,
            time_remaining,
            dt_new);

      /*
       * All finer levels are synchronized with this level now.
       * If appropriate, we regrid finer levels and re-synchronize
       * levels as needed.  Otherwise, we reset time-dependent data and
       * re-synchronize levels as needed.  Note that the regridding
       * process resets the data on each level involved in the regridding.
       */

      if (atRegridPoint(level_number)) {

         if (!lastLevelStep(level_number)
             || !coarserLevelRegridsToo(level_number)) {
#ifdef DEBUG_TIMES
            tbox::plog << "\nRegridding from level number = "
                       << level_number << std::endl;
#endif
            /*
             * Reset time dependent data.  If the gridding algorithm uses
             * time integration for error estimation, it will have already
             * reset time dependent data on all levels regridded, so
             * only reset data on levels that are not regridded.  If the
             * gridding algorithm does not used time integration, reset data
             * on all levels.
             */
            if (d_gridding_algorithm->getTagAndInitializeStrategy()->
                usesTimeIntegration(d_step_level[0], d_integrator_time)) {
               if (!d_patch_hierarchy->
                   levelCanBeRefined(finest_level_number)) {
                  d_refine_level_integrator->resetTimeDependentData(
                     d_patch_hierarchy->getPatchLevel(finest_level_number),
                     new_level_time,
                     d_patch_hierarchy->levelCanBeRefined(finest_level_number));
               }
            } else {
               for (int ln = level_number; ln <= finest_level_number; ++ln) {
                  d_refine_level_integrator->resetTimeDependentData(
                     d_patch_hierarchy->getPatchLevel(ln),
                     d_level_sim_time[ln],
                     d_patch_hierarchy->levelCanBeRefined(ln));
               }
            }

            d_last_finest_level = finest_level_number;

            /*
             * Regrid finer levels.  If the error estimation procedure
             * uses time integration (e.g. Richardson extrapolation) then
             * we must supply the oldest time at which data is stored.
             *
             * If the error coarsen ratio is two, data will be stored
             * from the previous timestep (at d_level_old_time).  If the
             * error coarsen ratio is three, data will be stored
             * from two previous timesteps (at d_level_old_old_time).
             *
             * If we are not using time integration, the oldest time
             * information should not be used, so it is set to NaNs
             * to throw an assertion if it is accessed.
             */

            std::vector<double> regrid_start_time;
            if (!d_gridding_algorithm->getTagAndInitializeStrategy()->
                usesTimeIntegration(d_step_level[0], d_integrator_time)) {

               int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
               regrid_start_time.resize(max_levels);
               int array_size = static_cast<int>(regrid_start_time.size());
               for (int i = 0; i < array_size; ++i) {
                  regrid_start_time[i] = 0.;
               }

            } else {

               if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio() ==
                   2) {
                  regrid_start_time = d_level_old_time;
               } else if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio()
                          == 3) {
                  regrid_start_time = d_level_old_old_time;
               } else {
                  TBOX_ERROR(
                     d_object_name << ": the supplied gridding "
                                   << "algorithm uses an error coarsen ratio of "
                                   << d_gridding_algorithm->
                     getTagAndInitializeStrategy()->getErrorCoarsenRatio()
                                   << " which is not supported in this class"
                                   << std::endl);
               }

            }

            d_gridding_algorithm->
            regridAllFinerLevels(
               level_number,
               d_tag_buffer,
               d_step_level[0],
               d_level_sim_time[level_number],
               regrid_start_time,
               (coarsest_sync_level >= level_number));

            d_just_regridded = true;

            if (level_number < d_patch_hierarchy->getFinestLevelNumber()) {
#ifdef DEBUG_TIMES
               tbox::plog << "\nSynchronizing levels after regrid : "
                          << level_number << " to "
                          << d_patch_hierarchy->getFinestLevelNumber()
                          << std::endl;
#endif

               // "false" argument: const bool initial_time = false;
               d_refine_level_integrator->
               synchronizeNewLevels(d_patch_hierarchy,
                  level_number,
                  d_patch_hierarchy->getFinestLevelNumber(),
                  d_level_sim_time[level_number],
                  false);
            }

         }

      } else {

         if (!lastLevelStep(level_number) || (level_number == 0)) {
            d_refine_level_integrator->resetTimeDependentData(
               patch_level,
               d_level_sim_time[level_number],
               d_patch_hierarchy->levelCanBeRefined(level_number));
         }

         if (d_just_regridded) {
#ifdef DEBUG_TIMES
            tbox::plog << "\nSynchronizing levels after regrid : "
                       << level_number << " to "
                       << d_patch_hierarchy->getFinestLevelNumber()
                       << std::endl;
#endif
            // "false" argument: const bool initial_time = false;
            d_refine_level_integrator->
            synchronizeNewLevels(d_patch_hierarchy,
               level_number,
               d_patch_hierarchy->getFinestLevelNumber(),
               d_level_sim_time[level_number],
               false);
         }

      }

   }

}

/*
 *************************************************************************
 *
 * Advance data on all hierarchy levels through specified time increment
 * using the same timesteps on each level.  Synchronization and
 * regridding are performed as needed.
 *
 *************************************************************************
 */

double
TimeRefinementIntegrator::advanceForSynchronizedTimestepping(
   const double end_time)
{
   TBOX_ASSERT(end_time >= d_integrator_time);

   double dt = end_time - d_integrator_time;

   int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
   double dt_new = tbox::MathUtilities<double>::getMax();

   d_level_0_advanced = false;
   d_hierarchy_advanced = false;

   int level_num;
   for (level_num = 0; level_num <= finest_level_number; ++level_num) {

      std::shared_ptr<hier::PatchLevel> patch_level(
         d_patch_hierarchy->getPatchLevel(level_num));

      if (level_num > 0) {
         d_step_level[level_num] = 1;
         d_max_steps_level[level_num] = 1;
      }
      d_dt_max_level[level_num] = dt;
      d_level_sim_time[level_num] = d_integrator_time;

#ifdef DEBUG_TIMES
      tbox::plog << "\nAdvancing level number = " << level_num << std::endl;
      tbox::plog << "step number = " << d_step_level[level_num] << std::endl;
      tbox::plog << "max steps = " << d_max_steps_level[level_num] << std::endl;
      tbox::plog << "current time = " << d_integrator_time << std::endl;
      tbox::plog << "dt used = " << dt << std::endl;
      tbox::plog << "new level time = " << d_integrator_time + dt << std::endl;
#endif

      if (d_barrier_and_time) {
         t_advance_level->barrierAndStart();
      }
      // "true" argument: bool first_step = true;
      // "false" argument: bool last_step = false;
      RANGE_PUSH("advanceLevel", 1);
      double dt_next_level =
         d_refine_level_integrator->advanceLevel(patch_level,
            d_patch_hierarchy,
            d_integrator_time,
            d_integrator_time + dt,
            true,
            false);
      RANGE_POP;

      if (level_num == 0) {
         d_level_0_advanced = true;
      }

      if (d_barrier_and_time) {
         t_advance_level->stop();
      }

      dt_new = tbox::MathUtilities<double>::Min(dt_new, dt_next_level);
   }

   dt_new = tbox::MathUtilities<double>::Min(dt_new, d_grow_dt * dt);

   for (level_num = 0; level_num <= finest_level_number; ++level_num) {
      d_dt_max_level[level_num] = d_dt_actual_level[level_num] = dt_new;
   }

   d_integrator_time += dt;
   ++d_step_level[0];
   d_hierarchy_advanced = true;

   int coarse_level_number = 0;

   SAMRAI_CALI_MARK_BEGIN("syncLevels");
   if (finest_level_number > 0) {

#ifdef DEBUG_TIMES
      tbox::plog << "\nSynchronizing levels " << coarse_level_number
                 << " to " << finest_level_number << std::endl;
#endif
      std::vector<double> old_times(finest_level_number + 1,
                                    d_integrator_time - dt);

      d_refine_level_integrator->standardLevelSynchronization(
         d_patch_hierarchy,
         0,
         finest_level_number,
         d_integrator_time,
         old_times);
   }
   SAMRAI_CALI_MARK_END("syncLevels");
   /*
    * Store the time from the previous step(s) in the "old_time"
    * array.  This information may be used during the re-gridding
    * process, if time integration is used during error estimation.
    */
   for (level_num = 0; level_num <= finest_level_number; ++level_num) {
      d_level_old_old_time[level_num] = d_level_old_time[level_num];
      d_level_old_time[level_num] = d_level_sim_time[level_num];
   }

   /*
    * Are we ready to re-grid??
    */
   bool regrid_now = (d_step_level[0] % d_regrid_interval[0] == 0);

   RANGE_PUSH("regrid", 1);
   if (!regrid_now) {

      /*
       * We are not re-gridding, so all we need to do is simply
       * reset the time dependent data on all the levels to
       * prepare for the next advance.
       */
      for (int ln = 0; ln <= finest_level_number; ++ln) {
         d_refine_level_integrator->resetTimeDependentData(
            d_patch_hierarchy->getPatchLevel(ln),
            d_integrator_time,
            d_patch_hierarchy->levelCanBeRefined(ln));
      }

   } else {

      /*
       * We are re-gridding.  Reset the time dependent data.  If
       * the error estimation procedure does NOT use time
       * integration (e.g. gradient detection) then simply go
       * through all levels of the hierarchy and reset.
       *
       * If the error estimation does use time integration
       * (e.g. Richardson extrapolation) the time dependent
       * data will be accessed during the tagging phase.
       * Hence, we only reset time dependent data on levels
       * that are not tagged (i.e. the finest level) and
       * defer resetting time dependent data on levels that
       * are tagged to the level integrator performing the
       * re-gridding operations. NOTE:  One should assure
       * this is properly done in the level integrator.
       */

      if (d_gridding_algorithm->getTagAndInitializeStrategy()->
          usesTimeIntegration(d_step_level[0], d_integrator_time)) {
         if (!d_patch_hierarchy->levelCanBeRefined(finest_level_number)) {
            d_refine_level_integrator->resetTimeDependentData(
               d_patch_hierarchy->getPatchLevel(finest_level_number),
               d_integrator_time,
               d_patch_hierarchy->levelCanBeRefined(finest_level_number));
         }
      } else {
         for (int ln = 0; ln <= finest_level_number; ++ln) {
            d_refine_level_integrator->resetTimeDependentData(
               d_patch_hierarchy->getPatchLevel(ln),
               d_integrator_time,
               d_patch_hierarchy->levelCanBeRefined(ln));
         }
      }

      /*
       * Regrid all levels, from coarsest to finest.  If the error
       * estimation procedure uses time integration (e.g. Richardson
       * extrapolation) then we must supply the oldest time at which
       * data is stored.
       *
       * If the error coarsen ratio is two, data will be stored
       * from the previous timestep (at d_level_old_time).  If the
       * error coarsen ratio is three, data will be stored
       * from two previous timesteps (at d_level_old_old_time).
       *
       * If we are not using time integration, the oldest time
       * information should not be used, so it is set to NaNs
       * to throw an assertion if it is accessed.
       */

      std::vector<double> regrid_start_time;
      if (!d_gridding_algorithm->getTagAndInitializeStrategy()->
          usesTimeIntegration(d_step_level[0], d_integrator_time)) {

         int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
         regrid_start_time.resize(max_levels);
         int array_size = static_cast<int>(regrid_start_time.size());
         for (int i = 0; i < array_size; ++i) {
            regrid_start_time[i] = 0.;
         }

      } else {

         if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio() == 2) {
            regrid_start_time = d_level_old_time;
         } else if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio() ==
                    3) {
            regrid_start_time = d_level_old_old_time;
         } else {
            TBOX_ERROR(
               d_object_name << ": the supplied gridding "
                             << "algorithm uses an error coarsen ratio of "
                             << d_gridding_algorithm->
               getTagAndInitializeStrategy()->getErrorCoarsenRatio()
                             << " which is not supported in this class"
                             << std::endl);
         }

      }

      d_gridding_algorithm->
      regridAllFinerLevels(
         coarse_level_number,
         d_tag_buffer,
         d_step_level[0],
         d_integrator_time,
         regrid_start_time);

      /*
       * Synchronize data on new levels.
       */
      SAMRAI_CALI_MARK_BEGIN("syncPostRegrid");
      if (d_patch_hierarchy->getFinestLevelNumber() > 0) {


#ifdef DEBUG_TIMES
         tbox::plog << "\nSynchronizing levels after regrid : "
                    << coarse_level_number << " to "
                    << d_patch_hierarchy->getFinestLevelNumber() << std::endl;
#endif
         const bool initial_time = false;
         d_refine_level_integrator->synchronizeNewLevels(
            d_patch_hierarchy,
            coarse_level_number,
            d_patch_hierarchy->getFinestLevelNumber(),
            d_integrator_time,
            initial_time);
      }
      SAMRAI_CALI_MARK_END("syncPostRegrid");
   }
   RANGE_POP;
   return dt_new;
}

/*
 *************************************************************************
 *
 * Compute time step data for given level and estimate the number of
 * time steps needed in current step sequence on that level.  The
 * boolean return value indicates whether the next step taken will
 * be the last in the step sequence on the level.  The outline of the
 * procedure is as follows:
 *
 *    1) Determine maximum time increment allowed, permitting time step
 *       growth if appropriate.
 *    2) If time steps remain in current step sequence:
 *       a) Estimate the number of timesteps left in the sequence based
 *          on the time remaining.
 *       b) Possibly increase the number of time steps remaining so
 *          that total step count remains consistent with regridding
 *          sequence.  Note that the total number of steps must be an
 *          integer multiple of the regrid interval for the level.
 *       c) Compute time increment by dividing number of steps left
 *          into remaining time.
 *    3) If no steps remain on the level in the current sequence, the
 *       next time increment is set to the current maximum increment.
 *
 *************************************************************************
 */

bool
TimeRefinementIntegrator::findNextDtAndStepsRemaining(
   const int level_number,
   const double time_remaining,
   const double dt_bound)
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(time_remaining >= 0.0);
   TBOX_ASSERT(dt_bound >= 0.0);

   /*
    * Grow time increment from previous if possible, but time step size
    * cannot be greater than dt_bound. also, we cannot take any more
    * steps if we have exceeded number allowable.
    */

   d_dt_max_level[level_number] =
      tbox::MathUtilities<double>::Min(dt_bound,
         d_dt_max_level[level_number] * d_grow_dt);

   if (stepsRemaining(level_number)) {

      /*
       * If we have not exceeded the max number of steps, the max number
       * of steps can be adjusted according to the amount of time remaining
       * of the level.
       */

      int number_steps_remaining =
         int(time_remaining / d_dt_max_level[level_number]);

      double dt_temp =
         d_dt_max_level[level_number] * double(number_steps_remaining);

      if (time_remaining - dt_temp >
          sqrt(tbox::MathUtilities<double>::getEpsilon()) * time_remaining) {
         ++number_steps_remaining;
      }

      if (level_number > 0) {
         d_max_steps_level[level_number] =
            d_step_level[level_number] + number_steps_remaining;
      }

      /*
       * If we are not on the coarsest hierarchy level, there must be
       * an integer multiple of the regrid interval of steps in the
       * step sequence.
       */

      if (level_number > 0
          && d_patch_hierarchy->levelCanBeRefined(level_number)) {

         d_max_steps_level[level_number] =
            tbox::MathUtilities<int>::Max(d_max_steps_level[level_number],
               d_regrid_interval[level_number]);

         int number_regrids =
            d_max_steps_level[level_number] / d_regrid_interval[level_number];

         if (d_max_steps_level[level_number]
             % d_regrid_interval[level_number]) ++number_regrids;

         d_max_steps_level[level_number] =
            number_regrids * d_regrid_interval[level_number];

         if (d_step_level[level_number] >= d_max_steps_level[level_number]) {
            TBOX_ERROR(
               d_object_name << ":  "
                             << "no steps left to divide remaining time ...\n"
                             << "level_number = " << level_number
                             << std::endl
                             << "time_remaining = " << time_remaining
                             << "\ndt_bound = " << dt_bound
                             << "\nnumber_steps_remaining = "
                             << number_steps_remaining
                             << std::endl);
         }
      }

      /*
       * Adjust current time increment so that amount of time remaining
       * on the level is divided approximately evenly among the time
       * steps remaining.
       */

      if (level_number == 0) {
         d_dt_actual_level[level_number] = time_remaining;
      } else {
         d_dt_actual_level[level_number] =
            time_remaining / double(d_max_steps_level[level_number]
                                    - d_step_level[level_number]);
      }

   } else {
      d_dt_actual_level[level_number] = d_dt_max_level[level_number];
   }

   if (level_number == 0) {
      return true;
   } else {
      return (d_max_steps_level[level_number]
              - d_step_level[level_number]) <= 1;
   }

}

/*
 *************************************************************************
 *
 * Return true if the level can be remeshed at the current step.
 * regrid step interval.  Otherwise, return false.
 *
 *************************************************************************
 */

bool
TimeRefinementIntegrator::atRegridPoint(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));

   int step_number;
   if (level_number == 0) {
      // If the entire hierarchy has advanced then so has level 0.
      // d_step_level[0] will have advanced by one step and reflects level 0's
      // step number.  However, if the hierarchy has not yet advanced then
      // we're in the middle of a recursive hierarchy advancement and
      // d_step_level[0] has not yet been incremented.  However, level 0 will
      // have advanced by 1 step.  So d_step_level[0] + 1 reflects level 0's
      // step number.
      if (d_hierarchy_advanced) {
         step_number = d_step_level[0];
      } else {
         step_number = d_step_level[0] + 1;
      }
   } else {
      step_number = d_step_level[level_number];
   }

   return (step_number > 0)
          && d_patch_hierarchy->levelCanBeRefined(level_number)
          && (step_number % d_regrid_interval[level_number] == 0);
}

/*
 *************************************************************************
 *
 * Return true if the specified level can be regridded and the the next
 * coarser level can be regridded too.  Otherwise, false is returned.
 *
 *************************************************************************
 */

bool
TimeRefinementIntegrator::coarserLevelRegridsToo(
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
   return (level_number > 0) ? atRegridPoint(level_number - 1) : false;
}

/*
 *************************************************************************
 *************************************************************************
 */

void
TimeRefinementIntegrator::setRegridInterval(
   const int regrid_interval)
{
   TBOX_ASSERT(!d_use_refined_timestepping);
   int array_size = static_cast<int>(d_regrid_interval.size());
   for (int i = 0; i < array_size; ++i) {
      d_regrid_interval[i] = regrid_interval;
   }

   tbox::plog << "TimeRefinementIntegrator::setRegridInterval setting regrid intervals:";
   for (size_t i = 0; i < d_regrid_interval.size(); ++i) {
      tbox::plog << "  [" << i << "]=" << d_regrid_interval[i];
   }
   tbox::plog << "\n";
}

/*
 *************************************************************************
 *
 * Print all data member for TimeRefinementIntegrator object.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::printClassData(
   std::ostream& os) const
{
   os << "\nTimeRefinementIntegrator::printClassData..." << std::endl;
   os << "\nTimeRefinementIntegrator: this = "
      << (TimeRefinementIntegrator *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_integrator_time = " << d_integrator_time << "\n"
      << "d_start_time = " << d_start_time << "\n"
      << "d_end_time = " << d_end_time << "\n"
      << "d_integrator_step = " << d_step_level[0] << "\n"
      << "d_max_integrator_steps = " << d_max_steps_level[0] << "\n"
      << "d_grow_dt = " << d_grow_dt << std::endl;
   os << "d_just_regridded = " << d_just_regridded << std::endl;
   os << "d_last_finest_level = " << d_last_finest_level << std::endl;
   os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << std::endl;
   os << "d_refine_level_integrator = "
      << d_refine_level_integrator.get() << std::endl;
   os << "d_gridding_algorithm = "
      << d_gridding_algorithm.get() << std::endl;

   const int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
   for (int level_number = 0; level_number < max_levels; ++level_number) {
      printDataForLevel(os, level_number);
   }
}

/*
 *************************************************************************
 *
 * Print all level-specific data for TimeRefinementIntegrator object.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::printDataForLevel(
   std::ostream& os,
   const int level_number) const
{
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
   os << "\nTimeRefinementIntegrator::printDataForLevel..." << std::endl;
   os << "\nd_level_sim_time[" << level_number << "] = "
      << d_level_sim_time[level_number] << std::endl;
   os << "\nd_level_old_time[" << level_number << "] = "
      << d_level_old_time[level_number] << std::endl;
   os << "\nd_level_old_old_time[" << level_number << "] = "
      << d_level_old_old_time[level_number] << std::endl;
   os << "\nd_dt_max_level[" << level_number << "] = "
      << d_dt_max_level[level_number] << std::endl;
   os << "\nd_dt_actual_level[" << level_number << "] = "
      << d_dt_actual_level[level_number] << std::endl;
   os << "\nd_step_level[" << level_number << "] = "
      << d_step_level[level_number] << std::endl;
   os << "\nd_max_steps_level[" << level_number << "] = "
      << d_max_steps_level[level_number] << std::endl;
}

/*
 *************************************************************************
 *
 * Write the class version number and data members to restart database object.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION",
      ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION);

   restart_db->putDouble("start_time", d_start_time);
   restart_db->putDouble("end_time", d_end_time);
   restart_db->putDouble("grow_dt", d_grow_dt);
   restart_db->putInteger("max_integrator_steps", d_max_steps_level[0]);
   restart_db->putIntegerVector("regrid_interval", d_regrid_interval);
   restart_db->putIntegerVector("tag_buffer", d_tag_buffer);
   restart_db->putBool("DEV_barrier_and_time", d_barrier_and_time);
   restart_db->putDouble("d_integrator_time", d_integrator_time);
   restart_db->putInteger("d_integrator_step", d_step_level[0]);
   restart_db->putInteger("d_last_finest_level", d_last_finest_level);
   restart_db->putDoubleVector("d_dt_max_level", d_dt_max_level);
   restart_db->putDoubleVector("d_dt_actual_level", d_dt_actual_level);
}

/*
 *************************************************************************
 *
 * If simulation is not from restart, read in all data members from
 * the input database.  Otherwise, only override end_time, grow_dt
 * max_integrator_steps, and tag_buffer from the input database.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (!is_from_restart && !input_db) {
      TBOX_ERROR(": TimeRefinementIntegrator::getFromInput()\n"
         << "no input database supplied" << std::endl);
   }

   if (!is_from_restart) {
      /*
       * If not from restart, read in all data members from input database.
       */

      if (!d_use_refined_timestepping) {
         int regrid_interval =
            input_db->getIntegerWithDefault("regrid_interval", 1);
         if (!(regrid_interval >= 1)) {
            INPUT_RANGE_ERROR("regrid_interval");
         }
         setRegridInterval(regrid_interval);
      } else if (input_db->keyExists("regrid_interval")) {
         TBOX_WARNING("TimeRefinementIntegrator::getFromInput() warning...\n"
            << "regrid_interval input parameter not applicable with\n"
            << "refined timestepping and will be ignored." << std::endl);
      }

      d_start_time = input_db->getDouble("start_time");
      if (!(d_start_time >= 0)) {
         INPUT_RANGE_ERROR("start_time");
      }

      d_end_time = input_db->getDouble("end_time");
      if (!(d_end_time >= d_start_time)) {
         INPUT_RANGE_ERROR("end_time");
      }

      d_grow_dt = input_db->getDoubleWithDefault("grow_dt", 1.0);
      if (!(d_grow_dt > 0)) {
         INPUT_RANGE_ERROR("grow_dt");
      }

      d_max_steps_level[0] = input_db->getInteger("max_integrator_steps");
      if (!(d_max_steps_level[0] >= 0)) {
         INPUT_RANGE_ERROR("max_integrator_steps");
      }

      if (input_db->keyExists("tag_buffer")) {
         d_tag_buffer = input_db->getIntegerVector("tag_buffer");
         if (static_cast<int>(d_tag_buffer.size()) <
             (d_patch_hierarchy->getMaxNumberOfLevels() - 1)) {
            int tsize = static_cast<int>(d_tag_buffer.size());
            d_tag_buffer.resize(
               d_patch_hierarchy->getMaxNumberOfLevels() - 1,
               d_tag_buffer[tsize - 1]);
         }
      } else {
         int level_number;

         d_tag_buffer.resize(d_patch_hierarchy->getMaxNumberOfLevels());
         for (level_number = 0;
              level_number < d_patch_hierarchy->getMaxNumberOfLevels();
              ++level_number) {
            d_tag_buffer[level_number] = d_regrid_interval[level_number];
         }

         TBOX_WARNING("TimeRefinementIntegrator::getFromInput() warning...\n"
            << "Key data `tag_buffer' not found in input.  "
            << "Default values used.  See class header for details."
            << std::endl);
      }

      d_barrier_and_time =
         input_db->getBoolWithDefault("DEV_barrier_and_time", false);
   } else if (input_db) {
      bool read_on_restart =
         input_db->getBoolWithDefault("read_on_restart", false);

      if (read_on_restart) {
         if (!d_use_refined_timestepping) {
            int regrid_interval =
               input_db->getIntegerWithDefault("regrid_interval", 1);
            if (regrid_interval < 1) {
               TBOX_ERROR("TimeRefinementIntegrator::getFromInput() error...\n"
                  << "regrid_interval must be >=1." << std::endl);
            }
            setRegridInterval(regrid_interval);
         } else if (input_db->keyExists("regrid_interval")) {
            TBOX_WARNING("TimeRefinementIntegrator::getFromInput() warning...\n"
               << "regrid_interval input parameter not applicable with\n"
               << "refined timestepping and will be ignored." << std::endl);
         }

         if (input_db->keyExists("start_time")) {
            double tmp = input_db->getDouble("start_time");
            if (tmp != d_start_time) {
               TBOX_WARNING("TimeRefinementIntegrator::getFromInput warning...\n"
                  << "start_time may not be changed on restart." << std::endl);
            }
         }

         d_end_time = input_db->getDoubleWithDefault("end_time", d_end_time);
         if (d_end_time < d_start_time) {
            TBOX_ERROR("TimeRefinementIntegrator::getFromInput() error...\n"
               << "end_time must be >= start_time." << std::endl);
         }

         d_grow_dt = input_db->getDoubleWithDefault("grow_dt", d_grow_dt);
         if (d_grow_dt <= 0) {
            TBOX_ERROR("TimeRefinementIntegrator::getFromInput() error...\n"
               << "grow_dt must be > 0." << std::endl);
         }

         d_max_steps_level[0] =
            input_db->getIntegerWithDefault("max_integrator_steps",
               d_max_steps_level[0]);
         if (d_max_steps_level[0] < 0) {
            TBOX_ERROR("TimeRefinementIntegrator::getFromInput() error...\n"
               << "max_integrator_steps must be >= 0." << std::endl);
         } else if (d_max_steps_level[0] < d_step_level[0]) {
            TBOX_ERROR("TimeRefinementIntegrator::getFromInput() error...\n"
               << "max_integrator_steps must be >= current integrator step."
               << std::endl);
         }

         if (input_db->keyExists("tag_buffer")) {
            d_tag_buffer = input_db->getIntegerVector("tag_buffer");
            if (static_cast<int>(d_tag_buffer.size()) <
                (d_patch_hierarchy->getMaxNumberOfLevels() - 1)) {
               int tsize = static_cast<int>(d_tag_buffer.size());
               d_tag_buffer.resize(
                  d_patch_hierarchy->getMaxNumberOfLevels() - 1,
                  d_tag_buffer[tsize - 1]);
            }
         }

         d_barrier_and_time =
            input_db->getBoolWithDefault("DEV_barrier_and_time",
               d_barrier_and_time);
      }
   }
}

/*
 *************************************************************************
 *
 * Gets the database in the restart root database that corresponds to
 * the object name.  This method then checks to make sure that the class
 * version number and the restart version number are the same.  If they
 * are, then reads in the objects data members from the restart
 * database.
 *
 * Data read from restart database: d_start_time, d_end_time, d_grow_dt,
 * d_max_steps_level[0], d_regrid_interval, d_tag_buffer,
 * d_step_level[0], d_dt_max_level, d_dt_actual_level.
 *
 *************************************************************************
 */

void
TimeRefinementIntegrator::getFromRestart()
{

   std::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file." << std::endl);
   }
   std::shared_ptr<tbox::Database> db(
      restart_db->getDatabase(d_object_name));

   int ver = db->getInteger("ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION");
   if (ver != ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version."
                       << std::endl);
   }

   d_start_time = db->getDouble("start_time");
   d_end_time = db->getDouble("end_time");
   d_grow_dt = db->getDouble("grow_dt");
   d_max_steps_level[0] = db->getInteger("max_integrator_steps");
   d_regrid_interval = db->getIntegerVector("regrid_interval");
   d_tag_buffer = db->getIntegerVector("tag_buffer");
   d_barrier_and_time = db->getBool("DEV_barrier_and_time");
   d_integrator_time = db->getDouble("d_integrator_time");
   d_step_level[0] = db->getInteger("d_integrator_step");
   d_last_finest_level = db->getInteger("d_last_finest_level");
   d_dt_max_level = db->getDoubleVector("d_dt_max_level");
   d_dt_actual_level = db->getDoubleVector("d_dt_actual_level");
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TimeRefinementIntegrator::initializeCallback()
{
   t_initialize_hier = tbox::TimerManager::getManager()->
      getTimer("algs::TimeRefinementIntegrator::initializeHierarchy()");
   t_advance_hier = tbox::TimerManager::getManager()->
      getTimer("algs::TimeRefinementIntegrator::advanceHierarchy()");
   t_advance_level = tbox::TimerManager::getManager()->
      getTimer("algs::TimeRefinementIntegrator::advance_level");
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TimeRefinementIntegrator::finalizeCallback()
{
   t_initialize_hier.reset();
   t_advance_hier.reset();
   t_advance_level.reset();
}

}
}
