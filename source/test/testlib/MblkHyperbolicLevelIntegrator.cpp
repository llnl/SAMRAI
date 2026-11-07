/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Integration routines for single level in AMR hierarchy
 *                (basic hyperbolic systems)
 *
 ************************************************************************/

#include "MblkHyperbolicLevelIntegrator.h"

#include <stdlib.h>
#include <fstream>
#include <string>

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

//#define RECORD_STATS
#undef RECORD_STATS
#ifdef RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in flux
 * synchronization process between hierarchy levels.
 *
 *************************************************************************
 */

extern "C" {
// in upfluxsum.m4:

void SAMRAI_F77_FUNC(upfluxsumface2d0, UPFLUXSUMFACE2D0) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumface2d1, UPFLUXSUMFACE2D1) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumside2d0, UPFLUXSUMSIDE2D0) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumside2d1, UPFLUXSUMSIDE2D1) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, double *);

void SAMRAI_F77_FUNC(upfluxsumface3d0, UPFLUXSUMFACE3D0) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumface3d1, UPFLUXSUMFACE3D1) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumface3d2, UPFLUXSUMFACE3D2) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumside3d0, UPFLUXSUMSIDE3D0) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumside3d1, UPFLUXSUMSIDE3D1) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
void SAMRAI_F77_FUNC(upfluxsumside3d2, UPFLUXSUMSIDE3D2) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&,
   const double *, double *);
}

#define ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION (3)

using namespace SAMRAI;
using namespace algs;

/*
 *************************************************************************
 *
 * This constructor sets the HyperbolicPatchStrategy pointer and
 * initializes integration parameters to default values.  Communication
 * algorithms are created here too.  Other data members are read in
 * from the input database or from the restart database corresponding
 * to the specified object_name.
 *
 *************************************************************************
 */

MblkHyperbolicLevelIntegrator::MblkHyperbolicLevelIntegrator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<tbox::Database> input_db,
   MblkHyperbolicPatchStrategy* patch_strategy,
   const std::shared_ptr<hier::PatchHierarchy>& mblk_hierarchy,
   const bool use_time_refinement):
   d_patch_strategy(patch_strategy),
   d_object_name(object_name),
   d_dim(dim),
   d_use_time_refinement(use_time_refinement),
   d_cfl(tbox::MathUtilities<double>::getSignalingNaN()),
   d_cfl_init(tbox::MathUtilities<double>::getSignalingNaN()),
   d_lag_dt_computation(true),
   d_use_ghosts_for_dt(false),
   d_flux_is_face(true),
   d_flux_face_registered(false),
   d_flux_side_registered(false),
   d_mblk_bdry_fill_advance(new xfer::RefineAlgorithm()),
   d_mblk_bdry_fill_advance_new(new xfer::RefineAlgorithm()),
   d_mblk_bdry_fill_advance_old(new xfer::RefineAlgorithm()),
   d_mblk_coarsen_fluxsum(new xfer::CoarsenAlgorithm(dim)),
   d_mblk_coarsen_sync_data(new xfer::CoarsenAlgorithm(dim)),
   d_mblk_sync_initial_data(new xfer::CoarsenAlgorithm(dim)),
   d_coarsen_rich_extrap_init(new xfer::CoarsenAlgorithm(dim)),
   d_coarsen_rich_extrap_final(new xfer::CoarsenAlgorithm(dim)),
   d_mblk_fill_new_level(new xfer::RefineAlgorithm()),
   d_number_time_data_levels(2),
   d_scratch(hier::VariableDatabase::getDatabase()->getContext("SCRATCH")),
   d_current(hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
   d_new(hier::VariableDatabase::getDatabase()->getContext("NEW")),
   d_plot_context(
      hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
   d_have_flux_on_level_zero(false),
   d_distinguish_mpi_reduction_costs(false),
   t_advance_bdry_fill_comm(tbox::TimerManager::getManager()->
                            getTimer("algs::MblkHyperbolicLevelIntegrator::advance_bdry_fill_comm")),
   t_error_bdry_fill_create(tbox::TimerManager::getManager()->
                            getTimer("algs::MblkHyperbolicLevelIntegrator::error_bdry_fill_create")),
   t_error_bdry_fill_comm(tbox::TimerManager::getManager()->
                          getTimer("algs::MblkHyperbolicLevelIntegrator::error_bdry_fill_comm")),
   t_mpi_reductions(tbox::TimerManager::getManager()->
                    getTimer("algs::MblkHyperbolicLevelIntegrator::mpi_reductions")),
   t_initialize_level_data(tbox::TimerManager::getManager()->
                           getTimer("algs::MblkHyperbolicLevelIntegrator::initializeLevelData()")),
   t_fill_new_level_create(tbox::TimerManager::getManager()->
                           getTimer("algs::MblkHyperbolicLevelIntegrator::fill_new_level_create")),
   t_fill_new_level_comm(tbox::TimerManager::getManager()->
                         getTimer("algs::MblkHyperbolicLevelIntegrator::fill_new_level_comm")),
   t_advance_bdry_fill_create(tbox::TimerManager::getManager()->
                              getTimer(
                                 "algs::MblkHyperbolicLevelIntegrator::advance_bdry_fill_create")),
   t_new_advance_bdry_fill_create(tbox::TimerManager::getManager()->
                                  getTimer(
                                     "algs::MblkHyperbolicLevelIntegrator::new_advance_bdry_fill_create")),
   t_apply_gradient_detector(tbox::TimerManager::getManager()->
                             getTimer(
                                "algs::MblkHyperbolicLevelIntegrator::applyGradientDetector()")),
   t_coarsen_rich_extrap(tbox::TimerManager::getManager()->
                         getTimer("algs::MblkHyperbolicLevelIntegrator::coarsen_rich_extrap")),
   t_get_level_dt(tbox::TimerManager::getManager()->
                  getTimer("algs::MblkHyperbolicLevelIntegrator::getLevelDt()")),
   t_get_level_dt_sync(tbox::TimerManager::getManager()->
                       getTimer("algs::MblkHyperbolicLevelIntegrator::getLevelDt()_sync")),
   t_advance_level(tbox::TimerManager::getManager()->
                   getTimer("algs::MblkHyperbolicLevelIntegrator::advanceLevel()")),
   t_new_advance_bdry_fill_comm(tbox::TimerManager::getManager()->
                                getTimer(
                                   "algs::MblkHyperbolicLevelIntegrator::new_advance_bdry_fill_comm")),
   t_patch_num_kernel(tbox::TimerManager::getManager()->
                      getTimer("algs::MblkHyperbolicLevelIntegrator::patch_numerical_kernels")),
   t_advance_level_sync(tbox::TimerManager::getManager()->
                        getTimer("algs::MblkHyperbolicLevelIntegrator::advanceLevel()_sync")),
   t_std_level_sync(tbox::TimerManager::getManager()->
                    getTimer(
                       "algs::MblkHyperbolicLevelIntegrator::standardLevelSynchronization()")),
   t_sync_new_levels(tbox::TimerManager::getManager()->
                     getTimer("algs::MblkHyperbolicLevelIntegrator::synchronizeNewLevels()"))
{
   NULL_USE(mblk_hierarchy);

   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(patch_strategy != 0);

   tbox::RestartManager::getManager()->
   registerRestartItem(d_object_name, this);

   /*
    * Initialize object with data read from the input and restart databases.
    */

   bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, from_restart);
}

/*
 *************************************************************************
 *
 * Destructor tells the tbox::RestartManager to remove this object from
 * the list of restart items.
 *
 *************************************************************************
 */
MblkHyperbolicLevelIntegrator::~MblkHyperbolicLevelIntegrator()
{
   tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/*
 *************************************************************************
 *
 * Initialize integration data on all patches on level.  This process
 * is used at the start of the simulation to set the initial hierarchy
 * data and after adaptive regridding.  In the second case, the old
 * level pointer points to the level that existed in the hierarchy
 * before regridding.  This pointer may be null, in which case it is
 * ignored.  If it is non-null, then data is copied from the old level
 * to the new level before the old level is discarded.
 *
 * Note that we also allocate flux storage for the coarsest AMR
 * hierarchy level here (i.e., level 0).  The time step sequence on
 * level 0 is dictated by the user code; so to avoid any memory
 * management errors, flux storage on level 0 persists as long as the
 * level does.
 *
 *************************************************************************
 */
void MblkHyperbolicLevelIntegrator::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0) &&
      (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(!old_level || (level_number == old_level->getLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

   t_initialize_level_data->start();

   std::shared_ptr<hier::PatchLevel> mblk_level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Allocate storage needed to initialize level and fill data
    * from coarser levels in AMR hierarchy, potentially. Since
    * time gets set when we allocate data, re-stamp it to current
    * time if we don't need to allocate.
    */
   if (allocate_data) {
      mblk_level->allocatePatchData(d_new_patch_init_data, init_data_time);
      mblk_level->allocatePatchData(d_old_time_dep_data, init_data_time);
   } else {
      mblk_level->setTime(init_data_time, d_new_patch_init_data);
   }

   /*
    * Create schedules for filling new level and fill data.
    */

   if ((level_number > 0) || old_level) {
      t_fill_new_level_create->start();

      std::shared_ptr<xfer::RefineSchedule> sched(
         d_mblk_fill_new_level->createSchedule(mblk_level,
            old_level,
            level_number - 1,
            hierarchy,
            d_patch_strategy));
      t_fill_new_level_create->stop();

      d_patch_strategy->setDataContext(d_scratch);

      t_fill_new_level_comm->start();
      sched->fillData(init_data_time);
      t_fill_new_level_comm->stop();

      d_patch_strategy->clearDataContext();
   }

   if ((d_number_time_data_levels == 3) && can_be_refined) {

      hier::VariableDatabase* variable_db =
         hier::VariableDatabase::getDatabase();

      for (hier::PatchLevel::iterator mi(mblk_level->begin());
           mi != mblk_level->end(); ++mi) {

         std::list<std::shared_ptr<hier::Variable> >::iterator time_dep_var =
            d_time_dep_variables.begin();
         while (time_dep_var != d_time_dep_variables.end()) {
            int old_indx =
               variable_db->mapVariableAndContextToIndex(*time_dep_var,
                  d_old);
            int cur_indx =
               variable_db->mapVariableAndContextToIndex(*time_dep_var,
                  d_current);

            (*mi)->setPatchData(old_indx, (*mi)->getPatchData(cur_indx));

            ++time_dep_var;
         }
      }

   } // loop over patches

   /*
    * Initialize data on patch interiors.
    */
   d_patch_strategy->setDataContext(d_current);

   for (hier::PatchLevel::iterator mi(mblk_level->begin());
        mi != mblk_level->end(); ++mi) {

      (*mi)->allocatePatchData(d_temp_var_scratch_data, init_data_time);

      d_patch_strategy->initializeDataOnPatch(**mi,
         init_data_time,
         initial_time);

      (*mi)->deallocatePatchData(d_temp_var_scratch_data);

   } // loop over patches

   d_patch_strategy->clearDataContext();

   //d_mblk_fill_new_level.reset();

   t_initialize_level_data->stop();
}

/*
 *************************************************************************
 *
 * Reset hierarchy configuration information where the range of new
 * hierarchy levels is specified.   The information updated involves
 * the cached communication schedules maintained by the algorithm.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::resetHierarchyConfiguration(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0) &&
      (coarsest_level <= finest_level) &&
      (finest_level <= hierarchy->getFinestLevelNumber()));
   for (int ln0 = 0; ln0 <= finest_level; ++ln0) {
      TBOX_ASSERT(hierarchy->getPatchLevel(ln0));
   }
#else
   NULL_USE(finest_level);
#endif

   int finest_hiera_level = hierarchy->getFinestLevelNumber();

   d_mblk_bdry_sched_advance.resize(finest_hiera_level + 1);
   d_mblk_bdry_sched_advance_new.resize(finest_hiera_level + 1);

   for (int ln = coarsest_level; ln <= finest_hiera_level; ++ln) {
      std::shared_ptr<hier::PatchLevel> mblk_level(
         hierarchy->getPatchLevel(ln));

      t_advance_bdry_fill_create->start();
      d_mblk_bdry_sched_advance[ln] =
         d_mblk_bdry_fill_advance->createSchedule(mblk_level,
            ln - 1,
            hierarchy,
            d_patch_strategy);
      t_advance_bdry_fill_create->stop();

      if (!d_lag_dt_computation && d_use_ghosts_for_dt) {
         t_new_advance_bdry_fill_create->start();
         d_mblk_bdry_sched_advance_new[ln] =
            d_mblk_bdry_fill_advance_new->createSchedule(mblk_level,
               ln - 1,
               hierarchy,
               d_patch_strategy);
         t_new_advance_bdry_fill_create->stop();
      }

   }

}

/*
 *************************************************************************
 *
 * Call patch routines to tag cells near large gradients.
 * These cells will be refined.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

   t_apply_gradient_detector->start();

   std::shared_ptr<hier::PatchLevel> mblk_level(
      hierarchy->getPatchLevel(level_number));

   mblk_level->allocatePatchData(d_saved_var_scratch_data, error_data_time);

   d_patch_strategy->setDataContext(d_scratch);

   t_error_bdry_fill_comm->start();
   d_mblk_bdry_sched_advance[level_number]->fillData(error_data_time);
   t_error_bdry_fill_comm->stop();

   for (hier::PatchLevel::iterator mi(mblk_level->begin());
        mi != mblk_level->end(); ++mi) {

      d_patch_strategy->
      tagGradientDetectorCells(**mi,
         error_data_time,
         initial_time,
         tag_index,
         uses_richardson_extrapolation_too);

   } // loop over patches

   d_patch_strategy->clearDataContext();

   mblk_level->deallocatePatchData(d_saved_var_scratch_data);

   t_apply_gradient_detector->stop();

}

/*
 *************************************************************************
 *
 * The Richardson extrapolation algorithm requires a coarsened version
 * of the level on which error estiamtion is performed.  This routine
 * is used to coarsen data from a level in the AMR hierarchy to some
 * coarsened version of it.  Note that this routine will be called twice
 * The init_coarse_level boolean argument indicates whether data is
 * set on the coarse level by coarsening the "old" time level solution
 * or by coarsening the "new" solution on the fine level (i.e., after
 * it has been advanced).
 *
 * The contexts used for coarsening old data depends on the number of
 * time levels.  We always want to use data at the oldest time on the
 * fine level, coarsened to the CURRENT context on the coarse level.
 * Thus, if the problem uses two time levels, we coarsen data from
 * CURRENT on fine level (since CURRENT is the oldest time maintained)
 * to CURRENT on the coarse level.  If the problem uses three time
 * levels, we coarsen from OLD on the fine level (since OLD is the
 * time maintained) to CURRENT on the coarse level.
 *
 * When the boolean is false, indicating we are operating at the new
 * time, we coarsen the time advanced solution at the NEW context on
 * the fine level to the NEW context on the coarse level so that they
 * may be compared later.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::coarsenDataForRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::PatchLevel>& coarse_level,
   const double coarsen_data_time,
   const bool before_advance)
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(hierarchy);
   NULL_USE(level_number);
#endif

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(coarse_level);

   t_coarsen_rich_extrap->start();

//   std::shared_ptr<hier::PatchLevel> level(
//      hierarchy->getPatchLevel(level_number));

   if (before_advance) {

      coarse_level->allocatePatchData(d_new_patch_init_data,
         coarsen_data_time);

      if (d_number_time_data_levels == 3) {
         d_patch_strategy->setDataContext(d_old);
      } else {
         d_patch_strategy->setDataContext(d_current);
      }

      TBOX_ERROR("Incomplete DLBG code.");
//      d_coarsen_rich_extrap_init->
//         createSchedule(coarse_level, level, d_patch_strategy)->
//            coarsenData();

      d_patch_strategy->clearDataContext();

   } else {

      coarse_level->allocatePatchData(d_new_time_dep_data,
         coarsen_data_time);

      d_patch_strategy->setDataContext(d_new);

      TBOX_ERROR("Incomplete DLBG code.");
//      d_coarsen_rich_extrap_final->
//         createSchedule(coarse_level, level, d_patch_strategy)->
//         coarsenData();

      d_patch_strategy->clearDataContext();

   }

   t_coarsen_rich_extrap->stop();

}

/*
 *************************************************************************
 *
 * Call patch routines to tag cells for refinement using Richardson
 * extrapolation.    Richardson extrapolation requires two copies of
 * the solution to compare.  The NEW context holds the solution
 * computed on the fine level and coarsened, whereas the CURRENT
 * context holds the solution integrated on the coarse level after
 * coarsening the initial data from the fine level.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::applyRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double error_data_time,
   const int tag_index,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_time,
   const bool uses_gradient_detector_too)
{
   TBOX_ASSERT(level);

   /*
    * Compare solutions computed on level (stored in NEW context) and on
    * the coarser level (stored in CURR context) on the patches of the
    * coarser level.  The patch strategy implements the compare operations
    * performed on each patch.
    */

   int error_level_number =
      level->getNextCoarserHierarchyLevelNumber() + 1;

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      d_patch_strategy->
      tagRichardsonExtrapolationCells(*patch,
         error_level_number,
         d_new,                                     //  finer context
         d_current,                                 //  coarser context
         error_data_time,
         deltat,
         error_coarsen_ratio,
         initial_time,
         tag_index,
         uses_gradient_detector_too);
   }

}

/*
 *************************************************************************
 *
 * Initialize level integrator by:
 *
 *   (1) Setting the number of time data levels based on needs of
 *       the gridding algorithm
 *   (2) Invoking variable registration in patch strategy.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::initializeLevelIntegrator(
   const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg)
{
   TBOX_ASSERT(gridding_alg);

   d_number_time_data_levels = 2;

   if ((gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() < 1) ||
       (gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() > 3)) {
      TBOX_ERROR("MblkHyperbolicLevelIntegrator::initializeLevelIntegrator "
         << "error...\n" << "   object name = " << d_object_name
         << "   gridding algorithm has bad error coarsen ratio" << std::endl);
   }

   if ((gridding_alg->getTagAndInitializeStrategy()->everUsesTimeIntegration()) &&
       (gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() == 3)) {
      d_number_time_data_levels = 3;
      d_old = hier::VariableDatabase::getDatabase()->getContext("OLD");
   }

   d_patch_strategy->registerModelVariables(this);

   d_patch_strategy->setupLoadBalancer(this,
      gridding_alg.get());
}

/*
 *************************************************************************
 *
 * Invoke dt calculation routines in patch strategy and take a min
 * over all patches on the level.  The result will be the max of the
 * next timestep on the level. If the boolean recompute_dt is true,
 * the max timestep on the level will be computed.  If it is false,
 * the method will simply access the latest dt stored in the time
 * refinement integrator.
 *
 *************************************************************************
 */

double
MblkHyperbolicLevelIntegrator::getLevelDt(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double dt_time,
   const bool initial_time)
{
   TBOX_ASSERT(level);

   t_get_level_dt->start();

   double dt = tbox::MathUtilities<double>::getMax();

   if (!d_use_ghosts_for_dt) {

      d_patch_strategy->setDataContext(d_current);

      for (hier::PatchLevel::iterator mi(level->begin());
           mi != level->end(); ++mi) {

         (*mi)->allocatePatchData(d_temp_var_scratch_data, dt_time);

         double patch_dt;
         patch_dt = d_patch_strategy->
            computeStableDtOnPatch(**mi,
               initial_time,
               dt_time);

         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);

         (*mi)->deallocatePatchData(d_temp_var_scratch_data);

      } // loop over patches

      d_patch_strategy->clearDataContext();

   } else {

      level->allocatePatchData(d_saved_var_scratch_data, dt_time);

      d_patch_strategy->setDataContext(d_scratch);

      t_advance_bdry_fill_comm->start();
      d_mblk_bdry_sched_advance[level->getLevelNumber()]->fillData(dt_time);
      t_advance_bdry_fill_comm->stop();

      for (hier::PatchLevel::iterator mi(level->begin());
           mi != level->end(); ++mi) {

         (*mi)->allocatePatchData(d_temp_var_scratch_data, dt_time);

         double patch_dt;
         patch_dt = d_patch_strategy->
            computeStableDtOnPatch(**mi,
               initial_time,
               dt_time);

         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);

         (*mi)->deallocatePatchData(d_temp_var_scratch_data);

      } // loop over patches

      d_patch_strategy->clearDataContext();

      /*
       * Copy data from scratch to current and de-allocate scratch storage.
       * This may be excessive here, but seems necessary if the
       * computation of dt affects the state of the problem solution.
       * Also, this getLevelDt() routine is called at initialization only
       * in most cases.
       */

      copyTimeDependentData(level, d_scratch, d_current);

      level->deallocatePatchData(d_saved_var_scratch_data);

   }

   t_get_level_dt_sync->start();

   if (d_distinguish_mpi_reduction_costs) {
      tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
      t_get_level_dt_sync->stop();
      t_mpi_reductions->start();
   }

   /*
    * The level time increment is a global min over all patches.
    */

   double global_dt = dt;
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_dt, 1, MPI_MIN);
   }
   global_dt *= tbox::MathUtilities<double>::Min(d_cfl_init, d_cfl);

   if (d_distinguish_mpi_reduction_costs) {
      t_mpi_reductions->stop();
   } else {
      t_get_level_dt_sync->stop();
   }

   t_get_level_dt->stop();

   return global_dt;

}

/*
 *************************************************************************
 *
 * For the standard explicit integration algorithm for hyperbolic
 * conservation laws, the fine time increment is the coarse increment
 * divided by the maximum mesh ratio (independent of level number).
 *
 *************************************************************************
 */

double
MblkHyperbolicLevelIntegrator::getMaxFinerLevelDt(
   const int finer_level_number,
   const double coarse_dt,
   const hier::IntVector& ratio)
{
   NULL_USE(finer_level_number);
   TBOX_ASSERT(ratio.min() > 0);
   return coarse_dt / static_cast<double>(ratio.max());
}

/*
 *************************************************************************
 *
 * Integrate data on all patches in patch level from current time
 * to new time (new_time) using a single time step.  Before the advance
 * can occur, proper ghost cell information is obtained for all patches
 * on the level.  Then, local patches are advanced sequentially in the
 * loop over patches.  The details of the routine are as follows:
 *
 *  0) Allocate storage for new time level data. Also, allocate
 *     necessary FLUX and flux integral storage if needed
 *     (i.e., if regrid_advance is false, first_step is true, and
 *     coarser or finer level than current level exists in hierarchy.)
 *
 *  1) Scratch space is filled so that, for each patch, interior data
 *     and ghost cell bdry data correspond to specified time.
 *
 *  1a) Call user routines to pre-process advance data, if needed.
 *
 *  2) Compute explicit fluxes in scratch space using data on
 *     patch + ghosts at given time.
 *
 *  3) Apply conservative difference in scratch space to advance patch
 *     interior data to time = new_time.
 *
 *  3a) Call user routines to post-process advance data, if needed.
 *
 *  4) Compute next stable time increment for subsequent level advances:
 *
 *     4a) If (d_lag_dt_computation == true) {
 *            DO NOT RECOMPUTE characteristic data after advancing
 *            data on patch. Use characteristic data corresponding
 *            to current time level, computed prior to flux computation,
 *            in dt calculation.
 *            If (d_use_ghosts_for_dt == true)
 *               - Compute dt using data on patch+ghosts at time.
 *            Else
 *               - Compute dt using data on patch interior ONLY.
 *         }
 *
 *     4b) Copy data from scratch space patch interior to new data
 *         storage for patch (i.e., at time = new_time).
 *
 *     4a) If (d_lag_dt_computation == false) {
 *            RECOMPUTE characteristic data after advancing data on
 *            patch. Use characteristic data corresponding to new time
 *            level in dt calculation.
 *            If (d_use_ghosts_for_dt == true)
 *               - Refill scratch space with new interior patch data
 *                 and ghost cell bdry data correspond to new time.
 *                 (NOTE: This requires a new boundary schedule.)
 *               - Compute dt using data on patch+ghosts at new_time.
 *            Else
 *               - Compute dt using data on patch interior ONLY.
 *                 (using patch interior data at new_time)
 *         }
 *
 *  5) If (ln > 0), update flux integrals by adding patch bdry FLUXes
 *     to flux sums.
 *
 * Important Notes:
 *    1) In order to advance finer levels (if they exist), both old
 *       and new data for each patch on the level must be maintained.
 *    2) If the timestep is the first in the timestep loop on the level
 *       (indicated by first_step), then time interpolation is
 *       is unnecessary to fill ghost cells from the next coarser level.
 *    3) The new dt is not calculated if regrid_advance is true.
 *       If this is the case, it is assumed that the results of the
 *       advance and the timestep calculation will be discarded
 *       (e.g., during regridding, or initialization).  Also, allocation
 *       and post-processing of FLUX/flux integral data is not performed
 *       in this case.
 *
 *************************************************************************
 */

double
MblkHyperbolicLevelIntegrator::advanceLevel(
   const std::shared_ptr<hier::PatchLevel>& level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const double current_time,
   const double new_time,
   const bool first_step,
   const bool last_step,
   const bool regrid_advance)
{
   TBOX_ASSERT(level);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(current_time <= new_time);

#ifdef RECORD_STATS
   std::shared_ptr<tbox::Statistic> num_boxes_l0(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberBoxesL0", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_boxes_l1(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberBoxesL1", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_boxes_l2(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberBoxesL2", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_boxes_l3(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberBoxesL3", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_gridcells_l0(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberGridcellsL0", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_gridcells_l1(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberGridcellsL1", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_gridcells_l2(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberGridcellsL2", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> num_gridcells_l3(
      tbox::Statistician::getStatistician()->
      getStatistic("NumberGridcellsL3", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> timestamp_l0(
      tbox::Statistician::getStatistician()->
      getStatistic("TimeStampL0", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> timestamp_l1(
      tbox::Statistician::getStatistician()->
      getStatistic("TimeStampL1", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> timestamp_l2(
      tbox::Statistician::getStatistician()->
      getStatistic("TimeStampL2", "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> timestamp_l3(
      tbox::Statistician::getStatistician()->
      getStatistic("TimeStampL3", "PROC_STAT"));

   int level_num = level->getLevelNumber();

   /*
    * Record number of gridcells on each patch.  Note that patch
    * stat requires a seq number to be identified.
    */
   double level_gridcells = 0.;
   double level_local_patches = 0.;
   // to count total gridcells on level
   //hier::BoxList boxes = level->getBoxes();
   //for (hier::BoxList::Iterator i(boxes); i; ++i) {
   //   level_gridcells += itr().size();
   //}
   // to count gridcells on this processor

   for (hier::PatchLevel::iterator mi(level->begin());
        mi != level->end(); ++mi) {
      level_gridcells += (*mi)->getBox().size();
      level_local_patches += 1.0;
   } // loop over patches

   if (level_num == 0) {
      num_boxes_l0->recordProcStat(level_local_patches);
      num_gridcells_l0->recordProcStat(level_gridcells);
      timestamp_l0->recordProcStat(current_time);
   }
   if (level_num == 1) {
      num_boxes_l1->recordProcStat(level_local_patches);
      num_gridcells_l1->recordProcStat(level_gridcells);
      timestamp_l1->recordProcStat(current_time);
   }
   if (level_num == 2) {
      num_boxes_l2->recordProcStat(level_local_patches);
      num_gridcells_l2->recordProcStat(level_gridcells);
      timestamp_l2->recordProcStat(current_time);
   }
   if (level_num == 3) {
      num_boxes_l3->recordProcStat(level_local_patches);
      num_gridcells_l3->recordProcStat(level_gridcells);
      timestamp_l3->recordProcStat(current_time);
   }
#endif

   t_advance_level->start();

   const int level_number = level->getLevelNumber();
   const double dt = new_time - current_time;

   /*
    * (1) Allocate data needed for advancing level.
    * (2) Generate temporary communication schedule to fill ghost
    *     cells, if needed.
    * (3) Fill ghost cell data.
    * (4) Process flux storage before the advance.
    */

   level->allocatePatchData(d_new_time_dep_data, new_time);
   level->allocatePatchData(d_saved_var_scratch_data, current_time);

   std::shared_ptr<xfer::RefineSchedule> mblk_fill_schedule;

   const bool in_hierarchy = level->inHierarchy();

   if (!in_hierarchy) {
      t_error_bdry_fill_create->start();
      if (d_number_time_data_levels == 3) {
         mblk_fill_schedule = d_mblk_bdry_fill_advance_old->
            createSchedule(level,
               level->getLevelNumber() - 1,
               hierarchy,
               d_patch_strategy);
      } else {
         mblk_fill_schedule = d_mblk_bdry_fill_advance->
            createSchedule(level,
               level->getLevelNumber() - 1,
               hierarchy,
               d_patch_strategy);
      }
      t_error_bdry_fill_create->stop();
   } else {
      mblk_fill_schedule = d_mblk_bdry_sched_advance[level_number];
   }

   d_patch_strategy->setDataContext(d_scratch);
   if (regrid_advance) {
      t_error_bdry_fill_comm->start();
   } else {
      t_advance_bdry_fill_comm->start();
   }
   mblk_fill_schedule->fillData(current_time);
   if (regrid_advance) {
      t_error_bdry_fill_comm->stop();
   } else {
      t_advance_bdry_fill_comm->stop();
   }

   d_patch_strategy->clearDataContext();
   mblk_fill_schedule.reset();

   preprocessFluxData(level,
      current_time,
      new_time,
      regrid_advance,
      first_step,
      last_step);

   /*
    * (5) Call user-routine to pre-process state data, if needed.
    * (6) Advance solution on all level patches (scratch storage).
    * (7) Copy new solution to from scratch to new storage.
    * (8) Call user-routine to post-process state data, if needed.
    */
   t_patch_num_kernel->start();
   d_patch_strategy->preprocessAdvanceLevelState(level,
      current_time,
      dt,
      first_step,
      last_step,
      regrid_advance);
   t_patch_num_kernel->stop();

   d_patch_strategy->setDataContext(d_scratch);
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      patch->allocatePatchData(d_temp_var_scratch_data, current_time);

      t_patch_num_kernel->start();
      d_patch_strategy->computeFluxesOnPatch(*patch,
         current_time,
         dt);
      t_patch_num_kernel->stop();

      bool at_syncronization = false;

      t_patch_num_kernel->start();
      d_patch_strategy->conservativeDifferenceOnPatch(*patch,
         current_time,
         dt,
         at_syncronization);
      t_patch_num_kernel->stop();

      patch->deallocatePatchData(d_temp_var_scratch_data);
   }
   d_patch_strategy->clearDataContext();

   level->setTime(new_time, d_saved_var_scratch_data);
   level->setTime(new_time, d_flux_var_data);

   copyTimeDependentData(level, d_scratch, d_new);

   t_patch_num_kernel->start();
   d_patch_strategy->postprocessAdvanceLevelState(level,
      current_time,
      dt,
      first_step,
      last_step,
      regrid_advance);
   t_patch_num_kernel->stop();

   /*
    * (9) If the level advance is for regridding, we compute the next timestep:
    *
    * (a) If the dt computation is lagged (i.e., we use pre-advance data
    *     to compute timestep), we reset scratch space on patch interiors
    *     if needed.  Then, we set the strategy context to current or scratch
    *     depending on whether ghost values are used to compute dt.
    * (b) If the dt computation is not lagged (i.e., we use advanced data
    *     to compute timestep), we refill scratch space, including ghost
    *     data with new solution values if needed.  Then, we set the strategy
    *     context to new or scratch depending on whether ghost values are
    *     used to compute dt.
    * (c) Then, we loop over patches and compute the dt on each patch.
    */

   double dt_next = tbox::MathUtilities<double>::getMax();

   if (!regrid_advance) {

      if (d_lag_dt_computation) {

         if (d_use_ghosts_for_dt) {
            d_patch_strategy->setDataContext(d_scratch);
            copyTimeDependentData(level, d_current, d_scratch);
         } else {
            d_patch_strategy->setDataContext(d_current);
         }

      } else {

         if (d_use_ghosts_for_dt) {

            if (!d_mblk_bdry_sched_advance_new[level_number]) {
               TBOX_ERROR(
                  d_object_name << ":  "
                                << "Attempt to fill new ghost data for timestep"
                                << "computation, but schedule not defined." << std::endl);
            }

            d_patch_strategy->setDataContext(d_scratch);
            t_new_advance_bdry_fill_comm->start();
            d_mblk_bdry_sched_advance_new[level_number]->fillData(new_time);
            t_new_advance_bdry_fill_comm->stop();

         } else {
            d_patch_strategy->setDataContext(d_new);
         }

      }

      for (hier::PatchLevel::iterator mi(level->begin());
           mi != level->end(); ++mi) {

         (*mi)->allocatePatchData(d_temp_var_scratch_data, new_time);
         // "false" argument indicates "initial_time" is false.
         t_patch_num_kernel->start();
         double patch_dt =
            d_patch_strategy->computeStableDtOnPatch(**mi,
               false,
               new_time);
         t_patch_num_kernel->stop();

         dt_next = tbox::MathUtilities<double>::Min(dt_next, patch_dt);

         (*mi)->deallocatePatchData(d_temp_var_scratch_data);

      } // loop over patches

      d_patch_strategy->clearDataContext();

   } // !regrid_advance

   level->deallocatePatchData(d_saved_var_scratch_data);

   postprocessFluxData(level,
      regrid_advance,
      first_step,
      last_step);

   t_advance_level->stop();

   t_advance_level_sync->start();

   if (d_distinguish_mpi_reduction_costs) {
      tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
      t_advance_level_sync->stop();
      t_mpi_reductions->start();
   }

   double next_dt = dt_next;
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&next_dt, 1, MPI_MIN);
   }
   next_dt *= d_cfl;

   if (d_distinguish_mpi_reduction_costs) {
      t_mpi_reductions->stop();
   } else {
      t_advance_level_sync->stop();
   }

   return next_dt;
}

/*
 *************************************************************************
 *
 * Synchronize data between patch levels according to the standard
 * hyperbolic flux correction algorithm.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::standardLevelSynchronization(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level,
   const double sync_time,
   const double old_time)
{
   std::vector<double> old_times(finest_level - coarsest_level + 1);
   for (int i = coarsest_level; i <= finest_level; ++i) {
      old_times[i] = old_time;
   }
   standardLevelSynchronization(hierarchy, coarsest_level, finest_level,
      sync_time, old_times);
}

void
MblkHyperbolicLevelIntegrator::standardLevelSynchronization(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level,
   const double sync_time,
   const std::vector<double>& old_times)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level < finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(static_cast<int>(old_times.size()) >= finest_level);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int ln = coarsest_level; ln < finest_level; ++ln) {
      TBOX_ASSERT(hierarchy->getPatchLevel(ln));
      TBOX_ASSERT(sync_time >= old_times[ln]);
   }
#endif
   TBOX_ASSERT(hierarchy->getPatchLevel(finest_level));

   t_std_level_sync->start();

   for (int fine_ln = finest_level; fine_ln > coarsest_level; --fine_ln) {
      const int coarse_ln = fine_ln - 1;

      TBOX_ASSERT(sync_time >= old_times[coarse_ln]);

      std::shared_ptr<hier::PatchLevel> mblk_fine_level(
         hierarchy->getPatchLevel(fine_ln));
      std::shared_ptr<hier::PatchLevel> mblk_coarse_level(
         hierarchy->getPatchLevel(coarse_ln));

      synchronizeLevelWithCoarser(mblk_fine_level,
         mblk_coarse_level,
         sync_time,
         old_times[coarse_ln]);

      mblk_fine_level->deallocatePatchData(d_fluxsum_data);
      mblk_fine_level->deallocatePatchData(d_flux_var_data);

      if (coarse_ln > coarsest_level) {
         mblk_coarse_level->deallocatePatchData(d_flux_var_data);
      } else {
         if (coarsest_level == 0) {
            mblk_coarse_level->deallocatePatchData(d_flux_var_data);
            d_have_flux_on_level_zero = false;
         }
      }

   }

   t_std_level_sync->stop();

}

/*
 *************************************************************************
 *
 * Coarsen current solution data from finest hierarchy level specified
 * down through the coarsest hierarchy level specified, if initial_time
 * is true (i.e., hierarchy is being constructed at initial simulation
 * time).  After data is coarsened, the user's initialization routine
 * is called to reset data (as needed by the application) before
 * that solution is further coarsened to the next coarser level in the
 * hierarchy.  If initial_time is false, then this routine does nothing
 * In that case, interpolation of data from coarser levels is sufficient
 * to set data on new levels in the hierarchy during regridding.
 *
 * NOTE: The fact that this routine does nothing when called at any
 *       time later than when the AMR hierarchy is constructed initially
 *        may need to change at some point based on application needs.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::synchronizeNewLevels(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level,
   const double sync_time,
   const bool initial_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level < finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int ln = coarsest_level; ln <= finest_level; ++ln) {
      TBOX_ASSERT(hierarchy->getPatchLevel(ln));
   }
#endif

   std::shared_ptr<tbox::Timer> t_sync_initial_create(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::sync_initial_create"));
   std::shared_ptr<tbox::Timer> t_sync_initial_comm(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::sync_initial_comm"));

   t_sync_new_levels->start();

   if (initial_time) {

      d_patch_strategy->setDataContext(d_current);

      for (int fine_ln = finest_level; fine_ln > coarsest_level; --fine_ln) {
         const int coarse_ln = fine_ln - 1;

         std::shared_ptr<hier::PatchLevel> fine_level(
            hierarchy->getPatchLevel(fine_ln));

         std::shared_ptr<hier::PatchLevel> coarse_level(
            hierarchy->getPatchLevel(coarse_ln));

         if (d_do_coarsening) {
            t_sync_initial_create->start();
            std::shared_ptr<xfer::CoarsenSchedule> sched(
               d_mblk_sync_initial_data->createSchedule(coarse_level,
                  fine_level,
                  d_patch_strategy));
            t_sync_initial_create->stop();

            t_sync_initial_comm->start();
            sched->coarsenData();
            t_sync_initial_comm->stop();
         }

         for (hier::PatchLevel::iterator mi(coarse_level->begin());
              mi != coarse_level->end(); ++mi) {

            (*mi)->allocatePatchData(d_temp_var_scratch_data, sync_time);

            d_patch_strategy->initializeDataOnPatch(**mi,
               sync_time,
               initial_time);
            (*mi)->deallocatePatchData(d_temp_var_scratch_data);

         }
      }

      d_patch_strategy->clearDataContext();

   } // if (initial_time)

   t_sync_new_levels->stop();

}

/*
 *************************************************************************
 *
 * Synchronize data between coarse and fine patch levels according to
 * the standard hyperbolic flux correction algorithm.  The steps of
 * the algorithm are:
 *
 *    (1) Replace coarse time-space flux integrals at coarse-fine
 *        boundaries with time-space flux integrals computed on fine
 *        level.
 *    (2) Repeat conservative difference on coarse level with corrected
 *        fluxes.
 *    (3) Conservatively coarsen solution on interior of fine level to
 *        coarse level.
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::synchronizeLevelWithCoarser(
   const std::shared_ptr<hier::PatchLevel> mblk_fine_level,
   const std::shared_ptr<hier::PatchLevel> mblk_coarse_level,
   const double sync_time,
   const double coarse_sim_time)
{
   TBOX_ASSERT(mblk_fine_level);
   TBOX_ASSERT(mblk_coarse_level);
   TBOX_ASSERT(mblk_coarse_level->getLevelNumber() ==
      (mblk_fine_level->getLevelNumber() - 1));

   std::shared_ptr<tbox::Timer> t_coarsen_fluxsum_create(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::coarsen_fluxsum_create"));
   std::shared_ptr<tbox::Timer> t_coarsen_fluxsum_comm(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::coarsen_fluxsum_comm"));
   std::shared_ptr<tbox::Timer> t_coarsen_sync_create(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::coarsen_sync_create"));
   std::shared_ptr<tbox::Timer> t_coarsen_sync_comm(
      tbox::TimerManager::getManager()->
      getTimer("algs::MblkHyperbolicLevelIntegrator::coarsen_sync_comm"));

   /*
    * Coarsen flux integrals around fine patch boundaries to coarser level
    * and replace coarse flux information where appropriate.  NULL patch
    * model is passed in to avoid over complicating coarsen process;
    * i.e. patch model is not needed in coarsening of flux integrals.
    */

   std::shared_ptr<xfer::CoarsenSchedule> sched;
   if (d_do_coarsening) {
      t_coarsen_fluxsum_create->start();
      sched = d_mblk_coarsen_fluxsum->createSchedule(mblk_coarse_level,
            mblk_fine_level,
            0);
      t_coarsen_fluxsum_create->stop();

      d_patch_strategy->setDataContext(d_current);
      t_coarsen_fluxsum_comm->start();
      sched->coarsenData();
      t_coarsen_fluxsum_comm->stop();
      d_patch_strategy->clearDataContext();
   }

   /*
    * Repeat conservative difference on coarser level.
    */
   mblk_coarse_level->allocatePatchData(d_saved_var_scratch_data,
      coarse_sim_time);
   mblk_coarse_level->setTime(coarse_sim_time, d_flux_var_data);

   d_patch_strategy->setDataContext(d_scratch);
   t_advance_bdry_fill_comm->start();
   d_mblk_bdry_sched_advance[mblk_coarse_level->getLevelNumber()]->
   fillData(coarse_sim_time);
   t_advance_bdry_fill_comm->stop();

   const double reflux_dt = sync_time - coarse_sim_time;

   for (hier::PatchLevel::iterator mi(mblk_coarse_level->begin());
        mi != mblk_coarse_level->end(); ++mi) {

      (*mi)->allocatePatchData(d_temp_var_scratch_data, coarse_sim_time);

      bool at_syncronization = true;
      d_patch_strategy->conservativeDifferenceOnPatch(**mi,
         coarse_sim_time,
         reflux_dt,
         at_syncronization);
      (*mi)->deallocatePatchData(d_temp_var_scratch_data);

   } // loop over patches

   d_patch_strategy->clearDataContext();

   copyTimeDependentData(mblk_coarse_level, d_scratch, d_new);

   mblk_coarse_level->deallocatePatchData(d_saved_var_scratch_data);

   /*
    * Coarsen time-dependent data from fine patch interiors to coarse patches.
    */

   if (d_do_coarsening) {
      t_coarsen_sync_create->start();
      sched = d_mblk_coarsen_sync_data->createSchedule(mblk_coarse_level,
            mblk_fine_level,
            d_patch_strategy);
      t_coarsen_sync_create->stop();

      d_patch_strategy->setDataContext(d_new);

      t_coarsen_sync_comm->start();
      sched->coarsenData();
      t_coarsen_sync_comm->stop();

      d_patch_strategy->clearDataContext();
   }
}

/*
 *************************************************************************
 *
 * Reset time-dependent data on patch level by replacing current data
 * with new.  The boolean argument is used for odd refinement ratios
 * (in particular 3 used in certain applications).
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::resetTimeDependentData(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double new_time,
   const bool can_be_refined)
{
   TBOX_ASSERT(level);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   double cur_time = 0.;

   for (hier::PatchLevel::iterator mi(level->begin());
        mi != level->end(); ++mi) {

      std::list<std::shared_ptr<hier::Variable> >::iterator time_dep_var =
         d_time_dep_variables.begin();
      while (time_dep_var != d_time_dep_variables.end()) {

         int cur_indx =
            variable_db->mapVariableAndContextToIndex(*time_dep_var,
               d_current);
         int new_indx =
            variable_db->mapVariableAndContextToIndex(*time_dep_var,
               d_new);

         cur_time = (*mi)->getPatchData(cur_indx)->getTime();

         if (can_be_refined && d_number_time_data_levels == 3) {

            int old_indx =
               variable_db->mapVariableAndContextToIndex(*time_dep_var,
                  d_old);

            (*mi)->setPatchData(old_indx, (*mi)->getPatchData(cur_indx));

            (*mi)->setPatchData(cur_indx, (*mi)->getPatchData(new_indx));

         } else {

            if (d_number_time_data_levels == 3) {

               int old_indx =
                  variable_db->mapVariableAndContextToIndex(*time_dep_var,
                     d_old);

               (*mi)->setPatchData(old_indx, (*mi)->getPatchData(cur_indx));

            }

            (*mi)->setPatchData(cur_indx, (*mi)->getPatchData(new_indx));

         }

         (*mi)->deallocatePatchData(new_indx);

         ++time_dep_var;

      }

   } // loop over patches

   level->setTime(new_time, d_new_patch_init_data);

   if (d_number_time_data_levels == 3) {
      level->setTime(cur_time, d_old_time_dep_data);
   }

}

/*
 *************************************************************************
 *
 * Discard new data on level.  This is used primarily to reset patch
 * data after error estimation (e.g., Richardson extrapolation.)
 *
 *************************************************************************
 */

void
MblkHyperbolicLevelIntegrator::resetDataToPreadvanceState(
   const std::shared_ptr<hier::PatchLevel>& level)
{
   TBOX_ASSERT(level);

   /*
    * De-allocate new context
    */
   level->deallocatePatchData(d_new_time_dep_data);

}

/*
 *************************************************************************
 *
 * Register given variable with algorithm according to specified
 * algorithm role (i.e., HYP_VAR_TYPE).  Assignment of descriptor
 * indices to variable lists, component selectors, and communication
 * algorithms takes place here.  The different cases are:
 *
 * TIME_DEP:
 *            The number of factories depends on the number of time
 *            levels of data that must be stored on patches to satisfy
 *            regridding reqs.  Currently, there are two possibilities:
 *
 *            (1) If the coarsen ratios between levels are even, the
 *                error coarsening ratio will be two and so only two
 *                time levels of data must be maintained on every level
 *                but the finest as usual.
 *
 *            (2) If the coarsen ratios between levels are three, and
 *                time integration is used during regridding (e.g., Rich-
 *                ardson extrapolation), then three time levels of data
 *                must be maintained on every level but the finest so
 *                that error estimation can be executed properly.
 *
 *            In case (1), three factories are needed:
 *                         SCRATCH, CURRENT, NEW.
 *            In case (2), four factories are needed:
 *                         SCRATCH, OLD, CURRENT, NEW.
 *
 *            SCRATCH index is added to d_saved_var_scratch_data.
 *            CURRENT index is added to d_new_patch_init_data.
 *            NEW index is added to d_new_time_dep_data.
 *
 * INPUT:
 *            Only one time level of data is maintained and once values
 *            are set on patches, they do not change in time.
 *
 *            Two factories are needed: SCRATCH, CURRENT.
 *
 *            SCRATCH index is added to d_saved_var_scratch_data.
 *            CURRENT index is added to d_new_patch_init_data.
 *
 * NO_FILL:
 *            Only one time level of data is stored and no scratch space
 *            is used.  Data may be set and manipulated at will in user
 *            routines.  Data (including ghost values) is never touched
 *            outside of user routines.
 *
 *            Two factories are needed: CURRENT, SCRATCH.
 *
 *            CURRENT index is added to d_new_patch_init_data.
 *            SCRATCH index is needed only for temporary work space to
 *            fill new patch levels.
 *
 * FLUX:
 *            One factory is needed: SCRATCH.
 *
 *            SCRATCH index is added to d_flux_var_data.
 *
 *            Additionally, a variable for flux integral data is created
 *            for each FLUX variable. It has a single factory, SCRATCH,
 *            which is added to d_fluxsum_data.
 *
 * TEMPORARY:
 *            One factory needed: SCRATCH.
 *            SCRATCH index is added to d_temp_var_scratch_data.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::registerVariable(
   const std::shared_ptr<hier::Variable> var,
   const hier::IntVector ghosts,
   const HYP_VAR_TYPE h_v_type,
   const std::shared_ptr<hier::CoarsenOperator> coarsen_op,
   const std::shared_ptr<hier::RefineOperator> refine_op,
   const std::shared_ptr<hier::TimeInterpolateOperator> time_int)
{
   TBOX_ASSERT(var);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   const hier::IntVector zero_ghosts(d_dim, 0);

   d_all_variables.push_back(var);

   switch (h_v_type) {

      case TIME_DEP: {

         //TBOX_ASSERT(refine_op);
         //TBOX_ASSERT(coarsen_op);
         //TBOX_ASSERT(time_int);

         d_time_dep_variables.push_back(var);

         int cur_id = variable_db->registerVariableAndContext(var,
               d_current,
               zero_ghosts);
         int new_id = variable_db->registerVariableAndContext(var,
               d_new,
               zero_ghosts);
         int scr_id = variable_db->registerVariableAndContext(var,
               d_scratch,
               ghosts);

         d_saved_var_scratch_data.setFlag(scr_id);

         d_new_patch_init_data.setFlag(cur_id);

         d_new_time_dep_data.setFlag(new_id);

         /*
          * Register variable and context needed for restart.
          */
         hier::PatchDataRestartManager::getManager()->
         registerPatchDataForRestart(cur_id);

         /*
          * Set boundary fill schedules for time-dependent variable.
          * If time interpolation operator is non-NULL, regular advance
          * bdry fill algorithm will time interpolate between current and
          * new data on coarser levels, and fill from current data on
          * same level.  New advance bdry fill algorithm will time interpolate
          * between current and new data on coarser levels, and fill from new
          * data on same level.  If time interpolation operator is NULL,
          * regular and new bdry fill algorithms will use current and new
          * data, respectively.
          */

         d_mblk_bdry_fill_advance->registerRefine(
            scr_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);
         d_mblk_bdry_fill_advance_new->registerRefine(
            scr_id, new_id, cur_id, new_id, scr_id, refine_op, time_int);
         d_mblk_fill_new_level->registerRefine(
            cur_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);

         /*
          * For data synchronization between levels, the coarsen algorithm
          * will coarsen new data on finer level to new data on coarser.
          * Recall that coarser level data pointers will not be reset until
          * after synchronization so we always coarsen to new
          * (see synchronizeLevelWithCoarser routine).
          */

         d_mblk_coarsen_sync_data->registerCoarsen(new_id, new_id, coarsen_op);

         d_mblk_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);

         /*
          * Coarsen operations used in Richardson extrapolation.  The init
          * initializes data on coarser level, before the coarse level
          * advance.  If two time levels are used, coarsening occurs between
          * the CURRENT context on both levels.  If three levels are used,
          * coarsening occurs between the OLD context on the fine level and
          * the CURRENT context on the coarse level.  The final coarsen
          * algorithm coarsens data after it has been advanced on the fine
          * level to the NEW context on the coarser level.
          */

         if (d_number_time_data_levels == 3) {

            int old_id = variable_db->registerVariableAndContext(var,
                  d_old,
                  zero_ghosts);
            d_old_time_dep_data.setFlag(old_id);

            d_mblk_bdry_fill_advance_old->registerRefine(
               scr_id, cur_id, old_id, new_id, scr_id, refine_op, time_int);

            d_coarsen_rich_extrap_init->
            registerCoarsen(cur_id, old_id, coarsen_op);

         } else {

            d_coarsen_rich_extrap_init->
            registerCoarsen(cur_id, cur_id, coarsen_op);
         }

         d_coarsen_rich_extrap_final->
         registerCoarsen(new_id, new_id, coarsen_op);

         break;
      }

      case INPUT: {

         //         TBOX_ASSERT(refine_op);
         //         TBOX_ASSERT(coarsen_op);

         int cur_id = variable_db->registerVariableAndContext(var,
               d_current,
               zero_ghosts);
         int scr_id = variable_db->registerVariableAndContext(var,
               d_scratch,
               ghosts);

         d_saved_var_scratch_data.setFlag(scr_id);

         d_new_patch_init_data.setFlag(cur_id);

         /*
          * Register variable and context needed for restart.
          */
         hier::PatchDataRestartManager::getManager()->
         registerPatchDataForRestart(cur_id);

         /*
          * Bdry algorithms for input variables will fill from current only.
          */

         d_mblk_bdry_fill_advance->registerRefine(
            scr_id, cur_id, scr_id, refine_op);
         d_mblk_bdry_fill_advance_new->registerRefine(
            scr_id, cur_id, scr_id, refine_op);
         d_mblk_fill_new_level->registerRefine(
            cur_id, cur_id, scr_id, refine_op);

         /*
          * At initialization, it may be necessary to coarsen INPUT data
          * up through the hierarchy so that all levels are consistent.
          */

         d_mblk_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);

         /*
          * Coarsen operation for setting initial data on coarser level
          * in the Richardson extrapolation algorithm.
          */

         d_coarsen_rich_extrap_init->
         registerCoarsen(cur_id, cur_id, coarsen_op);

         break;
      }

      case NO_FILL: {

         TBOX_ASSERT(refine_op);
         TBOX_ASSERT(coarsen_op);
         int cur_id = variable_db->registerVariableAndContext(var,
               d_current,
               ghosts);

         int scr_id = variable_db->registerVariableAndContext(var,
               d_scratch,
               ghosts);

         d_new_patch_init_data.setFlag(cur_id);

         /*
          * Register variable and context needed for restart.
          */
         hier::PatchDataRestartManager::getManager()->
         registerPatchDataForRestart(cur_id);

         d_mblk_fill_new_level->registerRefine(
            cur_id, cur_id, scr_id, refine_op);

         /*
          * Coarsen operation for setting initial data on coarser level
          * in the Richardson extrapolation algorithm.
          */

         d_coarsen_rich_extrap_init->
         registerCoarsen(cur_id, cur_id, coarsen_op);

         break;
      }

      case FLUX: {

//         TBOX_ASSERT(coarsen_op);
         /*
          * Note that we force all flux variables to hold double precision
          * data and be face- or side-centered.  Also, for each flux variable,
          * a corresponding "fluxsum" variable is created to manage
          * synchronization of data betweeen patch levels in the hierarchy.
          */
         const std::shared_ptr<pdat::FaceVariable<double> > face_var(
            std::dynamic_pointer_cast<pdat::FaceVariable<double>,
                                        hier::Variable>(var));
         const std::shared_ptr<pdat::SideVariable<double> > side_var(
            std::dynamic_pointer_cast<pdat::SideVariable<double>,
                                        hier::Variable>(var));

         if (face_var) {
            if (d_flux_side_registered) {
               TBOX_ERROR(
                  d_object_name << ":  "
                                << "Attempt to register FaceVariable when "
                                << "SideVariable already registered."
                                << std::endl);
            }

            d_flux_is_face = true;

         } else if (side_var) {
            if (d_flux_face_registered) {
               TBOX_ERROR(
                  d_object_name << ":  "
                                << "Attempt to register SideVariable when "
                                << "FaceVariable already registered."
                                << std::endl);
            }

            d_flux_is_face = false;

         } else {
            TBOX_ERROR(
               d_object_name << ":  "
                             << "Flux is neither face- or side-centered." << std::endl);
         }

         d_flux_variables.push_back(var);

         int scr_id = variable_db->registerVariableAndContext(var,
               d_scratch,
               ghosts);

         d_flux_var_data.setFlag(scr_id);

         std::string var_name = var->getName();
         std::string fs_suffix = "_fluxsum";
         std::string fsum_name = var_name;
         fsum_name += fs_suffix;

         std::shared_ptr<hier::Variable> fluxsum;

         if (d_flux_is_face) {
            std::shared_ptr<pdat::FaceDataFactory<double> > fdf(
               SAMRAI_SHARED_PTR_CAST<pdat::FaceDataFactory<double>,
                          hier::PatchDataFactory>(var->getPatchDataFactory()));
            TBOX_ASSERT(fdf);
            fluxsum.reset(new pdat::OuterfaceVariable<double>(
                  d_dim,
                  fsum_name,
                  fdf->getDepth()));
            d_flux_face_registered = true;
         } else {
            std::shared_ptr<pdat::SideDataFactory<double> > sdf(
               SAMRAI_SHARED_PTR_CAST<pdat::SideDataFactory<double>,
                          hier::PatchDataFactory>(var->getPatchDataFactory()));
            TBOX_ASSERT(sdf);
            fluxsum.reset(new pdat::OutersideVariable<double>(
                  d_dim,
                  fsum_name,
                  sdf->getDepth()));
            d_flux_side_registered = true;
         }

         d_fluxsum_variables.push_back(fluxsum);

         int fs_id = variable_db->registerVariableAndContext(fluxsum,
               d_scratch,
               zero_ghosts);

         d_fluxsum_data.setFlag(fs_id);

         d_mblk_coarsen_fluxsum->registerCoarsen(scr_id, fs_id, coarsen_op);

         break;
      }

      case TEMPORARY: {

         int scr_id = variable_db->registerVariableAndContext(var,
               d_scratch,
               ghosts);

         d_temp_var_scratch_data.setFlag(scr_id);

         break;
      }

      default: {

         TBOX_ERROR(
            d_object_name << ":  "
                          << "unknown HYP_VAR_TYPE = " << h_v_type
                          << std::endl);

      }

   }
}

/*
 *************************************************************************
 *
 * Process FLUX and FLUX INTEGRAL data before integration on the level.
 *
 * We allocate FLUX storage if appropriate.
 *
 * If the advance is not temporary, we also zero out the FLUX INTEGRALS
 * on the first step of any level finer than level zero.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::preprocessFluxData(
   const std::shared_ptr<hier::PatchLevel> mblk_level,
   const double cur_time,
   const double new_time,
   const bool regrid_advance,
   const bool first_step,
   const bool last_step)
{
   NULL_USE(cur_time);
   NULL_USE(last_step);

   TBOX_ASSERT(mblk_level);

   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   const int level_number = mblk_level->getLevelNumber();

   if (!regrid_advance) {
      if (((level_number > 0) && first_step) ||
          ((level_number == 0) && !d_have_flux_on_level_zero)) {
         mblk_level->allocatePatchData(d_flux_var_data, new_time);
         if (level_number == 0) {
            d_have_flux_on_level_zero = true;
         }
      }
   } else {
      if (first_step) {
         mblk_level->allocatePatchData(d_flux_var_data, new_time);
      }
   }

   if (!regrid_advance && (level_number > 0)) {

      if (first_step) {

         mblk_level->allocatePatchData(d_fluxsum_data, new_time);

         for (hier::PatchLevel::iterator mi(mblk_level->begin());
              mi != mblk_level->end(); ++mi) {

            std::list<std::shared_ptr<hier::Variable> >::iterator fs_var =
               d_fluxsum_variables.begin();

            while (fs_var != d_fluxsum_variables.end()) {
               int fsum_id =
                  variable_db->mapVariableAndContextToIndex(*fs_var,
                     d_scratch);

               if (d_flux_is_face) {
                  std::shared_ptr<pdat::OuterfaceData<double> > fsum_data(
                     SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<double>, hier::PatchData>(
                        (*mi)->getPatchData(fsum_id)));

                  TBOX_ASSERT(fsum_data);
                  fsum_data->fillAll(0.0);
               } else {
                  std::shared_ptr<pdat::OutersideData<double> > fsum_data(
                     SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
                        (*mi)->getPatchData(fsum_id)));

                  TBOX_ASSERT(fsum_data);

                  fsum_data->fillAll(0.0);
               }

               ++fs_var;
            }
         } // loop over patches

      } else {
         mblk_level->setTime(new_time, d_fluxsum_data);
      }

   } // if ( !regrid_advance && (level_number > 0) )

}

/*
 *************************************************************************
 *
 * Process FLUX and FLUX INTEGRAL data after advancing the solution on
 * the level.  During normal integration steps, the flux integrals are
 * updated for subsequent synchronization by adding FLUX values to
 * flux integrals.
 *
 * If the advance is not temporary (regular integration step):
 * 1) If the level is the finest in the hierarchy, FLUX data is
 *    deallocated.  It is not used during synchronization, and is only
 *    maintained if needed for the advance.
 *
 * 2) If the level is not the coarsest in the hierarchy, update the
 *    flux integrals for later synchronization by adding FLUX values to
 *    flux integrals.
 *
 * If the advance is temporary, deallocate the flux data if first step.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::postprocessFluxData(
   const std::shared_ptr<hier::PatchLevel> mblk_level,
   const bool regrid_advance,
   const bool first_step,
   const bool last_step)
{
   NULL_USE(last_step);

   TBOX_ASSERT(mblk_level);

   if (regrid_advance && first_step) {
      mblk_level->deallocatePatchData(d_flux_var_data);
   }

   if (!regrid_advance && (mblk_level->getLevelNumber() > 0)) {

      for (hier::PatchLevel::iterator mi(mblk_level->begin());
           mi != mblk_level->end(); ++mi) {

         std::list<std::shared_ptr<hier::Variable> >::iterator flux_var =
            d_flux_variables.begin();
         std::list<std::shared_ptr<hier::Variable> >::iterator fluxsum_var =
            d_fluxsum_variables.begin();

         const hier::Index& ilo = (*mi)->getBox().lower();
         const hier::Index& ihi = (*mi)->getBox().upper();

         while (flux_var != d_flux_variables.end()) {

            std::shared_ptr<hier::PatchData> flux_data(
               (*mi)->getPatchData(*flux_var, d_scratch));
            std::shared_ptr<hier::PatchData> fsum_data(
               (*mi)->getPatchData(*fluxsum_var, d_scratch));

            std::shared_ptr<pdat::FaceData<double> > fflux_data;
            std::shared_ptr<pdat::OuterfaceData<double> > ffsum_data;

            std::shared_ptr<pdat::SideData<double> > sflux_data;
            std::shared_ptr<pdat::OutersideData<double> > sfsum_data;

            int ddepth;
            hier::IntVector flux_ghosts(d_dim);

            if (d_flux_is_face) {
               fflux_data = SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>,
                                       hier::PatchData>(flux_data);
               ffsum_data = SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<double>,
                                       hier::PatchData>(fsum_data);

               TBOX_ASSERT(fflux_data && ffsum_data);
               TBOX_ASSERT(fflux_data->getDepth() == ffsum_data->getDepth());

               ddepth = fflux_data->getDepth();
               flux_ghosts = fflux_data->getGhostCellWidth();
            } else {
               sflux_data = SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>,
                                       hier::PatchData>(flux_data);
               sfsum_data = SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>,
                                       hier::PatchData>(fsum_data);

               TBOX_ASSERT(sflux_data && sfsum_data);
               TBOX_ASSERT(sflux_data->getDepth() == sfsum_data->getDepth());

               ddepth = sflux_data->getDepth();
               flux_ghosts = sflux_data->getGhostCellWidth();
            }

            for (int d = 0; d < ddepth; ++d) {
               // loop over lower and upper parts of outer face/side arrays
               for (int ifs = 0; ifs < 2; ++ifs) {
                  if (d_flux_is_face) {
                     if (d_dim == tbox::Dimension(2)) {
                        SAMRAI_F77_FUNC(upfluxsumface2d0, UPFLUXSUMFACE2D0) (
                           ilo(0), ilo(1), ihi(0), ihi(1),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           ifs,
                           fflux_data->getPointer(0, d),
                           ffsum_data->getPointer(0, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumface2d1, UPFLUXSUMFACE2D1) (
                           ilo(0), ilo(1), ihi(0), ihi(1),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           ifs,
                           fflux_data->getPointer(1, d),
                           ffsum_data->getPointer(1, ifs, d));
                     }
                     if (d_dim == tbox::Dimension(3)) {
                        SAMRAI_F77_FUNC(upfluxsumface3d0, UPFLUXSUMFACE3D0) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           fflux_data->getPointer(0, d),
                           ffsum_data->getPointer(0, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumface3d1, UPFLUXSUMFACE3D1) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           fflux_data->getPointer(1, d),
                           ffsum_data->getPointer(1, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumface3d2, UPFLUXSUMFACE3D2) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           fflux_data->getPointer(2, d),
                           ffsum_data->getPointer(2, ifs, d));
                     }
                  } else {
                     if (d_dim == tbox::Dimension(2)) {
                        SAMRAI_F77_FUNC(upfluxsumside2d0, UPFLUXSUMSIDE2D0) (
                           ilo(0), ilo(1), ihi(0), ihi(1),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           ifs,
                           sflux_data->getPointer(0, d),
                           sfsum_data->getPointer(0, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumside2d1, UPFLUXSUMSIDE2D1) (
                           ilo(0), ilo(1), ihi(0), ihi(1),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           ifs,
                           sflux_data->getPointer(1, d),
                           sfsum_data->getPointer(1, ifs, d));
                     }
                     if (d_dim == tbox::Dimension(3)) {
                        SAMRAI_F77_FUNC(upfluxsumside3d0, UPFLUXSUMSIDE3D0) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           sflux_data->getPointer(0, d),
                           sfsum_data->getPointer(0, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumside3d1, UPFLUXSUMSIDE3D1) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           sflux_data->getPointer(1, d),
                           sfsum_data->getPointer(1, ifs, d));
                        SAMRAI_F77_FUNC(upfluxsumside3d2, UPFLUXSUMSIDE3D2) (
                           ilo(0), ilo(1), ilo(2),
                           ihi(0), ihi(1), ihi(2),
                           flux_ghosts(0),
                           flux_ghosts(1),
                           flux_ghosts(2),
                           ifs,
                           sflux_data->getPointer(2, d),
                           sfsum_data->getPointer(2, ifs, d));
                     }
                  }  // if face operations vs. side operations
               }  // loop over lower and upper sides/faces
            }  // loop over depth

            ++flux_var;
            ++fluxsum_var;

         }  // loop over flux variables

      }  // loop over patches

   }  // if !regrid_advance and level number > 0 ....

}

/*
 *************************************************************************
 *
 * Copy time-dependent data from source to destination on level.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::copyTimeDependentData(
   const std::shared_ptr<hier::PatchLevel> level,
   const std::shared_ptr<hier::VariableContext> src_context,
   const std::shared_ptr<hier::VariableContext> dst_context)
{
   TBOX_ASSERT(level);
   TBOX_ASSERT(src_context);
   TBOX_ASSERT(src_context);

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      std::shared_ptr<hier::Patch> patch = *ip;

      std::list<std::shared_ptr<hier::Variable> >::iterator time_dep_var =
         d_time_dep_variables.begin();
      while (time_dep_var != d_time_dep_variables.end()) {
         std::shared_ptr<hier::PatchData> src_data =
            patch->getPatchData(*time_dep_var, src_context);
         std::shared_ptr<hier::PatchData> dst_data =
            patch->getPatchData(*time_dep_var, dst_context);

         TBOX_ASSERT(src_data);
         TBOX_ASSERT(dst_data);

         dst_data->copy(*src_data);
         ++time_dep_var;
      }

   }

}

/*
 *************************************************************************
 *
 * Print all class data for MblkHyperbolicLevelIntegrator object.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::printClassData(
   std::ostream& os) const
{
   os << "\nMblkHyperbolicLevelIntegrator::printClassData..." << std::endl;
   os << "MblkHyperbolicLevelIntegrator: this = "
      << (MblkHyperbolicLevelIntegrator *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_cfl = " << d_cfl << "\n"
      << "d_cfl_init = " << d_cfl_init << std::endl;
   os << "d_lag_dt_computation = " << d_lag_dt_computation << "\n"
      << "d_use_ghosts_for_dt = "
      << d_use_ghosts_for_dt << std::endl;
   os << "d_patch_strategy = "
      << (MblkHyperbolicPatchStrategy *)d_patch_strategy << std::endl;
   os
   << "NOTE: Not printing variable arrays, ComponentSelectors, communication schedules, etc."
   << std::endl;
}

/*
 *************************************************************************
 *
 * Writes out the class version number, d_cfl, d_cfl_init,
 * d_lag_dt_computation, and d_use_ghosts_for_dt to the restart database.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION",
      ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION);

   restart_db->putDouble("d_cfl", d_cfl);
   restart_db->putDouble("d_cfl_init", d_cfl_init);
   restart_db->putBool("d_lag_dt_computation", d_lag_dt_computation);
   restart_db->putBool("d_use_ghosts_for_dt", d_use_ghosts_for_dt);
   restart_db->putBool("d_do_coarsening", d_do_coarsening);
}

/*
 *************************************************************************
 *
 * Reads in cfl, cfl_init, lag_dt_computation, and
 * use_ghosts_to_compute_dt from the input database.
 * Note all restart values are overriden with values from the input
 * database.
 *
 *************************************************************************
 */

void MblkHyperbolicLevelIntegrator::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   if (input_db->keyExists("cfl")) {
      d_cfl = input_db->getDouble("cfl");
   } else {
      if (!is_from_restart) {
         d_cfl = input_db->getDoubleWithDefault("cfl", d_cfl);
      }
   }

   if (input_db->keyExists("cfl_init")) {
      d_cfl_init = input_db->getDouble("cfl_init");
   } else {
      if (!is_from_restart) {
         d_cfl_init = input_db->getDoubleWithDefault("cfl_init", d_cfl_init);
      }
   }

   if (input_db->keyExists("lag_dt_computation")) {
      d_lag_dt_computation = input_db->getBool("lag_dt_computation");
   } else {
      if (!is_from_restart) {
         d_lag_dt_computation =
            input_db->getBoolWithDefault("lag_dt_computation",
               d_lag_dt_computation);
      }
   }

   if (input_db->keyExists("use_ghosts_to_compute_dt")) {
      d_use_ghosts_for_dt = input_db->getBool("use_ghosts_to_compute_dt");
   } else {
      if (!is_from_restart) {
         d_use_ghosts_for_dt =
            input_db->getBoolWithDefault("use_ghosts_for_dt",
               d_use_ghosts_for_dt);
         TBOX_WARNING(
            d_object_name << ":  "
                          << "Key data `use_ghosts_to_compute_dt' not found in input."
                          << "  Using default value "
                          << d_use_ghosts_for_dt << std::endl);
      }
   }

   if (input_db->keyExists("distinguish_mpi_reduction_costs")) {
      d_distinguish_mpi_reduction_costs =
         input_db->getBool("distinguish_mpi_reduction_costs");
   }

   d_do_coarsening = input_db->getBoolWithDefault("do_coarsening", true);
}

/*
 *************************************************************************
 *
 * First, gets the database corresponding to the object_name from the
 * restart file.   If this database exists, this method checks to make
 * sure that the version number of the class matches the version number
 * of the restart file.  If they match, then d_cfl, d_cfl_init,
 * d_lag_dt_computation, and d_use_ghosts_to_compute_dt are read from
 * restart database.
 * Note all restart values can be overriden with values from the input
 * database.
 *
 *************************************************************************
 */
void MblkHyperbolicLevelIntegrator::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file" << std::endl);
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION");
   if (ver != ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
                               << "Restart file version different "
                               << "than class version." << std::endl);
   }

   d_cfl = db->getDouble("d_cfl");
   d_cfl_init = db->getDouble("d_cfl_init");
   d_lag_dt_computation = db->getBool("d_lag_dt_computation");
   d_use_ghosts_for_dt = db->getBool("d_use_ghosts_for_dt");
   d_do_coarsening = db->getBool("d_do_coarsening");
}

/*
 *************************************************************************
 *
 * Utility routines to retrieve variable contexts used by integrator.
 *
 *************************************************************************
 */

std::shared_ptr<hier::VariableContext>
MblkHyperbolicLevelIntegrator::getCurrentContext() const
{
   return d_current;
}

std::shared_ptr<hier::VariableContext>
MblkHyperbolicLevelIntegrator::getNewContext() const
{
   return d_new;
}

std::shared_ptr<hier::VariableContext>
MblkHyperbolicLevelIntegrator::getOldContext() const
{
   return d_old;
}

std::shared_ptr<hier::VariableContext>
MblkHyperbolicLevelIntegrator::getScratchContext() const
{
   return d_scratch;
}

std::shared_ptr<hier::VariableContext>
MblkHyperbolicLevelIntegrator::getPlotContext() const
{
   return d_plot_context;
}

bool
MblkHyperbolicLevelIntegrator::usingRefinedTimestepping() const
{
   return d_use_time_refinement;
}
