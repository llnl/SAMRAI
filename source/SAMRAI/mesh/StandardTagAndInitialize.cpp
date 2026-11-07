/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for performing cell-tagging and initializing
 *                a new level.
 *
 ************************************************************************/
#include "SAMRAI/mesh/StandardTagAndInitialize.h"

#include "SAMRAI/pdat/CellIntegerConstantRefine.h"
#include "SAMRAI/pdat/CellConstantRefine.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Collectives.h"

#include <stdio.h>
#include <memory>

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in Richardson
 * extrapolation algorithm to coarsen tagged cells from fine to coarse
 * level.
 *
 *************************************************************************
 */

extern "C" {
#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif
// in coarsentags1d.f:
void SAMRAI_F77_FUNC(coarsentags1d, COARSENTAGS1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);
// in coarsentags2d.f:
void SAMRAI_F77_FUNC(coarsentags2d, COARSENTAGS2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in coarsentags3d.f:
void SAMRAI_F77_FUNC(coarsentags3d, COARSENTAGS3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const int *, int *);
}

#define DEBUG_TIMES
//#undef DEBUG_TIMES

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * Static function that computes greatest common divisor.
 *
 *************************************************************************
 */

static int
GCD(
   const int a,
   const int b);

/*
 *************************************************************************
 *
 * Constructors and destructor for StandardTagAndInitialize.
 *
 *************************************************************************
 */

StandardTagAndInitialize::StandardTagAndInitialize(
   const std::string& object_name,
   StandardTagAndInitStrategy* tag_strategy,
   const std::shared_ptr<tbox::Database>& input_db):
   TagAndInitializeStrategy(object_name),
   d_tag_strategy(tag_strategy),
   d_error_coarsen_ratio(1),
   d_use_cycle_criteria(false),
   d_use_time_criteria(false),
   d_ever_uses_richardson_extrapolation(false),
   d_ever_uses_gradient_detector(false),
   d_ever_uses_refine_boxes(false),
   d_boxes_changed(false),
   d_old_cycle(-1)
{
   TBOX_ASSERT(!object_name.empty());

   getFromInput(input_db);

   d_cur_time_criteria = d_time_criteria.begin();
   d_cur_cycle_criteria = d_cycle_criteria.begin();

   /*
    * If the user wishes to only use the REFINE_BOXES tagging option,
    * the registered strategy class may be null.  In order to use
    * the GRADIENT_DETECTOR or RICHARDSON_EXTRAPOLATION options, the
    * registered StandardTagAndInitStrategy must be non-NULL.
    */
   if ((d_ever_uses_richardson_extrapolation || d_ever_uses_gradient_detector) &&
       tag_strategy == 0) {
      TBOX_ERROR(
         getObjectName() << ":constructor "
                         << "\nThe supplied implementation of the "
                         << "\nStandardTagAndInitStrategy is NULL.  It must be"
                         << "\nnon-NULL to use the GRADIENT_DETECTOR or"
                         << "\nRICHARDSON_EXTRAPOLATION tagging options."
                         << std::endl);
   }
}

StandardTagAndInitialize::~StandardTagAndInitialize()
{
}

/*
 *************************************************************************
 *
 * Pass requests to initialize level data, reset hierarchy information,
 * and apply an application-specific gradient detector to
 * the subclass of the StandardTagAndInitStrategX data member.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::initializeLevelData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(!old_level || level_number == old_level->getLevelNumber());
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

   if (d_tag_strategy != 0) {
      d_tag_strategy->initializeLevelData(hierarchy,
         level_number,
         init_data_time,
         can_be_refined,
         initial_time,
         old_level,
         allocate_data);
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
}

/*
 *************************************************************************
 *
 * Reset hierarchy configuration information where the range of new
 * hierarchy levels is specified.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::resetHierarchyConfiguration(
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

   if (d_tag_strategy != 0) {
      d_tag_strategy->resetHierarchyConfiguration(hierarchy,
         coarsest_level,
         finest_level);
   }

}

/*
 *************************************************************************
 *
 * Tag cells on level where refinement should occur.   The method can
 * tag cells using either of three options:
 *
 *    1) Richardson extrapolation
 *    2) gradient detection
 *    3) user supplied refine boxes.
 *
 * These options may be used individually or in combination.  If used in
 * combination,  it is IMPORTANT TO PRESERVE THE ORDER of the calls
 * (Richardson extrapolation 1st, gradient detection 2nd, user-supplied
 * refine boxes 3rd) in this method because users may have logic in
 * their code to compare how cells are tagged and changing the order
 * could destroy this logic.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::tagCellsForRefinement(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int regrid_cycle,
   const double regrid_time,
   const int tag_index,
   const bool initial_time,
   const bool coarsest_sync_level,
   const bool can_be_refined,
   const double regrid_start_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(tag_index >= 0);

   bool usesRichExtrap =
      usesRichardsonExtrapolation(regrid_cycle, regrid_time);

   if (usesRichExtrap) {
      tagCellsUsingRichardsonExtrapolation(hierarchy,
         level_number,
         regrid_cycle,
         regrid_time,
         regrid_start_time,
         tag_index,
         initial_time,
         coarsest_sync_level,
         can_be_refined);
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   if (usesGradientDetector(regrid_cycle, regrid_time)) {

      TBOX_ASSERT(d_tag_strategy != 0);

      d_tag_strategy->applyGradientDetector(hierarchy,
         level_number,
         regrid_time,
         tag_index,
         initial_time,
         usesRichExtrap);
   }

   /*
    * If user-supplied refine boxes are to be used, get refine box information
    * from the TagAndInitializeStrategy class, from which this class is
    * derived.
    */
   if (usesRefineBoxes(regrid_cycle, regrid_time)) {

      hier::BoxContainer refine_boxes;
      getUserSuppliedRefineBoxes(
         refine_boxes,
         level_number,
         regrid_cycle,
         regrid_time);

      std::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(level_number));

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<pdat::CellData<int> > tag_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
               patch->getPatchData(tag_index)));

         TBOX_ASSERT(tag_data);

         const hier::Box& tag_data_box = tag_data->getBox();
         const hier::BlockId& tag_data_block_id = tag_data_box.getBlockId();
         for (hier::BoxContainerSingleBlockIterator ib = refine_boxes.begin(tag_data_block_id);
              ib != refine_boxes.end(tag_data_block_id); ++ib) {
            hier::Box intersection = *ib * tag_data_box;
            if (!(intersection.empty())) {
               tag_data->fill(1, intersection);
            }
         }
      }
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }

}

/*
 *************************************************************************
 *
 * The Richardson extrapolation error estimation tags cells according
 * to differences in the solution computed on two different levels of
 * refinement.  The preprocessRichardsonExtrapolation method advanced
 * data on a COARSENED VERSION of the hierarchy level where regridding
 * is applied.  This method advances data on the level itself and
 * compares the solutions on the level and the coarsened version of the
 * level at the advanced time.
 *
 * The steps are summarized as follows:
 *
 * 0) Advance data on the level
 *    0a) if (initial_time) {
 *           advance level for ErrorCoarsenRatio steps with time
 *           increment dt
 *        } else {
 *           advance level by 1 step (see discussion
 *           under 2b for reasons for the difference).
 *        }
 *        NOTE: The "first_step" argument in the
 *        tag_strategy->advanceLevel() method sets the conditions of
 *        the advance.  The following conditions define the state of
 *        "first_step":
 *           - first_step is always true at the initial time.
 *           - at subsequent times, it is true only if:
 *                    level <=  finest level that has not been regridded
 *              -AND- the level is the coarsest used in flux syncs
 *    0b) reset the time dependent data for all but the last step.  The
 *        reason we don't do it on the last step is that we are going to
 *        use the allocated space for step 4.
 *
 * 1) Coarsen data computed on the level to the coarsened version of the
 *    level.  Apply the
 *    tag_strategy->coarsenDataForRichardsonExtrapolation() method again.
 *    This time, before_advance set to false since data has already been
 *    advanced on both levels.
 *
 * 2) Allocate tags and apply the
 *    tag_strategy->applyRichardsonExtrapolation() method.  This method
 *    sets tags on the COARSENED VERSION of the level according to the
 *    criteria set in the tag_strategy.
 *
 * 3) Refine tags from the coarsened version of the level to the level
 *    where tagging is performed.
 *
 * 4) Put data back into the state it was before calling the RE routine.
 *    4a) Apply the d_tag_strategy->resetDataToPreAdvanceState() to take
 *        care of resetting data after the last advance step (see 3b
 *        discussion).
 *    4b) Reset the timestamp of the data if we are at the initial time
 *        by calling the initializeLevelData() method with
 *        "allocate_data"
 *        set to false. We have already allocated and operated on the
 *        initialized data, we just want to reset the timestamp of it.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::tagCellsUsingRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int regrid_cycle,
   const double regrid_time,
   const double regrid_start_time,
   const int tag_index,
   const bool initial_time,
   const bool coarsest_sync_level,
   const bool can_be_refined)
{
   TBOX_ASSERT(regrid_start_time <= regrid_time);
   TBOX_ASSERT(d_tag_strategy != 0);

   const tbox::Dimension& dim = hierarchy->getDim();

   std::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Determine the level timestep.  If the error coarsen ratio is 2
    * (i.e. even refinement ratio) only a single step has been taken
    * between the regrid_start_time and regrid_time, so the timestep is
    * the difference between the regrid_time and the regrid_start_time.
    * If the error coarsen ratio is 3 (i.e. factor 3 refine ratio)
    * then two timesteps have been taken between the regrid_start_time
    * and the regrid_time, so the timestep will be one-half the difference
    * between the regrid_start_time and the regrid_time.
    */
   double dt = (regrid_time - regrid_start_time)
      / static_cast<double>(d_error_coarsen_ratio - 1);

   /*
    * Determine number of advance steps for time integration on the level.
    */
   int n_steps;
   if (initial_time) {
      n_steps = d_error_coarsen_ratio;
   } else {
      n_steps = 1;
      d_tag_strategy->resetTimeDependentData(patch_level,
         regrid_time,
         can_be_refined);
   }

   /*
    * Arguments to the advanceLevel() method
    * are set as follows:
    *   first_step - if the level is NOT the coarsest sync level (that is
    *                the coarsest level to synchronize at the regridding
    *                time), first step is true.  Otherwise, it is false.
    *   last_step - true or false: depends on step count
    *   regrid_advance - true: this is a time-dependent regrid advance
    */

   bool first_step = !coarsest_sync_level;
   bool regrid_advance = true;

   double start_time = regrid_time;
   double end_time = 0.0;
   for (int step_cnt = 0; step_cnt < n_steps; ++step_cnt) {

      end_time = start_time + dt;
      bool last_step = (step_cnt == (n_steps - 1));

#ifdef DEBUG_TIMES
      tbox::plog << "\nAdvancing Data on level in Rich. Extrap" << std::endl;
      tbox::plog << "level number = " << patch_level->getLevelNumber()
                 << std::endl;
      tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
                 << std::endl;
      tbox::plog << "start time = " << start_time << std::endl;
      tbox::plog << "end time = " << end_time << std::endl;
      tbox::plog << "first step? = " << first_step << std::endl;
      tbox::plog << "last step? = " << last_step << std::endl;
      tbox::plog << "regrid advance? = " << regrid_advance << std::endl;
#endif

      (void)d_tag_strategy->advanceLevel(patch_level,
         hierarchy,
         start_time,
         end_time,
         first_step,
         last_step,
         regrid_advance);

      if (step_cnt < (n_steps - 1)) {
         d_tag_strategy->resetTimeDependentData(patch_level,
            end_time,
            can_be_refined);
      }

      start_time = end_time;
   }

   std::shared_ptr<hier::PatchLevel> coarser_level(
      d_rich_extrap_coarsened_levels[level_number]);

   /*
    * Coarsen data from hierarchy level to coarser level.
    */
   bool before_advance = false;
   d_tag_strategy->coarsenDataForRichardsonExtrapolation(hierarchy,
      level_number,
      coarser_level,
      end_time,
      before_advance);
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   coarser_level->allocatePatchData(tag_index, end_time);
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   /*
    * Coarsen tags from level to coarser level.
    */
   hier::IntVector coarsen_ratio(dim, d_error_coarsen_ratio);
   for (hier::PatchLevel::iterator ip(coarser_level->begin());
        ip != coarser_level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& coarse_patch = *ip;
      std::shared_ptr<hier::Patch> fine_patch(
         patch_level->getPatch(coarse_patch->getGlobalId()));
      std::shared_ptr<pdat::CellData<int> > ftags(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            fine_patch->getPatchData(tag_index)));
      std::shared_ptr<pdat::CellData<int> > ctags(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            coarse_patch->getPatchData(tag_index)));

      TBOX_ASSERT(ftags);
      TBOX_ASSERT(ctags);
      TBOX_ASSERT(ctags->getDepth() == ftags->getDepth());

      const hier::Index& filo = ftags->getGhostBox().lower();
      const hier::Index& fihi = ftags->getGhostBox().upper();
      const hier::Index& cilo = ctags->getGhostBox().lower();
      const hier::Index& cihi = ctags->getGhostBox().upper();

      const hier::Index& ifirstc = coarse_patch->getBox().lower();
      const hier::Index& ilastc = coarse_patch->getBox().upper();

      for (int d = 0; d < ctags->getDepth(); ++d) {
         if (dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(coarsentags1d, COARSENTAGS1D) (ifirstc(0), ilastc(0),
               filo(0), fihi(0),
               cilo(0), cihi(0),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));
         } else if ((dim == tbox::Dimension(2))) {
            SAMRAI_F77_FUNC(coarsentags2d, COARSENTAGS2D) (ifirstc(0), ifirstc(1),
               ilastc(0), ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));
         } else if ((dim == tbox::Dimension(3))) {
            SAMRAI_F77_FUNC(coarsentags3d, COARSENTAGS3D) (ifirstc(0), ifirstc(1),
               ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));

         } else {
            TBOX_ERROR("StandardTagAndInitialize error...\n"
               << "DIM > 3 not supported." << std::endl);

         }

      }
   }

   /*
    * Tag cells on coarser level.
    */
   d_tag_strategy->applyRichardsonExtrapolation(coarser_level,
      end_time,
      tag_index,
      dt,
      d_error_coarsen_ratio,
      initial_time,
      usesGradientDetector(regrid_cycle, regrid_time));

   /*
    * Refine tags from coarser level to level.
    */
   //pdat::CellIntegerConstantRefine copytags;
   pdat::CellConstantRefine<int> copytags;
   for (hier::PatchLevel::iterator ip(coarser_level->begin());
        ip != coarser_level->end(); ++ip) {

      const std::shared_ptr<hier::Patch>& coarse_patch = *ip;
      std::shared_ptr<hier::Patch> fine_patch(
         patch_level->getPatch(coarse_patch->getGlobalId()));
      copytags.refine(*fine_patch, *coarse_patch,
         tag_index, tag_index,
         fine_patch->getBox(), coarsen_ratio);
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   /*
    * Final cleanup.  Reset data to initial state before entering this routine.
    */
   d_tag_strategy->resetDataToPreadvanceState(patch_level);

   if (initial_time) {
      bool allocate_data = false;
      initializeLevelData(hierarchy, level_number, regrid_time,
         can_be_refined, initial_time,
         std::shared_ptr<hier::PatchLevel>(),
         allocate_data);
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
}

/*
 *************************************************************************
 *
 * Preprocess data before cell tagging, if appropriate.  For the options
 * provided in this class, only Richardson extrapolation requires
 * any pre-processing.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::preprocessErrorEstimation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int cycle,
   const double regrid_time,
   const double regrid_start_time,
   const bool initial_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

   if (usesRichardsonExtrapolation(cycle, regrid_time)) {
      preprocessRichardsonExtrapolation(hierarchy,
         level_number,
         regrid_time,
         regrid_start_time,
         initial_time);
   }
}

/*
 *************************************************************************
 *
 * The preprocess method for Richardson extrapolation error estimation
 * creates a coarsened version of a level in the hierarchy and advances
 * data on the coarsened level to a prescribed time.
 *
 * The steps are summarized as follows:
 *
 * 0) Create a coarser version of the patch level where tagging is
 *    being performed. The coarser level is coarsened by the
 *    "error coarsen ratio", which is the greatest common divisor of the
 *    refinement ratio (e.g. GCD of ratio 4 refinement would be 2).
 *
 * 1) Initialize data on the coarser level by applying the
 *    tag_strategy->coarsenDataForRichardsonExtrapolation() method.
 *    Note that "before_advance" is set true in this call since we have
 *    not yet advanced data on the coarser level.
 *
 * 2) Advance data on the coarsened version of the level:
 *    2a) get timestep (dt) for level where tagging is performed
 *    2b) if (initial_time) {
 *          Advance coarse level by ErrorCoarsenRatio*dt
 *        } else {
 *          Advance coarse level by dt
 *        }
 *    2c) reset the time dependent data on the coarser level by calling
 *        the tag_strategy's resetTimeDependentData() function.
 *
 * The constructed coarsened levels are stored and used in the
 * tagCellsUsingRichardsonExtrapolation() method.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::preprocessRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double regrid_time,
   const double regrid_start_time,
   const bool initial_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(regrid_start_time <= regrid_time);
   TBOX_ASSERT(d_tag_strategy != 0);

   const tbox::Dimension& dim = hierarchy->getDim();

   std::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Determine the level timestep.  If the error coarsen ratio is 2
    * (i.e. even refinement ratio) only a single step has been taken
    * between the regrid_start_time and regrid_time, so the timestep is
    * the difference between the regrid_time and the regrid_start_time.
    * If the error coarsen ratio is 3 (i.e. factor 3 refine ratio)
    * then two timesteps have been taken between the regrid_start_time
    * and the regrid_time, so the timestep will be one-half the difference
    * between the regrid_start_time and the regrid_time.
    */
   double dt = (regrid_time - regrid_start_time)
      / static_cast<double>(d_error_coarsen_ratio - 1);

   /*
    * Determine start and end times for integration on the coarsened level.
    *
    * At the initial time, the start time is the regrid_time and the end
    * time is the regrid_time + error_coarsen_ratio*dt.
    *
    * In subsequent to the initial time, the start time is the
    * regrid_start_time and the end time is the regrid_time + dt.
    */
   double coarse_start_time;
   double coarse_end_time;

   if (initial_time) {
      coarse_start_time = regrid_time;
      coarse_end_time =
         regrid_time + dt * d_error_coarsen_ratio;
   } else {
      coarse_start_time = regrid_start_time;
      coarse_end_time = regrid_time + dt;
   }

#ifdef DEBUG_TIMES
   tbox::plog << "\nRegridding on level using Rich. Extrap" << std::endl;
   tbox::plog << "level number = " << patch_level->getLevelNumber()
              << std::endl;
   tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
              << std::endl;
   tbox::plog << "coarsened hier level = " << level_number - 1 << std::endl;
   tbox::plog << "coarsened start time = " << coarse_start_time << std::endl;
   tbox::plog << "coarsened end time = " << coarse_end_time << std::endl;
#endif

   /*
    * Generate coarsened version and initialize data on it.  If coarsened
    * level aligns with next coarsened level in hierarchy, set level number
    * so user routines can use this information.
    */

   std::shared_ptr<hier::PatchLevel> coarsened_level(
      std::make_shared<hier::PatchLevel>(dim));
   hier::IntVector coarsen_ratio(dim, d_error_coarsen_ratio);
   coarsened_level->setCoarsenedPatchLevel(patch_level, coarsen_ratio);

   if ((level_number > 0)
       && (hierarchy->getPatchLevel(level_number - 1)->getRatioToLevelZero() ==
           coarsened_level->getRatioToLevelZero())) {
      coarsened_level->setLevelNumber(level_number - 1);
   }
   coarsened_level->setNextCoarserHierarchyLevelNumber(level_number - 1);

   /*
    * Generate Connector patch_level<==>coarsened_level and
    * coarsened_level--->coarsened_level.  To support recursive data
    * transfer, these Connectors should have widths equivalent to
    * patch_level<==>patch_level.
    */
   const hier::IntVector level_to_level_width =
      hierarchy->getRequiredConnectorWidth(level_number,
         level_number);

   const hier::Connector& level_to_level =
      patch_level->findConnector(*patch_level,
         level_to_level_width,
         hier::CONNECTOR_IMPLICIT_CREATION_RULE);

   std::shared_ptr<hier::Connector> tmp_coarsened(
      std::make_shared<hier::Connector>(level_to_level));
   tmp_coarsened->setBase(*coarsened_level->getBoxLevel());
   tmp_coarsened->setHead(*coarsened_level->getBoxLevel());
   tmp_coarsened->setWidth(
      hier::IntVector::ceilingDivide(level_to_level_width, coarsen_ratio), true);
   tmp_coarsened->coarsenLocalNeighbors(coarsen_ratio);
   tmp_coarsened->setTranspose(0, false);

   std::shared_ptr<hier::Connector> level_to_coarsened(
      std::make_shared<hier::Connector>(*tmp_coarsened));
   level_to_coarsened->setBase(*patch_level->getBoxLevel());
   level_to_coarsened->setHead(*coarsened_level->getBoxLevel());
   level_to_coarsened->setWidth(level_to_level_width, true);
   level_to_coarsened->setTranspose(0, false);

   coarsened_level->cacheConnector(tmp_coarsened);

   if (level_number > 0) {
      /*
       * Get Connectors coarsened<==>coarser, which are used for recursive
       * refinement filling of the coarsened level's ghosts.
       */
      std::shared_ptr<hier::Connector> coarser_to_coarsened;
      std::shared_ptr<hier::PatchLevel> coarser_level(
         hierarchy->getPatchLevel(level_number - 1));
      const hier::Connector& coarser_to_level =
         coarser_level->findConnectorWithTranspose(*patch_level,
            hierarchy->getRequiredConnectorWidth(
               level_number - 1, level_number, true),
            hierarchy->getRequiredConnectorWidth(
               level_number, level_number - 1, true),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE);

      // level_to_coarsened only needs its transpose set temporarily for this
      // call to bridge.  When it is cached later we don't want to also cache
      // its transpose which is why the tranpose is set to a null shared_ptr a
      // few lines below.
      hier::Connector coarsened_to_level(level_to_level);
      coarsened_to_level.setBase(*coarsened_level->getBoxLevel());
      coarsened_to_level.setHead(*patch_level->getBoxLevel());
      coarsened_to_level.setWidth(
         hier::IntVector::ceilingDivide(level_to_level_width, coarsen_ratio),
         true);
      level_to_coarsened->setTranspose(&coarsened_to_level, false);

      hier::OverlapConnectorAlgorithm oca;
      oca.bridge(coarser_to_coarsened,
         coarser_to_level,
         *level_to_coarsened,
         true);
      coarser_level->cacheConnector(coarser_to_coarsened);

      level_to_coarsened->setTranspose(0, false);
   }
   patch_level->cacheConnector(level_to_coarsened);

   bool before_advance = true;
   d_tag_strategy->coarsenDataForRichardsonExtrapolation(hierarchy,
      level_number,
      coarsened_level,
      coarse_start_time,
      before_advance);

   /*
    * Advance data on coarsened level.  Arguments to the advanceLevel() method
    * are set as follows:
    *   first_step - true: this is the first step on the coarsened level
    *                so it is necessary to do any required
    *                setup in the advance method.
    *   last_step - true: only one step will occur on the the coarsened level
    *   regrid_advance - true: this is a time-dependent regrid advance
    */
   bool first_step = true;
   bool last_step = true;
   bool regrid_advance = true;

#ifdef DEBUG_TIMES
   tbox::plog << "\nAdvancing Data on coarsened in Rich. Extrap" << std::endl;
   tbox::plog << "level number = " << coarsened_level->getLevelNumber()
              << std::endl;
   tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
              << std::endl;
   tbox::plog << "start time = " << coarse_start_time << std::endl;
   tbox::plog << "end time = " << coarse_end_time << std::endl;
   tbox::plog << "first step? = " << first_step << std::endl;
   tbox::plog << "last step? = " << last_step << std::endl;
   tbox::plog << "regrid advance? = " << regrid_advance << std::endl;
#endif

   (void)d_tag_strategy->advanceLevel(coarsened_level,
      hierarchy,
      coarse_start_time,
      coarse_end_time,
      first_step,
      last_step,
      regrid_advance);

   /*
    * Reset data on the coarse level.  Since this level does not reside
    * in the hierarchy, it cannot be refined (and no special storage
    * manipulation is needed; so we set this to false.
    */
   bool level_can_be_refined = false;
   d_tag_strategy->resetTimeDependentData(coarsened_level,
      coarse_end_time,
      level_can_be_refined);

   /*
    * Add the constructed coarsened level to the array of maintained
    * coarsened levels for Richardson extrapolation.
    */
   if (static_cast<int>(d_rich_extrap_coarsened_levels.size()) < level_number + 1) {
      d_rich_extrap_coarsened_levels.resize(level_number + 1);
   }

   d_rich_extrap_coarsened_levels[level_number] = coarsened_level;

}

/*
 *************************************************************************
 *
 * Boxes on the coarsest level must be able to be coarsened by the
 * error coarsen ratio to apply Richardson Extrapolation.  This method
 * simply checks that this is the case.
 *
 *************************************************************************
 */

bool
StandardTagAndInitialize::coarsestLevelBoxesOK(
   const hier::BoxContainer& boxes) const
{
   TBOX_ASSERT(!boxes.empty());

   bool boxes_ok = true;
   if (everUsesRichardsonExtrapolation()) {

      const tbox::Dimension& dim = boxes.begin()->getDim();

      for (hier::BoxContainer::const_iterator ib = boxes.begin();
           ib != boxes.end(); ++ib) {
         hier::IntVector n_cells = ib->numberCells();
         for (int i = 0; i < dim.getValue(); ++i) {
            int error_coarsen_ratio = getErrorCoarsenRatio();
            if (!((n_cells(i) % error_coarsen_ratio) == 0)) {
               tbox::perr << "Bad domain box: " << *ib << std::endl;
               TBOX_WARNING(
                  getObjectName() << "At least one box on the\n"
                                  << "coarsest level could not be coarsened by the ratio: "
                                  << error_coarsen_ratio << std::endl);
               boxes_ok = false;
            }
         }
      }
   }
   return boxes_ok;

}

/*
 *************************************************************************
 *
 * Compute Error coarsen ratio for Richardson extrapolation. For a given
 * level, the error coarsen ratio should be the greatest common divisor
 * (GCD) of the refinement ratio applied to the level.  This value
 * should generally be 2 or 3 (e.g. refinement ratio=2 gives GCD=2;
 * rr=3 gives GCD=3; rr=4 gives GCD=2; etc.).
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::checkCoarsenRatios(
   const std::vector<hier::IntVector>& ratio_to_coarser)
{
   if (everUsesRichardsonExtrapolation()) {
      const tbox::Dimension& dim = ratio_to_coarser[1].getDim();

      /*
       * Compute GCD on first coordinate direction of level 1
       */
      int error_coarsen_ratio = 0;
      int gcd_level1 = ratio_to_coarser[1](0,0);
      if ((gcd_level1 % 2) == 0) {
         error_coarsen_ratio = 2;
      } else if ((gcd_level1 % 3) == 0) {
         error_coarsen_ratio = 3;
      } else {
         TBOX_ERROR("Unable to perform Richardson extrapolation algorithm "
            << "with ratio_to_coarser[1](0) = " << gcd_level1 << std::endl);
      }

      /*
       * Iterate through levels and check the coarsen ratios to make sure the
       * error coarsen ratios computed in every coordinate direction on every
       * level are between the supported 2 or 3, and that the error coarsen
       * ratios are constant over the hierarchy.
       */
      for (int ln = 1; ln < static_cast<int>(ratio_to_coarser.size()); ++ln) {

         for (int d = 0; d < dim.getValue(); ++d) {
            int gcd = GCD(error_coarsen_ratio, ratio_to_coarser[ln](0,d));
            if ((gcd % error_coarsen_ratio) != 0) {
               gcd = ratio_to_coarser[ln](0,d);
               TBOX_ERROR(
                  getObjectName() << "\n"
                                  << "Unable to perform Richardson extrapolation because\n"
                                  << "the error coarsen ratio computed from the\n"
                                  << "ratio_to_coarser entries is not constant across all\n"
                                  << "levels, in all coordinate directions, of the hierarchy. In\n"
                                  << "order to use Richardson extrapolation, the minimum\n"
                                  << "divisor (> 1) of all the ratio_to_coarser entries must\n"
                                  << "be 2 -or- 3:\n"
                                  << "   level 1(0): minimum divisor: "
                                  << error_coarsen_ratio
                                  << "\n   level " << ln << "(" << d
                                  << "):"
                                  << ": ratio_to_coarser = " << gcd
                                  << std::endl);
            }
         }
      }

      d_error_coarsen_ratio = error_coarsen_ratio;

   }

}

/*
 *************************************************************************
 * Returns true if there is ever a tagging crtieria which advances the
 * solution data in the regridding process.
 *************************************************************************
 */
bool
StandardTagAndInitialize::everUsesTimeIntegration() const
{
   return everUsesRichardsonExtrapolation();
}

/*
 *************************************************************************
 * Returns true if the tagging crtieria for the supplied cycle/time
 * advances the solution data in the regridding process.
 *************************************************************************
 */
bool
StandardTagAndInitialize::usesTimeIntegration(
   int cycle,
   double time)
{
   return usesRichardsonExtrapolation(cycle, time);
}

/*
 *************************************************************************
 * Returns true if there is ever a Richardson extrapolation tagging crtieria.
 *************************************************************************
 */
bool
StandardTagAndInitialize::everUsesRichardsonExtrapolation() const
{
   return d_ever_uses_richardson_extrapolation;
}

/*
 *************************************************************************
 * Returns true if there is a Richardson extrapolation tagging crtieria
 * for the supplied cycle/time.
 *************************************************************************
 */
bool
StandardTagAndInitialize::usesRichardsonExtrapolation(
   int cycle,
   double time)
{
   TBOX_ASSERT(!d_use_cycle_criteria || !d_use_time_criteria);

   bool result = false;

   setCurrentTaggingCriteria(cycle, time);
   if (d_use_cycle_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_cycle_criteria->d_tag_criteria.begin();
           i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
            result = true;
            break;
         }
      }
   } else if (d_use_time_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_time_criteria->d_tag_criteria.begin();
           i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
            result = true;
            break;
         }
      }
   }
   return result;
}

/*
 *************************************************************************
 * Returns true if there is ever a gradient detector tagging crtieria.
 *************************************************************************
 */
bool
StandardTagAndInitialize::everUsesGradientDetector() const
{
   return d_ever_uses_gradient_detector;
}

/*
 *************************************************************************
 * Returns true if there is a gradient detector tagging crtieria for the
 * supplied cycle/time.
 *************************************************************************
 */
bool
StandardTagAndInitialize::usesGradientDetector(
   int cycle,
   double time)
{
   TBOX_ASSERT(!d_use_cycle_criteria || !d_use_time_criteria);

   bool result = false;

   setCurrentTaggingCriteria(cycle, time);
   if (d_use_cycle_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_cycle_criteria->d_tag_criteria.begin();
           i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "GRADIENT_DETECTOR") {
            result = true;
            break;
         }
      }
   } else if (d_use_time_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_time_criteria->d_tag_criteria.begin();
           i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "GRADIENT_DETECTOR") {
            result = true;
            break;
         }
      }
   }
   return result;
}

/*
 *************************************************************************
 * Returns true if there is ever a refine boxes tagging crtieria.
 *************************************************************************
 */
bool
StandardTagAndInitialize::everUsesRefineBoxes() const
{
   return d_ever_uses_refine_boxes;
}

/*
 *************************************************************************
 * Returns true if there is a refine boxes tagging crtieria for the
 * supplied cycle/time.
 *************************************************************************
 */
bool
StandardTagAndInitialize::usesRefineBoxes(
   int cycle,
   double time)
{
   TBOX_ASSERT(!d_use_cycle_criteria || !d_use_time_criteria);

   bool result = false;

   setCurrentTaggingCriteria(cycle, time);
   if (d_use_cycle_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_cycle_criteria->d_tag_criteria.begin();
           i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "REFINE_BOXES") {
            result = true;
            break;
         }
      }
   } else if (d_use_time_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_time_criteria->d_tag_criteria.begin();
           i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "REFINE_BOXES") {
            result = true;
            break;
         }
      }
   }
   return result;
}

int
StandardTagAndInitialize::getErrorCoarsenRatio() const
{
   return d_error_coarsen_ratio;
}

bool
StandardTagAndInitialize::refineUserBoxInputOnly(
   int cycle,
   double time)
{
   bool use_only_refine_boxes = false;
   if (usesRefineBoxes(cycle, time)) {
      use_only_refine_boxes = true;
      if (usesGradientDetector(cycle, time) ||
          usesRichardsonExtrapolation(cycle, time)) {
         use_only_refine_boxes = false;
      }
   }
   return use_only_refine_boxes;
}

/*
 *************************************************************************
 *
 * Read cell tagging option and, if required, specified refinement boxes.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{
   /*
    * If no input database is provided, no criteria is set to tag cells
    * so cell-tagging will not occur.  Print a warning to indicate if
    * this is the case.
    */
   if (input_db) {

      /*
       * Allow shortcut syntax.
       */
      if (input_db->keyExists("tagging_method")) {
         // Shortcut syntax read here.

         CycleTagCriteria this_cycle_crit;
         TagCriteria this_tag_crit;
         this_cycle_crit.d_cycle = 0;
         std::string tagging_method = input_db->getString("tagging_method");
         if (!(tagging_method == "RICHARDSON_EXTRAPOLATION" ||
               tagging_method == "GRADIENT_DETECTOR" ||
               tagging_method == "REFINE_BOXES" ||
               tagging_method == "NONE")) {
            INPUT_VALUE_ERROR("tagging_method");
         }
         this_tag_crit.d_tagging_method = tagging_method;
         if (tagging_method == "REFINE_BOXES") {
            d_ever_uses_refine_boxes = true;
            std::vector<std::string> level_keys = input_db->getAllKeys();
            int n_level_keys = static_cast<int>(level_keys.size());
            if (n_level_keys <= 1) {
               TBOX_ERROR(
                  getObjectName() << "::getFromInput \n"
                                  << "No refine boxes supplied."
                                  << std::endl);
            }
            for (int k = 0; k < n_level_keys; ++k) {
               hier::BoxContainer level_boxes;
               if (level_keys[k] == "tagging_method") {
                  continue;
               }
               if (level_keys[k].find("level_") != 0) {
                  TBOX_ERROR(
                     getObjectName() << "::getFromInput \n"
                                     << "Invalid syntax for level refine boxes."
                                     << std::endl);
               }
               int level = atoi(level_keys[k].substr(6).c_str());
               std::shared_ptr<tbox::Database> level_db(
                  input_db->getDatabase(level_keys[k]));
               std::vector<std::string> block_keys = level_db->getAllKeys();
               int n_block_keys = static_cast<int>(block_keys.size());
               if (n_block_keys <= 0) {
                  TBOX_ERROR(
                     getObjectName() << "::getFromInput \n"
                                     << "No refine boxes supplied."
                                     << std::endl);
               }
               if ((n_block_keys == 1) && (block_keys[0] == "boxes")) {
                  std::vector<tbox::DatabaseBox> db_box_vector =
                     level_db->getDatabaseBoxVector("boxes");
                  hier::BoxContainer boxes(db_box_vector);
                  for (hier::BoxContainer::iterator b = boxes.begin();
                       b != boxes.end(); ++b) {
                     b->setBlockId(hier::BlockId(0));
                  }
                  level_boxes.spliceBack(boxes);
               } else {
                  for (int l = 0; l < n_block_keys; ++l) {
                     if (block_keys[l].find("block_") != 0) {
                        TBOX_ERROR(
                           getObjectName() << "::getFromInput \n"
                                           << "Invalid syntax for level refine boxes."
                                           << std::endl);
                     }
                     int block = atoi(block_keys[l].substr(6).c_str());
                     std::shared_ptr<tbox::Database> block_db(
                        level_db->getDatabase(block_keys[l]));
                     std::vector<tbox::DatabaseBox> db_box_vector =
                        block_db->getDatabaseBoxVector("boxes");
                     hier::BoxContainer boxes(db_box_vector);
                     for (hier::BoxContainer::iterator b = boxes.begin();
                          b != boxes.end(); ++b) {
                        b->setBlockId(hier::BlockId(block));
                     }
                     level_boxes.spliceBack(boxes);
                     boxes.clear();
                  }
               }
               this_tag_crit.d_level_refine_boxes.insert(
                  std::make_pair(level, level_boxes));
               level_boxes.clear();
            }
         } else if (tagging_method == "RICHARDSON_EXTRAPOLATION") {
            d_ever_uses_richardson_extrapolation = true;
         } else if (tagging_method == "GRADIENT_DETECTOR") {
            d_ever_uses_gradient_detector = true;
         } else {
            TBOX_WARNING(
               getObjectName() << "::getFromInput\n"
                               << "NO METHOD IS SPECIFIED TO TAG CELLS FOR\n"
                               << "REFINEMENT so no tagging is performed."
                               << std::endl);
         }
         this_cycle_crit.d_tag_criteria.push_back(this_tag_crit);
         d_cycle_criteria.insert(this_cycle_crit);
      } else {
         // Complete syntax read here.

         // Read the set of criteria for each cycle or time having tagging
         // criteria.
         int n_at_commands = static_cast<int>(input_db->getAllKeys().size());
         for (int i = 0; i < n_at_commands; ++i) {
            std::string at_name = "at_" + tbox::Utilities::intToString(i);
            if (!input_db->keyExists(at_name)) {
               TBOX_ERROR(
                  getObjectName() << "::getFromInput \n"
                                  << "Missing tagging criteria " << i << "."
                                  << std::endl);
            }
            std::shared_ptr<tbox::Database> at_db(
               input_db->getDatabase(at_name));

            // Read information specific to a cycle or a time criteria.
            CycleTagCriteria this_cycle_crit;
            TimeTagCriteria this_time_crit;
            bool is_cycle = false;
            if (at_db->keyExists("cycle")) {
               /*
                * Read cycle tagging criteria.
                */
               int cycle = at_db->getInteger("cycle");
               if (!(cycle >= 0)) {
                  INPUT_RANGE_ERROR("cycle");
               }
               this_cycle_crit.d_cycle = cycle;
               is_cycle = true;
            } else if (at_db->keyExists("time")) {
               /*
                * Read time tagging criteria.
                */
               double time = at_db->getDouble("time");
               if (!(time >= 0)) {
                  INPUT_RANGE_ERROR("time");
               }
               this_time_crit.d_time = time;
            } else {
               TBOX_ERROR(
                  getObjectName() << "::getFromInput \n"
                                  << "Invalid tagging period, must be 'cycle' or 'time'."
                                  << std::endl);
            }

            /*
             * Read info common to cycle and time tagging criteria.
             */
            int n_tag_keys = static_cast<int>(at_db->getAllKeys().size()) - 1;
            if (n_tag_keys <= 0) {
               TBOX_ERROR(
                  getObjectName() << "::getFromInput \n"
                                  << "No tagging criteria supplied for tag "
                                  << i << "."
                                  << std::endl);
            }

            // Read the set of criteria at this cycle or time.
            for (int j = 0; j < n_tag_keys; ++j) {
               TagCriteria this_tag_crit;
               std::string tag_name = "tag_" + tbox::Utilities::intToString(j);
               if (!at_db->keyExists(tag_name)) {
                  TBOX_ERROR(
                     getObjectName() << "::getFromInput \n"
                                     << "Missing tag criteria." << std::endl);
               }
               std::shared_ptr<tbox::Database> this_tag_db(
                  at_db->getDatabase(tag_name));
               std::string tagging_method =
                  this_tag_db->getString("tagging_method");
               if (tagging_method != "RICHARDSON_EXTRAPOLATION" &&
                   tagging_method != "GRADIENT_DETECTOR" &&
                   tagging_method != "REFINE_BOXES" &&
                   tagging_method != "NONE") {
                  TBOX_ERROR(
                     getObjectName() << "::getFromInput \n"
                                     << "Invalid tagging_method supplied."
                                     << std::endl);
               }
               this_tag_crit.d_tagging_method = tagging_method;

               // If REFINE_BOXES is the method then the boxes for each level
               // need to be read.
               if (tagging_method == "REFINE_BOXES") {
                  d_ever_uses_refine_boxes = true;
                  std::vector<std::string> level_keys = this_tag_db->getAllKeys();
                  int n_level_keys = static_cast<int>(level_keys.size());
                  if (n_level_keys <= 1) {
                     TBOX_ERROR(
                        getObjectName() << "::getFromInput \n"
                                        << "No refine boxes supplied."
                                        << std::endl);
                  }

                  // For each level specified, read the refine boxes.
                  for (int k = 0; k < n_level_keys; ++k) {
                     hier::BoxContainer level_boxes;
                     if (level_keys[k] == "tagging_method") {
                        continue;
                     }
                     if (level_keys[k].find("level_") != 0) {
                        TBOX_ERROR(
                           getObjectName() << "::getFromInput \n"
                                           << "Invalid syntax for level refine boxes."
                                           << std::endl);
                     }
                     int level = atoi(level_keys[k].substr(6).c_str());
                     std::shared_ptr<tbox::Database> level_db(
                        this_tag_db->getDatabase(level_keys[k]));
                     std::vector<std::string> block_keys =
                        level_db->getAllKeys();
                     int n_block_keys = static_cast<int>(block_keys.size());
                     if (n_block_keys <= 0) {
                        TBOX_ERROR(
                           getObjectName() << "::getFromInput \n"
                                           << "No refine boxes supplied."
                                           << std::endl);
                     }
                     if ((n_block_keys == 1) && (block_keys[0] == "boxes")) {
                        std::vector<tbox::DatabaseBox> db_box_vector =
                           level_db->getDatabaseBoxVector("boxes");
                        hier::BoxContainer boxes(db_box_vector);
                        for (hier::BoxContainer::iterator b = boxes.begin();
                             b != boxes.end(); ++b) {
                           b->setBlockId(hier::BlockId(0));
                        }
                        level_boxes.spliceBack(boxes);
                     } else {
                        for (int l = 0; l < n_block_keys; ++l) {
                           if (block_keys[l].find("block_") != 0) {
                              TBOX_ERROR(
                                 getObjectName() << "::getFromInput \n"
                                                 << "Invalid syntax for level refine boxes."
                                                 << std::endl);
                           }
                           int block = atoi(block_keys[l].substr(6).c_str());
                           std::shared_ptr<tbox::Database> block_db(
                              level_db->getDatabase(block_keys[l]));
                           std::vector<tbox::DatabaseBox> db_box_vector =
                              block_db->getDatabaseBoxVector("boxes");
                           hier::BoxContainer boxes(db_box_vector);
                           for (hier::BoxContainer::iterator b = boxes.begin();
                                b != boxes.end(); ++b) {
                              b->setBlockId(hier::BlockId(block));
                           }
                           level_boxes.spliceBack(boxes);
                           boxes.clear();
                        }
                     }
                     this_tag_crit.d_level_refine_boxes.insert(
                        std::make_pair(level, level_boxes));
                     level_boxes.clear();
                  }
               } else if (tagging_method == "RICHARDSON_EXTRAPOLATION") {
                  d_ever_uses_richardson_extrapolation = true;
               } else if (tagging_method == "GRADIENT_DETECTOR") {
                  d_ever_uses_gradient_detector = true;
               } else if (n_tag_keys != 1) {
                  TBOX_ERROR(
                     getObjectName() << "::getFromInput \n"
                                     << "A tagging method of NONE has been specified \n"
                                     << "along with other tagging methods for the \n."
                                     << "same cycle or time.  This is contradictory."
                                     << std::endl);
               }
               // Add the tagging criteria to the vector of criteria.
               if (is_cycle) {
                  this_cycle_crit.d_tag_criteria.push_back(this_tag_crit);
               } else {
                  this_time_crit.d_tag_criteria.push_back(this_tag_crit);
               }

               // Clear the tagging criteria's boxes for the next criteria in
               // the set.
               this_tag_crit.d_level_refine_boxes.clear();
            }

            // Add the collection of tagging criteria for this cycle or time to
            // the appropriate set.
            if (is_cycle) {
               d_cycle_criteria.insert(this_cycle_crit);
            } else {
               d_time_criteria.insert(this_time_crit);
            }

            // Clear the collection of tagging criteria for the next cycle or
            // time.
            this_cycle_crit.d_tag_criteria.clear();
         }
      }

      if (d_cycle_criteria.empty() && d_time_criteria.empty()) {
         TBOX_WARNING(
            getObjectName() << "::getFromInput\n"
                            << "NO METHOD IS SPECIFIED TO TAG CELLS FOR\n"
                            << "REFINEMENT so no tagging is performed."
                            << std::endl);
      }
   } else {
      TBOX_WARNING(
         getObjectName() << "::getFromInput\n"
                         << "no input database specified - NO METHOD IS SPECIFIED TO TAG\n"
                         << "CELLS FOR REFINEMENT so no tagging is performed."
                         << std::endl);
   }

}

/*
 *************************************************************************
 *
 * Sets refine boxes for case where refine region is specified by the
 * user.  The bool return value specifies whether or not the refine
 * boxes have been reset from the last time the method was called
 * (true = they are reset, false = they have NOT changed).
 *
 * Note that if any method which invokes tagging is performed there
 * is always potential that the boxes have changed so this method will
 * always return true in this case.
 *
 *************************************************************************
 */
bool
StandardTagAndInitialize::getUserSuppliedRefineBoxes(
   hier::BoxContainer& refine_boxes,
   const int level_num,
   const int cycle,
   const double time)
{
   TBOX_ASSERT(refine_boxes.empty());
   TBOX_ASSERT(level_num >= 0);
   TBOX_ASSERT(time >= 0.);
   TBOX_ASSERT(!d_use_cycle_criteria || !d_use_time_criteria);

   setCurrentTaggingCriteria(cycle, time);

   /*
    * Get the boxes to refine on the supplied level depending on whether we're
    * using a scycle or a time criteria.
    */
   if (d_use_cycle_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_cycle_criteria->d_tag_criteria.begin();
           i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "REFINE_BOXES") {
            std::map<int, hier::BoxContainer>::const_iterator boxes =
               i->d_level_refine_boxes.find(level_num);
            if (boxes != i->d_level_refine_boxes.end()) {
               refine_boxes = boxes->second;
            }
            break;
         }
      }
   } else if (d_use_time_criteria) {
      for (std::vector<TagCriteria>::const_iterator i = d_cur_time_criteria->d_tag_criteria.begin();
           i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "REFINE_BOXES") {
            std::map<int, hier::BoxContainer>::const_iterator boxes =
               i->d_level_refine_boxes.find(level_num);
            if (boxes != i->d_level_refine_boxes.end()) {
               refine_boxes = boxes->second;
            }
            break;
         }
      }
   }

   /*
    * If the user has requested their own particular set of refine
    * boxes (i.e. by calling resetRefineBoxes()), overwrite any previously
    * determined refine boxes with their requested set.
    */
   bool use_reset = false;
   if (static_cast<int>(d_refine_boxes_reset.size()) > level_num) {
      if (d_refine_boxes_reset[level_num]) {
         use_reset = true;
         refine_boxes = d_reset_refine_boxes[level_num];
      }
   }

   if (refine_boxes.empty()) {
      TBOX_WARNING(
         getObjectName() << ": getRefineBoxes\n"
                         << "No refine boxes specified on level " << level_num
                         << ".\n No refinement will be performed."
                         << std::endl);
   }

   /*
    * If we have reset the boxes or the boxes have changed due to a change in
    * the tagging criteria from the last cycle/time then return "true".
    */
   bool modified_refine_boxes = false;
   if (use_reset || d_boxes_changed) {
      modified_refine_boxes = true;
   }

   return modified_refine_boxes;

}

/*
 *************************************************************************
 *
 * Resets refine boxes for specified level.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::resetRefineBoxes(
   const hier::BoxContainer& refine_boxes,
   const int level_num)
{
   TBOX_ASSERT(level_num >= 0);

   int i = static_cast<int>(d_reset_refine_boxes.size());
   if (i <= level_num) {
      d_reset_refine_boxes.resize(level_num + 1);
      d_refine_boxes_reset.resize(level_num + 1);
      for ( ; i < static_cast<int>(d_reset_refine_boxes.size()); ++i) {
         d_refine_boxes_reset[i] = false;
      }
   }

   d_refine_boxes_reset[level_num] = true;
   d_reset_refine_boxes[level_num] = refine_boxes;

}

void
StandardTagAndInitialize::turnOnRefineBoxes(
   double time)
{
   TimeTagCriteria search_for;
   search_for.d_time = time;
   std::set<TimeTagCriteria, time_tag_criteria_less>::iterator existing =
      d_time_criteria.find(search_for);
   if (existing == d_time_criteria.end()) {
      TimeTagCriteria this_time_crit;
      TagCriteria this_tag_crit;
      this_tag_crit.d_tagging_method = "REFINE_BOXES";
      this_time_crit.d_time = time;
      this_time_crit.d_tag_criteria.push_back(this_tag_crit);
      d_cur_time_criteria = d_time_criteria.insert(this_time_crit).first;
   } else {
      bool refine_boxes_already_on = false;
      for (std::vector<TagCriteria>::const_iterator i = existing->d_tag_criteria.begin();
           i != existing->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "REFINE_BOXES") {
            refine_boxes_already_on = true;
            break;
         }
      }
      if (!refine_boxes_already_on) {
         TagCriteria this_tag_crit;
         this_tag_crit.d_tagging_method = "REFINE_BOXES";
         std::vector<TagCriteria>& this_tag_criteria =
            const_cast<std::vector<TagCriteria>&>(existing->d_tag_criteria);
         this_tag_criteria.push_back(this_tag_crit);
      }
   }
}

void
StandardTagAndInitialize::turnOffRefineBoxes(
   double time)
{
   for (std::set<TimeTagCriteria, time_tag_criteria_less>::iterator i = d_time_criteria.begin();
        i != d_time_criteria.end(); ++i) {
      if (i->d_time <= time) {
         std::vector<TagCriteria>& tag_crits =
            const_cast<std::vector<TagCriteria>&>(i->d_tag_criteria);
         std::vector<TagCriteria>::iterator j = tag_crits.begin();
         while (j != tag_crits.end()) {
            if (j->d_tagging_method == "REFINE_BOXES") {
               tag_crits.erase(j);
            } else {
               ++j;
            }
         }
         break;
      }
   }
}

void
StandardTagAndInitialize::turnOnGradientDetector(
   double time)
{
   if (!d_tag_strategy) {
      TBOX_ERROR("StandardTagAndInitialize::turnOnGradientDetector\n"
         << "A tagging strategy must be defined if graient detector is used.\n");
   }

   TimeTagCriteria search_for;
   search_for.d_time = time;
   std::set<TimeTagCriteria, time_tag_criteria_less>::iterator existing =
      d_time_criteria.find(search_for);
   if (existing == d_time_criteria.end()) {
      TimeTagCriteria this_time_crit;
      TagCriteria this_tag_crit;
      this_tag_crit.d_tagging_method = "GRADIENT_DETECTOR";
      this_time_crit.d_time = time;
      this_time_crit.d_tag_criteria.push_back(this_tag_crit);
      d_cur_time_criteria = d_time_criteria.insert(this_time_crit).first;
   } else {
      bool grad_detect_already_on = false;
      for (std::vector<TagCriteria>::const_iterator i = existing->d_tag_criteria.begin();
           i != existing->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "GRADIENT_DETECTOR") {
            grad_detect_already_on = true;
            break;
         }
      }
      if (!grad_detect_already_on) {
         TagCriteria this_tag_crit;
         this_tag_crit.d_tagging_method = "GRADIENT_DETECTOR";
         std::vector<TagCriteria>& this_tag_criteria =
            const_cast<std::vector<TagCriteria>&>(existing->d_tag_criteria);
         this_tag_criteria.push_back(this_tag_crit);
      }
   }
}

void
StandardTagAndInitialize::turnOffGradientDetector(
   double time)
{
   for (std::set<TimeTagCriteria, time_tag_criteria_less>::iterator i = d_time_criteria.begin();
        i != d_time_criteria.end(); ++i) {
      if (i->d_time <= time) {
         std::vector<TagCriteria>& tag_crits =
            const_cast<std::vector<TagCriteria>&>(i->d_tag_criteria);
         std::vector<TagCriteria>::iterator j = tag_crits.begin();
         while (j != tag_crits.end()) {
            if (j->d_tagging_method == "GRADIENT_DETECTOR") {
               tag_crits.erase(j);
            } else {
               ++j;
            }
         }
         break;
      }
   }
}

void
StandardTagAndInitialize::turnOnRichardsonExtrapolation(
   double time)
{
   if (!d_tag_strategy) {
      TBOX_ERROR("StandardTagAndInitialize::turnOnRichardsonExtrapolation\n"
         << "A tagging strategy must be defined if\n"
         << "Richardson extrapolation is used.\n");
   }

   TimeTagCriteria search_for;
   search_for.d_time = time;
   std::set<TimeTagCriteria, time_tag_criteria_less>::iterator existing =
      d_time_criteria.find(search_for);
   if (existing == d_time_criteria.end()) {
      TimeTagCriteria this_time_crit;
      TagCriteria this_tag_crit;
      this_tag_crit.d_tagging_method = "RICHARDSON_EXTRAPOLATION";
      this_time_crit.d_time = time;
      this_time_crit.d_tag_criteria.push_back(this_tag_crit);
      d_cur_time_criteria = d_time_criteria.insert(this_time_crit).first;
   } else {
      bool rich_extrap_already_on = false;
      for (std::vector<TagCriteria>::const_iterator i = existing->d_tag_criteria.begin();
           i != existing->d_tag_criteria.end(); ++i) {
         if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
            rich_extrap_already_on = true;
            break;
         }
      }
      if (!rich_extrap_already_on) {
         TagCriteria this_tag_crit;
         this_tag_crit.d_tagging_method = "RICHARDSON_EXTRAPOLATION";
         std::vector<TagCriteria>& this_tag_criteria =
            const_cast<std::vector<TagCriteria>&>(existing->d_tag_criteria);
         this_tag_criteria.push_back(this_tag_crit);
      }
   }
}

void
StandardTagAndInitialize::turnOffRichardsonExtrapolation(
   double time)
{
   for (std::set<TimeTagCriteria, time_tag_criteria_less>::iterator i = d_time_criteria.begin();
        i != d_time_criteria.end(); ++i) {
      if (i->d_time <= time) {
         std::vector<TagCriteria>& tag_crits =
            const_cast<std::vector<TagCriteria>&>(i->d_tag_criteria);
         std::vector<TagCriteria>::iterator j = tag_crits.begin();
         while (j != tag_crits.end()) {
            if (j->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
               tag_crits.erase(j);
            } else {
               ++j;
            }
         }
         break;
      }
   }
}

/*
 *************************************************************************
 * Find the tagging criteria for the supplied time and cycle.  If there
 * are both cycle and time criteria for this cycle/time combination the
 * time criteria has precedence.
 *************************************************************************
 */
void
StandardTagAndInitialize::setCurrentTaggingCriteria(
   int cycle,
   double time)
{
   TBOX_ASSERT(!d_use_cycle_criteria || !d_use_time_criteria);

   // Only change tagging criteria if the problem has advanced to a new
   // cycle.  We specifically do not change tagging criteria level by level
   // due solely to changes in level simulation time.
   if (d_old_cycle != cycle) {
      d_old_cycle = cycle;

      // Cache information about the old (current) tagging criteria so that
      // it can be determined if the tagged boxes have changed due to a change
      // in tagging criteria.
      bool old_use_cycle_criteria = d_use_cycle_criteria;
      bool old_use_time_criteria = d_use_time_criteria;
      std::set<CycleTagCriteria, cycle_tag_criteria_less>::iterator old_cur_cycle_criteria =
         d_cur_cycle_criteria;
      std::set<TimeTagCriteria, time_tag_criteria_less>::iterator old_cur_time_criteria =
         d_cur_time_criteria;
      bool old_use_re = false;
      bool old_use_gd = false;
      bool old_use_rb = false;
      if (d_use_cycle_criteria) {
         for (std::vector<TagCriteria>::const_iterator i =
                 d_cur_cycle_criteria->d_tag_criteria.begin();
              i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
            if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
               old_use_re = true;
            } else if (i->d_tagging_method == "GRADIENT_DETECTOR") {
               old_use_gd = true;
            } else if (i->d_tagging_method == "REFINE_BOXES") {
               old_use_rb = true;
            }
         }
      } else if (d_use_time_criteria) {
         for (std::vector<TagCriteria>::const_iterator i =
                 d_cur_time_criteria->d_tag_criteria.begin();
              i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
            if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
               old_use_re = true;
            } else if (i->d_tagging_method == "GRADIENT_DETECTOR") {
               old_use_gd = true;
            } else if (i->d_tagging_method == "REFINE_BOXES") {
               old_use_rb = true;
            }
         }
      }

      // Now set the new tagging criteria.  If both cycle and time criteria
      // are valid at the current cycle/time, the time criteria takes
      // precedence.
      if (!d_cycle_criteria.empty()) {
         // Find the cycle criteria for the supplied cycle, if any.
         for ( ; d_cur_cycle_criteria != d_cycle_criteria.end();
               ++d_cur_cycle_criteria) {
            std::set<CycleTagCriteria, cycle_tag_criteria_less>::iterator next_cycle_criteria =
               d_cur_cycle_criteria;
            ++next_cycle_criteria;
            if (next_cycle_criteria == d_cycle_criteria.end()) {
               // This is the last cycle criteria so see if it's valid at this
               // cycle and quit.
               if (cycle >= d_cur_cycle_criteria->d_cycle) {
                  d_use_cycle_criteria = true;
                  d_use_time_criteria = false;
               }
               break;
            } else if (cycle >= d_cur_cycle_criteria->d_cycle &&
                       cycle < next_cycle_criteria->d_cycle) {
               // This is not the last cycle criteria so see if it's valid and
               // the next one isn't at this cycle.  If so then quit.
               d_use_cycle_criteria = true;
               d_use_time_criteria = false;
               break;
            } else if (cycle < d_cur_cycle_criteria->d_cycle) {
               // If this cycle criteria isn't valid yet then none of the
               // others will be either so just quit looking.
               break;
            }
         }
      }

      if (!d_time_criteria.empty()) {
         // Find the time criteria for the supplied time, if any.
         for ( ; d_cur_time_criteria != d_time_criteria.end();
               ++d_cur_time_criteria) {
            std::set<TimeTagCriteria, time_tag_criteria_less>::iterator next_time_criteria =
               d_cur_time_criteria;
            ++next_time_criteria;
            if (next_time_criteria == d_time_criteria.end()) {
               // This is the last time criteria so see if it's valid at
               // this time and quit.
               if (time >= d_cur_time_criteria->d_time) {
                  d_use_time_criteria = true;
                  d_use_cycle_criteria = false;
               }
               break;
            } else if (time >= d_cur_time_criteria->d_time &&
                       time < next_time_criteria->d_time) {
               // This is not the last time criteria so see if it's valid
               // and the next one isn't at this time.  If so then quit.
               d_use_time_criteria = true;
               d_use_cycle_criteria = false;
               break;
            } else if (time < d_cur_time_criteria->d_time) {
               // If this time criteria isn't valid yet then none of the
               // others will be either so just quit looking.
               break;
            }
         }
      }

      // Now summarize information about the new tagging criteria for
      // comparison with the information cached at the beginning of this
      // method about the old tagging criteria.
      bool new_use_re = false;
      bool new_use_gd = false;
      bool new_use_rb = false;
      if (d_use_cycle_criteria) {
         for (std::vector<TagCriteria>::const_iterator i =
                 d_cur_cycle_criteria->d_tag_criteria.begin();
              i != d_cur_cycle_criteria->d_tag_criteria.end(); ++i) {
            if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
               new_use_re = true;
            } else if (i->d_tagging_method == "GRADIENT_DETECTOR") {
               new_use_gd = true;
            } else if (i->d_tagging_method == "REFINE_BOXES") {
               new_use_rb = true;
            }
         }
      } else if (d_use_time_criteria) {
         for (std::vector<TagCriteria>::const_iterator i =
                 d_cur_time_criteria->d_tag_criteria.begin();
              i != d_cur_time_criteria->d_tag_criteria.end(); ++i) {
            if (i->d_tagging_method == "RICHARDSON_EXTRAPOLATION") {
               new_use_re = true;
            } else if (i->d_tagging_method == "GRADIENT_DETECTOR") {
               new_use_gd = true;
            } else if (i->d_tagging_method == "REFINE_BOXES") {
               new_use_rb = true;
            }
         }
      }

      // Compare the old and new tagging criteria to determine if the tagged
      // boxes changed.
      if ((old_use_re != new_use_re) || (old_use_gd != new_use_gd) ||
          (old_use_rb != new_use_rb)) {
         // If one of the tagging methods which was used is now not used or
         // vice-versa, then the tagged boxes have changed.
         d_boxes_changed = true;
      } else {
         // The tagging methods are the same as they were last cycle.
         if (new_use_re || new_use_gd) {
            // If we're using either Richardson extrapolation or gradient
            // detector we must assume that the boxes have changed.
            d_boxes_changed = true;
         } else if ((old_use_cycle_criteria != d_use_cycle_criteria) ||
                    (old_use_time_criteria != d_use_time_criteria)) {
            // Assume that switching from a cycle to a time criteria or a time
            // to a cycle criteria results in different boxes.
            d_boxes_changed = true;
         } else if (d_use_cycle_criteria) {
            // We're still using a cycle tagging criteria.  If we're using a
            // different cycle tagging criteria assume that the boxes have
            // changed.
            if (d_cur_cycle_criteria != old_cur_cycle_criteria) {
               d_boxes_changed = true;
            } else {
               d_boxes_changed = false;
            }
         } else if (d_use_time_criteria) {
            // We're still using a time tagging criteria.  If we're using a
            // different time tagging criteria assume that the boxes have
            // changed.
            if (d_cur_time_criteria != old_cur_time_criteria) {
               d_boxes_changed = true;
            } else {
               d_boxes_changed = false;
            }
         }
      }
   }
}

void
StandardTagAndInitialize::processHierarchyBeforeAddingNewLevel(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::BoxLevel>& new_box_level)
{
   d_tag_strategy->processHierarchyBeforeAddingNewLevel(hierarchy,
      level_number,
      new_box_level);
}

void
StandardTagAndInitialize::processLevelBeforeRemoval(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::PatchLevel>& old_level)
{
   d_tag_strategy->processLevelBeforeRemoval(hierarchy,
      level_number,
      old_level);
}


/*
 *************************************************************************
 * Pass to StandardTagAndInitStrategy to check user tags on a tagged
 * level.  Only calls into the strategy if using gradient detector or
 * Richardson extrapolation at the curren time and cycle.
 *************************************************************************
 */
void
StandardTagAndInitialize::checkUserTagData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int regrid_cycle,
   const double regrid_time,
   const int tag_index)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(tag_index >= 0);

   if (usesGradientDetector(regrid_cycle, regrid_time) ||
       usesRichardsonExtrapolation(regrid_cycle, regrid_time)) {

      TBOX_ASSERT(d_tag_strategy != 0);

      d_tag_strategy->checkUserTagData(hierarchy,
         level_number,
         tag_index);

   }
}

/*
 *************************************************************************
 * Pass to StandardTagAndInitStrategy to check saved tags on a new level.
 *************************************************************************
 */
void
StandardTagAndInitialize::checkNewLevelTagData(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const int tag_index)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(tag_index >= 0);
   TBOX_ASSERT(d_tag_strategy != 0);

   d_tag_strategy->checkNewLevelTagData(hierarchy,
      level_number,
      tag_index);

}




static int GCD(
   const int a,
   const int b)
{
   int at = tbox::MathUtilities<int>::Min(a, b);
   int bt = tbox::MathUtilities<int>::Max(a, b);

   if (at == 0 || bt == 0) return bt;

   at = (at > 0 ? at : -at);
   bt = (bt > 0 ? bt : -bt);

   int r0 = at;
   int r1 = at;
   int r2 = bt;

   while (!(r2 == 0)) {
      r0 = r1;
      r1 = r2;
      int q = r0 / r1;
      r2 = r0 - r1 * q;
   }

   return r1;
}

}
}
