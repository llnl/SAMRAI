/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Coarsening schedule for data transfer between AMR levels
 *
 ************************************************************************/
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"
#include "SAMRAI/xfer/CoarsenCopyTransaction.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"

#include <vector>

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Initialization for static data members.
 *
 *************************************************************************
 */

std::string CoarsenSchedule::s_schedule_generation_method = "DLBG";
bool CoarsenSchedule::s_extra_debug = false;
bool CoarsenSchedule::s_barrier_and_time = false;
bool CoarsenSchedule::s_read_static_input = false;

std::shared_ptr<tbox::Timer> CoarsenSchedule::t_coarsen_schedule;
std::shared_ptr<tbox::Timer> CoarsenSchedule::t_coarsen_data;
std::shared_ptr<tbox::Timer> CoarsenSchedule::t_gen_sched_n_squared;
std::shared_ptr<tbox::Timer> CoarsenSchedule::t_gen_sched_dlbg;
std::shared_ptr<tbox::Timer> CoarsenSchedule::t_coarse_data_fill;

tbox::StartupShutdownManager::Handler
CoarsenSchedule::s_initialize_finalize_handler(
   CoarsenSchedule::initializeCallback,
   0,
   0,
   CoarsenSchedule::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 * ************************************************************************
 *
 * Static function to set box intersection algorithm for schedules.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::setScheduleGenerationMethod(
   const std::string& method)
{
   if (!((method == "ORIG_NSQUARED") ||
         (method == "DLBG"))) {
      TBOX_ERROR("CoarsenSchedule::setScheduleGenerationMethod\n"
         << "Given method std::string "
         << method << " is invalid.\n Options are\n"
         << "'ORIG_NSQUARED' and 'DLBG'."
         << std::endl);
   }

   s_schedule_generation_method = method;
}

/*
 * ************************************************************************
 *
 * Create a coarsening schedule that transfers data from the source
 * patch data components of the fine level into the destination patch
 * data components of the coarse level.  If the coarsening operators
 * require data in ghost cells on the source level, then those ghost
 * cells must be filled before this call.
 *
 * ************************************************************************
 */

CoarsenSchedule::CoarsenSchedule(
   const std::shared_ptr<hier::PatchLevel>& crse_level,
   const std::shared_ptr<hier::PatchLevel>& fine_level,
   const std::shared_ptr<CoarsenClasses>& coarsen_classes,
   const std::shared_ptr<CoarsenTransactionFactory>& transaction_factory,
   CoarsenPatchStrategy* patch_strategy,
   bool fill_coarse_data):
   d_number_coarsen_items(0),
   d_coarsen_items(0),
   d_crse_level(crse_level),
   d_fine_level(fine_level),
   d_coarsen_patch_strategy(patch_strategy),
   d_transaction_factory(transaction_factory),
   d_ratio_between_levels(crse_level->getDim(),
                          0,
                          crse_level->getGridGeometry()->getNumberBlocks()),
   d_fill_coarse_data(fill_coarse_data)
{
   TBOX_ASSERT(crse_level);
   TBOX_ASSERT(fine_level);
   TBOX_ASSERT(coarsen_classes);
   TBOX_ASSERT(transaction_factory);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*crse_level, *fine_level);

   getFromInput();

   if (s_barrier_and_time) {
      t_coarsen_schedule->start();
   }

   if (s_extra_debug) {
      tbox::plog << "CoarsenSchedule::CoarsenSchedule " << this << " entered" << std::endl;
   }

   const tbox::Dimension& dim(crse_level->getDim());

   /*
    * Compute ratio between fine and coarse levels and then check for
    * correctness.
    */

   const hier::IntVector& fine(d_fine_level->getRatioToLevelZero());
   const hier::IntVector& crse(d_crse_level->getRatioToLevelZero());

   const size_t nblocks = d_crse_level->getGridGeometry()->getNumberBlocks();

   for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
      for (unsigned int i = 0; i < dim.getValue(); ++i) {
         if (fine(b,i) > 1) {
            d_ratio_between_levels(b,i) = fine(b,i) / crse(b,i);
         } else {
            d_ratio_between_levels(b,i) = tbox::MathUtilities<int>::Abs(crse(
                     b,i) / fine(b,i));
         }
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (dim > tbox::Dimension(1)) {
      for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
         for (unsigned int i = 0; i < dim.getValue(); ++i) {
            if (d_ratio_between_levels(b,i)
                * d_ratio_between_levels(b,(i + 1) % dim.getValue()) < 0) {
               TBOX_ASSERT((d_ratio_between_levels(b,i) == 1) ||
                  (d_ratio_between_levels(b,(i + 1) % dim.getValue()) == 1));
            }
         }
      }
   }
#endif

   setCoarsenItems(coarsen_classes);
   initialCheckCoarsenClassItems();

   /*
    * Set up refine schedules to transfer coarsened data and to fill temporary
    * coarse level data before coarsening operations, if needed.  Then,
    * generate communication schedules to transfer data.
    */

   setupRefineAlgorithm();

   generateSchedule();

   if (s_extra_debug) {
      tbox::plog << "CoarsenSchedule::CoarsenSchedule " << this << " returning" << std::endl;
   }

   if (s_barrier_and_time) {
      t_coarsen_schedule->barrierAndStop();
   }
}

/*
 * ************************************************************************
 *
 * The destructor for the coarsen schedule class implicitly deallocates
 * all of the data associated with the communication schedule.
 *
 * ************************************************************************
 */

CoarsenSchedule::~CoarsenSchedule()
{
   clearCoarsenItems();
   delete[] d_coarsen_items;
}

/*
 * ***********************************************************************
 *
 * Read static member data from input database once.
 *
 * ***********************************************************************
 */
void
CoarsenSchedule::getFromInput()
{
   if (!s_read_static_input) {
      s_read_static_input = true;
      std::shared_ptr<tbox::Database> idb(
         tbox::InputManager::getInputDatabase());
      if (idb && idb->isDatabase("CoarsenSchedule")) {
         std::shared_ptr<tbox::Database> csdb(
            idb->getDatabase("CoarsenSchedule"));
         s_extra_debug = csdb->getBoolWithDefault("DEV_extra_debug", s_extra_debug);
         s_barrier_and_time =
            csdb->getBoolWithDefault("DEV_barrier_and_time", s_barrier_and_time);
      }
   }
}

/*
 * ***********************************************************************
 *
 * Reset schedule with new set of coarsen items.
 *
 * ***********************************************************************
 */

void
CoarsenSchedule::reset(
   const std::shared_ptr<CoarsenClasses>& coarsen_classes)
{
   TBOX_ASSERT(coarsen_classes);

   setCoarsenItems(coarsen_classes);

   setupRefineAlgorithm();

   if (d_fill_coarse_data) {
      t_coarse_data_fill->start();
      d_precoarsen_refine_algorithm->resetSchedule(
         d_precoarsen_refine_schedule);
      t_coarse_data_fill->stop();
   }

}

/*
 * ************************************************************************
 *
 * Execute the stored communication schedule that copies data into the
 * the destination patch data components of the destination level from
 * the source patch data components of the source level.  The steps
 * to the algorithm are as follows:
 *
 *      (1) Allocate the source space on the temporary patch level.
 *      (2) Coarsen the data from the fine patch level to the temporary
 *          patch level (local operation).
 *      (3) Copy data from the source space of the temporary patch
 *          level into the destination space of the destination patch
 *          level (requires interprocessor communication).
 *      (4) Deallocate the source space on the temporary patch level.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::coarsenData() const
{
   if (s_extra_debug) {
      tbox::plog << "CoarsenSchedule::coarsenData " << this << " entered" << std::endl;
   }
   if (s_barrier_and_time) {
      t_coarsen_data->barrierAndStart();
   }

   /*
    * Allocate the source data space on the temporary patch level.
    * We do not know the current time, so set it to zero.  It should
    * not matter, since the copy routines do not require that
    * the time markers match.
    */

   d_temp_crse_level->allocatePatchData(d_sources, 0.0);

   if (d_fill_coarse_data) {
      t_coarse_data_fill->start();
      d_precoarsen_refine_schedule->fillData(0.0);
      t_coarse_data_fill->stop();
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }

   /*
    * Coarsen the data from the sources on the fine data level into the
    * sources on the temporary data level
    */

   coarsenSourceData(d_coarsen_patch_strategy);

   /*
    * Copy data from the source interiors of the temporary patch level
    * into the destination interiors of the destination patch level.
    */

   d_schedule->communicate();

   /*
    * Deallocate the source data in the temporary patch level.
    */

   d_temp_crse_level->deallocatePatchData(d_sources);

   if (s_extra_debug) {
      tbox::plog << "CoarsenSchedule::coarsenData " << this << " returning" << std::endl;
   }

   if (s_barrier_and_time) {
      t_coarsen_data->stop();
   }
}

/*
 * ************************************************************************
 *
 * Generate the temporary coarse level by coarsening the fine patch
 * level boxes.  Note that no patch data components are allocated until
 * they are needed during the coarsening operation.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::generateTemporaryLevel()
{
   const tbox::Dimension& dim(d_crse_level->getDim());

   d_temp_crse_level.reset(new hier::PatchLevel(dim));
   d_temp_crse_level->setCoarsenedPatchLevel(d_fine_level,
      d_ratio_between_levels);
   d_temp_crse_level->setLevelNumber(d_crse_level->getLevelNumber());
   d_temp_crse_level->setNextCoarserHierarchyLevelNumber(
      d_crse_level->getLevelNumber());

   const hier::IntVector max_ghosts = getMaxGhostsToGrow();
   hier::IntVector min_width(dim);
   TBOX_ASSERT(d_crse_level->getBoxLevel()->getRefinementRatio() != 0);
   if (d_crse_level->getBoxLevel()->getRefinementRatio() >
       hier::IntVector::getZero(dim)) {
      min_width =
         (d_fine_level->getBoxLevel()->getRefinementRatio()
          / d_crse_level->getBoxLevel()->getRefinementRatio())
         * max_ghosts;
   } else {
      TBOX_ASSERT(d_fine_level->getBoxLevel()->getRefinementRatio() >=
         hier::IntVector::getOne(dim));
      min_width =
         (-d_crse_level->getBoxLevel()->getRefinementRatio()
          / d_fine_level->getBoxLevel()->getRefinementRatio())
         * max_ghosts;
   }

   const hier::IntVector transpose_width =
      hier::Connector::convertHeadWidthToBase(
         d_crse_level->getBoxLevel()->getRefinementRatio(),
         d_fine_level->getBoxLevel()->getRefinementRatio(),
         min_width);

   const hier::Connector& coarse_to_fine =
      d_crse_level->findConnectorWithTranspose(*d_fine_level,
         transpose_width,
         min_width,
         hier::CONNECTOR_IMPLICIT_CREATION_RULE,
         true);

   /*
    * Generate temporary BoxLevel and Connectors.
    */

   /*
    * Compute d_coarse_to_temp and its transpose.
    *
    * We use the fact that d_temp_crse_level patches are numbered just
    * like the fine level patches.  The Connectors between coarse and
    * temp are very similar to those between coarse and fine.
    */
   d_coarse_to_temp.reset(new hier::Connector(coarse_to_fine));
   d_coarse_to_temp->setBase(*d_crse_level->getBoxLevel());
   d_coarse_to_temp->setHead(*d_temp_crse_level->getBoxLevel());
   d_coarse_to_temp->setWidth(coarse_to_fine.getConnectorWidth(), true);
   d_coarse_to_temp->coarsenLocalNeighbors(d_ratio_between_levels);
   /*
    * temp_to_coarse is a Connector from a coarsened version of fine to
    * coarse.  Therefore it has the same neighborhoods as fine_to_coarse
    * but it's base, head and width are different.  So first assign
    * fine_to_coarse to temp_to_coarse which will properly set the
    * neighborhoods.  Then initialize it with the proper base/head/width
    * keeping the neighborhoods that we just set.
    */
   hier::Connector* temp_to_coarse =
      new hier::Connector(coarse_to_fine.getTranspose());
   temp_to_coarse->setBase(*d_temp_crse_level->getBoxLevel());
   temp_to_coarse->setHead(coarse_to_fine.getBase());
   temp_to_coarse->setWidth(coarse_to_fine.getConnectorWidth(), true);
   const hier::IntVector& one_vector(hier::IntVector::getOne(dim));
   d_coarse_to_temp->shrinkWidth(one_vector);
   temp_to_coarse->shrinkWidth(one_vector);
   d_coarse_to_temp->setTranspose(temp_to_coarse, true);
}

/*
 * ************************************************************************
 *
 * Set up refine algorithms to transfer coarsened data and to fill
 * temporary coarse level before performing coarsening operations,
 * if necessary.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::setupRefineAlgorithm()
{
   if (d_fill_coarse_data) {
      t_coarse_data_fill->barrierAndStart();

      d_precoarsen_refine_algorithm.reset(new RefineAlgorithm());

      for (size_t ici = 0; ici < d_number_coarsen_items; ++ici) {
         const int src_id = d_coarsen_items[ici]->d_src;
         d_precoarsen_refine_algorithm->registerRefine(src_id,
            src_id,
            src_id,
            std::shared_ptr<hier::RefineOperator>());
      }

      t_coarse_data_fill->stop();
   }

}

/*
 * ************************************************************************
 *
 * Generate communication schedule that copies source patch data
 * from the temporary level into the destination patch data of the
 * destination (coarse) level.  The source and destination
 * spaces may be the same.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::generateSchedule()
{

   /*
    * Set up coarsened version of fine level for temporary data storage.
    * Next, create refine algorithm if needed to fill temporary coarse
    * level before coarsen operations occur.  Then, create empty schedule
    * that will hold transactions for moving data.  Finally, generate
    * schedule based on chosen generation method.
    */
   generateTemporaryLevel();

   if (d_fill_coarse_data) {
      t_coarse_data_fill->barrierAndStart();
      d_precoarsen_refine_schedule =
         d_precoarsen_refine_algorithm->createSchedule(d_temp_crse_level,
            d_crse_level, 0);
      t_coarse_data_fill->stop();
   }

   d_schedule.reset(new tbox::Schedule());
   d_schedule->setTimerPrefix("xfer::CoarsenSchedule");

   if (s_schedule_generation_method == "ORIG_NSQUARED") {

      generateScheduleNSquared();

   } else if (s_schedule_generation_method == "DLBG") {

      generateScheduleDLBG();

   } else {

      TBOX_ERROR("Internal CoarsenSchedule error..."
         << "\n unrecognized schedule generation option: "
         << s_schedule_generation_method << std::endl);

   }

}

/*
 *************************************************************************
 *
 * This version of the schedule generation procedure uses the original
 * SAMRAI N^2 algorithms to construct communication schedules.  Here,
 * we loop over all of the patches on the source and destination levels.
 * check to see whether source or destination is local to this processor.
 * If not, then skip over schedule construction operations.
 *
 *************************************************************************
 */

void
CoarsenSchedule::generateScheduleNSquared()
{

   t_gen_sched_n_squared->start();

   const int dst_npatches = d_crse_level->getGlobalNumberOfPatches();
   const int src_npatches = d_temp_crse_level->getGlobalNumberOfPatches();

   const hier::ProcessorMapping& dst_mapping =
      d_crse_level->getProcessorMapping();
   const hier::ProcessorMapping& src_mapping =
      d_temp_crse_level->getProcessorMapping();

   hier::BoxContainer::const_iterator crse_itr_dp =
      d_crse_level->getBoxes().begin();
   for (int dp = 0; dp < dst_npatches; ++dp, ++crse_itr_dp) {

      const hier::Box dst_box(*crse_itr_dp,
                              hier::LocalId(dp),
                              dst_mapping.getProcessorAssignment(dp));

      hier::BoxContainer::const_iterator crse_itr_sp =
         d_temp_crse_level->getBoxes().begin();
      for (int sp = 0; sp < src_npatches; ++sp, ++crse_itr_sp) {

         const hier::Box src_box(
            *crse_itr_sp,
            hier::LocalId(sp),
            src_mapping.getProcessorAssignment(sp));

         if (dst_mapping.isMappingLocal(dp)
             || src_mapping.isMappingLocal(sp)) {

            constructScheduleTransactions(d_crse_level, dst_box,
               d_temp_crse_level, src_box);

         }  // if either source or destination patch is local

      } // loop over source patches

   } // loop over destination patches

   t_gen_sched_n_squared->stop();

}

/*
 *************************************************************************
 *************************************************************************
 */
void
CoarsenSchedule::generateScheduleDLBG()
{
   t_gen_sched_dlbg->start();

   /*
    * Construct sending transactions for local src Boxes.
    */
   /*
    * Restructure the temp_to_coarse edge data to arange neighbors by the
    * coarse boxes, as required to match the transaction ordering on the
    * receiving processors.  At the same time, shift temp-coarse pairs to
    * make the coarse shifts zero.
    */
   FullNeighborhoodSet temp_eto_coarse_bycoarse;
   restructureNeighborhoodSetsByDstNodes(temp_eto_coarse_bycoarse,
      d_coarse_to_temp->getTranspose());

   for (FullNeighborhoodSet::const_iterator ei = temp_eto_coarse_bycoarse.begin();
        ei != temp_eto_coarse_bycoarse.end(); ++ei) {

      /*
       * coarse_box can be remote (by definition of FullNeighborhoodSet).
       * local_temp_boxes are the local source boxes that contribute data
       * to box.
       */
      const hier::Box& coarse_box = ei->first;
      const hier::BoxContainer& local_temp_boxes = ei->second;
      TBOX_ASSERT(!coarse_box.isPeriodicImage());

      /*
       * Construct transactions for data going from local source boxes
       * to remote coarse boxes.
       */
      for (hier::BoxContainer::const_iterator ni =
              local_temp_boxes.begin();
           ni != local_temp_boxes.end(); ++ni) {
         const hier::Box& temp_box = *ni;
         if (temp_box.getOwnerRank() ==
             coarse_box.getOwnerRank()) {
            /*
             * Disregard local coarse_box to avoid duplicating same
             * transactions created by the second loop below.
             */
            continue;
         }
         constructScheduleTransactions(d_crse_level,
            coarse_box,
            d_temp_crse_level,
            temp_box);
      }

   }

   /*
    * Construct receiving transactions for local dst boxes.
    */
   const hier::BoxLevel& coarse_box_level = *d_crse_level->getBoxLevel();
   for (hier::Connector::ConstNeighborhoodIterator ei = d_coarse_to_temp->begin();
        ei != d_coarse_to_temp->end(); ++ei) {

      const hier::BoxId& dst_gid = *ei;
      const hier::Box& dst_box =
         *coarse_box_level.getBoxStrict(dst_gid);

      for (hier::Connector::ConstNeighborIterator ni = d_coarse_to_temp->begin(ei);
           ni != d_coarse_to_temp->end(ei); ++ni) {
         const hier::Box& src_box = *ni;

         constructScheduleTransactions(d_crse_level,
            dst_box,
            d_temp_crse_level,
            src_box);

      }

   }

   t_gen_sched_dlbg->stop();

}

/*
 ***********************************************************************
 * This method does 2 important things to the src_to_dst edges:
 *
 * 1. It puts the edge data in dst-major order so the src owners can
 * easily loop through the dst-src edges in the same order that dst
 * owners see them.  Transactions must have the same order on the
 * sending and receiving processors.
 *
 * 2. It shifts periodic image dst boxes back to the zero-shift position,
 * and applies a similar shift to src boxes so that the overlap is
 * unchanged.  The constructScheduleTransactions method requires all
 * shifts to be absorbed in the src box.
 ***********************************************************************
 */
void
CoarsenSchedule::restructureNeighborhoodSetsByDstNodes(
   FullNeighborhoodSet& full_inverted_edges,
   const hier::Connector& src_to_dst) const
{
   const tbox::Dimension& dim(d_crse_level->getDim());

   const hier::BoxLevel& src_box_level = src_to_dst.getBase();
   const hier::IntVector& src_ratio(src_to_dst.getBase().getRefinementRatio());
   const hier::IntVector& dst_ratio(src_to_dst.getHead().getRefinementRatio());

   const hier::PeriodicShiftCatalog& shift_catalog =
      src_box_level.getGridGeometry()->getPeriodicShiftCatalog();

   /*
    * These are the counterparts to shifted dst boxes and unshifted src boxes.
    */
   hier::Box shifted_box(dim), unshifted_nabr(dim);
   full_inverted_edges.clear();
   for (hier::Connector::ConstNeighborhoodIterator ci = src_to_dst.begin();
        ci != src_to_dst.end();
        ++ci) {
      const hier::Box& box = *src_box_level.getBoxStrict(*ci);
      for (hier::Connector::ConstNeighborIterator na = src_to_dst.begin(ci);
           na != src_to_dst.end(ci); ++na) {
         const hier::Box& nabr = *na;
         if (nabr.isPeriodicImage()) {
            shifted_box.initialize(
               box,
               shift_catalog.getOppositeShiftNumber(nabr.getPeriodicId()),
               src_ratio,
               shift_catalog);
            unshifted_nabr.initialize(
               nabr,
               shift_catalog.getZeroShiftNumber(),
               dst_ratio,
               shift_catalog);

            full_inverted_edges[unshifted_nabr].insert(shifted_box);
         } else {
            full_inverted_edges[nabr].insert(box);
         }
      }
   }
}

/*
 **************************************************************************
 * Calculate the max ghost cell width to grow boxes to check for
 * overlaps.  Given in source (fine) level's index space.
 **************************************************************************
 */

hier::IntVector
CoarsenSchedule::getMaxGhostsToGrow() const
{
   const tbox::Dimension& dim(d_crse_level->getDim());
   std::shared_ptr<hier::PatchDescriptor> pd(
      d_temp_crse_level->getPatchDescriptor());

   /*
    * Box, face and side elements of adjacent cells overlap even though
    * the cells do not overlap.  Therefore, we always grow at least one
    * cell catch overlaps of box, face and side elements.
    */
   hier::IntVector gcw(dim, 1);

   for (size_t ici = 0; ici < d_number_coarsen_items; ++ici) {

      const int src_id = d_coarsen_items[ici]->d_src;
      gcw.max(pd->getPatchDataFactory(src_id)->getGhostCellWidth());

      hier::IntVector gcw1 = d_coarsen_items[ici]->d_gcw_to_coarsen;
      if (d_coarsen_items[ici]->d_opcoarsen) {
         gcw1 += d_coarsen_items[ici]->d_opcoarsen->getStencilWidth(dim);
      }
      gcw.max(gcw1);
   }

   return gcw;
}

/*
 *************************************************************************
 *
 * Private utility function that constructs schedule transactions that
 * move data from source patch on source level to destination patch
 * on destination level.
 *
 *************************************************************************
 */

void
CoarsenSchedule::constructScheduleTransactions(
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const hier::Box& dst_box,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const hier::Box& src_box)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);

   const tbox::Dimension& dim(d_crse_level->getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(dim,
      *dst_level,
      *src_level,
      dst_box,
      src_box);

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));
   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   if (s_extra_debug) {
      tbox::plog << "CoarsenSchedule::constructScheduleTransactions:"
                 << "\n  src: L" << src_level->getLevelNumber()
                 << "R" << src_level->getRatioToLevelZero()
                 << " / " << src_box << ""
                 << "\n  dst: L" << dst_level->getLevelNumber()
                 << "R" << dst_level->getRatioToLevelZero()
                 << " / " << dst_box
                 << std::endl;
   }

   std::shared_ptr<hier::PatchDescriptor> dst_patch_descriptor(
      dst_level->getPatchDescriptor());
   std::shared_ptr<hier::PatchDescriptor> src_patch_descriptor(
      src_level->getPatchDescriptor());

   const int num_equiv_classes =
      d_coarsen_classes->getNumberOfEquivalenceClasses();

   const hier::PeriodicShiftCatalog& shift_catalog =
      src_level->getGridGeometry()->getPeriodicShiftCatalog();

   /*
    * Calculate the shift and the shifted source box.
    */
   hier::IntVector src_shift(dim, 0);
   hier::IntVector dst_shift(dim, 0);
   hier::Box unshifted_src_box = src_box;
   hier::Box unshifted_dst_box = dst_box;
   const hier::BlockId& dst_block_id = dst_box.getBlockId();
   const hier::BlockId& src_block_id = src_box.getBlockId();
   if (src_box.isPeriodicImage()) {
      TBOX_ASSERT(!dst_box.isPeriodicImage());
      src_shift = shift_catalog.shiftNumberToShiftDistance(
            src_box.getPeriodicId());
      src_shift *= src_level->getRatioToLevelZero();
      unshifted_src_box.shift(-src_shift);
   }
   if (dst_box.isPeriodicImage()) {
      TBOX_ASSERT(!src_box.isPeriodicImage());
      dst_shift = shift_catalog.shiftNumberToShiftDistance(
            dst_box.getPeriodicId());
      dst_shift *= dst_level->getRatioToLevelZero();
      unshifted_dst_box.shift(-dst_shift);
   }

   /*
    * Transformation initialized to src_shift with no rotation.
    * It will never be modified in single-block runs, nor in multiblock runs
    * when src_box and dst_box are on the same block.
    */
   hier::Transformation transformation(src_shift);

   /*
    * When src_box and dst_box are on different blocks
    * transformed_src_box is a representation of the source box in the
    * destination coordinate system.
    *
    * For all other cases, transformed_src_box is simply a copy of the
    * box from src_box.
    */
   hier::Box transformed_src_box(src_box);

   /*
    * When needed, transform the source box and determine if src and
    * dst touch at an enhance connectivity singularity.
    */
   if (src_block_id != dst_block_id) {

      std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
         d_crse_level->getGridGeometry());

      hier::Transformation::RotationIdentifier rotation =
         grid_geometry->getRotationIdentifier(dst_block_id,
            src_block_id);
      hier::IntVector offset(
         grid_geometry->getOffset(dst_block_id, src_block_id, d_crse_level->getLevelNumber()));

      transformation = hier::Transformation(rotation, offset,
            src_block_id, dst_block_id);
      transformation.transform(transformed_src_box);

#ifdef DEBUG_CHECK_ASSERTIONS
      if (grid_geometry->areSingularityNeighbors(dst_block_id, src_block_id)) {
         for (int nc = 0; nc < num_equiv_classes; ++nc) {
            const CoarsenClasses::Data& rep_item =
               d_coarsen_classes->getClassRepresentative(nc);

            TBOX_ASSERT(rep_item.d_var_fill_pattern->getStencilWidth() ==
               hier::IntVector::getZero(dim));

         }
      }
#endif
   }

   const int num_coarsen_items = d_coarsen_classes->getNumberOfCoarsenItems();
   std::vector<std::shared_ptr<tbox::Transaction> > transactions(
      num_coarsen_items);

   for (int nc = 0; nc < num_equiv_classes; ++nc) {

      if (s_extra_debug) {
         tbox::plog << " equivalent class " << nc << "/" << num_equiv_classes << std::endl;
      }
      const CoarsenClasses::Data& rep_item =
         d_coarsen_classes->getClassRepresentative(nc);

      const int rep_item_dst_id = rep_item.d_dst;
      const int rep_item_src_id = rep_item.d_src;

      std::shared_ptr<hier::PatchDataFactory> src_pdf(
         src_patch_descriptor->getPatchDataFactory(rep_item_src_id));
      std::shared_ptr<hier::PatchDataFactory> dst_pdf(
         dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id));

      const hier::IntVector& dst_gcw(dst_pdf->getGhostCellWidth());

      hier::Box dst_fill_box(unshifted_dst_box);
      dst_fill_box.grow(dst_gcw);

      hier::Box test_mask(dst_fill_box * transformed_src_box);
      if ((dst_gcw == constant_zero_intvector) &&
          dst_pdf->dataLivesOnPatchBorder() &&
          test_mask.empty()) {
         test_mask = dst_fill_box;
         test_mask.grow(constant_one_intvector);
         test_mask = test_mask * transformed_src_box;
      }
      hier::Box src_mask(test_mask);
      transformation.inverseTransform(src_mask);

      if (s_extra_debug) {
         tbox::plog << " dst_gcw = " << dst_gcw
                    << "\n dst_fill_box = " << dst_fill_box
                    << "\n test_mask = " << test_mask
                    << "\n src_mask (before += test_mask) = " << src_mask
                    << std::endl;
      }

      if (!src_mask.empty()) {
         // What does this block do?  Need comments!
         test_mask = unshifted_src_box;
         test_mask.grow(
            hier::IntVector::min(
               rep_item.d_gcw_to_coarsen,
               src_pdf->getGhostCellWidth()));
         src_mask += test_mask;
      }

      if (s_extra_debug) {
         tbox::plog << "\n src_mask (after += test_mask) = " << src_mask
                    << std::endl;
      }

      std::shared_ptr<hier::BoxOverlap> overlap(
         rep_item.d_var_fill_pattern->calculateOverlap(
            *dst_pdf->getBoxGeometry(unshifted_dst_box),
            *src_pdf->getBoxGeometry(unshifted_src_box),
            dst_box,
            src_mask,
            dst_fill_box,
            true, transformation));

      if (!overlap) {
         TBOX_ERROR("Internal CoarsenSchedule error..."
            << "\n Overlap is NULL for "
            << "\n src box = " << src_box
            << "\n dst box = " << dst_box
            << "\n src mask = " << src_mask << std::endl);
      }
      if (s_extra_debug) {
         tbox::plog << " Overlap:\n" << std::endl;
         overlap->print(tbox::plog);
      }

      if (!overlap->isOverlapEmpty()) {
         if (s_extra_debug) {
            tbox::plog << " Overlap FINITE." << std::endl;
         }
         for (std::list<int>::iterator l(d_coarsen_classes->getIterator(nc));
              l != d_coarsen_classes->getIteratorEnd(nc); ++l) {
            const CoarsenClasses::Data& item =
               d_coarsen_classes->getCoarsenItem(*l);
            TBOX_ASSERT(item.d_class_index == nc);
            TBOX_ASSERT(item.d_tag == *l);
            TBOX_ASSERT(&item == d_coarsen_items[*l]);

            const int citem_count = item.d_tag;
            transactions[citem_count] =
               d_transaction_factory->allocate(dst_level,
                  src_level,
                  overlap,
                  dst_box,
                  src_box,
                  d_coarsen_items,
                  citem_count);
         }
      } else {
         if (s_extra_debug) {
            tbox::plog << " Overlap empty." << std::endl;
         }
      }

   }  // iterate over all coarsen equivalence classes

   for (int i = 0; i < num_coarsen_items; ++i) {
      if (transactions[i]) {
         d_schedule->appendTransaction(transactions[i]);
      }
   }
}

/*
 * ************************************************************************
 *
 * Coarsen data from the source space on the fine patch level into the
 * source space on the coarse temporary patch level.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::coarsenSourceData(
   CoarsenPatchStrategy* patch_strategy) const
{
   if (patch_strategy) {
      patch_strategy->preprocessCoarsenLevel(
         *d_temp_crse_level,
         *d_fine_level);
#if defined(HAVE_RAJA)
      if (patch_strategy->needSynchronize()) {
         tbox::parallel_synchronize();
      }
#endif
   }

   /*
    * Loop over all local patches (fine and temp have the same mapping)
    */

   for (hier::PatchLevel::iterator p(d_fine_level->begin());
        p != d_fine_level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& fine_patch = *p;
      std::shared_ptr<hier::Patch> temp_patch(
         d_temp_crse_level->getPatch(fine_patch->getGlobalId()));

      const hier::Box& box = temp_patch->getBox();
      const hier::BlockId& block_id = box.getBlockId();
      hier::IntVector block_ratio(
         d_ratio_between_levels.getBlockVector(block_id));
      /*
       * Coarsen the fine space onto the temporary coarse space
       */
      if (patch_strategy) {
         patch_strategy->preprocessCoarsen(*temp_patch,
            *fine_patch, box, block_ratio);
#if defined(HAVE_RAJA)
         if (patch_strategy->needSynchronize()) {
            tbox::parallel_synchronize();
         }
#endif
      }

      for (size_t ici = 0; ici < d_number_coarsen_items; ++ici) {
         const CoarsenClasses::Data * const crs_item =
            d_coarsen_items[ici];
         if (crs_item->d_opcoarsen) {
            const int source_id = crs_item->d_src;
            crs_item->d_opcoarsen->coarsen(*temp_patch, *fine_patch,
               source_id, source_id,
               box, block_ratio);
         }
      }

      if (patch_strategy) {
         patch_strategy->setPostCoarsenSyncFlag();
#if defined(HAVE_RAJA)
         if (patch_strategy->needSynchronize()) {
            tbox::parallel_synchronize();
         }
#endif

         patch_strategy->postprocessCoarsen(*temp_patch,
            *fine_patch,
            box,
            block_ratio);
      }
   }

#if defined(HAVE_RAJA)
   if (!patch_strategy || patch_strategy->needSynchronize()) {
      tbox::parallel_synchronize();
   }
#endif

   if (patch_strategy) {
      patch_strategy->postprocessCoarsenLevel(
         *d_temp_crse_level,
         *d_fine_level);
#if defined(HAVE_RAJA)
      if (patch_strategy->needSynchronize()) {
         tbox::parallel_synchronize();
      }
#endif
   }

}

/*
 * ***********************************************************************
 *
 * Private utility function to set up local array of coarsen items.
 *
 * ***********************************************************************
 */

void
CoarsenSchedule::setCoarsenItems(
   const std::shared_ptr<CoarsenClasses>& coarsen_classes)
{

   clearCoarsenItems();

   d_coarsen_classes = coarsen_classes;

   d_number_coarsen_items = d_coarsen_classes->getNumberOfCoarsenItems();

   /*
    * Determine total number of coarsen items and set state of
    * component selector used to manage storage on temporary level.
    */
   d_sources.clrAllFlags();

   for (unsigned int nc = 0; nc < d_number_coarsen_items; ++nc) {
      const CoarsenClasses::Data& item = d_coarsen_classes->getCoarsenItem(nc);
      d_sources.setFlag(item.d_src);
   }

   /*
    * Allocate and initialize array of coarsen items.
    */

   if (!d_coarsen_items) {
      d_coarsen_items =
         new const CoarsenClasses::Data *[d_number_coarsen_items];
   }

   int ircount = 0;
   for (unsigned int nc = 0; nc < d_number_coarsen_items; ++nc) {
      d_coarsen_classes->getCoarsenItem(nc).d_tag = ircount;
      d_coarsen_items[ircount] = &(d_coarsen_classes->getCoarsenItem(nc));
      ++ircount;
   }

}

/*
 * ***********************************************************************
 *
 * Private utility function to check coarsen items in initial setup to
 * see whether source and destination patch data components have
 * sufficient ghost cell widths to satisfy the "ghost width to coarsen"
 * functionality described in the CoarsenAlgorithm class header.
 * Specifically, the destination data must have a ghost cell width at
 * least as large as the ghost cell width to coarsen.  The source data
 * must have a ghost cell width at least as large as the ghost cell
 * width to coarsen refined to the source (finer) level index space.
 * Other checks are also performed here by calling the
 * CoarsenClasses::itemIsValid() routine.
 *
 * ***********************************************************************
 */

void
CoarsenSchedule::initialCheckCoarsenClassItems() const
{
   const tbox::Dimension& dim(d_crse_level->getDim());

   std::shared_ptr<hier::PatchDescriptor> pd(
      d_crse_level->getPatchDescriptor());

   hier::IntVector user_gcw(dim, 0);
   if (d_coarsen_patch_strategy) {
      user_gcw = d_coarsen_patch_strategy->getCoarsenOpStencilWidth(dim);
   }

   for (size_t ici = 0; ici < d_number_coarsen_items; ++ici) {

      const CoarsenClasses::Data * const crs_item = d_coarsen_items[ici];

#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_coarsen_classes->itemIsValid(*crs_item, pd)) {
#endif

      const int dst_id = crs_item->d_dst;
      const int src_id = crs_item->d_src;

      const hier::IntVector& dst_gcw(
         pd->getPatchDataFactory(dst_id)->getGhostCellWidth());
      const hier::IntVector& src_gcw(
         pd->getPatchDataFactory(src_id)->getGhostCellWidth());

      if (crs_item->d_gcw_to_coarsen > dst_gcw) {
         TBOX_ERROR("Bad data given to CoarsenSchedule...\n"
            << "`ghost cell width to coarsen' specified in\n"
            << "registration of `Destination' patch data "
            << pd->mapIndexToName(dst_id)
            << " with CoarsenAlgorithm\n"
            << " is larger than ghost cell width of data \n"
            << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
            << "\n data ghost cell width = " << dst_gcw << std::endl);
      }

      if (d_ratio_between_levels * crs_item->d_gcw_to_coarsen > src_gcw) {
         TBOX_ERROR("Bad data given to CoarsenSchedule...\n"
            << "`Source' patch data " << pd->mapIndexToName(src_id)
            << " has ghost cell width too small to support the\n"
            << "`ghost cell width to coarsen' specified in"
            << " registration with CoarsenAlgorithm\n"
            << "data ghost cell width = " << src_gcw
            << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
            << "\nratio between levels = " << d_ratio_between_levels
            << "\n Thus, data ghost width must be >= "
            << d_ratio_between_levels * crs_item->d_gcw_to_coarsen
            << std::endl);
      }

      if (user_gcw > src_gcw) {
         TBOX_ERROR("Bad data given to CoarsenSchedule...\n"
            << "User supplied coarsen stencil width = "
            << user_gcw
            << "\nis larger than ghost cell width of `Source'\n"
            << "patch data " << pd->mapIndexToName(src_id)
            << " , which is " << src_gcw << std::endl);
      }

#ifdef DEBUG_CHECK_ASSERTIONS
   }
#endif

   }

}

/*
 * ************************************************************************
 *
 * Private utility function to clear array of coarsen items.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::clearCoarsenItems()
{
   if (d_coarsen_items) {
      for (size_t ici = 0; ici < d_number_coarsen_items; ++ici) {
         d_coarsen_items[ici] = 0;
      }
      d_number_coarsen_items = 0;
   }
}

/*
 **************************************************************************
 **************************************************************************
 */

void
CoarsenSchedule::setDeterministicUnpackOrderingFlag(bool flag)
{
   if (d_schedule) {
      d_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
   if (d_precoarsen_refine_schedule) {
      d_precoarsen_refine_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
}

void
CoarsenSchedule::setScheduleOpsStrategy(tbox::ScheduleOpsStrategy* strategy)
{
   if (d_schedule) {
      d_schedule->setScheduleOpsStrategy(strategy);
   }
   if (d_precoarsen_refine_schedule) {
      d_precoarsen_refine_schedule->setScheduleOpsStrategy(strategy);
   }
}

/*
 * ************************************************************************
 *
 * Print coarsen schedule data to the specified output stream.
 *
 * ************************************************************************
 */

void
CoarsenSchedule::printClassData(
   std::ostream& stream) const
{
   stream << "CoarsenSchedule::printClassData()" << std::endl;
   stream << "---------------------------------------" << std::endl;
   stream << "s_schedule_generation_method = "
          << s_schedule_generation_method << std::endl;
   stream << "d_fill_coarse_data = " << d_fill_coarse_data << std::endl;

   d_coarsen_classes->printClassData(stream);

   d_schedule->printClassData(stream);

   if (d_fill_coarse_data) {
      stream
      << "Printing pre-coarsen refine algorithm that fills data before coarsening...\n";
      d_precoarsen_refine_algorithm->printClassData(stream);
      stream
      << "Printing pre-coarsen refine schedule that fills data before coarsening...\n";
      d_precoarsen_refine_schedule->printClassData(stream);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
CoarsenSchedule::initializeCallback()
{
   t_coarsen_schedule = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule::CoarsenSchedule()");
   t_coarsen_data = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule::coarsenData()");
   t_gen_sched_n_squared = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule::generateScheduleNSquared()");
   t_gen_sched_dlbg = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule::generateScheduleDLBG()");
   t_coarse_data_fill = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule::coarse_data_fill");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */

void
CoarsenSchedule::finalizeCallback()
{
   t_coarsen_schedule.reset();
   t_coarsen_data.reset();
   t_gen_sched_n_squared.reset();
   t_gen_sched_dlbg.reset();
   t_coarse_data_fill.reset();
}

}
}
