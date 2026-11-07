/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Refine schedule for data transfer between AMR levels
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/xfer/RefineCopyTransaction.h"
#include "SAMRAI/xfer/RefineScheduleConnectorWidthRequestor.h"
#include "SAMRAI/xfer/RefineTimeTransaction.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/NVTXUtilities.h"
#include "SAMRAI/tbox/Collectives.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

bool RefineSchedule::s_extra_debug = false;
bool RefineSchedule::s_barrier_and_time = false;
bool RefineSchedule::s_read_static_input = false;

std::shared_ptr<tbox::Timer> RefineSchedule::t_refine_schedule;
std::shared_ptr<tbox::Timer> RefineSchedule::t_fill_data;
std::shared_ptr<tbox::Timer> RefineSchedule::t_fill_data_nonrecursive;
std::shared_ptr<tbox::Timer> RefineSchedule::t_fill_data_recursive;
std::shared_ptr<tbox::Timer> RefineSchedule::t_fill_physical_boundaries;
std::shared_ptr<tbox::Timer> RefineSchedule::t_fill_singularity_boundaries;
std::shared_ptr<tbox::Timer> RefineSchedule::t_refine_scratch_data;
std::shared_ptr<tbox::Timer> RefineSchedule::t_finish_sched_const;
std::shared_ptr<tbox::Timer> RefineSchedule::t_finish_sched_const_recurse;
std::shared_ptr<tbox::Timer> RefineSchedule::t_gen_comm_sched;
std::shared_ptr<tbox::Timer> RefineSchedule::t_shear;
std::shared_ptr<tbox::Timer> RefineSchedule::t_get_global_box_count;
std::shared_ptr<tbox::Timer> RefineSchedule::t_coarse_shear;
std::shared_ptr<tbox::Timer> RefineSchedule::t_setup_coarse_interp_box_level;
std::shared_ptr<tbox::Timer> RefineSchedule::t_bridge_coarse_interp_hiercoarse;
std::shared_ptr<tbox::Timer> RefineSchedule::t_bridge_dst_hiercoarse;
std::shared_ptr<tbox::Timer> RefineSchedule::t_invert_edges;
std::shared_ptr<tbox::Timer> RefineSchedule::t_construct_send_trans;
std::shared_ptr<tbox::Timer> RefineSchedule::t_construct_recv_trans;

tbox::StartupShutdownManager::Handler
RefineSchedule::s_initialize_finalize_handler(
   RefineSchedule::initializeCallback,
   0,
   0,
   RefineSchedule::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 **************************************************************************
 *
 * Create a refine schedule that copies data from the source level into
 * the destination level on the components represented by the refine
 * classes.  Ony data on the intersection of the two levels will be
 * copied. It is assumed that the index spaces of the source and
 * destination levels are "consistent"; i.e., they represent the same
 * grid resolution.  The levels do not have to be part of the same
 * AMR patch hierarchy, however.
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& dst_level_fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   const std::shared_ptr<RefineClasses>& refine_classes,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_number_refine_items(0),
   d_refine_items(0),
   d_dst_level(dst_level),
   d_src_level(src_level),
   d_refine_patch_strategy(patch_strategy),
   d_singularity_patch_strategy(dynamic_cast<SingularityPatchStrategy *>(patch_strategy)),
   d_transaction_factory(transaction_factory),
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_force_boundary_fill(false),
   d_num_periodic_directions(0),
   d_periodic_shift(dst_level->getDim()),
   d_coarse_priority_level_schedule(std::make_shared<tbox::Schedule>()),
   d_fine_priority_level_schedule(std::make_shared<tbox::Schedule>()),
   d_encon_level(std::make_shared<hier::PatchLevel>(dst_level->getDim())),
   d_dst_to_src(0),
   d_max_fill_boxes(0),
   d_dst_level_fill_pattern(dst_level_fill_pattern),
   d_top_refine_schedule(this),
   d_internal_allocated(false)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(refine_classes);
   TBOX_ASSERT(transaction_factory);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
#endif

   getFromInput();

   if (s_barrier_and_time) {
      t_refine_schedule->barrierAndStart();
   }

   if (d_dst_level->getGridGeometry()->getNumberOfBlockSingularities() > 0 &&
       !d_singularity_patch_strategy && d_refine_patch_strategy) {
      TBOX_ERROR("RefineSchedule: Schedules for meshes with singularities\n"
         << "requires a SingularityPatchStrategy implementation along\n"
         << "with the RefinePatchStrategy.  To do this,\n"
         << "inherit SinglarityPatchStrategy with the user class\n"
         << "that inherited RefinePatchStrategy and implement\n"
         << "the SingularityPatchStrategy pure virtual methods.");
   }

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_domain_is_one_box.resize(
      d_dst_level->getGridGeometry()->getNumberBlocks(), false);

   d_coarse_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule_fill");
   d_fine_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule_fill");

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */
   initializeDomainAndGhostInformation();

   hier::IntVector min_connector_width(getMinConnectorWidth());
   if (!d_dst_level_fill_pattern->fillingCoarseFineGhosts()) {
      min_connector_width = hier::IntVector::getZero(dst_level->getDim());
   }

   d_dst_to_src = &d_dst_level->findConnectorWithTranspose(*d_src_level,
         min_connector_width,
         hier::Connector::convertHeadWidthToBase(
            d_src_level->getBoxLevel()->getRefinementRatio(),
            d_dst_level->getBoxLevel()->getRefinementRatio(),
            min_connector_width),
         hier::CONNECTOR_IMPLICIT_CREATION_RULE,
         true);
   hier::Connector& src_to_dst = d_dst_to_src->getTranspose();

   TBOX_ASSERT(d_dst_to_src->getBase() == *d_dst_level->getBoxLevel());
   TBOX_ASSERT(src_to_dst.getHead() == *d_dst_level->getBoxLevel());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dst_level_fill_pattern->fillingCoarseFineGhosts()) {
      TBOX_ASSERT(d_dst_to_src->getConnectorWidth() >= d_max_scratch_gcw);
      TBOX_ASSERT(d_dst_to_src->getConnectorWidth() >= d_boundary_fill_ghost_width);
   } 
#endif

   if (s_extra_debug) {
      /*
       * This check may be redundant because
       * PersistentOverlapConnectors should already guarantee
       * completeness.
       */
      d_dst_to_src->assertOverlapCorrectness();
      src_to_dst.assertOverlapCorrectness();
   }

   /*
    * Create fill_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will attempt to fill.
    */

   std::shared_ptr<hier::BoxLevel> fill_box_level;
   std::shared_ptr<hier::Connector> dst_to_fill;
   hier::BoxNeighborhoodCollection dst_to_fill_on_src_proc;
   setDefaultFillBoxLevel(
      fill_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc);

   /*
    * Generation the communication transactions that will move data from
    * the source to the destination.  generateCommunicationSchedule will
    * initialize the "unused" objects with information about the parts of
    * fill_box_level that cannot be filled from the source level.
    * They are unused because this RefineSchedule constructor creates
    * schedules that do not do anything to fill the parts of the
    * destination that can't be filled directly from the source.
    */
   std::shared_ptr<hier::BoxLevel> unused_unfilled_box_level;
   std::shared_ptr<hier::Connector> unused_dst_to_unfilled;
   std::shared_ptr<hier::BoxLevel> unused_unfilled_encon_box_level;
   std::shared_ptr<hier::Connector> unused_encon_to_unfilled_encon;
   bool create_transactions = true;
   generateCommunicationSchedule(
      unused_unfilled_box_level,
      unused_dst_to_unfilled,
      unused_unfilled_encon_box_level,
      unused_encon_to_unfilled_encon,
      *dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_refinement,
      create_transactions);

   if (d_coarse_interp_level) {
      computeRefineOverlaps(d_refine_overlaps,
         d_dst_level,
         d_coarse_interp_level,
         d_dst_to_coarse_interp->getTranspose(),
         *d_coarse_interp_to_unfilled);
   }

   if (d_coarse_interp_encon_level) {
      computeRefineOverlaps(d_encon_refine_overlaps,
         d_encon_level,
         d_coarse_interp_encon_level,
         d_encon_to_coarse_interp_encon->getTranspose(),
         *d_coarse_interp_encon_to_unfilled_encon);
   }

   if (s_barrier_and_time) {
      t_refine_schedule->barrierAndStop();
   }
}

/*
 **************************************************************************
 *
 * Create a refine schedule that copies data from the source level into
 * the destination level on the components represented by the refine
 * classes.  If portions of the destination level remain unfilled, then
 * the algorithm recursively fills those unfilled portions from coarser
 * levels in the AMR hierarchy.  It is assumed that the index spaces of
 * the source and destination levels are "consistent"; i.e., they
 * represent the same grid resolution.  Also, the next coarser level
 * integer argument must be the number of level in the specified
 * hierarchy representing the next coarser level of mesh resolution to
 * the destination level.
 *
 * IMPORTANT NOTES: The source level may be NULL, in which case the
 * destination level will be filled only using data interpolated from
 * coarser levels in the AMR hierarchy.  The hierarchy may be NULL only
 * if the next coarser level is -1 (that is, there is no coarser level).
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   const std::shared_ptr<PatchLevelFillPattern>& dst_level_fill_pattern,
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   int next_coarser_ln,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::shared_ptr<RefineClasses>& refine_classes,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement):
   d_number_refine_items(0),
   d_refine_items(0),
   d_dst_level(dst_level),
   d_src_level(src_level),
   d_refine_patch_strategy(patch_strategy),
   d_singularity_patch_strategy(dynamic_cast<SingularityPatchStrategy *>(patch_strategy)),
   d_transaction_factory(transaction_factory),
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_force_boundary_fill(false),
   d_num_periodic_directions(0),
   d_periodic_shift(dst_level->getDim()),
   d_encon_level(std::make_shared<hier::PatchLevel>(dst_level->getDim())),
   d_dst_to_src(0),
   d_max_fill_boxes(0),
   d_dst_level_fill_pattern(dst_level_fill_pattern),
   d_top_refine_schedule(this),
   d_internal_allocated(false)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT((next_coarser_ln == -1) || hierarchy);
   TBOX_ASSERT(refine_classes);
   TBOX_ASSERT(transaction_factory);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (src_level) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
   }
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
   }
#endif

   getFromInput();

   if (s_barrier_and_time) {
      t_refine_schedule->barrierAndStart();
   }

   const tbox::Dimension& dim(dst_level->getDim());

   if (dst_level->getGridGeometry()->getNumberOfBlockSingularities() > 0 &&
       !d_singularity_patch_strategy && d_refine_patch_strategy) {
      TBOX_ERROR("RefineSchedule: Schedules for meshes with singularities\n"
         << "requires a SingularityPatchStrategy implementation along\n"
         << "with the RefinePatchStrategy.  To do this,\n"
         << "inherit SinglarityPatchStrategy with the user class\n"
         << "that inherited RefinePatchStrategy and implement\n"
         << "the SingularityPatchStrategy pure virtual methods.");
   }

   setRefineItems(refine_classes);
   initialCheckRefineClassItems();

   d_domain_is_one_box.resize(
      dst_level->getGridGeometry()->getNumberBlocks(), false);

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */
   initializeDomainAndGhostInformation();

   hier::IntVector min_connector_width(getMinConnectorWidth());

   if (d_src_level &&
       d_src_level->getRatioToLevelZero() != d_dst_level->getRatioToLevelZero()) {
      if (d_src_level->getRatioToLevelZero() >= d_dst_level->getRatioToLevelZero()) {
         const hier::IntVector src_dst_ratio =
            d_src_level->getRatioToLevelZero() / d_dst_level->getRatioToLevelZero();
         if (d_dst_level->getRatioToLevelZero() * src_dst_ratio !=
             d_src_level->getRatioToLevelZero()) {
            TBOX_ERROR("RefineSchedule::RefineSchedule error: source and destination\n"
               << "levels must be a simple refinement of one another.\n"
               << "src resolution: " << d_src_level->getRatioToLevelZero() << "\n"
               << "dst resolution: " << d_dst_level->getRatioToLevelZero());
         }
         min_connector_width *= src_dst_ratio;
      } else if (d_src_level->getRatioToLevelZero() <= d_dst_level->getRatioToLevelZero()) {
         TBOX_ERROR("RefineSchedule:RefineSchedule error: We are not currently\n"
            << "supporting RefineSchedules with the source level finer\n"
            << "than the destination level.");
      } else {
         TBOX_ERROR("RefineSchedule::RefineSchedule error: src level may not be\n"
            << "coarser than dst level in one direction and finer in another.\n"
            << "src resolution: " << d_src_level->getRatioToLevelZero() << "\n"
            << "dst resolution: " << d_dst_level->getRatioToLevelZero());
      }
   }

   if (next_coarser_ln >= 0) {
      RefineScheduleConnectorWidthRequestor rscwr;

      if (hierarchy->getNumberOfLevels() > next_coarser_ln + 1) {
         if (d_dst_level->getRatioToLevelZero() !=
             hierarchy->getPatchLevel(next_coarser_ln + 1)->getRatioToLevelZero()) {
            hier::IntVector expansion_ratio =
               hierarchy->getPatchLevel(next_coarser_ln+1)->getRatioToLevelZero() / d_dst_level->getRatioToLevelZero();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT( expansion_ratio * d_dst_level->getRatioToLevelZero() == hierarchy->getPatchLevel(next_coarser_ln+1)->getRatioToLevelZero() );
            // All values in expansion_ratio must be identical.
            TBOX_ASSERT( hier::IntVector(dim,expansion_ratio(0,0),expansion_ratio.getNumBlocks()) == expansion_ratio );
#endif
            rscwr.setGhostCellWidthFactor(expansion_ratio(0,0));
         }
      }
      rscwr.computeRequiredFineConnectorWidthsForRecursiveRefinement(
         d_fine_connector_widths,
         min_connector_width,
         d_max_stencil_width,
         *hierarchy,
         next_coarser_ln + 1);
   }

   std::shared_ptr<hier::Connector> dummy_connector(
      std::make_shared<hier::Connector>(dim));

   if (d_src_level) {

      hier::IntVector transpose_min_connector_width =
         hier::Connector::convertHeadWidthToBase(
            d_src_level->getBoxLevel()->getRefinementRatio(),
            dst_level->getBoxLevel()->getRefinementRatio(),
            min_connector_width);
      d_dst_to_src = &dst_level->findConnectorWithTranspose(*d_src_level,
            min_connector_width,
            transpose_min_connector_width,
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      TBOX_ASSERT(d_dst_to_src->getBase() == *dst_level->getBoxLevel());
      TBOX_ASSERT(d_dst_to_src->getTranspose().getHead() == *dst_level->getBoxLevel());
      TBOX_ASSERT(d_dst_to_src->getConnectorWidth() >= d_max_scratch_gcw);
      TBOX_ASSERT(d_dst_to_src->getConnectorWidth() >= d_boundary_fill_ghost_width);
   } else {
      dummy_connector->setTranspose(dummy_connector.get(), false);
      d_dst_to_src = dummy_connector.get();
   }

   /*
    * Create fill_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will fill.
    */

   std::shared_ptr<hier::BoxLevel> fill_box_level;
   std::shared_ptr<hier::Connector> dst_to_fill;
   hier::BoxNeighborhoodCollection dst_to_fill_on_src_proc;

   setDefaultFillBoxLevel(
      fill_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc);

   const bool skip_first_generate_schedule =
      !d_dst_level_fill_pattern->doesSourceLevelCommunicateToDestination();

   const hier::IntVector dummy_intvector(dim, -1);

   /*
    * finishScheduleConstruction sets up all transactions to communicate
    * data from source to destination, and sets up recursive schedules to
    * fill whatever cannot be filled by the source.
    */
   int errf = finishScheduleConstruction(
         next_coarser_ln,
         hierarchy,
         dummy_intvector,
         *dst_to_fill,
         dst_to_fill_on_src_proc,
         use_time_refinement,
         skip_first_generate_schedule);
   if (errf) {
      tbox::perr
      << "Internal error in RefineSchedule constructor..."
      << "\n dst_to_fill:\n" << dst_to_fill->format("\tDF->", 2)
      << "\n dst:\n" << d_dst_level->getBoxLevel()->format("\tD->", 2)
      << std::endl;
      TBOX_ERROR("Top RefineSchedule constructor aborting due to above error.");
      return;
   }

   /*
    * Compute the BoxOverlap objects that will be used to refine the
    * data from coarser levels onto the destination.
    */
   if (d_coarse_interp_schedule) {
      computeRefineOverlaps(d_refine_overlaps,
         d_dst_level,
         d_coarse_interp_level,
         d_dst_to_coarse_interp->getTranspose(),
         *d_coarse_interp_to_unfilled);
   }

   if (d_coarse_interp_encon_schedule) {
      computeRefineOverlaps(d_encon_refine_overlaps,
         d_encon_level,
         d_coarse_interp_encon_level,
         d_encon_to_coarse_interp_encon->getTranspose(),
         *d_coarse_interp_encon_to_unfilled_encon);
   }

   if (s_barrier_and_time) {
      t_refine_schedule->barrierAndStop();
   }
}

/*
 **************************************************************************
 *
 * This private constructor is used to create internal schedules that
 * fill internal levels that are used as coarse levels in refinement
 * operations.
 *
 **************************************************************************
 */

RefineSchedule::RefineSchedule(
   int& errf,
   const std::shared_ptr<hier::PatchLevel>& dst_level,
   const std::shared_ptr<hier::PatchLevel>& src_level,
   int next_coarser_ln,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::Connector& dst_to_src,
   const hier::IntVector& src_growth_to_nest_dst,
   const std::shared_ptr<RefineClasses>& refine_classes,
   const std::shared_ptr<RefineTransactionFactory>& transaction_factory,
   RefinePatchStrategy* patch_strategy,
   const RefineSchedule* top_refine_schedule):
   d_number_refine_items(0),
   d_refine_items(0),
   d_dst_level(dst_level),
   d_src_level(src_level),
   d_refine_patch_strategy(patch_strategy),
   d_singularity_patch_strategy(dynamic_cast<SingularityPatchStrategy *>(patch_strategy)),
   d_transaction_factory(transaction_factory),
   d_max_stencil_width(dst_level->getDim()),
   d_max_scratch_gcw(dst_level->getDim()),
   d_boundary_fill_ghost_width(dst_level->getDim()),
   d_force_boundary_fill(false),
   d_domain_is_one_box(dst_level->getGridGeometry()->getNumberBlocks(), false),
   d_num_periodic_directions(0),
   d_periodic_shift(dst_level->getDim()),
   d_encon_level(std::make_shared<hier::PatchLevel>(dst_level->getDim())),
   d_dst_to_src(&dst_to_src),
   d_max_fill_boxes(0),
   d_dst_level_fill_pattern(std::make_shared<PatchLevelFullFillPattern>()),
   d_top_refine_schedule(top_refine_schedule),
   d_internal_allocated(false)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT((next_coarser_ln == -1) || hierarchy);
   TBOX_ASSERT(dst_to_src.hasTranspose());
   TBOX_ASSERT(refine_classes);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
   if (hierarchy) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
   }
#endif

   // Don't time this constructor because it's recursive.

   getFromInput();

   /*
    * Initial values; some will change in setup operations.
    * Note that we do not check refine items here, since this
    * constructor is private and called recursively (i.e., the
    * items have been checked already).
    */

   setRefineItems(refine_classes);

   /*
    * Initialize destination level, ghost cell widths,
    * and domain information data members.
    */
   initializeDomainAndGhostInformation();

   hier::Connector& src_to_dst = d_dst_to_src->getTranspose();

   TBOX_ASSERT(d_dst_to_src->getBase() == *d_dst_level->getBoxLevel());
   TBOX_ASSERT(src_to_dst.getHead() == *d_dst_level->getBoxLevel());

   if (s_extra_debug) {
      src_to_dst.assertOverlapCorrectness(false, true, true);
      d_dst_to_src->assertOverlapCorrectness(false, true, true);
   }

   /*
    * Create fill_box_level, representing all parts of the
    * destination level, including ghost regions if desired, that this
    * schedule will fill.  Here, the destination is always a coarse interpolation
    * level constructed by coarsening another RefineSchedule's unfilled
    * boxes.  As the destination will be used as a coarse level in a
    * refinement operation, the fill_box_level will be the boxes
    * of the destination level grown by the maximum interplation stencil
    * width.
    */

   std::shared_ptr<hier::BoxLevel> fill_box_level;
   std::shared_ptr<hier::Connector> dst_to_fill;
   hier::BoxNeighborhoodCollection dst_to_fill_on_src_proc;
   setDefaultFillBoxLevel(
      fill_box_level,
      dst_to_fill,
      dst_to_fill_on_src_proc);

   bool use_time_refinement = true;

   /*
    * finishScheduleConstruction sets up all transactions to communicate
    * data from source to destination, and sets up recursive schedules to
    * fill whatever cannot be filled by the source.
    */

   errf = finishScheduleConstruction(
         next_coarser_ln,
         hierarchy,
         src_growth_to_nest_dst,
         *dst_to_fill,
         dst_to_fill_on_src_proc,
         use_time_refinement);
   if (errf) {
      tbox::perr
      << "Internal error in private RefineSchedule constructor..."
      << "\n next_coarser_ln: " << next_coarser_ln
      << "\n dst_to_fill:\n" << dst_to_fill->format("\tDF->", 2)
      << std::endl;
      return;
   }

   /*
    * Compute the BoxOverlap objects that will be used to refine the
    * data from coarser levels onto the destination.
    */

   if (d_coarse_interp_schedule) {
      computeRefineOverlaps(d_refine_overlaps,
         d_dst_level,
         d_coarse_interp_level,
         d_dst_to_coarse_interp->getTranspose(),
         *d_coarse_interp_to_unfilled);
   }

   if (d_coarse_interp_encon_schedule) {
      computeRefineOverlaps(d_encon_refine_overlaps,
         d_encon_level,
         d_coarse_interp_encon_level,
         d_encon_to_coarse_interp_encon->getTranspose(),
         *d_coarse_interp_encon_to_unfilled_encon);
   }

}

/*
 **************************************************************************
 *
 * The destructor for the refine schedule class implicitly deallocates
 * all of the data associated with the communication schedule.
 *
 **************************************************************************
 */

RefineSchedule::~RefineSchedule()
{
   clearRefineItems();
   delete[] d_refine_items;
   if (d_internal_allocated) {
      deallocateInternalData();
   }
}

/*
 *************************************************************************
 *
 * Read static member data from input database once.
 *
 ************************************************************************
 */
void
RefineSchedule::getFromInput()
{
   if (!s_read_static_input) {
      s_read_static_input = true;
      std::shared_ptr<tbox::Database> idb(
         tbox::InputManager::getInputDatabase());
      if (idb && idb->isDatabase("RefineSchedule")) {
         std::shared_ptr<tbox::Database> rsdb(
            idb->getDatabase("RefineSchedule"));
         s_extra_debug = rsdb->getBoolWithDefault("DEV_extra_debug", false);
         s_barrier_and_time =
            rsdb->getBoolWithDefault("DEV_barrier_and_time", false);
      }
   }
}

/*
 *************************************************************************
 *
 * Reset schedule with new set of refine items.
 *
 ************************************************************************
 */

void
RefineSchedule::reset(
   const std::shared_ptr<RefineClasses>& refine_classes)
{
   TBOX_ASSERT(refine_classes);

   if (d_internal_allocated) {
      deallocateInternalData();
   }

   setRefineItems(refine_classes);
   if (d_coarse_interp_schedule) {
      d_coarse_interp_schedule->reset(refine_classes);
   }
   if (d_coarse_interp_encon_schedule) {
      d_coarse_interp_encon_schedule->reset(refine_classes);
   }
}

/*
 ************************************************************************
 * Construct transactions for schedule and set up recursive schedules if
 * needed.
 *
 * Generate communication schedules to transfer data from src to
 * fillboxes associated with dst boxes.  What parts cannot be filled
 * from the src becomes the "unfilled" boxes.  If no source, all fill
 * boxes become "unfilled" boxes.  We also construct unfilled boxes at
 * enhanced connectivity block boundaries.
 *
 * If there are any unfilled boxes, we coarsen them to create a
 * coarse interpolation level and set up a recursive schedule for filling the
 * coarse interpolation level.  The idea is to interpolate data from the
 * coarse interpolation level to fill the unfilled boxes.
 ************************************************************************
 */

int
RefineSchedule::finishScheduleConstruction(
   int next_coarser_ln,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& src_growth_to_nest_dst,
   const hier::Connector& dst_to_fill,
   const hier::BoxNeighborhoodCollection& dst_to_fill_on_src_proc,
   bool use_time_interpolation,
   bool skip_generate_schedule)
{
   if (s_barrier_and_time) {
      t_finish_sched_const->barrierAndStart();
   }
   TBOX_ASSERT(d_dst_to_src);
   TBOX_ASSERT(d_dst_to_src->hasTranspose());
   TBOX_ASSERT((next_coarser_ln == -1) || hierarchy);

   // Get data that will be used below.

   const tbox::Dimension& dim(hierarchy->getDim());

   hier::BoxLevelConnectorUtils edge_utils;
   hier::OverlapConnectorAlgorithm oca;
   oca.setTimerPrefix("xfer::RefineSchedule_build");

   if (d_src_level) {
      // Should never have a source without connection from destination.
      TBOX_ASSERT(d_dst_to_src->isFinalized());
   }

   d_coarse_priority_level_schedule.reset(new tbox::Schedule());
   d_fine_priority_level_schedule.reset(new tbox::Schedule());

   d_coarse_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule_fill");
   d_fine_priority_level_schedule->setTimerPrefix("xfer::RefineSchedule_fill");

   /*
    * Generate the schedule for filling the boxes in dst_to_fill.
    * Any portions of the fill boxes that cannot be filled from
    * the source is placed in d_unfilled_box_level.
    *
    * If the source is not given or skip_generate_schedule==true,
    * the schedule generation degenates to turning all the fill boxes
    * into unfilled boxes.
    */

   std::shared_ptr<hier::Connector> dst_to_unfilled;
   const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry(
      d_dst_level->getGridGeometry());
   std::shared_ptr<hier::Connector> encon_to_unfilled_encon;

   bool create_transactions = true;
   if (!d_src_level || skip_generate_schedule) {
      create_transactions = false;
   }

   generateCommunicationSchedule(
      d_unfilled_box_level,
      dst_to_unfilled,
      d_unfilled_encon_box_level,
      encon_to_unfilled_encon,
      dst_to_fill,
      dst_to_fill_on_src_proc,
      use_time_interpolation,
      create_transactions);

   /*
    * d_unfilled_box_level may include ghost cells that lie
    * outside the physical domain.  These parts must be removed if
    * they are not at periodic boundaries.  They will be filled
    * through a user call-back method.
    */

   shearUnfilledBoxesOutsideNonperiodicBoundaries(
      *d_unfilled_box_level,
      *dst_to_unfilled,
      hierarchy);

   t_get_global_box_count->barrierAndStart();

   const bool need_to_fill =
      (d_unfilled_box_level->getGlobalNumberOfBoxes() > 0);

   const bool need_to_fill_encon =
      grid_geometry->hasEnhancedConnectivity() &&
      (d_unfilled_encon_box_level->getGlobalNumberOfBoxes() > 0);

   t_get_global_box_count->stop();

   /*
    * If there remain boxes to be filled from coarser levels, then set
    * up data for recursive schedule generation:
    *
    * 1. Generate a coarse interpolation BoxLevel
    * (coarse_interp_box_level) by coarsening the unfilled boxes.
    *
    * 2. Connect coarse_interp_box_level to the next coarser level on
    * the hierarchy.
    *
    * 3: Construct the coarse interpolation PatchLevel (d_coarse_interp_level)
    * and construct d_coarse_interp_schedule to fill d_coarse_interp_level.
    * The coarser level on the hierarchy will be the source for filling
    * d_coarse_interp_level, which is why we need step 2..
    *
    * The idea is that once d_coarse_interp_level is filled, we can refine its
    * data to fill the current unfilled boxes.
    */

   if (need_to_fill) {

      t_finish_sched_const_recurse->start();

      makeNodeCenteredUnfilledBoxLevel(*d_unfilled_box_level,
         *dst_to_unfilled);

      /*
       * If there are no coarser levels in the hierarchy or the
       * hierarchy is null, then throw an error.  Something is messed
       * up someplace and code execution cannot proceed.
       */
      if (next_coarser_ln < 0) {
         tbox::perr
         << "Internal error in RefineSchedule::finishScheduleConstruction..."
         << "\n In finishScheduleConstruction() -- "
         << "\n No coarser levels...will not fill from coarser."
         << "\n next_coarser_ln: " << next_coarser_ln
         << "\n src_growth_to_nest_dst: " << src_growth_to_nest_dst
         << "\n dst_to_unfilled:\n" << dst_to_unfilled->format("\tDU->", 2)
         << "\n d_unfilled_box_level:\n" << d_unfilled_box_level->format("\tUF->", 2)
         << std::endl;
         return 1;
      } else {
         if (!hierarchy) {
            tbox::perr
            << "Internal RefineSchedule error..."
            << "\n In finishScheduleConstruction() -- "
            << "\n Need to fill from coarser hierarchy level and \n"
            << "hierarchy is unavailable." << std::endl;
            return 2;
         }
      }

      /*
       * hiercoarse is the coarse level on the hierarchy.  It is to be
       * differentiated from the coarse interpolation (coarse_interp) level,
       * which is at the same resolution and level number but is not on the
       * hierarchy.
       */
      const std::shared_ptr<hier::PatchLevel> hiercoarse_level(
         hierarchy->getPatchLevel(next_coarser_ln));

      const hier::BoxLevel& hiercoarse_box_level(
         *hiercoarse_level->getBoxLevel());

      /*
       * Ratio to the next coarser level in the hierarchy.
       */
      const hier::IntVector dst_hiercoarse_ratio(
         d_dst_level->getRatioToLevelZero()
         / hiercoarse_level->getRatioToLevelZero());

      /*
       * Set up the coarse interpolation BoxLevel and also set up
       * d_dst_to_coarse_interp, its transpose and
       * d_coarse_interp_to_unfilled.  These
       * Connectors are easily generated using dst_to_unfilled.
       */

      std::shared_ptr<hier::BoxLevel> coarse_interp_box_level;
      setupCoarseInterpBoxLevel(
         coarse_interp_box_level,
         d_dst_to_coarse_interp,
         d_coarse_interp_to_unfilled,
         hiercoarse_box_level,
         *dst_to_unfilled);

      /*
       * Create the coarse interpolation PatchLevel and connect its
       * BoxLevel (the next recursion's dst) to the hiercoarse
       * BoxLevel (the next recursion's src).
       */

      std::shared_ptr<hier::Connector> coarse_interp_to_hiercoarse;

      createCoarseInterpPatchLevel(
         d_coarse_interp_level,
         coarse_interp_box_level,
         coarse_interp_to_hiercoarse,
         next_coarser_ln,
         hierarchy,
         *d_dst_to_src,
         *d_dst_to_coarse_interp,
         d_dst_level);

      /*
       * Compute how much hiercoarse would have to grow to nest coarse_interp, a
       * required parameter in the private constructor.
       *
       * If dst is a coarse interpolation level (generated by RefineSchedule),
       * we have the info to compute the growth.  If not, we make some
       * assumptions about where dst came from in order to determine
       * how its fill boxes nest in hiercoarse.
       */
      hier::IntVector hiercoarse_growth_to_nest_coarse_interp(
         hier::IntVector::getZero(dim));
      const bool dst_is_coarse_interp_level = this != d_top_refine_schedule;
      if (dst_is_coarse_interp_level) {
         /*
          * Assume that src barely nests in hiercoarse.  (In most
          * places, it nests by a margin equal to the nesting buffer,
          * but we don't count on that because the nesting buffer is
          * zero at physical boundaries.)  To nest dst, hiercoarse
          * has to grow as much as the src does, plus the ghost width
          * of the fill.
          *
          * REMARK: We may in fact be able to count on the nesting
          * buffer because extending boxes to physical boundaries do
          * not create any extra relationships.  However, we don't
          * currently have access to the size of the nesting buffer.
          */
         hiercoarse_growth_to_nest_coarse_interp =
            src_growth_to_nest_dst + dst_to_fill.getConnectorWidth();
      } else {
         /*
          * dst may be:
          * 1. The hierarchy level just finer than level number next_coarser_ln.
          * 2. A level that nests in level number next_coarser_ln:
          *    a. A new level generated by GriddingAlgorithm.
          *    b. The hierarchy level just finer than level number next_coarser_ln,
          *       coarsened for Richardson extrapolation.
          * In any case, dst should nest in hiercoarse.  Furthermore, it does
          * not grow when coarsened into the hiercoarse index space.
          * To nest dst and its fill boxes, hiercoarse just has to grow by
          * the ghost width of the fill.
          */
         hiercoarse_growth_to_nest_coarse_interp =
            dst_to_fill.getConnectorWidth();
      }
      hiercoarse_growth_to_nest_coarse_interp.ceilingDivide(
         dst_hiercoarse_ratio);

      t_finish_sched_const_recurse->stop();

      /*
       * We now have all the data for building the coarse interpolation
       * schedule using the private constructor.
       *
       * We need to make sure that the coarse schedule uses
       * BoxGeometryVariableFillPattern, so that it fills all needed
       * parts of d_coarse_interp_level
       */
      std::shared_ptr<BoxGeometryVariableFillPattern> bg_fill_pattern(
         std::make_shared<BoxGeometryVariableFillPattern>());

      std::shared_ptr<RefineClasses> coarse_schedule_refine_classes(
         std::make_shared<RefineClasses>());

      const int num_refine_items =
         d_refine_classes->getNumberOfRefineItems();

      for (int nd = 0; nd < num_refine_items; ++nd) {
         RefineClasses::Data item = d_refine_classes->getRefineItem(nd);
         item.d_var_fill_pattern = bg_fill_pattern;
         coarse_schedule_refine_classes->insertEquivalenceClassItem(item);
      }

      if (t_finish_sched_const->isRunning()) {
         t_finish_sched_const->stop();
      }

      int errf;
      d_coarse_interp_schedule.reset(new RefineSchedule(errf,
            d_coarse_interp_level,
            hiercoarse_level,
            next_coarser_ln - 1,
            hierarchy,
            *coarse_interp_to_hiercoarse,
            hiercoarse_growth_to_nest_coarse_interp,
            coarse_schedule_refine_classes,
            d_transaction_factory,
            d_refine_patch_strategy,
            d_top_refine_schedule));
      if (errf) {
         tbox::perr
         << "In finishScheduleConstruction after failure to generate d_coarse_interp_schedule:"
         << "\n next_coarser_ln: " << next_coarser_ln
         << "\n src_growth_to_nest_dst: " << src_growth_to_nest_dst
         << "\n coarse_interp_level:\n" << d_coarse_interp_level->getBoxLevel()->format("\tCI->", 2)
         << "\n coarse_interp_to_hiercoarse:\n" << coarse_interp_to_hiercoarse->format("\tCHC->", 2)
         << std::endl;
         return errf;
      }

   } else {
      if (t_finish_sched_const->isRunning()) {
         t_finish_sched_const->stop();
      }
   }

   if (need_to_fill_encon) {

      /*
       * Create schedule to fill unfilled boxes at enhanced connectivity
       */

      const std::shared_ptr<hier::PatchLevel> hiercoarse_level(
         hierarchy->getPatchLevel(next_coarser_ln));

      createEnconFillSchedule(
         hierarchy,
         hiercoarse_level,
         src_growth_to_nest_dst,
         *encon_to_unfilled_encon);
   }

   return 0;
}

/*
 ***********************************************************************
 * Create schedule for filling unfilled boxes at enhanced connectivity
 *
 * d_coarse_interp_encon_level is created by coarsening d_unfilled_encon_level.
 * d_coarse_interp_encon_schedule is created to communicate data from the
 * hierarchy to fill d_coarse_interp_encon_level.
 ***********************************************************************
 */
void
RefineSchedule::createEnconFillSchedule(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::shared_ptr<hier::PatchLevel>& hiercoarse_level,
   const hier::IntVector& src_growth_to_nest_dst,
   const hier::Connector& encon_to_unfilled_encon)
{
   const hier::BoxLevel& hiercoarse_box_level(
      *hiercoarse_level->getBoxLevel());

   const int next_coarser_ln = hiercoarse_level->getLevelNumber();

   /*
    * Ratio to the next coarser level in the hierarchy.
    */
   const hier::IntVector dst_hiercoarse_ratio(
      d_dst_level->getRatioToLevelZero()
      / hiercoarse_level->getRatioToLevelZero());

   std::shared_ptr<hier::BoxLevel> coarse_interp_encon_box_level;

   setupCoarseInterpBoxLevel(
      coarse_interp_encon_box_level,
      d_encon_to_coarse_interp_encon,
      d_coarse_interp_encon_to_unfilled_encon,
      hiercoarse_box_level,
      encon_to_unfilled_encon);

   std::shared_ptr<hier::Connector> coarse_interp_encon_to_hiercoarse;

   createCoarseInterpPatchLevel(
      d_coarse_interp_encon_level,
      coarse_interp_encon_box_level,
      coarse_interp_encon_to_hiercoarse,
      next_coarser_ln,
      hierarchy,
      *d_encon_to_src,
      *d_encon_to_coarse_interp_encon,
      d_encon_level);

   /*
    * Compute this nesting value the same as for coarse_interp
    */
   hier::IntVector hiercoarse_growth_to_nest_coarse_interp_encon(
      hier::IntVector::getZero(hiercoarse_level->getDim()));
   const bool dst_is_coarse_interp_level = this != d_top_refine_schedule;
   if (dst_is_coarse_interp_level) {
      hiercoarse_growth_to_nest_coarse_interp_encon =
         src_growth_to_nest_dst + encon_to_unfilled_encon.getConnectorWidth();
      hiercoarse_growth_to_nest_coarse_interp_encon.ceilingDivide(dst_hiercoarse_ratio);
   } else {
      hiercoarse_growth_to_nest_coarse_interp_encon =
         encon_to_unfilled_encon.getConnectorWidth();
      hiercoarse_growth_to_nest_coarse_interp_encon.ceilingDivide(dst_hiercoarse_ratio);
   }

   /*
    * We need to make sure that the coarse schedule uses
    * BoxGeometryVariableFillPattern, so that it fills all needed parts of
    * d_coarse_interp_encon_level
    */
   std::shared_ptr<BoxGeometryVariableFillPattern> bg_fill_pattern(
      std::make_shared<BoxGeometryVariableFillPattern>());

   std::shared_ptr<RefineClasses> coarse_schedule_refine_classes(
      std::make_shared<RefineClasses>());

   const int num_refine_items =
      d_refine_classes->getNumberOfRefineItems();

   for (int nd = 0; nd < num_refine_items; ++nd) {
      RefineClasses::Data item = d_refine_classes->getRefineItem(nd);
      item.d_var_fill_pattern = bg_fill_pattern;
      coarse_schedule_refine_classes->insertEquivalenceClassItem(item);
   }

   /*
    * Schedule to fill d_coarse_interp_encon_level
    */
   int errf;
   d_coarse_interp_encon_schedule.reset(new RefineSchedule(
         errf,
         d_coarse_interp_encon_level,
         hiercoarse_level,
         next_coarser_ln - 1,
         hierarchy,
         *coarse_interp_encon_to_hiercoarse,
         hiercoarse_growth_to_nest_coarse_interp_encon,
         coarse_schedule_refine_classes,
         d_transaction_factory,
         d_refine_patch_strategy,
         d_top_refine_schedule));
   if (errf) {
      TBOX_ERROR("RefineSchedule constructor aborting due to above errors.");
   }

}

/*
 **************************************************************************
 * Shear off parts of unfilled boxes that lie outside non-periodic
 * domain boundaries.
 *
 * Update the overlap Connector dst_to_unfilled.
 *
 * If domain is fully periodic, this is a no-op.
 **************************************************************************
 */

void
RefineSchedule::shearUnfilledBoxesOutsideNonperiodicBoundaries(
   hier::BoxLevel& unfilled,
   hier::Connector& dst_to_unfilled,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   const tbox::Dimension& dim(unfilled.getDim());

   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   if (fully_periodic) {
      /*
       * Bypass shearing for fully_periodic domains, where it would be a
       * no-op anyway.
       */
      return;
   }

   t_shear->start();

   const hier::Connector& unfilled_to_periodic_domain(
      d_unfilled_box_level->findConnector(
         hierarchy->getDomainBoxLevel(),
         dst_to_unfilled.getConnectorWidth(),
         hier::CONNECTOR_CREATE));

   std::shared_ptr<hier::MappingConnector> unfilled_to_sheared;
   std::shared_ptr<hier::BoxLevel> sheared_box_level;

   hier::BoxLevelConnectorUtils edge_utils;
   edge_utils.setTimerPrefix("xfer::RefineSchedule_build");
   edge_utils.computeInternalParts(
      sheared_box_level,
      unfilled_to_sheared,
      unfilled_to_periodic_domain,
      hier::IntVector::getZero(dim));

   hier::MappingConnectorAlgorithm mca;
   mca.setTimerPrefix("xfer::RefineSchedule_build");
   mca.setBarrierBeforeCommunication(false); // Next modify needs not communicate.
   mca.modify(dst_to_unfilled,
      *unfilled_to_sheared,
      d_unfilled_box_level.get());
   dst_to_unfilled.eraseEmptyNeighborSets();

   t_shear->stop();
}

/*
 **************************************************************************
 * Make the node-centered unfilled box level.
 **************************************************************************
 */

void
RefineSchedule::makeNodeCenteredUnfilledBoxLevel(
   const hier::BoxLevel& unfilled_box_level,
   const hier::Connector& dst_to_unfilled)
{
   const tbox::Dimension& dim = unfilled_box_level.getDim();

   d_unfilled_node_box_level.reset(new hier::BoxLevel(
         unfilled_box_level.getRefinementRatio(),
         unfilled_box_level.getGridGeometry(),
         unfilled_box_level.getMPI()));

   d_unfilled_to_unfilled_node.reset(new hier::Connector(
         unfilled_box_level,
         *d_unfilled_node_box_level,
         hier::IntVector::getZero(unfilled_box_level.getDim())));

   const hier::BoxLevel& dst_box_level = dst_to_unfilled.getBase();

   hier::Connector const* dst_to_src;
   if (d_src_level) {
      dst_to_src = d_dst_to_src;
   } else {
      dst_to_src = 0;
   }

   hier::LocalId last_unfilled_node_local_id(-1);
   const std::shared_ptr<const hier::BaseGridGeometry>& grid_geometry(
      dst_to_unfilled.getBase().getGridGeometry());

   for (hier::Connector::ConstNeighborhoodIterator cf = dst_to_unfilled.begin();
        cf != dst_to_unfilled.end(); ++cf) {

      const hier::BoxId& dst_box_id(*cf);
      const hier::Box& dst_box = *dst_box_level.getBox(dst_box_id);
      const hier::BlockId& dst_block = dst_box.getBlockId();
      hier::Box dst_node_box(dst_box);
      dst_node_box.setUpper(
         dst_node_box.upper() + hier::IntVector::getOne(dim));

      for (hier::Connector::ConstNeighborIterator ni = dst_to_unfilled.begin(cf);
           ni != dst_to_unfilled.end(cf); ++ni) {

         const hier::Box& unfilled_box = *ni;

         hier::BoxContainer unfilled_node_boxes(unfilled_box);
         unfilled_node_boxes.begin()->setUpper(
            unfilled_node_boxes.begin()->upper() + hier::IntVector::getOne(dim));

         if (d_src_level && dst_to_src->hasNeighborSet(dst_box_id)) {

            hier::Connector::ConstNeighborhoodIterator dst_src_itr =
               dst_to_src->findLocal(dst_box_id);

            for (hier::Connector::ConstNeighborIterator ds =
                    dst_to_src->begin(dst_src_itr);
                 ds != dst_to_src->end(dst_src_itr); ++ds) {

               const hier::Box& src_box = *ds;
               const hier::BlockId& src_block = src_box.getBlockId();
               hier::Box src_node_box(src_box);
               if (src_block == dst_block) {
                  src_node_box.setUpper(
                     src_node_box.upper() + hier::IntVector::getOne(dim));
               } else if (grid_geometry->areNeighbors(src_block, dst_block)) {
                  grid_geometry->transformBox(
                     src_node_box,
                     d_dst_level->getLevelNumber(),
                     dst_block,
                     src_block);
                  src_node_box.setUpper(
                     src_node_box.upper() + hier::IntVector::getOne(dim));
               }
               if (!(src_node_box * dst_node_box).empty()) {
                  unfilled_node_boxes.removeIntersections(src_node_box);
               }
            }
         }

         unfilled_node_boxes.coalesce();

         hier::Connector::NeighborhoodIterator unfilled_itr =
            d_unfilled_to_unfilled_node->
            makeEmptyLocalNeighborhood(unfilled_box.getBoxId());

         for (hier::BoxContainer::iterator ui = unfilled_node_boxes.begin();
              ui != unfilled_node_boxes.end(); ++ui) {

            hier::Box unfilled_nodal(*ui,
                                     ++last_unfilled_node_local_id,
                                     dst_box.getOwnerRank());
            TBOX_ASSERT(unfilled_nodal.getBlockId() == dst_block);

            d_unfilled_node_box_level->addBoxWithoutUpdate(unfilled_nodal);

            d_unfilled_to_unfilled_node->insertLocalNeighbor(
               unfilled_nodal,
               unfilled_itr);
         }

      }
   }

   d_unfilled_node_box_level->finalize();

}

/*
 **************************************************************************
 * Set up the coarse interpolation BoxLevel.  The coarse_interp_box_level
 * represents parts of the fill boxes that cannot be filled by
 * d_src_level (d_unfilled_box_level).  It will be filled by
 * appealing to coarser levels in the hierarchy.
 *
 * - Start with the d_unfilled_box_level.
 * - Coarsen unfilled boxes
 * - Shear off parts outside non-periodic boundaries or extend to boundary,
 *   if needed for conforming to gridding restrictions at boundary.
 *
 * Build Connectors d_dst_to_coarse_interp, its transpose and
 * d_coarse_interp_to_unfilled.  The ghost data to be filled on the coarse
 * interpolation level will be the max stencil width.
 *
 * We set d_dst_to_coarse_interp's width big enough so each dst Box, grown by
 * this width, nests its potential coarse interpolation Boxes.  Note that
 * d_dst_to_coarse_interp is incomplete because each dst Box only has edges to
 * coarse interpolation Boxes it generated.
 **************************************************************************
 */
void
RefineSchedule::setupCoarseInterpBoxLevel(
   std::shared_ptr<hier::BoxLevel>& coarse_interp_box_level,
   std::shared_ptr<hier::Connector>& dst_to_coarse_interp,
   std::shared_ptr<hier::Connector>& coarse_interp_to_unfilled,
   const hier::BoxLevel& hiercoarse_box_level,
   const hier::Connector& dst_to_unfilled)
{
   t_setup_coarse_interp_box_level->start();

   const tbox::Dimension& dim(dst_to_unfilled.getBase().getDim());

   const hier::IntVector dst_hiercoarse_ratio(
      d_dst_level->getRatioToLevelZero()
      / hiercoarse_box_level.getRefinementRatio());

   const bool fully_periodic = d_num_periodic_directions == dim.getValue();

   const std::shared_ptr<const hier::BaseGridGeometry>& grid_geometry(
      dst_to_unfilled.getBase().getGridGeometry());

   const size_t nblocks = grid_geometry->getNumberBlocks();

   hier::IntVector big_grow_vector(dim, 0);
   if (d_num_periodic_directions > 0) {
      for (int dir = 0; dir < dim.getValue(); ++dir) {
         if (d_periodic_shift(dir)) {
            big_grow_vector(dir) = BIG_GHOST_CELL_WIDTH;
         }
      }
   }
   hier::IntVector multi_big_grow(big_grow_vector, nblocks);

   std::vector<hier::BoxContainer> coarser_physical_domain(nblocks);
   std::vector<hier::BoxContainer> coarser_shear_domain(nblocks);
   std::vector<bool> do_coarse_shearing(nblocks);
   for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
      grid_geometry->computePhysicalDomain(
         coarser_physical_domain[b],
         hiercoarse_box_level.getRefinementRatio(),
         hier::BlockId(b));

      do_coarse_shearing[b] = (!fully_periodic && !d_domain_is_one_box[b]);

      if (do_coarse_shearing[b]) {
         t_coarse_shear->start();

         coarser_shear_domain[b] = coarser_physical_domain[b];

         if (d_num_periodic_directions > 0) {
            coarser_shear_domain[b].grow(multi_big_grow);
         }

         coarser_shear_domain[b].unorder();
         coarser_shear_domain[b].simplify();

         t_coarse_shear->stop();
      }
   }

   coarse_interp_box_level.reset(new hier::BoxLevel(
         hiercoarse_box_level.getRefinementRatio(),
         grid_geometry,
         hiercoarse_box_level.getMPI()));

   std::shared_ptr<hier::BoxLevel> nbr_blk_fill_box_level(
      new hier::BoxLevel(
         d_dst_level->getRatioToLevelZero(),
         grid_geometry,
         hiercoarse_box_level.getMPI())); 

   /*
    * Width of dst-->coarse_interp is
    *
    * - width of dst-->fill, but rounded up so it extends
    *   the growth of coarse_interp caused by coarsening unfilled.
    *
    * - extended by the stencil width, where coarse_interp has its ghost data.
    *
    * This width states that each dst box sees all of its
    * coarse interpolation boxes, including the ghost cells in the
    * coarse interpolation boxes.
    */
   const hier::IntVector dst_to_coarse_interp_width =
      (hier::IntVector::ceilingDivide(dst_to_unfilled.getConnectorWidth(),
          dst_hiercoarse_ratio) + d_max_stencil_width)
      * dst_hiercoarse_ratio;

   dst_to_coarse_interp.reset(new hier::Connector(dst_to_unfilled.getBase(),
         *coarse_interp_box_level,
         dst_to_coarse_interp_width));

   coarse_interp_to_unfilled.reset(new hier::Connector(
         *coarse_interp_box_level,
         dst_to_unfilled.getHead(),
         hier::IntVector::getZero(dim)));

   d_coarse_interp_to_nbr_fill.reset(new hier::Connector(
         *coarse_interp_box_level,
         *nbr_blk_fill_box_level,
         hier::IntVector::getZero(dim)));

   /*
    * This loop builds up coarse_interp_box_level.  It also builds up
    * the neighborhood sets of dst_to_coarse_interp and
    * coarse_interp_to_unfilled using simple associations in dst_to_unfilled.
    */
   for (hier::Connector::ConstNeighborhoodIterator ei = dst_to_unfilled.begin();
        ei != dst_to_unfilled.end(); ++ei) {

      const hier::BoxId& dst_box_id = *ei;

      hier::Connector::NeighborhoodIterator dst_base_box_itr =
         dst_to_coarse_interp->findLocal(dst_box_id);
      bool has_base_box = dst_base_box_itr != dst_to_coarse_interp->end();

      for (hier::Connector::ConstNeighborIterator ni = dst_to_unfilled.begin(ei);
           ni != dst_to_unfilled.end(ei); ++ni) {

         /*
          * Build the coarse interpolation box by coarsening the unfilled box.
          * The box is sheared at the physical boundary or extended to
          * the physical boundary to conform to normal gridding
          * restrictions.
          */
         const hier::Box& unfilled_box = *ni;
         const hier::BlockId& dst_blk_id = unfilled_box.getBlockId();
         const hier::BlockId::block_t& dst_blk = dst_blk_id.getBlockValue();

         hier::Box coarse_interp_box(unfilled_box);
         hier::BoxContainer coarse_interp_boxes; 
         if (nblocks == 1 || grid_geometry->hasIsotropicRatios()) {
            coarse_interp_boxes.pushBack(coarse_interp_box);
         } else {
            hier::BoxUtilities::growAndAdjustAcrossBlockBoundary(
               coarse_interp_boxes,
               coarse_interp_box,
               grid_geometry,
               d_dst_level->getBoxLevel()->getRefinementRatio(),
               hier::IntVector::getOne(dim),
               hier::IntVector::getZero(dim),
               false,
               false);
         }

         coarse_interp_boxes.coarsen(dst_hiercoarse_ratio);

         std::vector<hier::BoxContainer> sheared_coarse_interp_boxes(nblocks);

         if (nblocks > 1) {
            sheared_coarse_interp_boxes.resize(nblocks);
            for (hier::BoxContainer::iterator ci = coarse_interp_boxes.begin();
                 ci != coarse_interp_boxes.end(); ++ci) {
               const hier::BlockId& cblock_id = ci->getBlockId(); 
               sheared_coarse_interp_boxes[cblock_id.getBlockValue()].
                  pushBack(*ci);
            }
         } else {
            sheared_coarse_interp_boxes[0].spliceBack(coarse_interp_boxes);
         }

         bool sheared_boxes_exist = false;
         for (hier::BlockId::block_t blk = 0; blk < nblocks; ++blk) {
            if (do_coarse_shearing[blk] &&
                !sheared_coarse_interp_boxes[blk].empty() &&
                (d_dst_level->patchTouchesRegularBoundary(dst_box_id))) {
               sheared_coarse_interp_boxes[blk].intersectBoxes(
                  coarser_shear_domain[blk]);
               sheared_coarse_interp_boxes[blk].simplify();
            }

            if (!sheared_coarse_interp_boxes[blk].empty()) { 
               (void)hier::BoxUtilities::extendBoxesToDomainBoundary(
                  sheared_coarse_interp_boxes[blk],
                  coarser_physical_domain[blk],
                  d_max_stencil_width);
               sheared_boxes_exist = true;
            }
         }

         if (sheared_boxes_exist) {

            if (!has_base_box) {
               dst_base_box_itr =
                  dst_to_coarse_interp->makeEmptyLocalNeighborhood(dst_box_id);
               has_base_box = true;
            }

            for (hier::BlockId::block_t blk = dst_blk;
                 blk < nblocks + dst_blk; ++blk) {
               hier::BlockId::block_t cur_blk =
                  static_cast<hier::BlockId::block_t>(blk % nblocks);
               for (hier::BoxContainer::iterator bi =
                    sheared_coarse_interp_boxes[cur_blk].begin();
                    bi != sheared_coarse_interp_boxes[cur_blk].end(); ++bi) {
                  const hier::Box& coarse_interp_level_box =
                     *coarse_interp_box_level->addBox(*bi, bi->getBlockId());
                  dst_to_coarse_interp->insertLocalNeighbor(
                     coarse_interp_level_box,
                     dst_base_box_itr);

                  coarse_interp_to_unfilled->insertLocalNeighbor(
                     unfilled_box,
                     coarse_interp_level_box.getBoxId());

                  if (bi->getBlockId() != dst_blk_id &&
                      !grid_geometry->hasIsotropicRatios()) {
                     hier::Box fill_box(coarse_interp_level_box);
                     fill_box.refine(dst_hiercoarse_ratio);
                     nbr_blk_fill_box_level->addBoxWithoutUpdate(fill_box);
                     d_coarse_interp_to_nbr_fill->insertLocalNeighbor(
                        fill_box,
                        coarse_interp_level_box.getBoxId());
                  } 
               }
            }
         }
      }
   }

   nbr_blk_fill_box_level->finalize();
   if (nblocks > 1 && !grid_geometry->hasIsotropicRatios()) {
      d_nbr_blk_fill_level.reset(new hier::PatchLevel(
         nbr_blk_fill_box_level,
         d_dst_level->getGridGeometry(),
         d_dst_level->getPatchDescriptor()));
      d_nbr_blk_fill_level->setLevelNumber(d_dst_level->getLevelNumber());
      d_nbr_blk_fill_level->setNextCoarserHierarchyLevelNumber(
         d_dst_level->getLevelNumber()-1);

      d_nbr_blk_fill_level->getGridGeometry()->
         adjustMultiblockPatchLevelBoundaries(*d_nbr_blk_fill_level);
   } 

   /*
    * Get the transpose of dst_to_coarse_interp, which is simple to compute
    * because we know the edges are all local.
    */
   dst_to_coarse_interp->setTranspose(
      dst_to_coarse_interp->createLocalTranspose(),
      true);

   t_setup_coarse_interp_box_level->stop();
}

/*
 **************************************************************************
 * Create the coarse interpolation PatchLevel and build the Connectors
 * between the it and the coarser level on the hierarchy.
 *
 * The coarser level on the hierarchy at the same resolution as the
 * coarse interpolation level and will be used as the source for
 * filling the coarse interpolation level.
 **************************************************************************
 */
void
RefineSchedule::createCoarseInterpPatchLevel(
   std::shared_ptr<hier::PatchLevel>& coarse_interp_level,
   std::shared_ptr<hier::BoxLevel>& coarse_interp_box_level,
   std::shared_ptr<hier::Connector>& coarse_interp_to_hiercoarse,
   const int next_coarser_ln,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::Connector& dst_to_src,
   const hier::Connector& dst_to_coarse_interp,
   const std::shared_ptr<hier::PatchLevel>& dst_level)
{
   TBOX_ASSERT(dst_to_src.hasTranspose());
   TBOX_ASSERT(dst_to_coarse_interp.hasTranspose());

   const tbox::Dimension& dim(hierarchy->getDim());

   hier::OverlapConnectorAlgorithm oca;
   oca.setTimerPrefix("xfer::RefineSchedule_build");
   hier::BoxLevelConnectorUtils edge_utils;

   const std::shared_ptr<hier::PatchLevel> hiercoarse_level(
      hierarchy->getPatchLevel(next_coarser_ln));

   const hier::BoxLevel& hiercoarse_box_level(
      *hiercoarse_level->getBoxLevel());

   const hier::IntVector& zero_vec(hier::IntVector::getZero(dim));

   /*
    * To compute coarse_interp<==>hiercoarse, we will perform this bridge:
    * coarse_interp<==>dst<==>hiercoarse.
    *
    * First, we need dst<==>hiercoarse.
    *
    * If dst is a coarse interpolation level, we know that src is on the
    * hierarchy.  Get the cached src<==>hiercoarse and bridge
    * dst<==>src<==>hiercoarse.
    *
    * If dst is not a coarse interpolation level, we expect that what
    * ever generated dst also computed and cached the Connectors for
    * us to look up.  If not, then the PersistentOverlapConnectors
    * look-up mechanism can generate the Connectors at the cost of
    * scalability.
    */

   // Required hiercoarse<==>dst width for recursion:
   const hier::IntVector& hiercoarse_to_dst_width(
      d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln]);
   const hier::IntVector dst_to_hiercoarse_width(
      hier::Connector::convertHeadWidthToBase(
         dst_level->getBoxLevel()->getRefinementRatio(),
         hiercoarse_box_level.getRefinementRatio(),
         hiercoarse_to_dst_width));

   const hier::Connector* dst_to_hiercoarse = 0;
   std::shared_ptr<hier::Connector> bridged_dst_to_hiercoarse;
   std::shared_ptr<hier::Connector> bridged_hiercoarse_to_dst;

   /*
    * dst_to_hiercoarse and its transpose point to the
    * hiercoarse<==>dst Connector we will use.  The exact Connectors
    * will be either the ones cached in PersistentOverlapConnectors,
    * or if not there, the ones computed by the
    * dst<==>src<==>hiercoarse bridge.
    *
    * We expect to find cached hiercoarse<==>dst if:
    * - dst is on the hierarchy
    * - someone had the foresight to compute and cache it (such as the
    * coarsening step of the Richardson extrapolation).
    *
    * If there is no cached hiercoarse<==>dst, we'll bridge
    * dst<==>src<==>hiercoarse for it.  Therefore, we require src to
    * be on the hierarchy and satisfying certain nesting requirement
    * and having cached src<==>hiercoarse.
    *
    * An unscalable alternative of last resort would be to have the
    * PersistentOverlapConnectors compute the missing
    * dst<==>hiercoarse.  This is not implemented but could be in the
    * future.
    */

   bool has_cached_connectors =
      dst_level->hasConnector(*hiercoarse_level, dst_to_hiercoarse_width);
   has_cached_connectors = has_cached_connectors &&
      hiercoarse_level->hasConnector(*dst_level, hiercoarse_to_dst_width);

   if (has_cached_connectors) {

      dst_to_hiercoarse =
         &dst_level->findConnectorWithTranspose(*hiercoarse_level,
            dst_to_hiercoarse_width,
            hiercoarse_to_dst_width,
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

   } else {

      const hier::BoxLevel& src_box_level =
         *(d_src_level->getBoxLevel());
      if (hierarchy->getBoxLevel(next_coarser_ln + 1).get() !=
          &src_box_level) {
         TBOX_ERROR("Missing dst<==>hiercoarse connector and\n"
            << "src is not from hierarchy.  RefineSchedule cannot\n"
            << "continue because there is no way to connect\n"
            << "the destination to the hierarchy.");
      }

      const hier::BoxLevel& dst_box_level = dst_to_src.getBase();

      const hier::IntVector& hiercoarse_to_src_width =
         d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln];
      hier::IntVector src_to_hiercoarse_width =
         hiercoarse_to_src_width * d_src_level->getRatioToCoarserLevel();

      /*
       * The computation of src_to_hiercoarse_width puts it in the
       * resolution of next_coarser_ln+1, but that is not necessarily
       * the resolution of src (for example, src may be a Richardson
       * extrapolation temporary level), so we have to adjust.
       */
      if (d_src_level->getBoxLevel()->getRefinementRatio() <=
          hierarchy->getBoxLevel(next_coarser_ln + 1)->getRefinementRatio()) {
         src_to_hiercoarse_width *=
            d_src_level->getBoxLevel()->getRefinementRatio();
         src_to_hiercoarse_width /= hierarchy->getBoxLevel(
               next_coarser_ln + 1)->getRefinementRatio();
      } else if (d_src_level->getBoxLevel()->getRefinementRatio()
                 >=
                 hierarchy->getBoxLevel(next_coarser_ln
                    + 1)->getRefinementRatio()) {
         src_to_hiercoarse_width *= hierarchy->getBoxLevel(
               next_coarser_ln + 1)->getRefinementRatio();
         src_to_hiercoarse_width /=
            d_src_level->getBoxLevel()->getRefinementRatio();
      }

      const hier::Connector& src_to_hiercoarse =
         d_src_level->findConnectorWithTranspose(
            *hiercoarse_level,
            src_to_hiercoarse_width,
            hiercoarse_to_src_width,
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      const hier::IntVector& src_to_dst_width =
         dst_to_src.getTranspose().getConnectorWidth();
      const hier::IntVector& dst_to_src_width = dst_to_src.getConnectorWidth();

      bool do_simple_bridge;
      if (src_to_dst_width > zero_vec && dst_to_src_width > zero_vec) {
         do_simple_bridge = true;
      } else if (d_dst_level->getGridGeometry()->getNumberBlocks() == 1 &&
                 src_to_dst_width != zero_vec &&
                 dst_to_src_width != zero_vec) {
         /*
          * If a single-block domain is flat (width 1) in a direction with
          * zero connector width, we can still do the simple bridge for
          * bridged_dst_to_hiercoarse.
          */ 
         do_simple_bridge = true;
         hier::BoxContainer phys_domain(
            d_dst_level->getGridGeometry()->getPhysicalDomain());
         hier::Box bounding_box(phys_domain.getBoundingBox());
         hier::IntVector domain_size(bounding_box.numberCells());
         for (unsigned int d = 0; d < dim.getValue(); ++d) {
            if (src_to_dst_width[d] == 0 && domain_size[d] != 1) {
               do_simple_bridge = false;
            }
            if (dst_to_src_width[d] == 0 && domain_size[d] != 1) {
               do_simple_bridge = false;
            }
            if (!do_simple_bridge) {
               break;
            }
         }
      } else {
         do_simple_bridge = false;
      }
      
      if (do_simple_bridge) {

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->barrierAndStart();
         }

         oca.bridge(
            bridged_dst_to_hiercoarse,
            dst_to_src,
            src_to_hiercoarse,
            d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln],
            true);

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->stop();
         }
      } else {
         /*
          * This should be entered when setting up recursive schedules for
          * enhanced connectivity.
          */

         int level_number = d_src_level->getLevelNumber();

         const hier::Connector& found_dst_to_src =
            dst_box_level.findConnectorWithTranspose(
               src_box_level,
               hierarchy->getRequiredConnectorWidth(level_number, level_number),
               hierarchy->getRequiredConnectorWidth(level_number, level_number),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         const hier::Connector& found_src_to_hiercoarse =
            d_src_level->findConnectorWithTranspose(
               *hiercoarse_level,
               hierarchy->getRequiredConnectorWidth(level_number, level_number - 1),
               hierarchy->getRequiredConnectorWidth(level_number - 1, level_number),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->barrierAndStart();
         }

         oca.bridge(
            bridged_dst_to_hiercoarse,
            found_dst_to_src,
            found_src_to_hiercoarse,
            d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln],
            true);

         if (s_barrier_and_time) {
            t_bridge_dst_hiercoarse->stop();
         }

      }

      bridged_dst_to_hiercoarse->removePeriodicRelationships();
      bridged_dst_to_hiercoarse->getTranspose().removePeriodicRelationships();

      dst_to_hiercoarse = bridged_dst_to_hiercoarse.get();

   } /* !has_cached_connectors */

   hier::Connector& hiercoarse_to_dst = dst_to_hiercoarse->getTranspose();

   /*
    * Compute coarse_interp<==>hiercoarse by bridging
    * coarse_interp<==>dst<==>hiercoarse.
    */

   if (s_barrier_and_time) {
      t_bridge_coarse_interp_hiercoarse->barrierAndStart();
   }

   oca.bridge(
      coarse_interp_to_hiercoarse,
      dst_to_coarse_interp.getTranspose(),
      *dst_to_hiercoarse,
      hiercoarse_to_dst.getConnectorWidth(),
      true);
   hier::Connector& hiercoarse_to_coarse_interp =
      coarse_interp_to_hiercoarse->getTranspose();
   TBOX_ASSERT(coarse_interp_to_hiercoarse->getConnectorWidth() >= d_max_stencil_width);
   TBOX_ASSERT(hiercoarse_to_coarse_interp.getConnectorWidth() >= d_max_stencil_width);
   if (s_barrier_and_time) {
      t_bridge_coarse_interp_hiercoarse->stop();
   }

   if (d_num_periodic_directions > 0) {
      /*
       * Remove periodic relationships added by bridging.  Some of
       * them may be extraneous because coarse_interp may have parts
       * outside the domain boundary.  Then add periodic images for
       * coarse_interp and periodic relationships in
       * coarse_interp<==>hiercoarse.  Connector
       * hiercoarse--->hiercoarse must be as wide as hiercoarse--->dst
       * to make sure it sees everything dst sees.  Because dst sees
       * all of coarse_interp, this guarantees that
       * hiercoarse<==>coarse_interp does not miss any periodic
       * relationships.
       */
      coarse_interp_to_hiercoarse->removePeriodicRelationships();
      hiercoarse_to_coarse_interp.removePeriodicRelationships();

      const hier::Connector& hiercoarse_to_hiercoarse =
         hiercoarse_level->findConnector(*hiercoarse_level,
            hiercoarse_to_dst.getConnectorWidth(),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);
      edge_utils.addPeriodicImagesAndRelationships(
         *coarse_interp_box_level,
         *coarse_interp_to_hiercoarse,
         hierarchy->getGridGeometry()->getDomainSearchTree(),
         hiercoarse_to_hiercoarse);

   }

   if (s_extra_debug) {
      sanityCheckCoarseInterpAndHiercoarseLevels(
         next_coarser_ln,
         hierarchy,
         *coarse_interp_to_hiercoarse);
   }

   /*
    * Construct the coarse interpolation PatchLevel and reset
    * coarse_interp<==>hiercoarse connectors to use the PatchLevel's
    * BoxLevel.
    */

   coarse_interp_level.reset(new hier::PatchLevel(
         coarse_interp_box_level,
         hiercoarse_level->getGridGeometry(),
         hiercoarse_level->getPatchDescriptor()));
   coarse_interp_level->setLevelNumber(next_coarser_ln);
   coarse_interp_level->setNextCoarserHierarchyLevelNumber(next_coarser_ln - 1);

   if (hiercoarse_level->getGridGeometry()->getNumberBlocks() > 1) {
      hiercoarse_level->getGridGeometry()->
      adjustMultiblockPatchLevelBoundaries(*coarse_interp_level);
   }

}

/*
 **************************************************************************
 * Check that the Connectors between the coarse interpolation and
 * hiercoarse levels are transposes and that the coarse interpolation
 * BoxLevel sufficiently nests inside the hiercoarse.
 **************************************************************************
 */
void
RefineSchedule::sanityCheckCoarseInterpAndHiercoarseLevels(
   const int next_coarser_ln,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::Connector& coarse_interp_to_hiercoarse)
{
   TBOX_ASSERT(coarse_interp_to_hiercoarse.hasTranspose());

   const hier::Connector& hiercoarse_to_coarse_interp =
      coarse_interp_to_hiercoarse.getTranspose();

   /*
    * coarse_interp_to_hiercoarse and hiercoarse_to_coarse_interp should be
    * proper transposes.
    */
   size_t err1 = coarse_interp_to_hiercoarse.checkTransposeCorrectness(
         hiercoarse_to_coarse_interp);
   if (err1) tbox::perr
      << "coarse_interp_to_hiercoarse failed transpose correctness."
      << std::endl;
   size_t err2 = hiercoarse_to_coarse_interp.checkTransposeCorrectness(
         coarse_interp_to_hiercoarse);
   if (err2) tbox::perr
      << "hiercoarse_to_coarse_interp failed transpose correctness."
      << std::endl;

   hier::OverlapConnectorAlgorithm oca;
   oca.setTimerPrefix("xfer::RefineSchedule_build");

   hier::IntVector multi_max_stencil(d_max_stencil_width);

   /*
    * To work properly, we must ensure that
    * coarse_interp^d_max_stencil_width nests in
    * hiercoarse^fine_connector_width.
    *
    * The nesting guarantees that the Connectors sees enough of its
    * surroundings to generate all the necessary relationships in
    * further RefineSchedule recursions.  We know that
    * coarse_interp_to_hiercoarse and hiercoarse_to_coarse_interp are not
    * complete, but this nesting guarantees we still have enough relationships
    * to avoid missing any relationships when we use these Connectors in
    * bridge operations.
    */
   std::shared_ptr<hier::Connector> wider_coarse_interp_to_hiercoarse;
   oca.findOverlaps(wider_coarse_interp_to_hiercoarse,
      coarse_interp_to_hiercoarse.getBase(),
      coarse_interp_to_hiercoarse.getHead(),
      d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln] - multi_max_stencil);

   std::shared_ptr<hier::BoxLevel> external;
   std::shared_ptr<hier::MappingConnector> coarse_interp_to_external;
   hier::BoxLevelConnectorUtils edge_utils;
   edge_utils.computeExternalParts(
      external,
      coarse_interp_to_external,
      *wider_coarse_interp_to_hiercoarse,
      d_top_refine_schedule->d_fine_connector_widths[next_coarser_ln] - multi_max_stencil,
      hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());
   coarse_interp_to_external->eraseEmptyNeighborSets();

   int err3 = coarse_interp_to_external->getGlobalNumberOfRelationships();
   if (err3) {
      tbox::perr << "Some parts of coarse_interp lies outside of where we\n"
                 << "guarantee support for recursive RefineSchedule.\n"
                 << coarse_interp_to_external->format("SE: ", 2);
   }

   if (err1 || err2 || err3) {
      TBOX_ERROR(
         "coarse_interp<==>hiercoarse have problems as reported above.\n"
         << "coarse_interp:\n" << coarse_interp_to_hiercoarse.getBase().format("SUPP->", 2)
         << "hiercoarse:\n" << coarse_interp_to_hiercoarse.getHead().format("HCRS->", 2)
         << "dst_to_coarse_interp:\n" << d_dst_to_coarse_interp->format("DS->", 2)
         << "coarse_interp_to_hiercoarse:\n" << coarse_interp_to_hiercoarse.format("SH->", 2)
         << "hiercoarse_to_coarse_interp:\n" << hiercoarse_to_coarse_interp.format("HS->", 2));
   }

}

/*
 **************************************************************************
 *
 * Execute the stored communication schedule that copies data into the
 * destination component of the destination level.  The algorithm is as
 * follows:
 *
 * (1) Allocate scratch space on the destination level if it does
 *    not exist.
 * (2) Call the recursive fill function that will fill the scratch
 *    space of the destination level.
 * (3) If the scratch and destination spaces are not the same,
 *    then copy locally from the scratch into the destination.
 * (4) Deallocate any previously allocated data.
 *
 **************************************************************************
 */

void
RefineSchedule::fillData(
   double fill_time,
   bool do_physical_boundary_fill) const
{
  RANGE_PUSH("fillData", 1);
   if (s_barrier_and_time) {
      t_fill_data->barrierAndStart();
   }

   t_fill_data_nonrecursive->start();

   if (d_internal_allocated) {
      setInternalDataTime(fill_time);
   }

   /*
    * Set the refine items and time for all transactions.  These items will
    * be shared by all transaction objects in the communication schedule.
    */

   d_transaction_factory->setTransactionTime(fill_time);

   /*
    * Check whether scratch data needs to be allocated on the destination
    * level.  Keep track of those allocated components so that they may be
    * deallocated later.
    */

   hier::ComponentSelector allocate_vector;
   allocateScratchSpace(allocate_vector, d_dst_level, fill_time);

   hier::ComponentSelector encon_allocate_vector;
   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      allocateScratchSpace(encon_allocate_vector, d_encon_level, fill_time);
   }

   hier::ComponentSelector nbr_fill_scratch_vector;
   hier::ComponentSelector nbr_fill_dst_vector;
   if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 && 
       d_nbr_blk_fill_level.get()) {
      allocateScratchSpace(nbr_fill_scratch_vector,
                           d_nbr_blk_fill_level,
                           fill_time);
      allocateDestinationSpace(nbr_fill_dst_vector,
                               d_nbr_blk_fill_level,
                               fill_time);
   }

   /*
    * Begin the recursive algorithm that fills from coarser, fills from
    * same, and then fills physical boundaries.
    */

   t_fill_data_nonrecursive->stop();
   t_fill_data_recursive->start();
   recursiveFill(fill_time, do_physical_boundary_fill);
   t_fill_data_recursive->stop();
   t_fill_data_nonrecursive->start();

   /*
    * Copy the scratch space of the destination level to the destination
    * space.
    */

   bool copied = copyScratchToDestination();
#if defined(HAVE_RAJA)
   if (copied) {
      tbox::parallel_synchronize();
   }
#else
   NULL_USE(copied);
#endif

   /*
    * Deallocate any allocated scratch space on the destination level.
    */

   d_dst_level->deallocatePatchData(allocate_vector);

   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      d_encon_level->deallocatePatchData(encon_allocate_vector);
   }
   if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
       d_nbr_blk_fill_level.get()) {
      d_nbr_blk_fill_level->deallocatePatchData(nbr_fill_scratch_vector);
      d_nbr_blk_fill_level->deallocatePatchData(nbr_fill_dst_vector);
   }

   t_fill_data_nonrecursive->stop();

   if (s_barrier_and_time) {
      t_fill_data->stop();
   }
   RANGE_POP;
}

/*
 **************************************************************************
 *
 * Recursively fill the required fill boxes on the destination level.
 * The scratch component of the level will be filled.  The algorithm
 * is as follows:
 *
 * (1) If we need to get data from a coarser level, then:
 *   (a) allocate scratch data on the coarser level
 *   (b) recursively call this routine to get the data
 *   (c) refine data from the coarse level into this level
 *   (d) deallocate the scratch data on the coarser level
 * (2) Copy data from the same level of refinement
 * (3) Copy data from the physical boundary conditions
 *
 **************************************************************************
 */

void
RefineSchedule::recursiveFill(
   double fill_time,
   bool do_physical_boundary_fill) const
{
   /*
    * Copy data from the source interiors of the source level into the ghost
    * cells and interiors of the scratch space on the destination level
    * for data where coarse data takes priority on level boundaries.
    */
   d_coarse_priority_level_schedule->communicate();

   /*
    * If there is a coarser schedule stored in this object, then we will
    * need to get data from a coarser grid level.
    */

   if (d_coarse_interp_schedule) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::ComponentSelector allocate_vector;
      hier::ComponentSelector work_allocate_vector;
      allocateScratchSpace(allocate_vector, d_coarse_interp_level, fill_time);
      allocateWorkSpace(work_allocate_vector, d_coarse_interp_level, fill_time);

      hier::ComponentSelector encon_allocate_vector;
      hier::ComponentSelector encon_work_allocate_vector;
      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         allocateScratchSpace(encon_allocate_vector,
            d_coarse_interp_schedule->d_encon_level,
            fill_time);
         allocateWorkSpace(encon_work_allocate_vector,
            d_coarse_interp_schedule->d_encon_level,
            fill_time);
      }

      hier::ComponentSelector nbr_blk_fill_allocate_vector;
      hier::ComponentSelector nbr_blk_fill_work_allocate_vector;
      if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
          d_coarse_interp_schedule->d_nbr_blk_fill_level.get()) {
         allocateScratchSpace(nbr_blk_fill_allocate_vector,
            d_coarse_interp_schedule->d_nbr_blk_fill_level, fill_time);
         allocateWorkSpace(nbr_blk_fill_work_allocate_vector,
            d_coarse_interp_schedule->d_nbr_blk_fill_level, fill_time);
      }

      hier::ComponentSelector nbr_blk_fill_scratch_vector;
      hier::ComponentSelector nbr_blk_fill_work_vector;
      hier::ComponentSelector nbr_blk_fill_dst_vector;
      if (d_nbr_blk_fill_level.get()) {
         allocateScratchSpace(nbr_blk_fill_scratch_vector,
            d_nbr_blk_fill_level, fill_time);
         allocateScratchSpace(nbr_blk_fill_work_vector,
            d_nbr_blk_fill_level, fill_time);
         allocateDestinationSpace(nbr_blk_fill_dst_vector,
            d_nbr_blk_fill_level, fill_time);
      }

      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */
      d_coarse_interp_schedule->recursiveFill(fill_time,
         do_physical_boundary_fill);

#if defined(HAVE_RAJA)
//      TODO: This sync probably isn't necessary, but keep an eye on it. 
//      tbox::parallel_synchronize();
#endif

      /*
       * d_coarse_interp_level should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData(d_dst_level,
         d_coarse_interp_level,
         d_dst_to_coarse_interp->getTranspose(),
         *d_coarse_interp_to_unfilled,
         d_refine_overlaps);


      /*
       * Deallocate the scratch data from the coarse grid.
       */

      d_coarse_interp_level->deallocatePatchData(allocate_vector);
      d_coarse_interp_level->deallocatePatchData(work_allocate_vector);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_schedule->d_encon_level->deallocatePatchData(
            encon_allocate_vector);
         d_coarse_interp_schedule->d_encon_level->deallocatePatchData(
            encon_work_allocate_vector);
      }

      if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
          d_coarse_interp_schedule->d_nbr_blk_fill_level.get()) {
         d_coarse_interp_schedule->d_nbr_blk_fill_level->deallocatePatchData(
            nbr_blk_fill_allocate_vector);
         d_coarse_interp_schedule->d_nbr_blk_fill_level->deallocatePatchData(
            nbr_blk_fill_work_allocate_vector);
      }

      if (d_nbr_blk_fill_level.get()) {
         d_nbr_blk_fill_level->deallocatePatchData(
            nbr_blk_fill_scratch_vector);
         d_nbr_blk_fill_level->deallocatePatchData(
            nbr_blk_fill_work_vector);
         d_nbr_blk_fill_level->deallocatePatchData(
            nbr_blk_fill_dst_vector);
      }

   }

   if (d_coarse_interp_encon_schedule) {

      /*
       * Allocate data on the coarser level and keep track of the allocated
       * components so that they may be deallocated later.
       */

      hier::ComponentSelector allocate_vector;
      hier::ComponentSelector work_allocate_vector;
      allocateScratchSpace(allocate_vector,
                           d_coarse_interp_encon_level,
                           fill_time);
      allocateWorkSpace(work_allocate_vector,
                           d_coarse_interp_encon_level,
                           fill_time);

      hier::ComponentSelector encon_allocate_vector;
      hier::ComponentSelector encon_work_allocate_vector;
      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         allocateScratchSpace(encon_allocate_vector,
            d_coarse_interp_encon_schedule->d_encon_level,
            fill_time);
         allocateWorkSpace(encon_work_allocate_vector,
            d_coarse_interp_encon_schedule->d_encon_level,
            fill_time);
      }

      /*
       * Recursively call the fill routine to fill the required coarse fill
       * boxes on the coarser level.
       */

      d_coarse_interp_encon_schedule->recursiveFill(fill_time,
         do_physical_boundary_fill);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      /*
       * d_coarse_interp_encon_level should now be filled.  Now interpolate
       * data from the coarse grid into the fine grid.
       */

      refineScratchData(d_encon_level,
         d_coarse_interp_encon_level,
         d_encon_to_coarse_interp_encon->getTranspose(),
         *d_coarse_interp_encon_to_unfilled_encon,
         d_encon_refine_overlaps);

      /*
       * Deallocate the scratch data from the coarse grid.
       */

      d_coarse_interp_encon_level->deallocatePatchData(allocate_vector);
      d_coarse_interp_encon_level->deallocatePatchData(work_allocate_vector);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_encon_schedule->d_encon_level->deallocatePatchData(
            encon_allocate_vector);
         d_coarse_interp_encon_schedule->d_encon_level->deallocatePatchData(
            encon_work_allocate_vector);
      }

   }

   /*
    * Copy data from the source interiors of the source level into the ghost
    * cells and interiors of the scratch space on the destination level
    * for data where fine data takes priority on level boundaries.
    */
   d_fine_priority_level_schedule->communicate();

   /*
    * Fill the physical boundaries of the scratch space on the destination
    * level.
    */

   if (do_physical_boundary_fill || d_force_boundary_fill) {
      fillPhysicalBoundaries(fill_time);
   }

   if (d_dst_level->getGridGeometry()->getNumberOfBlockSingularities() > 0) {
      fillSingularityBoundaries(fill_time);
   }
}

/*
 **************************************************************************
 *
 * Fill the physical boundaries of the specified level with data.
 *
 **************************************************************************
 */

void
RefineSchedule::fillPhysicalBoundaries(
   double fill_time) const
{
   TBOX_ASSERT(d_dst_level);
   t_fill_physical_boundaries->start();

   d_dst_level->setBoundaryBoxes();

   if (d_refine_patch_strategy) {
#if defined(HAVE_RAJA)
      bool bdry_is_filled = false; 
#endif 
      for (hier::PatchLevel::iterator p(d_dst_level->begin());
           p != d_dst_level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch(*p);
         if (patch->getPatchGeometry()->intersectsPhysicalBoundary()) {
            d_refine_patch_strategy->
            setPhysicalBoundaryConditions(*patch,
               fill_time,
               d_boundary_fill_ghost_width);
#if defined(HAVE_RAJA)
            bdry_is_filled = true;
#endif
         }
      }
#if defined(HAVE_RAJA) 
      if (bdry_is_filled && d_refine_patch_strategy->needSynchronize()) {
         tbox::parallel_synchronize();
      }
#endif

   }
   t_fill_physical_boundaries->stop();
}

/*
 ********************************************************************
 *
 * Call patch strategy method for filling ghost regions around
 * enhanced connectivity block boundaries.
 *
 ********************************************************************
 */
void
RefineSchedule::fillSingularityBoundaries(
   double fill_time) const
{
   TBOX_ASSERT(d_dst_level);
   t_fill_singularity_boundaries->start();

   NULL_USE(fill_time);

   std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   if (grid_geometry->getNumberBlocks() > 1) {

      d_dst_level->setBoundaryBoxes();

      const tbox::Dimension& dim(d_dst_level->getDim());

      const hier::IntVector& ratio = d_dst_level->getRatioToLevelZero();

      if (d_singularity_patch_strategy) {
#if defined(HAVE_RAJA)
         bool sing_is_filled = false;
#endif
         for (hier::BlockId::block_t bn = 0; bn < grid_geometry->getNumberBlocks(); ++bn) {

            hier::BlockId block_id(bn);

            const hier::BoxContainer& sing_boxes =
               grid_geometry->getSingularityBoxContainer(block_id);
            for (hier::BoxContainer::const_iterator sb = sing_boxes.begin();
                 sb != sing_boxes.end(); ++sb) {

               hier::Box singularity(*sb);
               singularity.refine(ratio);

               const hier::BoxContainer& level_boxes(
                  d_dst_level->getBoxLevel()->getBoxes());
               hier::BoxContainerSingleBlockIterator dst_local_iter(
                  level_boxes.begin(block_id));

               for ( ; dst_local_iter != level_boxes.end(block_id);
                     ++dst_local_iter) {

                  const hier::BoxId& box_id = dst_local_iter->getBoxId();

                  std::shared_ptr<hier::Patch> patch(
                     d_dst_level->getPatch(box_id));
                  std::shared_ptr<hier::PatchGeometry> pgeom(
                     patch->getPatchGeometry());

                  const std::vector<hier::BoundaryBox>& nboxes =
                     pgeom->getNodeBoundaries();

                  if (nboxes.size()) {
                     for (int bb = 0;
                          bb < static_cast<int>(nboxes.size()); ++bb) {
                        hier::Box intersection((nboxes[bb].getBox())
                                               * singularity);
                        if (!(intersection.empty())) {
                           hier::Box fill_box(
                              pgeom->getBoundaryFillBox(nboxes[bb],
                                 patch->getBox(),
                                 d_boundary_fill_ghost_width));

                           if (!(fill_box.empty())) {
                              d_singularity_patch_strategy->
                              fillSingularityBoundaryConditions(
                                 *patch, *d_encon_level,
                                 d_dst_to_encon,
                                 fill_box, nboxes[bb],
                                 grid_geometry);
#if defined(HAVE_RAJA)
                              sing_is_filled = true;
#endif
                           }
                        }
                     }
                  }

                  if (dim == tbox::Dimension(3)) {
                     const std::vector<hier::BoundaryBox>& eboxes =
                        pgeom->getEdgeBoundaries();

                     if (eboxes.size()) {
                        for (int bb = 0;
                             bb < static_cast<int>(eboxes.size()); ++bb) {
                           hier::Box intersection(
                              (eboxes[bb].getBox()) * singularity);
                           if (!(intersection.empty())) {
                              hier::Box fill_box(
                                 pgeom->getBoundaryFillBox(eboxes[bb],
                                    patch->getBox(),
                                    d_boundary_fill_ghost_width));

                              if (!(fill_box.empty())) {
                                 d_singularity_patch_strategy->
                                 fillSingularityBoundaryConditions(
                                    *patch, *d_encon_level,
                                    d_dst_to_encon,
                                    fill_box, eboxes[bb],
                                    grid_geometry);
#if defined(HAVE_RAJA)
                                 sing_is_filled = true;
#endif
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
#if defined(HAVE_RAJA) 
         if (sing_is_filled && d_singularity_patch_strategy->needSynchronize()) {
            tbox::parallel_synchronize();
         }
#endif
      }
   }
   t_fill_singularity_boundaries->stop();
}

/*
 **************************************************************************
 *
 * Check whether the scratch data needs to be allocated on the specified
 * level.  Keep track of those allocated components so that they may be
 * deallocated later.
 *
 **************************************************************************
 */

void
RefineSchedule::allocateScratchSpace(
   hier::ComponentSelector& allocate_vector,
   const std::shared_ptr<hier::PatchLevel>& level,
   double fill_time) const
{
   TBOX_ASSERT(level);

   allocate_vector.clrAllFlags();

   hier::ComponentSelector preprocess_vector;

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      if (!level->checkAllocated(scratch_id)) {
         allocate_vector.setFlag(scratch_id);
      }
      preprocess_vector.setFlag(scratch_id);

   }

   level->allocatePatchData(allocate_vector, fill_time);

   if (d_transaction_factory) {
      d_transaction_factory->preprocessScratchSpace(level,
         fill_time,
         preprocess_vector);
   }
}

void
RefineSchedule::allocateDestinationSpace(
   hier::ComponentSelector& allocate_vector,
   const std::shared_ptr<hier::PatchLevel>& level,
   double fill_time) const
{
   TBOX_ASSERT(level);

   allocate_vector.clrAllFlags();

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      const int dst_id = d_refine_items[iri]->d_dst;
      if (!level->checkAllocated(dst_id)) {
         allocate_vector.setFlag(dst_id);
      }
   }

   level->allocatePatchData(allocate_vector, fill_time);

}

void
RefineSchedule::allocateWorkSpace(
   hier::ComponentSelector& allocate_vector,
   const std::shared_ptr<hier::PatchLevel>& level,
   double fill_time) const
{
   TBOX_ASSERT(level);

   allocate_vector.clrAllFlags();

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      for (std::vector<int>::const_iterator it =
               d_refine_items[iri]->d_work.begin();
            it != d_refine_items[iri]->d_work.end();
            ++it)
      {
         const int work_id = *it;
         if (!level->checkAllocated(work_id)) {
            allocate_vector.setFlag(work_id);
         }
      }
   }

   level->allocatePatchData(allocate_vector, fill_time);

}


/*
 **************************************************************************
 *
 * Copy data from the scratch space of the specified level into the
 * destination space.  If the two spaces are the same, then no copy
 * is performed.
 *
 **************************************************************************
 */

bool
RefineSchedule::copyScratchToDestination() const
{
   TBOX_ASSERT(d_dst_level);
   bool copied = false;

   for (hier::PatchLevel::iterator p(d_dst_level->begin());
        p != d_dst_level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch(*p);

      for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
         const int src_id = d_refine_items[iri]->d_scratch;
         const int dst_id = d_refine_items[iri]->d_dst;
         if (src_id != dst_id) {
            TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(patch->
                  getPatchData(dst_id)->getTime(),
                  patch->getPatchData(src_id)->getTime()));
            patch->getPatchData(dst_id)->copy(*patch->getPatchData(src_id));
            copied = true;
         }
      }
   }

   return copied;
}

/*
 **************************************************************************
 *
 * Refine data from the coarse level into the fine level on the provided
 * fill box regions.  All operations are performed on the scratch space.
 *
 **************************************************************************
 */

void
RefineSchedule::refineScratchData(
   const std::shared_ptr<hier::PatchLevel>& fine_level,
   const std::shared_ptr<hier::PatchLevel>& coarse_level,
   const hier::Connector& coarse_to_fine,
   const hier::Connector& coarse_to_unfilled,
   const std::vector<std::vector<std::shared_ptr<hier::BoxOverlap> > >&
   overlaps) const
{
   t_refine_scratch_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   bool is_encon = (fine_level == d_encon_level);
#endif
   int nbr_blk_copies = 0;

   if (d_refine_patch_strategy) {
      d_refine_patch_strategy->preprocessRefineLevel(
         *fine_level,
         *coarse_level,
         coarse_to_fine,
         coarse_to_unfilled,
         overlaps,
         d_refine_items);
#if defined(HAVE_RAJA)
      if (d_refine_patch_strategy->needSynchronize()) {
         tbox::parallel_synchronize();
      }
#endif
   }

   const hier::IntVector ratio(fine_level->getRatioToLevelZero()
                               / coarse_level->getRatioToLevelZero());

   /*
    * Loop over all the coarse patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   for (int pi = 0; pi < coarse_level->getLocalNumberOfPatches(); ++pi) {
      const hier::Box& crse_box = coarse_level->getPatch(pi)->getBox();
      const hier::BoxId& crse_box_id = crse_box.getBoxId();

      hier::Connector::ConstNeighborhoodIterator dst_nabrs =
         coarse_to_fine.find(crse_box_id);
      const hier::Box& dst_box = *coarse_to_fine.begin(dst_nabrs);
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each crse_box can point back to just one dst_box.
       * All other boxes in dst_nabrs must be a periodic image of
       * the same dst_box.
       */
      for (hier::Connector::ConstNeighborIterator na = coarse_to_fine.begin(dst_nabrs);
           na != coarse_to_fine.end(dst_nabrs); ++na) {
         TBOX_ASSERT(na->isPeriodicImage() ||
            na == coarse_to_fine.begin(dst_nabrs));
         TBOX_ASSERT(na->getGlobalId() == dst_box.getGlobalId());
      }
#endif
      std::shared_ptr<hier::Patch> fine_patch(fine_level->getPatch(
                                                   dst_box.getGlobalId()));
      std::shared_ptr<hier::Patch> crse_patch(coarse_level->getPatch(
                                                   crse_box.getGlobalId()));

      const hier::BlockId& crse_blk_id = crse_patch->getBox().getBlockId();
      hier::IntVector local_ratio(ratio.getBlockVector(crse_blk_id));
      if (fine_patch->getBox().getBlockId() == crse_blk_id) {

         TBOX_ASSERT(coarse_to_unfilled.numLocalNeighbors(crse_box.getBoxId()) == 1);
         hier::Connector::ConstNeighborhoodIterator unfilled_nabrs =
            coarse_to_unfilled.find(crse_box.getBoxId());
         const hier::Box& unfilled_nabr =
            *coarse_to_unfilled.begin(unfilled_nabrs);
         hier::BoxContainer fill_boxes(unfilled_nabr);

         if (d_refine_patch_strategy) {
            d_refine_patch_strategy->preprocessRefineBoxes(*fine_patch,
               *crse_patch,
               fill_boxes,
               local_ratio);
#if defined(HAVE_RAJA)
            if (d_refine_patch_strategy->needSynchronize()) {
               tbox::parallel_synchronize();
            }
#endif
         }

         for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
            const RefineClasses::Data * const ref_item = d_refine_items[iri];
            if (ref_item->d_oprefine) {

               std::shared_ptr<hier::BoxOverlap> refine_overlap =
                  (overlaps[pi])[ref_item->d_class_index];

               const int scratch_id = ref_item->d_scratch;

               ref_item->d_oprefine->refine(*fine_patch, *crse_patch,
                  scratch_id, scratch_id,
                  *refine_overlap, local_ratio);
            }
         }

         if (d_refine_patch_strategy) {
            d_refine_patch_strategy->setPostRefineSyncFlag();
#if defined(HAVE_RAJA)
            if (d_refine_patch_strategy->needSynchronize()) {
               tbox::parallel_synchronize();
            }
#endif
            d_refine_patch_strategy->postprocessRefineBoxes(*fine_patch,
               *crse_patch,
               fill_boxes,
               local_ratio);
#if defined(HAVE_RAJA)
            if (d_refine_patch_strategy->needSynchronize()) {
               tbox::parallel_synchronize();
            }
#endif
         }

      } else {
         /*
          * This section is only entered when filling ghost regions in
          * blocks neighboring the fine patch, and there is anisotropic
          * refinement so that the refinement ratio on the neighboring block
          * may be different from the ratio on the fine patch's block.
          */

         TBOX_ASSERT(!is_encon);
         TBOX_ASSERT(!d_dst_level->getGridGeometry()->hasIsotropicRatios());
         TBOX_ASSERT(d_coarse_interp_to_nbr_fill->numLocalNeighbors(crse_box.getBoxId()) == 1);
         hier::Connector::ConstNeighborhoodIterator unfilled_nabrs =
            d_coarse_interp_to_nbr_fill->find(crse_box.getBoxId());
         const hier::Box& unfilled_nabr =
            *d_coarse_interp_to_nbr_fill->begin(unfilled_nabrs);
         hier::BoxContainer fill_boxes(unfilled_nabr);

         const hier::BoxId& unfilled_id = unfilled_nabr.getBoxId();

         /*
          * The refinement operation interpolates data onto nbr_fill_patch.
          */

         std::shared_ptr<hier::Patch> nbr_fill_patch(
            d_nbr_blk_fill_level->getPatch(unfilled_id));

         if (d_refine_patch_strategy) {
            d_refine_patch_strategy->preprocessRefineBoxes(*nbr_fill_patch,
               *crse_patch,
               fill_boxes,
               local_ratio);
#if defined(HAVE_RAJA)
            if (d_refine_patch_strategy->needSynchronize()) {
               tbox::parallel_synchronize();
            }
#endif
         }

         for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
            const RefineClasses::Data * const ref_item = d_refine_items[iri];

            if (ref_item->d_oprefine) {

               std::shared_ptr<hier::BoxOverlap> refine_overlap =
                  (overlaps[pi])[ref_item->d_class_index];

               const int scratch_id = ref_item->d_scratch;

               ref_item->d_oprefine->refine(*nbr_fill_patch, *crse_patch,
                  scratch_id, scratch_id,
                  *refine_overlap, local_ratio);

            }
         }
#if defined(HAVE_RAJA)
         tbox::parallel_synchronize();
#endif

         if (d_refine_patch_strategy) {
            d_refine_patch_strategy->postprocessRefineBoxes(*nbr_fill_patch,
               *crse_patch,
               fill_boxes,
               local_ratio);
#if defined(HAVE_RAJA)
            if (d_refine_patch_strategy->needSynchronize()) {
               tbox::parallel_synchronize();
            }
#endif
         }


         /*
          * Post-interpolation loop to copy data from nbr_fill_patch to
          * fine_patch.
          */

         for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
            const RefineClasses::Data * const ref_item = d_refine_items[iri];

            if (ref_item->d_oprefine) {
               std::shared_ptr<hier::BoxOverlap> nbr_copy_overlap =
                  (d_nbr_blk_copy_overlaps[nbr_blk_copies])[ref_item->d_class_index];

               const int scratch_id = ref_item->d_scratch;

               fine_patch->getPatchData(scratch_id)->copy(
                  *nbr_fill_patch->getPatchData(scratch_id), *nbr_copy_overlap);

            }
         }

         ++nbr_blk_copies;
      }
   }

#if defined(HAVE_RAJA)
   if (!d_refine_patch_strategy || d_refine_patch_strategy->needSynchronize()) {
      tbox::parallel_synchronize();
   }
#endif

   if (d_refine_patch_strategy) {
      d_refine_patch_strategy->postprocessRefineLevel(
         *fine_level,
         *coarse_level,
         coarse_to_fine,
         coarse_to_unfilled);
#if defined(HAVE_RAJA)
      if (d_refine_patch_strategy->needSynchronize()) {
         tbox::parallel_synchronize();
      }
#endif

   }

   t_refine_scratch_data->stop();
}

/*
 **************************************************************************
 *
 * Compute the overlaps defining where the data will be refined into the
 * fine level.  All operations are performed on the scratch space.
 *
 **************************************************************************
 */
void
RefineSchedule::computeRefineOverlaps(
   std::vector<std::vector<std::shared_ptr<hier::BoxOverlap> > >& overlaps,
   const std::shared_ptr<hier::PatchLevel>& fine_level,
   const std::shared_ptr<hier::PatchLevel>& coarse_level,
   const hier::Connector& coarse_to_fine,
   const hier::Connector& coarse_to_unfilled)
{
   bool is_encon = (fine_level == d_encon_level);

   std::shared_ptr<BoxGeometryVariableFillPattern> bg_fill_pattern(
      std::make_shared<BoxGeometryVariableFillPattern>());

   std::shared_ptr<hier::PatchDescriptor> fine_patch_descriptor(
      fine_level->getPatchDescriptor());

   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

   TBOX_ASSERT(overlaps.empty());
   overlaps.reserve(coarse_level->getLocalNumberOfPatches());

   std::shared_ptr<hier::BaseGridGeometry> grid_geom(
      d_dst_level->getGridGeometry());

   /*
    * Loop over all the coarse patches and find the corresponding
    * destination patch and destination fill boxes.
    */

   for (hier::PatchLevel::iterator crse_itr(coarse_level->begin());
        crse_itr != coarse_level->end(); ++crse_itr) {
      const hier::Box& coarse_box = crse_itr->getBox();
      hier::Connector::ConstNeighborhoodIterator fine_nabrs =
         coarse_to_fine.find(coarse_box.getBoxId());
      hier::Connector::ConstNeighborIterator na =
         coarse_to_fine.begin(fine_nabrs);
      const hier::Box& fine_box = *na;
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Each coarse_box can point back to just one fine_box.
       * All other boxes in fine_nabrs must be a periodic image of
       * the same fine_box.
       */
      for ( ; na != coarse_to_fine.end(fine_nabrs); ++na) {
         TBOX_ASSERT(na->isPeriodicImage() ||
            na == coarse_to_fine.begin(fine_nabrs));
         TBOX_ASSERT(na->getGlobalId() == fine_box.getGlobalId());
      }
#endif
      std::shared_ptr<hier::Patch> fine_patch(fine_level->getPatch(
                                                   fine_box.getBoxId()));

      TBOX_ASSERT(coarse_to_unfilled.numLocalNeighbors(
            coarse_box.getBoxId()) == 1);
      hier::Connector::ConstNeighborhoodIterator unfilled_nabrs =
         coarse_to_unfilled.find(coarse_box.getBoxId());
      const hier::Box& unfilled_nabr = *coarse_to_unfilled.begin(
            unfilled_nabrs);

      if (unfilled_nabr.getBlockId() == coarse_box.getBlockId()) {
         hier::BoxContainer fill_boxes(unfilled_nabr);
         fill_boxes.intersectBoxes(
            hier::Box::refine(coarse_box, coarse_to_fine.getRatio()));
         const hier::BoxId& unfilled_id = unfilled_nabr.getBoxId();

         hier::BoxContainer node_fill_boxes;
         if (!is_encon) {
            d_unfilled_to_unfilled_node->getNeighborBoxes(unfilled_id,
               node_fill_boxes);
         }

         /*
          * The refine overlap will cover only  the fine fill box regions.
          * of index space.  Note that we restrict the interpolation range
          * to the intersection of the fine fill box and the ghost box of
          * the scratch data component (i.e., the destination of the
          * interpolation).  This is needed for the case where data
          * components treated by the schedule have different ghost
          * cell widths since the fill boxes are generated using the
          * maximum ghost cell width.
          */
         overlaps.push_back(std::vector<std::shared_ptr<hier::BoxOverlap> >(0));
         std::vector<std::shared_ptr<hier::BoxOverlap> >& refine_overlaps(
            overlaps.back());
         refine_overlaps.resize(num_equiv_classes);
         for (int ne = 0; ne < num_equiv_classes; ++ne) {

            const RefineClasses::Data& rep_item =
               d_refine_classes->getClassRepresentative(ne);

            if (rep_item.d_oprefine) {

               const int scratch_id = rep_item.d_scratch;
               std::shared_ptr<hier::PatchDataFactory> fine_pdf(
                  fine_patch_descriptor->getPatchDataFactory(scratch_id));
               const hier::IntVector& fine_scratch_gcw =
                  fine_pdf->getGhostCellWidth();
               hier::Box scratch_space(fine_patch->getBox());
               scratch_space.grow(fine_scratch_gcw);

               if (!is_encon) {
                  refine_overlaps[ne] =
                     rep_item.d_var_fill_pattern->computeFillBoxesOverlap(
                        fill_boxes,
                        node_fill_boxes,
                        fine_patch->getBox(),
                        scratch_space,
                        *fine_pdf);

               } else {
                  refine_overlaps[ne] =
                     bg_fill_pattern->computeFillBoxesOverlap(
                        fill_boxes,
                        node_fill_boxes,
                        fine_patch->getBox(),
                        scratch_space,
                        *fine_pdf);
               }
            }
         }
      } else {
         /*
          * This section only gets entered if fine patch touches a block
          * boundary and we need to fill ghost region across that block
          * boundary, and there is anisotropic refinement so that the
          * refinement ratio on the neighboring block may be different from
          * the ratio on the fine patch's block.
          */
         TBOX_ASSERT(!is_encon);
         TBOX_ASSERT(!grid_geom->hasIsotropicRatios()); 
         TBOX_ASSERT(d_coarse_interp_to_nbr_fill->numLocalNeighbors(coarse_box.getBoxId()) == 1);

         /*
          * Find the patch on d_nbr_blk_fill_level where data will
          * be interpolated, and the fill boxes to ensure that data is only
          * interpolated from space covered by the coarse box.
          */
         hier::Connector::ConstNeighborhoodIterator fill_nabrs =
            d_coarse_interp_to_nbr_fill->find(coarse_box.getBoxId());
         const hier::Box& fill_nabr =
            *d_coarse_interp_to_nbr_fill->begin(fill_nabrs);

         hier::BoxContainer fill_boxes(fill_nabr);
         fill_boxes.intersectBoxes(hier::Box::refine(coarse_box,
                                      coarse_to_fine.getRatio()));
         const hier::BoxId& fill_id = fill_nabr.getBoxId();

         std::shared_ptr<hier::Patch> nbr_fill_patch(
            d_nbr_blk_fill_level->getPatch(fill_id));

         /*
          * Access the Transformation so that dst_fill_box can be computed
          * in terms of the fine patch's block's index space.
          */

         hier::BaseGridGeometry::ConstNeighborIterator nbr_itr =
            grid_geom->find(fine_box.getBlockId(), coarse_box.getBlockId()); 
         TBOX_ASSERT(nbr_itr != grid_geom->end(fine_box.getBlockId()));

         const hier::Transformation& transformation =
            (*nbr_itr).getTransformation(fine_level->getLevelNumber());

         hier::Box dst_fill_box(fill_nabr);
         transformation.transform(dst_fill_box);

         hier::IntVector fill_ghost_width(
            this == d_top_refine_schedule ?
            d_boundary_fill_ghost_width : d_max_stencil_width);

         fill_ghost_width.max(hier::IntVector::getOne(fine_box.getDim()));
         hier::Box fill_fine_box(fine_box);
         fill_fine_box.grow(fill_ghost_width);
         dst_fill_box *= fill_fine_box;

         overlaps.push_back(
            std::vector<std::shared_ptr<hier::BoxOverlap> >(0));
         d_nbr_blk_copy_overlaps.push_back(
            std::vector<std::shared_ptr<hier::BoxOverlap> >(0));
         std::vector<std::shared_ptr<hier::BoxOverlap> >& refine_overlaps(
            overlaps.back());
         std::vector<std::shared_ptr<hier::BoxOverlap> >& nbr_overlaps(
            d_nbr_blk_copy_overlaps.back());
         refine_overlaps.resize(num_equiv_classes);
         nbr_overlaps.resize(num_equiv_classes);

         hier::BoxContainer node_fill_boxes; //Unneeded here, so empty.

         /*
          * The overlaps computed here are refine_overlaps, to control
          * the interpolation of data onto nbr_fill_patch, and nbr_overlaps,
          * to control local copies from nbr_fill_patch to fine_patch.
          */
         for (int ne = 0; ne < num_equiv_classes; ++ne) {

            const RefineClasses::Data& rep_item =
               d_refine_classes->getClassRepresentative(ne);

            if (rep_item.d_oprefine) {

               const int scratch_id = rep_item.d_scratch;
               std::shared_ptr<hier::PatchDataFactory> fine_pdf(
                  fine_patch_descriptor->getPatchDataFactory(scratch_id));
               hier::Box scratch_space(nbr_fill_patch->getBox());

               refine_overlaps[ne] =
                  bg_fill_pattern->computeFillBoxesOverlap(
                     fill_boxes,
                     node_fill_boxes,
                     nbr_fill_patch->getBox(),
                     scratch_space,
                     *fine_pdf);

               nbr_overlaps[ne] =
                  rep_item.d_var_fill_pattern->calculateOverlap(
                     *fine_pdf->getBoxGeometry(fine_patch->getBox()),
                     *fine_pdf->getBoxGeometry(nbr_fill_patch->getBox()),
                     fine_patch->getBox(),
                     fill_nabr,
                     dst_fill_box,
                     true,
                     transformation);

            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Generate communication schedule routine creates transactions to move
 * data from interiors of the source space on the source level into the
 * specified fill box regions of the destination level.
 *
 * The resulting transactions will only fill the regions of intersection
 * between the fill_boxes and destination level boxes.  The remaining
 * box regions are returned in unfilled_boxes.
 *
 *********************************************************************
 */

void
RefineSchedule::generateCommunicationSchedule(
   std::shared_ptr<hier::BoxLevel>& unfilled_box_level,
   std::shared_ptr<hier::Connector>& dst_to_unfilled,
   std::shared_ptr<hier::BoxLevel>& unfilled_encon_box_level,
   std::shared_ptr<hier::Connector>& encon_to_unfilled_encon,
   const hier::Connector& dst_to_fill,
   const hier::BoxNeighborhoodCollection& dst_to_fill_on_src_proc,
   const bool use_time_interpolation,
   const bool create_transactions)
{
   t_gen_comm_sched->start();
   TBOX_ASSERT(d_dst_to_src);
   TBOX_ASSERT(d_dst_to_src->hasTranspose());

   std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   if (s_extra_debug) {
      if (d_dst_to_src->isFinalized()) {
         hier::Connector& src_to_dst = d_dst_to_src->getTranspose();
         d_dst_to_src->assertTransposeCorrectness(src_to_dst);
         src_to_dst.assertTransposeCorrectness(*d_dst_to_src);
      }
   }

   const hier::BoxLevel& dst_box_level = dst_to_fill.getBase();

   unfilled_box_level.reset(new hier::BoxLevel(
         dst_box_level.getRefinementRatio(),
         dst_box_level.getGridGeometry(),
         dst_box_level.getMPI()));

   dst_to_unfilled.reset(new hier::Connector(
         dst_box_level,
         *unfilled_box_level,
         dst_to_fill.getConnectorWidth()));

   if (grid_geometry->hasEnhancedConnectivity()) {
      unfilled_encon_box_level.reset(new hier::BoxLevel(
            dst_box_level.getRefinementRatio(),
            dst_box_level.getGridGeometry(),
            dst_box_level.getMPI()));

      encon_to_unfilled_encon.reset(new hier::Connector(
            *(d_encon_level->getBoxLevel()),
            *unfilled_encon_box_level,
            dst_to_fill.getConnectorWidth()));
   }

   if (create_transactions) {

      /*
       * Reorder d_dst_to_src's transpose's edge data to arrange neighbors by
       * the dst boxes, as required to match the transaction ordering
       * on the receiving processors.
       */
      FullNeighborhoodSet src_to_dst_edges_bydst;
      t_invert_edges->start();
      reorderNeighborhoodSetsByDstNodes(src_to_dst_edges_bydst);
      t_invert_edges->stop();

      t_construct_send_trans->start();

      /*
       * Construct transactions with local source and remote destination.
       */
      for (FullNeighborhoodSet::const_iterator
           ei = src_to_dst_edges_bydst.begin();
           ei != src_to_dst_edges_bydst.end(); ++ei) {

         /*
          * dst_box can be remote (by definition of FullNeighborhoodSet).
          * local_src_boxes are the local source boxes that
          * contribute data to dst_box.
          */
         const hier::Box& dst_box = ei->first;
         const hier::BoxContainer& local_src_boxes = ei->second;
         TBOX_ASSERT(!dst_box.isPeriodicImage());

         hier::BoxNeighborhoodCollection::ConstIterator dst_fill_iter =
            dst_to_fill_on_src_proc.find(dst_box.getBoxId());
         if (dst_fill_iter == dst_to_fill_on_src_proc.end()) {
            /*
             * Missing fill boxes should indicate that the dst box
             * has no fill box.  One way this is possible is for
             * d_dst_level_fill_pattern to be of type PatchLevelBorderFillPattern
             * and for dst_box to be away from level borders.
             */
            continue;
         }

         int num_nbrs = dst_to_fill_on_src_proc.numNeighbors(dst_fill_iter);
         hier::BoxNeighborhoodCollection::ConstNeighborIterator nbrs_begin =
            dst_to_fill_on_src_proc.begin(dst_fill_iter);
         hier::BoxNeighborhoodCollection::ConstNeighborIterator nbrs_end =
            dst_to_fill_on_src_proc.end(dst_fill_iter);
         for (hier::BoxContainer::const_iterator ni = local_src_boxes.begin();
              ni != local_src_boxes.end(); ++ni) {
            const hier::Box& src_box = *ni;

            if (src_box.getOwnerRank() != dst_box.getOwnerRank()) {

               constructScheduleTransactions(
                  num_nbrs,
                  nbrs_begin,
                  nbrs_end,
                  dst_box,
                  src_box,
                  use_time_interpolation);

            }
         }

      } // end send transactions loop

      t_construct_send_trans->stop();

   } // if create_transactions

   hier::LocalId last_unfilled_local_id(-1);

   t_construct_recv_trans->start();
   for (hier::Connector::ConstNeighborhoodIterator cf = dst_to_fill.begin();
        cf != dst_to_fill.end(); ++cf) {

      const hier::BoxId& dst_box_id(*cf);
      const hier::Box& dst_box = *dst_box_level.getBox(dst_box_id);
      const hier::BlockId& dst_block_id = dst_box.getBlockId();

      hier::BoxContainer fill_boxes_list;
      for (hier::Connector::ConstNeighborIterator bi = dst_to_fill.begin(cf);
           bi != dst_to_fill.end(cf); ++bi) {
         fill_boxes_list.pushBack(*bi);
      }

      /*
       * If the destination box touches enhanced connectivity,
       * encon_fill_boxes will hold the portion of fill_boxes_list that lies
       * in the enhanced connectivity ghost region of the destination.
       * Otherwise, encon_fill_boxes will be empty.
       */
      hier::BoxContainer encon_fill_boxes;
      if (grid_geometry->hasEnhancedConnectivity()) {
         findEnconFillBoxes(encon_fill_boxes,
            fill_boxes_list,
            dst_block_id);
      }

      /*
       * unfilled_boxes_for_dst starts out containing the fill boxes
       * for the current dst_box.  As transactions are created,
       * the space covered by those transactions will be removed from this
       * list, and whatever is left cannot be filled from the source level.
       *
       * The boxes in encon_fill_boxes, if any, are removed from
       * unfilled_boxes_for_dst, because unfilled boxes at enhanced
       * connectivity are handled separately at a later step.
       */
      hier::BoxContainer unfilled_boxes_for_dst(fill_boxes_list);
      unfilled_boxes_for_dst.removeIntersections(encon_fill_boxes);

      if (create_transactions) {

         hier::Connector::ConstNeighborhoodIterator dst_to_src_iter =
            d_dst_to_src->findLocal(dst_box_id);

         if (dst_to_src_iter != d_dst_to_src->end()) {

            int num_nbrs = dst_to_fill.numLocalNeighbors(*cf);
            hier::Connector::ConstNeighborIterator nbrs_begin =
               dst_to_fill.begin(cf);
            hier::Connector::ConstNeighborIterator nbrs_end =
               dst_to_fill.end(cf);
            for (hier::Connector::ConstNeighborIterator
                 na = d_dst_to_src->begin(dst_to_src_iter);
                 na != d_dst_to_src->end(dst_to_src_iter); ++na) {

               const hier::Box& src_box = *na;
               const hier::BlockId& src_block_id = src_box.getBlockId();

               /*
                * Remove the source box from the unfilled boxes list.  If
                * necessary, the source box is transformed to the destination
                * coordinate system.
                *
                * The removal of the source box is skipped if the source
                * and destination are enhanced connectivity (singularity)
                * neighbors, since the handling of unfilled boxes in this case
                * is a separate step.
                */
               if (src_block_id != dst_block_id) {

                  if (!grid_geometry->areSingularityNeighbors(dst_block_id,
                         src_block_id)) {

                     hier::Box transformed_src_box(src_box);

                     grid_geometry->transformBox(
                        transformed_src_box,
                        d_dst_level->getLevelNumber(),
                        dst_block_id,
                        src_block_id);

                     unfilled_boxes_for_dst.removeIntersections(
                        transformed_src_box);

                  }
               } else {
                  unfilled_boxes_for_dst.removeIntersections(src_box);
               }
               constructScheduleTransactions(
                  num_nbrs,
                  nbrs_begin,
                  nbrs_end,
                  dst_box,
                  src_box,
                  use_time_interpolation);

            }
         }
      }

      if (grid_geometry->hasEnhancedConnectivity() &&
          !unfilled_boxes_for_dst.empty()) {
         hier::BoxContainer fixed_unfilled_boxes(unfilled_boxes_for_dst);
         hier::BoxContainer dst_block_domain;
         grid_geometry->computePhysicalDomain(
            dst_block_domain,
            d_dst_level->getRatioToLevelZero(),
            dst_box.getBlockId());
         fixed_unfilled_boxes.intersectBoxes(dst_block_domain);

         for (hier::BoxContainer::iterator bi = unfilled_boxes_for_dst.begin();
              bi != unfilled_boxes_for_dst.end(); ++bi) {

            for (hier::BaseGridGeometry::ConstNeighborIterator ni =
                    grid_geometry->begin(dst_box.getBlockId());
                 ni != grid_geometry->end(dst_box.getBlockId()); ++ni) {

               hier::BoxContainer transformed_domain(
                  (*ni).getTransformedDomain());

               transformed_domain.refine(d_dst_level->getRatioToLevelZero());

               transformed_domain.intersectBoxes(*bi);

               if (!transformed_domain.empty()) {
                  fixed_unfilled_boxes.spliceBack(transformed_domain);
               }
            }
         }
         unfilled_boxes_for_dst.clear();
         unfilled_boxes_for_dst.spliceBack(fixed_unfilled_boxes);
      }

      /*
       * Add mapping information to unfilled boxes and add them to
       * containers for the level.
       */
      if (!unfilled_boxes_for_dst.empty()) {

         unfilled_boxes_for_dst.coalesce();

         hier::Connector::NeighborhoodIterator base_box_itr =
            dst_to_unfilled->makeEmptyLocalNeighborhood(dst_box_id);
         for (hier::BoxContainer::iterator bi = unfilled_boxes_for_dst.begin();
              bi != unfilled_boxes_for_dst.end(); ++bi) {

            hier::Box unfilled_box(*bi,
                                   ++last_unfilled_local_id,
                                   dst_box.getOwnerRank());
            TBOX_ASSERT(unfilled_box.getBlockId() == dst_block_id);

            unfilled_box_level->addBoxWithoutUpdate(unfilled_box);
            dst_to_unfilled->insertLocalNeighbor(
               unfilled_box,
               base_box_itr);

         }
      }

      /*
       * Call private method to handle unfilled boxes where the destination
       * box touches enhanced connectivity.
       */
      if (!encon_fill_boxes.empty()) {
         findEnconUnfilledBoxes(
            unfilled_encon_box_level,
            encon_to_unfilled_encon,
            last_unfilled_local_id,
            dst_box,
            encon_fill_boxes);
      }

   } // End receive/copy transactions loop
   unfilled_box_level->finalize();
   if (grid_geometry->hasEnhancedConnectivity()) {
      unfilled_encon_box_level->finalize();
   }
   t_construct_recv_trans->stop();

   t_gen_comm_sched->stop();
}

/*
 ***********************************************************************
 * Given a BlockId and a list of fill boxes, populate encon_fill_boxes
 * with the portion of the fill boxes that lie across an enhanced
 * connectivity boundary from the destination block.
 ***********************************************************************
 */
void
RefineSchedule::findEnconFillBoxes(
   hier::BoxContainer& encon_fill_boxes,
   const hier::BoxContainer& fill_boxes_list,
   const hier::BlockId& dst_block_id)
{
   TBOX_ASSERT(encon_fill_boxes.empty());

   std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   hier::BoxContainer encon_neighbor_list;
   for (hier::BaseGridGeometry::ConstNeighborIterator ni =
           grid_geometry->begin(dst_block_id);
        ni != grid_geometry->end(dst_block_id); ++ni) {

      const hier::BaseGridGeometry::Neighbor& nbr = *ni;
      if (nbr.isSingularity()) {
         hier::BoxContainer transformed_domain(
            nbr.getTransformedDomain());
         encon_fill_boxes.spliceFront(transformed_domain);
      }

   }

   encon_fill_boxes.refine(d_dst_level->getRatioToLevelZero());

   encon_fill_boxes.intersectBoxes(fill_boxes_list);

   encon_fill_boxes.coalesce();
}

/*
 ***********************************************************************
 * Given a list of fill boxes for a destination box at enhanced
 * connectivity, determine which fill boxes cannot be filled from the
 * source level and add those unfilled boxes to the appropriate containers.
 ***********************************************************************
 */
void
RefineSchedule::findEnconUnfilledBoxes(
   const std::shared_ptr<hier::BoxLevel>& unfilled_encon_box_level,
   const std::shared_ptr<hier::Connector>& encon_to_unfilled_encon,
   hier::LocalId& last_unfilled_local_id,
   const hier::Box& dst_box,
   const hier::BoxContainer& encon_fill_boxes)
{
   TBOX_ASSERT(d_dst_to_src);

   std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());

   const hier::BoxId& dst_box_id = dst_box.getBoxId();
   const hier::BlockId& dst_block_id = dst_box.getBlockId();

   /*
    * map container will hold unfilled boxes for each block that is
    * an enhanced connectivity neighbor of the destination box.
    *
    * To start, each entry in the map will have the intersection of
    * encon_fill_boxes and the neighboring block's domain, represented in
    * the destination coordinate system.
    */
   std::map<hier::BlockId, hier::BoxContainer> unfilled_encon_nbr_boxes;

   for (hier::BaseGridGeometry::ConstNeighborIterator ni =
           grid_geometry->begin(dst_block_id);
        ni != grid_geometry->end(dst_block_id); ++ni) {

      const hier::BaseGridGeometry::Neighbor& nbr = *ni;
      if (nbr.isSingularity()) {
         const hier::BlockId nbr_block_id(nbr.getBlockId());

         hier::BoxContainer neighbor_boxes(nbr.getTransformedDomain());
         neighbor_boxes.refine(d_dst_level->getRatioToLevelZero());
         neighbor_boxes.intersectBoxes(encon_fill_boxes);
         unfilled_encon_nbr_boxes[nbr_block_id].spliceFront(neighbor_boxes);
      }
   }

   if (d_src_level) {
      /*
       * If there are overlapping source boxes found in d_dst_to_src,
       * and if those source boxes lie across enhanced connectivity from
       * the destination box, then we remove the source box from the
       * source block's entry in the unfilled_encon_nbr_boxes map container.
       */
      hier::Connector::ConstNeighborhoodIterator dst_to_src_iter =
         d_dst_to_src->findLocal(dst_box_id);

      if (dst_to_src_iter != d_dst_to_src->end()) {

         /*
          * If at enhanced connectivity, remove source box from container of
          * unfilled boxes
          */
         for (hier::Connector::ConstNeighborIterator na = d_dst_to_src->begin(dst_to_src_iter);
              na != d_dst_to_src->end(dst_to_src_iter); ++na) {

            const hier::Box& src_box = *na;
            const hier::BlockId& src_block_id = src_box.getBlockId();

            if (src_block_id != dst_block_id) {

               if (grid_geometry->areSingularityNeighbors(dst_block_id,
                      src_block_id)) {

                  hier::Box transformed_src_box(src_box);
                  grid_geometry->transformBox(transformed_src_box,
                     d_dst_level->getLevelNumber(),
                     dst_block_id,
                     src_block_id);

                  unfilled_encon_nbr_boxes[src_block_id].removeIntersections(
                     transformed_src_box);

               }
            }
         }
      }
   }

   /*
    * If any unfilled boxes remain for any block in unfilled_encon_nbr_boxes,
    * they need to be added to the output containers.
    */
   if (unfilled_encon_nbr_boxes.size() > 0) {

      hier::Connector::ConstNeighborhoodIterator find_encon_nabrs =
         d_dst_to_encon->findLocal(dst_box_id);

      if (find_encon_nabrs != d_dst_to_encon->end()) {

         for (hier::Connector::ConstNeighborIterator de_iter =
                 d_dst_to_encon->begin(find_encon_nabrs);
              de_iter != d_dst_to_encon->end(find_encon_nabrs); ++de_iter) {

            const hier::BoxId& encon_box_id = de_iter->getBoxId();
            const hier::BlockId& nbr_block_id = de_iter->getBlockId();

            const hier::BoxContainer& unfilled_boxes =
               unfilled_encon_nbr_boxes[nbr_block_id];

            if (!unfilled_boxes.empty()) {

               /*
                * The unfilled boxes are, at this point, represented in
                * the destination coordinates.  Here they are transformed
                * into the neighboring block's coordinates before being
                * added to the output containers.
                */

               hier::Connector::NeighborhoodIterator base_box_itr =
                  encon_to_unfilled_encon->makeEmptyLocalNeighborhood(
                     encon_box_id);
               for (hier::BoxContainer::const_iterator bi = unfilled_boxes.begin();
                    bi != unfilled_boxes.end(); ++bi) {

                  hier::Box unfilled_box(*bi);
                  grid_geometry->transformBox(
                     unfilled_box,
                     d_dst_level->getLevelNumber(),
                     nbr_block_id,
                     dst_block_id);

                  hier::Box unfilled_encon_box(unfilled_box,
                                               ++last_unfilled_local_id,
                                               dst_box.getOwnerRank());

                  TBOX_ASSERT(unfilled_encon_box.getBlockId() == nbr_block_id);
                  unfilled_encon_box_level->addBox(unfilled_encon_box);

                  encon_to_unfilled_encon->insertLocalNeighbor(
                     unfilled_encon_box,
                     base_box_itr);

               }
            }
         }
      }
   }
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
 * 2. It shifts periodic image dst Boxes back to the zero-shift position,
 * and applies a similar shift to src Boxes so that the overlap is
 * unchanged.  The constructScheduleTransactions method requires all
 * shifts to be absorbed in the src Box.
 ***********************************************************************
 */
void
RefineSchedule::reorderNeighborhoodSetsByDstNodes(
   FullNeighborhoodSet& full_inverted_edges) const
{
   TBOX_ASSERT(d_dst_to_src);
   TBOX_ASSERT(d_dst_to_src->hasTranspose());

   hier::Connector& src_to_dst = d_dst_to_src->getTranspose();

   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::BoxLevel& src_box_level = src_to_dst.getBase();
   const hier::IntVector& src_ratio = src_box_level.getRefinementRatio();
   const hier::IntVector& dst_ratio = src_to_dst.getHead().getRefinementRatio();

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
      const hier::Box& src_box = *src_box_level.getBoxStrict(*ci);
      for (hier::Connector::ConstNeighborIterator na = src_to_dst.begin(ci);
           na != src_to_dst.end(ci); ++na) {
         const hier::Box& nabr = *na;
         if (nabr.isPeriodicImage()) {
            shifted_box.initialize(
               src_box,
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
            full_inverted_edges[nabr].insert(src_box);
         }
      }
   }
}

/*
 ****************************************************************
 * Set the fill boxes and the dst--->fill Connector.
 * Initialize this same information on the src processes as well.
 *
 * Compute the fill boxes of the dst level according to the
 * d_dst_level_fill_pattern.  Initialize the fill BoxLevel from
 * the those fill boxes.
 *
 * Set the edges between dst and fill box_levels (trivial,
 * since fill BoxLevel is related to dst BoxLevel).  Set
 * edges between src and fill box_levels (by copying the edges
 * between src and dst.
 *
 * If d_dst_level_fill_pattern can compute the relationship between
 * dst boxes and fill boxes on the src processes, have it do that.
 * Otherwise, perform communication to get that information onto the
 * src processes.
 *
 * The Connectors dst_to_src and src_to_src (in the arguments)
 * are required only for fill_patterns PatchLevelBorderFillPattern and
 * and PatchLevelBorderAndInteriorFillPattern
 ****************************************************************
 */

void
RefineSchedule::setDefaultFillBoxLevel(
   std::shared_ptr<hier::BoxLevel>& fill_box_level,
   std::shared_ptr<hier::Connector>& dst_to_fill,
   hier::BoxNeighborhoodCollection& dst_to_fill_on_src_proc)
{
   TBOX_ASSERT(d_dst_to_src);
   TBOX_ASSERT(d_dst_to_src->hasTranspose());

   /*
    * If the destination variable of any registered communication
    * requires that coarse data take precedence over fine data at
    * coarse-fine boundaries, make sure that the ghost cell width is
    * at least one in each direction.  This is necessary so that
    * generateCommunicationSchedule() determines that boundary data for
    * the level needs to be transferred from the next coarser level
    * (i.e. the interiors of the fill boxes overlap the next coarser
    * level in a nontrivial way).
    */

   const tbox::Dimension& dim(d_dst_level->getDim());
   const hier::BoxLevel& dst_box_level(*d_dst_level->getBoxLevel());
   const hier::IntVector& fill_ghost_width(
      this == d_top_refine_schedule ?
      d_boundary_fill_ghost_width : d_max_stencil_width);

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim, dst_box_level, fill_ghost_width);

   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   bool need_nontrivial_ghosts = false;
   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      if (!(d_refine_items[iri]->d_fine_bdry_reps_var)) {
         need_nontrivial_ghosts = true;
      }
   }

   hier::IntVector fill_gcw = fill_ghost_width;
   if (need_nontrivial_ghosts) {
      fill_gcw.max(constant_one_intvector);
   }

   // New data computed here:

   /*
    * d_max_fill_boxes is the max number of fill boxes
    * for any dst box.  d_max_fill_boxes is nominally 1.
    * If more is required (based on fill_pattern),
    * it is incremented below.
    */
   d_max_fill_boxes = 1;

   /*
    * Generate fill boxes (grown dst boxes) and
    * get edges between dst and fill box_levels.
    * These are all local edges by definition.
    * Note that since fill boxes are simply grown
    * versions of dst_box (at most), edges between
    * fill and src box_levels would have gcw of zero.
    *
    * Technically speaking, edges between dst and
    * fill are not complete.  Once the dst box is
    * grown to make the fill box, the fill box may
    * intersect other boxes in the dst box_level.  We
    * ignore all these overlaps because they are not
    * needed by the algorithm.
    */

   /*
    * Note that dst_to_fill is NOT complete.
    * We care only about the intersection between a dst_box
    * and the fill_box it created.  In fact, a dst_box
    * may also intersect the fill_box of other nearby
    * destination boxes.
    */

   d_dst_level_fill_pattern->computeFillBoxesAndNeighborhoodSets(
      fill_box_level,
      dst_to_fill,
      dst_box_level,
      fill_gcw,
      d_data_on_patch_border_flag);

   d_max_fill_boxes = tbox::MathUtilities<int>::Max(
         d_max_fill_boxes,
         d_dst_level_fill_pattern->getMaxFillBoxes());

   if (d_src_level) {
      if (d_dst_level_fill_pattern->needsToCommunicateDestinationFillBoxes()) {

         /*
          * This part assumes src-dst ratio is one.
          * Should be modified if the assumption does not hold.
          */
         TBOX_ASSERT(d_dst_to_src->getRatio() == 1);

         /*
          * For these fill_pattern, the src owner could not compute fill boxes
          * for all its dst neighbors using local data, so the dst owners
          * must send this information.
          */
         communicateFillBoxes(dst_to_fill_on_src_proc, *dst_to_fill);

      } else {
         d_dst_level_fill_pattern->computeDestinationFillBoxesOnSourceProc(
            dst_to_fill_on_src_proc,
            dst_box_level,
            d_dst_to_src->getTranspose(),
            fill_gcw);
      }
   }

   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      createEnconLevel(fill_gcw);
   }
}

/*
 ****************************************************************
 * Setup the level representing ghost regions of destination patches
 * that touch enhanced connectivity.
 *****************************************************************
 */
void
RefineSchedule::createEnconLevel(const hier::IntVector& fill_gcw)
{
   const tbox::Dimension& dim = fill_gcw.getDim();

   std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
      d_dst_level->getGridGeometry());
   const size_t num_blocks = grid_geometry->getNumberBlocks();

   /*
    * Create encon_box_level and associated Connectors.
    *
    * Where destination patches have ghost regions across enhanced connectivity
    * boundaries, data communicated by this schedule will not be written
    * directly into those ghost regions, but rather into patches on
    * d_encon_level.  This level, once filled with data, will be provided
    * to SingularityPatchStrategy's fillSingularityBoundaryConditions.
    */
   std::shared_ptr<hier::BoxLevel> encon_box_level(
      std::make_shared<hier::BoxLevel>(d_dst_level->getRatioToLevelZero(),
                                         grid_geometry));

   d_dst_to_encon.reset(new hier::Connector(dim));
   d_dst_to_encon->setBase(*(d_dst_level->getBoxLevel()));

   hier::IntVector encon_gcw(
      hier::IntVector::max(fill_gcw, hier::IntVector::getOne(dim)));

   if (num_blocks > 1) {

      hier::LocalId encon_local_id(0);

      /*
       * Loop over blocks to find destination patches on each block
       * that touch enhanced connectivity.
       */
      for (hier::BlockId::block_t bn = 0; bn < num_blocks; ++bn) {

         hier::BlockId block_id(bn);

         /*
          * Test to see if there are any local destination boxes on this
          * block.  Move on to next block if not.
          */
         const hier::BoxContainer& level_boxes(
            d_dst_level->getBoxLevel()->getBoxes());
         hier::BoxContainerSingleBlockIterator dst_test_iter(
            level_boxes.begin(block_id));

         if (dst_test_iter == level_boxes.end(block_id)) {
            continue;
         }

         /*
          * An empty singularity box list would mean the block touches no
          * enhanced connectivy boundaries, so no need to go further with
          * this block.
          */
         const hier::BoxContainer& sing_boxes =
            grid_geometry->getSingularityBoxContainer(block_id);

         if (!sing_boxes.empty()) {

            /*
             * Loop over neighboring blocks and find the ones that are
             * singularity neighbors.
             */
            for (hier::BaseGridGeometry::ConstNeighborIterator ni =
                    grid_geometry->begin(block_id);
                 ni != grid_geometry->end(block_id); ++ni) {

               const hier::BaseGridGeometry::Neighbor& nbr = *ni;
               if (nbr.isSingularity()) {

                  const hier::BlockId& nbr_id = nbr.getBlockId();

                  /*
                   * Get the transformation from neighbor block to dst
                   * block, and get a representation of the neighbor block
                   * domain in coordinate system of dst block.
                   */
                  hier::Transformation::RotationIdentifier rotation =
                     nbr.getRotationIdentifier();
                  const hier::IntVector& offset =
                     nbr.getShift(d_dst_level->getLevelNumber());

                  hier::Transformation transformation(rotation, offset,
                                                      nbr_id, block_id);

                  hier::BoxContainer trans_neighbor_list;
                  grid_geometry->getTransformedBlock(trans_neighbor_list,
                     block_id,
                     nbr_id);
                  trans_neighbor_list.refine(
                     d_dst_level->getRatioToLevelZero());

                  /*
                   * Loop over dst boxes for this block and
                   * determine if they touch the current neighbor block
                   * at enhanced connectivity.
                   */
                  hier::BoxContainerSingleBlockIterator dst_local_iter(
                     level_boxes.begin(block_id));

                  for ( ; dst_local_iter != level_boxes.end(block_id);
                        ++dst_local_iter) {

                     const hier::BoxId& box_id = dst_local_iter->getBoxId();

                     std::shared_ptr<hier::Patch> patch(
                        d_dst_level->getPatch(box_id));
                     std::shared_ptr<hier::PatchGeometry> pgeom(
                        patch->getPatchGeometry());

                     if (pgeom->getTouchesRegularBoundary()) {

                        /*
                         * Grow the patch box and intersect with the
                         * representation of the neighboring block.
                         * A non-empty intersection indicates there is
                         * a ghost region across an enhanced connectivity
                         * boundary.
                         */
                        hier::BoxContainer encon_test_list(patch->getBox());
                        encon_test_list.grow(encon_gcw);
                        encon_test_list.intersectBoxes(trans_neighbor_list);

                        if (!encon_test_list.empty()) {

                           encon_test_list.coalesce();
                           TBOX_ASSERT(encon_test_list.size() == 1);

                           /*
                            * Transform the boxes representing the ghost
                            * region back to the neighbor block's
                            * coordinate system, and create a Box
                            * to be added to encon_box_level.
                            */

                           hier::Connector::NeighborhoodIterator base_box_itr =
                              d_dst_to_encon->makeEmptyLocalNeighborhood(
                                 box_id);
                           for (hier::BoxContainer::iterator bi = encon_test_list.begin();
                                bi != encon_test_list.end(); ++bi) {
                              hier::Box encon_box(*bi);

                              transformation.inverseTransform(encon_box);

                              /*
                               * If a Box at this location on the
                               * same neighbor block and on the same processor
                               * already exists in encon_box_level,
                               * do not create another.
                               */
                              hier::Box actual_encon_box(dim);
                              if (!encon_box_level->getSpatiallyEqualBox(
                                     encon_box, nbr_id, actual_encon_box)) {
                                 actual_encon_box = hier::Box(encon_box,
                                       ++encon_local_id,
                                       box_id.getOwnerRank());
                                 TBOX_ASSERT(actual_encon_box.getBlockId() ==
                                    nbr_id);
                                 encon_box_level->addBoxWithoutUpdate(
                                    actual_encon_box);
                              }

                              /*
                               * Add to the neighborhood set for the
                               * d_dst_to_encon connector.
                               */
                              d_dst_to_encon->insertLocalNeighbor(
                                 actual_encon_box,
                                 base_box_itr);

                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   /*
    * Finalize encon_box_level create d_encon_level and more associated
    * Connectors.
    */
   encon_box_level->finalize();

   d_encon_level.reset(new hier::PatchLevel(encon_box_level,
         grid_geometry,
         d_dst_level->getPatchDescriptor(),
         std::shared_ptr<hier::PatchFactory>(),
         true));

   d_encon_level->setLevelNumber(d_dst_level->getLevelNumber());

   d_dst_to_encon->setHead(*(d_encon_level->getBoxLevel()));
   d_dst_to_encon->setWidth(encon_gcw, true);

   if (d_src_level) {
      const hier::Connector& dst_to_src =
         d_dst_level->findConnectorWithTranspose(*d_src_level,
            encon_gcw,
            encon_gcw,
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      // d_dst_to_encon only needs its transpose set temporarily as the
      // transpose is only used in this call to bridge.  We do not want to
      // store d_dst_to_encon's transpose after this point which is why it is
      // deleted a few lines later.
      hier::Connector* encon_to_dst = d_dst_to_encon->createLocalTranspose();
      d_dst_to_encon->setTranspose(encon_to_dst, false);

      hier::OverlapConnectorAlgorithm oca;
      oca.setTimerPrefix("xfer::RefineSchedule_build");

      oca.bridge(
         d_encon_to_src,
         *encon_to_dst,
         dst_to_src,
         hier::IntVector::getZero(dim),
         true);

      d_dst_to_encon->setTranspose(0, false);
      delete encon_to_dst;
   }
}

/*
 **************************************************************************
 * Calculate the maximum ghost cell width of all destination patch data
 * components.
 *
 * It is possible for dst_to_fill_on_src_proc to omit a visible dst box
 * if the dst box has no fill boxes.
 **************************************************************************
 */
void
RefineSchedule::communicateFillBoxes(
   hier::BoxNeighborhoodCollection& dst_to_fill_on_src_proc,
   const hier::Connector& dst_to_fill)
{
   TBOX_ASSERT(d_dst_to_src);
   TBOX_ASSERT(d_dst_to_src->hasTranspose());

   hier::Connector& src_to_dst = d_dst_to_src->getTranspose();

   const tbox::Dimension& dim(d_dst_level->getDim());

   int rank = dst_to_fill.getBase().getMPI().getRank();

   std::set<int> src_owners;
   d_dst_to_src->getLocalOwners(src_owners);
   src_owners.erase(rank);

   std::set<int> dst_owners;
   src_to_dst.getLocalOwners(dst_owners);
   dst_owners.erase(rank);

   std::map<int, std::vector<int> > send_mesgs;
   std::map<int, std::vector<int> > recv_mesgs;

   int num_incoming_comms = static_cast<int>(dst_owners.size());
   int num_comms = static_cast<int>(src_owners.size()) + num_incoming_comms;
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* comms = new tbox::AsyncCommPeer<int>[num_comms];

   // Set up all comms.
   // We want to post the first receive for the processor with rank immediatly
   // less that this processor's.  Set mesg_number to the index of the comm
   // whose processor's receive we want to post first.
   int mesg_number = 0;
   int counter = 0;
   for (std::set<int>::const_iterator di = dst_owners.begin();
        di != dst_owners.end(); ++di) {
      comms[counter].initialize(&comm_stage);
      comms[counter].setPeerRank(*di);
      // Reuse communicator.  Assuming it has no outstanding messages!
      comms[counter].setMPI(src_to_dst.getMPI());
      comms[counter].setMPITag(0, 1);
      if (*di < rank) {
         ++mesg_number;
      }
      ++counter;
   }
   mesg_number = mesg_number > 0 ? mesg_number - 1 : num_incoming_comms - 1;

   // Set iterator si corresponds to comms[mesg_number].
   std::set<int>::const_iterator si = num_comms == 0 ?
      dst_owners.end() :
      dst_owners.find(comms[mesg_number].getPeerRank());

   /*
    * Post receives for fill boxes from dst owners.
    */
   for (counter = 0; counter < num_incoming_comms;
        ++counter, --si, --mesg_number) {
      comms[mesg_number].beginRecv();
      if (comms[mesg_number].isDone()) {
         comms[mesg_number].pushToCompletionQueue();
      }

      if (si == dst_owners.begin()) {
         // Continue loop at the opposite end.
         si = dst_owners.end();
         mesg_number = num_incoming_comms;
      }
   }

   /*
    * Pack fill boxes and send messages to src owners.
    */
   // Pack messages.
   std::vector<int> tmp_mesg;
   hier::BoxContainer tmp_fill_boxes;
   for (hier::Connector::ConstNeighborhoodIterator ei = dst_to_fill.begin();
        ei != dst_to_fill.end(); ++ei) {
      const hier::BoxId& dst_box_id = *ei;
      /*
       * Pack dst_box_id's fill box info into tmp_mesg.
       * - dst_box_id's LocalId
       * - number of fill neighbors
       * - fill neighbors (could just send box and save 2 ints)
       */
      tmp_mesg.clear();
      tmp_mesg.reserve(3 + dst_to_fill.numLocalNeighbors(dst_box_id)
         * hier::Box::commBufferSize(dim));
      tmp_mesg.insert(tmp_mesg.end(), 3, 0);
      tmp_mesg[0] = dst_box_id.getLocalId().getValue();
      tmp_mesg[1] = -1;
      tmp_mesg[2] = static_cast<int>(dst_to_fill.numLocalNeighbors(dst_box_id));
      tmp_fill_boxes.clear();
      for (hier::Connector::ConstNeighborIterator na = dst_to_fill.begin(ei);
           na != dst_to_fill.end(ei); ++na) {
         tmp_mesg.insert(tmp_mesg.end(), hier::Box::commBufferSize(dim), 0);
         na->putToIntBuffer(&tmp_mesg[tmp_mesg.size()
                                      - hier::Box::commBufferSize(dim)]);
         tmp_fill_boxes.insert(*na);
      }
      // Append tmp_mesg to buffers for sending to src owners.
      hier::Connector::ConstNeighborhoodIterator di =
         d_dst_to_src->findLocal(dst_box_id);
      if (di != d_dst_to_src->end()) {
         std::set<int> tmp_owners;
         d_dst_to_src->getLocalOwners(di, tmp_owners);
         for (std::set<int>::const_iterator so = tmp_owners.begin();
              so != tmp_owners.end(); ++so) {
            const int& src_owner = *so;
            if (src_owner == dst_box_id.getOwnerRank()) {
               dst_to_fill_on_src_proc.insert(
                  dst_box_id,
                  tmp_fill_boxes);
            } else {
               std::vector<int>& send_mesg = send_mesgs[src_owner];
               send_mesg.insert(send_mesg.end(),
                  tmp_mesg.begin(),
                  tmp_mesg.end());
            }
         }
      }
   }

   // Send messages.
   // We want the first send to be to the processor whose rank is immediately
   // greater that the sender's.  Set si to the processor whose send we want to
   // post first.
   si = src_owners.lower_bound(rank + 1);
   for (mesg_number = num_incoming_comms;
        mesg_number < num_comms;
        ++mesg_number, ++si) {
      if (si == src_owners.end()) {
         si = src_owners.begin();
      }
      comms[mesg_number].initialize(&comm_stage);
      comms[mesg_number].setPeerRank(*si);
      // Reuse communicator.  Assuming it has no outstanding messages!
      comms[mesg_number].setMPI(src_to_dst.getMPI());
      comms[mesg_number].setMPITag(0, 1);
      comms[mesg_number].beginSend(&send_mesgs[*si][0],
         send_mesgs[*si].size());
      if (comms[mesg_number].isDone()) {
         comms[mesg_number].pushToCompletionQueue();
      }
   }

   /*
    * Complete communication and unpack messages.
    */
   while (comm_stage.hasCompletedMembers() || comm_stage.advanceSome()) {
      tbox::AsyncCommPeer<int>* peer =
         CPP_CAST<tbox::AsyncCommPeer<int> *>(comm_stage.popCompletionQueue());
      TBOX_ASSERT(peer != 0);
      if (peer < comms + src_owners.size()) {
         // This is a receive.  Unpack it.  (Otherwise, ignore send completion.)
         const int* ptr = peer->getRecvData();
         while (ptr != peer->getRecvData() + peer->getRecvSize()) {
            const hier::BoxId distributed_id(hier::LocalId(ptr[0]),
                                             peer->getPeerRank());
            const unsigned int num_fill_boxes = ptr[2];
            ptr += 3;
            d_max_fill_boxes = tbox::MathUtilities<int>::Max(
                  d_max_fill_boxes,
                  num_fill_boxes);
            hier::BoxNeighborhoodCollection::Iterator fill_boxes_iter =
               dst_to_fill_on_src_proc.insert(distributed_id).first;
            for (size_t ii = 0; ii < num_fill_boxes; ++ii) {
               hier::Box tmp_dst_box(dim);
               tmp_dst_box.getFromIntBuffer(ptr);
               ptr += hier::Box::commBufferSize(dim);
               dst_to_fill_on_src_proc.insert(
                  fill_boxes_iter,
                  tmp_dst_box);
            }
         }
         TBOX_ASSERT(ptr == peer->getRecvData() + peer->getRecvSize());
      }
   }

   delete[] comms;
}

/*
 **************************************************************************
 * Get whether there is data living on patch borders.
 **************************************************************************
 */

bool
RefineSchedule::getDataOnPatchBorderFlag() const
{
   bool rval = false;
   std::shared_ptr<hier::PatchDescriptor> pd(
      d_dst_level->getPatchDescriptor());

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      const int dst_id = d_refine_items[iri]->d_dst;
      const hier::PatchDataFactory& pdf = *pd->getPatchDataFactory(dst_id);
      if (pdf.dataLivesOnPatchBorder()) {
         rval = true;
         break;
      }
   }

   return rval;
}

/*
 **************************************************************************
 * Calculate the maximum ghost cell width of all destination patch data
 * components for the purpose of determining overlap criterion.
 *
 * When data lives on the border, increment width to create an overlap
 * in the cell index space.
 **************************************************************************
 */

hier::IntVector
RefineSchedule::getMinConnectorWidth() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector width(dim, 0);
   std::shared_ptr<hier::PatchDescriptor> pd(
      d_dst_level->getPatchDescriptor());

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {

      const int dst_id = d_refine_items[iri]->d_dst;
      const hier::PatchDataFactory& dst_pdf = *pd->getPatchDataFactory(dst_id);

      if (dst_pdf.dataLivesOnPatchBorder()) {
         width.max(dst_pdf.getGhostCellWidth() + 1);
      } else {
         width.max(dst_pdf.getGhostCellWidth());
      }

      const int scratch_id = d_refine_items[iri]->d_scratch;
      const hier::PatchDataFactory& scratch_pdf = *pd->getPatchDataFactory(scratch_id);

      if (scratch_pdf.dataLivesOnPatchBorder()) {
         width.max(scratch_pdf.getGhostCellWidth() + 1);
      } else {
         width.max(scratch_pdf.getGhostCellWidth());
      }

   }

   return width;
}

/*
 **************************************************************************
 * Calculate the maximum ghost cell width of all destination patch data
 * components.
 **************************************************************************
 */

hier::IntVector
RefineSchedule::getMaxDestinationGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   std::shared_ptr<hier::PatchDescriptor> pd(
      d_dst_level->getPatchDescriptor());

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      const int dst_id = d_refine_items[iri]->d_dst;
      gcw.max(pd->getPatchDataFactory(dst_id)->getGhostCellWidth());
   }

   return gcw;
}

/*
 **************************************************************************
 *
 * Calculate the maximum ghost cell width of all scratch patch data
 * components.
 *
 **************************************************************************
 */

hier::IntVector
RefineSchedule::getMaxScratchGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   std::shared_ptr<hier::PatchDescriptor> pd(
      d_dst_level->getPatchDescriptor());

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      const int scratch_id = d_refine_items[iri]->d_scratch;
      gcw.max(pd->getPatchDataFactory(scratch_id)->getGhostCellWidth());
   }

   return gcw;
}

/*
 **************************************************************************
 *
 * Calculate the maximum ghost cell width required for all stencils.
 *
 **************************************************************************
 */

hier::IntVector
RefineSchedule::getMaxStencilGhosts() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   hier::IntVector gcw(dim, 0);
   if (d_refine_patch_strategy) {
      gcw = d_refine_patch_strategy->getRefineOpStencilWidth(dim);
   }

   for (size_t iri = 0; iri < d_number_refine_items; ++iri) {
      if (d_refine_items[iri]->d_oprefine) {
         gcw.max(d_refine_items[iri]->d_oprefine->getStencilWidth(dim));
      }
   }

   return gcw;
}

/*
 *************************************************************************
 *
 * Private utility function that constructs schedule transactions that
 * move data from source patch on source level to destination patch
 * on destination level on regions defined by list of fil boxes.
 *
 *************************************************************************
 */

void
RefineSchedule::constructScheduleTransactions(
   int num_nbrs,
   hier::BoxNeighborhoodCollection::ConstNeighborIterator& nbrs_begin,
   hier::BoxNeighborhoodCollection::ConstNeighborIterator& nbrs_end,
   const hier::Box& dst_box,
   const hier::Box& src_box,
   bool use_time_interpolation)
{
   TBOX_ASSERT(d_dst_level);
   TBOX_ASSERT(d_src_level);
   TBOX_ASSERT(!dst_box.isPeriodicImage()); // src absorbs the shift, if any.

   const tbox::Dimension& dim(d_dst_level->getDim());
   const int my_rank = d_dst_level->getBoxLevel()->getMPI().getRank();

   const hier::IntVector& constant_zero_intvector(hier::IntVector::getZero(dim));
   const hier::IntVector& constant_one_intvector(hier::IntVector::getOne(dim));

   if (s_extra_debug) {
      tbox::plog << "RefineSchedule::constructScheduleTransactions: " << use_time_interpolation
                 << std::endl;
      tbox::plog << "  src: L" << d_src_level->getLevelNumber()
                 << "R" << d_src_level->getRatioToLevelZero()
                 << " / " << src_box << std::endl;
      tbox::plog << "  dst: L" << d_dst_level->getLevelNumber()
                 << "R" << d_dst_level->getRatioToLevelZero()
                 << " / " << dst_box << std::endl;
      tbox::plog << "  fill_boxes (" << num_nbrs << ")";
      for (hier::BoxNeighborhoodCollection::ConstNeighborIterator bi = nbrs_begin;
           bi != nbrs_end; ++bi) {
         tbox::plog << " " << *bi;
      }
      tbox::plog << std::endl;
   }

   /*
    * d_src_masks and d_overlaps exist only in this method, but are
    * class members instead of temporaries so that they don't
    * have to be reallocated every time this method is used.
    */
   int max_overlap_array_size = d_max_fill_boxes;

   if (d_src_masks.size() < max_overlap_array_size) {
      for (int i = d_src_masks.size();
           i < max_overlap_array_size; ++i) {
         d_src_masks.emplaceBack(dim);
      }
   }

   if (static_cast<int>(d_overlaps.size()) < max_overlap_array_size) {
      d_overlaps.clear();
      d_overlaps.resize(max_overlap_array_size);
   }

   std::shared_ptr<hier::PatchDescriptor> dst_patch_descriptor(
      d_dst_level->getPatchDescriptor());
   std::shared_ptr<hier::PatchDescriptor> src_patch_descriptor(
      d_src_level->getPatchDescriptor());

   const bool same_patch =
      (d_dst_level == d_src_level &&
       dst_box.getGlobalId() == src_box.getGlobalId());
   const bool same_patch_no_shift =
      (same_patch && !src_box.isPeriodicImage());

   const int num_equiv_classes =
      d_refine_classes->getNumberOfEquivalenceClasses();

   const hier::PeriodicShiftCatalog& shift_catalog =
      d_dst_level->getGridGeometry()->getPeriodicShiftCatalog();

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
      const hier::IntVector& src_ratio =
         d_src_level->getRatioToLevelZero();
      src_shift = shift_catalog.shiftNumberToShiftDistance(
            src_box.getPeriodicId());
      src_shift = (src_ratio > constant_zero_intvector) ?
         (src_shift * src_ratio) :
         hier::IntVector::ceilingDivide(src_shift, -src_ratio);
      unshifted_src_box.shift(-src_shift);
   }
   if (dst_box.isPeriodicImage()) {
      TBOX_ASSERT(!src_box.isPeriodicImage());
      const hier::IntVector& dst_ratio =
         d_dst_level->getRatioToLevelZero();
      dst_shift = shift_catalog.shiftNumberToShiftDistance(
            dst_box.getPeriodicId());
      dst_shift = (dst_ratio > constant_zero_intvector) ?
         (dst_shift * dst_ratio) :
         hier::IntVector::ceilingDivide(dst_shift, -dst_ratio);
      unshifted_dst_box.shift(-dst_shift);
   }
   if (s_extra_debug) {
      tbox::plog << "  src_shift: " << src_shift
                 << " unshifted_src_box " << unshifted_src_box << std::endl;
      tbox::plog << "  dst_shift: " << dst_shift
                 << " unshifted_dst_box " << unshifted_dst_box << std::endl;
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
   bool is_singularity = false;
   if (src_block_id != dst_block_id) {

      std::shared_ptr<hier::BaseGridGeometry> grid_geometry(
         d_dst_level->getGridGeometry());

      hier::Transformation::RotationIdentifier rotation =
         grid_geometry->getRotationIdentifier(dst_block_id,
            src_block_id);
      hier::IntVector offset(
         grid_geometry->getOffset(dst_block_id, src_block_id, d_dst_level->getLevelNumber()));

      is_singularity = grid_geometry->areSingularityNeighbors(dst_block_id,
            src_block_id);

      transformation = hier::Transformation(rotation, offset,
            src_block_id, dst_block_id);
      transformation.transform(transformed_src_box);
   }

   /*
    * For any case except when src and dst touch at enhanced connectivity,
    * the transactions use d_dst_level as the destination level and
    * dst_box as the destination box.  For the enhanced connectivity
    * case, the destination level becomes d_encon_level and the destination
    * box is a member of d_encon_level.
    */
   std::shared_ptr<hier::PatchLevel> transaction_dst_level;
   hier::Box transaction_dst_box(dim);
   if (is_singularity) {

      /*
       * Determination of transaction_dst_box is done differently
       * depending of if dst_box is local.  When it is local
       * (regardless of whether source is also local), the appropriate
       * Box to assign to transaction_dst_box can be
       * found from the d_dst_to_encon connector.
       *
       * If the destination is not local, then the source is, and
       * d_encon_to_src->getTranspose() is searched to fined the right member
       * of d_encon_level to assign to transaction_dst_box.
       */
      if (dst_box.getOwnerRank() == my_rank) {

         hier::Connector::ConstNeighborhoodIterator ei =
            d_dst_to_encon->findLocal(dst_box.getBoxId());

         if (ei == d_dst_to_encon->end()) {
            return;
         }

         const hier::BlockId& src_block_id = src_box.getBlockId();

         for (hier::Connector::ConstNeighborIterator en = d_dst_to_encon->begin(ei);
              en != d_dst_to_encon->end(ei); ++en) {

            if (src_block_id == en->getBlockId()) {
               TBOX_ASSERT(transaction_dst_box.empty());
               transaction_dst_box = *en;
            }
         }

      } else {

         TBOX_ASSERT(src_box.getOwnerRank() == my_rank);

         hier::Connector& src_to_encon = d_encon_to_src->getTranspose();
         hier::Connector::ConstNeighborhoodIterator ei =
            src_to_encon.findLocal(src_box.getBoxId());

         if (ei == src_to_encon.end()) {
            return;
         }

         hier::IntVector test_gcw(
            hier::IntVector::max(d_boundary_fill_ghost_width,
               d_max_stencil_width));
         test_gcw.max(hier::IntVector::getOne(dim));

         hier::Box test_dst_box(dst_box);
         test_dst_box.grow(test_gcw);

         hier::BoxContainer encon_nbr_choices;
         for (hier::Connector::ConstNeighborIterator ni = src_to_encon.begin(ei);
              ni != src_to_encon.end(ei); ++ni) {
            if (ni->getOwnerRank() == dst_box.getOwnerRank()) {
               hier::Box encon_box(*ni);
               transformation.transform(encon_box);

               if (test_dst_box.contains(encon_box)) {
                  encon_nbr_choices.pushBack(*ni);
               }
            }
         }
         if (encon_nbr_choices.empty()) {
            transaction_dst_box.setEmpty();
            transaction_dst_box.setBlockId(src_box.getBlockId());
         } else {
            size_t max_nbr_size = 0;
            for (hier::BoxContainer::iterator en = encon_nbr_choices.begin();
                 en != encon_nbr_choices.end(); ++en) {
               const size_t box_size = en->size();
               if (box_size > max_nbr_size) {
                  max_nbr_size = box_size;
                  transaction_dst_box = *en;
               }
            }
         }
      }

      transaction_dst_level = d_encon_level;

   } else {
      /*
       * All cases except for handling enhance connectivity neighbors
       * go here to do a simple assignment.
       */
      transaction_dst_level = d_dst_level;
      transaction_dst_box = dst_box;
   }

   for (int nc = 0; nc < num_equiv_classes; ++nc) {

      const RefineClasses::Data& rep_item =
         d_refine_classes->getClassRepresentative(nc);

      const int rep_item_dst_id = rep_item.d_scratch;
      const int rep_item_src_id = rep_item.d_src;

      std::shared_ptr<hier::PatchDataFactory> src_pdf(
         src_patch_descriptor->getPatchDataFactory(rep_item_src_id));
      std::shared_ptr<hier::PatchDataFactory> dst_pdf(
         dst_patch_descriptor->getPatchDataFactory(rep_item_dst_id));

      const hier::IntVector& dst_gcw = dst_pdf->getGhostCellWidth();

      hier::BoxContainer::iterator box_itr = d_src_masks.begin();
      int box_num = 0;
      for (hier::BoxNeighborhoodCollection::ConstNeighborIterator bi = nbrs_begin;
           bi != nbrs_end; ++bi) {

         const hier::Box& fill_box = *bi;

         /*
          * Get the patch data factories and calculate the overlap.
          * Note that we restrict the domain of the time interpolation
          * to the intersection of the fill box and the ghost box of
          * the destination patch data component.  This is needed for
          * the case where the schedule treats data components with
          * different ghost cell widths since the fill boxes are
          * generated using the largest ghost width.
          */

         hier::Box dst_fill_box(dst_box);
         dst_fill_box.grow(dst_gcw);
         dst_fill_box = dst_fill_box * fill_box;

         std::shared_ptr<hier::BoxOverlap> overlap;
         hier::Box src_mask(dim);
         if (!is_singularity) {

            /*
             * Create overlap for normal cases (all but enhanced connectivity).
             */

            hier::Box test_mask(dst_fill_box * transformed_src_box);
            if (test_mask.empty() && dst_pdf->dataLivesOnPatchBorder()) {

               test_mask = dst_fill_box;
               test_mask.grow(constant_one_intvector);
               test_mask = test_mask * transformed_src_box;

            }

            src_mask = test_mask;
            transformation.inverseTransform(src_mask);

            if (!src_mask.empty()) {
               overlap =
                  rep_item.d_var_fill_pattern->calculateOverlap(
                     *dst_pdf->getBoxGeometry(unshifted_dst_box),
                     *src_pdf->getBoxGeometry(unshifted_src_box),
                     dst_box,
                     src_mask,
                     fill_box,
                     true, transformation);
            }
         } else {

            /*
             * Create overlap for enhanced connectivity.  This overlap
             * will be used in transaction from d_src_level to d_encon_level.
             */

            hier::Box test_mask(dst_fill_box);
            transformation.inverseTransform(test_mask);
            test_mask = test_mask * src_box;
            if (test_mask.empty() && dst_pdf->dataLivesOnPatchBorder()) {

               test_mask = dst_fill_box;
               test_mask.grow(constant_one_intvector);
               transformation.inverseTransform(test_mask);
               test_mask = test_mask * src_box;

            }

            src_mask = test_mask;
            hier::Box transformed_fill_box(fill_box);
            hier::Box transformed_dst_box(dst_box);
            transformation.inverseTransform(transformed_fill_box);
            transformation.inverseTransform(transformed_dst_box);

            if (!src_mask.empty()) {
               overlap =
                  rep_item.d_var_fill_pattern->calculateOverlap(
                     *dst_pdf->getBoxGeometry(transaction_dst_box),
                     *src_pdf->getBoxGeometry(unshifted_src_box),
                     transformed_dst_box,
                     src_mask,
                     transformed_fill_box,
                     true, hier::Transformation(hier::IntVector::getZero(dim)));
            }
         }

#ifdef DEBUG_CHECK_ASSERTIONS
         if (!overlap && !src_mask.empty()) {
            TBOX_ERROR("Internal RefineSchedule error..."
               << "\n Overlap is NULL for "
               << "\n src box = " << src_box
               << "\n dst box = " << dst_box
               << "\n src mask = " << src_mask << std::endl);
         }
#endif

         if (s_extra_debug) {
            tbox::plog << "  overlap: ";
            overlap->print(tbox::plog);
         }
         *box_itr = src_mask;
         d_overlaps[box_num] = overlap;
         ++box_num;
         ++box_itr;

      }

      /*
       * Iterate over components in refine description list
       */
      for (std::list<int>::iterator l(d_refine_classes->getIterator(nc));
           l != d_refine_classes->getIteratorEnd(nc); ++l) {
         const RefineClasses::Data& item = d_refine_classes->getRefineItem(*l);
         TBOX_ASSERT(item.d_class_index == nc);
         TBOX_ASSERT(&item == d_refine_items[*l]);
         TBOX_ASSERT(item.d_tag == *l);

         const int dst_id = item.d_scratch;
         const int src_id = item.d_src;

         /*
          * If the src and dst patches, levels, and components are the
          * same, and there is no shift, the data exchange is unnecessary.
          */
         if (!same_patch_no_shift || (dst_id != src_id)) {

            /*
             * Iterate over the fill boxes and create transactions
             * for each box that has a non-empty overlap.
             */
            hier::BoxContainer::iterator itr = d_src_masks.begin();
            for (int i = 0; i < box_num; ++i, ++itr) {

               /*
                * If overlap is not empty, then add the transaction
                * to the appropriate communication schedule.
                * There are two schedules depending on whether
                * coarse or fine data takes precedence at
                * coarse-fine boundaries for communications
                * where the destination variable quantity
                * has data residing on the boundary.
                * There are two types of transactions depending on
                * whether we use time interpolation.
                */

               if (d_overlaps[i] && !d_overlaps[i]->isOverlapEmpty()) {

                  std::shared_ptr<tbox::Transaction> transaction;

                  if (d_transaction_factory) {

                     transaction =
                        d_transaction_factory->allocate(transaction_dst_level,
                           d_src_level,
                           d_overlaps[i],
                           transaction_dst_box,
                           src_box,
                           d_refine_items,
                           item.d_tag,
                           *itr,
                           (use_time_interpolation && item.d_time_interpolate));
                  } else if (use_time_interpolation &&
                             item.d_time_interpolate) {

                     transaction.reset(new RefineTimeTransaction(
                           transaction_dst_level, d_src_level,
                           d_overlaps[i],
                           transaction_dst_box, src_box,
                           *itr,
                           d_refine_items,
                           item.d_tag));

                  } else {  // no time interpolation

                     transaction.reset(new RefineCopyTransaction(
                           transaction_dst_level, d_src_level,
                           d_overlaps[i],
                           transaction_dst_box, src_box,
                           d_refine_items,
                           item.d_tag));

                  }  // time interpolation conditional

                  if (item.d_fine_bdry_reps_var) {
                     if (same_patch) {
                        d_fine_priority_level_schedule->addTransaction(
                           transaction);
                     } else {
                        d_fine_priority_level_schedule->appendTransaction(
                           transaction);
                     }
                  } else {
                     if (same_patch) {
                        d_coarse_priority_level_schedule->addTransaction(
                           transaction);
                     } else {
                        d_coarse_priority_level_schedule->appendTransaction(
                           transaction);
                     }
                  }

               }  // if overlap not empty

            }  // iterate over fill_boxes

         }

      }  // iterate over refine components in equivalence class

   }  // iterate over refine equivalence classes
}

/*
 *************************************************************************
 *
 * Private member function to initialize data members for hierarchy info.
 *
 *************************************************************************
 */

void
RefineSchedule::initializeDomainAndGhostInformation()
{

   const tbox::Dimension& dim(d_dst_level->getDim());

   d_max_scratch_gcw = getMaxScratchGhosts();
   d_max_stencil_width = getMaxStencilGhosts();

   if (d_top_refine_schedule == this) {
      // Not a recursive schedule.
      d_boundary_fill_ghost_width = getMaxDestinationGhosts();
      d_force_boundary_fill = false;
   } else {
      // Recursive schedule.
      d_boundary_fill_ghost_width = d_max_stencil_width;
      d_force_boundary_fill = (d_boundary_fill_ghost_width.max() > 0);
   }

   d_data_on_patch_border_flag = getDataOnPatchBorderFlag();

   std::shared_ptr<hier::BaseGridGeometry> grid_geom(
      d_dst_level->getGridGeometry());
   const hier::IntVector& ratio_to_level_zero =
      d_dst_level->getRatioToLevelZero();

   for (hier::BlockId::block_t b = 0; b < grid_geom->getNumberBlocks(); ++b) {
      d_domain_is_one_box[b] = grid_geom->getDomainIsSingleBox(hier::BlockId(b));
   }

   if (grid_geom->getNumberBlocks() > 1) {
      d_periodic_shift = hier::IntVector::getZero(dim);
   } else {
      d_periodic_shift =
         grid_geom->getPeriodicShift(
            ratio_to_level_zero);
   }

   d_num_periodic_directions = 0;
   for (int d = 0; d < dim.getValue(); ++d) {
      if (d_periodic_shift(d)) {
         ++d_num_periodic_directions;
      }
   }

}

/*
 *************************************************************************
 *
 * Private utility function to set up local array of refine items.
 *
 *************************************************************************
 */

void
RefineSchedule::setRefineItems(
   const std::shared_ptr<RefineClasses>& refine_classes)
{

   clearRefineItems();
   TBOX_ASSERT(d_number_refine_items == 0);

   d_refine_classes = refine_classes;

   d_number_refine_items = d_refine_classes->getNumberOfRefineItems();

   if (!d_refine_items) {
      d_refine_items =
         new const RefineClasses::Data *[d_number_refine_items];
   }

   for (int nd = 0; nd < static_cast<int>(d_number_refine_items); ++nd) {
      d_refine_classes->getRefineItem(nd).d_tag = nd;
      d_refine_items[nd] = &(d_refine_classes->getRefineItem(nd));
   }

}

/*
 *************************************************************************
 *
 * Private utility function used during initial schedule set up to
 * check whether patch data entries have proper number of ghost cells.
 * In particular, each scratch data entry must have at least as many
 * ghost cells as the user-defined refine operator stencil width.
 * Other checks are performed in the
 * RefineClasses::itemIsValid() routine.
 *
 *************************************************************************
 */

void
RefineSchedule::initialCheckRefineClassItems() const
{
   const tbox::Dimension& dim(d_dst_level->getDim());

   const hier::IntVector& constant_zero_intvector(
      hier::IntVector::getZero(dim));

   std::shared_ptr<hier::PatchDescriptor> pd(
      d_dst_level->getPatchDescriptor());

   hier::IntVector user_gcw(constant_zero_intvector);
   if (d_refine_patch_strategy) {
      user_gcw = d_refine_patch_strategy->getRefineOpStencilWidth(dim);
   }

   if (user_gcw > constant_zero_intvector) {

      for (size_t iri = 0; iri < d_number_refine_items; ++iri) {

         const RefineClasses::Data * const ref_item = d_refine_items[iri];

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_refine_classes->itemIsValid(*ref_item, pd)) {
#endif

         const int scratch = ref_item->d_scratch;
         const hier::IntVector& scratch_gcw(pd->getPatchDataFactory(scratch)->
                                            getGhostCellWidth());

         if (user_gcw > scratch_gcw) {
            TBOX_ERROR("Bad data given to RefineSchedule...\n"
               << "User supplied interpolation stencil width = "
               << user_gcw
               << "\nis larger than ghost cell width of `Scratch'\n"
               << "patch data " << pd->mapIndexToName(scratch)
               << " , which is " << scratch_gcw << std::endl);
         }

#ifdef DEBUG_CHECK_ASSERTIONS
      }
#endif

      }

   }

}

/*
 **************************************************************************
 *
 * Private utility function to clear array of refine items.
 *
 **************************************************************************
 */

void
RefineSchedule::clearRefineItems()
{
   if (d_refine_items) {
      for (size_t ici = 0; ici < d_number_refine_items; ++ici) {
         d_refine_items[ici] = 0;
      }
      d_number_refine_items = 0;
   }
}

/*
 **************************************************************************
 **************************************************************************
 */

void
RefineSchedule::setDeterministicUnpackOrderingFlag(bool flag)
{
   if (d_coarse_priority_level_schedule) {
      d_coarse_priority_level_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
   if (d_fine_priority_level_schedule) {
      d_fine_priority_level_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
   if (d_coarse_interp_schedule) {
      d_coarse_interp_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
   if (d_coarse_interp_encon_schedule) {
      d_coarse_interp_encon_schedule->setDeterministicUnpackOrderingFlag(flag);
   }
}

void
RefineSchedule::setScheduleOpsStrategy(tbox::ScheduleOpsStrategy* strategy)
{
   if (d_coarse_priority_level_schedule) {
      d_coarse_priority_level_schedule->setScheduleOpsStrategy(strategy);
   }
   if (d_fine_priority_level_schedule) {
      d_fine_priority_level_schedule->setScheduleOpsStrategy(strategy);
   }
   if (d_coarse_interp_schedule) {
      d_coarse_interp_schedule->setScheduleOpsStrategy(strategy);
   }
   if (d_coarse_interp_encon_schedule) {
      d_coarse_interp_encon_schedule->setScheduleOpsStrategy(strategy);
   }
}


/*
 **************************************************************************
 *
 * Allocate internal data
 *
 **************************************************************************
 */

void
RefineSchedule::allocateInternalData(double fill_time)
{
   if (d_internal_allocated) {
      return;
   }

   allocateScratchSpace(d_dst_scratch_vector, d_dst_level, fill_time);

   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      allocateScratchSpace(d_encon_scratch_vector, d_encon_level, fill_time);
   }

   if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
       d_nbr_blk_fill_level.get()) {
      allocateScratchSpace(d_nbr_fill_scratch_vector,
                           d_nbr_blk_fill_level,
                           fill_time);
      allocateDestinationSpace(d_nbr_fill_dst_vector,
                               d_nbr_blk_fill_level,
                               fill_time);
   }

   if (d_coarse_interp_schedule) {
      allocateScratchSpace(d_coarse_scratch_vector,
         d_coarse_interp_level,
         fill_time);
      allocateWorkSpace(d_coarse_work_vector,
         d_coarse_interp_level,
         fill_time);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         allocateScratchSpace(d_coarse_encon_scratch_vector,
            d_coarse_interp_schedule->d_encon_level,
            fill_time);
         allocateWorkSpace(d_coarse_encon_work_vector,
            d_coarse_interp_schedule->d_encon_level,
            fill_time);
      }

      if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
          d_coarse_interp_schedule->d_nbr_blk_fill_level.get()) {
         allocateScratchSpace(d_coarse_nbr_fill_scratch_vector,
            d_coarse_interp_schedule->d_nbr_blk_fill_level, fill_time);
         allocateWorkSpace(d_coarse_nbr_fill_work_vector,
            d_coarse_interp_schedule->d_nbr_blk_fill_level, fill_time);
      }

      if (d_nbr_blk_fill_level.get()) {
         allocateWorkSpace(d_nbr_fill_work_vector,
            d_nbr_blk_fill_level, fill_time);
      }

      d_coarse_interp_schedule->allocateInternalData(fill_time);
   }

   if (d_coarse_interp_encon_schedule) {

      allocateScratchSpace(d_coarse_interp_encon_scratch_vector,
                           d_coarse_interp_encon_level,
                           fill_time);
      allocateWorkSpace(d_coarse_interp_encon_work_vector,
                        d_coarse_interp_encon_level,
                        fill_time);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         allocateScratchSpace(d_coarse_encon_encon_scratch_vector,
            d_coarse_interp_encon_schedule->d_encon_level,
            fill_time);
         allocateWorkSpace(d_coarse_encon_encon_work_vector,
            d_coarse_interp_encon_schedule->d_encon_level,
            fill_time);
      }

      d_coarse_interp_encon_schedule->allocateInternalData(fill_time);
   }
   d_internal_allocated = true;
}

/*
 **************************************************************************
 *
 * Set timestamps on internal data
 * 
 **************************************************************************
 */

void
RefineSchedule::setInternalDataTime(double fill_time) const
{
   d_dst_level->setTime(fill_time, d_dst_scratch_vector);

   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      d_encon_level->setTime(fill_time, d_encon_scratch_vector);
   }

   if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
       d_nbr_blk_fill_level.get()) {
      d_nbr_blk_fill_level->setTime(fill_time, d_nbr_fill_scratch_vector);
      d_nbr_blk_fill_level->setTime(fill_time, d_nbr_fill_dst_vector);
   }

   if (d_coarse_interp_schedule) {
      d_coarse_interp_level->setTime(fill_time, d_coarse_scratch_vector);
      d_coarse_interp_level->setTime(fill_time, d_coarse_work_vector);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_schedule->d_encon_level->setTime(fill_time,
            d_coarse_encon_scratch_vector);
         d_coarse_interp_schedule->d_encon_level->setTime(fill_time,
            d_coarse_encon_work_vector);
      }

      if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
          d_coarse_interp_schedule->d_nbr_blk_fill_level.get()) {
         d_coarse_interp_schedule->d_nbr_blk_fill_level->setTime(fill_time,
            d_coarse_nbr_fill_scratch_vector); 
         d_coarse_interp_schedule->d_nbr_blk_fill_level->setTime(fill_time,
            d_coarse_nbr_fill_work_vector); 
      }

      if (d_nbr_blk_fill_level.get()) {
         d_nbr_blk_fill_level->setTime(fill_time, d_nbr_fill_work_vector);
      }

      d_coarse_interp_schedule->setInternalDataTime(fill_time);
   }

   if (d_coarse_interp_encon_schedule) {

      d_coarse_interp_encon_level->setTime(fill_time,
         d_coarse_interp_encon_scratch_vector);
      d_coarse_interp_encon_level->setTime(fill_time,
         d_coarse_interp_encon_work_vector);

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_encon_schedule->d_encon_level->setTime(fill_time,
            d_coarse_encon_encon_scratch_vector);
         d_coarse_interp_encon_schedule->d_encon_level->setTime(fill_time,
            d_coarse_encon_encon_work_vector);
      }

      d_coarse_interp_encon_schedule->setInternalDataTime(fill_time);
   }
}

/*
 **************************************************************************
 *
 * Deallocate internal data
 *
 **************************************************************************
 */

void
RefineSchedule::deallocateInternalData()
{

   d_dst_level->deallocatePatchData(d_dst_scratch_vector);
   d_dst_scratch_vector.clrAllFlags();

   if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
      d_encon_level->deallocatePatchData(d_encon_scratch_vector);
      d_encon_scratch_vector.clrAllFlags();
   }

   if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
       d_nbr_blk_fill_level.get()) {
      d_nbr_blk_fill_level->deallocatePatchData(d_nbr_fill_scratch_vector);
      d_nbr_blk_fill_level->deallocatePatchData(d_nbr_fill_dst_vector);
      d_nbr_fill_scratch_vector.clrAllFlags();
      d_nbr_fill_dst_vector.clrAllFlags();
   }

   if (d_coarse_interp_schedule) {
      d_coarse_interp_level->deallocatePatchData(d_coarse_scratch_vector);
      d_coarse_interp_level->deallocatePatchData(d_coarse_work_vector);
      d_coarse_scratch_vector.clrAllFlags();
      d_coarse_work_vector.clrAllFlags();

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_schedule->d_encon_level->deallocatePatchData(
            d_coarse_encon_scratch_vector);
         d_coarse_interp_schedule->d_encon_level->deallocatePatchData(
            d_coarse_encon_work_vector);
         d_coarse_encon_scratch_vector.clrAllFlags();
         d_coarse_encon_work_vector.clrAllFlags();
      }

      if (d_dst_level->getGridGeometry()->getNumberBlocks() > 1 &&
          d_coarse_interp_schedule->d_nbr_blk_fill_level.get()) {
         d_coarse_interp_schedule->d_nbr_blk_fill_level->deallocatePatchData(
            d_coarse_nbr_fill_scratch_vector);
         d_coarse_interp_schedule->d_nbr_blk_fill_level->deallocatePatchData(
            d_coarse_nbr_fill_work_vector);
         d_coarse_nbr_fill_scratch_vector.clrAllFlags(); 
         d_coarse_nbr_fill_work_vector.clrAllFlags(); 
      }

      if (d_nbr_blk_fill_level.get()) {
         d_nbr_blk_fill_level->deallocatePatchData(d_nbr_fill_work_vector);
         d_nbr_fill_work_vector.clrAllFlags();
      }

      d_coarse_interp_schedule->deallocateInternalData();
   }

   if (d_coarse_interp_encon_schedule) {

      d_coarse_interp_encon_level->deallocatePatchData(
         d_coarse_interp_encon_scratch_vector);
      d_coarse_interp_encon_level->deallocatePatchData(
         d_coarse_interp_encon_work_vector);
      d_coarse_interp_encon_scratch_vector.clrAllFlags();
      d_coarse_interp_encon_work_vector.clrAllFlags();

      if (d_dst_level->getGridGeometry()->hasEnhancedConnectivity()) {
         d_coarse_interp_encon_schedule->d_encon_level->deallocatePatchData(
            d_coarse_encon_encon_scratch_vector);
         d_coarse_interp_encon_schedule->d_encon_level->deallocatePatchData(
            d_coarse_encon_encon_work_vector);
         d_coarse_encon_encon_scratch_vector.clrAllFlags();
         d_coarse_encon_encon_work_vector.clrAllFlags();
      }

      d_coarse_interp_encon_schedule->deallocateInternalData();
   }

   d_internal_allocated = false;
}


/*
 **************************************************************************
 *
 * Print class data to the specified output stream.
 *
 **************************************************************************
 */

void
RefineSchedule::printClassData(
   std::ostream& stream) const
{
   stream << "RefineSchedule::printClassData()\n";
   stream << "--------------------------------------\n";

   d_refine_classes->printClassData(stream);

   stream << "Printing coarse priority refine schedule...\n";
   d_coarse_priority_level_schedule->printClassData(stream);

   stream << "Printing fine priority refine schedule...\n";
   d_fine_priority_level_schedule->printClassData(stream);

   if (d_coarse_interp_schedule) {
      stream << "Printing coarse interpolation refine schedule...\n";
      d_coarse_interp_schedule->printClassData(stream);
   }
}

/*
 **************************************************************************
 **************************************************************************
 */

void
RefineSchedule::initializeCallback()
{
   t_refine_schedule = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::RefineSchedule()");
   t_fill_data = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillData()");
   t_fill_data_nonrecursive = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillData()_nonrecursive");
   t_fill_data_recursive = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillData()_recursive");
   t_fill_physical_boundaries = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillPhysicalBoundaries()");
   t_fill_singularity_boundaries = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::fillSingularityBoundaries()");
   t_refine_scratch_data = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::refineScratchData()");
   t_finish_sched_const = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finishScheduleConstruction()");
   t_finish_sched_const_recurse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finishScheduleConstruction()_recurse");
   t_gen_comm_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generateCommunicationSchedule()");
   t_shear = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_shear");
   t_get_global_box_count = tbox::TimerManager::getManager()->
      getTimer(
         "xfer::RefineSchedule::finish...()_get_global_box_count");
   t_coarse_shear = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_coarse_shear");
   t_setup_coarse_interp_box_level = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::setupCoarseInterpBoxLevel()");
   t_bridge_coarse_interp_hiercoarse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_bridge_coarse_interp_hiercoarse");
   t_bridge_dst_hiercoarse = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::finish...()_bridge_dst_hiercoarse");
   t_invert_edges = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_invert_edges");
   t_construct_send_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_construct_send_trans");
   t_construct_recv_trans = tbox::TimerManager::getManager()->
      getTimer("xfer::RefineSchedule::generate...()_construct_recv_trans");

}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
void
RefineSchedule::finalizeCallback()
{
   t_fill_data.reset();
   t_fill_data_nonrecursive.reset();
   t_fill_data_recursive.reset();
   t_refine_scratch_data.reset();
   t_finish_sched_const.reset();
   t_finish_sched_const_recurse.reset();
   t_gen_comm_sched.reset();
   t_shear.reset();
   t_get_global_box_count.reset();
   t_coarse_shear.reset();
   t_setup_coarse_interp_box_level.reset();
   t_bridge_coarse_interp_hiercoarse.reset();
   t_bridge_dst_hiercoarse.reset();
   t_invert_edges.reset();
   t_construct_send_trans.reset();
   t_construct_recv_trans.reset();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
