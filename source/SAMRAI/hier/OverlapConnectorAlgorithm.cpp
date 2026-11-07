/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with overlap Connectors.
 *
 ************************************************************************/
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/AssumedPartition.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerUtils.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <iterator>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const std::string OverlapConnectorAlgorithm::s_default_timer_prefix(
   "hier::OverlapConnectorAlgorithm");
std::map<std::string,
         OverlapConnectorAlgorithm::TimerStruct> OverlapConnectorAlgorithm::s_static_timers;
char OverlapConnectorAlgorithm::s_ignore_external_timer_prefix('n');

char OverlapConnectorAlgorithm::s_print_steps = '\0';

int OverlapConnectorAlgorithm::s_operation_mpi_tag = 0;
/*
 * Do we even need to use different tags each time we bridge???
 * Unique tags were used to help debug, but the methods may work
 * with reused tags anyway.
 */

tbox::StartupShutdownManager::Handler
OverlapConnectorAlgorithm::s_initialize_finalize_handler(
   OverlapConnectorAlgorithm::initializeCallback,
   0,
   0,
   OverlapConnectorAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

OverlapConnectorAlgorithm::OverlapConnectorAlgorithm():
   d_mpi(MPI_COMM_NULL),
   d_mpi_is_exclusive(false),
   d_object_timers(0),
   d_print_steps(s_print_steps == 'y'),
   d_barrier_before_communication(false),
   d_sanity_check_method_preconditions(false),
   d_sanity_check_method_postconditions(false)
{
   getFromInput();
   setTimerPrefix(s_default_timer_prefix);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

OverlapConnectorAlgorithm::~OverlapConnectorAlgorithm()
{
   if (d_mpi_is_exclusive) {
      d_mpi.freeCommunicator();
      d_mpi_is_exclusive = false;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void OverlapConnectorAlgorithm::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& mpi,
   bool make_duplicate)
{
   if (d_mpi_is_exclusive) {
      d_mpi.freeCommunicator();
      d_mpi_is_exclusive = false;
   }
   if (make_duplicate) {
      d_mpi.dupCommunicator(mpi);
      d_mpi_is_exclusive = true;
   } else {
      d_mpi = mpi;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::getFromInput()
{
   if (s_print_steps == '\0') {
      s_print_steps = 'n';
      if (tbox::InputManager::inputDatabaseExists()) {
         std::shared_ptr<tbox::Database> idb(
            tbox::InputManager::getInputDatabase());
         if (idb->isDatabase("OverlapConnectorAlgorithm")) {
            std::shared_ptr<tbox::Database> oca_db(
               idb->getDatabase("OverlapConnectorAlgorithm"));
            s_print_steps =
               oca_db->getCharWithDefault("DEV_print_bridge_steps", 'n');
            if (!(s_print_steps == 'n' || s_print_steps == 'y')) {
               INPUT_VALUE_ERROR("DEV_print_bridge_steps");
            }
            s_ignore_external_timer_prefix =
               oca_db->getCharWithDefault("DEV_ignore_external_timer_prefix",
                  'n');
            if (!(s_ignore_external_timer_prefix == 'n' ||
                  s_ignore_external_timer_prefix == 'y')) {
               INPUT_VALUE_ERROR("DEV_ignore_external_timer_prefix");
            }
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::extractNeighbors(
   NeighborSet& neighbors,
   const Connector& connector,
   const BoxId& box_id,
   const IntVector& width) const
{
   const tbox::Dimension& dim(width.getDim());

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!(width <= connector.getConnectorWidth())) {
      TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot provide\n"
         << "neighbors for a wider ghost cell width that used to initialize it.\n");
   }
   if (connector.getParallelState() != BoxLevel::GLOBALIZED &&
       box_id.getOwnerRank() != connector.getMPI().getRank()) {
      TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot get\n"
         << "neighbor data for a remote box unless in GLOBALIZED mode.\n");
   }
   if (!connector.getBase().hasBox(box_id)) {
      std::string dbgbord;
      TBOX_ERROR(
         "\nOverlapConnectorAlgorithm::extractNeighbors: box_id " << box_id
         << " is not in the base of the box_level.\n"
         << "base:\n" << connector.getBase().format(dbgbord, 2)
         << "head:\n" << connector.getHead().format(dbgbord, 2)
         << "connector:\n" << connector.format(dbgbord, 2));
   }
#endif

   /*
    * Temporarily disable extracting neighbors for remote boxes.  This
    * method functionality is not much used and prrobably should be
    * removed.
    */
   TBOX_ASSERT(box_id.getOwnerRank() == connector.getMPI().getRank());

   const std::shared_ptr<const BaseGridGeometry>& grid_geom(
      connector.getBase().getGridGeometry());

   const Box& box(*connector.getBase().getBox(Box(dim, box_id)));
   Connector::ConstNeighborhoodIterator ins =
      connector.findLocal(box_id);
   neighbors.clear();
   if (ins != connector.end()) {
      if (width == connector.getConnectorWidth()) {
         for (Connector::ConstNeighborIterator ni = connector.begin(ins);
              ni != connector.end(ins); ++ni) {
            neighbors.insert(neighbors.end(), *ni);
         }
      }
      else {
         BoxContainer grown_boxes;
         if (grid_geom->getNumberBlocks() == 1 ||
             grid_geom->hasIsotropicRatios()) {
            Box grown_box = box;
            grown_box.grow(width);
            if (connector.getHeadCoarserFlag() == false) {
               grown_box.refine(connector.getRatio());
            }
            grown_boxes.pushBack(grown_box);
         } else {
            bool head_is_finer = false;
            bool base_is_finer = false;
            if (connector.getRatio() != 1) {
               if (connector.getHeadCoarserFlag()) {
                  base_is_finer = true;
               } else {
                  head_is_finer = true;
               }
            }
            BoxUtilities::growAndAdjustAcrossBlockBoundary(
               grown_boxes,
               box,
               grid_geom,
               connector.getBase().getRefinementRatio(),
               connector.getRatio(),
               width, 
               head_is_finer,
               base_is_finer);
         }
         for (BoxContainer::iterator g_itr = grown_boxes.begin();
              g_itr != grown_boxes.end(); ++g_itr) {

            const Box& g_box = *g_itr;
            const BlockId& g_block_id = g_box.getBlockId();
            for (Connector::ConstNeighborIterator ni = connector.begin(ins);
                 ni != connector.end(ins); ++ni) {
               const Box& neighbor(*ni);
               Box nabr_box(neighbor);
               bool do_intersect = true;
               if (neighbor.getBlockId() != g_block_id) {
                  do_intersect =
                     grid_geom->transformBox(nabr_box,
                        connector.getHead().getRefinementRatio(),
                        g_block_id,
                        neighbor.getBlockId());
               }
               if (do_intersect) {
                  if (g_box.intersects(nabr_box)) {
                     neighbors.insert(neighbors.end(), neighbor);
                  }
               }
            }
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::extractNeighbors(
   Connector& other,
   const Connector& connector,
   const IntVector& width) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!(width <= connector.getConnectorWidth())) {
      TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot provide\n"
         << "neighbors for a wider ghost cell width that used to initialize it.\n");
   }
#endif

   other.clearNeighborhoods();
   for (Connector::ConstNeighborhoodIterator ni = connector.begin();
        ni != connector.end(); ++ni) {

      const BoxId& box_id = *ni;
      Connector::NeighborhoodIterator base_box_itr =
         other.makeEmptyLocalNeighborhood(box_id);

      const tbox::Dimension& dim(width.getDim());

#ifdef DEBUG_CHECK_ASSERTIONS
      if (connector.getParallelState() != BoxLevel::GLOBALIZED &&
          box_id.getOwnerRank() != connector.getMPI().getRank()) {
         TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot get\n"
            << "neighbor data for a remote box unless in GLOBALIZED mode.\n");
      }
      if (!connector.getBase().hasBox(box_id)) {
         std::string dbgbord;
         TBOX_ERROR(
            "\nOverlapConnectorAlgorithm::extractNeighbors: box_id " << box_id
                                                                     <<
            " is not in the base of the box_level.\n"
            << "base:\n" << connector.getBase().format(dbgbord, 2)
            << "head:\n" << connector.getHead().format(dbgbord, 2)
            << "connector:\n" << connector.format(dbgbord, 2));
      }
#endif

      /*
       * Temporarily disable extracting neighbors for remote boxes.  This
       * method functionality is not much used and prrobably should be
       * removed.
       */
      TBOX_ASSERT(box_id.getOwnerRank() == connector.getMPI().getRank());

      const std::shared_ptr<const BaseGridGeometry>& grid_geom(
         connector.getBase().getGridGeometry());

      const Box& box = *connector.getBase().getBox(Box(dim, box_id));

      if (width == connector.getConnectorWidth()) {
         for (Connector::ConstNeighborIterator si = connector.begin(ni);
              si != connector.end(ni); ++si) {
            other.insertLocalNeighbor(*si, base_box_itr);
         }
      } else {
         BoxContainer grown_boxes;
         if (grid_geom->getNumberBlocks() == 1 ||
             grid_geom->hasIsotropicRatios()) {

            Box grown_box = box;
            grown_box.grow(width);
            if (connector.getHeadCoarserFlag()) {
               grown_box.coarsen(connector.getRatio());
            } else {
               grown_box.refine(connector.getRatio());
            }
            grown_boxes.pushBack(grown_box);
         } else {
            bool head_is_finer = false;
            bool base_is_finer = false;
            if (connector.getRatio() != 1) {
               if (connector.getHeadCoarserFlag()) {
                  base_is_finer = true;
               } else {
                  head_is_finer = true;
               }
            }
            BoxUtilities::growAndAdjustAcrossBlockBoundary(
               grown_boxes,
               box,
               grid_geom,
               connector.getBase().getRefinementRatio(),
               connector.getRatio(),
               width,
               head_is_finer,
               base_is_finer);
         }
         for (BoxContainer::iterator g_itr = grown_boxes.begin();
              g_itr != grown_boxes.end(); ++g_itr) {
            const Box& g_box = *g_itr;
            const BlockId& g_block_id = g_box.getBlockId();
            for (Connector::ConstNeighborIterator si = connector.begin(ni);
                 si != connector.end(ni); ++si) {
               const Box& neighbor = *si;
               Box nabr_box(neighbor);
               bool do_intersect = true;
               if (neighbor.getBlockId() != g_block_id) {
                  do_intersect = grid_geom->transformBox(nabr_box,
                     connector.getHead().getRefinementRatio(),
                     g_block_id,
                     neighbor.getBlockId());
               }
               if (do_intersect) {
                  if (g_box.intersects(nabr_box)) {
                     other.insertLocalNeighbor(neighbor, base_box_itr);
                  }
               }
            }
         }
      }
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps(
   std::shared_ptr<Connector>& connector,
   const BoxLevel& base_box_level,
   const BoxLevel& head_box_level,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state,
   const bool ignore_self_overlap) const
{
   connector.reset(new Connector(base_box_level,
      head_box_level,
      base_width,
      parallel_state));
   findOverlaps(*connector,
      head_box_level.getGlobalizedVersion(),
      ignore_self_overlap);
   if (&base_box_level == &head_box_level) {
      connector->setTranspose(connector.get(), false);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlapsWithTranspose(
   std::shared_ptr<Connector>& connector,
   const BoxLevel& base_box_level,
   const BoxLevel& head_box_level,
   const IntVector& base_width,
   const IntVector& transpose_base_width,
   const BoxLevel::ParallelState parallel_state,
   const bool ignore_self_overlap) const
{
   findOverlaps(connector,
      base_box_level,
      head_box_level,
      base_width,
      parallel_state,
      ignore_self_overlap);
   if (&base_box_level != &head_box_level) {
      Connector* transpose = new Connector(head_box_level,
         base_box_level,
         transpose_base_width,
         parallel_state);
      findOverlaps(*transpose, ignore_self_overlap);
      connector->setTranspose(transpose, true);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps(
   Connector& connector,
   const bool ignore_self_overlap) const
{
   findOverlaps(connector,
      connector.getHead().getGlobalizedVersion(),
      ignore_self_overlap);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps(
   Connector& connector,
   const BoxLevel& globalized_head,
   const bool ignore_self_overlap) const
{
   d_object_timers->t_find_overlaps_rbbt->start();
   connector.findOverlaps_rbbt(globalized_head,
      ignore_self_overlap,
      d_sanity_check_method_postconditions);
   d_object_timers->t_find_overlaps_rbbt->stop();
}

/*
 ***********************************************************************
 * Find overlaps using the assumed partition algorithm.  In SAMRAI
 * terms, we create a center BoxLevel from the assumed partition,
 * connect the base and head BoxLevels to it, then bridge across the
 * assumed partition center to get base<==>head.
 *
 * 1. Get bounding boxes for base or head, which ever has smaller
 * bounding boxes.
 *
 * 2. Create a reasonably balanced center BoxLevel from an
 * AssumedPartition of the bounding boxes.
 *
 * 3. Populate head--->center and base--->center.
 *
 * 4. Get transposes center--->head and center--->base.  This is
 * implemented in the Connector's computeTransposeOf() method.
 *
 * 5. Bridge base<==>center<==>head.
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::findOverlaps_assumedPartition(
   Connector& conn) const
{
   if (d_print_steps) {
      tbox::plog << "OverlapConnectorAlgorithm::findOverlaps_assumedPartition: entered.\n";
   }
   d_object_timers->t_find_overlaps_assumed_partition->start();

   const BoxLevel& base = conn.getBase();
   const BoxLevel& head = conn.getHead();
   const BoxContainer& head_boxes = head.getBoxes();
   const BoxContainer& base_boxes = base.getBoxes();

   const IntVector& width_in_base_resolution = conn.getConnectorWidth();
   const IntVector width_in_head_resolution =
      Connector::convertHeadWidthToBase(head.getRefinementRatio(),
         base.getRefinementRatio(),
         width_in_base_resolution);

   const tbox::Dimension& dim = base.getDim();
   const tbox::SAMRAI_MPI& mpi = d_mpi.hasNullCommunicator() ? base.getMPI() : d_mpi;
   const std::shared_ptr<const BaseGridGeometry>& geom = base.getGridGeometry();
   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("OverlapConnectorAlgorithm::findOverlaps_assumedPartition: not starting\n"
         << "clean of receivable MPI messages.");
   }

   if (d_sanity_check_method_preconditions) {
      if (!d_mpi.hasNullCommunicator() && !d_mpi.isCongruentWith(base.getMPI())) {
         TBOX_ERROR("OverlapConnectorAlgorithm::findOverlaps_assumedPartition input error: Input\n"
            << "has SAMRAI_MPI that is incongruent with OverlapConnectorAlgorithm's.\n"
            << "See OverlapConnectorAlgorithm::setSAMRAI_MPI.\n");
      }
   }

   d_object_timers->t_find_overlaps_assumed_partition_get_ap->barrierAndStart();

   /*
    * Set up center BoxLevel.  We can use either the base or head to
    * construct the center.  We choose the smaller one because we
    * don't need to cover the bigger region.  There are no overlaps
    * away from the smaller BoxLevel anyway.
    *
    * As an optimization, try to make an AssumedPartition with about
    * 1 partition per rank, or fewer if the base or head have fewer boxes.
    * This avoids having more many more assumed partitions than boxes,
    * which can scale poorly.
    */
   BoxContainer base_bounding_boxes, head_bounding_boxes;
   size_t base_bounding_cell_count = 0, head_bounding_cell_count = 0;
   for (BlockId::block_t bn = 0; bn < geom->getNumberBlocks(); ++bn) {
      BlockId block_id(bn);
      base_bounding_boxes.push_back(base.getGlobalBoundingBox(block_id));
      head_bounding_boxes.push_back(head.getGlobalBoundingBox(block_id));
      base_bounding_cell_count += base_bounding_boxes.back().size();
      head_bounding_cell_count += head_bounding_boxes.back().size();
   }
   size_t num_parts = head.getGlobalNumberOfBoxes() < base.getGlobalNumberOfBoxes() ?
      head.getGlobalNumberOfBoxes() : base.getGlobalNumberOfBoxes();
   num_parts = static_cast<size_t>(mpi.getSize()) < num_parts ?
      static_cast<size_t>(mpi.getSize()) : num_parts;
   const AssumedPartition center_ap(
      head_bounding_cell_count < base_bounding_cell_count ? head_bounding_boxes : base_bounding_boxes,
      0, mpi.getSize(), 0, static_cast<double>(num_parts)/mpi.getSize() );
   base_bounding_boxes.clear();
   head_bounding_boxes.clear();

   BoxContainer center_boxes;
   center_ap.getAllBoxes(center_boxes, mpi.getRank());
   const IntVector& center_refinement_ratio = head_bounding_cell_count < base_bounding_cell_count ?
      head.getRefinementRatio() : base.getRefinementRatio();
   const BoxLevel center( center_boxes, center_refinement_ratio, geom, mpi );

   d_object_timers->t_find_overlaps_assumed_partition_get_ap->barrierAndStop();

   d_object_timers->t_find_overlaps_assumed_partition_connect_to_ap->barrierAndStart();

   // Set up base<==>center
   Connector base_to_center(base, center, width_in_base_resolution);
   BoxContainer base_boxes_mod(base_boxes);
   base_boxes_mod.grow(width_in_base_resolution);
   if (base.getRefinementRatio() != center_refinement_ratio) {
      if (base.getRefinementRatio() >= center_refinement_ratio) {
         base_boxes_mod.refine(conn.getRatio());
      } else {
         base_boxes_mod.coarsen(conn.getRatio());
      }
   }
   for (BoxContainer::const_iterator bi = base_boxes_mod.begin(); bi != base_boxes_mod.end();
        ++bi) {
      const Box& compare_box = *bi;
      BoxContainer neighbors;
      center_ap.findOverlaps(neighbors, compare_box, *geom, center_refinement_ratio);
      base_to_center.insertNeighbors(neighbors, bi->getBoxId());
   }
   base_boxes_mod.clear();

   // Set up head<==>center
   Connector head_to_center(head, center, width_in_head_resolution);
   BoxContainer head_boxes_mod(head_boxes);
   head_boxes_mod.grow(width_in_head_resolution);
   if (head.getRefinementRatio() != center_refinement_ratio) {
      if (head.getRefinementRatio() >= center_refinement_ratio) {
         head_boxes_mod.refine(conn.getRatio());
      } else {
         head_boxes_mod.coarsen(conn.getRatio());
      }
   }
   for (BoxContainer::const_iterator bi = head_boxes_mod.begin(); bi != head_boxes_mod.end();
        ++bi) {
      Box compare_box = *bi;
      BoxContainer neighbors;
      center_ap.findOverlaps(neighbors, compare_box, *geom, center_refinement_ratio);
      head_to_center.insertNeighbors(neighbors, bi->getBoxId());
   }
   head_boxes_mod.clear();

   d_object_timers->t_find_overlaps_assumed_partition_connect_to_ap->barrierAndStop();

   d_object_timers->t_find_overlaps_assumed_partition_transpose->start();
   if (d_print_steps) {
      tbox::plog
      << "OverlapConnectorAlgorithm::findOverlaps_assumedPartition: getting transposes.\n";
   }
   Connector center_to_base(dim), center_to_head(dim);
   center_to_base.computeTransposeOf(base_to_center, mpi);
   center_to_head.computeTransposeOf(head_to_center, mpi);
   d_object_timers->t_find_overlaps_assumed_partition_transpose->barrierAndStop();

   // Bridge for base<==>head
   base_to_center.setTranspose(&center_to_base, false);
   head_to_center.setTranspose(&center_to_head, false);
   const IntVector center_growth_to_nest_base(
      dim,
      head_bounding_cell_count <
      base_bounding_cell_count ? tbox::MathUtilities<int>::getMax() : 0);
   const IntVector center_growth_to_nest_head(
      dim,
      head_bounding_cell_count < base_bounding_cell_count ? 0 : tbox::MathUtilities<int>::getMax());
   std::shared_ptr<Connector> tmp_conn;
   if (d_print_steps) {
      tbox::plog << "OverlapConnectorAlgorithm::findOverlaps_assumedPartition: bridging.\n";
   }
   bridgeWithNesting(tmp_conn,
      base_to_center,
      center_to_head,
      center_growth_to_nest_base,
      center_growth_to_nest_head,
      IntVector(dim, -1),
      false);
   conn.clear();
   conn.setBase(tmp_conn->getBase());
   conn.setHead(tmp_conn->getHead());
   conn.setWidth(tmp_conn->getConnectorWidth(), true);
   for (Connector::NeighborhoodIterator ni = tmp_conn->begin(); ni != tmp_conn->end(); ++ni) {
      for (Connector::NeighborIterator na = tmp_conn->begin(ni); na != tmp_conn->end(ni); ++na) {
         conn.insertLocalNeighbor(*na, *ni);
      }
   }

   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("OverlapConnectorAlgorithm::findOverlaps_assumedPartition: not finishing\n"
         << "clean of receivable MPI messages.");
   }

   d_object_timers->t_find_overlaps_assumed_partition->stop();
   if (d_print_steps) {
      tbox::plog << "OverlapConnectorAlgorithm::findOverlaps_assumedPartition: leaving.\n";
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridgeWithNesting(
   std::shared_ptr<Connector>& west_to_east,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const IntVector& cent_growth_to_nest_west,
   const IntVector& cent_growth_to_nest_east,
   const IntVector& connector_width_limit,
   bool compute_transpose) const
{
   const tbox::SAMRAI_MPI& mpi =
      d_mpi.hasNullCommunicator() ? west_to_cent.getBase().getMPI() : d_mpi;

   if (d_barrier_before_communication) {
      mpi.Barrier();
   }
   d_object_timers->t_bridge->start();

   TBOX_ASSERT(west_to_cent.hasTranspose());
   TBOX_ASSERT(cent_to_east.hasTranspose());
   const tbox::Dimension& dim(connector_width_limit.getDim());
   IntVector west_to_east_width(dim);
   IntVector east_to_west_width(dim);
   std::set<int> incoming_ranks, outgoing_ranks;
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   privateBridge_prologue(
      west_to_cent,
      cent_to_east,
      cent_to_east.getTranspose(),
      west_to_cent.getTranspose(),
      (cent_growth_to_nest_west(0) >= 0),
      cent_growth_to_nest_west,
      (cent_growth_to_nest_east(0) >= 0),
      cent_growth_to_nest_east,
      connector_width_limit,
      compute_transpose,
      west_to_east_width,
      east_to_west_width,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   Connector* east_to_west = nullptr;
   west_to_east.reset(new Connector(west_to_cent.getBase(),
      cent_to_east.getHead(),
      west_to_east_width));
   if (compute_transpose) {
      east_to_west = new Connector(cent_to_east.getHead(),
         west_to_cent.getBase(),
         east_to_west_width);
   }
   privateBridge(
      *west_to_east,
      east_to_west,
      cent_to_east,
      compute_transpose,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   if (compute_transpose) {
      west_to_east->setTranspose(east_to_west, true);
   } else if (&west_to_east->getHead() == &west_to_east->getBase()) {
      west_to_east->setTranspose(west_to_east.get(), false);
   }

   d_object_timers->t_bridge->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   std::shared_ptr<Connector>& west_to_east,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const IntVector& connector_width_limit,
   bool compute_transpose) const
{
   const tbox::SAMRAI_MPI& mpi =
      d_mpi.hasNullCommunicator() ? west_to_cent.getBase().getMPI() : d_mpi;

   if (d_barrier_before_communication) {
      mpi.Barrier();
   }
   d_object_timers->t_bridge->start();

   TBOX_ASSERT(west_to_cent.hasTranspose());
   TBOX_ASSERT(cent_to_east.hasTranspose());
   const tbox::Dimension& dim(connector_width_limit.getDim());
   const IntVector& zero_vector(IntVector::getZero(dim));
   IntVector west_to_east_width(dim);
   IntVector east_to_west_width(dim);
   std::set<int> incoming_ranks, outgoing_ranks;
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   privateBridge_prologue(
      west_to_cent,
      cent_to_east,
      cent_to_east.getTranspose(),
      west_to_cent.getTranspose(),
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit,
      compute_transpose,
      west_to_east_width,
      east_to_west_width,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   Connector* east_to_west = nullptr;
   west_to_east.reset(new Connector(west_to_cent.getBase(),
      cent_to_east.getHead(),
      west_to_east_width));
   if (compute_transpose) {
      east_to_west = new Connector(cent_to_east.getHead(),
         west_to_cent.getBase(),
         east_to_west_width);
   }
   privateBridge(
      *west_to_east,
      east_to_west,
      cent_to_east,
      compute_transpose,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   if (compute_transpose) {
      west_to_east->setTranspose(east_to_west, true);
   } else if (&west_to_east->getHead() == &west_to_east->getBase()) {
      west_to_east->setTranspose(west_to_east.get(), false);
   }

   d_object_timers->t_bridge->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   std::shared_ptr<Connector>& west_to_east,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   bool compute_transpose) const
{
   const tbox::SAMRAI_MPI& mpi =
      d_mpi.hasNullCommunicator() ? west_to_cent.getBase().getMPI() : d_mpi;

   if (d_barrier_before_communication) {
      mpi.Barrier();
   }
   d_object_timers->t_bridge->start();

   TBOX_ASSERT(west_to_cent.hasTranspose());
   TBOX_ASSERT(cent_to_east.hasTranspose());
   const tbox::Dimension& dim(cent_to_east.getConnectorWidth().getDim());
   const IntVector& zero_vector(IntVector::getZero(dim));
   const IntVector connector_width_limit(dim, -1); // No user-imposed limit.
   IntVector west_to_east_width(dim);
   IntVector east_to_west_width(dim);
   std::set<int> incoming_ranks, outgoing_ranks;
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   privateBridge_prologue(
      west_to_cent,
      cent_to_east,
      cent_to_east.getTranspose(),
      west_to_cent.getTranspose(),
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit,
      compute_transpose,
      west_to_east_width,
      east_to_west_width,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   Connector* east_to_west = nullptr;
   west_to_east.reset(new Connector(west_to_cent.getBase(),
      cent_to_east.getHead(),
      west_to_east_width));
   if (compute_transpose) {
      east_to_west = new Connector(cent_to_east.getHead(),
      west_to_cent.getBase(),
      east_to_west_width);
   }
   privateBridge(
      *west_to_east,
      east_to_west,
      cent_to_east,
      compute_transpose,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   if (compute_transpose) {
      west_to_east->setTranspose(east_to_west, true);
   } else if (&west_to_east->getHead() == &west_to_east->getBase()) {
      west_to_east->setTranspose(west_to_east.get(), false);
   }

   d_object_timers->t_bridge->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_cent,
   const Connector& cent_to_east,
   const IntVector& connector_width_limit) const
{
   const tbox::SAMRAI_MPI& mpi =
      d_mpi.hasNullCommunicator() ? cent_to_east.getBase().getMPI() : d_mpi;

   if (d_barrier_before_communication) {
      mpi.Barrier();
   }
   d_object_timers->t_bridge->start();

   TBOX_ASSERT(west_to_cent.hasTranspose());
   TBOX_ASSERT(cent_to_east.hasTranspose());
   Connector& cent_to_west = west_to_cent.getTranspose();
   const tbox::Dimension& dim(connector_width_limit.getDim());
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   IntVector west_to_east_width(dim);
   IntVector east_to_west_width(dim);
   std::set<int> incoming_ranks, outgoing_ranks;
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   bool compute_transpose = &west_to_cent != &west_to_cent.getTranspose();
   privateBridge_prologue(
      west_to_cent,
      cent_to_east,
      cent_to_east.getTranspose(),
      west_to_cent.getTranspose(),
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit,
      compute_transpose,
      west_to_east_width,
      east_to_west_width,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);
   west_to_cent.clearNeighborhoods();
   west_to_cent.setBase(west_to_cent.getBase());
   west_to_cent.setHead(cent_to_east.getHead());
   west_to_cent.setWidth(west_to_east_width, true);
   if (compute_transpose) {
      cent_to_west.clearNeighborhoods();
      cent_to_west.setBase(cent_to_east.getHead());
      cent_to_west.setHead(west_to_cent.getBase());
      cent_to_west.setWidth(east_to_west_width, true);
   }
   privateBridge(
      west_to_cent,
      &cent_to_west,
      cent_to_east,
      compute_transpose,
      incoming_ranks,
      outgoing_ranks,
      visible_west_nabrs,
      visible_east_nabrs);

   d_object_timers->t_bridge->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_prologue(
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west,
   bool west_nesting_is_known,
   const IntVector& cent_growth_to_nest_west,
   bool east_nesting_is_known,
   const IntVector& cent_growth_to_nest_east,
   const IntVector& connector_width_limit,
   bool compute_transpose,
   IntVector& west_to_east_width,
   IntVector& east_to_west_width,
   std::set<int>& incoming_ranks,
   std::set<int>& outgoing_ranks,
   NeighborSet& visible_west_nabrs,
   NeighborSet& visible_east_nabrs) const
{

   if (d_print_steps) {
      std::string dbgbord("bridge->  ");
      tbox::plog
      << "bridge west:\n" << west_to_cent.getBase().format(dbgbord, 3)
      << "bridge east:\n" << east_to_cent.getBase().format(dbgbord, 3)
      << "bridge center:\n" << cent_to_west.getBase().format(dbgbord, 3)
      << "bridge west_to_cent:\n" << west_to_cent.format(dbgbord, 3)
      << "bridge cent_to_west:\n" << cent_to_west.format(dbgbord, 3)
      << "bridge cent_to_east:\n" << cent_to_east.format(dbgbord, 3)
      << "bridge east_to_cent:\n" << east_to_cent.format(dbgbord, 3);
   }

   privateBridge_checkParameters(
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west);

   const BoxLevel& cent = cent_to_west.getBase();
   const BoxLevel& west = cent_to_west.getHead();
   const BoxLevel& east = cent_to_east.getHead();
   const IntVector& cent_refinement_ratio = cent.getRefinementRatio();
   const IntVector& west_refinement_ratio = west.getRefinementRatio();
   const IntVector& east_refinement_ratio = east.getRefinementRatio();

   const tbox::Dimension& dim(connector_width_limit.getDim());

   const IntVector& zero_vector(IntVector::getZero(dim));

   const IntVector finest_refinement_ratio =
      IntVector::max(
         cent_refinement_ratio,
         IntVector::max(west_refinement_ratio, east_refinement_ratio));

   const size_t num_blocks = cent_refinement_ratio.getNumBlocks();

   IntVector width_limit(connector_width_limit);
   if (width_limit.getNumBlocks() == 1 && num_blocks != 1) {
      if (width_limit.max() == width_limit.min()) {
         width_limit = IntVector(width_limit, num_blocks);
      } else {
         TBOX_ERROR("Anisotropic width limit argument for bridge must be of size equal to the number of blocks." << std::endl);
      }
   }

   /*
    * Using the bridge theorem, compute the largest bridge width for
    * which we can guarantee discovering all the overlaps (when
    * nesting is satisfied).  If either the east or west
    * BoxLevel's nesting in the center is known, compute the
    * output width by the bridge theorem, and use the bigger one.  If
    * neither is known, we assume that both east and west nest in
    * center, and just to do something reasonable.
    */
   IntVector output_width1(zero_vector, num_blocks);
   IntVector output_width2(zero_vector, num_blocks);
   if (west_nesting_is_known || east_nesting_is_known) {
      if (west_nesting_is_known) {
         output_width1 =
            cent_to_east.getConnectorWidth() - cent_growth_to_nest_west;
      }
      if (east_nesting_is_known) {
         output_width2 =
            cent_to_west.getConnectorWidth() - cent_growth_to_nest_east;
      }
      if (!(output_width1 >= zero_vector || output_width2 >= zero_vector)) {
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge_prologue:\n"
            << "Useless nesting specifications!\n"
            << "Neither west nor east BoxLevel nest with enough\n"
            << "margin to guarantee finding all overlaps.\n"
            << "To ensure you understand completness is not guaranteed,\n"
            << "this is considered an error.  To proceed anyway and live\n"
            << "with potential incompleteness, use a bridge interface\n"
            << "that does not claim any nesting.  Or, you can specify\n"
            << "a different nesting claim (but don't enable sanity\n"
            << "checking, which will catch your fib).\n");
      }
   } else {
      output_width1 = cent_to_east.getConnectorWidth();
      output_width2 = cent_to_west.getConnectorWidth();
   }
   IntVector output_width_in_finest_refinement_ratio =
      IntVector::max(output_width1, output_width2) *
      finest_refinement_ratio / cent_refinement_ratio;

   /*
    * Reduce the output width to the user-specified width limit.  Note
    * that the width limit is specified in the coarser of the east and
    * west refinement ratios.
    */
   if (width_limit >= zero_vector) {
      const IntVector coarser_refinement_ratio =
         IntVector::min(west_refinement_ratio, east_refinement_ratio);
      const IntVector width_limit_in_finest_refinement_ratio(
         width_limit * finest_refinement_ratio / coarser_refinement_ratio);
      if (!(width_limit_in_finest_refinement_ratio <= output_width_in_finest_refinement_ratio)) {
         /*
          * If user specifies a width limit, he is probably assuming
          * that the bridge's allowable width is no smaller.  If that
          * is not the case, this method will not crash, but it will
          * give bad results that result in elusive bugs.  Therefore,
          * we catch it immediately.
          */
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge_prologue input error:\n"
            << "The given connector width limit, " << width_limit
            << " (" << width_limit_in_finest_refinement_ratio
            << " in finest index space)\n"
            << "is not <= the maximum width of the bridge, "
            << output_width_in_finest_refinement_ratio
            << " (in finest index space).");
      }
      output_width_in_finest_refinement_ratio.min(
         width_limit_in_finest_refinement_ratio);
   }

   west_to_east_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / west_refinement_ratio);
   east_to_west_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / east_refinement_ratio);

   const int rank = cent.getMPI().getRank();

   /*
    * Owners we have to exchange information with are the ones
    * owning east/west Boxes visible to the local process.
    */
   cent_to_west.getLocalOwners(outgoing_ranks);
   west_to_cent.getLocalOwners(incoming_ranks);
   if (compute_transpose && &cent_to_west != &cent_to_east) {
      cent_to_east.getLocalOwners(outgoing_ranks);
      east_to_cent.getLocalOwners(incoming_ranks);
   }
   outgoing_ranks.erase(rank);
   incoming_ranks.erase(rank);

   /*
    * Create BoxContainers which will later be used to initialize the search
    * trees for visible east and west neighbors:
    * visible_west_nabrs and visible_east_nabrs.
    */
   d_object_timers->t_bridge_discover_get_neighbors->start();
   cent_to_west.getLocalNeighbors(visible_west_nabrs);
   cent_to_east.getLocalNeighbors(visible_east_nabrs);
   d_object_timers->t_bridge_discover_get_neighbors->stop();
}

/*
 ***********************************************************************
 *
 *                           west to east
 *        (west box_level)  ------------------> (east box_level)
 *                       ^  <------------------ ^
 *                        \    east to west    /
 *                         \                  /
 *           center to west \                / center to east
 *                           \              /
 *                            \            /
 *                          (center box_level)
 *
 * Bridge operation is in two phases, discovery and
 * sharing.  The discovery phase loops through local
 * Boxes in the center and comparing the west and east neighbors
 * for overlaps.  Local overlaps are stored immediately.
 * Remote overlaps are placed in messages to be sent to appropriate
 * processors by the sharing phase.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge(
   Connector& west_to_east,
   Connector* east_to_west,
   const Connector& cent_to_east,
   bool compute_transpose,
   const std::set<int>& incoming_ranks,
   const std::set<int>& outgoing_ranks,
   NeighborSet& visible_west_nabrs,
   NeighborSet& visible_east_nabrs) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (compute_transpose) {
      const IntVector& west_refinement_ratio =
         west_to_east.getBase().getRefinementRatio();
      const IntVector& east_refinement_ratio =
         west_to_east.getHead().getRefinementRatio();
      if (west_refinement_ratio / east_refinement_ratio * east_refinement_ratio ==
          west_refinement_ratio ||
          east_refinement_ratio / west_refinement_ratio * west_refinement_ratio ==
          east_refinement_ratio) {
         /*
          * If it's possible to make west<==>east transposes, it
          * should happen.  The requirement is that one refinement ratio is
          * an IntVector times the other.
          */
         TBOX_ASSERT(west_to_east.isTransposeOf(*east_to_west));
         TBOX_ASSERT(east_to_west->isTransposeOf(west_to_east));
      }
   }
#endif


   /*
    * Set up communication mechanism and post receives.
    * Note that in all_comms, all the incoming_comm come
    * first, the outgoing_comm later.
    */

   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int> * all_comms(0);

   d_object_timers->t_bridge_share->start();
   d_object_timers->t_bridge_setup_comm->start();

   if (d_sanity_check_method_preconditions) {
      if (!d_mpi.hasNullCommunicator() &&
          !d_mpi.isCongruentWith(west_to_east.getBase().getMPI())) {
         TBOX_ERROR("OverlapConnectorAlgorithm::findOverlaps_assumedPartition input error: Input\n"
            << "has SAMRAI_MPI that is incongruent with OverlapConnectorAlgorithm's.\n"
            << "See OverlapConnectorAlgorithm::setSAMRAI_MPI.\n");
      }
   }

   const tbox::SAMRAI_MPI& mpi =
      d_mpi.hasNullCommunicator() ? west_to_east.getBase().getMPI() : d_mpi;
   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Errant message detected.");
   }

   s_operation_mpi_tag = 0;

   setupCommunication(
      all_comms,
      comm_stage,
      mpi,
      incoming_ranks,
      outgoing_ranks,
      d_object_timers->t_bridge_MPI_wait,
      s_operation_mpi_tag,
      d_print_steps);

   d_object_timers->t_bridge_setup_comm->stop();
   d_object_timers->t_bridge_share->stop();

   /*
    * Messages for other processors describing removed and added relationships.
    */
   std::map<int, std::vector<int> > send_mesgs;
   for (std::set<int>::const_iterator itr(outgoing_ranks.begin());
        itr != outgoing_ranks.end(); ++itr) {
      send_mesgs[*itr];
   }

   /*
    * First step: Remove neighbor data for Boxes that are
    * going away and cache information to be sent out.
    */
   privateBridge_removeAndCache(
      send_mesgs,
      west_to_east,
      east_to_west,
      cent_to_east);

   privateBridge_discoverAndSend(
      send_mesgs,
      west_to_east,
      east_to_west,
      incoming_ranks,
      outgoing_ranks,
      all_comms,
      visible_west_nabrs,
      visible_east_nabrs);

   d_object_timers->t_bridge_share->start();

   receiveAndUnpack(
      west_to_east,
      east_to_west,
      incoming_ranks,
      all_comms,
      comm_stage,
      d_object_timers->t_bridge_receive_and_unpack,
      d_print_steps);

   d_object_timers->t_bridge_share->stop();

   if (all_comms) {
      delete[] all_comms;
   }

   if (d_sanity_check_method_postconditions) {
      west_to_east.assertConsistencyWithBase();
      west_to_east.assertConsistencyWithHead();
      if (compute_transpose) {
         east_to_west->assertConsistencyWithBase();
         east_to_west->assertConsistencyWithHead();
         east_to_west->assertTransposeCorrectness(west_to_east, true);
      }
   }

   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Errant message detected.");
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_checkParameters(
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west) const
{
   const BoxLevel& cent = cent_to_west.getBase();

   /*
    * Ensure that Connectors incident to and from the center agree on
    * what the center is.  This can be an expensive (though still
    * scalable) check, so we only check in debug mode, unless debugging,
    * in which it can be independently enabled in any mode.
    */
   if (cent != cent_to_east.getBase() ||
       cent != east_to_cent.getHead() ||
       cent != west_to_cent.getHead()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "Given Connectors to base and head of bridge are not incident\n"
         << "from the same center in\n"
         << "OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "west_to_cent is  TO  " << &west_to_cent.getHead() << "\n"
         << "cent_to_east is FROM " << &cent_to_east.getBase() << "\n"
         << "east_to_cent is  TO  " << &east_to_cent.getHead() << "\n"
         << "cent_to_west is FROM " << &cent_to_west.getBase() << "\n"
         );
   }
   /*
    * Ensure that head and base box_levels in argument agree with
    * head and base in the object.
    */
   if (cent_to_west.getHead() != west_to_cent.getBase()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "Given Connectors to and from base of bridge do not refer\n"
         << "to the base of the bridge in\n"
         << "OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "west_to_cent is FROM " << &west_to_cent.getBase() << "\n"
         << "cent_to_west is  TO  " << &cent_to_west.getHead() << "\n"
         );
   }
   if (cent_to_east.getHead() != east_to_cent.getBase()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "Given Connectors to and from head of bridge do not refer\n"
         << "to the head of the bridge in\n"
         << "OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "east_to_cent is FROM " << &east_to_cent.getBase() << "\n"
         << "cent_to_east is  TO  " << &cent_to_east.getHead() << "\n"
         );
   }
   if (!west_to_cent.isTransposeOf(cent_to_west)) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "Given Connectors between base and center of bridge\n"
         << "are not transposes of each other.\n"
         << "See OverlapConnectorAlgorithm::isTransposeOf().\n"
         );
   }
   if (!east_to_cent.isTransposeOf(cent_to_east)) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::privateBridge_checkParameters:\n"
         << "Given Connectors between head and center of bridge\n"
         << "are not transposes of each other.\n"
         << "See OverlapConnectorAlgorithm::isTransposeOf().\n"
         );
   }

   // Expensive sanity checks:
   if (d_sanity_check_method_preconditions) {
      west_to_cent.assertConsistencyWithBase();
      west_to_cent.assertConsistencyWithHead();
      cent_to_east.assertConsistencyWithBase();
      cent_to_east.assertConsistencyWithHead();
      east_to_cent.assertConsistencyWithBase();
      east_to_cent.assertConsistencyWithHead();
      cent_to_west.assertConsistencyWithBase();
      cent_to_west.assertConsistencyWithHead();
      west_to_cent.assertTransposeCorrectness(cent_to_west);
      cent_to_west.assertTransposeCorrectness(west_to_cent);
      east_to_cent.assertTransposeCorrectness(cent_to_east);
      cent_to_east.assertTransposeCorrectness(east_to_cent);
   }
}

/*
 ***********************************************************************
 * Remove relationships from resulting overlap.  Cache outgoing
 * information in message buffers.
 *
 * TODO: This method is a no-op.  Is it a place holder for something?
 * Maybe it should be removed.  artf19532.
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_removeAndCache(
   std::map<int, std::vector<int> >& send_mesgs,
   Connector& overlap_connector,
   Connector* overlap_connector_transpose,
   const Connector& misc_connector) const
{
   d_object_timers->t_bridge_remove_and_cache->start();

   NULL_USE(send_mesgs);
   NULL_USE(overlap_connector);
   NULL_USE(overlap_connector_transpose);
   NULL_USE(misc_connector);
   /*
    * As the overlap relationships are empty to start there are never any
    * that need to be deleted.
    */
   d_object_timers->t_bridge_remove_and_cache->stop();
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_discoverAndSend(
   std::map<int, std::vector<int> >& send_mesgs,
   Connector& west_to_east,
   Connector* east_to_west,
   const std::set<int>& incoming_ranks,
   const std::set<int>& outgoing_ranks,
   tbox::AsyncCommPeer<int>* all_comms,
   NeighborSet& visible_west_nabrs,
   NeighborSet& visible_east_nabrs) const
{
   if (visible_west_nabrs.empty() && visible_east_nabrs.empty()) {
      return;
   }

   /*
    * Discover overlaps.  Overlaps are either locally stored or
    * packed into a message for sending.
    */

   d_object_timers->t_bridge_discover_and_send->start();

   if (d_print_steps) {
      tbox::plog << "Before building RBBTs:\n"
                 << "visible_west_nabrs:"
                 << visible_west_nabrs.format("\n  ")
                 << "visible_east_nabrs:"
                 << visible_east_nabrs.format("\n  ");
      tbox::plog << "\nincoming ranks (" << incoming_ranks.size() << "): ";
      copy(incoming_ranks.begin(), incoming_ranks.end(), std::ostream_iterator<int>(tbox::plog, " "));
      tbox::plog << "\noutgoing ranks (" << outgoing_ranks.size() << "): ";
      copy(outgoing_ranks.begin(), outgoing_ranks.end(), std::ostream_iterator<int>(tbox::plog, " "));
      tbox::plog << std::endl;
   }

   bool compute_transpose =
      (east_to_west != 0 && east_to_west != &west_to_east);

   const BoxLevel& east(west_to_east.getBase());
   const std::shared_ptr<const BaseGridGeometry>& grid_geometry(
      east.getGridGeometry());

   const tbox::SAMRAI_MPI& mpi = d_mpi.hasNullCommunicator() ? east.getMPI() : d_mpi;

   const tbox::Dimension& dim(east.getDim());
   const int rank = mpi.getRank();

   d_object_timers->t_bridge_discover_form_rbbt->start();
   const BoxContainer east_rbbt(visible_east_nabrs);
   east_rbbt.makeTree(grid_geometry.get());
   // Note: west_rbbt only needed when compute_transpose is true.
   BoxContainer empty_nabrs(true);
   const BoxContainer west_rbbt(
      compute_transpose ? visible_west_nabrs : empty_nabrs);
   west_rbbt.makeTree(grid_geometry.get());
   d_object_timers->t_bridge_discover_form_rbbt->stop();

   /*
    * Local process can find some neighbors for the (local and
    * remote) Boxes in visible_west_nabrs and visible_east_nabrs.
    * Separate this into 2 parts: discovery of remote Boxes which
    * may be threaded, and discovery of local Boxes can not be.
    * In either case we loop through the visible_west_nabrs and
    * compare each to visible_east_nabrs, looking for overlaps.
    * Then vice versa.  Since each of these NeighborSets is
    * ordered by processor owner first and we know each non-local
    * processor we can construct each non-local message in a
    * separate thread and then find and set all the local overlaps.
    *
    * To do this we first separate visible_west_nabrs into 2 groups
    * non-local and local neighbors.  Also do the same for
    * visible_east_nabrs.
    */
   bool ordered = true;
   NeighborSet visible_local_west_nabrs(ordered);
   NeighborSet visible_local_east_nabrs(ordered);
   const Box this_proc_start(dim, GlobalId(LocalId::getZero(), rank));
   NeighborSet::iterator west_ni =
      visible_west_nabrs.lowerBound(this_proc_start);
   NeighborSet::iterator east_ni =
      visible_east_nabrs.lowerBound(this_proc_start);
   while (west_ni != visible_west_nabrs.end() &&
          west_ni->getOwnerRank() == rank) {
      visible_local_west_nabrs.insert(*west_ni);
      visible_west_nabrs.erase(west_ni++);
   }
   while (east_ni != visible_east_nabrs.end() &&
          east_ni->getOwnerRank() == rank) {
      visible_local_east_nabrs.insert(*east_ni);
      visible_east_nabrs.erase(east_ni++);
   }

   // Discover all non-local overlaps.
   int i = 0;
   int imax = static_cast<int>(outgoing_ranks.size());
   std::vector<int> another_outgoing_ranks(outgoing_ranks.size());
   for (std::set<int>::const_iterator outgoing_ranks_itr(outgoing_ranks.begin());
        outgoing_ranks_itr != outgoing_ranks.end(); ++outgoing_ranks_itr) {
      another_outgoing_ranks[i++] = *outgoing_ranks_itr;
   }
#ifdef HAVE_OPENMP
#pragma omp parallel private(i) num_threads(4)
   {
#pragma omp for schedule(dynamic) nowait
#endif
   for (i = 0; i < imax; ++i) {
      BoxId outgoing_proc_start_id(
         LocalId::getZero(),
         another_outgoing_ranks[i]);
      Box outgoing_proc_start(dim, outgoing_proc_start_id);
      NeighborSet::const_iterator thread_west_ni =
         visible_west_nabrs.lowerBound(outgoing_proc_start);
      NeighborSet::const_iterator thread_east_ni =
         visible_east_nabrs.lowerBound(outgoing_proc_start);
      privateBridge_discover(
         send_mesgs[another_outgoing_ranks[i]],
         west_to_east,
         east_to_west,
         visible_west_nabrs,
         visible_east_nabrs,
         thread_west_ni,
         thread_east_ni,
         another_outgoing_ranks[i],
         east_rbbt,
         west_rbbt,
         dim,
         compute_transpose,
         rank);
   }
#ifdef HAVE_OPENMP
}
#endif

   /*
    * Send all non-local overlap messages.
    * As an optimization, send to the next higher ranked process first followed
    * by successively higher processes and finally looping around to process
    * 0 through the next lower ranked process.  This spreads out the sends more
    * evenly and prevents everyone from sending to the same processor (like
    * process 0) at the same time.
    */
   int num_outgoing_ranks = static_cast<int>(outgoing_ranks.size());
   int num_incoming_ranks = static_cast<int>(incoming_ranks.size());
   int num_comms = num_outgoing_ranks + num_incoming_ranks;
   std::set<int>::const_iterator outgoing_ranks_itr(
      outgoing_ranks.lower_bound(rank + 1));
   if (outgoing_ranks_itr == outgoing_ranks.end()) {
      outgoing_ranks_itr = outgoing_ranks.begin();
   }
   int comm_offset = num_incoming_ranks;
   for ( ; comm_offset < num_comms; ++comm_offset) {
      if (all_comms[comm_offset].getPeerRank() == *outgoing_ranks_itr) {
         break;
      }
   }
   TBOX_ASSERT(num_outgoing_ranks == 0 || comm_offset < num_comms);
   for (int outgoing_ranks_ctr = 0;
        outgoing_ranks_ctr < num_outgoing_ranks; ++outgoing_ranks_ctr) {
      std::vector<int>& send_mesg = send_mesgs[*outgoing_ranks_itr];
      tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[comm_offset];
      outgoing_comm.beginSend(
         &send_mesg[0],
         send_mesg.size());
      ++comm_offset;
      ++outgoing_ranks_itr;
      TBOX_ASSERT((outgoing_ranks_itr == outgoing_ranks.end()) ==
         (comm_offset == num_comms));
      if (outgoing_ranks_itr == outgoing_ranks.end()) {
         outgoing_ranks_itr = outgoing_ranks.begin();
      }
      if (comm_offset == num_comms) {
         comm_offset = num_incoming_ranks;
      }
      if (s_print_steps == 'y') {
         tbox::plog << "Sent to " << outgoing_comm.getPeerRank() << std::endl;
      }
   }

   // Discover all local overlaps and store them in the Connector(s).
   NeighborSet::const_iterator west_local_ni =
      visible_local_west_nabrs.lowerBound(this_proc_start);
   NeighborSet::const_iterator east_local_ni =
      visible_local_east_nabrs.lowerBound(this_proc_start);
   privateBridge_discover(
      send_mesgs[rank],
      west_to_east,
      east_to_west,
      visible_local_west_nabrs,
      visible_local_east_nabrs,
      west_local_ni,
      east_local_ni,
      rank,
      east_rbbt,
      west_rbbt,
      dim,
      compute_transpose,
      rank);

   d_object_timers->t_bridge_discover_and_send->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_discover(
   std::vector<int>& send_mesg,
   Connector& west_to_east,
   Connector* east_to_west,
   const NeighborSet& visible_west_nabrs,
   const NeighborSet& visible_east_nabrs,
   NeighborSet::const_iterator& west_ni,
   NeighborSet::const_iterator& east_ni,
   int curr_owner,
   const BoxContainer& east_rbbt,
   const BoxContainer& west_rbbt,
   const tbox::Dimension& dim,
   bool compute_transpose,
   int rank) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   // Owners that we have sent messages to.  Used for debugging.
   std::set<int> owners_sent_to;
#endif

   TBOX_ASSERT(owners_sent_to.find(curr_owner) == owners_sent_to.end());

   /*
    * Set up send_message to contain info discovered
    * locally but needed by curr_owner.
    *
    * Content of send_mesg:
    * - neighbor-removal section cached in send_mesg.
    * - offset to the reference section (see below)
    * - number of west boxes for which neighbors are found
    * - number of east boxes for which neighbors are found
    *   - id of west/east box
    *   - number of neighbors found for west/east box.
    *     - BoxId of neighbors found.
    *       Boxes of these found neighbors are given in the
    *       reference section of the message.
    * - reference section: all the Boxes referenced as
    *   neighbors (accumulated in referenced_west_nabrs
    *   and referenced_east_nabrs).
    *   - number of referenced west neighbors
    *   - number of referenced east neighbors
    *   - referenced west neighbors
    *   - referenced east neighbors
    *
    * The purpose of factoring out info on the neighbors referenced
    * is to reduce redundant data that can eat up lots of memory
    * and message passing bandwidth when there are lots of Boxes
    * with the same neighbors.
    */

   /*
    * The first section of the send_mesg is the remote neighbor-removal
    * section (computed above).
    */
   if (curr_owner != rank && send_mesg.empty()) {
      // No neighbor-removal data found for curr_owner.
      send_mesg.insert(send_mesg.end(), 0);
   }

   // Indices of certain positions in send_mesg.
   const int idx_offset_to_ref = static_cast<int>(send_mesg.size());
   const int idx_num_west_boxes = idx_offset_to_ref + 1;
   const int idx_num_east_boxes = idx_offset_to_ref + 2;
   send_mesg.insert(send_mesg.end(), 3, 0);

   // Boxes referenced in the message, used when adding ref section.
   BoxContainer referenced_west_nabrs;
   BoxContainer referenced_east_nabrs;

   if (s_print_steps == 'y') {
      tbox::plog << "Finding west --> east overlaps for owner "
                 << curr_owner << std::endl;
   }

   // Find neighbors for all west boxes owned by curr_owner.
   privateBridge_findOverlapsForOneProcess(
      curr_owner,
      visible_west_nabrs,
      west_ni,
      send_mesg,
      idx_num_west_boxes,
      west_to_east,
      referenced_east_nabrs,
      east_rbbt);

   // Find neighbors for all east boxes owned by curr_owner.
   if (compute_transpose) {
      if (s_print_steps == 'y') {
         tbox::plog << "Finding west <-- east overlaps for owner "
                    << curr_owner << std::endl;
      }
      privateBridge_findOverlapsForOneProcess(
         curr_owner,
         visible_east_nabrs,
         east_ni,
         send_mesg,
         idx_num_east_boxes,
         *east_to_west,
         referenced_west_nabrs,
         west_rbbt);
   }

   if (curr_owner != rank) {
      /*
       * If this discovery is off processor then the send message must be
       * filled with the referenced neighbors.
       */

      packReferencedNeighbors(
         send_mesg,
         idx_offset_to_ref,
         referenced_east_nabrs,
         referenced_west_nabrs,
         dim,
         s_print_steps == 'y');

#ifdef DEBUG_CHECK_ASSERTIONS
      owners_sent_to.insert(curr_owner);
#endif

   } // Block to send discoveries to curr_owner.
}

/*
 ***********************************************************************
 *
 * Find overlaps from visible_base_nabrs to head_rbbt.  Find only
 * overlaps for Boxes owned by owner_rank.
 *
 * On entry, base_ni points to the first Box in visible_base_nabrs
 * owned by owner_rank.  Increment base_ni past those Boxes
 * processed and remove them from visible_base_nabrs.
 *
 * Save local and semilocal overlaps in bridging_connector.  For
 * remote overlaps, pack in send_mesg, add head Box to
 * referenced_head_nabrs and increment
 * send_mesg[remote_box_counter_index].
 *
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge_findOverlapsForOneProcess(
   const int owner_rank,
   const NeighborSet& visible_base_nabrs,
   NeighborSet::const_iterator& base_ni,
   std::vector<int>& send_mesg,
   const int remote_box_counter_index,
   Connector& bridging_connector,
   NeighborSet& referenced_head_nabrs,
   const BoxContainer& head_rbbt) const
{
   const IntVector& head_refinement_ratio(
      bridging_connector.getHead().getRefinementRatio());

   const BaseGridGeometry& grid_geom = 
      *bridging_connector.getBase().getGridGeometry();

   bool refine_base = false;
   bool coarsen_base = false;
   if (bridging_connector.getHead().getRefinementRatio() ==
       bridging_connector.getBase().getRefinementRatio()) {
      /*
       * Don't do any coarsen/refine because head and base have same
       * refinement ratio.
       */
   } else if (bridging_connector.getHead().getRefinementRatio() <=
              bridging_connector.getBase().getRefinementRatio()) {
      coarsen_base = true;
   } else if (bridging_connector.getHead().getRefinementRatio() >=
              bridging_connector.getBase().getRefinementRatio()) {
      refine_base = true;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   else {
      TBOX_ERROR("Can't coarsen in one direction and refine in another");
   }
#endif

   const PeriodicShiftCatalog& shift_catalog =
      bridging_connector.getHead().getGridGeometry()->getPeriodicShiftCatalog(); 
   // Should be made a member to avoid repetitive alloc/dealloc.
   // Reserve in privateBridge and used here.
   BoxContainer found_nabrs, scratch_found_nabrs;

   while (base_ni != visible_base_nabrs.end() &&
          base_ni->getOwnerRank() == owner_rank) {
      const Box& visible_base_nabrs_box = *base_ni;
      if (d_print_steps) {
         tbox::plog << "Finding neighbors for non-periodic visible_base_nabrs_box "
                    << visible_base_nabrs_box << std::endl;
      }
      BoxContainer grown_boxes;
      if (grid_geom.getNumberBlocks() == 1 || grid_geom.hasIsotropicRatios()) {
         Box base_box = visible_base_nabrs_box;
         base_box.grow(bridging_connector.getConnectorWidth());
         if (refine_base) {
            base_box.refine(bridging_connector.getRatio());
         }
         else if (coarsen_base) {
            base_box.coarsen(bridging_connector.getRatio());
         }
         grown_boxes.pushBack(base_box);
      } else {
         BoxUtilities::growAndAdjustAcrossBlockBoundary(
            grown_boxes,
            visible_base_nabrs_box,
            bridging_connector.getBase().getGridGeometry(),
            bridging_connector.getBase().getRefinementRatio(),
            bridging_connector.getRatio(),
            bridging_connector.getConnectorWidth(),
            refine_base,
            coarsen_base);
      }

      found_nabrs.clear();
      for (BoxContainer::iterator g_itr = grown_boxes.begin();
           g_itr != grown_boxes.end(); ++g_itr) {

         head_rbbt.findOverlapBoxes(found_nabrs, *g_itr,
                                    head_refinement_ratio,
                                    true /* include singularity block neighbors */ );
      }
      if (d_print_steps) {
         tbox::plog << "Found " << found_nabrs.size() << " neighbors:";
         found_nabrs.print(tbox::plog);
         //BoxContainerUtils::recursivePrintBoxVector(found_nabrs, tbox::plog, "\n ");
         tbox::plog << std::endl;
      }
      if (!found_nabrs.empty()) {
         if (visible_base_nabrs_box.isPeriodicImage()) {
            privateBridge_unshiftOverlappingNeighbors(
               visible_base_nabrs_box,
               found_nabrs,
               scratch_found_nabrs,
               bridging_connector.getHead().getRefinementRatio(),
               shift_catalog);
         }
         if (owner_rank != bridging_connector.getMPI().getRank()) {
            // Pack up info for sending.
            ++send_mesg[remote_box_counter_index];
            const int subsize = 3
               + BoxId::commBufferSize() * static_cast<int>(found_nabrs.size());
            send_mesg.insert(send_mesg.end(), subsize, -1);
            int* submesg = &send_mesg[send_mesg.size() - subsize];
            *(submesg++) = visible_base_nabrs_box.getLocalId().getValue();
            *(submesg++) = static_cast<int>(
               visible_base_nabrs_box.getBlockId().getBlockValue());
            *(submesg++) = static_cast<int>(found_nabrs.size());
            for (BoxContainer::const_iterator na = found_nabrs.begin();
                 na != found_nabrs.end(); ++na) {
               const Box& head_nabr = *na;
               referenced_head_nabrs.insert(head_nabr);
               head_nabr.getBoxId().putToIntBuffer(submesg);
               submesg += BoxId::commBufferSize();
            }
         } else {
            // Save neighbor info locally.
            BoxId unshifted_base_box_id;
            if (!visible_base_nabrs_box.isPeriodicImage()) {
               unshifted_base_box_id = visible_base_nabrs_box.getBoxId();
            } else {
               unshifted_base_box_id.initialize(
                  visible_base_nabrs_box.getLocalId(),
                  visible_base_nabrs_box.getOwnerRank(),
                  PeriodicId::zero());
            }
            // Add found neighbors for visible_base_nabrs_box.
            if (!found_nabrs.empty()) {
               Connector::NeighborhoodIterator base_box_itr =
                  bridging_connector.makeEmptyLocalNeighborhood(
                     unshifted_base_box_id);
               for (BoxContainer::const_iterator na = found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  bridging_connector.insertLocalNeighbor(*na, base_box_itr);
               }
            }
         }
      }
      if (d_print_steps) {
         tbox::plog << "Erasing visible base nabr " << (*base_ni) << std::endl;
      }
      ++base_ni;
      if (s_print_steps == 'y') {
         if (base_ni == visible_base_nabrs.end()) {
            tbox::plog << "Next base nabr: end" << std::endl;
         } else {
            tbox::plog << "Next base nabr: " << *base_ni << std::endl;
         }
      }

   }
}

/*
 ***********************************************************************
 * Shift neighbors by amount equal and opposite of a Box's shift so that
 * they become neighbors of the unshifed box.  If this results in a
 * neighbor shift that is not in the shift catalog, discard the neighbor.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge_unshiftOverlappingNeighbors(
   const Box& box,
   BoxContainer& neighbors,
   BoxContainer& scratch_space,
   const IntVector& neighbor_refinement_ratio,
   const PeriodicShiftCatalog& shift_catalog) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, neighbor_refinement_ratio);

   scratch_space.clear();
//   scratch_space.reserve(neighbors.size());
   for (BoxContainer::iterator na = neighbors.begin();
        na != neighbors.end(); ++na) {
      Box& nabr = *na;
      IntVector sum_shift =
         shift_catalog.shiftNumberToShiftDistance(nabr.getPeriodicId())
         - shift_catalog.shiftNumberToShiftDistance(box.getPeriodicId());
      const PeriodicId new_shift_number =
         shift_catalog.shiftDistanceToShiftNumber(sum_shift);
      if (new_shift_number.getPeriodicValue() !=
          shift_catalog.getInvalidShiftNumber()) {
         nabr.initialize(nabr,
                         new_shift_number,
                         neighbor_refinement_ratio,
                         shift_catalog);
         scratch_space.pushBack(nabr);
      }
   }
   if (scratch_space.size() != neighbors.size()) {
      // We have discarded some neighbors due to invalid shift.
      neighbors.swap(scratch_space);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::initializeCallback()
{
   // Initialize timers with default prefix.
   getAllTimers(s_default_timer_prefix,
      s_static_timers[s_default_timer_prefix]);
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */

void
OverlapConnectorAlgorithm::finalizeCallback()
{
   s_static_timers.clear();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::string timer_prefix_used;
   if (s_ignore_external_timer_prefix == 'y') {
      timer_prefix_used = s_default_timer_prefix;
   } else {
      timer_prefix_used = timer_prefix;
   }
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix_used));
   if (ti == s_static_timers.end()) {
      d_object_timers = &s_static_timers[timer_prefix_used];
      getAllTimers(timer_prefix_used, *d_object_timers);
   } else {
      d_object_timers = &(ti->second);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_find_overlaps_rbbt = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::findOverlaps_rbbt()");

   timers.t_find_overlaps_assumed_partition = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::findOverlaps_assumedPartition()");
   timers.t_find_overlaps_assumed_partition_connect_to_ap = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::findOverlaps_assumedPartition()_connect_to_ap");
   timers.t_find_overlaps_assumed_partition_get_ap = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::findOverlaps_assumedPartition()_get_ap");
   timers.t_find_overlaps_assumed_partition_transpose = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::findOverlaps_assumedPartition()_transpose");

   timers.t_bridge = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateBridge()");
   timers.t_bridge_setup_comm = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::setupCommunication()");
   timers.t_bridge_remove_and_cache = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateBridge_removeAndCache()");
   timers.t_bridge_discover_and_send = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateBridge()_discoverAndSend");
   timers.t_bridge_discover_get_neighbors = tbox::TimerManager::getManager()->
      getTimer(
         timer_prefix + "::privateBridge()_discover_get_neighbors");
   timers.t_bridge_discover_form_rbbt = tbox::TimerManager::getManager()->
      getTimer(
         timer_prefix + "::privateBridge()_discover_form_rbbt");
   timers.t_bridge_share = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateBridge()_share");
   timers.t_bridge_receive_and_unpack = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::receiveAndUnpack");
   timers.t_bridge_MPI_wait = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateBridge()_MPI_wait");
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
