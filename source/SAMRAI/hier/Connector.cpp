/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of edges incident from a box_level of a distributed
 *                box graph.
 *
 ************************************************************************/
#include "SAMRAI/hier/Connector.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/ConnectorStatistics.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/CenteredRankTree.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <map>
#include <vector>
#include <set>
#include <algorithm>
//#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int Connector::HIER_CONNECTOR_VERSION = 0;

std::shared_ptr<tbox::Timer> Connector::t_acquire_remote_relationships;
std::shared_ptr<tbox::Timer> Connector::t_cache_global_reduced_data;
std::shared_ptr<tbox::Timer> Connector::t_find_overlaps_rbbt;

tbox::StartupShutdownManager::Handler
Connector::s_initialize_finalize_handler(
   Connector::initializeCallback,
   0,
   0,
   Connector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const tbox::Dimension& dim):
   d_base_handle(),
   d_head_handle(),
   d_base_width(dim),
   d_ratio(dim),
   d_ratio_is_exact(false),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_mpi(MPI_COMM_NULL),
   d_parallel_state(BoxLevel::DISTRIBUTED),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_transpose(0),
   d_owns_transpose(false)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const tbox::Dimension& dim,
   tbox::Database& restart_db):
   d_base_handle(),
   d_head_handle(),
   d_base_width(dim),
   d_ratio(dim),
   d_ratio_is_exact(false),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_mpi(MPI_COMM_NULL),
   d_parallel_state(BoxLevel::DISTRIBUTED),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_transpose(0),
   d_owns_transpose(false)
{
   getFromRestart(restart_db);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const Connector& other):
   d_base_handle(other.d_base_handle),
   d_head_handle(other.d_head_handle),
   d_base_width(other.d_base_width),
   d_ratio(other.d_ratio),
   d_ratio_is_exact(other.d_ratio_is_exact),
   d_head_coarser(other.d_head_coarser),
   d_relationships(other.d_relationships),
   d_global_relationships(other.d_global_relationships),
   d_mpi(other.d_mpi),
   d_parallel_state(other.d_parallel_state),
   d_finalized(other.d_finalized),
   d_global_number_of_neighbor_sets(other.d_global_number_of_neighbor_sets),
   d_global_number_of_relationships(other.d_global_number_of_relationships),
   d_global_data_up_to_date(other.d_global_data_up_to_date),
   d_transpose(other.d_transpose),
   d_owns_transpose(false)
{
   size_t num_blocks = 
      d_base_handle->getBoxLevel().getGridGeometry()->getNumberBlocks();

   if (d_base_width.getNumBlocks() == 1 && num_blocks != 1) {
      if (d_base_width.max() == d_base_width.min()) {
         d_base_width = IntVector(d_base_width, num_blocks);
      } else {
         TBOX_ERROR("Connector::Connector: anisotropic connector\n"
            << "width " << d_base_width << " must be \n"
            << "defined for " << num_blocks << " blocks.\n");
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const BoxLevel& base_box_level,
   const BoxLevel& head_box_level,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state):
   d_base_width(IntVector::getZero(base_width.getDim())),
   d_ratio(IntVector::getZero(base_width.getDim())),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_mpi(base_box_level.getMPI()),
   d_parallel_state(parallel_state),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(true),
   d_transpose(0),
   d_owns_transpose(false)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(base_box_level,
      head_box_level,
      base_width);
   IntVector tmp_base_width(base_width);
   size_t num_blocks = base_box_level.getGridGeometry()->getNumberBlocks();
   if (tmp_base_width.getNumBlocks() == 1 && num_blocks != 1) {
      if (tmp_base_width.max() == tmp_base_width.min()) {
         tmp_base_width = IntVector(tmp_base_width, num_blocks);
      } else {
         TBOX_ERROR("Connector::Connector: anisotropic connector\n"
            << "width " << base_width << " must be \n"
            << "defined for " << num_blocks << " blocks.\n");
      }
   }

   setBase(base_box_level);
   setHead(head_box_level);
   setWidth(tmp_base_width, true);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::~Connector()
{
   clear();
   if (d_transpose && d_owns_transpose && d_transpose != this) {
      delete d_transpose;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
Connector&
Connector::operator = (
   const Connector& rhs)
{
   if (this != &rhs) {
      d_base_handle = rhs.d_base_handle;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;
      d_global_number_of_neighbor_sets = rhs.d_global_number_of_neighbor_sets;
      d_head_handle = rhs.d_head_handle;
      d_global_number_of_relationships = rhs.d_global_number_of_relationships;
      d_relationships = rhs.d_relationships;
      d_global_relationships = rhs.d_global_relationships;
      d_mpi = rhs.d_mpi;
      d_base_width = rhs.d_base_width;
      d_ratio = rhs.d_ratio;
      d_head_coarser = rhs.d_head_coarser;
      d_parallel_state = rhs.d_parallel_state;
      d_finalized = rhs.d_finalized;
      d_transpose = rhs.d_transpose; // TODO: This leads to a memory error.
      d_owns_transpose = false;
   }
   return *this;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
Connector::operator == (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return true;
   }
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return true;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return false;
   }
   if (isFinalized() && !rhs.isFinalized()) {
      return false;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return false;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return false;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return false;
   }
   if (d_relationships != rhs.d_relationships) {
      return false;
   }

   return true;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
Connector::operator != (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return false;
   }
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return false;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return true;
   }
   if (isFinalized() && !rhs.isFinalized()) {
      return true;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return true;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return true;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return true;
   }
   if (d_relationships != rhs.d_relationships) {
      return true;
   }

   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::insertNeighbors(
   const BoxContainer& neighbors,
   const BoxId& base_box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       base_box.getOwnerRank() != getMPI().getRank()) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(base_box)) {
      TBOX_ERROR(
         "Exiting due to above reported error."
         << "Connector::insertNeighbors: Cannot access neighbors for\n"
         << "id " << base_box << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships.insert(base_box, neighbors);
   }
   if (base_box.getOwnerRank() == getMPI().getRank()) {
      d_relationships.insert(base_box, neighbors);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::eraseNeighbor(
   const Box& neighbor,
   const BoxId& box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       box_id.getOwnerRank() != getMPI().getRank()) {
      TBOX_ERROR("Connector::eraseNeighbor error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(box_id)) {
      TBOX_ERROR(
         "Connector::eraseNeighbors: Cannot access neighbors for\n"
         << "id " << box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships.erase(box_id, neighbor);
   }
   if (box_id.getOwnerRank() == getMPI().getRank()) {
      d_relationships.erase(box_id, neighbor);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::shrinkWidth(
   const IntVector& new_width)
{
   IntVector shrink_width(new_width);
   if (shrink_width.getNumBlocks() == 1 &&
       d_base_width.getNumBlocks() != 1) {
      if (shrink_width.max() == shrink_width.min()) {
         size_t new_size = d_base_width.getNumBlocks();
         shrink_width = IntVector(shrink_width, new_size);
      } else {
         TBOX_ERROR("Connector::shrinkWidth: new anisotropic connector\n"
            << "width " << shrink_width << " must be \n"
            << "defined for " << d_base_width.getNumBlocks() << " blocks.\n");
      }
   }

   TBOX_ASSERT(shrink_width.getNumBlocks() == d_base_width.getNumBlocks()); 
   if (!(shrink_width <= getConnectorWidth())) {
      TBOX_ERROR("Connector::shrinkWidth: new connector\n"
         << "width " << shrink_width << " involves an\n"
         << "enlargement of the current width "
         << getConnectorWidth());
   }
   else if (shrink_width == getConnectorWidth()) {
      // This is a no-op.
      return;
   }

   // Have not yet written this for GLOBALIZED mode.
   TBOX_ASSERT(getParallelState() == BoxLevel::DISTRIBUTED);

   /*
    * Remove overlaps that disappeared given the new GCW.
    * Swap out the overlaps, modify them then swap them back in.
    */

   const bool head_coarser = getHeadCoarserFlag();
   const bool base_coarser = !getHeadCoarserFlag() &&
      getBase().getRefinementRatio() != getHead().getRefinementRatio();

   const std::shared_ptr<const BaseGridGeometry>& grid_geom(
      getBase().getGridGeometry());

   if (grid_geom->getNumberBlocks() == 1 || grid_geom->hasIsotropicRatios()) {
      for (NeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& box_id = *ei;
         const Box& box = *getBase().getBoxStrict(box_id);
         Box box_box = box;
         box_box.grow(shrink_width);
         if (base_coarser) {
            box_box.refine(getRatio());
         }
         for (NeighborIterator na = begin(ei);
              na != end(ei); /* incremented in loop */) {
            const Box& nabr = *na;
            Box nabr_box = nabr;
            if (nabr.getBlockId() != box.getBlockId()) {
               grid_geom->transformBox(nabr_box,
                  getHead().getRefinementRatio(),
                  box.getBlockId(),
                  nabr.getBlockId());
            }
            if (head_coarser) {
               nabr_box.refine(getRatio());
            }
            ++na;
            if (!box_box.intersects(nabr_box)) {
               d_relationships.erase(ei, nabr);
            }
         }
      }
   } else {
      for (NeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& box_id = *ei;
         const Box& box = *getBase().getBoxStrict(box_id);
         BoxContainer grown_boxes;
         BoxUtilities::growAndAdjustAcrossBlockBoundary(grown_boxes,
            box,
            grid_geom,
            getBase().getRefinementRatio(),
            getRatio(),
            shrink_width,
            base_coarser,
            head_coarser);

         for (NeighborIterator na = begin(ei);
              na != end(ei); /* incremented in loop */) {

            const Box& nabr = *na;
            bool intersection = false;

            for (BoxContainer::iterator b_itr = grown_boxes.begin();
                 b_itr != grown_boxes.end(); ++b_itr) {

               if (nabr.getBlockId() == b_itr->getBlockId()) {
                  Box nabr_box = nabr;

                  if (head_coarser) {
                     nabr_box.refine(getRatio());
                  }
                  if (b_itr->intersects(nabr_box)) {
                     intersection = true;
                     break;
                  }
               }
            }
            ++na;
            if (!intersection) {
               d_relationships.erase(ei, nabr);
            }
         }
      }
   }

   d_base_width = shrink_width;
   return;
}

/*
 ***********************************************************************
 * Set this Connector to the transpose of another Connector.  The
 * other's base owners tell its head owners about relationships, and
 * head owners (base owners of this Connector) populate the this
 * Connector.
 *
 * This method uses the termination message technique of the assumed
 * partition algorithm.
 *
 * Two communication patterns are executed simultaneously: edge info
 * and termination messages.  For edge info, other's base box owners
 * send data to this's base box owners, who respond with
 * acknowledgement messages.  Termination messages let the processes
 * know when to stop checking for edge messages.  They are propagated
 * up then down a rank tree.  Upward messages inform processes that
 * their descendents have received all needed acknowledgements.
 * Downward messages inform processes that the entire tree completed
 * its acknowledgements, indicating that there are no messages in
 * transit and the process can stop.
 ***********************************************************************
 */
void
Connector::computeTransposeOf(const Connector& other,
                              const tbox::SAMRAI_MPI& mpi)
{
   *this = Connector(other.getHead(), other.getBase(),
         convertHeadWidthToBase(other.getHead().getRefinementRatio(),
            other.getBase().getRefinementRatio(),
            other.getConnectorWidth()));

   const tbox::SAMRAI_MPI& mpi1 = mpi.hasNullCommunicator() ? getBase().getMPI() : mpi;

   // Order locally visible edges by owners who need to know about them.
   typedef std::map<Box, BoxContainer, Box::id_less> FullNeighborhoodSet;
   FullNeighborhoodSet reordered_relationships;
   other.reorderRelationshipsByHead(reordered_relationships);

   /*
    * We receive different types of messages from different sources
    * without knowing which one is next, so we use the same MPI tag
    * but differentiate messages by embedding a type in each.
    */
   const int mpi_tag = 0;
   char edge_msg_type = 'e';
   char ack_msg_type = 'a';
   char upward_term_msg_type = 'u';
   char downward_term_msg_type = 'd';

   if (mpi1.hasReceivableMessage(0, MPI_ANY_SOURCE, mpi_tag)) {
      TBOX_ERROR("Connector::computeTransposeOf: not starting clean of receivable MPI messages.");
   }

   std::map<int, std::shared_ptr<tbox::MessageStream> > messages;
   std::vector<tbox::SAMRAI_MPI::Request> requests;
   tbox::SAMRAI_MPI::Status tmp_status;
   int mpi_err;

   // Send edge messages and remember to get receivers' acknowledgements.
   std::set<int> ack_needed;
   BoxContainer unshifted_head_nabrs;
   for (FullNeighborhoodSet::iterator rr = reordered_relationships.begin();
        rr != reordered_relationships.end(); ++rr) {

      const Box& base_box = rr->first;
      const BoxContainer& head_nabrs = rr->second;
      TBOX_ASSERT(!base_box.isPeriodicImage());

      /*
       * If base_box is local, store the neighbors.
       * Else, send neighbors to base_box's owner to store.
       */
      if (base_box.getOwnerRank() == mpi1.getRank()) {
         insertNeighbors(head_nabrs, base_box.getBoxId());
      } else {
         std::shared_ptr<tbox::MessageStream>& mstream = messages[base_box.getOwnerRank()];
         if (!mstream) {
            mstream.reset(new tbox::MessageStream);
            *mstream << edge_msg_type;
         }
         *mstream << base_box.getLocalId() << static_cast<size_t>(head_nabrs.size());
         for (BoxContainer::const_iterator bi = head_nabrs.begin(); bi != head_nabrs.end(); ++bi) {
            *mstream << *bi;
         }

         FullNeighborhoodSet::iterator nextrr = rr;
         ++nextrr;
         if (nextrr == reordered_relationships.end() ||
             nextrr->first.getOwnerRank() != base_box.getOwnerRank()) {
            requests.push_back(tbox::SAMRAI_MPI::Request());
            mpi_err = mpi1.Isend((void *)mstream->getBufferStart(),
                  static_cast<int>(mstream->getCurrentSize()), MPI_CHAR,
                  base_box.getOwnerRank(), mpi_tag,
                  &requests.back());
            TBOX_ASSERT(mpi_err == MPI_SUCCESS);
            ack_needed.insert(base_box.getOwnerRank());
         }
      }
   }

   // Data for propagating termination messages on the rank tree.
   tbox::CenteredRankTree rank_tree(mpi1);
   size_t child_term_needed = rank_tree.getNumberOfChildren();
   bool send_upward_term_msg = mpi1.getSize() > 1;

   if (ack_needed.empty() && child_term_needed == 0 && send_upward_term_msg) {
      // Leaves of the tree initiate upward termination message if no edge communication.
      requests.push_back(tbox::SAMRAI_MPI::Request());
      mpi_err = mpi1.Isend(&upward_term_msg_type, 1, MPI_CHAR,
            rank_tree.getParentRank(), mpi_tag,
            &requests.back());
      TBOX_ASSERT(mpi_err == MPI_SUCCESS);
      send_upward_term_msg = false;
   }

   /*
    * Receive edge messages and propgate termination messages: Both
    * communications must occur simultaneously.  Don't know where the
    * next message will come from, so must receive from any source.
    * Process messages based on the embedded msg_type.  Stop when
    * there are no edge messages are in transit, indicated by the
    * downward termination message.  Single process execution bypasses
    * communication by setting msg_type to downward termination.
    */
   int msg_length = 0;
   std::vector<char> recv_buffer;
   Box tmp_box(getBase().getDim());
   BoxContainer tmp_boxes;
   char msg_type = mpi1.getSize() == 1 ? downward_term_msg_type : char(0);

   while (msg_type != downward_term_msg_type) {

      mpi_err = mpi1.Probe(MPI_ANY_SOURCE, mpi_tag, &tmp_status);
      TBOX_ASSERT(mpi_err == MPI_SUCCESS);
      tbox::SAMRAI_MPI::Get_count(&tmp_status, MPI_CHAR, &msg_length);
      recv_buffer.resize(msg_length);
      mpi_err = mpi1.Recv(&recv_buffer[0], msg_length, MPI_CHAR,
            tmp_status.MPI_SOURCE, mpi_tag, &tmp_status);
      TBOX_ASSERT(mpi_err == MPI_SUCCESS);

      tbox::MessageStream mstream(recv_buffer.size(), tbox::MessageStream::Read,
                                  &recv_buffer[0], false);

      mstream >> msg_type;
      if (msg_type == edge_msg_type) {

         // Edge messages require acknowledgement and unpacking.
         requests.push_back(tbox::SAMRAI_MPI::Request());
         mpi_err = mpi1.Isend(static_cast<void *>(&ack_msg_type), 1, MPI_CHAR,
               tmp_status.MPI_SOURCE, mpi_tag, &requests.back());
         TBOX_ASSERT(mpi_err == MPI_SUCCESS);
         do {
            LocalId lid;
            size_t num_nabrs;
            mstream >> lid >> num_nabrs;
            for (size_t i = 0; i < num_nabrs; ++i) {
               mstream >> tmp_box;
               tmp_boxes.insert(tmp_box);
            }
            insertNeighbors(tmp_boxes, BoxId(lid, mpi1.getRank()));
            tmp_boxes.clear();
         } while (!mstream.endOfData());

      } else if (msg_type == ack_msg_type) {
         TBOX_ASSERT(ack_needed.find(tmp_status.MPI_SOURCE) != ack_needed.end());
         ack_needed.erase(tmp_status.MPI_SOURCE);
      } else if (msg_type == upward_term_msg_type) {
         TBOX_ASSERT(child_term_needed > 0);
         --child_term_needed;
      } else if (msg_type == downward_term_msg_type) {
         TBOX_ASSERT(child_term_needed == 0);
         // Propagate termination message downward.
         for (unsigned int ci = 0; ci < rank_tree.getNumberOfChildren(); ++ci) {
            requests.push_back(tbox::SAMRAI_MPI::Request());
            mpi_err = mpi1.Isend(&downward_term_msg_type, 1, MPI_CHAR,
                  rank_tree.getChildRank(ci), mpi_tag,
                  &requests.back());
            TBOX_ASSERT(mpi_err == MPI_SUCCESS);
         }
      } else {
         TBOX_ERROR(
            "Connector::computeTransposeOf: Library error: msg_type "
            << static_cast<int>(msg_type)
            <<
            " unrecognized,\npossibly due to receiving unrelated message.");
      }
      TBOX_ASSERT(mstream.endOfData());

      if (ack_needed.empty() && child_term_needed == 0 && send_upward_term_msg) {
         if (rank_tree.isRoot()) {
            // Initiate downward termination message.
            for (unsigned int ci = 0; ci < rank_tree.getNumberOfChildren(); ++ci) {
               requests.push_back(tbox::SAMRAI_MPI::Request());
               mpi_err = mpi1.Isend(&downward_term_msg_type, 1, MPI_CHAR,
                     rank_tree.getChildRank(ci), mpi_tag,
                     &requests.back());
               TBOX_ASSERT(mpi_err == MPI_SUCCESS);
            }
            msg_type = downward_term_msg_type;
         } else {
            // Propagate upward termination message.
            requests.push_back(tbox::SAMRAI_MPI::Request());
            mpi_err = mpi1.Isend(&upward_term_msg_type, 1, MPI_CHAR,
                  rank_tree.getParentRank(), mpi_tag,
                  &requests.back());
            TBOX_ASSERT(mpi_err == MPI_SUCCESS);
         }
         send_upward_term_msg = false;
      }

      recv_buffer.clear();
   }
   NULL_USE(mpi_err);

   if (!requests.empty()) {
      // Complete sends before allowing memory deallocation.
      std::vector<tbox::SAMRAI_MPI::Status> statuses(requests.size());
      tbox::SAMRAI_MPI::Waitall(static_cast<int>(requests.size()), &requests[0], &statuses[0]);
   }

   if (mpi1.hasReceivableMessage(0, MPI_ANY_SOURCE, mpi_tag)) {
      TBOX_ERROR("Connector::computeTransposeOf: not finishing clean of receivable MPI messages.");
   }
}

/*
 ***********************************************************************
 * This method does 2 important things with the edges:
 *
 * 1. It puts the edge data in head-major order so the base owners can
 * easily loop through the head--->base edges in the same order that
 * head owners see them.
 *
 * 2. It shifts periodic image head Boxes back to the zero-shift
 * position, and applies a similar shift to base Boxes so that the
 * overlap is unchanged.
 ***********************************************************************
 */
void
Connector::reorderRelationshipsByHead(
   std::map<Box, BoxContainer, Box::id_less>& relationships_by_head) const
{
   const tbox::Dimension& dim(getBase().getDim());

   const BoxLevel& base_box_level = getBase();
   const IntVector& base_ratio = getBase().getRefinementRatio();
   const IntVector& head_ratio = getHead().getRefinementRatio();

   const PeriodicShiftCatalog& shift_catalog =
      base_box_level.getGridGeometry()->getPeriodicShiftCatalog();

   Box shifted_box(dim), unshifted_nabr(dim);
   relationships_by_head.clear();
   for (Connector::ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {
      const Box& base_box = *base_box_level.getBoxStrict(*ci);
      for (Connector::ConstNeighborIterator na = begin(ci); na != end(ci); ++na) {
         const Box& nabr = *na;
         if (nabr.isPeriodicImage()) {
            shifted_box.initialize(
               base_box,
               shift_catalog.getOppositeShiftNumber(nabr.getPeriodicId()),
               base_ratio,
               shift_catalog);
            unshifted_nabr.initialize(
               nabr,
               shift_catalog.getZeroShiftNumber(),
               head_ratio,
               shift_catalog);
            relationships_by_head[unshifted_nabr].insert(shifted_box);
         } else {
            relationships_by_head[nabr].insert(base_box);
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods()
{
   tbox::SAMRAI_MPI mpi(getMPI());
   if (mpi.getSize() == 1) {
      // In single-proc mode, we already have all the relationships already.
      d_global_relationships = d_relationships;
      return;
   }

   t_acquire_remote_relationships->start();

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack relationships from all box_level relationship sets into a single message.
    * Note that each box_level relationship set object packs the size of its
    * sub-message into send_mesg.
    */
   acquireRemoteNeighborhoods_pack(send_mesg);
   int send_mesg_size = static_cast<int>(send_mesg.size());

   /*
    * Send and receive the data.
    */

   std::vector<int> recv_mesg_size(getMPI().getSize());
   mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(getMPI().getSize());
   int totl_size = 0;
   for (int n = 0; n < getMPI().getSize(); ++n) {
      proc_offset[n] = totl_size;
      totl_size += recv_mesg_size[n];
   }
   recv_mesg.resize(totl_size, tbox::MathUtilities<int>::getMax());
   mpi.Allgatherv(&send_mesg[0],
      send_mesg_size,
      MPI_INT,
      &recv_mesg[0],
      &recv_mesg_size[0],
      &proc_offset[0],
      MPI_INT);

   /*
    * Extract relationship info received from other processors.
    */
   acquireRemoteNeighborhoods_unpack(recv_mesg, proc_offset);

   t_acquire_remote_relationships->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods_pack(
   std::vector<int>& send_mesg) const
{
   const tbox::Dimension& dim = getBase().getDim();
   d_relationships.putToIntBuffer(send_mesg,
      dim,
      tbox::MathUtilities<int>::getMax());
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods_unpack(
   const std::vector<int>& recv_mesg,
   const std::vector<int>& proc_offset)
{
   const tbox::Dimension& dim = getBase().getDim();
   const int num_procs = getMPI().getSize();
   const int rank = getMPI().getRank();
   d_global_relationships = d_relationships;
   d_global_relationships.getFromIntBuffer(
      recv_mesg,
      proc_offset,
      dim,
      num_procs,
      rank);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::setParallelState(
   const BoxLevel::ParallelState parallel_state)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isFinalized()) {
      TBOX_ERROR(
         "Connector::setParallelState: Cannot change the parallel state of\n"
         << "an unfinalized Connector.  See Connector::finalizeContext()");
   }
#endif
   if (d_parallel_state == BoxLevel::DISTRIBUTED && parallel_state ==
       BoxLevel::GLOBALIZED) {
      acquireRemoteNeighborhoods();
   } else if (d_parallel_state == BoxLevel::GLOBALIZED && parallel_state ==
              BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::finalizeContext()
{
   TBOX_ASSERT(d_base_handle);
   TBOX_ASSERT(d_head_handle);

   const BoxLevel& base = d_base_handle->getBoxLevel();
   const BoxLevel& head = d_head_handle->getBoxLevel();
   const IntVector& baseRefinementRatio = base.getRefinementRatio();
   const IntVector& headRefinementRatio = head.getRefinementRatio();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!base.getMPI().isCongruentWith(head.getMPI())) {
      TBOX_ERROR("Connector::finalizeContext()\n"
         "base and head MPI communicators must be congruent.");
   }
#endif
   if (base.getGridGeometry() != head.getGridGeometry()) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Connector must be finalized with\n"
         << "BoxLevels using the same grid geometry.");
   }
   if (!(baseRefinementRatio >= headRefinementRatio ||
         baseRefinementRatio <= headRefinementRatio)) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Refinement ratio between base and head box_levels\n"
         << "cannot be mixed (bigger in some direction and\n"
         << "smaller in others).\n"
         << "Input base ratio = " << baseRefinementRatio
         << "\n"
         << "Input head ratio = " << headRefinementRatio
         << "\n");
   }
   if (d_parallel_state == BoxLevel::GLOBALIZED &&
       base.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR(
         "Connector::finalizeContext: base BoxLevel must be in\n"
         << "GLOBALIZED state for the Connector to be in\n"
         << "GLOBALIZED state.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   bool errf = false;
   for (NeighborhoodIterator ci = begin(); ci != end(); ++ci) {
      if (!base.hasBox(*ci)) {
         tbox::perr << "\nConnector::finalizeContext: NeighborhoodSet "
                    << "provided for non-existent box " << *ci
                    << "\n" << "Neighbors ("
                    << d_relationships.numNeighbors(*ci) << "):\n";
         for (NeighborIterator na = begin(ci); na != end(ci); ++na) {
            tbox::perr << (*na) << "\n";
         }
         errf = true;
      }
   }
   if (errf) {
      TBOX_ERROR(
         "Exiting due to errors."
         << "\nConnector::finalizeContext base box_level:\n"
         << base.format());
   }
#endif
   computeRatioInfo(
      baseRefinementRatio,
      headRefinementRatio,
      d_ratio,
      d_head_coarser,
      d_ratio_is_exact);

   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   } else {
      if (&d_relationships != &d_global_relationships) {
         d_global_relationships = d_relationships;
      }
   }

   // Erase remote relationships, if any, from d_relationships.
   // Note that we don't call getMPI here to get the rank as d_finalized isn't
   // set until after this step is complete.  It may be picky to insist that
   // d_finalized be set at the very end of the method but it's more correct.
   d_relationships.eraseNonLocalNeighborhoods(
      d_base_handle->getBoxLevel().getMPI().getRank());

   d_finalized = true;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setBase(
   const BoxLevel& new_base,
   bool finalize_context)
{
   if (!new_base.isInitialized()) {
      TBOX_ERROR("Connector::setBase():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_base_handle = new_base.getBoxLevelHandle();
   d_mpi = new_base.getMPI();
   if (finalize_context) {
      finalizeContext();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setHead(
   const BoxLevel& new_head,
   bool finalize_context)
{
   if (!new_head.isInitialized()) {
      TBOX_ERROR("Connector::setHead():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_head_handle = new_head.getBoxLevelHandle();
   if (finalize_context) {
      finalizeContext();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setWidth(
   const IntVector& new_width,
   bool finalize_context)
{
   if (!(new_width >= IntVector::getZero(new_width.getDim()))) {
      TBOX_ERROR("Connector::setWidth():\n"
         << "Invalid connector width: "
         << new_width << "\n");
   }

   d_finalized = false;
   d_base_width = new_width;
   size_t num_blocks =
      d_base_handle->getBoxLevel().getGridGeometry()->getNumberBlocks();

   if (d_base_width.getNumBlocks() == 1 && num_blocks != 1) {
      if (d_base_width.max() == d_base_width.min()) {
         d_base_width = IntVector(d_base_width, num_blocks);
      } else {
         TBOX_ERROR("Connector::setWidth: new anisotropic connector\n"
            << "width " << d_base_width << " must be \n"
            << "defined for " << num_blocks << " blocks.\n");
      }
   }

   if (finalize_context) {
      finalizeContext();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::computeRatioInfo(
   const IntVector& baseRefinementRatio,
   const IntVector& headRefinementRatio,
   IntVector& ratio,
   bool& head_coarser,
   bool& ratio_is_exact)
{
   if (baseRefinementRatio <= headRefinementRatio) {
      ratio = headRefinementRatio / baseRefinementRatio;
      head_coarser = false;
      ratio_is_exact = (ratio * baseRefinementRatio) == headRefinementRatio;
   } else {
      ratio = baseRefinementRatio / headRefinementRatio;
      head_coarser = true;
      ratio_is_exact = (ratio * headRefinementRatio) == baseRefinementRatio;
   }
   if (baseRefinementRatio * headRefinementRatio <
       IntVector::getZero(baseRefinementRatio.getDim())) {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -headRefinementRatio * baseRefinementRatio;
      ratio_is_exact = true;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::writeNeighborhoodsToErrorStream(
   const std::string& border) const
{
   const std::string indented_border = border + "  ";
   tbox::perr << border << "  " << d_relationships.numBoxNeighborhoods()
              << " neigborhoods:\n";
   for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
      tbox::perr << border << "  " << *ei << "\n";
      tbox::perr << border << "    Neighbors ("
                 << d_relationships.numNeighbors(ei) << "):\n";
      for (ConstNeighborIterator bi = begin(ei); bi != end(ei); ++bi) {
         const Box& box = *bi;
         tbox::perr << border << "    "
                    << box << "   "
                    << box.numberCells() << '\n';
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::writeNeighborhoodToStream(
   std::ostream& os,
   const BoxId& box_id) const
{
   const BoxNeighborhoodCollection& relationships = getRelations(box_id);
   BoxId non_per_id(box_id.getGlobalId(), PeriodicId::zero());
   ConstNeighborhoodIterator ei = relationships.find(non_per_id);
   if (ei == relationships.end()) {
      TBOX_ERROR("Connector::find: No neighbor set exists for\n"
         << "box " << box_id << ".\n");
   }
   for (ConstNeighborIterator bi = relationships.begin(ei);
        bi != relationships.end(ei); ++bi) {
      const Box& box = *bi;
      os << "    " << box << "   " << box.numberCells() << '\n';
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector *
Connector::createLocalTranspose() const
{
   const IntVector transpose_width = convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         getConnectorWidth());

   Connector* transpose = new Connector(getHead(), getBase(), transpose_width);
   doLocalTransposeWork(transpose);
   return transpose;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector *
Connector::createTranspose() const
{
   Connector* transpose =
      new Connector(getHead(),
         getBase(),
         convertHeadWidthToBase(getBase().getRefinementRatio(),
            getHead().getRefinementRatio(),
            getConnectorWidth()));

   doTransposeWork(transpose);
   return transpose;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::doLocalTransposeWork(
   Connector* transpose) const
{
   TBOX_ASSERT(transpose);
   TBOX_ASSERT(isTransposeOf(*transpose));

   const PeriodicShiftCatalog& shift_catalog =
      getBase().getGridGeometry()->getPeriodicShiftCatalog();

   for (ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {

      const BoxId& box_id = *ci;
      const BoxContainer::const_iterator ni = transpose->getHead().getBox(box_id);
      if (ni == transpose->getHead().getBoxes().end()) {
         TBOX_ERROR(
            "Connector::createLocalTranspose: box index\n"
            << box_id
            << " not found in local part of head box_level.\n"
            << "This means that the Connector data was not a\n"
            << "self-consistent local mapping.\n");
      }
      const Box& my_head_box = *ni;

      for (ConstNeighborIterator na = begin(ci); na != end(ci); ++na) {
         const Box& my_base_box = *na;
         if (my_base_box.getOwnerRank() != transpose->getMPI().getRank()) {
            TBOX_ERROR(
               "Connector::createLocalTranspose: base box "
               << my_head_box << "\n"
               << "has remote neighbor " << my_base_box
               << " which is disallowed.\n"
               << "Boxes must have only local neighbors in this method.");
         }
         if (my_base_box.isPeriodicImage()) {
            Box my_shifted_head_box(
               my_head_box,
               shift_catalog.getOppositeShiftNumber(
                  my_base_box.getPeriodicId()),
               transpose->getHead().getRefinementRatio(),
               shift_catalog);
            if (transpose->getHead().hasBox(my_shifted_head_box)) {
               BoxId base_non_per_id(
                  my_base_box.getGlobalId(),
                  PeriodicId::zero());
               transpose->d_relationships.insert(
                  base_non_per_id,
                  my_shifted_head_box);
            }
         } else {
            transpose->d_relationships.insert(
               my_base_box.getBoxId(),
               my_head_box);
         }
      }

   }

   if (0) {
      tbox::perr << "end of createLocalTranspose:\n"
                 << "base:\n" << transpose->getBase().format("BASE->", 3)
                 << "head:\n" << transpose->getHead().format("HEAD->", 3)
                 << "this:\n" << transpose->format("RRRR->", 3)
                 << "r:\n" << format("THIS->", 3)
                 << "Checking this transpose correctness:" << std::endl;
      transpose->assertTransposeCorrectness(*this, false);
      tbox::perr << "Checking r's transpose correctness:" << std::endl;
      assertTransposeCorrectness(*transpose, false);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::doTransposeWork(Connector* transpose) const
{
   TBOX_ASSERT(transpose);
   TBOX_ASSERT(isTransposeOf(*transpose));

   const tbox::Dimension dim(getBase().getDim());

   const Connector* globalized =
      (d_parallel_state == BoxLevel::GLOBALIZED) ?
      this : makeGlobalizedCopy(*this);

   const BoxLevel& globalized_base = getBase().getGlobalizedVersion();
   const BoxContainer& globalized_boxes = globalized_base.getGlobalBoxes();

   for (BoxNeighborhoodCollection::ConstIterator ni = globalized->d_global_relationships.begin();
        ni != globalized->d_global_relationships.end(); ++ni) {

      for (Connector::ConstNeighborIterator na = begin(ni); na != end(ni); ++na) {
         if (na->getOwnerRank() == globalized_base.getMPI().getRank()) {
            if (!na->isPeriodicImage()) {
               TBOX_ASSERT(getHead().hasBox(*na));
               transpose->insertLocalNeighbor(
                  *globalized_boxes.find(Box(dim, *ni)),
                  na->getBoxId());
            } else {
               // Need to do shifting.
               TBOX_ERROR("Unfinished Code!!!");
            }
         }
      }

   }

   if (globalized != this) {
      delete globalized;
      globalized = 0;
   }

   if (0) {
      TBOX_ASSERT(checkTransposeCorrectness(*transpose));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

bool
Connector::isTransposeOf(
   const Connector& other) const
{
   bool rval = false;
   if (d_base_handle == other.d_head_handle &&
       d_head_handle == other.d_base_handle) {
      if (d_head_coarser) {
         IntVector transpose_base_width = convertHeadWidthToBase(
               getHead().getRefinementRatio(),
               getBase().getRefinementRatio(),
               d_base_width);
         rval = other.d_base_width == transpose_base_width;
      } else {
         IntVector transpose_base_width = convertHeadWidthToBase(
               other.getHead().getRefinementRatio(),
               other.getBase().getRefinementRatio(),
               other.d_base_width);
         rval = d_base_width == transpose_base_width;
      }
   }
   return rval;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isFinalized());

   if (d_global_data_up_to_date) {
      return;
   }

   t_cache_global_reduced_data->barrierAndStart();

   tbox::SAMRAI_MPI mpi(getMPI());

   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_number_of_relationships =
         d_global_relationships.sumNumNeighbors();
      d_global_number_of_neighbor_sets =
         d_global_relationships.numBoxNeighborhoods();
   } else {
      if (mpi.getSize() > 1) {
         int tmpa[2], tmpb[2];
         tmpa[0] = getLocalNumberOfNeighborSets();
         tmpa[1] = getLocalNumberOfRelationships();

         TBOX_ASSERT(tmpa[0] >= 0);
         TBOX_ASSERT(tmpa[0] >= 0);

         mpi.Allreduce(tmpa,
            tmpb,                        // Better to use MPI_IN_PLACE, but not some MPI's do not support.
            2,
            MPI_INT,
            MPI_SUM);
         d_global_number_of_neighbor_sets = tmpb[0];
         d_global_number_of_relationships = tmpb[1];
      } else {
         d_global_number_of_neighbor_sets = getLocalNumberOfNeighborSets();
         d_global_number_of_relationships = getLocalNumberOfRelationships();
      }

      TBOX_ASSERT(d_global_number_of_neighbor_sets >= 0);
      TBOX_ASSERT(d_global_number_of_relationships >= 0);
   }

   d_global_data_up_to_date = true;

   t_cache_global_reduced_data->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

IntVector
Connector::convertHeadWidthToBase(
   const IntVector& base_refinement_ratio,
   const IntVector& head_refinement_ratio,
   const IntVector& head_width)
{
   if (!(base_refinement_ratio >= head_refinement_ratio ||
         base_refinement_ratio <= head_refinement_ratio)) {
      TBOX_ERROR("Connector::convertHeadWidthToBase:\n"
         << "head box_level must be either\n"
         << "finer or coarser than base.\n"
         << "Combined refinement and coarsening not allowed.");
   }

   tbox::Dimension dim(head_refinement_ratio.getDim());
   const size_t num_blocks = base_refinement_ratio.getNumBlocks();

   IntVector tmp_head_width(head_width);
   if (tmp_head_width.getNumBlocks() == 1 && num_blocks != 1) { 
      if (tmp_head_width.max() == tmp_head_width.min()) {
         tmp_head_width = IntVector(tmp_head_width, num_blocks);
      } else {
         TBOX_ERROR("Connector::convertHeadWidthToBase: anisotropic connector\n"
            << "width " << head_width << " must be \n"
            << "defined for " << num_blocks << " blocks.\n");
      }
   }

   IntVector ratio(dim); // Ratio between head and base.

   if ((head_refinement_ratio > IntVector::getZero(dim)) ==
       (base_refinement_ratio > IntVector::getZero(dim))) {
      // Same signs for both ratios -> simple to compute head-base ratio.
      if (base_refinement_ratio >= head_refinement_ratio) {
         ratio = base_refinement_ratio / head_refinement_ratio;
      } else {
         ratio = head_refinement_ratio / base_refinement_ratio;
      }
   } else {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -base_refinement_ratio * head_refinement_ratio;
   }
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));

   const IntVector base_width =
      (base_refinement_ratio >= head_refinement_ratio) ?
      (tmp_head_width * ratio) : IntVector::ceilingDivide(tmp_head_width, ratio);

   return base_width;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int detail_depth) const
{
   if (detail_depth < 0) {
      return;
   }

   if (!isFinalized()) {
      os << border << "Unfinalized.\n";
      return;
   }
   bool head_coarser = d_head_coarser;
   const IntVector head_width =
      convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         d_base_width);
   os << border << "Parallel state     : "
      << (getParallelState() == BoxLevel::DISTRIBUTED ? "DIST" : "GLOB")
      << '\n'
      << border << "Rank,nproc         : " << getMPI().getRank() << ", " << getMPI().getSize()
      << '\n'
      << border << "Base,head objects  :"
      << " ("
      << (d_base_handle == d_head_handle ? "same" : "different") << ") "
      << (void *)&d_base_handle->getBoxLevel() << ", "
      << (void *)&d_head_handle->getBoxLevel() << "\n"
      << border << "Base,head,/ ratios : "
      << getBase().getRefinementRatio() << ", "
      << getHead().getRefinementRatio() << ", "
      << d_ratio << (d_head_coarser ? " (head coarser)" : "") << '\n'
      << border << "Base,head widths   : " << d_base_width << ", "
      << head_width << '\n'
      << border << "Box count    : " << getBase().getLocalNumberOfBoxes()
      << " (" << getLocalNumberOfNeighborSets() << " with neighbor lists)\n"
   ;
   if (detail_depth > 0) {
      os << border << "Boxes with neighbors:\n";
      for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& box_id = *ei;
         BoxContainer::const_iterator ni = getBase().getBox(box_id);
         if (ni != getBase().getBoxes().end()) {
            os << border << "  "
               << (*ni) << "_"
               << (*ni).numberCells();
         } else {
            os << border << "  #"
               << box_id
               << ": INVALID DATA WARNING: No base box with this index!";
         }
         os << "  has " << numLocalNeighbors(box_id) << " neighbors:"
            << ((detail_depth > 1) ? "\n" : " ...\n");
         if (detail_depth > 1) {
            for (ConstNeighborIterator i_nabr = begin(ei);
                 i_nabr != end(ei); ++i_nabr) {
               os << border << "      "
                  << (*i_nabr) << "_"
                  << i_nabr->numberCells();
               if (ni != getBase().getBoxes().end()) {
                  Box ovlap = *i_nabr;
                  if (ni->getBlockId() != i_nabr->getBlockId()) {
                     d_base_handle->getBoxLevel().getGridGeometry()->
                     transformBox(
                        ovlap,
                        d_head_handle->getBoxLevel().getRefinementRatio(),
                        ni->getBlockId(),
                        i_nabr->getBlockId());
                  }
                  if (ovlap.getBlockId() != ni->getBlockId()) {
                     os << "\tov undefined (non-touching blocks.";
                  } else {
                     if (head_coarser) {
                        ovlap.refine(d_ratio);
                     } else if (d_ratio != 1) {
                        ovlap.coarsen(d_ratio);
                     }
                     Box ghost_box = (*ni);
                     ghost_box.grow(d_base_width);
                     ovlap = ovlap * ghost_box;
                     os << "\tov" << ovlap << "_" << ovlap.numberCells();
                  }
               }
               os << '\n';
            }
         }
      }
   }
}

/*
 ***********************************************************************
 * Outputter copy constructor
 ***********************************************************************
 */

Connector::Outputter::Outputter(
   const Connector::Outputter& other): 
   d_conn(other.d_conn),
   d_border(other.d_border),
   d_detail_depth(other.d_detail_depth),
   d_output_statistics(other.d_output_statistics)
{
}

/*
 ***********************************************************************
 * Construct a Connector Outputter with formatting parameters.
 ***********************************************************************
 */

Connector::Outputter::Outputter(
   const Connector& connector,
   const std::string& border,
   int detail_depth,
   bool output_statistics):
   d_conn(connector),
   d_border(border),
   d_detail_depth(detail_depth),
   d_output_statistics(output_statistics)
{
}

/*
 ***********************************************************************
 * Print out a Connector according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream&
operator << (
   std::ostream& os,
   const Connector::Outputter& format)
{
   if (format.d_output_statistics) {
      ConnectorStatistics cs(format.d_conn);
      cs.printNeighborStats(os, format.d_border);
   } else {
      format.d_conn.recursivePrint(os, format.d_border, format.d_detail_depth);
   }
   return os;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector *
Connector::makeGlobalizedCopy(
   const Connector& other) const
{
   // Prevent wasteful accidental use when this method is not needed.
   TBOX_ASSERT(other.getParallelState() != BoxLevel::GLOBALIZED);

   Connector* copy = new Connector(other);
   copy->setParallelState(BoxLevel::GLOBALIZED);
   return copy;
}

/*
 ***********************************************************************
 * Run checkTransposeCorrectness and assert that no errors are found.
 ***********************************************************************
 */

void
Connector::assertTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   size_t err_count =
      checkTransposeCorrectness(input_transpose, ignore_periodic_relationships);
   if (err_count) {
      TBOX_ERROR(
         "Connector::assertTransposeCorrectness:\n"
         << "Aborting with " << err_count << " transpose errors found:\n"
         << "this base:\n" << getBase().format("B:", 3)
         << "this head:\n" << getHead().format("H:", 3)
         << "this Connector:\n" << format("B->H:", 3)
         << "\ntranspose Connector:\n" << input_transpose.format("H->B:", 3));
   }
}

/*
 ***********************************************************************
 *
 * For every relationship in this, there should be reverse relationship in
 * transpose.
 *
 * This method does not check whether the Connectors are defined to
 * form logical transposes (based on their widths and their base and
 * head box_levels).  For that, see isTransposeOf().
 *
 ***********************************************************************
 */

size_t
Connector::checkTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   const tbox::Dimension dim(getBase().getDim());

   const Connector* transpose =
      (input_transpose.d_parallel_state == BoxLevel::GLOBALIZED) ?
      &input_transpose : makeGlobalizedCopy(input_transpose);

   const BoxLevel& head = getHead().getGlobalizedVersion();

   const PeriodicShiftCatalog& shift_catalog =
      head.getGridGeometry()->getPeriodicShiftCatalog();

   /*
    * Check for extraneous relationships.
    * For every relationship in this, there should be reverse relationship in transpose.
    */
   Box shifted_box(dim);   // Shifted version of an unshifted Box.
   Box unshifted_box(dim); // Unhifted version of a shifted Box.

   size_t err_count = 0;

   const BoxNeighborhoodCollection& tran_relationships =
      transpose->getGlobalNeighborhoodSets();
   for (ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {

      const BoxId& box_id = *ci;
      const Box& box = *getBase().getBox(box_id);

      size_t err_count_for_current_index = 0;

      for (ConstNeighborIterator ni = begin(ci); ni != end(ci); ++ni) {

         if (ignore_periodic_relationships && ni->isPeriodicImage()) {
            continue;
         }

         const Box& nabr = *ni;

         /*
          * Key for find in NeighborhoodSet must be non-periodic.
          */
         BoxId non_per_nabr_id(nabr.getGlobalId(),
                               PeriodicId::zero());

         ConstNeighborhoodIterator cn =
            tran_relationships.find(non_per_nabr_id);

         if (cn == tran_relationships.end()) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local box " << box
            << " has relationship to " << nabr
            << " but " << nabr << " has no relationship container.\n";
            ++err_count_for_current_index;
            continue;
         }

         TBOX_ASSERT(*cn == non_per_nabr_id);

         bool nabr_has_box;
         if (nabr.isPeriodicImage()) {
            shifted_box.initialize(
               box,
               shift_catalog.getOppositeShiftNumber(nabr.getPeriodicId()),
               getBase().getRefinementRatio(),
               shift_catalog);
            nabr_has_box =
               tran_relationships.hasNeighbor(cn, shifted_box);
         } else {
            nabr_has_box = tran_relationships.hasNeighbor(cn, box);
         }

         if (!nabr_has_box) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local box " << box;
            if (nabr.isPeriodicImage()) {
               tbox::perr << " (shifted version " << shifted_box << ")";
            }
            tbox::perr << " has relationship to " << nabr << " but "
            << nabr << " does not have the reverse relationship.\n"
            ;
            tbox::perr << "Neighbors of " << nabr << " are:\n";
            for (ConstNeighborIterator nj = tran_relationships.begin(cn);
                 nj != tran_relationships.end(cn); ++nj) {
               tbox::perr << "   " << *nj << std::endl;
            }
            ++err_count_for_current_index;
            continue;
         }

      }

      if (err_count_for_current_index > 0) {
         tbox::perr << "Box " << box << " had "
         << err_count_for_current_index
         << " errors.  Neighbors are:\n";
         for (ConstNeighborIterator nj = begin(ci); nj != end(ci); ++nj) {
            tbox::perr << "  " << *nj << std::endl;
         }
         err_count += err_count_for_current_index;
      }

   }

   /*
    * Check for missing relationships:
    * Transpose should not contain any relationship that does not correspond to
    * one in this.
    */

   for (ConstNeighborhoodIterator ci = tran_relationships.begin();
        ci != tran_relationships.end(); ++ci) {

      const BoxId& box_id = *ci;

      size_t err_count_for_current_index = 0;

      if (!head.hasBox(box_id)) {
         TBOX_ASSERT(head.hasBox(box_id));
      }
      const Box& head_box = *head.getBoxStrict(box_id);

      for (ConstNeighborIterator na = tran_relationships.begin(ci);
           na != tran_relationships.end(ci); ++na) {

         const Box nabr = *na;

         if (nabr.getOwnerRank() == getMPI().getRank()) {

            if (ignore_periodic_relationships && nabr.isPeriodicImage()) {
               continue;
            }

            if (!getBase().hasBox(nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head box " << head_box
               << " has neighbor " << nabr << "\n"
               << " but the neighbor does not exist "
               << "in the base box_level.\n";
               tbox::perr << "Neighbors of head box "
               << box_id << " are:\n";
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box& base_box = *getBase().getBoxStrict(nabr);

            /*
             * Non-periodic BoxId needed for NeighborhoodSet::find()
             */
            BoxId base_non_per_id(base_box.getGlobalId(), PeriodicId::zero());

            if (!d_relationships.isBaseBox(base_non_per_id)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head box " << head_box << "\n"
               << " has base box "
               << base_box << " as a neighbor.\n"
               << "But " << base_box
               << " has no neighbor container.\n";
               tbox::perr << "Neighbors of head box " << BoxId(box_id)
               << ":" << std::endl;
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box nabr_nabr(dim, box_id.getGlobalId(),
                                shift_catalog.getOppositeShiftNumber(
                                   base_box.getPeriodicId()));

            if (!d_relationships.hasNeighbor(base_non_per_id, nabr_nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head box " << head_box << "\n"
               << " has base box " << base_box
               << " as a neighbor.\n"
               << "But base box " << base_box
               << " does not have a box indexed "
               << nabr_nabr.getBoxId()
               << " in its neighbor list." << std::endl;
               tbox::perr << "Neighbors of head box " << nabr_nabr.getBoxId()
               << ":" << std::endl;
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               tbox::perr << "Neighbors of base box ";
               if (nabr.isPeriodicImage()) {
                  unshifted_box.initialize(
                     nabr,
                     shift_catalog.getZeroShiftNumber(),
                     getBase().getRefinementRatio(),
                     shift_catalog);
                  tbox::perr << unshifted_box;
               }
               tbox::perr << ":" << std::endl;
               ConstNeighborhoodIterator nabr_nabrs_ =
                  d_relationships.find(base_non_per_id);
               for (ConstNeighborIterator nj = begin(nabr_nabrs_);
                    nj != end(nabr_nabrs_); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

         }

      }

      if (err_count_for_current_index > 0) {
         err_count += err_count_for_current_index;
      }

   }

   if (transpose != &input_transpose) {
      delete transpose;
   }

   int global_err_count = static_cast<int>(err_count);
   if (getMPI().getSize() > 1) {
      getMPI().AllReduce(&global_err_count, 1, MPI_SUM);
   }

   return static_cast<size_t>(global_err_count);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t
Connector::checkConsistencyWithBase() const
{
   size_t num_errors = 0;
   for (ConstNeighborhoodIterator i_relationships = begin();
        i_relationships != end(); ++i_relationships) {
      const BoxId& box_id = *i_relationships;
      if (!getBase().hasBox(box_id)) {
         ++num_errors;
         tbox::plog << "ERROR->"
         << "Connector::checkConsistencyWithBase: Neighbor data given "
         << "\nfor box " << box_id
         << " but the box does not exist.\n";
      }
   }
   return num_errors;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::assertConsistencyWithBase() const
{
   if (checkConsistencyWithBase() > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithHead() found inconsistencies.\n"
         << "Base BoxLevel:\n" << getBase().format("base-> ", 3)
         << "Head BoxLevel:\n" << getHead().format("head-> ", 3)
         << "Connector:\n" << format("E-> ", 3));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::computeNeighborhoodDifferences(
   std::shared_ptr<Connector>& left_minus_right,
   const Connector& left,
   const Connector& right)
{
   if (0) {
      tbox::plog << "Computing relationship differences, a:\n" << left.format("A-> ")
      << "Computing relationship differences, b:\n" << right.format("B-> ");
   }
   left_minus_right.reset(new Connector(left.d_base_handle->getBoxLevel(),
         left.d_head_handle->getBoxLevel(),
         left.d_base_width,
         left.getParallelState()));

   for (ConstNeighborhoodIterator ai = left.begin(); ai != left.end(); ++ai) {

      const BoxId& box_id = *ai;

      ConstNeighborhoodIterator bi = right.findLocal(box_id);
      if (bi != right.end()) {
         // Remove bi from ai.  Put results in a_minus_b.

         /*
          * In theory, we should not have to resort to using std::set
          * in order to use set_difference, but our BoxContainer does
          * not implement all features necessary to use
          * set_difference.
          */
         std::set<Box, Box::id_less> anabrs(left.begin(ai), left.end(ai));
         std::set<Box, Box::id_less> bnabrs(right.begin(bi), right.end(bi));
         std::set<Box, Box::id_less> diff;
         std::insert_iterator<std::set<Box, Box::id_less> > ii(diff, diff.begin());
         std::set_difference(anabrs.begin(),
            anabrs.end(),
            bnabrs.begin(),
            bnabrs.end(),
            ii, Box::id_less());
         if (!diff.empty()) {
            NeighborhoodIterator base_box_itr =
               left_minus_right->makeEmptyLocalNeighborhood(box_id);
            for (std::set<Box, Box::id_less>::const_iterator ii = diff.begin();
                 ii != diff.end();
                 ++ii) {
               left_minus_right->insertLocalNeighbor(*ii, base_box_itr);
            }
         }
      } else if (left.numLocalNeighbors(box_id) != 0) {
         NeighborhoodIterator base_box_itr =
            left_minus_right->makeEmptyLocalNeighborhood(box_id);
         for (ConstNeighborIterator na = left.begin(ai);
              na != left.end(ai); ++na) {
            left_minus_right->insertLocalNeighbor(*na, base_box_itr);
         }
      }

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::assertConsistencyWithHead() const
{
   const size_t number_of_inconsistencies = checkConsistencyWithHead();
   if (number_of_inconsistencies > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithHead() found inconsistencies.\n"
         << "Base BoxLevel:\n" << getBase().format("base-> ", 3)
         << "Head BoxLevel:\n" << getHead().format("head-> ", 3)
         << "Connector:\n" << format("E-> ", 3));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t
Connector::checkConsistencyWithHead() const
{
   const BoxLevel& head_box_level = getHead().getGlobalizedVersion();

   TBOX_ASSERT(head_box_level.getParallelState() == BoxLevel::GLOBALIZED);

   const PeriodicShiftCatalog& shift_catalog =
      head_box_level.getGridGeometry()->getPeriodicShiftCatalog();

   const BoxContainer& head_boxes = head_box_level.getGlobalBoxes();

   size_t number_of_inconsistencies = 0;

   /*
    * For each neighbor in each neighbor list,
    * check that the neighbor is in the head_box_level.
    */

   for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {

      const BoxId& box_id = *ei;

      for (ConstNeighborIterator na = begin(ei); na != end(ei); ++na) {

         const Box& nabr = *na;
         const Box unshifted_nabr(
            nabr,
            PeriodicId::zero(),
            head_box_level.getRefinementRatio(),
            shift_catalog);

         BoxContainer::const_iterator na_in_head =
            head_boxes.find(unshifted_nabr);

         if (na_in_head == head_boxes.end()) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Neighbor list for box " << box_id << "\n"
            << "referenced nonexistent neighbor "
            << nabr << "\n";
            tbox::perr << "Neighbors of box " << box_id << ":\n";
            for (ConstNeighborIterator nb = begin(ei); nb != end(ei); ++nb) {
               tbox::perr << "   " << *nb << '\n';
            }
            ++number_of_inconsistencies;
            continue;
         }

         const Box& nabr_in_head = *na_in_head;
         if (!unshifted_nabr.isIdEqual(nabr_in_head) ||
             !unshifted_nabr.isSpatiallyEqual(nabr_in_head)) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Inconsistent box data at box "
            << box_id << "\n"
            << "Neighbor " << nabr << "(unshifted to "
            << unshifted_nabr << ") does not match "
            << "head box " << nabr_in_head
            << "\n";
            ++number_of_inconsistencies;
         }

      }
   }

   return number_of_inconsistencies;
}

/*
 ***********************************************************************
 * Checking is done as follows:
 *   - Rebuild the overlap containers using findOverlaps().
 *     Note that the rebuilt overlap set is complete.
 *   - Check the current overlap set against the rebuilt overlap set
 *     to find missing overlaps and extra overlaps.
 *
 * Currently, the rebuilt overlaps are rebuilt using findOverlaps().
 * Thus, it may be pointless to use this method as a check for that
 * method.
 ***********************************************************************
 */

void
Connector::findOverlapErrors(
   std::shared_ptr<Connector>& missing,
   std::shared_ptr<Connector>& extra,
   bool ignore_self_overlap) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!getBase().isInitialized() || !getHead().isInitialized()) {
      TBOX_ERROR(
         "Connector::findOverlapErrors: Cannot check overlaps\n"
         << "when base or head box_level is uninitialized.");
   }
#endif

   /*
    * Obtain a globalized version of the head for checking.
    */
   const BoxLevel& head = getHead().getGlobalizedVersion();

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Before checking on overlap errors, make sure the user gave a
    * valid Connector.
    *
    * Each neighbor set should correspond to a base box.
    *
    * Each referenced neighbor should exist in the head.
    */
   size_t num_base_consistency_errors = checkConsistencyWithBase();
   size_t num_head_consistency_errors = checkConsistencyWithHead();
   if (num_base_consistency_errors > 0) {
      tbox::perr
      << "Connector::findOverlapErrors: cannot check overlap errors\n"
      << "for inconsistent base data.\n";
   }
   if (num_head_consistency_errors > 0) {
      tbox::perr
      << "Connector::findOverlapErrors: cannot check overlap errors\n"
      << "for inconsistent head data.\n";
   }
   if (num_base_consistency_errors || num_head_consistency_errors) {
      TBOX_ERROR(
         "Connector::findOverlapErrors exiting due to\n"
         << "inconsistent data.\n"
         << "Base:\n" << getBase().format("B->", 2)
         << "Head:\n" << getHead().format("H->", 2)
         << "Connector:\n" << format("C->", 3));
   }
#endif

   /*
    * Rebuild the overlap Connector for checking.
    */
   Connector rebuilt(getBase(), getHead(), getConnectorWidth());
   rebuilt.findOverlaps_rbbt(head, ignore_self_overlap);

   /*
    * Check that the rebuilt overlaps match the existing overlaps.
    *
    * Currently, we use findOverlaps to rebuild the overlaps.
    * Thus, it may be pointless to use this method
    * as a check for that method.
    */
   computeNeighborhoodDifferences(extra, *this, rebuilt);
   TBOX_ASSERT(&extra->getBase() == &getBase());
   TBOX_ASSERT(&extra->getHead() == &getHead());
   computeNeighborhoodDifferences(missing, rebuilt, *this);
   TBOX_ASSERT(&missing->getBase() == &getBase());
   TBOX_ASSERT(&missing->getHead() == &getHead());
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::assertOverlapCorrectness(
   bool ignore_self_overlap,
   bool assert_completeness,
   bool ignore_periodic_images) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!getBase().isInitialized() || !getHead().isInitialized()) {
      TBOX_ERROR(
         "Connector::assertOverlapCorrectness: Cannot check overlaps\n"
         << "when base or head box_level is uninitialized.");

   }
#endif

   int local_error_count = checkOverlapCorrectness(ignore_self_overlap,
         assert_completeness,
         ignore_periodic_images);

   const tbox::SAMRAI_MPI& mpi(getMPI());
   int max_error_count = local_error_count;
   int rank_of_max = mpi.getRank();
   if (mpi.getSize() > 1) {
      IntIntStruct send, recv;
      send.rank = recv.rank = mpi.getRank();
      send.i = local_error_count;
      mpi.Allreduce(&send, &recv, 1, MPI_2INT, MPI_MAXLOC);
      max_error_count = recv.i;
      rank_of_max = recv.rank;
   }
   if (max_error_count > 0) {
      TBOX_ERROR(
         "Connector::assertOverlapCorrectness found missing and/or extra overlaps.\n"
         << "Error in connector, " << local_error_count
         << " local errors, "
         << max_error_count << " max errors on proc " << rank_of_max
         << ":\n"
         << format("E-> ")
         << "base:\n" << getBase().format("B-> ")
         << "head:\n" << getHead().format("H-> "));
   }
}

/*
 ***********************************************************************
 * Return number of missing and number of extra overlaps.
 ***********************************************************************
 */

int
Connector::checkOverlapCorrectness(
   bool ignore_self_overlap,
   bool assert_completeness,
   bool ignore_periodic_images) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!getBase().isInitialized() || !getHead().isInitialized()) {
      TBOX_ERROR(
         "Connector::checkOverlapCorrectness: Cannot check overlaps when\n"
         << "base or head box_level is uninitialized.");

   }
#endif
   TBOX_ASSERT(!hasPeriodicLocalNeighborhoodBaseBoxes());

   std::shared_ptr<Connector> missing, extra;
   findOverlapErrors(missing, extra, ignore_self_overlap);

   if (!assert_completeness) {
      // Disregard missing overlaps by resetting missing to empty.
      missing->clearNeighborhoods();
   } else if (ignore_periodic_images) {
      // Disregard missing overlaps if they are incident on a periodic box.
      missing->removePeriodicLocalNeighbors();
      missing->eraseEmptyNeighborSets();
   }

   const BoxId dummy_box_id;

   /*
    * Report the errors found, ordered by the Box where the
    * error appears.  In order to do this, we have to loop through
    * the neighborhoods of missing and extra at the same time.
    */

   Connector::ConstNeighborhoodIterator im = missing->begin();
   Connector::ConstNeighborhoodIterator ie = extra->begin();
   for ( ; im != missing->end() || ie != extra->end();
         /* incremented in loop */) {

      const BoxId& global_id_missing =
         im == missing->end() ? dummy_box_id : *im;
      const BoxId& global_id_extra =
         ie == extra->end() ? dummy_box_id : *ie;

      if (im != missing->end() && ie != extra->end() &&
          *im == *ie) {

         /*
          * im and ie are pointing at the same Box.  Report the
          * error for this Box.
          */

         const Box& box = *getBase().getBoxStrict(global_id_missing);
         tbox::perr << "Found " << missing->numLocalNeighbors(*im)
         << " missing and "
         << extra->numLocalNeighbors(*ie)
         << " extra overlaps for "
         << box << std::endl;
         Connector::ConstNeighborhoodIterator it = findLocal(global_id_missing);
         if (it == end()) {
            tbox::perr << "  Current Neighbors (no neighbor set)." << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
            << numLocalNeighbors(*it) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = begin(it);
                 na != end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         {
            tbox::perr << "  Missing Neighbors ("
            << missing->numLocalNeighbors(*im) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = missing->begin(im);
                 na != missing->end(im); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         {
            tbox::perr << "  Extra Neighbors ("
            << extra->numLocalNeighbors(*ie) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = extra->begin(ie);
                 na != extra->end(ie); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         ++im;
         ++ie;

      } else if ((ie == extra->end()) ||
                 (im != missing->end() && *im < *ie)) {

         /*
          * im goes before ie (or ie has reached the end).  Report the
          * errors for the Box at im.
          */

         const Box& box = *getBase().getBoxStrict(global_id_missing);
         tbox::perr << "Found " << missing->numLocalNeighbors(*im)
         << " missing overlaps for " << box << std::endl;
         Connector::ConstNeighborhoodIterator it = findLocal(global_id_missing);
         if (it == end()) {
            tbox::perr << "    Current Neighbors (no neighbor set)."
            << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
            << numLocalNeighbors(*it) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = begin(it);
                 na != end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         {
            tbox::perr << "  Missing Neighbors ("
            << missing->numLocalNeighbors(*im) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = missing->begin(im);
                 na != missing->end(im); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         ++im;
      } else if ((im == missing->end()) ||
                 (ie != extra->end() && *ie < *im)) {

         /*
          * ie goes before im (or im has reached the end).  Report the
          * errors for the Box at ie.
          */

         const Box& box = *getBase().getBoxStrict(
               global_id_extra);
         tbox::perr << "Found " << extra->numLocalNeighbors(*ie)
         << " extra overlaps for " << box << std::endl;
         Connector::ConstNeighborhoodIterator it = findLocal(global_id_extra);
         if (it == end()) {
            tbox::perr << "  Current Neighbors (no neighbor set)." << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
            << numLocalNeighbors(*it) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = begin(it);
                 na != end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         {
            tbox::perr << "  Extra Neighbors ("
            << extra->numLocalNeighbors(*ie) << "):"
            << std::endl;
            Box ghost_box = box;
            ghost_box.grow(getConnectorWidth());
            for (Connector::ConstNeighborIterator na = extra->begin(ie);
                 na != extra->end(ie); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (getHeadCoarserFlag()) {
                  nabr_box.refine(getRatio());
               } else if (getRatio() != 1) {
                  nabr_box.coarsen(getRatio());
               }
               if (nabr_box.getBlockId() != box.getBlockId()) {
                  getBase().getGridGeometry()->transformBox(nabr_box,
                     getBase().getRefinementRatio(),
                     box.getBlockId(),
                     nabr.getBlockId());
               }
               Box ovlap = nabr_box * ghost_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
               << "\tov" << ovlap << '_' << ovlap.numberCells()
               << std::endl;
            }
         }
         ++ie;
      }

   }

   return missing->getLocalNumberOfNeighborSets()
          + extra->getLocalNumberOfNeighborSets();
}

/*
 ***********************************************************************
 * ignore_self_overlap should be set to true only if
 * - the base and head box_levels represent the same box_level.
 *   Two different box_level objects may represent the same
 *   box_level if they are of the same refinement ratio.
 * - you want to ignore overlaps between a box and itself.
 ***********************************************************************
 */

void
Connector::findOverlaps_rbbt(
   const BoxLevel& head,
   bool ignore_self_overlap,
   bool sanity_check_method_postconditions)
{
   const tbox::Dimension dim(head.getDim());

   t_find_overlaps_rbbt->start();

   /*
    * Finding overlaps for this object, using
    * an externally provided head BoxLevel
    * meant to represent d_head.  We allow the
    * substitution of an external head because
    * we require the head is GLOBALIZED.  The
    * user may have a GLOBALIZED version already,
    * in which case we want to avoid the expense
    * of creating a temporary GLOBALIZED version.
    *
    * Global boxes provided by head are sorted in a BoxContainer
    * so they can be quickly searched to see which intersects the
    * boxes in this object.
    */
   if (head.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR("Connector::findOverlaps_rbbt() requires given head\n"
         << "to be GLOBALIZED.\n");
   }

   /*
    * The nomenclature "base" refers to the *this box_level
    * and "head" refer to the box_level in the argument.
    */
   const BoxLevel& base(getBase());

   /*
    * Determine relationship between base and head index spaces.
    */
   const bool head_is_finer =
      head.getRefinementRatio() >= base.getRefinementRatio() &&
      head.getRefinementRatio() != base.getRefinementRatio();
   const bool base_is_finer =
      base.getRefinementRatio() >= head.getRefinementRatio() &&
      base.getRefinementRatio() != head.getRefinementRatio();

   /*
    * Create single container of visible head boxes
    * to generate the search tree.
    */
   const BoxContainer& rbbt = head.getGlobalBoxes();
   rbbt.makeTree(head.getGridGeometry().get());

   /*
    * A neighbor of a Box would be discarded if
    * - ignore_self_overlap is true,
    * - the two are equal by comparison, and
    * - they are from box_levels with the same refinement ratio
    *   (we cannot compare box box_level pointers because that
    *   does not work when a box box_level is a temporary globalized object)
    */
   const bool discard_self_overlap =
      ignore_self_overlap &&
      (base.getRefinementRatio() == head.getRefinementRatio());

   /*
    * Discard current overlaps (if any).
    */
   clearNeighborhoods();

   /*
    * Use BoxTree to find local base Boxes intersecting head Boxes.
    */
   NeighborSet nabrs_for_box;
   const BoxContainer& base_boxes = base.getBoxes();
   for (RealBoxConstIterator ni(base_boxes.realBegin());
        ni != base_boxes.realEnd(); ++ni) {

      const Box& base_box = *ni;

      // Grow the base_box and put it in the head refinement ratio.
      Box box = base_box;
      BoxContainer grown_boxes;

      if (base.getGridGeometry()->getNumberBlocks() == 1 ||
          base.getGridGeometry()->hasIsotropicRatios()) {
         box.grow(getConnectorWidth());

         if (head_is_finer) {
            box.refine(getRatio());
         } else if (base_is_finer) {
            box.coarsen(getRatio());
         }
         grown_boxes.pushBack(box);
      } else {
         BoxUtilities::growAndAdjustAcrossBlockBoundary(grown_boxes,
            box,
            base.getGridGeometry(),
            base.getRefinementRatio(),
            getRatio(),
            getConnectorWidth(),
            head_is_finer,
            base_is_finer);
      }

      for (BoxContainer::iterator b_itr = grown_boxes.begin();
           b_itr != grown_boxes.end(); ++b_itr) {

         // Add found overlaps to neighbor set for box.
         rbbt.findOverlapBoxes(nabrs_for_box,
            *b_itr,
            head.getRefinementRatio(),
            true);
      }
      if (discard_self_overlap) {
         nabrs_for_box.order();
         nabrs_for_box.erase(base_box);
      }
      if (!nabrs_for_box.empty()) {
         insertNeighbors(nabrs_for_box, base_box.getBoxId());
         nabrs_for_box.clear();
      }

   }

   if (sanity_check_method_postconditions) {
      assertConsistencyWithBase();
      assertConsistencyWithHead();
   }

   t_find_overlaps_rbbt->stop();
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
