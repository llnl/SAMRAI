/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with mapping Connectors.
 *
 ************************************************************************/
#include "SAMRAI/hier/BaseConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace hier {

const int
BaseConnectorAlgorithm::BASE_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH = 1000;

/*
 ***********************************************************************
 * Receive messages and unpack info sent from other processes.
 ***********************************************************************
 */
void
BaseConnectorAlgorithm::setupCommunication(
   tbox::AsyncCommPeer<int> *& all_comms,
   tbox::AsyncCommStage& comm_stage,
   const tbox::SAMRAI_MPI& mpi,
   const std::set<int>& incoming_ranks,
   const std::set<int>& outgoing_ranks,
   const std::shared_ptr<tbox::Timer>& mpi_wait_timer,
   int& operation_mpi_tag,
   bool print_steps) const
{
   /*
    * Set up communication mechanism (and post receives).  We lump all
    * communication objects into one array, all_comms.  all_comms is
    * ordered with the incoming first and the outgoing afterward.
    */
   comm_stage.setCommunicationWaitTimer(mpi_wait_timer);
   const int n_comm = static_cast<int>(
         incoming_ranks.size() + outgoing_ranks.size());
   if (n_comm > 0) {
      all_comms = new tbox::AsyncCommPeer<int>[n_comm];
   }

   const int tag0 = ++operation_mpi_tag;
   const int tag1 = ++operation_mpi_tag;

   std::set<int>::const_iterator owneri;
   size_t comm_idx = 0;
   for (owneri = incoming_ranks.begin();
        owneri != incoming_ranks.end();
        ++owneri, ++comm_idx) {
      const int peer_rank = *owneri;
      tbox::AsyncCommPeer<int>& incoming_comm = all_comms[comm_idx];
      incoming_comm.initialize(&comm_stage);
      incoming_comm.setPeerRank(peer_rank);
      incoming_comm.setMPI(mpi);
      incoming_comm.setMPITag(tag0, tag1);
      incoming_comm.limitFirstDataLength(
         BASE_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH);
      incoming_comm.beginRecv();
      if (print_steps) {
         tbox::plog << "Receiving from " << incoming_comm.getPeerRank()
                    << std::endl;
      }
      if (incoming_comm.isDone()) {
         incoming_comm.pushToCompletionQueue();
      }
   }

   for (owneri = outgoing_ranks.begin();
        owneri != outgoing_ranks.end();
        ++owneri, ++comm_idx) {
      const int peer_rank = *owneri;
      tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[comm_idx];
      outgoing_comm.initialize(&comm_stage);
      outgoing_comm.setPeerRank(peer_rank);
      outgoing_comm.setMPI(mpi);
      outgoing_comm.setMPITag(tag0, tag1);
      outgoing_comm.limitFirstDataLength(
         BASE_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH);
      if (print_steps) {
         tbox::plog << "Sending to " << outgoing_comm.getPeerRank()
                    << std::endl;
      }
   }
}

/*
 ***********************************************************************
 * privateBridge/Modify_findOverlapsForOneProcess() cached some discovered
 * remote neighbors into send_mesg.  packReferencedNeighbors() packs the
 * message with this information.
 *
 * privateBridge/Modify_findOverlapsForOneProcess() placed neighbor data
 * in referenced_new_head_nabrs and referenced_new_base_nabrs rather than
 * directly into send_mesg.  This method packs the referenced neighbors.
 ***********************************************************************
 */
void
BaseConnectorAlgorithm::packReferencedNeighbors(
   std::vector<int>& send_mesg,
   int idx_offset_to_ref,
   const BoxContainer& referenced_new_head_nabrs,
   const BoxContainer& referenced_new_base_nabrs,
   const tbox::Dimension& dim,
   bool print_steps) const
{
   /*
    * Fill the messages's reference section with all the neighbors
    * that have been referenced.
    */
   const int offset = send_mesg[idx_offset_to_ref] =
         static_cast<int>(send_mesg.size());
   const int n_referenced_nabrs = static_cast<int>(
         referenced_new_head_nabrs.size() + referenced_new_base_nabrs.size());
   const int reference_section_size =
      2 + n_referenced_nabrs * Box::commBufferSize(dim);
   send_mesg.insert(
      send_mesg.end(),
      reference_section_size,
      -1);
   int* ptr = &send_mesg[offset];
   *(ptr++) = static_cast<int>(referenced_new_base_nabrs.size());
   *(ptr++) = static_cast<int>(referenced_new_head_nabrs.size());
   for (BoxContainer::const_iterator ni = referenced_new_base_nabrs.begin();
        ni != referenced_new_base_nabrs.end(); ++ni) {
      const Box& box = *ni;
      box.putToIntBuffer(ptr);
      ptr += Box::commBufferSize(dim);
   }
   for (BoxContainer::const_iterator ni = referenced_new_head_nabrs.begin();
        ni != referenced_new_head_nabrs.end(); ++ni) {
      const Box& box = *ni;
      box.putToIntBuffer(ptr);
      ptr += Box::commBufferSize(dim);
   }
   if (print_steps) {
      tbox::plog << "sending " << referenced_new_base_nabrs.size()
                 << " referenced_new_base_nabrs:\n"
                 << referenced_new_base_nabrs.format("  ") << std::endl
                 << "sending " << referenced_new_head_nabrs.size()
                 << " referenced_new_head_nabrs:\n"
                 << referenced_new_head_nabrs.format("  ") << std::endl;
   }

   TBOX_ASSERT(ptr == &send_mesg[send_mesg.size() - 1] + 1);
}

/*
 ***********************************************************************
 * Receive messages and unpack info sent from other processes.
 ***********************************************************************
 */
void
BaseConnectorAlgorithm::receiveAndUnpack(
   Connector& new_base_to_new_head,
   Connector* new_head_to_new_base,
   const std::set<int>& incoming_ranks,
   tbox::AsyncCommPeer<int>* all_comms,
   tbox::AsyncCommStage& comm_stage,
   const std::shared_ptr<tbox::Timer>& receive_and_unpack_timer,
   bool print_steps) const
{
   receive_and_unpack_timer->start();
   /*
    * Receive and unpack messages.
    */
   while (comm_stage.hasCompletedMembers() || comm_stage.advanceSome()) {

      tbox::AsyncCommPeer<int>* peer =
         CPP_CAST<tbox::AsyncCommPeer<int> *>(comm_stage.popCompletionQueue());

      TBOX_ASSERT(peer != 0);

      if ((size_t)(peer - all_comms) < incoming_ranks.size()) {
         // Receive from this peer.
         if (print_steps) {
            tbox::plog << "Received from " << peer->getPeerRank()
                       << std::endl;
         }
         unpackDiscoveryMessage(
            peer,
            new_base_to_new_head,
            new_head_to_new_base,
            print_steps);
      } else {
         // Sent to this peer.  No follow-up needed.
         if (print_steps) {
            tbox::plog << "Sent to " << peer->getPeerRank() << std::endl;
         }
      }

   }

   receive_and_unpack_timer->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BaseConnectorAlgorithm::unpackDiscoveryMessage(
   const tbox::AsyncCommPeer<int>* incoming_comm,
   Connector& new_base_to_new_head,
   Connector* new_head_to_new_base,
   bool print_steps) const
{
   const int sender = incoming_comm->getPeerRank();
   const int* ptr = incoming_comm->getRecvData();

#ifdef DEBUG_CHECK_ASSERTIONS
   const int msg_size = incoming_comm->getRecvSize();
   const int* ptr_end = ptr + msg_size;
#endif
   const int rank = new_base_to_new_head.getMPI().getRank();

   const tbox::Dimension& dim(new_base_to_new_head.getRatio().getDim());

   Box tmp_box(dim);

   const int box_com_buffer_size = Box::commBufferSize(dim);
   /*
    * Content of send_mesg, constructed largely in
    * privateBridge/Modify_findOverlapsForOneProcess() and
    * packReferencedNeighbors():
    *
    * -neighbor-removal section cached in neighbor_removal_mesg.
    * - offset to the reference section (see below)
    * - number of new_base boxes for which overlaps are found
    * - number of new_head boxes for which overlaps are found
    *   - id of new_base/new_head box
    *   - number of neighbors found for new_base/new_head box.
    *     - owner and local indices of neighbors found (unsorted).
    *       Boxes of found neighbors are given in the
    *       reference section of the message.
    * - reference section: all the boxes referenced as
    *   neighbors (accumulated in referenced_new_base_nabrs
    *   and referenced_new_head_nabrs).
    *   - number of referenced new_base neighbors
    *   - number of referenced new_head neighbors
    *   - referenced new_base neighbors
    *   - referenced new_head neighbors
    */

   // Unpack neighbor-removal section.
   // TODO: Get rid of 2 unused values, making sure adjustments in message sizes
   // is correct.
   const int num_removed_boxes = *(ptr++);
   for (int ii = 0; ii < num_removed_boxes; ++ii) {
      const LocalId id_gone(*(ptr++));
      // Skip unneeded value.
      ++ptr;
      const int number_affected = *(ptr++);
      const Box box_gone(dim, GlobalId(id_gone, sender));
      if (print_steps) {
         tbox::plog << "Box " << box_gone
                    << " removed, affecting " << number_affected
                    << " boxes." << std::endl;
      }
      for (int iii = 0; iii < number_affected; ++iii) {
         const LocalId id_affected(*(ptr++));
         // Skip unneeded block id.
         ++ptr;
         BoxId affected_nbrhd(id_affected, rank);
         if (print_steps) {
            tbox::plog << " Removing " << box_gone
                       << " from nabr list for " << id_affected
                       << std::endl;
         }
         TBOX_ASSERT(new_base_to_new_head.hasLocalNeighbor(
               affected_nbrhd,
               box_gone));
         new_base_to_new_head.eraseNeighbor(box_gone, affected_nbrhd);
      }
      TBOX_ASSERT(ptr != ptr_end);
   }

   // Get the referenced neighbor Boxes.
   bool ordered = true;
   BoxContainer referenced_new_base_nabrs(ordered);
   BoxContainer referenced_new_head_nabrs(ordered);
   const int offset = *(ptr++);
   const int n_new_base_boxes = *(ptr++);
   const int n_new_head_boxes = *(ptr++);
   const int* ref_box_ptr = incoming_comm->getRecvData() + offset;
   const int n_reference_new_base_boxes = *(ref_box_ptr++);
   const int n_reference_new_head_boxes = *(ref_box_ptr++);

   TBOX_ASSERT(new_head_to_new_base != 0 || n_new_head_boxes == 0);

#ifdef DEBUG_CHECK_ASSERTIONS
   const int correct_msg_size = offset
      + 2 /* counters of new_head and new_base reference boxes */
      + Box::commBufferSize(dim) * n_reference_new_base_boxes
      + Box::commBufferSize(dim) * n_reference_new_head_boxes;
   TBOX_ASSERT(msg_size == correct_msg_size);
#endif

   // Extract referenced boxes from message.
   for (int ii = 0; ii < n_reference_new_base_boxes; ++ii) {
      tmp_box.getFromIntBuffer(ref_box_ptr);
      referenced_new_base_nabrs.insert(
         referenced_new_base_nabrs.end(),
         tmp_box);
      ref_box_ptr += box_com_buffer_size;
   }
   for (int ii = 0; ii < n_reference_new_head_boxes; ++ii) {
      tmp_box.getFromIntBuffer(ref_box_ptr);
      referenced_new_head_nabrs.insert(
         referenced_new_head_nabrs.end(),
         tmp_box);
      ref_box_ptr += box_com_buffer_size;
   }
   TBOX_ASSERT(ref_box_ptr == ptr_end);

   if (print_steps) {
      tbox::plog << "received " << n_reference_new_base_boxes
                 << " referenced_new_base_nabrs:\n"
                 << referenced_new_base_nabrs.format("  ") << std::endl
                 << "received " << n_reference_new_head_boxes
                 << " referenced_new_head_nabrs:\n"
                 << referenced_new_head_nabrs.format("  ") << std::endl;
   }

   /*
    * Unpack neighbor data for new_head neighbors of new_base boxes
    * and new_base neighbors of new_head boxes.  The neighbor info
    * given includes only block and local index.  Refer to
    * reference data to get the box info.
    */
   for (int ii = 0; ii < n_new_base_boxes; ++ii) {
      const LocalId local_id(*(ptr++));
      const BlockId block_id(*(ptr++));
      const BoxId new_base_box_id(local_id, rank);
      const int n_new_head_nabrs_found = *(ptr++);
      // Add received neighbors to Box new_base_box_id.
      if (n_new_head_nabrs_found > 0) {
         Connector::NeighborhoodIterator base_box_itr =
            new_base_to_new_head.makeEmptyLocalNeighborhood(new_base_box_id);
         BoxId box_id;
         for (int j = 0; j < n_new_head_nabrs_found; ++j) {
            box_id.getFromIntBuffer(ptr);
            tmp_box.setId(box_id);
            ptr += BoxId::commBufferSize();
            BoxContainer::const_iterator na =
               referenced_new_head_nabrs.find(tmp_box);
            TBOX_ASSERT(na != referenced_new_head_nabrs.end());
            const Box& new_head_nabr = *na;
            new_base_to_new_head.insertLocalNeighbor(
               new_head_nabr,
               base_box_itr);
         }
      }
   }
   for (int ii = 0; ii < n_new_head_boxes; ++ii) {
      const LocalId local_id(*(ptr++));
      const BlockId block_id(*(ptr++));
      const BoxId new_head_box_id(local_id, rank);
      const int n_new_base_nabrs_found = *(ptr++);
      // Add received neighbors to Box new_head_box_id.
      if (n_new_base_nabrs_found > 0) {
         Connector::NeighborhoodIterator base_box_itr =
            new_head_to_new_base->makeEmptyLocalNeighborhood(new_head_box_id);
         BoxId box_id;
         for (int j = 0; j < n_new_base_nabrs_found; ++j) {
            box_id.getFromIntBuffer(ptr);
            tmp_box.setId(box_id);
            ptr += BoxId::commBufferSize();
            BoxContainer::const_iterator na =
               referenced_new_base_nabrs.find(tmp_box);
            TBOX_ASSERT(na != referenced_new_base_nabrs.end());
            const Box& new_base_nabr = *na;
            new_head_to_new_base->insertLocalNeighbor(
               new_base_nabr,
               base_box_itr);
         }
      }
   }
}

}
}
