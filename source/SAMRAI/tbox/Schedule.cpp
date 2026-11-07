/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Schedule of communication transactions between processors
 *
 ************************************************************************/
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Collectives.h"

#include <cstring>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

typedef std::list<std::shared_ptr<Transaction> >::iterator Iterator;
typedef std::list<std::shared_ptr<Transaction> >::const_iterator ConstIterator;

const int Schedule::s_default_first_tag = 0;
const int Schedule::s_default_second_tag = 1;
/*
 * TODO: Set the default first message length to the maximum value
 * possible without incurring any additional cost associated with the
 * MPI communication.  This parameter should be dependent on the MPI
 * implementation.
 */
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
const size_t Schedule::s_default_first_message_length = 0;
#else
const size_t Schedule::s_default_first_message_length = 1000;
#endif

const std::string Schedule::s_default_timer_prefix("tbox::Schedule");
std::map<std::string, Schedule::TimerStruct> Schedule::s_static_timers;
char Schedule::s_ignore_external_timer_prefix('\0');

StartupShutdownManager::Handler
Schedule::s_initialize_finalize_handler(
   Schedule::initializeCallback,
   0,
   0,
   0,
   StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */

Schedule::Schedule():
   d_com_stage(),
   d_mpi(SAMRAI_MPI::getSAMRAIWorld()),
   d_first_tag(s_default_first_tag),
   d_second_tag(s_default_second_tag),
   d_first_message_length(s_default_first_message_length),
   d_unpack_in_deterministic_order(false),
   d_object_timers(0)
{
   getFromInput();
   setTimerPrefix(s_default_timer_prefix);
}

/*
 *************************************************************************
 * Note that the destructor should not be called during a communication
 * phase.
 *************************************************************************
 */
Schedule::~Schedule()
{
   if (allocatedCommunicationObjects()) {
      TBOX_ERROR("Destructing a schedule while communication is pending\n"
         << "leads to lost messages.  Aborting.");
   }
}

/*
 *************************************************************************
 * Add a transaction to the head of a list of data transactions in
 * this schedule. The assignment of the transaction to a list depends
 * on the source and destination processors of the transaction.
 *************************************************************************
 */
void
Schedule::addTransaction(
   const std::shared_ptr<Transaction>& transaction)
{
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   std::shared_ptr<TransactionFuseable> fuseable_transaction{
      std::dynamic_pointer_cast<TransactionFuseable>(transaction)};

   if ((d_mpi.getRank() == src_id) && (d_mpi.getRank() == dst_id)) {
      if (fuseable_transaction) {
         if (!d_local_fusers) {
            d_local_fusers = StagedKernelFusers::getInstance();
         }
         fuseable_transaction->setKernelFuser(d_local_fusers);
         d_local_set_fuseable.push_front(fuseable_transaction);
      } else {
         d_local_set.push_front(transaction);
      }
   } else {
      if (d_mpi.getRank() == dst_id) {
         if (fuseable_transaction) {
            if (!d_recv_fusers) {
               d_recv_fusers = StagedKernelFusers::getInstance();
            }
            fuseable_transaction->setKernelFuser(d_recv_fusers);
            d_recv_sets_fuseable[src_id].push_front(fuseable_transaction);
         } else {
            d_recv_sets[src_id].push_front(transaction);
         }
      } else if (d_mpi.getRank() == src_id) {
         if (fuseable_transaction) {
            if (!d_send_fusers) {
               d_send_fusers = StagedKernelFusers::getInstance();
            }
            fuseable_transaction->setKernelFuser(d_send_fusers);
            d_send_sets_fuseable[dst_id].push_front(fuseable_transaction);
         } else {
            d_send_sets[dst_id].push_front(transaction);
         }
      }
   }
}

/*
 *************************************************************************
 * Append a transaction to the tail of a list of data transactions in
 * this schedule.  The assignment of the transaction to a list depends
 * on the source and destination processors of the transaction.
 *************************************************************************
 */
void
Schedule::appendTransaction(
   const std::shared_ptr<Transaction>& transaction)
{
   const int src_id = transaction->getSourceProcessor();
   const int dst_id = transaction->getDestinationProcessor();

   std::shared_ptr<TransactionFuseable> fuseable_transaction{
      std::dynamic_pointer_cast<TransactionFuseable>(transaction)};

   if ((d_mpi.getRank() == src_id) && (d_mpi.getRank() == dst_id)) {
      if (fuseable_transaction) {
         if (!d_local_fusers) {
            d_local_fusers = StagedKernelFusers::getInstance();
         }
         fuseable_transaction->setKernelFuser(d_local_fusers);
         d_local_set_fuseable.push_back(fuseable_transaction);
      } else {
         d_local_set.push_back(transaction);
      }
   } else {
      if (d_mpi.getRank() == dst_id) {
         if (fuseable_transaction) {
            if (!d_recv_fusers) {
               d_recv_fusers = StagedKernelFusers::getInstance();
            }
            fuseable_transaction->setKernelFuser(d_recv_fusers);
            d_recv_sets_fuseable[src_id].push_back(fuseable_transaction);
         } else {
            d_recv_sets[src_id].push_back(transaction);
         }
      } else if (d_mpi.getRank() == src_id) {
         if (fuseable_transaction) {
            if (!d_send_fusers) {
               d_send_fusers = StagedKernelFusers::getInstance();
            }
            fuseable_transaction->setKernelFuser(d_send_fusers);
            d_send_sets_fuseable[dst_id].push_back(transaction);
         } else {
            d_send_sets[dst_id].push_back(transaction);
         }
      }
   }
}

/*
 *************************************************************************
 * Access number of send transactions.
 *************************************************************************
 */
int
Schedule::getNumSendTransactions(
   const int rank) const
{
   int size = 0;
   TransactionSets::const_iterator mi = d_send_sets.find(rank);
   if (mi != d_send_sets.end()) {
      size += static_cast<int>(mi->second.size());
   }
   mi = d_send_sets_fuseable.find(rank);
   if (mi != d_send_sets_fuseable.end()) {
      size += static_cast<int>(mi->second.size());
   }

   return size;
}

/*
 *************************************************************************
 * Access number of receive transactions.
 *************************************************************************
 */
int
Schedule::getNumRecvTransactions(
   const int rank) const
{
   int size = 0;
   TransactionSets::const_iterator mi = d_recv_sets.find(rank);
   if (mi != d_recv_sets.end()) {
      size += static_cast<int>(mi->second.size());
   }
   mi = d_recv_sets_fuseable.find(rank);
   if (mi != d_recv_sets_fuseable.end()) {
      size += static_cast<int>(mi->second.size());
   }
   return size;
}

/*
 *************************************************************************
 * Perform the communication described by the schedule.
 *************************************************************************
 */
void
Schedule::communicate()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Schedule::communicate: Errant message detected before beginCommunication().");
   }
#endif

   d_object_timers->t_communicate->start();
   beginCommunication();
   finalizeCommunication();
   d_object_timers->t_communicate->stop();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Schedule::communicate: Errant message detected after finalizeCommunication().");
   }
#endif
}

/*
 *************************************************************************
 * Begin communication but do not wait for it to finish.  This routine
 * posts receives, and sends outgoing messages.  Since we do not wait
 * for message completion, use finalizeCommunication() to ensure that
 * communication has finished.
 *************************************************************************
 */
void
Schedule::beginCommunication()
{
   d_object_timers->t_begin_communication->start();
   if (d_ops_strategy) {
      d_ops_strategy->preCommunicate();
   }

   allocateCommunicationObjects();
   postReceives();
   postSends();
   d_object_timers->t_begin_communication->stop();
}

/*
 *************************************************************************
 * Perform the local data copies, complete receive operations and
 * unpack received data into their destinations.
 *************************************************************************
 */
void
Schedule::finalizeCommunication()
{
   d_object_timers->t_finalize_communication->start();
   performLocalCopies();
   processCompletedCommunications();
   deallocateCommunicationObjects();

   if (d_ops_strategy) {
      d_ops_strategy->postCommunicate();
   }

   d_object_timers->t_finalize_communication->stop();
}

/*
 *************************************************************************
 * Post receives.
 *
 * Where message lengths can be locally computed, use the correct
 * message lengths to avoid overheads due to unknown lengths.
 *************************************************************************
 */
void
Schedule::postReceives()
{
   if (d_recv_sets.empty() && d_recv_sets_fuseable.empty()) {
      /*
       * Short cut because some looping logic in this method assumes
       * non-empty d_recv_sets.
       */
      return;
   }


   int rank = d_mpi.getRank();

   /*
    * We loop through d_recv_sets starting with the highest rank that
    * is lower than local process.  We loop backwards, continuing at
    * the opposite end when we run out of sets.  This ordering is in
    * the reverse direction of the message send ordering so that a
    * send posted earlier is paired with a receive that is also posted
    * earlier.
    */
   for (CommMap::reverse_iterator comm_peer(d_recv_coms.lower_bound(rank));
        comm_peer != d_recv_coms.rend();
        ++comm_peer) {
      const int recv_rank = (*comm_peer).first;
      auto& comm = (*comm_peer).second;
      // Compute incoming message size, if possible.
      size_t byte_count = 0;
      bool can_estimate_incoming_message_size = true;

      for (const auto& t : d_recv_sets[recv_rank] ) {
         if (!t->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
            break;
         }
         byte_count += t->computeIncomingMessageSize();
      }

      for (const auto& t: d_recv_sets_fuseable[recv_rank]) {
         if (!t->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
            break;
         }
         byte_count += t->computeIncomingMessageSize();
      }

      // Set AsyncCommPeer to receive known message length.
      if (can_estimate_incoming_message_size) {
         comm->limitFirstDataLength(byte_count);
      }

      // Begin non-blocking receive operation.
      d_object_timers->t_post_receives->start();
      comm->beginRecv();
      if (comm->isDone()) {
         comm->pushToCompletionQueue();
      }
      d_object_timers->t_post_receives->stop();
   }

   CommMap::reverse_iterator stop(d_recv_coms.lower_bound(rank));
   for (CommMap::reverse_iterator comm_peer = d_recv_coms.rbegin(); comm_peer != stop; ++comm_peer) {
      const int recv_rank = (*comm_peer).first;
      auto& comm = (*comm_peer).second;
      // Compute incoming message size, if possible.
      size_t byte_count = 0;
      bool can_estimate_incoming_message_size = true;

      for (const auto& t : d_recv_sets[recv_rank] ) {
         if (!t->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
            break;
         }
         byte_count += t->computeIncomingMessageSize();
      }

      for (const auto& t: d_recv_sets_fuseable[recv_rank]) {
         if (!t->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
            break;
         }
         byte_count += t->computeIncomingMessageSize();
      }

      // Set AsyncCommPeer to receive known message length.
      if (can_estimate_incoming_message_size) {
         comm->limitFirstDataLength(byte_count);
      }

      // Begin non-blocking receive operation.
      d_object_timers->t_post_receives->start();
      comm->beginRecv();
      if (comm->isDone()) {
         comm->pushToCompletionQueue();
      }
      d_object_timers->t_post_receives->stop();
   }

}

/*
 *************************************************************************
 * Allocate the send buffer, pack the data, and initiate the message
 * sends.
 *************************************************************************
 */
void
Schedule::postSends()
{
   d_object_timers->t_post_sends->start();

   if (d_ops_strategy) {
      d_ops_strategy->prePack();
   }

   /*
    * We loop through d_send_sets starting with the first set with
    * rank higher than the local process, continuing at the opposite
    * end when we run out of sets.  This ordering tends to spread out
    * the communication traffic over the entire network to reduce the
    * potential network contention.
    */

   std::map<int, std::shared_ptr<MessageStream> > outgoing_streams;
   std::map<int, bool> defer_stream;

   bool defer_send = false;

   int rank = d_mpi.getRank();
   int mpi_size = d_mpi.getSize();
   int start_rank = (rank+1) % mpi_size;


   for (int ip = start_rank; ip != rank; ip = (ip+1) % mpi_size) {
      if (d_send_coms.find(ip) == d_send_coms.end()) {
         continue; 
      }
      const int peer_rank = ip;
      auto& comm = d_send_coms[ip];

      size_t byte_count = 0;
      bool can_estimate_incoming_message_size = true;
      for (const auto& transaction : d_send_sets[peer_rank]) {
         if (!transaction->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
         }
         byte_count += transaction->computeOutgoingMessageSize();
      }

      for (const auto& transaction : d_send_sets_fuseable[peer_rank]) {
         if (!transaction->canEstimateIncomingMessageSize()) {
            can_estimate_incoming_message_size = false;
         }
         byte_count += transaction->computeOutgoingMessageSize();
      }

      // Pack outgoing data into a message.
      outgoing_streams[peer_rank] = std::make_shared<MessageStream>(
         byte_count,
         MessageStream::Write,
         nullptr,
         true
#ifdef HAVE_UMPIRE
         , AllocatorDatabase::getDatabase()->getStreamAllocator()
#endif
         );

      MessageStream& outgoing_stream = *outgoing_streams[peer_rank];

      d_object_timers->t_pack_stream->start();

      bool have_fuseable = !(d_send_sets_fuseable[peer_rank].empty());
      bool have_nonfuseable = !(d_send_sets[peer_rank].empty());

      if (have_fuseable || have_nonfuseable) { 
         for (const auto& transaction : d_send_sets_fuseable[peer_rank]) {
            transaction->packStream(outgoing_stream);
         }
 
         for (const auto& transaction : d_send_sets[peer_rank]) {
            transaction->packStream(outgoing_stream);
         }
      }

      d_object_timers->t_pack_stream->stop();

      if (can_estimate_incoming_message_size) {
         // Receiver knows message size so set it exactly.
         comm->limitFirstDataLength(byte_count);
      }

      if (d_ops_strategy && !defer_send) {
         defer_send = d_ops_strategy->deferMessageSend();
      }

      if (!defer_send) {

#if defined(HAVE_RAJA)
         parallel_synchronize();
#endif

         // Begin non-blocking send operation.
         comm->beginSend(
            (const char *)outgoing_stream.getBufferStart(),
            outgoing_stream.getCurrentSize());
         defer_stream[peer_rank] = false;
         if (comm->isDone()) {
            comm->pushToCompletionQueue();
         }
      } else {
         defer_stream[peer_rank] = true;
      }

   }

#if defined(HAVE_RAJA)
   bool need_sync = true;
#endif

   if (d_ops_strategy) {
      d_ops_strategy->postPack();
#if defined(HAVE_RAJA)
      need_sync = d_ops_strategy->needSynchronize();
#endif
   }

   if (defer_send) {

#if defined(HAVE_RAJA)
      if (!d_ops_strategy || need_sync) {
         parallel_synchronize();
      }
#endif

      for (int ip = start_rank; ip != rank; ip = (ip+1) % mpi_size) {
         if (d_send_coms.find(ip) == d_send_coms.end()) {
            continue;
         }
         const int peer_rank = ip;
         auto& comm = d_send_coms[ip];

         MessageStream& outgoing_stream = *outgoing_streams[peer_rank];

         // Begin non-blocking send operation.
         if (defer_stream[peer_rank]) {
            comm->beginSend(
               (const char *)outgoing_stream.getBufferStart(),
               outgoing_stream.getCurrentSize());
            if (comm->isDone()) {
               comm->pushToCompletionQueue();
            }
         }
      }
   }

   d_object_timers->t_post_sends->stop();
}

/*
 *************************************************************************
 * Perform all of the local memory-to-memory copies for this processor.
 *************************************************************************
 */
void
Schedule::performLocalCopies()
{
   d_object_timers->t_local_copies->start();

   if (d_ops_strategy) {
      d_ops_strategy->preCopy();
   }

   for (const auto& local : d_local_set_fuseable) {
      local->copyLocalData();
   }

   for (const auto& local : d_local_set) {
      local->copyLocalData();
   }

   // need_sync initialized true unless both sets are empty.
#if defined(HAVE_RAJA) 
   bool need_sync = !(d_local_set_fuseable.empty() && d_local_set.empty());
#endif

   if (d_ops_strategy) {
      d_ops_strategy->postCopy();
      // d_ops_strategy may indicate sync no longer needed.
#if defined(HAVE_RAJA)
      need_sync = d_ops_strategy->needSynchronize();
#endif
   }

#if defined(HAVE_RAJA)
   if (need_sync) {
      parallel_synchronize();
   }
#endif

   d_object_timers->t_local_copies->stop();

}


/*
 *************************************************************************
 * Process completed operations as they come in.  Initially, completed
 * operations are placed in d_completed_comm.  Process these first,
 * then check for next set of completed operations.  Repeat until all
 * operations are completed.
 *
 * Once a receive is completed, put it in a MessageStream for unpacking.
 *************************************************************************
 */
void
Schedule::processCompletedCommunications()
{
   d_object_timers->t_process_incoming_messages->start();

   if (d_unpack_in_deterministic_order) {

      // Unpack in deterministic order.  Wait for receive as needed.
      // Deterministic order is lowest to highest recv rank

      for (auto& comms : d_recv_coms) {
         if (d_ops_strategy) {
            d_ops_strategy->preUnpack();
         }

         auto& completed_comm = comms.second;

         int sender = comms.first;
         TBOX_ASSERT(sender == completed_comm->getPeerRank());
         completed_comm->completeCurrentOperation();
         completed_comm->yankFromCompletionQueue();

         MessageStream incoming_stream(
            static_cast<size_t>(completed_comm->getRecvSize()) * sizeof(char),
            MessageStream::Read,
            completed_comm->getRecvData(),
            false /* don't use deep copy */
#ifdef HAVE_UMPIRE
            , AllocatorDatabase::getDatabase()->getStreamAllocator()
#endif
            );


         d_object_timers->t_unpack_stream->start();
         for (const auto& transaction : d_recv_sets_fuseable[sender]) {
            transaction->unpackStream(incoming_stream);
         }
         for (const auto& transaction : d_recv_sets[sender]) {
            transaction->unpackStream(incoming_stream);
         }

         // need_sync initialized true unless both sets are empty.
#if defined(HAVE_RAJA)
         bool need_sync = !(d_recv_sets_fuseable[sender].empty() &&
                            d_recv_sets[sender].empty());
#endif
         if (d_ops_strategy) {
            d_ops_strategy->postUnpack();
            // d_ops_strategy may indicate sync no longer needed.
#if defined(HAVE_RAJA)
            need_sync = d_ops_strategy->needSynchronize();
#endif
         }

#if defined(HAVE_RAJA)
         if (need_sync) {
            parallel_synchronize();
         }
#endif

         d_object_timers->t_unpack_stream->stop();
         completed_comm->clearRecvData();
      }

      // Complete sends.
      d_com_stage.advanceAll();
      while (d_com_stage.hasCompletedMembers()) {
         d_com_stage.popCompletionQueue();
      }

   } else {

      // Unpack in order of completed receives.

      if (d_ops_strategy) {
         d_ops_strategy->preUnpack();
      }

      bool have_fuseable = false;
      bool have_nonfuseable = false;
      while (d_com_stage.hasCompletedMembers() || d_com_stage.advanceSome()) {

         AsyncCommPeer<char>* completed_comm =
            CPP_CAST<AsyncCommPeer<char> *>(d_com_stage.popCompletionQueue());

         TBOX_ASSERT(completed_comm != 0);
         TBOX_ASSERT(completed_comm->isDone());
         if (!completed_comm->isSender()) {

            const int sender = completed_comm->getPeerRank();

            MessageStream incoming_stream(
               static_cast<size_t>(completed_comm->getRecvSize()) * sizeof(char),
               MessageStream::Read,
               completed_comm->getRecvData(),
               false /* don't use deep copy */
#ifdef HAVE_UMPIRE
               , AllocatorDatabase::getDatabase()->getStreamAllocator()
#endif
               );

            have_fuseable = have_fuseable ||
                            !(d_recv_sets_fuseable[sender].empty());

            have_nonfuseable = have_nonfuseable ||
                               !(d_recv_sets[sender].empty());

            d_object_timers->t_unpack_stream->start();
            for (const auto& transaction : d_recv_sets_fuseable[sender]) {
               transaction->unpackStream(incoming_stream);
            }

            for (const auto& transaction : d_recv_sets[sender]) {
               transaction->unpackStream(incoming_stream);
            }

            d_object_timers->t_unpack_stream->stop();
         } else {
            // No further action required for completed send.
         }
      }

      // need_sync initialized true unless there were no transactions found
      // in the above loop.
#if defined(HAVE_RAJA)
      bool need_sync = have_fuseable || have_nonfuseable;
#endif

      if (d_ops_strategy) {
         d_ops_strategy->postUnpack();
         // d_ops_strategy may indicate sync no longer needed.
#if defined(HAVE_RAJA)
         need_sync = d_ops_strategy->needSynchronize();
#endif
      }

#if defined(HAVE_RAJA)
      if (need_sync) {
         parallel_synchronize();
      }
#endif

      for (auto& comms : d_recv_coms) {
         auto& completed_comm = comms.second;
         completed_comm->clearRecvData();
      }

   }

   d_object_timers->t_process_incoming_messages->stop();
}

/*
 *************************************************************************
 * Allocate communication objects, set them up on the stage and get
 * them ready to send/receive.
 *************************************************************************
 */
void
Schedule::allocateCommunicationObjects()
{
   for (const auto& transaction : d_recv_sets) {
      int rank = transaction.first;

      auto peer = std::make_shared<AsyncCommPeer<char>>();
      peer->initialize(&d_com_stage);
      peer->setPeerRank(rank);
      peer->setMPITag(d_first_tag, d_second_tag);
      peer->setMPI(d_mpi);
      peer->limitFirstDataLength(d_first_message_length);
#ifdef HAVE_UMPIRE
      peer->setAllocator(AllocatorDatabase::getDatabase()->getStreamAllocator());
#endif
      d_recv_coms[rank] = peer;
   }

   for (const auto& transaction : d_recv_sets_fuseable) {
      int rank = transaction.first;
      
      if (d_recv_coms.find(rank) == d_recv_coms.end()) {
         auto peer = std::make_shared<AsyncCommPeer<char>>();
         peer->initialize(&d_com_stage);
         peer->setPeerRank(rank);
         peer->setMPITag(d_first_tag, d_second_tag);
         peer->setMPI(d_mpi);
         peer->limitFirstDataLength(d_first_message_length);
#ifdef HAVE_UMPIRE
         peer->setAllocator(AllocatorDatabase::getDatabase()->getStreamAllocator());
#endif
         d_recv_coms[rank] = peer;
      }
   }

   for (const auto& transaction : d_send_sets) {
      int rank = transaction.first;
      auto peer = std::make_shared<AsyncCommPeer<char>>();

      peer->initialize(&d_com_stage);
      peer->setPeerRank(rank);
      peer->setMPITag(d_first_tag, d_second_tag);
      peer->setMPI(d_mpi);
      peer->limitFirstDataLength(d_first_message_length);
#ifdef HAVE_UMPIRE
      peer->setAllocator(AllocatorDatabase::getDatabase()->getStreamAllocator());
#endif
      d_send_coms[rank] = peer;
   }

   for (const auto& transaction : d_send_sets_fuseable) {
      int rank = transaction.first;
      
      if (d_send_coms.find(rank) == d_send_coms.end()) {
         auto peer = std::make_shared<AsyncCommPeer<char>>();

         peer->initialize(&d_com_stage);
         peer->setPeerRank(rank);
         peer->setMPITag(d_first_tag, d_second_tag);
         peer->setMPI(d_mpi);
         peer->limitFirstDataLength(d_first_message_length);
#ifdef HAVE_UMPIRE
         peer->setAllocator(AllocatorDatabase::getDatabase()->getStreamAllocator());
#endif
         d_send_coms[rank] = peer;
      }
   }
}

/*
 *************************************************************************
 * Print class data to the specified output stream.
 *************************************************************************
 */
void
Schedule::printClassData(
   std::ostream& stream) const
{
   stream << "Schedule::printClassData()" << std::endl;
   stream << "-------------------------------" << std::endl;

   stream << "Number of sends: " << d_send_sets.size() << std::endl;
   stream << "Number of recvs: " << d_recv_sets.size() << std::endl;

   for (TransactionSets::const_iterator ss = d_send_sets.begin();
        ss != d_send_sets.end(); ++ss) {
      const std::list<std::shared_ptr<Transaction> >& send_set = ss->second;
      stream << "Send Set: " << ss->first << std::endl;
      for (ConstIterator send = send_set.begin();
           send != send_set.end(); ++send) {
         (*send)->printClassData(stream);
      }
   }

   for (TransactionSets::const_iterator rs = d_recv_sets.begin();
        rs != d_recv_sets.end(); ++rs) {
      const std::list<std::shared_ptr<Transaction> >& recv_set = rs->second;
      stream << "Recv Set: " << rs->first << std::endl;
      for (ConstIterator recv = recv_set.begin();
           recv != recv_set.end(); ++recv) {
         (*recv)->printClassData(stream);
      }
   }

   stream << "Local Set" << std::endl;
   for (ConstIterator local = d_local_set.begin();
        local != d_local_set.end(); ++local) {
      (*local)->printClassData(stream);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Schedule::getFromInput()
{
   /*
    * - set up debugging flags.
    */
   if (s_ignore_external_timer_prefix == '\0') {
      s_ignore_external_timer_prefix = 'n';
      if (InputManager::inputDatabaseExists()) {
         std::shared_ptr<Database> idb(
            InputManager::getInputDatabase());
         if (idb->isDatabase("Schedule")) {
            std::shared_ptr<Database> sched_db(
               idb->getDatabase("Schedule"));
            s_ignore_external_timer_prefix =
               sched_db->getCharWithDefault("DEV_ignore_external_timer_prefix",
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
Schedule::setTimerPrefix(
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
   d_com_stage.setCommunicationWaitTimer(d_object_timers->t_MPI_wait);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Schedule::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_communicate = TimerManager::getManager()->
      getTimer(timer_prefix + "::communicate()");
   timers.t_begin_communication = TimerManager::getManager()->
      getTimer(timer_prefix + "::beginCommunication()");
   timers.t_finalize_communication = TimerManager::getManager()->
      getTimer(timer_prefix + "::finalizeCommunication()");
   timers.t_post_receives = TimerManager::getManager()->
      getTimer(timer_prefix + "::postReceives()");
   timers.t_post_sends = TimerManager::getManager()->
      getTimer(timer_prefix + "::postSends()");
   timers.t_process_incoming_messages = TimerManager::getManager()->
      getTimer(timer_prefix + "::processIncomingMessages()");
   timers.t_MPI_wait = TimerManager::getManager()->
      getTimer(timer_prefix + "::MPI_wait");
   timers.t_pack_stream = TimerManager::getManager()->
      getTimer(timer_prefix + "::pack_stream");
   timers.t_unpack_stream = TimerManager::getManager()->
      getTimer(timer_prefix + "::unpack_stream");
   timers.t_local_copies = TimerManager::getManager()->
      getTimer(timer_prefix + "::performLocalCopies()");
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
