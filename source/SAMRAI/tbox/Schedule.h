/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Schedule of communication transactions between processors
 *
 ************************************************************************/
#ifndef included_tbox_Schedule
#define included_tbox_Schedule

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/ScheduleOpsStrategy.h"
#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/tbox/TransactionFuseable.h"
#include "SAMRAI/tbox/KernelFuser.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"

#include <iostream>
#include <map>
#include <list>
#include <memory>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class Schedule is used to construct and execute a set of
 * data communication transactions.  Each transaction represents some
 * data dependency and exchange between two processors, or locally
 * involving a single processor.
 *
 * Once a communication schedule is constructed, transactions are
 * provided to the schedule, using either the addTransaction() method
 * or the appendTransaction() method.  The schedule is then executed
 * forcing the communication, either interprocessor or local to occur.
 * The basic idea behind the schedule is that it enables the cost of
 * assembling communication dependencies and data transfers to be
 * amortized over many communication phases.
 *
 * Note that since the transactions are stored in lists, the "add" and
 * "append" mimick the semantics of the List class.  That is,
 * addTransaction() will put the transaction at the head of the list,
 * while appendTransaction() will put the transaction at the end of
 * the list.  This flexibility is provided for situations where the
 * order of transaction execution matters.  The transactions will be
 * executed in the order in which they appear in the list.
 *
 * @see Transaction
 */

class Schedule
{
public:
   /*!
    * @brief Create an empty schedule with no transactions.
    */
   Schedule();

   /*!
    * @brief The destructor deletes the schedule and all associated storage.
    *
    * Note that the schedule can not be deleted during a communication
    * phase; this will result in an assertion being thrown.
    *
    * @pre !allocatedCommunicationObjects()
    */
   ~Schedule();

   /*!
    * @brief Add a data transaction to the head of the list of transactions
    * in the schedule.
    *
    * The transaction must involve the local processor as either a
    * source or destination or both.  If the transaction does not
    * include the local processor, then the transaction is not placed
    * on the schedule.
    *
    * @param transaction  std::shared_ptr to transaction added to the
    * schedule.
    */
   void
   addTransaction(
      const std::shared_ptr<Transaction>& transaction);

   /*!
    * @brief Append a data transaction to the tail of the list of
    * transactions in the schedule.
    *
    * The transaction must involve the local processor as either a
    * source or destination or both.  If the transaction does not
    * include the local processor, then the transaction will not be
    * not placed on the schedule.
    *
    * @param transaction  std::shared_ptr to transaction appended to the
    * schedule.
    */
   void
   appendTransaction(
      const std::shared_ptr<Transaction>& transaction);

   /*!
    * @brief Return number of send transactions in the schedule.
    */
   int
   getNumSendTransactions(
      const int rank) const;

   /*!
    * @brief Return number of receive transactions in the schedule.
    */
   int
   getNumRecvTransactions(
      const int rank) const;

   /*!
    * @brief Return number of local transactions (non-communicating)
    * in the schedule.
    */
   int
   getNumLocalTransactions() const
   {
      return static_cast<int>(d_local_set.size());
   }

   /*!
    * @brief Set the MPI communicator used for communication.
    *
    * By default the communicator returned by
    * SAMRAI_MPI::getCommunicator() at the time the Schedule
    * constructor is called is used.  This method may be used to
    * override the default communicator.
    */
   void
   setMPI(
      const SAMRAI_MPI& mpi)
   {
      d_mpi = mpi;
   }

   /*!
    * @brief Specify MPI tag values to use in communication.
    *
    * If you are using MPI communicators that are not sufficiently
    * isolated from other communications, you can specify distinct
    * tags to avoid message mix-ups.  Up to two messages are sent from
    * each communicating pairs.  Specify two distinct tags.
    *
    * @pre first_tag >= 0
    * @pre second_tag >= 0
    */
   void
   setMPITag(
      const int first_tag,
      const int second_tag)
   {
      TBOX_ASSERT(first_tag >= 0);
      TBOX_ASSERT(second_tag >= 0);
      d_first_tag = first_tag;
      d_second_tag = second_tag;
   }

   /*!
    * @brief Specify the message length (in bytes) used in the first
    * message when the receiving processor cannot determine the
    * message length.
    *
    * A receiving transaction cannot always determine the length of data
    * it is receiving.  In this case, first a message of length
    * first_message_length is sent.  The first message includes the full
    * length of the data that is being sent so the receiver can
    * determine if a second message, with the remainder of the data,
    * is needed.  If the length of the data being sent is less than
    * first_message_length, only a single message is sent.
    *
    * Setting a first_message_length too small creates the need for many
    * second messages.  Setting it to an excessive length wastes memory.
    *
    * This two message protocol was designed based on the observed
    * platform-dependent MPI implementation feature where a message of
    * a small length takes the same time as sending an empty message.
    * A protocol of sending the message length first and then a second
    * message with the data does not exploit this property.  The
    * message protocol SAMRAI uses of sending some small amount of
    * data with the first message does exploit this property and will
    * save the cost of always communicating two messages for small
    * messages.  Currently SAMRAI does not attempt to automatically
    * determine this message length however that may be added in
    * future versions.
    *
    * first_message_length defaults to 1000.
    *
    * @param first_message_length length (in bytes) of first message
    *
    * @pre first_message_length > 0
    */
   void
   setFirstMessageLength(
      int first_message_length)
   {
      TBOX_ASSERT(first_message_length > 0);
      d_first_message_length = static_cast<size_t>(first_message_length);
   }

   /*!
    * @brief Perform the communication described by the schedule.
    *
    * This method is simply a <TT>beginCommunication()</TT> followed by
    * <TT>finalizeCommunication()</TT>.
    */
   void
   communicate();

   /*!
    * @brief Begin the communication process but do not deliver data to the
    * transaction objects.
    *
    * This method must be followed by a call to <TT>finalizeCommunication()</TT>
    * in order to complete the communication.
    */
   void
   beginCommunication();

   /*!
    * @brief Finish the communication and deliver the messages.
    *
    * This method must follow a communication begun by
    * <TT>beginCommunication()</TT> and communication is completed
    * when this method returns.
    */
   void
   finalizeCommunication();

   /*!
    * @brief Set whether to unpack messages in a deterministic order.
    *
    * By default message unpacking is ordered by receive time, which
    * is not deterministic.  If your results are dependent on unpack
    * ordering and you want deterministic results, set this flag to
    * true.
    *
    * @param [in] flag
    */
   void setDeterministicUnpackOrderingFlag(bool flag)
   {
      d_unpack_in_deterministic_order = flag;
   }

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named "tbox::Schedule::*",
    * where the third field is the specific steps performed
    * by the Schedule.  You can override the first two
    * fields with this method.  Conforming to the timer
    * naming convention, timer_prefix should have the form
    * "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   /*!
    * @brief Print class data to the specified output stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Returns true if the communication objects have been allocated.
    */
   bool
   allocatedCommunicationObjects()
   {
      return (d_recv_coms.size() > 0 && d_send_coms.size() > 0);
   }

   /*!
    * @brief Get the name of this object.
    */
   const std::string
   getObjectName() const
   {
      return "Schedule";
   }

   /*!
    * @brief Set a pointer to a ScheduleOpsStrategy object
    *
    * @param strategy   Pointer to a concrete instance of ScheduleOpsStrategy
    */
   void setScheduleOpsStrategy(ScheduleOpsStrategy* strategy)
   {
      d_ops_strategy = strategy;
   }

private:
   void
   allocateCommunicationObjects();
   void
   deallocateCommunicationObjects()
   {
      d_send_coms.clear();
      d_recv_coms.clear();
   }

   void
   postReceives();
   void
   postSends();
   void
   performLocalCopies();
   void
   processCompletedCommunications();
   void
   deallocateSendBuffers();

   Schedule(
      const Schedule&);                 // not implemented
   Schedule&
   operator = (
      const Schedule&);                 // not implemented

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      TimerStruct& timers(s_static_timers[s_default_timer_prefix]);
      getAllTimers(s_default_timer_prefix, timers);
   }

   /*!
    * @brief Read input data from input database and initialize class members.
    */
   void
   getFromInput();

   /*
    * @brief Transactions in this schedule.
    *
    * Three containers of transactions are maintained based on
    * source and destination.
    */
   typedef std::map<int, std::list<std::shared_ptr<Transaction> > > TransactionSets;
   TransactionSets d_send_sets;
   TransactionSets d_recv_sets;

   TransactionSets d_send_sets_fuseable;
   TransactionSets d_recv_sets_fuseable;

   StagedKernelFusers* d_send_fusers{nullptr};
   StagedKernelFusers* d_recv_fusers{nullptr};

   /*
    * @brief Transactions where the source and destination are the
    * local process.
    */
   std::list<std::shared_ptr<Transaction> > d_local_set;

   /*
    * @brief Fuseable transactions where the source and destination are the
    * local process.
    */
   std::list<std::shared_ptr<Transaction> > d_local_set_fuseable;

   StagedKernelFusers* d_local_fusers{nullptr};

   ScheduleOpsStrategy* d_ops_strategy{nullptr};

   //@{ @name High-level asynchronous messages passing objects

   /*!
    * @brief Peer-to-peer communication objects, one for each outgoing
    * and each incoming message.
    *
    * d_coms is typed for byte sending because our data is of
    * unknown mixed type.
    */
   using CommMap = std::map<int, std::shared_ptr<AsyncCommPeer<char>>>;
   CommMap d_send_coms;
   CommMap d_recv_coms;

   /*!
    * @brief Stage for advancing communication operations to
    * completion.
    */
   AsyncCommStage d_com_stage;

   //@}

   /*!
    * @brief MPI wrapper.
    *
    * See setMPI().
    */
   SAMRAI_MPI d_mpi;
   /*!
    * @brief Tag value for first message in schedule execution.
    *
    * See setMPITag().
    */
   int d_first_tag;
   /*!
    * @brief Tag value for second message in schedule execution.
    *
    * See setMPITag().
    */
   int d_second_tag;

   /*!
    * @brief If send message length is unknown, use
    * d_first_message_length for the first message's data size.
    */
   size_t d_first_message_length;

   /*!
    * @brief Whether to unpack messages in a deterministic order.
    *
    * @see setDeterministicUnpackOrderingFlag()
    */
   bool d_unpack_in_deterministic_order;

   static const int s_default_first_tag;
   static const int s_default_second_tag;
   static const size_t s_default_first_message_length;

   //@{
   //! @name Timer data for Schedule class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each Schedule object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      std::shared_ptr<Timer> t_communicate;
      std::shared_ptr<Timer> t_begin_communication;
      std::shared_ptr<Timer> t_finalize_communication;
      std::shared_ptr<Timer> t_post_receives;
      std::shared_ptr<Timer> t_post_sends;
      std::shared_ptr<Timer> t_process_incoming_messages;
      std::shared_ptr<Timer> t_MPI_wait;
      std::shared_ptr<Timer> t_pack_stream;
      std::shared_ptr<Timer> t_unpack_stream;
      std::shared_ptr<Timer> t_local_copies;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   static char s_ignore_external_timer_prefix;

   /*!
    * @brief Structure of timers in s_static_timers, matching this
    * object's timer prefix.
    */
   TimerStruct* d_object_timers;

   /*!
    * @brief Get all the timers defined in TimerStruct.  The timers
    * are named with the given prefix.
    */
   static void
   getAllTimers(
      const std::string& timer_prefix,
      TimerStruct& timers);

   //@}

   static StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
