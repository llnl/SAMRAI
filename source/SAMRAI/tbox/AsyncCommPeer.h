/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Staged peer-to-peer communication.
 *
 ************************************************************************/
#ifndef included_tbox_AsyncCommPeer
#define included_tbox_AsyncCommPeer

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <algorithm>
#include <cstring>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Supports point-to-point asynchronous communication
 * operations that can overlap other communcations using a
 * AsyncCommStage.
 *
 * To use this class:
 *
 * -# Construct an object with a stage or @c initialize() it with a
 *    stage.
 *
 * -# Set the MPI communicator and tag values for use in the MPI
 *    nonblocking message passing functions.
 *
 * -# Set the peer rank using setPeerRank().
 *
 * -# Call the beginSend() and beginRecv() functions.  Make sure
 *    there is a matching call on the peer processor.
 *
 * -# Use the stage to advance the message passing operations.  The
 *    stage tells you which members have completed their operations.
 *
 * -# Received data can be accessed using getRecvSize() and getRecvData().
 *
 * It is not neccesary to know how much data is being sent.  This
 * class will figure that out (using overhead data attached to the
 * actual user data being communicated).
 *
 * To check on sends and receives without using the stage, use
 * checkSend() and checkRecv() or proceedToNextWait().  To wait for
 * the current operation to complete, use completeCurrentOperation().
 */
template<class TYPE>
class AsyncCommPeer:public AsyncCommStage::Member
{

private:
   //! @brief Operations users would want to do.
   enum BaseOp { undefined,
                 send,
                 recv };
   //! @brief Tasks used to complete a base operation.
   enum TaskOp { send_start,
                 send_check,
                 recv_start,
                 recv_check0,
                 recv_check,
                 none };

public:
   /*!
    * @brief Default constructor does not set up anything.
    * You must initialize() the object before using it.
    */
   AsyncCommPeer();

   /*!
    * @brief Construct staged communication object.
    *
    * @see initialize().
    */
   explicit AsyncCommPeer(
      AsyncCommStage* stage,
      AsyncCommStage::Handler* handler = 0);

   /*!
    * @brief Destructor.
    *
    * @pre isDone()
    */
   virtual ~AsyncCommPeer();

   /*!
    * @brief Initialize the object.
    *
    * Attach self to the given stage and set the Handler.
    *
    * @param stage The stage handling communicaiton requests for the object.
    * @param handler Optional pointer to user-defined data.
    *
    * @pre isDone()
    */
   void
   initialize(
      AsyncCommStage* stage,
      AsyncCommStage::Handler* handler = 0);

   //@{
   //! @name Define the peer relationship.

   /*!
    * @brief Set the peer rank.
    *
    * We assume (without checking) that the peer processor does
    * a complementary call using this processor as its peer.
    * If this assumption is wrong, there will likely be
    * communication errors.
    *
    * @param peer_rank
    *
    * @pre isDone()
    */
   void
   setPeerRank(
      int peer_rank);

   /*!
    * @brief Rank of peer process.
    */
   int
   getPeerRank() const
   {
      return d_peer_rank;
   }

   /*!
    * @brief Limit the data length in first message of a communication.
    *
    * Sometimes, the receiver does not know how much data is being
    * sent.  In such cases, the receiver would assume the data is no
    * more than a size that the sender and receiver have agreed upon
    * before hand, @c max_first_data_len.  If the message is longer,
    * the sender sends it in two parts.  The first message contains no
    * more than @c max_first_data_len items and indicates the true
    * length of the data.  A second message is used if needed.
    *
    * The receiver must set a first data length no less than
    * the sender's or risk a fatal MPI "truncated message" error.
    *
    * The limit is set to one by default, meaning that any message
    * more than one item long would be sent in two parts.
    */
   void
   limitFirstDataLength(
      size_t max_first_data_len)
   {
      d_max_first_data_len = s_int_max > max_first_data_len ?
                             max_first_data_len : s_int_max;
   }
   //@}

   /*!
    * @brief Set the MPI tags used for communication.
    *
    * The @c primary_tag is used for every message.  If a message
    * must be sent in two parts (see limitFirstDataLength()),
    * the secondary tag is used for the follow-up message.  It is
    * advisable to use different values for the two tags to ensure
    * that the second message is not mistaken for the first.
    *
    * @attention This class is NOT (and cannot be) responsible for
    * ensuring that the MPI communicator and tag are sufficient to
    * select the correct messages.  Please specify appropriate values
    * for the MPI communicator and tag.  Very elusive bugs can occur
    * if incorrect messages are received.
    *
    * @param primary_tag
    * @param secondary_tag
    *
    * @pre isDone()
    */
   void
   setMPITag(
      const int primary_tag,
      const int secondary_tag);

   /*!
    * @brief Set the MPI object used for communications.
    *
    * @attention This class is NOT (and cannot be) responsible for
    * ensuring that the MPI object and tag are sufficient to
    * select the correct messages.  Please specify appropriate values
    * for the MPI object and tag.  Very elusive bugs can occur
    * if incorrect messages are received.  To be safe, it is best
    * to create a new communicator to avoid interference with other
    * communications within SAMRAI.
    *
    * @param mpi
    *
    * @pre isDone()
    */
   void
   setMPI(
      const SAMRAI_MPI& mpi);

   /*!
    * @brief Gets the MPI object used for communications.
    */
   const SAMRAI_MPI&
   getMPI() const
   {
      return d_mpi;
   }

#ifdef HAVE_UMPIRE
   /*!
    * Set Umpire allocator for internal buffer allocations.
    */
   void
   setAllocator(umpire::TypedAllocator<char> allocator)
   {
      d_allocator = allocator;
   }
#endif

   /*!
    * @brief Next task in a current communication operation.
    */
   TaskOp
   getNextTaskOp() const
   {
      return d_next_task_op;
   }

   /*!
    * @brief Operation being performed.
    */
   BaseOp
   getBaseOp() const
   {
      return d_base_op;
   }

   /*!
    * @brief Primary tag for the first message.
    */
   int
   getPrimaryTag() const
   {
      return d_tag0;
   }

   /*!
    * @brief Secondary tag for the first message.
    */
   int
   getSecondaryTag() const
   {
      return d_tag1;
   }

   //@{

   /*!
    * @name Communication methods
    */

   /*!
    * @brief Begin a send communication.
    *
    * This must be paired with a beginRecv() by the peer processor.
    *
    * On return, the data in @b buffer would have been copied so it is
    * safe to deallocate or modify @c buffer.
    *
    * @param buffer
    * @param size
    *
    * @param automatic_push_to_completion_queue Whether to automatically
    * push this member onto the completion queue of its stage (where you
    * can retrieve it using the stage's popFromCompletionQueue() method.
    *
    * @return Whether operation is completed.
    *
    * @pre getNextTaskOp() == none
    */
   bool
   beginSend(
      const TYPE* buffer,
      size_t size,
      bool automatic_push_to_completion_queue = false);

   /*!
    * @brief Check the current broadcast communication and complete
    * the broadcast if all MPI requests are fulfilled.
    *
    * If no communication is in progress, this call does nothing.
    *
    * @param automatic_push_to_completion_queue See beginSend().
    *
    * @return Whether operation is completed.
    *
    * @pre getBaseOp() == send
    */
   bool
   checkSend(
      bool automatic_push_to_completion_queue = false);

   /*!
    * @brief Begin a receive communication.
    *
    * This must be paired with a beginSend() by the peer processor.
    *
    * The received data is stored in an internal buffer accessible
    * through getRecvData().  After receiving, the internal buffer is
    * available until the object goes out of scope, another
    * communication operation is initiated or
    * deallocateInternalRecvBuffer() is called.
    *
    * The actual length of data received by the last receive operation
    * is returned by getRecvSize().
    *
    * @param automatic_push_to_completion_queue See beginSend().
    *
    * @return Whether operation is completed.
    *
    * @internal Once everthing is working, we should implement a
    * version that takes an external buffer.  This optimizes out the
    * need to copy from the internal to an external buffer.  But it
    * requires some thought about size checking, especially with
    * regards to how to use the limit on the first message size.
    *
    * @pre getNextTaskOp() == none
    */
   bool
   beginRecv(
      bool automatic_push_to_completion_queue = false);

   /*!
    * @brief Check the current receive communication and complete the
    * receive if the MPI request has been fulfilled.
    *
    * @return Whether operation is completed.
    *
    * @param automatic_push_to_completion_queue See beginSend().
    *
    * @pre getBaseOp() == recv
    */
   bool
   checkRecv(
      bool automatic_push_to_completion_queue = false);

   /*!
    * @brief Return the size of received data.
    *
    * @pre getBaseOp() == recv
    */
   int
   getRecvSize() const;

   /*!
    * @brief Get access to the received data.
    *
    * @pre getBaseOp() == recv
    */
   const TYPE *
   getRecvData() const;

   /*!
    * @brief Clear the receive buffer immediately.
    *
    * This is optional.  The received data remains inside the object
    * the next communicaiton operation starts, or until it is
    * explicitly freed.
    *
    * It is an error to clear the receive buffer in the middle of a
    * communication operation.
    *
    * @pre getNextTaksOp() == none
    */
   void
   clearRecvData();

   /*!
    * @brief Check the current communication and complete it if all
    * MPI requests are fulfilled.
    */
   bool
   proceedToNextWait();

   /*!
    * @brief Whether the last communication operation has finished.
    *
    * This means more than just whether there is outstanding MPI
    * requests such as that returned by hasOutstandingRequests().  The
    * communication is more complex, requiring up to two messages
    * and copying of the received message into the correct buffer.
    * When isDone() returns true, the operation is truly finished.
    */
   bool
   isDone() const;

   /*!
    * @brief Wait for the current operation to complete.
    */
   void
   completeCurrentOperation();

   /*!
    * @brief Whether the current (or last) operation was a send.
    */
   bool isSender() {
      return d_base_op == send;
   }

   /*!
    * @brief Whether the current (or last) operation was a receive.
    */
   bool isReceiver() {
      return d_base_op == recv;
   }

   //@}

   //@{
   //! @name Timers for MPI calls

   /*!
    * @brief Set the send-timer.
    *
    * Set the timer for non-blocking sends.
    * If the timer is null, revert to the default timer named
    * "tbox::AsyncCommPeer::MPI_Isend()".
    *
    * @param [in] send_timer
    */
   void
   setSendTimer(
      const std::shared_ptr<Timer>& send_timer);

   /*!
    * @brief Set the receive-timer.
    *
    * Set the timer for non-blocking receives.
    * If the timer is null, revert to the default timer named
    * "tbox::AsyncCommPeer::MPI_Irecv()".
    *
    * @param [in] recv_timer
    */
   void
   setRecvTimer(
      const std::shared_ptr<Timer>& recv_timer);

   /*!
    * @brief Set the wait-timer.
    *
    * Set the timer for blocking waits.
    * If the timer is null, revert to the default timer named
    * "tbox::AsyncCommPeer::MPI_Waitall()".
    *
    * This timer is used when in this class and is not the same
    * as the timer given to AsyncCommStage (unless the user sets
    * it up that way explicitly).
    *
    * @param [in] wait_timer
    */
   void
   setWaitTimer(
      const std::shared_ptr<Timer>& wait_timer);

   //@}

   /*!
    * @brief For use in debugging.
    */
   void
   logCurrentState(
      std::ostream& co) const;

private:

   /*
    * @brief Data structure for combining integer overhead data along with
    * user data TYPE in the same MPI message.
    */
   union FlexData {
      unsigned int d_uint;
      TYPE d_t;
      FlexData();
   };

   /*
    * @brief Assert that user-set MPI parameters are valid.
    *
    * @pre getPrimaryTag() >= 0 && getSecondaryTag() >= 0
    * @pre getMPI().getCommunicator() != MPI_COMM_NULL
    * @pre getPeerRank() >= 0
    * @pre getPeerRank() != getMPI().getRank() || SAMRAI_MPI::usingMPI()
    */
   void
   checkMPIParams();

   /*!
    * @brief Return the number of FlexData storage space to store a given
    * number of TYPE data.
    */
   size_t
   getNumberOfFlexData(
      size_t number_of_type_data) const;

   /*!
    * @brief Operation disallowed due to primitive internal memory management.
    */
   AsyncCommPeer(
      const AsyncCommPeer& r);

   /*!
    * @brief Operation disallowed by primitive internal memory management.
    */
   AsyncCommPeer&
   operator = (
      const AsyncCommPeer& r);

   void
   resetStatus(
      SAMRAI_MPI::Status& mpi_status)
   {
      mpi_status.MPI_TAG =
         mpi_status.MPI_SOURCE =
            mpi_status.MPI_ERROR = -1;
   }

   /*
    * @brief Resize the internal buffer to hold size FlexData unions.
    *
    * Note this is FlexData's not TYPE's.
    */
   void
   resizeBuffer(
      size_t size);

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Rank of peer process.
    *
    * If negative, peer has not been set.
    */
   int d_peer_rank;

   /*!
    * @brief Operation being performed.
    *
    * At the completion of the operation, d_base_op does NOT get
    * reset to undefined.  Its value is used in getRecvData() and
    * is also useful for debugging.
    */
   BaseOp d_base_op;

   /*!
    * @brief Next task in a current communication operation.
    *
    * If d_next_task_op is none, there is no current communication
    * operation (the last one is completed).
    */
   TaskOp d_next_task_op;

   /*!
    * @brief Max data size of first message when the receiver does
    * not know the message size.
    *
    * @see limitFirstDataLength().
    */
   size_t d_max_first_data_len;

   /*!
    * @brief The full count of data items before splitting up
    * to into the limit of the first chunk.
    */
   size_t d_full_count;

   /*!
    * @brief External data input buffer.
    *
    * Though we keep this around after beginSend() returns, in the
    * current implementation, it is not needed.  It is kept around
    * to help debugging.  After beginSend() returns, the pointer
    * is no longer guaranteed to be valid.
    */
   const TYPE* d_external_buf;

   /*!
    * @brief Internal buffer.
    *
    * Used for sends and receives.  This is typed on FlexData because
    * it must be able to store the templated type TYPE and some integer
    * for overhead data.
    */
   size_t d_internal_buf_size;
   FlexData* d_internal_buf = nullptr;

   /*!
    * @brief Small buffer for count information only
    *
    * Buffer used in instances where MPI messages will only contain
    * metadata information about the size of data and number of data
    * buffers that will be sent in future messages.
    */
   FlexData* d_count_buf = nullptr;

   /*!
    * @brief Buffer for the first receive
    */
   FlexData* d_first_recv_buf = nullptr;

   /*!
    * @brief The SAMRAI_MPI object
    */
   SAMRAI_MPI d_mpi;

   /*!
    * @brief Tag for the first message.
    */
   int d_tag0;
   /*!
    * @brief Tag for the second message.
    */
   int d_tag1;

   /*!
    * @brief Counter of the maximum number of sends that this object has sent.
    */
   unsigned int d_max_sends = 1;

   /*!
    * @brief Whether send completion should be reported when
    * AsyncCommPeer_DEBUG_OUTPUT is defined.
    *
    * This is non-essential data used in debugging.
    */
   bool d_report_send_completion[SAMRAI_MAX_COMM_BUFFERS];

   // Make some temporary variable statuses to avoid repetitious allocations.
   int d_mpi_err;

   /*
    * These are static constexpr values used to handle cases where a
    * buffer passed to beginSend is larger than the maximum value allowed by
    * the int type.  Due to the MPI API using int to indicate the size of
    * buffers, we cannot pass a buffer larger than the int maximum into
    * MPI send/recv calls.  We compute s_int_max as the size at which this
    * class does special operations to split beginSend's buffer into chunks
    * that are passed in a series of MPI sends.
    *
    * First we set s_prelim_max as a value just short of the true INT_MAX,
    * as a safety measure to not push all the way to the limit, and then
    * to set s_int_max, we adjust it to be a multiple of 128, which may be
    * helpful to align with byte sizes of typical data types.
    */
   static constexpr unsigned int s_prelim_max = MathUtilities<int>::getMax() - 256;
   static constexpr unsigned int s_int_max = s_prelim_max - (s_prelim_max % (128*sizeof(FlexData)));

   /*!
    * A constant used as a signal value to indicate that a single MPI
    * send/recv pair has communicated all needed data in one message.
    */
   static constexpr unsigned int s_onemsg_signal = MathUtilities<int>::getMax();

#ifdef HAVE_UMPIRE
   umpire::TypedAllocator<char> d_allocator;
#endif 

   std::shared_ptr<Timer> t_send_timer;
   std::shared_ptr<Timer> t_recv_timer;
   std::shared_ptr<Timer> t_wait_timer;

   static std::shared_ptr<Timer> t_default_send_timer;
   static std::shared_ptr<Timer> t_default_recv_timer;
   static std::shared_ptr<Timer> t_default_wait_timer;

   static StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#include "SAMRAI/tbox/AsyncCommPeer.cpp"

#endif  // included_tbox_AsyncCommPeer
