/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   All-to-one and one-to-all communication using a tree.
 *
 ************************************************************************/
#ifndef included_tbox_AsyncCommGroup
#define included_tbox_AsyncCommGroup

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Supports all-to-one and one-to-all asynchronous
 * communication operations within a given group of processes by
 * sending messages along the branches of a conceptual tree.
 *
 * This class was created to perform certain group communications
 * without using MPI global communications, which require creating new
 * MPI communicators (can be expensive) and does not support
 * asynchronous operations (until MPI-2).
 *
 * The supported communications are asynchronous in that you can start
 * one and wait for it or check back on it occassionally until it
 * completes.  Simultaneous Asynchronous operations of many groups can
 * be done by using a AsyncCommStage to allocate the groups and to
 * check for completed communications.
 *
 * Supported operations are currently broadcast, gather and sum
 * reduce.  Only integer data is supported.
 *
 * A tree is an acyclic graph in which a node at position pos has
 * nchild children, and the following positions for its
 *
 * - parent: (pos-1)/nchild
 * - first (oldest) child: pos*nchild+1
 * - last (youngest) child: (pos+1)*nchild
 *
 * For example, nchild=2 corresponds to a binary tree.
 *
 * Communication is done by sending messages toward the root (for
 * all-to-one operations) or leaves (for one-to-all operations).  For
 * the former, we receive data from the children and send to the
 * parent.  For the latter, we receive from the parent and send to the
 * children.  Thus every communication involves a receive and a send
 * (except at the root and leaf nodes of the tree).
 *
 * Using a tree generally gives better performance than having all
 * processes in the the tree communicate directly with the root
 * process.  Using MPI communicators corresponding to the groups may
 * faster than using this class, but the cost of creating MPI
 * communicators MAY be expensive.
 *
 * This class supports communication and uses MPI for message passing.
 * If MPI is disabled, the job of this class disappears and the class
 * is effectively empty.  The public interfaces still remain so the
 * class can compile, but the implementations are trivial.
 */
class AsyncCommGroup:public AsyncCommStage::Member
{

private:
   //! @brief Operations user would want to do.
   enum BaseOp { undefined,
                 gather,
                 bcast,
                 max_reduce,
                 min_reduce,
                 sum_reduce };
   //! @brief Tasks, executed in order, to complete a base operation.
   enum TaskOp { recv_start,
                 recv_check,
                 send_start,
                 send_check,
                 none };

public:
   /*!
    * @brief Default constructor does not set up anything.
    * You must initialize() the object before using it.
    */
   AsyncCommGroup();

   /*!
    * @brief Construct communication group.
    *
    * The number of children per node is flexible.
    *
    * @param nchild Number of children per tree node in the group,
    *        i.e., nchild=2 is a binary tree.
    * @param stage
    * @param handler
    *
    * @post nchild == numberOfRequests()
    */
   AsyncCommGroup(
      const size_t nchild,
      AsyncCommStage* stage,
      AsyncCommStage::Handler* handler = 0);

   /*!
    * @brief Destructor.
    *
    * @pre isDone()
    */
   virtual ~AsyncCommGroup();

   /*!
    * @brief Initialize the object.
    *
    * Attach self to the given stage and set the Handler.
    *
    * @param nchild Number of children per tree node in the group,
    * i.e., nchild=2 is a binary tree.
    *
    * @param stage The required stage used for completing non-blocking
    * message passing calls.
    *
    * @param handler Optional handler (see AsyncCommStage::Member).
    *
    * @pre isDone()
    */
   void
   initialize(
      const int nchild,
      AsyncCommStage* stage,
      AsyncCommStage::Handler* handler = 0);

   //@{
   //! @name Define the communication group

   /*!
    * @brief Setup the tree for the given group of processes.  The
    * root process is specified by its index in the group array.
    *
    * The root rank is specified by dereferencing @c group array with
    * @c root_index.
    *
    * @pre getNextTaskOp() == none
    */
   void
   setGroupAndRootIndex(
      const SAMRAI_MPI& mpi,
      const int* group_ranks,
      const int group_size,
      const int root_index);

   /*!
    * @brief Setup the group for the given group.  The root process is
    * specified by its rank.
    *
    * The rank of the root is root_rank, which must be one of the
    * ranks given in the group.
    */
   void
   setGroupAndRootRank(
      const SAMRAI_MPI& mpi,
      const int* group_ranks,
      const int group_size,
      const int root_rank);

   //@}

   /*!
    * @brief Set the MPI tag used for communication within the group.
    *
    * @attention This class is NOT (and cannot be) responsible for
    * ensuring that the MPI communicator and tag are sufficient to
    * select the correct messages.  Please specify appropriate values
    * for the MPI communicator and tag.  Very elusive bugs can occur
    * if incorrect messages are received.  To be safe, it is best to
    * create a new communicator to avoid interference with other
    * communications within SAMRAI.
    *
    * @pre isDone()
    */
   void
   setMPITag(
      const int mpi_tag);

   /*!
    * @brief Returns the MPI tag used for communication within the group.
    */
   int
   getMPITag() const
   {
      return d_mpi_tag;
   }

   /*!
    * @brief Returns the size of the group.
    */
   int
   getGroupSize() const
   {
      return d_group_size;
   }

   /*!
    * @brief Returns next task in a current communication operation.
    */
   TaskOp
   getNextTaskOp() const
   {
      return d_next_task_op;
   }

   /*!
    * @brief Returns operation being performed.
    */
   BaseOp
   getBaseOp() const
   {
      return d_base_op;
   }

   /*!
    * @brief Rank of parent process in the group.
    */
   int
   getParentRank() const
   {
      return d_parent_rank;
   }

   /*!
    * @brief Set whether to use native MPI collective function calls
    * when group includes all ranks in the MPI communicator.
    *
    * This option is off by default to avoid MPI lock-ups.  If you use
    * it, make sure all processors can get to the collective operation
    * to avoid lock-ups.
    *
    * @pre isDone()
    */
   void
   setUseMPICollectiveForFullGroups(
      bool use_mpi_collective = true);

   /*!
    * @brief Set whether sends to parents should be blocking.
    *
    * The default is to use blocking send to parent.  Because there is
    * just one parent, short messages can be buffered by MPI to
    * improve the performance of blocking sends.  Blocking sends need
    * not be checked for completion.
    */
   void
   setUseBlockingSendToParent(
      const bool flag)
   {
      d_use_blocking_send_to_parent = flag;
   }

   /*!
    * @brief Set whether sends to children should be blocking.
    *
    * The default is to use nonblocking send to children.
    * Nonblocking sends to children are generally appropriate
    * as there are multiple children.
    */
   void
   setUseBlockingSendToChildren(
      const bool flag)
   {
      d_use_blocking_send_to_children = flag;
   }

   //@{

   /*!
    * @name Communication methods
    */

   /*!
    * @brief Begin a broadcast communication.
    *
    * Root process of broadcast may send less data (smaller size) than
    * receivers of broadcast, in which case the missing data is
    * considered irrelevant by the root.
    *
    * If this method returns false, checkBcast() must be called until
    * it returns true before any change in object state is allowed.
    *
    * @return Whether operation is completed.
    *
    * @pre getNextTaskOp() == none
    */
   bool
   beginBcast(
      int* buffer,
      int size);

   /*!
    * @brief Check the current broadcast communication and complete
    * the broadcast if all MPI requests are fulfilled.
    *
    *
    * If no communication is in progress, this call does nothing.
    *
    * @return Whether operation is completed.
    *
    * @pre getBaseOp() == bcast
    *
    * @post (getParentRank() != -1) || (getNextTaskOp() != recv_check)
    */
   bool
   checkBcast();

   /*!
    * @brief Begin a gather communication.
    *
    * The gather operation mimics the results of MPI_Gather.
    *
    * Sending processes of gather may send less data (smaller size)
    * than receivers, in which case the missing data is considered
    * irrelevant by the sender.
    *
    * If this method returns false, checkGather() must be called until
    * it returns true before any change in object state is allowed.
    *
    * On non-root processes, buffer should contain the data to be
    * gathered.  On the root process, it should have enough space for
    * all the data from all the processes in the group.
    *
    * @param buffer Data to gather.
    *
    * @param size Number of items contributed by each process.  This
    * must be the same across processes.  However, the root's actual
    * usable buffer must be big enough to hold the gathered data (size
    * times the number of processes in the group).
    *
    * @return Whether operation is completed.
    *
    * @pre getNextTaskOp() == none
    */
   bool
   beginGather(
      int* buffer,
      int size);

   /*!
    * @brief Check the current gather communication and complete the
    * gather if all MPI requests are fulfilled.
    *
    * @return Whether operation is completed.
    *
    * @pre getBaseOp() == gather
    */
   bool
   checkGather();

   /*!
    * @brief Begin a sum reduce communication.
    *
    * Assume all messages are the same size.
    *
    * If this method returns false, checkSumReduce() must be called
    * until it returns true before any change in object state is
    * allowed.
    *
    * Buffer should contain the data to be gathered.
    *
    * @return Whether operation is completed.
    *
    * @pre getNextTaskOp() == none
    */
   bool
   beginSumReduce(
      int* buffer,
      int size);

   /*!
    * @brief Check the current sum reduce communication and complete
    * the sum reduce if all MPI requests are fulfilled.
    *
    * @return Whether operation is completed.
    */
   bool
   checkSumReduce()
   {
      return checkReduce();
   }

   /*!
    * @brief Check the current communication and complete it if all
    * MPI requests are fulfilled.
    */
   bool
   proceedToNextWait();

   /*!
    * @brief Whether the last communication operation has finished.
    *
    * This means more than just whether there is pending MPI requests
    * such as that returned by hasPendingRequests().  The communication
    * may be more complex, requiring several messages and copying of the
    * received message into the correct buffer.
    *
    * @pre (getNextTaskOp() != none) || !hasPendingRequests()
    */
   bool
   isDone() const;

   /*!
    * @brief Wait for the current operation to complete.
    */
   void
   completeCurrentOperation();

   //@}

   int
   getNumberOfChildren() const
   {
      return static_cast<int>(d_nchild);
   }

   void
   logCurrentState(
      std::ostream& co) const;

private:
   /*
    * @brief Assert that user-set MPI parameters are valid.
    *
    * @pre getMPITag() >= 0
    */
   void
   checkMPIParams();

   /*!
    * @brief Operation disallowed due to primitive internal memory management.
    */
   AsyncCommGroup(
      const AsyncCommGroup& r);

   /*!
    * @brief Operation disallowed by primitive internal memory management.
    */
   AsyncCommGroup&
   operator = (
      const AsyncCommGroup& r);

   /*
    * Use MPI collective function call to do communication.
    */
   bool
   bcastByMpiCollective()
   {
      d_mpi.Bcast(d_external_buf, d_external_size, MPI_INT, d_root_rank);
      d_next_task_op = none;
      return true;
   }

   /*
    * Use MPI collective function call to do communication.
    */
   bool
   gatherByMpiCollective();

   /*
    * Use MPI collective function call to do communication.
    */
   bool
   reduceByMpiCollective();

   struct ChildData {
      //! @brief Rank of child process in the group.
      int rank;
      //! @brief Number of descendants on each child branch.
      int size;
      ChildData();
   };

   /*!
    * @brief Begin a generic reduce communication.
    *
    * This method is the workhorse underneath the public reduce methods.
    *
    * If this method returns false, proceedToNextWait() must
    * be called until it returns true before any change
    * in object state is allowed.
    *
    * Buffer should contain the data to be gathered.
    *
    * @return Whether operation is completed.
    */
   bool
   beginReduce();

   /*!
    * @brief Check the current gather communication.
    *
    * This method is the workhorse underneath the public reduce methods.
    *
    * @return Whether operation is completed.
    *
    * @pre (getBaseOp() == max_reduce) || (getBaseOp() == min_reduce) ||
    *      (getBaseOp() == sum_reduce)
    *
    * @post (getParentRank() != -1) || (getNextTaskOp() != send_check)
    * @post (getNextTaskOp() == none) || (numberOfPendingRequests() > 0)
    */
   bool
   checkReduce();

   /*!
    * @brief Perform reduction on data that after it has been brought
    * to the local process.
    *
    * The exact reduce operation depends on the base operation.
    */
   void
   reduceData(
      int* output,
      const int* data) const;

   /*!
    * @brief Compute the data that depends on the group definition.
    *
    * Extract and compute parts and characteristics of the tree
    * relevant to the local process.
    */
   void
   computeDependentData(
      const int* group_ranks,
      const int group_size);

   void
   resetStatus()
   {
      d_mpi_status.MPI_TAG =
         d_mpi_status.MPI_SOURCE =
            d_mpi_status.MPI_ERROR = -1;
   }

   //@{
   /*!
    * @name Mappings between array indices, group positions and process ranks
    */

   /*
    * pos refers the process position in the group (where root position == 0)
    * idx refers the index of the process in the group
    * rank refers the process rank
    */

   /*!
    * @brief Convert the array index to the position.
    *
    * @pre (index >= 0) && (index < getGroupSize())
    */
   int
   toPosition(
      int index) const;
   /*!
    * @brief Convert the position to the array index.
    *
    * @pre (position >= 0) && (position < getGroupSize())
    */
   int
   toIndex(
      int position) const;

   /*!
    * @brief Compute the position of child child_id of a parent (whether
    * or not that child really exists).
    *
    * @param parent_pos Position of the parent in the group.
    * @param ic Index of the child.  (Zero coresponds to the first child.)
    *
    * @pre (parent_pos >= 0) && (parent_pos < getGroupSize())
    */
   int
   toChildPosition(
      int parent_pos,
      int ic) const;

   /*!
    * @brief Compute the oldest (lowest position) child position
    * of a given position (whether or not that child really exists).
    *
    * Same as toChildPosition( parent_pos, 0 );
    *
    * @pre (parent_pos >= 0) && (parent_pos < getGroupSize())
    */
   int
   toOldest(
      int parent_pos) const;
   /*!
    * @brief Compute the youngest (highest position) child position
    * of a given position (whether or not that child really exists).
    *
    * Same as toChildPosition( parent_pos, d_nchild-1 );
    *
    * @pre (parent_pos >= 0) && (parent_pos < getGroupSize())
    */
   int
   toYoungest(
      int parent_pos) const;

   //@}

   /*!
    * @brief Initialize static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      t_reduce_data = TimerManager::getManager()->
         getTimer("tbox::AsyncCommGroup::reduceData()");
      t_wait_all = TimerManager::getManager()->
         getTimer("tbox::AsyncCommGroup::mpi_wait_all");
   }

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
      t_reduce_data.reset();
      t_wait_all.reset();
   }

   /*!
    * @brief Number of children per node.
    */
   size_t d_nchild;

   /*!
    * @brief Index of the local process in d_group_ranks.
    *
    * The group is defined by an array of process ranks.
    * The index of a process refers to the index in the this array.
    * The "position" of a process in the group represents
    * the position in the group, where the nodes
    * are numbered sequentially, starting at zero for the root.
    *
    * We require that the root has position zero.
    * If the index of the root is zero, then positions are
    * identical to indices.  If not, then the index of the
    * root is swapped with zero to get positions.  These two
    * cases correspond respectively to following trivial and
    * nontrivial maps.
    *
    * @verbatim
    *
    *                 Trivial map           Nontrivial map
    * Parameter     d_root_idx == 0         d_root_idx > 0
    * ---------     ---------------         --------------
    *
    * index of root       0              d_root_idx
    *
    * index of          d_idx            d_idx
    * local process
    *
    * index of            p              p == 0 ? d_root_idx :
    * position p                         p == d_root_idx ? 0 :
    *                                    p
    *
    * position of         i              i == 0 ? d_root_idx :
    * index i                            i == d_root_idx ? 0 :
    *                                    i
    *
    * @endverbatim
    */
   int d_idx;

   /*!
    * @brief Index of the root process in d_group_ranks.
    */
   int d_root_idx;

   /*!
    * @brief Index of the root process in d_group_ranks.
    */
   int d_root_rank;

   int d_group_size;

   /*!
    * @brief Rank of parent process in the group.
    *
    * If negative, there is no parent (this is the root).
    * In send_start tasks, send only to children with valid ranks
    * (not -1).
    */
   int d_parent_rank;

   //! @brief Data on each child branch.
   ChildData* d_child_data;

   /*!
    * @brief Total of all branch sizes.
    */
   int d_branch_size_totl;

   /*!
    * @brief Operation being performed.
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
    * @brief External data input and output buffer.
    *
    * This provides the input and output for transfering data.
    * The expected size of the buffer depends on the communication.
    */
   int* d_external_buf;

   /*!
    * @brief Size of d_external_buf.
    */
   int d_external_size;

   /*!
    * @brief Internal buffer scratch space.
    *
    * Used for gather and reduce operations but not for
    * broadcast.  Not used when using MPI collective calls.
    */
   std::vector<int> d_internal_buf;

   int d_mpi_tag;
   SAMRAI_MPI d_mpi;

   bool d_use_mpi_collective_for_full_groups;
   bool d_use_blocking_send_to_children;
   bool d_use_blocking_send_to_parent;

   // Make some temporary variable statuses to avoid repetitious allocations.
   SAMRAI_MPI::Status d_mpi_status;

   int d_mpi_err;

   static std::shared_ptr<Timer> t_reduce_data;
   static std::shared_ptr<Timer> t_wait_all;

   static StartupShutdownManager::Handler s_initialize_finalize_handler;

#ifdef DEBUG_CHECK_ASSERTIONS
   /*!
    * @brief Array of process ranks in the group.
    *
    * It is possible to code this class without storing all the
    * ranks internally.  However, it is easier to debug if the
    * ranks are available.
    */
   std::vector<int> d_group_ranks;
#endif

};

}
}

#endif  // included_tbox_AsyncCommGroup
