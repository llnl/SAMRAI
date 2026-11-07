/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   All-to-one and one-to-all communication using a tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/AsyncCommGroup.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include STL_SSTREAM_HEADER_FILE

#ifdef OSTRINGSTREAM_TYPE_IS_BROKEN
#ifdef OSTRSTREAM_TYPE_IS_BROKEN
#error "Neither std::ostringstream nor std::ostrstream works"
#else
typedef std::ostringstream std::ostrstream
#endif
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

/*
 * This class uses a non-deterministic algorithm, which can be
 * very hard to debug.  To help debugging, we keep some special
 * debugging code that is activated when AsyncCommGroup_DEBUG_OUTPUT
 * is defined.
 */
// #define AsyncCommGroup_DEBUG_OUTPUT

std::shared_ptr<Timer> AsyncCommGroup::t_reduce_data;
std::shared_ptr<Timer> AsyncCommGroup::t_wait_all;

StartupShutdownManager::Handler
AsyncCommGroup::s_initialize_finalize_handler(
   AsyncCommGroup::initializeCallback,
   0,
   0,
   AsyncCommGroup::finalizeCallback,
   StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 * Construct a simple communication group that does not work
 * with a communication stage.  All parameters are set to reasonable
 * defaults or, if appropriate, invalid values.
 ***********************************************************************
 */
AsyncCommGroup::AsyncCommGroup():
   AsyncCommStage::Member(),
   d_nchild(MathUtilities<size_t>::getMax()),
   d_idx(-1),
   d_root_idx(-1),
   d_parent_rank(-1),
   d_child_data(0),
   d_branch_size_totl(-1),
   d_base_op(undefined),
   d_next_task_op(none),
   d_external_buf(0),
   d_external_size(0),
   d_internal_buf(),
   d_mpi_tag(-1),
   d_mpi(SAMRAI_MPI::getSAMRAIWorld()),
   d_use_mpi_collective_for_full_groups(false),
   d_use_blocking_send_to_children(false),
   d_use_blocking_send_to_parent(false)
#ifdef DEBUG_CHECK_ASSERTIONS
   ,
   d_group_ranks(0, true)
#endif
{
}

/*
 ***********************************************************************
 * Construct a simple communication group that does not work
 * with a communication stage.  All parameters are set to reasonable
 * defaults or, if appropriate, invalid values.
 ***********************************************************************
 */
AsyncCommGroup::AsyncCommGroup(
   const size_t nchild,
   AsyncCommStage* stage,
   AsyncCommStage::Handler* handler):
   AsyncCommStage::Member(nchild, stage, handler),
   d_nchild(nchild),
   d_idx(-1),
   d_root_idx(-1),
   d_parent_rank(-1),
   d_child_data(new ChildData[nchild]),
   d_branch_size_totl(-1),
   d_base_op(undefined),
   d_next_task_op(none),
   d_external_buf(0),
   d_external_size(0),
   d_internal_buf(),
   d_mpi_tag(-1),
   d_mpi(SAMRAI_MPI::getSAMRAIWorld()),
   d_use_mpi_collective_for_full_groups(false),
   d_use_blocking_send_to_children(false),
   d_use_blocking_send_to_parent(false)
#ifdef DEBUG_CHECK_ASSERTIONS
   ,
   d_group_ranks(0, true)
#endif
{
   TBOX_ASSERT(nchild == numberOfRequests());
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommGroup::AsyncCommGroup(
   const AsyncCommGroup& r):
   AsyncCommStage::Member(0, 0, 0),
   d_nchild(0),
   d_mpi(SAMRAI_MPI::getSAMRAIWorld())
{
   NULL_USE(r);
   TBOX_ERROR(
      "Copy constructor disallowed due to primitive internal memory management.");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommGroup&
AsyncCommGroup::operator = (
   const AsyncCommGroup& r) {
   NULL_USE(r);
   TBOX_ERROR(
      "Assignment operator disallowed due to primitive internal memory management.");
   return *this;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommGroup::~AsyncCommGroup()
{
   if (!isDone()) {
      TBOX_ERROR("Deallocating a group while communication is pending\n"
         << "leads to lost messages.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
   }
   delete[] d_child_data;
   d_child_data = 0;
}

/*
 ***********************************************************************
 * Initialize data as if constructed with the given arguments.
 ***********************************************************************
 */
void
AsyncCommGroup::initialize(
   const int nchild,
   AsyncCommStage* stage,
   AsyncCommStage::Handler* handler)
{
   if (!isDone()) {
      TBOX_ERROR("It is illegal to re-initialize a AsyncCommGroup\n"
         << "while it has current messages.\n");
   }
   attachStage(nchild, stage);
   setHandler(handler);
   d_nchild = nchild;
   delete[] d_child_data;
   d_idx = -1;
   d_root_idx = -1;
   d_parent_rank = -1;
   d_child_data = new ChildData[nchild];
   d_branch_size_totl = -1;
   d_base_op = undefined;
   d_next_task_op = none;
#ifdef DEBUG_CHECK_ASSERTIONS
   d_group_ranks.clear();
#endif
}

/*
 *********************************************************************
 * Check whether the current (or last) operation has completed.
 *********************************************************************
 */
bool
AsyncCommGroup::proceedToNextWait()
{
   switch (d_base_op) {
      case gather: return checkGather();

      case bcast: return checkBcast();

      case min_reduce:
      case max_reduce:
      case sum_reduce: return checkReduce();

      case undefined:
         TBOX_ERROR("There is no current operation to check.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
         break;
      default:
         TBOX_ERROR("Library error: attempt to use an operation that\n"
         << "has not been written yet"
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
   }
   return true;
}

/*
 ****************************************************************
 * Whether the last communicatoin operation has finished.
 ****************************************************************
 */
bool
AsyncCommGroup::isDone() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (getNextTaskOp() == none) {
      TBOX_ASSERT(!hasPendingRequests());    // Verify sane state.
   }
#endif
   return d_next_task_op == none;
}

/*
 *********************************************************************
 * Wait for current communication operation to complete.
 *
 * Wait for all requests to come in and call proceedToNextWait()
 * until all tasks of the communication operation is complete.
 *********************************************************************
 */
void
AsyncCommGroup::completeCurrentOperation()
{
   SAMRAI_MPI::Request * const req = getRequestPointer();
   SAMRAI_MPI::Status* mpi_stat = d_next_task_op == none ?
      0 : new SAMRAI_MPI::Status[d_nchild];

   while (d_next_task_op != none) {

      t_wait_all->start();
      int errf = SAMRAI_MPI::Waitall(static_cast<int>(d_nchild),
            req,
            mpi_stat);
      t_wait_all->stop();

      if (errf != MPI_SUCCESS) {
         TBOX_ERROR("Error in MPI_Waitall call.\n"
            << "mpi_communicator = " << d_mpi.getCommunicator()
            << "mpi_tag = " << d_mpi_tag);
      }

      proceedToNextWait();

   }

   if (mpi_stat != 0) {
      delete[] mpi_stat;
   }
}

/*
 ************************************************************************
 * Set internal parameters for performing the broadcast
 * and call checkBcast to perform the communication.
 ************************************************************************
 */
bool
AsyncCommGroup::beginBcast(
   int* buffer,
   int size)
{
   if (getNextTaskOp() != none) {
      TBOX_ERROR("Cannot begin communication while another is in progress."
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif
   d_external_buf = buffer;
   d_external_size = size;
   d_base_op = bcast;

   if (d_use_mpi_collective_for_full_groups && d_group_size == d_mpi.getSize()) {
      return bcastByMpiCollective();
   }

   d_next_task_op = recv_start;
   return checkBcast();
}

/*
 ************************************************************************
 * Broadcast is an one-to-all operation, so we receive from the
 * parent process and send to the children processes.
 ************************************************************************
 */
bool
AsyncCommGroup::checkBcast()
{
   if (getBaseOp() != bcast) {
      TBOX_ERROR("Cannot check nonexistent broadcast operation."
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
   }
   SAMRAI_MPI::Request * const req = getRequestPointer();
   size_t ic;
   int flag = 0;

   if (d_next_task_op != none) {
      bool task_entered = false;
      if (d_next_task_op == recv_start) {
         task_entered = true;
         if (d_parent_rank > -1) {
            d_mpi_err = d_mpi.Irecv(d_external_buf,
                  d_external_size,
                  MPI_INT,
                  d_parent_rank,
                  d_mpi_tag,
                  &req[0]);
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in MPI_Irecv."
                  << "mpi_communicator = " << d_mpi.getCommunicator()
                  << "mpi_tag = " << d_mpi_tag);
            }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            plog << "tag-" << d_mpi_tag
                 << " expecting " << d_external_size
                 << " from " << d_parent_rank
                 << " in checkBcast"
                 << std::endl;
#endif
         }
      }
      bool breakout = false;
      if (d_next_task_op == recv_start || d_next_task_op == recv_check) {
         task_entered = true;
         if (req[0] != MPI_REQUEST_NULL) {
            resetStatus();
            d_mpi_err = SAMRAI_MPI::Test(&req[0], &flag, &d_mpi_status);
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in MPI_Test.\n"
                  << "Error-in-status is "
                  << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                  << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                  << '\n'
                  << "mpi_communicator = " << d_mpi.getCommunicator()
                  << "mpi_tag = " << d_mpi_tag);
            }
            TBOX_ASSERT((req[0] == MPI_REQUEST_NULL) == (flag == 1));
            if (flag == 1) {
#ifdef DEBUG_CHECK_ASSERTIONS
               int count = -1;
               d_mpi_err = SAMRAI_MPI::Get_count(&d_mpi_status, MPI_INT, &count);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Get_count.\n"
                     << "Error-in-status is "
                     << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                     << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                     << '\n'
                     << "mpi_communicator = " << d_mpi.getCommunicator()
                     << "mpi_tag = " << d_mpi_tag);
               }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               plog << "tag-" << d_mpi_tag
                    << " received " << count
                    << " from " << d_mpi_status.MPI_SOURCE
                    << " in checkBcast"
                    << std::endl;
#endif
               TBOX_ASSERT(count <= d_external_size);
               TBOX_ASSERT(d_mpi_status.MPI_TAG == d_mpi_tag);
               TBOX_ASSERT(d_mpi_status.MPI_SOURCE == d_parent_rank);
#endif
            } else {
               d_next_task_op = recv_check;
               breakout = true;
            }
         }
      }

      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start)) {
         task_entered = true;
         for (ic = 0; ic < d_nchild; ++ic) {
            if (d_child_data[ic].rank >= 0) {
               if (d_use_blocking_send_to_children) {
                  d_mpi_err = d_mpi.Send(d_external_buf,
                        d_external_size,
                        MPI_INT,
                        d_child_data[ic].rank,
                        d_mpi_tag);
               } else {
                  d_mpi_err = d_mpi.Isend(d_external_buf,
                        d_external_size,
                        MPI_INT,
                        d_child_data[ic].rank,
                        d_mpi_tag,
                        &req[ic]);
               }
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in send."
                     << "mpi_communicator = " << d_mpi.getCommunicator()
                     << "mpi_tag = " << d_mpi_tag);
               }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               plog << "tag-" << d_mpi_tag
                    << " sending " << d_external_size
                    << " to " << d_child_data[ic].rank
                    << " in checkBcast"
                    << std::endl;
#endif
            }
         }
      }

      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start || d_next_task_op == send_check)) {
         task_entered = true;
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               resetStatus();
               d_mpi_err = SAMRAI_MPI::Test(&req[ic], &flag, &d_mpi_status);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Test.\n"
                     << "Error-in-status is "
                     << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                     << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                     << '\n'
                     << "mpi_communicator = " << d_mpi.getCommunicator()
                     << "mpi_tag = " << d_mpi_tag);
               }
               TBOX_ASSERT((req[ic] == MPI_REQUEST_NULL) == (flag == 1));
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               if (req[ic] == MPI_REQUEST_NULL) {
                  plog << "tag-" << d_mpi_tag
                       << " sent unknown size (MPI convention)"
                       << " to " << d_child_data[ic].rank
                       << " in checkBcast"
                       << std::endl;
               }
#endif
            }
         }
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               break;
            }
         }
         if (!breakout && ic < d_nchild) {
            d_next_task_op = send_check;
            breakout = true;
         }
         if (!breakout) { 
            d_next_task_op = none;
         }
      }

      if (!task_entered) {
         TBOX_ERROR("checkBcast is incompatible with current state."
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
      }
   }

   if (d_parent_rank == -1) {
      TBOX_ASSERT(getNextTaskOp() != recv_check);
   }

   return d_next_task_op == none;
}

/*
 ************************************************************************
 ************************************************************************
 */
bool
AsyncCommGroup::gatherByMpiCollective()
{
   if (d_mpi.getSize() > 1) {
      d_internal_buf.clear();
      d_internal_buf.insert(d_internal_buf.begin(),
         d_external_buf,
         d_external_buf + d_external_size);
      d_mpi.Gather(&d_internal_buf[0],
         d_external_size,
         MPI_INT,
         d_external_buf,
         d_external_size,
         MPI_INT,
         d_root_rank);
      d_internal_buf.clear();
   }
   d_next_task_op = none;
   return true;
}

/*
 ************************************************************************
 * Allocate enough memory internally to store all descendent data.
 * Place local process's contribution in the internal buffer.
 * Call checkGather to obtain data from descendants and send
 * to parent.
 ************************************************************************
 */
bool
AsyncCommGroup::beginGather(
   int* buffer,
   int size)
{

   if (getNextTaskOp() != none) {
      TBOX_ERROR("Cannot begin communication while another is in progress."
         << "mpi_communicator = " << d_mpi.getCommunicator()
         << "mpi_tag = " << d_mpi_tag);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif

   d_base_op = gather;
   d_external_buf = buffer;
   d_external_size = size;

   if (d_use_mpi_collective_for_full_groups && d_group_size == d_mpi.getSize()) {
      return gatherByMpiCollective();
   }

   /*
    * The internal buffer size is d_branch_size_totl+1 times the
    * message size.  There is one data block for each descendent,
    * plus one for the local contribution:
    *
    * |<------------------------ internal buffer ---------------------->|
    * |                                                                 |
    * |<--msg_size-->|<--msg_size-->| ... |<--msg_size-->|<--msg_size-->|
    * |              |              | ... |              |              |
    * |  recv from   |  recv from   | ... |  recv from   |    local     |
    * | descendant 0 | descendant 1 | ... | dsndt (nb-1) | contribution |
    *
    * (nb = d_branch_size_totl)
    *
    * The message size is the data buffer size, plus one integer
    * describing the index of the process contributing the data.
    */

   /*
    * Allocate enough space for data from this position plus all
    * descendant positions.  Each position contributes its external
    * data plus some data to help sort the final gathered data.
    */
   int per_proc_msg_size = (1 + d_external_size);
   d_internal_buf.clear();
   d_internal_buf.insert(d_internal_buf.end(),
      (d_branch_size_totl + 1) * per_proc_msg_size,
      0);

   /*
    * Add our contribution to the gathered data.
    */
   int* ptr = &d_internal_buf[0]
      + (d_branch_size_totl) * per_proc_msg_size;
   *(ptr++) = d_idx;
   int i;
   for (i = 0; i < d_external_size; ++i) {
      ptr[i] = d_external_buf[i];
   }

   d_next_task_op = recv_start;

   return checkGather();
}

/*
 ************************************************************************
 * Gather is an all-to-one operation, so we receive from the
 * children processes and send to the parent process.
 ************************************************************************
 */
bool
AsyncCommGroup::checkGather()
{
   if (getBaseOp() != gather) {
      TBOX_ERROR("Cannot check nonexistent gather operation\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }

   SAMRAI_MPI::Request * const req = getRequestPointer();
   int per_proc_msg_size = (1 + d_external_size);

   int i;
   size_t ic;
   int older_sibling_size;
   int flag = 0;

   if (d_next_task_op != none) {
      bool task_entered = false;
      if (d_next_task_op == recv_start) {
         task_entered = true;
         older_sibling_size = 0;
         for (ic = 0; ic < d_nchild; ++ic) {
            if (d_child_data[ic].rank >= 0) {
               /*
                * Child number ic exists.  We'll put its data after
                * its older siblings' data.
                */
               d_mpi_err = d_mpi.Irecv(&d_internal_buf[0]
                     + per_proc_msg_size * older_sibling_size,
                     d_child_data[ic].size * per_proc_msg_size,
                     MPI_INT,
                     d_child_data[ic].rank,
                     d_mpi_tag,
                     &req[ic]);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Irecv.\n"
                     << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
               }
               older_sibling_size += d_child_data[ic].size;
            }
         }
      }

      bool breakout = false;
      if (d_next_task_op == recv_start || d_next_task_op == recv_check) {
         task_entered = true;
         /*
          * Check all pending receives from the children.
          */
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               resetStatus();
               d_mpi_err = SAMRAI_MPI::Test(&req[ic], &flag, &d_mpi_status);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Test.\n"
                     << "Error-in-status is "
                     << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                     << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                     << '\n'
                     << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
               }
               TBOX_ASSERT((req[ic] == MPI_REQUEST_NULL) == (flag == 1));
#ifdef DEBUG_CHECK_ASSERTIONS
               if (flag == 1) {
                  TBOX_ASSERT(d_mpi_status.MPI_TAG == d_mpi_tag);
                  TBOX_ASSERT(d_mpi_status.MPI_SOURCE == d_child_data[ic].rank);
                  int count = -1;
                  SAMRAI_MPI::Get_count(&d_mpi_status, MPI_INT, &count);
#ifdef AsyncCommGroup_DEBUG_OUTPUT
                  plog << "tag-" << d_mpi_tag
                       << " received " << count
                       << " from " << d_mpi_status.MPI_SOURCE
                       << " in checkGather"
                       << std::endl;
#endif
                  if (count > d_child_data[ic].size * per_proc_msg_size) {
                     TBOX_ERROR("Message size bigger than expected from proc "
                        << d_child_data[ic].rank << "\n"
                        << "Expect "
                        << d_child_data[ic].size * per_proc_msg_size
                        << "\n"
                        << "Actual " << count << '\n'
                        << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                        << "mpi_tag = " << d_mpi_tag << '\n');
                  }
               }
#endif
            }
         }

         /*
          * If there are still pending requests, we cannot complete
          * the communication operation at this time.
          */
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               break;
            }
         }
         if (ic < d_nchild) {
            d_next_task_op = recv_check;
            breakout = true;
         }

         /*
          * At this point, all receives are completed.
          */
         if (!breakout && d_parent_rank < 0) {
            /*
             * The root of the gather (only the root!) transfers the
             * internal buffer into the external buffer, unshuffling
             * data in the process.
             */
            int n;
            for (n = 0; n < d_group_size; ++n) {
               int* ptr = &d_internal_buf[0]
                  + per_proc_msg_size * n;
               const int source_idx = *(ptr++);
               if (source_idx < 0 && source_idx >= d_group_size) {
                  TBOX_ERROR("Gathered data has out of range index "
                     << source_idx << '\n'
                     << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
               }
               int* dest_buf = d_external_buf
                  + d_external_size * (source_idx);
               for (i = 0; i < d_external_size; ++i, ++ptr) {
                  dest_buf[i] = *ptr;
               }
            }
         }
      }

      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start)) {
         task_entered = true;
         if (d_parent_rank >= 0) {
            if (d_use_blocking_send_to_parent) {
               d_mpi_err = d_mpi.Send(&d_internal_buf[0],
                     static_cast<int>(d_internal_buf.size()),
                     MPI_INT,
                     d_parent_rank,
                     d_mpi_tag);
            } else {
               d_mpi_err = d_mpi.Isend(&d_internal_buf[0],
                     static_cast<int>(d_internal_buf.size()),
                     MPI_INT,
                     d_parent_rank,
                     d_mpi_tag,
                     &req[0]);
            }
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in send.\n"
                  << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                  << "mpi_tag = " << d_mpi_tag << '\n');
            }
         }
      }

      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start || d_next_task_op == send_check)) {
         task_entered = true;
         if (req[0] != MPI_REQUEST_NULL) {
            resetStatus();
            d_mpi_err = SAMRAI_MPI::Test(&req[0], &flag, &d_mpi_status);
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in MPI_Test.\n"
                  << "Error-in-status is "
                  << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                  << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                  << '\n'
                  << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                  << "mpi_tag = " << d_mpi_tag << '\n');
            }
            TBOX_ASSERT((req[0] == MPI_REQUEST_NULL) == (flag == 1));
         }
         if (req[0] != MPI_REQUEST_NULL) {
            d_next_task_op = send_check;
         } else {
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            if (d_parent_rank > -1) {
               plog << "tag-" << d_mpi_tag
                    << " sent " << d_internal_buf.size()
                    << " to " << d_parent_rank
                    << " in checkGather"
                    << std::endl;
            }
#endif
            d_next_task_op = none;
         }
      }
      if (!task_entered) { 
         TBOX_ERROR("checkGather is incompatible with current state.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
      }
   }

   return d_next_task_op == none;
}

/*
 **********************************************************************
 * Set flag indicating a sum reduce and call the generic reduce method
 * to do the actual work.
 **********************************************************************
 */
bool
AsyncCommGroup::beginSumReduce(
   int* buffer,
   int size)
{
   if (getNextTaskOp() != none) {
      TBOX_ERROR("Cannot begin communication while another is in progress.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   checkMPIParams();
#endif
   d_base_op = sum_reduce;
   d_external_buf = buffer;
   d_external_size = size;
   return beginReduce();
}

/*
 ************************************************************************
 ************************************************************************
 */
bool
AsyncCommGroup::reduceByMpiCollective()
{
   SAMRAI_MPI::Op mpi_op =
      d_base_op == max_reduce ? MPI_MAX :
      d_base_op == min_reduce ? MPI_MIN :
      MPI_SUM;
   d_internal_buf.clear();
   d_internal_buf.insert(d_internal_buf.end(),
      d_external_size,
      0);
   d_mpi.Reduce(d_external_buf,
      &d_internal_buf[0],
      d_external_size,
      MPI_INT,
      mpi_op,
      d_root_rank);
   if (d_parent_rank < 0) {
      copy(d_internal_buf.begin(), d_internal_buf.end(), d_external_buf);
   }
   d_next_task_op = none;
   return true;
}

/*
 ************************************************************************
 * Allocate enough memory internally to store all children data.
 * Place local process's contribution in the internall buffer.
 * Call checkReduce to obtain data from children, reduce the
 * data and send result to parent..
 ************************************************************************
 */
bool
AsyncCommGroup::beginReduce()
{
   if (d_use_mpi_collective_for_full_groups && d_group_size == d_mpi.getSize()) {
      return reduceByMpiCollective();
   }

   int msg_size = d_external_size;
   /*
    * For reducing data, nc = number of actual children.  nc <= d_nchild.
    *
    * The internal buffer stores up to nc+1 times the message size:
    *
    * |<------------------------ internal buffer ---------------------->|
    * |                                                                 |
    * |<--msg_size-->|<--msg_size-->| ... |<--msg_size-->|<--msg_size-->|
    * |              |              | ... |              |              |
    * |  recv from   |  recv from   | ... |  recv from   |   send to    |
    * |   child 0    |   child 1    | ... | child (nc-1) |    parent    |
    * |   (if any)   |   (if any)   | ... |   (if any)   |   (if any)   |
    *
    * If a data block is not needed, it is not created and the
    * following blocks shift over.  For non-root processes,
    * the reduced data is placed in the "send to parent" section
    * so it can be passed up the tree.
    * For the root process, reduced data is placed directly
    * into d_external_buf.
    */

   /*
    * Compute the number of actual children, nc.  Note that
    * d_nchild is just the upper limit to the number of actual
    * children.
    */
   const int oldest_pos = toOldest(toPosition(d_idx));
   int limit_pos = toYoungest(toPosition(d_idx)) + 1;
   limit_pos = MathUtilities<int>::Min(limit_pos, d_group_size);
   const int n_children = limit_pos > oldest_pos ? limit_pos - oldest_pos : 0;

   d_internal_buf.clear();
   d_internal_buf.insert(d_internal_buf.end(),
      msg_size * (n_children + (d_parent_rank > -1)),
      0);

   if (d_parent_rank > -1) {
      int* ptr = &d_internal_buf[0] + d_internal_buf.size() - msg_size;
      for (int i = 0; i < d_external_size; ++i) {
         ptr[i] = d_external_buf[i];
      }
   }

   d_next_task_op = recv_start;

   return checkReduce();
}

/*
 ************************************************************************
 * Reduction is an all-to-one operation, so we receive from
 * the children and send to the parent.  The transfer data
 * toward the root process is the same for all types of
 * reductions.  After each child's data is received, the
 * specific data reduction arithmetic is done locally.
 * The result of the local reduction is sent to the parent.
 ************************************************************************
 */
bool
AsyncCommGroup::checkReduce()
{
   if (!(getBaseOp() == max_reduce ||
         getBaseOp() == min_reduce ||
         getBaseOp() == sum_reduce)) {
      TBOX_ERROR("Cannot check nonexistent reduce operation.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }

   SAMRAI_MPI::Request * const req = getRequestPointer();
   int msg_size = d_external_size;

   size_t ic;
   int flag = 0;
   if (d_next_task_op != none) {
      bool task_entered = false;
      if (d_next_task_op == recv_start) {
         task_entered = true;
         for (ic = 0; ic < d_nchild; ++ic) {
            if (d_child_data[ic].rank >= 0) {
               d_mpi_err = d_mpi.Irecv(&d_internal_buf[0] + ic * msg_size,
                     msg_size,
                     MPI_INT,
                     d_child_data[ic].rank,
                     d_mpi_tag,
                     &req[ic]);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Irecv.\n"
                     << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
               }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
               plog << "tag-" << d_mpi_tag
                    << " expecting " << msg_size
                    << " from " << d_child_data[ic].rank
                    << " in checkReduce"
                    << std::endl;
#endif
            }
         }
      }

      bool breakout = false;
      if (d_next_task_op == recv_start || d_next_task_op == recv_check) {
         task_entered = true;
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               resetStatus();
               d_mpi_err = SAMRAI_MPI::Test(&req[ic], &flag, &d_mpi_status);
               if (d_mpi_err != MPI_SUCCESS) {
                  TBOX_ERROR("Error in MPI_Test.\n"
                     << "Error-in-status is "
                     << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                     << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                     << '\n'
                     << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
               }
               TBOX_ASSERT((req[ic] == MPI_REQUEST_NULL) == (flag == 1));
               if (flag == 1) {
#ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(d_mpi_status.MPI_TAG == d_mpi_tag);
                  TBOX_ASSERT(d_mpi_status.MPI_SOURCE == d_child_data[ic].rank);
                  int count = -1;
                  SAMRAI_MPI::Get_count(&d_mpi_status, MPI_INT, &count);
#ifdef AsyncCommGroup_DEBUG_OUTPUT
                  plog << " child-" << ic << " tag-" << d_mpi_tag
                       << " received " << count
                       << " from " << d_mpi_status.MPI_SOURCE
                       << " in checkReduce"
                       << std::endl;
#endif
                  if (count != msg_size) {
                     TBOX_ERROR(
                        "Did not get the expected message size from proc "
                        << d_child_data[ic].rank << "\n"
                        << "Expect " << msg_size
                        << "\n"
                        << "Actual " << count << '\n'
                        << "mpi_communicator = "
                        << d_mpi.getCommunicator() << '\n'
                        << "mpi_tag = "
                        << d_mpi_tag << '\n');
                  }
#endif
               } else {
#ifdef AsyncCommGroup_DEBUG_OUTPUT
                  plog << " child-" << ic << " tag-" << d_mpi_tag
                       << " still waiting for proc " << d_child_data[ic].rank
                       << " in checkReduce"
                       << std::endl;
                  TBOX_ASSERT(req[ic] != MPI_REQUEST_NULL);
                  TBOX_ASSERT(numberOfPendingRequests() > 0);
#endif
               }
            }
         }
         for (ic = 0; ic < d_nchild; ++ic) {
            if (req[ic] != MPI_REQUEST_NULL) {
               break;
            }
         }
         if (ic < d_nchild) {
            d_next_task_op = recv_check;
            breakout = true;
         }

         if (!breakout) {
            int* local_data = d_parent_rank < 0 ? d_external_buf :
               &d_internal_buf[0] + d_internal_buf.size() - msg_size;
            t_reduce_data->start();
            for (ic = 0; ic < d_nchild; ++ic) {
               if (d_child_data[ic].rank > -1) {
                  int* child_data = &d_internal_buf[0] + ic * msg_size;
                  reduceData(local_data, child_data);
               }
            }
            t_reduce_data->stop();
         }
      }
      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start)) {
         task_entered = true;
         if (d_parent_rank >= 0) {
            int* ptr = &d_internal_buf[0]
               + d_internal_buf.size() - msg_size;
            if (d_use_blocking_send_to_parent) {
               d_mpi_err = d_mpi.Send(ptr,
                     msg_size,
                     MPI_INT,
                     d_parent_rank,
                     d_mpi_tag);
            } else {
               d_mpi_err = d_mpi.Isend(ptr,
                     msg_size,
                     MPI_INT,
                     d_parent_rank,
                     d_mpi_tag,
                     &req[0]);
            }
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in send.\n"
                  << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                  << "mpi_tag = " << d_mpi_tag << '\n');
            }
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            plog << "tag-" << d_mpi_tag
                 << " sending " << msg_size
                 << " to " << d_parent_rank
                 << " in checkReduce"
                 << std::endl;
#endif
         }
      }

      if (!breakout &&
            (d_next_task_op == recv_start || d_next_task_op == recv_check ||
             d_next_task_op == send_start || d_next_task_op == send_check)) {
         task_entered = true;
         if (req[0] != MPI_REQUEST_NULL) {
            resetStatus();
            d_mpi_err = SAMRAI_MPI::Test(&req[0], &flag, &d_mpi_status);
            if (d_mpi_err != MPI_SUCCESS) {
               TBOX_ERROR("Error in MPI_Test.\n"
                  << "Error-in-status is "
                  << (d_mpi_err == MPI_ERR_IN_STATUS) << '\n'
                  << "MPI_ERROR value is " << d_mpi_status.MPI_ERROR
                  << '\n'
                  << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
                  << "mpi_tag = " << d_mpi_tag << '\n');
            }
            TBOX_ASSERT((req[0] == MPI_REQUEST_NULL) == (flag == 1));
         }
         if (req[0] != MPI_REQUEST_NULL) {
#ifdef AsyncCommGroup_DEBUG_OUTPUT
            plog << "tag-" << d_mpi_tag
                 << " sent unknown size (MPI convention)"
                 << " to " << d_parent_rank
                 << " in checkReduce"
                 << std::endl;
#endif
            d_next_task_op = send_check;
         } else {
            d_next_task_op = none;
         }
      }
      if (!task_entered) {
         TBOX_ERROR("checkReduce is incompatible with current state.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
      }
   }

   if (getParentRank() == -1) {
      TBOX_ASSERT(getNextTaskOp() != send_check);
   }

   if (getNextTaskOp() != none) {
      TBOX_ASSERT(numberOfPendingRequests() > 0);
   }

   return d_next_task_op == none;

}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommGroup::setGroupAndRootRank(
   const SAMRAI_MPI& mpi,
   const int* group_ranks,
   const int group_size,
   const int root_rank)
{
   int i;
   for (i = 0; i < group_size; ++i) {
      if (group_ranks[i] == root_rank) {
         break;
      }
   }
   if (i == group_size) {
      TBOX_ERROR(
         "New root " << root_rank << " is not in the group.\n"
                     << "mpi_communicator = " << d_mpi.getCommunicator()
                     << '\n'
                     << "mpi_tag = " << d_mpi_tag << '\n');
   }
   setGroupAndRootIndex(mpi, group_ranks, group_size, i);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommGroup::setGroupAndRootIndex(
   const SAMRAI_MPI& mpi,
   const int* group_ranks,
   const int group_size,
   const int root_index)
{
   if (getNextTaskOp() != none) {
      TBOX_ERROR("AsyncCommGroup::setGroupAndRootIndex:\n"
         << "Changing group while a communication is occuring can\n"
         << "corrupt data and so is disallowed.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }

   d_mpi = mpi;

   // Set the index for local and root processes.
   d_group_size = group_size;
   d_idx = -1;
   for (int i = 0; i < d_group_size; ++i) {
      if (group_ranks[i] == d_mpi.getRank()) {
         d_idx = i;
         break;
      }
   }
   d_root_idx = root_index;
   d_root_rank = group_ranks[d_root_idx];

#ifdef DEBUG_CHECK_ASSERTIONS
   // Set d_group_ranks and do some sanity checks.
   if (d_group_size > d_mpi.getSize()) {
      TBOX_ERROR("AsyncCommGroup::setGroupAndRootIndex:\n"
         << "Group size (" << d_group_size << ") must not be greater than the size of\n"
         << "the MPI communicator group (" << d_mpi.getSize() << ").\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');

   }
   TBOX_ASSERT(d_group_size > 0);
   if (d_idx == -1) {
      TBOX_ERROR(
         "AsyncCommGroup::setGroupAndRootIndex:\n"
         << "The local process (" << d_mpi.getRank()
         << ") MUST be in the communication group.\n"
         << "mpi_communicator = "
         << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }
   d_group_ranks.clear();
   d_group_ranks.insert(d_group_ranks.end(), d_group_size, -1);
   int dup = 0;
   for (int i = 0; i < d_group_size; ++i) {
      if (group_ranks[i] < 0 || group_ranks[i] >= d_mpi.getSize()) {
         TBOX_ERROR(
            "AsyncCommGroup::setGroupAndRootIndex:\n"
            << "Rank " << group_ranks[i] << " is not in the current\n"
            << "MPI communicator.\n"
            << "mpi_communicator = "
            << d_mpi.getCommunicator() << '\n'
            << "mpi_tag = " << d_mpi_tag
            << '\n');
      }
      if (group_ranks[i] == d_mpi.getRank()) ++dup;
      d_group_ranks[i] = group_ranks[i];
   }
   if (dup != 1) {
      TBOX_ERROR("AsyncCommGroup::setGroupAndRootIndex:\n"
         << "The local process must appear exactly once in the group.\n"
         << "It appeared " << dup << " times.\n"
         << "mpi_communicator = " << d_mpi.getCommunicator() << '\n'
         << "mpi_tag = " << d_mpi_tag << '\n');
   }
#endif

   computeDependentData(group_ranks, group_size);
}

/*
 ***************************************************************************
 * Perform reduction on data that after it has been brought to the local
 * process.
 ***************************************************************************
 */
void
AsyncCommGroup::reduceData(
   int* output,
   const int* data) const
{
   int i;
   switch (d_base_op) {
      case max_reduce:
         for (i = 0; i < d_external_size; ++i) {
            if (output[i] < data[i]) output[i] = data[i];
         }
         break;
      case min_reduce:
         for (i = 0; i < d_external_size; ++i) {
            if (output[i] > data[i]) output[i] = data[i];
         }
         break;
      case sum_reduce:
         for (i = 0; i < d_external_size; ++i) {
            output[i] = output[i] + data[i];
         }
         break;
      default:
         TBOX_ERROR("Library error: d_base_op is somehow corrupted. ");
   }
}

/*
 ***************************************************************************
 * Compute data that is dependent on the group, root and local processes.
 ***************************************************************************
 */
void
AsyncCommGroup::computeDependentData(
   const int* group_ranks,
   const int group_size)
{
   NULL_USE(group_size);
   /*
    * Compute number of descendants in each child branch and in all branches.
    * To find the number of descendants in each branch find the oldest
    * and youngest descendants for each generation.  Add up contribution
    * of all descendant generations.
    */
   d_branch_size_totl = 0;
   unsigned int ic;
   for (ic = 0; ic < d_nchild; ++ic) {
      d_child_data[ic].size = 0;
      int pos_of_child_ic = toChildPosition(toPosition(d_idx), ic);
      int oldest = pos_of_child_ic;
      int yngest = pos_of_child_ic;
      while (oldest < d_group_size) {
         d_child_data[ic].size +=
            (yngest < d_group_size ? yngest : d_group_size - 1) - oldest + 1;
         oldest = toOldest(oldest);
         if (yngest < d_group_size) yngest = toYoungest(yngest);
      }
      d_branch_size_totl += d_child_data[ic].size;
   }

   /*
    * Find ranks of parent and children using specific arithmetic
    * relationships between parent and children positions.
    */
   const int pos = toPosition(d_idx);
   if (pos > 0) {
      d_parent_rank = group_ranks[toIndex((pos - 1)
                                     / static_cast<int>(d_nchild))];
   } else d_parent_rank = -1;
   for (ic = 0; ic < d_nchild; ++ic) {
      const int pos_of_child_ic = toChildPosition(pos, ic);
      if (pos_of_child_ic < d_group_size) {
         d_child_data[ic].rank = group_ranks[toIndex(pos_of_child_ic)];
      } else d_child_data[ic].rank = -1;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommGroup::logCurrentState(
   std::ostream& co) const
{
   SAMRAI_MPI::Request * const req = getRequestPointer();
   co << "State=" << 10 * d_base_op + d_next_task_op
      << "  tag=" << d_mpi_tag
      << "  communicator=" << d_mpi.getCommunicator()
      << "  extern. buff=" << d_external_buf
      << "  size=" << d_external_size
      << "  parent=" << d_parent_rank
      << "  root rank=" << d_root_rank
      << "  use_mpi_collective_for_full_groups="
      << d_use_mpi_collective_for_full_groups
   ;
   for (size_t i = 0; i < d_nchild; ++i) {
      co << "  [" << i << "]=" << d_child_data[i].rank
         << " (" << req[i] << ')';
   }
   co << '\n';
}

/*
 ****************************************************************
 ****************************************************************
 */
void
AsyncCommGroup::setMPITag(
   const int mpi_tag)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI tag is not allowed\n"
         << "during pending communications");
   }
   d_mpi_tag = mpi_tag;
}

/*
 ****************************************************************
 ****************************************************************
 */
void
AsyncCommGroup::setUseMPICollectiveForFullGroups(
   bool use_mpi_collective)
{
   if (!isDone()) {
      TBOX_ERROR("Resetting the MPI collective option is not allowed\n"
         << "during pending communications");
   }
   d_use_mpi_collective_for_full_groups = use_mpi_collective;
}

/*
 ****************************************************************
 ****************************************************************
 */
int
AsyncCommGroup::toPosition(
   int index) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (index < 0 || index >= getGroupSize()) {
      TBOX_ERROR(
         "Invalid index " << index << "\n"
                          << "should be in [0," << d_group_size - 1
                          << "].");
   }
#endif
   const int position = index == 0 ? d_root_idx :
      index == d_root_idx ? 0 :
      index;
   return position;
}

/*
 ****************************************************************
 ****************************************************************
 */
int
AsyncCommGroup::toIndex(
   int position) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (position < 0 || position >= getGroupSize()) {
      TBOX_ERROR(
         "Invalid parent position " << position
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   const int index = position == 0 ? d_root_idx :
      position == d_root_idx ? 0 :
      position;
   return index;
}

/*
 ****************************************************************
 ****************************************************************
 */
int
AsyncCommGroup::toChildPosition(
   int parent_pos,
   int child) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= getGroupSize()) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return parent_pos * static_cast<int>(d_nchild) + 1 + child;
}

/*
 ****************************************************************
 ****************************************************************
 */
int
AsyncCommGroup::toOldest(
   int parent_pos) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= getGroupSize()) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return parent_pos * static_cast<int>(d_nchild) + 1;
}

/*
 ****************************************************************
 ****************************************************************
 */
int
AsyncCommGroup::toYoungest(
   int parent_pos) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (parent_pos < 0 || parent_pos >= getGroupSize()) {
      TBOX_ERROR(
         "Invalid parent position " << parent_pos
                                    << " should be in [" << 0 << ','
                                    << d_group_size - 1 << "].");
   }
#endif
   return (parent_pos + 1) * static_cast<int>(d_nchild);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommGroup::checkMPIParams()
{
   if (getMPITag() < 0) {
      TBOX_ERROR("AsyncCommGroup: Invalid MPI tag value "
         << d_mpi_tag << "\nUse setMPITag() to set it.");
   }
}

AsyncCommGroup::ChildData::ChildData():
   rank(-1),
   size(-1)
{
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
