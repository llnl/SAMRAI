/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Support for coordinating multiple asynchronous communications
 *
 ************************************************************************/
#include "SAMRAI/tbox/AsyncCommStage.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include STL_SSTREAM_HEADER_FILE

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

// #define AsyncCommStage_ExtraDebug

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommStage::AsyncCommStage():
   d_members(),
   d_member_count(0),
   d_req(0),
   d_stat(0),
   d_req_to_member(0),
   d_member_to_req(1, 0)
{
}

/*
 ****************************************************************
 * Make sure the remaining Members are not in the middle
 * of a communication.  Deallocate allocated Members.
 ****************************************************************
 */
AsyncCommStage::~AsyncCommStage()
{
   for (size_t i = 0; i < d_members.size(); ++i) {
      if (d_members[i] != 0) {
         /*
          * Found an undeallocated Member.  Make sure it does not
          * have oustanding requests and deallocate it.
          */
         if (d_members[i]->hasPendingRequests()) {
            TBOX_ERROR("Destructing a stage while some Members\n"
               << "have pending communication leads to\n"
               << "abandoned MPI messages.  Member number "
               << i << "\n"
               << "is not yet done." << std::endl);
         }
         d_members[i] = 0;
      }
   }
}

/*
 ****************************************************************
 * Stage a new communication stage member.  Make sure there
 * is enough space in the arrays, reallocating space if needed.
 * Allocate the new Member and set the internal arrays for it.
 ****************************************************************
 */
void
AsyncCommStage::privateStageMember(
   Member* member,
   size_t nreq)
{
   TBOX_ASSERT(!member->hasStage());   // Double stage not allowed.
#ifdef DEBUG_CHECK_ASSERTIONS
   if (nreq < 1) {
      TBOX_ERROR("Each Member on a stage must have at least one request.\n");
   }
   assertDataConsistency();
#endif

   /*
    * Allocate space at the end of the current arrays for the Member
    * and its needed requests.
    *
    * Recall that d_member_to_req[d_members.size()] is always the
    * current number of requests currently used by the stage.  So
    * d_member_to_req[d_members.size()+1] is set to the new number of
    * requests used by the stage.
    */
   d_members.push_back(member);
   ++d_member_count;

   const size_t cur_total_request_count = d_member_to_req[d_members.size() - 1];
   const size_t new_total_request_count = cur_total_request_count + nreq;
   d_member_to_req.push_back(new_total_request_count);

   d_req.resize(new_total_request_count, MPI_REQUEST_NULL);
   d_stat.resize(new_total_request_count);
   d_req_to_member.resize(new_total_request_count, d_members.size() - 1);
   for (size_t i = cur_total_request_count; i < new_total_request_count; ++i) {
      d_stat[i].MPI_TAG = d_stat[i].MPI_SOURCE = d_stat[i].MPI_ERROR = -1;
   }

   member->d_stage = this;
   member->d_nreq = nreq;
   member->d_index_on_stage = d_members.size() - 1;

#ifdef DEBUG_CHECK_ASSERTIONS
   assertDataConsistency();
#endif
}

/*
 *******************************************************************
 * Remove mutual reference between Member and the stage.
 * 1. Confirm that the Member is currently on the stage
 *    (or else it is an illegal operation).
 * 2. Remove mutual references.
 * 3. Reduce the stage data size by skimming unused space
 *    off the ends of arrays, if possible.
 *******************************************************************
 */
void
AsyncCommStage::privateDestageMember(
   Member* member)
{
   if (member->hasPendingRequests()) {
      TBOX_ERROR("Cannot clear a Member with pending communications.\n"
         << "It would corrupt message passing algorithms.\n");
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assertDataConsistency();
#endif

   if (getMember(member->d_index_on_stage) != member) {
      /*
       * Member was not staged with this AsyncCommStage.  Since staging
       * and destaging are private methods, there must be some logic
       * bug in the library.
       */
      TBOX_ERROR("Library error: An AsyncCommStage cannot destage a Member\n"
         << "that was not staged with it." << std::endl);
   }

   d_members[member->d_index_on_stage] = 0;
   --d_member_count;

   /*
    * Remove the ends of the arrays as much as possible without
    * shifting arrays.  (This is only possible if member is at the
    * end.)
    */
   size_t min_required_len = d_members.size();
   while (min_required_len > 0 && d_members[min_required_len - 1] == 0) {
      --min_required_len;
   }
   if (min_required_len != d_members.size()) {
      d_members.resize(min_required_len, 0);
      d_member_to_req.resize(d_members.size() + 1,
         size_t(MathUtilities<int>::getMax()));

      const size_t new_num_req = d_member_to_req[d_member_to_req.size() - 1];
      d_req_to_member.resize(new_num_req,
         size_t(MathUtilities<int>::getMax()));
      d_req.resize(new_num_req, MPI_REQUEST_NULL);
   }

   member->d_nreq = member->d_index_on_stage = size_t(
            MathUtilities<int>::getMax());
   member->d_stage = 0;
   member->d_handler = 0;

#ifdef DEBUG_CHECK_ASSERTIONS
   assertDataConsistency();
#endif
}

/*
 ****************************************************************
 ****************************************************************
 */
void
AsyncCommStage::assertDataConsistency() const
{
   if (d_members.size() + 1 != d_member_to_req.size()) {
      TBOX_ERROR("d_members.size()=" << d_members.size()
                                     << "+1 is not d_member_to_req.size()="
                                     << d_member_to_req.size() << std::endl);
   }
   if (d_member_to_req[d_member_to_req.size() - 1] != d_req.size()) {
      TBOX_ERROR("d_member_to_req's last entry is bad." << std::endl);
   }

   if (d_members.empty()) {
      return;
   }

   for (size_t i = 0; i < d_members.size() - 1; ++i) {
      if (d_member_to_req[i] >= d_member_to_req[i + 1]) {
         TBOX_ERROR("d_member_to_req out of order at i=" << i << std::endl);
      }
      if (d_members[i] != 0) {
         if (d_members[i]->d_nreq !=
             d_member_to_req[i + 1] - d_member_to_req[i]) {
            TBOX_ERROR("d_members[" << i << "] has bad d_nreq="
                                    << d_members[i]->d_nreq << ", stage value is "
                                    << (d_member_to_req[i + 1] - d_member_to_req[i])
                                    << std::endl);
         }
      }
   }

   size_t member_index = 0;
   size_t number_of_requests = 0;
   for (size_t i = 0; i < d_req_to_member.size(); ++i) {
      if (d_req_to_member[i] == member_index) {
         ++number_of_requests;
      } else {
         if (d_members[member_index] != 0) {
            if (d_members[member_index]->d_nreq != number_of_requests) {
               TBOX_ERROR("d_members[" << member_index << "]->d_nreq is "
                                       << d_members[member_index]->d_nreq
                                       << " while stage claims it should have "
                                       << number_of_requests << " requests"
                                       << std::endl);
            }
         }
         member_index = d_req_to_member[i];
         number_of_requests = 1;
      }
   }
   if (d_members[member_index]->d_nreq != number_of_requests) {
      TBOX_ERROR("d_members[" << member_index << "]->d_nreq is "
                              << d_members[member_index]->d_nreq
                              << " while stage claims it should have "
                              << number_of_requests << " requests"
                              << std::endl);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
AsyncCommStage::privatePushToCompletionQueue(Member& member)
{
   TBOX_ASSERT(member.isDone());
   /*
    * Push member onto d_completed_members, but be sure to avoid
    * duplicating.
    */
   std::list<size_t>::iterator li;
   for (li = d_completed_members.begin();
        li != d_completed_members.end(); ++li) {
      if (member.d_index_on_stage == *li) break;
   }
   if (li == d_completed_members.end()) {
      d_completed_members.push_back(member.d_index_on_stage);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
AsyncCommStage::privateYankFromCompletionQueue(Member& member)
{
   /*
    * Remove member from d_completed_members.
    */
   std::list<size_t>::iterator li;
   for (li = d_completed_members.begin();
        li != d_completed_members.end(); ++li) {
      if (member.d_index_on_stage == *li) {
         d_completed_members.erase(li);
         break;
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommStage::Member *
AsyncCommStage::popCompletionQueue()
{
   if (!hasCompletedMembers()) {
      TBOX_ERROR("AsyncCommStage::popCompletionQueue(): There is no\n"
         << "completed member.  You cannot call this method\n"
         << "when hasCompletedMembers()." << std::endl);
   }
   if (!firstCompletedMember()->isDone()) {
      TBOX_ERROR("AsyncCommStage::popCompletionQueue error:\n"
         << "You asked for a completed AsyncCommStage Member\n"
         << "but its stage has changed to pending since the\n"
         << "stage last identified it as being completed.\n"
         << "This is likely caused by some code re-using the\n"
         << "Member for another operation before poping it\n"
         << "using this method." << std::endl);
   }
   Member* completed = d_members[d_completed_members.front()];
   d_completed_members.pop_front();
   return completed;
}

/*
 ******************************************************************
 * Advance all communication Members by calling advanceSome
 * repeatedly until all outstanding communications are complete.
 ******************************************************************
 */
bool
AsyncCommStage::advanceAll()
{
   while (hasPendingRequests()) {
      advanceSome();
   }
   return !d_completed_members.empty();
}

/*
 ******************************************************************
 * Advance one or more communication Members by using
 * MPI_Waitsome to complete one or more communication requests.
 *
 * Get one or more completed requests and check their communication
 * Members to see if any Member finished its communication operation.
 * If at least one Member finished its communication, return those that
 * finished.  If no Member has finished, repeat until at least one has.
 ******************************************************************
 */
bool
AsyncCommStage::advanceSome()
{
   if (!SAMRAI_MPI::usingMPI()) {
      return false;
   }

   if (d_members.empty()) {
      // Short cut for an empty stage.
      return false;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (unsigned int i = static_cast<unsigned int>(d_member_to_req[d_members.size()]);
        i < d_req.size();
        ++i) {
      if (d_req[i] != MPI_REQUEST_NULL)
         TBOX_WARNING("non-null request above d_n_req." << std::endl);
   }
#endif

   std::vector<int> index(static_cast<int>(d_member_to_req[d_members.size()]));
   std::vector<SAMRAI_MPI::Status> stat(
      static_cast<int>(d_member_to_req[d_members.size()]));

   size_t n_member_completed = 0;
   int n_req_completed = 0;

   do {

      int errf;
      if (d_communication_timer) d_communication_timer->start();
      errf = SAMRAI_MPI::Waitsome(
            static_cast<int>(d_member_to_req[d_members.size()]),
            &d_req[0],
            &n_req_completed,
            &index[0],
            &stat[0]);
      if (d_communication_timer) d_communication_timer->stop();
#ifdef DEBUG_CHECK_ASSERTIONS
      if (n_req_completed <= 0) {
         /*
          * Undocumented feature of some MPI_Waitsome implementations:
          * MPI_Waitsome sets n_req_completed to a negative number
          * if all the input requests are MPI_REQUEST_NULL.
          */
         for (size_t i = 0; i < d_member_to_req[d_members.size()]; ++i) {
            if (d_req[i] != MPI_REQUEST_NULL) {
               TBOX_ERROR("Library error in AsyncCommStage::advanceSome:\n"
                  << "errf = " << errf << '\n'
                  << "MPI_SUCCESS = " << MPI_SUCCESS << '\n'
                  << "MPI_ERR_IN_STATUS = " << MPI_ERR_IN_STATUS << '\n'
                  << "MPI_REQUEST_NULL = " << MPI_REQUEST_NULL << '\n'
                  << "number of requests = "
                  << d_member_to_req[d_members.size()] << '\n'
                  << "d_req.size() = " << d_req.size() << '\n'
                  << "n_req_completed = " << n_req_completed << '\n'
                  << "i = " << i << '\n'
                  );
            }
         }
         for (unsigned int i = static_cast<unsigned int>(d_member_to_req[d_members.size()]);
              i < d_req.size();
              ++i) {
            if (d_req[i] != MPI_REQUEST_NULL)
               TBOX_WARNING("non-null request above d_n_reg." << std::endl);
         }
      }
      if (n_req_completed == 0) {
         TBOX_ASSERT(!hasPendingRequests());
      }
#endif
      if (errf != MPI_SUCCESS) {
         TBOX_ERROR("Error in MPI_Waitsome call.\n"
            << "Error-in-status is "
            << (errf == MPI_ERR_IN_STATUS)
            << '\n');
      }

      /*
       * Construct array of Members with at least one completed
       * request.
       */
      // Number of Members to check with at least one completed request.
      unsigned int n_check_member = 0;

      for (int iout = 0; iout < n_req_completed; ++iout) {

         // Save status of completed request.
         d_stat[index[iout]] = stat[iout];
         /*
          * Change index from request index to Member index.
          * If the Member index is not a duplicate, add it to
          * the list of Members to check (which is actually
          * the same list) and increase n_check_member.
          */
         index[iout] = static_cast<int>(d_req_to_member[index[iout]]);
#ifdef AsyncCommStage_ExtraDebug
         plog << "AsyncCommStage::advanceSome completed:"
              << " tag-" << stat[iout].MPI_TAG
              << " source-" << stat[iout].MPI_SOURCE
              << " for member index " << index[iout]
              << std::endl;
#endif
         unsigned int i;
         for (i = 0; i < n_check_member; ++i) {
            if (index[i] == index[iout]) break;
         }
         if (i == n_check_member) {
            index[n_check_member++] = index[iout];
         }
      }

      /*
       * Check the Members whose requests completed and count up the
       * Members that completed all their communication tasks.
       */
      for (unsigned int imember = 0; imember < n_check_member; ++imember) {
         Member& memberi = *d_members[index[imember]];
         TBOX_ASSERT(!memberi.isDone());
         bool memberi_done = memberi.proceedToNextWait();
#ifdef AsyncCommStage_ExtraDebug
         plog
         << "AsyncCommStage::advanceSome proceedToNextWait for member:"
         << memberi.d_index_on_stage
         << " completion=" << memberi_done
         << std::endl;
#endif
         if (memberi_done) {
            ++n_member_completed;
            TBOX_ASSERT(!memberi.hasPendingRequests());
            privatePushToCompletionQueue(memberi);
         }
      }

   } while (n_req_completed > 0 && n_member_completed == 0);

   return !d_completed_members.empty();
}

/*
 ****************************************************************
 * Advance one communication Member by using MPI_Waitany
 * to complete one requests.
 *
 * Get one completed request and check its communication Member to see
 * if the Member finished its communication operation.  If it finished
 * its operation, return the member.  If not, repeat until one Member
 * has completed its communication operation.
 ****************************************************************
 */
bool
AsyncCommStage::advanceAny()
{
   if (!SAMRAI_MPI::usingMPI()) {
      return false;
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (unsigned int i = static_cast<unsigned int>(d_member_to_req[d_members.size()]);
        i < d_req.size();
        ++i) {
      if (d_req[i] != MPI_REQUEST_NULL)
         TBOX_WARNING("non-null request above d_n_reg." << std::endl);
   }
#endif

   int ireq = MPI_UNDEFINED;
   int member_index_on_stage = -1;
   bool member_done;

   do {

      SAMRAI_MPI::Status mpi_stat;
      int errf;
      if (d_communication_timer) d_communication_timer->start();
      errf = SAMRAI_MPI::Waitany(
            static_cast<int>(d_member_to_req[d_members.size()]),
            &d_req[0],
            &ireq,
            &mpi_stat);
      if (d_communication_timer) d_communication_timer->stop();
      if (errf != MPI_SUCCESS) {
         TBOX_ERROR("Error in MPI_Waitany call.\n"
            << "Error-in-status is "
            << (errf == MPI_ERR_IN_STATUS) << '\n'
            << "MPI_ERROR value is " << mpi_stat.MPI_ERROR
            << '\n');
      }

      if (ireq == MPI_UNDEFINED) {
         // All input requests are completed even before waiting.
         break;
      }

      TBOX_ASSERT(ireq >= 0 &&
         ireq < static_cast<int>(d_member_to_req[d_members.size()]));

      d_stat[ireq] = mpi_stat;
      member_index_on_stage = static_cast<int>(d_req_to_member[ireq]);

      TBOX_ASSERT(member_index_on_stage >= 0 &&
         member_index_on_stage < static_cast<int>(d_members.size()));
      TBOX_ASSERT(d_members[member_index_on_stage] != 0);

      Member* member = d_members[member_index_on_stage];
      /*
       * Member member_index_on_stage had a request completed.
       * See if all of its requests are now completed.
       * If so, run proceedToNextWait() to see if the member is done.
       * (The member may not be done because it may initiate
       * follow-up communication tasks.)
       * Exit the do loop when a Member is actually done
       * with all of its communication tasks.
       */
      const size_t init_req = d_member_to_req[member->d_index_on_stage];
      const size_t term_req = d_member_to_req[member->d_index_on_stage + 1];
      size_t i;
      for (i = init_req; i < term_req; ++i) {
         if (d_req[i] != MPI_REQUEST_NULL) {
            break;
         }
      }
      if (i == term_req) {
         /*
          * Member is done with at least its current communication
          * requests.  Follow-up communication may be launched,
          * so check the return value of proceedToNextWait().
          */
         member_done = d_members[member_index_on_stage]->proceedToNextWait();
      } else {
         member_done = false;
         TBOX_ASSERT(d_members[member_index_on_stage]->hasPendingRequests());
      }

   } while (member_done == false);

   if (member_index_on_stage >= 0) {
      privatePushToCompletionQueue(*d_members[member_index_on_stage]);
   }

   return !d_completed_members.empty();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t
AsyncCommStage::numberOfRequests(
   size_t index_on_stage) const
{
   TBOX_ASSERT(index_on_stage < numManagedMembers());
   TBOX_ASSERT(getMember(index_on_stage) != 0);

   const int init_req = static_cast<int>(d_member_to_req[index_on_stage]);
   const int term_req = static_cast<int>(d_member_to_req[index_on_stage + 1]);
   return term_req - init_req;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
AsyncCommStage::hasPendingRequests() const
{
   bool hasPending = false;
   for (size_t ireq = 0; ireq < d_member_to_req[d_members.size()]; ++ireq) {
      if (d_req[ireq] != MPI_REQUEST_NULL) {
         hasPending = true;
         break;
      }
   }
   return hasPending;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t
AsyncCommStage::numberOfPendingRequests() const
{
   size_t npending = 0;
   for (size_t ireq = 0; ireq < d_member_to_req[d_members.size()]; ++ireq) {
      if (d_req[ireq] != MPI_REQUEST_NULL) {
         ++npending;
      }
   }
   return npending;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t
AsyncCommStage::numberOfPendingMembers() const
{
   size_t nmember = 0;
   for (size_t imember = 0; imember < d_members.size(); ++imember) {
      if (d_members[imember] != 0 &&
          d_members[imember]->hasPendingRequests()) {
         ++nmember;
      }
   }
   return nmember;
}

/*
 ****************************************************************
 * Return the request pointer for a communication Member
 * allocated by this object.
 ****************************************************************
 */
SAMRAI_MPI::Request *
AsyncCommStage::lookupRequestPointer(
   const size_t imember) const
{
   TBOX_ASSERT(imember < numManagedMembers());
   TBOX_ASSERT(getMember(imember) != 0);
   return &d_req[d_member_to_req[imember]];
}

/*
 ****************************************************************
 * Return the status pointer for a communication Member
 * allocated by this object.
 ****************************************************************
 */
SAMRAI_MPI::Status *
AsyncCommStage::lookupStatusPointer(
   const size_t imember) const
{
   TBOX_ASSERT(imember < numManagedMembers());
   TBOX_ASSERT(getMember(imember) != 0);
   return &d_stat[d_member_to_req[imember]];
}

/*
 ****************************************************************
 * Implementations for AsyncCommStage::Member.
 ****************************************************************
 */

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommStage::Member::Member(
   const size_t nreq,
   AsyncCommStage* stage,
   AsyncCommStage::Handler* handler):
   d_stage(0),
   d_nreq(size_t(MathUtilities<int>::getMax())),
   d_index_on_stage(size_t(MathUtilities<int>::getMax())),
   d_handler(handler)
{
   stage->privateStageMember(this, nreq);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommStage::Member::Member():
   d_stage(0),
   d_nreq(size_t(MathUtilities<int>::getMax())),
   d_index_on_stage(size_t(MathUtilities<int>::getMax())),
   d_handler(0)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
AsyncCommStage::Member::~Member()
{
   if (hasPendingRequests()) {
      TBOX_ERROR("Cannot deallocate a Member with pending communications.\n"
         << "It would corrupt message passing algorithms.\n");
   }
   if (d_stage != 0) {
      d_stage->privateDestageMember(this);
   }
   d_handler = 0;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommStage::Member::attachStage(
   const size_t nreq,
   AsyncCommStage* stage)
{
   if (d_stage != 0) {
      // Deregister from current stage.
      d_stage->privateDestageMember(this);
   }
   if (stage != 0) {
      // Register with new stage, if any.
      stage->privateStageMember(this, nreq);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommStage::Member::detachStage()
{
   if (d_stage != 0) {
      // Deregister from current stage.
      d_stage->privateDestageMember(this);
   }
   d_nreq = 0;
   d_stage = 0;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
AsyncCommStage::Member::hasPendingRequests() const
{
   if (d_stage == 0) {
      return false;
   } else {
      SAMRAI_MPI::Request* req = getRequestPointer();
      for (size_t i = 0; i < d_nreq; ++i) {
         if (req[i] != MPI_REQUEST_NULL) return true;
      }
   }
   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t
AsyncCommStage::Member::numberOfPendingRequests() const
{
   size_t npending = 0;
   if (d_stage != 0) {
      SAMRAI_MPI::Request* req = getRequestPointer();
      for (size_t i = 0; i < d_nreq; ++i) {
         if (req[i] != MPI_REQUEST_NULL) ++npending;
      }
   }
   return npending;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_MPI::Request *
AsyncCommStage::Member::getRequestPointer() const
{
   if (!hasStage()) {
      TBOX_ERROR("AssyncCommStage::Member::getRequestPointer():\n"
         << "Empty stage encountered!\n"
         << "This probably means the stage Member has not been placed on a stage.\n"
         << "See documentation for the Membmber's concrete implementation for how\n"
         << "to place the Member on the stage.");

   }
   return d_stage->lookupRequestPointer(d_index_on_stage);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_MPI::Status *
AsyncCommStage::Member::getStatusPointer() const
{
   if (!hasStage()) {
      TBOX_ERROR("AssyncCommStage::Member::getStatusPointer():\n"
         << "Empty stage encountered!\n"
         << "This probably means the stage Member has not been placed on a stage.\n"
         << "See documentation for the Membmber's concrete implementation for how\n"
         << "to place the Member on the stage.");

   }
   return d_stage->lookupStatusPointer(d_index_on_stage);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
AsyncCommStage::Member::pushToCompletionQueue()
{
   if (!isDone()) {
      TBOX_ERROR("AsyncCommStage::Member::pushToCompletionQueue error:\n"
         << "This method may not be called by Members that have not\n"
         << "completed their operation (and returns true from isDone()."
         << std::endl);
   }
   d_stage->privatePushToCompletionQueue(*this);
}

/*
 ****************************************************************
 * Implementations for AsyncCommStage::Member.
 ****************************************************************
 */

AsyncCommStage::Handler::~Handler()
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
