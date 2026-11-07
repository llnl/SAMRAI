/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Support for coordinating multiple asynchronous communications
 *
 ************************************************************************/
#ifndef included_tbox_AsyncCommStage
#define included_tbox_AsyncCommStage

#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Timer.h"

#include <vector>
#include <list>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Stage multiple non-blocking MPI communications so that codes
 * waiting for them to complete can advance as their respective MPI
 * communication operations are completed.
 *
 * An AsyncCommStage object manages multiple non-blocking MPI
 * communications carried out by AsyncCommStage::Member objects.
 * These classes factor out multiple non-blocking communication codes.
 *
 * By being staged together, Members can overlap their communication
 * waits conveniently.  The stage determines which Member's
 * communication requests are completed.  User code then follows up
 * the communication with whatever was waiting for that operation to
 * complete.
 *
 * The Member class defines the interfaces required to work with an
 * AsyncCommStage.  Communication operations can be a simple
 * non-blocking MPI communication or a complex sequence of MPI
 * communications and local computations.  See Member subclasses for
 * examples.
 *
 * Use this class when you:
 * - Have multiple communications, each performed
 *   by a different and independent AsyncCommStage::Member object.
 * - Want the multiple communications to proceed
 *   independently and asynchronously.
 *
 * Each AsyncCommStage is a registry of Member objects on the stage.
 * To overlap the Members' communications, the stage manages space for
 * their SAMRAI_MPI::Request's and SAMRAI_MPI::Status'es.  The
 * requests are allocated such that a single MPI_Waitany or
 * MPI_Waitsome represents all Members registered with the stage.
 * Thus the communications performed by the Member objects can
 * complete in the order allowed by the arrival of MPI messages.  The
 * exact order is NOT deterministic!
 *
 * To advance the communication operation of any of the allocated
 * Member objects, use advanceAny() or advanceSome().  In general,
 * advanceSome() has better performance than advanceAny() because it
 * avoids the "starvation" problem.  See the MPI documentation for a
 * discussion of starvation.
 *
 * This class supports communication and uses MPI for message passing.
 * If MPI is disabled, the job of this class disappears and the class
 * is effectively empty, except for managing the registration of
 * Member objects.  The public interfaces still remain so the class
 * can compile, but the implementations are trivial.
 */
class AsyncCommStage
{

public:
   /*!
    * @brief Optional object that can be attached to each Member,
    * for use in determining what to do after the Member completes
    * its communication.
    *
    * This class is empty because it serves only to be subclassed.
    *
    * @see Member::setHandler().
    */
   struct Handler {
      virtual ~Handler();
   };

   /*!
    * @brief Something on a stage, using MPI requests and statuses
    * provided by the stage.
    *
    * This class defines the interfaces required to participate
    * in a AsyncCommStage.  Subclasses should implement communication
    * operations using the MPI requests and statuses returned by
    * getRequestPointer() and getStatusPointer().
    */
   struct Member {
      /*!
       * @brief Default constructor, requiring a follow-up call to
       * attachStage() to properly stage the Member.
       */
      Member();

      /*!
       * @brief Initializing constructor.
       *
       * Same as default construction plus attachStage(), plus
       * setHandler().
       */
      Member(
         const size_t nreq,
         AsyncCommStage* stage,
         Handler* handler);

      /*!
       * @brief Destructor detach the Member from its stage.
       * Memory allocated by the stage to support the Member
       * will be recycled.
       *
       * @pre !hasPendingRequests()
       */
      virtual ~Member();

      /*!
       * @brief Get the number of requests for the stage Member.
       *
       * This is the number used to initialize the object.  If the
       * Member has not been initialized to a stage, zero is returned.
       */
      size_t
      numberOfRequests() const
      {
         return d_nreq;
      }

      /*!
       * @brief Return the number of pending SAMRAI_MPI::Request
       * objects for the Member.
       */
      size_t
      numberOfPendingRequests() const;

      /*!
       * @brief Return whether the Member has some pending communication
       * requests.
       */
      bool
      hasPendingRequests() const;

      /*!
       * @brief Check whether entire operation is done.
       *
       * This should not just check the MPI requests--that is the job
       * for hasPendingRequests().  It should check the entire
       * communication operation (which may have follow-up
       * communications and computations).  If MPI communications are
       * completed but have not been followed through by a
       * completeCurrentOperation(), it should return false.
       */
      virtual bool
      isDone() const = 0;

      /*!
       * @brief If current MPI requests have completed, proceed to the
       * next communication wait (or to completion if no more
       * communication is needed).
       */
      virtual bool
      proceedToNextWait() = 0;

      /*!
       * @brief Complete entire operation, including waiting for
       * messages perform follow-up operations and waiting for
       * follow-ups to complete.
       *
       * After this method completes, isDone() should return true.
       *
       * This method should assume that all needed communications can
       * finish within this method.  Otherwise the method may hang.
       */
      virtual void
      completeCurrentOperation() = 0;

      /*!
       * @brief Set a user-defined Handler.
       *
       * A Handler is just some user-defined data.
       *
       * When the stage returns completed Members, the Handler makes
       * it easier to determine what was waiting for the Member to
       * complete.  If you do not need to find the handler this way,
       * you do not need to set it.  Typically, Handlers are only
       * needed if a stage is simultaneously coordinating multiple
       * Members per processor.
       *
       * By design, a non-const Handler is returned, even by a const
       * Member.  This is because the Handler is expected to follow
       * the communication operation with something non-trivial and
       * possibly state-changing.
       */
      void
      setHandler(
         Handler* handler)
      {
         d_handler = handler;
      }

      /*!
       * @brief Regurgitate the handler from setHandler().
       */
      Handler *
      getHandler() const
      {
         return d_handler;
      }

      /*!
       * @brief Push this onto its stage's list of completed Members.
       *
       * This causes the member to be returned by a call to
       * AsyncCommStage::popCompletionQueue(), eventually.
       *
       * @pre isDone()
       */
      void
      pushToCompletionQueue();

      /*!
       * @brief Yank this member from the stage's list of completed
       * Members so it would not be returned by
       * AsyncCommStage::popCompletionQueue().
       */
      void
      yankFromCompletionQueue()
      {
         d_stage->privateYankFromCompletionQueue(*this);
      }

      /*!
       * @brief Returns true if there is a stage associated with this member.
       */
      bool
      hasStage() const
      {
         return d_stage != 0;
      }

protected:
      /*!
       * @brief Associate the member with a stage.
       *
       * Specify the stage to use and the number of SAMRAI_MPI::Request
       * objects needed from the stage.
       *
       * @param[in] nreq Number of SAMRAI_MPI::Requests needed by the member.
       *
       * @param[in] stage The AsyncCommStage to attach to.
       */
      void
      attachStage(
         const size_t nreq,
         AsyncCommStage* stage);

      /*!
       * @brief Disassociate the member from its stage (undo attachStage()).
       */
      void
      detachStage();

      /*!
       * @brief Return the MPI requests provided by the stage for the
       * Member.
       *
       * To work on the stage, you must use these requests in your
       * non-blocking MPI calls.  The requests are a part of an array
       * of requests for all Members on the stage.  Use only
       * getNumberOfRequests() requests.  Terrible bugs will strike if
       * you modify more than what was allocated for the Member.
       *
       * @b Important: Do not save the pointer returned by this
       * method.  The stage's dynamic memory actions may render old
       * pointers invalid.
       *
       * @pre hasStage()
       */
      SAMRAI_MPI::Request *
      getRequestPointer() const;

      /*!
       * @brief Return MPI statuses for this Member of the stage.
       *
       * To work on the stage, you must use these statuses in your MPI
       * calls.  The status attributes are set by the stage only for
       * requests completed by the stage.  The statuses are a part of
       * an array of statuses for all Members on the stage.  Use only
       * getNumberOfRequests() statuses.  Terrible bugs will strike if
       * you modify more than what was allocated for the Member.
       *
       * @b Important: Do not save the pointer returned by this
       * method.  The stage's dynamic memory actions may render old
       * pointers invalid.
       *
       * @pre hasStage()
       */
      SAMRAI_MPI::Status *
      getStatusPointer() const;

      /*!
       * @brief Return the number of requests for this stage Member.
       */
      size_t getNumberOfRequests() const {
         return d_nreq;
      }

private:
      /*!
       * @brief Member is an integral part of the stage code and the
       * stage will set the Member's internal data.
       */
      friend class AsyncCommStage;

      // Unimplemented copy constructor.
      Member(
         const Member& other);

      // Unimplemented assignment operator.
      Member&
      operator = (
         const Member& rhs);

      /*!
       * @brief The stage this Member belongs to.
       *
       * This stage provides SAMRAI_MPI::Request and SAMRAI_MPI::Status
       * objects for the Member.
       */
      AsyncCommStage* d_stage;

      /*!
       * @brief Number of requests reserved on the stage, valid only
       * when d_stage is set.
       *
       * This parameter is always the same as
       * d_stage->numberOfRequests(d_index_on_stage), and exists only
       * for convenience.  Set by the stage when Member is staged or
       * destaged.
       */
      size_t d_nreq;

      /*!
       * @brief Member's index within the stage, valid only when
       * d_stage is set.
       *
       * Set by the stage when Member is staged or destaged.
       */
      size_t d_index_on_stage;

      /*!
       * @brief Pointer to the Member's Handler.
       *
       * @see getHandler()
       */
      Handler* d_handler;
   };

   /*!
    * @brief Construct a stage that may begin allocating and
    * managing Members.
    */
   AsyncCommStage();

   /*!
    * @brief Deallocate Members remaining in the stage and all
    * internal data used to manage Members.
    */
   ~AsyncCommStage();

   /*!
    * @brief Advance to completion one Member (any Member) that is
    * currently waiting for communication to complete.
    *
    * The completed Member is accessible through popCompletionQueue().
    *
    * This method uses MPI_Waitany, which may be prone to starvation.
    * When a process is potentially receiving multiple messages from
    * another processor, it is better to use advanceSome(), which uses
    * MPI_Waitsome, which avoids starvation.
    *
    * @return True if there are still completed Members in the completion queue.
    */
   bool
   advanceAny();

   /*!
    * @brief Advance to completion one or more Members (any Members)
    * that are currently waiting for communication to complete.
    *
    * The completed Members are accessible through popCompletionQueue().
    *
    * @return True if there are still completed Members in the completion queue.
    */
   bool
   advanceSome();

   /*!
    * @brief Advance all pending communications to operations
    *
    * The completed Members are accessible through popCompletionQueue().
    *
    * @return True if there are still completed Members in the completion queue.
    */
   bool
   advanceAll();

   /*
    * @brief Number of completed stage Members in the completion
    * queue.
    *
    * Members that completed their communication operation through
    * advanceAny(), advanceSome() or advanceAll(), and members manually
    * pushed onto the queue by pushToCompletionQueue(), can be accessed
    * through popCompletionQueue().  This method gives the size of
    * that queue.
    */
   size_t
   numberOfCompletedMembers() const
   {
      return d_completed_members.size();
   }

   /*
    * @brief Returns true if there are completed stage Members in the completion
    * queue.
    *
    * Members that completed their communication operation through
    * advanceAny(), advanceSome() or advanceAll() can be accessed
    * through popCompletionQueue().  This method tells if there are
    * any members in that queue.
    */
   bool
   hasCompletedMembers() const
   {
      return !d_completed_members.empty();
   }

   /*!
    * @brief Returns the first completed Member.
    */
   const Member *
   firstCompletedMember() const
   {
      return d_members[d_completed_members.front()];
   }

   /*!
    * @brief Returns the ith managed Member.
    *
    * @param i
    *
    * @pre i < numManagedMembers()
    */
   const Member *
   getMember(size_t i) const
   {
      return d_members[i];
   }

   /*
    * @brief Return a Member that has fully completed its
    * communication operation.
    *
    * This method provides the accessor for stage Members what have
    * completed their communication operations.  The queue is
    * populated by advanceAny(), advanceSome() and advanceAll().  You
    * can also push Members onto the queue using the Member's
    * pushToCompletionQueue().
    *
    * @pre hasCompletedMembers()
    * @pre firstCompletedMember()->isDone()
    */
   Member *
   popCompletionQueue();

   /*!
    * @brief Clear the internal completion queue.
    *
    * @see popCompletionQueue(), Member::pushToCompletionQueue().
    */
   void
   clearCompletionQueue()
   {
      d_completed_members.clear();
   }

   /*!
    * @brief Set optional timer for timing communication wait.
    *
    * The timer is used only for the MPI communication wait.  Timing
    * does not include time the stage Members use to process the
    * messages (the time spent in Member::proceedToNextWait()).
    *
    * If not set, or NULL is given, none will be used.
    */
   void
   setCommunicationWaitTimer(
      const std::shared_ptr<Timer>& communication_timer)
   {
      d_communication_timer = communication_timer;
   }

   /*!
    * @brief Get the number of Members on this stage.
    */
   size_t
   numberOfMembers() const
   {
      return d_member_count;
   }

   /*!
    * @brief Returns the number of staged and destaged Members.
    */
   size_t
   numManagedMembers() const
   {
      return d_members.size();
   }

   /*!
    * @brief Return whether the stage has any pending communication
    * requests.
    */
   bool
   hasPendingRequests() const;

   /*!
    * @brief Return the number of Members that have pending
    * communication.
    */
   size_t
   numberOfPendingMembers() const;

   /*!
    * @brief Return the number of pending SAMRAI_MPI::Request
    * objects on the stage.
    */
   size_t
   numberOfPendingRequests() const;

private:
   /*!
    * @brief Member is a friend so it can access private look-up and
    * stage/destage methods.  This avoids making those private methods
    * public.  This friendship is safe because Member is an integral
    * part of the stage code.
    */
   friend struct Member;

   // Unimplemented copy constructor.
   AsyncCommStage(
      const AsyncCommStage& other);

   // Unimplemented assignment operator.
   AsyncCommStage&
   operator = (
      const AsyncCommStage& rhs);

   //@{
   //! @name Private methods to be called only by Members of the stage.

   /*!
    * @brief Set up a Member to work this stage, initializing mutual
    * references between the stage and the member.
    *
    * @param member Member to work this stage.
    * @param nreq Number of requests needed on the stage.
    *
    * @pre !member->hasStage()
    * @pre nreq >= 1
    */
   void
   privateStageMember(
      Member* member,
      size_t nreq);

   /*!
    * @brief Remove a member from the stage, clearing mutual
    * references between the stage and the member.
    *
    * @param member Member removed from the stage.
    *
    * @pre !member->hasPendingRequests()
    * @pre getMember(member->d_index_on_stage) == member
    */
   void
   privateDestageMember(
      Member* member);

   /*!
    * @brief Get the number of requests for the given Member index.
    *
    * @param index_on_stage Member index.
    *
    * @pre index_on_stage < numManagedMembers()
    * @pre getMember(index_on_stage) != 0
    */
   size_t
   numberOfRequests(
      size_t index_on_stage) const;

   /*!
    * @brief Assert internal data consistency.
    */
   void
   assertDataConsistency() const;

   /*!
    * @brief Lookup and return the request pointer from the stage for
    * the given Member.
    *
    * The given Member MUST have been allocated by the stage.
    * The number of requests that the Member may use is the
    * same as the number originally requested.
    *
    * The pointer is NOT guaranteed to be the same for the life
    * of the Member, as a stage may rearange the array of
    * SAMRAI_MPI::Request objects.  However, the pointer is valid
    * until the next call to privateStageMember().
    *
    * This is private because only friend class Member should
    * use it.
    *
    * @param index_on_stage Member index.
    *
    * @pre index_on_stage < numManagedMembers()
    * @pre getMember(index_on_stage) != 0
    */
   SAMRAI_MPI::Request *
   lookupRequestPointer(
      const size_t index_on_stage) const;

   /*!
    * @brief Look up and return the status pointer from the stage
    * for the given Member.  (Works similarly to lookupRequestPointer().)
    *
    * @param index_on_stage Member index.
    *
    * @pre index_on_stage < numManagedMembers()
    * @pre getMember(index_on_stage) != 0
    */
   SAMRAI_MPI::Status *
   lookupStatusPointer(
      const size_t index_on_stage) const;

   /*!
    * @brief Push the given Member onto the stage's list of completed
    * Members.
    *
    * @param member
    *
    * @pre member.isDone()
    */
   void
   privatePushToCompletionQueue(
      Member& member);

   /*!
    * @brief Yank the given Member from the stage's list of completed
    * Members.
    */
   void
   privateYankFromCompletionQueue(
      Member& member);

   //@}

   /*!
    * @brief Members managed on this stage.
    *
    * Includes destaged Members that are still occupying space on the
    * stage because they are not at the end of the vector and cannot
    * be removed from the vector.
    */
   std::vector<Member *> d_members;

   /*!
    * @brief Number of members.
    *
    * Not necessarily the same as numManagedMembers(), because d_members
    * may have unused slots left behind by destaged members.
    */
   size_t d_member_count;

   /*!
    * @brief SAMRAI_MPI::Request objects used by the Members.
    *
    * This is mutable because the const method getRequestPointer()
    * needs to return a non-const SAMRAI_MPI::Request pointer.  The
    * pointer must be non-const for use in MPI calls.
    * getRequestPointer() should be const because no Member should
    * require a non-const stage just to get the request allocated for
    * it.
    */
   mutable std::vector<SAMRAI_MPI::Request> d_req;

   /*!
    * @brief SAMRAI_MPI::Status objects corresponding to requests on
    * the stage.
    */
   mutable std::vector<SAMRAI_MPI::Status> d_stat;

   //!@brief Map from request index to Member index.
   std::vector<size_t> d_req_to_member;

   /*!
    * @brief Map from Member index to (the Member's first) request index.
    *
    * Provides the index where Member i's requests starts
    * (d_member_to_req[i]) and the number of request it has
    * (d_member_to_req[i+1]-d_member_to_req[i]).  This vector is
    * always one longer than d_members.  The extra item,
    * d_member_to_req[numManagedMembers()], is the total number of
    * requests that the stage has allocated.
    */
   std::vector<size_t> d_member_to_req;

   /*!
    * @brief Members completed by by the stage.
    *
    * Member whose operations are completed through advanceAny(),
    * advanceSome() or advanceAll() are tracked here for access
    * through popCompletionQueue().  This list does NOT include
    * Members completed on their own, because users are expected
    * to track those by themselves.
    */
   std::list<size_t> d_completed_members;

   /*!
    * @brief Members who has completed their operations.
    */

   /*!
    * @brief Timer for communicaiton wait, set by
    * setCommunicationWaitTimer().
    */
   std::shared_ptr<Timer> d_communication_timer;

};

}
}

#endif  // included_tbox_AsyncCommStage
