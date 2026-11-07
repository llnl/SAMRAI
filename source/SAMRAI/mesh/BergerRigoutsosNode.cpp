/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Node in asynchronous Berger-Rigoutsos tree
 *
 ************************************************************************/
#include <cstring>
#include <algorithm>

#include "SAMRAI/mesh/BergerRigoutsosNode.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int BergerRigoutsosNode::BAD_INTEGER = -9999999;

/*
 *******************************************************************
 * Construct root node for a single block.
 *******************************************************************
 */
BergerRigoutsosNode::BergerRigoutsosNode(
   BergerRigoutsos* common,
   const hier::Box& box):
   d_pos(1),
   d_common(common),
   d_parent(0),
   d_lft_child(0),
   d_rht_child(0),
   d_box(box),
   d_group(0),
   d_min_box_size(common->d_min_box.getBlockVector(box.getBlockId())),
   d_min_cell_request(1),
   d_mpi_tag(-1),
   d_overlap(tbox::MathUtilities<size_t>::getMax()),
   d_box_acceptance(undetermined),
   d_box_iterator(hier::BoxContainer().end()),
   d_wait_phase(to_be_launched),
   d_send_msg(),
   d_recv_msg(),
   d_comm_group(0),
   d_generation(1),
   d_n_cont(0)
{
   if (box.empty()) {
      TBOX_ERROR("BergerRigoutsosNode: Library error: constructing\n"
         << "root node with an empty box.");
   }

   d_common->incNumNodesConstructed();
   d_common->incNumNodesExisting();
   if (d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
      d_common->incNumNodesOwned();
   }

   if (d_common->d_max_generation < d_generation) {
      d_common->d_max_generation = d_generation;
   }

   if (d_common->d_mpi.getRank() == 0) {
      claimMPITag();
   } else {
      d_mpi_tag = 0;
   }

   /*
    * Set the processor group.
    */
   d_group.resize(d_common->d_mpi.getSize(), BAD_INTEGER);
   for (unsigned int i = 0; i < d_group.size(); ++i) {
      d_group[i] = i;
   }

   if (d_common->d_log_node_history) {
      d_common->writeCounters();
      tbox::plog << "Construct root " << d_generation << ':' << d_pos
                 << ' ' << d_box
                 << ".\n";
   }
}

/*
 *******************************************************************
 * Construct non-root node of the tree.
 *******************************************************************
 */
BergerRigoutsosNode::BergerRigoutsosNode(
   BergerRigoutsos* common_params,
   BergerRigoutsosNode* parent,
   const int child_number):
   d_pos((parent->d_pos > 0 && parent->d_pos <
          tbox::MathUtilities<int>::getMax() / 2) ?
         2 * parent->d_pos + child_number :
         (child_number == 0 ? -1 : -2)),
   d_common(common_params),
   d_parent(parent),
   d_lft_child(0),
   d_rht_child(0),
   d_box(common_params->getDim()),
   d_group(0),
   d_min_box_size(d_parent->d_min_box_size),
   d_mpi_tag(-1),
   d_overlap(tbox::MathUtilities<size_t>::getMax()),
   d_box_acceptance(undetermined),
   d_box_iterator(hier::BoxContainer().end()),
   d_wait_phase(for_data_only),
   d_send_msg(),
   d_recv_msg(),
   d_comm_group(0),
   d_generation(d_parent->d_generation + 1),
   d_n_cont(0)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parent->d_pos >= 0 && d_pos < 0) {
      TBOX_WARNING("Too many generations for node identification.\n"
         << "The node id cannot be increased any further.\n"
         << "This affects only the node id, which is only\n"
         << "used for analysis and debugging and does not\n"
         << "affect the algorithm.\n"
         << "Last valid node id is " << d_parent->d_pos << '\n');
   }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   d_box_iterator = hier::BoxContainer().end();
#endif

   d_common->incNumNodesConstructed();
   d_common->incNumNodesExisting();

   if (d_common->d_max_generation < d_generation) {
      d_common->d_max_generation = d_generation;
   }

   if (d_common->d_log_node_history) {
      d_common->writeCounters();
      tbox::plog << "Construct " << d_generation << ':' << d_pos
                 << ", child of "
                 << d_parent->d_generation << ':' << d_parent->d_pos
                 << "   " << d_parent->d_box
                 << ".\n";
   }
}

BergerRigoutsosNode::~BergerRigoutsosNode()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Forbid deleting a node that is running because there may
    * be pending communication (by the node or its children).
    * Note that this is NOT an extra restriction over the
    * recursive implementation.
    */
   if (d_wait_phase != for_data_only &&
       d_wait_phase != to_be_launched &&
       d_wait_phase != completed) {
      TBOX_ERROR("Should not delete a node that is currently running\n"
         << "the Berger-Rigoutsos algorithm because there\n"
         << "may be pending communications." << std::endl);
   }
#endif

   if (d_comm_group != 0) {
      if (!d_comm_group->isDone()) {
         TBOX_ERROR("Library error: Destructing a node with an unfinished\n"
            << "communication tree is bad because it leaves\n"
            << "pending MPI messages." << std::endl);
      }
      delete d_comm_group;
      d_comm_group = 0;
   }

   if (d_common->d_log_node_history) {
      d_common->writeCounters();
      tbox::plog << "Destruct " << d_generation << ':' << d_pos
                 << "  " << d_box
                 << ".\n";
   }

   d_common->decNumNodesExisting();

   d_wait_phase = deallocated;
}

/*
 ********************************************************************
 * This method looks messy, but it is just the BR agorithm,
 * with multiple pause and continue points implemented by
 * the goto and labels.  Each pause point is accompanied by
 * a line setting d_wait_phase so that the algorithm can
 * continue where it left off when this method is called again.
 * The BR algorithm is not completed until this method returns
 * the WaitPhase value "completed".
 ********************************************************************
 */
BergerRigoutsosNode::WaitPhase
BergerRigoutsosNode::continueAlgorithm()
{
   d_common->d_object_timers->t_continue_algorithm->start();
   ++d_n_cont;

   TBOX_ASSERT(d_parent == 0 || d_parent->d_wait_phase != completed);
   TBOX_ASSERT(inRelaunchQueue(this) == d_common->d_relaunch_queue.end());

   /*
    * Skip right to where we left off,
    * which is specified by the wait phase variable.
    */
   switch (d_wait_phase) {
      case for_data_only:
         TBOX_ERROR("Library error: Attempt to execute data-only node."
         << std::endl);
      case to_be_launched:
         goto TO_BE_LAUNCHED;
      case reduce_histogram:
         goto REDUCE_HISTOGRAM;
      case bcast_acceptability:
         goto BCAST_ACCEPTABILITY;
      case gather_grouping_criteria:
         goto GATHER_GROUPING_CRITERIA;
      case bcast_child_groups:
         goto BCAST_CHILD_GROUPS;
      case run_children:
         goto RUN_CHILDREN;
      case bcast_to_dropouts:
         goto BCAST_TO_DROPOUTS;
      case completed:
         TBOX_ERROR("Library error: Senseless continuation of completed node."
         << std::endl);
         break;
      default:
         TBOX_ERROR("Library error: Nonexistent phase." << std::endl);
   }

   bool sub_completed;

   /*
    * Delegated tasks: Major tasks are delegated to private methods.
    * These methods may check whether the process is the owner or
    * just a contributor and fork appropriately.  The communication
    * checking tasks return whether communication is completed, but
    * they do NOT change the d_wait_phase variable, which is done
    * in this function.
    */

TO_BE_LAUNCHED:

   d_common->incNumNodesActive();

   if (d_common->d_log_node_history) {
      d_common->writeCounters();
      tbox::plog << "Commence " << d_generation << ':' << d_pos
                 << "  " << d_box
                 << "  accept=" << d_box_acceptance
                 << "  ovlap=" << d_overlap
                 << "  owner=" << d_box.getOwnerRank()
                 << "  gsize=" << d_group.size()
                 << ".\n";
   }

   if (d_parent == 0 || d_overlap > 0 || d_common->d_mpi.getRank() == d_box.getOwnerRank()) {

      TBOX_ASSERT(inGroup(d_group));

      // Set up communication group for operations in participating group.
      d_comm_group = new tbox::AsyncCommGroup(
            computeCommunicationTreeDegree(static_cast<int>(d_group.size())),
            &d_common->d_comm_stage,
            this);
      d_comm_group->setUseBlockingSendToParent(false);
      d_comm_group->setGroupAndRootRank(d_common->d_mpi,
         &d_group[0], static_cast<int>(d_group.size()), d_box.getOwnerRank());
      if (d_parent == 0) {
         /*
          * For the global group, MPI collective functions are presumably
          * faster than the peer-to-peer collective implementation in
          * AsyncCommGroup.
          *
          * Enable this mode only for the root node.  Child nodes are
          * not guaranteed to execute the communication operation at
          * the same point on all processors (even if all proccessors
          * participate).
          */
         d_comm_group->setUseMPICollectiveForFullGroups(true);
      }

      d_common->d_object_timers->t_local_tasks->start();
      makeLocalTagHistogram();
      d_common->d_object_timers->t_local_tasks->stop();

      if (d_group.size() > 1) {
         d_common->d_object_timers->t_reduce_histogram->start();
         reduceHistogram_start();
         d_common->incNumNodesCommWait();
REDUCE_HISTOGRAM:
         if (!d_common->d_object_timers->t_reduce_histogram->isRunning())
            d_common->d_object_timers->t_reduce_histogram->start();
         if (d_common->d_algo_advance_mode == BergerRigoutsos::SYNCHRONOUS) {
            d_comm_group->completeCurrentOperation();
         }
         sub_completed = reduceHistogram_check();
         d_common->d_object_timers->t_reduce_histogram->stop();
         if (!sub_completed) {
            d_wait_phase = reduce_histogram;
            goto RETURN;
         }
         d_common->decNumNodesCommWait();
      }

      if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
         /*
          * The owner node saves the tag count.  Participant nodes get
          * tag count from broadcastAcceptability().  This data is just for
          * analysis (not required) and I expect it to have trivial cost.
          */
         int narrowest_dir = 0;
         for (int d = 0; d < d_common->getDim().getValue(); ++d) {
            if (d_histogram[d].size() < d_histogram[narrowest_dir].size())
               narrowest_dir = d;
         }
         d_num_tags = 0;
         for (size_t i = 0; i < d_histogram[narrowest_dir].size(); ++i) {
            d_num_tags += d_histogram[narrowest_dir][i];
         }

         /*
          * If this is the root node, d_num_tags is the total tag count
          * in all nodes.
          */
         if (d_parent == 0) {
            d_common->d_num_tags_in_all_nodes += d_num_tags;
         }
      }

      if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
         d_common->d_object_timers->t_local_tasks->start();
         computeMinimalBoundingBoxForTags();
         acceptOrSplitBox();
         d_common->d_object_timers->t_local_tasks->stop();
         TBOX_ASSERT(boxAccepted() || boxRejected() ||
            (boxHasNoTag() && d_parent == 0));
         if (!boxHasNoTag()) {
            /*
             * A box_level node is created even if box is not acceptable,
             * so that the children can reference its local index in case
             * the box is later accepted based on the combined tolerance
             * of the children.  The node would be erased later if
             * it is not finally accepted.
             */
            createBox();
         }
      }

      if (d_group.size() > 1) {
         d_common->d_object_timers->t_bcast_acceptability->start();
         broadcastAcceptability_start();
         d_common->incNumNodesCommWait();
BCAST_ACCEPTABILITY:
         if (!d_common->d_object_timers->t_bcast_acceptability->isRunning())
            d_common->d_object_timers->t_bcast_acceptability->start();
         if (d_common->d_algo_advance_mode == BergerRigoutsos::SYNCHRONOUS) {
            d_comm_group->completeCurrentOperation();
         }
         sub_completed = broadcastAcceptability_check();
         d_common->d_object_timers->t_bcast_acceptability->stop();
         if (!sub_completed) {
            d_wait_phase = bcast_acceptability;
            goto RETURN;
         }
         d_common->decNumNodesCommWait();
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
         TBOX_ASSERT(d_box_acceptance == accepted_by_calculation ||
            d_box_acceptance == rejected_by_calculation ||
            d_box_acceptance == hasnotag_by_owner);
      } else {
         TBOX_ASSERT(d_box_acceptance == accepted_by_owner ||
            d_box_acceptance == rejected_by_owner ||
            d_box_acceptance == hasnotag_by_owner);
      }
#endif

      /*
       * If this is the root node, d_num_tags is the total tag count
       * in all nodes.
       */
      if (d_parent == 0 && d_common->d_mpi.getRank() != d_box.getOwnerRank()) {
         d_common->d_num_tags_in_all_nodes += d_num_tags;
      }

      if (boxRejected()) {

         /*
          * Compute children groups and owners without assuming
          * entire mesh structure is known locally.
          */
         d_common->d_object_timers->t_local_tasks->start();
         countOverlapWithLocalPatches();
         d_common->d_object_timers->t_local_tasks->stop();

         if (d_group.size() > 1) {
            d_common->d_object_timers->t_gather_grouping_criteria->start();
            gatherGroupingCriteria_start();
            d_common->incNumNodesCommWait();
GATHER_GROUPING_CRITERIA:
            if (!d_common->d_object_timers->t_gather_grouping_criteria->isRunning())
               d_common->d_object_timers->t_gather_grouping_criteria->start();
            if (d_common->d_algo_advance_mode == BergerRigoutsos::SYNCHRONOUS) {
               d_comm_group->completeCurrentOperation();
            }
            sub_completed = gatherGroupingCriteria_check();
            d_common->d_object_timers->t_gather_grouping_criteria->stop();
            if (!sub_completed) {
               d_wait_phase = gather_grouping_criteria;
               goto RETURN;
            }
            d_common->decNumNodesCommWait();
         }

         if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
            d_common->d_object_timers->t_local_tasks->start();
            formChildGroups();
            d_common->d_object_timers->t_local_tasks->stop();
         }

         if (d_group.size() > 1) {
            d_common->d_object_timers->t_bcast_child_groups->start();
            broadcastChildGroups_start();
            d_common->incNumNodesCommWait();
BCAST_CHILD_GROUPS:
            if (!d_common->d_object_timers->t_bcast_child_groups->isRunning())
               d_common->d_object_timers->t_bcast_child_groups->start();
            if (d_common->d_algo_advance_mode == BergerRigoutsos::SYNCHRONOUS) {
               d_comm_group->completeCurrentOperation();
            }
            sub_completed = broadcastChildGroups_check();
            d_common->d_object_timers->t_bcast_child_groups->stop();
            if (!sub_completed) {
               d_wait_phase = bcast_child_groups;
               goto RETURN;
            }
            d_common->decNumNodesCommWait();
         }

         if (d_lft_child->d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
            d_common->incNumNodesOwned();
         }
         if (d_rht_child->d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
            d_common->incNumNodesOwned();
         }

         runChildren_start();
RUN_CHILDREN:
         sub_completed = runChildren_check();
         if (!sub_completed) {
            d_wait_phase = run_children;
            goto RETURN;
         }
      } else if (boxAccepted()) {
         if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
            ++(d_common->d_num_boxes_generated);
         }
      } else {
         // Box has no tag.
      }

      // All done with communication within participating group.
      delete d_comm_group;
      d_comm_group = 0;

   } else {
      /*
       * This process is not in the group that decides on the box for
       * this node.
       */
      TBOX_ASSERT(!inGroup(d_group));
   }

   if (d_parent == 0) {
      /*
       * Compute relationships and set up relationship sharing data.
       * This is usually done by a node's parent in the
       * runChildren_check() method because only the
       * parent can know if the node's box will be
       * kept or recombined with the sibling.
       * But the root node must do this itself because it has no parent.
       */
      if (d_common->d_compute_relationships > 0 && boxAccepted()) {
         computeNewNeighborhoodSets();
      }
   }

   TBOX_ASSERT(d_lft_child == 0);
   TBOX_ASSERT(d_rht_child == 0);
   // TBOX_ASSERT( ! inRelaunchQueue(this) );
   TBOX_ASSERT(inRelaunchQueue(this) == d_common->d_relaunch_queue.end());

   /*
    * Broadcast the result to dropouts.
    * Dropout processes are those that participated in the
    * parent but not in this node.  They need the
    * result to perform combined efficiency check for the
    * parent.
    *
    * Processes that should participate in the dropout broadcast
    * are the dropouts (processes with zero overlap) and the owner.
    *
    * Broadcast to dropouts is only needed if:
    *
    *    - In multi-owner mode and relationship-computing mode.
    *      In single-owner mode, only the original owner needs
    *      the final result, and it participates everywhere,
    *      so there is no need for this phase.
    *      When computing relationships, participant processors must
    *      know results to do recombination check, to determine
    *      if parent box is preferred.
    *
    *    - This is NOT the root node.  The root node
    *      has no parent and no corresponding dropout group.
    *
    *    - Dropout group is not empty.  Number of dropouts
    *      is the difference between parent group size and this
    *      group size.
    */
   if (d_overlap == 0 || d_common->d_mpi.getRank() == d_box.getOwnerRank()) {

      if ((d_common->d_owner_mode != BergerRigoutsos::SINGLE_OWNER ||
           d_common->d_compute_relationships > 0) &&
          d_parent != 0 &&
          d_parent->d_group.size() > d_group.size()) {

         d_common->d_object_timers->t_bcast_to_dropouts->start();
         {
            // Create the communication group for the dropouts.
            BergerRigoutsos::VectorOfInts dropouts(0);
            d_common->d_object_timers->t_local_tasks->start();
            computeDropoutGroup(d_parent->d_group,
               d_group,
               dropouts,
               d_box.getOwnerRank());
            d_comm_group = new tbox::AsyncCommGroup(
                  computeCommunicationTreeDegree(
                     static_cast<int>(d_group.size())),
                  &d_common->d_comm_stage,
                  this);
            d_comm_group->setUseBlockingSendToParent(false);
            d_comm_group->setGroupAndRootIndex(d_common->d_mpi,
               &dropouts[0], static_cast<int>(dropouts.size()), 0);
            d_common->d_object_timers->t_local_tasks->stop();
         }

         broadcastToDropouts_start();
         d_common->incNumNodesCommWait();
BCAST_TO_DROPOUTS:
         if (!d_common->d_object_timers->t_bcast_to_dropouts->isRunning())
            d_common->d_object_timers->t_bcast_to_dropouts->start();
         sub_completed = broadcastToDropouts_check();
         d_common->d_object_timers->t_bcast_to_dropouts->stop();

         if (!sub_completed) {
            d_wait_phase = bcast_to_dropouts;
            goto RETURN;
         }

         d_common->decNumNodesCommWait();

         if (d_common->d_log_node_history && d_common->d_mpi.getRank() != d_box.getOwnerRank()) {
            d_common->writeCounters();
            tbox::plog << "DO Recv " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << "  accept=" << d_box_acceptance
                       << ".\n";
         }

         delete d_comm_group;
         d_comm_group = 0;
      }
   }

   d_wait_phase = completed;

   if (d_comm_group != 0) {
      // No further communication.  Deallocate the communication group.
      delete d_comm_group;
      d_comm_group = 0;
   }

   TBOX_ASSERT(d_common->d_num_nodes_owned >= 0);

   // Adjust counters.
   d_common->decNumNodesActive();
   d_common->incNumNodesCompleted();
   if (d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
      d_common->decNumNodesOwned();
   }
   d_common->incNumContinues(d_n_cont);

   if (d_common->d_log_node_history) {
      d_common->writeCounters();
      tbox::plog << "Complete " << d_generation << ':' << d_pos
                 << "  " << d_box
                 << "  accept=" << d_box_acceptance
                 << ".\n";
   }

   /*
    * Recall that a tree node waiting for its children
    * is not placed in the relaunch queue (because it is
    * pointless to relaunch it until the children are completed).
    * Therefore, to eventually continue that node, its last
    * child to complete must put it on the queue.  If this node
    * and its sibling are completed, put the parent on the FRONT
    * queue to be checked immediately (required for synchronous
    * mode).
    */
   if (d_parent != 0 &&
       d_parent->d_lft_child->d_wait_phase == completed &&
       d_parent->d_rht_child->d_wait_phase == completed) {
      TBOX_ASSERT(d_parent->d_wait_phase == run_children);
      // TBOX_ASSERT( ! inRelaunchQueue(d_parent) );
      TBOX_ASSERT(inRelaunchQueue(d_parent) == d_common->d_relaunch_queue.end());
      // d_common->d_relaunch_queue.push_front(d_parent);
      d_common->prependQueue(d_parent);
      if (d_common->d_log_node_history) {
         d_common->writeCounters();
         tbox::plog << "Parent " << d_parent->d_generation << ':'
                    << d_parent->d_pos
                    << " awoken by last child of "
                    << d_parent->d_lft_child->d_generation << ':'
                    << d_parent->d_lft_child->d_pos
                    << ", "
                    << d_parent->d_rht_child->d_generation << ':'
                    << d_parent->d_rht_child->d_pos
                    << " queue size " << d_common->d_relaunch_queue.size()
                    << ".\n";
      }
   }

RETURN:

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_wait_phase != completed && d_wait_phase != run_children) {
      TBOX_ASSERT(!d_comm_group->isDone());
      TBOX_ASSERT(d_common->d_comm_stage.hasPendingRequests());
   }
   if (d_wait_phase == run_children) {
      // TBOX_ASSERT( ! d_relaunch_queue.empty() );
      TBOX_ASSERT(!d_common->d_relaunch_queue.empty());
   }
#endif

   d_common->d_object_timers->t_continue_algorithm->stop();

   return d_wait_phase;
}

void
BergerRigoutsosNode::runChildren_start()
{
   /*
    * Children were created to store temporary data
    * and determine participation. Now, run them.
    */

   /*
    * Should only be here if box is rejected based on calculation.
    */
   TBOX_ASSERT(d_box_acceptance == rejected_by_calculation ||
      d_box_acceptance == rejected_by_owner);

   d_lft_child->d_wait_phase = to_be_launched;
   d_rht_child->d_wait_phase = to_be_launched;

   /*
    * Queue the children so they get executed.
    * Put them at the front so that in synchronous
    * mode, they can complete first before moving
    * to another task (important in synchronous mode).
    * It also does not hurt to put children at the
    * front of the queue because they have
    * immediate computation (compute histogram)
    * to perform.  Put the left child in front
    * of the right to more closely match the
    * progression of the recursive BR (not essential).
    */
   // d_common->d_relaunch_queue.push_front(d_rht_child);
   // d_common->d_relaunch_queue.push_front(d_lft_child);
   d_common->prependQueue(d_rht_child, d_lft_child);
}

/*
 ********************************************************************
 * Check for combined tolerance.
 * If both children accepted their boxes without further splitting
 * but their combined efficiency is not good enough to make
 * the splitting worth accepting, use the current box instead
 * of the children boxes.  Otherwise, use the children boxes.
 ********************************************************************
 */
bool
BergerRigoutsosNode::runChildren_check()
{
   if (d_lft_child->d_wait_phase != completed ||
       d_rht_child->d_wait_phase != completed) {
      return false;
   }

   const double combine_reduction =
      double(d_lft_child->d_box.size()
             + d_rht_child->d_box.size()) / static_cast<double>(d_box.size());
   if (d_lft_child->boxAccepted() &&
       d_rht_child->boxAccepted() &&
       d_box.numberCells() <= d_common->d_max_box_size &&
       (combine_reduction >= d_common->getCombineEfficiency(d_common->d_level_number))) {

      // Discard childrens' graph nodes in favor of recombination.

      d_box_acceptance = accepted_by_recombination;

      if (d_common->d_log_node_history) {
         d_common->writeCounters();
         tbox::plog << "Recombine " << d_generation << ':' << d_pos
                    << " insufficient reduction of " << combine_reduction
                    << "  " << d_box
                    << " <= " << d_lft_child->d_box
                    << " + " << d_rht_child->d_box
                    << "  " << "accept=" << d_box_acceptance
                    << ".\n";
      }

      if (d_lft_child->d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
         d_lft_child->eraseBox();
         d_lft_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->d_num_boxes_generated);
      }

      if (d_rht_child->d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
         d_rht_child->eraseBox();
         d_rht_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->d_num_boxes_generated);
      }

      if (d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
         ++(d_common->d_num_boxes_generated);
      }

   } else {

      // Accept childrens' results, discarding graph node.

      if (d_box.getOwnerRank() == d_common->d_mpi.getRank()) {
         eraseBox();
      }
      if (d_common->d_compute_relationships > 0) {
         if (d_lft_child->boxAccepted() &&
             d_lft_child->d_box_acceptance != accepted_by_dropout_bcast) {
            d_lft_child->computeNewNeighborhoodSets();
         }
         if (d_rht_child->boxAccepted() &&
             d_rht_child->d_box_acceptance != accepted_by_dropout_bcast) {
            d_rht_child->computeNewNeighborhoodSets();
         }
         if (d_common->d_log_node_history) {
            d_common->writeCounters();
            tbox::plog << "Discard " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << " => " << d_lft_child->d_box
                       << " + " << d_rht_child->d_box
                       << "  " << "accept=" << d_box_acceptance
                       << ".\n";
         }
      }

   }

   /*
    * No longer need children nodes after this point.
    */
   delete d_lft_child;
   delete d_rht_child;
   d_lft_child = 0;
   d_rht_child = 0;

   return true;
}

/*
 ********************************************************************
 *
 * Asynchronous methods: these methods have _start and _check
 * suffices.  They involve initiating some task and checking
 * whether that task is completed.
 *
 ********************************************************************
 */

void
BergerRigoutsosNode::reduceHistogram_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + reduce_histogram_tag);
   const int hist_size = getHistogramBufferSize(d_box);
   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      d_recv_msg.resize(hist_size, BAD_INTEGER);
      putHistogramToBuffer(&d_recv_msg[0]);
      d_comm_group->beginSumReduce(&d_recv_msg[0], hist_size);
   } else {
      d_send_msg.resize(hist_size, BAD_INTEGER);
      putHistogramToBuffer(&d_send_msg[0]);
      d_comm_group->beginSumReduce(&d_send_msg[0], hist_size);
   }
}

bool
BergerRigoutsosNode::reduceHistogram_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->proceedToNextWait();
   if (d_comm_group->isDone() && d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      getHistogramFromBuffer(&d_recv_msg[0]);
   }
   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::broadcastAcceptability_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + bcast_acceptability_tag);
   /*
    * Items communicated:
    * - local index of node
    * - whether box is accepted
    * - in case box is accepted:
    *   . box (which may have been trimmed to minimal tag bounding box)
    * - in case box is rejected:
    *   . left/right child boxes
    *   . left/right child MPI tags
    */

   const int buffer_size = 1          // Number of tags in candidate
      + 1                             // Acceptability flag.
      + 1                             // Local index of node.
      + getDim().getValue() * 2       // Box.
      + getDim().getValue() * 4       // Children boxes.
      + 2                             // Children MPI tags
   ;

   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      TBOX_ASSERT(d_box_acceptance == rejected_by_calculation ||
         d_box_acceptance == accepted_by_calculation ||
         (d_parent == 0 && d_box_acceptance == hasnotag_by_owner));
      d_send_msg.resize(buffer_size, BAD_INTEGER);
      int* ptr = &d_send_msg[0];
      *(ptr++) = d_num_tags;
      *(ptr++) = d_box_acceptance >= 0 ?
         d_box_acceptance + 2 /* indicate remote decision */ :
         d_box_acceptance;
      if (!boxHasNoTag()) {
         *(ptr++) = d_box.getLocalId().getValue();
         ptr = putBoxToBuffer(d_box, ptr);
         if (boxRejected()) {
            ptr = putBoxToBuffer(d_lft_child->d_box, ptr);
            ptr = putBoxToBuffer(d_rht_child->d_box, ptr);
            *(ptr++) = d_lft_child->d_mpi_tag;
            *(ptr++) = d_rht_child->d_mpi_tag;
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      else {
         // This may not be needed now that the messages are in vector<int>.
         // Suppress memory check warnings about uninitialized data.
         for (size_t c = ptr - (&d_send_msg[0]); c < d_send_msg.size(); ++c) {
            d_send_msg[c] = -1;
         }
      }
#endif
      d_comm_group->beginBcast(&d_send_msg[0], buffer_size);
   } else {
      d_recv_msg.resize(buffer_size, BAD_INTEGER);
      d_comm_group->beginBcast(&d_recv_msg[0], buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastAcceptability_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->checkBcast();
   if (d_comm_group->isDone() && d_common->d_mpi.getRank() != d_box.getOwnerRank()) {

      int* ptr = &d_recv_msg[0];

      d_num_tags = *(ptr++);

      d_box_acceptance = intToBoxAcceptance(*(ptr++));
      TBOX_ASSERT(boxAccepted() || boxRejected() ||
         (boxHasNoTag() && d_parent == 0));
      if (!boxHasNoTag()) {
         const hier::LocalId accepted_box_local_id(*(ptr++));
         ptr = getBoxFromBuffer(d_box, ptr);
         d_box.initialize(d_box, accepted_box_local_id, d_box.getOwnerRank());   // Reset local id.
         /*
          * Do not check for min_box violation in root node.  That
          * check should be done outside of this class in order to
          * have flexibility regarding how to handle it.
          */
         TBOX_ASSERT(d_parent == 0 || d_box.numberCells() >= d_min_box_size);
      }

      if (boxRejected()) {

         /*
          * The owner formed its children earlier so it can
          * use their parameters while determining which to run.
          * Contributors create the children when the receive
          * the d_box_acceptance flag indicates that further
          * branching is required.
          */
         d_lft_child = new BergerRigoutsosNode(d_common, this, 0);
         d_rht_child = new BergerRigoutsosNode(d_common, this, 1);

         ptr = getBoxFromBuffer(d_lft_child->d_box, ptr);
         ptr = getBoxFromBuffer(d_rht_child->d_box, ptr);
         d_lft_child->d_box.setBlockId(d_box.getBlockId());
         d_rht_child->d_box.setBlockId(d_box.getBlockId());

         d_lft_child->d_mpi_tag = *(ptr++);
         d_rht_child->d_mpi_tag = *(ptr++);

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_box.numberCells() >= d_min_box_size) {
            TBOX_ASSERT(d_lft_child->d_box.numberCells() >= d_min_box_size);
            TBOX_ASSERT(d_rht_child->d_box.numberCells() >= d_min_box_size);
         }
#endif
         TBOX_ASSERT(d_lft_child->d_mpi_tag > -1);
         TBOX_ASSERT(d_rht_child->d_mpi_tag > -1);
         if (d_common->d_log_node_history) {
            d_common->writeCounters();
            tbox::plog << "Rm Split " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << " => " << d_lft_child->d_box
                       << " + " << d_rht_child->d_box
                       << ".\n";
         }

      } else {
         if (d_common->d_log_node_history) {
            d_common->writeCounters();
            tbox::plog << "Rm Accepted " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << "  accept=" << d_box_acceptance << ".\n";
         }
      }
   }
   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::gatherGroupingCriteria_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + gather_grouping_criteria_tag);

   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      d_recv_msg.resize(4 * d_group.size(), BAD_INTEGER);
      d_comm_group->beginGather(&d_recv_msg[0], 4);
   } else {
      d_send_msg.resize(4, BAD_INTEGER);
      // TODO: Change message buffers to MessageStream to avoid limiting d_overlap to integer size.
      d_send_msg[0] = static_cast<int>(d_lft_child->d_overlap);
      d_send_msg[1] = static_cast<int>(d_rht_child->d_overlap);
      // Use negative burden measures for uniformity of criteria comparison.
      d_send_msg[2] = -d_common->d_num_nodes_owned;
      d_send_msg[3] = -d_common->d_num_nodes_active;
      d_comm_group->beginGather(&d_send_msg[0], 4);
   }
}

void
BergerRigoutsosNode::broadcastChildGroups_start()
{
   if (d_group.size() == 1) {
      return;
   }
   /*
    * Items communicated:
    * - left/right owner
    * - left/right group
    */
   d_comm_group->setMPITag(d_mpi_tag + bcast_child_groups_tag);

   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {

      /*
       * When d_parent == 0, use d_comm_group's MPI collective call option.
       * The option uses MPI_Bcast, which requires the buffer size is the same
       * on all processors.  When this is not the case, use the child group
       * sizes to save memory and possibly improve performance.
       */
      const int buffer_size = 2                // Left/right owners.
         + 2                                   // Left/right group sizes.
         + (d_parent == 0 ? static_cast<int>(d_group.size())
            : static_cast<int>(d_lft_child->d_group.size()))    // Left group.
         + (d_parent == 0 ? static_cast<int>(d_group.size())
            : static_cast<int>(d_rht_child->d_group.size()))    // Right group.
      ;

      d_send_msg.resize(buffer_size, BAD_INTEGER);
      int* ptr = &d_send_msg[0];

      *(ptr++) = d_lft_child->d_box.getOwnerRank();
      *(ptr++) = static_cast<int>(d_lft_child->d_group.size());
      for (size_t i = 0; i < d_lft_child->d_group.size(); ++i) {
         *(ptr++) = d_lft_child->d_group[i];
      }
      *(ptr++) = d_rht_child->d_box.getOwnerRank();
      *(ptr++) = static_cast<int>(d_rht_child->d_group.size());
      for (size_t i = 0; i < d_rht_child->d_group.size(); ++i) {
         *(ptr++) = d_rht_child->d_group[i];
      }
      if (d_parent == 0) {
         // Initialize unused data to avoid warnings and weird numbers.
         for (size_t i =
                 (d_lft_child->d_group.size() + d_rht_child->d_group.size());
              i < 2 * d_group.size(); ++i) {
            *(ptr++) = -1;
         }
      }

      d_comm_group->beginBcast(&d_send_msg[0], buffer_size);
   } else {
      const int buffer_size = 2                // Left/right owners.
         + 2                                   // Left/right group sizes.
         + 2 * static_cast<int>(d_group.size())   // Left/right groups.
      ;
      d_recv_msg.resize(buffer_size, BAD_INTEGER);

      d_comm_group->beginBcast(&d_recv_msg[0], buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastChildGroups_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->checkBcast();
   if (d_comm_group->isDone() && d_common->d_mpi.getRank() != d_box.getOwnerRank()) {

      int* ptr = &d_recv_msg[0];

      int lft_owner = *(ptr++);
      d_lft_child->d_group.resize(*(ptr++), BAD_INTEGER);
      for (size_t i = 0; i < d_lft_child->d_group.size(); ++i) {
         d_lft_child->d_group[i] = *(ptr++);
      }
      int rht_owner = *(ptr++);
      d_rht_child->d_group.resize(*(ptr++), BAD_INTEGER);
      for (size_t i = 0; i < d_rht_child->d_group.size(); ++i) {
         d_rht_child->d_group[i] = *(ptr++);
      }

      d_lft_child->d_box.initialize(d_lft_child->d_box,
         d_lft_child->d_box.getLocalId(),
         lft_owner);
      d_rht_child->d_box.initialize(d_rht_child->d_box,
         d_rht_child->d_box.getLocalId(),
         rht_owner);
      TBOX_ASSERT(d_lft_child->d_box.getOwnerRank() >= 0);
      TBOX_ASSERT(d_lft_child->d_group.size() > 0);
      TBOX_ASSERT((d_lft_child->d_overlap > 0) ==
         inGroup(d_lft_child->d_group));
      TBOX_ASSERT(d_rht_child->d_box.getOwnerRank() >= 0);
      TBOX_ASSERT(d_rht_child->d_group.size() > 0);
      TBOX_ASSERT((d_rht_child->d_overlap > 0) ==
         inGroup(d_rht_child->d_group));

   }

   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::broadcastToDropouts_start()
{
   TBOX_ASSERT(d_common->d_mpi.getRank() == d_box.getOwnerRank() || d_overlap == 0);
   d_comm_group->setMPITag(d_mpi_tag + bcast_to_dropouts_tag);

   const int buffer_size = 1      // d_box_acceptance
      + 1                         // local index of graph node
      + d_common->getDim().getValue() * 2   // d_box (in case it got reduced)
   ;
   d_send_msg.clear();
   d_recv_msg.clear();
   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      d_send_msg.resize(buffer_size, BAD_INTEGER);
      d_send_msg[0] = d_box_acceptance;
      d_send_msg[1] = d_box.getLocalId().getValue();
      putBoxToBuffer(d_box, &d_send_msg[2]);
      d_comm_group->beginBcast(&d_send_msg[0],
         buffer_size);
   } else {
      d_recv_msg.resize(buffer_size, BAD_INTEGER);
      d_comm_group->beginBcast(&d_recv_msg[0],
         buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastToDropouts_check()
{
   TBOX_ASSERT(d_common->d_mpi.getRank() == d_box.getOwnerRank() || d_overlap == 0);
   d_comm_group->checkBcast();
   if (d_comm_group->isDone()) {
      if (d_common->d_mpi.getRank() != d_box.getOwnerRank()) {
         /*
          * We check for the case of the box having no tags,
          * to keeps things explicit and help detect bugs.
          * But in fact, having no tags is impossible
          * in the broadcastToDropout step, because it is
          * only possible for the root node,
          * which has no dropout group.
          */
         TBOX_ASSERT(d_recv_msg[0] >= 0);

         d_box_acceptance = intToBoxAcceptance((d_recv_msg[0] % 2)
               + rejected_by_dropout_bcast);
         const hier::LocalId accepted_box_local_id(d_recv_msg[1]);
         getBoxFromBuffer(d_box, &d_recv_msg[2]);
         /*
          * Do not check for min_box violation in root node.  That
          * check should be done outside of this class in order to
          * have flexibility regarding how to handle it.
          */
         TBOX_ASSERT(d_parent == 0 || d_box.numberCells() >= d_min_box_size);
         d_box.initialize( d_box, accepted_box_local_id, d_box.getOwnerRank() ); // Reset local id.
      }
   }
   return d_comm_group->isDone();
}

/*
 ********************************************************************
 * Utility computations using local data.
 ********************************************************************
 */

void
BergerRigoutsosNode::makeLocalTagHistogram()
{
   RANGE_PUSH("make-histogram", 4);
   d_common->d_object_timers->t_local_histogram->start();

   /*
    * Compute the histogram size and allocate space for it.
    */
   for (tbox::Dimension::dir_t d = 0; d < d_common->getDim().getValue(); ++d) {
      TBOX_ASSERT(d_box.numberCells(d) > 0);
      d_histogram[d].clear();
      d_histogram[d].insert(d_histogram[d].end(), d_box.numberCells(d), 0);
   }

   /*
    * Accumulate tag counts in the histogram variable.
    */
   const hier::PatchLevel& tag_level = *d_common->d_tag_level;
   for (hier::PatchLevel::iterator ip(tag_level.begin());
        ip != tag_level.end(); ++ip) {
      hier::Patch& patch = **ip;

      const hier::BlockId& block_id = patch.getBox().getBlockId();

      if (block_id == d_box.getBlockId()) {
         const hier::Box intersection = patch.getBox() * d_box;
         const hier::Index& lower = d_box.lower();

         if (!(intersection.empty())) {

            std::shared_ptr<pdat::CellData<int> > tag_data_(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                  patch.getPatchData(d_common->d_tag_data_index)));

            TBOX_ASSERT(tag_data_);

            pdat::CellData<int>& tag_data = *tag_data_;
#if defined(HAVE_RAJA)
            tbox::parallel_synchronize();
#endif

            pdat::CellIterator ciend(pdat::CellGeometry::end(intersection));
            for (pdat::CellIterator ci(pdat::CellGeometry::begin(intersection));
                 ci != ciend; ++ci) {
               if (tag_data(*ci) == d_common->d_tag_val) {
                  const hier::Index& idx = *ci;
                  for (int d = 0; d < d_common->getDim().getValue(); ++d) {
                     ++(d_histogram[d][idx(d) - lower(d)]);
                  }
               }
            }
         }
      }
   }
   d_common->d_object_timers->t_local_histogram->stop();
   RANGE_POP;
}

/*
 ********************************************************************
 * Change d_box to that of the minimal bounding box for tags.
 * If d_box is changed, reduce d_histogram to new d_box.
 ********************************************************************
 */
void
BergerRigoutsosNode::computeMinimalBoundingBoxForTags()
{
   TBOX_ASSERT(!d_box.empty());

   hier::Index new_lower = d_box.lower();
   hier::Index new_upper = d_box.upper();

   const hier::IntVector& min_box = d_min_box_size;
   hier::IntVector box_size = d_box.numberCells();

   /*
    * Bring the lower side of the box up past untagged index planes.
    * Bring the upper side of the box down past untagged index planes.
    * Do not make the box smaller than the min_box requirement.
    */
   size_t num_cells = box_size.getProduct();

   d_min_cell_request = (d_common->d_min_cell_request / d_common->d_ratio.getProduct(d_box.getBlockId()));

   for (tbox::Dimension::dir_t d = 0; d < d_common->getDim().getValue(); ++d) {
      TBOX_ASSERT(d_histogram[d].size() != 0);
      int* histogram_beg = &d_histogram[d][0];
      int* histogram_end = histogram_beg + d_box.numberCells(d) - 1;
      while (*histogram_beg == 0 &&
             box_size(d) > min_box(d) && num_cells > d_min_cell_request) {
         ++new_lower(d);
         ++histogram_beg;
         num_cells /= box_size(d);
         --box_size(d);
         num_cells *= box_size(d);
      }
      while (*histogram_end == 0 &&
             box_size(d) > min_box(d) && num_cells > d_min_cell_request) {
         --new_upper(d);
         --histogram_end;
         num_cells /= box_size(d);
         --box_size(d);
         num_cells *= box_size(d);
      }
   }

   const hier::Box new_box(new_lower, new_upper, d_box.getBlockId());
   const hier::IntVector new_size = new_box.numberCells();

   if (!new_box.isSpatiallyEqual(d_box)) {
      /*
       * Do not check for min_box violation in root node.  That
       * check should be done outside of this class in order to
       * have flexibility regarding how to handle it.
       */
      TBOX_ASSERT(d_parent == 0 || new_box.numberCells() >= min_box);
      /*
       * Save tagged part of the current histogram and reset the box.
       * Is this step really required?  No, we can just keep the
       * shift in a hier::IntVector and adjust.
       */
      for (tbox::Dimension::dir_t d = 0; d < d_common->getDim().getValue(); ++d) {
         VectorOfInts& h = d_histogram[d];
         const int shift = new_lower(d) - d_box.lower() (d);
         if (shift > 0) {
            int i;
            for (i = 0; i < new_size(d); ++i) {
               h[i] = h[i + shift];
            }
         }
         h.resize(new_size(d), BAD_INTEGER);
      }
      if (d_common->d_log_node_history) {
         d_common->writeCounters();
         tbox::plog << "Shrunken " << d_generation << ':' << d_pos
                    << "  " << d_box << " -> " << new_box
                    << ".\n";
      }
      d_box.initialize(new_box, d_box.getLocalId(), d_box.getOwnerRank());
   }
}

/*
 *********************************************************************
 * Accept the box or split it, setting d_box_acceptance accordingly.
 *********************************************************************
 */
void
BergerRigoutsosNode::acceptOrSplitBox()
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_box.getOwnerRank() != d_common->d_mpi.getRank()) {
      TBOX_ERROR("Only the owner can determine\n"
         "whether to accept or split a box.\n");
   }
#endif
   TBOX_ASSERT(d_box_acceptance == undetermined);

   const hier::IntVector boxdims(d_box.numberCells());
   const hier::IntVector oversize(boxdims - d_common->d_max_box_size);

   /*
    * Box d_box is acceptable if
    * - it has a high enough fraction of tagged cells, or
    * - it cannot be split without breaking the minimum
    *   box requirement, or
    *
    * If d_box has no tags:
    * - set d_box_acceptance = hasnotag_by_owner
    * If accepting d_box:
    * - set d_box_acceptance = accepted_by_calculation
    * If rejecting d_box:
    * - set d_box_acceptance = rejected_by_calculation
    * - create left and right children
    * - set children boxes
    * - claim MPI tags for communication by children nodes
    */

   if (d_box_acceptance == undetermined) {
      if (oversize <= hier::IntVector::getZero(d_common->getDim())) {
         /*
          * See if d_box should be accepted based on efficiency.
          */
         int num_tagged = 0;
         for (size_t i = 0; i < d_histogram[0].size(); ++i) {
            num_tagged += d_histogram[0][i];
         }
         size_t boxsize = d_box.size();
         double efficiency = (boxsize == 0 ? 1.e0 :
                              ((double)num_tagged) / static_cast<double>(boxsize));

         if (d_common->d_max_tags_owned < num_tagged) {
            d_common->d_max_tags_owned = num_tagged;
         }

         if (efficiency >= d_common->getEfficiencyTolerance(d_common->d_level_number)) {
            d_box_acceptance = accepted_by_calculation;
            if (d_common->d_log_node_history) {
               d_common->writeCounters();
               tbox::plog << "Accepted " << d_generation << ':' << d_pos
                          << "  " << d_box << " by sufficient efficiency of " << efficiency
                          << "  accept=" << d_box_acceptance << ".\n";
            }
         } else if (num_tagged == 0) {
            // No tags!  This should be caught at the root.
            TBOX_ASSERT(d_parent == 0);
            d_box_acceptance = hasnotag_by_owner;
            if (d_common->d_log_node_history) {
               d_common->writeCounters();
               tbox::plog << "HasNoTag " << d_generation << ':' << d_pos
                          << "  " << d_box
                          << ".\n";
            }
         }
      }
   }

   /*
    * If d_box cannot be split without violating min_size, it should
    * be accepted.
    *
    * If cut_margin is negative in any direction, we cannot cut d_box
    * across that direction without violating min_box.
    */
   hier::IntVector min_size(d_min_box_size);
   min_size.max( d_common->d_min_box_size_from_cutting );
   const hier::IntVector cut_margin = boxdims - min_size * 2;

   if (d_box_acceptance == undetermined) {
      if (cut_margin < hier::IntVector::getZero(d_common->getDim())) {
         d_box_acceptance = accepted_by_calculation;
      } else if (d_box.size() <= d_min_cell_request * 2) {
         d_box_acceptance = accepted_by_calculation;
      }
   }

   hier::IntVector sorted_margins(d_common->getDim());

   if (d_box_acceptance == undetermined) {
      /*
       * Sort the bounding box directions from largest to smallest cut
       * margin.  If there are multiple cuttable directions, we will
       * favor the direction with the greatest cut_margin.
       */
      for (tbox::Dimension::dir_t dim = 0; dim < d_common->getDim().getValue(); ++dim) {
         sorted_margins(dim) = dim;
      }
      for (tbox::Dimension::dir_t d0 = 0; d0 < d_common->getDim().getValue() - 1; ++d0) {
         for (tbox::Dimension::dir_t d1 = static_cast<tbox::Dimension::dir_t>(d0 + 1);
              d1 < d_common->getDim().getValue();
              ++d1) {
            if (cut_margin(sorted_margins(d0)) <
                cut_margin(sorted_margins(d1))) {
               int tmp_dim = sorted_margins(d0);
               sorted_margins(d0) = sorted_margins(d1);
               sorted_margins(d1) = tmp_dim;
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (tbox::Dimension::dir_t dim = 0; dim < d_common->getDim().getValue() - 1; ++dim) {
         TBOX_ASSERT(cut_margin(sorted_margins(dim)) >=
            cut_margin(sorted_margins(dim + 1)));
      }
#endif
   }

   const int max_margin_dir = sorted_margins(0);
   const int min_margin_dir = sorted_margins(d_common->getDim().getValue() - 1);

   int num_cuttable_dim = 0;

   if (d_box_acceptance == undetermined) {
      /*
       * Determine number of coordinate directions that are cuttable
       * according to the cut_margin.
       */
      for (num_cuttable_dim = 0; num_cuttable_dim < d_common->getDim().getValue();
           ++num_cuttable_dim) {
         if (cut_margin(sorted_margins(num_cuttable_dim)) < 0) {
            break;
         }
      }
      TBOX_ASSERT(num_cuttable_dim > 0);   // We already accounted for un-cuttable case before this point.
   }

   if (d_box_acceptance == undetermined) {

      /*
       * Attempt to split box at a zero interior point in the
       * histogram.  Check each cuttable direction, from
       * largest to smallest, until zero point found.
       */

      int cut_lo, cut_hi;
      int cut_pt = -(tbox::MathUtilities<int>::getMax());
      tbox::Dimension::dir_t cut_dir = 0;
      tbox::Dimension::dir_t dir = 0;
      const hier::Index& box_lo(d_box.lower());
      const hier::Index& box_hi(d_box.upper());
      hier::Index lft_hi(box_hi);
      hier::Index rht_lo(box_lo);

      for (dir = 0; dir < d_common->getDim().getValue(); ++dir) {
         cut_dir = static_cast<tbox::Dimension::dir_t>(sorted_margins(dir));
         if (cut_margin(cut_dir) < 0) {
            continue;  // This direction is too small to cut.
         }
         if (findZeroCutSwath(cut_lo, cut_hi, cut_dir)) {
            // Split bound box at cut_pt; cut_dir is splitting direction.
            TBOX_ASSERT(cut_hi - cut_lo >= 0);
            lft_hi(cut_dir) = cut_lo - 1;
            rht_lo(cut_dir) = cut_hi + 1;
            if (d_common->d_log_node_history) {
               d_common->writeCounters();
               tbox::plog << "HoleCut " << d_generation << ':' << d_pos
                          << "  " << d_box << " d=" << cut_dir
                          << " at " << cut_lo << '-' << cut_hi
                          << ".\n";
            }
            break;
         }
      }

      /*
       * If no zero point found, try inflection cut.
       */

      if (dir == d_common->getDim().getValue()) {

         /*
          * inflection_cut_threshold_ar specifies the mininum box
          * thickness that can be cut, as a ratio to the thinnest box
          * direction.  If the box doesn't have any direction thick
          * enough, then it has a reasonable aspect ratio, so we can
          * cut it in any direction.
          *
          * Degenerate values of inflection_cut_threshold_ar:
          *
          * 1: cut any direction except the thinnest.
          *
          * (0,1) and huge values: cut any direction.
          *
          * 0: Not a degenerate case but a special case meaning cut
          * only the thickest direction.  This leads to more cubic
          * boxes but can miss feature edges aligned across other
          * directions.
          *
          * Experiments show that a value of 4 works well.
          */
         int max_box_length_to_leave = boxdims(max_margin_dir) - 1;
         if (d_common->d_inflection_cut_threshold_ar > 0.0) {
            max_box_length_to_leave =
               static_cast<int>(0.5 + boxdims(min_margin_dir)
                                * d_common->d_inflection_cut_threshold_ar);
            if (max_box_length_to_leave >= boxdims(max_margin_dir)) {
               /*
                * Box aspect ratio is not too bad. Disable preference
                * for cutting longer dirs.
                */
               max_box_length_to_leave = 0;
            }
         }

         int inflection = -1;
         for (tbox::Dimension::dir_t d = 0; d < d_common->getDim().getValue(); ++d) {
            if (cut_margin(d) < 0 || boxdims(d) <= max_box_length_to_leave) {
               continue;  // Direction d is too small to cut.
            }
            int try_cut_pt, try_inflection;
            cutAtInflection(try_cut_pt, try_inflection, d);
            if (inflection < try_inflection ||
                (inflection == try_inflection && cut_margin(d) > cut_margin(cut_dir))) {
               cut_dir = d;
               cut_pt = try_cut_pt;
               inflection = try_inflection;
            }
         }
         TBOX_ASSERT(cut_dir < d_common->getDim().getValue());

         // Split bound box at cut_pt; cut_dir is splitting direction.
         lft_hi(cut_dir) = cut_pt - 1;
         rht_lo(cut_dir) = cut_pt;
         if (d_common->d_log_node_history) {
            d_common->writeCounters();
            tbox::plog << "LapCut " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << " d=" << cut_dir << " at " << cut_pt
                       << ".\n";
         }
      }

      /*
       * The owner forms its children now so it can use their
       * parameters while determining which to run.
       * Contributors create the children when they receive
       * the d_box_acceptance flag from the owner.
       */
      d_lft_child = new BergerRigoutsosNode(d_common, this, 0);
      d_rht_child = new BergerRigoutsosNode(d_common, this, 1);

      d_lft_child->d_box = hier::Box(box_lo, lft_hi, d_box.getBlockId());
      d_rht_child->d_box = hier::Box(rht_lo, box_hi, d_box.getBlockId());
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_box.numberCells() >= d_min_box_size) {
         TBOX_ASSERT(d_lft_child->d_box.numberCells() >= d_min_box_size);
         TBOX_ASSERT(d_rht_child->d_box.numberCells() >= d_min_box_size);
      }
#endif

      d_lft_child->claimMPITag();
      d_rht_child->claimMPITag();

      d_box_acceptance = rejected_by_calculation;

      if (d_common->d_log_node_history) {
         d_common->writeCounters();
         tbox::plog << "Lc Split "
                    << d_generation << ':' << d_pos << "  " << d_box
                    << " => " << d_lft_child->d_generation << ':'
                    << d_lft_child->d_pos << d_lft_child->d_box
                    << " + " << d_rht_child->d_generation << ':'
                    << d_rht_child->d_pos << d_rht_child->d_box
                    << ".\n";
      }

   }
}

/*
 ********************************************************************
 *
 * Attempt to find a range with zero histogram value near the
 * middle of d_box in the given coordinate direction.
 * Note that the hole is kept more than a minimium distance from
 * the endpoints of of the index interval.
 *
 * Note that it is assumed that box indices are cell indices.
 *
 * If a hole is found, cut_lo and cut_hi are set to the
 * range of zero tag cells.
 *
 * Optimization note: There seems to be no reason to look for a single
 * zero swath.  If there are multiple zero swaths in the signature,
 * why not cut through them all and produce multiple children?  We may
 * have to change d_lft_child and d_rgt_child to d_children[].  If
 * we don't cut all the signatures we see, we'd just force the
 * children to recompose those signatures themselves.  Making multiple
 * cuts in the box can significantly reduce the amount of data
 * communication in the children nodes.  The downside is that we have
 * to spend more time finding these cuts before we can notify the
 * children.  However, if we don't find the cuts when we have the
 * chance, the children would have to spend their time looking.
 *
 *
 ********************************************************************
 */

bool
BergerRigoutsosNode::findZeroCutSwath(
   int& cut_lo,
   int& cut_hi,
   const tbox::Dimension::dir_t dim)
{
   const int lo = d_box.lower(dim);
   const int hi = d_box.upper(dim);
   // Compute the limit for the swath.
   const int cut_lo_lim = lo + d_min_box_size(dim);
   const int cut_hi_lim = hi - d_min_box_size(dim);

   /*
    * Start in the middle of the box.
    * Move cut_lo down and cut_hi up until a hole is found.
    * Keep moving in same direction of the hole until the
    * other side of the hole is found.  The two planes form
    * the widest cut possible at the hole.
    */
   cut_lo = cut_hi = (lo + hi) / 2;
   while ((cut_lo >= cut_lo_lim) && (cut_hi <= cut_hi_lim)) {
      if (d_histogram[dim][cut_lo - lo] == 0) {
         /* The narrow cut is at cut_lo.  Initialize the cut swath here
          * and move cut_lo down until the far side the hole is found.
          */
         cut_hi = cut_lo;
         while (((cut_lo > cut_lo_lim)) &&
                (d_histogram[dim][cut_lo - lo - 1] == 0)) {
            --cut_lo;
         }
         TBOX_ASSERT(cut_hi >= cut_lo);
         TBOX_ASSERT(cut_lo - lo >= d_min_box_size(dim));
         TBOX_ASSERT(hi - cut_hi >= d_min_box_size(dim));
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int i = cut_lo; i <= cut_hi; ++i) {
            TBOX_ASSERT(d_histogram[dim][i - lo] == 0);
         }
#endif
         return true;
      }
      if (d_histogram[dim][cut_hi - lo] == 0) {
         /* The narrow cut is at cut_hi.  Initialize the cut swath here
          * and move cut_hi up until the far side the hole is found.
          */
         cut_lo = cut_hi;
         while (((cut_hi < cut_hi_lim)) &&
                (d_histogram[dim][cut_hi - lo + 1] == 0)) {
            ++cut_hi;
         }
         TBOX_ASSERT(cut_hi >= cut_lo);
         TBOX_ASSERT(cut_lo - lo >= d_min_box_size(dim));
         TBOX_ASSERT(hi - cut_hi >= d_min_box_size(dim));
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int i = cut_lo; i <= cut_hi; ++i) {
            TBOX_ASSERT(d_histogram[dim][i - lo] == 0);
         }
#endif
         return true;
      }
      --cut_lo;
      ++cut_hi;
   }

   return false;
}

/*
 ***********************************************************************
 *
 * Attempt to find a point in the given coordinate direction near an
 * inflection point in the histogram for that direction. Note that the
 * cut point is kept more than a minimium distance from the endpoints
 * of the index interval (lo, hi).  Also, the box must have at least
 * three cells along a side to apply the inflection test.  If no
 * inflection point is found, the mid-point of the interval is
 * returned as the cut point.
 *
 * Note that it is assumed that box indices are cell indices.
 *
 ***********************************************************************
 */

void
BergerRigoutsosNode::cutAtInflection(
   int& cut_pt,
   int& inflection,
   const tbox::Dimension::dir_t dim)
{
   /*
    * New implementation prefers and possibly restricts the inflection
    * cut to the center part of the box.
    *
    * The cuts refer to face indices, not cell indices.
    *
    * Note that we work in the index space centered on the box's lower
    * cell and add the box lower cell index at the end.
    */

   const VectorOfInts& hist = d_histogram[dim];
   const unsigned int hist_size = static_cast<int>(hist.size());
   TBOX_ASSERT(d_box.upper() (dim) - d_box.lower() (dim) + 1 == static_cast<int>(hist_size));
   TBOX_ASSERT(hist_size >= 2);

   /*
    * Inflection cut requires at least 4 cells of histogram, so it can
    * compute an inflection value.  Without 4 cells, we just cut
    * across the largest change in the histogram.
    */
   if (hist_size < 4) {
      cut_pt = 1;
      for (unsigned int i = 2; i < hist_size; ++i) {
         if (tbox::MathUtilities<int>::Abs(hist[cut_pt] - hist[cut_pt - 1]) <
             tbox::MathUtilities<int>::Abs(hist[i] - hist[i - 1])) {
            cut_pt = i;
         }
      }
      cut_pt += d_box.lower() (dim);
      inflection = 0;  // Not necessarily an inflection point.
      return;
   }

   const int min_box_size =
      tbox::MathUtilities<int>::Max( d_min_box_size(dim),
                                     d_common->d_min_box_size_from_cutting(dim) );

   const int box_lo = 0;
   const int box_hi = hist_size - 1;
   const int max_dist_from_center =
      int(d_common->d_max_inflection_cut_from_center * hist_size / 2);
   const int box_mid = (box_lo + box_hi + 1) / 2;

   const int cut_lo_lim = tbox::MathUtilities<int>::Max(
         box_lo + min_box_size, box_mid - max_dist_from_center);

   const int cut_hi_lim = tbox::MathUtilities<int>::Min(
         box_hi - min_box_size + 1, box_mid + max_dist_from_center);

   /*
    * Initial cut point and differences between the Laplacian on
    * either side of it.  We want to cut where the difference between
    * the two Laplacians is greatest and they have oposite signs.
    */
   cut_pt = box_mid;
   inflection =
      (hist[cut_pt - 1] - 2 * hist[cut_pt] + hist[cut_pt + 1])
      - (hist[cut_pt - 2] - 2 * hist[cut_pt - 1] + hist[cut_pt]);
   inflection = tbox::MathUtilities<int>::Abs(inflection);

   int cut_lo = box_mid - 1;
   int cut_hi = box_mid + 1;

   while (cut_lo > cut_lo_lim || cut_hi < cut_hi_lim) {
      if (cut_lo > cut_lo_lim) {
         const int la = (hist[cut_lo - 1] - 2 * hist[cut_lo] + hist[cut_lo + 1]);
         const int lb = (hist[cut_lo - 2] - 2 * hist[cut_lo - 1] + hist[cut_lo]);
         if (la * lb <= 0) {
            const int try_inflection = tbox::MathUtilities<int>::Abs(la - lb);
            if (try_inflection > inflection) {
               cut_pt = cut_lo;
               inflection = try_inflection;
            }
         }
      }
      if (cut_hi < cut_hi_lim) {
         const int la = (hist[cut_hi - 1] - 2 * hist[cut_hi] + hist[cut_hi + 1]);
         const int lb = (hist[cut_hi - 2] - 2 * hist[cut_hi - 1] + hist[cut_hi]);
         if (la * lb <= 0) {
            const int try_inflection = tbox::MathUtilities<int>::Abs(la - lb);
            if (try_inflection > inflection) {
               cut_pt = cut_hi;
               inflection = try_inflection;
            }
         }
      }
      --cut_lo;
      ++cut_hi;
   }

   cut_pt += d_box.lower() (dim);
}

/*
 ********************************************************************
 * Create a DLBG Box in d_new_box_level,
 * where the output boxes of the algorithm is saved.
 *
 * Only the owner should create the box_level node this way.
 * Other processes build box_level node using data from owner.
 *
 * TODO: this should be renamed putBoxInBoxLevel.
 ********************************************************************
 */
void
BergerRigoutsosNode::createBox()
{
   TBOX_ASSERT(d_common->d_mpi.getRank() == d_box.getOwnerRank());
   hier::LocalId last_index =
      d_common->d_new_box_level->getBoxes().empty() ? hier::LocalId(-1) :
      d_common->d_new_box_level->getBoxes().back().getLocalId();

   hier::Box new_box(d_box, last_index + 1, d_common->d_mpi.getRank());
   TBOX_ASSERT(new_box.getBlockId() == d_box.getBlockId());
   d_common->d_new_box_level->addBoxWithoutUpdate(new_box);
   d_box_iterator = d_common->d_new_box_level->getBox(new_box);

   TBOX_ASSERT(d_box_iterator->isSpatiallyEqual(d_box));
   d_box = *d_box_iterator;
}

/*
 ********************************************************************
 * Discard the Box.  On the owner, this Box is a part of
 * d_new_box_level where it must be removed.  On
 * contributors the Box can just be ignored.  To prevent bugs,
 * the node and its iterator are set to unusable values.
 ********************************************************************
 */
void
BergerRigoutsosNode::eraseBox()
{
   if (d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      d_common->d_new_box_level->eraseBoxWithoutUpdate(
         *d_box_iterator);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   d_box_iterator = hier::BoxContainer().end();
#endif
}

void
BergerRigoutsosNode::countOverlapWithLocalPatches()
{
   /*
    * Count overlaps for the left and right sides.
    *
    * Remove the child if it has zero overlap.
    */
   hier::Box lft_grown_box = d_lft_child->d_box;
   lft_grown_box.grow(d_common->d_tag_to_new_width);
   hier::Box rht_grown_box = d_rht_child->d_box;
   rht_grown_box.grow(d_common->d_tag_to_new_width);
   size_t& lft_overlap = d_lft_child->d_overlap;
   size_t& rht_overlap = d_rht_child->d_overlap;
   lft_overlap = rht_overlap = 0;

   const hier::PatchLevel& tag_level = *d_common->d_tag_level;
   for (hier::PatchLevel::iterator ip(tag_level.begin());
        ip != tag_level.end(); ++ip) {

      const hier::Box& patch_box = (*ip)->getBox();
      const hier::BlockId& block_id = patch_box.getBlockId();

      if (block_id == d_box.getBlockId()) {

         hier::Box lft_intersection = patch_box * lft_grown_box;
         lft_overlap += lft_intersection.size();

         hier::Box rht_intersection = patch_box * rht_grown_box;
         rht_overlap += rht_intersection.size();

      } else {

         hier::Box transform_box(patch_box);
         bool transformed =
            d_common->d_tag_level->getGridGeometry()->transformBox(
               transform_box,
               d_common->d_tag_level->getLevelNumber(),
               d_box.getBlockId(),
               block_id);

         if (transformed) {
            hier::Box lft_intersection = transform_box * lft_grown_box;
            lft_overlap += lft_intersection.size();

            hier::Box rht_intersection = transform_box * rht_grown_box;
            rht_overlap += rht_intersection.size();
         }
      }
   }
}

/*
 *************************************************************************
 * Child groups are subsets of current group.  Each child group
 * includes processes owning patches that overlap the box of that child.
 * The overlap data has been gathered in d_recv_msg.
 * See gatherGroupingCriteria_start() for the format of the message.
 *************************************************************************
 */
void
BergerRigoutsosNode::formChildGroups()
{
   /*
    * Form child groups and determine owners from data gathered
    * in the gather_overlap_counts phase.
    */
   if (d_group.size() == 1) {
      // Short cut for trivial groups.
      d_lft_child->d_group.resize(1, BAD_INTEGER);
      d_rht_child->d_group.resize(1, BAD_INTEGER);
      d_lft_child->d_group[0] = d_group[0];
      d_rht_child->d_group[0] = d_group[0];
      d_lft_child->d_box.initialize(d_lft_child->d_box,
         d_lft_child->d_box.getLocalId(),
         d_box.getOwnerRank());
      d_rht_child->d_box.initialize(d_rht_child->d_box,
         d_rht_child->d_box.getLocalId(),
         d_box.getOwnerRank());
      return;
   }

   d_lft_child->d_group.resize(d_group.size(), BAD_INTEGER);
   d_rht_child->d_group.resize(d_group.size(), BAD_INTEGER);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Only owner process should be here.
    */
   if (d_common->d_mpi.getRank() != d_box.getOwnerRank()) {
      TBOX_ERROR("Library error!" << std::endl);
   }
#endif
   TBOX_ASSERT(d_recv_msg.size() == 4 * d_group.size());

   int* lft_overlap = &d_recv_msg[0];
   int* rht_overlap = &d_recv_msg[1];

   const int imyself = findOwnerInGroup(d_common->d_mpi.getRank(), d_group);
   // TODO: Use MessageStream instead of int* for buffer to pass values greater than int.
   lft_overlap[imyself * 4] = static_cast<int>(d_lft_child->d_overlap);
   rht_overlap[imyself * 4] = static_cast<int>(d_rht_child->d_overlap);

   int* lft_criteria = 0;
   int* rht_criteria = 0;
   switch (d_common->d_owner_mode) {
      case BergerRigoutsos::SINGLE_OWNER:
         lft_criteria = &d_recv_msg[0];
         rht_criteria = &d_recv_msg[1];
         lft_criteria[imyself * 4] = tbox::MathUtilities<int>::getMax();
         rht_criteria[imyself * 4] = tbox::MathUtilities<int>::getMax();
         break;
      case BergerRigoutsos::MOST_OVERLAP:
         lft_criteria = &d_recv_msg[0];
         rht_criteria = &d_recv_msg[1];
         lft_criteria[imyself * 4] = static_cast<int>(d_lft_child->d_overlap);
         rht_criteria[imyself * 4] = static_cast<int>(d_rht_child->d_overlap);
         break;
      case BergerRigoutsos::FEWEST_OWNED:
         lft_criteria = &d_recv_msg[2];
         rht_criteria = &d_recv_msg[2];
         lft_criteria[imyself * 4] = -d_common->d_num_nodes_owned;
         rht_criteria[imyself * 4] = -d_common->d_num_nodes_owned;
         break;
      case BergerRigoutsos::LEAST_ACTIVE:
         lft_criteria = &d_recv_msg[3];
         rht_criteria = &d_recv_msg[3];
         lft_criteria[imyself * 4] = -d_common->d_num_nodes_active;
         rht_criteria[imyself * 4] = -d_common->d_num_nodes_active;
         break;
      default:
         TBOX_ERROR("LIBRARY error" << std::endl);
         break;
   }

   int n_lft = 0;
   int n_rht = 0;

   int lft_owner_score = tbox::MathUtilities<int>::getMin();
   int rht_owner_score = tbox::MathUtilities<int>::getMin();

   /*
    * Loop through the group to see which process should participate
    * on the left/right sides.  Also see which process should be the
    * owner of the left/right sides.  For efficiency in some searches
    * through d_groups, make sure that d_group is ordered.
    */
   int lft_owner = -1;
   int rht_owner = -1;
   for (unsigned int i = 0; i < d_group.size(); ++i) {
      int i4 = i * 4;
      if (lft_overlap[i4] != 0) {
         d_lft_child->d_group[n_lft++] = d_group[i];
         if (lft_criteria[i4] > lft_owner_score) {
            lft_owner = d_group[i];
            lft_owner_score = lft_criteria[i4];
         }
      }
      if (rht_overlap[i4] != 0) {
         d_rht_child->d_group[n_rht++] = d_group[i];
         if (rht_criteria[i4] > rht_owner_score) {
            rht_owner = d_group[i];
            rht_owner_score = rht_criteria[i4];
         }
      }
   }
   d_lft_child->d_box.initialize(d_lft_child->d_box,
      d_lft_child->d_box.getLocalId(),
      lft_owner);
   d_rht_child->d_box.initialize(d_rht_child->d_box,
      d_rht_child->d_box.getLocalId(),
      rht_owner);

   d_lft_child->d_group.resize(n_lft, BAD_INTEGER);
   d_rht_child->d_group.resize(n_rht, BAD_INTEGER);

   // Recall that only the owner should execute this code.
   TBOX_ASSERT(d_lft_child->d_box.getOwnerRank() >= 0);
   TBOX_ASSERT(d_lft_child->d_group.size() > 0);
   TBOX_ASSERT(d_lft_child->d_group.size() <= d_group.size());
   TBOX_ASSERT(d_common->d_owner_mode == BergerRigoutsos::SINGLE_OWNER ||
      ((d_lft_child->d_overlap == 0) !=
       inGroup(d_lft_child->d_group)));
   TBOX_ASSERT(d_rht_child->d_box.getOwnerRank() >= 0);
   TBOX_ASSERT(d_rht_child->d_group.size() > 0);
   TBOX_ASSERT(d_rht_child->d_group.size() <= d_group.size());
   TBOX_ASSERT(d_common->d_owner_mode == BergerRigoutsos::SINGLE_OWNER ||
      ((d_rht_child->d_overlap == 0) !=
       inGroup(d_rht_child->d_group)));
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_common->d_owner_mode == BergerRigoutsos::SINGLE_OWNER) {
      TBOX_ASSERT(inGroup(d_lft_child->d_group, d_box.getOwnerRank()));
      TBOX_ASSERT(inGroup(d_rht_child->d_group, d_box.getOwnerRank()));
   }
   for (size_t i = 0; i < d_group.size(); ++i) {
      TBOX_ASSERT(i == 0 || d_group[i] > d_group[i - 1]);
      TBOX_ASSERT((lft_overlap[i * 4] > 0 ||
                   (d_group[i] == d_lft_child->d_box.getOwnerRank()))
         == inGroup(d_lft_child->d_group, d_group[i]));
      TBOX_ASSERT((rht_overlap[i * 4] > 0 ||
                   (d_group[i] == d_rht_child->d_box.getOwnerRank()))
         == inGroup(d_rht_child->d_group, d_group[i]));
   }
#endif
}

/*
 *************************************************************************
 *
 * Compute overlaps between the new graph node and nodes on
 * the tagged level, saving that data in the form of relationships.
 *
 * Note that the relationship data may be duplicated in two objects.
 * - tag_to_new stores the relationships organized around each node
 *   in the tagged level.  For each node on the tagged level,
 *   we store a container of neighbors on the new box_level.
 * - new_to_tag stores the relationships organized around each NEW node.
 *   For each new node we store a container of neighbors on the
 *   tagged level.
 *
 * If compute_relationships > 0, we store tag_to_new.
 *
 * If compute_relationships > 1, we also compute new_to_tag.
 * The data in new_to_tag are
 * computed by the particant processes but eventually stored on the
 * owners of the new nodes, so their computation requires caching
 * the relationship data in relationship_messages for sending to the appropriate
 * processes later.
 *
 *************************************************************************
 */
void
BergerRigoutsosNode::computeNewNeighborhoodSets()
{
   d_common->d_object_timers->t_compute_new_neighborhood_sets->start();
   TBOX_ASSERT(d_common->d_compute_relationships > 0);
   TBOX_ASSERT(boxAccepted());
   TBOX_ASSERT(d_box_acceptance != accepted_by_dropout_bcast);
   /*
    * Do not check for min_box violation in root node.  That
    * check should be done outside of this class in order to
    * have flexibility regarding how to handle it.
    */
   TBOX_ASSERT(d_parent == 0 || d_box.numberCells() >= d_min_box_size);
   /*
    * We should not compute nabrs if we got the node
    * by a dropout broadcast because we already know
    * there is no overlap!
    */
   TBOX_ASSERT(d_box_acceptance != accepted_by_dropout_bcast);

   // Create an expanded box for intersection check.
   hier::Box grown_box = d_box;
   grown_box.grow(d_common->d_tag_to_new_width);
   hier::BoxContainer grown_boxes;
   if (d_common->d_tag_level->getGridGeometry()->getNumberBlocks() == 1 ||
       d_common->d_tag_level->getGridGeometry()->hasIsotropicRatios()) {
      grown_boxes.pushBack(d_box);
      grown_boxes.grow(d_common->d_tag_to_new_width);
   } else {
      hier::BoxUtilities::growAndAdjustAcrossBlockBoundary(
         grown_boxes,
         d_box,
         d_common->d_tag_level->getGridGeometry(),
         d_common->d_tag_to_new->getBase().getRefinementRatio(),
         d_common->d_tag_to_new->getRatio(),
         d_common->d_tag_to_new_width,
         false,
         false);
   }

   /*
    * On the owner process, we store the neighbors of the new node.
    * This data is NOT required on other processes.
    */
   bool on_owner_process = d_common->d_mpi.getRank() == d_box.getOwnerRank();
   if (on_owner_process) {
      d_common->d_tag_to_new->getTranspose().makeEmptyLocalNeighborhood(d_box.getBoxId());
   }

   // Data to send to owner regarding new relationships found by local process.
   VectorOfInts* relationship_message = 0;
   if (d_common->d_compute_relationships > 1 && d_common->d_mpi.getRank() !=
       d_box.getOwnerRank()) {
      /*
       * Will have to send to owner the relationships found locally for
       * d_box.
       * Label the id of the new node and the (yet unknown) number
       * of relationship found for it.
       *
       * The message to be sent to owner is appended the following
       * data:
       * - index of new node
       * - number of relationships found for the new node
       * - index of nodes on the tagged level overlapping new node.
       */
      relationship_message = &d_common->d_relationship_messages[d_box.getOwnerRank()];
      relationship_message->insert(relationship_message->end(), d_box.getLocalId().getValue());
      relationship_message->insert(relationship_message->end(), 0);
   }

   const int index_of_counter =
      (relationship_message != 0 ? static_cast<int>(relationship_message->size()) : 0) - 1;
   const int ints_per_node = hier::Box::commBufferSize(d_common->getDim());

   const hier::BoxContainer& tag_boxes = d_common->d_tag_level->getBoxLevel()->getBoxes();

   for (hier::RealBoxConstIterator ni(tag_boxes.realBegin());
        ni != tag_boxes.realEnd(); ++ni) {

      const hier::Box& tag_box = *ni;

      bool intersection = false;
      for (hier::BoxContainer::const_iterator b_itr = grown_boxes.begin();
           !intersection && b_itr != grown_boxes.end(); ++b_itr) {
         if (tag_box.getBlockId() == b_itr->getBlockId()) {
            intersection = tag_box.intersects(*b_itr);
         } else {
            hier::Box transform_box(tag_box);
            bool transformed =
               d_common->d_tag_level->getGridGeometry()->transformBox(
                  transform_box,
                  d_common->d_tag_level->getRatioToLevelZero(),
                  b_itr->getBlockId(),
                  tag_box.getBlockId());
            if (transformed) {
               intersection = transform_box.intersects(*b_itr);
            }
         }
      }

      if (intersection) {

         // Add d_box as a neighbor of tag_box.
         d_common->d_tag_to_new->insertLocalNeighbor(d_box,
            tag_box.getBoxId());

         if (on_owner_process) {
            // Owner adds tag_box as a neighbor of d_box.
            d_common->d_tag_to_new->getTranspose().insertLocalNeighbor(tag_box,
               d_box.getBoxId());
         }

         if (relationship_message != 0) {
            /* Non-owners put found relationship in the message
             * to (eventually) send to owner.
             */
            relationship_message->insert(relationship_message->end(), ints_per_node, 0);
            int* ptr = &(*relationship_message)[relationship_message->size() - ints_per_node];
            tag_box.putToIntBuffer(ptr);
            ++(*relationship_message)[index_of_counter];
         }
      }
   }

   if (d_common->d_compute_relationships > 1 &&
       d_common->d_mpi.getRank() == d_box.getOwnerRank()) {
      /*
       * If box was accepted, the owner should remember
       * which process will be sending relationship data.
       * Update the list of relationship senders to make sure
       * it includes all processes in the group.
       * We use this list in
       * BergerRigoutsos::shareNewNeighborhoodSetsWithOwners to
       * tell us which processors are sending us new relationships.
       * The relationship senders are the participants of the group.
       */
      d_common->d_relationship_senders.insert(d_group.begin(), d_group.end());
   }

   d_common->d_object_timers->t_compute_new_neighborhood_sets->stop();
}

/*
 ********************************************************************
 * Utility methods.
 ********************************************************************
 */

int *
BergerRigoutsosNode::putHistogramToBuffer(
   int* buffer)
{
   int dim_val = d_common->getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dim_val; ++d) {
      d_histogram[d].resize(d_box.numberCells(d), BAD_INTEGER);
      memcpy(buffer,
         &d_histogram[d][0],
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

int *
BergerRigoutsosNode::getHistogramFromBuffer(
   int* buffer)
{
   tbox::Dimension::dir_t dim_val = d_common->getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dim_val; ++d) {
      TBOX_ASSERT((int)d_histogram[d].size() == d_box.numberCells(d));
      // d_histogram[d].resizeArray( d_box.numberCells(d) );
      memcpy(&d_histogram[d][0],
         buffer,
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

int *
BergerRigoutsosNode::putBoxToBuffer(
   const hier::Box& box,
   int* buffer) const
{
   const hier::Index& l = box.lower();
   const hier::Index& u = box.upper();
   int dim_val = d_common->getDim().getValue();
   for (int d = 0; d < dim_val; ++d) {
      *(buffer++) = l(d);
      *(buffer++) = u(d);
   }
   return buffer;
}

int *
BergerRigoutsosNode::getBoxFromBuffer(
   hier::Box& box,
   int* buffer) const
{
   int dim_val = d_common->getDim().getValue();
   for (int d = 0; d < dim_val; ++d) {
      box.setLower(static_cast<hier::Box::dir_t>(d), *(buffer++));
      box.setUpper(static_cast<hier::Box::dir_t>(d), *(buffer++));
   }
   return buffer;
}

/*
 ***********************************************************************
 * Put in dropouts things that are in main_group but
 * not in sub_group.
 *
 * Assume that sub_group is a subset of elements in main_group.
 * Assume that sub_group and main_group are sorted in ascending order.
 *
 * Assume add_root is NOT in the dropout and add it anyway.
 ***********************************************************************
 */
void
BergerRigoutsosNode::computeDropoutGroup(
   const VectorOfInts& main_group,
   const VectorOfInts& sub_group,
   VectorOfInts& dropout_group,
   int add_root) const
{
   TBOX_ASSERT(main_group.size() >= sub_group.size());

   dropout_group.resize(main_group.size(), BAD_INTEGER);

   size_t i, j, k = 0;
   dropout_group[k++] = add_root;
   for (i = 0, j = 0; i < main_group.size(); ++i) {
      if (main_group[i] != sub_group[j]) {
         dropout_group[k++] = main_group[i];
      } else {
         ++j;
         if (j == sub_group.size()) {
            // No more in the sub_group so the rest of main_group
            // goes in dropout_group.
            for (i = i + 1; i < main_group.size(); ++i, ++k) {
               dropout_group[k] = main_group[i];
            }
         }
      }
   }

   TBOX_ASSERT(j = sub_group.size());
   dropout_group.resize(k, BAD_INTEGER);
}

/*
 **********************************************************************
 * Claim a unique tag from the processor's available tag pool.
 * Check that the pool is not overused.
 **********************************************************************
 */
void
BergerRigoutsosNode::claimMPITag()
{
   /*
    * Each node should claim no more than one MPI tag
    * so make sure it does not already have one.
    */
   TBOX_ASSERT(d_mpi_tag < 0);

   d_mpi_tag = d_common->d_available_mpi_tag;
   d_common->d_available_mpi_tag = d_mpi_tag + total_phase_tags;
   if (d_mpi_tag + total_phase_tags - 1 >
       d_common->d_tag_upper_bound / (d_common->d_mpi.getSize())
       * (d_common->d_mpi.getRank() + 1)) {
      /*
       * Each process is alloted tag_upper_bound/(d_common->d_mpi.getSize())
       * tag values.  If it needs more than this, it will encroach
       * on the tag pool of the next process and may lead to using
       * non-unique tags.
       */
      TBOX_ERROR("Out of MPI tag values need to ensure that\n"
         << "messages are properly differentiated."
         << "\nd_mpi_tag = " << d_mpi_tag
         << "\ntag_upper_bound = " << d_common->d_tag_upper_bound
         << "\nmber of nodes = " << d_common->d_mpi.getSize()
         << "\nmax tag required = " << d_mpi_tag + total_phase_tags - 1
         << "\nmax tag available = "
         << d_common->d_tag_upper_bound / (d_common->d_mpi.getSize())
         * (d_common->d_mpi.getRank() + 1)
         << std::endl);
      /*
       * It is probably safe to recycle tags if we run out of MPI tags.
       * This is not implemented because thus far, there is no need for it.
       * Recycling is starting over from the initial tag set aside for the
       * local process.  To make sure that recycled tags are not still
       * in use, we should claim a new (or recycled) tag for the dropout
       * broadcast phase.  This is because descendant nodes may recycle
       * the current claimed tag before this phase starts.  All other
       * phases are not interupted by descendant communications, so we
       * are assured that their tag is not doubly claimed.
       */
   }
}

/*
 **********************************************************************
 * Convert an integer value to BoxAcceptance.
 * This is needed because the compiler cannot
 * cast an integer to an enum type.
 **********************************************************************
 */
BergerRigoutsosNode::BoxAcceptance
BergerRigoutsosNode::intToBoxAcceptance(
   int i) const
{
   switch (i) {
      case undetermined: return undetermined;

      case hasnotag_by_owner: return hasnotag_by_owner;

      case rejected_by_calculation: return rejected_by_calculation;

      case accepted_by_calculation: return accepted_by_calculation;

      case rejected_by_owner: return rejected_by_owner;

      case accepted_by_owner: return accepted_by_owner;

      case rejected_by_recombination: return rejected_by_recombination;

      case accepted_by_recombination: return accepted_by_recombination;

      case rejected_by_dropout_bcast: return rejected_by_dropout_bcast;

      case accepted_by_dropout_bcast: return accepted_by_dropout_bcast;

      default:
         TBOX_ERROR("Library error: bad BoxAcceptance data of " << i << ".\n");
   }
   return undetermined;
}

/*
 **********************************************************************
 **********************************************************************
 */
void
BergerRigoutsosNode::printNodeState(
   std::ostream& co) const
{
   co << d_generation << ':' << d_pos << '=' << d_box
      << "  o=" << d_box.getOwnerRank() << ',' << (d_common->d_mpi.getRank() == d_box.getOwnerRank())
      << "  a=" << d_box_acceptance
      << "  w=" << d_wait_phase << '/' << bool(d_comm_group)
      << (d_comm_group ? d_comm_group->isDone() : true)
      << "  t=" << d_num_tags;
   if (d_lft_child) {
      co << "  l=" << d_lft_child->d_generation << ':' << d_lft_child->d_pos
         << '=' << d_lft_child->d_box;
   }
   if (d_rht_child) {
      co << "  r=" << d_rht_child->d_generation << ':' << d_rht_child->d_pos
         << '=' << d_rht_child->d_box;
   }
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
