/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerTree_C
#define included_mesh_CascadePartitionerTree_C

#include "SAMRAI/mesh/CascadePartitionerTree.h"
#include "SAMRAI/mesh/CascadePartitioner.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 * Construct the tree for the given CascadePartitioner by constructing
 * the root and recursively constructing its children.
 *
 * Only children relevant to local process is actually constructed.
 * These are the groups containing the local processes and their
 * sibling groups.
 */
CascadePartitionerTree::CascadePartitionerTree(
   const CascadePartitioner& partitioner):
   d_common(&partitioner),
   d_gen_num(0),

   d_begin(0),
   d_end(partitioner.d_mpi.getSize()),

   d_parent(0),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(partitioner.d_local_load->getSumLoad()),
   d_obligation(d_common->d_global_work_avg * d_common->d_mpi.getSize()),
   d_group_may_supply(false)
{
   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::root constructor: entered generation "
                 << d_gen_num << "  ranks " << d_begin << '-' << d_end << std::endl;
   }

   makeChildren();

   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::root constructor: leaving" << std::endl;
      printClassData(tbox::plog, "\t");
   }
}

/*
 * Construct a child group.  Child group has either the lower half or
 * the upper half of the parent group, indicated by group_position.
 * If parent group has odd number of processes, the extra process is
 * placed in the upper group.
 *
 * Assign contacts by pairing this process with the process having the
 * same relative rank in the sibling group.  If upper group has an
 * extra rank, pair it with the last rank in lower group (giving last
 * rank in lower group two contacts).
 */
CascadePartitionerTree::CascadePartitionerTree(
   CascadePartitionerTree& parent,
   Position group_position):
   d_common(parent.d_common),
   d_gen_num(1 + parent.d_gen_num),

   d_begin(parent.d_begin),
   d_end(parent.d_end),

   d_parent(&parent),
   d_near(0),
   d_far(0),
   d_leaf(0),

   d_work(-1.0),
   d_obligation(-1.0),
   d_group_may_supply(false)
{
   if (d_common->d_print_steps && d_common->d_print_child_steps) {
      tbox::plog << d_common->d_object_name << "::non-root constructor: entered generation "
                 << d_gen_num << "  parent ranks " << d_begin << '-' << d_end
                 << "  position " << group_position << std::endl;
   }

   d_children[0] = d_children[1] = 0;
   d_contact[0] = d_contact[1] = -1;
   d_process_may_supply[0] = d_process_may_supply[1] = false;

   const int upper_begin = (d_parent->d_begin + d_parent->d_end) / 2;

   if (group_position == Lower) {
      d_end = upper_begin;
      const int relative_rank = d_common->d_mpi.getRank() - d_parent->d_begin;
      d_contact[0] = relative_rank + upper_begin;

      if ((d_parent->d_end - d_parent->d_begin) % 2 &&
          d_common->d_mpi.getRank() == upper_begin - 1) {
         d_contact[1] = 1 + d_contact[0];
      }
   } else {
      d_begin = upper_begin;
      const int relative_rank = d_common->d_mpi.getRank() - upper_begin;
      d_contact[0] = tbox::MathUtilities<int>::Min(
            relative_rank + d_parent->d_begin, upper_begin - 1);
   }

   d_obligation = d_common->d_global_work_avg * (d_end - d_begin);

   if (containsRank(d_common->d_mpi.getRank())) {
      makeChildren();
   }

   if (d_common->d_print_steps && d_common->d_print_child_steps) {
      tbox::plog << d_common->d_object_name << "::non-root constructor: leaving" << std::endl;
      printClassData(tbox::plog, "\t");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::makeChildren()
{
   if (size() > 1) {

      d_children[0] = new CascadePartitionerTree(*this, Lower);
      d_children[1] = new CascadePartitionerTree(*this, Upper);

      const bool in_upper_branch = d_children[1]->containsRank(d_common->d_mpi.getRank());
      d_near = d_children[in_upper_branch];
      d_far = d_children[!in_upper_branch];

      d_leaf = d_near->d_leaf;
   } else {
      if (d_begin == d_common->d_mpi.getRank()) {
         d_leaf = this;
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
CascadePartitionerTree::~CascadePartitionerTree()
{
   if (d_children[0]) delete d_children[0];
   if (d_children[1]) delete d_children[1];
   d_children[0] = 0;
   d_children[1] = 0;
   d_near = 0;
   d_far = 0;
   d_leaf = 0;
   d_common = 0;
}

/*
 *************************************************************************
 * This method contains the looping structure required to balance all
 * groups.  The outer loop (top_group) cycles from the largest group
 * to the smallest, balancing the two halves of each group.  The inner
 * loop (current_group) cycles from the leaf to the current top_group,
 * combining groups to build a sufficient global picture of top_group,
 * then balance its two halves.
 *************************************************************************
 */
void CascadePartitionerTree::distributeLoad()
{
   d_common->t_distribute_load->start();

   const double connector_update_interval = computeConnectorUpdateInterval();

   const int tree_depth = CascadePartitioner::lgInt(d_common->d_mpi.getSize());

   CascadePartitionerTree* top_group = this;

   for (int top_gen_number = 0; top_gen_number < tree_depth; ++top_gen_number) {

      /*
       * This block combines and balances child branches.  Each time through:
       * - Compute group weight by combining descendant group weights.
       * - Reset obligations of processes based on group weight.
       * - Balance the two children group.
       *
       * For certain top groups, update Connectors.  All nodes must
       * participate in updating Connectors (or the step will hang
       * when diagnostic barriers are used).
       */
      TBOX_ASSERT(top_group->d_gen_num == top_gen_number);

      if (top_group != d_leaf) {

         if (d_common->d_print_steps) {
            tbox::plog << d_common->d_object_name << "::distributeLoad balancing outer top_group "
                       << top_group->d_gen_num
                       << "  with exact local_load=" << d_common->d_local_load->getSumLoad()
                       << std::endl;
            tbox::plog << "\ttop_group:" << std::endl;
            top_group->printClassData(tbox::plog, "\t");
            tbox::plog << "\tchild 0:" << std::endl;
            top_group->d_children[0]->printClassData(tbox::plog, "\t");
            tbox::plog << "\tchild 1:" << std::endl;
            top_group->d_children[1]->printClassData(tbox::plog, "\t");
         }

         d_leaf->recomputeLeafData();

         /*
          * Loop from leaf's parent toward current_group, combining
          * group weights.
          */
         for ( CascadePartitionerTree *current_group = d_leaf->d_parent;
               current_group != 0 && current_group->d_near != top_group;
               current_group = current_group->d_parent ) {

            current_group->combineChildren();
            if (d_common->d_print_steps && d_common->d_print_child_steps) {
               tbox::plog << d_common->d_object_name << "::distributeLoad outer top_group "
                          << top_group->d_gen_num << "  combined generation "
                          << current_group->d_gen_num << ".  All d_work values are exact."
                          << std::endl;
               tbox::plog << "\tcurrent_group:" << std::endl;
               current_group->printClassData(tbox::plog, "\t");
               tbox::plog << "\tchild 0:" << std::endl;
               current_group->d_children[0]->printClassData(tbox::plog, "\t");
               tbox::plog << "\tchild 1:" << std::endl;
               current_group->d_children[1]->printClassData(tbox::plog, "\t");
            }

         } // Inner loop, current_group

         /*
          * Non-root groups may have average weights different than the
          * global average.  Reset the obligation to the group average
          * if option is on.
          */
         if ( d_common->d_reset_obligations && top_group->d_gen_num != 0 ) {
            const double old_obligation = top_group->d_obligation;
            top_group->resetObligation( top_group->d_work/static_cast<double>(top_group->size()) );
            if ( d_common->d_print_steps ) {
               tbox::plog << d_common->d_object_name << "::distributeLoad generation "
                          << top_group->d_gen_num << " reset obligation from "
                          << old_obligation << " to " << top_group->d_obligation
                          << std::endl;
            }
         }

         /*
          * Balance between children of the top_group.
          */
         top_group->balanceChildren();
         if ( d_common->d_print_steps ) {
            tbox::plog << d_common->d_object_name << "::distributeLoad outer top_group "
                       << top_group->d_gen_num << "  shuffled generation "
                       << top_group->d_gen_num << ".  d_work is exact, but childrens' are estimates."
                       << std::endl;
            tbox::plog << "\ttop_group:" << std::endl;
            top_group->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 0:" << std::endl;
            top_group->d_children[0]->printClassData( tbox::plog, "\t" );
            tbox::plog << "\tchild 1:" << std::endl;
            top_group->d_children[1]->printClassData( tbox::plog, "\t" );
         }

         if ( d_common->d_print_steps ) {
            tbox::plog << d_common->d_object_name << "::distributeLoad completed inner loop for generation "
                       << top_group->d_gen_num << std::endl;
         }

      }

      /*
       * Update Connectors at appropriate intervals or if this is the
       * last time through the loop.
       */
      if ( static_cast<int>(top_gen_number/connector_update_interval) !=
           static_cast<int>((top_gen_number+1)/connector_update_interval) ||
           (top_gen_number == tree_depth-2 && !d_common->d_pparams->usingVouchers()) ) {
         if ( d_common->d_print_steps ) {
            tbox::plog << d_common->d_object_name << "::distributeLoad updating Connectors after balancing generation "
                       << top_group->d_gen_num << std::endl;
         }
         d_common->t_distribute_load->stop();
         d_common->updateConnectors();
         d_common->t_distribute_load->start();
         if (top_group != d_leaf->d_parent) {
            d_common->d_local_load->insertAllWithExistingLoads(
               d_common->d_balance_box_level->getBoxes());
         }
      }

      top_group = top_group->d_near;

   } // Outer loop, top_group

   d_common->t_distribute_load->stop();
}

/*
 *************************************************************************
 * Combine near and far children data (using communication) to compute
 * work-related data for this group.  Update work-related values in
 * this group and its far child.
 *
 * After combining children, d_work, d_near->d_work and d_far->d_work
 * are exact.  d_near->d_work and d_far->d_work remain exact until
 * balanceChildren(), when we just estimate transfers by remote ranks
 * in this group.  d_work remains exact longer--until the next bigger
 * group calls balanceChildren().
 *
 * NOTE: This method can probably be re-organized as if there's only
 * one contact.  We only need to exchange data with the second contact
 * if the the far group may supply and the near group has a deficit.
 *************************************************************************
 */
void
CascadePartitionerTree::combineChildren()
{
   d_common->t_combine_children->start();

   tbox::MessageStream send_msg;
   send_msg << d_near->d_work << d_near->d_group_may_supply << d_near->d_process_may_supply[0];

   for (int i = 0; i < 2; ++i) {
      if (d_near->d_contact[i] >= 0) {

         d_common->d_comm_peer[i].setPeerRank(d_near->d_contact[i]);
         d_common->d_comm_peer[i].setMPITag(CascadePartitionerTree_TAG_InfoExchange0,
            CascadePartitionerTree_TAG_InfoExchange1);
         d_common->d_comm_peer[i].limitFirstDataLength(send_msg.getCurrentSize());
         d_common->d_comm_peer[i].beginRecv(true);

         d_common->d_comm_peer[2 + i].setPeerRank(d_near->d_contact[i]);
         d_common->d_comm_peer[2 + i].setMPITag(CascadePartitionerTree_TAG_InfoExchange0,
            CascadePartitionerTree_TAG_InfoExchange1);
         d_common->d_comm_peer[2 + i].limitFirstDataLength(send_msg.getCurrentSize());
         d_common->d_comm_peer[2 + i].beginSend(static_cast<const char *>(send_msg.getBufferStart()),
            send_msg.getCurrentSize(), true);

      }
   }

   d_far->d_work = 0.0;
   d_far->d_group_may_supply = true;
   while (d_common->d_comm_stage.numberOfCompletedMembers() > 0 ||
          d_common->d_comm_stage.advanceAny()) {

      tbox::AsyncCommPeer<char>* completed = static_cast<tbox::AsyncCommPeer<char> *>(
            d_common->d_comm_stage.popCompletionQueue());

      const int i = static_cast<int>(completed - d_common->d_comm_peer);
      TBOX_ASSERT(i >= 0 && i < 4);
      if (i < 2) {
         // This was a receive.
         tbox::MessageStream recv_msg(completed->getRecvSize(),
                                      tbox::MessageStream::Read,
                                      completed->getRecvData(),
                                      false);
         recv_msg >> d_far->d_work >> d_far->d_group_may_supply >> d_far->d_process_may_supply[i];
      }
   }
   TBOX_ASSERT(d_common->d_comm_stage.numberOfPendingMembers() == 0);

   d_work = d_children[0]->d_work + d_children[1]->d_work;
   d_group_may_supply = estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol();

   // If process still may supply for near child, it may supply for this group.
   d_process_may_supply[0] = d_near->d_process_may_supply[0];

   d_common->t_combine_children->stop();
}

/*
 *************************************************************************
 * If one child has a positive surplus and the other has a negative
 * surplus, the former supplies work to the latter.  Amount supplied
 * is ideally the minimum of the supplier's surplus and the
 * requestor's deficit.  (Actual ammounts are affected by load cutting
 * restrictions.)  Note that only children groups will be balanced
 * (not all descendents).
 *
 * This method records estimates of the work changes to the groups it
 * knows about.  It doesn't record the actual work changes because
 * that happens remotely.  Each process in the supply group send a
 * message received by its contact(s) on the requesting group.  The
 * messages has the actual work to be transfered.
 *
 * Note: This method balances only the children groups (the lower
 * child and upper), not all descendents.  To balance the near grand
 * child, call d_near->balanceChildren().
 *************************************************************************
 */
void
CascadePartitionerTree::balanceChildren()
{
   d_common->t_balance_children->start();

   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::balanceChildren: entered" << std::endl;
   }

   TBOX_ASSERT(d_common->d_shipment->empty());

   if (d_near->estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol() &&
       d_far->estimatedSurplus() < -d_common->d_pparams->getLoadComparisonTol()) {
      // Outgoing work, from near child to far child.

      if (d_near->d_process_may_supply[0]) {

         d_common->t_supply_work->start();
         double work_supplied = d_near->supplyWork(-d_far->estimatedSurplus(), d_near->d_contact[0]);
         d_common->t_supply_work->stop();

         // Record work taken by the far child.
         d_far->d_work += work_supplied;
         d_far->d_group_may_supply = d_far->d_process_may_supply[0] =
               d_far->d_process_may_supply[1] = false;

         if (d_common->d_print_steps) {
            tbox::plog << d_common->d_object_name << "::balanceChildren:"
                       << "  record outgoing shipment of " << work_supplied
                       << " from our half to far half.  Send to " << d_near->d_contact[0]
                       << std::endl;
         }

         TBOX_ASSERT(d_near->d_contact[0] >= 0);
         sendShipment(d_near->d_contact[0]); // If 2 contacts in far group, send to the first one only.
      }
   } else if (d_far->estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol() &&
              d_near->estimatedSurplus() < -d_common->d_pparams->getLoadComparisonTol()) {
      // Incoming work, from far child to near child.

      /*
       * Even when a group may supply, it will not supply to a
       * redundant_demand.  It will supply to its first contact
       * everything it wants to send to the group.  A redundant_demand
       * is one from the second contact in the same group.  Local
       * demand is redundant if the near group is bigger than the far
       * group and local process is the last rank in its group.
       */
      const bool redundant_demand = d_near->size() > d_far->size() &&
         d_common->d_mpi.getRank() == d_end - 1;

      if (d_far->d_group_may_supply && !redundant_demand) {
         if (d_far->d_process_may_supply[0]) {
            d_common->d_comm_peer[0].setPeerRank(d_near->d_contact[0]);
            d_common->d_comm_peer[0].setMPITag(CascadePartitionerTree_TAG_LoadTransfer0,
               CascadePartitionerTree_TAG_LoadTransfer1);
            d_common->d_comm_peer[0].beginRecv(true);
            if (d_common->d_print_steps) {
               tbox::plog << d_common->d_object_name << "::balanceChildren:"
                          << "  expecting shipment from first contact " << d_near->d_contact[0]
                          << "  redundant_demand=" << redundant_demand
                          << "  d_near->size()=" << d_near->size()
                          << "  d_far->size()=" << d_far->size()
                          << "  d_end=" << d_end
                          << "  d_common->d_mpi.getRank()=" << d_common->d_mpi.getRank()
                          << std::endl;
            }
         }
         if (d_far->d_process_may_supply[1]) {
            d_common->d_comm_peer[1].setPeerRank(d_near->d_contact[1]);
            d_common->d_comm_peer[1].setMPITag(CascadePartitionerTree_TAG_LoadTransfer0,
               CascadePartitionerTree_TAG_LoadTransfer1);
            d_common->d_comm_peer[1].beginRecv(true);
            if (d_common->d_print_steps) {
               tbox::plog << d_common->d_object_name << "::balanceChildren:"
                          << "  expecting shipment from second contact " << d_near->d_contact[1]
                          << "  redundant_demand=" << redundant_demand
                          << "  d_near->size()=" << d_near->size()
                          << "  d_far->size()=" << d_far->size()
                          << "  d_end=" << d_end
                          << "  d_common->d_mpi.getRank()=" << d_common->d_mpi.getRank()
                          << std::endl;
            }
         }
      }

      d_common->t_supply_work->start();
      double work_supplied = d_far->supplyWork(-d_near->estimatedSurplus(), d_common->d_mpi.getRank());
      d_common->t_supply_work->stop();

      // Record work taken by near child group.
      d_near->d_work += work_supplied;
      d_near->d_group_may_supply = d_near->d_process_may_supply[0] = false;

      if (d_common->d_print_steps) {
         tbox::plog << d_common->d_object_name << "::balanceChildren:"
                    << "  recorded incoming shipment of " << work_supplied
                    << " from far half to our half." << std::endl;
      }

      if (d_far->d_process_may_supply[0] || d_far->d_process_may_supply[1]) {
         receiveAndUnpackSuppliedLoad();
      }
   } else {
      if (d_common->d_print_steps) {
         tbox::plog << d_common->d_object_name
                    << "::balanceChildren: not supplying or demanding" << std::endl;
      }
   }

   // Complete the load send, if there was any.
   d_common->d_comm_stage.advanceAll();
   while (d_common->d_comm_stage.numberOfCompletedMembers() > 0) {
      d_common->d_comm_stage.popCompletionQueue();
   }

   d_common->t_balance_children->stop();

   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::balanceChildren: leaving" << std::endl;
   }
}

/*
 *************************************************************************
 * Supply work_requested to another group requesting work.  Any load
 * supplied by local process will be sent to the designated taker, a
 * process in the requesting group.  Give priority to supplies closest
 * to the taker in rank space.
 *
 * 1. If this is a leaf, remove load from d_common->d_local_load
 *    and put it in d_common->d_shipment.
 * 2. Else:
 *    A: Remove load from the child group closer to the taker.
 *    B: If step A didn't supply enough, remove some from the other child.
 *
 * This method recurses into descendent groups, estimating changes
 * needed to supply the work_requested.  Estimation is necessary
 * because most of the changes take place remotely.  If the recursion
 * reaches the leaf group containing the local process, it sets aside
 * the work that the local process supplies.
 *************************************************************************
 */
double
CascadePartitionerTree::supplyWork(double work_requested, int taker)
{
   TBOX_ASSERT(work_requested > 0.0);
   TBOX_ASSERT(!containsRank(taker));
   TBOX_ASSERT(containsRank(d_common->d_mpi.getRank()) || d_children[0] == 0);   // Only near groups should store children.
   TBOX_ASSERT(d_group_may_supply ==
      (estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol()));

   if (d_common->d_print_steps && d_common->d_print_child_steps) {
      tbox::plog << d_common->d_object_name << "::supplyWork generation "
                 << d_gen_num << " [" << d_begin << ',' << d_end << ')'
                 << " attempting to supply " << work_requested << " to " << taker
                 << std::endl;
   }

   double est_work_supplied = 0.0; // Estimate of work supplied by this group.

   if (d_group_may_supply) {

      const double allowed_supply = tbox::MathUtilities<double>::Min(
            work_requested, d_common->d_limit_supply_to_surplus ? estimatedSurplus() : d_work);
      TBOX_ASSERT(allowed_supply >= 0.0);

      if (d_children[0] != 0) {
         // This is a near group but not a leaf: Recursively supply load from children.
         const int priority = taker >= d_begin;
         est_work_supplied = d_children[priority]->supplyWork(allowed_supply, taker);
         if (est_work_supplied < allowed_supply) {
            est_work_supplied +=
               d_children[!priority]->supplyWork(allowed_supply - est_work_supplied, taker);
         }
      } else {
         // This is a leaf and/or a far group.  No children, and no recursion.
         est_work_supplied = allowed_supply;

         if (containsRank(d_common->d_mpi.getRank())) {
            // This is a near leaf group: apportion the load shipment.
            TBOX_ASSERT(size() == 1);
            const double tolerance = d_common->d_flexible_load_tol * d_common->d_global_work_avg;
            d_common->d_shipment->adjustLoad(
               *d_common->d_local_load,
               est_work_supplied,
               est_work_supplied - tolerance,
               est_work_supplied + tolerance);

            if (d_common->d_print_steps) {
               tbox::plog << d_common->d_object_name << "::supplyWork giving to " << taker << ": ";
               d_common->d_shipment->recursivePrint();
               tbox::plog << d_common->d_object_name << "::supplyWork keeping: ";
               d_common->d_local_load->recursivePrint();
               if (d_common->d_shipment->getSumLoad() > est_work_supplied + tolerance ||
                   d_common->d_shipment->getSumLoad() < est_work_supplied - tolerance) {
                  tbox::plog << "  shipment missed range: target shipment "
                             << d_common->d_shipment->getSumLoad() << " / " << est_work_supplied
                             << " [" << est_work_supplied - tolerance << ','
                             << est_work_supplied + tolerance << "] by "
                             << d_common->d_shipment->getSumLoad() - est_work_supplied
                             << " units.  kept: "
                             << d_common->d_local_load->getSumLoad()
                             << std::endl;
               }
            }

            d_process_may_supply[0] =
               d_common->d_local_load->getSumLoad() - d_obligation >
               d_common->d_pparams->getLoadComparisonTol();
         }
      }

      d_group_may_supply = estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol();

      d_work -= est_work_supplied;

   } // d_group_may_supply

   if (d_common->d_print_steps && d_common->d_print_child_steps) {
      tbox::plog << d_common->d_object_name << "::supplyWork generation "
                 << d_gen_num << " [" << d_begin << ',' << d_end << ')'
                 << " supplied estimated " << est_work_supplied << " to " << taker
                 << std::endl;
   }
   return est_work_supplied;
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::sendShipment(int taker)
{
   d_common->t_send_shipment->start();

   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::sendMyShipment: sending to " << taker << ' ';
      d_common->d_shipment->recursivePrint(tbox::plog, "", 0);
      tbox::plog << " leaving d_local_load with ";
      d_common->d_local_load->recursivePrint(tbox::plog, "", 0);
      tbox::plog << std::endl;
   }
   tbox::MessageStream msg;
   msg << *d_common->d_shipment;
   d_common->d_comm_peer[0].setPeerRank(taker);
   d_common->d_comm_peer[0].setMPITag(CascadePartitionerTree_TAG_LoadTransfer0,
      CascadePartitionerTree_TAG_LoadTransfer1);
   d_common->d_comm_peer[0].beginSend(static_cast<const char *>(msg.getBufferStart()),
      msg.getCurrentSize(), true);
   d_common->d_shipment->clear();

   d_common->t_send_shipment->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::receiveAndUnpackSuppliedLoad()
{
   d_common->t_receive_and_unpack_supplied_load->start();

   while (d_common->d_comm_stage.numberOfCompletedMembers() > 0 ||
          d_common->d_comm_stage.advanceAny()) {
      tbox::AsyncCommStage::Member* completed = d_common->d_comm_stage.popCompletionQueue();
      tbox::AsyncCommPeer<char>* comm_peer = static_cast<tbox::AsyncCommPeer<char> *>(completed);
      tbox::MessageStream recv_msg(comm_peer->getRecvSize(),
                                   tbox::MessageStream::Read,
                                   comm_peer->getRecvData(),
                                   true);
      recv_msg >> *d_common->d_shipment;
      d_common->d_local_load->insertAll(*d_common->d_shipment);
      if (d_common->d_print_steps) {
         tbox::plog << d_common->d_object_name << "::receiveAndUnpackSuppliedLoad: received ";
         d_common->d_shipment->recursivePrint(tbox::plog, "", 0);
         tbox::plog << " from process " << comm_peer->getPeerRank()
                    << " and updated d_local_load to ";
         d_common->d_local_load->recursivePrint(tbox::plog, "", 0);
         tbox::plog << std::endl;
      }
      d_common->d_shipment->clear();
   }

   d_common->t_receive_and_unpack_supplied_load->stop();
}

/*
 *************************************************************************
 * Should only be called for leaves.  For bigger groups, communication
 * is required and is determined as part of work exchange.
 *************************************************************************
 */
void
CascadePartitionerTree::recomputeLeafData()
{
   TBOX_ASSERT(this == d_leaf);   // Should only be called for leaves.
   d_work = d_common->d_local_load->getSumLoad();
   d_group_may_supply = d_process_may_supply[0] =
         estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol();
}

/*
 *************************************************************************
 * Reset the obligation of this group and its descendents based on the
 * given avg_load.
 *
 * This method is needed after a group balances its children, and the
 * children are stuck with their load, which may be higher than the
 * global average.  Because we are stuck with this work, we try to
 * distribute it evenly rather than have one process absorb all the
 * extra work.
 *************************************************************************
 */
void
CascadePartitionerTree::resetObligation(double avg_load)
{
   d_obligation = avg_load * static_cast<double>(size());
   d_group_may_supply = estimatedSurplus() > d_common->d_pparams->getLoadComparisonTol();

   if (d_children[0]) {
      d_children[0]->resetObligation(avg_load);
      d_children[1]->resetObligation(avg_load);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
double
CascadePartitionerTree::computeConnectorUpdateInterval() const
{
   const double fanout_size = d_common->d_global_work_avg >
      d_common->d_pparams->getLoadComparisonTol() ?
      d_common->d_local_work_max / d_common->d_global_work_avg : 1.0;

   /*
    * Set default number_of_updates to 1, because if fanout_size == 1.0,
    * then the computed number_of_updates will be 0, which would cause
    * a divide by 0.
    */
   int number_of_updates = 1;
   if (fanout_size > 1.0 && d_common->d_max_spread_procs > 1) {
      number_of_updates =
         static_cast<int>(ceil(log(fanout_size) / log(static_cast<double>(d_common->d_max_spread_procs))));
   }
   TBOX_ASSERT(number_of_updates > 0);

   const int tree_depth = CascadePartitioner::lgInt(d_common->d_mpi.getSize());
   const double update_interval = static_cast<double>(tree_depth) / number_of_updates;
   if (d_common->d_print_steps) {
      tbox::plog << d_common->d_object_name << "::computeConnectorUpdateInterval"
                 << "  max_spread_procs=" << d_common->d_max_spread_procs
                 << "  fanout_size=" << fanout_size
                 << "  number_of_updates=" << number_of_updates
                 << "  update_interval=" << update_interval
                 << std::endl;
   }
   return update_interval;
}

/*
 *************************************************************************
 *************************************************************************
 */
void
CascadePartitionerTree::printClassData(std::ostream& co, const std::string& border) const
{
   const std::string indent(border + std::string(d_gen_num, ' ') + std::string(d_gen_num, ' '));
   const int cycle_num = CascadePartitioner::lgInt(d_common->d_mpi.getSize()) - d_gen_num;
   co << indent << "gen_num=" << d_gen_num << "  cycle=" << cycle_num
      << "  [" << d_begin << ',' << d_end << ")  group_size=" << d_end - d_begin
      << "  local leaf=" << (this == d_leaf) << "  this=" << this
      << "  near=" << d_near << "  far=" << d_far
      << '\n' << indent
      << "contact=" << d_contact[0] << ',' << d_contact[1]
      << "  work=" << d_work << '/' << d_obligation << " (" << size() * d_common->d_global_work_avg
      << ")  estimated surplus=" << estimatedSurplus()
      << " (" << d_work - (size() * d_common->d_global_work_avg)
      << ")  local_load=" << d_common->d_local_load->getSumLoad()
      << '\n' << indent
      << "group_may_supply=" << d_group_may_supply
      << "  process_may_supply=" << d_process_may_supply[0] << ',' << d_process_may_supply[1]
      << '\n';
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

#endif
