/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_BoxTransitSet_C
#define included_mesh_BoxTransitSet_C

#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int BoxTransitSet::BoxTransitSet_EDGETAG0;
const int BoxTransitSet::BoxTransitSet_EDGETAG1;
const int BoxTransitSet::BoxTransitSet_FIRSTDATALEN;

const std::string BoxTransitSet::s_default_timer_prefix("mesh::BoxTransitSet");
std::map<std::string, BoxTransitSet::TimerStruct> BoxTransitSet::s_static_timers;

tbox::StartupShutdownManager::Handler
BoxTransitSet::s_initialize_finalize_handler(
   BoxTransitSet::initializeCallback,
   0,
   0,
   BoxTransitSet::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */
BoxTransitSet::BoxTransitSet(
   const PartitioningParams& pparams):
   TransitLoad(),
   d_set(),
   d_sumload(0.0),
   d_sumsize(0.0),
   d_pparams(&pparams),
   d_box_breaker(pparams),
   d_print_steps(false),
   d_print_pop_steps(false),
   d_print_swap_steps(false),
   d_print_break_steps(false),
   d_print_edge_steps(false),
   d_object_timers(0)
{
   getFromInput();
   setTimerPrefix(s_default_timer_prefix);
   d_box_breaker.setPrintBreakSteps(d_print_break_steps);
}

/*
 *************************************************************************
 *************************************************************************
 */
BoxTransitSet::BoxTransitSet(
   const BoxTransitSet& other,
   bool copy_load):
   TransitLoad(other),
   d_set(),
   d_sumload(0.0),
   d_sumsize(0.0),
   d_pparams(other.d_pparams),
   d_box_breaker(other.d_box_breaker),
   d_print_steps(other.d_print_steps),
   d_print_pop_steps(other.d_print_pop_steps),
   d_print_swap_steps(other.d_print_swap_steps),
   d_print_break_steps(other.d_print_break_steps),
   d_print_edge_steps(other.d_print_edge_steps),
   d_object_timers(other.d_object_timers)
{
   if (copy_load) {
      d_set = other.d_set;
      d_sumload = other.d_sumload;
      d_sumsize = other.d_sumsize;
   }
   d_box_breaker.setPrintBreakSteps(d_print_break_steps);
}

/*
 *************************************************************************
 * Initialize sets to a new (empty) container but retains current
 * supplemental data such as control and diagnostic parameters.
 *************************************************************************
 */
void BoxTransitSet::initialize()
{
   d_set.clear();
   d_sumload = 0.0;
   d_sumsize = 0.0;
}

/*
 *************************************************************************
 * Allocate a new object exactly like this, but empty.
 *************************************************************************
 */
BoxTransitSet *BoxTransitSet::clone() const
{
   BoxTransitSet* new_object = new BoxTransitSet(*this, false);
   return new_object;
}

/*
 *************************************************************************
 *************************************************************************
 */
void BoxTransitSet::insertAll(const hier::BoxContainer& other)
{
   size_t old_size = d_set.size();
   for (hier::BoxContainer::const_iterator bi = other.begin(); bi != other.end(); ++bi) {
      BoxInTransit new_box(*bi);
      d_set.insert(new_box);
      d_sumload += new_box.getLoad();
      d_sumsize += new_box.getSize();
   }
   if (d_set.size() != old_size + other.size()) {
      TBOX_ERROR("BoxTransitSet's insertAll currently can't weed out duplicates.");
   }
}

void BoxTransitSet::insertAllWithArtificialMinimum(
   const hier::BoxContainer& other,
   double minimum_load)
{
   size_t old_size = d_set.size();
   for (hier::BoxContainer::const_iterator bi = other.begin(); bi != other.end(); ++bi) {
      BoxInTransit new_box(*bi);
      new_box.setLoad(
         tbox::MathUtilities<double>::Max(new_box.getLoad(), minimum_load));
      d_set.insert(new_box);
      d_sumload += new_box.getLoad();
      d_sumsize += new_box.getSize();
   }
   if (d_set.size() != old_size + other.size()) {
      TBOX_ERROR("BoxTransitSet's insertAllWithArtificialMinimum currently can't weed out duplicates.");
   }
}


/*
 *************************************************************************
 *************************************************************************
 */
void BoxTransitSet::insertAll(TransitLoad& other_transit_load)
{
   const BoxTransitSet& other = recastTransitLoad(other_transit_load);
   size_t old_size = d_set.size();
   d_set.insert(other.d_set.begin(), other.d_set.end());
   d_sumload += other.d_sumload;
   d_sumsize += other.d_sumsize;
   if (d_set.size() != old_size + other.size()) {
      TBOX_ERROR("BoxTransitSet's insertAll currently can't weed out duplicates.");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void BoxTransitSet::insertAllWithExistingLoads(
   const hier::BoxContainer& other)
{
   std::set<BoxInTransit, BoxInTransitMoreLoad> tmp_set;

   for (iterator si = begin(); si != end(); ++si) {
      hier::BoxContainer::const_iterator itr = other.find(si->getBox());
      if (itr != other.end()) {
         BoxInTransit new_box(*itr);
         new_box.setLoad(si->getLoad());
         new_box.setCornerWeights(si->getCornerWeights());
         tmp_set.insert(new_box);
      } else {
         TBOX_ERROR("BoxTransitSet::insertAllWithExistingLoads requires that the BoxContainer input contains the same boxes as this BoxTransitSet");
      }
   }
   d_set.swap(tmp_set);

}

/*
 *************************************************************************
 *************************************************************************
 */
void BoxTransitSet::setWorkload(
   const hier::PatchLevel& patch_level,
   const int work_data_id)
{
   /*
    * Set the workload for all the BoxInTransit members of d_set based
    * on the data represented by the work_data_id.  Since we cannot
    * change the members of a set in place, we construct a temporary set
    * and then swap.
    */
   std::set<BoxInTransit, BoxInTransitMoreLoad> tmp_set;
   LoadType sumload = 0.0;
   for (iterator si = begin(); si != end(); ++si) {
      const hier::BoxId& box_id = si->getBox().getBoxId();
      const std::shared_ptr<hier::Patch>& patch =
         patch_level.getPatch(box_id);
      BoxInTransit new_transit_box(*si);
      std::vector<double> corner_weights;
      new_transit_box.setLoad(
         BalanceUtilities::computeNonUniformWorkloadOnCorners(corner_weights,
            patch,
            work_data_id,
            patch->getBox()));
      new_transit_box.setCornerWeights(corner_weights);
      sumload += new_transit_box.getLoad();
      tmp_set.insert(new_transit_box);
   }
   d_set.swap(tmp_set);
   d_sumload = sumload;
}

/*
 *************************************************************************
 *************************************************************************
 */
size_t BoxTransitSet::getNumberOfItems() const
{
   return size();
}

/*
 *************************************************************************
 *************************************************************************
 */
size_t BoxTransitSet::getNumberOfOriginatingProcesses() const
{
   std::set<int> originating_procs;
   for (const_iterator si = begin(); si != end(); ++si) {
      originating_procs.insert(si->getOrigBox().getOwnerRank());
   }
   return originating_procs.size();
}

/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and put edges in balanced<==>unbalanced Connector).
 *
 * We can generate balanced--->unbalanced edges for all boxes because
 * we have their origin info.  If the box originated locally, we can
 * generate the unbalanced--->balanced edge for them as well.
 * However, we can't generate these edges for boxes originating
 * remotely.  They are generated in
 * constructSemilocalUnbalancedToBalanced, which uses communication.
 */
void
BoxTransitSet::assignToLocalAndPopulateMaps(
   hier::BoxLevel& balanced_box_level,
   hier::MappingConnector& balanced_to_unbalanced,
   hier::MappingConnector& unbalanced_to_balanced,
   double flexible_load_tol,
   const tbox::SAMRAI_MPI& alt_mpi)
{
   NULL_USE(flexible_load_tol);

   d_object_timers->t_assign_to_local_process_and_populate_maps->start();

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::assignToLocalAndPopulateMaps: entered." << std::endl;
   }

   assignToLocal(balanced_box_level, unbalanced_to_balanced.getBase(), flexible_load_tol);
   populateMaps(balanced_to_unbalanced, unbalanced_to_balanced, alt_mpi);

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::assignToLocalAndPopulateMaps: exiting." << std::endl;
   }

   d_object_timers->t_assign_to_local_process_and_populate_maps->stop();
}

/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level).
 */
void
BoxTransitSet::assignToLocal(
   hier::BoxLevel& balanced_box_level,
   const hier::BoxLevel& unbalanced_box_level,
   double flexible_load_tol,
   const tbox::SAMRAI_MPI& alt_mpi)
{
   NULL_USE(flexible_load_tol);
   NULL_USE(alt_mpi);
   /*
    * Reassign contents to local process, assigning IDs that don't
    * conflict with current Boxes.
    */
   hier::SequentialLocalIdGenerator id_gen(
      unbalanced_box_level.getLastLocalId());
   reassignOwnership(id_gen, balanced_box_level.getMPI().getRank());

   putInBoxLevel(balanced_box_level);
}

/*
 *************************************************************************
 * We can generate balanced--->unbalanced edges for all boxes because
 * we have their origin info.  If the box originated locally, we can
 * generate the unbalanced--->balanced edge for them as well.
 * However, we can't generate these edges for boxes originating
 * remotely.  They are generated in
 * constructSemilocalUnbalancedToBalanced, which uses communication.
 */
void
BoxTransitSet::populateMaps(
   hier::MappingConnector& balanced_to_unbalanced,
   hier::MappingConnector& unbalanced_to_balanced,
   const tbox::SAMRAI_MPI& alt_mpi) const
{
   d_object_timers->t_populate_maps->start();

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::populateMaps: entered." << std::endl;
   }

   generateLocalBasedMapEdges(unbalanced_to_balanced, balanced_to_unbalanced);

   constructSemilocalUnbalancedToBalanced(
      unbalanced_to_balanced,
      alt_mpi.getCommunicator() == MPI_COMM_NULL ?
      unbalanced_to_balanced.getBase().getMPI() : alt_mpi);

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::populateMaps: exiting." << std::endl;
   }

   d_object_timers->t_populate_maps->stop();
}

/*
 *************************************************************************
 * Communicate semilocal relationships in unbalanced--->balanced
 * Connectors.  These relationships must be represented by this
 * object.  Semilocal means the local process owns either d_box or
 * getOrigBox() (not both!) of each item in this BoxTransitSet.  The
 * owner of the other doesn't have this data, so this method does the
 * necessary P2P communication to set up the transpose edges.
 *
 * Each process already knows the data in its BoxTransitSet,
 * obviously.  The idea is to acquire relevant data from other
 * processes.
 *************************************************************************
 */
void
BoxTransitSet::constructSemilocalUnbalancedToBalanced(
   hier::MappingConnector& unbalanced_to_balanced,
   const tbox::SAMRAI_MPI& mpi) const
{
   d_object_timers->t_construct_semilocal->start();

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::constructSemilocalUnbalancedToBalanced: entered."
                 << std::endl;
   }

   const hier::BoxLevel& unbalanced_box_level = unbalanced_to_balanced.getBase();
   const hier::BoxLevel& balanced_box_level = unbalanced_to_balanced.getHead();

   size_t num_cells_imported = 0;

   // Pack the imported boxes into buffers by their original owners.
   d_object_timers->t_pack_edge->start();
   std::map<int, std::shared_ptr<tbox::MessageStream> > outgoing_messages;
   for (const_iterator bi = begin(); bi != end(); ++bi) {
      const BoxInTransit& bit = *bi;
      TBOX_ASSERT(bit.getBox().getOwnerRank() == mpi.getRank());
      if (bit.getOrigBox().getOwnerRank() == mpi.getRank()) {
         // Not imported.
         continue;
      }
      num_cells_imported += bit.getBox().size();
      std::shared_ptr<tbox::MessageStream>& mstream =
         outgoing_messages[bit.getOrigBox().getOwnerRank()];
      if (!mstream) {
         mstream.reset(new tbox::MessageStream);
      }
      bit.putToMessageStream(*mstream);
   }
   d_object_timers->t_pack_edge->stop();

   /*
    * Send outgoing_messages.  Optimization for mitigating contention:
    * Start by sending to the first recipient with a rank higher than
    * the local rank.
    */

   std::map<int, std::shared_ptr<tbox::MessageStream> >::iterator recip_itr =
      outgoing_messages.upper_bound(mpi.getRank());
   if (recip_itr == outgoing_messages.end()) {
      recip_itr = outgoing_messages.begin();
   }

   size_t outgoing_messages_size = static_cast<int>(outgoing_messages.size());
   std::vector<tbox::SAMRAI_MPI::Request>
   send_requests(outgoing_messages_size, MPI_REQUEST_NULL);

   d_object_timers->t_construct_semilocal_send_edges->start();
   for (size_t send_number = 0; send_number < outgoing_messages_size; ++send_number) {

      int recipient = recip_itr->first;
      tbox::MessageStream& mstream = *recip_itr->second;

      if (d_print_edge_steps) {
         tbox::plog << "Accounting for cells on proc " << recipient << std::endl;
      }

      mpi.Isend(
         (void *)(mstream.getBufferStart()),
         static_cast<int>(mstream.getCurrentSize()),
         MPI_CHAR,
         recipient,
         BoxTransitSet_EDGETAG0,
         &send_requests[send_number]);

      ++recip_itr;
      if (recip_itr == outgoing_messages.end()) {
         recip_itr = outgoing_messages.begin();
      }

   }
   d_object_timers->t_construct_semilocal_send_edges->stop();

   TBOX_ASSERT(unbalanced_box_level.getLocalNumberOfCells() + num_cells_imported
      >= balanced_box_level.getLocalNumberOfCells());

   size_t num_unaccounted_cells =
      unbalanced_box_level.getLocalNumberOfCells()
      - balanced_box_level.getLocalNumberOfCells()
      + num_cells_imported;

   if (d_print_edge_steps) {
      tbox::plog << num_unaccounted_cells << " unaccounted cells." << std::endl;
   }

   /*
    * Receive info about exported cells from processes that now own
    * those cells.  Receive until all cells are accounted for.
    * This gives us all missing semilocal unbalanced--->balanced.
    */

   std::vector<char> incoming_message;
   BoxInTransit balanced_box_in_transit(unbalanced_box_level.getDim());

   while (num_unaccounted_cells > 0) {

      d_object_timers->t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Status status;
      mpi.Probe(MPI_ANY_SOURCE, BoxTransitSet_EDGETAG0, &status);

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count(&status, MPI_CHAR, &count);
      incoming_message.resize(count, '\0');

      mpi.Recv(
         static_cast<void *>(&incoming_message[0]),
         count,
         MPI_CHAR,
         source,
         BoxTransitSet_EDGETAG0,
         &status);
      d_object_timers->t_construct_semilocal_comm_wait->stop();

      tbox::MessageStream msg(incoming_message.size(),
                              tbox::MessageStream::Read,
                              static_cast<void *>(&incoming_message[0]),
                              false);
      const size_t old_count = num_unaccounted_cells;
      d_object_timers->t_unpack_edge->start();
      while (!msg.endOfData()) {
         balanced_box_in_transit.getFromMessageStream(msg);
         TBOX_ASSERT(balanced_box_in_transit.getBox().size() <= num_unaccounted_cells);
         unbalanced_to_balanced.insertLocalNeighbor(
            balanced_box_in_transit.getBox(),
            balanced_box_in_transit.getOrigBox().getBoxId());
         TBOX_ASSERT(num_unaccounted_cells >= balanced_box_in_transit.getBox().size());
         num_unaccounted_cells -= balanced_box_in_transit.getBox().size();
      }
      d_object_timers->t_unpack_edge->stop();

      if (d_print_edge_steps) {
         tbox::plog << "Process " << source << " accounted for "
                    << (old_count - num_unaccounted_cells) << " cells, leaving "
                    << num_unaccounted_cells << " unaccounted." << std::endl;
      }

      incoming_message.clear();
   }
   TBOX_ASSERT(num_unaccounted_cells == 0);

   // Wait for the sends to complete before clearing outgoing_messages.
   if (send_requests.size() > 0) {
      std::vector<tbox::SAMRAI_MPI::Status> status(send_requests.size());
      d_object_timers->t_construct_semilocal_comm_wait->start();
      tbox::SAMRAI_MPI::Waitall(
         static_cast<int>(send_requests.size()),
         &send_requests[0],
         &status[0]);
      d_object_timers->t_construct_semilocal_comm_wait->stop();
      outgoing_messages.clear();
   }

   if (d_print_steps || d_print_edge_steps) {
      tbox::plog << "BoxTransitSet::constructSemilocalUnbalancedToBalanced: exiting."
                 << std::endl;
   }

   d_object_timers->t_construct_semilocal->stop();
}

/*
 *************************************************************************
 * Reassign the boxes to the new owner.  Any box that isn't already
 * owned by the new owner or doesn't have a valid LocalId, is given one
 * by the SequentialLocalIdGenerator.
 *************************************************************************
 */
void
BoxTransitSet::reassignOwnership(
   hier::SequentialLocalIdGenerator& id_gen,
   int new_owner_rank)
{
   std::set<BoxInTransit, BoxInTransitMoreLoad> tmp_set;

   for (const_iterator bi = begin(); bi != end(); ++bi) {
      if (bi->getOwnerRank() != new_owner_rank ||
          !bi->getLocalId().isValid()) {
         BoxInTransit reassigned_box(
            *bi, bi->getBox(), new_owner_rank, id_gen.nextValue());
         tmp_set.insert(tmp_set.end(), reassigned_box);
      } else {
         tmp_set.insert(tmp_set.end(), *bi);
      }
   }
   d_set.swap(tmp_set);

}

/*
 *************************************************************************
 * Put all local d_box into a BoxLevel.
 * Each d_box must have a valid BoxId.
 *************************************************************************
 */
void
BoxTransitSet::putInBoxLevel(
   hier::BoxLevel& box_level) const
{
   for (iterator ni = begin(); ni != end(); ++ni) {
      TBOX_ASSERT(ni->getBox().getBoxId().isValid());
      if (ni->getBox().getOwnerRank() == box_level.getMPI().getRank()) {
         box_level.addBox(ni->getBox());
      }
   }
}

/*
 *************************************************************************
 * Generate all d_box<==>getOrigBox() mapping edges, except for those
 * that cannot be set up without communication.  These semilocal edges
 * have either a remote d_box or a remote getOrigBox().
 *
 * Each d_box must have a valid BoxId.
 *************************************************************************
 */
void
BoxTransitSet::generateLocalBasedMapEdges(
   hier::MappingConnector& unbalanced_to_balanced,
   hier::MappingConnector& balanced_to_unbalanced) const
{

   tbox::SAMRAI_MPI mpi = balanced_to_unbalanced.getBase().getMPI();

   for (iterator ni = begin(); ni != end(); ++ni) {

      const BoxInTransit& added_box = *ni;

      if (!added_box.isOriginal()) {
         // ID changed means mapping needed, but store only for local boxes.

         if (added_box.getBox().getOwnerRank() == mpi.getRank()) {
            balanced_to_unbalanced.insertLocalNeighbor(
               added_box.getOrigBox(),
               added_box.getBox().getBoxId());
         }

         if (added_box.getOrigBox().getOwnerRank() == mpi.getRank()) {
            unbalanced_to_balanced.insertLocalNeighbor(
               added_box.getBox(),
               added_box.getOrigBox().getBoxId());
         }

      }

   }
}

/*
 *************************************************************************
 *
 * This method adjusts the load in this BoxTransitSet by
 * moving work between it (main_bin) and a holding_bin.  It tries to bring
 * main_bin's load to the specified ideal_load.
 *
 * The high_load and low_load define an acceptable range around the
 * ideal_load.  As soon as the main load falls in this range, no
 * further change is tried, even if it may bring the load closer to
 * the ideal.
 *
 * This method makes a best effort and returns the amount of load
 * moved.  It can move BoxInTransits between given sets and, if needed,
 * break some BoxInTransits up to move part of the work.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Return amount of load moved to main_bin from hold_bin.  Negative
 * amount means load moved from main_bin to hold_bin.
 *
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoad(
   TransitLoad& transit_load_hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   BoxTransitSet& main_bin(*this);
   BoxTransitSet& hold_bin(recastTransitLoad(transit_load_hold_bin));

   if (d_print_steps) {
      tbox::plog << "  BoxTransitSet::adjustLoad attempting to bring main load from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " or within [" << low_load << ", " << high_load << "]."
                 << std::endl;
   }
   TBOX_ASSERT(low_load <= ideal_load);
   TBOX_ASSERT(high_load >= ideal_load);

   LoadType actual_transfer = 0;

   if ((main_bin.empty() && ideal_load <= 0) ||
       (hold_bin.empty() && main_bin.getSumLoad() < ideal_load)) {
      return actual_transfer;
   }

   d_object_timers->t_adjust_load->start();

   actual_transfer = adjustLoadByPopping(
         hold_bin,
         ideal_load,
         low_load,
         high_load);

   if (d_print_steps) {
      double balance_penalty = computeBalancePenalty(
            (main_bin.getSumLoad() - ideal_load));
      tbox::plog << "  Balance penalty after adjustLoadByPopping = "
                 << balance_penalty
                 << ", needs " << (ideal_load - main_bin.getSumLoad())
                 << " more with " << main_bin.size() << " main_bin and "
                 << hold_bin.size() << " hold_bin Boxes remaining."
                 << "\n  main_bin now has " << main_bin.getSumLoad()
                 << " in " << main_bin.size() << " boxes."
                 << std::endl;
   }

   /*
    * The algorithm cycles through a do-loop.  Each time around, we
    * try to swap some BoxInTransit between main_bin and hold_bin
    * until we have main_bin's load in [low_load,high_load] or we
    * cannot improve the actual_transfer any further.  Then, we try
    * breaking up a BoxInTransit to improve the results.  If we break
    * some BoxInTransit, we generate some more swapping options that
    * were not there before, so we loop back to try swapping again.
    *
    * If a break phase does not break any Box (and does not generate
    * more swap options), the loop will stop making changes.  We break
    * the loop at that point (and whenever we get main_bin's load in
    * the correct range).  We also break out if there is no improvement,
    * which can happen when the swp and break steps undo each other's
    * work (due to round-off errors).
    *
    * TODO: This should be a while loop.  We don't need to enter it
    * if already in range.
    */
   do {

      const LoadType old_distance_to_ideal = ideal_load - main_bin.getSumLoad();

      /*
       * Try to balance load through swapping.
       */
      LoadType swap_transfer = adjustLoadBySwapping(
            hold_bin,
            ideal_load,
            low_load,
            high_load);

      actual_transfer += swap_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
               (main_bin.getSumLoad() - ideal_load));
         tbox::plog << "  Balance penalty after adjustLoadBySwapping = "
                    << balance_penalty
                    << ", needs " << (ideal_load - main_bin.getSumLoad())
                    << " more with " << main_bin.size() << " main_bin and "
                    << hold_bin.size() << " hold_bin Boxes remaining."
                    << "\n  main_bin now has " << main_bin.getSumLoad()
                    << " in " << main_bin.size() << " boxes."
                    << std::endl;
      }

      // Skip breaking if already in range.
      if (main_bin.getSumLoad() <= high_load && main_bin.getSumLoad() >= low_load) break;

      /*
       * Skip breaking if adding/subtracting the min load overshoots the range and worsens distance to range.
       */
      if (tbox::MathUtilities<double>::Abs(main_bin.getSumLoad() - 0.5 * (high_load + low_load)) <=
          0.5 * d_pparams->getMinBoxSizeProduct()) {
         break;
      }

      if (getAllowBoxBreaking()) {
         /*
          * Assuming that we did the best we could, swapping
          * some BoxInTransit without breaking any, we now break up a Box
          * in the overloaded side for partial transfer to the
          * underloaded side.
          */
         LoadType brk_transfer = adjustLoadByBreaking(
               hold_bin,
               ideal_load,
               low_load,
               high_load);
         actual_transfer += brk_transfer;

         if (d_print_steps) {
            double balance_penalty = computeBalancePenalty(
                  (main_bin.getSumLoad() - ideal_load));
            tbox::plog << "  Balance penalty after adjustLoadByBreaking = "
                       << balance_penalty
                       << ", needs " << (ideal_load - main_bin.getSumLoad())
                       << " more with " << main_bin.size() << " main_bin and "
                       << hold_bin.size() << " hold_bin Boxes remaining."
                       << "\n  main_bin now has " << main_bin.getSumLoad()
                       << " in " << main_bin.size() << " boxes."
                       << std::endl;
         }
         if (brk_transfer == 0) {
            /*
             * If no box can be broken to improve the actual_transfer,
             * there is nothing further we can do.  The swap phase, tried
             * before the break phase, also generated no transfer, so
             * there's no point trying again.  Break out now to save
             * retrying the swap phase.
             */
            if (d_print_steps) {
               tbox::plog << "  adjustLoad stopping due to unsuccessful break."
                          << std::endl;
            }
            break;
         }
      }

      LoadType improvement =
         tbox::MathUtilities<double>::Abs(old_distance_to_ideal
            - (ideal_load - main_bin.getSumLoad()));
      if (improvement < d_pparams->getLoadComparisonTol()) {
         break;
      }

      /*
       * Now that we have broken up a Box, redo this loop to
       * see if swapping can produce a better result.
       */
   } while ((main_bin.getSumLoad() >= high_load) ||
            (main_bin.getSumLoad() <= low_load));

   if (d_print_steps) {
      const LoadType point_miss = main_bin.getSumLoad() - ideal_load;
      const LoadType range_miss =
         main_bin.getSumLoad() > high_load ? main_bin.getSumLoad() - high_load :
         main_bin.getSumLoad() < low_load ? low_load - main_bin.getSumLoad() : 0;
      tbox::plog << "  adjustLoad point_miss=" << point_miss
                 << "  range_miss="
                 << (range_miss > 0 ? " " : "") // Add space if missed range
                 << (range_miss > 0.5
          * static_cast<double>(d_pparams->getMinBoxSize().getProduct()) ? " " : "")                             // Add space if missed range by a lot
                 << range_miss
                 << "  " << main_bin.getSumLoad() << '/'
                 << ideal_load << " [" << low_load << ',' << high_load << ']'
                 << std::endl;
   }

   d_object_timers->t_adjust_load->stop();

   return actual_transfer;
}

/*
 *************************************************************************
 * Attempt bring main_bin to within a specific load range by moving
 * one box to/from it from/to hold_bin.  This method is allowed to break
 * the box and move parts of it.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadByBreaking(
   BoxTransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   LoadType actual_transfer = 0;

   if (d_print_steps) {
      tbox::plog << "  Attempting to bring this bin from "
                 << getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by breaking."
                 << std::endl;
   }

   if (getSumLoad() > high_load) {
      if (d_print_steps) {
         tbox::plog << "  BoxTransitSet::adjustLoadByBreaking reversing direction."
                    << std::endl;
      }
      // The logic below does not handle bi-directional transfers, so handle it here.
      actual_transfer = -hold_bin.adjustLoadByBreaking(
            *this,
            hold_bin.getSumLoad() - (ideal_load - getSumLoad()),
            hold_bin.getSumLoad() - (high_load - getSumLoad()),
            hold_bin.getSumLoad() - (low_load - getSumLoad()));
      return actual_transfer;
   }

   BoxTransitSet& main_bin(*this);

   TBOX_ASSERT(low_load <= ideal_load);
   TBOX_ASSERT(ideal_load <= high_load);
   TBOX_ASSERT(main_bin.getSumLoad() <= high_load);

   TBOX_ASSERT(main_bin.size() + hold_bin.size() > 0);

   d_object_timers->t_shift_loads_by_breaking->start();

   const LoadType ideal_transfer = ideal_load - main_bin.getSumLoad();
   const LoadType high_transfer = high_load - main_bin.getSumLoad();
   const LoadType low_transfer = low_load - main_bin.getSumLoad();

   if (d_print_steps) {
      tbox::plog << "    adjustLoadByBreaking asked to break off "
                 << ideal_transfer << " [" << low_transfer << ','
                 << high_transfer << "] from one of " << hold_bin.size()
                 << " Boxes to add to set of " << main_bin.size()
                 << " Boxes."
                 << std::endl;
   }

   // Data for the best cutting results so far:
   hier::BoxContainer breakoff;
   hier::BoxContainer leftover;
   double breakoff_amt = 0.0;
   BoxInTransit breakbox(d_pparams->getMinBoxSize().getDim());

   int break_acceptance_flags[4] = { 0, 0, 0, 0 };
   int& found_breakage = break_acceptance_flags[2];

   /*
    * Find best box to break.  Loop in reverse because smaller boxes
    * are cheaper to analyze for bad cuts.
    */
   for (reverse_iterator si = hold_bin.rbegin(); si != hold_bin.rend(); ++si) {

      /*
       * Skip boxes smaller than ideal_transfer.  If we called
       * adjustLoadBySwapping before entering this method, there
       * should not be any such boxes.
       */
      if (si->getLoad() < ideal_transfer) {
         continue;
      }

      const BoxInTransit& candidate = *si;

      if (d_print_steps) {
         tbox::plog << "    Considering break candidate " << candidate
                    << std::endl;
      }

      hier::BoxContainer trial_breakoff;
      hier::BoxContainer trial_leftover;
      double trial_breakoff_amt;

      d_box_breaker.breakOffLoad(
         trial_breakoff,
         trial_leftover,
         trial_breakoff_amt,
         candidate.getBox(),
         candidate.getLoad(),
         candidate.getCornerWeights(),
         ideal_transfer,
         low_transfer,
         high_transfer,
         getThresholdWidth());

      if (!trial_breakoff.empty()) {

         const bool accept_break = BalanceUtilities::compareLoads(
               break_acceptance_flags, breakoff_amt, trial_breakoff_amt,
               ideal_transfer, low_transfer, high_transfer, *d_pparams);
         if (d_print_break_steps) {
            tbox::plog << "    adjustLoadByBreaking sees potential to replace "
                       << candidate << " with "
                       << trial_breakoff.size() << " breakoff Boxes and "
                       << trial_leftover.size() << " leftover Boxes."
                       << "  break amount = " << trial_breakoff_amt
                       << "\n    Break evaluation:"
                       << "  " << break_acceptance_flags[0]
                       << "  " << break_acceptance_flags[1]
                       << "  " << break_acceptance_flags[2]
                       << "  " << break_acceptance_flags[3]
                       << std::endl;
         }

         if (accept_break) {
            breakbox = candidate;
            breakoff_amt = trial_breakoff_amt;
            breakoff.swap(trial_breakoff);
            leftover.swap(trial_leftover);
            if (break_acceptance_flags[0] == 1) {
               // We are in the [low,high] range.  That is sufficient.
               break;
            }
         }

      } else {
         if (d_print_break_steps) {
            tbox::plog << "    Break step could not break " << ideal_transfer
                       << " from hold_bin box " << candidate
                       << std::endl;
         }
      }

   }

   if (found_breakage == 1) {

      int work_data_id = d_pparams->getWorkloadDataId();

      /*
       * Remove the chosen candidate.  Put its breakoff parts
       * in main_bin and its leftover parts back into hold_bin.
       */
      hold_bin.erase(breakbox);
      size_t breakoff_size = breakoff.getTotalSizeOfBoxes();
      for (hier::BoxContainer::const_iterator bi = breakoff.begin();
           bi != breakoff.end();
           ++bi) {
         /*
          * The breakoff load (breakoff_amt) is apportioned proportionally
          * according to box size to the boxes in the breakoff container.
          * No corner weight information is stored in the resulting
          * BoxInTransits.
          *
          * When uniform load balancing is being used, all the loads should
          * be equal to the box zone counts.
          */
         BoxInTransit give_box_in_transit(
            breakbox,
            *bi,
            breakbox.getOwnerRank(),
            hier::LocalId::getInvalidId());

         if (work_data_id >= 0 && d_pparams->usingVouchers()) {
            const hier::BoxId& orig_box_id =
               give_box_in_transit.getOrigBox().getBoxId();
            const std::shared_ptr<hier::Patch>& patch =
               d_pparams->getWorkloadPatchLevel().getPatch(orig_box_id);

            std::vector<double> corner_weights;
            give_box_in_transit.setLoad(
               BalanceUtilities::computeNonUniformWorkloadOnCorners(
                  corner_weights,
                  patch,
                  work_data_id,
                  give_box_in_transit.getBox()));
            give_box_in_transit.setCornerWeights(corner_weights);
         } else {
            double load_frac = static_cast<double>(bi->size()) /
                               static_cast<double>(breakoff_size);
            give_box_in_transit.setLoad(load_frac * breakoff_amt);
            give_box_in_transit.setCornerWeights(std::vector<double>(0));
         }
         main_bin.insert(give_box_in_transit);
         actual_transfer += give_box_in_transit.getLoad();
      }
      LoadType leftover_amt = breakbox.getLoad() - 
                              static_cast<LoadType>(breakoff_amt);
      size_t leftover_size = leftover.getTotalSizeOfBoxes();
      for (hier::BoxContainer::const_iterator bi = leftover.begin();
           bi != leftover.end();
           ++bi) {
         /*
          * The leftover load (origial load minus breakoff_amt) is
          * aportioned proportionally according to box size to the boxes
          * in the leftover container.  No corner weight information is
          * stored in the resulting BoxInTransits.
          */
         BoxInTransit keep_box_in_transit(
            breakbox,
            *bi,
            breakbox.getOwnerRank(),
            hier::LocalId::getInvalidId());
         if (work_data_id >= 0 && d_pparams->usingVouchers()) {
            const hier::BoxId& orig_box_id =
               keep_box_in_transit.getOrigBox().getBoxId();
            const std::shared_ptr<hier::Patch>& patch =
               d_pparams->getWorkloadPatchLevel().getPatch(orig_box_id);

            std::vector<double> corner_weights;
            keep_box_in_transit.setLoad(
               BalanceUtilities::computeNonUniformWorkloadOnCorners(
                  corner_weights,
                  patch,
                  work_data_id,
                  keep_box_in_transit.getBox()));
            keep_box_in_transit.setCornerWeights(corner_weights);
         } else {
            double load_frac = static_cast<double>(bi->size()) /
                               static_cast<double>(leftover_size);

            keep_box_in_transit.setLoad(load_frac * leftover_amt);
            keep_box_in_transit.setCornerWeights(std::vector<double>(0));
         }
         hold_bin.insert(keep_box_in_transit);
      }
   }

   d_object_timers->t_shift_loads_by_breaking->stop();
   return actual_transfer;
}

/*
 *************************************************************************
 * Attempt to adjust the load of a main_bin by swapping boxes with
 * a hold_bin.
 *
 * Transfering a BoxInTransit from one BoxTransitSet to another
 * is considered a degenerate "swap" (a BoxInTransit is
 * swapped for nothing) handled by this function.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of main_bin
 * ideal_transfer < 0 means to raise the load of hold_bin
 * The iterative do loop may overshoot the ideal_transfer
 * and may have to swap to shift some of the load
 * back.
 *
 * Return amount of load transfered.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadBySwapping(
   BoxTransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   TBOX_ASSERT(high_load >= ideal_load);
   TBOX_ASSERT(low_load <= ideal_load);

   d_object_timers->t_adjust_load_by_swapping->start();

   BoxTransitSet& main_bin(*this);

   if (d_print_steps) {
      tbox::plog << "  Attempting to bring main_bin from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by swapping."
                 << std::endl;
   }

   bool found_swap;

   LoadType actual_transfer = 0;

   do {

      /*
       * Ammount we seek to transfer from hi to lo
       * (the "ideal" for this particular iteration).
       * Unlike ideal_transfer and actual_transfer, this quantity is positive.
       */
      LoadType rem_transfer = main_bin.getSumLoad() - ideal_load;
      LoadType low_transfer = main_bin.getSumLoad() - high_load;
      LoadType high_transfer = main_bin.getSumLoad() - low_load;
      if (d_print_swap_steps) {
         tbox::plog << "    Swap progress: " << main_bin.getSumLoad()
                    << " / " << ideal_load << " remaining transfer = "
                    << rem_transfer << " [" << low_transfer << ','
                    << high_transfer << ']' << std::endl;
      }

      LoadType swap_transfer;
      found_swap = swapLoadPair(
            main_bin,
            hold_bin,
            swap_transfer,
            rem_transfer,
            low_transfer,
            high_transfer);
      swap_transfer = -swap_transfer;

      if (found_swap) {
         actual_transfer += swap_transfer;
      }

   } while (found_swap &&
            (main_bin.getSumLoad() < low_load || main_bin.getSumLoad() > high_load));

   if (d_print_swap_steps) {
      tbox::plog << "  Final balance for adjustLoadBySwapping: "
                 << main_bin.getSumLoad() << " / " << ideal_load
                 << "  Off by " << (main_bin.getSumLoad() - ideal_load)
                 << std::endl;
   }

   d_object_timers->t_adjust_load_by_swapping->stop();

   return actual_transfer;
}

/*
 *************************************************************************
 * Attempt to adjust the load of a main_bin by popping the biggest boxes
 * from a source bin and moving them to a destination bin.
 *
 * This method should give results similar to adjustLoadBySwapping,
 * but when the boxes are small in comparison to the load changed, it
 * should be faster.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of main_bin
 * ideal_transfer < 0 means to raise the load of hold_bin
 *
 * Return amount of load transfered.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadByPopping(
   BoxTransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   TBOX_ASSERT(high_load >= ideal_load);
   TBOX_ASSERT(low_load <= ideal_load);

   d_object_timers->t_adjust_load_by_popping->start();

   BoxTransitSet& main_bin(*this);

   /*
    * Logic in this method assumes positive transfer from hold_bin
    * (the source) to main_bin (the destination).  When transfering
    * the other way, switch the roles of main_bin and hold_bin.
    */
   BoxTransitSet* src = &hold_bin;
   BoxTransitSet* dst = &main_bin;
   LoadType dst_ideal_load = ideal_load;
   LoadType dst_low_load = low_load;
   LoadType dst_high_load = high_load;

   if (main_bin.getSumLoad() > ideal_load) {

      dst_ideal_load = hold_bin.getSumLoad() + (main_bin.getSumLoad() - ideal_load);
      dst_low_load = hold_bin.getSumLoad() + (main_bin.getSumLoad() - high_load);
      dst_high_load = hold_bin.getSumLoad() + (main_bin.getSumLoad() - low_load);

      src = &main_bin;
      dst = &hold_bin;

   }

   if (d_print_pop_steps) {
      tbox::plog << "  Attempting to bring main_bin from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by popping."
                 << std::endl;
   }

   LoadType actual_transfer = 0;
   int acceptance_flags[4] = { 0, 0, 0, 0 };

   size_t num_boxes_popped = 0;

   while (!src->empty()) {

      const BoxInTransit& candidate_box = *src->begin();

      bool improved = BalanceUtilities::compareLoads(
            acceptance_flags, dst->getSumLoad(),
            dst->getSumLoad() + candidate_box.getLoad(),
            dst_ideal_load, dst_low_load, dst_high_load, *d_pparams);

      if (improved) {

         if (d_print_pop_steps) {
            tbox::plog << "    adjustLoadByPopping pop #" << num_boxes_popped
                       << ", " << candidate_box;
         }

         actual_transfer += candidate_box.getLoad();
         dst->insert(candidate_box);
         src->erase(src->begin());
         ++num_boxes_popped;

         if (d_print_pop_steps) {
            tbox::plog << ", main_bin load is " << main_bin.getSumLoad() << std::endl;
         }
      }
      if ((dst->getSumLoad() >= dst_low_load && dst->getSumLoad() <= high_load) ||
          !improved) {
         /*
          * TODO: Popping a box is so inexpensive that we should
          * really continue until !improved, instead of breaking out
          * right when the dst is within the desired range.  The only
          * argument for exiting early here is to keep as much load
          * off the transit lines as possible to avoid bandwith
          * limitations.  Experience shows that bandwidth is
          * challenged by number of box sources rather than number of
          * boxes.
          */
         break;
      }
   }

   if (d_print_pop_steps) {
      tbox::plog << "  Final result in adjustLoadByPopping: "
                 << main_bin.getSumLoad() << " / " << ideal_load
                 << "  Off by " << (main_bin.getSumLoad() - ideal_load)
                 << ".  " << num_boxes_popped << " boxes popped."
                 << std::endl;
   }

   d_object_timers->t_adjust_load_by_popping->stop();

   return actual_transfer;
}

/*
 *************************************************************************
 * Find a BoxInTransit in src and a BoxInTransit in dst which when
 * swapped results in shifting close to ideal_shift from src to dst.
 * Make the swap.  Return whether a swap pair was found.
 *************************************************************************
 */
bool
BoxTransitSet::swapLoadPair(
   BoxTransitSet& src,
   BoxTransitSet& dst,
   LoadType& actual_transfer,
   LoadType ideal_transfer,
   LoadType low_transfer,
   LoadType high_transfer) const
{
   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      bool rval = swapLoadPair(
            dst,
            src,
            actual_transfer,
            -ideal_transfer,
            -high_transfer,
            -low_transfer);
      actual_transfer = -actual_transfer;
      return rval;
   }

   d_object_timers->t_find_swap_pair->start();

   if (d_print_swap_steps) {
      tbox::plog << "    swapLoadPair looking for transfer of "
                 << ideal_transfer
                 << " between " << src.size() << "-box src and "
                 << dst.size() << "-box dst." << std::endl;
      tbox::plog << "      src (" << src.size() << "):" << std::endl;
      if (src.size() < 10) {
         for (iterator si = src.begin(); si != src.end(); ++si) {
            tbox::plog << "        " << *si << std::endl;
         }
      }
      tbox::plog << "      dst (" << dst.size() << "):" << std::endl;
      if (dst.size() < 10) {
         for (iterator si = dst.begin(); si != dst.end(); ++si) {
            tbox::plog << "        " << *si << std::endl;
         }
      }
   }

   /*
    * Look for two swap options.  The "high side" option would
    * transfer at least ideal_transfer.  The "low side" option would
    * transfer up to ideal_transfer.
    *
    * Each option is defined by a box from src and a box from dst,
    * designated by the iterators src_hiside, dst_hiside, src_loside
    * and dst_loside.  src_hiside points to the box in the src for the
    * high-side transfer, and similarly for dst_hiside.  src_loside
    * points to the box in the src for the low-side transfer, and
    * similarly for dst_loside.
    *
    * Note that in the degenerate case, the dst box does not exist,
    * and the swap degenerates to moving a box from the src to the
    * dst.
    *
    * Compute the balance_penalty if high and low were swapped.  Keep
    * looking until we find the pair giving the lowest balance_penalty
    * on swapping.
    *
    * isrc and idst point to the current best pair to swap.
    *
    * src_test and dst_test are trial pairs to check to see if we can
    * improve on new_balance_penalty.
    *
    * We will look for two "best" pairs:
    *
    * TODO: This method was originally written to compute the best
    * hiside and loside options separately and compare them at the
    * end.  That separation may not be needded anymore.  It may be
    * possible to simplify this method by keeping only the best option
    * at any time.
    */

   // Initialization indicating no swap pair found yet.
   iterator src_hiside = src.end();
   iterator dst_hiside = dst.end();
   iterator src_loside = src.end();
   iterator dst_loside = dst.end();

   // A dummy BoxInTransit for set searches.
   hier::Box dummy_box(d_pparams->getMinBoxSize().getDim());
   BoxInTransit dummy_search_target(d_pparams->getMinBoxSize().getDim());

   // Difference between swap results and ideal, >= 0
   LoadType hiside_transfer = 0.0;
   LoadType loside_transfer = 0.0;

   int loside_acceptance_flags[4] = { 0, 0, 0, 0 };
   int hiside_acceptance_flags[4] = { 0, 0, 0, 0 };

   if (dst.empty()) {
      /*
       * There is no dst BoxInTransit, so the swap would
       * degnerate to moving a box from src to dst.  Find
       * the best src BoxInTransit to move.
       */
      dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
      dummy_search_target.setLoad(ideal_transfer);
      const iterator src_test = src.lower_bound(dummy_search_target);

      if (d_print_swap_steps) {
         tbox::plog << "  swapLoadPair with empty dst: ";
      }

      if (src_test != src.begin()) {
         iterator src_test1 = src_test;
         --src_test1;
         if (BalanceUtilities::compareLoads(
                hiside_acceptance_flags, hiside_transfer,
                src_test1->getLoad(), ideal_transfer,
                low_transfer, high_transfer, *d_pparams)) {
            src_hiside = src_test1;
            hiside_transfer = src_hiside->getLoad();
            if (d_print_swap_steps) {
               tbox::plog << "  hi src: " << (*src_hiside)
                          << " with transfer " << src_hiside->getLoad()
                          << ", off by " << hiside_transfer - ideal_transfer
                          << ", acceptance_flags=" << hiside_acceptance_flags[0]
                          << ',' << hiside_acceptance_flags[1]
                          << ',' << hiside_acceptance_flags[2]
                          << ',' << hiside_acceptance_flags[3];
            }
         }
      }
      if (src_test != src.end()) {
         if (BalanceUtilities::compareLoads(
                loside_acceptance_flags, loside_transfer,
                src_test->getLoad(), ideal_transfer,
                low_transfer, high_transfer, *d_pparams)) {
            src_loside = src_test;
            loside_transfer = src_loside->getLoad();
            if (d_print_swap_steps) {
               tbox::plog << "  lo src: " << (*src_loside)
                          << " with transfer " << src_loside->getLoad()
                          << ", off by " << loside_transfer - ideal_transfer
                          << ", acceptance_flags=" << loside_acceptance_flags[0]
                          << ',' << loside_acceptance_flags[1]
                          << ',' << loside_acceptance_flags[2]
                          << ',' << loside_acceptance_flags[3];
            }
         }
      }
      if (d_print_swap_steps) {
         tbox::plog << std::endl;
      }

   } else {

      /*
       * Start search through src beginning with the box whose load
       * exceeds the biggest dst box by at least ideal_transfer.
       */
      dummy_search_target = *dst.begin();
      dummy_search_target.setLoad(dummy_search_target.getLoad() + ideal_transfer);
      iterator src_beg = src.lower_bound(dummy_search_target);

      for (iterator src_test = src_beg; src_test != src.end(); ++src_test) {

         /*
          * Set dst_test pointing to where we should start looking in dst.
          * Look for a load less than the load of src_test by
          * ideal_transfer.
          */
         dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
         dummy_search_target.setLoad(tbox::MathUtilities<LoadType>::Max(
               src_test->getLoad() - ideal_transfer,
               0));
         iterator dst_test = dst.lower_bound(dummy_search_target);

         if (dst_test != dst.end()) {

            /*
             * lower_bound returned dst_test that would transfer >=
             * ideal_transfer when swapped with src_test.  Check
             * transfererence between src_test and dst_test for the
             * high-side transfer.  Also check the next smaller box in
             * dst for the low-side transfer.
             */

            BalanceUtilities::compareLoads(
               hiside_acceptance_flags, hiside_transfer,
               src_test->getLoad() - dst_test->getLoad(),
               ideal_transfer, low_transfer, high_transfer, *d_pparams);

            if (hiside_acceptance_flags[2] == 1) {
               src_hiside = src_test;
               dst_hiside = dst_test;
               hiside_transfer = src_hiside->getLoad() - dst_hiside->getLoad();
               if (d_print_swap_steps) {
                  tbox::plog << "    new hi-swap pair: " << (*src_hiside)
                             << " & " << (*dst_hiside) << " with transfer "
                             << hiside_transfer
                             << " missing by " << hiside_transfer - ideal_transfer
                             << std::endl;
               }
            }

            if (dst_test != dst.begin()) {
               --dst_test; // Now, src_test and dst_test transferer by *less* than ideal_transfer.

               BalanceUtilities::compareLoads(
                  loside_acceptance_flags, loside_transfer,
                  src_test->getLoad() - dst_test->getLoad(),
                  ideal_transfer, low_transfer, high_transfer, *d_pparams);

               if (loside_acceptance_flags[2] == 1) {
                  src_loside = src_test;
                  dst_loside = dst_test;
                  loside_transfer = src_loside->getLoad() - dst_loside->getLoad();
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap pair: " << (*src_loside)
                                << " & " << (*dst_loside) << " with transfer "
                                << loside_transfer
                                << " missing by " << loside_transfer - ideal_transfer
                                << std::endl;
                  }
               }
            }

         } else {

            /*
             * The ideal dst to swap is smaller than the smallest dst
             * box.  So the only choice is swapping src_test for nothing.
             * Chech this against the current high- and low-side choices.
             */
            if (src_test->getLoad() > ideal_transfer) {
               // Moving src_test to dst is moving too much--hiside.

               BalanceUtilities::compareLoads(
                  hiside_acceptance_flags, hiside_transfer,
                  src_test->getLoad(), ideal_transfer,
                  low_transfer, high_transfer, *d_pparams);

               if (hiside_acceptance_flags[2] == 1) {
                  src_hiside = src_test;
                  dst_hiside = dst.end();
                  hiside_transfer = src_hiside->getLoad();
                  if (d_print_swap_steps) {
                     tbox::plog << "    new hi-swap source: " << (*src_hiside)
                                << " & " << "no dst" << " with transfer "
                                << (src_hiside->getLoad())
                                << " missing by " << hiside_transfer - ideal_transfer
                                << std::endl;
                  }
               }
            } else {
               // Moving src_test to dst is moving (just right or) too little--loside.

               BalanceUtilities::compareLoads(
                  loside_acceptance_flags, loside_transfer,
                  src_test->getLoad(), ideal_transfer,
                  low_transfer, high_transfer, *d_pparams);

               if (loside_acceptance_flags[2] == 1) {
                  src_loside = src_test;
                  dst_loside = dst.end();
                  loside_transfer = src_loside->getLoad();
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap source: " << (*src_loside)
                                << " & " << "no dst" << " with transfer "
                                << (src_loside->getLoad())
                                << " missing by " << loside_transfer - ideal_transfer
                                << std::endl;
                  }
               }
               /*
                * Break out of the loop early because there is no
                * point checking smaller src boxes.
                */
               break;
            }
         }

         if ((low_transfer <= loside_transfer && loside_transfer <= high_transfer) ||
             (low_transfer <= hiside_transfer && hiside_transfer <= high_transfer)) {
            // Found a transfer satisfying the range.  Stop searching.
            break;
         }

      }

   }

   if (d_print_swap_steps) {
      double balance_penalty_current = static_cast<double>(ideal_transfer);
      double balance_penalty_loside = static_cast<double>(loside_transfer - ideal_transfer);
      double balance_penalty_hiside = static_cast<double>(hiside_transfer - ideal_transfer);
      tbox::plog.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      tbox::plog.precision(8);
      tbox::plog << "    Swap candidates give penalties (unswap,lo,hi): "
                 << balance_penalty_current << " , " << balance_penalty_loside
                 << " , " << balance_penalty_hiside << std::endl;
   }

   bool found_swap = false;
   iterator isrc = src.end();
   iterator idst = dst.end();
   actual_transfer = 0;

   if (BalanceUtilities::compareLoads(
          hiside_acceptance_flags, actual_transfer,
          hiside_transfer, ideal_transfer,
          low_transfer, high_transfer, *d_pparams)) {
      isrc = src_hiside;
      idst = dst_hiside;
      actual_transfer = hiside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking hiside." << std::endl;
      }
   }

   if (BalanceUtilities::compareLoads(
          loside_acceptance_flags, actual_transfer,
          loside_transfer, ideal_transfer,
          low_transfer, high_transfer, *d_pparams)) {
      isrc = src_loside;
      idst = dst_loside;
      actual_transfer = loside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking loside." << std::endl;
      }
   }

   if (found_swap) {

      // We can improve balance_penalty by swapping isrc with idst.
      if (d_print_swap_steps) {
         tbox::plog << "    Swapping " << actual_transfer << " units using ";
         if (isrc != src.end()) tbox::plog << *isrc;
         else tbox::plog << "X";
         tbox::plog << " <--> ";
         if (idst != dst.end()) tbox::plog << *idst;
         else tbox::plog << "X";
         tbox::plog << std::endl;
      }

      if (isrc != src.end()) {
         dst.insert(*isrc);
         src.erase(isrc);
      }
      if (idst != dst.end()) {
         src.insert(*idst);
         dst.erase(idst);
      }

   } else {
      if (d_print_swap_steps) {
         if (isrc == src.end()) {
            tbox::plog << "    Cannot find swap pair for " << ideal_transfer
                       << " units." << std::endl;
         } else {
            tbox::plog << "    Keeping original (no swap)." << std::endl;
         }
      }
   }

   d_object_timers->t_find_swap_pair->stop();
   return found_swap;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::setPrintFlags(
   bool steps, bool pop_steps, bool swap_steps, bool break_steps, bool edge_steps)
{
   d_print_steps = steps;
   d_print_pop_steps = pop_steps;
   d_print_swap_steps = swap_steps;
   d_print_break_steps = break_steps;
   d_print_edge_steps = edge_steps;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix));
   if (ti == s_static_timers.end()) {
      d_object_timers = &s_static_timers[timer_prefix];
      getAllTimers(timer_prefix, *d_object_timers);
   } else {
      d_object_timers = &(ti->second);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_adjust_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoad()");
   timers.t_adjust_load_by_popping = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoadByPopping()");
   timers.t_adjust_load_by_swapping = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoadBySwapping()");
   timers.t_shift_loads_by_breaking = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoadByBreaking()");
   timers.t_find_swap_pair = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::swapLoadPair()");

   timers.t_assign_to_local_process_and_populate_maps = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::assignToLocalAndPopulateMaps()");

   timers.t_populate_maps = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::populateMaps()");
   timers.t_construct_semilocal = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::constructSemilocalUnbalancedToBalanced()");
   timers.t_construct_semilocal_comm_wait = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::constructSemilocalUnbalancedToBalanced()_comm_wait");
   timers.t_construct_semilocal_send_edges = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::constructSemilocalUnbalancedToBalanced()_send_edges");

   timers.t_pack_edge = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::pack_edge");
   timers.t_unpack_edge = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::unpack_edge");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::putToMessageStream(tbox::MessageStream& msg) const
{
   msg << size();
   for (const_iterator ni = begin(); ni != end(); ++ni) {
      const BoxInTransit& box_in_transit = *ni;
      box_in_transit.putToMessageStream(msg);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::getFromMessageStream(tbox::MessageStream& msg)
{
   /*
    * As we pull each BoxInTransit out, give it a new id that reflects
    * its new owner.
    */
   size_t num_boxes = 0;
   msg >> num_boxes;
   BoxInTransit received_box(d_pparams->getDim());
   for (size_t i = 0; i < num_boxes; ++i) {
      received_box.getFromMessageStream(msg);
      insert(received_box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   co << border << getSumLoad() << " units in " << size() << " boxes";
   if (detail_depth > 0) {
      size_t count = 0;
      co << ":\n";
      for (BoxTransitSet::const_iterator bi = begin();
           bi != end() && count < 10; ++bi, ++count) {
         tbox::plog << border << "    " << *bi << '\n';
      }
   }
}

/*
 *************************************************************************
 * Look for an input database called "BoxTransitSet" and read
 * parameters if it exists.
 *************************************************************************
 */

void
BoxTransitSet::getFromInput()
{
   if (!tbox::InputManager::inputDatabaseExists()) return;

   std::shared_ptr<tbox::Database> input_db = tbox::InputManager::getInputDatabase();

   if (input_db->isDatabase("BoxTransitSet")) {

      std::shared_ptr<tbox::Database> my_db = input_db->getDatabase("BoxTransitSet");

      d_print_steps = my_db->getBoolWithDefault("DEV_print_steps", d_print_steps);
      d_print_break_steps =
         my_db->getBoolWithDefault("DEV_print_break_steps", d_print_break_steps);
      d_print_pop_steps =
         my_db->getBoolWithDefault("DEV_print_pop_steps", d_print_pop_steps);
      d_print_swap_steps =
         my_db->getBoolWithDefault("DEV_print_swap_steps", d_print_swap_steps);
      d_print_edge_steps =
         my_db->getBoolWithDefault("DEV_print_edge_steps", d_print_edge_steps);

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

#endif
