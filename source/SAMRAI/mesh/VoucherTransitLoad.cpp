/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TransitLoad by vouchers.
 *
 ************************************************************************/

#ifndef included_mesh_VoucherTransitLoad_C
#define included_mesh_VoucherTransitLoad_C

#include "SAMRAI/mesh/VoucherTransitLoad.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int VoucherTransitLoad::VoucherTransitLoad_DEMANDTAG;
const int VoucherTransitLoad::VoucherTransitLoad_SUPPLYTAG;

const std::string VoucherTransitLoad::s_default_timer_prefix("mesh::VoucherTransitLoad");
std::map<std::string, VoucherTransitLoad::TimerStruct> VoucherTransitLoad::s_static_timers;

tbox::StartupShutdownManager::Handler
VoucherTransitLoad::s_initialize_finalize_handler(
   VoucherTransitLoad::initializeCallback,
   0,
   0,
   VoucherTransitLoad::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *************************************************************************
 */
VoucherTransitLoad::VoucherTransitLoad(
   const PartitioningParams& pparams):
   TransitLoad(),
   d_vouchers(),
   d_sumload(0),
   d_pparams(&pparams),
   d_partition_work_supply_recursively(true),
   d_flexible_load_tol(0.0),
   d_reserve(pparams),
   d_print_steps(false),
   d_print_edge_steps(false),
   d_object_timers(0)
{
   getFromInput();
   setTimerPrefix(s_default_timer_prefix);
}

/*
 *************************************************************************
 *************************************************************************
 */
VoucherTransitLoad::VoucherTransitLoad(
   const VoucherTransitLoad& other,
   bool copy_load):
   TransitLoad(other),
   d_vouchers(),
   d_sumload(0),
   d_pparams(other.d_pparams),
   d_partition_work_supply_recursively(other.d_partition_work_supply_recursively),
   d_flexible_load_tol(other.d_flexible_load_tol),
   d_reserve(*other.d_pparams),
   d_print_steps(other.d_print_steps),
   d_print_edge_steps(other.d_print_edge_steps),
   d_object_timers(other.d_object_timers)
{
   if (copy_load) {
      d_vouchers = other.d_vouchers;
      d_sumload = other.d_sumload;
   }
}

/*
 *************************************************************************
 * Initialize sets to a new (empty) container but retains current
 * supplemental data such as control and diagnostic parameters.
 *************************************************************************
 */
void VoucherTransitLoad::initialize()
{
   d_vouchers.clear();
   d_sumload = 0.0;
}

/*
 *************************************************************************
 * Allocate a new object exactly like this, but empty.
 *************************************************************************
 */
VoucherTransitLoad *VoucherTransitLoad::clone() const
{
   VoucherTransitLoad* new_object = new VoucherTransitLoad(*this, false);
   return new_object;
}

/*
 *************************************************************************
 *************************************************************************
 */
void VoucherTransitLoad::insertAll(const hier::BoxContainer& other)
{
   for (hier::BoxContainer::const_iterator bi = other.begin(); bi != other.end(); ++bi) {
      insertCombine(Voucher(LoadType(bi->size()), bi->getOwnerRank()));
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void VoucherTransitLoad::insertAllWithArtificialMinimum(
   const hier::BoxContainer& other,
   double minimum_load)
{
   for (hier::BoxContainer::const_iterator bi = other.begin(); bi != other.end(); ++bi) {
      LoadType load = static_cast<LoadType>(bi->size());
      if (load < minimum_load) {
         load = static_cast<LoadType>(minimum_load);
      }
      insertCombine(Voucher(load, bi->getOwnerRank()));
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void VoucherTransitLoad::insertAll(TransitLoad& other_transit_load)
{
   const VoucherTransitLoad& other = recastTransitLoad(other_transit_load);
   for (const_iterator si = other.d_vouchers.begin();
        si != other.d_vouchers.end(); ++si) {
      TBOX_ASSERT(si->d_load >= d_pparams->getLoadComparisonTol());
      insertCombine(*si);
   }
}


/*
 *************************************************************************
 *************************************************************************
 */
void VoucherTransitLoad::setWorkload(
   const hier::PatchLevel& patch_level,
   const int work_data_id)
{
   d_reserve.clear();
   d_reserve.setAllowBoxBreaking(getAllowBoxBreaking());
   d_reserve.setThresholdWidth(getThresholdWidth());

   /*
    */
   clear();
   LoadType sumload = 0.0;
   for (hier::PatchLevel::iterator pi = patch_level.begin();
        pi != patch_level.end(); ++pi) {
      const std::shared_ptr<hier::Patch>& patch = *pi;
      const hier::BoxId& box_id = patch->getBox().getBoxId();
      BoxInTransit new_transit_box(patch->getBox());
      std::vector<double> corner_weights;
      double load =
         BalanceUtilities::computeNonUniformWorkloadOnCorners(corner_weights,
            patch,
            work_data_id,
            new_transit_box.getBox());
      const double artificial_minimum = d_pparams->getArtificialMinimumLoad();
      if (load < artificial_minimum) {
         load = artificial_minimum;
      }
      new_transit_box.setLoad(load);
      new_transit_box.setCornerWeights(corner_weights);
      sumload += new_transit_box.getLoad();
      d_reserve.insert(new_transit_box);
      insertCombine(Voucher(LoadType(new_transit_box.getLoad()),
                    box_id.getOwnerRank()));
   }
   d_sumload = sumload;
}


/*
 *************************************************************************
 *************************************************************************
 */
size_t VoucherTransitLoad::getNumberOfItems() const
{
   return d_vouchers.size();
}

/*
 *************************************************************************
 *************************************************************************
 */
size_t VoucherTransitLoad::getNumberOfOriginatingProcesses() const
{
   return d_vouchers.size();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::putToMessageStream(tbox::MessageStream& msg) const
{
   msg << d_vouchers.size();
   for (const_iterator ni = d_vouchers.begin(); ni != d_vouchers.end(); ++ni) {
      msg << *ni;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::getFromMessageStream(tbox::MessageStream& msg)
{
   size_t num_vouchers = 0;
   msg >> num_vouchers;
   Voucher v;
   for (size_t i = 0; i < num_vouchers; ++i) {
      msg >> v;
      insert(v);
   }
}

/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level).
 *************************************************************************
 */
void
VoucherTransitLoad::assignToLocal(
   hier::BoxLevel& balanced_box_level,
   const hier::BoxLevel& unbalanced_box_level,
   double flexible_load_tol,
   const tbox::SAMRAI_MPI& alt_mpi)
{
   d_object_timers->t_assign_to_local->start();
   // Delegate to assignToLocalAndPopulateMap but give it no map to work with.
   assignToLocalAndPopulateMaps(
      balanced_box_level,
      0, 0,
      unbalanced_box_level,
      flexible_load_tol,
      alt_mpi);
   d_object_timers->t_assign_to_local->stop();
}

/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and populate balanced<==>unbalanced).
 */
void
VoucherTransitLoad::assignToLocalAndPopulateMaps(
   hier::BoxLevel& balanced_box_level,
   hier::MappingConnector& balanced_to_unbalanced,
   hier::MappingConnector& unbalanced_to_balanced,
   double flexible_load_tol,
   const tbox::SAMRAI_MPI& alt_mpi)
{
   // Delegate to the more general assignToLocalAndPopulateMap.
   assignToLocalAndPopulateMaps(
      balanced_box_level,
      &balanced_to_unbalanced,
      &unbalanced_to_balanced,
      unbalanced_to_balanced.getBase(),
      flexible_load_tol,
      alt_mpi);
}

/*
 *************************************************************************
 * Assign boxes to local process (put them in the balanced_box_level
 * and populate balanced<==>unbalanced).
 *
 * This method does two things:
 * - Voucher redeemers request and receive work vouchers they hold.
 * - Voucher fulfillers receive and fulfill redemption requests.
 * The code is writen to let each process be both redeemers and
 * fulfillers.  Logic should drop through correctly on processes
 * that plays just one role.
 *
 * There are four major steps, organized to overlap communication.
 * 1. Request work for vouchers to be redeemed.
 * 2. Receive redemption requests.
 * 3. Fulfill redemption requests.
 * 4. Receive work for redeemed vouchers.
 *
 * To combine communicating vouchers and map edges, this method does
 * both at the same time.  If the maps are omitted (NULL), populating
 * maps is omitted.
 */
void
VoucherTransitLoad::assignToLocalAndPopulateMaps(
   hier::BoxLevel& balanced_box_level,
   hier::MappingConnector* balanced_to_unbalanced,
   hier::MappingConnector* unbalanced_to_balanced,
   const hier::BoxLevel& unbalanced_box_level,
   double flexible_load_tol,
   const tbox::SAMRAI_MPI& alt_mpi)
{
   d_object_timers->t_assign_to_local_and_populate_maps->start();

   d_flexible_load_tol = flexible_load_tol;

   if (d_print_steps) {
      tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps: entered."
                 << "\npparams: " << *d_pparams
                 << std::endl;
   }

   const tbox::SAMRAI_MPI& mpi = alt_mpi.hasNullCommunicator() ?
      unbalanced_box_level.getMPI() : alt_mpi;

   /*
    * unaccounted_work is amount of work we started with, minus what
    * we can still acount for (the part we still hold).  The rest have
    * been sent off in the form of vouchers and we don't know where
    * they ended up.
    */
   LoadType original_work;
   if (d_reserve.empty()) {
      const double artificial_minimum = d_pparams->getArtificialMinimumLoad();
      if (artificial_minimum > 1.0) {
         original_work = 0.0;
         const hier::BoxContainer& boxes = unbalanced_box_level.getBoxes();
         for (hier::BoxContainer::const_iterator bi = boxes.begin();
              bi != boxes.end(); ++bi) {
            double load = static_cast<double>(bi->size());
            if (load < artificial_minimum) {
               load = artificial_minimum;
            }
            original_work += static_cast<LoadType>(load);
         }
      } else {
         original_work = LoadType(unbalanced_box_level.getLocalNumberOfCells());
      }
   } else {
      original_work = d_reserve.getSumLoad();
   }
   LoadType unaccounted_work = original_work
      - findVoucher(mpi.getRank()).d_load;

   if (d_print_edge_steps) {
      tbox::plog << "unaccounted work is " << unaccounted_work << std::endl;
   }

   // 1. Send work demands for vouchers we want to redeem.

   std::map<int, VoucherRedemption> redemptions_to_request;
   std::map<int, VoucherRedemption> redemptions_to_fulfill;

   hier::LocalId local_id_offset = unbalanced_box_level.getLastLocalId();
   const hier::LocalId local_id_inc(static_cast<int>(d_vouchers.size()));

   for (const_iterator si = d_vouchers.begin(); si != d_vouchers.end(); ++si) {
      hier::SequentialLocalIdGenerator id_gen(++local_id_offset, local_id_inc);
      if (si->d_issuer_rank != mpi.getRank()) {
         if (d_print_edge_steps) {
            tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                       << " sending demand for voucher " << *si << '.'
                       << std::endl;
         }
         redemptions_to_request[si->d_issuer_rank].sendWorkDemand(
            si, id_gen, mpi);
      } else {
         // Locally fulfilled.  Place in redemptions_to_fulfill for processing.
         redemptions_to_fulfill[mpi.getRank()].setLocalRedemption(si, id_gen, mpi);
      }
   }

   // Set up the reserve for fulfilling incoming redemption requests.
   if (d_reserve.empty()) {
      d_reserve.clear();
      d_reserve.setAllowBoxBreaking(getAllowBoxBreaking());
      d_reserve.setThresholdWidth(getThresholdWidth());
      const double artificial_minimum = d_pparams->getArtificialMinimumLoad();
      if (artificial_minimum > 1.0) {
         d_reserve.insertAllWithArtificialMinimum(
            unbalanced_box_level.getBoxes(),
            artificial_minimum);
      } else {
         d_reserve.insertAll(unbalanced_box_level.getBoxes());
      }
   }
   if (d_print_edge_steps) {
      tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                 << " reserve before redemption steps: "
                 << d_reserve.format();
      tbox::plog << std::endl;
   }

   // 2. Receive work demands for voucher we generated but can't account for.

   /*
    * If there is unaccounted work, then some processes must have our
    * voucher.  Reveive their demands for work until we've accounted
    * for everything.  Don't supply work until all demands are
    * received, because that leads to non-deterministic results.
    */

   while (unaccounted_work > d_pparams->getLoadComparisonTol()) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe(MPI_ANY_SOURCE, VoucherTransitLoad_DEMANDTAG, &status);

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count(&status, MPI_CHAR, &count);

      VoucherRedemption& vr = redemptions_to_fulfill[source];
      vr.recvWorkDemand(source, count, mpi);
      unaccounted_work -= vr.d_voucher.d_load;

      if (d_print_edge_steps) {
         tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                    << " received demand from " << source << " for voucher "
                    << vr.d_voucher << " leaving unaccounted_work = " << unaccounted_work
                    << std::endl;
      }

      TBOX_ASSERT(unaccounted_work >= -d_pparams->getLoadComparisonTol());

   }

   // 3. Supply work according to received demands.

   if (d_partition_work_supply_recursively) {
      if (!redemptions_to_fulfill.empty()) {
         recursiveSendWorkSupply(
            redemptions_to_fulfill.begin(),
            redemptions_to_fulfill.end(),
            d_reserve);

         // Done fulfilling redemptions.  The rest stay on local process.
         for (std::map<int, VoucherRedemption>::const_iterator mi = redemptions_to_fulfill.begin();
              mi != redemptions_to_fulfill.end(); ++mi) {
            mi->second.d_box_shipment->putInBoxLevel(balanced_box_level);
            if (unbalanced_to_balanced) {
               mi->second.d_box_shipment->generateLocalBasedMapEdges(
                  *unbalanced_to_balanced,
                  *balanced_to_unbalanced);
            }
         }
      }
   } else {
      for (std::map<int, VoucherRedemption>::iterator mi = redemptions_to_fulfill.begin();
           mi != redemptions_to_fulfill.end(); ++mi) {

         VoucherRedemption& vr = mi->second;
         if (vr.d_demander_rank != mpi.getRank()) {
            vr.sendWorkSupply(d_reserve, flexible_load_tol, *d_pparams, false);
            if (d_print_edge_steps) {
               tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                          << " sent supply to " << mi->first << " for voucher "
                          << vr.d_voucher << ": "
                          << vr.d_box_shipment->format()
                          << std::endl;
            }

            if (unbalanced_to_balanced) {
               vr.d_box_shipment->generateLocalBasedMapEdges(
                  *unbalanced_to_balanced,
                  *balanced_to_unbalanced);
            }
         } else {
            vr.fulfillLocalRedemption(d_reserve, *d_pparams, false);
            vr.d_box_shipment->putInBoxLevel(balanced_box_level);
            if (unbalanced_to_balanced) {
               vr.d_box_shipment->generateLocalBasedMapEdges(
                  *unbalanced_to_balanced,
                  *balanced_to_unbalanced);
            }
         }
      }

      // Anything left in d_reserve is kept locally.
      hier::SequentialLocalIdGenerator id_gen(unbalanced_box_level.getLastLocalId(), local_id_inc);
      d_reserve.reassignOwnership(
         id_gen,
         balanced_box_level.getMPI().getRank());
   }

   if (d_print_edge_steps) {
      tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                 << " reserve after sending work supplies: "
                 << d_reserve.format()
                 << std::endl;
   }

   d_reserve.putInBoxLevel(balanced_box_level);
   if (unbalanced_to_balanced) {
      d_reserve.generateLocalBasedMapEdges(
         *unbalanced_to_balanced,
         *balanced_to_unbalanced);
   }

   // 4. Receive work according to the demands we sent.

   while (!redemptions_to_request.empty()) {

      tbox::SAMRAI_MPI::Status status;
      mpi.Probe(MPI_ANY_SOURCE, VoucherTransitLoad_SUPPLYTAG, &status);

      int source = status.MPI_SOURCE;
      int count = -1;
      tbox::SAMRAI_MPI::Get_count(&status, MPI_CHAR, &count);

      VoucherRedemption& vr = redemptions_to_request[source];
      vr.recvWorkSupply(count, *d_pparams);

      vr.d_box_shipment->putInBoxLevel(balanced_box_level);
      if (unbalanced_to_balanced) {
         vr.d_box_shipment->generateLocalBasedMapEdges(
            *unbalanced_to_balanced,
            *balanced_to_unbalanced);
      }

      if (d_print_edge_steps) {
         tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps:"
                    << " received supply from rank " << source << " for " << vr.d_voucher
                    << ": " << vr.d_box_shipment->format()
                    << std::endl;
      }

      redemptions_to_request.erase(source);
   }

   if (d_print_steps) {
      tbox::plog << "VoucherTransitLoad::assignToLocalAndPopulateMaps: exiting." << std::endl;
   }

   d_flexible_load_tol = 0.0;

   d_object_timers->t_assign_to_local_and_populate_maps->stop();
}

/*
 *************************************************************************
 * Alternative option to recursively partition work supply.
 * This version tries to avoid cutting small amounts out of big
 * boxes, which unavoidablly generates slivers.
 *
 * This method splits the reserve into a left half and a right half.
 * It gives the left side to demanders in first half of [begin,end) and
 * the right side to the second half.
 *************************************************************************
 */
void VoucherTransitLoad::recursiveSendWorkSupply(
   const std::map<int, VoucherRedemption>::iterator& begin,
   const std::map<int, VoucherRedemption>::iterator& end,
   BoxTransitSet& reserve)
{
   std::map<int, VoucherRedemption>::iterator left_begin = begin;
   std::map<int, VoucherRedemption>::iterator left_end = begin;
   ++left_end;

   std::map<int, VoucherRedemption>::iterator right_end = end;
   std::map<int, VoucherRedemption>::iterator right_begin = end;
   --right_begin;

   if (right_begin == left_begin) {
      if (begin->second.d_demander_rank != begin->second.d_mpi.getRank()) {
         begin->second.sendWorkSupply(reserve, d_flexible_load_tol, *d_pparams, true);
      } else {
         begin->second.fulfillLocalRedemption(reserve, *d_pparams, true);
      }

      if (d_print_edge_steps) {
         tbox::plog << "VoucherTransitLoad::recursiveSendWorkSupply:"
                    << " sent supply to " << begin->first << " for voucher "
                    << begin->second.d_voucher << ": "
                    << begin->second.d_box_shipment->format();
         tbox::plog << std::endl;
      }

   } else {

      double left_load = left_begin->second.d_voucher.d_load;
      double right_load = right_begin->second.d_voucher.d_load;

      // Find midpoint, where left side ends and right side begins.
      while (left_end != right_begin) {
         left_load += left_end->second.d_voucher.d_load;
         ++left_end;
         if (left_end != right_begin) {
            --right_begin;
            right_load += right_begin->second.d_voucher.d_load;
         }
      }

      BoxTransitSet left_reserve(reserve, false);
      double upper_lim = left_load * (1 + d_flexible_load_tol);
      double lower_lim = left_load * (1 - d_flexible_load_tol);
      lower_lim = tbox::MathUtilities<double>::Max(
            lower_lim, left_load - right_load * d_flexible_load_tol);
      upper_lim = tbox::MathUtilities<double>::Min(
            upper_lim, left_load + right_load * d_flexible_load_tol);
      if (d_print_edge_steps) {
         tbox::plog << "VoucherTransitLoad::recursiveSendWorkSupply: splitting reserve:"
                    << reserve.format()
                    << "to get left_load of " << left_load << " [" << lower_lim << ',' << upper_lim
                    << " and right_load of " << right_load << "\n";
      }

      left_reserve.adjustLoad(reserve, left_load, lower_lim, upper_lim);

      if (d_print_edge_steps) {
         tbox::plog << "VoucherTransitLoad::recursiveSendWorkSupply split reserve into left: "
                    << left_reserve.format() << "and right: " << reserve.format()
                    << std::endl;
      }

      recursiveSendWorkSupply(left_begin, left_end, left_reserve);
      recursiveSendWorkSupply(right_begin, right_end, reserve);

      reserve.insertAll(left_reserve);

   }
}

/*
 *************************************************************************
 * Send a demand for work to the voucher issuer, starting the voucher
 * redemption sequence of steps.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::sendWorkDemand(
   const VoucherTransitLoad::const_iterator& voucher,
   const hier::SequentialLocalIdGenerator& id_gen,
   const tbox::SAMRAI_MPI& mpi)
{
   d_voucher = *voucher;
   d_demander_rank = mpi.getRank();
   d_id_gen = id_gen;
   d_mpi = mpi;

   d_msg = std::make_shared<tbox::MessageStream>();
   (*d_msg) << d_id_gen << d_voucher;

   d_mpi.Isend(
      (void *)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_voucher.d_issuer_rank,
      VoucherTransitLoad_DEMANDTAG,
      &d_mpi_request);
}

/*
 *************************************************************************
 * Receive a demand for work from a process holding a locally issued
 * voucher.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::recvWorkDemand(
   int demander_rank,
   int message_length,
   const tbox::SAMRAI_MPI& mpi)
{
   d_mpi = mpi;
   d_demander_rank = demander_rank;

   std::vector<char> incoming_message(message_length);
   tbox::SAMRAI_MPI::Status status;
   d_mpi.Recv(
      static_cast<void *>(&incoming_message[0]),
      message_length,
      MPI_CHAR,
      d_demander_rank,
      VoucherTransitLoad_DEMANDTAG,
      &status);

   d_msg = std::make_shared<tbox::MessageStream>(
         message_length, tbox::MessageStream::Read,
         static_cast<void *>(&incoming_message[0]), false);
   (*d_msg) >> d_id_gen >> d_voucher;
   TBOX_ASSERT(d_msg->endOfData());
   TBOX_ASSERT(d_voucher.d_issuer_rank == mpi.getRank());

   d_msg.reset();
}

/*
 *************************************************************************
 * Send a supply of work to the demander to fulfill a voucher.  The
 * work is taken from a reserve.  Save mapping edges incident from
 * local boxes.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::sendWorkSupply(
   BoxTransitSet& reserve,
   double flexible_load_tol,
   const PartitioningParams& pparams,
   bool send_all)
{
   d_pparams = &pparams;
   d_box_shipment = std::make_shared<BoxTransitSet>(pparams);
   d_box_shipment->setAllowBoxBreaking(reserve.getAllowBoxBreaking());
   d_box_shipment->setThresholdWidth(reserve.getThresholdWidth());
   if (send_all) {
      d_box_shipment->swap(reserve);
   } else {
      double upper_lim = d_voucher.d_load * (1 + flexible_load_tol);
      double lower_lim = d_voucher.d_load * (1 - flexible_load_tol);
      lower_lim = tbox::MathUtilities<double>::Max(
            lower_lim,
            d_voucher.d_load - (reserve.getSumLoad() - d_voucher.d_load) * flexible_load_tol);
      d_box_shipment->adjustLoad(reserve,
         d_voucher.d_load,
         lower_lim,
         upper_lim);
   }
   d_box_shipment->reassignOwnership(d_id_gen, d_demander_rank);

   d_msg = std::make_shared<tbox::MessageStream>();
   d_box_shipment->putToMessageStream(*d_msg);

   d_mpi.Isend(
      (void *)(d_msg->getBufferStart()),
      static_cast<int>(d_msg->getCurrentSize()),
      MPI_CHAR,
      d_demander_rank,
      VoucherTransitLoad_SUPPLYTAG,
      &d_mpi_request);
}

/*
 *************************************************************************
 * Receive work supply from the issuer of the voucher.  Save mapping
 * edges incident from local boxes.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::recvWorkSupply(
   int message_length,
   const PartitioningParams& pparams)
{
   std::vector<char> incoming_message(message_length);
   tbox::SAMRAI_MPI::Status status;
   d_mpi.Recv(static_cast<void *>(&incoming_message[0]),
      message_length,
      MPI_CHAR,
      d_voucher.d_issuer_rank,
      VoucherTransitLoad_SUPPLYTAG,
      &status);

   d_msg = std::make_shared<tbox::MessageStream>(
         message_length, tbox::MessageStream::Read,
         static_cast<void *>(&incoming_message[0]), false);
   d_box_shipment = std::make_shared<BoxTransitSet>(pparams);
   d_box_shipment->getFromMessageStream(*d_msg);
   d_msg.reset();
}

/*
 *************************************************************************
 * Set a demand for a local redemption.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::setLocalRedemption(
   const VoucherTransitLoad::const_iterator& voucher,
   const hier::SequentialLocalIdGenerator& id_gen,
   const tbox::SAMRAI_MPI& mpi)
{
   d_voucher = *voucher;
   d_demander_rank = mpi.getRank();
   d_id_gen = id_gen;
   d_mpi = mpi;
}

/*
 *************************************************************************
 * Supply work to fulfill local redemption.
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::fulfillLocalRedemption(
   BoxTransitSet& reserve,
   const PartitioningParams& pparams,
   bool all)
{
   d_pparams = &pparams;
   d_box_shipment = std::shared_ptr<BoxTransitSet>(reserve.clone());
   if (all) {
      d_box_shipment->swap(reserve);
   } else {
      d_box_shipment->adjustLoad(reserve,
         d_voucher.d_load,
         d_voucher.d_load,
         d_voucher.d_load);
   }
   d_box_shipment->reassignOwnership(d_id_gen, d_demander_rank);
}

/*
 *************************************************************************
 *************************************************************************
 */
void VoucherTransitLoad::VoucherRedemption::finishSendRequest()
{
   if (d_mpi_request != MPI_REQUEST_NULL) {
      tbox::SAMRAI_MPI::Status status;
      tbox::SAMRAI_MPI::Wait(&d_mpi_request, &status);
   }
   TBOX_ASSERT(d_mpi_request == MPI_REQUEST_NULL);
}

/*
 *************************************************************************
 * Adjust the VoucherTransitLoad by moving work between it (main_bin)
 * and another (hold_bin).  Try to bring the load to the specified
 * ideal.
 *
 * Move Vouchers between given bins and, if needed, break some Vouchers
 * up to move part of the work.
 *
 * The high_load and low_load define an acceptable range around the
 * ideal_load.  Because vouchers can be cut to any arbitrary amount
 * (unlike boxes), it only uses the ideal_load.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Return amount of load moved to this object from hold_bin.  Negative
 * amount means load moved from this object to hold_bin.
 *************************************************************************
 */
VoucherTransitLoad::LoadType
VoucherTransitLoad::adjustLoad(
   TransitLoad& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   VoucherTransitLoad& main_bin(*this);

   if (d_print_steps) {
      tbox::plog << "  adjustLoad attempting to bring main load from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " or within [" << low_load << ", " << high_load << "]."
                 << std::endl;
   }
   TBOX_ASSERT(low_load <= ideal_load);
   TBOX_ASSERT(high_load >= ideal_load);

   LoadType actual_transfer = 0;

   if ((main_bin.empty() && ideal_load <= d_pparams->getLoadComparisonTol()) ||
       (hold_bin.empty() && main_bin.getSumLoad() < ideal_load
        + d_pparams->getLoadComparisonTol())) {
      return actual_transfer;
   }

   d_object_timers->t_adjust_load->start();

   const LoadType change = ideal_load - main_bin.getSumLoad();

   if (change > 0) {
      // Move load to main_bin.
      actual_transfer =
         raiseDstLoad(recastTransitLoad(hold_bin),
            main_bin,
            ideal_load);
   } else if (change < 0) {
      // Move load to hold_bin.
      actual_transfer =
         -raiseDstLoad(main_bin,
            recastTransitLoad(hold_bin),
            hold_bin.getSumLoad() - change);
   }

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
 * Raise load of dst container by shifing load from src.
 *************************************************************************
 */
VoucherTransitLoad::LoadType
VoucherTransitLoad::raiseDstLoad(
   VoucherTransitLoad& src,
   VoucherTransitLoad& dst,
   LoadType ideal_dst_load)
{
   TBOX_ASSERT(ideal_dst_load >= dst.getSumLoad());

   if (src.empty()) {
      return 0;
      // No-op empty-container cases is not handled below.
   }

   d_object_timers->t_raise_dst_load->start();

   /*
    * Decide whether to take work from the beginning or the end of
    * src, whichever is closer to the dst.  This is a minor
    * optimization to better preserve locality when the issuer rank
    * ranges of src and dst don't overlap much.
    */
   bool take_from_src_end = true;
   if (!dst.empty()) {

      const LoadType gap_at_src_end =
         tbox::MathUtilities<double>::Abs(
            src.rbegin()->d_issuer_rank - dst.begin()->d_issuer_rank);

      const LoadType gap_at_src_begin =
         tbox::MathUtilities<double>::Abs(
            dst.rbegin()->d_issuer_rank - src.begin()->d_issuer_rank);

      take_from_src_end = gap_at_src_end < gap_at_src_begin;
   }

   LoadType old_dst_load = dst.getSumLoad();

   int num_src_vouchers = src.getNumberOfItems();
   int num_iterations = -1;
   while (dst.getSumLoad() < ideal_dst_load - dst.d_pparams->getLoadComparisonTol() &&
          !src.empty()) {
      ++num_iterations;
      LoadType give_load = ideal_dst_load - dst.getSumLoad();

      iterator src_itr;
      if (take_from_src_end) {
         if (num_iterations < num_src_vouchers) {
            src_itr = src.end();
            do {
               --src_itr;
            } while (((give_load < src_itr->d_load && 2.0*give_load > src_itr->d_load) || (src_itr->d_load < give_load && (give_load - src_itr->d_load) < 0.5*give_load)) && src_itr != src.begin());
         } else {
            src_itr = src.end();
            --src_itr;
         }
      } else {
         src_itr = src.begin();
      }

      Voucher free_voucher = *src_itr;
      src.erase(src_itr);

      if (free_voucher.d_load <
          (ideal_dst_load - dst.getSumLoad() + dst.d_pparams->getLoadComparisonTol())) {
         dst.insert(free_voucher);
      } else {
         Voucher partial_voucher((ideal_dst_load - dst.getSumLoad()), free_voucher);
         dst.insert(partial_voucher);
         src.insert(free_voucher);
         /*
          * Breaking not strictly needed, except rounding error may
          * cause adding another partial voucher from same issuer.
          */
         break;
      }
   }

   d_object_timers->t_raise_dst_load->stop();
   return dst.getSumLoad() - old_dst_load;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
VoucherTransitLoad::eraseIssuer(int issuer_rank)
{
   Voucher tmp_voucher(0.0, issuer_rank);
   const_iterator vi = d_vouchers.lower_bound(tmp_voucher);
   if (vi != d_vouchers.end() && vi->d_issuer_rank == issuer_rank) {
#ifdef DEBUG_CHECK_ASSERTIONS
      const_iterator vi1 = vi;
      ++vi1;
      TBOX_ASSERT(vi1 == d_vouchers.end() || vi1->d_issuer_rank != issuer_rank);
#endif
      d_sumload -= vi->d_load;
      d_vouchers.erase(vi);
      return true;
   }
   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
VoucherTransitLoad::Voucher
VoucherTransitLoad::findVoucher(int issuer_rank) const
{
   Voucher tmp_voucher(0.0, issuer_rank);
   const_iterator vi = d_vouchers.lower_bound(tmp_voucher);
   if (vi != d_vouchers.end() && vi->d_issuer_rank == issuer_rank) {
      tmp_voucher.d_load = vi->d_load;
   }
   return tmp_voucher;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
VoucherTransitLoad::Voucher
VoucherTransitLoad::yankVoucher(int issuer_rank)
{
   Voucher tmp_voucher(0.0, issuer_rank);
   const_iterator vi = d_vouchers.lower_bound(tmp_voucher);
   if (vi != d_vouchers.end() && vi->d_issuer_rank == issuer_rank) {
      tmp_voucher.d_load = vi->d_load;
      erase(vi);
   }
   return tmp_voucher;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::setTimerPrefix(
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
   d_reserve.setTimerPrefix(timer_prefix + "::redemption");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_adjust_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::adjustLoad()");

   timers.t_raise_dst_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::raiseDstLoad()");

   timers.t_assign_to_local = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::assignToLocal()");

   timers.t_assign_to_local_and_populate_maps = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::assignToLocalAndPopulateMaps()");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
VoucherTransitLoad::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   co << getSumLoad() << " units in " << size() << " vouchers";
   if (detail_depth > 0) {
      size_t count = 0;
      co << ":";
      for (VoucherTransitLoad::const_iterator vi = begin();
           vi != end() && count < 50; ++vi, ++count) {
         co << border << "  " << *vi;
      }
      co << '\n';
   } else {
      co << ".\n";
   }
   co.flush();
}

/*
 *************************************************************************
 * Look for an input database called "VoucherTransitLoad" and read
 * parameters if it exists.
 *************************************************************************
 */
void
VoucherTransitLoad::getFromInput()
{
   if (!tbox::InputManager::inputDatabaseExists()) return;

   std::shared_ptr<tbox::Database> input_db = tbox::InputManager::getInputDatabase();

   if (input_db->isDatabase("VoucherTransitLoad")) {

      std::shared_ptr<tbox::Database> my_db = input_db->getDatabase("VoucherTransitLoad");

      d_partition_work_supply_recursively = my_db->getBoolWithDefault(
            "DEV_partition_work_supply_recursively", d_partition_work_supply_recursively);

      d_print_steps = my_db->getBoolWithDefault(
            "DEV_print_steps", d_print_steps);

      d_print_edge_steps = my_db->getBoolWithDefault(
            "DEV_print_edge_steps", d_print_edge_steps);

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
