/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_VoucherTransitLoad
#define included_mesh_VoucherTransitLoad

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/SequentialLocalIdGenerator.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <set>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Implementation of TransitLoad, representing the load with a
 * set of vouchers, each of which has a work value and the rank of the
 * process that issued the voucher.
 *
 * As a container, this class is nearly identical to
 * std::set<Voucher,VoucherMoreLoad>.
 *
 * Terminology: To avoid confusing the sending and receiving of
 * messages, the sending and receiving of work use the terms supply
 * and demand.  During redemption of vouchers, the holder of a voucher
 * demands the work.  The issuer of that voucher supplies the work.
 * Both send and receive messages to accomplish this.
 *
 * @internal
 * There are 2 private nested classes.  Voucher encapsulates voucher info.
 * VoucherRedemption encapsulates the process of redeeming a voucher.
 */

class VoucherTransitLoad:public TransitLoad
{

public:
   typedef double LoadType;

   //! @name Constructor
   VoucherTransitLoad(
      const PartitioningParams& pparams);

   //! @name Copy constructor
   VoucherTransitLoad(
      const VoucherTransitLoad& other,
      bool copy_load = true);

   //@{
   //! @name TransitLoad abstract interfaces

   //! @copydoc TransitLoad::clone()
   VoucherTransitLoad *
   clone() const;

   //! @copydoc TransitLoad::initialize()
   void
   initialize();

   //! @copydoc TransitLoad::getSumLoad()
   LoadType getSumLoad() const {
      return d_sumload;
   }

   //! @copydoc TransitLoad::insertAll( const hier::BoxContainer & )
   void
   insertAll(
      const hier::BoxContainer& box_container);

   //! @copydoc TransitLoad::insertAllWithArtificialMinimum(
   //!           const hier::BoxContainer&, double)
   void
   insertAllWithArtificialMinimum(
      const hier::BoxContainer& box_container,
      double minimum_load);

   //! @copydoc TransitLoad::insertAll( TransitLoad & )
   void
   insertAll(
      TransitLoad& other);

   //! @copydoc TransitLoad::insertAllWithExisitngLoads( const hier::BoxContainer & )
   void
   insertAllWithExistingLoads(
      const hier::BoxContainer& box_container)
   {
      clear();
      insertAll(box_container);
   }

   //! @copydoc TransitLoad::setWorkload( const hier::PatchLevel&, const int )
   void
   setWorkload(
      const hier::PatchLevel& patch_level,
      const int work_data_id);

   //! @copydoc TransitLoad::getNumberOfItems()
   size_t
   getNumberOfItems() const;

   //! @copydoc TransitLoad::getNumberOfOriginatingProcesses()
   size_t
   getNumberOfOriginatingProcesses() const;

   //! @copydoc TransitLoad::putToMessageStream()
   void
   putToMessageStream(
      tbox::MessageStream& msg) const;

   //! @copydoc TransitLoad::getFromMessageStream()
   void
   getFromMessageStream(
      tbox::MessageStream& msg);

   /*!
    * @copydoc TransitLoad::adjustLoad()
    */
   LoadType
   adjustLoad(
      TransitLoad& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load);

   /*!
    * @copydoc TransitLoad::assignToLocal()
    *
    * This method uses communication to redeem vouchers.
    */
   void
   assignToLocal(
      hier::BoxLevel& balanced_box_level,
      const hier::BoxLevel& unbalanced_box_level,
      double flexible_load_tol = 0.0,
      const tbox::SAMRAI_MPI& alt_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL));

   /*!
    * @copydoc TransitLoad::assignToLocalAndPopulateMaps()
    *
    * This method uses communication to redeem vouchers and up the map.
    */
   void
   assignToLocalAndPopulateMaps(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector& balanced_to_unbalanced,
      hier::MappingConnector& unbalanced_to_balanced,
      double flexible_load_tol = 0.0,
      const tbox::SAMRAI_MPI& alt_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL));

   //@}

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named "mesh::VoucherTransitLoad::*",
    * where the third field is the specific steps performed
    * by the Schedule.  You can override the first two
    * fields with this method.  Conforming to the timer
    * naming convention, timer_prefix should have the form
    * "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   void
   recursivePrint(
      std::ostream& co = tbox::plog,
      const std::string& border = std::string(),
      int detail_depth = 1) const;

private:
   /*!
    * @brief Real implementation of assignToLocalAndPopulateMaps, with
    * slightly modified interface.
    */
   void
   assignToLocalAndPopulateMaps(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector* balanced_to_unbalanced,
      hier::MappingConnector* unbalanced_to_balanced,
      const hier::BoxLevel& unbalanced_box_level,
      double flexible_load_tol,
      const tbox::SAMRAI_MPI& alt_mpi);

   //! @brief MPI tag for demand communication.
   static const int VoucherTransitLoad_DEMANDTAG = 3;
   //! @brief MPI tag for supply communication.
   static const int VoucherTransitLoad_SUPPLYTAG = 4;

   //! @brief Voucher.
   struct Voucher {
      //! @brief Default constructor sets zero value and invalid issuer.
      Voucher():
         d_issuer_rank(tbox::SAMRAI_MPI::getInvalidRank()),
         d_load(0.0) {
      }
      //! @brief Initializing constructor sets voucher value and issuer.
      Voucher(
         const LoadType& load,
         int issuer_rank):
         d_issuer_rank(issuer_rank),
         d_load(load) {
      }
      //! @brief Construct Voucher by combining two vouchers from the same issuer.
      Voucher(
         const Voucher& a,
         const Voucher& b):
         d_issuer_rank(a.d_issuer_rank),
         d_load(a.d_load + b.d_load) {
         if (a.d_issuer_rank != b.d_issuer_rank) {
            TBOX_ERROR("VoucherTransitLoad: Cannot combine vouchers from different issuers.");
         }
      }
      //@ @brief Construct Voucher by taking value from an existing Voucher.
      Voucher(
         LoadType load,
         Voucher& src):
         d_issuer_rank(src.d_issuer_rank),
         d_load(load <= src.d_load ? load : src.d_load) {
         src.d_load -= d_load;
      }
      friend std::ostream& operator << (std::ostream& co, const Voucher& v) {
         co << v.d_issuer_rank << '|' << v.d_load;
         return co;
      }
      friend tbox::MessageStream& operator << (tbox::MessageStream& ms, const Voucher& v) {
         ms << v.d_issuer_rank << v.d_load;
         return ms;
      }
      friend tbox::MessageStream& operator >> (tbox::MessageStream& ms, Voucher& v) {
         ms >> v.d_issuer_rank >> v.d_load;
         return ms;
      }
      /*!
       * @brief Adjust load by taking work from or giving work to
       * another Voucher.
       *
       * Similar to the interface defined in TransitLoad but working
       * with individual vouchers instead of containers.
       */
      LoadType
      adjustLoad(
         Voucher& other,
         LoadType ideal_load);
      int d_issuer_rank;
      LoadType d_load;
   };

   //! @brief Comparison functor for sorting Vouchers by issuer rank and load.
   struct VoucherRankCompare {
      bool operator () (const Voucher& a, const Voucher& b) const {
         if (a.d_issuer_rank != b.d_issuer_rank) {
            return a.d_issuer_rank < b.d_issuer_rank;
         }
         return a.d_load < b.d_load;
      }
   };

   //@{
   //! @name Interfaces like the C++ standard stl::set, to help readability.
   typedef std::set<Voucher, VoucherRankCompare> RankOrderedVouchers;
   typedef RankOrderedVouchers::iterator iterator;
   typedef RankOrderedVouchers::const_iterator const_iterator;
   typedef RankOrderedVouchers::reverse_iterator reverse_iterator;
   typedef RankOrderedVouchers::key_type key_type;
   typedef RankOrderedVouchers::value_type value_type;
   iterator begin() {
      return d_vouchers.begin();
   }
   iterator end() {
      return d_vouchers.end();
   }
   const_iterator begin() const {
      return d_vouchers.begin();
   }
   const_iterator end() const {
      return d_vouchers.end();
   }
   reverse_iterator rbegin() {
      return d_vouchers.rbegin();
   }
   reverse_iterator rend() {
      return d_vouchers.rend();
   }
   size_t size() const {
      return d_vouchers.size();
   }
   bool empty() const {
      return d_vouchers.empty();
   }
   void clear() {
      d_sumload = 0;
      d_vouchers.clear();
   }
   std::pair<iterator, bool> insert(const Voucher& v) {
      TBOX_ASSERT(v.d_load >= d_pparams->getLoadComparisonTol());
      iterator itr = d_vouchers.lower_bound(Voucher(0, v.d_issuer_rank));
      if (itr != d_vouchers.end() &&
          itr->d_issuer_rank == v.d_issuer_rank) {
         TBOX_ERROR(
            "Cannot insert Voucher " << v
                                     << ".\nExisting voucher " << *itr
                                     << " is from the same issuer."
                                     << "\nTo combine the vouchers, use insertCombine().");
      }
      itr = d_vouchers.insert(itr, v);
      d_sumload += v.d_load;
      checkDupes();
      return std::pair<iterator, bool>(itr, true);
   }
   size_t erase(const Voucher& v) {
      iterator vi = d_vouchers.lower_bound(v);
      if (vi != d_vouchers.end()) {
         d_sumload -= vi->d_load;
         erase(vi);
         return 1;
      }
      return 0;
   }
   void erase(iterator pos) {
      d_sumload -= pos->d_load;
#ifdef DEBUG_CHECK_ASSERTIONS
      size_t old_size = d_vouchers.size();
#endif
      d_vouchers.erase(pos);
      TBOX_ASSERT(d_vouchers.size() == old_size - 1);
   }
   void swap(VoucherTransitLoad& other) {
      const LoadType tmpload = d_sumload;
      d_sumload = other.d_sumload;
      other.d_sumload = tmpload;
      d_vouchers.swap(other.d_vouchers);
   }
   //@}

   //! @brief Encapsulates voucher redemption for both demander and supplier.
   struct VoucherRedemption {
      VoucherRedemption():
         d_demander_rank(tbox::SAMRAI_MPI::getInvalidRank()),
         d_pparams(0),
         d_mpi(MPI_COMM_NULL),
         d_mpi_request(MPI_REQUEST_NULL) {
      }
      ~VoucherRedemption() {
         finishSendRequest();
         d_pparams = 0;
      }

      //@{
      //! @name Demanding and supplying work based on a voucher.
      void
      sendWorkDemand(
         const VoucherTransitLoad::const_iterator& voucher,
         const hier::SequentialLocalIdGenerator& id_gen,
         const tbox::SAMRAI_MPI& mpi);

      void
      recvWorkDemand(
         int demander_rank,
         int message_length,
         const tbox::SAMRAI_MPI& mpi);

      void
      sendWorkSupply(
         BoxTransitSet& reserve,
         double flexible_load_tol,
         const PartitioningParams& pparams,
         bool send_all);

      void
      recvWorkSupply(
         int message_length,
         const PartitioningParams& pparams);

      void
      setLocalRedemption(
         const VoucherTransitLoad::const_iterator& voucher,
         const hier::SequentialLocalIdGenerator& id_gen,
         const tbox::SAMRAI_MPI& mpi);

      void
      fulfillLocalRedemption(
         BoxTransitSet& reserve,
         const PartitioningParams& pparams,
         bool all);
      //@}

      void
      takeWorkFromReserve(
         BoxTransitSet& work,
         BoxTransitSet& reserve);

      void
      finishSendRequest();

      Voucher d_voucher;
      int d_demander_rank;
      //! @brief Demander-specified LocalId generator to avoid ID clashes.
      hier::SequentialLocalIdGenerator d_id_gen;
      //! @brief Shipment of work, as boxes, sent or received.
      std::shared_ptr<BoxTransitSet> d_box_shipment;
      const PartitioningParams* d_pparams;

      std::shared_ptr<tbox::MessageStream> d_msg;
      tbox::SAMRAI_MPI d_mpi;
      tbox::SAMRAI_MPI::Request d_mpi_request;
   };

   /*!
    * @brief Insert voucher, combining with existing voucher from same issuer.
    */
   iterator insertCombine(const Voucher& v) {
      iterator itr = d_vouchers.lower_bound(Voucher(0.0, v.d_issuer_rank));
      if (itr == d_vouchers.end() ||
          v.d_issuer_rank != itr->d_issuer_rank) {
         // Safe to insert.
         TBOX_ASSERT(v.d_load >= d_pparams->getLoadComparisonTol());
         itr = d_vouchers.insert(itr, v);
      } else {
         // Create combined voucher to replace existing one.
         Voucher combined(*itr, v);
         TBOX_ASSERT(combined.d_load >= d_pparams->getLoadComparisonTol());
         d_vouchers.erase(itr++);
         itr = d_vouchers.insert(itr, combined);
      }
      d_sumload += v.d_load;
      checkDupes();
      return itr;
   }

   /*!
    * @brief Erase voucher issued by the given process.
    *
    * @return Whether there was a Voucher to be erased.
    */
   bool
   eraseIssuer(
      int issuer_rank);

   //! @brief Sanity check enforcing no-duplicate-issuer rule.
   void checkDupes() const {
      for (const_iterator i = begin(); i != end(); ++i) {
         const_iterator i1 = i;
         ++i1;
         if (i1 != end() && i1->d_issuer_rank == i->d_issuer_rank) {
         }
      }
   }

   //! @brief Sanity check catching extremely small vouchers.
   void checkSmallVouchers() const {
      for (const_iterator i = begin(); i != end(); ++i) {
         if (i->d_load < d_pparams->getLoadComparisonTol()) {
            TBOX_ERROR("Voucher " << *i << " is smaller than tolerance "
                                  << d_pparams->getLoadComparisonTol());
         }
      }
   }

   /*!
    * @brief Raise load of dst container by shifing load from src.
    */
   LoadType
   raiseDstLoad(
      VoucherTransitLoad& src,
      VoucherTransitLoad& dst,
      LoadType ideal_dst_load);

   /*!
    * @brief Assign a reserve to a set of VoucherRedemption.
    *
    * Alternative option to recursively partition work supply.
    *
    * On return, work assignments will be reflected in reserve.
    */
   void
   recursiveSendWorkSupply(
      const std::map<int, VoucherRedemption>::iterator& begin,
      const std::map<int, VoucherRedemption>::iterator& end,
      BoxTransitSet& reserve);

   /*!
    * @brief Re-cast a TransitLoad object to a VoucherTransitLoad.
    */
   const VoucherTransitLoad& recastTransitLoad(const TransitLoad& transit_load) {
      const VoucherTransitLoad* ptr = static_cast<const VoucherTransitLoad *>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }

   /*!
    * @brief Re-cast a TransitLoad object to a VoucherTransitLoad.
    */
   VoucherTransitLoad& recastTransitLoad(TransitLoad& transit_load) {
      VoucherTransitLoad* ptr = static_cast<VoucherTransitLoad *>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }

   /*!
    * @brief Return the Voucher issued by the given process.  If
    * Voucher isn't there, return zero-valued Voucher.
    */
   Voucher
   findVoucher(
      int issuer_rank) const;

   /*!
    * @brief Yank out and return the Voucher issued by the given
    * process.  If Voucher isn't there, return zero-valued Voucher.
    */
   Voucher
   yankVoucher(
      int issuer_rank);

   /*!
    * @brief Look for an input database called "VoucherTransitLoad"
    * and read parameters if it exists.
    */
   void
   getFromInput();

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void initializeCallback() {
      TimerStruct& timers(s_static_timers[s_default_timer_prefix]);
      getAllTimers(s_default_timer_prefix, timers);
   }

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void finalizeCallback() {
      s_static_timers.clear();
   }

   //! @brief Work load, sorted by issuer rank.
   RankOrderedVouchers d_vouchers;

   LoadType d_sumload;

   const PartitioningParams* d_pparams;

   bool d_partition_work_supply_recursively;

   double d_flexible_load_tol;

   //! @brief Reserve load container used during redemption phase.
   BoxTransitSet d_reserve;

   //@{
   //! @name Debugging stuff.
   bool d_print_steps;
   bool d_print_edge_steps;
   //@}

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
      std::shared_ptr<tbox::Timer> t_adjust_load;
      std::shared_ptr<tbox::Timer> t_raise_dst_load;
      std::shared_ptr<tbox::Timer> t_assign_to_local;
      std::shared_ptr<tbox::Timer> t_assign_to_local_and_populate_maps;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

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

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};

}
}

#endif
