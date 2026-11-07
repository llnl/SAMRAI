/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_BoxTransitSet
#define included_mesh_BoxTransitSet

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/SequentialLocalIdGenerator.h"
#include "SAMRAI/mesh/BalanceBoxBreaker.h"
#include "SAMRAI/mesh/BoxInTransit.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <set>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Implementation of TransitLoad, representing the load with a
 * set of boxes, each of which represents a load and knows the origin
 * of its load.
 *
 * As a container, this class is nearly identical to
 * std::set<BoxInTransit,BoxInTransitMoreLoad>.
 */

class BoxTransitSet:public TransitLoad
{

public:
   typedef double LoadType;

   //! @name Constructor
   BoxTransitSet(
      const PartitioningParams& pparams);

   /*!
    * @name Copy constructor
    *
    * The content may be omitted from the copy, using the flag copy_load.
    */
   BoxTransitSet(
      const BoxTransitSet& other,
      bool copy_load = true);

   //@{
   //! @name TransitLoad abstract interfaces

   //! @copydoc TransitLoad::clone()
   BoxTransitSet *
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

   /*! @brief Insert all boxes while applying a minimum load value
    *
    * The minimum load value is assigned to any boxes whose computed load
    * is smaller than the given value.  This allows for boxes with small
    * true loads to be treated as having a larger fixed load value, if that
    * behavior is desired by the calling code.
    *
    * @param other  box container with boxes to inserted into this object.
    * @param minimum_load  artificial minimum load value
    */
   void
   insertAllWithArtificialMinimum(
      const hier::BoxContainer& box_container,
      double minimum_load);

   //! @copydoc TransitLoad::insertAll( TransitLoad & )
   void
   insertAll(
      TransitLoad& other);

   //! @copydoc TransitLoad::insertAllWithExistingLoads( const hier::BoxContainer & )
   void
   insertAllWithExistingLoads(
      const hier::BoxContainer& box_container);

   /*!
    * @brief Set workloads in the members of this BoxTransitSet.
    *
    * The BoxInTransit members of this BoxTransitSet will have their
    * workload values set according to the data on the given patch level
    * associated with the given data id
    *
    * @param[in]  patch_level  Level with boxes that are the same as those
    *                          held in this BoxTransitSet
    * @param[in]  work_data_id Data id for the workload data that exists on
    *                          the patch level
    */
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
    * This method uses communication to populate the map.
    */
   void
   assignToLocalAndPopulateMaps(
      hier::BoxLevel& balanced_box_level,
      hier::MappingConnector& balanced_to_unbalanced,
      hier::MappingConnector& unbalanced_to_balanced,
      double flexible_load_tol = 0.0,
      const tbox::SAMRAI_MPI& alt_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL));

   //@}

   /*
    * @brief Reassign the boxes to the new owner.
    *
    * Any box that isn't already owned by the new owner or doesn't
    * have a valid LocalId, is given one by the
    * SequentialLocalIdGenerator.
    */
   void
   reassignOwnership(
      hier::SequentialLocalIdGenerator& id_gen,
      int new_owner_rank);

   /*!
    * @brief Put local Boxes into a BoxLevel.
    */
   void
   putInBoxLevel(
      hier::BoxLevel& box_level) const;

   /*!
    * @brief Generate unbalanced<==>balanced edges incident from
    * local boxes.
    *
    * This method does no communication.  Semilocal edges not incident
    * from remote Boxes are not communicated.
    */
   void
   generateLocalBasedMapEdges(
      hier::MappingConnector& unbalanced_to_balanced,
      hier::MappingConnector& balanced_to_unbalanced) const;

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named "mesh::BoxTransitSet::*",
    * where the third field is the specific steps performed
    * by the Schedule.  You can override the first two
    * fields with this method.  Conforming to the timer
    * naming convention, timer_prefix should have the form
    * "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   /*!
    * @brief Set print flags for individual object.
    */
   void
   setPrintFlags(
      bool steps,
      bool pop_steps,
      bool swap_steps,
      bool break_steps,
      bool edge_steps);

   void
   recursivePrint(
      std::ostream& co = tbox::plog,
      const std::string& border = std::string(),
      int detail_depth = 1) const;

   /*!
    * @brief Intermediary between BoxTransitSet and output streams,
    * adding ability to control the output.  See
    * BoxTransitSet::format().
    */
   class Outputter
   {

      //! @brief Insert a BoxTransitSet to the stream according to Outputter settings.
      friend std::ostream&
      operator << (std::ostream& s,
                   const Outputter& f)
      {
         f.d_boxes.recursivePrint(s, f.d_border, f.d_detail_depth);
         return s;
      }

private:
      friend class BoxTransitSet;
      /*!
       * @brief Copy constructor
       */
      Outputter(
         const Outputter& other);

      /*!
       * @brief Construct the Outputter with a BoxTransitSet and the
       * parameters needed to output the BoxTransitSet to a stream.
       */
      Outputter(
         const BoxTransitSet& boxes,
         const std::string& border,
         int detail_depth = 2):
         d_boxes(boxes),
         d_border(border),
         d_detail_depth(detail_depth)
      {
      }

      void
      operator = (
         const Outputter& rhs);               // Unimplemented private.
      const BoxTransitSet& d_boxes;
      const std::string d_border;
      const int d_detail_depth;
   };

   /*!
    * @brief Return a object to that can format the BoxTransitSet for
    * inserting into output streams.
    *
    * Usage example (printing with a tab indentation):
    * @verbatim
    *    cout << "my boxes:\n" << boxes.format("\t") << endl;
    * @endverbatim
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 2) const
   {
      return Outputter(*this, border, detail_depth);
   }

   //@}

private:
   void
   populateMaps(
      hier::MappingConnector& balanced_to_unbalanced,
      hier::MappingConnector& unbalanced_to_balanced,
      const tbox::SAMRAI_MPI& alt_mpi) const;

   static const int BoxTransitSet_EDGETAG0 = 3;
   static const int BoxTransitSet_EDGETAG1 = 4;
   static const int BoxTransitSet_FIRSTDATALEN = 1000;

   /*!
    * @brief Comparison functor for sorting BoxInTransit from more to
    * less loads.
    *
    * Ties are broken by BlockId, then lexical comparison of the box's
    * lower corner, then lexical comparison of the upper corner, then
    * orginator BoxId.  The comparison should not use the box's BoxId
    * because some boxes may not have valid ones.
    */
   struct BoxInTransitMoreLoad {
      bool operator () (
         const BoxInTransit& a,
         const BoxInTransit& b) const {
         if (tbox::MathUtilities<double>::Abs(a.getLoad() - b.getLoad()) > 1.0e-10) {
            return a.getLoad() > b.getLoad();
         }
         if (a.getBox().getBlockId() != b.getBox().getBlockId()) {
            return a.getBox().getBlockId() < b.getBox().getBlockId();
         }
         if (a.getBox().lower() != b.getBox().lower()) {
            return lexicalIndexLessThan(a.getBox().lower(), b.getBox().lower());
         }
         if (a.getBox().upper() != b.getBox().upper()) {
            return lexicalIndexLessThan(a.getBox().upper(), b.getBox().upper());
         }
         return a.getOrigBox().getBoxId() < b.getOrigBox().getBoxId();
      }
private:
      bool lexicalIndexLessThan(const hier::Index& a,
                                const hier::Index& b) const {
         for (hier::Index::dir_t i = 0; i < a.getDim().getValue(); ++i) {
            if (a(i) != b(i)) return a(i) < b(i);
         }
         return false;
      }
   };

public:
   //@{
   //! @name Interfaces like the C++ standard stl::set, to help readability.
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::iterator iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::const_iterator const_iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::reverse_iterator reverse_iterator;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::key_type key_type;
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad>::value_type value_type;
   iterator begin() {
      return d_set.begin();
   }
   iterator end() {
      return d_set.end();
   }
   const_iterator begin() const {
      return d_set.begin();
   }
   const_iterator end() const {
      return d_set.end();
   }
   reverse_iterator rbegin() const {
      return d_set.rbegin();
   }
   reverse_iterator rend() const {
      return d_set.rend();
   }
   size_t size() const {
      return d_set.size();
   }
   std::pair<iterator, bool> insert(const value_type& x) {
      std::pair<iterator, bool> rval = d_set.insert(x);
      if (rval.second) d_sumload += x.getLoad();
      if (rval.second) d_sumsize += x.getSize();
      return rval;
   }
   void erase(iterator pos) {
      d_sumload -= pos->getLoad();
      d_sumsize -= pos->getSize();
      d_set.erase(pos);
   }
   size_t erase(const key_type& k) {
      const size_t num_erased = d_set.erase(k);
      if (num_erased) d_sumload -= k.getLoad();
      if (num_erased) d_sumsize -= k.getSize();
      return num_erased;
   }
   bool empty() const {
      return d_set.empty();
   }
   void clear() {
      d_sumload = 0;
      d_sumsize = 0;
      d_set.clear();
   }
   void swap(BoxTransitSet& other) {
      const LoadType tl = d_sumload;
      const LoadType ts = d_sumsize;
      d_sumload = other.d_sumload;
      d_sumsize = other.d_sumsize;
      other.d_sumload = tl;
      other.d_sumsize = ts;
      d_set.swap(other.d_set);
   }
   iterator lower_bound(const key_type& k) const {
      return d_set.lower_bound(k);
   }
   iterator upper_bound(const key_type& k) const {
      return d_set.upper_bound(k);
   }
   //@}

private:
   //@{ @name Load adjustment methods

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving the
    * biggest between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadByPopping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load);

   /*!
    * @brief Adjust the load in this BoxTransitSet by swapping boxes
    * between it and another BoxTransitSet.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadBySwapping(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load);

   /*!
    * @brief Adjust the load in this BoxTransitSet by moving work
    * between it and another BoxTransitSet.  One box may be broken
    * up to have a part of its load moved.
    *
    * @param[in,out] hold_bin Holding bin for reserve load.
    *
    * @param[in] ideal_load The load that this bin should have.
    *
    * @param[in] low_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @param[in] high_load Return when this bin's load is in the range
    * [low_load,high_load]
    *
    * @return Net load added to this BoxTransitSet.  If negative, load
    * decreased.
    */
   LoadType
   adjustLoadByBreaking(
      BoxTransitSet& hold_bin,
      LoadType ideal_load,
      LoadType low_load,
      LoadType high_load);

   /*!
    * @brief Find a BoxInTransit in each of the source and destination
    * containers that, when swapped, effects a transfer of the given
    * amount of work from the source to the destination.  Swap the boxes.
    *
    * @param [in,out] src
    *
    * @param [in,out] dst
    *
    * @param actual_transfer [out] Amount of work transfered from src to
    * dst.
    *
    * @param ideal_transfer [in] Amount of work to be transfered from
    * src to dst.
    *
    * @param low_transfer
    *
    * @param high_transfer
    */
   bool
   swapLoadPair(
      BoxTransitSet& src,
      BoxTransitSet& dst,
      LoadType& actual_transfer,
      LoadType ideal_transfer,
      LoadType low_transfer,
      LoadType high_transfer) const;

   //@}

   /*!
    * @brief Communicate semilocal relationships to load donors and
    * fill in unbalanced--->balanced Connectors.
    *
    * These relationships must be represented by this object.
    * Semilocal means the local process owns either Box or
    * the orginal box (not both!) of each item in this BoxTransitSet.
    *
    * Recall that semi-local relationships are those where the base
    * and head boxes are owned by different processes.  These edges
    * require communication to set up.
    *
    * @param [out] unbalanced_to_balanced
    */
   void
   constructSemilocalUnbalancedToBalanced(
      hier::MappingConnector& unbalanced_to_balanced,
      const tbox::SAMRAI_MPI& mpi) const;

   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   const BoxTransitSet& recastTransitLoad(const TransitLoad& transit_load) {
      const BoxTransitSet* ptr = static_cast<const BoxTransitSet *>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }

   /*!
    * @brief Re-cast a TransitLoad object to a BoxTransitSet.
    */
   BoxTransitSet& recastTransitLoad(TransitLoad& transit_load) {
      BoxTransitSet* ptr = static_cast<BoxTransitSet *>(&transit_load);
      TBOX_ASSERT(ptr);
      return *ptr;
   }

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

   //! @brief Compute the load for a Box.
   double computeLoad(const hier::Box& box) const {
      return static_cast<double>(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      return static_cast<double>((box * restriction).size());
   }

   /*!
    * @brief Look for an input database called "BoxTransitSet" and
    * read parameters if it exists.
    */
   void
   getFromInput();

   //! @brief Balance penalty is proportional to imbalance.
   double computeBalancePenalty(double imbalance) const {
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   std::set<BoxInTransit, BoxInTransitMoreLoad> d_set;
   LoadType d_sumload;
   LoadType d_sumsize;

   const PartitioningParams* d_pparams;

   BalanceBoxBreaker d_box_breaker;

   //@{
   //! @name Debugging stuff.
   bool d_print_steps;
   bool d_print_pop_steps;
   bool d_print_swap_steps;
   bool d_print_break_steps;
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
      std::shared_ptr<tbox::Timer> t_adjust_load_by_popping;
      std::shared_ptr<tbox::Timer> t_adjust_load_by_swapping;
      std::shared_ptr<tbox::Timer> t_shift_loads_by_breaking;
      std::shared_ptr<tbox::Timer> t_find_swap_pair;
      std::shared_ptr<tbox::Timer> t_assign_to_local_process_and_populate_maps;
      std::shared_ptr<tbox::Timer> t_populate_maps;
      std::shared_ptr<tbox::Timer> t_construct_semilocal;
      std::shared_ptr<tbox::Timer> t_construct_semilocal_comm_wait;
      std::shared_ptr<tbox::Timer> t_construct_semilocal_send_edges;
      std::shared_ptr<tbox::Timer> t_pack_edge;
      std::shared_ptr<tbox::Timer> t_unpack_edge;
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
