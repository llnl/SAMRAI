/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with overlap Connectors.
 *
 ************************************************************************/
#ifndef included_hier_OverlapConnectorAlgorithm
#define included_hier_OverlapConnectorAlgorithm

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BaseConnectorAlgorithm.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <map>
#include <set>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Algorithms for working Connectors whose neighbor data
 * represents overlaps.
 *
 * An overlap Connector is one in which neighbors represent a pair of
 * overlapping Boxes.  If a base Box grown by the
 * Connector width overlaps a head Box, the head Box is a
 * neighbor of the base Box.  This class implements some
 * functions for working with a overlap Connectors.
 *
 * OverlapConnectorAlgorithm objects create, check and operate on overlap
 * Connectors.
 */
class OverlapConnectorAlgorithm:public BaseConnectorAlgorithm
{

public:
   /*!
    * @brief Constructor.
    */
   OverlapConnectorAlgorithm();

   /*!
    * @brief Destructor.
    */
   virtual ~OverlapConnectorAlgorithm();

   /*!
    * @brief Read extra debugging flag from input database.
    */
   void
   getFromInput();

   /*!
    * @brief Create overlap Connector then discover and add overlaps from base
    * to head to it.
    *
    * The Connector's neighbor information is modified.
    *
    * If the Connector's head is not GLOBALIZED, a copy is made and
    * globalized.  Once a globalized head is obtained, this method
    * simply calls findOverlaps(const BoxLevel &globalized_head).
    *
    * @param[in,out] connector
    * @param[in] base_box_level
    * @param[in] head_box_level
    * @param[in] base_width
    * @param[in] parallel_state
    * @param[in] ignore_self_overlap
    */
   void
   findOverlaps(
      std::shared_ptr<Connector>& connector,
      const BoxLevel& base_box_level,
      const BoxLevel& head_box_level,
      const IntVector& base_width,
      const BoxLevel::ParallelState parallel_state = BoxLevel::DISTRIBUTED,
      const bool ignore_self_overlap = false) const;

   /*!
    * @brief Create overlap Connector with its transpose then discover and add
    * overlaps from base to head to it and overlaps from head to base to
    * transpose.
    *
    * The Connector's neighbor information is modified.
    *
    * If the Connector's head is not GLOBALIZED, a copy is made and
    * globalized.  Once a globalized head is obtained, this method
    * simply calls findOverlaps(const BoxLevel &globalized_head).
    *
    * @param[in,out] connector
    * @param[in] base_box_level
    * @param[in] head_box_level
    * @param[in] base_width
    * @param[in] transpose_base_width
    * @param[in] parallel_state
    * @param[in] ignore_self_overlap
    */
   void
   findOverlapsWithTranspose(
      std::shared_ptr<Connector>& connector,
      const BoxLevel& base_box_level,
      const BoxLevel& head_box_level,
      const IntVector& base_width,
      const IntVector& transpose_base_width,
      const BoxLevel::ParallelState parallel_state = BoxLevel::DISTRIBUTED,
      const bool ignore_self_overlap = false) const;

   /*!
    * @brief Discover and add overlaps from base to head for an
    * overlap Connector.
    *
    * The Connector's neighbor information is modified.
    *
    * If the Connector's head is not GLOBALIZED, a copy is made and
    * globalized.  Once a globalized head is obtained, this method
    * simply calls findOverlaps(const BoxLevel &globalized_head).
    *
    * @param[in,out] connector
    * @param[in] ignore_self_overlap
    */
   void
   findOverlaps(
      Connector& connector,
      const bool ignore_self_overlap = false) const;

   /*
    * @brief Discover and add overlaps from base to
    * (a globalized version of the) head for an overlap Connector.
    *
    * A Box's overlap with another Box is disregarded if:
    * @li @c ignore_self_overlap is true and
    * @li the two Boxes have the same BoxId and.
    * @li the Connectors's base and head have the same refinement ratio.
    *
    * @param[in,out] connector
    * @param[in] globalized_head
    * @param[in] ignore_self_overlap
    */
   void
   findOverlaps(
      Connector& connector,
      const BoxLevel& globalized_head,
      const bool ignore_self_overlap = false) const;

   /*
    * @brief Populate overlap relationships for the given Connector by
    * using the assumed partition algorithm to find overlaps.
    *
    * For the assumed partition algorithm, see Allison Baker's paper.
    */
   void
   findOverlaps_assumedPartition(
      Connector& connector) const;

   /*!
    * @brief For a given Connector, get the subset of overlapping neighbors
    * defined by the given Connector width.
    *
    * The difference between extractNeighbors() and
    * Connector::getNeighbors() is that extractNeighbors() extracts
    * the subset of the currently stored neighbors associated with the
    * given Connector width.  getNeighbors() returns the entire
    * neighbor set, which corresponds to the Connector's width.  The
    * specified Connector width must be less than or equal to the
    * Connector's width.
    *
    * @param[out] neighbors
    * @param[in] connector
    * @param[in] box_id
    * @param[in] width
    *
    * @pre width <= connector.getConnectorWidth()
    * @pre (connector.getParallelState() == BoxLevel::GLOBALIZED) ||
    *      (box_id.getOwnerRank() == connector.getMPI().getRank())
    * @pre connector.getBase().hasBox(box_id)
    */
   void
   extractNeighbors(
      Connector::NeighborSet& neighbors,
      const Connector& connector,
      const BoxId& box_id,
      const IntVector& width) const;

   /*!
    * @brief Like extractNeighbors above except that it computes all
    * neighborhoods of connector placing them into other.
    *
    * @param[out] other
    * @param[in] connector
    * @param[in] width
    *
    * @pre width <= connector.getConnectorWidth()
    * @pre for the box_id of each neighborhood base box in connector,
    *      (connector.getParallelState() == BoxLevel::GLOBALIZED) ||
    *      (box_id.getOwnerRank() == connector.getMPI().getRank())
    * @pre for the box_id of each neighborhood base box in connector,
    *      connector.getBase().hasBox(box_id)
    */
   void
   extractNeighbors(
      Connector& other,
      const Connector& connector,
      const IntVector& width) const;

   /*!
    * @brief Compute the overlap Connectors between BoxLevels
    * efficiently by using information from existing overlap
    * Connectors.
    *
    * Let east, west and center be BoxLevels.  Given
    * overlap Connectors between center and east and between center and
    * west, compute overlap Connectors between east and west.
    *
    * @code
    *
    *               Input:                            Output:
    *
    *                                               west to east
    *      (west)              (east)      (west) ---------------> (east)
    *        \ ^                ^ /               <---------------
    *         \ \center  center/ /                  east to west
    *          \ \ to     to  / /
    *       west\ \west  east/ /east
    *        to  \ \        / /  to
    *       center\ \      / /center
    *              v \    / v
    *               (center)
    *
    * @endcode
    *
    * The "bridge" is the Connectors between @c west and @c east.
    * Bridging is an algorithm for finding their overlaps using
    * existing overlaps incident from @c center.  This is more
    * efficient than findOverlaps() because it does not require
    * acquiring, storing or searching globalized data.
    *
    * Preconditions:
    *
    * - Four input Connectors refer to the center BoxLevel, two
    * refer to the west, and two refer to the east.  The center, east
    * and west BoxLevels must be the same regardless of what
    * Connector is used to get them.  For example, <tt>
    * &west_to_center.getTranspose().getHead() == &west_to_center.getBase() </tt> must
    * be true.
    *
    * Postconditions:
    *
    * @li @c west and @c east, as referenced by the output Connectors
    * are heads of west_to_center's transpose and @c center_to_east.
    *
    * @li Widths of the output Connectors will be either @c
    * center_to_east's width reduced by @c center_growth_to_nest_west
    * or west_to_center's transpose's width reduced by @c
    * center_growth_to_nest_east, which ever is greater.  Output
    * Connector widths are still limited by the @c
    * connector_width_limit argument.  All comparisons are done after
    * the appropriate index space conversions, of course.
    *
    * @note
    *
    * @li Although @c east, @c west and @c center are not explicit
    * parameters, they are specified implicitly through the input
    * Connectors.
    *
    * @li The bridge operation works in the degenerate case where @c
    * east and @c west are the same object.  In that case, @c
    * west_to_east and its transpose should also be the same object.
    *
    * @li Bridging finds as many overlaps as it can, given the inputs,
    * but can only guarantee that the output Connectors are complete
    * overlap Connectors if there is nesting as specified by the
    * arguments @c center_growth_to_nest_west and/or @c
    * center_growth_to_nest_east.  Nesting is not checked, because
    * checking requires an iterative non-scalable computation.
    * Nesting is best determined by the code that created the
    * BoxLevels involved.  Bridging is meant to be fast and
    * scalable.
    *
    * @li Periodic relationships are automatically generated if a head
    * Box in its shifted position overlaps the grown base Box.
    * However, because neighbors may be remote, we cannot scalably
    * verify that the periodic neighbors actually exist in the heads
    * of the output Connectors.  Spurious periodic relationships may
    * be generated when the periodic image simply has not be added to
    * the head or when parts of east or west lie outside the domain
    * extents.  They can be discarded using
    * Connector::removePeriodicRelationships and properly regenerated
    * using BoxLevelConnectorUtils.
    *
    * @param[out] west_to_east
    * @param[in] west_to_center
    * @param[in] center_to_east
    *
    * @param center_growth_to_nest_west The amount by which the center
    * BoxLevel must grow to nest the west BoxLevel.
    * Bridging guarantees completeness if the width of @b
    * center_to_east exceeds this amount.  If unknown, set to negative
    * value if unknown so it won't be considered when computing the
    * output Connector widths.
    *
    * @param center_growth_to_nest_east The amount by which the center
    * BoxLevel must grow to nest the east BoxLevel.
    * Bridging guarantees completeness if the width of @b
    * west_to_center's transpose exceeds this amount.  If unknown, set to
    * negative value if unknown so it won't be considered when computing the
    * output Connector widths.
    *
    * @param connector_width_limit specifies the maximum Connector
    * width to compute overlaps for.  The connector_width should be in
    * the coarser of the east and west indices.  If
    * connector_width_limit is negative, use the default
    * connector_width, which is the larger of the west_to_center's transpose
    * and center_to_east Connectors', coarsened into the coarser of east
    * and west indices.
    *
    * @param compute_transpose true if west_to_east's transpose should be
    * computed
    *
    * @pre west_to_cent.hasTranspose()
    * @pre cent_to_east.hasTranspose()
    */
   void
   bridgeWithNesting(
      std::shared_ptr<Connector>& west_to_east,
      const Connector& west_to_center,
      const Connector& center_to_east,
      const IntVector& center_growth_to_nest_west,
      const IntVector& center_growth_to_nest_east,
      const IntVector& connector_width_limit,
      bool compute_transpose) const;

   /*!
    * @brief A version of bridge without any guarantee of nesting.
    *
    * The east and west BoxLevels are assumed
    * to nest in the center BoxLevel.  If they do not,
    * the results are not guaranteed to be complete.
    *
    * The output Connector widths are the greater of the widths of @c
    * center_to_east and @c center_to_west (converted into the proper
    * index space, of course).
    *
    * @param[out] west_to_east
    * @param[in] west_to_center
    * @param[in] center_to_east
    * @param[in] connector_width_limit
    * @param compute_transpose true if west_to_east's transpose should be
    * computed
    *
    * @pre west_to_cent.hasTranspose()
    * @pre cent_to_east.hasTranspose()
    *
    * @see bridgeWithNesting( Connector& west_to_east, const Connector& west_to_center, const Connector& center_to_east, const IntVector& center_growth_to_nest_west, const IntVector& center_growth_to_nest_east, const IntVector& connector_width_limit, bool compute_transpose) const;
    */
   void
   bridge(
      std::shared_ptr<Connector>& west_to_east,
      const Connector& west_to_center,
      const Connector& center_to_east,
      const IntVector& connector_width_limit,
      bool compute_transpose) const;

   /*!
    * @brief A version of bridge without limiting the connector_width of the
    * result.
    *
    * @param[out] west_to_east
    * @param[in] west_to_center
    * @param[in] center_to_east
    * @param compute_transpose true if west_to_east's transpose should be
    * computed
    *
    * @pre west_to_cent.hasTranspose()
    * @pre cent_to_east.hasTranspose()
    *
    * @see bridge( Connector& west_to_east, const Connector& west_to_center, const Connector& center_to_east, const IntVector& connector_width_limit, bool compute_transpose) const;
    */
   void
   bridge(
      std::shared_ptr<Connector>& west_to_east,
      const Connector& west_to_center,
      const Connector& center_to_east,
      bool compute_transpose) const;

   /*!
    * @brief A version of bridge without any guarantee of nesting in which
    * an input connector and its transpose are modified to form the resulting
    * bridge connectors.
    *
    * The east and west BoxLevels are assumed
    * to nest in the center BoxLevel.  If they do not,
    * the results are not guaranteed to be complete.
    *
    * The output Connector widths are the greater of the widths of @c
    * center_to_east and west_to_center's transpose (converted into the proper
    * index space, of course).
    *
    * @param[in,out] west_to_center
    * @param[in] center_to_east
    * @param[in] connector_width_limit
    *
    * @pre west_to_cent.hasTranspose()
    * @pre cent_to_east.hasTranspose()
    *
    * @see bridgeWithNesting( Connector& west_to_east, const Connector& west_to_center, const Connector& center_to_east, const IntVector& center_growth_to_nest_west, const IntVector& center_growth_to_nest_east, const IntVector& connector_width_limit, bool compute_transpose) const;
    */
   void
   bridge(
      Connector& west_to_center,
      const Connector& center_to_east,
      const IntVector& connector_width_limit) const;

   /*!
    * @brief Set whether to barrier before potential major
    * communication.
    *
    * This developer feature makes sure all processes start major
    * operations at the same time so that timers do not include the
    * time waiting for slower processes to get to the starting point.
    *
    * @param[in] do_barrier
    */
   void
   setBarrierBeforeCommunication(
      bool do_barrier)
   {
      d_barrier_before_communication = do_barrier;
   }

   /*!
    * @brief When @c do_check is true, turn on sanity checks for input
    * parameters.
    *
    * @note
    * Sanity checks occur at the beginning of certain methods only.
    * The checks are expensive and meant mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPreconditions(
      bool do_check)
   {
      d_sanity_check_method_preconditions = do_check;
   }

   /*!
    * @brief When @c do_check is true, turn on sanity checks for outputs
    *
    * @note
    * Sanity checks occur at the end of certain methods only.
    * The checks are expensive and meant mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPostconditions(
      bool do_check)
   {
      d_sanity_check_method_postconditions = do_check;
   }

   /*!
    * @brief Set the SAMRAI_MPI to use.
    *
    * If set, communication will use the specified SAMRAI_MPI instead
    * of the SAMRAI_MPI from BoxLevels.  This protects communication
    * operations from accidentally interacting with unrelated
    * communications, but it limits operations to work only with
    * metadata objects with comptatible (congruent) SAMRAI_MPI
    * objects.
    *
    * If make_duplicate is true, the specified SAMRAI_MPI will be
    * duplicated for exclusise use.  The duplicate will be freed upon
    * object destruction.
    */
   void
   setSAMRAI_MPI(
      const tbox::SAMRAI_MPI& mpi,
      bool make_duplicate = true);

   /*!
    * @brief When @c print_steps is true, print what the code is
    * doing.
    *
    * @note Step printing may be expensive and and is meant mainly for
    * debugging.
    *
    * @param[in] print_steps
    */
   void
   setPrintSteps(
      bool print_steps)
   {
      d_print_steps = print_steps;
   }

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named
    * "hier::OverlapConnectorAlgorithm::*", where the third field is
    * the specific steps performed by the MappingConnectorAlgorithm.
    * You can override the first two fields with this method.
    * Conforming to the timer naming convention, timer_prefix should
    * have the form "*::*".
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

   /*!
    * @brief Get the name of this object.
    */
   const std::string
   getObjectName() const
   {
      return "OverlapConnectorAlgorithm";
   }

private:
   // Internal shorthand.
   typedef Connector::NeighborSet NeighborSet;

   void
   privateBridge_prologue(
      const Connector& west_to_cent,
      const Connector& cent_to_east,
      const Connector& east_to_cent,
      const Connector& cent_to_west,
      bool west_nesting_is_known,
      const IntVector& cent_growth_to_nest_west,
      bool east_nesting_is_known,
      const IntVector& cent_growth_to_nest_east,
      const IntVector& connector_width_limit,
      bool compute_transpose,
      IntVector& west_to_east_width,
      IntVector& east_to_west_width,
      std::set<int>& incoming_ranks,
      std::set<int>& outgoing_ranks,
      NeighborSet& visible_west_nabrs,
      NeighborSet& visible_east_nabrs) const;

   /*!
    * @brief This is where the bridge algorithm is implemented.
    * All public bridge interfaces which do not perform an "in place" bridge
    * call this method underneath.
    *
    * @param west_to_east
    *
    * @param east_to_west
    *
    * @param cent_to_east
    *
    * @param compute_transpose true if east_to_west should be computed
    *
    * @param incoming_ranks
    *
    * @param outgoing_ranks
    *
    * @param visible_west_nabrs
    *
    * @param visible_east_nabrs
    */
   void
   privateBridge(
      Connector& west_to_east,
      Connector* east_to_west,
      const Connector& cent_to_east,
      bool compute_transpose,
      const std::set<int>& incoming_ranks,
      const std::set<int>& outgoing_ranks,
      NeighborSet& visible_west_nabrs,
      NeighborSet& visible_east_nabrs) const;

   /*
    * @brief Perform checks on the arguments of bridge.
    */
   void
   privateBridge_checkParameters(
      const Connector& west_to_cent,
      const Connector& cent_to_east,
      const Connector& east_to_cent,
      const Connector& cent_to_west) const;

   /*!
    * @brief Relationship removal part of overlap algorithm, caching
    * outgoing information in message buffers.
    */
   void
   privateBridge_removeAndCache(
      std::map<int, std::vector<int> >& send_mesgs,
      Connector& overlap_connector,
      Connector* overlap_connector_transpose,
      const Connector& misc_connector) const;

   /*!
    * @brief Find all relationships in the Connector(s) to be computed and send
    * outgoing information.
    */
   void
   privateBridge_discoverAndSend(
      std::map<int, std::vector<int> >& send_mesgs,
      Connector& west_to_east,
      Connector* east_to_west,
      const std::set<int>& incoming_ranks,
      const std::set<int>& outgoing_ranks,
      tbox::AsyncCommPeer<int>* all_comms,
      NeighborSet& visible_west_nabrs,
      NeighborSet& visible_east_nabrs) const;

   /*!
    * @brief Find all relationships in the Connector(s) to be computed.
    */
   void
   privateBridge_discover(
      std::vector<int>& send_mesg,
      Connector& west_to_east,
      Connector* east_to_west,
      const NeighborSet& visible_west_nabrs,
      const NeighborSet& visible_east_nabrs,
      NeighborSet::const_iterator& west_ni,
      NeighborSet::const_iterator& east_ni,
      int curr_owner,
      const BoxContainer& east_rbbt,
      const BoxContainer& west_rbbt,
      const tbox::Dimension& dim,
      bool compute_transpose,
      int rank) const;

   /*!
    * @brief Find overlap and save in bridging connector or pack
    * into send message, used in privateBridge().
    */
   void
   privateBridge_findOverlapsForOneProcess(
      const int curr_owner,
      const NeighborSet& visible_base_nabrs,
      NeighborSet::const_iterator& base_ni,
      std::vector<int>& send_mesg,
      const int remote_box_counter_index,
      Connector& bridging_connector,
      NeighborSet& referenced_head_nabrs,
      const BoxContainer& head_rbbt) const;

   /*!
    * @brief Utility used in privateBridge()
    */
   void
   privateBridge_unshiftOverlappingNeighbors(
      const Box& box,
      BoxContainer& neighbors,
      BoxContainer& scratch_space,
      const IntVector& neighbor_refinement_ratio,
      const PeriodicShiftCatalog& shift_catalog) const;

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   //! @brief SAMRAI_MPI for internal communications.
   tbox::SAMRAI_MPI d_mpi;

   //! @brief Whether d_mpi was duplicated for exclusive use.
   bool d_mpi_is_exclusive;

   // Extra checks independent of optimization/debug.
   static char s_print_steps;

   /*!
    * @brief Tag to use (and increment) at begining of operations that
    * require nearest-neighbor communication, to aid in eliminating
    * mixing of messages from different internal operations.
    *
    * Using unique tags was important to the early development of this
    * class but may not be necessary anymore.
    */
   static int s_operation_mpi_tag;


   //@{
   //! @name Timer data for this class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      std::shared_ptr<tbox::Timer> t_find_overlaps_rbbt;
      std::shared_ptr<tbox::Timer> t_find_overlaps_assumed_partition;
      std::shared_ptr<tbox::Timer> t_find_overlaps_assumed_partition_get_ap;
      std::shared_ptr<tbox::Timer> t_find_overlaps_assumed_partition_connect_to_ap;
      std::shared_ptr<tbox::Timer> t_find_overlaps_assumed_partition_transpose;

      std::shared_ptr<tbox::Timer> t_bridge;
      std::shared_ptr<tbox::Timer> t_bridge_setup_comm;
      std::shared_ptr<tbox::Timer> t_bridge_remove_and_cache;
      std::shared_ptr<tbox::Timer> t_bridge_discover_and_send;
      std::shared_ptr<tbox::Timer> t_bridge_discover_get_neighbors;
      std::shared_ptr<tbox::Timer> t_bridge_discover_form_rbbt;
      std::shared_ptr<tbox::Timer> t_bridge_share;
      std::shared_ptr<tbox::Timer> t_bridge_receive_and_unpack;
      std::shared_ptr<tbox::Timer> t_bridge_MPI_wait;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   static char s_ignore_external_timer_prefix;

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

   bool d_print_steps;
   bool d_barrier_before_communication;
   bool d_sanity_check_method_preconditions;
   bool d_sanity_check_method_postconditions;

   static tbox::StartupShutdownManager::Handler s_initialize_finalize_handler;

};

}
}

#endif // included_hier_OverlapConnectorAlgorithm
