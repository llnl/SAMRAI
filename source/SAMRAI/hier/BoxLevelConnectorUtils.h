/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utilities for working on DLBG edges.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevelConnectorUtils
#define included_hier_BoxLevelConnectorUtils

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/BoxLevel.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for common operating on BoxLevels.
 *
 * Objects of this class can be set to perform certain sanity checks
 * on the pre and post conditions of the methods.  See
 * setSanityCheckMethodPreconditions() and
 * setSanityCheckMethodPostconditions().
 */
class BoxLevelConnectorUtils
{

public:
   /*!
    * @brief Default constructor.
    *
    * By default, sanity checks are disabled.  To enable them, see
    * setSanityCheckMethodPreconditions() and
    * setSanityCheckMethodPostconditions().
    */
   BoxLevelConnectorUtils();

   /*!
    * @brief Destructor.
    */
   ~BoxLevelConnectorUtils() = default;

   /*!
    * @brief Set whether to run expensive sanity checks on input parameters.
    *
    * Mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPreconditions(
      bool do_check)
   {
      d_sanity_check_precond = do_check;
   }

   /*!
    * @brief Set whether to run expensive sanity checks on output parameters.
    *
    * Mainly for debugging.
    *
    * @param[in] do_check
    */
   void
   setSanityCheckMethodPostconditions(
      bool do_check)
   {
      d_sanity_check_postcond = do_check;
   }

   //@{

   //! @name Comparing boxes of two BoxLevels

   /*!
    * @brief Given an overlap Connector, determine the extent to which
    * the Connector's base nests in its head.
    *
    * This method returns true if the base, grown by @c base_swell,
    * nests inside the head, grown by @c head_swell, by a margin of @c
    * head_nesting_margin.  @c base_swell and @c head_swell should be
    * non-negative and specified in the base and head index spaces,
    * respectively.  @c head_nesting_margin should be in the head
    * index space.
    *
    * The Connector width must be at least the sum of the @c
    * base_swell and the (appropriately converted) @c head_swell and
    * @c head_nesting_margin.  We require and assume without verifying
    * that the Connector is complete.
    *
    * If the domain is given, non-nesting parts outside of the domain
    * are disregarded.
    *
    * @param[out] locally_nests Whether the local parts of the base
    * nests in the head.  This output may vary among the processes.
    *
    * @param[in] connector
    *
    * @param[in] base_swell the amount that the base is grown by, given in the
    * base index space and non-negative
    *
    * @param[in] head_swell the amount that the head is grown by, given in the
    * head index space and non-negative
    *
    * @param[in] head_nesting_margin given in the head index space.
    *
    * @param[in] domain Domain description, in reference index space,
    * in search tree format.
    *
    * @return True if the given base BoxLevel nests in the head,
    * otherwise False.
    * @pre (connector.getBase().getDim() == base_swell.getDim()) &&
    *      (connector.getBase().getDim() == head_nesting_margin.getDim())
    * @pre connector.isFinalized()
    * @pre base_swell >= IntVector::getZero(connector.getBase().getDim())
    * @pre head_nesting_margin >= IntVector::getZero(connector.getBase().getDim())
    */
   bool
   baseNestsInHead(
      bool* locally_nests,
      const Connector& connector,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_nesting_margin,
      const BoxContainer* domain = 0) const;

   /*!
    * @brief Given base and head BoxLevels, determine the extent
    * to which the base nests in the head.
    *
    * This method is similar to the version taking a Connector instead
    * of the base and head BoxLevels, except that it will use
    * the base BoxLevel's PersistentOverlapConnectors object to
    * get the base--->head Connector.  If such a Connector does not
    * exist, the PersistentOverlapConnectors object will create it, an
    * unscalable operation possibly requiring collective
    * communication.
    *
    * @param domain Domain description, in reference index space, in
    * search tree format.
    *
    * @return Whether the given base BoxLevel nests in the head.
    *
    * @param[out] locally_nests Whether the local parts of the base
    * nests in the head.  This output may vary among the processes.
    *
    * @param[in] base
    *
    * @param[in] head
    *
    * @param[in] base_swell the amount that the base is grown by, given in the
    * base index space and non-negative
    *
    * @param[in] head_swell the amount that the head is grown by, given in the
    * head index space and non-negative
    *
    * @param[in] head_margin given in the head index space.
    *
    * @param[in] domain Domain description, in reference index space,
    * in search tree format.
    *
    * @return Whether the given base BoxLevel nests in the head.
    *
    * @pre head_nesting_margin.getDim() == base_swell.getDim()
    * @pre head.getMPI() == base.getMPI()
    * @pre base_swell >= IntVector::getZero(head.getDim())
    * @pre head_swell >= IntVector::getZero(head.getDim())
    * @pre head_nesting_margin >= IntVector::getZero(head.getDim())
    * @pre (head.getRefinementRatio() <= base.getRefinementRatio()) ||
    *      (head.getRefinementRatio() >= base.getRefinementRatio()) ||
    */
   bool
   baseNestsInHead(
      bool* locally_nests,
      const BoxLevel& base,
      const BoxLevel& head,
      const IntVector& base_swell,
      const IntVector& head_swell,
      const IntVector& head_margin,
      const BoxContainer* domain = 0) const;

   /*!
    * @brief Compute the parts of one BoxLevel that are external
    * to another BoxLevel.
    *
    * Compare an input BoxLevel to a "reference" BoxLevel.
    * Compute the parts of the input that are external to the
    * reference.  Build the "external" BoxLevel representing the
    * external parts.  Build a MappingConnector with the input as its
    * base and the external as its head.
    *
    * A partially external input cell (possible when input is coarser
    * than reference) is considered to be external.
    *
    * For the purpose of defining what is external, the reference
    * level can be grown by nesting_width before comparing.  This
    * feature can be used to determing which parts of the input does
    * not nest in the reference comparison.  A negative growth
    * indicates shrinking of the reference level.
    *
    * This method does not require any communication.
    *
    * @param[out] external  The existing state will be discarded.
    *
    * @param[out] input_to_external  The existing state will be
    * discarded.
    *
    * @param[in] input_to_reference Overlap Connector from input to
    * reference BoxLevel.  The width of input_to_reference must be
    * at least one and at least the absolute value of nesting_width.
    *
    * @param[in] nesting_width Growth of the reference BoxLevel for
    * the purpose of comparing to input.  Must be in resolution of
    * input BoxLevel.  Must be either non-negative or non-positive but
    * not mixed.  If any width is negative, then input_to_reference
    * must have a Connector width of at least 1.  Otherwise, an error
    * is thrown (because correct results cannot be guaranteed).
    *
    * @param[in] domain The domain representation, without periodic
    * images, in search tree form.  These boxes should be in the
    * reference index space.  If domain is given, do not shrink the
    * reference BoxLevel where it touches the domain boundary.
    *
    * @post input_to_external.isLocal()
    */
   void
   computeExternalParts(
      std::shared_ptr<BoxLevel>& external,
      std::shared_ptr<MappingConnector>& input_to_external,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const BoxContainer& domain = BoxContainer()) const
   {
      d_object_timers->t_compute_external_parts->start();
      computeInternalOrExternalParts(
         external,
         input_to_external,
         'e',
         input_to_reference,
         nesting_width,
         domain);
      d_object_timers->t_compute_external_parts->stop();
   }

   /*!
    * @brief Compute the parts of one BoxLevel that are internal
    * to another BoxLevel.
    *
    * Compare an input BoxLevel to a "reference" box_level.
    * Identify parts of the input that are internal to the reference
    * BoxLevel, and store the internal parts in a
    * BoxLevel.  Set up a MappingConnector between the input
    * and its internal parts.
    *
    * A partially internal input cell (possible when input is coarser
    * than reference) is considered to be internal.
    *
    * For the purpose of defining what is external, the reference
    * level can be grown by nesting_width before comparing.  This
    * feature can be used to determing which parts of the input does
    * not nest in the reference comparison.  A negative growth
    * indicates shrinking of the reference level.
    *
    * This method does not require any communication.
    *
    * @param[out] internal  The existing state will be discarded.
    *
    * @param[out] input_to_internal  The existing state will be
    * discarded.
    *
    * @param[in] input_to_reference Overlap Connector from input to
    * reference BoxLevel.  The width of input_to_reference must be
    * at least one and at least the absolute value of nesting_width.
    *
    * @param[in] nesting_width Growth of the reference BoxLevel for
    * the purpose of comparing to input.  Must be in resolution of
    * input BoxLevel.  Must be either non-negative or non-positive but
    * not mixed.  If any width is negative, then input_to_reference
    * must have a Connector width of at least 1.  Otherwise, an error
    * is thrown (because correct results cannot be guaranteed).
    *
    * @param[in] domain The domain representation, without periodic
    * images, in search tree form.  These boxes should be in the
    * reference index space.  If domain is given, do not shrink the
    * reference BoxLevel where it touches the domain boundary.
    *
    * @post input_to_internal.isLocal()
    */
   void
   computeInternalParts(
      std::shared_ptr<BoxLevel>& internal,
      std::shared_ptr<MappingConnector>& input_to_internal,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const BoxContainer& domain = BoxContainer()) const
   {
      d_object_timers->t_compute_internal_parts->start();
      computeInternalOrExternalParts(
         internal,
         input_to_internal,
         'i',
         input_to_reference,
         nesting_width,
         domain);
      d_object_timers->t_compute_internal_parts->stop();
   }

   /*!
    * @brief Compute parts of a level that do not intersect with another
    *        level
    *
    * Given the Connector input_to_takeaway, this method computes the parts
    * of the base level 'input' that do not intersect the head level
    * 'takeaway'.  The results are stored in the BoxLevel 'remainder',
    * with the input_to_remainder Connector storing the overlap
    * relationships between 'input' and 'remainder'.
    *
    * @param[out] remainder  The non-intersecting parts of input level
    * @param[out] input_to_remainder  Overlap connector with width zero
    * @param[in] input_to_takeaway    Overlap connector, may have any width
    *
    * @pre  input_to_takeaway is a complete overlap connector
    * @post input_to_remainder.isLocal()
    */
   void
   computeNonIntersectingParts(
      std::shared_ptr<BoxLevel>& remainder,
      std::shared_ptr<Connector>& input_to_remainder,
      const Connector& input_to_takeaway) const;

   //@}

   /*!
    * @brief Given a set of Boxes, compute its boundary as a set
    * of boxes located just outside it.
    *
    * @param[out] boundary Boundary boxes, sorted into BoxContainers
    * according to the BlockId.
    *
    * @param[in] refinement_ratio Refinement ratio of boxes.
    *
    * @param[in] grid_geometry
    *
    * @param simplify_boundary_boxes Whether to simplify the boundary
    * boxes after computing them.
    */
   void
   computeBoxesAroundBoundary(
      BoxContainer& boundary,
      const IntVector& refinement_ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geometry,
      const bool simplify_boundary_boxes = true) const;

   //@{

   //! @name Setting up common MappingConnectors

   /*
    * @brief Sort the Boxes in BoxLevel and make a MappingConnector
    * from the unsorted BoxLevel to the sorted one.
    * The sorting can renumber the LocalIndices of the Boxes
    * or put the Boxes in spatial ordering, or both.
    *
    * The MappingConnector created is local (no Box is mapped to a new
    * owner).
    *
    * If @c sort_boxes_by_corner is true, the map will reorder
    * local Boxes by their box corners.  This is useful for
    * making random box ordering into something deterministic.
    *
    * If @c sequentialize_global_indices is true, determine the lowest
    * index the local processor should use for the output
    * BoxLevel so that the global set of Boxes have
    * sequential indices.  (This requires communication.)  If false,
    * each processor start numbering with LocalIndex zero.  If true,
    * @c initial_sequential_index can specify the first index of first
    * Box of the lowest rank processor.
    *
    * For more information on MappingConnectors, see
    * MappingConnectorAlgorithm.
    *
    * @param[out] sorted_box_level Sorted version of the input
    * unsorted_box_level.
    *
    * @param[out] output_map Mapping from @c unsorted_box_level
    * to @c sorted_box_level.
    *
    * @param[in] unsorted_box_level
    *
    * @param[in] sort_boxes_by_corner Whether to sort local
    * Boxes by their indices to make their ordering solely a
    * function of box positions.
    *
    * @param[in] sequentialize_global_indices Whether to renumber the
    * LocalIndices into a globally sequential numbering.
    *
    * @param[in] initial_sequential_index The first index of first
    * Box of the lowest rank process.  This parameter is
    * disregarded when not globally sequentializing the indices.
    */
   void
   makeSortingMap(
      std::shared_ptr<BoxLevel>& sorted_box_level,
      std::shared_ptr<MappingConnector>& output_map,
      const BoxLevel& unsorted_box_level,
      bool sort_boxes_by_corner = true,
      bool sequentialize_global_indices = true,
      LocalId initial_sequential_index = LocalId::getZero()) const;

   /*
    * @brief Given a mapping from an original BoxLevel to parts
    * to be removed (rejected), construct the remainder BoxLevel
    * and the mapping from the original to a remainder.
    *
    * This method does no communication.
    *
    * @see MappingConnectorAlgorithm.
    *
    * @param[out] remainder The new BoxLevel resulting from
    * removing the rejected parts from the original BoxLevel.
    *
    * @param[out] orig_to_remainder The output mapping.  This is a
    * local map.
    *
    * @param[in] orig_to_rejections Mapping from original
    * BoxLevel to its parts that should be be removed.  This
    * must be a local map.
    *
    * @pre orig_to_rejection.isLocal()
    */
   void
   makeRemainderMap(
      std::shared_ptr<BoxLevel>& remainder,
      std::shared_ptr<MappingConnector>& orig_to_remainder,
      const MappingConnector& orig_to_rejections) const;

   //@}

   //@{

   //! @name Adding periodic images

   /*!
    * @brief Add periodic images to a BoxLevel.
    *
    * This method is a no-op in the case of non-periodic domains.
    *
    * Any periodic image within a certain distance of the domain is
    * added, Those farther out are not added.  The threshold distance
    * is @c threshold_distance.
    *
    * @param[in,out] box_level BoxLevel subject to the
    * addition of periodic Boxes.
    *
    * @param[in] domain_search_tree Domain description in the reference
    * index space.  This tree must NOT include periodic images.
    *
    * @param[in] threshold_distance
    */
   void
   addPeriodicImages(
      BoxLevel& box_level,
      const BoxContainer& domain_search_tree,
      const IntVector& threshold_distance) const;

   /*!
    * @brief Add periodic images to a BoxLevel and add new
    * relationships to the periodic images.
    *
    * This method is a no-op in the case of non-periodic domains.
    *
    * Any periodic images within a certain distance of the domain is
    * added, but the rest are not added.  The threshold distance is
    * the width of the Connector @c box_level_to_anchor.
    *
    * This method updates the overlap Connectors between the
    * BoxLevel getting new periodic Boxes and an "anchor"
    * BoxLevel.  New periodic overlap relationships generated
    * are added to the overlap Connector @c
    * box_level_to_anchor.  If you don't need to have a
    * Connector updated, use addPeriodicImages() instead of this
    * method.
    *
    * Preconditions: box_level<==>anchor must be transpose
    * overlap Connectors (or communications may hang).
    * anchor--->anchor must be a complete overlap Connector.
    *
    * @param[in,out] box_level BoxLevel subject to the
    * addition of periodic Boxes.
    *
    * @param[in,out] box_level_to_anchor Overlap Connector to
    * be updated with new relationships.
    *
    * @param[in] domain_search_tree Domain description in the
    * reference index space.  This tree must NOT include periodic
    * images.
    *
    * @param[in] anchor_to_anchor Self overlap Connector for anchor
    * BoxLevel.  Must be a complete overlap Connector with periodic
    * relationships.  To guarantee complete periodic relationships,
    * the anchor, grown by the width of this Connector, must nest
    * box_level grown by the width of box_level_to_anchor.
    */
   void
   addPeriodicImagesAndRelationships(
      BoxLevel& box_level,
      Connector& box_level_to_anchor,
      const BoxContainer& domain_search_tree,
      const Connector& anchor_to_anchor) const;

   //@}

   /*!
    * @brief Setup names of timers.
    *
    * By default, timers are named
    * "hier::BoxLevelConnectorUtils::*", where the third field is
    * the specific steps performed by the BoxLevelConnectorUtils.
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
      return "BoxLevelConnectorUtils";
   }

private:
   /*!
    * @brief Delegated work of computeInternalParts and computeExternalParts.
    *
    * @post input_to_parts.isLocal()
    */
   void
   computeInternalOrExternalParts(
      std::shared_ptr<BoxLevel>& parts,
      std::shared_ptr<MappingConnector>& input_to_parts,
      char internal_or_external,
      const Connector& input_to_reference,
      const IntVector& nesting_width,
      const BoxContainer& domain) const;

   /*!
    * @brief Call-back function to sort boxes.
    */
   static int
   qsortBoxCompare(
      const void* v,
      const void* w);

   /*!
    * @brief Read extra debugging flag from input database.
    */
   void
   getFromInput();

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

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
      std::shared_ptr<tbox::Timer> t_make_sorting_map;
      std::shared_ptr<tbox::Timer> t_compute_boxes_around_boundary;
      std::shared_ptr<tbox::Timer> t_compute_boxes_around_boundary_singularity;
      std::shared_ptr<tbox::Timer> t_compute_boxes_around_boundary_simplify;
      std::shared_ptr<tbox::Timer> t_compute_external_parts;
      std::shared_ptr<tbox::Timer> t_compute_internal_parts;
      std::shared_ptr<tbox::Timer> t_compute_internal_or_external_parts;
      std::shared_ptr<tbox::Timer> t_compute_internal_or_external_parts_manip_reference;
      std::shared_ptr<tbox::Timer> t_compute_internal_or_external_parts_simplify;
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

   bool d_sanity_check_precond;
   bool d_sanity_check_postcond;

   static tbox::StartupShutdownManager::Handler s_initialize_handler;

};

}
}

#endif  // included_hier_BoxLevelConnectorUtils
