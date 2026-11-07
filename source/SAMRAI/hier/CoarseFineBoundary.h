/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   For describing coarse-fine boundary interfaces
 *
 ************************************************************************/

#ifndef included_hier_CoarseFineBoundary
#define included_hier_CoarseFineBoundary

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <vector>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utility class to construct and maintain a description of the
 * coarse-fine boundary between a patch level and a coarser level.
 *
 * A coarse-fine boundary box is a BoundaryBox object, but it is generated
 * differently than a typical boundary box maintained by a patch geometry
 * object.  A boundary box serving as a coarse-fine boundary box describes part
 * of the boundary of a given patch with its next coarser AMR hierarchy level.
 * It does not intersect any other patch on the same level, nor does it lie
 * on a physical domain boundary, except where the physical boundary is
 * periodic and the appropriate continuation of that boundary is part of a
 * coarser patch level.
 *
 * The coarse-fine boundary is typically created from two adjacent
 * hierarchy levels, but the description lives on (refers to the index
 * space of) the finer level.  Since the coarse-fine boundary
 * describes the boundary to the next coarser level, the coarsest
 * level (level zero) has no coarse-fine boundary.
 *
 * Each CoarseFineBoundary object corresponds to one level, so to
 * represent a entire hierarchy, one would need an array or list of
 * such objects.
 */

class CoarseFineBoundary
{
public:
   /*!
    * @brief Construct a CoarseFineBoundary object with no boundary boxes.
    *
    * @param[in] dim Dimension
    */
   explicit CoarseFineBoundary(
      const tbox::Dimension& dim);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] rhs
    */
   CoarseFineBoundary(
      const CoarseFineBoundary& rhs);

   /*!
    * @brief Construct a CoarseFineBoundary object for the specified
    * level in the given patch hierarchy.
    *
    * @note If level number is zero, the coarse-fine boundary will be empty.
    *
    * @param[in] hierarchy
    * @param[in] level_num
    * @param[in] max_ghost_width The ghost width determines the extent of the
    *                            boundary boxes along the level domain boundary,
    *                            similar to regular domain boundary boxes.  Note
    *                            that as in the case of regular boundary boxes,
    *                            each box will always be one cell wide in the
    *                            direction perpendicular to the patch boundary.
    *
    * @pre max_ghost_width > IntVector(max_ghost_width.getDim(), -1)
    */
   CoarseFineBoundary(
      const PatchHierarchy& hierarchy,
      int level_num,
      const IntVector& max_ghost_width);

   /*!
    * @brief Construct a CoarseFineBoundary object for a specified level.
    *
    * The coarse-fine boundary will be computed using the physical domain
    * as the reference coarser level.  The physical domain is provided to this
    * method as the 'head' level of the box_level_to_domain Connector.
    *
    * @note If the level covers the entire physical domain, the coarse-fine
    * boundary will be empty.
    *
    * @param[in] level
    * @param[in] box_level_to_domain
    * @param[in] box_level_to_self
    * @param[in] max_ghost_width The ghost width determines the extent of the
    *                            boundary boxes along the level domain boundary,
    *                            similar to regular domain boundary boxes.  Note
    *                            that as in the case of regular boundary boxes,
    *                            each box will always be one cell wide in the
    *                            direction perpendicular to the patch boundary.
    *
    * @pre max_ghost_width > IntVector(max_ghost_width.getDim(), -1)
    */
   CoarseFineBoundary(
      const PatchLevel& level,
      const Connector& box_level_to_domain,
      const Connector& box_level_to_self,
      const IntVector& max_ghost_width);

   /*!
    * @brief Destructor.
    */
   ~CoarseFineBoundary() = default;

   /*!
    * @brief Clear all boundary data.
    */
   void
   clear()
   {
      d_boundary_boxes.clear();
   }

   //@{
   /*!
    * @name Functions to get the computed coarse-fine boundaries.
    */

   /*!
    * @brief Get a vector of boundary boxes of a given type
    * for a specified patch.
    *
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param[in] global_id
    * @param[in] boundary_type Codimension of boundaries.
    * @param[in] block_id     Defaults to 0 for the single block case
    *
    * @pre d_initialize[block_id.getBlockValue()]
    */
   const std::vector<BoundaryBox>&
   getBoundaries(
      const GlobalId& global_id,
      const int boundary_type,
      const BlockId& block_id = BlockId::zero()) const;

   /*!
    * @brief Get a vector of node boundary boxes for a specified patch.
    *
    * @see BoundaryBox for more information.
    *
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param[in] global_id
    * @param[in] block_id     Defaults to 0 for the single block case
    */
   const std::vector<BoundaryBox>&
   getNodeBoundaries(
      const GlobalId& global_id,
      const BlockId& block_id = BlockId::zero()) const
   {
      return getBoundaries(global_id, getDim().getValue(), block_id);
   }

   /*!
    * @brief Get a vector of edge boundary boxes for a specified patch.
    *
    * @see BoundaryBox for more information.
    *
    * Note that edge boxes are only meaningful if the dimension is > 1.
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param[in] global_id
    * @param[in] block_id     Defaults to 0 for the single block case
    *
    * @pre getDim().getValue() >= 2
    */
   const std::vector<BoundaryBox>&
   getEdgeBoundaries(
      const GlobalId& global_id,
      const BlockId& block_id = BlockId::zero()) const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (getDim().getValue() < 2) {
         TBOX_ERROR("CoarseFineBoundary::getEdgeBoundaries():  There are\n"
            << "no edge boundaries in " << d_dim << "d.\n");
      }
#endif
      return getBoundaries(global_id, getDim().getValue() - 1, block_id);
   }

   /*!
    * @brief Get a vector of face boundary boxes for a specified patch.
    *
    * @see BoundaryBox for more information.
    *
    * Note that face boxes are only meaningful if the dimension is > 2.
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param[in] global_id
    * @param[in] block_id     Defaults to 0 for the single block case
    *
    * @pre getDim().getValue() >= 3
    */
   const std::vector<BoundaryBox>&
   getFaceBoundaries(
      const GlobalId& global_id,
      const BlockId& block_id = BlockId::zero()) const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (getDim().getValue() < 3) {
         TBOX_ERROR("CoarseFineBoundary::getFaceBoundaries():  There are\n"
            << "no face boundaries in " << d_dim << "d.\n");
      }
#endif
      return getBoundaries(global_id, getDim().getValue() - 2, block_id);
   }

   //@}

   /*!
    * @brief Print out class data.
    *
    * @param[in] os Output stream
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    */
   CoarseFineBoundary&
   operator = (
      const CoarseFineBoundary& rhs)
   {
      d_initialized = rhs.d_initialized;
      d_boundary_boxes = rhs.d_boundary_boxes;
      return *this;
   }

   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   /* Don't allow default ctor */
   CoarseFineBoundary();

   /*!
    * @brief Compute a CoarseFineBoundary object for a specified level.
    *
    * The coarse-fine boundary will be computed using the physical domain
    * as the reference coarser level.  The physical domain is provided to this
    * method as the 'head' level of the box_level_to_domain Connector.
    *
    * @note If the level covers the entire physical domain, the coarse-fine
    * boundary will be empty.
    *
    * @param[in] level
    * @param[in] level_to_domain Connector from level to physical domain level
    * @param[in] level_to_level  Connector from level to itself
    * @param[in] max_ghost_width The ghost width determines the extent of the
    *                            boundary boxes along the level domain boundary,
    *                            similar to regular domain boundary boxes.  Note
    *                            that as in the case of regular boundary boxes,
    *                            each box will always be one cell wide in the
    *                            direction perpendicular to the patch boundary.
    *
    * @pre getDim() == max_ghost_width.getDim()
    */
   void
   computeFromLevel(
      const PatchLevel& level,
      const Connector& level_to_domain,
      const Connector& level_to_level,
      const IntVector& max_ghost_width);

   /*!
    * @brief Compute a CoarseFineBoundary object for a specified level.
    *
    * The coarse-fine boundary will be computed using the physical domain
    * as the reference coarser level.
    *
    * @note If the level covers the entire physical domain, the coarse-fine
    * boundary will be empty.
    *
    * @param[in] level
    * @param[in] level_to_domain Connector from level to physical domain level
    * @param[in] level_to_level  Connector from level to itself
    * @param[in] max_ghost_width The ghost width determines the extent of the
    *                            boundary boxes along the level domain boundary,
    *                            similar to regular domain boundary boxes.  Note
    *                            that as in the case of regular boundary boxes,
    *                            each box will always be one cell wide in the
    *                            direction perpendicular to the patch boundary.
    *
    * @pre getDim() == max_ghost_width.getDim()
    */
   void
   computeFromMultiblockLevel(
      const PatchLevel& level,
      const Connector& level_to_domain,
      const Connector& level_to_level,
      const IntVector& max_ghost_width);

   /*!
    * @brief Dimension of the object.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Whether the boundary boxes have been computed.  One entry on the
    * vector for each block in a multiblock grid.
    */
   std::vector<bool> d_initialized;

   /*!
    * @brief Patch boundary boxes describing the coarse-fine boundary.
    *
    * Map each local patch to a PatchBoundaries object.
    */
   std::map<BoxId, PatchBoundaries> d_boundary_boxes;

};

}
}

#endif
