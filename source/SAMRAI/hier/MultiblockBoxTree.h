/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of Boxes for overlap searches.
 *
 ************************************************************************/

#ifndef included_hier_MultiblockBoxTree
#define included_hier_MultiblockBoxTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxTree.h"

#include <vector>
#include <map>
#include <memory>

namespace SAMRAI {
namespace hier {

class BaseGridGeometry;
class BoxContainer;

/*!
 * @brief Utility sorting Boxes into tree-like form for finding
 * box overlaps.  Boxes are sorted by BlockId and then a BoxTree is constructed
 * for each block in a BoxContainer.
 *
 * Except for a destructor needed by shared_ptr, the entire interface is
 * private.  This class is intended to be only used by BoxContainer, which
 * is made a friend class.
 */

class MultiblockBoxTree
{

   friend class BoxContainer;

public:
   /*!
    * @brief Destructor.
    */
   ~MultiblockBoxTree() = default;

private:
   /*!
    * @brief Constructs a MultiblockBoxTree from set of Boxes.
    *
    * @param[in] boxes  No empty boxes are allowed.
    *
    * @param[in] grid_geometry GridGeometry associated with boxes in tree.
    *
    * @param[in] min_number Split up sets of boxes while the number of
    * boxes in a subset is greater than this value.  Setting to a
    * larger value tends to make tree building faster but tree
    * searching slower, and vice versa.  @b Default: 10
    *
    * @pre for each box in boxes, box.getBlockId().isValid()
    */
   MultiblockBoxTree(
      const BoxContainer& boxes,
      const BaseGridGeometry* grid_geometry,
      const int min_number = 10);

   /*!
    * Default constructor is unimplemented and should not be used.
    */
   MultiblockBoxTree() = delete;

   /*!
    * @brief Return whether the tree contains any Boxes with the
    * given BlockId.
    */
   bool
   hasBoxInBlock(
      const BlockId& block_id) const
   {
      return d_single_block_trees.find(block_id) !=
             d_single_block_trees.end();
   }

   /*!
    * @brief Return the number of blocks represented in this tree.
    */
   int
   getNumberBlocksInTree() const
   {
      return static_cast<int>(d_single_block_trees.size());
   }

   const BaseGridGeometry *
   getGridGeometry() const
   {
      return d_grid_geometry;
   }

   /*!
    * @brief Reset to uninitialized state.
    *
    * Uninitialized trees can be initialized using generateTree().
    */
   void
   clear()
   {
      d_single_block_trees.clear();
   }

   //@{

   //! @name Overlap checks

   /*
    * @brief Whether the given box has an overlap with Boxes in the
    * tree.
    *
    * @param[in] box The box must have the same BlockId as all Boxes in the
    * tree.
    *
    * @pre getNumberBlocksInTree() == 1
    */
   bool
   hasOverlap(
      const Box& box) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * A pointer to every Box in the tree that intersects with the
    * box argument will be added to the overlap_boxes output vector.  The
    * vector is not sorted in any way.
    *
    * This only works if the tree represents Boxes all having the same BlockId
    * as the argument box.
    *
    * @param[out] overlap_boxes
    *
    * @param[in] box
    *
    * @pre getNumberBlocksInTree() == 1
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * Uses refinement ratio and grid geometry to handle intersections
    * across block boundaries if needed.
    *
    * A pointer to every Box in the tree that intersects with the
    * box argument will be added to the overlap_boxes output vector.  The
    * vector is not sorted in any way.
    *
    * @param[out]  overlap_boxes
    *
    * @param[in]  box
    *
    * @param[in]  refinement_ratio  All boxes in this BoxContainer
    * are assumed to exist in index space that has this refinement ratio
    * relative to the coarse-level domain stored in the grid geometry.
    *
    * @param[in]  include_singularity_block_neighbors  If true, intersections
    * with neighboring blocks that touch only across an enhanced connectivity
    * singularity will be added to output.  If false, those intersections are
    * ignored.
    *
    * @pre (getGridGeometry().getDim() == box.getDim()) &&
    *      (getGridGeometry().getDim() == refinement_ratio.getDim())
    * @pre (box.getBlockId().getBlockValue() >= 0) &&
    *      (box.getBlockId().getBlockValue() < getGridGeometry()->getNumberBlocks())
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * Every Box in this tree that intersects with the box argument
    * will be copied to the overlap_boxes output container.  The output
    * container will retain the same ordered/unordered state that it had
    * prior to being passed into this method.
    *
    * This only works if the tree represents Boxes all having the same BlockId
    * as the argument box.
    *
    * @param[out] overlap_boxes
    *
    * @param[in] box
    *
    * @pre getNumberBlocksInTree() == 1
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * Uses refinement ratio and grid geometry to handle intersections
    * across block boundaries if needed.
    *
    * Every Box in this tree that intersects with the box argument
    * will be copied to the overlap_boxes output container.  The output
    * container will retain the same ordered/unordered state that it had
    * prior to being passed into this method.
    *
    * @param[out]  overlap_boxes
    *
    * @param[in]  box
    *
    * @param[in]  refinement_ratio  All boxes in this BoxContainer
    * are assumed to exist in index space that has this refinement ratio
    * relative to the coarse-level domain stored in the grid geometry.
    *
    * @param[in]  include_singularity_block_neighbors  If true, intersections
    * with neighboring blocks that touch only across an enhanced connectivity
    * singularity will be added to output.  If false, those intersections are
    * ignored.
    *
    * @pre (getGridGeometry().getDim() == box.getDim()) &&
    *      (getGridGeometry().getDim() == refinement_ratio.getDim())
    * @pre (box.getBlockId().getBlockValue() >= 0) &&
    *      (box.getBlockId().getBlockValue() < getGridGeometry()->getNumberBlocks())
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   //@}

private:
   /*!
    * @brief Container of single-block BoxTrees.
    *
    * For each BlockId represented in the tree, there is
    * an entry in this container.
    */
   std::map<BlockId, std::shared_ptr<BoxTree> > d_single_block_trees;

   const BaseGridGeometry* d_grid_geometry;
};

}
}

#endif
