/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Multiblock binary trees of Boxes for overlap searches.
 *
 ************************************************************************/
#include "SAMRAI/hier/MultiblockBoxTree.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BaseGridGeometry.h"


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 * Constructor
 *************************************************************************
 */
MultiblockBoxTree::MultiblockBoxTree(
   const BoxContainer& boxes,
   const BaseGridGeometry* grid_geometry,
   int min_number):
   d_grid_geometry(grid_geometry)
{
   /*
    * Group Boxes by their BlockId and
    * create a tree for each BlockId.
    */
   std::map<BlockId, std::list<const Box *> > single_block_boxes;
   for (BoxContainer::const_iterator bi = boxes.begin();
        bi != boxes.end(); ++bi) {
      TBOX_ASSERT((*bi).getBlockId().isValid());
      const BlockId& block_id = (*bi).getBlockId();
      single_block_boxes[block_id].push_back(&(*bi));
   }

   for (std::map<BlockId, std::list<const Box *> >::iterator blocki =
           single_block_boxes.begin();
        blocki != single_block_boxes.end(); ++blocki) {

      d_single_block_trees[blocki->first].reset(new BoxTree(blocki->second,
            min_number));
   }
}


/*
 *************************************************************************
 * Tell if any Box in the tree intersects the given Box
 *************************************************************************
 */
bool
MultiblockBoxTree::hasOverlap(
   const Box& box) const
{
   if (getNumberBlocksInTree() != 1) {
      TBOX_ERROR("Single block version of hasOverlap called on search tree with multiple blocks.");
   }

   if (hasBoxInBlock(box.getBlockId())) {
      return d_single_block_trees.begin()->second->hasOverlap(box);
   } else {
      return false;
   }
}

/*
 **************************************************************************
 * Fills the vector with pointers to Boxes that intersect the arguement
 **************************************************************************
 */
void
MultiblockBoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*d_grid_geometry, box, refinement_ratio);

   const BlockId& block_id = box.getBlockId();
   TBOX_ASSERT(block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, std::shared_ptr<BoxTree> >::const_iterator blocki(d_single_block_trees.find(
                                                                            block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second->findOverlapBoxes(overlap_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   for (BaseGridGeometry::ConstNeighborIterator ni =
           d_grid_geometry->begin(block_id);
        ni != d_grid_geometry->end(block_id); ++ni) {

      const BaseGridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second->findOverlapBoxes(overlap_boxes, transformed_box);

   }
}

/*
 **************************************************************************
 * Fills the container with pointers to Boxes that intersect the arguement
 **************************************************************************
 */
void
MultiblockBoxTree::findOverlapBoxes(
   BoxContainer& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*d_grid_geometry, box, refinement_ratio);

   const BlockId& block_id = box.getBlockId();
   TBOX_ASSERT(block_id.getBlockValue() < d_grid_geometry->getNumberBlocks());

   /*
    * Search in the index space of block_id for overlaps.
    */

   std::map<BlockId, std::shared_ptr<BoxTree> >::const_iterator blocki(d_single_block_trees.find(
                                                                            block_id));

   if (blocki != d_single_block_trees.end()) {
      blocki->second->findOverlapBoxes(overlap_boxes, box);
   }

   /*
    * Search in the index spaces neighboring block_id for overlaps.
    */

   for (BaseGridGeometry::ConstNeighborIterator ni =
           d_grid_geometry->begin(block_id);
        ni != d_grid_geometry->end(block_id); ++ni) {

      const BaseGridGeometry::Neighbor& neighbor(*ni);

      if (!include_singularity_block_neighbors && neighbor.isSingularity()) {
         continue;
      }

      const BlockId neighbor_block_id(neighbor.getBlockId());

      blocki = d_single_block_trees.find(neighbor_block_id);

      if (blocki == d_single_block_trees.end()) {
         continue;
      }

      Box transformed_box(box);

      d_grid_geometry->transformBox(transformed_box,
         refinement_ratio,
         neighbor_block_id,
         block_id);

      blocki->second->findOverlapBoxes(overlap_boxes, transformed_box);

   }
}

/*
 **************************************************************************
 * Fills the container with pointers to Boxes that intersect the arguement
 **************************************************************************
 */
void
MultiblockBoxTree::findOverlapBoxes(
   BoxContainer& overlap_boxes,
   const Box& box) const
{
   if (getNumberBlocksInTree() != 1) {
      TBOX_ERROR(
         "Single block version of findOverlapBoxes called on search tree with multiple blocks.");
   }

   if (hasBoxInBlock(box.getBlockId())) {
      d_single_block_trees.begin()->second->findOverlapBoxes(
         overlap_boxes,
         box);
   }
}

/*
 **************************************************************************
 * Fills the vector with pointers to Boxes that intersect the arguement
 **************************************************************************
 */
void
MultiblockBoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_boxes,
   const Box& box) const
{
   if (getNumberBlocksInTree() != 1) {
      TBOX_ERROR(
         "Single block version of findOverlapBoxes called on search tree with multiple blocks.");
   }

   if (hasBoxInBlock(box.getBlockId())) {
      d_single_block_trees.begin()->second->findOverlapBoxes(
         overlap_boxes,
         box);
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
