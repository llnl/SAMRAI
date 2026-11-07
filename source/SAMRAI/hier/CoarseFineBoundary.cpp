/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   For describing coarse-fine boundary interfaces
 *
 ************************************************************************/
#include "SAMRAI/hier/CoarseFineBoundary.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/PatchLevel.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

CoarseFineBoundary::CoarseFineBoundary(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_initialized(1, false)
{
}

CoarseFineBoundary::CoarseFineBoundary(
   const CoarseFineBoundary& rhs):
   d_dim(rhs.d_dim),
   d_initialized(1, false),
   d_boundary_boxes(rhs.d_boundary_boxes)
{
   /*
    * This needs to be written this way since STL vector for bools
    * causes uninitialized memory reads because it is poorly implemented.
    * So the vector is initialized and then copied.
    */
   d_initialized = rhs.d_initialized;
}

CoarseFineBoundary::CoarseFineBoundary(
   const PatchHierarchy& hierarchy,
   int level_num,
   const IntVector& max_ghost_width):
   d_dim(max_ghost_width.getDim()),
   d_initialized(hierarchy.getGridGeometry()->getNumberBlocks(), false)
{
   TBOX_ASSERT(max_ghost_width > IntVector(d_dim, -1));

   size_t number_blocks = hierarchy.getGridGeometry()->getNumberBlocks();
   const PatchLevel& level = *hierarchy.getPatchLevel(level_num);
   const IntVector& ratio_to_zero =
      level.getRatioToLevelZero();

   IntVector connector_width(max_ghost_width, number_blocks);
   connector_width.max(IntVector::getOne(d_dim));
   const Connector& level_to_level = level.findConnector(level,
         connector_width,
         CONNECTOR_CREATE);
   if (level_num != 0) {
      for (BlockId::block_t b = 0; b < connector_width.getNumBlocks(); ++b) {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (connector_width(b,d) % ratio_to_zero(b,d) != 0) {
               connector_width(b,d) = 
                  (connector_width(b,d) / ratio_to_zero(b,d)) +
                  ratio_to_zero(b,d);
            }
         }
      }
   }

   const Connector& level_to_domain =
      level.getBoxLevel()->findConnector(hierarchy.getDomainBoxLevel(),
         connector_width,
         CONNECTOR_CREATE);

   if (hierarchy.getGridGeometry()->getNumberBlocks() == 1) {
      computeFromLevel(
         level,
         level_to_domain,
         level_to_level,
         max_ghost_width);
   } else {
      computeFromMultiblockLevel(
         level,
         level_to_domain,
         level_to_level,
         max_ghost_width);
   }

}

CoarseFineBoundary::CoarseFineBoundary(
   const PatchLevel& level,
   const Connector& level_to_domain,
   const Connector& level_to_level,
   const IntVector& max_ghost_width):
   d_dim(max_ghost_width.getDim()),
   d_initialized(1, false)
{
   TBOX_ASSERT(max_ghost_width > IntVector(d_dim, -1));

   computeFromLevel(level,
      level_to_domain,
      level_to_level,
      max_ghost_width);

}

/*
 ************************************************************************
 * Use grid_geometry.computeBoundaryGeometry function,
 * setting up the arguments in a way that will generate
 * the coarse-fine boundary (instead of the domain boundary).
 ************************************************************************
 */
void
CoarseFineBoundary::computeFromLevel(
   const PatchLevel& level,
   const Connector& level_to_domain,
   const Connector& level_to_level,
   const IntVector& max_ghost_width)
{
//this is single block version.
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, max_ghost_width);

   clear();

   const BoxLevel& box_level = *level.getBoxLevel();
   const IntVector& ratio = level.getRatioToLevelZero();

   std::shared_ptr<BaseGridGeometry> grid_geometry(level.getGridGeometry());

   /*
    * Get the domain's periodic shift.
    */
   const IntVector periodic_shift(grid_geometry->getPeriodicShift(ratio));

   bool is_periodic = false;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      is_periodic = is_periodic || periodic_shift(i);
   }

   /*
    * Here we add some boxes outside of non-periodic boundaries to the
    * adjusted level.  For each patch that touches a regular boundary,
    * grow the patch box (and any periodic images of the patch box) by
    * the max ghost width.  Remove intersections with the periodic
    * adjusted physical domain.  Add what remains to the adjusted level.
    *
    * This will ensure that ensuing call to create boundary boxes will not
    * create boundary boxes at the locations where the level touches a
    * non-periodic physical boundary, but only where there is a coarse-fine
    * interface in the domain interior (A periodic boundary is considered
    * part of the domain interior for this purpose).
    */

   /*
    * Build a fake domain in fake_domain_list.
    *
    * The fake domain should be such that when fed to computeBoundaryBoxesOnLevel,
    * the coarse-fine boundaries are computed rather than the physical boundary.
    * computeBoundaryBoxesOnLevel defines boundaries of a patch to be the
    * parts of the grown patch box that lie outside the "domain".  So se make
    * the fake domain be everywhere there is NOT a coarse-fine boundary--or
    * everywhere there IS a physical boundary or a fine-boundary.
    */
   std::vector<BoxContainer> fake_domain(1);
   BoxContainer& fake_domain_list = fake_domain[0];

   // Every box should connect to the domain box_level.

   TBOX_ASSERT(level_to_domain.getLocalNumberOfNeighborSets() ==
      static_cast<int>(box_level.getLocalNumberOfBoxes()));

   // Add physical boundaries to the fake domain.
   IntVector physical_grow_width(max_ghost_width);
   physical_grow_width.max(IntVector::getOne(d_dim));
   for (Connector::ConstNeighborhoodIterator ei = level_to_domain.begin();
        ei != level_to_domain.end(); ++ei) {
      const Box& box = *box_level.getBoxStrict(*ei);
      BoxContainer refined_domain_nabrs;
      for (Connector::ConstNeighborIterator ni = level_to_domain.begin(ei);
           ni != level_to_domain.end(ei); ++ni) {
         refined_domain_nabrs.insert(refined_domain_nabrs.end(), *ni);
      }
      refined_domain_nabrs.refine(ratio);
      Box grow_box(box);
      grow_box.grow(physical_grow_width);
      BoxContainer physical_boundary_portion(grow_box);
      physical_boundary_portion.removeIntersections(refined_domain_nabrs);
      fake_domain_list.spliceBack(physical_boundary_portion);
   }

   // Add fine-fine boundaries to the fake domain.
#ifdef DEBUG_CHECK_ASSERTIONS
   for (BlockId::block_t b = 0; b < grid_geometry->getNumberBlocks(); ++b) {
      TBOX_ASSERT(level_to_level.getConnectorWidth()
         >= IntVector::getOne(d_dim));
   }
#endif
   BoxContainer level_neighbors;
   level_neighbors.order();
   level_to_level.getLocalNeighbors(level_neighbors);
   level_neighbors.unorder();
   fake_domain_list.spliceBack(level_neighbors);

   /*
    * Call BaseGridGeometry::computeBoundaryBoxesOnLevel with arguments contrived
    * such that they give the coarse-fine boundaries instead of the domain
    * boundaries.  The basic algorithm used by
    * BaseGridGeometry::computeBoundaryBoxesOnLevel is
    * 1. grow boxes by ghost width
    * 2. remove intersection with domain
    * 3. reorganize and classify resulting boxes
    *
    * This is how we get BaseGridGeometry::computeBoundaryGeometry to
    * compute the coarse-fine boundary instead of the physical boundary.
    *
    * Since we handle the periodic boundaries ourselves, do not treat
    * them differently from regular boundaries.  State that all boundaries
    * are non-periodic boundaries.
    *
    * Send the periodic-adjusted level boxes as the domain for the
    * remove-intersection-with-domain operation.  This causes that
    * operation to remove non-coarse-fine (that is, fine-fine) boxes
    * along the periodic boundaries, leaving the coarse-fine boundary
    * boxes.
    *
    * Send the periodic-adjusted domain for the limit-domain intersect
    * operation.  This removes the boundaries that are on the non-periodic
    * boundaries, which is what we want because there is no possibility
    * of a coarse-fine boundary there.
    */
   bool do_all_patches = true;
   const IntVector use_periodic_shift(d_dim, 0);
   grid_geometry->computeBoundaryBoxesOnLevel(
      d_boundary_boxes,
      level,
      use_periodic_shift,
      max_ghost_width,
      fake_domain,
      do_all_patches);

   d_initialized[0] = true;
}

void
CoarseFineBoundary::computeFromMultiblockLevel(
   const PatchLevel& level,
   const Connector& level_to_domain,
   const Connector& level_to_level,
   const IntVector& max_ghost_width)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, max_ghost_width);

   clear();

   const BoxLevel& box_level = *level.getBoxLevel();
   const IntVector& ratio = level.getRatioToLevelZero();
   const int level_number = level.getLevelNumber();

   /*
    * Get the number of blocks from the grid geometry.
    */
   std::shared_ptr<BaseGridGeometry> grid_geometry(level.getGridGeometry());
   size_t nblocks = grid_geometry->getNumberBlocks();

   std::vector<BoxContainer> fake_domain(nblocks);

   // Every box should connect to the domain box_level.
   TBOX_ASSERT(level_to_domain.getLocalNumberOfNeighborSets() ==
      static_cast<int>(box_level.getLocalNumberOfBoxes()));

   // Add physical boundaries to the fake domain.
   IntVector physical_grow_width(max_ghost_width);
   physical_grow_width.max(IntVector::getOne(d_dim));
   for (Connector::ConstNeighborhoodIterator ei = level_to_domain.begin();
        ei != level_to_domain.end(); ++ei) {
      const Box& box = *box_level.getBoxStrict(*ei);
      const BlockId& block_id = box.getBlockId();
      BoxContainer& fake_domain_list = fake_domain[block_id.getBlockValue()];

      BoxContainer refined_domain_nabrs;
      for (Connector::ConstNeighborIterator ni = level_to_domain.begin(ei);
           ni != level_to_domain.end(ei); ++ni) {
         if (block_id != ni->getBlockId()) {
            Box transform_box(*ni);
            grid_geometry->transformBox(transform_box,
                                        0,
                                        block_id,
                                        ni->getBlockId());
            refined_domain_nabrs.pushBack(transform_box);
         } else {
            refined_domain_nabrs.pushBack(*ni);
         }
      }
      refined_domain_nabrs.refine(ratio);
      BoxContainer physical_boundary_portion(box);
      physical_boundary_portion.grow(physical_grow_width);
      physical_boundary_portion.removeIntersections(refined_domain_nabrs);
      fake_domain_list.spliceBack(physical_boundary_portion);
   }

   for (Connector::ConstNeighborhoodIterator ei = level_to_level.begin();
        ei != level_to_level.end(); ++ei) {
      const Box& box = *box_level.getBoxStrict(*ei);
      const BlockId& block_id = box.getBlockId();
      BoxContainer& fake_domain_list = fake_domain[block_id.getBlockValue()];

      BoxContainer level_nabrs;
      for (Connector::ConstNeighborIterator ni = level_to_level.begin(ei);
           ni != level_to_level.end(ei); ++ni) {
         if (block_id != ni->getBlockId()) {
            Box transform_box(*ni);
            grid_geometry->transformBox(transform_box,
                                        level_number,
                                        block_id,
                                        ni->getBlockId());
            level_nabrs.pushBack(transform_box);
         } else {
            level_nabrs.pushBack(*ni);
         }
      }
      fake_domain_list.spliceBack(level_nabrs);
   }

   for (BlockId::block_t i = 0; i < nblocks; ++i) {
      d_initialized[i] = true;
   }

   d_boundary_boxes.clear();

   /*
    * Call BaseGridGeometry::computeBoundaryGeometry with arguments contrived
    * such that they give the coarse-fine boundaries instead of the
    * domain boundaries.  The basic algorithm used by
    * BaseGridGeometry::computeBoundaryGeometry is
    * 1. grow boxes by ghost width
    * 2. remove intersection with domain
    * 3. reorganize and classify resulting boxes
    *
    * This is how we get BaseGridGeometry::computeBoundaryGeometry to
    * compute the coarse-fine boundary instead of the physical boundary.
    *
    * Send the adjusted level boxes as the domain for the
    * remove-intersection-with-domain operation.  This causes that
    * operation to remove non-coarse-fine (that is, fine-fine) boxes
    * along the periodic boundaries, leaving the coarse-fine boundary
    * boxes.
    *
    * Send the adjusted domain for the limit-domain intersect
    * operation.  This removes the boundaries that are on the physical
    * boundaries, which is what we want because there is no possibility
    * of a coarse-fine boundary there.
    */
   bool do_all_patches = true;
   IntVector use_periodic_shift(d_dim, 0);
   grid_geometry->computeBoundaryBoxesOnLevel(
      d_boundary_boxes,
      level,
      use_periodic_shift,
      max_ghost_width,
      fake_domain,
      do_all_patches);

}

const std::vector<BoundaryBox>&
CoarseFineBoundary::getBoundaries(
   const GlobalId& global_id,
   const int boundary_type,
   const BlockId& block_id) const
{
   const BlockId::block_t& block_num = block_id.getBlockValue();
   if (!d_initialized[block_num]) {
      TBOX_ERROR("The boundary boxes have not been computed.");
   }

   BoxId box_id(global_id);
   std::map<BoxId, PatchBoundaries>::const_iterator
      mi = d_boundary_boxes.find(box_id);
   TBOX_ASSERT(mi != d_boundary_boxes.end());
   return (*mi).second[boundary_type - 1];
}

void
CoarseFineBoundary::printClassData(
   std::ostream& os) const {
   os << "\nCoarseFineBoundary::printClassData...";
   for (std::map<BoxId, PatchBoundaries>::const_iterator
        mi = d_boundary_boxes.begin(); mi != d_boundary_boxes.end(); ++mi) {
      os << "\n         patch " << (*mi).first;
      for (unsigned int btype = 0; btype < d_dim.getValue(); ++btype) {
         os << "\n                type " << btype;
         const std::vector<BoundaryBox>& array_of_boxes = (*mi).second[btype];
         int num_boxes = static_cast<int>(array_of_boxes.size());
         int bn;
         for (bn = 0; bn < num_boxes; ++bn) {
            os << "\n                           box "
               << bn << "/" << num_boxes << ":";
            os << array_of_boxes[bn].getBox();
         }
      }
   }
   os << "\n";
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
