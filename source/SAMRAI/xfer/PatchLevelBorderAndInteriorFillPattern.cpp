/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelBorderAndInteriorFillPattern::PatchLevelBorderAndInteriorFillPattern():
   d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchLevelBorderAndInteriorFillPattern::~PatchLevelBorderAndInteriorFillPattern()
{
}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */
void
PatchLevelBorderAndInteriorFillPattern::computeFillBoxesAndNeighborhoodSets(
   std::shared_ptr<hier::BoxLevel>& fill_box_level,
   std::shared_ptr<hier::Connector>& dst_to_fill,
   const hier::BoxLevel& dst_box_level,
   const hier::IntVector& fill_ghost_width,
   bool data_on_patch_border)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_box_level, fill_ghost_width);

   fill_box_level.reset(new hier::BoxLevel(
         dst_box_level.getRefinementRatio(),
         dst_box_level.getGridGeometry(),
         dst_box_level.getMPI()));

   dst_to_fill.reset(new hier::Connector(dst_box_level,
         *fill_box_level,
         fill_ghost_width));

   const hier::BoxContainer& dst_boxes = dst_box_level.getBoxes();

   hier::IntVector dst_to_dst_width(fill_ghost_width);
   if (data_on_patch_border) {
      dst_to_dst_width += hier::IntVector::getOne(fill_ghost_width.getDim());
   }

   const int dst_level_num = dst_box_level.getGridGeometry()->
      getEquivalentLevelNumber(dst_box_level.getRefinementRatio());

   const hier::Connector& dst_to_dst =
      dst_box_level.findConnector(dst_box_level,
         dst_to_dst_width,
         hier::CONNECTOR_IMPLICIT_CREATION_RULE,
         true);

   /*
    * Grow each patch box and remove the level from it, except the
    * patch box itself.  (Do not fill ghost cells that are normally
    * filled by same box_level.  Do fill ghost cells that are
    * normally filled by coarser box_level.)
    */
   hier::LocalId last_id = dst_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_boxes.realBegin());
        ni != dst_boxes.realEnd(); ++ni) {

      const hier::BoxId& gid(ni->getBoxId());
      const hier::Box& dst_box = *dst_box_level.getBox(gid);
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_box, fill_ghost_width));
      hier::Connector::ConstNeighborhoodIterator nabrs =
         dst_to_dst.find(dst_box.getBoxId());

      for (hier::Connector::ConstNeighborIterator na = dst_to_dst.begin(nabrs);
           na != dst_to_dst.end(nabrs); ++na) {
         if (!ni->isSpatiallyEqual(*na)) {
            if (dst_box.getBlockId() == na->getBlockId()) {
               fill_boxes.removeIntersections(*na);
            } else {

               std::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
                  dst_box_level.getGridGeometry());

               const hier::BlockId& dst_block_id = dst_box.getBlockId();
               const hier::BlockId& nbr_block_id = na->getBlockId();

               TBOX_ASSERT(grid_geometry->areNeighbors(dst_block_id,
                     nbr_block_id));

               hier::Transformation::RotationIdentifier rotation =
                  grid_geometry->getRotationIdentifier(dst_block_id,
                     nbr_block_id);
               hier::IntVector offset(
                  grid_geometry->getOffset(dst_block_id, nbr_block_id, dst_level_num));

               hier::Transformation transformation(rotation, offset,
                                                   nbr_block_id, dst_block_id);

               hier::Box nbr_box(*na);
               transformation.transform(nbr_box);

               fill_boxes.removeIntersections(nbr_box);

            }
         }
      }

      if (!fill_boxes.empty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());

         hier::Connector::NeighborhoodIterator base_box_itr =
            dst_to_fill->makeEmptyLocalNeighborhood(gid);
         for (hier::BoxContainer::iterator li = fill_boxes.begin();
              li != fill_boxes.end(); ++li) {
            hier::Box fill_box(*li,
                               ++last_id,
                               dst_box.getOwnerRank());
            TBOX_ASSERT(fill_box.getBlockId() == dst_box.getBlockId());
            fill_box_level->addBoxWithoutUpdate(fill_box);
            dst_to_fill->insertLocalNeighbor(fill_box, base_box_itr);
         }
      }
   }
   fill_box_level->finalize();
}

void
PatchLevelBorderAndInteriorFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_fill_boxes_on_src_proc);
   NULL_USE(dst_box_level);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   if (!needsToCommunicateDestinationFillBoxes()) {
      TBOX_ERROR(
         "PatchLevelBorderAndInteriorFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelBorderAndInteriorFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelBorderAndInteriorFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return true;
}

int
PatchLevelBorderAndInteriorFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

bool
PatchLevelBorderAndInteriorFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelBorderAndInteriorFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

}
}
