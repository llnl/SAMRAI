/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/Box.h"
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

PatchLevelBorderFillPattern::PatchLevelBorderFillPattern():
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

PatchLevelBorderFillPattern::~PatchLevelBorderFillPattern()
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
PatchLevelBorderFillPattern::computeFillBoxesAndNeighborhoodSets(
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

   const int dst_level_num = dst_box_level.getGridGeometry()->
      getEquivalentLevelNumber(dst_box_level.getRefinementRatio());

   hier::IntVector dst_to_dst_width(fill_ghost_width);
   if (data_on_patch_border) {
      dst_to_dst_width += hier::IntVector::getOne(fill_ghost_width.getDim());
   }

   const hier::Connector& dst_to_dst =
      dst_box_level.findConnector(dst_box_level,
         dst_to_dst_width,
         hier::CONNECTOR_IMPLICIT_CREATION_RULE,
         true);

   /*
    * To get the level border, grow each patch box and remove
    * the level from it.
    */
   hier::LocalId last_id = dst_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_boxes.realBegin());
        ni != dst_boxes.realEnd(); ++ni) {
      const hier::Box& dst_box = *ni;
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_box, fill_ghost_width));
      hier::Connector::ConstNeighborhoodIterator nabrs =
         dst_to_dst.find(dst_box.getBoxId());
      for (hier::Connector::ConstNeighborIterator na = dst_to_dst.begin(nabrs);
           na != dst_to_dst.end(nabrs); ++na) {
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

      if (!fill_boxes.empty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());
         hier::Connector::NeighborhoodIterator base_box_itr =
            dst_to_fill->makeEmptyLocalNeighborhood(dst_box.getBoxId());
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
PatchLevelBorderFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_box_level);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   NULL_USE(dst_fill_boxes_on_src_proc);
   if (!needsToCommunicateDestinationFillBoxes()) {
      TBOX_ERROR(
         "PatchLevelBorderFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelBorderFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelBorderFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return false;
}

bool
PatchLevelBorderFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelBorderFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

int
PatchLevelBorderFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
