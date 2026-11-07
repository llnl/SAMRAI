/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.h"
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

PatchLevelEnhancedFillPattern::PatchLevelEnhancedFillPattern():
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

PatchLevelEnhancedFillPattern::~PatchLevelEnhancedFillPattern()
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
PatchLevelEnhancedFillPattern::computeFillBoxesAndNeighborhoodSets(
   std::shared_ptr<hier::BoxLevel>& fill_box_level,
   std::shared_ptr<hier::Connector>& dst_to_fill,
   const hier::BoxLevel& dst_box_level,
   const hier::IntVector& fill_ghost_width,
   bool data_on_patch_border)
{
   NULL_USE(data_on_patch_border);
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_box_level, fill_ghost_width);

   fill_box_level.reset(new hier::BoxLevel(
         dst_box_level.getRefinementRatio(),
         dst_box_level.getGridGeometry(),
         dst_box_level.getMPI()));

   dst_to_fill.reset(new hier::Connector(dst_box_level,
         *fill_box_level,
         fill_ghost_width));

   std::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
      dst_box_level.getGridGeometry());

   const hier::BoxContainer& dst_boxes = dst_box_level.getBoxes();

   hier::LocalId last_id = dst_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_boxes.realBegin());
        ni != dst_boxes.realEnd(); ++ni) {
      const hier::Box& dst_box = *ni;
      const hier::BoxId& dst_box_id = dst_box.getBoxId();
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_box, fill_ghost_width));

      hier::BoxContainer constructed_fill_boxes;

      hier::Connector::NeighborhoodIterator base_box_itr =
         dst_to_fill->findLocal(dst_box_id);
      bool has_base_box = base_box_itr != dst_to_fill->end();

      for (hier::BaseGridGeometry::ConstNeighborIterator ni =
              grid_geometry->begin(dst_box.getBlockId());
           ni != grid_geometry->end(dst_box.getBlockId()); ++ni) {
         const hier::BaseGridGeometry::Neighbor& nbr = *ni;
         if (nbr.isSingularity()) {

            hier::BoxContainer encon_boxes(nbr.getTransformedDomain());
            encon_boxes.refine(dst_box_level.getRefinementRatio());
            encon_boxes.intersectBoxes(fill_boxes);
            encon_boxes.removeIntersections(constructed_fill_boxes);

            if (encon_boxes.size()) {

               if (!has_base_box) {
                  base_box_itr = dst_to_fill->makeEmptyLocalNeighborhood(
                        dst_box_id);
                  has_base_box = true;
               }
               for (hier::BoxContainer::iterator ei = encon_boxes.begin();
                    ei != encon_boxes.end(); ++ei) {

                  hier::Box fill_box(
                     *ei,
                     ++last_id,
                     dst_box.getOwnerRank());

                  TBOX_ASSERT(fill_box.getBlockId() == dst_box.getBlockId());

                  fill_box_level->addBoxWithoutUpdate(fill_box);

                  dst_to_fill->insertLocalNeighbor(
                     fill_box,
                     base_box_itr);

                  constructed_fill_boxes.pushBack(*ei);
               }
            }
         }
      }

      d_max_fill_boxes = tbox::MathUtilities<int>::Max(
            d_max_fill_boxes,
            constructed_fill_boxes.size());
   }
   fill_box_level->finalize();
}

void
PatchLevelEnhancedFillPattern::computeDestinationFillBoxesOnSourceProc(
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
         "PatchLevelEnhancedFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelEnhancedFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelEnhancedFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return false;
}

bool
PatchLevelEnhancedFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelEnhancedFillPattern::fillingEnhancedConnectivityOnly() const
{
   return true;
}

int
PatchLevelEnhancedFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
