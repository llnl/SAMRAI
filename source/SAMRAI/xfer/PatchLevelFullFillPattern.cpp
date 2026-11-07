/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelFullFillPattern::PatchLevelFullFillPattern():
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

PatchLevelFullFillPattern::~PatchLevelFullFillPattern()
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
PatchLevelFullFillPattern::computeFillBoxesAndNeighborhoodSets(
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

   const hier::BoxContainer& dst_boxes = dst_box_level.getBoxes();

   for (hier::RealBoxConstIterator ni(dst_boxes.realBegin());
        ni != dst_boxes.realEnd(); ++ni) {
      const hier::Box& dst_box = *ni;
      hier::Box fill_box(dst_box);
      fill_box.grow(fill_ghost_width);
      fill_box_level->addBoxWithoutUpdate(fill_box);
      dst_to_fill->insertLocalNeighbor(fill_box, dst_box.getBoxId());
      TBOX_ASSERT(dst_to_fill->numLocalNeighbors(dst_box.getBoxId()) == 1);
   }
   fill_box_level->finalize();
}

/*
 *************************************************************************
 *
 * computeDestinationFillBoxesOnSourceProc
 *
 *************************************************************************
 */

void
PatchLevelFullFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_box_level, fill_ghost_width);

   const tbox::Dimension& dim(fill_ghost_width.getDim());
   const hier::IntVector& ratio(dst_box_level.getRefinementRatio());

   bool is_periodic = false;
   if (dst_box_level.getGridGeometry()->getPeriodicShift(ratio) !=
       hier::IntVector::getZero(dim)) {
      is_periodic = true;
   }

   /*
    * src_to_dst initialized only when there is a src box_level.
    * Without the src box_level, we do not need to compute
    * dst_fill_boxes_on_src_proc.
    *
    * For PatchLevelFullFillPattern, the src owner can compute fill boxes
    * for all its dst neighbors using local data.  This info is
    * stored in dst_fill_boxes_on_src_proc.
    */
   bool ordered = true;
   hier::BoxContainer all_dst_nabrs(ordered);
   if (is_periodic) {
      hier::BoxContainer tmp_nabrs(ordered);
      src_to_dst.getLocalNeighbors(tmp_nabrs);
      tmp_nabrs.unshiftPeriodicImageBoxes(
         all_dst_nabrs,
         dst_box_level.getRefinementRatio(),
         dst_box_level.getGridGeometry()->getPeriodicShiftCatalog());
   } else {
      src_to_dst.getLocalNeighbors(all_dst_nabrs);
   }
   for (hier::BoxContainer::const_iterator na = all_dst_nabrs.begin();
        na != all_dst_nabrs.end(); ++na) {
      FillSet::Iterator dst_fill_boxes_iter =
         dst_fill_boxes_on_src_proc.insert(na->getBoxId()).first;
      hier::Box fill_box(*na);
      fill_box.grow(fill_ghost_width);
      dst_fill_boxes_on_src_proc.insert(dst_fill_boxes_iter, fill_box);
      d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
            static_cast<int>(dst_fill_boxes_on_src_proc.numNeighbors(
                                dst_fill_boxes_iter)));
   }
}

bool
PatchLevelFullFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return false;
}

bool
PatchLevelFullFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return true;
}

bool
PatchLevelFullFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelFullFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

int
PatchLevelFullFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
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
