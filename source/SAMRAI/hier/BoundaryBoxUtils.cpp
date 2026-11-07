/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Generic utilities for boundary box calculus.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {

BoundaryBoxUtils::BoundaryBoxUtils(
   const BoundaryBox& bbox):
   d_bbox(bbox),
   d_outward(bbox.getDim(), 0)
{
   computeOutwardShift();
}

void
BoundaryBoxUtils::computeOutwardShift()
{

   const tbox::Dimension& dim(d_bbox.getDim());
   /*
    * Note that d_outward contains information that is redundant
    * with respect to the boundary box.  The values of d_outward
    * depends strictly the location of the boundary box.
    */
   const int lidx = d_bbox.getLocationIndex();

   switch (d_bbox.getBoundaryType()) {

      // Note: number of non-zero in d_outward is the same as boundary type.

      case 1:
      {
         int i = lidx / 2;
         d_outward(i) = lidx % 2 == 0 ? -1 /* lower side */ : 1 /* upper side */;
      }
      break;

      case 2:
      {
         if (dim.getValue() == 2) {
            d_outward(0) = lidx % 2 == 0 ? -1 : 1;
            d_outward(1) = lidx / 2 == 0 ? -1 : 1;
         } else if (dim.getValue() == 3) {
            const int dir = 2 - (lidx / 4);
            const int rem = lidx % 4;
            if (dir == 0) {
               // Nonzero in dirs 1 and 2.
               d_outward(1) = rem % 2 == 0 ? -1 : 1;
               d_outward(2) = rem / 2 == 0 ? -1 : 1;
            } else if (dir == 1) {
               // Nonzero in dirs 0 and 2.
               d_outward(0) = rem % 2 == 0 ? -1 : 1;
               d_outward(2) = rem / 2 == 0 ? -1 : 1;
            } else {
               // Nonzero in dirs 0 and 1.
               d_outward(0) = rem % 2 == 0 ? -1 : 1;
               d_outward(1) = rem / 2 == 0 ? -1 : 1;
            }
         } else {
            TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
               << "boundary direction for " << d_bbox.getBox());
         }
      }
      break;

      case 3:
      {
         if (dim.getValue() == 3) {
            d_outward(0) = lidx % 2 == 0 ? -1 : 1;
            d_outward(1) = (lidx % 4) / 2 == 0 ? -1 : 1;
            d_outward(2) = lidx / 4 == 0 ? -1 : 1;
         } else {
            TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
               << "boundary direction for " << d_bbox.getBox());
         }
      }
      break;

      default:
         TBOX_ERROR("BoundaryBoxUtils cannot compute\n"
         << "boundary direction for type "
         << d_bbox.getBoundaryType() << " in " << dim << "D");
         break;
   }
}

void
BoundaryBoxUtils::stretchBoxToGhostWidth(
   Box& box,
   const IntVector& ghost_cell_width) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(d_bbox, box);

   const tbox::Dimension& dim(d_bbox.getDim());

   TBOX_ASSERT(ghost_cell_width >= IntVector::getZero(dim));

   box = d_bbox.getBox();
   for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
      /*
       * If gcw along direction d is > 1, stretch it out to that width.
       * If gcw a long direction d is 0, shrink the box down to nothing
       * in that direction.
       */
      if (d_outward(d) == -1) {
         if (ghost_cell_width(d) > 1) box.growLower(d, ghost_cell_width(d) - 1);
         else box.setLower(d, box.upper(d) - (ghost_cell_width(d) - 1));
      } else if (d_outward(d) == 1) {
         if (ghost_cell_width(d) > 1) box.growUpper(d, ghost_cell_width(d) - 1);
         else box.setUpper(d, box.lower(d) + (ghost_cell_width(d) - 1));
      }
   }
}

void
BoundaryBoxUtils::extendBoxOutward(
   Box& box,
   const IntVector& extension) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(d_bbox, box);

   const tbox::Dimension& dim(d_bbox.getDim());

   box = d_bbox.getBox();
   for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
      if (d_outward(d) == -1) box.growLower(d, extension(d));
      else if (d_outward(d) == 1) box.growUpper(d, extension(d));
   }
}

/*
 ************************************************************************
 * Make surface box on boundary using standard boundary box
 ************************************************************************
 */

Box
BoundaryBoxUtils::getSurfaceBoxFromBoundaryBox() const
{
   if (d_bbox.getBoundaryType() != 1) {
      TBOX_ERROR("BoundaryBoxUtils::getSurfaceBoxFromBoundaryBox\n"
         << "called with improper boundary box\n");
   }
   Box side_index_box = d_bbox.getBox();
   int location_index = d_bbox.getLocationIndex();
   if (location_index % 2 == 0) {
      /*
       * On the min index side, the face indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      side_index_box.shift(static_cast<tbox::Dimension::dir_t>(location_index / 2), 1);
   }
   return side_index_box;
}

/*
 ************************************************************************
 * Trim a boundary box so it does not stick out past the corners of a
 * patch.
 ************************************************************************
 */

BoundaryBox
BoundaryBoxUtils::trimBoundaryBox(
   const Box& limit_box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(d_bbox, limit_box);

   const tbox::Dimension& dim(d_bbox.getDim());

   TBOX_ASSERT(d_bbox.getBoundaryType() < dim.getValue());

   const Box& bbox = d_bbox.getBox();
   const Index& plo = limit_box.lower();
   const Index& pup = limit_box.upper();
   const Index& blo = bbox.lower();
   const Index& bup = bbox.upper();
   Index newlo(dim), newup(dim);

   if (d_bbox.getBoundaryType() == 1) {
      /*
       * Loop through directions.
       * Preserve box size in direction normal to boundary.
       * Trim box size in direction parallel to boundary.
       */
      const int boundary_normal = d_bbox.getLocationIndex() / 2;
      int d;
      for (d = 0; d < dim.getValue(); ++d) {
         if (d == boundary_normal) {
            newlo(d) = blo(d);
            newup(d) = bup(d);
         } else {
            // On min side, use max between boundary and patch boxes.
            newlo(d) = tbox::MathUtilities<int>::Max(blo(d), plo(d));
            // On max side, use min between boundary and patch boxes.
            newup(d) = tbox::MathUtilities<int>::Min(bup(d), pup(d));
         }
      }
   } else if (d_bbox.getBoundaryType() == 2) {
      /*
       * Loop through directions.
       * Preserve box size in direction normal to boundary.
       * Trim box size in direction parallel to boundary.
       */
      const int boundary_dir = 4 - (d_bbox.getLocationIndex() / 4);
      int d;
      for (d = 0; d < dim.getValue(); ++d) {
         if (d == boundary_dir) {
            // On min side, use max between boundary and patch boxes.
            newlo(d) = tbox::MathUtilities<int>::Max(blo(d), plo(d));
            // On max side, use min between boundary and patch boxes.
            newup(d) = tbox::MathUtilities<int>::Min(bup(d), pup(d));
         } else {
            newlo(d) = blo(d);
            newup(d) = bup(d);
         }
      }
   }

   const Box newbox(newlo, newup, d_bbox.getBox().getBlockId());
   const BoundaryBox newbbox(newbox,
                             d_bbox.getBoundaryType(),
                             d_bbox.getLocationIndex());

   return newbbox;
}

}
}
