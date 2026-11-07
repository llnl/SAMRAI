/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a cell geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

CellGeometry::CellGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

CellGeometry::~CellGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two cell centered box
 * geometries.  The calculateOverlap() checks whether both arguments are
 * cell geometries; if so, it computes the intersection.  If not, then
 * it calls calculateOverlap() on the source object (if retry is true)
 * to allow the source a chance to calculate the intersection.  See the
 * hier::BoxGeometry base class for more information about the protocol.
 * A pointer to null is returned if the intersection cannot be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
CellGeometry::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const bool retry,
   const hier::BoxContainer& dst_restrict_boxes) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(d_box, src_mask);

   const CellGeometry* t_dst =
      dynamic_cast<const CellGeometry *>(&dst_geometry);
   const CellGeometry* t_src =
      dynamic_cast<const CellGeometry *>(&src_geometry);

   std::shared_ptr<hier::BoxOverlap> over;
   if ((t_src != 0) && (t_dst != 0)) {
      over = doOverlap(*t_dst, *t_src, src_mask, fill_box, overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if (retry) {
      over = src_geometry.calculateOverlap(dst_geometry, src_geometry,
            src_mask, fill_box, overwrite_interior,
            transformation, false,
            dst_restrict_boxes);
   }
   return over;
}

/*
 *************************************************************************
 *
 * Compute the boxes that will be used to construct an overlap object
 *
 *************************************************************************
 */

void
CellGeometry::computeDestinationBoxes(
   hier::BoxContainer& dst_boxes,
   const CellGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(d_box, src_mask);

   // Translate the source box and grow the destination box by the ghost cells

   hier::Box src_box(src_geometry.d_box);
   src_box.grow(src_geometry.d_ghosts);
   src_box = src_box * src_mask;
   transformation.transform(src_box);
   hier::Box dst_ghost(d_box);
   dst_ghost.grow(d_ghosts);

   // Convert the boxes into cell space and compute the intersection

   const hier::Box together(dst_ghost * src_box * fill_box);

   if (!together.empty()) {
      if (!overwrite_interior) {
         dst_boxes.removeIntersections(together, d_box);
      } else {
         dst_boxes.pushBack(together);
      }
   }

   if (!dst_boxes.empty() && !dst_restrict_boxes.empty()) {
      dst_boxes.intersectBoxes(dst_restrict_boxes);
   }
}

/*
 *************************************************************************
 *
 * Set up a CellOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
CellGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   // Create the cell overlap data object using the boxes and source shift
   return std::make_shared<CellOverlap>(boxes, transformation);

}

/*
 *************************************************************************
 *
 * Transform a CellIndex
 *
 *************************************************************************
 */

void
CellGeometry::transform(
   CellIndex& index,
   const hier::Transformation& transformation)
{
   const tbox::Dimension& dim = index.getDim();

   if (transformation.getRotation() == hier::Transformation::NO_ROTATE &&
       transformation.getOffset() == hier::IntVector::getZero(dim)) {
      return;
   }

   const hier::Transformation::RotationIdentifier& rotation =
      transformation.getRotation();
   if (dim.getValue() == 1) {
      const int rotation_num = static_cast<int>(rotation);
      if (rotation_num > 1) {
         TBOX_ERROR("CellGeometry::transform invalid 1D RotationIdentifier.");
      }

      if (rotation_num) {
         CellIndex tmp_index(index);
         index(0) = -tmp_index(0) - 1;
      }
   } else if (dim.getValue() == 2) {
      const int rotation_num = static_cast<int>(rotation);
      if (rotation_num > 3) {
         TBOX_ERROR("CellGeometry::transform invalid 2D RotationIdentifier.");
      }

      if (rotation_num) {
         CellIndex tmp_index(dim);
         for (int r = 0; r < rotation_num; ++r) {
            tmp_index = index;
            index(0) = tmp_index(1);
            index(1) = -tmp_index(0) - 1;
         }
      }

   } else if (dim.getValue() == 3) {
      switch (rotation) {

         case hier::Transformation::NO_ROTATE:
            break;

         case hier::Transformation::KUP_IUP_JUP:
            rotateAboutAxis(index, 0, 3);
            rotateAboutAxis(index, 2, 3);
            break;

         case hier::Transformation::JUP_KUP_IUP:
            rotateAboutAxis(index, 1, 1);
            rotateAboutAxis(index, 2, 1);
            break;

         case hier::Transformation::IDOWN_KUP_JUP:
            rotateAboutAxis(index, 1, 2);
            rotateAboutAxis(index, 0, 3);
            break;

         case hier::Transformation::KUP_JUP_IDOWN:
            rotateAboutAxis(index, 1, 3);
            break;

         case hier::Transformation::JUP_IDOWN_KUP:
            rotateAboutAxis(index, 2, 1);
            break;

         case hier::Transformation::KDOWN_JUP_IUP:
            rotateAboutAxis(index, 1, 1);
            break;

         case hier::Transformation::IUP_KDOWN_JUP:
            rotateAboutAxis(index, 0, 3);
            break;

         case hier::Transformation::JUP_IUP_KDOWN:
            rotateAboutAxis(index, 0, 2);
            rotateAboutAxis(index, 2, 3);
            break;

         case hier::Transformation::KDOWN_IDOWN_JUP:
            rotateAboutAxis(index, 0, 3);
            rotateAboutAxis(index, 2, 1);
            break;

         case hier::Transformation::IDOWN_JUP_KDOWN:
            rotateAboutAxis(index, 1, 2);
            break;

         case hier::Transformation::JUP_KDOWN_IDOWN:
            rotateAboutAxis(index, 0, 3);
            rotateAboutAxis(index, 1, 3);
            break;

         case hier::Transformation::JDOWN_IUP_KUP:
            rotateAboutAxis(index, 2, 3);
            break;

         case hier::Transformation::IUP_KUP_JDOWN:
            rotateAboutAxis(index, 0, 1);
            break;

         case hier::Transformation::KUP_JDOWN_IUP:
            rotateAboutAxis(index, 0, 2);
            rotateAboutAxis(index, 1, 1);
            break;

         case hier::Transformation::JDOWN_KUP_IDOWN:
            rotateAboutAxis(index, 0, 1);
            rotateAboutAxis(index, 1, 3);
            break;

         case hier::Transformation::IDOWN_JDOWN_KUP:
            rotateAboutAxis(index, 0, 2);
            rotateAboutAxis(index, 1, 2);
            break;

         case hier::Transformation::KUP_IDOWN_JDOWN:
            rotateAboutAxis(index, 0, 1);
            rotateAboutAxis(index, 2, 1);
            break;

         case hier::Transformation::JDOWN_KDOWN_IUP:
            rotateAboutAxis(index, 0, 3);
            rotateAboutAxis(index, 1, 1);
            break;

         case hier::Transformation::KDOWN_IUP_JDOWN:
            rotateAboutAxis(index, 0, 1);
            rotateAboutAxis(index, 2, 3);
            break;

         case hier::Transformation::IUP_JDOWN_KDOWN:
            rotateAboutAxis(index, 0, 2);
            break;

         case hier::Transformation::JDOWN_IDOWN_KDOWN:
            rotateAboutAxis(index, 0, 2);
            rotateAboutAxis(index, 2, 1);
            break;

         case hier::Transformation::KDOWN_JDOWN_IDOWN:
            rotateAboutAxis(index, 0, 2);
            rotateAboutAxis(index, 1, 3);
            break;

         case hier::Transformation::IDOWN_KDOWN_JDOWN:
            rotateAboutAxis(index, 1, 2);
            rotateAboutAxis(index, 0, 1);
            break;

         default:
            TBOX_ERROR("CellGeometry::transform invalid 3D RotationIdentifier.");
      }

   } else {
      TBOX_ERROR("CellGeometry::transform implemented for 2D and 3d only.");
   }

   index += transformation.getOffset();
}

void
CellGeometry::rotateAboutAxis(CellIndex& index,
                              const int axis,
                              const int num_rotations)
{
   const tbox::Dimension& dim = index.getDim();
   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   CellIndex tmp_index(dim);
   for (int j = 0; j < num_rotations; ++j) {
      tmp_index = index;
      index(a) = tmp_index(b);
      index(b) = -tmp_index(a) - 1;
   }
}

CellIterator
CellGeometry::begin(
   const hier::Box& box)
{
   return CellIterator(box, true);
}

CellIterator
CellGeometry::end(
   const hier::Box& box)
{
   return CellIterator(box, false);
}

}
}
