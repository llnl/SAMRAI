/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a side geometry object given the box, ghost cell width, and
 * direction information.
 *
 *************************************************************************
 */

SideGeometry::SideGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts,
   const hier::IntVector& directions):
   d_box(box),
   d_ghosts(ghosts),
   d_directions(directions)

{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);
}

/*
 *************************************************************************
 *
 * Create a side geometry object given the box and ghost cell width
 *
 *************************************************************************
 */

SideGeometry::SideGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts),
   d_directions(hier::IntVector::getOne(ghosts.getDim()))
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

SideGeometry::~SideGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two side centered box
 * geometries.  The calculateOverlap() checks whether both arguments are
 * side geometries; if so, it compuates the intersection.  If not, then
 * it calls calculateOverlap() on the source object (if retry is true)
 * to allow the source a chance to calculate the intersection.  See the
 * hier::BoxGeometry base class for more information about the protocol.
 * A pointer to null is returned if the intersection cannot be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
SideGeometry::calculateOverlap(
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

   const SideGeometry* t_dst =
      dynamic_cast<const SideGeometry *>(&dst_geometry);
   const SideGeometry* t_src =
      dynamic_cast<const SideGeometry *>(&src_geometry);

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
SideGeometry::computeDestinationBoxes(
   std::vector<hier::BoxContainer>& dst_boxes,
   const SideGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes) const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);
#endif
   TBOX_ASSERT(
      getDirectionVector() == src_geometry.getDirectionVector());

   const tbox::Dimension& dim(src_mask.getDim());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap
   hier::Box src_shift(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   transformation.transform(src_shift);
   hier::Box dst_ghost(d_box);
   dst_ghost.grow(d_ghosts);

   // Compute the intersection (if any) for each of the side directions
   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector)
      * hier::Box::grow(dst_ghost, one_vector));

   if (!quick_check.empty()) {

      const hier::IntVector& dirs = src_geometry.getDirectionVector();
      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
         if (dirs(d)) {
            const hier::Box dst_side(toSideBox(dst_ghost, d));
            const hier::Box src_side(toSideBox(src_shift, d));
            const hier::Box fill_side(toSideBox(fill_box, d));
            const hier::Box together(dst_side * src_side * fill_side);
            if (!together.empty()) {
               if (!overwrite_interior) {
                  const hier::Box int_side(toSideBox(d_box, d));
                  dst_boxes[d].removeIntersections(together, int_side);
               } else {
                  dst_boxes[d].pushBack(together);
               }
            }  // if (!together.empty())
         } // if (dirs(d))

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer side_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               side_restrict_boxes.pushBack(toSideBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(side_restrict_boxes);
         }

      }  // loop over dim && dirs(d)

   }  // if (!quick_check.empty())
}

/*
 *************************************************************************
 *
 * Convert an AMR-index space hier::Box into a side-index space box by a
 * increasing the index size by one in the axis direction.
 *
 *************************************************************************
 */

hier::Box
SideGeometry::toSideBox(
   const hier::Box& box,
   tbox::Dimension::dir_t side_normal)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT((side_normal < dim.getValue()));

   hier::Box side_box(dim);

   if (!box.empty()) {
      side_box = box;
      side_box.setUpper(side_normal, side_box.upper(side_normal) + 1);
   }

   return side_box;
}

/*
 *************************************************************************
 *
 * Compute the overlap between two side centered boxes.  The algorithm
 * is fairly straight-forward.  First, we perform a quick-and-dirty
 * intersection to see if the boxes might overlap.  If that intersection
 * is not empty, then we need to do a better job calculating the overlap
 * for each direction.  Note that the AMR index space boxes must be
 * shifted into the side centered space before we calculate the proper
 * intersections.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
SideGeometry::doOverlap(
   const SideGeometry& dst_geometry,
   const SideGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);
#endif
   TBOX_ASSERT(
      dst_geometry.getDirectionVector() == src_geometry.getDirectionVector());

   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   dst_geometry.computeDestinationBoxes(dst_boxes,
      src_geometry,
      src_mask,
      fill_box,
      overwrite_interior,
      transformation,
      dst_restrict_boxes);

   // Create the side overlap data object using the boxes and source shift

   return std::make_shared<SideOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Set up a SideOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
SideGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   const tbox::Dimension& dim(transformation.getOffset().getDim());
   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
         hier::Box side_box(SideGeometry::toSideBox(*b, d));
         dst_boxes[d].pushBack(side_box);
      }
   }

   // Create the side overlap data object using the boxes and source shift
   return std::make_shared<SideOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Transform a box
 *
 *************************************************************************
 */

void
SideGeometry::transform(
   hier::Box& box,
   int& normal_direction,
   const hier::Transformation& transformation)
{

   const tbox::Dimension& dim = box.getDim();

   if (transformation.getRotation() == hier::Transformation::NO_ROTATE &&
       transformation.getOffset() == hier::IntVector::getZero(dim)) {
      return;
   }

   if (!box.empty()) {
      const hier::Transformation::RotationIdentifier rotation =
         transformation.getRotation();

      if (rotation == hier::Transformation::NO_ROTATE) {

         transformation.transform(box);

      } else {

         box.setUpper(static_cast<hier::Box::dir_t>(normal_direction),
            box.upper(static_cast<hier::Box::dir_t>(normal_direction)) - 1);
         transformation.transform(box);
         if (dim.getValue() == 2) {
            const int rotation_num = static_cast<int>(rotation);
            if (rotation_num % 2) {
               normal_direction = (normal_direction + 1) % 2;
            }
         } else if (dim.getValue() == 3) {

            if (normal_direction == 0) {

               switch (rotation) {

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:

                     normal_direction = 0;
                     break;

                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:

                     normal_direction = 1;
                     break;

                  default:

                     normal_direction = 2;
                     break;

               }

            } else if (normal_direction == 1) {

               switch (rotation) {
                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:

                     normal_direction = 0;
                     break;

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:

                     normal_direction = 1;
                     break;

                  default:

                     normal_direction = 2;
                     break;
               }

            } else if (normal_direction == 2) {

               switch (rotation) {
                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:

                     normal_direction = 0;
                     break;

                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:

                     normal_direction = 1;
                     break;

                  default:

                     normal_direction = 2;
                     break;

               }
            }
         }

         box.setUpper(static_cast<hier::Box::dir_t>(normal_direction),
            box.upper(static_cast<hier::Box::dir_t>(normal_direction)) + 1);
      }
   }
}

/*
 *************************************************************************
 *
 * Transform a SideIndex
 *
 *************************************************************************
 */

void
SideGeometry::transform(
   SideIndex& index,
   const hier::Transformation& transformation)
{
   const tbox::Dimension& dim = index.getDim();

   if (transformation.getRotation() == hier::Transformation::NO_ROTATE &&
       transformation.getOffset() == hier::IntVector::getZero(dim)) {
      return;
   }

   const hier::Transformation::RotationIdentifier& rotation =
      transformation.getRotation();

   const int normal_direction = index.getAxis();

   for (int i = 0; i < dim.getValue(); ++i) {
      if (i != normal_direction && index(i) >= 0) {
         ++index(i);
      }
   }

   int new_normal_direction = normal_direction;
   if (dim.getValue() == 2) {
      const int rotation_num = static_cast<int>(rotation);

      TBOX_ASSERT(rotation_num <= 3);

      if (rotation_num) {

         SideIndex tmp_index(dim);
         for (int r = 0; r < rotation_num; ++r) {
            tmp_index = index;
            index(0) = tmp_index(1);
            index(1) = -tmp_index(0);
         }

         new_normal_direction = (normal_direction + rotation_num) % 2;

         index.setAxis(new_normal_direction);
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
            TBOX_ERROR("SideGeometry::transform invalid 3D RotationIdentifier.");
      }
      new_normal_direction = index.getAxis();
   }

   for (int i = 0; i < dim.getValue(); ++i) {
      if (i != new_normal_direction && index(i) > 0) {
         --index(i);
      }
   }

   index += transformation.getOffset();
}

void
SideGeometry::rotateAboutAxis(SideIndex& index,
                              const tbox::Dimension::dir_t axis,
                              const int num_rotations)
{
   const tbox::Dimension& dim = index.getDim();
   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   SideIndex tmp_index(dim);
   for (int j = 0; j < num_rotations; ++j) {
      tmp_index = index;
      index(a) = tmp_index(b);
      index(b) = -tmp_index(a);
   }

   int new_normal_direction = index.getAxis();
   if (new_normal_direction != axis) {
      for (int j = 0; j < num_rotations; ++j) {
         new_normal_direction = new_normal_direction == a ? b : a;
      }
   }
   index.setAxis(new_normal_direction);
}

SideIterator
SideGeometry::begin(
   const hier::Box& box,
   tbox::Dimension::dir_t axis)
{
   return SideIterator(box, axis, true);
}

SideIterator
SideGeometry::end(
   const hier::Box& box,
   tbox::Dimension::dir_t axis)
{
   return SideIterator(box, axis, false);
}

}
}
