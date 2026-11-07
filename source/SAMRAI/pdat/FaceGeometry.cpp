/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a face geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

FaceGeometry::FaceGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

FaceGeometry::~FaceGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two face centered box
 * geometries.  The calculateOverlap() checks whether both arguments are
 * face geometries; if so, it compuates the intersection.  If not, then
 * it calls calculateOverlap() on the source object (if retry is true)
 * to allow the source a chance to calculate the intersection.  See the
 * hier::BoxGeometry base class for more information about the protocol.
 * A pointer to null is returned if the intersection cannot be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
FaceGeometry::calculateOverlap(
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

   const FaceGeometry* t_dst =
      dynamic_cast<const FaceGeometry *>(&dst_geometry);
   const FaceGeometry* t_src =
      dynamic_cast<const FaceGeometry *>(&src_geometry);

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
 * Convert an AMR-index space hier::Box into a face-index space box by a
 * cyclic shift of indices.
 *
 *************************************************************************
 */

hier::Box
FaceGeometry::toFaceBox(
   const hier::Box& box,
   tbox::Dimension::dir_t face_normal)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT((face_normal < dim.getValue()));

   hier::Box face_box(dim);

   if (!box.empty()) {
      const tbox::Dimension::dir_t x = face_normal;
      face_box.setLower(0, box.lower(x));
      face_box.setUpper(0, box.upper(x) + 1);
      for (tbox::Dimension::dir_t i = 1; i < dim.getValue(); ++i) {
         const tbox::Dimension::dir_t y =
            static_cast<tbox::Dimension::dir_t>((face_normal + i) % dim.getValue());
         face_box.setLower(i, box.lower(y));
         face_box.setUpper(i, box.upper(y));
      }
      face_box.setBlockId(box.getBlockId());
   }

   return face_box;
}

/*
 *************************************************************************
 *
 * Compute the boxes that will be used to construct an overlap object
 *
 *************************************************************************
 */

void
FaceGeometry::computeDestinationBoxes(
   std::vector<hier::BoxContainer>& dst_boxes,
   const FaceGeometry& src_geometry,
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

   const tbox::Dimension& dim(src_mask.getDim());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(d_box, d_ghosts));

   // Compute the intersection (if any) for each of the face directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {
      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
         const hier::Box dst_face(toFaceBox(dst_ghost, d));
         const hier::Box src_face(toFaceBox(src_shift, d));
         const hier::Box fill_face(toFaceBox(fill_box, d));
         const hier::Box together(dst_face * src_face * fill_face);
         if (!together.empty()) {
            if (!overwrite_interior) {
               const hier::Box int_face(toFaceBox(d_box, d));
               dst_boxes[d].removeIntersections(together, int_face);
            } else {
               dst_boxes[d].pushBack(together);
            }
         }  // if (!together.empty())

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer face_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               face_restrict_boxes.pushBack(toFaceBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(face_restrict_boxes);
         }
      }  // loop over dim
   }  // !quick_check.empty()
}

/*
 *************************************************************************
 *
 * Compute the overlap between two face centered boxes.  The algorithm
 * is fairly straight-forward.  First, we perform a quick-and-dirty
 * intersection to see if the boxes might overlap.  If that intersection
 * is not empty, then we need to do a better job calculating the overlap
 * for each direction.  Note that the AMR index space boxes must be
 * shifted into the face centered space before we calculate the proper
 * intersections.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
FaceGeometry::doOverlap(
   const FaceGeometry& dst_geometry,
   const FaceGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   dst_geometry.computeDestinationBoxes(dst_boxes,
      src_geometry,
      src_mask,
      fill_box,
      overwrite_interior,
      transformation,
      dst_restrict_boxes);

   // Create the face overlap data object using the boxes and source shift

   return std::make_shared<FaceOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Set up a FaceOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
FaceGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   const tbox::Dimension& dim(transformation.getOffset().getDim());
   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {
         hier::Box face_box(FaceGeometry::toFaceBox(*b, d));
         dst_boxes[d].pushBack(face_box);
      }
   }

   // Create the face overlap data object using the boxes and source shift
   return std::make_shared<FaceOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Transform a box
 *
 *************************************************************************
 */

void
FaceGeometry::transform(
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

      hier::Box cell_box(dim);
      for (int d = 0; d < dim.getValue(); ++d) {
         int cell_dim = (normal_direction + d) % dim.getValue();
         cell_box.setLower(static_cast<hier::Box::dir_t>(cell_dim),
            box.lower(static_cast<hier::Box::dir_t>(d)));
         cell_box.setUpper(static_cast<hier::Box::dir_t>(cell_dim),
            box.upper(static_cast<hier::Box::dir_t>(d)));
      }
      cell_box.setUpper(static_cast<hier::Box::dir_t>(normal_direction),
         cell_box.upper(static_cast<hier::Box::dir_t>(normal_direction)) - 1);
      cell_box.setBlockId(box.getBlockId());
      transformation.transform(cell_box);

      if (rotation != hier::Transformation::NO_ROTATE) {

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
      }

      for (int d = 0; d < dim.getValue(); ++d) {
         int cell_dim = (normal_direction + d) % dim.getValue();
         box.setLower(static_cast<hier::Box::dir_t>(d),
            cell_box.lower(static_cast<hier::Box::dir_t>(cell_dim)));
         box.setUpper(static_cast<hier::Box::dir_t>(d),
            cell_box.upper(static_cast<hier::Box::dir_t>(cell_dim)));
      }

      box.setUpper(0, box.upper(0) + 1);
      box.setBlockId(cell_box.getBlockId());
   }
}

/*
 *************************************************************************
 *
 * Transform a FaceIndex
 *
 *************************************************************************
 */

void
FaceGeometry::transform(
   FaceIndex& index,
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

   FaceIndex rotate_index(dim);
   for (int d = 0; d < dim.getValue(); ++d) {
      int rotate_dim = (normal_direction + d) % dim.getValue();
      rotate_index(rotate_dim) = index(d);
      if (d != 0 && rotate_index(rotate_dim) >= 0) {
         ++rotate_index(rotate_dim);
      }
   }
   rotate_index.setAxis(normal_direction);

   int new_normal_direction = normal_direction;
   if (dim.getValue() == 2) {
      const int rotation_num = static_cast<int>(rotation);

      TBOX_ASSERT(rotation_num <= 3);

      if (rotation_num) {

         hier::Index tmp_index(dim);
         for (int r = 0; r < rotation_num; ++r) {
            tmp_index = rotate_index;
            rotate_index(0) = tmp_index(1);
            rotate_index(1) = -tmp_index(0);
         }

         new_normal_direction = (normal_direction + rotation_num) % 2;

         index.setAxis(new_normal_direction);
      }
   } else if (dim.getValue() == 3) {

      switch (rotation) {

         case hier::Transformation::NO_ROTATE:
            break;

         case hier::Transformation::KUP_IUP_JUP:
            rotateAboutAxis(rotate_index, 0, 3);
            rotateAboutAxis(rotate_index, 2, 3);
            break;

         case hier::Transformation::JUP_KUP_IUP:
            rotateAboutAxis(rotate_index, 1, 1);
            rotateAboutAxis(rotate_index, 2, 1);
            break;

         case hier::Transformation::IDOWN_KUP_JUP:
            rotateAboutAxis(rotate_index, 1, 2);
            rotateAboutAxis(rotate_index, 0, 3);
            break;

         case hier::Transformation::KUP_JUP_IDOWN:
            rotateAboutAxis(rotate_index, 1, 3);
            break;

         case hier::Transformation::JUP_IDOWN_KUP:
            rotateAboutAxis(rotate_index, 2, 1);
            break;

         case hier::Transformation::KDOWN_JUP_IUP:
            rotateAboutAxis(rotate_index, 1, 1);
            break;

         case hier::Transformation::IUP_KDOWN_JUP:
            rotateAboutAxis(rotate_index, 0, 3);
            break;

         case hier::Transformation::JUP_IUP_KDOWN:
            rotateAboutAxis(rotate_index, 0, 2);
            rotateAboutAxis(rotate_index, 2, 3);
            break;

         case hier::Transformation::KDOWN_IDOWN_JUP:
            rotateAboutAxis(rotate_index, 0, 3);
            rotateAboutAxis(rotate_index, 2, 1);
            break;

         case hier::Transformation::IDOWN_JUP_KDOWN:
            rotateAboutAxis(rotate_index, 1, 2);
            break;

         case hier::Transformation::JUP_KDOWN_IDOWN:
            rotateAboutAxis(rotate_index, 0, 3);
            rotateAboutAxis(rotate_index, 1, 3);
            break;

         case hier::Transformation::JDOWN_IUP_KUP:
            rotateAboutAxis(rotate_index, 2, 3);
            break;

         case hier::Transformation::IUP_KUP_JDOWN:
            rotateAboutAxis(rotate_index, 0, 1);
            break;

         case hier::Transformation::KUP_JDOWN_IUP:
            rotateAboutAxis(rotate_index, 0, 2);
            rotateAboutAxis(rotate_index, 1, 1);
            break;

         case hier::Transformation::JDOWN_KUP_IDOWN:
            rotateAboutAxis(rotate_index, 0, 1);
            rotateAboutAxis(rotate_index, 1, 3);
            break;

         case hier::Transformation::IDOWN_JDOWN_KUP:
            rotateAboutAxis(rotate_index, 0, 2);
            rotateAboutAxis(rotate_index, 1, 2);
            break;

         case hier::Transformation::KUP_IDOWN_JDOWN:
            rotateAboutAxis(rotate_index, 0, 1);
            rotateAboutAxis(rotate_index, 2, 1);
            break;

         case hier::Transformation::JDOWN_KDOWN_IUP:
            rotateAboutAxis(rotate_index, 0, 3);
            rotateAboutAxis(rotate_index, 1, 1);
            break;

         case hier::Transformation::KDOWN_IUP_JDOWN:
            rotateAboutAxis(rotate_index, 0, 1);
            rotateAboutAxis(rotate_index, 2, 3);
            break;

         case hier::Transformation::IUP_JDOWN_KDOWN:
            rotateAboutAxis(rotate_index, 0, 2);
            break;

         case hier::Transformation::JDOWN_IDOWN_KDOWN:
            rotateAboutAxis(rotate_index, 0, 2);
            rotateAboutAxis(rotate_index, 2, 1);
            break;

         case hier::Transformation::KDOWN_JDOWN_IDOWN:
            rotateAboutAxis(rotate_index, 0, 2);
            rotateAboutAxis(rotate_index, 1, 3);
            break;

         case hier::Transformation::IDOWN_KDOWN_JDOWN:
            rotateAboutAxis(rotate_index, 1, 2);
            rotateAboutAxis(rotate_index, 0, 1);
            break;

         default:
            TBOX_ERROR("FaceGeometry::transform invalid 3D RotationIdentifier.");
      }
      new_normal_direction = rotate_index.getAxis();
   }

   for (int d = 0; d < dim.getValue(); ++d) {
      if (d != new_normal_direction && rotate_index(d) > 0) {
         --rotate_index(d);
      }
   }

   rotate_index += transformation.getOffset();
   for (int d = 0; d < dim.getValue(); ++d) {
      int rotate_dim = (new_normal_direction + d) % dim.getValue();
      index(d) = rotate_index(rotate_dim);
   }

   index.setAxis(new_normal_direction);

}

void
FaceGeometry::rotateAboutAxis(FaceIndex& index,
                              const int axis,
                              const int num_rotations)
{
   const tbox::Dimension& dim = index.getDim();
   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   FaceIndex tmp_index(dim);
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

FaceIterator
FaceGeometry::begin(
   const hier::Box& box,
   tbox::Dimension::dir_t axis)
{
   return FaceIterator(box, axis, true);
}

FaceIterator
FaceGeometry::end(
   const hier::Box& box,
   tbox::Dimension::dir_t axis)
{
   return FaceIterator(box, axis, false);
}

}
}
