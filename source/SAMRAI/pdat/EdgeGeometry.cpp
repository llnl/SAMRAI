/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/EdgeIterator.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a edge geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

EdgeGeometry::EdgeGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

EdgeGeometry::~EdgeGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two edge centered box
 * geometries.  The calculateOverlap() checks whether both arguments are
 * edge geometries; if so, it compuates the intersection.  If not, then
 * it calls calculateOverlap() on the source object (if retry is true)
 * to allow the source a chance to calculate the intersection.  See the
 * hier::BoxGeometry base class for more information about the protocol.
 * A pointer to null is returned if the intersection cannot be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
EdgeGeometry::calculateOverlap(
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

   const EdgeGeometry* t_dst =
      dynamic_cast<const EdgeGeometry *>(&dst_geometry);
   const EdgeGeometry* t_src =
      dynamic_cast<const EdgeGeometry *>(&src_geometry);

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
EdgeGeometry::computeDestinationBoxes(
   std::vector<hier::BoxContainer>& dst_boxes,
   const EdgeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes) const
{
   const tbox::Dimension& dim(src_mask.getDim());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   hier::Box src_shift(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   transformation.transform(src_shift);
   hier::Box dst_ghost(d_box);
   dst_ghost.grow(d_ghosts);

   // Compute the intersection (if any) for each of the edge directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector)
      * hier::Box::grow(dst_ghost, one_vector));

   if (!quick_check.empty()) {

      for (int d = 0; d < dim.getValue(); ++d) {

         const hier::Box dst_edge(toEdgeBox(dst_ghost, d));
         const hier::Box src_edge(toEdgeBox(src_shift, d));
         const hier::Box fill_edge(toEdgeBox(fill_box, d));
         const hier::Box together(dst_edge * src_edge * fill_edge);

         if (!together.empty()) {

            if (!overwrite_interior) {
               const hier::Box int_edge(toEdgeBox(d_box, d));
               dst_boxes[d].removeIntersections(together, int_edge);
            } else {
               dst_boxes[d].pushBack(together);
            }

         }  // if (!together.empty())

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer edge_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               edge_restrict_boxes.pushBack(toEdgeBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(edge_restrict_boxes);
         }
      }  // loop over dim

   }  // if (!quick_check.empty())

}

/*
 *************************************************************************
 *
 * Convert an AMR-index space hier::Box into a edge-index space box by a
 * cyclic shift of indices.
 *
 *************************************************************************
 */

hier::Box
EdgeGeometry::toEdgeBox(
   const hier::Box& box,
   int axis)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(0 <= axis && axis < dim.getValue());

   hier::Box edge_box(dim);

   if (!box.empty()) {
      edge_box = box;
      for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
         if (axis != i) {
            edge_box.setUpper(i, edge_box.upper(i) + 1);
         }
      }
   }

   return edge_box;
}

/*
 *************************************************************************
 *
 * Compute the overlap between two edge centered boxes.  The algorithm
 * is fairly straight-forward.  First, we perform a quick-and-dirty
 * intersection to see if the boxes might overlap.  If that intersection
 * is not empty, then we need to do a better job calculating the overlap
 * for each direction.  Note that the AMR index space boxes must be
 * shifted into the edge centered space before we calculate the proper
 * intersections.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
EdgeGeometry::doOverlap(
   const EdgeGeometry& dst_geometry,
   const EdgeGeometry& src_geometry,
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

   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Set up a EdgeOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
EdgeGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   const tbox::Dimension& dim(transformation.getOffset().getDim());
   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      for (int d = 0; d < dim.getValue(); ++d) {
         hier::Box edge_box(EdgeGeometry::toEdgeBox(*b, d));
         dst_boxes[d].pushBack(edge_box);
      }
   }

   // Create the edge overlap data object using the boxes and source shift
   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Transform a box
 *
 *************************************************************************
 */

void
EdgeGeometry::transform(
   hier::Box& box,
   int& axis_direction,
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

         for (int d = 0; d < dim.getValue(); ++d) {
            if (d != axis_direction) {
               box.setUpper(static_cast<hier::Box::dir_t>(d),
                  box.upper(static_cast<hier::Box::dir_t>(d)) - 1);
            }
         }
         transformation.transform(box);
         if (dim.getValue() == 2) {
            const int rotation_num = static_cast<int>(rotation);
            if (rotation_num % 2) {
               axis_direction = (axis_direction + 1) % 2;
            }
         } else if (dim.getValue() == 3) {

            if (axis_direction == 0) {

               switch (rotation) {

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:

                     axis_direction = 0;
                     break;

                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:

                     axis_direction = 1;
                     break;

                  default:

                     axis_direction = 2;
                     break;

               }

            } else if (axis_direction == 1) {

               switch (rotation) {
                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::JUP_IDOWN_KUP:
                  case hier::Transformation::JUP_IUP_KDOWN:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::JDOWN_IUP_KUP:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::JDOWN_IDOWN_KDOWN:

                     axis_direction = 0;
                     break;

                  case hier::Transformation::IUP_JUP_KUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::IDOWN_JUP_KDOWN:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::IDOWN_JDOWN_KUP:
                  case hier::Transformation::IUP_JDOWN_KDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:

                     axis_direction = 1;
                     break;

                  default:

                     axis_direction = 2;
                     break;
               }

            } else if (axis_direction == 2) {

               switch (rotation) {
                  case hier::Transformation::KUP_IUP_JUP:
                  case hier::Transformation::KUP_JUP_IDOWN:
                  case hier::Transformation::KDOWN_JUP_IUP:
                  case hier::Transformation::KDOWN_IDOWN_JUP:
                  case hier::Transformation::KUP_JDOWN_IUP:
                  case hier::Transformation::KUP_IDOWN_JDOWN:
                  case hier::Transformation::KDOWN_IUP_JDOWN:
                  case hier::Transformation::KDOWN_JDOWN_IDOWN:

                     axis_direction = 0;
                     break;

                  case hier::Transformation::JUP_KUP_IUP:
                  case hier::Transformation::IDOWN_KUP_JUP:
                  case hier::Transformation::IUP_KDOWN_JUP:
                  case hier::Transformation::JUP_KDOWN_IDOWN:
                  case hier::Transformation::IUP_KUP_JDOWN:
                  case hier::Transformation::JDOWN_KUP_IDOWN:
                  case hier::Transformation::JDOWN_KDOWN_IUP:
                  case hier::Transformation::IDOWN_KDOWN_JDOWN:

                     axis_direction = 1;
                     break;

                  default:

                     axis_direction = 2;
                     break;

               }
            }
         }

         for (int d = 0; d < dim.getValue(); ++d) {
            if (d != axis_direction) {
               box.setUpper(static_cast<hier::Box::dir_t>(d),
                  box.upper(static_cast<hier::Box::dir_t>(d)) + 1);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Transform an EdgeIndex
 *
 *************************************************************************
 */

void
EdgeGeometry::transform(
   EdgeIndex& index,
   const hier::Transformation& transformation)
{
   const tbox::Dimension& dim = index.getDim();

   if (transformation.getRotation() == hier::Transformation::NO_ROTATE &&
       transformation.getOffset() == hier::IntVector::getZero(dim)) {
      return;
   }

   const hier::Transformation::RotationIdentifier& rotation =
      transformation.getRotation();

   const int axis_direction = index.getAxis();

   for (int i = 0; i < dim.getValue(); ++i) {
      if (i == axis_direction && index(i) >= 0) {
         ++index(i);
      }
   }

   int new_axis_direction = axis_direction;

   if (dim.getValue() == 2) {
      const int rotation_num = static_cast<int>(rotation);

      TBOX_ASSERT(rotation_num <= 3);

      if (rotation_num) {

         EdgeIndex tmp_index(dim);
         for (int r = 0; r < rotation_num; ++r) {
            tmp_index = index;
            index(0) = tmp_index(1);
            index(1) = -tmp_index(0);
         }

         new_axis_direction = (axis_direction + rotation_num) % 2;

         index.setAxis(new_axis_direction);
      }
   } else {

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
            TBOX_ERROR("EdgeGeometry::transform invalid 3D RotationIdentifier.");

      }

      new_axis_direction = index.getAxis();

   }

   for (int i = 0; i < dim.getValue(); ++i) {
      if (i == new_axis_direction && index(i) > 0) {
         --index(i);
      }
   }

   index += transformation.getOffset();
}

void
EdgeGeometry::rotateAboutAxis(EdgeIndex& index,
                              const int axis,
                              const int num_rotations)
{
   const tbox::Dimension& dim = index.getDim();
   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   EdgeIndex tmp_index(dim);
   for (int j = 0; j < num_rotations; ++j) {
      tmp_index = index;
      index(a) = tmp_index(b);
      index(b) = -tmp_index(a);
   }

   int new_axis_direction = index.getAxis();
   if (new_axis_direction != axis) {
      for (int j = 0; j < num_rotations; ++j) {
         new_axis_direction = new_axis_direction == a ? b : a;
      }
   }
   index.setAxis(new_axis_direction);
}

EdgeIterator
EdgeGeometry::begin(
   const hier::Box& box,
   int axis)
{
   return EdgeIterator(box, axis, true);
}

EdgeIterator
EdgeGeometry::end(
   const hier::Box& box,
   int axis)
{
   return EdgeIterator(box, axis, false);
}

}
}
