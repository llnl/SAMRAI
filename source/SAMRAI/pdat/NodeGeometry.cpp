/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a node geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

NodeGeometry::NodeGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

NodeGeometry::~NodeGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two node centered box
 * geometries.  The calculateOverlap() checks whether both arguments are
 * node geometries; if so, it computes the intersection.  If not, then
 * it calls calculateOverlap() on the source object (if retry is true)
 * to allow the source a chance to calculate the intersection.  See the
 * hier::BoxGeometry base class for more information about the
 * protocol. A pointer to null is returned if the intersection cannot
 * be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
NodeGeometry::calculateOverlap(
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

   const NodeGeometry* t_dst =
      dynamic_cast<const NodeGeometry *>(&dst_geometry);
   const NodeGeometry* t_src =
      dynamic_cast<const NodeGeometry *>(&src_geometry);

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
NodeGeometry::computeDestinationBoxes(
   hier::BoxContainer& dst_boxes,
   const NodeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes) const
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   const hier::IntVector& src_offset(transformation.getOffset());
#endif
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);

   // Translate the source box and grow the destination box by the ghost cells

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(d_box, d_ghosts));

   // Convert the boxes into node space and compute the intersection

   const hier::Box dst_node(toNodeBox(dst_ghost));
   const hier::Box src_node(toNodeBox(src_shift));
   const hier::Box fill_node(toNodeBox(fill_box));
   const hier::Box together(dst_node * src_node * fill_node);

   if (!together.empty()) {
      if (!overwrite_interior) {
         const hier::Box int_node(toNodeBox(d_box));
         dst_boxes.removeIntersections(together, int_node);
      } else {
         dst_boxes.pushBack(together);
      }
   }

   if (!dst_restrict_boxes.empty() && !dst_boxes.empty()) {
      hier::BoxContainer node_restrict_boxes;
      for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
           b != dst_restrict_boxes.end(); ++b) {
         node_restrict_boxes.pushBack(toNodeBox(*b));
      }
      dst_boxes.intersectBoxes(node_restrict_boxes);
   }
}

/*
 *************************************************************************
 *
 * Set up a NodeOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
NodeGeometry::setUpOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation) const
{
   hier::BoxContainer dst_boxes;

   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      hier::Box node_box(NodeGeometry::toNodeBox(*b));
      dst_boxes.pushBack(node_box);
   }

   // Create the node overlap data object using the boxes and source shift
   return std::make_shared<NodeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Transform a box
 *
 *************************************************************************
 */

void
NodeGeometry::transform(
   hier::Box& box,
   const hier::Transformation& transformation)
{
   if (transformation.getRotation() == hier::Transformation::NO_ROTATE &&
       transformation.getOffset() == hier::IntVector::getZero(box.getDim())) {
      return;
   }

   if (!box.empty()) {
      box.setUpper(box.upper() - hier::IntVector::getOne(box.getDim()));
      transformation.transform(box);
      box.setUpper(box.upper() + hier::IntVector::getOne(box.getDim()));
   }
}

/*
 *************************************************************************
 *
 * Transform a NodeIndex
 *
 *************************************************************************
 */

void
NodeGeometry::transform(
   NodeIndex& index,
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
         TBOX_ERROR("NodeGeometry::transform invalid 1D RotationIdentifier.");
      }

      if (rotation_num) {
         NodeIndex tmp_index(index);
         index(0) = -tmp_index(0);
      }
   } else if (dim.getValue() == 2) {
      const int rotation_num = static_cast<int>(rotation);
      if (rotation_num > 3) {
         TBOX_ERROR("NodeGeometry::transform invalid 2D RotationIdentifier.");
      }

      if (rotation_num) {
         NodeIndex tmp_index(dim);
         for (int r = 0; r < rotation_num; ++r) {
            tmp_index = index;
            index(0) = tmp_index(1);
            index(1) = -tmp_index(0);
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
            TBOX_ERROR("NodeGeometry::transform invalid 3D RotationIdentifier.");
      }
   } else {
      TBOX_ERROR("NodeGeometry::transform implemented for 2D and 3d only.");
   }

   index += transformation.getOffset();
}

void
NodeGeometry::rotateAboutAxis(NodeIndex& index,
                              const int axis,
                              const int num_rotations)
{
   const tbox::Dimension& dim = index.getDim();
   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   NodeIndex tmp_index(dim);
   for (int j = 0; j < num_rotations; ++j) {
      tmp_index = index;
      index(a) = tmp_index(b);
      index(b) = -tmp_index(a);
   }
}

NodeIterator
NodeGeometry::begin(
   const hier::Box& box)
{
   return NodeIterator(box, true);
}

NodeIterator
NodeGeometry::end(
   const hier::Box& box)
{
   return NodeIterator(box, false);
}

}
}
