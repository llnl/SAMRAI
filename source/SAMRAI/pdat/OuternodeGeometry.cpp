/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/OuternodeGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a outernode geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

OuternodeGeometry::OuternodeGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)

{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

OuternodeGeometry::~OuternodeGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two outernode centered
 * box geometries.  The calculateOverlap() checks whether both arguments
 * are outernode geometries; if so, it computes the intersection.  If
 * not, then it calls calculateOverlap() on the source object (if retry
 * is true) to allow the source a chance to calculate the intersection.
 * See the hier::BoxGeometry base class for more information about the
 * protocol.  A pointer to null is returned if the intersection canot be
 * computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuternodeGeometry::calculateOverlap(
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

   const NodeGeometry* t_dst_node =
      dynamic_cast<const NodeGeometry *>(&dst_geometry);
   const OuternodeGeometry* t_dst_onode =
      dynamic_cast<const OuternodeGeometry *>(&dst_geometry);
   const NodeGeometry* t_src_node =
      dynamic_cast<const NodeGeometry *>(&src_geometry);
   const OuternodeGeometry* t_src_onode =
      dynamic_cast<const OuternodeGeometry *>(&src_geometry);

   std::shared_ptr<hier::BoxOverlap> over;
   if ((t_src_onode != 0) && (t_dst_node != 0)) {
      over = doOverlap(*t_dst_node, *t_src_onode, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if ((t_dst_onode != 0) && (t_src_node != 0)) {
      over = doOverlap(*t_dst_onode, *t_src_node, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if ((t_src_onode != 0) && (t_dst_onode != 0)) {
      over = doOverlap(*t_dst_onode, *t_src_onode, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if (retry) {
      over = src_geometry.calculateOverlap(
            dst_geometry, src_geometry, src_mask, fill_box,
            overwrite_interior, transformation, false, dst_restrict_boxes);
   }
   return over;
}

/*
 *************************************************************************
 *
 * Compute the overlap between a node geometry destination box and an
 * outernode geometry source box.  The intersection algorithm is similar
 * the node geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuternodeGeometry::doOverlap(
   const NodeGeometry& dst_geometry,
   const OuternodeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);

   const tbox::Dimension& dim(src_mask.getDim());

   hier::BoxContainer dst_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_box_shifted(src_box);
   transformation.transform(src_box_shifted);
   const hier::Box dst_box(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   const hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_box));
   const hier::Box src_node_box(NodeGeometry::toNodeBox(src_box_shifted));

   // Compute the intersection (if any) for each of the side directions

   if (dst_node_box.intersects(src_node_box)) {

      const hier::Box msk_node_box(
         NodeGeometry::toNodeBox(hier::Box::shift(src_mask, src_offset)));
      const hier::Box fill_node_box(
         NodeGeometry::toNodeBox(fill_box));

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         hier::Box trimmed_src_node_box = src_node_box;
         for (tbox::Dimension::dir_t dh = static_cast<tbox::Dimension::dir_t>(d + 1);
              dh < dim.getValue();
              ++dh) {
            /*
             * For directions higher than d, narrow the box down to avoid
             * representing edge and corner nodes multiple times.
             */
            trimmed_src_node_box.setLower(dh,
               trimmed_src_node_box.lower(dh) + 1);
            trimmed_src_node_box.setUpper(dh,
               trimmed_src_node_box.upper(dh) - 1);
         }

         // Add lower side intersection (if any) to the box list
         hier::Box low_node_box(trimmed_src_node_box);
         low_node_box.setUpper(d, low_node_box.lower(d));

         hier::Box low_overlap(low_node_box * msk_node_box * dst_node_box
                               * fill_node_box);
         if (!low_overlap.empty()) {
            dst_boxes.pushBack(low_overlap);
         }

         // Add upper side intersection (if any) to the box list
         hier::Box hig_node_box(trimmed_src_node_box);
         hig_node_box.setLower(d, hig_node_box.upper(d));

         hier::Box hig_overlap(hig_node_box * msk_node_box * dst_node_box
                               * fill_node_box);
         if (!hig_overlap.empty()) {
            dst_boxes.pushBack(hig_overlap);
         }

         // Take away the interior if over_write interior is not set

         if (!overwrite_interior) {
            dst_boxes.removeIntersections(
               NodeGeometry::toNodeBox(dst_geometry.getBox()));
         }

      }  // loop over dim

      if (!dst_restrict_boxes.empty() && !dst_boxes.empty()) {
         hier::BoxContainer node_restrict_boxes;
         for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
              b != dst_restrict_boxes.end(); ++b) {
            node_restrict_boxes.pushBack(NodeGeometry::toNodeBox(*b));
         }
         dst_boxes.intersectBoxes(node_restrict_boxes);
      }

   }  // src and dst boxes intersect

   // Create the outernode overlap data object using the boxes and source shift

   return std::make_shared<NodeOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Compute the overlap between an outernode geometry destination box and a
 * node geometry source box.  The intersection algorithm is similar
 * the node geometry algorithm except that only the borders of the dest
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuternodeGeometry::doOverlap(
   const OuternodeGeometry& dst_geometry,
   const NodeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);

   const tbox::Dimension& dim(src_mask.getDim());

   hier::BoxContainer src_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.getBox(),
         src_geometry.getGhosts()) * src_mask);
   hier::Box src_box_shifted(src_box);
   transformation.transform(src_box_shifted);
   const hier::Box dst_box(
      hier::Box::grow(dst_geometry.d_box, dst_geometry.d_ghosts));

   const hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_box));
   const hier::Box src_node_box(NodeGeometry::toNodeBox(src_box_shifted));

   // Compute the intersection (if any) for each of the side directions

   if (dst_node_box.intersects(src_node_box)) {

      const hier::Box msk_node_box(
         NodeGeometry::toNodeBox(hier::Box::shift(src_mask, src_offset)));
      const hier::Box fill_node_box(
         NodeGeometry::toNodeBox(fill_box));

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         hier::Box trimmed_dst_node_box(dst_node_box * fill_node_box);
         for (tbox::Dimension::dir_t dh = static_cast<tbox::Dimension::dir_t>(d + 1);
              dh < dim.getValue();
              ++dh) {
            /*
             * For directions higher than d, narrow the box down to avoid
             * representing edge and corner nodes multiple times.
             */
            trimmed_dst_node_box.setLower(dh,
               trimmed_dst_node_box.lower(dh) + 1);
            trimmed_dst_node_box.setUpper(dh,
               trimmed_dst_node_box.upper(dh) - 1);
         }

         // Add lower side intersection (if any) to the box list
         hier::Box low_node_box(trimmed_dst_node_box);
         low_node_box.setUpper(d, low_node_box.lower(d));

         hier::Box low_overlap(low_node_box * msk_node_box * src_node_box);
         if (!low_overlap.empty()) {
            src_boxes.pushBack(low_overlap);
         }

         // Add upper side intersection (if any) to the box list
         hier::Box hig_node_box(trimmed_dst_node_box);
         hig_node_box.setLower(d, hig_node_box.upper(d));

         hier::Box hig_overlap(hig_node_box * msk_node_box * src_node_box);
         if (!hig_overlap.empty()) {
            src_boxes.pushBack(hig_overlap);
         }

         // Take away the interior of over_write interior is not set

         if (!overwrite_interior) {
            src_boxes.removeIntersections(
               NodeGeometry::toNodeBox(dst_geometry.getBox()));
         }

      }  // loop over dim

      if (!dst_restrict_boxes.empty() && !src_boxes.empty()) {
         hier::BoxContainer node_restrict_boxes;
         for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
              b != dst_restrict_boxes.end(); ++b) {
            node_restrict_boxes.pushBack(NodeGeometry::toNodeBox(*b));
         }
         src_boxes.intersectBoxes(node_restrict_boxes);
      }

   }  // src and dst boxes intersect

   // Create the side overlap data object using the boxes and source shift

   return std::make_shared<NodeOverlap>(src_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Compute the overlap between an outernode geometry destination box and an
 * outernode geometry source box.  The intersection algorithm is similar
 * the node geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuternodeGeometry::doOverlap(
   const OuternodeGeometry& dst_geometry,
   const OuternodeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const hier::IntVector& src_offset = transformation.getOffset();
   TBOX_ASSERT_OBJDIM_EQUALITY2(src_mask, src_offset);

   const tbox::Dimension& dim(src_mask.getDim());

   hier::BoxContainer dst_boxes;

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box =
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask;
   hier::Box src_box_shifted(src_box);
   transformation.transform(src_box_shifted);
   const hier::Box dst_box =
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts());

   const hier::Box dst_node_box = NodeGeometry::toNodeBox(dst_box);
   const hier::Box src_node_box = NodeGeometry::toNodeBox(src_box_shifted);

   // Compute the intersection (if any) for each of the side directions

   if (dst_node_box.intersects(src_node_box)) {

      const hier::Box msk_node_box =
         NodeGeometry::toNodeBox(hier::Box::shift(src_mask, src_offset));
      const hier::Box fill_node_box =
         NodeGeometry::toNodeBox(fill_box);

      tbox::Dimension::dir_t dst_d, src_d;

      for (dst_d = 0; dst_d < dim.getValue(); ++dst_d) {

         hier::Box trimmed_dst_node_box(dst_node_box * fill_node_box);
         for (tbox::Dimension::dir_t dh = static_cast<tbox::Dimension::dir_t>(dst_d + 1);
              dh < dim.getValue();
              ++dh) {
            trimmed_dst_node_box.setLower(dh,
               trimmed_dst_node_box.lower(dh) + 1);
            trimmed_dst_node_box.setUpper(dh,
               trimmed_dst_node_box.upper(dh) - 1);
         }

         hier::Box lo_dst_node_box = trimmed_dst_node_box;
         lo_dst_node_box.setUpper(dst_d, lo_dst_node_box.lower(dst_d));

         hier::Box hi_dst_node_box = trimmed_dst_node_box;
         hi_dst_node_box.setLower(dst_d, hi_dst_node_box.upper(dst_d));

         for (src_d = 0; src_d < dim.getValue(); ++src_d) {

            hier::Box trimmed_src_node_box = src_node_box;
            for (tbox::Dimension::dir_t dh = static_cast<tbox::Dimension::dir_t>(src_d + 1);
                 dh < dim.getValue();
                 ++dh) {
               trimmed_src_node_box.setLower(dh,
                  trimmed_src_node_box.lower(dh) + 1);
               trimmed_src_node_box.setUpper(dh,
                  trimmed_src_node_box.upper(dh) - 1);
            }

            hier::Box lo_src_node_box = trimmed_src_node_box;
            lo_src_node_box.setUpper(src_d, lo_src_node_box.lower(src_d));

            hier::Box hi_src_node_box = trimmed_src_node_box;
            hi_src_node_box.setLower(src_d, hi_src_node_box.upper(src_d));

            hier::Box lo_lo_box(
               lo_src_node_box * msk_node_box * lo_dst_node_box);
            if (!lo_lo_box.empty()) {
               dst_boxes.pushBack(lo_lo_box);
            }

            hier::Box hi_lo_box(
               hi_src_node_box * msk_node_box * lo_dst_node_box);
            if (!hi_lo_box.empty()) {
               dst_boxes.pushBack(hi_lo_box);
            }

            hier::Box lo_hi_box(
               lo_src_node_box * msk_node_box * hi_dst_node_box);
            if (!lo_hi_box.empty()) {
               dst_boxes.pushBack(lo_hi_box);
            }

            hier::Box hi_hi_box(
               hi_src_node_box * msk_node_box * hi_dst_node_box);
            if (!hi_hi_box.empty()) {
               dst_boxes.pushBack(hi_hi_box);
            }

            // Take away the interior of over_write interior is not set

            if (!overwrite_interior) {
               dst_boxes.removeIntersections(
                  NodeGeometry::toNodeBox(dst_geometry.d_box));
            }

         }  // loop over src dim

      }  // loop over dst dim

      if (!dst_restrict_boxes.empty() && !dst_boxes.empty()) {
         hier::BoxContainer node_restrict_boxes;
         for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
              b != dst_restrict_boxes.end(); ++b) {
            node_restrict_boxes.pushBack(NodeGeometry::toNodeBox(*b));
         }
         dst_boxes.intersectBoxes(node_restrict_boxes);
      }

   }  // if src and dst boxes intersect

   // Create the side overlap data object using the boxes and source shift

   return std::make_shared<NodeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Set up a NodeOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
OuternodeGeometry::setUpOverlap(
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

}
}
