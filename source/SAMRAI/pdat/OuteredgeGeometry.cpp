/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Box geometry information for outeredge centered objects
 *
 ************************************************************************/
#include "SAMRAI/pdat/OuteredgeGeometry.h"

#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create an outeredge geometry object given box and ghost cell width.
 *
 *************************************************************************
 */

OuteredgeGeometry::OuteredgeGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)

{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

OuteredgeGeometry::~OuteredgeGeometry()
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
OuteredgeGeometry::calculateOverlap(
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

   const EdgeGeometry* t_dst_edge =
      dynamic_cast<const EdgeGeometry *>(&dst_geometry);
   const OuteredgeGeometry* t_dst_oedge =
      dynamic_cast<const OuteredgeGeometry *>(&dst_geometry);
   const OuteredgeGeometry* t_src =
      dynamic_cast<const OuteredgeGeometry *>(&src_geometry);

   std::shared_ptr<hier::BoxOverlap> over;

   if ((t_src != 0) && (t_dst_edge != 0)) {
      over = doOverlap(*t_dst_edge, *t_src, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if ((t_src != 0) && (t_dst_oedge != 0)) {
      over = doOverlap(*t_dst_oedge, *t_src, src_mask, fill_box,
            overwrite_interior,
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
 * Compute the overlap between an edge and an outeredge centered boxes.
 * The algorithm is similar to the standard edge intersection algorithm
 * except we operate only on the boundaries of the source box.
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
OuteredgeGeometry::doOverlap(
   const EdgeGeometry& dst_geometry,
   const OuteredgeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_box_shifted(src_box);
   transformation.transform(src_box_shifted);
   const hier::Box dst_box(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the edge directions

   const hier::IntVector one_vector(dim, 1);

   bool quick_boxes_intersect =
      (hier::Box::grow(src_box_shifted, one_vector)).intersects(
         hier::Box::grow(dst_box, one_vector));
   if (quick_boxes_intersect) {

      for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {

         const hier::Box dst_edge_box(
            EdgeGeometry::toEdgeBox(dst_box, axis));
         const hier::Box src_edge_box(
            EdgeGeometry::toEdgeBox(src_box_shifted, axis));

         bool boxes_intersect = dst_edge_box.intersects(src_edge_box);

         if (boxes_intersect) {

            const hier::Box fill_edge_box(
               EdgeGeometry::toEdgeBox(fill_box, axis));

            for (tbox::Dimension::dir_t face_normal = 0;
                 face_normal < dim.getValue();
                 ++face_normal) {

               if (face_normal != axis) {

                  for (tbox::Dimension::dir_t side = 0; side < 2; ++side) {
                     hier::Box outeredge_src_box(
                        toOuteredgeBox(src_box_shifted,
                           axis,
                           face_normal,
                           side));
                     hier::Box overlap_box(
                        outeredge_src_box * dst_edge_box * fill_edge_box);
                     if (!overlap_box.empty()) {
                        dst_boxes[axis].pushBack(
                           outeredge_src_box * dst_edge_box * fill_edge_box);
                     }
                  }

               }  // data is not defined when face_normal == axis

            }  // iterate over face normal directions

            if (!overwrite_interior) {
               const hier::Box interior_edges(
                  EdgeGeometry::toEdgeBox(dst_geometry.getBox(),
                     axis));
               dst_boxes[axis].removeIntersections(interior_edges);
            }

         }  // if source and destination edge boxes overlap in axis direction

         if (!dst_restrict_boxes.empty() && !dst_boxes[axis].empty()) {
            hier::BoxContainer edge_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               edge_restrict_boxes.pushBack(EdgeGeometry::toEdgeBox(*b, axis));
            }
            dst_boxes[axis].intersectBoxes(edge_restrict_boxes);
         }

      }  // iterate over axis directions

   }  // if quick check passes

   // Create the edge overlap data object using the boxes and source shift
   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Compute the overlap between two outeredge centered boxes.
 * The algorithm is similar to the standard edge intersection algorithm
 * except we operate only on the boundaries of the source box.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuteredgeGeometry::doOverlap(
   const OuteredgeGeometry& dst_geometry,
   const OuteredgeGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{

   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_box_shifted(src_box);
   transformation.transform(src_box_shifted);
   const hier::Box dst_box(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the edge directions

   const hier::IntVector one_vector(dim, 1);

   bool quick_boxes_intersect =
      (hier::Box::grow(src_box_shifted, one_vector)).intersects(
         hier::Box::grow(dst_box, one_vector));
   if (quick_boxes_intersect) {

      for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {

         const hier::Box dst_edge_box(
            EdgeGeometry::toEdgeBox(dst_box, axis));
         const hier::Box src_edge_box(
            EdgeGeometry::toEdgeBox(src_box_shifted, axis));

         bool boxes_intersect = dst_edge_box.intersects(src_edge_box);

         if (boxes_intersect) {

            const hier::Box fill_edge_box(
               EdgeGeometry::toEdgeBox(fill_box, axis));

            for (tbox::Dimension::dir_t src_face_normal = 0;
                 src_face_normal < dim.getValue();
                 ++src_face_normal) {

               if (src_face_normal != axis) {

                  hier::Box outeredge_src_box_lo(toOuteredgeBox(
                                                    src_box_shifted,
                                                    axis,
                                                    src_face_normal,
                                                    0));
                  hier::Box outeredge_src_box_up(toOuteredgeBox(
                                                    src_box_shifted,
                                                    axis,
                                                    src_face_normal,
                                                    1));

                  for (tbox::Dimension::dir_t dst_face_normal = 0;
                       dst_face_normal < dim.getValue();
                       ++dst_face_normal) {

                     if (dst_face_normal != axis) {

                        hier::Box outeredge_dst_box_lo(toOuteredgeBox(dst_box,
                                                          axis,
                                                          dst_face_normal,
                                                          0));
                        hier::Box outeredge_dst_box_up(toOuteredgeBox(dst_box,
                                                          axis,
                                                          dst_face_normal,
                                                          1));

                        outeredge_dst_box_lo =
                           outeredge_dst_box_lo * fill_edge_box;
                        outeredge_dst_box_up =
                           outeredge_dst_box_up * fill_edge_box;

                        hier::Box lo_lo_box(
                           outeredge_src_box_lo * outeredge_dst_box_lo);
                        if (!lo_lo_box.empty()) {
                           dst_boxes[axis].pushBack(lo_lo_box);
                        }

                        hier::Box lo_up_box(
                           outeredge_src_box_lo * outeredge_dst_box_up);
                        if (!lo_up_box.empty()) {
                           dst_boxes[axis].pushBack(lo_up_box);
                        }

                        hier::Box up_lo_box(
                           outeredge_src_box_up * outeredge_dst_box_lo);
                        if (!up_lo_box.empty()) {
                           dst_boxes[axis].pushBack(up_lo_box);
                        }

                        hier::Box up_up_box(
                           outeredge_src_box_up * outeredge_dst_box_up);
                        if (!up_up_box.empty()) {
                           dst_boxes[axis].pushBack(up_up_box);
                        }

                     }  // dst data undefined when dst_face_normal == axis

                  }  // iterate over dst face normal directions

               }  // src data undefined when src_face_normal == axis

            }  // iterate over src face normal directions

         }  // if source and destination edge boxes overlap in axis direction

         if (!overwrite_interior) {
            const hier::Box interior_edges(
               EdgeGeometry::toEdgeBox(dst_geometry.getBox(),
                  axis));
            dst_boxes[axis].removeIntersections(interior_edges);
         }

         if (!dst_restrict_boxes.empty() && !dst_boxes[axis].empty()) {
            hier::BoxContainer edge_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               edge_restrict_boxes.pushBack(EdgeGeometry::toEdgeBox(*b, axis));
            }
            dst_boxes[axis].intersectBoxes(edge_restrict_boxes);
         }

      }  // iterate over axis directions

   }  // if quick check passes

   // Create the edge overlap data object using the boxes and source shift
   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Convert an AMR-index space hier::Box into a edge-index space box
 * for an outeredge region.
 *
 *************************************************************************
 */

hier::Box
OuteredgeGeometry::toOuteredgeBox(
   const hier::Box& box,
   tbox::Dimension::dir_t axis,
   tbox::Dimension::dir_t face_normal,
   int side)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(axis < dim.getValue());
   TBOX_ASSERT(face_normal < dim.getValue());
   TBOX_ASSERT(face_normal != axis);
   TBOX_ASSERT(side == 0 || side == 1);

   hier::Box oedge_box(dim);

   /*
    * If data is defined (i.e., face_normal != axis), then
    *    1) Make an edge box for the given axis.
    *    2) Trim box as needed to avoid redundant edge indices
    *       for different face normal directions.
    *    3) Restrict box to lower or upper face for given
    *       face normal direction.
    */

   if ((face_normal != axis) && !box.empty()) {

      oedge_box = EdgeGeometry::toEdgeBox(box, axis);

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         if (d != axis) {    // do not trim in axis direction

            for (tbox::Dimension::dir_t dh = static_cast<tbox::Dimension::dir_t>(d + 1);
                 dh < dim.getValue();
                 ++dh) {                                                                                              // trim higher directions

               if (dh != axis && dh != face_normal) {
                  // do not trim in axis or face_normal direction

                  oedge_box.setLower(dh, oedge_box.lower(dh) + 1);
                  oedge_box.setUpper(dh, oedge_box.upper(dh) - 1);

               }

            }

         }

      }

      if (side == 0) {   // lower side in face normal direction
         oedge_box.setUpper(face_normal, oedge_box.lower(face_normal));
      } else {  // side == 1; upper side in face normal direction
         oedge_box.setLower(face_normal, oedge_box.upper(face_normal));
      }

   }

   return oedge_box;
}

/*
 *************************************************************************
 *
 * Set up a EdgeOverlap oject using the given boxes and offset
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
OuteredgeGeometry::setUpOverlap(
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

   // Create the edge overlap data object using the boxes and transformation
   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);

}

}
}
