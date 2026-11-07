/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/OuterfaceGeometry.h"
#include "SAMRAI/pdat/FaceGeometry.h"
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

OuterfaceGeometry::OuterfaceGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);
}

OuterfaceGeometry::~OuterfaceGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two outerface centered
 * box geometries.  The calculateOverlap() checks whether both arguments
 * are outerface geometries; if so, it compuates the intersection.  If
 * not, then it calls calculateOverlap() on the source object (if retry
 * is true) to allow the source a chance to calculate the intersection.
 * See the hier::BoxGeometry base class for more information about
 * the protocol.  A pointer to null is returned if the intersection
 * cannot be computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuterfaceGeometry::calculateOverlap(
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

   const FaceGeometry* t_dst_face =
      dynamic_cast<const FaceGeometry *>(&dst_geometry);
   const OuterfaceGeometry* t_dst_oface =
      dynamic_cast<const OuterfaceGeometry *>(&dst_geometry);
   const OuterfaceGeometry* t_src =
      dynamic_cast<const OuterfaceGeometry *>(&src_geometry);

   std::shared_ptr<hier::BoxOverlap> over;

   if ((t_src != 0) && (t_dst_face != 0)) {
      over = doOverlap(*t_dst_face, *t_src, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if ((t_src != 0) && (t_dst_oface != 0)) {
      over = doOverlap(*t_dst_oface, *t_src, src_mask, fill_box,
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
 * Compute the overlap between a face geometry destination box and an
 * outerface geometry source box.  The intersection algorithm is similar
 * the face geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuterfaceGeometry::doOverlap(
   const FaceGeometry& dst_geometry,
   const OuterfaceGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   hier::Box src_box(src_geometry.d_box);
   src_box.grow(src_geometry.d_ghosts);
   src_box = src_box * src_mask;
   transformation.transform(src_box);
   hier::Box dst_ghost(dst_geometry.getBox());
   dst_ghost.grow(dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the face directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_box, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      hier::Box mask_shift(src_mask);
      transformation.transform(mask_shift);

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         const hier::Box msk_face(
            FaceGeometry::toFaceBox(mask_shift, d));
         const hier::Box dst_face(
            FaceGeometry::toFaceBox(dst_ghost, d));
         const hier::Box src_face(
            FaceGeometry::toFaceBox(src_box, d));
         const hier::Box fill_face(
            FaceGeometry::toFaceBox(fill_box, d));

         const hier::Box together(dst_face * src_face * fill_face);

         if (!together.empty()) {

            // Add lower face intersection (if any) to the box list
            hier::Box low_face(src_face);
            low_face.setUpper(0, low_face.lower(0));  //+ghosts;

            hier::Box low_overlap(low_face * msk_face * dst_face);
            if (!low_overlap.empty()) {
               dst_boxes[d].pushBack(low_overlap);
            }

            // Add upper face intersection (if any) to the box list
            hier::Box hig_face(src_face);
            hig_face.setLower(0, hig_face.upper(0));  //-ghosts;

            hier::Box hig_overlap(hig_face * msk_face * dst_face);
            if (!hig_overlap.empty()) {
               dst_boxes[d].pushBack(hig_overlap);
            }

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  FaceGeometry::toFaceBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer face_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               face_restrict_boxes.pushBack(FaceGeometry::toFaceBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(face_restrict_boxes);
         }

      }  // loop over dim

   } // if (!quick_check.empty())

   // Create the face overlap data object using the boxes and source shift

   return std::make_shared<FaceOverlap>(dst_boxes, transformation);
}
/*
 *************************************************************************
 *
 * Compute the overlap between a face geometry destination box and an
 * outerface geometry source box.  The intersection algorithm is similar
 * the face geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OuterfaceGeometry::doOverlap(
   const OuterfaceGeometry& dst_geometry,
   const OuterfaceGeometry& src_geometry,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation,
   const hier::BoxContainer& dst_restrict_boxes)
{
   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   hier::Box src_box(src_geometry.d_box);
   src_box.grow(src_geometry.d_ghosts);
   src_box = src_box * src_mask;
   transformation.transform(src_box);
   hier::Box dst_ghost(dst_geometry.getBox());
   dst_ghost.grow(dst_geometry.getGhosts());

   // Compute the intersection (if any) for each of the face directions

   const hier::IntVector one_vector(dim, 1);

   const hier::Box quick_check(
      hier::Box::grow(src_box, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      hier::Box mask_shift(src_mask);
      transformation.transform(mask_shift);

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         const hier::Box dst_face(
            FaceGeometry::toFaceBox(dst_geometry.getBox(), d));
         const hier::Box src_face(
            FaceGeometry::toFaceBox(src_box, d));
         const hier::Box fill_face(
            FaceGeometry::toFaceBox(fill_box, d));

         const hier::Box together(dst_face * src_face * fill_face);

         if (!together.empty()) {

            const hier::Box msk_face(
               FaceGeometry::toFaceBox(mask_shift, d));

            hier::Box low_dst_face(dst_face);
            low_dst_face.setUpper(0, low_dst_face.lower(0));
            hier::Box hig_dst_face(dst_face);
            hig_dst_face.setLower(0, hig_dst_face.upper(0));

            // Add lower face intersection (if any) to the box list
            hier::Box low_src_face(src_face);
            low_src_face.setUpper(0, low_src_face.lower(0));

            hier::Box low_low_overlap(low_src_face * msk_face * low_dst_face);
            if (!low_low_overlap.empty()) {
               dst_boxes[d].pushBack(low_low_overlap);
            }

            hier::Box low_hig_overlap(low_src_face * msk_face * hig_dst_face);
            if (!low_hig_overlap.empty()) {
               dst_boxes[d].pushBack(low_hig_overlap);
            }

            // Add upper face intersection (if any) to the box list
            hier::Box hig_src_face(src_face);
            hig_src_face.setLower(0, hig_src_face.upper(0));  //-ghosts;

            hier::Box hig_low_overlap(hig_src_face * msk_face * low_dst_face);
            if (!hig_low_overlap.empty()) {
               dst_boxes[d].pushBack(hig_low_overlap);
            }

            hier::Box hig_hig_overlap(hig_src_face * msk_face * hig_dst_face);
            if (!hig_hig_overlap.empty()) {
               dst_boxes[d].pushBack(hig_hig_overlap);
            }

            // Take away the interior of over_write interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  FaceGeometry::toFaceBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer face_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               face_restrict_boxes.pushBack(FaceGeometry::toFaceBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(face_restrict_boxes);
         }

         dst_boxes[d].coalesce();

      }  // loop over dim

   } // if (!quick_check.empty())

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
OuterfaceGeometry::setUpOverlap(
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

}
}
