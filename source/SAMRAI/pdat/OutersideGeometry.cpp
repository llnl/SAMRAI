/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/OutersideGeometry.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Create a side geometry object given the box and ghost cell width.
 *
 *************************************************************************
 */

OutersideGeometry::OutersideGeometry(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   d_box(box),
   d_ghosts(ghosts)

{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(ghosts.min() >= 0);

}

OutersideGeometry::~OutersideGeometry()
{
}

/*
 *************************************************************************
 *
 * Attempt to calculate the intersection between two outerside centered
 * box geometries.  The calculateOverlap() checks whether both arguments
 * are outerside geometries; if so, it compuates the intersection.  If
 * not, then it calls calculateOverlap() on the source object (if retry
 * is true) to allow the source a chance to calculate the intersection.
 * See the hier::BoxGeometry base class for more information about the
 * protocol.  A pointer to null is returned if the intersection canot be
 * computed.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OutersideGeometry::calculateOverlap(
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

   const SideGeometry* t_dst_side =
      dynamic_cast<const SideGeometry *>(&dst_geometry);
   const OutersideGeometry* t_dst_oside =
      dynamic_cast<const OutersideGeometry *>(&dst_geometry);
   const OutersideGeometry* t_src =
      dynamic_cast<const OutersideGeometry *>(&src_geometry);

   std::shared_ptr<hier::BoxOverlap> over;

   if ((t_src != 0) && (t_dst_side != 0)) {
      over = doOverlap(*t_dst_side, *t_src, src_mask, fill_box,
            overwrite_interior,
            transformation, dst_restrict_boxes);
   } else if ((t_src != 0) && (t_dst_oside != 0)) {
      over = doOverlap(*t_dst_oside, *t_src, src_mask, fill_box,
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
 * Compute the overlap between a side geometry destination box and an
 * outerside geometry source box.  The intersection algorithm is similar
 * the side geometry algorithm except that only the borders of source
 * are used in the intersection computation.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OutersideGeometry::doOverlap(
   const SideGeometry& dst_geometry,
   const OutersideGeometry& src_geometry,
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

   const tbox::Dimension& dim(src_mask.getDim());

   TBOX_ASSERT(dst_geometry.getDirectionVector() == hier::IntVector::getOne(dim));

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the side directions

   const hier::IntVector& one_vector = hier::IntVector::getOne(dim);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      hier::Box mask_shift(src_mask);
      transformation.transform(mask_shift);

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         const hier::Box dst_side(
            SideGeometry::toSideBox(dst_ghost, d));
         const hier::Box src_side(
            SideGeometry::toSideBox(src_shift, d));

         const hier::Box together(dst_side * src_side);

         if (!together.empty()) {

            const hier::Box msk_side(
               SideGeometry::toSideBox(mask_shift, d));

            const hier::Box fill_side(
               SideGeometry::toSideBox(fill_box, d));

            // Add lower side intersection (if any) to the box list
            hier::Box low_side(src_side);
            low_side.setUpper(d, low_side.lower(d)); //+ghosts;

            hier::Box low_overlap(low_side * msk_side * dst_side * fill_side);
            if (!low_overlap.empty()) {
               dst_boxes[d].pushBack(low_overlap);
            }

            // Add upper side intersection (if any) to the box list
            hier::Box hig_side(src_side);
            hig_side.setLower(d, hig_side.upper(d)); //-ghosts;

            hier::Box hig_overlap(hig_side * msk_side * dst_side * fill_side);
            if (!hig_overlap.empty()) {
               dst_boxes[d].pushBack(hig_overlap);
            }

            // Take away the interior if overwrite_interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  SideGeometry::toSideBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())

         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer side_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               side_restrict_boxes.pushBack(SideGeometry::toSideBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(side_restrict_boxes);
         }

      }  // loop over dim

   } // if (!quick_check.empty())

   // Create the side overlap data object using the boxes and source shift

   return std::make_shared<SideOverlap>(dst_boxes, transformation);
}

/*
 *************************************************************************
 *
 * Compute the overlap between two outerside centered boxes.
 * The algorithm is similar to the standard side intersection algorithm
 * except we operate only on the boundaries of the source box.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
OutersideGeometry::doOverlap(
   const OutersideGeometry& dst_geometry,
   const OutersideGeometry& src_geometry,
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

   const tbox::Dimension& dim(src_mask.getDim());

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   // Perform a quick-and-dirty intersection to see if the boxes might overlap

   const hier::Box src_box(
      hier::Box::grow(src_geometry.d_box, src_geometry.d_ghosts) * src_mask);
   hier::Box src_shift(src_box);
   transformation.transform(src_shift);
   const hier::Box dst_ghost(
      hier::Box::grow(dst_geometry.getBox(), dst_geometry.getGhosts()));

   // Compute the intersection (if any) for each of the side directions

   const hier::IntVector& one_vector = hier::IntVector::getOne(dim);

   const hier::Box quick_check(
      hier::Box::grow(src_shift, one_vector) * hier::Box::grow(dst_ghost,
         one_vector));

   if (!quick_check.empty()) {

      hier::Box mask_shift(src_mask);
      transformation.transform(mask_shift);

      for (tbox::Dimension::dir_t d = 0; d < dim.getValue(); ++d) {

         const hier::Box dst_side(
            SideGeometry::toSideBox(dst_geometry.getBox(), d));
         const hier::Box src_side(
            SideGeometry::toSideBox(src_shift, d));

         const hier::Box together(dst_side * src_side);

         if (!together.empty()) {

            const hier::Box msk_side(
               SideGeometry::toSideBox(mask_shift, d));

            const hier::Box fill_side(
               SideGeometry::toSideBox(fill_box, d));

            hier::Box low_dst_side(dst_side);
            low_dst_side.setUpper(d, low_dst_side.lower(d));
            hier::Box hig_dst_side(dst_side);
            hig_dst_side.setLower(d, hig_dst_side.upper(d));

            // Add lower side intersection (if any) to the box list
            hier::Box low_src_side(src_side);
            low_src_side.setUpper(d, low_src_side.lower(d));

            hier::Box low_low_overlap(low_src_side * msk_side
                                      * low_dst_side * fill_side);
            if (!low_low_overlap.empty()) {
               dst_boxes[d].pushBack(low_low_overlap);
            }

            hier::Box low_hig_overlap(low_src_side * msk_side
                                      * hig_dst_side * fill_side);
            if (!low_hig_overlap.empty()) {
               dst_boxes[d].pushBack(low_hig_overlap);
            }

            // Add upper side intersection (if any) to the box list
            hier::Box hig_src_side(src_side);
            hig_src_side.setLower(d, hig_src_side.upper(d));

            hier::Box hig_low_overlap(hig_src_side * msk_side
                                      * low_dst_side * fill_side);
            if (!hig_low_overlap.empty()) {
               dst_boxes[d].pushBack(hig_low_overlap);
            }

            hier::Box hig_hig_overlap(hig_src_side * msk_side
                                      * hig_dst_side * fill_side);
            if (!hig_hig_overlap.empty()) {
               dst_boxes[d].pushBack(hig_hig_overlap);
            }

            // Take away the interior if overwrite_interior is not set
            if (!overwrite_interior) {
               dst_boxes[d].removeIntersections(
                  SideGeometry::toSideBox(dst_geometry.getBox(), d));
            }

         }  // if (!together.empty())
         if (!dst_restrict_boxes.empty() && !dst_boxes[d].empty()) {
            hier::BoxContainer side_restrict_boxes;
            for (hier::BoxContainer::const_iterator b = dst_restrict_boxes.begin();
                 b != dst_restrict_boxes.end(); ++b) {
               side_restrict_boxes.pushBack(SideGeometry::toSideBox(*b, d));
            }
            dst_boxes[d].intersectBoxes(side_restrict_boxes);
         }

         dst_boxes[d].coalesce();

      }  // loop over dim

   } // if (!quick_check.empty())

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
OutersideGeometry::setUpOverlap(
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

}
}
