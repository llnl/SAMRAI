/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class that fills patch interiors only
 *
 ************************************************************************/
#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

const std::string PatchInteriorVariableFillPattern::s_name_id =
   "PATCH_INTERIOR_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Default contructor only sets the string name identifier
 *
 *************************************************************************
 */

PatchInteriorVariableFillPattern::PatchInteriorVariableFillPattern(
   const tbox::Dimension& dim):
   d_dim(dim)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchInteriorVariableFillPattern::~PatchInteriorVariableFillPattern()
{
}

/*
 *************************************************************************
 *
 * Calculate the overlap using the implemented calculateOverlap() method
 * for the destination geometry.
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
PatchInteriorVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   NULL_USE(dst_patch_box);
   NULL_USE(overwrite_interior);
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);
   hier::BoxContainer dst_restrict_boxes(dst_patch_box);
   return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
      true, transformation,
      dst_restrict_boxes);
}

/*
 *************************************************************************
 *
 * Compute BoxOverlap that specifies data to be filled by refinement
 * operator.
 *
 *************************************************************************
 */
std::shared_ptr<hier::BoxOverlap>
PatchInteriorVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::BoxContainer& unfilled_node_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(unfilled_node_boxes);

   /*
    * For this case, the overlap is simply the intersection of
    * fill_boxes, data_box, and patch_box.
    */
   hier::Transformation transformation(
      hier::IntVector::getZero(patch_box.getDim()));

   hier::BoxContainer overlap_boxes(fill_boxes);
   overlap_boxes.intersectBoxes(data_box * patch_box);

   return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes,
      transformation);
}

const hier::IntVector&
PatchInteriorVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getZero(d_dim);
}

const std::string&
PatchInteriorVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

}
}
