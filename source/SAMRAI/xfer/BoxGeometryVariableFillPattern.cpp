/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace xfer {

const std::string BoxGeometryVariableFillPattern::s_name_id =
   "BOX_GEOMETRY_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Default contructor only sets the string name identifier
 *
 *************************************************************************
 */

BoxGeometryVariableFillPattern::BoxGeometryVariableFillPattern()
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

BoxGeometryVariableFillPattern::~BoxGeometryVariableFillPattern()
{
}

/*
 *************************************************************************
 *
 * getStencilWidth() throws an error if called.  Only overridding
 * versions of this method in concrete subclasses should be called.
 *
 *************************************************************************
 */
const hier::IntVector&
BoxGeometryVariableFillPattern::getStencilWidth()
{
   TBOX_ERROR(
      "BoxGeometryVariableFillPattern::getStencilWidth() should not be\n"
      << "called.  This pattern creates overlaps based on\n"
      << "the BoxGeometry objects and is not restricted to a\n"
      << "specific stencil.\n");

   /*
    * Dummy return value that will never get reached.
    */
   return hier::IntVector::getZero(tbox::Dimension(1));
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
BoxGeometryVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::BoxContainer& node_fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(node_fill_boxes);

   /*
    * For this (default) case, the overlap is simply the intersection of
    * fill_boxes and data_box.
    */
   hier::Transformation transformation(
      hier::IntVector::getZero(patch_box.getDim()));

   hier::BoxContainer overlap_boxes(fill_boxes);
   overlap_boxes.intersectBoxes(data_box);

   return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes,
      transformation);
}

}
}
