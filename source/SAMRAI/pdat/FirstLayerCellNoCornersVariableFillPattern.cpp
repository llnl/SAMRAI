/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/pdat/FirstLayerCellNoCornersVariableFillPattern.h"

#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

const std::string FirstLayerCellNoCornersVariableFillPattern::s_name_id =
   "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

FirstLayerCellNoCornersVariableFillPattern::
FirstLayerCellNoCornersVariableFillPattern(
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

FirstLayerCellNoCornersVariableFillPattern::~
FirstLayerCellNoCornersVariableFillPattern()
{
}

/*
 *************************************************************************
 *
 * Calculate the overlap according to the desired pattern
 *
 *************************************************************************
 */

std::shared_ptr<hier::BoxOverlap>
FirstLayerCellNoCornersVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, dst_patch_box);

   return dst_geometry.calculateOverlap(src_geometry,
      src_mask,
      fill_box,
      overwrite_interior,
      transformation,
      stencil_boxes);
}

/*
 *************************************************************************
 *
 * Return the stencil width (1)
 *
 *************************************************************************
 */

const hier::IntVector&
FirstLayerCellNoCornersVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getOne(d_dim);
}

/*
 *************************************************************************
 *
 * Return the string name identifier
 *
 *************************************************************************
 */

const std::string&
FirstLayerCellNoCornersVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

/*
 *************************************************************************
 *
 * Compute the boxes for the stencil around a given patch box
 *
 *************************************************************************
 */
void
FirstLayerCellNoCornersVariableFillPattern::computeStencilBoxes(
   hier::BoxContainer& stencil_boxes,
   const hier::Box& dst_box) const
{
   TBOX_ASSERT(stencil_boxes.size() == 0);

   const tbox::Dimension& dim = dst_box.getDim();

   for (unsigned short i = 0; i < dim.getValue(); ++i) {
      hier::Box low_box(dst_box);
      low_box.setLower(i, dst_box.lower(i) - 1);
      low_box.setUpper(i, low_box.lower(i));
      stencil_boxes.pushFront(low_box);

      hier::Box high_box(dst_box);
      high_box.setLower(i, dst_box.upper(i) + 1);
      high_box.setUpper(i, high_box.lower(i));
      stencil_boxes.pushFront(high_box);
   }
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
FirstLayerCellNoCornersVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::BoxContainer& node_fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(pdf);
   NULL_USE(node_fill_boxes);

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, patch_box);

   hier::BoxContainer overlap_boxes(fill_boxes);
   overlap_boxes.intersectBoxes(data_box);
   overlap_boxes.intersectBoxes(stencil_boxes);

   return std::make_shared<CellOverlap>(
             overlap_boxes,
             hier::Transformation(hier::IntVector::getZero(patch_box.getDim())));
}

}
}
