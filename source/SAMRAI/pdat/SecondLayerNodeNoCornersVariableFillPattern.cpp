/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/pdat/SecondLayerNodeNoCornersVariableFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

const std::string SecondLayerNodeNoCornersVariableFillPattern::s_name_id =
   "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

SecondLayerNodeNoCornersVariableFillPattern::
SecondLayerNodeNoCornersVariableFillPattern(
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

SecondLayerNodeNoCornersVariableFillPattern::~
SecondLayerNodeNoCornersVariableFillPattern()
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
SecondLayerNodeNoCornersVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);
   NULL_USE(overwrite_interior);

   hier::BoxContainer dst_boxes;

   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_patch_box));
   hier::Box src_node_mask(NodeGeometry::toNodeBox(src_mask));

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, dst_patch_box);

   const NodeGeometry* t_dst =
      CPP_CAST<const NodeGeometry *>(&dst_geometry);
   const NodeGeometry* t_src =
      CPP_CAST<const NodeGeometry *>(&src_geometry);

   TBOX_ASSERT(t_dst);
   TBOX_ASSERT(t_src);

   t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box,
      false, transformation);

   dst_boxes.intersectBoxes(stencil_boxes);

   return std::make_shared<NodeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Return the stencil width (1)
 *
 *************************************************************************
 */

const hier::IntVector&
SecondLayerNodeNoCornersVariableFillPattern::getStencilWidth()
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
SecondLayerNodeNoCornersVariableFillPattern::getPatternName() const
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
SecondLayerNodeNoCornersVariableFillPattern::computeStencilBoxes(
   hier::BoxContainer& stencil_boxes,
   const hier::Box& dst_box) const
{
   TBOX_ASSERT(stencil_boxes.size() == 0);

   const tbox::Dimension& dim = dst_box.getDim();
   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_box));

   for (unsigned short i = 0; i < dim.getValue(); ++i) {
      hier::Box low_box(dst_node_box);
      low_box.setLower(i, dst_node_box.lower(i) - 1);
      low_box.setUpper(i, low_box.lower(i));
      stencil_boxes.pushFront(low_box);

      hier::Box high_box(dst_node_box);
      high_box.setLower(i, dst_node_box.upper(i) + 1);
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
SecondLayerNodeNoCornersVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::BoxContainer& node_fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(pdf);
   const tbox::Dimension& dim = patch_box.getDim();

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, patch_box);

   hier::BoxContainer overlap_boxes(fill_boxes);

   /*
    * This is the equivalent of converting every box in overlap_boxes
    * to a node centering, which must be done before intersecting with
    * stencil_boxes, which is node-centered.
    */
   for (hier::BoxContainer::iterator b = overlap_boxes.begin();
        b != overlap_boxes.end(); ++b) {
      b->growUpper(hier::IntVector::getOne(dim));
   }

   overlap_boxes.intersectBoxes(NodeGeometry::toNodeBox(data_box));

   overlap_boxes.intersectBoxes(stencil_boxes);
   overlap_boxes.intersectBoxes(node_fill_boxes);

   overlap_boxes.coalesce();

   return std::make_shared<NodeOverlap>(
             overlap_boxes,
             hier::Transformation(hier::IntVector::getZero(dim)));
}

}
}
