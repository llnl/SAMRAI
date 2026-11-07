/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/pdat/SecondLayerNodeVariableFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

const std::string SecondLayerNodeVariableFillPattern::s_name_id =
   "SECOND_LAYER_NODE_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

SecondLayerNodeVariableFillPattern::SecondLayerNodeVariableFillPattern(
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

SecondLayerNodeVariableFillPattern::~SecondLayerNodeVariableFillPattern()
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
SecondLayerNodeVariableFillPattern::calculateOverlap(
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

   const tbox::Dimension dim(dst_patch_box.getDim());

   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_patch_box));
   hier::Box src_node_mask(NodeGeometry::toNodeBox(src_mask));

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, dst_patch_box);

   hier::BoxContainer dst_boxes;

   const NodeGeometry* t_dst = CPP_CAST<const NodeGeometry *>(&dst_geometry);
   const NodeGeometry* t_src = CPP_CAST<const NodeGeometry *>(&src_geometry);

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
SecondLayerNodeVariableFillPattern::getStencilWidth()
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
SecondLayerNodeVariableFillPattern::getPatternName() const
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
SecondLayerNodeVariableFillPattern::computeStencilBoxes(
   hier::BoxContainer& stencil_boxes,
   const hier::Box& dst_box) const
{
   TBOX_ASSERT(stencil_boxes.size() == 0);

   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_box));

   hier::Box ghost_box(dst_node_box);
   ghost_box.grow(hier::IntVector::getOne(dst_box.getDim()));
   stencil_boxes.removeIntersections(ghost_box, dst_node_box);
   stencil_boxes.coalesce();
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
SecondLayerNodeVariableFillPattern::computeFillBoxesOverlap(
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
      b->growUpper(hier::IntVector::getOne(patch_box.getDim()));
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
