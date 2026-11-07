/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/
#include "SAMRAI/pdat/FirstLayerEdgeVariableFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

const std::string FirstLayerEdgeVariableFillPattern::s_name_id =
   "FIRST_LAYER_EDGE_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

FirstLayerEdgeVariableFillPattern::FirstLayerEdgeVariableFillPattern(
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

FirstLayerEdgeVariableFillPattern::~FirstLayerEdgeVariableFillPattern()
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
FirstLayerEdgeVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

   const tbox::Dimension& dim = dst_patch_box.getDim();
   std::vector<hier::BoxContainer> stencil_boxes(dim.getValue());
   computeStencilBoxes(stencil_boxes, dst_patch_box);

   std::vector<hier::BoxContainer> dst_boxes(dim.getValue());

   const EdgeGeometry* t_dst =
      dynamic_cast<const EdgeGeometry *>(&dst_geometry);
   const EdgeGeometry* t_src =
      dynamic_cast<const EdgeGeometry *>(&src_geometry);

   TBOX_ASSERT(t_dst);
   TBOX_ASSERT(t_src);

   t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box,
      overwrite_interior, transformation);

   for (int d = 0; d < dim.getValue(); ++d) {
      dst_boxes[d].intersectBoxes(stencil_boxes[d]);
   }

   return std::make_shared<EdgeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Return the stencil width (0)
 *
 *************************************************************************
 */

const hier::IntVector&
FirstLayerEdgeVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getZero(d_dim);
}

/*
 *************************************************************************
 *
 * Return the string name identifier
 *
 *************************************************************************
 */

const std::string&
FirstLayerEdgeVariableFillPattern::getPatternName() const
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
FirstLayerEdgeVariableFillPattern::computeStencilBoxes(
   std::vector<hier::BoxContainer>& stencil_boxes,
   const hier::Box& dst_box) const
{
   const tbox::Dimension& dim = dst_box.getDim();
   TBOX_ASSERT(static_cast<int>(stencil_boxes.size()) == dim.getValue());

   for (int d = 0; d < dim.getValue(); ++d) {
      hier::Box dst_edge_box(EdgeGeometry::toEdgeBox(dst_box, d));
      hier::Box interior_edge_box(dst_edge_box);
      hier::IntVector shrink_vector(dim, -1);
      shrink_vector[d] = 0;
      interior_edge_box.grow(shrink_vector);

      stencil_boxes[d].removeIntersections(dst_edge_box, interior_edge_box);
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
FirstLayerEdgeVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::BoxContainer& node_fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(pdf);
   NULL_USE(node_fill_boxes);

   const tbox::Dimension& dim = patch_box.getDim();

   std::vector<hier::BoxContainer> stencil_boxes(dim.getValue());
   computeStencilBoxes(stencil_boxes, patch_box);

   std::vector<hier::BoxContainer> overlap_boxes(dim.getValue());
   for (int d = 0; d < dim.getValue(); ++d) {

      /*
       * This is the equivalent of converting every box in overlap_boxes
       * to a edge centering, which must be done before intersecting with
       * stencil_boxes, which is edge-centered.
       */
      for (hier::BoxContainer::const_iterator b = fill_boxes.begin();
           b != fill_boxes.end(); ++b) {
         overlap_boxes[d].pushBack(EdgeGeometry::toEdgeBox(*b, d));
      }

      overlap_boxes[d].intersectBoxes(EdgeGeometry::toEdgeBox(data_box, d));

      overlap_boxes[d].intersectBoxes(stencil_boxes[d]);

      /*
       * We need to coalesce the boxes to prevent redundant edges in the
       * overlap, which can produce erroneous results accumulation
       * communication.
       */

      overlap_boxes[d].coalesce();
   }

   return std::make_shared<EdgeOverlap>(
             overlap_boxes,
             hier::Transformation(hier::IntVector::getZero(dim)));
}

}
}
