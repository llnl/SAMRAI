/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/EdgeOverlap.h"

#include "SAMRAI/pdat/EdgeGeometry.h"

namespace SAMRAI {
namespace pdat {

EdgeOverlap::EdgeOverlap(
   const std::vector<hier::BoxContainer>& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(true),
   d_transformation(transformation)
{
   const tbox::Dimension dim(transformation.getOffset().getDim());
   d_dst_boxes.resize(boxes.size());

   for (int d = 0; d < static_cast<int>(boxes.size()); ++d) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].empty()) d_is_overlap_empty = false;
   }
}

EdgeOverlap::~EdgeOverlap()
{
}

bool
EdgeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
EdgeOverlap::getDestinationBoxContainer(
   const int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < static_cast<int>(d_dst_boxes.size())));

   return d_dst_boxes[axis];
}

void
EdgeOverlap::getSourceBoxContainer(hier::BoxContainer& src_boxes,
                                   int& axis_direction) const
{
   TBOX_ASSERT(src_boxes.empty());
   TBOX_ASSERT(axis_direction >= 0 &&
      axis_direction < static_cast<int>(d_dst_boxes.size()));

   src_boxes = d_dst_boxes[axis_direction];
   int transform_direction = axis_direction;
   if (!src_boxes.empty()) {
      hier::Transformation inverse_transform =
         d_transformation.getInverseTransformation();
      for (hier::BoxContainer::iterator bi = src_boxes.begin();
           bi != src_boxes.end(); ++bi) {
         transform_direction = axis_direction;
         EdgeGeometry::transform(*bi,
            transform_direction,
            inverse_transform);
      }
   }

   axis_direction = transform_direction;

   TBOX_ASSERT(axis_direction >= 0 &&
      axis_direction < static_cast<int>(d_dst_boxes.size()));

}

const hier::IntVector&
EdgeOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
EdgeOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
