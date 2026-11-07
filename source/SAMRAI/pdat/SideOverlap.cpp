/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/SideOverlap.h"

#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace pdat {

SideOverlap::SideOverlap(
   const std::vector<hier::BoxContainer>& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(true),
   d_transformation(transformation)
{
   const tbox::Dimension& dim = d_transformation.getOffset().getDim();
   d_dst_boxes.resize(boxes.size());

   for (int d = 0; d < dim.getValue(); ++d) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].empty()) d_is_overlap_empty = false;
   }
}

SideOverlap::~SideOverlap()
{
}

bool
SideOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
SideOverlap::getDestinationBoxContainer(
   const int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < static_cast<int>(d_dst_boxes.size())));

   return d_dst_boxes[axis];
}

void
SideOverlap::getSourceBoxContainer(hier::BoxContainer& src_boxes,
                                   int& normal_direction) const
{
   TBOX_ASSERT(src_boxes.empty());
   TBOX_ASSERT(normal_direction >= 0 &&
      normal_direction < static_cast<int>(d_dst_boxes.size()));

   src_boxes = d_dst_boxes[normal_direction];
   int transform_normal = normal_direction;
   if (!src_boxes.empty()) {
      hier::Transformation inverse_transform =
         d_transformation.getInverseTransformation();
      for (hier::BoxContainer::iterator bi = src_boxes.begin();
           bi != src_boxes.end(); ++bi) {
         if (d_transformation.getRotation() == 0) {
            bi->setUpper(static_cast<tbox::Dimension::dir_t>(normal_direction),
               bi->upper(static_cast<tbox::Dimension::dir_t>(normal_direction)) - 1);
            inverse_transform.transform(*bi);
            bi->setUpper(static_cast<tbox::Dimension::dir_t>(normal_direction),
               bi->upper(static_cast<tbox::Dimension::dir_t>(normal_direction)) + 1);
         } else {
            transform_normal = normal_direction;
            SideGeometry::transform(*bi,
               transform_normal,
               inverse_transform);
         }
      }
   }
   normal_direction = transform_normal;

   TBOX_ASSERT(normal_direction >= 0 &&
      normal_direction < static_cast<int>(d_dst_boxes.size()));

}

const hier::IntVector&
SideOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
SideOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
