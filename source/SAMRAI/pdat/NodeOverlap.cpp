/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/NodeOverlap.h"

namespace SAMRAI {
namespace pdat {

NodeOverlap::NodeOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.empty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

NodeOverlap::~NodeOverlap()
{
}

bool
NodeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
NodeOverlap::getDestinationBoxContainer() const
{
   return d_dst_boxes;
}

void
NodeOverlap::getSourceBoxContainer(hier::BoxContainer& src_boxes) const
{
   TBOX_ASSERT(src_boxes.empty());

   src_boxes = d_dst_boxes;
   if (!src_boxes.empty()) {
      const tbox::Dimension& dim = src_boxes.front().getDim();
      for (hier::BoxContainer::iterator bi = src_boxes.begin();
           bi != src_boxes.end(); ++bi) {
         bi->setUpper(bi->upper() - hier::IntVector::getOne(dim));
         d_transformation.inverseTransform(*bi);
         bi->setUpper(bi->upper() + hier::IntVector::getOne(dim));
      }
   }
}

const hier::IntVector&
NodeOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
NodeOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
