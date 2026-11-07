/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellOverlap.h"

#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace pdat {

CellOverlap::CellOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.empty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

CellOverlap::~CellOverlap()
{
}

bool
CellOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
CellOverlap::getDestinationBoxContainer() const
{
   return d_dst_boxes;
}

void
CellOverlap::getSourceBoxContainer(hier::BoxContainer& src_boxes) const
{
   TBOX_ASSERT(src_boxes.empty());

   src_boxes = d_dst_boxes;
   if (!src_boxes.empty()) {
      for (hier::BoxContainer::iterator bi = src_boxes.begin();
           bi != src_boxes.end(); ++bi) {
         d_transformation.inverseTransform(*bi);
      }
   }
}

const hier::IntVector&
CellOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
CellOverlap::getTransformation() const
{
   return d_transformation;
}

void
CellOverlap::print(
   std::ostream& os) const
{
   os << "CellOverlap boxes:";
   for (hier::BoxContainer::const_iterator b = d_dst_boxes.begin();
        b != d_dst_boxes.end(); ++b) {
      const hier::Box& box = *b;
      os << "  " << box << std::endl;
   }
}

}
}
