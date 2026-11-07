/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsReal_C
#define included_math_PatchSideDataOpsReal_C

#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataOpsReal<TYPE>::PatchSideDataOpsReal()
{
}

template<class TYPE>
PatchSideDataOpsReal<TYPE>::~PatchSideDataOpsReal()
{
}

/*
 *************************************************************************
 *
 * General templated operations for real side-centered patch data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchSideDataOpsReal<TYPE>::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::SideData<TYPE> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<TYPE>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::SideData<TYPE> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<TYPE>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector());
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

template<class TYPE>
void
PatchSideDataOpsReal<TYPE>::printData(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<class TYPE>
void
PatchSideDataOpsReal<TYPE>::copyData(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         (dst->getArrayData(d)).copy(src->getArrayData(d), side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataOpsReal<TYPE>::setToScalar(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
#endif
