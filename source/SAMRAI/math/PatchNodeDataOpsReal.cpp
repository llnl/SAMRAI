/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real node-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataOpsReal_C
#define included_math_PatchNodeDataOpsReal_C

#include "SAMRAI/math/PatchNodeDataOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataOpsReal<TYPE>::PatchNodeDataOpsReal()
{
}

#if 0
/*
 * This was moved into the header due to what looks like bug in the
 * XLC compiler.
 */

template<class TYPE>
PatchNodeDataOpsReal<TYPE>::~PatchNodeDataOpsReal()
{
}
#endif

/*
 *************************************************************************
 *
 * General templated operations for real node-centered patch data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchNodeDataOpsReal<TYPE>::swapData(
   const std::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   std::shared_ptr<pdat::NodeData<TYPE> > d1(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
         patch->getPatchData(data1_id)));
   std::shared_ptr<pdat::NodeData<TYPE> > d2(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<TYPE>, hier::PatchData>(
         patch->getPatchData(data2_id)));

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

template<class TYPE>
void
PatchNodeDataOpsReal<TYPE>::printData(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
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
PatchNodeDataOpsReal<TYPE>::copyData(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   (dst->getArrayData()).copy(src->getArrayData(), node_box);
}

template<class TYPE>
void
PatchNodeDataOpsReal<TYPE>::setToScalar(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
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
