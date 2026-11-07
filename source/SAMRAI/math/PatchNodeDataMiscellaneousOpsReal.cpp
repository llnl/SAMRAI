/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real node-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataMiscellaneousOpsReal_C
#define included_math_PatchNodeDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchNodeDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataMiscellaneousOpsReal<TYPE>::PatchNodeDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchNodeDataMiscellaneousOpsReal<TYPE>::~PatchNodeDataMiscellaneousOpsReal()
{
}

/*
 *************************************************************************
 *
 * Templated miscellaneous opertions for real node-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
int
PatchNodeDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   int retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
            data2->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.computeConstrProdPosWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
void
PatchNodeDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      d_array_ops.compareToScalar(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         node_box);
   } else {
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         cvol->getArrayData(),
         node_box);
   }
}

template<class TYPE>
int
PatchNodeDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);

   int retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
            src->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.testReciprocalWithControlVolume(
            dst->getArrayData(),
            src->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchNodeDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const std::shared_ptr<pdat::NodeData<TYPE> >& numer,
   const std::shared_ptr<pdat::NodeData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         node_box);
   return retval;
}

template<class TYPE>
TYPE
PatchNodeDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const std::shared_ptr<pdat::NodeData<TYPE> >& numer,
   const std::shared_ptr<pdat::NodeData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         node_box);
   return retval;
}

}
}
#endif
