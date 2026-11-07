/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated node-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataBasicOps_C
#define included_math_PatchNodeDataBasicOps_C

#include "SAMRAI/math/PatchNodeDataBasicOps.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataBasicOps<TYPE>::PatchNodeDataBasicOps()
{
}

template<class TYPE>
PatchNodeDataBasicOps<TYPE>::~PatchNodeDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * Generic basic templated operations for node-centered patch data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::scale(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.scale(dst->getArrayData(),
      alpha, src->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::addScalar(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.addScalar(dst->getArrayData(),
      src->getArrayData(), alpha,
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::add(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.add(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::subtract(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.subtract(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::multiply(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.multiply(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::divide(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.divide(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::reciprocal(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.reciprocal(dst->getArrayData(),
      src->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::linearSum(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const TYPE& beta,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.linearSum(dst->getArrayData(),
      alpha, src1->getArrayData(),
      beta, src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::axpy(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.axpy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::axmy(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.axmy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      node_box);
}

template<class TYPE>
TYPE
PatchNodeDataBasicOps<TYPE>::min(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.min(data->getArrayData(), node_box);
}

template<class TYPE>
TYPE
PatchNodeDataBasicOps<TYPE>::max(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.max(data->getArrayData(), node_box);
}

template<class TYPE>
void
PatchNodeDataBasicOps<TYPE>::setRandomValues(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.setRandomValues(dst->getArrayData(),
      width, low, node_box);
}

}
}
#endif
