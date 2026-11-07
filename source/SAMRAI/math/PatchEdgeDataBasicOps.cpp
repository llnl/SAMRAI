/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated edge-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataBasicOps_C
#define included_math_PatchEdgeDataBasicOps_C

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/math/PatchEdgeDataBasicOps.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchEdgeDataBasicOps<TYPE>::PatchEdgeDataBasicOps()
{
}

template<class TYPE>
PatchEdgeDataBasicOps<TYPE>::~PatchEdgeDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * General basic templated operations for edge data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::scale(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
         alpha, src->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::addScalar(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
         src->getArrayData(d), alpha,
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::add(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::subtract(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::multiply(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::divide(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::reciprocal(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
         src->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::linearSum(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const TYPE& beta,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         beta, src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::axpy(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::axmy(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
void
PatchEdgeDataBasicOps<TYPE>::setRandomValues(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   int dimVal = box.getDim().getValue();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
         width, low, edge_box);
   }
}

template<class TYPE>
TYPE
PatchEdgeDataBasicOps<TYPE>::min(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
            minval, d_array_ops.min(data->getArrayData(d), edge_box));
   }
   return minval;
}

template<class TYPE>
TYPE
PatchEdgeDataBasicOps<TYPE>::max(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (int d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
            maxval, d_array_ops.max(data->getArrayData(d), edge_box));
   }
   return maxval;
}

}
}
#endif
