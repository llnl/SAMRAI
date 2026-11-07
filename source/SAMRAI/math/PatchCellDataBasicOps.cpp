/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated cell-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataBasicOps_C
#define included_math_PatchCellDataBasicOps_C

#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchCellDataBasicOps<TYPE>::PatchCellDataBasicOps()
{
}

template<class TYPE>
PatchCellDataBasicOps<TYPE>::~PatchCellDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * Generic templated basic operations for cell-centered patch data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::scale(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::CellData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   d_array_ops.scale(dst->getArrayData(),
      alpha, src->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::addScalar(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   d_array_ops.addScalar(dst->getArrayData(),
      src->getArrayData(), alpha,
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::add(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.add(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::subtract(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.subtract(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::multiply(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.multiply(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::divide(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.divide(dst->getArrayData(),
      src1->getArrayData(), src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::reciprocal(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   d_array_ops.reciprocal(dst->getArrayData(),
      src->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::linearSum(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const TYPE& beta,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.linearSum(dst->getArrayData(),
      alpha, src1->getArrayData(),
      beta, src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::axpy(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.axpy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::axmy(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::CellData<TYPE> >& src1,
   const std::shared_ptr<pdat::CellData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   d_array_ops.axmy(dst->getArrayData(),
      alpha, src1->getArrayData(),
      src2->getArrayData(),
      box);
}

template<class TYPE>
TYPE
PatchCellDataBasicOps<TYPE>::min(
   const std::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   return d_array_ops.min(data->getArrayData(), box);
}

template<class TYPE>
TYPE
PatchCellDataBasicOps<TYPE>::max(
   const std::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   return d_array_ops.max(data->getArrayData(), box);
}

template<class TYPE>
void
PatchCellDataBasicOps<TYPE>::setRandomValues(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   d_array_ops.setRandomValues(dst->getArrayData(),
      width, low, box);
}

}
}
#endif
