/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated side-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataBasicOps_C
#define included_math_PatchSideDataBasicOps_C

#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataBasicOps<TYPE>::PatchSideDataBasicOps()
{
}

template<class TYPE>
PatchSideDataBasicOps<TYPE>::~PatchSideDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * General basic templated opertions for side data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::scale(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.scale(dst->getArrayData(d),
            alpha, src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::addScalar(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.addScalar(dst->getArrayData(d),
            src->getArrayData(d), alpha,
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::add(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.add(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::subtract(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.subtract(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::multiply(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.multiply(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::divide(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.divide(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::reciprocal(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.reciprocal(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::linearSum(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const TYPE& beta,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.linearSum(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            beta, src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::axpy(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axpy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::axmy(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::SideData<TYPE> >& src1,
   const std::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axmy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::setRandomValues(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.setRandomValues(dst->getArrayData(d),
            width, low, side_box);
      }
   }
}

template<class TYPE>
TYPE
PatchSideDataBasicOps<TYPE>::min(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         minval = tbox::MathUtilities<TYPE>::Min(
               minval, d_array_ops.min(data->getArrayData(d), side_box));
      }
   }
   return minval;
}

template<class TYPE>
TYPE
PatchSideDataBasicOps<TYPE>::max(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         maxval = tbox::MathUtilities<TYPE>::Max(
               maxval, d_array_ops.max(data->getArrayData(d), side_box));
      }
   }
   return maxval;
}

}
}
#endif
