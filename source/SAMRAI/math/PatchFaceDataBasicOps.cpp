/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated face-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataBasicOps_C
#define included_math_PatchFaceDataBasicOps_C

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/math/PatchFaceDataBasicOps.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataBasicOps<TYPE>::PatchFaceDataBasicOps()
{
}

template<class TYPE>
PatchFaceDataBasicOps<TYPE>::~PatchFaceDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * General basic templated opertions for face data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::scale(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.scale(dst->getArrayData(d),
         alpha, src->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::addScalar(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.addScalar(dst->getArrayData(d),
         src->getArrayData(d), alpha,
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::add(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.add(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::subtract(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.subtract(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::multiply(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.multiply(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::divide(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.divide(dst->getArrayData(d),
         src1->getArrayData(d), src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::reciprocal(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.reciprocal(dst->getArrayData(d),
         src->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::linearSum(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const TYPE& beta,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.linearSum(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         beta, src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::axpy(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.axpy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::axmy(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.axmy(dst->getArrayData(d),
         alpha, src1->getArrayData(d),
         src2->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
void
PatchFaceDataBasicOps<TYPE>::setRandomValues(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, box);

   int dimVal = dst->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.setRandomValues(dst->getArrayData(d),
         width, low, face_box);
   }
}

template<class TYPE>
TYPE
PatchFaceDataBasicOps<TYPE>::min(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      minval = tbox::MathUtilities<TYPE>::Min(
            minval, d_array_ops.min(data->getArrayData(d), face_box));
   }
   return minval;
}

template<class TYPE>
TYPE
PatchFaceDataBasicOps<TYPE>::max(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      maxval = tbox::MathUtilities<TYPE>::Max(
            maxval, d_array_ops.max(data->getArrayData(d), face_box));
   }
   return maxval;
}

}
}
#endif
