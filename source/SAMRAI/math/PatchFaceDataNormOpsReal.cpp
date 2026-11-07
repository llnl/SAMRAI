/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataNormOpsReal_C
#define included_math_PatchFaceDataNormOpsReal_C

#include "SAMRAI/math/PatchFaceDataNormOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataNormOpsReal<TYPE>::PatchFaceDataNormOpsReal()
{
}

template<class TYPE>
PatchFaceDataNormOpsReal<TYPE>::~PatchFaceDataNormOpsReal()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
PatchFaceDataNormOpsReal<TYPE>::numberOfEntries(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   size_t retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box dbox = pdat::FaceGeometry::toFaceBox(ibox, d);
      retval += (dbox.size() * data->getDepth());
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real face-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::sumControlVolumes(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const std::shared_ptr<pdat::FaceData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
            cvol->getArrayData(d),
            face_box);
   }
   return retval;
}

template<class TYPE>
void
PatchFaceDataNormOpsReal<TYPE>::abs(
   const std::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const std::shared_ptr<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         face_box);
   }
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::L1Norm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), face_box);
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
      }
   }
   return retval;
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::L2Norm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), face_box);
         retval += aval * aval;
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
               data->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::weightedL2Norm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const std::shared_ptr<pdat::FaceData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
               weight->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
               data->getArrayData(d),
               weight->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::RMSNorm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::weightedRMSNorm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const std::shared_ptr<pdat::FaceData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchFaceDataNormOpsReal<TYPE>::maxNorm(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNorm(data->getArrayData(d), face_box));
      }
   } else {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNormWithControlVolume(
                  data->getArrayData(d),
                  cvol->getArrayData(d),
                  face_box));
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchFaceDataNormOpsReal<TYPE>::dot(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data1,
   const std::shared_ptr<pdat::FaceData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   tbox::Dimension::dir_t dimVal = data1->getDim().getValue();

   TYPE retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
               data2->getArrayData(d),
               face_box);
      }
   } else {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
               data1->getArrayData(d),
               data2->getArrayData(d),
               cvol->getArrayData(d),
               face_box);
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchFaceDataNormOpsReal<TYPE>::integral(
   const std::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& vol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   TYPE retval = 0.0;

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      retval += d_array_ops.integral(data->getArrayData(d),
            vol->getArrayData(d),
            face_box);
   }

   return retval;
}

}
}
#endif
