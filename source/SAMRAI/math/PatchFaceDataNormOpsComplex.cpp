/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex face-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchFaceDataNormOpsComplex.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cmath>

namespace SAMRAI {
namespace math {

PatchFaceDataNormOpsComplex::PatchFaceDataNormOpsComplex()
{
}

PatchFaceDataNormOpsComplex::~PatchFaceDataNormOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int
PatchFaceDataNormOpsComplex::numberOfEntries(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();
   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const int data_depth = data->getDepth();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      retval += static_cast<int>((pdat::FaceGeometry::toFaceBox(ibox, d).size()) * data_depth);
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex face-centered data.
 *
 *************************************************************************
 */

double
PatchFaceDataNormOpsComplex::sumControlVolumes(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const std::shared_ptr<pdat::FaceData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();
   double retval = 0.0;
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
            cvol->getArrayData(d),
            pdat::FaceGeometry::toFaceBox(box, d));
   }
   return retval;
}

void
PatchFaceDataNormOpsComplex::abs(
   const std::shared_ptr<pdat::FaceData<double> >& dst,
   const std::shared_ptr<pdat::FaceData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         pdat::FaceGeometry::toFaceBox(box, d));
   }
}

double
PatchFaceDataNormOpsComplex::L1Norm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

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

double
PatchFaceDataNormOpsComplex::L2Norm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

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

double
PatchFaceDataNormOpsComplex::weightedL2Norm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const std::shared_ptr<pdat::FaceData<dcomplex> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

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

double
PatchFaceDataNormOpsComplex::RMSNorm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
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

double
PatchFaceDataNormOpsComplex::weightedRMSNorm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const std::shared_ptr<pdat::FaceData<dcomplex> >& weight,
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

double
PatchFaceDataNormOpsComplex::maxNorm(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

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
                  data->getArrayData(d), cvol->getArrayData(d), face_box));
      }
   }
   return retval;
}

dcomplex
PatchFaceDataNormOpsComplex::dot(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data1,
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);
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

dcomplex
PatchFaceDataNormOpsComplex::integral(
   const std::shared_ptr<pdat::FaceData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::FaceData<double> >& vol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();
   dcomplex retval = dcomplex(0.0, 0.0);
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      retval += d_array_ops.integral(data->getArrayData(d),
            vol->getArrayData(d),
            pdat::FaceGeometry::toFaceBox(box, d));
   }
   return retval;
}

}
}
