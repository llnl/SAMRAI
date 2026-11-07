/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex side-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchSideDataNormOpsComplex.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cmath>

namespace SAMRAI {
namespace math {

PatchSideDataNormOpsComplex::PatchSideDataNormOpsComplex()
{
}

PatchSideDataNormOpsComplex::~PatchSideDataNormOpsComplex()
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
PatchSideDataNormOpsComplex::numberOfEntries(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();
   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   const int data_depth = data->getDepth();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         retval +=
            static_cast<int>((pdat::SideGeometry::toSideBox(ibox, d).size()) * data_depth);
      }
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex side-centered data.
 *
 *************************************************************************
 */

double
PatchSideDataNormOpsComplex::sumControlVolumes(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const std::shared_ptr<pdat::SideData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, cvol->getDirectionVector()));

   int dimVal = box.getDim().getValue();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
               cvol->getArrayData(d),
               pdat::SideGeometry::toSideBox(box, d));
      }
   }
   return retval;
}

void
PatchSideDataNormOpsComplex::abs(
   const std::shared_ptr<pdat::SideData<double> >& dst,
   const std::shared_ptr<pdat::SideData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

double
PatchSideDataNormOpsComplex::L1Norm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1Norm(data->getArrayData(d), side_box);
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

double
PatchSideDataNormOpsComplex::L2Norm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2Norm(data->getArrayData(d), side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2NormWithControlVolume(
                  data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

double
PatchSideDataNormOpsComplex::weightedL2Norm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const std::shared_ptr<pdat::SideData<dcomplex> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, weight->getDirectionVector()));

   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                  weight->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2NormWithControlVolume(
                  data->getArrayData(d),
                  weight->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

double
PatchSideDataNormOpsComplex::RMSNorm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
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
PatchSideDataNormOpsComplex::weightedRMSNorm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const std::shared_ptr<pdat::SideData<dcomplex> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
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
PatchSideDataNormOpsComplex::maxNorm(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNorm(data->getArrayData(d), side_box));
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNormWithControlVolume(
                     data->getArrayData(d),
                     cvol->getArrayData(d), side_box));
         }
      }
   }
   return retval;
}

dcomplex
PatchSideDataNormOpsComplex::dot(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data1,
   const std::shared_ptr<pdat::SideData<dcomplex> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector());

   int dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);
   const hier::IntVector& directions = data1->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dot(data1->getArrayData(d),
                  data2->getArrayData(d),
                  side_box);
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dotWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

dcomplex
PatchSideDataNormOpsComplex::integral(
   const std::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& vol) const
{
   TBOX_ASSERT(data);

   int dimVal = box.getDim().getValue();
   dcomplex retval = dcomplex(0.0, 0.0);
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, vol->getDirectionVector()));

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         retval += d_array_ops.integral(
               data->getArrayData(d),
               vol->getArrayData(d),
               pdat::SideGeometry::toSideBox(box, d));
      }
   }
   return retval;
}

}
}
