/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataNormOpsReal_C
#define included_math_PatchSideDataNormOpsReal_C

#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

#include <cmath>

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataNormOpsReal<TYPE>::PatchSideDataNormOpsReal()
{
}

template<class TYPE>
PatchSideDataNormOpsReal<TYPE>::~PatchSideDataNormOpsReal()
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
PatchSideDataNormOpsReal<TYPE>::numberOfEntries(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   size_t retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box dbox = pdat::SideGeometry::toSideBox(ibox, d);
         retval += (dbox.size() * data->getDepth());
      }
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real side-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::sumControlVolumes(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const std::shared_ptr<pdat::SideData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, cvol->getDirectionVector()));

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
               cvol->getArrayData(d),
               side_box);
      }
   }
   return retval;
}

template<class TYPE>
void
PatchSideDataNormOpsReal<TYPE>::abs(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::L1Norm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

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
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

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

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::L2Norm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

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
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

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

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::weightedL2Norm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const std::shared_ptr<pdat::SideData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

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
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

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

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::RMSNorm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
// SGS

   TBOX_ASSERT(data);

   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::weightedRMSNorm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const std::shared_ptr<pdat::SideData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::maxNorm(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

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
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNormWithControlVolume(
                     data->getArrayData(d),
                     cvol->getArrayData(d),
                     side_box));
         }
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchSideDataNormOpsReal<TYPE>::dot(
   const std::shared_ptr<pdat::SideData<TYPE> >& data1,
   const std::shared_ptr<pdat::SideData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector());

   tbox::Dimension::dir_t dimVal = data1->getDim().getValue();

   TYPE retval = 0.0;
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
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data1, *cvol);
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

template<class TYPE>
TYPE
PatchSideDataNormOpsReal<TYPE>::integral(
   const std::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& vol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT(vol);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, box, *vol);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   TYPE retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, vol->getDirectionVector()));

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.integral(
               data->getArrayData(d),
               vol->getArrayData(d),
               side_box);
      }
   }

   return retval;
}

}
}
#endif
