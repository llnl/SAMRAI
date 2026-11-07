/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real edge-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataNormOpsReal_C
#define included_math_PatchEdgeDataNormOpsReal_C

#include "SAMRAI/math/PatchEdgeDataNormOpsReal.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/EdgeGeometry.h"

#include <cmath>

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchEdgeDataNormOpsReal<TYPE>::PatchEdgeDataNormOpsReal()
{
}

template<class TYPE>
PatchEdgeDataNormOpsReal<TYPE>::~PatchEdgeDataNormOpsReal()
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
PatchEdgeDataNormOpsReal<TYPE>::numberOfEntries(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   size_t retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box dbox = pdat::EdgeGeometry::toEdgeBox(ibox, d);
      retval += (dbox.size() * data->getDepth());
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real edge-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchEdgeDataNormOpsReal<TYPE>::sumControlVolumes(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
            cvol->getArrayData(d),
            edge_box);
   }
   return retval;
}

template<class TYPE>
void
PatchEdgeDataNormOpsReal<TYPE>::abs(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& dst,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   tbox::Dimension::dir_t dimVal = box.getDim().getValue();

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         edge_box);
   }
}

template<class TYPE>
double
PatchEdgeDataNormOpsReal<TYPE>::L1Norm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), edge_box);
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

template<class TYPE>
double
PatchEdgeDataNormOpsReal<TYPE>::L2Norm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
               data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchEdgeDataNormOpsReal<TYPE>::weightedL2Norm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
               weight->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
               data->getArrayData(d),
               weight->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchEdgeDataNormOpsReal<TYPE>::RMSNorm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
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
PatchEdgeDataNormOpsReal<TYPE>::weightedRMSNorm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
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
PatchEdgeDataNormOpsReal<TYPE>::maxNorm(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNorm(data->getArrayData(d), edge_box));
      }
   } else {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNormWithControlVolume(
                  data->getArrayData(d),
                  cvol->getArrayData(d),
                  edge_box));
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchEdgeDataNormOpsReal<TYPE>::dot(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data1,
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   tbox::Dimension::dir_t dimVal = data1->getDim().getValue();

   TYPE retval = 0.0;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
               data2->getArrayData(d),
               edge_box);
      }
   } else {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
               data1->getArrayData(d),
               data2->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchEdgeDataNormOpsReal<TYPE>::integral(
   const std::shared_ptr<pdat::EdgeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::EdgeData<double> >& vol) const
{
   TBOX_ASSERT(data);

   tbox::Dimension::dir_t dimVal = data->getDim().getValue();

   TYPE retval = 0.0;

   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box side_box = pdat::EdgeGeometry::toEdgeBox(box, d);
      retval += d_array_ops.integral(
            data->getArrayData(d),
            vol->getArrayData(d),
            side_box);
   }

   return retval;
}

}
}
#endif
