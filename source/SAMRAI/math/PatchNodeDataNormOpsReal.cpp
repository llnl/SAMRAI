/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real node-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchNodeDataNormOpsReal_C
#define included_math_PatchNodeDataNormOpsReal_C

#include "SAMRAI/math/PatchNodeDataNormOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

#include <cmath>

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchNodeDataNormOpsReal<TYPE>::PatchNodeDataNormOpsReal()
{
}

template<class TYPE>
PatchNodeDataNormOpsReal<TYPE>::~PatchNodeDataNormOpsReal()
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
PatchNodeDataNormOpsReal<TYPE>::numberOfEntries(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   const hier::Box ibox =
      pdat::NodeGeometry::toNodeBox(box * data->getGhostBox());
   size_t retval = ibox.size() * data->getDepth();
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real node-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchNodeDataNormOpsReal<TYPE>::sumControlVolumes(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const std::shared_ptr<pdat::NodeData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      node_box);
}

template<class TYPE>
void
PatchNodeDataNormOpsReal<TYPE>::abs(
   const std::shared_ptr<pdat::NodeData<TYPE> >& dst,
   const std::shared_ptr<pdat::NodeData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      node_box);
}

template<class TYPE>
double
PatchNodeDataNormOpsReal<TYPE>::L1Norm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.L1Norm(data->getArrayData(), node_box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
double
PatchNodeDataNormOpsReal<TYPE>::L2Norm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.L2Norm(data->getArrayData(), node_box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
double
PatchNodeDataNormOpsReal<TYPE>::weightedL2Norm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const std::shared_ptr<pdat::NodeData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
            weight->getArrayData(),
            node_box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.weightedL2NormWithControlVolume(
            data->getArrayData(),
            weight->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
double
PatchNodeDataNormOpsReal<TYPE>::RMSNorm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
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
PatchNodeDataNormOpsReal<TYPE>::weightedRMSNorm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const std::shared_ptr<pdat::NodeData<TYPE> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
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
PatchNodeDataNormOpsReal<TYPE>::maxNorm(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   double retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.maxNorm(data->getArrayData(), node_box);
   } else {
      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchNodeDataNormOpsReal<TYPE>::dot(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data1,
   const std::shared_ptr<pdat::NodeData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);
   if (!cvol) {
      retval = d_array_ops.dot(data1->getArrayData(),
            data2->getArrayData(),
            node_box);
   } else {
      retval = d_array_ops.dotWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            node_box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchNodeDataNormOpsReal<TYPE>::integral(
   const std::shared_ptr<pdat::NodeData<TYPE> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& vol) const
{
   TBOX_ASSERT(data);

   TYPE retval;
   const hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);

   retval = d_array_ops.integral(
         data->getArrayData(),
         vol->getArrayData(),
         node_box);

   return retval;
}

}
}
#endif
