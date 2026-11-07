/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex node-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchNodeDataNormOpsComplex.h"
#include "SAMRAI/pdat/NodeGeometry.h"

#include <cmath>

namespace SAMRAI {
namespace math {

PatchNodeDataNormOpsComplex::PatchNodeDataNormOpsComplex()
{
}

PatchNodeDataNormOpsComplex::~PatchNodeDataNormOpsComplex()
{
}

double
PatchNodeDataNormOpsComplex::L1Norm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
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

double
PatchNodeDataNormOpsComplex::L2Norm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
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

double
PatchNodeDataNormOpsComplex::weightedL2Norm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const std::shared_ptr<pdat::NodeData<dcomplex> >& weight,
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

double
PatchNodeDataNormOpsComplex::RMSNorm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
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

double
PatchNodeDataNormOpsComplex::weightedRMSNorm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
   const std::shared_ptr<pdat::NodeData<dcomplex> >& weight,
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

double
PatchNodeDataNormOpsComplex::maxNorm(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data,
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

dcomplex
PatchNodeDataNormOpsComplex::dot(
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data1,
   const std::shared_ptr<pdat::NodeData<dcomplex> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::NodeData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   dcomplex retval;
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

}
}
