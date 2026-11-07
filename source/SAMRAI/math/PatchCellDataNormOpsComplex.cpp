/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex cell-centered patch data.
 *
 ************************************************************************/
#include "SAMRAI/math/PatchCellDataNormOpsComplex.h"

#include <cmath>

namespace SAMRAI {
namespace math {

PatchCellDataNormOpsComplex::PatchCellDataNormOpsComplex()
{
}

PatchCellDataNormOpsComplex::~PatchCellDataNormOpsComplex()
{
}

double
PatchCellDataNormOpsComplex::L1Norm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.L1Norm(data->getArrayData(), box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double
PatchCellDataNormOpsComplex::L2Norm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.L2Norm(data->getArrayData(), box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double
PatchCellDataNormOpsComplex::weightedL2Norm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const std::shared_ptr<pdat::CellData<dcomplex> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
            weight->getArrayData(),
            box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.weightedL2NormWithControlVolume(
            data->getArrayData(),
            weight->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

double
PatchCellDataNormOpsComplex::RMSNorm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double
PatchCellDataNormOpsComplex::weightedRMSNorm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const std::shared_ptr<pdat::CellData<dcomplex> >& weight,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data, *weight, box);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

double
PatchCellDataNormOpsComplex::maxNorm(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*data, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.maxNorm(data->getArrayData(), box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data, *cvol);

      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

dcomplex
PatchCellDataNormOpsComplex::dot(
   const std::shared_ptr<pdat::CellData<dcomplex> >& data1,
   const std::shared_ptr<pdat::CellData<dcomplex> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data1, *data2, box);

   dcomplex retval;
   if (!cvol) {
      retval = d_array_ops.dot(data1->getArrayData(),
            data2->getArrayData(),
            box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data1, *cvol);

      retval = d_array_ops.dotWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

}
}
