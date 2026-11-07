/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real cell-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataMiscellaneousOpsReal_C
#define included_math_PatchCellDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchCellDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchCellDataMiscellaneousOpsReal<TYPE>::PatchCellDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchCellDataMiscellaneousOpsReal<TYPE>::~PatchCellDataMiscellaneousOpsReal()
{
}

/*
 *************************************************************************
 *
 * Templated miscellaneous operations for real cell-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
int
PatchCellDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const std::shared_ptr<pdat::CellData<TYPE> >& data1,
   const std::shared_ptr<pdat::CellData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*data1, *data2, box);

   int retval;
   if (!cvol) {
      retval = d_array_ops.computeConstrProdPos(data1->getArrayData(),
            data2->getArrayData(),
            box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*data1, *cvol);

      retval = d_array_ops.computeConstrProdPosWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
void
PatchCellDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   if (!cvol) {
      d_array_ops.compareToScalar(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, *cvol);
      d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(),
         src->getArrayData(),
         alpha,
         cvol->getArrayData(),
         box);
   }
}

template<class TYPE>
int
PatchCellDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   const std::shared_ptr<pdat::CellData<TYPE> >& dst,
   const std::shared_ptr<pdat::CellData<TYPE> >& src,
   const hier::Box& box,
   const std::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*dst, *src, box);

   int retval;
   if (!cvol) {
      retval = d_array_ops.testReciprocal(dst->getArrayData(),
            src->getArrayData(),
            box);
   } else {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*dst, *cvol);

      retval = d_array_ops.testReciprocalWithControlVolume(
            dst->getArrayData(),
            src->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchCellDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const std::shared_ptr<pdat::CellData<TYPE> >& numer,
   const std::shared_ptr<pdat::CellData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*numer, *denom, box);

   TYPE retval;
   retval = d_array_ops.maxPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         box);
   return retval;
}

template<class TYPE>
TYPE
PatchCellDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const std::shared_ptr<pdat::CellData<TYPE> >& numer,
   const std::shared_ptr<pdat::CellData<TYPE> >& denom,
   const hier::Box& box) const
{

   TBOX_ASSERT(numer && denom);
   TBOX_ASSERT_OBJDIM_EQUALITY3(*numer, *denom, box);

   TYPE retval;
   retval = d_array_ops.minPointwiseDivide(numer->getArrayData(),
         denom->getArrayData(),
         box);
   return retval;
}

}
}
#endif
