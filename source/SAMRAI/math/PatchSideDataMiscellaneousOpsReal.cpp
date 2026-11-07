/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real side-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataMiscellaneousOpsReal_C
#define included_math_PatchSideDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchSideDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataMiscellaneousOpsReal<TYPE>::PatchSideDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchSideDataMiscellaneousOpsReal<TYPE>::~PatchSideDataMiscellaneousOpsReal()
{
}

/*
 *************************************************************************
 *
 * Templated miscellaneous opertions for real side-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
int
PatchSideDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const std::shared_ptr<pdat::SideData<TYPE> >& data1,
   const std::shared_ptr<pdat::SideData<TYPE> >& data2,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector());

   int retval = 1;
   tbox::Dimension::dir_t dimVal = data1->getDim().getValue();

   const hier::IntVector& directions = data1->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<int>::Min(retval,
                  d_array_ops.computeConstrProdPos(
                     data1->getArrayData(d),
                     data2->getArrayData(d),
                     side_box));
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<int>::Min(retval,
                  d_array_ops.computeConstrProdPosWithControlVolume(
                     data1->getArrayData(d),
                     data2->getArrayData(d),
                     cvol->getArrayData(d),
                     side_box));
         }
      }
   }
   return retval;
}

template<class TYPE>
void
PatchSideDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());

   tbox::Dimension::dir_t dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            d_array_ops.compareToScalar(dst->getArrayData(d),
               src->getArrayData(d),
               alpha,
               side_box);
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
               src->getArrayData(d),
               alpha,
               cvol->getArrayData(d),
               side_box);
         }
      }
   }
}

template<class TYPE>
int
PatchSideDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   const std::shared_ptr<pdat::SideData<TYPE> >& dst,
   const std::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box,
   const std::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());

   tbox::Dimension::dir_t dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   int retval = 1;
   if (!cvol) {
      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<int>::Min(retval,
                  d_array_ops.testReciprocal(
                     dst->getArrayData(d),
                     src->getArrayData(d),
                     side_box));
         }
      }
   } else {
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<int>::Min(retval,
                  d_array_ops.testReciprocalWithControlVolume(
                     dst->getArrayData(d),
                     src->getArrayData(d),
                     cvol->getArrayData(d),
                     side_box));
         }
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchSideDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const std::shared_ptr<pdat::SideData<TYPE> >& numer,
   const std::shared_ptr<pdat::SideData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   tbox::Dimension::dir_t dimVal = numer->getDim().getValue();

   TYPE retval = 0.0;
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box side_box =
         pdat::SideGeometry::toSideBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            side_box);
      retval = tbox::MathUtilities<TYPE>::Max(retval, dirval);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchSideDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const std::shared_ptr<pdat::SideData<TYPE> >& numer,
   const std::shared_ptr<pdat::SideData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   tbox::Dimension::dir_t dimVal = numer->getDim().getValue();

   TYPE retval = 0.0;
   for (tbox::Dimension::dir_t d = 0; d < dimVal; ++d) {
      const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            side_box);
      retval = tbox::MathUtilities<TYPE>::Min(retval, dirval);
   }
   return retval;
}

}
}
#endif
