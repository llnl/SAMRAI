/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Miscellaneous templated operations for real array data
 *
 ************************************************************************/

#ifndef included_math_ArrayDataMiscellaneousOpsReal_C
#define included_math_ArrayDataMiscellaneousOpsReal_C

#include "SAMRAI/math/ArrayDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
ArrayDataMiscellaneousOpsReal<TYPE>::ArrayDataMiscellaneousOpsReal()
{
}

template<class TYPE>
ArrayDataMiscellaneousOpsReal<TYPE>::~ArrayDataMiscellaneousOpsReal()
{
}

/*
 *************************************************************************
 *
 * General templated miscellaneous operations for array data.
 *
 *************************************************************************
 */

template<class TYPE>
int
ArrayDataMiscellaneousOpsReal<TYPE>::computeConstrProdPosWithControlVolume(
   const pdat::ArrayData<TYPE>& data1,
   const pdat::ArrayData<TYPE>& data2,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(data1, data2, cvol, box);
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());

   tbox::Dimension::dir_t dimVal = data1.getDim().getValue();

   int test = 1;

   const hier::Box d1_box = data1.getBox();
   const hier::Box d2_box = data2.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d1_box * d2_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth = data1.getDepth();
      const int cvdepth = cvol.getDepth();
      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d1_w[SAMRAI::MAX_DIM_VAL];
      int d2_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d1_offset = data1.getOffset();
      const size_t d2_offset = data2.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d1_begin = d1_box.offset(ibox.lower());
      size_t d2_begin = d2_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      const TYPE* dd1 = data1.getPointer();
      const TYPE* dd2 = data2.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         int d1_counter = static_cast<int>(d1_begin);
         int d2_counter = static_cast<int>(d2_begin);
         int cv_counter = static_cast<int>(cv_begin);

         int d1_b[SAMRAI::MAX_DIM_VAL];
         int d2_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (cvd[cv_counter + i0] > 0.0) {
                  if (tbox::MathUtilities<TYPE>::Abs(dd2[d2_counter + i0]) >
                      0.0
                      && (dd1[d1_counter + i0] * dd2[d2_counter + i0] <= 0.0)
                      ) {
                     test = 0;
                  }
               }
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
                  cv_step *= cv_w[k];
               }
               d1_counter = d1_b[dim_jump - 1] + d1_step;
               d2_counter = d2_b[dim_jump - 1] + d2_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
                  cv_b[m] = cv_counter;
               }
            }

         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
         cv_begin += cv_offset;
      }

   }

   return test;
}

template<class TYPE>
int
ArrayDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const pdat::ArrayData<TYPE>& data1,
   const pdat::ArrayData<TYPE>& data2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data1, data2, box);
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());

   tbox::Dimension::dir_t dimVal = data1.getDim().getValue();

   int test = 1;

   const hier::Box d1_box = data1.getBox();
   const hier::Box d2_box = data2.getBox();
   const hier::Box ibox = box * d1_box * d2_box;

   if (!ibox.empty()) {
      const int ddepth = data1.getDepth();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d1_w[SAMRAI::MAX_DIM_VAL];
      int d2_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d1_w[i] = d1_box.numberCells(i);
         d2_w[i] = d2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d1_offset = data1.getOffset();
      const size_t d2_offset = data2.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d1_begin = d1_box.offset(ibox.lower());
      size_t d2_begin = d2_box.offset(ibox.lower());

      const TYPE* dd1 = data1.getPointer();
      const TYPE* dd2 = data2.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         int d1_counter = static_cast<int>(d1_begin);
         int d2_counter = static_cast<int>(d2_begin);

         int d1_b[SAMRAI::MAX_DIM_VAL];
         int d2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d1_b[nd] = d1_counter;
            d2_b[nd] = d2_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (tbox::MathUtilities<TYPE>::Abs(dd2[d2_counter + i0]) > 0.0
                   && (dd1[d1_counter + i0] * dd2[d2_counter + i0] <= 0.0)) {
                  test = 0;
               }
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d1_step = 1;
               int d2_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d1_step *= d1_w[k];
                  d2_step *= d2_w[k];
               }
               d1_counter = d1_b[dim_jump - 1] + d1_step;
               d2_counter = d2_b[dim_jump - 1] + d2_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d1_b[m] = d1_counter;
                  d2_b[m] = d2_counter;
               }
            }

         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
      }

   }

   return test;
}

template<class TYPE>
void
ArrayDataMiscellaneousOpsReal<TYPE>::compareToScalarWithControlVolume(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const TYPE& alpha,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src, cvol, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   const hier::Box d_box = dst.getBox();
   const hier::Box s_box = src.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * s_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth = dst.getDepth();
      const int cvdepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int s_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = dst.getOffset();
      const size_t s_offset = src.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t s_begin = s_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         int d_counter = static_cast<int>(d_begin);
         int s_counter = static_cast<int>(s_begin);
         int cv_counter = static_cast<int>(cv_begin);

         int d_b[SAMRAI::MAX_DIM_VAL];
         int s_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
            cv_b[nd] = cv_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {

               if (cvd[cv_counter + i0] > 0.0) {
                  dd[d_counter + i0] = (
                        (tbox::MathUtilities<TYPE>::Abs(sd[s_counter + i0]) >=
                         alpha)
                        ? 1.0F : 0.0F);
               }
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               s_counter = s_b[dim_jump - 1] + s_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
                  cv_b[m] = cv_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
         cv_begin += cv_offset;
      }

   }
}

template<class TYPE>
void
ArrayDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   const hier::Box d_box = dst.getBox();
   const hier::Box s_box = src.getBox();
   const hier::Box ibox = box * d_box * s_box;

   if (!ibox.empty()) {

      const int ddepth = dst.getDepth();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int s_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = dst.getOffset();
      const size_t s_offset = src.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t s_begin = s_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t s_counter = s_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int s_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            s_b[nd] = static_cast<int>(s_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[d_counter + i0] = (
                     (tbox::MathUtilities<TYPE>::Abs(sd[s_counter + i0]) >=
                      alpha)
                     ? 1.0F : 0.0F);
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               s_counter = s_b[dim_jump - 1] + s_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  s_b[m] = static_cast<int>(s_counter);
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
      }

   }
}

template<class TYPE>
int
ArrayDataMiscellaneousOpsReal<TYPE>::testReciprocalWithControlVolume(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src, cvol, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   int test = 1;

   const hier::Box d_box = dst.getBox();
   const hier::Box s_box = src.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * s_box * cv_box;

   if (!ibox.empty()) {
      const int ddepth = dst.getDepth();
      const int cvdepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int s_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = dst.getOffset();
      const size_t s_offset = src.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t s_begin = s_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();
      const double* cvd = cvol.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t s_counter = s_begin;
         size_t cv_counter = cv_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int s_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            s_b[nd] = static_cast<int>(s_counter);
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (cvd[cv_counter + i0] > 0.0) {
                  if (sd[s_counter + i0] == 0.0) {
                     test = 0;
                     dd[d_counter + i0] = 0.0;
                  } else {
                     dd[d_counter + i0] = 1.0F / sd[s_counter + i0];
                  }
               }
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }
            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               s_counter = s_b[dim_jump - 1] + s_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  s_b[m] = static_cast<int>(s_counter);
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
         cv_begin += cv_offset;
      }

   }

   return test;
}

template<class TYPE>
int
ArrayDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const hier::Box& box) const
{
// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   int test = 1;

   const hier::Box d_box = dst.getBox();
   const hier::Box s_box = src.getBox();
   const hier::Box ibox = box * d_box * s_box;

   if (!ibox.empty()) {
      const int ddepth = dst.getDepth();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int s_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         s_w[i] = s_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = dst.getOffset();
      const size_t s_offset = src.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t s_begin = s_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         int d_counter = static_cast<int>(d_begin);
         int s_counter = static_cast<int>(s_begin);

         int d_b[SAMRAI::MAX_DIM_VAL];
         int s_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = d_counter;
            s_b[nd] = s_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (sd[s_counter + i0] == 0.0) {
                  test = 0;
                  dd[d_counter + i0] = 0.0F;
               } else {
                  dd[d_counter + i0] = 1.0F / sd[s_counter + i0];
               }
            }

            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }
            if (dim_jump > 0) {
               int d_step = 1;
               int s_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  s_step *= s_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               s_counter = s_b[dim_jump - 1] + s_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = d_counter;
                  s_b[m] = s_counter;
               }
            }
         }

         d_begin += d_offset;
         s_begin += s_offset;
      }

   }

   return test;
}

template<class TYPE>
TYPE
ArrayDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const pdat::ArrayData<TYPE>& numer,
   const pdat::ArrayData<TYPE>& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(numer, denom, box);
   TBOX_ASSERT(denom.getDepth() == numer.getDepth());

   tbox::Dimension::dir_t dimVal = numer.getDim().getValue();

   TYPE max = 0.0, quot;

   const hier::Box n_box = numer.getBox();
   const hier::Box d_box = denom.getBox();
   const hier::Box ibox = box * d_box * n_box;

   if (!ibox.empty()) {
      const int ddepth = denom.getDepth();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int n_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         n_w[i] = n_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t n_offset = numer.getOffset();
      const size_t d_offset = denom.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t n_begin = n_box.offset(ibox.lower());
      size_t d_begin = d_box.offset(ibox.lower());

      const TYPE* nd = numer.getPointer();
      const TYPE* dd = denom.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         size_t n_counter = n_begin;
         size_t d_counter = d_begin;

         int n_b[SAMRAI::MAX_DIM_VAL];
         int d_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nm = 0; nm < dimVal; ++nm) {
            n_b[nm] = static_cast<int>(n_counter);
            d_b[nm] = static_cast<int>(d_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (dd[d_counter + i0] == 0.0) {
                  quot = tbox::MathUtilities<TYPE>::Abs(nd[n_counter + i0]);
               } else {
                  quot = tbox::MathUtilities<TYPE>::Abs(nd[n_counter + i0]
                        / dd[d_counter + i0]);
               }
               if (max < quot) max = quot;
            }
            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int n_step = 1;
               int d_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  n_step *= n_w[k];
                  d_step *= d_w[k];
               }
               n_counter = n_b[dim_jump - 1] + n_step;
               d_counter = d_b[dim_jump - 1] + d_step;

               for (int m = 0; m < dim_jump; ++m) {
                  n_b[m] = static_cast<int>(n_counter);
                  d_b[m] = static_cast<int>(d_counter);
               }
            }
         }

         n_begin += n_offset;
         d_begin += d_offset;
      }

   }

   return max;
}

template<class TYPE>
TYPE
ArrayDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const pdat::ArrayData<TYPE>& numer,
   const pdat::ArrayData<TYPE>& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(numer, denom, box);
   TBOX_ASSERT(denom.getDepth() == numer.getDepth());

   tbox::Dimension::dir_t dimVal = numer.getDim().getValue();

   TYPE min = tbox::MathUtilities<TYPE>::getMax();
   TYPE quot = tbox::MathUtilities<TYPE>::getMax();

   const hier::Box n_box = numer.getBox();
   const hier::Box d_box = denom.getBox();
   const hier::Box ibox = box * d_box * n_box;

   if (!ibox.empty()) {
      const int ddepth = denom.getDepth();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int n_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         n_w[i] = n_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t n_offset = numer.getOffset();
      const size_t d_offset = denom.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t n_begin = n_box.offset(ibox.lower());
      size_t d_begin = d_box.offset(ibox.lower());

      const TYPE* nd = numer.getPointer();
      const TYPE* dd = denom.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         size_t n_counter = n_begin;
         size_t d_counter = d_begin;

         int n_b[SAMRAI::MAX_DIM_VAL];
         int d_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nm = 0; nm < dimVal; ++nm) {
            n_b[nm] = static_cast<int>(n_counter);
            d_b[nm] = static_cast<int>(d_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (dd[d_counter + i0] != 0.0) {
                  quot = nd[n_counter + i0] / dd[d_counter + i0];
               }
               if (quot < min) min = quot;
            }
            int dim_jump = 0;

            for (tbox::Dimension::dir_t j = 1; j < dimVal; ++j) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int n_step = 1;
               int d_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  n_step *= n_w[k];
                  d_step *= d_w[k];
               }
               n_counter = n_b[dim_jump - 1] + n_step;
               d_counter = d_b[dim_jump - 1] + d_step;

               for (int m = 0; m < dim_jump; ++m) {
                  n_b[m] = static_cast<int>(n_counter);
                  d_b[m] = static_cast<int>(d_counter);
               }
            }
         }

         n_begin += n_offset;
         d_begin += d_offset;
      }

   }

   return min;
}

}
}
#endif
