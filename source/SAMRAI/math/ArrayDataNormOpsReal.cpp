/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated array data norm operations.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataNormOpsReal_C
#define included_math_ArrayDataNormOpsReal_C

#include "SAMRAI/math/ArrayDataNormOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cmath>

namespace SAMRAI {
namespace math {

template<class TYPE>
ArrayDataNormOpsReal<TYPE>::ArrayDataNormOpsReal()
{
}

template<class TYPE>
ArrayDataNormOpsReal<TYPE>::~ArrayDataNormOpsReal()
{
}

/*
 *************************************************************************
 *
 * Generic class template norm operations for real array data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayDataNormOpsReal<TYPE>::abs(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   const hier::Box dst_box = dst.getBox();
   const hier::Box src_box = src.getBox();
   const hier::Box ibox = box * dst_box * src_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = static_cast<int>(ibox.numberCells(i));
         dst_w[i] = static_cast<int>(dst_box.numberCells(i));
         src_w[i] = static_cast<int>(src_box.numberCells(i));
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src_offset = src.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src_begin = src_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();

      const unsigned int ddepth = dst.getDepth();
      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src_counter = src_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src_b[nd] = static_cast<int>(src_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] =
                  tbox::MathUtilities<TYPE>::Abs(sd[src_counter + i0]);
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
               int dst_step = 1;
               int src_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  dst_step *= dst_w[k];
                  src_step *= src_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src_counter = src_b[dim_jump - 1] + src_step;

               for (int m = 0; m < dim_jump; ++m) {
                  dst_b[m] = static_cast<int>(dst_counter);
                  src_b[m] = static_cast<int>(src_counter);
               }
            }
         }

         dst_begin += dst_offset;
         src_begin += src_offset;

      }
   }
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::sumControlVolumes(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, cvol, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double sum = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * cv_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t cv_offset = cvol.getOffset();

      size_t cv_begin = cv_box.offset(ibox.lower());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      const unsigned int ddepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == data.getDepth()) || (ddepth == 1));

      const double* cvd = cvol.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t cv_counter = cv_begin;

         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               sum += cvd[cv_counter + i0];
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
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  cv_step *= cv_w[k];
               }
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }
         cv_begin += cv_offset;
      }

      if (ddepth != data.getDepth()) sum *= data.getDepth();

   }

   return sum;
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::L1NormWithControlVolume(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, cvol, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double l1norm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * cv_box;

   if (!ibox.empty()) {
      const unsigned int ddepth = data.getDepth();
      const unsigned int cvdepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      size_t d_begin = d_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      const TYPE* dd = data.getPointer();
      const double* cvd = cvol.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t cv_counter = cv_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               l1norm += tbox::MathUtilities<TYPE>::Abs(dd[d_counter + i0])
                  * cvd[cv_counter + i0];
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
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }

         d_begin += d_offset;
         cv_begin += cv_offset;
      }
   }

   return l1norm;
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::L1Norm(
   const pdat::ArrayData<TYPE>& data,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double l1norm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();

      size_t d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      const unsigned int ddepth = data.getDepth();

      const TYPE* dd = data.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               l1norm += tbox::MathUtilities<TYPE>::Abs(dd[d_counter + i0]);
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
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
               }
            }
         }

         d_begin += d_offset;
      }
   }

   return l1norm;
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::L2NormWithControlVolume(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, cvol, box);

   return sqrt(ArrayDataNormOpsReal<TYPE>::
      dotWithControlVolume(data, data, cvol, box));
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::L2Norm(
   const pdat::ArrayData<TYPE>& data,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);

   return sqrt(ArrayDataNormOpsReal<TYPE>::dot(data, data, box));
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::weightedL2NormWithControlVolume(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<TYPE>& weight,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(data, weight, cvol, box);
   TBOX_ASSERT(data.getDepth() == weight.getDepth());

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double wl2norm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box w_box = weight.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * w_box * cv_box;

   if (!ibox.empty()) {
      const unsigned int ddepth = data.getDepth();
      const unsigned int cvdepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int w_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         w_w[i] = w_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();
      const size_t w_offset = weight.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t w_begin = w_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      const TYPE* dd = data.getPointer();
      const TYPE* wd = weight.getPointer();
      const double* cvd = cvol.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t w_counter = w_begin;
         size_t cv_counter = cv_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int w_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            w_b[nd] = static_cast<int>(w_counter);
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               TYPE val = dd[d_counter + i0] * wd[w_counter + i0];
               wl2norm += val * val * cvd[cv_counter + i0];
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
               int w_step = 1;
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  w_step *= w_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               w_counter = w_b[dim_jump - 1] + w_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  w_b[m] = static_cast<int>(w_counter);
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }

         d_begin += d_offset;
         w_begin += w_offset;
         cv_begin += cv_offset;
      }
   }

   return sqrt(wl2norm);
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::weightedL2Norm(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<TYPE>& weight,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, weight, box);
   TBOX_ASSERT(data.getDepth() == weight.getDepth());

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double wl2norm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box w_box = weight.getBox();
   const hier::Box ibox = box * d_box * w_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int w_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         w_w[i] = w_box.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();
      const size_t w_offset = weight.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t w_begin = w_box.offset(ibox.lower());

      const TYPE* dd = data.getPointer();
      const TYPE* wd = weight.getPointer();

      const unsigned int ddepth = data.getDepth();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t w_counter = w_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int w_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            w_b[nd] = static_cast<int>(w_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {
            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               TYPE val = dd[d_counter + i0] * wd[w_counter + i0];
               wl2norm += val * val;
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
               int w_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  w_step *= w_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               w_counter = w_b[dim_jump - 1] + w_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  w_b[m] = static_cast<int>(w_counter);
               }
            }
         }

         d_begin += d_offset;
         w_begin += w_offset;
      }
   }

   return sqrt(wl2norm);
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::maxNormWithControlVolume(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, cvol, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double maxnorm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d_box * cv_box;

   if (!ibox.empty()) {
      const unsigned int ddepth = data.getDepth();
      const unsigned int cvdepth = cvol.getDepth();

      TBOX_ASSERT((ddepth == cvdepth) || (cvdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int cv_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         cv_w[i] = cv_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();
      const size_t cv_offset = ((cvdepth == 1) ? 0 : cvol.getOffset());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());
      size_t cv_begin = cv_box.offset(ibox.lower());

      const TYPE* dd = data.getPointer();
      const double* cvd = cvol.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;
         size_t cv_counter = cv_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               if (cvd[cv_counter + i0] > 0.0) {
                  maxnorm =
                     tbox::MathUtilities<double>::Max(
                        maxnorm,
                        tbox::MathUtilities<TYPE>::Abs(dd[d_counter + i0]));
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
               int cv_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  cv_step *= cv_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               cv_counter = cv_b[dim_jump - 1] + cv_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }

         d_begin += d_offset;
         cv_begin += cv_offset;
      }
   }

   return maxnorm;
}

template<class TYPE>
double
ArrayDataNormOpsReal<TYPE>::maxNorm(
   const pdat::ArrayData<TYPE>& data,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   double maxnorm = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box ibox = box * d_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t d_offset = data.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t d_begin = d_box.offset(ibox.lower());

      const TYPE* dd = data.getPointer();

      const unsigned int ddepth = data.getDepth();
      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               maxnorm = tbox::MathUtilities<double>::Max(
                     maxnorm,
                     tbox::MathUtilities<TYPE>::Abs(dd[d_counter + i0]));
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
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = static_cast<int>(d_counter);
               }
            }
         }

         d_begin += d_offset;
      }
   }

   return maxnorm;
}

template<class TYPE>
TYPE
ArrayDataNormOpsReal<TYPE>::dotWithControlVolume(
   const pdat::ArrayData<TYPE>& data1,
   const pdat::ArrayData<TYPE>& data2,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(data1, data2, cvol, box);
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   tbox::Dimension::dir_t dimVal = data1.getDim().getValue();

   TYPE dprod = 0.0;

   const hier::Box d1_box = data1.getBox();
   const hier::Box d2_box = data2.getBox();
   const hier::Box cv_box = cvol.getBox();
   const hier::Box ibox = box * d1_box * d2_box * cv_box;

   if (!ibox.empty()) {
      const unsigned int d1depth = data1.getDepth();
      const unsigned int cvdepth = cvol.getDepth();

      TBOX_ASSERT(d1depth == data2.getDepth());
      TBOX_ASSERT((d1depth == cvdepth) || (cvdepth == 1));

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

      for (unsigned int d = 0; d < d1depth; ++d) {

         size_t d1_counter = d1_begin;
         size_t d2_counter = d2_begin;
         size_t cv_counter = cv_begin;

         int d1_b[SAMRAI::MAX_DIM_VAL];
         int d2_b[SAMRAI::MAX_DIM_VAL];
         int cv_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d1_b[nd] = static_cast<int>(d1_counter);
            d2_b[nd] = static_cast<int>(d2_counter);
            cv_b[nd] = static_cast<int>(cv_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dprod += static_cast<TYPE>(dd1[d1_counter + i0] * dd2[d2_counter + i0]
                                          * cvd[cv_counter + i0]);
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
                  d1_b[m] = static_cast<int>(d1_counter);
                  d2_b[m] = static_cast<int>(d2_counter);
                  cv_b[m] = static_cast<int>(cv_counter);
               }
            }
         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;
         cv_begin += cv_offset;

      }

   }

   return dprod;
}

template<class TYPE>
TYPE
ArrayDataNormOpsReal<TYPE>::dot(
   const pdat::ArrayData<TYPE>& data1,
   const pdat::ArrayData<TYPE>& data2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data1, data2, box);
   TBOX_ASSERT(data1.getDepth() == data2.getDepth());

   tbox::Dimension::dir_t dimVal = data1.getDim().getValue();

   TYPE dprod = 0.0;

   const hier::Box d1_box = data1.getBox();
   const hier::Box d2_box = data2.getBox();
   const hier::Box ibox = box * d1_box * d2_box;

   if (!ibox.empty()) {
      const unsigned int d1depth = data1.getDepth();

      TBOX_ASSERT(d1depth == data2.getDepth());

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

      for (unsigned int d = 0; d < d1depth; ++d) {

         size_t d1_counter = d1_begin;
         size_t d2_counter = d2_begin;

         int d1_b[SAMRAI::MAX_DIM_VAL];
         int d2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d1_b[nd] = static_cast<int>(d1_counter);
            d2_b[nd] = static_cast<int>(d2_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dprod += dd1[d1_counter + i0] * dd2[d2_counter + i0];
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
                  d1_b[m] = static_cast<int>(d1_counter);
                  d2_b[m] = static_cast<int>(d2_counter);
               }
            }
         }

         d1_begin += d1_offset;
         d2_begin += d2_offset;

      }

   }

   return dprod;
}

template<class TYPE>
TYPE
ArrayDataNormOpsReal<TYPE>::integral(
   const pdat::ArrayData<TYPE>& data,
   const pdat::ArrayData<double>& vol,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(data, vol, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();

   TYPE integral = 0.0;

   const hier::Box d_box = data.getBox();
   const hier::Box v_box = vol.getBox();
   const hier::Box ibox = box * d_box * v_box;

   if (!ibox.empty()) {
      const int ddepth = data.getDepth();
      const int vdepth = vol.getDepth();

      TBOX_ASSERT((ddepth == vdepth) || (vdepth == 1));

      int box_w[SAMRAI::MAX_DIM_VAL];
      int d_w[SAMRAI::MAX_DIM_VAL];
      int v_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         d_w[i] = d_box.numberCells(i);
         v_w[i] = v_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int d_offset = static_cast<int>(data.getOffset());
      const int v_offset =
         ((vdepth == 1) ? 0 : static_cast<int>(vol.getOffset()));

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      int d_begin = static_cast<int>(d_box.offset(ibox.lower()));
      int v_begin = static_cast<int>(v_box.offset(ibox.lower()));

      const TYPE* dd = data.getPointer();
      const double* vd = vol.getPointer();

      for (int d = 0; d < ddepth; ++d) {

         int d_counter = d_begin;
         int v_counter = v_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         int v_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = d_counter;
            v_b[nd] = v_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               integral += dd[d_counter + i0] * vd[v_counter + i0];
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
               int v_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  d_step *= d_w[k];
                  v_step *= v_w[k];
               }
               d_counter = d_b[dim_jump - 1] + d_step;
               v_counter = v_b[dim_jump - 1] + v_step;

               for (int m = 0; m < dim_jump; ++m) {
                  d_b[m] = d_counter;
                  v_b[m] = v_counter;
               }
            }
         }

         d_begin += d_offset;
         v_begin += v_offset;
      }
   }

   return integral;
}

}
}
#endif
