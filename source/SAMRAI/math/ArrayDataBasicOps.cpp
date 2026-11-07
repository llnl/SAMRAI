/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated opertions for array data.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataBasicOps_C
#define included_math_ArrayDataBasicOps_C

#include "SAMRAI/math/ArrayDataBasicOps.h"

#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
ArrayDataBasicOps<TYPE>::ArrayDataBasicOps()
{
}

template<class TYPE>
ArrayDataBasicOps<TYPE>::~ArrayDataBasicOps()
{
}

/*
 *************************************************************************
 *
 * General templated operations for array data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::scale(
   pdat::ArrayData<TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<TYPE>& src,
   const hier::Box& box) const
{
// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.fillAll(alpha, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      dst.copy(src, box);
   } else {
      const unsigned int ddepth = dst.getDepth();

      TBOX_ASSERT(ddepth == src.getDepth());

      const hier::Box dst_box = dst.getBox();
      const hier::Box src_box = src.getBox();
      const hier::Box ibox = box * dst_box * src_box;

      if (!ibox.empty()) {

         int box_w[SAMRAI::MAX_DIM_VAL];
         int dst_w[SAMRAI::MAX_DIM_VAL];
         int src_w[SAMRAI::MAX_DIM_VAL];
         int dim_counter[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src_w[i] = src_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const size_t dst_offset = dst.getOffset();
         const size_t src_offset = src.getOffset();

         const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

         size_t dst_begin = dst_box.offset(ibox.lower());
         size_t src_begin = src_box.offset(ibox.lower());

         TYPE* dd = dst.getPointer();
         const TYPE* sd = src.getPointer();

         for (unsigned int d = 0; d < ddepth; ++d) {

            size_t dst_counter = dst_begin;
            size_t src_counter = src_begin;

            int dst_b[SAMRAI::MAX_DIM_VAL];
            int src_b[SAMRAI::MAX_DIM_VAL];
            for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
               dst_b[nd] = static_cast<int>(dst_counter);
               src_b[nd] = static_cast<int>(src_counter);
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; ++nb) {

               for (int i0 = 0; i0 < box_w[0]; ++i0) {
                  dd[dst_counter + i0] = alpha * sd[src_counter + i0];
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
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::addScalar(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.copy(src, box);
   } else {
      const unsigned int ddepth = dst.getDepth();

      TBOX_ASSERT(ddepth == src.getDepth());

      const hier::Box dst_box = dst.getBox();
      const hier::Box src_box = src.getBox();
      const hier::Box ibox = box * dst_box * src_box;

      if (!ibox.empty()) {

         int box_w[SAMRAI::MAX_DIM_VAL];
         int dst_w[SAMRAI::MAX_DIM_VAL];
         int src_w[SAMRAI::MAX_DIM_VAL];
         int dim_counter[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src_w[i] = src_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const size_t dst_offset = dst.getOffset();
         const size_t src_offset = src.getOffset();

         const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

         size_t dst_begin = dst_box.offset(ibox.lower());
         size_t src_begin = src_box.offset(ibox.lower());

         TYPE* dd = dst.getPointer();
         const TYPE* sd = src.getPointer();

         for (unsigned int d = 0; d < ddepth; ++d) {

            size_t dst_counter = dst_begin;
            size_t src_counter = src_begin;

            int dst_b[SAMRAI::MAX_DIM_VAL];
            int src_b[SAMRAI::MAX_DIM_VAL];
            for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
               dst_b[nd] = static_cast<int>(dst_counter);
               src_b[nd] = static_cast<int>(src_counter);
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; ++nb) {

               for (int i0 = 0; i0 < box_w[0]; ++i0) {
                  dd[dst_counter + i0] = alpha + sd[src_counter + i0];
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
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::add(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

   const hier::Box dst_box = dst.getBox();
   const hier::Box src1_box = src1.getBox();
   const hier::Box src2_box = src2.getBox();
   const hier::Box ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src1_w[SAMRAI::MAX_DIM_VAL];
      int src2_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src1_offset = src1.getOffset();
      const size_t src2_offset = src2.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src1_begin = src1_box.offset(ibox.lower());
      size_t src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src1_counter = src1_begin;
         size_t src2_counter = src2_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src1_b[SAMRAI::MAX_DIM_VAL];
         int src2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src1_b[nd] = static_cast<int>(src1_counter);
            src2_b[nd] = static_cast<int>(src2_counter);
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] = s1d[src1_counter + i0]
                  + s2d[src2_counter + i0];
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
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src1_counter = src1_b[dim_jump - 1] + src1_step;
               src2_counter = src2_b[dim_jump - 1] + src2_step;

               for (int m = 0; m < dim_jump; ++m) {
                  dst_b[m] = static_cast<int>(dst_counter);
                  src1_b[m] = static_cast<int>(src1_counter);
                  src2_b[m] = static_cast<int>(src2_counter);
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::subtract(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

   const hier::Box dst_box = dst.getBox();
   const hier::Box src1_box = src1.getBox();
   const hier::Box src2_box = src2.getBox();
   const hier::Box ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src1_w[SAMRAI::MAX_DIM_VAL];
      int src2_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src1_offset = src1.getOffset();
      const size_t src2_offset = src2.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src1_begin = src1_box.offset(ibox.lower());
      size_t src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src1_counter = src1_begin;
         size_t src2_counter = src2_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src1_b[SAMRAI::MAX_DIM_VAL];
         int src2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src1_b[nd] = static_cast<int>(src1_counter);
            src2_b[nd] = static_cast<int>(src2_counter);
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] = s1d[src1_counter + i0]
                  - s2d[src2_counter + i0];
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
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src1_counter = src1_b[dim_jump - 1] + src1_step;
               src2_counter = src2_b[dim_jump - 1] + src2_step;

               for (int m = 0; m < dim_jump; ++m) {
                  dst_b[m] = static_cast<int>(dst_counter);
                  src1_b[m] = static_cast<int>(src1_counter);
                  src2_b[m] = static_cast<int>(src2_counter);
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::multiply(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();
   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

   const hier::Box dst_box = dst.getBox();
   const hier::Box src1_box = src1.getBox();
   const hier::Box src2_box = src2.getBox();
   const hier::Box ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src1_w[SAMRAI::MAX_DIM_VAL];
      int src2_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src1_offset = src1.getOffset();
      const size_t src2_offset = src2.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src1_begin = src1_box.offset(ibox.lower());
      size_t src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src1_counter = src1_begin;
         size_t src2_counter = src2_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src1_b[SAMRAI::MAX_DIM_VAL];
         int src2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src1_b[nd] = static_cast<int>(src1_counter);
            src2_b[nd] = static_cast<int>(src2_counter);
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] = s1d[src1_counter + i0]
                  * s2d[src2_counter + i0];
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
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src1_counter = src1_b[dim_jump - 1] + src1_step;
               src2_counter = src2_b[dim_jump - 1] + src2_step;

               for (int m = 0; m < dim_jump; ++m) {
                  dst_b[m] = static_cast<int>(dst_counter);
                  src1_b[m] = static_cast<int>(src1_counter);
                  src2_b[m] = static_cast<int>(src2_counter);
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::divide(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();
   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

   const hier::Box dst_box = dst.getBox();
   const hier::Box src1_box = src1.getBox();
   const hier::Box src2_box = src2.getBox();
   const hier::Box ibox = box * dst_box * src1_box * src2_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src1_w[SAMRAI::MAX_DIM_VAL];
      int src2_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src1_w[i] = src1_box.numberCells(i);
         src2_w[i] = src2_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src1_offset = src1.getOffset();
      const size_t src2_offset = src2.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src1_begin = src1_box.offset(ibox.lower());
      size_t src2_begin = src2_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* s1d = src1.getPointer();
      const TYPE* s2d = src2.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src1_counter = src1_begin;
         size_t src2_counter = src2_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src1_b[SAMRAI::MAX_DIM_VAL];
         int src2_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src1_b[nd] = static_cast<int>(src1_counter);
            src2_b[nd] = static_cast<int>(src2_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] = s1d[src1_counter + i0]
                  / s2d[src2_counter + i0];
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
               int src1_step = 1;
               int src2_step = 1;
               for (int k = 0; k < dim_jump; ++k) {
                  dst_step *= dst_w[k];
                  src1_step *= src1_w[k];
                  src2_step *= src2_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src1_counter = src1_b[dim_jump - 1] + src1_step;
               src2_counter = src2_b[dim_jump - 1] + src2_step;

               for (int m = 0; m < dim_jump; ++m) {
                  dst_b[m] = static_cast<int>(dst_counter);
                  src1_b[m] = static_cast<int>(src1_counter);
                  src2_b[m] = static_cast<int>(src2_counter);
               }

            }
         }

         dst_begin += dst_offset;
         src1_begin += src1_offset;
         src2_begin += src2_offset;

      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::reciprocal(
   pdat::ArrayData<TYPE>& dst,
   const pdat::ArrayData<TYPE>& src,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(dst, src, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();
   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src.getDepth());

   const hier::Box dst_box = dst.getBox();
   const hier::Box src_box = src.getBox();
   const hier::Box ibox = box * dst_box * src_box;

   if (!ibox.empty()) {

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int src_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src_w[i] = src_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const size_t dst_offset = dst.getOffset();
      const size_t src_offset = src.getOffset();

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      size_t dst_begin = dst_box.offset(ibox.lower());
      size_t src_begin = src_box.offset(ibox.lower());

      TYPE* dd = dst.getPointer();
      const TYPE* sd = src.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t dst_counter = dst_begin;
         size_t src_counter = src_begin;

         int dst_b[SAMRAI::MAX_DIM_VAL];
         int src_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            dst_b[nd] = static_cast<int>(dst_counter);
            src_b[nd] = static_cast<int>(src_counter);
         }

         /*
          * Loop over each contiguous block of data.
          */
         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[dst_counter + i0] =
                  tbox::MathUtilities<TYPE>::getOne() / sd[src_counter + i0];
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
void
ArrayDataBasicOps<TYPE>::linearSum(
   pdat::ArrayData<TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<TYPE>& src1,
   const TYPE& beta,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
// Ignore Intel warning about floating point comparisons
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();
   const unsigned int ddepth = dst.getDepth();

   TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      if (beta == tbox::MathUtilities<TYPE>::getZero()) {
         dst.fillAll(tbox::MathUtilities<TYPE>::getZero(), box);
      } else {
         scale(dst, beta, src2, box);
      }
   } else if (beta == tbox::MathUtilities<TYPE>::getZero()) {
      scale(dst, alpha, src1, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      axpy(dst, beta, src2, src1, box);
   } else if (beta == tbox::MathUtilities<TYPE>::getOne()) {
      axpy(dst, alpha, src1, src2, box);
   } else if (alpha == -tbox::MathUtilities<TYPE>::getOne()) {
      axmy(dst, beta, src2, src1, box);
   } else if (beta == -tbox::MathUtilities<TYPE>::getOne()) {
      axmy(dst, alpha, src1, src2, box);
   } else {

      const hier::Box dst_box = dst.getBox();
      const hier::Box src1_box = src1.getBox();
      const hier::Box src2_box = src2.getBox();
      const hier::Box ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[SAMRAI::MAX_DIM_VAL];
         int dst_w[SAMRAI::MAX_DIM_VAL];
         int src1_w[SAMRAI::MAX_DIM_VAL];
         int src2_w[SAMRAI::MAX_DIM_VAL];
         int dim_counter[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const size_t dst_offset = dst.getOffset();
         const size_t src1_offset = src1.getOffset();
         const size_t src2_offset = src2.getOffset();

         const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

         size_t dst_begin = dst_box.offset(ibox.lower());
         size_t src1_begin = src1_box.offset(ibox.lower());
         size_t src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (unsigned int d = 0; d < ddepth; ++d) {

            size_t dst_counter = dst_begin;
            size_t src1_counter = src1_begin;
            size_t src2_counter = src2_begin;

            int dst_b[SAMRAI::MAX_DIM_VAL];
            int src1_b[SAMRAI::MAX_DIM_VAL];
            int src2_b[SAMRAI::MAX_DIM_VAL];
            for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
               dst_b[nd] = static_cast<int>(dst_counter);
               src1_b[nd] = static_cast<int>(src1_counter);
               src2_b[nd] = static_cast<int>(src2_counter);
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; ++nb) {

               for (int i0 = 0; i0 < box_w[0]; ++i0) {
                  dd[dst_counter + i0] = alpha * s1d[src1_counter + i0]
                     + beta * s2d[src2_counter + i0];
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
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; ++k) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump - 1] + dst_step;
                  src1_counter = src1_b[dim_jump - 1] + src1_step;
                  src2_counter = src2_b[dim_jump - 1] + src2_step;

                  for (int m = 0; m < dim_jump; ++m) {
                     dst_b[m] = static_cast<int>(dst_counter);
                     src1_b[m] = static_cast<int>(src1_counter);
                     src2_b[m] = static_cast<int>(src2_counter);
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::axpy(
   pdat::ArrayData<TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      dst.copy(src2, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      add(dst, src1, src2, box);
   } else if (alpha == -tbox::MathUtilities<TYPE>::getOne()) {
      subtract(dst, src2, src1, box);
   } else {
      const unsigned int ddepth = dst.getDepth();

      TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

      const hier::Box dst_box = dst.getBox();
      const hier::Box src1_box = src1.getBox();
      const hier::Box src2_box = src2.getBox();
      const hier::Box ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[SAMRAI::MAX_DIM_VAL];
         int dst_w[SAMRAI::MAX_DIM_VAL];
         int src1_w[SAMRAI::MAX_DIM_VAL];
         int src2_w[SAMRAI::MAX_DIM_VAL];
         int dim_counter[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const size_t dst_offset = dst.getOffset();
         const size_t src1_offset = src1.getOffset();
         const size_t src2_offset = src2.getOffset();

         const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

         size_t dst_begin = dst_box.offset(ibox.lower());
         size_t src1_begin = src1_box.offset(ibox.lower());
         size_t src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (unsigned int d = 0; d < ddepth; ++d) {

            size_t dst_counter = dst_begin;
            size_t src1_counter = src1_begin;
            size_t src2_counter = src2_begin;

            int dst_b[SAMRAI::MAX_DIM_VAL];
            int src1_b[SAMRAI::MAX_DIM_VAL];
            int src2_b[SAMRAI::MAX_DIM_VAL];
            for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
               dst_b[nd] = static_cast<int>(dst_counter);
               src1_b[nd] = static_cast<int>(src1_counter);
               src2_b[nd] = static_cast<int>(src2_counter);
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; ++nb) {

               for (int i0 = 0; i0 < box_w[0]; ++i0) {
                  dd[dst_counter + i0] = alpha * s1d[src1_counter + i0]
                     + s2d[src2_counter + i0];
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
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; ++k) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump - 1] + dst_step;
                  src1_counter = src1_b[dim_jump - 1] + src1_step;
                  src2_counter = src2_b[dim_jump - 1] + src2_step;

                  for (int m = 0; m < dim_jump; ++m) {
                     dst_b[m] = static_cast<int>(dst_counter);
                     src1_b[m] = static_cast<int>(src1_counter);
                     src2_b[m] = static_cast<int>(src2_counter);
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::axmy(
   pdat::ArrayData<TYPE>& dst,
   const TYPE& alpha,
   const pdat::ArrayData<TYPE>& src1,
   const pdat::ArrayData<TYPE>& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src1, src2, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();

   if (alpha == tbox::MathUtilities<TYPE>::getZero()) {
      scale(dst, -tbox::MathUtilities<TYPE>::getOne(), src2, box);
   } else if (alpha == tbox::MathUtilities<TYPE>::getOne()) {
      subtract(dst, src1, src2, box);
   } else {
      const unsigned int ddepth = dst.getDepth();

      TBOX_ASSERT(ddepth == src1.getDepth() && ddepth == src2.getDepth());

      const hier::Box dst_box = dst.getBox();
      const hier::Box src1_box = src1.getBox();
      const hier::Box src2_box = src2.getBox();
      const hier::Box ibox = box * dst_box * src1_box * src2_box;

      if (!ibox.empty()) {

         int box_w[SAMRAI::MAX_DIM_VAL];
         int dst_w[SAMRAI::MAX_DIM_VAL];
         int src1_w[SAMRAI::MAX_DIM_VAL];
         int src2_w[SAMRAI::MAX_DIM_VAL];
         int dim_counter[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t i = 0; i < dimVal; ++i) {
            box_w[i] = ibox.numberCells(i);
            dst_w[i] = dst_box.numberCells(i);
            src1_w[i] = src1_box.numberCells(i);
            src2_w[i] = src2_box.numberCells(i);
            dim_counter[i] = 0;
         }

         const size_t dst_offset = dst.getOffset();
         const size_t src1_offset = src1.getOffset();
         const size_t src2_offset = src2.getOffset();

         const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

         size_t dst_begin = dst_box.offset(ibox.lower());
         size_t src1_begin = src1_box.offset(ibox.lower());
         size_t src2_begin = src2_box.offset(ibox.lower());

         TYPE* dd = dst.getPointer();
         const TYPE* s1d = src1.getPointer();
         const TYPE* s2d = src2.getPointer();

         for (unsigned int d = 0; d < ddepth; ++d) {

            size_t dst_counter = dst_begin;
            size_t src1_counter = src1_begin;
            size_t src2_counter = src2_begin;

            int dst_b[SAMRAI::MAX_DIM_VAL];
            int src1_b[SAMRAI::MAX_DIM_VAL];
            int src2_b[SAMRAI::MAX_DIM_VAL];
            for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
               dst_b[nd] = static_cast<int>(dst_counter);
               src1_b[nd] = static_cast<int>(src1_counter);
               src2_b[nd] = static_cast<int>(src2_counter);
            }

            /*
             * Loop over each contiguous block of data.
             */
            for (int nb = 0; nb < num_d0_blocks; ++nb) {

               for (int i0 = 0; i0 < box_w[0]; ++i0) {
                  dd[dst_counter + i0] = alpha * s1d[src1_counter + i0]
                     - s2d[src2_counter + i0];
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
                  int src1_step = 1;
                  int src2_step = 1;
                  for (int k = 0; k < dim_jump; ++k) {
                     dst_step *= dst_w[k];
                     src1_step *= src1_w[k];
                     src2_step *= src2_w[k];
                  }
                  dst_counter = dst_b[dim_jump - 1] + dst_step;
                  src1_counter = src1_b[dim_jump - 1] + src1_step;
                  src2_counter = src2_b[dim_jump - 1] + src2_step;

                  for (int m = 0; m < dim_jump; ++m) {
                     dst_b[m] = static_cast<int>(dst_counter);
                     src1_b[m] = static_cast<int>(src1_counter);
                     src2_b[m] = static_cast<int>(src2_counter);
                  }

               }
            }

            dst_begin += dst_offset;
            src1_begin += src1_offset;
            src2_begin += src2_offset;

         }
      }
   }
}

template<class TYPE>
TYPE
ArrayDataBasicOps<TYPE>::min(
   const pdat::ArrayData<TYPE>& data,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

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
               minval =
                  tbox::MathUtilities<TYPE>::Min(minval, dd[d_counter + i0]);
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

   return minval;
}

template<class TYPE>
TYPE
ArrayDataBasicOps<TYPE>::max(
   const pdat::ArrayData<TYPE>& data,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);

   tbox::Dimension::dir_t dimVal = data.getDim().getValue();
   TYPE maxval = -(tbox::MathUtilities<TYPE>::getMax());

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
               maxval =
                  tbox::MathUtilities<TYPE>::Max(maxval, dd[d_counter + i0]);
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

   return maxval;
}

template<class TYPE>
void
ArrayDataBasicOps<TYPE>::setRandomValues(
   pdat::ArrayData<TYPE>& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(dst, box);

   tbox::Dimension::dir_t dimVal = dst.getDim().getValue();
   const hier::Box d_box = dst.getBox();
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

      const size_t d_offset = dst.getOffset();

      size_t d_begin = d_box.offset(ibox.lower());

      const int num_d0_blocks = static_cast<int>(ibox.size() / box_w[0]);

      const unsigned int ddepth = dst.getDepth();

      TYPE* dd = dst.getPointer();

      for (unsigned int d = 0; d < ddepth; ++d) {

         size_t d_counter = d_begin;

         int d_b[SAMRAI::MAX_DIM_VAL];
         for (tbox::Dimension::dir_t nd = 0; nd < dimVal; ++nd) {
            d_b[nd] = static_cast<int>(d_counter);
         }

         for (int nb = 0; nb < num_d0_blocks; ++nb) {

            for (int i0 = 0; i0 < box_w[0]; ++i0) {
               dd[d_counter + i0] = tbox::MathUtilities<TYPE>::Rand(low, width);
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
}

}
}
#endif
