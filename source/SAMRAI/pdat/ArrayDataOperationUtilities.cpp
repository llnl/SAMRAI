/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated array data looping operations supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataOperationUtilities_C
#define included_pdat_ArrayDataOperationUtilities_C

#include "SAMRAI/pdat/ArrayDataOperationUtilities.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/ForAll.h"
#include "SAMRAI/pdat/SumOperation.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"
#include "SAMRAI/tbox/Utilities.h"

#include <array>

namespace SAMRAI
{
namespace pdat
{

/*
 * Note on the usage of RAJA in this implementation:
 *
 * To support the allocation of ArrayData on either the host or device when
 * running on GPU architectures, each RAJA kernel is implemented twice in
 * an if/else code block.  The looping structure host_parallel_for_all()
 * is used when ArrayData is allocated on the host and guarantees that
 * the loop will execute on the host.  When the GPU devices is available and
 * ArrayData is allocated on the device, parallel_for_all() is used to
 * launch the kernels on the device.
 */

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination array data objects and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

template <class TYPE, class OP>
void ArrayDataOperationUtilities<TYPE, OP>::doArrayDataOperationOnBox(
    ArrayData<TYPE>& dst,
    const ArrayData<TYPE>& src,
    const hier::Box& opbox,
    const hier::IntVector& src_shift,
    unsigned int dst_start_depth,
    unsigned int src_start_depth,
    unsigned int num_depth,
    const OP& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src, opbox, src_shift);
   TBOX_ASSERT((dst_start_depth + num_depth <= dst.getDepth()));
   TBOX_ASSERT((src_start_depth + num_depth <= src.getDepth()));

   const tbox::Dimension& dim(dst.getDim());

   TYPE* const dst_ptr = dst.getPointer();
   const TYPE* const src_ptr = src.getPointer();

#if !defined(HAVE_RAJA)
   const hier::Box& dst_box(dst.getBox());
   const hier::Box& src_box(src.getBox());

   std::array<int,SAMRAI::MAX_DIM_VAL> box_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dst_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> src_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dim_counter = {};
   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dst_w[i] = dst_box.numberCells(i);
      src_w[i] = src_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dst_offset = dst.getOffset();
   const size_t src_offset = src.getOffset();

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dst_begin, src_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dst_begin = dst_box.offset(opbox.lower()) + dst_start_depth * dst_offset;
   size_t src_begin = src_box.offset(opbox.lower() - src_shift) + src_start_depth * src_offset;

#else
   NULL_USE(src_ptr);
   NULL_USE(dst_ptr);
   NULL_USE(src_shift);
   NULL_USE(dst_start_depth);
   NULL_USE(src_start_depth);

   bool src_on_host = src.dataOnHost();
   bool dst_on_host = dst.dataOnHost();
   TBOX_ASSERT(src_on_host == dst_on_host);
   bool on_host = (src_on_host && dst_on_host);
#endif

#if defined(HAVE_RAJA)
   bool use_fuser = dst.useFuser();
   tbox::KernelFuser* fuser = use_fuser ?
      tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;
#endif

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < num_depth; ++d) {

#if defined(HAVE_RAJA)

      switch (dim.getValue()) {
         case 1: {
            auto dest = get_view<1>(dst, d);
            auto source = get_const_view<1>(src, d);
            const int shift_i = src_shift[0];
            auto s2 = source.shift({{shift_i}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i) {
                  op(dest(i), s2(i));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i) {
                  op(dest(i), s2(i));
               });
            }
         } break;

         case 2: {
            auto dest = get_view<2>(dst, d);
            auto source = get_const_view<2>(src, d);
            const int shift_i = src_shift[0];
            const int shift_j = src_shift[1];
            auto s2 = source.shift({{shift_i, shift_j}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j) {
                  op(dest(i, j), s2(i, j));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i, int j) {
                  op(dest(i, j), s2(i, j));
               });
            }
         } break;

         case 3: {
            auto dest = get_view<3>(dst, d);
            auto source = get_const_view<3>(src, d);
            const int shift_i = src_shift[0];
            const int shift_j = src_shift[1];
            const int shift_k = src_shift[2];
            auto s2 = source.shift({{shift_i, shift_j, shift_k}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j, int k) {
                  op(dest(i, j, k), s2(i, j, k));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  op(dest(i, j, k), s2(i, j, k));
               });
            }
         } break;

         default:
            tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
            TBOX_ERROR("Aborting in " __FILE__);
            break;
      }
#else

      size_t dst_counter = dst_begin;
      size_t src_counter = src_begin;

      std::array<size_t,SAMRAI::MAX_DIM_VAL> dst_b = {};
      std::array<size_t,SAMRAI::MAX_DIM_VAL> src_b = {};
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dst_b[nd] = dst_counter;
         src_b[nd] = src_counter;
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is copied, calculate the
          * beginning array index for the next block.
          */

         for (tbox::Dimension::dir_t j = 1; j < dim.getValue(); ++j) {
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
               dst_b[m] = dst_counter;
               src_b[m] = src_counter;
            }

         }  // if dim_jump > 0

      }  // nb loop over contiguous data blocks

      /*
       * After copy is complete on a full box for one depth index,
       * advance by the offset values.
       */

      dst_begin += dst_offset;
      src_begin += src_offset;
#endif

   }  // d loop over depth indices
}  // end doArrayDataOperationOnBox

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination data pointers and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

template <class TYPE, class OP>
void ArrayDataOperationUtilities<TYPE, OP>::doArrayDataBufferOperationOnBox(
    const ArrayData<TYPE>& arraydata,
    const TYPE* buffer,
    const hier::Box& opbox,
    bool src_is_buffer,
    const OP& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(arraydata, opbox);
   TBOX_ASSERT(buffer != 0);
   TBOX_ASSERT(opbox.isSpatiallyEqual((opbox * arraydata.getBox())));

   const tbox::Dimension& dim(arraydata.getDim());

   TYPE* const dst_ptr =
       (src_is_buffer ? const_cast<TYPE*>(arraydata.getPointer())
                      : const_cast<TYPE*>(buffer));
   const TYPE* const src_ptr =
       (src_is_buffer ? buffer : arraydata.getPointer());

   const unsigned int array_d_depth = arraydata.getDepth();

#if !defined(HAVE_RAJA)
   const hier::Box& array_d_box(arraydata.getBox());

   std::array<int,SAMRAI::MAX_DIM_VAL> box_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dat_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dim_counter = {};
   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dat_w[i] = array_d_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dat_offset = arraydata.getOffset();
   const size_t buf_offset = box_w[0];

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dat_begin, buf_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dat_begin = array_d_box.offset(opbox.lower());
   size_t buf_begin = 0;
#else
   bool on_host = arraydata.dataOnHost();
#endif

#if defined(HAVE_RAJA)
   bool use_fuser = arraydata.useFuser();
   tbox::KernelFuser* fuser = use_fuser ?
      tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;
#endif

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < array_d_depth; ++d) {

#if defined(HAVE_RAJA)
      const hier::Box& dst_box = src_is_buffer ? arraydata.getBox() : opbox;
      const hier::Box& src_box = src_is_buffer ? opbox : arraydata.getBox();

      const int dst_offset = dst_box.size() * d;
      const int src_offset = src_box.size() * d;

      switch (dim.getValue()) {
         case 1: {
            typename pdat::ArrayData<TYPE>::template View<1> dest(dst_ptr + dst_offset, dst_box);
            typename pdat::ArrayData<TYPE>::template ConstView<1> source(src_ptr + src_offset, src_box);

            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i) {
                  op(dest(i), source(i));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i) {
                  op(dest(i), source(i));
               });
            }
         } break;

         case 2: {
            typename pdat::ArrayData<TYPE>::template View<2> dest(dst_ptr + dst_offset, dst_box);
            typename pdat::ArrayData<TYPE>::template ConstView<2> source(src_ptr + src_offset, src_box);
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j) {
                  op(dest(i, j), source(i, j));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i, int j) {
                  op(dest(i, j), source(i, j));
               });
            }
         } break;

         case 3: {
            typename pdat::ArrayData<TYPE>::template View<3> dest(dst_ptr + dst_offset, dst_box);
            typename pdat::ArrayData<TYPE>::template ConstView<3> source(src_ptr + src_offset, src_box);
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j, int k) {
                  op(dest(i, j, k), source(i, j, k));
               });
            } else {
               hier::parallel_for_all(fuser, opbox, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  op(dest(i, j, k), source(i, j, k));
               });
            }
         } break;

         default:
            tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
            TBOX_ERROR("Aborting in " __FILE__);
            break;
      }

#else  // !HAVE_RAJA
      size_t dat_counter = dat_begin;
      size_t buf_counter = buf_begin;

      size_t& dst_counter = (src_is_buffer ? dat_counter : buf_counter);
      size_t& src_counter = (src_is_buffer ? buf_counter : dat_counter);

      std::array<int,SAMRAI::MAX_DIM_VAL> dat_b = {};
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dat_b[nd] = static_cast<int>(dat_counter);
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is packed, calculate the
          * beginning array index for the next block.
          */

         for (int j = 1; j < dim.getValue(); ++j) {
            if (dim_counter[j] < box_w[j] - 1) {
               ++dim_counter[j];
               dim_jump = j;
               break;
            } else {
               dim_counter[j] = 0;
            }
         }

         if (dim_jump > 0) {

            int dat_step = 1;
            for (int k = 0; k < dim_jump; ++k) {
               dat_step *= dat_w[k];
            }
            dat_counter = dat_b[dim_jump - 1] + dat_step;

            for (int m = 0; m < dim_jump; ++m) {
               dat_b[m] = static_cast<int>(dat_counter);
            }

         }  // if dim_jump > 0

         buf_counter += buf_offset;

      }  // nb loop over contiguous data blocks

      /*
       * After packing is complete on a full box for one depth index,
       * advance by the offset value.
       */

      dat_begin += dat_offset;
      buf_begin = buf_counter;
#endif

   }  // d loop over depth indices
}  // end doArrayDataBufferOperationOnBox

// specialize for dcomplex sum
template <>
inline void ArrayDataOperationUtilities<dcomplex,SumOperation<dcomplex> >::doArrayDataOperationOnBox(
    ArrayData<dcomplex>& dst,
    const ArrayData<dcomplex>& src,
    const hier::Box& opbox,
    const hier::IntVector& src_shift,
    unsigned int dst_start_depth,
    unsigned int src_start_depth,
    unsigned int num_depth,
    const SumOperation<dcomplex>& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst, src, opbox, src_shift);
   TBOX_ASSERT((dst_start_depth + num_depth <= dst.getDepth()));
   TBOX_ASSERT((src_start_depth + num_depth <= src.getDepth()));

   const tbox::Dimension& dim(dst.getDim());

   dcomplex* const dst_ptr = dst.getPointer();
   const dcomplex* const src_ptr = src.getPointer();

#if !defined(HAVE_RAJA)
   const hier::Box& dst_box(dst.getBox());
   const hier::Box& src_box(src.getBox());

   std::array<int,SAMRAI::MAX_DIM_VAL> box_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dst_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> src_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dim_counter = {};

   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dst_w[i] = dst_box.numberCells(i);
      src_w[i] = src_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dst_offset = dst.getOffset();
   const size_t src_offset = src.getOffset();

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dst_begin, src_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dst_begin = dst_box.offset(opbox.lower()) + dst_start_depth * dst_offset;
   size_t src_begin = src_box.offset(opbox.lower() - src_shift) + src_start_depth * src_offset;

#else
   NULL_USE(src_ptr);
   NULL_USE(dst_ptr);
   NULL_USE(src_shift);
   NULL_USE(dst_start_depth);
   NULL_USE(src_start_depth);
   NULL_USE(op);

   bool src_on_host = src.dataOnHost();
   bool dst_on_host = dst.dataOnHost();
   TBOX_ASSERT(src_on_host == dst_on_host);
   bool on_host = (src_on_host && dst_on_host);
#endif

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < num_depth; ++d) {

#if defined(HAVE_RAJA)
      SumOperation<double> sumop_dbl;
      switch (dim.getValue()) {
         case 1: {
            auto dest = get_view<1>(dst, d);
            auto source = get_const_view<1>(src, d);
            const int shift_i = src_shift[0];
            auto s2 = source.shift({{shift_i}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            }
         } break;

         case 2: {
            auto dest = get_view<2>(dst, d);
            auto source = get_const_view<2>(src, d);
            const int shift_i = src_shift[0];
            const int shift_j = src_shift[1];
            auto s2 = source.shift({{shift_i, shift_j}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i,j))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i,j))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i, int j) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i,j))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i,j))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            }
         } break;

         case 3: {
            auto dest = get_view<3>(dst, d);
            auto source = get_const_view<3>(src, d);
            const int shift_i = src_shift[0];
            const int shift_j = src_shift[1];
            const int shift_k = src_shift[2];
            auto s2 = source.shift({{shift_i, shift_j, shift_k}});
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j, int k) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j,k))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j,k))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i,j,k))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i,j,k))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j,k))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j,k))[1];
                  const double &s2_real = reinterpret_cast<const double(&)[2]>(s2(i,j,k))[0];
                  const double &s2_imag = reinterpret_cast<const double(&)[2]>(s2(i,j,k))[1];
                  sumop_dbl(dest_real, s2_real);
                  sumop_dbl(dest_imag, s2_imag);
               });
            }
         } break;

         default:
            tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
            TBOX_ERROR("Aborting in " __FILE__);
            break;
      }

#else

      size_t dst_counter = dst_begin;
      size_t src_counter = src_begin;

      std::array<size_t,SAMRAI::MAX_DIM_VAL> dst_b = {};
      std::array<size_t,SAMRAI::MAX_DIM_VAL> src_b = {};
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dst_b[nd] = dst_counter;
         src_b[nd] = src_counter;
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is copied, calculate the
          * beginning array index for the next block.
          */

         for (tbox::Dimension::dir_t j = 1; j < dim.getValue(); ++j) {
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
               dst_b[m] = dst_counter;
               src_b[m] = src_counter;
            }

         }  // if dim_jump > 0

      }  // nb loop over contiguous data blocks

      /*
       * After copy is complete on a full box for one depth index,
       * advance by the offset values.
       */

      dst_begin += dst_offset;
      src_begin += src_offset;
#endif

   }  // d loop over depth indices
}  // end doArrayDataOperationOnBox

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination data pointers and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

// specialize for dcomplex sum
template <>
inline void ArrayDataOperationUtilities<dcomplex, SumOperation<dcomplex> >::doArrayDataBufferOperationOnBox(
    const ArrayData<dcomplex>& arraydata,
    const dcomplex* buffer,
    const hier::Box& opbox,
    bool src_is_buffer,
    const SumOperation<dcomplex>& op)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(arraydata, opbox);
   TBOX_ASSERT(buffer != 0);
   TBOX_ASSERT(opbox.isSpatiallyEqual((opbox * arraydata.getBox())));

   const tbox::Dimension& dim(arraydata.getDim());


   dcomplex* const dst_ptr =
       (src_is_buffer ? const_cast<dcomplex*>(arraydata.getPointer())
                      : const_cast<dcomplex*>(buffer));
   const dcomplex* const src_ptr =
       (src_is_buffer ? buffer : arraydata.getPointer());

   const unsigned int array_d_depth = arraydata.getDepth();

#if !defined(HAVE_RAJA)
   const hier::Box& array_d_box(arraydata.getBox());

   std::array<int,SAMRAI::MAX_DIM_VAL> box_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dat_w = {};
   std::array<int,SAMRAI::MAX_DIM_VAL> dim_counter = {};
   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      box_w[i] = opbox.numberCells(i);
      dat_w[i] = array_d_box.numberCells(i);
      dim_counter[i] = 0;
   }

   const size_t dat_offset = arraydata.getOffset();
   const size_t buf_offset = box_w[0];

   /*
    * Data on the opbox can be decomposed into a set of
    * contiguous array sections representing data in a straight line
    * in the 0 coordinate direction.
    *
    * num_d0_blocks is the number of such array sections.
    * dat_begin, buf_begin are the array indices for the first
    * data items in each array section to be copied.
    */

   const int num_d0_blocks = static_cast<int>(opbox.size() / box_w[0]);

   size_t dat_begin = array_d_box.offset(opbox.lower());
   size_t buf_begin = 0;
#else
   NULL_USE(op);
   bool on_host = arraydata.dataOnHost();
#endif

   /*
    * Loop over the depth sections of the data arrays.
    */

   for (unsigned int d = 0; d < array_d_depth; ++d) {

#if defined(HAVE_RAJA)
      SumOperation<double> sumop_dbl;
      const hier::Box& dst_box = src_is_buffer ? arraydata.getBox() : opbox;
      const hier::Box& src_box = src_is_buffer ? opbox : arraydata.getBox();

      const int dst_offset = dst_box.size() * d;
      const int src_offset = src_box.size() * d;

      switch (dim.getValue()) {
         case 1: {
            typename pdat::ArrayData<double>::template View<1> dest((double*)(dst_ptr + dst_offset), dst_box);
            typename pdat::ArrayData<double>::template View<1> source((double*)(src_ptr + src_offset), src_box);
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i))[1];
                     sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i))[1];
                     sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            }
         } break;

         case 2: {
            typename pdat::ArrayData<double>::template View<2> dest((double*)(dst_ptr + dst_offset), dst_box);
            typename pdat::ArrayData<double>::template View<2> source((double*)src_ptr + src_offset, src_box);
            if (on_host) {
               hier::host_parallel_for_all(opbox, [=] (int i, int j) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i,j))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i,j))[1];
                  sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i, int j) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i,j))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i,j))[1];
                  sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            }
         } break;

         case 3: {
            typename pdat::ArrayData<double>::template View<3> dest((double*)(dst_ptr + dst_offset), dst_box);
            typename pdat::ArrayData<double>::template View<3> source((double*)(src_ptr + src_offset), src_box);
            if (on_host) { 
               hier::host_parallel_for_all(opbox, [=] (int i, int j, int k) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j,k))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j,k))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i,j,k))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i,j,k))[1];
                  sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            } else {
               hier::parallel_for_all(opbox, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  double &dest_real = reinterpret_cast<double(&)[2]>(dest(i,j,k))[0];
                  double &dest_imag = reinterpret_cast<double(&)[2]>(dest(i,j,k))[1];
                  double &source_real = reinterpret_cast<double(&)[2]>(source(i,j,k))[0];
                  double &source_imag = reinterpret_cast<double(&)[2]>(source(i,j,k))[1];
                  sumop_dbl(dest_real, source_real);
                  sumop_dbl(dest_imag, source_imag);
               });
            }
         } break;

         default:
            tbox::perr << "Dimension > 3 not supported with RAJA parallel kernels" << std::endl;
            TBOX_ERROR("Aborting in " __FILE__);
            break;
      }
#else  // !HAVE_RAJA
      size_t dat_counter = dat_begin;
      size_t buf_counter = buf_begin;

      size_t& dst_counter = (src_is_buffer ? dat_counter : buf_counter);
      size_t& src_counter = (src_is_buffer ? buf_counter : dat_counter);

      std::array<int,SAMRAI::MAX_DIM_VAL> dat_b = {};
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dat_b[nd] = static_cast<int>(dat_counter);
      }

      /*
       * Loop over each contiguous block of data.
       */

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            op(dst_ptr[dst_counter + i0], src_ptr[src_counter + i0]);
         }
         int dim_jump = 0;

         /*
          * After each contiguous block is packed, calculate the
          * beginning array index for the next block.
          */

         for (int j = 1; j < dim.getValue(); ++j) {
            if (dim_counter[j] < box_w[j] - 1) {
               ++dim_counter[j];
               dim_jump = j;
               break;
            } else {
               dim_counter[j] = 0;
            }
         }

         if (dim_jump > 0) {

            int dat_step = 1;
            for (int k = 0; k < dim_jump; ++k) {
               dat_step *= dat_w[k];
            }
            dat_counter = dat_b[dim_jump - 1] + dat_step;

            for (int m = 0; m < dim_jump; ++m) {
               dat_b[m] = static_cast<int>(dat_counter);
            }

         }  // if dim_jump > 0

         buf_counter += buf_offset;

      }  // nb loop over contiguous data blocks

      /*
       * After packing is complete on a full box for one depth index,
       * advance by the offset value.
       */

      dat_begin += dat_offset;
      buf_begin = buf_counter;
#endif

   }  // d loop over depth indices

}  // end doArrayDataBufferOperationOnBox

}  // namespace pdat
}  // namespace SAMRAI

#endif
