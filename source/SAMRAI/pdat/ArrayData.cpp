/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayData_C
#define included_pdat_ArrayData_C

#include "SAMRAI/tbox/KernelFuser.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/ForAll.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/ArrayDataOperationUtilities.h"
#include "SAMRAI/pdat/CopyOperation.h"
#include "SAMRAI/pdat/SumOperation.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

#ifdef HAVE_UMPIRE
#include "umpire/ResourceManager.hpp"
#endif

#include <array>
#include <utility>


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 * Note on the usage of RAJA in this implementation:
 *
 * To support the allocation of ArrayData on either the host or device when
 * running on GPU architectures, each RAJA kernel is implemented twice in
 * an if/else code block.  The looping structure host_parallel_for_all()
 * is used when ArrayData is allocated on the host and guarantees that
 * the loop will execute on the host.  When the GPU device is available and
 * ArrayData is allocated on the device, parallel_for_all() is used to
 * launch the kernels on the device.
 */

namespace SAMRAI
{
namespace pdat
{

template <class TYPE>
const int ArrayData<TYPE>::PDAT_ARRAYDATA_VERSION = 1;

template <class TYPE>
bool ArrayData<TYPE>::canEstimateStreamSizeFromBox()
{
   if ((typeid(TYPE) == typeid(bool)) || (typeid(TYPE) == typeid(char)) || (typeid(TYPE) == typeid(double)) || (typeid(TYPE) == typeid(float)) || (typeid(TYPE) == typeid(int)) || (typeid(TYPE) == typeid(dcomplex))) {
      return true;
   } else {
      return false;
   }
}

template <class TYPE>
size_t
ArrayData<TYPE>::getSizeOfData(
    const hier::Box& box,
    unsigned int depth)
{
   return tbox::MemoryUtilities::align(box.size() * depth * sizeof(TYPE));
}

/*
 *************************************************************************
 *
 * The main constructor allocates data for the given box and depth.  It
 * does not initialize the memory.  The destructor automatically
 * deallocates memory via the array destructor.
 *
 *************************************************************************
 */

template <class TYPE>
ArrayData<TYPE>::ArrayData(
    const hier::Box& box,
    unsigned int depth) : d_depth(depth),
                          d_offset(box.size()),
                          d_box(box)
#if defined(HAVE_UMPIRE)
                          ,
                          d_allocator(umpire::ResourceManager::getInstance().getAllocator("samrai::data_allocator")),
                          d_array(d_allocator.allocate(d_depth * d_offset)),
#else
                          ,
                          d_array(d_depth * d_offset),
#endif
                          d_on_host(true),
                          d_use_fuser(false)
{
   TBOX_ASSERT(depth > 0);

#if defined (HAVE_UMPIRE)
   tbox::ResourceAllocator data_allocator =
      umpire::ResourceManager::getInstance().getAllocator("samrai::data_allocator");
   if (data_allocator.getPlatform() != umpire::Platform::host) {
      d_on_host = false;
   }
#endif

#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

template<class TYPE>
ArrayData<TYPE>::ArrayData(
   const hier::Box& box,
   unsigned int depth,
   tbox::ResourceAllocator allocator):
   d_depth(depth),
   d_offset(box.size()),
   d_box(box),
#if defined(HAVE_UMPIRE)
   d_allocator(allocator),
   d_array(d_allocator.allocate(d_depth * d_offset)),
#else
   d_array(d_depth * d_offset),
#endif
   d_on_host(true),
   d_use_fuser(false)
{
#ifndef HAVE_UMPIRE
   NULL_USE(allocator);
#endif

   TBOX_ASSERT(depth > 0);
#if defined (HAVE_UMPIRE)
   if (allocator.getPlatform() != umpire::Platform::host) {
      d_on_host = false;
   }
#endif

#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

template <class TYPE>
ArrayData<TYPE>::~ArrayData()
{
#if defined(HAVE_UMPIRE)
   d_allocator.deallocate(d_array, d_depth * d_offset);
#endif
}


template <class TYPE>
bool ArrayData<TYPE>::isInitialized() const
{
   return d_depth * d_offset > 0;
}

template <class TYPE>
const hier::Box&
ArrayData<TYPE>::getBox() const
{
   return d_box;
}

template <class TYPE>
unsigned int
ArrayData<TYPE>::getDepth() const
{
   return d_depth;
}

template <class TYPE>
size_t
ArrayData<TYPE>::getOffset() const
{
   return d_offset;
}

template <class TYPE>
size_t
ArrayData<TYPE>::getIndex(
    const hier::Index& i,
    unsigned int d) const
{
   TBOX_ASSERT((d < d_depth));

   size_t index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index < d_depth * d_offset));

   return index;
}

template <class TYPE>
TYPE* ArrayData<TYPE>::getPointer(
    unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   return &d_array[d * d_offset];
}

template <class TYPE>
const TYPE*
ArrayData<TYPE>::getPointer(
    unsigned int d) const
{
   TBOX_ASSERT((d < d_depth));

   return &d_array[d * d_offset];
}

#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename ArrayData<TYPE>::template View<DIM>
ArrayData<TYPE>::getView(
    int depth)
{
   return ArrayData<TYPE>::View<DIM>(getPointer(depth), getBox());
}

template <class TYPE>
template <int DIM>
typename ArrayData<TYPE>::template ConstView<DIM>
ArrayData<TYPE>::getConstView(
    int depth) const
{
   return ArrayData<TYPE>::ConstView<DIM>(getPointer(depth), getBox());
}
#endif

template <class TYPE>
TYPE& ArrayData<TYPE>::operator()(
    const hier::Index& i,
    unsigned int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((d < d_depth));

   return d_array[getIndex(i, d)];
}

template <class TYPE>
const TYPE&
ArrayData<TYPE>::operator()(
    const hier::Index& i,
    unsigned int d) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((d < d_depth));

   return d_array[getIndex(i, d)];
}

/*
 *************************************************************************
 *
 * Copy data between two array data objects on a specified box domain.
 * Don't use C++ indexing member functions, since compilers are probably
 * too stupid to do strength reduction on the loops to get performance.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast copy
 * of all data.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::copy(
    const ArrayData<TYPE>& src,
    const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      TYPE* const dst_ptr = &d_array[0];
      const TYPE* const src_ptr = &src.d_array[0];
      const size_t n = d_offset * d_depth;
#if defined(HAVE_RAJA)
      tbox::KernelFuser* fuser = d_use_fuser ?
         tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;
      if (d_on_host) {
         hier::host_parallel_for_all(0, n, [=] (int i) {
            copyop(dst_ptr[i], src_ptr[i]);
         });
      } else {
         hier::parallel_for_all(fuser, 0, n, [=] SAMRAI_HOST_DEVICE(int i) {
            copyop(dst_ptr[i], src_ptr[i]);
         });
      } 
#else
      for (size_t i = 0; i < n; ++i) {
         copyop(dst_ptr[i], src_ptr[i]);
      }
#endif

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       copyop);
      }
   }
}

/*
 *************************************************************************
 *
 * Copy data from source ArrayData object to this (destination)
 * ArrayData object on given box domain.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::copy(
    const ArrayData<TYPE>& src,
    const hier::Box& box,
    const hier::IntVector& src_shift)
{

   if (src_shift == hier::IntVector::getZero(box.getDim())) {

      copy(src, box);

   } else {

      const hier::Box copybox =
          box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         CopyOperation<TYPE> copyop;

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       copyop);
      }
   }
}

template <class TYPE>
void ArrayData<TYPE>::copy(
    const ArrayData<TYPE>& src,
    const hier::Box& box,
    const hier::Transformation& transformation)
{
   if (transformation.getRotation() == hier::Transformation::NO_ROTATE && transformation.getOffset() == hier::IntVector::getZero(box.getDim()) && transformation.getBeginBlock() == transformation.getEndBlock()) {

      copy(src, box);

   } else {

      hier::Box transformed_src(src.d_box);
      transformation.transform(transformed_src);
      const hier::Box copybox(
          box * d_box * transformed_src);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         CopyOperation<TYPE> copyop;

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       transformation.getOffset(),
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       copyop);
      }
   }
}

/*
 *************************************************************************
 *
 * Copy over the boxlist by calling the single-box copy for each box in
 * the boxlist.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::copy(
    const ArrayData<TYPE>& src,
    const hier::BoxContainer& boxes,
    const hier::IntVector& src_shift)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      copy(src, *b, src_shift);
   }
}

template <class TYPE>
void ArrayData<TYPE>::copy(
    const ArrayData<TYPE>& src,
    const hier::BoxContainer& boxes,
    const hier::Transformation& transformation)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      copy(src, *b, transformation);
   }
}

/*
 *************************************************************************
 *
 * Copy data between two array data objects on a specified box domain.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast copy
 * of all data.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::copyDepth(
    unsigned int dst_depth,
    const ArrayData<TYPE>& src,
    unsigned int src_depth,
    const hier::Box& box)
{
   TBOX_ASSERT((dst_depth <= d_depth));
   TBOX_ASSERT((src_depth <= src.d_depth));

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_box.isSpatiallyEqual(src.d_box)) && (box.isSpatiallyEqual(d_box))) {

      TYPE* const dst_ptr = &d_array[0];
      const TYPE* const src_ptr = &src.d_array[0];

      TYPE* const dst_ptr_d = dst_ptr + dst_depth * d_offset;
      const TYPE* const src_ptr_d = src_ptr + src_depth * d_offset;


#if defined(HAVE_RAJA)
      tbox::KernelFuser* fuser = d_use_fuser ?
         tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;

      if (d_on_host) {
         hier::host_parallel_for_all(0, d_offset, [=] (int i) {
            copyop(dst_ptr_d[i], src_ptr_d[i]);
         });
      } else {
         hier::parallel_for_all(fuser, 0, d_offset, [=] SAMRAI_HOST_DEVICE(int i) {
            copyop(dst_ptr_d[i], src_ptr_d[i]);
         });
      }
#else
      for (size_t i = 0; i < d_offset; ++i) {
         copyop(dst_ptr_d[i], src_ptr_d[i]);
      }
#endif
   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = dst_depth;
         const unsigned int src_start_depth = src_depth;
         const unsigned int num_depth = 1;
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       copyop);
      }
   }
}

/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on given box region.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast sum
 * on all data rather than performing explicit looping operations.
 *
 *************************************************************************
 */
template <class TYPE>
void ArrayData<TYPE>::sum(
    const ArrayData<TYPE>& src,
    const hier::Box& box)
{

   SumOperation<TYPE> sumop;

   /*
    * Do a fast copy and add if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      TYPE* const dst_ptr = &d_array[0];
      const TYPE* const src_ptr = &src.d_array[0];
      const size_t n = d_offset * d_depth;

#if defined(HAVE_RAJA)
      if (d_on_host) {
         hier::host_parallel_for_all(0, n, [=] (int i) {
            sumop(dst_ptr[i], src_ptr[i]);
         });
      } else { 
         hier::parallel_for_all(0, n, [=] SAMRAI_HOST_DEVICE(int i) {
            sumop(dst_ptr[i], src_ptr[i]);
         });
      }
#else
      for (size_t i = 0; i < n; ++i) {
         sumop(dst_ptr[i], src_ptr[i]);
      }
#endif

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       sumop);
      }
   }
}

////////////////////////////////
// specialization for dcomplex
template <>
inline void ArrayData<dcomplex>::sum(
    const ArrayData<dcomplex>& src,
    const hier::Box& box)
{

   SumOperation<dcomplex> sumop;

   /*
    * Do a fast copy and add if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      dcomplex* const dst_ptr = &d_array[0];
      const dcomplex* const src_ptr = &src.d_array[0];
      const size_t n = d_offset * d_depth;

#if defined(HAVE_RAJA)
      SumOperation<double> sumop_dbl;
      if (d_on_host) {
         hier::host_parallel_for_all(0, n, [=] (int i) {   
            double &dst_ptr_real =
               reinterpret_cast<double(&)[2]>(dst_ptr[i])[0];
            double &dst_ptr_imag =
               reinterpret_cast<double(&)[2]>(dst_ptr[i])[1];
            const double &src_ptr_real =
               reinterpret_cast<const double(&)[2]>(src_ptr[i])[0];
            const double &src_ptr_imag =
               reinterpret_cast<const double(&)[2]>(src_ptr[i])[1];

            sumop_dbl(dst_ptr_real, src_ptr_real);
            sumop_dbl(dst_ptr_imag, src_ptr_imag);
         });
      } else {
         hier::parallel_for_all(0, n, [=] SAMRAI_HOST_DEVICE(int i) {   
            double &dst_ptr_real = 
               reinterpret_cast<double(&)[2]>(dst_ptr[i])[0];
            double &dst_ptr_imag = 
               reinterpret_cast<double(&)[2]>(dst_ptr[i])[1];
            const double &src_ptr_real =
               reinterpret_cast<const double(&)[2]>(src_ptr[i])[0];
            const double &src_ptr_imag = 
               reinterpret_cast<const double(&)[2]>(src_ptr[i])[1];

            sumop_dbl(dst_ptr_real, src_ptr_real);
            sumop_dbl(dst_ptr_imag, src_ptr_imag);
         });
      }
#else
      for (size_t i = 0; i < n; ++i) {
         sumop(dst_ptr[i], src_ptr[i]);
      }
#endif

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<dcomplex, SumOperation<dcomplex> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       sumop);
      }
   }
}


/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on region described by given box and offset.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::sum(
    const ArrayData<TYPE>& src,
    const hier::Box& box,
    const hier::IntVector& src_shift)
{

   if (src_shift == hier::IntVector::getZero(box.getDim())) {

      sum(src, box);

   } else {

      const hier::Box copybox =
          box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         SumOperation<TYPE> sumop;

         ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
             doArrayDataOperationOnBox(*this,
                                       src,
                                       copybox,
                                       src_shift,
                                       dst_start_depth,
                                       src_start_depth,
                                       num_depth,
                                       sumop);
      }
   }
}

/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on regions described by given boxes and offset.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::sum(
    const ArrayData<TYPE>& src,
    const hier::BoxContainer& boxes,
    const hier::IntVector& src_shift)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      sum(src, *b, src_shift);
   }
}

template <class TYPE>
size_t
ArrayData<TYPE>::getDataStreamSize(
    const hier::BoxContainer& boxes,
    const hier::IntVector& source_shift) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(source_shift);
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source_shift);

   const size_t nelements = boxes.getTotalSizeOfBoxes();

   if (typeid(TYPE) == typeid(bool)) {
      return tbox::MessageStream::getSizeof<bool>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(char)) {
      return tbox::MessageStream::getSizeof<char>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(dcomplex)) {
      return tbox::MessageStream::getSizeof<dcomplex>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(double)) {
      return tbox::MessageStream::getSizeof<double>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(float)) {
      return tbox::MessageStream::getSizeof<float>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(int)) {
      return tbox::MessageStream::getSizeof<int>(d_depth * nelements);
   }

   TBOX_ERROR("ArrayData::getDataStreamSize() -- Invalid type" << std::endl);
   return 0;
}

/*
 *************************************************************************
 *
 * Pack data into the message stream.  Both packing routines add one
 * level of copy into a temporary buffer to reduce the number of calls
 * to the abstract stream packing routines.  These definitions will only
 * work for the standard built-in types of bool, char, double, float,
 * and int.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::packStream(
    tbox::MessageStream& stream,
    const hier::Box& dest_box,
    const hier::IntVector& src_shift) const
{
   RANGE_PUSH("ArrayData::pack", 2);

   const size_t size = d_depth * dest_box.size();

   TYPE* buffer = stream.getWriteBuffer<TYPE>(size);
   packBuffer(&buffer[0], hier::Box::shift(dest_box, -src_shift));

   // std::vector< TYPE, umpire::TypedAllocator<TYPE> > buffer(size);
   // stream.pack(&buffer[0], size);

   RANGE_POP
}

template <class TYPE>
void ArrayData<TYPE>::packStream(
    tbox::MessageStream& stream,
    const hier::BoxContainer& dest_boxes,
    const hier::IntVector& src_shift) const
{
   RANGE_PUSH("ArrayData::pack", 2);

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   TYPE* buffer = stream.getWriteBuffer<TYPE>(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      packBuffer(&buffer[ptr], hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   // stream.pack(&buffer[0], size);

   RANGE_POP
}

template <class TYPE>
void ArrayData<TYPE>::packStream(
    tbox::MessageStream& stream,
    const hier::Box& dest_box,
    const hier::Transformation& transformation) const
{
   RANGE_PUSH("ArrayData::pack", 2);

   const size_t size = d_depth * dest_box.size();
   TYPE* buffer = stream.getWriteBuffer<TYPE>(size);

   hier::Box pack_box(dest_box);
   transformation.inverseTransform(pack_box);
   packBuffer(&buffer[0], pack_box);
   //      hier::Box::shift(dest_box, -src_shift));

   //stream.pack(&buffer[0], size);

   RANGE_POP
}

template <class TYPE>
void ArrayData<TYPE>::packStream(
    tbox::MessageStream& stream,
    const hier::BoxContainer& dest_boxes,
    const hier::Transformation& transformation) const
{
   RANGE_PUSH("ArrayData::pack", 2);


   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   TYPE* buffer = stream.getWriteBuffer<TYPE>(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      hier::Box pack_box(*b);
      transformation.inverseTransform(pack_box);
      packBuffer(&buffer[ptr], pack_box);
      //         hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   // stream.pack(&buffer[0], size);

   RANGE_POP
}

/*
 *************************************************************************
 *
 * Unpack data from the message stream.  Both unpacking routines add one
 * level of copy into a temporary buffer to reduce the number of calls
 * to the abstract stream packing routines.  These definitions will only
 * work for the standard built-in types of bool, char, double, float,
 * and int.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::unpackStream(
    tbox::MessageStream& stream,
    const hier::Box& dest_box,
    const hier::IntVector& src_shift)
{
   RANGE_PUSH("ArrayData::unpack", 2);

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_box.size();
   const TYPE* buffer = stream.getReadBuffer<TYPE>(size);

   unpackBuffer(&buffer[0], dest_box);

   RANGE_POP;
}

template <class TYPE>
void ArrayData<TYPE>::unpackStream(
    tbox::MessageStream& stream,
    const hier::BoxContainer& dest_boxes,
    const hier::IntVector& src_shift)
{
   RANGE_PUSH("ArrayData::unpack", 2);

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   const TYPE* buffer = stream.getReadBuffer<TYPE>(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      unpackBuffer(&buffer[ptr], *b);
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   RANGE_POP;
}

/*
 *************************************************************************
 *
 * Unpack data from the message stream and add to this array data object.
 * Both unpacking routines add one level of copy into a temporary buffer
 * to reduce the number of calls to the abstract stream packing routines.
 * These definitions will only work for the standard built-in types of
 * bool, char, double, float, and int.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::unpackStreamAndSum(
    tbox::MessageStream& stream,
    const hier::Box& dest_box,
    const hier::IntVector& src_shift)
{
   NULL_USE(src_shift);

   const size_t size = d_depth * dest_box.size();
   const TYPE* buffer = stream.getReadBuffer<TYPE>(size);

   unpackBufferAndSum(&buffer[0], dest_box);
}

template <class TYPE>
void ArrayData<TYPE>::unpackStreamAndSum(
    tbox::MessageStream& stream,
    const hier::BoxContainer& dest_boxes,
    const hier::IntVector& src_shift)
{
   NULL_USE(src_shift);

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   const TYPE* buffer = stream.getReadBuffer<TYPE>(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      unpackBufferAndSum(&buffer[ptr], *b);
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);
}

/*
 *************************************************************************
 *
 * Fill all or portions of the array with the specified data value.
 * The templated TYPE must define the assignment operator.
 *
 *************************************************************************
 */
template <class TYPE>
void ArrayData<TYPE>::fillAll(
    const TYPE& t)
{
   if (!d_box.empty()) {
      TYPE* ptr = &d_array[0];
      const size_t n = d_depth * d_offset;
#if defined(HAVE_RAJA)
      tbox::KernelFuser* fuser = d_use_fuser ?
         tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;

      if (d_on_host) {
         hier::host_parallel_for_all(0, n, [=] (int i) {
            ptr[i] = t;
         });
      } else {
         hier::parallel_for_all(fuser, 0, n, [=] SAMRAI_HOST_DEVICE(int i) {
            ptr[i] = t;
         });
      }
#if defined(DEBUG_INITIALIZE_UNDEFINED)      
      tbox::parallel_synchronize();
#endif
#else
      for (size_t i = 0; i < n; ++i) {
         ptr[i] = t;
      }
#endif
   }

}

template <class TYPE>
void ArrayData<TYPE>::fillAll(
    const TYPE& t,
    const hier::Box& box)
{
   for (tbox::Dimension::dir_t d = 0; d < d_depth; ++d) {
      fill(t, box, d);
   }
}

template <class TYPE>
void ArrayData<TYPE>::fill(
    const TYPE& t,
    const unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   TYPE* ptr = &d_array[d * d_offset];
   const size_t n = d_offset;
   if (!d_box.empty()) {
#if defined(HAVE_RAJA)
      tbox::KernelFuser* fuser = d_use_fuser ?
         tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;

      if (d_on_host) {
         hier::host_parallel_for_all(0, n, [=] (int i) {
            ptr[i] = t;
         });
      } else {
         hier::parallel_for_all(fuser, 0, n, [=] SAMRAI_HOST_DEVICE(int i) {
            ptr[i] = t;
         });
      }
#else
      for (size_t i = 0; i < n; ++i) {
         ptr[i] = t;
      }
#endif
   }

}


template <class TYPE>
void ArrayData<TYPE>::fill(
    const TYPE& t,
    const hier::Box& box,
    const unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   const hier::Box ispace = d_box * box;

   if (!ispace.empty()) {
#if defined(HAVE_RAJA)
      tbox::KernelFuser* fuser = d_use_fuser ?
         tbox::StagedKernelFusers::getInstance()->getFuser(0) : nullptr;

      switch (ispace.getDim().getValue()) {
         case 1: {
            auto data = getView<1>(d);
            if (d_on_host) {
               hier::host_parallel_for_all(ispace, [=] (int i) {
                  data(i) = t;
               });
            } else {
               hier::parallel_for_all(fuser, ispace, [=] SAMRAI_HOST_DEVICE(int i) {
                  data(i) = t;
               });
            }
            break;
         }
         case 2: {
            auto data = getView<2>(d);
            if (d_on_host) {
               hier::host_parallel_for_all(ispace, [=] (int i, int j) {
                  data(i,j) = t;
               });
            } else {
               hier::parallel_for_all(fuser, ispace, [=] SAMRAI_HOST_DEVICE(int i, int j) {
                  data(i,j) = t;
               });
            }
            break;
         }
         case 3: {
            auto data = getView<3>(d);
            if (d_on_host) {
               hier::host_parallel_for_all(ispace, [=] (int i, int j, int k) {
                  data(i,j,k) = t;
               });
            } else {
               hier::parallel_for_all(fuser, ispace, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  data(i,j,k) = t;
               });
            }
            break;
         }
         default:
            TBOX_ERROR("hier::parallel_for_all undefined for dim > 3" << std::endl);
      }
#else
      const tbox::Dimension& dim = box.getDim();

      std::array<int,SAMRAI::MAX_DIM_VAL> box_w = {};
      std::array<int,SAMRAI::MAX_DIM_VAL> dst_w = {};
      std::array<int,SAMRAI::MAX_DIM_VAL> dim_counter = {};
      for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
         box_w[i] = ispace.numberCells(i);
         dst_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int num_d0_blocks = static_cast<int>(ispace.size() / box_w[0]);

      size_t dst_counter = d_box.offset(ispace.lower()) + d * d_offset;

      std::array<size_t,SAMRAI::MAX_DIM_VAL> dst_b = {};
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dst_b[nd] = dst_counter;
      }

      TYPE* const dst_ptr = &d_array[0];

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            dst_ptr[dst_counter + i0] = t;
         }
         int dim_jump = 0;

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
            int dst_step = 1;
            for (int k = 0; k < dim_jump; ++k) {
               dst_step *= dst_w[k];
            }
            dst_counter = dst_b[dim_jump - 1] + dst_step;

            for (int m = 0; m < dim_jump; ++m) {
               dst_b[m] = dst_counter;
            }
         }
      }
#endif
   }
}

template <class TYPE>
void ArrayData<TYPE>::fillSequential(
    const TYPE& t,
    const hier::Box& box,
    const unsigned int d)
{
   const tbox::Dimension& dim = box.getDim();
   if (dim.getValue() > 3) {
      fill(t, box, d);
   } else {
      hier::Box fill_box(box * d_box);
      const hier::Index& b_lower = fill_box.lower();
      const hier::Index b_upper = fill_box.upper();
      int b_ilo = b_lower(0);
      int b_ihi = b_upper(0);
      int b_jlo = 0;
      int b_jhi = 0;
      int b_klo = 0;
      int b_khi = 0;
      if (dim.getValue() > 1) {
         b_jlo = b_lower(1);
         b_jhi = b_upper(1);
      }
      if (dim.getValue() == 3) {
         b_klo = b_lower(2);
         b_khi = b_upper(2);
      }

      const hier::Index& dbox_lower = d_box.lower();
      const hier::Index& dbox_upper = d_box.upper();
      int iwidth = dbox_upper(0)-dbox_lower(0) + 1;
      int jwidth = 0;
      int dbox_ilo = dbox_lower(0);
      int dbox_jlo = 0;
      int dbox_klo = 0;
      if (dim.getValue() > 1) {
         jwidth = dbox_upper(1) - dbox_lower(1) + 1;
         dbox_jlo = dbox_lower(1);
      }
      if (dim.getValue() == 3) {
         dbox_klo = dbox_lower(2);
      }

      for (int k = b_klo; k <= b_khi; ++k) { 
         int koffset = (d * d_offset) + (k - dbox_klo) * iwidth * jwidth; 
         for (int j = b_jlo; j <= b_jhi; ++j) {
            int joffset = (j - dbox_jlo) * iwidth;
            for (int i = b_ilo; i <= b_ihi; ++i) {
               int ioffset = i - dbox_ilo;
               d_array[ioffset+joffset+koffset] = t;
            }
         }
      }
   }

}


/*
 *************************************************************************
 *
 * Checks to make sure that class and restart file version numbers are
 * equal.  If so, reads in d_depth, d_offset, and d_box from the
 * database.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::getFromRestart(
    const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   int ver = restart_db->getInteger("PDAT_ARRAYDATA_VERSION");
   if (ver != PDAT_ARRAYDATA_VERSION) {
      TBOX_ERROR("ArrayData::getFromRestart error...\n"
                 << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");
   d_offset = restart_db->getInteger("d_offset");
   d_box = restart_db->getDatabaseBox("d_box");

#if defined(HAVE_UMPIRE)
   std::vector<TYPE> temp;
   restart_db->getVector("d_array", temp);
   std::copy(temp.begin(), temp.end(), d_array);
#else
   restart_db->getVector("d_array", d_array);
#endif
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth, d_offset, d_box, and
 * d_array to the restart database.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::putToRestart(
    const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("PDAT_ARRAYDATA_VERSION", PDAT_ARRAYDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);
   restart_db->putInteger("d_offset", static_cast<int>(d_offset));
   restart_db->putDatabaseBox("d_box", d_box);

#if defined(HAVE_UMPIRE)
   restart_db->putVector("d_array", std::vector<TYPE>(d_array, d_array + d_depth * d_offset));
#else
   restart_db->putVector("d_array", d_array);
#endif
}

template <class TYPE>
const tbox::Dimension&
ArrayData<TYPE>::getDim() const
{
   return d_box.getDim();
}

template <class TYPE>
bool ArrayData<TYPE>::isValid()
{
   return !d_box.empty();
}

/*
 *************************************************************************
 *
 * Set all array data to undefined values.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::undefineData()
{
   fillAll(tbox::MathUtilities<TYPE>::getSignalingNaN());
}

/*
 *************************************************************************
 *
 * Private member functions to pack and unpack data on the specified box
 * (for all components) into/from the buffer.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::packBuffer(
    TYPE* buffer,
    const hier::Box& box) const
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = false;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
       doArrayDataBufferOperationOnBox(*this,
                                       buffer,
                                       box,
                                       src_is_buffer,
                                       copyop);
}

template <class TYPE>
void ArrayData<TYPE>::unpackBuffer(
    const TYPE* buffer,
    const hier::Box& box)
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = true;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
       doArrayDataBufferOperationOnBox(*this,
                                       buffer,
                                       box,
                                       src_is_buffer,
                                       copyop);
}

/*
 *************************************************************************
 *
 * Private member function to unpack data on the specified box
 * (all components) from the buffer and add to this array data object.
 *
 *************************************************************************
 */

template <class TYPE>
void ArrayData<TYPE>::unpackBufferAndSum(
    const TYPE* buffer,
    const hier::Box& box)
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = true;

   SumOperation<TYPE> sumop;

   ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
       doArrayDataBufferOperationOnBox(*this,
                                       buffer,
                                       box,
                                       src_is_buffer,
                                       sumop);
}

#if defined(HAVE_RAJA)
template <int DIM, typename DATA, typename... Args>
typename DATA::template View<DIM> get_view(std::shared_ptr<hier::PatchData> src, Args&&... args)
{
   auto ptr = SAMRAI_SHARED_PTR_CAST<DATA, hier::PatchData>(src);
   return ptr->template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename DATA, typename... Args>
typename DATA::template ConstView<DIM> get_const_view(const std::shared_ptr<hier::PatchData> src, Args&&... args)
{
   auto ptr = SAMRAI_SHARED_PTR_CAST<const DATA, const hier::PatchData>(src);
   return ptr->template getConstView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename ArrayData<TYPE>::template View<DIM> get_view(ArrayData<TYPE>& data, Args&&... args)
{
   return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename ArrayData<TYPE>::template ConstView<DIM> get_const_view(const ArrayData<TYPE>& data, Args&&... args)
{
   return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif


}  // namespace pdat
}  // namespace SAMRAI

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
