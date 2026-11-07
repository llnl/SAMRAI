/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayData
#define included_pdat_ArrayData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayDataIterator.h"
#include "SAMRAI/pdat/ArrayView.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <typeinfo>
#include <vector>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
class OuteredgeData;
template<class TYPE>
class OuternodeData;
template<class TYPE>
class SideData;

/*!
 * @brief Class ArrayData<TYPE> is a basic templated array structure defined
 * over the index space of a box (with a specified depth) that provides
 * the support for the various standard array-based patch data subclasses.
 *
 * The data storage is in (i,...,k,d) order, where i,...,k indicates
 * spatial indices and the d indicates the component at that location.
 * Memory allocation is in column-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.  Note that a number of
 * functions only work for standard built-in types (bool, char, double,
 * float, and int).  To use this class with other user-defined types,
 * many of these functions will need to be specialized, especially those
 * that deal with message packing and unpacking.
 */

template<class TYPE>
class ArrayData
{
public:
   /*!
    * Static member function that returns true when the amount of buffer space
    * in a message stream can be estimated from box only.  For built-in types
    * (bool, char, double, float, int, and dcomplex), this routine returns
    * true.  For other data types (template paramters) that may require special
    * handling, a different implementation must be provided.
    */
   static bool
   canEstimateStreamSizeFromBox();

   /*!
    * Static member function that returns the amount of memory space needed to
    * store data of given depth on a box.
    *
    * Note that this function is only defined for the standard data types:
    * bool, char, double, float, int, and dcomplex.  It must be provided for
    * other template parameter types.
    *
    * @return size_t value indicating the amount of memory space needed for the
    *         data.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the array data index region of interest.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    */
   static size_t
   getSizeOfData(
      const hier::Box& box,
      unsigned int depth);

   /*!
    * Construct an array data object.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    *
    * @pre depth > 0
    */
   ArrayData(
      const hier::Box& box,
      unsigned int depth);

   /*!
    * Construct an array data object using a ResourceAllocator.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    * @param allocator A ResourceAllocator
    *
    * @pre depth > 0
    */
   ArrayData(
      const hier::Box& box,
      unsigned int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * The destructor for an array data object releases all memory allocated
    * for the array elements.
    */
   ~ArrayData();

   /*!
    * @brief Returns true when the array has been properly initialized
    * and storage has been allocated; otherwise, return false.
    *
    * Note: Only arrays that have been initialized can do anything useful.
    */
   bool
   isInitialized() const;

   /*!
    * Set the array data to an ``undefined'' state appropriate for the data
    * type. For example, for float and double, this means setting data to
    * signaling NaNs that cause a floating point assertion when used in a
    * numerical expression without being set to valid values.
    */
   void
   undefineData();

   /*!
    * Return the box over which the array is defined.
    */
   const hier::Box&
   getBox() const;

   /*!
    * Return the depth (e.g., the number of data values at each spatial
    * location) of this array.
    */
   unsigned int
   getDepth() const;

   /*!
    * Return the offset (e.g., the number of data values for each
    * depth component) of this array.
    */
   size_t
   getOffset() const;

   /*!
    * Get a non-const pointer to the beginning of the given depth
    * component of this data array.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   TYPE *
   getPointer(
      const unsigned int d = 0);

   /*!
    * Get a const pointer to the beginning of the given depth
    * component of this data array.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   const TYPE *
   getPointer(
      const unsigned int d = 0) const;

#if defined(HAVE_RAJA)
   template<int DIM>
   using View = pdat::ArrayView<DIM, TYPE>;

   template<int DIM>
   using ConstView = pdat::ArrayView<DIM, const TYPE>;

   /*!
    * @brief Get an ArrayView that can access the array for RAJA looping.
    */
   template <int DIM>
   View<DIM>
   getView(
      int depth = 0);

   /*!
    * @brief Get a const ArrayView that can access the array for RAJA looping.
    */
   template <int DIM>
   ConstView<DIM>
   getConstView(
      int depth = 0) const;
#endif

   /*!
    * Return reference to value in this array associated with the given
    * box index and depth component.
    *
    * @pre getDim() == i.getDim()
    * @pre (d >= 0) && (d < getDepth())
    */
   TYPE&
   operator () (
      const hier::Index& i,
      const unsigned int d);

   /*!
    * Return const reference to value in this array associated with the given
    * box index and depth component.
    *
    * @pre getDim() == i.getDim()
    * @pre (d >= 0) && (d < getDepth())
    */
   const TYPE&
   operator () (
      const hier::Index& i,
      const unsigned int d) const;

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy
    *              operation.
    *              Note: the box is in either the source or destination index
    *                    space (which are assumed to be the same).
    *
    * @pre (getDim() == src.getDim()) && (getDim() == box.getDim())
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box
    * and shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy
    *              operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array
    *              data object.
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::IntVector& src_shift);

   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::Transformation& transformation);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box
    * and shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the copy
    *              operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array
    *              data object.
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift);

   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation);

   /*!
    * Copy given source depth of source array data object to given destination
    * depth of this array data object on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param dst_depth Integer depth of destination array.
    * @param src       Const reference to source array data object.
    * @param src_depth Integer depth of source array.
    * @param box       Const reference to box object describing the spatial
    *                  extents of the index space region over which to perform
    *                  the copy operation.
    *                  Note: the box is in either the source or destination
    *                        index space (which are assumed to be the same).
    *
    * @pre (0 <= dst_depth) && (dst_depth <= getDepth()))
    * @pre (0 <= src_depth) && (src_depth <= src.getDepth())
    */
   void
   copyDepth(
      unsigned int dst_depth,
      const ArrayData<TYPE>& src,
      unsigned int src_depth,
      const hier::Box& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum
    *              operation.
    *              Note: the box is in either the source or destination index
    *                    space (which are assumed to be the same).
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::Box& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box
    * and shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum
    *              operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array
    *              data object.
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::IntVector& src_shift);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box
    * and shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the sum
    *              operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array
    *              data object.
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift);

   /*!
    * Calculate the number of bytes needed to stream the data living
    * in the specified box domains.  This routine is only defined for
    * the built-in types of bool, char, double, float, int, and dcomplex.  For
    * all other types, a specialized implementation must be provided.
    *
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions of interest.
    *              Note: the boxes are assumed to be in the index space of this
    *              array data object.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this array data object.
    *              Note: this argument is currently ignored.
    *
    * @pre (getDim() == src_shift.getDim())
    */
   size_t
   getDataStreamSize(
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift) const;

   /*!
    * Pack data living on the specified index region into the stream.
    *
    * Note that this routine assumes that the given box region must
    * be shifted to be consistent with the source (this) array box region.
    *
    * @param stream Reference to stream into which to pack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              box into the index space region of this (source) array data
    *              object.
    *
    * Note: The shifted box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the shifted box is not contained in the index space of this
    * array.
    */
   void
   packStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_shift) const;

   void
   packStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::Transformation& src_shift) const;

   /*!
    * Pack data living on the specified index regions into the stream.
    *
    * Note that this routine assumes that the given box regions must
    * be shifted to be consistent with the source (this) array box region.
    *
    * @param stream Reference to stream into which to pack data.
    * @param dest_boxes Const reference to boxes describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this (source) array
    *              data object.
    *
    * Note: The shifted boxes must lie completely within the index space of
    * this array.  If compiled with assertions enabled, the routine will abort
    * if the shifted boxes are not contained in the index space of this array.
    */
   void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_shift) const;

   void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::Transformation& transformation) const;

   /*!
    * Unpack data from the stream into the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream into the index regions specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial
    *              extents of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array
    *              data object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if some box is not contained in the index space of this array.
    */
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream and add to the array in the index region
    * specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void
   unpackStreamAndSum(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream and ad to the array in the index region
    * specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial
    *              extents of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array
    *              data object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array.  If compiled with assertions enabled, the routine will abort if
    * some box is not contained in the index space of this array.
    */
   void
   unpackStreamAndSum(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_offset);

   /*!
    * Fill all array values with value t.
    */
   void
   fillAll(
      const TYPE& t);

   /*!
    * Fill all array values within the box with value t.
    */
   void
   fillAll(
      const TYPE& t,
      const hier::Box& box);

   /*!
    * Fill all array values associated with depth component d with the value t.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   void
   fill(
      const TYPE& t,
      const unsigned int d = 0);

   /*!
    * Fill all array values associated with depth component d
    * within the box with the value t.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   void
   fill(
      const TYPE& t,
      const hier::Box& box,
      const unsigned int d = 0);

   /*!
    * @brief Fill all arrray values associated with depth component d
    * with the box with the value t
    *
    * This is equivalent to the fill() method with the same signature, with an
    * implementation that never uses parallel threading.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   void
   fillSequential(
      const TYPE& t,
      const hier::Box& box,
      const unsigned int d = 0);

   /*!
    * Check to make sure that the class version and restart file
    * version are equal.  If so, read in data from restart database.
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * Write out array data object data to restart database.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

   /**
    * @brief Returns true if the array is valid.
    */
   bool
   isValid();

   /*!
    * @brief Returns true if the data is on the CPU host.
    */
   bool dataOnHost() const
   {
      return d_on_host;
   }

   void startKernelFuser()
   {
      d_use_fuser = true;
   }

   void stopKernelFuser()
   {
      d_use_fuser = false;
   }

   bool useFuser() const
   {
      return d_use_fuser;
   }

   /*!
    * The array data iterator iterates over the elements of a box
    * associated with an ArrayData object.  This typedef is
    * convenient link to the ArrayDataIterator class.
    */
   typedef ArrayDataIterator iterator;

private:
   // Unimplemented default constructor.
   ArrayData();

   // Unimplemented copy constructor.
   ArrayData(
      const ArrayData&);

   // Unimplemented assignment operator.
   ArrayData&
   operator = (
      const ArrayData&);

   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_ARRAYDATA_VERSION;

   /*
    * Private member functions to pack/unpack data to/from buffer.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void
   packBuffer(
      TYPE* buffer,
      const hier::Box& box) const;
   void
   unpackBuffer(
      const TYPE* buffer,
      const hier::Box& box);

   /*
    * Private member functions to unpack data from buffer and add to
    * this array data object.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void
   unpackBufferAndSum(
      const TYPE* buffer,
      const hier::Box& box);

   /*!
    * @brief Compte index into d_array for data at index i and depth d.
    *
    * @param i Index of data
    * @param d depth of data
    *
    * @pre (d >= 0) && (d < getDepth())
    *
    * @post (index >= 0) && (index < getDepth() * getOffset())
    */
   size_t
   getIndex(
      const hier::Index& i,
      unsigned int d) const;

   unsigned int d_depth;
   size_t d_offset;
   hier::Box d_box;
#if defined(HAVE_UMPIRE)
   umpire::TypedAllocator<TYPE> d_allocator;
   TYPE* d_array;
#else
   std::vector<TYPE> d_array;
#endif

   bool d_on_host;
   bool d_use_fuser;
};

#if defined(HAVE_RAJA)
template<int DIM, typename DATA, typename... Args>
typename DATA::template View<DIM> get_view(std::shared_ptr<hier::PatchData> src, Args&&... args);

template<int DIM, typename DATA, typename... Args>
typename DATA::template ConstView<DIM> get_const_view(const std::shared_ptr<hier::PatchData> src, Args&&... args);

template<int DIM, typename TYPE, typename... Args>
typename ArrayData<TYPE>::template View<DIM> get_view(ArrayData<TYPE>& data, Args&&... args);

template<int DIM, typename TYPE, typename... Args>
typename ArrayData<TYPE>::template ConstView<DIM> get_const_view(const ArrayData<TYPE>& data, Args&&... args);
#endif

}
}

#include "SAMRAI/pdat/ArrayData.cpp"

#endif
