/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated edge centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_EdgeData
#define included_pdat_EdgeData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/pdat/EdgeIterator.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class EdgeData<TYPE> provides an implementation for data defined
 * at cell edges on AMR patches.  It is derived from the hier::PatchData
 * interface common to all SAMRAI patch data types.  Given a CELL-centered
 * AMR index space box, an edge data object represents data of some template
 * TYPE and depth on the edges of the cells in the box.  Here, depth
 * indicates the number of data values at each edge index location.  The
 * EdgeGeometry class provides the translation between the standard SAMRAI
 * cell-centered AMR index space and edge-centered data.
 *
 * Edge data is stored in DIM arrays, each of which contains the
 * data for the edges tangent to a corresponding coordinate direction.
 * Within each array, data is stored in (i,...,k,d) order, where i,...,k
 * indicates a spatial index and the d indicates the component depth at
 * that locaion.  Memory allocation is in column-major ordering (e.g., Fortran
 * style) so that the leftmost index runs fastest in memory.
 * For example, a three-dimensional edge data object created over a
 * CELL-centered AMR index space [l0:u0,l1:u1,l2:u2] allocates three data
 * arrays sized as follows:
 * \verbatim
 *
 * axis 0
 *   [ l0 : u0 ,
 *     l1 : u1+1 ,
 *     l2 : u2+1 , d ]   ,
 *
 * axis 1
 *   [ l0 : u0+1 ,
 *     l1 : u1 ,
 *     l2 : u2+1 , d ]   ,
 *
 * axis 2
 *   [ l0 : u0+1 ,
 *     l1 : u1+1 ,
 *     l2 : u2 , d ]   ,
 *
 * \endverbatim
 * Here the axis directions 0, 1, 2 can be thought of as the x, y, and z
 * edge directions, respectively, and d is the depth index (i.e., number
 * of values at each face index location).  Other spatial dimensions are
 * represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see EdgeDataFactory
 * @see EdgeIndex
 * @see EdgeIterator
 * @see EdgeGeometry
 */

template<class TYPE>
class EdgeData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent edge-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            edge data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the edge data will be allocated.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre depth > 0
    */
   static size_t
   getSizeOfData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts);

   /*!
    * @brief The constructor for an edge data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            edge data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the edge data will be allocated.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   EdgeData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts);

   /*!
    * @brief The constructor for an edge data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            edge data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the edge data will be allocated.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   EdgeData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief The virtual destructor for an edge data object.
    */
   virtual ~EdgeData();

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular axis and
    * depth component of the edge centered array.
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE *
   getPointer(
      int axis,
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular axis and
    * depth component of the edge centered array.
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE *
   getPointer(
      int axis,
      int depth = 0) const;

#if defined(HAVE_RAJA)
  template <int DIM>
  using View = pdat::ArrayView<DIM, TYPE>;

  template <int DIM>
  using ConstView = pdat::ArrayView<DIM, const TYPE>;

   /*!
    * @brief Get an ArrayView that can access the array for RAJA looping.
    */
  template <int DIM>
  View<DIM> getView(int axis, int depth = 0);

   /*!
    * @brief Get a const ArrayView that can access the array for RAJA looping.
    */
  template <int DIM>
  ConstView<DIM> getConstView(int axis, int depth = 0) const;
#endif

   /*!
    * @brief Return a reference to the data entry corresponding
    * to a given edge index and depth.
    *
    * @pre getDim() == i.getDim()
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE&
   operator () (
      const EdgeIndex& i,
      int depth = 0);

   /*!
    * @brief Return a const reference to the data entry corresponding
    * to a given edge index and depth.
    *
    * @pre getDim() == i.getDim()
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE&
   operator () (
      const EdgeIndex& i,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for the
    * given axis of the edge centered data object.
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    */
   ArrayData<TYPE>&
   getArrayData(
      int axis);

   /*!
    * @brief Return a const reference to the array data object for the
    * given axis of the edge centered data object.
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    */
   const ArrayData<TYPE>&
   getArrayData(
      int axis) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * an EdgeData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    *
    * @pre getDim() == src.getDim()
    */
   virtual void
   copy(
      const hier::PatchData& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * an EdgeData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<EdgeData<TYPE> *>(&dst) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be EdgeData of the same DIM and TYPE
    * and the overlap must be a EdgeOverlap of the same DIM. If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be EdgeData of the same DIM and TYPE
    * and the overlap must be a EdgeOverlap of the same DIM.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<EdgeData<TYPE> *>(&dst) != 0
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given CELL-centered AMR index box.
    *
    * @pre (getDim() == src.getDim()) && (getDim() == box.getDim())
    */
   void
   copyOnBox(
      const EdgeData<TYPE>& src,
      const hier::Box& box);

   /*!
    * @brief Fast copy (i.e., source and this edge data objects are
    * defined over the same box) to this destination edge data object
    * from the given source edge data object at the specified depths.
    *
    * @pre getDim() == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const EdgeData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box).
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual bool
   canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Return the number of bytes needed to stream the data
    * in this patch data object lying in the specified box overlap
    * region.
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    *
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * EdgeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * EdgeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);


   /*!
    * @brief Fill all values at depth d with the value t.
    *
    * @pre (d >= 0) && (d < getDepth())
    */
   void
   fill(
      const TYPE& t,
      int d = 0);

   /*!
    * @brief Fill all values at depth d within the box with the value t.
    *
    * @pre getDim() == box.getDim()
    * @pre (d >= 0) && (d < getDepth())
    */
   void
   fill(
      const TYPE& t,
      const hier::Box& box,
      int d = 0);

   /*!
    * @brief Fill all depth components with value t.
    */
   void
   fillAll(
      const TYPE& t);

   /*!
    * @brief Fill all depth components within the box with value t.
    *
    * @pre getDim() == box.getDim()
    */
   void
   fillAll(
      const TYPE& t,
      const hier::Box& box);

   /*!
    * @brief Print all edge data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void
   print(
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all edge data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    */
   void
   print(
      const hier::Box& box,
      int depth,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all edge centered data values for specified axis
    * residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (axis >= 0) && (axis < getDim().getValue())
    */
   void
   printAxis(
      int axis,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all edge centered data values for specified axis
    * residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    * @pre (axis >= 0) && (axis < getDim().getValue())
    */
   void
   printAxis(
      int axis,
      const hier::Box& box,
      int depth,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * Check that class version and restart file version are equal.  If so,
    * read data members from the restart database.
    *
    * @pre restart_db
    */
   virtual void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * Write out the class version number and other data members to
    * the restart database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * The edge iterator iterates over the elements on one axis of an edge
    * centered box geometry.  This typedef is a convenience for
    * using the EdgeIterator class.
    */
   typedef EdgeIterator iterator;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_EDGEDATA_VERSION;

   // Unimplemented copy constructor
   EdgeData(
      const EdgeData&);

   // Unimplemented assignment operator
   EdgeData&
   operator = (
      const EdgeData&);

   void
   copyWithRotation(
      const EdgeData<TYPE>& src,
      const EdgeOverlap& overlap);

   void
   packWithRotation(
      tbox::MessageStream& stream,
      const EdgeOverlap& overlap) const;

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> > d_data[SAMRAI::MAX_DIM_VAL];

};

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename EdgeData<TYPE>::template View<DIM> get_view(EdgeData<TYPE>& data,
                                                     Args&&... args);

template <int DIM, typename TYPE, typename... Args>
typename EdgeData<TYPE>::template ConstView<DIM> get_const_view(
    const EdgeData<TYPE>& data,
    Args&&... args);
#endif


}
}

#include "SAMRAI/pdat/EdgeData.cpp"

#endif
