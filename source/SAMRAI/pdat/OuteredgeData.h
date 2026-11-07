/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outeredge centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeData
#define included_pdat_OuteredgeData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/pdat/OuteredgeGeometry.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OuteredgeData<TYPE> provides an implementation for data defined
 * at cell edges on the boundaries of AMR patches.  It is derived from the
 * hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outeredge data object represents
 * data of some template TYPE and depth on the cell edges on the boundary
 * of the box.  Here, depth indicates the number of data values at each edge
 * index location.  The OuteredgedgeGeometry class provides the translation
 * between the standard SAMRAI cell-centered AMR index space and
 * outeredge-centered data.
 *
 * Outeredge data is stored in 2*DIM*(DIM-2) arrays, each of which contains
 * data associated with edge indices in a coordinate axis direction, an outward
 * pointing face normal direction, and an upper or lower box face in the
 * face normal direction.  The data layout in the outernode data arrays matches
 * the corresponding array sections provided by the node data implementation.
 * Note that outeredge data is NOT defined when the axis and face normal are
 * equal.  This is consistent with the edge data representation.  Where an edge
 * index falls on more than one box face (patch boundary edges and corners),
 * the outeredge data value belongs to only one data array so that there are no
 * redundant data values. Specifically, when DIM > 2, outeredge data boxes are
 * "trimmed" for each axis edge direction so that each edge index that lives on
 * more than one face on the box boundary will be associated with the largest
 * face normal direction and only that face.  Within each array, data is stored
 * in (i,...,k,d) order, where i,...,k indicates a spatial index and the
 * d indicates the component depth at that location.  Memory allocation is
 * in column-major ordering (e.g., Fortran style) so that the leftmost
 * index runs fastest in memory.
 *
 * To illustrate the outeredge data layout, in particular the "box trimming"
 * that prevents redundant data values, we describe the data for a
 * three-dimensional outeredge data object instantiated over a
 * box [l0:u0,l1:u1,l2:u2] in the standard SAMRAI cell-centered AMR index
 * space.
 * Note: no boxes are trimmed when DIM < 3.
 *
 * \verbatim
 *
 *    a = edge axis (corresponds to the standard EdgeData axis)
 *    f = face normal dir
 *    s = lower/upper face
 *    d = data depth (i.e., number of values at each edge index point).
 *
 *    Here axis and face normal values 0, 1, and 2 can be thought of as X, Y, Z
 *    respectively.
 *
 *        (a,f,s)
 *    axis 0:
 *        (0,0,[0,1])  DATA IS NOT DEFINED when a == f
 *        (0,1,0)      [l0:u0, l1:l1,     l2+1:u2,   d]
 *        (0,1,1)      [l0:u0, u1+1:u1+1, l2+1:u2,   d]
 *        (0,2,0)      [l0:u0, l1:u1+1,   l2:l2,     d]
 *        (0,2,1)      [l0:u0, l1:u1+1,   u2+1:u2+1, d]
 *        Note: Edge indices are duplicated at the intersection of faces normal
 *              to directions 1 and 2.
 *              So boxes for face normal direction 1 are trimmed in direction 2
 *              so that edge indices shared with faces in direction 2 only
 *              appear in face normal direction 2 boxes.
 *
 *    axis 1:
 *        (1,0,0)      [l0:l0,     l1:u1, l2+1:u2,   d]
 *        (1,0,1)      [u0+1:u0+1, l1:u1, l2+1:u2,   d]
 *        (1,1,[0,1])  DATA IS NOT DEFINED when a == f
 *        (1,2,0)      [l0:u0+1,   l1:u1, l2:l2,     d]
 *        (1,2,1)      [l0:u0+1,   l1:u1, u2+1:u2+1, d]
 *        Note: Edge indices are duplicated at the intersection of faces normal
 *              to directions 0 and 2.
 *              So boxes for face normal direction 0 are trimmed in direction 2
 *              so that edge indices shared with faces in direction 2 only
 *              appear in face normal direction 2 boxes.
 *
 *    axis 2:
 *        (2,0,0)      [l0:l0,     l1+1:u1,   l2:u2, d]
 *        (2,0,1)      [u0+1:u0+1, l1+1:u1,   l2:u2, d]
 *        (2,1,0)      [l0:u0+1,   l1:l1,     l2:u2, d]
 *        (2,1,1)      [l0:u0+1,   u1+1:u1+1, l2:u2, d]
 *        (2,2,[0,1])  DATA IS NOT DEFINED when a == f
 *        Note: Edge indices are duplicated at the intersection of faces normal
 *              to directions 0 and 1.
 *              So boxes for face normal direction 0 are trimmed in direction 1
 *              so that edge indices shared with faces in direction 1 only
 *              appear in face normal direction 1 boxes.
 *
 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see OuteredgeDataFactory
 * @see OuteredgeGeometry
 * @see EdgeIterator
 * @see EdgeIndex
 */

template<class TYPE>
class OuteredgeData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outeredge-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outeredge data object will be created.
    *            Note: the ghost cell width is assumed to be zero.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    *
    * @pre depth > 0
    */
   static size_t
   getSizeOfData(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Constructor for an outeredge data object.
    *
    * Note: Outeredge data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outeredge data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    *
    * @pre depth > 0
    */
   OuteredgeData(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Constructor for an outeredge data object.
    *
    * Note: Outeredge data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outeredge data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data. 
    *
    *
    * @pre depth > 0
    */
   OuteredgeData(
      const hier::Box& box,
      int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Virtual destructor for a outeredge data object.
    */
   virtual ~OuteredgeData();

   /*!
    * @brief Return the depth (i.e., the number of data values for
    * each spatial location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Returns true if outeredge data exists for the given
    * axis and face normal direction; false otherwise.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    */
   bool
   dataExists(
      int axis,
      int face_normal) const;

   /*!
    * @brief Return the box of valid edge indices for
    *        outeredge data.  Note: the returned box
    *        will reside in the @em edge index space.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   hier::Box
   getDataBox(
      int axis,
      int face_normal,
      int side);

   /*!
    * @brief Get a pointer to the beginning of a particular axis,
    * face normal, side, and depth component of the outeredge centered
    * array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE *
   getPointer(
      int axis,
      int face_normal,
      int side,
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular axis,
    * face normal, side, and depth component of the outeredge centered
    * array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE *
   getPointer(
      int axis,
      int face_normal,
      int side,
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
   View<DIM> getView(int axis, int face_normal, int side, int depth = 0);

   /*!
    * @brief Get a const ArrayView that can access the array for RAJA looping.
    */
   template <int DIM>
   ConstView<DIM> getConstView(int axis, int face_normal, int side, int depth = 0) const;
#endif

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given edge index and depth.
    *
    * @param i const reference to EdgeIndex, @em MUST be
    *          an index on the outeredge of the box.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
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
    * @brief Return a const reference to data entry corresponding
    * to a given edge index and depth.
    *
    * @param i const reference to EdgeIndex, @em MUST be
    *          an index on the outeredge of the box.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
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
    * @brief Return a reference to the array data object for
    * given axis, face normal, and side index of the outeredge
    * centered array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   ArrayData<TYPE>&
   getArrayData(
      int axis,
      int face_normal,
      int side);

   /*!
    * @brief Return a const reference to the array data object for
    * given axis, face normal, and side index of the outeredge
    * centered array.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   const ArrayData<TYPE>&
   getArrayData(
      int axis,
      int face_normal,
      int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * either EdgeData or OuteredgeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const EdgeData<TYPE> *>(&src) != 0 ||
    *      dynamic_cast<const OuteredgeData<TYPE> *>(&src) != 0
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
    * either EdgeData or OuteredgeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<EdgeData<TYPE> *>(&dst) != 0 ||
    *      dynamic_cast<OuteredgeData<TYPE> *>(&dst) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be either EdgeData or OuteredgeData
    * of the same DIM and TYPE and the overlap must be an EdgeOverlap
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    * @pre dynamic_cast<const EdgeData<TYPE> *>(&src) != 0 ||
    *      dynamic_cast<const OuteredgeData<TYPE> *>(&src) != 0
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be either EdgeData or OuteredgeData
    * of the same DIM and TYPE and the overlap must be an EdgeOverlap
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    * @pre dynamic_cast<EdgeData<TYPE> *>(&dst) != 0 ||
    *      dynamic_cast<OuteredgeData<TYPE> *>(&dst) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes edge and outeredge data objects are
    * defined over the same box) from the given edge source data object to
    * this destination outeredge data object at the specified depths.
    *
    * @pre getDim() == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const EdgeData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes edge and outeredge data objects are
    * defined over the same box) to the given edge destination data object
    * from this source outeredge data object at the specified depths.
    *
    * @pre getDim() == dst.getDim()
    */
   void
   copyDepth2(
      int dst_depth,
      EdgeData<TYPE>& dst,
      int src_depth) const;

   /*!
    * @brief Add data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be OuteredgeData of the same DIM and
    * TYPE and the overlap must be an EdgeOverlap of the same DIM.
    * If not, then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    * @pre dynamic_cast<const OuteredgeData<TYPE> *>(&src) != 0
    */
   virtual void
   sum(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

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
    * @brief Unpack data from stream and add into this patch data object
    * over the specified box overlap region.
    *
    * @pre dynamic_cast<const EdgeOverlap *>(&overlap) != 0
    */
   virtual void
   unpackStreamAndSum(
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
    * @brief Print all outeredge data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point
    *        numbers, and the default is 6 decimal places floats.  For other
    *        types, this value is ignored.
    *
    * @pre getDim() == box.getDim()
    */
   void
   print(
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outeredge data values at the given array depth in
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
    *        decimal places for double and complex floating point
    *        numbers, and the default is 6 decimal places floats.  For other
    *        types, this value is ignored.
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
    * @brief Print all outeredge centered data values for specified axis,
    * face_normal, and side residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point
    *        numbers, and the default is 6 decimal places floats.  For other
    *        types, this value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisSide(
      int axis,
      int face_normal,
      int side,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outeredge centered data values for specified axis,
    * face_normal, side, and depth residing in the specified box.
    *
    * @param axis  integer edge data coordinate axis,
    *              must satisfy 0 <= axis < DIM
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to edge index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point
    *        numbers, and the default is 6 decimal places floats.  For other
    *        types, this value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisSide(
      int axis,
      int face_normal,
      int side,
      const hier::Box& box,
      int depth,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Check that class version and restart file version are equal.
    * If so, read data members from the restart database.
    *
    * @pre restart_db
    */
   virtual void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * @brief Write out the class version number and other data members to
    * the restart database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_OUTEREDGEDATA_VERSION;

   // Unimplemented copy constructor
   OuteredgeData(
      const OuteredgeData&);

   // Unimplemented assignment operator
   OuteredgeData&
   operator = (
      const OuteredgeData&);

   //@{
   //! @name Internal implementations for data copy operations.
   void
   copyFromEdge(
      const EdgeData<TYPE>& src);
   void
   copyFromEdge(
      const EdgeData<TYPE>& src,
      const EdgeOverlap& overlap);
   void
   copyToEdge(
      EdgeData<TYPE>& dst) const;
   void
   copyToEdge(
      EdgeData<TYPE>& dst,
      const EdgeOverlap& overlap) const;
   void
   copyFromOuteredge(
      const OuteredgeData<TYPE>& src);
   void
   copyFromOuteredge(
      const OuteredgeData<TYPE>& src,
      const EdgeOverlap& overlap);
   void
   copyToOuteredge(
      OuteredgeData<TYPE>& dst) const;
   void
   copyToOuteredge(
      OuteredgeData<TYPE>& dst,
      const EdgeOverlap& overlap) const;
   //@}

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> >
   d_data[SAMRAI::MAX_DIM_VAL][SAMRAI::MAX_DIM_VAL][2];
};

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OuteredgeData<TYPE>::template View<DIM> get_view(OuteredgeData<TYPE>& data,
                                                     Args&&... args);

template <int DIM, typename TYPE, typename... Args>
typename OuteredgeData<TYPE>::template ConstView<DIM> get_const_view(
    const OuteredgeData<TYPE>& data,
    Args&&... args);
#endif

}
}

#include "SAMRAI/pdat/OuteredgeData.cpp"

#endif
