/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outerside centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OutersideData
#define included_pdat_OutersideData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OutersideData<TYPE> provides an implementation for data defined
 * at cell sides (faces) on the boundaries of AMR patches.  It is derived from
 * the hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outerside data object represents
 * data of some template TYPE and depth on the cell sides (faces) on the
 * boundary of the box.  Here, depth indicates the number of data values at
 * each face index location.  The OuteredgsideGeometry class provides the
 * translation between the standard SAMRAI cell-centered AMR index space and
 * outerside-centered data.
 *
 * Outerside data is stored in 2*DIM arrays, each of which contains data
 * associated with side (face) indices normal to a coordinate axis direction
 * and an upper or lower box side (face) in the face normal direction.
 * The data layout in the outerside data arrays matches the corresponding array
 * sections provided by the side data implementation.  Also, in each of array,
 * memory allocation is in column-major ordering (e.g., Fortran style) so that
 * the leftmost index runs fastest in memory.  For example, a three-dimensional
 * outerside data object created over a CELL-centered AMR index space
 * [l0:u0,l1:u1,l2:u2] allocates six data arrays sized as follows:
 * \verbatim
 *
 * face normal 0:
 *   lower face      [ l0:l0     , l1:u1     , l2:u2     , d ]
 *   upper face      [ u0+1:u0+1 , l1:u1     , l2:u2     , d ]
 *
 * face normal 1:
 *   lower face      [ l0:u0     , l1:l1     , l2:u2     , d ]
 *   upper face      [ l0:u0     , u1+1:u1+1 , l2:u2     , d ]
 *
 * face normal 2:
 *   lower face      [ l0:u0     , l1:u1     , l2:l2     , d ]
 *   upper face      [ l0:u0     , l1:u1     , u2+1:u2+1 , d ]
 *
 * \endverbatim
 * Here the face normal directions 0, 1, 2 can be thought of as the x, y, and z
 * face normal directions, respectively, and d is the depth index (i.e., number
 * of values at each face index location).  Other spatial dimensions are
 * represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * IMPORTANT: The OuterfaceData<TYPE> class provides the same storage
 * as this outerside data class, except that the coordinate directions of the
 * individual arrays are permuted; i.e., OuterfaceData is consistent
 * with the FaceData implementation.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see OutersideDataFactory
 * @see OutersideGeometry
 * @see SideIterator
 * @see SideIndex
 */

template<class TYPE>
class OutersideData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outerside-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerside data object will be created.
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
    * @brief Constructor for an outerside data object.
    *
    * Note: Outerside data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerside data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    *
    * @pre depth > 0
    */
   OutersideData(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Constructor for an outerside data object.
    *
    * Note: Outerside data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerside data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data.
    *
    * @pre depth > 0
    */
   OutersideData(
      const hier::Box& box,
      int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Virtual destructor for a outerside data object.
    */
   virtual ~OutersideData();

   /*!
    * @brief Return the depth (i.e., the number of data values for
    * each spatial location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular
    * side normal, side, and depth component of the outerside centered
    * array.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param depth integer depth component
    *
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE *
   getPointer(
      int side_normal,
      int side,
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular
    * side normal, side, and depth component of the outerside centered
    * array.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param depth integer depth component
    *
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE *
   getPointer(
      int side_normal,
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
   View<DIM> getView(int side_normal, int side, int depth = 0);

   /*!
    * @brief Get an ArrayView that can access the array for RAJA looping.
    */
   template <int DIM>
   ConstView<DIM> getConstView(int side_normal, int side, int depth = 0) const;
#endif

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given side index, side location, and depth.
    *
    * @param i const reference to SideIndex, @em MUST be
    *          an index on the outerside of the box.
    * @param side  integer (lower/upper location of outerside data)
    * @param depth integer depth component
    *
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE&
   operator () (
      const SideIndex& i,
      int side,
      int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given side index, side location, and depth.
    *
    * @param i const reference to SideIndex, @em MUST be
    *          an index on the outerside of the box.
    * @param side  integer (lower/upper location of outerside data)
    * @param depth integer depth component
    *
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE&
   operator () (
      const SideIndex& i,
      int side,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * side normal and side location of the outerside centered array.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    *
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   ArrayData<TYPE>&
   getArrayData(
      int side_normal,
      int side);

   /*!
    * @brief Return a const reference to the array data object for
    * side normal and side location of the outerside centered array.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    *
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   const ArrayData<TYPE>&
   getArrayData(
      int side_normal,
      int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * SideData the same DIM and TYPE.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const SideData<TYPE> *>(&src) != 0
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
    * SideData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<SideData<TYPE> *>(&dst) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * IMPORTANT: this routine is @b not @b yet @b implemented!
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be SideData of the same DIM
    * and TYPE and the overlap must be a SideOverlap of the same
    * DIM.  If not, then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<SideData<TYPE> *>(&dst) != 0
    * @pre dynamic_cast<const SideOverlap *>(&overlap) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes side and outerside data objects are
    * defined over the same box) from the given side source data object to
    * this destination outerside data object at the specified depths.
    *
    * @pre getDim == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const SideData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes side and outerside data objects are
    * defined over the same box) to the given side destination data object
    * from this source outerside data object at the specified depths.
    *
    * @pre getDim == dst.getDim()
    */
   void
   copyDepth2(
      int dst_depth,
      SideData<TYPE>& dst,
      int src_depth) const;

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
    * @pre dynamic_cast<const SideOverlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * SideOverlap of the same DIM.
    *
    * @pre dynamic_cast<const SideOverlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * SideOverlap of the same DIM.
    *
    * @pre dynamic_cast<const SideOverlap *>(&overlap) != 0
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
    * @brief Print all outerside data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    */
   void
   print(
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outerside data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param depth integer depth component
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
    * @brief Print all outerside centered data values for specified
    * side_normal and side location residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisSide(
      tbox::Dimension::dir_t side_normal,
      int side,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outerside centered data values for specified
    * side_normal, side location, and depth residing in the specified box.
    *
    * @param side_normal  integer side normal direction for data
    * @param side integer lower (0) or upper (1) side of outerside
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to side index space.
    * @param depth integer depth component
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    * @pre (side_normal >= 0) && (side_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisSide(
      tbox::Dimension::dir_t side_normal,
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
   static const int PDAT_OUTERSIDEDATA_VERSION;

   // Unimplemented copy constructor
   OutersideData(
      const OutersideData&);

   // Unimplemented assignment operator
   OutersideData&
   operator = (
      const OutersideData&);

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> > d_data[SAMRAI::MAX_DIM_VAL][2];
};

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OutersideData<TYPE>::template View<DIM> get_view(OutersideData<TYPE>& data,
                                                     Args&&... args);

template <int DIM, typename TYPE, typename... Args>
typename OutersideData<TYPE>::template ConstView<DIM> get_const_view(
    const OutersideData<TYPE>& data,
    Args&&... args);
#endif


}
}

#include "SAMRAI/pdat/OutersideData.cpp"

#endif
