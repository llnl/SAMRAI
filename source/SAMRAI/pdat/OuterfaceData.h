/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outerface centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceData
#define included_pdat_OuterfaceData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OuterfaceData<TYPE> provides an implementation for data defined
 * at cell faces on the boundaries of AMR patches.  It is derived from the
 * hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outerface data object represents
 * data of some template TYPE and depth on the cell faces on the boundary
 * of the box.  Here, depth indicates the number of data values at each face
 * index location.  The OuteredgfaceGeometry class provides the translation
 * between the standard SAMRAI cell-centered AMR index space and
 * outerface-centered data.
 *
 * Outerface data is stored in 2*DIM arrays, each of which contains data
 * associated with face indices normal to a coordinate axis direction and an
 * upper or lower box face in the face normal direction.  The data layout in
 * the outerface data arrays matches the corresponding array sections provided
 * by the face data implementation.  Specifically, within each array, the data
 * indices are cyclically permuted to be consistent with the FaceData<TYPE>
 * implementation.  Also, in each of array, memory allocation is in
 * column-major ordering (e.g., Fortran style) so that the leftmost index runs
 * fastest in memory.  For example, a three-dimensional outerface data object
 * created over a CELL-centered AMR index space [l0:u0,l1:u1,l2:u2] allocates
 * six data arrays sized as follows:
 * \verbatim
 *
 * face normal 0:
 *   lower face      [ l0:l0     , l1:u1 , l2:u2 , d ]
 *   upper face      [ u0+1:u0+1 , l1:u1 , l2:u2 , d ]
 *
 * face normal 1:
 *   lower face      [ l1:l1     , l2:u2 , l0:u0 , d ]
 *   upper face      [ u1+1:u1+1 , l2:u2 , l0:u0 , d ]
 *
 * face normal 2:
 *   lower face      [ l2:l2     , l0:u0 , l1:u1 , d ]
 *   upper face      [ u2+1:u2+1 , l0:u0 , l1:u1 , d ]
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
 * IMPORTANT: The OutersideData<TYPE> class provides the same storage
 * as this outerface data class, except that the coordinate directions of the
 * individual arrays are not permuted; i.e., OutersideData is consistent
 * with the SideData implementation.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see OuterfaceDataFactory
 * @see OuterfaceGeometry
 * @see FaceIterator
 * @see FaceIndex
 */

template<class TYPE>
class OuterfaceData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outerface-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerface data object will be created.
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
    * @brief Constructor for an outerface data object.
    *
    * Note: Outerface data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerface data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    *
    * @pre depth > 0
    */
   OuterfaceData(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Constructor for an outerface data object.
    *
    * Note: Outerface data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outerface data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data. 
    *
    * @pre depth > 0
    */
   OuterfaceData(
      const hier::Box& box,
      int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Virtual destructor for a outerface data object.
    */
   virtual ~OuterfaceData();

   /*!
    * @brief Return the depth (i.e., the number of data values for
    * each spatial location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular
    * face normal, side, and depth component of the outerface centered
    * array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    * @param depth integer depth component
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE *
   getPointer(
      int face_normal,
      int side,
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular
    * face normal, side location, and depth component of the outerface
    * centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    * @param depth integer depth component
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE *
   getPointer(
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
   View<DIM> getView(int face_normal, int side, int depth = 0);
 
   /*!
    * @brief Get a const ArrayView that can access the array for RAJA looping.
    */
   template <int DIM>
   ConstView<DIM> getConstView(int face_normal, int side, int depth = 0) const;
#endif

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given face index, side location, and depth.
    *
    * @param i const reference to FaceIndex, @em MUST be
    *          an index on the outerface of the box.
    * @param side  integer (lower/upper location of outerface data)
    * @param depth integer depth component
    *
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE&
   operator () (
      const FaceIndex& i,
      int side,
      int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given face index, side location, and depth.
    *
    * @param i const reference to FaceIndex, @em MUST be
    *          an index on the outerface of the box.
    * @param side  integer (lower/upper location of outerface data)
    * @param depth integer depth component
    *
    * @pre (i.getAxis() >= 0) && (i.getAxis() < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE&
   operator () (
      const FaceIndex& i,
      int side,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * face normal and side location of the outerface centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   ArrayData<TYPE>&
   getArrayData(
      int face_normal,
      int side);

   /*!
    * @brief Return a const reference to the array data object for
    * face normal and side location of the outerface centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   const ArrayData<TYPE>&
   getArrayData(
      int face_normal,
      int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * FaceData the same DIM and TYPE.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const FaceData<TYPE> *>(&src) != 0
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
    * FaceData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<FaceData<TYPE> *>(&dst) != 0
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
    * Currently, destination data must be FaceData of the same DIM
    * and TYPE and the overlap must be a FaceOverlap of the same
    * DIM.  If not, then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<FaceData<TYPE> *>(&dst) != 0
    * @pre dynamic_cast<const FaceOverlap *>(&overlap) != 0
    * @pre overlap.getTransformation().getRotation() == hier::Transformation::NO_ROTATE
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes face and outerface data objects are
    * defined over the same box) from the given face source data object to
    * this destination outerface data object at the specified depths.
    *
    * @pre getDim() == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const FaceData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes face and outerface data objects are
    * defined over the same box) to the given face destination data object
    * from this source outerface data object at the specified depths.
    *
    * @pre getDim() == dst.getDim()
    */
   void
   copyDepth2(
      int dst_depth,
      FaceData<TYPE>& dst,
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
    * @pre dynamic_cast<const FaceOverlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * FaceOverlap of the same DIM.
    *
    * @pre dynamic_cast<const FaceOverlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap,
      tbox::KernelFuser& fuser) const
   {
      NULL_USE(fuser);
      packStream(stream, overlap);
   }

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * FaceOverlap of the same DIM.
    *
    * @pre dynamic_cast<const FaceOverlap *>(&overlap) != 0
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap,
      tbox::KernelFuser& fuser)
   {
      NULL_USE(fuser);
      unpackStream(stream, overlap);
   }

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
    * @brief Print all outerface data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
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
    * @brief Print all outerface data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
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
    * @brief Print all outerface centered data values for specified
    * face_normal and side residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default is 12
    *        decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    *
    * @pre getDim() == box.getDim()
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisFace(
      tbox::Dimension::dir_t face_normal,
      int side,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outerface centered data values for specified
    * face_normal, side, and depth residing in the specified box.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outerface
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
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
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisFace(
      tbox::Dimension::dir_t face_normal,
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
   static const int PDAT_OUTERFACEDATA_VERSION;

   // Unimplemented copy constructor
   OuterfaceData(
      const OuterfaceData&);

   // Unimplemented assignment operator
   OuterfaceData&
   operator = (
      const OuterfaceData&);

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> > d_data[SAMRAI::MAX_DIM_VAL][2];
};

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OuterfaceData<TYPE>::template View<DIM> get_view(OuterfaceData<TYPE>& data,
                                                     Args&&... args);

template <int DIM, typename TYPE, typename... Args>
typename OuterfaceData<TYPE>::template ConstView<DIM> get_const_view(
    const OuterfaceData<TYPE>& data,
    Args&&... args);
#endif


}
}

#include "SAMRAI/pdat/OuterfaceData.cpp"

#endif
