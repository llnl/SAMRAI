/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outernode centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeData
#define included_pdat_OuternodeData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class OuternodeData<TYPE> provides an implementation for data defined
 * at cell nodes on the boundaries of AMR patches.  It is derived from the
 * hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outernode data object represents
 * data of some template TYPE and depth on the cell nodes on the boundary
 * of the box.  Here, depth indicates the number of data values at each node
 * index location.  The OuternodnodeGeometry class provides the translation
 * between the standard SAMRAI cell-centered AMR index space and
 * outernode-centered data.
 *
 * Outernode data is stored in 2*DIM arrays, each of which contains data
 * associated with node indices on an upper or lower box face in some
 * corrdinate direction.  The data layout in the outernode data arrays matches
 * the corresponding array sections provided by the node data implementation.
 * Where a node index falls on more than one box face (patch boundary edges and
 * corners), the outernode data value belongs to only one data array so that
 * there are no redundant data values.  Specifically, when DIM > 1, outernode
 * data boxes are "trimmed" so that each node index that lives on more than one
 * face on the box boundary will be associated with the face of the largest
 * coordinate direction and only that face.  Within each array, data is stored
 * in (i,...,k,d) order, where i,...,k indicates a spatial index and the d
 * indicates the component depth at that location.  Memory allocation is
 * in column-major ordering (e.g., Fortran style) so that the leftmost
 * index runs fastest in memory.
 *
 * To illustrate the outernode data layout, in particular the "box trimming"
 * that prevents redundant data values, we describe the data for a
 * three-dimensional outernode data object instantiated over a box
 * [l0:u0,l1:u1,l2:u2] in the standard SAMRAI cell-centered AMR index space.
 *
 * \verbatim
 *
 *    Here face normal directions 0, 1, and 2 can be thought of as X, Y, Z
 *    respectively, and d is the data depth.
 *
 *    face normal 0:
 *        lower    [ l0 : l0     , l1+1 : u1   , l2+1 : u2 , d ]
 *        upper    [ u0+1 : u0+1 , l1+1 : u1   , l2+1 : u2 , d ]
 *        Note: Boxes are trimmed at edges intersecting faces with
 *              normal directions 1 and 2 so that node indices shared
 *              with those faces appear in data arrays associated with
 *              higher direction faces.
 *
 *    face normal 1:
 *        lower    [ l0 : u0+1   , l1 : l1     , l2+1 : u2 , d ]
 *        upper    [ l0 : u0+1   , u1+1 : u1+1 , l2+1 : u2 , d ]
 *        Note: Boxes are trimmed at edges intersecting faces with
 *              normal direction 2 so that node indices shared
 *              with those faces appear in data arrays associated with
 *              higher direction faces.
 *
 *    face normal 2:
 *        lower    [ l0 : u0+1   , l1 : u1+1   , l2 : l2     , d ]
 *        upper    [ l0 : u0+1   , l1 : u1+1   , u2+1 : u2+1 , d ]
 *        Note: Boxes are not trimmed.
 *
 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see OuternodeDataFactory
 * @see OuternodeGeometry
 * @see NodeIterator
 * @see NodeIndex
 */

template<class TYPE>
class OuternodeData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outernode-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outernode data object will be created.
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
    * @brief Constructor for an outernode data object.
    *
    * Note: Outernode data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outernode data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    *
    * @pre depth > 0
    */
   OuternodeData(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Constructor for an outernode data object.
    *
    * Note: Outernode data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outernode data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data. 
    *
    * @pre depth > 0
    */
   OuternodeData(
      const hier::Box& box,
      int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Virtual destructor for a outernode data object.
    */
   virtual ~OuternodeData();

   /*!
    * @brief Return the depth (e.g., the number of components at each spatial
    * location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Returns true if outernode data exists for the given
    * face normal direction; false otherwise.
    *
    * @param face_normal  integer face normal direction for data
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    */
   bool
   dataExists(
      int face_normal) const;

   /*!
    * @brief Return the box of valid node indices for
    *        outernode data.  Note: the returned box
    *        will reside in the @em node index space.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   hier::Box
   getDataBox(
      tbox::Dimension::dir_t face_normal,
      int side);

   /*!
    * @brief Get a pointer to the beginning of a particular face normal,
    * side, and depth component of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outernode
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
    * @brief Get a const pointer to the beginning of a particular face
    * normal, side, and depth component of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outernode
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
    * to a given node index and depth.
    *
    * @param i const reference to NodeIndex, @em MUST be
    *          an index on the outernode of the box.
    * @param depth integer depth component
    *
    * @pre getDim() == i.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE&
   operator () (
      const NodeIndex& i,
      int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given node index and depth.
    *
    * @param i const reference to NodeIndex, @em MUST be
    *          an index on the outernode of the box.
    * @param depth integer depth component
    *
    * @pre getDim() == i.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE&
   operator () (
      const NodeIndex& i,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * face normal, and side index of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue()));
    * @pre (side == 0) || (side == 1)
    */
   ArrayData<TYPE>&
   getArrayData(
      int face_normal,
      int side);

   /*!
    * @brief Return a const reference to the array data object for
    * face normal, and side index of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    *
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue()));
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
    * either NodeData or OuternodeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    * @pre (dynamic_cast<const NodeData<TYPE> *>(&src) != 0) ||
    *      (dynamic_cast<const OuternodeData<TYPE> *>(&src) != 0)
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
    * either NodeData or OuternodeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre (dynamic_cast<NodeData<TYPE> *>(&dst) != 0) ||
    *      (dynamic_cast<OuternodeData<TYPE> *>(&dst) != 0)
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be either NodeData or OuternodeData
    * of the same DIM and TYPE and the overlap must be an NodeOverlap
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    * @pre (dynamic_cast<const NodeData<TYPE> *>(&src) != 0) ||
    *      (dynamic_cast<const OuternodeData<TYPE> *>(&src) != 0)
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be either NodeData or OuternodeData
    * of the same DIM and TYPE and the overlap must be an NodeOverlap
    * of the same DIM.  If not, then an unrecoverable error
    * results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    * @pre (dynamic_cast<NodeData<TYPE> *>(&dst) != 0) ||
    *      (dynamic_cast<OuternodeData<TYPE> *>(&dst) != 0)
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes node and outernode data objects are
    * defined over the same box) from the given node source data object to
    * this destination outernode data object at the specified depths.
    *
    * @pre getDim() == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const NodeData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes node and outernode data objects are
    * defined over the same box) to the given node destination data object
    * from this source outernode data object at the specified depths.
    *
    * @pre getDim() == dst.getDim()
    */
   void
   copyDepth2(
      int dst_depth,
      NodeData<TYPE>& dst,
      int src_depth) const;

   /*!
    * @brief Add data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be OuternodeData of the same DIM and
    * TYPE and the overlap must be an EdgeOverlap of the same DIM.
    * If not, then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    * @pre dynamic_cast<const OuternodeData<TYPE> *>(&src) != 0
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
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * NodeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region.  The overlap must be an
    * NodeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Unpack data from stream and add into this patch data object
    * over the specified box overlap region.  The overlap must be an
    * NodeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
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
    * @brief Print all outernode data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
    * @brief Print all outernode data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
    * @brief Print all outernode centered data values for specified
    * face_normal and side residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
   printAxisSide(
      int face_normal,
      int side,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all outernode centered data values for specified
    * face_normal, side, and depth residing in the specified box.
    *
    * @param face_normal  integer face normal direction for data
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
    * @pre (face_normal >= 0) && (face_normal < getDim().getValue())
    * @pre (side == 0) || (side == 1)
    */
   void
   printAxisSide(
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
   static const int PDAT_OUTERNODEDATA_VERSION;

   // Unimplemented copy constructor
   OuternodeData(
      const OuternodeData&);

   // Unimplemented assignment operator
   OuternodeData&
   operator = (
      const OuternodeData&);

   //@
   //! @name Internal implementations of data copy operations.
   void
   copyFromNode(
      const NodeData<TYPE>& src);
   void
   copyFromNode(
      const NodeData<TYPE>& src,
      const NodeOverlap& overlap);
   void
   copyToNode(
      NodeData<TYPE>& dst) const;
   void
   copyToNode(
      NodeData<TYPE>& dst,
      const NodeOverlap& overlap) const;
   void
   copyFromOuternode(
      const OuternodeData<TYPE>& src);
   void
   copyFromOuternode(
      const OuternodeData<TYPE>& src,
      const NodeOverlap& overlap);
   void
   copyToOuternode(
      OuternodeData<TYPE>& dst) const;
   void
   copyToOuternode(
      OuternodeData<TYPE>& dst,
      const NodeOverlap& overlap) const;
   //@

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> > d_data[SAMRAI::MAX_DIM_VAL][2];
};


#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OuternodeData<TYPE>::template View<DIM> get_view(OuternodeData<TYPE>& data,
                                                     Args&&... args);

template <int DIM, typename TYPE, typename... Args>
typename OuternodeData<TYPE>::template ConstView<DIM> get_const_view(
    const OuternodeData<TYPE>& data,
    Args&&... args);
#endif


}
}

#include "SAMRAI/pdat/OuternodeData.cpp"

#endif
