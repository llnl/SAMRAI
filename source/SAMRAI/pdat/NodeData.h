/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated node centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_NodeData
#define included_pdat_NodeData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/ArrayView.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class NodeData<TYPE> provides an implementation for data defined
 * at nodes on AMR patches.  It is derived from the hier::PatchData
 * interface common to all SAMRAI patch data types.  Given a CELL-centered
 * AMR index space box, a node data object represents data of some template
 * TYPE and depth at the nodes of the cells in the box.  Here, depth indicates
 * the number of data values at each node index location.  The NodeGeometry
 * class provides the translation between the standard SAMRAI cell-centered
 * AMR index space and node-centered data.
 *
 * A node data array is stored in (i,...,k,d) order, where i,...,k indicates
 * spatial indices and the d indicates the component depth at that locaion.
 * Memory allocation is in column-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.  For example, a
 * three-dimensional node data object created over a CELL-centered
 * AMR index space box [l0:u0,l1:u1,l2:u2] allocates a data array
 * sized as
 * \verbatim
 *
 *   [ l0 : u0+1 ,
 *     l1 : u1+1,
 *     l2 : u2+1 , d ]
 *
 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see ArrayData
 * @see hier::PatchData
 * @see NodeDataFactory
 * @see NodeIndex
 * @see NodeIterator
 * @see NodeGeometry
 */

template<class TYPE>
class NodeData:public hier::PatchData
{
public:
   /*!
    * @brief  Calculate the amount of memory needed to represent node-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            node data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the node data will be allocated.
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
    * @brief The constructor for a node data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            node data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the node data will be allocated.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   NodeData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts);

   /*!
    * @brief The constructor for an node data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            node data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *               of the ghost cell region around the box over which
    *               the node data will be allocated.
    * @param allocator A ResourceAllocator to manage the allocation of the
    *                  underlying data.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   NodeData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts,
      tbox::ResourceAllocator allocator);



   /*!
    * @brief The virtual destructor for a node data object.
    */
   virtual ~NodeData();

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular depth
    * component of the node centered array.
    *
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE *
   getPointer(
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular depth
    * component of the node centered array.
    *
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE *
   getPointer(
      int depth = 0) const;

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
    * @brief Return a reference to the data entry corresponding
    * to a given node index and depth.
    *
    * @pre getDim() == i.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    */
   TYPE&
   operator () (
      const NodeIndex& i,
      int depth = 0);

   /*!
    * @brief Return a const reference to the data entry corresponding
    * to a given node index and depth.
    *
    * @pre getDim() == i.getDim()
    * @pre (depth >= 0) && (depth < getDepth())
    */
   const TYPE&
   operator () (
      const NodeIndex& i,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object of
    * the node centered array.
    */
   ArrayData<TYPE>&
   getArrayData();

   /*!
    * @brief Return a const reference to the array data object of
    * the node centered array.
    */
   const ArrayData<TYPE>&
   getArrayData() const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * a NodeData of the same DIM and TYPE.  If not, then an unrecoverable
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
    * a NodeData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<NodeData<TYPE> *>(&dst) != 0
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be NodeData of the same DIM and TYPE
    * and the overlap must be a NodeOverlap of the same DIM.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == src.getDim()
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   virtual void
   copyFuseable(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap)
   {
       d_data->startKernelFuser();
       copy(src, overlap);
       d_data->stopKernelFuser();
   }

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be NodeData of the same DIM and TYPE
    * and the overlap must be a NodeOverlap of the same DIM.  If not,
    * then an unrecoverable error results.
    *
    * @pre getDim() == dst.getDim()
    * @pre dynamic_cast<NodeData<TYPE> *>(&dst) != 0
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
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
      const NodeData<TYPE>& src,
      const hier::Box& box);

   /*!
    * @brief Fast copy (i.e., source and this node data objects are
    * defined over the same box) from the given node source data object to
    * this destination node data object at the specified depths.
    *
    * @pre getDim() == src.getDim()
    */
   void
   copyDepth(
      int dst_depth,
      const NodeData<TYPE>& src,
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
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be a
    * NodeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   virtual void
   packStreamFuseable(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const
   {
      d_data->startKernelFuser();
      packStream(stream, overlap);
      d_data->stopKernelFuser();
   }

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be a
    * NodeOverlap of the same DIM.
    *
    * @pre dynamic_cast<const NodeOverlap *>(&overlap) != 0
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   virtual void
   unpackStreamFuseable(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap)
   {
      d_data->startKernelFuser();
      unpackStream(stream, overlap);
      d_data->stopKernelFuser();
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
    * @brief Print all node data values residing in the specified box.
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
    * @brief Print all node data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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

#ifdef SAMRAI_HAVE_CONDUIT
   /*!
    * @brief Put data into a conduit node for the blueprint format
    *
    * @param domain_node   Node holding blueprint data for one patch
    * @param field_name    Name of this field
    * @param topology_name Identifier of the topology for this field
    * @param depth         Use the data at this depth
    */
   void
   putBlueprintField(
      conduit::Node& domain_node,
      const std::string& field_name,
      const std::string& topology_name,
      int depth = 0);
#endif

   /*!
    * The node iterator iterates over the elements of a node
    * centered box geometry.  This typedef is a convenience for
    * using the NodeIterator class.
    */
   typedef NodeIterator iterator;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_NODEDATA_VERSION;

   // Unimplemented copy constructor
   NodeData(
      const NodeData&);

   // Unimplemented assignment operator
   NodeData&
   operator = (
      const NodeData&);

   void
   copyWithRotation(
      const NodeData<TYPE>& src,
      const NodeOverlap& overlap);

   void
   packWithRotation(
      tbox::MessageStream& stream,
      const NodeOverlap& overlap) const;

   int d_depth;

   std::shared_ptr<ArrayData<TYPE> > d_data;

};

#if defined(HAVE_RAJA)
template<int DIM, typename TYPE, typename... Args>
typename NodeData<TYPE>::template View<DIM> get_view(NodeData<TYPE>& data, Args&&... args);

template<int DIM, typename TYPE, typename... Args>
typename NodeData<TYPE>::template ConstView<DIM> get_const_view(const NodeData<TYPE>& data, Args&&... args);
#endif

}
}

#include "SAMRAI/pdat/NodeData.cpp"

#endif
