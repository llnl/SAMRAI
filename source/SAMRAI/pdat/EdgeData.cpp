/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated edge centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_EdgeData_C
#define included_pdat_EdgeData_C

#include "SAMRAI/pdat/EdgeData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int EdgeData<TYPE>::PDAT_EDGEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for edge data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
EdgeData<TYPE>::EdgeData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);

   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::Box edge_box =
         EdgeGeometry::toEdgeBox(getGhostBox(), d);
      d_data[d].reset(new ArrayData<TYPE>(edge_box, depth));
   }
}

template<class TYPE>
EdgeData<TYPE>::EdgeData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts,
   tbox::ResourceAllocator allocator):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);

   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::Box edge_box =
         EdgeGeometry::toEdgeBox(getGhostBox(), d);
      d_data[d].reset(new ArrayData<TYPE>(edge_box, depth, allocator));
   }
}

template<class TYPE>
EdgeData<TYPE>::~EdgeData()
{
}

template<class TYPE>
int
EdgeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE *
EdgeData<TYPE>::getPointer(
   int axis,
   int depth)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis]->getPointer(depth);
}

template<class TYPE>
const TYPE *
EdgeData<TYPE>::getPointer(
   int axis,
   int depth) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis]->getPointer(depth);
}

#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename EdgeData<TYPE>::template View<DIM> EdgeData<TYPE>::getView(
    int axis,
    int depth)
{
  const hier::Box edge_box =
      EdgeGeometry::toEdgeBox(getGhostBox(), axis);
  return EdgeData<TYPE>::View<DIM>(getPointer(axis, depth), edge_box);
}

template <class TYPE>
template <int DIM>
typename EdgeData<TYPE>::template ConstView<DIM> EdgeData<TYPE>::getConstView(
    int axis,
    int depth) const
{
  const hier::Box edge_box =
      EdgeGeometry::toEdgeBox(getGhostBox(), axis);
  return EdgeData<TYPE>::ConstView<DIM>(getPointer(axis, depth),
                                        edge_box);
}
#endif

template<class TYPE>
TYPE&
EdgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*(d_data[axis]))(i, depth);
}

template<class TYPE>
const TYPE&
EdgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*(d_data[axis]))(i, depth);
}

template<class TYPE>
ArrayData<TYPE>&
EdgeData<TYPE>::getArrayData(
   int axis)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   return *(d_data[axis]);
}

template<class TYPE>
const ArrayData<TYPE>&
EdgeData<TYPE>::getArrayData(
   int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   return *(d_data[axis]);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two edge centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const EdgeData<TYPE>* t_src = dynamic_cast<const EdgeData<TYPE> *>(&src);
   if (t_src == 0) {
      src.copy2(*this);
   } else {
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::Box box = d_data[d]->getBox() * t_src->d_data[d]->getBox();
         if (!box.empty()) {
            d_data[d]->copy(*(t_src->d_data[d]), box);
         }
      }
   }
}

template<class TYPE>
void
EdgeData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   EdgeData<TYPE>* t_dst = CPP_CAST<EdgeData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != 0);

   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::Box box = d_data[d]->getBox() * t_dst->d_data[d]->getBox();
      if (!box.empty()) {
         t_dst->d_data[d]->copy(*(d_data[d]), box);
      }
   }
}

/*
 *************************************************************************
 *
 * Copy data from the source into the destination according to the
 * overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const EdgeData<TYPE>* t_src = dynamic_cast<const EdgeData<TYPE> *>(&src);

   const EdgeOverlap* t_overlap = dynamic_cast<const EdgeOverlap *>(&overlap);

   if ((t_src == 0) || (t_overlap == 0)) {
      src.copy2(*this, overlap);
   } else {
      if (t_overlap->getTransformation().getRotation() ==
          hier::Transformation::NO_ROTATE) {

         const hier::Transformation& transformation =
            t_overlap->getTransformation();
         for (int d = 0; d < getDim().getValue(); ++d) {
            const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
            d_data[d]->copy(*(t_src->d_data[d]), box_list, transformation);
         }
      } else {
         copyWithRotation(*t_src, *t_overlap);
      }
   }
}

template<class TYPE>
void
EdgeData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   EdgeData<TYPE>* t_dst = CPP_CAST<EdgeData<TYPE> *>(&dst);
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != 0);
   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      const hier::IntVector& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
         t_dst->d_data[d]->copy(*(d_data[d]), box_list, src_offset);
      }
   } else {
      t_dst->copyWithRotation(*this, *t_overlap);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::copyOnBox(
   const EdgeData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      const hier::Box edge_box = EdgeGeometry::toEdgeBox(box, axis);
      d_data[axis]->copy(src.getArrayData(axis), edge_box);
   }

}

template<class TYPE>
void
EdgeData<TYPE>::copyWithRotation(
   const EdgeData<TYPE>& src,
   const EdgeOverlap& overlap)
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(src.getDim());
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Box rotatebox(src.getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   for (int i = 0; i < dim.getValue(); ++i) {
      const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

      hier::Box edge_rotatebox(EdgeGeometry::toEdgeBox(rotatebox, i));

      for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
           bi != overlap_boxes.end(); ++bi) {
         const hier::Box& overlap_box = *bi;

         const hier::Box copybox(edge_rotatebox * overlap_box);

         if (!copybox.empty()) {
            const int depth = ((getDepth() < src.getDepth()) ?
                               getDepth() : src.getDepth());

            hier::Box::iterator ciend(copybox.end());
            for (hier::Box::iterator ci(copybox.begin()); ci != ciend; ++ci) {

               EdgeIndex dst_index(*ci, 0, 0);
               dst_index.setAxis(i);
               EdgeIndex src_index(dst_index);
               EdgeGeometry::transform(src_index, back_trans);

               for (int d = 0; d < depth; ++d) {
                  (*this)(dst_index, d) = src(src_index, d);
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from an edge data object to this edge data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::copyDepth(
   int dst_depth,
   const EdgeData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::Box box = d_data[d]->getBox() * src.d_data[d]->getBox();
      if (!box.empty()) {
         d_data[d]->copyDepth(dst_depth, *(src.d_data[d]), src_depth, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the buffer space needed to pack/unpack messages on the box
 * region using the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
bool
EdgeData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
EdgeData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   size_t size = 0;
   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      size += d_data[d]->getDataStreamSize(
            t_overlap->getDestinationBoxContainer(d),
            offset);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Pack/unpack data into/out of the message streams using the index
 * space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      const hier::Transformation& transformation = t_overlap->getTransformation();
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
         if (!boxes.empty()) {
            d_data[d]->packStream(stream, boxes, transformation);
         }
      }
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const EdgeOverlap& overlap) const
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(getDim());
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Box rotatebox(getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   const int depth = getDepth();

   for (int i = 0; i < dim.getValue(); ++i) {
      const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

      const size_t size = depth * overlap_boxes.getTotalSizeOfBoxes();
      std::vector<TYPE> buffer(size);

      hier::Box edge_rotatebox(EdgeGeometry::toEdgeBox(rotatebox, i));

      int buf_count = 0;
      for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
           bi != overlap_boxes.end(); ++bi) {
         const hier::Box& overlap_box = *bi;

         const hier::Box copybox(edge_rotatebox * overlap_box);

         if (!copybox.empty()) {

            for (int d = 0; d < depth; ++d) {

               hier::Box::iterator ciend(copybox.end());
               for (hier::Box::iterator ci(copybox.begin());
                    ci != ciend; ++ci) {

                  EdgeIndex src_index(*ci, 0, 0);
                  src_index.setAxis(i);
                  EdgeGeometry::transform(src_index, back_trans);

                  buffer[buf_count] = (*this)(src_index, d);
                  ++buf_count;
               }
            }
         }
      }
      stream.pack(&buffer[0], size);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      if (!boxes.empty()) {
         d_data[d]->unpackStream(stream, boxes, offset);
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  edge centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
EdgeData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   for (int d = 0; d < box.getDim().getValue(); ++d) {
      const hier::Box edge_box = EdgeGeometry::toEdgeBox(ghost_box, d);
      size += ArrayData<TYPE>::getSizeOfData(edge_box, depth);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Fill the edge centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fill(t, d);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fill(t, EdgeGeometry::toEdgeBox(box, i), d);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fillAll(t);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fillAll(t, EdgeGeometry::toEdgeBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print edge centered data.  Note:  makes call to specialized printAxis
 * routine in EdgeDataSpecialized.cpp
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {
      os << "Array axis = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      os << "Array axis = " << axis << std::endl;
      printAxis(axis, box, depth, os, prec);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::printAxis(
   int axis,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<class TYPE>
void
EdgeData<TYPE>::printAxis(
   int axis,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   os.precision(prec);
   EdgeIterator iend(EdgeGeometry::end(box, axis));
   for (EdgeIterator i(EdgeGeometry::begin(box, axis)); i != iend; ++i) {
      os << "array" << *i << " = " << (*(d_data[axis]))(*i, depth)
         << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in the d_depth data member to the database.  Then tell
 * d_data to read itself in from the restart database.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_EDGEDATA_VERSION");
   if (ver != PDAT_EDGEDATA_VERSION) {
      TBOX_ERROR("EdgeData<getDim()>::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");

   std::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); ++i) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = restart_db->getDatabase(array_name);
      d_data[i]->getFromRestart(array_database);
   }
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth data member to the
 * database.  Then tells d_data to write itself to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
EdgeData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_EDGEDATA_VERSION", PDAT_EDGEDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);

   std::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); ++i) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i);
      array_database = restart_db->putDatabase(array_name);
      d_data[i]->putToRestart(array_database);
   }
}


#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename EdgeData<TYPE>::template View<DIM> get_view(EdgeData<TYPE>& data,
                                                     Args&&... args)
{
  return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename EdgeData<TYPE>::template ConstView<DIM> get_const_view(
    const EdgeData<TYPE>& data,
    Args&&... args)
{
  return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#endif
