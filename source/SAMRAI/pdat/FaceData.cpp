/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated face centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_FaceData_C
#define included_pdat_FaceData_C

#include "SAMRAI/pdat/FaceData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int FaceData<TYPE>::PDAT_FACEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for face data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
FaceData<TYPE>::FaceData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::Box face = FaceGeometry::toFaceBox(getGhostBox(), d);
      d_data[d].reset(new ArrayData<TYPE>(face, depth));
   }
}

template <class TYPE>
FaceData<TYPE>::FaceData(const hier::Box& box,
                         int depth,
                         const hier::IntVector& ghosts,
                         tbox::ResourceAllocator allocator)
    : hier::PatchData(box, ghosts), d_depth(depth)
{
  TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
  TBOX_ASSERT(depth > 0);
  TBOX_ASSERT(ghosts.min() >= 0);

  for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
    const hier::Box face = FaceGeometry::toFaceBox(getGhostBox(), d);
    d_data[d].reset(new ArrayData<TYPE>(face, depth, allocator));
  }
}

template<class TYPE>
FaceData<TYPE>::~FaceData()
{
}

template<class TYPE>
int
FaceData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE *
FaceData<TYPE>::getPointer(
   int face_normal,
   int depth)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal]->getPointer(depth);
}

template<class TYPE>
const TYPE *
FaceData<TYPE>::getPointer(
   int face_normal,
   int depth) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal]->getPointer(depth);
}

#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename FaceData<TYPE>::template View<DIM> FaceData<TYPE>::getView(
    int face_normal,
    int depth)
{
  const hier::Box face_box =
      FaceGeometry::toFaceBox(getGhostBox(), face_normal);
  return FaceData<TYPE>::View<DIM>(getPointer(face_normal, depth), face_box);
}

template <class TYPE>
template <int DIM>
typename FaceData<TYPE>::template ConstView<DIM> FaceData<TYPE>::getConstView(
    int face_normal,
    int depth) const
{
  const hier::Box face_box =
      FaceGeometry::toFaceBox(getGhostBox(), face_normal);
  return FaceData<TYPE>::ConstView<DIM>(getPointer(face_normal, depth),
                                        face_box);
}
#endif

template<class TYPE>
TYPE&
FaceData<TYPE>::operator () (
   const FaceIndex& i,
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
FaceData<TYPE>::operator () (
   const FaceIndex& i,
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
FaceData<TYPE>::getArrayData(
   int face_normal)
{

   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return *(d_data[face_normal]);
}

template<class TYPE>
const ArrayData<TYPE>&
FaceData<TYPE>::getArrayData(
   int face_normal) const
{

   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return *(d_data[face_normal]);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two face centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
FaceData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const FaceData<TYPE>* t_src = dynamic_cast<const FaceData<TYPE> *>(&src);

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
FaceData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   FaceData<TYPE>* t_dst = CPP_CAST<FaceData<TYPE> *>(&dst);

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
FaceData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const FaceData<TYPE>* t_src = dynamic_cast<const FaceData<TYPE> *>(&src);
   const FaceOverlap* t_overlap = dynamic_cast<const FaceOverlap *>(&overlap);

   if ((t_src == 0) || (t_overlap == 0)) {
      src.copy2(*this, overlap);
   } else {
      if (t_overlap->getTransformation().getRotation() ==
          hier::Transformation::NO_ROTATE) {
         const hier::IntVector& src_offset = t_overlap->getSourceOffset();
         for (int d = 0; d < getDim().getValue(); ++d) {
            hier::IntVector face_offset(src_offset);
            if (d > 0) {
               for (int i = 0; i < getDim().getValue(); ++i) {
                  face_offset(i) = src_offset((d + i) % getDim().getValue());
               }
            }
            hier::Transformation transform(hier::Transformation::NO_ROTATE,
                                           face_offset,
                                           t_src->getBox().getBlockId(),
                                           getBox().getBlockId());

            const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
            d_data[d]->copy(*(t_src->d_data[d]), box_list, transform);
         }
      } else {
         copyWithRotation(*t_src, *t_overlap);
      }
   }
}

template<class TYPE>
void
FaceData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   FaceData<TYPE>* t_dst = CPP_CAST<FaceData<TYPE> *>(&dst);

   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != 0);
   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      const hier::IntVector& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); ++d) {
         hier::IntVector face_offset(src_offset);
         if (d > 0) {
            for (int i = 0; i < getDim().getValue(); ++i) {
               face_offset(i) = src_offset((d + i) % getDim().getValue());
            }
         }
         hier::Transformation transform(hier::Transformation::NO_ROTATE,
                                        face_offset,
                                        getBox().getBlockId(),
                                        t_dst->getBox().getBlockId());

         const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
         t_dst->d_data[d]->copy(*(d_data[d]), box_list, transform);
      }
   } else {
      t_dst->copyWithRotation(*this, *t_overlap);
   }
}

template<class TYPE>
void
FaceData<TYPE>::copyOnBox(
   const FaceData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);

   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {
      const hier::Box face_box = FaceGeometry::toFaceBox(box, axis);
      d_data[axis]->copy(src.getArrayData(axis), face_box);
   }

}

template<class TYPE>
void
FaceData<TYPE>::copyWithRotation(
   const FaceData<TYPE>& src,
   const FaceOverlap& overlap)
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

   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

      hier::Box face_rotatebox(FaceGeometry::toFaceBox(rotatebox, i));

      for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
           bi != overlap_boxes.end(); ++bi) {
         const hier::Box& overlap_box = *bi;

         const hier::Box copybox(face_rotatebox * overlap_box);

         if (!copybox.empty()) {
            const int depth = ((getDepth() < src.getDepth()) ?
                               getDepth() : src.getDepth());

            hier::Box::iterator ciend(copybox.end());
            for (hier::Box::iterator ci(copybox.begin());
                 ci != ciend; ++ci) {

               FaceIndex dst_index(*ci, 0, 0);
               dst_index.setAxis(i);
               FaceIndex src_index(dst_index);
               FaceGeometry::transform(src_index, back_trans);

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
 * Perform a fast copy between two arrays at the
 * specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
FaceData<TYPE>::copyDepth(
   int dst_depth,
   const FaceData<TYPE>& src,
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
FaceData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
FaceData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   size_t size = 0;
   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); ++i) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      size += d_data[d]->getDataStreamSize(t_overlap->getDestinationBoxContainer(d),
            face_offset);
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
FaceData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::Transformation& transformation =
      t_overlap->getTransformation();
   if (transformation.getRotation() ==
       hier::Transformation::NO_ROTATE) {
      const hier::IntVector& offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); ++d) {
         hier::IntVector face_offset(offset);
         if (d > 0) {
            for (int i = 0; i < getDim().getValue(); ++i) {
               face_offset(i) = offset((d + i) % getDim().getValue());
            }
         }
         hier::Transformation transform(hier::Transformation::NO_ROTATE,
                                        face_offset,
                                        transformation.getBeginBlock(),
                                        transformation.getEndBlock());

         const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
         if (!boxes.empty()) {
            d_data[d]->packStream(stream, boxes, transform);
         }
      }
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
FaceData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const FaceOverlap& overlap) const
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

   for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
      const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

      const size_t size = depth * overlap_boxes.getTotalSizeOfBoxes();
      std::vector<TYPE> buffer(size);

      hier::Box face_rotatebox(FaceGeometry::toFaceBox(rotatebox, i));

      int buf_count = 0;
      for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
           bi != overlap_boxes.end(); ++bi) {
         const hier::Box& overlap_box = *bi;

         const hier::Box copybox(face_rotatebox * overlap_box);

         if (!copybox.empty()) {

            for (int d = 0; d < depth; ++d) {

               hier::Box::iterator ciend(copybox.end());
               for (hier::Box::iterator ci(copybox.begin());
                    ci != ciend; ++ci) {

                  FaceIndex src_index(*ci, 0, 0);
                  src_index.setAxis(i);
                  FaceGeometry::transform(src_index, back_trans);

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
FaceData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); ++i) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }

      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      if (!boxes.empty()) {
         d_data[d]->unpackStream(stream, boxes, face_offset);
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  face centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
FaceData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);

   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   for (tbox::Dimension::dir_t d = 0; d < box.getDim().getValue(); ++d) {
      const hier::Box face_box = FaceGeometry::toFaceBox(ghost_box, d);
      size += ArrayData<TYPE>::getSizeOfData(face_box, depth);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Fill the face centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
FaceData<TYPE>::fill(
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
FaceData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fill(t, FaceGeometry::toFaceBox(box, i), d);
   }
}

template<class TYPE>
void
FaceData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fillAll(t);
   }
}

template<class TYPE>
void
FaceData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (tbox::Dimension::dir_t i = 0; i < getDim().getValue(); ++i) {
      d_data[i]->fillAll(t, FaceGeometry::toFaceBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print face centered data.  Note:  makes call to specialized printAxis
 * routine in FaceDataSpecialized.cpp
 *
 *************************************************************************
 */

template<class TYPE>
void
FaceData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<class TYPE>
void
FaceData<TYPE>::print(
   const hier::Box& box,
   int d,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      os << "Array face normal = " << axis << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<class TYPE>
void
FaceData<TYPE>::printAxis(
   tbox::Dimension::dir_t face_normal,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((face_normal < getDim().getValue()));

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      printAxis(face_normal, box, d, os, prec);
   }
}

template<class TYPE>
void
FaceData<TYPE>::printAxis(
   tbox::Dimension::dir_t face_normal,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((face_normal < getDim().getValue()));

   os.precision(prec);
   FaceIterator iend(FaceGeometry::end(box, face_normal));
   for (FaceIterator i(FaceGeometry::begin(box, face_normal)); i != iend; ++i) {
      os << "array" << *i << " = "
         << (*(d_data[face_normal]))(*i, depth) << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in the d_depth data member from the restart database.  Then tells
 * d_data to read itself in from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
FaceData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_FACEDATA_VERSION");
   if (ver != PDAT_FACEDATA_VERSION) {
      TBOX_ERROR("FaceData<getDim()>::getFromRestart error...\n"
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
FaceData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_FACEDATA_VERSION", PDAT_FACEDATA_VERSION);

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
typename FaceData<TYPE>::template View<DIM> get_view(FaceData<TYPE>& data,
                                                     Args&&... args)
{
  return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename FaceData<TYPE>::template ConstView<DIM> get_const_view(
    const FaceData<TYPE>& data,
    Args&&... args)
{
  return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#endif
