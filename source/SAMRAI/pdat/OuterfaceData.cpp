/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outerface centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceData_C
#define included_pdat_OuterfaceData_C

#include "SAMRAI/pdat/OuterfaceData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/OuterfaceGeometry.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OuterfaceData<TYPE>::PDAT_OUTERFACEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outerface data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OuterfaceData<TYPE>::OuterfaceData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{

   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box facebox = FaceGeometry::toFaceBox(ghosts, d);
      hier::Box outerfacebox = facebox;
      outerfacebox.setUpper(0, facebox.lower(0));
      d_data[d][0].reset(new ArrayData<TYPE>(outerfacebox, depth));
      outerfacebox.setLower(0, facebox.upper(0));
      outerfacebox.setUpper(0, facebox.upper(0));
      d_data[d][1].reset(new ArrayData<TYPE>(outerfacebox, depth));
   }
}

template<class TYPE>
OuterfaceData<TYPE>::OuterfaceData(
   const hier::Box& box,
   int depth,
   tbox::ResourceAllocator allocator):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{

   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box facebox = FaceGeometry::toFaceBox(ghosts, d);
      hier::Box outerfacebox = facebox;
      outerfacebox.setUpper(0, facebox.lower(0));
      d_data[d][0].reset(new ArrayData<TYPE>(outerfacebox, depth,allocator));
      outerfacebox.setLower(0, facebox.upper(0));
      outerfacebox.setUpper(0, facebox.upper(0));
      d_data[d][1].reset(new ArrayData<TYPE>(outerfacebox, depth,allocator));
   }
}

template<class TYPE>
OuterfaceData<TYPE>::~OuterfaceData()
{
}

template<class TYPE>
int
OuterfaceData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE *
OuterfaceData<TYPE>::getPointer(
   int face_normal,
   int side,
   int d)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_data[face_normal][side]->getPointer(d);
}

template<class TYPE>
const TYPE *
OuterfaceData<TYPE>::getPointer(
   int face_normal,
   int side,
   int d) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_data[face_normal][side]->getPointer(d);
}

#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename OuterfaceData<TYPE>::template View<DIM> OuterfaceData<TYPE>::getView(
    int face_normal,
    int side,
    int depth)
{
  const hier::Box& box = getGhostBox();
  hier::Box outerfacebox = box;
  if(DIM < 3) { // we can get away with not transposing since the slices are 1d
     if(face_normal == 0) {
       if(side == 0) {
          outerfacebox.setUpper(0, box.lower(0));
       } else if (side==1) {
          outerfacebox.setLower(0, box.upper(0));
          outerfacebox.setUpper(0, box.upper(0));
       }
     } else if(face_normal == 1) {
       if(side == 0 ) {
          outerfacebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outerfacebox.setLower(1, box.upper(1));
          outerfacebox.setUpper(1, box.upper(1));
       }
     } // else face_normal == 1
  } else if(DIM == 3) {

    if(face_normal == 0) {
       if(side == 0) {
          outerfacebox.setUpper(0, box.lower(0));
       } else if (side == 1) {
          outerfacebox.setLower(0, box.upper(0));
          outerfacebox.setUpper(0, box.upper(0));
       }
    }
    else if(face_normal == 1) {
       outerfacebox.setLower(0, box.lower(2));
       outerfacebox.setLower(1, box.lower(1));
       outerfacebox.setLower(2, box.lower(0));
       outerfacebox.setUpper(0, box.upper(2));
       outerfacebox.setUpper(1, box.upper(1));
       outerfacebox.setUpper(2, box.upper(0));
       if(side == 0) {
          outerfacebox.setUpper(1, outerfacebox.lower(1));
       } else if (side == 1) {
          outerfacebox.setLower(1, outerfacebox.upper(1));
       }   
    } else if(face_normal == 2) {
       if(side == 0) {
          outerfacebox.setUpper(2, box.lower(2));
       } else if (side == 1) {
          outerfacebox.setLower(2, box.upper(2));
          outerfacebox.setUpper(2, box.upper(2));
       }
    }
    
  } // if DIM == 3
  return OuterfaceData<TYPE>::View<DIM>(getPointer(face_normal,side, depth), outerfacebox);
}

template <class TYPE>
template <int DIM>
typename OuterfaceData<TYPE>::template ConstView<DIM> OuterfaceData<TYPE>::getConstView(
    int face_normal,
    int side,
    int depth) const
{
  const hier::Box& box = getGhostBox();
  hier::Box outerfacebox = box;
  if(DIM < 3) { // we can get away with not transposing since the slices are 1d
     if(face_normal == 0) {
       if(side == 0) {
          outerfacebox.setUpper(0, box.lower(0));
       } else if (side==1) {
          outerfacebox.setLower(0, box.upper(0));
          outerfacebox.setUpper(0, box.upper(0));
       }
     } else if(face_normal == 1) {
       if(side == 0 ) {
          outerfacebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outerfacebox.setLower(1, box.upper(1));
          outerfacebox.setUpper(1, box.upper(1));
       }
     } // else face_normal == 1
  } else if(DIM == 3) {

    if(face_normal == 0) {
       if(side == 0) {
          outerfacebox.setUpper(0, box.lower(0));
       } else if (side == 1) {
          outerfacebox.setLower(0, box.upper(0));
          outerfacebox.setUpper(0, box.upper(0));
       }
    } else if(face_normal == 1) {
       outerfacebox.setLower(0, box.lower(2));
       outerfacebox.setLower(1, box.lower(1));
       outerfacebox.setLower(2, box.lower(0));
       outerfacebox.setUpper(0, box.upper(2));
       outerfacebox.setUpper(1, box.upper(1));
       outerfacebox.setUpper(2, box.upper(0));
       if(side == 0) {
          outerfacebox.setUpper(1, outerfacebox.lower(1));
       } else if (side == 1) {
          outerfacebox.setLower(1, outerfacebox.upper(1));
       }   
    } else if(face_normal == 2) {
       if(side == 0) {
          outerfacebox.setUpper(2, box.lower(2));
       } else if (side == 1) {
          outerfacebox.setLower(2, box.upper(2));
          outerfacebox.setUpper(2, box.upper(2));
       }
    }
  } // if DIM == 3
  return OuterfaceData<TYPE>::ConstView<DIM>(getPointer(face_normal, side, depth),
                                        outerfacebox);
}
#endif

template<class TYPE>
ArrayData<TYPE>&
OuterfaceData<TYPE>::getArrayData(
   int face_normal,
   int side)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[face_normal][side]);
}

template<class TYPE>
const ArrayData<TYPE>&
OuterfaceData<TYPE>::getArrayData(
   int face_normal,
   int side) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[face_normal][side]);
}

template<class TYPE>
TYPE&
OuterfaceData<TYPE>::operator () (
   const FaceIndex& i,
   int side,
   int depth)
{
   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*(d_data[axis][side]))(i, depth);
}

template<class TYPE>
const TYPE&
OuterfaceData<TYPE>::operator () (
   const FaceIndex& i,
   int side,
   int depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*(d_data[axis][side]))(i, depth);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outerface patch data type (source) and
 * a face patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const FaceData<TYPE> * const t_src =
      CPP_CAST<const FaceData<TYPE> *>(&src);

   TBOX_ASSERT(t_src != 0);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      const ArrayData<TYPE>& face_array = t_src->getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         ArrayData<TYPE>& oface_array = *(d_data[axis][loc]);
         oface_array.copy(face_array, oface_array.getBox());
      }
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   FaceData<TYPE>* t_dst = CPP_CAST<FaceData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != 0);

   for (int d = 0; d < getDim().getValue(); ++d) {
      t_dst->getArrayData(d).copy(*(d_data[d][0]), d_data[d][0]->getBox());
      t_dst->getArrayData(d).copy(*(d_data[d][1]), d_data[d][1]->getBox());
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
OuterfaceData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const OuterfaceData<TYPE>* t_oface_src =
      dynamic_cast<const OuterfaceData<TYPE> *>(&src);
   const FaceData<TYPE>* t_face_src =
      dynamic_cast<const FaceData<TYPE> *>(&src);

   TBOX_ASSERT(t_oface_src == 0 || t_face_src == 0);
   TBOX_ASSERT(t_oface_src != 0 || t_face_src != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   if (t_oface_src != 0) {
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& box_list =
            t_overlap->getDestinationBoxContainer(d);
         d_data[d][0]->copy(t_oface_src->getArrayData(d, 0), box_list, src_offset);
         d_data[d][0]->copy(t_oface_src->getArrayData(d, 1), box_list, src_offset);
         d_data[d][1]->copy(t_oface_src->getArrayData(d, 0), box_list, src_offset);
         d_data[d][1]->copy(t_oface_src->getArrayData(d, 1), box_list, src_offset);
      }
   } else if (t_face_src != 0) {
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& box_list =
            t_overlap->getDestinationBoxContainer(d);
         d_data[d][0]->copy(t_face_src->getArrayData(d), box_list, src_offset);
         d_data[d][1]->copy(t_face_src->getArrayData(d), box_list, src_offset);
      }
   } else {
      TBOX_ERROR("OuterfaceData<TYPE>::copy error...\n"
         << " : Cannot copy from type other than FaceData or OuterfaceData " << std::endl);
   }

}

template<class TYPE>
void
OuterfaceData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   FaceData<TYPE>* t_dst = CPP_CAST<FaceData<TYPE> *>(&dst);
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != 0);
   TBOX_ASSERT(t_overlap != 0);

   const hier::Transformation& transformation = t_overlap->getTransformation();
   TBOX_ASSERT(transformation.getRotation() == hier::Transformation::NO_ROTATE);

   const hier::IntVector& src_offset = transformation.getOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      hier::IntVector face_offset(src_offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); ++i) {
            face_offset(i) = src_offset((d + i) % getDim().getValue());
         }
      }
      hier::Transformation face_transform(hier::Transformation::NO_ROTATE,
                                          face_offset,
                                          getBox().getBlockId(),
                                          t_dst->getBox().getBlockId());

      const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
      t_dst->getArrayData(d).copy(*(d_data[d][0]), box_list, face_transform);
      t_dst->getArrayData(d).copy(*(d_data[d][1]), box_list, face_transform);
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a face data object to this outerface data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copyDepth(
   int dst_depth,
   const FaceData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      const ArrayData<TYPE>& src_face_array = src.getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         ArrayData<TYPE>& dst_oface_array = *(d_data[axis][loc]);
         dst_oface_array.copyDepth(dst_depth,
            src_face_array,
            src_depth,
            dst_oface_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy to a face data object from this outerface data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copyDepth2(
   int dst_depth,
   FaceData<TYPE>& dst,
   int src_depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      ArrayData<TYPE>& dst_face_array = dst.getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         const ArrayData<TYPE>& src_oface_array = *(d_data[axis][loc]);
         dst_face_array.copyDepth(dst_depth,
            src_oface_array,
            src_depth,
            src_oface_array.getBox());
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
OuterfaceData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
OuterfaceData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   size_t size = 0;
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxlist = t_overlap->getDestinationBoxContainer(d);
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); ++i) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      size += d_data[d][0]->getDataStreamSize(boxlist, face_offset);
      size += d_data[d][1]->getDataStreamSize(boxlist, face_offset);
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
OuterfaceData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);

      if (!boxes.empty()) {
         hier::IntVector face_offset(offset);
         if (d > 0) {
            for (int i = 0; i < getDim().getValue(); ++i) {
               face_offset(i) = offset((d + i) % getDim().getValue());
            }
         }

         hier::Transformation face_transform(hier::Transformation::NO_ROTATE,
                                             face_offset,
                                             getBox().getBlockId(),
                                             boxes.begin()->getBlockId());

         for (hier::BoxContainer::const_iterator b = boxes.begin();
              b != boxes.end(); ++b) {
            hier::Box src_box(*b);
            face_transform.inverseTransform(src_box);
            for (int f = 0; f < 2; ++f) {
               hier::Box intersect(src_box * d_data[d][f]->getBox());
               if (!intersect.empty()) {
                  face_transform.transform(intersect);
                  d_data[d][f]->packStream(stream,
                     intersect,
                     face_transform);
               }
            }
         }
      }

   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const FaceOverlap* t_overlap = CPP_CAST<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); ++i) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }

      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {
         for (int f = 0; f < 2; ++f) {
            const hier::Box intersect = (*b) * d_data[d][f]->getBox();
            if (!intersect.empty()) {
               d_data[d][f]->unpackStream(stream, intersect, face_offset);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outerface centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OuterfaceData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   for (tbox::Dimension::dir_t d = 0; d < box.getDim().getValue(); ++d) {
      hier::Box lower = FaceGeometry::toFaceBox(box, d);
      hier::Box upper = FaceGeometry::toFaceBox(box, d);
      lower.setUpper(d, box.lower(d));
      upper.setLower(d, box.upper(d));
      size += ArrayData<TYPE>::getSizeOfData(lower, depth);
      size += ArrayData<TYPE>::getSizeOfData(upper, depth);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Fill the outerface centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fill(t, d);
      d_data[i][1]->fill(t, d);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fill(t, FaceGeometry::toFaceBox(box, i), d);
      d_data[i][1]->fill(t, FaceGeometry::toFaceBox(box, i), d);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fillAll(t);
      d_data[i][1]->fillAll(t);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fillAll(t, FaceGeometry::toFaceBox(box, i));
      d_data[i][1]->fillAll(t, FaceGeometry::toFaceBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print routines for outerface centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int d = 0; d < d_depth; ++d) {
      print(box, d, os, prec);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (tbox::Dimension::dir_t face_normal = 0;
        face_normal < getDim().getValue(); ++face_normal) {
      os << "Array face normal = " << face_normal << std::endl;
      for (int side = 0; side < 2; ++side) {
         os << "side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisFace(face_normal, side, box, depth, os, prec);
      }
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::printAxisFace(
   tbox::Dimension::dir_t face_normal,
   int side,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT(face_normal < getDim().getValue());
   TBOX_ASSERT((side == 0) || (side == 1));

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      printAxisFace(face_normal, side, box, d, os, prec);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::printAxisFace(
   tbox::Dimension::dir_t face_normal,
   int side,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT(face_normal < getDim().getValue());
   TBOX_ASSERT((side == 0) || (side == 1));

   const hier::Box facebox =
      FaceGeometry::toFaceBox(box, face_normal);
   const hier::Box region =
      facebox * d_data[face_normal][side]->getBox();
   os.precision(prec);
   hier::Box::iterator iend(region.end());
   for (hier::Box::iterator i(region.begin()); i != iend; ++i) {
      os << "array" << *i << " = "
         << (*(d_data[face_normal][side]))(*i, depth) << std::endl;
      os << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in d_depth from the restart database.  Then has each item in d_data
 * read in its data from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_OUTERFACEDATA_VERSION");
   if (ver != PDAT_OUTERFACEDATA_VERSION) {
      TBOX_ERROR(
         "OuterfaceData<getDim()>::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");

   std::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); ++i) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = restart_db->getDatabase(array_name);
      d_data[i][0]->getFromRestart(array_database);

      array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
      array_database = restart_db->getDatabase(array_name);
      d_data[i][1]->getFromRestart(array_database);
   }
}

/*
 *************************************************************************
 *
 * Writes out class version number, d_depth to the restart database.
 * Then has each item in d_data write out its data to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{

   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_OUTERFACEDATA_VERSION",
      PDAT_OUTERFACEDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);

   std::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); ++i) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = restart_db->putDatabase(array_name);
      d_data[i][0]->putToRestart(array_database);

      array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
      array_database = restart_db->putDatabase(array_name);
      d_data[i][1]->putToRestart(array_database);
   }
}

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OuterfaceData<TYPE>::template View<DIM> get_view(OuterfaceData<TYPE>& data,
                                                     Args&&... args)
{
  return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename OuterfaceData<TYPE>::template ConstView<DIM> get_const_view(
    const OuterfaceData<TYPE>& data,
    Args&&... args)
{
  return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#endif
