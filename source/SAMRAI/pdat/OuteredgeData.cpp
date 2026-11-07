/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outeredge centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeData_C
#define included_pdat_OuteredgeData_C

#include "SAMRAI/pdat/OuteredgeData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OuteredgeData<TYPE>::PDAT_OUTEREDGEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outeredge data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OuteredgeData<TYPE>::OuteredgeData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {

      for (tbox::Dimension::dir_t face_normal = 0; face_normal < getDim().getValue();
           ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               hier::Box oedge_data_box =
                  OuteredgeGeometry::toOuteredgeBox(getGhostBox(),
                     axis,
                     face_normal,
                     side);

               d_data[axis][face_normal][side].reset(
                  new ArrayData<TYPE>(oedge_data_box, depth));

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
OuteredgeData<TYPE>::OuteredgeData(
   const hier::Box& box,
   int depth,
   tbox::ResourceAllocator allocator):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {

      for (tbox::Dimension::dir_t face_normal = 0; face_normal < getDim().getValue();
           ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               hier::Box oedge_data_box =
                  OuteredgeGeometry::toOuteredgeBox(getGhostBox(),
                     axis,
                     face_normal,
                     side);

               d_data[axis][face_normal][side].reset(
                  new ArrayData<TYPE>(oedge_data_box, depth, allocator));

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
OuteredgeData<TYPE>::~OuteredgeData()
{
}

template<class TYPE>
int
OuteredgeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
bool
OuteredgeData<TYPE>::dataExists(
   int axis,
   int face_normal) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[axis][face_normal][0]->isInitialized();
}

template<class TYPE>
TYPE *
OuteredgeData<TYPE>::getPointer(
   int axis,
   int face_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][face_normal][side]->getPointer(depth);
}

template<class TYPE>
const TYPE *
OuteredgeData<TYPE>::getPointer(
   int axis,
   int face_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][face_normal][side]->getPointer(depth);
}

#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename OuteredgeData<TYPE>::template View<DIM> OuteredgeData<TYPE>::getView(
   int axis,
   int face_normal,
   int side,
   int depth)
{
   ArrayData<TYPE>& array_data = getArrayData(axis, face_normal, side);
   return array_data.getView(depth);
}

template <class TYPE>
template <int DIM>
typename OuteredgeData<TYPE>::template ConstView<DIM>
OuteredgeData<TYPE>::getConstView(
   int axis,
   int face_normal,
   int side,
   int depth) const
{
   const ArrayData<TYPE>& array_data = getArrayData(axis, face_normal, side);
   return array_data.getConstView(depth);
}
#endif

template<class TYPE>
TYPE&
OuteredgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

      if (face_normal != axis) {

         for (int side = 0; side < 2; ++side) {

            if (d_data[axis][face_normal][side]->getBox().contains(i)) {
               return (*(d_data[axis][face_normal][side]))(i, depth);
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
      " edge index " << i << " with axis = " << axis << std::endl);
   return (*(d_data[0][0][0]))(i, depth);
}

template<class TYPE>
const TYPE&
OuteredgeData<TYPE>::operator () (
   const EdgeIndex& i,
   int depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

      if (face_normal != axis) {

         for (int side = 0; side < 2; ++side) {

            if (d_data[axis][face_normal][side]->getBox().contains(i)) {
               return (*(d_data[axis][face_normal][side]))(i, depth);
            }

         }  // iterate over lower/upper sides

      }  // data is undefined when axis == face_normal

   }  // iterate over face normal directions

   TBOX_ERROR("Attempt to access OuteredgeData value with bad index"
      " edge index " << i << " with axis = " << axis << std::endl);
   return (*(d_data[0][0][0]))(i, depth);
}

template<class TYPE>
ArrayData<TYPE>&
OuteredgeData<TYPE>::getArrayData(
   int axis,
   int face_normal,
   int side)
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[axis][face_normal][side]);
}

template<class TYPE>
const ArrayData<TYPE>&
OuteredgeData<TYPE>::getArrayData(
   int axis,
   int face_normal,
   int side) const
{
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[axis][face_normal][side]);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outeredge patch data type (source) and
 * a edge patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const EdgeData<TYPE> * const t_edge_src =
      dynamic_cast<const EdgeData<TYPE> *>(&src);
   const OuteredgeData<TYPE> * const t_oedge_src =
      dynamic_cast<const OuteredgeData<TYPE> *>(&src);

   if (t_edge_src != 0) {
      copyFromEdge(*t_edge_src);
   } else if (t_oedge_src != 0) {
      copyFromOuteredge(*t_oedge_src);
   } else {
      TBOX_ERROR("OuteredgeData<getDim()>::copy error!\n"
         << "Can copy only from EdgeData<TYPE> or "
         << "OuteredgeData<TYPE> of the same "
         << "getDim() and TYPE.");
   }

}

template<class TYPE>
void
OuteredgeData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   EdgeData<TYPE>* t_edge_dst = dynamic_cast<EdgeData<TYPE> *>(&dst);
   OuteredgeData<TYPE>* t_oedge_dst = dynamic_cast<OuteredgeData<TYPE> *>(&dst);

   if (t_edge_dst != 0) {
      copyToEdge(*t_edge_dst);
   } else if (t_oedge_dst != 0) {
      copyToOuteredge(*t_oedge_dst);
   } else {
      TBOX_ERROR("OuteredgeData<getDim()>::copy2 error!\n"
         << "Can copy only to EdgeData<TYPE> or "
         << "OuteredgeData<TYPE> of the same "
         << "getDim() and TYPE.");
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
OuteredgeData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const EdgeData<TYPE>* t_edge_src =
      dynamic_cast<const EdgeData<TYPE> *>(&src);
   const OuteredgeData<TYPE>* t_oedge_src =
      dynamic_cast<const OuteredgeData<TYPE> *>(&src);

   if (t_edge_src != 0) {
      copyFromEdge(*t_edge_src, *t_overlap);
   } else if (t_oedge_src != 0) {
      copyFromOuteredge(*t_oedge_src, *t_overlap);
   } else {
      TBOX_ERROR("OuternodeData<getDim()>::copy error!\n"
         << "Can copy only from EdgeData<TYPE> or "
         << "OuteredgeData<TYPE> of the same "
         << "getDim() and TYPE.");
   }

}

template<class TYPE>
void
OuteredgeData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   EdgeData<TYPE>* t_edge_dst = dynamic_cast<EdgeData<TYPE> *>(&dst);
   OuteredgeData<TYPE>* t_oedge_dst = dynamic_cast<OuteredgeData<TYPE> *>(&dst);

   if (t_edge_dst != 0) {
      copyToEdge(*t_edge_dst, *t_overlap);
   } else if (t_oedge_dst != 0) {
      copyToOuteredge(*t_oedge_dst, *t_overlap);
   } else {
      TBOX_ERROR("OuternodeData<getDim()>::copy2 error!\n"
         << "Can copy only to EdgeData<TYPE> or "
         << "OuteredgeData<TYPE> of the same "
         << "getDim() and TYPE.");
   }

}

/*
 *************************************************************************
 *
 * Perform a fast copy from an edge data object to this outeredge data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::copyDepth(
   int dst_depth,
   const EdgeData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const ArrayData<TYPE>& src_edge_array = src.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               ArrayData<TYPE>& dst_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_oedge_array.copyDepth(dst_depth,
                  src_edge_array,
                  src_depth,
                  dst_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
 *************************************************************************
 *
 * Perform a fast copy to an edge data object from this outeredge data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::copyDepth2(
   int dst_depth,
   EdgeData<TYPE>& dst,
   int src_depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      ArrayData<TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               const ArrayData<TYPE>& src_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_edge_array.copyDepth(dst_depth,
                  src_oedge_array,
                  src_depth,
                  src_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
 *************************************************************************
 *
 * Add source data to the destination according to overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::sum(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const OuteredgeData<TYPE>* t_oedge_src =
      CPP_CAST<const OuteredgeData<TYPE> *>(&src);

   TBOX_ASSERT(t_oedge_src != 0);

   // NOTE:  We assume this operation is only needed to
   //        copy and add data to another outeredge data
   //        object.  If we ever need to provide this for edge
   //        data or other flavors of the copy operation, we
   //        should refactor the routine similar to the way
   //        the regular copy operations are implemented.
   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& box_list =
         t_overlap->getDestinationBoxContainer(axis);

      for (int src_face_normal = 0;
           src_face_normal < getDim().getValue();
           ++src_face_normal) {

         if (src_face_normal != axis) {

            for (int src_side = 0; src_side < 2; ++src_side) {

               if (t_oedge_src->d_data[axis][src_face_normal][src_side]->
                   isInitialized()) {

                  const ArrayData<TYPE>& src_array =
                     *(t_oedge_src->d_data[axis][src_face_normal][src_side]);

                  for (int dst_face_normal = 0;
                       dst_face_normal < getDim().getValue(); ++dst_face_normal) {

                     if (dst_face_normal != axis) {

                        for (int dst_side = 0; dst_side < 2; ++dst_side) {

                           if (d_data[axis][dst_face_normal][dst_side]->
                               isInitialized()) {

                              d_data[axis][dst_face_normal][dst_side]->
                              sum(src_array,
                                 box_list,
                                 src_offset);

                           }  // if dst data array is initialized

                        }  // iterate over dst lower/upper sides

                     }  // dst data is undefined when axis == face_normal

                  }  // iterate over dst face normal directions

               }  // if src data array is initialized

            }  // iterate over src lower/upper sides

         }  // src data is undefined when axis == face_normal

      }  // iterate over src face normal directions

   }  // iterate over axis directions

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
OuteredgeData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
OuteredgeData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   size_t size = 0;

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   for (tbox::Dimension::dir_t axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& boxlist =
         t_overlap->getDestinationBoxContainer(axis);

      for (tbox::Dimension::dir_t face_normal = 0; face_normal < getDim().getValue();
           ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               if (d_data[axis][face_normal][side]->isInitialized()) {

                  size += d_data[axis][face_normal][side]->
                     getDataStreamSize(boxlist, src_offset);

               }   // if data arrays is initialized

            }   // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

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
OuteredgeData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& dst_boxes =
         t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator dst_box = dst_boxes.begin();
           dst_box != dst_boxes.end(); ++dst_box) {

         const hier::Box src_box = hier::Box::shift(*dst_box,
               -src_offset);

         for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

            if (face_normal != axis) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box intersection =
                     src_box * d_data[axis][face_normal][side]->getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side]->
                     packStream(stream,
                        hier::Box::shift(intersection,
                           src_offset),
                        src_offset);

                  } // if intersection non-empty

               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& dst_boxes =
         t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator dst_box = dst_boxes.begin();
           dst_box != dst_boxes.end(); ++dst_box) {

         for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

            if (face_normal != axis) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box intersection =
                     (*dst_box) * d_data[axis][face_normal][side]->getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side]->unpackStream(stream,
                        intersection,
                        src_offset);

                  } // if intersection non-empty

               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

/*
 *************************************************************************
 *
 * Unpack data from the message stream and add to this outeredge data
 * object using the index space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& dst_boxes =
         t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator dst_box = dst_boxes.begin();
           dst_box != dst_boxes.end(); ++dst_box) {

         for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

            if (face_normal != axis) {

               for (int side = 0; side < 2; ++side) {

                  const hier::Box intersection =
                     (*dst_box) * d_data[axis][face_normal][side]->getBox();

                  if (!intersection.empty()) {

                     d_data[axis][face_normal][side]->
                     unpackStreamAndSum(stream,
                        intersection,
                        src_offset);

                  } // if intersection non-empty

               }  // iterate over lower/upper sides

            }  // data is undefined when axis == face_normal

         }  // iterate over face normal directions

      }  // iterate over overlap boxes

   }  // iterate over axis directions

}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outeredge centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OuteredgeData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;

   for (tbox::Dimension::dir_t axis = 0; axis < box.getDim().getValue(); ++axis) {

      for (tbox::Dimension::dir_t face_normal = 0;
           face_normal < box.getDim().getValue();
           ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               hier::Box oedge_data_box =
                  OuteredgeGeometry::toOuteredgeBox(box,
                     axis,
                     face_normal,
                     side);

               size +=
                  ArrayData<TYPE>::getSizeOfData(oedge_data_box, depth);

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

   return size;
}

/*
 *************************************************************************
 *
 * Compute the box of valid edge indices given values of
 * direction and side designating the set of data indices.
 *
 *************************************************************************
 */

template<class TYPE>
hier::Box
OuteredgeData<TYPE>::getDataBox(
   int axis,
   int face_normal,
   int side)
{
   return OuteredgeGeometry::toOuteredgeBox(getGhostBox(),
      axis,
      face_normal,
      side);
}

/*
 *************************************************************************
 *
 * Fill the outeredge centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               if (d_data[axis][face_normal][side]->isInitialized()) {
                  d_data[axis][face_normal][side]->fill(t, d);
               }

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      hier::Box databox = EdgeGeometry::toEdgeBox(box, axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               if (d_data[axis][face_normal][side]->isInitialized()) {
                  d_data[axis][face_normal][side]->fill(t, databox, d);
               }

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::fillAll(
   const TYPE& t)
{

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               if (d_data[axis][face_normal][side]->isInitialized()) {
                  d_data[axis][face_normal][side]->fillAll(t);
               }

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      hier::Box databox = EdgeGeometry::toEdgeBox(box, axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               if (d_data[axis][face_normal][side]->isInitialized()) {
                  d_data[axis][face_normal][side]->fillAll(t, databox);
               }

            }   // iterate over lower/upper sides

         }   // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outeredge patch data type (source) and
 * a edge patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::copyFromEdge(
   const EdgeData<TYPE>& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const ArrayData<TYPE>& src_edge_array = src.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               ArrayData<TYPE>& dst_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_oedge_array.copy(src_edge_array,
                  dst_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::copyToEdge(
   EdgeData<TYPE>& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      ArrayData<TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               const ArrayData<TYPE>& src_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_edge_array.copy(src_oedge_array,
                  src_oedge_array.getBox());

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }   // iterate over axis directions

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
OuteredgeData<TYPE>::copyFromEdge(
   const EdgeData<TYPE>& src,
   const EdgeOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const hier::IntVector& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer(axis);
      const ArrayData<TYPE>& src_edge_array = src.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               ArrayData<TYPE>& dst_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_oedge_array.copy(src_edge_array,
                  box_list,
                  src_offset);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::copyToEdge(
   EdgeData<TYPE>& dst,
   const EdgeOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   const hier::IntVector& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer(axis);
      ArrayData<TYPE>& dst_edge_array = dst.getArrayData(axis);

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               const ArrayData<TYPE>& src_oedge_array =
                  *(d_data[axis][face_normal][side]);

               dst_edge_array.copy(src_oedge_array,
                  box_list,
                  src_offset);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }   // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::copyFromOuteredge(
   const OuteredgeData<TYPE>& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int src_face_normal = 0;
           src_face_normal < getDim().getValue();
           ++src_face_normal) {

         if (src_face_normal != axis) {

            for (int src_side = 0; src_side < 2; ++src_side) {

               const ArrayData<TYPE>& src_oedge_array =
                  *(src.d_data[axis][src_face_normal][src_side]);

               for (int dst_face_normal = 0;
                    dst_face_normal < getDim().getValue();
                    ++dst_face_normal) {

                  if (dst_face_normal != axis) {

                     for (int dst_side = 0; dst_side < 2; ++dst_side) {

                        ArrayData<TYPE>& dst_oedge_array =
                           *(d_data[axis][dst_face_normal][dst_side]);

                        dst_oedge_array.copy(src_oedge_array,
                           dst_oedge_array.getBox());

                     }  // iterate over dst lower/upper sides

                  }  // dst data is undefined when axis == face_normal

               }  // iterate over dst face normal directions

            }  // iterate over src lower/upper sides

         }  // src data is undefined when axis == face_normal

      }  // iterate over src face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::copyFromOuteredge(
   const OuteredgeData<TYPE>& src,
   const EdgeOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const hier::IntVector& src_offset = overlap.getSourceOffset();

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      const hier::BoxContainer& box_list =
         overlap.getDestinationBoxContainer(axis);

      for (int src_face_normal = 0;
           src_face_normal < getDim().getValue();
           ++src_face_normal) {

         if (src_face_normal != axis) {

            for (int src_side = 0; src_side < 2; ++src_side) {

               const ArrayData<TYPE>& src_oedge_array =
                  *(src.d_data[axis][src_face_normal][src_side]);

               for (int dst_face_normal = 0;
                    dst_face_normal < getDim().getValue();
                    ++dst_face_normal) {

                  if (dst_face_normal != axis) {

                     for (int dst_side = 0; dst_side < 2; ++dst_side) {

                        ArrayData<TYPE>& dst_oedge_array =
                           *(d_data[axis][dst_face_normal][dst_side]);

                        dst_oedge_array.copy(src_oedge_array,
                           box_list,
                           src_offset);

                     }  // iterate over dst lower/upper sides

                  }  // dst data is undefined when axis == face_normal

               }  // iterate over dst face normal directions

            }  // iterate over src lower/upper sides

         }  // src data is undefined when axis == face_normal

      }  // iterate over src face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::copyToOuteredge(
   OuteredgeData<TYPE>& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   dst.copyFromOuteredge(*this);
}

template<class TYPE>
void
OuteredgeData<TYPE>::copyToOuteredge(
   OuteredgeData<TYPE>& dst,
   const EdgeOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   dst.copyFromOuteredge(*this, overlap);
}

/*
 *************************************************************************
 *
 * Print routines for outeredge centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::print(
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
OuteredgeData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         os << "Array axis, face normal = "
            << axis << "," << face_normal << std::endl;

         for (int side = 0; side < 2; ++side) {

            os << "side  = "
               << ((side == 0) ? "lower" : "upper") << std::endl;

            printAxisSide(axis, face_normal, side,
               box, depth, os, prec);

         }  // iterate over lower/upper sides

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

template<class TYPE>
void
OuteredgeData<TYPE>::printAxisSide(
   int axis,
   int face_normal,
   int side,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(axis, face_normal, side,
         box, d, os, prec);
   }

}

template<class TYPE>
void
OuteredgeData<TYPE>::printAxisSide(
   int axis,
   int face_normal,
   int side,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   NULL_USE(prec);

   if (axis == face_normal) {

      os << "array data undefined" << std::endl;

   } else {

      const hier::Box edgebox = EdgeGeometry::toEdgeBox(box, axis);
      const hier::Box region =
         edgebox * d_data[axis][face_normal][side]->getBox();
      hier::Box::iterator iiend(region.end());
      for (hier::Box::iterator ii(region.begin()); ii != iiend; ++ii) {
         os << "array" << *ii << " = "
            << (*(d_data[axis][face_normal][side]))(*ii, depth) << std::endl;
         os << std::flush;
      }

   }

}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.
 * If so, reads in d_depth from the restart database.
 * Then has each item in d_data read in its data from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuteredgeData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_OUTEREDGEDATA_VERSION");
   if (ver != PDAT_OUTEREDGEDATA_VERSION) {
      TBOX_ERROR(
         "OuteredgeData<getDim()>::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");

   std::shared_ptr<tbox::Database> array_database;

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {
               std::string array_name = "d_data" + tbox::Utilities::intToString(
                     axis)
                  + "_" + tbox::Utilities::intToString(face_normal) + "_"
                  + tbox::Utilities::intToString(side);
               array_database = restart_db->getDatabase(array_name);
               d_data[axis][face_normal][side]->getFromRestart(array_database);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

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
OuteredgeData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_OUTEREDGEDATA_VERSION",
      PDAT_OUTEREDGEDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);

   std::shared_ptr<tbox::Database> array_database;

   for (int axis = 0; axis < getDim().getValue(); ++axis) {

      for (int face_normal = 0; face_normal < getDim().getValue(); ++face_normal) {

         if (face_normal != axis) {

            for (int side = 0; side < 2; ++side) {

               std::string array_name = "d_data" + tbox::Utilities::intToString(
                     axis)
                  + "_" + tbox::Utilities::intToString(face_normal) + "_"
                  + tbox::Utilities::intToString(side);
               array_database = restart_db->putDatabase(array_name);
               d_data[axis][face_normal][side]->putToRestart(
                  array_database);

            }  // iterate over lower/upper sides

         }  // data is undefined when axis == face_normal

      }  // iterate over face normal directions

   }  // iterate over axis directions

}

#if defined(HAVE_RAJA)
template <int DIM, typename TYPE, typename... Args>
typename OuteredgeData<TYPE>::template View<DIM> get_view(OuteredgeData<TYPE>& data,
                                                     Args&&... args)
{
  return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename OuteredgeData<TYPE>::template ConstView<DIM> get_const_view(
    const OuteredgeData<TYPE>& data,
    Args&&... args)
{
  return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif


}
}

#endif
