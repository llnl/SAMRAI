/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated outerside centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OutersideData_C
#define included_pdat_OutersideData_C

#include "SAMRAI/pdat/OutersideData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/tbox/Utilities.h"
#include <stdio.h>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OutersideData<TYPE>::PDAT_OUTERSIDEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outerside data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OutersideData<TYPE>::OutersideData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box sidebox = SideGeometry::toSideBox(ghosts, d);
      hier::Box outersidebox = sidebox;
      outersidebox.setUpper(d, sidebox.lower(d));
      d_data[d][0].reset(new ArrayData<TYPE>(outersidebox, depth));
      outersidebox.setLower(d, sidebox.upper(d));
      outersidebox.setUpper(d, sidebox.upper(d));
      d_data[d][1].reset(new ArrayData<TYPE>(outersidebox, depth));
   }
}

template<class TYPE>
OutersideData<TYPE>::OutersideData(
   const hier::Box& box,
   int depth,
   tbox::ResourceAllocator allocator):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{

   TBOX_ASSERT(depth > 0);

   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box sidebox = SideGeometry::toSideBox(ghosts, d);
      hier::Box outersidebox = sidebox;
      outersidebox.setUpper(d, sidebox.lower(d));
      d_data[d][0].reset(new ArrayData<TYPE>(outersidebox, depth, allocator));
      outersidebox.setLower(d, sidebox.upper(d));
      outersidebox.setUpper(d, sidebox.upper(d));
      d_data[d][1].reset(new ArrayData<TYPE>(outersidebox, depth, allocator));
   }
}

template<class TYPE>
OutersideData<TYPE>::~OutersideData()
{
}

template<class TYPE>
int
OutersideData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE *
OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side]->getPointer(depth);
}

template<class TYPE>
const TYPE *
OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side]->getPointer(depth);
}


#if defined(HAVE_RAJA)
template <class TYPE>
template <int DIM>
typename OutersideData<TYPE>::template View<DIM> OutersideData<TYPE>::getView(
    int side_normal,
    int side,
    int depth)
{
  const hier::Box& box = getGhostBox();
  hier::Box outersidebox = box;
  if(DIM < 3) { // we can get away with not transposing since the slices are 1d
     if(side_normal == 0) {
       if(side == 0) {
          outersidebox.setUpper(0, box.lower(0));
       } else if (side==1) {
          outersidebox.setLower(0, box.upper(0));
          outersidebox.setUpper(0, box.upper(0));
       }
     } else if(side_normal == 1) {
       if(side == 0 ) {
          outersidebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outersidebox.setLower(1, box.upper(1));
          outersidebox.setUpper(1, box.upper(1));
       }
     } // else side_normal == 1
  } else if(DIM == 3) {

    if(side_normal == 0) {
       if(side == 0) {
          outersidebox.setUpper(0, box.lower(0));
       } else if (side == 1) {
          outersidebox.setLower(0, box.upper(0));
          outersidebox.setUpper(0, box.upper(0));
       }
    }
    else if(side_normal == 1) {
       if(side == 0) {
          outersidebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outersidebox.setLower(1, box.upper(1));
          outersidebox.setUpper(1, box.upper(1));
       }   
    } else if(side_normal == 2) {
       if(side == 0) {
          outersidebox.setUpper(2, box.lower(2));
       } else if (side == 1) {
          outersidebox.setLower(2, box.upper(2));
          outersidebox.setUpper(2, box.upper(2));
       }
    }
    
  } // if DIM == 3
  return OutersideData<TYPE>::View<DIM>(getPointer(side_normal,side, depth), outersidebox);
}

template <class TYPE>
template <int DIM>
typename OutersideData<TYPE>::template ConstView<DIM> OutersideData<TYPE>::getConstView(
    int side_normal,
    int side,
    int depth) const
{
  const hier::Box& box = getGhostBox();
  hier::Box outersidebox = box;
  if(DIM < 3) { // we can get away with not transposing since the slices are 1d
     if(side_normal == 0) {
       if(side == 0) {
          outersidebox.setUpper(0, box.lower(0));
       } else if (side==1) {
          outersidebox.setLower(0, box.upper(0));
          outersidebox.setUpper(0, box.upper(0));
       }
     } else if(side_normal == 1) {
       if(side == 0 ) {
          outersidebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outersidebox.setLower(1, box.upper(1));
          outersidebox.setUpper(1, box.upper(1));
       }
     } // else side_normal == 1
  } else if(DIM == 3) {

    if(side_normal == 0) {
       if(side == 0) {
          outersidebox.setUpper(0, box.lower(0));
       } else if (side == 1) {
          outersidebox.setLower(0, box.upper(0));
          outersidebox.setUpper(0, box.upper(0));
       }
    } else if(side_normal == 1) {
       if(side == 0) {
          outersidebox.setUpper(1, box.lower(1));
       } else if (side == 1) {
          outersidebox.setLower(1, box.upper(1));
          outersidebox.setUpper(1, box.upper(1));
       }   
    } else if(side_normal == 2) {
       if(side == 0) {
          outersidebox.setUpper(2, box.lower(2));
       } else if (side == 1) {
          outersidebox.setLower(2, box.upper(2));
          outersidebox.setUpper(2, box.upper(2));
       }
    }
  } // if DIM == 3
  return OutersideData<TYPE>::ConstView<DIM>(getPointer(side_normal, side, depth),
                                        outersidebox);
}
#endif


template<class TYPE>
ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[side_normal][side]);
}

template<class TYPE>
const ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return *(d_data[side_normal][side]);
}

template<class TYPE>
TYPE&
OutersideData<TYPE>::operator () (
   const SideIndex& i,
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
OutersideData<TYPE>::operator () (
   const SideIndex& i,
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
 * Perform a fast copy between an outerside patch data type (source) and
 * a side patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const SideData<TYPE> * const t_src = CPP_CAST<const SideData<TYPE> *>(&src);

   TBOX_ASSERT(t_src != 0);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      const ArrayData<TYPE>& side_array = t_src->getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         ArrayData<TYPE>& oside_array = *(d_data[axis][loc]);
         oside_array.copy(side_array, oside_array.getBox());
      }
   }

}

template<class TYPE>
void
OutersideData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   SideData<TYPE>* t_dst = CPP_CAST<SideData<TYPE> *>(&dst);

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
OutersideData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{

   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const SideOverlap* t_overlap = CPP_CAST<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const OutersideData<TYPE>* t_oside_src =
      dynamic_cast<const OutersideData<TYPE> *>(&src);
   const SideData<TYPE>* t_side_src =
      dynamic_cast<const SideData<TYPE> *>(&src);

   TBOX_ASSERT(t_oside_src == 0 || t_side_src == 0);
   TBOX_ASSERT(t_oside_src != 0 || t_side_src != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   if (t_oside_src != 0) {
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& box_list =
            t_overlap->getDestinationBoxContainer(d);
         d_data[d][0]->copy(t_oside_src->getArrayData(d, 0), box_list, src_offset);
         d_data[d][0]->copy(t_oside_src->getArrayData(d, 1), box_list, src_offset);
         d_data[d][1]->copy(t_oside_src->getArrayData(d, 0), box_list, src_offset);
         d_data[d][1]->copy(t_oside_src->getArrayData(d, 1), box_list, src_offset);
      }
   } else if (t_side_src != 0) {
      for (int d = 0; d < getDim().getValue(); ++d) {
         const hier::BoxContainer& box_list =
            t_overlap->getDestinationBoxContainer(d);
         d_data[d][0]->copy(t_side_src->getArrayData(d), box_list, src_offset);
         d_data[d][1]->copy(t_side_src->getArrayData(d), box_list, src_offset);
      }
   } else {
      TBOX_ERROR("OutersideData<TYPE>::copy error...\n"
         << " : Cannot copy from type other than SideData or OutersideData " << std::endl);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   SideData<TYPE>* t_dst = CPP_CAST<SideData<TYPE> *>(&dst);
   const SideOverlap* t_overlap = CPP_CAST<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != 0);
   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
      t_dst->getArrayData(d).copy(*(d_data[d][0]), box_list, src_offset);
      t_dst->getArrayData(d).copy(*(d_data[d][1]), box_list, src_offset);
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a side data object to this outerside data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copyDepth(
   int dst_depth,
   const SideData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      const ArrayData<TYPE>& src_side_array = src.getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         ArrayData<TYPE>& dst_oside_array = *(d_data[axis][loc]);
         dst_oside_array.copyDepth(dst_depth,
            src_side_array,
            src_depth,
            dst_oside_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy to a side data object from this outerside data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copyDepth2(
   int dst_depth,
   SideData<TYPE>& dst,
   int src_depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); ++axis) {
      ArrayData<TYPE>& dst_side_array = dst.getArrayData(axis);
      for (int loc = 0; loc < 2; ++loc) {
         const ArrayData<TYPE>& src_oside_array = *(d_data[axis][loc]);
         dst_side_array.copyDepth(dst_depth,
            src_oside_array,
            src_depth,
            src_oside_array.getBox());
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
OutersideData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
OutersideData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap = CPP_CAST<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   size_t size = 0;
   for (tbox::Dimension::dir_t d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxlist = t_overlap->getDestinationBoxContainer(d);
      size += d_data[d][0]->getDataStreamSize(boxlist, src_offset);
      size += d_data[d][1]->getDataStreamSize(boxlist, src_offset);
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
OutersideData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap = CPP_CAST<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {
         const hier::Box src_box = hier::Box::shift(*b, -src_offset);
         for (int f = 0; f < 2; ++f) {
            const hier::Box intersect = src_box * d_data[d][f]->getBox();
            if (!intersect.empty()) {
               d_data[d][f]->packStream(stream,
                  hier::Box::shift(intersect, src_offset),
                  src_offset);
            }
         }
      }
   }
}

template<class TYPE>
void
OutersideData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const SideOverlap* t_overlap = CPP_CAST<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); ++d) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {
         for (int f = 0; f < 2; ++f) {
            const hier::Box intersect = (*b) * d_data[d][f]->getBox();
            if (!intersect.empty()) {
               d_data[d][f]->unpackStream(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outerside centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OutersideData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   for (tbox::Dimension::dir_t d = 0; d < box.getDim().getValue(); ++d) {
      hier::Box lower = SideGeometry::toSideBox(box, d);
      hier::Box upper = SideGeometry::toSideBox(box, d);
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
 * Fill the outerside centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::fill(
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
OutersideData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fill(t, SideGeometry::toSideBox(box, i), d);
      d_data[i][1]->fill(t, SideGeometry::toSideBox(box, i), d);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fillAll(t);
      d_data[i][1]->fillAll(t);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int i = 0; i < getDim().getValue(); ++i) {
      d_data[i][0]->fillAll(t, SideGeometry::toSideBox(box, i));
      d_data[i][1]->fillAll(t, SideGeometry::toSideBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print routines for outerside centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::print(
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
OutersideData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (tbox::Dimension::dir_t side_normal = 0;
        side_normal < getDim().getValue(); ++side_normal) {
      os << "Array side normal  = " << side_normal << std::endl;
      for (int side = 0; side < 2; ++side) {
         os << "side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisSide(side_normal, side, box, depth, os, prec);
      }
   }
}

template<class TYPE>
void
OutersideData<TYPE>::printAxisSide(
   tbox::Dimension::dir_t side_normal,
   int side,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT(side_normal < getDim().getValue());
   TBOX_ASSERT((side == 0) || (side == 1));

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(side_normal, side, box, d, os, prec);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::printAxisSide(
   tbox::Dimension::dir_t side_normal,
   int side,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT(side_normal < getDim().getValue());
   TBOX_ASSERT((side == 0) || (side == 1));

   const hier::Box sidebox =
      SideGeometry::toSideBox(box, side_normal);
   const hier::Box region =
      sidebox * d_data[side_normal][side]->getBox();
   os.precision(prec);
   hier::Box::iterator iend(region.end());
   for (hier::Box::iterator i(region.begin()); i != iend; ++i) {
      os << "array" << *i << " = "
         << (*(d_data[side_normal][side]))(*i, depth) << std::endl;
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
OutersideData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_OUTERSIDEDATA_VERSION");
   if (ver != PDAT_OUTERSIDEDATA_VERSION) {
      TBOX_ERROR("OutersideData<TYPE>::getFromRestart error...\n"
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
OutersideData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_OUTERSIDEDATA_VERSION",
      PDAT_OUTERSIDEDATA_VERSION);

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
typename OutersideData<TYPE>::template View<DIM> get_view(OutersideData<TYPE>& data,
                                                     Args&&... args)
{
  return data.template getView<DIM>(std::forward<Args>(args)...);
}

template <int DIM, typename TYPE, typename... Args>
typename OutersideData<TYPE>::template ConstView<DIM> get_const_view(
    const OutersideData<TYPE>& data,
    Args&&... args)
{
  return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#endif
