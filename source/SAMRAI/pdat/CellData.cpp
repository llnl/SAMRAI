/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated cell centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_CellData_C
#define included_pdat_CellData_C

#include "SAMRAI/pdat/CellData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int CellData<TYPE>::PDAT_CELLDATA_VERSION = 1;

template<class TYPE>
std::shared_ptr<tbox::Timer> CellData<TYPE>::t_copy;

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a cell centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
CellData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);

   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   return ArrayData<TYPE>::getSizeOfData(ghost_box, depth);
}

/*
 *************************************************************************
 *
 * Constructor and destructor for cell data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
CellData<TYPE>::CellData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   t_copy = tbox::TimerManager::getManager()->
      getTimer("pdat::CellData::copy");

   d_data.reset(new ArrayData<TYPE>(getGhostBox(), depth));
}

template<class TYPE>
CellData<TYPE>::CellData(
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

   t_copy = tbox::TimerManager::getManager()->
      getTimer("pdat::CellData::copy");

   d_data.reset(new ArrayData<TYPE>(getGhostBox(), depth, allocator));
}

template<class TYPE>
CellData<TYPE>::~CellData()
{
}

template<class TYPE>
int
CellData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE *
CellData<TYPE>::getPointer(
   int depth)
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data->getPointer(depth);
}

template<class TYPE>
const TYPE *
CellData<TYPE>::getPointer(
   int depth) const
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data->getPointer(depth);
}

#if defined(HAVE_RAJA)

template<class TYPE>
template<int DIM>
typename CellData<TYPE>::template View<DIM>
CellData<TYPE>::getView(int depth)
{
   return CellData<TYPE>::View<DIM>(getPointer(depth), getGhostBox());
}

template<class TYPE>
template<int DIM>
typename CellData<TYPE>::template ConstView<DIM>
CellData<TYPE>::getConstView(int depth) const
{
   return CellData<TYPE>::ConstView<DIM>(getPointer(depth), getGhostBox());
}

#endif

template<class TYPE>
TYPE&
CellData<TYPE>::operator () (
   const CellIndex& i,
   int depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*d_data)(i, depth);
}

template<class TYPE>
const TYPE&
CellData<TYPE>::operator () (
   const CellIndex& i,
   int depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*d_data)(i, depth);
}

template<class TYPE>
ArrayData<TYPE>&
CellData<TYPE>::getArrayData()
{
   return *d_data;
}

template<class TYPE>
const ArrayData<TYPE>&
CellData<TYPE>::getArrayData() const
{
   return *d_data;
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two cell centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
CellData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*d_data, src);

   const CellData<TYPE>* t_src = dynamic_cast<const CellData<TYPE> *>(&src);
   if (t_src == 0) {
      src.copy2(*this);
   } else {
      const hier::Box box = d_data->getBox() * t_src->d_data->getBox();
      if (!box.empty()) {
         d_data->copy(*(t_src->d_data), box);
      }
   }
}

template<class TYPE>
void
CellData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*d_data, dst);

   CellData<TYPE>* t_dst = CPP_CAST<CellData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != 0);

   const hier::Box box = d_data->getBox() * t_dst->d_data->getBox();
   if (!box.empty()) {
      t_dst->d_data->copy(*d_data, box);
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
CellData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   t_copy->start();
   const CellData<TYPE>* t_src = dynamic_cast<const CellData<TYPE> *>(&src);

   const CellOverlap* t_overlap = dynamic_cast<const CellOverlap *>(&overlap);

   if ((t_src == 0) || (t_overlap == 0)) {
      src.copy2(*this, overlap);
   } else {
      if (t_overlap->getTransformation().getRotation() ==
          hier::Transformation::NO_ROTATE) {
         d_data->copy(*(t_src->d_data),
            t_overlap->getDestinationBoxContainer(),
            t_overlap->getTransformation());
      } else {
         copyWithRotation(*t_src, *t_overlap);
      }
   }
   t_copy->stop();
}

template<class TYPE>
void
CellData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   CellData<TYPE>* t_dst = CPP_CAST<CellData<TYPE> *>(&dst);
   const CellOverlap* t_overlap = CPP_CAST<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != 0);
   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      t_dst->d_data->copy(*d_data,
         t_overlap->getDestinationBoxContainer(),
         t_overlap->getTransformation());
   } else {
      t_dst->copyWithRotation(*this, *t_overlap);
   }
}

template<class TYPE>
void
CellData<TYPE>::copyOnBox(
   const CellData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);
   d_data->copy(src.getArrayData(), box);
}

template<class TYPE>
void
CellData<TYPE>::copyWithRotation(
   const CellData<TYPE>& src,
   const CellOverlap& overlap)
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(src.getDim());
   const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer();
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   hier::Box rotatebox(src.getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(rotatebox * overlap_box);

      if (!copybox.empty()) {
         const int depth = ((getDepth() < src.getDepth()) ?
                            getDepth() : src.getDepth());

         CellData<double>::iterator ciend(CellGeometry::end(copybox));
         for (CellData<double>::iterator ci(CellGeometry::begin(copybox));
              ci != ciend; ++ci) {

            const CellIndex& dst_index = *ci;
            CellIndex src_index(dst_index);
            hier::Transformation::rotateIndex(src_index, back_rotate);
            src_index += back_shift;

            for (int d = 0; d < depth; ++d) {
               (*d_data)(dst_index, d) = (*src.d_data)(src_index, d);
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
CellData<TYPE>::copyDepth(
   int dst_depth,
   const CellData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*d_data, src);

   const hier::Box box = d_data->getBox() * src.d_data->getBox();
   if (!box.empty()) {
      d_data->copyDepth(dst_depth, *(src.d_data), src_depth, box);
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
CellData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
CellData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const CellOverlap* t_overlap = CPP_CAST<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   return d_data->getDataStreamSize(
      t_overlap->getDestinationBoxContainer(),
      t_overlap->getSourceOffset());
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
CellData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const CellOverlap* t_overlap = CPP_CAST<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {
      d_data->packStream(stream, t_overlap->getDestinationBoxContainer(),
         t_overlap->getTransformation());
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
CellData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const CellOverlap& overlap) const
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(getDim());
   const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer();
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   hier::Box rotatebox(getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   const int depth = getDepth();

   const size_t size = depth * overlap_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   int i = 0;
   for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(rotatebox * overlap_box);

      if (!copybox.empty()) {

         for (int d = 0; d < depth; ++d) {
            CellData<double>::iterator ciend(CellGeometry::end(copybox));
            for (CellData<double>::iterator ci(CellGeometry::begin(copybox));
                 ci != ciend; ++ci) {

               CellIndex src_index(*ci);
               hier::Transformation::rotateIndex(src_index, back_rotate);
               src_index += back_shift;

               buffer[i] = (*d_data)(src_index, d);
               ++i;
            }
         }
      }
   }

   stream.pack(&buffer[0], size);
}

template<class TYPE>
void
CellData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const CellOverlap* t_overlap = CPP_CAST<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   d_data->unpackStream(stream, t_overlap->getDestinationBoxContainer(),
      t_overlap->getSourceOffset());
}

/*
 *************************************************************************
 *                                                                       *
 * Add source data to the destination according to overlap descriptor.   *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::sum(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const CellData<TYPE>* t_onode_src =
      dynamic_cast<const CellData<TYPE> *>(&src);

   // NOTE:  We assume this operation is only needed to
   //        copy and add data to another cell data
   //        object.  If we ever need to provide this for node
   //        data or other flavors of the copy operation, we
   //        should refactor the routine similar to the way
   //        the regular copy operations are implemented.
   if (t_onode_src == 0) {
      TBOX_ERROR("CellData<dim>::sum error!\n"
         << "Can copy and add only from CellData<TYPE> "
         << "of the same dim and TYPE.");
   } else {

      const hier::IntVector& src_offset(t_overlap->getSourceOffset());
      const hier::BoxContainer& box_container(
         t_overlap->getDestinationBoxContainer());
      const ArrayData<TYPE>& src_array(*t_onode_src->d_data);
      if (d_data->isInitialized()) {
         d_data->sum(src_array, box_container, src_offset);
      }
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Unpack data from the message stream and add to this cell data         *
 * object using the index space in the overlap descriptor.               *
 *                                                                       *
 *************************************************************************
 */

template<class TYPE>
void CellData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& dst_boxes(
      t_overlap->getDestinationBoxContainer());
   const hier::IntVector& src_offset(t_overlap->getSourceOffset());
   for (hier::BoxContainer::const_iterator dst_box(dst_boxes.begin());
        dst_box != dst_boxes.end(); ++dst_box) {
      const hier::Box intersect(*dst_box * d_data->getBox());
      if (!intersect.empty()) {
         d_data->unpackStreamAndSum(stream, intersect, src_offset);
      }
   }
}

template<class TYPE>
void
CellData<TYPE>::fill(
   const TYPE& t,
   int d)
{

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data->fill(t, d);
}

template<class TYPE>
void
CellData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data->fill(t, box, d);
}

template<class TYPE>
void
CellData<TYPE>::fillAll(
   const TYPE& t)
{
   d_data->fillAll(t);
}

template<class TYPE>
void
CellData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   d_data->fillAll(t, box);
}

#ifdef SAMRAI_HAVE_CONDUIT
template<class TYPE>
void
CellData<TYPE>::putBlueprintField(
   conduit::Node& domain_node,
   const std::string& field_name,
   const std::string& topology_name,
   int depth)
{
   size_t data_size = getGhostBox().size();
   domain_node["fields"][field_name]["values"].set_external(
      getPointer(depth), data_size);
   domain_node["fields"][field_name]["association"].set_string("element");
   domain_node["fields"][field_name]["type"].set_string("scalar");
   domain_node["fields"][field_name]["topology"].set_string(topology_name);
}
#endif

/*
 *************************************************************************
 *
 * Print cell centered data.  Note:  makes call to specialized print
 * routine in CellDataSpecialized.cpp
 *
 *************************************************************************
 */

template<class TYPE>
void
CellData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   for (int d = 0; d < d_depth; ++d) {
      os << "Array depth = " << d << std::endl;
      print(box, d, os, prec);
   }
}

template<class TYPE>
void
CellData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   os.precision(prec);
   CellIterator iend(CellGeometry::end(box));
   for (CellIterator i(CellGeometry::begin(box)); i != iend; ++i) {
      os << "array" << *i << " = "
         << (*d_data)(*i, depth) << std::endl << std::flush;
      os << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in the d_depth data member to the restart database.  Then tells
 * d_data to read itself in from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
CellData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{

   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_CELLDATA_VERSION");
   if (ver != PDAT_CELLDATA_VERSION) {
      TBOX_ERROR("CellData<TYPE>::getFromRestart error...\n"
         << "Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");

   d_data->getFromRestart(restart_db->getDatabase("d_data"));
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth data member to the
 * restart database.  Then tells d_data to write itself to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
CellData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_CELLDATA_VERSION", PDAT_CELLDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);

   d_data->putToRestart(restart_db->putDatabase("d_data"));
}

#if defined(HAVE_RAJA)
template<int DIM, typename TYPE, typename... Args>
typename CellData<TYPE>::template View<DIM> get_view(CellData<TYPE>& data, Args&&... args)
{
   return data.template getView<DIM>(std::forward<Args>(args)...);
}

template<int DIM, typename TYPE, typename... Args>
typename CellData<TYPE>::template ConstView<DIM> get_const_view(const CellData<TYPE>& data, Args&&... args)
{
   return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
