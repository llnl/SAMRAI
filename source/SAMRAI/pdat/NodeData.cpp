/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated node centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_NodeData_C
#define included_pdat_NodeData_C

#include "SAMRAI/pdat/NodeData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int NodeData<TYPE>::PDAT_NODEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for node data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
NodeData<TYPE>::NodeData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   const hier::Box node = NodeGeometry::toNodeBox(getGhostBox());
   d_data.reset(new ArrayData<TYPE>(node, depth));
}

template<class TYPE>
NodeData<TYPE>::NodeData(
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

   const hier::Box node = NodeGeometry::toNodeBox(getGhostBox());
   d_data.reset(new ArrayData<TYPE>(node, depth, allocator));
}

template<class TYPE>
NodeData<TYPE>::~NodeData()
{
}

template<class TYPE>
int
NodeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
ArrayData<TYPE>&
NodeData<TYPE>::getArrayData()
{
   return *d_data;
}

template<class TYPE>
const ArrayData<TYPE>&
NodeData<TYPE>::getArrayData() const
{
   return *d_data;
}

template<class TYPE>
TYPE *
NodeData<TYPE>::getPointer(
   int depth)
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data->getPointer(depth);
}

template<class TYPE>
const TYPE *
NodeData<TYPE>::getPointer(
   int depth) const
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data->getPointer(depth);
}

#if defined(HAVE_RAJA)
template<class TYPE>
template<int DIM>
typename NodeData<TYPE>::template View<DIM>
NodeData<TYPE>::getView(int depth)
{
   const hier::Box node_box = NodeGeometry::toNodeBox(getGhostBox());
   return NodeData<TYPE>::View<DIM>(getPointer(depth), node_box);
}

template<class TYPE>
template<int DIM>
typename NodeData<TYPE>::template ConstView<DIM>
NodeData<TYPE>::getConstView(int depth) const
{
   const hier::Box node_box = NodeGeometry::toNodeBox(getGhostBox());
   return NodeData<TYPE>::ConstView<DIM>(getPointer(depth), node_box);
}
#endif

template<class TYPE>
TYPE&
NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*d_data)(i, depth);
}

template<class TYPE>
const TYPE&
NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return (*d_data)(i, depth);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two node centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const NodeData<TYPE>* t_src = dynamic_cast<const NodeData<TYPE> *>(&src);
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
NodeData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   NodeData<TYPE>* t_dst = CPP_CAST<NodeData<TYPE> *>(&dst);

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
NodeData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const NodeData<TYPE>* t_src = dynamic_cast<const NodeData<TYPE> *>(&src);
   const NodeOverlap* t_overlap = dynamic_cast<const NodeOverlap *>(&overlap);

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
}

template<class TYPE>
void
NodeData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   NodeData<TYPE>* t_dst = CPP_CAST<NodeData<TYPE> *>(&dst);
   const NodeOverlap* t_overlap = CPP_CAST<const NodeOverlap *>(&overlap);

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
NodeData<TYPE>::copyOnBox(
   const NodeData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);
   const hier::Box node_box = NodeGeometry::toNodeBox(box);
   d_data->copy(src.getArrayData(), node_box);
}

template<class TYPE>
void
NodeData<TYPE>::copyWithRotation(
   const NodeData<TYPE>& src,
   const NodeOverlap& overlap)
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

   hier::Box node_rotatebox(NodeGeometry::toNodeBox(rotatebox));

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   node_rotatebox.getBlockId(),
                                   getBox().getBlockId());

   for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(node_rotatebox * overlap_box);

      if (!copybox.empty()) {
         const int depth = ((getDepth() < src.getDepth()) ?
                            getDepth() : src.getDepth());

         hier::Box::iterator ciend(copybox.end());
         for (hier::Box::iterator ci(copybox.begin()); ci != ciend; ++ci) {

            NodeIndex dst_index(*ci, hier::IntVector::getZero(dim));
            NodeIndex src_index(dst_index);
            NodeGeometry::transform(src_index, back_trans);

            for (int d = 0; d < depth; ++d) {
               (*d_data)(dst_index, d) = (*(src.d_data))(src_index, d);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a node data object to this node data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::copyDepth(
   int dst_depth,
   const NodeData<TYPE>& src,
   int src_depth)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

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
NodeData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
size_t
NodeData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap = CPP_CAST<const NodeOverlap *>(&overlap);

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
NodeData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap = CPP_CAST<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {
      d_data->packStream(stream,
         t_overlap->getDestinationBoxContainer(),
         t_overlap->getTransformation());
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
NodeData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const NodeOverlap& overlap) const
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

   hier::Box node_rotatebox(NodeGeometry::toNodeBox(rotatebox));

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   const int depth = getDepth();

   const size_t size = depth * overlap_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   int i = 0;
   for (hier::BoxContainer::const_iterator bi = overlap_boxes.begin();
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(node_rotatebox * overlap_box);

      if (!copybox.empty()) {

         for (int d = 0; d < depth; ++d) {
            hier::Box::iterator ciend(copybox.end());
            for (hier::Box::iterator ci(copybox.begin()); ci != ciend; ++ci) {

               NodeIndex src_index(*ci, hier::IntVector::getZero(dim));
               NodeGeometry::transform(src_index, back_trans);

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
NodeData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const NodeOverlap* t_overlap = CPP_CAST<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != 0);

   d_data->unpackStream(stream,
      t_overlap->getDestinationBoxContainer(),
      t_overlap->getSourceOffset());
}

template<class TYPE>
void
NodeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data->fill(t, d);
}

template<class TYPE>
void
NodeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data->fill(t, NodeGeometry::toNodeBox(box), d);
}

template<class TYPE>
void
NodeData<TYPE>::fillAll(
   const TYPE& t)
{
   d_data->fillAll(t);
}

template<class TYPE>
void
NodeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   d_data->fillAll(t, NodeGeometry::toNodeBox(box));
}

#ifdef SAMRAI_HAVE_CONDUIT
template<class TYPE>
void
NodeData<TYPE>::putBlueprintField(
   conduit::Node& domain_node,
   const std::string& field_name,
   const std::string& topology_name,
   int depth)
{
   size_t data_size = getGhostBox().size();
   domain_node["fields"][field_name]["values"].set_external(
      getPointer(depth), data_size);
   domain_node["fields"][field_name]["association"].set_string("element");
   domain_node["fields"][field_name]["type"].set_string("vertex");
   domain_node["fields"][field_name]["topology"].set_string(topology_name);
}
#endif


/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a node centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
NodeData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
   TBOX_ASSERT(depth > 0);

   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   const hier::Box node_box = NodeGeometry::toNodeBox(ghost_box);
   return ArrayData<TYPE>::getSizeOfData(node_box, depth);
}

/*
 *************************************************************************
 *
 * Print node-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::print(
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
NodeData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   os.precision(prec);
   NodeIterator iend(NodeGeometry::end(box));
   for (NodeIterator i(NodeGeometry::begin(box)); i != iend; ++i) {
      os << "array" << *i << " = "
         << (*d_data)(*i, depth) << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks to make sure that the class version and restart file
 * version are equal.  If so, reads in d_depth and has d_data
 * retrieve its own data from the restart database.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_NODEDATA_VERSION");
   if (ver != PDAT_NODEDATA_VERSION) {
      TBOX_ERROR("NodeData<TYPE>::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");

   d_data->getFromRestart(restart_db->getDatabase("d_data"));
}

/*
 *************************************************************************
 *
 * Writes out the class version number and d_depth, Then has d_data
 * write its own data to the restart database.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_NODEDATA_VERSION", PDAT_NODEDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);

   d_data->putToRestart(restart_db->putDatabase("d_data"));
}

#if defined(HAVE_RAJA)
template<int DIM, typename TYPE, typename... Args>
typename NodeData<TYPE>::template View<DIM> get_view(NodeData<TYPE>& data, Args&&... args)
{
   return data.template getView<DIM>(std::forward<Args>(args)...);
}

template<int DIM, typename TYPE, typename... Args>
typename NodeData<TYPE>::template ConstView<DIM> get_const_view(const NodeData<TYPE>& data, Args&&... args)
{
   return data.template getConstView<DIM>(std::forward<Args>(args)...);
}
#endif

}
}

#endif
