/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_IndexData_C
#define included_pdat_IndexData_C

#include "SAMRAI/pdat/IndexData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE, class BOX_GEOMETRY>
const int IndexData<TYPE, BOX_GEOMETRY>::PDAT_INDEXDATA_VERSION = 1;

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>::IndexDataNode(
   const hier::Index& index,
   const size_t offset,
   TYPE& t,
   IndexDataNode<TYPE, BOX_GEOMETRY>* n,
   IndexDataNode<TYPE, BOX_GEOMETRY>* p):
   d_index(index),
   d_offset(offset),
   d_item(&t),
   d_next(n),
   d_prev(p)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>::~IndexDataNode()
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::getNode()
{
   return *d_node;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexData<TYPE, BOX_GEOMETRY>& index_data,
   bool begin):
   d_index_data(const_cast<IndexData<TYPE, BOX_GEOMETRY> *>(&index_data)),
   d_node(begin ? d_index_data->d_list_head : 0)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   IndexData<TYPE, BOX_GEOMETRY>* index_data,
   IndexDataNode<TYPE, BOX_GEOMETRY>* node):
   d_index_data(index_data),
   d_node(node)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter):
   d_index_data(iter.d_index_data),
   d_node(iter.d_node)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::~IndexIterator()
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator = (
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
TYPE&
IndexIterator<TYPE, BOX_GEOMETRY>::operator * ()
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const TYPE&
IndexIterator<TYPE, BOX_GEOMETRY>::operator * () const
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const hier::Index&
IndexIterator<TYPE, BOX_GEOMETRY>::getIndex() const
{
   return d_node->d_index;
}

template<class TYPE, class BOX_GEOMETRY>
TYPE *
IndexIterator<TYPE, BOX_GEOMETRY>::operator -> ()
{
   return d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const TYPE *
IndexIterator<TYPE, BOX_GEOMETRY>::operator -> () const
{
   return d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator ++ ()
{
   if (d_node) {
      d_node = d_node->d_next;
   }
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::operator ++ (
   int)
{
   IndexIterator<TYPE, BOX_GEOMETRY> tmp = *this;
   if (d_node) {
      d_node = d_node->d_next;
   }
   return tmp;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator -- ()
{
   if (d_node) {
      d_node = d_node->d_prev;
   }
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::operator -- (
   int)
{
   IndexIterator<TYPE, BOX_GEOMETRY> tmp = *this;
   if (d_node) {
      d_node = d_node->d_prev;
   }
   return tmp;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexIterator<TYPE, BOX_GEOMETRY>::operator == (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node == i.d_node;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexIterator<TYPE, BOX_GEOMETRY>::operator != (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node != i.d_node;
}

/*
 *************************************************************************
 *
 * The constructor for the irregular grid object simply initializes the
 * irregular data list to be null (this is done implicitly).
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
IndexData<TYPE, BOX_GEOMETRY>::IndexData(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_dim(box.getDim()),
   d_data(hier::PatchData::getGhostBox().size()),
   d_list_head(0),
   d_list_tail(0),
   d_number_items(0)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);
}

template<class TYPE, class BOX_GEOMETRY>
IndexData<TYPE, BOX_GEOMETRY>::~IndexData()
{
   removeAllItems();
}

/*
 *************************************************************************
 *
 * Copy into dst where src overlaps on interiors.
 *
 *************************************************************************
 */
template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const IndexData<TYPE, BOX_GEOMETRY>* t_src =
      CPP_CAST<const IndexData<TYPE, BOX_GEOMETRY> *>(&src);

   TBOX_ASSERT(t_src != 0);

   const hier::Box& src_ghost_box = t_src->getGhostBox();
   removeInsideBox(src_ghost_box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*t_src, false);
   for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*t_src, true);
        s != send;
        ++s) {
      if (getGhostBox().contains(s.getNode().d_index)) {
         appendItem(s.getNode().d_index, *(s.getNode().d_item));
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   dst.copy(*this);
}

/*
 *************************************************************************
 *
 * Copy data from the source into the destination according to the
 * overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const IndexData<TYPE, BOX_GEOMETRY>* t_src =
      CPP_CAST<const IndexData<TYPE, BOX_GEOMETRY> *>(&src);
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);

   TBOX_ASSERT(t_src != 0);
   TBOX_ASSERT(t_overlap != 0);

   const hier::IntVector& src_offset(t_overlap->getSourceOffset());
   const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer();
   const hier::Box& src_ghost_box = t_src->getGhostBox();

   for (hier::BoxContainer::const_iterator b = box_list.begin();
        b != box_list.end(); ++b) {
      const hier::Box& dst_box = *b;
      const hier::Box src_box(hier::Box::shift(*b, -src_offset));
      removeInsideBox(dst_box);
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*t_src, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*t_src, true);
           s != send;
           ++s) {
         if (src_box.contains(s.getNode().d_index)) {
            TYPE new_item;
            new_item.copySourceItem(
               s.getNode().d_index,
               src_offset,
               *(t_src->d_data[src_ghost_box.offset(s.getNode().d_index)]->
                 d_item));
            appendItem(s.getNode().d_index + src_offset, new_item);
         }
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);

   dst.copy(*this, overlap);
}

/*
 *************************************************************************
 *
 * Calculate the buffer space needed to pack/unpack messages on the box
 * region using the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::canEstimateStreamSizeFromBox() const
{
   return false;
}

template<class TYPE, class BOX_GEOMETRY>
size_t
IndexData<TYPE, BOX_GEOMETRY>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != 0);

   size_t bytes = 0;
   int num_items = 0;
   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*b, -(t_overlap->getSourceOffset()));
      hier::Box::iterator indexend(box.end());
      for (hier::Box::iterator index(box.begin());
           index != indexend; ++index) {
         TYPE* item = getItem(*index);
         if (item) {
            ++num_items;
            bytes += item->getDataStreamSize();
         }
      }
   }
   const size_t index_size = d_dim.getValue() * tbox::MessageStream::getSizeof<int>();
   bytes += (num_items * index_size + tbox::MessageStream::getSizeof<int>());
   return bytes;
}

/*
 *************************************************************************
 *
 * Pack/unpack data into/out of the message streams using the index
 * space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   int num_items = 0;
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*b, -(t_overlap->getSourceOffset()));
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*this, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*this, true);
           s != send; ++s) {
         if (box.contains(s.getNode().d_index)) {
            ++num_items;
         }
      }
   }

   stream << num_items;

   for (hier::BoxContainer::const_iterator c = boxes.begin();
        c != boxes.end(); ++c) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*c, -(t_overlap->getSourceOffset()));
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator tend(*this, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator t(*this, true);
           t != tend; ++t) {
         if (box.contains(t.getNode().d_index)) {
            TYPE* item = &(*t);
            TBOX_ASSERT(item != 0);

            int index_buf[SAMRAI::MAX_DIM_VAL];
            for (int i = 0; i < d_dim.getValue(); ++i) {
               index_buf[i] = t.getNode().d_index(i);
            }
            stream.pack(index_buf, d_dim.getValue());
            item->packStream(stream);
         }
      }
   }

}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != 0);

   int num_items;
   stream >> num_items;

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      removeInsideBox(*b);
   }

   int i;
   TYPE* items = 0;
   if (num_items > 0) {
      items = new TYPE[num_items];
   }
   for (i = 0; i < num_items; ++i) {
      int index_buf[SAMRAI::MAX_DIM_VAL];
      stream.unpack(index_buf, d_dim.getValue());
      hier::Index index(d_dim);
      for (int j = 0; j < d_dim.getValue(); ++j) {
         index(j) = index_buf[j];
      }
      (items + i)->unpackStream(stream, t_overlap->getSourceOffset());
      addItem(index + (t_overlap->getSourceOffset()), items[i]);
   }
   if (items) {
      delete[] items;
   }
}

/*
 *************************************************************************
 *
 * List manipulation stuff.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != 0);

   *new_item = item;
   addItemToList(index, offset, *new_item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   appendItemToList(index, offset, *item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != 0);

   *new_item = item;
   addItemToList(index, offset, *new_item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   addItemToList(index, offset, *item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAddItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != 0);

   *new_item = item;

   if (node == 0) {

      addItemToList(index, offset, *new_item);

   } else {
      delete node->d_item;

      node->d_item = new_item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAddItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   if (node == 0) {

      addItemToList(index, offset, *item);

   } else {

      delete node->d_item;

      node->d_item = item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAppendItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != 0);

   *new_item = item;

   if (node == 0) {

      appendItemToList(index, offset, *new_item);

   } else {
      delete node->d_item;

      node->d_item = new_item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAppendItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   if (node == 0) {

      appendItemToList(index, offset, *item);

   } else {

      delete node->d_item;

      node->d_item = item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeItem(
   const hier::Index& index)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   size_t offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   removeItem(offset);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeItem(
   const size_t offset)
{
   TBOX_ASSERT(offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TBOX_ASSERT(node);

   removeNodeFromList(node);

   delete node->d_item;
   delete node;

   d_data[offset] = 0;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItemToList(
   const hier::Index& index,
   const size_t offset,
   TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   IndexDataNode<TYPE, BOX_GEOMETRY>* new_node =
      new IndexDataNode<TYPE, BOX_GEOMETRY>(index,
                                            offset,
                                            item,
                                            d_list_head,
                                            0);

   if (d_list_head) {
      d_list_head->d_prev = new_node;
   }

   d_list_head = new_node;

   if (!d_list_tail) {
      d_list_tail = new_node;
   }

   d_data[offset] = new_node;

   ++d_number_items;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItemToList(
   const hier::Index& index,
   const size_t offset,
   TYPE& item)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   IndexDataNode<TYPE, BOX_GEOMETRY>* new_node =
      new IndexDataNode<TYPE, BOX_GEOMETRY>(index,
                                            offset,
                                            item,
                                            0,
                                            d_list_tail);

   if (d_list_tail) {
      d_list_tail->d_next = new_node;
   }

   d_list_tail = new_node;

   if (!d_list_head) {
      d_list_head = new_node;
   }

   d_data[offset] = new_node;

   ++d_number_items;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeNodeFromList(
   IndexDataNode<TYPE, BOX_GEOMETRY>* node)
{
   if ((d_list_head == node) && (d_list_tail == node)) {
      d_list_head = d_list_tail = 0;

   } else if (d_list_head == node) {
      d_list_head = node->d_next;
      node->d_next->d_prev = 0;

   } else if (d_list_tail == node) {
      d_list_tail = node->d_prev;
      node->d_prev->d_next = 0;

   } else {
      node->d_next->d_prev = node->d_prev;
      node->d_prev->d_next = node->d_next;
   }

   d_data[node->d_offset] = 0;

   --d_number_items;
}

template<class TYPE, class BOX_GEOMETRY>
size_t
IndexData<TYPE, BOX_GEOMETRY>::getNumberOfItems() const
{
   return d_number_items;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeInsideBox(
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator l(*this, true);
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator lend(*this, false);

   while (l != lend) {
      if (box.contains(l.getNode().d_index)) {
         hier::Index index(l.getNode().d_index);
         ++l;
         removeItem(index);
      } else {
         ++l;
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeOutsideBox(
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator l(*this, true);
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator lend(*this, false);

   while (l != lend) {
      if (!box.contains(l.getNode().d_index)) {
         hier::Index index(l.getNode().d_index);
         ++l;
         removeItem(index);
      } else {
         ++l;
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeGhostItems()
{
   removeOutsideBox(hier::PatchData::getBox());
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeAllItems()
{
   removeInsideBox(hier::PatchData::getGhostBox());
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::isElement(
   const hier::Index& index) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   return d_data[hier::PatchData::getGhostBox().offset(index)] != 0;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::isElement(
   size_t offset) const
{
   return d_data[offset] != 0;
}

/*
 *************************************************************************
 *
 * Just checks to make sure that the class version is the same
 * as the restart file version number.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   int ver = restart_db->getInteger("PDAT_INDEXDATA_VERSION");
   if (ver != PDAT_INDEXDATA_VERSION) {
      TBOX_ERROR("IndexData::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   int item_count = 0;
   bool item_found = true;

   do {
      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(
            item_count,
            6);

      if (restart_db->isDatabase(index_keyword)) {

         std::shared_ptr<tbox::Database> item_db(
            restart_db->getDatabase(index_keyword));

         std::vector<int> index_array =
            item_db->getIntegerVector(index_keyword);
         hier::Index index(d_dim);
         for (int j = 0; j < d_dim.getValue(); ++j) {
            index(j) = index_array[j];
         }

         TYPE item;
         item.getFromRestart(item_db);

         appendItem(index, item);

      } else {
         item_found = false;
      }

      ++item_count;

   } while (item_found);

}

/*
 *************************************************************************
 *
 * Just writes out the class version number to the restart database.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   restart_db->putInteger("PDAT_INDEXDATA_VERSION", PDAT_INDEXDATA_VERSION);

   int item_count = 0;
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*this, false);
   for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*this, true);
        s != send; ++s) {

      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(
            item_count,
            6);
      hier::Index index = s.getNode().d_index;
      std::vector<int> index_array(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); ++i) {
         index_array[i] = index(i);
      }

      std::shared_ptr<tbox::Database> item_db(
         restart_db->putDatabase(index_keyword));

      item_db->putIntegerArray(index_keyword,
         &index_array[0],
         static_cast<int>(index_array.size()));

      TYPE* item = getItem(index);

      item->putToRestart(item_db);

      ++item_count;
   }
}

template<class TYPE, class BOX_GEOMETRY>
TYPE *
IndexData<TYPE, BOX_GEOMETRY>::getItem(
   const hier::Index& index) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, index);

   TYPE* item;
   if (!isElement(index)) {
      item = 0;
   } else {
      item = d_data[hier::PatchData::getGhostBox().offset(index)]->d_item;
   }

   return item;
}

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
