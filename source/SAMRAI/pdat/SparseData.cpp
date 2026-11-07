/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   pdat
 *
 ************************************************************************/
#ifndef included_pdat_SparseData_C
#define included_pdat_SparseData_C

#include "SAMRAI/pdat/SparseData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/PIO.h"

#include <stdexcept>
#include <algorithm>
#include <cctype>

#ifdef __GNUC__
#if __GNUC__ == 4 && __GNUC_MINOR__ == 1
#define GNUC_VERSION_412 412
#endif
#endif

namespace SAMRAI {
namespace pdat  {

/**********************************************************************
 * PDAT_SPARSEDATA_VERSION
 *********************************************************************/
template<typename BOX_GEOMETRY>
const int
SparseData<BOX_GEOMETRY>::PDAT_SPARSEDATA_VERSION = 2;

/**********************************************************************
 * INVALID_ID
 *********************************************************************/
template<typename BOX_GEOMETRY>
const int
SparseData<BOX_GEOMETRY>::INVALID_ID = 1;

/**********************************************************************
* Constructor
**********************************************************************/
template<typename BOX_GEOMETRY>
SparseData<BOX_GEOMETRY>::Attributes::Attributes(
   const int dsize,
   const int isize):
   d_dbl_attrs(dsize),
   d_int_attrs(isize)
{
}

/**********************************************************************
* Copy ctor
**********************************************************************/
template<typename BOX_GEOMETRY>
SparseData<BOX_GEOMETRY>::Attributes::Attributes(
   const Attributes& other):
   d_dbl_attrs(other.d_dbl_attrs),
   d_int_attrs(other.d_int_attrs)
{
}

/**********************************************************************
* Dtor
**********************************************************************/
template<typename BOX_GEOMETRY>
SparseData<BOX_GEOMETRY>::Attributes::~Attributes()
{
}

/**********************************************************************
* Assignment
**********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::Attributes&
SparseData<BOX_GEOMETRY>::Attributes::operator = (
   const typename SparseData<BOX_GEOMETRY>::Attributes& rhs)
{
   if (this != &rhs) {
      d_dbl_attrs = rhs.d_dbl_attrs;
      d_int_attrs = rhs.d_int_attrs;
   }
   return *this;
}

/**********************************************************************
 * add(dvals, ivals)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::Attributes::add(
   const double* dvals, const int* ivals)
{
   std::copy(dvals, dvals + d_dbl_attrs.size(), d_dbl_attrs.begin());
   std::copy(ivals, ivals + d_int_attrs.size(), d_int_attrs.begin());
}

template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::Attributes::add(
   const std::vector<double>&dvals, const std::vector<int>&ivals)
{
   std::copy(dvals.begin(), dvals.end(), d_dbl_attrs.begin());
   std::copy(ivals.begin(), ivals.end(), d_int_attrs.begin());
}

/**********************************************************************
* non-modifying operations
**********************************************************************/
template<typename BOX_GEOMETRY>
const double *
SparseData<BOX_GEOMETRY>::Attributes::getDoubleAttributes() const {
   return &d_dbl_attrs[0];
}

template<typename BOX_GEOMETRY>
const int *
SparseData<BOX_GEOMETRY>::Attributes::getIntAttributes() const {
   return &d_int_attrs[0];
}

/**********************************************************************
 * access operators
 *********************************************************************/
template<typename BOX_GEOMETRY>
double&
SparseData<BOX_GEOMETRY>::Attributes::operator [] (const DoubleAttributeId& id)
{
   return d_dbl_attrs[id()];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
const double&
SparseData<BOX_GEOMETRY>::Attributes::operator [] (const DoubleAttributeId& id) const
{
   return d_dbl_attrs[id()];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
int&
SparseData<BOX_GEOMETRY>::Attributes::operator [] (const IntegerAttributeId& id)
{
   return d_int_attrs[id()];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
const int&
SparseData<BOX_GEOMETRY>::Attributes::operator [] (const IntegerAttributeId& id) const
{
   return d_int_attrs[id()];
}

/**********************************************************************
*
**********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::Attributes::operator == (
   const Attributes& rhs) const
{
   return d_dbl_attrs == rhs.d_dbl_attrs &&
          d_int_attrs == rhs.d_int_attrs;
}

/**********************************************************************
 * Attribute Class implementations.
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::Attributes::printAttributes(
   std::ostream & out) const
{
   dbl_iterator diter(d_dbl_attrs.begin()), diterend(d_dbl_attrs.end());
   int_iterator iiter(d_int_attrs.begin()), iiterend(d_int_attrs.end());

   out << "Double Attributes ( ";
   for ( ; diter != diterend; ++diter) {
      out.precision(6);
      out << *diter << " ";
   }
   out << ")" << std::endl;

   out << "Integer Attributes ( ";
   for ( ; iiter != iiterend; ++iiter) {
      out << *iiter << " ";
   }
   out << ")" << std::endl;
}

/**********************************************************************
 * hasher
 *********************************************************************/
template<typename BOX_GEOMETRY>
std::size_t
SparseData<BOX_GEOMETRY>::index_hash::operator () (
   const hier::Index& index) const
{
   std::size_t seed = 0;
   int dim = index.getDim().getValue();
   SparseData<BOX_GEOMETRY>::hash_combine(seed, dim);
   for (int i = 0; i < dim; ++i) {
      SparseData<BOX_GEOMETRY>::hash_combine(seed, index[i]);
      SparseData<BOX_GEOMETRY>::hash_combine(seed, index[i]);
   }
   return seed;
}
/**********************************************************************
 * SparseData
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseData<BOX_GEOMETRY>::SparseData(
   const hier::Box& box,
   const hier::IntVector& ghosts,
   const std::vector<std::string>& dbl_names,
   const std::vector<std::string>& int_names):
   hier::PatchData(box, ghosts),
   d_dim(box.getDim()),
   d_dbl_attr_size(static_cast<int>(dbl_names.size())),
   d_int_attr_size(static_cast<int>(int_names.size()))
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, ghosts);

   std::vector<std::string>::const_iterator name_iter;
   int val(0);
   std::string tmp;
   for (name_iter = dbl_names.begin(); name_iter != dbl_names.end();
        ++name_iter, ++val) {
      tmp = *name_iter;
      to_lower(tmp);
      d_dbl_names.insert(std::make_pair(tmp, DoubleAttributeId(val)));
   }

   val = 0;
   for (name_iter = int_names.begin(); name_iter != int_names.end();
        ++name_iter, ++val) {
      tmp = *name_iter;
      to_lower(tmp);
      d_int_names.insert(std::make_pair(tmp, IntegerAttributeId(val)));
   }

}

/**********************************************************************
 * d'tor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseData<BOX_GEOMETRY>::~SparseData()
{
}

/**********************************************************************
 * copy(src)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::copy(
   const hier::PatchData& src)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const SparseData<BOX_GEOMETRY>* tmp_src =
      CPP_CAST<const SparseData<BOX_GEOMETRY> *>(&src);

   TBOX_ASSERT(tmp_src != 0);
   const hier::Box& src_ghost_box = tmp_src->getGhostBox();
   _removeInsideBox(src_ghost_box);

   typename IndexMap::const_iterator src_index_map_iterator =
      tmp_src->d_index_to_attribute_map.begin();

   for ( ; src_index_map_iterator != tmp_src->d_index_to_attribute_map.end();
         ++src_index_map_iterator) {

      if (getGhostBox().contains(src_index_map_iterator->first)) {
         _add(src_index_map_iterator);
      }

   }
}

/**********************************************************************
 * copy2(dst)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::copy2(
   hier::PatchData& dst) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);
   dst.copy(*this);
}

/**********************************************************************
 * copy(src, overlap)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, src);

   const SparseData<BOX_GEOMETRY>* tmp_src =
      CPP_CAST<const SparseData<BOX_GEOMETRY> *>(&src);

   const typename BOX_GEOMETRY::Overlap * tmp_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);

   TBOX_ASSERT(tmp_src != 0);
   TBOX_ASSERT(tmp_overlap != 0);

   const hier::IntVector& src_offset(tmp_overlap->getSourceOffset());
   const hier::BoxContainer& box_list = tmp_overlap->getDestinationBoxContainer();
   const hier::Box& src_ghost_box = tmp_src->getGhostBox();

   for (hier::BoxContainer::const_iterator overlap_box = box_list.begin();
        overlap_box != box_list.end(); ++overlap_box) {

      const hier::Box& dst_box = *overlap_box;
      const hier::Box src_box(hier::Box::shift(*overlap_box, -src_offset));
      _removeInsideBox(dst_box);

      typename IndexMap::const_iterator src_index_map_iter =
         tmp_src->d_index_to_attribute_map.begin();

      for ( ; src_index_map_iter != tmp_src->d_index_to_attribute_map.end();
            ++src_index_map_iter) {

         if (src_ghost_box.contains(src_index_map_iter->first)) {

            hier::Index idx = src_index_map_iter->first + src_offset;
            d_index_to_attribute_map.insert(
               std::make_pair(idx, src_index_map_iter->second));

         } // if (src_ghost_box.contains(...
      } // for (; src_index_map_iter != ...
   } // for (hier::BoxContainer::const_iterator overlap_box = ...
}

/**********************************************************************
 * copy2(dst, overlap)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, dst);
   dst.copy(*this, overlap);
}

/**********************************************************************
 * canEstimateStreamSizeFromBox()
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::canEstimateStreamSizeFromBox() const
{
   return false;
}

/**********************************************************************
 * getDataStreamSize(overlap)
 *********************************************************************/
template<typename BOX_GEOMETRY>
size_t
SparseData<BOX_GEOMETRY>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * tmp_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);

   TBOX_ASSERT(tmp_overlap != 0);

   size_t bytes = 0;
   int num_items = 0;
   int num_attributes = 0;
   //int exc_count = 0;
   const hier::BoxContainer& boxes = tmp_overlap->getDestinationBoxContainer();

   // first count up the number of items that we'll need to deal
   // with
   for (hier::BoxContainer::const_iterator overlap_box = boxes.begin();
        overlap_box != boxes.end(); ++overlap_box) {

      const hier::Box& box = hier::PatchData::getBox()
         * hier::Box::shift(*overlap_box, -(tmp_overlap->getSourceOffset()));

      typename IndexMap::const_iterator iter = d_index_to_attribute_map.begin();
      typename IndexMap::const_iterator iend = d_index_to_attribute_map.end();

      for ( ; iter != iend; ++iter) {
         if (box.contains(iter->first)) {
            ++num_items;
            num_attributes += static_cast<int>(iter->second.size());
         }
      }
   }

   // an int for the number of items.  This value is always packed
   // into the MessageStream, even if it is zero.
   bytes += tbox::MessageStream::getSizeof<int>();

   if (num_items > 0) {

      // an int for the each of the attribute list sizes.
      bytes += num_items * tbox::MessageStream::getSizeof<int>();

      // an int for the number of double attribute names
      bytes += tbox::MessageStream::getSizeof<int>();

      typename DoubleAttrNameMap::const_iterator dnames =
         d_dbl_names.begin();

      for ( ; dnames != d_dbl_names.end(); ++dnames) {

         // two ints to store the key name size and the value (id)
         bytes += tbox::MessageStream::getSizeof<int>() * 2;

         // and a char each for the actual key
         bytes += (static_cast<int>(dnames->first.size())
                   * tbox::MessageStream::getSizeof<char>());

      } // for (; dname ....

      // and do the same for the integer attribute names
      bytes += tbox::MessageStream::getSizeof<int>();
      typename IntAttrNameMap::const_iterator inames = d_int_names.begin();

      for ( ; inames != d_int_names.end(); ++inames) {

         // two ints to store the key name size and the value (id)
         bytes += tbox::MessageStream::getSizeof<int>() * 2;

         // and a char eachy the actual key
         bytes += (static_cast<int>(inames->first.size())
                   * tbox::MessageStream::getSizeof<char>());

      } // for ( ; inames ...

      // record the size of the Indexes
      const size_t index_size =
         d_dim.getValue() * tbox::MessageStream::getSizeof<int>();
      bytes += (num_items * index_size + tbox::MessageStream::getSizeof<int>());

      // calculate the size of the attributes values
      bytes += num_attributes
         * d_dbl_attr_size
         * tbox::MessageStream::getSizeof<double>();

      bytes += num_attributes
         * d_int_attr_size
         * tbox::MessageStream::getSizeof<int>();

   } // if (num_item > 0)
   return bytes;
}

/**********************************************************************
 * packStream(stream, overlap)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * tmp_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(tmp_overlap != 0);

   // Calculate the number of matching items
   const hier::BoxContainer& boxes = tmp_overlap->getDestinationBoxContainer();

   int num_items = 0;

   for (hier::BoxContainer::const_iterator overlap_box = boxes.begin();
        overlap_box != boxes.end(); ++overlap_box) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*overlap_box, -(tmp_overlap->getSourceOffset()));

      typename IndexMap::const_iterator iter = d_index_to_attribute_map.begin();
      typename IndexMap::const_iterator iend = d_index_to_attribute_map.end();

      for ( ; iter != iend; ++iter) {
         if (box.contains(iter->first)) {
            ++num_items;
         }
      }
   }

   // pack number of total matching items
   stream << num_items;

   int key_size = 0;
   int num_dbl_attrs = static_cast<int>(d_dbl_names.size());
   int num_int_attrs = static_cast<int>(d_int_names.size());

   if (num_items > 0) {
      // pack the double keys first
      // start with the number of attributes
      stream << num_dbl_attrs;
      typename DoubleAttrNameMap::const_iterator dbl_name_iter =
         d_dbl_names.begin();

      // then pack the key-value pairs for the double attribute names.
      for ( ; dbl_name_iter != d_dbl_names.end(); ++dbl_name_iter) {

         // key size plus the key itself and it's mapped value
         key_size = static_cast<int>(dbl_name_iter->first.size());
         stream << key_size;
         std::string key(dbl_name_iter->first);
         for (int i = 0; i < key_size; ++i) {
            stream.pack<char>(&key[i], 1);
         }
         stream << dbl_name_iter->second();
      }

      // pack the int keys next
      // start with the number of integer attributes
      stream << num_int_attrs;
      typename IntAttrNameMap::const_iterator int_name_iter =
         d_int_names.begin();

      // then pack the key-value pairs for the integer attribute names
      for ( ; int_name_iter != d_int_names.end(); ++int_name_iter) {

         // key size plus the key itself and it's mapped value
         key_size = static_cast<int>(int_name_iter->first.size());
         stream << key_size;
         std::string key(int_name_iter->first);
         for (int i = 0; i < key_size; ++i) {
            stream.pack<char>(&key[i], 1);
         }
         stream << int_name_iter->second();
      }
   }

   // pack the individual items
   for (hier::BoxContainer::const_iterator overlap_box = boxes.begin();
        overlap_box != boxes.end(); ++overlap_box) {

      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*overlap_box, -(tmp_overlap->getSourceOffset()));

      typename IndexMap::const_iterator index_map_iter =
         d_index_to_attribute_map.begin();

      typename IndexMap::const_iterator index_map_iend =
         d_index_to_attribute_map.end();

      for ( ; index_map_iter != index_map_iend; ++index_map_iter) {

         if (box.contains(index_map_iter->first)) {

            // first pack the Index
            std::unique_ptr<int[]> index_buf = std::make_unique<int[]>(d_dim.getValue());

            for (int i = 0; i < d_dim.getValue(); ++i) {
               index_buf[i] = index_map_iter->first(i);
            }

            stream.pack<int>(index_buf.get(), d_dim.getValue());

            // pack the number of attributes
            int num_list_items =
               static_cast<int>(index_map_iter->second.size());

            stream << num_list_items;

            //typename std::list<Attributes>::const_iterator list_iter =
            typename SparseData<BOX_GEOMETRY>::AttributeList::const_iterator
            list_iter = index_map_iter->second.begin();
            //typename std::list<Attributes>::const_iterator list_end =
            typename SparseData<BOX_GEOMETRY>::AttributeList::const_iterator
            list_end = index_map_iter->second.end();

            for ( ; list_iter != list_end; ++list_iter) {
               Attributes tmp = *list_iter;
               for (int i = 0; i < d_dbl_attr_size; ++i) {
                  stream.pack<double>(&tmp[DoubleAttributeId(i)], 1);
               }
               for (int i = 0; i < d_int_attr_size; ++i) {
                  stream.pack<int>(&tmp[IntegerAttributeId(i)], 1);
               }
            }
         } //  if (box.contains(...
      } // for (; index_map_iter
   } // for (hier::BoxContainer::const_iterator overlap_box = ...
}

/**********************************************************************
 * unpackStream(stream, overlap)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const typename BOX_GEOMETRY::Overlap * tmp_overlap =
      CPP_CAST<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(tmp_overlap != 0);

   int num_items;
   // unpack total number of items
   stream >> num_items;

   int num_dbl_attrs = 0;
   int num_int_attrs = 0;
   //unpack the keys
   d_dbl_names.clear();

   if (num_items > 0) {
      int key_size = 0;
      int value = 0;
      // double keys first starting with the number of attribute names
      stream >> num_dbl_attrs;

      // then unpack the actual key-value pairs for the dbl attribute names
      for (int i = 0; i < num_dbl_attrs; ++i) {
         stream >> key_size;

         // unpack<char> wants to add a delete character (^?) when
         // unpacking more than one character at a time, so we'll unpack a
         // single character at a time.
         char c;
         std::string tmp;
         for (int j = 0; j < key_size; ++j) {
            stream.unpack<char>(&c, 1);
            tmp += c;
         }
         stream >> value;
         d_dbl_names.insert(std::make_pair(tmp,
               DoubleAttributeId(value)));
      }
      // unpack the int keys next, starting with the total number of attrs
      stream >> num_int_attrs;

      // then unpack the key-value pairs for the integer attribute names.
      d_int_names.clear();
      for (int i = 0; i < num_int_attrs; ++i) {
         stream >> key_size;
         char c;
         std::string tmp;
         for (int j = 0; j < key_size; ++j) {
            stream.unpack<char>(&c, 1);
            tmp += c;
         }
         stream >> value;
         d_int_names.insert(std::make_pair(tmp,
               IntegerAttributeId(value)));
      }

   }

   const hier::BoxContainer& boxes = tmp_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator overlap_box = boxes.begin();
        overlap_box != boxes.end(); ++overlap_box) {

      _removeInsideBox(*overlap_box);
   }

   // finally unpack the individual items.
   for (int i = 0; i < num_items; ++i) {

      int num_attrs = 0;
      // Unpack the Index
      std::unique_ptr<int[]> index_buf = std::make_unique<int[]>(d_dim.getValue());
      stream.unpack<int>(index_buf.get(), d_dim.getValue());

      hier::Index index(d_dim);
      for (int j = 0; j < d_dim.getValue(); ++j) {
         index(j) = index_buf[j];
      }

      iterator map_iter = registerIndex(index);

      // unpack the number of attributes
      stream >> num_attrs;
      std::unique_ptr<double[]> dvals = std::make_unique<double[]>(d_dbl_attr_size);
      std::unique_ptr<int[]> ivals = std::make_unique<int[]>(d_int_attr_size);
      for (int count = 0; count < num_attrs; ++count) {
         stream.unpack<double>(dvals.get(), d_dbl_attr_size);
         stream.unpack<int>(ivals.get(), d_int_attr_size);
         map_iter.insert(dvals.get(), ivals.get());
      }
   }
}

/**********************************************************************
 * getFromRestart(database)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::getFromRestart(restart_db);

   // get and check the version
   int ver = restart_db->getInteger("PDAT_SPARSEDATA_VERSION");
   if (ver != PDAT_SPARSEDATA_VERSION) {
      TBOX_ERROR("SparseData::getFromRestart error...\n"
         << " : Restart file version is "
         << "different than class version" << std::endl);
   }

   // Number of elements in this SparseData object
   int count = restart_db->getInteger("sparse_data_count");

   // number of double node item attributes
   TBOX_ASSERT(restart_db->getInteger("dbl_attr_item_count") ==
      d_dbl_attr_size);

   // get the registered double keys and their associated id.
   std::string* keys = 0;
   int* ids = 0;
   if (d_dbl_attr_size > 0) {
      keys = new std::string[d_dbl_attr_size];
      ids = new int[d_dbl_attr_size];
   }
   restart_db->getStringArray("sparse_data_dbl_keys", keys, d_dbl_attr_size);
   restart_db->getIntegerArray("sparse_data_dbl_ids", ids, d_dbl_attr_size);

   d_dbl_names.clear();
   for (int i = 0; i < d_dbl_attr_size; ++i) {
      d_dbl_names.insert(
         std::make_pair(keys[i], DoubleAttributeId(ids[i])));
   }

   if (d_dbl_attr_size > 0) {
      delete[] keys;
      keys = 0;
      delete[] ids;
      ids = 0;
   }

   // number of double node item attributes
   TBOX_ASSERT(restart_db->getInteger("int_attr_item_count") ==
      d_int_attr_size);

   // get the registered integer keys and their associated id.
   if (d_int_attr_size > 0) {
      keys = new std::string[d_int_attr_size];
      ids = new int[d_int_attr_size];
   }
   restart_db->getStringArray("sparse_data_int_keys", keys, d_int_attr_size);
   restart_db->getIntegerArray("sparse_data_int_ids", ids, d_int_attr_size);

   d_int_names.clear();
   for (int i = 0; i < d_int_attr_size; ++i) {
      d_int_names.insert(
         std::make_pair(keys[i], IntegerAttributeId(ids[i])));
   }

   if (d_int_attr_size > 0) {
      delete[] keys;
      keys = 0;
      delete[] ids;
      ids = 0;
   }

   // get the data for each node in this sparse data object
   for (int curr_item = 0; curr_item < count; ++curr_item) {

      std::string index_keyword =
         "attr_index_data_" + tbox::Utilities::intToString(curr_item, 6);

      // get the next item
      if (restart_db->isDatabase(index_keyword)) {
         std::shared_ptr<tbox::Database> item_db(
            restart_db->getDatabase(index_keyword));

         // unpack the index
         std::vector<int> index_array =
            item_db->getIntegerVector(index_keyword);
         hier::Index index(d_dim);
         for (int j = 0; j < d_dim.getValue(); ++j) {
            index(j) = index_array[j];
         }

         // register the new Index so that we can add the attributes
         // to its list.
         iterator new_item = registerIndex(index);

         // get the list size.
         std::string list_size_keyword = "attr_list_size_"
            + tbox::Utilities::intToString(curr_item, 6);
         int list_size = item_db->getInteger(list_size_keyword);

         // all the values for each list are packed together
         std::string dvalues_keyword = "attr_dbl_values_"
            + tbox::Utilities::intToString(curr_item, 6);
         std::string ivalues_keyword = "attr_int_values_"
            + tbox::Utilities::intToString(curr_item, 6);

         int dbl_ary_size = d_dbl_attr_size * list_size;
	 std::unique_ptr<double[]> dvalues =
            std::make_unique<double[]>(static_cast<unsigned int>(dbl_ary_size));

         int int_ary_size = d_int_attr_size * list_size;
         std::unique_ptr<int[]> ivalues = 
            std::make_unique<int[]>(static_cast<unsigned int>(int_ary_size));

         item_db->getDoubleArray(dvalues_keyword, dvalues.get(),
            (dbl_ary_size));
         item_db->getIntegerArray(ivalues_keyword, ivalues.get(),
            (int_ary_size));

         int doffset(0), ioffset(0);
         for (int curr_list = 0; curr_list < list_size; ++curr_list) {

            // unpack the double attributes
            std::vector<double> dbl_ary(d_dbl_attr_size);

            for (int i = doffset; i < doffset + d_dbl_attr_size; ++i) {
               dbl_ary[i - doffset] = dvalues[i];
            }

            // unpack the integer attributes
            std::vector<int> int_ary(d_int_attr_size);

            for (int k = ioffset; k < ioffset + d_int_attr_size; ++k) {
               int_ary[k - ioffset] = ivalues[k];
            }

            new_item.insert(dbl_ary, int_ary);
            doffset += d_dbl_attr_size;
            ioffset += d_int_attr_size;
         } // for (int curr_list ...
      } // if (restart_db->isDatabase(...
      else {
         TBOX_ERROR("SparseData::getFromRestart error...\n"
            << " : Restart database missing data for attribute index "
            << index_keyword << std::endl);
      }
   } // for (int curr_item = ...
}

/**********************************************************************
 * putToRestart(database)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   hier::PatchData::putToRestart(restart_db);

   // record the version
   restart_db->putInteger("PDAT_SPARSEDATA_VERSION", PDAT_SPARSEDATA_VERSION);

   // record the number of sparse data elements
   restart_db->putInteger("sparse_data_count",
      static_cast<int>(d_index_to_attribute_map.size()));
   restart_db->putInteger("dbl_attr_item_count", d_dbl_attr_size);

   // record the keys for the attributes
   std::string* keys = 0;
   int* ids = 0;
   if (d_dbl_attr_size > 0) {
      keys = new std::string[d_dbl_attr_size];
      ids = new int[d_dbl_attr_size];
   }
   typename DoubleAttrNameMap::const_iterator
   dbl_name_iter = d_dbl_names.begin();

   for (int i = 0; i < d_dbl_attr_size; ++i, ++dbl_name_iter) {
      keys[i] = dbl_name_iter->first;
      ids[i] = dbl_name_iter->second();
   }

   restart_db->putStringArray("sparse_data_dbl_keys", keys, d_dbl_attr_size);
   restart_db->putIntegerArray("sparse_data_dbl_ids", ids, d_dbl_attr_size);

   if (d_dbl_attr_size > 0) {
      delete[] keys;
      keys = 0;
      delete[] ids;
      ids = 0;
   }

   // record the keys for the attributes
   restart_db->putInteger("int_attr_item_count", d_int_attr_size);

   if (d_int_attr_size > 0) {
      keys = new std::string[d_int_attr_size];
      ids = new int[d_int_attr_size];
   }
   typename IntAttrNameMap::const_iterator int_name_iter = d_int_names.begin();

   for (int i = 0; i < d_int_attr_size; ++i, ++int_name_iter) {
      keys[i] = int_name_iter->first;
      ids[i] = int_name_iter->second();
   }

   restart_db->putStringArray("sparse_data_int_keys", keys, d_int_attr_size);
   restart_db->putIntegerArray("sparse_data_int_ids", ids, d_int_attr_size);

   if (d_int_attr_size > 0) {
      delete[] keys;
      keys = 0;
      delete[] ids;
      ids = 0;
   }

   // record the actual data for each element
   int curr_item(0);
   typename IndexMap::iterator index_iter =
      const_cast<IndexMap&>(d_index_to_attribute_map).begin();
   typename IndexMap::iterator index_iter_end =
      const_cast<IndexMap&>(d_index_to_attribute_map).end();

   for ( ; index_iter != index_iter_end; ++index_iter) {

      std::string index_keyword =
         "attr_index_data_" + tbox::Utilities::intToString(curr_item, 6);

      // First deal with the Index
      const hier::Index& index = index_iter->first;
      std::vector<int> index_array(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); ++i) {
         index_array[i] = index(i);
      }

      std::shared_ptr<tbox::Database> item_db(
         restart_db->putDatabase(index_keyword));

      item_db->putIntegerVector(index_keyword, index_array);

      // Next get the node and record the double attribute data
      typename SparseData<BOX_GEOMETRY>::AttributeIterator attributes(
         index_iter->second, index_iter->second.begin());

      typename SparseData<BOX_GEOMETRY>::AttributeIterator attr_end(
         index_iter->second, index_iter->second.end());

      int list_size = static_cast<int>(index_iter->second.size());
      std::string list_size_keyword = "attr_list_size_"
         + tbox::Utilities::intToString(curr_item, 6);
      item_db->putInteger(list_size_keyword, list_size);

      std::unique_ptr<double[]> dvalues =
         std::make_unique<double[]>(static_cast<unsigned int>(d_dbl_attr_size * list_size));
      std::unique_ptr<int[]> ivalues =
         std::make_unique<int[]>(static_cast<unsigned int>(d_int_attr_size * list_size));

      // pack all the data together.
      int doffset(0), ioffset(0);
      for ( ; attributes != attr_end; ++attributes) {

         // Record the double attribute data
         for (int i = doffset; i < doffset + d_dbl_attr_size; ++i) {
            dvalues[i] = attributes[DoubleAttributeId(i - doffset)];
         }

         // Record the integer attribute data
         for (int i = ioffset; i < ioffset + d_int_attr_size; ++i) {
            ivalues[i] = attributes[IntegerAttributeId(i - ioffset)];
         }
         doffset += d_dbl_attr_size;
         ioffset += d_int_attr_size;
      }

      std::string dvalues_keyword = "attr_dbl_values_"
         + tbox::Utilities::intToString(curr_item, 6);
      item_db->putDoubleArray(dvalues_keyword, dvalues.get(),
         (d_dbl_attr_size * list_size));

      std::string ivalues_keyword = "attr_int_values_"
         + tbox::Utilities::intToString(curr_item, 6);
      item_db->putIntegerArray(ivalues_keyword, ivalues.get(),
         (d_int_attr_size * list_size));

      ++curr_item;
   }
}

/**********************************************************************
 * getDblAttributeId(name)
 *********************************************************************/
template<typename BOX_GEOMETRY>
const DoubleAttributeId
SparseData<BOX_GEOMETRY>::getDblAttributeId(
   const std::string& attribute) const
{
   DoubleAttributeId id(-1);
   auto iter = d_dbl_names.find(attribute);
   if (iter != d_dbl_names.end()) {
      id = iter->second;
   }
   return id;
}

/**********************************************************************
 * getIntAttributeId(name)
 *********************************************************************/
template<typename BOX_GEOMETRY>
const IntegerAttributeId
SparseData<BOX_GEOMETRY>::getIntAttributeId(
   const std::string& attribute) const
{
   IntegerAttributeId id(-1);
   auto iter = d_int_names.find(attribute);
   if (iter != d_int_names.end()) {
      id = iter->second;
   }
   return id;
}

/**********************************************************************
 * empty()
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::empty()
{
   return d_index_to_attribute_map.empty();
}

/**********************************************************************
 * registerIndex(index)
 *********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::iterator
SparseData<BOX_GEOMETRY>::registerIndex(
   const hier::Index& index)
{
   std::pair<typename IndexMap::iterator, bool> result =
      d_index_to_attribute_map.insert(std::make_pair(
            index, typename SparseData<BOX_GEOMETRY>::AttributeList()));

   return SparseDataIterator<BOX_GEOMETRY>(*this, result.first);
}

/**********************************************************************
 * remove(index)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::remove(iterator& iterToRemove)
{
   d_index_to_attribute_map.erase(iterToRemove.d_iterator++);
}

/**********************************************************************
 * clear()
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::clear()
{
   d_index_to_attribute_map.clear();
}

/**********************************************************************
 * size()
 *********************************************************************/
template<typename BOX_GEOMETRY>
int
SparseData<BOX_GEOMETRY>::size()
{
   return static_cast<int>(d_index_to_attribute_map.size());
}

/**********************************************************************
 * isValid(double_id)
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::isValid(
   const DoubleAttributeId& id) const
{
   return (id() >= 0) && (id() < d_dbl_attr_size);
}

/**********************************************************************
 * isValid(int_id)
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::isValid(
   const IntegerAttributeId& id) const
{
   return (id() >= 0) && (id() < d_dbl_attr_size);
}

/**********************************************************************
 * begin()
 *********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::iterator
SparseData<BOX_GEOMETRY>::begin()
{
   return SparseDataIterator<BOX_GEOMETRY>(this);
}

/**********************************************************************
 * end()
 *********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::iterator
SparseData<BOX_GEOMETRY>::end()
{
   return SparseDataIterator<BOX_GEOMETRY>(
             *this, d_index_to_attribute_map.end());
}

/**********************************************************************
 * begin(index)
 *********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::AttributeIterator
SparseData<BOX_GEOMETRY>::begin(
   const hier::Index& index)
{
   return SparseDataAttributeIterator<BOX_GEOMETRY>(
             d_index_to_attribute_map[index],
             d_index_to_attribute_map[index].begin());
}

/**********************************************************************
 * end(index)
 *********************************************************************/
template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::AttributeIterator
SparseData<BOX_GEOMETRY>::end(
   const hier::Index& index)
{
   return SparseDataAttributeIterator<BOX_GEOMETRY>(
             d_index_to_attribute_map[index], d_index_to_attribute_map[index].end());
}

/**********************************************************************
 * printNames()
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::printNames(std::ostream& out) const
{
   typename std::unordered_map<
      std::string, DoubleAttributeId>::const_iterator
   dbl_iter = d_dbl_names.begin();

   typename std::unordered_map<
      std::string, IntegerAttributeId>::const_iterator
   int_iter = d_int_names.begin();

   for ( ; dbl_iter != d_dbl_names.end(); ++dbl_iter) {
      out << dbl_iter->first << ": "
          << dbl_iter->second() << std::endl;
   }

   for ( ; int_iter != d_int_names.end(); ++int_iter) {
      out << int_iter->first << ": "
          << int_iter->second() << std::endl;
   }
}

/**********************************************************************
 * printAttributes()
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::printAttributes(
   std::ostream& out) const
{
   typename IndexMap::const_iterator i = d_index_to_attribute_map.begin();
   typename SparseData<BOX_GEOMETRY>::AttributeList::const_iterator list_iter;
   ///typename std::list<Attributes>::const_iterator list_iter;
   for ( ; i != d_index_to_attribute_map.end(); ++i) {
      out << "Index: " << i->first << std::endl;
      out << "====" << i->second.size() << " items ====" << std::endl;
      for (list_iter = i->second.begin(); list_iter != i->second.end();
           ++list_iter) {
         list_iter->printAttributes(out);
      }
      out << std::endl;
   }
}

/**********************************************************************
* equality
**********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::operator == (
   const SparseData<BOX_GEOMETRY>& other) const
{
   return d_index_to_attribute_map == other.d_index_to_attribute_map;
}

/**********************************************************************
* inequality
**********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseData<BOX_GEOMETRY>::operator != (
   const SparseData<BOX_GEOMETRY>& other) const
{
   return d_index_to_attribute_map != other.d_index_to_attribute_map;
}

template<typename BOX_GEOMETRY>
typename SparseData<BOX_GEOMETRY>::AttributeList
& SparseData<BOX_GEOMETRY>::_get(
   const hier::Index & index) const
{
   typename SparseData<BOX_GEOMETRY>::AttributeList * list = 0;
   try
   {
      list = &d_index_to_attribute_map.at(index);
   } catch (std::out_of_range e) {
      TBOX_ASSERT_MSG(list != 0,
         "The index was not found in this sparse data object");
   }
   return *list;
}

/**********************************************************************
 * _add(IndexMap::const_iterator)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::_add(
   const typename IndexMap::const_iterator& item_to_add)
{
   registerIndex(item_to_add->first);
   d_index_to_attribute_map[item_to_add->first] = item_to_add->second;
}

/**********************************************************************
 * _removeInsideBox
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseData<BOX_GEOMETRY>::_removeInsideBox(
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   typename IndexMap::iterator index_map_iter =
      d_index_to_attribute_map.begin();

   while (index_map_iter != d_index_to_attribute_map.end()) {
      if (box.contains(index_map_iter->first)) {
         // erase returns the next iterator in the list prior to erasure
         index_map_iter = d_index_to_attribute_map.erase(index_map_iter);
      } else {
         ++index_map_iter;
      }
   }
}

template <typename BOX_GEOMETRY>
template <class T>
inline void 
SparseData<BOX_GEOMETRY>::hash_combine(std::size_t& seed, const T& v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <typename BOX_GEOMETRY>
template <class T>
inline void
SparseData<BOX_GEOMETRY>::to_lower(T& input)
{
  std::transform(input.begin(), input.end(), input.begin(), ::tolower);
}

/**********************************************************************
 * SparseDataIterator methods
 *********************************************************************/

/**********************************************************************
 * ctor's
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator():
   d_data(0)
{
}

/**********************************************************************
 * ctor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>& sparse_data):
   d_data(&sparse_data)
{
   d_iterator = d_data->d_index_to_attribute_map.begin();
}

/**********************************************************************
 * ctor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>* sparse_data):
   d_data(sparse_data)
{
   d_iterator = d_data->d_index_to_attribute_map.begin();
}

/**********************************************************************
 * private c'tor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>& sparse_data,
   typename SparseData<BOX_GEOMETRY>::IndexMap::iterator iterator):
   d_data(&sparse_data),
   d_iterator(iterator)
{
}

/**********************************************************************
 * copy ctor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   const SparseDataIterator<BOX_GEOMETRY>& other)
{
   if (this != &other) {
      d_data = other.d_data;
      d_iterator = other.d_iterator;
   }
}

/**********************************************************************
 * d'tor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::~SparseDataIterator()
{
   d_iterator = d_data->d_index_to_attribute_map.end();
   d_data = 0;
}

/**********************************************************************
 * assignment operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>&
SparseDataIterator<BOX_GEOMETRY>::operator = (
   const SparseDataIterator<BOX_GEOMETRY>& rhs)
{
   if (this != &rhs) {
      d_data = 0;
      d_data = rhs.d_data;
      d_iterator = rhs.d_iterator;
   }
   return *this;
}

/**********************************************************************
 * operator==
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseDataIterator<BOX_GEOMETRY>::operator == (
   const SparseDataIterator<BOX_GEOMETRY>& rhs) const
{
   return d_iterator == rhs.d_iterator;
}

/**********************************************************************
 * inequality
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseDataIterator<BOX_GEOMETRY>::operator != (
   const SparseDataIterator<BOX_GEOMETRY>& rhs) const
{
   return !this->operator == (rhs);
}

/**********************************************************************
 * pre-increment operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>&
SparseDataIterator<BOX_GEOMETRY>::operator ++ ()
{
   ++d_iterator;
   return *this;
}

/**********************************************************************
 * post-increment operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>
SparseDataIterator<BOX_GEOMETRY>::operator ++ (int)
{
   SparseDataIterator<BOX_GEOMETRY> tmp = *this;
   ++d_iterator;
   return tmp;
}

/**********************************************************************
* getIndex()
**********************************************************************/
template<typename BOX_GEOMETRY>
const hier::Index&
SparseDataIterator<BOX_GEOMETRY>::getIndex() const
{
   return d_iterator->first;
}

/**********************************************************************
 * insert(double*, int*)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataIterator<BOX_GEOMETRY>::insert(
   const double* dvalues, const int* ivalues)
{
   typename SparseData<BOX_GEOMETRY>::Attributes
   tmp(d_data->d_dbl_attr_size, d_data->d_int_attr_size);
   tmp.add(dvalues, ivalues);
   d_iterator->second.push_back(tmp);
}

/**********************************************************************
 * insert(std::vector<double>&, std::vector<int>&)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataIterator<BOX_GEOMETRY>::insert(
   const std::vector<double>& dvalues,
   const std::vector<int>& ivalues)
{
   typename SparseData<BOX_GEOMETRY>::Attributes
   tmp(d_data->d_dbl_attr_size, d_data->d_int_attr_size);
   tmp.add(dvalues, ivalues);
   d_iterator->second.push_back(tmp);
}

/**********************************************************************
 * equals(Iterator)
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseDataIterator<BOX_GEOMETRY>::equals(
   const SparseDataIterator<BOX_GEOMETRY>& rhs) const
{
   bool success = d_iterator->first == rhs.d_iterator->first;
   success = success && (d_iterator->second == rhs.d_iterator->second);
   return success;
}

/**********************************************************************
* move(toIndex)
**********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataIterator<BOX_GEOMETRY>::move(
   const hier::Index& toIndex)
{
   // ensure that the index exists.
   TBOX_ASSERT(d_data->d_index_to_attribute_map->count(toIndex) == 1);

   AttributeList& list = d_data->d_index_to_attribute_map[toIndex];
   AttributeList& this_list = d_iterator->second;
   list.insert(list.end(), this_list.begin(), this_list.end());

   // NOTE:  d_iterator++ increments so that it refers to the next element but
   // it will yield a copy of the original value.  Thus, d_iterator does not
   // refer to the element that is removed when erase() is called.
   d_data->d_index_to_attribute_map.erase(d_iterator++);
}

/**********************************************************************
 * insert(Attributes)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataIterator<BOX_GEOMETRY>::_insert(
   const typename SparseData<BOX_GEOMETRY>::Attributes& attributes)
{
   d_iterator->second.push_back(attributes);
}

/**********************************************************************
 * printIterator(ostream)
 *********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataIterator<BOX_GEOMETRY>::printIterator(
   std::ostream& out) const
{
   typename AttributeList::const_iterator iter = d_iterator->second.begin();

   out << "Index: " << d_iterator->first << std::endl;
   for ( ; iter != d_iterator->second.end(); ++iter) {
      iter->printAttributes(out);
   }
}

/**********************************************************************
 * SparseDataAttributeIterator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>::SparseDataAttributeIterator(
   const SparseData<BOX_GEOMETRY>& sparse_data,
   const hier::Index& index):
   d_list(sparse_data._get(index)),
   d_list_iterator(sparse_data._get(index).begin())
{
}

/**********************************************************************
 * copy c'tor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>::SparseDataAttributeIterator(
   const SparseDataAttributeIterator<BOX_GEOMETRY>& other)
{
   if (this != &other) {
      d_list = other.d_list;
      d_list_iterator = other.d_list_iterator;
   }
}

/**********************************************************************
 * private c'tor
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>::SparseDataAttributeIterator(
   const AttributeList& attributes,
   const typename AttributeList::iterator& iterator):
   d_list(attributes),
   d_list_iterator(iterator)
{
}

/**********************************************************************
 * equality operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseDataAttributeIterator<BOX_GEOMETRY>::operator == (
   const SparseDataAttributeIterator<BOX_GEOMETRY>& rhs) const
{
   return d_list_iterator == rhs.d_list_iterator;
}

/**********************************************************************
 * inequality operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
bool
SparseDataAttributeIterator<BOX_GEOMETRY>::operator != (
   const SparseDataAttributeIterator<BOX_GEOMETRY>& rhs) const
{
   return !this->operator == (rhs);
}

/**********************************************************************
 * pre-increment operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator ++ ()
{
   ++d_list_iterator;
   return *this;
}

/**********************************************************************
 * post-increment operator
 *********************************************************************/
template<typename BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>
SparseDataAttributeIterator<BOX_GEOMETRY>::operator ++ (int)
{
   SparseDataAttributeIterator<BOX_GEOMETRY> tmp = *this;
   ++d_list_iterator;
   return tmp;
}

/**********************************************************************
 * access operators
 *********************************************************************/
template<typename BOX_GEOMETRY>
double&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator [] (
   const DoubleAttributeId& id)
{
   return (*d_list_iterator)[id];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
const double&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator [] (
   const DoubleAttributeId& id) const
{
   return (*d_list_iterator)[id];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
int&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator [] (
   const IntegerAttributeId& id)
{
   return (*d_list_iterator)[id];
}

/**********************************************************************
 *
 *********************************************************************/
template<typename BOX_GEOMETRY>
const int&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator [] (
   const IntegerAttributeId& id) const
{
   return (*d_list_iterator)[id];
}

/**********************************************************************
* printAttribute(out)
**********************************************************************/
template<typename BOX_GEOMETRY>
void
SparseDataAttributeIterator<BOX_GEOMETRY>::printAttribute(
   std::ostream& out) const
{
   d_list_iterator->printAttributes(out);
}

/**********************************************************************
* output operators for SparseDataIterator and SparseDataAttributeIterator
**********************************************************************/
template<typename BOX_GEOMETRY>
std::ostream& operator << (std::ostream& out,
                           SparseDataIterator<BOX_GEOMETRY>& sparse_data_iterator)
{
   sparse_data_iterator.printIterator(out);
   return out;
}

template<typename BOX_GEOMETRY>
std::ostream& operator << (std::ostream& out,
                           SparseDataAttributeIterator<BOX_GEOMETRY>& attr_iterator)
{
   attr_iterator.printAttribute(out);
   return out;
}


} // namespace pdat
} // namespace SAMRAI
#endif
