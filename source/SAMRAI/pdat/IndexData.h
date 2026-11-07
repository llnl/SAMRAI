/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_IndexData
#define included_pdat_IndexData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Database.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace pdat {

template<class TYPE, class BOX_GEOMETRY>
class IndexDataNode;
template<class TYPE, class BOX_GEOMETRY>
class IndexIterator;

/**
 * IndexData is used for storing sparse data.  The iteration over the
 * data preserves insertion order similar to a linked list, hence data
 * may be inserted at the front (AddItem) or tail (AppendItem).
 *
 * For example, this class is used to represent embedded * boundary
 * features as a regular patch data type using the BoundaryCell class
 * as the template type.  The iteration ordering property is used to
 * visit the boundary cells in a well defined way.
 *
 * Additional insertion methods are provided for optimization in
 * specific use cases.  An extra constructor call can be avoided if
 * the "Pointer" versions of Add/AppendItem are used.   These methods
 * assume that the items are being "given" to IndexData thus IndexData
 * will delete them when needed.  The items must be dynamically
 * allocated for this to work correctly, E.G.
 *
 * index_data.addItemPointer(someIndex, new ItemToInsert());
 *
 * The "replace" versions of the insert routines may be used when one
 * wishes to preserve the existing iteration ordering.  If an item
 * already exists the replace insert will delete the old item and put
 * the item provided in it's place (inserting it into that location in
 * the linked list).
 *
 * If one does not care about the iteration order and is frequently
 * updating the items stored, the "replace" versions are significantly
 * faster.
 *
 * The template parameter TYPE * defines the storage at each index
 * location.  IndexData is derived from * hier::PatchData.
 *
 * The data type TYPE must define the following five methods which are
 * require by this class:
 *
 *    - \b - Default constructor (taking no arguments).
 *    - \b - Assignment operator; i.e.,  TYPE\& operator=(const TYPE\& rhs)
 *    - \b - Copy; copySourceItem(const hier::Index\& index,
 *                                const hier::IntVector\& src_offset,
 *                                const TYPE\& src_item)
 *    - \b - Return size of data; size_t getDataStreamSize()
 *    - \b - Pack data into message stream; i.e.,
 *             packStream(MessageStream\& stream,
 *    - \b - Unpack data from message stream; i.e.,
 *             unpackStream(MessageStream\& stream,
 *             const hier::IntVector\& offset)
 *    - \b - Write to restart;
 *             putToRestart(std::shared_ptr<tbox::Database>\& restart_db)
 *    - \b - Retrieve from restart;
 *             getFromRestart(std::shared_ptr<tbox::Database>\& restart_db)
 *
 * The BOX_GEOMETRY template parameter defines the geometry.   BOX_GEOMETRY must
 * have a nested class name Overlap that implements he following methods:
 *
 *   - \b - getSourceOffset
 *   - \b - getDestinationBoxList
 *
 *
 * More information about the templated TYPE is provided in the IndexData
 * README file.
 *
 * IndexData objects are created by the IndexDataFactory
 * factory object just as all other patch data types.
 *
 * @see IndexData
 * @see hier::PatchData
 * @see IndexDataFactory
 */

template<class TYPE, class BOX_GEOMETRY>
class IndexData:public hier::PatchData
{
public:
   /**
    * Define the iterator.
    */
   typedef IndexIterator<TYPE, BOX_GEOMETRY> iterator;

   /**
    * The constructor for an IndexData object.  The box describes the interior
    * of the index space and the ghosts vector describes the ghost nodes in
    * each coordinate direction.
    *
    * @pre box.getDim() == ghosts.getDim()
    */
   IndexData(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /**
    * The virtual destructor for an IndexData object.
    */
   virtual ~IndexData();

   /**
    * A fast copy between the source and destination.  All data is copied
    * from the source into the destination where there is overlap in the
    * index space.
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const IndexData<TYPE, BOX_GEOMETRY> *>(&src) != 0
    */
   virtual void
   copy(
      const hier::PatchData& src);

   /**
    * A fast copy between the source and destination.  All data is copied
    * from the source into the destination where there is overlap in the
    * index space.
    *
    * @pre getDim() == dst.getDim()
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description should have been computed
    * previously from computeIntersection().
    *
    * @pre getDim() == src.getDim()
    * @pre dynamic_cast<const IndexData<TYPE, BOX_GEOMETRY> *>(&src) != 0
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description should have been computed
    * previously from computeIntersection().
    *
    * @pre getDim() == dst.getDim()
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /**
    * Determines whether the hier::PatchData subclass can estinate the necessary
    * stream size using only index space information.
    */
   virtual bool
   canEstimateStreamSizeFromBox() const;

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   virtual size_t
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /**
    * Pack data lying on the specified index set into the output stream.
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /**
    * Unpack data from the message stream into the specified index set.
    *
    * @pre dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap) != 0
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   /**
    * Add a new item to the tail of the irregular index set.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   appendItem(
      const hier::Index& index,
      const TYPE& item);

   /**
    * Add a pointer to a new item to the tail of the irregular index
    * set.  IndexData will delete the item when it is no longer needed
    * by IndexData.  Due to this behavior item must be dynamically
    * created (e.g. new) so that it may be deleted.
    *
    * NOTE: This is an optimization to avoid an extra constructor
    * call.  It should be used with caution, the caller MUST NOT
    * delete the referenced item.  Think of this as giving up control
    * of the item to IndexData.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   appendItemPointer(
      const hier::Index& index,
      TYPE * item);

   /**
    * Add a new item to the head of the irregular index set
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   addItem(
      const hier::Index& index,
      const TYPE& item);

   /**
    * Add a pointer to a new item to the head of the irregular index
    * set.  IndexData will delete the item when it is no longer needed
    * by IndexData.  Due to this behavior item must be dynamically
    * created (e.g. new) so that it may be deleted.
    *
    * NOTE: This is an optimization to avoid an extra constructor
    * call.  It should be used with caution, the caller MUST NOT
    * delete the referenced item.  Think of this as giving up control
    * of the item to IndexData.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */

   void
   addItemPointer(
      const hier::Index& index,
      TYPE * item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the orginal insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to addItem.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   replaceAddItem(
      const hier::Index& index,
      const TYPE& item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the original insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to addItemPointer.
    *
    * See addItemPointer for additional comments on pointer semantics.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   replaceAddItemPointer(
      const hier::Index& index,
      TYPE * item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the orginal insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to appendItem.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   replaceAppendItem(
      const hier::Index& index,
      const TYPE& item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index and replace it with a new item.
    *
    * This preserves the iteration order of the original insertions.
    *
    * If an item does not already exist at index this is equivelent
    * to appendItemPointer.
    *
    * See addItemPointer for additional comments on pointer semantics.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    */
   void
   replaceAppendItemPointer(
      const hier::Index& index,
      TYPE * item);

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the specified hier::Index.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    * @pre (hier::PatchData::getGhostBox().offset(index) >= 0) &&
    *      (hier::PatchData::getGhostBox().offset(index) <= hier::PatchData::getGhostBox().size())
    *
    */
   void
   removeItem(
      const hier::Index& index);

   /**
    * Return the number of data items (i.e. the number of indices) in
    * the index data list.
    *
    * @pre offset >= 0 && offset <= hier::PatchData::getGhostBox().size()
    * @pre d_data[offset] != 0
    */
   size_t
   getNumberOfItems() const;

   /**
    * Remove (deallocate) any items in the irregular index set located in
    * the index space of the hier::Box.
    *
    * @pre getDim() == box.getDim()
    */
   void
   removeInsideBox(
      const hier::Box& box);

   /**
    * Remove (deallocate) any items in the irregular index set located
    * outside of the index space of the hier::Box.
    *
    * @pre getDim() == box.getDim()
    */
   void
   removeOutsideBox(
      const hier::Box& box);

   /**
    * Remove (deallocate) the items in the irregular index set located in
    * the ghost region of the patch.
    */
   void
   removeGhostItems();

   /**
    * Remove (deallocate) all items in the irregular index set.
    */
   void
   removeAllItems();

   /**
    * Returns true if there is an element of the irregular index set at
    * the specified hier::Index.
    *
    * @pre getDim() == index.getDim()
    * @pre hier::PatchData::getGhostBox().contains(index)
    */
   bool
   isElement(
      const hier::Index& index) const;

   /**
    * Given an index, return a pointer to the item located at that index.
    * If there is no item at the index, null is returned.
    *
    * @pre getDim() == index.getDim()
    */
   TYPE *
   getItem(
      const hier::Index& index) const;

   /**
    * Check to make sure that the class version number is the same
    * as the restart file version number.
    *
    * @pre restart_db
    */
   virtual void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /**
    * Write out the class version number to the restart database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

private:
   friend class IndexIterator<TYPE, BOX_GEOMETRY>;

   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_INDEXDATA_VERSION;

   /**
    * Returns true if element exists at offset
    */
   bool
   isElement(
      size_t offset) const;

   /**
    * Remove (deallocate) the item in the irregular index set located at
    * the offset.
    *
    * NOTE: This is for optimization to avoid computing
    * offset repeatedly.
    */
   void
   removeItem(
      const size_t offset);

   /**
    * Internal routine to append item to the linked list
    * representation.
    *
    * NOTE: Offset is not strictly necessary but was include to avoid
    * computing it repeatedly.
    *
    * @pre getDim() == index.getDim()
    */
   void
   addItemToList(
      const hier::Index& index,
      const size_t offset,
      TYPE& item);

   /**
    * Internal routine to add item to the linked list
    * representation.
    *
    * NOTE: Offset is not strictly necessary but was include to avoid
    * computing it repeatedly.
    *
    * @pre getDim() == index.getDim()
    */
   void
   appendItemToList(
      const hier::Index& index,
      const size_t offset,
      TYPE& item);

   /**
    * Remove the specified node from the linked list.
    */
   void
   removeNodeFromList(
      IndexDataNode<TYPE, BOX_GEOMETRY> * node);

   // Unimplemented copy constructor
   IndexData(
      const IndexData&);

   // Unimplemented assignment operator
   IndexData&
   operator = (
      const IndexData&);

   const tbox::Dimension d_dim;

   std::vector<IndexDataNode<TYPE, BOX_GEOMETRY> *> d_data;

   /*
    * Doublely linked list of nodes.
    */
   IndexDataNode<TYPE, BOX_GEOMETRY>* d_list_head;
   IndexDataNode<TYPE, BOX_GEOMETRY>* d_list_tail;
   int d_number_items;
};

/**
 * Class IndexDataNode holds items for the linked list.  This class should
 * be defined inside the IndexData class, but nested template classes can
 * cause some compilers to barf.  List nodes should never be seen by the
 * user of a list.
 *
 */
template<class TYPE, class BOX_GEOMETRY>
class IndexDataNode
{
public:
   friend class IndexData<TYPE, BOX_GEOMETRY>;
   friend class IndexIterator<TYPE, BOX_GEOMETRY>;

   IndexDataNode(
      const hier::Index & index,
      const size_t d_offset,
      TYPE & t,
      IndexDataNode<TYPE, BOX_GEOMETRY>* n,
      IndexDataNode<TYPE, BOX_GEOMETRY>* p);

   virtual ~IndexDataNode();

private:
   // Unimplemented default constructor.
   IndexDataNode();

   hier::Index d_index;
   size_t d_offset;
   TYPE* d_item;

   IndexDataNode<TYPE, BOX_GEOMETRY>* d_next;
   IndexDataNode<TYPE, BOX_GEOMETRY>* d_prev;
};

/**
 * Class IndexIterator is the iterator associated with the IndexData
 * This class provides methods for stepping through the
 * list that contains the irregular index set.  The user should
 * access this class through the name IndexData<TYPE, BOX_GEOMETRY>::iterator.
 *
 * This iterator should be used as follows:
 * \verbatim
 * IndexData<TYPE, BOX_GEOMETRY> data;
 * ...
 * IndexData<TYPE, BOX_GEOMETRY>::iterator iterend(data, false);
 * for (IndexData<TYPE, BOX_GEOMETRY>::iterator iter(data, true); iter != iterend; ++iter) {
 *    ... = *iter;
 * }
 * \endverbatim
 *
 * @see IndexData
 * @see IndexIterator
 */

template<class TYPE, class BOX_GEOMETRY>
class IndexIterator
{
public:
   /**
    * Constructor for the index list iterator.  The iterator will iterate
    * over the irregular index set of the argument data object.
    */
   IndexIterator(
      const IndexData<TYPE, BOX_GEOMETRY>& data,
      bool begin);

   /**
    * Copy constructor for the index list iterator.
    */
   IndexIterator(
      const IndexIterator<TYPE, BOX_GEOMETRY>& iterator);

   /**
    * Assignment operator for the index list iterator.
    */
   IndexIterator&
   operator = (
      const IndexIterator<TYPE, BOX_GEOMETRY>& iterator);

   /**
    * Destructor for the index list iterator.
    */
   ~IndexIterator();

   /**
    * Return a reference to the current item in the irregular index set.
    */
   TYPE&
   operator * ();

   /**
    * Return a const reference to the current item in the irregular
    * index set.
    */
   const TYPE&
   operator * () const;

   /**
    * Return a pointer to the current item in the irregular index set.
    */
   TYPE *
   operator -> ();

   /**
    * Return a const pointer to the current item in the irregular
    * index set.
    */
   const TYPE *
   operator -> () const;

   /**
    * Return the index of the current item in the irregular index set
    */
   const hier::Index&
   getIndex() const;

   /**
    * Pre-increment the iterator to point to the next item in the index set.
    */
   IndexIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next item in the index set.
    */
   IndexIterator
   operator ++ (
      int);

   /**
    * Pre-decrement the iterator to point to the previous item in the index
    * set.
    */
   IndexIterator&
   operator -- ();

   /**
    * Post-decrement the iterator to point to the previous item in the index
    * set.
    */
   IndexIterator
   operator -- (
      int);

   /**
    * Test two iterators for equality (pointing to the same item).
    */
   bool
   operator == (
      const IndexIterator<TYPE, BOX_GEOMETRY>& iterator) const;

   /**
    * Test two iterators for inequality (pointing to different items).
    */
   bool
   operator != (
      const IndexIterator<TYPE, BOX_GEOMETRY>& iterator) const;

private:
   friend class IndexData<TYPE, BOX_GEOMETRY>;

   IndexIterator(
      IndexData<TYPE, BOX_GEOMETRY>* index_data,
      IndexDataNode<TYPE, BOX_GEOMETRY>* node);

   IndexDataNode<TYPE, BOX_GEOMETRY>&
   getNode();

   IndexData<TYPE, BOX_GEOMETRY>* d_index_data;

   IndexDataNode<TYPE, BOX_GEOMETRY>* d_node;
};

}
}

#include "SAMRAI/pdat/IndexData.cpp"

#endif
