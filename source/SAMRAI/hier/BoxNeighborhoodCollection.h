/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A class describing the adjacency of Boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxNeighborhoodCollection
#define included_hier_BoxNeighborhoodCollection

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include <map>
#include <set>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Given a Box in a base BoxLevel, the Boxes in a head BoxLevel which
 * are adjacent to the base Box are its neighbors and are said to form the
 * neighborhood of the base Box.  This class describes the neighborhoods of a
 * collection of base Boxes.  Each base Box in the collection has a
 * neighborhood of adjacent head Boxes.
 */
class BoxNeighborhoodCollection
{
   friend class Iterator;
   friend class ConstIterator;

private:
   // Strict weak ordering for pointers to Boxes.
   struct box_ptr_less {
      bool
      operator () (const Box* box0, const Box* box1) const
      {
         return box0->getBoxId() < box1->getBoxId();
      }
   };

   // Strict weak ordering for Boxes.
   struct box_less {
      bool
      operator () (const Box& box0, const Box& box1) const
      {
         return box0.getBoxId() < box1.getBoxId();
      }
   };

   // Strict weak ordering for pointers to BoxIds.
   struct box_id_ptr_less {
      bool
      operator () (const BoxId* id0, const BoxId* id1) const
      {
         return *id0 < *id1;
      }
   };

   // Strict weak ordering for BoxIds.
   struct box_id_less {
      bool
      operator () (const BoxId& id0, const BoxId& id1) const
      {
         return id0 < id1;
      }
   };

   // type aliases
   using BaseBoxPool = std::set<BoxId, box_id_less>;
   using BaseBoxPoolItr = BaseBoxPool::iterator;
   using HeadBoxPool = std::set<Box, box_less>;
   using HeadBoxLinkCt = std::map<const Box*, int, box_ptr_less>;
   using Neighborhood = std::set<const Box*, box_ptr_less>;
   using NeighborhoodItr = Neighborhood::iterator;
   using NeighborhoodConstItr = Neighborhood::const_iterator;
   using AdjList = std::map<const BoxId*, Neighborhood, box_id_ptr_less>;
   using AdjListItr = AdjList::iterator;
   using AdjListConstItr = AdjList::const_iterator;

   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_BOX_NBRHD_COLLECTION_VERSION;

   /*!
    * @brief The pool of BoxIds of base Boxes.
    */
   BaseBoxPool d_base_boxes;

   /*!
    * @brief The pool of head Boxes.
    */
   HeadBoxPool d_nbrs;

   /*!
    * @brief The links between members of the pool of base Box BoxIds and
    * members of the pool of head Boxes.  The links are directed from base
    * to head hence a given base BoxId is linked to arbitrarily many head
    * Boxes.
    */
   AdjList d_adj_list;

   /*!
    * @brief The number of incident links for each member of the pool of
    * head Boxes.
    */
   HeadBoxLinkCt d_nbr_link_ct;

public:
   // Constructors.

   /*!
    * @brief Constructs an empty object. There are not yet any base Boxes
    * whose neighborhoods are represented by this object.
    */
   BoxNeighborhoodCollection() = default;

   /*!
    * @brief Constructs a collection of empty neighborhoods for each base
    * box in base_boxes.
    *
    * @param base_boxes Base boxes whose neighborhoods will be represented
    * by this object.
    */
   BoxNeighborhoodCollection(
      const BoxContainer& base_boxes);

   /*!
    * @brief Copy constructor.
    *
    * @param other
    */
   BoxNeighborhoodCollection(
      const BoxNeighborhoodCollection& other);

   /*!
    * @brief Assignment operator.
    *
    * @param rhs
    */
   BoxNeighborhoodCollection&
   operator = (
      const BoxNeighborhoodCollection& rhs);

   // Destructor

   /*!
    * @brief Destructor
    */
   ~BoxNeighborhoodCollection() = default;

   // Operators

   /*!
    * @brief Determine if two collections are equivalent.
    *
    * @param rhs
    */
   bool
   operator == (
      const BoxNeighborhoodCollection& rhs) const;

   /*!
    * @brief Determine if two collections are not equivalent.
    *
    * @param rhs
    */
   bool
   operator != (
      const BoxNeighborhoodCollection& rhs) const;

   //@{
   /*!
    * @name Iteration
    */

   class Iterator;
   class ConstNeighborIterator;
   /*!
    * @brief An iterator over the base Boxes of the neighborhoods in a const
    * BoxNeighborhoodCollection.  The interface does not allow modification
    * of the base Boxes.
    */
   class ConstIterator
   {
      friend class BoxNeighborhoodCollection;
      friend class Connector;
      friend class ConstNeighborIterator;

public:
      // Constructors.

      /*!
       * @brief Constructs an iterator over the base Boxes in the
       * supplied collection.
       *
       * @param nbrhds
       * @param from_start
       */
      ConstIterator(
         const BoxNeighborhoodCollection& nbrhds,
         bool from_start = true);

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      ConstIterator(const ConstIterator& other) = default;

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      ConstIterator(
         const Iterator& other);

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      ConstIterator&
      operator = (
         const ConstIterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_itr = rhs.d_itr;
         d_base_boxes_itr = rhs.d_base_boxes_itr;
         return *this;
      }

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      ConstIterator&
      operator = (
         const Iterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_itr = rhs.d_itr;
         d_base_boxes_itr = rhs.d_base_boxes_itr;
         return *this;
      }

      // Destructor

      /*!
       * @brief Performs necessary deletion.
       */
      ~ConstIterator() = default;

      // Operators

      /*!
       * @brief Extracts the BoxId of the base Box of the current
       * neighborhood in the iteration.
       */
      const BoxId&
      operator * () const
      {
         return *(d_itr->first);
      }

      /*!
       * @brief Extracts a pointer to the BoxId of the base Box of the
       * current neighborhood in the iteration.
       */
      const BoxId *
      operator -> () const
      {
         return d_itr->first;
      }

      /*!
       * @brief Post-increment iterator to point to BoxId of the base Box
       * of next neighborhood in the collection.
       */
      ConstIterator
      operator ++ (
         int)
      {
         // Go to the next base Box.
         ConstIterator tmp = *this;
         if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
            ++d_base_boxes_itr;
            ++d_itr;
         }
         return tmp;
      }

      /*!
       * @brief Pre-increment iterator to point to BoxId of the base Box
       * of next neighborhood in the collection.
       */
      ConstIterator&
      operator ++ ()
      {
         // Go to the next base Box.
         if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
            ++d_base_boxes_itr;
            ++d_itr;
         }
         return *this;
      }

      /*!
       * @brief Determine if two iterators are equivalent.
       *
       * @param rhs
       */
      bool
      operator == (
         const ConstIterator& rhs) const
      {
         return d_collection == rhs.d_collection &&
                d_itr == rhs.d_itr &&
                d_base_boxes_itr == rhs.d_base_boxes_itr;
      }

      /*!
       * @brief Determine if two iterators are not equivalent.
       *
       * @param rhs
       */
      bool
      operator != (
         const ConstIterator& rhs) const
      {
         return !(*this == rhs);
      }

private:
      // Default constructor does not exist.
      ConstIterator();

      /*!
       * @brief Constructs an iterator pointing to a specific base Box in
       * nbrhds.  Should only be called by BoxNeighborhoodCollection.
       *
       * @param nbrhds
       *
       * @param itr
       */
      ConstIterator(
         const BoxNeighborhoodCollection& nbrhds,
         AdjListConstItr itr);

      const BoxNeighborhoodCollection* d_collection;

      AdjListConstItr d_itr;

      BaseBoxPoolItr d_base_boxes_itr;
   };

   class NeighborIterator;

   /*!
    * @brief An iterator over the base Boxes of the neighborhoods in a
    * BoxNeighborhoodCollection.  The interface does not allow modification
    * of the base Boxes.
    */
   class Iterator
   {
      friend class BoxNeighborhoodCollection;
      friend class Connector;
      friend class ConstIterator;
      friend class NeighborIterator;

public:
      // Constructors.

      /*!
       * @brief Constructs an iterator over the base Boxes in the
       * supplied collection.
       *
       * @param nbrhds
       * @param from_start
       */
      Iterator(
         BoxNeighborhoodCollection& nbrhds,
         bool from_start = true);

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      Iterator(const Iterator& other) = default;

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      Iterator&
      operator = (
         const Iterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_itr = rhs.d_itr;
         d_base_boxes_itr = rhs.d_base_boxes_itr;
         return *this;
      }

      // Destructor

      /*!
       * @brief Performs necessary deletion.
       */
      ~Iterator() = default;

      // Operators

      /*!
       * @brief Extracts the BoxId of the base Box of the current
       * neighborhood in the iteration.
       */
      const BoxId&
      operator * () const
      {
         return *(d_itr->first);
      }

      /*!
       * @brief Extracts a pointer to the BoxId of the base Box of the
       * current neighborhood in the iteration.
       */
      const BoxId *
      operator -> () const
      {
         return d_itr->first;
      }

      /*!
       * @brief Post-increment iterator to point to BoxId of the base Box
       * of next neighborhood in the collection.
       */
      Iterator
      operator ++ (
         int)
      {
         // Go to the next base Box.
         Iterator tmp = *this;
         if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
            ++d_base_boxes_itr;
            ++d_itr;
         }
         return tmp;
      }

      /*!
       * @brief Pre-increment iterator to point to BoxId of the base Box
       * of next neighborhood in the collection.
       */
      Iterator&
      operator ++ ()
      {
         // Go to the next base Box.
         if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
            ++d_base_boxes_itr;
            ++d_itr;
         }
         return *this;
      }

      /*!
       * @brief Determine if two iterators are equivalent.
       *
       * @param rhs
       */
      bool
      operator == (
         const Iterator& rhs) const
      {
         return d_collection == rhs.d_collection &&
                d_itr == rhs.d_itr &&
                d_base_boxes_itr == rhs.d_base_boxes_itr;
      }

      /*!
       * @brief Determine if two iterators are not equivalent.
       *
       * @param rhs
       */
      bool
      operator != (
         const Iterator& rhs) const
      {
         return !(*this == rhs);
      }

private:
      // Default constructor does not exist.
      Iterator();

      /*!
       * @brief Constructs an iterator pointing to a specific base Box in
       * nbrhds.  Should only be called by BoxNeighborhoodCollection.
       *
       * @param nbrhds
       *
       * @param itr
       */
      Iterator(
         BoxNeighborhoodCollection& nbrhds,
         AdjListItr itr);

      const BoxNeighborhoodCollection* d_collection;

      AdjListItr d_itr;

      BaseBoxPoolItr d_base_boxes_itr;
   };

   /*!
    * @brief An iterator over the neighbors in the neighborhood of a base
    * Box in a const BoxNeighborhoodCollection.  The interface does not
    * allow modification of the neighbors.
    */
   class ConstNeighborIterator
   {
      friend class BoxNeighborhoodCollection;

public:
      typedef std::forward_iterator_tag iterator_category;
      typedef Box value_type;
      typedef std::ptrdiff_t difference_type;
      typedef Box * pointer;
      typedef Box& reference;

      // Constructors

      /*!
       * @brief Constructs an iterator over the neighbors of the base Box
       * pointed to by the supplied Iterator in the supplied collection
       * of neighborhoods.
       *
       * @param base_box_itr
       * @param from_start If true constructs an iterator pointing to the
       * first neighbor.  Otherwise constructs an iterator pointing one
       * past the last neighbor.
       */
      ConstNeighborIterator(
         const ConstIterator& base_box_itr,
         bool from_start = true);

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      ConstNeighborIterator(const ConstNeighborIterator& other) = default;

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      ConstNeighborIterator(
         const NeighborIterator& other);

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      ConstNeighborIterator&
      operator = (
         const ConstNeighborIterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_base_box = rhs.d_base_box;
         d_itr = rhs.d_itr;
         return *this;
      }

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      ConstNeighborIterator&
      operator = (
         const NeighborIterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_base_box = rhs.d_base_box;
         d_itr = rhs.d_itr;
         return *this;
      }

      // Destructor

      /*!
       * @brief Performs necessary deletion.
       */
      ~ConstNeighborIterator() = default;

      // Operators

      /*!
       * @brief Extract the Box which is the current neighbor in the
       * iteration of the neighborhood of the base Box.
       */
      const Box&
      operator * () const
      {
         return *(*d_itr);
      }

      /*!
       * @brief Extracts a pointer to the Box which is current neighbor
       * in the iteration of the neighborhood of the base Box.
       */
      const Box *
      operator -> () const
      {
         return *d_itr;
      }

      /*!
       * @brief Post-increment iterator to point to the Box which is the
       * next neighbor of the base Box.
       */
      ConstNeighborIterator
      operator ++ (
         int)
      {
         ConstNeighborIterator tmp = *this;
         if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
            ++d_itr;
         }
         return tmp;
      }

      /*!
       * @brief Pre-increment iterator to point to the Box which is the
       * next neighbor of the base Box.
       */
      ConstNeighborIterator&
      operator ++ ()
      {
         if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
            ++d_itr;
         }
         return *this;
      }

      /*!
       * @brief Determine if two iterators are equivalent.
       *
       * @param rhs
       */
      bool
      operator == (
         const ConstNeighborIterator& rhs) const
      {
         return d_collection == rhs.d_collection &&
                d_base_box == rhs.d_base_box &&
                d_itr == rhs.d_itr;
      }

      /*!
       * @brief Determine if two iterators are not equivalent.
       *
       * @param rhs
       */
      bool
      operator != (
         const ConstNeighborIterator& rhs) const
      {
         return !(*this == rhs);
      }

private:
      // Default constructor does not exist.
      ConstNeighborIterator();

      const BoxNeighborhoodCollection* d_collection;

      const BoxId* d_base_box;

      NeighborhoodConstItr d_itr;
   };

   /*!
    * @brief An iterator over the neighbors in the neighborhood of a base
    * Box in a BoxNeighborhoodCollection.  The interface does not allow
    * modification of the neighbors.
    */
   class NeighborIterator
   {
      friend class BoxNeighborhoodCollection;
      friend class ConstNeighborIterator;

public:
      // Constructors

      /*!
       * @brief Constructs an iterator over the neighbors of the base Box
       * pointed to by the supplied Iterator in the supplied collection
       * of neighborhoods.
       *
       * @param base_box_itr
       * @param from_start If true constructs an iterator pointing to the
       * first neighbor.  Otherwise constructs an iterator pointing one
       * past the last neighbor.
       */
      NeighborIterator(
         Iterator& base_box_itr,
         bool from_start = true);

      /*!
       * @brief Copy constructor.
       *
       * @param other
       */
      NeighborIterator(const NeighborIterator& other) = default;

      /*!
       * @brief Assignment operator.
       *
       * @param rhs
       */
      NeighborIterator&
      operator = (
         const NeighborIterator& rhs)
      {
         d_collection = rhs.d_collection;
         d_base_box = rhs.d_base_box;
         d_itr = rhs.d_itr;
         return *this;
      }

      // Destructor

      /*!
       * @brief Performs necessary deletion.
       */
      ~NeighborIterator() = default;

      // Operators

      /*!
       * @brief Extract the Box which is the current neighbor in the
       * iteration of the neighborhood of the base Box.
       */
      const Box&
      operator * () const
      {
         return *(*d_itr);
      }

      /*!
       * @brief Extracts a pointer to the Box which is current neighbor
       * in the iteration of the neighborhood of the base Box.
       */
      const Box *
      operator -> () const
      {
         return *d_itr;
      }

      /*!
       * @brief Post-increment iterator to point to the Box which is the
       * next neighbor of the base Box.
       */
      NeighborIterator
      operator ++ (
         int)
      {
         NeighborIterator tmp = *this;
         if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
            ++d_itr;
         }
         return tmp;
      }

      /*!
       * @brief Pre-increment iterator to point to the Box which is the
       * next neighbor of the base Box.
       */
      NeighborIterator&
      operator ++ ()
      {
         if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
            ++d_itr;
         }
         return *this;
      }

      /*!
       * @brief Determine if two iterators are equivalent.
       *
       * @param rhs
       */
      bool
      operator == (
         const NeighborIterator& rhs) const
      {
         return d_collection == rhs.d_collection &&
                d_base_box == rhs.d_base_box &&
                d_itr == rhs.d_itr;
      }

      /*!
       * @brief Determine if two iterators are not equivalent.
       *
       * @param rhs
       */
      bool
      operator != (
         const NeighborIterator& rhs) const
      {
         return !(*this == rhs);
      }

private:
      // Default constructor does not exist.
      NeighborIterator();

      const BoxNeighborhoodCollection* d_collection;

      const BoxId* d_base_box;

      NeighborhoodItr d_itr;
   };

   /*!
    * @brief Returns an iterator pointing to the beginning of the collection
    * of neighborhoods.
    */
   Iterator
   begin()
   {
      return Iterator(*this);
   }

   /*!
    * @brief Returns an iterator pointing to the beginning of the collection
    * of neighborhoods.
    */
   ConstIterator
   begin() const
   {
      return ConstIterator(*this);
   }

   /*!
    * @brief Returns an iterator pointing just past the end of the
    * collection of neighborhoods.
    */
   Iterator
   end()
   {
      return Iterator(*this, false);
   }

   /*!
    * @brief Returns an iterator pointing just past the end of the
    * collection of neighborhoods.
    */
   ConstIterator
   end() const
   {
      return ConstIterator(*this, false);
   }

   /*!
    * @brief Returns an iterator pointing to the first neighbor in the
    * neighborhood of the base Box with the supplied BoxId.
    *
    * @param base_box_id
    */
   NeighborIterator
   begin(
      const BoxId& base_box_id)
   {
      Iterator itr(find(base_box_id));
      return begin(itr);
   }

   /*!
    * @brief Returns an iterator pointing to the first neighbor in the
    * neighborhood of the base Box with the supplied BoxId.
    *
    * @param base_box_id
    */
   ConstNeighborIterator
   begin(
      const BoxId& base_box_id) const
   {
      return begin(find(base_box_id));
   }

   /*!
    * @brief Returns an iterator pointing to the first neighbor in the
    * neighborhood of the base Box pointed to by base_box_itr.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   NeighborIterator
   begin(
      Iterator& base_box_itr)
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return NeighborIterator(base_box_itr);
   }

   /*!
    * @brief Returns an iterator pointing to the first neighbor in the
    * neighborhood of the base Box pointed to by base_box_itr.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   ConstNeighborIterator
   begin(
      const ConstIterator& base_box_itr) const
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return ConstNeighborIterator(base_box_itr);
   }

   /*!
    * @brief Returns an iterator pointing just past the last neighbor in the
    * neighborhood of the base Box with the supplied BoxId.
    *
    * @param base_box_id
    */
   NeighborIterator
   end(
      const BoxId& base_box_id)
   {
      Iterator itr(find(base_box_id));
      return end(itr);
   }

   /*!
    * @brief Returns an iterator pointing just past the last neighbor in the
    * neighborhood of the base Box with the supplied BoxId.
    *
    * @param base_box_id
    */
   ConstNeighborIterator
   end(
      const BoxId& base_box_id) const
   {
      return end(find(base_box_id));
   }

   /*!
    * @brief Returns an iterator pointing just past the last neighbor in the
    * neighborhood of the base Box pointed to by base_box_itr.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   NeighborIterator
   end(
      Iterator& base_box_itr)
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return NeighborIterator(base_box_itr, false);
   }

   /*!
    * @brief Returns an iterator pointing just past the last neighbor in the
    * neighborhood of the base Box pointed to by base_box_itr.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   ConstNeighborIterator
   end(
      const ConstIterator& base_box_itr) const
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return ConstNeighborIterator(base_box_itr, false);
   }

   //@}

   //@{
   /*!
    * @name Lookup
    */

   /*!
    * @brief Returns an iterator pointing to the base Box with the supplied
    * BoxId.  If no base Box's BoxId is base_box_id this method returns
    * end().
    *
    * @param base_box_id
    */
   ConstIterator
   find(
      const BoxId& base_box_id) const
   {
      BaseBoxPoolItr base_boxes_itr = d_base_boxes.find(base_box_id);
      if (base_boxes_itr == d_base_boxes.end()) {
         return end();
      } else {
         return ConstIterator(*this, d_adj_list.find(&(*base_boxes_itr)));
      }
   }

   /*!
    * @brief Returns an iterator pointing to the base Box with the supplied
    * BoxId.  If no base Box's BoxId is base_box_id this method returns
    * end().
    *
    * @param base_box_id
    */
   Iterator
   find(
      const BoxId& base_box_id)
   {
      BaseBoxPoolItr base_boxes_itr = d_base_boxes.find(base_box_id);
      if (base_boxes_itr == d_base_boxes.end()) {
         return end();
      } else {
         return Iterator(*this, d_adj_list.find(&(*base_boxes_itr)));
      }
   }

   //@}

   // Typedefs
   typedef std::pair<Iterator, bool> InsertRetType;

   //@{
   /*!
    * @name State queries
    */

   /*!
    * @brief Returns true if the number of box neighborhoods == 0.
    */
   bool
   empty() const
   {
      return d_base_boxes.empty();
   }

   /*!
    * @brief Returns the number of box neighborhoods.
    */
   int
   numBoxNeighborhoods() const
   {
      return static_cast<int>(d_base_boxes.size());
   }

   /*!
    * @brief Returns true if the neighborhood of the base Box with the
    * supplied BoxId is empty.
    *
    * @param base_box_id
    */
   bool
   emptyBoxNeighborhood(
      const BoxId& base_box_id) const
   {
      return emptyBoxNeighborhood(find(base_box_id));
   }

   /*!
    * @brief Returns true if the neighborhood of the base Box pointed to by
    * base_box_itr is empty.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   bool
   emptyBoxNeighborhood(
      const ConstIterator& base_box_itr) const
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return base_box_itr.d_itr->second.empty();
   }

   /*!
    * @brief Returns the number of neighbors in the neighborhood of the base
    * Box with the supplied BoxId.
    *
    * @param base_box_id
    */
   int
   numNeighbors(
      const BoxId& base_box_id) const
   {
      return numNeighbors(find(base_box_id));
   }

   /*!
    * @brief Returns the number of neighbors in the neighborhood of the base
    * Box pointed to by base_box_itr.
    *
    * @param base_box_itr
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   int
   numNeighbors(
      const ConstIterator& base_box_itr) const
   {
      TBOX_ASSERT(base_box_itr.d_collection == this);
      TBOX_ASSERT(base_box_itr != end());
      return static_cast<int>(base_box_itr.d_itr->second.size());
   }

   /*!
    * @brief Returns the number of neighbors in all neighborhoods.
    */
   int
   sumNumNeighbors() const;

   /*!
    * @brief Returns true if nbr is a neighbor of the base Box with the
    * supplied BoxId.
    *
    * @param base_box_id
    * @param nbr
    */
   bool
   hasNeighbor(
      const BoxId& base_box_id,
      const Box& nbr) const
   {
      return hasNeighbor(find(base_box_id), nbr);
   }

   /*!
    * @brief Returns true if nbr is a neighbor of the base Box pointed to by
    * base_box_itr.
    *
    * @param base_box_itr
    * @param nbr
    */
   bool
   hasNeighbor(
      const ConstIterator& base_box_itr,
      const Box& nbr) const;

   /*!
    * @brief Returns true if the neighborhood of the base Box with the
    * supplied BoxId is the same in this and other.
    *
    * @param base_box_id
    * @param other
    */
   bool
   neighborhoodEqual(
      const BoxId& base_box_id,
      const BoxNeighborhoodCollection& other) const;

   /*!
    * @brief Returns true if base_box_id is the BoxId of the base Box of a
    * neighborhood held by this object.
    *
    * @param base_box_id
    */
   bool
   isBaseBox(
      const BoxId& base_box_id) const
   {
      return find(base_box_id) != end();
   }

   /*!
    * @brief Returns true if all neighbors of all base Boxes are owned by
    * the processor owning this object.
    *
    * @note Currently, this is only called (via Connector::isLocal) from
    * within a small number of TBOX_ASSERTs so it is only a sanity check and
    * is not called in optimized code.  If the method becomes
    * algorithmically important an optimization would be to store a boolean
    * to indicate if there are any non-local neighbors which would avoid the
    * search that this method currently does.
    *
    * @param rank The rank of the process owning this object.
    */
   bool
   isLocal(
      int rank) const;

   //@}

   /*!
    * @brief Insert the rank of the processor owning each neighbor in each
    * neighborhood into the supplied set.
    *
    * @param owners
    */
   void
   getOwners(
      std::set<int>& owners) const;

   /*!
    * @brief Insert the rank of the processor owning each neighbor in each
    * neighborhood into the supplied set.
    *
    * @param itr
    * @param owners
    */
   void
   getOwners(
      ConstIterator& itr,
      std::set<int>& owners) const;

   //@{
   /*!
    * @name Neighborhood editing
    */

   /*!
    * @brief Inserts a new neighbor into the neighborhood of the base Box
    * with the supplied BoxId.
    *
    * @param base_box_id The BoxId of the neighborhood base Box.
    *
    * @param new_nbr The new neighbor of base_box_id.
    *
    * @return An Iterator pointing to the base Box with the supplied BoxId.
    */
   Iterator
   insert(
      const BoxId& base_box_id,
      const Box& new_nbr)
   {
      Iterator base_boxes_itr = insert(base_box_id).first;
      insert(base_boxes_itr, new_nbr);
      return base_boxes_itr;
   }

   /*!
    * @brief Inserts a new neighbor into the neighborhood of the base Box
    * pointed to by base_box_itr.
    *
    * @note base_box_itr must point to a valid base Box or this function
    * can not work.  Unlike the other versions of insert, this version has
    * no return value.  The base Box must already exist and the Iterator
    * already points to it so returning an Iterator or bool has no value.
    *
    * @param base_box_itr Iterator pointing to the base Box.
    *
    * @param new_nbr The new neighbor of the base Box.
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   void
   insert(
      Iterator& base_box_itr,
      const Box& new_nbr);

   /*!
    * @brief Inserts new neighbors into the neighborhood of the base Box
    * with the supplied BoxId.
    *
    * @param base_box_id The BoxId of the base Box.
    *
    * @param new_nbrs The new neighbors of the base Box.
    *
    * @return An Iterator pointing to the base Box with the supplied BoxId.
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   Iterator
   insert(
      const BoxId& base_box_id,
      const BoxContainer& new_nbrs)
   {
      Iterator base_boxes_itr = insert(base_box_id).first;
      insert(base_boxes_itr, new_nbrs);
      return base_boxes_itr;
   }

   /*!
    * @brief Inserts new neighbors into the neighborhood of the base Box
    * pointed to by base_box_itr.
    *
    * @note base_box_itr must point to a valid base Box or this function
    * can not work.  Unlike the other versions of insert, this version has
    * no return value.  The base Box must already exist and the Iterator
    * already points to it so returning an Iterator or bool has no value.
    *
    * @param base_box_itr Iterator pointing to the base Box.
    *
    * @param new_nbrs The new neighbors of the base Box.
    */
   void
   insert(
      Iterator& base_box_itr,
      const BoxContainer& new_nbrs);

   /*!
    * @brief Erases a neighbor from the neighborhood of the base Box with
    * the supplied BoxId.
    *
    * @param base_box_id The BoxId of the base Box.
    *
    * @param nbr The neighbor of the base Box to be erased.
    */
   void
   erase(
      const BoxId& base_box_id,
      const Box& nbr)
   {
      Iterator itr(find(base_box_id));
      erase(itr, nbr);
   }

   /*!
    * @brief Erases a neighbor from the neighborhood of the base Box pointed
    * to by base_box_itr.
    *
    * @param base_box_itr An iterator pointing to the base Box.
    *
    * @param nbr The neighbor of the base Box to be erased.
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    * @pre d_nbrs.find(nbr) != d_nbrs.end()
    */
   void
   erase(
      Iterator& base_box_itr,
      const Box& nbr);

   /*!
    * @brief Erases neighbors from the neighborhood of the base Box with the
    * supplied BoxId.
    *
    * @param base_box_id The BoxId of the base Box.
    *
    * @param nbrs The neighbors of base Box to be erased.
    */
   void
   erase(
      const BoxId& base_box_id,
      const BoxContainer& nbrs)
   {
      Iterator itr(find(base_box_id));
      erase(itr, nbrs);
   }

   /*!
    * @brief Erases neighbors from the neighborhood of the base Box pointed
    * to by base_box_itr.
    *
    * @param base_box_itr An iterator pointing to the base Box.
    *
    * @param nbrs The neighbors of base Box to be erased.
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    * @pre for each box in nbrs, d_nbrs.find(box) != d_nbrs.end()
    */
   void
   erase(
      Iterator& base_box_itr,
      const BoxContainer& nbrs);

   /*!
    * @brief Inserts a base Box with an empty neighborhood.  If the base Box
    * does not exist in this object then this function returns true.
    * Otherwise it returns false and takes no action.
    *
    * @param new_base_box The new base_box of an empty neighborhood.
    *
    * @return A pair whose first member is an Iterator and second member is
    * a bool.  If new_base_box was already a base Box then the Iterator
    * points to the already existing base Box's neighborhood and the bool is
    * false.  If new_base_box did not already exist then the Iterator points
    * to the newly inserted base Box's neighborhood and the bool is true.
    */
   InsertRetType
   insert(
      const BoxId& new_base_box);

   /*!
    * @brief Erases the neighbors of the base Box with the supplied BoxId
    * including the base Box itself.
    *
    * @param base_box_id The BoxId of the base Box whose neighbors are to be
    * erased.
    */
   void
   erase(
      const BoxId& base_box_id)
   {
      Iterator itr(find(base_box_id));
      erase(itr);
   }

   /*!
    * @brief Erases the neighbors of the base Box pointed to by base_box_itr
    * including base Box iself.
    *
    * @param base_box_itr Iterator pointing to the base box whose neighbors
    * are to be erased.
    *
    * @pre base_box_itr.d_collection == this
    * @pre base_box_itr != end()
    */
   void
   erase(
      Iterator& base_box_itr);

   /*!
    * @brief Erases the neighbors of the base Boxes pointed to in the range
    * [first_base_box_itr, last_base_box_itr) including the base Boxes
    * themselves.
    *
    * @param first_base_box_itr Iterator pointing to the first base Box
    * whose neighbors are to be erased.
    *
    * @param last_base_box_itr Iterator pointing one past the last base Box
    * whose neighbors are to be erased.
    */
   void
   erase(
      Iterator& first_base_box_itr,
      Iterator& last_base_box_itr);

   /*!
    * @brief For all base Boxes not owned by the process owning this object
    * this method erases the base Box and its neighbors.
    *
    * @param rank the rank of the process owning this object.
    */
   void
   eraseNonLocalNeighborhoods(
      int rank);

   /*!
    * @brief Erases all base Boxes having no neighbors.
    */
   void
   eraseEmptyNeighborhoods();

   /*!
    * @brief Erases all neighbors in all neighborhoods which are periodic
    * image boxes.
    */
   void
   erasePeriodicNeighbors();

   /*!
    * @brief Erases all contents so empty() == true.
    */
   void
   clear();

   //@}

   //@{
   /*!
    * @name Coarsen, refine, grow.
    */

   /*!
    * @brief Coarsens the neighbors of each base Box by the given ratio.
    *
    * @param ratio
    */
   void
   coarsenNeighbors(
      const IntVector& ratio);

   /*!
    * @brief Refines the neighbors of each base Box by the given ratio.
    *
    * @param ratio
    */
   void
   refineNeighbors(
      const IntVector& ratio);

   /*!
    * @brief Grows the neighbors of each base Box by the given amount.
    *
    * @param growth
    */
   void
   growNeighbors(
      const IntVector& growth);

   //@}

   //@{
   /*!
    * @name Neighborhood member extraction
    * Currently, the way some algorithms are implemented these are needed
    * but we may find that it is not necessary.
    */

   /*!
    * @brief Fill the supplied BoxContainer with the neighbors from all the
    * neighborhoods in this object.  The container has no notion of the
    * neighborhoods to which its contents belong.
    *
    * @param neighbors
    */
   void
   getNeighbors(
      BoxContainer& neighbors) const;

   /*!
    * @brief Fill the supplied BoxContainer with the neighbors having the
    * specified block id from all the neighborhoods in this object.  The
    * container has no notion of the neighborhoods to which its contents
    * belong.
    *
    * @param neighbors
    *
    * @param block_id
    */
   void
   getNeighbors(
      BoxContainer& neighbors,
      const BlockId& block_id) const;

   /*!
    * @brief Fill the supplied map with the neighbors from all the
    * neighborhoods in this object by block id.  The container has no notion
    * of the neighborhoods to which its contents belong.
    *
    * @param neighbors
    */
   void
   getNeighbors(
      std::map<BlockId, BoxContainer>& neighbors) const;

   /*!
    * @brief Fill the supplied BoxContainer with the neighbors of the base
    * Box with the supplied BoxId.
    *
    * @param base_box_id
    * @param neighbors
    */
   void
   getNeighbors(
      const BoxId& base_box_id,
      BoxContainer& neighbors) const;

   /*!
    * @brief Place any periodic neighbors from each neighborhood into the
    * supplied BoxContainer.
    *
    * @param result
    */
   void
   getPeriodicNeighbors(
      BoxContainer& result) const;

   //@}

   //@{
   /*!
    * @name Communication packing/unpacking
    */

   /*!
    * @brief Load an integer communication buffer with the data from this
    * object.
    *
    * @param send_mesg The integer communication buffer
    *
    * @param dim
    *
    * @param buff_init Initializer for newly allocated buffer data.
    */
   void
   putToIntBuffer(
      std::vector<int>& send_mesg,
      const tbox::Dimension& dim,
      int buff_init) const;

   /*!
    * @brief Populate object based on information contained in an integer
    * communication buffer.
    *
    * @param recv_mesg The integer communication buffer
    *
    * @param proc_offset Offset of beginning of message stream in recv_mesg
    * for each process.
    *
    * @param dim
    *
    * @param num_proc
    *
    * @param rank
    */
   void
   getFromIntBuffer(
      const std::vector<int>& recv_mesg,
      const std::vector<int>& proc_offset,
      const tbox::Dimension& dim,
      int num_proc,
      int rank);

   //@}

   //@{
   /*!
    * @name IO--these are only called from the mblktree test.
    */

   /*!
    * @brief Writes the neighborhood information to the supplied restart
    * database.
    *
    * @param restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * @brief Constructs the neighborhoods from the supplied restart
    * database.
    *
    * @param restart_db
    */
   void
   getFromRestart(
      tbox::Database& restart_db);

   //@}
};

}
}

#endif // included_hier_BoxNeighborhoodCollection
