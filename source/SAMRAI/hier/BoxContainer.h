/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_BoxContainer
#define included_hier_BoxContainer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MultiblockBoxTree.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <list>
#include <set>
#include <vector>

namespace SAMRAI {
namespace hier {

class BoxContainerSingleBlockIterator;
class BoxContainerSingleOwnerIterator;
class RealBoxConstIterator;

/*!
 * @brief A container for Boxes.
 *
 * The BoxContainer stores a collection of Boxes and provides methods for
 * access and manipulation of that collection.
 *
 * A BoxContainer exists in either an "ordered" or "unordered" state.
 * The ordered state means that the Boxes have been ordered according to
 * the comparison operators defined in the BoxId class, while the unordered
 * state has no meaningful ordering of the Boxes besides the sequence that the
 * Boxes were added to the container.  Additionally the ordered state
 * requires that all Boxes in the container have a valid and unique BoxId,
 * while there is no such restriction for unordered containers.
 *
 * An ordered container can always have its state switched to unordered by
 * a call to the unorder() method.  An unordered container can also have its
 * state switched to ordered by a call to the order() method, but only under
 * certain conditions specified below in the comments for order().
 *
 * Certain methods in this class can only be called on ordered containers
 * while others can only be called on unordered containers.  Violating these
 * restrictions will result in a run-time error.
 *
 * Regardless of unordered/unordered state, all Boxes within a BoxContainer
 * must be of the same Dimension.  If a new Box added to a container has
 * a different Dimension than the Boxes already in the container, an assertion
 * failure will occur.
 *
 * An option exists to create an internal search tree representation based on
 * the spatial coordinates of the Boxes in the container.  This option can be
 * used to reduce the cost of searching operations in the methods
 * removeIntersections(), intersectBoxes(), findOverlapBoxes() and
 * hasOverlap().  This option is invoked by calling the method makeTree().
 * This option should only be used in cases where the listed search methods
 * will be called multiple times on the same unchanging BoxContainer.  The
 * cost of building the tree representaion is O(N(log(N)), while the tree
 * reduces the cost of the search operations to O(log(N)) rather than O(N),
 * thus it is adviseable to only use the tree representation when the
 * reduction in search cost is expected to outweigh the increased cost of
 * building the tree.
 *
 * Constructing the tree represenation via makeTree() will change nothing
 * about the Boxes stored in the container, nor will it change the
 * ordered/unordered state of the container.
 *
 * @see BoxId
 */
class BoxContainer
{
   friend class BoxContainerIterator;
   friend class BoxContainerConstIterator;

public:
   typedef const Box value_type;
   class BoxContainerIterator;

   /*!
    * @brief A immutable iterator over the boxes in a BoxContainer.
    *
    * If iterating over an ordered BoxContainer, then iteration will follow
    * the sequence of the BoxId-based ordering of the container.  If iterating
    * over an unordered BoxContainer, the sequence of the iteration will be
    * based on how the members of the container were added.
    *
    * @see BoxContainer
    * @see BoxContainerIterator
    */
   class BoxContainerConstIterator
   {
      friend class BoxContainer;
      friend class BoxContainerIterator;

public:
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef const Box value_type;
      typedef std::ptrdiff_t difference_type;
      typedef const Box * pointer;
      typedef const Box& reference;

      /*!
       * @brief Copy constructor.
       *
       * @param[in] other
       */
      BoxContainerConstIterator(
         const BoxContainerConstIterator& other);

      /*!
       * @brief Copy constructor to copy mutable iterator to an immutable
       * iterator.
       */
      BoxContainerConstIterator(
         const BoxContainerIterator& other);

      /*!
       * @brief Assignment operator.
       *
       * @param[in] rhs
       */
      BoxContainerConstIterator&
      operator = (
         const BoxContainerConstIterator& rhs)
      {
         if (this != &rhs) {
            d_ordered = rhs.d_ordered;
            if (d_ordered) {
               d_set_iter = rhs.d_set_iter;
            } else {
               d_list_iter = rhs.d_list_iter;
            }
         }
         return *this;
      }

      /*!
       * @brief The destructor releases all storage.
       */
      ~BoxContainerConstIterator() = default;

      /*!
       * @brief Get box corresponding to iterator's position in container.
       *
       * @return An immutable reference to the current Box in the iteration.
       */
      const Box&
      operator * () const
      {
         return d_ordered ? **d_set_iter : *d_list_iter;
      }

      /*!
       * @brief Get pointer to box at iterator's position in container.
       *
       * @return Const pointer to the current box.
       */
      const Box *
      operator -> () const
      {
         return d_ordered ? *d_set_iter : &(*d_list_iter);
      }

      /*!
       * @brief Post-increment iterator to point to next box in the container.
       *
       * @return Iterator at the position in the container before the
       * increment.
       */
      BoxContainerConstIterator
      operator ++ (
         int)
      {
         BoxContainerConstIterator return_iter(*this);
         if (d_ordered) {
            ++d_set_iter;
         } else {
            ++d_list_iter;
         }
         return return_iter;
      }

      /*!
       * @brief Pre-increment iterator to point to next box in the container.
       *
       * @return Reference to iterator at the position in the container after
       * the increment.
       */
      const BoxContainerConstIterator&
      operator ++ ()
      {
         if (d_ordered) {
            ++d_set_iter;
         } else {
            ++d_list_iter;
         }
         return *this;
      }

      /*!
       * @brief Post-decrement iterator to point to next box in the container.
       *
       * @return Iterator at the position in the container before the
       * decrement.
       */
      BoxContainerConstIterator
      operator -- (
         int)
      {
         BoxContainerConstIterator return_iter(*this);
         if (d_ordered) {
            --d_set_iter;
         } else {
            --d_list_iter;
         }
         return return_iter;
      }

      /*!
       * @brief Pre-decrement iterator to point to next box in the container.
       *
       * @return Reference to iterator at the position in the container after
       * the decrement.
       */
      const BoxContainerConstIterator&
      operator -- ()
      {
         if (d_ordered) {
            --d_set_iter;
         } else {
            --d_list_iter;
         }
         return *this;
      }

      /*!
       * @brief Equality operator.
       *
       * @return true if both iterators point to the same box.
       *
       * @param[in] other
       */
      bool
      operator == (
         const BoxContainerConstIterator& other) const
      {
         return d_ordered ? d_set_iter == other.d_set_iter :
                d_list_iter == other.d_list_iter;
      }

      /*!
       * @brief Inequality operator.
       *
       * @return true if both iterators point to different boxes.
       *
       * @param[in] other
       */
      bool
      operator != (
         const BoxContainerConstIterator& other) const
      {
         return d_ordered ? d_set_iter != other.d_set_iter :
                d_list_iter != other.d_list_iter;
      }

private:
      /*!
       * @brief Default constructor is defined but accessible only by friends.
       */
      BoxContainerConstIterator();

      /*!
       * @brief Constructor for the BoxContainerConstIterator.
       *
       * The iterator will point to the beginning or the end of the argument
       * container, depending on the from_start argument
       *
       * @param[in] container The container whose members are iterated.
       * @param[in] from_start true if iteration starts at beginning of
       * container.
       */
      explicit BoxContainerConstIterator(
         const BoxContainer& container,
         bool from_start = true);

      /*
       * Underlying iterator to be used when unordered.
       */
      std::list<Box>::const_iterator d_list_iter;

      /*
       * Underlying iterator to be used when ordered.
       */
      std::set<Box *, Box::id_less>::const_iterator d_set_iter;

      bool d_ordered;
   };

   /*!
    * @brief A mutable iterator over the boxes in a BoxContainer.
    *
    * If iterating over an ordered BoxContainer, then iteration will follow
    * the sequence of the BoxId-based ordering of the container.  If iterating
    * over an unordered BoxContainer, the sequence of the iteration will be
    * based on how the members of the container were added.
    *
    * @see BoxContainer
    * @see BoxContainerConstIterator
    */
   class BoxContainerIterator
   {
      friend class BoxContainer;
      friend class BoxContainerConstIterator;

public:
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef Box value_type;
      typedef std::ptrdiff_t difference_type;
      typedef Box * pointer;
      typedef Box& reference;

      /*!
       * @brief Copy constructor.
       *
       * @param[in] other
       */
      BoxContainerIterator(
         const BoxContainerIterator& other);

      /*!
       * @brief Assignment operator.
       *
       * @param[in] rhs
       */
      BoxContainerIterator&
      operator = (
         const BoxContainerIterator& rhs)
      {
         if (this != &rhs) {
            d_ordered = rhs.d_ordered;
            if (d_ordered) {
               d_set_iter = rhs.d_set_iter;
            } else {
               d_list_iter = rhs.d_list_iter;
            }
         }
         return *this;
      }

      /*!
       * @brief The destructor releases all storage.
       */
      ~BoxContainerIterator() = default;

      /*!
       * @brief Get box corresponding to iterator's position in container.
       *
       * @return A mutable reference to the current Box in the iteration.
       */
      Box&
      operator * () const
      {
         return d_ordered ? **d_set_iter : *d_list_iter;
      }

      /*!
       * @brief Get pointer to box at iterator's position in container.
       *
       * @return Pointer to the current box.
       */
      Box *
      operator -> () const
      {
         return d_ordered ? *d_set_iter : &(*d_list_iter);
      }

      /*!
       * @brief Post-increment iterator to point to next box in the container.
       *
       * @return Iterator at the position in the container before the
       * increment.
       */
      BoxContainerIterator
      operator ++ (
         int)
      {
         BoxContainerIterator return_iter(*this);
         if (d_ordered) {
            ++d_set_iter;
         } else {
            ++d_list_iter;
         }
         return return_iter;
      }

      /*!
       * @brief Pre-increment iterator to point to next box in the container.
       *
       * @return Reference to iterator at the position in the container after
       * the increment.
       */
      const BoxContainerIterator&
      operator ++ ()
      {
         if (d_ordered) {
            ++d_set_iter;
         } else {
            ++d_list_iter;
         }
         return *this;
      }

      /*!
       * @brief Post-decrement iterator to point to next box in the container.
       *
       * @return Iterator at the position in the container before the
       * decrement.
       */
      BoxContainerIterator
      operator -- (
         int)
      {
         BoxContainerIterator return_iter(*this);
         if (d_ordered) {
            --d_set_iter;
         } else {
            --d_list_iter;
         }
         return return_iter;
      }

      /*!
       * @brief Pre-decrement iterator to point to next box in the container.
       *
       * @return Reference to iterator at the position in the container after
       * the decrement.
       */
      const BoxContainerIterator&
      operator -- ()
      {
         if (d_ordered) {
            --d_set_iter;
         } else {
            --d_list_iter;
         }
         return *this;
      }

      /*!
       * @brief Equality operators
       *
       * @return true if both iterators point to the same box.
       *
       * @param[in] other
       */
      bool
      operator == (
         const BoxContainerIterator& other) const
      {
         return d_ordered ? d_set_iter == other.d_set_iter :
                d_list_iter == other.d_list_iter;
      }

      bool
      operator == (
         const BoxContainerConstIterator& other) const
      {
         return d_ordered ? d_set_iter == other.d_set_iter :
                d_list_iter == other.d_list_iter;
      }

      /*!
       * @brief Inequality operators.
       *
       * @return true if both iterators point to different boxes.
       *
       * @param[in] other
       */
      bool
      operator != (
         const BoxContainerIterator& other) const
      {
         return d_ordered ? d_set_iter != other.d_set_iter :
                d_list_iter != other.d_list_iter;
      }

      bool
      operator != (
         const BoxContainerConstIterator& other) const
      {
         return d_ordered ? d_set_iter != other.d_set_iter :
                d_list_iter != other.d_list_iter;
      }

private:
      /*!
       * @brief Default constructor is defined but accessible only by friends.
       */
      BoxContainerIterator();

      /*!
       * @brief Constructor for the BoxContainerIterator.
       *
       * The iterator will point to the beginning or the end of the argument
       * container, depending on the from_start argument
       *
       * @param[in] container The container whose members are iterated.
       * @param[in] from_start true if iteration starts at beginning of
       * container.
       */
      explicit BoxContainerIterator(
         BoxContainer& container,
         bool from_start = true);

      /*
       * Underlying iterator to be used when unordered.
       */
      std::list<Box>::iterator d_list_iter;

      /*
       * Underlying iterator to be used when ordered.
       */
      std::set<Box *, Box::id_less>::iterator d_set_iter;

      bool d_ordered;

   };

   /*!
    * @brief The iterator for class BoxContainer.
    */
   typedef BoxContainerIterator iterator;

   /*!
    * @brief The const iterator for class BoxContainer.
    */
   typedef BoxContainerConstIterator const_iterator;

   //@{ @name Constructors, Destructors, Assignment

   /*!
    * @brief Default constructor creates empty container in unordered state.
    */
   BoxContainer();

   /*!
    * @brief Creates empty container in state determined by boolean
    *
    * @param[in] ordered  Container will be ordered if true, unordered if false.
    */
   explicit BoxContainer(
      const bool ordered);

   /*!
    * @brief Move constructor.
    *
    * Constructs a BoxContainer by taking ownership of the internal
    * storage and state from another BoxContainer.
    *
    * @param[in] other  BoxContainer to move from.
    */
   BoxContainer(
      BoxContainer&& other) noexcept
   : d_list(std::move(other.d_list)),
     d_set(std::move(other.d_set)),
     d_ordered(other.d_ordered),
     d_tree(std::move(other.d_tree))
   {
      // No additional work needed.
   }

   /*!
    * @brief Move assignment operator.
    *
    * Replaces the contents of this BoxContainer by taking ownership of
    * the internal storage and state from another BoxContainer.
    *
    * @param[in] other  BoxContainer to move from.
    *
    * @return Reference to this BoxContainer.
    */
   BoxContainer&
   operator =(
      BoxContainer&& other) noexcept
   {
      if (this != &other) {
         d_list = std::move(other.d_list);
         d_set = std::move(other.d_set);
         d_ordered = other.d_ordered;
         d_tree = std::move(other.d_tree);
      }
      return *this;
   }
   /*!
    * @brief Create container containing members from another container.
    *
    * Members in the range [first, last) are copied to new container.
    *
    * @param[in] first
    * @param[in] last
    * @param[in] ordered  Container will be ordered if true, unordered if false.
    */
   BoxContainer(
      const_iterator first,
      const_iterator last,
      const bool ordered = false);

   /*!
    * @brief Create a container with 1 box.
    *
    * @param[in] box  Box to copy into new container.
    * @param[in] ordered  Container will be ordered if true, unordered if false.
    */
   explicit BoxContainer(
      const Box& box,
      const bool ordered = false);

   /*!
    * @brief Copy constructor from another BoxContainer.
    *
    * All boxes and the ordered/unordered state will be copied to the new
    * BoxContainer.
    *
    * @param[in] other
    */
   BoxContainer(
      const BoxContainer& other);

   /*!
    * @brief Copy constructor from an array of tbox::DatabaseBox objects.
    *
    * The new BoxContainer will be unordered.
    *
    * @param[in] other
    */
   explicit BoxContainer(
      const std::vector<tbox::DatabaseBox>& other);

   /*!
    * @brief Constructor that copies only Boxes having the given BlockId
    * from the other container.
    *
    * The unordered or ordered state will be the same as that of the argument
    * container.
    *
    * @param[in] other
    * @param[in] block_id
    */
   BoxContainer(
      const BoxContainer& other,
      const BlockId& block_id);

   /*!
    * @brief Assignment from other BoxContainer.
    *
    * All boxes and the ordered/unordered state will be copied to the
    * assigned BoxContainer.  Any previous state of the assigned
    * BoxContainer is discarded.
    *
    * @param[in] rhs
    */
   BoxContainer&
   operator = (
      const BoxContainer& rhs);

   /*!
    * @brief Assignment from an array of tbox::DatabaseBox objects.
    *
    * The assigned BoxContainer will be unordered.  Any previous state of the
    * assigned BoxContainer is discarded.
    *
    * @param[in] rhs
    */
   BoxContainer&
   operator = (
      const std::vector<tbox::DatabaseBox>& rhs);

   /*!
    * @brief The destructor releases all storage.
    */
   ~BoxContainer() = default;

   //@}

   //@{ @name Methods that may be called on ordered or unordered BoxContainers.

   /*!
    * @brief Return the number of boxes in the container.
    *
    * @return The number of boxes in the container.
    */
   int
   size() const
   {
      if (!d_ordered) {
         return static_cast<int>(d_list.size());
      } else {
         return static_cast<int>(d_set.size());
      }
   }

   /*!
    * @brief Returns true if there are no boxes in the container.
    *
    * This version follows the naming standards used in STL.
    *
    * @return True if the container is empty.
    */
   bool
   empty() const
   {
      return d_list.empty();
   }

   /*!
    * @brief Return a const_iterator pointing to the start of the container.
    *
    * @return An immutable iterator pointing to the first box.
    */
   const_iterator
   begin() const
   {
      return const_iterator(*this);
   }

   /*!
    * @brief Return a const_iterator pointing to the end of the container.
    *
    * @return An immutable iterator pointing beyond the last box.
    */
   const_iterator
   end() const
   {
      return const_iterator(*this, false);
   }

   /*!
    * @brief Return an iterator pointing to the start of the container.
    *
    * @return A mutable iterator pointing to the first box.
    */
   iterator
   begin()
   {
      return iterator(*this);
   }

   /*!
    * @brief Return an iterator pointing to the end of the container.
    *
    * @return A mutable iterator pointing beyond the last box.
    */
   iterator
   end()
   {
      return iterator(*this, false);
   }

   /*!
    * @brief Return a BoxContainerSingleBlockIterator pointing to the first
    * box in the container with the given BlockId.
    *
    * @param block_id The BlockId of the boxes we want.
    *
    * @return A mutable iterator pointing to the first box with the given
    * BlockId.
    */
   BoxContainerSingleBlockIterator
   begin(
      const BlockId& block_id) const;

   /*!
    * @brief Return a BoxContainerSingleBlockIterator pointing to the end of
    * the container.
    *
    * @param block_id The BlockId of the boxes we want.
    *
    * @return A mutable iterator pointing beyond the last box.
    */
   BoxContainerSingleBlockIterator
   end(
      const BlockId& block_id) const;

   /*!
    * @brief Return a BoxContainerSingleOwnerIterator pointing to the first
    * box in the container having the given owner.
    *
    * @param owner_rank The processaor whose boxes we want.
    *
    * @return A mutable iterator pointing to the first box having the given
    * owner.
    */
   BoxContainerSingleOwnerIterator
   begin(
      const int& owner_rank) const;

   /*!
    * @brief Return a BoxContainerSingleOwnerIterator pointing to the end of
    * the container.
    *
    * @param owner_rank The processaor whose boxes we want.
    *
    * @return A mutable iterator pointing beyond the last box.
    */
   BoxContainerSingleOwnerIterator
   end(
      const int& owner_rank) const;

   /*!
    * @brief Return a RealBoxConstIterator pointing to the first real
    * (non-periodic) box in the container.
    *
    * @return A mutable iterator pointing to the first non-periodic box.
    */
   RealBoxConstIterator
   realBegin() const;

   /*!
    * @brief Return a RealBoxConstIterator pointing to the end of the
    * container.
    *
    * @return A mutable iterator pointing beyond the last box.
    */
   RealBoxConstIterator
   realEnd() const;

   /*!
    * @brief Returns the first element in the container.
    *
    * @return A const reference to the first Box in the container.
    */
   const Box&
   front() const
   {
      return d_ordered ? **(d_set.begin()) : d_list.front();
   }

   /*!
    * @brief Returns the first element in the container.
    *
    * @return A const reference to the last Box in the container.
    */
   const Box&
   back() const
   {
      return d_ordered ? **(d_set.rbegin()) : d_list.back();
   }

   /*!
    * @brief Remove the member of the container pointed to by "iter".
    *
    * Can be called on ordered or unordered containers.
    *
    * @param[in] iter
    */
   void
   erase(
      iterator iter);

   /*!
    * @brief Remove the members of the container in the range [first, last).
    *
    * Can be called on ordered or unordered containers.
    *
    * @param[in] first
    * @param[in] last
    */
   void
   erase(
      iterator first,
      iterator last);

   /*!
    * @brief Removes all the members of the container.
    *
    * Can be called on unordered or unordered containers.  Sets the state to
    * unordered.
    */
   void
   clear()
   {
      d_list.clear();
      d_set.clear();
      d_ordered = false;
      d_tree.reset();
   }

   /*!
    * @brief  Swap all contents and state with another BoxContainer.
    *
    * This container and other container exchange all member Boxes and
    * ordered/unordered state.
    *
    * @param[in,out] other  Other container for swap.
    */
   void
   swap(
      BoxContainer& other)
   {
      d_list.swap(other.d_list);
      d_set.swap(other.d_set);
      bool other_set_created = other.d_ordered;
      other.d_ordered = d_ordered;
      d_ordered = other_set_created;
      d_tree.swap(other.d_tree);
   }

   /*!
    * @brief  Get all of the ranks that own Boxes in this container
    *
    * The rank of every member of this container is inserted into the set.
    *
    * @param[out] owners
    *
    * @pre isOrdered() || for each box in container has valid BoxId
    */
   void
   getOwners(
      std::set<int>& owners) const;

   /*!
    * @brief Grow boxes in the container by the specified ghost cell width.
    *
    * @param[in] ghosts
    */
   void
   grow(
      const IntVector& ghosts);

   /*!
    * @brief Shift boxes in the container by the specified offset.
    *
    * @param[in] offset
    */
   void
   shift(
      const IntVector& offset);

   /*!
    * @brief Refine boxes in container by the specified refinement ratio.
    *
    * @param[in] ratio
    */
   void
   refine(
      const IntVector& ratio);

   /*!
    * @brief Coarsen boxes in container by the specified coarsening ratio.
    *
    * @param[in] ratio
    */
   void
   coarsen(
      const IntVector& ratio);

   /*!
    * @brief Count total number of indices in the boxes in the container.
    *
    * @return Total number of indices of all boxes in the container.
    */
   size_t
   getTotalSizeOfBoxes() const;

   /*!
    * @brief Determine if "idx" lies within bounds of boxes in container.
    *
    * Only boxes of the given BlockId will be checked.
    *
    * @return true if idx lies within bounds of boxes in container.
    *
    * @param[in] idx
    * @param[in] block_id
    *
    * @pre each box in container must have a valid BlockId
    */
   bool
   contains(
      const Index& idx,
      const BlockId& block_id) const;

   /*!
    * @brief  Returns the bounding box for all the boxes in the container.
    *
    * @pre !empty()
    * @pre each Box in container has same BlockId
    */
   Box
   getBoundingBox() const;

   /*!
    * @brief  Returns the bounding box for all the boxes in the container
    *         having the given BlockId.
    *
    * @param[in] block_id
    *
    * @pre !empty()
    */
   Box
   getBoundingBox(
      const BlockId& block_id) const;

   /*!
    * @brief Check for non-empty intersection among boxes in container.
    *
    * @return Returns true if there exists any non-empty intersection among
    * the boxes in the container.
    */
   bool
   boxesIntersect() const;

   //@}

   //@{ @name Methods to change or query ordered/unordered state

   /*!
    * @brief Changes state of this container to ordered.
    *
    * If called on a container that is already ordered, nothing changes.
    *
    * @pre each box in container must have valid and unique BoxId
    */
   void
   order();

   /*!
    * @brief Changes state of this container to unordered.
    *
    * This method can be called on any container.
    */
   void
   unorder();

   /*!
    * @brief Return whether this container is ordered.
    *
    * @return  True if ordered, false if unordered.
    */
   constexpr bool
   isOrdered() const
   {
      return d_ordered;
   }

   //@}

   //@{ Methods that may only be called on unordered containers.

   /*!
    * @brief Adds "item" to the "front" of the container.
    *
    * Makes "item" the member of the container that will be returned by
    * front() in an unordered container.
    *
    * @param[in] item
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    */
   void
   pushFront(
      const Box& item)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!empty()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(front(), item);
      }
#endif
      if (!d_ordered) {
         d_list.push_front(item);
      } else {
         TBOX_ERROR("Attempted pushFront on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Adds "item" to the "end" of the container.
    *
    * Makes "item" the member of the container that will be returned by
    * back() in an unordered container.
    *
    * @param[in] item
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    */
   void
   pushBack(
      const Box& item)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!empty()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(front(), item);
      }
#endif
      if (!d_ordered) {
         d_list.push_back(item);
      } else {
         TBOX_ERROR("Attempted pushBack on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Emplaces constructed Box at the front of the container.
    *
    * Forwards arguments to the Box constructor and makes the
    * new Box the member that will be returned by front() in an
    * unordered container.
    *
    * When assertions are enabled and a Box is emplaced into a non-empty
    * container, there will be an assertion error if the dimension of
    * the new Box is not equal to the dimension of the Box at the
    * back of the container.
    *
    * @tparam Args  Valid Box constructor arguments, forwarded to the
    *               constructor
    *
    * @pre !isOrdered()
    */
   template <typename... Args>
   void
   emplaceFront(Args&&... args)
   {
      if (!d_ordered) {
         d_list.emplace_front(std::forward<Args>(args)...);

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_list.size() > 1) {
            TBOX_ASSERT_OBJDIM_EQUALITY2(back(), d_list.front());
         }
#endif
      } else {
         TBOX_ERROR("Attempted emplaceFront on an ordered BoxContainer" << std::endl);
      }

      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Emplaces constructed Box at the end of the container.
    *
    * Forwards arguments to the Box constructor and makes the
    * new Box the member that will be returned by back() in an
    * unordered container.
    *
    * When assertions are enabled and a Box is emplaced into a non-empty
    * container, there will be an assertion error if the dimension of
    * the new Box is not equal to the dimension of the Box at the
    * front of the container.
    *
    * @tparam Args  Valid Box constructor arguments, forwarded to the
    *               constructor
    *
    * @pre !isOrdered()
    */
   template <typename... Args>
   void
   emplaceBack(Args&&... args)
   {
      if (!d_ordered) {
         d_list.emplace_back(std::forward<Args>(args)...);

#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_list.size() > 1) {
            TBOX_ASSERT_OBJDIM_EQUALITY2(front(), d_list.back());
         }
#endif
      } else {
         TBOX_ERROR("Attempted emplaceBack on an ordered BoxContainer" << std::endl);
      }

      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief STL-named version of pushFront().
    */
   void
   push_front(
      const Box& item)
   {
      pushFront(item);
   }

   /*!
    * @brief STL-named version of pushBack().
    */
   void
   push_back(
      const Box& item)
   {
      pushBack(item);
   }

   /*!
    * @brief STL-named version of emplaceFront().
    */
   void
   emplace_front(
      const Box& item)
   {
      pushFront(item);
   }

   /*!
    * @brief STL-named version of emplaceBack().
    */
   void
   emplace_back(
      const Box& item)
   {
      emplaceBack(item);
   }

   /*!
    * @brief Add "item" to specific place in the container.
    *
    * Places "item" immediately before the member of the container pointed
    * to by "iter" in an unordered container.
    *
    * @param[in] iter Location to add item before.
    * @param[in] item Box to add to container.
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    */
   void
   insertBefore(
      iterator iter,
      const Box& item)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!empty()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(front(), item);
      }
#endif
      if (!d_ordered) {
         d_list.insert(iter.d_list_iter, item);
      } else {
         TBOX_ERROR("Attempted insertBefore on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Add "item" to specific place in the container.
    *
    * Places "item" immediately after the member of the container pointed
    * to by "iter" in an unordered container.
    *
    * @param[in] iter Location to add item after.
    * @param[in] item Box to add to container.
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    */
   void
   insertAfter(
      iterator iter,
      const Box& item)
   {
      if (!d_ordered) {
         iterator tmp = iter;
         ++tmp;
         if (tmp == end()) {
            pushBack(item);
         } else {
            insertBefore(tmp, item);
         }
      } else {
         TBOX_ERROR("Attempted insertAfter called on ordered BoxContainer." << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Prepends the Boxes in "boxes" to this BoxContainer.
    *
    * "boxes" will be empty following this operation.
    *
    * @param[in] boxes
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    *
    * @post boxes.empty()
    */
   void
   spliceFront(
      BoxContainer& boxes)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!empty() && !boxes.empty()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(front(), boxes.front());
      }
#endif
      if (!d_ordered) {
         d_list.splice(begin().d_list_iter, boxes.d_list);
      } else {
         TBOX_ERROR("Attempted spliceFront on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Appends the Boxes in "boxes" to this BoxContainer.
    *
    * "boxes" will be empty following this operation.
    *
    * @param[in] boxes
    *
    * @pre empty() || (front().getDim() == item.getDim())
    * @pre !isOrdered()
    *
    * @post boxes.empty()
    */
   void
   spliceBack(
      BoxContainer& boxes)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!empty() && !boxes.empty()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(front(), boxes.front());
      }
#endif
      if (!d_ordered) {
         boxes.spliceFront(*this);
         d_list.swap(boxes.d_list);
      } else {
         TBOX_ERROR("Attempted spliceBack on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Remove the first member of the unordered container.
    *
    * @pre !isOrdered()
    */
   void
   popFront()
   {
      if (!d_ordered) {
         d_list.pop_front();
      } else {
         TBOX_ERROR("Attempted popFront on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Remove the last member of the unordered container.
    *
    * @pre !isOrdered()
    */
   void
   popBack()
   {
      if (!d_ordered) {
         d_list.pop_back();
      } else {
         TBOX_ERROR("Attempted popBack on an ordered BoxContainer" << std::endl);
      }
      if (d_tree) {
         d_tree.reset();
      }
   }

   /*!
    * @brief Place the boxes in the container into a canonical ordering.
    *
    * The canonical ordering for boxes is defined such that boxes that lie
    * next to each other in higher directions are coalesced together before
    * boxes that lie next to each other in lower directions.  This ordering
    * provides a standard representation that can be used to compare box
    * containers.  The canonical ordering also does not allow any overlap
    * between the boxes in the container.  This routine is potentially
    * expensive, since the running time is \f$O(N^2)\f$ for N boxes.  None
    * of the domain calculus routines call simplify(); all calls to simplify
    * the boxes must be explicit.  Note that this routine is distinct from
    * coalesce(), which is not guaranteed to produce a canonical ordering.
    *
    * @pre !isOrdered()
    * @pre empty() || all Boxes in container have same BlockId
    */
   void
   simplify();

   /*!
    * @brief Combine any boxes in the container which may be coalesced.
    *
    * Two boxes may be coalesced if their union is a box (recall that boxes
    * are not closed under index set unions).  Empty boxes in the container
    * are removed during this process.  Note that this is potentially an
    * expensive calculation (e.g., it will require \f$(N-1)!\f$ box
    * comparisons for a box container with \f$N\f$ boxes in the worst
    * possible case).  So this routine should be used sparingly.  Also note
    * that this routine is different than simplify() since it does not
    * produce a canonical ordering.  In particular, this routine processes
    * the boxes in the order in which they appear in the container, rather
    * than attempting to coalesce boxes along specific coordinate directions
    * before others.
    *
    * @pre !isOrdered()
    * @pre empty() || all Boxes in container have same BlockId
    */
   void
   coalesce();

   /*!
    * @brief Rotate boxes in container according to a RotationIdentifier
    *
    * @note Works only in 2D or 3D.
    *
    * @param[in] rotation_ident
    *
    * @pre each Box in container has dim 2 or 3
    * @pre each Box in container has same BlockId
    */
   void
   rotate(
      const Transformation::RotationIdentifier rotation_ident);

   /*!
    * @brief Remove from each box the portions that intersect takeaway.
    *
    * This operation can be thought of as a set difference defined over the
    * abstract AMR box index space.  Performing the set difference will
    * require \f$O(N)\f$ time for a container with \f$N\f$ boxes.  For each
    * box, b, in this container this operation computes b-(b^takeaway) where
    * '^' indicates intersection.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    *
    * @pre !isOrdered()
    */
   void
   removeIntersections(
      const Box& takeaway);

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * For each box, b, in this container and for each box, t, in takeaway
    * this operation computes b-(b^t) where '^' indicates intersection.
    *
    * This only works if all boxes in this BoxContainer and the takeaway
    * BoxContainer have the same BlockId.  An error will occur otherwise.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    *
    * @pre !isOrdered()
    */
   void
   removeIntersections(
      const BoxContainer& takeaway);

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * Uses refinement ratio and grid geometry to handle intersections
    * across block boundaries if needed.
    *
    * @param[in] refinement_ratio  All boxes in this BoxContainer
    * are assumed to exist in index space that has this refinement ratio
    * relative to the coarse-level domain.
    *
    * @param[in] takeaway  The boxes to take away from this BoxContainer.  An
    * error will occur if makeTree with a non-null BaseGridGeometry argument
    * has not been previously called on this container.
    *
    * @param[in] include_singularity_block_neighbors  If true, intersections
    * with neighboring blocks that touch only across an enhanced connectivity
    * singularity will be removed.  If false, those intersections are ignored.
    *
    * @pre !isOrdered()
    * @pre takeaway.hasTree()
    */
   void
   removeIntersections(
      const IntVector& refinement_ratio,
      const BoxContainer& takeaway,
      const bool include_singularity_block_neighbors = false);

   /*!
    * @brief Remove from box the portions intersecting takeaway.
    *
    * This is special version for the case where the container is empty
    * initially.  Upon completion this container contains the result of the
    * removal from box of the intersection of box with takeaway.  If the
    * boxes do not intersect, box is simply added to this container.  This
    * routine is primarily suited for applications which are looking only
    * for the intersection of two boxes.  This operation computes
    * box-(box^takeaway) where '^' indicates intersection.
    *
    * @param[in] box
    * @param[in] takeaway
    *
    * @pre !isOrdered()
    * @pre empty()
    * @pre box.getBlockId() == takeaway.getBlockId()
    */
   void
   removeIntersections(
      const Box& box,
      const Box& takeaway);

   /*!
    * @brief Keep the intersection of the container's boxes and keep.
    *
    * Performing the intersection will require \f$O(N)\f$ time for a
    * container with \f$N\f$ boxes.  The complement of removeIntersections.
    *
    * @param[in] keep
    *
    * @pre !isOrdered()
    */
   void
   intersectBoxes(
      const Box& keep);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * Intersect the boxes in the current container against the boxes in the
    * specified container.  The intersection calculation will require
    * \f$O(N^2)\f$ time for containers with \f$N\f$ boxes.  The complement
    * of removeIntersections.
    *
    * This only works if all boxes in this BoxContainer and the keep
    * BoxContainer have the same BlockId.  An error will occur otherwise.
    *
    * @param[in] keep
    *
    * @pre !isOrdered()
    */
   void
   intersectBoxes(
      const BoxContainer& keep);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * Uses refinement ratio and grid geometry to handle intersections
    * across block boundaries if needed.
    *
    * @param[in]  refinement_ratio  All boxes in this BoxContainer
    * are assumed to exist in index space that has this refinement ratio
    * relative to the coarse-level domain.
    *
    * @param[in] keep  The boxes to intersect with this BoxContainer.  An
    * error will occur if makeTree with a non-null BaseGridGeometry argument
    * has not been previously called on this container.
    *
    * @param[in] include_singularity_block_neighbors  If true, intersections
    * with neighboring blocks that touch only across an enhanced connectivity
    * singularity will be kept.  If false, those intersections are ignored.
    *
    * @pre !isOrdered()
    * @pre keep.hasTree()
    */
   void
   intersectBoxes(
      const IntVector& refinement_ratio,
      const BoxContainer& keep,
      bool include_singularity_block_neighbors = false);

   //@}

   //@{ @name Ordered insert methods

   /*!
    * The insert methods are used to add Boxes to ordered containers.  They
    * may be called on an unordered container only if the size of the unordered
    * container is zero.  If called on such an empty unordered container,
    * the state of the container will be changed to ordered.  A run-time error
    * will occur if called on a non-empty unordered container.
    */

   /*!
    * @brief  Insert a single Box.
    *
    * The Box will be added to the container unless the container already
    * contains a Box with the same BoxId.  If a Box with the same BoxId
    * does already exist in the contianer, the container will not be changed.
    *
    * @return  True if the container did not already have a Box with the
    *          same BoxId, false otherwise.
    *
    * @param[in]  box Box to attempt to insert into the container.
    *
    * @pre box.getBoxId().isValid()
    * @pre empty() || (front().getDim() == box.getDim())
    * @pre empty() || isOrdered()
    */
   bool
   insert(
      const Box& box);

   /*!
    * @brief  Insert a single Box.
    *
    * The Box will be added to the container unless the container already
    * contains a Box with the same BoxId.  If a Box with the same BoxId
    * does already exist in the contianer, the container will not be changed.
    *
    * This version of insert includes an iterator argument pointing somewhere
    * in this container.  This iterator indicates a position in the ordered
    * container where the search for the proper place to insert the given
    * Box will begin.
    *
    * The iterator argument does not determine the place the Box will end up
    * in the ordered container, as that is always determined by BoxId; it
    * is intended only to provide a means of optimization when the calling
    * code knows something about the ordering of the container.
    *
    * @return  iterator pointing to the newly-added Box if the container
    *          did not already have a Box with the same BoxId.  If the
    *          container did have a Box with the same BoxId, the returned
    *          iterator points to that Box.
    *
    * @param[in] position  Location to begin searching for place to insert Box
    * @param[in] box       Box to attempt to insert into the container
    *
    * @pre box.getBoxId().isValid()
    * @pre box.getBlockId() != BlockId::invalidId()
    * @pre empty() || (front().getDim() == box.getDim())
    * @pre empty() || isOrdered()
    */
   iterator
   insert(
      iterator position,
      const Box& box);

   /*!
    * @brief  Insert all Boxes within a range.
    *
    * Boxes in the range [first, last) are added to the ordered container, as
    * long as they do not have a BoxId matching that of a Box already in the
    * container.
    *
    * @param[in] first
    * @param[in] last
    *
    * @pre empty() || isOrdered()
    * @pre for each box in [first, last), box.getBoxId().isValid() &&
    *      (empty || front().getDim() == box.getDim())
    */
   void
   insert(
      const_iterator first,
      const_iterator last);

   //@}

   //@{ @name Methods that may only be called on an ordered container

   /*!
    * @brief  Find a box in an ordered container.
    *
    * Search for a Box having the same BoxId as the given box argument.  This
    * may only be called on an ordered container.
    *
    * @return  If a Box with the same BoxId as the argument is found,
    *          the iterator points to that Box in this container, otherwise
    *          end() for this container is returned.
    *
    * @param[in]  box  Box serving as key for the find operation.  Only
    *                  its BoxId is compared to members of this container.
    *
    * @pre isOrdered()
    */
   iterator
   find(
      const Box& box) const
   {
      if (!d_ordered) {
         TBOX_ERROR("find attempted on unordered BoxContainer." << std::endl);
      }
      iterator iter;
      iter.d_set_iter = d_set.find(const_cast<Box *>(&box));
      iter.d_ordered = true;
      return iter;
   }

   /*!
    * @brief  Get lower bound iterator for a given Box.
    *
    * This may only be called on an ordered container.
    *
    * @return  iterator pointing to the first member of this container
    *          with a BoxId value greater than or equal to the BoxId of
    *          the argument Box.
    *
    * @param[in]  box  Box serving as key for the lower bound search.
    *
    * @pre isOrdered()
    */
   iterator
   lowerBound(
      const Box& box) const
   {
      if (!d_ordered) {
         TBOX_ERROR("lowerBound attempted on unordered BoxContainer." << std::endl);
      }
      iterator iter;
      iter.d_set_iter = d_set.lower_bound(const_cast<Box *>(&box));
      iter.d_ordered = true;
      return iter;
   }

   /*!
    * @brief  Get upper bound iterator for a given Box.
    *
    * @return  iterator pointing to the first member of this container
    *          with a BoxId value greater than the BoxId of the argument
    *          Box.  Will return end() if there are no members with a greater
    *          BoxId value.
    *
    * @param[in]  box  Box serving as key for the upper bound search.
    *
    * @pre isOrdered()
    */
   iterator
   upperBound(
      const Box& box) const
   {
      if (!d_ordered) {
         TBOX_ERROR("upperBound attempted on unordered BoxContainer." << std::endl);
      }
      iterator iter;
      iter.d_set_iter = d_set.upper_bound(const_cast<Box *>(&box));
      iter.d_ordered = true;
      return iter;
   }

   /*!
    * @brief  Erase a Box from the container.
    *
    * If a member of the container has the same BoxId as the argument Box,
    * it will be erased from the container.  If no such member is found,
    * the container is unchanged.
    *
    * @return  1 if a Box is erased, 0 otherwise.
    *
    * @param[in]  box  Box serving as key to find a Box to be erased.
    *
    * @pre isOrdered()
    */
   int
   erase(
      const Box& box);

   // The following may only be called on ordered containers.

   /*!
    * @brief Copy the members of this BoxContainer into two vector<Box>
    * objects, one containing real Boxes and one containing their
    * periodic images.
    *
    * Put the results in the output vectors.  For flexibility and
    * efficiency, the output containers are NOT cleared first, so users
    * may want to clear them before calling this method.
    *
    * @param[out] real_box_vector
    * @param[out] periodic_image_box_vector
    * @param[in] shift_catalog PeriodicShiftCatalog object that maps the
    * PeriodicId to a specific shift.
    *
    * @pre isOrdered()
    */
   void
   separatePeriodicImages(
      std::vector<Box>& real_box_vector,
      std::vector<Box>& periodic_image_box_vector,
      const PeriodicShiftCatalog& shift_catalog) const;

   /*!
    * @brief  Any members of this container that are periodic images will
    *         be erased.
    *
    * @pre empty() || isOrdered()
    */
   void
   removePeriodicImageBoxes();

   /*!
    * @brief  Place unshifted versions of Boxes into a BoxContainer.
    *
    * For all members of this container that are periodic images, create
    * an unshifted copy the member box and add insert it the output container.
    * Additionally, insert all members of the container that are not
    * periodic images to the output container.
    *
    * For flexibility and efficiency, the output container is NOT cleared
    * first, so users may want to clear it before calling this method.
    *
    * @param[out] output_boxes
    *
    * @param[in] refinement_ratio Refinement ratio where the boxes live.
    *
    * @param[in] shift_catalog PeriodicShiftCatalog object that maps the
    * PeriodicId to a specific shift.
    *
    * @pre isOrdered()
    */
   void
   unshiftPeriodicImageBoxes(
      BoxContainer& output_boxes,
      const IntVector& refinement_ratio,
      const PeriodicShiftCatalog& shift_catalog) const;

   //@}

   /*!
    * @brief  Equality operator
    *
    * If the container is ordered, then this checks if all boxes have
    * identical BoxIds and are spatially equal to the boxes in rhs.  If
    * unordered, then this checks only the spatial equality of the boxes.
    *
    * Addtionally to be considered equal, both containers must contain the
    * same number of boxes and have the same ordered/unordered state.
    *
    * @return  true only if the containers meet the equality conditions.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const BoxContainer& rhs) const;

   /*!
    * @brief  Check for equality of BoxIds in ordered containers.
    *
    * An error will occur if this operator is called with either BoxContainer
    * in unordered state.
    *
    * @return  true if both containers are the same size and the boxes
    *          have identical BoxIds.
    *
    * @param[in] other
    * @pre isOrdered() && other.isOrdered()
    *
    */
   bool
   isIdEqual(
      const BoxContainer& other) const;

   /*!
    * @brief  Check for spatial equality of all boxes
    *
    * @return  true if both containers are the same size and the boxes
    *          are spatially equal (same extents and same BlockId).
    *
    * @param[in] other
    */
   bool
   isSpatiallyEqual(
      const BoxContainer& other) const;

   /*!
    * @brief  Inequality operator
    *
    * @return  Return true if operator== would return false.
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const BoxContainer& rhs) const
   {
      return !(*this == rhs);
   }

   //@{ @name I/O

   /*!
    * @brief Write the BoxContainer to a restart database.
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * @brief Read the BoxContainer from a restart database.
    */
   void
   getFromRestart(
      tbox::Database& restart_db);

   /*!
    * @brief Conversion from BoxContainer to std::vector<tbox::DatabaseBox>.
    */
   operator std::vector<tbox::DatabaseBox>() const;

   /*!
    * @brief Print each box in the container to the specified output stream.
    *
    * @param[in] os
    * @param[in] border
    */
   void
   print(
      std::ostream& os = tbox::plog,
      const std::string& border = std::string()) const;

   /*!
    * @brief Intermediary between BoxContainer and output streams,
    * adding ability to control the output.  See
    * BoxContainer::format().
    */
   class Outputter
   {

      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);

private:
      friend class BoxContainer;

      /*!
       * @brief Copy constructor
       */
      Outputter(
         const Outputter& other);

      /*!
       * @brief Construct the Outputter with a BoxContainer and the
       * parameters needed to output the BoxContainer to a stream.
       */
      Outputter(
         const BoxContainer& boxes,
         const std::string& border,
         int detail_depth = 2);

      Outputter&
      operator = (
         const Outputter& rhs);               // Unimplemented private.

      const BoxContainer& d_set;

      const std::string d_border;

      const int d_detail_depth;
   };

   /*!
    * @brief Return a object to that can format the BoxContainer for
    * inserting into output streams.
    *
    * Usage example (printing with a tab indentation):
    * @verbatim
    *    cout << "my boxes:\n" << boxes.format("\t") << endl;
    * @endverbatim
    *
    * @param[in] border Left border of the output
    *
    * @param[in] detail_depth How much detail to print.
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 2) const;

   //@}

   /*!
    * @brief Create a search tree representation of the boxes in this container.
    *
    * If the size of this container is greater than the min_number
    * argument, then an internal BoxTree representation of the boxes will
    * be created and held by the container.
    *
    * This method may only be used if all members of this container have the
    * same BlockId.  An assertion failure will occur if this condition is not
    * met.
    *
    * The building of the tree is an O(N log(N)) operation, but it reduces
    * the cost of the search methods findOverlapBoxes() and hasOverlap() to
    * O(log(N)), rather than O(N).  The tree representation is intended for
    * when these search methods are called multiple times on the same
    * BoxContainer which is not changing, so that the benefit of more
    * efficient search should outweigh the cost of building the tree.
    *
    * The min_number argument is used to indicate a container size below
    * which there is no benefit from building the tree, so no tree is created
    * when the size is less than or equal to this value.
    *
    * If a tree representation has been created via this method, and then
    * any other BoxContainer method is called that changes the container
    * by adding or removing boxes, or by changing the spatial coordinates of
    * the boxes, the tree representation is destroyed.
    *
    * A non-null pointer to a BaseGridGeometry must be provided if this
    * container is going to be used in any of the methods that handle
    * multiblock transformations.  If this container is used only in a
    * single-block context, no BaseGridGeometry argument is necessary.
    *
    * @note The grid_geometry argument is required for multiblock.  It must
    * be the GridGeometry from which the Boxes stored in the container came
    * from.
    *
    * @param[in]  grid_geometry  To handle multiblock transformations if
    *                            needed.
    * @param[in]  min_number
    *
    * @pre min_number > 0
    */
   void
   makeTree(
      const BaseGridGeometry* grid_geometry = 0,
      const int min_number = 10) const;

   /*!
    * @brief Query if the search tree representation exists.
    */
   bool
   hasTree() const
   {
      return d_tree.get() != 0;
   }

   /*!
    * @brief Query if this BoxContainer contains any Box with the given
    * BlockId.
    *
    * This is an efficient query if the the tree representation of this
    * container has been constructed.  If not, it requires a linear search
    * over the entire container and may not be efficient.
    */
   bool
   hasBoxInBlock(
      const BlockId& block_id) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * Every Box in this BoxContainer that intersects with the box argument
    * will be added to the overlap_boxes output container.  The output
    * container will retain the same ordered/unordered state that it had
    * prior to being passed into this method.
    *
    * If this method is used multiple times on the same BoxContainer, it
    * is recommended for efficiency's sake to call makeTree() on this
    * BoxContainer before calling this method.
    *
    * @param[out] overlap_boxes
    *
    * @param[in] box
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * A pointer to every Box in this BoxContainer that intersects with the
    * box argument will be copied to the overlap_boxes output vector.  The
    * vector is not sorted in any way.
    *
    * If this method is used multiple times on the same BoxContainer, it
    * is recommended for efficiency's sake to call makeTree() on this
    * BoxContainer before calling this method.
    *
    * @param[out] overlap_boxes
    *
    * @param[in] box
    *
    * @pre !hasTree() || (d_tree->getNumberBlocksInTree() == 1)
    */
   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box) const;

   /*!
    * @brief Find all boxes that intersect with a given box.
    *
    * Uses refinement ratio and grid geometry to handle intersections
    * across block boundaries if needed.  The makeTree method with a non-null
    * BaseGridGeometry pointer must be called on this container before calling
    * this version of findOverlapBoxes.
    *
    * Every Box in this BoxContainer that intersects with the box argument
    * will be copied to the overlap_boxes output container.  The output
    * container will retain the same ordered/unordered state that it had
    * prior to being passed into this method.
    *
    * @param[out]  overlap_boxes
    *
    * @param[in]  box
    *
    * @param[in]  refinement_ratio  All boxes in this BoxContainer
    * are assumed to exist in index space that has this refinement ratio
    * relative to the coarse-level domain.
    *
    * @param[in]  include_singularity_block_neighbors  If true, intersections
    * with neighboring blocks that touch only across an enhanced connectivity
    * singularity will be added to output.  If false, those intersections are
    * ignored.
    *
    * @pre hasTree()
    */
   void
   findOverlapBoxes(
      BoxContainer& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   void
   findOverlapBoxes(
      std::vector<const Box *>& overlap_boxes,
      const Box& box,
      const IntVector& refinement_ratio,
      bool include_singularity_block_neighbors = false) const;

   /*!
    * @brief Determine if a given box intersects with the BoxContainer.
    *
    * If this method is used multiple times on the same BoxContainer, it
    * is recommended for efficiency's sake to call makeTree() on this
    * BoxContainer before calling this method.
    *
    * This only works if all boxes in this BoxContainer have the same BlockId,
    * and the argument box also has that same BlockId.  An error will occur if
    * these conditions are not met.
    *
    * @return  True if box intersects with any member of the BoxContainer,
    *          false otherwise.
    *
    * @param[in] box
    */
   bool
   hasOverlap(
      const Box& box) const;

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_BOX_CONTAINER_VERSION;

   /*!
    * @brief Remove from each box portions intersecting boxes in takeaway.
    *
    * MultiblockBoxTree has an efficient overlap search method so this
    * version of removeIntersection is relatively fast.
    * For each box, b, in this container and for each box, t, in takeaway
    * this operation computes b-(b^t) where '^' indicates intersection.
    *
    * @param[in] takeaway What to exclude from each box in the container.
    *
    * @pre !isOrdered()
    */
   void
   removeIntersections(
      const MultiblockBoxTree& takeaway);

   /*!
    * @brief Keep the intersection of the container's boxes and keep's boxes
    *
    * MultiblockBoxTree has an efficient overlap search method so this
    * version of intersectBoxes is relatively fast.  The complement of
    * removeIntersections.
    *
    * @param[in] keep
    *
    * @pre !isOrdered()
    */
   void
   intersectBoxes(
      const MultiblockBoxTree& keep);

   /*!
    * @brief Break up bursty against solid and adds the pieces to container.
    *
    * The bursting is done on directions 0 through dimension-1, starting
    * with lowest directions first to try to maintain the canonical
    * representation for the bursted domains.
    *
    * @param[in] bursty
    * @param[in] solid
    * @param[in] direction
    *
    * @pre bursty.getDim() == solid.getDim()
    * @pre dimension <= bursty.getDim().getValue()
    */
   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int direction);

   /*!
    * @brief Break up bursty against solid and adds the pieces to container
    * starting at location pointed to by itr.
    *
    * The bursting is done on directions 0 through dimension-1, starting
    * with lowest directions first to try to maintain the canonical
    * representation for the bursted domains.
    *
    * @param[in] bursty
    * @param[in] solid
    * @param[in] direction
    * @param[in] itr
    *
    * @pre bursty.getDim() == solid.getDim()
    * @pre dimension <= bursty.getDim().getValue()
    */
   void
   burstBoxes(
      const Box& bursty,
      const Box& solid,
      const int direction,
      iterator& itr);

   /*!
    * @brief Remove from each box in the sublist of this container defined
    * by sublist_start and sublist_end portions intersecting takeaway.
    *
    * @param[in] takeaway
    * @param[in] sublist_start
    * @param[in] sublist_end
    * @param[in] insertion_pt Where to put new boxes created by this
    * operation.
    *
    * @pre !isOrdered()
    */
   void
   removeIntersectionsFromSublist(
      const Box& takeaway,
      iterator& sublist_start,
      iterator& sublist_end,
      iterator& insertion_pt);

   /*!
    * List that provides the internal storage for the member Boxes.
    */
   std::list<Box> d_list;

   /*!
    * Set of Box* used for ordered containers.  Each Box* in the set
    * points to a member of d_list.
    */
   std::set<Box *, Box::id_less> d_set;

   bool d_ordered;

   mutable std::shared_ptr<MultiblockBoxTree> d_tree;
};

}
}

#endif
