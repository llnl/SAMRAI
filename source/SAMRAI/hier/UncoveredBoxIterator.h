/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_UncoveredBoxIterator
#define included_hier_UncoveredBoxIterator

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BoxContainer.h"

#include <utility>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

class Patch;
class PatchHierarchy;
class FlattenedHierarchy;

/*!
 * @brief An iterator over the uncovered Boxes in a hierarhcy.
 *
 * This iterator points to a pair consisting of the current uncovered Box and
 * the Patch with which it is associated.  Note that in the case of overlapping
 * Patches in which the overlap is uncovered, the overlapping region will
 * appear multiple times in the iteration.  The number of appearances is equal
 * to the number of Patches that overlap that region.
 */
class UncoveredBoxIterator
{
   friend class PatchHierarchy;
   friend class FlattenedHierarchy;

public:
   /*!
    * @brief Copy constructor.
    *
    * @param[in] other The iterator being copied.
    */
   UncoveredBoxIterator(
      const UncoveredBoxIterator& other);

   /*!
    * Destructor.
    */
   ~UncoveredBoxIterator();

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs The right hand side of the assignment.
    */
   UncoveredBoxIterator&
   operator = (
      const UncoveredBoxIterator& rhs);

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const std::pair<std::shared_ptr<Patch>, Box>&
   operator * () const;

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const std::pair<std::shared_ptr<Patch>, Box> *
   operator -> () const;

   /*!
    * @brief Equality comparison.
    *
    * @param[in] rhs The right hand side of the comparison.
    */
   bool
   operator == (
      const UncoveredBoxIterator& rhs) const;

   /*!
    * @brief Inequality comparison.
    *
    * @param[in] rhs The right hand side of the comparison.
    */
   bool
   operator != (
      const UncoveredBoxIterator& rhs) const;

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   UncoveredBoxIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increments it and returns the
    * saved iterator.
    */
   UncoveredBoxIterator
   operator ++ (
      int);

private:
   /*!
    * @brief Unimplemented default constructor.
    */
   UncoveredBoxIterator();

   /*!
    * @brief Constructor.
    *
    * @param[in] hierarchy The hierarchy over which the iteration will occur.
    * @param[in] begin If true the iteration starts from the beginning.
    */
   UncoveredBoxIterator(
      const PatchHierarchy* hierarchy,
      bool begin);

   UncoveredBoxIterator(
      const FlattenedHierarchy* hierarchy,
      bool begin);

   /*!
    * @brief Private method to do work common to both pre and post increments.
    */
   void
   incrementIterator();

   /*!
    * @brief Private method to find first uncovered box
    */
   void
   findFirstUncoveredBox();

   /*!
    * @brief Set the item based on the current patch id and uncovered box.
    */
   void
   setIteratorItem();

   /* The PatchHierarchy on which this iterator operates. */
   const PatchHierarchy* d_hierarchy;
   const FlattenedHierarchy* d_flattened_hierarchy;

   bool d_allocated_flattened_hierarchy;

   /* The current level in the PatchHierarchy. */
   int d_level_num;

   /* The id of the current patch */
   BoxId d_current_patch_id;

   /* The iterator over the uncovered boxes for the current patch. */
   BoxContainer::const_iterator d_uncovered_boxes_itr;

   /* The iterator at the end of the uncovered boxes for the current patch. */
   BoxContainer::const_iterator d_uncovered_boxes_itr_end;

   /* The current item in the iteration. */
   std::pair<std::shared_ptr<Patch>, Box>* d_item;

   /* The number of the finest level in the hierarchy. */
   int d_finest_level_num;

};

}
}

#endif
