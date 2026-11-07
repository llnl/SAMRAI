/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleOwnerIterator
#define included_hier_BoxContainerSingleOwnerIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief BoxContainer iterator picking items with a specified
 * owner rank.
 *
 * This iterator runs through all Boxes in a BoxContainer that
 * has the given owner rank.  The iterator runs through the
 * Boxes in the order they appear in the BoxContainer, skipping
 * over Boxes that do not have the specified owner rank.
 */
class BoxContainerSingleOwnerIterator
{
   friend class BoxContainer;

public:

   /*!
    * @brief Copy constructor
    */
   BoxContainerSingleOwnerIterator(
      const BoxContainerSingleOwnerIterator& other) = default;

   //! @brief Destructor.
   ~BoxContainerSingleOwnerIterator() = default;

   /*!
    * @brief Assignment operator.
    */
   BoxContainerSingleOwnerIterator&
   operator = (
      const BoxContainerSingleOwnerIterator& r) = default;

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const Box&
   operator * () const
   {
      return *d_iter;
   }

   /*!
    * @brief Dereference operator mimicking a pointer dereference.
    */
   const Box *
   operator -> () const
   {
      return &(*d_iter);
   }

   /*!
    * @brief Equality comparison.
    */
   bool
   operator == (
      const BoxContainerSingleOwnerIterator& r) const
   {
      return d_boxes == r.d_boxes &&
             d_owner_rank == r.d_owner_rank &&
             d_iter == r.d_iter;
   }

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const BoxContainerSingleOwnerIterator& r) const
   {
      return d_boxes != r.d_boxes ||
             d_owner_rank != r.d_owner_rank ||
             d_iter != r.d_iter;
   }

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   BoxContainerSingleOwnerIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   BoxContainerSingleOwnerIterator
   operator ++ (
      int);

private:
   /*!
    * @brief Constructor
    *
    * @param [in] container
    * @param [in] owner_rank
    * @param [in] begin
    */
   BoxContainerSingleOwnerIterator(
      const BoxContainer& container,
      const int& owner_rank,
      bool begin);

   /*!
    * @brief BoxContainer being iterated through.
    */
   const BoxContainer* d_boxes;

   /*!
    * @brief The owner_rank.
    */
   int d_owner_rank;

   /*!
    * @brief The iterator.
    */
   BoxContainer::const_iterator d_iter;

};

}
}

#endif  // included_hier_BoxContainerSingleOwnerIterator
