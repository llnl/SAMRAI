/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleBlockIterator
#define included_hier_BoxContainerSingleBlockIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief BoxContainer iterator picking items with a specified
 * BlockId.
 *
 * This iterator runs through all Boxes in a BoxContainer that
 * has the given BlockId.  The iterator runs through the Boxes
 * in the order they appear in the BoxContainer, skipping over
 * Boxes that do not have the specified owner rank.
 */
class BoxContainerSingleBlockIterator
{
   friend class BoxContainer;

public:

   /*!
    * @brief Copy constructor
    */
   BoxContainerSingleBlockIterator(
      const BoxContainerSingleBlockIterator& other) = default;

   //! @brief Destructor.
   ~BoxContainerSingleBlockIterator() = default;

   /*!
    * @brief Assignment operator.
    */
   BoxContainerSingleBlockIterator&
   operator = (
      const BoxContainerSingleBlockIterator& r) = default;

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
      const BoxContainerSingleBlockIterator& r) const
   {
      return d_boxes == r.d_boxes &&
             d_block_id == r.d_block_id &&
             d_iter == r.d_iter;
   }

   /*!
    * @brief Inequality comparison.
    */
   bool
   operator != (
      const BoxContainerSingleBlockIterator& r) const
   {
      return d_boxes != r.d_boxes ||
             d_block_id != r.d_block_id ||
             d_iter != r.d_iter;
   }

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increment the iterator and returns the incremented
    * state.
    */
   BoxContainerSingleBlockIterator&
   operator ++ ();

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the iterator, increment it and returns the
    * saved iterator.
    */
   BoxContainerSingleBlockIterator
   operator ++ (
      int);

   /*!
    * @brief Returns the number of BoxContainers being iterated through.
    */
   int
   count() const;

private:
   /*!
    * @brief Constructor
    *
    * @param [in] container
    * @param [in] block_id
    * @param [in] begin
    */
   BoxContainerSingleBlockIterator(
      const BoxContainer& container,
      const BlockId& block_id,
      bool begin);

   /*!
    * @brief BoxContainer being iterated through.
    */
   const BoxContainer* d_boxes;

   /*!
    * @brief The BlockId.
    */
   BlockId d_block_id;

   /*!
    * @brief The iterator.
    */
   BoxContainer::const_iterator d_iter;

};

}
}

#endif  // included_hier_BoxContainerSingleBlockIterator
