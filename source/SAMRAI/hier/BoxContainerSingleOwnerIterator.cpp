/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"

namespace SAMRAI {
namespace hier {

BoxContainerSingleOwnerIterator::BoxContainerSingleOwnerIterator(
   const BoxContainer& boxes,
   const int& owner_rank,
   bool begin):
   d_boxes(&boxes),
   d_owner_rank(owner_rank),
   d_iter(begin ? d_boxes->begin() : d_boxes->end())
{
   if (begin) {
      while (d_iter != d_boxes->end() &&
             d_iter->getOwnerRank() != d_owner_rank) {
         ++d_iter;
      }
   }
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator&
BoxContainerSingleOwnerIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator
BoxContainerSingleOwnerIterator::operator ++ (
   int)
{
   BoxContainerSingleOwnerIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return saved;
}

}
}
