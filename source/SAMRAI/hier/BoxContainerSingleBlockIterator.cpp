/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"

namespace SAMRAI {
namespace hier {

BoxContainerSingleBlockIterator::BoxContainerSingleBlockIterator(
   const BoxContainer& boxes,
   const BlockId& block_id,
   bool begin):
   d_boxes(&boxes),
   d_block_id(block_id),
   d_iter(begin ? d_boxes->begin() : d_boxes->end())
{
   if (begin) {
      while (d_iter != d_boxes->end() &&
             d_iter->getBlockId() != d_block_id) {
         ++d_iter;
      }
   }
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

BoxContainerSingleBlockIterator&
BoxContainerSingleBlockIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxContainerSingleBlockIterator
BoxContainerSingleBlockIterator::operator ++ (
   int)
{
   BoxContainerSingleBlockIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_boxes->end() &&
            d_iter->getBlockId() != d_block_id);
   return saved;
}

int
BoxContainerSingleBlockIterator::count() const
{
   int ct = 0;
   BoxContainerSingleBlockIterator iter(d_boxes->begin(d_block_id));
   while (iter != d_boxes->end(d_block_id)) {
      ++ct;
      ++iter;
   }
   return ct;
}

}
}
