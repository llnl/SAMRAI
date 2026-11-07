/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real Boxes in a BoxContainer.
 *
 ************************************************************************/
#include "SAMRAI/hier/RealBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

RealBoxConstIterator::RealBoxConstIterator(
   const BoxContainer& boxes,
   bool begin):
   d_boxes(&boxes),
   d_ni(begin ? d_boxes->begin() : d_boxes->end())
{
   if (begin) {
      while (d_ni != d_boxes->end() && d_ni->isPeriodicImage()) {
         ++d_ni;
      }
   }
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator&
RealBoxConstIterator::operator ++ ()
{
   do {
      ++d_ni;
   } while (d_ni != d_boxes->end() && d_ni->isPeriodicImage());
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator
RealBoxConstIterator::operator ++ (
   int)
{
   RealBoxConstIterator saved = *this;
   do {
      ++d_ni;
   } while (d_ni != d_boxes->end() && d_ni->isPeriodicImage());
   return saved;
}

}
}
