/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Iterator for cell centered patch data types
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellIterator.h"

namespace SAMRAI {
namespace pdat {

CellIterator::CellIterator(
   const hier::Box& box,
   bool begin):
   d_index(box.lower()),
   d_box(box)
{
   if (!d_box.empty() && !begin) {
      d_index(d_box.getDim().getValue() - 1) =
         d_box.upper(static_cast<tbox::Dimension::dir_t>(d_box.getDim().getValue() - 1)) + 1;
   }
}

CellIterator::CellIterator(
   const CellIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

CellIterator::~CellIterator()
{
}

CellIterator&
CellIterator::operator ++ ()
{
   ++d_index(0);
   for (tbox::Dimension::dir_t i = 0; i < d_box.getDim().getValue() - 1; ++i) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         ++d_index(i + 1);
      } else {
         break;
      }
   }
   return *this;
}

CellIterator
CellIterator::operator ++ (
   int)
{
   CellIterator tmp = *this;
   ++d_index(0);
   for (tbox::Dimension::dir_t i = 0; i < d_box.getDim().getValue() - 1; ++i) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         ++d_index(i + 1);
      } else {
         break;
      }
   }
   return tmp;
}

}
}
