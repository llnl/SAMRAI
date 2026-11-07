/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/SideIndex.h"

namespace SAMRAI {
namespace pdat {

SideIndex::SideIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

SideIndex::SideIndex(
   const hier::Index& rhs,
   const int axis,
   const int side):
   hier::Index(rhs),
   d_axis(axis)
{
   (*this)(d_axis) += side;
}

SideIndex::SideIndex(
   const SideIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

SideIndex::~SideIndex()
{
}

hier::Index
SideIndex::toCell(
   const int side) const
{
   const tbox::Dimension& dim(getDim());

   hier::Index index(dim);

   for (int i = 0; i < dim.getValue(); ++i) {
      index(i) = (*this)(i);
   }

   index(d_axis) += (side - 1);

   return index;
}

}
}
