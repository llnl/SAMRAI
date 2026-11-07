/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/EdgeIndex.h"

namespace SAMRAI {
namespace pdat {

EdgeIndex::EdgeIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

EdgeIndex::EdgeIndex(
   const hier::Index& rhs,
   const int axis,
   const int edge):
   hier::Index(rhs),
   d_axis(axis)
{
   if (getDim() > tbox::Dimension(1)) {
      (*this)((d_axis + 1) % getDim().getValue()) += edge % 2;
   }
   for (int j = 2; j < getDim().getValue(); ++j) {
      (*this)((d_axis + j) % getDim().getValue()) += (edge / (1 << (j - 1))) % 2;
   }
}

EdgeIndex::EdgeIndex(
   const EdgeIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

EdgeIndex::~EdgeIndex()
{
}

hier::Index
EdgeIndex::toCell(
   const int edge) const
{
   const tbox::Dimension& dim(getDim());
   hier::Index index(dim);

   for (int i = 0; i < dim.getValue(); ++i) {
      index(i) = (*this)(i);
   }

   if (dim > tbox::Dimension(1)) {
      index((d_axis + 1) % dim.getValue()) += ((edge % 2) - 1);
   }
   for (int j = 2; j < dim.getValue(); ++j) {
      index((d_axis + j) % dim.getValue()) += (((edge / (2 << (j - 1))) % 2) - 1);
   }
   return index;
}

}
}
