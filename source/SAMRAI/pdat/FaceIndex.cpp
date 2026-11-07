/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/FaceIndex.h"

namespace SAMRAI {
namespace pdat {

FaceIndex::FaceIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

FaceIndex::FaceIndex(
   const hier::Index& rhs,
   const int axis,
   const int face):
   hier::Index(rhs),
   d_axis(axis)
{
   (*this)(0) = rhs(d_axis) + face;
   for (int i = 1; i < getDim().getValue(); ++i) {
      (*this)(i) = rhs((d_axis + i) % getDim().getValue());
   }
}

FaceIndex::FaceIndex(
   const FaceIndex& rhs):
   hier::Index(rhs),
   d_axis(rhs.d_axis)
{
}

FaceIndex::~FaceIndex()
{
}

hier::Index
FaceIndex::toCell(
   const int face) const
{
   hier::Index index(getDim());
   index(d_axis) = (*this)(0) + face - 1;
   for (int i = 1; i < getDim().getValue(); ++i) {
      index((d_axis + i) % getDim().getValue()) = (*this)(i);
   }
   return index;
}

}
}
