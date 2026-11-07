/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellIndex.h"

namespace SAMRAI {
namespace pdat {

CellIndex::CellIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
}

CellIndex::CellIndex(
   const hier::Index& rhs):hier::Index(rhs)
{
}

CellIndex::CellIndex(
   const CellIndex& rhs):hier::Index(rhs)
{
}

CellIndex::~CellIndex()
{
}

}
}
