/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Describes boundaries for a patch
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchBoundaries.h"

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Constructor leaves the arrays empty.
 *
 *************************************************************************
 */
PatchBoundaries::PatchBoundaries(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_array_of_bboxes(dim.getValue())
{
}

} // SAMRAI namespace
} // hier namespace
