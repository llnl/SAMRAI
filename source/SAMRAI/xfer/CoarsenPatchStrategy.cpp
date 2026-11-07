/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data.
 *
 ************************************************************************/
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"

namespace SAMRAI {
namespace xfer {

CoarsenPatchStrategy::CoarsenPatchStrategy()
{
   registerObject();
}

CoarsenPatchStrategy::~CoarsenPatchStrategy()
{
}

/*
 *************************************************************************
 * Compute the max coarsen stencil width from all constructed
 * coarsen patch strategies.
 *************************************************************************
 */
hier::IntVector
CoarsenPatchStrategy::getMaxCoarsenOpStencilWidth(
   const tbox::Dimension& dim)
{
   hier::IntVector max_width(dim, 0);

   std::set<CoarsenPatchStrategy *>& current_objects =
      CoarsenPatchStrategy::getCurrentObjects();
   for (std::set<CoarsenPatchStrategy *>::const_iterator
        si = current_objects.begin(); si != current_objects.end(); ++si) {
      const CoarsenPatchStrategy* strategy = *si;
      max_width.max(strategy->getCoarsenOpStencilWidth(dim));
   }

   return max_width;
}

}
}
