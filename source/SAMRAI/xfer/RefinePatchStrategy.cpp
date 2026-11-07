/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

RefinePatchStrategy::RefinePatchStrategy()
{
   registerObject();
}

RefinePatchStrategy::~RefinePatchStrategy()
{
   unregisterObject();
}

/*
 *************************************************************************
 * Compute the max refine stencil width from all constructed
 * refine patch strategies.
 *************************************************************************
 */
hier::IntVector
RefinePatchStrategy::getMaxRefineOpStencilWidth(
   const tbox::Dimension& dim)
{
   hier::IntVector max_width(dim, 0);

   std::set<RefinePatchStrategy *>& current_objects =
      RefinePatchStrategy::getCurrentObjects();
   for (std::set<RefinePatchStrategy *>::const_iterator
        si = current_objects.begin(); si != current_objects.end(); ++si) {
      const RefinePatchStrategy* strategy = *si;
      max_width.max(strategy->getRefineOpStencilWidth(dim));
   }

   return max_width;
}

}
}
