/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch classes
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchFactory.h"


namespace SAMRAI {
namespace hier {

std::shared_ptr<Patch>
PatchFactory::allocate(
   const Box& box_level_box,
   const std::shared_ptr<PatchDescriptor>& descriptor) const
{
   return std::make_shared<Patch>(box_level_box, descriptor);
}

}
}
