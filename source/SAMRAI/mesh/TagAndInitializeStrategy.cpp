/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for params, tagging, init for gridding.
 *
 ************************************************************************/
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>

namespace SAMRAI {
namespace mesh {

TagAndInitializeStrategy::TagAndInitializeStrategy(
   const std::string& object_name):
   d_object_name(object_name)
{
}

TagAndInitializeStrategy::~TagAndInitializeStrategy()
{
}

}
}
