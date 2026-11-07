/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*
 ******************************************************************************
 * Stream-insert operator.
 ******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const BoxId& r)
{
   co << r.d_global_id.getOwnerRank()
   << '#' << r.d_global_id.getLocalId()
   << '/' << r.d_periodic_id;
   return co;
}

}
}
