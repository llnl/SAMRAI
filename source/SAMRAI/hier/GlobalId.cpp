/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Globally unique identifier that can be locally determined.
 *
 ************************************************************************/
#include "SAMRAI/hier/GlobalId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

std::ostream&
operator << (
   std::ostream& co,
   const GlobalId& r)
{
   co << r.d_owner_rank << '#' << r.d_local_id;
   return co;
}

}
}
