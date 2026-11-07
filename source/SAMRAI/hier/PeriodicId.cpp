/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Periodic shift identifier in periodic domain.
 *
 ************************************************************************/
#include "SAMRAI/hier/PeriodicId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const PeriodicId PeriodicId::s_invalid_id(s_invalid_val);
const PeriodicId PeriodicId::s_zero_id(s_zero_val);


/*
 ******************************************************************************
 ******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const PeriodicId& r)
{
   co << r.d_value;
   return co;
}

}
}
