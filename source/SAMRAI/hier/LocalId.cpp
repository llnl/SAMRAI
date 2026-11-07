/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Generic identifier used on a single process.
 *
 ************************************************************************/
#include "SAMRAI/hier/LocalId.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const LocalId LocalId::s_invalid_id(s_invalid_val);
const LocalId LocalId::s_zero_id(s_zero_val);

/*
 *******************************************************************************
 *******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const LocalId& r)
{
   if (r.isValid()) {
      co << r.d_value;
   } else {
      co << 'X';
   }
   return co;
}

}
}
