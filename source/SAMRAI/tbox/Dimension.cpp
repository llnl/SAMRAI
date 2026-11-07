/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Dimension class for abstracting dimension
 *
 ************************************************************************/
#include "SAMRAI/tbox/Dimension.h"

namespace SAMRAI {
namespace tbox {

std::ostream&
operator << (
   std::ostream& s,
   const Dimension& dim)
{
   s << dim.getValue() << 'D';
   return s;
}

}
}
