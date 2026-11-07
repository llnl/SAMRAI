/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class that describes intersections between AMR boxes
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxOverlap.h"

namespace SAMRAI {
namespace hier {

void
BoxOverlap::print(
   std::ostream& os) const
{
   os << "print() method not implemented for this overlap type" << std::endl;
}

}
}
