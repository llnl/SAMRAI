/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory abstract base class for creating patch data objects
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchDataFactory.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

PatchDataFactory::PatchDataFactory(
   const IntVector& ghosts):
   d_ghosts(ghosts)
{
   TBOX_ASSERT(ghosts.min() >= 0);
}

/**********************************************************************
* Default implementation
**********************************************************************/

MultiblockDataTranslator *
PatchDataFactory::getMultiblockDataTranslator()
{
   return 0;
}

}
}
