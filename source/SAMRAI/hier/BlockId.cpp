/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const BlockId BlockId::s_invalid_id(s_invalid_val);
const BlockId BlockId::s_zero_id(s_zero_val);


}
}
