/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to user-defined operations used in FAC solve.
 *
 ************************************************************************/
#include "SAMRAI/solv/FACOperatorStrategy.h"

namespace SAMRAI {
namespace solv {

FACOperatorStrategy::FACOperatorStrategy()
{
}

FACOperatorStrategy::~FACOperatorStrategy()
{
}

void
FACOperatorStrategy::deallocateOperatorState()
{
}

}
}
