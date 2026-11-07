/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   base class for fuseable schedule transactions
 *
 ************************************************************************/


#include "SAMRAI/tbox/TransactionFuseable.h"

namespace SAMRAI {
namespace tbox {

TransactionFuseable::TransactionFuseable()
{
}

TransactionFuseable::~TransactionFuseable()
{
}

void
TransactionFuseable::setKernelFuser(StagedKernelFusers* fuser)
{
   d_fuser = fuser;
}

StagedKernelFusers*
TransactionFuseable::getKernelFuser()
{
   return d_fuser;
}

}
}
