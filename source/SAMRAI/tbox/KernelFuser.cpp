/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Kernel fuser
 *
 ************************************************************************/

#include "SAMRAI/tbox/KernelFuser.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"

namespace SAMRAI {
namespace tbox {


KernelFuser::KernelFuser() :
#ifdef HAVE_RAJA
   d_workpool(AllocatorDatabase::getDatabase()->getKernelFuserAllocator()),
   d_workgroup(d_workpool.instantiate()),
   d_worksite(d_workgroup.run()),
#endif
   d_launched(false)
{
}

KernelFuser::~KernelFuser()
{
}

/*
void
KernelFuser::initialize()
{
}
*/

}
}

