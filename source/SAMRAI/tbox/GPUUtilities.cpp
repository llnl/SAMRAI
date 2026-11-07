/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for GPUs
 *
 ************************************************************************/

#include "SAMRAI/tbox/GPUUtilities.h"


namespace SAMRAI {
namespace tbox {

bool GPUUtilities::s_using_gpu =
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
   true;
#else
   false;
#endif

void GPUUtilities::setUsingGPU(bool using_gpu)
{
   s_using_gpu = using_gpu;
}



}
}

