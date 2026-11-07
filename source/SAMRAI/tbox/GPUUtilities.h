/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for GPU
 *
 ************************************************************************/

#ifndef included_tbox_GPUUtilities
#define included_tbox_GPUUtilities

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_RAJA)
#include "RAJA/RAJA.hpp"
#endif

namespace SAMRAI {
namespace tbox {


/*!
 * Utility functions to support running on GPUs.
 */

struct GPUUtilities {

/*!
 * @brief Manually set flag telling whether the run is using GPUs.
 *
 * Usually the default value of the flag will be correct. This may be used
 * to set the flag to false if running a GPU-enabled build without access to
 * GPUs.
 */
static void
setUsingGPU(bool using_gpu);

static bool
isUsingGPU()
{
   return s_using_gpu;
}

/*!
 * @brief Synchronizes GPU threads.
 *
 * May be called after RAJA-CUDA loops to ensure data modified inside
 * kernels is ready for next operation.
 */
static void
parallel_synchronize()
{
#if defined(HAVE_CUDA) && defined(HAVE_RAJA)
   if (s_using_gpu) {
      RAJA::synchronize<RAJA::cuda_synchronize>();
   }
#elif defined(HAVE_HIP) && defined(HAVE_RAJA)
   if (s_using_gpu) {
      RAJA::synchronize<RAJA::hip_synchronize>();
   }
#endif       
}

private:

   /*!
    * Flag indicating wheter the run is using GPUs. This has a default of
    * true when SAMRAI is compiled for GPUs, false otherwise.
    */
   static bool s_using_gpu;

};

}
}

#endif
