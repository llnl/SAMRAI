/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Kernel fuser
 *
 ************************************************************************/

#ifndef included_tbox_KernelFuser
#define included_tbox_KernelFuser

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_UMPIRE
#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"
#endif

#include "SAMRAI/tbox/ExecutionPolicy.h"
#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#ifdef HAVE_RAJA
#include "RAJA/RAJA.hpp"
#endif


#ifdef HAVE_UMPIRE
#include "umpire/Allocator.hpp"
#include "umpire/TypedAllocator.hpp"
#endif

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class KernelFuser is a class that is used to encapsulate RAJA
 * features for fusing the launch of multiple RAJA loop kernels into a
 * single kernel launch on the GPU.
 *
 * The fuser works by adding any number of lambda functions that can run as
 * independent kernels to a pool using the enqueue method, then calling the
 * launch method to execute all of the kernels with one launch.  To re-use
 * a KernelFuser object after a launch, a host-device synchronization
 * and then the cleanup method must be called.
 */

class KernelFuser
{
public:

   /*!
    * @brief Constructor
    */
   KernelFuser();

   /*!
    * @brief Destructor
    */
   virtual ~KernelFuser();

   /*!
    * @brief Enqueue a kernel for kernel fusion.
    *
    * @param  begin   Beginning loop index value for the threaded kernel
    * @param  end     Ending loop index value for the threaded kernel
    * @param  kernel  The kernel containing the lambda function.
    */
#ifdef HAVE_RAJA
   template<typename Kernel>
   void enqueue(int begin, int end, Kernel&& kernel) {
      if (d_launched) {
         TBOX_ERROR("KernelFuser Error: Cannont enqueue until cleanup called after previous launch.");
      }

      d_workpool.enqueue(RAJA::RangeSegment(begin, end), std::forward<Kernel>(kernel));
   }
#endif

   /*!
    * @brief Launch all enqueued kernels
    *
    * When running on GPU systems all enqueued kernels will run concurrently.
    */
   void launch()
   {
      if (d_launched) {
         TBOX_ERROR("KernelFuser Error: This KernelFuser already launched.");
      }
 
#ifdef HAVE_RAJA
      if (d_workpool.num_loops() > 0) {
         d_workgroup = d_workpool.instantiate();
         d_worksite = d_workgroup.run();
      }
      d_launched = true;
 #endif
   }

   /*!
    * @brief Clean up the kernel fuser
    *
    * This clears all of the internal data in this class.  If this object
    * has previously called launch, a host-device synchronization should
    * be called before cleanup to avoid race conditions.  After cleanup,
    * this object can be re-used with more kernels to be enqueued and launched.
    */
   void cleanup()
   {
#ifdef HAVE_RAJA
      d_workpool.clear();
      d_workgroup.clear();
      d_worksite.clear();
      d_launched = false;
#endif
   }

   /*!
    * @brief Returns true if this object has launched and has not
    * subsequently called cleanup
    */
   bool launched() const
   {
      return d_launched;
   }

   /*!
    * @brief Returns true if this objects has enqueued kernels
    */
   bool enqueued() const
   { 
#ifdef HAVE_RAJA
      return (d_workpool.num_loops() > 0);
#else
      return false;
#endif
   }


private:

#ifdef HAVE_UMPIRE
   using Allocator = umpire::TypedAllocator<char>;
#else
   using Allocator = ResourceAllocator;
#endif

#ifdef HAVE_RAJA
   using Policy = typename tbox::detail::policy_traits< tbox::policy::parallel >::WorkGroupPolicy;
   using WorkPool  = RAJA::WorkPool <Policy, int, RAJA::xargs<>, Allocator>;
   using WorkGroup = RAJA::WorkGroup<Policy, int, RAJA::xargs<>, Allocator>;
   using WorkSite  = RAJA::WorkSite <Policy, int, RAJA::xargs<>, Allocator>;
#endif

#ifdef HAVE_RAJA
   WorkPool d_workpool;
   WorkGroup d_workgroup;
   WorkSite d_worksite;
#endif

   bool d_launched;
};

}
}

#endif 
