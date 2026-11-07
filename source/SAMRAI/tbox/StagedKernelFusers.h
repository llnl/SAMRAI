/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton container of kernel fusers
 *
 ************************************************************************/

#ifndef included_tbox_StagedKernelFusers
#define included_tbox_StagedKernelFusers

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/KernelFuser.h"

#ifdef HAVE_RAJA
#include "RAJA/RAJA.hpp"
#endif

#include <map>

namespace SAMRAI {
namespace tbox {
/*!
 * @brief Class StagedKernelFusers is a singleton class that serves as
 * a container of KernelFuser objects that can be executed in sequence. 
 *
 * The singleton object is accessed by the getInstance() method, and the
 * KernelFuser objects are ordered by the integer value used in the
 * getFuser() method, so that their operations can be executed in ordered
 * "stages".
 *
 * A typical usage pattern would be a loop where there is a series of
 * kernels which are dependent on each other in the sequence they appear in
 * the loop, but the kernels for each iteration are independent and
 * thus could benefit from usage of KernelFuser. This class allows kernels
 * to be fused in different KernelFuser instances while maintaining the
 * required order of execution.
 *
 *   for loop over patches
 *      staged_fusers->enqueue(0, begin, end, kernel);
 *      staged_fusers->enqueue(1, begin, end, kernel);
 *      staged_fusers->enqueue(2, begin, end, kernel);
 *   end loop
 *
 *   staged_fusers->launchAndCleanup();
 *
 * In the above pseudocode, launchAndCleanup will launch the fused kernel
 * from each stage with a synchronization and cleanup to ensure that the
 * object can be used again.
 */

class StagedKernelFusers
{
public:

   /*!
    * @brief Get a pointer to the singleton instance of this object.
    */
   static StagedKernelFusers* getInstance();

   /*!
    * @brief Enqueue a kernel for fusion.
    *
    * @param stage   The stage value to control order of kernel fuser launches
    * @param begin   Beginning loop index value for the threaded kernel
    * @param end     Ending loop index value for the threaded kernel
    * @param kernel  The kernel containing the lambda function.
    */
#ifdef HAVE_RAJA
   template<typename Kernel>
   void enqueue(int stage, int begin, int end, Kernel&& kernel) {
      d_kernel_fusers[stage].enqueue(begin, end, kernel);
      d_enqueued = true;
   }
#endif

   /*!
    * @brief Launch all enqueued kernels in all stages
    *
    * The KernelFuser for each stage is launched in order.  Note that there
    * are no synchronizations between launches, so this should only be
    * called if all kernels are independent and there is no danger of
    * race conditions.
    */
   void launch()
   {
      for (auto& fuser : d_kernel_fusers) {
         fuser.second.launch();
         d_launched = (d_launched || fuser.second.launched()); 
      }
   }

   /*!
    * @brief Launch all enqueued kernels in all stages
    *
    * The KernelFuser for each stage is launched in order, with a 
    * synchronization after each launch, so that the kernels from one stage
    * will be complete before the launch of kernels from the next stage in
    * order.
    */
   void launchWithSynchronize()
   {
      for (auto& fuser : d_kernel_fusers) {
         fuser.second.launch();
#ifdef HAVE_RAJA
         if (fuser.second.launched()) {
            tbox::parallel_synchronize();
         }
#endif
         d_launched = (d_launched || fuser.second.launched()); 
      }
   }

   /*!
    * @brief Clean up the staged kernel fusers
    *
    * This clears all of the internal data in the KernelFusers. If this object
    * has previously called launch without a synchronization, a host-device
    * synchronization should be called before cleanup to avoid race conditions.
    * After cleanup, this object can be re-used with more kernels to be
    * enqueued and launched.
    */
   void cleanup()
   {
      for (auto& fuser : d_kernel_fusers) {
         fuser.second.cleanup();
      }
      d_enqueued = false;
      d_launched = false;
   }

   /*!
    * @brief Get a pointer to the KernelFuser object for a given stage.
    */
   KernelFuser* getFuser(int stage)
   {
      return &d_kernel_fusers[stage];
   }

   /*!
    * @brief Clear the KernelFuser object for the given stage.
    *
    * All data for this stage, including any enqueued kernels, will be
    * erased.  Be sure to only call this when there are no remaining
    * operations needed for this stage.
    */
   void clearKernelFuser(int stage)
   {
      d_kernel_fusers.erase(stage);
   }

   /*!
    * @brief Clear all staged kernel fusers.
    *
    * This returns the singleton StagedKernelFusers object to its
    * post-initialization state.
    */
   void clearAllFusers()
   {
      d_kernel_fusers.clear();
      d_launched = false;
      d_enqueued = false;
   }

   /*!
    * @brief Return true if this object has launched kernels and has not
    * cleaned up.
    */
   bool launched()
   {
      return d_launched;
   }

   /*!
    * @brief Return true if this object has any kernels that have been
    * enqueued but not launched.
    */
   bool enqueued()
   {
      return d_enqueued;
   }

   /*!
    * @brief Initialize the singleton object.
    */
   void initialize();

   /*!
    * @brief Launch all kernel fuser in staged order and clean up.
    *
    * A host-device synchronization will occur after the launch of each
    * staged fuser.
    */
   void launchAndCleanup()
   {
      launchWithSynchronize();
      cleanup();
   }

protected:

   /*!
    * @brief Protected constructor for singleton class.
    */
   StagedKernelFusers();

   /*!
    * @brief Protected destructor for singleton class.
    */
   virtual ~StagedKernelFusers();

private:

   static void startupCallback();
   static void shutdownCallback();

   static StagedKernelFusers* s_staged_kernel_fusers_instance;

   static StartupShutdownManager::Handler
   s_startup_handler;

   std::map<int, KernelFuser> d_kernel_fusers;

   bool d_launched = false;
   bool d_enqueued = false;

};

}
}

#endif
