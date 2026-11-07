/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An restart manager singleton class
 *
 ************************************************************************/

#ifndef included_hier_PatchDataRestartManager
#define included_hier_PatchDataRestartManager

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

/**
 * Class PatchDataRestartManager handles registration of PatchData for restart.
 */

class PatchDataRestartManager
{
public:
   /**
    * Return a pointer to the single instance of the PatchDataRestartManager.
    * All access to the PatchDataRestartManager object is through getManager().
    *
    * Note that when the manager is accessed for the first time, the
    * Singleton instance is registered with the StartupShutdownManager
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    */
   static PatchDataRestartManager *
   getManager();

   /**
    * @brief Check whether given patch data index is registered with the
    * manager for restart.
    *
    * @param[in]  index  Integer patch data index to check.
    *
    * @return Boolean true if the patch data with the given index
    *         is registered for restart; otherwise false.
    *
    */
   bool
   isPatchDataRegisteredForRestart(
      int index) const
   {
      return d_patchdata_restart_table.isSet(index);
   }

   /**
    * Registers a patch data index for restart.
    *
    * @param[in]  index  Integer patch data index to check.
    */
   void
   registerPatchDataForRestart(
      int index)
   {
      d_patchdata_restart_table.setFlag(index);
   }

   /**
    * Unregisters a patch data index for restart.
    *
    * @param[in]  index  Integer patch data index to check.
    */
   void
   unregisterPatchDataForRestart(
      int index)
   {
      d_patchdata_restart_table.clrFlag(index);
   }

   /**
    * Returns true if the patch data components selected by "selected" are
    * identical to the components registered in the manager and vice-versa.
    *
    * @param[in]  selected  The patch data components that have been selected.
    */
   bool
   registeredPatchDataMatches(const ComponentSelector& selected)
   {
      return selected == d_patchdata_restart_table;
   }

private:
   /**
    * The constructor for PatchDataRestartManager is private.
    * Consistent with the definition of a Singleton class, only the
    * manager object has access to the constructor for the class.
    *
    * The constructor for PatchDataRestartManager initializes the root
    * data base to a NullDatabase and sets the restart flag to false.
    */
   PatchDataRestartManager() = default;

   /**
    * The destructor for the restart manager is protected, since only the
    * singleton class and subclasses may destroy the manager objects.
    */
   ~PatchDataRestartManager() = default;

   // Unimplemented copy constructor.
   PatchDataRestartManager(
      const PatchDataRestartManager& other) = delete;

   // Unimplemented assignment operator.
   PatchDataRestartManager&
   operator = (
      const PatchDataRestartManager& rhs) = delete;

   /**
    * Deallocate the restart manager instance.  It is not necessary to call
    * this routine at program termination, since it is automatically called
    * by the StartupShutdownManager class.
    */
   static void
   shutdownCallback();

   static PatchDataRestartManager* s_manager_instance;

   /*
    * ComponentSelector holds bits that determine which patch data items need
    * to be written to the restart database.  The bit in position j corresponds
    * to the patch data associated with the j-th index of the patch descriptor
    * object.
    */
   ComponentSelector d_patchdata_restart_table;

   static tbox::StartupShutdownManager::Handler s_shutdown_handler;
};

}
}

#endif
