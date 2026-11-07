/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An restart manager singleton class
 *
 ************************************************************************/

#include "SAMRAI/hier/PatchDataRestartManager.h"

namespace SAMRAI {
namespace hier {

PatchDataRestartManager * PatchDataRestartManager::s_manager_instance = 0;

tbox::StartupShutdownManager::Handler
PatchDataRestartManager::s_shutdown_handler(
   0,
   0,
   PatchDataRestartManager::shutdownCallback,
   0,
   tbox::StartupShutdownManager::priorityRestartManager);

/*
 *************************************************************************
 *
 * Basic singleton classes to create, set, and destroy the manager
 * instance.
 *
 *************************************************************************
 */

PatchDataRestartManager *
PatchDataRestartManager::getManager()
{
   if (!s_manager_instance) {
      s_manager_instance = new PatchDataRestartManager;
   }
   return s_manager_instance;
}

void
PatchDataRestartManager::shutdownCallback()
{
   if (s_manager_instance) {
      delete s_manager_instance;
      s_manager_instance = 0;
   }
}

}
}
