/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI class to manage package startup and shutdown
 *
 ************************************************************************/

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <new>

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Set static member to set maximum number of patch data entries,
 * statistics, and timers supported by the code.
 * This value is used to set the sizes of certain array containers
 * used in the SAMRAI library.  It is set here so that it may be
 * resized if necessary from a single access point (i.e., via the
 * SAMRAIManager) if the default size (set here) are insufficient.
 *
 * This value can be changed by calling the static functions:
 * SAMRAIManager::setMaxNumberPatchDataEntries(),
 *
 *************************************************************************
 */

bool SAMRAIManager::s_initialized = false;
bool SAMRAIManager::s_started = false;
int SAMRAIManager::s_max_patch_data_entries = 256;

/*
 *************************************************************************
 *
 * Initialize the SAMRAI package.  This routine performs the following
 * tasks:
 *
 * (1) Initialize the SAMRAI MPI package
 * (2) Initialize the parallel I/O routines
 * (3) Set up IEEE assertion handlers
 * (4) Set new handler so that an error message is printed if new fails.
 *
 *************************************************************************
 */

static void badnew()
{
   TBOX_ERROR("operator `new' failed -- program abort!" << std::endl);
}

void
SAMRAIManager::initialize(
   bool initialize_IEEE_assertion_handlers)
{
   TBOX_ASSERT(!s_initialized);

   PIO::initialize();

   if (initialize_IEEE_assertion_handlers) {
      IEEE::setupFloatingPointExceptionHandlers();
   }

#ifndef LACKS_PROPER_MEMORY_HANDLER
   std::set_new_handler(badnew);
#endif

   StartupShutdownManager::initialize();

   s_initialized = true;
}

void
SAMRAIManager::setMaxNumberPatchDataEntries(
   int maxnum)
{
   s_max_patch_data_entries = MathUtilities<int>::Max(maxnum,
         s_max_patch_data_entries);
}

}
}
