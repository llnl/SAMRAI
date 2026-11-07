/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Reference counting class for Array
 *
 ************************************************************************/

#ifndef included_tbox_ReferenceCounter
#define included_tbox_ReferenceCounter

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <cstddef>

namespace SAMRAI {
namespace tbox {

/**
 * Class ReferenceCounter manages the shared reference counter and
 * arena resources used by Array.  It uses a local
 * free pool of objects to speed memory allocation and deallocation.  The
 * locally cached free pool can be freed by calling freeCachedCopies().
 *
 * {\b Do not subclass!}  Changing the size of a ReferenceCounter
 * object will cause my simple memory allocation mechanism to break in
 * horrible and disgusting ways.
 *
 * @see Array
 */

class ReferenceCounter
{
public:
   /**
    * Create a ReferenceCounter.
    * The number of references is set to one.
    */
   ReferenceCounter();

   /**
    * Destructor for ReferenceCounter.  The destructor releases
    * the managed memory arena if its count has gone to zero.
    */
   ~ReferenceCounter();

   /**
    * Decrement the number of references.  True is returned if the
    * reference count has gone to zero; false otherwise.
    */
   bool
   deleteReference()
   {
      return --d_references == 0;
   }

   /**
    * Increment the number of references.
    */
   void
   addReference()
   {
      ++d_references;
   }

   /**
    * Class-specific operator new.  Data is allocated off of an
    * internal free list to speed memory allocation.
    */
   void *
   operator new (
      size_t bytes);

   /**
    * Class-specific operator delete.  Freed data is returned to
    * an internal free list for re-use by operator new.
    */
   void
   operator delete (
      void* what);

   /**
    * Returns true if ReferenceCounter class has been be stopped.
    * This method is used only for a debugging check assert in
    * the pointer class and should not normally be used.
    */
   static bool
   isFinalized()
   {
      return s_is_finalized;
   }

private:
   ReferenceCounter(
      const ReferenceCounter&);                 // not implemented
   ReferenceCounter&
   operator = (
      const ReferenceCounter&);                 // not implemented

   /**
    * Release the memory for all currently cached ReferenceCounter
    * copies.
    */
   static void
   finalizeCallback();

   int d_references;

   // Next element on free list
   ReferenceCounter* d_next;

   // Free list of ReferenceCounter objects
   static ReferenceCounter* s_free_list;

   static StartupShutdownManager::Handler s_handler;

   static bool s_is_finalized;
};

}
}

#endif
