/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Reference counting class for Array
 *
 ************************************************************************/

#include "SAMRAI/tbox/ReferenceCounter.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <cstdlib>
#include <cstdio>

namespace SAMRAI {
namespace tbox {

ReferenceCounter * ReferenceCounter::s_free_list = 0;
bool ReferenceCounter::s_is_finalized = false;

StartupShutdownManager::Handler
ReferenceCounter::s_handler(
   0,
   0,
   0,
   ReferenceCounter::finalizeCallback,
   StartupShutdownManager::priorityReferenceCounter);

ReferenceCounter::ReferenceCounter()
{
   d_references = 1;
   d_next = 0;
}

ReferenceCounter::~ReferenceCounter()
{
   if ((d_next) && (--d_next->d_references == 0)) {
      delete d_next;
   }
}

void *
ReferenceCounter::operator new (
   size_t bytes)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
   /* Since low level class; tbox Utilities may not function here */
   assert(!ReferenceCounter::isFinalized());
#endif

   if (s_free_list) {
      ReferenceCounter* node = s_free_list;
      s_free_list = s_free_list->d_next;
      return node;
   } else {
      return ::operator new (
                bytes);
   }
}

void
ReferenceCounter::operator delete (
   void* what)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
   /* Since low level class; tbox Utilities may not function here */
   assert(!ReferenceCounter::isFinalized());
#endif

   ReferenceCounter* node = (ReferenceCounter *)what;
   node->d_next = s_free_list;
   s_free_list = node;
}

void
ReferenceCounter::finalizeCallback()
{
   while (s_free_list) {
      void * byebye = s_free_list;
      s_free_list = s_free_list->d_next
      ;
      ::operator delete (
         byebye);
   }

   s_is_finalized = true;
}

}
}
