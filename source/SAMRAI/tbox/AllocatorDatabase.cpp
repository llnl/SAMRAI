/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton database for managing Umpire allocators 
 *
 ************************************************************************/

#include "SAMRAI/tbox/AllocatorDatabase.h"

#ifdef HAVE_UMPIRE
#include "umpire/strategy/QuickPool.hpp"
#include "umpire/ResourceManager.hpp"

#include "umpire/resource/MemoryResourceTypes.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"
#endif

#include <iostream>

namespace SAMRAI {
namespace tbox {

AllocatorDatabase * AllocatorDatabase::s_allocator_database_instance(0);

StartupShutdownManager::Handler
AllocatorDatabase::s_startup_handler(
    0,
    AllocatorDatabase::startupCallback,
    0,
    0,
    tbox::StartupShutdownManager::priorityArenaManager);

void
AllocatorDatabase::startupCallback()
{
  AllocatorDatabase::getDatabase()->initialize();
}

void
AllocatorDatabase::shutdownCallback()
{
   if (s_allocator_database_instance) {
      delete s_allocator_database_instance;
   }
   s_allocator_database_instance = 0;
}

AllocatorDatabase *
AllocatorDatabase::getDatabase()
{
   if (!s_allocator_database_instance) {
      s_allocator_database_instance = new AllocatorDatabase();
   }
   return s_allocator_database_instance;
}

AllocatorDatabase::~AllocatorDatabase()
{
}

void
AllocatorDatabase::initialize()
{
#if defined(HAVE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  if (!rm.isAllocator("samrai::data_allocator")) {
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    // Internal pool for allocations
#if 0
    auto allocator = rm.makeAllocator<umpire::strategy::AllocationAdvisor>(
        "internal::samrai::um_allocation_advisor",
        rm.getAllocator(umpire::resource::Unified),
        // Set preferred location to GPU
        "SET_PREFERRED_LOCATION");
#endif

#if defined(USE_DEVICE_ALLOCATOR)
    auto allocator = rm.getAllocator(umpire::resource::Device);
#else
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#endif

#else 
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif

    rm.makeAllocator<umpire::strategy::QuickPool>("samrai::data_allocator", allocator);

  }

  if (!rm.isAllocator("samrai::tag_allocator")) {
    rm.makeAllocator<umpire::strategy::QuickPool>("samrai::tag_allocator",
        rm.getAllocator(umpire::resource::Host));
  }

  if (!rm.isAllocator("samrai::stream_allocator")) {
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
#if defined(USE_DEVICE_ALLOCATOR)
    auto allocator = rm.getAllocator(umpire::resource::Device);
#else
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#endif
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif

    rm.makeAllocator<umpire::strategy::QuickPool>("samrai::stream_allocator", allocator);
  }

  if (!rm.isAllocator("samrai::fuser_allocator")) {
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif

    rm.makeAllocator<umpire::strategy::QuickPool>("samrai::fuser_allocator", allocator);
  }

  if (!rm.isAllocator("samrai::temporary_data_allocator")) {
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
#if defined(USE_DEVICE_ALLOCATOR)
    auto allocator = rm.getAllocator(umpire::resource::Device);
#else
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#endif
#else
    auto allocator = rm.getAllocator(umpire::resource::Host);
#endif
    rm.makeAllocator<umpire::strategy::QuickPool>("samrai::temporary_data_allocator", allocator);
  }

#endif
}

ResourceAllocator
AllocatorDatabase::getDevicePool()
{
#if defined(HAVE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("samrai::temporary_data_allocator");
#else
  return ResourceAllocator();
#endif
}


#if defined(HAVE_UMPIRE)
umpire::TypedAllocator<char>
AllocatorDatabase::getStreamAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return umpire::TypedAllocator<char>(rm.getAllocator("samrai::stream_allocator"));
}

umpire::TypedAllocator<char>
AllocatorDatabase::getKernelFuserAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return umpire::TypedAllocator<char>(rm.getAllocator("samrai::fuser_allocator"));
}

umpire::TypedAllocator<char>
AllocatorDatabase::getInternalHostAllocator()
{
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return umpire::TypedAllocator<char>(rm.getAllocator(umpire::resource::Host));
}
#endif

ResourceAllocator
AllocatorDatabase::getTagAllocator()
{
#if defined(HAVE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("samrai::tag_allocator");
#else
  return ResourceAllocator();
#endif
}

ResourceAllocator
AllocatorDatabase::getDefaultAllocator()
{
#if defined(HAVE_UMPIRE)
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

  return rm.getAllocator("samrai::data_allocator");
#else
  return ResourceAllocator();
#endif
}


}
}

