/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A wrapper for Umpire allocators.
 *
 ************************************************************************/

#ifndef included_tbox_ResourceAllocator
#define included_tbox_ResourceAllocator

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_UMPIRE
#include "umpire/Allocator.hpp"
#endif

namespace SAMRAI {
namespace tbox {

/**
 * ResourceAllocator a type alias for umpire::Allocator that enables API
 * consistency when SAMRAI is built with or without the Resource library.
 * When Umpire is available ResourceAllocator is an alias for umpire::Allocator,
 * so calling codes can pass in an umpire::Allocator anywhere that
 * ResourceAllocator is required in the SAMRAI API.  If Umpire is not
 * available, ResourceAllocator is an empty struct.
 */

#ifdef HAVE_UMPIRE

using ResourceAllocator = umpire::Allocator;

#else

struct ResourceAllocator
{
};

#endif

}
}

#endif
