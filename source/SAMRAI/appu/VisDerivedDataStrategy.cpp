/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for writing user-defined data to VisIt
 *
 ************************************************************************/
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace appu {

VisDerivedDataStrategy::VisDerivedDataStrategy()
{
}

VisDerivedDataStrategy::~VisDerivedDataStrategy()
{
}

bool
VisDerivedDataStrategy::packMixedDerivedDataIntoDoubleBuffer(
   double* buffer,
   std::vector<double>& mixbuffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index) const
{
   NULL_USE(buffer);
   NULL_USE(mixbuffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(variable_name);
   NULL_USE(depth_index);
   TBOX_ERROR("VisDerivedDataStrategy::"
      << "packMixedDerivedDataIntoDoubleBuffer()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes nothing) is executed" << std::endl);
   return 0;
}

}
}
