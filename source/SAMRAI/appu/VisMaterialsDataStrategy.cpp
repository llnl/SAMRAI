/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for writing material related data to VisIt
 *                dump file
 *
 ************************************************************************/
#include "SAMRAI/appu/VisMaterialsDataStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace appu {

VisMaterialsDataStrategy::VisMaterialsDataStrategy()
{
}

VisMaterialsDataStrategy::~VisMaterialsDataStrategy()
{
}

int
VisMaterialsDataStrategy::packMaterialFractionsIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& material_name) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(material_name);
   TBOX_ERROR("VisMaterialsDataStrategy::"
      << "packMaterialFractionsIntoDoubleBuffer()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes nothing) is executed" << std::endl);
   return 0;
}

int
VisMaterialsDataStrategy::packMaterialFractionsIntoSparseBuffers(
   int* mat_list,
   std::vector<int>& mix_zones,
   std::vector<int>& mix_mat,
   std::vector<double>& vol_fracs,
   std::vector<int>& next_mat,
   const hier::Patch& patch,
   const hier::Box& region) const
{
   NULL_USE(mat_list);
   NULL_USE(mix_zones);
   NULL_USE(mix_mat);
   NULL_USE(vol_fracs);
   NULL_USE(next_mat);
   NULL_USE(patch);
   NULL_USE(region);
   TBOX_ERROR("VisMaterialsDataStrategy::"
      << "packSparseMaterialFractionsIntoDoubleBuffer()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes nothing) is executed" << std::endl);
   return 0;
}

int
VisMaterialsDataStrategy::packSpeciesFractionsIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& material_name,
   const std::string& species_name) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(material_name);
   NULL_USE(species_name);
   TBOX_ERROR("VisMaterialsDataStrategy::"
      << "packSpeciesFractionsIntoDoubleBuffer()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes nothing) is executed" << std::endl);
   return 0;
}

}
}
