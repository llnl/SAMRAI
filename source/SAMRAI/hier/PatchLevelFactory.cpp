/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchLevelFactory.h"


namespace SAMRAI {
namespace hier {

std::shared_ptr<PatchLevel>
PatchLevelFactory::allocate(
   const BoxLevel& box_level,
   const std::shared_ptr<BaseGridGeometry>& grid_geometry,
   const std::shared_ptr<PatchDescriptor>& descriptor,
   const std::shared_ptr<PatchFactory>& factory) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box_level, *grid_geometry);
   std::shared_ptr<PatchLevel> pl(
      std::make_shared<PatchLevel>(
         box_level,
         grid_geometry,
         descriptor,
         factory));
   return pl;
}

std::shared_ptr<PatchLevel>
PatchLevelFactory::allocate(
   const std::shared_ptr<BoxLevel> box_level,
   const std::shared_ptr<BaseGridGeometry>& grid_geometry,
   const std::shared_ptr<PatchDescriptor>& descriptor,
   const std::shared_ptr<PatchFactory>& factory) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*box_level, *grid_geometry);
   std::shared_ptr<PatchLevel> pl(
      std::make_shared<PatchLevel>(
         box_level,
         grid_geometry,
         descriptor,
         factory));
   return pl;
}

std::shared_ptr<PatchLevel>
PatchLevelFactory::allocate(
   const std::shared_ptr<tbox::Database>& database,
   const std::shared_ptr<BaseGridGeometry>& grid_geometry,
   const std::shared_ptr<PatchDescriptor>& descriptor,
   const std::shared_ptr<PatchFactory>& factory,
   const bool defer_boundary_box_creation) const
{
   std::shared_ptr<PatchLevel> pl(
      std::make_shared<PatchLevel>(
         database,
         grid_geometry,
         descriptor,
         factory,
         defer_boundary_box_creation));
   return pl;
}

}
}
