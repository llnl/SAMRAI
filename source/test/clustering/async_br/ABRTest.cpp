/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   ABRTest class implementation
 *
 ************************************************************************/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ABRTest.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"

#include <iomanip>

using namespace SAMRAI;

ABRTest::ABRTest(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
   std::shared_ptr<tbox::Database> database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(patch_hierarchy),
   d_sine_wall(object_name + ":tagger",
               d_dim,
               database->getDatabaseWithDefault("sine_tagger", std::shared_ptr<tbox::Database>()))
{
   d_sine_wall.resetHierarchyConfiguration(d_hierarchy, 0, 0);
}

ABRTest::~ABRTest()
{
}

mesh::StandardTagAndInitStrategy *ABRTest::getStandardTagAndInitObject()
{
   return &d_sine_wall;
}

/*
 * Deallocate patch data allocated by this class.
 */
void ABRTest::computeHierarchyData(
   hier::PatchHierarchy& hierarchy,
   double time)
{
   d_sine_wall.computeHierarchyData(hierarchy, time);
}

/*
 * Deallocate patch data allocated by this class.
 */
void ABRTest::deallocatePatchData(
   hier::PatchHierarchy& hierarchy)
{
   d_sine_wall.deallocatePatchData(hierarchy);
}

/*
 * Deallocate patch data allocated by this class.
 */
void ABRTest::deallocatePatchData(
   hier::PatchLevel& level)
{
   d_sine_wall.deallocatePatchData(level);
}

#ifdef HAVE_HDF5
int ABRTest::registerVariablesWithPlotter(
   std::shared_ptr<appu::VisItDataWriter> writer)
{
   if (writer)
      d_sine_wall.registerVariablesWithPlotter(*writer);
   return 0;
}
#endif

bool ABRTest::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   if (variable_name == "Patch level number") {
      double pln = patch.getPatchLevelNumber();
      size_t i, size = region.size();
      for (i = 0; i < size; ++i) buffer[i] = pln;
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        << "ABRTest::packDerivedPatchDataIntoDoubleBuffer");
   }

   return false;
}
