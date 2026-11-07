/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   DLBGTest class implementation
 *
 ************************************************************************/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "DLBGTest.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"

#include <iomanip>

using namespace SAMRAI;


DLBGTest::DLBGTest(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
   std::shared_ptr<tbox::Database> database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(patch_hierarchy),
   d_sine_wall(object_name + ":tagger",
               d_dim,
               database->getDatabaseWithDefault("sine_tagger", std::shared_ptr<tbox::Database>())),
   d_time(0.5)
{
   d_sine_wall.resetHierarchyConfiguration(patch_hierarchy, 0, 0);
}

DLBGTest::~DLBGTest()
{
}

mesh::StandardTagAndInitStrategy *DLBGTest::getStandardTagAndInitObject()
{
   return &d_sine_wall;
}

/*
 * Deallocate patch data allocated by this class.
 */
void DLBGTest::computeHierarchyData(
   hier::PatchHierarchy& hierarchy,
   double time)
{
   d_sine_wall.computeHierarchyData(hierarchy, time);
}

/*
 * Deallocate patch data allocated by this class.
 */
void DLBGTest::deallocatePatchData(
   hier::PatchHierarchy& hierarchy)
{
   d_sine_wall.deallocatePatchData(hierarchy);
}

/*
 * Deallocate patch data allocated by this class.
 */
void DLBGTest::deallocatePatchData(
   hier::PatchLevel& level)
{
   d_sine_wall.deallocatePatchData(level);
}

#ifdef HAVE_HDF5
int DLBGTest::registerVariablesWithPlotter(
   std::shared_ptr<appu::VisItDataWriter> writer)
{
   if (writer) {
      d_sine_wall.registerVariablesWithPlotter(*writer);
   }
   return 0;
}
#endif

bool DLBGTest::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(variable_name);
   NULL_USE(depth_id);
   NULL_USE(simulation_time);
   TBOX_ERROR("Should not be here.  This object didn't register any derived plot variables.");
   return true;
}
