/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   TreeLoadBalancer test.
 *
 ************************************************************************/
#include "DerivedVisOwnerData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

using namespace SAMRAI;

DerivedVisOwnerData::DerivedVisOwnerData()
{
}

DerivedVisOwnerData::~DerivedVisOwnerData()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
#ifdef HAVE_HDF5
int DerivedVisOwnerData::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   writer.registerDerivedPlotQuantity("Owner", "SCALAR", this);
   return 0;
}
#endif

bool DerivedVisOwnerData::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(patch);
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   if (variable_name == "Owner") {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      double owner = mpi.getRank();
      size_t i, size = region.size();
      for (i = 0; i < size; ++i) buffer[i] = owner;
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        <<
         "DerivedVisOwnerData::packDerivedPatchDataIntoDoubleBuffer");
   }

   return true;
}
