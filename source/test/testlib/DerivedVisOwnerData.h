/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   TreeLoadBalancer test.
 *
 ************************************************************************/
#ifndef included_DerivedVisOwnerData
#define included_DerivedVisOwnerData

#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"

#include "SAMRAI/tbox/Database.h"

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"


using namespace SAMRAI;

/*!
 * @brief Write owner rank using VisDerivedDataStrategy.
 */
class DerivedVisOwnerData:
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   DerivedVisOwnerData();

   ~DerivedVisOwnerData();

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer);
#endif

   //@{ @name SAMRAI::appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time) const;

   //@}

private:
};

#endif  // included_DerivedVisOwnerData
