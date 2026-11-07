/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for cell-centered patch data
 *
 ************************************************************************/

#ifndef included_pdat_CellDataTest
#define included_pdat_CellDataTest

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"

#include "PatchDataTestStrategy.h"

#include <string>
#include <memory>



namespace SAMRAI {

class CommTester;

using CELL_KERNEL_TYPE = double; // float | double | dcomplex

/**
 * Class CellDataTest provides routines to test communication operations
 * for cell-centered patch data on an AMR patch hierarchy.
 *
 * Required input keys and data types:
 *
 *
 *
 *
 *   Double values that define linear function initial data to test refine
 *   operations (Ax + By + Cz + D = f(x,y,z), where f(x,y,z) is the value
 *   assigned to each array value at initialization and against which
 *   linear interpolation is tested:
 *
 *    Acoef, Dcoef always required.
 *    If (dim > 1), Bcoef is needed.
 *    If (dim > 2), Ccoef is needed.
 *
 *
 *
 *
 *
 * See PatchDataTestStrategy header file comments for variable and
 * refinement input data description.
 */

class CellDataTest:public PatchDataTestStrategy
{
public:
   /**
    * The constructor initializes variable data arrays to zero length.
    */
   CellDataTest(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> main_input_db,
      bool do_refine,
      bool do_coarsen,
      const std::string& refine_option);

   /**
    * Virtual destructor for CellDataTest.
    */
   ~CellDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   void
   setPhysicalBoundaryConditions(
      const hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw_to_fill) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   void
   registerVariables(
      CommTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void
   initializeDataOnPatch(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number,
      char src_or_dst);

   /**
    * Function for checking results of communication operations.
    */
   bool
   verifyResults(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

   void
   setDataIds(std::list<int>& data_ids);

   bool
   verifyCompositeBoundaryData(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int data_id,
      int level_number,
      const std::vector<std::shared_ptr<hier::PatchData> >& bdry_data);

#ifdef SAMRAI_HAVE_CONDUIT
   void addFields(
      conduit::Node& node,
      int domain_id,
      const std::shared_ptr<hier::Patch>& patch);
#endif

private:
   /**
    * Function for reading test data from input file.
    */
   void
   readTestInput(
      std::shared_ptr<tbox::Database> db);

   /**
    * Set linear function data for testing interpolation
    */
   void
   setLinearData(
      std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   /**
    * Set conservative linear function data for testing coarsening
    */
   void
   setConservativeData(
      std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
      const hier::Box& box,
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number) const;

   /**
    * Set periodic linear function data for testing interpolation in
    * periodic domains.
    */
   void
   setPeriodicData(
      std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   void
   checkPatchInteriorData(
      const std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> >& data,
      const hier::Box& interior,
      const hier::Patch& patch) const;

   const tbox::Dimension d_dim;

   /*
    * Object string identifier for error reporting
    */
   std::string d_object_name;

   /*
    * Data members specific to this cell data test.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_cart_grid_geometry;

   /*
    * Data members specific to this cell data test.
    */
   double d_Acoef;
   double d_Bcoef;
   double d_Ccoef;
   double d_Dcoef;

   bool d_do_refine;
   bool d_do_coarsen;
   std::string d_refine_option;
   int d_finest_level_number;

   std::vector<std::shared_ptr<hier::Variable> > d_variables;

};
} // Namespace SAMRAI
#endif
