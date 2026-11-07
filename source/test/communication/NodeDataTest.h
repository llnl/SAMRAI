/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data
 *
 ************************************************************************/

#ifndef included_pdat_NodeDataTest
#define included_pdat_NodeDataTest

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/Patch.h"
#include "PatchDataTestStrategy.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "SAMRAI/hier/Variable.h"

#include <memory>



namespace SAMRAI {

class CommTester;

using NODE_KERNEL_TYPE = double; // use float | double | dcomplex
/**
 * Class NodeDataTest provides routines to test communication operations
 * for node-centered patch data on an AMR patch hierarchy.
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


class NodeDataTest:public PatchDataTestStrategy
{
public:
   /**
    * The constructor initializes variable data arrays to zero length.
    */
   NodeDataTest(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> main_input_db,
      bool do_refine,
      bool do_coarsen,
      const std::string& refine_option);

   /**
    * Virtual destructor for NodeDataTest.
    */
   ~NodeDataTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined coarsen and refine operations.
    */
   virtual void
   setPhysicalBoundaryConditions(
      const hier::Patch& patch,
      const double time,
      const hier::IntVector&) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost node widths,
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
      std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   /**
    * Set periodic linear function data for testing interpolation in
    * periodic domains.
    */
   void
   setPeriodicData(
      std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > data,
      const hier::Box& box,
      const hier::Patch& patch) const;

   void
   checkPatchInteriorData(
      const std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> >& data,
      const hier::Box& interior,
      const hier::Patch& patch) const;

   const tbox::Dimension d_dim;

   /*
    * Object std::string identifier for error reporting
    */
   std::string d_object_name;

   /*
    * Data members specific to this node data test.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_cart_grid_geometry;

   /*
    * Data members specific to this node data test.
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

}
#endif
