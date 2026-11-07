/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for edge-centered patch data
 *
 ************************************************************************/

#ifndef included_SideMultiblockTest
#define included_SideMultiblockTest

#include "SAMRAI/SAMRAI_config.h"

#include "PatchMultiblockTestStrategy.h"

using namespace SAMRAI;

/**
 * Class SideMultiblockTest provides routines to test communication operations
 * for edge-centered patch data on an AMR patch hierarchy.
 *
 * See PatchMultiblockTestStrategy header file comments for variable and
 * refinement input data description.
 */

class SideMultiblockTest:public PatchMultiblockTestStrategy
{
public:
   /**
    * The constructor initializes variable data arrays to zero length.
    */
   SideMultiblockTest(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> main_input_db,
      const std::string& refine_option);

   /**
    * Virtual destructor for SideMultiblockTest.
    */
   virtual ~SideMultiblockTest();

   /**
    * User-supplied boundary conditions.  Note that we do not implement
    * user-defined refine operations.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw_to_fill) const;

   void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      std::shared_ptr<const hier::Connector> dst_to_encon,
      const hier::Box& fill_box,
      const hier::BoundaryBox& bbox,
      const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry);

   /**
    * This function is called from the MultiblockTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * refine operations) with the MultiblockTester object, which
    * manages the variable storage.
    */
   void
   registerVariables(
      MultiblockTester* commtest);

   /**
    * Function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void
   initializeDataOnPatch(
      hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const hier::BlockId& block_id,
      char src_or_dst);

   /**
    * Function for tagging cells on each patch to refine.
    */
   void
   tagCellsToRefine(
      hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number,
      int tag_index);

   /**
    * Function for checking results of communication operations.
    */
   bool
   verifyResults(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      const int level_number,
      const hier::BlockId& block_id);

   ///
   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const std::shared_ptr<hier::VariableContext>& context,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;

private:
   /**
    * Function for reading test data from input file.
    */
   void
   readTestInput(
      std::shared_ptr<tbox::Database> db);

   /*
    * Object string identifier for error reporting
    */
   std::string d_object_name;

   const tbox::Dimension d_dim;

   std::string d_refine_option;
   int d_finest_level_number;

   std::vector<std::shared_ptr<hier::Variable> > d_variables;

};

#endif
