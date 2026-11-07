/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for patch data test operations.
 *
 ************************************************************************/

#ifndef included_hier_PatchDataTestStrategy
#define included_hier_PatchDataTestStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <memory>


namespace SAMRAI {

class CommTester;

/**
 * Class PatchDataTestStrategy defines an interface for testing specific
 * patch data transfer operations on individual patches when using
 * the CommTester class.  This base class provides two member functions
 * for reading test input information.  These are:
 *
 * readVariableInput(): This function reads in a collection of databases,
 * each of which contains parameters for a single variable.  The names of
 * the sub-databases are arbitrary, but must be distinct.  Each variable
 * sub-database has the following input keys:
 *
 *
 *
 *    - \b  name         variable name std::string (required)
 *    - \b  depth        optional variable depth (default = 1)
 *    - \b  src_ghosts   optional comm source ghost width (default = 0,0,0)
 *    - \b  dst_ghosts   optional comm dest ghost width (default = 0,0,0)
 *    - \b  coarsen_operator   opt. coarsen op name (default = "NO_COARSEN")
 *    - \b  refine_operator    opt. refine op name (default = "NO_REFINE")
 *
 *
 *
 * readRefinementInput(): This function reads in a collection of box
 * arrays, each of which describes the region to refine on each level.
 * The key names of the box arrays are arbitrary, but must be distinct.
 * For example,
 *
 *
 *
 *    - \b  level0_boxes   boxes to refine on level zero.
 *    - \b  level1_boxes   boxes to refine on level one.
 *    - \b  ...            etc...
 *
 *
 *
 *
 * The pure virtual functions in this class that must be provided by
 * concrete test subclass:
 * \begin{enumerate}
 *    - [registerVariables(...)] register variables with CommTester.
 *    - [initializeDataOnPatch(...)] set patch data on new patch.
 *    - [verifyResults(...)] check results of communication operations.
 * \end{enumerate}
 *
 * The following virtual functions are given default non-operations in this
 * class so that concrete test subclass can either implement them to test
 * specific functionality or simply ignore.  They are pure virtual in the
 * coarsen and refine patch strategy classes:
 * \begin{enumerate}
 *    - [setPhysicalBoundaryConditions(...)]
 *    - [preprocessRefine(...)]
 *    - [postprocessRefine(...)]
 *    - [preprocessCoarsen(...)]
 *    - [postprocessCoarsen(...)]
 * \end{enumerate}
 */

class PatchDataTestStrategy
{
public:
   /**
    * The constructor initializes variable data arrays to zero length.
    */
   PatchDataTestStrategy(
      const tbox::Dimension& dim);

   /**
    * Virtual destructor for PatchDataTestStrategy.
    */
   virtual ~PatchDataTestStrategy();

   /**
    * Grid geometry access operations.
    */
   void setGridGeometry(
      std::shared_ptr<geom::CartesianGridGeometry> grid_geom)
   {
      TBOX_ASSERT(grid_geom);
      d_grid_geometry = grid_geom;
   }

   ///
   std::shared_ptr<geom::CartesianGridGeometry> getGridGeometry() const
   {
      return d_grid_geometry;
   }

   /**
    * Utility functions for managing patch data context.
    */
   void setDataContext(
      std::shared_ptr<hier::VariableContext> context)
   {
      TBOX_ASSERT(context);
      d_data_context = context;
   }

   ///
   std::shared_ptr<hier::VariableContext> getDataContext() const
   {
      return d_data_context;
   }

   void clearDataContext()
   {
      d_data_context.reset();
   }

   /**
    * Read variable parameters from input database.
    */
   void
   readVariableInput(
      std::shared_ptr<tbox::Database> db);

   /**
    * Virtual functions in interface to user-supplied boundary conditions,
    * coarsen and refine operations.
    */
   virtual void
   setPhysicalBoundaryConditions(
      const hier::Patch& patch,
      const double time,
      const hier::IntVector& gcw) const;

   ///
   virtual void
   preprocessRefine(
      const hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;

   ///
   virtual void
   postprocessRefine(
      const hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;

   ///
   virtual void
   preprocessCoarsen(
      const hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;

   ///
   virtual void
   postprocessCoarsen(
      const hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;

   /**
    * This function is called from the CommTester constructor.  Its
    * purpose is to register variables used in the patch data test
    * and appropriate communication parameters (ghost cell widths,
    * coarsen/refine operations) with the CommTester object, which
    * manages the variable storage.
    */
   virtual void
   registerVariables(
      CommTester* commtest) = 0;

   /**
    * Virtual function for setting data on new patch in hierarchy.
    *
    * @param src_or_dst Flag set to 's' for source or 'd' for destination
    *        to indicate variables to set data for.
    */
   virtual void
   initializeDataOnPatch(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number,
      char src_or_dst) = 0;
   /**
    * Virtual function for checking results of communication operations.
    *
    * @returns Whether test passed.
    */
   virtual bool
   verifyResults(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number) = 0;

   virtual void
   setDataIds(std::list<int>& data_ids)
   {
      NULL_USE(data_ids);
   }

   virtual bool
   verifyCompositeBoundaryData(
      const hier::Patch& patch,
      const std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int data_id,
      int level_number,
      const std::vector<std::shared_ptr<hier::PatchData> >& bdry_data)
   {
      NULL_USE(patch);
      NULL_USE(hierarchy);
      NULL_USE(data_id);
      NULL_USE(level_number);
      NULL_USE(bdry_data);
      return true;
   }

protected:
   const tbox::Dimension d_dim;
   /*
    * Vectors of information read from input file describing test variables
    */
   std::vector<std::string> d_variable_src_name;
   std::vector<std::string> d_variable_dst_name;
   std::vector<int> d_variable_depth;
   std::vector<hier::IntVector> d_variable_src_ghosts;
   std::vector<hier::IntVector> d_variable_dst_ghosts;
   std::vector<std::string> d_variable_coarsen_op;
   std::vector<std::string> d_variable_refine_op;

private:
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   std::shared_ptr<hier::VariableContext> d_data_context;

};

}
#endif
