/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class to test usage of boundary utilities
 *
 ************************************************************************/

#ifndef included_BoundaryDataTester
#define included_BoundaryDataTester

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <vector>
#include <memory>

using namespace SAMRAI;

class BoundaryDataTester:
   public xfer::RefinePatchStrategy,
   public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor reads variable data from input database.
    */
   BoundaryDataTester(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   /**
    * Virtual destructor for BoundaryDataTester.
    */
   virtual ~BoundaryDataTester();

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class RefinePatchStrategy.  It sets the boundary
    * conditions for the variables.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);

   /**
    * The next three functions are dummy implementations of the pure
    * virtual functions declared in the RefinePatchStrategy base class.
    * They are not needed for this example since we only have one level
    * in the hierarchy.
    */
   hier::IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector(dim, 0);
   }

   void preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   void postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void
   readDirichletBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads NEUMANN
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void
   readNeumannBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * Set data on patch interiors on given level in hierarchy.
    */
   void
   initializeDataOnPatchInteriors(
      std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

   /**
    * Run boundary tests for given level in hierarchy and return integer
    * number of test failures.
    */
   int
   runBoundaryTest(
      std::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

   /**
    * Print all class data members to given output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension& getDim() const
   {
      return d_dim;
   }

private:
   /*
    * The object name is used for error/warning reporting.
    */
   std::string d_object_name;

   const tbox::Dimension d_dim;

   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   /*
    * Vectors of information read from input file describing test variables
    */
   std::vector<std::string> d_variable_name;
   std::vector<int> d_variable_depth;
   std::vector<hier::IntVector> d_variable_num_ghosts;
   std::vector<std::vector<double> > d_variable_interior_values;

   /*
    * Items used to manage variables and data in test program.
    */
   std::vector<std::shared_ptr<hier::Variable> > d_variables;
   std::shared_ptr<hier::VariableContext> d_variable_context;
   hier::ComponentSelector d_patch_data_components;

   /*
    * Vectors of information read from input file for boundary conditions
    */
   std::vector<int> d_master_bdry_edge_conds;
   std::vector<int> d_scalar_bdry_edge_conds;
   std::vector<int> d_vector_bdry_edge_conds;

   std::vector<int> d_master_bdry_node_conds;
   std::vector<int> d_scalar_bdry_node_conds;
   std::vector<int> d_vector_bdry_node_conds;

   std::vector<int> d_master_bdry_face_conds; // Used only in 3D
   std::vector<int> d_scalar_bdry_face_conds; // Used only in 3D
   std::vector<int> d_vector_bdry_face_conds; // Used only in 3D

   std::vector<int> d_node_bdry_edge; // Used only in 2D
   std::vector<int> d_edge_bdry_face; // Used only in 3D
   std::vector<int> d_node_bdry_face; // Used only in 3D

   std::vector<std::vector<double> > d_variable_bc_values;

   int d_fail_count;

   /*
    * Private functions to perform tasks for boundary testing.
    */
   void
   readVariableInputAndMakeVariables(
      std::shared_ptr<tbox::Database> db);
   void
   readBoundaryDataInput(
      std::shared_ptr<tbox::Database> db);
   void
   readBoundaryDataStateEntry(
      std::shared_ptr<tbox::Database> db,
      std::string& db_name,
      int bdry_location_index);
   void
   setBoundaryDataDefaults();
   void
   postprocessBoundaryInput();
   void
   checkBoundaryData(
      int btype,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_check);

};

#endif
