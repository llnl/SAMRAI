/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in Heat equation ex.
 *
 ************************************************************************/

#ifndef included_ConvDiffXD
#define included_ConvDiffXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/algs/MethodOfLinesIntegrator.h"
#include "SAMRAI/algs/MethodOfLinesPatchStrategy.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Serializable.h"
#include <string>
#define included_String
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include <vector>
#include <memory>

/**
 * The ConvDiff class provides numerical routines for a sample problem
 * which illustrates use of AMR for solution of a system of ODEs.
 * This class is derived from the algs::MethodOfLinesPatchStrategy
 * and provides implementations of the virtual functions declared in that
 * class.  Other member functions are specific to this application.  Most
 * member functions in ConvDiff provide numerical routines that apply to
 * individual patches in an AMR hierarchy.
 *
 * The convection-diffusion equation is  du/dt + div(a*u) = mu * div^2(u)
 * + gamma, where "u" is a scalar-valued function and "a", "mu" and "gamma"
 * are constant vectors.  Time integration of this equation is performed
 * using the algs::MethodOfLinesIntegrator.  The PDE is cast as a set of ODEs
 * (i.e. du/dt = F(u) where F(u) = -div(a*u) + mu*div^2(u) + gamma).
 *
 * The primary numerical quantities are "u" and "F(u)", defined in the
 * code as "primitive_vars" and "function_eval", respectively.  All
 * other variables are temporary quantities used in the numerical routines.
 */

#define NEQU (1)  // depth of u

using namespace SAMRAI;

class ConvDiff:
   public tbox::Serializable,
   public algs::MethodOfLinesPatchStrategy,
   public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor for ConvDiff sets default model parameters to
    * initialize the object. It creates variables that represent
    * the state of the solution, and initializes pertinent private
    * data members.  It also registers the object with the
    * tbox::RestartManager.
    *
    * After setting the default values, the routine calls
    * calls getFromRestart() if this is a restart case.  It next
    * calls getfromInput() to read values from the input database,
    * potentially overriding those in the restart file.
    */
   ConvDiff(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   /**
    * The destructor for ConvDiff.
    */
   ~ConvDiff();

   ///
   ///  The following routines:
   ///
   ///      registerModelVariables(),
   ///      initializeDataOnPatch(),
   ///      computeStableDtOnPatch(),
   ///      singleStep(),
   ///      tagGradientDetectorCells(),
   ///      preprocessRefine(),
   ///      postprocessRefine(),
   ///      preprocessCoarsen(),
   ///      postprocessCoarsen()
   ///
   ///  are concrete implementations of functions declared in the
   ///  algs::MethodOfLinesPatchStrategy abstract base class.
   ///

   /**
    * Register the variables with algs::MethodOfLinesIntegrator.  This
    * registration defines the ways in which data will be manipulated on
    * patches.  Two variable types are available; SOLN and RHS.  For
    * instance, in the solution of du/dt = F(u), u is of SOLN type and
    * F(u) is RHS type.
    *
    * @see algs::MethodOfLinesIntegrator.
    */
   void
   registerModelVariables(
      algs::MethodOfLinesIntegrator* integrator);

   /**
    * Set the data on the patch interior to some initial values,
    * depending on the input parameters and numerical routines.
    * If the "initial_time" flag is false, indicating that the
    * routine is called after a regridding stepa the routine does nothing.
    */
   void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double time,
      const bool initial_time) const;

   /**
    * Compute the stable time increment for a patch using a CFL-based
    * criteria.  Return computed dt.
    */
   double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const double time) const;

   /**
    * Perform a single step of Runge-Kutta routine.  That is, an nth-order
    * RK scheme will perform n sub-iterations at each timestep to integrate
    * over time dt.  The singleStep routine performs one of these
    * sub-iterations.
    */
   void
   singleStep(
      hier::Patch& patch,
      const double dt,
      const double alpha_1,
      const double alpha_2,
      const double beta) const;

   /**
    * Tag cells which need refinement.
    */
   void
   tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_index,
      const bool uses_richardson_extrapolation_too);

   ///
   ///  The following routines:
   ///
   ///      setPhysicalBoundaryConditions()
   ///
   ///  are concrete implementations of functions declared in the
   ///  RefinePatchStrategy abstract base class.
   ///

   /**
    * Set the data in ghost cells corresponding to the physical domain
    * boundary.  Specific boundary conditions are determined by information
    * specified in input file and numerical routines.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector&
      ghost_width_to_fill);

   //@{
   //! @name Required implementations of MethodOfLinesPatchStrategy pure virtuals.

   hier::IntVector
   getRefineOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector::getZero(dim);
   }

   void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   hier::IntVector
   getCoarsenOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector::getZero(dim);
   }

   void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   //@}

   /**
    * Writes state of ConvDiff object to the specified restart database.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * boundary state values from the given database with the
    * given name std::string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary
    * condition applies.
    */
   void
   readDirichletBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   void
   readNeumannBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt
    * visualization tool.
    */
#ifdef HAVE_HDF5
   void
   registerVisItDataWriter(
      std::shared_ptr<appu::VisItDataWriter> viz_writer);
#endif

   /**
    * Prints all class data members, if assertion is thrown.
    */
   void
   printClassData(
      std::ostream& os) const;

private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    *
    * An assertion results if the database pointer is null.
    */
   virtual void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);

   virtual void
   getFromRestart();

   void
   readStateDataEntry(
      std::shared_ptr<tbox::Database> db,
      const std::string& db_name,
      int array_indx,
      std::vector<double>& uval);

   /*
    * Private member function to check correctness of boundary data.
    */
   void
   checkBoundaryData(
      int btype,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const std::vector<int>& scalar_bconds) const;

   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   std::string d_object_name;

   /*
    * Dimension of problem.
    */
   tbox::Dimension d_dim;

   /*
    * std::shared_ptr to the grid geometry object used (Cartesian) to setup
    * initial data and to set physical boundary conditions.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

#ifdef HAVE_HDF5
   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif

   tbox::ResourceAllocator d_allocator;

   /*
    * std::shared_ptrs to variables.  d_primitive_vars - [u]
    *                                   d_function_eval  - [F(u)]
    */
   std::shared_ptr<pdat::CellVariable<double> > d_primitive_vars;
   std::shared_ptr<pdat::CellVariable<double> > d_function_eval;

   /*
    * Convection-diffusion equation constant vectors
    */
   double d_convection_coeff[SAMRAI::MAX_DIM_VAL];
   double d_diffusion_coeff;
   double d_source_coeff;

   /*
    *  Parameters for numerical method:
    *
    *    d_cfl ................ CFL condition for timestepping.
    *
    *    d_tolerance .......... Tolerance used for tagging cells - if
    *                           value[N] > d_tolerance[n] (where N is
    *                           between 0 and NEQU-1) cell is tagged.
    *
    *    d_nghosts ............ number of ghost cells for cell-centered
    *                           variables
    *
    */
   double d_cfl;
   double d_tolerance[NEQU];
   hier::IntVector d_nghosts;
   hier::IntVector d_zero_ghosts;

   /*
    * Indicator for problem type and initial conditions
    */
   std::string d_data_problem;
   int d_data_problem_int;

   /*
    * Input for SPHERE problem
    */
   double d_radius;
   double d_center[SAMRAI::MAX_DIM_VAL];
   double d_val_inside[NEQU];
   double d_val_outside[NEQU];

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET, NEUMANN
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   std::vector<int> d_scalar_bdry_edge_conds;
   std::vector<int> d_scalar_bdry_node_conds;
   std::vector<int> d_scalar_bdry_face_conds; // 3D use only.

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
   std::vector<int> d_node_bdry_edge; // 2D use only.
   std::vector<int> d_edge_bdry_face; // 3D use only.
   std::vector<int> d_node_bdry_face; // 3D use only.

   /*
    * Vectors of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
   std::vector<double> d_bdry_edge_val; // 2D use only.
   std::vector<double> d_bdry_face_val; // 3D use only.

};

#endif
