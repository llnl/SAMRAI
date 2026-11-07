/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in Euler equation ex.
 *
 ************************************************************************/

#ifndef included_EulerXD
#define included_EulerXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MessageStream.h"

#include <string>
#include <vector>
#include <memory>


/**
 * The Euler class provides routines for a sample application code that
 * solves the Euler equations of gas dynamics.  This code illustrates the
 * manner in which a code employing the standard Berger/Oliger AMR algorithm
 * for explicit hydrodynamics can be used in the SAMRAI framework.
 * This class is derived from the algs::HyperbolicPatchStrategy abstract base
 * class which defines the bulk of the interface between the hyperbolic
 * intergration algorithm provided by SAMRAI and the numerical routines
 * specific to Euler.  In particular, this class provides routines which
 * maybe applied to any patch in an AMR patch hierarchy.
 *
 * The numerical routines model the Euler equations of gas dynamics with
 * explicit timestepping and a second-order unsplit Godunov method.
 * The primary numerical quantities are density, velocity, and pressure.
 */

using namespace SAMRAI;

class Euler:
   public tbox::Serializable,
   public algs::HyperbolicPatchStrategy,
   public appu::BoundaryUtilityStrategy,
   public appu::VisDerivedDataStrategy
{
public:
   /**
    * The constructor for Euler sets default parameters for the
    * Euler model.  Specifically, it allocates the variables that represent
    * the state of the solution.  The constructor also registers this
    * object for restart with the restart manager using the object name.
    *
    * After default values are set, this routine calls getFromRestart()
    * if execution from a restart file is specified.  Finally,
    * getFromInput() is called to read values from the given input
    * database (potentially overriding those found in the restart file).
    */
   Euler(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   /**
    * The destructor for Euler does nothing.
    */
   ~Euler();

   ///
   ///  The following routines:
   ///
   ///      registerModelVariables(),
   ///      initializeDataOnPatch(),
   ///      computeStableDtOnPatch(),
   ///      computeFluxesOnPatch(),
   ///      conservativeDifferenceOnPatch(),
   ///      tagGradientDetectorCells(),
   ///      tagRichardsonExtrapolationCells()
   ///
   ///  are concrete implementations of functions declared in the
   ///  algs::HyperbolicPatchStrategy abstract base class.
   ///

   /**
    * Register Euler model variables with algs::HyperbolicLevelIntegrator
    * according to variable registration function provided by the integrator.
    * In other words, variables are registered according to their role
    * in the integration process (e.g., time-dependent, flux, etc.).
    * This routine also registers variables for plotting with the
    * Vis writer.
    */
   void
   registerModelVariables(
      algs::HyperbolicLevelIntegrator* integrator);

   /**
    * Set up parameters in the load balancer object (owned by the gridding
    * algorithm) if needed.  The Euler model allows non-uniform load balancing
    * to be used based on the input file parameter called
    * "use_nonuniform_workload".  The default case is to use uniform
    * load balancing (i.e., use_nonuniform_workload == false).  For
    * illustrative and testing purposes, when non-uniform load balancing is
    * turned on, a weight of one will be applied to every grid cell.  This
    * should produce an identical patch configuration to the uniform load
    * balance case.
    */
   void
   setupLoadBalancer(
      algs::HyperbolicLevelIntegrator* integrator,
      mesh::GriddingAlgorithm* gridding_algorithm);

   /**
    * Set the data on the patch interior to some initial values,
    * depending on the input parameters and numerical routines.
    * If the "initial_time" flag is false, indicating that the
    * routine is called after a regridding step, the routine does nothing.
    */
   void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double data_time,
      const bool initial_time);

   /**
    * Compute the stable time increment for patch using a CFL
    * condition and return the computed dt.
    */
   double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const bool initial_time,
      const double dt_time);

   /**
    * Compute time integral of fluxes to be used in conservative difference
    * for patch integration.  When d_dim == tbox::Dimension(3)), this function calls either
    * compute3DFluxesWithCornerTransport1(), or
    * compute3DFluxesWithCornerTransport2() depending on which
    * transverse flux correction option that is specified in input.
    * The conservative difference used to update the integrated quantities
    * is implemented in the conservativeDifferenceOnPatch() routine.
    */
   void
   computeFluxesOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt);

   /**
    * Update Euler solution variables by performing a conservative
    * difference with the fluxes calculated in computeFluxesOnPatch().
    */
   void
   conservativeDifferenceOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt,
      bool at_syncronization);

   /**
    * Tag cells for refinement using gradient detector.
    */
   void
   tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_indx,
      const bool uses_richardson_extrapolation_too);

   /**
    * Tag cells for refinement using Richardson extrapolation.
    */
   void
   tagRichardsonExtrapolationCells(
      hier::Patch& patch,
      const int error_level_number,
      const std::shared_ptr<hier::VariableContext>& coarsened_fine,
      const std::shared_ptr<hier::VariableContext>& advanced_coarse,
      const double regrid_time,
      const double deltat,
      const int error_coarsen_ratio,
      const bool initial_error,
      const int tag_index,
      const bool uses_gradient_detector_too);

   //@{
   //! @name Required implementations of HyperbolicPatchStrategy pure virtuals.

   ///
   ///  The following routines:
   ///
   ///      setPhysicalBoundaryConditions(),
   ///      getRefineOpStencilWidth(),
   ///      preprocessRefine()
   ///      postprocessRefine()
   ///
   ///  are concrete implementations of functions declared in the
   ///  RefinePatchStrategy abstract base class.
   ///

   /**
    * Set the data in ghost cells corresponding to physical boundary
    * conditions.  Specific boundary conditions are determined by
    * information specified in input file and numerical routines.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector&
      ghost_width_to_fill);

   /**
    * Return stencil width of conservative linear interpolation operations.
    */
   hier::IntVector
   getRefineOpStencilWidth(const tbox::Dimension& dim) const {
      return hier::IntVector::getOne(dim);
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

   /**
    * Refine velocity and pressure from coarse patch to fine patch
    * so that momentum and total energy are conserved.
    */
   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   ///
   ///  The following routines:
   ///
   ///      getCoarsenOpStencilWidth(),
   ///      preprocessCoarsen()
   ///      postprocessCoarsen()
   ///
   ///  are concrete implementations of functions declared in the
   ///  CoarsenPatchStrategy abstract base class.  They are trivial
   ///  because this class doesn't do any pre/postprocessCoarsen.
   ///

   /**
    * Return stencil width of conservative averaging operations.
    */
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

   /**
    * Coarsen velocity and pressure from coarse patch to fine patch
    * so that momentum and total energy are conserved.
    */
   void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   //@}

   /**
    * Write state of Euler object to the given database for restart.
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
    * given name string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary
    * condition applies.
    */
   void
   readDirichletBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It is a blank implementation
    * for the purposes of this class.
    */
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
    * This routine is a concrete implementation of the virtual function
    * in the base class appu::VisDerivedDataStrategy.  It computes derived
    * plot quantities registered with the VisIt data
    * writers from data  that is maintained on each patch in the
    * hierarchy.  In particular, it writes the plot quantity
    * identified by the string variable name to the specified
    * double buffer on the patch in the given region.  The depth_id
    * integer argument indicates which entry in the "depth" of the
    * vector is being written; for a scalar quantity, this may be
    * ignored.  For a vector quantity, it may be used to compute
    * the quantity at the particular depth (e.g. mom[depth_id] =
    * rho * vel[depth_id]).  The boolean return value specifies
    * whether or not derived data exists on the patch.  Generally,
    * this will be TRUE.  If the derived data does NOT exist on
    * the patch, return FALSE.
    */
   bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time) const;

   ///
   ///  The following routines are specific to the Euler class and
   ///  are not declared in any base class.
   ///

   /**
    * Reset physical boundary values in special cases, such as when
    * using symmetric (i.e., reflective) boundary conditions.
    */
   void
   boundaryReset(
      hier::Patch& patch,
      pdat::FaceData<double>& traced_left,
      pdat::FaceData<double>& traced_right) const;

   /**
    * Print all data members for Euler class.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*
    * Dump data in intersection of 1-dimensional "pencil box" to file
    * with given name.  The direction corresponds to the axis of the
    * pencil box in the domain.  Data dumped by this routine is
    * readable by Matlab.
    */
   void
   writeData1dPencil(
      const std::shared_ptr<hier::Patch> patch,
      const hier::Box& pencil_box,
      const tbox::Dimension::dir_t idir,
      std::ostream& file);

private:
   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.   If the boolean flag "is_from_restart"
    * is true when reading from input, some restart values may be
    * overridden by those in the input file.
    *
    * An assertion results if the database pointer is null.
    */
   void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);
   void
   getFromRestart();

   void
   readStateDataEntry(
      std::shared_ptr<tbox::Database> db,
      const std::string& db_name,
      int array_indx,
      std::vector<double>& density,
      std::vector<double>& velocity,
      std::vector<double>& pressure);

   /*
    * Private member function to check correctness of boundary data.
    */
   void
   checkBoundaryData(
      int btype,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const std::vector<int>& scalar_bconds,
      const std::vector<int>& vector_bconds) const;

   /*
    * Three-dimensional flux computation routines corresponding to
    * either of the two transverse flux correction options.  These
    * routines are called from the computeFluxesOnPatch() function.
    */
   void
   compute3DFluxesWithCornerTransport1(
      hier::Patch& patch,
      const double dt);
   void
   compute3DFluxesWithCornerTransport2(
      hier::Patch& patch,
      const double dt);

   /*
    * The object name is used for error/warning reporting and also as a
    * string label for restart database entries.
    */
   std::string d_object_name;

   /*
    * We cache pointers to the grid geometry and Vis data writers
    * to set up initial data, set physical boundary conditions,
    * and register plot variables.  We also cache a pointer to the
    * plot context passed to the variable registration routine.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

#ifdef HAVE_HDF5
   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif
   std::shared_ptr<hier::VariableContext> d_plot_context;

   /*
    * Problem dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * Data items used for nonuniform load balance, if used.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_workload_variable;
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   /*
    * Euler solution state is represented by "primitive" variables,
    * density, velocity, and pressure.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_density;
   std::shared_ptr<pdat::CellVariable<double> > d_velocity;
   std::shared_ptr<pdat::CellVariable<double> > d_pressure;

   /*
    * std::shared_ptr to flux variable vector  - [frho, fu, fp]
    */
   std::shared_ptr<pdat::FaceVariable<double> > d_flux;

   /*
    * Ratio of specific heats for ideal gas.
    */
   double d_gamma;

   /*
    *  Parameters for numerical method:
    *
    *    d_riemann_solve ....... Riemann solver used in flux calculation
    *
    *    d_godunov_order ....... order of Godunov slopes (1, 2, or 4)
    *
    *    d_corner_transport .... type of finite difference approximation
    *                            for 3d transverse flux correction
    *
    *    d_nghosts ............. number of ghost cells for cell-centered
    *                            and face/side-centered variables
    *
    *    d_fluxghosts .......... number of ghost cells for fluxes
    *
    */
   std::string d_riemann_solve;
   int d_riemann_solve_int;
   int d_godunov_order;
   std::string d_corner_transport;
   hier::IntVector d_nghosts;
   hier::IntVector d_fluxghosts;

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
   double d_density_inside;
   double d_velocity_inside[SAMRAI::MAX_DIM_VAL];
   double d_pressure_inside;
   double d_density_outside;
   double d_velocity_outside[SAMRAI::MAX_DIM_VAL];
   double d_pressure_outside;

   /*
    * Input for PIECEWISE_CONSTANT_*  and STEP problems
    */
   int d_number_of_intervals;
   std::vector<double> d_front_position;
   std::vector<double> d_interval_density;
   std::vector<double> d_interval_velocity;
   std::vector<double> d_interval_pressure;

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   std::vector<int> d_master_bdry_edge_conds;
   std::vector<int> d_master_bdry_node_conds;
   std::vector<int> d_master_bdry_face_conds; // Used in 3D only.

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
   std::vector<int> d_scalar_bdry_edge_conds;
   std::vector<int> d_vector_bdry_edge_conds;

   std::vector<int> d_scalar_bdry_node_conds;
   std::vector<int> d_vector_bdry_node_conds;

   std::vector<int> d_scalar_bdry_face_conds; // Used in 3D only.
   std::vector<int> d_vector_bdry_face_conds; // Used in 3D only.

   std::vector<int> d_node_bdry_edge; // Used in 2D only.
   std::vector<int> d_edge_bdry_face; // Used in 3D only.
   std::vector<int> d_node_bdry_face; // Used in 3D only.

   /*
    * Vectors of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
   std::vector<double> d_bdry_edge_density; // Used in 2D only.
   std::vector<double> d_bdry_edge_velocity; // Used in 2D only.
   std::vector<double> d_bdry_edge_pressure; // Used in 2D only.
   std::vector<double> d_bdry_face_density; // Used in 3D only.
   std::vector<double> d_bdry_face_velocity; // Used in 3D only.
   std::vector<double> d_bdry_face_pressure; // Used in 3D only.

   /*
    * Refinement criteria parameters for gradient detector and
    * Richardson extrapolation.
    */
   std::vector<std::string> d_refinement_criteria;
   std::vector<double> d_density_dev_tol;
   std::vector<double> d_density_dev;
   std::vector<double> d_density_dev_time_max;
   std::vector<double> d_density_dev_time_min;
   std::vector<double> d_density_grad_tol;
   std::vector<double> d_density_grad_time_max;
   std::vector<double> d_density_grad_time_min;
   std::vector<double> d_density_shock_onset;
   std::vector<double> d_density_shock_tol;
   std::vector<double> d_density_shock_time_max;
   std::vector<double> d_density_shock_time_min;
   std::vector<double> d_density_rich_tol;
   std::vector<double> d_density_rich_time_max;
   std::vector<double> d_density_rich_time_min;
   std::vector<double> d_pressure_dev_tol;
   std::vector<double> d_pressure_dev;
   std::vector<double> d_pressure_dev_time_max;
   std::vector<double> d_pressure_dev_time_min;
   std::vector<double> d_pressure_grad_tol;
   std::vector<double> d_pressure_grad_time_max;
   std::vector<double> d_pressure_grad_time_min;
   std::vector<double> d_pressure_shock_onset;
   std::vector<double> d_pressure_shock_tol;
   std::vector<double> d_pressure_shock_time_max;
   std::vector<double> d_pressure_shock_time_min;
   std::vector<double> d_pressure_rich_tol;
   std::vector<double> d_pressure_rich_time_max;
   std::vector<double> d_pressure_rich_time_min;

   /*
    * Timers.
    */
   static std::shared_ptr<tbox::Timer> t_init;
   static std::shared_ptr<tbox::Timer> t_compute_dt;
   static std::shared_ptr<tbox::Timer> t_compute_fluxes;
   static std::shared_ptr<tbox::Timer> t_conservdiff;
   static std::shared_ptr<tbox::Timer> t_setphysbcs;
   static std::shared_ptr<tbox::Timer> t_taggradient;

};

#endif
