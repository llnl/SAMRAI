/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in linear advection ex.
 *
 ************************************************************************/

#ifndef included_MblkLinAdvXD
#define included_MblkLinAdvXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"

#include "test/testlib/SkeletonCellDoubleConservativeLinearRefine.h"
#include "test/testlib/SkeletonCellDoubleWeightedAverage.h"
#include "SkeletonOutersideDoubleWeightedAverage.h"
#include <string>
#include <vector>
#include <memory>
#define included_String
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisItDataWriter.h"

// Local classes used for this application
#include "MblkGeometry.h"
#include "test/testlib/MblkHyperbolicLevelIntegrator.h"
#include "test/testlib/MblkHyperbolicPatchStrategy.h"


/**
 * The MblkLinAdv class provides routines for a sample application code that
 * solves linear advection.  This code illustrates the manner in which
 * a code employing the standard Berger/Oliger AMR algorithm for
 * explicit hydrodynamics can be used in the SAMRAI framework.
 * This class is derived from the algs::HyperbolicPatchStrategy abstract base
 * class which defines the bulk of the interface between the hyperbolic
 * intergration algorithm provided by SAMRAI and the numerical routines
 * specific to linear advection.  In particular, this class provides routines
 * which maybe applied to any patch in an AMR patch hierarchy.
 *
 * The linear advection problem is simply du/dt + div(a*u) = 0, where
 * "u" is a scalar-valued function and "a" is a constant vector.  The
 * primary numerical quantities are "uval" and "flux", which represent
 * "u" and "a*u", respectively.  All other variables are temporary
 * quantities used in the numerical routines.  The numerical routines
 * use explicit timestepping and a second-order unsplit Godunov method.
 */

using namespace SAMRAI;

class MblkLinAdv:
   public tbox::Serializable,
   public MblkHyperbolicPatchStrategy,
   public xfer::SingularityPatchStrategy,
   public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor for MblkLinAdv sets default parameters for the linear
    * advection model.  Specifically, it creates variables that represent
    * the state of the solution.  The constructor also registers this
    * object for restart with the restart manager using the object name.
    *
    * After default values are set, this routine calls getFromRestart()
    * if execution from a restart file is specified.  Finally,
    * getFromInput() is called to read values from the given input
    * database (potentially overriding those found in the restart file).
    */
   MblkLinAdv(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<hier::BaseGridGeometry>& grid_geoms);

   /**
    * The destructor for MblkLinAdv does nothing.
    */
   ~MblkLinAdv();

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
    * Register MblkLinAdv model variables with algs::HyperbolicLevelIntegrator
    * according to variable registration function provided by the integrator.
    * In other words, variables are registered according to their role
    * in the integration process (e.g., time-dependent, flux, etc.).
    * This routine also registers variables for plotting with the
    * Vis writer.
    */
   void
   registerModelVariables(
      MblkHyperbolicLevelIntegrator* integrator);

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
    * Update linear advection solution variables by performing a conservative
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
      const int tag_indexindx,
      const bool uses_richardson_extrapolation_too);

   ///
   ///  The following routines:
   ///
   ///      postprocessRefine()
   ///      setPhysicalBoundaryConditions()
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
    * Refine operations for cell data.
    */
   void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);
   void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   /**
    * Coarsen operations for cell data.
    */
   void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   /**
    * Fill the singularity conditions for the multi-block case
    */
   void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      std::shared_ptr<const hier::Connector> dst_to_encon,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box,
      const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry);

   /**
    * Build mapped grid on patch
    */
   void
   setMappedGridOnPatch(
      const hier::Patch& patch,
      const int level_number,
      const int block_number);

   /**
    * Write state of MblkLinAdv object to the given database for restart.
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

   void
   readNeumannBoundaryDataEntry(
      const std::shared_ptr<tbox::Database>& db,
      std::string& db_name,
      int bdry_location_index)
   {
      NULL_USE(db);
      NULL_USE(db_name);
      NULL_USE(bdry_location_index);
   }

   hier::IntVector
   getMultiblockRefineOpStencilWidth() const;
   hier::IntVector
   getMultiblockCoarsenOpStencilWidth();

#ifdef HAVE_HDF5
   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt
    * visualization tool.
    */
   void
   registerVisItDataWriter(
      std::shared_ptr<appu::VisItDataWriter> viz_writer);
#endif

   /**
    * Print all data members for MblkLinAdv class.
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

   const tbox::Dimension d_dim;

   /*
    * We cache pointers to the grid geometry and VisIt data writer
    * object to set up initial data, set physical boundary conditions,
    * and register plot variables.
    */
   std::shared_ptr<hier::BaseGridGeometry> d_grid_geometry;
#ifdef HAVE_HDF5
   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif

   /*
    * Data items used for nonuniform load balance, if used.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_workload_variable;
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   /**
    * std::shared_ptr to state variable vector - [u]
    */
   std::shared_ptr<pdat::CellVariable<double> > d_uval;

   /**
    * std::shared_ptr to flux variable vector  - [F]
    */
   std::shared_ptr<pdat::FaceVariable<double> > d_flux;

   /**
    * std::shared_ptr to grid - [xyz]
    */
   std::shared_ptr<pdat::NodeVariable<double> > d_xyz;

   /**
    * linear advection velocity vector
    */
   double d_advection_velocity[SAMRAI::MAX_DIM_VAL];

   /*
    *  Parameters for numerical method:
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
   int d_godunov_order;
   std::string d_corner_transport;
   hier::IntVector d_nghosts;
   hier::IntVector d_fluxghosts;
   hier::IntVector d_nodeghosts;

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
   double d_uval_inside;
   double d_uval_outside;

   /*
    * Input for FRONT problem
    */
   int d_number_of_intervals;
   std::vector<double> d_front_position;
   std::vector<double> d_interval_uval;

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   std::vector<int> d_scalar_bdry_edge_conds;
   std::vector<int> d_scalar_bdry_node_conds;
   std::vector<int> d_scalar_bdry_face_conds; // 3D only

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
   std::vector<int> d_node_bdry_edge; // 2D only
   std::vector<int> d_edge_bdry_face; // 3D only
   std::vector<int> d_node_bdry_face; // 3D only

   /*
    * Vectors of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
   std::vector<double> d_bdry_edge_uval;
   std::vector<double> d_bdry_face_uval;

   /*
    * Input for Sine problem initialization
    */
   double d_amplitude;
   double d_frequency[SAMRAI::MAX_DIM_VAL];

   /*
    * Refinement criteria parameters for gradient detector and
    * Richardson extrapolation.
    */
   std::vector<std::string> d_refinement_criteria;
   std::vector<double> d_dev_tol;
   std::vector<double> d_dev;
   std::vector<double> d_dev_time_max;
   std::vector<double> d_dev_time_min;
   std::vector<double> d_grad_tol;
   std::vector<double> d_grad_time_max;
   std::vector<double> d_grad_time_min;
   std::vector<double> d_shock_onset;
   std::vector<double> d_shock_tol;
   std::vector<double> d_shock_time_max;
   std::vector<double> d_shock_time_min;
   std::vector<double> d_rich_tol;
   std::vector<double> d_rich_time_max;
   std::vector<double> d_rich_time_min;

   /*
    * This class stores geometry information used for constructing the
    * mapped multiblock hierarchy
    */
   MblkGeometry* d_mblk_geometry;

   /*
    * Operators to be used with GridGeometry
    */
   std::shared_ptr<SkeletonCellDoubleConservativeLinearRefine>
   d_cell_cons_linear_refine_op;
   std::shared_ptr<SkeletonCellDoubleWeightedAverage> d_cell_cons_coarsen_op;
   std::shared_ptr<hier::TimeInterpolateOperator> d_cell_time_interp_op;
   std::shared_ptr<SkeletonOutersideDoubleWeightedAverage>
   d_side_cons_coarsen_op;

};

#endif
