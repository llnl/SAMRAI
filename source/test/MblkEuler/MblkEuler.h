/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in linear advection ex.
 *
 ************************************************************************/

#ifndef included_MblkEulerXD
#define included_MblkEulerXD

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include <string>
#include <vector>
#define included_String

#include "MblkGeometry.h"
#include "test/testlib/MblkHyperbolicLevelIntegrator.h"
#include "test/testlib/MblkHyperbolicPatchStrategy.h"

// ----------------------------------------------------------------------

using namespace SAMRAI;

class MblkEuler:
   public tbox::Serializable,
   public MblkHyperbolicPatchStrategy,
   public xfer::SingularityPatchStrategy
{
public:
   //
   // the constructor and destructor
   //
   MblkEuler(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<hier::BaseGridGeometry>& grid_geom);

   ~MblkEuler();

   //
   // register with the framework
   //
   void
   registerModelVariables(
      MblkHyperbolicLevelIntegrator* integrator);

   //
   // set the patch initial conditions
   //
   void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double data_time,
      const bool initial_time);

   //
   // Compute the stable time increment for patch using a CFL
   // condition and return the computed dt.
   //
   double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const bool initial_time,
      const double dt_time);

   //
   // compute the state extrema for debugging
   //
   void
   testPatchExtrema(
      hier::Patch& patch,
      const char* pos);

   //
   // compute the fluxes and the initial update in a timestep
   //
   void
   computeFluxesOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt);

   //
   // update the state (currently only for refluxing)
   //
   void
   conservativeDifferenceOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt,
      bool at_syncronization);

   //
   // Tag cells for refinement using gradient detector.
   //
   void
   tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_indexindx,
      const bool uses_richardson_extrapolation_too);

   //
   //  The following routines:
   //
   //      postprocessRefine()
   //      setPhysicalBoundaryConditions()
   //
   //  are concrete implementations of functions declared in the
   //  RefinePatchStrategy abstract base class.
   //

   //
   // mark the zones to track what zones are being filled
   //
   void
   markPhysicalBoundaryConditions(
      hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill);

   //
   // set the data in the physical ghost zones
   //
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector&
      ghost_width_to_fill);

   //
   // Refine operations for cell data.
   //
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

   //
   // Coarsen operations for cell data.
   //
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
      const hier::Patch& patch);

   //
   // build the volume on a mapped grid
   //
   void
   setVolumeOnPatch(
      const hier::Patch& patch);

   /**
    * Write state of MblkEuler object to the given database for restart.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

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
    * Print all data members for MblkEuler class.
    */
   void
   printClassData(
      std::ostream& os) const;

private:
   /*
    * This private member function reads data from input.  If the boolean flag
    * is true all input must be present in input database.
    *
    * An assertion results if the database pointer is null.
    */
   void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);

   void
   getFromRestart();

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

   //
   // Data items used for nonuniform load balance, if used.
   //
   int d_workload_data_id;
   bool d_use_nonuniform_workload;

   //
   // =============== State and Variable definitions (private) ================
   //

   //
   // std::shared_ptr to state variable vector - [state]
   //
   int d_nState;  // depth of the state vector
   std::shared_ptr<pdat::CellVariable<double> > d_state;

   //
   // std::shared_ptr to cell volume - [v]
   //
   std::shared_ptr<pdat::CellVariable<double> > d_vol;

   //
   // std::shared_ptr to flux variable vector  - [F]
   //
   std::shared_ptr<pdat::SideVariable<double> > d_flux;

   //
   // std::shared_ptr to grid - [xyz]
   //
   std::shared_ptr<pdat::NodeVariable<double> > d_xyz;

   //
   // =========================== Initial Conditions (private) ================
   //

   /// center of the sphere or revolution origin
   double d_center[SAMRAI::MAX_DIM_VAL];

   /// revolution axis
   double d_axis[SAMRAI::MAX_DIM_VAL];

   /// revolution radius and pos on axis of radius
   std::vector<std::vector<double> > d_rev_rad;
   std::vector<std::vector<double> > d_rev_axis;

   ///
   /// Rayleigh Taylor Shock tube experiments
   ///
   double d_dt_ampl;
   std::vector<double> d_amn;
   std::vector<double> d_m_mode;
   std::vector<double> d_n_mode;
   std::vector<double> d_phiy;
   std::vector<double> d_phiz;

   ///
   /// input for all the geometries
   ///

   //
   // linear advection velocity vector for unit test
   //
   int d_advection_test;      // run the linear advection unit test
   int d_advection_vel_type;  // type of velocity to use
   double d_advection_velocity[SAMRAI::MAX_DIM_VAL];

   //
   // sizing of zonal, flux, and nodal ghosts
   //
   hier::IntVector d_nghosts;
   hier::IntVector d_fluxghosts;
   hier::IntVector d_nodeghosts;

   //
   // Indicator for problem type and initial conditions
   //
   std::string d_data_problem;

   //
   // region initialization inputs
   //
   int d_number_of_regions;
   std::vector<double> d_front_position;

   //
   // array of initial conditions and their names [region][state]
   //
   std::vector<std::vector<double> > d_state_ic;
   std::vector<std::string> d_state_names;

   //
   // This class stores geometry information used for constructing the
   // mapped multiblock hierarchy
   //
   MblkGeometry* d_mblk_geometry;

   /// the bound on the index space for the current block
   int d_dom_current_bounds[6];

   /// the number of boxes needed to describe the index space
   /// for the current block
   int d_dom_current_nboxes;

   /// the blocks bounding the current patch
   hier::BlockId::block_t d_dom_local_blocks[6];

   //
   // ====================== Refinement Data (private) =======================
   //

   std::vector<std::string> d_refinement_criteria;

   /// history variable gradient tagging tolerance
   std::vector<std::vector<double> > d_state_grad_tol;
   std::vector<std::string> d_state_grad_names;
   std::vector<int> d_state_grad_id;

   //
   // ==================== Boundary Conditions (private) ======================
   //

   /// factors for the boundary conditions
   std::vector<int> d_wall_factors;

   //
   // Operators to be used with GridGeometry
   //
   std::shared_ptr<hier::TimeInterpolateOperator> d_cell_time_interp_op;
};

#endif
