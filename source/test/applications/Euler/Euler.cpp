/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for Euler equations SAMRAI example
 *
 ************************************************************************/

#include "Euler.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (0)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

// External definitions for Fortran numerical routines
#include "EulerFort.h"

// Number of entries in state vector (d_dim velocity comps + pressure + density)
#define NEQU (d_dim.getValue() + 2)

// Number of ghosts cells used for each variable quantity.
#define CELLG (4)
#define FACEG (4)
#define FLUXG (1)

// defines for initialization
#define PIECEWISE_CONSTANT_X (10)
#define PIECEWISE_CONSTANT_Y (11)
#define PIECEWISE_CONSTANT_Z (12)
#define SPHERE (40)
#define STEP (80)

// defines for Riemann solver used in Godunov flux calculation
#define APPROX_RIEM_SOLVE (20)   // Colella-Glaz approx Riemann solver
#define EXACT_RIEM_SOLVE (21)    // Exact Riemann solver
#define HLLC_RIEM_SOLVE (22)     // Harten, Lax, van Leer approx Riemann solver

// defines for cell tagging routines
#define RICHARDSON_NEWLY_TAGGED (-10)
#define RICHARDSON_ALREADY_TAGGED (-11)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

// Version of Euler restart file data
#define EULER_VERSION (3)

std::shared_ptr<tbox::Timer> Euler::t_init;
std::shared_ptr<tbox::Timer> Euler::t_compute_dt;
std::shared_ptr<tbox::Timer> Euler::t_compute_fluxes;
std::shared_ptr<tbox::Timer> Euler::t_conservdiff;
std::shared_ptr<tbox::Timer> Euler::t_setphysbcs;
std::shared_ptr<tbox::Timer> Euler::t_taggradient;

/*
 *************************************************************************
 *
 * The constructor for Euler class sets data members to defualt values,
 * creates variables that define the solution state for the Euler
 * equations.
 *
 * After default values are set, this routine calls getFromRestart()
 * if execution from a restart file is specified.  Finally,
 * getFromInput() is called to read values from the given input
 * database (potentially overriding those found in the restart file).
 *
 *************************************************************************
 */

Euler::Euler(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<geom::CartesianGridGeometry> grid_geom):
   algs::HyperbolicPatchStrategy(),
   d_object_name(object_name),
   d_grid_geometry(grid_geom),
   d_dim(dim),
   d_use_nonuniform_workload(false),
   d_density(new pdat::CellVariable<double>(dim, "density", 1)),
   d_velocity(new pdat::CellVariable<double>(
                 dim, "velocity", d_dim.getValue())),
   d_pressure(new pdat::CellVariable<double>(dim, "pressure", 1)),
   d_flux(new pdat::FaceVariable<double>(dim, "flux", NEQU)),
   d_gamma(1.4),
   // specific heat ratio for ideal diatomic gas (e.g., air)
   d_riemann_solve("APPROX_RIEM_SOLVE"),
   d_godunov_order(1),
   d_corner_transport("CORNER_TRANSPORT_1"),
   d_nghosts(hier::IntVector(dim, CELLG)),
   d_fluxghosts(hier::IntVector(dim, FLUXG)),
   d_radius(tbox::MathUtilities<double>::getSignalingNaN()),
   d_density_inside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_pressure_inside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_density_outside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_pressure_outside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_number_of_intervals(0)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(grid_geom);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   if (!t_init) {
      t_init = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::initializeDataOnPatch()");
      t_compute_dt = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::computeStableDtOnPatch()");
      t_compute_fluxes = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::computeFluxesOnPatch()");
      t_conservdiff = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::conservativeDifferenceOnPatch()");
      t_setphysbcs = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::setPhysicalBoundaryConditions()");
      t_taggradient = tbox::TimerManager::getManager()->
         getTimer("apps::Euler::tagGradientDetectorCells()");
   }

   TBOX_ASSERT(CELLG == FACEG);

   /*
    * Defaults for problem type and initial data
    */

   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center,
      d_dim.getValue());
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_velocity_inside,
      d_dim.getValue());
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_velocity_outside,
      d_dim.getValue());

   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

   if (d_dim == tbox::Dimension(2)) {
      d_master_bdry_edge_conds.resize(NUM_2D_EDGES);
      d_scalar_bdry_edge_conds.resize(NUM_2D_EDGES);
      d_vector_bdry_edge_conds.resize(NUM_2D_EDGES);
      for (int ei = 0; ei < NUM_2D_EDGES; ++ei) {
         d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_node_conds.resize(NUM_2D_NODES);
      d_scalar_bdry_node_conds.resize(NUM_2D_NODES);
      d_vector_bdry_node_conds.resize(NUM_2D_NODES);
      d_node_bdry_edge.resize(NUM_2D_NODES);

      for (int ni = 0; ni < NUM_2D_NODES; ++ni) {
         d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
      }

      d_bdry_edge_density.resize(NUM_2D_EDGES);
      d_bdry_edge_velocity.resize(NUM_2D_EDGES * d_dim.getValue());
      d_bdry_edge_pressure.resize(NUM_2D_EDGES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_edge_density);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_edge_velocity);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_edge_pressure);
   }
   if (d_dim == tbox::Dimension(3)) {
      d_master_bdry_face_conds.resize(NUM_3D_FACES);
      d_scalar_bdry_face_conds.resize(NUM_3D_FACES);
      d_vector_bdry_face_conds.resize(NUM_3D_FACES);
      for (int fi = 0; fi < NUM_3D_FACES; ++fi) {
         d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
         d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
         d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_edge_conds.resize(NUM_3D_EDGES);
      d_scalar_bdry_edge_conds.resize(NUM_3D_EDGES);
      d_vector_bdry_edge_conds.resize(NUM_3D_EDGES);
      d_edge_bdry_face.resize(NUM_3D_EDGES);
      for (int ei = 0; ei < NUM_3D_EDGES; ++ei) {
         d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_node_conds.resize(NUM_3D_NODES);
      d_scalar_bdry_node_conds.resize(NUM_3D_NODES);
      d_vector_bdry_node_conds.resize(NUM_3D_NODES);
      d_node_bdry_face.resize(NUM_3D_NODES);

      for (int ni = 0; ni < NUM_3D_NODES; ++ni) {
         d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
      }

      d_bdry_face_density.resize(NUM_3D_FACES);
      d_bdry_face_velocity.resize(NUM_3D_FACES * d_dim.getValue());
      d_bdry_face_pressure.resize(NUM_3D_FACES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_face_density);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_face_velocity);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(
         d_bdry_face_pressure);
   }

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   /*
    * Set problem data to values read from input/restart.
    */

   if (d_riemann_solve == "APPROX_RIEM_SOLVE") {
      d_riemann_solve_int = APPROX_RIEM_SOLVE;
   } else if (d_riemann_solve == "EXACT_RIEM_SOLVE") {
      d_riemann_solve_int = EXACT_RIEM_SOLVE;
   } else if (d_riemann_solve == "HLLC_RIEM_SOLVE") {
      d_riemann_solve_int = HLLC_RIEM_SOLVE;
   } else {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Unknown d_riemann_solve string = "
                       << d_riemann_solve
                       << " encountered in constructor" << std::endl);
   }

   if (d_data_problem == "PIECEWISE_CONSTANT_X") {
      d_data_problem_int = PIECEWISE_CONSTANT_X;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
      d_data_problem_int = PIECEWISE_CONSTANT_Y;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
      d_data_problem_int = PIECEWISE_CONSTANT_Z;
   } else if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
   } else if (d_data_problem == "STEP") {
      d_data_problem_int = STEP;
   } else {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Unknown d_data_problem string = "
                       << d_data_problem
                       << " encountered in constructor" << std::endl);
   }

   /*
    * Postprocess boundary data from input/restart values.
    */
   if (d_dim == tbox::Dimension(2)) {
      for (int i = 0; i < NUM_2D_EDGES; ++i) {
         d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
         d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

         if (d_master_bdry_edge_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::FLOW;
         }
      }

      for (int i = 0; i < NUM_2D_NODES; ++i) {
         d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
         d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

         if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }

         if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_edge[i] =
               appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                  i, d_master_bdry_node_conds[i]);
         }
      }
   }
   if (d_dim == tbox::Dimension(3)) {
      for (int i = 0; i < NUM_3D_FACES; ++i) {
         d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
         d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];

         if (d_master_bdry_face_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_face_conds[i] = BdryCond::FLOW;
         }
      }

      for (int i = 0; i < NUM_3D_EDGES; ++i) {
         d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
         d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

         if (d_master_bdry_edge_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_edge_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::YFLOW;
         }
         if (d_master_bdry_edge_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::ZFLOW;
         }

         if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
            d_edge_bdry_face[i] =
               appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                  i, d_master_bdry_edge_conds[i]);
         }
      }

      for (int i = 0; i < NUM_3D_NODES; ++i) {
         d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
         d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

         if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::ZFLOW;
         }

         if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_face[i] =
               appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                  i, d_master_bdry_node_conds[i]);
         }
      }

   }

   SAMRAI_F77_FUNC(stufprobc, STUFPROBC) (APPROX_RIEM_SOLVE, EXACT_RIEM_SOLVE,
      HLLC_RIEM_SOLVE,
      PIECEWISE_CONSTANT_X, PIECEWISE_CONSTANT_Y, PIECEWISE_CONSTANT_Z,
      SPHERE, STEP,
      CELLG, FACEG, FLUXG);
}

/*
 *************************************************************************
 *
 * Empty destructor for Euler class.
 *
 *************************************************************************
 */

Euler::~Euler()
{
   t_init.reset();
   t_compute_dt.reset();
   t_compute_fluxes.reset();
   t_conservdiff.reset();
   t_setphysbcs.reset();
   t_taggradient.reset();
}

/*
 *************************************************************************
 *
 * Register density, velocity, pressure (i.e., solution state variables),
 * and flux variables with hyperbolic integrator that manages storage
 * for those quantities.  Also, register plot data with the vis tool.
 *
 * Note that density coarsening/refining uses standard conservative
 * operations provided in SAMRAI library.   Velocity and pressure
 * are not conserved.  The Euler code provides operations to coarsen/
 * refine momentum and total energy conservatively.  Velocity and
 * pressure are calculated from the conserved quantities.
 *
 *************************************************************************
 */

void Euler::registerModelVariables(
   algs::HyperbolicLevelIntegrator* integrator)
{
   TBOX_ASSERT(integrator != 0);
   TBOX_ASSERT(CELLG == FACEG);

   integrator->registerVariable(d_density, d_nghosts,
      algs::HyperbolicLevelIntegrator::TIME_DEP,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(d_velocity, d_nghosts,
      algs::HyperbolicLevelIntegrator::TIME_DEP,
      d_grid_geometry,
      "USER_DEFINED_COARSEN",
      "USER_DEFINED_REFINE");

   integrator->registerVariable(d_pressure, d_nghosts,
      algs::HyperbolicLevelIntegrator::TIME_DEP,
      d_grid_geometry,
      "USER_DEFINED_COARSEN",
      "USER_DEFINED_REFINE");

   integrator->registerVariable(d_flux, d_fluxghosts,
      algs::HyperbolicLevelIntegrator::FLUX,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "NO_REFINE");

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

   d_plot_context = integrator->getPlotContext();

#ifdef HAVE_HDF5
   if (d_visit_writer) {
      d_visit_writer->registerPlotQuantity("Density",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_density, d_plot_context));

      d_visit_writer->registerPlotQuantity("Velocity",
         "VECTOR",
         vardb->mapVariableAndContextToIndex(
            d_velocity, d_plot_context));

      d_visit_writer->registerPlotQuantity("Pressure",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_pressure, d_plot_context));

      d_visit_writer->registerDerivedPlotQuantity("Total Energy",
         "SCALAR",
         this);
      d_visit_writer->registerDerivedPlotQuantity("Momentum",
         "VECTOR",
         this);
   }

   if (!d_visit_writer) {
      TBOX_WARNING(d_object_name << ": registerModelVariables()\n"
                                 << "VisIt data writer was not registered\n"
                                 << "Consequently, no plot data will\n"
                                 << "be written." << std::endl);
   }
#endif

}

/*
 *************************************************************************
 *
 * Set up parameters for nonuniform load balancing, if used.
 *
 *************************************************************************
 */

void Euler::setupLoadBalancer(
   algs::HyperbolicLevelIntegrator* integrator,
   mesh::GriddingAlgorithm* gridding_algorithm)
{
   NULL_USE(integrator);

   const hier::IntVector& zero_vec = hier::IntVector::getZero(d_dim);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
   hier::PatchDataRestartManager* pdrm =
      hier::PatchDataRestartManager::getManager();

   if (d_use_nonuniform_workload && gridding_algorithm) {
      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         std::dynamic_pointer_cast<mesh::TreeLoadBalancer, mesh::LoadBalanceStrategy>(
            gridding_algorithm->getLoadBalanceStrategy()));

      if (load_balancer) {
         d_workload_variable.reset(
            new pdat::CellVariable<double>(
               d_dim,
               "workload_variable",
               1));
         d_workload_data_id =
            vardb->registerVariableAndContext(d_workload_variable,
               vardb->getContext("WORKLOAD"),
               zero_vec);
         load_balancer->setWorkloadPatchDataIndex(d_workload_data_id);
         pdrm->registerPatchDataForRestart(d_workload_data_id);
      } else {
         TBOX_WARNING(
            d_object_name << ": "
                          << "  Unknown load balancer used in gridding algorithm."
                          << "  Ignoring request for nonuniform load balancing." << std::endl);
         d_use_nonuniform_workload = false;
      }
   } else {
      d_use_nonuniform_workload = false;
   }

}

/*
 *************************************************************************
 *
 * Set initial data for solution variables on patch interior.
 * This routine is called whenever a new patch is introduced to the
 * AMR patch hierarchy.  Note that the routine does nothing unless
 * we are at the initial time.  In all other cases, conservative
 * interpolation from coarser levels and copies from patches at the
 * same mesh resolution are sufficient to set data.
 *
 *************************************************************************
 */

void Euler::initializeDataOnPatch(
   hier::Patch& patch,
   const double data_time,
   const bool initial_time)
{
   NULL_USE(data_time);

   t_init->start();

   if (initial_time) {

      const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(pgeom);
      const double* dx = pgeom->getDx();
      const double* xlo = pgeom->getXLower();
      const double* xhi = pgeom->getXUpper();

      std::shared_ptr<pdat::CellData<double> > density(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_density, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > velocity(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_velocity, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > pressure(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_pressure, getDataContext())));

      TBOX_ASSERT(density);
      TBOX_ASSERT(velocity);
      TBOX_ASSERT(pressure);

      const hier::IntVector& ghost_cells = density->getGhostCellWidth();

      TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
      TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();

      if (d_data_problem == "SPHERE") {

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(eulerinitsphere2d, EULERINITSPHERE2d) (d_data_problem_int,
               dx, xlo, xhi,
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ghost_cells(0),
               ghost_cells(1),
               d_gamma,
               density->getPointer(),
               velocity->getPointer(),
               pressure->getPointer(),
               d_density_inside,
               d_velocity_inside,
               d_pressure_inside,
               d_density_outside,
               d_velocity_outside,
               d_pressure_outside,
               d_center, d_radius);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(eulerinitsphere3d, EULERINITSPHERE3d) (d_data_problem_int,
               dx, xlo, xhi,
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               ghost_cells(0),
               ghost_cells(1),
               ghost_cells(2),
               d_gamma,
               density->getPointer(),
               velocity->getPointer(),
               pressure->getPointer(),
               d_density_inside,
               d_velocity_inside,
               d_pressure_inside,
               d_density_outside,
               d_velocity_outside,
               d_pressure_outside,
               d_center, d_radius);
         }

      } else {

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(eulerinit2d, EULERINIT2D) (d_data_problem_int,
               dx, xlo, xhi,
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ghost_cells(0),
               ghost_cells(1),
               d_gamma,
               density->getPointer(),
               velocity->getPointer(),
               pressure->getPointer(),
               d_number_of_intervals,
               &d_front_position[0],
               &d_interval_density[0],
               &d_interval_velocity[0],
               &d_interval_pressure[0]);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(eulerinit3d, EULERINIT3D) (d_data_problem_int,
               dx, xlo, xhi,
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               ghost_cells(0),
               ghost_cells(1),
               ghost_cells(2),
               d_gamma,
               density->getPointer(),
               velocity->getPointer(),
               pressure->getPointer(),
               d_number_of_intervals,
               &d_front_position[0],
               &d_interval_density[0],
               &d_interval_velocity[0],
               &d_interval_pressure[0]);
         }

      }

   }

   if (d_use_nonuniform_workload) {
      if (!patch.checkAllocated(d_workload_data_id)) {
         patch.allocatePatchData(d_workload_data_id);
      }
      std::shared_ptr<pdat::CellData<double> > workload_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_workload_data_id)));
      TBOX_ASSERT(workload_data);
      workload_data->fillAll(1.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }

   t_init->stop();

}

/*
 *************************************************************************
 *
 * Compute stable time increment for patch.  Return this value.
 *
 *************************************************************************
 */

double Euler::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time,
   const double dt_time)
{
   NULL_USE(initial_time);
   NULL_USE(dt_time);

   t_compute_dt->start();

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   const std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, getDataContext())));
   const std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, getDataContext())));
   const std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, getDataContext())));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);

   const hier::IntVector& ghost_cells = density->getGhostCellWidth();

   TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
   TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);

   double stabdt = 0.;
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(stabledt2d, STABLEDT2D) (dx,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ghost_cells(0),
         ghost_cells(1),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         stabdt);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(stabledt3d, STABLEDT3D) (dx,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         ghost_cells(0),
         ghost_cells(1),
         ghost_cells(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         stabdt);
   }

   t_compute_dt->stop();
   return stabdt;
}

/*
 *************************************************************************
 *
 * Compute time integral of numerical fluxes for finite difference
 * at each cell face on patch.  When d_dim == tbox::Dimension(3)), there are two options
 * for the transverse flux correction.  Otherwise, there is only one.
 *
 *************************************************************************
 */

void Euler::computeFluxesOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt)
{
   NULL_USE(time);
   t_compute_fluxes->start();

   if (d_dim == tbox::Dimension(3)) {
      if (d_corner_transport == "CORNER_TRANSPORT_2") {
         compute3DFluxesWithCornerTransport2(patch, dt);
      } else {
         compute3DFluxesWithCornerTransport1(patch, dt);
      }
   }

   if (d_dim == tbox::Dimension(2)) {

      TBOX_ASSERT(CELLG == FACEG);

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(patch_geom);
      const double* dx = patch_geom->getDx();

      hier::Box pbox = patch.getBox();
      const hier::Index ifirst = pbox.lower();
      const hier::Index ilast = pbox.upper();

      std::shared_ptr<pdat::CellData<double> > density(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_density, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > velocity(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_velocity, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > pressure(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_pressure, getDataContext())));
      std::shared_ptr<pdat::FaceData<double> > flux(
         SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
            patch.getPatchData(d_flux, getDataContext())));

      TBOX_ASSERT(density);
      TBOX_ASSERT(velocity);
      TBOX_ASSERT(pressure);
      TBOX_ASSERT(flux);
      TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
      TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
      TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
      TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

      /*
       * Allocate patch data for temporaries local to this routine.
       */
      pdat::FaceData<double> traced_left(pbox, NEQU, d_nghosts);
      pdat::FaceData<double> traced_right(pbox, NEQU, d_nghosts);
      pdat::CellData<double> sound_speed(pbox, 1, d_nghosts);

      /*
       *  Initialize traced states (w^R and w^L) with proper cell-centered values.
       */
      SAMRAI_F77_FUNC(inittraceflux2d, INITTRACEFLUX2D) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_right.getPointer(0),
         traced_right.getPointer(1),
         flux->getPointer(0),
         flux->getPointer(1));

      /*
       * If Godunov method requires slopes with order greater than one, perform
       * characteristic tracing to compute higher-order slopes.
       */
      if (d_godunov_order > 1) {

         /*
          * Prepare temporary data for characteristic tracing.
          */
         int Mcells = 0;
         for (tbox::Dimension::dir_t k = 0; k < d_dim.getValue(); ++k) {
            Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
         }

         // Face-centered temporary arrays
         std::vector<double> ttedgslp((2 * FACEG + 1 + Mcells) * NEQU);
         std::vector<double> ttraclft((2 * FACEG + 1 + Mcells) * NEQU);
         std::vector<double> ttracrgt((2 * FACEG + 1 + Mcells) * NEQU);

         // Cell-centered temporary arrays
         std::vector<double> ttsound((2 * CELLG + Mcells));
         std::vector<double> ttcelslp((2 * CELLG + Mcells) * NEQU);

         /*
          * Compute local sound speed in each computational cell.
          */
         SAMRAI_F77_FUNC(computesound2d, COMPUTESOUND2D) (ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            d_gamma,
            density->getPointer(),
            velocity->getPointer(),
            pressure->getPointer(),
            sound_speed.getPointer());

         /*
          *  Apply characteristic tracing to compute initial estimate of
          *  traces w^L and w^R at faces.
          *  Inputs: sound_speed, w^L, w^R (traced_left/right)
          *  Output: w^L, w^R
          */
         SAMRAI_F77_FUNC(chartracing2d0, CHARTRACING2D0) (dt,
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            Mcells, dx[0], d_gamma, d_godunov_order,
            sound_speed.getPointer(),
            traced_left.getPointer(0),
            traced_right.getPointer(0),
            &ttcelslp[0],
            &ttedgslp[0],
            &ttsound[0],
            &ttraclft[0],
            &ttracrgt[0]);

         SAMRAI_F77_FUNC(chartracing2d1, CHARTRACING2D1) (dt,
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            Mcells, dx[1], d_gamma, d_godunov_order,
            sound_speed.getPointer(),
            traced_left.getPointer(1),
            traced_right.getPointer(1),
            &ttcelslp[0],
            &ttedgslp[0],
            &ttsound[0],
            &ttraclft[0],
            &ttracrgt[0]);

      }  // if (d_godunov_order > 1) ...

// SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,1,dx, to get artificial viscosity
// SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,0,dx, to get NO artificial viscosity

      /*
       *  Compute preliminary fluxes at faces by solving approximate
       *  Riemann problem using the trace states computed so far.
       *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
       *  Output: F (flux)
       */
      SAMRAI_F77_FUNC(fluxcalculation2d, FLUXCALCULATION2D) (dt, 1, 0, dx,
         ifirst(0), ilast(0), ifirst(1), ilast(1),
         d_gamma,
         d_riemann_solve_int,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         flux->getPointer(0),
         flux->getPointer(1),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_right.getPointer(0),
         traced_right.getPointer(1));

      /*
       *  Update trace states at cell faces with transverse correction applied.
       *  Inputs: F (flux)
       *  Output: w^L, w^R (traced_left/right)
       */
      SAMRAI_F77_FUNC(fluxcorrec, FLUXCORREC) (dt,
         ifirst(0), ilast(0), ifirst(1), ilast(1),
         dx, d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         flux->getPointer(0),
         flux->getPointer(1),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_right.getPointer(0),
         traced_right.getPointer(1));

      boundaryReset(patch, traced_left, traced_right);

      /*
       *  Re-compute fluxes with updated trace states.
       *  Inputs: w^L, w^R (traced_left/right)
       *  Output: F (flux)
       */
      SAMRAI_F77_FUNC(fluxcalculation2d, FLUXCALCULATION2D) (dt, 0, 0, dx,
         ifirst(0), ilast(0), ifirst(1), ilast(1),
         d_gamma,
         d_riemann_solve_int,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         flux->getPointer(0),
         flux->getPointer(1),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_right.getPointer(0),
         traced_right.getPointer(1));

   }

   t_compute_fluxes->stop();
}

/*
 *************************************************************************
 *
 * Compute numerical approximations to flux terms using an extension
 * to three dimensions of Collella's corner transport upwind approach.
 * I.E. input value corner_transport = "CORNER_TRANSPORT_1"
 *
 *************************************************************************
 */

void Euler::compute3DFluxesWithCornerTransport1(
   hier::Patch& patch,
   const double dt)
{
   TBOX_ASSERT(CELLG == FACEG);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, getDataContext())));
   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<double> traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<double> traced_right(pbox, NEQU, d_nghosts);
   pdat::CellData<double> sound_speed(pbox, 1, d_nghosts);
   pdat::FaceData<double> temp_flux(pbox, NEQU, d_fluxghosts);
   pdat::FaceData<double> temp_traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<double> temp_traced_right(pbox, NEQU, d_nghosts);

   /*
    *  Initialize traced states (w^R and w^L) with proper cell-centered values.
    */
   SAMRAI_F77_FUNC(inittraceflux3d, INITTRACEFLUX3D) (ifirst(0), ilast(0),
      ifirst(1), ilast(1),
      ifirst(2), ilast(2),
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2));

   /*
    * If Godunov method requires slopes with order greater than one, perform
    * characteristic tracing to compute higher-order slopes.
    */
   if (d_godunov_order > 1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (tbox::Dimension::dir_t k = 0; k < d_dim.getValue(); ++k) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      std::vector<double> ttedgslp((2 * FACEG + 1 + Mcells) * NEQU);
      std::vector<double> ttraclft((2 * FACEG + 1 + Mcells) * NEQU);
      std::vector<double> ttracrgt((2 * FACEG + 1 + Mcells) * NEQU);

      // Cell-centered temporary arrays
      std::vector<double> ttsound((2 * CELLG + Mcells));
      std::vector<double> ttcelslp((2 * CELLG + Mcells) * NEQU);

      /*
       * Compute local sound speed in each computational cell.
       */
      SAMRAI_F77_FUNC(computesound3d, COMPUTESOUND3D) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         sound_speed.getPointer());

      /*
       *  Apply characteristic tracing to compute initial estimate of
       *  traces w^L and w^R at faces.
       *  Inputs: sound_speed, w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      SAMRAI_F77_FUNC(chartracing3d0, CHARTRACING3D0) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[0], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(0),
         traced_right.getPointer(0),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d1, CHARTRACING3D1) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[1], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(1),
         traced_right.getPointer(1),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d2, CHARTRACING3D2) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[2], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(2),
         traced_right.getPointer(2),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

   }  // if (d_godunov_order > 1) ...

   /*
    *  Compute preliminary fluxes at faces by solving approximate
    *  Riemann problem using the trace states computed so far.
    *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,1,dx,  to do artificial viscosity
//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,0,dx,  to do NO artificial viscosity
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

   /*
    *  Re-compute face traces to include one set of correction terms with
    *  transverse flux differences.  Store result in temporary vectors
    *  (i.e. temp_traced_left/right).
    *  Inputs: F (flux), P, rho, v, w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   SAMRAI_F77_FUNC(fluxcorrec2d, FLUXCORREC2D) (dt,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      dx, d_gamma, 1,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2),
      temp_traced_left.getPointer(0),
      temp_traced_left.getPointer(1),
      temp_traced_left.getPointer(2),
      temp_traced_right.getPointer(0),
      temp_traced_right.getPointer(1),
      temp_traced_right.getPointer(2));

   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Compute fluxes with partially-corrected trace states.  Store
    *  result in temporary flux vector.
    *  Inputs: P, rho, v, temp_traced_left/right
    *  Output: temp_flux
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 0, 1, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      temp_flux.getPointer(0),
      temp_flux.getPointer(1),
      temp_flux.getPointer(2),
      temp_traced_left.getPointer(0),
      temp_traced_left.getPointer(1),
      temp_traced_left.getPointer(2),
      temp_traced_right.getPointer(0),
      temp_traced_right.getPointer(1),
      temp_traced_right.getPointer(2));
   /*
    *  Compute face traces with other transverse correction flux
    *  difference terms included.  Store result in temporary vectors
    *  (i.e. temp_traced_left/right).
    *  Inputs: F (flux), P, rho, v, w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   SAMRAI_F77_FUNC(fluxcorrec2d, FLUXCORREC2D) (dt,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      dx, d_gamma, -1,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2),
      temp_traced_left.getPointer(0),
      temp_traced_left.getPointer(1),
      temp_traced_left.getPointer(2),
      temp_traced_right.getPointer(0),
      temp_traced_right.getPointer(1),
      temp_traced_right.getPointer(2));

   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Compute final predicted fluxes with both sets of transverse flux
    *  differences included.  Store the result in regular flux vector.
    *  NOTE:  the fact that we store  these fluxes in the regular (i.e.
    *  not temporary) flux vector does NOT indicate this is the final result.
    *  Rather, the flux vector is used as a convenient storage location.
    *  Inputs: P, rho, v, temp_traced_left/right
    *  Output: flux
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      temp_traced_left.getPointer(0),
      temp_traced_left.getPointer(1),
      temp_traced_left.getPointer(2),
      temp_traced_right.getPointer(0),
      temp_traced_right.getPointer(1),
      temp_traced_right.getPointer(2));
   /*
    *  Compute the final trace state vectors at cell faces using transverse
    *  differences of final predicted fluxes.  Store the result in w^L
    *  (traced_left) and w^R (traced_right) vectors.
    *  Inputs: temp_flux, flux
    *  Output: w^L, w^R (traced_left/right)
    */
   SAMRAI_F77_FUNC(fluxcorrec3d, FLUXCORREC3D) (dt,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      dx, d_gamma,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      temp_flux.getPointer(0),
      temp_flux.getPointer(1),
      temp_flux.getPointer(2),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

   /*
    *  Final flux calculation using corrected trace states.
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 0, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));
}

/*
 *************************************************************************
 *
 * Compute numerical approximations to flux terms using John
 * Trangenstein's interpretation of the three-dimensional version of
 * Collella's corner transport upwind approach.
 * I.E. input value corner_transport = "CORNER_TRANSPORT_2"
 *
 *************************************************************************
 */

void Euler::compute3DFluxesWithCornerTransport2(
   hier::Patch& patch,
   const double dt)
{
   TBOX_ASSERT(CELLG == FACEG);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, getDataContext())));
   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<double> traced_left(pbox, NEQU, d_nghosts);
   pdat::FaceData<double> traced_right(pbox, NEQU, d_nghosts);
   pdat::CellData<double> sound_speed(pbox, 1, d_nghosts);
   pdat::FaceData<double> temp_flux(pbox, NEQU, d_fluxghosts);
   pdat::CellData<double> third_state(pbox, NEQU, d_nghosts);

   /*
    *  Initialize traced states (w^R and w^L) with proper cell-centered values.
    */
   SAMRAI_F77_FUNC(inittraceflux3d, INITTRACEFLUX3D) (ifirst(0), ilast(0),
      ifirst(1), ilast(1),
      ifirst(2), ilast(2),
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2));

   /*
    *  Compute fluxes at faces by solving approximate Riemann problem
    *  using initial trace states.
    *  Inputs: P, rho, v, w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,1,dx,  to do artificial viscosity
//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,0,dx,  to do NO artificial viscosity
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 1, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

   /*
    * If Godunov method requires slopes with order greater than one, perform
    * characteristic tracing to compute higher-order slopes.
    */
   if (d_godunov_order > 1) {

      /*
       * Prepare temporary data for characteristic tracing.
       */
      int Mcells = 0;
      for (tbox::Dimension::dir_t k = 0; k < d_dim.getValue(); ++k) {
         Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
      }

      // Face-centered temporary arrays
      std::vector<double> ttedgslp((2 * FACEG + 1 + Mcells) * NEQU);
      std::vector<double> ttraclft((2 * FACEG + 1 + Mcells) * NEQU);
      std::vector<double> ttracrgt((2 * FACEG + 1 + Mcells) * NEQU);

      // Cell-centered temporary arrays
      std::vector<double> ttsound((2 * CELLG + Mcells));
      std::vector<double> ttcelslp((2 * CELLG + Mcells) * NEQU);

      /*
       * Compute local sound speed in each computational cell.
       */
      SAMRAI_F77_FUNC(computesound3d, COMPUTESOUND3D) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         sound_speed.getPointer());

      /*
       *  Apply characteristic tracing to update traces w^L and w^R at faces.
       *  Inputs: sound_speed, w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      SAMRAI_F77_FUNC(chartracing3d0, CHARTRACING3D0) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[0], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(0),
         traced_right.getPointer(0),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d1, CHARTRACING3D1) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[1], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(1),
         traced_right.getPointer(1),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d2, CHARTRACING3D2) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[2], d_gamma, d_godunov_order,
         sound_speed.getPointer(),
         traced_left.getPointer(2),
         traced_right.getPointer(2),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttsound[0],
         &ttraclft[0],
         &ttracrgt[0]);

   } // if (d_godunov_order > 1) ...

   for (int idir = 0; idir < d_dim.getValue(); ++idir) {

      /*
       *    Approximate traces at cell centers (in idir direction);
       * i.e.,  "1/3 state".
       *    Inputs:  F (flux), rho, v, P
       *    Output:  third_state
       */
      SAMRAI_F77_FUNC(onethirdstate, ONETHIRDSTATE) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         flux->getPointer(0),
         flux->getPointer(1),
         flux->getPointer(2),
         third_state.getPointer());

      /*
       *    Compute fluxes using 1/3 state traces, in the two directions OTHER
       *    than idir.
       *    Inputs:  third_state, rho, v, P
       *    Output:  temp_flux (only directions other than idir are modified)
       */
      SAMRAI_F77_FUNC(fluxthird, FLUXTHIRD) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_gamma,
         d_riemann_solve_int,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         third_state.getPointer(),
         temp_flux.getPointer(0),
         temp_flux.getPointer(1),
         temp_flux.getPointer(2));
      /*
       *    Compute transverse corrections for the traces in the two
       *    directions (OTHER than idir) using the flux differences
       *    computed in those directions.
       *    Inputs:  temp_flux, rho, v, P
       *    Output:  w^L, w^R (traced_left/right)
       */
      SAMRAI_F77_FUNC(fluxcorrecjt, FLUXCORRECJT) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer(),
         temp_flux.getPointer(0),
         temp_flux.getPointer(1),
         temp_flux.getPointer(2),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_left.getPointer(2),
         traced_right.getPointer(0),
         traced_right.getPointer(1),
         traced_right.getPointer(2));

   }  // loop over directions...

   boundaryReset(patch, traced_left, traced_right);

   /*
    *  Final flux calculation using corrected trace states.
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 0, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_gamma,
      d_riemann_solve_int,
      density->getPointer(),
      velocity->getPointer(),
      pressure->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

}

/*
 *************************************************************************
 *
 * Update Euler solution variables by performing a conservative
 * difference with the fluxes calculated in computeFluxesOnPatch().
 * Although, "primitive" variables are maintained (i.e., density,
 * velocity, pressure), "conserved" variables (i.e., density,
 * momentum, total energy) are conserved.
 *
 *************************************************************************
 */

void Euler::conservativeDifferenceOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt,
   bool at_syncronization)
{
   NULL_USE(time);
   NULL_USE(dt);
   NULL_USE(at_syncronization);

   t_conservdiff->start();

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, getDataContext())));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(density->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(velocity->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(pressure->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(consdiff2d, CONSDIFF2D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         dx,
         flux->getPointer(0),
         flux->getPointer(1),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer());
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(consdiff3d, CONSDIFF3D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         ifirst(2), ilast(2), dx,
         flux->getPointer(0),
         flux->getPointer(1),
         flux->getPointer(2),
         d_gamma,
         density->getPointer(),
         velocity->getPointer(),
         pressure->getPointer());
   }

   t_conservdiff->stop();

}

/*
 *************************************************************************
 *
 * Reset physical boundary values for special cases, such as those
 * involving reflective boundary conditions and when the "STEP"
 * problem is run.
 *
 *************************************************************************
 */

void Euler::boundaryReset(
   hier::Patch& patch,
   pdat::FaceData<double>& traced_left,
   pdat::FaceData<double>& traced_right) const
{
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   int idir;
   bool bdry_cell = true;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   hier::BoxContainer domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes,
      patch_geom->getRatio(),
      hier::BlockId::zero());
   const double* dx = patch_geom->getDx();
   const double* xpatchhi = patch_geom->getXUpper();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   pdat::CellIndex icell(ifirst);
   hier::BoxContainer bdrybox;
   hier::Index ibfirst = ifirst;
   hier::Index iblast = ilast;
   int bdry_case = 0;

   for (idir = 0; idir < d_dim.getValue(); ++idir) {
      ibfirst(idir) = ifirst(idir) - 1;
      iblast(idir) = ifirst(idir) - 1;
      bdrybox.pushBack(hier::Box(ibfirst, iblast, hier::BlockId(0)));

      ibfirst(idir) = ilast(idir) + 1;
      iblast(idir) = ilast(idir) + 1;
      bdrybox.pushBack(hier::Box(ibfirst, iblast, hier::BlockId(0)));
   }

   hier::BoxContainer::iterator ib = bdrybox.begin();
   for (idir = 0; idir < d_dim.getValue(); ++idir) {
      int bside = 2 * idir;
      if (d_dim == tbox::Dimension(2)) {
         bdry_case = d_master_bdry_edge_conds[bside];
      }
      if (d_dim == tbox::Dimension(3)) {
         bdry_case = d_master_bdry_face_conds[bside];
      }
      if (bdry_case == BdryCond::REFLECT) {
         pdat::CellIterator icend(pdat::CellGeometry::end(*ib));
         for (pdat::CellIterator ic(pdat::CellGeometry::begin(*ib));
              ic != icend; ++ic) {
            for (hier::BoxContainer::iterator domain_boxes_itr =
                    domain_boxes.begin();
                 domain_boxes_itr != domain_boxes.end();
                 ++domain_boxes_itr) {
               if (domain_boxes_itr->contains(*ic))
                  bdry_cell = false;
            }
            if (bdry_cell) {
               pdat::FaceIndex sidein = pdat::FaceIndex(*ic, idir, 1);
               (traced_left)(sidein, 0) = (traced_right)(sidein, 0);
            }
         }
      }
      ++ib;

      int bnode = 2 * idir + 1;
      if (d_dim == tbox::Dimension(2)) {
         bdry_case = d_master_bdry_edge_conds[bnode];
      }
      if (d_dim == tbox::Dimension(3)) {
         bdry_case = d_master_bdry_face_conds[bnode];
      }
// BEGIN SIMPLE-MINDED FIX FOR STEP PROBLEM
      if ((d_data_problem == "STEP") && (bnode == 1) &&
          (tbox::MathUtilities<double>::Abs(xpatchhi[0] - xdomainhi[0]) < dx[0])) {
         bdry_case = BdryCond::FLOW;
      }
// END SIMPLE-MINDED FIX FOR STEP PROBLEM
      if (bdry_case == BdryCond::REFLECT) {
         pdat::CellIterator icend(pdat::CellGeometry::end(*ib));
         for (pdat::CellIterator ic(pdat::CellGeometry::begin(*ib));
              ic != icend; ++ic) {
            for (hier::BoxContainer::iterator domain_boxes_itr =
                    domain_boxes.begin();
                 domain_boxes_itr != domain_boxes.end();
                 ++domain_boxes_itr) {
               if (domain_boxes_itr->contains(*ic))
                  bdry_cell = false;
            }
            if (bdry_cell) {
               pdat::FaceIndex sidein = pdat::FaceIndex(*ic, idir, 0);
               (traced_right)(sidein, 0) = (traced_left)(sidein, 0);
            }
         }
      }
      ++ib;
   }
}

/*
 *************************************************************************
 *
 * Refine velocity and pressure by conservatively refining
 * momentum and total energy.
 *
 *************************************************************************
 */

void Euler::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{

   std::shared_ptr<pdat::CellData<double> > cdensity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvelocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cpressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_pressure, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > fdensity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvelocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fpressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_pressure, getDataContext())));

   TBOX_ASSERT(cdensity);
   TBOX_ASSERT(cvelocity);
   TBOX_ASSERT(cpressure);
   TBOX_ASSERT(fdensity);
   TBOX_ASSERT(fvelocity);
   TBOX_ASSERT(fpressure);

#ifdef DEBUG_CHECK_ASSERTIONS
   hier::IntVector gccheck = cdensity->getGhostCellWidth();
   TBOX_ASSERT(cvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   TBOX_ASSERT(fvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Box cgbox(cdensity->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdensity->getGhostBox().lower();
   const hier::Index fihi = fdensity->getGhostBox().upper();

   const std::shared_ptr<geom::CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<geom::CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));
   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   const hier::IntVector cons_ghosts(d_dim, 1);
   pdat::CellData<double> conserved(coarse_box, 1, cons_ghosts);

   const hier::IntVector tmp_ghosts(d_dim, 0);

   double* diff0 = new double[coarse_box.numberCells(0) + 1];
   pdat::CellData<double> slope0(coarse_box, 1, tmp_ghosts);

   double* diff1 = new double[coarse_box.numberCells(1) + 1];
   pdat::CellData<double> slope1(coarse_box, 1, tmp_ghosts);

   double* diff2 = d_dim ==
      tbox::Dimension(3) ? new double[coarse_box.numberCells(2) + 1] : 0;
   pdat::CellData<double> slope2(coarse_box, 1, tmp_ghosts);

   if (d_dim == tbox::Dimension(2)) {
      pdat::CellData<double> flat0(coarse_box, 1, tmp_ghosts);
      pdat::CellData<double> flat1(coarse_box, 1, tmp_ghosts);
      int mc = cihi(0) - cilo(0) + 1;
      mc = tbox::MathUtilities<int>::Max(mc, cihi(1) - cilo(1) + 1);
      double* tflat = new double[mc];
      double* tflat2 = new double[mc];
      double* tsound = new double[mc];
      double* tdensc = new double[mc];
      double* tpresc = new double[mc];
      double* tvelc = new double[mc];
      SAMRAI_F77_FUNC(conservlinint2d, CONSERVLININT2D) (ifirstc(0), ifirstc(1),
         ilastc(0), ilastc(1),                                                              /* input */
         ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
         cilo(0), cilo(1), cihi(0), cihi(1),
         filo(0), filo(1), fihi(0), fihi(1),
         &ratio[0],
         cgeom->getDx(),
         fgeom->getDx(),
         d_gamma,
         cdensity->getPointer(),
         fdensity->getPointer(),
         cvelocity->getPointer(),
         cpressure->getPointer(),
         fvelocity->getPointer(),                                 /* output */
         fpressure->getPointer(),
         conserved.getPointer(),                                 /* temporaries */
         tflat, tflat2, tsound, mc,
         tdensc, tpresc, tvelc,
         flat0.getPointer(),
         flat1.getPointer(),
         diff0, slope0.getPointer(),
         diff1, slope1.getPointer());
      delete[] tflat;
      delete[] tflat2;
      delete[] tsound;
      delete[] tdensc;
      delete[] tpresc;
      delete[] tvelc;
   } else if (d_dim == tbox::Dimension(3)) {
      pdat::CellData<double> flat0(coarse_box, 1, tmp_ghosts);
      pdat::CellData<double> flat1(coarse_box, 1, tmp_ghosts);
      pdat::CellData<double> flat2(coarse_box, 1, tmp_ghosts);
      int mc = cihi(0) - cilo(0) + 1;
      mc = tbox::MathUtilities<int>::Max(mc, cihi(1) - cilo(1) + 1);
      mc = tbox::MathUtilities<int>::Max(mc, cihi(2) - cilo(2) + 1);
      double* tflat = new double[mc];
      double* tflat2 = new double[mc];
      double* tsound = new double[mc];
      double* tdensc = new double[mc];
      double* tpresc = new double[mc];
      double* tvelc = new double[mc];
      SAMRAI_F77_FUNC(conservlinint3d, CONSERVLININT3D) (ifirstc(0), ifirstc(1),
         ifirstc(2),                                                                      /* input */
         ilastc(0), ilastc(1), ilastc(2),
         ifirstf(0), ifirstf(1), ifirstf(2),
         ilastf(0), ilastf(1), ilastf(2),
         cilo(0), cilo(1), cilo(2), cihi(0), cihi(1), cihi(2),
         filo(0), filo(1), filo(2), fihi(0), fihi(1), fihi(2),
         &ratio[0],
         cgeom->getDx(),
         fgeom->getDx(),
         d_gamma,
         cdensity->getPointer(),
         fdensity->getPointer(),
         cvelocity->getPointer(),
         cpressure->getPointer(),
         fvelocity->getPointer(),                                 /* output */
         fpressure->getPointer(),
         conserved.getPointer(),                                 /* temporaries */
         tflat, tflat2, tsound, mc,
         tdensc, tpresc, tvelc,
         flat0.getPointer(),
         flat1.getPointer(),
         flat2.getPointer(),
         diff0, slope0.getPointer(),
         diff1, slope1.getPointer(),
         diff2, slope2.getPointer());
      delete[] tflat;
      delete[] tflat2;
      delete[] tsound;
      delete[] tdensc;
      delete[] tpresc;
      delete[] tvelc;
   }

   delete[] diff0;
   delete[] diff1;
   if (d_dim == tbox::Dimension(3)) {
      delete[] diff2;
   }

}

/*
 *************************************************************************
 *
 * Coarsen velocity and pressure by conservatively coarsening
 * momentum and total energy.
 *
 *************************************************************************
 */

void Euler::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{

   std::shared_ptr<pdat::CellData<double> > fdensity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvelocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fpressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_pressure, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cdensity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvelocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cpressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_pressure, getDataContext())));

   TBOX_ASSERT(cdensity);
   TBOX_ASSERT(cvelocity);
   TBOX_ASSERT(cpressure);
   TBOX_ASSERT(fdensity);
   TBOX_ASSERT(fvelocity);
   TBOX_ASSERT(fpressure);

#ifdef DEBUG_CHECK_ASSERTIONS
   hier::IntVector gccheck = cdensity->getGhostCellWidth();
   TBOX_ASSERT(cvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(cpressure->getGhostCellWidth() == gccheck);

   gccheck = fdensity->getGhostCellWidth();
   TBOX_ASSERT(fvelocity->getGhostCellWidth() == gccheck);
   TBOX_ASSERT(fpressure->getGhostCellWidth() == gccheck);
#endif

   const hier::Index filo = fdensity->getGhostBox().lower();
   const hier::Index fihi = fdensity->getGhostBox().upper();
   const hier::Index cilo = cdensity->getGhostBox().lower();
   const hier::Index cihi = cdensity->getGhostBox().upper();

   const std::shared_ptr<geom::CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));
   const std::shared_ptr<geom::CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   TBOX_ASSERT(fgeom);
   TBOX_ASSERT(cgeom);

   const hier::Box fine_box = hier::Box::refine(coarse_box, ratio);
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   const hier::IntVector cons_ghosts(d_dim, 0);
   pdat::CellData<double> conserved(fine_box, 1, cons_ghosts);

   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(conservavg2d, CONSERVAVG2D) (ifirstf(0), ifirstf(1), ilastf(0),
         ilastf(1),                                                                   /* input */
         ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
         filo(0), filo(1), fihi(0), fihi(1),
         cilo(0), cilo(1), cihi(0), cihi(1),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         d_gamma,
         fdensity->getPointer(),
         cdensity->getPointer(),
         fvelocity->getPointer(),
         fpressure->getPointer(),
         cvelocity->getPointer(),                               /* output */
         cpressure->getPointer(),
         conserved.getPointer());                              /* temporary */
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(conservavg3d, CONSERVAVG3D) (ifirstf(0), ifirstf(1), ifirstf(2),       /* input */
         ilastf(0), ilastf(1), ilastf(2),
         ifirstc(0), ifirstc(1), ifirstc(2),
         ilastc(0), ilastc(1), ilastc(2),
         filo(0), filo(1), filo(2), fihi(0), fihi(1), fihi(2),
         cilo(0), cilo(1), cilo(2), cihi(0), cihi(1), cihi(2),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         d_gamma,
         fdensity->getPointer(),
         cdensity->getPointer(),
         fvelocity->getPointer(),
         fpressure->getPointer(),
         cvelocity->getPointer(),                               /* output */
         cpressure->getPointer(),
         conserved.getPointer());
   }

}

/*
 *************************************************************************
 *
 * Set the data in ghost cells corresponding to physical boundary
 * conditions.  Note that boundary geometry configuration information
 * (i.e., faces, edges, and nodes) is obtained from the patch geometry
 * object owned by the patch.
 *
 *************************************************************************
 */

void Euler::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);
   t_setphysbcs->start();

   std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, getDataContext())));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);
#ifdef DEBUG_CHECK_ASSERTIONS
   hier::IntVector ghost_cells = density->getGhostCellWidth();
   TBOX_ASSERT(velocity->getGhostCellWidth() == ghost_cells);
   TBOX_ASSERT(pressure->getGhostCellWidth() == ghost_cells);
#endif

   if (d_dim == tbox::Dimension(2)) {

      /*
       * Set boundary conditions for cells corresponding to patch edges.
       *
       * Note: We apply a simple-minded adjustment for the "STEP" problem
       *       so that the right edge of the domain gets (out)FLOW conditions
       *       whereas the right edge at the step gets REFLECT condtions (from input),
       */
      std::vector<int> tmp_edge_scalar_bcond(NUM_2D_EDGES);
      std::vector<int> tmp_edge_vector_bcond(NUM_2D_EDGES);
      for (int i = 0; i < NUM_2D_EDGES; ++i) {
         tmp_edge_scalar_bcond[i] = d_scalar_bdry_edge_conds[i];
         tmp_edge_vector_bcond[i] = d_vector_bdry_edge_conds[i];
      }

      if (d_data_problem == "STEP") {

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch.getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();
         const double* xpatchhi = patch_geom->getXUpper();
         const double* xdomainhi = d_grid_geometry->getXUpper();

         if (tbox::MathUtilities<double>::Abs(xpatchhi[0] - xdomainhi[0]) <
             dx[0]) {
            tmp_edge_scalar_bcond[BdryLoc::XHI] = BdryCond::FLOW;
            tmp_edge_vector_bcond[BdryLoc::XHI] = BdryCond::FLOW;
         }

      }

      appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("density", density,
         patch,
         ghost_width_to_fill,
         tmp_edge_scalar_bcond,
         d_bdry_edge_density);
      appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("velocity", velocity,
         patch,
         ghost_width_to_fill,
         tmp_edge_vector_bcond,
         d_bdry_edge_velocity);
      appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("pressure", pressure,
         patch,
         ghost_width_to_fill,
         tmp_edge_scalar_bcond,
         d_bdry_edge_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill,
         tmp_edge_scalar_bcond, tmp_edge_vector_bcond);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("density", density,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_edge_density);
      appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("velocity", velocity,
         patch,
         ghost_width_to_fill,
         d_vector_bdry_node_conds,
         d_bdry_edge_velocity);
      appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("pressure", pressure,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_edge_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds, d_vector_bdry_node_conds);
#endif
#endif

   }

   if (d_dim == tbox::Dimension(3)) {

      /*
       *  Set boundary conditions for cells corresponding to patch faces.
       */

      appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("density", density,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_face_conds,
         d_bdry_face_density);
      appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("velocity", velocity,
         patch,
         ghost_width_to_fill,
         d_vector_bdry_face_conds,
         d_bdry_face_velocity);
      appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("pressure", pressure,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_face_conds,
         d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_face_conds, d_vector_bdry_face_conds);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch edges.
       */

      appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("density", density,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_face_density);
      appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("velocity", velocity,
         patch,
         ghost_width_to_fill,
         d_vector_bdry_edge_conds,
         d_bdry_face_velocity);
      appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("pressure", pressure,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds, d_vector_bdry_edge_conds);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("density", density,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_face_density);
      appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("velocity", velocity,
         patch,
         ghost_width_to_fill,
         d_vector_bdry_node_conds,
         d_bdry_face_velocity);
      appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("pressure", pressure,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_face_pressure);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds, d_scalar_bdry_node_conds);
#endif
#endif

   }

   t_setphysbcs->stop();
}

/*
 *************************************************************************
 *
 * Tag cells for refinement using gradient detector.  Tagging criteria
 * defined in input.
 *
 *************************************************************************
 */

void Euler::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{
   t_taggradient->start();

   const int error_level_number = patch.getPatchLevelNumber();

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<int> > tags(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(tag_indx)));
   TBOX_ASSERT(tags);

   hier::Box pbox = patch.getBox();
   hier::BoxContainer domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes,
      patch_geom->getRatio(),
      hier::BlockId::zero());
   /*
    * Construct domain bounding box
    */
   hier::Box domain(d_dim);
   for (hier::BoxContainer::iterator i = domain_boxes.begin();
        i != domain_boxes.end(); ++i) {
      domain += *i;
   }

   const hier::Index domfirst = domain.lower();
   const hier::Index domlast = domain.upper();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   /*
    * Create a set of temporary tags and set to untagged value.
    */

   std::shared_ptr<pdat::CellData<int> > temp_tags(
      new pdat::CellData<int>(
         pbox,
         1,
         d_nghosts));
   temp_tags->fillAll(FALSE);
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

   if (d_dim == tbox::Dimension(2)) {
      /*
       * Problem specific criteria for step case.
       */
      if (initial_error && d_data_problem == "STEP") {
         if (error_level_number < 2) {
            hier::Box tagbox(hier::Index(9, 0), hier::Index(9, 3), hier::BlockId(0));
            if (error_level_number == 1) {
               tagbox.refine(hier::IntVector(d_dim, 2));
#if defined(HAVE_RAJA)
               tbox::parallel_synchronize();
#endif
            }
            hier::Box ibox = pbox * tagbox;

            pdat::CellIterator itcend(pdat::CellGeometry::end(ibox));
            for (pdat::CellIterator itc(pdat::CellGeometry::begin(ibox));
                 itc != itcend; ++itc) {
               (*temp_tags)(*itc, 0) = TRUE;
            }
         }
      }
   }

   /*
    * Possible tagging criteria includes
    *    DENSITY_DEVIATION, DENSITY_GRADIENT, DENSITY_SHOCK
    *    PRESSURE_DEVIATION, PRESSURE_GRADIENT, PRESSURE_SHOCK
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0;
        ncrit < static_cast<int>(d_refinement_criteria.size()); ++ncrit) {

      std::string ref = d_refinement_criteria[ncrit];
      std::shared_ptr<pdat::CellData<double> > var;
      int size = 0;
      double tol = 0.;
      double onset = 0.;
      double dev = 0.;
      bool time_allowed = false;

      if (ref == "DENSITY_DEVIATION") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_density, getDataContext()));
         size = static_cast<int>(d_density_dev_tol.size());
         tol = ((error_level_number < size)
                ? d_density_dev_tol[error_level_number]
                : d_density_dev_tol[size - 1]);
         size = static_cast<int>(d_density_dev.size());
         dev = ((error_level_number < size)
                ? d_density_dev[error_level_number]
                : d_density_dev[size - 1]);
         size = static_cast<int>(d_density_dev_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_density_dev_time_min[error_level_number]
                            : d_density_dev_time_min[size - 1]);
         size = static_cast<int>(d_density_dev_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_density_dev_time_max[error_level_number]
                            : d_density_dev_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "DENSITY_GRADIENT") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_density, getDataContext()));
         size = static_cast<int>(d_density_grad_tol.size());
         tol = ((error_level_number < size)
                ? d_density_grad_tol[error_level_number]
                : d_density_grad_tol[size - 1]);
         size = static_cast<int>(d_density_grad_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_density_grad_time_min[error_level_number]
                            : d_density_grad_time_min[size - 1]);
         size = static_cast<int>(d_density_grad_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_density_grad_time_max[error_level_number]
                            : d_density_grad_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "DENSITY_SHOCK") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_density, getDataContext()));
         size = static_cast<int>(d_density_shock_tol.size());
         tol = ((error_level_number < size)
                ? d_density_shock_tol[error_level_number]
                : d_density_shock_tol[size - 1]);
         size = static_cast<int>(d_density_shock_onset.size());
         onset = ((error_level_number < size)
                  ? d_density_shock_onset[error_level_number]
                  : d_density_shock_onset[size - 1]);
         size = static_cast<int>(d_density_shock_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_density_shock_time_min[error_level_number]
                            : d_density_shock_time_min[size - 1]);
         size = static_cast<int>(d_density_shock_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_density_shock_time_max[error_level_number]
                            : d_density_shock_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "PRESSURE_DEVIATION") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_pressure, getDataContext()));
         size = static_cast<int>(d_pressure_dev_tol.size());
         tol = ((error_level_number < size)
                ? d_pressure_dev_tol[error_level_number]
                : d_pressure_dev_tol[size - 1]);
         size = static_cast<int>(d_pressure_dev.size());
         dev = ((error_level_number < size)
                ? d_pressure_dev[error_level_number]
                : d_pressure_dev[size - 1]);
         size = static_cast<int>(d_pressure_dev_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_pressure_dev_time_min[error_level_number]
                            : d_pressure_dev_time_min[size - 1]);
         size = static_cast<int>(d_pressure_dev_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_pressure_dev_time_max[error_level_number]
                            : d_pressure_dev_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "PRESSURE_GRADIENT") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_pressure, getDataContext()));
         size = static_cast<int>(d_pressure_grad_tol.size());
         tol = ((error_level_number < size)
                ? d_pressure_grad_tol[error_level_number]
                : d_pressure_grad_tol[size - 1]);
         size = static_cast<int>(d_pressure_grad_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_pressure_grad_time_min[error_level_number]
                            : d_pressure_grad_time_min[size - 1]);
         size = static_cast<int>(d_pressure_grad_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_pressure_grad_time_max[error_level_number]
                            : d_pressure_grad_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "PRESSURE_SHOCK") {
         var = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_pressure, getDataContext()));
         size = static_cast<int>(d_pressure_shock_tol.size());
         tol = ((error_level_number < size)
                ? d_pressure_shock_tol[error_level_number]
                : d_pressure_shock_tol[size - 1]);
         size = static_cast<int>(d_pressure_shock_onset.size());
         onset = ((error_level_number < size)
                  ? d_pressure_shock_onset[error_level_number]
                  : d_pressure_shock_onset[size - 1]);
         size = static_cast<int>(d_pressure_shock_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_pressure_shock_time_min[error_level_number]
                            : d_pressure_shock_time_min[size - 1]);
         size = static_cast<int>(d_pressure_shock_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_pressure_shock_time_max[error_level_number]
                            : d_pressure_shock_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (time_allowed) {

         TBOX_ASSERT(var);

         hier::IntVector vghost = var->getGhostCellWidth();
         hier::IntVector tagghost = tags->getGhostCellWidth();

         if (ref == "DENSITY_DEVIATION" || ref == "PRESSURE_DEVIATION") {

            /*
             * Check for tags that have already been set in a previous
             * step.  Do NOT consider values tagged with value
             * RICHARDSON_NEWLY_TAGGED since these were set most recently
             * by Richardson extrapolation.
             */
            pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
                 ic != icend; ++ic) {
               double locden = tol;
               int tag_val = (*tags)(*ic, 0);
               if (tag_val) {
                  if (tag_val != RICHARDSON_NEWLY_TAGGED) {
                     locden *= 0.75;
                  }
               }
               if (tbox::MathUtilities<double>::Abs((*var)(*ic) - dev) >
                   locden) {
                  (*temp_tags)(*ic, 0) = TRUE;
               }
            }
         } else if (ref == "DENSITY_GRADIENT" || ref == "PRESSURE_GRADIENT") {
            if (d_dim == tbox::Dimension(2)) {
               SAMRAI_F77_FUNC(detectgrad2d, DETECTGRAD2D) (ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  vghost(0), tagghost(0), d_nghosts(0),
                  vghost(1), tagghost(1), d_nghosts(1),
                  dx,
                  tol,
                  TRUE, FALSE,
                  var->getPointer(),
                  tags->getPointer(), temp_tags->getPointer());
            } else if (d_dim == tbox::Dimension(3)) {
               SAMRAI_F77_FUNC(detectgrad3d, DETECTGRAD3D) (ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  vghost(0), tagghost(0), d_nghosts(0),
                  vghost(1), tagghost(1), d_nghosts(1),
                  vghost(2), tagghost(2), d_nghosts(2),
                  dx,
                  tol,
                  TRUE, FALSE,
                  var->getPointer(),
                  tags->getPointer(), temp_tags->getPointer());
            }
         } else if (ref == "DENSITY_SHOCK" || ref == "PRESSURE_SHOCK") {
            if (d_dim == tbox::Dimension(2)) {
               SAMRAI_F77_FUNC(detectshock2d, DETECTSHOCK2D) (ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  vghost(0), tagghost(0), d_nghosts(0),
                  vghost(1), tagghost(1), d_nghosts(1),
                  dx,
                  tol,
                  onset,
                  TRUE, FALSE,
                  var->getPointer(),
                  tags->getPointer(), temp_tags->getPointer());
            } else if (d_dim == tbox::Dimension(3)) {
               SAMRAI_F77_FUNC(detectshock3d, DETECTSHOCK3D) (ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  vghost(0), tagghost(0), d_nghosts(0),
                  vghost(1), tagghost(1), d_nghosts(1),
                  vghost(2), tagghost(2), d_nghosts(2),
                  dx,
                  tol,
                  onset,
                  TRUE, FALSE,
                  var->getPointer(),
                  tags->getPointer(), temp_tags->getPointer());
            }
         }

      }  // if time_allowed

   }  // loop over criteria

   /*
    * Adjust temp_tags from those tags set in Richardson extrapolation.
    */
   if (uses_richardson_extrapolation_too) {
      pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
      for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
           ic != icend; ++ic) {
         if ((*tags)(*ic, 0) == RICHARDSON_ALREADY_TAGGED ||
             (*tags)(*ic, 0) == RICHARDSON_NEWLY_TAGGED) {
            (*temp_tags)(*ic, 0) = TRUE;
         }
      }
   }

   /*
    * Update tags
    */
   pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
   for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
        ic != icend; ++ic) {
      (*tags)(*ic, 0) = (*temp_tags)(*ic, 0);
   }

   t_taggradient->stop();
}

/*
 *************************************************************************
 *
 * Tag cells for refinement using Richardson extrapolation.  Criteria
 * defined in input.
 *
 *************************************************************************
 */

void Euler::tagRichardsonExtrapolationCells(
   hier::Patch& patch,
   const int error_level_number,
   const std::shared_ptr<hier::VariableContext>& coarsened_fine,
   const std::shared_ptr<hier::VariableContext>& advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   NULL_USE(initial_error);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* xdomainlo = d_grid_geometry->getXLower();
   const double* xdomainhi = d_grid_geometry->getXUpper();

   hier::Box pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<int> > tags(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(tag_index)));
   TBOX_ASSERT(tags);

   /*
    * Possible tagging criteria includes
    *    DENSITY_RICHARDSON, PRESSURE_RICHARDSON
    * The criteria is specified over a time interval.
    *
    * Loop over criteria provided and check to make sure we are in the
    * specified time interval.  If so, apply appropriate tagging for
    * the level.
    */
   for (int ncrit = 0;
        ncrit < static_cast<int>(d_refinement_criteria.size()); ++ncrit) {

      std::string ref = d_refinement_criteria[ncrit];
      std::shared_ptr<pdat::CellData<double> > coarsened_fine_var;
      std::shared_ptr<pdat::CellData<double> > advanced_coarse_var;
      int size = 0;
      double tol = 0.;
      bool time_allowed = false;

      if (ref == "DENSITY_RICHARDSON") {
         coarsened_fine_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_density, coarsened_fine));
         advanced_coarse_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_density, advanced_coarse));
         size = static_cast<int>(d_density_rich_tol.size());
         tol = ((error_level_number < size)
                ? d_density_rich_tol[error_level_number]
                : d_density_rich_tol[size - 1]);
         size = static_cast<int>(d_density_rich_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_density_rich_time_min[error_level_number]
                            : d_density_rich_time_min[size - 1]);
         size = static_cast<int>(d_density_rich_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_density_rich_time_max[error_level_number]
                            : d_density_rich_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      } else if (ref == "PRESSURE_RICHARDSON") {
         coarsened_fine_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_pressure, coarsened_fine));
         advanced_coarse_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_pressure, advanced_coarse));
         size = static_cast<int>(d_pressure_rich_tol.size());
         tol = ((error_level_number < size)
                ? d_pressure_rich_tol[error_level_number]
                : d_pressure_rich_tol[size - 1]);
         size = static_cast<int>(d_pressure_rich_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_pressure_rich_time_min[error_level_number]
                            : d_pressure_rich_time_min[size - 1]);
         size = static_cast<int>(d_pressure_rich_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_pressure_rich_time_max[error_level_number]
                            : d_pressure_rich_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);
      }

      if (time_allowed) {

         TBOX_ASSERT(coarsened_fine_var);
         TBOX_ASSERT(advanced_coarse_var);

         if (ref == "DENSITY_RICHARDSON" || ref == "PRESSURE_RICHARDSON") {

            /*
             * We tag wherever the global error > specified tolerance.
             * The estimated global error is the
             * local truncation error * the approximate number of steps
             * used in the simulation.  Approximate the number of steps as:
             *
             *       steps = L / (s*deltat)
             * where
             *       L = length of problem domain
             *       s = wave speed
             *       delta t = timestep on current level
             *
             * Compute max wave speed from delta t.  This presumes that
             * deltat was computed as deltat = dx/s_max.  We have deltat
             * and dx, so back out s_max from this.
             */

            const double* dx = patch_geom->getDx();
            double max_dx = 0.;
            double max_length = 0.;
            for (int idir = 0; idir < d_dim.getValue(); ++idir) {
               max_dx = tbox::MathUtilities<double>::Max(max_dx, dx[idir]);
               double length = xdomainhi[idir] - xdomainlo[idir];
               max_length =
                  tbox::MathUtilities<double>::Max(max_length, length);
            }
            double max_wave_speed = max_dx / deltat;
            double steps = max_length / (max_wave_speed * deltat);

            /*
             * Tag cells where |w_c - w_f| * (r^n -1) * steps
             *
             * where
             *       w_c = soln on coarse level (pressure_crse)
             *       w_f = soln on fine level (pressure_fine)
             *       r   = error coarsen ratio
             *       n   = spatial order of scheme (1st or 2nd depending
             *             on whether Godunov order is 1st or 2nd/4th)
             */
            int order = 1;
            if (d_godunov_order > 1) order = 2;
            double r = error_coarsen_ratio;
            double rnminus1 = pow(r, order) - 1;

            double diff = 0.;
            double error = 0.;

            pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
                 ic != icend; ++ic) {

               /*
                * Compute error norm
                */
               diff = (*advanced_coarse_var)(*ic, 0)
                  - (*coarsened_fine_var)(*ic, 0);
               error =
                  tbox::MathUtilities<double>::Abs(diff) * rnminus1 * steps;

               /*
                * Tag cell if error > prescribed threshold. Since we are
                * operating on the actual tag values (not temporary ones)
                * distinguish here tags that were previously set before
                * coming into this routine and those that are set here.
                *     RICHARDSON_ALREADY_TAGGED - tagged before coming
                *                                 into this method.
                *     RICHARDSON_NEWLY_TAGGED - newly tagged in this method
                */
               if (error > tol) {
                  if ((*tags)(*ic, 0)) {
                     (*tags)(*ic, 0) = RICHARDSON_ALREADY_TAGGED;
                  } else {
                     (*tags)(*ic, 0) = RICHARDSON_NEWLY_TAGGED;
                  }
               }
            }

         }

      } // time_allowed

   } // loop over refinement criteria

   /*
    * If we are NOT performing gradient detector (i.e. only
    * doing Richardson extrapolation) set tags marked in this method
    * to TRUE and all others false.  Otherwise, leave tags set to the
    * RICHARDSON_ALREADY_TAGGED and RICHARDSON_NEWLY_TAGGED as we may
    * use this information in the gradient detector.
    */
   if (!uses_gradient_detector_too) {
      pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
      for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
           ic != icend; ++ic) {
         if ((*tags)(*ic, 0) == RICHARDSON_ALREADY_TAGGED ||
             (*tags)(*ic, 0) == RICHARDSON_NEWLY_TAGGED) {
            (*tags)(*ic, 0) = TRUE;
         } else {
            (*tags)(*ic, 0) = FALSE;
         }
      }
   }

}

/*
 *************************************************************************
 *
 * Register VisIt data writer to write data to plot files that may
 * be postprocessed by the VisIt tool.
 *
 *************************************************************************
 */

#ifdef HAVE_HDF5
void Euler::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);
   d_visit_writer = viz_writer;
}
#endif

/*
 *************************************************************************
 *
 * Pack "total energy" and "momentum" (derived Vis plot quantities)
 * for the patch into a double precision buffer.
 *
 *************************************************************************
 */

bool Euler::packDerivedDataIntoDoubleBuffer(
   double* dbuffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(simulation_time);

   TBOX_ASSERT((region * patch.getBox()).isSpatiallyEqual(region));

   bool data_on_patch = FALSE;

   std::shared_ptr<pdat::CellData<double> > density(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_density, d_plot_context)));
   std::shared_ptr<pdat::CellData<double> > velocity(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_velocity, d_plot_context)));
   std::shared_ptr<pdat::CellData<double> > pressure(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_pressure, d_plot_context)));

   TBOX_ASSERT(density);
   TBOX_ASSERT(velocity);
   TBOX_ASSERT(pressure);
   TBOX_ASSERT(density->getGhostBox().isSpatiallyEqual(patch.getBox()));
   TBOX_ASSERT(velocity->getGhostBox().isSpatiallyEqual(patch.getBox()));
   TBOX_ASSERT(pressure->getGhostBox().isSpatiallyEqual(patch.getBox()));

   const hier::Box& data_box = density->getGhostBox();
   const int box_w0 = region.numberCells(0);
   const int dat_w0 = data_box.numberCells(0);
   const int box_w1 = region.numberCells(1);
   const int dat_w1 = d_dim >
      tbox::Dimension(2) ? data_box.numberCells(1) : tbox::MathUtilities<int>::getMax();
   const int box_w2 = d_dim >
      tbox::Dimension(2) ? region.numberCells(2) : tbox::MathUtilities<int>::getMax();

   if (variable_name == "Total Energy") {
      const double * const dens = density->getPointer();
      const double * const xvel = velocity->getPointer(0);
      const double * const yvel = velocity->getPointer(1);
      const double * const zvel = d_dim > tbox::Dimension(2) ? velocity->getPointer(2) : 0;
      const double * const pres = pressure->getPointer();

      double valinv = 1.0 / (d_gamma - 1.0);
      int buf_b1 = 0;
      size_t dat_b2 = data_box.offset(region.lower());

      if (d_dim > tbox::Dimension(2)) {
         for (int i2 = 0; i2 < box_w2; ++i2) {
            size_t dat_b1 = dat_b2;
            for (int i1 = 0; i1 < box_w1; ++i1) {
               for (int i0 = 0; i0 < box_w0; ++i0) {
                  size_t dat_indx = dat_b1 + i0;
                  double v2norm = pow(xvel[dat_indx], 2.0)
                     + pow(yvel[dat_indx], 2.0)
                     + pow(zvel[dat_indx], 2.0)
                  ;
                  double rho = dens[dat_indx];
                  double int_energy = 0.0;
                  if (rho > 0.0) {
                     int_energy = valinv * pres[dat_indx] / dens[dat_indx];
                  }
                  dbuffer[buf_b1 + i0] =
                     dens[dat_indx] * (0.5 * v2norm + int_energy);
               }
               dat_b1 += dat_w0;
               buf_b1 += box_w0;
            }
            dat_b2 += dat_w1 * dat_w0;
         }
      }

      if (d_dim == tbox::Dimension(2)) {
         size_t dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; ++i1) {
            for (int i0 = 0; i0 < box_w0; ++i0) {
               size_t dat_indx = dat_b1 + i0;
               double v2norm = pow(xvel[dat_indx], 2.0)
                  + pow(yvel[dat_indx], 2.0)
               ;
               double rho = dens[dat_indx];
               double int_energy = 0.0;
               if (rho > 0.0) {
                  int_energy = valinv * pres[dat_indx] / dens[dat_indx];
               }
               dbuffer[buf_b1 + i0] =
                  dens[dat_indx] * (0.5 * v2norm + int_energy);
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
      }

      data_on_patch = TRUE;

   } else if (variable_name == "Momentum") {
      TBOX_ASSERT(depth_id < d_dim.getValue());

      const double * const dens = density->getPointer();
      const double * const vel = velocity->getPointer(depth_id);
      int buf_b1 = 0;
      size_t dat_b2 = data_box.offset(region.lower());

      if (d_dim == tbox::Dimension(2)) {
         size_t dat_b1 = dat_b2;
         for (int i1 = 0; i1 < box_w1; ++i1) {
            for (int i0 = 0; i0 < box_w0; ++i0) {
               size_t dat_indx = dat_b1 + i0;
               dbuffer[buf_b1 + i0] = dens[dat_indx] * vel[dat_indx];
            }
            dat_b1 += dat_w0;
            buf_b1 += box_w0;
         }
      }
      if (d_dim == tbox::Dimension(3)) {
         for (int i2 = 0; i2 < box_w2; ++i2) {
            size_t dat_b1 = dat_b2;
            for (int i1 = 0; i1 < box_w1; ++i1) {
               for (int i0 = 0; i0 < box_w0; ++i0) {
                  size_t dat_indx = dat_b1 + i0;
                  dbuffer[buf_b1 + i0] = dens[dat_indx] * vel[dat_indx];
               }
               dat_b1 += dat_w0;
               buf_b1 += box_w0;
            }
            dat_b2 += dat_w1 * dat_w0;
         }
      }

      data_on_patch = TRUE;

   } else {
      TBOX_ERROR("Euler::packDerivedDataIntoDoubleBuffer()"
         << "\n    unknown variable_name " << variable_name << "\n");
   }

   return data_on_patch;

}

/*
 *************************************************************************
 *
 * Write 1d data intersection of patch and pencil box to file
 * with given name for plotting with Matlab.
 *
 *************************************************************************
 */

void Euler::writeData1dPencil(
   const std::shared_ptr<hier::Patch> patch,
   const hier::Box& pencil_box,
   const tbox::Dimension::dir_t idir,
   std::ostream& file)
{

   const hier::Box& patchbox = patch->getBox();
   const hier::Box box = pencil_box * patchbox;

   if (!box.empty()) {

      std::shared_ptr<pdat::CellData<double> > density(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_density, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > velocity(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_velocity, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > pressure(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(d_pressure, getDataContext())));

      TBOX_ASSERT(density);
      TBOX_ASSERT(velocity);
      TBOX_ASSERT(pressure);

      const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));
      TBOX_ASSERT(pgeom);
      const double* dx = pgeom->getDx();
      const double* xlo = pgeom->getXLower();

      const double cell_center = xlo[idir]
         + (double(box.lower(idir)
                   - patchbox.lower(idir))
            + 0.5) * dx[idir];

      double valinv = 1.0 / (d_gamma - 1.0);

      int ccount = 0;
      pdat::CellIterator icend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator ic(pdat::CellGeometry::begin(box));
           ic != icend; ++ic) {
         file << cell_center + ccount * dx[idir] << " ";
         ++ccount;

         double rho = (*density)(*ic, 0);
         double vel = (*velocity)(*ic, idir);
         double p = (*pressure)(*ic, 0);

         double mom = rho * vel;
         double eint = 0.0;
         if (rho > 0.0) {
            eint = valinv * (p / rho);
         }
         double etot = rho * (0.5 * vel * vel + eint);

         /*
          * Write out conserved quantities.
          */
         file << rho << " ";
         file << mom << " ";
         file << etot << " ";

         /*
          * Write out "primitive" quantities and internal energy.
          */
         file << p << " ";
         file << vel << " ";
         file << eint << " ";

         file << std::endl;
      }

   }

}

/*
 *************************************************************************
 *
 * Write all class data members to specified output stream.
 *
 *************************************************************************
 */

void Euler::printClassData(
   std::ostream& os) const
{
   int j, k;

   os << "\nEuler::printClassData..." << std::endl;
   os << "Euler: this = " << (Euler *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_grid_geometry = "
      << d_grid_geometry.get() << std::endl;

   os << "Parameters for physical problem ..." << std::endl;
   os << "   d_gamma = " << d_gamma << std::endl;

   os << "Numerical method description and ghost sizes..." << std::endl;
   os << "   d_riemann_solve = " << d_riemann_solve << std::endl;
   os << "   d_riemann_solve_int = " << d_riemann_solve_int << std::endl;
   os << "   d_godunov_order = " << d_godunov_order << std::endl;
   os << "   d_corner_transport = " << d_corner_transport << std::endl;
   os << "   d_nghosts = " << d_nghosts << std::endl;
   os << "   d_fluxghosts = " << d_fluxghosts << std::endl;

   os << "Problem description and initial data..." << std::endl;
   os << "   d_data_problem = " << d_data_problem << std::endl;
   os << "   d_data_problem_int = " << d_data_problem_int << std::endl;

   os << "       d_radius = " << d_radius << std::endl;
   os << "       d_center = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_center[j] << " ";
   os << std::endl;
   os << "       d_density_inside = " << d_density_inside << std::endl;
   os << "       d_velocity_inside = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_velocity_inside[j] << " ";
   os << std::endl;
   os << "       d_pressure_inside = " << d_pressure_inside << std::endl;
   os << "       d_density_outside = " << d_density_outside << std::endl;
   os << "       d_velocity_outside = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_velocity_outside[j] << " ";
   os << std::endl;
   os << "       d_pressure_outside = " << d_pressure_outside << std::endl;

   os << "       d_number_of_intervals = " << d_number_of_intervals << std::endl;
   os << "       d_front_position = ";
   for (k = 0; k < d_number_of_intervals - 1; ++k) {
      os << d_front_position[k] << "  ";
   }
   os << std::endl;
   os << "       d_interval_density = " << std::endl;
   for (k = 0; k < d_number_of_intervals; ++k) {
      os << "            " << d_interval_density[k] << std::endl;
   }
   os << "       d_interval_velocity = " << std::endl;
   for (k = 0; k < d_number_of_intervals; ++k) {
      os << "            ";
      for (j = 0; j < d_dim.getValue(); ++j) {
         os << d_interval_velocity[k * d_dim.getValue() + j] << "  ";
      }
      os << std::endl;
   }
   os << "       d_interval_pressure = " << std::endl;
   for (k = 0; k < d_number_of_intervals; ++k) {
      os << "            " << d_interval_pressure[k] << std::endl;
   }

   os << "   Boundary condition data " << std::endl;

   if (d_dim == tbox::Dimension(2)) {
      for (j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); ++j) {
         os << "\n       d_master_bdry_edge_conds[" << j << "] = "
            << d_master_bdry_edge_conds[j] << std::endl;
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << std::endl;
         os << "       d_vector_bdry_edge_conds[" << j << "] = "
            << d_vector_bdry_edge_conds[j] << std::endl;
         if (d_master_bdry_edge_conds[j] == BdryCond::DIRICHLET) {
            os << "         d_bdry_edge_density[" << j << "] = "
               << d_bdry_edge_density[j] << std::endl;
            os << "         d_bdry_edge_velocity[" << j << "] = "
               << d_bdry_edge_velocity[j * d_dim.getValue() + 0] << " , "
               << d_bdry_edge_velocity[j * d_dim.getValue() + 1] << std::endl;
            os << "         d_bdry_edge_pressure[" << j << "] = "
               << d_bdry_edge_pressure[j] << std::endl;
         }
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); ++j) {
         os << "\n       d_master_bdry_node_conds[" << j << "] = "
            << d_master_bdry_node_conds[j] << std::endl;
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << std::endl;
         os << "       d_vector_bdry_node_conds[" << j << "] = "
            << d_vector_bdry_node_conds[j] << std::endl;
         os << "       d_node_bdry_edge[" << j << "] = "
            << d_node_bdry_edge[j] << std::endl;
      }
   }
   if (d_dim == tbox::Dimension(3)) {
      for (j = 0; j < static_cast<int>(d_master_bdry_face_conds.size()); ++j) {
         os << "\n       d_master_bdry_face_conds[" << j << "] = "
            << d_master_bdry_face_conds[j] << std::endl;
         os << "       d_scalar_bdry_face_conds[" << j << "] = "
            << d_scalar_bdry_face_conds[j] << std::endl;
         os << "       d_vector_bdry_face_conds[" << j << "] = "
            << d_vector_bdry_face_conds[j] << std::endl;
         if (d_master_bdry_face_conds[j] == BdryCond::DIRICHLET) {
            os << "         d_bdry_face_density[" << j << "] = "
               << d_bdry_face_density[j] << std::endl;
            os << "         d_bdry_face_velocity[" << j << "] = "
               << d_bdry_face_velocity[j * d_dim.getValue() + 0] << " , "
               << d_bdry_face_velocity[j * d_dim.getValue() + 1] << " , "
               << d_bdry_face_velocity[j * d_dim.getValue() + 2] << std::endl;
            os << "         d_bdry_face_pressure[" << j << "] = "
               << d_bdry_face_pressure[j] << std::endl;
         }
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_master_bdry_edge_conds.size()); ++j) {
         os << "\n       d_master_bdry_edge_conds[" << j << "] = "
            << d_master_bdry_edge_conds[j] << std::endl;
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << std::endl;
         os << "       d_vector_bdry_edge_conds[" << j << "] = "
            << d_vector_bdry_edge_conds[j] << std::endl;
         os << "       d_edge_bdry_face[" << j << "] = "
            << d_edge_bdry_face[j] << std::endl;
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_master_bdry_node_conds.size()); ++j) {
         os << "\n       d_master_bdry_node_conds[" << j << "] = "
            << d_master_bdry_node_conds[j] << std::endl;
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << std::endl;
         os << "       d_vector_bdry_node_conds[" << j << "] = "
            << d_vector_bdry_node_conds[j] << std::endl;
         os << "       d_node_bdry_face[" << j << "] = "
            << d_node_bdry_face[j] << std::endl;
      }
   }

   os << "   Refinement criteria parameters " << std::endl;

   for (j = 0; j < static_cast<int>(d_refinement_criteria.size()); ++j) {
      os << "       d_refinement_criteria[" << j << "] = "
         << d_refinement_criteria[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_dev_tol.size()); ++j) {
      os << "       d_density_dev_tol[" << j << "] = "
         << d_density_dev_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_dev.size()); ++j) {
      os << "       d_density_dev[" << j << "] = "
         << d_density_dev[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_dev_time_max.size()); ++j) {
      os << "       d_density_dev_time_max[" << j << "] = "
         << d_density_dev_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_dev_time_min.size()); ++j) {
      os << "       d_density_dev_time_min[" << j << "] = "
         << d_density_dev_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_grad_tol.size()); ++j) {
      os << "       d_density_grad_tol[" << j << "] = "
         << d_density_grad_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_grad_time_max.size()); ++j) {
      os << "       d_density_grad_time_max[" << j << "] = "
         << d_density_grad_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_grad_time_min.size()); ++j) {
      os << "       d_density_grad_time_min[" << j << "] = "
         << d_density_grad_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_shock_onset.size()); ++j) {
      os << "       d_density_shock_onset[" << j << "] = "
         << d_density_shock_onset[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_shock_tol.size()); ++j) {
      os << "       d_density_shock_tol[" << j << "] = "
         << d_density_shock_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_shock_time_max.size()); ++j) {
      os << "       d_density_shock_time_max[" << j << "] = "
         << d_density_shock_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_shock_time_min.size()); ++j) {
      os << "       d_density_shock_time_min[" << j << "] = "
         << d_density_shock_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_rich_tol.size()); ++j) {
      os << "       d_density_rich_tol[" << j << "] = "
         << d_density_rich_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_rich_time_max.size()); ++j) {
      os << "       d_density_rich_time_max[" << j << "] = "
         << d_density_rich_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_density_rich_time_min.size()); ++j) {
      os << "       d_density_rich_time_min[" << j << "] = "
         << d_density_rich_time_min[j] << std::endl;
   }
   os << std::endl;

   for (j = 0; j < static_cast<int>(d_pressure_dev_tol.size()); ++j) {
      os << "       d_pressure_dev_tol[" << j << "] = "
         << d_pressure_dev_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_dev.size()); ++j) {
      os << "       d_pressure_dev[" << j << "] = "
         << d_pressure_dev[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_dev_time_max.size()); ++j) {
      os << "       d_pressure_dev_time_max[" << j << "] = "
         << d_pressure_dev_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_dev_time_min.size()); ++j) {
      os << "       d_pressure_dev_time_min[" << j << "] = "
         << d_pressure_dev_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_grad_tol.size()); ++j) {
      os << "       d_pressure_grad_tol[" << j << "] = "
         << d_pressure_grad_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_grad_time_max.size()); ++j) {
      os << "       d_pressure_grad_time_max[" << j << "] = "
         << d_pressure_grad_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_grad_time_min.size()); ++j) {
      os << "       d_pressure_grad_time_min[" << j << "] = "
         << d_pressure_grad_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_shock_onset.size()); ++j) {
      os << "       d_pressure_shock_onset[" << j << "] = "
         << d_pressure_shock_onset[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_shock_tol.size()); ++j) {
      os << "       d_pressure_shock_tol[" << j << "] = "
         << d_pressure_shock_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_shock_time_max.size()); ++j) {
      os << "       d_pressure_shock_time_max[" << j << "] = "
         << d_pressure_shock_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_shock_time_min.size()); ++j) {
      os << "       d_pressure_shock_time_min[" << j << "] = "
         << d_pressure_shock_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_rich_tol.size()); ++j) {
      os << "       d_pressure_rich_tol[" << j << "] = "
         << d_pressure_rich_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_rich_time_max.size()); ++j) {
      os << "       d_pressure_rich_time_max[" << j << "] = "
         << d_pressure_rich_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_pressure_rich_time_min.size()); ++j) {
      os << "       d_pressure_rich_time_min[" << j << "] = "
         << d_pressure_rich_time_min[j] << std::endl;
   }
   os << std::endl;

}

/*
 *************************************************************************
 *
 * Read data members from input.  Note all values set from restart
 * can be overridden by values in the input database.
 *
 *************************************************************************
 */

void Euler::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   /*
    * Note: if we are restarting, then we only allow nonuniform
    * workload to be used if nonuniform workload was used originally.
    */
   if (!is_from_restart) {
      d_use_nonuniform_workload =
         input_db->getBoolWithDefault("use_nonuniform_workload",
            d_use_nonuniform_workload);
   } else {
      if (d_use_nonuniform_workload) {
         d_use_nonuniform_workload =
            input_db->getBool("use_nonuniform_workload");
      }
   }

   if (!is_from_restart) {
      d_gamma = input_db->getDoubleWithDefault("gamma", d_gamma);
   }

   if (input_db->keyExists("riemann_solve")) {
      d_riemann_solve = input_db->getString("riemann_solve");
      if ((d_riemann_solve != "APPROX_RIEM_SOLVE") &&
          (d_riemann_solve != "EXACT_RIEM_SOLVE") &&
          (d_riemann_solve != "HLLC_RIEM_SOLVE")) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`riemann_solve' in input must be either string "
                          << "'APPROX_RIEM_SOLVE', 'EXACT_RIEM_SOLVE', "
                          << "'HLLC_RIEM_SOLVE'." << std::endl);

      }
   } else {
      d_riemann_solve = input_db->getStringWithDefault("d_riemann_solve",
            d_riemann_solve);
   }

   if (input_db->keyExists("godunov_order")) {
      d_godunov_order = input_db->getInteger("godunov_order");
      if ((d_godunov_order != 1) &&
          (d_godunov_order != 2) &&
          (d_godunov_order != 4)) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`godunov_order' in input must be 1, 2, or 4." << std::endl);

      }
   } else {
      d_godunov_order = input_db->getIntegerWithDefault("d_godunov_order",
            d_godunov_order);
   }

   if (input_db->keyExists("corner_transport")) {
      d_corner_transport = input_db->getString("corner_transport");
      if ((d_corner_transport != "CORNER_TRANSPORT_1") &&
          (d_corner_transport != "CORNER_TRANSPORT_2")) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`corner_transport' in input must be either string"
                          << " 'CORNER_TRANSPORT_1' or 'CORNER_TRANSPORT_2'." << std::endl);
      }
   } else {
      d_corner_transport = input_db->getStringWithDefault("corner_transport",
            d_corner_transport);
   }

   if (input_db->keyExists("Refinement_data")) {
      std::shared_ptr<tbox::Database> refine_db(
         input_db->getDatabase("Refinement_data"));
      std::vector<std::string> refinement_keys = refine_db->getAllKeys();
      int num_keys = static_cast<int>(refinement_keys.size());

      if (refine_db->keyExists("refine_criteria")) {
         d_refinement_criteria = refine_db->getStringVector("refine_criteria");
      } else {
         TBOX_WARNING(
            d_object_name << ": "
                          << "No key `refine_criteria' found in data for"
                          << " RefinementData. No refinement will occur." << std::endl);
      }

      std::vector<std::string> ref_keys_defined(num_keys);
      int def_key_cnt = 0;
      std::shared_ptr<tbox::Database> error_db;
      for (int i = 0; i < num_keys; ++i) {

         std::string error_key = refinement_keys[i];
         error_db.reset();

         if (!(error_key == "refine_criteria")) {

            if (!(error_key == "DENSITY_DEVIATION" ||
                  error_key == "DENSITY_GRADIENT" ||
                  error_key == "DENSITY_SHOCK" ||
                  error_key == "DENSITY_RICHARDSON" ||
                  error_key == "PRESSURE_DEVIATION" ||
                  error_key == "PRESSURE_GRADIENT" ||
                  error_key == "PRESSURE_SHOCK" ||
                  error_key == "PRESSURE_RICHARDSON")) {
               TBOX_ERROR(
                  d_object_name << ": "
                                << "Unknown refinement criteria: "
                                << error_key
                                << "\nin input." << std::endl);
            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               ++def_key_cnt;
            }

            if (error_db && error_key == "DENSITY_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_density_dev_tol = error_db->getDoubleVector("dev_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `dev_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("density_dev")) {
                  d_density_dev = error_db->getDoubleVector("density_dev");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `density_dev' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_density_dev_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_density_dev_time_max.resize(1);
                  d_density_dev_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_density_dev_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_density_dev_time_min.resize(1);
                  d_density_dev_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "DENSITY_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_density_grad_tol = error_db->getDoubleVector("grad_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `grad_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_density_grad_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_density_grad_time_max.resize(1);
                  d_density_grad_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_density_grad_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_density_grad_time_min.resize(1);
                  d_density_grad_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "DENSITY_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_density_shock_onset =
                     error_db->getDoubleVector("shock_onset");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_onset' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_density_shock_tol = error_db->getDoubleVector("shock_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_density_shock_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_density_shock_time_max.resize(1);
                  d_density_shock_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_density_shock_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_density_shock_time_min.resize(1);
                  d_density_shock_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "DENSITY_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_density_rich_tol = error_db->getDoubleVector("rich_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `rich_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_density_rich_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_density_rich_time_max.resize(1);
                  d_density_rich_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_density_rich_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_density_rich_time_min.resize(1);
                  d_density_rich_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "PRESSURE_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_pressure_dev_tol = error_db->getDoubleVector("dev_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `dev_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("pressure_dev")) {
                  d_pressure_dev = error_db->getDoubleVector("pressure_dev");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `pressure_dev' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_pressure_dev_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_pressure_dev_time_max.resize(1);
                  d_pressure_dev_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_pressure_dev_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_pressure_dev_time_min.resize(1);
                  d_pressure_dev_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "PRESSURE_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_pressure_grad_tol = error_db->getDoubleVector("grad_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `grad_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_pressure_grad_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_pressure_grad_time_max.resize(1);
                  d_pressure_grad_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_pressure_grad_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_pressure_grad_time_min.resize(1);
                  d_pressure_grad_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "PRESSURE_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_pressure_shock_onset =
                     error_db->getDoubleVector("shock_onset");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_onset' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_pressure_shock_tol =
                     error_db->getDoubleVector("shock_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_pressure_shock_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_pressure_shock_time_max.resize(1);
                  d_pressure_shock_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_pressure_shock_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_pressure_shock_time_min.resize(1);
                  d_pressure_shock_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "PRESSURE_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_pressure_rich_tol = error_db->getDoubleVector("rich_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `rich_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_pressure_rich_time_max =
                     error_db->getDoubleVector("time_max");
               } else {
                  d_pressure_rich_time_max.resize(1);
                  d_pressure_rich_time_max[0] =
                     tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_pressure_rich_time_min =
                     error_db->getDoubleVector("time_min");
               } else {
                  d_pressure_rich_time_min.resize(1);
                  d_pressure_rich_time_min[0] = 0.;
               }

            }

         }

      } // loop over refine criteria

      /*
       * Check that input is found for each string identifier in key list.
       */
      for (int k0 = 0;
           k0 < static_cast<int>(d_refinement_criteria.size()); ++k0) {
         std::string use_key = d_refinement_criteria[k0];
         bool key_found = false;
         for (int k1 = 0; k1 < def_key_cnt; ++k1) {
            std::string def_key = ref_keys_defined[k1];
            if (def_key == use_key) key_found = true;
         }

         if (!key_found) {
            TBOX_ERROR(d_object_name << ": "
                                     << "No input found for specified refine criteria: "
                                     << d_refinement_criteria[k0] << std::endl);
         }
      }

   } // if "Refinement_data" db entry exists

   if (!is_from_restart) {

      if (input_db->keyExists("data_problem")) {
         d_data_problem = input_db->getString("data_problem");
      } else {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`data_problem' value not found in input."
                          << std::endl);
      }

      if (!input_db->keyExists("Initial_data")) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "No `Initial_data' database found in input." << std::endl);
      }
      std::shared_ptr<tbox::Database> init_data_db(
         input_db->getDatabase("Initial_data"));

      bool found_problem_data = false;

      if (d_data_problem == "SPHERE") {

         if (init_data_db->keyExists("radius")) {
            d_radius = init_data_db->getDouble("radius");
         } else {
            TBOX_ERROR(
               d_object_name << ": "
                             << "`radius' input required for SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("center")) {
            init_data_db->getDoubleArray("center", d_center, d_dim.getValue());
         } else {
            TBOX_ERROR(
               d_object_name << ": "
                             << "`center' input required for SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("density_inside")) {
            d_density_inside = init_data_db->getDouble("density_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`density_inside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("velocity_inside")) {
            init_data_db->getDoubleArray("velocity_inside",
               d_velocity_inside, d_dim.getValue());
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`velocity_inside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("pressure_inside")) {
            d_pressure_inside = init_data_db->getDouble("pressure_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`pressure_inside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("density_outside")) {
            d_density_outside = init_data_db->getDouble("density_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`density_outside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("velocity_outside")) {
            init_data_db->getDoubleArray("velocity_outside",
               d_velocity_outside, d_dim.getValue());
         } else {
            TBOX_ERROR(
               d_object_name << ": "
                             << "`velocity_outside' input required for "
                             << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("pressure_outside")) {
            d_pressure_outside = init_data_db->getDouble("pressure_outside");
         } else {
            TBOX_ERROR(
               d_object_name << ": "
                             << "`pressure_outside' input required for "
                             << "SPHERE problem." << std::endl);
         }

         found_problem_data = true;

      }

      if (!found_problem_data &&
          ((d_data_problem == "PIECEWISE_CONSTANT_X") ||
           (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
           (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
           (d_data_problem == "STEP"))) {

         int idir = 0;
         if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
            idir = 1;
         }

         if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
            if (d_dim < tbox::Dimension(3)) {
               TBOX_ERROR(
                  d_object_name << ": `PIECEWISE_CONSTANT_Z' "
                                << "problem invalid in 2 dimensions."
                                << std::endl);
            }
            idir = 2;
         }

         std::vector<std::string> init_data_keys = init_data_db->getAllKeys();

         if (init_data_db->keyExists("front_position")) {
            d_front_position = init_data_db->getDoubleVector("front_position");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`front_position' input required for "
                                     << "PIECEWISE_CONSTANT_* problem." << std::endl);
         }

         d_number_of_intervals =
            tbox::MathUtilities<int>::Min(static_cast<int>(d_front_position.size()) + 1,
               static_cast<int>(init_data_keys.size()) - 1);

         d_front_position.resize(static_cast<int>(d_front_position.size()) + 1);
         d_front_position[static_cast<int>(d_front_position.size()) - 1] =
            d_grid_geometry->getXUpper()[idir];

         d_interval_density.resize(d_number_of_intervals);
         d_interval_velocity.resize(d_number_of_intervals * d_dim.getValue());
         d_interval_pressure.resize(d_number_of_intervals);

         int i = 0;
         int nkey = 0;
         bool found_interval_data = false;

         while (!found_interval_data
                && (i < d_number_of_intervals)
                && (nkey < static_cast<int>(init_data_keys.size()))) {

            if (!(init_data_keys[nkey] == "front_position")) {

               std::shared_ptr<tbox::Database> interval_db(
                  init_data_db->getDatabase(init_data_keys[nkey]));

               readStateDataEntry(interval_db,
                  init_data_keys[nkey],
                  i,
                  d_interval_density,
                  d_interval_velocity,
                  d_interval_pressure);

               ++i;

               found_interval_data = (i == d_number_of_intervals);

            }

            ++nkey;

         }

         if (!found_interval_data) {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Insufficient interval data given in input"
                             << " for PIECEWISE_CONSTANT_* or STEP problem." << std::endl);
         }

         found_problem_data = true;

      }

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
                                  << "`Initial_data' database found in input."
                                  << " But bad data supplied." << std::endl);
      }

   } // if !is_from_restart read in problem data

   const hier::IntVector& one_vec = hier::IntVector::getOne(d_dim);
   hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);
   int num_per_dirs = 0;
   for (int id = 0; id < d_dim.getValue(); ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < d_dim.getValue()) {

      if (input_db->keyExists("Boundary_data")) {

         std::shared_ptr<tbox::Database> bdry_db = input_db->getDatabase(
               "Boundary_data");

         if (d_dim == tbox::Dimension(2)) {
            appu::CartesianBoundaryUtilities2::getFromInput(this,
               bdry_db,
               d_master_bdry_edge_conds,
               d_master_bdry_node_conds,
               periodic);
         }
         if (d_dim == tbox::Dimension(3)) {
            appu::CartesianBoundaryUtilities3::getFromInput(this,
               bdry_db,
               d_master_bdry_face_conds,
               d_master_bdry_edge_conds,
               d_master_bdry_node_conds,
               periodic);
         }

      } else {
         TBOX_ERROR(
            d_object_name << ": "
                          << "Key data `Boundary_data' not found in input. " << std::endl);
      }

   }

}

/*
 *************************************************************************
 *
 * Routines to put/get data members to/from from restart database.
 *
 *************************************************************************
 */

void Euler::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("EULER_VERSION", EULER_VERSION);

   restart_db->putDouble("d_gamma", d_gamma);

   restart_db->putString("d_riemann_solve", d_riemann_solve);
   restart_db->putInteger("d_godunov_order", d_godunov_order);
   restart_db->putString("d_corner_transport", d_corner_transport);
   restart_db->putIntegerArray("d_nghosts", &d_nghosts[0], d_dim.getValue());
   restart_db->putIntegerArray("d_fluxghosts",
      &d_fluxghosts[0],
      d_dim.getValue());

   restart_db->putString("d_data_problem", d_data_problem);

   if (d_data_problem == "SPHERE") {
      restart_db->putDouble("d_radius", d_radius);
      restart_db->putDoubleArray("d_center", d_center, d_dim.getValue());
      restart_db->putDouble("d_density_inside", d_density_inside);
      restart_db->putDoubleArray("d_velocity_inside",
         d_velocity_inside,
         d_dim.getValue());
      restart_db->putDouble("d_pressure_inside", d_pressure_inside);
      restart_db->putDouble("d_density_outside", d_density_outside);
      restart_db->putDoubleArray("d_velocity_outside",
         d_velocity_outside,
         d_dim.getValue());
      restart_db->putDouble("d_pressure_outside", d_pressure_outside);
   }

   if ((d_data_problem == "PIECEWISE_CONSTANT_X") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
       (d_data_problem == "STEP")) {
      restart_db->putInteger("d_number_of_intervals", d_number_of_intervals);
      if (d_number_of_intervals > 0) {
         restart_db->putDoubleVector("d_front_position", d_front_position);
         restart_db->putDoubleVector("d_interval_density", d_interval_density);
         restart_db->putDoubleVector("d_interval_velocity", d_interval_velocity);
         restart_db->putDoubleVector("d_interval_pressure", d_interval_pressure);
      }
   }

   restart_db->putIntegerVector("d_master_bdry_edge_conds",
      d_master_bdry_edge_conds);
   restart_db->putIntegerVector("d_master_bdry_node_conds",
      d_master_bdry_node_conds);

   if (d_dim == tbox::Dimension(2)) {
      restart_db->putDoubleVector("d_bdry_edge_density",
         d_bdry_edge_density);
      restart_db->putDoubleVector("d_bdry_edge_velocity",
         d_bdry_edge_velocity);
      restart_db->putDoubleVector("d_bdry_edge_pressure",
         d_bdry_edge_pressure);
   }
   if (d_dim == tbox::Dimension(3)) {
      restart_db->putIntegerVector("d_master_bdry_face_conds",
         d_master_bdry_face_conds);

      restart_db->putDoubleVector("d_bdry_face_density",
         d_bdry_face_density);
      restart_db->putDoubleVector("d_bdry_face_velocity",
         d_bdry_face_velocity);
      restart_db->putDoubleVector("d_bdry_face_pressure",
         d_bdry_face_pressure);
   }

   if (d_refinement_criteria.size() > 0) {
      restart_db->putStringVector("d_refinement_criteria",
         d_refinement_criteria);
   }
   for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); ++i) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         restart_db->putDoubleVector("d_density_dev_tol",
            d_density_dev_tol);
         restart_db->putDoubleVector("d_density_dev",
            d_density_dev);
         restart_db->putDoubleVector("d_density_dev_time_max",
            d_density_dev_time_max);
         restart_db->putDoubleVector("d_density_dev_time_min",
            d_density_dev_time_min);

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         restart_db->putDoubleVector("d_density_grad_tol",
            d_density_grad_tol);
         restart_db->putDoubleVector("d_density_grad_time_max",
            d_density_grad_time_max);
         restart_db->putDoubleVector("d_density_grad_time_min",
            d_density_grad_time_min);

      } else if (d_refinement_criteria[i] == "DENSITY_SHOCK") {

         restart_db->putDoubleVector("d_density_shock_onset",
            d_density_shock_onset);
         restart_db->putDoubleVector("d_density_shock_tol",
            d_density_shock_tol);
         restart_db->putDoubleVector("d_density_shock_time_max",
            d_density_shock_time_max);
         restart_db->putDoubleVector("d_density_shock_time_min",
            d_density_shock_time_min);

      } else if (d_refinement_criteria[i] == "DENSITY_RICHARDSON") {

         restart_db->putDoubleVector("d_density_rich_tol",
            d_density_rich_tol);
         restart_db->putDoubleVector("d_density_rich_time_max",
            d_density_rich_time_max);
         restart_db->putDoubleVector("d_density_rich_time_min",
            d_density_rich_time_min);

      } else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         restart_db->putDoubleVector("d_pressure_dev_tol",
            d_pressure_dev_tol);
         restart_db->putDoubleVector("d_pressure_dev",
            d_pressure_dev);
         restart_db->putDoubleVector("d_pressure_dev_time_max",
            d_pressure_dev_time_max);
         restart_db->putDoubleVector("d_pressure_dev_time_min",
            d_pressure_dev_time_min);

      } else if (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         restart_db->putDoubleVector("d_pressure_grad_tol",
            d_pressure_grad_tol);
         restart_db->putDoubleVector("d_pressure_grad_time_max",
            d_pressure_grad_time_max);
         restart_db->putDoubleVector("d_pressure_grad_time_min",
            d_pressure_grad_time_min);

      } else if (d_refinement_criteria[i] == "PRESSURE_SHOCK") {

         restart_db->putDoubleVector("d_pressure_shock_onset",
            d_pressure_shock_onset);
         restart_db->putDoubleVector("d_pressure_shock_tol",
            d_pressure_shock_tol);
         restart_db->putDoubleVector("d_pressure_shock_time_max",
            d_pressure_shock_time_max);
         restart_db->putDoubleVector("d_pressure_shock_time_min",
            d_pressure_shock_time_min);

      } else if (d_refinement_criteria[i] == "PRESSURE_RICHARDSON") {

         restart_db->putDoubleVector("d_pressure_rich_tol",
            d_pressure_rich_tol);
         restart_db->putDoubleVector("d_pressure_rich_time_max",
            d_pressure_rich_time_max);
         restart_db->putDoubleVector("d_pressure_rich_time_min",
            d_pressure_rich_time_min);
      }

   }

}

void Euler::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file." << std::endl);
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("EULER_VERSION");
   if (ver != EULER_VERSION) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Restart file version different than class version." << std::endl);
   }

   d_gamma = db->getDouble("d_gamma");

   d_riemann_solve = db->getString("d_riemann_solve");
   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = &d_nghosts[0];
   db->getIntegerArray("d_nghosts", tmp_nghosts, d_dim.getValue());
   for (int i = 0; i < d_dim.getValue(); ++i) {
      if (d_nghosts(i) != CELLG) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "Key data `d_nghosts' in restart file != CELLG." << std::endl);
      }
   }
   int* tmp_fluxghosts = &d_fluxghosts[0];
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, d_dim.getValue());
   for (int i = 0; i < d_dim.getValue(); ++i) {
      if (d_fluxghosts(i) != FLUXG) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "Key data `d_fluxghosts' in restart file != FLUXG." << std::endl);
      }
   }

   d_data_problem = db->getString("d_data_problem");

   if (d_data_problem == "SPHERE") {
      d_radius = db->getDouble("d_radius");
      db->getDoubleArray("d_center", d_center, d_dim.getValue());
      d_density_inside = db->getDouble("d_density_inside");
      db->getDoubleArray("d_velocity_inside", d_velocity_inside, d_dim.getValue());
      d_pressure_inside = db->getDouble("d_pressure_inside");
      d_density_outside = db->getDouble("d_density_outside");
      db->getDoubleArray("d_velocity_outside", d_velocity_outside, d_dim.getValue());
      d_pressure_outside = db->getDouble("d_pressure_outside");
   }

   if ((d_data_problem == "PIECEWISE_CONSTANT_X") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
       (d_data_problem == "STEP")) {
      d_number_of_intervals = db->getInteger("d_number_of_intervals");
      if (d_number_of_intervals > 0) {
         d_front_position = db->getDoubleVector("d_front_position");
         d_interval_density = db->getDoubleVector("d_interval_density");
         d_interval_velocity = db->getDoubleVector("d_interval_velocity");
         d_interval_pressure = db->getDoubleVector("d_interval_pressure");
      }
   }

   d_master_bdry_edge_conds = db->getIntegerVector("d_master_bdry_edge_conds");
   d_master_bdry_node_conds = db->getIntegerVector("d_master_bdry_node_conds");

   if (d_dim == tbox::Dimension(2)) {
      d_bdry_edge_density = db->getDoubleVector("d_bdry_edge_density");
      d_bdry_edge_velocity = db->getDoubleVector("d_bdry_edge_velocity");
      d_bdry_edge_pressure = db->getDoubleVector("d_bdry_edge_pressure");
   }
   if (d_dim == tbox::Dimension(3)) {
      d_master_bdry_face_conds =
         db->getIntegerVector("d_master_bdry_face_conds");

      d_bdry_face_density = db->getDoubleVector("d_bdry_face_density");
      d_bdry_face_velocity = db->getDoubleVector("d_bdry_face_velocity");
      d_bdry_face_pressure = db->getDoubleVector("d_bdry_face_pressure");
   }

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringVector("d_refinement_criteria");
   }

   for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); ++i) {

      if (d_refinement_criteria[i] == "DENSITY_DEVIATION") {

         d_density_dev_tol =
            db->getDoubleVector("d_density_dev_tol");
         d_density_dev =
            db->getDoubleVector("d_density_dev");
         d_density_dev_time_max =
            db->getDoubleVector("d_density_dev_time_max");
         d_density_dev_time_min =
            db->getDoubleVector("d_density_dev_time_min");

      } else if (d_refinement_criteria[i] == "DENSITY_GRADIENT") {

         d_density_grad_tol =
            db->getDoubleVector("d_density_grad_tol");
         d_density_grad_time_max =
            db->getDoubleVector("d_density_grad_time_max");
         d_density_grad_time_min =
            db->getDoubleVector("d_density_grad_time_min");

      } else if (d_refinement_criteria[i] == "DENSITY_SHOCK") {

         d_density_shock_onset =
            db->getDoubleVector("d_density_shock_onset");
         d_density_shock_tol =
            db->getDoubleVector("d_density_shock_tol");
         d_density_shock_time_max =
            db->getDoubleVector("d_density_shock_time_max");
         d_density_shock_time_min =
            db->getDoubleVector("d_density_shock_time_min");

      } else if (d_refinement_criteria[i] == "DENSITY_RICHARDSON") {

         d_density_rich_tol =
            db->getDoubleVector("d_density_rich_tol");
         d_density_rich_time_max =
            db->getDoubleVector("d_density_rich_time_max");
         d_density_rich_time_min =
            db->getDoubleVector("d_density_rich_time_min");

      } else if (d_refinement_criteria[i] == "PRESSURE_DEVIATION") {

         d_pressure_dev_tol =
            db->getDoubleVector("d_pressure_dev_tol");
         d_pressure_dev =
            db->getDoubleVector("d_pressure_dev");
         d_pressure_dev_time_max =
            db->getDoubleVector("d_pressure_dev_time_max");
         d_pressure_dev_time_min =
            db->getDoubleVector("d_pressure_dev_time_min");

      } else if (d_refinement_criteria[i] == "PRESSURE_GRADIENT") {

         d_pressure_grad_tol =
            db->getDoubleVector("d_pressure_grad_tol");
         d_pressure_grad_time_max =
            db->getDoubleVector("d_pressure_grad_time_max");
         d_pressure_grad_time_min =
            db->getDoubleVector("d_pressure_grad_time_min");

      } else if (d_refinement_criteria[i] == "PRESSURE_SHOCK") {

         d_pressure_shock_onset =
            db->getDoubleVector("d_pressure_shock_onset");
         d_pressure_shock_tol =
            db->getDoubleVector("d_pressure_shock_tol");
         d_pressure_shock_time_max =
            db->getDoubleVector("d_pressure_shock_time_max");
         d_pressure_shock_time_min =
            db->getDoubleVector("d_pressure_shock_time_min");

      } else if (d_refinement_criteria[i] == "PRESSURE_RICHARDSON") {

         d_pressure_rich_tol =
            db->getDoubleVector("d_pressure_rich_tol");
         d_pressure_rich_time_max =
            db->getDoubleVector("d_pressure_rich_time_max");
         d_pressure_rich_time_min =
            db->getDoubleVector("d_pressure_rich_time_min");

      }

   }

}

/*
 *************************************************************************
 *
 * Routines to read boundary data from input database.
 *
 *************************************************************************
 */

void Euler::readDirichletBoundaryDataEntry(
   const std::shared_ptr<tbox::Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());

   if (d_dim == tbox::Dimension(2)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_edge_density,
         d_bdry_edge_velocity,
         d_bdry_edge_pressure);
   }
   if (d_dim == tbox::Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_density,
         d_bdry_face_velocity,
         d_bdry_face_pressure);
   }
}

void Euler::readNeumannBoundaryDataEntry(
   const std::shared_ptr<tbox::Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   NULL_USE(db);
   NULL_USE(db_name);
   NULL_USE(bdry_location_index);
}

void Euler::readStateDataEntry(
   std::shared_ptr<tbox::Database> db,
   const std::string& db_name,
   int array_indx,
   std::vector<double>& density,
   std::vector<double>& velocity,
   std::vector<double>& pressure)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(static_cast<int>(density.size()) > array_indx);
   TBOX_ASSERT(static_cast<int>(velocity.size()) > array_indx * d_dim.getValue());
   TBOX_ASSERT(static_cast<int>(pressure.size()) > array_indx);

   if (db->keyExists("density")) {
      density[array_indx] = db->getDouble("density");
   } else {
      TBOX_ERROR(d_object_name << ": "
                               << "`density' entry missing from " << db_name
                               << " input database. " << std::endl);
   }
   if (db->keyExists("velocity")) {
      std::vector<double> tmp_vel = db->getDoubleVector("velocity");
      if (static_cast<int>(tmp_vel.size()) < d_dim.getValue()) {
         TBOX_ERROR(d_object_name << ": "
                                  << "Insufficient number `velocity' values"
                                  << " given in " << db_name
                                  << " input database." << std::endl);
      }
      for (int iv = 0; iv < d_dim.getValue(); ++iv) {
         velocity[array_indx * d_dim.getValue() + iv] = tmp_vel[iv];
      }
   } else {
      TBOX_ERROR(d_object_name << ": "
                               << "`velocity' entry missing from " << db_name
                               << " input database. " << std::endl);
   }
   if (db->keyExists("pressure")) {
      pressure[array_indx] = db->getDouble("pressure");
   } else {
      TBOX_ERROR(d_object_name << ": "
                               << "`pressure' entry missing from " << db_name
                               << " input database. " << std::endl);
   }

}

/*
 *************************************************************************
 *
 * Routine to check boundary data when debugging.
 *
 *************************************************************************
 */

void Euler::checkBoundaryData(
   int btype,
   const hier::Patch& patch,
   const hier::IntVector& ghost_width_to_check,
   const std::vector<int>& scalar_bconds,
   const std::vector<int>& vector_bconds) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim == tbox::Dimension(2)) {
      TBOX_ASSERT(btype == Bdry::EDGE2D || btype == Bdry::NODE2D);
   }
   if (d_dim == tbox::Dimension(3)) {
      TBOX_ASSERT(btype == Bdry::FACE3D || btype == Bdry::EDGE3D || btype == Bdry::NODE3D);
   }
#endif

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const std::vector<hier::BoundaryBox>& bdry_boxes =
      pgeom->getCodimensionBoundaries(btype);

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   for (int i = 0; i < static_cast<int>(bdry_boxes.size()); ++i) {
      hier::BoundaryBox bbox = bdry_boxes[i];
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
      int bloc = bbox.getLocationIndex();

      int bscalarcase = 0;
      int bvelocitycase = 0;
      int refbdryloc = 0;
      if (d_dim == tbox::Dimension(2)) {
         if (btype == Bdry::EDGE2D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_2D_EDGES);
            TBOX_ASSERT(static_cast<int>(vector_bconds.size()) ==
               NUM_2D_EDGES);

            bscalarcase = scalar_bconds[bloc];
            bvelocitycase = vector_bconds[bloc];
            refbdryloc = bloc;
         } else { // btype == Bdry::NODE2D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_2D_NODES);
            TBOX_ASSERT(static_cast<int>(vector_bconds.size()) ==
               NUM_2D_NODES);

            bscalarcase = scalar_bconds[bloc];
            bvelocitycase = vector_bconds[bloc];
            refbdryloc = d_node_bdry_edge[bloc];
         }
      }
      if (d_dim == tbox::Dimension(3)) {
         if (btype == Bdry::FACE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_FACES);
            TBOX_ASSERT(static_cast<int>(vector_bconds.size()) ==
               NUM_3D_FACES);

            bscalarcase = scalar_bconds[bloc];
            bvelocitycase = vector_bconds[bloc];
            refbdryloc = bloc;
         } else if (btype == Bdry::EDGE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_EDGES);
            TBOX_ASSERT(static_cast<int>(vector_bconds.size()) ==
               NUM_3D_EDGES);

            bscalarcase = scalar_bconds[bloc];
            bvelocitycase = vector_bconds[bloc];
            refbdryloc = d_edge_bdry_face[bloc];
         } else { // btype == Bdry::NODE3D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_NODES);
            TBOX_ASSERT(static_cast<int>(vector_bconds.size()) ==
               NUM_3D_NODES);

            bscalarcase = scalar_bconds[bloc];
            bvelocitycase = vector_bconds[bloc];
            refbdryloc = d_node_bdry_face[bloc];
         }
      }

      int num_bad_values = 0;

      if (d_dim == tbox::Dimension(2)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities2::checkBdryData(
               d_density->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_density,
                  getDataContext()), 0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_edge_density[refbdryloc]);
      }
      if (d_dim == tbox::Dimension(3)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities3::checkBdryData(
               d_density->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_density,
                  getDataContext()), 0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_face_density[refbdryloc]);
      }
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nEuler Boundary Test FAILED: \n"
                    << "     " << num_bad_values
                    << " bad DENSITY values found for\n"
                    << "     boundary type " << btype << " at location "
                    << bloc << std::endl;
      }
#endif

      if (d_dim == tbox::Dimension(2)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities2::checkBdryData(
               d_pressure->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_pressure, getDataContext()),
               0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_edge_density[refbdryloc]);
      }
      if (d_dim == tbox::Dimension(3)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities3::checkBdryData(
               d_pressure->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_pressure, getDataContext()),
               0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_face_density[refbdryloc]);
      }
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nEuler Boundary Test FAILED: \n"
                    << "     " << num_bad_values
                    << " bad PRESSURE values found for\n"
                    << "     boundary type " << btype << " at location "
                    << bloc << std::endl;
      }
#endif

      for (int idir = 0; idir < d_dim.getValue(); ++idir) {

         int vbcase = bscalarcase;
         if (d_dim == tbox::Dimension(2)) {
            if (btype == Bdry::EDGE2D) {
               if ((idir == 0 && (bloc == BdryLoc::XLO || bloc == BdryLoc::XHI)) ||
                   (idir == 1 && (bloc == BdryLoc::YLO || bloc == BdryLoc::YHI))) {
                  vbcase = bvelocitycase;
               }
            } else if (btype == Bdry::NODE2D) {
               if ((idir == 0 && bvelocitycase == BdryCond::XREFLECT) ||
                   (idir == 1 && bvelocitycase == BdryCond::YREFLECT)) {
                  vbcase = bvelocitycase;
               }
            }
         }
         if (d_dim == tbox::Dimension(3)) {
            if (btype == Bdry::FACE3D) {
               if ((idir == 0 && (bloc == BdryLoc::XLO || bloc == BdryLoc::XHI)) ||
                   (idir == 1 && (bloc == BdryLoc::YLO || bloc == BdryLoc::YHI)) ||
                   (idir == 2 && (bloc == BdryLoc::ZLO || bloc == BdryLoc::ZHI))) {
                  vbcase = bvelocitycase;
               }
            } else if (btype == Bdry::EDGE3D || btype == Bdry::NODE3D) {
               if ((idir == 0 && bvelocitycase == BdryCond::XREFLECT) ||
                   (idir == 1 && bvelocitycase == BdryCond::YREFLECT) ||
                   (idir == 2 && bvelocitycase == BdryCond::ZREFLECT)) {
                  vbcase = bvelocitycase;
               }
            }
         }

         if (d_dim == tbox::Dimension(2)) {
            num_bad_values =
               appu::CartesianBoundaryUtilities2::checkBdryData(
                  d_velocity->getName(),
                  patch,
                  vdb->mapVariableAndContextToIndex(d_velocity, getDataContext()),
                  idir,
                  ghost_width_to_check,
                  bbox,
                  vbcase,
                  d_bdry_edge_velocity[refbdryloc * d_dim.getValue() + idir]);
         }
         if (d_dim == tbox::Dimension(3)) {
            num_bad_values =
               appu::CartesianBoundaryUtilities3::checkBdryData(
                  d_velocity->getName(),
                  patch,
                  vdb->mapVariableAndContextToIndex(d_velocity, getDataContext()),
                  idir,
                  ghost_width_to_check,
                  bbox,
                  vbcase,
                  d_bdry_face_velocity[refbdryloc * d_dim.getValue() + idir]);
         }
#if (TESTING == 1)
         if (num_bad_values > 0) {
            tbox::perr << "\nEuler Boundary Test FAILED: \n"
                       << "     " << num_bad_values
                       << " bad VELOCITY values found in direction " << idir
                       << " for\n"
                       << "     boundary type " << btype << " at location "
                       << bloc << std::endl;
         }
#endif
      }

   }

}
