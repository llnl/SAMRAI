/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in linear advection ex.
 *
 ************************************************************************/
#include "LinAdv.h"

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
#include <strstream.h>
#endif
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (0)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

// External definitions for Fortran numerical routines
#include "LinAdvFort.h"

// Number of ghosts cells used for each variable quantity
#define CELLG (4)
#define FACEG (4)
#define FLUXG (1)

// defines for initialization
#define PIECEWISE_CONSTANT_X (10)
#define PIECEWISE_CONSTANT_Y (11)
#define PIECEWISE_CONSTANT_Z (12)
#define SINE_CONSTANT_X (20)
#define SINE_CONSTANT_Y (21)
#define SINE_CONSTANT_Z (22)
#define SPHERE (40)

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

// Version of LinAdv restart file data
#define LINADV_VERSION (3)

/*
 *************************************************************************
 *
 * The constructor for LinAdv class sets data members to defualt values,
 * creates variables that define the solution state for the linear
 * advection equation.
 *
 * After default values are set, this routine calls getFromRestart()
 * if execution from a restart file is specified.  Finally,
 * getFromInput() is called to read values from the given input
 * database (potentially overriding those found in the restart file).
 *
 *************************************************************************
 */

LinAdv::LinAdv(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
   const std::shared_ptr<MeshGenerationStrategy>& sine_wall):
   algs::HyperbolicPatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_mesh_gen(sine_wall),
   d_grid_geometry(grid_geom),
   d_use_nonuniform_workload(false),
   d_uval(new pdat::CellVariable<double>(dim, "uval", 1)),
   d_flux(new pdat::FaceVariable<double>(dim, "flux", 1)),
   d_godunov_order(1),
   d_corner_transport("CORNER_TRANSPORT_1"),
   d_nghosts(hier::IntVector(dim, CELLG)),
   d_fluxghosts(hier::IntVector(dim, FLUXG))
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(grid_geom);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   if (!t_init) {
      t_init = tbox::TimerManager::getManager()->
         getTimer("apps::LinAdv::initializeDataOnPatch()");
      t_init_first_time = tbox::TimerManager::getManager()->
         getTimer("apps::LinAdv::initializeDataOnPatch()_initial_time");
      t_analytical_tag = tbox::TimerManager::getManager()->
         getTimer("apps::LinAdv::analytical_tag");
   }

   TBOX_ASSERT(CELLG == FACEG);

   /*
    * Defaults for problem type and initial data.
    */

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(stufprobc2d, STUFPROBC2D) (PIECEWISE_CONSTANT_X,
         PIECEWISE_CONSTANT_Y, PIECEWISE_CONSTANT_Z,
         SINE_CONSTANT_X, SINE_CONSTANT_Y, SINE_CONSTANT_Z, SPHERE,
         CELLG, FACEG, FLUXG);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(stufprobc3d, STUFPROBC3D) (PIECEWISE_CONSTANT_X,
         PIECEWISE_CONSTANT_Y, PIECEWISE_CONSTANT_Z,
         SINE_CONSTANT_X, SINE_CONSTANT_Y, SINE_CONSTANT_Z, SPHERE,
         CELLG, FACEG, FLUXG);
   }

}

/*
 *************************************************************************
 *
 * Empty destructor for LinAdv class.
 *
 *************************************************************************
 */

LinAdv::~LinAdv() {
}

/*
 *************************************************************************
 *
 * Register conserved variable (u) (i.e., solution state variable) and
 * flux variable with hyperbolic integrator that manages storage for
 * those quantities.  Also, register plot data with VisIt
 *
 *************************************************************************
 */

void LinAdv::registerModelVariables(
   algs::HyperbolicLevelIntegrator* integrator)
{

   TBOX_ASSERT(integrator != 0);
   TBOX_ASSERT(CELLG == FACEG);

   integrator->registerVariable(d_uval, d_nghosts,
      algs::HyperbolicLevelIntegrator::TIME_DEP,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "CONSERVATIVE_LINEAR_REFINE");

   integrator->registerVariable(d_flux, d_fluxghosts,
      algs::HyperbolicLevelIntegrator::FLUX,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "NO_REFINE");

}

/*
 *************************************************************************
 *
 * Set up parameters for nonuniform load balancing, if used.
 *
 *************************************************************************
 */

void LinAdv::setupLoadBalancer(
   algs::HyperbolicLevelIntegrator* integrator,
   mesh::GriddingAlgorithm* gridding_algorithm)
{
   NULL_USE(integrator);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();
   hier::PatchDataRestartManager* pdrm =
      hier::PatchDataRestartManager::getManager();

   if (d_use_nonuniform_workload && gridding_algorithm) {
      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         std::dynamic_pointer_cast<mesh::TreeLoadBalancer, mesh::LoadBalanceStrategy>(
            gridding_algorithm->getLoadBalanceStrategy()));

      if (load_balancer) {
         d_workload_variable.reset(new pdat::CellVariable<double>(
               d_dim,
               "workload_variable",
               1));
         d_workload_data_id =
            vardb->registerVariableAndContext(d_workload_variable,
               vardb->getContext("WORKLOAD"),
               hier::IntVector(d_dim, 0));
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
void LinAdv::initializeDataOnPatch(
   hier::Patch& patch,
   const double data_time,
   const bool initial_time)
{
   t_init->start();

   if (initial_time) {

      t_init_first_time->start();
      const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(pgeom);

      std::shared_ptr<pdat::CellData<double> > uval(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_uval, getDataContext())));
      TBOX_ASSERT(uval);
      uval->setTime(data_time);

      d_mesh_gen->computePatchData(
         patch,
         uval.get(), 0,
         patch.getBox());

      t_init_first_time->stop();
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

double LinAdv::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time,
   const double dt_time)
{
   NULL_USE(initial_time);
   NULL_USE(dt_time);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));

   TBOX_ASSERT(uval);

   hier::IntVector ghost_cells = uval->getGhostCellWidth();

   double stabdt;
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(stabledt2d, STABLEDT2D) (dx,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ghost_cells(0),
         ghost_cells(1),
         d_advection_velocity,
         uval->getPointer(),
         stabdt);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(stabledt3d, STABLEDT3D) (dx,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         ghost_cells(0),
         ghost_cells(1),
         ghost_cells(2),
         d_advection_velocity,
         uval->getPointer(),
         stabdt);
   } else {
      TBOX_ERROR("Only 2D or 3D allowed in LinAdv::computeStableDtOnPatch");
      stabdt = 0;
   }

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

void LinAdv::computeFluxesOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt)
{
   NULL_USE(time);

   if (d_dim == tbox::Dimension(3)) {

      if (d_corner_transport == "CORNER_TRANSPORT_2") {
         compute3DFluxesWithCornerTransport2(patch, dt);
      } else {
         compute3DFluxesWithCornerTransport1(patch, dt);
      }

   }

   if (d_dim < tbox::Dimension(3)) {

      TBOX_ASSERT(CELLG == FACEG);

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(patch_geom);
      const double* dx = patch_geom->getDx();

      hier::Box pbox = patch.getBox();
      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();

      std::shared_ptr<pdat::CellData<double> > uval(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_uval, getDataContext())));
      std::shared_ptr<pdat::FaceData<double> > flux(
         SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
            patch.getPatchData(d_flux, getDataContext())));

      /*
       * Verify that the integrator providing the context correctly
       * created it, and that the ghost cell width associated with the
       * context matches the ghosts defined in this class...
       */
      TBOX_ASSERT(uval);
      TBOX_ASSERT(flux);
      TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
      TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

      /*
       * Allocate patch data for temporaries local to this routine.
       */
      pdat::FaceData<double> traced_left(pbox, 1, d_nghosts);
      pdat::FaceData<double> traced_right(pbox, 1, d_nghosts);

      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(inittraceflux2d, INITTRACEFLUX2D) (ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            uval->getPointer(),
            traced_left.getPointer(0),
            traced_left.getPointer(1),
            traced_right.getPointer(0),
            traced_right.getPointer(1),
            flux->getPointer(0),
            flux->getPointer(1)
            );
      }

      if (d_godunov_order > 1) {

         /*
          * Prepare temporary data for characteristic tracing.
          */
         int Mcells = 0;
         for (tbox::Dimension::dir_t k = 0; k < d_dim.getValue(); ++k) {
            Mcells = tbox::MathUtilities<int>::Max(Mcells, pbox.numberCells(k));
         }

// Face-centered temporary arrays
         std::vector<double> ttedgslp(2 * FACEG + 1 + Mcells);
         std::vector<double> ttraclft(2 * FACEG + 1 + Mcells);
         std::vector<double> ttracrgt(2 * FACEG + 1 + Mcells);

// Cell-centered temporary arrays
         std::vector<double> ttcelslp(2 * CELLG + Mcells);

/*
 *  Apply characteristic tracing to compute initial estimate of
 *  traces w^L and w^R at faces.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: w^L, w^R
 */
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(chartracing2d0, CHARTRACING2D0) (dt,
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               Mcells, dx[0], d_advection_velocity[0], d_godunov_order,
               uval->getPointer(),
               traced_left.getPointer(0),
               traced_right.getPointer(0),
               &ttcelslp[0],
               &ttedgslp[0],
               &ttraclft[0],
               &ttracrgt[0]);
         }

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(chartracing2d1, CHARTRACING2D1) (dt,
               ifirst(0), ilast(0), ifirst(1), ilast(1),
               Mcells, dx[1], d_advection_velocity[1], d_godunov_order,
               uval->getPointer(),
               traced_left.getPointer(1),
               traced_right.getPointer(1),
               &ttcelslp[0],
               &ttedgslp[0],
               &ttraclft[0],
               &ttracrgt[0]);
         }

      }  // if (d_godunov_order > 1) ...

      if (d_dim == tbox::Dimension(2)) {
/*
 *  Compute fluxes at faces using the face states computed so far.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: F (flux)
 */
// SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,1,dx, to get artificial viscosity
// SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,0,dx, to get NO artificial viscosity

         SAMRAI_F77_FUNC(fluxcalculation2d, FLUXCALCULATION2D) (dt, 1, 0, dx,
            ifirst(0), ilast(0), ifirst(1), ilast(1),
            d_advection_velocity,
            uval->getPointer(),
            flux->getPointer(0),
            flux->getPointer(1),
            traced_left.getPointer(0),
            traced_left.getPointer(1),
            traced_right.getPointer(0),
            traced_right.getPointer(1));

/*
 *  Re-compute traces at cell faces with transverse correction applied.
 *  Inputs: F (flux)
 *  Output: w^L, w^R (traced_left/right)
 */
         SAMRAI_F77_FUNC(fluxcorrec2d, FLUXCORREC2D) (dt, ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            dx, d_advection_velocity,
            uval->getPointer(),
            flux->getPointer(0),
            flux->getPointer(1),
            traced_left.getPointer(0),
            traced_left.getPointer(1),
            traced_right.getPointer(0),
            traced_right.getPointer(1));

/*
 *  Re-compute fluxes with updated traces.
 *  Inputs: w^L, w^R (traced_left/right)
 *  Output: F (flux)
 */
         SAMRAI_F77_FUNC(fluxcalculation2d, FLUXCALCULATION2D) (dt, 0, 0, dx,
            ifirst(0), ilast(0), ifirst(1), ilast(1),
            d_advection_velocity,
            uval->getPointer(),
            flux->getPointer(0),
            flux->getPointer(1),
            traced_left.getPointer(0),
            traced_left.getPointer(1),
            traced_right.getPointer(0),
            traced_right.getPointer(1));

      }

//     tbox::plog << "flux values: option1...." << std::endl;
//     flux->print(pbox, tbox::plog);
   }
}

/*
 *************************************************************************
 *
 * Compute numerical approximations to flux terms using an extension
 * to three dimensions of Collella's corner transport upwind approach.
 * I.E. input value corner_transport = CORNER_TRANSPORT_1
 *
 *************************************************************************
 */
void LinAdv::compute3DFluxesWithCornerTransport1(
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
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   TBOX_ASSERT(uval);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::FaceData<double> temp_traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<double> temp_traced_right(pbox, 1, d_nghosts);

   SAMRAI_F77_FUNC(inittraceflux3d, INITTRACEFLUX3D) (
      ifirst(0), ilast(0),
      ifirst(1), ilast(1),
      ifirst(2), ilast(2),
      uval->getPointer(),
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
      std::vector<double> ttedgslp(2 * FACEG + 1 + Mcells);
      std::vector<double> ttraclft(2 * FACEG + 1 + Mcells);
      std::vector<double> ttracrgt(2 * FACEG + 1 + Mcells);

      // Cell-centered temporary arrays
      std::vector<double> ttcelslp(2 * CELLG + Mcells);

      /*
       *  Apply characteristic tracing to compute initial estimate of
       *  traces w^L and w^R at faces.
       *  Inputs: w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      SAMRAI_F77_FUNC(chartracing3d0, CHARTRACING3D0) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[0], d_advection_velocity[0], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(0),
         traced_right.getPointer(0),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d1, CHARTRACING3D1) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[1], d_advection_velocity[1], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(1),
         traced_right.getPointer(1),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d2, CHARTRACING3D2) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[2], d_advection_velocity[2], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(2),
         traced_right.getPointer(2),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);
   }

   /*
    *  Compute preliminary fluxes at faces using the face states computed
    *  so far.
    *  Inputs: w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */

//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,1,dx,  to do artificial viscosity
//  SAMRAI_F77_FUNC(fluxcalculation,FLUXCALCULATION)(dt,*,*,0,dx,  to do NO artificial viscosity
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_advection_velocity,
      uval->getPointer(),
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
    *  Inputs: F (flux), w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   SAMRAI_F77_FUNC(fluxcorrec3d2d, FLUXCORREC3D2D) (
      dt, ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      dx, d_advection_velocity, 1,
      uval->getPointer(),
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

   /*
    *  Compute fluxes with partially-corrected trace states.  Store result in
    *  temporary flux vector.
    *  Inputs: temp_traced_left/right
    *  Output: temp_flux
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 0, 1, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_advection_velocity,
      uval->getPointer(),
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
    *  Inputs: F (flux), w^L, w^R (traced_left/right)
    *  Output: temp_traced_left/right
    */
   SAMRAI_F77_FUNC(fluxcorrec3d2d, FLUXCORREC3D2D) (dt, ifirst(0), ilast(0), ifirst(1),
      ilast(1), ifirst(2), ilast(2),
      dx, d_advection_velocity, -1,
      uval->getPointer(),
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

   /*
    *  Compute final predicted fluxes with both sets of transverse flux
    *  differences included.  Store the result in regular flux vector.
    *  NOTE:  the fact that we store  these fluxes in the regular (i.e.
    *  not temporary) flux vector does NOT indicate this is the final result.
    *  Rather, the flux vector is used as a convenient storage location.
    *  Inputs: temp_traced_left/right
    *  Output: flux
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_advection_velocity,
      uval->getPointer(),
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
    *  Compute the final trace state vectors at cell faces, using transverse
    *  differences of final predicted fluxes.  Store result w^L
    *  (traced_left) and w^R (traced_right) vectors.
    *  Inputs: temp_flux, flux
    *  Output: w^L, w^R (traced_left/right)
    */
   SAMRAI_F77_FUNC(fluxcorrec3d3d, FLUXCORREC3D3D) (dt, ifirst(0), ilast(0), ifirst(1),
      ilast(1), ifirst(2), ilast(2),
      dx, d_advection_velocity,
      uval->getPointer(),
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
      d_advection_velocity,
      uval->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

//     tbox::plog << "flux values: option1...." << std::endl;
//     flux->print(pbox, tbox::plog);

}

/*
 *************************************************************************
 *
 * Compute numerical approximations to flux terms using John
 * Trangenstein's interpretation of the three-dimensional version of
 * Collella's corner transport upwind approach.
 * I.E. input value corner_transport = CORNER_TRANSPORT_2
 *
 *************************************************************************
 */
void LinAdv::compute3DFluxesWithCornerTransport2(
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
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   TBOX_ASSERT(uval);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   /*
    * Allocate patch data for temporaries local to this routine.
    */
   pdat::FaceData<double> traced_left(pbox, 1, d_nghosts);
   pdat::FaceData<double> traced_right(pbox, 1, d_nghosts);
   pdat::FaceData<double> temp_flux(pbox, 1, d_fluxghosts);
   pdat::CellData<double> third_state(pbox, 1, d_nghosts);

   /*
    *  Initialize trace fluxes (w^R and w^L) with cell-centered values.
    */
   SAMRAI_F77_FUNC(inittraceflux3d, INITTRACEFLUX3D) (ifirst(0), ilast(0),
      ifirst(1), ilast(1),
      ifirst(2), ilast(2),
      uval->getPointer(),
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
    *  Compute preliminary fluxes at faces using the face states computed
    *  so far.
    *  Inputs: w^L, w^R (traced_left/right)
    *  Output: F (flux)
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 1, 1, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_advection_velocity,
      uval->getPointer(),
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
      std::vector<double> ttedgslp(2 * FACEG + 1 + Mcells);
      std::vector<double> ttraclft(2 * FACEG + 1 + Mcells);
      std::vector<double> ttracrgt(2 * FACEG + 1 + Mcells);

      // Cell-centered temporary arrays
      std::vector<double> ttcelslp(2 * CELLG + Mcells);

      /*
       *  Apply characteristic tracing to update traces w^L and
       *  w^R at faces.
       *  Inputs: w^L, w^R (traced_left/right)
       *  Output: w^L, w^R
       */
      SAMRAI_F77_FUNC(chartracing3d0, CHARTRACING3D0) (dt,
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[0], d_advection_velocity[0], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(0),
         traced_right.getPointer(0),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d1, CHARTRACING3D1) (dt,
         ifirst(0), ilast(0), ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         Mcells, dx[1], d_advection_velocity[1], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(1),
         traced_right.getPointer(1),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);

      SAMRAI_F77_FUNC(chartracing3d2, CHARTRACING3D2) (dt,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         Mcells, dx[2], d_advection_velocity[2], d_godunov_order,
         uval->getPointer(),
         traced_left.getPointer(2),
         traced_right.getPointer(2),
         &ttcelslp[0],
         &ttedgslp[0],
         &ttraclft[0],
         &ttracrgt[0]);

   } //  if (d_godunov_order > 1) ...

   for (int idir = 0; idir < d_dim.getValue(); ++idir) {

      /*
       *    Approximate traces at cell centers (in idir direction) - denoted
       *    1/3 state.
       *    Inputs:  F (flux)
       *    Output:  third_state
       */
      SAMRAI_F77_FUNC(onethirdstate3d, ONETHIRDSTATE3D) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_advection_velocity,
         uval->getPointer(),
         flux->getPointer(0),
         flux->getPointer(1),
         flux->getPointer(2),
         third_state.getPointer());
      /*
       *    Compute fluxes using 1/3 state traces, in the two directions OTHER
       *    than idir.
       *    Inputs:  third_state
       *    Output:  temp_flux (only two directions (i.e. those other than idir)
       *             are modified)
       */
      SAMRAI_F77_FUNC(fluxthird3d, FLUXTHIRD3D) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_advection_velocity,
         uval->getPointer(),
         third_state.getPointer(),
         temp_flux.getPointer(0),
         temp_flux.getPointer(1),
         temp_flux.getPointer(2));

      /*
       *    Compute transverse corrections for the traces in the two directions
       *    (OTHER than idir) using the differenced fluxes computed in those
       *    directions.
       *    Inputs:  temp_flux
       *    Output:  w^L, w^R (traced_left/right)
       */
      SAMRAI_F77_FUNC(fluxcorrecjt3d, FLUXCORRECJT3D) (dt, dx, idir,
         ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
         d_advection_velocity,
         uval->getPointer(),
         temp_flux.getPointer(0),
         temp_flux.getPointer(1),
         temp_flux.getPointer(2),
         traced_left.getPointer(0),
         traced_left.getPointer(1),
         traced_left.getPointer(2),
         traced_right.getPointer(0),
         traced_right.getPointer(1),
         traced_right.getPointer(2));

   } // loop over directions...

   /*
    *  Final flux calculation using corrected trace states.
    *  Inputs:  w^L, w^R (traced_left/right)
    *  Output:  F (flux)
    */
   SAMRAI_F77_FUNC(fluxcalculation3d, FLUXCALCULATION3D) (dt, 0, 0, 0, dx,
      ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
      d_advection_velocity,
      uval->getPointer(),
      flux->getPointer(0),
      flux->getPointer(1),
      flux->getPointer(2),
      traced_left.getPointer(0),
      traced_left.getPointer(1),
      traced_left.getPointer(2),
      traced_right.getPointer(0),
      traced_right.getPointer(1),
      traced_right.getPointer(2));

//     tbox::plog << "flux values: option2...." << std::endl;
//     flux->print(pbox, tbox::plog);
}

/*
 *************************************************************************
 *
 * Update solution variables by performing a conservative
 * difference with the fluxes calculated in computeFluxesOnPatch().
 *
 *************************************************************************
 */

void LinAdv::conservativeDifferenceOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt,
   bool at_syncronization)
{
   NULL_USE(time);
   NULL_USE(dt);
   NULL_USE(at_syncronization);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::FaceData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   TBOX_ASSERT(uval);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);

   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(consdiff2d, CONSDIFF2D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         dx,
         flux->getPointer(0),
         flux->getPointer(1),
         d_advection_velocity,
         uval->getPointer());
   }
   if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(consdiff3d, CONSDIFF3D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         ifirst(2), ilast(2), dx,
         flux->getPointer(0),
         flux->getPointer(1),
         flux->getPointer(2),
         d_advection_velocity,
         uval->getPointer());
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

void LinAdv::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));

   TBOX_ASSERT(uval);
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(uval->getTime() == fill_time);

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));

   for (int codim = 1; codim <= patch.getDim().getValue(); ++codim) {

      const std::vector<hier::BoundaryBox>& boundary_boxes =
         pgeom->getCodimensionBoundaries(codim);

      for (int bn = 0; bn < static_cast<int>(boundary_boxes.size()); ++bn) {

         const hier::Box fill_box =
            pgeom->getBoundaryFillBox(boundary_boxes[bn],
               patch.getBox(),
               ghost_width_to_fill);

         d_mesh_gen->computePatchData(patch, uval.get(), 0, fill_box);

      }

   }

}

/*
 *************************************************************************
 *
 * Tag cells for refinement using Richardson extrapolation.  Criteria
 * defined in input.
 *
 *************************************************************************
 */
void LinAdv::tagRichardsonExtrapolationCells(
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

   hier::Box pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<int> > tags(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(tag_index)));
   TBOX_ASSERT(tags);

   /*
    * Possible tagging criteria includes
    *    UVAL_RICHARDSON
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
      int size;
      double tol;
      bool time_allowed;

      if (ref == "UVAL_RICHARDSON") {
         coarsened_fine_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_uval, coarsened_fine));
         advanced_coarse_var =
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_uval, advanced_coarse));
         size = static_cast<int>(d_rich_tol.size());
         tol = ((error_level_number < size)
                ? d_rich_tol[error_level_number]
                : d_rich_tol[size - 1]);
         size = static_cast<int>(d_rich_time_min.size());
         double time_min = ((error_level_number < size)
                            ? d_rich_time_min[error_level_number]
                            : d_rich_time_min[size - 1]);
         size = static_cast<int>(d_rich_time_max.size());
         double time_max = ((error_level_number < size)
                            ? d_rich_time_max[error_level_number]
                            : d_rich_time_max[size - 1]);
         time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

         if (time_allowed) {

            TBOX_ASSERT(coarsened_fine_var);
            TBOX_ASSERT(advanced_coarse_var);

            /*
             * We tag wherever the global error > specified tolerance
             * (i.e. d_rich_tol).  The estimated global error is the
             * local truncation error * the approximate number of steps
             * used in the simulation.  Approximate the number of steps as:
             *
             *       steps = L / (s*deltat)
             * where
             *       L = length of problem domain
             *       s = wave speed
             *       delta t = timestep on current level
             *
             */
            const double* xdomainlo = d_grid_geometry->getXLower();
            const double* xdomainhi = d_grid_geometry->getXUpper();
            double max_length = 0.;
            double max_wave_speed = 0.;
            for (int idir = 0; idir < d_dim.getValue(); ++idir) {
               double length = xdomainhi[idir] - xdomainlo[idir];
               if (length > max_length) max_length = length;

               double wave_speed = d_advection_velocity[idir];
               if (wave_speed > max_wave_speed) max_wave_speed = wave_speed;
            }

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
                *
                */
               if (error > tol) {
                  if ((*tags)(*ic, 0)) {
                     (*tags)(*ic, 0) = RICHARDSON_ALREADY_TAGGED;
                  } else {
                     (*tags)(*ic, 0) = RICHARDSON_NEWLY_TAGGED;
                  }
               }

            }

         } // time_allowed

      } // if UVAL_RICHARDSON

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
 * Tag cells for refinement using gradient detector.  Tagging criteria
 * defined in input.
 *
 *************************************************************************
 */

void LinAdv::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(initial_error);

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
   TBOX_ASSERT(tags->getTime() == regrid_time);

   hier::Box pbox = patch.getBox();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   hier::Index ict(d_dim);

   int not_refine_tag_val = FALSE;
   int refine_tag_val = TRUE;

   /*
    * Create a set of temporary tags and set to untagged value.
    */
   std::shared_ptr<pdat::CellData<int> > temp_tags(
      new pdat::CellData<int>(pbox, 1, d_nghosts));
   temp_tags->setTime(regrid_time);
   temp_tags->fillAll(not_refine_tag_val);

   if (d_mesh_gen) {
      t_analytical_tag->start();
      d_mesh_gen->computePatchData(patch,
         0,
         tags.get(),
         patch.getBox());
      t_analytical_tag->stop();
   } else {
      /*
       * Possible tagging criteria includes
       *    UVAL_DEVIATION, UVAL_GRADIENT, UVAL_SHOCK
       * The criteria is specified over a time interval.
       *
       * Loop over criteria provided and check to make sure we are in the
       * specified time interval.  If so, apply appropriate tagging for
       * the level.
       */
      for (int ncrit = 0;
           ncrit < static_cast<int>(d_refinement_criteria.size()); ++ncrit) {

         std::string ref = d_refinement_criteria[ncrit];
         std::shared_ptr<pdat::CellData<double> > var(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(d_uval, getDataContext())));

         TBOX_ASSERT(var);

         hier::IntVector vghost(var->getGhostCellWidth());
         hier::IntVector tagghost(tags->getGhostCellWidth());

         int size = 0;
         double tol = 0.;
         double onset = 0.;
         bool time_allowed = false;

         if (ref == "UVAL_DEVIATION") {
            size = static_cast<int>(d_dev_tol.size());
            tol = ((error_level_number < size)
                   ? d_dev_tol[error_level_number]
                   : d_dev_tol[size - 1]);
            size = static_cast<int>(d_dev.size());
            double dev = ((error_level_number < size)
                          ? d_dev[error_level_number]
                          : d_dev[size - 1]);
            size = static_cast<int>(d_dev_time_min.size());
            double time_min = ((error_level_number < size)
                               ? d_dev_time_min[error_level_number]
                               : d_dev_time_min[size - 1]);
            size = static_cast<int>(d_dev_time_max.size());
            double time_max = ((error_level_number < size)
                               ? d_dev_time_max[error_level_number]
                               : d_dev_time_max[size - 1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed) {

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
                     (*temp_tags)(*ic, 0) = refine_tag_val;
                  }
               }
            }
         }

         if (ref == "UVAL_GRADIENT") {
            size = static_cast<int>(d_grad_tol.size());
            tol = ((error_level_number < size)
                   ? d_grad_tol[error_level_number]
                   : d_grad_tol[size - 1]);
            size = static_cast<int>(d_grad_time_min.size());
            double time_min = ((error_level_number < size)
                               ? d_grad_time_min[error_level_number]
                               : d_grad_time_min[size - 1]);
            size = static_cast<int>(d_grad_time_max.size());
            double time_max = ((error_level_number < size)
                               ? d_grad_time_max[error_level_number]
                               : d_grad_time_max[size - 1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed) {

               if (d_dim == tbox::Dimension(2)) {
                  SAMRAI_F77_FUNC(detectgrad2d, DETECTGRAD2D) (
                     ifirst(0), ilast(0), ifirst(1), ilast(1),
                     vghost(0), tagghost(0), d_nghosts(0),
                     vghost(1), tagghost(1), d_nghosts(1),
                     dx,
                     tol,
                     refine_tag_val, not_refine_tag_val,
                     var->getPointer(),
                     tags->getPointer(), temp_tags->getPointer());
               } else if (d_dim == tbox::Dimension(3)) {
                  SAMRAI_F77_FUNC(detectgrad3d, DETECTGRAD3D) (
                     ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2),
                     ilast(2),
                     vghost(0), tagghost(0), d_nghosts(0),
                     vghost(1), tagghost(1), d_nghosts(1),
                     vghost(2), tagghost(2), d_nghosts(2),
                     dx,
                     tol,
                     refine_tag_val, not_refine_tag_val,
                     var->getPointer(),
                     tags->getPointer(), temp_tags->getPointer());
               }
            }

         }

         if (ref == "UVAL_SHOCK") {
            size = static_cast<int>(d_shock_tol.size());
            tol = ((error_level_number < size)
                   ? d_shock_tol[error_level_number]
                   : d_shock_tol[size - 1]);
            size = static_cast<int>(d_shock_onset.size());
            onset = ((error_level_number < size)
                     ? d_shock_onset[error_level_number]
                     : d_shock_onset[size - 1]);
            size = static_cast<int>(d_shock_time_min.size());
            double time_min = ((error_level_number < size)
                               ? d_shock_time_min[error_level_number]
                               : d_shock_time_min[size - 1]);
            size = static_cast<int>(d_shock_time_max.size());
            double time_max = ((error_level_number < size)
                               ? d_shock_time_max[error_level_number]
                               : d_shock_time_max[size - 1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed) {

               if (d_dim == tbox::Dimension(2)) {
                  SAMRAI_F77_FUNC(detectshock2d, DETECTSHOCK2D) (
                     ifirst(0), ilast(0), ifirst(1), ilast(1),
                     vghost(0), tagghost(0), d_nghosts(0),
                     vghost(1), tagghost(1), d_nghosts(1),
                     dx,
                     tol,
                     onset,
                     refine_tag_val, not_refine_tag_val,
                     var->getPointer(),
                     tags->getPointer(), temp_tags->getPointer());
               } else if (d_dim == tbox::Dimension(3)) {
                  SAMRAI_F77_FUNC(detectshock3d, DETECTSHOCK3D) (
                     ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2),
                     ilast(2),
                     vghost(0), tagghost(0), d_nghosts(0),
                     vghost(1), tagghost(1), d_nghosts(1),
                     vghost(2), tagghost(2), d_nghosts(2),
                     dx,
                     tol,
                     onset,
                     refine_tag_val, not_refine_tag_val,
                     var->getPointer(),
                     tags->getPointer(), temp_tags->getPointer());
               }
            }

         }

      }  // loop over criteria

      /*
       * Adjust temp_tags from those tags set in Richardson extrapolation.
       * Here, we just reset any tags that were set in Richardson extrapolation
       * to be the designated "refine_tag_val".
       */
      if (uses_richardson_extrapolation_too) {
         pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
              ic != icend; ++ic) {
            if ((*tags)(*ic, 0) == RICHARDSON_ALREADY_TAGGED ||
                (*tags)(*ic, 0) == RICHARDSON_NEWLY_TAGGED) {
               (*temp_tags)(*ic, 0) = refine_tag_val;
            }
         }
      }

      /*
       * Update tags.
       */
      pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
      for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
           ic != icend; ++ic) {
         (*tags)(*ic, 0) = (*temp_tags)(*ic, 0);
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
void LinAdv::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);

   d_mesh_gen->registerVariablesWithPlotter(*viz_writer);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

#ifdef HAVE_HDF5
   if (viz_writer) {
      viz_writer->
      registerPlotQuantity("U",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_uval, vardb->getContext("CURRENT")));
   }
#endif
}
#endif

bool LinAdv::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id,
   double simulation_time) const
{
   NULL_USE(buffer);
   NULL_USE(patch);
   NULL_USE(region);
   NULL_USE(variable_name);
   NULL_USE(depth_id);
   NULL_USE(simulation_time);
   TBOX_ERROR("Should not be here.  This object didn't register any derived plot variables.");
   return true;
}

/*
 *************************************************************************
 *
 * Write LinAdv object state to specified stream.
 *
 *************************************************************************
 */

void LinAdv::printClassData(
   std::ostream& os) const
{
   int j;

   os << "\nLinAdv::printClassData..." << std::endl;
   os << "LinAdv: this = " << (LinAdv *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_grid_geometry = "
      << d_grid_geometry.get() << std::endl;

   os << "Parameters for numerical method ..." << std::endl;
   os << "   d_advection_velocity = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_advection_velocity[j] << " ";
   os << std::endl;
   os << "   d_godunov_order = " << d_godunov_order << std::endl;
   os << "   d_corner_transport = " << d_corner_transport << std::endl;
   os << "   d_nghosts = " << d_nghosts << std::endl;
   os << "   d_fluxghosts = " << d_fluxghosts << std::endl;
   os << "   Boundary condition data " << std::endl;

   os << "   Refinement criteria parameters " << std::endl;

   for (j = 0; j < static_cast<int>(d_refinement_criteria.size()); ++j) {
      os << "       d_refinement_criteria[" << j << "] = "
         << d_refinement_criteria[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_dev_tol.size()); ++j) {
      os << "       d_dev_tol[" << j << "] = "
         << d_dev_tol[j] << std::endl;
   }
   for (j = 0; j < static_cast<int>(d_dev.size()); ++j) {
      os << "       d_dev[" << j << "] = "
         << d_dev[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_dev_time_max.size()); ++j) {
      os << "       d_dev_time_max[" << j << "] = "
         << d_dev_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_dev_time_min.size()); ++j) {
      os << "       d_dev_time_min[" << j << "] = "
         << d_dev_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_grad_tol.size()); ++j) {
      os << "       d_grad_tol[" << j << "] = "
         << d_grad_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_grad_time_max.size()); ++j) {
      os << "       d_grad_time_max[" << j << "] = "
         << d_grad_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_grad_time_min.size()); ++j) {
      os << "       d_grad_time_min[" << j << "] = "
         << d_grad_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_shock_onset.size()); ++j) {
      os << "       d_shock_onset[" << j << "] = "
         << d_shock_onset[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_shock_tol.size()); ++j) {
      os << "       d_shock_tol[" << j << "] = "
         << d_shock_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_shock_time_max.size()); ++j) {
      os << "       d_shock_time_max[" << j << "] = "
         << d_shock_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_shock_time_min.size()); ++j) {
      os << "       d_shock_time_min[" << j << "] = "
         << d_shock_time_min[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_rich_tol.size()); ++j) {
      os << "       d_rich_tol[" << j << "] = "
         << d_rich_tol[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_rich_time_max.size()); ++j) {
      os << "       d_rich_time_max[" << j << "] = "
         << d_rich_time_max[j] << std::endl;
   }
   os << std::endl;
   for (j = 0; j < static_cast<int>(d_rich_time_min.size()); ++j) {
      os << "       d_rich_time_min[" << j << "] = "
         << d_rich_time_min[j] << std::endl;
   }
   os << std::endl;

}

/*
 *************************************************************************
 *
 * Read data members from input.  All values set from restart can be
 * overridden by values in the input database.
 *
 *************************************************************************
 */
void LinAdv::getFromInput(
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

   if (input_db->keyExists("advection_velocity")) {
      input_db->getDoubleArray("advection_velocity",
         d_advection_velocity, d_dim.getValue());
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `advection_velocity' not found in input.");
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
         d_refinement_criteria =
            refine_db->getStringVector("refine_criteria");
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

            if (!(error_key == "UVAL_DEVIATION" ||
                  error_key == "UVAL_GRADIENT" ||
                  error_key == "UVAL_SHOCK" ||
                  error_key == "UVAL_RICHARDSON")) {
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

            if (error_db && error_key == "UVAL_DEVIATION") {

               if (error_db->keyExists("dev_tol")) {
                  d_dev_tol = error_db->getDoubleVector("dev_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `dev_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("uval_dev")) {
                  d_dev = error_db->getDoubleVector("uval_dev");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `uval_dev' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_dev_time_max = error_db->getDoubleVector("time_max");
               } else {
                  d_dev_time_max.resize(1);
                  d_dev_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_dev_time_min = error_db->getDoubleVector("time_min");
               } else {
                  d_dev_time_min.resize(1);
                  d_dev_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "UVAL_GRADIENT") {

               if (error_db->keyExists("grad_tol")) {
                  d_grad_tol = error_db->getDoubleVector("grad_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `grad_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_grad_time_max = error_db->getDoubleVector("time_max");
               } else {
                  d_grad_time_max.resize(1);
                  d_grad_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_grad_time_min = error_db->getDoubleVector("time_min");
               } else {
                  d_grad_time_min.resize(1);
                  d_grad_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "UVAL_SHOCK") {

               if (error_db->keyExists("shock_onset")) {
                  d_shock_onset = error_db->getDoubleVector("shock_onset");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_onset' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("shock_tol")) {
                  d_shock_tol = error_db->getDoubleVector("shock_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `shock_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_shock_time_max = error_db->getDoubleVector("time_max");
               } else {
                  d_shock_time_max.resize(1);
                  d_shock_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_shock_time_min = error_db->getDoubleVector("time_min");
               } else {
                  d_shock_time_min.resize(1);
                  d_shock_time_min[0] = 0.;
               }

            }

            if (error_db && error_key == "UVAL_RICHARDSON") {

               if (error_db->keyExists("rich_tol")) {
                  d_rich_tol = error_db->getDoubleVector("rich_tol");
               } else {
                  TBOX_ERROR(
                     d_object_name << ": "
                                   << "No key `rich_tol' found in data for "
                                   << error_key << std::endl);
               }

               if (error_db->keyExists("time_max")) {
                  d_rich_time_max = error_db->getDoubleVector("time_max");
               } else {
                  d_rich_time_max.resize(1);
                  d_rich_time_max[0] = tbox::MathUtilities<double>::getMax();
               }

               if (error_db->keyExists("time_min")) {
                  d_rich_time_min = error_db->getDoubleVector("time_min");
               } else {
                  d_rich_time_min.resize(1);
                  d_rich_time_min[0] = 0.;
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

   } // refine db entry exists

   hier::IntVector periodic(
      d_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   int num_per_dirs = 0;
   for (int id = 0; id < d_dim.getValue(); ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

}

/*
 *************************************************************************
 *
 * Routines to put/get data members to/from restart database.
 *
 *************************************************************************
 */

void LinAdv::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("LINADV_VERSION", LINADV_VERSION);

   restart_db->putDoubleArray("d_advection_velocity",
      d_advection_velocity,
      d_dim.getValue());

   restart_db->putInteger("d_godunov_order", d_godunov_order);
   restart_db->putString("d_corner_transport", d_corner_transport);
   restart_db->putIntegerArray("d_nghosts", &d_nghosts[0], d_dim.getValue());
   restart_db->putIntegerArray("d_fluxghosts",
      &d_fluxghosts[0],
      d_dim.getValue());

   if (d_refinement_criteria.size() > 0) {
      restart_db->putStringVector("d_refinement_criteria",
         d_refinement_criteria);
   }
   for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); ++i) {

      if (d_refinement_criteria[i] == "UVAL_DEVIATION") {
         restart_db->putDoubleVector("d_dev_tol", d_dev_tol);
         restart_db->putDoubleVector("d_dev", d_dev);
         restart_db->putDoubleVector("d_dev_time_max", d_dev_time_max);
         restart_db->putDoubleVector("d_dev_time_min", d_dev_time_min);
      } else if (d_refinement_criteria[i] == "UVAL_GRADIENT") {
         restart_db->putDoubleVector("d_grad_tol", d_grad_tol);
         restart_db->putDoubleVector("d_grad_time_max", d_grad_time_max);
         restart_db->putDoubleVector("d_grad_time_min", d_grad_time_min);
      } else if (d_refinement_criteria[i] == "UVAL_SHOCK") {
         restart_db->putDoubleVector("d_shock_onset", d_shock_onset);
         restart_db->putDoubleVector("d_shock_tol", d_shock_tol);
         restart_db->putDoubleVector("d_shock_time_max", d_shock_time_max);
         restart_db->putDoubleVector("d_shock_time_min", d_shock_time_min);
      } else if (d_refinement_criteria[i] == "UVAL_RICHARDSON") {
         restart_db->putDoubleVector("d_rich_tol", d_rich_tol);
         restart_db->putDoubleVector("d_rich_time_max", d_rich_time_max);
         restart_db->putDoubleVector("d_rich_time_min", d_rich_time_min);
      }

   }

}

/*
 *************************************************************************
 *
 *    Access class information from restart database.
 *
 *************************************************************************
 */
void LinAdv::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file.");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("LINADV_VERSION");
   if (ver != LINADV_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   db->getDoubleArray("d_advection_velocity", d_advection_velocity, d_dim.getValue());

   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

   int* tmp_nghosts = &d_nghosts[0];
   db->getIntegerArray("d_nghosts", tmp_nghosts, d_dim.getValue());
   if (!(d_nghosts == CELLG)) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Key data `d_nghosts' in restart file != CELLG." << std::endl);
   }
   int* tmp_fluxghosts = &d_fluxghosts[0];
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, d_dim.getValue());
   if (!(d_fluxghosts == FLUXG)) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Key data `d_fluxghosts' in restart file != FLUXG." << std::endl);
   }

   if (db->keyExists("d_refinement_criteria")) {
      d_refinement_criteria = db->getStringVector("d_refinement_criteria");
   }
   for (int i = 0; i < static_cast<int>(d_refinement_criteria.size()); ++i) {

      if (d_refinement_criteria[i] == "UVAL_DEVIATION") {
         d_dev_tol = db->getDoubleVector("d_dev_tol");
         d_dev_time_max = db->getDoubleVector("d_dev_time_max");
         d_dev_time_min = db->getDoubleVector("d_dev_time_min");
      } else if (d_refinement_criteria[i] == "UVAL_GRADIENT") {
         d_grad_tol = db->getDoubleVector("d_grad_tol");
         d_grad_time_max = db->getDoubleVector("d_grad_time_max");
         d_grad_time_min = db->getDoubleVector("d_grad_time_min");
      } else if (d_refinement_criteria[i] == "UVAL_SHOCK") {
         d_shock_onset = db->getDoubleVector("d_shock_onset");
         d_shock_tol = db->getDoubleVector("d_shock_tol");
         d_shock_time_max = db->getDoubleVector("d_shock_time_max");
         d_shock_time_min = db->getDoubleVector("d_shock_time_min");
      } else if (d_refinement_criteria[i] == "UVAL_RICHARDSON") {
         d_rich_tol = db->getDoubleVector("d_rich_tol");
         d_rich_time_max = db->getDoubleVector("d_rich_time_max");
         d_rich_time_min = db->getDoubleVector("d_rich_time_min");
      }

   }

}

void LinAdv::readStateDataEntry(
   std::shared_ptr<tbox::Database> db,
   const std::string& db_name,
   int array_indx,
   std::vector<double>& uval)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(static_cast<int>(uval.size()) > array_indx);

   if (db->keyExists("uval")) {
      uval[array_indx] = db->getDouble("uval");
   } else {
      TBOX_ERROR(d_object_name << ": "
                               << "`uval' entry missing from " << db_name
                               << " input database. " << std::endl);
   }

}
