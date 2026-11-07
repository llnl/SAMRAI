/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in linear advection ex.
 *
 ************************************************************************/
#include "MblkLinAdv.h"

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

#include <vector>


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeDoubleInjection.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (0)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "test/testlib/SkeletonBoundaryUtilities2.h"
#include "test/testlib/SkeletonBoundaryUtilities3.h"

// Depth of the advected variable
#define DEPTH (1)

// Number of ghosts cells used for each variable quantity
#define CELLG (1)
#define FLUXG (0)
#define NODEG (0)

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

// Version of MblkLinAdv restart file data
#define MBLKLINADV_VERSION (3)

//
// some extra defines for C code
//
#define real8 double
#define POLY3(i, j, k, imin, jmin, kmin, nx, \
              nxny) ((i - imin) + (j - jmin) * (nx) + (k - kmin) * (nxny))
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

//
// inline geometry functions
//
#include "test/testlib/GeomUtilsAMR.h"

/*
 *************************************************************************
 *
 * The constructor for MblkLinAdv class sets data members to defualt values,
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

MblkLinAdv::MblkLinAdv(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<hier::BaseGridGeometry>& grid_geom):
   MblkHyperbolicPatchStrategy(dim),
   d_object_name(object_name),
   d_dim(dim),
   d_grid_geometry(grid_geom),
   d_use_nonuniform_workload(false),
   d_uval(new pdat::CellVariable<double>(dim, "uval", DEPTH)),
   d_vol(new pdat::CellVariable<double>(dim, "vol", 1)),
   d_flux(new pdat::SideVariable<double>(dim, "flux",
                                         hier::IntVector::getOne(dim), 1)),
   d_xyz(new pdat::NodeVariable<double>(dim, "xyz", dim.getValue())),
   d_dx_set(false),
   d_godunov_order(1),
   d_corner_transport("CORNER_TRANSPORT_1"),
   d_nghosts(hier::IntVector(dim, CELLG)),
   d_fluxghosts(hier::IntVector(dim, FLUXG)),
   d_nodeghosts(hier::IntVector(dim, NODEG)),
   d_data_problem_int(tbox::MathUtilities<int>::getMax()),
   d_radius(tbox::MathUtilities<double>::getSignalingNaN()),
   d_uval_inside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_uval_outside(tbox::MathUtilities<double>::getSignalingNaN()),
   d_number_of_intervals(0),
   d_amplitude(0.)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   /*
    * Setup MblkGeometry object to manage construction of mapped grids
    */
   d_mblk_geometry = new MblkGeometry("MblkGeometry",
         d_dim,
         input_db,
         grid_geom->getNumberBlocks());


   // SPHERE problem...
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, d_dim.getValue());

   // SINE problem
   for (int k = 0; k < d_dim.getValue(); ++k) d_frequency[k] = 0.;

   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

   if (d_dim == tbox::Dimension(2)) {
      d_scalar_bdry_edge_conds.resize(NUM_2D_EDGES);
      for (int ei = 0; ei < NUM_2D_EDGES; ++ei) {
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      }

      d_scalar_bdry_node_conds.resize(NUM_2D_NODES);
      d_node_bdry_edge.resize(NUM_2D_NODES);

      for (int ni = 0; ni < NUM_2D_NODES; ++ni) {
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
      }

      d_bdry_edge_uval.resize(NUM_2D_EDGES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_uval);
   } else if (d_dim == tbox::Dimension(3)) {
      d_scalar_bdry_face_conds.resize(NUM_3D_FACES);
      for (int fi = 0; fi < NUM_3D_FACES; ++fi) {
         d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      }

      d_scalar_bdry_edge_conds.resize(NUM_3D_EDGES);
      d_edge_bdry_face.resize(NUM_3D_EDGES);
      for (int ei = 0; ei < NUM_3D_EDGES; ++ei) {
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
      }

      d_scalar_bdry_node_conds.resize(NUM_3D_NODES);
      d_node_bdry_face.resize(NUM_3D_NODES);

      for (int ni = 0; ni < NUM_3D_NODES; ++ni) {
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
      }

      d_bdry_face_uval.resize(NUM_3D_FACES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_uval);
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

   if (d_data_problem == "PIECEWISE_CONSTANT_X") {
      d_data_problem_int = PIECEWISE_CONSTANT_X;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
      d_data_problem_int = PIECEWISE_CONSTANT_Y;
   } else if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
      d_data_problem_int = PIECEWISE_CONSTANT_Z;
   } else if (d_data_problem == "SINE_CONSTANT_X") {
      d_data_problem_int = SINE_CONSTANT_X;
   } else if (d_data_problem == "SINE_CONSTANT_Y") {
      d_data_problem_int = SINE_CONSTANT_Y;
   } else if (d_data_problem == "SINE_CONSTANT_Z") {
      d_data_problem_int = SINE_CONSTANT_Z;
   } else if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
   } else {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Unknown d_data_problem string = "
                       << d_data_problem
                       << " encountered in constructor" << std::endl);
   }

   /*
    * Postprocess boundary data from input/restart values.  Note: scalar
    * quantity in this problem cannot have reflective boundary conditions
    * so we reset them to FLOW.
    */
   if (d_dim == tbox::Dimension(2)) {
      for (int i = 0; i < NUM_2D_EDGES; ++i) {
         if (d_scalar_bdry_edge_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::FLOW;
         }
      }

      for (int i = 0; i < NUM_2D_NODES; ++i) {
         if (d_scalar_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_scalar_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }

         if (d_scalar_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_edge[i] =
               SkeletonBoundaryUtilities2::getEdgeLocationForNodeBdry(
                  i, d_scalar_bdry_node_conds[i]);
         }
      }
   } else if (d_dim == tbox::Dimension(3)) {
      for (int i = 0; i < NUM_3D_FACES; ++i) {
         if (d_scalar_bdry_face_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_face_conds[i] = BdryCond::FLOW;
         }
      }

      for (int i = 0; i < NUM_3D_EDGES; ++i) {
         if (d_scalar_bdry_edge_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::XFLOW;
         }
         if (d_scalar_bdry_edge_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::YFLOW;
         }
         if (d_scalar_bdry_edge_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::ZFLOW;
         }

         if (d_scalar_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
            d_edge_bdry_face[i] =
               SkeletonBoundaryUtilities3::getFaceLocationForEdgeBdry(
                  i, d_scalar_bdry_edge_conds[i]);
         }
      }

      for (int i = 0; i < NUM_3D_NODES; ++i) {
         if (d_scalar_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_scalar_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }
         if (d_scalar_bdry_node_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::ZFLOW;
         }

         if (d_scalar_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_face[i] =
               SkeletonBoundaryUtilities3::getFaceLocationForNodeBdry(
                  i, d_scalar_bdry_node_conds[i]);
         }
      }

   }

}

/*
 *************************************************************************
 *
 * Empty destructor for MblkLinAdv class.
 *
 *************************************************************************
 */

MblkLinAdv::~MblkLinAdv() {
   if (d_mblk_geometry) delete d_mblk_geometry;
}

/*
 *************************************************************************
 *
 * Register conserved variable (u) (i.e., solution state variable) and
 * flux variable with hyperbolic integrator that manages storage for
 * those quantities.  Also, register plot data with VisIt.
 *
 *************************************************************************
 */

void MblkLinAdv::registerModelVariables(
   MblkHyperbolicLevelIntegrator* integrator)
{
   TBOX_ASSERT(integrator != 0);

   d_cell_cons_linear_refine_op.reset(
      new SkeletonCellDoubleConservativeLinearRefine(d_dim));
   d_cell_cons_coarsen_op.reset(
      new SkeletonCellDoubleWeightedAverage(d_dim));
   d_cell_time_interp_op.reset(
      new pdat::CellDoubleLinearTimeInterpolateOp());
   d_side_cons_coarsen_op.reset(
      new SkeletonOutersideDoubleWeightedAverage(d_dim));

   // Note that the Node linear refine operator is null for this case,
   // which is OK because the only node data is the grid coordinates,
   // which we explicitly set on any new patch
   std::shared_ptr<hier::RefineOperator> node_linear_refine_op;
   std::shared_ptr<pdat::NodeDoubleInjection> node_cons_coarsen_op(
      new pdat::NodeDoubleInjection());

   integrator->registerVariable(d_uval, d_nghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      d_cell_cons_coarsen_op,
      d_cell_cons_linear_refine_op,
      d_cell_time_interp_op);

   integrator->registerVariable(d_vol, d_nghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      d_cell_cons_coarsen_op,
      d_cell_cons_linear_refine_op,
      d_cell_time_interp_op);

   integrator->registerVariable(d_flux, d_fluxghosts,
      MblkHyperbolicLevelIntegrator::FLUX,
      d_side_cons_coarsen_op);

   std::shared_ptr<hier::TimeInterpolateOperator> node_time_interp_op(
      new pdat::NodeDoubleLinearTimeInterpolateOp());
   integrator->registerVariable(d_xyz, d_nodeghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      node_cons_coarsen_op,
      node_linear_refine_op,
      node_time_interp_op);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

#ifdef HAVE_HDF5
   if (d_visit_writer) {
      d_visit_writer->
      registerPlotQuantity("U",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_uval, integrator->getPlotContext()));
      d_visit_writer->
      registerPlotQuantity("vol",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_vol, integrator->getPlotContext()));
      d_visit_writer->
      registerNodeCoordinates(
         vardb->mapVariableAndContextToIndex(
            d_xyz, integrator->getPlotContext()));
   }

   if (!d_visit_writer) {
      TBOX_WARNING(
         d_object_name << ": registerModelVariables()"
                       << "\nVisIt data writer was"
                       << "\nregistered.  Consequently, no plot data will"
                       << "\nbe written." << std::endl);
   }
#endif

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
void MblkLinAdv::initializeDataOnPatch(
   hier::Patch& patch,
   const double data_time,
   const bool initial_time)
{
   NULL_USE(data_time);
   /*
    * Build the mapped grid on the patch.
    */
   hier::BlockId::block_t block_number =
      patch.getBox().getBlockId().getBlockValue();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   /*
    * Set the dx in the operators
    */
   double dx[SAMRAI::MAX_DIM_VAL];
   d_mblk_geometry->getDx(level_number, dx);
   d_dx_set = true;

   d_cell_cons_linear_refine_op->setDx(level_number, dx);
   d_cell_cons_coarsen_op->setDx(level_number, dx);
   d_side_cons_coarsen_op->setDx(level_number, dx);

   if (initial_time) {

      std::shared_ptr<pdat::CellData<double> > uval(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_uval, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > vol(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_vol, getDataContext())));
      std::shared_ptr<pdat::NodeData<double> > xyz(
         SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
            patch.getPatchData(d_xyz, getDataContext())));

      TBOX_ASSERT(uval);
      TBOX_ASSERT(vol);
      TBOX_ASSERT(xyz);
      TBOX_ASSERT(uval->getGhostCellWidth() == vol->getGhostCellWidth());

      hier::IntVector uval_ghosts = uval->getGhostCellWidth();
      hier::IntVector xyz_ghosts = xyz->getGhostCellWidth();

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();

      int imin = ifirst(0) - uval_ghosts(0);
      int imax = ilast(0) + uval_ghosts(0);
      int jmin = ifirst(1) - uval_ghosts(1);
      int jmax = ilast(1) + uval_ghosts(1);
      int kmin = ifirst(2) - uval_ghosts(2);
      //int kmax = ilast(2)  + uval_ghosts(2);
      int nx = imax - imin + 1;
      int ny = jmax - jmin + 1;
      // int nz   = kmax - kmin + 1;
      int nxny = nx * ny;
      // int nel  = nx*ny*nz;

      int nd_imin = ifirst(0) - xyz_ghosts(0);
      int nd_imax = ilast(0) + 1 + xyz_ghosts(0);
      int nd_jmin = ifirst(1) - xyz_ghosts(1);
      int nd_jmax = ilast(1) + 1 + xyz_ghosts(1);
      int nd_kmin = ifirst(2) - xyz_ghosts(2);
      //int nd_kmax = ilast(2)  + 1 + xyz_ghosts(2);
      int nd_nx = nd_imax - nd_imin + 1;
      int nd_ny = nd_jmax - nd_jmin + 1;
      // int nd_nz   = nd_kmax - nd_kmin + 1;
      int nd_nxny = nd_nx * nd_ny;
      // int nd_nel  = nd_nx*nd_ny*nd_nz;

      //
      // get the pointers
      //
      double* psi = uval->getPointer();
      double* cvol = vol->getPointer();

      double* x = xyz->getPointer(0);
      double* y = xyz->getPointer(1);
      double* z = xyz->getPointer(2);

      //
      // compute the source due to the upwind method
      //
      for (int k = ifirst(2); k <= ilast(2); ++k) {
         for (int j = ifirst(1); j <= ilast(1); ++j) {
            for (int i = ifirst(0); i <= ilast(0); ++i) {

               int cind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

               int n1 = POLY3(i,
                     j,
                     k,
                     nd_imin,
                     nd_jmin,
                     nd_kmin,
                     nd_nx,
                     nd_nxny);
               int n2 = n1 + 1;
               int n3 = n1 + 1 + nd_nx;
               int n4 = n1 + nd_nx;

               int n5 = n1 + nd_nxny;
               int n6 = n1 + nd_nxny + 1;
               int n7 = n1 + nd_nxny + 1 + nd_nx;
               int n8 = n1 + nd_nxny + nd_nx;

               real8 vol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                     x[n5], x[n6], x[n7], x[n8],

                     y[n1], y[n2], y[n3], y[n4],
                     y[n5], y[n6], y[n7], y[n8],

                     z[n1], z[n2], z[n3], z[n4],
                     z[n5], z[n6], z[n7], z[n8]);

               cvol[cind] = vol;

               real8 xmid = 0.125 * (x[n1] + x[n2] + x[n3] + x[n4]
                                     + x[n5] + x[n6] + x[n7] + x[n8]);

               real8 ymid = 0.125 * (y[n1] + y[n2] + y[n3] + y[n4]
                                     + y[n5] + y[n6] + y[n7] + y[n8]);

               real8 zmid = 0.125 * (z[n1] + z[n2] + z[n3] + z[n4]
                                     + z[n5] + z[n6] + z[n7] + z[n8]);

               real8 xc = xmid - d_center[0];
               real8 yc = ymid - d_center[1];
               real8 zc = zmid - d_center[2];

               real8 radsq = xc * xc + yc * yc + zc * zc;

               bool inside = (d_radius * d_radius) > radsq;

               if (inside) {
                  psi[cind] = d_uval_inside;
               } else {
                  psi[cind] = d_uval_outside;
               }

            }
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
   }

}

//
// ==================================================================
//    extra code
// ==================================================================
//

/*
 *************************************************************************
 *
 * Compute stable time increment for patch.  Return this value.
 *
 *************************************************************************
 */

double MblkLinAdv::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time,
   const double dt_time)
{
   NULL_USE(initial_time);
   NULL_USE(dt_time);

   /*
    * Build the mapped grid on the patch.
    */
   hier::BlockId::block_t block_number =
      patch.getBox().getBlockId().getBlockValue();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > vol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_vol, getDataContext())));
   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));

   TBOX_ASSERT(uval);
   TBOX_ASSERT(vol);
   TBOX_ASSERT(xyz);
   TBOX_ASSERT(uval->getGhostCellWidth() == vol->getGhostCellWidth());

   hier::IntVector uval_ghosts = uval->getGhostCellWidth();
   hier::IntVector xyz_ghosts = xyz->getGhostCellWidth();

   /*
    * Adjust advection velocity based on rotation of block
    */
   //int block_number = getBlockNumber();

   int imin = ifirst(0) - uval_ghosts(0);
   int imax = ilast(0) + uval_ghosts(0);
   int jmin = ifirst(1) - uval_ghosts(1);
   int jmax = ilast(1) + uval_ghosts(1);
   int kmin = ifirst(2) - uval_ghosts(2);
   //int kmax = ilast(2)  + uval_ghosts(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   // int nz   = kmax - kmin + 1;
   int nxny = nx * ny;
   // int nel  = nx*ny*nz;

   int nd_imin = ifirst(0) - xyz_ghosts(0);
   int nd_imax = ilast(0) + 1 + xyz_ghosts(0);
   int nd_jmin = ifirst(1) - xyz_ghosts(1);
   int nd_jmax = ilast(1) + 1 + xyz_ghosts(1);
   int nd_kmin = ifirst(2) - xyz_ghosts(2);
   //int nd_kmax = ilast(2)  + 1 + xyz_ghosts(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   // int nd_nz   = nd_kmax - nd_kmin + 1;
   int nd_nxny = nd_nx * nd_ny;
   // int nd_nel  = nd_nx*nd_ny*nd_nz;

   //int rot = d_mblk_geometry->getBlockRotation(block_number);

   real8 u = d_advection_velocity[0];
   real8 v = d_advection_velocity[1];
   real8 w = d_advection_velocity[2];

   //
   // compute the source due to the upwind method
   //
   double stabdt = 1.e20;

   double* cvol = vol->getPointer();

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int cind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

            //   have set up the normal so that it always points in the positive, i, j, and k directions

            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;
            int n7 = n1 + nd_nxny + 1 + nd_nx;
            int n8 = n1 + nd_nxny + nd_nx;

            real8 vol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                  x[n5], x[n6], x[n7], x[n8],

                  y[n1], y[n2], y[n3], y[n4],
                  y[n5], y[n6], y[n7], y[n8],

                  z[n1], z[n2], z[n3], z[n4],
                  z[n5], z[n6], z[n7], z[n8]);

            cvol[cind] = vol;

            if (vol < 0.) {
               TBOX_ERROR("Error:  negative volume computed in UpwindVolume");
            }

            real8 xx[8];
            real8 yy[8];
            real8 zz[8];

            xx[0] = x[n1];
            xx[1] = x[n2];
            xx[2] = x[n3];
            xx[3] = x[n4];
            xx[4] = x[n5];
            xx[5] = x[n6];
            xx[6] = x[n7];
            xx[7] = x[n8];

            yy[0] = y[n1];
            yy[1] = y[n2];
            yy[2] = y[n3];
            yy[3] = y[n4];
            yy[4] = y[n5];
            yy[5] = y[n6];
            yy[6] = y[n7];
            yy[7] = y[n8];

            zz[0] = z[n1];
            zz[1] = z[n2];
            zz[2] = z[n3];
            zz[3] = z[n4];
            zz[4] = z[n5];
            zz[5] = z[n6];
            zz[6] = z[n7];
            zz[7] = z[n8];

            real8 arealg = UpwindCharacteristicLength(xx, yy, zz, vol);

            real8 uu = MAX(u, MAX(v, w));

            stabdt = MIN(stabdt, arealg / uu);

         }
      }
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

void MblkLinAdv::computeFluxesOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt)
{
   NULL_USE(time);

   //return;

   hier::BlockId::block_t block_number =
      patch.getBox().getBlockId().getBlockValue();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > uval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > vol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::SideData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));

   /*
    * Verify that the integrator providing the context correctly
    * created it, and that the ghost cell width associated with the
    * context matches the ghosts defined in this class...
    */
   TBOX_ASSERT(uval);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(vol);
   TBOX_ASSERT(xyz);
   TBOX_ASSERT(uval->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(vol->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(xyz->getGhostCellWidth() == d_nodeghosts);

   //
   // ------------------------------- spliced in code ----------------------------
   //

   int fx_imn = ifirst(0) - d_fluxghosts(0);
   int fx_imx = ilast(0) + 1 + d_fluxghosts(0);
   int fx_jmn = ifirst(1) - d_fluxghosts(1);
   int fx_jmx = ilast(1) + d_fluxghosts(1);
   int fx_kmn = ifirst(2) - d_fluxghosts(2);
   //int fx_kmx  = ilast(2)      + d_fluxghosts(2);
   int fx_nx = fx_imx - fx_imn + 1;
   int fx_ny = fx_jmx - fx_jmn + 1;
   //int fx_nz   = fx_kmx - fx_kmn + 1;
   int fx_nxny = fx_nx * fx_ny;
   //int fx_nel  = fx_nx*fx_ny*fx_nz;

   int fy_imn = ifirst(0) - d_fluxghosts(0);
   int fy_imx = ilast(0) + d_fluxghosts(0);
   int fy_jmn = ifirst(1) - d_fluxghosts(1);
   int fy_jmx = ilast(1) + 1 + d_fluxghosts(1);
   int fy_kmn = ifirst(2) - d_fluxghosts(2);
   //int fy_kmx  = ilast(2)      + d_fluxghosts(2);
   int fy_nx = fy_imx - fy_imn + 1;
   int fy_ny = fy_jmx - fy_jmn + 1;
   //int fy_nz   = fy_kmx - fy_kmn + 1;
   int fy_nxny = fy_nx * fy_ny;
   //int fy_nel  = fy_nx*fy_ny*fy_nz;

   int fz_imn = ifirst(0) - d_fluxghosts(0);
   int fz_imx = ilast(0) + d_fluxghosts(0);
   int fz_jmn = ifirst(1) - d_fluxghosts(1);
   int fz_jmx = ilast(1) + d_fluxghosts(1);
   int fz_kmn = ifirst(2) - d_fluxghosts(2);
   //int fz_kmx  = ilast(2)   + 1 + d_fluxghosts(2);
   int fz_nx = fz_imx - fz_imn + 1;
   int fz_ny = fz_jmx - fz_jmn + 1;
   //int fz_nz   = fz_kmx - fz_kmn + 1;
   int fz_nxny = fz_nx * fz_ny;
   //int fz_nel  = fz_nx*fz_ny*fz_nz;

   int imin = ifirst(0) - d_nghosts(0);
   int imax = ilast(0) + d_nghosts(0);
   int jmin = ifirst(1) - d_nghosts(1);
   int jmax = ilast(1) + d_nghosts(1);
   int kmin = ifirst(2) - d_nghosts(2);
   //int kmax = ilast(2)  + d_nghosts(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   // int nz   = kmax - kmin + 1;
   int nxny = nx * ny;
   // int nel  = nx*ny*nz;

   int nd_imin = ifirst(0) - d_nodeghosts(0);
   int nd_imax = ilast(0) + 1 + d_nodeghosts(0);
   int nd_jmin = ifirst(1) - d_nodeghosts(1);
   int nd_jmax = ilast(1) + 1 + d_nodeghosts(1);
   int nd_kmin = ifirst(2) - d_nodeghosts(2);
   //int nd_kmax = ilast(2)  + 1 + d_nodeghosts(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   // int nd_nz   = nd_kmax - nd_kmin + 1;
   int nd_nxny = nd_nx * nd_ny;
   // int nd_nel  = nd_nx*nd_ny*nd_nz;

   //
   // get the pointers
   //
   double* psi = uval->getPointer();

   double* cvol = vol->getPointer();

   double* fx = flux->getPointer(0);
   double* fy = flux->getPointer(1);
   double* fz = flux->getPointer(2);

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   real8 u = -d_advection_velocity[0];
   real8 v = -d_advection_velocity[1];
   real8 w = -d_advection_velocity[2];

   //
   // compute the source due to the upwind method
   //
   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0) + 1; ++i) {

            int ifx = POLY3(i, j, k, fx_imn, fx_jmn, fx_kmn, fx_nx, fx_nxny);

            // --------- get the neighbors
            int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
            int ib = ind - 1;

            // ---------- use a righthand rule in evaluating face data, n1 -- n4 go around the bottom plane
            //   of the element, n5-n8 go around the top plane of the element.  note, I
            //   have set up the normal so that it always points in the positive, i, j, and k directions
            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n4 = n1 + nd_nx;
            int n5 = n1 + nd_nxny;
            int n8 = n1 + nd_nxny + nd_nx;

            fx[ifx] = UpwindFlux(x[n1], x[n4], x[n8], x[n5],  // 1 - 4 - 8 - 5
                  y[n1], y[n4], y[n8], y[n5],
                  z[n1], z[n4], z[n8], z[n5],
                  u, v, w,
                  psi[ib], psi[ind]);
         }
      }
   }

   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1) + 1; ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int ify = POLY3(i, j, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny);

            // --------- get the neighbors
            int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
            int jb = ind - nx;

            // ---------- use a righthand rule in evaluating face data, n1 -- n4 go around the bottom plane
            //   of the element, n5-n8 go around the top plane of the element.  note, I
            //   have set up the normal so that it always points in the positive, i, j, and k directions
            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;

            fy[ify] = UpwindFlux(x[n1], x[n5], x[n6], x[n2],  // 1 - 5 - 6 - 2
                  y[n1], y[n5], y[n6], y[n2],
                  z[n1], z[n5], z[n6], z[n2],
                  u, v, w,
                  psi[jb], psi[ind]);
         }
      }
   }

   for (int k = ifirst(2); k <= ilast(2) + 1; ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int ifz = POLY3(i, j, k, fz_imn, fz_jmn, fz_kmn, fz_nx, fz_nxny);

            // --------- get the neighbors
            int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
            int kb = ind - nxny;

            // ---------- use a righthand rule in evaluating face data, n1 -- n4 go around the bottom plane
            //   of the element, n5-n8 go around the top plane of the element.  note, I
            //   have set up the normal so that it always points in the positive, i, j, and k directions
            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            fz[ifz] = UpwindFlux(x[n1], x[n2], x[n3], x[n4],  // 1 - 2 - 3 - 4
                  y[n1], y[n2], y[n3], y[n4],
                  z[n1], z[n2], z[n3], z[n4],
                  u, v, w,
                  psi[kb], psi[ind]);
         }
      }
   }

   //
   // compute the source due to the upwind method
   //
   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

            int ib = POLY3(i, j, k, fx_imn, fx_jmn, fx_kmn, fx_nx, fx_nxny);
            int ie = POLY3(i + 1, j, k, fx_imn, fx_jmn, fx_kmn, fx_nx, fx_nxny);

            int jb = POLY3(i, j, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny);
            int je = POLY3(i, j + 1, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny);

            int kb = POLY3(i, j, k, fz_imn, fz_jmn, fz_kmn, fz_nx, fz_nxny);
            int ke = POLY3(i, j, k + 1, fz_imn, fz_jmn, fz_kmn, fz_nx, fz_nxny);

            //   have set up the normal so that it always points in the positive, i, j, and k directions
            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;
            int n7 = n1 + nd_nxny + 1 + nd_nx;
            int n8 = n1 + nd_nxny + nd_nx;

            real8 vol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                  x[n5], x[n6], x[n7], x[n8],

                  y[n1], y[n2], y[n3], y[n4],
                  y[n5], y[n6], y[n7], y[n8],

                  z[n1], z[n2], z[n3], z[n4],
                  z[n5], z[n6], z[n7], z[n8]);

            cvol[ind] = vol;

            psi[ind] -= dt * ((fx[ie] - fx[ib])
                              + (fy[je] - fy[jb])
                              + (fz[ke] - fz[kb])) / vol;

         }
      }
   }

}

/*
 *************************************************************************
 *
 * Update solution variables by performing a conservative
 * difference with the fluxes calculated in computeFluxesOnPatch().
 *
 *************************************************************************
 */

void MblkLinAdv::conservativeDifferenceOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt,
   bool at_syncronization)
{
   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(dt);
   NULL_USE(at_syncronization);
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

void MblkLinAdv::setPhysicalBoundaryConditions(
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

   if (d_dim == tbox::Dimension(2)) {

      /*
       * Set boundary conditions for cells corresponding to patch edges.
       */
      SkeletonBoundaryUtilities2::
      fillEdgeBoundaryData("uval", uval,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_edge_uval);

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      SkeletonBoundaryUtilities2::
      fillNodeBoundaryData("uval", uval,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_edge_uval);

   } // d_dim == tbox::Dimension(2))
   else if (d_dim == tbox::Dimension(3)) {

      /*
       *  Set boundary conditions for cells corresponding to patch faces.
       */

      SkeletonBoundaryUtilities3::
      fillFaceBoundaryData("uval", uval,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_face_conds,
         d_bdry_face_uval);

      /*
       *  Set boundary conditions for cells corresponding to patch edges.
       */

      SkeletonBoundaryUtilities3::
      fillEdgeBoundaryData("uval", uval,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_face_uval);

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      SkeletonBoundaryUtilities3::
      fillNodeBoundaryData("uval", uval,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_face_uval);
   } // d_dim == tbox::Dimension(3))

}

/*
 *************************************************************************
 *
 * Refine operations
 *
 *************************************************************************
 */

void MblkLinAdv::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine_box);
   NULL_USE(ratio);

   int xyz_id = hier::VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   int fln = fine.getPatchLevelNumber();
   int cln = coarse.getPatchLevelNumber();
   if (fln < 0) {
      fln = cln + 1;
      if (!fine.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessRefine()"
                                  << "\nfine xyz data not allocated" << std::endl);
      }
   }
   if (cln < 0) {
      cln = fln - 1;
      if (!coarse.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessRefine()"
                                  << "\ncoarse xyz data not allocated" << std::endl);
      }
   }
   hier::BlockId::block_t block_number =
      coarse.getBox().getBlockId().getBlockValue();
   setMappedGridOnPatch(coarse, cln, block_number);
   setMappedGridOnPatch(fine, fln, block_number);

}

void MblkLinAdv::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{

   std::shared_ptr<pdat::CellData<double> > cuval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > fuval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_vol, getDataContext())));

   TBOX_ASSERT(cuval);
   TBOX_ASSERT(fuval);
   TBOX_ASSERT(cvol);
   TBOX_ASSERT(fvol);
   TBOX_ASSERT(cuval->getDepth() == fuval->getDepth());

   //
   // get the needed geometry and box information
   //
   const hier::Box cgbox(cuval->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fuval->getGhostBox().lower();
   const hier::Index fihi = fuval->getGhostBox().upper();

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   int flev = fine.getPatchLevelNumber();
   int clev = coarse.getPatchLevelNumber();

   tbox::plog << "--------------------- start postprocessRefineData ("
              << flev << "," << clev << ")" << std::endl;
   tbox::plog << "flevel     = " << flev << std::endl;
   tbox::plog << "clevel     = " << clev << std::endl;
   tbox::plog << "fine       = " << fine.getBox() << std::endl;
   tbox::plog << "coarse     = " << coarse.getBox() << std::endl;
   tbox::plog << "fine_box   = " << fine_box << std::endl;
   tbox::plog << "coarse_box = " << coarse_box << std::endl;

   tbox::plog << "filo = " << filo << ", fihi = " << fihi << std::endl;
   tbox::plog << "cilo = " << cilo << ", cihi = " << cihi << std::endl;

   //
   // setup work variables
   //
   const hier::IntVector tmp_ghosts(d_dim, 0);
   const int numInterp = 1;
   pdat::CellData<double> val_vals(coarse_box, numInterp, tmp_ghosts);
   pdat::CellData<double> slope0_vals(coarse_box, numInterp, tmp_ghosts);
   pdat::CellData<double> slope1_vals(coarse_box, numInterp, tmp_ghosts);
   pdat::CellData<double> slope2_vals(coarse_box, numInterp, tmp_ghosts);
   double* val = val_vals.getPointer();
   double* slope0 = slope0_vals.getPointer();
   double* slope1 = slope1_vals.getPointer();
   double* slope2 = slope2_vals.getPointer();

   int imin = ifirstc(0);  // the box of coarse elements being refined
   int imax = ilastc(0);   // the work data is sized to this box
   int jmin = ifirstc(1);
   int jmax = ilastc(1);
   int kmin = ifirstc(2);
   int kmax = ilastc(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nz = kmax - kmin + 1;
   int nxny = nx * ny;
   int nel = nx * ny * nz;

   int cimin = cilo(0);  // the coarse data bounds
   int cimax = cihi(0);
   int cjmin = cilo(1);
   int cjmax = cihi(1);
   int ckmin = cilo(2);
   int ckmax = cihi(2);
   int cnx = cimax - cimin + 1;
   int cny = cjmax - cjmin + 1;
   int cnz = ckmax - ckmin + 1;
   int cnxny = cnx * cny;
   int cnel = cnx * cny * cnz;

   int fimin = filo(0);  // the fine data bounds
   int fimax = fihi(0);
   int fjmin = filo(1);
   int fjmax = fihi(1);
   int fkmin = filo(2);
   int fkmax = fihi(2);
   int fnx = fimax - fimin + 1;
   int fny = fjmax - fjmin + 1;
   int fnz = fkmax - fkmin + 1;
   int fnxny = fnx * fny;
   int fnel = fnx * fny * fnz;

   double rat0 = ratio[0];
   double rat1 = ratio[1];
   double rat2 = ratio[2];
   double fact = 2.0;  // xi varies from -1 to 1

   double* cvolume = cvol->getPointer();
   double* cdata = cuval->getPointer();
   double* fdata = fuval->getPointer();

   //
   // ================================= history variable refinement ====================
   //
   for (int n = 0; n < DEPTH; ++n) {

      for (int l = 0; l < nel; ++l) {       // default slopes are zero
         slope0[l] = 0.0;  // this yields piecewise constant interpolation
         slope1[l] = 0.0;  // and makes a handy initializer
         slope2[l] = 0.0;
      }

      for (int k = ifirstc(2); k <= ilastc(2); ++k) {
         for (int j = ifirstc(1); j <= ilastc(1); ++j) {
            for (int i = ifirstc(0); i <= ilastc(0); ++i) {

               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
               int chind =
                  POLY3(i, j, k, cimin, cjmin, ckmin, cnx, cnxny) + n * cnel;

               int wind = POLY3(i, j, k, cimin, cjmin, ckmin, cnx, cnxny);
               double w_i = cvolume[wind];
               double w_im = cvolume[wind - 1];
               double w_ip = cvolume[wind + 1];
               double w_jm = cvolume[wind - cnx];
               double w_jp = cvolume[wind + cnx];
               double w_km = cvolume[wind - cnxny];
               double w_kp = cvolume[wind + cnxny];

               int im1 = chind - 1;
               int ip1 = chind + 1;
               int jm1 = chind - cnx;
               int jp1 = chind + cnx;
               int km1 = chind - cnxny;
               int kp1 = chind + cnxny;

               double aii = cdata[chind];

               double aip = cdata[ip1];
               double aim = cdata[im1];
               double ajp = cdata[jp1];
               double ajm = cdata[jm1];
               double akp = cdata[kp1];
               double akm = cdata[km1];

               val[ind] = aii;

               TBOX_ASSERT(ind >= 0);   // debug assertions
               TBOX_ASSERT(ind < nel);

               my_slopes(aii, aip, aim, ajp, ajm, akp, akm,
                  w_i, w_ip, w_im, w_jp, w_jm, w_kp, w_km,
                  slope0[ind],
                  slope1[ind],
                  slope2[ind]);

            }
         }
      }

      //
      // compute the interpolated data from the cached slopes, looping
      // over the fine zones
      //
      for (int k = ifirstf(2); k <= ilastf(2); ++k) {
         for (int j = ifirstf(1); j <= ilastf(1); ++j) {
            for (int i = ifirstf(0); i <= ilastf(0); ++i) {

               int fhind =
                  POLY3(i, j, k, fimin, fjmin, fkmin, fnx, fnxny) + n * fnel;

               double ric = (double(i) + 0.5) / rat0 - 0.5;
               double rjc = (double(j) + 0.5) / rat1 - 0.5;
               double rkc = (double(k) + 0.5) / rat2 - 0.5;

               int ic = (int)(ric + (ric >= 0 ? 0.5 : -0.5));
               int jc = (int)(rjc + (rjc >= 0 ? 0.5 : -0.5));
               int kc = (int)(rkc + (rkc >= 0 ? 0.5 : -0.5));

               double ldx = fact * (ric - ic);
               double ldy = fact * (rjc - jc);
               double ldz = fact * (rkc - kc);

               int ind = POLY3(ic, jc, kc, imin, jmin, kmin, nx, nxny);  // work pos

               TBOX_ASSERT(0 <= ind);
               TBOX_ASSERT(ind < nel);

               fdata[fhind] =
                  cdata[ind]
                  + ldx * slope0[ind]
                  + ldy * slope1[ind]
                  + ldz * slope2[ind];

            }
         }
      } // end of i,j,k loops for finding the fine history variables

   } // end of history loop

   tbox::plog << "--------------------- end postprocessRefine" << std::endl;

}

/*
 *************************************************************************
 *
 * Coarsen operations
 *
 *************************************************************************
 */
void MblkLinAdv::preprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   NULL_USE(coarse_box);
   NULL_USE(ratio);

   int xyz_id = hier::VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   int fln = fine.getPatchLevelNumber();
   int cln = coarse.getPatchLevelNumber();
   if (fln < 0) {
      fln = cln + 1;
      if (!fine.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessCoarsen()"
                                  << "\nfine xyz data not allocated" << std::endl);
      }
   }
   if (cln < 0) {
      cln = fln - 1;
      if (!coarse.checkAllocated(xyz_id)) {
         TBOX_ERROR(d_object_name << ":preprocessCoarsen()"
                                  << "\ncoarse xyz data not allocated" << std::endl);
      }
   }
   hier::BlockId::block_t block_number =
      coarse.getBox().getBlockId().getBlockValue();
   setMappedGridOnPatch(coarse, cln, block_number);
   setMappedGridOnPatch(fine, fln, block_number);

}

//
// the coarsening function
//
void MblkLinAdv::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{

   std::shared_ptr<pdat::CellData<double> > cuval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > fuval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_uval, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_vol, getDataContext())));

   TBOX_ASSERT(cuval);
   TBOX_ASSERT(cvol);
   TBOX_ASSERT(fuval);
   TBOX_ASSERT(fvol);
   TBOX_ASSERT(cuval->getDepth() == fuval->getDepth());

   //
   // box and geometry information
   //
   const hier::Index filo = fuval->getGhostBox().lower();
   const hier::Index fihi = fuval->getGhostBox().upper();
   const hier::Index cilo = cuval->getGhostBox().lower();
   const hier::Index cihi = cuval->getGhostBox().upper();

   const hier::Box fine_box = hier::Box::refine(coarse_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();  // coarse basis
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();    // fine data
   const hier::Index ilastf = fine_box.upper();

   int flev = fine.getPatchLevelNumber();
   int clev = coarse.getPatchLevelNumber();

   tbox::plog << "--------------------- start postprocessCoarsenData ("
              << flev << "," << clev << ")" << std::endl;
   tbox::plog << "flevel     = " << flev << std::endl;
   tbox::plog << "clevel     = " << clev << std::endl;
   tbox::plog << "fine       = " << fine.getBox() << std::endl;
   tbox::plog << "coarse     = " << coarse.getBox() << std::endl;
   tbox::plog << "fine_box   = " << fine_box << std::endl;
   tbox::plog << "coarse_box = " << coarse_box << std::endl;

   tbox::plog << "filo = " << filo << ", fihi = " << fihi << std::endl;
   tbox::plog << "cilo = " << cilo << ", cihi = " << cihi << std::endl;

   //
   // work variables
   //

   int cimin = cilo(0);
   int cimax = cihi(0);
   int cjmin = cilo(1);
   int cjmax = cihi(1);
   int ckmin = cilo(2);
   int ckmax = cihi(2);
   int cnx = cimax - cimin + 1;
   int cny = cjmax - cjmin + 1;
   int cnz = ckmax - ckmin + 1;
   int cnxny = cnx * cny;
   int cnel = cnx * cny * cnz;

   int fimin = filo(0);
   int fimax = fihi(0);
   int fjmin = filo(1);
   int fjmax = fihi(1);
   int fkmin = filo(2);
   int fkmax = fihi(2);
   int fnx = fimax - fimin + 1;
   int fny = fjmax - fjmin + 1;
   int fnz = fkmax - fkmin + 1;
   int fnxny = fnx * fny;
   int fnel = fnx * fny * fnz;

   double rat0 = (double)(ratio(0));
   double rat1 = (double)(ratio(1));
   double rat2 = (double)(ratio(2));

   double* cvolume = cvol->getPointer();
   double* fvolume = fvol->getPointer();
   double* cdata = cuval->getPointer();
   double* fdata = fuval->getPointer();

   //
   // average the data
   //
   for (int n = 0; n < DEPTH; ++n) {

      //
      // zero out the underlying coarse data to serve as a counter
      //
      for (int k = ifirstc(2); k <= ilastc(2); ++k) {    // loop over the coarse zones
         for (int j = ifirstc(1); j <= ilastc(1); ++j) {
            for (int i = ifirstc(0); i <= ilastc(0); ++i) {

               int chind =
                  POLY3(i, j, k, cimin, cjmin, ckmin, cnx, cnxny) + n * cnel;

               cdata[chind] = 0.0;
            }
         }
      }

      //
      // compute the interpolated data from the cached slopes
      //

      for (int k = ifirstf(2); k <= ilastf(2); ++k) {    // loop over the coarse zones
         for (int j = ifirstf(1); j <= ilastf(1); ++j) {
            for (int i = ifirstf(0); i <= ilastf(0); ++i) {

               int vol_ind = POLY3(i, j, k, fimin, fjmin, fkmin, fnx, fnxny);
               int fhind = vol_ind + n * fnel;

               real8 ric = (double(i) + 0.5) / rat0 - 0.5;
               real8 rjc = (double(j) + 0.5) / rat1 - 0.5;
               real8 rkc = (double(k) + 0.5) / rat2 - 0.5;

               int ic = (int)(ric + (ric >= 0 ? 0.5 : -0.5));        // a round operation
               int jc = (int)(rjc + (rjc >= 0 ? 0.5 : -0.5));        // shift up and truncate if ic > 0
               int kc = (int)(rkc + (rkc >= 0 ? 0.5 : -0.5));        // shift down and truncate if ic < 0

               int chind =
                  POLY3(ic, jc, kc, cimin, cjmin, ckmin, cnx, cnxny) + n * cnel;  // post + history offset

               TBOX_ASSERT(cimin <= ic && ic <= cimax);
               TBOX_ASSERT(cjmin <= jc && jc <= cjmax);
               TBOX_ASSERT(ckmin <= kc && kc <= ckmax);

               real8 fmass = fvolume[vol_ind];

               cdata[chind] += fmass * fdata[fhind];  // sum of fine extensives now
            }
         }
      }

      //
      // normalize the completed sum by converting back from extensive to intensive
      //
      for (int k = ifirstc(2); k <= ilastc(2); ++k) {    // loop over the coarse zones
         for (int j = ifirstc(1); j <= ilastc(1); ++j) {
            for (int i = ifirstc(0); i <= ilastc(0); ++i) {

               int vol_ind = POLY3(i, j, k, cimin, cjmin, ckmin, cnx, cnxny);
               int chind = vol_ind + n * cnel;

               real8 cmass = cvolume[vol_ind];

               cdata[chind] /= cmass + 1.e-80;
            }
         }
      }

   }
   tbox::plog << "--------------------- end postprocessCoarsen" << std::endl;
}

// -------------------------------------------------------------------

/*
 *************************************************************************
 *
 * Tag cells for refinement using gradient detector.  Tagging criteria
 * defined in input.
 *
 *************************************************************************
 */

//
// tag cells for refinement
//
void MblkLinAdv::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(regrid_time);
   NULL_USE(uses_richardson_extrapolation_too);

   NULL_USE(initial_error);

   hier::BlockId::block_t block_number =
      patch.getBox().getBlockId().getBlockValue();
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch, level_number, block_number);

   //
   // get geometry data
   //
   const int error_level_number = patch.getPatchLevelNumber();

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));
   TBOX_ASSERT(xyz);
   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();
   int level = patch.getPatchLevelNumber();

   tbox::plog << "--------------------- start tagGradientCells (" << level
              << ")" << std::endl;
   tbox::plog << "level  = " << level << std::endl;
   tbox::plog << "box    = " << patch.getBox() << std::endl;

   std::shared_ptr<pdat::CellData<int> > tags(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(tag_indx)));
   std::shared_ptr<pdat::CellData<double> > var(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_uval, getDataContext())));
   TBOX_ASSERT(tags);
   TBOX_ASSERT(var);

   //
   // Create a set of temporary tags and set to untagged value.
   //
   std::shared_ptr<pdat::CellData<int> > temp_tags(
      new pdat::CellData<int>(pbox, 1, d_nghosts));
   temp_tags->fillAll(FALSE);

   hier::IntVector tag_ghost = tags->getGhostCellWidth();

   hier::IntVector nghost_cells = xyz->getGhostCellWidth();
   int nd_imin = ifirst(0) - nghost_cells(0);
   int nd_imax = ilast(0) + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1) - nghost_cells(1);
   int nd_jmax = ilast(1) + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2) - nghost_cells(2);
   //int nd_kmax = ilast(2)  + 1 + nghost_cells(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   //int nd_nz   = nd_kmax - nd_kmin + 1;
   int nd_nxny = nd_nx * nd_ny;
   //int nd_nel  = nd_nx*nd_ny*nd_nz;

   hier::IntVector v_ghost = var->getGhostCellWidth();     // has ghost zones
   int imin = ifirst(0) - v_ghost(0); // the polynomial for the field
   int imax = ilast(0) + v_ghost(0);
   int jmin = ifirst(1) - v_ghost(1);
   int jmax = ilast(1) + v_ghost(1);
   int kmin = ifirst(2) - v_ghost(2);
   //int kmax = ilast(2)  + v_ghost(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   // int nz   = kmax - kmin + 1;
   int nxny = nx * ny;
   //int nel  = nx*ny*nz;

   hier::IntVector temp_tag_ghost = temp_tags->getGhostCellWidth();
   int imn = ifirst(0) - temp_tag_ghost(0);  // the polynomial for temp_tags
   int imx = ilast(0) + temp_tag_ghost(0);
   int jmn = ifirst(1) - temp_tag_ghost(1);
   int jmx = ilast(1) + temp_tag_ghost(1);
   int kmn = ifirst(2) - temp_tag_ghost(2);
   //int kmx = ilast(2)  + temp_tag_ghost(2);
   int tnx = imx - imn + 1;
   int tny = jmx - jmn + 1;
   //int tnz   = kmx - kmn + 1;
   int tnxny = tnx * tny;
   //int tnel  = tnx*tny*tnz;

   double* lvar = var->getPointer();
   int* ltags = temp_tags->getPointer();
   double dv_x[3];
   double dv_xi[3];
   double xfact = 0.5;

   //
   // Possible tagging criteria includes
   //    DENSITY_DEVIATION, DENSITY_GRADIENT, DENSITY_SHOCK
   //    PRESSURE_DEVIATION, PRESSURE_GRADIENT, PRESSURE_SHOCK
   // The criteria is specified over a time interval.
   //
   // Loop over criteria provided and check to make sure we are in the
   // specified time interval.  If so, apply appropriate tagging for
   // the level.
   //
   for (int ncrit = 0;
        ncrit < static_cast<int>(d_refinement_criteria.size()); ++ncrit) {

      std::string ref = d_refinement_criteria[ncrit];
      int size = 0;
      double tol = 0.;

      if (ref == "UVAL_GRADIENT") {
         size = static_cast<int>(d_grad_tol.size());  // max depth of gradient tolerance
         tol = ((error_level_number < size)    // find the tolerance
                ? d_grad_tol[error_level_number]
                : d_grad_tol[size - 1]);

         for (int k = ifirst(2); k <= ilast(2); ++k) {
            for (int j = ifirst(1); j <= ilast(1); ++j) {
               for (int i = ifirst(0); i <= ilast(0); ++i) {

                  int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
                  int tind = POLY3(i, j, k, imn, jmn, kmn, tnx, tnxny);

                  int ib = ind - 1;
                  int ie = ind + 1;
                  int jb = ind - nx;
                  int je = ind + nx;
                  int kb = ind - nxny;
                  int ke = ind + nxny;

                  //
                  // vector is now a max gradient in xi, eta, zeta
                  //
                  dv_xi[0] = xfact
                     * MAX(fabs(lvar[ind] - lvar[ib]),
                        fabs(lvar[ind] - lvar[ie]));
                  dv_xi[1] = xfact
                     * MAX(fabs(lvar[ind] - lvar[jb]),
                        fabs(lvar[ind] - lvar[je]));
                  dv_xi[2] = xfact
                     * MAX(fabs(lvar[ind] - lvar[kb]),
                        fabs(lvar[ind] - lvar[ke]));

                  //
                  // convert vector to max gradine in x, y, z
                  //
                  int n1 = POLY3(i,
                        j,
                        k,
                        nd_imin,
                        nd_jmin,
                        nd_kmin,
                        nd_nx,
                        nd_nxny);                                               // -1, -1, -1
                  int n2 = n1 + 1;                                              //  1, -1, -1
                  int n3 = n1 + 1 + nd_nx;                                      //  1,  1, -1
                  int n4 = n1 + nd_nx;                                          // -1,  1, -1

                  int n5 = n1 + nd_nxny;                                        // -1, -1,  1
                  int n6 = n1 + nd_nxny + 1;                                 //  1, -1,  1
                  int n7 = n1 + nd_nxny + 1 + nd_nx;                         //  1,  1,  1
                  int n8 = n1 + nd_nxny + nd_nx;                                     // -1,  1,  1

                  // ------------------------------------------------ x

                  real8 x1 = 0.25 * (x[n1] + x[n4] + x[n5] + x[n8]);  // xi
                  real8 x2 = 0.25 * (x[n2] + x[n3] + x[n6] + x[n7]);

                  real8 x3 = 0.25 * (x[n1] + x[n2] + x[n5] + x[n6]);  // eta
                  real8 x4 = 0.25 * (x[n3] + x[n4] + x[n7] + x[n8]);

                  real8 x5 = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4]);  // zeta
                  real8 x6 = 0.25 * (x[n5] + x[n6] + x[n7] + x[n8]);

                  // ------------------------------------------------ y

                  real8 y1 = 0.25 * (y[n1] + y[n4] + y[n5] + y[n8]);
                  real8 y2 = 0.25 * (y[n2] + y[n3] + y[n6] + y[n7]);

                  real8 y3 = 0.25 * (y[n1] + y[n2] + y[n5] + y[n6]);
                  real8 y4 = 0.25 * (y[n3] + y[n4] + y[n7] + y[n8]);

                  real8 y5 = 0.25 * (y[n1] + y[n2] + y[n3] + y[n4]);
                  real8 y6 = 0.25 * (y[n5] + y[n6] + y[n7] + y[n8]);

                  // ------------------------------------------------ z

                  real8 z1 = 0.25 * (z[n1] + z[n4] + z[n5] + z[n8]);
                  real8 z2 = 0.25 * (z[n2] + z[n3] + z[n6] + z[n7]);

                  real8 z3 = 0.25 * (z[n1] + z[n2] + z[n5] + z[n6]);
                  real8 z4 = 0.25 * (z[n3] + z[n4] + z[n7] + z[n8]);

                  real8 z5 = 0.25 * (z[n1] + z[n2] + z[n3] + z[n4]);
                  real8 z6 = 0.25 * (z[n5] + z[n6] + z[n7] + z[n8]);

                  //
                  // the components of the matrices that we want to invert
                  //
                  real8 dx_xi = 0.5 * (x2 - x1);
                  real8 dy_xi = 0.5 * (y2 - y1);
                  real8 dz_xi = 0.5 * (z2 - z1);

                  real8 dx_eta = 0.5 * (x4 - x3);
                  real8 dy_eta = 0.5 * (y4 - y3);
                  real8 dz_eta = 0.5 * (z4 - z3);

                  real8 dx_zeta = 0.5 * (x6 - x5);
                  real8 dy_zeta = 0.5 * (y6 - y5);
                  real8 dz_zeta = 0.5 * (z6 - z5);

                  //
                  // invert dx/dxi as in dx/dxi d/dx = d/dxi, note this
                  // is the transpose of the above matrix M, also via
                  // Kramer's rule
                  //
                  real8 detMt = (dx_xi * dy_eta * dz_zeta
                                 + dx_eta * dy_zeta * dz_xi
                                 + dx_zeta * dy_xi * dz_eta
                                 - dx_zeta * dy_eta * dz_xi
                                 - dx_eta * dy_xi * dz_zeta
                                 - dx_xi * dy_zeta * dz_eta);

                  real8 detC11 = dy_eta * dz_zeta - dz_eta * dy_zeta;
                  real8 detC21 = dx_eta * dz_zeta - dz_eta * dx_zeta;
                  real8 detC31 = dx_eta * dy_zeta - dy_eta * dx_zeta;

                  real8 detC12 = dy_xi * dz_zeta - dz_xi * dy_zeta;
                  real8 detC22 = dx_xi * dz_zeta - dz_xi * dx_zeta;
                  real8 detC32 = dx_xi * dy_zeta - dy_xi * dx_zeta;

                  real8 detC13 = dy_xi * dz_eta - dz_xi * dy_eta;
                  real8 detC23 = dx_xi * dz_eta - dz_xi * dx_eta;
                  real8 detC33 = dx_xi * dy_eta - dy_xi + dx_eta;

                  // -------------------

                  real8 b11 = detC11 / detMt;
                  real8 b21 = detC21 / detMt;
                  real8 b31 = detC31 / detMt;

                  real8 b12 = detC12 / detMt;
                  real8 b22 = detC22 / detMt;
                  real8 b32 = detC32 / detMt;

                  real8 b13 = detC13 / detMt;
                  real8 b23 = detC23 / detMt;
                  real8 b33 = detC33 / detMt;

                  //
                  // determine the maximum gradient in x, y and z (nice orthonormal basis)
                  //
                  dv_x[0] = b11 * dv_xi[0] + b12 * dv_xi[1] + b13 * dv_xi[2];
                  dv_x[1] = b21 * dv_xi[0] + b22 * dv_xi[1] + b23 * dv_xi[2];
                  dv_x[2] = b31 * dv_xi[0] + b32 * dv_xi[1] + b33 * dv_xi[2];

                  double vmax = MAX(dv_x[0], MAX(dv_x[1], dv_x[2]));

                  if (vmax > tol)
                     ltags[tind] = TRUE;
               }
            }
         }

      } // criteria = UVAL_GRADIENT

      if (ref == "USER_DEFINED") {

         /*
          * For user-defined, access refine box data from the MblkGeometry
          * class.
          */
         hier::BlockId::block_t block_number =
            patch.getBox().getBlockId().getBlockValue();
         int level_number = patch.getPatchLevelNumber();
         hier::BoxContainer refine_boxes;
         if (d_mblk_geometry->getRefineBoxes(refine_boxes,
                block_number,
                level_number)) {
            for (hier::BoxContainer::iterator b = refine_boxes.begin();
                 b != refine_boxes.end(); ++b) {
               hier::Box intersect = pbox * (*b);
               if (!intersect.empty()) {
                  temp_tags->fill(TRUE, intersect);
               }
            }
         }

      } // criteria = USER_DEFINED

   }  // loop over criteria

   //
   // Update tags
   //
   pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
   for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
        ic != icend; ++ic) {
      (*tags)(*ic, 0) = (*temp_tags)(*ic, 0);
   }

   tbox::plog << "--------------------- end tagGradientCells" << std::endl;

}

/*
 *************************************************************************
 *                                                                       *
 * Fill the singularity conditions for the multi-block case
 *                                                                       *
 *************************************************************************
 */
void MblkLinAdv::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   std::shared_ptr<const hier::Connector> dst_to_encon,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{

   NULL_USE(patch);
   NULL_USE(encon_level);
   NULL_USE(dst_to_encon);
   NULL_USE(fill_box);
   NULL_USE(boundary_box);
   NULL_USE(grid_geometry);

}

/*
 *************************************************************************
 *                                                                       *
 * Private method to build XYZ coordinates on a patch                    *
 *                                                                       *
 *************************************************************************
 */
void MblkLinAdv::setMappedGridOnPatch(
   const hier::Patch& patch,
   const int level_number,
   const hier::BlockId::block_t block_number)
{
   TBOX_ASSERT(level_number >= 0);

   // compute level domain
   const std::shared_ptr<hier::PatchGeometry> patch_geom(
      patch.getPatchGeometry());
   hier::IntVector ratio = patch_geom->getRatio();
   hier::BoxContainer domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, ratio,
      hier::BlockId(block_number));
   int num_domain_boxes = domain_boxes.size();

   if (num_domain_boxes > 1) {
      TBOX_ERROR("Sorry, cannot handle non-rectangular domains..." << std::endl);
   }

   int xyz_id = hier::VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   d_mblk_geometry->buildGridOnPatch(patch,
      domain_boxes.front(),
      xyz_id,
      level_number,
      block_number);
}

/*
 *************************************************************************
 *                                                                       *
 * Register VisIt data writer to write data to plot files that may       *
 * be postprocessed by the VisIt tool.                                   *
 *                                                                       *
 *************************************************************************
 */

#ifdef HAVE_HDF5
void MblkLinAdv::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);
   d_visit_writer = viz_writer;
}
#endif

/*
 *************************************************************************
 *                                                                       *
 * Write MblkLinAdv object state to specified stream.                        *
 *                                                                       *
 *************************************************************************
 */

void MblkLinAdv::printClassData(
   std::ostream& os) const
{
   int j, k;

   os << "\nMblkLinAdv::printClassData..." << std::endl;
   os << "MblkLinAdv: this = " << (MblkLinAdv *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_grid_geometry = " << std::endl;
//   for (j=0; j < d_grid_geometry.getSize(); ++j) {
//      os << (*((std::shared_ptr<geom::GridGeometry >)(d_grid_geometry[j]))) << std::endl;
//   }

   os << "Parameters for numerical method ..." << std::endl;
   os << "   d_advection_velocity = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_advection_velocity[j] << " ";
   os << std::endl;
   os << "   d_godunov_order = " << d_godunov_order << std::endl;
   os << "   d_corner_transport = " << d_corner_transport << std::endl;
   os << "   d_nghosts = " << d_nghosts << std::endl;
   os << "   d_fluxghosts = " << d_fluxghosts << std::endl;

   os << "Problem description and initial data..." << std::endl;
   os << "   d_data_problem = " << d_data_problem << std::endl;
   os << "   d_data_problem_int = " << d_data_problem << std::endl;

   os << "       d_radius = " << d_radius << std::endl;
   os << "       d_center = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_center[j] << " ";
   os << std::endl;
   os << "       d_uval_inside = " << d_uval_inside << std::endl;
   os << "       d_uval_outside = " << d_uval_outside << std::endl;

   os << "       d_number_of_intervals = " << d_number_of_intervals << std::endl;
   os << "       d_front_position = ";
   for (k = 0; k < d_number_of_intervals - 1; ++k) {
      os << d_front_position[k] << "  ";
   }
   os << std::endl;
   os << "       d_interval_uval = " << std::endl;
   for (k = 0; k < d_number_of_intervals; ++k) {
      os << "            " << d_interval_uval[k] << std::endl;
   }
   os << "   Boundary condition data " << std::endl;

   if (d_dim == tbox::Dimension(2)) {
      for (j = 0; j < static_cast<int>(d_scalar_bdry_edge_conds.size()); ++j) {
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << std::endl;
         if (d_scalar_bdry_edge_conds[j] == BdryCond::DIRICHLET) {
            os << "         d_bdry_edge_uval[" << j << "] = "
               << d_bdry_edge_uval[j] << std::endl;
         }
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_scalar_bdry_node_conds.size()); ++j) {
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << std::endl;
         os << "       d_node_bdry_edge[" << j << "] = "
            << d_node_bdry_edge[j] << std::endl;
      }
   } else if (d_dim == tbox::Dimension(3)) {
      for (j = 0; j < static_cast<int>(d_scalar_bdry_face_conds.size()); ++j) {
         os << "       d_scalar_bdry_face_conds[" << j << "] = "
            << d_scalar_bdry_face_conds[j] << std::endl;
         if (d_scalar_bdry_face_conds[j] == BdryCond::DIRICHLET) {
            os << "         d_bdry_face_uval[" << j << "] = "
               << d_bdry_face_uval[j] << std::endl;
         }
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_scalar_bdry_edge_conds.size()); ++j) {
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << std::endl;
         os << "       d_edge_bdry_face[" << j << "] = "
            << d_edge_bdry_face[j] << std::endl;
      }
      os << std::endl;
      for (j = 0; j < static_cast<int>(d_scalar_bdry_node_conds.size()); ++j) {
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << std::endl;
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
 *                                                                       *
 * Read data members from input.  All values set from restart can be     *
 * overridden by values in the input database.                           *
 *                                                                       *
 *************************************************************************
 */
void MblkLinAdv::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   std::shared_ptr<tbox::Database> db(input_db->getDatabase("MblkLinAdv"));

   /*
    * Note: if we are restarting, then we only allow nonuniform
    * workload to be used if nonuniform workload was used originally.
    */
   if (!is_from_restart) {
      d_use_nonuniform_workload =
         db->getBoolWithDefault("use_nonuniform_workload",
            d_use_nonuniform_workload);
   } else {
      if (d_use_nonuniform_workload) {
         d_use_nonuniform_workload =
            db->getBool("use_nonuniform_workload");
      }
   }

   if (db->keyExists("advection_velocity")) {
      db->getDoubleArray("advection_velocity",
         d_advection_velocity, d_dim.getValue());
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `advection_velocity' not found in input.");
   }

   if (db->keyExists("godunov_order")) {
      d_godunov_order = db->getInteger("godunov_order");
      if ((d_godunov_order != 1) &&
          (d_godunov_order != 2) &&
          (d_godunov_order != 4)) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`godunov_order' in input must be 1, 2, or 4." << std::endl);
      }
   } else {
      d_godunov_order = db->getIntegerWithDefault("d_godunov_order",
            d_godunov_order);
   }

   if (db->keyExists("corner_transport")) {
      d_corner_transport = db->getString("corner_transport");
      if ((d_corner_transport != "CORNER_TRANSPORT_1") &&
          (d_corner_transport != "CORNER_TRANSPORT_2")) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`corner_transport' in input must be either string"
                          << " 'CORNER_TRANSPORT_1' or 'CORNER_TRANSPORT_2'." << std::endl);
      }
   } else {
      d_corner_transport = db->getStringWithDefault("corner_transport",
            d_corner_transport);
   }

   if (db->keyExists("Refinement_data")) {
      std::shared_ptr<tbox::Database> refine_db(
         db->getDatabase("Refinement_data"));
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

   } // refine db entry exists

   if (!is_from_restart) {

      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`data_problem' value not found in input."
                          << std::endl);
      }

      if (!db->keyExists("Initial_data")) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "No `Initial_data' database found in input." << std::endl);
      }
      std::shared_ptr<tbox::Database> init_data_db(
         db->getDatabase("Initial_data"));

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
         if (init_data_db->keyExists("uval_inside")) {
            d_uval_inside = init_data_db->getDouble("uval_inside");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`uval_inside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("uval_outside")) {
            d_uval_outside = init_data_db->getDouble("uval_outside");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`uval_outside' input required for "
                                     << "SPHERE problem." << std::endl);
         }

         found_problem_data = true;

      }

      if (!found_problem_data && (
             (d_data_problem == "PIECEWISE_CONSTANT_X") ||
             (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
             (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
             (d_data_problem == "SINE_CONSTANT_X") ||
             (d_data_problem == "SINE_CONSTANT_Y") ||
             (d_data_problem == "SINE_CONSTANT_Z"))) {

         if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
            if (d_dim < tbox::Dimension(2)) {
               TBOX_ERROR(
                  d_object_name << ": `PIECEWISE_CONSTANT_Y' "
                                << "problem invalid in 1 dimension."
                                << std::endl);
            }
         }

         if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
            if (d_dim < tbox::Dimension(3)) {
               TBOX_ERROR(
                  d_object_name << ": `PIECEWISE_CONSTANT_Z' "
                                << "problem invalid in 1 or 2 dimensions." << std::endl);
            }
         }

         std::vector<std::string> init_data_keys = init_data_db->getAllKeys();

         if (init_data_db->keyExists("front_position")) {
            d_front_position = init_data_db->getDoubleVector("front_position");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`front_position' input required for "
                                     << d_data_problem << " problem." << std::endl);
         }

         d_number_of_intervals =
            tbox::MathUtilities<int>::Min(static_cast<int>(d_front_position.size()) + 1,
               static_cast<int>(init_data_keys.size()) - 1);

         d_interval_uval.resize(d_number_of_intervals);

         int i = 0;
         int nkey = 0;
         bool found_interval_data = false;

         while (!found_interval_data
                && (i < d_number_of_intervals)
                && (nkey < static_cast<int>(init_data_keys.size()))) {

            if (!(init_data_keys[nkey] == "front_position")) {

               std::shared_ptr<tbox::Database> interval_db(
                  init_data_db->getDatabase(init_data_keys[nkey]));

               if (interval_db->keyExists("uval")) {
                  d_interval_uval[i] = interval_db->getDouble("uval");
               } else {
                  TBOX_ERROR(d_object_name << ": "
                                           << "`uval' data missing in input for key = "
                                           << init_data_keys[nkey] << std::endl);
               }
               ++i;

               found_interval_data = (i == d_number_of_intervals);

            }

            ++nkey;

         }

         if ((d_data_problem == "SINE_CONSTANT_X") ||
             (d_data_problem == "SINE_CONSTANT_Y") ||
             (d_data_problem == "SINE_CONSTANT_Z")) {
            if (init_data_db->keyExists("amplitude")) {
               d_amplitude = init_data_db->getDouble("amplitude");
            }
            if (init_data_db->keyExists("frequency")) {
               init_data_db->getDoubleArray("frequency", d_frequency, d_dim.getValue());
            } else {
               TBOX_ERROR(
                  d_object_name << ": "
                                << "`frequency' input required for SINE problem." << std::endl);
            }
         }

         if (!found_interval_data) {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Insufficient interval data given in input"
                             << " for PIECEWISE_CONSTANT_*problem."
                             << std::endl);
         }

         found_problem_data = true;
      }

      if (!found_problem_data) {
         TBOX_ERROR(d_object_name << ": "
                                  << "`Initial_data' database found in input."
                                  << " But bad data supplied." << std::endl);
      }

   } // if !is_from_restart read in problem data

   hier::IntVector periodic = d_grid_geometry->getPeriodicShift(
         hier::IntVector(d_dim, 1));
   int num_per_dirs = 0;
   for (int id = 0; id < d_dim.getValue(); ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   /*
    * If there are multiple blocks, periodicity is not currently supported.
    */
   if ((d_grid_geometry->getNumberBlocks() > 1) && (num_per_dirs > 0)) {
      TBOX_ERROR(d_object_name << ": cannot have periodic BCs when there"
                               << "\nare multiple blocks." << std::endl);
   }

   if (db->keyExists("Boundary_data")) {

      std::shared_ptr<tbox::Database> bdry_db(
         db->getDatabase("Boundary_data"));

      if (d_dim == tbox::Dimension(2)) {
         SkeletonBoundaryUtilities2::getFromInput(this,
            bdry_db,
            d_scalar_bdry_edge_conds,
            d_scalar_bdry_node_conds,
            periodic);
      } else if (d_dim == tbox::Dimension(3)) {
         SkeletonBoundaryUtilities3::getFromInput(this,
            bdry_db,
            d_scalar_bdry_face_conds,
            d_scalar_bdry_edge_conds,
            d_scalar_bdry_node_conds,
            periodic);
      }

   } else {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Key data `Boundary_data' not found in input. " << std::endl);
   }

}

/*
 *************************************************************************
 *                                                                       *
 * Routines to put/get data members to/from restart database.            *
 *                                                                       *
 *************************************************************************
 */

void MblkLinAdv::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("MBLKLINADV_VERSION", MBLKLINADV_VERSION);

   restart_db->putDoubleArray("d_advection_velocity",
      d_advection_velocity,
      d_dim.getValue());

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
      restart_db->putDouble("d_uval_inside", d_uval_inside);
      restart_db->putDouble("d_uval_outside", d_uval_outside);
   }

   if ((d_data_problem == "PIECEWISE_CONSTANT_X") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
       (d_data_problem == "SINE_CONSTANT_X") ||
       (d_data_problem == "SINE_CONSTANT_Y") ||
       (d_data_problem == "SINE_CONSTANT_Z")) {
      restart_db->putInteger("d_number_of_intervals", d_number_of_intervals);
      if (d_number_of_intervals > 0) {
         restart_db->putDoubleVector("d_front_position", d_front_position);
         restart_db->putDoubleVector("d_interval_uval", d_interval_uval);
      }
   }

   restart_db->putIntegerVector("d_scalar_bdry_edge_conds",
      d_scalar_bdry_edge_conds);
   restart_db->putIntegerVector("d_scalar_bdry_node_conds",
      d_scalar_bdry_node_conds);

   if (d_dim == tbox::Dimension(2)) {
      restart_db->putDoubleVector("d_bdry_edge_uval", d_bdry_edge_uval);
   } else if (d_dim == tbox::Dimension(3)) {
      restart_db->putIntegerVector("d_scalar_bdry_face_conds",
         d_scalar_bdry_face_conds);
      restart_db->putDoubleVector("d_bdry_face_uval", d_bdry_face_uval);
   }

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
 *                                                                       *
 *    Access class information from restart database.                    *
 *                                                                       *
 *************************************************************************
 */
void MblkLinAdv::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file.");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("MBLKLINADV_VERSION");
   if (ver != MBLKLINADV_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   db->getDoubleArray("d_advection_velocity", d_advection_velocity, d_dim.getValue());

   d_godunov_order = db->getInteger("d_godunov_order");
   d_corner_transport = db->getString("d_corner_transport");

#if 0
   int* tmp_nghosts = d_nghosts;
   db->getIntegerArray("d_nghosts", tmp_nghosts, d_dim);
   if (!(d_nghosts == CELLG)) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Key data `d_nghosts' in restart file != CELLG." << std::endl);
   }
   int* tmp_fluxghosts = d_fluxghosts;
   db->getIntegerArray("d_fluxghosts", tmp_fluxghosts, d_dim);
   if (!(d_fluxghosts == FLUXG)) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Key data `d_fluxghosts' in restart file != FLUXG." << std::endl);
   }
#endif

   d_data_problem = db->getString("d_data_problem");

   if (d_data_problem == "SPHERE") {
      d_data_problem_int = SPHERE;
      d_radius = db->getDouble("d_radius");
      db->getDoubleArray("d_center", d_center, d_dim.getValue());
      d_uval_inside = db->getDouble("d_uval_inside");
      d_uval_outside = db->getDouble("d_uval_outside");
   }

   if ((d_data_problem == "PIECEWISE_CONSTANT_X") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Y") ||
       (d_data_problem == "PIECEWISE_CONSTANT_Z") ||
       (d_data_problem == "SINE_CONSTANT_X") ||
       (d_data_problem == "SINE_CONSTANT_Y") ||
       (d_data_problem == "SINE_CONSTANT_Z")) {
      d_number_of_intervals = db->getInteger("d_number_of_intervals");
      if (d_number_of_intervals > 0) {
         d_front_position = db->getDoubleVector("d_front_position");
         d_interval_uval = db->getDoubleVector("d_interval_uval");
      }
   }

   d_scalar_bdry_edge_conds = db->getIntegerVector("d_scalar_bdry_edge_conds");
   d_scalar_bdry_node_conds = db->getIntegerVector("d_scalar_bdry_node_conds");

   if (d_dim == tbox::Dimension(2)) {
      d_bdry_edge_uval = db->getDoubleVector("d_bdry_edge_uval");
   } else if (d_dim == tbox::Dimension(3)) {
      d_scalar_bdry_face_conds =
         db->getIntegerVector("d_scalar_bdry_face_conds");

      d_bdry_face_uval = db->getDoubleVector("d_bdry_face_uval");
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

/*
 *************************************************************************
 *                                                                       *
 * Routines to read boundary data from input database.                   *
 *                                                                       *
 *************************************************************************
 */

void MblkLinAdv::readDirichletBoundaryDataEntry(
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
         d_bdry_edge_uval);
   } else if (d_dim == tbox::Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_uval);
   }
}

void MblkLinAdv::readStateDataEntry(
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

void MblkLinAdv::readNeumannBoundaryDataEntry(
   const std::shared_ptr<tbox::Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   NULL_USE(db);
   NULL_USE(db_name);
   NULL_USE(bdry_location_index);
}

/*
 *************************************************************************
 *                                                                       *
 * Routine to check boundary data when debugging.                        *
 *                                                                       *
 *************************************************************************
 */

void MblkLinAdv::checkBoundaryData(
   int btype,
   const hier::Patch& patch,
   const hier::IntVector& ghost_width_to_check,
   const std::vector<int>& scalar_bconds) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim == tbox::Dimension(2)) {
      TBOX_ASSERT(btype == Bdry::EDGE2D ||
         btype == Bdry::NODE2D);
   } else if (d_dim == tbox::Dimension(3)) {
      TBOX_ASSERT(btype == Bdry::FACE3D ||
         btype == Bdry::EDGE3D ||
         btype == Bdry::NODE3D);
   }
#endif

   const std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());
   const std::vector<hier::BoundaryBox>& bdry_boxes =
      pgeom->getCodimensionBoundaries(btype);

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   for (int i = 0; i < static_cast<int>(bdry_boxes.size()); ++i) {
      hier::BoundaryBox bbox = bdry_boxes[i];
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
      int bloc = bbox.getLocationIndex();

      int bscalarcase = 0, refbdryloc = 0;
      if (d_dim == tbox::Dimension(2)) {
         if (btype == Bdry::EDGE2D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) == NUM_2D_EDGES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = bloc;
         } else { // btype == Bdry::NODE2D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) == NUM_2D_NODES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_node_bdry_edge[bloc];
         }
      } else if (d_dim == tbox::Dimension(3)) {
         if (btype == Bdry::FACE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) == NUM_3D_FACES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = bloc;
         } else if (btype == Bdry::EDGE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) == NUM_3D_EDGES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_edge_bdry_face[bloc];
         } else { // btype == Bdry::NODE3D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) == NUM_3D_NODES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_node_bdry_face[bloc];
         }
      }

#if (TESTING == 1)
      int num_bad_values = 0;
#endif

      if (d_dim == tbox::Dimension(2)) {
#if (TESTING == 1)
         num_bad_values =
#endif
         SkeletonBoundaryUtilities2::checkBdryData(
            d_uval->getName(),
            patch,
            vdb->mapVariableAndContextToIndex(d_uval, getDataContext()), 0,
            ghost_width_to_check,
            bbox,
            bscalarcase,
            d_bdry_edge_uval[refbdryloc]);
      } else if (d_dim == tbox::Dimension(3)) {
#if (TESTING == 1)
         num_bad_values =
#endif
         SkeletonBoundaryUtilities3::checkBdryData(
            d_uval->getName(),
            patch,
            vdb->mapVariableAndContextToIndex(d_uval, getDataContext()), 0,
            ghost_width_to_check,
            bbox,
            bscalarcase,
            d_bdry_face_uval[refbdryloc]);
      }
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nMblkLinAdv Boundary Test FAILED: \n"
                    << "     " << num_bad_values
                    << " bad UVAL values found for\n"
                    << "     boundary type " << btype << " at location "
                    << bloc << std::endl;
      }
#endif

   }

}

hier::IntVector MblkLinAdv::getMultiblockRefineOpStencilWidth() const
{
   return hier::IntVector(d_dim, 1);
}

hier::IntVector MblkLinAdv::getMultiblockCoarsenOpStencilWidth()
{
   return hier::IntVector(d_dim, 0);
}
