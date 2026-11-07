/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Template for a multiblock AMR Euler code
 *
 ************************************************************************/
#include "MblkEuler.h"

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
#include <vector>

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

// Number of ghosts cells used for each variable quantity
#define CELLG (1)
#define FLUXG (0)
#define NODEG (1)

// defines for cell tagging routines
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

// Version of MblkEuler restart file data
#define MBLKEULER_VERSION (3)

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

// ================================= MblkEuler::Initialization =============================

/*
 *************************************************************************
 *
 * The constructor for MblkEuler class sets data members to defualt values,
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

MblkEuler::MblkEuler(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<hier::BaseGridGeometry>& grid_geom):
   MblkHyperbolicPatchStrategy(dim),
   d_object_name(object_name),
   d_dim(dim),
   d_grid_geometry(grid_geom),
   d_use_nonuniform_workload(false),
   d_nghosts(hier::IntVector(d_dim, CELLG)),
   d_fluxghosts(hier::IntVector(d_dim, FLUXG)),
   d_nodeghosts(hier::IntVector(d_dim, NODEG))
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   //
   // Setup MblkGeometry object to manage construction of mapped grids
   //
   d_mblk_geometry = new MblkGeometry("MblkGeometry",
         d_dim,
         input_db,
         grid_geom->getNumberBlocks());

   std::shared_ptr<tbox::Database> mbe_db(
      input_db->getDatabase("MblkEuler"));

   //
   // zero out initial condition variables
   //
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, d_dim.getValue());
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_axis, d_dim.getValue());

   //
   // Initialize object with data read from given input/restart databases.
   //
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   //
   // quantities that define state of linear advection problem.
   //
   d_state.reset(new pdat::CellVariable<double>(dim, "state", d_nState));
   d_vol.reset(new pdat::CellVariable<double>(dim, "vol", 1));
   d_flux.reset(new pdat::SideVariable<double>(dim, "flux",
         hier::IntVector::getOne(dim), d_nState));
   d_xyz.reset(new pdat::NodeVariable<double>(dim, "xyz", d_dim.getValue()));

   //
   // drop the region layout as a table
   //
   tbox::plog << "region layout follows:" << std::endl;

   tbox::plog << "field";
   for (int ir = 0; ir < d_number_of_regions; ++ir)
      tbox::plog << "\t" << ir;
   tbox::plog << std::endl;

   for (int ii = 0; ii < d_number_of_regions * 10; ++ii)
      tbox::plog << "-";
   tbox::plog << std::endl;

   for (int istate = 0; istate < d_nState; ++istate) {
      tbox::plog << d_state_names[istate];
      for (int ir = 0; ir < d_number_of_regions; ++ir)
         tbox::plog << "\t" << d_state_ic[ir][istate];
      tbox::plog << std::endl;
   }
}

/*
 *************************************************************************
 *
 * Empty destructor for MblkEuler class.
 *
 *************************************************************************
 */

MblkEuler::~MblkEuler() {
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

void MblkEuler::registerModelVariables(
   MblkHyperbolicLevelIntegrator* integrator)
{
   TBOX_ASSERT(integrator != 0);

   //
   // zonal data and its fluxes
   //
   d_cell_time_interp_op.reset(
      new pdat::CellDoubleLinearTimeInterpolateOp());

   integrator->registerVariable(d_state, d_nghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      std::shared_ptr<hier::CoarsenOperator>(),
      std::shared_ptr<hier::RefineOperator>(),
      d_cell_time_interp_op);

   integrator->registerVariable(d_vol, d_nghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      std::shared_ptr<hier::CoarsenOperator>(),
      std::shared_ptr<hier::RefineOperator>(),
      d_cell_time_interp_op);

   integrator->registerVariable(d_flux, d_fluxghosts,
      MblkHyperbolicLevelIntegrator::FLUX,
      std::shared_ptr<hier::CoarsenOperator>());

   //
   // The nodal position data
   //
   std::shared_ptr<hier::TimeInterpolateOperator> node_time_interp_op(
      new pdat::NodeDoubleLinearTimeInterpolateOp());

   integrator->registerVariable(d_xyz, d_nodeghosts,
      MblkHyperbolicLevelIntegrator::TIME_DEP,
      std::shared_ptr<hier::CoarsenOperator>(),
      std::shared_ptr<hier::RefineOperator>(),
      node_time_interp_op);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

#ifdef HAVE_HDF5
   if (!d_visit_writer) {
      TBOX_WARNING(
         d_object_name << ": registerModelVariables()"
                       << "\nVisIt data writer was"
                       << "\nregistered.  Consequently, no plot data will"
                       << "\nbe written." << std::endl);
   }

   for (int n = 0; n < d_nState; ++n) {
      std::string vname = d_state_names[n];
      d_visit_writer->registerPlotQuantity(vname, "SCALAR",
         vardb->mapVariableAndContextToIndex(d_state,
            integrator->getPlotContext()),
         n);
   }

   d_visit_writer->registerPlotQuantity("vol", "SCALAR",
      vardb->mapVariableAndContextToIndex(d_vol, integrator->getPlotContext()));

   d_visit_writer->registerNodeCoordinates(vardb->mapVariableAndContextToIndex(
         d_xyz, integrator->getPlotContext()));
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
void MblkEuler::initializeDataOnPatch(
   hier::Patch& patch,
   const double data_time,
   const bool initial_time)
{
   NULL_USE(data_time);

   //
   // Build the mapped grid on the patch.
   //
   setMappedGridOnPatch(patch);

   if (initial_time) {

      std::shared_ptr<pdat::CellData<double> > state(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_state, getDataContext())));
      std::shared_ptr<pdat::CellData<double> > vol(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_vol, getDataContext())));
      std::shared_ptr<pdat::NodeData<double> > xyz(
         SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
            patch.getPatchData(d_xyz, getDataContext())));

      TBOX_ASSERT(state);
      TBOX_ASSERT(vol);
      TBOX_ASSERT(xyz);
      TBOX_ASSERT(state->getGhostCellWidth() == vol->getGhostCellWidth());

      hier::IntVector state_ghosts = state->getGhostCellWidth();
      hier::IntVector xyz_ghosts = xyz->getGhostCellWidth();

      const hier::Index ifirst = patch.getBox().lower();
      const hier::Index ilast = patch.getBox().upper();

      int imin = ifirst(0) - state_ghosts(0);
      int imax = ilast(0) + state_ghosts(0);
      int jmin = ifirst(1) - state_ghosts(1);
      int jmax = ilast(1) + state_ghosts(1);
      int kmin = ifirst(2) - state_ghosts(2);
      int kmax = ilast(2) + state_ghosts(2);
      int nx = imax - imin + 1;
      int ny = jmax - jmin + 1;
      int nxny = nx * ny;

      int nd_imin = ifirst(0) - xyz_ghosts(0);
      int nd_imax = ilast(0) + 1 + xyz_ghosts(0);
      int nd_jmin = ifirst(1) - xyz_ghosts(1);
      int nd_jmax = ilast(1) + 1 + xyz_ghosts(1);
      int nd_kmin = ifirst(2) - xyz_ghosts(2);
      int nd_nx = nd_imax - nd_imin + 1;
      int nd_ny = nd_jmax - nd_jmin + 1;
      int nd_nxny = nd_nx * nd_ny;

      //
      // get the pointers
      //
      double* cvol = vol->getPointer();

      double* x = xyz->getPointer(0);
      double* y = xyz->getPointer(1);
      double* z = xyz->getPointer(2);

      hier::IntVector ghost_cells = state->getGhostCellWidth();
      pdat::CellData<double> elemCoords(patch.getBox(), 3, ghost_cells);         // storage for the average of the element coordinates
      pdat::CellData<int> region_ids_data(patch.getBox(), 1, ghost_cells);       // storage for the slopes

      int* region_ids = region_ids_data.getPointer();
      double* xc = elemCoords.getPointer(0);
      double* yc = elemCoords.getPointer(1);
      double* zc = elemCoords.getPointer(2);

      //
      //  ---------------- compute the element coordinates
      //
      for (int k = kmin; k <= kmax; ++k) {
         for (int j = jmin; j <= jmax; ++j) {
            for (int i = imin; i <= imax; ++i) {
               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

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

               double xavg = 0.125 * (x[n1] + x[n2] + x[n3] + x[n4]
                                      + x[n5] + x[n6] + x[n7] + x[n8]);

               double yavg = 0.125 * (y[n1] + y[n2] + y[n3] + y[n4]
                                      + y[n5] + y[n6] + y[n7] + y[n8]);

               double zavg = 0.125 * (z[n1] + z[n2] + z[n3] + z[n4]
                                      + z[n5] + z[n6] + z[n7] + z[n8]);
               xc[ind] = xavg;
               yc[ind] = yavg;
               zc[ind] = zavg;
            }
         }
      }

      //
      //  ---------------- compute the element volume
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

               double lvol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                     x[n5], x[n6], x[n7], x[n8],

                     y[n1], y[n2], y[n3], y[n4],
                     y[n5], y[n6], y[n7], y[n8],

                     z[n1], z[n2], z[n3], z[n4],
                     z[n5], z[n6], z[n7], z[n8]);

               cvol[cind] = lvol;
            }
         }
      }

      //
      //  ---------------- process the different initialization regions
      //
      if (d_data_problem == "REVOLUTION") {

         for (int m = 0; m < d_number_of_regions; ++m) {    // loop over the regions and shape in data

            std::vector<double>& lrad = d_rev_rad[m];
            std::vector<double>& laxis = d_rev_axis[m];
            int naxis = static_cast<int>(laxis.size());

            for (int k = kmin; k <= kmax; ++k) {
               for (int j = jmin; j <= jmax; ++j) {
                  for (int i = imin; i <= imax; ++i) {
                     int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

                     double x2 = zc[ind] - d_center[2];
                     double x1 = yc[ind] - d_center[1];
                     double x0 = xc[ind] - d_center[0];

                     double Dz = x0 * d_axis[0] + x1 * d_axis[1] + x2
                        * d_axis[2];
                     double r2 = x0 * x0 + x1 * x1 + x2 * x2;
                     double Dr = sqrt(r2 - Dz * Dz);

                     if (laxis[0] <= Dz && Dz <= laxis[naxis - 1]) { // test to see if we are contained

                        int lpos = 0;
                        while (Dz > laxis[lpos])
                           ++lpos;

                        double a =
                           (Dz
                            - laxis[lpos - 1]) / (laxis[lpos] - laxis[lpos - 1]);
                        double lr = (a) * lrad[lpos]
                           + (1.0 - a) * lrad[lpos - 1];

                        if (Dr <= lr) { // if we are within the radius set the region id
                           region_ids[ind] = m;
                        }

                     } // laxis

                  }
               }
            } // k

         } // end of region loop
      }
      //
      //  the spherical initialization
      //
      else if (d_data_problem == "SPHERE") {
         double* front = &d_front_position[0];
         for (int k = kmin; k <= kmax; ++k) {
            for (int j = jmin; j <= jmax; ++j) {
               for (int i = imin; i <= imax; ++i) {
                  int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

                  double x2 = zc[ind] - d_center[2];
                  double x1 = yc[ind] - d_center[1];
                  double x0 = xc[ind] - d_center[0];

                  double rad2 = sqrt(x0 * x0 + x1 * x1);
                  // double phi  = atan2(rad2,x2);
                  double rad3 = sqrt(rad2 * rad2 + x2 * x2);

                  int ifr = 0;  // find the region we draw from (ifr=0 is always origin)
                  while (rad3 > front[ifr + 1]) {
                     ++ifr;
                  }
                  region_ids[ind] = ifr;
               }
            }
         }
      }
      //
      //  the planar initialization
      //
      else if ((d_data_problem == "PIECEWISE_CONSTANT_X")
               || (d_data_problem == "PIECEWISE_CONSTANT_Y")
               || (d_data_problem == "PIECEWISE_CONSTANT_Z")) {

         double* front = &d_front_position[0];
         double* xx = xc;
         if (d_data_problem == "PIECEWISE_CONSTANT_Y") {
            xx = yc;
         }
         if (d_data_problem == "PIECEWISE_CONSTANT_Z") {
            xx = zc;
         }

         for (int k = kmin; k <= kmax; ++k) {
            for (int j = jmin; j <= jmax; ++j) {
               for (int i = imin; i <= imax; ++i) {
                  int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

                  int ifr = 0;   // the geion we are in
                  while (xx[ind] > front[ifr + 1]) {
                     ++ifr;
                  }
                  region_ids[ind] = ifr;
               }
            }
         }
      }
      //
      //  the planar initialization
      //
      else if ((d_data_problem == "RT_SHOCK_TUBE")) {

         double* front = &d_front_position[0];
         double shock_pos = front[1]; // the shock front
         double front_pos = front[2]; // the sinusoidal perturbation between the two fluids

         double dt_ampl = d_dt_ampl;
         int nmodes = static_cast<int>(d_amn.size());
         double* amn = &d_amn[0];
         double* n_mode = &d_n_mode[0];
         double* m_mode = &d_m_mode[0];
         double* phiy = &d_phiy[0];
         double* phiz = &d_phiz[0];

         // ... this is a cartesian problem by definition
         const double* xdlo = 0;  // d_cart_xlo[0][0];
         const double* xdhi = 0;  // d_cart_xhi[0][0];
         TBOX_ASSERT(0);  // xdlo, xdhi wrong

         double l_y = xdhi[1] - xdlo[1];      // width of the domain in y and z
         double l_z = xdhi[2] - xdlo[2];
         double lpi = 3.14159265358979310862446895044;

         for (int k = kmin; k <= kmax; ++k) {
            for (int j = jmin; j <= jmax; ++j) {
               for (int i = imin; i <= imax; ++i) {
                  int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

                  double lpert = 0.0;
                  double ly = yc[ind];
                  double lz = zc[ind];
                  for (int m = 0; m < nmodes; ++m) {
                     double lphiy = 2.0 * lpi * n_mode[m] * ly / l_y + phiy[m];
                     double lphiz = 2.0 * lpi * m_mode[m] * lz / l_z + phiz[m];
                     double cy = cos(lphiy);
                     double cz = cos(lphiz);
                     lpert += amn[m] * cy * cz;
                  }
                  lpert *= dt_ampl;

                  int ifr = 0;
                  if (xc[ind] > shock_pos)
                     ifr = 1;
                  if (xc[ind] > front_pos + lpert)
                     ifr = 2;

                  region_ids[ind] = ifr;
               }
            }
         }
      } else if (d_data_problem == "BLIP") {
         TBOX_ASSERT(0);
      }

      //
      // ---------------- state vector
      //
      int depth = state->getDepth();

      for (int idepth = 0; idepth < depth; ++idepth) {
         double* psi = state->getPointer(idepth);

         for (int k = ifirst(2); k <= ilast(2); ++k) {
            for (int j = ifirst(1); j <= ilast(1); ++j) {
               for (int i = ifirst(0); i <= ilast(0); ++i) {

                  int cind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
#if 0
                  // psi[cind] = (double) i;  // set values to coordinate ids as debug tool
                  // psi[cind] = (double) j;
                  // psi[cind] = (double) k;
#endif
#if 1
                  int ireg = region_ids[cind];
                  psi[cind] = d_state_ic[ireg][idepth];
#endif
               }
            }
         }

      } // end of depth loop

   } // end of initial time if test

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

// ================================= MblkEuler::Integration =============================

/*
 *************************************************************************
 *
 * Compute stable time increment for patch.  Return this value.
 *
 *************************************************************************
 */

double MblkEuler::computeStableDtOnPatch(
   hier::Patch& patch,
   const bool initial_time,
   const double dt_time)
{
   NULL_USE(initial_time);
   NULL_USE(dt_time);

   //
   // Build the mapped grid on the patch.
   //
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   tbox::plog << "--------------------- start stableDtOnPatch on patch ("
              << level_number << ")" << std::endl;
   tbox::plog << "level = " << level_number << std::endl;
   tbox::plog << "box   = " << patch.getBox() << std::endl;

   //
   // get the cell data and their bounds
   //
   std::shared_ptr<pdat::CellData<double> > state(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));
   hier::IntVector state_ghosts = state->getGhostCellWidth();
   hier::IntVector xyz_ghosts = xyz->getGhostCellWidth();

   pdat::CellData<double> Aii(patch.getBox(), 9, hier::IntVector(d_dim, 0));

   TBOX_ASSERT(state);
   TBOX_ASSERT(xyz);

   int imin = ifirst(0);
   int imax = ilast(0);
   int jmin = ifirst(1);
   int jmax = ilast(1);
   int kmin = ifirst(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nxny = nx * ny;

   int nd_imin = ifirst(0) - xyz_ghosts(0);
   int nd_imax = ilast(0) + 1 + xyz_ghosts(0);
   int nd_jmin = ifirst(1) - xyz_ghosts(1);
   int nd_jmax = ilast(1) + 1 + xyz_ghosts(1);
   int nd_kmin = ifirst(2) - xyz_ghosts(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   int pind = 0;
   double* a11 = Aii.getPointer(pind);
   ++pind;
   double* a12 = Aii.getPointer(pind);
   ++pind;
   double* a13 = Aii.getPointer(pind);
   ++pind;
   double* a21 = Aii.getPointer(pind);
   ++pind;
   double* a22 = Aii.getPointer(pind);
   ++pind;
   double* a23 = Aii.getPointer(pind);
   ++pind;
   double* a31 = Aii.getPointer(pind);
   ++pind;
   double* a32 = Aii.getPointer(pind);
   ++pind;
   double* a33 = Aii.getPointer(pind);
   ++pind;

   //
   // compute direction cosines
   //
   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;
            int n7 = n1 + nd_nxny + 1 + nd_nx;
            int n8 = n1 + nd_nxny + nd_nx;

            // ------------------------------------------------ x

            double x1 = 0.25 * (x[n1] + x[n4] + x[n5] + x[n8]);  // xi
            double x2 = 0.25 * (x[n2] + x[n3] + x[n6] + x[n7]);

            double x3 = 0.25 * (x[n1] + x[n2] + x[n5] + x[n6]);  // eta
            double x4 = 0.25 * (x[n3] + x[n4] + x[n7] + x[n8]);

            double x5 = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4]);  // zeta
            double x6 = 0.25 * (x[n5] + x[n6] + x[n7] + x[n8]);

            // ------------------------------------------------ y

            double y1 = 0.25 * (y[n1] + y[n4] + y[n5] + y[n8]);
            double y2 = 0.25 * (y[n2] + y[n3] + y[n6] + y[n7]);

            double y3 = 0.25 * (y[n1] + y[n2] + y[n5] + y[n6]);
            double y4 = 0.25 * (y[n3] + y[n4] + y[n7] + y[n8]);

            double y5 = 0.25 * (y[n1] + y[n2] + y[n3] + y[n4]);
            double y6 = 0.25 * (y[n5] + y[n6] + y[n7] + y[n8]);

            // ------------------------------------------------ z

            double z1 = 0.25 * (z[n1] + z[n4] + z[n5] + z[n8]);
            double z2 = 0.25 * (z[n2] + z[n3] + z[n6] + z[n7]);

            double z3 = 0.25 * (z[n1] + z[n2] + z[n5] + z[n6]);
            double z4 = 0.25 * (z[n3] + z[n4] + z[n7] + z[n8]);

            double z5 = 0.25 * (z[n1] + z[n2] + z[n3] + z[n4]);
            double z6 = 0.25 * (z[n5] + z[n6] + z[n7] + z[n8]);

            //
            // the components of the matrices that we want to invert
            //
            double dx_xi = 0.5 * (x2 - x1);
            double dy_xi = 0.5 * (y2 - y1);
            double dz_xi = 0.5 * (z2 - z1);

            double dx_eta = 0.5 * (x4 - x3);
            double dy_eta = 0.5 * (y4 - y3);
            double dz_eta = 0.5 * (z4 - z3);

            double dx_zeta = 0.5 * (x6 - x5);
            double dy_zeta = 0.5 * (y6 - y5);
            double dz_zeta = 0.5 * (z6 - z5);

            //
            // invert M = dx/dxi to find the matrix needed to convert
            // displacements in x into displacements in xi (dxi/dx)
            // via kramer's rule
            //
            double detM = (dx_xi * dy_eta * dz_zeta
                           + dx_eta * dy_zeta * dz_xi
                           + dx_zeta * dy_xi * dz_eta
                           - dx_zeta * dy_eta * dz_xi
                           - dx_eta * dy_xi * dz_zeta
                           - dx_xi * dy_zeta * dz_eta);

            double detB11 = dy_eta * dz_zeta - dy_zeta * dz_eta;
            double detB21 = dy_zeta * dz_xi - dy_xi * dz_zeta;
            double detB31 = dy_xi * dz_eta - dy_eta * dz_xi;

            double detB12 = dx_zeta * dz_eta - dx_eta * dz_zeta;
            double detB22 = dx_xi * dz_zeta - dx_zeta * dz_xi;
            double detB32 = dx_eta * dz_xi - dx_xi * dz_eta;

            double detB13 = dx_eta * dy_zeta - dx_zeta * dy_eta;
            double detB23 = dx_zeta * dy_xi - dx_xi * dy_zeta;
            double detB33 = dx_xi * dy_eta - dx_eta * dy_xi;

            a11[ind] = detB11 / detM;
            a21[ind] = detB21 / detM;
            a31[ind] = detB31 / detM;

            a12[ind] = detB12 / detM;
            a22[ind] = detB22 / detM;
            a32[ind] = detB32 / detM;

            a13[ind] = detB13 / detM;
            a23[ind] = detB23 / detM;
            a33[ind] = detB33 / detM;
         }
      }
   }

   //
   // print out patch extrema
   //
   testPatchExtrema(patch, "before timestep evaluation");

   //
   // the timestep evaluation
   //
   double stabdt = 1.e20;
   if (d_advection_test) {
      // ------------------------------------- for linear advection
      double u = d_advection_velocity[0];
      double v = d_advection_velocity[1];
      double w = d_advection_velocity[2];

      for (int k = ifirst(2); k <= ilast(2); ++k) {
         for (int j = ifirst(1); j <= ilast(1); ++j) {
            for (int i = ifirst(0); i <= ilast(0); ++i) {
               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

               double uxi = a11[ind] * u + a12[ind] * v + a13[ind] * w;  // max parametric signal speed
               double ueta = a21[ind] * u + a22[ind] * v + a23[ind] * w;
               double uzeta = a31[ind] * u + a32[ind] * v + a33[ind] * w;

               uxi = fabs(uxi);
               ueta = fabs(ueta);
               uzeta = fabs(uzeta);

               double dtxi = 2.0 / MAX(1.e-80, uxi);
               double dteta = 2.0 / MAX(1.e-80, ueta);
               double dtzeta = 2.0 / MAX(1.e-80, uzeta);

               double ldelt = MIN(dtxi, MIN(dteta, dtzeta));
               stabdt = MIN(stabdt, ldelt);

            }
         }
      }
   } else {
      TBOX_ASSERT(d_advection_test);
   }

   //
   // process the timestep constraints
   //
   double dt_fixed = -1;
   double returned_dt = stabdt;
   if (dt_fixed > 0.0)
      returned_dt = dt_fixed;

   tbox::plog << "stabdt      = " << stabdt << std::endl;
   tbox::plog << "dt_fixed    = " << dt_fixed << std::endl;
   tbox::plog << "returned_dt = " << returned_dt << std::endl;

   tbox::plog << "--------------------- end stableDtOnPatch on patch" << std::endl;

   return returned_dt;
}

// ---------------------------------------------------------------------------

//
// Test the extrema on the patch
//
void MblkEuler::testPatchExtrema(
   hier::Patch& patch,
   const char* pos)
{
   std::shared_ptr<pdat::CellData<double> > state(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_state, getDataContext())));
   TBOX_ASSERT(state);
   hier::IntVector state_ghosts = state->getGhostCellWidth();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   int imgn = ifirst(0) - state_ghosts(0);
   int imgx = ilast(0) + state_ghosts(0);
   int jmgn = ifirst(1) - state_ghosts(1);
   int jmgx = ilast(1) + state_ghosts(1);
   int kmgn = ifirst(2) - state_ghosts(2);
   int nxg = imgx - imgn + 1;
   int nyg = jmgx - jmgn + 1;
   int nxnyg = nxg * nyg;

   //
   //  compute field extrema to tag along with this and print out
   //
   double* psi_min = new double[d_nState];
   double* psi_max = new double[d_nState];
   for (int ii = 0; ii < d_nState; ++ii) {
      psi_max[ii] = -1.e80;
      psi_min[ii] = 1.e80;
   }

   for (int ii = 0; ii < d_nState; ++ii) {
      double* lstate = state->getPointer(ii);

      for (int k = ifirst(2); k <= ilast(2); ++k) {     // just loop over interior elements
         for (int j = ifirst(1); j <= ilast(1); ++j) {
            for (int i = ifirst(0); i <= ilast(0); ++i) {

               int gind = POLY3(i, j, k, imgn, jmgn, kmgn, nxg, nxnyg);

               //
               // some extra bounds checks for sanity
               //
               psi_max[ii] = MAX(psi_max[ii], lstate[gind]);
               psi_min[ii] = MIN(psi_min[ii], lstate[gind]);
            }
         }
      }
   }

   tbox::plog << std::endl << "extrema for the state follow " << pos
              << " (min,max) = " << std::endl;
   for (int ii = 0; ii < d_nState; ++ii) {
      tbox::plog << d_state_names[ii] << " (min,max) = ";
      tbox::plog << psi_min[ii] << " " << psi_max[ii] << std::endl;
   }

   delete[] psi_min;
   delete[] psi_max;
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

void MblkEuler::computeFluxesOnPatch(
   hier::Patch& patch,
   const double time,
   const double dt)
{
   //
   // process the SAMRAI data
   //
   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch);

   hier::Box pbox = patch.getBox();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > state(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > vol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_vol, getDataContext())));
   std::shared_ptr<pdat::SideData<double> > flux(
      SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
         patch.getPatchData(d_flux, getDataContext())));
   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));

   TBOX_ASSERT(state);
   TBOX_ASSERT(flux);
   TBOX_ASSERT(vol);
   TBOX_ASSERT(xyz);
   TBOX_ASSERT(state->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(vol->getGhostCellWidth() == d_nghosts);
   TBOX_ASSERT(flux->getGhostCellWidth() == d_fluxghosts);
   TBOX_ASSERT(xyz->getGhostCellWidth() == d_nodeghosts);

   tbox::plog << "--------------------- start computeFluxesOnPatch on patch (";
   tbox::plog << level_number << ")" << std::endl;
   tbox::plog << "TIMESTEP for level = " << level_number << ", ";
   tbox::plog << "dt    = " << dt << ", stime  = " << time << ", ftime = "
              << time + dt << std::endl;
   tbox::plog << "box   = " << pbox << std::endl;

   //
   // ------------------------------- the upwind bounds ----------------------------
   //

   int fx_imn = ifirst(0) - d_fluxghosts(0);
   int fx_imx = ilast(0) + 1 + d_fluxghosts(0);
   int fx_jmn = ifirst(1) - d_fluxghosts(1);
   int fx_jmx = ilast(1) + d_fluxghosts(1);
   int fx_kmn = ifirst(2) - d_fluxghosts(2);
   int fx_nx = fx_imx - fx_imn + 1;
   int fx_ny = fx_jmx - fx_jmn + 1;
   int fx_nxny = fx_nx * fx_ny;

   int fy_imn = ifirst(0) - d_fluxghosts(0);
   int fy_imx = ilast(0) + d_fluxghosts(0);
   int fy_jmn = ifirst(1) - d_fluxghosts(1);
   int fy_jmx = ilast(1) + 1 + d_fluxghosts(1);
   int fy_kmn = ifirst(2) - d_fluxghosts(2);
   int fy_nx = fy_imx - fy_imn + 1;
   int fy_ny = fy_jmx - fy_jmn + 1;
   int fy_nxny = fy_nx * fy_ny;

   int fz_imn = ifirst(0) - d_fluxghosts(0);
   int fz_imx = ilast(0) + d_fluxghosts(0);
   int fz_jmn = ifirst(1) - d_fluxghosts(1);
   int fz_jmx = ilast(1) + d_fluxghosts(1);
   int fz_kmn = ifirst(2) - d_fluxghosts(2);
   int fz_nx = fz_imx - fz_imn + 1;
   int fz_ny = fz_jmx - fz_jmn + 1;
   int fz_nxny = fz_nx * fz_ny;

   int imin = ifirst(0) - d_nghosts(0);
   int imax = ilast(0) + d_nghosts(0);
   int jmin = ifirst(1) - d_nghosts(1);
   int jmax = ilast(1) + d_nghosts(1);
   int kmin = ifirst(2) - d_nghosts(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nxny = nx * ny;

   int nd_imin = ifirst(0) - d_nodeghosts(0);
   int nd_imax = ilast(0) + 1 + d_nodeghosts(0);
   int nd_jmin = ifirst(1) - d_nodeghosts(1);
   int nd_jmax = ilast(1) + 1 + d_nodeghosts(1);
   int nd_kmin = ifirst(2) - d_nodeghosts(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   //
   // get the pointers
   //
   double* cvol = vol->getPointer();

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   double u = d_advection_velocity[0];
   double v = d_advection_velocity[1];
   double w = d_advection_velocity[2];
   double u0 = sqrt(u * u + v * v + w * w);

   // 0, cartesian, 1 R translation, 2 rigid body theta rotation, 4, rigid body phi rotation
   int VEL_TYPE = d_advection_vel_type;

   //
   // note on areas, we set the area vector Avector,
   //
   //     Avector = -0.5*( x3 - x1 ) cross product ( x4 - x2 )
   //
   // where the x1 .. x4 are the right hand rule circulation for the face
   //

   int depth = state->getDepth();

   for (int idepth = 0; idepth < depth; ++idepth) {

      double* psi = state->getPointer(idepth); // assumed single depth here !!!!

      double* fx = flux->getPointer(0, idepth);
      double* fy = flux->getPointer(1, idepth);
      double* fz = flux->getPointer(2, idepth);

      //
      // compute the fluxes for the upwind method
      //
      for (int k = ifirst(2); k <= ilast(2); ++k) {
         for (int j = ifirst(1); j <= ilast(1); ++j) {
            for (int i = ifirst(0); i <= ilast(0) + 1; ++i) {

               // --------- get the indices
               int ifx = POLY3(i, j, k, fx_imn, fx_jmn, fx_kmn, fx_nx, fx_nxny);

               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
               int ib = ind - 1;

               int n1 = POLY3(i,
                     j,
                     k,
                     nd_imin,
                     nd_jmin,
                     nd_kmin,
                     nd_nx,
                     nd_nxny);
               int n4 = n1 + nd_nx;
               int n5 = n1 + nd_nxny;
               int n8 = n1 + nd_nxny + nd_nx;

               // --------- get the positions // 1 - 4 - 8 - 5
               double x1 = x[n1];
               double y1 = y[n1];
               double z1 = z[n1];
               double x2 = x[n4];
               double y2 = y[n4];
               double z2 = z[n4];
               double x3 = x[n8];
               double y3 = y[n8];
               double z3 = z[n8];
               double x4 = x[n5];
               double y4 = y[n5];
               double z4 = z[n5];

               double xm = 0.25 * (x1 + x2 + x3 + x4);
               double ym = 0.25 * (y1 + y2 + y3 + y4);
               double zm = 0.25 * (z1 + z2 + z3 + z4);

               double xn1 = sqrt(xm * xm + ym * ym);
               double R = sqrt(xm * xm + ym * ym + zm * zm);

               // --------- compute the flux
               double u1 = u0 * xm / R;
               double v1 = u0 * ym / R;
               double w1 = u0 * zm / R;

               double dx31 = x3 - x1;
               double dx42 = x4 - x2;

               double dy31 = y3 - y1;
               double dy42 = y4 - y2;

               double dz31 = z3 - z1;
               double dz42 = z4 - z2;

               double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
               double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
               double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

               double Audotn = 0.0;
               if (VEL_TYPE == 0) {
                  Audotn = Ax * u + Ay * v + Az * w;
               }
               if (VEL_TYPE == 1) {
                  Audotn = Ax * u1 + Ay * v1 + Az * w1;
               }
               if (VEL_TYPE == 2) {
                  double u2 = -u0 * ym;
                  double v2 = u0 * xm;
                  double w2 = 0.0;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }
               if (VEL_TYPE == 3) {
                  double cosphi = xn1 / R;
                  double sinphi =
                     sqrt(1 - xn1 * xn1 / (R * R)) * (zm > 0.0 ? 1.0 : -1.0);
                  double u2 = -u0 * xm * sinphi / cosphi;
                  double v2 = -u0 * ym * sinphi / cosphi;
                  double w2 = u0 * cosphi / sinphi;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }

               double lflux = (Audotn > 0.0 ? psi[ib] : psi[ind]) * Audotn;

               fx[ifx] = lflux;  // (Ax,Ay,Az) vector points towards ind
            }
         }
      }

      for (int k = ifirst(2); k <= ilast(2); ++k) {
         for (int j = ifirst(1); j <= ilast(1) + 1; ++j) {
            for (int i = ifirst(0); i <= ilast(0); ++i) {

               // --------- get the indices
               int ify = POLY3(i, j, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny);

               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
               int jb = ind - nx;

               int n1 = POLY3(i,
                     j,
                     k,
                     nd_imin,
                     nd_jmin,
                     nd_kmin,
                     nd_nx,
                     nd_nxny);
               int n2 = n1 + 1;
               int n5 = n1 + nd_nxny;
               int n6 = n1 + nd_nxny + 1;

               // --------- get the positions // 1 - 5 - 6 - 2
               double x1 = x[n1];
               double y1 = y[n1];
               double z1 = z[n1];
               double x2 = x[n5];
               double y2 = y[n5];
               double z2 = z[n5];
               double x3 = x[n6];
               double y3 = y[n6];
               double z3 = z[n6];
               double x4 = x[n2];
               double y4 = y[n2];
               double z4 = z[n2];

               double xm = 0.25 * (x1 + x2 + x3 + x4);
               double ym = 0.25 * (y1 + y2 + y3 + y4);
               double zm = 0.25 * (z1 + z2 + z3 + z4);

               double xn1 = sqrt(xm * xm + ym * ym);
               double R = sqrt(xm * xm + ym * ym + zm * zm);

               // --------- compute the flux
               double u1 = u0 * xm / R;
               double v1 = u0 * ym / R;
               double w1 = u0 * zm / R;

               double dx31 = x3 - x1;
               double dx42 = x4 - x2;

               double dy31 = y3 - y1;
               double dy42 = y4 - y2;

               double dz31 = z3 - z1;
               double dz42 = z4 - z2;

               double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
               double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
               double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

               double Audotn = 0.0;
               if (VEL_TYPE == 0) {
                  Audotn = Ax * u + Ay * v + Az * w;
               }
               if (VEL_TYPE == 1) {
                  Audotn = Ax * u1 + Ay * v1 + Az * w1;
               }
               if (VEL_TYPE == 2) {
                  double u2 = -u0 * ym;
                  double v2 = u0 * xm;
                  double w2 = 0.0;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }
               if (VEL_TYPE == 3) {
                  double cosphi = xn1 / R;
                  double sinphi =
                     sqrt(1 - xn1 * xn1 / (R * R)) * (zm > 0.0 ? 1.0 : -1.0);
                  double u2 = -u0 * xm * sinphi / cosphi;
                  double v2 = -u0 * ym * sinphi / cosphi;
                  double w2 = u0 * cosphi / sinphi;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }

               double lflux = (Audotn > 0.0 ? psi[jb] : psi[ind]) * Audotn;

               fy[ify] = lflux;  // (Ax,Ay,Az) vector points towards ind
            }
         }
      }

      for (int k = ifirst(2); k <= ilast(2) + 1; ++k) {
         for (int j = ifirst(1); j <= ilast(1); ++j) {
            for (int i = ifirst(0); i <= ilast(0); ++i) {

               // --------- get the indices
               int ifz = POLY3(i, j, k, fz_imn, fz_jmn, fz_kmn, fz_nx, fz_nxny);

               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
               int kb = ind - nxny;

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

               // --------- get the positions // 1 - 2 - 3 - 4
               double x1 = x[n1];
               double y1 = y[n1];
               double z1 = z[n1];
               double x2 = x[n2];
               double y2 = y[n2];
               double z2 = z[n2];
               double x3 = x[n3];
               double y3 = y[n3];
               double z3 = z[n3];
               double x4 = x[n4];
               double y4 = y[n4];
               double z4 = z[n4];

               double xm = 0.25 * (x1 + x2 + x3 + x4);
               double ym = 0.25 * (y1 + y2 + y3 + y4);
               double zm = 0.25 * (z1 + z2 + z3 + z4);

               double xn1 = sqrt(xm * xm + ym * ym);
               double R = sqrt(xm * xm + ym * ym + zm * zm);

               // --------- compute the flux
               double u1 = u0 * xm / R;
               double v1 = u0 * ym / R;
               double w1 = u0 * zm / R;

               double dx31 = x3 - x1;
               double dx42 = x4 - x2;

               double dy31 = y3 - y1;
               double dy42 = y4 - y2;

               double dz31 = z3 - z1;
               double dz42 = z4 - z2;

               double Ax = -0.5 * (dy42 * dz31 - dz42 * dy31);
               double Ay = -0.5 * (dz42 * dx31 - dx42 * dz31);
               double Az = -0.5 * (dx42 * dy31 - dy42 * dx31);

               double Audotn = 0.0;
               if (VEL_TYPE == 0) {
                  Audotn = Ax * u + Ay * v + Az * w;
               }
               if (VEL_TYPE == 1) {
                  Audotn = Ax * u1 + Ay * v1 + Az * w1;
               }
               if (VEL_TYPE == 2) {
                  double u2 = -u0 * ym;
                  double v2 = u0 * xm;
                  double w2 = 0.0;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }
               if (VEL_TYPE == 3) {
                  double cosphi = xn1 / R;
                  double sinphi =
                     sqrt(1 - xn1 * xn1 / (R * R)) * (zm > 0.0 ? 1.0 : -1.0);
                  double u2 = -u0 * xm * sinphi / cosphi;
                  double v2 = -u0 * ym * sinphi / cosphi;
                  double w2 = u0 * cosphi / sinphi;
                  Audotn = Ax * u2 + Ay * v2 + Az * w2;
               }

               double lflux = (Audotn > 0.0 ? psi[kb] : psi[ind]) * Audotn;

               fz[ifz] = lflux;  // (Ax,Ay,Az) vector points towards ind
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
               int ie = POLY3(i + 1,
                     j,
                     k,
                     fx_imn,
                     fx_jmn,
                     fx_kmn,
                     fx_nx,
                     fx_nxny);

               int jb = POLY3(i, j, k, fy_imn, fy_jmn, fy_kmn, fy_nx, fy_nxny);
               int je = POLY3(i,
                     j + 1,
                     k,
                     fy_imn,
                     fy_jmn,
                     fy_kmn,
                     fy_nx,
                     fy_nxny);

               int kb = POLY3(i, j, k, fz_imn, fz_jmn, fz_kmn, fz_nx, fz_nxny);
               int ke = POLY3(i,
                     j,
                     k + 1,
                     fz_imn,
                     fz_jmn,
                     fz_kmn,
                     fz_nx,
                     fz_nxny);

               //   have set up the normal so that it always points in the positive, i, j, and k directions
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

               double lvol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                     x[n5], x[n6], x[n7], x[n8],

                     y[n1], y[n2], y[n3], y[n4],
                     y[n5], y[n6], y[n7], y[n8],

                     z[n1], z[n2], z[n3], z[n4],
                     z[n5], z[n6], z[n7], z[n8]);

               cvol[ind] = lvol;

               psi[ind] -= dt * ((fx[ie] - fx[ib])
                                 + (fy[je] - fy[jb])
                                 + (fz[ke] - fz[kb])) / lvol;

            }
         }
      }

   } // end of field loop

   tbox::plog << "--------------------- end computeFluxesOnPatch on patch"
              << std::endl;
}

/*
 *************************************************************************
 *
 * Update solution variables by performing a conservative
 * difference with the fluxes calculated in computeFluxesOnPatch().
 *
 *************************************************************************
 */

void MblkEuler::conservativeDifferenceOnPatch(
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
 * this routine initializes all the zonal data that one needs to fill to dummy
 * values as a check
 *
 *************************************************************************
 */

void MblkEuler::markPhysicalBoundaryConditions(
   hier::Patch& patch,
   const hier::IntVector& ghost_width_to_fill)
{
   std::shared_ptr<pdat::CellData<double> > state(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_state, getDataContext())));
   TBOX_ASSERT(state);

   //
   // the domain and its ghost box
   //
   const hier::Box& interior = patch.getBox();

   const hier::Box& ghost_box = state->getGhostBox();
   const hier::IntVector& ghost_cells = state->getGhostCellWidth();

   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_width_to_fill);

   int bc_types[3] = { Bdry::FACE3D, Bdry::EDGE3D, Bdry::NODE3D };

   int imin = ghost_box.lower(0);
   int imax = ghost_box.upper(0);
   int jmin = ghost_box.lower(1);
   int jmax = ghost_box.upper(1);
   int kmin = ghost_box.lower(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nxny = nx * ny;

   const std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());

   for (int ii = 0; ii < 3; ++ii) {

      const std::vector<hier::BoundaryBox>& bc_bdry =
         pgeom->getCodimensionBoundaries(bc_types[ii]);

      for (int jj = 0; jj < static_cast<int>(bc_bdry.size()); ++jj) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(bc_bdry[jj],
               interior,
               gcw_to_fill);
         int l_imin = fill_box.lower(0);
         int l_jmin = fill_box.lower(1);
         int l_kmin = fill_box.lower(2);

         int l_imax = fill_box.upper(0);
         int l_jmax = fill_box.upper(1);
         int l_kmax = fill_box.upper(2);

         int nd = state->getDepth();

         for (int n = 0; n < nd; ++n) {
            double* sptr = state->getPointer(n);

            for (int k = l_kmin; k <= l_kmax; ++k) {
               for (int j = l_jmin; j <= l_jmax; ++j) {
                  for (int i = l_imin; i <= l_imax; ++i) {

                     int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

                     sptr[ind] = ii;

                  }
               }
            }
         } // end of marking loops

      } // loop over boxes in each type

   } // loop over box types

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

extern "C" {
void SAMRAI_F77_FUNC(bcmultiblock, BCMULTIBLOCK) (
   const int *, const int *, const int *,
   const int *, const int *, const int *,
   double *,
   const int *, const int *, const int *,
   const int *, const int *, const int *,
   double *,
   const int* ncomp,
   const int* dlo, const int* dhi,
   const int* glo, const int* ghi,
   const int* lo, const int* hi,
   const int* dir, const int* side);
}

void MblkEuler::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);

   const hier::BlockId::block_t block_number =
      patch.getBox().getBlockId().getBlockValue();

   markPhysicalBoundaryConditions(patch, ghost_width_to_fill);

   std::shared_ptr<pdat::NodeData<double> > position(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > state(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_state, getDataContext())));
   TBOX_ASSERT(position);
   TBOX_ASSERT(state);

   //
   // the patch and its ghost box
   //
   const hier::Box& patch_box = patch.getBox();
   const hier::Box& ghost_box = state->getGhostBox();
   const hier::Box& nghost_box = position->getGhostBox();

   double* state_ptr = state->getPointer();
   int state_depth = state->getDepth();
   double* position_ptr = position->getPointer();

   //
   // the index space of this block and its neighbors
   //
   const std::shared_ptr<hier::PatchGeometry> patch_geom(
      patch.getPatchGeometry());
   const hier::IntVector ratio = patch_geom->getRatio();
   hier::BoxContainer domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, ratio,
      patch.getBox().getBlockId());

   const hier::IntVector& periodic =
      d_grid_geometry->getPeriodicShift(
         hier::IntVector(d_dim, 1));

   const hier::Box& domain_box = domain_boxes.front();
   // domain_box.refine(patch_geom->getRatio());

   d_mblk_geometry->buildLocalBlocks(patch_box,
      domain_box,
      patch.getBox().getBlockId().getBlockValue(),
      d_dom_local_blocks);
   //
   // loop over the directions, filling in boundary conditions where needed
   // note that here we have a check for is this really a physical boundary or is it
   // just a periodic boundary condition, or is it just an internal block boundary
   //
   int imin = ghost_box.lower(0);
   int imax = ghost_box.upper(0);
   int jmin = ghost_box.lower(1);
   int jmax = ghost_box.upper(1);
   int kmin = ghost_box.lower(2);
   int kmax = ghost_box.upper(2);

   int nd_imin = nghost_box.lower(0);
   int nd_imax = nghost_box.upper(0);
   int nd_jmin = nghost_box.lower(1);
   int nd_jmax = nghost_box.upper(1);
   int nd_kmin = nghost_box.lower(2);
   int nd_kmax = nghost_box.upper(2);

   for (tbox::Dimension::dir_t dir = 0; dir < d_dim.getValue(); ++dir) {
      if (!periodic(dir)) {

         if ((ghost_box.lower(dir) < domain_box.lower(dir)) &&
             (d_dom_local_blocks[dir] == block_number)) {
            int iside = 0;
            int intdir = dir;
            SAMRAI_F77_FUNC(bcmultiblock, BCMULTIBLOCK) (
               &nd_imin, &nd_imax, &nd_jmin, &nd_jmax, &nd_kmin, &nd_kmax,
               position_ptr,
               &imin, &imax, &jmin, &jmax, &kmin, &kmax, state_ptr,
               &state_depth,
               &domain_box.lower()[0],
               &domain_box.upper()[0],
               &ghost_box.lower()[0],
               &ghost_box.upper()[0],
               &patch_box.lower()[0],
               &patch_box.upper()[0],
               &intdir, &iside);
         }

         if ((ghost_box.upper(dir) > domain_box.upper(dir)) &&
             (d_dom_local_blocks[d_dim.getValue() + dir] == block_number)) {
            int iside = 1;
            int intdir = dir;
            SAMRAI_F77_FUNC(bcmultiblock, BCMULTIBLOCK) (
               &nd_imin, &nd_imax, &nd_jmin, &nd_jmax, &nd_kmin, &nd_kmax,
               position_ptr,
               &imin, &imax, &jmin, &jmax, &kmin, &kmax, state_ptr,
               &state_depth,
               &domain_box.lower()[0],
               &domain_box.upper()[0],
               &ghost_box.lower()[0],
               &ghost_box.upper()[0],
               &patch_box.lower()[0],
               &patch_box.upper()[0],
               &intdir, &iside);
         }

      } // end of periodic check

   } // end of direction loop

}

/*
 *************************************************************************
 *
 * Refine operations
 *
 *************************************************************************
 */

void MblkEuler::preprocessRefine(
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
   setMappedGridOnPatch(coarse);
   setMappedGridOnPatch(fine);
   setVolumeOnPatch(coarse);
   setVolumeOnPatch(fine);
}

// ---------------------------------------------------------------------------

//
// the refinement operator
//
void MblkEuler::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,                                 // where the fine data is needed
   const hier::IntVector& ratio)
{
   std::shared_ptr<pdat::CellData<double> > cstate(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > fstate(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_vol, getDataContext())));

   TBOX_ASSERT(cstate);
   TBOX_ASSERT(fstate);
   TBOX_ASSERT(cvol);
   TBOX_ASSERT(fvol);
   TBOX_ASSERT(cstate->getDepth() == fstate->getDepth());

   int depth = cstate->getDepth();

   //
   // get the boxes and bounds
   //

   // ... the bounds of the data
   const hier::Box cgbox(cstate->getGhostBox());
   const hier::Box fgbox(fstate->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fgbox.lower();
   const hier::Index fihi = fgbox.upper();

   // ... the bounds we actually work on
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
   tbox::plog << "clevel     = " << clev << std::endl << std::endl;
   tbox::plog << "fine_box         = " << fine_box << std::endl;
   tbox::plog << "fine_patch_box   = " << fine.getBox() << std::endl;
   tbox::plog << "fine_ghost_box   = " << fgbox << std::endl << std::endl;

   tbox::plog << "coarse_box       = " << coarse_box << std::endl;
   tbox::plog << "coarse_patch_box = " << coarse.getBox() << std::endl;
   tbox::plog << "coarse_ghost_box = " << coarse_box << std::endl;

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

   //
   // setup coarse strides
   //
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
   int cnx = cimax - cimin + 1;
   int cny = cjmax - cjmin + 1;
   int cnxny = cnx * cny;

   int fimin = filo(0);  // the fine data bounds
   int fimax = fihi(0);
   int fjmin = filo(1);
   int fjmax = fihi(1);
   int fkmin = filo(2);
   int fnx = fimax - fimin + 1;
   int fny = fjmax - fjmin + 1;
   int fnxny = fnx * fny;

   double rat0 = ratio[0];
   double rat1 = ratio[1];
   double rat2 = ratio[2];
   double fact = 2.0;  // xi varies from -1 to 1

   //
   // ================================= state variable refinement ====================
   //
   double* cdata = 0; // keeps pointers around till end of loop
   double* fdata = 0;

   for (int n = 0; n < depth; ++n) {

      cdata = cstate->getPointer(n);
      fdata = fstate->getPointer(n);

      for (int l = 0; l < nel; ++l) {       // default slopes are zero
         slope0[l] = 0.0;  // this yields piecewise constant interpolation
         slope1[l] = 0.0;  // and makes a handy initializer
         slope2[l] = 0.0;
      }

      for (int k = ifirstc(2); k <= ilastc(2); ++k) {
         for (int j = ifirstc(1); j <= ilastc(1); ++j) {
            for (int i = ifirstc(0); i <= ilastc(0); ++i) {

               int ind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);
               int cind = POLY3(i, j, k, cimin, cjmin, ckmin, cnx, cnxny);

               double aii = cdata[cind];
               val[ind] = aii;

               TBOX_ASSERT(ind >= 0);   // debug assertions
               TBOX_ASSERT(ind < nel);
#if 0 // turn to zero for simple interp
               int im1 = cind - 1;
               int ip1 = cind + 1;
               int jm1 = cind - cnx;
               int jp1 = cind + cnx;
               int km1 = cind - cnxny;
               int kp1 = cind + cnxny;

               double w_i = cvolume[cind];
               double w_ip = cvolume[ip1];
               double w_im = cvolume[im1];
               double w_jp = cvolume[jp1];
               double w_jm = cvolume[jm1];
               double w_kp = cvolume[kp1];
               double w_km = cvolume[km1];

               double aip = cdata[ip1];
               double aim = cdata[im1];
               double ajp = cdata[jp1];
               double ajm = cdata[jm1];
               double akp = cdata[kp1];
               double akm = cdata[km1];

               my_slopes(aii, aip, aim, ajp, ajm, akp, akm,
                  w_i, w_ip, w_im, w_jp, w_jm, w_kp, w_km,
                  slope0[ind],
                  slope1[ind],
                  slope2[ind]);
#endif
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

               int find = POLY3(i, j, k, fimin, fjmin, fkmin, fnx, fnxny);

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

               fdata[find] =
                  val[ind]
                  + ldx * slope0[ind]
                  + ldy * slope1[ind]
                  + ldz * slope2[ind];
            }
         }
      } // end of i,j,k loops for finding the fine state variables

   } // end of state loop

   tbox::plog << "--------------------- end postprocessRefine" << std::endl;
}

/*
 *************************************************************************
 *
 * Coarsen operations
 *
 *************************************************************************
 */
void MblkEuler::preprocessCoarsen(
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
   setMappedGridOnPatch(coarse);
   setMappedGridOnPatch(fine);
   setVolumeOnPatch(coarse);
   setVolumeOnPatch(fine);
}

// ---------------------------------------------------------------------------

//
// the coarsening function
//
void MblkEuler::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   std::shared_ptr<pdat::CellData<double> > cstate(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > cvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::CellData<double> > fstate(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_state, getDataContext())));
   std::shared_ptr<pdat::CellData<double> > fvol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(d_vol, getDataContext())));

   TBOX_ASSERT(cstate);
   TBOX_ASSERT(cvol);
   TBOX_ASSERT(fstate);
   TBOX_ASSERT(fvol);
   TBOX_ASSERT(cstate->getDepth() == fstate->getDepth());

   int depth = cstate->getDepth();

   //
   // box and geometry information
   //
   const hier::Index filo = fstate->getGhostBox().lower();
   const hier::Index fihi = fstate->getGhostBox().upper();
   const hier::Index cilo = cstate->getGhostBox().lower();
   const hier::Index cihi = cstate->getGhostBox().upper();

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
   double* cdata = cstate->getPointer();
   double* fdata = fstate->getPointer();

   //
   // average the data
   //
   for (int n = 0; n < depth; ++n) {

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

               double ric = (double(i) + 0.5) / rat0 - 0.5;
               double rjc = (double(j) + 0.5) / rat1 - 0.5;
               double rkc = (double(k) + 0.5) / rat2 - 0.5;

               int ic = (int)(ric + (ric >= 0 ? 0.5 : -0.5));        // a round operation
               int jc = (int)(rjc + (rjc >= 0 ? 0.5 : -0.5));        // shift up and truncate if ic > 0
               int kc = (int)(rkc + (rkc >= 0 ? 0.5 : -0.5));        // shift down and truncate if ic < 0

               int chind =
                  POLY3(ic, jc, kc, cimin, cjmin, ckmin, cnx, cnxny) + n * cnel;  // post + state offset

               TBOX_ASSERT(cimin <= ic && ic <= cimax);
               TBOX_ASSERT(cjmin <= jc && jc <= cjmax);
               TBOX_ASSERT(ckmin <= kc && kc <= ckmax);

               double fmass = fvolume[vol_ind];

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

               double cmass = cvolume[vol_ind];

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
void MblkEuler::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_indx,
   const bool uses_richardson_extrapolation_too)
{

   NULL_USE(initial_error);
   NULL_USE(regrid_time);
   NULL_USE(uses_richardson_extrapolation_too);

   int level_number = patch.getPatchLevelNumber();
   setMappedGridOnPatch(patch);

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
         patch.getPatchData(d_state, getDataContext())));
   TBOX_ASSERT(tags);
   TBOX_ASSERT(var);

   //
   // Create a set of temporary tags and set to untagged value.
   //
   std::shared_ptr<pdat::CellData<int> > temp_tags(
      new pdat::CellData<int>(
         pbox,
         1,
         d_nghosts));
   temp_tags->fillAll(FALSE);

   hier::IntVector tag_ghost = tags->getGhostCellWidth();

   hier::IntVector nghost_cells = xyz->getGhostCellWidth();
   int nd_imin = ifirst(0) - nghost_cells(0);
   int nd_imax = ilast(0) + 1 + nghost_cells(0);
   int nd_jmin = ifirst(1) - nghost_cells(1);
   int nd_jmax = ilast(1) + 1 + nghost_cells(1);
   int nd_kmin = ifirst(2) - nghost_cells(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   hier::IntVector v_ghost = var->getGhostCellWidth();     // has ghost zones
   int imin = ifirst(0) - v_ghost(0); // the polynomial for the field
   int imax = ilast(0) + v_ghost(0);
   int jmin = ifirst(1) - v_ghost(1);
   int jmax = ilast(1) + v_ghost(1);
   int kmin = ifirst(2) - v_ghost(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nxny = nx * ny;

   hier::IntVector temp_tag_ghost = temp_tags->getGhostCellWidth();
   int imn = ifirst(0) - temp_tag_ghost(0);  // the polynomial for temp_tags
   int imx = ilast(0) + temp_tag_ghost(0);
   int jmn = ifirst(1) - temp_tag_ghost(1);
   int jmx = ilast(1) + temp_tag_ghost(1);
   int kmn = ifirst(2) - temp_tag_ghost(2);
   int tnx = imx - imn + 1;
   int tny = jmx - jmn + 1;
   int tnxny = tnx * tny;

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

      TBOX_ASSERT(var);

      std::string ref = d_refinement_criteria[ncrit];

      if (ref == "GRADIENT") {
         int nStateLocal = static_cast<int>(d_state_grad_names.size());
         for (int id = 0; id < nStateLocal; ++id) {

            double* lvar = var->getPointer(d_state_grad_id[id]);

            int size = static_cast<int>(d_state_grad_tol[id].size());  // max depth of gradient tolerance
            double tol = ((error_level_number < size)    // find the tolerance
                          ? d_state_grad_tol[id][error_level_number]
                          : d_state_grad_tol[id][size - 1]);

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
                     int n6 = n1 + nd_nxny + 1;                              //  1, -1,  1
                     int n7 = n1 + nd_nxny + 1 + nd_nx;                      //  1,  1,  1
                     int n8 = n1 + nd_nxny + nd_nx;                                  // -1,  1,  1

                     // ------------------------------------------------ x

                     double x1 = 0.25 * (x[n1] + x[n4] + x[n5] + x[n8]);  // xi
                     double x2 = 0.25 * (x[n2] + x[n3] + x[n6] + x[n7]);

                     double x3 = 0.25 * (x[n1] + x[n2] + x[n5] + x[n6]);  // eta
                     double x4 = 0.25 * (x[n3] + x[n4] + x[n7] + x[n8]);

                     double x5 = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4]);  // zeta
                     double x6 = 0.25 * (x[n5] + x[n6] + x[n7] + x[n8]);

                     // ------------------------------------------------ y

                     double y1 = 0.25 * (y[n1] + y[n4] + y[n5] + y[n8]);
                     double y2 = 0.25 * (y[n2] + y[n3] + y[n6] + y[n7]);

                     double y3 = 0.25 * (y[n1] + y[n2] + y[n5] + y[n6]);
                     double y4 = 0.25 * (y[n3] + y[n4] + y[n7] + y[n8]);

                     double y5 = 0.25 * (y[n1] + y[n2] + y[n3] + y[n4]);
                     double y6 = 0.25 * (y[n5] + y[n6] + y[n7] + y[n8]);

                     // ------------------------------------------------ z

                     double z1 = 0.25 * (z[n1] + z[n4] + z[n5] + z[n8]);
                     double z2 = 0.25 * (z[n2] + z[n3] + z[n6] + z[n7]);

                     double z3 = 0.25 * (z[n1] + z[n2] + z[n5] + z[n6]);
                     double z4 = 0.25 * (z[n3] + z[n4] + z[n7] + z[n8]);

                     double z5 = 0.25 * (z[n1] + z[n2] + z[n3] + z[n4]);
                     double z6 = 0.25 * (z[n5] + z[n6] + z[n7] + z[n8]);

                     //
                     // the components of the matrices that we want to invert
                     //
                     double dx_xi = 0.5 * (x2 - x1);
                     double dy_xi = 0.5 * (y2 - y1);
                     double dz_xi = 0.5 * (z2 - z1);

                     double dx_eta = 0.5 * (x4 - x3);
                     double dy_eta = 0.5 * (y4 - y3);
                     double dz_eta = 0.5 * (z4 - z3);

                     double dx_zeta = 0.5 * (x6 - x5);
                     double dy_zeta = 0.5 * (y6 - y5);
                     double dz_zeta = 0.5 * (z6 - z5);

                     //
                     // invert dx/dxi as in dx/dxi d/dx = d/dxi, note this
                     // is the transpose of the above matrix M, also via
                     // Kramer's rule
                     //
                     double detMt = (dx_xi * dy_eta * dz_zeta
                                     + dx_eta * dy_zeta * dz_xi
                                     + dx_zeta * dy_xi * dz_eta
                                     - dx_zeta * dy_eta * dz_xi
                                     - dx_eta * dy_xi * dz_zeta
                                     - dx_xi * dy_zeta * dz_eta);

                     double detC11 = dy_eta * dz_zeta - dz_eta * dy_zeta;
                     double detC21 = dx_eta * dz_zeta - dz_eta * dx_zeta;
                     double detC31 = dx_eta * dy_zeta - dy_eta * dx_zeta;

                     double detC12 = dy_xi * dz_zeta - dz_xi * dy_zeta;
                     double detC22 = dx_xi * dz_zeta - dz_xi * dx_zeta;
                     double detC32 = dx_xi * dy_zeta - dy_xi * dx_zeta;

                     double detC13 = dy_xi * dz_eta - dz_xi * dy_eta;
                     double detC23 = dx_xi * dz_eta - dz_xi * dx_eta;
                     double detC33 = dx_xi * dy_eta - dy_xi + dx_eta;

                     // -------------------

                     double b11 = detC11 / detMt;
                     double b21 = detC21 / detMt;
                     double b31 = detC31 / detMt;

                     double b12 = detC12 / detMt;
                     double b22 = detC22 / detMt;
                     double b32 = detC32 / detMt;

                     double b13 = detC13 / detMt;
                     double b23 = detC23 / detMt;
                     double b33 = detC33 / detMt;

                     //
                     // determine the maximum gradient in x, y and z (nice orthonormal basis)
                     //
                     dv_x[0] = b11 * dv_xi[0] + b12 * dv_xi[1] + b13 * dv_xi[2];
                     dv_x[1] = b21 * dv_xi[0] + b22 * dv_xi[1] + b23 * dv_xi[2];
                     dv_x[2] = b31 * dv_xi[0] + b32 * dv_xi[1] + b33 * dv_xi[2];

                     double vmax = MAX(dv_x[0], MAX(dv_x[1], dv_x[2]));

                     if (vmax > tol)
                        ltags[tind] = TRUE;

                  } // i, j, k loops
               }
            }

         } // end of state position loop
      } // criteria = STATE_GRADIENT

      //
      // For user-defined fixed refinement, access refine box data from the MblkGeometry
      // class.
      //
      if (ref == "USER_DEFINED") {

         hier::BoxContainer refine_boxes;
         if (d_mblk_geometry->getRefineBoxes(refine_boxes,
                patch.getBox().getBlockId().getBlockValue(),
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
void MblkEuler::fillSingularityBoundaryConditions(
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
void MblkEuler::setMappedGridOnPatch(
   const hier::Patch& patch)
{

   //
   // compute level domain
   //
   const std::shared_ptr<hier::PatchGeometry> patch_geom(
      patch.getPatchGeometry());
   hier::IntVector ratio = patch_geom->getRatio();
   hier::BoxContainer domain_boxes;
   d_grid_geometry->computePhysicalDomain(domain_boxes, ratio,
      patch.getBox().getBlockId());

   //
   // statistics on the level domain
   //
   d_dom_current_nboxes = domain_boxes.size();

   hier::BoxContainer::iterator itr = domain_boxes.begin();
   d_dom_current_bounds[0] = itr->lower(0);
   d_dom_current_bounds[1] = itr->lower(1);
   d_dom_current_bounds[2] = itr->lower(2);
   d_dom_current_bounds[3] = itr->upper(0);
   d_dom_current_bounds[4] = itr->upper(1);
   d_dom_current_bounds[5] = itr->upper(2);
   ++itr;

   for (int i = 1; i < d_dom_current_nboxes; ++i, ++itr) {
      d_dom_current_bounds[0] = MIN(d_dom_current_bounds[0], itr->lower(0));
      d_dom_current_bounds[1] = MIN(d_dom_current_bounds[1], itr->lower(1));
      d_dom_current_bounds[2] = MIN(d_dom_current_bounds[2], itr->lower(2));

      d_dom_current_bounds[3] = MAX(d_dom_current_bounds[3], itr->upper(0));
      d_dom_current_bounds[4] = MAX(d_dom_current_bounds[4], itr->upper(1));
      d_dom_current_bounds[5] = MAX(d_dom_current_bounds[5], itr->upper(2));
   }

   //
   // now build the mesh
   //
   int xyz_id = hier::VariableDatabase::getDatabase()->
      mapVariableAndContextToIndex(d_xyz, getDataContext());

   d_mblk_geometry->buildGridOnPatch(patch,
      domain_boxes.front(),
      xyz_id,
      patch.getBox().getBlockId().getBlockValue(),
      d_dom_local_blocks);
}

/*
 *************************************************************************
 *                                                                       *
 * Private method to build the volume on a patch                         *
 *                                                                       *
 *************************************************************************
 */
void MblkEuler::setVolumeOnPatch(
   const hier::Patch& patch)
{
   std::shared_ptr<pdat::CellData<double> > vol(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_vol, getDataContext())));

   std::shared_ptr<pdat::NodeData<double> > xyz(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<double>, hier::PatchData>(
         patch.getPatchData(d_xyz, getDataContext())));
   TBOX_ASSERT(vol);
   TBOX_ASSERT(xyz);

   hier::IntVector vol_ghosts = vol->getGhostCellWidth();
   hier::IntVector xyz_ghosts = xyz->getGhostCellWidth();

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   int imin = ifirst(0) - vol_ghosts(0);
   int imax = ilast(0) + vol_ghosts(0);
   int jmin = ifirst(1) - vol_ghosts(1);
   int jmax = ilast(1) + vol_ghosts(1);
   int kmin = ifirst(2) - vol_ghosts(2);
   int nx = imax - imin + 1;
   int ny = jmax - jmin + 1;
   int nxny = nx * ny;

   int nd_imin = ifirst(0) - xyz_ghosts(0);
   int nd_imax = ilast(0) + 1 + xyz_ghosts(0);
   int nd_jmin = ifirst(1) - xyz_ghosts(1);
   int nd_jmax = ilast(1) + 1 + xyz_ghosts(1);
   int nd_kmin = ifirst(2) - xyz_ghosts(2);
   int nd_nx = nd_imax - nd_imin + 1;
   int nd_ny = nd_jmax - nd_jmin + 1;
   int nd_nxny = nd_nx * nd_ny;

   double* cvol = vol->getPointer();

   double* x = xyz->getPointer(0);
   double* y = xyz->getPointer(1);
   double* z = xyz->getPointer(2);

   for (int k = ifirst(2); k <= ilast(2); ++k) {
      for (int j = ifirst(1); j <= ilast(1); ++j) {
         for (int i = ifirst(0); i <= ilast(0); ++i) {

            int cind = POLY3(i, j, k, imin, jmin, kmin, nx, nxny);

            int n1 = POLY3(i, j, k, nd_imin, nd_jmin, nd_kmin, nd_nx, nd_nxny);
            int n2 = n1 + 1;
            int n3 = n1 + 1 + nd_nx;
            int n4 = n1 + nd_nx;

            int n5 = n1 + nd_nxny;
            int n6 = n1 + nd_nxny + 1;
            int n7 = n1 + nd_nxny + 1 + nd_nx;
            int n8 = n1 + nd_nxny + nd_nx;

            double lvol = UpwindVolume(x[n1], x[n2], x[n3], x[n4],
                  x[n5], x[n6], x[n7], x[n8],

                  y[n1], y[n2], y[n3], y[n4],
                  y[n5], y[n6], y[n7], y[n8],

                  z[n1], z[n2], z[n3], z[n4],
                  z[n5], z[n6], z[n7], z[n8]);

            cvol[cind] = lvol;
         }
      }
   }
}

// ================================= MblkEuler::Visualization and IO =============================

/*
 *************************************************************************
 *                                                                       *
 * Register VisIt data writer to write data to plot files that may       *
 * be postprocessed by the VisIt tool.                                   *
 *                                                                       *
 *************************************************************************
 */

#ifdef HAVE_HDF5
void MblkEuler::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);
   d_visit_writer = viz_writer;
}
#endif

/*
 *************************************************************************
 *                                                                       *
 * Write MblkEuler object state to specified stream.                        *
 *                                                                       *
 *************************************************************************
 */

void MblkEuler::printClassData(
   std::ostream& os) const
{
   os << "\nMblkEuler::printClassData..." << std::endl;
   os << "MblkEuler: this = " << (MblkEuler *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_grid_geometry = " << std::endl;
   for (hier::BlockId::block_t j = 0; j < d_grid_geometry->getNumberBlocks(); ++j) {
//      os << (geom::GridGeometry*)d_grid_geometry[j] << std::endl;
   }

   // ----------------------------------------------

   os << "Parameters for numerical method ..." << std::endl;
   os << "   d_advection_velocity = ";
   for (tbox::Dimension::dir_t j = 0; j < d_dim.getValue(); ++j) os << d_advection_velocity[j] << " ";
   os << std::endl;

   os << "   d_nghosts    = " << d_nghosts << std::endl;
   os << "   d_fluxghosts = " << d_fluxghosts << std::endl;

   os << "Problem description and initial data..." << std::endl;
   os << "   d_data_problem = " << d_data_problem << std::endl;
}

/*
 *************************************************************************
 *                                                                       *
 * Read data members from input.  All values set from restart can be     *
 * overridden by values in the input database.                           *
 *                                                                       *
 *************************************************************************
 */
void MblkEuler::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   std::shared_ptr<tbox::Database> db(input_db->getDatabase("MblkEuler"));

   //
   // --------------- load balancing inputs
   //
   // Note: if we are restarting, then we only allow nonuniform
   // workload to be used if nonuniform workload was used originally.
   //
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

   //
   // --------------- initialize boundary condition factors
   //
   if (db->keyExists("wall_factors")) {
      d_wall_factors = db->getIntegerVector("wall_factors");
   } else {
      d_wall_factors.resize(6);
      for (int i = 0; i < 6; ++i) d_wall_factors[i] = 0;
   }

   //
   // --------------- process the linear advection test ---------------------
   //
   d_advection_test = 1;
   d_advection_velocity[0] = d_advection_velocity[1] =
         d_advection_velocity[2] = FLT_MAX;
   d_advection_vel_type = 0;
   if (db->keyExists("advection_velocity")) {
      db->getDoubleArray("advection_velocity",
         d_advection_velocity, d_dim.getValue());
      d_advection_vel_type = db->getInteger("advection_vel_type");
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `advection_velocity' not found in input.");
   }

   //
   // --------------- The state names inputs ---------------------
   //
   if (d_advection_test) {
      if (db->keyExists("state_names")) {
         d_state_names = db->getStringVector("state_names");
         d_nState = static_cast<int>(d_state_names.size());
      } else {
         TBOX_ERROR("missing 'state_names' input for sizing the state" << std::endl);
      }
   } else {
      TBOX_ASSERT(d_advection_test);
   }

   //
   // --------------- region Initialization inputs ---------------------
   //
   if (!is_from_restart) {

      if (db->keyExists("data_problem")) {
         d_data_problem = db->getString("data_problem");
      } else {
         TBOX_ERROR(
            d_object_name << ": "
                          << "`data_problem' value not found in input."
                          << std::endl);
      }

      int problem_1d = 1;

      //
      //  axis of revolution inputs
      //
      if (d_data_problem == "REVOLUTION") {

         if (db->keyExists("center")) {
            db->getDoubleArray("center", d_center, d_dim.getValue());
         } else {
            TBOX_ERROR(
               "`center' input required for REVOLUTION problem." << std::endl);
         }

         if (db->keyExists("axis")) {
            db->getDoubleArray("axis", d_axis, d_dim.getValue());
         } else {
            TBOX_ERROR("`axis' input required for REVOLUTION problem." << std::endl);
         }

         // normalize the axis to a unit vector
         double anorm = sqrt(
               d_axis[0] * d_axis[0] + d_axis[1] * d_axis[1] + d_axis[2]
               * d_axis[2]);
         d_axis[0] /= anorm;
         d_axis[1] /= anorm;
         d_axis[2] /= anorm;

         d_rev_rad.resize(d_number_of_regions);
         d_rev_axis.resize(d_number_of_regions);

         for (int i = 0; i < d_number_of_regions; ++i) {

            char tmp[20];
            snprintf(tmp, 20, "region_%d", i + 1);  //
            std::string lkey = tmp;
            std::shared_ptr<tbox::Database> region_db(
               db->getDatabase(lkey));

            d_rev_rad[i] = region_db->getDoubleVector("radius");
            d_rev_axis[i] = region_db->getDoubleVector("axis");

         }

         problem_1d = 0;
      }

      //
      //  Spherical inputs
      //
      if (d_data_problem == "SPHERE") {

         if (db->keyExists("center")) {
            db->getDoubleArray("center", d_center, d_dim.getValue());
         } else {
            TBOX_ERROR("`center' input required for SPHERE problem." << std::endl);
         }
      }

      //
      //  Rayleigh tayler shock tube inputs
      //
      if (d_data_problem == "RT_SHOCK_TUBE") {

         if (db->keyExists("amn")) {
            d_dt_ampl = db->getDouble("ampl");
            d_amn = db->getDoubleVector("amn");
            d_m_mode = db->getDoubleVector("m_mode");
            d_n_mode = db->getDoubleVector("n_mode");
            d_phiy = db->getDoubleVector("phiy");
            d_phiz = db->getDoubleVector("phiz");
         } else {
            TBOX_ERROR("missing input for RT_SHOCK_TUBE problem." << std::endl);
         }
      }

      //
      //  shared inputs for the 1d style problems
      //
      if (problem_1d) {
         if (db->keyExists("front_position")) {
            d_front_position = db->getDoubleVector("front_position");
            d_number_of_regions = static_cast<int>(d_front_position.size()) - 1;
            TBOX_ASSERT(d_number_of_regions > 0);
         } else {
            TBOX_ERROR("Missing`front_position' input required" << std::endl);
         }
      }

      //
      // the state data entries for the initial conditions
      //
//      int llen = d_number_of_regions*d_nState;
//     //double *tmp = new double[llen];
//    std::vector<double> tmp(llen);
//   for ( int ii = 0 ; ii < llen ; ++ii ) {
//       tmp[ii] = FLT_MAX;
//     }

      d_state_ic.resize(d_number_of_regions);
      //d_state_ic  = new double *[d_number_of_regions];
      for (int iReg = 0; iReg < d_number_of_regions; ++iReg) {
         d_state_ic[iReg].resize(d_nState);
         //d_state_ic[iReg] = tmp[d_nState*iReg];
      }

      //
      // pull in the data for each region
      //
      if (d_advection_test) {
         if (db->keyExists("state_data")) {
            std::shared_ptr<tbox::Database> state_db(
               db->getDatabase("state_data"));
            std::vector<double> lpsi;
            for (int iState = 0; iState < d_nState; ++iState) {
               lpsi = state_db->getDoubleVector(d_state_names[iState]);
               for (int iReg = 0; iReg < d_number_of_regions; ++iReg) {
                  d_state_ic[iReg][iState] = lpsi[iReg];
               }
            }
         } else {
            TBOX_ERROR(
               "missing 'state_data' input for initial conditions" << std::endl);
         }
      } else {
         TBOX_ASSERT(d_advection_test);
      }

   } // if !is_from_restart read in problem data

   //
   //  --------------- refinement criteria inputs
   //
   if (db->keyExists("Refinement_data")) {
      std::shared_ptr<tbox::Database> refine_db = db->getDatabase(
            "Refinement_data");
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

            //
            // allow only valid refinement criteria as the remaining keys
            //
            if (!(error_key == "GRADIENT")) {

               TBOX_ERROR(
                  "Unknown refinement criteria: " << error_key
                                                  << " in input."
                                                  << std::endl);

            } else {
               error_db = refine_db->getDatabase(error_key);
               ref_keys_defined[def_key_cnt] = error_key;
               ++def_key_cnt;
            }

            //
            // process the specific keys
            //
            if (error_db && error_key == "GRADIENT") {

               d_state_grad_names = error_db->getStringVector("names");
               int nStateLocal = static_cast<int>(d_state_grad_names.size());

               d_state_grad_tol.resize(nStateLocal);
               d_state_grad_id.resize(nStateLocal);

               for (int id = 0; id < nStateLocal; ++id) {
                  std::string grad_name = d_state_grad_names[id];

                  // ... the index needed
                  d_state_grad_id[id] = -1;
                  bool found = false;
                  for (int idj = 0; idj < d_nState && !found; ++idj) {
                     if (grad_name == d_state_names[idj]) {
                        found = true;
                        d_state_grad_id[id] = idj;
                     }
                  }

                  // ... the tolerance array needed
                  if (error_db->keyExists(grad_name)) {
                     d_state_grad_tol[id] =
                        error_db->getDoubleVector(grad_name);
                  } else {
                     TBOX_ERROR(
                        "No tolerance array " << grad_name
                                              << "found for gradient detector" << std::endl);
                  }

               }
            }

         } // refine criteria if test

      } // loop over refine criteria

   } // refine db entry exists

   //
   // --------------- boundary condition inputs ---------------------
   //
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
}

/*
 *************************************************************************
 *                                                                       *
 * Routines to put/get data members to/from restart database.            *
 *                                                                       *
 *************************************************************************
 */

void MblkEuler::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("MBLKEULER_VERSION", MBLKEULER_VERSION);

   restart_db->putDoubleArray("d_advection_velocity",
      d_advection_velocity,
      d_dim.getValue());

   restart_db->putIntegerArray("d_nghosts", &d_nghosts[0], d_dim.getValue());
   restart_db->putIntegerArray("d_fluxghosts",
      &d_fluxghosts[0],
      d_dim.getValue());

   restart_db->putString("d_data_problem", d_data_problem);
}

/*
 *************************************************************************
 *                                                                       *
 *    Access class information from restart database.                    *
 *                                                                       *
 *************************************************************************
 */
void MblkEuler::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file.");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("MBLKEULER_VERSION");
   if (ver != MBLKEULER_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   d_data_problem = db->getString("d_data_problem");
}

/*
 *************************************************************************
 *
 * Routine to check boundary data when debugging.
 *
 *************************************************************************
 */

void MblkEuler::checkBoundaryData(
   int btype,
   const hier::Patch& patch,
   const hier::IntVector& ghost_width_to_check,
   const std::vector<int>& scalar_bconds) const
{
   NULL_USE(btype);
   NULL_USE(patch);
   NULL_USE(ghost_width_to_check);
   NULL_USE(scalar_bconds);
}

hier::IntVector MblkEuler::getMultiblockRefineOpStencilWidth() const
{
   return hier::IntVector(d_dim, 1);
}

hier::IntVector MblkEuler::getMultiblockCoarsenOpStencilWidth()
{
   return hier::IntVector(d_dim, 0);
}
