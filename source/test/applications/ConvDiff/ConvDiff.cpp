/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for single patch in convection
 *                diffusion example.
 *
 ************************************************************************/

#include "ConvDiff.h"

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <float.h>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (1)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

// external definitions for Fortran numerical routines
#include "ConvDiffFort.h"

// defines for initialization
#define SPHERE (40)

// Number of ghosts cells used for each variable quantity.
#define CELLG (1)

// Define class version number
#define CONV_DIFF_VERSION (2)

/*
 *************************************************************************
 *
 * The constructor for ConvDiff class sets data members to defualt
 * values, creates variables that define the solution state for the
 * convection diffusion equation.
 *
 * After default values are set, this routine calls getFromRestart()
 * if execution from a restart file is specified.  Finally,
 * getFromInput() is called to read values from the given input
 * database (potentially overriding those found in the restart file).
 *
 *************************************************************************
 */

ConvDiff::ConvDiff(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db,
   std::shared_ptr<geom::CartesianGridGeometry> grid_geom):
   algs::MethodOfLinesPatchStrategy::MethodOfLinesPatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_grid_geometry(grid_geom),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_primitive_vars(new pdat::CellVariable<double>(dim, "primitive_vars", d_allocator)),
   d_function_eval(new pdat::CellVariable<double>(dim, "function_eval", d_allocator)),
   d_diffusion_coeff(1.),
   d_source_coeff(0.),
   d_cfl(0.9),
   d_nghosts(dim, 1),
   d_zero_ghosts(dim, 0),
   d_radius(tbox::MathUtilities<double>::getSignalingNaN())
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(grid_geom);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   int k;

   /*
    * *hier::Variable quantities that define state of convection diffusion
    * problem.
    *
    *     dU/dt + alpha dU/dx = mu d^2U/dx^2 + gamma
    *
    *     U     = primitive variable(s)
    *     F(U)  = function evaluation
    *     alpha = convection coefficient
    *     mu    = diffusion coefficient
    *     gamma = source coefficient
    */
   for (k = 0; k < d_dim.getValue(); ++k) d_convection_coeff[k] = 0.;

   // Physics parameters
   for (k = 0; k < NEQU; ++k) d_tolerance[k] = 0.;

   /*
    * Defaults for problem type and initial data.  Set initial
    * data to NaNs so we make sure input has set it to appropriate
    * problem.
    */
   d_data_problem = tbox::MathUtilities<char>::getMax(),

   // SPHERE problem...
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_center, d_dim.getValue());
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_val_inside, NEQU);
   tbox::MathUtilities<double>::setArrayToSignalingNaN(d_val_inside, NEQU);

   /*
    * Boundary condition initialization.
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

      d_bdry_edge_val.resize(NUM_2D_EDGES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_val);
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

      d_bdry_face_val.resize(NUM_3D_FACES);
      tbox::MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_val);
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

   if (d_data_problem == "SPHERE") {
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
               appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
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
               appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
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
               appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                  i, d_scalar_bdry_node_conds[i]);
         }
      }

   }

}

/*************************************************/
ConvDiff::~ConvDiff()
{
}

/*
 *************************************************************************
 *
 * Register variables with the MOL integrator.  Since the integrator
 * that manages storage for these quantities, this is how it finds
 * out about them.
 *
 *
 *************************************************************************
 */
void ConvDiff::registerModelVariables(
   algs::MethodOfLinesIntegrator* integrator)
{
   /*
    * Only two types of variables are used by the integrator - SOLN
    * and RHS.  The primitive variables are of type SOLN while
    * the function evaluation is of type RHS.
    *
    * SOLN needs two contexts - current (without ghosts) and
    * scratch (with ghosts). Current is maintained between timesteps,
    * while scratch is created and destroyed within a timestep.
    */

   integrator->registerVariable(d_primitive_vars, d_nghosts,
      algs::MethodOfLinesIntegrator::SOLN,
      d_grid_geometry,
      "CONSERVATIVE_COARSEN",
      "LINEAR_REFINE");

   /*
    * RHS needs only one context - scratch (with ghosts).  It is used
    * to store communicated ghost information within the timestep.
    * The function evaluation is not communicated across levels,
    * it is recomputed using the interpolated primitive variables.
    * Hence, we don't need to define a coarsen or refine operator.
    */
   integrator->registerVariable(d_function_eval, d_nghosts,
      algs::MethodOfLinesIntegrator::RHS,
      d_grid_geometry,
      "NO_COARSEN",
      "NO_REFINE");

   /*
    * Loop over primitive variables and register each with the
    * data writer.
    */
   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

   int prim_var_id = vardb->mapVariableAndContextToIndex(
         d_primitive_vars, getInteriorContext());

   std::string dump_name = "Primitive Var #";
   const int size = static_cast<int>(dump_name.length()) + 16;
   char* buffer = new char[size];

   for (int n = 0; n < NEQU; ++n) {
      snprintf(buffer, size, "%s%01d", dump_name.c_str(), n);
      std::string variable_name(buffer);
#ifdef HAVE_HDF5
      if (d_visit_writer) {
         d_visit_writer->
         registerPlotQuantity(variable_name, "SCALAR",
            prim_var_id, n);
      }
      if (!d_visit_writer) {
         TBOX_WARNING(
            d_object_name << ": registerModelVariables()\n"
                          << "VisIt data writer was not registered.\n"
                          << "Consequently, no plot data will\n"
                          << "be written." << std::endl);
      }
#endif

   }

   delete[] buffer;

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
void ConvDiff::initializeDataOnPatch(
   hier::Patch& patch,
   const double time,
   const bool initial_time) const
{
   NULL_USE(time);

   if (initial_time) {

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(patch_geom);

      const double* dx = patch_geom->getDx();
      const double* xlo = patch_geom->getXLower();
      const double* xhi = patch_geom->getXUpper();

      std::shared_ptr<pdat::CellData<double> > primitive_vars(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_primitive_vars, getInteriorContext())));

      TBOX_ASSERT(primitive_vars);

      hier::IntVector ghost_cells = primitive_vars->getGhostCellWidth();

      const hier::Box& pbox = patch.getBox();
      const hier::Index ifirst = pbox.lower();
      const hier::Index ilast = pbox.upper();

      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(initsphere2d, INITSPHERE2D) (dx, xlo, xhi,
            ifirst(0), ilast(0), ifirst(1), ilast(1),
            ghost_cells(0), ghost_cells(1),
            primitive_vars->getPointer(),
            d_val_inside,
            d_val_outside,
            d_center, d_radius,
            NEQU);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(initsphere3d, INITSPHERE3D) (dx, xlo, xhi,
            ifirst(0), ilast(0), ifirst(1), ilast(1),
            ifirst(2), ilast(2),
            ghost_cells(0), ghost_cells(1),
            ghost_cells(2),
            primitive_vars->getPointer(),
            d_val_inside,
            d_val_outside,
            d_center, d_radius,
            NEQU);
      }

      // tbox::plog << "Level:" << patch.getPatchLevelNumber() << "\n" << std::endl;
      // tbox::plog << "Patch:" << std::endl;
      // primitive_vars->print(pbox,tbox::plog);
   }
}

/*
 *************************************************************************
 *
 * Compute stable time increment for patch.  Return this value.
 * (See Hirsch, Vol 1, pp 448 for description of stability analysis)
 *
 *************************************************************************
 */
double ConvDiff::computeStableDtOnPatch(
   hier::Patch& patch,
   const double time) const
{
   NULL_USE(time);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   double stabdt = 0.0;

//  Use the condition defined on pg. 448 of Hirsch, Vol. 1.
//     for system du/dt = beta*u_xx
//           0 <= beta*dt/(dx**2) <= 1
//     assume Beta=d_cfl for this case.

   if (d_dim == tbox::Dimension(2)) {
      stabdt = d_cfl * ((*dx) * (*dx));
   } else if (d_dim == tbox::Dimension(3)) {
      stabdt = d_cfl * ((*dx) * (*dx) * (*dx));
   }

//   Alternatively, one could use a fortran function here if you want
//   something more complex.
//
//   FORT_STABLE_DT(stabdt
//                  ifirst(0),ilast(0),
//                  ifirst(1),ilast(1),
//                  ifirst(2),ilast(2),
//                  d_params,
//                  density->getPointer(),
//                  velocity->getPointer(),
//                  pressure->getPointer(),
//                  stabdt);
   return stabdt;

}

/*
 *************************************************************************
 *
 * Perform a single Runge-Kutta sub-iteration using the passed-in
 * alpha.
 *
 *************************************************************************
 */
void ConvDiff::singleStep(
   hier::Patch& patch,
   const double dt,
   const double alpha_1,
   const double alpha_2,
   const double beta) const
{

   std::shared_ptr<pdat::CellData<double> > prim_var_updated(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_primitive_vars, getInteriorWithGhostsContext())));

   std::shared_ptr<pdat::CellData<double> > prim_var_fixed(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_primitive_vars, getInteriorContext())));

   std::shared_ptr<pdat::CellData<double> > function_eval(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_function_eval, getInteriorContext())));
   TBOX_ASSERT(prim_var_updated);
   TBOX_ASSERT(prim_var_fixed);
   TBOX_ASSERT(function_eval);

   const hier::Box& pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

//      tbox::plog << "----primitive_var_current" << std::endl;
//      prim_var_current->print(prim_var_current->getGhostBox(),tbox::plog);
//      tbox::plog << "----primitive_var_scratch" << std::endl;
//      prim_var_scratch->print(prim_var_scratch->getGhostBox(),tbox::plog);
//
// Evaluate Right hand side F(prim_var_scratch)
//
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(computerhs2d, COMPUTERHS2D) (ifirst(0), ilast(0), ifirst(1),
         ilast(1),
         d_nghosts(0), d_nghosts(1),
         dx,
         d_convection_coeff,
         d_diffusion_coeff,
         d_source_coeff,
         prim_var_updated->getPointer(),
         function_eval->getPointer(),
         NEQU);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(computerhs3d, COMPUTERHS3D) (ifirst(0), ilast(0), ifirst(1),
         ilast(1),
         ifirst(2), ilast(2),
         d_nghosts(0), d_nghosts(1),
         d_nghosts(2),
         dx,
         d_convection_coeff,
         d_diffusion_coeff,
         d_source_coeff,
         prim_var_updated->getPointer(),
         function_eval->getPointer(),
         NEQU);
   }

//    tbox::plog << "Function Evaluation" << std::endl;
//    function_eval->print(function_eval->getBox());
//
// Take RK step
//
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(rkstep2d, RKSTEP2D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         d_nghosts(0), d_nghosts(1),
         dt, alpha_1, alpha_2, beta,
         d_convection_coeff,
         d_diffusion_coeff,
         d_source_coeff,
         prim_var_updated->getPointer(),
         prim_var_fixed->getPointer(),
         function_eval->getPointer(),
         NEQU);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(rkstep3d, RKSTEP3D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         d_nghosts(0), d_nghosts(1),
         d_nghosts(2),
         dt, alpha_1, alpha_2, beta,
         d_convection_coeff,
         d_diffusion_coeff,
         d_source_coeff,
         prim_var_updated->getPointer(),
         prim_var_fixed->getPointer(),
         function_eval->getPointer(),
         NEQU);
   }
//        tbox::plog << "----prim_var_scratch after RK step" << std::endl;
//        prim_var_scratch->print(prim_var_scratch->getGhostBox(),tbox::plog);

}

/*
 *************************************************************************
 *
 *  Cell tagging routine - tag cells that require refinement based on
 *  a provided condition.
 *
 *************************************************************************
 */
void ConvDiff::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_richardson_extrapolation_too);

   std::shared_ptr<pdat::CellData<int> > tags(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(tag_index)));
   std::shared_ptr<pdat::CellData<double> > primitive_vars(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_primitive_vars, getInteriorWithGhostsContext())));
   TBOX_ASSERT(tags);
   TBOX_ASSERT(primitive_vars);

   const hier::Box& pbox = patch.getBox();
   const hier::Index ifirst = pbox.lower();
   const hier::Index ilast = pbox.upper();

   const hier::IntVector var_ghosts = primitive_vars->getGhostCellWidth();

   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(tagcells2d, TAGCELLS2D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         var_ghosts(0), var_ghosts(1),
         tags->getPointer(),
         primitive_vars->getPointer(),
         true,
         d_tolerance,
         NEQU);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(tagcells3d, TAGCELLS3D) (ifirst(0), ilast(0), ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         var_ghosts(0), var_ghosts(1),
         var_ghosts(2),
         tags->getPointer(),
         primitive_vars->getPointer(),
         true,
         d_tolerance,
         NEQU);
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
void ConvDiff::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);

   std::shared_ptr<pdat::CellData<double> > primitive_vars(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(d_primitive_vars, getInteriorWithGhostsContext())));

   TBOX_ASSERT(primitive_vars);
   TBOX_ASSERT(primitive_vars->getGhostCellWidth() == d_nghosts);

   if (d_dim == tbox::Dimension(2)) {

      /*
       * Set boundary conditions for cells corresponding to patch edges.
       */
      appu::CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("primitive_vars", primitive_vars,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_edge_val);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      appu::CartesianBoundaryUtilities2::
      fillNodeBoundaryData("primitive_vars", primitive_vars,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_edge_val);

#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds);
#endif
#endif

   } // d_dim == tbox::Dimension(2))
   else if (d_dim == tbox::Dimension(3)) {

      /*
       *  Set boundary conditions for cells corresponding to patch faces.
       */
      appu::CartesianBoundaryUtilities3::
      fillFaceBoundaryData("primitive_vars", primitive_vars,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_face_conds,
         d_bdry_face_val);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_face_conds);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch edges.
       */

      appu::CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("primitive_vars", primitive_vars,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_face_val);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_edge_conds);
#endif
#endif

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      appu::CartesianBoundaryUtilities3::
      fillNodeBoundaryData("primitive_vars", primitive_vars,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_face_val);
#ifdef DEBUG_CHECK_ASSERTIONS
#if CHECK_BDRY_DATA
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill,
         d_scalar_bdry_node_conds);
#endif
#endif

   } // d_dim == tbox::Dimension(3))

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
void ConvDiff::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);
   d_visit_writer = viz_writer;
}
#endif

/*
 *************************************************************************
 *
 * Prints class data - writes out info in class if assertion is thrown
 *
 *************************************************************************
 */

void ConvDiff::printClassData(
   std::ostream& os) const
{
   fflush(stdout);
   int j;

   os << "ptr ConvDiff = " << (ConvDiff *)this << std::endl;
   os << "ptr grid geometry = "
      << d_grid_geometry.get() << std::endl;

   os << "Coefficients..." << std::endl;
   for (j = 0; j < d_dim.getValue(); ++j) os << "d_convection_coeff[" << j << "] = "
                                             << d_convection_coeff[j] << std::endl;
   os << "d_diffusion_coeff = " << d_diffusion_coeff << std::endl;
   os << "d_source_coeff = " << d_source_coeff << std::endl;

   os << "Problem description and initial data..." << std::endl;
   os << "   d_data_problem = " << d_data_problem << std::endl;

   os << "       d_radius = " << d_radius << std::endl;
   os << "       d_center = ";
   for (j = 0; j < d_dim.getValue(); ++j) os << d_center[j] << " ";
   os << std::endl;
   os << "       d_val_inside = ";
   for (j = 0; j < NEQU; ++j) os << d_val_inside[j] << " ";
   os << std::endl;
   os << "       d_val_outside = ";
   for (j = 0; j < NEQU; ++j) os << d_val_outside[j] << " ";
   os << std::endl;

   os << "Boundary Condition data..." << std::endl;
   if (d_dim == tbox::Dimension(2)) {
      for (j = 0; j < static_cast<int>(d_scalar_bdry_edge_conds.size()); ++j) {
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << std::endl;
         if (d_scalar_bdry_edge_conds[j] == BdryCond::DIRICHLET) {
            os << "         d_bdry_edge_val[" << j << "] = "
               << d_bdry_edge_val[j] << std::endl;
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
            os << "         d_bdry_face_val[" << j << "] = "
               << d_bdry_face_val[j] << std::endl;
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

   os << "         d_nghosts = " << d_nghosts << std::endl;
   os << "         d_zero_ghosts = " << d_zero_ghosts << std::endl;
   os << "         d_cfl = " << d_cfl << std::endl;

}

/*
 *************************************************************************
 *
 *
 *************************************************************************
 */
void ConvDiff::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   if (input_db->keyExists("convection_coeff")) {
      input_db->getDoubleArray("convection_coeff",
         d_convection_coeff, d_dim.getValue());
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `convection_coeff' not found in input.");
   }
   if (input_db->keyExists("diffusion_coeff")) {
      d_diffusion_coeff = input_db->getDouble("diffusion_coeff");
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `diffusion_coeff' not found in input.");
   }
   if (input_db->keyExists("source_coeff")) {
      d_source_coeff = input_db->getDouble("source_coeff");
   } else {
      TBOX_ERROR(d_object_name << ":  "
                               << "Key data `source_coeff' not found in input.");
   }

   d_cfl = input_db->getDoubleWithDefault("cfl", d_cfl);

   if (input_db->keyExists("cell_tagging_tolerance")) {
      input_db->getDoubleArray("cell_tagging_tolerance",
         d_tolerance, NEQU);
   } else {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Key data `cell_tagging_tolerance' not found in input.");
   }

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
         if (init_data_db->keyExists("val_inside")) {
            init_data_db->getDoubleArray("val_inside", d_val_inside, NEQU);
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "val_inside' input required for "
                                     << "SPHERE problem." << std::endl);
         }
         if (init_data_db->keyExists("val_outside")) {
            init_data_db->getDoubleArray("val_outside", d_val_outside, NEQU);
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "`val_outside' input required for "
                                     << "SPHERE problem." << std::endl);
         }

         found_problem_data = true;

      }

      if (!found_problem_data) {
         TBOX_ERROR(
            d_object_name << ": "
                          << "Bad data in `Initial_data' database."
                          << std::endl);
      }

      const hier::IntVector& one_vec(hier::IntVector::getOne(d_dim));
      hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);
      int num_per_dirs = 0;
      for (int id = 0; id < d_dim.getValue(); ++id) {
         if (periodic(id)) ++num_per_dirs;
      }

      if (input_db->keyExists("Boundary_data")) {
         std::shared_ptr<tbox::Database> boundary_db(
            input_db->getDatabase("Boundary_data"));

         if (d_dim == tbox::Dimension(2)) {
            appu::CartesianBoundaryUtilities2::getFromInput(this,
               boundary_db,
               d_scalar_bdry_edge_conds,
               d_scalar_bdry_node_conds,
               periodic);
         } else if (d_dim == tbox::Dimension(3)) {
            appu::CartesianBoundaryUtilities3::getFromInput(this,
               boundary_db,
               d_scalar_bdry_face_conds,
               d_scalar_bdry_edge_conds,
               d_scalar_bdry_node_conds,
               periodic);
         }

      } else {
         TBOX_WARNING(
            d_object_name << ": "
                          << "Key data `Boundary_data' not found in input. "
                          << "Using default FLOW boundary conditions." << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 * Routines to put/get data members to/from restart database.
 *
 *************************************************************************
 */

void ConvDiff::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("CONV_DIFF_VERSION", CONV_DIFF_VERSION);

   restart_db->putDouble("d_diffusion_coeff", d_diffusion_coeff);
   restart_db->putDoubleArray("d_convection_coeff",
      d_convection_coeff,
      d_dim.getValue());
   restart_db->putDouble("d_source_coeff", d_source_coeff);
   restart_db->putIntegerArray("d_nghosts", &d_nghosts[0], d_dim.getValue());

   restart_db->putString("d_data_problem", d_data_problem);

   restart_db->putDouble("d_radius", d_radius);
   restart_db->putDoubleArray("d_center", d_center, d_dim.getValue());
   restart_db->putDoubleArray("d_val_inside", d_val_inside, NEQU);
   restart_db->putDoubleArray("d_val_outside", d_val_outside, NEQU);

   restart_db->putDouble("d_cfl", d_cfl);

   restart_db->putIntegerVector("d_scalar_bdry_edge_conds",
      d_scalar_bdry_edge_conds);
   restart_db->putIntegerVector("d_scalar_bdry_node_conds",
      d_scalar_bdry_node_conds);

   if (d_dim == tbox::Dimension(2)) {
      restart_db->putDoubleVector("d_bdry_edge_val", d_bdry_edge_val);
   } else if (d_dim == tbox::Dimension(3)) {
      restart_db->putIntegerVector("d_scalar_bdry_face_conds",
         d_scalar_bdry_face_conds);
      restart_db->putDoubleVector("d_bdry_face_val", d_bdry_face_val);
   }

}

/*
 *************************************************************************
 *
 *
 *************************************************************************
 */
void ConvDiff::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in the restart file.");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("CONV_DIFF_VERSION");
   if (ver != CONV_DIFF_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   d_diffusion_coeff = db->getDouble("d_diffusion_coeff");
   db->getDoubleArray("d_convection_coeff", d_convection_coeff, d_dim.getValue());
   d_source_coeff = db->getDouble("d_source_coeff");
   db->getIntegerArray("d_nghosts", &d_nghosts[0], d_dim.getValue());

   d_data_problem = static_cast<char>(db->getInteger("d_data_problem"));

   d_radius = db->getDouble("d_radius");
   db->getDoubleArray("d_center", d_center, d_dim.getValue());
   db->getDoubleArray("d_val_inside", d_val_inside, NEQU);
   db->getDoubleArray("d_val_outside", d_val_outside, NEQU);

   d_cfl = db->getDouble("d_cfl");

   d_scalar_bdry_edge_conds = db->getIntegerVector("d_scalar_bdry_edge_conds");
   d_scalar_bdry_node_conds = db->getIntegerVector("d_scalar_bdry_node_conds");

   if (d_dim == tbox::Dimension(2)) {
      d_bdry_edge_val = db->getDoubleVector("d_bdry_edge_val");
   } else if (d_dim == tbox::Dimension(3)) {
      d_scalar_bdry_face_conds =
         db->getIntegerVector("d_scalar_bdry_face_conds");

      d_bdry_face_val = db->getDoubleVector("d_bdry_face_val");
   }

}

/*
 *************************************************************************
 *
 * Routines to read boundary data from input database.
 *
 *************************************************************************
 */

void ConvDiff::readDirichletBoundaryDataEntry(
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
         d_bdry_edge_val);
   } else if (d_dim == tbox::Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_val);
   }
}

void ConvDiff::readNeumannBoundaryDataEntry(
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
         d_bdry_edge_val);
   } else if (d_dim == tbox::Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_val);
   }
}

void ConvDiff::readStateDataEntry(
   std::shared_ptr<tbox::Database> db,
   const std::string& db_name,
   int array_indx,
   std::vector<double>& val)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(array_indx >= 0);
   TBOX_ASSERT(static_cast<int>(val.size()) > array_indx);

   if (db->keyExists("val")) {
      val[array_indx] = db->getDouble("val");
   } else {
      TBOX_ERROR(d_object_name << ": "
                               << "`val' entry missing from " << db_name
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

void ConvDiff::checkBoundaryData(
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

      int bscalarcase = 0, refbdryloc = 0;
      if (d_dim == tbox::Dimension(2)) {
         if (btype == Bdry::EDGE2D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_2D_EDGES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = bloc;
         } else { // btype == Bdry::NODE2D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_2D_NODES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_node_bdry_edge[bloc];
         }
      } else if (d_dim == tbox::Dimension(3)) {
         if (btype == Bdry::FACE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_FACES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = bloc;
         } else if (btype == Bdry::EDGE3D) {
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_EDGES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_edge_bdry_face[bloc];
         } else { // btype == Bdry::NODE3D
            TBOX_ASSERT(static_cast<int>(scalar_bconds.size()) ==
               NUM_3D_NODES);
            bscalarcase = scalar_bconds[bloc];
            refbdryloc = d_node_bdry_face[bloc];
         }
      }

      int num_bad_values = 0;

      if (d_dim == tbox::Dimension(2)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities2::checkBdryData(
               d_primitive_vars->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_primitive_vars,
                  getInteriorWithGhostsContext()),
               0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_edge_val[refbdryloc]);
      } else if (d_dim == tbox::Dimension(3)) {
         num_bad_values =
            appu::CartesianBoundaryUtilities3::checkBdryData(
               d_primitive_vars->getName(),
               patch,
               vdb->mapVariableAndContextToIndex(d_primitive_vars,
                  getInteriorWithGhostsContext()),
               0,
               ghost_width_to_check,
               bbox,
               bscalarcase,
               d_bdry_face_val[refbdryloc]);
      }
#if (TESTING == 1)
      if (num_bad_values > 0) {
         tbox::perr << "\nConvDiff Boundary Test FAILED: \n"
                    << "     " << num_bad_values
                    << " bad VAL values found for\n"
                    << "     boundary type " << btype << " at location "
                    << bloc << std::endl;
      }
#endif

   }

}
