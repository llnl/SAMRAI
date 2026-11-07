/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Poisson using FAC
 *
 ************************************************************************/
#include "SAMRAI/solv/CellPoissonFACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/CellPoissonHypreSolver.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"


namespace SAMRAI {
namespace solv {

std::shared_ptr<pdat::CellVariable<double> >
CellPoissonFACOps::s_cell_scratch_var[SAMRAI::MAX_DIM_VAL];

std::shared_ptr<pdat::SideVariable<double> >
CellPoissonFACOps::s_flux_scratch_var[SAMRAI::MAX_DIM_VAL];

std::shared_ptr<pdat::OutersideVariable<double> >
CellPoissonFACOps::s_oflux_scratch_var[SAMRAI::MAX_DIM_VAL];

tbox::StartupShutdownManager::Handler
CellPoissonFACOps::s_finalize_handler(
   0,
   0,
   0,
   CellPoissonFACOps::finalizeCallback,
   tbox::StartupShutdownManager::priorityVariables);

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

void SAMRAI_F77_FUNC(compfluxvardc2d, COMPFLUXVARDC2D) (
   double* xflux,
   double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const int* dcgi,
   const int* dcgj,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void SAMRAI_F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D) (
   double* xflux,
   double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double& diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void SAMRAI_F77_FUNC(rbgswithfluxmaxvardcvarsf2d, RBGSWITHFLUXMAXVARDCVARSF2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const int* dcgi,
   const int* dcgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxcondcvarsf2d, RBGSWITHFLUXMAXCONDCVARSF2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double& dc,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf2d, RBGSWITHFLUXMAXVARDCCONSF2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const int* dcgi,
   const int* dcgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const double& scalar_field,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf2d, RBGSWITHFLUXMAXCONDCCONSF2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double& dc,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const double& scalar_field,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(compresvarsca2d, COMPRESVARSCA2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void SAMRAI_F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const double& scalar_field,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* dx);
void SAMRAI_F77_FUNC(ewingfixfluxvardc2d, EWINGFIXFLUXVARDC2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const int* dcgi,
   const int* dcgj,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);
void SAMRAI_F77_FUNC(ewingfixfluxcondc2d, EWINGFIXFLUXCONDC2D) (
   const double* xflux,
   const double* yflux,
   const int* fluxgi,
   const int* fluxgj,
   const double& diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);

void SAMRAI_F77_FUNC(compfluxvardc3d, COMPFLUXVARDC3D) (
   double* xflux,
   double* yflux,
   double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const double* zdiff_coef,
   const int* dcgi,
   const int* dcgj,
   const int* dcgk,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx);
void SAMRAI_F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D) (
   double* xflux,
   double* yflux,
   double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double& diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx);
void SAMRAI_F77_FUNC(rbgswithfluxmaxvardcvarsf3d, RBGSWITHFLUXMAXVARDCVARSF3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const double* zdiff_coef,
   const int* dcgi,
   const int* dcgj,
   const int* dcgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const int* scalar_field_gk,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxcondcvarsf3d, RBGSWITHFLUXMAXCONDCVARSF3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double& dc,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const int* scalar_field_gk,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf3d, RBGSWITHFLUXMAXVARDCCONSF3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const double* zdiff_coef,
   const int* dcgi,
   const int* dcgj,
   const int* dcgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   const double& scalar_field,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf3d, RBGSWITHFLUXMAXCONDCCONSF3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double& dc,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   const double& scalar_field,
   double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx,
   const int* offset,
   const double* maxres);
void SAMRAI_F77_FUNC(compresvarsca3d, COMPRESVARSCA3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const int* residualgk,
   const double* scalar_field,
   const int* scalar_field_gi,
   const int* scalar_field_gj,
   const int* scalar_field_gk,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx);
void SAMRAI_F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* rhs,
   const int* rhsgi,
   const int* rhsgj,
   const int* rhsgk,
   double* residual,
   const int* residualgi,
   const int* residualgj,
   const int* residualgk,
   const double& scalar_field,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* dx);
void SAMRAI_F77_FUNC(ewingfixfluxvardc3d, EWINGFIXFLUXVARDC3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double* xdiff_coef,
   const double* ydiff_coef,
   const double* zdiff_coef,
   const int* dcgi,
   const int* dcgj,
   const int* dcgk,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);
void SAMRAI_F77_FUNC(ewingfixfluxcondc3d, EWINGFIXFLUXCONDC3D) (
   const double* xflux,
   const double* yflux,
   const double* zflux,
   const int* fluxgi,
   const int* fluxgj,
   const int* fluxgk,
   const double& diff_coef,
   const double* soln,
   const int* solngi,
   const int* solngj,
   const int* solngk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const int* location_index,
   const int* ratio_to_coarser,
   const int* blower,
   const int* bupper,
   const double* dx);

}

/*
 ********************************************************************
 * Constructor.
 ********************************************************************
 */
#ifdef HAVE_HYPRE
CellPoissonFACOps::CellPoissonFACOps(
   const std::shared_ptr<CellPoissonHypreSolver>& hypre_solver,
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(object_name),
   d_ln_min(-1),
   d_ln_max(-1),
   d_poisson_spec(object_name + "::Poisson specs"),
   d_coarse_solver_choice("hypre"),
   d_cf_discretization("Ewing"),
   d_prolongation_method("CONSTANT_REFINE"),
   d_coarse_solver_tolerance(1.e-10),
   d_coarse_solver_max_iterations(20),
   d_residual_tolerance_during_smoothing(-1.0),
   d_flux_id(-1),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_hypre_solver(hypre_solver),
   d_physical_bc_coef(0),
   d_context(hier::VariableDatabase::getDatabase()->getContext(
                object_name + "::PRIVATE_CONTEXT")),
   d_cell_scratch_id(-1),
   d_flux_scratch_id(-1),
   d_oflux_scratch_id(-1),
   d_bc_helper(dim,
               d_object_name + "::bc helper"),
   d_enable_logging(false)
{
   buildObject(input_db);
}
#else
CellPoissonFACOps::CellPoissonFACOps(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(object_name),
   d_ln_min(-1),
   d_ln_max(-1),
   d_poisson_spec(object_name + "::Poisson specs"),
   d_coarse_solver_choice("redblack"),
   d_cf_discretization("Ewing"),
   d_prolongation_method("CONSTANT_REFINE"),
   d_coarse_solver_tolerance(1.e-8),
   d_coarse_solver_max_iterations(500),
   d_residual_tolerance_during_smoothing(-1.0),
   d_flux_id(-1),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_physical_bc_coef(0),
   d_context(hier::VariableDatabase::getDatabase()->getContext(
                object_name + "::PRIVATE_CONTEXT")),
   d_cell_scratch_id(-1),
   d_flux_scratch_id(-1),
   d_oflux_scratch_id(-1),
   d_bc_helper(dim,
               d_object_name + "::bc helper"),
   d_enable_logging(false)
{
   buildObject(input_db);
}
#endif

CellPoissonFACOps::~CellPoissonFACOps()
{
}

void
CellPoissonFACOps::buildObject(
   const std::shared_ptr<tbox::Database>& input_db)
{
   if (d_dim == tbox::Dimension(1) || d_dim > tbox::Dimension(3)) {
      TBOX_ERROR("CellPoissonFACOps : DIM == 1 or > 3 not implemented yet.\n");
   }

   t_restrict_solution = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::restrictSolution()");
   t_restrict_residual = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::restrictResidual()");
   t_prolong = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::prolongErrorAndCorrect()");
   t_smooth_error = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::smoothError()");
   t_solve_coarsest = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::solveCoarsestLevel()");
   t_compute_composite_residual = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::computeCompositeResidualOnLevel()");
   t_compute_residual_norm = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonFACOps::computeResidualNorm()");

   if (!s_cell_scratch_var[d_dim.getValue() - 1]) {
      TBOX_ASSERT(!s_cell_scratch_var[d_dim.getValue() - 1]);

      std::ostringstream ss;
      ss << "CellPoissonFACOps::private_cell_scratch" << d_dim.getValue();
      s_cell_scratch_var[d_dim.getValue() - 1].reset(
         new pdat::CellVariable<double>(d_dim, ss.str(), d_allocator));
      ss.str("");
      ss << "CellPoissonFACOps::private_flux_scratch" << d_dim.getValue();
      s_flux_scratch_var[d_dim.getValue() - 1].reset(
         new pdat::SideVariable<double>(d_dim, ss.str(),
                                        hier::IntVector::getOne(d_dim),
                                        d_allocator));
      ss.str("");
      ss << "CellPoissonFACOps::private_oflux_scratch" << d_dim.getValue();
      s_oflux_scratch_var[d_dim.getValue() - 1].reset(
         new pdat::OutersideVariable<double>(d_dim, ss.str(),
                                             d_allocator));
   }

   /*
    * Some variables initialized by default are overriden by input.
    */
   getFromInput(input_db);

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
   d_cell_scratch_id = vdb->
      registerVariableAndContext(s_cell_scratch_var[d_dim.getValue() - 1],
         d_context,
         hier::IntVector::getOne(d_dim));
   d_flux_scratch_id = vdb->
      registerVariableAndContext(s_flux_scratch_var[d_dim.getValue() - 1],
         d_context,
         hier::IntVector::getZero(d_dim));
   d_oflux_scratch_id = vdb->
      registerVariableAndContext(s_oflux_scratch_var[d_dim.getValue() - 1],
         d_context,
         hier::IntVector::getZero(d_dim));

   /*
    * Check input validity and correctness.
    */
   checkInputPatchDataIndices();
}

/*
 ************************************************************************
 * Read input parameters from database.
 ************************************************************************
 */
void
CellPoissonFACOps::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{
   if (input_db) {
      d_coarse_solver_choice =
         input_db->getStringWithDefault("coarse_solver_choice",
            d_coarse_solver_choice);
      if (!(d_coarse_solver_choice == "hypre" ||
            d_coarse_solver_choice == "redblack" ||
            d_coarse_solver_choice == "jacobi")) {
         INPUT_VALUE_ERROR("coarse_solver_choice");
      }

      d_coarse_solver_tolerance =
         input_db->getDoubleWithDefault("coarse_solver_tolerance",
            d_coarse_solver_tolerance);
      if (!(d_coarse_solver_tolerance > 0)) {
         INPUT_RANGE_ERROR("coarse_solver_tolerance");
      }

      d_coarse_solver_max_iterations =
         input_db->getIntegerWithDefault("coarse_solver_max_iterations",
            d_coarse_solver_max_iterations);
      if (!(d_coarse_solver_max_iterations >= 1)) {
         INPUT_RANGE_ERROR("coarse_solver_max_iterations");
      }

      d_cf_discretization =
         input_db->getStringWithDefault("cf_discretization", "Ewing");
      if (!(d_cf_discretization == "Ewing" ||
            d_cf_discretization == "CONSTANT_REFINE" ||
            d_cf_discretization == "LINEAR_REFINE" ||
            d_cf_discretization == "CONSERVATIVE_LINEAR_REFINE")) {
         INPUT_VALUE_ERROR("cf_discretization");
      }

      d_prolongation_method =
         input_db->getStringWithDefault("prolongation_method",
            "CONSTANT_REFINE");
      if (!(d_prolongation_method == "CONSTANT_REFINE" ||
            d_prolongation_method == "LINEAR_REFINE" ||
            d_prolongation_method == "CONSERVATIVE_LINEAR_REFINE")) {
         INPUT_VALUE_ERROR("prolongation_method");
      }

      d_enable_logging = input_db->getBoolWithDefault("enable_logging", false);
   }
}

/*
 ************************************************************************
 * FACOperatorStrategy virtual initializeOperatorState function.
 *
 * Set internal variables to correspond to the solution passed in.
 * Look up transfer operators.
 ************************************************************************
 */

void
CellPoissonFACOps::initializeOperatorState(
   const SAMRAIVectorReal<double>& solution,
   const SAMRAIVectorReal<double>& rhs)
{
   deallocateOperatorState();
   int ln;
   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   d_hierarchy = solution.getPatchHierarchy();
   d_ln_min = solution.getCoarsestLevelNumber();
   d_ln_max = solution.getFinestLevelNumber();
   d_hopscell.reset(new math::HierarchyCellDataOpsReal<double>(d_hierarchy,
         d_ln_min,
         d_ln_max));
   d_hopsside.reset(new math::HierarchySideDataOpsReal<double>(d_hierarchy,
         d_ln_min,
         d_ln_max));

#ifdef DEBUG_CHECK_ASSERTIONS

   if (d_physical_bc_coef == 0) {
      /*
       * It's an error not to have bc object set.
       * Note that the bc object cannot be passed in through
       * the argument because the interface is inherited.
       */
      TBOX_ERROR(
         d_object_name << ": No physical bc object in\n"
                       << "CellPoissonFACOps::initializeOperatorState\n"
                       << "You must use "
                       << "CellPoissonFACOps::setPhysicalBcCoefObject\n"
                       << "to set one before calling initializeOperatorState\n");
   }

   if (solution.getNumberOfComponents() != 1) {
      TBOX_WARNING(d_object_name
         << ": Solution vector has multiple components.\n"
         << "Solver is for component 0 only.\n");
   }
   if (rhs.getNumberOfComponents() != 1) {
      TBOX_WARNING(d_object_name
         << ": RHS vector has multiple components.\n"
         << "Solver is for component 0 only.\n");
   }

   /*
    * Make sure that solution and rhs data
    *   are of correct type
    *   are allocated
    *   has sufficient ghost width
    */
   std::shared_ptr<hier::Variable> var;
   {
      vdb->mapIndexToVariable(rhs.getComponentDescriptorIndex(0),
         var);
      if (!var) {
         TBOX_ERROR(d_object_name << ": RHS component does not\n"
                                  << "correspond to a variable.\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
         SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(var));
      TBOX_ASSERT(cell_var);
   }
   {
      vdb->mapIndexToVariable(solution.getComponentDescriptorIndex(0),
         var);
      if (!var) {
         TBOX_ERROR(d_object_name << ": Solution component does not\n"
                                  << "correspond to a variable.\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
         SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(var));
      TBOX_ASSERT(cell_var);
   }
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level_ptr(
         d_hierarchy->getPatchLevel(ln));
      hier::PatchLevel& level = *level_ptr;
      for (hier::PatchLevel::iterator pi(level.begin());
           pi != level.end(); ++pi) {
         hier::Patch& patch = **pi;
         std::shared_ptr<hier::PatchData> fd(
            patch.getPatchData(rhs.getComponentDescriptorIndex(0)));
         if (fd) {
            /*
             * Some data checks can only be done if the data already exists.
             */
            std::shared_ptr<pdat::CellData<double> > cd(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(fd));
            TBOX_ASSERT(cd);
            if (cd->getDepth() > 1) {
               TBOX_WARNING(d_object_name
                  << ": RHS data has multiple depths.\n"
                  << "Solver is for depth 0 only.\n");
            }
         }
         std::shared_ptr<hier::PatchData> ud(
            patch.getPatchData(solution.getComponentDescriptorIndex(0)));
         if (ud) {
            /*
             * Some data checks can only be done if the data already exists.
             */
            std::shared_ptr<pdat::CellData<double> > cd(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(ud));
            TBOX_ASSERT(cd);
            if (cd->getDepth() > 1) {
               TBOX_WARNING(d_object_name
                  << ": Solution data has multiple depths.\n"
                  << "Solver is for depth 0 only.\n");
            }
            if (cd->getGhostCellWidth() < hier::IntVector::getOne(d_dim)) {
               TBOX_ERROR(d_object_name
                  << ": Solution data has insufficient ghost width\n");
            }
         }
      }
   }

   /*
    * Solution and rhs must have some similar properties.
    */
   if (rhs.getPatchHierarchy() != d_hierarchy
       || rhs.getCoarsestLevelNumber() != d_ln_min
       || rhs.getFinestLevelNumber() != d_ln_max) {
      TBOX_ERROR(d_object_name << ": solution and rhs do not have\n"
                               << "the same set of patch levels.\n");
   }

#endif

   /*
    * Initialize the coarse-fine boundary description for the
    * hierarchy.
    */
   d_cf_boundary.resize(d_hierarchy->getNumberOfLevels());

   hier::IntVector max_gcw(d_dim, 1);
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_cf_boundary[ln].reset(
         new hier::CoarseFineBoundary(*d_hierarchy,
            ln,
            max_gcw));
   }
#ifdef HAVE_HYPRE
   if (d_coarse_solver_choice == "hypre") {
      d_hypre_solver->initializeSolverState(d_hierarchy, d_ln_min);
      /*
       * Share the boundary condition object with the hypre solver
       * to make sure that boundary condition settings are consistent
       * between the two objects.
       */
      d_hypre_solver->setPhysicalBcCoefObject(d_physical_bc_coef);
      d_hypre_solver->setMatrixCoefficients(d_poisson_spec);
   }
#endif

   /*
    * Get the transfer operators.
    * Flux coarsening is conservative.
    * Cell (solution, error, etc) coarsening is conservative.
    * Cell refinement from same level is constant refinement.
    * Cell refinement from coarser level is chosen by the
    *   choice of coarse-fine discretization, d_cf_discretization,
    *   which should be set to either "Ewing" or one of the
    *   acceptable strings for looking up the refine operator.
    */
   std::shared_ptr<geom::CartesianGridGeometry> geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
         d_hierarchy->getGridGeometry()));
   TBOX_ASSERT(geometry);
   std::shared_ptr<hier::Variable> variable;

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_prolongation_refine_operator =
      geometry->lookupRefineOperator(variable,
         d_prolongation_method);

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_urestriction_coarsen_operator =
      d_rrestriction_coarsen_operator =
         geometry->lookupCoarsenOperator(variable,
            "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_oflux_scratch_id, variable);
   d_flux_coarsen_operator =
      geometry->lookupCoarsenOperator(variable,
         "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_refine_operator =
      geometry->lookupRefineOperator(variable,
         d_cf_discretization == "Ewing" ?
         "CONSTANT_REFINE" : d_cf_discretization);

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_nocoarse_refine_operator =
      geometry->lookupRefineOperator(variable,
         "CONSTANT_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_prolongation_refine_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find prolongation refine operator" << std::endl);
   }
   if (!d_urestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find restriction coarsening operator" << std::endl);
   }
   if (!d_rrestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find restriction coarsening operator" << std::endl);
   }
   if (!d_flux_coarsen_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find flux coarsening operator" << std::endl);
   }
   if (!d_ghostfill_refine_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find ghost filling refinement operator" << std::endl);
   }
   if (!d_ghostfill_nocoarse_refine_operator) {
      TBOX_ERROR(d_object_name
         << ": Cannot find ghost filling refinement operator" << std::endl);
   }
#endif

   for (ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->
      allocatePatchData(d_oflux_scratch_id);
   }

   /*
    * Make space for saving communication schedules.
    * There is no need to delete the old schedules first
    * because we have deallocated the solver state above.
    */
   d_prolongation_refine_schedules.resize(d_ln_max + 1);
   d_ghostfill_refine_schedules.resize(d_ln_max + 1);
   d_ghostfill_nocoarse_refine_schedules.resize(d_ln_max + 1);
   d_urestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_rrestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_flux_coarsen_schedules.resize(d_ln_max + 1);

   d_prolongation_refine_algorithm.reset(
      new xfer::RefineAlgorithm());
   d_urestriction_coarsen_algorithm.reset(
      new xfer::CoarsenAlgorithm(d_dim));
   d_rrestriction_coarsen_algorithm.reset(
      new xfer::CoarsenAlgorithm(d_dim));
   d_flux_coarsen_algorithm.reset(
      new xfer::CoarsenAlgorithm(d_dim));
   d_ghostfill_refine_algorithm.reset(
      new xfer::RefineAlgorithm());
   d_ghostfill_nocoarse_refine_algorithm.reset(
      new xfer::RefineAlgorithm());

   d_prolongation_refine_algorithm->registerRefine(
      d_cell_scratch_id,
      solution.getComponentDescriptorIndex(0),
      d_cell_scratch_id,
      d_prolongation_refine_operator);
   d_urestriction_coarsen_algorithm->registerCoarsen(
      solution.getComponentDescriptorIndex(0),
      solution.getComponentDescriptorIndex(0),
      d_urestriction_coarsen_operator);
   d_rrestriction_coarsen_algorithm->registerCoarsen(
      rhs.getComponentDescriptorIndex(0),
      rhs.getComponentDescriptorIndex(0),
      d_rrestriction_coarsen_operator);
   d_ghostfill_refine_algorithm->registerRefine(
      solution.getComponentDescriptorIndex(0),
      solution.getComponentDescriptorIndex(0),
      solution.getComponentDescriptorIndex(0),
      d_ghostfill_refine_operator);
   d_flux_coarsen_algorithm->registerCoarsen(
      ((d_flux_id != -1) ? d_flux_id : d_flux_scratch_id),
      d_oflux_scratch_id,
      d_flux_coarsen_operator);
   d_ghostfill_nocoarse_refine_algorithm->registerRefine(
      solution.getComponentDescriptorIndex(0),
      solution.getComponentDescriptorIndex(0),
      solution.getComponentDescriptorIndex(0),
      d_ghostfill_nocoarse_refine_operator);

   for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln) {

      std::shared_ptr<xfer::PatchLevelFullFillPattern> fill_pattern(
         std::make_shared<xfer::PatchLevelFullFillPattern>());
      d_prolongation_refine_schedules[dest_ln] =
         d_prolongation_refine_algorithm->
         createSchedule(fill_pattern,
            d_hierarchy->getPatchLevel(dest_ln),
            std::shared_ptr<hier::PatchLevel>(),
            dest_ln - 1,
            d_hierarchy,
            &d_bc_helper);
      if (!d_prolongation_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name
            << ": Cannot create a refine schedule for prolongation!\n");
      }
      d_ghostfill_refine_schedules[dest_ln] =
         d_ghostfill_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dest_ln),
            dest_ln - 1,
            d_hierarchy,
            &d_bc_helper);
      if (!d_ghostfill_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name
            << ": Cannot create a refine schedule for ghost filling!\n");
      }
      d_ghostfill_nocoarse_refine_schedules[dest_ln] =
         d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dest_ln),
            &d_bc_helper);
      if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
         TBOX_ERROR(
            d_object_name
            << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
      }
   }
   for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln) {
      d_urestriction_coarsen_schedules[dest_ln] =
         d_urestriction_coarsen_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dest_ln),
            d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_urestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name
            << ": Cannot create a coarsen schedule for U restriction!\n");
      }
      d_rrestriction_coarsen_schedules[dest_ln] =
         d_rrestriction_coarsen_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dest_ln),
            d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_rrestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name
            << ": Cannot create a coarsen schedule for R restriction!\n");
      }
      d_flux_coarsen_schedules[dest_ln] =
         d_flux_coarsen_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dest_ln),
            d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_flux_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name
            << ": Cannot create a coarsen schedule for flux transfer!\n");
      }
   }
   d_ghostfill_nocoarse_refine_schedules[d_ln_min] =
      d_ghostfill_nocoarse_refine_algorithm->createSchedule(
         d_hierarchy->getPatchLevel(d_ln_min),
         &d_bc_helper);
   if (!d_ghostfill_nocoarse_refine_schedules[d_ln_min]) {
      TBOX_ERROR(
         d_object_name
         << ": Cannot create a refine schedule for ghost filling on bottom level!\n");
   }
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual deallocateOperatorState
 * function.  Deallocate internal hierarchy-dependent data.
 * State is allocated iff hierarchy is set.
 ********************************************************************
 */

void
CellPoissonFACOps::deallocateOperatorState()
{
   if (d_hierarchy) {
      int ln;
      for (ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
         d_hierarchy->getPatchLevel(ln)->
         deallocatePatchData(d_oflux_scratch_id);
      }
      d_cf_boundary.resize(0);
#ifdef HAVE_HYPRE
      d_hypre_solver->deallocateSolverState();
#endif
      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;

      d_prolongation_refine_algorithm.reset();
      d_prolongation_refine_schedules.clear();

      d_urestriction_coarsen_algorithm.reset();
      d_urestriction_coarsen_schedules.clear();

      d_rrestriction_coarsen_algorithm.reset();
      d_rrestriction_coarsen_schedules.clear();

      d_flux_coarsen_algorithm.reset();
      d_flux_coarsen_schedules.clear();

      d_ghostfill_refine_algorithm.reset();
      d_ghostfill_refine_schedules.clear();

      d_ghostfill_nocoarse_refine_algorithm.reset();
      d_ghostfill_nocoarse_refine_schedules.clear();

   }
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual postprocessOneCycle function.
 ********************************************************************
 */

void
CellPoissonFACOps::postprocessOneCycle(
   int fac_cycle_num,
   const SAMRAIVectorReal<double>& current_soln,
   const SAMRAIVectorReal<double>& residual)
{
   NULL_USE(current_soln);
   NULL_USE(residual);

   if (d_enable_logging) {
      if (d_preconditioner) {
         /*
          * Output convergence progress.  This is probably only appropriate
          * if the solver is NOT being used as a preconditioner.
          */
         double avg_factor, final_factor;
         d_preconditioner->getConvergenceFactors(avg_factor, final_factor);
         tbox::plog
         << "iter=" << std::setw(4) << fac_cycle_num
         << " resid=" << d_preconditioner->getResidualNorm()
         << " net conv=" << d_preconditioner->getNetConvergenceFactor()
         << " final conv=" << d_preconditioner->getNetConvergenceFactor()
         << " avg conv=" << d_preconditioner->getAvgConvergenceFactor()
         << std::endl;
      }
   }
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual restrictSolution function.
 * After restricting solution, update ghost cells of the affected
 * level.
 ********************************************************************
 */

void
CellPoissonFACOps::restrictSolution(
   const SAMRAIVectorReal<double>& s,
   SAMRAIVectorReal<double>& d,
   int dest_ln)
{
   t_restrict_solution->start();

   xeqScheduleURestriction(d.getComponentDescriptorIndex(0),
      s.getComponentDescriptorIndex(0),
      dest_ln);

   d_bc_helper.setHomogeneousBc(false);
   d_bc_helper.setTargetDataId(d.getComponentDescriptorIndex(0));

   if (dest_ln == d_ln_min) {
      xeqScheduleGhostFillNoCoarse(d.getComponentDescriptorIndex(0),
         dest_ln);
   } else {
      xeqScheduleGhostFill(d.getComponentDescriptorIndex(0),
         dest_ln);
   }

   t_restrict_solution->stop();
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual restrictresidual function.
 ********************************************************************
 */

void
CellPoissonFACOps::restrictResidual(
   const SAMRAIVectorReal<double>& s,
   SAMRAIVectorReal<double>& d,
   int dest_ln)
{

   t_restrict_residual->start();

   xeqScheduleRRestriction(d.getComponentDescriptorIndex(0),
      s.getComponentDescriptorIndex(0),
      dest_ln);

   t_restrict_residual->stop();
}

/*
 ***********************************************************************
 * FACOperatorStrategy virtual prolongErrorAndCorrect function.
 * After the prolongation, we set the physical boundary condition
 * for the correction, which is zero.  Other ghost cell values,
 * which are preset to zero, need not be set.
 ***********************************************************************
 */

void
CellPoissonFACOps::prolongErrorAndCorrect(
   const SAMRAIVectorReal<double>& s,
   SAMRAIVectorReal<double>& d,
   int dest_ln)
{
   t_prolong->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (s.getPatchHierarchy() != d_hierarchy
       || d.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
         "internal state hierarchy." << std::endl);
   }
#endif

   std::shared_ptr<hier::PatchLevel> fine_level(
      d_hierarchy->getPatchLevel(dest_ln));

   /*
    * Data is prolonged into the scratch space corresponding
    * to index d_cell_scratch_id and allocated here.
    */
   fine_level->allocatePatchData(d_cell_scratch_id);

   /*
    * Refine solution into scratch space to fill the fine level
    * interior in the scratch space, then use that refined data
    * to correct the fine level error.
    */
   d_bc_helper.setTargetDataId(d_cell_scratch_id);
   d_bc_helper.setHomogeneousBc(true);
   const int src_index = s.getComponentDescriptorIndex(0);
   xeqScheduleProlongation(d_cell_scratch_id,
      src_index,
      d_cell_scratch_id,
      dest_ln);

   /*
    * Add the refined error in the scratch space
    * to the error currently residing in the destination level.
    */
   math::HierarchyCellDataOpsReal<double>
   hierarchy_math_ops(d_hierarchy, dest_ln, dest_ln);
   const int dst_index = d.getComponentDescriptorIndex(0);
   hierarchy_math_ops.add(dst_index, dst_index, d_cell_scratch_id);

   fine_level->deallocatePatchData(d_cell_scratch_id);

   t_prolong->stop();

}

/*
 ********************************************************************
 ********************************************************************
 */

void
CellPoissonFACOps::smoothError(
   SAMRAIVectorReal<double>& data,
   const SAMRAIVectorReal<double>& residual,
   int ln,
   int num_sweeps)
{

   t_smooth_error->start();

   checkInputPatchDataIndices();
   smoothErrorByRedBlack(data,
      residual,
      ln,
      num_sweeps,
      d_residual_tolerance_during_smoothing);

   t_smooth_error->stop();
}

/*
 ********************************************************************
 * Workhorse function to smooth error using red-black
 * Gauss-Seidel iterations.
 ********************************************************************
 */

void
CellPoissonFACOps::smoothErrorByRedBlack(
   SAMRAIVectorReal<double>& data,
   const SAMRAIVectorReal<double>& residual,
   int ln,
   int num_sweeps,
   double residual_tolerance)
{

   checkInputPatchDataIndices();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (data.getPatchHierarchy() != d_hierarchy
       || residual.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
         "internal hierarchy." << std::endl);
   }
#endif
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(ln));

   const int data_id = data.getComponentDescriptorIndex(0);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   d_bc_helper.setTargetDataId(data_id);
   d_bc_helper.setHomogeneousBc(true);
   xeqScheduleGhostFillNoCoarse(data_id, ln);

   if (ln > d_ln_min) {
      /*
       * Perform a one-time transfer of data from coarser level,
       * to fill ghost boundaries that will not change through
       * the smoothing loop.
       */
      xeqScheduleGhostFill(data_id, ln);
   }

   /*
    * Smooth the number of sweeps specified or until
    * the convergence is satisfactory.
    */
   int isweep;
   double red_maxres, blk_maxres, maxres = 0;
   red_maxres = blk_maxres = residual_tolerance + 1;
   /*
    * Instead of checking residual convergence globally,
    * we check the not_converged flag.  This avoids possible
    * round-off errors affecting different processes differently,
    * leading to disagreement on whether to continue smoothing.
    */
   int not_converged = 1;
   for (isweep = 0; isweep < num_sweeps && not_converged; ++isweep) {
      red_maxres = blk_maxres = 0;

      // Red sweep.
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      for (hier::PatchLevel::iterator pi(level->begin());
           pi != level->end(); ++pi) {
         const std::shared_ptr<hier::Patch>& patch = *pi;

         bool deallocate_flux_data_when_done = false;
         if (flux_id == d_flux_scratch_id) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if (!patch->checkAllocated(flux_id)) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         std::shared_ptr<pdat::CellData<double> > err_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               data.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > residual_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               residual.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::SideData<double> > flux_data(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(flux_id)));

         TBOX_ASSERT(err_data);
         TBOX_ASSERT(residual_data);
         TBOX_ASSERT(flux_data);

         computeFluxOnPatch(
            *patch,
            level->getRatioToCoarserLevel(),
            *err_data,
            *flux_data);

         redOrBlackSmoothingOnPatch(*patch,
            *flux_data,
            *residual_data,
            *err_data,
            'r',
            &red_maxres);

         if (deallocate_flux_data_when_done) {
            patch->deallocatePatchData(flux_id);
         }
      }        // End patch number *pi
      xeqScheduleGhostFillNoCoarse(data_id, ln);

      // Black sweep.
      for (hier::PatchLevel::iterator pi(level->begin());
           pi != level->end(); ++pi) {
         const std::shared_ptr<hier::Patch>& patch = *pi;

         bool deallocate_flux_data_when_done = false;
         if (flux_id == d_flux_scratch_id) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if (!patch->checkAllocated(flux_id)) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         std::shared_ptr<pdat::CellData<double> > err_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               data.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > residual_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               residual.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::SideData<double> > flux_data(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData(flux_id)));

         TBOX_ASSERT(err_data);
         TBOX_ASSERT(residual_data);
         TBOX_ASSERT(flux_data);

         computeFluxOnPatch(
            *patch,
            level->getRatioToCoarserLevel(),
            *err_data,
            *flux_data);

         redOrBlackSmoothingOnPatch(*patch,
            *flux_data,
            *residual_data,
            *err_data,
            'b',
            &blk_maxres);

         if (deallocate_flux_data_when_done) {
            patch->deallocatePatchData(flux_id);
         }
      }        // End patch number *pi
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      if (residual_tolerance >= 0.0) {
         /*
          * Check for early end of sweeps due to convergence
          * only if it is numerically possible (user gave a
          * non negative value for residual tolerance).
          */
         maxres = tbox::MathUtilities<double>::Max(red_maxres, blk_maxres);
         not_converged = maxres > residual_tolerance;
         const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
         if (mpi.getSize() > 1) {
            mpi.AllReduce(&not_converged, 1, MPI_MAX);
         }
      }
   }        // End sweep number isweep
   if (d_enable_logging) tbox::plog
      << d_object_name << " RBGS smoothing maxres = " << maxres << "\n"
      << "  after " << isweep << " sweeps.\n";

}

/*
 ********************************************************************
 * Fix flux on coarse-fine boundaries computed from a
 * constant-refine interpolation of coarse level data.
 ********************************************************************
 */

void
CellPoissonFACOps::ewingFixFlux(
   const hier::Patch& patch,
   const pdat::CellData<double>& soln_data,
   pdat::SideData<double>& flux_data,
   const hier::IntVector& ratio_to_coarser) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(d_dim, patch, soln_data, flux_data,
      ratio_to_coarser);

   const int patch_ln = patch.getPatchLevelNumber();
   const hier::GlobalId id = patch.getGlobalId();
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();
   const hier::Box& patch_box(patch.getBox());
   const hier::Index& plower = patch_box.lower();
   const hier::Index& pupper = patch_box.upper();

   hier::IntVector block_ratio(ratio_to_coarser);
   if (block_ratio.getNumBlocks() != 1) {
      block_ratio = hier::IntVector(d_dim);
      hier::BlockId::block_t b = patch_box.getBlockId().getBlockValue();
      for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
         block_ratio[d] = ratio_to_coarser(b,d);
      }
   }

   const std::vector<hier::BoundaryBox>& bboxes =
      d_cf_boundary[patch_ln]->getBoundaries(id, 1);
   int bn, nboxes = static_cast<int>(bboxes.size());

   if (d_poisson_spec.dIsVariable()) {

      std::shared_ptr<pdat::SideData<double> > diffcoef_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_poisson_spec.getDPatchDataId())));

      TBOX_ASSERT(diffcoef_data);

      for (bn = 0; bn < nboxes; ++bn) {
         const hier::BoundaryBox& boundary_box = bboxes[bn];

         TBOX_ASSERT(boundary_box.getBoundaryType() == 1);

         const hier::Box& bdry_box = boundary_box.getBox();
         const hier::Index& blower = bdry_box.lower();
         const hier::Index& bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(ewingfixfluxvardc2d, EWINGFIXFLUXVARDC2D) (
               flux_data.getPointer(0), flux_data.getPointer(1),
               &flux_data.getGhostCellWidth()[0],
               &flux_data.getGhostCellWidth()[1],
               diffcoef_data->getPointer(0), diffcoef_data->getPointer(1),
               &diffcoef_data->getGhostCellWidth()[0],
               &diffcoef_data->getGhostCellWidth()[1],
               soln_data.getPointer(),
               &soln_data.getGhostCellWidth()[0],
               &soln_data.getGhostCellWidth()[1],
               &plower[0], &pupper[0], &plower[1], &pupper[1],
               &location_index,
               &block_ratio[0],
               &blower[0], &bupper[0],
               dx);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(ewingfixfluxvardc3d, EWINGFIXFLUXVARDC3D) (
               flux_data.getPointer(0),
               flux_data.getPointer(1),
               flux_data.getPointer(2),
               &flux_data.getGhostCellWidth()[0],
               &flux_data.getGhostCellWidth()[1],
               &flux_data.getGhostCellWidth()[2],
               diffcoef_data->getPointer(0),
               diffcoef_data->getPointer(1),
               diffcoef_data->getPointer(2),
               &diffcoef_data->getGhostCellWidth()[0],
               &diffcoef_data->getGhostCellWidth()[1],
               &diffcoef_data->getGhostCellWidth()[2],
               soln_data.getPointer(),
               &soln_data.getGhostCellWidth()[0],
               &soln_data.getGhostCellWidth()[1],
               &soln_data.getGhostCellWidth()[2],
               &plower[0], &pupper[0],
               &plower[1], &pupper[1],
               &plower[2], &pupper[2],
               &location_index,
               &block_ratio[0],
               &blower[0], &bupper[0],
               dx);
         } else {
            TBOX_ERROR("CellPoissonFACOps : DIM > 3 not supported" << std::endl);
         }

      }
   } else {

      const double diffcoef_constant = d_poisson_spec.getDConstant();

      for (bn = 0; bn < nboxes; ++bn) {
         const hier::BoundaryBox& boundary_box = bboxes[bn];

         TBOX_ASSERT(boundary_box.getBoundaryType() == 1);

         const hier::Box& bdry_box = boundary_box.getBox();
         const hier::Index& blower = bdry_box.lower();
         const hier::Index& bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(ewingfixfluxcondc2d, EWINGFIXFLUXCONDC2D) (
               flux_data.getPointer(0), flux_data.getPointer(1),
               &flux_data.getGhostCellWidth()[0],
               &flux_data.getGhostCellWidth()[1],
               diffcoef_constant,
               soln_data.getPointer(),
               &soln_data.getGhostCellWidth()[0],
               &soln_data.getGhostCellWidth()[1],
               &plower[0], &pupper[0],
               &plower[1], &pupper[1],
               &location_index,
               &block_ratio[0],
               &blower[0], &bupper[0],
               dx);
         } else if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(ewingfixfluxcondc3d, EWINGFIXFLUXCONDC3D) (
               flux_data.getPointer(0),
               flux_data.getPointer(1),
               flux_data.getPointer(2),
               &flux_data.getGhostCellWidth()[0],
               &flux_data.getGhostCellWidth()[1],
               &flux_data.getGhostCellWidth()[2],
               diffcoef_constant,
               soln_data.getPointer(),
               &soln_data.getGhostCellWidth()[0],
               &soln_data.getGhostCellWidth()[1],
               &soln_data.getGhostCellWidth()[2],
               &plower[0], &pupper[0],
               &plower[1], &pupper[1],
               &plower[2], &pupper[2],
               &location_index,
               &block_ratio[0],
               &blower[0], &bupper[0],
               dx);
         }
      }
   }
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual solveCoarsestLevel
 * function
 ********************************************************************
 */

int
CellPoissonFACOps::solveCoarsestLevel(
   SAMRAIVectorReal<double>& data,
   const SAMRAIVectorReal<double>& residual,
   int coarsest_ln)
{
   t_solve_coarsest->start();

   checkInputPatchDataIndices();

   int return_value = 0;

   if (d_coarse_solver_choice == "jacobi") {
      d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
      smoothError(data,
         residual,
         coarsest_ln,
         d_coarse_solver_max_iterations);
      d_residual_tolerance_during_smoothing = -1.0;
   } else if (d_coarse_solver_choice == "redblack") {
      d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
      smoothError(data,
         residual,
         coarsest_ln,
         d_coarse_solver_max_iterations);
      d_residual_tolerance_during_smoothing = -1.0;
   } else if (d_coarse_solver_choice == "hypre") {
#ifndef HAVE_HYPRE
      TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
                               << d_coarse_solver_choice
                               << "' unavailable in "
                               << "scapCellPoissonOps::solveCoarsestLevel."
                               << std::endl);
#else
      return_value = solveCoarsestLevel_HYPRE(data, residual, coarsest_ln);
#endif
   } else {
      TBOX_ERROR(
         d_object_name << ": Bad coarse level solver choice '"
                       << d_coarse_solver_choice
                       << "' in scapCellPoissonOps::solveCoarsestLevel."
                       << std::endl);
   }

   xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
      coarsest_ln);

   t_solve_coarsest->stop();

   return return_value;
}

#ifdef HAVE_HYPRE
/*
 ********************************************************************
 * Solve coarsest level using Hypre
 * We only solve for the error, so we always use homogeneous bc.
 ********************************************************************
 */

int
CellPoissonFACOps::solveCoarsestLevel_HYPRE(
   SAMRAIVectorReal<double>& data,
   const SAMRAIVectorReal<double>& residual,
   int coarsest_ln)
{
   NULL_USE(coarsest_ln);

#ifndef HAVE_HYPRE
   TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
                            << d_coarse_solver_choice
                            << "' unavailable in "
                            << "CellPoissonFACOps::solveCoarsestLevel_HYPRE."
                            << std::endl);

   return 0;

#else

   checkInputPatchDataIndices();
   d_hypre_solver->setStoppingCriteria(d_coarse_solver_max_iterations,
      d_coarse_solver_tolerance);
   const int solver_ret =
      d_hypre_solver->solveSystem(
         data.getComponentDescriptorIndex(0),
         residual.getComponentDescriptorIndex(0),
         true);
   /*
    * Present data on the solve.
    * The Hypre solver returns 0 if converged.
    */
   if (d_enable_logging) tbox::plog
      << d_object_name << " Hypre solve " << (solver_ret ? "" : "NOT ")
      << "converged\n"
      << "\titerations: " << d_hypre_solver->getNumberOfIterations() << "\n"
      << "\tresidual: " << d_hypre_solver->getRelativeResidualNorm() << "\n";

   return !solver_ret;

#endif

}
#endif

/*
 ********************************************************************
 * FACOperatorStrategy virtual
 * computeCompositeResidualOnLevel function
 ********************************************************************
 */

void
CellPoissonFACOps::computeCompositeResidualOnLevel(
   SAMRAIVectorReal<double>& residual,
   const SAMRAIVectorReal<double>& solution,
   const SAMRAIVectorReal<double>& rhs,
   int ln,
   bool error_equation_indicator)
{
   t_compute_composite_residual->start();

   checkInputPatchDataIndices();
#ifdef DEBUG_CHECK_ASSERTIONS
   if (residual.getPatchHierarchy() != d_hierarchy
       || solution.getPatchHierarchy() != d_hierarchy
       || rhs.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
         "internal hierarchy." << std::endl);
   }
#endif
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(ln));

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */
   const int soln_id = solution.getComponentDescriptorIndex(0);
   d_bc_helper.setTargetDataId(soln_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. Residual data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute residual data from flux.
    */

   /* S1. Fill solution ghost data. */
   {
      if (ln > d_ln_min) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(soln_id, ln);
      } else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(soln_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done = false;
   if (flux_id == d_flux_scratch_id) {
      if (!level->checkAllocated(flux_id)) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      const std::shared_ptr<hier::Patch>& patch = *pi;

      std::shared_ptr<pdat::CellData<double> > soln_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch->getPatchData(flux_id)));

      TBOX_ASSERT(soln_data);
      TBOX_ASSERT(flux_data);

      computeFluxOnPatch(
         *patch,
         level->getRatioToCoarserLevel(),
         *soln_data,
         *flux_data);

   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if (ln < d_ln_max) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Compute residual on patches in level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      const std::shared_ptr<hier::Patch>& patch = *pi;
      std::shared_ptr<pdat::CellData<double> > soln_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > rhs_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            rhs.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > residual_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            residual.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch->getPatchData(flux_id)));

      TBOX_ASSERT(soln_data);
      TBOX_ASSERT(rhs_data);
      TBOX_ASSERT(residual_data);
      TBOX_ASSERT(flux_data);

      computeResidualOnPatch(*patch,
         *flux_data,
         *soln_data,
         *rhs_data,
         *residual_data);

      if (ln > d_ln_min) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute residual"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         std::shared_ptr<pdat::OutersideData<double> > oflux_data(
            SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
               patch->getPatchData(d_oflux_scratch_id)));

         TBOX_ASSERT(oflux_data);

         oflux_data->copy(*flux_data);
      }
   }

   if (deallocate_flux_data_when_done) {
      level->deallocatePatchData(flux_id);
   }

   t_compute_composite_residual->stop();
}

/*
 ********************************************************************
 * FACOperatorStrategy virtual computeResidualNorm
 * function
 ********************************************************************
 */

double
CellPoissonFACOps::computeResidualNorm(
   const SAMRAIVectorReal<double>& residual,
   int fine_ln,
   int coarse_ln)
{

   if (coarse_ln != residual.getCoarsestLevelNumber() ||
       fine_ln != residual.getFinestLevelNumber()) {
      TBOX_ERROR("CellPoissonFACOps::computeResidualNorm() is not\n"
         << "set up to compute residual except on the range of\n"
         << "levels defining the vector.\n");
   }
   t_compute_residual_norm->start();
   /*
    * The residual vector was cloned from vectors that has
    * the proper weights associated with them, so we do not
    * have to explicitly weight the residuals.
    *
    * maxNorm: not good to use because Hypre's norm does not
    *   correspond to it.  Also maybe too sensitive to spikes.
    * L2Norm: maybe good.  But does not correspond to the
    *   scale of the quantity.
    * L1Norm: maybe good.  Correspond to scale of quantity,
    *   but may be too insensitive to spikes.
    * RMSNorm: maybe good.
    */
   double norm = residual.RMSNorm();
   t_compute_residual_norm->stop();
   return norm;
}

/*
 ********************************************************************
 * Compute the vector weight and put it at a specified patch data
 * index.
 ********************************************************************
 */

void
CellPoissonFACOps::computeVectorWeights(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int weight_id,
   int coarsest_ln,
   int finest_ln) const
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

   if (coarsest_ln == -1) coarsest_ln = 0;
   if (finest_ln == -1) finest_ln = hierarchy->getFinestLevelNumber();
   if (finest_ln < coarsest_ln) {
      TBOX_ERROR(d_object_name
         << ": Illegal level number range.  finest_ln < coarsest_ln."
         << std::endl);
   }

   int ln;
   for (ln = finest_ln; ln >= coarsest_ln; --ln) {

      /*
       * On every level, first assign cell volume to vector weight.
       */

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));

         TBOX_ASSERT(patch_geometry);

         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (d_dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
         }

         if (d_dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
         }

         std::shared_ptr<pdat::CellData<double> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(weight_id)));
         TBOX_ASSERT(w);
         w->fillAll(cell_vol);
      }

      /*
       * On all but the finest level, assign 0 to vector
       * weight to cells covered by finer cells.
       */

      if (ln < finest_ln) {

         /*
          * First get the boxes that describe index space of the next finer
          * level and coarsen them to describe corresponding index space
          * at this level.
          */

         std::shared_ptr<hier::PatchLevel> next_finer_level(
            hierarchy->getPatchLevel(ln + 1));
         hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
         hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
         coarsen_ratio /= level->getRatioToLevelZero();
         coarsened_boxes.coarsen(coarsen_ratio);

         /*
          * Then set vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for (hier::PatchLevel::iterator p(level->begin());
              p != level->end(); ++p) {

            const std::shared_ptr<hier::Patch>& patch = *p;
            for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
                 i != coarsened_boxes.end(); ++i) {

               hier::Box intersection = *i * (patch->getBox());
               if (!intersection.empty()) {
                  std::shared_ptr<pdat::CellData<double> > w(
                     SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(weight_id)));
                  TBOX_ASSERT(w);
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
         }  // loop over patches in level
      }  // all levels except finest
   }  // loop over levels
}

/*
 ********************************************************************
 * Check the validity and correctness of input data for this class.
 ********************************************************************
 */

void
CellPoissonFACOps::checkInputPatchDataIndices() const
{
   /*
    * Check input validity and correctness.
    */
   hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());

   if (!d_poisson_spec.dIsConstant()
       && d_poisson_spec.getDPatchDataId() != -1) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_poisson_spec.getDPatchDataId(), var);
      std::shared_ptr<pdat::SideVariable<double> > diffcoef_var(
         SAMRAI_SHARED_PTR_CAST<pdat::SideVariable<double>, hier::Variable>(var));

      TBOX_ASSERT(diffcoef_var);
   }

   if (!d_poisson_spec.cIsConstant() && !d_poisson_spec.cIsZero()) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_poisson_spec.getCPatchDataId(), var);
      std::shared_ptr<pdat::CellVariable<double> > scalar_field_var(
         SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(var));

      TBOX_ASSERT(scalar_field_var);
   }

   if (d_flux_id != -1) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_flux_id, var);
      std::shared_ptr<pdat::SideVariable<double> > flux_var(
         SAMRAI_SHARED_PTR_CAST<pdat::SideVariable<double>, hier::Variable>(var));

      TBOX_ASSERT(flux_var);
   }

}

/*
 *******************************************************************
 *
 * AMR-unaware patch-centered computational kernels.
 *
 *******************************************************************
 */

void
CellPoissonFACOps::computeFluxOnPatch(
   const hier::Patch& patch,
   const hier::IntVector& ratio_to_coarser_level,
   const pdat::CellData<double>& w_data,
   pdat::SideData<double>& Dgradw_data) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(d_dim, patch, ratio_to_coarser_level,
      w_data, Dgradw_data);
   TBOX_ASSERT(patch.inHierarchy());
   TBOX_ASSERT(w_data.getGhostCellWidth() >=
      hier::IntVector::getOne(ratio_to_coarser_level.getDim()));

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const hier::Box& box = patch.getBox();
   const int* lower = &box.lower()[0];
   const int* upper = &box.upper()[0];
   const double* dx = patch_geom->getDx();

   if (d_poisson_spec.dIsConstant()) {
      double D_value = d_poisson_spec.getDConstant();
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D) (
            Dgradw_data.getPointer(0),
            Dgradw_data.getPointer(1),
            &Dgradw_data.getGhostCellWidth()[0],
            &Dgradw_data.getGhostCellWidth()[1],
            D_value,
            w_data.getPointer(),
            &w_data.getGhostCellWidth()[0],
            &w_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D) (
            Dgradw_data.getPointer(0),
            Dgradw_data.getPointer(1),
            Dgradw_data.getPointer(2),
            &Dgradw_data.getGhostCellWidth()[0],
            &Dgradw_data.getGhostCellWidth()[1],
            &Dgradw_data.getGhostCellWidth()[2],
            D_value,
            w_data.getPointer(),
            &w_data.getGhostCellWidth()[0],
            &w_data.getGhostCellWidth()[1],
            &w_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx);
      }
   } else {
      std::shared_ptr<pdat::SideData<double> > D_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_poisson_spec.getDPatchDataId())));
      TBOX_ASSERT(D_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(compfluxvardc2d, COMPFLUXVARDC2D) (
            Dgradw_data.getPointer(0),
            Dgradw_data.getPointer(1),
            &Dgradw_data.getGhostCellWidth()[0],
            &Dgradw_data.getGhostCellWidth()[1],
            D_data->getPointer(0),
            D_data->getPointer(1),
            &D_data->getGhostCellWidth()[0],
            &D_data->getGhostCellWidth()[1],
            w_data.getPointer(),
            &w_data.getGhostCellWidth()[0],
            &w_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx);
      }
      if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(compfluxvardc3d, COMPFLUXVARDC3D) (
            Dgradw_data.getPointer(0),
            Dgradw_data.getPointer(1),
            Dgradw_data.getPointer(2),
            &Dgradw_data.getGhostCellWidth()[0],
            &Dgradw_data.getGhostCellWidth()[1],
            &Dgradw_data.getGhostCellWidth()[2],
            D_data->getPointer(0),
            D_data->getPointer(1),
            D_data->getPointer(2),
            &D_data->getGhostCellWidth()[0],
            &D_data->getGhostCellWidth()[1],
            &D_data->getGhostCellWidth()[2],
            w_data.getPointer(),
            &w_data.getGhostCellWidth()[0],
            &w_data.getGhostCellWidth()[1],
            &w_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx);
      }
   }

   const int patch_ln = patch.getPatchLevelNumber();

   if (d_cf_discretization == "Ewing" && patch_ln > d_ln_min) {
      ewingFixFlux(patch,
         w_data,
         Dgradw_data,
         ratio_to_coarser_level);
   }

}

void
CellPoissonFACOps::computeResidualOnPatch(
   const hier::Patch& patch,
   const pdat::SideData<double>& flux_data,
   const pdat::CellData<double>& soln_data,
   const pdat::CellData<double>& rhs_data,
   pdat::CellData<double>& residual_data) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY5(d_dim, patch, flux_data, soln_data,
      rhs_data, residual_data);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const hier::Box& box = patch.getBox();
   const int* lower = &box.lower()[0];
   const int* upper = &box.upper()[0];
   const double* dx = patch_geom->getDx();

   double scalar_field_constant;
   if (d_poisson_spec.cIsVariable()) {
      std::shared_ptr<pdat::CellData<double> > scalar_field_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_poisson_spec.getCPatchDataId())));
      TBOX_ASSERT(scalar_field_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(compresvarsca2d, COMPRESVARSCA2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0], &lower[1], &upper[1],
            dx);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(compresvarsca3d, COMPRESVARSCA3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            &residual_data.getGhostCellWidth()[2],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            &scalar_field_data->getGhostCellWidth()[2],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
            dx);
      }
   } else if (d_poisson_spec.cIsConstant()) {
      scalar_field_constant = d_poisson_spec.getCConstant();
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0], &lower[1], &upper[1],
            dx);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            &residual_data.getGhostCellWidth()[2],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
            dx);
      }
   } else {
      scalar_field_constant = 0.0;
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(compresconsca2d, COMPRESCONSCA2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0], &lower[1], &upper[1],
            dx);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(compresconsca3d, COMPRESCONSCA3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            residual_data.getPointer(),
            &residual_data.getGhostCellWidth()[0],
            &residual_data.getGhostCellWidth()[1],
            &residual_data.getGhostCellWidth()[2],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2],
            dx);
      }
   }
}

void
CellPoissonFACOps::redOrBlackSmoothingOnPatch(
   const hier::Patch& patch,
   const pdat::SideData<double>& flux_data,
   const pdat::CellData<double>& rhs_data,
   pdat::CellData<double>& soln_data,
   char red_or_black,
   double* p_maxres) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(d_dim, patch, flux_data, soln_data,
      rhs_data);
   TBOX_ASSERT(red_or_black == 'r' || red_or_black == 'b');

   const int offset = red_or_black == 'r' ? 0 : 1;
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const hier::Box& box = patch.getBox();
   const int* lower = &box.lower()[0];
   const int* upper = &box.upper()[0];
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > scalar_field_data;
   double scalar_field_constant;
   std::shared_ptr<pdat::SideData<double> > diffcoef_data;
   double diffcoef_constant;

   if (d_poisson_spec.cIsVariable()) {
      scalar_field_data = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(d_poisson_spec.getCPatchDataId()));
   } else if (d_poisson_spec.cIsConstant()) {
      scalar_field_constant = d_poisson_spec.getCConstant();
   } else {
      scalar_field_constant = 0.0;
   }
   if (d_poisson_spec.dIsVariable()) {
      diffcoef_data = SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_poisson_spec.getDPatchDataId()));
   } else {
      diffcoef_constant = d_poisson_spec.getDConstant();
   }

   double maxres = 0.0;
   if (d_poisson_spec.dIsVariable() && d_poisson_spec.cIsVariable()) {
      TBOX_ASSERT(scalar_field_data);
      TBOX_ASSERT(diffcoef_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcvarsf2d, RBGSWITHFLUXMAXVARDCVARSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcvarsf3d, RBGSWITHFLUXMAXVARDCVARSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            diffcoef_data->getPointer(2),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            &diffcoef_data->getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            &scalar_field_data->getGhostCellWidth()[2],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   } else if (d_poisson_spec.dIsVariable() && d_poisson_spec.cIsConstant()) {
      TBOX_ASSERT(diffcoef_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf2d, RBGSWITHFLUXMAXVARDCCONSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf3d, RBGSWITHFLUXMAXVARDCCONSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            diffcoef_data->getPointer(2),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            &diffcoef_data->getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   } else if (d_poisson_spec.dIsVariable() && d_poisson_spec.cIsZero()) {
      TBOX_ASSERT(diffcoef_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf2d, RBGSWITHFLUXMAXVARDCCONSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxvardcconsf3d, RBGSWITHFLUXMAXVARDCCONSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_data->getPointer(0),
            diffcoef_data->getPointer(1),
            diffcoef_data->getPointer(2),
            &diffcoef_data->getGhostCellWidth()[0],
            &diffcoef_data->getGhostCellWidth()[1],
            &diffcoef_data->getGhostCellWidth()[2],
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   } else if (!d_poisson_spec.dIsVariable() && d_poisson_spec.cIsVariable()) {
      TBOX_ASSERT(scalar_field_data);
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcvarsf2d, RBGSWITHFLUXMAXCONDCVARSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcvarsf3d, RBGSWITHFLUXMAXCONDCVARSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            scalar_field_data->getPointer(),
            &scalar_field_data->getGhostCellWidth()[0],
            &scalar_field_data->getGhostCellWidth()[1],
            &scalar_field_data->getGhostCellWidth()[2],
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   } else if (!d_poisson_spec.dIsVariable() && d_poisson_spec.cIsConstant()) {
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf2d, RBGSWITHFLUXMAXCONDCCONSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf3d, RBGSWITHFLUXMAXCONDCCONSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            scalar_field_constant,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   } else if (!d_poisson_spec.dIsVariable() && d_poisson_spec.cIsZero()) {
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf2d, RBGSWITHFLUXMAXCONDCCONSF2D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            dx,
            &offset, &maxres);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(rbgswithfluxmaxcondcconsf3d, RBGSWITHFLUXMAXCONDCCONSF3D) (
            flux_data.getPointer(0),
            flux_data.getPointer(1),
            flux_data.getPointer(2),
            &flux_data.getGhostCellWidth()[0],
            &flux_data.getGhostCellWidth()[1],
            &flux_data.getGhostCellWidth()[2],
            diffcoef_constant,
            rhs_data.getPointer(),
            &rhs_data.getGhostCellWidth()[0],
            &rhs_data.getGhostCellWidth()[1],
            &rhs_data.getGhostCellWidth()[2],
            0.0,
            soln_data.getPointer(),
            &soln_data.getGhostCellWidth()[0],
            &soln_data.getGhostCellWidth()[1],
            &soln_data.getGhostCellWidth()[2],
            &lower[0], &upper[0],
            &lower[1], &upper[1],
            &lower[2], &upper[2],
            dx,
            &offset, &maxres);
      }
   }

   *p_maxres = maxres;
}

void
CellPoissonFACOps::xeqScheduleProlongation(
   int dst_id,
   int src_id,
   int scr_id,
   int dest_ln)
{
   if (!d_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }
   xfer::RefineAlgorithm refiner;
   refiner.registerRefine(dst_id,
      src_id,
      scr_id,
      d_prolongation_refine_operator);
   refiner.resetSchedule(d_prolongation_refine_schedules[dest_ln]);
   d_prolongation_refine_schedules[dest_ln]->fillData(0.0);
   d_prolongation_refine_algorithm->resetSchedule(
      d_prolongation_refine_schedules[dest_ln]);
}

void
CellPoissonFACOps::xeqScheduleURestriction(
   int dst_id,
   int src_id,
   int dest_ln)
{
   if (!d_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }

   xfer::CoarsenAlgorithm coarsener(d_dim);
   coarsener.registerCoarsen(dst_id,
      src_id,
      d_urestriction_coarsen_operator);
   coarsener.resetSchedule(d_urestriction_coarsen_schedules[dest_ln]);
   d_urestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_urestriction_coarsen_algorithm->resetSchedule(
      d_urestriction_coarsen_schedules[dest_ln]);
}

void
CellPoissonFACOps::xeqScheduleRRestriction(
   int dst_id,
   int src_id,
   int dest_ln)
{
   if (!d_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }

   xfer::CoarsenAlgorithm coarsener(d_dim);
   coarsener.registerCoarsen(dst_id,
      src_id,
      d_rrestriction_coarsen_operator);
   coarsener.resetSchedule(d_rrestriction_coarsen_schedules[dest_ln]);
   d_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_rrestriction_coarsen_algorithm->resetSchedule(
      d_rrestriction_coarsen_schedules[dest_ln]);
}

void
CellPoissonFACOps::xeqScheduleFluxCoarsen(
   int dst_id,
   int src_id,
   int dest_ln)
{
   if (!d_flux_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }

   xfer::CoarsenAlgorithm coarsener(d_dim);
   coarsener.registerCoarsen(dst_id,
      src_id,
      d_flux_coarsen_operator);

   coarsener.resetSchedule(d_flux_coarsen_schedules[dest_ln]);
   d_flux_coarsen_schedules[dest_ln]->coarsenData();
   d_flux_coarsen_algorithm->resetSchedule(d_flux_coarsen_schedules[dest_ln]);
}

void
CellPoissonFACOps::xeqScheduleGhostFill(
   int dst_id,
   int dest_ln)
{
   if (!d_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }
   xfer::RefineAlgorithm refiner;
   refiner.
   registerRefine(dst_id,
      dst_id,
      dst_id,
      d_ghostfill_refine_operator);
   refiner.
   resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
   d_ghostfill_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_refine_algorithm->resetSchedule(
      d_ghostfill_refine_schedules[dest_ln]);
}

void
CellPoissonFACOps::xeqScheduleGhostFillNoCoarse(
   int dst_id,
   int dest_ln)
{
   if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found." << std::endl);
   }
   xfer::RefineAlgorithm refiner;
   refiner.
   registerRefine(dst_id,
      dst_id,
      dst_id,
      d_ghostfill_nocoarse_refine_operator);
   refiner.
   resetSchedule(d_ghostfill_nocoarse_refine_schedules[dest_ln]);
   d_ghostfill_nocoarse_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_nocoarse_refine_algorithm->resetSchedule(
      d_ghostfill_nocoarse_refine_schedules[dest_ln]);
}

void
CellPoissonFACOps::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_cell_scratch_var[d].reset();
      s_flux_scratch_var[d].reset();
      s_oflux_scratch_var[d].reset();
   }
}

}
}
