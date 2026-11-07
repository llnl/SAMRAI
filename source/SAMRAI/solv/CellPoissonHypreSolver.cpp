/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Hypre solver interface for diffusion-like elliptic problems.
 *
 ************************************************************************/
#include "SAMRAI/solv/CellPoissonHypreSolver.h"

#ifdef HAVE_HYPRE

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/math/ArrayDataBasicOps.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <cstdlib>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

void SAMRAI_F77_FUNC(compdiagvariablec2d, COMPDIAGVARIABLEC2D) (
   double* diag,
   const double* c,
   const double* offdiagi,
   const double* offdiagj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(compdiagscalarc2d, COMPDIAGSCALARC2D) (
   double* diag,
   const double* c,
   const double* offdiagi,
   const double* offdiagj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(compdiagzeroc2d, COMPDIAGZEROC2D) (
   double* diag,
   const double* offdiagi,
   const double* offdiagj,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(adjbdry2d, ADJBDRY2D) (
   double* diag,
   const double* offdiagi,
   const double* offdiagj,
   const int* pifirst, const int* pilast,
   const int* pjfirst, const int* pjlast,
   const double* acoef,
   const double* bcoef,
   const int* aifirst, const int* ailast,
   const int* ajfirst, const int* ajlast,
   const double* Ak0,
   const int* kifirst, const int* kilast,
   const int* kjfirst, const int* kjlast,
   const int* lower, const int* upper,
   const int* location,
   const double* h);
void SAMRAI_F77_FUNC(adjbdryconstoffdiags2d, ADJBDRYCONSTOFFDIAGS2D) (
   double* diag,
   const double* offdiag,
   const int* pifirst,
   const int* pilast,
   const int* pjfirst,
   const int* pjlast,
   const double* acoef,
   const int* aifirst,
   const int* ailast,
   const int* ajfirst,
   const int* ajlast,
   const double* Ak0,
   const int* kifirst,
   const int* kilast,
   const int* kjfirst,
   const int* kjlast,
   const int* lower, const int* upper,
   const int* location,
   const double* h);
void SAMRAI_F77_FUNC(adjustrhs2d, ADJUSTRHS2D) (double* rhs,
   const int* rifirst,
   const int* rilast,
   const int* rjfirst,
   const int* rjlast,
   const double* Ak0,
   const int* kifirst,
   const int* kilast,
   const int* kjfirst,
   const int* kjlast,
   const double* gcoef,
   const int* aifirst,
   const int* ailast,
   const int* ajfirst,
   const int* ajlast,
   const int* lower, const int* upper,
   const int* location);

void SAMRAI_F77_FUNC(compdiagvariablec3d, COMPDIAGVARIABLEC3D) (
   double* diag,
   const double* c,
   const double* offdiagi,
   const double* offdiagj,
   const double* offdiagk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(compdiagscalarc3d, COMPDIAGSCALARC3D) (
   double* diag,
   const double* c,
   const double* offdiagi,
   const double* offdiagj,
   const double* offdiagk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(compdiagzeroc3d, COMPDIAGZEROC3D) (
   double* diag,
   const double* offdiagi,
   const double* offdiagj,
   const double* offdiagk,
   const int* ifirst,
   const int* ilast,
   const int* jfirst,
   const int* jlast,
   const int* kfirst,
   const int* klast,
   const double* cscale,
   const double* dscale);
void SAMRAI_F77_FUNC(adjbdry3d, ADJBDRY3D) (
   double* diag,
   const double* offdiagi,
   const double* offdiagj,
   const double* offdiagk,
   const int* pifirst,
   const int* pilast,
   const int* pjfirst,
   const int* pjlast,
   const int* pkfirst,
   const int* pklast,
   const double* acoef,
   const double* bcoef,
   const int* aifirst,
   const int* ailast,
   const int* ajfirst,
   const int* ajlast,
   const int* akfirst,
   const int* aklast,
   const double* Ak0,
   const int* kifirst,
   const int* kilast,
   const int* kjfirst,
   const int* kjlast,
   const int* kkfirst,
   const int* kklast,
   const int* lower, const int* upper,
   const int* location,
   const double* h);
void SAMRAI_F77_FUNC(adjbdryconstoffdiags3d, ADJBDRYCONSTOFFDIAGS3D) (
   double* diag,
   const double* offdiag,
   const int* pifirst,
   const int* pilast,
   const int* pjfirst,
   const int* pjlast,
   const int* pkfirst,
   const int* pklast,
   const double* acoef,
   const int* aifirst,
   const int* ailast,
   const int* ajfirst,
   const int* ajlast,
   const int* akfirst,
   const int* aklast,
   const double* Ak0,
   const int* kifirst,
   const int* kilast,
   const int* kjfirst,
   const int* kjlast,
   const int* kkfirst,
   const int* kklast,
   const int* lower, const int* upper,
   const int* location,
   const double* h);
void SAMRAI_F77_FUNC(adjustrhs3d, ADJUSTRHS3D) (double* rhs,
   const int* rifirst,
   const int* rilast,
   const int* rjfirst,
   const int* rjlast,
   const int* rkfirst,
   const int* rklast,
   const double* Ak0,
   const int* kifirst,
   const int* kilast,
   const int* kjfirst,
   const int* kjlast,
   const int* kkfirst,
   const int* kklast,
   const double* gcoef,
   const int* aifirst,
   const int* ailast,
   const int* ajfirst,
   const int* ajlast,
   const int* akfirst,
   const int* aklast,
   const int* lower, const int* upper,
   const int* location);

}

namespace SAMRAI {
namespace solv {

std::shared_ptr<pdat::OutersideVariable<double> >
CellPoissonHypreSolver::s_Ak0_var[SAMRAI::MAX_DIM_VAL];

tbox::StartupShutdownManager::Handler CellPoissonHypreSolver::s_finalize_handler(
   0,
   0,
   0,
   CellPoissonHypreSolver::finalizeCallback,
   tbox::StartupShutdownManager::priorityVariables);

/*
 *************************************************************************
 * Constructor
 *************************************************************************
 */

CellPoissonHypreSolver::CellPoissonHypreSolver(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(object_name),
   d_ln(-1),
   d_context(hier::VariableDatabase::getDatabase()->getContext(
                object_name + "::context")),
   d_physical_bc_coef_strategy(&d_physical_bc_simple_case),
   d_physical_bc_simple_case(dim, d_object_name + "::simple bc"),
   d_cf_bc_coef(dim, object_name + "::coarse-fine bc coefs"),
   d_Ak0_id(-1),
   d_soln_depth(0),
   d_rhs_depth(0),
   d_max_iterations(10),
   d_relative_residual_tol(1e-10),
   d_number_iterations(-1),
   d_num_pre_relax_steps(1),
   d_num_post_relax_steps(1),
   d_relative_residual_norm(-1.0),
   d_use_smg(true),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_grid(0),
   d_stencil(0),
   d_matrix(0),
   d_linear_rhs(0),
   d_linear_sol(0),
   d_mg_data(0),
   d_print_solver_info(false)
{
   if (d_dim == tbox::Dimension(1) || d_dim > tbox::Dimension(3)) {
      TBOX_ERROR(" CellPoissonHypreSolver : DIM == 1 or > 3 not implemented");
   }

   t_solve_system = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonHypreSolver::solveSystem()");
   t_set_matrix_coefficients = tbox::TimerManager::getManager()->
      getTimer("solv::CellPoissonHypreSolver::setMatrixCoefficients()");

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
   if (!s_Ak0_var[d_dim.getValue() - 1]) {
      s_Ak0_var[d_dim.getValue() - 1].reset(
         new pdat::OutersideVariable<double>(d_dim, d_object_name + "::Ak0",
                                             d_allocator,
                                             1));
   }

   d_Ak0_id =
      vdb->registerVariableAndContext(s_Ak0_var[d_dim.getValue() - 1],
         d_context,
         hier::IntVector::getZero(d_dim));

   getFromInput(input_db);
}

/*
 ********************************************************************
 * Set state from database
 ********************************************************************
 */

void
CellPoissonHypreSolver::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{
   if (input_db) {
      d_print_solver_info =
         input_db->getBoolWithDefault("print_solver_info", false);

      d_max_iterations = input_db->getIntegerWithDefault("max_iterations", 10);
      if (!(d_max_iterations >= 1)) {
         INPUT_RANGE_ERROR("max_iterations");
      }

      d_relative_residual_tol =
         input_db->getDoubleWithDefault("relative_residual_tol", 1.0e-10);
      if (!(d_relative_residual_tol > 0)) {
         INPUT_RANGE_ERROR("relative_residual_tol");
      }

      if (input_db->isDouble("residual_tol")) {
         TBOX_ERROR("CellPoissonHypreSolver getFromInput error...\n"
            << "The parameter 'residual_tol' has been replaced\n"
            << "by 'relative_residual_tol' to be more descriptive.\n"
            << "Please change the parameter name in the input database.");
      }

      d_num_pre_relax_steps =
         input_db->getIntegerWithDefault("num_pre_relax_steps", 1);
      if (!(d_num_pre_relax_steps >= 0)) {
         INPUT_RANGE_ERROR("num_pre_relax_steps");
      }

      d_num_post_relax_steps =
         input_db->getIntegerWithDefault("num_post_relax_steps", 1);
      if (!(d_num_post_relax_steps >= 0)) {
         INPUT_RANGE_ERROR("num_post_relax_steps");
      }

      d_use_smg = input_db->getBoolWithDefault("use_smg", true);
   }
}

/*
 ********************************************************************
 * Initialize internal data for a given hierarchy level
 * After setting internal data, propagate the information
 * to the major algorithm objects.  Allocate data for
 * storing boundary condition-dependent quantities for
 * adding to souce term before solving.
 ********************************************************************
 */

void
CellPoissonHypreSolver::initializeSolverState(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int ln)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);

   deallocateSolverState();

   d_hierarchy = hierarchy;
   d_ln = ln;

   hier::IntVector max_gcw(d_dim, 1);
   d_cf_boundary.reset(
      new hier::CoarseFineBoundary(*d_hierarchy, d_ln, max_gcw));

   d_physical_bc_simple_case.setHierarchy(d_hierarchy, d_ln, d_ln);

   d_number_iterations = -1;
   d_relative_residual_norm = -1.0;

   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(d_ln));
   level->allocatePatchData(d_Ak0_id);
   allocateHypreData();
}

/*
 ********************************************************************
 * Deallocate data initialized by initializeSolverState
 ********************************************************************
 */

void
CellPoissonHypreSolver::deallocateSolverState()
{
   if (!d_hierarchy) {
      return;
   }

   d_cf_boundary->clear();
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(d_ln));
   level->deallocatePatchData(d_Ak0_id);
   deallocateHypreData();
   d_hierarchy.reset();
   d_ln = -1;
}

/*
 *************************************************************************
 *
 * Allocate the HYPRE data structures that depend only on the level
 * and will not change (grid, stencil, matrix, and vectors).
 *
 *************************************************************************
 */
void
CellPoissonHypreSolver::allocateHypreData()
{
   tbox::SAMRAI_MPI::Comm communicator =
      d_hierarchy->getMPI().getCommunicator();

   /*
    * Set up the grid data - only set grid data for local boxes
    */

   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(d_ln));
   std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
         d_hierarchy->getGridGeometry()));

   TBOX_ASSERT(grid_geometry);

   const hier::IntVector ratio = level->getRatioToLevelZero();
   hier::IntVector periodic_shift =
      grid_geometry->getPeriodicShift(ratio);

   int periodic_flag[SAMRAI::MAX_DIM_VAL];
   tbox::Dimension::dir_t d;
   bool is_periodic = false;
   for (d = 0; d < d_dim.getValue(); ++d) {
      periodic_flag[d] = periodic_shift[d] != 0;
      is_periodic = is_periodic || periodic_flag[d];
   }

   HYPRE_StructGridCreate(communicator, d_dim.getValue(), &d_grid);
   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const hier::Box& box = (*p)->getBox();
      hier::Index lower = box.lower();
      hier::Index upper = box.upper();
      HYPRE_StructGridSetExtents(d_grid, &lower[0], &upper[0]);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (is_periodic) {
      const hier::BoxContainer& level_domain =
         level->getPhysicalDomain(hier::BlockId::zero());
      hier::Box domain_bound(level_domain.front());
      for (hier::BoxContainer::const_iterator i = level_domain.begin();
           i != level_domain.end(); ++i) {
         domain_bound.setLower(
            hier::Index::min(domain_bound.lower(), i->lower()));
         domain_bound.setUpper(
            hier::Index::min(domain_bound.upper(), i->upper()));
      }
      for (d = 0; d < d_dim.getValue(); ++d) {
         if (periodic_flag[d] == true) {
            int tmpi = 1;
            unsigned int p_of_two;
            for (p_of_two = 0; p_of_two < 8 * sizeof(p_of_two) - 1;
                 ++p_of_two) {
               if (tmpi == domain_bound.numberCells(d)) {
                  break;
               }
               if (tmpi > domain_bound.numberCells(d)) {
                  TBOX_ERROR(
                     d_object_name << ": Hypre currently requires\n"
                                   << "that grid size in periodic directions be\n"
                                   << "powers of two.  (This requirement may go\n"
                                   << "away in future versions of hypre.)\n"
                                   << "Size problem in direction "
                                   << d << "\n"
                                   << "Domain bound is "
                                   << domain_bound << ",\n"
                                   << "Size of "
                                   << domain_bound.numberCells() << "\n");
               }
               tmpi = tmpi ? tmpi << 1 : 1;
            }
         }
      }
   }
#endif

   HYPRE_StructGridSetPeriodic(d_grid, &periodic_shift[0]);
   HYPRE_StructGridAssemble(d_grid);

   {
      /*
       * Allocate stencil data and set stencil offsets
       */

      if (d_dim == tbox::Dimension(1)) {
         const int stencil_size = 2;
         int stencil_offsets[2][1] = {
            { -1 }, { 0 }
         };
         HYPRE_StructStencilCreate(d_dim.getValue(), stencil_size, &d_stencil);
         for (int s = 0; s < stencil_size; ++s) {
            HYPRE_StructStencilSetElement(d_stencil, s,
               stencil_offsets[s]);
         }
      } else if (d_dim == tbox::Dimension(2)) {
         const int stencil_size = 3;
         int stencil_offsets[3][2] = {
            { -1, 0 }, { 0, -1 }, { 0, 0 }
         };
         HYPRE_StructStencilCreate(d_dim.getValue(), stencil_size, &d_stencil);
         for (int s = 0; s < stencil_size; ++s) {
            HYPRE_StructStencilSetElement(d_stencil, s,
               stencil_offsets[s]);
         }
      } else if (d_dim == tbox::Dimension(3)) {
         const int stencil_size = 4;
         int stencil_offsets[4][3] = {
            { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }, { 0, 0, 0 }
         };
         HYPRE_StructStencilCreate(d_dim.getValue(), stencil_size, &d_stencil);
         for (int s = 0; s < stencil_size; ++s) {
            HYPRE_StructStencilSetElement(d_stencil, s,
               stencil_offsets[s]);
         }
      }
   }

   {
      int full_ghosts1[2 * 3] = { 1, 1, 0, 0, 0, 0 };
      int no_ghosts1[2 * 3] = { 0, 0, 0, 0, 0, 0 };

      int full_ghosts2[2 * 3] = { 1, 1, 1, 1, 0, 0 };
      int no_ghosts2[2 * 3] = { 0, 0, 0, 0, 0, 0 };

      int full_ghosts3[2 * 3] = { 1, 1, 1, 1, 1, 1 };
      int no_ghosts3[2 * 3] = { 0, 0, 0, 0, 0, 0 };

      /*
       * Allocate the structured matrix
       */

      int* full_ghosts = 0;
      int* no_ghosts = 0;

      if (d_dim == tbox::Dimension(1)) {
         full_ghosts = full_ghosts1;
         no_ghosts = no_ghosts1;
      } else if (d_dim == tbox::Dimension(2)) {
         full_ghosts = full_ghosts2;
         no_ghosts = no_ghosts2;
      } else if (d_dim == tbox::Dimension(3)) {
         full_ghosts = full_ghosts3;
         no_ghosts = no_ghosts3;
      } else {
         TBOX_ERROR(
            "CellPoissonHypreSolver does not yet support dimension " << d_dim);
      }

      HYPRE_StructMatrixCreate(communicator,
         d_grid,
         d_stencil,
         &d_matrix);
      HYPRE_StructMatrixSetNumGhost(d_matrix, full_ghosts);
      HYPRE_StructMatrixSetSymmetric(d_matrix, 1);
      HYPRE_StructMatrixInitialize(d_matrix);

      HYPRE_StructVectorCreate(communicator,
         d_grid,
         &d_linear_rhs);
      HYPRE_StructVectorSetNumGhost(d_linear_rhs, no_ghosts);
      HYPRE_StructVectorInitialize(d_linear_rhs);

      HYPRE_StructVectorCreate(communicator,
         d_grid,
         &d_linear_sol);
      HYPRE_StructVectorSetNumGhost(d_linear_sol, full_ghosts);
      HYPRE_StructVectorInitialize(d_linear_sol);
   }
}

/*
 *************************************************************************
 *
 * The destructor deallocates solver data.
 *
 *************************************************************************
 */

CellPoissonHypreSolver::~CellPoissonHypreSolver()
{
   deallocateHypreData();

   if (d_hierarchy) {
      d_hierarchy->getPatchLevel(0)->deallocatePatchData(d_Ak0_id);
   }
   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();
   vdb->removePatchDataIndex(d_Ak0_id);
}

/*
 *************************************************************************
 *
 * Deallocate HYPRE data and solver.  HYPRE requires that we
 * check whether HYPRE has already deallocated this data.
 * Note that the HYPRE solver, d_mg_data, was created at
 * the end of setMatrixCoefficients.
 *
 *************************************************************************
 */

void
CellPoissonHypreSolver::deallocateHypreData()
{
   if (d_stencil) {
      HYPRE_StructStencilDestroy(d_stencil);
      d_stencil = 0;
   }
   if (d_grid) {
      HYPRE_StructGridDestroy(d_grid);
      d_grid = 0;
   }
   if (d_matrix) {
      HYPRE_StructMatrixDestroy(d_matrix);
      d_matrix = 0;
   }
   if (d_linear_rhs) {
      HYPRE_StructVectorDestroy(d_linear_rhs);
      d_linear_rhs = 0;
   }
   if (d_linear_sol) {
      HYPRE_StructVectorDestroy(d_linear_sol);
      d_linear_sol = 0;
   }
   destroyHypreSolver();
}

/*
 *************************************************************************
 *
 * Copy data into the HYPRE vector structures.
 *
 *************************************************************************
 */

void
CellPoissonHypreSolver::copyToHypre(
   HYPRE_StructVector vector,
   pdat::CellData<double>& src,
   int depth,
   const hier::Box& box)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, src, box);

   hier::Index lower(box.lower());
   hier::Index upper(box.upper());

   if (src.getGhostBox().isSpatiallyEqual(box)) {
      HYPRE_StructVectorSetBoxValues(
         vector, &lower[0], &upper[0], src.getPointer(depth)); 
   } else {
      pdat::CellData<double> tmp(box, 1, hier::IntVector::getZero(d_dim),
                                 d_allocator);

      tmp.copyDepth(0, src, depth);
      HYPRE_StructVectorSetBoxValues(
         vector, &lower[0], &upper[0], tmp.getPointer());
   } 
}

/*
 *************************************************************************
 *
 * Copy data out of the HYPRE vector structures.
 *
 *************************************************************************
 */

void
CellPoissonHypreSolver::copyFromHypre(
   pdat::CellData<double>& dst,
   int depth,
   HYPRE_StructVector vector,
   const hier::Box box)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, dst, box);

   hier::Index lower(box.lower());
   hier::Index upper(box.upper());
   if (dst.getGhostBox().isSpatiallyEqual(box)) {
      HYPRE_StructVectorGetBoxValues(
         vector, &lower[0], &upper[0], dst.getPointer(depth));
   } else {
      pdat::CellData<double> tmp(box, 1, hier::IntVector::getZero(d_dim),
                                 d_allocator);

      HYPRE_StructVectorGetBoxValues(
         vector, &lower[0], &upper[0], tmp.getPointer());
      dst.copyDepth(depth,tmp,0);
   }
}

/*
 *************************************************************************
 *
 * Set the matrix coefficients for the linear system.
 * The matrix coefficients are dependent on the problem
 * specification described by the PoissonSpecificiations
 * object and by the boundary condition.
 *
 *************************************************************************
 */

void
CellPoissonHypreSolver::setMatrixCoefficients(
   const PoissonSpecifications& spec)
{
   if (d_physical_bc_coef_strategy == 0) {
      TBOX_ERROR(
         d_object_name << ": No BC coefficient strategy object!\n"
                       << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                       << "to specify the boundary conidition.  Do it before\n"
                       << "calling setMatrixCoefficients.");
   }

   t_set_matrix_coefficients->start();

   std::shared_ptr<pdat::CellData<double> > C_data;
   std::shared_ptr<pdat::SideData<double> > D_data;

   /*
    * Some computations can be done using high-level math objects.
    * Define the math objects.
    */
   math::ArrayDataBasicOps<double> array_math;
   math::PatchSideDataBasicOps<double> patch_side_math;

   /*
    * The value of the ghost cell based on the Robin boundary condition
    * can be written as the sum of a constant, k0, plus a multiple of the
    * internal cell value, k1*ui.  k1*ui depends on the value of u so it
    * contributes to the product Au,
    * while the constant k0 contributes the right hand side f.
    * We save Ak0 = A*k0(a) to add to f when solving.
    * We assume unit g here because we will multiply it in just before
    * solving, thus allowing everything that does not affect A to change
    * from solve to solve.
    */
   std::shared_ptr<pdat::OutersideData<double> > Ak0;

   /*
    * Loop over patches and set matrix entries for each patch.
    */
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(d_ln));
   const hier::IntVector no_ghosts(d_dim, 0);
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {

      hier::Patch& patch = **pi;

      std::shared_ptr<geom::CartesianPatchGeometry> pg(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));

      TBOX_ASSERT(pg);

      const double* h = pg->getDx();

      const hier::Box patch_box = patch.getBox();

      if (!spec.cIsZero() && !spec.cIsConstant()) {
         C_data = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData(spec.getCPatchDataId()));
         TBOX_ASSERT(C_data);
      }

      if (!spec.dIsConstant()) {
         D_data = SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch.getPatchData(spec.getDPatchDataId()));
         TBOX_ASSERT(D_data);
      }

      Ak0 = SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
            patch.getPatchData(d_Ak0_id));
      TBOX_ASSERT(Ak0);

      Ak0->fillAll(0.0);

      pdat::CellData<double> diagonal(patch_box, 1, no_ghosts,
                                      d_allocator);

      /*
       * Set diagonals to zero so we can accumulate to it.
       * Accumulation is used at boundaries to shift weights
       * for ghost cells onto the diagonal.
       */
      diagonal.fillAll(0.0);

      /*
       * Storage for off-diagonal entries,
       * which can be variable or constant.
       */
      pdat::SideData<double> off_diagonal(patch_box, 1, no_ghosts,
                                          d_allocator);
      /*
       * Compute all off-diagonal entries with no regard to BCs.
       * These off-diagonal entries are simply D/(h*h), according
       * to our central difference formula.
       */
      if (spec.dIsConstant()) {
         for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
            double dhh = spec.getDConstant() / (h[i] * h[i]);
            pdat::ArrayData<double>& off_diag_array(off_diagonal.getArrayData(i));
            off_diag_array.fill(dhh);
         }
      } else {
         for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
            hier::Box sbox(patch_box);
            sbox.growUpper(i, 1);
            array_math.scale(off_diagonal.getArrayData(i),
               1.0 / (h[i] * h[i]),
               D_data->getArrayData(i),
               sbox);
         }
      }

      /*
       * Compute diagonal entries using off-diagonal contributions.
       */
      if (spec.cIsZero()) {
         computeDiagonalEntries(diagonal,
            off_diagonal,
            patch_box);
      } else if (spec.cIsConstant()) {
         computeDiagonalEntries(diagonal,
            spec.getCConstant(),
            off_diagonal,
            patch_box);
      } else {
         computeDiagonalEntries(diagonal,
            *C_data,
            off_diagonal,
            patch_box);
      }

      /*
       * Walk physical domain boundaries and adjust off-diagonals
       * before computation of diagonal entries.
       * The exterior cell's value is
       * uo = ( h*gamma + ui*(beta-h*alpha/2) )/( beta+h*alpha/2 )
       *   = k0 + k1*ui
       * where k0 = h*gamma/( beta+h*alpha/2 )
       * k1 = ( beta-h*alpha/2 )/( beta+h*alpha/2 )
       * Split coupling between interior-exterior cells
       * into two parts: interior-interior coupling (k1)
       * and rhs contribution (k0).
       */
      {
         const std::vector<hier::BoundaryBox>& surface_boxes =
            pg->getCodimensionBoundaries(1);
         const int n_bdry_boxes = static_cast<int>(surface_boxes.size());
         for (int n = 0; n < n_bdry_boxes; ++n) {

            const hier::BoundaryBox& boundary_box = surface_boxes[n];
            if (boundary_box.getBoundaryType() != 1) {
               TBOX_ERROR(
                  d_object_name << ": Illegal boundary type in "
                                << "CellPoissonHypreSolver::setMatrixCoefficients\n");
            }
            const hier::BoundaryBoxUtils bbu(boundary_box);
            const int location_index = boundary_box.getLocationIndex();
            const hier::BoundaryBox trimmed_boundary_box =
               bbu.trimBoundaryBox(patch.getBox());
            const hier::Box bccoef_box =
               bbu.getSurfaceBoxFromBoundaryBox();
            std::shared_ptr<pdat::ArrayData<double> > acoef_data(
               std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1,
                                                          d_allocator));
            std::shared_ptr<pdat::ArrayData<double> > bcoef_data(
               std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1,
                                                          d_allocator));
            std::shared_ptr<pdat::ArrayData<double> > gcoef_data;
            static const double fill_time = 0.0;
            d_physical_bc_coef_strategy->setBcCoefs(acoef_data,
               bcoef_data,
               gcoef_data,
               d_physical_bc_variable,
               patch,
               boundary_box,
               fill_time);
            pdat::ArrayData<double>& Ak0_data =
               Ak0->getArrayData(location_index / 2,
                  location_index % 2);
            adjustBoundaryEntries(diagonal,
               off_diagonal,
               patch_box,
               *acoef_data,
               *bcoef_data,
               bccoef_box,
               Ak0_data,
               trimmed_boundary_box,
               h);
         }
      }
      /*
       * Walk coarse-fine boundaries and adjust off-diagonals
       * according data in ghost cells.
       */
      if (d_ln > 0) {
         /*
          * There are potentially coarse-fine boundaries to deal with.
          */

         std::vector<hier::BoundaryBox> empty_vector(0,
                                                     hier::BoundaryBox(d_dim));
         const std::vector<hier::BoundaryBox>& surface_boxes =
            d_dim == tbox::Dimension(2) ? d_cf_boundary->getEdgeBoundaries(pi->getGlobalId()) :
            (d_dim ==
             tbox::Dimension(3) ? d_cf_boundary->getFaceBoundaries(pi->getGlobalId()) :
             empty_vector);

         const int n_bdry_boxes = static_cast<int>(surface_boxes.size());

         for (int n = 0; n < n_bdry_boxes; ++n) {

            const hier::BoundaryBox& boundary_box = surface_boxes[n];
            if (boundary_box.getBoundaryType() != 1) {
               TBOX_ERROR(
                  d_object_name << ": Illegal boundary type in "
                                << "CellPoissonHypreSolver::setMatrixCoefficients\n");
            }
            const int location_index = boundary_box.getLocationIndex();
            const hier::BoundaryBoxUtils bbu(boundary_box);
            const hier::BoundaryBox trimmed_boundary_box =
               bbu.trimBoundaryBox(patch.getBox());
            const hier::Box bccoef_box =
               bbu.getSurfaceBoxFromBoundaryBox();
            std::shared_ptr<pdat::ArrayData<double> > acoef_data(
               std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1,
                                                          d_allocator));
            std::shared_ptr<pdat::ArrayData<double> > bcoef_data(
               std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1,
                                                          d_allocator));
            std::shared_ptr<pdat::ArrayData<double> > gcoef_data;
            static const double fill_time = 0.0;
            /*
             * Reset invalid ghost data id to help detect use in setBcCoefs.
             */
            d_cf_bc_coef.setGhostDataId(-1, hier::IntVector::getZero(d_dim));
            d_cf_bc_coef.setBcCoefs(acoef_data,
               bcoef_data,
               gcoef_data,
               d_coarsefine_bc_variable,
               patch,
               boundary_box,
               fill_time);
            pdat::ArrayData<double>& Ak0_data =
               Ak0->getArrayData(location_index / 2,
                  location_index % 2);
            adjustBoundaryEntries(diagonal,
               off_diagonal,
               patch_box,
               *acoef_data,
               *bcoef_data,
               bccoef_box,
               Ak0_data,
               trimmed_boundary_box,
               h);
         }
      }

      /*
       * Copy matrix entries to HYPRE matrix structure.  Note that
       * we translate our temporary diagonal/off-diagonal storage into the
       * HYPRE symmetric storage scheme for the stencil specified earlier.
       */
      const int stencil_size = d_dim.getValue() + 1;
      int stencil_indices[stencil_size];
      double* mat_entries =
         new double[stencil_size*patch_box.size()];

      for (int i = 0; i < stencil_size; ++i) stencil_indices[i] = i;

      unsigned int offset=0;
      unsigned int count=0;
      unsigned int offsetx=0;
      double* offdx=off_diagonal.getPointer(0);
      unsigned int offsety=0;
      double* offdy= (d_dim > tbox::Dimension(1)) ?
                     off_diagonal.getPointer(1) : 0;
      unsigned int offsetz=0;
      double* offdz= (d_dim > tbox::Dimension(2)) ?
                     off_diagonal.getPointer(2) : 0;
      double* diag=diagonal.getPointer();

      hier::Index lower(patch_box.lower());
      hier::Index upper(patch_box.upper());
      int ilower[3]={0,0,0};
      int iupper[3]={0,0,0};
      for(int d=0;d<d_dim.getValue();d++){
         ilower[d]=lower[d];
         iupper[d]=upper[d];
      }
      for(int k=ilower[2];k<=iupper[2];k++){
         for(int j=ilower[1];j<=iupper[1];j++){
            for(int i=ilower[0];i<=iupper[0];i++){
               mat_entries[offset++] = offdx[offsetx++];
               if (d_dim > tbox::Dimension(1))
                  mat_entries[offset++] = offdy[offsety++];
               if (d_dim > tbox::Dimension(2))
                  mat_entries[offset++] = offdz[offsetz++];
               mat_entries[offset++] = diag[count++];
            }
            offsetx++;
         }
         if (d_dim > tbox::Dimension(2))
            offsety+=patch_box.numberCells(0);
      }
      HYPRE_StructMatrixSetBoxValues(d_matrix, &ilower[0], &iupper[0],
                                     stencil_size, stencil_indices,
                                     mat_entries);
      delete[] mat_entries;
   } // end patch loop

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat_bA.out", d_matrix, 1);
   }

   HYPRE_StructMatrixAssemble(d_matrix);

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat_aA.out", d_matrix, 1);
   }

   t_set_matrix_coefficients->stop();

   setupHypreSolver();
}

/*
 **********************************************************************
 * Add g*A*k0(a) from physical boundaries to rhs.
 * This operation is done for physical as well as cf boundaries,
 * so it is placed in a function.
 **********************************************************************
 */

void
CellPoissonHypreSolver::add_gAk0_toRhs(
   const hier::Patch& patch,
   const std::vector<hier::BoundaryBox>& bdry_boxes,
   const RobinBcCoefStrategy* robin_bc_coef,
   pdat::CellData<double>& rhs)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, patch, rhs);

   /*
    * g*A*k0(a) is the storage for adjustments to be made to the rhs
    * when we solve. This is the value of the weight of the ghost cell
    * value for the interior cell, times k0.  It is independent of u,
    * and so is moved to the rhs.  Before solving, g*A*k0(a) is added
    * to rhs.
    */
   std::shared_ptr<pdat::OutersideData<double> > Ak0(
      SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
         patch.getPatchData(d_Ak0_id)));

   TBOX_ASSERT(Ak0);

   const int n_bdry_boxes = static_cast<int>(bdry_boxes.size());
   for (int n = 0; n < n_bdry_boxes; ++n) {

      const hier::BoundaryBox& boundary_box = bdry_boxes[n];
#ifdef DEBUG_CHECK_ASSERTIONS
      if (boundary_box.getBoundaryType() != 1) {
         TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                                  << "CellPoissonHypreSolver::add_gAk0_toRhs\n");
      }
#endif
      const int location_index = boundary_box.getLocationIndex();
      const hier::BoundaryBoxUtils bbu(boundary_box);
      const hier::BoundaryBox trimmed_boundary_box =
         bbu.trimBoundaryBox(patch.getBox());
      const hier::Index& lower = trimmed_boundary_box.getBox().lower();
      const hier::Index& upper = trimmed_boundary_box.getBox().upper();
      const hier::Box& rhsbox = rhs.getArrayData().getBox();
      const hier::Box& Ak0box = Ak0->getArrayData(location_index / 2,
            location_index % 2).getBox();
      const hier::Box bccoef_box = bbu.getSurfaceBoxFromBoundaryBox();
      std::shared_ptr<pdat::ArrayData<double> > acoef_data;
      std::shared_ptr<pdat::ArrayData<double> > bcoef_data;
      std::shared_ptr<pdat::ArrayData<double> > gcoef_data(
         std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1,
                                                    d_allocator));
      static const double fill_time = 0.0;
      robin_bc_coef->setBcCoefs(acoef_data,
         bcoef_data,
         gcoef_data,
         d_physical_bc_variable,
         patch,
         boundary_box,
         fill_time);
      /*
       * Nomenclature for indices: cel=first-cell, gho=ghost,
       * beg=beginning, end=ending.
       */
      if (d_dim == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(adjustrhs2d, ADJUSTRHS2D) (rhs.getPointer(d_rhs_depth),
            &rhsbox.lower()[0],
            &rhsbox.upper()[0],
            &rhsbox.lower()[1],
            &rhsbox.upper()[1],
            Ak0->getPointer(location_index / 2, location_index % 2),
            &Ak0box.lower()[0],
            &Ak0box.upper()[0],
            &Ak0box.lower()[1],
            &Ak0box.upper()[1],
            gcoef_data->getPointer(),
            &bccoef_box.lower()[0],
            &bccoef_box.upper()[0],
            &bccoef_box.lower()[1],
            &bccoef_box.upper()[1],
            &lower[0], &upper[0],
            &location_index);
      } else if (d_dim == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(adjustrhs3d, ADJUSTRHS3D) (rhs.getPointer(d_rhs_depth),
            &rhsbox.lower()[0],
            &rhsbox.upper()[0],
            &rhsbox.lower()[1],
            &rhsbox.upper()[1],
            &rhsbox.lower()[2],
            &rhsbox.upper()[2],
            Ak0->getPointer(location_index / 2, location_index % 2),
            &Ak0box.lower()[0],
            &Ak0box.upper()[0],
            &Ak0box.lower()[1],
            &Ak0box.upper()[1],
            &Ak0box.lower()[2],
            &Ak0box.upper()[2],
            gcoef_data->getPointer(),
            &bccoef_box.lower()[0],
            &bccoef_box.upper()[0],
            &bccoef_box.lower()[1],
            &bccoef_box.upper()[1],
            &bccoef_box.lower()[2],
            &bccoef_box.upper()[2],
            &lower[0], &upper[0],
            &location_index);
      }
   }
}

/*
 *************************************************************************
 * Create the hypre solver and set it according to the current state.
 *************************************************************************
 */
void
CellPoissonHypreSolver::setupHypreSolver()
{
   TBOX_ASSERT(d_mg_data == 0);

   tbox::SAMRAI_MPI::Comm communicator =
      d_hierarchy->getMPI().getCommunicator();

   if (d_use_smg) {
      HYPRE_StructSMGCreate(communicator, &d_mg_data);
      HYPRE_StructSMGSetMemoryUse(d_mg_data, 0);
      HYPRE_StructSMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructSMGSetTol(d_mg_data, d_relative_residual_tol);
      HYPRE_StructSMGSetLogging(d_mg_data, 1);
      HYPRE_StructSMGSetNumPreRelax(d_mg_data,
         d_num_pre_relax_steps);
      HYPRE_StructSMGSetNumPostRelax(d_mg_data,
         d_num_post_relax_steps);
      HYPRE_StructSMGSetup(d_mg_data,
         d_matrix,
         d_linear_rhs,
         d_linear_sol);
   } else {
      HYPRE_StructPFMGCreate(communicator, &d_mg_data);
      HYPRE_StructPFMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructPFMGSetTol(d_mg_data, d_relative_residual_tol);
      HYPRE_StructPFMGSetLogging(d_mg_data, 1);
      HYPRE_StructPFMGSetNumPreRelax(d_mg_data,
         d_num_pre_relax_steps);
      HYPRE_StructPFMGSetNumPostRelax(d_mg_data,
         d_num_post_relax_steps);
      HYPRE_StructPFMGSetup(d_mg_data,
         d_matrix,
         d_linear_rhs,
         d_linear_sol);
   }
}

void
CellPoissonHypreSolver::destroyHypreSolver()
{
   if (d_mg_data != 0) {
      if (d_use_smg) {
         HYPRE_StructSMGDestroy(d_mg_data);
      } else {
         HYPRE_StructPFMGDestroy(d_mg_data);
      }
      d_mg_data = 0;
   }
}

/*
 *************************************************************************
 *
 * Solve the linear system.  This routine assumes that the boundary
 * conditions and the matrix coefficients have been specified.
 *
 *************************************************************************
 */

int
CellPoissonHypreSolver::solveSystem(
   const int u,
   const int f,
   bool homogeneous_bc,
   bool initial_zero)
{
   if (d_physical_bc_coef_strategy == 0) {
      TBOX_ERROR(
         d_object_name << ": No BC coefficient strategy object!\n"
                       << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                       << "to specify the boundary conidition.  Do it before\n"
                       << "calling solveSystem.");
   }
   // Tracer t("CellPoissonHypreSolver::solveSystem");

   t_solve_system->start();

   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(d_ln));
   TBOX_ASSERT(u >= 0);
   TBOX_ASSERT(
      u < level->getPatchDescriptor()->getMaxNumberRegisteredComponents());
   TBOX_ASSERT(f >= 0);
   TBOX_ASSERT(
      f < level->getPatchDescriptor()->getMaxNumberRegisteredComponents());

   if (d_physical_bc_coef_strategy == &d_physical_bc_simple_case) {
      /*
       * If we are using the simple bc implementation, the final piece
       * of information it requires is the Dirichlet boundary value
       * set in the ghost cells.  Now that we have the ghost cell data,
       * we can complete the boundary condition setup.
       */
      d_physical_bc_simple_case.cacheDirichletData(u);
   }

   /*
    * Modify right-hand-side to account for boundary conditions and
    * Modify right-hand-side to account for boundary conditions and
    * copy solution and right-hand-side to HYPRE structures.
    */

   const hier::IntVector no_ghosts(d_dim, 0);
   const hier::IntVector ghosts(d_dim, 1);

   /*
    * At coarse-fine boundaries, we expect ghost cells to have correct
    * values to be used in our bc, so u provides the ghost cell data.
    * Assume that the user only provided data for the immediate first
    * ghost cell, so pass zero for the number of extensions fillable.
    */
   d_cf_bc_coef.setGhostDataId(u, hier::IntVector::getZero(d_dim));

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      const hier::Box box = patch->getBox();

      /*
       * Set up variable data needed to prepare linear system solver.
       */
      std::shared_ptr<pdat::CellData<double> > u_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(u)));
      TBOX_ASSERT(u_data_);
      pdat::CellData<double>& u_data = *u_data_;
      pdat::CellData<double> rhs_data(box, 1, no_ghosts,
                                      d_allocator);

      /*
       * Copy rhs and solution from the hierarchy into HYPRE structures.
       * For rhs, add in the contribution from boundary conditions, if
       * needed.  If boundary condition is homogenous, this only adds
       * zero, so we skip it.
       */
      if (!initial_zero) { 
         copyToHypre(d_linear_sol, u_data, d_soln_depth, box);
      } else {
         HYPRE_StructVectorSetConstantValues(d_linear_sol, 0.0);
      }
      rhs_data.copy(*(patch->getPatchData(f)));
      if (!homogeneous_bc) {
         /*
          * Add g*A*k0(a) from physical and coarse-fine boundaries to rhs.
          */
         add_gAk0_toRhs(*patch,
            patch->getPatchGeometry()->getCodimensionBoundaries(1),
            d_physical_bc_coef_strategy,
            rhs_data);
         add_gAk0_toRhs(*patch,
            d_cf_boundary->getBoundaries(patch->getGlobalId(), 1),
            &d_cf_bc_coef,
            rhs_data);
      }
      copyToHypre(d_linear_rhs, rhs_data, d_rhs_depth, box);

   } // end patch loop

   /*
    * Reset invalid ghost data id to help detect erroneous further use.
    */
   d_cf_bc_coef.setGhostDataId(-1, hier::IntVector::getZero(d_dim));

   /*
    * Finish assembly of the vectors
    */
   HYPRE_StructVectorAssemble(d_linear_sol);

   HYPRE_StructVectorAssemble(d_linear_rhs);

   /*
    * Solve the system - zero means convergence
    * Solve takes the same arguments as Setup
    */

   if (d_print_solver_info) {
      HYPRE_StructVectorPrint("sol0.out", d_linear_sol, 1);
      HYPRE_StructMatrixPrint("mat0.out", d_matrix, 1);
      HYPRE_StructVectorPrint("rhs.out", d_linear_rhs, 1);
   }

   if (d_use_smg) {
      // HYPRE_StructSMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructSMGSetTol(d_mg_data, d_relative_residual_tol);
      /* converge = */ HYPRE_StructSMGSolve(d_mg_data,
         d_matrix,
         d_linear_rhs,
         d_linear_sol);
   } else {
      // HYPRE_StructPFMGSetMaxIter(d_mg_data, d_max_iterations);
      HYPRE_StructPFMGSetTol(d_mg_data, d_relative_residual_tol);
      /* converge = */ HYPRE_StructPFMGSolve(d_mg_data,
         d_matrix,
         d_linear_rhs,
         d_linear_sol);
   }

   if (d_print_solver_info) {
      HYPRE_StructMatrixPrint("mat.out", d_matrix, 1);
      HYPRE_StructVectorPrint("sol.out", d_linear_sol, 1);
   }

   if (d_use_smg) {
      HYPRE_StructSMGGetNumIterations(d_mg_data,
         &d_number_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(d_mg_data,
         &d_relative_residual_norm);
   } else {
      HYPRE_StructPFMGGetNumIterations(d_mg_data,
         &d_number_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(d_mg_data,
         &d_relative_residual_norm);
   }

   /*
    * Pull the solution vector out of the HYPRE structures
    */
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;
      std::shared_ptr<pdat::CellData<double> > u_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch->getPatchData(u)));

      TBOX_ASSERT(u_data_);

      pdat::CellData<double>& u_data = *u_data_;
      copyFromHypre(u_data,
         d_soln_depth,
         d_linear_sol,
         patch->getBox());
   }

   t_solve_system->stop();

   return d_relative_residual_norm <= d_relative_residual_tol;
}

void
CellPoissonHypreSolver::computeDiagonalEntries(
   pdat::CellData<double>& diagonal,
   const pdat::CellData<double>& C_data,
   const pdat::SideData<double>& off_diagonal,
   const hier::Box& patch_box)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(d_dim,
      diagonal,
      C_data,
      off_diagonal,
      patch_box);

   const hier::Index& patch_lo = patch_box.lower();
   const hier::Index& patch_up = patch_box.upper();
   const double c = 1.0, d = 1.0;
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(compdiagvariablec2d, COMPDIAGVARIABLEC2D) (diagonal.getPointer(),
         C_data.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &c, &d);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(compdiagvariablec3d, COMPDIAGVARIABLEC3D) (diagonal.getPointer(),
         C_data.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         off_diagonal.getPointer(2),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &patch_lo[2], &patch_up[2],
         &c, &d);
   }
}

void
CellPoissonHypreSolver::computeDiagonalEntries(
   pdat::CellData<double>& diagonal,
   const double C,
   const pdat::SideData<double>& off_diagonal,
   const hier::Box& patch_box)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(d_dim, diagonal, off_diagonal, patch_box);

   const hier::Index& patch_lo = patch_box.lower();
   const hier::Index& patch_up = patch_box.upper();
   const double c = 1.0, d = 1.0;
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(compdiagscalarc2d, COMPDIAGSCALARC2D) (diagonal.getPointer(),
         &C,
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &c, &d);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(compdiagscalarc3d, COMPDIAGSCALARC3D) (diagonal.getPointer(),
         &C,
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         off_diagonal.getPointer(2),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &patch_lo[2], &patch_up[2],
         &c, &d);
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
         << "DIM > 3 not supported." << std::endl);
   }
}

void
CellPoissonHypreSolver::computeDiagonalEntries(
   pdat::CellData<double>& diagonal,
   const pdat::SideData<double>& off_diagonal,
   const hier::Box& patch_box)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(d_dim, diagonal, off_diagonal, patch_box);

   const hier::Index& patch_lo = patch_box.lower();
   const hier::Index& patch_up = patch_box.upper();
   const double c = 1.0, d = 1.0;
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(compdiagzeroc2d, COMPDIAGZEROC2D) (diagonal.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &c, &d);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(compdiagzeroc3d, COMPDIAGZEROC3D) (diagonal.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         off_diagonal.getPointer(2),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &patch_lo[2], &patch_up[2],
         &c, &d);
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
         << "DIM > 3 not supported." << std::endl);
   }
}

void
CellPoissonHypreSolver::adjustBoundaryEntries(
   pdat::CellData<double>& diagonal,
   const pdat::SideData<double>& off_diagonal,
   const hier::Box& patch_box,
   const pdat::ArrayData<double>& acoef_data,
   const pdat::ArrayData<double>& bcoef_data,
   const hier::Box bccoef_box,
   pdat::ArrayData<double>& Ak0_data,
   const hier::BoundaryBox& trimmed_boundary_box,
   const double h[SAMRAI::MAX_DIM_VAL])
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY8(d_dim, diagonal, off_diagonal, patch_box,
      acoef_data, bcoef_data,
      bccoef_box, Ak0_data, trimmed_boundary_box);

   const hier::Index& patch_lo = patch_box.lower();
   const hier::Index& patch_up = patch_box.upper();
   const int location_index = trimmed_boundary_box.getLocationIndex();
   const hier::Index& lower = trimmed_boundary_box.getBox().lower();
   const hier::Index& upper = trimmed_boundary_box.getBox().upper();
   const hier::Box& Ak0_box = Ak0_data.getBox();
   if (d_dim == tbox::Dimension(2)) {
      SAMRAI_F77_FUNC(adjbdry2d, ADJBDRY2D) (diagonal.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         acoef_data.getPointer(),
         bcoef_data.getPointer(),
         &bccoef_box.lower()[0],
         &bccoef_box.upper()[0],
         &bccoef_box.lower()[1],
         &bccoef_box.upper()[1],
         Ak0_data.getPointer(),
         &Ak0_box.lower()[0],
         &Ak0_box.upper()[0],
         &Ak0_box.lower()[1],
         &Ak0_box.upper()[1],
         &lower[0], &upper[0],
         &location_index, h);
   } else if (d_dim == tbox::Dimension(3)) {
      SAMRAI_F77_FUNC(adjbdry3d, ADJBDRY3D) (diagonal.getPointer(),
         off_diagonal.getPointer(0),
         off_diagonal.getPointer(1),
         off_diagonal.getPointer(2),
         &patch_lo[0], &patch_up[0],
         &patch_lo[1], &patch_up[1],
         &patch_lo[2], &patch_up[2],
         acoef_data.getPointer(),
         bcoef_data.getPointer(),
         &bccoef_box.lower()[0],
         &bccoef_box.upper()[0],
         &bccoef_box.lower()[1],
         &bccoef_box.upper()[1],
         &bccoef_box.lower()[2],
         &bccoef_box.upper()[2],
         Ak0_data.getPointer(),
         &Ak0_box.lower()[0],
         &Ak0_box.upper()[0],
         &Ak0_box.lower()[1],
         &Ak0_box.upper()[1],
         &Ak0_box.lower()[2],
         &Ak0_box.upper()[2],
         &lower[0], &upper[0],
         &location_index, h);
   } else {
      TBOX_ERROR("CellPoissonHypreSolver error...\n"
         << "DIM > 3 not supported." << std::endl);
   }
}

void
CellPoissonHypreSolver::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_Ak0_var[d].reset();
   }
}

}
}

#endif
