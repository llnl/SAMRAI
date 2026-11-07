/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:
 *
 ************************************************************************/

#include "CVODEModel.h"

#if defined(HAVE_SUNDIALS) && defined(HAVE_HYPRE)

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

//integer constants for boundary conditions
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary dat
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"

#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

// Define class version number
#define CVODE_MODEL_VERSION (1)

// This is used in the cell tagging routine.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

extern "C" {
void SAMRAI_F77_FUNC(comprhs2d, COMPRHS2D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *,
   const double *, const double *,
   double *);
#ifdef USE_FAC_PRECONDITIONER
void SAMRAI_F77_FUNC(setneufluxvalues2d, SETNEUFLUXVALUES2D) (
   const int&, const int&,
   const int&, const int&,
   const int *, const double *,
   int *, int *, int *, int *,
   double *, double *, double *, double *);
#endif
void SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double *,
   const double *, const double *, const double *,
   double *);
#ifdef USE_FAC_PRECONDITIONER
void SAMRAI_F77_FUNC(setneufluxvalues3d, SETNEUFLUXVALUES3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *,
   int *, int *, int *, int *, int *, int *,
   double *, double *, double *, double *, double *, double *);
#endif
}

/*************************************************************************
 *
 * Constructor and Destructor for CVODEModel class.
 *
 ************************************************************************/

CVODEModel::CVODEModel(
   const std::string& object_name,
   const Dimension& dim,
   std::shared_ptr<CellPoissonFACSolver> fac_solver,
   std::shared_ptr<Database> input_db,
   std::shared_ptr<CartesianGridGeometry> grid_geom):
   RefinePatchStrategy(),
   CoarsenPatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_allocator(tbox::AllocatorDatabase::getDatabase()->getDefaultAllocator()),
   d_soln_var(new CellVariable<double>(dim, "soln", d_allocator)),
   d_FAC_solver(fac_solver),
   d_grid_geometry(grid_geom)
{
   /*
    * set up variables and contexts
    */
   VariableDatabase* variable_db = VariableDatabase::getDatabase();

   d_cur_cxt = variable_db->getContext("CURRENT");
   d_scr_cxt = variable_db->getContext("SCRATCH");

   d_soln_cur_id = variable_db->registerVariableAndContext(d_soln_var,
         d_cur_cxt,
         IntVector(d_dim, 0));
   d_soln_scr_id = variable_db->registerVariableAndContext(d_soln_var,
         d_scr_cxt,
         IntVector(d_dim, 1));
#ifdef USE_FAC_PRECONDITIONER
   d_diff_var.reset(new SideVariable<double>(d_dim, "diffusion",
         hier::IntVector::getOne(d_dim), d_allocator));

   d_diff_id = variable_db->registerVariableAndContext(d_diff_var,
         d_cur_cxt,
         IntVector(d_dim, 0));

   /*
    * Set default values for preconditioner.
    */
   d_use_neumann_bcs = false;

   d_current_soln_time = 0.;

#endif

   /*
    * Print solver data.
    */
   d_print_solver_info = false;

   /*
    * Counters.
    */
   d_number_rhs_eval = 0;
   d_number_precond_setup = 0;
   d_number_precond_solve = 0;

   /*
    * Boundary condition initialization.
    */
   if (d_dim == Dimension(2)) {
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
      MathUtilities<double>::setVectorToSignalingNaN(d_bdry_edge_val);
   } else if (d_dim == Dimension(3)) {
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
      MathUtilities<double>::setVectorToSignalingNaN(d_bdry_face_val);
   }

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

#ifdef USE_FAC_PRECONDITIONER
   /*
    * Construct outerface variable to hold boundary flags and Neumann fluxes.
    */
   if (d_use_neumann_bcs) {
      d_flag_var.reset(new OuterfaceVariable<int>(d_dim, "bdryflag",
                                                  d_allocator));
      d_flag_id = variable_db->registerVariableAndContext(d_flag_var,
            d_cur_cxt,
            IntVector(d_dim, 0));
      d_neuf_var.reset(new OuterfaceVariable<double>(d_dim, "neuflux",
                                                     d_allocator));
      d_neuf_id = variable_db->registerVariableAndContext(d_neuf_var,
            d_cur_cxt,
            IntVector(d_dim, 0));
   } else {
      d_flag_id = -1;
      d_neuf_id = -1;
   }

   /*
    * Set boundary types for FAC preconditioner.
    *  bdry_types holds a flag where 0 = dirichlet, 1 = neumann
    */
   if (d_dim == Dimension(2)) {
      for (int i = 0; i < NUM_2D_EDGES; ++i) {
         d_bdry_types[i] = 0;
         if (d_scalar_bdry_edge_conds[i] == BdryCond::DIRICHLET) d_bdry_types[i] = 0;
         if (d_scalar_bdry_edge_conds[i] == BdryCond::NEUMANN) d_bdry_types[i] = 1;
      }
   } else if (d_dim == Dimension(3)) {
      for (int i = 0; i < NUM_3D_FACES; ++i) {
         d_bdry_types[i] = 0;
         if (d_scalar_bdry_face_conds[i] == BdryCond::DIRICHLET) d_bdry_types[i] = 0;
         if (d_scalar_bdry_face_conds[i] == BdryCond::NEUMANN) d_bdry_types[i] = 1;
      }
   }
#endif

   /*
    * Postprocess boundary data from input/restart values.  Note: scalar
    * quantity in this problem cannot have reflective boundary conditions
    * so we reset them to BdryCond::FLOW.
    */
   if (d_dim == Dimension(2)) {
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
               CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                  i, d_scalar_bdry_node_conds[i]);
         }
      }
   } else if (d_dim == Dimension(3)) {
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
               CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
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
               CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                  i, d_scalar_bdry_node_conds[i]);
         }
      }

   }

}

CVODEModel::~CVODEModel()
{
   std::shared_ptr<SAMRAIVectorReal<double> > soln_samvect =
      Sundials_SAMRAIVector::getSAMRAIVector(d_solution_vector);
   Sundials_SAMRAIVector::destroySundialsVector(d_solution_vector);

   soln_samvect->freeVectorComponents();
   soln_samvect.reset();

   // if (d_level_solver_allocated) delete d_level_solver;
   // d_level_solver_allocated = false;
   // if (d_FAC_solver_allocated) delete d_FAC_solver;
   // d_FAC_solver_allocated = false;

}

/*************************************************************************
 *
 * Methods inherited from mesh::StandardTagAndInitStrategy.
 *
 ************************************************************************/

void
CVODEModel::initializeLevelData(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(time);
   NULL_USE(old_level);
   NULL_USE(allocate_data);

   // This method is empty because initialization is taken care of
   // by the setInitialConditions() method below.  If there is any
   // data that is not managed inside the SAMRAI CVODESolver class
   // but that must be set on the level, do it here.

}

void
CVODEModel::resetHierarchyConfiguration(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   // This method is empty because this example does not exercise the
   // situation when the grid changes, so it effectively is never called.
   // This is a subject for future work...
}

/*
 *************************************************************************
 *
 * Cell tagging and patch level data initialization routines declared
 * in the GradientDetectorStrategy interface.  They are used to
 * construct the hierarchy initially.
 *
 *************************************************************************
 */

void
CVODEModel::applyGradientDetector(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(time);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);

   std::shared_ptr<PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<Patch>& patch = *p;

      std::shared_ptr<CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<CellData<int>, PatchData>(
            patch->getPatchData(tag_index)));
      TBOX_ASSERT(tag_data);

      // dumb implementation that tags all cells.
      tag_data->fillAll(TRUE);
   }
}

/*
 *************************************************************************
 *
 * Methods inherited from RefinePatchStrategy.
 *
 ***********************************************************************
 */

void
CVODEModel::setPhysicalBoundaryConditions(
   Patch& patch,
   const double time,
   const IntVector& ghost_width_to_fill)
{
   NULL_USE(time);

   std::shared_ptr<CellData<double> > soln_data(
      SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
         patch.getPatchData(d_soln_scr_id)));

   TBOX_ASSERT(soln_data);

   IntVector ghost_cells(soln_data->getGhostCellWidth());

   if (d_dim == Dimension(2)) {

      /*
       * Set boundary conditions for cells corresponding to patch edges.
       */
      CartesianBoundaryUtilities2::
      fillEdgeBoundaryData("soln_data", soln_data,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_edge_val);

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */

      CartesianBoundaryUtilities2::
      fillNodeBoundaryData("soln_data", soln_data,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_edge_val);

   } else if (d_dim == Dimension(3)) {

      /*
       *  Set boundary conditions for cells corresponding to patch faces.
       */
      CartesianBoundaryUtilities3::
      fillFaceBoundaryData("soln_data", soln_data,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_face_conds,
         d_bdry_face_val);

      /*
       *  Set boundary conditions for cells corresponding to patch edges.
       */
      CartesianBoundaryUtilities3::
      fillEdgeBoundaryData("soln_data", soln_data,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_edge_conds,
         d_bdry_face_val);

      /*
       *  Set boundary conditions for cells corresponding to patch nodes.
       */
      CartesianBoundaryUtilities3::
      fillNodeBoundaryData("soln_data", soln_data,
         patch,
         ghost_width_to_fill,
         d_scalar_bdry_node_conds,
         d_bdry_face_val);

   }

//    plog << "----Boundary Conditions "  << std::endl;
//    soln_data->print(soln_data->getGhostBox());

}

void
CVODEModel::preprocessRefine(
   Patch& fine,
   const Patch& coarse,
   const Box& fine_box,
   const IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void
CVODEModel::postprocessRefine(
   Patch& fine,
   const Patch& coarse,
   const Box& fine_box,
   const IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

/*************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy.
 *
 ************************************************************************/

void
CVODEModel::preprocessCoarsen(
   Patch& coarse,
   const Patch& fine,
   const Box& coarse_box,
   const IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void
CVODEModel::postprocessCoarsen(
   Patch& coarse,
   const Patch& fine,
   const Box& coarse_box,
   const IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunction
 *
 ************************************************************************/

int
CVODEModel::evaluateRHSFunction(
   double time,
   SundialsAbstractVector* y,
   SundialsAbstractVector* y_dot)
{
   /*
    * Convert Sundials vectors to SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > y_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y));
   std::shared_ptr<SAMRAIVectorReal<double> > y_dot_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y_dot));

   std::shared_ptr<PatchHierarchy> hierarchy(y_samvect->getPatchHierarchy());

   /*
    * Compute max norm of solution vector.
    */
   //std::shared_ptr<HierarchyDataOpsReal<double> > hierops(
   //   new HierarchyCellDataOpsReal<double>(hierarchy));
   //double max_norm = hierops->maxNorm(y_samvect->
   //                                   getComponentDescriptorIndex(0));

   if (d_print_solver_info) {
      pout << "\t\tEval RHS: "
           << "\n   \t\t\ttime = " << time
           << "\n   \t\t\ty_maxnorm = " << y_samvect->maxNorm()
           << std::endl;
   }

   /*
    * Allocate scratch space and fill ghost cells in the solution vector
    * 1) Create a refine algorithm
    * 2) Register with the algorithm the current & scratch space, along
    *    with a refine operator.
    * 3) Use the refine algorithm to construct a refine schedule
    * 4) Use the refine schedule to fill data on fine level.
    */
   std::shared_ptr<RefineAlgorithm> bdry_fill_alg(
      new RefineAlgorithm());
   std::shared_ptr<RefineOperator> refine_op(d_grid_geometry->
                                               lookupRefineOperator(d_soln_var,
                                                  "CONSERVATIVE_LINEAR_REFINE"));
   bdry_fill_alg->registerRefine(d_soln_scr_id,  // dest
      y_samvect->
      getComponentDescriptorIndex(0),                            // src
      d_soln_scr_id,                            // scratch
      refine_op);

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));
      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }

      // Note:  a pointer to "this" tells the refine schedule to invoke
      // the setPhysicalBCs defined in this class.
      std::shared_ptr<RefineSchedule> bdry_fill_alg_schedule(
         bdry_fill_alg->createSchedule(level,
            ln - 1,
            hierarchy,
            this));

      bdry_fill_alg_schedule->fillData(time);
   }

   /*
    * Step through the levels and compute rhs
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
         const std::shared_ptr<Patch>& patch = *ip;

         std::shared_ptr<CellData<double> > y(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_soln_scr_id)));
         std::shared_ptr<SideData<double> > diff(
            SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
               patch->getPatchData(d_diff_id)));
         std::shared_ptr<CellData<double> > rhs(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(y_dot_samvect->getComponentDescriptorIndex(0))));
         TBOX_ASSERT(y);
         TBOX_ASSERT(diff);
         TBOX_ASSERT(rhs);

         const Index ifirst(patch->getBox().lower());
         const Index ilast(patch->getBox().upper());

         const std::shared_ptr<CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();

         IntVector ghost_cells(y->getGhostCellWidth());

         /*
          * 1 eqn radiation diffusion
          */
         if (d_dim == Dimension(2)) {
            SAMRAI_F77_FUNC(comprhs2d, COMPRHS2D) (
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ghost_cells(0), ghost_cells(1),
               dx,
               y->getPointer(),
               diff->getPointer(0),
               diff->getPointer(1),
               rhs->getPointer());
         } else if (d_dim == Dimension(3)) {
            SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D) (
               ifirst(0), ilast(0),
               ifirst(1), ilast(1),
               ifirst(2), ilast(2),
               ghost_cells(0), ghost_cells(1),
               ghost_cells(2),
               dx,
               y->getPointer(),
               diff->getPointer(0),
               diff->getPointer(1),
               diff->getPointer(2),
               rhs->getPointer());
         }

      } // loop over patches
   } // loop over levels

   /*
    * Deallocate scratch space.
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(d_soln_scr_id);
   }

   /*
    * record current time and increment counter for number of RHS
    * evaluations.
    */
   d_current_soln_time = time;
   ++d_number_rhs_eval;

   return 0;
}

/*
 *****************************************************************
 *
 * Set up FAC preconditioner for Jacobian system.  Here we
 * use the FAC hierarchy solver in SAMRAI which automatically sets
 * up the composite grid system and uses hypre as a solver on each
 * level.
 *
 *****************************************************************
 */

int CVODEModel::CVSpgmrPrecondSet(
   double t,
   SundialsAbstractVector* y,
   SundialsAbstractVector* fy,
   int jok,
   int* jcurPtr,
   double gamma)
{
#ifndef USE_FAC_PRECONDITIONER
   NULL_USE(t);
   NULL_USE(y);
   NULL_USE(gamma);
#endif
   NULL_USE(fy);
   NULL_USE(jok);
   NULL_USE(jcurPtr);

#ifdef USE_FAC_PRECONDITIONER

   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > y_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y));

   std::shared_ptr<PatchHierarchy> hierarchy(
      y_samvect->getPatchHierarchy());

   int y_indx = y_samvect->getComponentDescriptorIndex(0);

   /*
    * Construct refine algorithm to fill boundaries of solution vector
    */
   RefineAlgorithm fill_soln_vector_bounds;
   std::shared_ptr<RefineOperator> refine_op(d_grid_geometry->
                                               lookupRefineOperator(d_soln_var,
                                                  "CONSERVATIVE_LINEAR_REFINE"));
   fill_soln_vector_bounds.registerRefine(d_soln_scr_id,
      y_samvect->getComponentDescriptorIndex(0),
      d_soln_scr_id,
      refine_op);

   /*
    * Construct coarsen algorithm to fill interiors on coarser levels
    * with solution on finer level.
    */
   CoarsenAlgorithm fill_soln_interior_on_coarser(d_dim);
   std::shared_ptr<CoarsenOperator> coarsen_op(d_grid_geometry->
                                                 lookupCoarsenOperator(d_soln_var,
                                                    "CONSERVATIVE_COARSEN"));

   fill_soln_interior_on_coarser.registerCoarsen(y_indx,
      y_indx,
      coarsen_op);

   /*
    * Step through levels - largest to smallest
    */
   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<PatchLevel> level(
         hierarchy->getPatchLevel(amr_level));

      std::shared_ptr<RefineSchedule> fill_soln_vector_bounds_sched =
         fill_soln_vector_bounds.createSchedule(level,
            amr_level - 1,
            hierarchy,
            this);

      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }
      fill_soln_vector_bounds_sched->fillData(t);

      /*
       * Construct a coarsen schedule for all levels larger than coarsest,
       * and fill interiors of solution vector on coarser levels using fine
       * data.
       */
      if (amr_level > 0) {
         std::shared_ptr<PatchLevel> coarser_level(
            hierarchy->getPatchLevel(amr_level - 1));

         std::shared_ptr<CoarsenSchedule> fill_soln_interior_on_coarser_sched(
            fill_soln_interior_on_coarser.createSchedule(coarser_level,
               level));

         fill_soln_interior_on_coarser_sched->coarsenData();
      }

      for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
         const std::shared_ptr<Patch>& patch = *p;

         const Index ifirst(patch->getBox().lower());
         const Index ilast(patch->getBox().upper());

         std::shared_ptr<SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
               patch->getPatchData(d_diff_id)));
         TBOX_ASSERT(diffusion);

         diffusion->fillAll(1.0);

         TBOX_ASSERT((t - d_current_soln_time) >= 0.);

         /*
          * Set Neumann fluxes and flag array (if desired)
          */
         if (d_use_neumann_bcs) {

            std::shared_ptr<OuterfaceData<int> > flag_data(
               SAMRAI_SHARED_PTR_CAST<OuterfaceData<int>, PatchData>(
                  patch->getPatchData(d_flag_id)));
            std::shared_ptr<OuterfaceData<double> > neuf_data(
               SAMRAI_SHARED_PTR_CAST<OuterfaceData<double>, PatchData>(
                  patch->getPatchData(d_neuf_id)));
            TBOX_ASSERT(flag_data);
            TBOX_ASSERT(neuf_data);

            /*
             * Outerface data access:
             *    neuf_data->getPointer(axis,face);
             * where axis specifies X, Y, or Z (0,1,2 respectively)
             * and face specifies lower or upper (0,1 respectively)
             */

            if (d_dim == Dimension(2)) {
               SAMRAI_F77_FUNC(setneufluxvalues2d, SETNEUFLUXVALUES2D) (
                  ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  d_bdry_types,
                  &d_bdry_edge_val[0],
                  flag_data->getPointer(0, 0), // x lower
                  flag_data->getPointer(0, 1), // x upper
                  flag_data->getPointer(1, 0), // y lower
                  flag_data->getPointer(1, 1), // y upper
                  neuf_data->getPointer(0, 0), // x lower
                  neuf_data->getPointer(0, 1), // x upper
                  neuf_data->getPointer(1, 0), // y lower
                  neuf_data->getPointer(1, 1)); // y upper
            } else if (d_dim == Dimension(3)) {
               SAMRAI_F77_FUNC(setneufluxvalues3d, SETNEUFLUXVALUES3D) (
                  ifirst(0), ilast(0),
                  ifirst(1), ilast(1),
                  ifirst(2), ilast(2),
                  d_bdry_types,
                  &d_bdry_face_val[0],
                  flag_data->getPointer(0, 0), // x lower
                  flag_data->getPointer(0, 1), // x upper
                  flag_data->getPointer(1, 0), // y lower
                  flag_data->getPointer(1, 1), // y upper
                  flag_data->getPointer(2, 0), // z lower
                  flag_data->getPointer(2, 1), // z lower
                  neuf_data->getPointer(0, 0), // x lower
                  neuf_data->getPointer(0, 1), // x upper
                  neuf_data->getPointer(1, 0), // y lower
                  neuf_data->getPointer(1, 1), // y upper
                  neuf_data->getPointer(2, 0), // z lower
                  neuf_data->getPointer(2, 1)); // z upper
            }
         }

      } // patch loop

      level->deallocatePatchData(d_soln_scr_id);

   } // level loop

   /*
    * Set boundaries.  The "bdry_types" array holds a set of integers
    * where 0 = dirichlet and 1 = neumann boundary conditions.
    */
   if (d_use_neumann_bcs) {
      d_FAC_solver->setBoundaries("Mixed", d_neuf_id, d_flag_id, d_bdry_types);
   } else {
      d_FAC_solver->setBoundaries("Dirichlet");
   }

   d_FAC_solver->setCConstant(1.0 / gamma);
   d_FAC_solver->setDPatchDataId(d_diff_id);

   /*
    * increment counter for number of precond setup calls
    */
   ++d_number_precond_setup;

#endif
   /*
    * We return 0 or 1 here - 0 if it passes, 1 if it fails.  For now,
    * just be optimistic and return 0. Eventually we should add some
    * assertion handling above to set what this value should be.
    */
   return 0;
}

/*
 *************************************************************************
 *
 * Apply preconditioner where right-hand-side is "r" and "z" is the
 * solution.   This routine assumes that the preconditioner setup call
 * has already been invoked.  Return 0 if preconditioner fails;
 * return 1 otherwise.
 *
 *************************************************************************
 */

int CVODEModel::CVSpgmrPrecondSolve(
   double t,
   SundialsAbstractVector* y,
   SundialsAbstractVector* fy,
   SundialsAbstractVector* r,
   SundialsAbstractVector* z,
   double gamma,
   double delta,
   int lr)
{
   NULL_USE(y);
   NULL_USE(fy);
#ifndef USE_FAC_PRECONDITIONER
   NULL_USE(gamma);
#endif
   NULL_USE(delta);
   NULL_USE(lr);

#ifdef USE_FAC_PRECONDITIONER

   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > r_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(r));
   std::shared_ptr<SAMRAIVectorReal<double> > z_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(z));

   int ret_val = 0;

   std::shared_ptr<PatchHierarchy> hierarchy(
      r_samvect->getPatchHierarchy());

   int r_indx = r_samvect->getComponentDescriptorIndex(0);
   int z_indx = z_samvect->getComponentDescriptorIndex(0);
   /******************************************************************
    *
    * We need to supply to the FAC solver a "version" of the z vector
    * that contains ghost cells.  The operations below allocate
    * on the patches a scratch context of the solution vector z and
    * fill it with z vector data
    *
    *****************************************************************/

   /*
    * Construct a communication schedule which will fill ghosts of
    * soln_scratch with z vector data (z -> soln_scratch).
    */
   RefineAlgorithm fill_z_vector_bounds;
   std::shared_ptr<RefineOperator> refine_op(d_grid_geometry->
                                               lookupRefineOperator(d_soln_var,
                                                  "CONSERVATIVE_LINEAR_REFINE"));
   fill_z_vector_bounds.registerRefine(d_soln_scr_id,
      z_indx,
      d_soln_scr_id,
      refine_op);

   /*
    * Set initial guess for z (if applicable) and copy z data into the
    * solution scratch context.
    */
   int ln;
   for (ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      if (!level->checkAllocated(d_soln_scr_id)) {
         level->allocatePatchData(d_soln_scr_id);
      }

      for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

         const std::shared_ptr<Patch>& patch = *p;

         std::shared_ptr<CellData<double> > z_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(z_indx)));
         TBOX_ASSERT(z_data);

         /*
          * Set initial guess for z here.
          */
         z_data->fillAll(0.);

         /*
          * Scale RHS by 1/gamma
          */
         PatchCellDataOpsReal<double> math_ops;
         std::shared_ptr<CellData<double> > r_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(r_indx)));
         TBOX_ASSERT(r_data);
         math_ops.scale(r_data, 1.0 / gamma, r_data, r_data->getBox());

         /*
          * Copy interior data from z vector to soln_scratch
          */
         std::shared_ptr<CellData<double> > z_scr_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_soln_scr_id)));
         TBOX_ASSERT(z_scr_data);
         z_scr_data->copy(*z_data);
      }

      /*
       * Fill ghost boundaries of soln_scratch.
       * Construct a schedule for each level, from the algorithm
       * constructed above.
       */

      std::shared_ptr<RefineSchedule> fill_z_vector_bounds_sched(
         fill_z_vector_bounds.createSchedule(level,
            ln - 1,
            hierarchy,
            this));

      fill_z_vector_bounds_sched->fillData(t);

   }

   /******************************************************************
   *
   * Apply the FAC solver.  It solves the system Az=r with the
   * format "solveSystem(z, r)". A was constructed in the precondSetup()
   * method.
   *
   ******************************************************************/

   if (d_print_solver_info) {
      pout << "\t\tBefore FAC Solve (Az=r): "
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm()
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm()
           << "\n   \t\t\tr_l2norm = " << r_samvect->L2Norm()
           << "\n   \t\t\tr_maxnorm = " << r_samvect->maxNorm()
           << std::endl;
   }
   /*
    * Set paramemters in the FAC solver.  It solves the system Az=r.
    * Here we supply the max norm of r in order to scale the
    * residual (i.e. residual = Az - r) to properly scale the convergence
    * error.
    */

   const int coarsest_solve_ln = 0;
   const int finest_solve_ln = 0;
   /*
    * Note: I don't know why we are only solving on level 0 here.
    * When upgrading to the new FAC solver from the old, I noticed
    * that the old solver only solved on level 0.  BTNG.
    */
   bool converge = d_FAC_solver->solveSystem(d_soln_scr_id,
         r_indx,
         hierarchy,
         coarsest_solve_ln,
         finest_solve_ln);

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver->getConvergenceFactors(avg_convergence, final_convergence);
      pout << "   \t\t\tFinal Residual Norm: "
           << d_FAC_solver->getResidualNorm() << std::endl;
      pout << "   \t\t\tFinal Convergence Error: "
           << final_convergence << std::endl;
      pout << "   \t\t\tFinal Convergence Rate: "
           << avg_convergence << std::endl;
   }

   /******************************************************************
   *
   * The FAC solver has computed a solution to z but it is stored
   * in the soln_scratch data space.  Copy it from soln_scratch back
   * into the z vector.
   *
   ******************************************************************/
   for (ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
         const std::shared_ptr<Patch>& patch = *p;

         std::shared_ptr<CellData<double> > soln_scratch(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_soln_scr_id)));
         std::shared_ptr<CellData<double> > z(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(z_indx)));
         TBOX_ASSERT(soln_scratch);
         TBOX_ASSERT(z);

         z->copy(*soln_scratch);
      }

   }

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver->getConvergenceFactors(avg_convergence, final_convergence);
      pout << "\t\tAfter FAC Solve (Az=r): "
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm()
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm()
           << "\n   \t\t\tResidual Norm: " << d_FAC_solver->getResidualNorm()
           << "\n   \t\t\tConvergence Error: " << final_convergence
           << std::endl;
   }

   if (converge != true) {
      ret_val = 1;
   }

   /*
    * Increment counter for number of precond solves
    */
   ++d_number_precond_solve;

   return ret_val;

#else

   return 0;

#endif

}

/*************************************************************************
 *
 * Methods specific to CVODEModel class.
 *
 ************************************************************************/

void
CVODEModel::setupSolutionVector(
   std::shared_ptr<PatchHierarchy> hierarchy)
{
   /* create SAMRAIVector */
   std::shared_ptr<SAMRAIVectorReal<double> > soln_samvect(
      new SAMRAIVectorReal<double>(
         "solution",
         hierarchy,
         0,
         hierarchy->getFinestLevelNumber()));
   soln_samvect->addComponent(d_soln_var, d_soln_cur_id);

   /* allocate memory for vectors. */
   soln_samvect->allocateVectorData();

   /* create SundialsAbstractVector */
   d_solution_vector =
      Sundials_SAMRAIVector::createSundialsVector(soln_samvect);

#ifdef USE_FAC_PRECONDITIONER
   /*
    * Allocate memory for preconditioner variables.
    */

   const int nlevels = hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < nlevels; ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));
      TBOX_ASSERT(level);
      level->allocatePatchData(d_diff_id);
      if (d_use_neumann_bcs) {
         level->allocatePatchData(d_flag_id);
         level->allocatePatchData(d_neuf_id);
      }

   }
#endif

}

SundialsAbstractVector *
CVODEModel::getSolutionVector(
   void)
{
   return d_solution_vector;
}

/*
 *************************************************************************
 *
 * Set initial conditions for CVODE solver
 *
 *************************************************************************
 */
void
CVODEModel::setInitialConditions(
   SundialsAbstractVector* soln_init)
{
   std::shared_ptr<SAMRAIVectorReal<double> > soln_init_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(soln_init));

   std::shared_ptr<PatchHierarchy> hierarchy(
      soln_init_samvect->getPatchHierarchy());

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (int cn = 0; cn < soln_init_samvect->getNumberOfComponents(); ++cn) {
         for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
            const std::shared_ptr<Patch>& patch = *p;

            /*
             * Set initial conditions for y
             */
            std::shared_ptr<CellData<double> > y_init(
               SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
                  soln_init_samvect->getComponentPatchData(cn, *patch)));
            TBOX_ASSERT(y_init);
            y_init->fillAll(d_initial_value);

            /*
             * Set initial diffusion coeff values.
             * NOTE: in a "real" application, the diffusion coefficient is
             * some function of y.  Here, we just do a simple minded
             * approach and set it to 1.
             */
            std::shared_ptr<SideData<double> > diffusion(
               SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
                  patch->getPatchData(d_diff_id)));
            TBOX_ASSERT(diffusion);

            diffusion->fillAll(1.0);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Return array of program counters.  Currently, the array holds the
 * following entries:
 *    1) number of RHS evaluations
 *    2) number of precond setup calls
 *    3) number of precond solve calls
 * More counters may be added, as desired.
 *
 *************************************************************************
 */
void
CVODEModel::getCounters(
   std::vector<int>& counters)
{
   counters.resize(3);
   counters[0] = d_number_rhs_eval;
   counters[1] = d_number_precond_setup;
   counters[2] = d_number_precond_solve;
}

/*
 *************************************************************************
 *
 * Get data from input database.
 *
 *************************************************************************
 */
void
CVODEModel::getFromInput(
   std::shared_ptr<Database> input_db,
   bool is_from_restart)
{
   NULL_USE(is_from_restart);

   d_initial_value = input_db->getDoubleWithDefault("initial_value", 0.0);

   IntVector periodic(d_grid_geometry->getPeriodicShift(IntVector(d_dim,
                            1)));
   int num_per_dirs = 0;
   for (int id = 0; id < d_dim.getValue(); ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (input_db->keyExists("Boundary_data")) {
      std::shared_ptr<Database> boundary_db(
         input_db->getDatabase("Boundary_data"));

      if (d_dim == Dimension(2)) {
         CartesianBoundaryUtilities2::getFromInput(this,
            boundary_db,
            d_scalar_bdry_edge_conds,
            d_scalar_bdry_node_conds,
            periodic);
      } else if (d_dim == Dimension(3)) {
         CartesianBoundaryUtilities3::getFromInput(this,
            boundary_db,
            d_scalar_bdry_face_conds,
            d_scalar_bdry_edge_conds,
            d_scalar_bdry_node_conds,
            periodic);
      }

   } else {
      TBOX_WARNING(
         d_object_name << ": "
                       << "Key data `Boundary_data' not found in input. " << std::endl);
   }

#ifdef USE_FAC_PRECONDITIONER
   d_use_neumann_bcs =
      input_db->getBoolWithDefault("use_neumann_bcs", d_use_neumann_bcs);
   d_print_solver_info =
      input_db->getBoolWithDefault("print_solver_info", d_print_solver_info);
#endif

}

/*
 *************************************************************************
 *
 * Write data to  restart database.
 *
 *************************************************************************
 */
void CVODEModel::putToRestart(
   const std::shared_ptr<Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("CVODE_MODEL_VERSION", CVODE_MODEL_VERSION);

   restart_db->putDouble("d_initial_value", d_initial_value);

   restart_db->putIntegerVector("d_scalar_bdry_edge_conds",
      d_scalar_bdry_edge_conds);
   restart_db->putIntegerVector("d_scalar_bdry_node_conds",
      d_scalar_bdry_node_conds);

   if (d_dim == Dimension(2)) {
      restart_db->putDoubleVector("d_bdry_edge_val", d_bdry_edge_val);
   } else if (d_dim == Dimension(3)) {
      restart_db->putIntegerVector("d_scalar_bdry_face_conds",
         d_scalar_bdry_face_conds);
      restart_db->putDoubleVector("d_bdry_face_val", d_bdry_face_val);
   }

}

/*
 *************************************************************************
 *
 * Read data from restart database.
 *
 *************************************************************************
 */
void CVODEModel::getFromRestart()
{

   std::shared_ptr<Database> root_db(
      RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in the restart file.");
   }
   std::shared_ptr<Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("CVODE_MODEL_VERSION");
   if (ver != CVODE_MODEL_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   d_initial_value = db->getDouble("d_initial_value");

   d_scalar_bdry_edge_conds = db->getIntegerVector("d_scalar_bdry_edge_conds");
   d_scalar_bdry_node_conds = db->getIntegerVector("d_scalar_bdry_node_conds");

   if (d_dim == Dimension(2)) {
      d_bdry_edge_val = db->getDoubleVector("d_bdry_edge_val");
   } else if (d_dim == Dimension(3)) {
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

void CVODEModel::readDirichletBoundaryDataEntry(
   const std::shared_ptr<Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());

   if (d_dim == Dimension(2)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_edge_val);
   } else if (d_dim == Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_val);
   }
}

void CVODEModel::readNeumannBoundaryDataEntry(
   const std::shared_ptr<Database>& db,
   std::string& db_name,
   int bdry_location_index)
{
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());

   if (d_dim == Dimension(2)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_edge_val);
   } else if (d_dim == Dimension(3)) {
      readStateDataEntry(db,
         db_name,
         bdry_location_index,
         d_bdry_face_val);
   }
}

void CVODEModel::readStateDataEntry(
   std::shared_ptr<Database> db,
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
 * Prints class data - writes out info in class if assertion is thrown
 *
 *************************************************************************
 */

void CVODEModel::printClassData(
   std::ostream& os) const
{
   fflush(stdout);
   int j;

   os << "ptr CVODEModel = " << (CVODEModel *)this << std::endl;

   os << "d_object_name = " << d_object_name << std::endl;

   os << "d_soln_cur_id = " << d_soln_cur_id << std::endl;
   os << "d_soln_scr_id = " << d_soln_scr_id << std::endl;

   os << "d_initial_value = " << d_initial_value << std::endl;

   os << "Boundary Condition data..." << std::endl;
   if (d_dim == Dimension(2)) {
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
   } else if (d_dim == Dimension(3)) {
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

}

void CVODEModel::setPrintSolverInfo(
   const bool info)
{
   d_print_solver_info = info;
}
#endif // HAVE_SUNDIALS
