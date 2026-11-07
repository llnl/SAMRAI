/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Example demonstrating use of CVODE vectors.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#if !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)

/*
 *************************************************************************
 * If the library is not compiled with CVODE, print an error.
 * If we're running autotests, skip the error and compile an empty
 * class.
 *************************************************************************
 */
#if (TESTING != 1)
#error "This example requires SAMRAI be compiled with CVODE -and- HYPRE."
#endif

#else

/*
 * Header file for base classes.
 */
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/solv/CVODEAbstractFunctions.h"

/*
 * Header file for SAMRAI classes referenced in this class.
 */
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

#define USE_FAC_PRECONDITIONER
// comment out line below to invoke preconditioning
// #undef USE_FAC_PRECONDITIONER

#ifdef USE_FAC_PRECONDITIONER
#include "SAMRAI/solv/CellPoissonFACSolver.h"
#endif
/*
 * Header files for CVODE wrapper classes
 */
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"


#include <vector>
#include <memory>

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace xfer;
using namespace pdat;
using namespace math;
using namespace mesh;
using namespace geom;
using namespace solv;
using namespace appu;

/**
 * The cvode_Model class tests the CVODE-SAMRAI interface using
 * two problems: (1) y' = k * d/dx (dy/dx) and (2) y' = y.
 *
 * The choice of which problem to solve and other input parameters
 * are specified through the input database.
 *
 * Input Parameters:
 *
 *
 *
 *
 *    - \b Problem_type
 *       1 for diffusion equation, 2 for y' = y.  By default, the heat
 *       equation is solved.
 *
 *    - \b Diffusion_coefficient
 *       specifies the diffusion coefficient to use when the
 *       has been specified that the diffusion equation will be
 *       solved.
 *
 *    - \b Initial_condition_type
 *       0 for constant initial conditions, 1 for sinusoidal initial
 *       conditions
 *
 *    - \b Initial_value
 *       specifies the initial value to be used for all grid points
 *       when constant initial conditions is specified
 *
 *    - \b Boundary_value
 *       specifies what value should be used for the dirichlet
 *       boundary conditions applied to this problem.
 *
 */

class CVODEModel:
   public StandardTagAndInitStrategy,
   public RefinePatchStrategy,
   public CoarsenPatchStrategy,
   public BoundaryUtilityStrategy,
   public CVODEAbstractFunctions
{
public:
   /**
    * Default constructor for CVODEModel.
    */
   CVODEModel(
      const std::string& object_name,
      const Dimension& dim,
      std::shared_ptr<CellPoissonFACSolver> fac_solver,
      std::shared_ptr<Database> input_db,
      std::shared_ptr<CartesianGridGeometry> grid_geom);

   /**
    * Empty destructor for CVODEModel.
    */
   virtual ~CVODEModel();

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy.
 *
 ************************************************************************/

   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old
    * level.  The specific operations that occur when initializing level
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
   virtual void
   initializeLevelData(
      const std::shared_ptr<PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<PatchLevel>& old_level =
         std::shared_ptr<PatchLevel>(),
      const bool allocate_data = true);

   /**
    * After hierarchy levels have changed and data has been initialized on
    * the new levels, this routine can be used to reset any information
    * needed by the solution method that is particular to the hierarchy
    * configuration.  For example, the solution procedure may cache
    * communication schedules to amortize the cost of data movement on the
    * AMR patch hierarchy.  This function will be called by the gridding
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must
    * make the solution consistent across multiple levels after the hierarchy
    * is changed, this process may be invoked by this routine.  Of course the
    * details of these processes are determined by the particular solution
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void
   resetHierarchyConfiguration(
      const std::shared_ptr<PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Set tags to the specified tag value where refinement of the given
    * level should occur using the user-supplied gradient detector.  The
    * value "tag_index" is the index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument indicates
    * whether cells are being tagged on the level for the first time;
    * i.e., when the hierarchy is initially constructed.  If it is false,
    * it should be assumed that cells are being tagged at some later time
    * after the patch hierarchy was initially constructed.  This information
    * is provided since the application of the error estimator may be
    * different in each of those circumstances.
    */
   virtual void
   applyGradientDetector(
      const std::shared_ptr<PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * Option to output solver info.  Set to true to turn on, false to
    * turn off.
    */
   void
   setPrintSolverInfo(
      const bool info);

/*************************************************************************
 *
 * Methods inherited from RefinePatchStrategy.
 *
 ************************************************************************/

   /**
    * Set the data at patch boundaries corresponding to the physical domain
    * boundary.  The specific boundary conditions are determined by the user.
    */
   virtual void
   setPhysicalBoundaryConditions(
      Patch& patch,
      const double time,
      const IntVector& ghost_width_to_fill);

   /**
    * Perform user-defined refining operations.  This member function
    * is called before the other refining operators.  The preprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the preprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void
   preprocessRefine(
      Patch& fine,
      const Patch& coarse,
      const Box& fine_box,
      const IntVector& ratio);

   /**
    * Perform user-defined refining operations.  This member function
    * is called after the other refining operators.  The postprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the postprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void
   postprocessRefine(
      Patch& fine,
      const Patch& coarse,
      const Box& fine_box,
      const IntVector& ratio);

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual IntVector getRefineOpStencilWidth(const Dimension& dim) const
   {
      return IntVector(dim, 0);
   }

/*************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy.
 *
 ************************************************************************/

   /**
    * Perform user-defined coarsening operations.  This member function
    * is called before the other coarsening operators.  The preprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void
   preprocessCoarsen(
      Patch& coarse,
      const Patch& fine,
      const Box& coarse_box,
      const IntVector& ratio);

   /**
    * Perform user-defined coarsening operations.  This member function
    * is called after the other coarsening operators.  The postprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void
   postprocessCoarsen(
      Patch& coarse,
      const Patch& fine,
      const Box& coarse_box,
      const IntVector& ratio);

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual IntVector getCoarsenOpStencilWidth(const Dimension& dim) const
   {
      return IntVector(dim, 0);
   }

   /*!
    * @brief Return the dimension of this object.
    */
   const Dimension& getDim() const
   {
      return d_dim;
   }

/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunctions
 *
 ************************************************************************/

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *
    *
    *
    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *
    *
    *
    *
    * IMPORTANT: This function must not modify the vector y. (KTC??)
    */
   virtual int
   evaluateRHSFunction(
      double time,
      SundialsAbstractVector* y,
      SundialsAbstractVector* y_dot);

   virtual int
   CVSpgmrPrecondSet(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* fy,
      int jok,
      int* jcurPtr,
      double gamma);

   virtual int
   CVSpgmrPrecondSolve(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* fy,
      SundialsAbstractVector* r,
      SundialsAbstractVector* z,
      double gamma,
      double delta,
      int lr);

   virtual int
   applyProjection(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* corr,
      double epsProj,
      SundialsAbstractVector* err)
   {
      NULL_USE(t);
      NULL_USE(y);
      NULL_USE(corr);
      NULL_USE(epsProj);
      NULL_USE(err);

      return 0;
   }

   virtual int evaluateJTimesRHSFunction(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* y_dot)
   {
      NULL_USE(t);
      NULL_USE(y);
      NULL_USE(y_dot);

      return 0;
   }


/*************************************************************************
 *
 * Methods particular to CVODEModel class.
 *
 ************************************************************************/

   /**
    * Set up solution vector.
    */
   void
   setupSolutionVector(
      std::shared_ptr<PatchHierarchy> hierarchy);

   /**
    * Get pointer to the solution vector.
    */
   SundialsAbstractVector *
   getSolutionVector(
      void);

   /**
    * Set initial conditions for problem.
    */
   void
   setInitialConditions(
      SundialsAbstractVector* y_init);

   /**
    * Return array of program counters.
    */
   void
   getCounters(
      std::vector<int>& counters);

   /**
    * Writes state of CVODEModel object to the specified restart database.
    *
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void
   putToRestart(
      const std::shared_ptr<Database>& restart_db) const;

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * and NEUMANN boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face (in 3D) or edge (in 2D) to which the boundary
    * condition applies.
    */
   void
   readDirichletBoundaryDataEntry(
      const std::shared_ptr<Database>& db,
      std::string& db_name,
      int bdry_location_index);

   void
   readNeumannBoundaryDataEntry(
      const std::shared_ptr<Database>& db,
      std::string& db_name,
      int bdry_location_index);

   /**
    * Prints all class data members, if assertion is thrown.
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
   virtual void
   getFromInput(
      std::shared_ptr<Database> input_db,
      bool is_from_restart);

   virtual void
   getFromRestart();

   void
   readStateDataEntry(
      std::shared_ptr<Database> db,
      const std::string& db_name,
      int array_indx,
      std::vector<double>& uval);

   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   std::string d_object_name;

   const Dimension d_dim;

   /*
    * Pointer to solution vector
    */
   SundialsAbstractVector* d_solution_vector;

   tbox::ResourceAllocator d_allocator;

   /*
    * Variables
    */
   std::shared_ptr<CellVariable<double> > d_soln_var;

   /*
    * Variable Contexts
    */
   std::shared_ptr<VariableContext> d_cur_cxt;
   std::shared_ptr<VariableContext> d_scr_cxt;

   /*
    * Patch Data ids
    */
   int d_soln_cur_id;
   int d_soln_scr_id;

#ifdef USE_FAC_PRECONDITIONER
   std::shared_ptr<SideVariable<double> > d_diff_var;
   std::shared_ptr<OuterfaceVariable<int> > d_flag_var;
   std::shared_ptr<OuterfaceVariable<double> > d_neuf_var;

   int d_diff_id;
   int d_flag_id;
   int d_neuf_id;
   int d_bdry_types[2 * MAX_DIM_VAL];

   std::shared_ptr<CellPoissonFACSolver> d_FAC_solver;
   bool d_FAC_solver_allocated;
   bool d_level_solver_allocated;
   bool d_use_neumann_bcs;

   double d_current_soln_time;
#endif

   /*
    * Print CVODE solver information
    */
   bool d_print_solver_info;

   /*
    * Grid geometry
    */
   std::shared_ptr<CartesianGridGeometry> d_grid_geometry;

   /*
    * Initial value
    */
   double d_initial_value;

   /*
    * Program counters
    *   1 - number of RHS evaluations
    *   2 - number of precond setups
    *   3 - number of precond solves
    */
   int d_number_rhs_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   /*
    * Boundary condition cases and boundary values.
    * Options are: FLOW, REFLECT, DIRICHLET, NEUMANN
    * and variants for nodes and edges.
    *
    * Input file values are read into these arrays.
    */
   std::vector<int> d_scalar_bdry_edge_conds;
   std::vector<int> d_scalar_bdry_node_conds;
   std::vector<int> d_scalar_bdry_face_conds; // Only used for 3D.

   /*
    * Boundary condition cases for scalar and vector (i.e., depth > 1)
    * variables.  These are post-processed input values and are passed
    * to the boundary routines.
    */
   std::vector<int> d_node_bdry_edge; // Only used for 2D.
   std::vector<int> d_edge_bdry_face; // Only used for 3D.
   std::vector<int> d_node_bdry_face; // Only used for 3D.

   /*
    * Arrays of face (3d) or edge (2d) boundary values for DIRICHLET case.
    */
   std::vector<double> d_bdry_edge_val; // Only used for 2D
   std::vector<double> d_bdry_face_val; // Only used for 3D

};
#endif // HAVE_SUNDIALS
