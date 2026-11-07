/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class containing numerical routines for modified Bratu problem
 *
 ************************************************************************/

#ifndef included_ModifiedBratuProblem
#define included_ModifiedBratuProblem

#include "SAMRAI/SAMRAI_config.h"

#if !defined(HAVE_PETSC) || !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)

/*
 *************************************************************************
 * If the library is not compiled with PETSC -and- KINSOL, print an error.
 * If we're running autotests, skip the error and compile an empty
 * class.
 *************************************************************************
 */
#if (TESTING != 1)
#error \
   "This example requires SAMRAI be compiled with KINSOL, PETSC, and HYPRE."
#endif

#else

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/solv/CellPoissonFACSolver.h"
#include "SAMRAI/solv/KINSOLAbstractFunctions.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/PETScAbstractVectorReal.h"
#include "SAMRAI/solv/SNESAbstractFunctions.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/algs/ImplicitEquationStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <vector>
#include <memory>

using namespace SAMRAI;
using namespace xfer;
using namespace hier;

/**
 * Class ModifiedBratuProblem class provides operations needed to solve
 *
 *     du/dt = div( D(x,t)*grad(u) ) + lambda * exp(u) + f(u,x,t)
 *
 * using implicit time integration and either KINSOL or PETSc to solve the
 * nonlinear system at each step.  Specifically, it provides operations
 * needed by the algs::ImplicitIntegrator class as well as those defined by the
 * interfaces to KINSOL and PETSc; i.e., KINSOLAbstractFunctions and
 * SNESAbstractFunctions respectively.
 *
 * This example is implemented only for 2D, 2:1 refinement ratios only.
 */

class ModifiedBratuProblem:
   public algs::ImplicitEquationStrategy,
   public mesh::StandardTagAndInitStrategy,
   public solv::SNESAbstractFunctions,
   public solv::KINSOLAbstractFunctions,
   public RefinePatchStrategy,
   public CoarsenPatchStrategy,
   public tbox::Serializable
{
public:
   /**
    * Constructor for ModifiedBratuProblem class creates problem variables
    * to represent the solution and other quantities on the patch hierarchy.
    * It initializes data members to default values and sets others based
    * on input and/or restart values.  The constructor also sets up algorithms
    * for communicating data between patches on the hierarchy.
    */
   ModifiedBratuProblem(
      const std::string& object_name,
      const tbox::Dimension& dim,
      const std::shared_ptr<solv::CellPoissonFACSolver> fac_solver,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
      std::shared_ptr<appu::VisItDataWriter> visit_data_writer =
         std::shared_ptr<appu::VisItDataWriter>());

   /**
    * Destructor for ModifiedBratuProblem class does nothing.
    */
   ~ModifiedBratuProblem();

   //@{
   /*!
    * @name Implicit integrator interfaces
    */

   /**
    * Set vector weights on the hierarchy.
    */
   void
   setVectorWeights(
      std::shared_ptr<hier::PatchHierarchy> hierarchy);

   /**
    * Set the nonlinear solution vector so that the new solution data is
    * solved for when the nonlinear solver advances the solution.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   void
   setupSolutionVector(
      const std::shared_ptr<solv::SAMRAIVectorReal<double> >& solution);

   /**
    * Return time increment for advancing the solution at the first timestep.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   double
   getInitialDt();

   /**
    * Return the next time increment through which to advance the solution.
    * The good_solution is the value returned by a call to checkNewSolution(),
    * which determines whether the computed solution is acceptable or not.
    * The integer solver_retcode is the return code generated by the
    * nonlinear solver.   This value must be interpreted in a manner
    * consistant with the solver in use.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   double
   getNextDt(
      const bool good_solution,
      const int solver_retcode);

   /**
    * Set the initial guess for the time advanced solution at the start
    * of the nonlinear iteration.  The boolean argument first_step
    * indicates whether we are at the first step on the current hierarchy
    * configuration.  This is true when the hierarchy is constructed
    * initially and after regridding.  In these cases, setting the initial
    * iterate using extrapolation, for example, may not be possible.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   void
   setInitialGuess(
      const bool first_step,
      const double current_time,
      const double current_dt,
      const double old_dt);

   /**
    * Check the computed solution and return true if it is acceptable;
    * otherwise return false.  The integer solver_retcode is the return
    * code generated by the nonlinear solver.  This value must be
    * interpreted in a manner consistent with the solver in use.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   bool
   checkNewSolution(
      const int solver_retcode);

   /**
    * Update solution storage and dependent quantities after computing an
    * acceptable time advanced solution.   The new_time value is the new
    * solution time.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
   void
   updateSolution(
      const double new_time);

   //@}

   //@{
   /*!
    * @name Functions overloaded from mesh::StandardTagAndInitStrategy.
    */

   virtual void
   initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

   void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);
   //@}

   //@{
   /*!
    * @name Interface functions overloaded from KINSOLAbstractFunctions.
    */

   void
   evaluateNonlinearFunction(
      solv::SundialsAbstractVector* soln,
      solv::SundialsAbstractVector* fval);

   int
   precondSetup(
      solv::SundialsAbstractVector* soln,
      solv::SundialsAbstractVector* soln_scale,
      solv::SundialsAbstractVector* fval,
      solv::SundialsAbstractVector* fval_scale,
      int& num_feval);

   int
   precondSolve(
      solv::SundialsAbstractVector* soln,
      solv::SundialsAbstractVector* soln_scale,
      solv::SundialsAbstractVector* fval,
      solv::SundialsAbstractVector* fval_scale,
      solv::SundialsAbstractVector* rhs,
      int& num_feval);

   int
   jacobianTimesVector(
      solv::SundialsAbstractVector* vector,
      solv::SundialsAbstractVector* product,
      const bool soln_changed,
      solv::SundialsAbstractVector* soln);

   //@}

   //@{
   /*!
    * @name Interface functions overloaded from SNESAbstractFunctions.
    */

   int
   evaluateNonlinearFunction(
      Vec xcur,
      Vec fcur);

   int
   evaluateJacobian(
      Vec x);

   int
   jacobianTimesVector(
      Vec xin,
      Vec xout);

   int
   setupPreconditioner(
      Vec x);

   int
   applyPreconditioner(
      Vec r,
      Vec z);
   //@}

   /*!
    * @brief Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from RefinePatchStrategy.
    */

   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& ghost_width_to_fill);

   //@{
   /*!
    * @name Empty functions for applying user-defined data refine operations
    */

   /*
    * These are overloaded from RefinePatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   void postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   hier::IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

   //@}

   //@{
   /*!
    * @name Empty functions for applying user-defined data coarsen operations
    */
   /*
    * These are overloaded from CoarsenPatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   void postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   hier::IntVector getCoarsenOpStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension& getDim() const
   {
      return d_dim;
   }

   //@}

   /**
    * Write data members to given database for restart.
    *
    * Inherited from tbox::Serializable.
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Write class data to given output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

private:
   /*
    * Functions to read data from input database. If the boolean
    * flag is true, all data members must be present in input.
    *
    * An assertion results if the database pointer is null.
    */
   void
   getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);

   /*
    * Functions for fixing up flux computations along coarse/fine interfaces
    * when ghost cells are filled with CONSTANT_REFINE refinement operators.
    */

   void
   getLevelEdges(
      hier::BoxContainer& boxes,
      std::shared_ptr<hier::Patch> patch,
      std::shared_ptr<hier::PatchLevel> level,
      const tbox::Dimension::dir_t dim,
      const int face);

   void
   correctLevelFlux(
      std::shared_ptr<hier::PatchLevel> level);

   void
   correctPatchFlux(
      std::shared_ptr<hier::PatchLevel> level,
      std::shared_ptr<hier::Patch> patch,
      std::shared_ptr<pdat::CellData<double> > u);

   //@{
   /*!
    * @name Numerical routines specific to modified Bratu problem
    */
   /*
    * These are needed by the nonlinear solvers.
    * They are called by the interface
    * routines after the vectors and other data has been appropriately
    * unwrapped so that these routines are solver-independent.
    */

   void
   evaluateBratuFunction(
      std::shared_ptr<solv::SAMRAIVectorReal<double> > x,
      std::shared_ptr<solv::SAMRAIVectorReal<double> > f);

   /*!
    * @brief Compute A(x)*x.
    *
    * The A(x) used is the one computed in evaluateBratuJacobian()
    * and stored at d_jacobian_a_id and d_jacobian_b_id.
    */
   int
   jacobianTimesVector(
      std::shared_ptr<solv::SAMRAIVectorReal<double> > vector,
      std::shared_ptr<solv::SAMRAIVectorReal<double> > product);

   void
   setupBratuPreconditioner(
      std::shared_ptr<solv::SAMRAIVectorReal<double> > x);

   int
   applyBratuPreconditioner(
      std::shared_ptr<solv::SAMRAIVectorReal<double> > r,
      std::shared_ptr<solv::SAMRAIVectorReal<double> > z);

   /*!
    * @brief Recompute the jacobian A(x).
    *
    * The diagonal of A(x) is placed at d_jacobian_a_id.
    * The off-diagonals are not computed here, because they are
    * independent of x, and it is easier to not explicitly compute them.
    */
   void
   evaluateBratuJacobian(
      std::shared_ptr<solv::SAMRAIVectorReal<double> > x);

   //@}

   /*
    * The object name is used as a handle to databases stored in
    * restart files and for error reporting purposes.
    */
   std::string d_object_name;

   /*
    * Dimension of the problem.
    */
   const tbox::Dimension d_dim;

   /*
    * We cache a pointer to the grid geometry object to set up initial
    * data and set physical boundary conditions.
    */
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   /*
    * Parameters read from input.
    */

   double d_lambda;      // factor multiplying exponential term
   double d_input_dt;    // time increment

   tbox::ResourceAllocator d_allocator;

   /*
    * *hier::Variable data management.
    *
    * Contexts are labels to describe the way variables are used.
    */
   std::shared_ptr<hier::VariableContext> d_current;
   std::shared_ptr<hier::VariableContext> d_new;
   std::shared_ptr<hier::VariableContext> d_scratch;

   /*
    * Variables for the discrete problem; see comments above class constructor.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_solution;
   std::shared_ptr<pdat::CellVariable<double> > d_source_term;
   std::shared_ptr<pdat::CellVariable<double> > d_exponential_term;
   std::shared_ptr<pdat::SideVariable<double> > d_diffusion_coef;
   std::shared_ptr<pdat::SideVariable<double> > d_flux;
   std::shared_ptr<pdat::OutersideVariable<double> > d_coarse_fine_flux;

   /*
    * For storing Jacobian A(x) stuff and computing Jacobian-vector
    * multiply A(x)*v.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_jacobian_a;
   std::shared_ptr<pdat::FaceVariable<double> > d_jacobian_b;
   int d_jacobian_a_id;
   int d_jacobian_b_id;
   hier::ComponentSelector d_jacobian_data;

   /*
    * For storing Jacobian A(x) stuff in setting up and applying
    * the preconditioner A(x)*z=r.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_precond_a;
   std::shared_ptr<pdat::FaceVariable<double> > d_precond_b;
   int d_precond_a_id;
   int d_precond_b_id;
   hier::ComponentSelector d_precond_data;

   int d_soln_scratch_id;
   int d_flux_id;
   int d_coarse_fine_flux_id;
   int d_function_id;

   hier::ComponentSelector d_problem_data;
   hier::ComponentSelector d_new_patch_problem_data;

   hier::IntVector d_nghosts;

   /*
    * The nonlinear solution process requires a solution vector; we cache
    * a pointer to it here.  A variable is used to define weights for the
    * solution vector entries on a composite grid.
    */
   std::shared_ptr<solv::SAMRAIVectorReal<double> > d_solution_vector;

   std::shared_ptr<pdat::CellVariable<double> > d_weight;

   int d_weight_id;

   /*
    * Communication algorithms and schedules used for filling ghost cells
    * and moving data between levels.
    * Schedules stored in arrays are indexed by the destination
    * level number in the transfer.  They are cached to save the cost
    * of generating them multiple times for the same hierarchy
    * configuration.
    */
   RefineAlgorithm d_fill_new_level;
   RefineAlgorithm d_soln_fill;
   std::vector<std::shared_ptr<RefineSchedule> > d_soln_fill_schedule;
   CoarsenAlgorithm d_flux_coarsen;
   std::vector<std::shared_ptr<CoarsenSchedule> > d_flux_coarsen_schedule;
   CoarsenAlgorithm d_soln_coarsen;
   std::vector<std::shared_ptr<CoarsenSchedule> > d_soln_coarsen_schedule;
   CoarsenAlgorithm d_scratch_soln_coarsen;
   std::vector<std::shared_ptr<CoarsenSchedule> > d_scratch_soln_coarsen_schedule;

   std::shared_ptr<RefineOperator> d_soln_refine_op;
   std::shared_ptr<CoarsenOperator> d_soln_coarsen_op;

   /*
    * Current solution time and time increment used in the solution process.
    * New time is current time + current dt.
    */
   double d_current_time;
   double d_new_time;
   double d_current_dt;

   /*
    * Preconditioner and parameters used for Jacobian system.
    *
    * The FAC solver manages the composite grid solution procedure.
    * The Poisson level strategy solves the problem on each level
    * in the hierarchy.
    */
   bool d_use_old_solver;
   std::shared_ptr<solv::CellPoissonFACSolver> d_FAC_solver;

   int d_max_precond_its;
   double d_precond_tol;

   static std::shared_ptr<tbox::Timer> s_copy_timer;
   static std::shared_ptr<tbox::Timer> s_pc_timer;

};

#endif
#endif
