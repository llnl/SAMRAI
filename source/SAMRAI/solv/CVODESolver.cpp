/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   C++ Wrapper class for CVODE solver package
 *
 ************************************************************************/

#include "SAMRAI/solv/CVODESolver.h"

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
namespace solv {

const int CVODESolver::STAT_OUTPUT_BUFFER_SIZE = 256;

/*
 *************************************************************************
 *
 * CVODESolver constructor and destructor.
 *
 *************************************************************************
 */
CVODESolver::CVODESolver(
   const std::string& object_name,
   CVODEAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(my_functions != 0);

   d_object_name = object_name;
   d_cvode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;

   d_solution_vector = 0;

   /*
    * Set default parameters to safe values or to CVODE/CVSpgmr defaults.
    */

   /*
    * CVODE memory record and log file.
    */
   d_cvode_mem = 0;
   d_linear_solver = 0;
   d_cvode_log_file = 0;
   d_cvode_log_file_name = "cvode.log";

   /*
    * ODE parameters.
    */
   d_t_0 = 0.0;
   d_user_t_f = 0.0;
   d_actual_t_f = 0.0;
   d_ic_vector = 0;

   /*
    * ODE integration parameters.
    */

   setLinearMultistepMethod(CV_BDF);
   setRelativeTolerance(0.0);
   setAbsoluteTolerance(0.0);
   d_absolute_tolerance_vector = 0;
   setSteppingMethod(CV_NORMAL);

   d_max_order = -1;
   d_max_num_internal_steps = -1;
   d_max_num_warnings = -1;
   d_init_step_size = -1;
   d_max_step_size = -1;
   d_min_step_size = -1;

   /*
    * CVSpgmr parameters.
    *
    * Note that when the maximum krylov dimension and CVSpgmr
    * tolerance scale factor are set to 0, CVSpgmr uses its
    * internal default values.  These are described in the header for
    * this class.
    */
   setPreconditioningType(PREC_NONE);
   setGramSchmidtType(MODIFIED_GS);
   setMaxKrylovDimension(0);
   setCVSpgmrToleranceScaleFactor(0);

   d_CVODE_needs_initialization = true;
   d_uses_projectionfn = false;
   d_uses_jtimesrhsfn = false;
}

CVODESolver::~CVODESolver()
{
   if (d_cvode_log_file) {
      fclose(d_cvode_log_file);
   }
   if (d_cvode_mem) {
      CVodeFree(&d_cvode_mem);
   }
   if (d_linear_solver) {
      SUNLinSolFree(d_linear_solver);
   }
}

/*
 *************************************************************************
 *
 * Functions to initialize linear solver and reset CVODE structure.
 *
 *************************************************************************
 */

void
CVODESolver::initializeCVODE()
{
   TBOX_ASSERT(d_solution_vector != 0);

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   if (d_CVODE_needs_initialization) {

      /*
       * Set CVODE log file.
       */
      if (d_cvode_log_file) {
         fclose(d_cvode_log_file);
      }
      d_cvode_log_file = fopen(d_cvode_log_file_name.c_str(), "w");

      /*
       * Make sure that either the relative tolerance or the
       * absolute tolerance has been set to a nonzero value.
       */
      bool tolerance_error = false;
      if (d_use_scalar_absolute_tolerance) {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_scalar == 0.0)) {
            tolerance_error = true;
         }
      } else {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_vector->maxNorm() == 0.0)) {
            tolerance_error = true;
         }
      }

      if (tolerance_error && d_cvode_log_file) {
         fprintf(d_cvode_log_file,
            "%s: Both relative and absolute tolerance have value 0.0",
            d_object_name.c_str());
      }

      /*
       * CVODE function pointer.
       */
      CVRhsFn RHSFunc = CVODESolver::CVODERHSFuncEval;

      /*
       * Free previously allocated CVode memory.  Note that the
       * CVReInit() function is not used since the d_neq variable
       * might have been changed from the previous initializeCVODE()
       * call.
       */
      if (d_cvode_mem) CVodeFree(&d_cvode_mem);

      /*
       * Allocate main memory for CVODE package.
       */

      d_cvode_mem = CVodeCreate(d_linear_multistep_method);

      int ierr = CVodeInit(d_cvode_mem,
            RHSFunc,
            d_t_0,
            d_ic_vector->getNVector());
      CVODE_SAMRAI_ERROR(ierr);

      ierr = CVodeSetUserData(d_cvode_mem, this);
      CVODE_SAMRAI_ERROR(ierr);

      ierr = CVodeSStolerances(d_cvode_mem,
                               d_relative_tolerance,
                               d_absolute_tolerance_scalar);

      d_linear_solver = SUNSPGMR(d_solution_vector->getNVector(),
                                 d_precondition_type,
                                 d_max_krylov_dim);

      ierr = CVSpilsSetLinearSolver(d_cvode_mem, d_linear_solver);
      CVODE_SAMRAI_ERROR(ierr);

      if (!(d_max_order < 1)) {
         ierr = CVodeSetMaxOrd(d_cvode_mem, d_max_order);
         CVODE_SAMRAI_ERROR(ierr);
      }

      /*
       * Setup CVSpgmr function pointers.
       */
      CVSpilsPrecSetupFn precond_set = 0;
      CVSpilsPrecSolveFn precond_solve = 0;

      if (d_uses_preconditioner) {
         precond_set = CVODESolver::CVSpgmrPrecondSet;
         precond_solve = CVODESolver::CVSpgmrPrecondSolve;
         CVSpilsSetPreconditioner(d_cvode_mem, precond_set,
            precond_solve);
      }

      if (!(d_max_num_internal_steps < 0)) {
         ierr = CVodeSetMaxNumSteps(d_cvode_mem, d_max_num_internal_steps);
         CVODE_SAMRAI_ERROR(ierr);
      }

      if (!(d_max_num_warnings < 0)) {
         ierr = CVodeSetMaxHnilWarns(d_cvode_mem, d_max_num_warnings);
         CVODE_SAMRAI_ERROR(ierr);
      }

      if (!(d_init_step_size < 0)) {
         ierr = CVodeSetInitStep(d_cvode_mem, d_init_step_size);
         CVODE_SAMRAI_ERROR(ierr);
      }

      if (!(d_max_step_size < 0)) {
         ierr = CVodeSetMaxStep(d_cvode_mem, d_max_step_size);
         CVODE_SAMRAI_ERROR(ierr);
      }

      if (!(d_min_step_size < 0)) {
         ierr = CVodeSetMinStep(d_cvode_mem, d_min_step_size);
         CVODE_SAMRAI_ERROR(ierr);
      }

      if (d_uses_projectionfn) {
         CVProjFn proj_fn = CVODESolver::CVODEProjEval;
         ierr = CVodeSetProjFn(d_cvode_mem , proj_fn);
      }

      if (d_uses_jtimesrhsfn) {
         CVRhsFn jtimesrhs_fn = CVODESolver::CVODEJTimesRHSFuncEval;
         ierr = CVodeSetJacTimesRhsFn(d_cvode_mem , jtimesrhs_fn);
      }

   } // if no need to initialize CVODE, function does nothing

   d_CVODE_needs_initialization = false;
}

/*
 *************************************************************************
 *
 * Access methods for CVODE statistics.
 *
 *************************************************************************
 */

void
CVODESolver::printCVODEStatistics(
   std::ostream& os) const
{

   char buf[STAT_OUTPUT_BUFFER_SIZE];

   os << "\nCVODESolver: CVODE statistics... " << std::endl;

   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "lenrw           = %5d     leniw            = %5d\n",
      getCVODEMemoryUsageForDoubles(),
      getCVODEMemoryUsageForIntegers());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "nst             = %5d     nfe              = %5d\n",
      getNumberOfInternalStepsTaken(),
      getNumberOfRHSFunctionCalls());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "nni             = %5d     nsetups          = %5d\n",
      getNumberOfNewtonIterations(),
      getNumberOfLinearSolverSetupCalls());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "netf            = %5d     ncfn             = %5d\n",
      getNumberOfLocalErrorTestFailures(),
      getNumberOfNonlinearConvergenceFailures());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "qu              = %5d     qcur             = %5d\n",
      getOrderUsedDuringLastInternalStep(),
      getOrderToBeUsedDuringNextInternalStep());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "\nhu              = %e      hcur             = %e\n",
      getStepSizeForLastInternalStep(),
      getStepSizeForNextInternalStep());
   os << buf;
   snprintf(buf, STAT_OUTPUT_BUFFER_SIZE, "tcur            = %e      tolsf            = %e\n",
      getCurrentInternalValueOfIndependentVariable(),
      getCVODESuggestedToleranceScalingFactor());
   os << buf;
}

/*
 *************************************************************************
 *
 * Access methods for CVSpgmr statistics.
 *
 *************************************************************************
 */

void
CVODESolver::printCVSpgmrStatistics(
   std::ostream& os) const
{
   os << "CVODESolver: CVSpgmr statistics... " << std::endl;

   os << "spgmr_lrw       = "
      << tbox::Utilities::intToString(getCVSpgmrMemoryUsageForDoubles(), 5)
      << "     spgmr_liw        = "
      << tbox::Utilities::intToString(getCVSpgmrMemoryUsageForIntegers(),
      5) << std::endl;

   os << "nli             = "
      << tbox::Utilities::intToString(getNumberOfLinearIterations(), 5)
      << "     ncfl             = "
      << tbox::Utilities::intToString(
      getNumberOfLinearConvergenceFailures(), 5) << std::endl;

   os << "npe             = "
      << tbox::Utilities::intToString(
      getNumberOfPreconditionerEvaluations(), 5)
      << "     nps              = "
      << tbox::Utilities::intToString(getNumberOfPrecondSolveCalls(),
      5) << std::endl;
}

/*
 *************************************************************************
 *
 * Print CVODESolver object data to given output stream.
 *
 *************************************************************************
 */
void
CVODESolver::printClassData(
   std::ostream& os) const
{
   os << "\nCVODESolver object data members..." << std::endl;
   os << "Object name = "
      << d_object_name << std::endl;

   os << "this = " << (CVODESolver *)this << std::endl;
   os << "d_solution_vector = "
      << (SundialsAbstractVector *)d_solution_vector << std::endl;

   os << "d_CVODE_functions = "
      << (CVODEAbstractFunctions *)d_cvode_functions << std::endl;

   os << "&d_cvode_mem = " << d_cvode_mem << std::endl;
   os << "d_cvode_log_file = " << (FILE *)d_cvode_log_file << std::endl;
   os << "d_cvode_log_file_name = " << d_cvode_log_file_name << std::endl;

   os << std::endl;
   os << "CVODE parameters..." << std::endl;
   os << "d_t_0 = "
      << d_t_0 << std::endl;
   os << "d_ic_vector = "
      << (SundialsAbstractVector *)d_ic_vector << std::endl;

   os << "d_linear_multistep_method = "
      << d_linear_multistep_method << std::endl;
   os << "d_relative_tolerance = "
      << d_relative_tolerance << std::endl;
   os << "d_use_scalar_absolute_tolerance = ";
   if (d_use_scalar_absolute_tolerance) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }
   os << "d_absolute_tolerance_scalar = "
      << d_absolute_tolerance_scalar << std::endl;
   os << "d_absolute_tolerance_vector= " << std::endl;
   d_absolute_tolerance_vector->printVector();

   os << "Optional CVODE inputs (see CVODE docs for details):"
      << std::endl;

   os << "maximum linear multistep method order = "
      << d_max_order << std::endl;
   os << "maximum number of internal steps = "
      << d_max_num_internal_steps << std::endl;
   os << "maximum number of nil internal step warnings = "
      << d_max_num_warnings << std::endl;

   os << "initial step size = "
      << d_init_step_size << std::endl;
   os << "maximum absolute value of step size = "
      << d_max_step_size << std::endl;
   os << "minimum absolute value of step size = "
      << d_min_step_size << std::endl;
   os << "last step size = "
      << getStepSizeForLastInternalStep() << std::endl;
   os << "...end of CVODE parameters\n" << std::endl;

   os << std::endl;
   os << "CVSpgmr parameters..." << std::endl;
   os << "d_precondition_type = "
      << d_precondition_type << std::endl;
   os << "d_gram_schmidt_type = "
      << d_gram_schmidt_type << std::endl;
   os << "d_max_krylov_dim = "
      << d_max_krylov_dim << std::endl;
   os << "d_tol_scale_factor = "
      << d_tol_scale_factor << std::endl;
   os << "...end of CVSpgmr parameters\n" << std::endl;

   os << "d_CVODE_needs_initialization = ";
   if (d_CVODE_needs_initialization) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }

   os << "...end of CVODESolver object data members\n" << std::endl;
}

}
}

#endif
