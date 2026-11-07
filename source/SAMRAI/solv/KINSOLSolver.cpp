/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   C++ Wrapper class for KINSOL solver package
 *
 ************************************************************************/

#include "SAMRAI/solv/KINSOLSolver.h"

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * KINSOLSolver constructor and destructor.
 *
 *************************************************************************
 */
KINSOLSolver::KINSOLSolver(
   const std::string& object_name,
   KINSOLAbstractFunctions* my_functions,
   const int uses_preconditioner,
   const int uses_jac_times_vector):
   d_object_name(object_name),
   d_solution_vector(0),
   d_KINSOL_functions(my_functions),
   d_uses_preconditioner(uses_preconditioner),
   d_uses_jac_times_vector(uses_jac_times_vector),
   d_kin_mem(0),
   d_kinsol_log_file(0),
   d_kinsol_log_file_name("kinsol.log"),
   d_soln_scale(0),
   d_my_soln_scale_vector(false),
   d_fval_scale(0),
   d_my_fval_scale_vector(false),
   d_constraints(0),
   d_linear_solver(0),
   d_KINSOL_needs_initialization(true),
   d_krylov_dimension(15),
   d_max_restarts(0),
   d_max_solves_no_set(10),
   d_max_iter(200),
   d_max_newton_step(-1.0),
   d_global_strategy(KIN_NONE),
   d_residual_tol(-1.0),
   d_step_tol(-1.0),
   d_maxsub(5),
   d_no_initial_setup(0),
   d_no_residual_monitoring(0),
   d_omega_min(0.00001),
   d_omega_max(0.9),
   d_omega(0.0),
   d_no_min_eps(0),
   d_max_beta_fails(10),
   d_eta_choice(KIN_ETACONSTANT),
   d_eta_constant(0.1),
   d_eta_gamma(0.9),
   d_eta_alpha(2.0),
   d_relative_function_error(-1.0),
   d_print_level(0)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(my_functions != 0);
}

void
KINSOLSolver::freeInternalVectors() {

   if (d_my_soln_scale_vector && d_my_fval_scale_vector && d_soln_scale) {
      d_soln_scale->freeVector();
      d_soln_scale = 0;
      d_fval_scale = 0;
      d_my_soln_scale_vector = false;
      d_my_fval_scale_vector = false;
   }

   if (d_my_soln_scale_vector && d_soln_scale) {
      d_soln_scale->freeVector();
      d_soln_scale = 0;
      d_my_soln_scale_vector = false;
   }

   if (d_my_fval_scale_vector && d_fval_scale) {
      d_fval_scale->freeVector();
      d_fval_scale = 0;
      d_my_fval_scale_vector = false;
   }
}

KINSOLSolver::~KINSOLSolver()
{

   freeInternalVectors();

   if (d_kinsol_log_file) {
      fclose(d_kinsol_log_file);
   }

   if (d_kin_mem) {
      KINFree(&d_kin_mem);
   }

   if (d_linear_solver) {
      SUNLinSolFree(d_linear_solver);
   }
}

/*
 *************************************************************************
 *
 * Functions to initialize nonlinear solver and reset KINSOL structure.
 *
 *************************************************************************
 */

void
KINSOLSolver::initialize(
   SundialsAbstractVector* solution,
   SundialsAbstractVector* uscale,
   SundialsAbstractVector* fscale)
{
   TBOX_ASSERT(solution != 0);

   d_solution_vector = solution;

   // Free previously allocated scaling vectors if
   // KINSOLSolver allocated them.
   freeInternalVectors();

   // If user is providing scaling vectors use them
   // otherwise allocate them.
   if (uscale) {
      if (d_my_soln_scale_vector && d_soln_scale) {
         d_soln_scale->freeVector();
      }
      d_soln_scale = uscale;
      d_my_soln_scale_vector = false;
   }

   if (fscale) {
      if (d_my_fval_scale_vector && d_fval_scale) {
         d_fval_scale->freeVector();
      }
      d_fval_scale = fscale;
      d_my_fval_scale_vector = false;
   }

   // Initialize KINSOL.
   d_KINSOL_needs_initialization = true;

   initializeKINSOL();
}

void
KINSOLSolver::initializeKINSOL()
{
   TBOX_ASSERT(d_solution_vector != 0);

   if (d_KINSOL_needs_initialization) {

      if (d_kinsol_log_file) {
         fclose(d_kinsol_log_file);
      }

      d_kinsol_log_file = fopen(d_kinsol_log_file_name.c_str(), "w");

      /*
       * KINSOL function pointers.
       */

      KINSpilsPrecSetupFn precond_set = 0;
      KINSpilsPrecSolveFn precond_solve = 0;
      KINSpilsJacTimesVecFn jac_times_vec = 0;

      if (d_uses_preconditioner) {
         precond_set = KINSOLSolver::KINSOLPrecondSet;
         precond_solve = KINSOLSolver::KINSOLPrecondSolve;
      } else {
         precond_set = 0;
         precond_solve = 0;
      }

      if (d_uses_jac_times_vector) {
         jac_times_vec = KINSOLSolver::KINSOLJacobianTimesVector;
      } else {
         jac_times_vec = 0;
      }

      if (d_kin_mem) KINFree(&d_kin_mem);

      /*
       * Initialize KINSOL structures and set options
       */

      d_kin_mem = KINCreate();

      int ierr = KINInit(d_kin_mem,
            KINSOLSolver::KINSOLFuncEval,
            d_solution_vector->getNVector());
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetUserData(d_kin_mem, this);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetInfoFile(d_kin_mem, d_kinsol_log_file);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaForm(d_kin_mem, d_eta_choice);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaConstValue(d_kin_mem, d_eta_constant);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetEtaParams(d_kin_mem, d_eta_gamma, d_eta_alpha);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxSetupCalls(d_kin_mem, d_max_solves_no_set);
      KINSOL_SAMRAI_ERROR(ierr);

      /*
       * Initialize KINSOL memory record.
       */
      d_linear_solver = SUNSPGMR(d_solution_vector->getNVector(),
                                 PREC_RIGHT,
                                 d_krylov_dimension);

      ierr = KINSpilsSetLinearSolver(d_kin_mem, d_linear_solver);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = SUNSPGMRSetMaxRestarts(d_linear_solver, d_max_restarts);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSpilsSetPreconditioner(d_kin_mem,
            precond_set,
            precond_solve);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSpilsSetJacTimesVecFn(d_kin_mem,
            jac_times_vec);
      KINSOL_SAMRAI_ERROR(ierr);

      if (!(d_residual_tol < 0)) {
         ierr = KINSetFuncNormTol(d_kin_mem, d_residual_tol);
         KINSOL_SAMRAI_ERROR(ierr);
      }

      if (!(d_step_tol < 0)) {
         ierr = KINSetScaledStepTol(d_kin_mem, d_step_tol);
         KINSOL_SAMRAI_ERROR(ierr);
      }

      ierr = KINSetConstraints(d_kin_mem,
            (d_constraints != 0) ? d_constraints->getNVector() : 0);
      KINSOL_SAMRAI_ERROR(ierr);

      // Keep default unless user specifies one.
      if (!(d_max_newton_step < 0)) {
         ierr = KINSetMaxNewtonStep(d_kin_mem, d_max_newton_step);
         KINSOL_SAMRAI_ERROR(ierr);
      }

      if (!(d_relative_function_error < 0)) {
         ierr = KINSetRelErrFunc(d_kin_mem, d_relative_function_error);
         KINSOL_SAMRAI_ERROR(ierr);
      }

      ierr = KINSetPrintLevel(d_kin_mem, d_print_level);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetConstraints(d_kin_mem,
            d_constraints == 0 ? 0 : d_constraints->getNVector());
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNumMaxIters(d_kin_mem, d_max_iter);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoInitSetup(d_kin_mem, d_no_initial_setup);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoResMon(d_kin_mem, d_no_residual_monitoring);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxSubSetupCalls(d_kin_mem, d_maxsub);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetResMonParams(d_kin_mem, d_omega_min, d_omega_max);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetResMonConstValue(d_kin_mem, d_omega);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetNoMinEps(d_kin_mem, d_no_min_eps);
      KINSOL_SAMRAI_ERROR(ierr);

      ierr = KINSetMaxBetaFails(d_kin_mem, d_max_beta_fails);
      KINSOL_SAMRAI_ERROR(ierr);

   } // if no need to initialize KINSOL, function does nothing

   d_KINSOL_needs_initialization = false;
}

/*
 *************************************************************************
 *
 * Solve nonlinear system; re-initialize KINSOL solver, if necessary.
 *
 *************************************************************************
 */
int
KINSOLSolver::solve()
{
   initializeKINSOL();

   /*
    * If scaling vectors are not provided, we make defaults here.
    */
   if (!d_soln_scale) {
      d_soln_scale = d_solution_vector->makeNewVector();
      d_soln_scale->setToScalar(1.0);
      d_my_soln_scale_vector = true;

      if (!d_fval_scale) {
         d_fval_scale = d_soln_scale;
         d_my_fval_scale_vector = true;
      }
   }

   if (!d_fval_scale) {
      d_fval_scale = d_solution_vector->makeNewVector();
      d_fval_scale->setToScalar(1.0);
      d_my_fval_scale_vector = true;
   }

   /*
    * See kinsol.h header file for definition of return types.
    */

   int retval = KINSol(d_kin_mem,
         d_solution_vector->getNVector(),
         d_global_strategy,
         d_soln_scale->getNVector(),
         d_fval_scale->getNVector());

   return retval;

}

/*
 *************************************************************************
 *
 * Setting KINSOL log file name and print flag for KINSOL statistics.
 *
 *************************************************************************
 */

void
KINSOLSolver::setLogFileData(
   const std::string& log_fname,
   const int flag)
{
   TBOX_ASSERT(flag >= 0 && flag <= 3);
   if (!(log_fname == d_kinsol_log_file_name)) {
      if (!log_fname.empty()) {
         d_kinsol_log_file_name = log_fname;
      }
   }
   d_print_level = flag;
   d_KINSOL_needs_initialization = true;
}

/*
 *************************************************************************
 *
 * Print KINSOLSolver object data to given output stream.
 *
 *************************************************************************
 */
void
KINSOLSolver::printClassData(
   std::ostream& os) const
{
   os << "\nKINSOLSolver object data members..." << std::endl;
   os << "this = " << (KINSOLSolver *)this << std::endl;
   os << "d_solution_vector = "
      << (SundialsAbstractVector *)d_solution_vector << std::endl;
   os << "d_soln_scale = "
      << (SundialsAbstractVector *)d_soln_scale << std::endl;
   os << "d_fval_scale = "
      << (SundialsAbstractVector *)d_fval_scale << std::endl;
   os << "d_my_soln_scale_vector = " << d_my_soln_scale_vector << std::endl;
   os << "d_my_fval_scale_vector = " << d_my_fval_scale_vector << std::endl;
   os << "d_constraints = " << (SundialsAbstractVector *)d_constraints
      << std::endl;

   os << "d_KINSOL_functions = "
      << (KINSOLAbstractFunctions *)d_KINSOL_functions << std::endl;

   os << "d_uses_preconditioner = " << d_uses_preconditioner << std::endl;
   os << "d_uses_jac_times_vector = " << d_uses_jac_times_vector << std::endl;

   os << "d_kin_mem = " << d_kin_mem << std::endl;
   os << "d_kinsol_log_file = " << (FILE *)d_kinsol_log_file << std::endl;
   os << "d_kinsol_log_file_name = " << d_kinsol_log_file_name << std::endl;

   os << "d_krylov_dimension = " << d_krylov_dimension << std::endl;
   os << "d_max_restarts = " << d_max_restarts << std::endl;
   os << "d_max_solves_no_set = " << d_max_solves_no_set << std::endl;
   os << "d_global_strategy = " << d_global_strategy << std::endl;
   os << "d_residual_tol = " << d_residual_tol << std::endl;
   os << "d_step_tol = " << d_step_tol << std::endl;

   // SGS add missing output

   os << "...end of KINSOLSolver object data members\n" << std::endl;

}

}
}

#endif
