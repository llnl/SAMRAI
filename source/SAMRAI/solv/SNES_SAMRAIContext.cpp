/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Wrapper for SNES solver for use in a SAMRAI-based application.
 *
 ************************************************************************/
#include "SAMRAI/solv/SNES_SAMRAIContext.h"

#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/tbox/RestartManager.h"

#ifdef HAVE_PETSC

namespace SAMRAI {
namespace solv {

const int SNES_SAMRAIContext::SOLV_SNES_SAMRAI_CONTEXT_VERSION = 1;

/*
 *************************************************************************
 *
 * Static member functions that provide linkage with PETSc/SNES package.
 * See header file for SNESAbstractFunctions for more information.
 *
 *************************************************************************
 */

int
SNES_SAMRAIContext::SNESJacobianSet(
   SNES snes,
   Vec x,
   Mat A,
   Mat B,
   void* ctx)
{
   NULL_USE(snes);
   NULL_USE(B);
   int retval = 0;
   if (((SNES_SAMRAIContext *)ctx)->getUsesExplicitJacobian()) {
      retval =
         ((SNES_SAMRAIContext *)ctx)->getSNESFunctions()->
         evaluateJacobian(x);
   } else {
      int ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      PETSC_SAMRAI_ERROR(ierr);
      ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
      PETSC_SAMRAI_ERROR(ierr);
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Constructor and destructor for SNES_SAMRAIContext.  The
 * constructor sets default values for data members, then overrides
 * them with values read from input or restart.  The destructor destroys
 * the SNES object.
 *
 *************************************************************************
 */
SNES_SAMRAIContext::SNES_SAMRAIContext(
   const std::string& object_name,
   SNESAbstractFunctions* my_functions,
   const std::shared_ptr<tbox::Database>& input_db):
   d_object_name(object_name),
   d_context_needs_initialization(true),
   d_SNES_solver(0),
   d_krylov_solver(0),
   d_jacobian(0),
   d_preconditioner(0),
   d_solution_vector(0),
   d_residual_vector(0),
   d_SNES_functions(my_functions),
   d_uses_preconditioner(true),
   d_uses_explicit_jacobian(true),
   d_maximum_nonlinear_iterations(PETSC_DEFAULT),
   d_maximum_function_evals(PETSC_DEFAULT),
   d_absolute_tolerance(PETSC_DEFAULT),
   d_relative_tolerance(PETSC_DEFAULT),
   d_step_tolerance(PETSC_DEFAULT),
   d_forcing_term_strategy("CONSTANT"),
   d_forcing_term_flag(PETSC_DEFAULT),
   d_constant_forcing_term(PETSC_DEFAULT),
   d_initial_forcing_term(PETSC_DEFAULT),
   d_maximum_forcing_term(PETSC_DEFAULT),
   d_EW_choice2_alpha(PETSC_DEFAULT),
   d_EW_choice2_gamma(PETSC_DEFAULT),
   d_EW_safeguard_exponent(PETSC_DEFAULT),
   d_EW_safeguard_disable_threshold(PETSC_DEFAULT),
   d_SNES_completion_code(SNES_CONVERGED_ITERATING),
   d_linear_solver_absolute_tolerance(PETSC_DEFAULT),
   d_linear_solver_divergence_tolerance(PETSC_DEFAULT),
   d_maximum_linear_iterations(PETSC_DEFAULT),
   d_maximum_gmres_krylov_dimension(PETSC_DEFAULT),
   d_differencing_parameter_strategy(MATMFFD_WP),
   d_function_evaluation_error(PETSC_DEFAULT),
   d_nonlinear_iterations(0)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(my_functions != 0);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
      this);

   /*
    * Initialize members with data read from the input and restart
    * databases.  Note that PETSc object parameters are set in
    * initialize().
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

}

SNES_SAMRAIContext::~SNES_SAMRAIContext()
{
   if (d_solution_vector) {
      PETSc_SAMRAIVectorReal<double>::destroyPETScVector(
         d_solution_vector);
   }

   if (d_residual_vector) {
      PETSc_SAMRAIVectorReal<double>::destroyPETScVector(
         d_residual_vector);
   }

   destroyPetscObjects();
}

/*
 *************************************************************************
 *
 * Routines to initialize PETSc/SNES solver and solve nonlinear system.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::initialize(
   const std::shared_ptr<SAMRAIVectorReal<double> >& solution)
{
   TBOX_ASSERT(solution);

   /*
    * Set up vectors for solution and nonlinear residual.
    */

   d_solution_vector =
      PETSc_SAMRAIVectorReal<double>::createPETScVector(solution);

   std::shared_ptr<SAMRAIVectorReal<double> > residual(
      solution->cloneVector("residual"));
   residual->allocateVectorData();
   d_residual_vector =
      PETSc_SAMRAIVectorReal<double>::createPETScVector(residual);

   createPetscObjects();
   initializePetscObjects();
}

/*
 *************************************************************************
 *
 * Reset the state of the nonlinear solver.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::resetSolver(
   const int coarsest_level,
   const int finest_level)
{
   std::shared_ptr<SAMRAIVectorReal<double> > solution_vector(
      PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(d_solution_vector));
   solution_vector->deallocateVectorData();
   solution_vector->resetLevels(coarsest_level, finest_level);
   solution_vector->allocateVectorData();

   std::shared_ptr<SAMRAIVectorReal<double> > residual_vector(
      PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(d_residual_vector));
   residual_vector->deallocateVectorData();
   residual_vector->resetLevels(coarsest_level, finest_level);
   residual_vector->allocateVectorData();

   destroyPetscObjects();
   createPetscObjects();
   initializePetscObjects();
}

/*
 *************************************************************************
 *
 * Solve the nonlinear system.
 *
 *************************************************************************
 */
int
SNES_SAMRAIContext::solve()
{
   int ierr;

   if (d_context_needs_initialization) initializePetscObjects();

   Vec initial_guess;

   ierr = VecDuplicate(d_solution_vector, &initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecSet(initial_guess, 0.0);
   PETSC_SAMRAI_ERROR(ierr);
   ierr = SNESSolve(d_SNES_solver,
         initial_guess,
         d_solution_vector);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetIterationNumber(d_SNES_solver,
         &d_nonlinear_iterations);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetConvergedReason(d_SNES_solver,
         &d_SNES_completion_code);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = VecDestroy(&initial_guess);
   PETSC_SAMRAI_ERROR(ierr);

   return ((int)d_SNES_completion_code > 0) ? 1 : 0;
}

/*
 *************************************************************************
 *
 *  Report the reason for termination of nonlinear iterations.  SNES
 *  return codes are translated here, and a message is placed in the
 *  specified output stream.  Test only on relevant completion codes.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::reportCompletionCode(
   std::ostream& os) const
{
   switch ((int)d_SNES_completion_code) {
      case SNES_CONVERGED_FNORM_ABS:
         os << " Fnorm less than specified absolute tolerance.\n";
         break;
      case SNES_CONVERGED_FNORM_RELATIVE:
         os << " Fnorm less than specified relative tolerance.\n";
         break;
      case SNES_CONVERGED_SNORM_RELATIVE:
         os << " Step size less than specified tolerance.\n";
         break;
      case SNES_DIVERGED_FUNCTION_COUNT:
         os << " Maximum function evaluation count exceeded.\n";
         break;
      case SNES_DIVERGED_FNORM_NAN:
         os << " Norm of F is NAN.\n";
         break;
      case SNES_DIVERGED_MAX_IT:
         os << " Maximum nonlinear iteration count exceeded.\n";
         break;
      case SNES_DIVERGED_LINE_SEARCH:
         os << " Failure in linesearch procedure.\n";
         break;
      default:
         os << " Inappropriate completion code reported.\n";
         break;
   }
}

/*
 *************************************************************************
 *
 * Create needed Petsc objects and cache a pointer to them.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::createPetscObjects()
{
   int ierr = 0;
   NULL_USE(ierr);

   /*
    * Create the nonlinear solver, specify linesearch backtracking,
    * and register method for nonlinear residual evaluation.
    */
   ierr = SNESCreate(PETSC_COMM_SELF, &d_SNES_solver);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetType(d_SNES_solver, SNESNEWTONLS);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESSetFunction(d_SNES_solver,
         d_residual_vector,
         SNES_SAMRAIContext::SNESFuncEval,
         (void *)this);
   PETSC_SAMRAI_ERROR(ierr);
   /*
    * Cache the linear solver object, as well as the wrapped Krylov
    * solver and preconditioner.
    */
//   ierr = SNESGetSLES(d_SNES_solver,
//                      &d_SLES_solver);
//                      PETSC_SAMRAI_ERROR(ierr);

//   ierr = SLESGetKSP(d_SLES_solver,
//                     &d_krylov_solver);
//                     PETSC_SAMRAI_ERROR(ierr);

   ierr = SNESGetKSP(d_SNES_solver, &d_krylov_solver);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPSetPCSide(d_krylov_solver, PC_RIGHT);
   PETSC_SAMRAI_ERROR(ierr);

   ierr = KSPGetPC(d_krylov_solver, &d_preconditioner);
   PETSC_SAMRAI_ERROR(ierr);

}

/*
 *************************************************************************
 *
 * Initialize the state of cached Petsc objects from cached information.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::initializePetscObjects()
{
   int ierr = 0;
   NULL_USE(ierr);

   /*
    * Set tolerances in nonlinear solver.  Also set parameters if
    * the Jacobian-free option has been selected.
    */
   ierr = SNESSetTolerances(d_SNES_solver,
         d_absolute_tolerance,
         d_relative_tolerance,
         d_step_tolerance,
         d_maximum_nonlinear_iterations,
         d_maximum_function_evals);
   PETSC_SAMRAI_ERROR(ierr);

   if (!(d_forcing_term_strategy == "CONSTANT")) {

      ierr = SNESKSPSetUseEW(d_SNES_solver, PETSC_TRUE);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = SNESKSPSetParametersEW(d_SNES_solver,
            d_forcing_term_flag,
            d_initial_forcing_term,
            d_maximum_forcing_term,
            d_EW_choice2_gamma,
            d_EW_choice2_alpha,
            d_EW_safeguard_exponent,
            d_EW_safeguard_disable_threshold);
      PETSC_SAMRAI_ERROR(ierr);
   }

   /*
    * Create data structures needed for Jacobian.  This is done
    * here in case an application toggles use of an explicit
    * Jacobian within a run.
    *
    * First delete any Jacobian object that already has been created.
    */
   if (d_jacobian) {
      MatDestroy(&d_jacobian);
   }
   if (d_uses_explicit_jacobian) {

      ierr = MatCreateShell(PETSC_COMM_SELF,
            0,                   // dummy number of local rows
            0,                   // dummy number of local columns
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            (void *)this,
            &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = MatShellSetOperation(d_jacobian,
            MATOP_MULT,
            (void (*)()) SNES_SAMRAIContext::
            SNESJacobianTimesVector);
      PETSC_SAMRAI_ERROR(ierr);

   } else {

      ierr = MatCreateSNESMF(d_SNES_solver,
            &d_jacobian);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = MatMFFDSetType(
            d_jacobian,
            (MatMFFDType)d_differencing_parameter_strategy.c_str());

      ierr = MatMFFDSetFunctionError(d_jacobian,
            d_function_evaluation_error);
   }

   /*
    * Register method for setting up Jacobian; this is the same
    * for both options.
    *
    * N.B.  In principle, the second Mat argument should not
    * be the same as the first Mat argument.  However we
    * restrict to either no preconditioner, or a shell
    * preconditioner; in these circumstances that seems to
    * cause no problem, since the shell preconditioner provides
    * its own setup method.
    */
   ierr = SNESSetJacobian(d_SNES_solver,
         d_jacobian,
         d_jacobian,
         SNES_SAMRAIContext::SNESJacobianSet,
         (void *)this);
   PETSC_SAMRAI_ERROR(ierr);

   /*
    * Initialize the Krylov solver object.  This includes setting the
    * type of Krylov method that is used and tolerances used by the
    * method.
    */
   ierr = KSPSetType(d_krylov_solver, (KSPType)d_linear_solver_type.c_str());
   PETSC_SAMRAI_ERROR(ierr);

   if (d_linear_solver_type == "gmres") {

      ierr = KSPGMRESSetRestart(
            d_krylov_solver,
            d_maximum_gmres_krylov_dimension);
      PETSC_SAMRAI_ERROR(ierr);

      if (d_gmres_orthogonalization_algorithm == "modifiedgramschmidt") {

         ierr = KSPGMRESSetOrthogonalization(
               d_krylov_solver,
               KSPGMRESModifiedGramSchmidtOrthogonalization);
         PETSC_SAMRAI_ERROR(ierr);

      } else if (d_gmres_orthogonalization_algorithm ==
                 "gmres_cgs_refine_ifneeded") {

         ierr = KSPGMRESSetCGSRefinementType(
               d_krylov_solver,
               KSP_GMRES_CGS_REFINE_IFNEEDED);
         PETSC_SAMRAI_ERROR(ierr);
      } else if (d_gmres_orthogonalization_algorithm ==
                 "gmres_cgs_refine_always") {

         ierr = KSPGMRESSetCGSRefinementType(
               d_krylov_solver,
               KSP_GMRES_CGS_REFINE_ALWAYS);
         PETSC_SAMRAI_ERROR(ierr);
      }
   }

   if (d_forcing_term_strategy == "CONSTANT") {

      ierr = KSPSetTolerances(d_krylov_solver,
            d_constant_forcing_term,
            d_linear_solver_absolute_tolerance,
            d_linear_solver_divergence_tolerance,
            d_maximum_linear_iterations);
      PETSC_SAMRAI_ERROR(ierr);
   }

   /*
    * Initialize the precondtioner.  Only shell PCs are supported.
    * For these, register the methods used to set up and apply
    * the preconditioner.
    */
   if (d_uses_preconditioner) {

      std::string pc_type = "shell";
      ierr = PCSetType(d_preconditioner, (PCType)pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetSetUp(d_preconditioner,
            SNES_SAMRAIContext::SNESsetupPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetContext(d_preconditioner, this);
      PETSC_SAMRAI_ERROR(ierr);

      ierr = PCShellSetApply(d_preconditioner,
            SNES_SAMRAIContext::SNESapplyPreconditioner);
      PETSC_SAMRAI_ERROR(ierr);

   } else {

      std::string pc_type = "none";
      ierr = PCSetType(d_preconditioner, (PCType)pc_type.c_str());
      PETSC_SAMRAI_ERROR(ierr);

   }

   d_context_needs_initialization = false;
}

/*
 *************************************************************************
 *
 * Destroy cached Petsc objects.
 *
 *************************************************************************
 */
void
SNES_SAMRAIContext::destroyPetscObjects()
{
   if (d_jacobian) {
      MatDestroy(&d_jacobian);
      d_jacobian = 0;
   }

   if (d_SNES_solver) {
      SNESDestroy(&d_SNES_solver);
//     if (d_SLES_solver) d_SLES_solver = 0;
      if (d_preconditioner) d_preconditioner = 0;
      if (d_krylov_solver) d_krylov_solver = 0;
   }
}

/*
 *************************************************************************
 *
 * Read parameters from input that are cached in this object.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (input_db) {
      if (!is_from_restart) {

         d_maximum_nonlinear_iterations =
            input_db->getIntegerWithDefault("maximum_nonlinear_iterations", PETSC_DEFAULT);

         d_maximum_function_evals =
            input_db->getIntegerWithDefault("maximum_function_evals", PETSC_DEFAULT);

         d_uses_preconditioner =
            input_db->getBoolWithDefault("uses_preconditioner", true);

         d_uses_explicit_jacobian =
            input_db->getBoolWithDefault("uses_explicit_jacobian", true);

         d_absolute_tolerance =
            input_db->getDoubleWithDefault("absolute_tolerance", PETSC_DEFAULT);

         d_relative_tolerance =
            input_db->getDoubleWithDefault("relative_tolerance", PETSC_DEFAULT);

         d_step_tolerance =
            input_db->getDoubleWithDefault("step_tolerance", PETSC_DEFAULT);

         d_forcing_term_strategy =
            input_db->getStringWithDefault("forcing_term_strategy", "CONSTANT");
         if (d_forcing_term_strategy == "EWCHOICE1") {
            d_forcing_term_flag = 1;
         } else if (d_forcing_term_strategy == "EWCHOICE2") {
            d_forcing_term_flag = 2;
         } else if (!(d_forcing_term_strategy == "CONSTANT")) {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Key data `forcing_term_strategy' = "
                             << d_forcing_term_strategy
                             << " in input not recognized.");
         }

         d_constant_forcing_term =
            input_db->getDoubleWithDefault("constant_forcing_term", PETSC_DEFAULT);

         d_initial_forcing_term =
            input_db->getDoubleWithDefault("initial_forcing_term", PETSC_DEFAULT);

         d_maximum_forcing_term =
            input_db->getDoubleWithDefault("maximum_forcing_term", PETSC_DEFAULT);

         d_EW_choice2_alpha =
            input_db->getDoubleWithDefault("EW_choice2_alpha", PETSC_DEFAULT);

         d_EW_choice2_gamma =
            input_db->getDoubleWithDefault("EW_choice2_gamma", PETSC_DEFAULT);

         d_EW_safeguard_exponent =
            input_db->getDoubleWithDefault("EW_safeguard_exponent", PETSC_DEFAULT);

         d_EW_safeguard_disable_threshold =
            input_db->getDoubleWithDefault("EW_safeguard_disable_threshold", PETSC_DEFAULT);

         d_linear_solver_type =
            input_db->getStringWithDefault("linear_solver_type", "");

         d_linear_solver_absolute_tolerance =
            input_db->getDoubleWithDefault("linear_solver_absolute_tolerance", PETSC_DEFAULT);

         d_linear_solver_divergence_tolerance =
            input_db->getDoubleWithDefault("linear_solver_divergence_tolerance", PETSC_DEFAULT);

         d_maximum_linear_iterations =
            input_db->getIntegerWithDefault("maximum_linear_iterations", PETSC_DEFAULT);

         d_maximum_gmres_krylov_dimension =
            input_db->getIntegerWithDefault("maximum_gmres_krylov_dimension", PETSC_DEFAULT);

         d_gmres_orthogonalization_algorithm =
            input_db->getStringWithDefault("gmres_orthogonalization_algorithm", "");

         d_differencing_parameter_strategy =
            input_db->getStringWithDefault("differencing_parameter_strategy", MATMFFD_WP);
         if (!(d_differencing_parameter_strategy == MATMFFD_WP ||
               d_differencing_parameter_strategy == MATMFFD_DS)) {
            INPUT_VALUE_ERROR("differencing_parameter_strategy");
         }

         d_function_evaluation_error =
            input_db->getDoubleWithDefault("function_evaluation_error", PETSC_DEFAULT);
      } else {
         bool read_on_restart =
            input_db->getBoolWithDefault("read_on_restart", false);
         if (!read_on_restart) {
            return;
         }

         d_maximum_nonlinear_iterations =
            input_db->getIntegerWithDefault("maximum_nonlinear_iterations",
               d_maximum_nonlinear_iterations);
         d_maximum_function_evals =
            input_db->getIntegerWithDefault("maximum_function_evals",
               d_maximum_function_evals);
         d_uses_preconditioner =
            input_db->getBoolWithDefault("uses_preconditioner",
               d_uses_preconditioner);
         d_uses_explicit_jacobian =
            input_db->getBoolWithDefault("uses_explicit_jacobian",
               d_uses_explicit_jacobian);
         d_absolute_tolerance =
            input_db->getDoubleWithDefault("absolute_tolerance",
               d_absolute_tolerance);
         d_relative_tolerance =
            input_db->getDoubleWithDefault("relative_tolerance",
               d_relative_tolerance);
         d_step_tolerance =
            input_db->getDoubleWithDefault("step_tolerance", d_step_tolerance);
         d_forcing_term_strategy =
            input_db->getStringWithDefault("forcing_term_strategy",
               d_forcing_term_strategy);
         if (d_forcing_term_strategy == "EWCHOICE1") {
            d_forcing_term_flag = 1;
         } else if (d_forcing_term_strategy == "EWCHOICE2") {
            d_forcing_term_flag = 2;
         } else if (!(d_forcing_term_strategy == "CONSTANT")) {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Key data `forcing_term_strategy' = "
                             << d_forcing_term_strategy
                             << " in input not recognized.");
         }
         d_constant_forcing_term =
            input_db->getDoubleWithDefault("constant_forcing_term",
               d_constant_forcing_term);
         d_initial_forcing_term =
            input_db->getDoubleWithDefault("initial_forcing_term",
               d_initial_forcing_term);
         d_maximum_forcing_term =
            input_db->getDoubleWithDefault("maximum_forcing_term",
               d_maximum_forcing_term);
         d_EW_choice2_alpha =
            input_db->getDoubleWithDefault("EW_choice2_alpha",
               d_EW_choice2_alpha);
         d_EW_choice2_gamma =
            input_db->getDoubleWithDefault("EW_choice2_gamma",
               d_EW_choice2_gamma);
         d_EW_safeguard_exponent =
            input_db->getDoubleWithDefault("EW_safeguard_exponent",
               d_EW_safeguard_exponent);
         d_EW_safeguard_disable_threshold =
            input_db->getDoubleWithDefault("EW_safeguard_disable_threshold",
               d_EW_safeguard_disable_threshold);
         d_linear_solver_type =
            input_db->getStringWithDefault("linear_solver_type",
               d_linear_solver_type);
         d_linear_solver_absolute_tolerance =
            input_db->getDoubleWithDefault("linear_solver_absolute_tolerance",
               d_linear_solver_absolute_tolerance);
         d_linear_solver_divergence_tolerance =
            input_db->getDoubleWithDefault("linear_solver_divergence_tolerance",
               d_linear_solver_divergence_tolerance);
         d_maximum_linear_iterations =
            input_db->getIntegerWithDefault("maximum_linear_iterations",
               d_maximum_linear_iterations);
         d_maximum_gmres_krylov_dimension =
            input_db->getIntegerWithDefault("maximum_gmres_krylov_dimension",
               d_maximum_gmres_krylov_dimension);
         d_gmres_orthogonalization_algorithm =
            input_db->getStringWithDefault("gmres_orthogonalization_algorithm",
               d_gmres_orthogonalization_algorithm);
         d_differencing_parameter_strategy =
            input_db->getStringWithDefault("differencing_parameter_strategy",
               d_differencing_parameter_strategy);
         if (d_differencing_parameter_strategy != MATMFFD_WP &&
             d_differencing_parameter_strategy != MATMFFD_DS) {
            TBOX_ERROR("SNES_SAMRAIContext::getFromInput error...\n"
               << "differencing_parameter_strategy must be \"wp\" or \"ds\"."
               << std::endl);
         }
         d_function_evaluation_error =
            input_db->getDoubleWithDefault("function_evaluation_error",
               d_function_evaluation_error);
      }
   }

}

/*
 *************************************************************************
 *
 * Routines to read/write from/to restart/database.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION");
   if (ver != SOLV_SNES_SAMRAI_CONTEXT_VERSION) {
      TBOX_ERROR(d_object_name << ":SNES_SAMRAIContext::getFromRestart() error ...\n "
                               << "Restart file version different "
                               << "than class version.");
   }

   d_uses_preconditioner = db->getBool("uses_preconditioner");
   d_uses_explicit_jacobian = db->getBool("uses_explicit_jacobian");

   d_maximum_nonlinear_iterations =
      db->getInteger("maximum_nonlinear_iterations");
   d_maximum_function_evals = db->getInteger("maximum_function_evals");

   d_absolute_tolerance = db->getDouble("absolute_tolerance");
   d_relative_tolerance = db->getDouble("relative_tolerance");
   d_step_tolerance = db->getDouble("step_tolerance");

   d_forcing_term_strategy = db->getString("forcing_term_strategy");
   d_forcing_term_flag = db->getInteger("d_forcing_term_flag");

   d_constant_forcing_term = db->getDouble("constant_forcing_term");
   d_initial_forcing_term = db->getDouble("initial_forcing_term");
   d_maximum_forcing_term = db->getDouble("maximum_forcing_term");
   d_EW_choice2_alpha = db->getDouble("EW_choice2_alpha");
   d_EW_choice2_gamma = db->getDouble("EW_choice2_gamma");
   d_EW_safeguard_exponent = db->getDouble("EW_safeguard_exponent");
   d_EW_safeguard_disable_threshold =
      db->getDouble("EW_safeguard_disable_threshold");

   d_linear_solver_type = db->getString("linear_solver_type");
   d_linear_solver_absolute_tolerance =
      db->getDouble("linear_solver_absolute_tolerance");
   d_linear_solver_divergence_tolerance =
      db->getDouble("linear_solver_divergence_tolerance");
   d_maximum_linear_iterations =
      db->getInteger("maximum_linear_iterations");

   d_maximum_gmres_krylov_dimension =
      db->getInteger("maximum_gmres_krylov_dimension");
   d_gmres_orthogonalization_algorithm =
      db->getString("gmres_orthogonalization_algorithm");

   d_function_evaluation_error = db->getDouble("function_evaluation_error");
   d_differencing_parameter_strategy =
      db->getString("differencing_parameter_strategy");

}

void
SNES_SAMRAIContext::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("SOLV_SNES_SAMRAI_CONTEXT_VERSION",
      SOLV_SNES_SAMRAI_CONTEXT_VERSION);

   restart_db->putBool("uses_preconditioner", d_uses_preconditioner);
   restart_db->putBool("uses_explicit_jacobian", d_uses_explicit_jacobian);

   restart_db->putInteger("maximum_nonlinear_iterations",
      d_maximum_nonlinear_iterations);
   restart_db->putInteger("maximum_function_evals",
      d_maximum_function_evals);

   restart_db->putDouble("absolute_tolerance", d_absolute_tolerance);
   restart_db->putDouble("relative_tolerance", d_relative_tolerance);
   restart_db->putDouble("step_tolerance", d_step_tolerance);

   restart_db->putString("forcing_term_strategy", d_forcing_term_strategy);
   restart_db->putInteger("d_forcing_term_flag", d_forcing_term_flag);

   restart_db->putDouble("constant_forcing_term", d_constant_forcing_term);
   restart_db->putDouble("initial_forcing_term", d_initial_forcing_term);
   restart_db->putDouble("maximum_forcing_term", d_maximum_forcing_term);
   restart_db->putDouble("EW_choice2_alpha", d_EW_choice2_alpha);
   restart_db->putDouble("EW_choice2_gamma", d_EW_choice2_gamma);
   restart_db->putDouble("EW_safeguard_exponent", d_EW_safeguard_exponent);
   restart_db->putDouble("EW_safeguard_disable_threshold",
      d_EW_safeguard_disable_threshold);

   restart_db->putString("linear_solver_type", d_linear_solver_type);
   restart_db->putDouble("linear_solver_absolute_tolerance",
      d_linear_solver_absolute_tolerance);
   restart_db->putDouble("linear_solver_divergence_tolerance",
      d_linear_solver_divergence_tolerance);
   restart_db->putInteger("maximum_linear_iterations",
      d_maximum_linear_iterations);

   restart_db->putInteger("maximum_gmres_krylov_dimension",
      d_maximum_gmres_krylov_dimension);
   restart_db->putString("gmres_orthogonalization_algorithm",
      d_gmres_orthogonalization_algorithm);

   restart_db->putDouble("function_evaluation_error",
      d_function_evaluation_error);
   restart_db->putString("differencing_parameter_strategy",
      d_differencing_parameter_strategy);

}

/*
 *************************************************************************
 *
 * Write all class data members to specified output stream.
 *
 *************************************************************************
 */

void
SNES_SAMRAIContext::printClassData(
   std::ostream& os) const
{
   os << "\nSNES_SAMRAIContext::printClassData..." << std::endl;
   os << "SNES_SAMRAIContext: this = "
      << (SNES_SAMRAIContext *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_SNES_functions = "
      << (SNESAbstractFunctions *)d_SNES_functions << std::endl;
   os << "d_SNES_solver = " << (SNES)d_SNES_solver << std::endl;
//   os << "d_SLES_solver = " << (SLES)d_SLES_solver << std::endl;
   os << "d_krylov_solver = " << (KSP)d_krylov_solver << std::endl;
   os << "d_jacobian = " << (Mat *)&d_jacobian << std::endl;
   os << "d_preconditioner = " << (PC *)&d_preconditioner << std::endl;

   os << "d_solution_vector = " << (Vec *)&d_solution_vector << std::endl;
   os << "d_residual_vector = " << (Vec *)&d_residual_vector << std::endl;

   os << "d_uses_preconditioner = " << d_uses_preconditioner << std::endl;
   os << "d_uses_explicit_jacobian = " << d_uses_explicit_jacobian << std::endl;

   os << "d_maximum_nonlinear_iterations = "
      << d_maximum_nonlinear_iterations << std::endl;
   os << "d_maximum_function_evals = " << d_maximum_function_evals << std::endl;

   os << "d_absolute_tolerance = " << d_absolute_tolerance << std::endl;
   os << "d_relative_tolerance = " << d_relative_tolerance << std::endl;
   os << "d_step_tolerance = " << d_step_tolerance << std::endl;

   os << "d_forcing_term_strategy = " << d_forcing_term_strategy << std::endl;
   os << "d_forcing_term_flag = " << d_forcing_term_flag << std::endl;

   os << "d_constant_forcing_term = " << d_constant_forcing_term << std::endl;
   os << "d_initial_forcing_term = " << d_initial_forcing_term << std::endl;
   os << "d_EW_choice2_alpha = " << d_EW_choice2_alpha << std::endl;
   os << "d_EW_choice2_gamma = " << d_EW_choice2_gamma << std::endl;
   os << "d_EW_safeguard_exponent = " << d_EW_safeguard_exponent << std::endl;
   os << "d_EW_safeguard_disable_threshold = "
      << d_EW_safeguard_disable_threshold << std::endl;

   os << "d_linear_solver_type = " << d_linear_solver_type << std::endl;
   os << "d_linear_solver_absolute_tolerance = "
      << d_linear_solver_absolute_tolerance << std::endl;
   os << "d_linear_solver_divergence_tolerance = "
      << d_linear_solver_divergence_tolerance << std::endl;
   os << "d_maximum_linear_iterations = "
      << d_maximum_linear_iterations << std::endl;

   os << "d_maximum_gmres_krylov_dimension = "
      << d_maximum_gmres_krylov_dimension << std::endl;
   os << "d_gmres_orthogonalization_algorithm = "
      << d_gmres_orthogonalization_algorithm << std::endl;

   os << "d_differencing_parameter_strategy = "
      << d_differencing_parameter_strategy << std::endl;
   os << "d_function_evaluation_error = "
      << d_function_evaluation_error << std::endl;
}

}
}

#endif
