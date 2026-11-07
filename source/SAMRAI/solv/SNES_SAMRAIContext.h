/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Wrapper for SNES solver for use in a SAMRAI-based application.
 *
 ************************************************************************/

#ifndef included_solv_SNES_SAMRAIContext
#define included_solv_SNES_SAMRAIContext

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
 ************************************************************************
 */

#ifdef HAVE_PETSC
#ifndef included_petsc_snes
#define included_petsc_snes
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif
#ifdef OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif
#include "petscsnes.h"
#endif

#include "SAMRAI/solv/NonlinearSolverStrategy.h"
#include "SAMRAI/solv/SNESAbstractFunctions.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace solv {

/*!
 * Class SNES_SAMRAIContext provides an interface to the SNES
 * nonlinear solver capabilities in PETSc to facilitate their use with
 * SAMRAI.  While PETSc is implemented in an object-based manner, this
 * class makes it easier to use PETSc's routines with SAMRAI data structures.
 * In particular, it hides from the user some of the messy details required
 * to link the C++ class library with PETSc which is written in C and to
 * override PETSc objects, like vectors, with those found in SAMRAI.
 *
 * This class declares five private static member functions to link
 * user-defined routines for nonlinear residual calculation, Jacobian
 * evaluation, preconditioner setup and solve, and Jacobian-vector product
 * operations.  The implementation of these functions is defined by the user
 * in a subclass of the abstract base class SNESAbstractFunctions.
 * The vector objects used within the solver are provided by the
 * PETSc_SAMRAIVectorReal wrapper class.
 *
 * <b> Input Parameters </b>
 *
 * If no parameters are read from input, PETSc defaults are used.  See the
 * PETSc documentation (http://www-unix.mcs.anl.gov/petsc/).
 * for more information on default parameters and SNES functionality.
 *
 * <b> Definitions: </b>
 *    - \b maximum_nonlinear_iterations
 *       maximum number of nonlinear iterations
 *
 *    - \b maximum_function_evals
 *       maximum number of nonlinear function evaluations
 *
 *    - \b uses_preconditioner
 *       whether or not a preconditioner is used in the solution of the
 *       Newton equations.
 *
 *    - \b uses_explicit_jacobian
 *       whether or not the user provides code to explicitly calculate
 *       Jacobian-vector products.
 *
 *    - \b absolute_tolerance
 *       absolute nonlinear convergence tolerance
 *
 *    - \b relative_tolerance
 *       relative nonlinear convergence tolerance
 *
 *    - \b step_tolerance
 *       minimum tolerance on change in solution norm between nonlinear
 *       iterations
 *
 *    - \b forcing_term_strategy
 *       forcing term choice for linear solvers within the inexact Newton
 *       method.
 *
 *    - \b constant_forcing_term
 *       constant relative convergence tolerance in Krylov solver
 *
 *    - \b initial_forcing_term
 *       initial relative convergence tolerance in Krylov solver
 *       (used in Eisenstat-Walker case).
 *
 *    - \b maximum_forcing_term
 *       maximum relative convergence tolerance in Krylov solver
 *       (used in Eisenstat-Walker case).
 *
 *    - \b EW_choice2_alpha
 *       power used in Eisenstat-Walker choice 2 relative convergence
 *       tolerance computation.
 *
 *    - \b EW_choice2_gamma
 *       multiplicative factor used in Eisenstat-Walker choice 2 relative
 *       convergence tolerance computation.
 *
 *    - \b EW_safeguard_exponent
 *       power for safeguard used in Eisenstat-Walker choice 2
 *
 *    - \b EW_safeguard_disable_threshold
 *       threshold for imposing safeguard in Eisenstat-Walker choice 2.
 *
 *    - \b linear_solver_type
 *       value for type of linear solver.
 *
 *    - \b linear_solver_absolute_tolerance
 *       absolute convergence tolerance in linear solver
 *
 *    - \b linear_solver_divergence_tolerance
 *       amount linear solver residual can increase before solver concludes
 *       method is diverging
 *
 *    - \b maximum_linear_iterations
 *       maximum number of linear solver iterations
 *
 *    - \b maximum_gmres_krylov_dimension
 *       maximum dimension of Krylov subspace before restarting.
 *       Valid only if GMRES is used as the linear solver.
 *
 *    - \b gmres_orthogonalization_algorithm
 *       algorithm used to incrementally construct the orthonormal basis of
 *       the Krylov subspace used by GMRES.  Valid only if GMRES is used as
 *       the linear solver
 *
 *    - \b differencing_parameter_strategy
 *       strategy used for computing the differencing parameter when
 *       Jacobian-vector products are approximated via finite differences
 *
 *    - \b function_evaluation_error
 *       square root of the estimated relative error in function evaluation
 *
 * All values read in from a restart database may be overriden by input
 * database values.  If no new input database value is given, the restart
 * database value is used.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>maximum_nonliner_iterations</td>
 *     <td>int</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>maximum_function_evals</td>
 *     <td>int</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>uses_preconditioner</td>
 *     <td>bool</td>
 *     <td>TRUE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>uses_explicit_jacobian</td>
 *     <td>bool</td>
 *     <td>TRUE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>absolute_tolerance</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>relative_tolerance</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>step_tolerance</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>forcing_term_strategy</td>
 *     <td>string</td>
 *     <td>"CONSTANT"</td>
 *     <td>"CONSTANT", "EWCHOICE1", "EWCHOICE2"</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>constant_forcing_term</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>initial_forcing_term</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>maximum_forcing_term</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>EW_choice2_alpha</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>EW_choice2_gamma</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>EW_safeguard_exponent</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>EW_safeguard_disable_threshold</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>linear_solver_type</td>
 *     <td>string</td>
 *     <td>""</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>linear_solver_absolute_tolerance</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>linear_solver_divergence_tolerance</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>maximum_linear_iterations</td>
 *     <td>int</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>maximum_gmres_krylov_dimension</td>
 *     <td>int</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>gmres_orthogonalization_algorithm</td>
 *     <td>string</td>
 *     <td>""</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>differencing_parameter_strategy</td>
 *     <td>string</td>
 *     <td>MATMFFD_WP</td>
 *     <td>MATMFFD_WP, MATFFD_DS</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>function_evaluation_error</td>
 *     <td>double</td>
 *     <td>PETSC_DEFAULT</td>
 *     <td>Refer to PETSc documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * A sample input file entry might look like:
 *
 * @code
 *    absolute_tolerance             = 10.e-10
 *    relative_tolerance             = 10.e-6
 *    step_tolerance                 = 10.e-8
 *    maximum_nonlinear_iterations   = 200
 *    forcing_term_strategy          = "EWCHOICE1"
 * @endcode
 *
 * Note that input values can also be set using accessor functions.
 * Values that are set via this mechanism will be cached both in the
 * solver context as well as in the corresponding PETSc object.  Thus
 * values changed on-the-fly will be written to restart.  These input
 * values can also be changed by directly accessing the corresponding
 * PETSc object and using native PETSc function calls; however such
 * settings/changes will NOT be cached in the solver context, and so
 * will not be written to restart.
 *
 * @see SNESAbstractFunctions
 * @see NonlinearSolverStrategy
 */

class SNES_SAMRAIContext:
   public NonlinearSolverStrategy,
   public tbox::Serializable
{
public:
   /*!
    * Constructor for SNES_SAMRAIContext allocates the SNES
    * object and initializes rudimentary state associated with
    * user-supplied solver components.  Then, it reads solver parameter
    * from input and restart which may override default values.
    *
    * @pre !object_name.empty()
    * @pre my_functions != 0
    */
   SNES_SAMRAIContext(
      const std::string& object_name,
      SNESAbstractFunctions* my_functions,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * Destructor for solve_SNES_SAMRAIContext destroys the SNES
    * and the PETSc solution vector wrapper.
    */
   ~SNES_SAMRAIContext();

   /*!
    * Return the PETSc nonlinear solver object.
    */
   SNES
   getSNESSolver() const
   {
      return d_SNES_solver;
   }

   /*!
    * Return pointer to object providing user-defined functions for SNES.
    */
   SNESAbstractFunctions *
   getSNESFunctions() const
   {
      return d_SNES_functions;
   }

   /*!
    * Return the PETSc linear solver object.
    */
//   SLES getSLESSolver() const;

   /*!
    * Return the PETSc Krylov solver object.
    */
   KSP
   getKrylovSolver() const
   {
      return d_krylov_solver;
   }

   /*!
    * Return the PETSc Mat object for the Jacobian.
    */
   Mat
   getJacobianMatrix() const
   {
      return d_jacobian;
   }

   /*!
    *  Get absolute tolerance for nonlinear solver.
    */
   double
   getAbsoluteTolerance() const
   {
      return d_absolute_tolerance;
   }

   /*!
    *  Set absolute tolerance for nonlinear solver.
    */
   void
   setAbsoluteTolerance(
      double abs_tol)
   {
      d_absolute_tolerance = abs_tol;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get relative tolerance for nonlinear solver.
    */
   double
   getRelativeTolerance() const
   {
      return d_relative_tolerance;
   }

   /*!
    *  Set relative tolerance for nonlinear solver.
    */
   void
   setRelativeTolerance(
      double rel_tol)
   {
      d_relative_tolerance = rel_tol;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get step tolerance for nonlinear solver.
    */
   double
   getStepTolerance() const
   {
      return d_step_tolerance;
   }

   /*!
    *  Set step tolerance for nonlinear solver.
    */
   void
   setStepTolerance(
      double step_tol)
   {
      d_step_tolerance = step_tol;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get maximum iterations for nonlinear solver.
    */
   int
   getMaxNonlinearIterations() const
   {
      return d_maximum_nonlinear_iterations;
   }

   /*!
    *  Set maximum iterations for nonlinear solver.
    */
   void
   setMaxNonlinearIterations(
      int max_nli)
   {
      d_maximum_nonlinear_iterations = max_nli;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get maximum function evaluations by nonlinear solver.
    */
   int
   getMaxFunctionEvaluations() const
   {
      return d_maximum_function_evals;
   }

   /*!
    *  Set maximum function evaluations in nonlinear solver.
    */
   void
   setMaxFunctionEvaluations(
      int max_feval)
   {
      d_maximum_function_evals = max_feval;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get strategy for forcing term.
    */
   std::string
   getForcingTermStrategy() const
   {
      return d_forcing_term_strategy;
   }

   /*!
    *  Set strategy for forcing term.
    *
    * @pre (strategy == "CONSTANT") || (strategy == "EWCHOICE1") ||
    *      (strategy == "EWCHOICE2")
    */
   void
   setForcingTermStrategy(
      std::string& strategy)
   {
      TBOX_ASSERT(strategy == "CONSTANT" ||
         strategy == "EWCHOICE1" ||
         strategy == "EWCHOICE2");
      d_forcing_term_strategy = strategy;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of constant forcing term.
    */
   double
   getConstantForcingTerm() const
   {
      return d_constant_forcing_term;
   }

   /*!
    *  Set value of constant forcing term.
    */
   void
   setConstantForcingTerm(
      double fixed_eta)
   {
      d_constant_forcing_term = fixed_eta;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of initial forcing term.
    */
   double
   getInitialForcingTerm() const
   {
      return d_initial_forcing_term;
   }

   /*!
    *  Set value of initial forcing term.
    */
   void
   setInitialForcingTerm(
      double initial_eta)
   {
      d_initial_forcing_term = initial_eta;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of maximum forcing term.
    */
   double
   getMaximumForcingTerm() const
   {
      return d_maximum_forcing_term;
   }

   /*!
    *  Set value of maximum forcing term.
    */
   void
   setMaximumForcingTerm(
      double max_eta)
   {
      d_maximum_forcing_term = max_eta;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of exponent in Eisenstat-Walker Choice 2 forcing term.
    */
   double
   getEWChoice2Exponent() const
   {
      return d_EW_choice2_alpha;
   }

   /*!
    *  Set value of exponent in Eisenstat-Walker Choice 2 forcing term.
    */
   void
   setEWChoice2Exponent(
      double alpha)
   {
      d_EW_choice2_alpha = alpha;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of exponent in Eisenstat-Walker Choice 2 safeguard.
    */
   double
   getEWChoice2SafeguardExponent() const
   {
      return d_EW_safeguard_exponent;
   }

   /*!
    *  Set value of exponent in Eisenstat-Walker Choice 2 safeguard.
    */
   void
   setEWChoice2SafeguardExponent(
      double beta)
   {
      d_EW_safeguard_exponent = beta;
      d_context_needs_initialization = true;
   }

   /*!
    * Get value of factor used to scale Eisenstat-Walker Choice 2
    * forcing term.
    */
   double
   getEWChoice2ScaleFactor() const
   {
      return d_EW_choice2_gamma;
   }

   /*!
    * Set value of factor used to scale Eisenstat-Walker Choice 2
    * forcing term.
    */
   void
   setEWChoice2ScaleFactor(
      double gamma)
   {
      d_EW_choice2_gamma = gamma;
      d_context_needs_initialization = true;
   }

   /*!
    *  Get value of threshold to disable safeguard in Eisenstat-Walker
    *  forcing terms.
    */
   double
   getEWSafeguardThreshold() const
   {
      return d_EW_safeguard_disable_threshold;
   }

   /*!
    *  Set value of threshold to disable safeguard in Eisenstat-Walker
    *  forcing terms.
    */
   void
   setEWSafeguardThreshold(
      double threshold)
   {
      d_EW_safeguard_disable_threshold = threshold;
      d_context_needs_initialization = true;
   }

   /*!
    * Get type of linear solver.
    */
   std::string
   getLinearSolverType() const
   {
      return d_linear_solver_type;
   }

   /*!
    * Set type of linear solver.
    */
   void
   setLinearSolverType(
      std::string& type)
   {
      d_linear_solver_type = type;
      d_context_needs_initialization = true;
   }

   /*!
    * Get whether a preconditioner is used.
    */
   bool
   getUsesPreconditioner() const
   {
      return d_uses_preconditioner;
   }

   /*!
    * Set whether to use a preconditioner.
    */
   void
   setUsesPreconditioner(
      bool uses_preconditioner)
   {
      d_uses_preconditioner = uses_preconditioner;
      d_context_needs_initialization = true;
   }

   /*!
    * Get absolute tolerance for linear solver.
    */
   double
   getLinearSolverAbsoluteTolerance() const
   {
      return d_linear_solver_absolute_tolerance;
   }

   /*!
    * Set absolute tolerance for linear solver.
    */
   void
   setLinearSolverAbsoluteTolerance(
      double abs_tol)
   {
      d_linear_solver_absolute_tolerance = abs_tol;
      d_context_needs_initialization = true;
   }

   /*!
    * Get divergence tolerance for linear solver.
    */
   double
   getLinearSolverDivergenceTolerance() const
   {
      return d_linear_solver_divergence_tolerance;
   }

   /*!
    * Set divergence tolerance for linear solver.
    */
   void
   setLinearSolverDivergenceTolerance(
      double div_tol)
   {
      d_linear_solver_divergence_tolerance = div_tol;
      d_context_needs_initialization = true;
   }

   /*!
    * Get maximum linear iterations for linear solver.
    */
   int
   getMaximumLinearIterations() const
   {
      return d_maximum_linear_iterations;
   }

   /*!
    * Set maximum linear iterations for linear solver.
    */
   void
   setMaximumLinearIterations(
      int max_li)
   {
      d_maximum_linear_iterations = max_li;
      d_context_needs_initialization = true;
   }

   /*!
    * Get maximum Krylov dimension in GMRES linear solver.
    */
   int
   getMaximumGMRESKrylovDimension() const
   {
      return d_maximum_gmres_krylov_dimension;
   }

   /*!
    * Set maximum Krylov dimension in GMRES linear solver.
    */
   void
   setMaximumGMRESKrylovDimension(
      int d)
   {
      d_maximum_gmres_krylov_dimension = d;
      d_context_needs_initialization = true;
   }

   /*!
    * Get orthogonalization method used in GMRES linear solver.
    */
   std::string
   getGMRESOrthogonalizationMethod() const
   {
      return d_gmres_orthogonalization_algorithm;
   }

   /*!
    * Set orthogonalization method used in GMRES linear solver.
    */
   void
   setGMRESOrthogonalizationMethod(
      std::string& method)
   {
      d_gmres_orthogonalization_algorithm = method;
      d_context_needs_initialization = true;
   }

   /*!
    * Get whether a method for explicit Jacobian-vector products is provided.
    */
   bool
   getUsesExplicitJacobian() const
   {
      return d_uses_explicit_jacobian;
   }

   /*!
    * Set whether a method for explicit Jacobian-vector products is provided.
    */
   void
   setUsesExplicitJacobian(
      bool use_jac)
   {
      d_uses_explicit_jacobian = use_jac;
      d_context_needs_initialization = true;
   }

   /*!
    * Get method for computing differencing parameter.
    */
   std::string
   getDifferencingParameterMethod() const
   {
      return d_differencing_parameter_strategy;
   }

   /*!
    * Set method for computing differencing parameter.
    */
   void
   setDifferencingParameterMethod(
      std::string& method)
   {
      d_differencing_parameter_strategy = method;
      d_context_needs_initialization = true;
   }

   /*!
    * Get estimate of error in function evaluation.
    */
   double
   getFunctionEvaluationError() const
   {
      return d_function_evaluation_error;
   }

   /*!
    * Set estimate of error in function evaluation.
    */
   void
   setFunctionEvaluationError(
      double evaluation_error)
   {
      d_function_evaluation_error = evaluation_error;
      d_context_needs_initialization = true;
   }

   /*!
    * Initialize the state of the SNES solver based on vector argument
    * representing the solution of the nonlinear system.  In general, this
    * routine must be called before the solve() routine is invoked.
    *
    * @pre solution
    */
   void
   initialize(
      const std::shared_ptr<SAMRAIVectorReal<double> >& solution);

   /*!
    *  Reset the state of the nonlinear solver after regridding.
    */
   void
   resetSolver(
      const int coarsest_level,
      const int finest_level);

   /*!
    * Solve the nonlinear problem.  In general, the initialize() routine
    * must be called before this solve function to set up the solver.
    * Returns 1 if successful, 0 otherwise.
    */
   int
   solve();

   /*!
    * Obtain number of nonlinear iterations.
    */
   int
   getNumberOfNonlinearIterations() const
   {
      return d_nonlinear_iterations;
   }

   /*!
    * Obtain total number of linear iterations accumulated over all
    * nonlinear iterations.
    */
   int
   getTotalNumberOfLinearIterations() const
   {
      int linear_itns;
      int ierr = SNESGetLinearSolveIterations(d_SNES_solver, &linear_itns);
      PETSC_SAMRAI_ERROR(ierr);
      return linear_itns;
   }

   /*!
    * Report reason for termination.
    */
   void
   reportCompletionCode(
      std::ostream& os = tbox::plog) const;

   /*!
    * Write solver parameters to restart database matching object name.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * Print out all members of integrator instance to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int SOLV_SNES_SAMRAI_CONTEXT_VERSION;

   /*
    * Static member functions for linkage with PETSc routines.
    * See header file for SNESAbstractFunctions for more information.
    */
   static int
   SNESFuncEval(
      SNES snes,                            // SNES context
      Vec x,                                // input vector
      Vec f,                                // residual vector
      void* ctx)                            // user-defined context
   {
      NULL_USE(snes);
      ((SNES_SAMRAIContext *)ctx)->getSNESFunctions()->
      evaluateNonlinearFunction(x, f);
      return 0;
   }

   static int
   SNESJacobianSet(
      SNES snes,                                      // SNES context
      Vec x,                                          // input vector
      Mat A,                                          // Jacobian matrix
      Mat B,                                          // precond matrix
      void* ctx);                                     // user-defined context

   static int
   SNESJacobianTimesVector(
      Mat M,                                     //  Jacobian matrix
      Vec xin,                                   //  input vector
      Vec xout)                                  //  output vector
   {
      void* ctx;
      int ierr = MatShellGetContext(M, &ctx);
      PETSC_SAMRAI_ERROR(ierr);
      return ((SNES_SAMRAIContext *)ctx)->
             getSNESFunctions()->jacobianTimesVector(xin, xout);
   }

   static int
   SNESsetupPreconditioner(
      PC pc)
   {
      int ierr = 0;
      Vec current_solution;
      void* ctx;
      PCShellGetContext(pc, &ctx);
      ierr = SNESGetSolution(((SNES_SAMRAIContext *)ctx)->getSNESSolver(),
            &current_solution);
      PETSC_SAMRAI_ERROR(ierr);
      return ((SNES_SAMRAIContext *)ctx)->
             getSNESFunctions()->setupPreconditioner(current_solution);
   }

   static int
   SNESapplyPreconditioner(
      PC pc,
      Vec xin,                                    // input vector
      Vec xout)                                   // output vector
   {
      void* ctx;
      PCShellGetContext(pc, &ctx);
      return ((SNES_SAMRAIContext *)ctx)->
             getSNESFunctions()->applyPreconditioner(xin, xout);
   }

   /*
    * Create and cache needed Petsc objects.
    */
   void
   createPetscObjects();

   /*
    * Initialize cached Petsc objects.
    */
   void
   initializePetscObjects();

   /*
    * Destroy cached Petsc objects.
    */
   void
   destroyPetscObjects();

   /*
    * Read input values from given database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /*!
    * Read solver parameters from restart database matching object name.
    */
   void
   getFromRestart();

   /*!
    * Internal state parameters:
    */
   std::string d_object_name;
   bool d_context_needs_initialization;

   /*
    * PETSc solver and preconditioner objects:
    *
    * d_SNES_solver ..... PETSc nonlinear solver object.
    *
    * d_linear_solver ... PETSc linear solver object, cached here so that
    *                     users may manipulate it through the
    *                     interface of this class.
    *
    * d_krylov_solver ... PETSc Krylov solver context.
    *
    * d_jacobian ........ PETSc matrix object, cached here so that users
    *                     may specify Jacobian operations through the
    *                     interface of this class without having to know
    *                     about PETSc matrix shells or matrix-free matrices.
    *
    * d_preconditioner... PETSc preconditioner object, cached here so that
    *                     users may specify operations through the
    *                     interface of this class without having to know
    *                     about PETSc preconditioner shells.
    */

   SNES d_SNES_solver;
   KSP d_krylov_solver;
   Mat d_jacobian;
   PC d_preconditioner;

   /*
    * Solution and residual vectors for nonlinear system.
    */

   Vec d_solution_vector;
   Vec d_residual_vector;

   /*
    * Pointer to object which provides user-supplied functions to SNES.
    */
   SNESAbstractFunctions* d_SNES_functions;

   /*
    * Boolean flags used during SNES initialization to provide correct
    * static function linkage with PETSc.
    */
   bool d_uses_preconditioner;
   bool d_uses_explicit_jacobian;

   /*
    * SNES state data maintained here for input/restart capabilities.
    */

   // Nonlinear solver parameters:

   int d_maximum_nonlinear_iterations;
   int d_maximum_function_evals;

   double d_absolute_tolerance;
   double d_relative_tolerance;
   double d_step_tolerance;

   std::string d_forcing_term_strategy;  // string is for input
   int d_forcing_term_flag;         // int is for passing choice to PETSc

   double d_constant_forcing_term;
   double d_initial_forcing_term;
   double d_maximum_forcing_term;
   double d_EW_choice2_alpha;
   double d_EW_choice2_gamma;
   double d_EW_safeguard_exponent;
   double d_EW_safeguard_disable_threshold;

   SNESConvergedReason d_SNES_completion_code;

   // Linear solver parameters:

   std::string d_linear_solver_type;
   double d_linear_solver_absolute_tolerance;
   double d_linear_solver_divergence_tolerance;
   int d_maximum_linear_iterations;

   int d_maximum_gmres_krylov_dimension;
   std::string d_gmres_orthogonalization_algorithm;

   // "Matrix-free" parameters:

   std::string d_differencing_parameter_strategy;
   double d_function_evaluation_error;

   // Output parameters:

   int d_nonlinear_iterations;
};

}
}

#endif
#endif
