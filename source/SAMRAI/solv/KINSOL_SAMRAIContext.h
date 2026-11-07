/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   KINSOL solver for use within a SAMRAI-based application.
 *
 ************************************************************************/

#ifndef included_solv_KINSOL_SAMRAIContext
#define included_solv_KINSOL_SAMRAIContext

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT KINSOL
 ************************************************************************
 */
#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/NonlinearSolverStrategy.h"
#include "SAMRAI/solv/KINSOLSolver.h"
#include "SAMRAI/solv/KINSOLAbstractFunctions.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"

namespace SAMRAI {
namespace solv {

/*!
 * @brief Wraps the KINSOLSolver C++ wrapper class so that
 * KINSOL may be used in applications that require a nonlinear solver.
 *
 * Class KINSOL_SAMRAIContext wraps the KINSOLSolver
 * C++ wrapper class so that KINSOL may be used in applications that
 * require a nonlinear solver.  The class KINSOLSolver does
 * not depend on SAMRAI.  Making it derived from the SAMRAI nonlinear solver
 * interface would require the KINSOL C++ wrapper to depend on SAMRAI.
 *
 * Important note:  This class can only create a KINSOL C++ wrapper instance,
 * initialize it in a rudimentary way, and invoke the solution process.
 * All other interaction with the nonlinear solver (i.e., setting parameters,
 * retrieving solver statistics, etc.) must be done directly with the
 * KINSOL wrapper accessible via the routine getKINSOLSolver().
 * Alternatively, solver parameters may be set up at initialization time
 * using the SAMRAI input database.
 *
 * <b> Input Parameters </b>
 *
 * If no parameters are read from input, KINSOL defaults are used.  See KINSOL
 * documentation for default information.
 *
 * <b> Definitions: </b>
 *  - @b residual_stop_tolerance
 *     stopping tolerance on norm of scaled residual
 *
 *  - @b max_nonlinear_iterations
 *     maximum number of nonlinear iterations
 *
 *  - @b max_krylov_dimension
 *     maximum dimension of Krylov space
 *
 *  - @b global_newton_strategy
 *     globalization strategy
 *
 *  - @b max_newton_step
 *     maximum allowable Newton step
 *
 *  - @b nonlinear_step_tolerance
 *     stopping tolerance on maximum entry in scaled Newton step
 *
 *  - @b relative_function_error
 *     relative error in function evaluation
 *
 *  - @b linear_convergence_test
 *     linear solver convergence tolerance
 *
 *  - @b max_subsetup_calls
 *     number of nonlinear iterations between checks by the nonlinear residual
 *     monitoring algorithm (specifies lenght of subinterval)
 *     NOTE: should be a multiple of max_solves_no_precond_setup
 *
 *  - @b residual_monitoring_params
 *     values of omega_min and omega_max scalars used by nonlinear residual
 *     monitoring algorithm
 *
 *  - @b residual_monitoring_constant
 *     constant value used by residual monitoring algorithm. If omega=0, then
 *     it is estimated using omega_min and omega_max.
 *
 *  - @b no_min_eps
 *     control whether or not the value * of eps is bounded below by
 *     0.01*fnormtol. FALSE = "constrain value of eps by setting to the
 *     following: eps = MAX{0.01*fnormtol, eps}" TRUE = "do
 *     notconstrain value of eps".
 *
 *  - @b eisenstat_walker_params
 *     Eisenstat-Walker choice 2; i.e., the values are given as ETAALPHA,
 *     followed by ETAGAMMA.  Note: the values only apply when linear
 *     convergence test is set to ETACHOICE2.
 *
 *  - @b linear_solver_constant_tolerance
 *     constant linear solver relative tolerance
 *     Note: value only applies when convergence test is set to ETACONSTANT.
 *
 *  - @b max_solves_no_precond_setup
 *     number of nonlinear steps separating successive calls to preconditioner
 *     setup routine
 *
 *  - @b max_linear_solve_restarts
 *     maximum number of linear solver restarts allowed
 *
 *  - @b KINSOL_log_filename
 *     name of KINSOL log file
 *
 *  - @b KINSOL_print_flag
 *     KINSOL log file print options
 *
 *  - @b uses_preconditioner
 *     indicates whether a preconditioner is supplied
 *
 *  - @b uses_jac_times_vector
 *     indicates whether an analytic Jacobian-vector product is supplied
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
 *     <td>residual_stop_tolerance</td>
 *     <td>double</td>
 *     <td>-1.0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_nonlinear_iterations</td>
 *     <td>int</td>
 *     <td>200</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_krylov_dimension</td>
 *     <td>int</td>
 *     <td>10</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>global_newton_strategy</td>
 *     <td>int</td>
 *     <td>0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_newton_step</td>
 *     <td>double</td>
 *     <td>-1.0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>nonlinear_step_tolerance</td>
 *     <td>double</td>
 *     <td>-1.0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>relative_function_error</td>
 *     <td>double</td>
 *     <td>-1.0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>linear_convergence_test</td>
 *     <td>int</td>
 *     <td>3</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_subsetup_calls</td>
 *     <td>int</td>
 *     <td>5</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>residual_monitoring_params</td>
 *     <td>double[2]</td>
 *     <td>[0.00001, 0.9]</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>residual_monitoring_constant</td>
 *     <td>double</td>
 *     <td>0.0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>no_min_eps</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>eisenstat_walker_params</td>
 *     <td>double[2]</td>
 *     <td>[2.0, 0.9]</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>linear_solver_constant_tolerance</td>
 *     <td>double</td>
 *     <td>0.1</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_solves_no_precond_setup</td>
 *     <td>int</td>
 *     <td>10</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_linear_solve_restarts</td>
 *     <td>int</td>
 *     <td>0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>KINSOL_log_filename</td>
 *     <td>string</td>
 *     <td>""</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>KINSOL_print_flag</td>
 *     <td>int</td>
 *     <td>0</td>
 *     <td>Refer to KINSOL documentation</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>uses_preconditioner</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>uses_jac_times_vector</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * A sample input file entry might look like:
 *
 * @code
 *   residual_stop_tolerance  =  10.e-6
 *   max_nonlinear_iterations =  200
 *   max_newton_step          =  0.1
 *   KINSOL_log_filename      =  "mylogfile"
 *   KINSOL_print_flag        =  3   // print all output KINSOL has to offer
 * @endcode
 *
 * @see NonlinearSolverStrategy
 */
class KINSOL_SAMRAIContext:
   public NonlinearSolverStrategy,
   public tbox::Serializable
{
public:
   /**
    * Constructor for algs::KINSOL_SAMRAIContext allocates the KINSOL
    * C++ wrapper object and initializes rudimentary state associated
    * with user-supplied solver components.  Then, it reads solver parameter
    * from input and restart which may override default values.
    *
    * @pre !object_name.empty()
    * @pre my_functions != 0
    */
   KINSOL_SAMRAIContext(
      const std::string& object_name,
      KINSOLAbstractFunctions* my_functions,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /**
    * Destructor for algs::KINSOL_SAMRAIContext destroys the KINSOL
    * C++ wrapper object and the KINSOL solution vector wrapper.
    */
   ~KINSOL_SAMRAIContext();

   /**
    * Initialize the state of KINSOL based on vector argument representing
    * the solution of the nonlinear system.  In general, this routine must
    * be called before the solve() routine is invoked.
    *
    * @pre solution
    */
   void
   initialize(
      const std::shared_ptr<SAMRAIVectorReal<double> >& solution);

   /**
    * Solve the nonlinear problem and return and integer value defined by
    * KINSOL.  A return value of 1 indicates success (i.e., KINSOL_SUCCESS).
    * Consult the KINSOL documentation, KINSOL header file kinsol.h, or the
    * header file for the class KINSOLSolver for more information
    * about KINSOL return codes.  In general, the initialize() routine must
    * be called before this solve function to set up the solver.
    */
   int
   solve();

   /**
    * Return pointer to KINSOL solver C++ wrapper object.
    */
   KINSOLSolver *
   getKINSOLSolver()
   {
      return d_KINSOL_solver;
   }

   /**
    * Read input parameters from given database.
    *
    * @param[in] input_db
    *
    * @param[in] is_from_restart
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /**
    * Retrieve solver parameters from restart database matching object name.
    *
    * When assertion checking is active, an unrecoverable assertion
    * will result if a restart database matching the object name does not
    * exist, or if the class version number does not match that in restart.
    */
   void
   getFromRestart();

   /**
    * Retrieve solver parameters from restart database matching object name.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Print out all members of integrator instance to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
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
   static const int SOLV_KINSOL_SAMRAI_CONTEXT_VERSION;

   std::string d_object_name;

   /*
    * KINSOL nonlinear solver object and solution vector for nonlinear system.
    */

   KINSOLSolver* d_KINSOL_solver;
   SundialsAbstractVector* d_solution_vector;

   /*
    * KINSOL state data maintained here for input/restart capabilities.
    */

   double d_residual_stop_tolerance;
   int d_max_nonlinear_iterations;
   int d_max_krylov_dimension;
   int d_global_newton_strategy;
   double d_max_newton_step;
   double d_nonlinear_step_tolerance;
   double d_relative_function_error;
   int d_linear_convergence_test;
   int d_max_subsetup_calls;
   double d_residual_monitoring_params[2];
   double d_residual_monitoring_constant;
   double d_eisenstat_walker_params[2];
   double d_linear_solver_constant_tolerance;
   int d_max_solves_no_precond_setup;
   int d_max_linear_solve_restarts;
   std::string d_KINSOL_log_filename;
   int d_KINSOL_print_flag;
   bool d_no_min_eps;
   bool d_uses_preconditioner;
   bool d_uses_jac_times_vector;
};

}
}

#endif
#endif
