/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Wrapper class for KINSOL solver function calls and data
 *
 ************************************************************************/

#ifndef included_solv_KINSOLSolver
#define included_solv_KINSOLSolver

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT KINSOL
 ************************************************************************
 */
#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/KINSOLAbstractFunctions.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"

extern "C" {
#include "kinsol/kinsol.h"
#include "sunlinsol/sunlinsol_spgmr.h"
}

#include "kinsol/kinsol_spils.h"

#include <string>

#ifndef LACKS_SSTREAM
#define KINSOL_SAMRAI_ERROR(ierr)                                   \
   do {                                                             \
      if (ierr != KIN_SUCCESS) {                                    \
         std::ostringstream tboxos;                                 \
         tbox::Utilities::abort(                                    \
            tboxos.str().c_str(), __FILE__, __LINE__);              \
      }                                                             \
   } while (0)
#else
#define KINSOL_SAMRAI_ERROR(ierr)                                   \
   do {                                                             \
      if (ierr != KIN_SUCCESS) {                                    \
         std::ostrstream tboxos;                                    \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);  \
      }                                                             \
   } while (0)
#endif

#define SABSVEC_CAST(v) \
   (static_cast<SundialsAbstractVector *>(v \
                                          -> \
                                          content))

namespace SAMRAI {
namespace solv {

/**
 * Class KINSOLSolver serves as a C++ wrapper for the KINSOL nonlinear
 * algebraic equation solver package and its data structures.  It is intended
 * to be sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares four private static member functions to link
 * user-defined routines for nonlinear residual calculation, preconditioner
 * setup and solve, and Jacobian-vector product.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base class
 * KINSOLAbstractFunctions.  The vector objects used within the solver
 * are given in a subclass of the abstract class SundialsAbstractVector.
 * The SundialsAbstractVector class defines the vector kernel operations
 * required by the KINSOL package so that they may be easily supplied
 * by a user who opts not to use the vector kernel supplied by the KINSOL
 * package.
 *
 * Note that this class provides no input or restart capabilities and
 * relies on KINSOL for output reporting.  When using KINSOL in an
 * application using SAMRAI, it is straightforward to include this
 * functionality in the entity using this solver class.
 *
 * KINSOL was developed in the Center for Applied Scientific Computing (CASC)
 * at Lawrence Livermore National Laboratory (LLNL).  For more information
 * about KINSOL and a complete description of the operations and data
 * structures used by this class, see A.G. Taylor and A.C. Hindmarsh,
 * "User documentation for KINSOL, a nonlinear solver for sequential and
 * parallel computers", UCRL-ID-131185, Lawrence Livermore National
 * Laboratory, 1998.
 *
 * @see KINSOLAbstractFunctions
 * @see SundialsAbstractVector
 */

class KINSOLSolver
{
public:
   /**
    * Constructor for KINSOLSolver sets default KINSOL parameters
    * and initializes the solver package with user-supplied functions.  Solver
    * parameters may be changed later using member functions described
    * below.  The integer flags indicate whether user-supplied preconditioner
    * and Jacobian-vector product function should be used.  Zero indicates
    * no user function; otherwise, user function will be used by the nonlinear
    * solver.
    *
    * Important note:  The solution vector is not passed into the constructor.
    * Before the solver can be used, the initialize() function must be called.
    *
    * @pre !object_name.empty()
    * @pre my_functions != 0
    */
   KINSOLSolver(
      const std::string& object_name,
      KINSOLAbstractFunctions* my_functions,
      const int uses_preconditioner,
      const int uses_jac_times_vector);

   /**
    * Virtual destructor for KINSOLSolver.
    */
   virtual ~KINSOLSolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * KINSOL.  This routine must be called before the solver can be used.
    *
    * Optionally set the scaling vectors used by KINSOL to scale
    * either nonlinear solution vector or nonlinear residual vector.
    * The elements of the scaling vectors must be positive.  In either
    * case, the scaling vector should be defined so that the vector
    * formed by taking the element-wise product of the
    * solution/residual vector and scaling vector has all elements
    * roughly the same magnitude when the solution vector IS/IS NOT
    * NEAR a root of the nonlinear function.
    *
    * See KINSOL documentation for more information.
    *
    * @pre solution != 0
    */
   void
   initialize(
      SundialsAbstractVector* solution,
      SundialsAbstractVector* uscale = 0,
      SundialsAbstractVector* fscale = 0);

   /**
    * Solve nonlinear problem and return integer termination code defined
    * by KINSOL.  The default return value is KINSOL_SUCCESS (= 1)
    * indicating success.  Return values which indicate non-recoverable
    * nonlinear solver behavior are KINSOL_NO_MEM (= -1),
    * KINSOL_INPUT_ERROR (= -2), and KINSOL_LSOLV_NO_MEM (= -3).
    * Return values PRECONDSET_FAILURE (= 9), and PRECONDSOLVE_FAILURE (= 10)
    * generally indicate non-recoverable behavior in the preconditioner.
    * See kinsol.h header file for more information about return values.
    *
    * If KINSOL requires re-initialization, it is automatically done before
    * the solve.  This may be required if any of the KINSOL data parameters
    * have changed since the last call to the solver.
    */
   int
   solve();

   /**
    * Accessory function for setting KINSOL output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    *
    * KINSOL printing options are:
    *
    *
    *
    * - \b 0 {no statistics printed}
    * - \b 1 {output iteration count, residual norm, number function calls}
    * - \b 2 {same as 1, but with statistics on globalization process}
    * - \b 3 {same as 2, but with more Krylov iteration statistics}
    *
    *
    *
    * The default is no output (i.e., 0).  If the file name string is empty
    * the default file name "kinsol.log" is used.
    *
    * See KINSOL documentation for more information.
    *
    * @pre flag >= 0 && flag <= 3
    */
   void
   setLogFileData(
      const std::string& log_fname,
      const int flag);

   /**
    * Accessory functions for passing user-defined function information
    * to KINSOL.
    *
    * my_functions is a pointer to the abstract function subclass object
    * that defines the residual calculation and preconditioner functions.
    *
    * uses_preconditioner turns user preconditioner on or off.
    *
    * uses_jac_times_vector turns user Jacobian-vector product on or off.
    *
    * Flags use "TRUE"/"FALSE" values defined in KINSOL.  See KINSOL
    * documentation for more information.
    *
    * @pre my_functions != 0
    */
   void
   setKINSOLFunctions(
      KINSOLAbstractFunctions* my_functions,
      const int uses_preconditioner,
      const int uses_jac_times_vector)
   {
      TBOX_ASSERT(my_functions != 0);
      d_KINSOL_functions = my_functions;
      d_uses_preconditioner = uses_preconditioner;
      d_uses_jac_times_vector = uses_jac_times_vector;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setPreconditioner(
      const int uses_preconditioner)
   {
      d_uses_preconditioner = uses_preconditioner;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setJacobianTimesVector(
      const int uses_jac_times_vector)
   {
      d_uses_jac_times_vector = uses_jac_times_vector;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Return pointer to object that provides user-defined functions for KINSOL.
    */
   KINSOLAbstractFunctions *
   getKINSOLFunctions() const
   {
      return d_KINSOL_functions;
   }

   /**
    * Set constraints on nonlinear solution.  By default the constraint
    * vector is null.
    *
    * The constraints are applied in KINSOL as follows:
    *
    *
    *
    * - \b {if constraints[i] > 0.0, then the constraint is solution[i]>0.0}
    * - \b {if constraints[i] < 0.0, then the constraint is solution[i]<0.0}
    * - \b {if constraints[i] = 0.0, then no constraint on solution[i]}
    *
    *
    *
    *
    * See KINSOL documentation for more information.
    */
   void
   setConstraintVector(
      SundialsAbstractVector* constraints)
   {
      d_constraints = constraints;
   }

   /**
    * Accessory functions for setting nonlinear solver parameters.
    * Parameters and default values are:
    *
    * Residual stopping tolerance is tolerarnce on max_norm(fscale * residual),
    * where product of vectors is another vector each element of which is
    * the product of the corresponding entries in the original vectors.
    * The default is \f$machine_epsilon^(1/3)\f$.
    *
    * Default maximum nonlinear iterations is 200.
    *
    * Default maximum Krylov dimension is 1.
    *
    * Options for global Newton method are: INEXACT_NEWTON = 0, LINESEARCH = 1.
    * The default is INEXACT_NEWTON.
    *
    * Default maximum Newton step is 1000*max(norm(uscale*u_0), norm(uscale)),
    * where u_0 is the initial guess at the solution.
    *
    * Default scaled step tolerarnce between successive nonlinear iterates is
    * \f$machine_epsilon^(2/3)\f$.
    *
    * Default relative error for nonlinear function is set to machine_epsilon.
    *
    * Scalar update constraint value restricts update of solution to
    * del(u)/u < constraint_value.  Here, vector ratio is another vector
    * each element of which is the ratio of the corresponding entries in
    * the original vectors.  The default is no constraint.
    *
    * See KINSOL documentation for more information.
    *
    * @pre (tol == -1.0) || (tol >= 0.0)
    */
   void
   setResidualStoppingTolerance(
      const double tol)
   {
      TBOX_ASSERT(tol == -1.0 || tol >= 0.0);
      d_residual_tol = tol;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setMaxIterations(
      const int maxits)
   {
      TBOX_ASSERT(maxits >= 0);
      d_max_iter = maxits;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setMaxKrylovDimension(
      const int kdim)
   {
      TBOX_ASSERT(kdim >= 0);
      d_krylov_dimension = kdim;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setGlobalStrategy(
      const int global)
   {
      TBOX_ASSERT(global == KIN_NONE || global == KIN_LINESEARCH);
      d_global_strategy = global;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setMaxNewtonStep(
      const double maxstep)
   {
      TBOX_ASSERT(maxstep == -1.0 || maxstep > 0.0);
      d_max_newton_step = maxstep;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setNonlinearStepTolerance(
      const double tol)
   {
      TBOX_ASSERT(tol == -1.0 || tol >= 0.0);
      d_step_tol = tol;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setRelativeFunctionError(
      const double reserr)
   {
      TBOX_ASSERT(reserr == -1.0 || reserr > 0.0);
      d_relative_function_error = reserr;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Accessory functions for setting convergence tests for inner linear
    * solvers within an inexact Newton method.  In general, the linear
    * solver attempts to produce a step p, satisfying:
    * norm(F(u) + J(u)*p) <= (eta + u_round)*norm(F(u)), where the norm
    * is a scaled L2-norm.
    *
    * The convergence test indicates the value for eta; options are:
    *
    *
    *
    * - \b 0 == ETACHOICE1{Choice 1 of Eisenstat and Walker}
    * - \b 1 == ETACHOICE2{Choice 2 of Eisenstat and Walker}
    * - \b 2 == ETACONSTANT{use constant value for eta}.
    *
    *
    *
    * The default option is ETACONSTANT.
    *
    * The default constant value for eta is 0.1.
    *
    * For choice ETACHOICE2, alpha = 2.0 and gamma = 0.9 are defaults.
    *
    * See KINSOL documentation for more information.
    *
    * @pre (conv == KIN_ETACONSTANT) || (conv == KIN_ETACHOICE1) ||
    *      (conv == KIN_ETACHOICE2)
    */
   void
   setLinearSolverConvergenceTest(
      const int conv)
   {
      TBOX_ASSERT(conv == KIN_ETACONSTANT || conv == KIN_ETACHOICE1 ||
         conv == KIN_ETACHOICE2);
      d_eta_choice = conv;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setLinearSolverConstantTolerance(
      const double tol)
   {
      TBOX_ASSERT(tol >= 0.0);
      //
      d_eta_constant = tol;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setEisenstatWalkerParameters(
      const double alpha,
      const double gamma)
   {
      TBOX_ASSERT(alpha >= 0.0);
      TBOX_ASSERT(gamma >= 0.0);
      // sgs
      d_eta_alpha = alpha;
      d_eta_gamma = gamma;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setMaxStepsWithNoPrecondSetup(
      const int maxsolv)
   {
      TBOX_ASSERT(maxsolv > 0);
      d_max_solves_no_set = maxsolv;
      d_KINSOL_needs_initialization = true;
   }

   ///
   void
   setMaxLinearSolveRestarts(
      const int restarts)
   {
      TBOX_ASSERT(restarts >= 0);
      d_max_restarts = restarts;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * The number of nonlinear iterations between checks by the
    * nonlinear residual monitoring algorithm (specifies lenght of
    * subinterval) NOTE: should be a multiple of
    * MaxStepsWithNoPrecondSetup
    */
   void
   setMaxSubSetupCalls(
      const int maxsub)
   {
      TBOX_ASSERT(maxsub >= 0);
      d_maxsub = maxsub;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Set values of omega_min and omega_max scalars used by nonlinear
    * residual monitoring algorithm.
    *
    *  Defaults is [0.00001 and 0.9]
    *
    * @pre omega_min >= 0
    * @pre omega_max >= 0
    * @pre omega_max >= omega_min
    */
   void
   setResidualMonitoringParams(
      const double omega_min,
      const double omega_max)
   {
      TBOX_ASSERT(omega_min >= 0);
      TBOX_ASSERT(omega_max >= 0);
      TBOX_ASSERT(omega_max >= omega_min);
      d_omega_min = omega_min;
      d_omega_max = omega_max;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Set constant value used by residual monitoring algorithm. If
    * omega=0, then it is estimated using omega_min and
    * omega_max.
    *
    * Default is 0.0.
    *
    * @pre omega >= 0
    */
   void
   setResidualMonitoringConstant(
      const double omega)
   {
      TBOX_ASSERT(omega >= 0);
      d_omega = omega;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Set flag controlling whether or not the value * of eps is
    * bounded below by 0.01*fnormtol.
    *
    *       FALSE constrain value of eps by setting to the following:
    *                eps = MAX{0.01*fnormtol, eps}
    *
    *       TRUE do notconstrain value of eps
    *
    * Default is FALSE
    */
   void
   setNoMinEps(
      const bool flag)
   {
      if (flag) {
         d_no_min_eps = 1;
      } else {
         d_no_min_eps = 0;
      }
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Set maximum number of beta condition failures in the line search algorithm.
    *
    * Default is [MXNBCF_DEFAULT] (defined in kinsol_impl.h)
    *
    * @pre max_beta_fails >= 0
    */
   void
   setMaxBetaFails(
      const int max_beta_fails)
   {
      TBOX_ASSERT(max_beta_fails >= 0);
      d_max_beta_fails = max_beta_fails;
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Flag controlling whether or not the KINSol routine makes an
    * initial call to the linearl solver setup routine.
    * Default is false.
    */
   void
   setNoInitialSetup(
      const bool flag)
   {
      if (flag) {
         d_no_initial_setup = 1;
      } else {
         d_no_initial_setup = 0;
      }
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Flag controlling whether or not the nonlinear residual
    * monitoring schemes is used to control Jacobian updating Default
    * is FALSE.
    */
   void
   setNoResidualMonitoring(
      const bool flag)
   {
      if (flag) {
         d_no_residual_monitoring = 1;
      } else {
         d_no_residual_monitoring = 0;
      }
      d_KINSOL_needs_initialization = true;
   }

   /**
    * Accessory functions to retrieve information fom KINSOL.
    *
    * See KINSOL documentation for more information.
    */
   int
   getTotalNumberOfNonlinearIterations() const
   {
      long int num;
      int ierr = KINGetNumNonlinSolvIters(d_kin_mem, &num);
      KINSOL_SAMRAI_ERROR(ierr);
      return static_cast<int>(num);
   }

   ///
   int
   getTotalNumberOfFunctionCalls() const
   {
      long int num;
      int ierr = KINGetNumFuncEvals(d_kin_mem, &num);
      KINSOL_SAMRAI_ERROR(ierr);
      return static_cast<int>(num);
   }

   ///
   int
   getTotalNumberOfBetaConditionFailures() const
   {
      long int num;
      int ierr = KINGetNumBetaCondFails(d_kin_mem, &num);
      KINSOL_SAMRAI_ERROR(ierr);
      return static_cast<int>(num);
   }

   ///
   int
   getTotalNumberOfBacktracks() const
   {
      long int num;
      int ierr = KINGetNumBacktrackOps(d_kin_mem, &num);
      KINSOL_SAMRAI_ERROR(ierr);
      return static_cast<int>(num);
   }

   ///
   double
   getScaledResidualNorm() const
   {
      realtype norm;
      int ierr = KINGetFuncNorm(d_kin_mem, &norm);
      KINSOL_SAMRAI_ERROR(ierr);
      return norm;
   }

   ///
   double
   getNewtonStepLength() const
   {
      realtype step_length;
      int ierr = KINGetStepLength(d_kin_mem, &step_length);
      KINSOL_SAMRAI_ERROR(ierr);
      return step_length;
   }

   /**
    * Print out all data members for this object.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*
    * Open KINSOL log file, allocate main memory for KINSOL and initialize
    * KINSOL memory record.  KINSOL is initialized based on current state
    * of solver parameter data members.  If any solver parameters have
    * changed since last initialization, this function will be automatically
    * invoked at next call to solver.
    */
   void
   initializeKINSOL();

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
    * Free internally allocated vectors.
    */
   void
   freeInternalVectors();

   /*
    * Static member functions for linkage with KINSOL routines.
    * See header file for KINSOLAbstractFunctions for more information.
    */
   static int
   KINSOLFuncEval(
      N_Vector soln,
      N_Vector fval,
      void* my_solver)
   {
      int success = 0;
      // SGS why no error condition?
      ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
      evaluateNonlinearFunction(SABSVEC_CAST(soln), SABSVEC_CAST(fval));
      return success;
   }

   static int
   KINSOLPrecondSet(
      N_Vector uu,
      N_Vector uscale,
      N_Vector fval,
      N_Vector fscale,
      void* my_solver)
   {
      ((KINSOLSolver *)my_solver)->initializeKINSOL();
      int num_feval = 0;
      int success = ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
         precondSetup(SABSVEC_CAST(uu),
            SABSVEC_CAST(uscale),
            SABSVEC_CAST(fval),
            SABSVEC_CAST(fscale),
            num_feval);
      return success;
   }

   static int
   KINSOLPrecondSolve(
      N_Vector uu,
      N_Vector uscale,
      N_Vector fval,
      N_Vector fscale,
      N_Vector vv,
      void* my_solver)

   {
      int num_feval = 0;
      int success = ((KINSOLSolver *)my_solver)->getKINSOLFunctions()->
         precondSolve(SABSVEC_CAST(uu),
            SABSVEC_CAST(uscale),
            SABSVEC_CAST(fval),
            SABSVEC_CAST(fscale),
            SABSVEC_CAST(vv),
            num_feval);
      return success;
   }

   static int
   KINSOLJacobianTimesVector(
      N_Vector v,
      N_Vector Jv,
      N_Vector uu,
      int* new_uu,
      void* my_solver)
   {
      bool soln_changed = true;
      if (*new_uu == 0) {
         soln_changed = false;
      }
      int success = ((KINSOLSolver *)my_solver)->
         getKINSOLFunctions()->jacobianTimesVector(SABSVEC_CAST(v),
            SABSVEC_CAST(Jv),
            soln_changed,
            SABSVEC_CAST(uu));
      return success;
   }

   std::string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the KINSOLSolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * KINSOL may need to be re-initialized (e.g., if Krylov dimension
    * changes, KINSOL must change its memeory record).  Then the
    * initializeKINSOL() member function will be invoked when
    * nonlinear solve function is called next.
    */

   /*
    * Nonlinear solution vector.
    */
   SundialsAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to KINSOL.
    */
   KINSOLAbstractFunctions* d_KINSOL_functions;

   /*
    * Boolean flags used during KINSOL initialization to provide correct
    * static function linkage with KINSOL package.
    */
   bool d_uses_preconditioner;
   bool d_uses_jac_times_vector;

   /*
    * KINSOL input and initialization parameters.
    */
   void* d_kin_mem;                             // KINSOL memory structure
   FILE* d_kinsol_log_file;                     // KINSOL message log file
   std::string d_kinsol_log_file_name;          // KINSOL log file name

   /*
    * Nonlinear solution and residual scaling vectors, and integer flags
    * to determine ownership of scaling vectors.
    */
   SundialsAbstractVector* d_soln_scale;
   bool d_my_soln_scale_vector;
   SundialsAbstractVector* d_fval_scale;
   bool d_my_fval_scale_vector;

   /*
    * Constraints on nonlinear solution vector.
    */
   SundialsAbstractVector* d_constraints;

   /*
    * Linear solver for preconditioning
    */
   SUNLinearSolver d_linear_solver;

   /*
    * Integer flag indicating whether KINSOL needs initialization
    * when solver is called.
    */
   int d_KINSOL_needs_initialization;

   /*
    * KINSOL nonlinear and linear solver parameters
    */
   int d_krylov_dimension;       // maximum krylov dimension
   int d_max_restarts;           // max. num. of linear solver restarts allowed
   int d_max_solves_no_set;      // max. num. of steps calling preconditioner
                                 // without resetting preconditioner

   int d_max_iter;               // maximum number of nonlinear iterations
   double d_max_newton_step;     // maximum scaled length of Newton step

   int d_global_strategy;        // globalization method for Newton steps.
   double d_residual_tol;        // stop tol. on scaled nonlinear residual
   double d_step_tol;            // stop tol. on consecutive step difference

   int d_maxsub;                 // number of nonlinear iterations
                                 // between checks by the nonlinear
                                 // residual monitoring alg

   int d_no_initial_setup;       // intitial setup
   int d_no_residual_monitoring; // residual monitoring to control Jacobian update

   double d_omega_min;           // residual monitoring alg params
   double d_omega_max;
   double d_omega;

   int d_no_min_eps;             //  eps is bounded

   int d_max_beta_fails;         // maximum number of beta condition failures

   // flag indicating which method to use to compute the value of the
   // eta coefficient used in the calculation of the linear solver
   // convergence tolerance:
   int d_eta_choice;

   // KINSetEtaConstValue constant value of eta - use with
   // KIN_ETACONSTANT option
   double d_eta_constant;

   // values of eta_gamma (egamma) and eta_alpha (ealpha) coefficients
   // use with KIN_ETACHOICE2
   double d_eta_gamma;
   double d_eta_alpha;

   // real scalar equal to realative error in computing F(u)
   double d_relative_function_error;

   // level of verbosity of output
   int d_print_level;

};

}
}

#endif
#endif
