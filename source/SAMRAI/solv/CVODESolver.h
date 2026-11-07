/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Wrapper class for CVODE solver function calls and data
 *
 ************************************************************************/

#ifndef included_solv_CVODESolver
#define included_solv_CVODESolver

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/CVODEAbstractFunctions.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"

// CVODE includes
extern "C" {
#include "cvode/cvode.h"
}

extern "C" {
#include "cvode/cvode_spils.h"
#include "sunlinsol/sunlinsol_spgmr.h"
}

#include <string>

namespace SAMRAI {
namespace solv {

#ifndef LACKS_SSTREAM
#define CVODE_SAMRAI_ERROR(ierr) \
   do {                                           \
      if (ierr != CV_SUCCESS) {                                                                 \
         std::ostringstream tboxos;                                                     \
         tbox::Utilities::abort( \
            tboxos.str().c_str(), __FILE__, __LINE__);      \
      }                                                                         \
   } while (0)
#else
#define CVODE_SAMRAI_ERROR(ierr) \
   do {                                           \
      if (ierr != CV_SUCCESS) {                                                                 \
         std::ostrstream tboxos;                                                        \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);              \
      }                                                                         \
   } while (0)
#endif

#define SABSVEC_CAST(v) \
   (static_cast<SundialsAbstractVector *>(v \
                                          -> \
                                          content))

/*!
 * @brief Class CVODESolver serves as a C++ wrapper for the CVODE
 * ordinary differential equation solver package.
 *
 * It is intended to be
 * sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares one private static member function to link the
 * user-defined routine for right-hand side function evaluation and
 * two private statice member functions to link the user-defined
 * preconditioner setup and solve routines.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base
 * class CVODEAbstractFunctions.  The vector objects used within the
 * solver are given in a subclass of the abstract class
 * SundialsAbstractVector. The SundialsAbstractVector
 * class defines the vector kernel operations required by the CVODE
 * package so that they may be easily supplied by a user who opts not
 * to use the vector kernel supplied by the CVODE package.  (It should be
 * noted that the vector kernel used by CVODE is the same as the one
 * used by the other packages in the Sundials of solvers).
 *
 * Note that this class provides no input or restart capabilities and
 * relies on CVODE for output reporting.
 *
 * CVODESolver Usage:
 *
 *
 *    -  In order to use the CVODESolver, the user must provide a
 *           concrete subclass of CVODEAbstractFunctions abstract
 *           base class which defines the evaluateRHSFunction(),
 *           CVSpgmrPrecondSet(), and CVSpgmrPrecondSolve() methods.
 *
 *    -  Solving a system of ODEs using this CVODE C++ interface
 *           requires four main stages.  First, a CVODESolver
 *           object is created with a user-specified name and
 *           CVODEAbstractFunctions object.  Second, the
 *           user must specify the integration parameters that s/he
 *           wishes to use.  Next, the user must call the CVODESolver
 *           method initialize(solution_vector) with the
 *           SundialsAbstractVector that s/he wants to put the solution
 *           in.  Finally, the solve() method is invoked to solve the
 *           system of ODEs to the specified value of the independent
 *           variable.
 *
 *    -  The following is a list of integration parameters that
 *           must be specified by the user before calling the solve()
 *           method:
 *
 *
 *            - Either relative or absolute tolerance must
 *                  be set - setRelativeTolerance(relative_tolerance),
 *                  setAbsoluteTolerance(absolute_tolerance)
 *
 *            - Initial value of independent variable -
 *                  setInitialValueOfIndependentVariable(init_time)
 *            - Final value of independent variable -
 *                  setFinalValueOfIndependentVariable(final_time
 *                      cvode_needs_initialization)
 *            - Initial condition vector -
 *                  setInitialConditionVector(ic_vector)
 *
 *
 *
 *    -  The following is a list of default values for integration
 *           parameters:
 *
 *
 *
 *           - @b Linear Multistep Method
 *                BDF
 *
 *           - @b Relative Tolerance
 *                0.0
 *
 *           - @b Scalar Absolute Tolerance
 *                0.0
 *
 *           - @b Vector Absolute Tolerance
 *                NULL
 *
 *           - @b Stepping Method
 *                NORMAL
 *
 *           - @b Maximum Order for Multistep Method
 *                12 for ADAMS, 5 for BDF
 *
 *           - @b Maximum Number of Internal Steps
 *                500
 *
 *           - @b Maximum Number of NIL Step Warnings
 *                10
 *
 *           - @b Initial Step Size
 *                determined by CVODE
 *
 *           - @b Maximum Absolute Value of Step Size
 *                infinity
 *
 *           - @b Minimum Absolute Value of Step Size
 *                0.0
 *
 *           - @b CVSpgmr Preconditioning Type
 *                NONE
 *
 *           - @b CVSpgmr Gram Schmidt Algorithm
 *                MODIFIED_GS
 *
 *           - @b CVSpgmr Maximum Krylov Dimension
 *                MIN(num_equations, CVSPGMR_MAXL=5)
 *
 *           - @b CVSpgmr Tolerance Scale Factor
 *                CVSPGMR_DELT = 0.05.
 *
 *
 *
 *
 *
 * CVODE was developed in the Center for Applied Scientific Computing (CASC)
 * at Lawrence Livermore National Laboratory (LLNL).  Many of the comments
 * in this class were taken verbatim from CVODE header files.  For more
 * information about CVODE and a complete description of the operations
 * and data structures used by this class, see S.D. Cohen and A.C. Hindmarsh,
 * "CVODE User Guide", UCRL-MA-118618, Lawrence Livermore National
 * Laboratory, 1994.
 *
 * @see CVODEAbstractFunctions
 * @see SundialsAbstractVector
 */

class CVODESolver
{
public:
   /**
    * Constructor for CVODESolver sets default CVODE parameters
    * and initializes the solver package with user-supplied functions
    * CVODESolver parameters may be changed later using member
    * functions described below.
    *
    * Notes:
    *
    *
    *
    *
    *    -
    *        The solution vector is not passed into the constructor.
    *        Before the solver can be used, the initialize() function must
    *        be called.
    *
    * @pre !object_name.empty()
    * @pre my_functions != 0
    */
   CVODESolver(
      const std::string& object_name,
      CVODEAbstractFunctions* my_functions,
      const bool uses_preconditioner);

   /**
    * Virtual destructor for CVODESolver closes the
    * CVODE log file and frees the memory allocated for the
    * CVODE memory record.
    */
   virtual ~CVODESolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * CVODE.  This routine must be called before the solver can be used.
    *
    * @pre solution != 0
    * @pre d_solution_vector == 0
    */
   void
   initialize(
      SundialsAbstractVector* solution)
   {
      TBOX_ASSERT(solution != 0);
      TBOX_ASSERT(d_solution_vector == 0);
      d_solution_vector = solution;
      d_CVODE_needs_initialization = true;
      initializeCVODE();
   }

   /**
    * Integrate ODE system specified t_f.  The integer return value is
    * a termination code defined by CVODE.  The following is a table
    * of termination codes and a brief description of their meanings.
    *
    * CVODE Termination Codes:
    *
    *
    *
    *
    *    - @b SUCCESS (=0)
    *        CVode succeeded.
    *
    *    - @b CVODE_NO_MEM (=-1)
    *        The cvode_mem argument was NULL.
    *
    *    - @b ILL_INPUT (=-2)
    *        One of the inputs to CVode is illegal. This
    *        includes the situation when a component of the
    *        error weight vectors becomes < 0 during
    *        internal time-stepping. The ILL_INPUT flag
    *        will also be returned if the linear solver
    *        routine CV--- (called by the user after
    *        calling CVodeMalloc) failed to set one of the
    *        linear solver-related fields in cvode_mem or
    *        if the linear solver's init routine failed. In
    *        any case, the user should see the printed
    *        error message for more details.
    *
    *    - @b TOO_MUCH_WORK (=-3)
    *        The solver took maxstep internal steps but
    *        could not reach t_f. The default value for
    *        mxstep is MXSTEP_DEFAULT = 500.
    *
    *    - @b TOO_MUCH_ACC (=-4)
    *        The solver could not satisfy the accuracy
    *        demanded by the user for some internal step.
    *
    *    - @b ERR_FAILURE (=-5)
    *        Error test failures occurred too many times
    *        (= MXNEF = 7) during one internal time step or
    *        occurred with |h| = hmin.
    *
    *    - @b CONV_FAILURE (=-6)
    *        Convergence test failures occurred too many
    *        times (= MXNCF = 10) during one internal time
    *        step or occurred with |h| = hmin.
    *
    *    - @b SETUP_FAILURE (=-7)
    *        The linear solver's setup routine failed in an
    *                 unrecoverable manner.
    *
    *    - @b SOLVE_FAILURE (=-8)
    *        The linear solver's solve routine failed in an
    *                 unrecoverable manner.
    *
    *
    *
    *
    *
    * See cvode.h header file for more information about return values.
    *
    * If CVODE or CVSpgmr requires re-initialization, it is
    * automatically done before the solve.  This may be required if any
    * of the CVODE or CVSpgmr data parameters have changed since the
    * last call to the solver.
    *
    * @pre d_user_t_f > d_t_0
    */
   int
   solve()
   {
      initializeCVODE();

      /*
       * Check to make sure that user specified final value for t
       * is greater than initial value for t.
       */
      TBOX_ASSERT(d_user_t_f > d_t_0);

      /*
       * See cvode.h header file for definition of return types.
       */
      int retval = CVode(d_cvode_mem,
            d_user_t_f,
            d_solution_vector->getNVector(),
            &d_actual_t_f,
            d_stepping_method);
      return retval;
   }

   /**
    * Accessor function for setting CVODE output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    */
   void
   setLogFileData(
      const std::string& log_fname = std::string())
   {
      if (log_fname != d_cvode_log_file_name) {
         if (log_fname.empty()) {
            d_cvode_log_file_name = "cvode.log";
         } else {
            d_cvode_log_file_name = log_fname;
         }
         d_CVODE_needs_initialization = true;
      }
   }

   /**
    * Set CVODESolver to use my_functions as the concrete subclass
    * of the CVODEAbstractFunctions class that defines the
    * right-hand side evaluation and preconditioner functions.  The
    * uses_preconditioner argument indicates whether or not the
    * the user has defined preconditioner routines in their concrete
    * subclass of the CVODEAbstractFunctions class.
    *
    * @pre my_functions != 0
    */
   void
   setCVODEFunctions(
      CVODEAbstractFunctions* my_functions,
      const bool uses_preconditioner)
   {
      TBOX_ASSERT(my_functions != 0);
      d_cvode_functions = my_functions;
      d_uses_preconditioner = uses_preconditioner;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Return pointer to object that provides user-defined functions for
    * CVODE and CVSpgmr.
    */
   CVODEAbstractFunctions *
   getCVODEFunctions() const
   {
      return d_cvode_functions;
   }

   // Methods for setting CVODE parameters.

   /**
    * Set linear multistep method.  The user can specify either
    * ADAMS or BDF (backward differentiation formula) methods
    * The BDF method is recommended  for stiff problems, and
    * the ADAMS method is recommended for nonstiff problems.
    *
    * Note: the enumeration constants ADAMS and BDF are defined in cvode.h.
    *
    * @pre (linear_multistep_method == CV_ADAMS) ||
    *      (linear_multistep_method == CV_BDF)
    */
   void
   setLinearMultistepMethod(
      int linear_multistep_method)
   {
      TBOX_ASSERT((linear_multistep_method == CV_ADAMS) ||
         (linear_multistep_method == CV_BDF));
      d_linear_multistep_method = linear_multistep_method;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the relative tolerance level.
    *
    * Note that pure absolute tolerance can be used by
    * setting the relative tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre relative_tolerance >= 0.0
    */
   void
   setRelativeTolerance(
      double relative_tolerance)
   {
      TBOX_ASSERT(relative_tolerance >= 0.0);
      d_relative_tolerance = relative_tolerance;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the scalar absolute tolerance level.
    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre absolute_tolerance >= 0.0
    */
   void
   setAbsoluteTolerance(
      double absolute_tolerance)
   {
      TBOX_ASSERT(absolute_tolerance >= 0.0);
      d_absolute_tolerance_scalar = absolute_tolerance;
      d_use_scalar_absolute_tolerance = true;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the vector absolute tolerance level.
    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre absolute_tolerance != 0
    * @pre absolute_tolerance->vecMin() >= 0.0
    */
   void
   setAbsoluteTolerance(
      SundialsAbstractVector* absolute_tolerance)
   {
      TBOX_ASSERT(absolute_tolerance != 0);
      TBOX_ASSERT(absolute_tolerance->vecMin() >= 0.0);
      d_absolute_tolerance_vector = absolute_tolerance;
      d_use_scalar_absolute_tolerance = false;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set stepping method to use for integration.  There are
    * stepping methods: NORMAL and ONE_STEP.  The NORMAL
    * method has the solver take internal steps until
    * it has reached or just passed the user specified t_f
    * parameter. The solver then interpolates in order to
    * return an approximate value of y(t_f). The ONE_STEP
    * option tells the solver to just take one internal step
    * and return the solution at the point reached by that
    * step.
    *
    * Note: the enumeration constants NORMAL and ONE_STEP are
    * defined in cvode.h.
    *
    * @pre (stepping_method == CV_NORMAL) || (stepping_method == CV_ONE_STEP)
    */
   void
   setSteppingMethod(
      int stepping_method)
   {
      TBOX_ASSERT((stepping_method == CV_NORMAL) ||
         (stepping_method == CV_ONE_STEP));
      d_stepping_method = stepping_method;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set initial value for independent variable.
    */
   void
   setInitialValueOfIndependentVariable(
      double t_0)
   {
      d_t_0 = t_0;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set final value for independent variable (i.e. the value of
    * independent variable to integrate the system to).  The boolean
    * argument specifies whether CVODE should be re-initialized (i.e.
    * on first step) or if we are taking subsequent steps in a
    * sequence, in which case it is not initialized.
    */
   void
   setFinalValueOfIndependentVariable(
      double t_f,
      bool cvode_needs_initialization)
   {
      d_user_t_f = t_f;
      d_CVODE_needs_initialization = cvode_needs_initialization;
   }

   /**
    * Set initial condition vector.
    *
    * @pre ic_vector != 0
    */
   void
   setInitialConditionVector(
      SundialsAbstractVector* ic_vector)
   {
      TBOX_ASSERT(ic_vector != 0);
      d_ic_vector = ic_vector;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set maximum order for the linear multistep method.
    * By default, this is set to 12 for ADAMS methods and 5 for BDF
    * methods.
    *
    * @pre max_order >= 0
    */
   void
   setMaximumLinearMultistepMethodOrder(
      int max_order)
   {
      TBOX_ASSERT(max_order >= 0);
      d_max_order = max_order;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set maximum number of internal steps to be taken by
    * the solver in its attempt to reach t_f.
    * By default, this is set to 500.
    *
    * @pre max_num_internal_steps >= 0
    */
   void
   setMaximumNumberOfInternalSteps(
      int max_num_internal_steps)
   {
      TBOX_ASSERT(max_num_internal_steps >= 0);
      d_max_num_internal_steps = max_num_internal_steps;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set maximum number of warning messages issued by the solver
    * that (t + h == t) on the next internal step.  By default,
    * this is set to 10.
    *
    * @pre max_num_warnings >= 0
    */
   void
   setMaximumNumberOfNilStepWarnings(
      int max_num_warnings)
   {
      TBOX_ASSERT(max_num_warnings >= 0);
      d_max_num_warnings = max_num_warnings;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set initial step size.
    *
    * @pre init_step_size >= 0.0
    */
   void
   setInitialStepSize(
      double init_step_size)
   {
      TBOX_ASSERT(init_step_size >= 0.0);
      d_init_step_size = init_step_size;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set maximum absolute value of step size allowed.
    * By default, there is no upper bound on the absolute value
    * of step size.
    *
    * @pre max_step_size >= 0.0
    */
   void
   setMaximumAbsoluteStepSize(
      double max_step_size)
   {
      TBOX_ASSERT(max_step_size >= 0.0);
      d_max_step_size = max_step_size;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set minimum absolute value of step size allowed.
    * By default, this is set to 0.0.
    *
    * @pre min_step_size >= 0.0
    */
   void
   setMinimumAbsoluteStepSize(
      double min_step_size)
   {
      TBOX_ASSERT(min_step_size >= 0.0);
      d_min_step_size = min_step_size;
      d_CVODE_needs_initialization = true;
   }

   // Methods for setting CVSpgmr parameters.

   /**
    * Set the preconditioning type to be used by CVSpgmr.
    * This must be one of the four enumeration constants
    * NONE, LEFT, RIGHT, or BOTH defined in iterativ.h.
    * These correspond to no preconditioning, left preconditioning only,
    * right preconditioning only, and both left and right
    * preconditioning, respectively.
    *
    * @pre (precondition_type == PREC_NONE) ||
    *      (precondition_type == PREC_LEFT) ||
    *      (precondition_type == PREC_RIGHT) ||
    *      (precondition_type == PREC_BOTH)
    */
   void
   setPreconditioningType(
      int precondition_type)
   {
      TBOX_ASSERT((precondition_type == PREC_NONE) ||
         (precondition_type == PREC_LEFT) ||
         (precondition_type == PREC_RIGHT) ||
         (precondition_type == PREC_BOTH));
      d_precondition_type = precondition_type;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the Gram-Schmidt orthogonalization type to be used by CVSpgmr.
    * This must be one of the two enumeration constants MODIFIED_GS
    * or CLASSICAL_GS defined in iterativ.h. These correspond to
    * using modified Gram-Schmidt and classical Gram-Schmidt, respectively.
    *
    * @pre (gs_type == CLASSICAL_GS) || (gs_type == MODIFIED_GS)
    */
   void
   setGramSchmidtType(
      int gs_type)
   {
      TBOX_ASSERT((gs_type == CLASSICAL_GS) || (gs_type == MODIFIED_GS));
      d_gram_schmidt_type = gs_type;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the maximum Krylov dimension to be used by CVSpgmr.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value MIN(num_equations, CVSPGMR_MAXL=5).
    *
    * @pre max_krylov_dim >= 0
    */
   void
   setMaxKrylovDimension(
      int max_krylov_dim)
   {
      TBOX_ASSERT(max_krylov_dim >= 0);
      d_max_krylov_dim = max_krylov_dim;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Set the factor by which the tolerance on the nonlinear
    * iteration is multiplied to get a tolerance on the linear iteration.
    * This is an optional input to the CVSPGMR solver. Pass 0 to
    * use the default value CVSPGMR_DELT = 0.05.
    *
    * @pre tol_scale_factor >= 0
    */
   void
   setCVSpgmrToleranceScaleFactor(
      double tol_scale_factor)
   {
      TBOX_ASSERT(tol_scale_factor >= 0);
      d_tol_scale_factor = tol_scale_factor;
      d_CVODE_needs_initialization = true;
   }

   /**
    * Enable use of a user -defined projection function. The uses_projectionfn
    * argument indicates whether or not the user has defined a projection
    * routines in their concrete subclass of the CVODEAbstractFunctions class.
    * Note this function must be called before initializeCVODE().
    */
   void
   setProjectionFunction(bool uses_projectionfn)
   {
      d_uses_projectionfn = uses_projectionfn;
   }

   /**
    * Enable use of a different RHS function in Jacobian -vector products.
    * The uses_jtimesrhsfn argument indicates whether or not the user has
    * defined a projection routine in their concrete subclass of the
    * CVODEAbstractFunctions class. Note this function must be called before
    * initializeCVODE().
    */
   void
   setJTimesRhsFunction(bool uses_jtimesrhsfn)
   {
      d_uses_jtimesrhsfn = uses_jtimesrhsfn;
   }

   /**
    * Get solution vector.
    */
   SundialsAbstractVector *
   getSolutionVector() const
   {
      return d_solution_vector;
   }

   /**
    * Get k-th derivative vector at the specified value of the
    * independent variable, t.  The integer return value is
    * return code the CVODE CVodeDky() function.  The following is a table
    * of termination codes and a brief description of their meanings.
    *
    * CVodeDky Return Codes:
    *
    *
    *
    *
    *    - @b OKAY (=0)
    *        CVodeDky succeeded.
    *
    *    - @b BAD_K (=-1)
    *
    *    - @b BAD_T (=-2)
    *
    *    - @b BAD_DKY (=-3)
    *
    *    - @b DKY_NO_MEM (=-4)
    *
    *
    *
    *
    *
    * Important Notes:
    *
    *
    *
    *
    *    -
    *       t must lie in the interval [t_cur - h, t_cur]
    *       where t_cur is the current internal time reached
    *       and h is the last internal step size successfully
    *       used by the solver.
    *
    *    -
    *       k may take on value 0, 1, . . . q where q is the order
    *       of the current linear multistep method being used.
    *
    *    -
    *       the dky vector must be allocated by the user.
    *
    *    -
    *       it is only leagal to call this method after a
    *       successful return from the solve() method.
    *
    *
    *
    *
    *
    */
   int
   getDkyVector(
      double t,
      int k,
      SundialsAbstractVector* dky) const
   {
      int return_code = CVodeGetDky(d_cvode_mem, t, k, dky->getNVector());
      return return_code;
   }

   /**
    * Get actual value of the independent variable that CVODE integrated
    * to (i.e. the value of t that actually corresponds to the solution
    * vector y).
    */
   double
   getActualFinalValueOfIndependentVariable() const
   {
      return d_actual_t_f;
   }

   /**
    * Print CVODE and CVSpgmr statistics.
    */
   void
   printStatistics(
      std::ostream& os) const
   {
      printCVODEStatistics(os);
      printCVSpgmrStatistics(os);
   }

   /**
    * Print CVODE statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    *
    *
    *
    *
    *    - @b lenrw
    *       size (in double words) of memory used for doubles
    *
    *    - @b leniw
    *       size (in integer words) of memory used for integers
    *
    *    - @b nst
    *       cumulative number of internal steps taken by solver
    *
    *    - @b nfe
    *       number of right-hand side function evaluations
    *
    *    - @b nni
    *       number of NEWTON iterations performed
    *
    *    - @b nsetups
    *       number of calls made to linear solver's setup routine
    *
    *    - @b netf
    *       number of local error test failures
    *
    *    - @b ncfn
    *       number of nonlinear convergence failures
    *
    *    - @b qu
    *       order used during the last internal step
    *
    *    - @b qcur
    *       order to be used on the next internal step
    *
    *    - @b hu
    *       step size for the last internal step
    *
    *    - @b hcur
    *       step size to be attempted on the next internal step
    *
    *    - @b tcur
    *       current internal value of t reached by the solver
    *
    *    - @b tolsf
    *       suggested tolerance scaling factor
    *
    *
    *
    *
    */
   void
   printCVODEStatistics(
      std::ostream& os) const;

   // CVODE optional return values.

   /**
    * Return the cumulative number of internal steps taken by
    * the solver.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfInternalStepsTaken() const
   {
      long int r;
      int ierr = CVodeGetNumSteps(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of calls to the right-hand side function.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfRHSFunctionCalls() const
   {
      long int r;
      int ierr = CVodeGetNumRhsEvals(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of calls made to linear solver setup
    * routines.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfLinearSolverSetupCalls() const
   {
      long int r;
      int ierr = CVodeGetNumLinSolvSetups(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of NEWTON iterations performed.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfNewtonIterations() const
   {
      long int r;
      int ierr = CVodeGetNumNonlinSolvIters(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of nonlinear convergence failures that have
    * occurred.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfNonlinearConvergenceFailures() const
   {
      long int r;
      int ierr = CVodeGetNumNonlinSolvConvFails(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of local error test failures.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getNumberOfLocalErrorTestFailures() const
   {
      long int r;
      int ierr = CVodeGetNumErrTestFails(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the order of the linear multistep method used during
    * the last internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getOrderUsedDuringLastInternalStep() const
   {
      int r;
      int ierr = CVodeGetLastOrder(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the order of the linear multistep method to be used during
    * the next internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getOrderToBeUsedDuringNextInternalStep() const
   {
      int r;
      int ierr = CVodeGetCurrentOrder(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the size (in LLNL_REAL words) of memory used
    * for LLNL_REALS.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getCVODEMemoryUsageForDoubles() const
   {
      long int r1;
      long int r2;
      int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r1);
   }

   /**
    * Return the size (in integer words) of memory used
    * for integers.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int
   getCVODEMemoryUsageForIntegers() const
   {
      long int r1;
      long int r2;
      int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r2);
   }

   /**
    * Return the step size for the last internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double
   getStepSizeForLastInternalStep() const
   {
      realtype r;
      int ierr = CVodeGetLastStep(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return r;
   }

   /**
    * Return the step size to be used in the next internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double
   getStepSizeForNextInternalStep() const
   {
      realtype r;
      int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return r;
   }

   /**
    * Return the current internal value of the independent
    * variable reached by the solver.
    *
    * Note: if the solver was not set to collect statistics,
    * the minimum double value (as defined in float.h) is
    * returned.
    */
   double
   getCurrentInternalValueOfIndependentVariable() const
   {
      realtype r;
      int ierr = CVodeGetCurrentStep(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return r;
   }

   /**
    * Return the suggested tolerance scaling factor.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double
   getCVODESuggestedToleranceScalingFactor() const
   {
      realtype r;
      int ierr = CVodeGetTolScaleFactor(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return r;
   }

   // CVSpgmr optional return values.

   /**
    * Print CVSpgmr statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    *
    *
    *
    *
    *    - @b spgmr_lrw
    *      size (in double words) of memory used for doubles
    *
    *    - @b spgmr_liw
    *       size (in integer words) of memory used for integers
    *
    *    - @b nli
    *       number of linear iterations
    *
    *    - @b ncfl
    *       number of linear convergence failures
    *
    *    - @b npe
    *       number of preconditioner evaluations
    *
    *    - @b nps
    *       number of calls to CVSpgmrPrecondSolve()
    *
    *
    *
    *
    */
   void
   printCVSpgmrStatistics(
      std::ostream& os) const;

   /**
    * Return the number of preconditioner evaluations.
    */
   int
   getNumberOfPreconditionerEvaluations() const
   {
      long int r;
      int ierr = CVSpilsGetNumPrecEvals(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of linear iterations.
    */
   int
   getNumberOfLinearIterations() const
   {
      long int r;
      int ierr = CVSpilsGetNumLinIters(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of CVSpgmrPrecondSolve() calls.
    */
   int
   getNumberOfPrecondSolveCalls() const
   {
      long int r;
      int ierr = CVSpilsGetNumPrecSolves(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of linear convergence failures.
    */
   int
   getNumberOfLinearConvergenceFailures() const
   {
      long int r;
      int ierr = CVSpilsGetNumConvFails(d_cvode_mem, &r);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the size (in double words) of memory used for doubles.
    */
   int
   getCVSpgmrMemoryUsageForDoubles() const
   {
      long int r1;
      long int r2;
      int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r1);
   }

   /**
    * Return the size (in integer words) of memory used for integers.
    */
   int
   getCVSpgmrMemoryUsageForIntegers() const
   {
      long int r1;
      long int r2;
      int ierr = CVodeGetWorkSpace(d_cvode_mem, &r1, &r2);
      CVODE_SAMRAI_ERROR(ierr);
      return static_cast<int>(r2);
   }

   /**
    * Print out all data members for this object.
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
    * Static integer constant describing the size of an output buffer of a
    * CVODE statistic.
    */
   static const int STAT_OUTPUT_BUFFER_SIZE;

   /*
    * Static member function for linkage with CVODE routines.
    */
   static int
   CVODERHSFuncEval(
      realtype t,
      N_Vector y,
      N_Vector y_dot,
      void* my_solver)
   {
      return ((CVODESolver *)my_solver)->getCVODEFunctions()->
             evaluateRHSFunction(t, SABSVEC_CAST(y), SABSVEC_CAST(y_dot));
   }

   /*
    * Static member functions for linkage with CVSpgmr routines.
    */
   static int
   CVSpgmrPrecondSet(
      realtype t,
      N_Vector y,
      N_Vector fy,
      int jok,
      booleantype* jcurPtr,
      realtype gamma,
      void* my_solver)
   {
      int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
         CVSpgmrPrecondSet(t,
            SABSVEC_CAST(y),
            SABSVEC_CAST(fy),
            jok,
            jcurPtr,
            gamma);
      return success;
   }

   static int
   CVSpgmrPrecondSolve(
      realtype t,
      N_Vector y,
      N_Vector fy,
      N_Vector r,
      N_Vector z,
      realtype gamma,
      realtype delta,
      int lr,
      void* my_solver)
   {
      int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
         CVSpgmrPrecondSolve(t,
            SABSVEC_CAST(y),
            SABSVEC_CAST(fy),
            SABSVEC_CAST(r),
            SABSVEC_CAST(z),
            gamma,
            delta,
            lr);
      return success;
   }

   static int
   CVODEProjEval(
      realtype t,
      N_Vector y,
      N_Vector corr ,
      realtype epsProj ,
      N_Vector err, void* my_solver)
   {
      int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
         applyProjection(t,
            SABSVEC_CAST(y),
            SABSVEC_CAST(corr),
            epsProj,
            SABSVEC_CAST (err));
      return success;
   }

   static int
   CVODEJTimesRHSFuncEval(
      realtype t,
      N_Vector y,
      N_Vector y_dot ,
      void* my_solver)
   {
      int success = ((CVODESolver *)my_solver)->getCVODEFunctions()->
         evaluateJTimesRHSFunction(t,
            SABSVEC_CAST(y),
            SABSVEC_CAST(y_dot));
      return success;
   }


   /*
    * Open CVODE log file, allocate main memory for CVODE and initialize
    * CVODE memory record.  CVODE is initialized based on current state
    * of solver parameter data members.  If any solver parameters have
    * changed since last initialization, this function will be automatically
    * invoked at next call to the solve() method.  Also, if NEWTON iteration
    * is specified, this method also initializes the CVSpgmr linear solver.
    *
    * @pre d_solution_vector != 0
    */
   void
   initializeCVODE();

   std::string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the CVODESolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * CVODE may need to be re-initialized (e.g., if the linear solver
    * changes, CVODE must change its memory record).  In this case,
    * the initializeCVODE() member function is invoked in the next
    * call to solve().
    */

   /*
    * Solution vector.
    */
   SundialsAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to CVODE
    * and CVSpgmr.
    */
   CVODEAbstractFunctions* d_cvode_functions;

   /*
    * CVODE memory record.
    */
   void* d_cvode_mem;                    // CVODE memory structure

   /*
    * Linear solver for preconditioning
    */
   SUNLinearSolver d_linear_solver;


   /*
    * CVODE log file information.
    */
   FILE* d_cvode_log_file;                   // CVODE message log file
   std::string d_cvode_log_file_name;        // CVODE log file name

   /*
    * ODE parameters.
    */
   double d_t_0;        // initial value for independent variable
   double d_user_t_f;   // user-specified final value for independent variable
   double d_actual_t_f; // actual final value of indep. variable after a step
   SundialsAbstractVector* d_ic_vector;

   /*
    * ODE integration parameters.
    */
   int d_linear_multistep_method;
   double d_relative_tolerance;
   bool d_use_scalar_absolute_tolerance;
   double d_absolute_tolerance_scalar;
   SundialsAbstractVector* d_absolute_tolerance_vector;
   int d_stepping_method;

   /*
    * Optional CVODE parameters.
    */
   int d_max_order;
   int d_max_num_internal_steps;
   int d_max_num_warnings;
   double d_init_step_size;
   double d_max_step_size;
   double d_min_step_size;
   /*
    * CVSpgmr parameters
    */
   int d_precondition_type;
   int d_gram_schmidt_type;
   int d_max_krylov_dim;
   double d_tol_scale_factor;

   /*
    * Boolean flag indicating whether CVODE needs initialization
    * when solver is called.
    */
   bool d_CVODE_needs_initialization;

   /*
    * Boolean flag indicating whether user-supplied preconditioner
    * routines are provided in the concrete subclass of
    * CVODEAbstractFunctions.
    */
   bool d_uses_preconditioner;

   /*
    * Boolean flag indicating whether a user-supplied projection
    * routine is provided in the concrete subclass of
    * CVODEAbstractFunctions.
    */
   bool d_uses_projectionfn;

   /*
    * Boolean flag indicating whether a different RHS
    * routine for Jacobian -vector products is provided
    * in the concrete subclass of CVODEAbstractFunctions.
    */ 
   bool d_uses_jtimesrhsfn;
};

}
}

#endif
#endif
