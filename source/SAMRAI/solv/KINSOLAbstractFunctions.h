/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to user-specified functions for KINSOL package
 *
 ************************************************************************/

#ifndef included_solv_KINSOLAbstractFunctions
#define included_solv_KINSOLAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT KINSOL
 ************************************************************************
 */
#ifdef HAVE_SUNDIALS

namespace SAMRAI {
namespace solv {

/**
 * Class KINSOLAbstractFunctions is an abstract base class that defines
 * an interface for user-supplied functions to be used with KINSOL via
 * the C++ wrapper class KINSOLSolver.  To use KINSOL with the
 * C++ wrapper one must derive a subclass of this base class and pass it
 * into the KINSOLSolver constructor.  The pure virtual member
 * functions in this interface are used by KINSOL during the nonlinear
 * system solution process.  The complete argument lists in the function
 * signatures defined by KINSOL have been preserved for the user-supplied
 * routines have been preserved for the most part.  In a few cases, some
 * arguments do not appear in the function signatures below since they
 * are superfluous via this interface.
 *
 * KINSOL only requires that the function evaluateNonlinearFunction()
 * be supplied.  The other virtual functions are optional in KINSOL.
 * Note that the use of the optional functions may be turned on and off
 * via boolean arguments to the constructor of the KINSOLSolver
 * class, or using the setKINSOLFunctions() member function of that class.
 *
 * @see KINSOLSolver
 * @see SundialsAbstractVector
 */

class KINSOLAbstractFunctions
{
public:
   /**
    * Uninteresting constructor destructor for KINSOLAbstractFunctions.
    */
   KINSOLAbstractFunctions();
   virtual ~KINSOLAbstractFunctions();

   /**
    * User-supplied nonlinear residual function evaluation.
    *
    * The function arguments are:
    *
    *
    *
    * - \b soln   (INPUT) {current iterate for the nonlinear solution.}
    * - \b fval   (OUTPUT){current value of residual function.}
    *
    *
    *
    *
    * IMPORTANT: This function must not modify the vector soln.
    */
   virtual void
   evaluateNonlinearFunction(
      SundialsAbstractVector* soln,
      SundialsAbstractVector* fval) = 0;

   /**
    * User-supplied preconditioner setup function.  The setup function
    * is called to provide matrix data for the subsequent call(s) to
    * precondSolve().  That is, this preconditioner setup function
    * is used to evaluate and preprocess any Jacobian-related data
    * needed by the preconditioner solve function.   The integer return
    * value is a flag indicating success if 0 is returned, and failure
    * otherwise.  If a non-zero value is returned, KINSOL stops.  Together
    * precondSetup() and precondSolve() form a right preconditoner for the
    * KINSOL Krylov solver.  This function will not be called prior to
    * every call of precondSolve(), but instead will be called only as
    * often as needed to achieve convergence within the Newton iteration.
    *
    * The function arguments are:
    *
    *
    *
    * - \b soln          (INPUT) {current iterate for the nonlinear solution}
    * - \b fval          (INPUT) {current values of the nonlinear residual}
    * - \b soln_scale   (INPUT) {diagonal entries of the nonlinear solution
    *                               scaling matrix}
    * - \b fval_scale   (INPUT) {diagonal entries of the nonlinear residual
    *                               scaling matrix}
    * - \b num_feval    (OUTPUT){number of nonlinear function evaluations
    *                               made to approximate the Jacobian, if any.
    *                               For example, if the routine evaluates the
    *                               function twice, num_feval is set to 2}
    *
    *
    *
    *
    * The scaling vectors are provided for possible use in approximating
    * Jacobian data; e.g., uing difference quotients.  The
    *
    * IMPORTANT: This function must not modify the vector arguments.
    */
   virtual int
   precondSetup(
      SundialsAbstractVector* soln,
      SundialsAbstractVector* soln_scale,
      SundialsAbstractVector* fval,
      SundialsAbstractVector* fval_scale,
      int& num_feval) = 0;

   /**
    * User-supplied preconditioner solve function.  This function must
    * solve \f$P x = r\f$, where \f$P\f$ is the right preconditioner matrix formed
    * by precondSetup().  The integer return value is a flag indicating
    * success if 0 is returned, and failure otherwise.  If a non-zero
    * value is returned, KINSOL stops.
    *
    * The function arguments are:
    *
    *
    *
    * - \b soln          (INPUT) {current iterate for the nonlinear solution}
    * - \b fval          (INPUT) {current iterate for the nonlinear residual}
    * - \b soln_scale   (INPUT) {diagonal entries of the nonlinear solution
    *                               scaling matrix}
    * - \b fval_scale   (INPUT) {diagonal entries of the nonlinear residual
    *                               scaling matrix}
    * - \b rhs           (OUTPUT){rhs-side (r) on input and must be set to
    *                              preconditioner solution (i.e., x) on output}
    * - \b num_feval    (OUTPUT){number of nonlinear function evaluations
    *                               made to approximate the Jacobian, if any.
    *                               For example, if the routine evaluates the
    *                               function twice, num_feval is set to 2}
    *
    *
    *
    *
    * IMPORTANT: This function must not modify soln, fval, or the scaling
    *            vectors.
    */
   virtual int
   precondSolve(
      SundialsAbstractVector* soln,
      SundialsAbstractVector* soln_scale,
      SundialsAbstractVector* fval,
      SundialsAbstractVector* fval_scale,
      SundialsAbstractVector* rhs,
      int& num_feval) = 0;

   /**
    * Optional user-supplied A times x routine, where A is an approximation
    * to the Jacobian matrix and v is some vector.  product = (A * vector)
    * is computed.
    *
    * The function arguments are:
    *
    *
    *
    * - \b vector       (INPUT) {the vector multiplied by the Jacobian}
    * - \b product      (OUTPUT){product of the Jacobian and vector; \f$A v\f$)}
    * - \b new_soln    (INPUT) {flag indicating whether solution has changed
    *                              since last call to this routine.
    *                              For example, if this routine computes and
    *                              saves the Jacobian, then the Jacobian does
    *                              not require computation if flag is false.}
    * - \b soln         (INPUT) {current iterate for the nonlinear solution}
    *
    *
    *
    *
    * IMPORTANT: This function must not modify soln vector.
    */
   virtual int
   jacobianTimesVector(
      SundialsAbstractVector* vector,
      SundialsAbstractVector* product,
      const bool soln_changed,
      SundialsAbstractVector* soln) = 0;

};

}
}
#endif

#endif
