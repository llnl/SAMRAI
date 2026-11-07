/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to user-specified functions for CVODE package
 *
 ************************************************************************/

#ifndef included_solv_CVODEAbstractFunctions
#define included_solv_CVODEAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
namespace solv {

/**
 * Class CVODEAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner
 * routines to be used with CVODE and CVSpgmr via the C++ wrapper
 * class CVODESolver.  To use CVODE with the C++ wrapper one must
 * derive a subclass of this base class and pass it into the CVODESolver
 * constructor.  The pure virtual member functions in this interface are
 * used by CVODE and CVSpgmr during the ODE integration process.  The
 * complete argument lists in the function signatures defined by CVODE
 * for the user-supplied routines have been preserved for the most part.
 * In a few cases, some arguments do not appear in the function signatures
 * below since they are superfluous via this interface.
 *
 * @see CVODESolver
 * @see SundialsAbstractVector
 */

class CVODEAbstractFunctions
{
public:
   /**
    * The constructor and destructor for CVODEAbstractFunctions
    * is empty.
    */
   CVODEAbstractFunctions();
   virtual ~CVODEAbstractFunctions();

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
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int
   evaluateRHSFunction(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* y_dot) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int
   CVSpgmrPrecondSet(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* fy,
      int jok,
      int* jcurPtr,
      double gamma) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int
   CVSpgmrPrecondSolve(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* fy,
      SundialsAbstractVector* r,
      SundialsAbstractVector* z,
      double gamma,
      double delta,
      int lr) = 0;

   /*!
    * @brief User-supplied function to project the current solution and
    * estimated error onto the constraint manifold.
    *
    * The function arguments are:
    *
    * - \b t (INPUT) {current value of the independent variable}
    * - \b y (INPUT) {current value of the dependent variable vector}
    * - \b corr (OUTPUT) {correction such that y+corr is on the constraint
    *                     manifold}
    * - \b epsProj (INPUT) {WRMS norm tolerance for a nonlinear solver
    *                       iteration}
    * - \b err (INPUT/OUTPUT) {input: unprojected error , output: projected
    *                          error}
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int
   applyProjection(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* corr,
      double epsProj,
      SundialsAbstractVector* err) = 0;

   /*!
    * @brief User-supplied right-hand side function used in evaluating finite
    * difference Jacobian -vector products.
    *
    * The function arguments are:
    *
    * - \b t (INPUT) {current value of the independent variable}
    * - \b y (INPUT) {current value of the dependent variable vector}
    * - \b ydot (OUTPUT) {current value of the derivative of y}
    * 
    * IMPORTANT: This function must not modify the vector y
    */
   virtual int evaluateJTimesRHSFunction(
      double t,
      SundialsAbstractVector* y,
      SundialsAbstractVector* y_dot) = 0;

};

}
}

#endif

#endif
