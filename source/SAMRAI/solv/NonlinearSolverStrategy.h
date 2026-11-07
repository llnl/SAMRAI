/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface between implicit integrator and nonlinear solver.
 *
 ************************************************************************/

#ifndef included_solv_NonlinearSolverStrategy
#define included_solv_NonlinearSolverStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/SAMRAIVectorReal.h"

#include <memory>

namespace SAMRAI {
namespace solv {

/**
 * @brief Abstract base class defining interface between an
 * algs::ImplicitIntegrator object and
 * a nonlinear solver used to advance the solution in time.
 *
 * The interface follows the Strategy design pattern.
 * The methods declared in the interface are provided in a
 * concrete solver derived from this base class.
 *
 * @see algs::ImplicitIntegrator
 */

class NonlinearSolverStrategy
{
public:
   /**
    * Empty constructor for algs::NonlinearSolverStrategy.
    */
   NonlinearSolverStrategy();

   /**
    * Empty constructor for algs::NonlinearSolverStrategy.
    */
   virtual ~NonlinearSolverStrategy();

   /**
    * Initialize the solver state.  The vector argument represents the
    * solution of the nonlinear system.  In general, this routine must
    * be called before the solve() routine is invoked.
    */
   virtual void
   initialize(
      const std::shared_ptr<SAMRAIVectorReal<double> >& solution) = 0;

   /**
    * Solve the nonlinear problem and return the integer code defined by the
    * particular nonlinear solver package in use (e.g., indicating success
    * or failure of solution process).  In general, the initialize() routine
    * must be called before this solve function.
    */
   virtual int
   solve() = 0;

};

}
}

#endif
