/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to user-defined operations used in FAC solve.
 *
 ************************************************************************/
#ifndef included_solv_FACOperatorStrategy
#define included_solv_FACOperatorStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/SAMRAIVectorReal.h"


namespace SAMRAI {
namespace solv {

/*!
 * @brief Defines a Strategy pattern interface to problem-specific
 * operations needed to implement the FAC preconditioner algorithm.
 *
 * The FACPreconditioner constructor accepts a concrete
 * implementation of this interface and calls the concrete implementations
 * of the virtual functions declared herein during the solution
 * process.
 *
 * All vector arguments in these interfaces are guaranteed to be
 * either the vectors given to in FACPreconditioner::solveSystem()
 * or FACPreconditioner::initializeSolverState()
 * or vectors cloned from them.
 *
 * @see FACPreconditioner
 */

class FACOperatorStrategy
{
public:
   /*!
    * @brief Empty constructor.
    */
   FACOperatorStrategy();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~FACOperatorStrategy();

   //@{
   /*!
    * @name Operator-dependent virtual methods
    */

   /*!
    * @brief Restrict the solution quantity to the specified level
    * from the next finer level.
    *
    * Restrict the residual data to level dest_ln in the destination
    * vector d, from level dest_ln+1 in the source vector s.
    *
    * Can assume:
    * -# dest_ln is not the finest level in the range being solved.
    * -# corresponding solution has been computed on level dest_ln+1.
    * -# the source and destination residual vectors may or may not
    *      be the same.  (This method must work in either case.)
    *
    * Upon return from this function, the solution on the refined region
    * of the coarse level will represent the coarsened version of the
    * fine solution in a manner that is consistent with the linear system
    * approximation on the composite grid.  This function must not change
    * the solution values anywhere except on level dest_ln of the destination
    * vector.
    *
    * The source and destination vectors may be the same.
    *
    * @param source source solution
    * @param dest destination solution
    * @param dest_ln destination level number
    */
   virtual void
   restrictSolution(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln) = 0;

   /*!
    * @brief Restrict the residual quantity to the specified level
    * from the next finer level.
    *
    * Restrict the residual data to level dest_ln in the destination
    * vector d, from level dest_ln+1 in the source vector s.
    *
    * Can assume:
    * -# dest_ln is not the finest level in the range being solved.
    * -# correspnding residual has been computed on level dest_ln+1.
    * -# the source and destination residual vectors may or may not
    *      be the same.  (This method must work in either case.)
    *
    * Upon return from this function, the residual on the refined region
    * of the coarse level will represent the coarsened version of the
    * fine residual in a manner that is consistent with the linear system
    * approximation on the composite grid.  This function must not change
    * the residual values anywhere except on level dest_ln of the destination
    * vector.
    *
    * The source and destination vectors may be the same.
    *
    * @param source source residual
    * @param dest destination residual
    * @param dest_ln destination level number
    */
   virtual void
   restrictResidual(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln) = 0;

   /*!
    * @brief Prolong the error quantity to the specified level
    * from the next coarser level and apply the correction to
    * the fine-level error.
    *
    * On the part of the coarse level that does @em not overlap
    * the fine level, the error is the corection to Au=f.
    *
    * On the part of the coarse level that @em does overlap the fine level,
    * the error is the corection to Ae=r of the fine level.
    *
    * This function should apply the coarse-level correction to the
    * fine level, that is @f[
    * e^{fine} \leftarrow e^{fine} + I^{fine}_{coarse} e^{coarse}
    * @f]
    *
    * @b Note: You probably have to store the refined error in a
    * temporary location before adding it to the current error.
    *
    * The array of boundary information contains a description of
    * the coarse-fine level boundary for each patch on the level;
    * the boundary information for patch N is obtained as the
    * N-th element in the array, coarse_fine_boundary[N].
    *
    * Upon return from this function,
    * the error on the fine level must represent the correction
    * to the solution on that level.
    * Also, this function must not change the error values on the
    * coarse level.
    *
    * The source and destination vectors may be the same.
    *
    * @param source source error vector
    * @param dest destination error vector
    * @param dest_ln destination level number of data transfer
    */
   virtual void
   prolongErrorAndCorrect(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln) = 0;

   /*!
    * @brief Perform a given number of relaxations on the error.
    *
    * Relax the residual equation Ae=r by applying the
    * given number of smoothing sweeps on the specified level.
    * The relaxation may ignore the possible existence of finer
    * levels on a given level.
    *
    * The array of boundary
    * information contains a description of the coarse-fine level boundary
    * for each patch on the level; the boundary information for patch N is
    * obtained as the N-th element in the array, coarse_fine_boundary[N].
    *
    * May assume:
    * - If intermediate data from level l+1 is needed (for example,
    * to match flux at coarse-fine boundaries), that data is already
    * computed and stored on level l+1.
    * - The error in the next finer level has been computed and stored
    * there.
    *
    * Steps for each iteration.
    * -# Fill ghost boundaries
    * -# Compute intermediate data (if needed) and coarsen intermediate
    * data stored in level l+1 (if needed).
    * -# Perform relaxation step (update e toward a better
    * approximation).
    *
    * Final step before leaving function.
    * - If needed, compute and store intermediate data for
    * next coarser level l-1.
    *
    * @param error error vector
    * @param residual residual vector
    * @param ln level number
    * @param num_sweeps number of sweeps
    */
   virtual void
   smoothError(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int ln,
      int num_sweeps) = 0;

   /*!
    * @brief Solve the residual equation Ae=r on the coarsest
    * level in the FAC iteration.
    *
    * Here e is the given error quantity and r is the
    * given residual quantity.  The array of boundary information contains a
    * description of the coarse-fine level boundary for each patch on the
    * level; the boundary information for patch N is obtained as the N-th
    * element in the array, coarse_fine_boundary[N].
    *
    * This routine must fill boundary values for given solution quantity
    * on all patches on the specified level before the solve is performed.
    *
    * @param error error vector
    * @param residual residual vector
    * @param coarsest_ln coarsest level number
    * @return 0 if solver converged to specified level, nonzero otherwise.
    *
    */
   virtual int
   solveCoarsestLevel(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int coarsest_ln) = 0;

   /*!
    * @brief Compute composite grid residual on a single level.
    *
    * For the specified level number ln,
    * compute the @em composite residual r=f-Au,
    * where f is the right hand side and u is the solution.
    * Note that the composite residual is not a one-level
    * residual.  It must take into account the composite grid
    * stencil around the coarse-fine grid interface.
    *
    * May assume:
    * - Composite residual on next finer level l+1,
    *   has been computed already.
    * - If any intermediately computed data is needed from
    *   level l+1, it has been done and stored on that level.
    * - Residual computations for the original equation and
    *   the error equations will not be intermingled within
    *   one FAC cycle.
    *
    * Steps:
    * -# Fill boundary ghosts.
    * -# If needed, coarsen intermediate data from level l+1.
    * -# Compute residual @f$ r^l \leftarrow f - A u^l @f$.
    *
    * Final step before leaving function:
    * - If any intermediately computed data is needed in at
    *   level l-1, it must be computed and stored before
    *   leaving this function.
    *
    * @b Important: Do not restrict residual from finer levels.
    * (However, you must write the function restrictResidual()
    * to do this.)
    *
    * @b Important: This function must also work when the
    * right-hand-side and the residual are identical.
    * In that case, it should effectively do @f$ r \leftarrow r - A u @f$.
    *
    * @param residual residual vector
    * @param solution solution vector
    * @param rhs source (right hand side) vector
    * @param ln level number
    * @param error_equation_indicator flag stating whether u is an error
    * vector or a solution vector
    */
   virtual void
   computeCompositeResidualOnLevel(
      SAMRAIVectorReal<double>& residual,
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs,
      int ln,
      bool error_equation_indicator) = 0;

   /*!
    * @brief Compute the norm of the residual quantity
    *
    * Compute norm of the given residual on the given range of
    * hierarchy levels.  The residual vector is computed already
    * and you should @b not change it.
    * The only purpose of this function to allow you to choose
    * how to define the norm.
    *
    * The norm value is used during the FAC iteration
    * to determine convergence of the composite grid linear system.
    *
    * Residual values that lie under a finer level should not be counted.
    *
    * @param residual residual vector
    * @param fine_ln finest level number
    * @param coarse_ln coarsest level number
    *
    * @return norm value of residual vector, which should be non-negative
    */
   virtual double
   computeResidualNorm(
      const SAMRAIVectorReal<double>& residual,
      int fine_ln,
      int coarse_ln) = 0;

   /*!
    * @brief Regular call back routine to be called after each FAC cycle.
    *
    * This function is called after each FAC cycle.
    * It allows you to monitor the progress and do other things.
    * You should @em not modify the solution vector in the argument.
    *
    * The default implementation does nothing.
    *
    * @param fac_cycle_num FAC cycle number completed
    * @param current_soln current solution
    * @param residual residual based on the current solution
    */
   virtual void
   postprocessOneCycle(
      int fac_cycle_num,
      const SAMRAIVectorReal<double>& current_soln,
      const SAMRAIVectorReal<double>& residual) = 0;

   /*!
    * @brief Compute hierarchy-dependent data if any is required
    *
    * This function is called when the hierarchy configuration changes.
    * If you maintain any hierarchy-dependent data in your implementation
    * (for example, caching communication schedules or computing
    * coarse-fine boundaries),
    * use this function to update that data.
    *
    * If you do not maintain such data, this function may be empty.
    *
    * Note that although the vector arguments given to other
    * methods in this class may not necessarily be the same
    * as those given to this method, there will be similarities,
    * including:
    * - hierarchy configuration (hierarchy pointer and level range)
    * - number, type and alignment of vector component data
    * - ghost cell width of data in the solution (or solution-like) vector
    *
    * @param solution solution vector u
    * @param rhs right hand side vector f
    *
    * The default implementation does nothing.
    */
   virtual void
   initializeOperatorState(
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs) = 0;

   /*!
    * @brief Remove all hierarchy-dependent data.
    *
    * Remove all hierarchy-dependent data set by initializeOperatorState().
    *
    * @see initializeOperatorState
    */
   virtual void
   deallocateOperatorState();

   //@}

};

}
}

#endif
