/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Poisson using FAC
 *
 ************************************************************************/
#ifndef included_solv_CellPoissonFACOps
#define included_solv_CellPoissonFACOps

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/CellPoissonHypreSolver.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/FACOperatorStrategy.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace solv {

/*!
 * @brief FAC operator class to solve Poisson's equation on a SAMR grid,
 * using cell-centered, second-order finite-volume method, with Robin
 * boundary conditions.
 *
 * This class provides operators that are used by the FAC
 * preconditioner FACPreconditioner.
 * It is used to solve the scalar Poisson's equation using a cell-centered
 * second-order finite-volume discretization.
 * It is designed to provide all operations specific to
 * the scalar Poisson's equation,
 * @f[ \nabla \cdot D \nabla u + C u = f @f]
 * (see PoissonSpecifications) where
 * - C, D and f are indpendent of u
 * - C is a cell-centered scalar field
 * - D is the @em diffusion @em coefficients, stored on faces
 * - f is a cell-centered scalar function
 *
 * You are left to provide the source function, initial guess, etc.,
 * by specifying them in specific forms.
 *
 * This class provides:
 * -# 5-point (second order), cell-centered stencil operations
 *    for the discrete Laplacian.
 * -# Red-black Gauss-Seidel smoothing.
 * -# Provisions for working Robin boundary conditions
 *    (see RobinBcCoefStrategy).
 *
 * This class is meant to provide the Poisson-specific operator
 * used by the FAC preconditioner, FACPreconditioner.
 * To use the preconditioner with this class, you will have to provide:
 * -# The solution vector SAMRAIVectorReal,
 *    with appropriate norm weighting for the cell-centered AMR mesh.
 *    This class provides the function computeVectorWeights()
 *    to help with computing the appropriate weights.
 *    Since this is for a scalar equation, only the first depth
 *    of the first component of the vectors are used.
 *    All other parts are ignored.
 * -# The source vector SAMRAIVectorReal for f.
 * -# A PoissonSpecifications objects to specify
 *    the cell-centered scalar field C and the side-centered
 *    diffusion coefficients D
 * -# The boundary condition specifications in terms of the coefficients
 *    @f$ \alpha @f$, @f$ \beta @f$ and @f$ \gamma @f$ in the
 *    Robin formula @f$  \alpha u + \beta u_n = \gamma @f$ applied on the
 *    boundary faces.  See RobinBcCoefStrategy.
 *
 * This class allocates and deallocates only its own scratch data.
 * Other data that it manipuates are passed in as function arguments.
 * Hence, it owns none of the solution vectors, error vectors,
 * diffusion coefficient data, or any such things.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b    coarse_solver_choice
 *
 *    - \b    coarse_solver_tolerance
 *
 *    - \b    coarse_solver_max_iterations
 *
 *    - \b    cf_discretization
 *
 *    - \b    prolongation_method
 *
 *    - \b    enable_logging
 *
 * <b> Details:</b> <br>
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
 *     <td>coarse_solver_choice</td>
 *     <td>string</td>
 *     <td>"hypre"</td>
 *     <td>"hypre", "redblack", "jacobi"</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coarse_solver_tolerance</td>
 *     <td>double</td>
 *     <td>1e-14</td>
 *     <td>>0.0</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coarse_solver_max_iterations</td>
 *     <td>int</td>
 *     <td>10</td>
 *     <td>>=1</td>
 *     <td>opt
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>cf_discretization</td>
 *     <td>string</td>
 *     <td>"Ewing"</td>
 *     <td>"Ewing", "CONSTANT_REFINE", "CONSERVATIVE_LINEAR_REFINE", "LINEAR_REFINE"</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>prolongation_method</td>
 *     <td>string</td>
 *     <td>"LINEAR_REFINE"</td>
 *     <td>"CONSTANT_REFINE", "CONSERVATIVE_LINEAR_REFINE", "LINEAR_REFINE"</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>enable_logging</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Not written to restart.  Value in input db used.</td>
 *   </tr>
 * </table>
 *
 */
class CellPoissonFACOps:
   public FACOperatorStrategy
{

public:
#ifdef HAVE_HYPRE
   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    * @param dim
    * @param object_name Object name
    * @param hypre_solver underlying hypre solver
    * @param input_db Input database
    */
   CellPoissonFACOps(
      const std::shared_ptr<CellPoissonHypreSolver>& hypre_solver,
      const tbox::Dimension& dim,
      const std::string& object_name = std::string(),
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());
#else
   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    * @param dim
    * @param object_name Object name
    * @param input_db Input database
    */
   CellPoissonFACOps(
      const tbox::Dimension& dim,
      const std::string& object_name = std::string(),
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());
#endif

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~CellPoissonFACOps();

   /*!
    * @brief Set the scalar Poisson equation specifications.
    */
   void
   setPoissonSpecifications(
      const PoissonSpecifications& spec)
   {
      d_poisson_spec = spec;
   }

   //@{
   /*!
    * @name Functions for setting patch data indices and coefficients
    */

   /*!
    * @brief Set the scratch patch data index for the flux.
    *
    * The use of this function is optional.
    * The patch data index should be a pdat::SideData<TYPE> type of variable.
    * If the flux id is -1 (the default initial value), scratch space
    * for the flux is allocated as needed and immediately deallocated
    * afterward, level by level.  If you have space preallocated for
    * flux and you would like that to be used, set flux id to the
    * patch data index of that space.
    */
   void
   setFluxId(
      int flux_id)
   {
      d_flux_id = flux_id;
#ifdef DEBUG_CHECK_ASSERTIONS
      checkInputPatchDataIndices();
#endif
   }

   //@}

   /*!
    * @brief Provide an implementation for getting the
    * physical bc coefficients
    *
    * If your solution is fixed at the physical boundary
    * ghost cell centers AND those cells have the correct
    * values before entering solveSystem(), you may use a
    * GhostCellRobinBcCoefs object.
    *
    * If your solution is @b not fixed at the ghost cell centers,
    * the ghost cell values will change as the interior
    * cell values change.  In those cases, the flexible
    * Robin boundary conditions are applied.  You must
    * call this function to provide the implementation for
    * determining the boundary condition coefficients.
    *
    * @param physical_bc_coef Pointer to an object that can
    *        set the Robin bc coefficients.
    */
   void
   setPhysicalBcCoefObject(
      const RobinBcCoefStrategy* physical_bc_coef)
   {
      d_physical_bc_coef = physical_bc_coef;
      d_bc_helper.setCoefImplementation(physical_bc_coef);
#ifdef HAVE_HYPRE
      d_hypre_solver->setPhysicalBcCoefObject(d_physical_bc_coef);
#endif
   }

   //@{

   /*!
    * @name Functions for checking validity and correctness of state.
    */

   /*!
    * @brief Check validity and correctness of input patch data indices.
    *
    * Descriptors checked:
    * -# Diffusion coefficient (see setDiffcoefId())
    * -# Flux (see setFluxId())
    * -# Source (see setScalarFieldId())
    */
   void
   checkInputPatchDataIndices() const;

   //@}

   /*!
    * @brief Set weight appropriate for computing vector norms.
    *
    * If you this function to set the weights used when you
    * SAMRAIVectorReal::addComponent, you can use the
    * vector norm functions of SAMRAIVectorReal, and
    * the weights will be used to blank out coarse grid
    * regions under fine grids.
    *
    * The weights computed are specific to the cell-centered
    * discretization used by this class.  The weight is equal
    * to the cell volume if the cell has not been refined,
    * and zero if it has.
    *
    * This function is state-independent.  All inputs are in
    * the argument list.
    *
    * @param hierarchy Hierarchy configuration to compute weights for
    * @param weight_id hier::Patch data index of the weight
    * @param coarsest_ln Coarsest level number.  Must be included
    *        in hierarchy.  Must not be greater than @c finest_ln.
    *        Default to 0.
    * @param finest_ln Finest level number.  Must be included
    *        in hierarchy.  Must not be less than @c coarsest_ln.
    *        Default to finest level in @c hierarchy.
    *
    * @pre hierarchy
    * @pre d_dim == hierarchy->getDim()
    * @pre finest_ln >= coarsest_ln
    */
   void
   computeVectorWeights(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int weight_id,
      int coarsest_ln = -1,
      int finest_ln = -1) const;

   /*!
    * @brief Set the FAC preconditioner that will be using this object.
    *
    * The FAC preconditioner is accessed to get convergence data during
    * the cycle postprocessing step.  It is optional.
    */
   void
   setPreconditioner(
      const FACPreconditioner* preconditioner)
   {
      d_preconditioner = preconditioner;
   }

   /*!
    * @brief function to compute flux, using general diffusion
    * coefficient data.
    *
    * Recall that this solver class discretizes the PDE
    * @f[ \nabla \cdot D \nabla u + C u = f @f] on an AMR grid.  This member
    * function allows users of this solver class to compute gradient
    * terms, @f[ D \nabla w @f], in their code in a manner consistent with the
    * solver discretization.   In particular, when solving PDE systems, it may
    * be necessary to discretize the gradient operator appearing in equations
    * not treated by the solver class in the same way as those treated by this
    * class.  These funtions allow users to do this easily.  The divergence
    * operator used in this solver is the standard sum of centered differences
    * involving flux terms on the cell sides computed by these routines.
    *
    * Note that the patch must exist on a level in an AMR hierarchy so that
    * the discretization can be computed properly at the coarse-fine interface.
    * Poisson coefficients C and D must exist on the patch, if they are variable.
    * Also, calling this function does not affect the internal solver state in any
    * way.  However, the solver must be fully initialized before it is called and care
    * should be exercised to pass arguments so that the solver solution quantity and
    * other internal solver quantities are not adversely affected.
    *
    * @param patch patch on which computation will take place
    * @param ratio_to_coarser_level refinement ratio from coarser level to level
    *                               on which patch lives; if current patch level
    *                               is level zero, this is ignored
    * @param w_data cell-centered data
    * @param Dgradw_data side-centered flux data (i.e., D (grad w))
    *
    * @pre (d_dim == patch.getDim()) &&
    *      (d_dim == ratio_to_coarser_level.getDim()) &&
    *      (d_dim == w_data.getDim()) && (d_dim == Dgradw_data.getDim())
    * @pre patch.inHierarchy()
    * @pre w_data.getGhostCellWidth() >= hier::IntVector::getOne(ratio_to_coarser_level.getDim())
    */
   void
   computeFluxOnPatch(
      const hier::Patch& patch,
      const hier::IntVector& ratio_to_coarser_level,
      const pdat::CellData<double>& w_data,
      pdat::SideData<double>& Dgradw_data) const;

   //@{ @name FACOperatorStrategy virtuals

   virtual void
   restrictSolution(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln);
   virtual void
   restrictResidual(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln);

   virtual void
   prolongErrorAndCorrect(
      const SAMRAIVectorReal<double>& source,
      SAMRAIVectorReal<double>& dest,
      int dest_ln);

   virtual void
   smoothError(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int ln,
      int num_sweeps);

   virtual int
   solveCoarsestLevel(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int coarsest_ln);

   virtual void
   computeCompositeResidualOnLevel(
      SAMRAIVectorReal<double>& residual,
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs,
      int ln,
      bool error_equation_indicator);

   virtual double
   computeResidualNorm(
      const SAMRAIVectorReal<double>& residual,
      int fine_ln,
      int coarse_ln);

   virtual void
   initializeOperatorState(
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs);

   virtual void
   deallocateOperatorState();

   virtual void
   postprocessOneCycle(
      int fac_cycle_num,
      const SAMRAIVectorReal<double>& current_soln,
      const SAMRAIVectorReal<double>& residual);

   //@}

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

protected:
   /*!
    * @brief Read parameters from input database.
    *
    * @param input_db Input Database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

private:
   // Internals of both constructors
   void
   buildObject(
      const std::shared_ptr<tbox::Database>& input_db);

   //@{
   /*!
    * @name Private workhorse functions.
    */

   /*!
    * @brief Red-black Gauss-Seidel error smoothing on a level.
    *
    * Smoothes on the residual equation @f$ Ae=r @f$ on a level.
    *
    * @param error error vector
    * @param residual residual vector
    * @param ln level number
    * @param num_sweeps number of sweeps
    * @param residual_tolerance the maximum residual considered to be
    *        converged
    *
    * @pre data.getPatchHierarchy() == d_hierarchy &&
    *      residual.getPatchHierarchy() == d_hierarchy
    */
   void
   smoothErrorByRedBlack(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int ln,
      int num_sweeps,
      double residual_tolerance = -1.0);

   /*!
    * @brief Solve the coarsest level using HYPRE
    */
   int
   solveCoarsestLevel_HYPRE(
      SAMRAIVectorReal<double>& error,
      const SAMRAIVectorReal<double>& residual,
      int ln);

   /*!
    * @brief Fix flux per Ewing's coarse-fine boundary treatment.
    *
    * Ewing's coarse-fine boundary treatment can be implemented
    * using a constant refinement into the fine-grid ghost boundary,
    * naively computing the flux using the constant-refined data then
    * fixing up the flux to correct the error.
    *
    * To use this function
    * -# you must use constant refinement to fill the fine level ghost cells
    * -# the flux must first be computed and stored
    *
    * @param patch patch
    * @param soln_data cell-centered solution data
    * @param flux_data side-centered flux data
    * @param ratio_to_coarser Refinement ratio to the next coarser level.
    *
    * @pre (d_dim == patch.getDim()) && (d_dim == soln_data.getDim()) &&
    *      (d_dim == flux_data.getDim()) &&
    *      (d_dim == ratio_to_coarser.getDim())
    */
   void
   ewingFixFlux(
      const hier::Patch& patch,
      const pdat::CellData<double>& soln_data,
      pdat::SideData<double>& flux_data,
      const hier::IntVector& ratio_to_coarser) const;

   /*!
    * @brief AMR-unaware function to compute residual on a single patch,
    * with variable scalar field.
    *
    * @param patch patch
    * @param flux_data side-centered flux data
    * @param soln_data cell-centered solution data
    * @param rhs_data cell-centered rhs data
    * @param residual_data cell-centered residual data
    *
    * @pre (d_dim = patch.getDim()) && (d_dim == flux_data.getDim()) &&
    *      (d_dim == soln_data.getDim()) && (d_dim == rhs_data.getDim()) &&
    *      (d_dim == residual_data.getDim())
    */
   void
   computeResidualOnPatch(
      const hier::Patch& patch,
      const pdat::SideData<double>& flux_data,
      const pdat::CellData<double>& soln_data,
      const pdat::CellData<double>& rhs_data,
      pdat::CellData<double>& residual_data) const;

   /*!
    * @brief AMR-unaware function to red or black smoothing on a single patch,
    * for variable diffusion coefficient and variable scalar field.
    *
    * @param patch patch
    * @param flux_data side-centered flux data
    * @param rhs_data cell-centered rhs data
    * @param soln_data cell-centered solution data
    * @param red_or_black red-black switch.  Set to 'r' or 'b'.
    * @param p_maxres max residual output.  Set to NULL to avoid computing.
    *
    * @pre (d_dim == patch.getDim()) && (d_dim == flux_data.getDim()) &&
    *      (d_dim == soln_data.getDim()) && (d_dim == rhs_data.getDim())
    * @pre (red_or_black == 'r') || (red_or_black == 'b')
    */
   void
   redOrBlackSmoothingOnPatch(
      const hier::Patch& patch,
      const pdat::SideData<double>& flux_data,
      const pdat::CellData<double>& rhs_data,
      pdat::CellData<double>& soln_data,
      char red_or_black,
      double* p_maxres = 0) const;

   //@}

   //@{ @name For executing, caching and resetting communication schedules.

   /*!
    * @brief Execute a refinement schedule
    * for prolonging cell data.
    *
    * General notes regarding internal objects for communication:
    * We maintain objects to support caching schedules to improve
    * efficiency.  Communication is needed in 5 distinct tasks.
    *   -# Prolongation
    *   -# Restriction
    *   -# Flux coarsening.  Changing the coarse grid flux to the
    *      composite grid flux by coarsening the fine grid flux
    *      at the coarse-fine boundaries.
    *   -# Fill boundary data from other patches in the same level
    *      and physical boundary condition.
    *   -# Fill boundary data from same level, coarser levels
    *      and physical boundary condition.
    *
    * For each task, we maintain a refine or coarsen operator,
    * and a array of communication schedules (one for each
    * destination level).
    *
    * The 5 member functions named @c xeqSchedule... execute
    * communication schedules appropriate for five specific tasks.
    * They use a cached schedule if possible or create and cache
    * a new schedule if needed.  These functions and the data
    * they manipulate are as follows:
    * <ol>
    *   <li> xeqScheduleProlongation():
    *        d_prolongation_refine_operator
    *        d_prolongation_refine_schedules
    *   <li> xeqScheduleURestriction():
    *        d_restriction_coarsen_operator,
    *        d_urestriction_coarsen_schedules.
    *   <li> xeqScheduleRRestriction():
    *        d_restriction_coarsen_operator,
    *        d_rrestriction_coarsen_schedules.
    *   <li> xeqScheduleFluxCoarsen():
    *        d_flux_coarsen_operator,
    *        d_flux_coarsen_schedules.
    *   <li> xeqScheduleGhostFill():
    *        d_ghostfill_refine_operator,
    *        d_ghostfill_refine_schedules.
    *   <li> xeqScheduleGhostFillNoCoarse():
    *        d_ghostfill_nocoarse_refine_operator,
    *        d_ghostfill_nocoarse_refine_schedules.
    * </ol>
    *
    * @return refinement schedule for prolongation
    *
    * @pre d_prolongation_refine_schedules[dest_ln]
    */
   void
   xeqScheduleProlongation(
      int dst_id,
      int src_id,
      int scr_id,
      int dest_ln);

   /*!
    * @brief Execute schedule for restricting solution to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for restriction
    *
    * @pre d_urestriction_coarsen_schedules[dest_ln]
    */
   void
   xeqScheduleURestriction(
      int dst_id,
      int src_id,
      int dest_ln);

   /*!
    * @brief Execute schedule for restricting residual to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for restriction
    *
    * @pre d_rrestriction_coarsen_schedules[dest_ln]
    */
   void
   xeqScheduleRRestriction(
      int dst_id,
      int src_id,
      int dest_ln);

   /*!
    * @brief Execute schedule for coarsening flux to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return coarsening schedule for setting composite grid flux at
    * coarse-fine boundaries.
    *
    * @pre d_flux_coarsen_schedules[dest_ln]
    */
   void
   xeqScheduleFluxCoarsen(
      int dst_id,
      int src_id,
      int dest_ln);

   /*!
    * @brief Execute schedule for filling ghosts on the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * @return refine schedule for filling ghost data from coarser level
    * and physical bc.
    *
    * @pre d_ghostfill_refine_schedules[dest_ln]
    */
   void
   xeqScheduleGhostFill(
      int dst_id,
      int dest_ln);

   /*!
    * @brief Execute schedule for filling ghosts on the specified
    * level or reregister an existing one.
    * This version does not get data from coarser levels.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * This function is used for the bottom solve level, since it does
    * not access data from any coarser level.  (Ghost data obtained
    * from coarser level must have been placed there before solve begins!)
    *
    * @return refine schedule for filling ghost data from same level
    * and physical bc.
    *
    * @pre d_ghostfill_nocoarse_refine_schedules[dest_ln]
    */
   void
   xeqScheduleGhostFillNoCoarse(
      int dst_id,
      int dest_ln);

   //@}

   //! @brief Return the patch data index for cell scratch data.
   int
   registerCellScratch() const;
   //! @brief Return the patch data index for flux scratch data.
   int
   registerFluxScratch() const;
   //! @brief Return the patch data index for outerflux scratch data.
   int
   registerOfluxScratch() const;

   //! @brief Free static variables at shutdown time.
   static void
   finalizeCallback();

   /*!
    * @brief Object dimension.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   //@{ @name Hierarchy-dependent objects.

   /*!
    * @brief Reference hierarchy
    *
    * This variable is non-null between the initializeOperatorState()
    * and deallocateOperatorState() calls.  It is not truly needed,
    * because the hierarchy is obtainable through variables in most
    * function argument lists.  We use it to enforce working on one
    * hierarchy at a time.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Coarsest level for solve.
    */
   int d_ln_min;

   /*!
    * @brief Finest level for solve.
    */
   int d_ln_max;

   /*!
    * @brief Description of coarse-fine boundaries.
    *
    * There is one coarse-fine boundary object for each level.
    * d_coarse_fine_boundary[i] is the description of
    * the coarse-fine boundary between level i and level i-1.
    * The coarse-fine boundary does not exist at the coarsest level,
    * although the hier::CoarseFineBoundary object still exists (it
    * should not contain any boxes).
    *
    * This array is initialized in initializeOperatorState() and
    * deallocated in deallocateOperatorState().  When allocated,
    * it is allocated for the index range [0,d_ln_max], though
    * the range [0,d_ln_min-1] is not used.  This is okay because
    * hier::CoarseFineBoundary is a light object before
    * it is set for a level.
    */
   std::vector<std::shared_ptr<hier::CoarseFineBoundary> > d_cf_boundary;

   //@}

   //@{
   /*!
    * @name Private state variables for solution process.
    */

   /*!
    * @brief Scalar Poisson equations specifications.
    * @see setPoissonSpecifications().
    */
   PoissonSpecifications d_poisson_spec;

   /*!
    * @brief Coarse level solver.
    */
   std::string d_coarse_solver_choice;

   /*!
    * @brief Coarse-fine discretization method.
    */
   std::string d_cf_discretization;

   /*!
    * @brief Coarse-fine discretization method.
    *
    * The name of the refinement operator used to prolong the
    * coarse grid correction.
    */
   std::string d_prolongation_method;

   /*!
    * @brief Tolerance specified to coarse solver
    */
   double d_coarse_solver_tolerance;

   /*!
    * @brief Coarse level solver iteration limit.
    */
   int d_coarse_solver_max_iterations;

   /*!
    * @brief Residual tolerance to govern smoothing.
    *
    * When we use one of the internal error smoothing functions
    * and want to terminate the smoothing sweeps at a certain
    * level of residual, this will be set to > 0.  If it is
    * < 0, the smoothing function effectively ignores it.
    *
    * This variable is needed because some coarse-level solver
    * simply runs the smoothing function until convergence.
    * It sets this variable to > 0, calls the smoothing function,
    * then resets it to < 0.
    */
   double d_residual_tolerance_during_smoothing;

   /*!
    * @brief Id of the flux.
    *
    * If set to -1, create and delete storage space on the fly.
    * Else, user has provided space for flux.
    *
    * @see setFluxId
    */
   int d_flux_id;

   /*!
    * Umpire Allocator for internal data allocations
    */
   tbox::ResourceAllocator d_allocator;

#ifdef HAVE_HYPRE
   /*!
    * @brief HYPRE coarse-level solver object.
    */
   std::shared_ptr<CellPoissonHypreSolver> d_hypre_solver;
#endif

   /*!
    * @brief Externally provided physical boundary condition object.
    *
    * see setPhysicalBcCoefObject()
    */
   const RobinBcCoefStrategy* d_physical_bc_coef;

   //@}

   //@{ @name Internal context and scratch data

   static std::shared_ptr<pdat::CellVariable<double> >
   s_cell_scratch_var[SAMRAI::MAX_DIM_VAL];

   static std::shared_ptr<pdat::SideVariable<double> >
   s_flux_scratch_var[SAMRAI::MAX_DIM_VAL];

   static std::shared_ptr<pdat::OutersideVariable<double> >
   s_oflux_scratch_var[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Default context of internally maintained hierarchy data.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   /*!
    * @brief ID of the solution-like scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::CellVariable<double> named
    * @c d_object_name+"::cell_scratch".
    * Scratch data is allocated and removed as needed
    * to reduce memory usage.
    */
   int d_cell_scratch_id;

   /*!
    * @brief ID of the side-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::SideVariable<double> named
    * @c d_object_name+"::flux_scratch".
    *
    * This data is allocated only as needed and deallocated
    * immediately after use.
    */
   int d_flux_scratch_id;

   /*!
    * @brief ID of the outerside-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::OutersideVariable<double> named
    * @c d_object_name+"::oflux_scratch".
    */
   int d_oflux_scratch_id;

   //@}

   //@{
   /*!
    * @name Various refine and coarsen objects used internally.
    */

   //! @brief Error prolongation (refinement) operator.
   std::shared_ptr<hier::RefineOperator> d_prolongation_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm> d_prolongation_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_prolongation_refine_schedules;

   //! @brief Solution restriction (coarsening) operator.
   std::shared_ptr<hier::CoarsenOperator> d_urestriction_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_urestriction_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
   d_urestriction_coarsen_schedules;

   //! @brief Residual restriction (coarsening) operator.
   std::shared_ptr<hier::CoarsenOperator> d_rrestriction_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_rrestriction_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
   d_rrestriction_coarsen_schedules;

   //! @brief Coarsen operator for outerflux-to-flux
   std::shared_ptr<hier::CoarsenOperator> d_flux_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_flux_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
   d_flux_coarsen_schedules;

   //! @brief Refine operator for cell-like data from coarser level.
   std::shared_ptr<hier::RefineOperator> d_ghostfill_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm> d_ghostfill_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_ghostfill_refine_schedules;

   //! @brief Refine operator for cell-like data from same level.
   std::shared_ptr<hier::RefineOperator>
   d_ghostfill_nocoarse_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm>
   d_ghostfill_nocoarse_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_ghostfill_nocoarse_refine_schedules;

   //@}

   /*!
    * @brief Utility object employed in setting ghost cells and providing
    * xfer::RefinePatchStrategy implementation.
    *
    * Since this class deals only in scalar variables having
    * Robin boundary conditions, we take advantage of the corresponding
    * implementation in CartesianRobinBcHelper.  Whenever
    * we need an implementation of xfer::RefinePatchStrategy,
    * this object is used.  Note that in the code, before we
    * use this object to set ghost cell values, directly or
    * indirectly by calling xfer::RefineSchedule::fillData(),
    * we must tell d_bc_helper the patch data index we want
    * to set and whether we are setting data with homogeneous
    * boundary condition.
    */
   CartesianRobinBcHelper d_bc_helper;

   //@{
   /*!
    * @name Non-essential objects used in outputs and debugging.
    */

   /*!
    * @brief Logging flag.
    */
   bool d_enable_logging;

   /*!
    * @brief Preconditioner using this object.
    *
    * This must remain a raw pointer.  Do not attempt to make this a smart
    * pointer.  The preconditioner and this class have data members referring
    * to each other.  If this is made into a smart pointer then things will not
    * tear themselves down properly leading to memory leaks unless there is
    * explicit action to reset the smart pointer.
    *
    * See setPreconditioner().
    */
   const FACPreconditioner* d_preconditioner;

   /*!
    * @brief Hierarchy cell operator used in debugging.
    */
   std::shared_ptr<math::HierarchyCellDataOpsReal<double> > d_hopscell;

   /*!
    * @brief Hierarchy side operator used in debugging.
    */
   std::shared_ptr<math::HierarchySideDataOpsReal<double> > d_hopsside;

   /*!
    * @brief Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_restrict_solution;
   std::shared_ptr<tbox::Timer> t_restrict_residual;
   std::shared_ptr<tbox::Timer> t_prolong;
   std::shared_ptr<tbox::Timer> t_smooth_error;
   std::shared_ptr<tbox::Timer> t_solve_coarsest;
   std::shared_ptr<tbox::Timer> t_compute_composite_residual;
   std::shared_ptr<tbox::Timer> t_compute_residual_norm;

   static tbox::StartupShutdownManager::Handler s_finalize_handler;
};

}
}

#endif // included_solv_CellPoissonFACOps
