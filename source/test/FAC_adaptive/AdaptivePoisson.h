/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AdaptivePoisson class declaration
 *
 ************************************************************************/
#ifndef included_AdaptivePoisson
#define included_AdaptivePoisson

#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/CellPoissonFACOps.h"
#include "SAMRAI/tbox/Dimension.h"
#include "PoissonSineSolution.h"
#include "PoissonPolynomialSolution.h"
#include "PoissonGaussianDiffcoefSolution.h"
#include "PoissonGaussianSolution.h"
#include "PoissonMultigaussianSolution.h"

#include <string>
#include <memory>

#include "SAMRAI/tbox/Database.h"

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/GhostCellRobinBcCoefs.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"


using namespace SAMRAI;

/*!
 * @brief Class to solve Poisson's equation on a SAMR grid.
 *
 * This class tests the FAC solver class solving
 * Poisson's equation on a SAMR grid.
 *
 * This class inherits and implements virtual functions from
 * mesh::StandardTagAndInitStrategy to initialize data
 * on the SAMR grid.
 */
class AdaptivePoisson:
   public mesh::StandardTagAndInitStrategy,
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    *
    * Requirements:
    * - the referenced objects
    *
    * Actions:
    * - Set up private member data
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    */
   AdaptivePoisson(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<solv::CellPoissonFACOps>& fac_ops,
      std::shared_ptr<solv::FACPreconditioner>& fac_precond,
      tbox::Database& database,
      /*! Log output stream */ std::ostream* log_stream = 0);

   //@{ @name mesh::StandardTagAndInitStrategy virtuals

public:
   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on the
    * grid.  Nevermind when it is called or where in the program that
    * happens.  All the information you need to initialize the grid
    * are in the arguments.
    *
    * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
    */
   virtual void
   initializeLevelData(
      /*! Hierarchy to initialize */
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Level to initialize */
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      /*! Whether level is being introduced for the first time */
      const bool initial_time,
      /*! Level to copy data from */
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true);

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ int coarsest_level,
      /*! Finest level */ int finest_level);

   virtual void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation);

   //@}

#if 0
   //@{
   /*!
    * @name Functions inherited from solv::CellPoissonFACOps
    */
   virtual void
   postprocessOneCycle(
      int iteration_num,
      const solv::SAMRAIVectorReal<double>& current_soln,
      const solv::SAMRAIVectorReal<double>& residual);
   //@}
#endif

   //@{ @name appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time) const;

   //@}

public:
   /*!
    * @brief Solve using FAC solver.
    *
    * Set up the linear algebra problem and use a
    * solv::FACPreconditioner object to solve it.
    *
    * @param hierarchy the hierarchy to solve on
    * @param initial_u how to set the initial guess for u.
    *       A std::string is used so the option "random" can be
    *       used.  If "random" is not used, set the std::string
    *       to a floating point number.
    */
   int
   solvePoisson(
      std::shared_ptr<hier::PatchHierarchy> hierarchy,
      std::string initial_u = std::string("0.0"));

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& visit_writer);
#endif

   /*!
    * @brief Compute the error of the current solution.
    *
    * Compute the @f$L_2@f$ and @f$L_\infty@f$ norms of the error,
    * for each level and over all levels.
    */
   int
   computeError(
      /*! hierarchy */ const hier::PatchHierarchy& hierarchy,
      /*! L2 norm */ double* l2norm,
      /*! L-inf norm */ double* linorm,
      /*! L2 norm on each level */ std::vector<double>& l2norms,
      /*! L-inf norm on each level */ std::vector<double>& linorms) const;

   /*!
    * @brief Compute error estimator (for adaption or plotting).
    *
    * Computes in the box defined by @c estimate_data.
    */
   void
   computeAdaptionEstimate(
      pdat::CellData<double>& estimate_data,
      const pdat::CellData<double>& soln_cell_data) const;

private:
   std::string d_name;
   const tbox::Dimension d_dim;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   //@{
   /*!
    * @name Major algorithm objects.
    */

   std::shared_ptr<solv::CellPoissonFACOps> d_fac_ops;

   std::shared_ptr<solv::FACPreconditioner> d_fac_preconditioner;

   //@}

   tbox::ResourceAllocator d_allocator;

   //@{

   /*!
    * @name Private state variables for solution.
    */

   /*!
    * @brief Context for persistent data.
    */
   std::shared_ptr<hier::VariableContext> d_context_persistent;

   /*!
    * @brief Context for scratch data.
    */
   std::shared_ptr<hier::VariableContext> d_context_scratch;

   /*!
    * @brief Diffusion coefficient.
    */
   std::shared_ptr<pdat::SideVariable<double> > d_diffcoef;

   /*!
    * @brief Flux.
    */
   std::shared_ptr<pdat::SideVariable<double> > d_flux;

   /*!
    * @brief Scalar solution of Poisson's equation.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_scalar;

   /*!
    * @brief Source for Poisson's equation.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_constant_source;

   /*!
    * @brief Linear source operator for linear system.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_ccoef;

   /*!
    * @brief Right hand side for linear system.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_rhs;

   /*!
    * @brief Exact solution.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_exact;

   /*!
    * @brief Residual.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_resid;

   /*!
    * @brief Vector weights.
    *
    * For cells not covered by a finer cell, the weight
    * is the volume.  For cells that are, the weight is zero.
    * This is used in computing norms on the AMR grid.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_weight;

   /*!
    * @brief Saved variable-context index.
    *
    * Because we refer to them often, variable-context indices are saved.
    * They are initialized in the constructor and never change.
    * Each index represents a variable-context pair in this class.
    * Thus the indices are @em not independent state variables.
    * If we had to set them, we need the associated variable and
    * context (and the variable database which manages the mapping).
    * See the hier::VariableDatabase class for more into.
    */
   int d_scalar_persistent, d_diffcoef_persistent,
       d_constant_source_persistent, d_weight_persistent,
       d_exact_persistent, d_rhs_scratch, d_resid_scratch,
       d_flux_persistent, d_ccoef_persistent;

   //@}

   //@{
private:
   /*!
    * @name Output streams.
    */
   /*!
    * @brief Log stream pointer.
    *
    * If set to NULL, no logging.
    */
   std::ostream* d_lstream;
   //@}

   //@{
   /*!
    * @name Miscellaneous.
    */
   std::string d_problem_name;
   //! @brief Poisson equation specifications.
   solv::PoissonSpecifications d_sps;
   //! @brief Things specific to the sinusoid solution
   PoissonSineSolution d_sine_solution;
   //! @brief Things specific to the Gaussian solution
   PoissonGaussianSolution d_gaussian_solution;
   //! @brief Things specific to the multi-Gaussian solution
   PoissonMultigaussianSolution d_multigaussian_solution;
   //! @brief Things specific to the polynomial solution
   PoissonPolynomialSolution d_polynomial_solution;
   //! @brief Things specific to the Gaussian coefficient solution
   PoissonGaussianDiffcoefSolution d_gaussian_diffcoef_solution;
   /*!
    * @brief Generic xfer::RefinePatchStrategy implementation for Robin bc.
    */
   solv::CartesianRobinBcHelper d_robin_refine_patch;
   /*!
    * @brief Physical bc coefficient strategy selecting one of the solutions'.
    */
   solv::RobinBcCoefStrategy* d_physical_bc_coef;
   //@}

   double d_adaption_threshold;

   //@{
private:
   /*!
    * @name Objects to help debugging.
    */
#ifdef HAVE_HDF5
   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif
   int d_finest_dbg_plot_ln;
   //@}

};

#endif  // included_AdaptivePoisson
