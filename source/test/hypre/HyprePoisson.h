/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Example user class for solving Poisson using Hypre.
 *
 ************************************************************************/
#ifndef included_HyprePoisson
#define included_HyprePoisson

#include "SAMRAI/SAMRAI_config.h"

#if !defined(HAVE_HYPRE)

/*
 *************************************************************************
 * If the library is not compiled with hypre, print an error.
 * If we're running autotests, skip the error and compile an empty
 * class.
 *************************************************************************
 */
#if (TESTING != 1)
#error "This example requires SAMRAI be compiled with hypre."
#endif

#else

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/solv/CellPoissonHypreSolver.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include <memory>


namespace SAMRAI {

/*!
 * @brief Class to solve a sample Poisson equation on a SAMR grid.
 *
 * This class demonstrates how use the HYPRE Poisson solver
 * class to solve Poisson's equation on a single level
 * within a hierarchy.
 *
 * We set up and solve the following problem:
 *
 *   2d: div(grad(u)) = -2 (pi^2) sin(pi x) sin(pi y)
 *
 *   3d: div(grad(u)) = -3 (pi^2) sin(pi x) sin(pi y) sin(pi z)
 *
 * which has the exact solution
 *
 *   2d: u = sin(pi x) sin(pi y)
 *
 *   3d: u = sin(pi x) sin(pi y) sin(pi z)
 *
 * This class inherits and implements virtual functions from
 * - mesh::StandardTagAndInitStrategy to initialize data
 *   on the SAMR grid.
 * - appu::VisDerivedDataStrategy to write out certain data
 *   in a vis file, such as the error of the solution.
 */
class HyprePoisson:
   public mesh::StandardTagAndInitStrategy,
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    *
    * @param object_name Ojbect name
    * @param dim
    * @param database
    */
   HyprePoisson(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<solv::CellPoissonHypreSolver>& hypre_solver,
      std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs);

   virtual ~HyprePoisson();

   //@{ @name mesh::StandardTagAndInitStrategy virtuals

   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on
    * the grid.  All the information needed to initialize the grid
    * are in the arguments.
    *
    * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
    */
   virtual void
   initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

   /*!
    * @brief Reset any internal hierarchy-dependent information.
    */
   virtual void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      int coarsest_level,
      int finest_level);

   //@}

   //@{ @name appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id,
      double simulation_time = 0.0) const;

   //@}

   /*!
    * @brief Solve using HYPRE Poisson solver
    *
    * Set up the linear algebra problem and use a
    * solv::CellPoissonHypreSolver object to solve it.
    * -# Set initial guess
    * -# Set boundary conditions
    * -# Specify Poisson equation parameters
    * -# Call solver
    *
    * @return whether solver converged
    */
   bool
   solvePoisson();

   /*!
    * @brief Set up external plotter to plot internal
    * data from this class.
    *
    * After calling this function, the external
    * data writer may be used to write the
    * visit file for this object.
    *
    * The internal hierarchy is used and must be
    * established before calling this function.
    * (This is commonly done by building a hierarchy
    * with the mesh::StandardTagAndInitStrategy virtual
    * functions implemented by this class.)
    *
    * @param visit_writer VisIt data writer
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& visit_writer) const;

private:
   std::string d_object_name;

   const tbox::Dimension d_dim;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   //@{
   /*!
    * @name Major algorithm objects.
    */

   /*!
    * @brief HYPRE poisson solver.
    */
   std::shared_ptr<solv::CellPoissonHypreSolver> d_poisson_hypre;

   /*!
    * @brief Boundary condition coefficient implementation.
    */
   std::shared_ptr<solv::LocationIndexRobinBcCoefs> d_bc_coefs;

   //@}

   //@{
private:
   /*!
    * @name Private state variables for solution.
    */

   /*!
    * @brief Context owned by this object.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   /*!
    * @brief Descriptor indices of internal data.
    *
    * These are initialized in the constructor and never change.
    */
   int d_comp_soln_id, d_exact_id, d_rhs_id;

   //@}

};

}

#endif
#endif  // included_HyprePoisson
