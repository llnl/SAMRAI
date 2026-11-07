/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Poisson solver
 *
 ************************************************************************/
#ifndef included_FACPoisson
#define included_FACPoisson

#include "SAMRAI/solv/CellPoissonFACSolver.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
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
 * This class demonstrates how use the FAC Poisson solver
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
 *
 * Inputs:  The only input parameter for this class is
 * "fac_poisson", the input database for the solv::CellPoissonFACSolver
 * object.  See the documentation for solv::CellPoissonFACSolver
 * for its input parameters.
 */
class FACPoisson:
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
    * @param fac_solver
    * @param bc_coefs
    */
   FACPoisson(
      const std::string& object_name,
      const tbox::Dimension& dim,
      const std::shared_ptr<solv::CellPoissonFACSolver>& fac_solver,
      const std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs);

   virtual ~FACPoisson();

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
      const std::shared_ptr<hier::PatchLevel>& old_level,
      const bool allocate_data);

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
      double simulation_time) const;

   //@}

   /*!
    * @brief Solve using HYPRE Poisson solver
    *
    * Set up the linear algebra problem and use a
    * solv::CellPoissonFACSolver object to solve it.
    * -# Set initial guess
    * -# Set boundary conditions
    * -# Specify Poisson equation parameters
    * -# Call solver
    */
   int
   solvePoisson();

#ifdef HAVE_HDF5
   /*!
    * @brief Set up external plotter to plot internal
    * data from this class.
    *
    * After calling this function, the external
    * data writer may be used to write the
    * viz file for this object.
    *
    * The internal hierarchy is used and must be
    * established before calling this function.
    * (This is commonly done by building a hierarchy
    * with the mesh::StandardTagAndInitStrategy virtual
    * functions implemented by this class.)
    *
    * @param viz_writer VisIt writer
    */
   int
   setupPlotter(
      appu::VisItDataWriter& plotter) const;
#endif

private:
   std::string d_object_name;

   const tbox::Dimension d_dim;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   //@{
   /*!
    * @name Major algorithm objects.
    */

   /*!
    * @brief FAC poisson solver.
    */
   std::shared_ptr<solv::CellPoissonFACSolver> d_poisson_fac_solver;

   /*!
    * @brief Boundary condition coefficient implementation.
    */
   std::shared_ptr<solv::LocationIndexRobinBcCoefs> d_bc_coefs;

   //@}

   //@{

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

#endif  // included_FACPoisson
