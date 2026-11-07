/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic method-of-lines time integration algorithm
 *
 ************************************************************************/

#ifndef included_algs_MethodOfLinesIntegrator
#define included_algs_MethodOfLinesIntegrator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/MethodOfLinesPatchStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace algs {

/*!
 * \brief Class MethodOfLinesIntegrator implements a spatially
 * adaptive version of the Strong Stability Preserving (SSP) Runge-Kutta
 * time integration algorithm.
 *
 * The original non-adaptive version of the algorithm is described in
 * S. Gottlieb, C.W. Shu, E. Tadmor, SIAM Review, Vol. 43, No. 1, pp. 89-112.
 * The advanceHierarchy() method integrates all levels of an AMR hierarchy
 * through a specified timestep.  See this method for details of the
 * time-stepping process.  Application-specific numerical routines that
 * are necessary for these operations are provided by the
 * MethodOfLinesPatchStrategy data member.  The collaboration between
 * this class and the patch strategy follows the the Strategy design pattern.
 * A concrete patch strategy object is derived from the base class to
 * provide those routines for a specific problem.
 *
 * This class is derived from the mesh::StandardTagAndInitStrategy abstract
 * base class which defines an interface for routines required by the
 * dynamic adaptive mesh refinement routines in the mesh::GriddingAlgorithm
 * class.  This collaboration also follows the Strategy design pattern.
 *
 * Initialization of an MethodOfLinesIntegrator object is performed
 * by first setting default values, then reading from input.  All input
 * values may override values read from restart.  Data read from input is
 * summarized as follows:
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b    alpha_1
 *    - \b    alpha_2
 *    - \b    beta <br>
 *       arrays of double values (length = order) specifying the coeffients
 *       used in the multi-step Strong Stability Preserving (SSP) Runge-Kutta
 *       algorithm.
 *
 * Note that when continuing from restart, the input parameters in the input
 * database override all values read in from the restart database.
 *
 * <b> Details: </b> <br>
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
 *      <td>alpha_1</td>
 *      <td>array of doubles</td>
 *      <td>[1.0, 0.75, 2.0/3.0]</td>
 *      <td>any doubles but no more than 3 of them</td>
 *      <td>opt</td>
 *      <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *      <td>alpha_2</td>
 *      <td>array of doubles</td>
 *      <td>[0.0, 0.25, 2.0/3.0]</td>
 *      <td>any doubles but no more than 3 of them</td>
 *      <td>opt</td>
 *      <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *      <td>beta</td>
 *      <td>array of doubles</td>
 *      <td>[1.0, 0.25, 2.0/3.0]</td>
 *      <td>any doubles but no more than 3 of them</td>
 *      <td>opt</td>
 *      <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * The following represents a sample input entry:
 *
 * @code
 *  MethodOfLinesIntegrator{
 *     alpha_1               = 1., 0.75, 0.33333
 *     alpha_2               = 0., 0.25, 0.66666
 *     beta                  = 1., 0.25, 0.66666
 *  }
 *  @endcode
 *
 * @see mesh::StandardTagAndInitStrategy
 */

class MethodOfLinesIntegrator:
   public tbox::Serializable,
   public mesh::StandardTagAndInitStrategy
{
public:
   /*!
    * Enumerated type for the different categories of variable
    * quantities allowed by the method of lines integration algorithm.
    * See registerVariable(...) function for more details.
    *
    * - \b SOLN           {Solution quantity for time-dependent ODE problem
    *                      solved by RK time-stepping algorithm.}
    * - \b RHS            {Right-hand-side of ODE problem solved;
    *                      i.e., du/dt = RHS.}
    *
    *
    *
    */
   enum MOL_VAR_TYPE { SOLN = 0,
                       RHS = 1 };

   /*!
    * The constructor for MethodOfLinesIntegrator configures the method
    * of lines integration algorithm with the concrete patch strategy object
    * (containing problem-specific numerical routines) and initializes
    * integration algorithm parameters provided in the specified input
    * database and in the restart database corresponding to the
    * specified object_name.
    *
    * @pre !object_name.empty()
    * @pre patch_strategy != 0
    */
   MethodOfLinesIntegrator(
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db,
      MethodOfLinesPatchStrategy* patch_strategy);

   /*!
    * The destructor for MethodOfLinesIntegrator unregisters
    * the integrator object with the restart manager.
    */
   virtual ~MethodOfLinesIntegrator();

   /*!
    * Initialize integrator by setting the number of time levels
    * of data needed based on specifications of the gridding algorithm.
    *
    * This routine also invokes variable registration in the patch strategy.
    *
    * @pre gridding_alg
    */
   void
   initializeIntegrator(
      const std::shared_ptr<mesh::GriddingAlgorithm>& gridding_alg);

   /*!
    * Return a suitable time increment over which to integrate the ODE
    * problem.  A minimum is taken over the increment computed on
    * each patch in the hierarchy.
    *
    * @pre hierarchy
    */
   double
   getTimestep(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const double time) const;

   /*!
    * Advance the solution through the specified dt, which is assumed
    * for the problem and state of the solution.  Advances all patches
    * in the hierarchy passed in.
    *
    * @pre hierarchy
    */
   void
   advanceHierarchy(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const double time,
      const double dt);

   /*!
    * Register variable quantity defined in the patch strategy with the
    * method of lines integrator which manipulates its storage.
    *
    * @pre variable
    * @pre transfer_geom
    * @pre variable->getDim() == ghosts.getDim()
    */
   void
   registerVariable(
      const std::shared_ptr<hier::Variable>& variable,
      const hier::IntVector& ghosts,
      const MOL_VAR_TYPE m_v_type,
      const std::shared_ptr<hier::BaseGridGeometry>& transfer_geom,
      const std::string& coarsen_name = std::string(),
      const std::string& refine_name = std::string());

   /*!
    * Print all data members of MethodOfLinesIntegrator object.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*!
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.  The old_level pointer corresponds to
    * the level that resided in the hierarchy before the level with the
    * specified number was introduced.  If the pointer is NULL, there was
    * no level in the hierarchy prior to the call and the level data is set
    * based on the user routines and the simulation time.  Otherwise, the
    * specified level replaces the old level and the new level receives data
    * from the old level appropriately before it is destroyed.
    *
    * Typically, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the data is to be set, the level number must
    * match that of the old level, if non-NULL.  If the old level is
    * non-NULL, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.  Then, user-supplied patch routines are called to further
    * initialize the data if needed.  The boolean argument after_regrid
    * is passed into the user's routines.
    *
    * The boolean argument initial_time indicates whether the integration
    * time corresponds to the initial simulation time.  If true, the level
    * should be initialized with initial simulation values.  Otherwise, it
    * should be assumed that the simulation time is at some point after the
    * start of the simulation.  This information is provided since the
    * initialization of the data on a patch may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    *
    * @pre hierarchy
    * @pre level_number >= 0
    * @pre hierarchy->getPatchLevel(level_number)
    * @pre !old_level || (level_number == old_level->getLevelNumber())
    * @pre !old_level || (hierarchy->getDim() == old_level->getDim())
    */
   void
   initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

#if !defined(__xlC__)
   using mesh::StandardTagAndInitStrategy::initializeLevelData;
#endif

   /*!
    * Reset cached communication schedules after the hierarchy has changed
    * (due to regridding, for example) and the data has been initialized on
    * the new levels.  The intent is that the cost of data movement on the
    * hierarchy will be amortized across multiple communication cycles,
    * if possible.  Note, that whenever this routine is called, communication
    * schedules are updated for every level finer than and including that
    * indexed by coarsest_level.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (coarsest_level <= finest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    */
   void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    * Set integer tags to "one" on the given level where refinement
    * of that level should occur using the user-supplied gradient detector.
    * The boolean argument initial_time is true when the level is being
    * subject to error estimation at initialization time.  If it is false,
    * the error estimation process is being invoked at some later time
    * after the AMR hierarchy was initially constructed.  The boolean argument
    * uses_richardson_extrapolation_too is true when Richardson
    * extrapolation error estimation is used in addition to the gradient
    * detector, and false otherwise.  This argument helps the user to
    * manage multiple regridding criteria.  This information
    * is passed along to the user's patch data tagging routines since the
    * application of the error estimator may be different in each of those
    * circumstances.
    *
    * @pre hierarchy
    * @pre hierarchy->getPatchLevel(level_number)
    */
   virtual void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /*!
    * Writes object state out to the given restart database.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int ALGS_METHOD_OF_LINES_INTEGRATOR_VERSION;

   /*
    * Copy all solution data from current context to scratch context.
    */
   void
   copyCurrentToScratch(
      const std::shared_ptr<hier::PatchLevel>& level) const;

   /*
    * Copy all solution data from scratch context to current context.
    */
   void
   copyScratchToCurrent(
      const std::shared_ptr<hier::PatchLevel>& level) const;

   /*
    * Reads in parameters from the input database.  All
    * values from the input file take precedence over values from the
    * restart file.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    *
    *    -Data is missing from the restart database or is inconsistent.
    *
    */
   void
   getFromRestart();

   /*
    * The object name is used as a handle to the database stored in
    * restart files and for error reporting purposes.
    */
   std::string d_object_name;

   /*
    * Order of the Runge-Kutta method, and array of alpha values used in
    * updating solution during multi-step process.
    */
   int d_order;
   std::vector<double> d_alpha_1;
   std::vector<double> d_alpha_2;
   std::vector<double> d_beta;

   /*
    * A pointer to the method of lines patch model that will perform
    * the patch-based numerical operations.
    */
   MethodOfLinesPatchStrategy* d_patch_strategy;

   /*
    * The communication algorithms and schedules are created and
    * maintained to manage inter-patch communication during AMR integration.
    * The algorithms are created in the class constructor.  They are
    * initialized when variables are registered with the integrator.
    */

   /*
    * The "advance" schedule is used prior to advancing a level and
    * prior to computing dt at initialization.  It must be reset each
    * time a level is regridded.  All ghosts are filled with current
    * data at specified time.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_advance;
   std::vector<std::shared_ptr<xfer::RefineSchedule> > d_bdry_sched_advance;

   /*
    * Algorithm for transferring data from coarse patch to fine patch
    * after a regrid.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_fill_after_regrid;

   /*
    * Algorithm for copying data from current context to scratch context,
    * on the same level.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_fill_before_tagging;

   /*
    * Algorithm and communication schedule for transferring data from
    * fine to coarse grid.
    */
   std::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> > d_coarsen_schedule;

   /*
    * This algorithm has two variable contexts.  The current context is the
    * solution state at the current simulation time.  These data values do
    * not need ghost cells.  The scratch context is the temporary solution
    * state during the time integration process.  These variables will require
    * ghost cell widths that depend on the spatial discretization.
    */
   std::shared_ptr<hier::VariableContext> d_current;
   std::shared_ptr<hier::VariableContext> d_scratch;

   std::list<std::shared_ptr<hier::Variable> > d_soln_variables;
   std::list<std::shared_ptr<hier::Variable> > d_rhs_variables;

   /*
    * The component selectors for current and scratch are used by the
    * algorithm to collectively allocate/deallocate data for variables.
    */
   hier::ComponentSelector d_current_data;
   hier::ComponentSelector d_scratch_data;
   hier::ComponentSelector d_rhs_data;

};

}
}

#endif
