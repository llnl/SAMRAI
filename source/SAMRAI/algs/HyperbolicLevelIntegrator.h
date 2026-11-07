/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Integration routines for single level in AMR hierarchy
 *                (basic hyperbolic systems)
 *
 ************************************************************************/

#ifndef included_algs_HyperbolicLevelIntegrator
#define included_algs_HyperbolicLevelIntegrator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/tbox/Timer.h"

#include <list>
#include <vector>
#include <memory>

#define HLI_RECORD_STATS
// #undef DGA_RECORD_STATS

#ifdef HLI_RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

#include <iostream>
#include <string>

namespace SAMRAI {
namespace algs {

/**
 * Class HyperbolicLevelIntegrator provides routines needed to
 * integrate a system of hyperbolic conservation laws on a structured
 * AMR patch hierarchy using local time refinement.  The routines include
 * initializing a level, advance a level, and synchronize levels in a
 * time-dependent AMR application.  The AMR timestepping algorithm that
 * cycles through the patch levels and calls these routines is provided by
 * the TimeRefinementIntegrator class.  Together, that hierarchy
 * integration class and this single level integration class produce the
 * common AMR algorithm due to Berger, Colella and Oliger
 * (see e.g., Berger and Colella, J. Comp. Phys. (82)1:64-84, 1989).
 * The operations performed on single patches on each level are implemented
 * in the user-defined, problem-specific class derived from the abstract
 * base class HyperbolicPatchStrategy.
 *
 * It is important to note that the variable contexts used by the concrete
 * patch strategy subclass must be consistent with those defined in this
 * class which manages the data for the variables.
 *
 * This class is derived from the abstract base class
 * TimeRefinementLevelStrategy, which defines routines needed by
 * the time refinement integrator.  There is an argument in the constructor
 * that determines whether this class will be used by the time
 * refinement integrator for refined timestepping or synchronized
 * timestepping.  The routines overloaded in
 * TimeRefinementLevelStrategy are: initializeLevelIntegrator(),
 * getLevelDt(), getMaxFinerLevelDt(), advanceLevel(),
 * standardLevelSynchronization(), synchronizeNewLevels(),
 * resetTimeDependentData(), and resetDataToPreadvanceState().
 * This class is also derived from mesh::StandardTagAndInitStrategy,
 * which defines routines needed by the gridding algorithm classes.  The
 * routines overloaded in mesh::StandardTagAndInitStrategy are:
 * initializeLevelData(), resetHierarchyConfiguration(),
 * applyGradientDetector(), applyRichardsonExtrapolation(), and
 * coarsenDataForRichardsonExtrapolation().
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b    cfl
 *       the CFL factor used for timestep selection (dt used = CFL * max dt)
 *
 *    - \b    cfl_init
 *       the CFL factor used for initial timestep
 *
 *    - \b    lag_dt_computation
 *       indicates whether dt is based on current solution or solution from
 *       previous step (possible optimization in communication for
 *       characteristic analysis)
 *
 *    - \b    use_ghosts_to_compute_dt
 *       indicates whether ghost data must be filled before timestep is
 *       computed on each patch (possible communication optimization)
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
 *     <td>cfl</td>
 *     <td>double</td>
 *     <td>none</td>
 *     <td>any double</td>
 *     <td>req</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>cfl_init</td>
 *     <td>double</td>
 *     <td>none</td>
 *     <td>any double</td>
 *     <td>req</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>lag_dt_computation</td>
 *     <td>bool</td>
 *     <td>TRUE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>use_ghosts_to_compute_dt</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * A sample input file entry might look like:
 *
 * @code
 *    cfl = 0.9
 *    cfl_init = 0.9
 *    lag_dt_computation = FALSE
 *    use_ghosts_to_compute_dt = TRUE
 * @endcode
 *
 * @see TimeRefinementIntegrator
 * @see mesh::StandardTagAndInitStrategy
 * @see HyperbolicPatchStrategy
 */

class HyperbolicLevelIntegrator:
   public TimeRefinementLevelStrategy,
   public mesh::StandardTagAndInitStrategy,
   public tbox::Serializable
{
public:
   /**
    * Enumerated type for the different ways in which variable storage
    * can be manipulated by the level integration algorithm.
    * See registerVariable(...) function for more details.
    *
    *
    *
    * - \b TIME_DEP      {Data that changes in time and needs more than one
    *                      time level to be stored.}
    * - \b INPUT          {Data that is set once and do not change during
    *                      the ghosts are never re-filled outside of
    *                      user-defined routines.}
    * - \b FLUX           {Face-centered double values used in conservative
    *                      difference and synchronization (i.e., refluxing)
    *                      process.  A corresponding variable to store flux
    *                      integral information is created for each FLUX
    *                      variable.}
    * - \b TEMPORARY      {Accessory values intended to live only for
    *                      computation on a single patch (i.e., they cannot
    *                      be assumed to exist between patch routine function
    *                      calls.)}
    *
    *
    *
    */
   enum HYP_VAR_TYPE { TIME_DEP = 0,
                       INPUT = 1,
                       NO_FILL = 2,
                       FLUX = 3,
                       TEMPORARY = 4 };

   /**
    * Constructor for HyperbolicLevelIntegrator initializes
    * integration parameters to default values and constructs standard
    * communication algorithms.  Other data members are read in from
    * the specified input database or the restart database corresponding
    * to the specified object_name.  This class is used by
    * the time refinement integrator for refined timestepping when the
    * use_time_refinement argument is true, and for synchronized
    * timestepping when the boolean is false.
    *
    * @pre !object_name.empty()
    * @pre patch_strategy != 0
    */
   HyperbolicLevelIntegrator(
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db,
      HyperbolicPatchStrategy* patch_strategy,
      bool use_time_refinement = true);

   /**
    * The destructor for HyperbolicLevelIntegrator unregisters
    * the integrator object with the restart manager.
    */
   virtual ~HyperbolicLevelIntegrator();

   /**
    * Initialize level integrator by by setting the number of time levels
    * of data needed based on specifications of the gridding algorithm.
    *
    * This routine also invokes variable registration in the patch strategy.
    *
    * @pre gridding_alg_strategy is actually a std::dynamic_pointer_cast<mesh::GriddingAlgorithm>
    */
   virtual void
   initializeLevelIntegrator(
      const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg_strategy);

   /**
    * Determine time increment to advance data on level and return that
    * value.  The double dt_time argument is the simulation time when
    * the routine is called.  The initial_time boolean is true if this
    * routine is called during hierarchy initialization (i.e., at the
    * initial simulation time).  Otherwise, it is false.  The
    * recompute_dt option specifies whether to compute the timestep using
    * the current level data or to return the value stored by the time
    * integrator. The default true setting means the timestep will be
    * computed if no value is supplied.
    *
    * @pre level
    */
   virtual double
   getLevelDt(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double dt_time,
      const bool initial_time);

   /**
    * Return the maximum allowable time increment for the level with
    * the specified level number based on the time increment for the
    * next coarser level and the mesh refinement ratio between the two
    * levels.  For the common explicit integration methods for hyperbolic
    * conservation laws (constrained by a CFL limit), the fine time increment
    * is typically the coarse increment divided by the refinement ratio.
    *
    * @pre each component of ratio_to_coarser is > 0
    */
   virtual double
   getMaxFinerLevelDt(
      const int finer_level_number,
      const double coarse_dt,
      const hier::IntVector& ratio_to_coarser);

   /**
    * Integrate data on all patches on the given patch level from current
    * time (current_time) to new time (new_time).  This routine is used
    * to advance the solution on each level in the hierarchy and during
    * time-dependent regridding procedures, such as Richardson extrapolation.
    * The boolean arguments are used to determine the state of the algorithm
    * and the data when the advance routine is called.   The first_step
    * and last_step indicate whether the step is the first or last in the
    * current timestep sequence on the level.  Typically, the current timestep
    * sequence means each step on the level between advance steps on a
    * coarser level in the hierarchy, if one exists.  The regrid_advance
    * value is true when the advance is called as part of a time-dependent
    * regridding procedure.  Usually when this happens, the results of the
    * colution advance will be discarded.  So, for example, when this is true
    * flux information is not maintained and flux integrals are not updated.
    * The final boolean argument indicates whether or not the level resides
    * in the hierarchy.  For example, during time-dependent regridding, such
    * as Richardson extrapolation, a temporary level that is not in the
    * hierarchy is created and advanced.  Then, a communication schedule
    * must be generated for the level before the advance begins.
    *
    * This routine is called at two different points during time integration:
    * during the regular time advance sequence, and possibly at the initial
    * simulation time.  The second call advances the solution on a coarser
    * level ahead in time to provide time-dependent boundary values for some
    * finer level when time-dependent regridding is used.  In the first case,
    * the values of the boolean flags are:
    *
    *
    *
    *    - \b  first_step
    *        = true for first step in level time step sequence; else, false.
    *    - \b  last_step
    *        = true for last step in level time step sequence; else, false.
    *    - \b  regrid_advance
    *        = false.
    *
    *
    *
    * In the second case, the values of the boolean flags are:
    *
    *
    *
    *    - \b  first_step
    *        = true.
    *    - \b  last_step
    *        = false.
    *    - \b  regrid_advance
    *        = true.
    *
    *
    *
    *
    * When time-dependent regridding (i.e., Richardson extrapolation) is
    * used, the routine is called from two different points in addition to
    * those described above: to advance a temporary level that is coarser
    * than the hierarchy level on which error estimation is performed, and
    * to advance the hierarchy level itself.  In the first case, the values of
    * the boolean flags are:
    *
    *
    *
    *    - \b  first_step
    *        = true.
    *    - \b  last_step
    *        = true.
    *    - \b  regrid_advance
    *        = true.
    *
    *
    *
    * In the second case, the values of the boolean flags are:
    *
    *
    *
    *    - \b  first_step
    *      (when regridding during time integration sequence)
    *        = true when the level is not coarsest level to synchronize
    *          immediately before the regridding process; else, false.
    *      (when generating initial hierarchy construction)
    *        = true, even though there may be multiple advance steps.
    *    - \b  last_step
    *        = true when the advance is the last in the Richardson
    *          extrapolation step sequence; else false.
    *    - \b  regrid_advance
    *        = true.
    *
    * @pre level
    * @pre hierarchy
    * @pre current_time <= new_time
    * @pre level->getDim() == hierarchy->getDim()
    */

   virtual double
   advanceLevel(
      const std::shared_ptr<hier::PatchLevel>& level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const double current_time,
      const double new_time,
      const bool first_step,
      const bool last_step,
      const bool regrid_advance = false);

   /**
    * Synchronize data between given patch levels in patch hierarchy
    * according to the standard hyperbolic AMR flux correction algorithm.
    * This routine synchronizes data between two levels at a time from
    * the level with index finest_level down to the level with index
    * coarsest_level.  The array of old time values are used in the
    * re-integration of the time-dependent data.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (coarsest_level < finest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    * @pre old_times.size() >= finest_level
    * @pre for each i in [coarsest_level, finest_level) hierarchy->getPatchLevel(i) && sync_time >= old_times[i]
    * @pre hierarchy->getPatchLevel(finest_level)
    */
   virtual void
   standardLevelSynchronization(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const std::vector<double>& old_times);

   /**
    * This overloaded version of standardLevelSynchronization implements
    * a routine used for synchronized timestepping.  Only a single
    * value for the old time is needed, since all levels would have the
    * same old time.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (coarsest_level < finest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    * @pre for each i in [coarsest_level, finest_level) hierarchy->getPatchLevel(i)
    * @pre sync_time >= old_time
    * @pre hierarchy->getPatchLevel(finest_level)
    */
   virtual void
   standardLevelSynchronization(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const double old_time);

   /**
    * Coarsen current solution data from finest hierarchy level specified
    * down through the coarsest hierarchy level specified, if initial_time
    * is true.  In this case, the hierarchy is being constructed at the
    * initial simulation time,  After data is coarsened, the application-
    * specific initialization routine is called to set data before that
    * solution is further coarsened to the next coarser level in the
    * hierarchy.  This operation makes the solution consistent between
    * coarser levels and finer levels that did not exist when the coarse
    * levels where created and initialized originally.
    *
    * When initial_time is false, this routine does nothing since the
    * standard hyperbolic AMR algorithm for conservation laws requires
    * no data synchronization after regridding beyond interpolation of
    * data from coarser levels in the hierarchy in some conservative fashion.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (coarsest_level < finest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    * @pre for each i in [coarsest_level, finest_level] hierarchy->getPatchLevel(i)
    */
   virtual void
   synchronizeNewLevels(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const bool initial_time);

   /**
    * Resets time-dependent data storage and update time for patch level.
    *
    * @pre level
    */
   virtual void
   resetTimeDependentData(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double new_time,
      const bool can_be_refined);

   /**
    * Deallocate all new simulation data on the given level.  This may
    * be necessary during regridding, or setting up levels initially.
    *
    * @pre level
    */
   virtual void
   resetDataToPreadvanceState(
      const std::shared_ptr<hier::PatchLevel>& level);

   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.  The old_level pointer corresponds to
    * the level that resided in the hierarchy before the level with the
    * specified number was introduced.  If the pointer is null, there was
    * no level in the hierarchy prior to the call and the level data is set
    * based on the user routines and the simulation time.  Otherwise, the
    * specified level replaces the old level and the new level receives data
    * from the old level appropriately before it is destroyed.
    *
    * Typically, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the data is to be set, the level number must
    * match that of the old level, if non-null.  If the old level is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.  Then, user-supplied patch routines are called to further
    * initialize the data if needed.  The boolean argument initial_time
    * is passed into the user's routines.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data on a patch may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest level allowed in the hierarchy.  This may or
    * may not affect the data initialization process depending on the problem.
    *
    * @pre hierarchy
    * @pre (level_number >= 0) &&
    *      (level_number <= hierarchy->getFinestLevelNumber())
    * @pre !old_level || (level_number == old_level->getLevelNumber())
    * @pre hierarchy->getPatchLevel(level_number))
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

#if !defined(__xlC__)
   using mesh::StandardTagAndInitStrategy::initializeLevelData;
#endif

   /**
    * Reset cached communication schedules after the hierarchy has changed
    * (due to regidding, for example) and the data has been initialized on
    * the new levels.  The intent is that the cost of data movement on the
    * hierarchy will be amortized across multiple communication cycles,
    * if possible.  The level numbers indicate the range of levels in the
    * hierarchy that have changed.  However, this routine updates
    * communication schedules every level finer than and including that
    * indexed by the coarsest level number given.
    *
    * @pre hierarchy
    * @pre (coarsest_level >= 0) && (coarsest_level <= finest_level) &&
    *      (finest_level <= hierarchy->getFinestLevelNumber())
    * @pre for all i <= finest_level hierarchy->getPatchLevel(i) is non-NULL
    */
   virtual void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Set integer tags to "one" in cells where refinement of the given
    * level should occur according to some user-supplied gradient criteria.
    * The double time argument is the regrid time.  The integer "tag_index"
    * argument is the patch descriptor index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument
    * initial_time indicates whether the level is being subject to refinement
    * at the initial simulation time.  If it is false, then the error
    * estimation process is being invoked at some later time after the AMR
    * hierarchy was initially constructed.  The boolean argument
    * uses_richardson_extrapolation_too is true when Richardson
    * extrapolation error estimation is used in addition to the gradient
    * detector, and false otherwise.  This argument helps the user to
    * manage multiple regridding criteria.   This information is passed along
    * to the user's patch tagging routines since the application of the
    * gradient detector may be different in each case.
    *
    * @pre hierarchy
    * @pre (level_number >= 0) &&
    *      (level_number <= hierarchy->getFinestLevelNumber())
    * @pre hierarchy->getPatchLevel(level_number)
    */
   virtual void
   applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   /**
    * Set integer tags to "one" where refinement onf the given
    * level should occur according to some user-supplied Richardson
    * extrapolation criteria.  The "error_data_time" argument is the
    * regrid time.  The "deltat" argument is the time increment to advance
    * the solution on the level to be refined.  Note that that level is
    * finer than the level in the argument list, in general.  The
    * ratio between the argument level and the actual hierarchy level
    * is given by the integer "coarsen ratio".
    *
    * The integer "tag_index" argument is the patch descriptor index of
    * the cell-centered integer tag array on each patch in the hierarchy.
    *
    * The boolean argument initial_time indicates whether the level is being
    * subject to refinement at the initial simulation time.  If it is false,
    * then the error estimation process is being invoked at some later time
    * after the AMR hierarchy was initially constructed.  Typically, this
    * information is passed to the user's patch tagging routines since the
    * application of the Richardson extrapolation process may be different
    * in each case.
    *
    * The boolean uses_gradient_detector_too is true when a gradient
    * detector procedure is used in addition to Richardson extrapolation,
    * and false otherwise.  This argument helps the user to manage multiple
    * regridding criteria.
    *
    * @pre level
    */
   virtual void
   applyRichardsonExtrapolation(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double error_data_time,
      const int tag_index,
      const double deltat,
      const int error_coarsen_ratio,
      const bool initial_time,
      const bool uses_gradient_detector_too);

   /**
    * Coarsen solution data from level to coarse_level for Richardson
    * extrapolation.  Note that this routine will be called twice during
    * the Richardson extrapolation error estimation process.  The
    * before_advance boolean argument indicates whether data is
    * set on the coarse level by coarsening the "old" time level solution
    * (i.e., before it has been advanced) or by coarsening the "new"
    * solution on the fine level (i.e., after it has been advanced).
    *
    * @pre hierarchy
    * @pre (level_number >= 0) &&
    *      (level_number <= hierarchy->getFinestLevelNumber())
    * @pre (hierarchy->getPatchLevel(level_number)
    * @pre coarse_level
    * @pre hierarchy->getDim() == coarse_level->getDim()
    */
   virtual void
   coarsenDataForRichardsonExtrapolation(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const std::shared_ptr<hier::PatchLevel>& coarse_level,
      const double coarsen_data_time,
      const bool before_advance);

   /**
    * Register a variable with the hyperbolic integration algorithm.  The
    * variable type must be one of the options defined by the enumerated
    * type defined above.  Typically, this routine is called from the
    * hyperbolic patch model when the variable registration process is
    * invoked by calling the function initializeLevelIntegrator() above.
    * In fact, that function should be called before this routine is called.
    *
    * @pre var
    * @pre transfer_geom
    * @pre ghosts.getDim() == var->getDim()
    */
   virtual void
   registerVariable(
      const std::shared_ptr<hier::Variable>& var,
      const hier::IntVector ghosts,
      const HYP_VAR_TYPE h_v_type,
      const std::shared_ptr<hier::BaseGridGeometry>& transfer_geom,
      const std::string& coarsen_name = std::string(),
      const std::string& refine_name = std::string());

   /**
    * Print class data representation for hyperbolic level integrator object.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Write out object state to the given restart database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Return pointer to "current" variable context used by integrator.
    * Current data corresponds to state data at the beginning of a
    * timestep, or when a new level is initialized.
    */
   std::shared_ptr<hier::VariableContext>
   getCurrentContext() const
   {
      return d_current;
   }

   /**
    * Return pointer to "new" variable context used by integrator.
    * New data corresponds to advanced state data at the end of a timestep.
    * The data is one timestep later than the "current" data.
    */
   std::shared_ptr<hier::VariableContext>
   getNewContext() const
   {
      return d_new;
   }

   /**
    * Return pointer to "old" variable context used by integrator.
    * Old data corresponds to an extra time level of state data used
    * for Richardson extrapolation error estimation.  The data is
    * one timestep earlier than the "current" data.
    *
    * Note that only in certain cases when using time-dependent error
    * estimation, such as Richardson extrapolation, is the returned
    * pointer will non-null.  See contructor for more information.
    */
   std::shared_ptr<hier::VariableContext>
   getOldContext() const
   {
      return d_old;
   }

   /**
    * Return pointer to "scratch" variable context used by integrator.
    * Scratch data typically corresponds to storage that user-routines
    * in the concrete HyperbolicPatchStrategy object manipulate;
    * in particular, scratch data contains ghost cells.
    */
   std::shared_ptr<hier::VariableContext>
   getScratchContext() const
   {
      return d_scratch;
   }

   /**
    * Return pointer to variable context used for plotting.  This
    * context corresponds to the data storage that should be written
    * to plot files.  Typically, this is the same as the "current" context.
    */
   std::shared_ptr<hier::VariableContext>
   getPlotContext() const
   {
      return d_plot_context;
   }

   /**
    * Return true if this class has been constructed to use refined
    * timestepping and false if it has been constructed to use
    * synchronized timestepping.
    */
   bool
   usingRefinedTimestepping() const
   {
      return d_use_time_refinement;
   }

   /*
    * Write out statistics recorded on numbers of cells and patches generated.
    */
   void
   printStatistics(
      std::ostream& s = tbox::plog) const;

   /**
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

protected:
   /**
    * Read values, indicated above, from given input database.  The boolean
    * argument is_from_restart should be set to true if the simulation
    * is beginning from restart.  Otherwise it should be set to false.
    *
    * @pre is_from_restart || input_db
    */
   virtual void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /**
    * Read object state from the restart file and initialize class
    * data members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *
    *
    *
    *    -
    *        The database corresponding to object_name is not found
    *        in the restart file.
    *
    *    -
    *        The class version number and restart version number do not
    *        match.
    *
    *
    *
    *
    *
    */
   virtual void
   getFromRestart();

   /*
    * Pre-process flux storage before advancing solution on level from
    * cur_time to new_time.  The boolean flags are used to determine
    * how flux and flux integral storage is allocated and initialized.
    * These are needed since the advanceLevel() routine is used for
    * both level integration and time-dependent error estimation.
    */
   virtual void
   preprocessFluxData(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double cur_time,
      const double new_time,
      const bool regrid_advance,
      const bool first_step,
      const bool last_step);

   /*
    * Post-process flux storage after advancing solution on level.
    * The boolean flag is used to determine how flux and flux integral
    * storage is copied and de-allocated. This is needed since the
    * advanceLevel() routine is used for both level integration and
    * time-dependent error estimation.
    */
   virtual void
   postprocessFluxData(
      const std::shared_ptr<hier::PatchLevel>& level,
      const bool regrid_advance,
      const bool first_step,
      const bool last_step);

   /*
    * Copy time-dependent data from source space to destination space.
    */
   virtual void
   copyTimeDependentData(
      const std::shared_ptr<hier::PatchLevel>& level,
      const std::shared_ptr<hier::VariableContext>& src_context,
      const std::shared_ptr<hier::VariableContext>& dst_context);

   /**
    * Apply the standard AMR hyperbolic flux synchronization process preserve
    * conservation properties in the solution between the fine level and the
    * coarse level.  The sync_time argument indicates the time at which
    * the solution on the two levels is being synchronized.  The variable
    * coarse_sim_time indicates the previous simulation time on the
    * coarse level (recall the conservative difference will be repeated on the
    * coarse level during the synchronization process).  After the
    * synchronization, the flux and flux integral data storage is reset on
    * the levels.
    *
    * @pre fine_level
    * @pre coarse_level
    * @pre coarse_level->getLevelNumber() == (fine_level->getLevelNumber()-1))
    * @pre fine_level->getDim() == coarse_level->getDim()
    * @pre sync_time > coarse_sim_time
    */
   virtual void
   synchronizeLevelWithCoarser(
      const std::shared_ptr<hier::PatchLevel>& fine,
      const std::shared_ptr<hier::PatchLevel>& coarse,
      const double sync_time,
      const double coarse_sim_time);

   /*!
    * @brief Check the tags on a tagged level.
    *
    * This is a callback method that is used to check the values held in
    * user tag PatchData.  The tag data will contain the tags created
    * by either a gradient detector or Richardson extrapolation as well as
    * any tags added internally by the GriddingAlgorithm (for
    * example, buffering).
    *
    * @param[in] hierarchy
    * @param[in] level_number  Level number of the tagged level
    * @param[in] tag_index     Patch data index for user tags
    */
   virtual void
   checkUserTagData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int tag_index) const;

   /*!
    * @brief Check the tags on a newly-created level.
    *
    * This is a callback method that allow for checking tag values that
    * have been saved on a new level that has been created during
    * initialization or regridding.  The tag values will be the values
    * of the user tags on the coarser level, constant-refined onto the
    * cells of the new level.
    *
    * @param[in] hierarchy
    * @param[in] level_number   Level number of the new level
    * @param[in] tag_index      Patch data index for the new tags.
    */
   virtual void
   checkNewLevelTagData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int tag_index) const;

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int ALGS_HYPERBOLIC_LEVEL_INTEGRATOR_VERSION;

   /*
    * Record statistics on how many patches and cells were generated.
    */
   void
   recordStatistics(
      const hier::PatchLevel& patch_level,
      double current_time);

   /*
    * The patch strategy supplies the application-specific operations
    * needed to treat data on patches in the AMR hierarchy.
    */
   HyperbolicPatchStrategy* d_patch_strategy;

   /*
    * The gridding algorithm is used to generate some DLBG data for
    * levels not on a hiearchy.
    */
   std::shared_ptr<mesh::GriddingAlgorithm> d_gridding_alg;

   /*
    * The object name is used as a handle to databases stored in
    * restart files and for error reporting purposes.
    */
   std::string d_object_name;
   bool d_use_time_refinement;

   /*
    * Courant-Friedrichs-Levy parameters for time increment selection.
    */
   double d_cfl;
   double d_cfl_init;

   /*
    * Boolean flags for algorithm variations during time integration.
    *
    * d_lag_dt_computation indicates when time increment is computed for
    *                      next step on a level.  A value of true means
    *                      that the current solution values will be used to
    *                      compute dt.  A value of false means that dt will
    *                      be computed after the current solution is advanced
    *                      and the new solution is used to compute dt. The
    *                      default value is true.
    *
    * d_use_ghosts_for_dt  indicates whether the time increment computation
    *                      on a patch requires ghost cell data (e.g., if
    *                      boundary conditions are needed).  This value must
    *                      be consistent with the numerical routines used
    *                      in the hyperbolic patch strategy object to
    *                      calculate the time step size.  The default is true.
    *
    * d_use_flux_correction indicates whether the synchronization step
    *                       will replace coarse fluxes with integrated
    *                       finer level fluxes and repeat the conservative
    *                       update.  The default is true.
    */
   bool d_lag_dt_computation;
   bool d_use_ghosts_for_dt;
   bool d_use_flux_correction;

   /*
    * Boolean flags for indicated whether face or side data types are
    * used for fluxes (choice is determined by numerical routines in
    * hyperbolic patch model).
    */
   bool d_flux_is_face;
   bool d_flux_face_registered;
   bool d_flux_side_registered;

/*
 * The following communication algorithms and schedules are created and
 * maintained to manage inter-patch communication during AMR integration.
 * The algorithms are created in the class constructor.  They are initialized
 * when variables are "registered" are registered with the integrator.
 */

   /*
    * The "advance" schedule is used prior to advancing a level and
    * prior to computing dt at initialization. It must be reset each
    * time a level is regridded. All ghosts are filled with TIME_DEP
    * and INPUT data at specified time. TIME_DEP data in patch interiors
    * will be filled with CURRENT_VAR values.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_advance;
   std::vector<std::shared_ptr<xfer::RefineSchedule> > d_bdry_sched_advance;

   /*
    * The "advance new" schedule can be used twice during a time integration
    * cycle. The first is when ghost cell data is required during the
    * conservative difference process (i.e., d_use_ghosts_for_cons_diff
    * is true).  If this is the case, ghosts must be refilled before the
    * conservative difference on a coarser level during the refluxing
    * process can take place.  See synchronizeLevelWithCoarser in class
    * HyperbolicLevelIntegrator second occurs when the dt calculation is
    * not lagged and the physical boundary conditions are needed to compute dt
    * (i.e., (!d_lag_dt_computation && d_use_ghosts_for_dt_computation)
    * is true).  In either case, all ghosts are filled with TIME_DEP and INPUT
    * data at specified time.  TIME_DEP data in patch interiors will be filled
    * with values corresponding to NEW descriptor indices.  See notes
    * accompanying HyperbolicLevelIntegrator::advanceLevel.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_advance_new;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
   d_bdry_sched_advance_new;

   /*
    * The "advance old" algorithm is used to fill ghosts using time
    * interpolated data from OLD_VAR and NEW_VAR on the coarser hierarchy
    * level.  It is currently only used for advancing data on a temporary
    * level during the Richardson extrapolation algorithm. Use of OLD_VAR
    * data is required only when three time levels are used
    * (i.e. d_number_time_data_levels=3).
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_advance_old;

   /*
    * Coarsen algorithms for conservative data synchronization
    * (e.g., flux correction or refluxing).
    */
   std::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_fluxsum;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_sync_data;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_sync_initial_data;

   /*
    * Coarsen algorithms for Richardson extrapolation.
    */
   std::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_rich_extrap_init;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_rich_extrap_final;

   /*
    * Algorithm for filling a new patch level in the hierarchy.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_fill_new_level;

   /*
    * Number of levels of time-dependent data that must be maintained
    * on each patch level.  This value is used to coordinate the needs
    * of the time integration and the regridding process with the
    * patch data types and descriptor indices.
    */
   int d_number_time_data_levels;

   /*
    * hier::Variable contexts and lists of variables used for data management.
    * The contexts are set in the constructor.   Note that they must
    * be consistent with those defined by the concrete subclass of
    * the HyperbolicPatchStrategy object.  The variable lists
    * and component selectors are set in the registerVariable() function.
    */

   std::shared_ptr<hier::VariableContext> d_scratch;
   std::shared_ptr<hier::VariableContext> d_current;
   std::shared_ptr<hier::VariableContext> d_new;
   std::shared_ptr<hier::VariableContext> d_old;
   std::shared_ptr<hier::VariableContext> d_plot_context;

   std::list<std::shared_ptr<hier::Variable> > d_all_variables;
   std::list<std::shared_ptr<hier::Variable> > d_time_dep_variables;
   std::list<std::shared_ptr<hier::Variable> > d_flux_variables;
   std::list<std::shared_ptr<hier::Variable> > d_fluxsum_variables;

   /*
    * SCRATCH descriptor indices for (non-TEMPORARY) variables
    * (i.e., TIME_DEP, INPUT, FLUX). Note that these are used
    * to create scratch space before ghost cells are filled
    * on level prior to advancing the data.
    */
   hier::ComponentSelector d_saved_var_scratch_data;

   /*
    * SCRATCH descriptor indices for TEMPORARY variables.  Note that
    * these are used to create scratch space on a patch-by-patch basis.
    */
   hier::ComponentSelector d_temp_var_scratch_data;

   /*
    * CURRENT descriptor indices for TIME_DEP, INPUT, NO_FILL
    * variables.  Note that these are used to create storage for quantities
    * when new patches are made (e.g., during hierachy initialization,
    * before error estimation during regridding, after regridding new
    * patch levels, etc.).
    */
   hier::ComponentSelector d_new_patch_init_data;

   /*
    * NEW descriptor indices for TIME_DEP variables.  Note that these
    * are used to create space for new data before patch level is advanced.
    */
   hier::ComponentSelector d_new_time_dep_data;

   /*
    * Descriptor indices for FLUX quantities and integrals of fluxes
    * (used to accumulate flux information around fine patch boundaries).
    * Also, a boolean flag to track flux storage on level 0.
    */
   hier::ComponentSelector d_flux_var_data;
   hier::ComponentSelector d_fluxsum_data;
   bool d_have_flux_on_level_zero;

   /*
    * OLD descriptor indices for TIME_DEP variables.  Note that
    * these are used only when three time levels of data are used.
    */
   hier::ComponentSelector d_old_time_dep_data;

   /*
    * Option to distinguish tbox::MPI reduction costs from load imbalances
    * when doing performance timings.
    */
   bool d_distinguish_mpi_reduction_costs;

   static bool s_barrier_after_error_bdry_fill_comm;

   /*!
    * @brief For diagnostics: whether to separate major advanceLevel sections with MPI barriers.
    */
   bool d_barrier_advance_level_sections;

   /*
    * Timers interspersed throughout the class.
    */
   static std::shared_ptr<tbox::Timer> t_advance_bdry_fill_comm;
   static std::shared_ptr<tbox::Timer> t_error_bdry_fill_create;
   static std::shared_ptr<tbox::Timer> t_error_bdry_fill_comm;
   static std::shared_ptr<tbox::Timer> t_mpi_reductions;
   static std::shared_ptr<tbox::Timer> t_initialize_level_data;
   static std::shared_ptr<tbox::Timer> t_init_level_create_sched;
   static std::shared_ptr<tbox::Timer> t_init_level_fill_data;
   static std::shared_ptr<tbox::Timer> t_init_level_fill_interior;
   static std::shared_ptr<tbox::Timer> t_advance_bdry_fill_create;
   static std::shared_ptr<tbox::Timer> t_new_advance_bdry_fill_create;
   static std::shared_ptr<tbox::Timer> t_apply_gradient_detector;
   static std::shared_ptr<tbox::Timer> t_tag_cells;
   static std::shared_ptr<tbox::Timer> t_coarsen_rich_extrap;
   static std::shared_ptr<tbox::Timer> t_get_level_dt;
   static std::shared_ptr<tbox::Timer> t_get_level_dt_sync;
   static std::shared_ptr<tbox::Timer> t_advance_level;
   static std::shared_ptr<tbox::Timer> t_advance_level_integrate;
   static std::shared_ptr<tbox::Timer> t_advance_level_pre_integrate;
   static std::shared_ptr<tbox::Timer> t_advance_level_post_integrate;
   static std::shared_ptr<tbox::Timer> t_advance_level_patch_loop;
   static std::shared_ptr<tbox::Timer> t_new_advance_bdry_fill_comm;
   static std::shared_ptr<tbox::Timer> t_patch_num_kernel;
   static std::shared_ptr<tbox::Timer> t_preprocess_flux_data;
   static std::shared_ptr<tbox::Timer> t_postprocess_flux_data;
   static std::shared_ptr<tbox::Timer> t_advance_level_sync;
   static std::shared_ptr<tbox::Timer> t_advance_level_compute_dt;
   static std::shared_ptr<tbox::Timer> t_copy_time_dependent_data;
   static std::shared_ptr<tbox::Timer> t_std_level_sync;
   static std::shared_ptr<tbox::Timer> t_sync_new_levels;
   static std::shared_ptr<tbox::Timer> t_barrier_after_error_bdry_fill_comm;
   static std::shared_ptr<tbox::Timer> t_sync_initial_comm;
   static std::shared_ptr<tbox::Timer> t_sync_initial_create;
   static std::shared_ptr<tbox::Timer> t_coarsen_fluxsum_create;
   static std::shared_ptr<tbox::Timer> t_coarsen_fluxsum_comm;
   static std::shared_ptr<tbox::Timer> t_coarsen_sync_create;
   static std::shared_ptr<tbox::Timer> t_coarsen_sync_comm;

#ifdef HLI_RECORD_STATS
   /*
    * Statistics on number of cells and patches generated.
    */
   static std::vector<std::shared_ptr<tbox::Statistic> > s_boxes_stat;
   static std::vector<std::shared_ptr<tbox::Statistic> > s_cells_stat;
   static std::vector<std::shared_ptr<tbox::Statistic> > s_timestamp_stat;
#endif

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*
    * Static initialization and cleanup handler.
    */
   static tbox::StartupShutdownManager::Handler
      s_initialize_handler;
};

}
}

#endif
