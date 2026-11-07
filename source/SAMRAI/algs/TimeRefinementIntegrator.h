/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Time integration manager for AMR with local time stepping.
 *
 ************************************************************************/

#ifndef included_algs_TimeRefinementIntegrator
#define included_algs_TimeRefinementIntegrator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/algs/TimeRefinementIntegratorConnectorWidthRequestor.h"
#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <iostream>
#include <memory>

namespace SAMRAI {
namespace algs {

/**
 * Class TimeRefinementIntegrator manages time integration over an
 * AMR patch hierarchy using local time refinement for finer hierarchy levels.
 * This class orchestrates hierarchy construction, data advancement and
 * synchronization, and the dynamic grid refinement processes.  The basic
 * ideas behind these algorithms are described in several sources on
 * structured adaptive mesh refinement.  See Berger and Colella, J. Comp.
 * Phys. (82)1:64-84, 1989 for an introduction to algorithm.  See Trangenstein,
 * SIAM J. Sci. Comput. 16(4):819-839, 1995, or Hornung, PhD thesis, Dept.
 * of Mathematics, Duke University, 1994 for further discussion.
 *
 * This class can be used in two different modes:  refined timestepping,
 * which divides they hierarchy's timestep into smaller timesteps on finer
 * levels, or synchronized timestepping, which advances all levels in a
 * hierarchy by the same timestep.  The mode that is used is determined
 * by querying the level integrator that is passed into the constructor.
 * One and only one mode can be used for each instantiation of
 * this class.
 *
 * The algorithm requires that integration steps on different levels are
 * interleaved since the time increment used on each level is determined
 * by the spatial resolution of the mesh on that level (e.g., CFL condition).
 * Generally, when using refined timestepping, coarser levels use larger time
 * increments than finer levels.  Thus, data must be synchronized between
 * levels and the dynamic regridding process must be coordinated with
 * the time stepping sequence.
 *
 * The routines in this class are implemented in a manner that is generic
 * with respect to the details of the level integration and regridding
 * methods, and the discrete equations being solved.  Thus, the class may
 * be employed for a variety of applications.  Upon construction, an object
 * of this class is configured with routines suitable for a given problem.
 * The TimeRefinementLevelStrategy data member supplies routines
 * for advancing levels and synchronizing data on different levels during
 * time integration.  The mesh::GriddingAlgorithmStrategy data member provides
 * operations that construct and dynamically reconfigure the patch hierarchy.
 * The collaboration between this class and each of those objects follows
 * the Strategy design pattern.
 *
 * Initialization begins by setting data on the coarsest AMR hierarchy
 * level.  Then, each successively finer level is created and initialized
 * after invoking the regridding procedures on the most recently initialized
 * level.  This process is performed until either a maximum number of levels
 * is reached or no further refinement is needed.
 *
 * Time integration is performed by invoking a recursive advance procedure on
 * the coarsest AMR hierarchy level.  On a level, data is integrated to
 * a given point using a sequence of integration steps, where the size of
 * each time increment depends on the problem being solved.  After each
 * step on a level, the next finer level (if it exists) is integrated to
 * the same time using a sequence of time steps appropriate for the level.
 * This class may dynamically adjust the time step sequence used on each
 * level during the data advance process depending on requirements of the
 * integrator and information about stable time step size.  Dynamic
 * mesh regridding is invoked during the integration process so that
 * time integration, data synchronization, and mesh movement are coordinated
 * properly.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *    - \b    regrid_interval
 *       when using synchronized timestepping, number of timesteps between each
 *        regrid of the hierarchy
 *
 *    - \b    start_time
 *       start time for the simulation.
 *
 *    - \b    end_time
 *       end time for the simulation.
 *
 *    - \b    grow_dt
 *       maximum factor by which each succesive time increment may grow
 *       (typically >= 1.0).
 *
 *    - \b    max_integrator_steps
 *       maximum number of timesteps performed on the coarsest hierarchy level
 *       during the simulation.
 *
 *    - \b    tag_buffer
 *       array of integer values (one for each level that may be refined)
 *       representing the number of cells by which tagged cells are buffered
 *       before clustering into boxes.
 *
 * Note that the input values for regrid_interval, end_time, grow_dt,
 * max_integrator_steps, and tag_buffer override values read in from restart.
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
 *     <td>regrid_interval</td>
 *     <td>int</td>
 *     <td>1</td>
 *     <td>>=1</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>start_time</td>
 *     <td>double</td>
 *     <td>none</td>
 *     <td>start_time >=0</td>
 *     <td>req</td>
 *     <td>May not be modified by input db on restart</td>
 *   </tr>
 *   <tr>
 *     <td>end_time</td>
 *     <td>double</td>
 *     <td>none</td>
 *     <td>end_time >= start_time</td>
 *     <td>req</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>grow_dt</td>
 *     <td>double</td>
 *     <td>1.0</td>
 *     <td>>0</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>max_integrator_steps</td>
 *     <td>int</td>
 *     <td>none</td>
 *     <td>>=0</td>
 *     <td>req</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>tag_buffer</td>
 *     <td>array of ints</td>
 *     <td>regrid_interval value for corresponding level</td>
 *     <td>all values >= 0</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * A sample input file entry might look like:
 *
 * @code
 *    start_time            = 0.e0      // initial simulation time
 *    end_time              = 10.e0     // final simulation time
 *    grow_dt               = 1.1e0     // growth factor for timesteps
 *    max_integrator_steps  = 50        // max number of simulation timesteps
 *    tag_buffer            = 1,1,1,1   // a max of 4 finer levels in hierarchy
 * @endcode
 *
 * @see TimeRefinementLevelStrategy
 * @see mesh::GriddingAlgorithmStrategy
 */

class TimeRefinementIntegrator:
   public tbox::Serializable
{
public:
   /**
    * The constructor for TimeRefinementIntegrator initializes the
    * time stepping parameters needed to integrate the levels in the AMR
    * hierarchy.   Some data is set to default values; others are read
    * from the specified input database and the restart database
    * corresponding to the specified object_name.  Consult top of
    * this header file for further details.
    *
    * Note that this object also invokes the variable creation and
    * registration process in the level strategy.
    *
    * @pre !object_name.empty()
    * @pre hierarchy
    * @pre level_integrator
    * @pre gridding_algorithm
    */
   TimeRefinementIntegrator(
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const std::shared_ptr<TimeRefinementLevelStrategy>& level_integrator,
      const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_algorithm);

   /**
    * The destructor for TimeRefinementIntegrator unregisters
    * the integrator object with the restart manager.
    */
   virtual ~TimeRefinementIntegrator();

   /*!
    * Set AMR patch hierarchy configuration and data at start of simulation.
    * If the run is begun from a restart file, the hierarchy and data
    * are read from the hierarchy database.  Otherwise, the hierarchy
    * and data are initialized by the gridding algorithm data member.
    * In this case, the coarsest level is constructed and initialized.
    * Then, error estimation is performed to determine if and where it
    * should be refined.  Successively finer levels are created and
    * initialized until the maximum allowable number of levels is achieved
    * or no further refinement is needed.  The double return value is the
    * time increment for the first data advance step on the coarsest
    * hierarchy level (i.e., level 0).
    *
    * This function assumes that the hierarchy exists, but that it contains
    * no patch levels, when it is called.  On return from this function, the
    * initial hierarchy configuration and simulation data is set properly for
    * the advanceHierarchy() function to be called.  In particular, on each
    * level constructed only the data needed for initialization exists.
    */
   double
   initializeHierarchy();

   /**
    * Advance each level in the hierarchy through the given time increment
    * and return an appropriate time increment for subsequent advances of the
    * coarsest hierarchy level (level 0).  The boolean argument indicates
    * whether the coarsest hierarchy level (i.e., level 0) should be load
    * balanced before the levels are advanced.  In general, the problem
    * domain (determined by the union of patches on level 0) does not change
    * once set.  However, the boolean flag here allows one to reconfigure
    * the patches on the coarsest level which constitute this union.  This
    * may be required depending on a dynamic change of the work load.
    * By default, the level will not be subject to load balancing.
    *
    * This function assumes that all data on each level in the hierarchy
    * has been set and that only the data need for initialization exists
    * on each level (as opposed to both current and new data, for example).
    * Upon return from this function, the simulation data on each hierarchy
    * levels is advanced through the time increment dt.  In addition, data on
    * all hierarchy levels has been synchronized so that it is consistent at
    * the new simulation time (where this synchronization process is defined
    * by the level strategy).  Thus, the data is set properly for any
    * subsequent calls to this function.
    *
    * @pre dt >= 0
    */
   double
   advanceHierarchy(
      const double dt,
      const bool rebalance_coarsest = false);

   /**
    * Return true if the current step count for the level indicates
    * that regridding should occur.  In particular, true is returned
    * if both the level allows refinement and the step count is an
    * integer multiple of the regrid step interval.
    * Otherwise, false is returned.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   bool
   atRegridPoint(
      const int level_number) const;

   /**
    * Return current integration time for coarsest hierarchy level.
    */
   double
   getIntegratorTime() const
   {
      return d_integrator_time;
   }

   /**
    * Return initial integration time.
    */
   double
   getStartTime() const
   {
      return d_start_time;
   }

   /**
    * Return final integration time.
    */
   double
   getEndTime() const
   {
      return d_end_time;
   }

   /**
    * Return integration step count for entire hierarchy
    * (i.e., number of steps taken by the coarsest level).
    */
   int
   getIntegratorStep() const
   {
      return d_step_level[0];
   }

   /**
    * Return maximum number of integration steps allowed for entire
    * hierarchy (i.e., steps allowed on coarsest level).
    */
   int
   getMaxIntegratorSteps() const
   {
      return d_max_steps_level[0];
   }

   /**
    * Return true if any steps remain in current step sequence on level
    * (i.e., before it will synchronize with some coarser level).
    * Return false otherwise.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   bool
   stepsRemaining(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      if (level_number == 0) {
         return !d_level_0_advanced;
      } else {
         return d_step_level[level_number] < d_max_steps_level[level_number];
      }
   }

   /**
    * Return true if any integration steps remain, false otherwise.
    */
   bool
   stepsRemaining() const
   {
      return d_step_level[0] < d_max_steps_level[0];
   }

   /**
    * Return current time increment used to advance level.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   getLevelDtActual(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      return d_dt_actual_level[level_number];
   }

   /**
    * Return maximum time increment currently allowed on level.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   getLevelDtMax(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      return d_dt_max_level[level_number];
   }

   /**
    * Return current simulation time for level.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   getLevelSimTime(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      return d_level_sim_time[level_number];
   }

   /**
    * Return step count for current integration sequence on level.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   int
   getLevelStep(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      if (level_number == 0) {
         return d_level_0_advanced ? 1 : 0;
      } else {
         return d_step_level[level_number];
      }
   }

   /**
    * Return maximum number of time steps allowed on level in
    * current integration step sequence.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   int
   getLevelMaxSteps(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      if (level_number == 0) {
         return 1;
      } else {
         return d_max_steps_level[level_number];
      }
   }

   /**
    * Return const pointer to patch hierarchy managed by integrator.
    */
   const std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const
   {
      return d_patch_hierarchy;
   }

   /**
    * Return pointer to level integrator.
    */
   std::shared_ptr<TimeRefinementLevelStrategy>
   getLevelIntegrator() const
   {
      return d_refine_level_integrator;
   }

   /**
    * Return pointer to gridding algorithm object.
    */
   std::shared_ptr<mesh::GriddingAlgorithmStrategy>
   getGriddingAlgorithm() const
   {
      return d_gridding_algorithm;
   }

   /**
    * Return true if current step on level is first in current step
    * sequence; otherwise return false.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   bool
   firstLevelStep(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      if (level_number == 0) {
         return !d_level_0_advanced;
      } else {
         return d_step_level[level_number] <= 0;
      }
   }

   /**
    * Return true if current step on level is last in current step
    * sequence; otherwise return false.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   bool
   lastLevelStep(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
         (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
      if (level_number == 0) {
         return d_level_0_advanced;
      } else {
         return d_step_level[level_number] >= d_max_steps_level[level_number];
      }
   }

   /**
    * Set the regrid interval to a new value.  This may only be used
    * when using synchronized timestepping.
    *
    * @pre !d_use_refined_timestepping
    */
   void
   setRegridInterval(
      const int regrid_interval);

   /**
    * Print data representation of this object to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Print time stepping data for a single level to given output stream.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   printDataForLevel(
      std::ostream& os,
      const int level_number) const;

   /**
    * Write object state out to the given restart database.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
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
   static const int ALGS_TIME_REFINEMENT_INTEGRATOR_VERSION;

   /*
    * Initialize data on given level.  If the level can be refined, a problem-
    * dependent error estimation procedure is invoked to determine whether
    * further refinement is needed.  If a new level is created, this function
    * is called recursively to initialize the next finest level.
    */
   void
   initializeRefinedTimesteppingLevelData(
      const int level_number);
   void
   initializeSynchronizedTimesteppingLevelData(
      const int level_number);

   /*
    * Advance the data on the level to the specified time using a dynamically
    * adjusted sequence of time increments.  If any finer levels exist
    * in the hierarchy when this function is called or are generated during
    * the regridding process, they will be advanced to the specified time
    * as well.  This function is recursive.  After a single timestep is
    * performed on each level, all finer levels are advanced to the new
    * integration time before another timestep is taken on the original level.
    */
   void
   advanceRecursivelyForRefinedTimestepping(
      const int level_number,
      const double end_time);

   double
   advanceForSynchronizedTimestepping(
      const double dt);

   /*
    * Determine the next stable time increment (dt) and adjust the step
    * sequence if necessary for the given level.  The computed dt will
    * be less than or equal to the specified bound and the time remaining.
    * In adjusting the step sequence, an attempt is made to partition the
    * remaining time interval into a sequence of equal time increments.
    * However, the total number of time steps in the sequence for the
    * level must satisfy any constraints placed on it by the regridding
    * procedure.
    */
   bool
   findNextDtAndStepsRemaining(
      const int level_number,
      const double time_remaining,
      const double dt_bound);

   /*
    * Return true if the this level can be regridded at the current step
    * and the next coarser level can be regridded too.  Otherwise,
    * false is returned.
    */
   bool
   coarserLevelRegridsToo(
      const int level_number) const;

   /*
    * Read input data from specified input database and initialize class
    * members.  The argument is_from_restart should be set to true if the
    * simulation is from restart.  Otherwise, it should be set to false.
    *
    * If the simulation is not from restart, read in start_time, end_time,
    * grow_dt, max_integrator_step, and possibly tag_buffer
    * from the database.
    *
    * If the simulation is from restart, then only read in end_time,
    * grow_dt, max_integrator_step and tag_buffer if they are
    * found in the input database.
    */
   virtual void
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
    */
   virtual void
   getFromRestart();

   /*
    * The object name is used as a handle to databases stored in
    * restart files and for error reporting purposes.
    */
   std::string d_object_name;

   /*
    * Pointers to the patch hierarchy, level integration and gridding
    * algorithm objects associated with this time integration object.
    * The level integrator defines operations for advancing data on
    * individual levels in the AMR patch hierarchy.  The gridding algorithm
    * provides grid generation and regridding routines for the AMR hierarchy.
    */
   std::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
   std::shared_ptr<TimeRefinementLevelStrategy> d_refine_level_integrator;
   std::shared_ptr<mesh::GriddingAlgorithmStrategy> d_gridding_algorithm;

   /*
    */
   bool d_use_refined_timestepping;

   /*
    * Integrator data read from input or set at initialization.
    */
   double d_start_time;
   double d_end_time;
   double d_grow_dt;

   /*
    * The regrid interval indicates the number of integration steps taken
    * on a level between successive invocations of the regridding process
    * on that level.  In general, this class enforces the constraint that
    * each synchronization time between two successive hierarchy levels
    * will always be a potential regrid point for the coarser of the two
    * levels.  Specifically, it sets the regrid interval for each level
    * to be the greatest common divisor between the entries in the grid
    * refinement ratio vector between the level and the next coarsest level
    * in the hierarchy.  The regrid interval for level 0 is set equal to
    * that for level 1.  In the future, users may be able to specify
    * this value in the input file.
    */
   std::vector<int> d_regrid_interval;

   /*
    * The tag buffer indicates the number of cells on each level by which
    * tagged cells will be buffered after they have selected for refinement.
    * These values are passed into the gridding algorithm routines during
    * hierarchy construction and regridding.  The tag buffer helps to
    * guarantee that refined cells near important features in the solution
    * will remain refined until the level is regridded next.
    *
    * Important note: these values may be specified in the input file.
    * If not, default values are set based on the regrid intervals.
    * However, if the user attempts to specify these values, care must
    * be taken to assure that improper tag buffering will not degrade the
    * calculation.
    */
   std::vector<int> d_tag_buffer;

   /*
    * Integrator data that evolves during time integration and maintains
    * the state of the timestep sequence over the levels in the AMR hierarchy.
    */
   double d_integrator_time;
   bool d_just_regridded;
   int d_last_finest_level;
   std::vector<double> d_level_old_old_time;
   std::vector<double> d_level_old_time;
   std::vector<double> d_level_sim_time;
   std::vector<double> d_dt_max_level;
   std::vector<double> d_dt_actual_level;
   std::vector<int> d_step_level;
   std::vector<int> d_max_steps_level;
   bool d_level_0_advanced;
   bool d_hierarchy_advanced;

   double d_dt;

   TimeRefinementIntegratorConnectorWidthRequestor d_connector_width_requestor;

   bool d_barrier_and_time;

   /*
    * tbox::Timer objects for performance measurement.
    */
   static std::shared_ptr<tbox::Timer> t_initialize_hier;
   static std::shared_ptr<tbox::Timer> t_advance_hier;
   static std::shared_ptr<tbox::Timer> t_advance_level;

   // The following are not implemented:
   TimeRefinementIntegrator(
      const TimeRefinementIntegrator&);
   TimeRefinementIntegrator&
   operator = (
      const TimeRefinementIntegrator&);

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
