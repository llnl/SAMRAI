/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to level routines for time-refinement integrator.
 *
 ************************************************************************/

#ifndef included_algs_TimeRefinementLevelStrategy
#define included_algs_TimeRefinementLevelStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Utilities.h"


#include <vector>
#include <memory>

namespace SAMRAI {
namespace algs {

/**
 * Class TimeRefinementLevelStrategy is an abstract base class that
 * defines the interface to level integration and synchronization routines
 * needed by the hierarchy integration class TimeRefinementIntegrator.
 * In particular, this class insulates the hierarchy integrator from the
 * routines that manipulate data on the hierarchy levels in a problem-specific
 * fashion.  When the AMR hierarchy integration and regridding sequence
 * provided by the class TimeRefinementIntegrator are appropriate
 * for some computational problem, a subclass of this base class can be
 * used to provide the necessary operations to the hierarchy integrator.
 * That is, a TimeRefinementIntegrator object may be configured with
 * a concrete implementation of this base class by passing the concrete
 * object into the to the time refinement integrator constructor.
 *
 * @see TimeRefinementIntegrator
 */

class TimeRefinementLevelStrategy
{
public:
   /**
    * Default constructor for TimeRefinementLevelStrategy.
    */
   TimeRefinementLevelStrategy();

   /**
    * Virtual destructor for TimeRefinementLevelStrategy.
    */
   virtual ~TimeRefinementLevelStrategy();

   /**
    * Initialize the state of the integrator that performs time
    * advances on the levels of the hierarchy.  Typically, this
    * involves setting up information, such as communication algorithms,
    * to manage variable storage. The pointer to the gridding algorithm
    * is provided so that the integrator may access information about
    * regridding procedures or the structure the hierarchy, which is not
    * yet created.
    */
   virtual void
   initializeLevelIntegrator(
      const std::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg) = 0;

   /**
    * Return appropriate time increment for given level in the patch
    * hierarchy.  This routine is called during the initial
    * generation of the AMR patch hierarchy and possibly during regridding
    * if time regridding uses a time advance.  It should be assumed that
    * the only data that exists on the level when this routine is called
    * is that which is needed to initialize the level.  The initial_time
    * boolean flag is true if this routine is called at the initial
    * simulation time (i.e., when hierarchy is generated for first time);
    * otherwise (e.g., at an advance step) it is false.
    *
    * The recompute_dt option specifies whether to compute
    * the timestep using the current level data or to return the value
    * stored by the time integrator. The default true setting means
    * the timestep will be computed if no value is supplied.
    */
   virtual double
   getLevelDt(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double dt_time,
      const bool initial_time) = 0;

   /**
    * Return the maximum allowable time increment for level in the hierarchy
    * with the given level number.  The point of this routine is to determine
    * the increment for that level based on the current time increment used
    * on the next coarser level.  The coarse dt is the current time step
    * size for the next coarser level.  The ratio is the mesh refinement
    * ratio between the two levels.
    *
    * If the concrete implentation of this class only supports synchronized
    * timestepping, this should return the time increment that is applicable on
    * all levels of the hierarchy.
    */
   virtual double
   getMaxFinerLevelDt(
      const int finer_level_number,
      const double coarse_dt,
      const hier::IntVector& ratio) = 0;

   /**
    * Advance data on all patches on specified patch level from current time
    * (current_time) to new time (new_time).   The boolean value first_step
    * indicates whether the advance step is the first in a time sequence
    * on the level.  The boolean value last_step indicates whether the
    * advance step is the last in a time sequence on the level.   Usually,
    * the timestep sequence refers to the steps taken to advance the solution
    * through the time interval of the most recent advance on the next coarser
    * level, if such a level exists.  The boolean regrid_advance is false
    * when the advance is part of the actual hierarchy integration process
    * and true when the advance is called during time-dependent regridding
    * (e.g., when using Richardson extrapolation).  The default value is false.
    * The last boolean argument is true when the level is in the hierarchy,
    * and false otherwise.  The default value is true.  Usually, this value
    * is false only during time-dependent regridding operations performed
    * on some temporary level; thus, a schedule must be generated for the
    * level before the advance can occur, for example.
    *
    * When this function is called, the level data required to begin the
    * advance must be allocated and be defined appropriately.  Typically,
    * this is equivalent to what is needed to initialize a new level after
    * regridding.  Upon exiting this routine, both current and new data may
    * exist on the level.  This data is needed until level synchronization
    * occurs, in general. Current and new data may be reset by calling
    * the member function resetTimeDependentData().
    *
    * This routine is called from two different points within the
    * TimeRefinementIntegrator class: during the regular time
    * advance sequence, and at the initial simulation time.  The second
    * call is made to advance the solution on a coarser level ahead in time to
    * provide time-dependent boundary values for some finer level when
    * time-dependent regridding is used.  In the first case, the values of
    * the boolean flags are:
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
    */
   virtual double
   advanceLevel(
      const std::shared_ptr<hier::PatchLevel>& level,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const double current_time,
      const double new_time,
      const bool first_step,
      const bool last_step,
      const bool regrid_advance = false) = 0;

   /**
    * Synchronize data on specified patch levels in AMR hierarchy at the
    * given synchronization time.  The array of time values provides the
    * previous integration time for each level involved in the synchronization.
    * In other words, (sync_time - old_times[ln]) is the most recent time
    * increment used to advance data on level ln.  These times are used when
    * the synchronization process requires re-integration of the data.  Note
    * that other synchronization routines are defined below for other points
    * in the hierarchy integration sequence.
    *
    * When this routine is called, both current and new data may exist on each
    * level involved in the synchronization.  The new data on each level
    * corresponds to the synchronization time.   Each entry in the array
    * of time values specifies the time to which the current data on each
    * level corresponds.  It is assumed that this routine will reset the
    * synchronized data on each level so that only the current data will
    * exist on each level when done.
    *
    * Note that this routine is distinct from the synchronizeNewLevels()
    * function below.  This routine is used to synchronize levels during
    * the time integration process.  The other routine is used to synchronize
    * new levels in the hierarchy, either at initialization time or after
    * regridding.
    *
    * In the case of the time refinement integrator using synchronized
    * timestepping, the old_times argument should be an array containing the
    * same time value for all levels, since all levels are advanced with the
    * same timestep.
    */
   virtual void
   standardLevelSynchronization(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const std::vector<double>& old_times) = 0;

   /**
    * Synchronize specified levels after regridding has occurred or during
    * initial construction of the AMR patch hierarchy.  Note that this
    * synchronization may be different than the standard time-dependent
    * synchronization (above in standardLevelSynchronization()) depending
    * on the level integration algorithm.
    *
    * Before this routine is called, all time-dependent data on all levels
    * involved in the synchronization has been reset.  Thus, this routine
    * must only synchronize the current data on each level.  On return from
    * this function, only the current data on each level must be present.
    */
   virtual void
   synchronizeNewLevels(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level,
      const double sync_time,
      const bool initial_time) = 0;

   /**
    * Reset time-dependent data storage for the specified patch level.  This
    * routine is called when the current level data is no longer needed
    * and it is appropriate to replace the current data with the new data
    * on the level, if such data exists.
    */
   virtual void
   resetTimeDependentData(
      const std::shared_ptr<hier::PatchLevel>& level,
      const double new_time,
      const bool can_be_refined) = 0;

   /**
    * Reset data on the patch level to state before time advance.
    * This is needed whenever the solution on a level is advanced beyond
    * the new time on a level.  For example, during time-dependent regridding
    * (Richardson extrapolation) or initialization, it is necessary to
    * do such an advance.  This routine is called to discard the new solution
    * data so that subsequent calls to advance are provided proper data at the
    * correct time.
    */
   virtual void
   resetDataToPreadvanceState(
      const std::shared_ptr<hier::PatchLevel>& level) = 0;

   /**
    * Return true if the implementation of this class is constructed
    * to use refined timestepping, and false otherwise.
    */
   virtual bool
   usingRefinedTimestepping() const = 0;

private:
};

}
}
#endif
