/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Timer class to track elapsed time in portions of a program.
 *
 ************************************************************************/

#ifndef included_tbox_Timer
#define included_tbox_Timer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Clock.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/PIO.h"

#include <string>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace tbox {

class TimerManager;

/**
 * Class Timer holds the exclusive and total start, stop, and elapsed
 * time for timers instrumented in SAMRAI.  Total time is simply the time
 * between calls to the start() and stop() functions.  Exclusive time is
 * applicable if there are nested timers called.
 *
 * System and user start and end times are stored as variables of type
 * clock_t, defined in the sys/times.h include file.  A detailed explanation
 * of the structures used to store system and user times is given in the
 * header for the Clock class. This routine simply accesses the functions
 * specified in that class.
 *
 * Wallclock time may be computed by the systems internal clocks which require
 * an object of type clock_t, or by SAMRAI_MPI::Wtime() if the code is linked
 * to MPI libraries.
 *
 * In addition to running or not running, a timer may be active or inactive.
 * An inactive timer is one that is created within a program but will never
 * be turned on or off because it is either not specified as active in
 * an input file or it was not explicitly made active by the user.  When
 * a timer is created, it is active by default.
 *
 * Note that the constructor is protected so that timer objects can only
 * be created by the TimerManager class.
 *
 * @see TimerManager
 */

class Timer
{
   friend class TimerManager;
public:
   /**
    * Empty destructor for Timer class.
    */
   ~Timer();

   /**
    * Return string name for timer.
    */
   const std::string&
   getName() const
   {
      return d_name;
   }

   /**
    * Start the timer if active.
    *
    * @pre !isActive() || !isRunning()
    */
   void
   start();

   /**
    * Stop the timer if active.
    *
    * @pre !isActive() || isRunning()
    */
   void
   stop();

   /**
    * If active, SAMRAI_MPI::getSAMRAIWorld().Barrier() then start the timer.
    */
   void
   barrierAndStart();

   /**
    * If active, SAMRAI_MPI::getSAMRAIWorld().Barrier() then stop the timer.
    */
   void
   barrierAndStop();

   /**
    * Start exclusive time.
    */
   void
   startExclusive();

   /**
    * Stop exclusive time.
    */
   void
   stopExclusive();

   /**
    * Reset the state of the timing information.
    */
   void
   reset();

   /**
    * Return total system time (between starts and stops)
    */
   double
   getTotalSystemTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_system_total / Clock::getClockCycle();

#else
      return 0.0;

#endif
   }

   /**
    * Return total user time
    */
   double
   getTotalUserTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_user_total / Clock::getClockCycle();

#else
      return 0.0;

#endif
   }

   /**
    * Return total wallclock time
    */
   double
   getTotalWallclockTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_wallclock_total;

#else
      return 0.0;

#endif
   }

   /**
    * Return max wallclock time
    */
   double
   getMaxWallclockTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_max_wallclock;

#else
      return 0.0;

#endif
   }

   /**
    * Return exclusive system time.
    */
   double
   getExclusiveSystemTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_system_exclusive / Clock::getClockCycle();

#else
      return 0.0;

#endif
   }

   /**
    * Return exclusive user time.
    */
   double
   getExclusiveUserTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_user_exclusive / Clock::getClockCycle();

#else
      return 0.0;

#endif
   }

   /**
    * Return exclusive wallclock time.
    */
   double
   getExclusiveWallclockTime() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_wallclock_exclusive;

#else
      return 0.0;

#endif
   }

   /**
    * Return true if the timer is active; false otherwise.
    */
   bool
   isActive() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_is_active;

#else
      return false;

#endif
   }

   /**
    * Return true if timer is running; false otherwise.
    */
   bool
   isRunning() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_is_running;

#else
      return false;

#endif
   }

   /**
    * Return number of accesses to start()-stop() functions for the
    * timer.
    */
   int
   getNumberAccesses() const
   {
#ifdef ENABLE_SAMRAI_TIMERS
      return d_accesses;

#else
      return 0;

#endif
   }

   /**
    * Compute load balance efficiency based on wallclock (non-exclusive)
    * time.
    */
   double
   computeLoadBalanceEfficiency();

   /**
    * Compute max wallclock time based on total (non-exclusive) time.
    */
   void
   computeMaxWallclock();

   /**
    * Write timer data members to restart database.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<Database>& restart_db) const;

   /**
    * Read restarted times from restart database.
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<Database>& restart_db);

protected:
   /**
    * The constructor for the Timer class sets timer name string
    * and integer identifiers, and initializes the timer state.
    */
   explicit Timer(
      const std::string& name);

   /*
    * Set this timer object to be a active or inactive.
    */
   void
   setActive(
      bool is_active)
   {
#ifdef ENABLE_SAMRAI_TIMERS
      d_is_active = is_active;
#else
      NULL_USE(is_active);
#endif
   }

   /**
    * Add Timer that running concurrently with this one.
    */
   void
   addConcurrentTimer(
      const Timer& timer)
   {
#ifdef ENABLE_SAMRAI_TIMERS
      if (!isConcurrentTimer(timer)) {
         d_concurrent_timers.push_back(&timer);
      }
#else
      NULL_USE(timer);
#endif
   }

   /**
    * Return if the timer is running concurrently with this one.
    */
   bool
   isConcurrentTimer(
      const Timer& timer) const;

private:
   // Unimplemented default constructor.
   Timer();

   // Unimplemented copy constructor.
   Timer(
      const Timer& other);

   // Unimplemented assignment operator.
   Timer&
   operator = (
      const Timer& rhs);

   /*
    * Class name, id, and concurrent timer flag.
    */
   std::string d_name;
   std::vector<const Timer *> d_concurrent_timers;

   bool d_is_running;
   bool d_is_active;

   /*
    *  Total times (non-exclusive)
    */
   double d_user_total;
   double d_system_total;
   double d_wallclock_total;

   /*
    *  Exclusive times
    */
   double d_user_exclusive;
   double d_system_exclusive;
   double d_wallclock_exclusive;

   /*
    *  Cross processor times (i.e. determined across processors)
    */
   double d_max_wallclock;

   /*
    *  Timestamps.  User and system times are stored as type clock_t.
    *  Wallclock time is also stored as clock_t unless the library has
    * been compiled with MPI.  In this case, the wall time is stored
    * as type double.
    */
   clock_t d_user_start_total;
   clock_t d_user_stop_total;
   clock_t d_system_start_total;
   clock_t d_system_stop_total;
   clock_t d_user_start_exclusive;
   clock_t d_user_stop_exclusive;
   clock_t d_system_start_exclusive;
   clock_t d_system_stop_exclusive;
   double d_wallclock_start_total;
   double d_wallclock_stop_total;
   double d_wallclock_start_exclusive;
   double d_wallclock_stop_exclusive;

   /*
    * Counter of number of times timers start/stop
    * are accessed.
    */
   int d_accesses;

   static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT;

   /*
    * Static integer constant describing this class's version number.
    */
   static const int TBOX_TIMER_VERSION;
};

}
}

#endif
