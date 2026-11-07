/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple utility class for interfacing with system clock
 *
 ************************************************************************/

#ifndef included_tbox_Clock
#define included_tbox_Clock

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#ifdef HAVE_CTIME
#include <ctime>
#endif

#ifdef SAMRAI_HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef SAMRAI_HAVE_UNISTD_H
#include <unistd.h>
#endif

namespace SAMRAI {
namespace tbox {

/**
 * Class Clock serves as a single point of access for system clock
 * information.  System and user time are computed via the POSIX compliant
 * times() function.  This is described on p. 137, Lewine, POSIX programmers
 * guide, 1992.  The methods and structs used in this utility are defined
 * in <sys/times.h>.  Start and end times are stored as variables of type
 * clock_t.  A clock_t value can be converted to seconds by dividing by
 * CLK_TCK (which is defined in <sys/times.h>).  Different systems may use
 * different CLK_TCK.  Time is accessed by calling the times() function which
 * takes as an argument a reference to an object of type struct tms.  This
 * object will record the system and user time (obj.tms_utime \&
 * obj.tms_stime) and will return the time since the system was started.
 *
 * The return value from the call to times() can be used to compute elapsed
 * wallclock time.  Alternatively, one can use SAMRAI_MPI::Wtime() if MPI
 * libraries are included.  Two methods are defined for accessing system time
 * - one that has a clock_t struct argument for wallclock time (the non-MPI
 * case) and one that has a double argument to record the value of
 * SAMRAI_MPI::Wtime().
 *
 * Computing user/system/wallclock time with the times() function is performed
 * as follows:
 * \verbatim
 *    struct tms buffer;
 *    clock_t wtime_start = times(&buffer);
 *    clock_t stime_start = buffer.tms_stime;
 *    clock_t utime_start = buffer.tms_utime;
 *     (do some computational work)
 *    clock_t wtime_stop  = times(&buffer);
 *    clock_t stime_stop  = buffer.tms_stime;
 *    clock_t utime_stop  = buffer.tms_utime;
 *    double wall_time   = double(wtime_stop-wtime_start)/double(CLK_TCK);
 *    double user_time   = double(utime_stop-utime_start)/double(CLK_TCK);
 *    double sys_time    = double(stime_stop-stime_start)/double(CLK_TCK);
 * \endverbatim
 *
 */

struct Clock {
   /**
    * Initialize system clock.  Argument must be in the "clock_t" format
    * which is a standard POSIX struct provided on most systems in the
    * <sys/times.h> include file. On Microsoft systems, it is provided in
    * <time.h>.
    */
   static void
   initialize(
      clock_t& clock)
   {
#ifdef SAMRAI_HAVE_SYS_TIMES_H
      clock = times(&s_tms_buffer);
#endif
   }

   /**
    * Initialize system clock, where clock is in double format.
    */
   static void
   initialize(
      double& clock)
   {
      clock = 0.0;
   }

   /**
    * Timestamp user, system, and walltime clocks.  Wallclock argument is in
    * double format since it will access wallclock times from
    * SAMRAI_MPI::Wtime() function.
    */
   static void
   timestamp(
      clock_t& user,
      clock_t& sys,
      double& wall)
   {
#ifdef SAMRAI_HAVE_SYS_TIMES_H
      s_null_clock_t = times(&s_tms_buffer);
      wall = SAMRAI_MPI::Wtime();
      sys = s_tms_buffer.tms_stime;
      user = s_tms_buffer.tms_utime;
#endif
   }

   /**
    * Returns clock cycle for the system.
    */
   static double
   getClockCycle()
   {
#ifdef _POSIX_VERSION
      double clock_cycle = double(sysconf(_SC_CLK_TCK));
#else
      double clock_cycle = 1.0;
#endif
      return clock_cycle;
   }

private:
#ifdef SAMRAI_HAVE_SYS_TIMES_H
   static struct tms s_tms_buffer;
#endif
   static clock_t s_null_clock_t;
};

}
}

#endif
