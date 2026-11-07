/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program to demonstrate/test timers.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdlib.h>

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Database.h"
#include "Foo.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>

// Simple code to check timer overhead
// Not part of test, but kept here in
// case it is useful for something
#undef CHECK_PTIMER
//#define CHECK_PTIMER

#ifdef CHECK_PTIMER
#include <sys/types.h>
#include <sys/time.h>

#include <memory>

class PTimer
{
public:
   PTimer():
      d_accesses(0),
      d_total_time(0.0),
      d_last_start_time(0.0) {
   }

   ~PTimer() {
   }

   void start()
   {
      ++d_accesses;

      static struct timeval tp;
      gettimeofday(&tp, (struct timezone *)0);
      d_last_start_time = static_cast<double>(tp.tv_sec)
         + (1.0e-6) * (tp.tv_usec);
   }

   void stop()
   {
      static struct timeval tp;
      gettimeofday(&tp, (struct timezone *)0);
      d_total_time += static_cast<double>(tp.tv_sec)
         + (1.0e-6) * (tp.tv_usec) - d_last_start_time;
   }

   int getNumAccesses() const {
      return d_accesses;
   }
   double getTotalTime() const {
      return d_total_time;
   }

private:
   PTimer(
      const PTimer&);
   PTimer&
   operator = (
      const PTimer&);

   int d_accesses;
   double d_total_time;
   double d_last_start_time;
};
#endif

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("Timer.log");

      std::string input_filename;
      std::string restart_dirname;
      int restore_num = 0;

      bool is_from_restart = false;

      if ((argc != 2) && (argc != 4)) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                    << "<restart dir> <restore number> [options]\n"
                    << "  options:\n"
                    << "  none at this time"
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
         return -1;
      } else {
         input_filename = argv[1];
         if (argc == 4) {
            restart_dirname = argv[2];
            restore_num = atoi(argv[3]);

            is_from_restart = true;
         }
      }

#ifndef HAVE_HDF5
      is_from_restart = false;
#endif

      int i;

      /*
       * Create an input database "input_db" and parse input file (specified
       * on the command line.
       */
      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(
         input_filename, input_db);

      /*
       * Retrieve "Main" section of the input database.  Read ntimes,
       * which is the number of times the functions are called, and
       * depth, which is the depth of the exclusive timer tree.
       */
      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      int ntimes = 1;
      if (main_db->keyExists("ntimes")) {
         ntimes = main_db->getInteger("ntimes");
      }

      int exclusive_tree_depth = 1;
      if (main_db->keyExists("exclusive_tree_depth")) {
         exclusive_tree_depth = main_db->getInteger("exclusive_tree_depth");
      }

      /*
       * Open the restart file and read information from file into the
       * restart database.
       */
      if (is_from_restart) {
         tbox::RestartManager::getManager()->
         openRestartFile(restart_dirname,
            restore_num,
            mpi.getSize());
      }

      std::shared_ptr<tbox::Database> restart_db(
         tbox::RestartManager::getManager()->getRootDatabase());
      NULL_USE(restart_db);

      /*
       * Create timer manager, reading timer list from input.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Add a timer "manually" (that is, not thru the input file).
       */
      std::shared_ptr<tbox::Timer> timer(
         tbox::TimerManager::getManager()->getTimer("apps::main::main"));
      timer->start();

      /*
       * We no longer need the restart file so lets close it.
       */
      tbox::RestartManager::getManager()->closeRestartFile();

      /*
       * Class Foo contains the functions we want to call.
       */
      Foo* foo = new Foo();

      /*
       * Check time to call function with timer name that is NOT
       * registered.  That is, time a NULL timer call.
       */
      std::shared_ptr<tbox::Timer> timer_off(
         tbox::TimerManager::getManager()->getTimer("apps::main::timer_off"));
      timer_off->start();
      for (i = 0; i < ntimes; ++i) {
         foo->timerOff();
      }
      timer_off->stop();

      /*
       * Check time to call function with timer name that IS
       * registered.
       */
      std::shared_ptr<tbox::Timer> timer_on(
         tbox::TimerManager::getManager()->getTimer("apps::main::timer_on"));
      std::shared_ptr<tbox::Timer> dummy_timer(
         tbox::TimerManager::getManager()->getTimer("apps::Foo::timerOn()"));
      NULL_USE(dummy_timer);
      timer_on->start();
      for (i = 0; i < ntimes; ++i) {
         foo->timerOn();
      }
      timer_on->stop();

      /*
       * Time to call tree-based set of exclusive timers. i.e.
       * Foo->zero() calls Foo->one(), which calls Foo->two(), ...
       * and so forth until we reach specified "exclusive_tree_depth.
       */
      std::shared_ptr<tbox::Timer> timer_excl(
         tbox::TimerManager::getManager()->getTimer("apps::main::exclusive_timer"));
      timer_excl->start();
      for (i = 0; i < ntimes; ++i) {
         foo->zero(exclusive_tree_depth);
      }
      timer_excl->stop();
      timer->stop();

      double eff = timer->computeLoadBalanceEfficiency();
      tbox::pout << "Load Balance eff: " << eff << "%" << std::endl;

//   tbox::TimerManager::getManager()->print(tbox::plog);

      tbox::TimerManager::getManager()->resetAllTimers();

      timer->start();

      /*
       * Check time to call function with timer name that is NOT
       * registered.  That is, time a NULL timer call.
       */
      timer_off->start();
      for (i = 0; i < ntimes; ++i) {
         foo->timerOff();
      }
      timer_off->stop();

      /*
       * Check time to call function with timer name that IS
       * registered.
       */
      timer_on->start();
      for (i = 0; i < ntimes; ++i) {
         foo->timerOn();
      }
      timer_on->stop();

      /*
       * Time to call tree-based set of exclusive timers. i.e.
       * Foo->zero() calls Foo->one(), which calls Foo->two(), ...
       * and so forth until we reach specified "exclusive_tree_depth.
       */
      timer_excl->start();
      for (i = 0; i < ntimes; ++i) {
         foo->zero(exclusive_tree_depth);
      }
      timer_excl->stop();
      timer->stop();

      /*
       * Check if we can allocate a large number of timers
       */
      const int max_timers = 575;

      std::shared_ptr<tbox::Timer> timers[max_timers];
      for (int timer_number = 0; timer_number < max_timers; ++timer_number) {

         std::string timer_name = "testcount-" + tbox::Utilities::intToString(
               timer_number,
               4);

         timers[timer_number] = tbox::TimerManager::getManager()->
            getTimer(timer_name);

         if (!timers[timer_number]) {
            TBOX_ERROR("Failed to allocate timer, max was " << timer_number);
         }
      }

      for (int timer_number = 0; timer_number < max_timers; ++timer_number) {
         timers[timer_number].reset();
      }

#ifdef CHECK_PTIMER

      const int nsleepsec = 1;
      const int testit = 3;

      std::shared_ptr<tbox::Timer> tarray1[testit];
      std::shared_ptr<tbox::Timer> tarray2[testit];

      for (int tcnt = 0; tcnt < testit; ++tcnt) {

         int mfactor = 1;
         for (int i = 0; i < tcnt; ++i) {
            mfactor *= 10;
         }
         std::string suffix(tbox::Utilities::intToString(mfactor, testit));

         std::string t1name("ttest1-" + suffix);
         tarray1[tcnt] = tbox::TimerManager::getManager()->getTimer(t1name,
               true);

         std::string t2name("ttest2-" + suffix);
         tarray2[tcnt] = tbox::TimerManager::getManager()->getTimer(t2name,
               true);
      }

      tbox::pout << "\n\nEstimate SAMRAI Timer overhead..." << std::endl;

      for (int tcnt = 0; tcnt < testit; ++tcnt) {

         int ntest = 1;
         for (int i = 0; i < tcnt; ++i) {
            ntest *= 10;
         }

         for (int i = 0; i < ntest; ++i) {
            tarray1[tcnt]->start();
            tarray2[tcnt]->start();
            sleep(nsleepsec);
            tarray2[tcnt]->stop();
            tarray1[tcnt]->stop();
         }

         double access_time1 = tarray1[tcnt]->getTotalWallclockTime();
         double access_time2 = tarray2[tcnt]->getTotalWallclockTime();
         tbox::pout.precision(16);
         tbox::pout << "\ntcnt, ntest, accesses = "
                    << tcnt << " , " << ntest << " , "
                    << tarray1[tcnt]->getNumberAccesses() << std::endl;
         tbox::pout << "access_time1 = " << access_time1 << std::endl;
         tbox::pout << "access_time2 = " << access_time2 << std::endl;
         double access_time12 = access_time1 - access_time2;
         double cost12 = access_time12 / static_cast<double>(ntest);
         tbox::pout << "Access time 12 overhead for each Timer: " << cost12
                    << " sec" << std::endl;
      }

      tbox::pout << "\n\nEstimate PTimer overhead..." << std::endl;

      PTimer Ptarray1[testit];
      PTimer Ptarray2[testit];

      for (int tcnt = 0; tcnt < testit; ++tcnt) {

         int ntest = 1;
         for (int i = 0; i < tcnt; ++i) {
            ntest *= 10;
         }

         for (int i = 0; i < ntest; ++i) {
            Ptarray1[tcnt].start();
            Ptarray2[tcnt].start();
            sleep(nsleepsec);
            Ptarray2[tcnt].stop();
            Ptarray1[tcnt].stop();
         }

         double access_time1 = Ptarray1[tcnt].getTotalTime();
         double access_time2 = Ptarray2[tcnt].getTotalTime();
         tbox::pout.precision(16);
         tbox::pout << "\ntcnt, ntest, accesses = "
                    << tcnt << " , " << ntest << " , "
                    << Ptarray1[tcnt].getNumAccesses() << std::endl;
         tbox::pout << "access_time1 = " << access_time1 << std::endl;
         tbox::pout << "access_time2 = " << access_time2 << std::endl;
         double access_time12 = access_time1 - access_time2;
         double cost12 = access_time12 / static_cast<double>(ntest);
         tbox::pout << "Access time 12 overhead for each Timer: " << cost12
                    << " sec" << std::endl;
      }

#endif

      tbox::TimerManager::getManager()->print(tbox::plog);

      /*
       * We're done.  Write the restart file.
       */
      std::string restart_write_dirname = "restart";
#ifdef HAVE_HDF5
      int timestep = 0;
      tbox::RestartManager::getManager()->writeRestartFile(
         restart_write_dirname,
         timestep);
#endif

      delete foo;

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  timertest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}
