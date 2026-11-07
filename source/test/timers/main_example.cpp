/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Example program to demonstrate timers.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <memory>

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      tbox::PIO::logAllNodes("Timer.log");

      std::string input_filename;

      if ((argc != 2)) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                    << "  options:\n"
                    << "  none at this time"
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
         return -1;
      }

      input_filename = argv[1];

      /*
       * Create an input database "input_db" and parse input file (specified
       * on the command line.
       */
      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Create timer manager, reading timer list from input.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Add a timer.  All timers are maintained by the tbox::TimerManager
       * and each timer is accessed by its name.  A static pointer is used
       * to avoid the name lookup cost each time it is called (its only
       * looked up the first time).
       */
      std::string name = "main::test";
      std::shared_ptr<tbox::Timer> timer(
         tbox::TimerManager::getManager()->getTimer(name));

      /*
       * Start timer.  If the timer name was not specified in the input
       * file, the timer start() and stop() functions will simply drop through
       * without recording the time.
       */
      timer->start();

      /*
       * Sleep for 10 sec.
       */
      sleep(10);

      /*
       * Stop the timer
       */
      timer->stop();

      /*
       * Print results to log file.
       */
      tbox::TimerManager::getManager()->print(tbox::plog);

   }

   /*
    * We're done.  Shut down application ...
    */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return 0;
}
