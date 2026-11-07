/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class to manage different timer objects used throughout the
 *                library.
 *
 ************************************************************************/

#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>

#ifndef ENABLE_SAMRAI_TIMERS
#ifdef __INTEL_COMPILER
#pragma warning (disable:869)
#endif
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

TimerManager * TimerManager::s_timer_manager_instance = 0;

int TimerManager::s_main_timer_identifier = -1;
int TimerManager::s_inactive_timer_identifier = -9999;

StartupShutdownManager::Handler
TimerManager::s_finalize_handler(
   0,
   0,
   0,
   TimerManager::finalizeCallback,
   StartupShutdownManager::priorityTimerManger);

/*
 *************************************************************************
 *
 * Static timer manager member functions.
 *
 *************************************************************************
 */

void
TimerManager::createManager(
   const std::shared_ptr<Database>& input_db)
{
   if (!s_timer_manager_instance) {
      s_timer_manager_instance = new TimerManager(input_db);

      /*
       * Compute the overheads
       */
      s_timer_manager_instance->computeOverheadConstants();

      s_timer_manager_instance->d_main_timer->start();
   } else {
      /*
       * Manager already exists so simply apply the rules
       * from the input file to the existing timers.
       */
      s_timer_manager_instance->getFromInput(input_db);
      s_timer_manager_instance->activateExistingTimers();
   }
}

TimerManager *
TimerManager::getManager()
{
   if (!s_timer_manager_instance) {
      TBOX_WARNING("TimerManager::getManager() is called before\n"
         << "createManager().  Creating the timer manager\n"
         << "(without using input database.)\n");
      createManager(std::shared_ptr<Database>());
   }

   return s_timer_manager_instance;
}

void
TimerManager::finalizeCallback()
{
   if (s_timer_manager_instance) {
      delete s_timer_manager_instance;
      s_timer_manager_instance = 0;
   }
}

void
TimerManager::registerSingletonSubclassInstance(
   TimerManager* subclass_instance)
{
   if (!s_timer_manager_instance) {
      s_timer_manager_instance = subclass_instance;
   } else {
      TBOX_ERROR("TimerManager internal error...\n"
         << "Attemptng to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Protected timer manager constructor and destructor.
 *
 *************************************************************************
 */

TimerManager::TimerManager(
   const std::shared_ptr<Database>& input_db)
#ifdef ENABLE_SAMRAI_TIMERS
   : d_timer_active_access_time(-9999.0),
   d_timer_inactive_access_time(-9999.0),
   d_main_timer(new Timer("TOTAL RUN TIME")),
   d_length_package_names(0),
   d_length_class_names(0),
   d_length_class_method_names(0),
   d_print_threshold(0.25),
   d_print_exclusive(false),
   d_print_total(true),
   d_print_processor(true),
   d_print_max(false),
   d_print_summed(false),
   d_print_user(false),
   d_print_sys(false),
   d_print_wall(true),
   d_print_percentage(true),
   d_print_concurrent(false),
   d_print_timer_overhead(false)
#endif
{
   /*
    * Create a timer that measures overall solution time
    */
#ifdef ENABLE_SAMRAI_TIMERS
   getFromInput(input_db);
#else
   NULL_USE(input_db);
#endif
}

TimerManager::~TimerManager()
{
#ifdef ENABLE_SAMRAI_TIMERS
   d_main_timer->stop();
   d_main_timer.reset();

   d_timers.clear();
   d_inactive_timers.clear();

   d_exclusive_timer_stack.clear();

   d_package_names.clear();
   d_class_names.clear();
   d_class_method_names.clear();
#endif
}

/*
 *************************************************************************
 *
 * Utility functions for creating timers, adding them to the manager
 * database, and checking whether a particular timer exists.
 *
 *    o checkTimerExistsInArray is private.  It returns true if a timer
 *      matching the name string exists in the array and returns false
 *      otherwise.  If such a timer exists, the pointer is set to that
 *      timer; otherwise the pointer is null.
 *
 *    o getTimer returns a timer with the given name.  It will be active
 *      if its name appears in the input file, or if it is added with
 *      a `true' argument.  Otherwise the timer will be inactive.
 *
 *    o checkTimerExists returns true if a timer matching the name
 *      string exists and false otherwise.  If such a timer exists,
 *      the pointer is set to that timer; otherwise the pointer is null.
 *
 *************************************************************************
 */

bool
TimerManager::checkTimerExistsInArray(
   std::shared_ptr<Timer>& timer,
   const std::string& name,
   const std::vector<std::shared_ptr<Timer> >& timer_array) const
{

   bool timer_found = false;
#ifdef ENABLE_SAMRAI_TIMERS

   timer.reset();
   if (!name.empty()) {
      for (size_t i = 0; i < timer_array.size(); ++i) {
         if (timer_array[i]->getName() == name) {
            timer_found = true;
            timer = timer_array[i];
            break;
         }
      }
   }
#else
   NULL_USE(timer);
   NULL_USE(name);
   NULL_USE(timer_array);
#endif
   return timer_found;
}

void
TimerManager::activateExistingTimers() {
#ifdef ENABLE_SAMRAI_TIMERS
   std::vector<std::shared_ptr<Timer> >::iterator it =
      d_inactive_timers.begin();
   while (it != d_inactive_timers.end()) {
      bool timer_active = checkTimerInNameLists((*it)->getName());
      if (timer_active) {
         (*it)->setActive(true);
         d_timers.push_back((*it));
         it = d_inactive_timers.erase(it);
      } else {
         ++it;
      }
   }
#endif
}

std::shared_ptr<Timer>
TimerManager::getTimer(
   const std::string& name,
   bool ignore_timer_input)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::shared_ptr<Timer> timer;

   TBOX_ASSERT(!name.empty());
   bool timer_active = true;
   if (!ignore_timer_input) {
      /*
       * Check if name is in either the d_package_names, d_class_names,
       * or d_class_method_names lists.  Add it if it is.
       */
      timer_active = checkTimerInNameLists(name);
   }

   /*
    * Add the timer to the appropriate array, if necessary.
    */
   if (timer_active) {
      if (!checkTimerExistsInArray(timer,
             name,
             d_timers)) {
         if (d_timers.size() == d_timers.capacity()) {
            d_timers.reserve(d_timers.size()
               + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
         }
         timer.reset(new Timer(name));
         d_timers.push_back(timer);
      }
   } else {
      if (!checkTimerExistsInArray(timer,
             name,
             d_inactive_timers)) {
         if (d_inactive_timers.size() == d_inactive_timers.capacity()) {
            d_inactive_timers.reserve(d_inactive_timers.size()
               + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
         }
         timer.reset(new Timer(name));
         timer->setActive(false);
         d_inactive_timers.push_back(timer);
      }
   }
   return timer;

#else
   // since timers aren't active - and we need to still provide
   // pseudo-timer functionality (i.e., a valid timer), we'll
   // create one on the fly, but not track it.
   NULL_USE(ignore_timer_input);
   std::shared_ptr<Timer> timer(new Timer(name));
   timer->setActive(false);
   return timer;

#endif
}

bool
TimerManager::checkTimerExists(
   std::shared_ptr<Timer>& timer,
   const std::string& name) const
{
#ifdef ENABLE_SAMRAI_TIMERS
   TBOX_ASSERT(!name.empty());

   bool timer_found = checkTimerExistsInArray(timer,
         name,
         d_timers);

   if (!timer_found) {
      timer_found = checkTimerExistsInArray(timer,
            name,
            d_inactive_timers);
   }

   return timer_found;

#else
   NULL_USE(timer);
   NULL_USE(name);

   return false;

#endif
}

/*
 *************************************************************************
 *
 * Utility functions to check whether timer is running and to reset
 * all active timers.
 *
 *************************************************************************
 */

bool
TimerManager::checkTimerRunning(
   const std::string& name) const
{
   bool is_running = false;
#ifdef ENABLE_SAMRAI_TIMERS

   TBOX_ASSERT(!name.empty());
   std::shared_ptr<Timer> timer;
   if (checkTimerExistsInArray(timer, name, d_timers)) {
      is_running = timer->isRunning();
   }
#else
   NULL_USE(name);
#endif
   return is_running;
}

void
TimerManager::resetAllTimers()
{
#ifdef ENABLE_SAMRAI_TIMERS
   d_main_timer->stop();
   d_main_timer->reset();
   d_main_timer->start();

   for (size_t i = 0; i < d_timers.size(); ++i) {
      d_timers[i]->reset();
   }
   for (size_t j = 0; j < d_inactive_timers.size(); ++j) {
      d_inactive_timers[j]->reset();
   }
#endif
}

/*
 *************************************************************************
 *
 * Protected start and stop routines for exclusive timer management.
 *
 *************************************************************************
 */

void
TimerManager::startTime(
   Timer* timer)
{
#ifdef ENABLE_SAMRAI_TIMERS
   TBOX_ASSERT(timer != 0);

   if (timer->isActive()) {
   }

   if (d_print_exclusive) {
      if (!d_exclusive_timer_stack.empty()) {
         ((Timer *)d_exclusive_timer_stack.front())->stopExclusive();
      }
      Timer* stack_timer = timer;
      d_exclusive_timer_stack.push_front(stack_timer);
      stack_timer->startExclusive();
   }

   if (d_print_concurrent) {
      for (size_t i = 0; i < d_timers.size(); ++i) {
         if ((d_timers[i].get() != timer) && d_timers[i]->isRunning()) {
            d_timers[i]->addConcurrentTimer(*d_timers[i]);
         }
      }
   }
#else
   NULL_USE(timer);
#endif
}

void
TimerManager::stopTime(
   Timer* timer)
{
#ifdef ENABLE_SAMRAI_TIMERS
   TBOX_ASSERT(timer != 0);

   if (d_print_exclusive) {
      timer->stopExclusive();
      if (!d_exclusive_timer_stack.empty()) {
         d_exclusive_timer_stack.pop_front();
         if (!d_exclusive_timer_stack.empty()) {
            ((Timer *)d_exclusive_timer_stack.front())->startExclusive();
         }
      }
   }
#else
   NULL_USE(timer);
#endif
}

/*
 *************************************************************************
 *
 * Parser to see if Timer has been registered or not.
 *
 *************************************************************************
 */

bool
TimerManager::checkTimerInNameLists(
   const std::string& copy)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string name = copy;
   /*
    * string::size_type is generally an int, but it may depend on vendor's
    * implementation of string class.  Just to be safe, we use the definition
    * from the string class (string::size_type).
    */
   std::string::size_type string_length, list_entry_length, position;

   /*
    * Specification of whether we will use the timer after comparing to
    * package, class, and class::method lists.
    */
   bool will_use_timer = false;

   /*
    * Start by evaluating the Timer's package.  If the list of packages
    * has length zero, we can skip this step. We can also skip this step
    * if the timer does not have two "::" in its name.
    *
    * The manager will assume the following with regard to the specified timer
    * name:  The name has...
    *  1) Two "::" - it is of the form Package::Class::method.  A user can
    *     turn on the timer by specifying the timers Package, Class, or
    *     Class::method combination in the input file.  This is the most
    *     versatile.
    *  2) One "::" - it has the form Class::method.  The timer is assumed not
    *     to be in a package.  A user can turn on the timer by specifying its
    *     Class or class::method combination.
    *  3) No "::" - it has the form Class.  The timer is can be turned on by
    *     entering the class name in the input file.
    */
   if (d_length_package_names > 0) {
      /*
       * See how many instances of "::" there are in the name.  If there
       * are at least two, the timer has a package which might be in the
       * package list.
       */
      int occurrences = 0;
      position = name.find("::");
      if (position < name.size()) {
         ++occurrences;
         std::string substring = name.substr(position + 2);
         position = substring.find("::");
         if (position < substring.size()) {
            ++occurrences;
         }
      }

      if (occurrences >= 2) {
         /*
          * The timer may be in the package list.  Parse to get its package
          * name and compare to the list entries.
          */
         position = name.find("::");
         std::string package = name.substr(0, position);

         /*
          * Iterate through package list and see if the timer's package
          * name is there.
          */
         bool package_exists = false;
         string_length = package.size();
         for (std::list<std::string>::iterator i = d_package_names.begin();
              i != d_package_names.end(); ++i) {
            list_entry_length = i->size();
            if (string_length == list_entry_length) {
               package_exists = (*i == package);
            }
            if (package_exists) {
               break;
            }
         }
         will_use_timer = package_exists;
      }
   }

   if (!will_use_timer) {

      /*
       * The timer's package has already been compared to the package list,
       * so we can parse it off to ease further evaluations of the name.
       */
      int occurrences = 0;
      position = name.find("::");
      if (position < name.size()) {
         ++occurrences;
         std::string substring = name.substr(position + 2);
         position = substring.find("::");
         if (position < substring.size()) {
            ++occurrences;
         }
      }
      if (occurrences >= 2) {
         position = name.find("::");
         name = name.substr(position + 2);
      }

      /*
       * See if Timer's class is in d_class_names.  If the list of classes
       * has length zero, we can skip this step.
       */
      if (d_length_class_names > 0) {

         /*
          * If a "::" exists in the timer name, parse off the part before
          * before it.  This is the class name.  If no "::" exists, assume
          * the timer name is a class name (don't do any parsing) and compare
          * it directly to the class list.
          */
         position = name.find("::");
         std::string class_name;
         if (position < name.size()) {
            class_name = name.substr(0, position);
         } else {
            class_name = name;
         }

         /*
          * Is class name dimensional?
          */
         string_length = class_name.size();
         std::string dim = class_name.substr(string_length - 1, 1);
         bool is_dimensional = false;
         std::string nondim_class_name;
         if (dim == "1" || dim == "2" || dim == "3") {
            is_dimensional = true;
            nondim_class_name = class_name.substr(0, string_length - 1);
         } else {
            nondim_class_name = class_name;
         }

         /*
          * See if class name is in class list.  Accelerate this by comparing
          * the length of the entries.  Do non-dimensional comparison first.
          */
         string_length = nondim_class_name.size();
         bool class_exists = false;
         for (std::list<std::string>::iterator i = d_class_names.begin();
              i != d_class_names.end(); ++i) {
            list_entry_length = i->size();
            if (string_length == list_entry_length) {
               class_exists = (*i == nondim_class_name);
            }
            if (class_exists) {
               break;
            }
         }

         /*
          * Now do dimensional comparison.
          */
         string_length = class_name.size();
         if (is_dimensional && !class_exists) {
            for (std::list<std::string>::iterator i = d_class_names.begin();
                 i != d_class_names.end(); ++i) {
               list_entry_length = i->size();
               if (string_length == list_entry_length) {
                  class_exists = (*i == class_name);
               }
               if (class_exists) {
                  break;
               }
            }
         }
         will_use_timer = class_exists;

      }

      /*
       * See if Timer's class::method name is in d_class_method_names list.
       *
       * If the list of class_method_names has length zero, we can skip
       * this step.  Also, if no "::" exists in the timer name then it
       * cannot be in the timer list so lets avoid the comparison.
       */
      position = name.find("::");
      occurrences = 0;
      if (position < name.size()) {
         occurrences = 1;
      }

      if (!will_use_timer && d_length_class_method_names > 0 &&
          occurrences > 0) {

         /*
          * Parse name before "::" - this is the class name.
          */
         position = name.find("::");
         std::string class_name = name.substr(0, position);

         /*
          * Is class name dimensional?
          */
         string_length = class_name.size();
         std::string dim = class_name.substr(string_length - 1, 1);
         bool is_dimensional = false;
         std::string nondim_name;
         if (dim == "1" || dim == "2" || dim == "3") {
            is_dimensional = true;
            std::string nondim_class_name = class_name.substr(0,
                  string_length - 1);
            std::string method_name = name.substr(position);
            nondim_name = nondim_class_name;
            nondim_name += method_name;

         } else {
            nondim_name = name;
         }

         /*
          * See if name is in class_method_names list.  Accelerate this by
          * comparing the length of the entries.  Do non-dimensional
          * comparison first.
          */
         bool class_method_exists = false;
         string_length = nondim_name.size();
         for (std::list<std::string>::iterator i = d_class_method_names.begin();
              i != d_class_method_names.end(); ++i) {
            list_entry_length = i->size();
            if (string_length == list_entry_length) {
               class_method_exists = (*i == nondim_name);
            }
            if (class_method_exists) {
               break;
            }
         }

         /*
          * Now do dimensional comparison.
          */
         if (is_dimensional && !class_method_exists) {
            string_length = name.size();
            for (std::list<std::string>::iterator i = d_class_method_names.begin();
                 i != d_class_method_names.end(); ++i) {
               list_entry_length = i->size();
               if (string_length == list_entry_length) {
                  class_method_exists = (*i == name);
               }
               if (class_method_exists) {
                  break;
               }
            }
         }

         will_use_timer = class_method_exists;

      }

   }

   return will_use_timer;

#else
   NULL_USE(copy);

   return false;

#endif
}

/*
 *************************************************************************
 *
 * Print timer information from each processor.
 *
 *************************************************************************
 */

void
TimerManager::print(
   std::ostream& os)
{
#ifdef ENABLE_SAMRAI_TIMERS
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   /*
    * There are 18 possible timer values that users may wish to look at.
    * (i.e. User/sys/wallclock time, Total or Exclusive, for individual
    * processor, max across all procs, or summed across all procs).
    * This method builds an array that holds these values and outputs
    * them in column format.
    */

   /*
    * First, stop the main_timer so we have an accurate measure of
    * overall run time so far.
    */
   d_main_timer->stop();

   /*
    * If we are doing max or sum operations, make sure timers are
    * consistent across processors.
    */
   if (d_print_summed || d_print_max) {
      checkConsistencyAcrossProcessors();
   }

   /*
    * Invoke arrays used to hold timer names, timer_values, and
    * max_processor_ids (i.e. processor holding the maximum time).
    */
   double(*timer_values)[18] = new double[d_timers.size() + 1][18];
   int(*max_processor_id)[2] = new int[d_timers.size() + 1][2];
   std::vector<std::string> timer_names(static_cast<int>(d_timers.size()) + 1);

   /*
    * Fill in timer_values and timer_names arrays, based on values of
    * d_print_total, d_print_exclusive,
    * d_print_user, d_print_sys, d_print_wallclock,
    * d_print_processor, d_print_summed, d_print_max.
    */
   buildTimerArrays(timer_values, max_processor_id, timer_names);

   /*
    * Now that we have built up the array, select output options.
    *
    * 1) If user has requested any two (or more) of {user, system, and
    *    walltime} then use the format:
    *    [Exclusive,Total]
    *    [Processor,Summed,Max]
    *       Name     User    System   Wall  (Max Processor)
    *
    * 2) If user chose just one of {User, system, or walltime}, then use the
    *    format:
    *    [Exclusive,Total]
    *       Name     Processor   Summed   Max  Max Processor
    *
    * 3) If user chose just one of {user, system, or walltime} and just one
    *    of {processor, summed, max} then use the format:
    *       Name     [Processor,Summed,Max]  Total   Exclusive
    *
    * If the user wants overhead stats, print those as:
    *   Timer Overhead:
    *     Timer Name    number calls     Estimated overhead
    *
    * If they want output of a concurrent timer tree, print this as:
    *   Concurrent Tree:
    *     Timer Name    names of timers called by it.
    */

   /*
    * Determine which case we are doing - #1, #2, or #3
    * #1 case - user requested any two (or more) of [user,system,wallclock]
    * #2 case - user requested one of [user,system,wallclock]
    * #3 case - user requested one of [user,system,wallclock] and one of
    *           [processor, summed, max]
    */
   bool case1 = false;
   bool case2 = false;
   bool case3 = false;
   if ((d_print_user && d_print_sys) ||
       (d_print_sys && d_print_wall) ||
       (d_print_wall && d_print_user)) {
      case1 = true;
   } else {
      case2 = true;
      if ((d_print_processor && !d_print_summed && !d_print_max) ||
          (!d_print_processor && d_print_summed && !d_print_max) ||
          (!d_print_processor && !d_print_summed && d_print_max)) {
         case2 = false;
         case3 = true;
      }
   }

   std::string table_title;
   std::vector<std::string> column_titles(4);
   int column_ids[3] = { 0, 0, 0 };
   int j, k;

   /*
    * Now print out the data
    */
   if (case1) {

      column_titles[0] = "";
      column_titles[1] = "";
      column_titles[2] = "";
      column_titles[3] = "Proc";

      if (d_print_user) {
         column_titles[0] = "User Time";
      } else {
         column_titles[0] = "";
      }
      if (d_print_sys) {
         column_titles[1] = "Sys Time";
      } else {
         column_titles[1] = "";
      }
      if (d_print_wall) {
         column_titles[2] = "Wall Time";
      } else {
         column_titles[2] = "";
      }

      for (k = 0; k < 2; ++k) {

         if ((k == 0 && d_print_exclusive) ||
             (k == 1 && d_print_total)) {

            for (j = 0; j < 3; ++j) {

               if ((j == 0 && d_print_processor) ||
                   (j == 1 && d_print_summed) ||
                   (j == 2 && d_print_max)) {

                  if (j == 0) {
#ifndef LACKS_SSTREAM
                     std::ostringstream out;
#endif
                     if (k == 0) {
#ifndef LACKS_SSTREAM
                        out << "EXCLUSIVE TIME \nPROCESSOR:"
                            << mpi.getRank();
                        table_title = out.str();
#else
                        table_title = "EXCLUSIVE TIME \nPROCESSOR:";
#endif
                        column_ids[0] = 0;
                        column_ids[1] = 1;
                        column_ids[2] = 2;
                     } else if (k == 1) {
#ifndef LACKS_SSTREAM
                        out << "TOTAL TIME \nPROCESSOR:"
                            << mpi.getRank();
                        table_title = out.str();
#else
                        table_title = "TOTAL TIME \nPROCESSOR:";
#endif
                        column_ids[0] = 9;
                        column_ids[1] = 10;
                        column_ids[2] = 11;
                     }
                     printTable(table_title,
                        column_titles,
                        timer_names,
                        column_ids,
                        timer_values,
                        os);
                  } else if (j == 1) {

                     if (k == 0) {
                        table_title =
                           "EXCLUSIVE TIME \nSUMMED ACROSS ALL PROCESSORS";
                        column_ids[0] = 3;
                        column_ids[1] = 4;
                        column_ids[2] = 5;
                     } else if (k == 1) {
                        table_title =
                           "TOTAL TIME \nSUMMED ACROSS ALL PROCESSORS:";
                        column_ids[0] = 12;
                        column_ids[1] = 13;
                        column_ids[2] = 14;
                     }
                     printTable(table_title,
                        column_titles,
                        timer_names,
                        column_ids,
                        timer_values,
                        os);

                  } else if (j == 2) {

                     int max_array_id = 0; // identifies which of the two
                                           // max_processor_id values to print
                     if (k == 0) {
                        table_title =
                           "EXCLUSIVE TIME \nMAX ACROSS ALL PROCESSORS";
                        column_ids[0] = 6;
                        column_ids[1] = 7;
                        column_ids[2] = 8;
                        max_array_id = 0;
                     } else if (k == 1) {
                        table_title =
                           "TOTAL TIME \nMAX ACROSS ALL PROCESSORS";
                        column_ids[0] = 15;
                        column_ids[1] = 16;
                        column_ids[2] = 17;
                        max_array_id = 1;
                     }
#if 0
                     printTable(table_title,
                        column_titles,
                        timer_names,
                        &max_processor_id[0][max_array_id],
                        column_ids,
                        timer_values,
                        os);
#else
                     printTable(table_title,
                        column_titles,
                        timer_names,
                        max_processor_id,
                        max_array_id,
                        column_ids,
                        timer_values,
                        os);
#endif
                  }

               } // if j
            } // for j
         } // if k
      } // for k
   } // if case 1

   if (case2) {

      for (k = 0; k < 2; ++k) {

         if ((k == 0 && d_print_exclusive) ||
             (k == 1 && d_print_total)) {

            int max_array_id = 0;
            std::string table_title_line_1;
            std::string table_title_line_2;
            if (k == 0) {
               table_title_line_1 = "EXCLUSVE \n";
               max_array_id = 0;
            } else if (k == 1) {
               table_title_line_1 = "TOTAL \n";
               max_array_id = 1;
            }
            if (d_print_user) {
               table_title_line_2 = "USER TIME";
            } else if (d_print_sys) {
               table_title_line_2 = "SYSTEM TIME";
            } else if (d_print_wall) {
               table_title_line_2 = "WALLCLOCK TIME";
            }
            table_title = table_title_line_1;
            table_title += table_title_line_2;

            column_titles[0] = "";
            column_titles[1] = "";
            column_titles[2] = "";
            column_titles[3] = "";
            if (d_print_processor) {
#ifndef LACKS_SSTREAM
               std::ostringstream out;
               out << "Proc: " << mpi.getRank();
               column_titles[0] = out.str();
#else
               column_titles[0] = "Proc: ";
#endif
            }
            if (d_print_summed) {
               column_titles[1] = "Summed";
            }
            if (d_print_max) {
               column_titles[2] = "Max";
               column_titles[3] = "Proc";
            }

            if (d_print_user) {
               if (k == 0) {
                  column_ids[0] = 0;
                  column_ids[1] = 3;
                  column_ids[2] = 6;
               } else if (k == 1) {
                  column_ids[0] = 9;
                  column_ids[1] = 12;
                  column_ids[2] = 15;
               }
            } else if (d_print_sys) {
               if (k == 0) {
                  column_ids[0] = 1;
                  column_ids[1] = 4;
                  column_ids[2] = 7;
               } else if (k == 1) {
                  column_ids[0] = 10;
                  column_ids[1] = 13;
                  column_ids[2] = 16;
               }
            } else if (d_print_wall) {
               if (k == 0) {
                  column_ids[0] = 2;
                  column_ids[1] = 5;
                  column_ids[2] = 8;
               } else if (k == 1) {
                  column_ids[0] = 11;
                  column_ids[1] = 14;
                  column_ids[2] = 17;
               }
            }
#if 0
            printTable(table_title,
               column_titles,
               timer_names,
               &max_processor_id[0][max_array_id],
               column_ids,
               timer_values,
               os);
#else
            printTable(table_title,
               column_titles,
               timer_names,
               max_processor_id,
               max_array_id,
               column_ids,
               timer_values,
               os);
#endif
         } // if k
      }  // for k
   } // if case2

   if (case3) {

      if (d_print_exclusive && !d_print_total) {
         column_titles[0] = "Exclusive";
         column_titles[1] = "";
      } else if (!d_print_exclusive && d_print_total) {
         column_titles[0] = "";
         column_titles[1] = "Total";
      } else if (d_print_exclusive && d_print_total) {
         column_titles[0] = "Exclusive";
         column_titles[1] = "Total";
      }
      column_titles[3] = "";

      column_ids[2] = 0;
      if (d_print_user) {
         if (d_print_processor) {
#ifndef LACKS_SSTREAM
            std::ostringstream out;
            out << "USER TIME \nPROCESSOR: " << mpi.getRank();
            table_title = out.str();
#else
            table_title = "USER TIME \nPROCESSOR: ";
#endif
            column_ids[0] = 0;
            column_ids[1] = 9;
         } else if (d_print_summed) {
            table_title = "USER TIME \nSUMMED ACROSS ALL PROCESSORS";
            column_ids[0] = 3;
            column_ids[1] = 12;
         } else if (d_print_max) {
            table_title = "USER TIME \nMAX ACROSS ALL PROCESSORS";
            column_ids[0] = 6;
            column_ids[1] = 15;
         }
      } else if (d_print_sys) {
         if (d_print_processor) {
#ifndef LACKS_SSTREAM
            std::ostringstream out;
            out << "SYSTEM TIME \nPROCESSOR: " << mpi.getRank();
            table_title = out.str();
#else
            table_title = "SYSTEM TIME \nPROCESSOR:";
#endif
            column_ids[0] = 1;
            column_ids[1] = 10;
         } else if (d_print_summed) {
            table_title = "SYSTEM TIME \nSUMMED ACROSS ALL PROCESSORS";
            column_ids[0] = 4;
            column_ids[1] = 13;
         } else if (d_print_max) {
            table_title = "SYSTEM TIME \nMAX ACROSS ALL PROCESSORS";
            column_ids[0] = 7;
            column_ids[1] = 16;
         }
      } else if (d_print_wall) {
         if (d_print_processor) {
#ifndef LACKS_SSTREAM
            std::ostringstream out;
            out << "WALLCLOCK TIME \nPROCESSOR: " << mpi.getRank();
            table_title = out.str();
#else
            table_title = "WALLCLOCK TIME \nPROCESSOR: ";
#endif
            column_ids[0] = 2;
            column_ids[1] = 11;
         } else if (d_print_summed) {
            table_title = "WALLCLOCK TIME \nSUMMED ACROSS ALL PROCESSORS";
            column_ids[0] = 5;
            column_ids[1] = 14;
         } else if (d_print_max) {
            table_title = "WALLCLOCK TIME \nMAX ACROSS ALL PROCESSORS";
            column_ids[0] = 8;
            column_ids[1] = 17;
         }
      }
      printTable(table_title,
         column_titles,
         timer_names,
         column_ids,
         timer_values,
         os);
   }

   /*
    * Print overhead stats - number of accesses and estimated cost
    * (estimated cost computed based on the number of accesses and
    * a fixed d_timer_active_access_time value).
    * Store the number of accesses in max_processor_id[0] and the estimated
    * cost in timer_values[0] and use the printTable method.
    */
   if (d_print_timer_overhead) {
      printOverhead(timer_names,
         timer_values,
         os);
   }

   /*
    * Print tree of concurrent timers.
    */
   if (d_print_concurrent) {
      printConcurrent(os);
   }

   delete[] timer_values;
   delete[] max_processor_id;
   /*
    * Lastly, restart the main_timer that we stopped at the beginning of
    * this routine
    */
   d_main_timer->start();
#else
   os << "Timers disabled\n";
#endif
}

void
TimerManager::printTable(
   const std::string& table_title,
   const std::vector<std::string>& column_titles,
   const std::vector<std::string>& timer_names,
   const int column_ids[],
   const double timer_values[][18],
   std::ostream& os)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
   std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
   std::string ascii_line = ascii_line1;
   ascii_line += ascii_line2;

   /*
    * By default, output in C++ is right justified with the setw()
    * option e.g. cout << "[" << setw(5) << 1 << "]" will output
    * [   1].  The line below makes it left justified, so the same line
    * will generate [1   ].  We us left justification because it is
    * more convenient to output columns of tables.
    */
   os.setf(std::ios::left);

   os << ascii_line << table_title << "\n";

   int i;

   /*
    * Determine maximum name length for formatting
    */
   int maxlen = 10;
   for (unsigned int n = 0; n < d_timers.size() + 1; ++n) {
      i = static_cast<int>(timer_names[n].size());
      if (i > maxlen) maxlen = i;
   }

   /*
    * Print table header.  If we are only printing the overall program
    * timer (i.e. d_num_timers = 0) with only d_print_processor,
    * d_print_total, and d_print_wall options being true (which
    * is the default case if the user doesn't add a "TimerManager"
    * section to the input file) then don't bother to print header as
    * it just clutters up the output.  Also, turn off percentages since
    * this doesn't mean anything with just one timer.
    */
   bool default_case = d_timers.size() == 0 && !d_print_exclusive &&
      !d_print_summed && !d_print_max &&
      !d_print_user && !d_print_sys;
   if (default_case) {
      d_print_percentage = false;
   } else {
      os << ascii_line
         << std::setw(maxlen + 3) << "Timer Name" << ' ';
      for (i = 0; i < 3; ++i) {
         if (!column_titles[i].empty()) {
            os << std::setw(15) << column_titles[i].c_str() << "  ";
         }
      }
      os << std::endl;
   }

   /*
    * Organize timers largest to smallest.  Apply this to the LAST NONZERO
    * column entry for the table by forming an ordering array - ordered_list
    * - that orders these values.
    */
   int last_nonzero_column = 0;
   for (i = 0; i < 3; ++i) {
      if (!column_titles[i].empty()) {
         last_nonzero_column = column_ids[i];
      }
   }
   int* ordered_list = new int[d_timers.size() + 1];
   buildOrderedList(timer_values,
      last_nonzero_column,
      ordered_list,
      static_cast<int>(d_timers.size()));

   /*
    * Tack on TOTAL TIME to end of ordered list
    */
   ordered_list[static_cast<int>(d_timers.size())] =
      static_cast<int>(d_timers.size());

   /*
    * Now output the rows of the table.
    */

   for (size_t k = 0; k < d_timers.size() + 1; ++k) {
      int n = ordered_list[k];

      /*
       * Check the print threshold to see if we should print this timer.
       */
      double frac = computePercentageDouble(
            timer_values[n][last_nonzero_column],
            timer_values[d_timers.size()][last_nonzero_column]);

      if (frac > d_print_threshold) {

         os << std::setw(maxlen + 3) << timer_names[n].c_str() << ' ';

         /*
          * Print column values
          */
         for (i = 0; i < 3; ++i) {

            /*
             * Print column values only title is non-null (i.e. not "")
             */
            if (!column_titles[i].empty()) {

               /*
                * Print percentages if requested.
                */
               int j = column_ids[i];

               if (d_print_percentage) {
#ifndef LACKS_SSTREAM
                  int perc = computePercentageInt(timer_values[n][j],
                        timer_values[d_timers.size()][j]);

                  std::ostringstream out;
                  out << timer_values[n][j] << " (" << perc << "%)";
                  os << std::setw(15) << out.str().c_str() << "  ";
#else
                  os << std::setw(15) << timer_values[n][j] << "  ";
#endif
               } else {
                  os << std::setw(15) << timer_values[n][j] << "  ";
               }

            } // if title is non-null

         } // loop over columns

         os << std::endl;

      } // if meets d_print_threshold condition

   } // loop over timers

   delete[] ordered_list;

   os << ascii_line << std::endl;
   os.setf(std::ios::right);
#else
   NULL_USE(table_title);
   NULL_USE(column_titles);
   NULL_USE(timer_names);
   NULL_USE(column_ids);
   NULL_USE(timer_values);
   NULL_USE(os);
#endif
}

void
TimerManager::printTable(
   const std::string& table_title,
   const std::vector<std::string>& column_titles,
   const std::vector<std::string>& timer_names,
   const int max_processor_id[][2],
   const int max_array_id,
   const int column_ids[],
   const double timer_values[][18],
   std::ostream& os)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
   std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
   std::string ascii_line = ascii_line1;
   ascii_line += ascii_line2;

   /*
    * Left-justify all output in this method.
    */
   os.setf(std::ios::left);

   os << ascii_line
      << table_title << "\n"
      << ascii_line;

   int i;

   /*
    * Determine maximum name length for formatting
    */
   int maxlen = 10;
   for (unsigned int n = 0; n < d_timers.size() + 1; ++n) {
      i = static_cast<int>(timer_names[n].size());
      if (i > maxlen) maxlen = i;
   }

   /*
    * Print table header
    */
   os << std::setw(maxlen + 3) << "Timer Name" << ' ';
   for (i = 0; i < 4; ++i) {
      if (!column_titles[i].empty()) {
         os << std::setw(15) << column_titles[i].c_str() << "  ";
      }
   }
   os << std::endl;

   /*
    * Organize timers largest to smallest.  Apply this to the LAST NONZERO
    * column entry for the table by forming an ordering array - ordered_list
    * - that orders these values.
    */
   int last_nonzero_column = 0;
   for (i = 0; i < 3; ++i) {
      if (!column_titles[i].empty()) {
         last_nonzero_column = column_ids[i];
      }
   }
   int* ordered_list = new int[d_timers.size() + 1];
   buildOrderedList(timer_values,
      last_nonzero_column,
      ordered_list,
      static_cast<int>(d_timers.size()));

   /*
    * Tack on TOTAL TIME to end of ordered list
    */
   ordered_list[static_cast<int>(d_timers.size())] =
      static_cast<int>(d_timers.size());

   /*
    * Now output the rows of the table.
    */
   for (size_t j = 0; j < d_timers.size() + 1; ++j) {
      unsigned int n = ordered_list[j];

      /*
       * Check the print threshold to see if we should print this timer.
       */
      double frac = computePercentageDouble(
            timer_values[n][last_nonzero_column],
            timer_values[d_timers.size()][last_nonzero_column]);

      if (frac > d_print_threshold) {

         os << std::setw(maxlen + 3) << timer_names[n].c_str() << ' ';

         /*
          * Print columns.
          */
         for (i = 0; i < 4; ++i) {

            /*
             * Print column values only title is non-null (i.e. not "")
             */
            if (!column_titles[i].empty()) {

               /*
                * Print percentages for columns 0-2
                */
               if (i < 3) {
                  int k = column_ids[i];

                  if (d_print_percentage) {
#ifndef LACKS_SSTREAM
                     int perc = computePercentageInt(timer_values[n][k],
                           timer_values[d_timers.size()][k]);
                     std::ostringstream out;
                     out << timer_values[n][k] << " (" << perc << "%)";
                     os << std::setw(15) << out.str().c_str() << "  ";
#else
                     os << std::setw(15) << timer_values[n][k] << "  ";
#endif
                  } else {
                     os << std::setw(15) << timer_values[n][k] << "  ";
                  }

               } else {

                  /*
                   * Print column 3 - the processor holding processor ID
                   * with max times (don't do for TOTAL TIME - this is
                   * meaningless since all processors are synchronized
                   * before and after this call).
                   */
                  if (n < d_timers.size()) {
                     os << std::setw(15) << max_processor_id[n][max_array_id];
                  }

               } // column 3

            }  // if column title is non-null

         } // loop over columns

         os << std::endl;

      } //  matches d_print_threshold conditions

   } // loop over timers

   delete[] ordered_list;

   os << ascii_line << std::endl;
   os.setf(std::ios::right);
#else
   NULL_USE(table_title);
   NULL_USE(column_titles);
   NULL_USE(timer_names);
   NULL_USE(max_processor_id);
   NULL_USE(max_array_id);
   NULL_USE(column_ids);
   NULL_USE(timer_values);
   NULL_USE(os);

#endif
}

void
TimerManager::printOverhead(
   const std::vector<std::string>& timer_names,
   const double timer_values[][18],
   std::ostream& os)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
   std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
   std::string ascii_line = ascii_line1;
   ascii_line += ascii_line2;

   /*
    * Left-justify all output in this method.
    */
   os.setf(std::ios::left);

   os << ascii_line
      << "TIMER OVERHEAD STATISTICS \n"
      << ascii_line;

   /*
    * Determine maximum name length for formatting
    */
   int maxlen = 10;
   for (unsigned int n = 0; n < d_timers.size(); ++n) {
      int i = static_cast<int>(timer_names[n].size());
      if (i > maxlen) maxlen = i;
   }

   /*
    * Print table header
    */
   os << std::setw(maxlen + 3) << "Timer Name"
      << std::setw(25) << "Number Accesses" << "  "
      << std::setw(25) << "Estimated Cost"
      << std::endl;

   /*
    * Compute totals: total number of REGISTERED accesses and total cost.
    * Total cost includes inactive timer costs.
    */
   int total_inactive_accesses = 0;
   for (size_t i = 0; i < d_inactive_timers.size(); ++i) {
      total_inactive_accesses += d_inactive_timers[i]->getNumberAccesses();
   }

   double est_cost = d_timer_inactive_access_time * total_inactive_accesses;
   double total_est_cost = est_cost;

   int total_accesses = 0;
   for (size_t n = 0; n < d_timers.size(); ++n) {
      total_accesses += d_timers[n]->getNumberAccesses();
   }
   est_cost = d_timer_active_access_time * total_accesses;

   /*
    * If we are keeping exclusive or concurrent times, each access costs
    * roughly four times as much.  Make this correction here...
    */
   if (d_print_exclusive || d_print_concurrent) {
      est_cost *= 4.;
   }
   total_est_cost += est_cost;

   /*
    * Output the rows of the table.  Start first with the inactive timers...
    */
   int num_accesses = total_inactive_accesses;
   est_cost = d_timer_inactive_access_time * num_accesses;
   int perc = computePercentageInt(est_cost, total_est_cost);

   os << std::setw(maxlen + 3) << "inactive timers"
      << std::setw(25) << num_accesses << "  ";
#ifndef LACKS_SSTREAM
   std::ostringstream out;
   out << est_cost << " (" << perc << "%)";
   os << std::setw(25) << out.str().c_str();
#else
   os << std::setw(25) << est_cost << " (" << perc << "%)";
#endif
   os << std::endl;

   /*
    * Now print the rest of the timers.  While we are cycling through them,
    * add up the total cost and print it at the end...
    */

   for (unsigned int n = 0; n < d_timers.size(); ++n) {

      num_accesses = d_timers[n]->getNumberAccesses();
      est_cost = d_timer_active_access_time * num_accesses;

      /*
       * If we are keeping exclusive or concurrent times, each access costs
       * roughly four times as much.  Make this correction here...
       */
      if (d_print_exclusive || d_print_concurrent) {
         est_cost *= 4.;
      }

      perc = computePercentageInt(est_cost, total_est_cost);

      os << std::setw(maxlen + 3) << timer_names[n].c_str()
         << std::setw(25) << num_accesses << "  ";
#ifndef LACKS_SSTREAM
      std::ostringstream out2;
      out2 << est_cost << " (" << perc << "%)";
      os << std::setw(25) << out2.str().c_str();
#else
      os << std::setw(25) << est_cost << " (" << perc << "%)";
#endif
      os << std::endl;
   }

   /*
    * Output the totals.
    */
   os << std::setw(maxlen + 3) << "TOTAL:"
      << std::setw(25) << total_accesses << "  "
      << std::setw(25) << total_est_cost
      << "\n" << std::endl;

   /*
    * Compare the total estimated cost with overall program wallclock time.
    * If it is a significant percentage (> 5%) print a warning
    */
   double perc_dbl = computePercentageDouble(total_est_cost,
         timer_values[d_timers.size()][11]);

   os << "Estimated Timer Costs as a percentage of overall Wallclock Time: "
      << perc_dbl << "% \n";
   if (perc_dbl > 5.) {
      os << "WARNING:  TIMERS ARE USING A SIGNIFICANT FRACTION OF RUN TIME"
         << std::endl;
   }

   os << ascii_line << std::endl;
   os.setf(std::ios::right);
#else
   NULL_USE(timer_names);
   NULL_USE(timer_values);
   NULL_USE(os);

#endif // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::printConcurrent(
   std::ostream& os)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
   std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
   std::string ascii_line = ascii_line1;
   ascii_line += ascii_line2;

   os << ascii_line
      << "CONCURRENT TIMERS\n"
      << ascii_line;

   /*
    * Determine maximum name length for formatting
    */
   int maxlen = 10;
   for (size_t n = 0; n < d_timers.size(); ++n) {
      int i = int((d_timers[n]->getName()).size());
      if (i > maxlen) maxlen = i;

   }

   /*
    * Print table header
    */
   os << std::setw(maxlen + 3) << "Timer Name"
      << std::setw(25) << "Nested Timers"
      << std::endl;

   /*
    * Output the rows of the table.
    */

   for (size_t n = 0; n < d_timers.size(); ++n) {

      os << std::setw(maxlen + 3) << d_timers[n]->getName().c_str();

      int count = 0;
      for (size_t i = 0; i < d_timers.size(); ++i) {
         if (d_timers[n]->isConcurrentTimer(*d_timers[i])) {
            ++count;
         }
      }
      if (count == 0) {
         os << std::setw(25) << "none " << std::endl;
      } else {
         /*
          * Format it like:    Timer Name      Concurrent Timer #1
          *                                    Concurrent Timer #2
          *                                    ...
          * Use "count" variable defined above to identify the first
          * line or subsequent lines.
          */
         count = 0;
         for (size_t j = 0; j < d_timers.size(); ++j) {
            if (d_timers[n]->isConcurrentTimer(*d_timers[j])) {
               if (count == 0) {
                  os << std::setw(25) << d_timers[j]->getName().c_str()
                     << std::endl;
               } else {
                  os << std::setw(maxlen + 3) << " "
                     << d_timers[j]->getName().c_str() << std::endl;
               }
               ++count;
            }
         }
      }

   }
   os << ascii_line << std::endl;
#else
   NULL_USE(os);
#endif // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::checkConsistencyAcrossProcessors()
{
#ifdef ENABLE_SAMRAI_TIMERS
   /*
    * Due to the difficulty of comparing strings using MPI calls,
    * we do a rough consistency check of
    * 1. the number of timers and
    * 2. the length of each timer name.
    *
    * Steps:
    * 1. Do global reductions to get the max number of timers
    *    and the max lengths of each timer name.
    * 2. Issue a warning if the number of timers is inconsistent.
    *    This inconsistency would be found on all processes
    *    except those with the biggest number of timers.
    * 3. Issue a warning for each individual timer if
    *    its name length is less than the max length of
    *    all timers at the same index in the timer manager.
    *    Even if the number of timers are consistent, this
    *    would find wrong timer orderings or inconsistent
    *    timer names, unless the errors are for timer names
    *    with identical lengths.
    * 4. Go global reductions to get the number of inconsistencies
    *    of other processes.  Turn off printing of sum and max
    *    if any processes has inconsistencies.
    *
    * In the future, we may want to convert the strings into
    * their MD5 signatures and compare those as integers.
    */

   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());

   unsigned int max_num_timers = static_cast<unsigned int>(d_timers.size());
   if (mpi.getSize() > 1) {
      int i = static_cast<int>(d_timers.size());
      mpi.Allreduce(&i, &max_num_timers, 1, MPI_INT, MPI_MAX);
   }

   std::vector<int> max_timer_lengths(max_num_timers);
   std::vector<int> rank_of_max(max_num_timers, mpi.getRank());

   for (unsigned int i = 0; i < max_num_timers; ++i) {
      max_timer_lengths[i] =
         i < static_cast<unsigned int>(d_timers.size())
         ? static_cast<int>(d_timers[i]->getName().size()) : 0;
   }

   if (mpi.getSize() > 1) {
      mpi.AllReduce(&max_timer_lengths[0],
         max_num_timers,
         MPI_MAXLOC,
         &rank_of_max[0]);
   }

   int inconsistency_count = 0;

   if (max_num_timers > d_timers.size()) {
      TBOX_WARNING("Timer selections across processors were determined to be"
         << "\ninconsistent.  This processor has only "
         << d_timers.size() << " while some has " << max_num_timers
         << ".\nThe consistency check"
         << "\nwill continue for this process, but checking only\n"
         << d_timers.size() << " timers."
         << "\nIt is not possible to print global"
         << "\nsummed or max timer information." << std::endl);
      ++inconsistency_count;
   }

   for (unsigned int i = 0; i < d_timers.size(); ++i) {
      if (max_timer_lengths[i] != int(d_timers[i]->getName().size())) {
         TBOX_WARNING("Timer[" << i << "]: " << d_timers[i]->getName()
                               << "\nis not consistent across all processors."
                               << "\nOther timer[" << i << "] has up to "
                               << max_timer_lengths[i] << " characters in their names."
                               << "\nIt is not possible to print global"
                               << "\nsummed or max timer information."
                               << std::endl);
         ++inconsistency_count;
      }
   }

   int max_inconsistency_count = inconsistency_count;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&inconsistency_count,
         &max_inconsistency_count,
         1,
         MPI_INT,
         MPI_MAX);
   }
   if (max_inconsistency_count > 0) {
      d_print_summed = false;
      d_print_max = false;
      if (inconsistency_count == 0) {
         TBOX_WARNING("Though this process found no timer inconsistencies,"
            << "\nother processes did.  It is not possible to print"
            << "\nglobal summed or max timer information." << std::endl);
      }
   }

   /*
    * NOTE:  It might be nice to someday add the capability to remove the
    * inconsistent timers and print the max/summed values of the
    * consistent ones.   Unfortunately, this is tough to implement.  If it
    * were just a matter of comparing timer names across processors it would be
    * easy. But with MPI, only ints and doubles can be exchanged across
    * processors so it is difficult to make string comparisons.
    * It is possible to compare the MD5 sum of the strings,
    * but that may make SAMRAI dependent on the MD5 library.
    */
#endif // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::buildTimerArrays(
   double timer_values[][18],
   int max_processor_id[][2],
   std::vector<std::string>& timer_names)
{
#ifdef ENABLE_SAMRAI_TIMERS
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   /*
    * timer_values - 2D array dimensioned [d_timers.size()][18]
    *     For each timer, there are 18 potential values which may be of
    *     interest.  This array stores them if they are requested.
    * max_processor_id - 2D array dimensioned [d_timers.size()][2]
    *     Holds the value of the processor that used the maximum amount
    *     of time.  [0] is for exclusive time, while [1] is for total time.
    */

   /*
    * Initialize arrays
    */
   for (unsigned int n = 0; n < d_timers.size() + 1; ++n) {
      timer_names[n] = "";
      max_processor_id[n][0] = 0;
      max_processor_id[n][1] = 0;
      for (int i = 0; i < 18; ++i) {
         timer_values[n][i] = 0.;
      }
   }

   /*
    * Build arrays.
    */
   for (unsigned int n = 0; n < d_timers.size(); ++n) {
      timer_names[n] = d_timers[n]->getName();

      /*
       *  Build timer_values[n][m] array:
       *    m = 0 :  processor exclusive user time
       *    m = 1 :  processor exclusive sys time
       *    m = 2 :  processor exclusive wall time
       *    m = 3 :  summed exclusive user time
       *    m = 4 :  summed exclusive sys time
       *    m = 5 :  summed exclusive wall time
       *    m = 6 :  max exclusive user time
       *    m = 7 :  max exclusive sys time
       *    m = 8 :  max exclusive wall time
       *    m = 9 :  processor total user time
       *    m = 10 :  processor total sys time
       *    m = 11 :  processor total wall time
       *    m = 12 :  summed total user time
       *    m = 13 :  summed total sys time
       *    m = 14 :  summed total wall time
       *    m = 15 :  max total user time
       *    m = 16 :  max total sys time
       *    m = 17 :  max total wall time
       */

      for (int k = 0; k < 2; ++k) {
         for (int j = 0; j < 3; ++j) {

            if ((k == 0 && d_print_exclusive) ||
                (k == 1 && d_print_total)) {

               if ((j == 0 && d_print_processor) ||
                   (j == 1 && d_print_summed) ||
                   (j == 2 && d_print_max)) {

                  if (k == 0 && j == 0) {
                     if (d_print_user) {
                        timer_values[n][0] =
                           d_timers[n]->getExclusiveUserTime();
                     }
                     if (d_print_sys) {
                        timer_values[n][1] =
                           d_timers[n]->getExclusiveSystemTime();
                     }
                     if (d_print_wall) {
                        timer_values[n][2] =
                           d_timers[n]->getExclusiveWallclockTime();
                     }
                  } else if (k == 0 && j == 1) {
                     if (d_print_user) {
                        timer_values[n][3] =
                           d_timers[n]->getExclusiveUserTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][3], 1, MPI_SUM);
                        }
                     }
                     if (d_print_sys) {
                        timer_values[n][3] =
                           d_timers[n]->getExclusiveSystemTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][4], 1, MPI_SUM);
                        }
                     }
                     if (d_print_wall) {
                        timer_values[n][3] =
                           d_timers[n]->getExclusiveWallclockTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][5], 1, MPI_SUM);
                        }
                     }
                  } else if (k == 0 && j == 2) {
                     if (d_print_user) {
                        double user_time =
                           d_timers[n]->getExclusiveUserTime();
                        if (mpi.getSize() > 1) {
                           mpi.Allreduce(&user_time, &timer_values[n][6], 1, MPI_DOUBLE, MPI_MAX);
                        }
                     }
                     if (d_print_sys) {
                        double sys_time =
                           d_timers[n]->getExclusiveSystemTime();
                        if (mpi.getSize() > 1) {
                           mpi.Allreduce(&sys_time, &timer_values[n][7], 1, MPI_DOUBLE, MPI_MAX);
                        }
                     }
                     if (d_print_wall) {
                        timer_values[n][8] = d_timers[n]->getExclusiveWallclockTime();
                        max_processor_id[n][0] = mpi.getRank();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][8],
                              1,
                              MPI_MAXLOC,
                              &max_processor_id[n][0]);
                        }
                     }

                  } else if (k == 1 && j == 0) {
                     if (d_print_user) {
                        timer_values[n][9] =
                           d_timers[n]->getTotalUserTime();
                     }
                     if (d_print_sys) {
                        timer_values[n][10] =
                           d_timers[n]->getTotalSystemTime();
                     }
                     if (d_print_wall) {
                        timer_values[n][11] =
                           d_timers[n]->getTotalWallclockTime();
                     }
                  } else if (k == 1 && j == 1) {
                     if (d_print_user) {
                        timer_values[n][12] =
                           d_timers[n]->getTotalUserTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][12], 1, MPI_SUM);
                        }
                     }
                     if (d_print_sys) {
                        timer_values[n][13] =
                           d_timers[n]->getTotalSystemTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][13], 1, MPI_SUM);
                        }
                     }
                     if (d_print_wall) {
                        timer_values[n][14] =
                           d_timers[n]->getTotalWallclockTime();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][14], 1, MPI_SUM);
                        }
                     }
                  } else if (k == 1 && j == 2) {
                     if (d_print_user) {
                        double user_time =
                           d_timers[n]->getTotalUserTime();
                        if (mpi.getSize() > 1) {
                           mpi.Allreduce(&user_time, &timer_values[n][15], 1, MPI_DOUBLE, MPI_MAX);
                        }
                     }
                     if (d_print_sys) {
                        double sys_time =
                           d_timers[n]->getTotalSystemTime();
                        if (mpi.getSize() > 1) {
                           mpi.Allreduce(&sys_time, &timer_values[n][16], 1, MPI_DOUBLE, MPI_MAX);
                        }
                     }
                     if (d_print_wall) {
                        timer_values[n][17] = d_timers[n]->getTotalWallclockTime();
                        max_processor_id[n][1] = mpi.getRank();
                        if (mpi.getSize() > 1) {
                           mpi.AllReduce(&timer_values[n][17], 1, MPI_MAXLOC,
                              &max_processor_id[n][1]);
                        }
                     }

                  }

               }   // if j
            }   // if k

         } // loop over j
      } // loop over k

   } // loop over n

   /*
    * Store main_timer data in timer_values[d_timers.size()][] location.  Max
    * time and exclusive time are not determined since these don't really
    * mean anything for an overall measurement of run time.
    */
   timer_names[static_cast<int>(d_timers.size())] = "TOTAL RUN TIME:";
   if (d_print_user) {
      double main_time = d_main_timer->getTotalUserTime();
      timer_values[d_timers.size()][0] = main_time;
      timer_values[d_timers.size()][3] = main_time;
      timer_values[d_timers.size()][6] = main_time;
      timer_values[d_timers.size()][9] = main_time;
      timer_values[d_timers.size()][12] = main_time;
      timer_values[d_timers.size()][15] = main_time;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][3], 1, MPI_DOUBLE, MPI_SUM);
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][12], 1, MPI_DOUBLE, MPI_SUM);
      }
   }
   if (d_print_sys) {
      double main_time = d_main_timer->getTotalSystemTime();
      timer_values[d_timers.size()][1] = main_time;
      timer_values[d_timers.size()][4] = main_time;
      timer_values[d_timers.size()][7] = main_time;
      timer_values[d_timers.size()][10] = main_time;
      timer_values[d_timers.size()][13] = main_time;
      timer_values[d_timers.size()][16] = main_time;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][4], 1, MPI_DOUBLE, MPI_SUM);
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][13], 1, MPI_DOUBLE, MPI_SUM);
      }
   }
   if (d_print_wall) {
      double main_time = d_main_timer->getTotalWallclockTime();
      timer_values[d_timers.size()][2] = main_time;
      timer_values[d_timers.size()][5] = main_time;
      timer_values[d_timers.size()][8] = main_time;
      timer_values[d_timers.size()][11] = main_time;
      timer_values[d_timers.size()][14] = main_time;
      timer_values[d_timers.size()][17] = main_time;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][5], 1, MPI_DOUBLE, MPI_SUM);
         mpi.Allreduce(&main_time, &timer_values[d_timers.size()][14], 1, MPI_DOUBLE, MPI_SUM);
      }
   }
#else
   NULL_USE(timer_values);
   NULL_USE(max_processor_id);
   NULL_USE(timer_names);
#endif // ENABLE_SAMRAI_TIMERS
}

/*
 *************************************************************************
 *
 * Build ordered_list which specifies order of timers - max to min.
 *
 *************************************************************************
 */
void
TimerManager::buildOrderedList(
   const double timer_values[][18],
   const int column,
   int index[],
   const int array_size)
{
#ifdef ENABLE_SAMRAI_TIMERS
   /*
    * initialize the arrays
    */
   std::vector<double> timer_vals(array_size);
   for (int i = 0; i < array_size; ++i) {
      index[i] = i;
      timer_vals[i] = timer_values[i][column];
   }

   /*
    * Do a quicksort on timer_values array to build index array
    * ordered_list.
    */
   quicksort(timer_vals, index, 0, array_size - 1);
#else
   NULL_USE(timer_values);
   NULL_USE(column);
   NULL_USE(index);
   NULL_USE(array_size);
#endif // ENABLE_SAMRAI_TIMERS
}

/*
 *************************************************************************
 *
 * Sort array a largest to smallest using quicksort algorithm.
 *
 *************************************************************************
 */
void
TimerManager::quicksort(
   const std::vector<double>& a,
   int index[],
   int lo,
   int hi)
{
#ifdef ENABLE_SAMRAI_TIMERS
   if (hi <= lo) return;

   /*
    * Put a[i] into position for i between lo and hi
    * (i.e. pivot point)
    */
   int i = lo - 1;
   int j = hi;
   double v = a[index[hi]];
   for ( ; ; ) {
      while (a[index[++i]] > v)
         NULL_STATEMENT;
      while (v > a[index[--j]]) {
         if (j == lo) break;
      }
      if (i >= j) break;

      // exchange i, j indices
      int temp = index[i];
      index[i] = index[j];
      index[j] = temp;
   }
   // exchange i, hi indices
   int temp = index[i];
   index[i] = index[hi];
   index[hi] = temp;

   quicksort(a, index, lo, i - 1);
   quicksort(a, index, i + 1, hi);
#else
   NULL_USE(a);
   NULL_USE(index);
   NULL_USE(lo);
   NULL_USE(hi);
#endif // ENABLE_SAMRAI_TIMERS
}

/*
 *************************************************************************
 *
 * Operation performed many times throughout print routines. We have
 * to have some safety checks to avoid divide-by-zero errors in some
 * cases.  Thus, I just made it a function.
 *
 *************************************************************************
 */
int
TimerManager::computePercentageInt(
   const double frac,
   const double tot)
{
   /*
    *  Put a cap on the percentage at 1000.  If tot = 0, this if
    *  test should catch it.
    */
   int perc = 0;
#ifdef ENABLE_SAMRAI_TIMERS
   if (tot > 0.1 * frac) {
      perc = int(frac / tot * 100.);
   } else {
      perc = 1000;
   }
#else
   NULL_USE(frac);
   NULL_USE(tot);
#endif
   return perc;
}

double
TimerManager::computePercentageDouble(
   const double frac,
   const double tot)
{
   /*
    *  Put a cap on the percentage at 1000.  If tot = 0, this if
    *  test should catch it.
    */
   double perc = 0;
#ifdef ENABLE_SAMRAI_TIMERS
   if (tot > 0.1 * frac) {
      perc = frac / tot * 100.;
   } else {
      perc = 1000;
   }
#else
   NULL_USE(frac);
   NULL_USE(tot);
#endif
   return perc;
}

/*
 *************************************************************************
 *
 * Private member function for processing input data.
 *
 *************************************************************************
 */

void
TimerManager::getFromInput(
   const std::shared_ptr<Database>& input_db)
{
#ifdef ENABLE_SAMRAI_TIMERS
   if (input_db) {

      d_print_exclusive =
         input_db->getBoolWithDefault("print_exclusive", false);

      d_print_total = input_db->getBoolWithDefault("print_total", true);

      d_print_processor =
         input_db->getBoolWithDefault("print_processor", true);

      d_print_max = input_db->getBoolWithDefault("print_max", false);

      d_print_summed = input_db->getBoolWithDefault("print_summed", false);

      d_print_user = input_db->getBoolWithDefault("print_user", false);

      d_print_sys = input_db->getBoolWithDefault("print_sys", false);

      d_print_wall = input_db->getBoolWithDefault("print_wall", true);

      d_print_percentage =
         input_db->getBoolWithDefault("print_percentage", true);

      d_print_concurrent =
         input_db->getBoolWithDefault("print_concurrent", false);

      d_print_timer_overhead =
         input_db->getBoolWithDefault("print_timer_overhead", false);

      d_print_threshold =
         input_db->getDoubleWithDefault("print_threshold", 0.25);

      std::vector<std::string> timer_list;
      if (input_db->keyExists("timer_list")) {
         timer_list = input_db->getStringVector("timer_list");
      }

      /*
       *  Step thru the input list and call addTimerToNameLists to add
       *  the input file entry to the d_package_names,
       *  d_class_names, and d_class_method_names lists.
       */
      for (int i = 0; i < static_cast<int>(timer_list.size()); ++i) {
         std::string entry = timer_list[i];
         addTimerToNameLists(entry);
      }
      d_length_package_names = static_cast<int>(d_package_names.size());
      d_length_class_names = static_cast<int>(d_class_names.size());
      d_length_class_method_names = static_cast<int>(d_class_method_names.size());

   }
#else
   NULL_USE(input_db);
#endif // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::addTimerToNameLists(
   const std::string& name)
{
#ifdef ENABLE_SAMRAI_TIMERS
   /*
    * Evaluate whether the name is a package, class, or class::method
    * combination.  This parser supports inputs of the form:
    *
    *    *::*::*         - ALL timers added
    *    Package::*::*   - "Package" added to package list.
    *    Class           - "Class" added to class list.
    *    *::Class        - "Class" added to class list.
    *    Class::*        - "Class" added to class list.
    *    *::Class::*     - "Class" added to class list.
    *    Package::Class::method  - "Class::method" put to class_method list
    *    Class::method   - "Class::method" put to class_method list
    */

   std::string::size_type position, string_length;

   /*
    *  Step thru the input list and form the d_package_names,
    *  d_class_names, and d_class_method_names lists.
    */

   if (!name.empty()) {    // Nested if #1

      std::string entry = name;

      /*
       *  Once we have determined whether the entry is a package,
       *  class, or class::method, use this bool to jump to the next
       *  loop entry.
       */
      bool determined_entry = false;

      /*
       * Check if its a wildcard entry - "*::*::*".  If so, add all package
       * names to the package name list.
       */
      position = entry.find("*::*::*");  // if not found, "position" runs off
                                         // end of entry so pos > entry.size()
      if (position < entry.size()) {
         d_package_names.push_front("algs");
         d_package_names.push_front("apps");
         d_package_names.push_front("appu");
         d_package_names.push_front("geom");
         d_package_names.push_front("hier");
         d_package_names.push_front("math");
         d_package_names.push_front("mesh");
         d_package_names.push_front("pdat");
         d_package_names.push_front("solv");
         d_package_names.push_front("tbox");
         d_package_names.push_front("xfer");
         determined_entry = true;
      }

      /*
       * Is it a package?  Look for "::*::*" string.  If its there,
       * parse it off and add the package to the package list.
       */
      if (!determined_entry) {
         position = entry.find("::*::*");
         if (position < entry.size()) {
            entry = entry.substr(0, position);
            d_package_names.push_front(entry);
            determined_entry = true;
         }
      }

      if (!determined_entry) {    // Nested if #2

         /*
          * Is it a class?  If it doesn't have any "::", it must be a class.
          */
         position = entry.find("::");
         if (position > entry.size()) {
            d_class_names.push_front(entry);
            determined_entry = true;
         }
         if (!determined_entry) {    // Nested if #3

            /*
             * At this point, we know the entry has a "::" but wasn't
             * identified as a package.  There are several options that
             * might make Foo a class entry:
             *  1) Foo::
             *  2) *::Foo::
             *  3) Package::Foo::
             *  4) *::Foo
             * Parse these as follows:  First, look for existence of "::*"
             * at the end of the entry.  This will identify the first 3
             * options.  Next look for existence of "*::" at front of the
             * string.  This will identify the fourth choice.
             *
             * Check for first three options...
             */
            string_length = entry.size();
            std::string substring = entry.substr(string_length - 3,
                  string_length);
            if (substring == "::*") {
               entry = entry.substr(0, string_length - 3);

               /*
                * If a preceeding "::" exists at the front of the entry
                * (i.e. option 2 and 3), parse off anything before it.
                */
               position = entry.find("::");
               if (position < entry.size()) {
                  entry = entry.substr(position + 2);
               }
               d_class_names.push_front(entry);
               determined_entry = true;
            }

            if (!determined_entry) {    // Nested if #4

               /*
                * Check for option 4.  The entry has a preceeding *::. Do not
                * accept case where there is a second "::" followed by anything
                * but "*", since this is a class::method combination.
                *
                */
               substring = entry.substr(0, 3);
               if (substring == "*::") {
                  entry = entry.substr(3);
                  position = entry.find("::");

                  /*
                   * There is no second "::".  Accept the entry as a class.
                   */
                  if (position > entry.size()) {
                     d_class_names.push_front(entry);
                     determined_entry = true;
                  } else {

                     /*
                      * There is a second "::".  See if it is followed by a
                      * "*".  If so, parse off the "::*" and accept entry as
                      * a class.  If not, let it be determined below to be a
                      * class::method entry.
                      */
                     string_length = entry.size();
                     substring = entry.substr(string_length - 1, string_length);
                     if (substring == "*") {
                        entry = entry.substr(0, string_length - 3);
                        d_class_names.push_front(entry);
                        determined_entry = true;
                     }
                  }
               }

               if (!determined_entry) {    // Nested if #5

                  /*
                   * The entry has not been identified as either a package or
                   * a class.  It must be a class::method combination. There
                   * are three options for entering class::method combinations:
                   *  1) Package::Foo::method
                   *  2) *::Foo::method
                   *  3) Foo::method
                   * We only want to maintain "Foo::method" in the package
                   * list.  Check first if there are two "::" in the entry.
                   * If there are, parse of whatever is in front of the
                   * first "::".  If not, just use the entry as is.
                   */
                  position = entry.find("::");
                  if (position < entry.size()) {
                     substring = entry.substr(position + 2);
                     position = substring.find("::");
                     if (position < substring.size()) {

                        /*
                         * There *are* two "::" so entry must contain a
                         * package.  Parse it off.
                         */
                        position = entry.find("::");
                        entry = entry.substr(position + 2);
                     }
                  }
                  d_class_method_names.push_front(entry);

               } // Nested if #5

            } // Nested if #4

         } // Nested if #3

      } // Nested if #2

   } // Nested if #1
#else
   NULL_USE(name);
#endif // ENABLE_SAMRAI_TIMERS
}

double
TimerManager::computeOverheadConstantActiveOrInactive(
   bool active)
{
#ifdef ENABLE_SAMRAI_TIMERS
   std::string outer_name("TimerManger::Outer");
   std::shared_ptr<Timer> outer_timer(
      TimerManager::getManager()->getTimer(outer_name, true));

   std::string inner_name("TimerMangerInner");
   std::shared_ptr<Timer> inner_timer(
      TimerManager::getManager()->getTimer(inner_name, active));

   const int ntest = 1000;
   for (int i = 0; i < ntest; ++i) {
      outer_timer->start();
      inner_timer->start();
      inner_timer->stop();
      outer_timer->stop();
   }

   return (outer_timer->getTotalWallclockTime()
           - inner_timer->getTotalWallclockTime()) / (static_cast<double>(ntest));

#else
   NULL_USE(active);

   return 0.0;

#endif // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::computeOverheadConstants()
{
#ifdef ENABLE_SAMRAI_TIMERS

   if (d_timer_active_access_time < 0.0) {

      clearArrays();
      d_timer_active_access_time = computeOverheadConstantActiveOrInactive(
            false);

      clearArrays();
      d_timer_inactive_access_time = computeOverheadConstantActiveOrInactive(
            true);

      clearArrays();
   }
#endif  // ENABLE_SAMRAI_TIMERS
}

void
TimerManager::clearArrays()
{
#ifdef ENABLE_SAMRAI_TIMERS
   /*
    * Create a timer that measures overall solution time.
    */
   d_main_timer.reset(new Timer("TOTAL RUN TIME"));

   d_timers.clear();
   d_inactive_timers.clear();

   d_exclusive_timer_stack.clear();
#endif // ENABLE_SAMRAI_TIMERS
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
