/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An restart manager singleton class
 *
 ************************************************************************/

#include <string>

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/HDFDatabaseFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/NullDatabase.h"
#include "SAMRAI/tbox/Parser.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

RestartManager * RestartManager::s_manager_instance = 0;

StartupShutdownManager::Handler
RestartManager::s_shutdown_handler(
   0,
   0,
   RestartManager::shutdownCallback,
   0,
   StartupShutdownManager::priorityRestartManager);

/*
 *************************************************************************
 *
 * Basic singleton classes to create, set, and destroy the manager
 * instance.
 *
 *************************************************************************
 */

RestartManager *
RestartManager::getManager()
{
   if (!s_manager_instance) {
      s_manager_instance = new RestartManager;
   }
   return s_manager_instance;
}

void
RestartManager::shutdownCallback()
{
   if (s_manager_instance) {
      s_manager_instance->clearRestartItems();
      delete s_manager_instance;
      s_manager_instance = 0;
   }
}

/*
 *************************************************************************
 *
 * The constructor and destructor are protected and call only be called
 * by the singleton class or its subclasses.
 *
 *************************************************************************
 */

RestartManager::RestartManager():
   d_database_root(std::make_shared<NullDatabase>()),
#ifdef HAVE_HDF5
   d_database_factory(std::make_shared<HDFDatabaseFactory>()),
#endif
   d_is_from_restart(false)
{
   clearRestartItems();
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */
RestartManager::~RestartManager()
{
}

/*
 *************************************************************************
 *
 * Mount restart_file to the empty database created in the
 * constructor and sets d_is_from_restart to true.
 * Return d_database_root.
 *
 *************************************************************************
 */

bool
RestartManager::openRestartFile(
   const std::string& root_dirname,
   const int restore_num,
   const int num_nodes)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int proc_num = mpi.getRank();

   /* create the intermediate parts of the full path name of restart file */
   std::string restore_buf = "/restore." + Utilities::intToString(
         restore_num,
         6);
   std::string nodes_buf = "/nodes." + Utilities::nodeToString(num_nodes);
   std::string proc_buf = "/proc." + Utilities::processorToString(
         proc_num);

   /* create full path name of restart file */
   std::string restart_filename = root_dirname + restore_buf
      + nodes_buf + proc_buf;

   bool open_successful = true;
   /* try to mount restart file */

   if (hasDatabaseFactory()) {

      std::shared_ptr<Database> database(d_database_factory->allocate(
                                              restart_filename));

      if (!database->open(restart_filename)) {
         TBOX_ERROR(
            "Error attempting to open restart file " << restart_filename
                                                     << "\n   No restart file for processor: "
                                                     << proc_num
                                                     << "\n   restart directory name = "
                                                     << root_dirname
                                                     << "\n   number of processors   = "
                                                     << num_nodes
                                                     << "\n   restore number         = "
                                                     << restore_num << std::endl);
         open_successful = false;
      } else {
         /* set d_database root and d_is_from_restart */
         d_database_root = database;
         d_is_from_restart = true;
      }
   } else {
      TBOX_ERROR("No DatabaseFactory supplied to RestartManager for opening "
         << restart_filename << std::endl);
   }

   return open_successful;
}

/*
 *************************************************************************
 *
 * Closes the restart file by unmounting d_database_root and setting it
 * to be a NullDatabase.
 *
 *************************************************************************
 */

void
RestartManager::closeRestartFile()
{
   if (d_database_root) {
      d_database_root->close();
      d_database_root.reset();
   }

   d_database_root.reset(new NullDatabase());
}

/*
 *************************************************************************
 *
 * Registers the object for restart by adding it to
 * d_restart_items_list.
 *
 *************************************************************************
 */
void
RestartManager::registerRestartItem(
   const std::string& name,
   Serializable* obj)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(obj != 0);

   /*
    * Run through list to see if there is another object registered
    * with the specified name.
    */
   std::list<RestartManager::RestartItem>::iterator iter =
      d_restart_items_list.begin();

   bool found_item = false;
   for ( ; !found_item && iter != d_restart_items_list.end(); ++iter) {
      found_item = (iter->name == name);
   }

   /*
    * If there are no other items registered with the specified name,
    * add the object to the restart list.  Otherwise, throw an
    * error.
    */
   if (!found_item) {
      RestartItem r_obj;
      r_obj.name = name;
      r_obj.obj = obj;

      d_restart_items_list.push_back(r_obj);

   } else {
      TBOX_ERROR("Register restart item error..."
         << "\n   Multiple objects with name `" << name << "' registered "
         << "with restart manager." << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Removes the object with the specified name from d_restart_items_list.
 *
 *************************************************************************
 */
void
RestartManager::unregisterRestartItem(
   const std::string& name)
{
   TBOX_ASSERT(!name.empty());

   std::list<RestartManager::RestartItem>::iterator iter =
      d_restart_items_list.begin();

   for ( ; iter != d_restart_items_list.end(); ++iter) {
      if (iter->name == name) {
         d_restart_items_list.erase(iter);
         break;
      }
   }
}

/*
 *************************************************************************
 *
 * Creates a new file with the given name and writes out the current
 * simulation state to the file by invoking the writeRestartFile()
 * method for all objects contained in d_restart_objects_list.
 *
 *************************************************************************
 */
void
RestartManager::writeRestartFile(
   const std::string& root_dirname,
   int restore_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   /* Create necessary directories and cd proper directory for writing */
   std::string restart_dirname = createDirs(root_dirname, restore_num);

   /* Create full path name of restart file */

   int proc_rank = mpi.getRank();

   std::string restart_filename_buf =
      "/proc." + Utilities::processorToString(proc_rank);

   std::string restart_filename = restart_dirname + restart_filename_buf;

   if (hasDatabaseFactory()) {

      std::shared_ptr<Database> new_restartDB(d_database_factory->allocate(
                                                   restart_filename));

      new_restartDB->create(restart_filename);

      writeRestartFile(new_restartDB);

      new_restartDB->close();

      new_restartDB.reset();
   } else {
      TBOX_ERROR(
         "No DatabaseFactory supplied to RestartManager for writeRestartFile "
         << restart_filename << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Write simulation state to supplied restart database.
 *
 *************************************************************************
 */
void
RestartManager::writeRestartFile(
   const std::shared_ptr<Database>& database)
{
   TBOX_ASSERT(database);

   std::list<RestartManager::RestartItem>::iterator i =
      d_restart_items_list.begin();
   for ( ; i != d_restart_items_list.end(); ++i) {
      std::shared_ptr<Database> obj_db(
         database->putDatabase(i->name));
      (i->obj)->putToRestart(obj_db);
   }
}

/*
 *************************************************************************
 *
 * Write simulation state to root database
 *
 *************************************************************************
 */
void
RestartManager::writeRestartToDatabase()
{
   if (hasRootDatabase()) {
      writeRestartFile(d_database_root);
   } else {
      TBOX_ERROR("writeRestartToDatabase has no database to write to"
         << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Creates the directory structure for the data files if they have not
 * already been created.
 *
 *************************************************************************
 */

std::string
RestartManager::createDirs(
   const std::string& root_dirname,
   int restore_num)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   int num_procs = mpi.getSize();

   std::string restore_buf = "/restore." + Utilities::intToString(
         restore_num,
         6);
   std::string nodes_buf = "/nodes." + Utilities::processorToString(
         num_procs);

   std::string full_dirname = root_dirname + restore_buf + nodes_buf;

   Utilities::recursiveMkdir(full_dirname);

   return full_dirname;
}

void
RestartManager::registerSingletonSubclassInstance(
   RestartManager* subclass_instance)
{
   if (!s_manager_instance) {
      s_manager_instance = subclass_instance;
   } else {
      TBOX_ERROR("RestartManager internal error...\n"
         << "Attemptng to set Singleton instance to subclass instance,"
         << "\n but Singleton instance already set." << std::endl);
   }
}

}
}
