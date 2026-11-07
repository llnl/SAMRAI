/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An restart manager singleton class
 *
 ************************************************************************/

#ifndef included_tbox_RestartManager
#define included_tbox_RestartManager

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/DatabaseFactory.h"

#include <string>
#include <list>
#include <memory>

namespace SAMRAI {
namespace tbox {

/**
 * Class RestartManager coordinates SAMRAI restart files (currently
 * implemented using the HDF database class) and the objects comprising
 * SAMRAI-based application code.  The manager class orchestrates opening
 * and closing the database, stores data to be written out for restart,
 * and writes out the restart data to the database.  Note that the restart
 * manager is a Singleton class that acts as a single point of control
 * for restart capabilities.  As such its constructor and destructor are
 * protected members; they are not to be called outside of this class.
 *
 * The general procedure for starting a simulation from a restart file
 * is as follows.
 *
 *
 * \li Open the restart file using openRestartFile("filename").
 * \li Get root of restart database using getRootDatabase().
 * \li Initialize simulation objects using the restart constructor
 *       for the objects.
 * \li Close the restart file using closeRestartFile().
 *
 * Technically, there is no need to close the restart file because this will
 * automatically be taken care of by the destructor for the database object.
 *
 * It is important to note in the initialization process, some objects
 * will need to be constructed in the "empty" state and filled in later
 * using some sort of getFromRestart() method.
 *
 * The process for writing out state to a restart file is somewhat more
 * complicated.  The following things need to be taken care of.
 *
 *
 * \li Each object that has state that needs to be saved for restart
 *       must be derived from the Serializable class (which
 *       responds to the putToRestart() method).
 * \li Any object that needs to save its state to the restart file
 *       must be registered with the restart manager using the
 *       registerRestartItem() method.   NOTE THAT NO TWO RESTARTABLE
 *       OBJECTS ARE ALLOWED TO HAVE THE SAME NAME STRING IDENTIFIER.
 * \li The patchdata to be written to restart need to be specified
 *       using the VariableDatabase::setPatchdataForRestart() method.
 *       This is usually taken care of by the numerical algorithm object.
 *
 * When all these items are accounted for, writing to the restart file
 * is accomplished using a writeRestartFile() method.  There are
 * two writeRestartFile() methods available.  One takes only
 * a restart directory name as an argument whereas the other takes
 * both a restart directory name and a restore number for its arguments.
 * See comments for member functions for more details.
 *
 * @see Database
 */

class RestartManager
{
public:
   /**
    * Return a pointer to the single instance of the restart manager.
    * All access to the restart manager object is through getManager().
    *
    * Note that when the manager is accessed for the first time, the
    * Singleton instance is registered with the StartupShutdownManager
    * class which destroys such objects at program completion.  Thus,
    * users of this class do not explicitly allocate or deallocate the
    * Singleton instance.
    */
   static RestartManager *
   getManager();

   /**
    * Returns true if the run is from a restart file (i.e. a restart file
    * has been opened from main()).  Returns false otherwise.
    */
   bool
   isFromRestart()
   {
      return d_is_from_restart;
   }

   /**
    * Attempts to mount, for reading, the restart file for the processor.
    * If there is no error opening the file, then the restart manager
    * mounts the restart file.
    * Returns true if open is successful; false otherwise.
    *
    * @pre hasDatabaseFactory()
    */
   bool
   openRestartFile(
      const std::string& root_dirname,
      const int restore_num,
      const int num_nodes);

   /**
    * Closes the restart file.
    */
   void
   closeRestartFile();

   /**
    * Returns a std::shared_ptr to the root of the database.
    */
   std::shared_ptr<Database>
   getRootDatabase()
   {
      return d_database_root;
   }

   /*!
    * @brief Returns true if the root of the database has been set.
    */
   bool
   hasRootDatabase()
   {
      return d_database_root.get();
   }

   /**
    * Sets the database for restore or dumps.
    *
    */
   void
   setRootDatabase(
      const std::shared_ptr<Database>& database)
   {
      if (!database) {
         d_database_root.reset();
         d_is_from_restart = false;
      } else {
         d_database_root = database;
         d_is_from_restart = true;
      }
   }

   /**
    * Sets the database for restore or dumps.
    *
    */
   void
   setDatabaseFactory(
      const std::shared_ptr<DatabaseFactory>& database_factory)
   {
      d_database_factory = database_factory;
   }

   /*!
    * @brief Returns true if the database for restore or dumps has been set.
    */
   bool
   hasDatabaseFactory()
   {
      return d_database_factory.get();
   }

   /**
    * Registers an object for restart with the given name.
    *
    * @pre !name.empty()
    * @pre obj != 0
    */
   void
   registerRestartItem(
      const std::string& name,
      Serializable* obj);

   /**
    * Removes the object with the specified name from the list of
    * restartable items.
    *
    * @pre !name.empty()
    */
   void
   unregisterRestartItem(
      const std::string& name);

   /**
    * Clear all restart items managed by the restart manager.
    */
   void
   clearRestartItems()
   {
      d_restart_items_list.clear();
   }

   /**
    * Write all objects registered to as restart objects to the
    * restart database.  The string argument is the name of the
    * root of restart directory.
    *
    * Note:  This method creates/uses a restart directory structure
    *    with 00000 as the restore number.
    */
   void
   writeRestartFile(
      const std::string& root_dirname)
   {
      writeRestartFile(root_dirname, 0);
   }

   /**
    * Write all objects registered to as restart objects to the
    * restart database.  The string argument is the name of the
    * root of restart directory.  The integer argument is the
    * identification number associated with the restart files generated.
    *
    * @pre hasDatabaseFactory()
    */
   void
   writeRestartFile(
      const std::string& root_dirname,
      const int restore_num);

   /**
    * Write all objects registered to as restart objects to the
    * restart database.
    *
    * @pre hasRootDatabase()
    */
   void
   writeRestartToDatabase();

protected:
   /**
    * The constructor for RestartManager is protected.
    * Consistent with the definition of a Singleton class, only the
    * manager object has access to the constructor for the class.
    *
    * The constructor for RestartManager initializes the root
    * data base to a NullDatabase and sets the restart flag to false.
    */
   RestartManager();

   /**
    * The destructor for the restart manager is protected, since only the
    * singleton class and subclasses may destroy the manager objects.
    */
   ~RestartManager();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    *
    * @pre !s_manager_instance
    */
   void
   registerSingletonSubclassInstance(
      RestartManager * subclass_instance);

private:
   // Unimplemented copy constructor.
   RestartManager(
      const RestartManager& other);

   // Unimplemented assignment operator.
   RestartManager&
   operator = (
      const RestartManager& rhs);

   /**
    * Write all objects registered to as restart objects to the
    * restart database.
    *
    * @pre database
    */
   void
   writeRestartFile(
      const std::shared_ptr<Database>& database);

   /*
    * Create the directory structure for the data files.
    * The directory structure created is
    *
    *   restart_dirname/
    *     restore.[restore number]/
    *       nodes.[number of processors]/
    *         proc.[processor number]
    */
   std::string
   createDirs(
      const std::string& root_dirname,
      int restore_num);

   struct RestartItem {
      std::string name;
      Serializable* obj;
   };

   /**
    * Deallocate the restart manager instance.  It is not necessary to call
    * this routine at program termination, since it is automatically called
    * by the StartupShutdownManager class.
    */
   static void
   shutdownCallback();

   static RestartManager* s_manager_instance;

   /*
    * list of objects registered to be written to the restart database
    */
   std::list<RestartManager::RestartItem> d_restart_items_list;

   std::shared_ptr<Database> d_database_root;

   /*
    * Database factory use to create new databases.
    * Defaults so HDFDatabaseFactory.
    */
   std::shared_ptr<DatabaseFactory> d_database_factory;

   bool d_is_from_restart;

   static StartupShutdownManager::Handler s_shutdown_handler;
};

}
}

#endif
