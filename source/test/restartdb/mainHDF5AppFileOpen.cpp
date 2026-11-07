/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Tests HDF database in SAMRAI
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"

#include <string>
#include <memory>

using namespace SAMRAI;

#include "database_tests.h"

class RestartTester:public tbox::Serializable
{
public:
   RestartTester()
   {
      tbox::RestartManager::getManager()->registerRestartItem("RestartTester",
         this);
   }

   virtual ~RestartTester() {
   }

   void putToRestart(
      const std::shared_ptr<tbox::Database>& db) const
   {
      writeTestData(db);
   }

   void getFromRestart()
   {
      std::shared_ptr<tbox::Database> root_db(
         tbox::RestartManager::getManager()->getRootDatabase());

      std::shared_ptr<tbox::Database> db;
      if (root_db->isDatabase("RestartTester")) {
         db = root_db->getDatabase("RestartTester");
      }

      readTestData(db);
   }

};

int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("HDF5test.log");

#ifdef HAVE_HDF5

      tbox::plog << "\n--- HDF5 database tests BEGIN ---" << std::endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      RestartTester hdf_tester;

      tbox::plog << "\n--- HDF5 write database tests BEGIN ---" << std::endl;

      setupTestData();

      std::shared_ptr<tbox::HDFDatabase> database(
         new tbox::HDFDatabase("SAMRAI Restart"));
      std::string name = "./restart." + tbox::Utilities::processorToString(
            mpi.getRank()) + ".hdf5";
      hid_t file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC,
            H5P_DEFAULT, H5P_DEFAULT);
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t hdf_group = H5Gcreate(file_id,
            "SAMRAIGroup",
            0,
            H5P_DEFAULT,
            H5P_DEFAULT);
#else
      hid_t hdf_group = H5Gcreate(file_id, "SAMRAIGroup", 0);
#endif
      database->attachToFile(hdf_group);

      restart_manager->setRootDatabase(database);

      restart_manager->writeRestartToDatabase();

      tbox::plog << "\n--- HDF5 write database tests END ---" << std::endl;

      tbox::plog << "\n--- HDF5 read database tests BEGIN ---" << std::endl;

      database->close();

      restart_manager->setRootDatabase(std::shared_ptr<tbox::Database>());

      H5Fclose(file_id);

      database.reset(new tbox::HDFDatabase("SAMRAI Restart"));
      file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hdf_group = H5Gopen(file_id, "SAMRAIGroup", H5P_DEFAULT);
#else
      hdf_group = H5Gopen(file_id, "SAMRAIGroup");
#endif
      database->attachToFile(hdf_group);

      restart_manager->setRootDatabase(database);

      hdf_tester.getFromRestart();

      database->close();

      restart_manager->setRootDatabase(std::shared_ptr<tbox::Database>());

      H5Fclose(file_id);

      tbox::plog << "\n--- HDF5 read database tests END ---" << std::endl;

      tbox::plog << "\n--- HDF5 database tests END ---" << std::endl;

#endif

      if (number_of_failures == 0) {
         tbox::pout << "\nPASSED:  HDF5" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return number_of_failures;

}
