/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for hierachy coarsen/refine tests.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include <string>
#include "SAMRAI/hier/VariableDatabase.h"

#include "HierarchyTester.h"

#include <memory>

using namespace SAMRAI;

/*
 ************************************************************************
 *
 *
 ************************************************************************
 */

int main(
   int argc,
   char* argv[])
{

   int fail_count;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      std::string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << std::endl);
      } else {
         input_filename = std::string(argv[1]);
      }

      tbox::plog << "\n Starting hierarchy refine/coarsen test..." << std::endl;
      tbox::plog << "Specified input file is: " << input_filename << std::endl;

      std::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "GlobalInputs" section of the input database and set
       * values accordingly.
       */

      if (input_db->keyExists("GlobalInputs")) {
         std::shared_ptr<tbox::Database> global_db(
            input_db->getDatabase("GlobalInputs"));
         if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
            bool flag = global_db->
               getBool("call_abort_in_serial_instead_of_exit");
            tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
         }
      }

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string log_file_name = "hierarchy_test.log";
      if (main_db->keyExists("log_file_name")) {
         log_file_name = main_db->getString("log_file_name");
      }
      bool log_all_nodes = false;
      if (main_db->keyExists("log_all_nodes")) {
         log_all_nodes = main_db->getBool("log_all_nodes");
      }
      if (log_all_nodes) {
         PIO::logAllNodes(log_file_name);
      } else {
         PIO::logOnlyNodeZero(log_file_name);
      }

      std::shared_ptr<HierarchyTester> hierarchy_tester(
         new HierarchyTester(
            "HierarchyTester",
            dim,
            input_db->getDatabase("HierarchyTest")));

      tbox::plog << "\nCreating initial patch hierarchy..." << std::endl;
      hierarchy_tester->setupInitialHierarchy(input_db);
      tbox::plog << "\nInitial patch hierarchy created..." << std::endl;

      tbox::plog << "\nInput file data is ...." << std::endl;
      input_db->printClassData(plog);

      tbox::plog << "\nVariable database..." << std::endl;
      hier::VariableDatabase::getDatabase()->printClassData(tbox::plog, false);

      tbox::plog << "\nPerforming refine/coarsen patch hierarchy test..."
                 << std::endl;
      fail_count = hierarchy_tester->runHierarchyTestAndVerify();
      tbox::plog << "\n Ending hierarchy refine/coarsen test..." << std::endl;

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  hierarchy tester" << std::endl;
      }
   }

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();

   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
