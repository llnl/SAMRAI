/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/
#include "SparseDataTester.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"

using namespace sam_test;
using namespace tbox;

int main(
   int argc,
   char* argv[])
{
   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();

   SAMRAI_MPI mpi(SAMRAI_MPI::getSAMRAIWorld());
   const Dimension dim(2);

   //if (argc < 2) {
   //   TBOX_ERROR("Usage: " << argv[0] << " [input file]");
   //}

   const std::string log_name = std::string("sparse_tester.")
      + Utilities::intToString(dim.getValue(), 1) + "d.log";

   bool is_restart = false;
   //if (argc == 3) {
   //   is restart = true;
   //   log_name + ".restart"
   //}

   PIO::logAllNodes(log_name);

   int fail_count = 0;

   // Need scope for tester object so it will be destroyed
   // and cleaned up before SAMRAI finalize is called.
   {
      SparseDataTester tester(dim);

      bool success = true;
      if (is_restart) {
         //success = tester.testRestart();
         if (success)
            tbox::plog << "PASSED: Restart Test" << std::endl;
         else {
            tbox::perr << "FAILED: Restart Test" << std::endl;
            ++fail_count;
         }

      } else {
         success = tester.testConstruction();
         if (success)
            tbox::plog << "PASSED: Test 1 (construction)" << std::endl;
         else {
            tbox::perr << "FAILED: construction" << std::endl;
            ++fail_count;
         }

         success = tester.testAdd();
         if (success)
            tbox::plog << "PASSED: Test 2 addItems" << std::endl;
         else {
            tbox::perr << "FAILED: addItems" << std::endl;
            ++fail_count;
         }

         success = tester.testRemove();
         if (success)
            tbox::plog << "PASSED: Test 3 removeItems" << std::endl;
         else {
            tbox::perr << "FAILED: remove items" << std::endl;
            ++fail_count;
         }

         success = tester.testCopy();
         if (success)
            tbox::plog << "PASSED: Test 4: copy" << std::endl;
         else {
            tbox::perr << "FAILED: copy items" << std::endl;
            ++fail_count;
         }

         success = tester.testCopy2();
         if (success)
            tbox::plog << "PASSED: Test 5: copy2 " << std::endl;
         else {
            tbox::perr << "FAILED: copy2 " << std::endl;
            ++fail_count;
         }

         success = tester.testPackStream();
         if (success)
            tbox::plog << "PASSED: Test 6: packStream" << std::endl;
         else {
            tbox::perr << "FAILED: packStream " << std::endl;
            ++fail_count;
         }

         success = tester.testPackStream(0);
         if (success)
            tbox::plog << "PASSED: Test 6a: packStream empty SparseData" << std::endl;
         else {
            tbox::perr << "FAILED: packStream empty SparseData " << std::endl;
            ++fail_count;
         }

         success = tester.testDatabaseInterface();
         if (success)
            tbox::plog << "PASSED: Test 7: database interface" << std::endl;
         else {
            tbox::perr << "FAILED: database interface " << std::endl;
            ++fail_count;
         }

         tester.testTiming();
      }
   }


   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  sparse data test" << std::endl;
   }

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return fail_count;
}
