/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test driver for the SAMRAI input database
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

//#include <stdlib.h>

#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <vector>
#include <memory>

using namespace SAMRAI;

int main(
   int argc,
   char** argv)
{
   /*
    * Initialize tbox::MPI and SAMRAI.  Enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   tbox::PIO::logOnlyNodeZero("inputdb.log");

   int fail_count = 0;

   {
      std::string input_filename = argv[1];

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

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

      /*******************************************************************
      * Test primitive types - int, float, double, bool, dcomplex,
      *                        std::string, box
      *******************************************************************/
      std::shared_ptr<tbox::Database> prim_type_db(
         input_db->getDatabase("PrimitiveTypes"));

      int i0_correct = 1;
      float f0_correct = 1.0;
      double d0_correct = 1.0;
      bool b0_correct = true;
      dcomplex c0_correct(1.0, 1.0);
      std::string s0_correct = "a string";

      int lower[NDIM];
      int upper[NDIM];
      for (int i = 0; i < NDIM; ++i) {
         lower[i] = 0;
         upper[i] = 9;
      }
      tbox::DatabaseBox box0_correct(tbox::Dimension(NDIM), lower, upper);

      int i0 = prim_type_db->getInteger("i0");
      float f0 = prim_type_db->getFloat("f0");
      double d0 = prim_type_db->getDouble("d0");
      bool b0 = prim_type_db->getBool("b0");
      dcomplex c0 = prim_type_db->getComplex("c0");
      std::string s0 = prim_type_db->getString("s0");
      tbox::DatabaseBox box0 = prim_type_db->getDatabaseBox("box0");

      if (i0 != i0_correct) {
         ++fail_count;
         tbox::perr << "Integer test #0 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<float>::equalEps(f0, f0_correct)) {
         ++fail_count;
         tbox::perr << "Float test #0 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<double>::equalEps(d0, d0_correct)) {
         ++fail_count;
         tbox::perr << "Double test #0 FAILED" << std::endl;
      }
      if (b0 != b0_correct) {
         ++fail_count;
         tbox::perr << "Bool test #0 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<dcomplex>::equalEps(c0, c0_correct)) {
         ++fail_count;
         tbox::perr << "Complex test #0 FAILED" << std::endl;
      }
      if (s0 != s0_correct) {
         ++fail_count;
         tbox::perr << "String test #0 FAILED" << std::endl;
      }
      if (!(box0 == box0_correct)) {
         ++fail_count;
         tbox::perr << "Box test #0 FAILED" << std::endl;
      }

      /*******************************************************************
       * Test Arrays of primitive types
       ******************************************************************/
      const int nsize = 5; // size of arrays

      /*
       * "Smart" arrays
       */
      std::shared_ptr<tbox::Database> smart_array_db(
         input_db->getDatabase("SmartArrays"));

      std::vector<int> i1_correct(5);
      std::vector<float> f1_correct(5);
      std::vector<double> d1_correct(5);
      std::vector<bool> b1_correct(5);
      std::vector<dcomplex> c1_correct(5);
      std::vector<std::string> s1_correct(5);
      std::vector<tbox::DatabaseBox> box1_correct(5);

      for (int i = 0; i < nsize; ++i) {
         i1_correct[i] = i0_correct;
         f1_correct[i] = f0_correct;
         d1_correct[i] = d0_correct;
         b1_correct[i] = b0_correct;
         c1_correct[i] = c0_correct;
         s1_correct[i] = s0_correct;
         box1_correct[i] = box0_correct;
      }

      std::vector<int> i1 = smart_array_db->getIntegerVector("i1");
      std::vector<float> f1 = smart_array_db->getFloatVector("f1");
      std::vector<double> d1 = smart_array_db->getDoubleVector("d1");
      std::vector<bool> b1 = smart_array_db->getBoolVector("b1");
      std::vector<dcomplex> c1 = smart_array_db->getComplexVector("c1");
      std::vector<std::string> s1 = smart_array_db->getStringVector("s1");
      std::vector<tbox::DatabaseBox> box1 =
         smart_array_db->getDatabaseBoxVector("box1");

      for (int i = 0; i < nsize; ++i) {
         if (i1[i] != i1_correct[i]) {
            ++fail_count;
            tbox::perr << "Integer test #1 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<float>::equalEps(f1[i], f1_correct[i])) {
            ++fail_count;
            tbox::perr << "Float test #1 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<double>::equalEps(d1[i], d1_correct[i])) {
            ++fail_count;
            tbox::perr << "Double test #1 FAILED" << std::endl;
         }
         if (b1[i] != b1_correct[i]) {
            ++fail_count;
            tbox::perr << "Bool test #1 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<dcomplex>::equalEps(c1[i], c1_correct[i])) {
            ++fail_count;
            tbox::perr << "Complex test #1 FAILED" << std::endl;
         }
         if (s1[i] != s1_correct[i]) {
            ++fail_count;
            tbox::perr << "String test #1 FAILED" << std::endl;
         }
         if (!(box1[i] == box1_correct[i])) {
            ++fail_count;
            tbox::perr << "Box test #1 FAILED" << std::endl;
         }
      }

      /*
       * Basic arrays (i.e. do not use the "smart" array construct)
       */
      std::shared_ptr<tbox::Database> basic_array_db(
         input_db->getDatabase("BasicArrays"));

      int i2_correct[nsize];
      float f2_correct[nsize];
      double d2_correct[nsize];
      bool b2_correct[nsize];
      dcomplex c2_correct[nsize];
      std::string s2_correct[nsize];
      tbox::DatabaseBox box2_correct[nsize];
      for (int i = 0; i < nsize; ++i) {
         i2_correct[i] = i0_correct;
         f2_correct[i] = f0_correct;
         d2_correct[i] = d0_correct;
         b2_correct[i] = b0_correct;
         c2_correct[i] = c0_correct;
         s2_correct[i] = s0_correct;
         box2_correct[i] = box0_correct;
      }

      int i2[nsize];
      float f2[nsize];
      double d2[nsize];
      bool b2[nsize];
      dcomplex c2[nsize];
      std::string s2[nsize];
      tbox::DatabaseBox box2[nsize];
      basic_array_db->getIntegerArray("i2", i2, nsize);
      basic_array_db->getFloatArray("f2", f2, nsize);
      basic_array_db->getDoubleArray("d2", d2, nsize);
      basic_array_db->getBoolArray("b2", b2, nsize);
      basic_array_db->getComplexArray("c2", c2, nsize);
      basic_array_db->getStringArray("s2", s2, nsize);
      basic_array_db->getDatabaseBoxArray("box2", box2, nsize);

      for (int i = 0; i < nsize; ++i) {
         if (i2[i] != i2_correct[i]) {
            ++fail_count;
            tbox::perr << "Integer test #2 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<float>::equalEps(f2[i], f2_correct[i])) {
            ++fail_count;
            tbox::perr << "Float test #2 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<double>::equalEps(d2[i], d2_correct[i])) {
            ++fail_count;
            tbox::perr << "Double test #2 FAILED" << std::endl;
         }
         if (b2[i] != b2_correct[i]) {
            ++fail_count;
            tbox::perr << "Bool test #2 FAILED" << std::endl;
         }
         if (!tbox::MathUtilities<dcomplex>::equalEps(c2[i], c2_correct[i])) {
            ++fail_count;
            tbox::perr << "Complex test #2 FAILED" << std::endl;
         }
         if (s2[i] != s2_correct[i]) {
            ++fail_count;
            tbox::perr << "String test #2 FAILED" << std::endl;
         }
         if (!(box2[i] == box2_correct[i])) {
            ++fail_count;
            tbox::perr << "Box test #2 FAILED" << std::endl;
         }
      }

      /*******************************************************************
       * Test "getWithDefault()" methods
       ******************************************************************/
      std::shared_ptr<tbox::Database> with_default_db(
         input_db->getDatabase("WithDefaultTypes"));

      int i3 = with_default_db->getIntegerWithDefault("i3", i0_correct);
      float f3 = with_default_db->getFloatWithDefault("f3", f0_correct);
      double d3 = with_default_db->getDoubleWithDefault("d3", d0_correct);
      bool b3 = with_default_db->getBoolWithDefault("b3", b0_correct);
      dcomplex c3 = with_default_db->getComplexWithDefault("c3", c0_correct);
      std::string s3 = with_default_db->getStringWithDefault("s3", s0_correct);
      tbox::DatabaseBox box3 = with_default_db->getDatabaseBoxWithDefault(
            "box3",
            box0_correct);

      if (i3 != i0_correct) {
         ++fail_count;
         tbox::perr << "Integer test #3 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<float>::equalEps(f3, f0_correct)) {
         ++fail_count;
         tbox::perr << "Float test #3 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<double>::equalEps(d3, d0_correct)) {
         ++fail_count;
         tbox::perr << "Double test #3 FAILED" << std::endl;
      }
      if (b3 != b0_correct) {
         ++fail_count;
         tbox::perr << "Bool test #3 FAILED" << std::endl;
      }
      if (!tbox::MathUtilities<dcomplex>::equalEps(c3, c0_correct)) {
         ++fail_count;
         tbox::perr << "Complex test #3 FAILED" << std::endl;
      }
      if (s3 != s0_correct) {
         ++fail_count;
         tbox::perr << "String test #3 FAILED" << std::endl;
      }
      if (!(box3 == box0_correct)) {
         ++fail_count;
         tbox::perr << "Box test #3 FAILED" << std::endl;
      }

      /*******************************************************************
       * Test replacing values in the database
       ******************************************************************/
      std::shared_ptr<tbox::Database> prim_type_db_new(
         input_db->getDatabase("PrimitiveTypes"));

      prim_type_db_new->putInteger("i0", i0_correct);
      prim_type_db_new->putFloat("f0", f0_correct);
      prim_type_db_new->putDouble("d0", d0_correct);
      prim_type_db_new->putBool("b0", b0_correct);
      prim_type_db_new->putComplex("c0", c0_correct);
      prim_type_db_new->putString("s0", s0_correct);
      prim_type_db_new->putDatabaseBox("box0", box0_correct);

      /*******************************************************************
       * Output contents of the input database
       ******************************************************************/
      tbox::plog << "Overall contents of input database..." << std::endl;
      input_db->printClassData(tbox::plog);

      tbox::plog << "\n\nUnused keys in the input database..." << std::endl;
      input_db->printUnusedKeys(tbox::plog);

      tbox::plog << "\n\nDefault keys in the input database..." << std::endl;
      input_db->printDefaultKeys(tbox::plog);

      input_db.reset();
   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  inputdb" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
