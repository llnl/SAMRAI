/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test variable database operations
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <string>
#include <memory>

using namespace SAMRAI;

int main(
   int argc,
   char* argv[]) {

   tbox::Dimension dim(2);

   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

// tbox::PIO::logOnlyNodeZero("main-var_db.log");
      tbox::PIO::logAllNodes("vdbtest.log");

      hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

      std::shared_ptr<hier::VariableContext> current_context(
         var_db->getContext("CURRENT"));

      hier::IntVector nghosts(dim, 4);
      hier::IntVector fluxghosts(dim, 1);
      hier::IntVector zero_ghosts(dim, 0);

      /* State variable */
      std::shared_ptr<pdat::CellVariable<double> > uval(
         new pdat::CellVariable<double>(dim, "uval", 1));

      /* Flux variable */
      std::shared_ptr<pdat::FaceVariable<double> > flux(
         new pdat::FaceVariable<double>(dim, "flux", 1));

      /* Register uval using ready made context
       * and by creating a new one on the fly
       */
      const int uval_current_id =
         var_db->registerVariableAndContext(uval,
            current_context,
            zero_ghosts);

      /* Register flux using pre-existing context ("CURRENT"), and by using
       * context that does not yet exist in the Database ("SCRATCH").
       */
// use of void eliminates compiler warning
// const int flux_current_id =
      (void)
      var_db->registerVariableAndContext(flux,
         current_context,
         zero_ghosts);

      const int flux_scratch_id =
         var_db->registerVariableAndContext(flux,
            var_db->getContext("SCRATCH"),
            nghosts);

      std::shared_ptr<pdat::OuterfaceVariable<double> > fluxsum(
         new pdat::OuterfaceVariable<double>(
            dim, "fluxsum", 1));

// use of void eliminates compiler warning
// const int fluxsum_current_id =
      (void)
      var_db->registerVariableAndContext(fluxsum,
         current_context,
         zero_ghosts);

      /*
       * Now, run the variable database through the ringer...
       */

      tbox::plog
      << "\n\nPrintout #1 of hier::Variable tbox::Database (after initial registration)..."
      << std::endl;
      var_db->printClassData(tbox::plog);

      // Test #1: Check Context functions...

      // Test #1a: hier::VariableDatabase::checkContextExists()
      tbox::plog
      << "Test #1a: hier::VariableDatabase::checkContextExists()..." << std::endl;
      if (!var_db->checkContextExists("SCRATCH")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #1a: hier::VariableDatabase::checkContextExists()\n"
         << "SCRATCH context added to var_db, but not found" << std::endl;
      }

      // Test #1b: hier::VariableDatabase::checkContextExists()
      tbox::plog
      << "Test #1b: hier::VariableDatabase::checkContextExists()..." << std::endl;
      if (!var_db->checkContextExists("CURRENT")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #1b: hier::VariableDatabase::checkContextExists()\n"
         << "CURRENT context added to var_db, but not found" << std::endl;
      }

      // Test #1c: hier::VariableDatabase::checkContextExists()
      tbox::plog
      << "Test #1c: hier::VariableDatabase::checkContextExists()..." << std::endl;
      if (var_db->checkContextExists("dummy")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #1c: hier::VariableDatabase::checkContextExists()\n"
         << "dummy context not added to var_db, but found" << std::endl;
      }

      // Adding dummy context to variable databse...
      /*
       * Although the dummy_ctxt is unused, we are checking for it in
       * the test.  So leave it in despite possible compiler warnings.
       */
      std::shared_ptr<hier::VariableContext> dummy_ctxt(
         var_db->getContext("dummy"));
      NULL_USE(dummy_ctxt);

      // Test #1d: hier::VariableDatabase::checkContextExists()
      tbox::plog
      << "Test #1d: hier::VariableDatabase::checkContextExists()..." << std::endl;
      if (!var_db->checkContextExists("dummy")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #1d: hier::VariableDatabase::checkContextExists()\n"
         << "dummy context added to var_db, but not found" << std::endl;
      }

      // Test #2,3: Check hier::Variable functions...

      // Test #2a: hier::VariableDatabase::checkVariableExists()
      tbox::plog
      << "Test #2a: hier::VariableDatabase::checkVariableExists()..." << std::endl;
      if (!var_db->checkVariableExists("uval")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #2a: hier::VariableDatabase::checkVariableExists()\n"
         << "uval variable added to var_db, but not found" << std::endl;
      }

      // Test #2b: hier::VariableDatabase::checkVariableExists()
      tbox::plog
      << "Test #2b: hier::VariableDatabase::checkVariableExists()..." << std::endl;
      if (!var_db->checkVariableExists("flux")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #2b: hier::VariableDatabase::checkVariableExists()\n"
         << "flux variable added to var_db, but not found" << std::endl;
      }

      // Test #2c: hier::VariableDatabase::checkVariableExists()
      tbox::plog
      << "Test #2c: hier::VariableDatabase::checkVariableExists()..." << std::endl;
      if (!var_db->checkVariableExists("fluxsum")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #2c: hier::VariableDatabase::checkVariableExists()\n"
         << "fluxsum variable added to var_db, but not found"
         << std::endl;
      }

      // Test #2d: hier::VariableDatabase::checkVariableExists()
      tbox::plog
      << "Test #2d: hier::VariableDatabase::checkVariableExists()..." << std::endl;
      if (var_db->checkVariableExists("dummy")) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #2d: hier::VariableDatabase::checkVariableExists()\n"
         << "dummy variable not added to var_db, but found" << std::endl;
      }

      // Test #3a: hier::VariableDatabase::getVariable()
      tbox::plog << "Test #3a: hier::VariableDatabase::getVariable()..."
                 << std::endl;
      std::shared_ptr<hier::Variable> tvar_uval(var_db->getVariable("uval"));
      if (!tvar_uval) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #3a: hier::VariableDatabase::getVariable()\n"
         << "uval variable added to var_db, but returning NULL"
         << std::endl;
      }

      // Test #3b: hier::VariableDatabase::getVariable()
      tbox::plog << "Test #3b: hier::VariableDatabase::getVariable()..."
                 << std::endl;
      std::shared_ptr<hier::Variable> tvar_flux(var_db->getVariable("flux"));
      if (!tvar_flux) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #3b: hier::VariableDatabase::getVariable()\n"
         << "flux variable added to var_db, but returning NULL"
         << std::endl;
      }

      // Test #3c: hier::VariableDatabase::getVariable()
      tbox::plog << "Test #3c: hier::VariableDatabase::getVariable()..."
                 << std::endl;
      std::shared_ptr<hier::Variable> tvar_fluxsum(
         var_db->getVariable("fluxsum"));
      if (!tvar_fluxsum) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #3c: hier::VariableDatabase::getVariable()\n"
         << "fluxsum variable added to var_db, but returning NULL"
         << std::endl;
      }

      // Test #3d: hier::VariableDatabase::getVariable()
      tbox::plog << "Test #3d: hier::VariableDatabase::getVariable()..."
                 << std::endl;
      //   tbox::perr << "Attempt to get variable named dummy..." << std::endl;
      std::shared_ptr<hier::Variable> tvar_dummy(
         var_db->getVariable("dummy"));
      if (tvar_dummy) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #3d: hier::VariableDatabase::getVariable()\n"
         << "dummy variable not added to var_db, but not returning NULL"
         << std::endl;
      }

      // Test #4: Check instance identifier assignments

      // Test #4a: hier::Variable::getInstanceIdentifier()
      tbox::plog << "Test #4a: hier::Variable::getInstanceIdentifier()..."
                 << std::endl;
      int uval_id = tvar_uval->getInstanceIdentifier();
      if (uval_id != 0) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #4a: hier::Variable::getInstanceIdentifier()\n"
         << "uval should have id = 0, but has id = " << uval_id
         << std::endl;
      }

      // Test #4b: hier::Variable::getInstanceIdentifier()
      tbox::plog << "Test #4b: hier::Variable::getInstanceIdentifier()..."
                 << std::endl;
      int flux_id = tvar_flux->getInstanceIdentifier();
      if (flux_id != 1) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #4b: hier::Variable::getInstanceIdentifier()\n"
         << "flux should have id = 1, but has id = " << flux_id
         << std::endl;
      }

      // Test #4c: hier::Variable::getInstanceIdentifier()
      tbox::plog << "Test #4c: hier::Variable::getInstanceIdentifier()..."
                 << std::endl;
      int fluxsum_id = tvar_fluxsum->getInstanceIdentifier();
      if (fluxsum_id != 2) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #4c: hier::Variable::getInstanceIdentifier()\n"
         << "fluxsum should have id = 2, but has id = "
         << fluxsum_id << std::endl;
      }

      // Test #5: Attempt to register (uval,CURRENT) again
      tbox::plog << "Test #5: Attempt to register (uval,CURRENT) again..."
                 << std::endl;
      std::shared_ptr<hier::VariableContext> tctxt_current(
         var_db->getContext("CURRENT"));
      hier::IntVector tzero_ghosts(dim, 0);
      int ti = var_db->registerVariableAndContext(
            tvar_uval, tctxt_current, tzero_ghosts);
      if (ti != 0) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #5: Re-registering a variable and context\n"
         << "Original id = 0 should be returned, but "
         << ti << "returned" << std::endl;
      }

      // Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()
      tbox::plog
      << "Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()..."
      << std::endl;
      ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_current);
      if (ti != 0) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #6a: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
         << "(uval,CURRENT) should be mapped to 0, but is mapped to "
         << ti << std::endl;
      }

      // Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()
      tbox::plog
      << "Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()..."
      << std::endl;
      tvar_uval = var_db->getVariable("uval");
      std::shared_ptr<hier::VariableContext> tctxt_scratch(
         var_db->getContext("SCRATCH"));
      ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_scratch);
      if (ti != -1) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #6b: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
         << "(uval,SCRATCH) should be mapped to -1, but is mapped to "
         << ti << std::endl;
      }

      // Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()
      tbox::plog
      << "Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()..."
      << std::endl;
      std::shared_ptr<pdat::CellVariable<double> > dummy_var(
         new pdat::CellVariable<double>(dim, "dummy", 3));
      tctxt_scratch = var_db->getContext("SCRATCH");
      ti = var_db->mapVariableAndContextToIndex(dummy_var, tctxt_scratch);
      if (ti != -1) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #6c: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
         << "(dummy,SCRATCH) should be mapped to -1, but is mapped to "
         << ti << std::endl;
      }

      // Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()
      tbox::plog
      << "Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()..."
      << std::endl;
      tvar_uval = var_db->getVariable("uval");
      std::shared_ptr<hier::VariableContext> tctxt_random(
         new hier::VariableContext("RANDOM"));
      ti = var_db->mapVariableAndContextToIndex(tvar_uval, tctxt_random);
      if (ti != -1) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #6d: hier::VariableDatabase::mapVariableAndContextToIndex()\n"
         << "(uval,RANDOM) should be mapped to -1, but is mapped to "
         << ti << std::endl;
      }

      // Test #7a: hier::VariableDatabase::mapIndexToVariableAndContext()
      tbox::plog
      << "Test #7a: hier::VariableDatabase::mapIndexToVariableAndContext()..."
      << std::endl;
      int search_id = 2;
      std::shared_ptr<hier::Variable> search_var;
      std::shared_ptr<hier::VariableContext> search_ctxt;
      std::string flux_variable("flux");
      std::string scratch_variable("SCRATCH");

      // searching for index = 2
      if (!var_db->mapIndexToVariableAndContext(
             search_id, search_var, search_ctxt)) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #7a.1: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
         << "Problem finding a (variable,context) pair for index = 2"
         << std::endl;
         if (search_var->getName() != flux_variable) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7a.2: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
            << "Returned var name should be \"flux\" but is "
            << search_var->getName() << std::endl;
         }
         if (search_ctxt->getName() != scratch_variable) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7a.3: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
            << "Returned context name should be \"SCRATCH\" but is "
            << search_ctxt->getName() << std::endl;
         }
/*      tbox::plog << "Var name = " << tvar->getName() << " = flux?, "
 *      << "Context name = " << tctxt->getName() << " = CURRENT?" << std::endl;
 */
/*   } else {
 *   tbox::plog << "Houston, we have a problem looking for index "
 *   << search_id << std::endl;
 */
      }

      // Test #7b: hier::VariableDatabase::mapIndexToVariableAndContext()
      tbox::plog
      << "Test #7b: hier::VariableDatabase::mapIndexToVariableAndContext()"
      << std::endl;
      search_id = 20;

      // searching for index = 20
      if (var_db->mapIndexToVariableAndContext(
             search_id, search_var, search_ctxt)) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #7b.1: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
         << "Something maps to index = 20 when nothing should.\n"
         << "Variable name: " << search_var->getName() << "\n"
         << "Context name: " << search_ctxt->getName()
         << std::endl;
      } else {

         if (search_var) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7b.2: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
            << "search_var should be NULL" << std::endl;
         }
         if (search_ctxt) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7b.3: hier::VariableDatabase::mapIndexToVariableAndContext()\n"
            << "search_ctxt should be NULL" << std::endl;
         }

      }

      // Test #7c: hier::VariableDatabase::mapIndexToVariable()
      tbox::plog
      << "Test #7c: hier::VariableDatabase::mapIndexToVariable()..." << std::endl;
      search_id = 2;

      // searching for index = 2
      if (!var_db->mapIndexToVariable(search_id, search_var)) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #7c.1: hier::VariableDatabase::mapIndexToVariable()\n"
         << "Problem finding a (variable) for index = 2"
         << std::endl;
         if (search_var->getName() != flux_variable) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7c.2: hier::VariableDatabase::mapIndexToVariable()\n"
            << "Returned var name should be \"flux\" but is "
            << search_var->getName() << std::endl;
         }
      }

      // Test #7d: hier::VariableDatabase::mapIndexToVariable()
      tbox::plog
      << "Test #7d: hier::VariableDatabase::mapIndexToVariable()..." << std::endl;
      search_id = 20;

      // searching for index = 20
      if (var_db->mapIndexToVariable(search_id, search_var)) {
         ++fail_count;
         tbox::perr
         << "FAILED: - Test #7d.1: hier::VariableDatabase::mapIndexToVariable()\n"
         << "Something maps to index = 20 when nothing should.\n"
         << "Variable name: " << search_var->getName() << "\n"
         << std::endl;
      } else {

         if (search_var) {
            ++fail_count;
            tbox::perr
            << "FAILED: - Test #7d.2: hier::VariableDatabase::mapIndexToVariable()\n"
            << "search_var should be NULL" << std::endl;
         }

      }

      // Testing new registration functions and ability of variable
      // database to add and delete descriptor indices...

      tbox::plog
      << "\nPrintout #2 of hier::Variable tbox::Database (after tests 1-7."
      << "\n Only difference with previous printout should be addition"
      << "\n of \"dummy\" context).....\n";
      var_db->printClassData(tbox::plog);

      // Test #8a: Checking mapping in variable database
      tbox::plog << "Test #8a: Checking mapping in variable database..."
                 << std::endl;
      if (!var_db->checkVariablePatchDataIndex(uval, uval_current_id)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8a: "
                    << "VariableDatabase:checkVariablePatchDataIndex()\n"
                    << "uval should be map to current id in patch descriptor"
                    << std::endl;
      }
      // Test #8b: Checking mapping in variable database
      tbox::plog << "Test #8b: Checking mapping in variable database..."
                 << std::endl;
      if (!var_db->checkVariablePatchDataIndexType(uval, uval_current_id)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8b: "
                    << "VariableDatabase:checkVariablePatchDataIndexType()\n"
                    << "uval should be map to current id in patch descriptor"
                    << std::endl;
      }

      // Test #8c: Checking restration of existing patch data index...
      tbox::plog
      << "Test #8c: Checking restration of existing patch data index..."
      << std::endl;
      int test_id = var_db->registerPatchDataIndex(uval, uval_current_id);

      if (test_id != uval_current_id) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8c: "
                    << "VariableDatabase:registerPatchDataIndex()\n"
                    << "re-registering current uval id should return same id"
                    << std::endl;
      }

      // Test #8d: Testing registration of new cloned factory to variable...
      tbox::plog
      << "Test #8d: Testing registration of new cloned factory to variable..."
      << std::endl;
      int new_id =
         var_db->registerClonedPatchDataIndex(uval, uval_current_id);

      if ((new_id < 0) || (new_id == uval_current_id)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8d: "
                    << "VariableDatabase:registerClonedPatchDataIndex()\n"
                    << "cloning current uval id return invalid id" << std::endl;
      }

      // Test #8e: Testing registration of new cloned factory to variable...
      tbox::plog
      << "Test #8e: Testing registration of new cloned factory to variable..."
      << std::endl;
      std::shared_ptr<hier::Variable> tvar;
      if (!var_db->mapIndexToVariable(new_id, tvar)
          || (tvar->getName() != "uval")) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8e: "
                    << "VariableDatabase:mapIndexToVariable()\n"
                    << "descriptor id = " << new_id
                    << " should now map to uval" << std::endl;
      }

      tbox::plog
      << "\nPrintout #3 of hier::Variable tbox::Database. (after tests 8a-8e"
      << "\n Descriptor index " << new_id
      << " should map to newly-created uval descriptor index)"
      << std::endl;
      var_db->printClassData(tbox::plog);

      // Test #8f: Testing removal of new variable-descriptor id mapping
      var_db->removePatchDataIndex(new_id);

      tbox::plog
      << "\nPrintout #4 of hier::Variable tbox::Database. (after test 8d"
      << "\n Descriptor index " << new_id
      << " should no longer be mapped to a variable)..." << std::endl;
      var_db->printClassData(tbox::plog);

      tvar.reset();
      if (var_db->mapIndexToVariable(new_id, tvar)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8f: "
                    << "VariableDatabase:removePatchDataIndex()\n"
                    << "descriptor id = " << new_id
                    << " should no longer map to uval variable" << std::endl;
      }

      // Test #8g-h: Testing whether inconsistent mapping is allowed...
      tbox::plog
      << "Test #8g-h: Testing whether inconsistent mapping is allowed..."
      << std::endl;
      if (!var_db->checkVariablePatchDataIndex(flux, flux_scratch_id)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8g: "
                    << "VariableDatabase:checkVariablePatchDataIndex()\n"
                    << "flux should be mapped to scratch flux id" << std::endl;
      }
      if (var_db->checkVariablePatchDataIndex(uval, flux_scratch_id)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8h: "
                    << "VariableDatabase:checkVariablePatchDataIndex()\n"
                    << "uval should not map to scratch flux id" << std::endl;
      }

      // Test #8i-j: Testing removal of existing descriptor id...
      tbox::plog
      << "Test #8i-j: Testing removal of existing descriptor id..." << std::endl;
      var_db->removePatchDataIndex(flux_scratch_id);

      int tindex =
         var_db->mapVariableAndContextToIndex(flux,
            var_db->getContext("SCRATCH"));
      if (tindex > -1) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8i: "
                    << "VariableDatabase:removePatchDataIndex()\n"
                    << "flux-SCRATCH mapping should no longer be in database"
                    << std::endl;
      }

      tvar.reset();
      if (var_db->mapIndexToVariable(flux_scratch_id, tvar)) {
         ++fail_count;
         tbox::perr << "FAILED: - Test #8j: "
                    << "VariableDatabase:mapIndexToVariable()\n"
                    << "descriptor id = " << flux_scratch_id
                    << " should no longer map to flux variable" << std::endl;
      }

      tbox::plog
      << "\nPrintout #5 of hier::Variable tbox::Database. (after tests 8e-j\n"
      << "Should be identical to first printout except for the "
      << "\naddition of \n\"dummy\" context and removal of "
      << "flux scratch variable-context mapping" << std::endl;
      var_db->printClassData(tbox::plog);

      tbox::plog
      << "\nPrintout #6 of hier::Variable tbox::Database after removal of \"flux\" "
      << "variable.";
      var_db->removeVariable("flux");
      var_db->printClassData(tbox::plog);

      /*
       * Tests that end in program abort...
       */
#if 0
      // Abort Test #1
      std::shared_ptr<pdat::CellVariable<double> > dummy(
         new pdat::CellVariable<double>("uval", 2));

      tbox::plog << "Attempt to add a different variable named uval."
                 << "This should bomb!!" << std::endl;
      var_db->addVariable(dummy);

      // Abort Test #2
      tbox::plog << "Attempt to register uval, CURRENT again w/ wrong ghosts."
                 << "This should bomb!!" << std::endl;
      tctxt = var_db->getContext("CURRENT");
      tvar = var_db->getVariable("uval");
      g = hier::IntVector(2);
      ti = var_db->registerVariableAndContext(tvar, tctxt, g);
      tbox::plog << "uval, CURRENT at index = " << ti << std::endl;

      // Abort Test #3
      tbox::plog << "Attempt to register uval with fake CURRENT context."
                 << "This should bomb!!" << std::endl;
      tctxt = new hier::VariableContext("CURRENT");
      tvar = var_db->getVariable("uval");
      g = hier::IntVector(0);
      ti = var_db->registerVariableAndContext(tvar, tctxt, g);
      tbox::plog << "uval, fake CURRENT at index = " << ti << std::endl;

      // Abort Test #4
      tbox::plog << "Attempt to map uval to descriptor id for flux."
                 << "This should bomb!!" << std::endl;
      var_db->addVariablePatchDataIndex(uval, flux_scratch_id);
#endif

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  vdbtest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
