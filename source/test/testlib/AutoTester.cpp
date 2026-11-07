/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class used for auto testing applications
 *
 ************************************************************************/

#include "AutoTester.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/FlattenedHierarchy.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/MathUtilities.h"

AutoTester::AutoTester(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> input_db):
   d_dim(dim),
   d_base_name("unnamed"),
   d_baseline_data("AutoTesterDatabase")
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_object_name = object_name;
   d_test_fluxes = false;
   d_test_iter_num = 10;
   d_output_correct = false;

   d_write_patch_boxes = false;
   d_read_patch_boxes = false;
   d_test_patch_boxes_at_steps.resize(0);
   d_test_patch_boxes_step_count = 0;

   getFromInput(input_db);

   std::string test_patch_boxes_filename = "test_inputs/";
#if defined(__xlC__)
#ifdef OPT_BUILD
   test_patch_boxes_filename += "xlC/";
#else
   test_patch_boxes_filename += "xlC_debug/";
#endif
#endif
   test_patch_boxes_filename += d_base_name + ".boxes";

   const std::string baseline_filename =
      test_patch_boxes_filename
      + "." + tbox::Utilities::nodeToString(mpi.getSize())
      + "." + tbox::Utilities::processorToString(mpi.getRank());

   if (d_read_patch_boxes) {
      if (!d_baseline_data.open(baseline_filename)) {
         TBOX_WARNING("Could not open baseline " << baseline_filename
            << ".baseline.  Skipping regression comparison.\n");
         d_read_patch_boxes = false;
      }
      if (d_read_patch_boxes && d_output_correct) {
         d_baseline_data.getDatabase()->printClassData(tbox::pout);
      }

   }

   if (d_write_patch_boxes) {
      d_baseline_data.create(baseline_filename);
   }

}

AutoTester::~AutoTester()
{
   d_baseline_data.close();
}

/*
 ******************************************************************
 *
 *  Method "evalTestData" compares the result of the run with
 *  the correct result for runs with the TimeRefinementIntegrator
 *  and HyperbolicLevelIntegrator.
 *
 ******************************************************************
 */
int AutoTester::evalTestData(
   int iter,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const std::shared_ptr<algs::TimeRefinementIntegrator> tri,
   const std::shared_ptr<algs::HyperbolicLevelIntegrator> hli,
   const std::shared_ptr<mesh::GriddingAlgorithm> ga)
{
   NULL_USE(ga);

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) {

      /*
       * set precision of output stream.
       */
      tbox::plog.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      std::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(nlevels));

      /*
       * Test 0: Time Refinement Integrator
       */
      double time = tri->getIntegratorTime();
      if (d_correct_result.size() > 0) {
         if (d_output_correct) {
            tbox::plog << "Test 0: Time Refinement Integrator "
                       << "\n   computed result: " << time;

            tbox::plog << "\n   specified result = "
                       << d_correct_result[0];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(time,
                d_correct_result[0])) {
            tbox::plog << "Test 0: Time Refinement check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 0 FAILED: Check Time Refinement Integrator"
                       << std::endl;
            ++num_failures;
         }
      }

      /*
       * Test 1: Time Refinement Integrator
       */
      double dt = tri->getLevelDtMax(nlevels);
      if (d_correct_result.size() > 1) {
         if (d_output_correct) {
            tbox::plog << "Test 1: Time Refinement Integrator "
                       << "\n   computed result: " << dt;
            tbox::plog << "\n   specified result = "
                       << d_correct_result[1];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[1])) {
            tbox::plog << "Test 1: Time Refinement check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 1 FAILED: Check Time Refinement Integrator"
                       << std::endl;
            ++num_failures;
         }
      }

      /*
       * Test 2: Hyperbolic Level Integrator
       */
      dt = hli->getLevelDt(level, time, false);
      if (d_correct_result.size() > 2) {
         if (d_output_correct) {
            tbox::plog << "Test 2: Hyperbolic Level Integrator "
                       << "\n   computed result: " << dt;

            tbox::plog << "\n   specified result = "
                       << d_correct_result[2];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[2])) {
            tbox::plog << "Test 2: Hyperbolic Level Int check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 2 FAILED: Check Hyperbolic Level Integrator"
                       << std::endl;
            ++num_failures;
         }
      }

      /*
       * Test 3: Gridding Algorithm
       */
      int n = hierarchy->getMaxNumberOfLevels();
      if (d_output_correct) {
         tbox::plog << "Test 3: Gridding Algorithm "
                    << "\n   computed result: " << n;
         tbox::plog << "\n   correct result = " << nlevels + 1;
         tbox::plog << std::endl;
      }
      if (n == (nlevels + 1)) {
         tbox::plog << "Test 3: Gridding Algorithm check successful"
                    << std::endl;
      } else {
         tbox::perr << "Test 3 FAILED: Check Gridding Algorithm" << std::endl;
         ++num_failures;
      }

      num_failures += testHierarchyNeighbors(hierarchy);
      num_failures += testFlattenedHierarchy(hierarchy);

   }

   if ((static_cast<int>(d_test_patch_boxes_at_steps.size()) >
        d_test_patch_boxes_step_count) &&
       (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] == iter)) {

      int num_levels = hierarchy->getNumberOfLevels();

      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_baseline_data.getDatabase()->printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         std::shared_ptr<tbox::Database> step_db(
            d_baseline_data.getDatabase()->getDatabase(step_name));

         /*
          * TODO: This check give false positives!!!!!
          * It writes the same file regardless of the number of processors.
          * We should be checking against base runs with the same number of processors,
          * compare different data.
          */
         for (int ln = 0; ln < num_levels; ++ln) {

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            std::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
               hierarchy->getGridGeometry());
            std::shared_ptr<tbox::Database> level_db(
               step_db->getDatabase(level_name));
            hier::BoxLevel correct_box_level(d_dim, *level_db, grid_geometry);

            num_failures += checkHierarchyBoxes(hierarchy,
                  ln,
                  correct_box_level,
                  iter);
         }

      }

      if (d_write_patch_boxes) {

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         std::shared_ptr<tbox::Database> step_db(
            d_baseline_data.getDatabase()->putDatabase(step_name));

         for (int ln = 0; ln < num_levels; ++ln) {
            std::shared_ptr<hier::PatchLevel> level(
               hierarchy->getPatchLevel(ln));

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            std::shared_ptr<tbox::Database> level_db(
               step_db->putDatabase(level_name));
            level->getBoxLevel()->putToRestart(level_db);
         }

         if (d_output_correct) {
            d_baseline_data.getDatabase()->printClassData(tbox::pout);
         }
      }

      ++d_test_patch_boxes_step_count;

   }

   return num_failures;
}

/*
 ******************************************************************
 *
 *  Method "evalTestData" compares the result of the run with
 *  the correct result for runs with the MethodOfLinesIntegrator.
 *
 ******************************************************************
 */
int AutoTester::evalTestData(
   int iter,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   double time,
   const std::shared_ptr<algs::MethodOfLinesIntegrator> mol,
   const std::shared_ptr<mesh::GriddingAlgorithm> ga)
{
   NULL_USE(ga);

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) {

      /*
       * set precision of output stream.
       */
      tbox::plog.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      std::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(nlevels));

      /*
       * Test 0: Time test
       */
      if (d_correct_result.size() > 0) {
         if (d_output_correct) {
            tbox::plog << "Test 0: Simulation Time: "
                       << "\n   computed result: " << time;
            tbox::plog << "\n   specified result = "
                       << d_correct_result[0];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(time,
                d_correct_result[0])) {
            tbox::plog << "Test 0: Simulation Time check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 0 FAILED: Simulation time incorrect"
                       << std::endl;
            ++num_failures;
         }
      }

      /*
       * Test 1: MethodOfLinesIntegrator
       */
      double dt = mol->getTimestep(hierarchy, time);
      if (d_correct_result.size() > 1) {
         if (d_output_correct) {
            tbox::plog << "Test 1: Method of Lines Integrator "
                       << "\n   computed result: " << dt;
            tbox::plog << "\n   specified result = "
                       << d_correct_result[1];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[1])) {
            tbox::plog << "Test 1: MOL Int check successful" << std::endl;
         } else {
            tbox::perr << "Test 1 FAILED: Check Method of Lines Integrator"
                       << std::endl;
            ++num_failures;
         }
      }

      /*
       * Test 2: Gridding Algorithm
       */
      int n = hierarchy->getMaxNumberOfLevels();
      if (d_output_correct) {
         tbox::plog << "Test 2: Gridding Algorithm "
                    << "\n   computed result: " << n;
         tbox::plog << "\n   correct result = " << nlevels + 1;
         tbox::plog << std::endl;
      }
      if (n == (nlevels + 1)) {
         tbox::plog << "Test 2: Gridding Alg check successful" << std::endl;
      } else {
         tbox::perr << "Test 2 FAILED: Check Gridding Algorithm" << std::endl;
         ++num_failures;
      }

      num_failures += testHierarchyNeighbors(hierarchy);
      num_failures += testFlattenedHierarchy(hierarchy);

   }

   if ((static_cast<int>(d_test_patch_boxes_at_steps.size()) > 0) &&
       (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] == iter)) {

      int num_levels = hierarchy->getNumberOfLevels();

      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_baseline_data.getDatabase()->printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         std::shared_ptr<tbox::Database> step_db(
            d_baseline_data.getDatabase()->getDatabase(step_name));

         for (int ln = 0; ln < num_levels; ++ln) {

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            std::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
               hierarchy->getGridGeometry());
            std::shared_ptr<tbox::Database> level_db(
               step_db->getDatabase(level_name));
            hier::BoxLevel correct_box_level(d_dim, *level_db, grid_geometry);

            num_failures += checkHierarchyBoxes(hierarchy,
                  ln,
                  correct_box_level,
                  iter);
         }

      }

      if (d_write_patch_boxes) {

         if (d_output_correct) {
            d_baseline_data.getDatabase()->printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         std::shared_ptr<tbox::Database> step_db(
            d_baseline_data.getDatabase()->putDatabase(step_name));

         for (int ln = 0; ln < num_levels; ++ln) {
            std::shared_ptr<hier::PatchLevel> level(
               hierarchy->getPatchLevel(ln));

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            std::shared_ptr<tbox::Database> level_db(
               step_db->putDatabase(level_name));
            level->getBoxLevel()->putToRestart(level_db);
         }

      }

      ++d_test_patch_boxes_step_count;

   }

   return num_failures;
}

int AutoTester::testHierarchyNeighbors(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   int num_failures = 0;
   int num_levels = hierarchy->getNumberOfLevels();
   hier::HierarchyNeighbors hier_nbrs(*hierarchy, 0, num_levels-1);
   for (int ln = 0; ln < num_levels; ++ln) {
      const std::shared_ptr<hier::PatchLevel>& current_level =
         hierarchy->getPatchLevel(ln);

      if (ln < num_levels-1) {

         std::shared_ptr<hier::PatchLevel> finer_level(
            hierarchy->getPatchLevel(ln+1));

         const hier::Connector& coarse_to_fine =
            current_level->findConnector(
               *finer_level,
               hier::IntVector::getOne(hierarchy->getDim()),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         for (hier::Connector::ConstNeighborhoodIterator cf =
              coarse_to_fine.begin(); cf != coarse_to_fine.end(); ++cf) {

            const hier::BoxId& crse_box_id(*cf);
            const hier::Box& crse_box =
               *(current_level->getBoxLevel()->getBox(crse_box_id));
            const hier::BoxContainer& finer_nbrs =
               hier_nbrs.getFinerLevelNeighbors(crse_box, ln);
            TBOX_ASSERT(finer_nbrs.isOrdered());

            for (hier::Connector::ConstNeighborIterator ni =
                 coarse_to_fine.begin(cf);
                 ni != coarse_to_fine.end(cf); ++ni) {

               const hier::Box& fine_box = *ni;
               if (finer_nbrs.find(fine_box) == finer_nbrs.end()) {
                  tbox::perr << "Test fine FAILED" << std::endl;
                  ++num_failures;
               }
            }
         }
      }

      if (ln > 0) {

         std::shared_ptr<hier::PatchLevel> coarser_level(
            hierarchy->getPatchLevel(ln-1));

         const hier::Connector& fine_to_coarse =
            current_level->findConnector(
               *coarser_level,
               hier::IntVector::getOne(hierarchy->getDim()),
               hier::CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         for (hier::Connector::ConstNeighborhoodIterator fc =
              fine_to_coarse.begin(); fc != fine_to_coarse.end(); ++fc) {

            const hier::BoxId& fine_box_id(*fc);
            const hier::Box& fine_box =
               *(current_level->getBoxLevel()->getBox(fine_box_id));
            const hier::BoxContainer& coarser_nbrs =
               hier_nbrs.getCoarserLevelNeighbors(fine_box, ln);
            TBOX_ASSERT(coarser_nbrs.isOrdered());

            for (hier::Connector::ConstNeighborIterator ni =
                 fine_to_coarse.begin(fc);
                 ni != fine_to_coarse.end(fc); ++ni) {

               const hier::Box& coarse_box = *ni;
               if (coarser_nbrs.find(coarse_box) == coarser_nbrs.end()) {
                  tbox::perr << "Test coarse FAILED" << std::endl;
                  ++num_failures;
               }
            }
         }
      }

      const hier::Connector& current_to_current =
         current_level->findConnector(
            *current_level,
            hier::IntVector::getOne(hierarchy->getDim()),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      for (hier::Connector::ConstNeighborhoodIterator cc =
           current_to_current.begin(); cc != current_to_current.end(); ++cc) {

         const hier::BoxId& current_box_id(*cc);
         const hier::Box& current_box =
            *(current_level->getBoxLevel()->getBox(current_box_id));
         const hier::BoxContainer& same_nbrs =
            hier_nbrs.getSameLevelNeighbors(current_box, ln);
         TBOX_ASSERT(same_nbrs.isOrdered() || same_nbrs.empty());

         for (hier::Connector::ConstNeighborIterator ni =
              current_to_current.begin(cc);
              ni != current_to_current.end(cc); ++ni) {

            const hier::Box& nbr_box = *ni;
            if (nbr_box.getBoxId() != current_box.getBoxId()) {
               if (same_nbrs.find(nbr_box) == same_nbrs.end()) {
                  tbox::perr << "Test same FAILED" << std::endl;
                  ++num_failures;
               }
            }
         }
      }
   }

   return num_failures;
}

int AutoTester::testFlattenedHierarchy(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   int num_failures = 0;
   int num_levels = hierarchy->getNumberOfLevels();
   hier::FlattenedHierarchy flat_hier(*hierarchy, 0, num_levels-1);

   const std::shared_ptr<hier::PatchLevel>& level_zero =
      hierarchy->getPatchLevel(0);

   double level_zero_size =
      static_cast<double>(level_zero->getBoxLevel()->getGlobalNumberOfCells());

   double local_size = 0.0;

   hier::UncoveredBoxIterator itr = flat_hier.beginUncovered();
   hier::UncoveredBoxIterator flat_end = flat_hier.endUncovered();
   for ( ; itr != flat_end; ++itr) {

      const hier::IntVector& ratio_to_zero =
         itr->first->getPatchGeometry()->getRatio();

      const hier::BlockId& block_id = itr->second.getBlockId();
      double refine_quotient =
         static_cast<double>(ratio_to_zero.getProduct(block_id));

      double cell_value = 1.0 / refine_quotient;
      local_size += (cell_value * static_cast<double>(itr->second.size()));

   }

   double global_flat_size = local_size;
   if (hierarchy->getMPI().AllReduce(&global_flat_size, 1, MPI_SUM) != MPI_SUCCESS) {
      tbox::perr << "FAILED: - AutoTester " << "\n"
                 << "MPI sum reduction failed" << std::endl;
      num_failures++;
   }

   if (tbox::MathUtilities<double>::Abs(global_flat_size-level_zero_size) >
       1.0e-8) {
      tbox::perr << "FAILED: - AutoTester " << "\n"
                 << "Flattened hierarchy size not equivalent \n"
                 << "to level zero size." << std::endl;
      num_failures++;
   }

   return num_failures;
}

/*
 ******************************************************************
 *
 *  Get test parameters from input.
 *
 ******************************************************************
 */

void AutoTester::getFromInput(
   std::shared_ptr<tbox::Database> input_db)
{
   std::shared_ptr<tbox::Database> tester_db(
      input_db->getDatabase(d_object_name));
   std::shared_ptr<tbox::Database> main_db(
      input_db->getDatabase("Main"));

   /*
    * Read testing parameters from testing_db
    */
   if (tester_db->keyExists("test_fluxes")) {
      d_test_fluxes = tester_db->getBool("test_fluxes");
   }

   if (tester_db->keyExists("test_iter_num")) {
      d_test_iter_num = tester_db->getInteger("test_iter_num");
   }

   if (tester_db->keyExists("write_patch_boxes")) {
      d_write_patch_boxes = tester_db->getBool("write_patch_boxes");
   }
   if (tester_db->keyExists("read_patch_boxes")) {
      d_read_patch_boxes = tester_db->getBool("read_patch_boxes");
   }
   if (d_read_patch_boxes && d_write_patch_boxes) {
      tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
                 << "Cannot 'read_patch_boxes' and 'write_patch_boxes' \n"
                 << "at the same time." << std::endl;
   }
   d_base_name = main_db->getStringWithDefault("base_name", d_base_name);
   if (d_read_patch_boxes || d_write_patch_boxes) {
      if (!tester_db->keyExists("test_patch_boxes_at_steps")) {
         tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
                    << "Must provide 'test_patch_boxes_at_steps' data."
                    << std::endl;
      } else {
         d_test_patch_boxes_at_steps =
            tester_db->getIntegerVector("test_patch_boxes_at_steps");
      }
   }

   if (d_test_fluxes) {

      /*
       * Read expected result for flux test...
       * Fluxes not verified in this routine.  Rather, we let it
       * write the result and do a "diff" within the script
       */

      tbox::pout << "Do a diff on the resulting *.dat file to verify result."
                 << std::endl;

   } else {

      /*
       * Read correct_result array for timestep test...
       */
      if (tester_db->keyExists("correct_result")) {
         d_correct_result = tester_db->getDoubleVector("correct_result");
      } else {
         TBOX_WARNING("main.cpp: TESTING is on but no `correct_result' array"
            << "is given in input file." << std::endl);
      }

      /* Specify whether to output "correct_result" result */

      if (tester_db->keyExists("output_correct")) {
         d_output_correct = tester_db->getBool("output_correct");
      }

   }

}

int AutoTester::checkHierarchyBoxes(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   const hier::BoxLevel& correct_box_level,
   int iter)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   const std::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));
   const hier::BoxLevel& box_level = *patch_level->getBoxLevel();

   const int local_exact_match = box_level == correct_box_level;

   int global_exact_match = local_exact_match;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_exact_match, 1, MPI_MIN);
   }

   /*
    * Check to make sure hierarchy's BoxLevel and
    * correct_box_level are identical.  If not, write an error
    * message.
    */

   int num_failures = 0;

   if (local_exact_match && global_exact_match) {
      tbox::plog << "Test 4: Level " << level_number
                 << " BoxLevel check successful for step " << iter
                 << std::endl << std::endl;
   } else {
      tbox::perr << "Test 4: FAILED: Level " << level_number
                 << " hier::BoxLevel configuration doesn't match at step " << iter
                 << std::endl << std::endl;

      hier::BoxContainer correct_minus_computed(
         correct_box_level.getGlobalizedVersion().getGlobalBoxes());
      correct_minus_computed.unorder();
      correct_minus_computed.removeIntersections(
         box_level.getGlobalizedVersion().getGlobalBoxes());

      hier::BoxContainer computed_minus_correct(
         box_level.getGlobalizedVersion().getGlobalBoxes());
      computed_minus_correct.unorder();
      computed_minus_correct.removeIntersections(
         correct_box_level.getGlobalizedVersion().getGlobalBoxes());

      tbox::plog << " global correct_box_level \\ box_level:\n"
                 << correct_minus_computed.format("\t");
      tbox::plog << " global box_level \\ correct_box_level:\n"
                 << computed_minus_correct.format("\t");
      tbox::plog << " correct_box_level:\n" << correct_box_level.format(" ", 2) << '\n'
                 << " box_level:\n" << box_level.format(" ", 2) << "\n\n";
      ++num_failures;
   }

   if (d_output_correct) {

      tbox::pout << "-------------------------------------------------------"
                 << std::endl;

      if (!local_exact_match) {
         tbox::pout << "LOCAL BOX LEVEL DOES NOT MATCH "
                    << "ON LEVEL: " << level_number << std::endl;
      }

      if (!global_exact_match) {
         tbox::pout << "GLOBAL BOX LEVEL DOES NOT MATCH "
                    << "ON LEVEL: " << level_number << std::endl;
      }
      tbox::pout << "BoxLevel: " << std::endl;
      box_level.recursivePrint(tbox::pout, "", 3);
      tbox::pout << "correct BoxLevel: " << std::endl;
      correct_box_level.recursivePrint(tbox::pout, "", 3);

      tbox::pout << "-------------------------------------------------------"
                 << std::endl << std::endl;

   }

   return num_failures;
}
