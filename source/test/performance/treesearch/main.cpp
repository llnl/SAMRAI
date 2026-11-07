/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Performance tests for tree searches.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxTree.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>
#include <random>
#include <vector>
#include <iomanip>

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 *
 * This is a performance test for the tree search algorithm
 * in BoxTree:
 *
 * 1. Generate a set of Boxes.
 *
 * 2. Sort the Boxes into trees using layerNodeTree.
 *
 * 3. Search for overlaps.
 *
 *************************************************************************
 */

typedef std::vector<hier::Box> BoxVec;

/*
 * Generate uniform boxes as specified in the database.
 */
void
generateBoxesUniform(
   const tbox::Dimension& dim,
   std::vector<hier::Box>& output,
   const std::shared_ptr<Database>& db);

int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */

   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   int fail_count = 0;

   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       * executable <input file name>
       */
      std::string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << std::endl);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(
         new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Set up the timer manager.
       */
      if (input_db->isDatabase("TimerManager")) {
         TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_file_name);
      } else {
         PIO::logOnlyNodeZero(log_file_name);
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      tbox::TimerManager * tm(tbox::TimerManager::getManager());
      const std::string dim_str(tbox::Utilities::intToString(dim.getValue()));
      std::shared_ptr<tbox::Timer> t_build_tree(
         tm->getTimer("apps::main::build_tree[" + dim_str + "]"));
      std::shared_ptr<tbox::Timer> t_search_tree_for_set(
         tm->getTimer("apps::main::search_tree_for_set[" + dim_str + "]"));
      std::shared_ptr<tbox::Timer> t_search_tree_for_vec(
         tm->getTimer("apps::main::search_tree_for_vec[" + dim_str + "]"));

      /*
       * Generate the boxes.
       */
      BoxVec boxes;
      generateBoxesUniform(dim,
         boxes,
         main_db->getDatabase("UniformBoxGen"));
      tbox::plog << "\n\n\nGenerated boxes (" << boxes.size() << "):\n";
      for (size_t i = 0; i < boxes.size(); ++i) {
         tbox::plog << '\t' << i << '\t' << boxes[i] << '\n';
         if (i > 20) {
            tbox::plog << "\t...\n";
            break;
         }
      }
      tbox::plog << "\n\n\n";

      /*
       * Compute bounding box.
       */
      hier::Box bounding_box(dim);
      for (BoxVec::iterator bi = boxes.begin(); bi != boxes.end(); ++bi) {
         bounding_box += *bi;
      }

      /*
       * Randomize the boxes.
       */
      bool randomize_order = main_db->getBoolWithDefault("randomize_order",
            false);
      if (randomize_order) {
         std::shuffle(boxes.begin(), boxes.end(),
            std::mt19937{std::random_device{}()});
      }

      /*
       * Scale up the number of boxes and time the sort and search for the
       * growing set of boxes.
       */
      size_t num_scale = (size_t)main_db->getIntegerWithDefault("num_scale", 1);
      for (unsigned int iscale = 0; iscale < num_scale; ++iscale) {

         if (iscale != 0) {
            /*
             * Scale up the box array.
             */
            tbox::Dimension::dir_t shift_dir =
               static_cast<tbox::Dimension::dir_t>((iscale - 1) % dim.getValue());
            /*
             * Shift distance is less than number of bounding boxes in shift_dir
             * in order to generate some non-trivial overlaps.
             */
            int shift_distance = int(0.91 * bounding_box.numberCells(shift_dir));

            const size_t old_size = boxes.size();
            boxes.insert(boxes.end(), boxes.begin(), boxes.end());
            for (size_t i = 0; i < old_size; ++i) {
               boxes[i].shift(shift_dir, shift_distance);
            }
            bounding_box.setUpper(shift_dir, shift_distance);
         }

         if (mpi.getRank() == 0) {
            tbox::pout << "Repetition " << iscale << std::endl;
         }
         tbox::plog << "Repetition " << iscale << " has "
                    << boxes.size() << " boxes bounded by "
                    << bounding_box << std::endl;

         /*
          * Generate the nodes from the boxes.
          */
         hier::BoxContainer nodes;
         for (hier::LocalId i(0); i < static_cast<int>(boxes.size()); ++i) {
            nodes.insert(nodes.end(),
               hier::Box(boxes[i.getValue()], i, 0));
         }
         const size_t node_count = nodes.size();

         /*
          * Grow the boxes for overlap search.
          */
         hier::IntVector growth(dim, 1);
         if (main_db->isInteger("growth")) {
            main_db->getIntegerArray("growth", &growth[0], dim.getValue());
         }
         BoxVec grown_boxes = boxes;
         if (growth != 1) {
            for (BoxVec::iterator bi = grown_boxes.begin();
                 bi != grown_boxes.end();
                 ++bi) {
               bi->grow(growth);
            }
         }

         /*
          * Reset timers and statistics.
          */
         tm->resetAllTimers();
         hier::BoxTree::resetStatistics(dim);

         /*
          * Build search tree.
          */
         t_build_tree->start();
         nodes.makeTree(0);
         t_build_tree->stop();

         /*
          * Search the tree.
          *
          * We test outputing in an unordered and ordered container.  The
          * can indicate the difference in performance due to sorting the
          * output for an ordered container.
          */
         hier::BoxContainer unordered_overlap;
         t_search_tree_for_set->start();
         for (BoxVec::iterator bi = grown_boxes.begin();
              bi != grown_boxes.end();
              ++bi) {
            unordered_overlap.clear();
            nodes.findOverlapBoxes(unordered_overlap, *bi);
         }
         t_search_tree_for_set->stop();

         hier::BoxContainer ordered_overlap;
         ordered_overlap.order();
         t_search_tree_for_vec->start();
         for (BoxVec::iterator bi = grown_boxes.begin();
              bi != grown_boxes.end();
              ++bi) {
            ordered_overlap.clear();
            nodes.findOverlapBoxes(ordered_overlap, *bi);
         }
         t_search_tree_for_vec->stop();

         /*
          * Output normalized timer to plog.
          */
         tbox::plog << "Timers for repetition " << iscale
                    << " (normalized by " << node_count << " nodes):\n";
         tbox::plog.precision(8);
         tbox::plog << t_build_tree->getName() << " = "
                    << t_build_tree->getTotalWallclockTime()
         / static_cast<double>(node_count)
                    << std::endl;
         tbox::plog << t_search_tree_for_set->getName() << " = "
                    << t_search_tree_for_set->getTotalWallclockTime()
         / static_cast<double>(node_count)
                    << std::endl;
         tbox::plog << t_search_tree_for_vec->getName() << " = "
                    << t_search_tree_for_vec->getTotalWallclockTime()
         / static_cast<double>(node_count)
                    << std::endl;

         /*
          * Log timer results and search tree statistics.
          */
         tbox::TimerManager::getManager()->print(tbox::plog);
         hier::BoxTree::printStatistics(dim);

         tbox::plog << "\n\n\n";

      }

      /*
       * Print input database again to fully show usage.
       */
      plog << "Input database after running..." << std::endl;
      input_db->printClassData(plog);

      tbox::pout << "\nPASSED:  Tree search" << std::endl;

      input_db.reset();
      main_db.reset();
      t_search_tree_for_set.reset();
      t_search_tree_for_vec.reset();

      /*
       * Exit properly by shutting down services in correct order.
       */
      tbox::plog << "\nShutting down..." << std::endl;

   }

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   return fail_count;
}

/*
 * Function to generate a uniform set of boxes.
 */
void generateBoxesUniform(
   const tbox::Dimension& dim,
   std::vector<hier::Box>& output,
   const std::shared_ptr<Database>& db)
{
   output.clear();

   hier::IntVector boxsize(dim, 1);
   if (db->isInteger("boxsize")) {
      db->getIntegerArray("boxsize", &boxsize[0], dim.getValue());
   } else {
      TBOX_ERROR("generateBoxesUniform() error...\n"
         << "    box size is absent.");
   }

   hier::IntVector boxrepeat(dim, 1);
   if (db->isInteger("boxrepeat")) {
      db->getIntegerArray("boxrepeat", &boxrepeat[0], dim.getValue());
   }

   /*
    * Create an array of boxes by repeating the given box.
    */
   hier::Index index(dim, 0);
   do {
      hier::Index lower(index * boxsize);
      hier::Index upper(lower + boxsize - 1);
      int& e = index(0);
      for (e = 0; e < boxrepeat(0); ++e) {
         lower(0) = e * boxsize(0);
         upper(0) = lower(0) + boxsize(0) - 1;
         output.insert(output.end(), hier::Box(lower, upper, hier::BlockId(0)));
      }
      for (int d = 0; d < dim.getValue(); ++d) {
         if (index(d) == boxrepeat(d) && d < dim.getValue() - 1) {
            index(d) = 0;
            ++index(d + 1);
         }
      }
   } while (index(dim.getValue() - 1) < boxrepeat(dim.getValue() - 1));
}
