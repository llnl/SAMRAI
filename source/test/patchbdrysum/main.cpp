/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for test of hierarchy sum
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Headers for basic SAMRAI objects
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"

// Headers for major algorithm/data structure objects
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"

// Header for application-specific algorithm/data structure object
#include "HierSumTest.h"

#include <memory>

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace geom;
using namespace mesh;

/**
 * This is the main program for an example case that tests SAMRAI
 * hierarchy sum classes for FE-type operations with node data,
 * and level edge sum operations with edge data.
 *
 * The main program constructs the various SAMRAI gridding objects
 * and performs the time-stepping loop.  The program should be
 * executed as:
 *
 *    executable <input file name>
 */

int main(
   int argc,
   char* argv[])
{

   int fail_count = 0;

   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

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
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                    << "[options]\n" << std::endl;
         tbox::SAMRAI_MPI::abort();
         return -1;
      } else {
         input_filename = argv[1];
      }

      tbox::plog << "input_filename = " << input_filename << std::endl;

      /****************************************************************
      *
      *  PROBLEM SETUP
      *
      ****************************************************************
      *
      *  Read data from input file and initialize SAMRAI classes
      *
      ****************************************************************/

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section of the input database.
       */

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      /*
       * Determine if we are doing node sum tests, edge sum tests,
       * or both.
       */
      bool do_node_sum = false;
      if (main_db->keyExists("do_node_sum")) {
         do_node_sum = main_db->getBool("do_node_sum");
      }
      bool do_edge_sum = false;
      if (main_db->keyExists("do_edge_sum")) {
         do_edge_sum = main_db->getBool("do_edge_sum");
      }

      int nsteps = 1;
      if (main_db->keyExists("nsteps")) {
         nsteps = main_db->getInteger("nsteps");
      }

      std::string log_file_name = "hiersumtest.log";
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

      std::string visit_dump_dirname = "visit_data";
      int visit_number_procs_per_file = 1;

      int visit_dump_interval =
         main_db->getIntegerWithDefault("visit_dump_interval", 0);

      if (visit_dump_interval > 0) {
         if (main_db->keyExists("visit_dump_dirname")) {
            visit_dump_dirname = main_db->getString("visit_dump_dirname");
         }
         if (main_db->keyExists("visit_number_procs_per_file")) {
            visit_number_procs_per_file =
               main_db->getInteger("visit_number_procs_per_file");
         }
      }

      /*
       * The grid geometry defines the grid type (e.g. cartesian, spherical,
       * etc.).  Because SAMRAI operates on block structured indices, it can
       * support any grid geometry that may be represented as an orthogonal
       * grid.
       */
      std::shared_ptr<CartesianGridGeometry> grid_geometry(
         new CartesianGridGeometry(dim,
            "CartesianGeometry",
            input_db->getDatabase("CartesianGeometry")));

      /*
       * The patch hierarchy defines the adaptive grid system.
       */
      std::shared_ptr<PatchHierarchy> patch_hierarchy(
         new PatchHierarchy(
            "PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

#ifdef HAVE_HDF5
      /*
       * Set up Visualization writer.
       */
      std::shared_ptr<appu::VisItDataWriter> visit_data_writer;
      if (visit_dump_interval > 0) {
         visit_data_writer.reset(
            new appu::VisItDataWriter(
               dim,
               "HierSumTest VisIt Writer",
               visit_dump_dirname,
               visit_number_procs_per_file));
      }
#endif

      /*
       * This is our problem class.  See the class header for comments on it.
       */
      HierSumTest* hier_sum_test = new HierSumTest(
            "HierSumTest",
            dim,
            input_db->getDatabase("HierSumTest")
#ifdef HAVE_HDF5
            , visit_data_writer
#endif
            );

      /*
       * The StandardTagAndInitialize class performs a variety of operations
       * with user-specified parameters related to adptive gridding.  For example,
       * it manages initialization of a level, cell tagging using a gradient
       * detector, and methods to reset data after the hierarchy has been
       * regridded.
       */
      std::shared_ptr<StandardTagAndInitialize> tag_and_init_ops(
         new StandardTagAndInitialize(
            "StandardTagAndInitialize",
            hier_sum_test,
            input_db->getDatabase("StandardTagAndInitialize")));

      /*
       * The gridding algorithm manages adaptive gridding.  It expects a
       * clustering scheme (i.e. how to cluster tagged-cells into patches),
       * and a load balance scheme to distribute work to processors.  In general
       * the baseline classes provided in SAMRAI should suffice for most
       * problems. It also requires a class that defines the particular tag
       * and initialization ops that correlate with the users problem.  For
       * this, we use the "tag_and_init_ops" above, which references our
       * "wave_eqn_model" problem class to define the user-specific operations.
       */
      std::shared_ptr<BergerRigoutsos> box_generator(
         new BergerRigoutsos(dim,
            input_db->getDatabase("BergerRigoutsos")));

      std::shared_ptr<TreeLoadBalancer> load_balancer(
         new TreeLoadBalancer(dim,
            "LoadBalancer",
            input_db->getDatabase("LoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      std::shared_ptr<GriddingAlgorithm> gridding_algorithm(
         new GriddingAlgorithm(
            patch_hierarchy,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            tag_and_init_ops,
            box_generator,
            load_balancer));

      /*
       * After creating all objects and initializing their state, we
       * print the input database and variable database contents to
       * the log file.
       */

      tbox::plog << "\nCheck input data and variables before simulation:"
                 << std::endl;
      tbox::plog << "Input database..." << std::endl;
      input_db->printClassData(plog);
      tbox::plog << "\nVariable database..." << std::endl;
      VariableDatabase::getDatabase()->printClassData(plog);

      /****************************************************************
      *
      *  INITIALIZE DATA ON PATCHES
      *
      ****************************************************************
      *
      *  Build patch hierarchy and initialize the data on the patches
      *  in the hierarchy.
      *  1) Create a "tag_buffer" for each level in the Hierarchy.
      *  2) Create the coarse (i.e. level 0) grid.
      *  3) Cycle through levels 1-max_levels, initializing data
      *     on each.  The makeFinerLevel method calls the error
      *     estimator (remember, it was registered with the
      *     gridding algorithm object) and tags cells for refinement
      *     as it generates patches on the finer levels.
      *
      ****************************************************************/

      double loop_time = 0.;
      int loop_cycle = 0;
      std::vector<int> tag_buffer_array(patch_hierarchy->getMaxNumberOfLevels());
      for (int il = 0; il < patch_hierarchy->getMaxNumberOfLevels(); ++il) {
         tag_buffer_array[il] = 1;
      }
      gridding_algorithm->makeCoarsestLevel(loop_time);

      bool done = false;
      bool initial_cycle = true;
      for (int ln = 0;
           patch_hierarchy->levelCanBeRefined(ln) && !done;
           ++ln) {
         gridding_algorithm->makeFinerLevel(
            tag_buffer_array[ln],
            initial_cycle,
            loop_cycle,
            loop_time);
         done = !(patch_hierarchy->finerLevelExists(ln));
      }

      tbox::plog
      << "************************************************************\n";
      tbox::plog
      << "************************* Hierarchy ************************\n";
      tbox::plog
      << "************************************************************\n";
      patch_hierarchy->recursivePrint(tbox::plog, "", 3);
      tbox::plog << "\n\n";

      int nlevels = patch_hierarchy->getNumberOfLevels();

      for (int pln = 0; pln <= patch_hierarchy->getFinestLevelNumber();
           ++pln) {
         std::shared_ptr<PatchLevel> level(
            patch_hierarchy->getPatchLevel(pln));

         tbox::plog << "\n PRINTING PATCHES ON LEVEL " << pln << std::endl;

         for (PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            tbox::plog << "patch # " << ip->getBox().getBoxId() << " : "
                       << ip->getBox() << std::endl;
         }
      }
      tbox::plog << std::endl;

      /*******************************************************************
       *
       * Test hier sum operation
       *
       *******************************************************************
       *
       *  1) Set node values (initial)
       *  2) Do hierarchy sum operation
       *  3) Check result
       *
       ******************************************************************/

      /*
       * Setup the sum operation(s)
       */
      if (do_node_sum) {
         hier_sum_test->setupOuternodeSum(patch_hierarchy);
      }
      if (do_edge_sum) {
         for (int ln = 0; ln < nlevels; ++ln) {
            hier_sum_test->setupOuteredgeSum(patch_hierarchy,
               ln);
         }
      }

      for (int i = 0; i < nsteps; ++i) {

         /*
          * In the process of constructing the hierarchy, we set cell values
          * to their proper weights.  Now go in and set the node/edge values.
          * (write data to VisIt once it is set).
          */
         if (do_node_sum) {
            fail_count += hier_sum_test->setInitialNodeValues(patch_hierarchy);
         }
         if (do_edge_sum) {
            for (int ln = 0; ln < nlevels; ++ln) {
               std::shared_ptr<PatchLevel> level(
                  patch_hierarchy->getPatchLevel(ln));
               fail_count += hier_sum_test->setInitialEdgeValues(level);
            }
         }

#ifdef HAVE_HDF5
         /*
          * Write the pre-summed cell/node data to VisIt
          */
         if ((visit_dump_interval > 0) &&
             ((i % visit_dump_interval) == 0)) {
            visit_data_writer->writePlotData(patch_hierarchy, i, loop_time);
         }
#endif

         /*
          * Perform the sum operation(s)
          */
         if (do_node_sum) {
            hier_sum_test->doOuternodeSum();
         }
         if (do_edge_sum) {
            for (int ln = 0; ln < nlevels; ++ln) {
               hier_sum_test->doOuteredgeSum(ln);
            }
         }
      }

      /*
       * Check result
       */
      if (do_node_sum) {
         fail_count += hier_sum_test->checkNodeResult(patch_hierarchy);
      }

      tbox::pout << "\n" << std::endl;

      if (do_edge_sum) {
         for (int ln = 0; ln < nlevels; ++ln) {
            std::shared_ptr<PatchLevel> level(
               patch_hierarchy->getPatchLevel(ln));
            fail_count += hier_sum_test->checkEdgeResult(level);
         }
      }

#ifdef HAVE_HDF5
      /*
       * Write the post-summed cell/node data to VisIt
       */
      if (visit_dump_interval > 0) {
         visit_data_writer->writePlotData(patch_hierarchy,
            nsteps + 1,
            loop_time);
      }
#endif

      /*
       * At conclusion of simulation, deallocate objects.
       */
#ifdef HAVE_HDF5
      visit_data_writer.reset();
#endif
      gridding_algorithm.reset();
      load_balancer.reset();
      box_generator.reset();
      tag_and_init_ops.reset();

      if (hier_sum_test) delete hier_sum_test;

      patch_hierarchy.reset();
      grid_geometry.reset();

      input_db.reset();
      main_db.reset();

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  patchbdrysum" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
