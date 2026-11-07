/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for SAMRAI Linear Advection example problem.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Header for application-specific algorithm/data structure object

#include "LinAdv.h"

// Headers for major algorithm/data structure objects

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"

// Headers for basic SAMRAI objects

#include "SAMRAI/hier/BlueprintUtils.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/ConduitDatabase.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"


// Classes for run-time plotting and autotesting.
//#define TESTING 1
#if (TESTING == 1)
#include "test/testlib/AutoTester.h"
#endif


#ifndef _MSC_VER
#include <unistd.h>
#endif

#ifdef SAMRAI_HAVE_CONDUIT
#include "conduit_blueprint.hpp"
#endif

#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <memory>

using namespace SAMRAI;

/************************************************************************
 *
 * This is the main program for an AMR solution of the linear advection
 * equation: du/dt + div(a*u) = 0, where  "u" is a scalar-valued
 * function and "a" is a constant vector.  This application program is
 * constructed by composing several algorithm objects found in the
 * SAMRAI library with a few that are specific to this application.
 * A brief description of these object follows.
 *
 * There are two main data containment objects.  These are:
 *
 *    hier::PatchHierarchy - A container for the AMR patch hierarchy and
 *       the data on the grid.
 *
 *    geom::CartesianGridGeometry - Defines and maintains the Cartesian
 *       coordinate system on the grid.  The hier::PatchHierarchy
 *       maintains a reference to this object.
 *
 * A single overarching algorithm object drives the time integration
 * and adaptive gridding processes:
 *
 *    algs::TimeRefinementIntegrator - Coordinates time integration and
 *       adaptive gridding procedures for the various levels
 *       in the AMR patch hierarchy.  Local time refinement is
 *       employed during hierarchy integration; i.e., finer
 *       levels are advanced using smaller time increments than
 *       coarser level.  Thus, this object also invokes data
 *       synchronization procedures which couple the solution on
 *       different patch hierarchy levels.
 *
 * The time refinement integrator is not specific to the numerical
 * methods used and the problem being solved.   It maintains references
 * to two other finer grain algorithmic objects that are more specific
 * to the problem at hand and with which it is configured when they are
 * passed into its constructor.   These finer grain algorithm objects
 * are:
 *
 *    algs::HyperbolicLevelIntegrator - Defines data management procedures
 *       for level integration, data synchronization between levels,
 *       and tagging cells for refinement.  These operations are
 *       tailored to explicit time integration algorithms used for
 *       hyperbolic systems of conservation laws, such as the Euler
 *       equations.  This integrator manages data for numerical
 *       routines that treat individual patches in the AMR patch
 *       hierarchy.  In this particular application, it maintains a
 *       pointer to the LinAdv object that defines variables and
 *       provides numerical routines for the linear advection problem.
 *
 *       LinAdv - Defines variables and numerical routines for the
 *          discrete linear advection equation on each patch in the
 *          AMR hierarchy.
 *
 *    mesh::GriddingAlgorithm - Drives the AMR patch hierarchy generation
 *       and regridding procedures.  This object maintains
 *       references to three other algorithmic objects with
 *       which it is configured when they are passed into its
 *       constructor.   They are:
 *
 *       mesh::BergerRigoutsos - Clusters cells tagged for refinement on a
 *          patch level into a collection of logically-rectangular
 *          box domains.
 *
 *       mesh::LoadBalancer - Processes the boxes generated by the
 *          mesh::BergerRigoutsos algorithm into a configuration from
 *          which patches are contructed.  The algorithm used in this
 *          class assumes a spatially-uniform workload distribution;
 *          thus, it attempts to produce a collection of boxes
 *          each of which contains the same number of cells.  The
 *          load balancer also assigns patches to processors.
 *
 *       mesh::StandardTagAndInitialize - Couples the gridding algorithm
 *          to the HyperbolicIntegrator. Selects cells for
 *          refinement based on either Gradient detection, Richardson
 *          extrapolation, or pre-defined Refine box region.  The
 *          object maintains a pointer to the algs::HyperbolicLevelIntegrator,
 *          which is passed into its constructor, for this purpose.
 *
 ************************************************************************
 */

/*
 *******************************************************************
 *
 * For each run, the input filename and restart information
 * (if needed) must be given on the command line.
 *
 *      For non-restarted case, command line is:
 *
 *          executable <input file name>
 *
 *      For restarted run, command line is:
 *
 *          executable <input file name> <restart directory> \
 *                     <restart number>
 *
 *******************************************************************
 */

int main(
   int argc,
   char* argv[])
{

   int num_failures = 0;
   const int number_of_runs = 2;

   /*
    * Initialize tbox::MPI.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();

   for (int run = 0; run < number_of_runs; ++run) {
      /*
       * Initialize SAMRAI, enable logging, and process command line.
       */
      tbox::SAMRAIManager::startup();
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      {
         std::string input_filename;
         std::string restart_read_dirname;
         int restore_num = 0;

         bool is_from_restart = false;

         if ((argc != 2) && (argc != 4)) {
            tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                       << "<restart dir> <restore number> [options]\n"
                       << "  options:\n"
                       << "  none at this time"
                       << std::endl;
            tbox::SAMRAI_MPI::abort();
            return -1;
         } else {
            input_filename = argv[1];
            if (argc == 4) {
               restart_read_dirname = argv[2];
               restore_num = atoi(argv[3]);

               is_from_restart = true;
            }
         }

         tbox::plog << "input_filename = " << input_filename << std::endl;
         tbox::plog << "restart_read_dirname = " << restart_read_dirname << std::endl;
         tbox::plog << "restore_num = " << restore_num << std::endl;

         /*
          * Create input database and parse all data in input file.
          */

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
#ifdef SGS
            if (global_db->keyExists("tag_clustering_method")) {
               std::string tag_clustering_method =
                  global_db->getString("tag_clustering_method");
               mesh::BergerRigoutsos::setClusteringOption(tag_clustering_method);
            }
#endif
            if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
               bool flag = global_db->
                  getBool("call_abort_in_serial_instead_of_exit");
               tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
            }
         }

         /*
          * Retrieve "Main" section of the input database.  First, read
          * dump information, which is used for writing plot files.
          * Second, if proper restart information was given on command
          * line, and the restart interval is non-zero, create a restart
          * database.
          */

         std::shared_ptr<tbox::Database> main_db(
            input_db->getDatabase("Main"));

         const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

         const std::string base_name =
            main_db->getStringWithDefault("base_name", "unnamed");

         const std::string log_filename =
            main_db->getStringWithDefault("log_filename", base_name + ".log");

         bool log_all_nodes = false;
         if (main_db->keyExists("log_all_nodes")) {
            log_all_nodes = main_db->getBool("log_all_nodes");
         }
         if (log_all_nodes) {
            tbox::PIO::logAllNodes(log_filename);
         } else {
            tbox::PIO::logOnlyNodeZero(log_filename);
         }

         int viz_dump_interval = 0;
         if (main_db->keyExists("viz_dump_interval")) {
            viz_dump_interval = main_db->getInteger("viz_dump_interval");
         }

         const std::string viz_dump_dirname =
            main_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");
         int visit_number_procs_per_file = 1;

         const bool viz_dump_data = (viz_dump_interval > 0);

         bool write_blueprint = false;
#ifdef SAMRAI_HAVE_CONDUIT
         write_blueprint =
            main_db->getBoolWithDefault("write_blueprint", false);
#endif

         int restart_interval = 0;
         if (main_db->keyExists("restart_interval")) {
            restart_interval = main_db->getInteger("restart_interval");
         }

         const std::string restart_write_dirname =
            main_db->getStringWithDefault("restart_write_dirname",
               base_name + ".restart");

         bool use_refined_timestepping = true;
         if (main_db->keyExists("timestepping")) {
            std::string timestepping_method = main_db->getString("timestepping");
            if (timestepping_method == "SYNCHRONIZED") {
               use_refined_timestepping = false;
            }
         }

         bool rebalance_coarsest = false;
         if (main_db->keyExists("rebalance_coarsest")) {
            rebalance_coarsest = main_db->getBool("rebalance_coarsest");
         }

#if (TESTING == 1) && !defined(HAVE_HDF5)
         /*
          * If we are autotesting on a system w/o HDF5, the read from
          * restart will result in an error.  We want this to happen
          * for users, so they know there is a problem with the restart,
          * but we don't want it to happen when autotesting.
          */
         is_from_restart = false;
         restart_interval = 0;
#endif

         const bool write_restart = (restart_interval > 0)
            && !(restart_write_dirname.empty());

         /*
          * Get the restart manager and root restart database.  If run is from
          * restart, open the restart file.
          */

         tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

         if (is_from_restart) {
            restart_manager->
            openRestartFile(restart_read_dirname, restore_num,
               mpi.getSize());
         }

         /*
          * Create major algorithm and data objects which comprise application.
          * Each object will be initialized either from input data or restart
          * files, or a combination of both.  Refer to each class constructor
          * for details.  For more information on the composition of objects
          * for this application, see comments at top of file.
          */

         std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
            new geom::CartesianGridGeometry(
               dim,
               "CartesianGeometry",
               input_db->getDatabase("CartesianGeometry")));

         std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
            new hier::PatchHierarchy(
               "PatchHierarchy",
               grid_geometry,
               input_db->getDatabase("PatchHierarchy")));

         LinAdv* linear_advection_model = new LinAdv(
               "LinAdv",
               dim,
               input_db->getDatabase("LinAdv"),
               grid_geometry);

         std::shared_ptr<algs::HyperbolicLevelIntegrator> hyp_level_integrator(
            new algs::HyperbolicLevelIntegrator(
               "HyperbolicLevelIntegrator",
               input_db->getDatabase("HyperbolicLevelIntegrator"),
               linear_advection_model,
               use_refined_timestepping));

         std::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
            new mesh::StandardTagAndInitialize(
               "StandardTagAndInitialize",
               hyp_level_integrator.get(),
               input_db->getDatabase("StandardTagAndInitialize")));

         std::shared_ptr<mesh::BergerRigoutsos> box_generator(
            new mesh::BergerRigoutsos(
               dim,
               input_db->getDatabaseWithDefault(
                  "BergerRigoutsos",
                  std::shared_ptr<tbox::Database>())));
         box_generator->useDuplicateMPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

         std::shared_ptr<mesh::CascadePartitioner> load_balancer(
            new mesh::CascadePartitioner(
               dim,
               "LoadBalancer",
               input_db->getDatabase("LoadBalancer")));
         load_balancer->setSAMRAI_MPI(
            tbox::SAMRAI_MPI::getSAMRAIWorld());

         std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
            new mesh::GriddingAlgorithm(
               patch_hierarchy,
               "GriddingAlgorithm",
               input_db->getDatabase("GriddingAlgorithm"),
               error_detector,
               box_generator,
               load_balancer));

         std::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
            new algs::TimeRefinementIntegrator(
               "TimeRefinementIntegrator",
               input_db->getDatabase("TimeRefinementIntegrator"),
               patch_hierarchy,
               hyp_level_integrator,
               gridding_algorithm));

         // VisItDataWriter is only present if HDF is available
#ifdef HAVE_HDF5
         std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
            new appu::VisItDataWriter(
               dim,
               "LinAdv VisIt Writer",
               viz_dump_dirname,
               visit_number_procs_per_file));
         linear_advection_model->
         registerVisItDataWriter(visit_data_writer);
#endif

         /*
          * Initialize hierarchy configuration and data on all patches.
          * Then, close restart file and write initial state for visualization.
          */

         double dt_now = time_integrator->initializeHierarchy();

         tbox::RestartManager::getManager()->closeRestartFile();

#if (TESTING == 1)
         /*
          * Create the autotesting object which will verify correctness
          * of the problem. If no automated testing is done, the object does
          * not get used.
          */
         AutoTester autotester("AutoTester", dim, input_db);
         std::cerr << "autotesting object initialized!!" << std::endl;
#endif

         /*
          * After creating all objects and initializing their state, we
          * print the input database and variable database contents
          * to the log file.
          */

         tbox::plog << "\nCheck input data and variables before simulation:"
                    << std::endl;
         tbox::plog << "Input database..." << std::endl;
         input_db->printClassData(tbox::plog);
         tbox::plog << "\nVariable database..." << std::endl;
         hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

         tbox::plog << "\nCheck Linear Advection data... " << std::endl;
         linear_advection_model->printClassData(tbox::plog);
         if (viz_dump_data &&
             time_integrator->getIntegratorStep() % viz_dump_interval == 0) {
#ifdef HAVE_HDF5
            visit_data_writer->writePlotData(
               patch_hierarchy,
               time_integrator->getIntegratorStep(),
               time_integrator->getIntegratorTime());
#endif
         }

         /*
          * Time step loop.  Note that the step count and integration
          * time are maintained by algs::TimeRefinementIntegrator.
          */

         double loop_time = time_integrator->getIntegratorTime();
         double loop_time_end = time_integrator->getEndTime();

         int iteration_num = time_integrator->getIntegratorStep();

#if (TESTING == 1)
         /*
          * If we are doing autotests, check result...
          */
         num_failures += autotester.evalTestData(iteration_num,
               patch_hierarchy,
               time_integrator,
               hyp_level_integrator,
               gridding_algorithm);
#endif

         while ((loop_time < loop_time_end) &&
                time_integrator->stepsRemaining()) {

            iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            tbox::pout << "At begining of timestep # " << iteration_num - 1
                       << std::endl;
            tbox::pout << "Simulation time is " << loop_time << std::endl;

            double dt_new = time_integrator->advanceHierarchy(
               dt_now,
               rebalance_coarsest);

            loop_time += dt_now;
            dt_now = dt_new;

            tbox::pout << "At end of timestep # " << iteration_num - 1 << std::endl;
            tbox::pout << "Simulation time is " << loop_time << std::endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

            /*
             * At specified intervals, write restart and visualization files.
             */
            if (write_restart) {

               if ((iteration_num % restart_interval) == 0) {
                  tbox::RestartManager::getManager()->
                  writeRestartFile(restart_write_dirname,
                     iteration_num);
               }
            }

            /*
             * At specified intervals, write out data files for plotting.
             */
            if (viz_dump_data) {
               if ((iteration_num % viz_dump_interval) == 0) {
#ifdef HAVE_HDF5
                  visit_data_writer->writePlotData(patch_hierarchy,
                     iteration_num,
                     loop_time);
#endif
               }
            }
            if (write_blueprint && viz_dump_interval) {
               if ((iteration_num % viz_dump_interval) == 0) {
#ifdef SAMRAI_HAVE_CONDUIT
                  std::shared_ptr<tbox::ConduitDatabase> conduit_db(
                     new tbox::ConduitDatabase("conduit_hierarchy"));
 
                  hier::BlueprintUtils bp_utils(linear_advection_model);
                  patch_hierarchy->makeBlueprintDatabase(conduit_db, bp_utils);

                  conduit::Node bp_node;
                  conduit_db->toConduitNode(bp_node);

                  std::vector<int> first_patch_id;
                  first_patch_id.push_back(0);

                  int num_levels = patch_hierarchy->getNumberOfLevels();
                  int patch_count = 0;
                  for (int i = 1; i < num_levels; ++i) {
                     patch_count +=
                        patch_hierarchy->getPatchLevel(i-1)->
                           getNumberOfPatches();
                     first_patch_id.push_back(patch_count);
                  }

                  int num_hier_patches = 0;
                  for (int i = 0; i < num_levels; ++i) {
                     const std::shared_ptr<hier::PatchLevel>& level =
                        patch_hierarchy->getPatchLevel(i);
                     num_hier_patches +=
                        level->getNumberOfPatches();

                     for (hier::PatchLevel::Iterator p(level->begin());
                          p != level->end(); ++p) {

                        const std::shared_ptr<hier::Patch>& patch = *p;
                        const hier::BoxId& box_id = patch->getBox().getBoxId();
                        const hier::LocalId& local_id = box_id.getLocalId();

                        int mesh_id = first_patch_id[i] + local_id.getValue();

                        linear_advection_model->addFields(
                           bp_node, mesh_id, patch); 
                     }
                  }

                  conduit::Node verify_info;
                  bool bp_verified =
                     conduit::blueprint::verify("mesh", bp_node, verify_info);

                  if (bp_verified) {
                     bp_utils.writeBlueprintMesh(
                        bp_node,
                        tbox::SAMRAI_MPI::getSAMRAIWorld(),
                        num_hier_patches,
                        "amr_mesh",
                        "LinAdvData",
                        "LinAdv.root",
                        "json");
                  } else {
                     num_failures += 1;
                  } 
#endif
               }
            }


#if (TESTING == 1)
            /*
             * If we are doing autotests, check result...
             */
            num_failures += autotester.evalTestData(iteration_num,
                  patch_hierarchy,
                  time_integrator,
                  hyp_level_integrator,
                  gridding_algorithm);
#endif

         }

         /*
          * At conclusion of simulation, deallocate objects.
          */
#ifdef HAVE_HDF5
         visit_data_writer.reset();
#endif

         time_integrator.reset();
         gridding_algorithm.reset();
         load_balancer.reset();
         box_generator.reset();
         error_detector.reset();
         hyp_level_integrator.reset();

         if (linear_advection_model) delete linear_advection_model;

         patch_hierarchy.reset();
         grid_geometry.reset();

         input_db.reset();
         main_db.reset();

      }
      tbox::SAMRAIManager::shutdown();
   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  LinAdv" << std::endl;
   }

   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}
