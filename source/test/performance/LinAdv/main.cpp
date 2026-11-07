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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <memory>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// Headers for basic SAMRAI objects

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

// Headers for major algorithm/data structure objects

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TileClustering.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#define RECORD_STATS
//#undef RECORD_STATS

// Header for application-specific algorithm/data structure object

#include "LinAdv.h"
#include "test/testlib/SinusoidalFrontGenerator.h"
#include "test/testlib/SphericalShellGenerator.h"
#include "test/testlib/MeshGenerationStrategy.h"


using namespace SAMRAI;

/*
 ************************************************************************
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
 * hier::PatchHierarchy - A container for the AMR patch hierarchy and
 *    the data on the grid.
 *
 * geom::CartesianGridGeometry - Defines and maintains the Cartesian
 *    coordinate system on the grid.  The hier::PatchHierarchy
 *    maintains a reference to this object.
 *
 * A single overarching algorithm object drives the time integration
 * and adaptive gridding processes:
 *
 * algs::TimeRefinementIntegrator - Coordinates time integration and
 *    adaptive gridding procedures for the various levels
 *    in the AMR patch hierarchy.  Local time refinement is
 *    employed during hierarchy integration; i.e., finer
 *    levels are advanced using smaller time increments than
 *    coarser level.  Thus, this object also invokes data
 *    synchronization procedures which couple the solution on
 *    different patch hierarchy levels.
 *
 * The time refinement integrator is not specific to the numerical
 * methods used and the problem being solved.   It maintains references
 * to two other finer grain algorithmic objects that are more specific
 * to the problem at hand and with which it is configured when they are
 * passed into its constructor.   These finer grain algorithm objects
 * are:
 *
 * algs::HyperbolicLevelIntegrator - Defines data management procedures
 *    for level integration, data synchronization between levels,
 *    and tagging cells for refinement.  These operations are
 *    tailored to explicit time integration algorithms used for
 *    hyperbolic systems of conservation laws, such as the Euler
 *    equations.  This integrator manages data for numerical
 *    routines that treat individual patches in the AMR patch
 *    hierarchy.  In this particular application, it maintains a
 *    pointer to the LinAdv object that defines variables and
 *    provides numerical routines for the linear advection problem.
 *
 *    LinAdv - Defines variables and numerical routines for the
 *       discrete linear advection equation on each patch in the
 *       AMR hierarchy.
 *
 * mesh::GriddingAlgorithm - Drives the AMR patch hierarchy generation
 *    and regridding procedures.  This object maintains
 *    references to three other algorithmic objects with
 *    which it is configured when they are passed into its
 *    constructor.   They are:
 *
 *    mesh::BergerRigoutsos - Clusters cells tagged for refinement on a
 *       patch level into a collection of logically-rectangular
 *       box domains.
 *
 *    mesh::LoadBalancer - Processes the boxes generated by the
 *       mesh::BergerRigoutsos algorithm into a configuration from
 *       which patches are contructed.  The algorithm used in this
 *       class assumes a spatially-uniform workload distribution;
 *       thus, it attempts to produce a collection of boxes
 *       each of which contains the same number of cells.  The
 *       load balancer also assigns patches to processors.
 *
 *    mesh::StandardTagAndInitialize - Couples the gridding algorithm
 *       to the HyperbolicIntegrator. Selects cells for
 *       refinement based on either Gradient detection, Richardson
 *       extrapolation, or pre-defined Refine box region.  The
 *       object maintains a pointer to the algs::HyperbolicLevelIntegrator,
 *       which is passed into its constructor, for this purpose.
 *
 ************************************************************************
 */

/*
 *******************************************************************
 *
 * For each run, the input filename and restart information
 * (if needed) must be given on the command line.
 *
 *   For non-restarted case, command line is:
 *
 *       executable <input file name>
 *
 *   For restarted run, command line is:
 *
 *       executable <input file name> <restart directory> \
 *                  <restart number>
 *
 *******************************************************************
 */

int main(
   int argc,
   char* argv[])
{

   using namespace tbox;

   /*
    * Initialize MPI and SAMRAI, enable logging, and process command line.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   int num_failures = 0;

   {
      std::string input_filename;
      std::string case_name;
      int scale_size = mpi.getSize();

      if ((argc != 2) && (argc != 3) && (argc != 4)) {
         tbox::pout << "USAGE:\n"
                    << argv[0] << " <input filename> "
                    << "or\n"
                    << argv[0] << " <input filename> <case name>"
                    << "or\n"
                    << argv[0] << " <input filename> <case name> <scale size> "
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
         return -1;
      } else {
         input_filename = argv[1];
         if (argc > 2) {
            case_name = argv[2];
         }
         if (argc > 3) {
            scale_size = atoi(argv[3]);
         }
      }

      pout << "input_filename = " << input_filename << std::endl;
      pout << "case_name = " << case_name << std::endl;
      pout << "scale_size = " << scale_size << std::endl;

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      std::shared_ptr<tbox::Timer> t_all =
         tbox::TimerManager::getManager()->getTimer("appu::main::all");
      t_all->start();

      std::shared_ptr<tbox::Timer> t_vis_writing(
         tbox::TimerManager::getManager()->getTimer("apps::Main::vis_writing"));

      /*
       * Retrieve "Main" section of the input database.  First, read
       * dump information, which is used for writing plot files.
       * Second, if proper restart information was given on command
       * line, and the restart interval is non-zero, create a restart
       * database.
       */

      std::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      bool use_scaled_input = main_db->getBoolWithDefault("use_scaled_input",
            true);

      std::string scaled_input_str =
         std::string("ScaledInput")
         + (use_scaled_input ? tbox::Utilities::intToString(scale_size) : std::string());
      std::shared_ptr<Database> scaled_input_db(input_db->getDatabase(scaled_input_str));

      std::string base_name = main_db->getStringWithDefault("base_name", "unnamed");

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      std::string base_name_ext = base_name;
      if (!case_name.empty()) {
         base_name_ext = base_name_ext + '-' + case_name;
      }
      base_name_ext = base_name_ext + '-' + tbox::Utilities::nodeToString(scale_size);
      tbox::pout << "Added case name (" << case_name << ") and nprocs ("
                 << mpi.getSize() << ") to base name -> '"
                 << base_name_ext << "'\n";

      /*
       * Logging.
       */
      std::string log_filename = base_name_ext + ".log";
      log_filename =
         main_db->getStringWithDefault("log_filename", base_name_ext + ".log");

      bool log_all_nodes = false;
      log_all_nodes =
         main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_filename);
      } else {
         PIO::logOnlyNodeZero(log_filename);
      }

      int viz_dump_interval = 0;
      if (main_db->keyExists("viz_dump_interval")) {
         viz_dump_interval = main_db->getInteger("viz_dump_interval");
      }

      const std::string viz_dump_dirname = main_db->getStringWithDefault(
            "viz_dump_dirname", base_name_ext + ".visit");
      int visit_number_procs_per_file = 1;

      const bool viz_dump_data = (viz_dump_interval > 0);

      int restart_interval = 0;
      if (main_db->keyExists("restart_interval")) {
         restart_interval = main_db->getInteger("restart_interval");
      }

      std::string restart_write_dirname;
      if (restart_interval > 0) {
         if (main_db->keyExists("restart_write_dirname")) {
            restart_write_dirname = main_db->getString("restart_write_dirname");
         } else {
            TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname'"
               << " not specifed in input file");
         }
      }

      bool use_refined_timestepping = true;
      if (main_db->keyExists("timestepping")) {
         std::string timestepping_method = main_db->getString("timestepping");
         if (timestepping_method == "SYNCHRONIZED") {
            use_refined_timestepping = false;
         }
      }

      const bool write_restart = (restart_interval > 0)
         && !(restart_write_dirname.empty());

      /*
       * Create major algorithm and data objects which comprise application.
       * Each object will be initialized either from input data or restart
       * files, or a combination of both.  Refer to each class constructor
       * for details.  For more information on the composition of objects
       * for this application, see comments at top of file.
       */

      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(dim,
            "CartesianGeometry",
            scaled_input_db->getDatabase("CartesianGeometry")));

      std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            "PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      std::shared_ptr<SinusoidalFrontGenerator> sine_wall;
      std::shared_ptr<SphericalShellGenerator> spherical_shell;
      std::shared_ptr<MeshGenerationStrategy> mesh_gen;

      if (input_db->isDatabase("SinusoidalFrontGenerator")) {
         sine_wall.reset(new SinusoidalFrontGenerator(
               "SinusoidalFrontGenerator", dim,
               input_db->getDatabase("SinusoidalFrontGenerator")));
         sine_wall->resetHierarchyConfiguration(
            patch_hierarchy, 0, patch_hierarchy->getMaxNumberOfLevels() - 1);
         mesh_gen = sine_wall;
      } else if (input_db->isDatabase("SphericalShellGenerator")) {
         spherical_shell.reset(new SphericalShellGenerator(
               "SphericalShellGenerator", dim,
               input_db->getDatabase("SphericalShellGenerator")));
         spherical_shell->resetHierarchyConfiguration(
            patch_hierarchy, 0, patch_hierarchy->getMaxNumberOfLevels() - 1);
         mesh_gen = spherical_shell;
      }

      LinAdv* linear_advection_model = new LinAdv(
            "LinAdv",
            dim,
            input_db->getDatabase("LinAdv"),
            grid_geometry,
            mesh_gen);

      std::shared_ptr<tbox::Database> hli_db(
         scaled_input_db->isDatabase("HyperbolicLevelIntegrator") ?
         scaled_input_db->getDatabase("HyperbolicLevelIntegrator") :
         input_db->getDatabase("HyperbolicLevelIntegrator"));
      std::shared_ptr<algs::HyperbolicLevelIntegrator> hyp_level_integrator(
         new algs::HyperbolicLevelIntegrator(
            "HyperbolicLevelIntegrator",
            hli_db,
            linear_advection_model, use_refined_timestepping));

      std::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
         new mesh::StandardTagAndInitialize(
            "StandardTagAndInitialize",
            hyp_level_integrator.get(),
            input_db->getDatabase("StandardTagAndInitialize")));

      // Set up the clustering.

      const std::string clustering_type =
         main_db->getStringWithDefault("clustering_type", "BergerRigoutsos");

      std::shared_ptr<mesh::BoxGeneratorStrategy> box_generator;

      if (clustering_type == "BergerRigoutsos") {

         std::shared_ptr<Database> abr_db(
            input_db->getDatabase("BergerRigoutsos"));
         std::shared_ptr<mesh::BoxGeneratorStrategy> berger_rigoutsos(
            new mesh::BergerRigoutsos(dim, abr_db));
         box_generator = berger_rigoutsos;

      } else if (clustering_type == "TileClustering") {

         std::shared_ptr<Database> tc_db(
            input_db->getDatabase("TileClustering"));
         std::shared_ptr<mesh::BoxGeneratorStrategy> tile_clustering(
            new mesh::TileClustering(dim, tc_db));
         box_generator = tile_clustering;

      }

      // Set up the load balancer.

      std::shared_ptr<mesh::LoadBalanceStrategy> load_balancer;
      std::shared_ptr<mesh::LoadBalanceStrategy> load_balancer0;

      const std::string load_balancer_type =
         main_db->getStringWithDefault("load_balancer_type", "TreeLoadBalancer");

      if (load_balancer_type == "TreeLoadBalancer") {

         std::shared_ptr<mesh::TreeLoadBalancer> tree_load_balancer(
            new mesh::TreeLoadBalancer(
               dim,
               "mesh::TreeLoadBalancer",
               input_db->getDatabase("TreeLoadBalancer"),
               std::shared_ptr<tbox::RankTreeStrategy>(new BalancedDepthFirstTree)));
         tree_load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

         std::shared_ptr<mesh::TreeLoadBalancer> tree_load_balancer0(
            new mesh::TreeLoadBalancer(
               dim,
               "mesh::TreeLoadBalancer0",
               input_db->getDatabase("TreeLoadBalancer"),
               std::shared_ptr<tbox::RankTreeStrategy>(new BalancedDepthFirstTree)));
         tree_load_balancer0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

         load_balancer = tree_load_balancer;
         load_balancer0 = tree_load_balancer0;
      } else if (load_balancer_type == "CascadePartitioner") {

         std::shared_ptr<mesh::CascadePartitioner> cascade_partitioner(
            new mesh::CascadePartitioner(
               dim,
               "mesh::CascadePartitioner",
               input_db->getDatabase("CascadePartitioner")));
         cascade_partitioner->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

         std::shared_ptr<mesh::CascadePartitioner> cascade_partitioner0(
            new mesh::CascadePartitioner(
               dim,
               "mesh::CascadePartitioner0",
               input_db->getDatabase("CascadePartitioner")));
         cascade_partitioner0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

         load_balancer = cascade_partitioner;
         load_balancer0 = cascade_partitioner0;
      } else if (load_balancer_type == "ChopAndPackLoadBalancer") {

         std::shared_ptr<mesh::ChopAndPackLoadBalancer> cap_load_balancer(
            new mesh::ChopAndPackLoadBalancer(
               dim,
               "mesh::ChopAndPackLoadBalancer",
               input_db->getDatabase("ChopAndPackLoadBalancer")));

         load_balancer = cap_load_balancer;

         /*
          * ChopAndPackLoadBalancer has trouble on L0 for some reason.
          * Work around by using the CascadePartitioner for L0.
          */
         std::shared_ptr<mesh::CascadePartitioner> cascade_partitioner0(
            new mesh::CascadePartitioner(
               dim,
               "mesh::CascadePartitioner0",
               input_db->getDatabase("CascadePartitioner")));
         cascade_partitioner0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
         load_balancer0 = cascade_partitioner0;
      }

      std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            patch_hierarchy,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_detector,
            box_generator,
            load_balancer,
            load_balancer0));

      std::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
         new algs::TimeRefinementIntegrator(
            "TimeRefinementIntegrator",
            input_db->getDatabase("TimeRefinementIntegrator"),
            patch_hierarchy,
            hyp_level_integrator,
            gridding_algorithm));

      /*
       * Initialize hierarchy configuration and data on all patches.
       * Then, close restart file and write initial state for visualization.
       */

      tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier(); // For timing.
      double dt_now = time_integrator->initializeHierarchy();

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

      RestartManager::getManager()->closeRestartFile();

      /*
       * After creating all objects and initializing their state, we
       * print the input database and variable database contents
       * to the log file.
       */

      if (mpi.getRank() == 0) {
         plog << "\nCheck input data and variables before simulation:" << std::endl;
         plog << "Input database..." << std::endl;
         input_db->printClassData(plog);
         plog << "\nVariable database..." << std::endl;
         hier::VariableDatabase::getDatabase()->printClassData(plog);

         plog << "\nCheck Linear Advection data... " << std::endl;
         linear_advection_model->printClassData(plog);
      }

      if (viz_dump_data) {
#ifdef HAVE_HDF5
         t_vis_writing->start();
         visit_data_writer->writePlotData(
            patch_hierarchy,
            time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime());
         t_vis_writing->stop();
#endif
      }

      /*
       * Time step loop.  Note that the step count and integration
       * time are maintained by algs::TimeRefinementIntegrator.
       */

      double loop_time = time_integrator->getIntegratorTime();
      double loop_time_end = time_integrator->getEndTime();

      int iteration_num = time_integrator->getIntegratorStep();

      while ((loop_time < loop_time_end) &&
             time_integrator->stepsRemaining()) {

         iteration_num = time_integrator->getIntegratorStep() + 1;

         pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
         pout << "At begining of timestep # " << iteration_num - 1 << std::endl;
         pout << "Simulation time is " << loop_time << std::endl;

         double dt_new = time_integrator->advanceHierarchy(dt_now);

         loop_time += dt_now;
         dt_now = dt_new;

         if (0) {
            /*
             * Logging can be very slow on I/O limited machines (such as
             * BlueGene).
             */
            plog << "Hierarchy summary:\n";
            patch_hierarchy->recursivePrint(plog, "H-> ", 1);
            plog << "PatchHierarchy summary:\n";
            patch_hierarchy->recursivePrint(plog, "L->", 1);
         }

         pout << "At end of timestep # " << iteration_num - 1 << std::endl;
         pout << "Simulation time is " << loop_time << std::endl;
         pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

         /*
          * At specified intervals, write restart and visualization files.
          */
         if (write_restart) {

            if ((iteration_num % restart_interval) == 0) {
               RestartManager::getManager()->
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
               t_vis_writing->start();
               visit_data_writer->writePlotData(patch_hierarchy,
                  iteration_num,
                  loop_time);
               t_vis_writing->stop();
#endif
            }
         }

#ifdef RECORD_STATS
         /*
          * Output statistics.
          */
         tbox::plog << "HyperbolicLevelIntegrator statistics:" << std::endl;
         hyp_level_integrator->printStatistics(tbox::plog);
         tbox::plog << "\nGriddingAlgorithm statistics:" << std::endl;
         gridding_algorithm->printStatistics(tbox::plog);
#endif
      } // End time-stepping loop.

      /*
       * Output timer results.
       */
      tbox::TimerManager::getManager()->print(tbox::plog);

      if (load_balancer_type == "TreeLoadBalancer") {
         /*
          * Output load balancing results for TreeLoadBalancer.
          */
         std::shared_ptr<mesh::TreeLoadBalancer> tree_load_balancer(
            SAMRAI_SHARED_PTR_CAST<mesh::TreeLoadBalancer, mesh::LoadBalanceStrategy>(
               load_balancer));
         TBOX_ASSERT(tree_load_balancer);
         tbox::plog << "\n\nLoad balancing results:\n";
         tree_load_balancer->printStatistics(tbox::plog);
      }
      if (load_balancer_type == "CascadePartitioner") {
         /*
          * Output load balancing results for CascadePartitioner.
          */
         std::shared_ptr<mesh::CascadePartitioner> cascade_partitioner(
            SAMRAI_SHARED_PTR_CAST<mesh::CascadePartitioner, mesh::LoadBalanceStrategy>(
               load_balancer));
         TBOX_ASSERT(cascade_partitioner);
         tbox::plog << "\n\nLoad balancing results:\n";
         cascade_partitioner->printStatistics(tbox::plog);
      }

      /*
       * Output box search results.
       */
      tbox::plog << "\n\nBox searching results:\n";
      hier::BoxTree::printStatistics(dim);

      t_all->stop();
      int size = tbox::SAMRAI_MPI::getSAMRAIWorld().getSize();
      if (tbox::SAMRAI_MPI::getSAMRAIWorld().getRank() == 0) {
         std::string timing_file =
            base_name + ".timing" + tbox::Utilities::intToString(size);
         FILE* fp = fopen(timing_file.c_str(), "w");
         fprintf(fp, "%f\n", t_all->getTotalWallclockTime());
         fclose(fp);
      }

      /*
       * At conclusion of simulation, deallocate objects.
       */

#ifdef HAVE_HDF5
      visit_data_writer.reset();
#endif

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

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  LinAdv" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}
