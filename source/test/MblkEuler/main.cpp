/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Run multiblock Euler AMR
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// Headers for basic SAMRAI objects

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/GridGeometry.h"

// Headers for major algorithm/data structure objects

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"

// Header for application-specific algorithm/data structure object

#include "test/testlib/MblkHyperbolicLevelIntegrator.h"
#include "MblkEuler.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace SAMRAI;

void
setupHierarchy(
   std::shared_ptr<tbox::Database> main_input_db,
   const tbox::Dimension& dim,
   std::shared_ptr<hier::BaseGridGeometry>& geometry,
   std::shared_ptr<hier::PatchHierarchy>& mblk_hierarchy);

//
// ===================================== The main code =======================
//

int main(
   int argc,
   char* argv[])
{
   //
   // initialize startup
   //
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

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

   //
   // fire up the log file
   //
   std::string log_file_name = "MblkEuler.log";
   bool log_all_nodes = false;
   if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
   } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
   }

#ifdef _OPENMP
   tbox::plog << "Compiled with OpenMP version " << _OPENMP
              << ".  Running with " << omp_get_max_threads() << " threads."
              << std::endl;
#else
   tbox::plog << "Compiled without OpenMP.\n";
#endif

   tbox::plog << "input_filename       = " << input_filename << std::endl;
   tbox::plog << "restart_read_dirname = " << restart_read_dirname << std::endl;
   tbox::plog << "restore_num          = " << restore_num << std::endl;

   //
   // Create input database and parse all data in input file.
   //

   std::shared_ptr<tbox::InputDatabase> input_db(
      new tbox::InputDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   tbox::plog << "---- done parsing input file" << std::endl << std::endl;

   //
   // Retrieve "GlobalInputs" section of the input database and set
   // values accordingly.
   //
   if (input_db->keyExists("GlobalInputs")) {
      std::shared_ptr<tbox::Database> global_db(
         input_db->getDatabase("GlobalInputs"));
//      if (global_db->keyExists("tag_clustering_method")) {
//       std::string tag_clustering_method =
//          global_db->getString("tag_clustering_method");
//       mesh::BergerRigoutsos::setClusteringOption(tag_clustering_method);
//     }
      if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
         bool flag = global_db->
            getBool("call_abort_in_serial_instead_of_exit");
         tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
      }
   }

   //
   // Retrieve "Main" section of the input database.  First, read dump
   // information, which is used for writing plot files.  Second, if
   // proper restart information was given on command line, and the
   // restart interval is non-zero, create a restart database.
   //
   std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

   const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

   //
   //..Initialize VisIt
   //
   int viz_dump_interval = 0;
   if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
   }

   std::string viz_dump_dirname = "";
   std::string visit_dump_dirname = "";
   int visit_number_procs_per_file = 1;

   if (viz_dump_interval > 0) {
      if (main_db->keyExists("viz_dump_dirname")) {
         viz_dump_dirname = main_db->getString("viz_dump_dirname");
      }
      visit_dump_dirname = viz_dump_dirname;
      if (viz_dump_dirname.empty()) {
         TBOX_ERROR("main(): "
            << "\nviz_dump_dirname is null ... "
            << "\nThis must be specified for use with VisIt"
            << std::endl);
      }
      if (main_db->keyExists("visit_number_procs_per_file")) {
         visit_number_procs_per_file =
            main_db->getInteger("visit_number_procs_per_file");
      }
   }

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

   //
   // Get restart manager and root restart database.  If run is from
   // restart, open the restart file.
   //
   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

   if (is_from_restart) {
      restart_manager->
      openRestartFile(restart_read_dirname, restore_num,
         mpi.getSize());
   }

   //
   // Setup the timer manager to trace timing statistics during execution
   // of the code.  The list of timers is given in the TimerManager
   // section of the input file.  Timing information is stored in the
   // restart file.  Timers will automatically be initialized to their
   // previous state if the run is restarted, unless they are explicitly
   // reset using the TimerManager::resetAllTimers() routine.
   //

   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   //
   // CREATE THE MULTIBLOCK HIERARCHY
   //
   std::shared_ptr<hier::PatchHierarchy> mblk_patch_hierarchy;
   std::shared_ptr<hier::BaseGridGeometry> geom;

   setupHierarchy(input_db,
      dim,
      geom,
      mblk_patch_hierarchy);

   //
   // -------------------- the patch operations --------------
   //
   MblkEuler* euler_model = new MblkEuler("MblkEuler",
         dim,
         input_db,
         geom);

   //
   // -------------------- the multiphase level operations --------------
   //
   std::shared_ptr<MblkHyperbolicLevelIntegrator> mblk_hyp_level_integrator(
      new MblkHyperbolicLevelIntegrator(
         "HyperbolicLevelIntegrator",
         dim,
         input_db->getDatabase("HyperbolicLevelIntegrator"),
         euler_model,
         mblk_patch_hierarchy,
         use_refined_timestepping));

   //
   // -------------------- the mesh refinement operations --------------
   //
   std::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
      new mesh::StandardTagAndInitialize(
         "StandardTagAndInitialize",
         mblk_hyp_level_integrator.get(),
         input_db->getDatabase("StandardTagAndInitialize")));

   std::shared_ptr<mesh::BergerRigoutsos> box_generator(
      new mesh::BergerRigoutsos(dim,
         input_db->getDatabase("BergerRigoutsos")));

   std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
      new mesh::TreeLoadBalancer(
         dim,
         "TreeLoadBalancer",
         input_db->getDatabase("TreeLoadBalancer"),
         std::shared_ptr<tbox::RankTreeStrategy>(new tbox::BalancedDepthFirstTree)));
   load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

   std::shared_ptr<mesh::GriddingAlgorithm> mblk_gridding_algorithm(
      new mesh::GriddingAlgorithm(
         mblk_patch_hierarchy,
         "GriddingAlgorithm",
         input_db->getDatabase("GriddingAlgorithm"),
         error_detector,
         box_generator,
         load_balancer,
         load_balancer));

   std::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
      new algs::TimeRefinementIntegrator(
         "TimeRefinementIntegrator",
         input_db->getDatabase("TimeRefinementIntegrator"),
         mblk_patch_hierarchy,
         mblk_hyp_level_integrator,
         mblk_gridding_algorithm));

#ifdef HAVE_HDF5
   //
   // ----------------------------- Set up Visualization writer(s).
   //
   bool is_multiblock = true;
   std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
      new appu::VisItDataWriter(
         dim,
         "MblkEuler VisIt Writer",
         visit_dump_dirname,
         visit_number_procs_per_file,
         is_multiblock));
   euler_model->
   registerVisItDataWriter(visit_data_writer);
#endif

   //
   // Initialize hierarchy configuration and data on all patches.
   // Then, close restart file and write initial state for visualization.
   //
   double dt_now = time_integrator->initializeHierarchy();

   tbox::RestartManager::getManager()->closeRestartFile();

   //
   // After creating all objects and initializing their state, we
   // print the input database and variable database contents
   // to the log file.
   //
   tbox::plog << "\nCheck input data and variables before simulation:" << std::endl;
   tbox::plog << "Input database..." << std::endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << std::endl;
   hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

   tbox::plog << "\nPrinting a summary of model input... " << std::endl;
   euler_model->printClassData(tbox::plog);

#ifdef HAVE_HDF5
   if (viz_dump_data) {
      visit_data_writer->writePlotData(
         mblk_patch_hierarchy,
         time_integrator->getIntegratorStep(),
         time_integrator->getIntegratorTime());
   }
#endif

   //
   // ==============================================================
   // Time step loop.  Note that the step count and integration
   // time are maintained by TimeRefinementIntegrator.
   // ==============================================================
   //

   double loop_time = time_integrator->getIntegratorTime();
   double loop_time_end = time_integrator->getEndTime();

   int iteration_num = time_integrator->getIntegratorStep();

   int old_log_style = 1;

   while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {

      iteration_num = time_integrator->getIntegratorStep() + 1;

      if (old_log_style) {
         tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
         tbox::pout << "At begining of timestep # " << iteration_num - 1
                    << std::endl;
         tbox::pout << "Simulation time is " << loop_time << std::endl;
      }

      //
      // advance the heirarchy a timestep
      //
      double dt_new = time_integrator->advanceHierarchy(dt_now);

      loop_time += dt_now;
      dt_now = dt_new;

      if (!old_log_style) {
         //
         // write out the timestep header
         //
         char my_line[256];
         snprintf(my_line, 256, "%4d time: %9.5e dt: %10.6e  ",
            iteration_num,
            loop_time,
            dt_new);

         tbox::pout << my_line << std::endl;
      }

      if (old_log_style) {
         tbox::pout << "At end of timestep # " << iteration_num - 1 << std::endl;
         tbox::pout << "Simulation time is " << loop_time << std::endl;
         tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      }

      //
      // At specified intervals, write restart files.
      //
      if (write_restart) {

         if ((iteration_num % restart_interval) == 0) {
            tbox::RestartManager::getManager()->
            writeRestartFile(restart_write_dirname,
               iteration_num);
         }
      }

#ifdef HAVE_HDF5
      //
      // At specified intervals, write out data files for plotting.
      //
      if (viz_dump_data) {
         if ((iteration_num % viz_dump_interval) == 0) {
            visit_data_writer->writePlotData(mblk_patch_hierarchy,
               iteration_num,
               loop_time);
         }
      }
#endif

   }   //-----------------------------------------------END TIME STEPPING LOOP

   //
   // Output timer results.
   //
   tbox::TimerManager::getManager()->print(tbox::plog);

   //
   // At conclusion of simulation, deallocate objects.
   //
#ifdef HAVE_HDF5
   visit_data_writer.reset();
#endif
   time_integrator.reset();
   mblk_gridding_algorithm.reset();
   load_balancer.reset();
   box_generator.reset();
   error_detector.reset();
   mblk_hyp_level_integrator.reset();

   if (euler_model) delete euler_model;

   mblk_patch_hierarchy.reset();
   geom.reset();

   input_db.reset();
   main_db.reset();

   tbox::pout << "\nPASSED:  MblkEuler" << std::endl;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}

// ----------------------------------------------------------------

//
// this function builds the skeleton grid geometry
//
void setupHierarchy(
   std::shared_ptr<tbox::Database> main_input_db,
   const tbox::Dimension& dim,
   std::shared_ptr<hier::BaseGridGeometry>& geometry,
   std::shared_ptr<hier::PatchHierarchy>& mblk_hierarchy)
{
   TBOX_ASSERT(main_input_db);

   std::shared_ptr<tbox::Database> mult_db(
      main_input_db->getDatabase("PatchHierarchy"));

   /*
    * Read the geometry information and build array of geometries
    */


   std::string geom_name("BlockGeometry");
   if (main_input_db->keyExists(geom_name)) {
      geometry.reset(
         new geom::GridGeometry(
            dim,
            geom_name,
            main_input_db->getDatabase(geom_name)));
   } else {
      TBOX_ERROR("main::setupHierarchy(): could not find entry `"
         << geom_name << "' in input.");
   }

   mblk_hierarchy.reset(
      new hier::PatchHierarchy(
         "PatchHierarchy",
         geometry,
         mult_db));

}
