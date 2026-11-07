/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for asynchronous BR implementation
 *
 ************************************************************************/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 * Headers for major algorithm/data structure objects from SAMRAI
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "ABRTest.h"

#include "test/testlib/get-input-filename.h"


#include <vector>

using namespace SAMRAI;

int main(
   int argc,
   char** argv)
{

   std::string input_filename;

   /*
    * Initialize MPI, process argv, and initialize SAMRAI
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   if (get_input_filename(&argc, argv, input_filename) == 1) {
      std::cout << "Usage: " << argv[0]
           << " <input file>."
           << std::endl;
      tbox::SAMRAI_MPI::finalize();
      return 0;
   }
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::plog << "Input file is " << input_filename << std::endl;

      std::string case_name;
      if (argc >= 2) {
         case_name = argv[1];
      }

      /*
       * Create input database and parse all data in input file into it.
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
         if (global_db->keyExists("call_abort_in_serial_instead_of_exit")) {
            bool flag = global_db->
               getBool("call_abort_in_serial_instead_of_exit");
            tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(flag);
         }
      }

      /*
       * Get the Main database part of the input database.
       * This database contains information relevant to main.
       */

      std::shared_ptr<tbox::Database> main_db(
         input_db->getDatabase("Main"));
      tbox::plog << "Main database:" << std::endl;
      main_db->printClassData(tbox::plog);

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      if (input_db->isDatabase("TimerManager")) {
         tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Base filename info.
       */

      std::string base_name = main_db->getStringWithDefault("base_name", "fp");

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      if (!case_name.empty()) {
         base_name = base_name + '-' + case_name;
      }
      if (mpi.getSize() > 1) {
         base_name = base_name + '-'
            + tbox::Utilities::intToString(mpi.getSize(), 5);
      }
      tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                 << mpi.getSize() << ") to base name -> '"
                 << base_name << "'\n";

      /*
       * Set the vis filename, defaults to base_name.
       */
      std::string vis_filename =
         main_db->getStringWithDefault("vis_filename", base_name);

      /*
       * Log file info.
       */

      std::string log_filename =
         main_db->getStringWithDefault("log_filename", base_name + ".log");
      bool log_all = false;
      log_all = main_db->getBoolWithDefault("log_all", log_all);
      if (log_all && mpi.getSize() > 1) {
         tbox::PIO::logAllNodes(log_filename);
      } else {
         tbox::PIO::logOnlyNodeZero(log_filename);
      }

      if (!case_name.empty()) {
         tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                    << mpi.getSize() << ") to base name -> '"
                    << base_name << "'\n";
      }
      tbox::plog << "Running on " << mpi.getSize()
                 << " processes.\n";

      int plot_step = main_db->getIntegerWithDefault("plot_step", 0);

      /*
       * Create a patch hierarchy for use later.
       * This object is a required input for these objects: abrtest.
       */
      /*
       * Create a grid geometry required for the
       * hier::PatchHierarchy object.
       */
      std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            dim,
            "CartesianGridGeometry",
            input_db->getDatabase("CartesianGridGeometry")));
      tbox::plog << "Grid Geometry:" << std::endl;
      grid_geometry->printClassData(tbox::plog);
      std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            "Patch Hierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      /*
       * Create the problem-specific object implementing the required
       * SAMRAI virtual functions.
       */
      tbox::plog << "Creating abrtest.\n";
      ABRTest abrtest("ABRTest",
                      dim,
                      patch_hierarchy,
                      input_db->getDatabase("ABRTest"));

      tbox::plog << "Creating box generator.\n";
      std::shared_ptr<mesh::BergerRigoutsos> new_br(
         new mesh::BergerRigoutsos(
            dim,
            input_db->isDatabase("BergerRigoutsos") ?
            input_db->getDatabase("BergerRigoutsos") :
            std::shared_ptr<tbox::Database>()));
      new_br->useDuplicateMPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      tbox::plog << "Creating grid algorithm.\n";
      /*
       * Create the tag-and-initializer, box-generator and load-balancer
       * object references required by the gridding_algorithm object.
       */
      std::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
         new mesh::StandardTagAndInitialize(
            "CellTaggingMethod",
            abrtest.getStandardTagAndInitObject(),
            input_db->getDatabase("StandardTagAndInitialize")));
      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "tree load balancer",
            input_db->getDatabase("TreeLoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      /*
       * Create the gridding algorithm used to generate the SAMR grid
       * and create the grid.
       */
      std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            patch_hierarchy,
            "Distributed Gridding Algorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            tag_and_initializer,
            new_br,
            load_balancer));
      tbox::plog << "Sistributed gridding algorithm:" << std::endl;
      gridding_algorithm->printClassData(tbox::plog);

      bool log_hierarchy = false;
      log_hierarchy = main_db->getBoolWithDefault("log_hierarchy",
            log_hierarchy);
      int num_steps = main_db->getIntegerWithDefault("num_steps", 0);

      /*
       * After setting up the problem and initializing the object states,
       * we print the input database and variable database contents
       * to the log file.
       */
      tbox::plog << "\nCheck input data:" << std::endl;
      tbox::plog << "Input database..." << std::endl;
      input_db->printClassData(tbox::plog);
      tbox::plog << "\nVariable database..." << std::endl;
      hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

      tbox::plog
      << "**********************************************************\n";
      tbox::plog << "Memory used before mesh generation:" << std::endl;
      tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
      tbox::plog
      << "**********************************************************\n";

      /*
       * Make the patch levels.
       */

      std::shared_ptr<tbox::Timer> t_generate_mesh(
         tbox::TimerManager::getManager()->
         getTimer("apps::main::generate_mesh"));
      t_generate_mesh->start();
      gridding_algorithm->makeCoarsestLevel(0.0);
      tbox::plog << "Memory used after creating level 0:" << std::endl;
      tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
      bool done = false;
      for (int ln = 0; patch_hierarchy->levelCanBeRefined(ln) && !done;
           ++ln) {
         tbox::plog << "Adding finer levels with ln = " << ln << std::endl;
         std::shared_ptr<hier::PatchLevel> level_(
            patch_hierarchy->getPatchLevel(ln));
         gridding_algorithm->makeFinerLevel(
            /* tag buffer size */ 0,
            /* whether initial cycle */ true,
            /* cycle */ 0,
            /* simulation time */ 0.0);
         tbox::plog << "Just added finer level " << ln << " -> " << ln + 1;
         if (patch_hierarchy->getNumberOfLevels() < ln + 2) {
            tbox::plog << " (no new level!)" << std::endl;
         } else {
            std::shared_ptr<hier::PatchLevel> finer_level_(
               patch_hierarchy->getPatchLevel(ln + 1));
            tbox::plog
            << " (" << level_->getNumberOfPatches()
            << " -> " << finer_level_->getNumberOfPatches()
            << " patches)"
            << std::endl;
         }
         done = !(patch_hierarchy->finerLevelExists(ln));

         tbox::plog << "Memory used after creating level " << ln + 1 << ":"
                    << std::endl;
         tbox::MemoryUtilities::printMemoryInfo(tbox::plog);

      }
      t_generate_mesh->stop();

      if (mpi.getRank() == 0) {
         tbox::plog << "Hierarchy generated:" << std::endl;
         patch_hierarchy->recursivePrint(tbox::plog, std::string("    "), 1);
      }
      if (log_hierarchy) {
         tbox::plog << "Hierarchy generated:" << std::endl;
         patch_hierarchy->recursivePrint(tbox::plog, std::string("H-> "), 3);
      }

#ifdef HAVE_HDF5
      /*
       * Write a plot file.
       */
      /* Get the output filename. */
      if (plot_step > 0) {
         const std::string visit_filename = vis_filename + ".visit";
         /* Create the VisIt data writer. */
         std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
            new appu::VisItDataWriter(
               dim,
               "VisIt Writer",
               visit_filename));
         /* Register variables with plotter. */
         abrtest.registerVariablesWithPlotter(visit_data_writer);
         /* Write the plot file. */
         visit_data_writer->writePlotData(patch_hierarchy, 0);
      }
#endif

      /*
       * Adapt the grid.
       */
      std::vector<int> tag_buffer(10);
      for (int i = 0; i < static_cast<int>(tag_buffer.size()); ++i) {
         tag_buffer[i] = 1;
      }

      for (int istep = 0; istep < num_steps; ++istep) {

         tbox::plog << "Adaption number " << istep << std::endl;

         // Recompute the front-dependent data at next time step.
         abrtest.computeHierarchyData(*patch_hierarchy,
            double(istep + 1));

         std::vector<double> regrid_start_time(patch_hierarchy->getMaxNumberOfLevels());
         for (int i = 0; i < static_cast<int>(regrid_start_time.size()); ++i) {
            regrid_start_time[i] = istep;
         }

         gridding_algorithm->regridAllFinerLevels(
            0,
            tag_buffer,
            istep + 1,
            double(istep + 1),
            regrid_start_time);

         if (mpi.getRank() == 0) {
            patch_hierarchy->recursivePrint(tbox::plog, std::string("    "), 1);
         }
         if (log_hierarchy) {
            tbox::plog << "Hierarchy adapted:" << std::endl;
            patch_hierarchy->recursivePrint(tbox::plog, std::string("H-> "), 3);
         }

         tbox::plog << "Memory used after adaption number " << istep << std::endl;
         tbox::MemoryUtilities::printMemoryInfo(tbox::plog);

#ifdef HAVE_HDF5
         if (plot_step > 0 && (istep + 1) % plot_step == 0) {
            const std::string visit_filename = vis_filename + ".visit";
            /* Create the VisIt data writer. */
            std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
               new appu::VisItDataWriter(
                  dim,
                  "VisIt Writer",
                  visit_filename));
            /* Register variables with plotter. */
            abrtest.registerVariablesWithPlotter(visit_data_writer);
            /* Write the plot file. */
            visit_data_writer->writePlotData(patch_hierarchy, istep + 1);
         }
#endif
      }

      tbox::TimerManager::getManager()->print(tbox::plog);

      tbox::pout << "\nPASSED:  async_br" << std::endl;

   }

   /*
    * Exit properly by shutting down services in correct order.
    */
   tbox::plog << "\nShutting down..." << std::endl;
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
