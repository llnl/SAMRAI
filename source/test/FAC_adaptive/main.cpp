/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Program for poisson solver on adaptive grid using FAC
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include IOMANIP_HEADER_FILE
#include <fstream>

#include "AdaptivePoisson.h"
#include "test/testlib/get-input-filename.h"

/*
 * Headers for basic SAMRAI objects used in this code.
 */
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 * Headers for major algorithm/data structure objects from SAMRAI
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/solv/FACPreconditioner.h"


#include <vector>
#include <string>
#include <memory>

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{

   std::string input_filename;

   /*
    * Initialize MPI, process argv, and initialize SAMRAI
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   if (get_input_filename(&argc, argv, input_filename) == 1) {
      tbox::pout << "Usage: " << argv[0] << " <input file>." << std::endl;
      tbox::SAMRAI_MPI::finalize();
      return 0;
   }
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   bool error_ok = false;

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      /*
       * These tests are set up to use hypre.
       * Do not run them without hypre.
       */
#ifdef HAVE_HYPRE
      tbox::pout << "Input file is " << input_filename << std::endl;

      std::string case_name;
      if (argc > 1) {
         case_name = argv[1];
      }

      /*
       * Create input database and parse all data in input file into it.
       */

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      if (input_db->isDatabase("TimerManager")) {
         tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Get the Main database part of the input database.
       * This database contains information relevant to main.
       */

      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      tbox::plog << "Main database:" << std::endl;
      main_db->printClassData(tbox::plog);

      /*
       * Base filename info.
       */

      std::string base_name =
         main_db->getStringWithDefault("base_name", "noname");

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      if (!case_name.empty()) {
         base_name = base_name + '-' + case_name;
      }
      base_name = base_name + '-'
         + tbox::Utilities::intToString(tbox::SAMRAI_MPI::getSAMRAIWorld().getSize(), 5);

      /*
       * Log file info.
       */
      {
         std::string log_filename =
            main_db->getStringWithDefault("log_filename", base_name + ".log");
         bool log_all =
            main_db->getBoolWithDefault("log_all", false);
         if (log_all)
            tbox::PIO::logAllNodes(log_filename);
         else
            tbox::PIO::logOnlyNodeZero(log_filename);
      }

      /*
       * Create a patch hierarchy for use later.
       * This object is a required input for these objects: adaptive_poisson.
       */
      /*
       * Create a grid geometry required for the patchHierarchy object.
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

      std::string adaptive_poisson_name = "AdaptivePoisson";
      std::string fac_ops_name =
         adaptive_poisson_name + ":scalar poisson operator";
      std::string fac_precond_name =
         "FAC preconditioner for Poisson's equation";
      std::string hypre_poisson_name = fac_ops_name + "::hypre_solver";

#ifdef HAVE_HYPRE
      std::shared_ptr<solv::CellPoissonHypreSolver> hypre_poisson(
         new solv::CellPoissonHypreSolver(
            dim,
            hypre_poisson_name,
            input_db->isDatabase("hypre_solver") ?
            input_db->getDatabase("hypre_solver") :
            std::shared_ptr<tbox::Database>()));

      std::shared_ptr<solv::CellPoissonFACOps> fac_ops(
         new solv::CellPoissonFACOps(
            hypre_poisson,
            dim,
            fac_ops_name,
            input_db->isDatabase("fac_ops") ?
            input_db->getDatabase("fac_ops") :
            std::shared_ptr<tbox::Database>()));
#else
      std::shared_ptr<solv::CellPoissonFACOps> fac_ops(
         new solv::CellPoissonFACOps(
            dim,
            fac_ops_name,
            input_db->isDatabase("fac_ops") ?
            input_db->getDatabase("fac_ops") :
            std::shared_ptr<tbox::Database>()));
#endif

      std::shared_ptr<solv::FACPreconditioner> fac_precond(
         new solv::FACPreconditioner(
            fac_precond_name,
            fac_ops,
            input_db->isDatabase("fac_precond") ?
            input_db->getDatabase("fac_precond") :
            std::shared_ptr<tbox::Database>()));

      AdaptivePoisson adaptive_poisson(adaptive_poisson_name,
                                       dim,
                                       fac_ops,
                                       fac_precond,
                                       *(input_db->getDatabase("AdaptivePoisson")),
                                       &tbox::plog);

      /*
       * Create the tag-and-initializer, box-generator and load-balancer
       * object references required by the gridding_algorithm object.
       */
      std::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
         new mesh::StandardTagAndInitialize(
            "CellTaggingMethod",
            &adaptive_poisson,
            input_db->getDatabase("StandardTagAndInitialize")));
      std::shared_ptr<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(
            dim,
            (input_db->isDatabase("BergerRigoutsos") ?
             input_db->getDatabase("BergerRigoutsos") :
             std::shared_ptr<tbox::Database>())));
      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "load balancer",
            input_db->getDatabase("TreeLoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      /*
       * Create the gridding algorithm used to generate the SAMR grid
       * and create the grid.
       */
      std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            patch_hierarchy,
            "Gridding Algorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            tag_and_initializer,
            box_generator,
            load_balancer));
      tbox::plog << "Gridding algorithm:" << std::endl;
      gridding_algorithm->printClassData(tbox::plog);
      /*
       * Make the coarse patch level.
       */
      gridding_algorithm->makeCoarsestLevel(0.0);

      int ln;

      /* Whether to plot */
      std::string vis_filename =
         main_db->getStringWithDefault("vis_filename", base_name);
      bool do_plot =
         main_db->getBoolWithDefault("do_plot", false);

      /*
       * After creating all objects and initializing their state,
       * we print the input database and variable database contents
       * to the log file.
       */
      tbox::plog << "\nCheck input data and variables before simulation:"
                 << std::endl;
      tbox::plog << "Input database..." << std::endl;
      input_db->printClassData(tbox::plog);
      tbox::plog << "\nVariable database..." << std::endl;
      hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

      tbox::plog << "\n\nFinal Hierarchy:\n";
      patch_hierarchy->recursivePrint(tbox::plog, "\t", 2);

      double target_l2norm = 1e-6;
      target_l2norm = main_db->getDoubleWithDefault("target_l2norm",
            target_l2norm);
      double l2norm, linorm;
      int max_adaptions = 1;
      max_adaptions = main_db->getIntegerWithDefault("max_adaptions",
            max_adaptions);
      int adaption_number = 0;
      bool done = false;
      do {
         /*
          * Solve.
          */
         tbox::pout.setf(std::ios::scientific);
         std::string initial_u =
            main_db->getStringWithDefault("initial_u", "0.0");
         adaptive_poisson.solvePoisson(patch_hierarchy,
            adaption_number ? std::string() : initial_u);
         std::vector<double> l2norms(patch_hierarchy->getNumberOfLevels());
         std::vector<double> linorms(patch_hierarchy->getNumberOfLevels());
         adaptive_poisson.computeError(*patch_hierarchy,
            &l2norm,
            &linorm,
            l2norms,
            linorms);
         error_ok = l2norm <= target_l2norm;
         tbox::plog << "Err " << (error_ok ? "" : "NOT ")
                    << "ok, err norm/target: "
                    << std::scientific << l2norm << '/' << std::scientific
                    << target_l2norm << std::endl;
         tbox::plog << "Err result after " << adaption_number
                    << " adaptions: \n"
                    << std::setw(15) << "l2: " << std::setw(10) << std::scientific << l2norm
                    << std::setw(15) << "li: " << std::setw(10) << std::scientific << linorm
                    << "\n";
         for (ln = 0; ln < patch_hierarchy->getNumberOfLevels(); ++ln) {
            tbox::plog << std::setw(10) << "l2[" << std::setw(2) << ln << "]: "
                       << std::setw(10) << std::scientific << l2norms[ln]
                       << std::setw(10) << "li[" << std::setw(2) << ln << "]: "
                       << std::setw(10) << std::scientific << linorms[ln]
                       << "\n";
         }

         /* Write the plot file. */
         if (do_plot) {
            std::shared_ptr<appu::VisItDataWriter> visit_writer(
               new appu::VisItDataWriter(
                  dim,
                  "VisIt Writer",
                  vis_filename + ".visit"));
            adaptive_poisson.registerVariablesWithPlotter(*visit_writer);
            visit_writer->writePlotData(patch_hierarchy,
               adaption_number);
            tbox::plog << "Wrote viz file " << vis_filename
                       << " for grid number "
                       << adaption_number << '\n';
         }

         /*
          * Done when max adaptions or convergence reached.
          */
         done = error_ok || (adaption_number >= max_adaptions);

         if (!done) {
            /*
             * Adapt grid.
             */
            ++adaption_number;
            tbox::plog << "Adaption number " << adaption_number << "\n";

            std::vector<int> tag_buffer(patch_hierarchy->getMaxNumberOfLevels());
            for (ln = 0; ln < static_cast<int>(tag_buffer.size()); ++ln) {
               tag_buffer[ln] = 1;
            }
            gridding_algorithm->regridAllFinerLevels(
               0,
               tag_buffer,
               0,
               0.0);
            tbox::plog << "Newly adapted hierarchy\n";
            patch_hierarchy->recursivePrint(tbox::plog, "    ", 1);
            if (0) {
               /* Write post-adapt viz file for debugging */
               std::shared_ptr<appu::VisItDataWriter> visit_writer(
                  new appu::VisItDataWriter(
                     dim,
                     "VisIt Writer",
                     "postadapt.visit"));
               adaptive_poisson.registerVariablesWithPlotter(*visit_writer);
               visit_writer->writePlotData(patch_hierarchy,
                  adaption_number - 1);
               tbox::plog << "Wrote viz file " << "postadapt.visit" << '\n';
            }
         }
      } while (!done);

      tbox::plog << "After " << adaption_number << "/" << max_adaptions
                 << " adaptions, residual is " << l2norm << "/"
                 << target_l2norm
                 << std::endl;

      tbox::TimerManager::getManager()->print(tbox::plog);

#else
      error_ok = true;
#endif

      if (error_ok) {
         tbox::pout << "\nPASSED:  FAC" << std::endl;
      } else {
         TBOX_ERROR("Failed to meet accuracy specifications.");
      }
   }

   /*
    * Exit properly by shutting down services in correct order.
    */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
