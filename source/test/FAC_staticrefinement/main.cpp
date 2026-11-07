/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for FAC Poisson example
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <memory>

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"

#include "FACPoisson.h"


using namespace SAMRAI;

/*
 ************************************************************************
 *
 * This is the driver program to demonstrate
 * how to use the FAC Poisson solver.
 *
 * We set up the simple problem
 *          u + div(grad(u)) = sin(x)*sin(y)
 * in the domain [0:1]x[0:1], with u=0 on the
 * boundary.
 *
 * FACPoisson is the primary object used to
 * set up and solve the system.  It maintains
 * the data for the computed solution u, the
 * exact solution, and the right hand side.
 *
 * The hierarchy created to solve this problem
 * has only one level.  (The FAC Poisson solver
 * is a single-level solver.)
 *
 *************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
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

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Set up the timer manager.
       */
      if (input_db->isDatabase("TimerManager")) {
         tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

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
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
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
            base_name + "CartesianGridGeometry",
            input_db->getDatabase("CartesianGridGeometry")));
      tbox::plog << "Cartesian Geometry:" << std::endl;
      grid_geometry->printClassData(tbox::plog);

      std::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            base_name + "::PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      std::string fac_poisson_name = base_name + "::FACPoisson";
      std::string fac_solver_name = fac_poisson_name + "::poisson_hypre";
      std::string fac_ops_name = fac_solver_name + "::fac_ops";
      std::string fac_precond_name = fac_solver_name + "::fac_precond";
      std::string hypre_poisson_name = fac_ops_name + "::hypre_solver";
      std::string bc_coefs_name = fac_poisson_name + "::bc_coefs";

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

      std::shared_ptr<solv::CellPoissonFACSolver> fac_solver(
         new solv::CellPoissonFACSolver(
            dim,
            fac_solver_name,
            fac_precond,
            fac_ops,
            input_db->isDatabase("fac_solver") ?
            input_db->getDatabase("fac_solver") :
            std::shared_ptr<tbox::Database>()));

      std::shared_ptr<solv::LocationIndexRobinBcCoefs> bc_coefs(
         new solv::LocationIndexRobinBcCoefs(
            dim,
            bc_coefs_name,
            input_db->isDatabase("bc_coefs") ?
            input_db->getDatabase("bc_coefs") :
            std::shared_ptr<tbox::Database>()));

      /*
       * The FACPoisson object is the main user object specific to the
       * problem being solved.  It provides the implementations for setting
       * up the grid and plotting data.  It also wraps up the solve
       * process that includes making the initial guess, specifying the
       * boundary conditions and call the solver.
       */
      FACPoisson fac_poisson(fac_poisson_name,
                             dim,
                             fac_solver,
                             bc_coefs);

      /*
       * Create the tag-and-initializer, box-generator and load-balancer
       * object references required by the gridding_algorithm object.
       */
      std::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
         new mesh::StandardTagAndInitialize(
            "CellTaggingMethod",
            &fac_poisson,
            input_db->getDatabase("StandardTagAndInitialize")));
      std::shared_ptr<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(dim));
      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "load balancer",
            std::shared_ptr<tbox::Database>()));
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
       * Make the coarsest patch level where we will be solving.
       */
      gridding_algorithm->makeCoarsestLevel(0.0);
      bool done = false;
      for (int lnum = 0;
           patch_hierarchy->levelCanBeRefined(lnum) && !done; lnum++) {
         tbox::plog << "Adding finner levels with lnum = " << lnum << std::endl;
         gridding_algorithm->makeFinerLevel(
            0,
            true,
            0,
            0.0);
         tbox::plog << "Just added finer levels with lnum = " << lnum << std::endl;
         done = !(patch_hierarchy->finerLevelExists(lnum));
      }

      /*
       * Set up the plotter for the hierarchy just created.
       * The FACPoisson object handles the data and has the
       * function setupExternalPlotter to register its data
       * with the plotter.
       */
#ifdef HAVE_HDF5
      std::string vis_filename =
         main_db->getStringWithDefault("vis_filename", base_name);
      std::shared_ptr<appu::VisItDataWriter> visit_writer(
         std::make_shared<appu::VisItDataWriter>(dim,
                                                   "VisIt Writer",
                                                   vis_filename + ".visit"));
      fac_poisson.setupPlotter(*visit_writer);
#endif

      /*
       * After creating all objects and initializing their state,
       * we print the input database and variable database contents
       * to the log file.
       */
      tbox::plog << "\nCheck input data and variables before simulation:"
                 << std::endl;
      tbox::plog << "Input database..." << std::endl;
      input_db->printClassData(tbox::plog);

      /*
       * Solve.
       */
      fac_poisson.solvePoisson();

#ifdef HAVE_HDF5
      /*
       * Plot.
       */
      visit_writer->writePlotData(patch_hierarchy, 0);
#endif

      /*
       * Deallocate objects when done.
       */

      tbox::TimerManager::getManager()->print(tbox::plog);
   }

   /*
    * This print is for the SAMRAI testing framework.  Passing here
    * means application ran.  A better test would actually test the
    * results.
    */
   tbox::pout << "\nPASSED:  FAC" << std::endl;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
