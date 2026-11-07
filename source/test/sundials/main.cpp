/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for testing Sundials/SAMRAI interface.
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

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/CVODESolver.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "CVODEModel.h"

using namespace SAMRAI;

/*
 * The cvode_test program is a general skeleton for using the
 * CVODESolver interface.  The main stages of this program are:
 *
 * (1)  Retrieving integration parameters from the input database.
 * (2)  Creating hierarchy, geometry, gridding, and CVODEModel
 *      objects.
 * (3)  Setting up the hierarchy configuration (grid configuration).
 * (4)  Setting the initial condition vector.
 * (5)  Creating a CVODESolver object.
 * (6)  Setting the integration parameters for CVODESolver.
 * (7)  Printing out to the log file the initial condition vector
 *      and computing some norm for checking purposes.
 * (8)  Solving the ODE system.
 * (9)  Printing out to the log file the solution vector produced
 *      by CVODE and computing some norms.
 * (10) Printing out the CVODE statistics.
 * (11) Cleaning up the memory allocated for the program.
 */

int main(
   int argc,
   char* argv[])
{

   /*
    * Initialize tbox::MPI and SAMRAI.  Enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      tbox::PIO::logAllNodes("cvode_test.log");

#if !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)
      tbox::pout << "Library compiled WITHOUT CVODE -and- HYPRE...\n"
                 << "SAMRAI was not configured with one, or both, of "
                 << "these packages.  Cannot run this example." << std::endl;
#else

      /*
       * Process command line arguments.
       */
      std::string input_filename;

      if (argc != 2) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << std::endl;
         exit(-1);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */
      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /**************************************************************************
      * Read input data and setup objects.
      **************************************************************************/

      /*
       * Retreive "Main" section of input db.
       */
      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      int max_order = main_db->getInteger("max_order");
      int max_internal_steps = main_db->getInteger("max_internal_steps");
      double init_time = main_db->getDouble("init_time");
      int init_cycle = main_db->getInteger("init_cycle");
      double print_interval = main_db->getDouble("print_interval");
      int num_print_intervals = main_db->getInteger("num_print_intervals");

      double relative_tolerance = main_db->getDouble("relative_tolerance");
      double absolute_tolerance = main_db->getDouble("absolute_tolerance");
      int stepping_method = main_db->getInteger("stepping_method");
      bool uses_preconditioning =
         main_db->getBoolWithDefault("uses_preconditioning", false);
      bool solution_logging =
         main_db->getBoolWithDefault("solution_logging", false);

      /*
       * Create geometry and hierarchy objects.
       */
      std::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            dim,
            "Geometry",
            input_db->getDatabase("Geometry")));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy(
            "Hierarchy",
            geometry,
            input_db->getDatabase("PatchHierarchy")));

      /*
       * Create gridding algorithm objects that will handle construction of
       * of the patch levels in the hierarchy.
       */

      std::string cvode_model_name = "CVODEModel";
      std::string fac_solver_name = cvode_model_name + ":FAC solver";
      std::string fac_ops_name = fac_solver_name + "::fac_ops";
      std::string fac_precond_name = fac_solver_name + "::fac_precond";
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

      std::shared_ptr<solv::CellPoissonFACSolver> fac_solver(
         new solv::CellPoissonFACSolver(
            dim,
            fac_solver_name,
            fac_precond,
            fac_ops,
            input_db->isDatabase("fac_solver") ?
            input_db->getDatabase("fac_solver") :
            std::shared_ptr<tbox::Database>()));

      std::shared_ptr<CVODEModel> cvode_model(
         new CVODEModel(
            cvode_model_name,
            dim,
            fac_solver,
            input_db->getDatabase("CVODEModel"),
            geometry));

      std::shared_ptr<mesh::StandardTagAndInitialize> error_est(
         new mesh::StandardTagAndInitialize(
            "StandardTagAndInitialize",
            cvode_model.get(),
            input_db->getDatabase("StandardTagAndInitialize")));

      std::shared_ptr<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(dim,
            input_db->getDatabase("BergerRigoutsos")));

      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "LoadBalancer",
            input_db->getDatabase("LoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            hierarchy,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_est,
            box_generator,
            load_balancer));

      /*
       * Setup hierarchy.
       */
      gridding_algorithm->makeCoarsestLevel(init_time);

      std::vector<int> tag_buffer_array(hierarchy->getMaxNumberOfLevels());
      for (int il = 0; il < hierarchy->getMaxNumberOfLevels(); ++il) {
         tag_buffer_array[il] = 1;
      }

      bool done = false;
      bool initial_cycle = true;
      for (int ln = 0; hierarchy->levelCanBeRefined(ln) && !done;
           ++ln) {
         gridding_algorithm->makeFinerLevel(
            tag_buffer_array[ln],
            initial_cycle,
            init_cycle,
            init_time);
         done = !(hierarchy->finerLevelExists(ln));
      }

      /*
       * Setup timer manager for profiling code.
       */
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      std::shared_ptr<tbox::Timer> t_cvode_solve(
         tbox::TimerManager::getManager()->
         getTimer("apps::main::cvode_solver"));
      std::shared_ptr<tbox::Timer> t_log_dump(
         tbox::TimerManager::getManager()->
         getTimer("apps::main::Solution log dump"));
      /*
       * Setup solution vector.
       */
      cvode_model->setupSolutionVector(hierarchy);
      SundialsAbstractVector* solution_vector =
         cvode_model->getSolutionVector();

      /*
       * Set initial conditions vector.
       */
      cvode_model->setInitialConditions(solution_vector);

      /**************************************************************************
      * Setup CVODESolver object.
      **************************************************************************/
      solv::CVODESolver* cvode_solver =
         new solv::CVODESolver("cvode_solver",
            cvode_model.get(),
            uses_preconditioning);

      size_t neq = 0;
      std::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      const hier::BoxContainer& level_0_boxes = level_zero->getBoxes();
      for (hier::BoxContainer::const_iterator i = level_0_boxes.begin();
           i != level_0_boxes.end(); ++i) {
         neq += i->size();
      }
      cvode_solver->setRelativeTolerance(relative_tolerance);
      cvode_solver->setAbsoluteTolerance(absolute_tolerance);
      cvode_solver->setMaximumNumberOfInternalSteps(max_internal_steps);
      cvode_solver->setSteppingMethod(stepping_method);
      cvode_solver->setMaximumLinearMultistepMethodOrder(max_order);
      if (uses_preconditioning) {
         cvode_solver->setPreconditioningType(PREC_LEFT);
      }

      cvode_solver->setInitialValueOfIndependentVariable(init_time);
      cvode_solver->setInitialConditionVector(solution_vector);
      cvode_solver->initialize(solution_vector);

      /*
       * Print initial vector (if solution logging is enabled)
       */
      std::shared_ptr<solv::SAMRAIVectorReal<double> > y_init(
         solv::Sundials_SAMRAIVector::getSAMRAIVector(solution_vector));

      if (solution_logging) {

         std::shared_ptr<hier::PatchHierarchy> init_hierarchy(
            y_init->getPatchHierarchy());

         tbox::pout << "Initial solution vector y() at initial time: " << std::endl;
         int ln;
         tbox::pout << "y(" << init_time << "): " << std::endl;
         for (ln = 0; ln < init_hierarchy->getNumberOfLevels(); ++ln) {
            std::shared_ptr<hier::PatchLevel> level(
               init_hierarchy->getPatchLevel(ln));
            tbox::plog << "level = " << ln << std::endl;

            for (hier::PatchLevel::iterator p(level->begin());
                 p != level->end(); ++p) {
               const std::shared_ptr<hier::Patch>& patch = *p;

               std::shared_ptr<CellData<double> > y_data(
                  SAMRAI_SHARED_PTR_CAST<CellData<double>, hier::PatchData>(
                     y_init->getComponentPatchData(0, *patch)));
               TBOX_ASSERT(y_data);
               y_data->print(y_data->getBox());
            }
         }
      }

      /*
       * Compute maxNorm and L1Norm of initial vector
       */
      if (solution_logging) {
         tbox::pout << "\n\nBefore solve..." << std::endl;
         tbox::pout << "Max Norm of y()= " << y_init->maxNorm() << std::endl;
         tbox::pout << "L1 Norm of y()= " << y_init->L1Norm() << std::endl;
         tbox::pout << "L2 Norm of y()= " << y_init->L2Norm() << std::endl;
      }

      /**************************************************************************
      * Start time-stepping.
      **************************************************************************/

      std::vector<double> time(num_print_intervals);
      std::vector<double> maxnorm(num_print_intervals);
      std::vector<double> l1norm(num_print_intervals);
      std::vector<double> l2norm(num_print_intervals);

      double final_time = init_time;
      int interval;
      for (interval = 1; interval <= num_print_intervals; ++interval) {

         /*
          * Set time interval
          */
         final_time += print_interval;
         cvode_solver->setFinalValueOfIndependentVariable(final_time, false);

         /*
          * Perform CVODE solve to the requested interval time.
          */
         t_cvode_solve->start();
         tbox::plog << "return code = " << cvode_solver->solve() << std::endl;
         t_cvode_solve->stop();
         double actual_time =
            cvode_solver->getActualFinalValueOfIndependentVariable();

         /*
          * Print statistics
          * Format:  time  max norm   l1 norm   l2 norm
          */
         std::shared_ptr<solv::SAMRAIVectorReal<double> > y_result(
            solv::Sundials_SAMRAIVector::getSAMRAIVector(solution_vector));
         std::shared_ptr<hier::PatchHierarchy> result_hierarchy(
            y_result->getPatchHierarchy());

         time[interval - 1] = actual_time;
         maxnorm[interval - 1] = y_result->maxNorm();
         l1norm[interval - 1] = y_result->L1Norm();
         l2norm[interval - 1] = y_result->L2Norm();

         if (solution_logging) {
            cvode_solver->printStatistics(tbox::pout);
         }

         /*
          * Write solution (if desired).
          */
         if (solution_logging) {
            tbox::plog << "y(" << final_time << "): " << std::endl << std::endl;
            t_log_dump->start();
            for (int ln = 0; ln < result_hierarchy->getNumberOfLevels();
                 ++ln) {
               std::shared_ptr<hier::PatchLevel> level(
                  result_hierarchy->getPatchLevel(ln));
               tbox::plog << "level = " << ln << std::endl;

               for (hier::PatchLevel::iterator p(level->begin());
                    p != level->end(); ++p) {
                  const std::shared_ptr<hier::Patch>& patch = *p;

                  std::shared_ptr<CellData<double> > y_data(
                     SAMRAI_SHARED_PTR_CAST<CellData<double>, hier::PatchData>(
                        y_result->getComponentPatchData(0, *patch)));
                  TBOX_ASSERT(y_data);
                  y_data->print(y_data->getBox());
               }
            }
            t_log_dump->stop();
         }
      } // end of timestep loop

      /*************************************************************************
       * Write summary information
       ************************************************************************/
      /*
       * Write CVODEModel stats
       */
      std::vector<int> counters;
      cvode_model->getCounters(counters);

#if (TESTING == 1)
      int correct_rhs_evals = main_db->getInteger("correct_rhs_evals");
      int correct_precond_setups = main_db->getInteger("correct_precond_setups");
      int correct_precond_solves = main_db->getInteger("correct_precond_solves");

      if (counters[0] == correct_rhs_evals) {
         tbox::plog << "Test 0: Number RHS evals CORRECT" << std::endl;
      } else {
         tbox::perr << "Test 0 FAILED: Number RHS evals INCORRECT" << std::endl;
         tbox::perr << "Correct Number RHS evals:  " << correct_rhs_evals
                    << std::endl;
         tbox::perr << "Number RHS evals computed: " << counters[0] << std::endl;
      }

      if (counters[1] == correct_precond_setups) {
         tbox::plog << "Test 1: Number precond setups CORRECT" << std::endl;
      } else {
         tbox::perr << "Test 1 FAILED: Number precond setups INCORRECT" << std::endl;
         tbox::perr << "Correct number precond setups:  "
                    << correct_precond_setups << std::endl;
         tbox::perr << "Number precond setups computed: " << counters[1]
                    << std::endl;
      }

      if (counters[2] == correct_precond_solves) {
         tbox::plog << "Test 2: Number precond solves CORRECT" << std::endl;
      } else {
         tbox::perr << "Test 2 FAILED: Number precond solves INCORRECT" << std::endl;
         tbox::perr << "Correct number precond solves:  "
                    << correct_precond_solves << std::endl;
         tbox::perr << "Number precond solves computed: " << counters[2]
                    << std::endl;
      }
#endif

      if (solution_logging) {
         tbox::plog << "\n\nEnd Timesteps - final time = " << final_time
                    << "\n\tTotal number of RHS evaluations = " << counters[0]
                    << "\n\tTotal number of precond setups = " << counters[1]
                    << "\n\tTotal number of precond solves = " << counters[2]
                    << std::endl;

         /*
          * Write out timestep sequence information
          */
         tbox::pout << "\n\nTimestep Summary of solution vector y()\n"
                    << "  time                   \t"
                    << "  Max Norm  \t"
                    << "  L1 Norm  \t"
                    << "  L2 Norm  " << std::endl;

         for (interval = 0; interval < num_print_intervals; ++interval) {
            tbox::pout.precision(18);
            tbox::pout << "  " << time[interval] << "  \t";
            tbox::pout.precision(6);
            tbox::pout << "  " << maxnorm[interval] << "  \t"
                       << "  " << l1norm[interval] << "  \t"
                       << "  " << l2norm[interval] << std::endl;
         }
      }

      /*
       * Write out timings
       */
#if (TESTING != 1)
      tbox::TimerManager::getManager()->print(tbox::pout);
#endif

      /*
       * Memory cleanup.
       */
      if (cvode_solver) delete cvode_solver;

      cvode_model.reset();
      gridding_algorithm.reset();
      error_est.reset();
      load_balancer.reset();
      box_generator.reset();
      hierarchy.reset();
      geometry.reset();

#endif // HAVE_SUNDIALS

      tbox::pout << "\nPASSED:  cvode" << std::endl;

   }

   /*
    * Shutdown SAMRAI and tbox::MPI.
    */
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
