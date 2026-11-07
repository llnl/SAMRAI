/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implicit time integration manager class for nonlinear problems.
 *
 ************************************************************************/
#include "SAMRAI/algs/ImplicitIntegrator.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace algs {

const int ImplicitIntegrator::ALGS_IMPLICIT_INTEGRATOR_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for ImplicitIntegrator.  The
 * constructor sets default values for data members, then overrides
 * them with values read from input or restart.  The destructor does
 * nothing interesting.
 *
 *************************************************************************
 */

ImplicitIntegrator::ImplicitIntegrator(
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db,
   ImplicitEquationStrategy* implicit_equations,
   solv::NonlinearSolverStrategy* nonlinear_solver,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy):
   d_object_name(object_name),
   d_implicit_equations(implicit_equations),
   d_nonlinear_solver(nonlinear_solver),
   d_patch_hierarchy(hierarchy),
   d_finest_level(-1),
   d_initial_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_final_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_current_time(tbox::MathUtilities<double>::getSignalingNaN()),
   d_current_dt(tbox::MathUtilities<double>::getSignalingNaN()),
   d_old_dt(tbox::MathUtilities<double>::getSignalingNaN()),
   d_integrator_step(0),
   d_max_integrator_steps(0)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(implicit_equations != 0);
   TBOX_ASSERT(nonlinear_solver != 0);
   TBOX_ASSERT(hierarchy);

   /*
    * Initialize object with data read from input and restart databases.
    */

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   d_current_dt = d_implicit_equations->getInitialDt();
   d_old_dt = 0.0;
   d_current_time = d_initial_time;

}

ImplicitIntegrator::~ImplicitIntegrator()
{
}

/*
 *************************************************************************
 *
 * Initialize integrator and nonlinear solver:
 *
 * (1) Create vector containing solution state advanced in time.
 *
 * (2) Equation class registers data components with solution vector.
 *
 * (3) Initialize nonlinear solver.
 *
 *************************************************************************
 */

void
ImplicitIntegrator::initialize()
{
   d_finest_level = d_patch_hierarchy->getFinestLevelNumber();

   d_solution_vector.reset(
      new solv::SAMRAIVectorReal<double>("solution_vector",
         d_patch_hierarchy,
         0, d_finest_level));

   d_implicit_equations->setupSolutionVector(d_solution_vector);

   if (d_solution_vector->getNumberOfComponents() == 0) {
      TBOX_ERROR("Solution vector has zero components." << std::endl);
   }

   d_nonlinear_solver->initialize(d_solution_vector);
}

/*
 *************************************************************************
 *
 * Integrate solution through given time increment:
 *
 * (1) If number of levels in hierarchy has changed since last advance
 *     due to regridding, the range of levels in the vectors is reset.
 *
 * (2) Construct initial guess at new solution by extrapolation.
 *
 * (3) Call the equation advance set up routine.
 *
 * (4) Compute the new solution using the nonlinear solver.
 *
 * (5) Return integer return code define by nonlinear solver.
 *
 *************************************************************************
 */

int
ImplicitIntegrator::advanceSolution(
   const double dt,
   const bool first_step)
{
   int retcode = tbox::MathUtilities<int>::getMax();

   if (stepsRemaining() && (d_current_time < d_final_time)) {

      d_current_dt = dt;

      const int finest_now = d_patch_hierarchy->getFinestLevelNumber();

      if (first_step && (finest_now != d_finest_level)) {

         d_finest_level = finest_now;

         d_solution_vector->resetLevels(0, d_finest_level);

         d_nonlinear_solver->initialize(d_solution_vector);

      }

      d_implicit_equations->setInitialGuess(first_step,
         d_current_time,
         d_current_dt,
         d_old_dt);

      retcode = d_nonlinear_solver->solve();

   }

   return retcode;
}

/*
 *************************************************************************
 *
 * Get next dt from user-supplied equation class.  Timestep selection
 * is generally based on whether the nonlinear solution iteration
 * converged and, if so, whether the solution meets some user-defined
 * criteria.  It is assumed that, before this routine is called, the
 * routine checkNewSolution() has been called.  The boolean argument
 * is the return value from that call.  The integer argument is
 * that which is returned by the particular nonlinear solver package
 * that generated the solution.
 *
 *************************************************************************
 */

double
ImplicitIntegrator::getNextDt(
   const bool good_solution,
   const int solver_retcode)
{
   double dt_next = d_implicit_equations->getNextDt(good_solution,
         solver_retcode);

   double global_dt_next = dt_next;
   const tbox::SAMRAI_MPI& mpi(d_patch_hierarchy->getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_dt_next, 1, MPI_MIN);
   }

   global_dt_next =
      tbox::MathUtilities<double>::Min(global_dt_next,
         d_final_time - d_current_time);

   return global_dt_next;
}

/*
 *************************************************************************
 *
 * Check whether time advanced solution is acceptable.  Note that the
 * user-supplied solution checking routine must interpret the integer
 * return code generated by the nonlinear solver correctly.
 *
 *************************************************************************
 */

bool
ImplicitIntegrator::checkNewSolution(
   const int solver_retcode) const
{
   bool good_solution =
      d_implicit_equations->checkNewSolution(solver_retcode);

   int good = (good_solution ? 1 : 0);
   int global_good = good;
   const tbox::SAMRAI_MPI& mpi(d_patch_hierarchy->getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_good, 1, MPI_MIN);
   }

   return global_good == 0 ? false : true;
}

/*
 *************************************************************************
 *
 * Assuming an acceptable time advanced solution is found, update
 * solution quantities and time information state of integrator.
 * Return the current simulation time.
 *
 *************************************************************************
 */

double
ImplicitIntegrator::updateSolution()
{
   d_current_time += d_current_dt;
   d_old_dt = d_current_dt;
   ++d_integrator_step;

   d_implicit_equations->updateSolution(d_current_time);

   return d_current_time;
}

/*
 *************************************************************************
 *
 * If simulation is not from restart, read data from input database.
 * Otherwise, override restart values for a subset of the data members
 * with those found in input.
 *
 *************************************************************************
 */

void
ImplicitIntegrator::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (!is_from_restart && !input_db) {
      TBOX_ERROR(": ImplicitIntegrator::getFromInput()\n"
         << "no input database supplied" << std::endl);
   }

   if (!is_from_restart) {

      d_initial_time = input_db->getDouble("initial_time");
      if (!(d_initial_time >= 0)) {
         INPUT_RANGE_ERROR("initial_time");
      }

      d_final_time = input_db->getDouble("final_time");
      if (!(d_final_time >= d_initial_time)) {
         INPUT_RANGE_ERROR("final_time");
      }

      d_max_integrator_steps = input_db->getInteger("max_integrator_steps");
      if (!(d_max_integrator_steps >= 0)) {
         INPUT_RANGE_ERROR("max_integrator_steps");
      }

   } else if (input_db) {
      bool read_on_restart =
         input_db->getBoolWithDefault("read_on_restart", false);

      if (read_on_restart) {
         if (input_db->keyExists("initial_time")) {
            double tmp = input_db->getDouble("initial_time");
            if (tmp != d_initial_time) {
               TBOX_WARNING("ImplicitIntegrator::getFromInput warning...\n"
                  << "initial_time may not be changed on restart."
                  << std::endl);
            }
         }

         d_final_time =
            input_db->getDoubleWithDefault("final_time", d_final_time);
         if (d_final_time < d_initial_time) {
            TBOX_ERROR("ImplicitIntegrator::getFromInput() error...\n"
               << "final_time must be >= initial_time " << std::endl);
         }

         d_max_integrator_steps =
            input_db->getIntegerWithDefault("max_integrator_steps",
               d_max_integrator_steps);
         if (d_max_integrator_steps < 0) {
            TBOX_ERROR("ImplicitIntegrator::getFromInput() error...\n"
               << "max_integrator_steps must be >= 0." << std::endl);
         } else if (d_max_integrator_steps < d_integrator_step) {
            TBOX_ERROR("ImplicitIntegrator::getFromInput() error...\n"
               << "max_integrator_steps must be >= current integrator step."
               << std::endl);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Write out class version number and data members to restart database.
 *
 *************************************************************************
 */

void
ImplicitIntegrator::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("ALGS_IMPLICIT_INTEGRATOR_VERSION",
      ALGS_IMPLICIT_INTEGRATOR_VERSION);

   restart_db->putDouble("initial_time", d_initial_time);
   restart_db->putDouble("final_time", d_final_time);

   restart_db->putInteger("d_integrator_step", d_integrator_step);
   restart_db->putInteger("max_integrator_steps", d_max_integrator_steps);

}

/*
 *************************************************************************
 *
 * Check to make sure that the version number of the class is that same
 * as the version number in the restart file.  If these values are equal
 * then read values for data members from the restart file.
 *
 *************************************************************************
 */

void
ImplicitIntegrator::getFromRestart()
{

   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file" << std::endl);
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("ALGS_IMPLICIT_INTEGRATOR_VERSION");
   if (ver != ALGS_IMPLICIT_INTEGRATOR_VERSION) {
      TBOX_ERROR(d_object_name << ":  "
                               << "Restart file version different "
                               << "than class version." << std::endl);
   }

   d_initial_time = db->getDouble("initial_time");
   d_final_time = db->getDouble("final_time");

   d_integrator_step = db->getInteger("d_integrator_step");
   d_max_integrator_steps = db->getInteger("max_integrator_steps");

}

/*
 *************************************************************************
 *
 * Print class data members to given output stream.
 *
 *************************************************************************
 */

void
ImplicitIntegrator::printClassData(
   std::ostream& os) const
{
   os << "\nImplicitIntegrator::printClassData..." << std::endl;
   os << "ImplicitIntegrator: this = "
      << (ImplicitIntegrator *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_implicit_equations = "
      << (ImplicitEquationStrategy *)d_implicit_equations << std::endl;
   os << "d_nonlinear_solver = "
      << (solv::NonlinearSolverStrategy *)d_nonlinear_solver << std::endl;
   os << "d_patch_hierarchy = "
      << d_patch_hierarchy.get() << std::endl;
   os << "d_solution_vector = "
      << d_solution_vector.get() << std::endl;

   os << "d_finest_level = " << d_finest_level << std::endl;
   os << "d_initial_time = " << d_initial_time << std::endl;
   os << "d_final_time = " << d_final_time << std::endl;
   os << "d_current_time = " << d_current_time << std::endl;
   os << "d_current_dt = " << d_current_dt << std::endl;
   os << "d_old_dt = " << d_old_dt << std::endl;
   os << "d_integrator_step = " << d_integrator_step << std::endl;
   os << "d_max_integrator_steps = " << d_max_integrator_steps << std::endl;
}

}
}
