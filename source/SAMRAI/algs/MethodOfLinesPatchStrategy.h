/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to application-specific patch functions in support
 *                Method of Lines integration algorithm
 *
 ************************************************************************/

#ifndef included_algs_MethodOfLinesPatchStrategy
#define included_algs_MethodOfLinesPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <memory>

namespace SAMRAI {
namespace algs {

class MethodOfLinesIntegrator;

/**
 * Class MethodOfLinesPatchStrategy is an abstract type defining the
 * interface for operations invoked during the integration routines defined
 * in the MethodOfLinesIntegrator class.  This class is derived from
 * the xfer::RefinePatchStrategy and xfer::CoarsenPatchStrategy abstract
 * base classes.  These base classes define the interfaces for user-defined
 * interlevel data refining and coarsening operations and the specification
 * of physical boundary conditions.
 *
 * @see MethodOfLinesIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

class MethodOfLinesPatchStrategy:
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy
{
public:
   /*!
    * Blank constructor for MethodOfLinesPatchStrategy.
    */
   MethodOfLinesPatchStrategy();

   /*!
    * Virtual destructor for MethodOfLinesPatchStrategy.
    */
   virtual ~MethodOfLinesPatchStrategy() = 0;

   /*!
    * Register variables specific to the problem to be solved with the
    * integrator using the registerVariable function.  This
    * defines the way data for each quantity will be manipulated on the
    * patches.  For more information, refer to
    * MethodOfLinesIntegrator::registerVariable.
    */
   virtual void
   registerModelVariables(
      MethodOfLinesIntegrator* integrator) = 0;

   /*!
    * Set the initial data on a patch interior (i.e., NO GHOST CELLS).
    * Setting "initial_time" true will initialize data at time. Setting it
    * false will interpolate data from the appropriate coarser patch.
    */
   virtual void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double time,
      const bool initial_time) const = 0;

   /*!
    * Compute the stable time increment for a patch.
    */
   virtual double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const double time) const = 0;

   /*!
    * Advance a single Runge Kutta step.
    *
    * @param patch patch that RK step is being applied
    * @param dt    timestep
    * @param alpha_1 first coefficient applied in the RK step
    * @param alpha_2 second coefficient
    * @param beta    third coefficient
    */
   virtual void
   singleStep(
      hier::Patch& patch,
      const double dt,
      const double alpha_1,
      const double alpha_2,
      const double beta) const = 0;

   /*!
    * Using a user-specified gradient detection scheme, determine cells which
    * have high gradients and, consequently, should be refined.
    */
   virtual void
   tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_index,
      const bool uses_richardson_extrapolation_too) = 0;

   /*!
    * Set user-defined boundary conditions at the physical domain boundary.
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill) = 0;

   /*!
    * The method of lines integrator controls the context for the data to
    * be used in the numerical routines implemented in the concrete patch
    * strategy. These accessor methods allow the patch strategy to access
    * the particular data contexts used in the integrator.
    *
    * Return pointer to data context with ghost cells.
    */
   std::shared_ptr<hier::VariableContext>
   getInteriorWithGhostsContext() const
   {
      return d_interior_with_ghosts;
   }

   /*!
    * Return pointer to data context with NO ghosts.
    */
   std::shared_ptr<hier::VariableContext>
   getInteriorContext() const
   {
      return d_interior;
   }

   /*!
    * Set pointer to data context with ghosts.
    */
   void
   setInteriorWithGhostsContext(
      const std::shared_ptr<hier::VariableContext>& context)
   {
      d_interior_with_ghosts = context;
   }

   /*!
    * Set pointer to data context with NO ghosts.
    */
   void
   setInteriorContext(
      const std::shared_ptr<hier::VariableContext>& context)
   {
      d_interior = context;
   }

private:
   std::shared_ptr<hier::VariableContext> d_interior_with_ghosts;
   std::shared_ptr<hier::VariableContext> d_interior;
};

}
}

#endif
