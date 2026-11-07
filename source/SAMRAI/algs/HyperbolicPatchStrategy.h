/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic integration scheme.
 *
 ************************************************************************/

#ifndef included_algs_HyperbolicPatchStrategy
#define included_algs_HyperbolicPatchStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"

#include <memory>

namespace SAMRAI {
namespace algs {

class HyperbolicLevelIntegrator;

/**
 * Class HyperbolicPatchStrategy is an abstract base class defining the
 * interface between an HyperbolicLevelIntegrator object and operations
 * applied to a single patch in a structured AMR hierarchy.   The operations
 * include patch initialization, dt calculation, flux computation,
 * conservative differencing, and error estimation.  This class is derived
 * from the xfer::RefinePatchStrategy and xfer::CoarsenPatchStrategy
 * abstract base classes.   These base classes provide the interface for
 * user-defined interlevel data refining and coarsening operations and the
 * specification of physical boundary conditions.  The functions
 * setPhysicalBoundaryConditions(), and pre/postprocessRefine() are
 * overloaded from the class xfer::RefinePatchStrategy.  The operations
 * pre/postprocessCoarsen() are overloaded from xfer::CoarsenPatchStrategy.
 * The pre/postprocessCoarsen/Refine() operations are given empty
 * implementations here so that the user does not need to proovide them if
 * the operations are not needed.
 *
 * It is important to recognize that for the concrete patch strategy subclass
 * and the HyperbolicLevelIntegrator to work together, the concrete
 * strategy must know which patch data to operate on.  The patch data storage
 * is manipulated by the level integrator.  The set/clearDataContext() methods
 * allow the integrator to inform the patch strategy of the correct data
 * context.  The concrete patch strategy subclass can access the appropriate
 * context via the getDataContext() method.
 *
 * @see HyperbolicLevelIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

class HyperbolicPatchStrategy:
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy
{
public:
   /**
    * Blank constructor for HyperbolicPatchStrategy.
    */
   HyperbolicPatchStrategy();

   /**
    * Virtual destructor for HyperbolicPatchStrategy.
    */
   virtual ~HyperbolicPatchStrategy();

   /**
    * Register specific variables needed in the numerical routines with the
    * hyperbolic level integrator using the registerVariable() function in that
    * class.  The integrator manipulates storage for the data and this
    * registration defines the way in which data for each quantity will
    * be manipulated on the patches.  Typically, the derived data quantities
    * for plotting are registered with a visualization data writer in this
    * routine as well, since the hyperbolic level integrator provides the
    * variable context for plotting (i.e., which data is available when a
    * plot file is generated).  The integrator pointer cannot be null in
    * most cases.
    *
    * The gridding algorithm pointer is provided so that patch data objects
    * may be registered with the load balancer object (owned by the gridding
    * algorithm) for non-uniform load balancing, if needed.
    */
   virtual void
   registerModelVariables(
      HyperbolicLevelIntegrator* integrator) = 0;

   /**
    * Set up parameters in the load balancer object (owned by the gridding
    * algorithm) if needed.  This function is called immediately after the
    * registerModelVariables() function is called.  The hyperbolic level
    * integrator pointer is provided so that the integrator can be used
    * to manage data for the load balancer if needed (e.g., when using
    * non-uniform load balancing).
    *
    * Note that this function is not pure virtual. It is given a
    * dummy implementation here so that users may ignore it when
    * inheriting from this class.
    */
   virtual void
   setupLoadBalancer(
      HyperbolicLevelIntegrator* integrator,
      mesh::GriddingAlgorithm* gridding_algorithm);

   /**
    * Set the initial data on a patch interior only.  Note that no ghost cells
    * need to be set in this routine regardless of whether the patch data
    * corresponding to the context requires ghost cells.  The data_time
    * is the simulation time when the routine is called.  The boolean
    * initial_time is true if the routine is called at the initial time
    * when the hierarchy is initially constructed, otherwise it is false.
    */
   virtual void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double data_time,
      const bool initial_time) = 0;

   /**
    * Compute the stable time increment for a patch on the level with the
    * given number.  The boolean flag initial_time is true if the routine
    * is called at the initial simulation time; otherwise it is false.
    * The double argument dt_time is the simulation time.
    */
   virtual double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const bool initial_time,
      const double dt_time) = 0;

   /**
    * Compute TIME INTEGRALS of fluxes to be used in conservative difference
    * for patch integration.  That is, it is assumed that this numerical
    * routine will compute the fluxes corresponding to the cell faces
    * multiplied by the time increment.  Typically, the numerical flux is
    * the normal flux at the cell face.  The flux integrals will be used in
    * the conservative difference that updates the conserved quantities.
    *
    * Note that the numerical routines in this method generally require
    * ghost cells.  Ghost cells data is filled before this routine is called.
    */
   virtual void
   computeFluxesOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt) = 0;

   /**
    * Update patch data with a conservative difference (approximating
    * the divergence theorem) using the flux integrals computed in
    * computeFluxesOnPatch() routine.  The boolean flag is true when this
    * routine is called during a flux synchronization step.   Otherwise,
    * it is false.  Note that the computeFluxesOnPatch() routine computes
    * TIME INTEGRALs of the numerical fluxes (e.g., they have been multiplied
    * by dt).  So the conservative difference routine should be consistent
    * with this.
    */
   virtual void
   conservativeDifferenceOnPatch(
      hier::Patch& patch,
      const double time,
      const double dt,
      bool at_syncronization) = 0;

   /**
    * This is an optional routine for user to process any application-specific
    * patch strategy data BEFORE patches are advanced on the given level.
    * This routine is called after patch boundary data is filled
    * (i.e., ghosts) and before computeFluxesOnPatch().  The arguments are:
    * level -- level that will be advanced, current_time -- current
    * integration time, dt -- current time increment, first_step -- boolean
    * flag that is true if advance is first in time step sequence on level
    * (i.e., previous advance step was on another level, false otherwise,
    * last_step -- boolean flag that is true if advance is last in time step
    * sequence on level (i.e., synchronization with coarser level will occur
    * immediately after this advance), regrid_advance -- boolean flag that
    * is true if the advance is during a regridding phase (i.e., the advance
    * is not used to integrate data on the hierarchy) in which case the
    * results of the advance will be discarded.
    *
    * Note that when this routine is called, the scratch data is filled on
    * all patches (i.e., ghost cells) and that data is the same as the
    * current level data on all patch interiors.  That is, both scratch and
    * current data correspond to current_time.
    *
    * Note that this function is not pure virtual. It is given a
    * dummy implementation here so that users may ignore it when
    * inheriting from this class.
    */
   virtual void
   preprocessAdvanceLevelState(
      const std::shared_ptr<hier::PatchLevel>& level,
      double current_time,
      double dt,
      bool first_step,
      bool last_step,
      bool regrid_advance);

   /**
    * This is an optional routine for user to process any application-specific
    * patch strategy data AFTER patches are advanced on the given level.
    * This routine is called after conservativeDifferenceOnPatch() is called
    * and before computeStableDtOnPatch().  The arguments are:
    * level -- level that will be advanced, current_time -- current
    * integration time, dt -- current time increment, first_step -- boolean
    * flag that is true if advance is first in time step sequence on level
    * (i.e., previous advance step was on another level, false otherwise,
    * last_step -- boolean flag that is true if advance is last in time step
    * sequence on level (i.e., synchronization with coarser
    * level will occur immediately after this advance), regrid_advance --
    * boolean flag that is true if the advance is during a regridding phase
    * (i.e., the advance is not used to integrate data on the hierarchy) in
    * which case the results of the advance will be discarded.
    *
    * Note that when this routine is called, the scratch data is filled on
    * all patches (i.e., ghost cells) and that data is the same as the
    * new level data on all patch interiors.  That is, both scratch and
    * new data correspond to current_time + dt on patch interiors.
    * The current data and ghost values correspond to the current_time.
    *
    * Note that this function is not pure virtual. It is given a
    * dummy implementation here so that users may ignore it when
    * inheriting from this class.
    */
   virtual void
   postprocessAdvanceLevelState(
      const std::shared_ptr<hier::PatchLevel>& level,
      double current_time,
      double dt,
      bool first_step,
      bool last_step,
      bool regrid_advance);

   /**
    * Tag cells on the given patch that require refinement based on
    * application-specific numerical quantities.  The tag index argument
    * indicates the index of the tag data on the patch data array.  The
    * boolean argument initial_error is true if tagging is being done at the
    * initial simulation time; otherwise, it is false.  The other boolean
    * flag uses_richardson_extrapolation_too is true when Richardson
    * extrapolation is used in addition to the gradient detector.  This flag
    * helps users manage multiple regridding criteria.
    *
    * Note that this function is not pure virtual. It is given a
    * dummy implementation here so that users may ignore it when
    * inheriting from this class.
    */
   virtual void
   tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_index,
      const bool uses_richardson_extrapolation_too);

   /**
    * Tag cells based from differences computed in the Richardson
    * extrapolation.  The Richardson
    * extrapolation algorithm creates a coarsened version of some hierarchy
    * patch level and advances data in time on both the coarsened patch
    * level and the hierarchy level. This routine takes the data resulting
    * from the advance on both the coarse and fine levels, compares them, and
    * tags cells according to the difference.
    * \verbatim
    *              (2)
    *      n+1 ^------->x finish      (1) advanced_coarse
    *          |        ^             (2) coarsened_fine
    *  time  n -        |
    *          ^        |(1)
    *          |        |
    *          <--------o start
    *        fine     coarse
    * \endverbatim
    *
    * The patch supplied to this routine is on the coarsened level.  However,
    * the error_level_number corresponds to the actual hierarchy level
    * from which it was coarsened.   Data resides on this patch in two
    * contexts - ``advanced_coarse'' and ``coarsened_fine''.  Advanced
    * coarse is data advanced on the coarsened version of the level, while
    * coarsened fine is the data advanced on the fine level and then
    * coarsened to the coarse level.  The regrid time and the time increment
    * are given for the actual hierarchy level.  The error coarsen ratio
    * argument is the ratio between the index spaces on the hierarchy level
    * and the coarsened hierarchy level.  The boolean flag ``initial_error''
    * is true when the error estimation is performed at the initial simulation
    * time; i.e., when the hierarchy levels are being constructed for the first
    * time.  The tag index argument is the index of the tag data on the patch
    * data array.  The other boolean flag uses_gradient_detector_too is
    * true when a gradient detector scheme is used in addition to Richardson
    * extrapolation.  This flag helps users manage multiple regridding
    * criteria.
    *
    * Note that this function is not pure virtual. It is given a
    * dummy implementation here so that users may ignore it when
    * inheriting from this class.
    */
   virtual void
   tagRichardsonExtrapolationCells(
      hier::Patch& patch,
      const int error_level_number,
      const std::shared_ptr<hier::VariableContext>& coarsened_fine,
      const std::shared_ptr<hier::VariableContext>& advanced_coarse,
      const double regrid_time,
      const double deltat,
      const int error_coarsen_ratio,
      const bool initial_error,
      const int tag_index,
      const bool uses_gradient_detector_too);

   /**
    * Set user-defined boundary conditions at the physical domain boundary.
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill) = 0;

   /*!
    * @brief Check the tags on a tagged patch.
    *
    * This virtual method provides an interface for a callback into
    * application code to check the values held in user tag PatchData
    * The tag data will contain the tags created by application code
    * plus any tags added internally by the GriddingAlgorithm (for
    * example, buffering).
    *
    * This is a virtual method with a no-op implementation provided, so that
    * users who do not need this callback need not implement anything.
    *
    * @param[in] patch
    * @param[in] tag_index     Patch data index for user tags
    */
   virtual void
   checkUserTagData(
      hier::Patch& patch,
      const int tag_index) const
   {
      NULL_USE(patch);
      NULL_USE(tag_index);
   }


   /*!
    * @brief Check the tags on a newly-created patch.
    *
    * This virtual method provides an interface for a callback into
    * application code to check the tag values that have been saved on
    * a new patch that has been created during intitialization or
    * regridding.  The tag values will be the values of the user tags
    * on the coarser level, constant refined onto the cells of the patch
    * on the new level.
    *
    * This is a virtual method with a no-op implementation provided, so that
    * users who do not need this callback need not implement anything.
    *
    * @param[in] patch
    * @param[in] tag_index      Patch data index for the new tags.
    */
   virtual void
   checkNewPatchTagData(
      hier::Patch& patch,
      const int tag_index) const
   {
      NULL_USE(patch);
      NULL_USE(tag_index);
   }

   /**
    * Return pointer to patch data context.
    */
   std::shared_ptr<hier::VariableContext>
   getDataContext() const
   {
      return d_data_context;
   }

   /**
    * The hyperbolic integrator controls the context for the data to be used
    * in the numerical routines implemented in the concrete patch strategy.
    * The setDataContext() allows the integrator to set the context for
    * data on a patch on which to operate.
    */
   void
   setDataContext(
      const std::shared_ptr<hier::VariableContext>& context)
   {
      d_data_context = context;
   }

   /**
    * The clearDataContext() routine resets the data context to be null.
    */
   void
   clearDataContext()
   {
      d_data_context.reset();
   }

private:
   std::shared_ptr<hier::VariableContext> d_data_context;
};

}
}

#endif
