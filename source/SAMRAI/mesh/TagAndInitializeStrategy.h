/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for params, tagging, init for gridding.
 *
 ************************************************************************/

#ifndef included_mesh_TagAndInitializeStrategy
#define included_mesh_TagAndInitializeStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/BoxLevel.h"

#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * Class TagAndInitializeStrategy is a base class that defines a
 * Strategy pattern interface for level initialization and cell tagging
 * routines that are needed by the adaptive meshing algorithms provided
 * by the class GriddingAlgorithm.  The class defines an interface to
 * construct refined regions based on a user-supplied set of boxes, but
 * its main role is to provide interfaces for level initialization and
 * cell tagging operations.
 *
 * All methods defined by this class with the exception of getObjectName are
 * abstract and must be supplied by a concrete sub-class of this base class.
 *
 * If user supplied refine boxes are used, they may be supplied through
 * input as defined by specific, concrete sub-classes of this base class.
 * Alternatively, they may be supplied through the "resetRefineBoxes()"
 * method.
 *
 * The virtual methods in this class may place constraints on the patch
 * hierarchy by the particular error estimation procedure in use.  Those
 * constraints and operations must be honored in the concrete subclass
 * implementations of these methods.  The constraints are discussed in
 * the method descriptions below.
 *
 * @see GriddingAlgorithm
 */

class TagAndInitializeStrategy
{
public:
   /*!
    * Empty constructor for TagAndInitializeStrategy.
    */
   TagAndInitializeStrategy(
      const std::string& object_name);

   /*!
    * Empty destructor for TagAndInitializeStrategy.
    */
   virtual ~TagAndInitializeStrategy();

   /*!
    * Return user supplied set of refine boxes for specified level number
    * and time.  The boolean return value specifies whether the boxes
    * have been reset from the last time this method was called.  If they
    * have been reset, it returns true.  If they are unchanged, it returns
    * false.
    */
   virtual bool
   getUserSuppliedRefineBoxes(
      hier::BoxContainer& refine_boxes,
      const int level_number,
      const int cycle,
      const double time) = 0;

   /*!
    * Reset the static refine boxes for the specified level number in the
    * hierarchy.  The level number must be greater than or equal to zero.
    */
   virtual void
   resetRefineBoxes(
      const hier::BoxContainer& refine_boxes,
      const int level_number) = 0;

   /*!
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.  The old_level pointer corresponds to
    * the level that resided in the hierarchy before the level with the
    * specified number was introduced.  If this pointer is null, there was
    * no level in the hierarchy prior to the call and the data on the new
    * level is set by interpolating data from coarser levels in the hierarchy.
    * Otherwise, the the new level is initialized by interpolating data from
    * coarser levels and copying data from the old level before it is
    * destroyed.
    *
    * The boolean argument initial_time indicates whether the integration
    * time corresponds to the initial simulation time.  If true, the level
    * should be initialized with initial simulation values.  Otherwise, it
    * should be assumed that the simulation time is at some point after the
    * start of the simulation.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  In any case, the double "time" value is the current
    * simulation time for the level.  The can_be_refined boolean argument
    * indicates whether the level is the finest allowable level in the
    * hierarchy.  This flag is included since data management on the finest
    * level may be different than other levels in the hierarchy in some cases.
    *
    * The last two (optional) arguments specify an old level from which the
    * data may be used to initialize data on this level, and a flag that
    * indicates whether data on the initialized level must first be allocated.
    * The allocate_data argument is used in cases where one wishes to
    * simply reset data to an initialized state on a level that has already
    * been allocated.
    */
   virtual void
   initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true) = 0;

   /*!
    * After hierarchy levels have changed and data has been initialized on
    * the new levels, this routine can be used to reset any information
    * needed by the solution method that is particular to the hierarchy
    * configuration.  For example, the solution procedure may cache
    * communication schedules to amortize the cost of data movement on the
    * AMR patch hierarchy.  This function will be called by the gridding
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must
    * make the solution consistent across multiple levels after the hierarchy
    * is changed, this process may be invoked by this routine.  Of course the
    * details of these processes are determined by the particular solution
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void
   resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level) = 0;

   /*!
    * Set integer tags to "one" on the given level to identify
    * where refinement of that level should occur.  The index is that of the
    * cell-centered integer tag array on each patch.  The boolean argument
    * initial_time indicates whether cells are being tagged at
    * initialization time, or at some later time during the calculation.
    * If it is false, it should be assumed that the error estimation process
    * is being invoked at some later time after the AMR hierarchy was
    * initially constructed.  This information is provided since application
    * of the error estimator may be different in each of those circumstances.
    *
    * The cell-tagging operation may use time advancement to determine
    * tag regions. The argument coarsest_sync_level provides information
    * for the tagging method to coordinate time advance with an integrator.
    * When time integration is used during regridding, this value is true
    * if the level is the coarsest level involved in level synchronization
    * immediately preceeding the regrid process; otherwise it is false.
    * If time advancement is not used, this argument are ignored.
    *
    * The boolean can_be_refined is used to coordinate data reset operations
    * with the time integrator when time-dependent regridding is used.  This
    * is provided since data may be managed differently on the finest hierarchy
    * level than on coarser levels.
    */
   virtual void
   tagCellsForRefinement(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int regrid_cycle,
      const double regrid_time,
      const int tag_index,
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined = true,
      const double regrid_start_time = 0.) = 0;

   /*!
    * Certain cases may require pre-processing of error estimation data
    * before tagging cells, which is handled by this method. For example,
    * Richardson extrapolation may require advances of data in time before
    * the error estimation procedure is implemented.
    *
    * The level number indicates the level in which pre-process steps
    * are applied, time is the time at which the operation is performed
    * (generally the regrid time), and the boolean argument indicates
    * whether the operation is performed at the initial time.
    */
   virtual void
   preprocessErrorEstimation(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int cycle,
      const double regrid_time,
      const double regrid_start_time,
      const bool initial_time) = 0;

   /*!
    * Return true if regridding process advances the data using some time
    * integration procedure at the supplied cycle or time; otherwise, return
    * false.
    */
   virtual bool
   usesTimeIntegration(
      int cycle,
      double time) = 0;

   /*!
    * Return true if regridding process advances the data using some time
    * integration procedure at any cycle or time; otherwise, return false.
    */
   virtual bool
   everUsesTimeIntegration() const = 0;

   /*!
    * Return true if boxes for coarsest hierarchy level are not appropriate
    * for gridding strategy.  Otherwise, return false.  If false is returned,
    * it is useful to provide a detailed explanatory message describing the
    * problems with the boxes.
    */
   virtual bool
   coarsestLevelBoxesOK(
      const hier::BoxContainer& boxes) const = 0;

   /*!
    * Return ratio by which level may be coarsened during the error
    * estimation process.  Generally, this is needed by the gridding
    * algorithm class so that the new patch levels that it constructs can
    * be coarsened properly (if needed) during the error estimation process.
    */
   virtual int
   getErrorCoarsenRatio() const = 0;

   /*!
    * Check ratios between hierarchy levels against any constraints that
    * may be required for the error estimation scheme.
    */
   virtual void
   checkCoarsenRatios(
      const std::vector<hier::IntVector>& ratio_to_coarser) = 0;

   /*!
    * Return whether refinement is being performed using ONLY
    * user-supplied refine boxes.  If any method is used that invokes
    * tagging, this will return false.
    */
   virtual bool
   refineUserBoxInputOnly(
      int cycle,
      double time) = 0;

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Process a hierarchy before swapping old and new levels during
    * regrid.
    *
    * During regrid, if user code needs to do any application-specific
    * operations on the PatchHierarchy before a new level is added or
    * an old level is swapped for a new level, this method provides a callback
    * for the user to define such operations.  The PatchHierarchy is provided
    * in its state with the old level, if it exists, still in place, while
    * new BoxLevel is also provided so that the user code can know the boxes
    * that will make up the new level.
    *
    * @param hierarchy The PatchHierarchy being modified.
    * @param level_number The number of the PatchLevel in hierarchy being
    *                     added or regridded.
    * @param new_box_level BoxLevel containing the boxes for the new level
    *
    */
   virtual void
   processHierarchyBeforeAddingNewLevel(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const std::shared_ptr<hier::BoxLevel>& new_box_level) = 0;

   /*!
    * @brief Process a level before it is removed from the hierarchy during
    * regrid.
    *
    * In some cases user code may wish to process a PatchLevel before it is
    * removed from the hierarchy.  For example, data may exist only on a given
    * PatchLevel such as the finest level.  If that level were to be removed
    * before this data is moved off of it then the data will be lost.  This
    * method is a user defined callback used by GriddingAlgorithm when a
    * PatchLevel is to be removed.  The callback performs any user actions on
    * the level about to be removed.  It is implemented by classes derived from
    * StandardTagAndInitStrategy.
    *
    * @param hierarchy The PatchHierarchy being modified.
    * @param level_number The number of the PatchLevel in hierarchy about to be
    *                     removed.
    * @param old_level The level in hierarchy about to be removed.
    *
    * @see GriddingAlgorithm
    * @see StandardTagAndInitStrategy
    */
   virtual void
   processLevelBeforeRemoval(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>()) = 0;

   /*!
    * @brief Check the tags on a tagged level.
    *
    * This virtual interface provides application code a callback that
    * allows for checking the values held in user tag PatchData.  The
    * tag data will contain the tags created by application code in
    * tagCellsForRefinement as well as any tags added internally by
    * the GriddingAlgorithm (for example, buffering).
    *
    * A no-op implementation is provided so that only applications that
    * want to use this method need to implement it.
    *
    * @param[in] hierarchy
    * @param[in] level_number  Level number of the tagged level
    * @param[in] regrid_cycle
    * @param[in] regrid_time
    * @param[in] tag_index     Patch data index for user tags
    */
   virtual void
   checkUserTagData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int regrid_cycle,
      const double regrid_time,
      const int tag_index)
   {
      NULL_USE(hierarchy);
      NULL_USE(level_number);
      NULL_USE(regrid_cycle);
      NULL_USE(regrid_time);
      NULL_USE(tag_index);
   } 

   /*!
    * @brief Check the tags on a newly-created level.
    *
    * This virtual interface provides application code a callback that
    * allow for checking tag values that have been saved on a new level
    * that has been created during initialization or regridding.  The
    * tag values will be the values of the user tags on the coarser level,
    * constant-refined onto the cells of the new level.
    *
    * A no-op implementation is provided so that only applications that
    * want to use this method need to implement it.
    *
    * @param[in] hierarchy
    * @param[in] level_number   Level number of the new level
    * @param[in] tag_index      Patch data index for the new tags.
    */
   virtual void
   checkNewLevelTagData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const int tag_index)
   {
      NULL_USE(hierarchy);
      NULL_USE(level_number);
      NULL_USE(tag_index);
   } 

private:
   std::string d_object_name;

   // The following are not implemented:
   TagAndInitializeStrategy(
      const TagAndInitializeStrategy&);
   TagAndInitializeStrategy&
   operator = (
      const TagAndInitializeStrategy&);

};

}
}

#endif
