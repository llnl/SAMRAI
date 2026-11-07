/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines.
 *
 ************************************************************************/

#ifndef included_mesh_GriddingAlgorithm
#define included_mesh_GriddingAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithmConnectorWidthRequestor.h"
#include "SAMRAI/mesh/MultiblockGriddingTagger.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#define GA_RECORD_STATS
// #undef GA_RECORD_STATS

#ifdef GA_RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class GriddingAlgorithm manages AMR patch hierarchy construction
 * operations in SAMRAI.  Specifically, it provides AMR patch hierarchy
 * generation and regridding routines that may be used with a variety
 * of AMR solution algorithms and application codes.
 *
 * The three main functions provided by this class are:
 *   - @b    makeCoarsestLevel()
 *      This routine constructs or repartitions
 *      the coarsest hierarchy level (level 0).
 *
 *   - @b    makeFinerLevel()
 *      This routine will attempt to add a new
 *      finest level to the hierarchy if the
 *      maximum number of levels allows it and
 *      cells on the current finest level are
 *      tagged for refinement.
 *
 *   - @b    regridAllFinerLevels()
 *      This routine will regrid all levels finer
 *      than some specified level based on cells
 *      that are tagged for refinement on each
 *      level finer than and including the given
 *      level.  This routine may add a new finest
 *      hierarchy level if the maximum number of
 *      levels allows it and cells on the current
 *      finest level are tagged for refinement.
 *      Levels may also be removed from the
 *      hierarchy if no cells are tagged.
 *
 *
 * These basic AMR operations are used to generate levels in
 * the AMR patch hierarchy at the beginning of a simulation, and regridding
 * collections of levels during an adaptive calculation.  More details are
 * found in the comments accompanying each member function below.
 *
 * Other objects passed to the
 * constructor provide the gridding algorithm with particular operations
 * needed during meshing operations.  Operations that tag cells for
 * refinement on a patch level and initialize data on new levels are
 * provided by the TagAndInitializeStrategy argument.  Operations that
 * cluster tagged cells into a boxes are provided by the
 * BoxGeneratorStrategy argument.  Routines that load balance patches
 * on a level are provided by the LoadBalanceStrategy constructor argument.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *   - \b    check_overflow_nesting
 *
 *   - \b    check_proper_nesting
 *
 *   - \b    check_nonrefined_tags
 *      controls how to resolve user-specified tags that violate proper
 *      nesting.  If a tag violates the nesting requirements, its location in
 *      index space will not be refined when creating the finer level.  This
 *      flag allows the user to determine what to do when this occurs. <br>
 *      Set to one of these characters: <br>
 *      @b 'i' - violating tags will be quietly disregarded. <br>
 *      @b 'w' - violating tags will cause a warning and be disregarded. <br>
 *      @b 'e' - violating tags will cause an unrecoverable assertion. <br>
 *      It is fastest to ignore non-nesting tags because no checking has to be
 *      done.
 *
 *   - \b    check_overlapping_patches
 *      controls checking for overlapping patches on a new level. <br>
 *      Set to one of these characters: <br>
 *      @b 'i' - there is no check for overlapping patches, and they will be
 *      quietly disregarded. <br>
 *      @b 'w' - overlapping patches will cause a warning and be
 *      disregarded. <br>
 *      @b 'e' - violating tags will cause an unrecoverable assertion. <br>
 *      The check for overlapping patches may be and should be bypassed by
 *      applications that can tolerate overlaps.  To prevent the creation of
 *      levels with overlapping patches, see the PatchHierarchy input flag
 *      "allow_patches_smaller_than_minimum_size_to_prevent_overlaps".
 *
 *   - \b    check_nonnesting_user_boxes
 *      controls how user-specified refinement boxes that violate proper
 *      nesting are handled. <br>
 *      Set to one of these characters: <br>
 *      @b 'i' - nesting violations will be quietly disregarded. <br>
 *      @b 'w' - nesting violations will cause a warning but the code will
 *      continue anyway. <br>
 *      @b 'e' - nesting violations will cause an unrecoverable assertion <br>
 *      We highly recommend making nesting violation an error.  The code may
 *      work anyway, but there are no guarantees.
 *
 *   - \b    sequentialize_patch_indices
 *      whether patch indices will be globally sequentialized.
 *      This is required for writing correct VisIt files, and
 *      due to the current VisIt requirement, this is currently true by
 *      default.  It will evetually be set back to false after we remove the
 *      VisIt requirement.
 *
 *   - \b    enforce_proper_nesting
 * 
 *   - \b    save_tag_data
 *      This is an option to save the tags that are used to create a new
 *      fine level in CellData on that level.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>check_overflow_nesting</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>check_proper_nesting</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>check_nonrefined_tags</td>
 *     <td>string</td>
 *     <td>"WARN"</td>
 *     <td>"WARN", "IGNORE", "ERROR"</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>check_overlapping_patches</td>
 *     <td>string</td>
 *     <td>"IGNORE"</td>
 *     <td>"WARN", "IGNORE", "ERROR"</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>check_nonnesting_user_boxes</td>
 *     <td>string</td>
 *     <td>"ERROR"</td>
 *     <td>"WARN", "IGNORE", "ERROR"</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>sequentialize_patch_indices</td>
 *     <td>bool</td>
 *     <td>TRUE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>enforce_proper_nesting</td>
 *     <td>bool</td>
 *     <td>TRUE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db may be overridden by input db</td>
 *   </tr>
 *   <tr>
 *     <td>save_tag_data</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE, FALSE</td>
 *     <td>opt</td>
 *     <td>Parameter read from restart db will not be overridden by input db</td>
 *   </tr>
 * </table>
 *
 * All values read in from a restart database may be overriden by input
 * database values.  If no new input database value is given, the restart
 * database value is used.
 *
 * @see TagAndInitializeStrategy
 * @see LoadBalanceStrategy
 * @see BoxGeneratorStrategy
 */

class GriddingAlgorithm:
   public GriddingAlgorithmStrategy,
   public tbox::Serializable
{
public:

   /*!
    * @brief static symbols to hold special tag values for user tags
    *
    * These symbols are used to hold special values that are added to the
    * user tag data when the GriddingAlgorithm adds tags for algorithmic
    * reasons.
    *
    * EXISTING_FINE_PATCHES_TAG_VAL:  This value is placed in all cells
    * of the tag level than intersect the existing next finest level.
    * This value is added before calling the application's tagging methods,
    * so that the application code can take into account the location
    * of an existing finer level if it is relevant to its tagging criteria.
    *
    * NEW_FINE_PATCHES_TAG_VAL:  This value is placed in all cells of the
    * tag level that intersect the level that is two levels finer, if such
    * a level exists.  This value is added after the application's tagging
    * methods are invoked, and is used to ensure that all finer levels of
    * a hierarchy are nested inside coarser levels.
    *
    * BUFFER_TAG_VAL:  This value is added to cells that are tagged due to
    * the buffering operation of the GriddingAlgorithm, according to the
    * tag buffer arguments that are passed into certain methods of this
    * class.
    */
   static const int EXISTING_FINE_PATCHES_TAG_VAL;
   static const int NEW_FINE_PATCHES_TAG_VAL;
   static const int BUFFER_TAG_VAL;

   /*!
    * @brief The constructor for GriddingAlgorithm configures the
    * gridding algorithm with the patch hierarchy and concrete algorithm
    * strategy objects in the argument list.
    *
    * Gridding parameters are initialized from values provided in the
    * specified input and in the restart database corresponding to the
    * specified object_name argument.
    *
    * @param[in] hierarchy The hierarchy that this GriddingAlgorithm will
    * work on.  The pointer is cached.  All hierarchy operations will
    * be on this hierarchy.
    *
    * @param[in] object_name For registering the object in the restart
    * database.
    *
    * @param[in] input_db
    *
    * @param[in] level_strategy
    *
    * @param[in] generator
    *
    * @param[in] balancer Load balancer
    *
    * @param[in] balancer_zero Special load balancer to use for level
    * zero.  If omitted, will use @c balancer instead.
    *
    * @pre hierarchy
    * @pre !object_name.empty()
    * @pre tag_init_strategy
    * @pre generator
    * @pre balancer
    */
   GriddingAlgorithm(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const std::string& object_name,
      const std::shared_ptr<tbox::Database>& input_db,
      const std::shared_ptr<TagAndInitializeStrategy>& level_strategy,
      const std::shared_ptr<BoxGeneratorStrategy>& generator,
      const std::shared_ptr<LoadBalanceStrategy>& balancer,
      const std::shared_ptr<LoadBalanceStrategy>& balancer_zero =
         std::shared_ptr<LoadBalanceStrategy>());

   /*!
    * @brief Destructor
    */
   virtual ~GriddingAlgorithm();

   /*!
    * @brief Create or rebalance the coarsest level.
    *
    * This is an implementation of interface defined in GriddingAlgorithmStrategy.
    *
    * This routine will attempt to construct the coarsest level in the AMR
    * patch hierarchy (i.e., level 0).  If level 0 does not already exist,
    * then the domain specification is checked against the constraints of
    * the grid generation procedures.  Recall that the domain specification
    * is maintained by the grid geometry object associated with the hierarchy.
    * Generally, an unrecoverable assertion will result if the constraints
    * are not satisfied.
    *
    * If level 0 already exists in the hierarchy, then the routine
    * will generate a new level zero by re-applying the load balancing
    * procedure to the existing level.  Data will be moved from the
    * old level to the new level and the pre-existing level 0 will be
    * discarded.  Note that this routine is different than the routine
    * makeFinerLevel() below, which is used to construct levels finer
    * than level zero.  In particular, this routine does not select
    * cells for refinement, whereas the other routine does.
    *
    * @param[in] level_time Simulation time.
    *
    * @pre d_hierarchy->getMaxNumberOfLevels() > 0
    */
   void
   makeCoarsestLevel(
      const double level_time);

   /*!
    * @brief Attempt to create a new level in the hierarchy finer than
    * the finest level currently residing in the hierarchy.
    *
    * This is an implementation of interface method
    * GriddingAlgorithmStrategy::makeFinerLevel().
    *
    * The tag buffer indicates the number of cells by which cells
    * selected for refinement should be buffered before new finer
    * level boxes are constructed.  All tagged cells should be refined
    * except where refinement would violate proper nesting.  The
    * buffer is meant to keep phenomena of interest on refined regions
    * of the mesh until adaptive regridding occurs next.  Callers of
    * this method should take into account how the simulation may
    * evolve before regridding occurs (e.g., number of timesteps
    * taken) when calculating the tag_buffer.
    *
    * @param[in] tag_buffer See above text.
    *
    * @param[in] initial_cycle See above text
    *
    * @param[in] cycle See above text.
    *
    * @param[in] level_time See above text..
    *
    * @param[in] regrid_start_time The simulation time when the
    * regridding operation began (this parameter is ignored except
    * when using Richardson extrapolation)
    *
    * @pre d_hierarchy
    * @pre d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber())
    * @pre tag_buffer >= 0
    */
   void
   makeFinerLevel(
      const int tag_buffer,
      const bool initial_cycle,
      const int cycle,
      const double level_time,
      const double regrid_start_time = 0.0);

   /*!
    * @brief Attempt to regrid each level in the PatchHierarchy
    * that is finer than the specified level.
    *
    * This method implements the virtual interface
    * GriddingAlgorithmStrategy::regridAllFinerLevels().
    *
    * Note that the current algorithm permits at most one new finest level
    * to be added to the hierarchy with each invocation of the regridding
    * process.  This constraint, though seemingly restrictive makes the
    * process of maintaining properly nested levels much easier.
    *
    * Note that the current algorithm permits at most one new finest level
    * to be added to the hierarchy with each call to this method.
    * This constraint, though seemingly restrictive makes the
    * process of maintaining properly nested levels much easier.
    *
    * @param[in] level_number Coarsest level on which cells will be
    * tagged for refinement
    *
    * @param[in] tag_buffer Size of buffer on each level around tagged
    * cells that will be covered by the next finer level
    *
    * @param[in] cycle Simulaition cycle when regridding occurs
    *
    * @param[in] level_time Simulation time of the level corresponding to
    *                       level_num when regridding occurs
    *
    * @param[in] regrid_start_time The simulation time when the
    * regridding operation began on each level (this parameter is
    * ignored except when using Richardson extrapolation)
    *
    * @param[in] level_is_coarsest_to_sync Level is the coarsest to sync
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= d_hierarchy->getFinestLevelNumber())
    * @pre d_hierarchy->getPatchLevel(level_number)
    * @pre tag_buffer.size() >= level_number + 1
    * @pre for each member, tb, of tag_buffer, tb >= 0
    */
   void
   regridAllFinerLevels(
      const int level_number,
      const std::vector<int>& tag_buffer,
      const int cycle,
      const double level_time,
      const std::vector<double>& regrid_start_time = std::vector<double>(),
      const bool level_is_coarsest_to_sync = true);

   /*!
    * @brief Return pointer to level gridding strategy data member.
    *
    * Access to this member is useful when an integrator implementation
    * needs to know if the error estimator uses time integration.
    *
    * @return pointer to TagAndInitializeStrategy data member.
    */
   virtual
   std::shared_ptr<TagAndInitializeStrategy>
   getTagAndInitializeStrategy() const;

   /*!
    * @brief Return pointer to load balance strategy data member.
    */
   virtual
   std::shared_ptr<LoadBalanceStrategy>
   getLoadBalanceStrategy() const;

   /*!
    * @brief Return pointer to load balance strategy specialized for
    * balancing level zero.
    *
    * @return pointer to load balance strategy specialized for the case
    * where one processor owns all the initial loads.
    */
   virtual
   std::shared_ptr<LoadBalanceStrategy>
   getLoadBalanceStrategyZero() const;

   /*!
    * @brief Return pointer to PatchHierarchy data member.
    */
   std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const;

   /*!
    * @brief Return patch data index for tag data saved in the hierarchy.
    */
   int
   getSavedTagPatchDataIndex() const
   {
      return d_saved_tag_indx;
   }

   /*!
    * @brief Print all data members of the class instance to given output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Write object state out to the given restart database.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*
    * @brief Write out statistics recorded on numbers of cells and patches generated.
    */
   void
   printStatistics(
      std::ostream& s = tbox::plog) const;

   /*!
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int ALGS_GRIDDING_ALGORITHM_VERSION;

   //! @brief Shorthand typedef.
   typedef hier::Connector::NeighborSet NeighborSet;

   /*!
    * @brief Read input data from specified database and initialize class members.
    *
    * The database pointer must be non-null.
    *
    * @param[in] input_db
    *
    * @param[in] is_from_restart
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db,
      bool is_from_restart);

   /*!
    * @brief Read object state from the restart file and initialize
    * class data members.
    *
    * The database from which the restart data is read is determined
    * by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *   -The database corresponding to object_name is not found
    *    in the restart file.
    *
    *   -The class version number and restart version number do not
    *    match.
    */
   void
   getFromRestart();

   /*!
    * @brief Recursively regrid the hierarchy level and all finer
    * levels in the hierarchy.
    *
    * This private member function is invoked by the
    * regridAllFinerLevels() routine.  It may invoke recursively
    * invoke itself.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= d_hierarchy->getFinestLevelNumber())
    * @pre d_hierarchy->getPatchLevel(level_number)
    * @pre (finest_level_not_regridded >= 0) &&
    *      (finest_level_not_regridded <= level_number)
    * @pre tag_buffer.size() >= level_number + 1
    * @pre for each member, tb, of tag_buffer, tb >= 0
    */
   void
   regridFinerLevel(
      const int level_number,
      const double regrid_time,
      const int regrid_cycle,
      const int finest_level_not_regridded,
      const bool level_is_coarsest_to_sync,
      const std::vector<int>& tag_buffer,
      const std::vector<double>& regrid_start_time = std::vector<double>(0));

   /*!
    * @brief Tagging stuff before recursive regrid, called from regridFinerLevel.
    */
   void
   regridFinerLevel_doTaggingBeforeRecursiveRegrid(
      const int tag_ln,
      const bool level_is_coarsest_sync_level,
      const std::vector<double>& regrid_start_time,
      const double regrid_time,
      const int regrid_cycle);

   /*!
    * @brief Tagging stuff after recursive regrid, called from regridFinerLevel.
    *
    * A side-effect of this method is setting the overlap Connectors
    * between the tag level (number tag_ln) and the finer level
    * (number tag_ln+2) if the finer level exists in the hierarchy.
    */
   void
   regridFinerLevel_doTaggingAfterRecursiveRegrid(
      std::shared_ptr<hier::Connector>& tag_to_finer,
      const int tag_ln,
      const std::vector<int>& tag_buffer,
      double regrid_time);

   /*!
    * @brief Given the metadata describing the new level, this method
    * creates and installs new PatchLevel in the hierarchy.
    *
    * @pre tag_to_new && tag_to_new->hasTranspose()
    * @pre new_box_level
    */
   void
   regridFinerLevel_createAndInstallNewLevel(
      const int tag_ln,
      const double regrid_time,
      std::shared_ptr<hier::Connector>& tag_to_new,
      std::shared_ptr<hier::BoxLevel> new_box_level);

   /*!
    * @brief Set all tags on a level to a given value.
    *
    * @param[in] tag_value
    *
    * @param[in] tag_level
    *
    * @param[in] tag_index
    *
    * @pre (tag_value == d_true_tag) || (tag_value == d_false_tag)
    * @pre tag_level
    * @pre (tag_index == d_tag_indx) || (tag_index == d_buf_tag_indx)
    */
   void
   fillTags(
      const int tag_value,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_index) const;

   /*!
    * @brief Set integer tags to specified value on intersection
    * between patch level and a BoxLevel
    *
    * The index value corresponds to the patch descriptor entry of the
    * cell-centered integer tag array.  The boolean flag indicates
    * whether tags are to be set on the regions corresponding to the
    * interior of the level only, if the tag arrays contain ghosts.
    *
    * @param[in] tag_value
    *
    * @param[in] tag_level
    *
    * @param[in] index
    *
    * @param[in] tag_level_to_fill_box_level Connector from the
    * level with the tags to the BoxLevel describing where to
    * fill.
    *
    * @param[in] interior_only
    *
    * @param[in] fill_box_growth How much to grow fill boxes before using
    * them to tag.  Must be in index space of level holding fill boxes.
    *
    * @param[in] preserve_existing_tags If set to true, only cells that
    * have a value of d_false_tag on input will be changed to tag_value. If
    * false, then all cells intersecting by the fill BoxLevel will be set to
    * tag_value. 
    *
    * @pre (tag_value == d_true_tag) || (tag_value == d_false_tag)
    * @pre tag_level
    * @pre (index == d_tag_indx || index == d_buf_tag_indx)
    * @pre tag_level_to_fill_box_level.getHeadCoarserFlag() == false
    */
   void
   fillTagsFromBoxLevel(
      const int tag_value,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const int index,
      const hier::Connector& tag_level_to_fill_box_level,
      const bool interior_only,
      const hier::IntVector& fill_box_growth,
      const bool preserve_existing_tags) const;

   /*!
    * @brief Enforce proper nesting.
    *
    * @param[in,out] new_box_level
    *
    * @param[in,out] tag_to_new
    *
    * @param[in] tag_ln
    */
   void
   enforceProperNesting(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      int tag_ln) const;

   /*!
    * @brief Make a map that, when applied to an improperly nested
    * BoxLevel, removes the nonnesting parts.
    *
    * @param[out] nested_box_level  The nesting parts of the
    * unnested BoxLevel.
    *
    * @param[out] unnested_to_nested  The mapping from the unnested
    * BoxLevel to the nested BoxLevel.
    *
    * @param[in] unnested_box_level
    *
    * @param[in] tag_to_unnested  Overlap Connector from the tag level
    * to unnested_box_level.
    *
    * @param[in] unnested_ln Level number of PatchLevel being
    * generated (one more than the tag level number).
    *
    * @param[in] oca
    *
    * @pre tag_to_unnested.hasTranspose()
    * @pre d_hierarchy->getDim() == unnested_box_level.getDim()
    */
   void
   makeProperNestingMap(
      std::shared_ptr<hier::BoxLevel>& nested_box_level,
      std::shared_ptr<hier::MappingConnector>& unnested_to_nested,
      const hier::BoxLevel& unnested_box_level,
      const hier::Connector& tag_to_unnested,
      const int unnested_ln,
      const hier::OverlapConnectorAlgorithm& oca) const;

   /*!
    * @brief Enforce overflow nesting.
    *
    * @param[in,out] new_box_level
    *
    * @param[in,out] tag_to_new
    */
   void
   enforceOverflowNesting(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new) const;

   /*!
    * @brief Make a map that, when applied to a BoxLevel that
    * extends outside of a reference BoxLevel, removes those
    * outside parts.
    *
    * @param[out] nested_box_level  The nesting parts of the unnested BoxLevel.
    *
    * @param[out] unnested_to_nested  The mapping from the unnested
    * BoxLevel to the nested BoxLevel.
    *
    * @param[in] unnested_box_level
    *
    * @param[in] unnested_to_reference
    *
    * @pre d_hierarchy->getDim() == unnested_box_level.getDim()
    */
   void
   makeOverflowNestingMap(
      std::shared_ptr<hier::BoxLevel>& nested_box_level,
      std::shared_ptr<hier::MappingConnector>& unnested_to_nested,
      const hier::BoxLevel& unnested_box_level,
      const hier::Connector& unnested_to_reference) const;

   /*!
    * @brief Make a map from a BoxLevel to parts of that BoxLevel
    * that violate proper nesting.
    *
    * @param[out] violator BoxLevel containing violating parts of candidate.
    *
    * @param[in] candidate_to_violator Connector between candidate and violator.
    *
    * @param[in] candidate BoxLevel being examined for nesting violation.
    *
    * @param[in] candidate_to_hierarchy Connector to box_level number
    *       tag_ln in the hierarchy.
    *
    * @param[in] tag_ln Level number of the level that candidate should nest in.
    *
    * @param[in] oca
    *
    * @pre candidate_to_hierarchy.hasTranspose()
    * @pre d_hierarchy->getDim() == candidate.getDim()
    * @pre candidate_to_hierarchy.getRatio() == hier::IntVector::getOne(d_hierarchy->getDim())
    * @pre candidate_to_hierarchy.getTranspose().getRatio() == hier::IntVector::getOne(d_hierarchy->getDim())
    */
   void
   computeNestingViolator(
      std::shared_ptr<hier::BoxLevel>& violator,
      std::shared_ptr<hier::MappingConnector>& candidate_to_violator,
      const hier::BoxLevel& candidate,
      const hier::Connector& candidate_to_hierarchy,
      const int tag_ln,
      const hier::OverlapConnectorAlgorithm& oca) const;

   /*!
    * @brief Extend Boxes to domain boundary if they they are too close.
    *
    * See hier::BoxUtilities::extendBoxToDomainBoundary().
    *
    * @param[in,out] new_box_level
    *
    * @param[in,out] tag_to_new
    *
    * @param[in] physical_domain_array  Vector holding domain for each block
    *
    * @param[in] extend_ghosts Extend the boxes to the boundary if
    * they less than extend_ghosts cells from the boundary.
    *
    * @pre tag_to_new.hasTranspose()
    */
   void
   extendBoxesToDomainBoundary(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      const std::vector<hier::BoxContainer>& physical_domain_array,
      const hier::IntVector& extend_ghosts) const;

   /*!
    * @brief Precompute data used to define proper nesting.
    *
    * Data is associated with level number ln, to be used for
    * constructing level number ln+1.
    *
    * @param[in] ln
    *
    * @param[in] oca
    *
    * @pre (d_base_ln >= 0) && (ln >= d_base_ln)
    * @pre (ln == d_base_ln) || (d_to_nesting_complement[ln - 1]->isFinalized())
    */
   void
   computeProperNestingData(
      const int ln,
      const hier::OverlapConnectorAlgorithm& oca);

   /*!
    * @brief Attempt to fix boxes that are too small by growing them
    * within the nesting-restricted domain.
    *
    * @param[in,out] new_box_level BoxLevel being formed
    * into the new new level.
    *
    * @param[in,out] tag_to_new  Connector to be updated.
    *
    * @param[in] min_size
    *
    * @param[in] tag_ln Level number of the tag level.
    *
    * @pre tag_no_new.hasTranspose()
    * @pre (d_hierarchy->getDim() == new_box_level.getDim()) &&
    *      (d_hierarchy->getDim() == min_size.getDim())
    */
   void
   growBoxesWithinNestingDomain(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      const hier::IntVector& min_size,
      const int tag_ln) const;

   /*!
    * @brief Refine a BoxLevel from the resolution of the tag
    * level to the resolution of the level being created.
    *
    * @param[in,out] new_box_level BoxLevel being formed
    * into the new new level.
    *
    * @param[in,out] tag_to_new  Connector to be updated.
    *
    * @param[in] ratio
    *
    * @pre tag_to_new.hasTranspose()
    */
   void
   refineNewBoxLevel(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      const hier::IntVector& ratio) const;

   /*!
    * @brief Renumber Boxes in a BoxLevel.
    *
    * This method renumbers Boxes to give them globally
    * sequential LocalIds and/or sorts the Boxes by their
    * coordinates.  Sequentializing the global indices numbers them
    * sequentially, like Patch numbering in the traditional SAMR
    * approach.  Sorting by box coordinates removes the randomness in
    * the ordering of boxes.
    *
    * @param[in,out] new_box_level
    *
    * @param[in,out] ref_to_new
    *
    * @param[in] sort_by_corners
    *
    * @param[in] sequentialize_global_indices
    *
    * @pre d_hierarchy->getDim() == new_box_level.getDim()
    */
   void
   renumberBoxes(
      hier::BoxLevel& new_box_level,
      hier::Connector* ref_to_new,
      bool sort_by_corners,
      bool sequentialize_global_indices) const;

   /*!
    * @brief Buffer each integer tag on patch level matching given tag
    * value with a border of matching tags.
    *
    * @pre (tag_value == d_true_tag) || (tag_value == d_false_tag)
    * @pre level
    * @pre buffer_size >= 0
    * @pre d_hierarchy->getDim() == level->getDim()
    */
   void
   bufferTagsOnLevel(
      const int tag_value,
      const std::shared_ptr<hier::PatchLevel>& level,
      const int buffer_size) const;

   /*!
    * @brief Set the new level boxes using information stored in a file.
    *
    * If cell tagging is not performed, the new level boxes may
    * be specified either from user input (by specifying "REFINE_BOXES"
    * input) or from output from a previous run stored in an
    * HDF file.  This method sets the "new_level_boxes" based on
    * the information in the file.  It also sets the boolean
    * parameter "remove_old_fine_level" which indicates whether
    * the level box configuration has changed and, consequently,
    * whether we need to remove the old level.
    *
    * @pre (level_number >= 0) &&
    *      (level_number <= d_hierarchy->getFinestLevelNumber())
    */
   void
   readLevelBoxes(
      std::shared_ptr<hier::BoxLevel>& new_box_level,
      std::shared_ptr<hier::Connector>& coarser_to_new,
      const int level_number,
      const double regrid_time,
      const int regrid_cycle,
      bool& remove_old_fine_level);

   /*!
    * @brief Given the number for the level where cells are tagged for
    * refinement, compute a BoxLevel from which a refinement of
    * the level may be constructed.
    *
    * It is assumed that the integer tags that identify cells for
    * refinement have been set on the level before this routine is
    * called.  At the end of this function, new_box_level will
    * represent the box extents of a new fine level on the given
    * block.
    *
    * @pre (tag_ln >= 0) && (tag_ln <= d_hierarchy->getFinestLevelNumber())
    * @pre d_hierarchy->getPatchLevel(tag_ln)
    */
   void
   findRefinementBoxes(
      std::shared_ptr<hier::BoxLevel>& new_box_level,
      std::shared_ptr<hier::Connector>& tag_to_new,
      const int tag_ln) const;

   /*!
    * @brief Set patch size and ghost cell information needed to
    * create new refinement levels.
    *
    * This method applies to levels that are being used to build new
    * finer levels (i.e. level_number is a coarser level in the
    * hierarchy) and to levels which are simply being reconstructed
    * (i.e. the same level in the hierarchy).  The boolean value
    * "for_building_finer" controls the logic for the two cases - in
    * the former case, it is true while in the latter case it is
    * false.
    *
    * When a finer level is being constructed, the maximum number of ghost
    * cells needed in any variable is used to compute the smallest patch
    * size allowed and the extent to which patches may be extended to touch
    * the physical boundary.  This avoids problems in setting ghost cell
    * data that may occur when ghost cell regions intersect the physical
    * boundary in strange ways.
    *
    * This routine sets the smallest and largest patch sizes for the specified
    * level, the smallest box to refine on the next coarser level, and the
    * number of cells that a patch may be extended to the boundary if it
    * sufficiently close to the boundary (extend_ghosts).
    *
    * @pre (d_hierarchy->getDim() == smallest_patch.getDim()) &&
    *      (d_hierarchy->getDim() == smallest_box_to_refine.getDim()) &&
    *      (d_hierarchy->getDim() == largest_patch.getDim()) &&
    *      (d_hierarchy->getDim() == extend_ghosts.getDim())
    * @pre (level_number >= 0) &&
    *      (level_number < d_hierarchy->getMaxNumberOfLevels())
    */
   void
   getGriddingParameters(
      hier::IntVector& smallest_patch,
      hier::IntVector& smallest_box_to_refine,
      hier::IntVector& largest_patch,
      hier::IntVector& extend_ghosts,
      const int level_number,
      const bool for_building_finer) const;

   /*!
    * @brief Compute d_tag_to_cluster_width.
    */
   void
   computeTagToClusterWidths();

   /*!
    * @brief Check domain boxes for violations of certain constraints.
    */
   void
   checkDomainBoxes(
      const hier::BoxContainer& domain_boxes) const;

   /*!
    * @brief Check for non-nesting user-specified boxes.
    */
   void
   checkNonnestingUserBoxes(
      const hier::Connector& new_to_tag,
      const hier::IntVector& nesting_buffer) const;

   /*!
    * @brief Check for boxes that are too close to the physical
    * boundary without touching it.
    *
    * Boxes close to the physical boundaries causes ghost boxes to
    * intersect the boundary in weird ways, so we disallow it.  This
    * method writes a warning describing each violation found.
    *
    * @param[in] box_level
    *
    * @param[in] extend_ghosts
    *
    * Return the number of violations found.
    */
   size_t
   checkBoundaryProximityViolation(
      const hier::BoxLevel& box_level,
      const hier::IntVector& extend_ghosts) const;

   /*!
    * @brief Check for boxes that are too close to the physical
    * boundary without touching it.
    *
    * Compute the allowable distance from boxes in new_box_level
    * to domain boundary and delegate the checking.
    */
   void
   checkBoundaryProximityViolation(
      const int tag_ln,
      const hier::BoxLevel& new_box_level) const;

   /*!
    * @brief Set tag data that will be only have the values expected by
    * the algorithms for box generation.
    *
    * Boolean tag data should have only values of d_true_tag or
    * d_false_tag.  This method interprets user tag data to set
    * boolean tag data to these values.
    *
    * Note that the use of the term "boolean" does not mean that the bool
    * type is used, but that the integer data will have only two possible
    * values.
    *
    * @param[in] tag_level  Level being tagged
    * @param[in] preserve_existing_tags  If set to true, any cells already
    * set to d_true_tag in the algorthmic tag data will be guaranteed to
    * remain so.  If false, no current values of the boolean tag data
    * will be guaranteed to be preserved.
    */
   void
   setBooleanTagData(
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      bool preserve_existing_tags) const;

   /*!
    * @brief Check for user tags that violate proper nesting.
    *
    * @param tag_level  Tag level
    *
    * @param tag_ln  Tag level number
    *
    * @param oca
    *
    * @pre d_hierarchy->getDim() == tag_level.getDim()
    */
   void
   checkNonrefinedTags(
      const hier::PatchLevel& tag_level,
      int tag_ln,
      const hier::OverlapConnectorAlgorithm& oca) const;

   /*!
    * @brief Reset data that handles tag buffering.
    *
    * Resets the tag buffering data so that it will be able to handle a
    * a tag buffer of the given size.
    *
    * @param tag_buffer  The size of the buffer
    */
   void
   resetTagBufferingData(
      const int tag_buffer);

   /*!
    * @brief Check for overlapping patches within a level when you
    * have the self Connector for the level.
    *
    * @param[in] box_level_to_self
    */
   void
   checkOverlappingPatches(
      const hier::Connector& box_level_to_self) const;

   /*!
    * @brief Check for overlapping patches within the level when you
    * do not have the self Connector for the level.
    */
   void
   checkOverlappingPatches(
      const hier::BoxLevel& box_level) const;

   /*!
    * @brief Check if the tag buffer data needs to be reset to a default value
    * at the conclusion of tagging operations.  This is true when the
    * tag buffer is larger than maximum ghost width for any
    * currently registered data.
    */
   bool
   needResetTagBuffer(const hier::IntVector& max_ghosts) const;

   /*!
    * @brief Warn if the domain is too small any periodic direction.
    */
   void
   warnIfDomainTooSmallInPeriodicDir() const;

   /*!
    * @brief Allocate timers.
    */
   void
   allocateTimers();

   /*!
    * @brief Assert that the given SAMRAI_MPI has no message waiting
    * to be received.
    */
   void
   assertNoMessageForCommunicator(
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * @brief Record statistics on how many patches and cells were generated.
    */
   void
   recordStatistics(
      double current_time);

   /*!
    * @brief The hierarchy that this GriddingAlgorithm works on.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Implementation registered with the hierarchy, telling the
    * hierarchy what width the GriddingAlgorithm will be requesting.
    */
   GriddingAlgorithmConnectorWidthRequestor d_connector_width_requestor;

   hier::IntVector d_buf_tag_ghosts;

   /*
    * The object name is used for error reporting and accessing
    * restart file information.
    */
   std::string d_object_name;

   /*
    * Data members that manage application-specific level initialization
    * and cell tagging, clustering of tagged cells into boxes, and load
    * balancing of patches to processors, respectively.
    */
   std::shared_ptr<TagAndInitializeStrategy> d_tag_init_strategy;
   std::shared_ptr<BoxGeneratorStrategy> d_box_generator;
   std::shared_ptr<LoadBalanceStrategy> d_load_balancer;
   std::shared_ptr<LoadBalanceStrategy> d_load_balancer0;

   /*
    * MultiblockGriddingTagger is the RefinePatchStrategy
    * implementation provided to the RefineSchedule d_bdry_sched_tags
    * inside GriddingAlgorithm.
    *
    * @see setMultiblockGriddingTagger().
    */
   MultiblockGriddingTagger * d_mb_tagger_strategy;

   /*!
    * @brief Variables and patch data indices for tags.
    *
    * Cell-centered integer variables use to tag cells for refinement.
    * The descriptor index d_user_tag_indx is used to obtain tag information
    * from user-defined routines on patch interiors.  The descriptor index
    * d_buf_tag_indx is used to buffer tags on patches that may be
    * distributed across processors.  d_boolean_tag_indx is used to put tag
    * values in a standard format that will be understood by box generator
    * implementations, and d_saved_tag_indx is used to preserve tag values
    * on new levels after gridding operations are completed.  Note that
    * d_boolean_tag does not use the bool type, but rather its associated
    * integer data is treated as a boolean with only two possible values.
    */
   std::shared_ptr<pdat::CellVariable<int> > d_user_tag;
   std::shared_ptr<pdat::CellVariable<int> > d_saved_tag;
   std::shared_ptr<pdat::CellVariable<int> > d_boolean_tag;
   std::shared_ptr<pdat::CellVariable<int> > d_buf_tag;
   int d_user_tag_indx;
   int d_saved_tag_indx;
   int d_boolean_tag_indx;
   int d_buf_tag_indx;

   std::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_tags;
   std::vector<std::shared_ptr<xfer::RefineSchedule> > d_bdry_sched_tags;

   /*!
    * @brief Refine algorithm and schedule for filling saved tag data on new
    * levels.
    */
   std::shared_ptr<xfer::RefineAlgorithm> d_fill_saved_tags;
   std::shared_ptr<xfer::RefineSchedule> d_saved_tags_sched;

   /*!
    * @brief Tag values defined in the constructor.
    * 
    * Tag values defined in the constructor and used to set and compare
    * integer tag values on the patch hierarchy.  Generally, these
    * these variables are easier to read in the code than integer constants,
    * such as 0 and 1.
    */
   const int d_true_tag;
   const int d_false_tag;
   const int d_from_fine_pretag;
   const int d_from_fine_tag;
   const int d_buffer_tag;

   /*!
    * @brief Finest level not changed during regrid.
    *
    * This member is temporary, used to coordinate with private methods
    * during mesh generation and regridding.  It has a value of -1 when
    * not being used.
    */
   int d_base_ln;

   /*!
    * @brief Connector widths to use when clustering.
    */
   std::vector<hier::IntVector> d_tag_to_cluster_width;

   /*
    * @brief When regridding level ln+1, the new level ln must not flow into
    * d_proper_nesting_complement[ln].
    *
    * Has length d_hierarchy->getMaxNumberOfLevels().  The objects are
    * initialized only during gridding/regridding.
    */
   std::vector<std::shared_ptr<hier::BoxLevel> > d_proper_nesting_complement;

   /*
    * @brief Connectors from the hierarchy to d_proper_nesting_complement.
    *
    * d_to_nesting_complement[ln] goes from level ln to
    * d_proper_nesting_complement[ln].
    */
   std::vector<std::shared_ptr<hier::Connector> > d_to_nesting_complement;

   /*!
    * @brief How to resolve user tags that violate nesting requirements.
    *
    * See input parameter "check_nonrefined_tags".
    */
   char d_check_nonrefined_tags;

   /*!
    * @brief Whether to check for overlapping patches on a level.
    *
    * See input parameter "check_overlapping_patches".
    */
   char d_check_overlapping_patches;

   /*!
    * @brief Whether to check for non-nesting user-specified refine
    * boxes.
    *
    * See input parameter "check_nonnesting_user_boxes".
    */
   char d_check_nonnesting_user_boxes;

   /*!
    * @brief Whether to check for user-specified refine boxes that
    * violate boundary proximity.
    *
    * See input parameter "check_boundary_proximity_violation".
    */
   char d_check_boundary_proximity_violation;

   /*!
    * @brief Whether to globally sequentialize patch indices on every level.
    */
   bool d_sequentialize_patch_indices;

   /*!
    * @brief Whether to log metadata statistics after generating a new
    * level.
    *
    * See input parameter description.
    */
   bool d_log_metadata_statistics;

   /*!
    * @brief OverlapConnectorAlgorithm object used for regrid.
    */
   mutable hier::OverlapConnectorAlgorithm d_oca;

   /*!
    * @brief MappingConnectorAlgorithm object used for regrid.
    */
   mutable hier::MappingConnectorAlgorithm d_mca;

   /*!
    * @brief BoxLevelConnectorUtils object used for regrid.
    */
   mutable hier::BoxLevelConnectorUtils d_blcu;

   /*!
    * @brief OverlapConnectorAlgorithm object used for initial mesh
    * construction and other one-time operations that are not expected
    * to scale well.
    */
   mutable hier::OverlapConnectorAlgorithm d_oca0;

   /*!
    * @brief MappingConnectorAlgorithm object used for initial mesh
    * construction and other one-time operations that are not expected
    * to scale well.
    */
   mutable hier::MappingConnectorAlgorithm d_mca0;

   /*!
    * @brief BoxLevelConnectorUtils object used for initial mesh
    * construction and other one-time operations that are not expected
    * to scale well.
    */
   mutable hier::BoxLevelConnectorUtils d_blcu0;

   /*
    * Switches for massaging boxes after clustering.  Should be on for
    * most AMR applications.  Turning off is mainly for debugging
    * purposes.
    */
   bool d_enforce_proper_nesting;
   bool d_extend_to_domain_boundary;
   bool d_load_balance;

   /*
    * Flag to save tag data on new fine levels.
    */
   bool d_save_tag_data; 

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_and_time;
   //@}

   /*
    * Timers interspersed throughout the class.
    */
   std::shared_ptr<tbox::Timer> t_find_domain_complement;
   std::shared_ptr<tbox::Timer> t_load_balance;
   std::shared_ptr<tbox::Timer> t_load_balance0;
   std::shared_ptr<tbox::Timer> t_bdry_fill_tags_create;
   std::shared_ptr<tbox::Timer> t_make_coarsest;
   std::shared_ptr<tbox::Timer> t_make_finer;
   std::shared_ptr<tbox::Timer> t_regrid_all_finer;
   std::shared_ptr<tbox::Timer> t_regrid_finer_do_tagging_before;
   std::shared_ptr<tbox::Timer> t_regrid_finer_do_tagging_after;
   std::shared_ptr<tbox::Timer> t_regrid_finer_create_and_install;
   std::shared_ptr<tbox::Timer> t_regrid_finer_create;
   std::shared_ptr<tbox::Timer> t_fill_tags;
   std::shared_ptr<tbox::Timer> t_tag_cells_for_refinement;
   std::shared_ptr<tbox::Timer> t_initialize_level_data;
   std::shared_ptr<tbox::Timer> t_buffer_tags;
   std::shared_ptr<tbox::Timer> t_bdry_fill_tags_comm;
   std::shared_ptr<tbox::Timer> t_second_finer_tagging;
   std::shared_ptr<tbox::Timer> t_find_refinement;
   std::shared_ptr<tbox::Timer> t_bridge_new_to_new;
   std::shared_ptr<tbox::Timer> t_find_new_to_new;
   std::shared_ptr<tbox::Timer> t_bridge_new_to_finer;
   std::shared_ptr<tbox::Timer> t_bridge_new_to_old;
   std::shared_ptr<tbox::Timer> t_find_boxes_containing_tags;
   std::shared_ptr<tbox::Timer> t_fix_zero_width_clustering;
   std::shared_ptr<tbox::Timer> t_enforce_proper_nesting;
   std::shared_ptr<tbox::Timer> t_compute_proper_nesting_data;
   std::shared_ptr<tbox::Timer> t_make_nesting_map;
   std::shared_ptr<tbox::Timer> t_use_nesting_map;
   std::shared_ptr<tbox::Timer> t_make_overflow_map;
   std::shared_ptr<tbox::Timer> t_use_overflow_map;
   std::shared_ptr<tbox::Timer> t_compute_external_parts;
   std::shared_ptr<tbox::Timer> t_compute_nesting_violator;
   std::shared_ptr<tbox::Timer> t_extend_to_domain_boundary;
   std::shared_ptr<tbox::Timer> t_extend_within_domain;
   std::shared_ptr<tbox::Timer> t_grow_boxes_within_domain;
   std::shared_ptr<tbox::Timer> t_renumber_boxes;
   std::shared_ptr<tbox::Timer> t_make_domain;
   std::shared_ptr<tbox::Timer> t_make_new;
   std::shared_ptr<tbox::Timer> t_process_error;
   std::shared_ptr<tbox::Timer> t_enforce_overflow_nesting;
   std::shared_ptr<tbox::Timer> t_reset_hier;

#ifdef GA_RECORD_STATS
   /*
    * Statistics on number of cells and patches generated.
    */
   std::vector<std::shared_ptr<tbox::Statistic> > d_boxes_stat;
   std::vector<std::shared_ptr<tbox::Statistic> > d_cells_stat;
   std::vector<std::shared_ptr<tbox::Statistic> > d_timestamp_stat;
#endif

   // The following are not yet implemented:
   GriddingAlgorithm(
      const GriddingAlgorithm&);
   GriddingAlgorithm&
   operator = (
      const GriddingAlgorithm&);

   // Verbose flags.
   bool d_check_overflow_nesting;
   bool d_check_proper_nesting;
   bool d_check_connectors;
   bool d_print_steps;

};

}
}

#endif
