/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Tile clustering algorithm.
 *
 ************************************************************************/
#ifndef included_mesh_TileClustering
#define included_mesh_TileClustering

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Tiled patch clustering algorithm.
 *
 * Tiling generates clusters of a predetermined tile size.  Tile size can
 * be different on different levels, and mixing tiled and untiled levels
 * is permitted.  However, tiling is most efficient when tile boundaries
 * coincide, which is achieved by setting Tc*R is divisible by Tf or vice
 * versa, where Tc is the tile size on a coarser level, Tf is the tile size
 * on the finer level and R is the refinement ratio between the two levels.
 * Be sure to use a compatible tile size in the partitioning object.
 *
 * The algorithm is described in the article "Advances in Patch-Based
 * Adaptive Mesh Refinement Scalability" submitted to JPDC.  Scaling
 * benchmark results are also in the article.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 *   - \b tile_size
 *   Tile size in the index space of the tag level.
 *
 *   - \b coalesce_boxes
 *   Whether to coalesce boxes after clustering.  This can lead to
 *   clusters that are bigger than specified tile size.
 *
 *   - \b coalesce_boxes_from_same_patch
 *   Whether to coalesce tiled-boxes that originate from the
 *   same tag patch.
 *   This can reduce number of tiles and lead to clusters bigger than
 *   tile size.
 *
 *   - \b allow_remote_tile_extent
 *   Whether tile may extend to remote tag patches.
 *   If false, tiles will be cut at process boundaries, resulting in
 *   completely local tiles.  If true, allow tiles to cross process
 *   boundaries where, resulting in less tile fragmentation.
 *   If false, clusters' extent can depend on how tag level
 *   is partitioned.
 *
 *   - \b DEV_debug_checks
 *   Whether to run expensive checks for debugging.
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
 *     <td>tile_size</td>
 *     <td>int[]</td>
 *     <td>all values are 8</td>
 *     <td>all values > 0</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coalesce_boxes</td>
 *     <td>bool</td>
 *     <td>true</td>
 *     <td>false/true</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coalesce_boxes_from_same_patch</td>
 *     <td>bool</td>
 *     <td>false</td>
 *     <td>false/true</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>allow_remote_tile_extent</td>
 *     <td>bool</td>
 *     <td>true</td>
 *     <td>false/true</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>DEV_debug_checks</td>
 *     <td>bool</td>
 *     <td>false</td>
 *     <td>false/true</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.
 * Defaults are listed in parenthesis:
 *
 * @internal DEV_print_steps (FALSE)
 * boolean
 */
class TileClustering:public BoxGeneratorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   explicit TileClustering(
      const tbox::Dimension& dim,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   virtual ~TileClustering();

   /*!
    * @brief Implement the BoxGeneratorStrategy interface
    * method of the same name.
    *
    * Create a set of boxes that covers all integer tags on
    * the patch level that match the specified tag value.
    * Each box will be at least as large as the given minimum
    * size and the tolerances will be met.
    */
   void
   findBoxesContainingTags(
      std::shared_ptr<hier::BoxLevel>& new_box_level,
      std::shared_ptr<hier::Connector>& tag_to_new,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::BoxContainer& bound_boxes,
      const hier::IntVector& min_box,
      const hier::IntVector& max_gcw);

   void
   setRatioToNewLevel(const hier::IntVector& ratio)
   {
      d_ratio = ratio;
   }

   /*!
    * @brief Setup names of timers.
    */
   void
   setTimerPrefix(
      const std::string& timer_prefix);

protected:
   /*!
    * @brief Read parameters from input database.
    *
    * @param input_db Input Database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

private:
   /*!
    * @brief Cluster, cutting off tiles at process boundaries.
    *
    * This is a special implementation for when we do now allow tiles
    * to cross process boundaries.
    */
   void
   clusterWithinProcessBoundaries(
      hier::BoxLevel& new_box_level,
      hier::Connector& tag_to_new,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const hier::BoxContainer& bound_boxes,
      int tag_data_index,
      int tag_val);

   /*!
    * @brief Create, populate and return a coarsened version of the
    * given tag data.
    *
    * The coarse cell values are set to tag_data if any corresponding
    * fine cell value is tag_value.  Otherwise, the coarse cell value
    * is set to zero.
    */
   std::shared_ptr<pdat::CellData<int> >
   makeCoarsenedTagData(
      const pdat::CellData<int>& tag_data,
      int tag_value) const;

   /*!
    * @brief Find tagged tiles in a single patch.
    */
   int
   findTilesContainingTags(
      hier::BoxContainer& tiles,
      const pdat::CellData<int>& tag_data,
      int tag_val,
      int first_tile_index);

   /*!
    * @brief Cluster tags into whole tiles.  The tiles are not cut up,
    * even where they cross process boundaries or level boundaries.
    */
   void
   clusterWholeTiles(
      hier::BoxLevel& new_box_level,
      std::shared_ptr<hier::Connector>& tag_to_new,
      int& local_tiles_have_remote_extent,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const hier::BoxContainer& bound_boxes,
      int tag_data_index,
      int tag_val);

   /*!
    * @brief Detect semilocal edges missing from the outputs of
    * clusterWholeTiles().
    */
   void
   detectSemilocalEdges(
      std::shared_ptr<hier::Connector>& tag_to_tile);

   /*!
    * @brief Remove duplicate tiles created when a tile crosses a
    * tag-level process boundary.
    */
   void
   removeDuplicateTiles(
      hier::BoxLevel& tile_box_level,
      hier::Connector& tag_to_tiles);

   /*
    * @brief Shear tiles at block boundaries so they don't cross the boundaries.
    */
   void
   shearTilesAtBlockBoundaries(
      hier::BoxLevel& tile_box_level,
      hier::Connector& tag_to_tiles);

   /*!
    * @brief Coalesce clusters (and update Connectors).
    */
   void
   coalesceClusters(
      hier::BoxLevel& tile_box_level,
      std::shared_ptr<hier::Connector>& tag_to_tile,
      int tiles_have_remote_extent);

   /*!
    * @brief Special version of coalesceClusters.
    * for use when tiles have no remote extent.
    */
   void
   coalesceClusters(
      hier::BoxLevel& tile_box_level,
      std::shared_ptr<hier::Connector>& tag_to_tile);

   /*!
    * @brief Recursive bi-section version of BoxContainer::coalesce,
    * having O(N lg N) expected complexity.
    */
   void
   coalesceBoxes(
      hier::BoxContainer &boxes );

   const tbox::Dimension d_dim;

   //! @brief Tile size constraint.
   hier::IntVector d_tile_size;

   /*!
    * @brief Whether to allow tiles to have remote extents.
    */
   bool d_allow_remote_tile_extent;

   /*!
    * @brief Whether to coalesce all local tiled-boxes after
    * clustering.
    *
    * This can reduce number of tiles and lead to clusters bigger than
    * tile size.
    */
   bool d_coalesce_boxes;

   /*!
    * @brief Whether to coalesce tiled-boxes that originate from the
    * same tag patch.
    *
    * This can reduce number of tiles and lead to clusters bigger than
    * tile size.
    */
   bool d_coalesce_boxes_from_same_patch;

   /*!
    * @brief Number of boxes at which to use recursive, instead of
    * simple, coalesce.
    */
   int d_recursive_coalesce_limit;

   hier::IntVector d_ratio;

   /*!
    * @brief Thread locker for modifying clustering outputs with multi-threads.
    */
   TBOX_omp_lock_t l_outputs;

   /*!
    * @brief Thread locker for modifying intermediate data with multi-threads.
    */
   TBOX_omp_lock_t l_interm;

   //@{
   //! @name Diagnostics and performance evaluation
   bool d_debug_checks;
   hier::OverlapConnectorAlgorithm d_oca;
   hier::MappingConnectorAlgorithm d_mca;
   bool d_log_cluster_summary;
   bool d_log_cluster;
   bool d_barrier_and_time;
   bool d_print_steps;
   //@}

   //@{
   //! @name Performance timer data for this class.

   /*
    * @brief Structure of timers used by this class.
    *
    * Each object can set its own timer names through
    * setTimerPrefix().  This leads to many timer look-ups.  Because
    * it is expensive to look up timers, this class caches the timers
    * that has been looked up.  Each TimerStruct stores the timers
    * corresponding to a prefix.
    */
   struct TimerStruct {
      std::shared_ptr<tbox::Timer> t_find_boxes_containing_tags;
      std::shared_ptr<tbox::Timer> t_cluster;
      std::shared_ptr<tbox::Timer> t_cluster_local;
      std::shared_ptr<tbox::Timer> t_coalesce;
      std::shared_ptr<tbox::Timer> t_coalesce_adjustment;
      std::shared_ptr<tbox::Timer> t_global_reductions;
      std::shared_ptr<tbox::Timer> t_cluster_setup;
      std::shared_ptr<tbox::Timer> t_cluster_wrapup;
   };

   //! @brief Default prefix for Timers.
   static const std::string s_default_timer_prefix;

   /*!
    * @brief Static container of timers that have been looked up.
    */
   static std::map<std::string, TimerStruct> s_static_timers;

   static int s_primary_mpi_tag;
   static int s_secondary_mpi_tag;
   static int s_first_data_length;

   /*!
    * @brief Structure of timers in s_static_timers, matching this
    * object's timer prefix.
    */
   TimerStruct* d_object_timers;

   //@}

};

}
}

#endif  // included_mesh_TileClustering
