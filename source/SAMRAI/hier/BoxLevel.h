/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of Boxes in the same "level".
 *
 ************************************************************************/
#ifndef included_hier_BoxLevel
#define included_hier_BoxLevel

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelHandle.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PersistentOverlapConnectors.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*
 *****************************************************************************
 * IMPORTANT
 * IF YOU ARE ADDING A NON-CONST METHOD IT MUST BE PREFACED WITH THE FOLLOWING
 * CODE AS LOCKED BOX_LEVELS MAY NOT BE MODIFIED IN ANY WAY
 * if (locked()) {
 *    TBOX_ERROR("BoxLevel::newMethodName(): operating on locked BoxLevel."
 *       << std::endl);
 * }
 *****************************************************************************
 */

/*!
 * @brief A distributed set of Box objects which reside in the
 * same index space.
 *
 *
 * TODO: Are we eliminating DLBG terminology?
 *
 * This class is a part of the distributed layered box graph (DLBG) for
 * managing SAMR meshes in parallel.  A BoxLevel is a set of
 * boxes in the same index space. Relationships (e.g., neighbor adjacency)
 * among boxes is contained in a Connector object. Also, each BoxLevel
 * has an refinement ratio vector describing the relationship of the
 * index space to that of a reference level in a patch hierarchy (typically
 * the coarsest level or level zero).
 *
 * Like a PatchLevel, a BoxLevel is a parallel object.  The
 * Boxes of a BoxLevel may be distributed across all the
 * processors in an MPI communicator and can be in one
 * of two parallel states:
 *
 * - @b DISTRIBUTED: Each MPI process knows only the Boxes in the set
 * that are "owned" by that process.  This is analogous to a PatchLevel
 * which owns only the Patches that reside on a process.
 *
 * - @b GLOBALIZED: All processes know all Boxes in the set.
 * This is analogous to PatchLevel BoxContainer state when it is globalized
 * (@see PatchLevel::getBoxes()).
 *
 * @par Performance notes
 * <ul>
 * <li> The parallel state is changed by calling setParallelState().  Going
 * from DISTRIBUTED to GLOBALIZED state is an expensive operation requiring
 * all-to-all communication.  Using this state can incur a significant
 * performance penalty.
 *
 * <li> The GLOBALIZED state requires more memory.
 *
 * <li> Transitioning from GLOBALIZED state to DISTRIBUTED state is
 * cheap.
 * </ul>
 *
 * @note
 * The general attributes of a BoxLevel are
 * <ul>
 * <li> the set of Box objects with unique BoxIds,
 * <li> the refinement ratio defining their index space, and
 * <li> the parallel state.
 * </ul>
 *
 * Box object uniqueness is based on the Box equality operator,
 * which compares owner MPI ranks and local indices.  Therefore,
 * a valid BoxLevel does not contain two Boxes with the same
 * owner and index.
 */
class BoxLevel
{

public:
   /*!
    * @brief Names of parallel states.
    */
   enum ParallelState { DISTRIBUTED, GLOBALIZED };

   /*!
    * @brief Construct a BoxLevel which will be initialized from the supplied
    * restart database.
    *
    * @param[in] dim
    * @param[in] restart_db
    * @param[in] grid_geom
    */
   BoxLevel(
      const tbox::Dimension& dim,
      tbox::Database& restart_db,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom);

   /*!
    * @brief Copy constructor.
    *
    * New object has the same parallel state as original.
    *
    * Persistent Connectors are not duplicated.  This
    * decision was based on expected usage, which is that
    * copies are either for short term usage or meant to
    * be changed in some way and will invalidate Connectors.
    *
    * @param[in] rhs
    */
   BoxLevel(
      const BoxLevel& rhs);

   /*!
    * @brief Constructs an empty, initialized object.
    *
    * @see addBox()
    * @see initialize()
    *
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   BoxLevel(
      const IntVector& ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Constructs a populated object.
    *
    * @see addBox()
    * @see initialize()
    *
    * @param[in] boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   BoxLevel(
      const BoxContainer& boxes,
      const IntVector& ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~BoxLevel();

   //@{
   //! @name Initialization and clearing methods

   /*!
    * @brief Initialize the BoxLevel
    *
    * The content and state of the object before calling this function
    * is discarded.
    *
    * @see addBox()
    * @see initialize(const BoxContainer&, const IntVector&, const tbox::SAMRAI_MPI&, const ParallelState)
    *
    * @param[in] boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    */
   void
   initialize(
      const BoxContainer& boxes,
      const IntVector& ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Initialize the BoxLevel.
    *
    * Similar to initialize(const BoxContainer&, const IntVector&, const tbox::SAMRAI_MPI&, const ParallelState), except that the @c boxes are mutable.
    *
    * The state of the object before calling this function is
    * discarded.  The Box content before calling this function
    * is returned via the @c boxes argument.
    *
    * @see initializePrivate()
    *
    * @param[in,out] boxes On input, this should contain the
    * Boxes to place in the BoxLevel.  On output, it
    * contains the Boxes that were in the BoxLevel before
    * the call.
    *
    * @param[in] boxes
    * @param[in] ratio
    * @param[in] grid_geom
    * @param[in] mpi
    * @param[in] parallel_state
    *
    * @pre &boxes != &getBoxes()
    */
   void
   swapInitialize(
      BoxContainer& boxes,
      const IntVector& ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Removes non-local boxes; computes bounding box, local number of
    * boxes, local number of cells, max/min box size.  To be called after all
    * Boxes in a BoxLevel have been added to indicate that the BoxLevel is
    * fully defined and ready to be used.
    */
   void
   finalize();

   /*!
    * @brief Returns True if the object has been initialized.
    */
   bool
   isInitialized() const
   {
      return d_ratio(0,0) != 0;
   }

   /*!
    * @brief Remove all the periodic image boxes in the BoxLevel.
    */
   void
   removePeriodicImageBoxes();

   /*!
    * @brief Clear the internal state of the BoxLevel.
    *
    * The BoxLevel will be in an uninitialized state
    * after a call to this method.
    */
   void
   clear();

   /*!
    * @brief Clear the globalized version and the persistent overlap
    * connectors for data consistency.
    *
    * Most of the time, this method is automatically called by methods
    * that know when some data is stale and needs to be cleared.
    * For example, adding a box makes the global number of
    * boxes stale.  However, sometimes it is necessary to call this
    * method manually.  For example, when only some processes add
    * boxes while others do not, resulting in some processes not
    * knowing that the global number of boxes is inconsistent.
    *
    * @param[in] isInvalid A flag indicating that boxes have been (or will
    *            be) removed, thus invalidating the handle.
    */
   void
   clearForBoxChanges(
      bool isInvalid = true)
   {
      if (locked()) {
         TBOX_ERROR("BoxLevel::clearForBoxChanges(): operating on locked BoxLevel."
            << std::endl);
      }
      deallocateGlobalizedVersion();
      clearPersistentOverlapConnectors();
      if (isInvalid) {
         /*
          * Box removal can lead on inconsistent Connectors holding on to
          * handle, so detach the handle.  Box addition does NOT lead to
          * such inconsistencies, so we can leave the handle alone in
          * those cases.
          */
         detachMyHandle();
      }
   }

   //@{

   //! @name Parallelism

   /*!
    * @brief Set the parallel state.
    *
    * This method is potentially expensive.
    * Acquiring remote Box information (when going
    * to GLOBALIZED mode) triggers all-gather communication.
    * More memory is required to store additional Boxes.
    *
    * Data not used by the new state gets deallocated.
    *
    * @param[in] parallel_state
    *
    * @pre isInitialized()
    */
   void
   setParallelState(
      const ParallelState parallel_state);

   /*!
    * @brief Returns the ParallelState of the object.
    */
   ParallelState
   getParallelState() const
   {
      return d_parallel_state;
   }

   /*!
    * @brief If global reduced data (global Box count, global
    * cell count and global bounding box) have not been updated,
    * compute and cache them (communication required).
    *
    * After this method is called, data requiring global reduction can
    * be accessed without further communications, until the object
    * changes.
    *
    * Sets d_global_data_up_to_date to true;
    *
    * @pre isInitialized()
    */
   void
   cacheGlobalReducedData() const;

   /*!
    * @brief Return the globalized version of the BoxLevel,
    * creating it if needed.
    *
    * If the BoxLevel is in globalized state, return @c *this.
    * If not, create and cache a globalized version (if necessary) and
    * return that.
    *
    * The cached version remains until it is removed by
    * deallocateGlobalizedVersion() or a method that can potentially
    * change the Boxes is called.  Note that globalizing and
    * globalized data is not scalable.  Use only when necessary.
    *
    * Obviously, when the globalized version must be created (when the
    * BoxLevel is in DISTRIBUTED state and there is no cached
    * version yet), all processes must make this call at the same
    * point.
    *
    * @pre isInitialized()
    *
    * @post d_globalized_version->getParallelState() == GLOBALIZED
    */
   const BoxLevel&
   getGlobalizedVersion() const;

   /*!
    * @brief Deallocate the internal globalized version of the
    * BoxLevel, if there is any.
    *
    * @pre (d_globalized_version == 0) ||
    *      (d_globalized_version->getParallelState() == GLOBALIZED)
    */
   void
   deallocateGlobalizedVersion() const
   {
      if (d_globalized_version) {
         TBOX_ASSERT(d_globalized_version->getParallelState() == GLOBALIZED);
         d_globalized_version.reset();
      }
   }

   /*!
    * @brief Returns the SAMRAI_MPI communicator over which the Boxes
    * are distributed.
    */
   const tbox::SAMRAI_MPI&
   getMPI() const
   {
      return d_mpi;
   }

   //@}

   /*!
    * @brief Assignment operator duplicates all internal data,
    * including parallel mode.
    *
    * Assignment is a modifying operation, causing the
    * PersistentOverlapConnectors to be cleared.
    *
    * Persistent Connectors are not duplicated.  This
    * decision was based on expected usage, which is
    * that copies are either for short term usage or meant to
    * be changed in some way (thus invalidating current
    * Connectors anyway).
    *
    * @see getPersistentOverlapConnectors()
    *
    * @param[in] rhs
    */
   BoxLevel&
   operator = (
      const BoxLevel& rhs);

   /*!
    * @brief Swap the contents of two BoxLevel objects.
    *
    * Swapping is a modifying operation, so the
    * PersistentOverlapConnectors of the operands are cleared.
    *
    * Persistent Connectors are not swapped.  This decision
    * was based on expected usage, which is that
    * copies are either for short term usage or meant to be
    * changed in some way (thus invalidating current
    * Connectors anyway).
    *
    * @param[in,out] level_a
    * @param[in,out] level_b
    *
    * @pre (&level_a == &level_b) ||
    *      !level_a.isInitialized() || !level_b.isInitialized() ||
    *      (level_a.getDim() == level_b.getDim())
    */
   static void
   swap(
      BoxLevel& level_a,
      BoxLevel& level_b);

   //@}

   /*!
    * @brief Equality comparison.
    *
    * All data required to initialize the object is compared, except
    * for the parallel state.  Thus equality here means just the local
    * parts are equal.  @b BEWARE!  This means that one processor may
    * see the equality differently from another.
    *
    * The cost for the comparison is on the order of the local
    * Box count.  An object may be compared to itself, an
    * efficient operation that always returns true.
    *
    * @param[in] rhs
    */
   bool
   operator == (
      const BoxLevel& rhs) const;

   /*!
    * @brief Inequality comparison.
    *
    * All data required to initialize the object is compared, except
    * for the parallel state.  Thus equality here means just the local
    * parts are equal.  @b BEWARE!  This means that one processor may
    * see the inequality differently from another.
    *
    * The cost for the comparison is on the order of the local
    * Box count.  However, an object may be compared to itself,
    * an efficient operation that always returns false.
    *
    * @param[in] rhs
    */
   bool
   operator != (
      const BoxLevel& rhs) const;

   //@{
   /*!
    * @name Accessors
    */

   /*!
    * @brief Returns the container of local Boxes.
    *
    * @par Important
    * The BoxContainer returned contains periodic image
    * Boxes (if any).  To iterate through real Boxes only, see
    * RealBoxConstIterator.
    *
    * You cannot directly modify the BoxContainer because it may
    * invalidate other internal data.  Use other methods for modifying
    * the BoxContainer.
    *
    * @note It is possible that one may wish to perform repeated searches on
    * the Boxes in the BoxContainer returned by this method.  As noted in
    * BoxContainer's documentation, it may be advantageous to call the makeTree
    * method on the BoxContainer in this case.  If you do, remember that the
    * GridGeometry passed to makeTree must be the GridGeometry held by this
    * BoxLevel.  Thus the proper use of makeTree with the result of a call to
    * this method will look like:
    * @verbatim
    *    box_level.getBoxes().makeTree(box_level.getGridGeometry());
    * @endverbatim
    *
    * @see getGlobalNumberOfBoxes()
    * @see getLocalNumberOfBoxes()
    *
    */
   const BoxContainer&
   getBoxes() const
   {
      return d_boxes;
   }

   /*!
    * @brief Returns the container of global Boxes.
    *
    * @note It is possible that one may wish to perform repeated searches on
    * the Boxes in the BoxContainer returned by this method.  As noted in
    * BoxContainer's documentation, it may be advantageous to call the makeTree
    * method on the BoxContainer in this case.  If you do, remember that the
    * GridGeometry passed to makeTree must be the GridGeometry held by this
    * BoxLevel.  Thus the proper use of makeTree with the result of a call to
    * this method will look like:
    * @verbatim
    *    box_level.getBoxes().makeTree(box_level.getGridGeometry());
    * @endverbatim
    */
   const BoxContainer&
   getGlobalBoxes() const
   {
      return d_global_boxes;
   }

   /*!
    * @brief Fill the container with the global Boxes.
    */
   void
   getGlobalBoxes(
      BoxContainer& global_boxes) const;

   /*
    * TODO: Why are the following two methods here?  Returning local id
    * information like this is dangerous in that it seems to imply that
    * they can be used as an integer range, or a count (last - first + 1).
    * Since the first method is not used and the second is used in a few
    * places, wouldn't it be better to just use the previous method?
    */
   /*!
    * @brief Returns the first LocalId, or one with a value of -1 if
    * no local Box exists.
    *
    * @pre isInitialized()
    */
   LocalId
   getFirstLocalId() const;

   /*!
    * @brief Returns the last LocalId, or one with a value of -1 if no
    * local Box exists.
    *
    * @pre isInitialized()
    */
   LocalId
   getLastLocalId() const;

   /*!
    * @brief Get const access to BoxLevel's refinement ratio
    * (with respect to a reference level).
    */
   const IntVector&
   getRefinementRatio() const
   {
      return d_ratio;
   }

   /*!
    * @brief Return local number of boxes.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   int
   getLocalNumberOfBoxes() const
   {
      TBOX_ASSERT(isInitialized());
      return d_local_number_of_boxes;
   }

   /*!
    * @brief Return number of boxes local to the given rank.
    *
    * Periodic image Boxes are excluded.
    *
    * @param[in] rank
    *
    * @pre isInitialized()
    * @pre (getParallelState() == GLOBALIZED) || (rank == getMPI().getRank())
    * @pre (rank >= 0) && (rank < getMPI().getSize())
    */
   int
   getLocalNumberOfBoxes(
      int rank) const;

   /*!
    * @brief Return global number of Boxes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   int
   getGlobalNumberOfBoxes() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_global_number_of_boxes;
   }

   /*!
    * @brief Return maximum number of Boxes over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   int
   getMaxNumberOfBoxes() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_max_number_of_boxes;
   }

   /*!
    * @brief Return maximum number of Boxes over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   int
   getMinNumberOfBoxes() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_min_number_of_boxes;
   }

   /*!
    * @brief Return local number of cells.
    *
    * Cells in periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   size_t
   getLocalNumberOfCells() const
   {
      TBOX_ASSERT(isInitialized());
      return d_local_number_of_cells;
   }

   /*!
    * @brief Return maximum number of cells over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   size_t
   getMaxNumberOfCells() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_max_number_of_cells;
   }

   /*!
    * @brief Return maximum number of cells over all processes.
    *
    * This requires a global reduction, if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   size_t
   getMinNumberOfCells() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_min_number_of_cells;
   }

   /*!
    * @brief Return number of cells local to the given rank.
    *
    * Cells in periodic image Boxes are excluded.
    *
    * @param[in] rank
    *
    * @pre isInitialized()
    * @pre (getParallelState() == GLOBALIZED) || (rank == getMPI().getRank())
    * @pre (rank >= 0) && (rank < getMPI().getSize())
    */
   size_t
   getLocalNumberOfCells(
      int rank) const;

   /*!
    * @brief Return global number of cells.
    *
    * This requires a global reduction if the global-reduced data has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * Cells in periodic image Boxes are excluded.
    *
    * @pre isInitialized()
    */
   size_t
   getGlobalNumberOfCells() const
   {
      TBOX_ASSERT(isInitialized());
      cacheGlobalReducedData();
      return d_global_number_of_cells;
   }

   /*!
    * @brief Return bounding box for local Boxes in a block.
    */
   const Box&
   getLocalBoundingBox(
      const BlockId& block_id) const
   {
      return d_local_bounding_box[block_id.getBlockValue()];
   }

   /*!
    * @brief Return bounding box for global Boxes in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const Box&
   getGlobalBoundingBox(
      const BlockId& block_id) const
   {
      cacheGlobalReducedData();
      return d_global_bounding_box[block_id.getBlockValue()];
   }

   /*!
    * @brief Return size of the largest local Box in a block.
    */
   const IntVector&
   getLocalMaxBoxSize(
      const BlockId& block_id) const
   {
      return d_local_max_box_size[block_id.getBlockValue()];
   }

   /*!
    * @brief Return size of the smallest local Box in a block.
    */
   const IntVector&
   getLocalMinBoxSize(
      const BlockId& block_id) const
   {
      return d_local_min_box_size[block_id.getBlockValue()];
   }

   /*!
    * @brief Return size of the largest Box globally in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const IntVector&
   getGlobalMaxBoxSize(
      const BlockId& block_id) const
   {
      cacheGlobalReducedData();
      return d_global_max_box_size[block_id.getBlockValue()];
   }

   /*!
    * @brief Return size of the smallest Box globally in a block.
    *
    * This requires a global reduction if the global bounding box has
    * not been computed and cached.  When communication is required,
    * all processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    */
   const IntVector&
   getGlobalMinBoxSize(
      const BlockId& block_id) const
   {
      cacheGlobalReducedData();
      return d_global_min_box_size[block_id.getBlockValue()];
   }

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_ratio.getDim();
   }

   /*!
    * @brief Return the grid geometry associated with this object.
    *
    * If object has never been initialized, return NULL pointer.
    */
   const std::shared_ptr<const BaseGridGeometry>&
   getGridGeometry() const
   {
      return d_grid_geometry;
   }

   //@}

   //@{

   //! @name Methods to modify all Boxes.
   /*
    */
   /*!
    * @brief Refine all Boxes of this BoxLevel by ratio placing result into
    * finer making finer's ratio final_ratio.
    *
    * @param[out] finer
    * @param[in] ratio
    * @param[in] final_ratio
    */
   void
   refineBoxes(
      BoxLevel& finer,
      const IntVector& ratio,
      const IntVector& final_ratio) const
   {
      finer.detachMyHandle();
      finer.d_globalized_version.reset();
      finer.d_boxes = d_boxes;
      finer.d_boxes.refine(ratio);
      finer.d_parallel_state = d_parallel_state;
      if (finer.d_parallel_state == GLOBALIZED) {
         finer.d_global_boxes = d_global_boxes;
         finer.d_global_boxes.refine(ratio);
      }
      finer.d_ratio = final_ratio;
      finer.computeLocalRedundantData();
   }

   /*!
    * @brief Coarsen all Boxes of this BoxLevel by ratio placing result into
    * coarser making coarser's ratio final_ratio.
    *
    * @param[out] coarser
    * @param[in] ratio
    * @param[in] final_ratio
    */
   void
   coarsenBoxes(
      BoxLevel& coarser,
      const IntVector& ratio,
      const IntVector& final_ratio) const
   {
      coarser.detachMyHandle();
      coarser.d_globalized_version.reset();
      coarser.d_boxes = d_boxes;
      coarser.d_boxes.coarsen(ratio);
      coarser.d_parallel_state = d_parallel_state;
      if (coarser.d_parallel_state == GLOBALIZED) {
         coarser.d_global_boxes = d_global_boxes;
         coarser.d_global_boxes.coarsen(ratio);
      }
      coarser.d_ratio = final_ratio;
      coarser.computeLocalRedundantData();
   }

   //@}

   //@{

   //! @name Individual Box methods.

   /*!
    * @brief Create new local Box from given Box and add it to this
    * level.
    *
    * The new Box will be assigned an unused local index, so the input
    * box need not have a valid one.  To be
    * efficient, no communication will be used.  Therefore, the state
    * must be distributed.
    *
    * The new Box will have a periodic shift number
    * corresponding to zero-shift.
    *
    * @note It is imperative that applications which call addBox also call
    * invalidateGlobalData.  It is possible for some processes to add Boxes and
    * for others to not.  Since the addBox method sets d_global_data_up_to_date
    * to false, some processes in this situation will have this flag set to
    * true and others will not.  This will result in a hang in
    * cacheGlobalReducedData when global data is accessed.  We could have
    * performed and allReduce of this flag to the top of cacheGlobalReducedData
    * but this would have added a costly call for every access of global data.
    *
    * @param[in] box
    * @param[in] block_id
    *
    * @return iterator to the new Box
    *
    * @pre getParallelState() == DISTRIBUTED
    * @pre (box.getBlockId() == BlockId::invalidId()) ||
    *      (box.getBlockId() == block_id)
    */
   BoxContainer::const_iterator
   addBox(
      const Box& box,
      const BlockId& block_id);

   /*!
    * @brief Add a Box to this level.
    *
    * @par CAUTION
    * To be efficient, no checks are made to make sure the
    * BoxLevel representation is consistent across all
    * processors.  Setting inconsistent data leads potentially
    * elusive bugs.
    *
    * @par Errors
    * It is an error to add a periodic image of a Box that is
    * not a part of the BoxLevel.
    *
    * It is an error to add any Box that already exists.
    *
    * @note It is imperative that applications which call addBox also call
    * invalidateGlobalData.  It is possible for some processes to add Boxes and
    * for others to not.  Since the addBox method sets d_global_data_up_to_date
    * to false, some processes in this situation will have this flag set to
    * true and others will not.  This will result in a hang in
    * cacheGlobalReducedData when global data is accessed.  We could have
    * performed and allReduce of this flag to the top of cacheGlobalReducedData
    * but this would have added a costly call for every access of global data.
    *
    * @param[in] box
    */
   void
   addBox(
      const Box& box);

   /*!
    * @brief Add a Box to this level without updating summary data such as
    * local number of boxes/cells, bounding box, max/min box size.  Meant to
    * be used during the construction of a BoxLevel as Boxes belonging to the
    * level are found.  finalize() should be called at the end of construction
    * making use of addBoxWithoutUpdate.
    *
    * @param[in] box
    */
   void
   addBoxWithoutUpdate(
      const Box& box)
   {
      if (locked()) {
         TBOX_ERROR("BoxLevel::addBoxWithoutUpdate(): operating on locked BoxLevel."
            << std::endl);
      }
      if (getParallelState() == GLOBALIZED) {
         d_global_boxes.insert(box);
      }
      d_boxes.insert(box);
   }

   /*!
    * @brief Insert given periodic image of an existing Box.
    *
    * Unlike adding a regular Box, it is OK to add a periodic
    * image Box that already exists.  However, that is a no-op.
    *
    * @par CAUTION
    * To be efficient, no checks are made to make sure the
    * BoxLevel representation is consistent across all
    * processors.  Setting inconsistent data leads to potentially
    * elusive bugs.
    *
    * @par Errors
    * It is an error to add a periodic image of a Box that does
    * not exist.
    *
    * TODO: Should we prevent this operation if persistent overlap
    * Connectors are attached to this object?
    *
    * @param[in] existing_box  An existing Box for reference.
    *      This Box must be in the BoxLevel.  The Box added
    *      is an image of the reference Box but shifted to another
    *      position.
    * @param[in] shift_number The valid shift number for the Box being
    *      added.  The shift amount is taken from the PeriodicShiftCatalog.
    *
    * @pre shift_number != getGridGeometry()->getPeriodicShiftCatalog().getZeroShiftNumber()
    */
   void
   addPeriodicBox(
      const Box& existing_box,
      const PeriodicId& shift_number);

   /*!
    * @brief Erase the existing Box specified by its iterator.
    *
    * The given iterator @em MUST be a valid iterator pointing to a
    * Box currently in this object.  After erasing, the iterator
    * is advanced to the next valid Box (or the end of its
    * BoxContainer).
    *
    * Erasing a Box also erases all of its periodic images.
    *
    * TODO: Should we prevent this operation if the object has
    * persistent overlap Connectors?
    *
    * @note It is imperative that applications which call eraseBox also call
    * invalidateGlobalData.  It is possible for some processes to erase Boxes
    * and for others to not.  Since the eraseBox method sets
    * d_global_data_up_to_date to false, some processes in this situation will
    * have this flag set to true and others will not.  This will result in a
    * hang in cacheGlobalReducedData when global data is accessed.  We could
    * have performed and allReduce of this flag to the top of
    * cacheGlobalReducedData but this would have erased a costly call for every
    * access of global data.
    *
    * @param[in] ibox The iterator of the Box to erase.
    *
    * @pre getParallelState() == DISTRIBUTED
    * @pre ibox == getBoxes().find(*ibox)
    */
   void
   eraseBox(
      BoxContainer::iterator& ibox);

   /*!
    * @brief Erase the Box matching the one given.
    *
    * The given Box @em MUST match a Box currently in this
    * object.  Matching means that the BoxId's match
    * (disregarding the spatial coordinates).
    *
    * Erasing a Box also erases all of its periodic images.
    *
    * TODO: Should we prevent this operation if the object has
    * persistent overlap Connectors?
    *
    * @note It is imperative that applications which call eraseBox also call
    * invalidateGlobalData.  It is possible for some processes to erase Boxes
    * and for others to not.  Since the eraseBox method sets
    * d_global_data_up_to_date to false, some processes in this situation will
    * have this flag set to true and others will not.  This will result in a
    * hang in cacheGlobalReducedData when global data is accessed.  We could
    * have performed and allReduce of this flag to the top of
    * cacheGlobalReducedData but this would have erased a costly call for every
    * access of global data.
    *
    * @param[in] box
    *
    * @pre getParallelState() == DISTRIBUTED
    * @pre getBoxes().find(box) != getBoxes().end()
    */
   void
   eraseBox(
      const Box& box);

   /*!
    * @brief Erases the Box matching the supplied Box from this level without
    * updating summary data such as local number of boxes/cells, bounding box,
    * max/min box size.  Meant to be used during the construction of a BoxLevel
    * as Boxes not belonging to the level are found.  finalize() should be
    * called at the end of construction making use of eraseBoxWithoutUpdate.
    *
    * @param[in] box
    */
   void
   eraseBoxWithoutUpdate(
      const Box& box)
   {
      if (locked()) {
         TBOX_ERROR("BoxLevel::eraseBoxWithoutUpdate(): operating on locked BoxLevel."
            << std::endl);
      }
      d_boxes.erase(box);
   }

   /*!
    * @brief Find the Box matching the one given.
    *
    * Only the BoxId matters in matching, so the actual Box can
    * be anything.
    *
    * If @c box is not a local Box, the state must be
    * GLOBALIZED.
    *
    * @param[in] box
    *
    * @return Iterator to the box, or @c
    * getBoxes(owner).end() if box does not exist in set.
    *
    * @pre (box.getOwnerRank() == getMPI().getRank()) ||
    *      (getParallelState() == GLOBALIZED)
    */
   BoxContainer::const_iterator
   getBox(
      const Box& box) const
   {
      if (box.getOwnerRank() == getMPI().getRank()) {
         return d_boxes.find(box);
      } else {
#ifdef DEBUG_CHECK_ASSERTIONS
         if (getParallelState() != GLOBALIZED) {
            TBOX_ERROR(
               "BoxLevel::getBox: cannot get remote box "
               << box << " without being in globalized state." << std::endl);
         }
#endif
         return d_global_boxes.find(box);
      }
   }

   /*!
    * @brief Find the Box specified by the given BoxId and
    * periodic shift.
    *
    * If @c box is not a local Box, the state must be
    * GLOBALIZED.
    * @param[in] box_id
    *
    * @return Iterator to the box, or @c
    * getBoxes(owner).end() if box does not exist in set.
    */
   BoxContainer::const_iterator
   getBox(
      const BoxId& box_id) const
   {
      const Box box(getDim(), box_id);
      return getBox(box);
   }

   /*
    * TODO: What is different about these "strict" methods compared to
    * the preceding ones?  I can't tell from the comments.
    */

   /*!
    * @brief Find the Box matching the one given.
    *
    * Only the BoxId matters in matching, so the actual Box can
    * be anything.
    *
    * If @c box is not owned by the local process, the state
    * must be GLOBALIZED.
    *
    * You cannot directly modify the BoxContainer because it may
    * invalidate other internal data.  Use other methods for modifying
    * the BoxContainer.
    *
    * @param[in] box
    *
    * @return Iterator to the box.
    *
    * @pre ((box.getOwnerRank() == getMPI().getRank()) &&
    *       (getBoxes().find(box) != getBoxes().end())) ||
    *      ((box.getOwnerRank() != getMPI().getRank()) &&
    *       (getParallelState() == GLOBALIZED) &&
    *       (getGlobalBoxes().find(box) != getGlobalBoxes().end()))
    */
   BoxContainer::const_iterator
   getBoxStrict(
      const Box& box) const;

   /*!
    * @brief Find the Box specified by the given BoxId
    *
    * You cannot directly modify the BoxContainer because it may
    * invalidate other internal data.  Use other methods for modifying
    * the BoxContainer.
    *
    * @param[in] box_id
    *
    * @return Iterator to the box.
    *
    * @pre (box_id.getOwnerRank() == getMPI().getRank()) ||
    *      (getParallelState() == GLOBALIZED)
    * @pre box with supplied BoxId exists in the BoxLevel
    */
   BoxContainer::const_iterator
   getBoxStrict(
      const BoxId& box_id) const;

   /*!
    * @brief Find the Box with the given BlockId which is spatially equal to
    * the supplied Box.
    *
    * @param[in] box_to_match
    * @param[in] block_id
    * @param[out] matching_box If there is a box with the supplied BlockId
    * spatially equal to box_to_match then this is set to that Box.
    *
    * @return true if a match is found.
    */
   bool
   getSpatiallyEqualBox(
      const Box& box_to_match,
      const BlockId& block_id,
      Box& matching_box) const;

   /*!
    * @brief Returns true when the object has a Box specified by the
    * BoxId.
    *
    * @param[in] box_id
    */
   bool
   hasBox(
      const BoxId& box_id) const
   {
      const Box box(getDim(), box_id);
      return hasBox(box);
   }

   /*!
    * @brief Returns true when the object has a Box consistent with all
    * of the arguments
    *
    * @param[in] global_id
    * @param[in] periodic_id
    */
   bool
   hasBox(
      const GlobalId& global_id,
      const PeriodicId& periodic_id) const
   {
      const Box box(
         getDim(),
         global_id,
         periodic_id);
      return hasBox(box);
   }

   /*!
    * @brief Returns true when this BoxLevel has a Box matching the
    * BoxId of the given box.
    *
    * @param[in] box
    *
    * @pre (box.getOwnerRank() == getMPI().getRank()) ||
    *      (getParallelState() == GLOBALIZED)
    */
   bool
   hasBox(
      const Box& box) const;

   //@}

   //@{
   /*!
    * @name IO support.
    */

   /*!
    * @brief Write the BoxLevel to a restart database.
    *
    * Write only local parts regardless of parallel state (to avoid
    * writing tons of repetitive data).
    *
    * @param[in,out] restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   //@}

   /*!
    * @brief Sets d_global_data_up_to_date to false.  Must be called after
    * calls to addBox or eraseBox.
    *
    * @see addBox()
    * @see eraseBox()
    */
   void
   invalidateGlobalData()
   {
      if (locked()) {
         TBOX_ERROR("BoxLevel::invalidateGlobalData(): operating on locked BoxLevel."
            << std::endl);
      }
      d_global_data_up_to_date = false;
   }

   /*!
    * @brief Deallocate persistent overlap Connectors, if there are any.
    */
   void
   clearPersistentOverlapConnectors()
   {
      if (d_persistent_overlap_connectors) {
         d_persistent_overlap_connectors->clear();
      }
   }

   /*!
    * @brief Find an overlap Connector with the given head and minimum
    * Connector width.  If the specified Connector is not found, take the
    * specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true) const
   {
      return getPersistentOverlapConnectors().findConnector(head,
         min_connector_width,
         not_found_action,
         exact_width_only);
   }

   /*!
    * @brief Find an overlap Connector with its transpose with the given head
    * and minimum Connector widths.  If the specified Connector is not found,
    * take the specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] transpose_min_connector_width Find the transpose overlap
    *      Connector satisfying this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      const IntVector& transpose_min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true) const
   {
      return getPersistentOverlapConnectors().findConnectorWithTranspose(head,
         min_connector_width,
         transpose_min_connector_width,
         not_found_action,
         exact_width_only);
   }

   /*!
    * @brief Create an overlap Connector, computing relationships by
    * globalizing data.
    *
    * The base will be this BoxLevel.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head This BoxLevel will be the head.
    * @param[in] connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width) const
   {
      return getPersistentOverlapConnectors().createConnector(head,
         connector_width);
   }

   /*!
    * @brief Create an overlap Connector with its transpose, computing
    * relationships by globalizing data.
    *
    * The base will be this BoxLevel.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head This BoxLevel will be the head.
    * @param[in] connector_width
    * @param[in] transpose_connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& connector_width,
      const IntVector& transpose_connector_width) const
   {
      return getPersistentOverlapConnectors().createConnectorWithTranspose(head,
         connector_width,
         transpose_connector_width);
   }

   /*!
    * @brief Cache the supplied overlap Connector and its transpose
    * if it exists.
    *
    * @param[in] connector
    *
    * @pre isInitialized()
    * @pre connector
    */
   void
   cacheConnector(
      std::shared_ptr<Connector>& connector) const
   {
      return getPersistentOverlapConnectors().cacheConnector(connector);
   }

   /*!
    * @brief Returns whether the object has overlap
    * Connectors with the given head and minimum Connector
    * width.
    *
    * TODO:  does the following comment mean that this must be called
    * before the call to findConnector?
    *
    * If this returns true, the Connector fitting the specification
    * exists and findConnector() will not throw an assertion.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    *
    * @return True if a Connector is found, otherwise false.
    */
   bool
   hasConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width) const
   {
      return getPersistentOverlapConnectors().hasConnector(head,
         min_connector_width);
   }

   /*
    * TODO: The following method is "not for general use" and indeed
    * is only used in two Connector classes.  Would it be better to
    * make the method private and make this class a friend of those?
    */

   /*!
    * @brief Get the handle with which Connectors
    * reference the BoxLevel instead of referencing the
    * BoxLevel itself.  Not for general use.
    *
    * Connectors referencing their base and head BoxLevels should
    * reference their handles instead of the BoxLevels themselves.
    * As long as the BoxLevel does not change in a way
    * that can invalidate Connector data, you can access
    * the BoxLevel from the BoxLevelHandle.
    *
    * If the BoxLevel go out of scope before the
    * Connector disconnects, this std::shared_ptr object will
    * stay around until all Connectors have disconnected.
    *
    * Operations that can invalidate Connector data are those
    * that remove information from the BoxLevel.  These
    * are:
    *
    * @li initialize()
    * @li swapInitialize()
    * @li swap()
    * @li clear()
    * @li operator=() (assignment) (Exception: assigning to
    *     self is a no-op, which does not invalidate Connector
    *     data.
    * @li eraseBox() (Note that adding a Box
    *     does not invalidate Connector data.)
    * @li going out of scope
    *
    * @see BoxLevelHandle.
    *
    * @return A std::shared_ptr to the BoxLevelHandle
    *
    * @pre !d_handle || (d_handle->d_box_level == this)
    */
   const std::shared_ptr<BoxLevelHandle>&
   getBoxLevelHandle() const
   {
      if (!d_handle) {
         /*
          * No handle yet.  Generate one attached to this object.
          */
         d_handle.reset(new BoxLevelHandle(this));
      }
      if (d_handle->d_box_level != this) {
         /*
          * Sanity check: The handle for this object should be attached
          * to this object.
          */
         TBOX_ERROR("Library error in BoxLevelHandle::getBoxLevel" << std::endl);
      }
      return d_handle;
   }

   /*!
    * @brief Effectively makes a non-const BoxLevel const.  Prevents any
    * non-const method from executing.
    */
   void
   lock()
   {
      d_locked = true;
   }

   /*!
    * @brief Returns true if the BoxLevel is locked.
    */
   bool
   locked()
   {
      return d_locked;
   }

   //@{

   /*!
    * @name Methods for outputs, error checking and debugging.
    */

   /*!
    * @brief Print Box info from this level
    *
    * @param[in,out] os The output stream
    * @param[in] border
    * @param[in] detail_depth
    */
   void
   recursivePrint(
      std::ostream& os,
      const std::string& border,
      int detail_depth = 2) const;

   /*!
    * @brief A class for outputting BoxLevel.
    *
    * To use, see BoxLevel::format() and BoxLevel::formatStatistics().
    *
    * This class simplifies the insertion of a BoxLevel into a
    * stream while letting the user control how the BoxLevel is
    * formatted for output.
    *
    * Each Outputter is a light-weight object constructed with a
    * BoxLevel and output parameters.  The Outputter is capable
    * of outputting its BoxLevel, formatted according to the
    * parameters.
    */
   class Outputter
   {
      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);
private:
      friend class BoxLevel;
      /*!
       * @brief Copy constructor
       */
      Outputter(
         const Outputter& other);

      /*!
       * @brief Construct the Outputter with a BoxLevel and the
       * parameters needed to output the BoxLevel to a stream.
       *
       * @param[in] box_level
       * @param[in] border
       * @param[in] detail_depth
       * @param[in] output_statistics
       */
      Outputter(
         const BoxLevel& box_level,
         const std::string& border,
         int detail_depth = 2,
         bool output_statistics = false);
      Outputter&
      operator = (
         const Outputter& r);               // Unimplemented private.
      const BoxLevel& d_level;
      const std::string d_border;
      const int d_detail_depth;
      const bool d_output_statistics;
   };

   /*!
    * @brief Return a object that can format the BoxLevel for
    * inserting into output streams.
    *
    * Usage example:
    * @code
    *    tbox::plog << "my box_level:\n" << box_level.format() << endl;
    * @endcode
    *
    * @param[in] border
    * @param[in] detail_depth
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 2) const;

   /*!
    * @brief Return a object that can format the BoxLevel for
    * inserting its global statistics into output streams.
    *
    * Usage example:
    * @code
    *    std::cout << "my box_level statistics:\n"
    *              << box_level.formatStatistics("  ") << std::endl;
    * @endcode
    *
    * @param[in] border
    */
   Outputter
   formatStatistics(
      const std::string& border = std::string()) const;

   //@}

private:
   friend class PersistentOverlapConnectors;

   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_BOX_LEVEL_VERSION;

   /*
    * Static integer constant describing the number of statistics in this
    * class.
    */
   static const int BOX_LEVEL_NUMBER_OF_STATS;

   /*
    * TODO: This same enum is defined in the Connector header.
    * If there is a common use for this, should it be defined in a
    * common location? Also, it seems to be used similarly to the
    * BAD_INTEGER #define in BergerRigoutsosNode.cpp.  Is the intent
    * really the same?
    */
   enum { BAD_INT = (1 << (8 * sizeof(int) - 2)) };

   /*
    * Unimplemented default constructor.
    */
   BoxLevel();

   /*
    * TODO: The comments for the following method use the phrase
    * "local redundant data" three times, but I still don't know
    * what that is!
    */
   /*!
    * @brief Recompute local redundant data.
    *
    * Local redundant data is usually updated immediately after their
    * dependencies change.  On certain occasions, we recompute all
    * the local redundant data.
    */
   void
   computeLocalRedundantData();

   //@{

   /*!
    * @brief Get and store info on remote Boxes.
    *
    * This requires global communication (all-gather).
    * Call acquireRemoteBoxes_pack to pack up messages.
    * Do an all-gather.  Call acquireRemoteBoxes_unpack
    * to unpack data from other processors.
    */
   void
   acquireRemoteBoxes();

   //! @brief Pack local Boxes into an integer array.
   void
   acquireRemoteBoxes_pack(
      std::vector<int>& send_mesg) const;

   /*!
    * @brief Unpack Boxes from an integer array into internal
    * storage.
    */
   void
   acquireRemoteBoxes_unpack(
      const std::vector<int>& recv_mesg,
      std::vector<int>& proc_offset);

   /*!
    * @brief Get and store info on remote Boxes for multiple
    * BoxLevel objects.
    *
    * This method combines communication for the multiple
    * box_levels to increase message passing efficiency.
    *
    * Note: This method is stateless (could be static).
    */
   void acquireRemoteBoxes(
      const int num_sets,
      BoxLevel * multiple_box_level[]);
   //@}

   /*!
    * @brief Get the collection of overlap Connectors dedicated to
    * provide overlap neighbors for this BoxLevel.
    *
    * The PersistentOverlapConnectors provides overlap neighbors for
    * this BoxLevel.  Its role is to create and manage
    * persistent overlap Connectors based at this BoxLevel and
    * persisting until the BoxLevel changes (so they should not
    * be set up until the BoxLevel is in its final state).  This
    * is the mechanism by which code that can efficiently generate the
    * overlap Connectors (usually the code that generated the
    * BoxLevel) provides overlap data to code using the
    * BoxLevel.  The PersistentOverlapConnectors are guaranteed
    * to be correct, so any changes to the BoxLevel will cause
    * current Connectors to be deallocated.
    *
    * @see PersistentOverlapConnectors for instructions on creating
    * the Connectors.
    */
   PersistentOverlapConnectors&
   getPersistentOverlapConnectors() const;

   /*!
    * @brief Detach this object from the handle it has been using.
    *
    * Postcondition: Objects that cached the handle would no longer
    * be able to access this BoxLevel by the handle.
    */
   void
   detachMyHandle()
   {
      if (d_handle) {
         clearPersistentOverlapConnectors();
         d_handle->detachMyBoxLevel();
         d_handle.reset();
      }
   }

   /*!
    * @brief Encapsulates functionality common to all initialization
    * functions.
    *
    * @pre getDim() == ratio.getDim()
    */
   void
   initializePrivate(
      const IntVector& ratio,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI::getSAMRAIWorld(),
      const ParallelState parallel_state = DISTRIBUTED);

   /*!
    * @brief Read the BoxLevel from a restart database.
    *
    * Put the BoxLevel in the DISTRIBUTED parallel state and
    * read only local parts.
    *
    * If the BoxLevel is initialized, use its SAMRAI_MPI object
    * and require its refinement ratio to match that in the database.
    * If the BoxLevel is uninitialized, it will be initialized
    * to use tbox::SAMRAI_MPI::getSAMRAIWorld() for the SAMRAI_MPI
    * object.  Note that these behaviors have not been extensively
    * discussed by the SAMRAI developers and may be subject to change.
    *
    * @param[in,out] restart_db
    * @param[in] grid_geom
    */
   void
   getFromRestart(
      tbox::Database& restart_db,
      const std::shared_ptr<const BaseGridGeometry>& grid_geom);

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      t_initialize_private = tbox::TimerManager::getManager()->
         getTimer("hier::BoxLevel::initializePrivate()");
      t_acquire_remote_boxes = tbox::TimerManager::getManager()->
         getTimer("hier::BoxLevel::acquireRemoteBoxes()");
      t_cache_global_reduced_data = tbox::TimerManager::getManager()->
         getTimer("hier::BoxLevel::cacheGlobalReducedData()");
   }

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
      t_initialize_private.reset();
      t_acquire_remote_boxes.reset();
      t_cache_global_reduced_data.reset();
   }

   /*!
    * @brief BoxLevel is a parallel object,
    * and this describes its MPI object.
    */
   tbox::SAMRAI_MPI d_mpi;

   /*!
    * @brief Locally-stored Boxes.
    *
    * This is always the container of local Boxes, regardless of
    * parallel mode.
    */
   BoxContainer d_boxes;

   /*!
    * @brief Locally-stored global Boxes (for GLOBALIZED mode).
    *
    * In DISTRIBUTED mode, this is empty.
    */
   BoxContainer d_global_boxes;

   /*
    * TODO: I certainly hope we are not using tests on whether the
    * ratio vector is zero to check whether we have an initialized object.
    */
   /*!
    * @brief Refinement ratio from a reference such as level 0.
    *
    * If d_ratio(0) == 0, the object is in uninitialized state.
    */
   IntVector d_ratio;

   /*!
    * @brief Local cell count, excluding periodic images.
    *
    * Unlike d_global_number_of_cells, this parameter is always current.
    */
   size_t d_local_number_of_cells;

   /*!
    * @brief Global cell count, excluding periodic images.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    */
   mutable size_t d_global_number_of_cells;

   /*!
    * @brief Local Box count, excluding periodic images.
    *
    * Unlike d_global_number_of_boxes, this parameter is always current.
    */
   int d_local_number_of_boxes;

   /*!
    * @brief Global box count, excluding periodic images.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    */
   mutable int d_global_number_of_boxes;

   //! @brief Global max box count on any proc, excluding periodic images.
   mutable int d_max_number_of_boxes;
   //! @brief Global min box count on any proc, excluding periodic images.
   mutable int d_min_number_of_boxes;
   //! @brief Global max cell count on any proc, excluding periodic images.
   mutable size_t d_max_number_of_cells;
   //! @brief Global min cell count on any proc, excluding periodic images.
   mutable size_t d_min_number_of_cells;

   //! @brief Max size of largest local box, one for each block.
   std::vector<IntVector> d_local_max_box_size;
   //! @brief Max size of largest box globally, one for each block.
   mutable std::vector<IntVector> d_global_max_box_size;
   //! @brief Min size of largest local box, one for each block.
   std::vector<IntVector> d_local_min_box_size;
   //! @brief Min size of largest box globally, one for each block.
   mutable std::vector<IntVector> d_global_min_box_size;

   /*!
    * @brief Bounding box of local Boxes, excluding periodic images.
    * One for each block.
    */
   std::vector<Box> d_local_bounding_box;

   /*!
    * @brief Whether d_local_bounding_box is up to date (or needs
    * recomputing.
    */
   bool d_local_bounding_box_up_to_date;

   /*!
    * @brief Bounding box of global Boxes, excluding periodic images.
    * One for each block.
    *
    * This is mutable because it depends on the Boxes and may be
    * saved by a const object if computed.
    */
   mutable std::vector<Box> d_global_bounding_box;

   /*!
    * @brief Whether globally reduced data is up to date or needs
    * recomputing using cacheGlobalReducedData().
    */
   mutable bool d_global_data_up_to_date;

   /*!
    * @brief State flag.
    *
    * Modified by setParallelState().
    */
   ParallelState d_parallel_state;

   /*!
    * @brief A globalized version of the BoxLevel.
    *
    * Initialized by getGlobalizedVersion().  Deallocated by
    * deallocateGlobalizedVersion().
    *
    * Like other redundant data, this is automatically removed if any
    * method that can potentially change the BoxLevel is called.
    *
    * This is mutable because it is redundant data and gets
    * automatically set as needed.
    */
   mutable std::unique_ptr<BoxLevel const> d_globalized_version;

   /*!
    * @brief Connectors managed by this BoxLevel,
    * providing overlap neighbor data across multiple
    * scopes.
    *
    * This is mutable so it can be allocated as needed (by
    * getPersistentOverlapConnectors()).  We can make it non-mutable
    * by always allocating the PersistentOverlapConnectors in the
    * constructor, but most BoxLevel won't need it at all.
    */
   mutable std::unique_ptr<PersistentOverlapConnectors> d_persistent_overlap_connectors;

   /*!
    * @brief A Handle for Connectors to reference this
    * BoxLevel, used to help prevent invalid Connector
    * data.
    *
    * Connectors reference the handle instead of the
    * BoxLevel directly.  When the BoxLevel
    * changes in a way that can invalidate Connector data,
    * it detaches its handle from itself.  A detached handle
    * tells Connectors that the BoxLevel has changed
    * in a way that can invalidate their data.
    *
    * Note: The automatic detaching mechanism prevents some
    * logic errors.  It cannot prevent incorrect Connector
    * data because correctness depends on the Connector's
    * intended usage.
    */
   mutable std::shared_ptr<BoxLevelHandle> d_handle;

   static std::shared_ptr<tbox::Timer> t_initialize_private;
   static std::shared_ptr<tbox::Timer> t_acquire_remote_boxes;
   static std::shared_ptr<tbox::Timer> t_cache_global_reduced_data;

   /*!
    * @brief std::shared_ptr to the grid geometry associated with this
    * object.
    */
   std::shared_ptr<const BaseGridGeometry> d_grid_geometry;

   bool d_locked;

   /*!
    * @brief A LocalId object with value of -1.
    */
   static const LocalId s_negative_one_local_id;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif  // included_hier_BoxLevel
