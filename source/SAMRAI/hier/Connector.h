/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of distributed box-graph relationships from one BoxLevel
 *                to another.
 *
 ************************************************************************/
#ifndef included_hier_Connector
#define included_hier_Connector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelHandle.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxNeighborhoodCollection.h"
#include "SAMRAI/tbox/Timer.h"

#include <set>
#include <string>
#include <vector>

namespace SAMRAI {
namespace hier {

class BoxLevelHandle;

/*!
 * @brief A container which holds relationship connections between two
 * BoxLevels.
 *
 * Connectors have a notion of a "base" and a "head", representing a
 * directional relationship between two BoxLevels.  The relationships
 * are a collection of Boxes in the head, pointed to by a Box in
 * the base.   The association between base and head relationships is
 * 1 .. 0-many.  That is, one Box in the base can be related to zero or
 * more Boxes (called its NeighborSet) in the head.
 *
 * @par Usage
 * Connections in a Connector can have three possible relationships:
 *
 * # A Box in the base has no related NeighborSet.  In this case,
 *   the Box in the base will exist as-is in the head.
 * # A Box in the base has NeighborSet which is empty.  In this case, the
 *   Box from the base will not exist in the head BoxLevel.
 * # A Box in the base has a corresponding NeighborSet which in
 *   non-empty.  In this case, the NeighborSet contains the set of Boxes
 *   to which the Box in the base is related.
 */

class Connector
{
public:
   /*!
    * @brief NeighborsSet is a clarifying typedef.
    */
   typedef BoxContainer NeighborSet;

   /*!
    * @brief Type of the iterator over neighborhoods.
    */
   typedef BoxNeighborhoodCollection::ConstIterator ConstNeighborhoodIterator;

   /*!
    * @brief Type of the iterator over neighborhoods.
    */
   typedef BoxNeighborhoodCollection::Iterator NeighborhoodIterator;

   /*!
    * @brief Type of the iterator over neighbors in a neighborhood.
    */
   typedef BoxNeighborhoodCollection::ConstNeighborIterator ConstNeighborIterator;

   /*!
    * @brief Type of the iterator over neighbors in a neighborhood.
    */
   typedef BoxNeighborhoodCollection::NeighborIterator NeighborIterator;

   /*!
    * @brief Creates an uninitialized Connector object in the
    * distributed state.
    *
    * @param dim The dimension of the head and base BoxLevels that
    * this object will eventually connect.
    *
    * @see setBase()
    * @see setHead()
    * @see setWidth()
    */
   explicit Connector(
      const tbox::Dimension& dim);

   /*!
    * @brief Creates a Connector which is initialized from a restart database.
    *
    * @param dim The dimension of the head and base BoxLevels that
    * this object will eventually connect.
    *
    * @param restart_db Restart Database written by a Connector.
    *
    * @see setBase()
    * @see setHead()
    * @see setWidth()
    */
   Connector(
      const tbox::Dimension& dim,
      tbox::Database& restart_db);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   Connector(
      const Connector& other);

   /*!
    * @brief Initialize a Connector with no defined relationships.
    *
    * The Connector's relationships are initialized to a dummy state.
    *
    * @param[in] base_box_level
    * @param[in] head_box_level
    * @param[in] base_width
    * @param[in] parallel_state
    *
    * @pre (base_box_level.getDim() == head_box_level.getDim()) &&
    *      (base_box_level.getDim() == base_width.getDim())
    */
   Connector(
      const BoxLevel& base_box_level,
      const BoxLevel& head_box_level,
      const IntVector& base_width,
      const BoxLevel::ParallelState parallel_state = BoxLevel::DISTRIBUTED);

   /*!
    * @brief Destructor.
    */
   virtual ~Connector();

   /*!
    * @brief Set this to the transpose of another Connector and
    * populate its edges with the other's transposed edges.
    *
    * This method uses communication to acquire the transpose edges.
    *
    * @param other [i]
    *
    * @param mpi SAMRAI_MPI to use for communication.  If omitted, use
    * the other.getBase().getMPI() by default.  If specified, must be
    * congruent with the default.
    */
   void
   computeTransposeOf(
      const Connector& other,
      const tbox::SAMRAI_MPI& mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL));

   /*!
    * @brief Transpose the visible relationships so that they point from
    * each visible head box to a set of local base boxes.
    */
   void
   reorderRelationshipsByHead(
      std::map<Box, BoxContainer, Box::id_less>& relationships_by_head) const;

   /*!
    * @brief Clear the Connector, putting it into an uninitialized state.
    */
   void
   clear()
   {
      if (d_base_handle) {
         d_relationships.clear();
         d_global_relationships.clear();
         d_mpi.setCommunicator(MPI_COMM_NULL);
         d_base_handle.reset();
         d_head_handle.reset();
         d_parallel_state = BoxLevel::DISTRIBUTED;
      }
   }

   /*!
    * @brief Clear the Connector's neighborhood relations.
    */
   void
   clearNeighborhoods()
   {
      d_relationships.clear();
      d_global_relationships.clear();
   }

   /*!
    * @brief Returns true if the object has been finalized
    */
   bool
   isFinalized() const
   {
      return d_finalized;
   }

   /*!
    * @brief Iterator pointing to the first neighborhood.
    */
   ConstNeighborhoodIterator
   begin() const
   {
      return d_relationships.begin();
   }

   /*!
    * @brief Iterator pointing to the first neighborhood.
    */
   NeighborhoodIterator
   begin()
   {
      return d_relationships.begin();
   }

   /*!
    * @brief Iterator pointing one past the last neighborhood.
    */
   ConstNeighborhoodIterator
   end() const
   {
      return d_relationships.end();
   }

   /*!
    * @brief Iterator pointing one past the last neighborhood.
    */
   NeighborhoodIterator
   end()
   {
      return d_relationships.end();
   }

   /*!
    * @brief Iterator pointing to the first neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   ConstNeighborIterator
   begin(
      const ConstNeighborhoodIterator& nbrhd) const
   {
      return nbrhd.d_collection->begin(nbrhd);
   }

   /*!
    * @brief Iterator pointing to the first neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   NeighborIterator
   begin(
      NeighborhoodIterator& nbrhd)
   {
      BoxNeighborhoodCollection* tmp =
         const_cast<BoxNeighborhoodCollection *>(nbrhd.d_collection);
      return tmp->begin(nbrhd);
   }

   /*!
    * @brief Iterator pointing one past the last neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   ConstNeighborIterator
   end(
      const ConstNeighborhoodIterator& nbrhd) const
   {
      return nbrhd.d_collection->end(nbrhd);
   }

   /*!
    * @brief Iterator pointing one past the last neighbor in nbrhd.
    *
    * @param nbrhd The neighborhood whose neighbors are to be iterated.
    */
   NeighborIterator
   end(
      NeighborhoodIterator& nbrhd)
   {
      BoxNeighborhoodCollection* tmp =
         const_cast<BoxNeighborhoodCollection *>(nbrhd.d_collection);
      return tmp->end(nbrhd);
   }

   /*!
    * @brief Returns an Iterator pointing to the neighborhood of box_id--
    * localized version.
    *
    * @param[in] box_id
    */
   ConstNeighborhoodIterator
   findLocal(
      const BoxId& box_id) const
   {
      BoxId non_per_id(box_id.getGlobalId(),
                       PeriodicId::zero());
      return d_relationships.find(non_per_id);
   }

   /*!
    * @brief Returns an Iterator pointing to the neighborhood of box_id--
    * localized version.
    *
    * @param[in] box_id
    */
   NeighborhoodIterator
   findLocal(
      const BoxId& box_id)
   {
      BoxId non_per_id(box_id.getGlobalId(),
                       PeriodicId::zero());
      return d_relationships.find(non_per_id);
   }

   /*!
    * @brief Returns an Iterator pointing to the neighborhood of box_id--
    * globalized version.
    *
    * @param[in] box_id
    */
   ConstNeighborhoodIterator
   find(
      const BoxId& box_id) const
   {
      const BoxNeighborhoodCollection& relationships = getRelations(box_id);
      BoxId non_per_id(box_id.getGlobalId(), PeriodicId::zero());
      ConstNeighborhoodIterator ei = relationships.find(non_per_id);
      if (ei == relationships.end()) {
         TBOX_ERROR("Connector::find: No neighbor set exists for\n"
            << "box " << box_id << ".\n");
      }
      return ei;
   }

   /*!
    * @brief Returns true if the local neighborhoods of this and other are the
    * same.
    *
    * @param[in] other
    */
   bool
   localNeighborhoodsEqual(
      const Connector& other) const
   {
      return d_relationships == other.d_relationships;
   }

   /*!
    * @brief Returns true if the neighborhood of the supplied BoxId of this
    * and other are the same.
    *
    * @param[in] box_id
    * @param[in] other
    */
   bool
   neighborhoodEqual(
      const BoxId& box_id,
      const Connector& other) const
   {
      const BoxNeighborhoodCollection& relationships = getRelations(box_id);
      const BoxNeighborhoodCollection& other_relationships =
         other.getRelations(box_id);
      BoxId non_per_id(box_id.getGlobalId(),
                       PeriodicId::zero());
      return relationships.neighborhoodEqual(box_id, other_relationships);
   }

   /*!
    * @brief Return true if a neighbor set exists for the specified
    * BoxId.
    *
    * @param[in] box_id
    */
   bool
   hasNeighborSet(
      const BoxId& box_id) const
   {
      const BoxNeighborhoodCollection& relationships = getRelations(box_id);
      BoxId non_per_id(box_id.getGlobalId(), PeriodicId::zero());
      ConstNeighborhoodIterator ei = relationships.find(non_per_id);
      return ei != relationships.end();
   }

   /*!
    * @brief Return true if the supplied box is in the neighborhood of the
    * supplied BoxId.
    *
    * @param[in] box_id
    * @param[in] neighbor
    *
    * @pre box_id.getOwnerRank() == getMPI().getRank()
    */
   bool
   hasLocalNeighbor(
      const BoxId& box_id,
      const Box& neighbor) const
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
      return d_relationships.hasNeighbor(box_id, neighbor);
   }

   /*!
    * @brief Return the neighbor set for the specified BoxId.
    *
    * @param[in] box_id
    * @param[out] nbr_boxes
    *
    * @pre hasNeighborSet(box_id)
    */
   void
   getNeighborBoxes(
      const BoxId& box_id,
      BoxContainer& nbr_boxes) const
   {
      TBOX_ASSERT(hasNeighborSet(box_id));

      const BoxNeighborhoodCollection& relationships = getRelations(box_id);
      BoxId non_per_id(box_id.getGlobalId(),
                       PeriodicId::zero());
      relationships.getNeighbors(non_per_id, nbr_boxes);
   }

   /*!
    * @brief Return all neighbors for all neighborhoods.
    *
    * @param[out] neighbors
    */
   void
   getLocalNeighbors(
      BoxContainer& neighbors) const
   {
      d_relationships.getNeighbors(neighbors);
   }

   /*!
    * @brief Return all neighbors for all neighborhoods segragated by BlockId.
    *
    * @param[out] neighbors
    */
   void
   getLocalNeighbors(
      std::map<BlockId, BoxContainer>& neighbors) const
   {
      d_relationships.getNeighbors(neighbors);
   }

   /*!
    * @brief Returns the number of neighbors in the neighborhood with the
    * supplied BoxId.
    *
    * @param[in] box_id
    *
    * @pre hasNeighborSet(box_id)
    */
   int
   numLocalNeighbors(
      const BoxId& box_id) const
   {
      TBOX_ASSERT(hasNeighborSet(box_id));
      BoxId non_per_id(box_id.getGlobalId(),
                       PeriodicId::zero());
      return d_relationships.numNeighbors(non_per_id);
   }

   /*!
    * @brief Returns the number of empty neighborhoods in the Connector.
    *
    * @return The number of empty neighborhoods in the Connector.
    */
   int
   numLocalEmptyNeighborhoods() const
   {
      int ct = 0;
      for (ConstNeighborhoodIterator itr = begin(); itr != end(); ++itr) {
         if (d_relationships.emptyBoxNeighborhood(itr)) {
            ++ct;
         }
      }
      return ct;
   }

   /*!
    * @brief Places the ranks of the processors owning all neighbors into
    * owners.
    *
    * @param[out] owners
    */
   void
   getLocalOwners(
      std::set<int>& owners) const
   {
      d_relationships.getOwners(owners);
   }

   /*!
    * @brief Places the ranks of the processors owning the neighbors of the Box
    * pointed to by base_boxes_itr into owners.
    *
    * @param[in] base_boxes_itr
    * @param[out] owners
    */
   void
   getLocalOwners(
      ConstNeighborhoodIterator& base_boxes_itr,
      std::set<int>& owners) const
   {
      d_relationships.getOwners(base_boxes_itr, owners);
   }

   //@{
   /*!
    * @name Algorithms for changing individual Box's neighbor data
    */

   /*!
    * @brief Insert additional neighbors for the specified base Box.
    *
    * @param[in] neighbors
    * @param[in] base_box
    *
    * @pre (getParallelState() == BoxLevel::GLOBALIZED) ||
    *      (base_box.getOwnerRank() == getMPI().getRank())
    * @pre getBase().hasBox(base_box)
    */
   void
   insertNeighbors(
      const BoxContainer& neighbors,
      const BoxId& base_box);

   /*!
    * @brief Erase neighbor of the specified BoxId.
    *
    * @param[in] neighbor
    * @param[in] box_id
    *
    * @pre (getParallelState() == BoxLevel::GLOBALIZED) ||
    *      (base_box.getOwnerRank() == getMPI().getRank())
    * @pre getBase().hasBox(box_id)
    */
   void
   eraseNeighbor(
      const Box& neighbor,
      const BoxId& box_id);

   /*!
    * @brief Adds a neighbor of the specified BoxId.
    *
    * @param[in] neighbor
    * @param[in] box_id
    *
    * @pre box_id.getOwnerRank() == getMPI().getRank()
    */
   void
   insertLocalNeighbor(
      const Box& neighbor,
      const BoxId& box_id)
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
      d_relationships.insert(box_id, neighbor);
   }

   /*!
    * @brief Adds a neighbor of the base box pointed to by base_box_itr.
    *
    * @param[in] neighbor
    * @param[in] base_box_itr
    *
    * @pre base_box_itr->getOwnerRank() == getMPI().getRank()
    */
   void
   insertLocalNeighbor(
      const Box& neighbor,
      NeighborhoodIterator& base_box_itr)
   {
      TBOX_ASSERT(base_box_itr->getOwnerRank() == d_mpi.getRank());
      d_relationships.insert(base_box_itr, neighbor);
   }

   /*!
    * @brief Erases the neighborhood of the specified BoxId.
    *
    * @param[in] box_id
    *
    * @pre box_id.getOwnerRank() == getMPI().getRank()
    */
   void
   eraseLocalNeighborhood(
      const BoxId& box_id)
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
      d_relationships.erase(box_id);
   }

   /*!
    * @brief Remove all the periodic relationships in the Connector.
    */
   void
   removePeriodicRelationships()
   {
      d_relationships.erasePeriodicNeighbors();
      if (d_parallel_state == BoxLevel::GLOBALIZED) {
         d_global_relationships.erasePeriodicNeighbors();
      }
   }

   /*!
    * @brief Remove all the periodic neighbors in all local neighborhoods.
    */
   void
   removePeriodicLocalNeighbors()
   {
      d_relationships.erasePeriodicNeighbors();
   }

   /*!
    * @brief Check for any base boxes which are periodic.
    *
    * @return true if any base box is periodic
    */
   bool
   hasPeriodicLocalNeighborhoodBaseBoxes() const
   {
      bool result = false;
      for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         if (ei->getPeriodicId().getPeriodicValue() != 0) {
            result = true;
            break;
         }
      }
      return result;
   }

   /*!
    * @brief Make an empty set of neighbors of the supplied box_id.
    *
    * @param[in] box_id
    *
    * @pre box_id.getOwnerRank() == getMPI().getRank()
    */
   NeighborhoodIterator
   makeEmptyLocalNeighborhood(
      const BoxId& box_id)
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
      return d_relationships.insert(box_id).first;
   }

   /*!
    * @brief Remove empty sets of neighbors.
    */
   void
   eraseEmptyNeighborSets()
   {
      d_relationships.eraseEmptyNeighborhoods();
      d_global_data_up_to_date = false;
   }

   /*!
    * @brief Returns true is the neighborhood of the supplied BoxId is empty.
    *
    * @param[in] box_id
    */
   bool
   isEmptyNeighborhood(
      const BoxId& box_id) const
   {
      return getRelations(box_id).emptyBoxNeighborhood(box_id);
   }

   /*!
    * @brief Coarsen all neighbors of this connector by ratio.
    *
    * @param[in] ratio
    */
   void
   coarsenLocalNeighbors(
      const IntVector& ratio)
   {
      d_relationships.coarsenNeighbors(ratio);
   }

   /*!
    * @brief Refine all neighbors of this connector by ratio.
    *
    * @param[in] ratio
    */
   void
   refineLocalNeighbors(
      const IntVector& ratio)
   {
      d_relationships.refineNeighbors(ratio);
   }

   /*!
    * @brief Grow all neighbors of this connector by growth.
    *
    * @param[in] growth
    */
   void
   growLocalNeighbors(
      const IntVector& growth)
   {
      d_relationships.growNeighbors(growth);
   }

   //@}

   /*!
    * @brief Enforces implicit class invariants and removes non-local neighbors
    * from relationships.
    *
    * To be called after modifying a Connector's context through setBase,
    * setHead, or setWidth methods.
    *
    * @pre d_base_handle && d_head_handle
    * @pre getBase().getGridGeometry() == getHead().getGridGeometry()
    * @pre (getBase().getRefinementRatio() >= getHead().getRefinementRatio()) ||
    *      (getBase().getRefinementRatio() <= getHead().getRefinementRatio())
    * @pre (getParallelState() == BoxLevel::DISTRIBUTED) ||
    *      (getBase().getParallelState() == BoxLevel::GLOBALIZED)
    */
   void
   finalizeContext();

   /*!
    * @brief Change the Connector base to new_base.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_base
    * @param finalize_context
    *
    * @pre new_base.isInitialized()
    */
   void
   setBase(
      const BoxLevel& new_base,
      bool finalize_context = false);

   /*!
    * @brief Return a reference to the base BoxLevel.
    *
    * @pre isFinalized()
    */
   const BoxLevel&
   getBase() const
   {
      TBOX_ASSERT(isFinalized());
      return d_base_handle->getBoxLevel();
   }

   /*!
    * @brief Change the Connector head to new_head.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_head
    * @param finalize_context
    *
    * @pre new_head.isInitialized()
    */
   void
   setHead(
      const BoxLevel& new_head,
      bool finalize_context = false);

   /*!
    * @brief Return a reference to the head BoxLevel.
    *
    * @pre isFinalized()
    */
   const BoxLevel&
   getHead() const
   {
      TBOX_ASSERT(isFinalized());
      return d_head_handle->getBoxLevel();
   }

   /*!
    * @brief Get the refinement ratio between the base and head
    * BoxLevels.
    *
    * The ratio is the same regardless of which is the coarser of the two.
    * Use getHeadCoarserFlag() to determine which is coarser.  If the ratio
    * cannot be represented by an IntVector, truncated.  @see ratioIsExact().
    *
    * @pre isFinalized()
    */
   const IntVector&
   getRatio() const
   {
      TBOX_ASSERT(isFinalized());
      return d_ratio;
   }

   /*!
    * @brief Whether the ratio given by getRatio() is exact.
    *
    * The ratio is exact if it can be represented by an IntVector.
    * @see getRatio().
    *
    * @pre isFinalized()
    */
   bool
   ratioIsExact() const
   {
      TBOX_ASSERT(isFinalized());
      return d_ratio_is_exact;
   }

   /*!
    * @brief Return true if head BoxLevel is coarser than base
    * BoxLevel.
    *
    * @pre isFinalized()
    */
   bool
   getHeadCoarserFlag() const
   {
      TBOX_ASSERT(isFinalized());
      return d_head_coarser;
   }

   /*!
    * @brief Return true if the Connector contains only relationships to local
    * Boxes.
    *
    * The check only applies to neighbors of local base Boxes,
    * so it is possible for the results to be different on different
    * processors.
    */
   bool
   isLocal() const
   {
      return d_relationships.isLocal(getMPI().getRank());
   }

   /*!
    * @brief Create and return this Connector's transpose, assuming that all
    * relationships are local (no remote neighbors).
    *
    * If any remote neighbor is found an unrecoverable assertion is
    * thrown.
    *
    * Non-periodic relationships in are simply reversed to get the transpose
    * relationship.  For each periodic relationship we create a periodic
    * relationship incident from the unshifted head neighbor to the shifted
    * base neighbor.  This is because all relationships must be incident from a
    * real (unshifted) Box.
    */
   virtual Connector *
   createLocalTranspose() const;

   /*!
    * @brief Create and return this Connector's transpose.
    *
    * Similar to createLocalTranspose(), but this method allows
    * non-local edges.  Global data is required, so this method
    * is not scalable.
    */
   virtual Connector *
   createTranspose() const;

   /*!
    * @brief Assignment operator
    */
   Connector&
   operator = (
      const Connector& rhs);

   //  TODO:  need to find out what the use case is for this, especially
   //  considering the caution statement.
   /*!
    * @brief Equality operator checks relationship data, Connector width and
    * equality of base and head Box pointers.
    *
    * @par CAUTION
    * Equality here means just the local parts are equal.
    * This means that one processor may see the equality differently
    * from another.
    *
    * The cost for the comparison is on the order of the local relationship
    * count.  However, an object may be compared to itself, an
    * efficient operation that always returns true.  When comparing
    * Connector objects, if you expect equality to hold, using the
    * same objects would improve performance.
    */
   bool
   operator == (
      const Connector& rhs) const;

   /*!
    * @brief Inequality operator checks the same data that equality
    * operator checks.
    *
    * @see operator==( const Connector &rhs );
    */
   bool
   operator != (
      const Connector& rhs) const;

   /*!
    * @brief Set the parallel distribution state.
    *
    * Before a Connector can be in a GLOBALIZED state, the base
    * BoxLevel given in setBase() must already be in
    * GLOBALIZED mode.  The base BoxLevel should remain in
    * GLOBALIZED mode for compatibility with the Connector.
    *
    * This method is not necessarily trivial.  More memory is required
    * to store additional relationships.
    *
    * For serial (one processor) runs, there is no difference between
    * the parallel states (except for the names), and there is no real
    * cost for switching parallel states.
    *
    * @param[in] parallel_state
    *
    * @pre isFinalized()
    */
   void
   setParallelState(
      const BoxLevel::ParallelState parallel_state);

   /*!
    * @brief Return the current parallel state.
    */
   BoxLevel::ParallelState
   getParallelState() const
   {
      return d_parallel_state;
   }

   /*!
    * @brief Returns the MPI communication object, which is always
    * that of the base BoxLevel.
    *
    * @pre isFinalized()
    */
   const tbox::SAMRAI_MPI&
   getMPI() const
   {
      TBOX_ASSERT(isFinalized());
      return getBase().getMPI();
   }

   /*!
    * @brief Change the Connector width to new_width.  If finalize_context is
    * true then this is the last atomic change being made and finalizeContext
    * should be called.
    *
    * @param new_width
    * @param finalize_context
    *
    * @pre new_width >= IntVector::getZero(new_width.getDim())
    */
   void
   setWidth(
      const IntVector& new_width,
      bool finalize_context = false);

   /*!
    * @brief Return the Connector width associated with the relationships.
    *
    * For overlap Connectors, an relationship exists between a base and head
    * Boxes if the base box, grown by this width,
    * overlaps the head box.  For mapping Connectors, the width
    * the amount that a pre-map box must grow to nest the post-map
    * boxes
    *
    * @pre isFinalized().
    */
   const IntVector&
   getConnectorWidth() const
   {
      TBOX_ASSERT(isFinalized());
      return d_base_width;
   }

   /*!
    * @brief Shrink the width of the connector modifying the proximity
    * relationships as needed.
    *
    * @param[in] new_width
    *
    * @pre new_width <= getConnectorWidth()
    * @pre getParallelState() == BoxLevel::DISTRIBUTED
    */
   void
   shrinkWidth(
      const IntVector& new_width);

   //@{
   /*!
    * @name For outputs, error checking and debugging.
    */

   /*
    * @brief output data
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
    * @brief Return true if two Connector objects are
    * transposes of each other.
    *
    * Each Connector represents a set of directed relationships incident from
    * its base to its head.  The transpose represent the relationships in the
    * opposite direction.  In order for two Connector objects to be transpose
    * of each other, their Connector widths and base refinement ratios must be
    * such that an relationship in one set also appears in the other set.
    * A transpose set must have
    * @li base and head BoxLevels reversed from the untransposed set.
    * @li the same Connector width, although it is described in the index
    *     space of a different base BoxLevel.
    *
    * @param[in] other
    */
   bool
   isTransposeOf(
      const Connector& other) const;

   /*!
    * @brief Given the Connector width in the head index space, convert
    * it to the base index space.
    *
    * This method is useful for computing Connector widths for
    * transpose Connectors.  It handles negative refinement ratios. By
    * SAMRAI convention, a refinement ratio of -N is interpreted as
    * 1/N.)
    *
    * This method is static because (1) it has nothing to do with an
    * existing Connector object, and (2) it is often used to compute a
    * Connector's initializing data.
    *
    * @param[in] base_refinement_ratio
    * @param[in] head_refinement_ratio
    * @param[in] head_width The connector width in the head index space.
    *
    * @return A copy of the connector width converted to the base index
    * space.
    *
    * @pre (base_refinement_ratio >= head_refinement_ratio) ||
    *      (base_refinement_ratio <= head_refinement_ratio)
    */
   static IntVector
   convertHeadWidthToBase(
      const IntVector& base_refinement_ratio,
      const IntVector& head_refinement_ratio,
      const IntVector& head_width);

   // TODO: refactor use of size_t as return type.  This could be
   // problematic.
   /*!
    * @brief Check for consistency between the relationship data and base
    * boxes, and return the number of consistency errors.
    *
    * Consistency stipulates that each neighbor list must correspond to
    * a base box.
    *
    * relationship consistency errors should be treated as fatal because many
    * operations assume consistency.
    */
   size_t
   checkConsistencyWithBase() const;

   /*!
    * @brief Run checkConsistencyWithBase().
    *
    * If any inconsistency is
    * found, write out diagnostic information and throw an
    * unrecoverable assertion is found.
    */
   void
   assertConsistencyWithBase() const;

   /*!
    * @brief Check that the neighbors specified by the relationships exist in
    * the head BoxLevel.
    *
    * If the head is not GLOBALIZED, a temporary copy is made and
    * globalized for checking, triggering communication.
    *
    * @return number of inconsistencies found.
    *
    * @pre getHead().getGlobalizedVersion().getParallelState() == BoxLevel::GLOBALIZED
    */

   size_t
   checkConsistencyWithHead() const;

   /*!
    * @brief Run checkConsistencyWithBase().  If any inconsistency is
    * found, write out diagnostic information and throw an
    * unrecoverable assertion.
    */
   void
   assertConsistencyWithHead() const;

   /*!
    * @brief Compute the differences between two relationship sets.
    *
    * Given Connectors @c left_connector and @c right_connector,
    * compute the relationships that are in @c left_connector but not in
    * @c right_connector.
    *
    * @param[out] left_minus_right
    * @param[in] left_connector
    * @param[in] right_connector
    */
   static void
   computeNeighborhoodDifferences(
      std::shared_ptr<Connector>& left_minus_right,
      const Connector& left_connector,
      const Connector& right_connector);

   /*!
    * @brief Returns true if the transpose of this Connector exists.
    */
   bool
   hasTranspose() const
   {
      return d_transpose;
   }

   /*!
    * @brief Returns the transpose of this Connector if it exists.
    *
    * @pre hasTranspose()
    */
   Connector&
   getTranspose() const
   {
      TBOX_ASSERT(hasTranspose());
      return *d_transpose;
   }

   /*!
    * @brief Sets this Connector's transpose and, if the transpose exists,
    * sets its transpose to this Connector.  If owns_transpose is true then
    * this object will delete the transpose during its deletion.
    *
    * @note If owns_transpose is false, then client code is responsible for
    * the deletion of the transpose.  Similarly, if owns_transpose is true
    * client code must never explicitly delete the transpose.
    *
    * @param[in] transpose
    * @param[in] owns_transpose
    */
   void
   setTranspose(
      Connector* transpose,
      bool owns_transpose)
   {
      if (d_transpose && d_owns_transpose &&
          (d_transpose != transpose) && (d_transpose != this)) {
         delete d_transpose;
      }
      d_transpose = transpose;
      d_owns_transpose = owns_transpose;
      if (d_transpose && d_transpose != this) {
         d_transpose->d_transpose = this;
         d_transpose->d_owns_transpose = false;

         if (d_ratio != IntVector::getOne(d_ratio.getDim())) {
            if ((d_ratio * d_base_width) != d_transpose->d_base_width &&
                d_base_width != (d_ratio * d_transpose->d_base_width)) {

               TBOX_ERROR("Connector::setTranspose: Base width for \n"
                  "this Connector and its transpose are inconsistent.\n");

            }
         } else {
            if (d_base_width != d_transpose->d_base_width) {

               TBOX_ERROR("Connector::setTranspose: Base width for \n"
                  "this Connector and its transpose are inconsistent.\n");

            }
         }
      }
   }

   /*!
    * @brief Check that the relationships are a correct transpose of another
    * Connector and return the number of erroneous relationships.
    *
    * For every relationship in this Connector, there should be a corresponding
    * relationship in the transpose Connector.  Any missing or extra
    * relationship constitutes an error.
    *
    * Errors found are written to perr.
    *
    * @param[in] transpose
    * @param[in] ignore_periodic_relationships
    *
    * @return Global number of errors in assuming that @c transpose is a
    * transpose of @c *this.
    */
   size_t
   checkTransposeCorrectness(
      const Connector& transpose,
      bool ignore_periodic_relationships = false) const;

   /*!
    * @brief Run checkTransposeCorrectness.  If any errors are found,
    * print out diagnostic information and throw an unrecoverable
    * assertion.
    *
    * @param[in] transpose
    * @param[in] ignore_periodic_relationships
    */
   void
   assertTransposeCorrectness(
      const Connector& transpose,
      const bool ignore_periodic_relationships = false) const;

   /*!
    * Check that overlap data is correct (represents overlaps).
    *
    * Checking is done as follows:
    *   - Find overlap errors using @c findOverlapErrors().
    *   - Report overlap errors to @c tbox::perr.
    *   - Return number of local errors.
    *
    * @note
    * This is an expensive operation (it uses @b findOverlapErrors())
    * and should only be used for debugging.
    *
    * @see findOverlaps()
    *
    * @param[in] ignore_self_overlap Ignore a box's overlap with itself
    * @param[in] assert_completeness If false, ignore missing overlaps. This
    *   will still look for overlaps that should not be there.
    * @param[in] ignore_periodic_images If true, do not require neighbors
    *   that are periodic images.
    *
    * @return Number of overlap errors found locally.
    *
    * @pre (getBase().isInitialized()) && (getHead().isInitialized())
    * @pre !hasPeriodicLocalNeighborhoodBaseBoxes()
    */
   int
   checkOverlapCorrectness(
      bool ignore_self_overlap = false,
      bool assert_completeness = true,
      bool ignore_periodic_images = false) const;

   /*!
    * @brief Assert overlap correctness.
    *
    * @par Assertions
    * if an error is found, the method will write out diagnostic information
    * and throw an an error on all processes.
    *
    * This is an expensive check.
    *
    * @see checkOverlapCorrectness().
    *
    * @param[in] ignore_self_overlap
    * @param[in] assert_completeness
    * @param[in] ignore_periodic_images
    *
    * @pre (getBase().isInitialized()) && (getHead().isInitialized())
    */
   void
   assertOverlapCorrectness(
      bool ignore_self_overlap = false,
      bool assert_completeness = true,
      bool ignore_periodic_images = false) const;

   /*!
    * @brief Find errors in overlap data of an overlap Connector.
    *
    * An error is either a missing overlap or an extra overlap.
    *
    * This is an expensive operation and should only be used for
    * debugging.
    *
    * @par Assertions
    * This version throws an assertion only if it finds inconsistent
    * Connector data.  Missing and extra overlaps are returned but do
    * not cause an assertion.
    *
    * @param[out] missing
    * @param[out] extra
    * @param[in] ignore_self_overlap
    *
    * @pre (getBase().isInitialized()) && (getHead().isInitialized())
    */
   void
   findOverlapErrors(
      std::shared_ptr<Connector>& missing,
      std::shared_ptr<Connector>& extra,
      bool ignore_self_overlap = false) const;

   //@}

   /*!
    * @brief Return local number of neighbor sets.
    */
   int
   getLocalNumberOfNeighborSets() const
   {
      return d_relationships.numBoxNeighborhoods();
   }

   /*!
    * @brief Return local number of relationships.
    */
   int
   getLocalNumberOfRelationships() const
   {
      return d_relationships.sumNumNeighbors();
   }

   /*!
    * @brief Return global number of neighbor sets.
    *
    * This requires a global sum reduction, if the global size has not
    * been computed and cached.  When communication is required, all
    * processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * @pre isFinalized()
    */
   int
   getGlobalNumberOfNeighborSets() const
   {
      TBOX_ASSERT(isFinalized());
      cacheGlobalReducedData();
      return d_global_number_of_neighbor_sets;
   }

   /*!
    * @brief Return global number of relationships.
    *
    * This requires a global sum reduction if the global size has not
    * been computed and cached.  When communication is required, all
    * processors must call this method.  To ensure that no
    * communication is needed, call cacheGlobalReducedData() first.
    *
    * @pre isFinalized()
    */
   int
   getGlobalNumberOfRelationships() const
   {
      TBOX_ASSERT(isFinalized());
      cacheGlobalReducedData();
      return d_global_number_of_relationships;
   }

   /*!
    * @brief If global reduced data (global number of relationships,
    * etc.) has not been updated, compute and cache them
    * (communication required).
    *
    * After this method is called, data requiring global reduction can
    * be accessed without further communications, until the object
    * changes.
    *
    * Sets d_global_data_up_to_date;
    *
    * @pre isFinalized()
    */
   void
   cacheGlobalReducedData() const;

   /*!
    * @brief Write the neighborhoods to a restart database.
    *
    * @param[in] restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const
   {
      d_relationships.putToRestart(restart_db);
   }

   /*!
    *
    * @brief Computes refinement ratio between head and base, whether that
    * ratio is exact and whether the head is coarser than the base.
    *
    * @param[in] baseRefinementRatio
    * @param[in] headRefinementRatio
    * @param[out] ratio
    * @param[out] head_coarser
    * @param[out] ratio_is_exact
    */
   static void
   computeRatioInfo(
      const IntVector& baseRefinementRatio,
      const IntVector& headRefinementRatio,
      IntVector& ratio,
      bool& head_coarser,
      bool& ratio_is_exact);

   /*!
    * @brief Writes the neighborhoods to tbox::perr.
    *
    * @param[in] border Left border of the output.
    */
   void
   writeNeighborhoodsToErrorStream(
      const std::string& border) const;

   /*!
    * @brief Writes the requested neighborhood to an output stream.
    *
    * @param[in] os
    * @param[in] box_id
    */
   void
   writeNeighborhoodToStream(
      std::ostream& os,
      const BoxId& box_id) const;

   /*!
    * @brief A class for outputting Connector.
    *
    * This class simplifies the insertion of a Connector into a stream
    * while letting the user control how the Connector is formatted
    * for output.
    *
    * Each Outputter is a light-weight object constructed with a
    * Connector and output parameters.  The Outputter is capable of
    * outputting its Connector, formatted according to the parameters.
    *
    * To use, @see Connector::format(), Connector::formatStatistics().
    */
   class Outputter
   {
      friend std::ostream&
      operator << (
         std::ostream& s,
         const Outputter& f);
private:
      friend class Connector;
      /*!
       * @brief Copy constructor
       */
      Outputter(
         const Outputter& other);
      /*!
       * @brief Construct the Outputter with a Connector and the
       * parameters needed to output the Connector to a stream.
       */
      Outputter(
         const Connector& connector,
         const std::string& border,
         int detail_depth = 2,
         bool output_statistics = false);
      Outputter&
      operator = (
         const Outputter& r);               // Unimplemented private.
      const Connector& d_conn;
      const std::string d_border;
      const int d_detail_depth;
      const bool d_output_statistics;
   };

   /*!
    * @brief Return an object that can format the Connector for
    * insertion into output streams.
    *
    * Usage example:
    * @code
    *    tbox::plog << "my connector:\n" << connector.format() << endl;
    * @endcode
    *
    * @param[in] border
    * @param[in] detail_depth
    */
   Outputter
   format(
      const std::string& border = std::string(),
      int detail_depth = 2) const
   {
      return Outputter(*this, border, detail_depth);
   }

   /*!
    * @brief Return an object that can format the Connector for
    * inserting its global statistics into output streams.
    *
    * Usage example:
    * @code
    *    cout << "my connector statistics:\n"
    *         << connector.formatStatistics("  ") << endl;
    * @endcode
    *
    * @param[in] border
    */
   Outputter
   formatStatistics(
      const std::string& border = std::string()) const
   {
      return Outputter(*this, border, 0, true);
   }

protected:
   /*!
    * @brief Method to do the work of createLocalTranspose.
    *
    * @param transpose
    *
    * @pre transpose
    */
   void
   doLocalTransposeWork(
      Connector* transpose) const;

   /*!
    * @brief Method to do the work of createTranspose.
    *
    * @param transpose
    *
    * @pre transpose
    */
   void
   doTransposeWork(
      Connector* transpose) const;

private:
   // To access findOverlaps_rbbt().
   friend class OverlapConnectorAlgorithm;

   /*
    * Static integer constant descibing class's version number.
    */
   static const int HIER_CONNECTOR_VERSION;

   //! @brief Data structure for MPI reductions.
   struct IntIntStruct { int i;
                         int rank;
   };

   /*
    * Uninitialized default constructor.
    */
   Connector();

   /*!
    * @brief Return the globalized relationship data.
    *
    * @pre getParallelState() == BoxLevel::GLOBALIZED
    */
   const BoxNeighborhoodCollection&
   getGlobalNeighborhoodSets() const
   {
      if (d_parallel_state == BoxLevel::DISTRIBUTED) {
         TBOX_ERROR("Global connectivity unavailable in DISTRIBUTED state." << std::endl);
      }
      return d_global_relationships;
   }

   /*!
    * @brief Return the relationships appropriate to the parallel state.
    *
    * @pre (getParallelState() == BoxLevel::GLOBALIZED) ||
    *      (box_id.getOwnerRank() == getMPI().getRank())
    */
   const BoxNeighborhoodCollection&
   getRelations(
      const BoxId& box_id) const
   {
#ifndef DEBUG_CHECK_ASSERTIONS
      NULL_USE(box_id);
#endif
      if (d_parallel_state == BoxLevel::DISTRIBUTED) {
         TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
      }
      const BoxNeighborhoodCollection& relationships =
         d_parallel_state == BoxLevel::DISTRIBUTED ?
         d_relationships : d_global_relationships;
      return relationships;
   }

   /*!
    * @brief Create a copy of a DISTRIBUTED Connector and
    * change its state to GLOBALIZED.
    *
    * The returned object should be deleted to prevent memory leaks.
    *
    * @pre other.getParallelState() != BoxLevel::GLOBALIZED
    */
   Connector *
   makeGlobalizedCopy(
      const Connector& other) const;

   /*!
    * @brief Get and store info on remote Boxes.
    *
    * This requires global communication (all gather).
    * Call acquireRemoteNeighborhoods_pack to pack up messages.
    * Do an all-gather.  Call acquireRemoteNeighborhoods_unpack
    * to unpack data from other processors.
    */
   void
   acquireRemoteNeighborhoods();

   //! @brief Pack local Boxes into an integer array.
   void
   acquireRemoteNeighborhoods_pack(
      std::vector<int>& send_mesg) const;

   //! @brief Unpack Boxes from an integer array into internal storage.
   void
   acquireRemoteNeighborhoods_unpack(
      const std::vector<int>& recv_mesg,
      const std::vector<int>& proc_offset);

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      t_acquire_remote_relationships = tbox::TimerManager::getManager()->
         getTimer("hier::Connector::acquireRemoteNeighborhoods()");
      t_cache_global_reduced_data = tbox::TimerManager::getManager()->
         getTimer("hier::Connector::cacheGlobalReducedData()");
      t_find_overlaps_rbbt = tbox::TimerManager::getManager()->
         getTimer("hier::Connector::findOverlaps_rbbt()");
   }

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
      t_acquire_remote_relationships.reset();
      t_cache_global_reduced_data.reset();
      t_find_overlaps_rbbt.reset();
   }

   /*!
    * @brief Read the neighborhoods from a restart database.
    *
    * @param[in] restart_db
    */
   void
   getFromRestart(
      tbox::Database& restart_db)
   {
      d_relationships.getFromRestart(restart_db);
   }

   /*!
    * @brief Discover and add overlaps from base and externally
    * provided head box_level.
    *
    * Relationships found are added to appropriate neighbor lists.  No overlap
    * is removed.  If existing overlaps are invalid, remove them first.
    *
    * The provided head should be a globalized version of
    * d_head_handle->getBoxLevel().  It should have the same
    * refinement ratio as d_head_handle->getBoxLevel().
    *
    * The ignore_self_overlap directs the method to not list
    * a box as its own neighbor.  This should be true only
    * when the head and tail objects represent the same box_level
    * (regardless of whether they are the same objects), and
    * you want to disregard self-overlaps.
    * Two boxes are considered the same if
    * - The boxes are equal by comparison (they have the same
    *   owner and the same indices), and
    * - They are from box_levels with the same refinement ratio.
    *
    * @pre head.getParallelState() == BoxLevel::GLOBALIZED
    */
   void
   findOverlaps_rbbt(
      const BoxLevel& head,
      bool ignore_self_overlap = false,
      bool sanity_check_method_postconditions = false);

   /*!
    * @brief Handle for access to the base BoxLevel.
    *
    * We don't use a pointer to the BoxLevel, because it would
    * become dangling when the BoxLevel goes out of scope.
    */
   std::shared_ptr<BoxLevelHandle> d_base_handle;

   /*!
    * @brief Handle for access to the base BoxLevel.
    *
    * We don't use a pointer to the BoxLevel, because it would
    * become dangling when the BoxLevel goes out of scope.
    */
   std::shared_ptr<BoxLevelHandle> d_head_handle;

   /*!
    * @brief Connector width for the base BoxLevel.
    *
    * This is the amount of growth applied to a box in the base BoxLevel
    * before checking if the box overlaps a box in the head BoxLevel.
    */
   IntVector d_base_width;

   /*!
    * @brief Refinement ratio between base and head.
    *
    * If d_head_coarser is false, the head is not coarser than
    * the base and this is the refinement ratio from base to head.
    * If d_head_coarser is true, this is the coarsen ratio
    * from base to head.
    *
    * This is redundant information.  You can compute it
    * from the base and head BoxLevels.
    */
   IntVector d_ratio;

   /*!
    * @brief Whether the ratio between the base and head
    * BoxLevel refinement ratios are exactly as given by
    * d_ratio.  It can only be exact if it can be represented as an
    * IntVector.
    */
   bool d_ratio_is_exact;

   /*!
    * @brief Whether the base BoxLevel is at a finer index space.
    *
    * When this is true, d_ratio is the refinement ratio going
    * from the head to the base.
    *
    * This is redundant information.  You can compute it
    * from the base and head BoxLevels.
    */
   bool d_head_coarser;

   /*!
    * @brief Neighbor data for local Boxes.
    */
   BoxNeighborhoodCollection d_relationships;

   /*!
    * @brief Neighbor data for global Boxes in GLOBALIZED mode.
    */
   BoxNeighborhoodCollection d_global_relationships;

   /*!
    * @brief SAMRAI_MPI object.
    *
    * This is a copy of the getBase().getMPI().  We maintain a copy to
    * allow continued limited functionality should the base detaches
    * itself.
    */
   tbox::SAMRAI_MPI d_mpi;

   /*!
    * @brief State flag.
    *
    * Modified by setParallelState().
    */
   BoxLevel::ParallelState d_parallel_state;

   /*!
    * @brief true when Container's context has been finalized--base, head and
    * width are all defined and consistent/valid.
    */
   bool d_finalized;

   /*!
    * @brief Number of NeighborSets in d_relationships globally.
    */
   mutable int d_global_number_of_neighbor_sets;

   /*!
    * @brief Number of relationships in d_relationships globally.
    */
   mutable int d_global_number_of_relationships;

   /*!
    * @brief Whether globally reduced data is up to date or needs
    * recomputing using cacheGlobalReducedData().
    */
   mutable bool d_global_data_up_to_date;

   Connector* d_transpose;

   bool d_owns_transpose;

   static std::shared_ptr<tbox::Timer> t_acquire_remote_relationships;
   static std::shared_ptr<tbox::Timer> t_cache_global_reduced_data;
   static std::shared_ptr<tbox::Timer> t_find_overlaps_rbbt;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif // included_hier_Connector
