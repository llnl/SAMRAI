/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a single box.
 *
 ************************************************************************/
#ifndef included_hier_AssumedPartitionBox
#define included_hier_AssumedPartitionBox

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/IntVector.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Compute an assumed partition of a box.  The assumed
 * partition should very fast to create and query and requires minimal
 * storage.
 *
 * An assumed partition avoids extreme imbalances, but its purpose is
 * not fine load balancing.
 */
class AssumedPartitionBox
{

public:
   /*!
    * @brief Construct AssumedPartition from a box.
    *
    * @param[in] unpartitioned_box Incoming unpartitioned box
    *
    * @param[in] rank_begin First rank
    *
    * @param[in] rank_end One past last rank
    *
    * @param[in] index_begin First index
    *
    * @param[in] avg_parts_per_rank See partition()
    *
    * @param[in] interleave See partition()
    */
   AssumedPartitionBox(
      const Box& unpartitioned_box,
      int rank_begin,
      int rank_end,
      int index_begin = 0,
      double avg_parts_per_rank = 1.0,
      bool interleave = false);

   /*!
    * @brief Constructor.
    */
   AssumedPartitionBox(
      const tbox::Dimension& dim);

   /*!
    * @brief Partition the given box, discarding the current state.
    *
    * The partition should degerate correctly if the box is empty, i.e.,
    * the partition size and count should be zero.
    *
    * @param[in] unpartitioned_box Incoming unpartitioned box
    * @param[in] rank_begin First rank
    * @param[in] rank_end One past last rank
    * @param[in] index_begin
    *
    * @param[in] avg_parts_per_rank Algorithm normally tries to get
    * one partition per rank.  This parameter is a request to change
    * that.
    *
    * @param[in] interleave Algorithm normally assign consecutive box
    * indices to a process.  This flag causes it to interleave
    * (round-robin) the box assignments.
    */
   void
   partition(
      const Box& unpartitioned_box,
      int rank_begin,
      int rank_end,
      int index_begin = 0,
      double avg_parts_per_rank = 1.0,
      bool interleave = false);

   /*!
    * @brief Destructor.
    */
   ~AssumedPartitionBox() = default;

   //! @brief Return the original unpartitioned box.
   const Box& getUnpartitionedBox() const {
      return d_box;
   }

   //! @brief Number of box partitions.
   int getNumberOfParts() const {
      return d_index_end - d_index_begin;
   }

   //! @brief Return the owner for a box.
   int
   getOwner(
      int box_index) const;

   //! @brief Return box for given index.
   Box
   getBox(
      int box_index) const;

   //! @brief Return box for given partition's position in the partition grid.
   Box
   getBox(
      const IntVector& position) const;

   //! @brief Get all partition boxes.
   void
   getAllBoxes(
      BoxContainer& all_boxes) const;

   //! @brief Get all partition boxes owned by a given rank.
   void
   getAllBoxes(
      BoxContainer& all_boxes,
      int rank) const;

   //! @brief Return index of first partition box.
   int begin() const {
      return d_index_begin;
   }

   //! @brief Return one past index of last partition box.
   int end() const {
      return d_index_end;
   }

   /*!
    * @brief Return index of first box in the contiguous index range
    * assigned to given rank.
    *
    * This method should not be used for objects when the ranks are
    * interleaved.
    *
    * If rank is lower than ranks in partitioning, return first index.
    * If rank is higher than ranks in partitioning, return one past
    * last index.
    */
   int
   beginOfRank(
      int rank) const;

   /*!
    * @brief Return one past index of last box in the contiguous index
    * range assigned to given rank.
    *
    * This method should not be used for objects when the ranks are
    * interleaved.
    *
    * If rank is lower than ranks in partitioning, return first index.
    * If rank is higher than ranks in partitioning, return one past
    * last index.
    */
   int
   endOfRank(
      int rank) const;

   /*!
    * @brief Find box partitions overlapping the given box.
    *
    * The search cost is proportional to number of overlapping boxes
    * found, NOT the total number of partitions.
    *
    * @param[out] overlapping_boxes
    * @param[in] box
    *
    * @return Whether any partitions are found.
    */
   bool
   findOverlaps(
      BoxContainer& overlapping_boxes,
      const Box& box) const;

   /*!
    * @brief Check the assumed partition for errors and
    * inconsistencies.  Write error diagnostics to plog.
    *
    * @return Number of errors found.  (Errors indicate
    * a bug in this class.)
    */
   size_t
   selfCheck() const;

   /*!
    * @brief Print info from this object
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

private:
   //! @brief Compute the partition lay-out.
   void
   computeLayout(
      double avg_parts_per_rank);

   //! @brief Compute rank assignment for the partition lay-out.
   void
   assignToRanks();

   /*!
    * @brief Compute rank assignment for the partition lay-out, using
    * contiguous index assignments.
    */
   void
   assignToRanks_contiguous();

   //! @brief Box partitioned.
   Box d_box;
   //! @brief First rank.
   int d_rank_begin;
   //! @brief One past last rank.
   int d_rank_end;
   //! @brief Index for first partition box.
   int d_index_begin;
   //! @brief One past index of last partition box.
   int d_index_end;

   //! @brief Size of each uniform partition.
   IntVector d_uniform_partition_size;
   //! @brief Number of partitions in each direction (size of partition grid).
   IntVector d_partition_grid_size;

   //! @brief Directions sorted from small to big, in d_partition_grid_size.
   IntVector d_major;
   //! @brief Stride of box index in each direction of the partition grid.
   IntVector d_index_stride;

   //! @brief Whether to interleave box assignments using round-robin.
   bool d_interleave;

   //@{
   //! @name Parameters for partition assignment in non-interleaved mode

   /*!
    * @brief Min (or max) parts per rank when there are more (fewer)
    * parts than ranks.
    */
   int d_parts_per_rank;

   //! @see assignToRanks_contiguous()
   int d_first_rank_with_1;

   //! @see assignToRanks_contiguous()
   int d_first_rank_with_0;

   //! @see assignToRanks_contiguous()
   int d_first_index_with_2;

   //! @see assignToRanks_contiguous()
   int d_first_index_with_1;

   //@}

};

}
}

#endif  // included_hier_AssumedPartitionBox
