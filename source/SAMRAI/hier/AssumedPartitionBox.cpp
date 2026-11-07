/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fast assumed partition for a box.
 *
 ************************************************************************/
#include "SAMRAI/hier/AssumedPartitionBox.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

/*
 ********************************************************************************
 ********************************************************************************
 */
AssumedPartitionBox::AssumedPartitionBox(
   const Box& box,
   int rank_begin,
   int rank_end,
   int index_begin,
   double avg_parts_per_rank,
   bool interleave):
   d_box(box),
   d_rank_begin(rank_begin),
   d_rank_end(rank_end),
   d_index_begin(index_begin),
   d_index_end(-1),
   d_uniform_partition_size(box.getDim()),
   d_partition_grid_size(box.getDim()),
   d_major(box.getDim()),
   d_index_stride(box.getDim()),
   d_interleave(interleave)
{
   computeLayout(avg_parts_per_rank);
   assignToRanks();
}

/*
 ********************************************************************************
 ********************************************************************************
 */
AssumedPartitionBox::AssumedPartitionBox(
   const tbox::Dimension& dim):
   d_box(dim),
   d_rank_begin(-1),
   d_rank_end(-1),
   d_index_begin(-1),
   d_index_end(-1),
   d_uniform_partition_size(dim),
   d_partition_grid_size(dim),
   d_major(dim),
   d_index_stride(dim),
   d_interleave(false)
{
}

/*
 ********************************************************************************
 * Partition the given box, discarding the current state.
 ********************************************************************************
 */
void
AssumedPartitionBox::partition(
   const Box& box,
   int rank_begin,
   int rank_end,
   int index_begin,
   double avg_parts_per_rank,
   bool interleave)
{
   d_box = box;
   d_rank_begin = rank_begin;
   d_rank_end = rank_end;
   d_index_begin = index_begin;
   d_interleave = interleave;
   computeLayout(avg_parts_per_rank);
   assignToRanks();
}

/*
 ********************************************************************************
 * Return index of first box assigned to given rank.
 ********************************************************************************
 */
int
AssumedPartitionBox::beginOfRank(int rank) const
{
   if (d_interleave) {
      TBOX_ERROR("AssumedPartitionBox::beginOfRank: You should not use beginOfRank()\n"
         << " or endOfRank() for an interleaved AsssumedPartitionBox object\n"
         << " because indices owned by the same rank are not contiguous.\n");
   }
   int index =
      rank < d_rank_begin ? d_index_begin :
      rank < d_first_rank_with_1 ?
      d_first_index_with_2 + (rank - d_rank_begin) * (1 + d_parts_per_rank) :
      rank < d_first_rank_with_0 ?
      d_first_index_with_1 + (rank - d_first_rank_with_1) * d_parts_per_rank :
      d_index_end;
   return index;
}

/*
 ********************************************************************************
 * Return one past index of last box assigned to given rank.
 ********************************************************************************
 */
int
AssumedPartitionBox::endOfRank(int rank) const
{
   if (d_interleave) {
      TBOX_ERROR("AssumedPartitionBox::beginOfRank: You should not use beginOfRank()\n"
         << " or endOfRank() for an interleaved AsssumedPartitionBox object\n"
         << " because indices owned by the same rank are not contiguous.\n");
   }
   int index =
      rank < d_rank_begin ? d_index_begin :
      rank < d_first_rank_with_1 ?
      d_first_index_with_2 + (1 + rank - d_rank_begin) * (1 + d_parts_per_rank) :
      rank < d_first_rank_with_0 ?
      d_first_index_with_1 + (1 + rank - d_first_rank_with_1) * d_parts_per_rank :
      d_index_end;
   return index;
}

/*
 ********************************************************************************
 * Compute the owner of the given index.
 ********************************************************************************
 */
int
AssumedPartitionBox::getOwner(int box_index) const
{
   TBOX_ASSERT(box_index >= d_index_begin);
   TBOX_ASSERT(box_index < d_index_begin + getNumberOfParts());

   int owner = tbox::SAMRAI_MPI::getInvalidRank();

   if (d_interleave) {
      if (box_index >= d_index_begin && box_index < d_index_end) {
         owner = d_rank_begin + (box_index - d_rank_begin) / (d_rank_end - d_rank_begin);
      }
   } else {
      if (box_index < d_index_begin || box_index >= d_index_end) {
         // Not an index in this object ==> invalid owner.
      } else if (box_index < d_first_index_with_1) {
         owner = d_rank_begin + (box_index - d_first_index_with_2) / (1 + d_parts_per_rank);
      } else if (box_index < d_index_end) {
         owner = d_first_rank_with_1 + (box_index - d_first_index_with_1) / d_parts_per_rank;
      }
   }
   return owner;
}

/*
 ********************************************************************************
 * Compute the partition box with the given index.
 ********************************************************************************
 */
Box
AssumedPartitionBox::getBox(int box_index) const
{
   TBOX_ASSERT(box_index >= d_index_begin);
   TBOX_ASSERT(box_index < d_index_begin + d_partition_grid_size.getProduct());

   // Set lower corner in partition grid resolution, based on the box_index.
   Box part(d_box);
   int box_index_diff = box_index - d_index_begin;
   for (int d = d_box.getDim().getValue() - 1; d >= 0; --d) {
      int dir = d_major[d];
      part.setLower(static_cast<tbox::Dimension::dir_t>(dir), box_index_diff / d_index_stride[dir]);
      box_index_diff -= part.lower()[dir] * d_index_stride[dir];
   }

   // Refine lower corner and set upper corner.
   part.setLower(part.lower() * d_uniform_partition_size);
   part.setLower(part.lower() + d_box.lower());
   part.setUpper(part.lower() + d_uniform_partition_size - IntVector::getOne(d_box.getDim()));
   part *= d_box;

   return Box(part, LocalId(box_index), getOwner(box_index));
}

/*
 ********************************************************************************
 * Compute the partition box with the given position in the grid of partitions.
 ********************************************************************************
 */
Box
AssumedPartitionBox::getBox(const IntVector& position) const
{
   TBOX_ASSERT(position >= IntVector::getZero(d_box.getDim()));
   TBOX_ASSERT(position < d_partition_grid_size);

   int box_index = d_index_begin;
   for (int d = 0; d < d_box.getDim().getValue(); ++d) {
      box_index += position[d] * d_index_stride[d];
   }
   const int owner = getOwner(box_index);
   const Index tmp_index(position);
   Box part(tmp_index,
            tmp_index,
            d_box.getBlockId(),
            LocalId(box_index),
            owner);
   part.refine(d_uniform_partition_size);
   part.shift(IntVector(d_box.lower()));
   part *= d_box;
   return part;
}

/*
 ********************************************************************************
 * Compute all partition boxes.
 ********************************************************************************
 */
void
AssumedPartitionBox::getAllBoxes(BoxContainer& all_boxes) const
{
   const int id_begin = begin();
   const int id_end = end();
   for (int id = id_begin; id < id_end; ++id) {
      all_boxes.push_back(getBox(id));
   }
}

/*
 ********************************************************************************
 * Compute all partition boxes owned by the given rank.
 ********************************************************************************
 */
void
AssumedPartitionBox::getAllBoxes(BoxContainer& all_boxes, int rank) const
{
   const int id_begin = beginOfRank(rank);
   const int id_end = endOfRank(rank);
   for (int id = id_begin; id < id_end; ++id) {
      all_boxes.push_back(getBox(id));
   }
}

/*
 ********************************************************************************
 * Find all partition boxes overlapping the given box.  Return whether
 * any boxes were found.
 ********************************************************************************
 */
bool
AssumedPartitionBox::findOverlaps(
   BoxContainer& overlapping_boxes,
   const Box& box) const
{
   Box coarsened_box = box;
   coarsened_box.coarsen(d_uniform_partition_size);
   coarsened_box *= Box(Index(IntVector::getZero(d_box.getDim())),
         Index(d_partition_grid_size - IntVector::getOne(d_box.getDim())),
         d_box.getBlockId());
   for (Box::iterator ci = coarsened_box.begin(); ci != coarsened_box.end(); ++ci) {
      overlapping_boxes.insert(getBox(IntVector(*ci)));
   }
   return !coarsened_box.empty();
}

/*
 ********************************************************************************
 * Check the assumed partition for errors and inconsistencies.  Write
 * error diagnostics to plog.
 *
 * Return number of errors found.  This class should prevent (or at
 * least catch) user errors, so any error found here indicates a bug in
 * the class.
 ********************************************************************************
 */
size_t
AssumedPartitionBox::selfCheck() const
{
   size_t nerr = 0;

   BoxContainer all_parts;
   for (int box_index = begin(); box_index != end(); ++box_index) {
      const Box box = getBox(box_index);
      if (box.getOwnerRank() == tbox::SAMRAI_MPI::getInvalidRank() ||
          box.getOwnerRank() < d_rank_begin || box.getOwnerRank() >= d_rank_end) {
         ++nerr;
         tbox::plog << "AssumedPartitionerBox::selfCheck(): Box "
                    << box << " has owner outside expected range of ["
                    << d_rank_begin << ',' << d_rank_end << ')'
                    << std::endl;
      }
      all_parts.pushBack(getBox(box_index));
   }
   all_parts.makeTree();

   // All parts should match boxes gotten through partition grid indices.
   const Box partition_grid(Index(IntVector::getZero(d_box.getDim())),
                            Index(d_partition_grid_size - IntVector::getOne(d_box.getDim())),
                            d_box.getBlockId());
   BoxContainer all_parts_by_grid;
   for (Box::iterator gi = partition_grid.begin(); gi != partition_grid.end(); ++gi) {
      all_parts_by_grid.pushBack(getBox(IntVector(*gi)));
   }
   // all_parts_by_grid may be in different order, so sort before comparing.
   all_parts.order();
   all_parts_by_grid.order();
   if (all_parts_by_grid != all_parts) {
      ++nerr;
      tbox::plog << "AssumedPartitionerBox::selfCheck(): Boxes gotten by\n"
                 << "index loop differs from boxes gotten by grid loop.\n"
                 << "Boxes by index loop:\n" << all_parts.format("\t") << '\n'
                 << "Boxes by grid loop:\n" << all_parts_by_grid.format("\t") << '\n';
   }
   all_parts.unorder();
   all_parts_by_grid.clear();

   // Parts should not overlap each other.
   BoxContainer tmp_boxes;
   for (BoxContainer::const_iterator bi = all_parts.begin(); bi != all_parts.end(); ++bi) {
      const Box& box = *bi;
      all_parts.findOverlapBoxes(tmp_boxes, box);
      tmp_boxes.order();
      if (!tmp_boxes.empty()) {
         BoxContainer::iterator self = tmp_boxes.find(box);
         if (self != tmp_boxes.end()) {
            tmp_boxes.erase(self);
         }
      }
      if (!tmp_boxes.empty()) {
         nerr += tmp_boxes.size();
         tbox::plog << "AssumedPartitionerBox::selfCheck(): Box "
                    << box << " unexpectedly overlaps these:\n"
                    << tmp_boxes.format("\t") << std::endl;
      }
      tmp_boxes.clear();
   }

   // Parts should cover no less than d_box.
   BoxContainer box_leftover(d_box);
   box_leftover.removeIntersections(all_parts);
   if (!box_leftover.empty()) {
      nerr += box_leftover.size();
      tbox::plog << "AssumedPartitionerBox::selfCheck(): Partitions cover less than box "
                 << d_box << "  Portions not covered by partitions:\n"
                 << box_leftover.format("\t") << std::endl;
   }

   // Parts should cover no more than d_box.
   BoxContainer parts_leftover = all_parts;
   parts_leftover.removeIntersections(d_box);
   if (!parts_leftover.empty()) {
      nerr += parts_leftover.size();
      tbox::plog << "AssumedPartitionerBox::selfCheck(): Partitions cover more than box "
                 << d_box << "  Portions outside the box:\n"
                 << parts_leftover.format("\t") << std::endl;
   }

   if (!d_interleave) {
      for (int rank = d_rank_begin; rank < d_rank_end; ++rank) {
         const int ibegin = beginOfRank(rank);
         const int iend = endOfRank(rank);
         if (rank > d_rank_begin) {
            if (ibegin < endOfRank(rank - 1)) {
               ++nerr;
               tbox::plog << "AssumedPartitionBox::selfCheck(): Index ranges overlap.\n"
                          << "Rank " << rank << " has [" << ibegin << ',' << iend << ")\n"
                          << "Rank " << rank - 1 << " has [" << beginOfRank(rank - 1)
                          << ',' << endOfRank(rank - 1) << ")\n";
            }
         }
         const int number_owned = iend - ibegin;
         const int expected_number =
            rank < d_first_rank_with_1 ? d_parts_per_rank + 1 :
            rank < d_first_rank_with_0 ? d_parts_per_rank : 0;
         if (number_owned != expected_number) {
            ++nerr;
            tbox::plog << "AssumedPartition::selfCheck(): Wrong parts count for rank " << rank
                       << ".\n"
                       << "It should own " << expected_number << ", but it owns " << number_owned
                       << ".\n";
         }
         for (int id = ibegin; id < iend; ++id) {
            const int owner = getOwner(id);
            if (owner != rank) {
               ++nerr;
               tbox::plog << "AssumedPartition::selfCheck(): Wrong owner for id " << id << ".\n"
                          << "Id " << id << " is owned by " << owner << " but rank " << rank
                          << " is assined the range (" << ibegin << ',' << iend << ").\n";
            }
         }
      }
   }

   return nerr;
}

/*
 ********************************************************************************
 * Compute the partition lay-out.  We use a grid of uniform sized
 * partitions whose union covers d_box and as little else as possible.
 */
void
AssumedPartitionBox::computeLayout(double avg_parts_per_rank)
{
   const IntVector box_size = d_box.numberCells();

   const int num_ranks = d_rank_end - d_rank_begin;

   if (box_size == 0) {
      d_uniform_partition_size = d_partition_grid_size = IntVector::getZero(d_box.getDim());
   } else {
      /*
       * Compute uniform partition size and how many partitions in each
       * direction.  There isn't one correct lay-out, but we try to avoid
       * excessive aspect ratios.
       */
      const int target_parts_count = static_cast<int>(num_ranks * avg_parts_per_rank + 0.5);
      d_uniform_partition_size = box_size;
      d_partition_grid_size = IntVector::getOne(d_box.getDim());
      long int parts_count = d_partition_grid_size.getProduct();
      IntVector num_parts_can_increase(d_box.getDim(), 1);
      IntVector sorter(d_box.getDim());
      while (parts_count < target_parts_count &&
             num_parts_can_increase != IntVector::getZero(d_box.getDim())) {
         sorter.sortIntVector(d_uniform_partition_size);
         int inc_dir = 0;
         for (inc_dir = d_box.getDim().getValue() - 1; inc_dir >= 0; --inc_dir) {
            if (num_parts_can_increase[sorter[inc_dir]]) break;
         }
         inc_dir = sorter[inc_dir];

         // Double partition grid size, unless it causes too many partitions.
         if (2 * parts_count > target_parts_count) {
            const long int cross_section = parts_count / d_partition_grid_size[inc_dir];
            d_partition_grid_size[inc_dir] =
               static_cast<int>((target_parts_count + cross_section - 1) / cross_section);
            parts_count = d_partition_grid_size.getProduct();
         } else {
            d_partition_grid_size[inc_dir] *= 2;
            parts_count *= 2;
         }

         d_uniform_partition_size = IntVector::ceilingDivide(box_size, d_partition_grid_size);
         num_parts_can_increase[inc_dir] = d_uniform_partition_size[inc_dir] > 1;
      }
      TBOX_ASSERT(parts_count == d_partition_grid_size.getProduct());

      // There can be partitions completele outside d_box.  Remove them.
      d_partition_grid_size = IntVector::ceilingDivide(box_size, d_uniform_partition_size);
   }

   d_index_end = d_index_begin + static_cast<int>(d_partition_grid_size.getProduct());

   d_major.sortIntVector(d_partition_grid_size);

   for (int d = 0; d < d_box.getDim().getValue(); ++d) {
      int dir = d_major[d];
      d_index_stride[dir] = 1;
      for (int d1 = d - 1; d1 >= 0; --d1) {
         d_index_stride[dir] *= d_partition_grid_size[d_major[d1]];
      }
   }
}

/*
 ********************************************************************************
 * Compute rank assignment for the partition lay-out.
 ********************************************************************************
 */
void
AssumedPartitionBox::assignToRanks()
{
   if (d_interleave) {
      // No special settings needed.
   } else {
      assignToRanks_contiguous();
   }
}

/*
 ********************************************************************************
 * Compute rank assignment for the partition lay-out where each rank in
 * [d_rank_begin,d_rank_end) gets a contiguous set of box indices in
 * [d_index_begin,d_index_end).
 *
 * Each rank has 0, d_parts_per_rank or 1+d_parts_per_rank partitions.
 * Lower ranks have more partitions than higher ranks do.
 ********************************************************************************
 */
void
AssumedPartitionBox::assignToRanks_contiguous()
{
   TBOX_ASSERT(d_index_end - d_index_begin == d_partition_grid_size.getProduct());

   /*
    * If there are more ranks than parts, the first getNumberOfParts()
    * ranks have 1 part each and the rest have none.  If there are
    * more parts than ranks, lower ranks have (d_parts_per_rank_1+1)
    * parts each and higher ranks have d_parts_per_rank each.  In
    * second case, index vs rank looks like this:
    *
    * index ^
    *       |
    *    i0 |                ......
    *       |             .
    *       |          .
    *       |       .
    *    i1 |    .
    *       |   .
    *       |  .
    *       | .
    *       |.
    *    i2 +-----------------------> rank
    *       r2   r1          r0
    *
    *    (r2,i2) = first rank with 1+d_parts_per_rank parts, its first index
    *    (r1,i1) = first rank with d_parts_per_rank parts, its first index
    *    (r0,i0) = first rank with 0 parts, d_index_end
    *    if r1==r2, no rank has 1+d_parts_per_rank parts (more ranks than parts)
    *    if r0==r1, no rank has d_parts_per_rank parts
    */
   if (d_index_end - d_index_begin <= d_rank_end - d_rank_begin) {
      d_parts_per_rank = 1;
      d_first_rank_with_1 = d_rank_begin;
      d_first_rank_with_0 = d_rank_begin + static_cast<int>(d_partition_grid_size.getProduct());
      d_first_index_with_1 = d_first_index_with_2 = d_index_begin;
   } else {
      d_parts_per_rank = (d_index_end - d_index_begin) / (d_rank_end - d_rank_begin);
      d_first_index_with_2 = d_index_begin;
      d_first_rank_with_1 = d_rank_begin
         + (d_index_end - d_index_begin) % (d_rank_end - d_rank_begin);
      d_first_index_with_1 = d_first_index_with_2
         + (1 + d_parts_per_rank) * (d_first_rank_with_1 - d_rank_begin);
      d_first_rank_with_0 = d_rank_end;
   }
}

/*
 ********************************************************************************
 ********************************************************************************
 */
void
AssumedPartitionBox::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   const char* to = "..";
   co << border << "d_box = " << d_box << "    "
      << d_box.numberCells() << "|" << d_box.size()
      << '\n' << border
      << "d_uniform_partition_size = " << d_uniform_partition_size
      << "  d_partition_grid_size = " << d_partition_grid_size
      << "  d_major = " << d_major
      << "  d_index_stride = " << d_index_stride
      << '\n' << border
      << "indices (" << d_partition_grid_size.getProduct() << "): "
      << d_index_begin << to << d_index_begin + d_partition_grid_size.getProduct() - 1
      << "    ranks (" << d_rank_end - d_rank_begin << "): "
      << d_rank_begin << to << d_rank_end - 1
      << "    parts_per_rank=" << d_parts_per_rank
      << "    interleave=" << d_interleave
      << '\n' << border
      << "ranks with " << d_parts_per_rank + 1 << " parts (" << d_first_rank_with_1
   - d_rank_begin << "): "
      << d_rank_begin << to << d_first_rank_with_1 - 1
      << '\n' << border
      << "ranks with " << d_parts_per_rank << " parts (" << d_first_rank_with_0
   - d_first_rank_with_1 << "): "
      << d_first_rank_with_1 << to << d_first_rank_with_0 - 1
      << '\n' << border
      << "ranks with 0 parts (" << d_rank_end - d_first_rank_with_0 << "): "
      << d_first_rank_with_0 << to << d_rank_end - 1
      << '\n' << border
      << "d_first_rank_with_1=" << d_first_rank_with_1
      << "    d_first_rank_with_0=" << d_first_rank_with_0
      << "    d_first_index_with_1=" << d_first_index_with_1
      << "    d_first_index_with_2=" << d_first_index_with_2
      << std::endl;
   if (detail_depth > 0) {
      BoxContainer parts;
      for (int box_index = begin(); box_index != end(); ++box_index) {
         parts.pushBack(getBox(box_index));
      }
      co << border << "Parts: " << parts.format(border);
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
