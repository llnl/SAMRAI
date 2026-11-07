/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A class to manage groups of processor ranks
 *
 ************************************************************************/

#include "SAMRAI/tbox/RankGroup.h"

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 * Default constructor creates RankGroup that uses all available ranks.
 ***********************************************************************
 */

RankGroup::RankGroup():
   d_min(-1),
   d_max(-1),
   d_ranks(0),
   d_storage(USING_ALL),
   d_samrai_mpi(SAMRAI_MPI::getSAMRAIWorld())
{
}

/*
 ***********************************************************************
 * Constructor that creates RankGroup that uses all available ranks in
 * a given communicator
 ***********************************************************************
 */
RankGroup::RankGroup(
   const SAMRAI_MPI& samrai_mpi):
   d_min(-1),
   d_max(-1),
   d_ranks(0),
   d_storage(USING_ALL),
   d_samrai_mpi(samrai_mpi)
{
}

/*
 ***********************************************************************
 * Constructor that takes a min and max rank.
 ***********************************************************************
 */

RankGroup::RankGroup(
   const int min,
   const int max,
   const SAMRAI_MPI& samrai_mpi):
   d_min(min),
   d_max(max),
   d_ranks(0),
   d_storage(USING_MIN_MAX),
   d_samrai_mpi(samrai_mpi)
{
   int nodes = 1;
   samrai_mpi.Comm_size(&nodes);

   TBOX_ASSERT(min >= 0);
   TBOX_ASSERT(max < nodes);
   TBOX_ASSERT(min <= max);

   /*
    * If min and max cover the full set of ranks, then switch to
    * "using_all" mode.
    */
   if (min == 0 && max == nodes - 1) {
      d_min = -1;
      d_max = -1;
      d_storage = USING_ALL;
   }
}

/*
 ***********************************************************************
 * Constructor that takes an array of ranks.
 ***********************************************************************
 */

RankGroup::RankGroup(
   const std::vector<int>& rank_group,
   const SAMRAI_MPI& samrai_mpi):
   d_min(-1),
   d_max(-1),
   d_ranks(rank_group),
   d_storage(USING_ARRAY),
   d_samrai_mpi(samrai_mpi)
{
   TBOX_ASSERT(rank_group.size() > 0);

   int nodes = 1;
   samrai_mpi.Comm_size(&nodes);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(static_cast<int>(rank_group.size()) <= nodes);

   /*
    * Check that each entry in the array has a unique value and is increasing
    * order
    */
   for (int i = 0; i < static_cast<int>(rank_group.size()); ++i) {
      TBOX_ASSERT(rank_group[i] >= 0);
      TBOX_ASSERT(rank_group[i] < nodes);
      if (i > 0) {
         TBOX_ASSERT(rank_group[i] > rank_group[i - 1]);
      }
   }
#endif

   /*
    * If array is the full set of ranks, then switch to "using_all" mode.
    */
   if (static_cast<int>(rank_group.size()) == nodes) {
      d_ranks.resize(0);
      d_storage = USING_ALL;
   }
}

/*
 ***********************************************************************
 * Copy constructor.
 ***********************************************************************
 */
RankGroup::RankGroup(
   const RankGroup& other):
   d_min(other.d_min),
   d_max(other.d_max),
   d_ranks(other.d_ranks),
   d_storage(other.d_storage),
   d_samrai_mpi(other.d_samrai_mpi)
{
}

/*
 ***********************************************************************
 * Destructor
 ***********************************************************************
 */

RankGroup::~RankGroup()
{
}

/*
 ***********************************************************************
 * Assignment operator.
 ***********************************************************************
 */
RankGroup&
RankGroup::operator = (
   const RankGroup& rhs)
{
   d_min = rhs.d_min;
   d_max = rhs.d_max;
   d_ranks = rhs.d_ranks;
   d_storage = rhs.d_storage;
   d_samrai_mpi = rhs.d_samrai_mpi;
   return *this;
}

/*
 ***********************************************************************
 * Determine if given rank is in the group.
 ***********************************************************************
 */

bool
RankGroup::isMember(
   const int rank) const
{
   bool is_member = false;

   switch (d_storage) {

      case USING_ALL:
         is_member = true;
         break;

      case USING_ARRAY:
      {
         int lo = 0;
         int hi = static_cast<int>(d_ranks.size()) - 1;
         while (!is_member && hi >= lo) {
            int i = (lo + hi) / 2;
            if (rank < d_ranks[i]) {
               hi = i - 1;
            } else if (rank > d_ranks[i]) {
               lo = i + 1;
            } else {
               is_member = true;
            }
         }
      }
      break;

      case USING_MIN_MAX:
         if (rank >= d_min && rank <= d_max) {
            is_member = true;
         }
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
   }

   return is_member;
}

/*
 ***********************************************************************
 * Size of RankGroup.
 ***********************************************************************
 */

int
RankGroup::size() const
{
   int size = -MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         d_samrai_mpi.Comm_size(&size);
         break;

      case USING_ARRAY:
         size = static_cast<int>(d_ranks.size());
         break;

      case USING_MIN_MAX:
         size = d_max - d_min + 1;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
         break;
   }

   return size;
}

/*
 ***********************************************************************
 * Get a rank using 1-to-1 mapping.
 ***********************************************************************
 */

int
RankGroup::getMappedRank(
   const int index) const
{
   TBOX_ASSERT(index >= 0 && index < size());

   int mapped_rank = -MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         mapped_rank = index;
         break;

      case USING_ARRAY:
         mapped_rank = d_ranks[index];
         break;

      case USING_MIN_MAX:
         mapped_rank = d_min + index;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
         break;
   }

   return mapped_rank;
}

/*
 ***********************************************************************
 * Get an integer identifier from 1-to-1 mapping.
 ***********************************************************************
 */

int
RankGroup::getMapIndex(
   const int rank) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   int nodes = 1;
   d_samrai_mpi.Comm_size(&nodes);
   TBOX_ASSERT(rank >= 0 && rank < nodes);
#endif

   TBOX_ASSERT(isMember(rank));

   int map_id = -MathUtilities<int>::getMax();

   switch (d_storage) {

      case USING_ALL:
         map_id = rank;
         break;

      case USING_ARRAY:
      {
         int lo = 0;
         int hi = static_cast<int>(d_ranks.size()) - 1;
         while (hi >= lo) {
            int i = (lo + hi) / 2;
            if (rank < d_ranks[i]) {
               hi = i - 1;
            } else if (rank > d_ranks[i]) {
               lo = i + 1;
            } else {
               map_id = i;
               break;
            }
         }
      }
      break;

      case USING_MIN_MAX:
         map_id = rank - d_min;
         break;

      default:
         TBOX_ERROR("RankGroup has not been set with a valid storage method");
         break;
   }

   return map_id;
}

}
}
