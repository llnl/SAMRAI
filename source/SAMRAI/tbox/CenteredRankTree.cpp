/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/CenteredRankTree.h"

#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

/*
 ****************************************************************
 ****************************************************************
 */
CenteredRankTree::CenteredRankTree():
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_root_rank(getInvalidRank()),
   d_num_children(0),
   d_make_first_rank_the_root(false)
{
}

/*
 ****************************************************************
 ****************************************************************
 */
CenteredRankTree::CenteredRankTree(
   const SAMRAI_MPI& mpi,
   bool make_first_rank_the_root):
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_root_rank(getInvalidRank()),
   d_num_children(0),
   d_make_first_rank_the_root(make_first_rank_the_root)
{
   setupTreeForContiguousRanks(0, mpi.getSize() - 1, mpi.getRank());
}

/*
 ****************************************************************
 ****************************************************************
 */
CenteredRankTree::CenteredRankTree(
   int first_rank,
   int last_rank,
   int my_rank,
   bool make_first_rank_the_root):
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_root_rank(getInvalidRank()),
   d_num_children(0),
   d_make_first_rank_the_root(make_first_rank_the_root)
{
   setupTreeForContiguousRanks(first_rank, last_rank, my_rank);
}

/*
 ****************************************************************
 ****************************************************************
 */
CenteredRankTree::~CenteredRankTree()
{
}

/*
 ****************************************************************
 * Set up the tree from a RankGroup.
 ****************************************************************
 */
void
CenteredRankTree::setupTree(
   const RankGroup& rank_group,
   int my_rank)
{
   TBOX_ASSERT(rank_group.isMember(my_rank));
   setupTreeForContiguousRanks(0, rank_group.size() - 1, rank_group.getMapIndex(my_rank));
}

/*
 ****************************************************************
 * Set up the tree for contiguous ranks.
 *
 * Determine the parent and children for the argument rank.  Start
 * with the range [first_rank,last_rank] and narrow it down by
 * moving in the direction of the rank until we find it.
 *
 * How to determine which direction to go:
 * Given contiguous ranks [rbeg, rend], the "midpoint" is
 * rmid=(rbeg+rend+1)/2.  Select the root and branches
 * for the range like this:
 *
 *            rmid
 *           /    \
 *  [rbeg,rmid)  (rmid,rend]
 *
 * The midpoint is the root for this range.  (Note that the midpoint is
 * rounded up, so that the left branch gets the odd rank.)
 *
 * Exception: If d_make_first_rank_the_root and we are at the root of
 * the tree, select the root and branches like this: rmid=(rbeg+rend)2,
 *
 *            rbeg
 *           /    \
 *  (rbeg,rmid]  (rmid,rend]
 ****************************************************************
 */
void
CenteredRankTree::setupTreeForContiguousRanks(
   int first_rank,
   int last_rank,
   int rank)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(first_rank <= rank);
   TBOX_ASSERT(rank <= last_rank);
#endif

   d_rank = rank;
   d_parent = d_children[0] = d_children[1] = getInvalidRank();
   d_generation = 0;
   d_child_number = getInvalidChildNumber();

   int rbeg = first_rank;
   int rend = last_rank;

   /*
    * Set up the root.  If it is d_rank, compute its children.
    */
   if (d_make_first_rank_the_root) {
      d_root_rank = rbeg;
      rbeg = rbeg + 1;
      const int rmid = (rbeg + rend) / 2;
      if (d_rank == d_root_rank) {
         if (rbeg <= rend) d_children[0] = (rbeg + rmid + 1) / 2;
         if (rend > rbeg) d_children[1] = (rmid + 1 + rend + 1) / 2;
      } else if (d_rank <= rmid) {
         rend = rmid;
         d_child_number = 0;
      } else { /* d_rank > rmid */
         rbeg = rmid + 1;
         d_child_number = 1;
      }
   } else {
      const int rmid = (rbeg + rend + 1) / 2;
      d_root_rank = rmid;
      if (d_rank == d_root_rank) {
         if (d_rank > rbeg) d_children[0] = (rbeg + rmid) / 2;
         if (d_rank < rend) d_children[1] = (rmid + 1 + rend + 1) / 2;
      }
   }

   if (d_rank != d_root_rank) {
      // Find d_rank's parent and children by walking toward it.
      d_parent = d_root_rank;
      ++d_generation;
      while (true) {
         const int rmid = (rbeg + rend + 1) / 2;
         if (d_rank == rmid) {
            if (d_rank > rbeg) d_children[0] = (rbeg + rmid) / 2;
            if (d_rank < rend) d_children[1] = (rmid + 1 + rend + 1) / 2;
            break;
         } else if (d_rank < rmid) {
            rend = rmid - 1;
            d_child_number = 0;
         } else { /* d_rank > rmid */
            rbeg = rmid + 1;
            d_child_number = 1;
         }
         d_parent = rmid;
         ++d_generation;
      }
   }

   // Count number of children.
   d_num_children = 0;
   for (int i = 0; i < 2; ++i) {
      if (d_children[i] != getInvalidRank()) {
         ++d_num_children;
      }
   }

}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
