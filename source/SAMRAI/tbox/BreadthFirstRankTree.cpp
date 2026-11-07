/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/BreadthFirstRankTree.h"

#include "SAMRAI/tbox/PIO.h"
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
BreadthFirstRankTree::BreadthFirstRankTree():
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_children(2, getInvalidRank()),
   d_num_children(0),
   d_child_number(getInvalidChildNumber()),
   d_generation(0),
   d_root_rank(getInvalidRank())
{
}

/*
 ****************************************************************
 ****************************************************************
 */
BreadthFirstRankTree::BreadthFirstRankTree(
   int first_rank,
   int last_rank,
   int my_rank,
   unsigned int degree):
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_children(degree, getInvalidRank()),
   d_num_children(0),
   d_child_number(getInvalidChildNumber()),
   d_generation(0),
   d_root_rank(getInvalidRank())
{
   setupTreeForContiguousRanks(first_rank, last_rank, my_rank);
}

/*
 ****************************************************************
 ****************************************************************
 */
BreadthFirstRankTree::~BreadthFirstRankTree()
{
}

/*
 ****************************************************************
 * Set up the tree from a RankGroup.
 ****************************************************************
 */
void
BreadthFirstRankTree::setupTree(
   const RankGroup& rank_group,
   int my_rank)
{
   TBOX_ASSERT(rank_group.isMember(my_rank));
   setupTreeForContiguousRanks(0, rank_group.size() - 1, rank_group.getMapIndex(my_rank));
}

/*
 ****************************************************************
 * Set up the tree for contiguous ranks.
 ****************************************************************
 */
void
BreadthFirstRankTree::setupTreeForContiguousRanks(
   int first_rank,
   int last_rank,
   int rank)
{
   TBOX_ASSERT(first_rank <= rank);
   TBOX_ASSERT(rank <= last_rank);

   const unsigned int degree = static_cast<unsigned int>(d_children.size());

   d_root_rank = first_rank;
   d_rank = rank;
   if (d_rank > first_rank) {
      d_parent = (d_rank - first_rank - 1) / degree;
      d_child_number = (d_rank - first_rank - 1) % degree;
   } else {
      d_parent = getInvalidRank();
      d_child_number = getInvalidChildNumber();
   }

   d_num_children = 0;
   for (d_num_children = 0; d_num_children < degree; ++d_num_children) {

      d_children[d_num_children] =
         first_rank + degree * (d_rank - first_rank) + 1 + d_num_children;

      if (d_children[d_num_children] > last_rank) {
         for (unsigned int i = d_num_children; i < degree; ++i) {
            d_children[i] = getInvalidRank();
         }
         break;
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
