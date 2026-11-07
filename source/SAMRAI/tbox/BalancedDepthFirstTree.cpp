/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"

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
BalancedDepthFirstTree::BalancedDepthFirstTree():
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_root_rank(getInvalidRank()),
   d_num_children(0),
   d_do_left_leaf_switch(true)
{
}

/*
 ****************************************************************
 ****************************************************************
 */
BalancedDepthFirstTree::BalancedDepthFirstTree(
   int first_rank,
   int last_rank,
   int my_rank,
   bool do_left_leaf_switch):
   d_rank(getInvalidRank()),
   d_parent(getInvalidRank()),
   d_root_rank(getInvalidRank()),
   d_num_children(0),
   d_do_left_leaf_switch(do_left_leaf_switch)
{
   setupTreeForContiguousRanks(first_rank, last_rank, my_rank);
}

/*
 ****************************************************************
 ****************************************************************
 */
BalancedDepthFirstTree::~BalancedDepthFirstTree()
{
}

/*
 ****************************************************************
 * Set up the tree from a RankGroup.
 ****************************************************************
 */
void
BalancedDepthFirstTree::setupTree(
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
BalancedDepthFirstTree::setupTreeForContiguousRanks(
   int first_rank,
   int last_rank,
   int rank)
{
   TBOX_ASSERT(first_rank <= rank);
   TBOX_ASSERT(rank <= last_rank);

   d_root_rank = first_rank;
   d_generation = 0;
   d_child_number = getInvalidChildNumber();

   int rbeg = first_rank;
   int rend = last_rank;
   int up = getInvalidRank();  // Temporary guess for parent.
   int upp = getInvalidRank(); // Temporary guess for grandparent.
   int cl, cr;         // Temporary guesses for children.
   unsigned int parents_child_number = getInvalidChildNumber();
   bool is_switchable = false;  // Part of a left-leaf switchable trio.

   size_t nr;           // Number of nodes on right branch
   size_t nl;           // Number of nodes on left branch

   while (1) {

      /*
       * Walk from root to leaf to find the position of rank, its
       * parent and its children.
       */

      int node = rbeg; // Node being examined
      size_t nrem = static_cast<size_t>(rend - rbeg);  // Number or nodes remaining, excluding node.

      nr = nrem / 2;      // Number on right branch
      nl = nrem - nr;     // Number on left branch

      /*
       * Both children are leaves => parent and children make
       * a switchable trio.
       */
      if (nrem == 2) {
         is_switchable = true;
      }

      cl = getInvalidRank();
      cr = getInvalidRank();
      if (nl > 0) cl = node + 1;        // left child
      if (nr > 0) cr = cl + static_cast<int>(nl);         // right child

      if (node == rank) break;
      else {
         TBOX_ASSERT(nl > 0);
         TBOX_ASSERT(cl != getInvalidRank());
         upp = up;
         up = node;
         ++d_generation;
         if (nr < 1 || rank < cr) {
            rbeg = cl;
            rend = cl + static_cast<int>(nl) - 1;
            parents_child_number = d_child_number;
            d_child_number = 0;
         } else {
            TBOX_ASSERT(nr > 0);
            TBOX_ASSERT(cr != getInvalidRank());
            rbeg = cr;
            rend = cr + static_cast<int>(nr) - 1;
            parents_child_number = d_child_number;
            d_child_number = 1;
         }
      }
   }

   const int gparent = upp;
   d_rank = rank;
   d_parent = up;
   d_children[0] = cl;
   d_children[1] = cr;
   d_num_children = 0;

   if (d_do_left_leaf_switch) {
      if (is_switchable) {
         /*
          * Trios of a parent and 2 leaf children are subject to
          * switching, in which the parent and left child switch
          * places:
          *
          * Before:    parent              After:   left
          *            /  \                        /  \
          *        left    right             parent    right
          */
         if (nl == 1) {
            // This is a parent in a left-leaf switchable.
            d_parent = cl;
            d_children[0] = getInvalidRank();
            d_children[1] = getInvalidRank();
            ++d_generation;
            d_child_number = 0;
         } else if (rank == d_parent + 1) {
            // This is a left child in a left-leaf switchable.
            d_children[0] = d_parent;
            d_parent = gparent;
            d_children[1] = rank + 1;
            --d_generation;
            d_child_number = parents_child_number;
         } else {
            // This is a right child in a left-leaf switchable.
            d_parent = d_parent + 1;
         }
         if (last_rank - first_rank + 1 == 3) {
            // Special case of exactly 3 ranks allows the root be switched.
            d_root_rank = first_rank + 1;
         }
      } else {
         /*
          * Rank is not in a switchable trio, but its children
          * may be.  Example:
          *
          * Before:      rank                   After:       rank
          *             /    \                              /    \
          *            /      \                            /      \
          *      rank+1        rank+4                rank+2        rank+5
          *     /      \      /      \              /      \      /      \
          *    /        \    /        \            /        \    /        \
          * rank+2   rank+3  rank+5   rank+6    rank+1   rank+3  rank+4   rank+6
          */
         if (nl == 3) {
            ++d_children[0];
         }
         if (nr == 3) {
            ++d_children[1];
         }
      }
   }

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
