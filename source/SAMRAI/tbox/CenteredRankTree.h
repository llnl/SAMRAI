/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_CenteredRankTree
#define included_tbox_CenteredRankTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/RankTreeStrategy.h"

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Implementation of RankTreeStrategy aranging ranks with the
 * root at the center of its descendent ranks.
 *
 * An example of a tree created is
 *
 * @verbatim
 *                     7
 *                    / \
 *                   /   \
 *                  /     \
 *                 3       11
 *                / \     / \
 *               1   5   9   13
 *              /|  /|  /|   |\
 *             0 2 4 6 8 10 12 14
 *
 * or with 0 at the root:
 *                     0
 *                    / \
 *                   /   \
 *                  /     \
 *                 4       11
 *                / \     / \
 *               2   6   9   13
 *              /|  /|  /|   |\
 *             1 3 5 7 8 10 12 14
 * @endverbatim
 *
 * The root of tree is the center, rounded up, of the rank in it left
 * and right branches.  The exception is the case where we force the
 * first rank to be the root of the entire tree (@see
 * makeFirstRankTheRoot()).
 *
 * Nodes that are close together
 * in the tree tends to be close together in natural ordering.  Without
 * knowing about the underlying message passing network structure, we
 * assume that close natural ordering usually means close together on
 * the network.  Thus nodes close together in the tree are also close
 * together on the network.  Thus, communication between nearest
 * neighbors in the tree tend to be faster.  If the weight of each edge
 * is the difference in ranks of its nodes, the sum of all weights in
 * a tree of N nodes is N*ln(N).
 *
 * The tree formed by this class has the property that each and every
 * subtree is composed nodes with contiguous natural ordering.  This
 * again benefits communication.
 */
class CenteredRankTree:public RankTreeStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   CenteredRankTree();

   /*!
    * @brief Initializing constructor.
    *
    * @see setupTree()
    *
    * @param[in] mpi
    * @param[in] make_first_rank_the_root See makeFirstRankTheRoot()
    */
   CenteredRankTree(
      const SAMRAI_MPI& mpi,
      bool make_first_rank_the_root = true);

   /*!
    * @brief Initializing constructor.
    *
    * @see setupTree()
    *
    * @param[in] first_rank
    * @param[in] last_rank
    * @param[in] rank
    * @param[in] make_first_rank_the_root See makeFirstRankTheRoot()
    */
   CenteredRankTree(
      int first_rank,
      int last_rank,
      int rank,
      bool make_first_rank_the_root = true);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~CenteredRankTree();

   /*!
    * @brief Set up the tree.
    *
    * Set up the tree for the processors in the given RankGroup.
    * Prepare to provide tree data for the given rank.
    *
    * @param[in] rank_group
    *
    * @param[in] my_rank The rank whose parent and children are
    * sought, usually the local process.
    */
   void
   setupTree(
      const RankGroup& rank_group,
      int my_rank);

   /*!
    * @brief Access the rank used to initialize.
    */
   int
   getRank() const {
      return d_rank;
   }

   /*!
    * @brief Access the parent rank.
    */
   int
   getParentRank() const {
      return d_parent;
   }

   /*!
    * @brief Access a child rank.
    */
   int
   getChildRank(unsigned int child_number) const {
      return (child_number < d_num_children) ?
             d_children[child_number] : getInvalidRank();
   }

   /*!
    * @brief Return the number of children.
    */
   unsigned int getNumberOfChildren() const {
      return d_num_children;
   }

   /*!
    * @brief Return the child number, or invalidChildNumber() if is
    * root of the tree.
    */
   unsigned int getChildNumber() const {
      return d_child_number;
   }

   /*!
    * @brief Return the degree of the tree (the maximum number of
    * children each node may have).
    */
   unsigned int getDegree() const {
      return 2;
   }

   /*!
    * @brief Return the generation number.
    */
   unsigned int getGenerationNumber() const {
      return d_generation;
   }

   /*!
    * @brief Return the rank of the root of the tree.
    */
   int getRootRank() const {
      return d_root_rank;
   }

   /*!
    * @brief Set whether to use the first rank in the range for the
    * root instead of the rank at the middle of the range.
    *
    * Default choice is false.  To change the choice, this call must
    * be made before setupTree().
    */
   void makeFirstRankTheRoot(bool make_first_rank_the_root) {
      TBOX_ASSERT(d_rank == getInvalidRank());
      d_make_first_rank_the_root = make_first_rank_the_root;
   }

private:
   // Unimplemented copy constructor.
   CenteredRankTree(
      const CenteredRankTree& other);

   // Unimplemented assignment operator.
   CenteredRankTree&
   operator = (
      const CenteredRankTree& rhs);

   /*!
    * @brief Set up the tree.
    *
    * Setting up has log complexity.
    *
    * @param[in] first_rank The first in a contiguous range of ranks in
    * the communication group.
    *
    * @param[in] last_rank The last in a contiguous range of ranks in
    * the communication group.
    *
    * @param[in] rank The rank whose parent and children are sought.
    */
   void
   setupTreeForContiguousRanks(
      int first_rank,
      int last_rank,
      int rank);

   /*!
    * @brief Initialized rank.
    */
   int d_rank;

   int d_parent;

   int d_children[2];

   int d_root_rank;

   unsigned int d_num_children;

   unsigned int d_child_number;

   unsigned int d_generation;

   bool d_make_first_rank_the_root;

};

}
}

#endif  // included_tbox_CenteredRankTree
