/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_BalancedDepthFirstTree
#define included_tbox_BalancedDepthFirstTree

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/RankTreeStrategy.h"

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Implementation of RankTreeStrategy aranging ranks in a
 * balanced, depth-first tree.
 *
 * An example of a tree created is
 *
 * @verbatim
 *                     0
 *                    / \
 *                   /   \
 *                  /     \
 *                 1       8
 *                / \     / \
 *               2   5   9   12
 *              /|  /|  / \  |\
 *             3 4 6 7 10 11 13...
 * @endverbatim
 *
 * The tree is as balanced as possible.  Nodes that are close together
 * in the tree tends to be close together in natural ordering.  Without
 * knowing about the underlying message passing network structure, we
 * assume that close natural ordering usually means close together on
 * the network.  Thus nodes close together in the tree are also close
 * together on the network.  Thus, communication between nearest
 * neighbors in the tree tend to be faster.
 *
 * The tree formed by this class has the property that each and every
 * subtree is composed nodes with contiguous natural ordering.  This
 * again benefits communication.
 */
class BalancedDepthFirstTree:public RankTreeStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   BalancedDepthFirstTree();

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] first_rank
    * @param[in] last_rank
    * @param[in] rank
    * @param[in] do_left_leaf_switch See setLeftLeafSwitching()
    */
   BalancedDepthFirstTree(
      int first_rank,
      int last_rank,
      int rank,
      bool do_left_leaf_switch = true);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~BalancedDepthFirstTree();

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
    *
    * @pre (first_rank <= rank) && (rank <= last_rank)
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
   getParentRank() const
   {
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
   unsigned int
   getNumberOfChildren() const {
      return d_num_children;
   }

   /*!
    * @brief Return the child number, or invalidChildNumber() if is
    * root of the tree.
    */
   unsigned int
   getChildNumber() const {
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
   unsigned int
   getGenerationNumber() const {
      return d_generation;
   }

   /*!
    * @brief Return the rank of the root of the tree.
    */
   int getRootRank() const {
      return d_root_rank;
   }

   /*!
    * @brief Whether to do left-leaf-switch, which puts all leaves
    * within one rank of their parents.
    *
    * Without left-leaf-switching, about half of the leaves are within
    * one rank of their parents.  The other half are within 2 ranks of
    * their parents.  On the downside, switching moves the leaves'
    * parents a little farther away from the grandparents.
    *
    * Default choice is true.  To change the choice, this call must
    * be made before setupTree().
    */
   void
   setLeftLeafSwitching(bool do_left_leaf_switch) {
      TBOX_ASSERT(d_rank == getInvalidRank());
      d_do_left_leaf_switch = do_left_leaf_switch;
   }

private:
   // Unimplemented copy constructor.
   BalancedDepthFirstTree(
      const BalancedDepthFirstTree& other);

   // Unimplemented assignment operator.
   BalancedDepthFirstTree&
   operator = (
      const BalancedDepthFirstTree& rhs);

   /*!
    * @brief Set up the tree.
    *
    * Setting up has log complexity.
    *
    * @param first_rank The first in a contiguous range of ranks in the
    * communication group.
    *
    * @param last_rank The last in a contiguous range of ranks in the
    * communication group.
    *
    * @param rank The rank whose parent and children are sought.
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

   bool d_do_left_leaf_switch;

};

}
}

#endif  // included_tbox_BalancedDepthFirstTree
