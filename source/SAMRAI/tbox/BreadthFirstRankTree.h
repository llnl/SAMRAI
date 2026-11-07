/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_BreadthFirstRankTree
#define included_tbox_BreadthFirstRankTree

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Implementation of RankTreeStrategy aranging ranks in a
 * breadth-first tree.
 *
 * An example of a binary tree created is
 *
 * @verbatim
 *                     0
 *                    / \
 *                   /   \
 *                  /     \
 *                 1       2
 *                / \     / \
 *               3   4   5   6
 *              /|  /|   |\   \
 *             7 8 9 10 11 12 13...
 * @endverbatim
 *
 * The degree of the tree (max number of children at each node) can be
 * set using setTreeDegree(), but that must be done before the tree is
 * set up.  The default is a binary tree (two children).
 *
 * A property of the tree is that the first rank, i0, is always the
 * root.  It is a d-degree tree.  The children of rank i are rank
 * d*(i-i0) and d*(i-i0)+1.  It's parent is rank (i-i0-1)/d.  The tree
 * is as balanced as possible.  Node ranks are close to their siblings
 * but not their parents.  This arangement does not map well to any
 * known communication networks, so this tree does not perform well.
 * It is included for comparison.
 */
class BreadthFirstRankTree:public RankTreeStrategy
{

public:
   typedef int LocalId;

   /*!
    * @brief Constructor.
    */
   BreadthFirstRankTree();

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] first_rank
    * @param[in] last_rank
    * @param[in] rank
    * @param[in] degree See setTreeDegree()
    */
   BreadthFirstRankTree(
      int first_rank,
      int last_rank,
      int rank,
      unsigned int degree = 2);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~BreadthFirstRankTree();

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
   getRank() const
   {
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
    *
    * @param [in] child_number
    * a binary tree.
    */
   int
   getChildRank(
      unsigned int child_number) const
   {
      return (child_number < d_num_children) ?
             d_children[child_number] : getInvalidRank();
   }

   unsigned int
   getNumberOfChildren() const
   {
      return d_num_children;
   }

   /*!
    * @brief Return the child number, or invalidChildNumber() if is
    * root of the tree.
    */
   unsigned int getChildNumber() const
   {
      return d_child_number;
   }

   /*!
    * @brief Return the degree of the tree (the maximum number of
    * children each node may have).
    */
   unsigned int getDegree() const
   {
      return static_cast<unsigned int>(d_children.size());
   }

   /*!
    * @brief Return the generation number.
    */
   unsigned int getGenerationNumber() const
   {
      return d_generation;
   }

   /*!
    * @brief Return the rank of the root of the tree.
    */
   int getRootRank() const
   {
      return d_root_rank;
   }

   /*!
    * @brief Set the degree (max number of children) of the tree.
    *
    * Default choice is 2 (binary tree).  To change the choice, this
    * call must be made before setupTree().
    */
   void setTreeDegree(unsigned int tree_degree)
   {
      TBOX_ASSERT(d_rank == getInvalidRank());
      d_children.resize(tree_degree, getInvalidRank());
   }

private:
   // Unimplemented copy constructor.
   BreadthFirstRankTree(
      const BreadthFirstRankTree& other);

   // Unimplemented assignment operator.
   BreadthFirstRankTree&
   operator = (
      const BreadthFirstRankTree& rhs);

   /*!
    * @brief Construct the tree.
    *
    * Setting up has complexity 1.
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
    *
    * @see setupTree();
    */
   int d_rank;

   int d_parent;

   /*!
    * @brief Children rank.  Length of this member is also the degree
    * of the tree.
    *
    * Number of valid child ranks is d_child_number.  The rest of the
    * entries should be invalid ranks.
    */
   std::vector<int> d_children;

   unsigned int d_num_children;

   unsigned int d_child_number;

   unsigned int d_generation;

   int d_root_rank;

};

}
}

#endif  // included_tbox_BreadthFirstRankTree
