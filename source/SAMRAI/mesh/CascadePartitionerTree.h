/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using a "cascade" algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitionerTree
#define included_mesh_CascadePartitionerTree

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

#include <memory>

namespace SAMRAI {
namespace mesh {

class CascadePartitioner;

/*!
 * @brief A binary-tree of process groups in the CascadePartitioner
 * algorithm.  This class is used internally in the CascadePartioner
 * class.  It should not be used otherwise.
 *
 * In cascade partitioner, the MPI ranks are recursively split into
 * groups, forming the nodes of a binary tree.  The root branch
 * contains all ranks.  The leaves are single-process groups.  Each
 * group is represented by a CascadePartitionerTree.
 *
 * @b Terminology: The root group includes all ranks.  It splits into
 * lower and upper groups (also called branches).  The lower branch
 * has the lower ranks.  If the group has an odd number of ranks, the
 * upper branch has one rank more than the lower.  The branch
 * containing the local process is also known as the "near branch";
 * the one not containing the local process is the "far branch".
 */
class CascadePartitionerTree
{

public:
   /*!
    * @brief Construct the tree for the given CascadePartitioner by
    * constructing the root and recursively constructing its children.
    */
   CascadePartitionerTree(
      const CascadePartitioner& partitioner);

   ~CascadePartitionerTree();

   /*!
    * @brief Distribute the load using the cascade algorithm.
    */
   void
   distributeLoad();

   void
   printClassData(
      std::ostream& co,
      const std::string& border) const;

   //! @brief Generation number (generation 0 contains all ranks).
   int generationNum() const {
      return d_gen_num;
   }

   //! @brief Size of group (number of processes in it).
   size_t size() const {
      return d_end - d_begin;
   }

   //! @brief Whether group contains a given rank.
   bool containsRank(int rank) const {
      return d_begin <= rank && rank < d_end;
   }

private:
   //! @brief Where a group falls in the next larger group.
   enum Position { Lower = 0, Upper = 1 };

   /*
    * Static integer constants.  Tags are for distinguishing messages
    * from different phases of the algorithm.
    */
   static const int CascadePartitionerTree_TAG_InfoExchange0 = 1000;
   static const int CascadePartitionerTree_TAG_InfoExchange1 = 1001;
   static const int CascadePartitionerTree_TAG_LoadTransfer0 = 1002;
   static const int CascadePartitionerTree_TAG_LoadTransfer1 = 1003;

   /*!
    * @brief Construct child node based on its position in the parent.
    *
    * @param parent
    * @param group_position Position of this group in its parent.
    */
   CascadePartitionerTree(
      CascadePartitionerTree& parent,
      Position group_position);

   /*!
    * @brief Allocate and set up the group's children, if any.
    */
   void
   makeChildren();

   /*!
    * @brief Combine near and far data for children group to compute
    * work-related data for this group.
    */
   void
   combineChildren();

   /*!
    * @brief Improve balance of the children of this group by
    * supplying load from overloaded child to underloaded child.
    */
   void
   balanceChildren();

   //! @brief Estimated surplus of the group.
   double estimatedSurplus() const {
      return d_work - d_obligation;
   }

   /*!
    * @brief Try to supply the requested amount of work by removing
    * it from this group, and return the (estimated) amount supplied.
    *
    * The return value is an estimate of the supplied work, because
    * the group almost always contains remote ranks whose exact
    * actions are not known.  Single-process groups set aside any work
    * it personally gives up.
    *
    * @param work_requested
    *
    * @param taker Representative of the group getting this work.
    *
    * @return Estimate of the amount supplied based on available work
    * and assuming perfect load cutting.
    */
   double
   supplyWork(
      double work_requested,
      int taker);

   void
   sendShipment(
      int taker);
   void
   receiveAndUnpackSuppliedLoad();

   //! @brief Recompute work-related data for a leaf group.
   void
   recomputeLeafData();

   /*!
    * @brief Reset obligation recursively for all descendents, based
    * on given average load.
    */
   void
   resetObligation(
      double avg_load);

   /*!
    * @brief Compute interval for updating Connector during load
    * distribution.
    */
   double
   computeConnectorUpdateInterval() const;

   //! @brief Data the main CascadePartitioner shares with all parts of the tree.
   const CascadePartitioner* d_common;

   //@{
   //! @brief Group specification

   //! @brief Generation number.  (Generation 0 contains all ranks.)
   int d_gen_num;

   //! @brief First rank in group.
   int d_begin;

   //! @brief One past last rank in group.
   int d_end;

   /*!
    * @brief Rank of contacts in sibling branch.
    *
    * Communication between sibling groups occurs between processes in
    * one group and their respective contacts in the sibling group.
    *
    * Most processes have just one contact in the sibling group.  The
    * only processes to have 2 contacts are in the lower sibling of a
    * parent that has an odd number of ranks.  In these cases, the
    * upper sibling has one more rank than the lower.  The last rank
    * of the lower sibling contacts the last 2 ranks in the upper.
    */
   int d_contact[2];

   //! @brief Parent group.
   CascadePartitionerTree* d_parent;

   /*!
    * @brief Lower and upper children branches.
    *
    * Children are allocated iff this is a near group and not a leaf.
    */
   CascadePartitionerTree* d_children[2];

   //! @brief Near child branch (branch containing local process).
   CascadePartitionerTree* d_near;

   //! @brief Far child branch (branch not containing local process).
   CascadePartitionerTree* d_far;

   //! @brief Group containing just the local process.
   CascadePartitionerTree* d_leaf;

   //@}

   //@{
   //! @name Work measures

   //! @brief Estimated or actual amount of work in this branch.
   double d_work;

   //! @brief Amount of work the group obligated to have.
   double d_obligation;

   //@}

   //@{
   //! @name For determining participation in certain group activities.

   //! @brief Whether this group may supply work to its sibling.
   bool d_group_may_supply;

   /*!
    * @brief Whether local process (if this is a near group) or
    * contact processes (if this is a far group) may supply load.
    *
    * If this is a near group, first value specifies whether the local
    * process may supply work and second value is unused.  If this is
    * a far group, the two values correspond to whether the two
    * contacts (d_contact) may supply work.
    */
   bool d_process_may_supply[2];

   //@}

};

}
}

#endif
