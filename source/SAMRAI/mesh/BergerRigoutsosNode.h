/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Node in asynchronous Berger-Rigoutsos tree
 *
 ************************************************************************/
#ifndef included_mesh_BergerRigoutsosNode
#define included_mesh_BergerRigoutsosNode

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/PatchLevel.h"

#include <set>
#include <list>
#include <vector>
#include <algorithm>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Node in the asynchronous Berger-Rigoutsos (BR) tree.
 * Do not directly use this class; for clustering, use BergerRigoutsos
 * instead.
 *
 * In mesh generation, the BR algorithm can be used to cluster
 * tagged cells into boxes.
 * This algorithm is described in Berger and Rigoutsos,
 * IEEE Trans. on Sys, Man, and Cyber (21)5:1278-1286.
 *
 * Each tree node is associated with a candidate box,
 * an owner process coordinating distributed computations on the box
 * and a group of processors participating in those computations.
 * Should the candidate box be one of the final output boxes,
 * the owner also owns the graph node associated with the box.
 *
 * To use this class:
 * -# Construct the root node, an object of type
 *    BergerRigoutsosNode.
 * -# Finetune the algorithm settings using the methods under
 *    "Algorithm settings".
 * -# Start clustering by calling clusterAndComputeRelationships().
 *
 * The 2 primary outputs of this implementation are:
 * -# A BoxLevel of Boxes containing input tags.  Each node
 *    corresponds to an output box.
 * -# Connector between the tag BoxLevel and the new BoxLevel.
 *
 * TODO:
 * -# Implement MOST_TAGS ownership option.  This may be an
 *    improvement over the MOST_OVERLAP and is easy to do
 *    because the number of local tags in the candidate box
 *    is already computed.
 */

class BergerRigoutsosNode:
   private tbox::AsyncCommStage::Handler
{

public:
   /*!
    * @brief Destructor.
    */
   virtual ~BergerRigoutsosNode();

private:
   /*
    * BergerRigoutsos and BergerRigoutsosNode are tightly coupled.
    * Technically, BergerRigoutsosNode can be made a private subclass
    * of BergerRigoutsos.  BergerRigoutsos has the common parts of the
    * data and algorithm.  BergerRigoutsosNode has the node-specific
    * parts.
    */
   friend class BergerRigoutsos;

   /*!
    * @brief Construct a root node for a single block.
    *
    * @param common_params  Parameters shares by all nodes in clustering
    * @param box            Global bounding box for a single block
    */
   BergerRigoutsosNode(
      BergerRigoutsos* common_params,
      const hier::Box& box);

   const tbox::Dimension& getDim() const {
      return d_box.getDim();
   }

   /*
    * Static integer constant defining value corresponding to a bad integer.
    */
   static const int BAD_INTEGER;

   /*!
    * @brief Shorthand for std::vector<int> for internal use.
    */
   typedef std::vector<int> VectorOfInts;

   /*!
    * @brief Construct a non-root node.
    *
    * This is private because the object requires setting up
    * after constructing.  Nodes constructed this way are
    * only meant for internal use by the recursion mechanism.
    */
   BergerRigoutsosNode(
      BergerRigoutsos* common_params,
      BergerRigoutsosNode* parent,
      const int child_number);

   /*!
    * @brief Names of algorithmic phases while outside of
    * continueAlgorithm().
    *
    * "For_data_only" phase is when the node is only used to
    * store data. If the node is to be executed, it enters the
    * "to_be_launched" phase.
    *
    * All names beginning with "reduce", "gather" or "bcast"
    * refer to communication phases, where control is
    * returned before the algorithm completes.
    *
    * The "children" phase does not explicitly contain communication,
    * but the children may perform communication.
    *
    * The "completed" phase is when the algorithm has run to completion.
    * This is where the recursive implementation would return.
    *
    * The "deallocated" phase is for debugging.  This phase is
    * set by the destructor, just to help find nodes that
    * are deallocated but somehow was referenced.
    */
   enum WaitPhase { for_data_only,
                    to_be_launched,
                    reduce_histogram,
                    bcast_acceptability,
                    gather_grouping_criteria,
                    bcast_child_groups,
                    run_children,
                    bcast_to_dropouts,
                    completed,
                    deallocated };

   /*!
    * @brief MPI tags identifying messages.
    *
    * Each message tag is the d_mpi_tag plus a PhaseTag.
    * Originally, there were different tags for different
    * communication phases, determined by d_mpi_tag plus
    * a PhaseTag.  But this is not really needed,
    * so all phases use the tag d_mpi_tag.  The PhaseTag
    * type is just here in case we have to go back to using
    * them.
    */
   enum PhaseTag { reduce_histogram_tag = 0,
                   bcast_acceptability_tag = 0,
                   gather_grouping_criteria_tag = 0,
                   bcast_child_groups_tag = 0,
                   bcast_to_dropouts_tag = 0,
                   total_phase_tags = 1 };

   /*!
    * @brief Continue the the BR algorithm.
    *
    * Parameters for finding boxes are internal.
    * They should be set in the constructor.
    *
    * In parallel, this the method may return before
    * algorithm is completed.  In serial, no communication
    * is done, so the algorithm IS completed when this
    * method returns.  The method is completed if it
    * returns WaitPhase::completed.  This method may
    * and @em should be called multiple times as long as
    * the algorithm has not completed.
    *
    * If this method returns before the algorithm is
    * complete, this object will have put itself on
    * the leaf queue to be checked for completion later.
    *
    * @return The communication phase currently running.
    *
    * @pre (d_parent == 0) || (d_parent->d_wait_phase != completed)
    * @pre inRelaunchQueue(this) == d_common->d_relaunch_queue.end()
    */
   WaitPhase
   continueAlgorithm();

   /*!
    * @brief Candidate box acceptance state.
    *
    * Note that accepted values are odd and rejected
    * and undetermined values are even!  See boxAccepted(),
    * boxRejected() and boxHasNoTag().
    *
    * It is not critical to have all values shown,
    * but the values help in debugging.
    *
    * Meaning of values:
    * - "hasnotag_by_owner": histogram is truly empty (after sum reduction).
    *   We don't accept the box, but we don't split it either.
    *   (This can only happen at the root node, as child
    *   boxes are guaranteed to have tags.)
    * - "(rejected|accepted)_by_calculation": decision by calculation
    *   on the owner process.
    * - "(rejected|accepted)_by_owner": decision by owner process,
    *   broadcast to participants.
    * - "(rejected|accepted)_by_recombination": decision by recombination
    *   on local process.
    * - "(rejected|accepted)_by_dropout_bcast": decision by participant group,
    *   broadcast
    *    to the dropout group.
    */
   enum BoxAcceptance { undetermined = -2,
                        hasnotag_by_owner = -1,
                        rejected_by_calculation = 0,
                        accepted_by_calculation = 1,
                        rejected_by_owner = 2,
                        accepted_by_owner = 3,
                        rejected_by_recombination = 4,
                        accepted_by_recombination = 5,
                        rejected_by_dropout_bcast = 6,
                        accepted_by_dropout_bcast = 7 };

   //@{
   //! @name Delegated tasks for various phases of running algorithm.
   void
   makeLocalTagHistogram();

   void
   reduceHistogram_start();

   bool
   reduceHistogram_check();

   void
   computeMinimalBoundingBoxForTags();

   void
   acceptOrSplitBox();

   void
   broadcastAcceptability_start();

   bool
   broadcastAcceptability_check();

   void
   countOverlapWithLocalPatches();

   void
   gatherGroupingCriteria_start();

   bool
   gatherGroupingCriteria_check()
   {
      if (d_group.size() == 1) {
         return true;
      }
      d_comm_group->checkGather();
      /*
       * Do nothing yet with the overlap data d_recv_msg.
       * We extract it in formChildGroups().
       */
      return d_comm_group->isDone();
   }

   //! @brief Form child groups from gathered overlap counts.
   // @pre d_common->d_rank == d_owner
   // @pre d_recv_msg.size() == 4 * d_group.size()
   void
   formChildGroups();

   //! @brief Form child groups from local copy of all level boxes.
   void
   broadcastChildGroups_start();

   bool
   broadcastChildGroups_check();

   void
   runChildren_start();

   bool
   runChildren_check();

   void
   broadcastToDropouts_start();

   bool
   broadcastToDropouts_check();

   void
   createBox();

   void
   eraseBox();

   //! @brief Compute new graph relationships touching local tag nodes.
   // @pre d_common->d_compute_relationships > 0
   // @pre boxAccepted()
   // @pre d_box_acceptance != accepted_by_dropout_bcast
   // @pre (d_parent == 0) || (d_box.numberCells() >= d_common->d_min_box)
   // @pre d_box_acceptance != accepted_by_dropout_bcast
   void
   computeNewNeighborhoodSets();
   //@}

   //@{
   //! @name Utilities for implementing algorithm

   //! @brief Find the index of the owner in the group.
   int
   findOwnerInGroup(
      int owner,
      const VectorOfInts& group) const
   {
      for (unsigned int i = 0; i < group.size(); ++i) {
         if (group[i] == owner) {
            return i;
         }
      }
      return -1;
   }

   //! @brief Claim a unique tag from process's available tag pool.
   // @pre d_mpi_tag < 0
   void
   claimMPITag();

   /*!
    * @brief Heuristically determine "best" tree degree for
    * communication group size.
    */
   int
   computeCommunicationTreeDegree(
      int group_size) const
   {
      int tree_deg = 2;
      int shifted_size = group_size >> 3;
      while (shifted_size > 0) {
         shifted_size >>= 3;
         ++tree_deg;
      }
      return tree_deg;
   }

   void
   computeGlobalTagDependentVariables();

   bool
   findZeroCutSwath(
      int& cut_lo,
      int& cut_hi,
      const tbox::Dimension::dir_t dim);

   void
   cutAtInflection(
      int& cut_pt,
      int& inflection,
      const tbox::Dimension::dir_t dim);

   int
   getHistogramBufferSize(
      const hier::Box& box) const
   {
      int size = box.numberCells(0);
      int dim_val = d_common->getDim().getValue();
      for (tbox::Dimension::dir_t d = 1; d < dim_val; ++d) {
         size += box.numberCells(d);
      }
      return size;
   }

   int *
   putHistogramToBuffer(
      int* buffer);

   int *
   getHistogramFromBuffer(
      int* buffer);

   int *
   putBoxToBuffer(
      const hier::Box& box,
      int* buffer) const;

   int *
   getBoxFromBuffer(
      hier::Box& box,
      int* buffer) const;

   //! @brief Compute list of non-participating processes.
   // @pre main_group.size() >= sub_group.size()
   void
   computeDropoutGroup(
      const VectorOfInts& main_group,
      const VectorOfInts& sub_group,
      VectorOfInts& dropouts,
      const int add_group) const;

   BoxAcceptance
   intToBoxAcceptance(
      int i) const;

   bool
   boxAccepted() const
   {
      return bool(d_box_acceptance >= 0 && d_box_acceptance % 2);
   }

   bool
   boxRejected() const
   {
      return bool(d_box_acceptance >= 0 && d_box_acceptance % 2 == 0);
   }

   bool
   boxHasNoTag() const
   {
      return bool(d_box_acceptance == -1);
   }
   //@}

   //@{
   //! @name Utilities to help analysis and debugging
   std::list<BergerRigoutsosNode *>::const_iterator
   inRelaunchQueue(
      BergerRigoutsosNode* node_ptr) const
   {
      std::list<BergerRigoutsosNode *>::const_iterator li =
         std::find(d_common->d_relaunch_queue.begin(),
            d_common->d_relaunch_queue.end(),
            node_ptr);
      return li;
   }

   bool
   inGroup(
      VectorOfInts& group,
      int rank = -1) const
   {
      if (rank < 0) {
         rank = d_common->d_mpi.getRank();
      }
      for (size_t i = 0; i < group.size(); ++i) {
         if (rank == group[i]) {
            return true;
         }
      }
      return false;
   }

   //! @name Developer's methods for analysis and debugging this class.
   void
   printNodeState(
      std::ostream& co) const;

   /*!
    * @brief Unique id in the binary tree.
    *
    * - To have succinct formula, the root node has d_pos of 1.
    * - Parent id is d_pos/2
    * - Left child id is 2*d_pos
    * - Right child id is 2*d_pos+1
    * - Generation number is ln(d_pos)
    *
    * This parameter is only used for debugging.
    *
    * The id of a node grows exponentially with each generation.
    * If the position in the binary tree is too big to be represented
    * by an integer, d_pos is set to -1 for a left child and -2 for a
    * right child.
    */
   const int d_pos;

   /*!
    * @brief Common parameters shared with descendents and ancestors.
    *
    * Only the root of the tree allocates the common parameters.
    * For all others, this pointer is set by the parent.
    */
   BergerRigoutsos* d_common;

   //@{
   /*!
    * @name Tree-related data
    */

   //! @brief Parent node (or NULL for the root node).
   BergerRigoutsosNode* d_parent;

   //! @brief Left child.
   BergerRigoutsosNode* d_lft_child;

   //! @brief Right child.
   BergerRigoutsosNode* d_rht_child;

   //@}

   //@{
   /*!
    * @name Data for one recursion of the BR algorithm
    */

   /*
    * These parameters are listed roughly in order of usage.
    */

   hier::Box d_box;

   /*!
    * @name Id of participating processes.
    */
   VectorOfInts d_group;

   /*!
    * Minimum size of a Box that d_box can potentially be chopped into.
    */
   hier::IntVector d_min_box_size;

   /*!
    * Requested minimum size for total cell count of an accepted box.
    * This minimum may not be enforced in all cases.
    */
   size_t d_min_cell_request;

   /*!
    * @brief MPI tag for message within a node.
    *
    * The tag is determined by on the process that owns the parent
    * when the parent decides to split its box.  The tags are broadcasted
    * along with the children boxes.
    */
   int d_mpi_tag;

   /*!
    * @brief Overlap count with d_box.
    */
   size_t d_overlap;

   /*!
    * @brief Whether and how box is accepted.
    *
    * @see BoxAcceptance.
    */
   BoxAcceptance d_box_acceptance;

   /*!
    * @brief Histogram for all directions of box d_box.
    *
    * If local process is owner, this is initially the
    * local histogram, then later, the reduced histogram.
    * If not, it is just the local histogram.
    */
   VectorOfInts d_histogram[SAMRAI::MAX_DIM_VAL];

   /*!
    * @brief Number of tags in the candidate box.
    */
   int d_num_tags;

   /*!
    * @brief Box iterator corresponding to an accepted box on
    * the owner.
    *
    * This is relevant only on the owner, where the d_box is
    * in a container.  On contributors, the graph node is non-local
    * and stands alone.
    */
   hier::BoxContainer::const_iterator d_box_iterator;

   /*!
    * @brief Name of wait phase when continueAlgorithm()
    * exits before completion.
    */
   WaitPhase d_wait_phase;

   //@}

   //@{
   /*!
    * @name Lower-level parameters for communication.
    */

   //! @brief Buffer for organizing outgoing data.
   VectorOfInts d_send_msg;
   //! @brief Buffer for organizing incoming data.
   VectorOfInts d_recv_msg;

   tbox::AsyncCommGroup* d_comm_group;
   //@}

   //@{
   //! @name Deubgging aid

   /*!
    * @brief Generation number.
    *
    * The generation number is the parent's generation number plus 1.
    * The root has generation number 1.
    */
   const int d_generation;

   //! @brief Number of times continueAlgorithm was called.
   int d_n_cont;

   //@}
};

}
}

#endif  // included_mesh_BergerRigoutsosNode
