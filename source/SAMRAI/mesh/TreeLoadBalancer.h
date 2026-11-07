/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer
#define included_mesh_TreeLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/CommGraphWriter.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/RankTreeStrategy.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Provides load balancing routines for AMR hierarchy by
 * implemementing the LoadBalancerStrategy.
 *
 * This class implements a tree-based load balancer.  The MPI
 * processes are arranged in a tree.  Work load is transmitted from
 * process to process along the edges of the tree.
 *
 * Currently, only uniform load balancing is supported.  Eventually,
 * non-uniform load balancing should be supported.  (Non-uniform load
 * balancing is supported by the CutAndPackLoadBalancer class.)
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 *   - \b flexible_load_tolerance
 *   Fraction of ideal load a process can
 *   take on in order to reduce box cutting and load movement.  Higher
 *   values often reduce partitioning time and box count but produce
 *   less balanced work loads.  Surplus work greater than this
 *   tolerance can still result due to other constraints, such as
 *   minimum box size.
 *
 *   - \b tile_size
 *   Tile size when using tile mode.  Tile mode restricts box cuts
 *   to tile boundaries.  Default is 1, which is equivalent to no restriction.
 *
 *   - \b max_spread_procs
 *   This parameter limits how many processes may receive the load of one
 *   process in a load distribution cycle.  If a process has too much
 *   initial load, this limit causes the load to distribute the load over
 *   multiple cycles.  It alleviates the bottle-neck of one process having
 *   to work with too many other processes in any cycle.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>flexible_load_tolerance</td>
 *     <td>double</td>
 *     <td>0.05</td>
 *     <td>0-1</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>tile_size</td>
 *     <td>IntVector</td>
 *     <td>1</td>
 *     <td>1-</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>max_spread_procs</td>
 *     <td>int</td>
 *     <td>500</td>
 *     <td> > 1</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.  Defaults listed
 * in parenthesis:
 *
 * @internal DEV_voucher_mode (false)
 * bool
 * Whether to use experimental voucher mode.
 *
 * @internal DEV_allow_box_breaking (true)
 * bool
 * Whether to allow box-breaking.  Set to false when boxes have
 * been pre-cut.
 *
 * @see LoadBalanceStrategy
 */

class TreeLoadBalancer:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if database provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined identifier used for diagnostic reports
    * and timer names.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    *
    * @param[in] rank_tree How to arange a contiguous range of MPI ranks
    * into a tree.  If omitted, we use a tbox::CenteredRankTree.
    *
    * @pre !name.empty()
    */
   TreeLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>(),
      const std::shared_ptr<tbox::RankTreeStrategy>& rank_tree =
         std::shared_ptr<tbox::RankTreeStrategy>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~TreeLoadBalancer();

   /*!
    * @brief Set the internal SAMRAI_MPI to a duplicate of the given
    * SAMRAI_MPI.
    *
    * The given SAMRAI_MPI must have a valid communicator.
    *
    * The given SAMRAI_MPI is duplicated for private use.  This
    * requires a global communication, so all processes in the
    * communicator must call it.  The advantage of a duplicate
    * communicator is that it ensures the communications for the
    * object won't accidentally interact with unrelated
    * communications.
    *
    * If the duplicate SAMRAI_MPI it is set, the TreeLoadBalancer will
    * only balance BoxLevels with congruent SAMRAI_MPI objects and
    * will use the duplicate SAMRAI_MPI for communications.
    * Otherwise, the SAMRAI_MPI of the BoxLevel will be used.  The
    * duplicate MPI communicator is freed when the object is
    * destructed, or freeMPICommunicator() is called.
    *
    * @pre samrai_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull
    */
   void
   setSAMRAI_MPI(
      const tbox::SAMRAI_MPI& samrai_mpi);

   /*!
    * @brief Free the internal MPI communicator, if any has been set.
    *
    * This is automatically done by the destructor, if needed.
    *
    * @see setSAMRAI_MPI().
    */
   void
   freeMPICommunicator();

   /*!
    * @brief Configure the load balancer to use the data stored
    * in the hierarchy at the specified descriptor index
    * for estimating the workload on each cell.
    *
    * Note: This method currently does not affect the results because
    * this class does not yet support uniform load balancing.
    *
    * @param data_id
    * Integer value of patch data identifier for workload
    * estimate on each cell.  An invalid value (i.e., < 0)
    * indicates that a spatially-uniform work estimate
    * will be used.  The default value is -1 (undefined)
    * implying the uniform work estimate.
    *
    * @param level_number
    * Optional integer number for level on which data id
    * is used.  If no value is given, the data will be
    * used for all levels.
    *
    * @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_id) is actually a  std::shared_ptr<pdat::CellDataFactory<double> >
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * @param[in] level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const;

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre !balance_to_anchor || balance_to_anchor->isTransposeOf(balance_to_anchor->getTranspose())
    * @pre (d_dim == balance_box_level.getDim()) &&
    *      (d_dim == min_size.getDim()) && (d_dim == max_size.getDim()) &&
    *      (d_dim == domain_box_level.getDim()) &&
    *      (d_dim == bad_interval.getDim()) && (d_dim == cut_factor.getDim())
    * @pre !hierarchy || (d_dim == hierarchy->getDim())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getSize())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getRank())
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * @brief Write out statistics recorded for the most recent load
    * balancing result.
    *
    * @param[in] output_stream
    */
   void
   printStatistics(
      std::ostream& output_stream = tbox::plog) const;

   /*!
    * @brief Enable or disable saving of tree data for diagnostics.
    *
    * @param [in] comm_graph_writer
    * External CommGraphWriter to save tree data to.
    * Use 0 to disable saving.
    */
   void
   setCommGraphWriter(
      const std::shared_ptr<tbox::CommGraphWriter>& comm_graph_writer)
   {
      d_comm_graph_writer = comm_graph_writer;
   }

   /*!
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   typedef double LoadType;

   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int TreeLoadBalancer_LOADTAG0 = 1;
   static const int TreeLoadBalancer_LOADTAG1 = 2;
   static const int TreeLoadBalancer_FIRSTDATALEN = 500;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   TreeLoadBalancer(
      const TreeLoadBalancer&);

   TreeLoadBalancer&
   operator = (
      const TreeLoadBalancer&);

   /*!
    * @brief Data to save for each sending/receiving process and the
    * branch at that process.
    *
    * Terminology: Parts of any tree may not be open to receiving work
    * from their parents because they have enough work already.  These
    * parts are not counted in the "effective" tree for the purpose of
    * distributing work on the branch.
    */
   class BranchData
   {

public:
      //! @brief Constructor.
      BranchData(
         const PartitioningParams& pparams,
         const TransitLoad& transit_load_prototype);
      //! @brief Copy constructor.
      BranchData(
         const BranchData& other);

      /*!
       * @brief Set the starting ideal, current and upper limit of the
       * load for the branch, which includes just the values from
       * local process.
       */
      void
      setStartingLoad(
         LoadType ideal,
         LoadType current,
         LoadType upperlimit);

      //! @brief Number of processes in branch.
      int numProcs() const
      {
         return d_num_procs;
      }
      //! @brief Number of processes in effective branch.
      int numProcsEffective() const
      {
         return d_eff_num_procs;
      }

      //@{
      //! @name Amount of work in branch, compared to various references.
      // surplus and deficit are current load compared to ideal.
      LoadType surplus() const
      {
         return d_branch_load_current - d_branch_load_ideal;
      }
      LoadType deficit() const
      {
         return d_branch_load_ideal - d_branch_load_current;
      }
      LoadType effSurplus() const
      {
         return d_eff_load_current - d_eff_load_ideal;
      }
      LoadType effDeficit() const
      {
         return d_eff_load_ideal - d_eff_load_current;
      }
      // excess and margin are current load compared to upper limit.
      LoadType excess() const
      {
         return d_branch_load_current - d_branch_load_upperlimit;
      }
      LoadType margin() const
      {
         return d_branch_load_upperlimit - d_branch_load_current;
      }
      LoadType effExcess() const
      {
         return d_eff_load_current - d_eff_load_upperlimit;
      }
      LoadType effMargin() const
      {
         return d_eff_load_upperlimit - d_eff_load_current;
      }
      //@}

      //! @brief Tell this tree to eventually ask for work from its parent.
      void setWantsWorkFromParent()
      {
         d_wants_work_from_parent = true;
      }

      //! @brief Get whether this branch want work from its parents.
      bool getWantsWorkFromParent() const
      {
         return d_wants_work_from_parent;
      }

      //@{
      //! @name Information on work shipped
      //! @brief Get amount of work shipped.
      LoadType getShipmentLoad() const
      {
         return d_shipment->getSumLoad();
      }
      //! @brief Get count of work shipped.
      size_t getShipmentPackageCount() const
      {
         return d_shipment->getNumberOfItems();
      }
      //! @brief Get count of originators of the work shipped.
      size_t getShipmentOriginatorCount() const
      {
         return d_shipment->getNumberOfOriginatingProcesses();
      }
      //@}

      //@{
      //! @name Methods supporting load import/export.
      /*!
       * @brief Adjust load to be sent away by taking work from or
       * dumping work into a reserve container.
       */
      LoadType
      adjustOutboundLoad(
         TransitLoad& reserve,
         LoadType ideal_load,
         LoadType low_load,
         LoadType high_load);

      //! @brief Move inbound load to the given reserve container.
      void
      moveInboundLoadToReserve(
         TransitLoad& reserve);

      /*!
       * @brief Incorporate child branch into this branch.
       */
      void
      incorporateChild(
         const BranchData& child);
      //@}

      //@{
      //! @name Packing/unpacking for communication up and down the tree.

      //! @brief Pack load/boxes for sending up to parent.
      void
      packDataToParent(
         tbox::MessageStream& msg) const;

      //! @brief Unpack load/boxes received from child.
      void
      unpackDataFromChild(
         tbox::MessageStream& msg);

      //! @brief Pack load/boxes for sending down to child.
      void
      packDataToChild(
         tbox::MessageStream& msg) const;

      //! @brief Unpack load/boxes received from parent.
      void
      unpackDataFromParentAndIncorporate(
         tbox::MessageStream& msg);

      //@}

      //! @brief Diagnostic printing.
      void
      recursivePrint(
         std::ostream& os,
         const std::string& border = std::string(),
         int detail_depth = 2) const;

      //! @brief Setup names of timers.
      void
      setTimerPrefix(
         const std::string& timer_prefix);

      //! @brief Whether to print steps for debugging.
      void setPrintSteps(bool print_steps)
      {
         d_print_steps = print_steps;
      }

private:
      /*!
       * @brief Number of processes in branch
       */
      int d_num_procs;

      /*!
       * @brief Current load in the branch, including local unassigned load.
       */
      LoadType d_branch_load_current;

      /*!
       * @brief Ideal load for the branch
       */
      LoadType d_branch_load_ideal;

      /*!
       * @brief Load the branch is willing to have, based on the load
       * tolerance and upper limits of children.
       */
      LoadType d_branch_load_upperlimit;

      /*!
       * @brief Number of processes in branch after pruning.
       */
      int d_eff_num_procs;

      /*!
       * @brief Current load in the effective branch.
       */
      LoadType d_eff_load_current;

      /*!
       * @brief Ideal load for the effective branch.
       */
      LoadType d_eff_load_ideal;

      /*!
       * @brief Load the effective branch is willing to have, which is
       * a sum of the upper limits of its effective children.
       */
      LoadType d_eff_load_upperlimit;

      /*!
       * @brief Load received or to be sent.
       *
       * If this object is for the local process, shipment is to or
       * from the process's *parent*.
       */
      std::shared_ptr<TransitLoad> d_shipment;

      /*!
       * @brief Whether branch expects its parent to send work down.
       */
      bool d_wants_work_from_parent;

      //! @brief Common partitioning parameters.
      const PartitioningParams* d_pparams;

      //@{
      //! @name Debugging and diagnostic data.
      std::shared_ptr<tbox::Timer> t_pack_load;
      std::shared_ptr<tbox::Timer> t_unpack_load;
      bool d_print_steps;
      //@}

   }; // BranchData declaration.

   /*!
    * @brief Check if there is any pending messages for the private
    * communication and throw an error if there is.
    */
   void
   assertNoMessageForPrivateCommunicator() const;

   /*!
    * Read parameters from input database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

   /*!
    * Utility functions to determine parameter values for level.
    */
   int
   getWorkloadDataId(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < static_cast<int>(d_workload_data_id.size()) ?
             d_workload_data_id[level_number] :
             d_master_workload_data_id;
   }

   /*!
    * @brief Compute the load for a Box.
    */
   double
   computeLoad(
      const hier::Box& box) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.
       */
      return double(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double
   computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.
       */
      return double((box * restriction).size());
   }

   /*
    * Count the local workload.
    */
   LoadType
   computeLocalLoad(
      const hier::BoxLevel& box_level) const;

   /*!
    * @brief Given an "unbalanced" BoxLevel, compute the BoxLevel that
    * is load-balanced within the given rank_group and update the
    * Connector between it and a reference BoxLevel.
    *
    * @pre !balance_to_reference || balance_to_reference->hasTranspose()
    * @pre d_dim == balance_box_level.getDim()
    */
   void
   loadBalanceWithinRankGroup(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_reference,
      const tbox::RankGroup& rank_group,
      const double group_sum_load) const;

   /*!
    * @brief Distribute load across the rank group using the tree
    * algorithm.
    *
    * Initial work is given in unbalanced_box_level.  Put the final
    * local work in balanced_work.
    */
   void
   distributeLoadAcrossRankGroup(
      TransitLoad& balanced_work,
      const hier::BoxLevel& unbalanced_box_level,
      const tbox::RankGroup& rank_group,
      double group_sum_load) const;

   /*!
    * @brief Compute surplus load per descendent who is still waiting
    * for load from parents.
    */
   LoadType
   computeSurplusPerEffectiveDescendent(
      const LoadType& unassigned_load,
      const LoadType& group_avg_load,
      const std::vector<BranchData>& child_branches,
      int first_child) const;

   /*!
    * @brief Create the cycle-based RankGroups the local process
    * belongs in.
    *
    * The RankGroup size increases exponentially with the cycle
    * number such that for the last cycle the rank group includes
    * all processes in d_mpi.
    *
    * @param [out] rank_group
    * @param [out] num_groups
    * @param [out] group_num
    * @param [in] cycle_fraction How far we are in the cycles.
    *   Value of 1 means the last cycle.
    */
   void
   createBalanceRankGroupBasedOnCycles(
      tbox::RankGroup& rank_group,
      int& num_groups,
      int& group_num,
      double cycle_fraction) const;

   /*!
    * @brief Set the AsyncCommPeer objects for this process to
    * communicate with its parent and children.
    *
    * @param [out] child_stage
    * @param [out] child_comms
    * @param [out] parent_stage
    * @param [out] parent_comm
    * @param [in] rank_group
    */
   void
   setupAsyncCommObjects(
      tbox::AsyncCommStage& child_stage,
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommStage& parent_stage,
      tbox::AsyncCommPeer<char> *& parent_comm,
      const tbox::RankGroup& rank_group) const;

   /*
    * @brief Undo the set-up done by setupAsyncCommObjects.
    *
    * @pre (d_mpi.getSize() != 1) || ((child_comms == 0) && (parent_comms == 0))
    */
   void
   destroyAsyncCommObjects(
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommPeer<char> *& parent_comm) const;

   /*!
    * @brief Set up timers for the object.
    */
   void
   setTimers();

   /*
    * Object dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   //! @brief Duplicated communicator object.  See setSAMRAI_MPI().
   mutable tbox::SAMRAI_MPI d_mpi;

   //! @brief Whether d_mpi is an internal duplicate.  See setSAMRAI_MPI().
   bool d_mpi_is_dupe;

   /*!
    * @brief Tile size, when restricting cuts to tile boundaries,
    * Set to 1 when not restricting.
    */
   hier::IntVector d_tile_size;

   //! @brief Max number of processes the a single process may spread load to per cycle.
   int d_max_spread_procs;

   //! @brief Whether to move load via vouchers.
   bool d_voucher_mode;

   //! @brief Whether to allow box breaking.
   bool d_allow_box_breaking;

   //! @brief How to arange a contiguous range of MPI ranks in a tree.
   const std::shared_ptr<tbox::RankTreeStrategy> d_rank_tree;

   /*!
    * @brief Utility to save data for communication graph output.
    */
   std::shared_ptr<tbox::CommGraphWriter> d_comm_graph_writer;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   std::vector<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Fraction of ideal load a process can accept over and above
    * the ideal.
    *
    * See input parameter "flexible_load_tolerance".
    */
   double d_flexible_load_tol;

   std::vector<double> d_artificial_minimum;

   /*!
    * @brief Metadata operations with timers set according to this object.
    */
   hier::MappingConnectorAlgorithm d_mca;

   //@{
   //! @name Data shared with private methods during balancing.
   mutable std::shared_ptr<PartitioningParams> d_pparams;
   mutable LoadType d_global_avg_load;
   //@}

   static const int s_default_data_id;

   //@{
   //! @name Used for evaluating peformance.

   bool d_barrier_before;
   bool d_barrier_after;

   /*!
    * @brief Whether to immediately report the results of the load
    * balancing cycles in the log files.
    */
   bool d_report_load_balance;

   /*!
    * @brief See "summarize_map" input parameter.
    */
   char d_summarize_map;

   /*
    * Performance timers.
    */
   std::shared_ptr<tbox::Timer> t_load_balance_box_level;
   std::shared_ptr<tbox::Timer> t_get_map;
   std::shared_ptr<tbox::Timer> t_use_map;
   std::shared_ptr<tbox::Timer> t_constrain_size;
   std::shared_ptr<tbox::Timer> t_distribute_load_across_rank_group;
   std::shared_ptr<tbox::Timer> t_compute_local_load;
   std::shared_ptr<tbox::Timer> t_compute_global_load;
   std::shared_ptr<tbox::Timer> t_compute_tree_load;
   std::vector<std::shared_ptr<tbox::Timer> > t_compute_tree_load_for_cycle;
   std::vector<std::shared_ptr<tbox::Timer> > t_load_balance_for_cycle;
   std::shared_ptr<tbox::Timer> t_send_load_to_children;
   std::shared_ptr<tbox::Timer> t_send_load_to_parent;
   std::shared_ptr<tbox::Timer> t_get_load_from_children;
   std::shared_ptr<tbox::Timer> t_get_load_from_parent;
   std::shared_ptr<tbox::Timer> t_post_load_distribution_barrier;
   std::shared_ptr<tbox::Timer> t_assign_to_local_and_populate_maps;
   std::shared_ptr<tbox::Timer> t_report_loads;
   std::shared_ptr<tbox::Timer> t_local_load_moves;
   std::shared_ptr<tbox::Timer> t_finish_sends;
   std::shared_ptr<tbox::Timer> t_barrier_before;
   std::shared_ptr<tbox::Timer> t_barrier_after;
   std::shared_ptr<tbox::Timer> t_child_send_wait;
   std::shared_ptr<tbox::Timer> t_child_recv_wait;
   std::shared_ptr<tbox::Timer> t_parent_send_wait;
   std::shared_ptr<tbox::Timer> t_parent_recv_wait;

   /*
    * Statistics on number of cells and patches generated.
    */
   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

   //@}

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_check_connectivity;
   char d_check_map;

};

}
}

#endif
