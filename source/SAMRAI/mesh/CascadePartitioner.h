/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner
#define included_mesh_CascadePartitioner

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/CascadePartitionerTree.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
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
 * implemementing the LoadBalancerStrategy using the cascade partitioning algorithm.
 *
 * The algorithm is described in the article "Advances in Patch-Based
 * Adaptive Mesh Refinement Scalability" submitted to JPDC.  Scaling
 * benchmark results are also in the article.
 *
 * This class can be used for both uniform or non-uniform load balancing.
 * To enable non-uniform load balancing, a call must be made to the method
 * setWorkloadPatchDataIndex to give this object a patch data id for
 * cell-centered workload data that must be set on the hierarchy outside of
 * this class.
 *
 * The default behavior of this class is to do uniform load balancing, treating
 * all cells of a level as having equal load value.
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
 *   process before updating Connectors.  If a process has too much
 *   initial load, this limit causes the Connector to be updated gradually,
 *   alleviating the bottle-neck of one process doing an excessive amount
 *   of Connector updates.
 *
 *   - \b use_vouchers
 *   Boolean parameter to turn on the optional voucher method for passing
 *   around workload during the cascade algorithm.  Note that non-uniform
 *   load balancing always uses the voucher method regardless of this
 *   parameter's value.
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
 *   <tr>
 *     <td>use_vouchers</td>
 *     <td>bool</td>
 *     <td>FALSE</td>
 *     <td>TRUE or FALSE</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @internal The following are developer inputs.  Defaults listed
 * in parenthesis:
 *
 * @internal DEV_reset_obligations (true)
 * bool
 * Whether to reset load obligations within groups that cannot change its load average.
 *
 * @internal DEV_limit_supply_to_surplus (true)
 * bool
 * Whether limit work a process can supply to its surplus.  The effects on partitioning
 * speed and quality are not yet known.
 *
 * @see LoadBalanceStrategy
 */

class CascadePartitioner:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if database provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined identifier used for error reporting
    * and timer names.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    *
    * @pre !name.empty()
    */
   CascadePartitioner(
      const tbox::Dimension& dim,
      const std::string& name,
      const std::shared_ptr<tbox::Database>& input_db =
         std::shared_ptr<tbox::Database>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~CascadePartitioner();

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
    * If the duplicate SAMRAI_MPI it is set, the CascadePartitioner will
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
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
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
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Configure the load balancer to use the data stored
    * in the hierarchy at the specified descriptor index
    * for estimating the workload on each cell.
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

private:
   typedef double LoadType;

   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int CascadePartitioner_LOADTAG0 = 1;
   static const int CascadePartitioner_LOADTAG1 = 2;
   static const int CascadePartitioner_FIRSTDATALEN = 500;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   CascadePartitioner(
      const CascadePartitioner&);

   void
   operator = (
      const CascadePartitioner&);

   /*
    * @brief Check if there is any pending messages for the private
    * communication and throw an error if there is.
    */
   void
   assertNoMessageForPrivateCommunicator() const;

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const std::shared_ptr<tbox::Database>& input_db);

   /*
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

   /*
    * Count the local workload.
    */
   LoadType
   computeLocalLoad(
      const hier::BoxLevel& box_level) const;

   /*
    * Compute the workload for the level based on a work function.
    */
   LoadType
   computeNonUniformWorkLoad(
      const hier::PatchLevel& patch_level) const;

   /*!
    * *@brief Implements the cascade partitioner algorithm.
    */
   void
   partitionByCascade(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_reference,
      bool use_vouchers = false) const;

   /*!
    * @brief Update Connectors balance_box_level<==>reference.
    */
   void
   updateConnectors() const;

   /*!
    * @brief Determine globally reduced work parameters.
    */
   void
   globalWorkReduction(
      LoadType local_work,
      bool has_any_load) const;

   //! @brief Compute log-base-2 of integer, rounded up.
   static int
   lgInt(
      int s);

   /*!
    * @brief Set up timers for the object.
    */
   void
   setTimers();

   /*
    * CascadePartitioner and CascadePartitionerTree are tightly
    * coupled.  CascadePartitioner has the common parts of the data
    * and algorithm.  CascadePartitionerTree has the group-specific
    * parts.  CascadePartitionerTree can be made a private subclass
    * of CascadePartitioner, but that would make a big file.
    */
   friend class CascadePartitionerTree;

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

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   std::vector<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Tile size, when restricting cuts to tile boundaries,
    * Set to 1 when not restricting.
    */
   hier::IntVector d_tile_size;

   /*!
    * @brief Max number of processes the a single process may spread
    * load before updating Connectors.
    */
   int d_max_spread_procs;

   /*!
    * @brief Whether to limit what a process can give to its surplus.
    */
   bool d_limit_supply_to_surplus;

   /*!
    * @brief Whether to reset load obligations within groups that
    * cannot change its load average.
    *
    * This option helps reduce imbalances but makes imbalances
    * caused by bugs more difficult to find.
    *
    * See input parameter "DEV_reset_obligations".
    */
   bool d_reset_obligations;

   /*!
    * @brief Fraction of ideal load a process can accept over and above
    * the ideal.
    *
    * See input parameter "flexible_load_tolerance".
    */
   double d_flexible_load_tol;

   std::vector<double> d_artificial_minimum;

   /*!
    * @brief Boolean to determine whether to use vouchers for transferring load.
    */
   bool d_use_vouchers;

   /*!
    * @brief Metadata operations with timers set according to this object.
    */
   hier::MappingConnectorAlgorithm d_mca;

   /*!
    * @brief Level holding workload data
    */
   mutable std::shared_ptr<hier::PatchLevel> d_workload_level;

   //@{
   //! @name Shared temporaries, used only when actively partitioning.
   mutable hier::BoxLevel* d_balance_box_level;
   mutable hier::Connector* d_balance_to_reference;
   mutable std::shared_ptr<PartitioningParams> d_pparams;
   mutable LoadType d_global_work_sum;
   mutable LoadType d_global_work_avg;
   mutable LoadType d_local_work_max;
   mutable size_t d_num_initial_owners;

   //! @brief Local load subject to change.
   mutable TransitLoad* d_local_load;
   //! @brief Load shipment for sending and receiving.
   mutable TransitLoad* d_shipment;
   //! @brief High-level communication stage.
   mutable tbox::AsyncCommStage d_comm_stage;
   //! @brief High-level peer-to-peer communication object (2 receives, 2 sends).
   mutable tbox::AsyncCommPeer<char> d_comm_peer[4];
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
   std::shared_ptr<tbox::Timer> t_assign_to_local_and_populate_maps;
   std::shared_ptr<tbox::Timer> t_use_map;
   std::shared_ptr<tbox::Timer> t_communication_wait;
   std::shared_ptr<tbox::Timer> t_distribute_load;
   std::shared_ptr<tbox::Timer> t_update_connectors;
   std::shared_ptr<tbox::Timer> t_global_work_reduction;
   std::shared_ptr<tbox::Timer> t_combine_children;
   std::shared_ptr<tbox::Timer> t_balance_children;
   std::shared_ptr<tbox::Timer> t_supply_work;
   std::shared_ptr<tbox::Timer> t_send_shipment;
   std::shared_ptr<tbox::Timer> t_receive_and_unpack_supplied_load;

   //@}

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_print_child_steps;
   char d_check_connectivity;
   char d_check_map;

   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

};

}
}

#endif
