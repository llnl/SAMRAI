/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer_C
#define included_mesh_TreeLoadBalancer_C

#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/mesh/VoucherTransitLoad.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/CenteredRankTree.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cmath>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG0;
const int TreeLoadBalancer::TreeLoadBalancer_LOADTAG1;
const int TreeLoadBalancer::TreeLoadBalancer_FIRSTDATALEN;

const int TreeLoadBalancer::s_default_data_id = -1;

/*
 *************************************************************************
 * TreeLoadBalancer constructor.
 *************************************************************************
 */

TreeLoadBalancer::TreeLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   const std::shared_ptr<tbox::Database>& input_db,
   const std::shared_ptr<tbox::RankTreeStrategy>& rank_tree):
   d_dim(dim),
   d_object_name(name),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mpi_is_dupe(false),
   d_tile_size(dim, 1),
   d_max_spread_procs(500),
   d_voucher_mode(false),
   d_allow_box_breaking(true),
   d_rank_tree(rank_tree ? rank_tree : std::shared_ptr<tbox::RankTreeStrategy>(new tbox::
                                                                                 CenteredRankTree)),
   d_comm_graph_writer(),
   d_master_workload_data_id(s_default_data_id),
   d_flexible_load_tol(0.05),
   d_artificial_minimum(1,1.0),
   d_mca(),
   // Performance evaluation.
   d_barrier_before(false),
   d_barrier_after(false),
   d_report_load_balance(false),
   d_summarize_map(false),
   d_print_steps(false),
   d_check_connectivity(false),
   d_check_map(false)
{
   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);
   setTimers();
   d_mca.setTimerPrefix(d_object_name);
}

/*
 *************************************************************************
 * TreeLoadBalancer constructor.
 *************************************************************************
 */

TreeLoadBalancer::~TreeLoadBalancer()
{
   freeMPICommunicator();
}

/*
 *************************************************************************
 * Accessory functions to get/set load balancing parameters.
 *************************************************************************
 */

bool
TreeLoadBalancer::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}

/*
 **************************************************************************
 **************************************************************************
 */
void
TreeLoadBalancer::setWorkloadPatchDataIndex(
   int data_id,
   int level_number)
{
   std::shared_ptr<pdat::CellDataFactory<double> > datafact(
      SAMRAI_SHARED_PTR_CAST<pdat::CellDataFactory<double>, hier::PatchDataFactory>(
         hier::VariableDatabase::getDatabase()->getPatchDescriptor()->
         getPatchDataFactory(data_id)));

   TBOX_ASSERT(datafact);

   if (level_number >= 0) {
      int asize = static_cast<int>(d_workload_data_id.size());
      if (asize < level_number + 1) {
         d_workload_data_id.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; ++i) {
            d_workload_data_id[i] = d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = data_id;
      }
   } else {
      d_master_workload_data_id = data_id;
      for (int ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ++ln) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}

/*
 *************************************************************************
 * This method implements the abstract LoadBalanceStrategy interface.
 *
 * This method does some preliminary setup then calls
 * loadBalanceWithinRankGroup to do the work.  The set-up includes
 * determining how many cycles to use to gradually balance a severely
 * unbalanced load.
 *
 * After load balancing, it enforces the maximum size restriction
 * by breaking up large boxes and update balance<==>reference again.
 *************************************************************************
 */
void
TreeLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_reference,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
   NULL_USE(domain_box_level);
#endif
   TBOX_ASSERT(!balance_to_reference || balance_to_reference->hasTranspose());
   TBOX_ASSERT(!balance_to_reference ||
      balance_to_reference->isTransposeOf(balance_to_reference->getTranspose()));
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY6(d_dim,
      balance_box_level,
      min_size,
      max_size,
      domain_box_level,
      bad_interval,
      cut_factor);
   
   if (hierarchy) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *hierarchy);
   }

   size_t minimum_cells = 1;
   double artificial_minimum = 1.0;
   if (hierarchy) {
      minimum_cells = hierarchy->getMinimumCellRequest(level_number);
      if (level_number < static_cast<int>(d_artificial_minimum.size())) {
         artificial_minimum = d_artificial_minimum[level_number];
      } else {
         artificial_minimum = d_artificial_minimum.back();
      }
      TBOX_ASSERT(artificial_minimum >= 0.0);
   }

   if (d_mpi_is_dupe) {
      /*
       * If user has set the duplicate communicator, make sure it is
       * compatible with the BoxLevel involved.
       */
      TBOX_ASSERT(d_mpi.getSize() == balance_box_level.getMPI().getSize());
      TBOX_ASSERT(d_mpi.getRank() == balance_box_level.getMPI().getRank());
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!d_mpi.isCongruentWith(balance_box_level.getMPI())) {
         TBOX_ERROR("TreeLoadBalancer::loadBalanceBoxLevel:\n"
            << "The input balance_box_level has a SAMRAI_MPI that is\n"
            << "not congruent with the one set with setSAMRAI_MPI().\n"
            << "You must use freeMPICommunicator() before balancing\n"
            << "a BoxLevel with an incongruent SAMRAI_MPI.");
      }
#endif
   } else {
      d_mpi = balance_box_level.getMPI();
   }

   if (d_print_steps) {
      tbox::plog << d_object_name << "::loadBalanceBoxLevel called with:"
                 << "\n  min_size = " << min_size
                 << "\n  max_size = " << max_size
                 << "\n  bad_interval = " << bad_interval
                 << "\n  cut_factor = " << cut_factor
                 << "\n  prebalance:\n"
                 << balance_box_level.format("  ", 2);
   }

   // Set effective_cut_factor to least common multiple of cut_factor and d_tile_size.
   const size_t nblocks = balance_box_level.getGridGeometry()->getNumberBlocks();
   hier::IntVector effective_cut_factor(cut_factor, nblocks);
   if ( d_tile_size != hier::IntVector::getOne(d_dim) ) {
      if (cut_factor.getNumBlocks() == 1) {
         for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
            for ( int d=0; d<d_dim.getValue(); ++d ) {
               while ( effective_cut_factor(b,d)/d_tile_size[d]*d_tile_size[d] != effective_cut_factor(b,d) ) {
                  effective_cut_factor(b,d) += cut_factor[d];
               }
            }
         }
      }  else {
         for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
            for ( int d=0; d<d_dim.getValue(); ++d ) {
               while ( effective_cut_factor(b,d)/d_tile_size[d]*d_tile_size[d] != effective_cut_factor(b,d) ) {
                  effective_cut_factor(b,d) += cut_factor(b,d);
               }
            }
         }
      }
      if (d_print_steps) {
         tbox::plog << d_object_name << "::loadBalanceBoxLevel effective_cut_factor = "
                    << effective_cut_factor << std::endl;
      }
   }

   /*
    * Periodic image Box should be ignored during load balancing
    * because they have no real work.  The load-balanced results
    * should contain no periodic images.
    *
    * To avoid need for special logic to skip periodic images while
    * load balancing, we just remove periodic images in the
    * balance_box_level and all periodic edges in
    * reference<==>balance.
    */

   balance_box_level.removePeriodicImageBoxes();
   if (balance_to_reference) {

      balance_to_reference->getTranspose().removePeriodicRelationships();
      balance_to_reference->getTranspose().setHead(balance_box_level, true);

      balance_to_reference->removePeriodicRelationships();
      balance_to_reference->setBase(balance_box_level, true);

   }

   if (d_barrier_before) {
      t_barrier_before->start();
      d_mpi.Barrier();
      t_barrier_before->stop();
   }

   if (!rank_group.containsAllRanks()) {
      BalanceUtilities::prebalanceBoxLevel(
         balance_box_level,
         balance_to_reference,
         rank_group);
   }

   t_load_balance_box_level->start();

   d_pparams = std::make_shared<PartitioningParams>(
         *balance_box_level.getGridGeometry(),
         balance_box_level.getRefinementRatio(),
         min_size, max_size, bad_interval, effective_cut_factor, minimum_cells,
         artificial_minimum,
         d_flexible_load_tol);

   /*
    * We expect the domain box_level to be in globalized state.
    */
   TBOX_ASSERT(
      domain_box_level.getParallelState() ==
      hier::BoxLevel::GLOBALIZED);

   LoadType local_load = computeLocalLoad(balance_box_level);

   LoadType max_local_load = local_load;

   LoadType global_sum_load = local_load;

   size_t nproc_with_initial_load =
      balance_box_level.getLocalNumberOfBoxes() > 0;

   /*
    * Determine the total load and number of processes that has any
    * initial load.
    */
   t_compute_global_load->start();
   if (d_mpi.getSize() > 1) {
      double dtmp[2], dtmp_sum[2], dtmp_max[2];

      dtmp[0] = local_load;
      dtmp[1] = static_cast<double>(nproc_with_initial_load);

      d_mpi.Allreduce(dtmp, dtmp_sum, 2, MPI_DOUBLE, MPI_SUM);
      global_sum_load = dtmp_sum[0];
      nproc_with_initial_load = (size_t)dtmp_sum[1];

      d_mpi.Allreduce(dtmp, dtmp_max, 1, MPI_DOUBLE, MPI_MAX);
      max_local_load = dtmp_max[0];

   }
   t_compute_global_load->stop();

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << d_object_name << "::loadBalanceBoxLevel"
                 << " max_local_load=" << max_local_load
                 << " global_sum_load=" << global_sum_load
                 << " (initially born on "
                 << nproc_with_initial_load << " procs) across all "
                 << d_mpi.getSize()
                 << " procs, averaging " << global_sum_load / d_mpi.getSize()
                 << " or " << pow(global_sum_load / d_mpi.getSize(), 1.0 / d_dim.getValue())
                 << "^" << d_dim << " per proc." << std::endl;
   }

   d_global_avg_load = global_sum_load / rank_group.size();

   /*
    * Compute how many balancing cycles to use based on severity of
    * imbalance, using formula
    * d_max_spread_procs^number_of_cycles >= fanout_size.
    *
    * The objective of balancing over multiple cycles is to avoid
    * unscalable performance in cases where just a few processes own
    * most of the initial load.  Each cycle spreads the load out a
    * little more.  By slowly spreading out the load, no process has
    * to set up unbalanced<==>balanced with number of relationships
    * that scales with the machine size.
    *
    * Exception: If given a RankGroup with less than all ranks, we
    * treat it as a specific user request to balance only within the
    * RankGroup and just use the RankGroup as is.  We are not set up
    * to support such a request and multi-cycling simultaneously.
    */
   const double fanout_size = d_global_avg_load > d_pparams->getLoadComparisonTol() ?
      max_local_load / d_global_avg_load : 1.0;
   const int number_of_cycles = !rank_group.containsAllRanks() ? 1 :
      int(ceil(log(fanout_size) / log(static_cast<double>(d_max_spread_procs))));
   if (d_print_steps) {
      tbox::plog << d_object_name << "::loadBalanceBoxLevel"
                 << " max_spread_procs=" << d_max_spread_procs
                 << " fanout_size=" << fanout_size
                 << " number_of_cycles=" << number_of_cycles
                 << std::endl;
   }

   /*
    * The icycle loop spreads out the work each time through.  If
    * using more than one cycle, only the last one tries to balance
    * across all processes.
    */

   for (int icycle = number_of_cycles - 1; icycle >= 0; --icycle) {

      // If not the first cycle, local_load needs updating.
      if (icycle != number_of_cycles - 1) {
         local_load = computeLocalLoad(balance_box_level);
      }

      if (d_report_load_balance) {
         // Debugging: check overall load balance at intermediate cycles.
         tbox::plog
         << d_object_name << "::loadBalanceBoxLevel results before cycle "
         << icycle << ":" << std::endl;
         BalanceUtilities::reduceAndReportLoadBalance(
            std::vector<double>(1, local_load),
            balance_box_level.getMPI());
      }

      const bool last_cycle = (icycle == 0);

      /*
       * Determine whether to use rank_group as is or subgroup it for
       * intemediate cycles.
       *
       * Set cycle_rank_group to either the input rank_group or a
       * subgroup that is a function of the cycle number.
       */
      int number_of_groups = 1;
      int group_num = 0;

      tbox::RankGroup cycle_rank_group(d_mpi);
      if (!last_cycle) {
         createBalanceRankGroupBasedOnCycles(
            cycle_rank_group,
            number_of_groups,
            group_num,
            1 - double(icycle) / number_of_cycles);
      }

      /*
       * Compute the group's load.
       */
      t_compute_tree_load->start();

      double group_sum_load;

      if (last_cycle) {

         group_sum_load = global_sum_load;

      } else {

         t_compute_tree_load_for_cycle[icycle]->start();

         /*
          * Use MPI's vector all-reduce to get individual group loads.
          * This gives more info than the process needs, but because the
          * number of groups << number of procs, it is still faster
          * (probably) than hand coded communication.
          */
         std::vector<double> group_loads(number_of_groups, 0.0);
         group_loads[group_num] = local_load;
         if (d_mpi.getSize() > 1) {
            d_mpi.AllReduce(&group_loads[0],
               static_cast<int>(group_loads.size()),
               MPI_SUM);
         }
         group_sum_load = group_loads[group_num];

         t_compute_tree_load_for_cycle[icycle]->stop();

      }

      t_compute_tree_load->stop();

      if (d_print_steps) {
         tbox::plog << d_object_name << "::loadBalanceBoxLevel"
                    << " cycle number=" << icycle
                    << " number_of_groups=" << number_of_groups
                    << " my group_num=" << group_num
                    << " my group size=" << cycle_rank_group.size()
                    << " my group_sum_load=" << group_sum_load
                    << std::endl;
      }

      // Run the tree load balancing algorithm.
      t_load_balance_for_cycle[icycle]->start();
      loadBalanceWithinRankGroup(
         balance_box_level,
         balance_to_reference,
         cycle_rank_group,
         group_sum_load);
      t_load_balance_for_cycle[icycle]->stop();

      if (d_barrier_after) {
         t_barrier_after->start();
         d_mpi.Barrier();
         t_barrier_after->stop();
      }

   }

   /*
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */

   hier::IntVector max_intvector(d_dim, tbox::MathUtilities<int>::getMax());
   if (max_size != max_intvector) {

      t_constrain_size->barrierAndStart();
      BalanceUtilities::constrainMaxBoxSizes(
         balance_box_level,
         balance_to_reference ? &balance_to_reference->getTranspose() : 0,
         *d_pparams);
      t_constrain_size->stop();

      if (d_print_steps) {
         tbox::plog << " TreeLoadBalancer completed constraining box sizes."
                    << "\n";
      }

   }

   /*
    * Finished load balancing.  Clean up and wrap up.
    */

   d_pparams.reset();

   t_load_balance_box_level->stop();

   local_load = computeLocalLoad(balance_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_box_level.getBoxes().size()));

   if (d_print_steps) {
      tbox::plog << "Post balanced:\n" << balance_box_level.format("", 2);
   }

   if (d_report_load_balance) {
      t_report_loads->start();
      tbox::plog
      << d_object_name << "::loadBalanceBoxLevel results after "
      << number_of_cycles << " cycles:" << std::endl;
      BalanceUtilities::reduceAndReportLoadBalance(
         std::vector<double>(1, local_load),
         balance_box_level.getMPI());
      t_report_loads->stop();
   }

   if (d_check_connectivity && balance_to_reference) {
      hier::Connector& reference_to_balance = balance_to_reference->getTranspose();
      tbox::plog << "TreeLoadBalancer checking balance-reference connectivity."
                 << std::endl;
      int errs = 0;
      if (reference_to_balance.checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in reference_to_balance!\n";
      }
      if (balance_to_reference->checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_reference!\n";
      }
      if (reference_to_balance.checkTransposeCorrectness(*balance_to_reference)) {
         ++errs;
         tbox::perr << "Error found in balance-reference transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found.\n"
            << "reference_box_level:\n" << reference_to_balance.getBase().format("", 2)
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "reference_to_balance:\n" << reference_to_balance.format("", 2)
            << "balance_to_reference:\n" << balance_to_reference->format("", 2));
      }
      tbox::plog << "TreeLoadBalancer checked balance-reference connectivity."
                 << std::endl;
   }

   if (d_barrier_after) {
      t_barrier_after->start();
      d_mpi.Barrier();
      t_barrier_after->stop();
   }

}

/*
 *************************************************************************
 * Given an "unbalanced" BoxLevel, load balance it within the given
 * RankGroup using the tree load balancing algorithm and update
 * Connectors.
 *************************************************************************
 */
void
TreeLoadBalancer::loadBalanceWithinRankGroup(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_reference,
   const tbox::RankGroup& rank_group,
   const double group_sum_load) const
{
   TBOX_ASSERT(!balance_to_reference || balance_to_reference->hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, balance_box_level);

   /*
    * Initialize empty balanced_box_level and mappings so they are
    * ready to be populated.
    */
   hier::BoxLevel balanced_box_level(
      balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI());
   hier::MappingConnector balanced_to_unbalanced(balanced_box_level,
         balance_box_level,
         hier::IntVector::getZero(d_dim));
   hier::MappingConnector unbalanced_to_balanced(balance_box_level,
         balanced_box_level,
         hier::IntVector::getZero(d_dim));
   unbalanced_to_balanced.setTranspose(&balanced_to_unbalanced, false);

   t_get_map->start();

   if (!rank_group.isMember(d_mpi.getRank())) {
      /*
       * If the local process is not a member of the RankGroup, it
       * does not participate in the work and just sets the output
       * objects to be locally empty.
       *
       * The following assert should be guaranteed by an earlier call
       * to prebalanceBoxLevel.  Having boxes without being in the
       * given rank group leads to undefined results.
       */
      TBOX_ASSERT(balance_box_level.getLocalNumberOfBoxes() == 0);

      t_post_load_distribution_barrier->start();
      d_mpi.Barrier();
      t_post_load_distribution_barrier->stop();

   } else {

      /*
       * Create a concrete TransitLoad container to hold the work
       * being redistributed.
       */
      std::shared_ptr<TransitLoad> balanced_work;
      if (d_voucher_mode) {
         balanced_work = std::make_shared<VoucherTransitLoad>(*d_pparams);
         balanced_work->setTimerPrefix(d_object_name + "::VoucherTransitLoad");
      } else {
         balanced_work = std::make_shared<BoxTransitSet>(*d_pparams);
         balanced_work->setTimerPrefix(d_object_name + "::BoxTransitSet");
      }

      distributeLoadAcrossRankGroup(
         *balanced_work,
         balance_box_level,
         rank_group,
         group_sum_load);

      t_post_load_distribution_barrier->start();
      d_mpi.Barrier();
      t_post_load_distribution_barrier->stop();

      if (d_print_steps) {
         tbox::plog << d_object_name
                    << "::loadBalanceWithinRankGroup constructing unbalanced<==>balanced.\n";
      }
      t_assign_to_local_and_populate_maps->start();
      balanced_work->assignToLocalAndPopulateMaps(
         balanced_box_level,
         balanced_to_unbalanced,
         unbalanced_to_balanced,
         d_flexible_load_tol,
         d_mpi);
      t_assign_to_local_and_populate_maps->stop();
      if (d_print_steps) {
         tbox::plog << d_object_name
                    <<
         "::loadBalanceWithinRankGroup finished constructing unbalanced<==>balanced.\n";
      }

   }

   t_get_map->stop();

   if (d_summarize_map) {
      tbox::plog << d_object_name << "::loadBalanceWithinRankGroup unbalanced--->balanced map:\n"
                 << unbalanced_to_balanced.format("\t", 0)
                 << "Map statistics:\n" << unbalanced_to_balanced.formatStatistics("\t")
                 << d_object_name << "::loadBalanceWithinRankGroup balanced--->unbalanced map:\n"
                 << balanced_to_unbalanced.format("\t", 0)
                 << "Map statistics:\n" << balanced_to_unbalanced.formatStatistics("\t")
                 << '\n';
   }

   if (d_check_map) {
      if (unbalanced_to_balanced.findMappingErrors() != 0) {
         TBOX_ERROR(
            d_object_name
            << "::loadBalanceWithinRankGroup Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            d_object_name << "::loadBalanceWithinRankGroup Transpose errors found!");
      }
   }

   if (d_summarize_map) {
      tbox::plog << d_object_name << "::loadBalanceWithinRankGroup: unbalanced--->balanced map:\n"
                 << unbalanced_to_balanced.format("\t", 0)
                 << "Map statistics:\n" << unbalanced_to_balanced.formatStatistics("\t")
                 << d_object_name << "::loadBalanceWithinRankGroup: balanced--->unbalanced map:\n"
                 << balanced_to_unbalanced.format("\t", 0)
                 << "Map statistics:\n" << balanced_to_unbalanced.formatStatistics("\t")
                 << '\n';
   }

   if (balance_to_reference && balance_to_reference->hasTranspose()) {
      t_use_map->barrierAndStart();
      d_mca.modify(
         balance_to_reference->getTranspose(),
         unbalanced_to_balanced,
         &balance_box_level,
         &balanced_box_level);
      t_use_map->barrierAndStop();
   } else {
      hier::BoxLevel::swap(balance_box_level, balanced_box_level);
   }

   if (d_print_steps) {
      tbox::plog << d_object_name << "::LoadBalanceWithinRankGroup: returning"
                 << std::endl;
   }
}

/*
 *************************************************************************
 * Distribute load on the tree and generate unbalanced<==>balanced
 * maps.
 *************************************************************************
 */
void
TreeLoadBalancer::distributeLoadAcrossRankGroup(
   TransitLoad& balanced_work,
   const hier::BoxLevel& unbalanced_box_level,
   const tbox::RankGroup& rank_group,
   double group_sum_load) const
{
   t_distribute_load_across_rank_group->start();

   TBOX_ASSERT(balanced_work.getNumberOfItems() == 0);
   TBOX_ASSERT(balanced_work.getSumLoad() <= d_pparams->getLoadComparisonTol());

   double group_avg_load = group_sum_load / rank_group.size();

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << d_object_name << "::LoadBalanceWithinRankGroup balancing "
                 << group_sum_load << " units in group of "
                 << rank_group.size() << " procs, averaging " << group_avg_load
                 << " or " << pow(group_avg_load, 1.0 / d_dim.getValue())
                 << "^" << d_dim << " per proc."
                 << "  Avg is " << group_avg_load
      / static_cast<double>(d_pparams->getMinBoxSize().getProduct())
                 << " times min size of " << d_pparams->getMinBoxSize()
                 << std::endl;
   }

   /*
    * Before the last cycle, it is possible for the group average load
    * to be below the global average, if the group just happens to have
    * underloaded processors.  However, there is no point in driving the
    * processor loads down below the global average just to have it
    * brought back up by the last cycle.  It just unnecessarily fragments
    * the boxes and costs more to do.  To prevent this, reset the group
    * average to the global average if it is below.
    */
   group_avg_load =
      tbox::MathUtilities<double>::Max(group_avg_load, d_global_avg_load);

   // Set parameters governing box breaking.
   balanced_work.setAllowBoxBreaking(d_allow_box_breaking);
   const double ideal_box_width = pow(group_avg_load, 1.0 / d_dim.getValue());
   balanced_work.setThresholdWidth(1.0 * ideal_box_width);
   if (d_print_steps) {
      tbox::plog << d_object_name << "::distributeLoadAcrossRankGroup: ideal_box_width = "
                 << ideal_box_width
                 << "\n  Set threshold width to " << balanced_work.getThresholdWidth()
                 << std::endl;
   }

   /*
    * Arrange the group ranks in a rank tree in order to get
    * the parent/children in the group.
    */
   d_rank_tree->setupTree(rank_group, d_mpi.getRank());
   if (d_print_steps) {
      // Write local part of tree to log.
      tbox::plog << "TreeLoadBalancer tree:\n"
                 << "  Root rank: " << d_rank_tree->getRootRank() << '\n'
                 << "  Child number: " << d_rank_tree->getChildNumber() << '\n'
                 << "  Generation number: " << d_rank_tree->getGenerationNumber() << '\n'
                 << "  Number of children: " << d_rank_tree->getNumberOfChildren() << '\n'
                 << "  Local relatives: "
                 << "  " << d_rank_tree->getParentRank() << " <- [" << d_rank_tree->getRank()
                 << "] -> {";
      for (unsigned int i = 0; i < d_rank_tree->getNumberOfChildren(); ++i) {
         tbox::plog << ' ' << d_rank_tree->getChildRank(i);
      }
      tbox::plog << " }" << std::endl;
   }

   const int num_children = d_rank_tree->getNumberOfChildren();

   /*
    * Communication objects for sending to/receiving from
    * parent/children: We could combine all of these AsyncCommStages
    * and most of the AsyncCommPeers, but we intentionally keep them
    * separate to aid performance analysis.
    */

   tbox::AsyncCommStage child_send_stage;
   tbox::AsyncCommPeer<char>* child_sends = 0;
   tbox::AsyncCommStage parent_send_stage;
   tbox::AsyncCommPeer<char>* parent_send = 0;

   setupAsyncCommObjects(
      child_send_stage,
      child_sends,
      parent_send_stage,
      parent_send,
      rank_group);
   child_send_stage.setCommunicationWaitTimer(t_child_send_wait);
   parent_send_stage.setCommunicationWaitTimer(t_parent_send_wait);

   tbox::AsyncCommStage child_recv_stage;
   tbox::AsyncCommPeer<char>* child_recvs = 0;
   tbox::AsyncCommStage parent_recv_stage;
   tbox::AsyncCommPeer<char>* parent_recv = 0;

   setupAsyncCommObjects(
      child_recv_stage,
      child_recvs,
      parent_recv_stage,
      parent_recv,
      rank_group);
   child_recv_stage.setCommunicationWaitTimer(t_child_recv_wait);
   parent_recv_stage.setCommunicationWaitTimer(t_parent_recv_wait);

   /*
    * Essential outline of the tree load balancing algorithm as implemented:
    *
    * 1. For each child of the local process:
    * Receive data from branch rooted at child (nodes in
    * branch, excess work, remaining work in branch, etc.).
    *
    * 2. Compute data for branch rooted at self by combining
    * local data with children branch data.
    *
    * 3. If parent exists:
    * Send branch info to parent, including work, if any.
    *
    * 4. If parent exists and we need more work:
    * Receive additional work from parent.
    *
    * 5. For each child who asked for more work:
    * Partition requested work amount and send.
    */

   /*
    * Step 1:
    *
    * Post receive for data from branch rooted at children.
    * We have to do some local setup, but post the receive
    * now to overlap communication.
    */
   t_get_load_from_children->start();
   for (int c = 0; c < num_children; ++c) {
      child_recvs[c].setRecvTimer(t_child_recv_wait);
      child_recvs[c].setWaitTimer(t_child_recv_wait);
      child_recvs[c].beginRecv();
      if (child_recvs[c].isDone()) {
         child_recvs[c].pushToCompletionQueue();
      }
   }
   t_get_load_from_children->stop();

   // State of the tree, as seen by local process.
   BranchData my_branch(*d_pparams, balanced_work);
   my_branch.setTimerPrefix(d_object_name);
   my_branch.setPrintSteps(d_print_steps);
   std::vector<BranchData> child_branches(num_children, my_branch);

   /*
    * Step 2, local part:
    *
    * unassigned is a container of loads that have been released by
    * their owners and have not yet been assigned to another.  First,
    * put all initial local work in unassigned.  Received loads are
    * placed here before determining whether to keep them or send them
    * to another part of the tree.
    *
    * unassigned is a reference to balanced_work, because at the end
    * of the algorithm, everything left unassigned is actually the
    * balanced load.
    */
   TransitLoad& unassigned(balanced_work);

   t_local_load_moves->start();
   unassigned.insertAll(unbalanced_box_level.getBoxes());
   t_local_load_moves->stop();

   my_branch.setStartingLoad(group_avg_load,
      unassigned.getSumLoad(),
      group_avg_load * (1 + d_flexible_load_tol));

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << "Initial local load: ";
      unassigned.recursivePrint(tbox::plog, "", 1);
      my_branch.recursivePrint(tbox::plog, "", 0);
      tbox::plog << std::endl;
   }

   /*
    * Step 2, remote part:
    *
    * Finish getting tree and load data from children, incorporate
    * its data and add any imported work to unassigned bin.
    */
   t_get_load_from_children->start();
   while (child_recv_stage.hasCompletedMembers() ||
          child_recv_stage.advanceSome()) {

      tbox::AsyncCommPeer<char>* child_recv =
         CPP_CAST<tbox::AsyncCommPeer<char> *>(child_recv_stage.popCompletionQueue());
      TBOX_ASSERT(child_recv != 0);
      TBOX_ASSERT(child_recv >= child_recvs);
      TBOX_ASSERT(child_recv < child_recvs + num_children);

      const int cindex = static_cast<int>(child_recv - child_recvs);
      TBOX_ASSERT(cindex >= 0 && cindex < num_children);

      // Extract data from the child and store in child_branches[cindex].
      tbox::MessageStream mstream(child_recv->getRecvSize(),
                                  tbox::MessageStream::Read,
                                  child_recv->getRecvData(),
                                  false);

      if (d_print_steps) {
         tbox::plog << "Unpacking from child "
                    << cindex << ':' << d_rank_tree->getChildRank(cindex) << ":\n";
      }
      child_branches[cindex].unpackDataFromChild(mstream);

      t_local_load_moves->start();
      my_branch.incorporateChild(child_branches[cindex]);
      child_branches[cindex].moveInboundLoadToReserve(unassigned);
      t_local_load_moves->stop();

   }

   // We should have received everything by this point.
   TBOX_ASSERT(!child_recv_stage.hasPendingRequests());

   size_t unassigned_highwater = unassigned.getNumberOfItems();

   /*
    * TODO: Maybe this should be deficit() or
    * max(deficit(),effDeficit()) instead of effDeficit().  The
    * argument for using deficit() is: if a branch ends up underloaded
    * (for any reason) its deficit should be kept within its parent's
    * branch.  The argument for using effDeficit() is: if part of the
    * branch took more work, the rest of the branch should take the
    * fair share because it can be done without aggravating overloads
    * and it would help keeping the parent branch's work within the
    * parent's branch (preserve locality).  Not using deficit() to
    * decide to get work from parent may be responsible for the 20%
    * unbalance observed at 1M procs.  There should be a corresponding
    * change to computeSurplusPerEffectiveDescendent().
    */
   if (my_branch.effDeficit() > 0 && !d_rank_tree->isRoot()) {
      my_branch.setWantsWorkFromParent();
   }

   if (d_print_steps) {
      tbox::plog << "Received children branch data." << std::endl;
      for (int c = 0; c < num_children; ++c) {
         tbox::plog << "Child "
                    << c << ':' << d_rank_tree->getChildRank(c)
                    << " branch:\n";
         child_branches[c].recursivePrint(tbox::plog, "  ");
      }
      tbox::plog << "Initial branch:\n";
      my_branch.recursivePrint(tbox::plog, "  ");
      tbox::plog << "unassigned: ";
      unassigned.recursivePrint(tbox::plog, "  ", 1);
   }

   t_get_load_from_children->stop();

   /*
    * Step 3:
    *
    * Send to parent (if any) branch info and excess work.
    */
   t_send_load_to_parent->start();
   if (parent_send != 0) {

      if (my_branch.effExcess() > 0) {

         /*
          * Try to send work in the range of [effective excess,surplus].
          * Sending less than effective excess would overload this branch.
          * Sending more than the surplus would overload its complement.
          *
          * It is possible to have effective excess > surplus.  In
          * these cases, don't send more than the surplus, because
          * that would overload the complement of this branch.  It is
          * better to overload this branch than to progressively push
          * unwanted load up the tree, making the root extremely
          * overloaded.  Keeping the overload in the branch also
          * maintains some data locality.
          */
         const LoadType export_load_low = tbox::MathUtilities<double>::Min(
               my_branch.effExcess(), my_branch.surplus());
         const LoadType export_load_high = my_branch.surplus();
         const LoadType export_load_ideal = my_branch.surplus();

         if (d_print_steps) {
            tbox::plog << "Pushing to parent rank "
                       << d_rank_tree->getParentRank() << ' ' << export_load_ideal
                       << " [" << export_load_low << ", " << export_load_high << ']'
                       << std::endl;
         }

         t_local_load_moves->start();
         my_branch.adjustOutboundLoad(
            unassigned,
            export_load_ideal,
            export_load_low,
            export_load_high);
         t_local_load_moves->stop();

      }

      tbox::MessageStream mstream;
      my_branch.packDataToParent(mstream);
      parent_send->setSendTimer(t_parent_send_wait);
      parent_send->setWaitTimer(t_parent_send_wait);
      parent_send->beginSend(static_cast<const char *>(mstream.getBufferStart()),
         mstream.getCurrentSize());

   }
   t_send_load_to_parent->stop();

   /*
    * Step 4:
    *
    * Get work from parent, if any.
    */
   if (my_branch.getWantsWorkFromParent()) {
      t_get_load_from_parent->start();

      parent_recv->setRecvTimer(t_parent_recv_wait);
      parent_recv->setWaitTimer(t_parent_recv_wait);

      parent_recv->beginRecv();
      parent_recv->completeCurrentOperation();

      if (d_print_steps) {
         tbox::plog << "Received from parent " << d_rank_tree->getParentRank()
                    << "... Unpacking." << std::endl;
      }

      tbox::MessageStream mstream(parent_recv->getRecvSize(),
                                  tbox::MessageStream::Read,
                                  parent_recv->getRecvData(),
                                  false);
      my_branch.unpackDataFromParentAndIncorporate(mstream);

      t_local_load_moves->start();
      my_branch.moveInboundLoadToReserve(unassigned);
      t_local_load_moves->stop();

      if (unassigned_highwater < unassigned.getNumberOfItems()) {
         unassigned_highwater = unassigned.getNumberOfItems();
      }

      t_get_load_from_parent->stop();
   } else {
      if (d_print_steps) {
         tbox::plog << "Did not request work from parent.\n";
      }
   }

   if (d_print_steps) {
      tbox::plog << "Postparent branch:\n";
      my_branch.recursivePrint(tbox::plog, "  ");
      tbox::plog << "unassigned: ";
      unassigned.recursivePrint(tbox::plog, "  ", 1);
   }

   /*
    * Step 5:
    *
    * Reassign and send work to each child that requested work.
    */

   t_send_load_to_children->start();

   for (int ichild = 0; ichild < num_children; ++ichild) {

      BranchData& recip_branch = child_branches[ichild];

      if (recip_branch.getWantsWorkFromParent()) {

         const LoadType surplus_per_eff_des =
            computeSurplusPerEffectiveDescendent(
               unassigned.getSumLoad(),
               group_avg_load,
               child_branches,
               ichild);

         const LoadType export_load_ideal = recip_branch.effDeficit()
            + (surplus_per_eff_des < 0.0 ? 0.0 :
               surplus_per_eff_des * recip_branch.numProcsEffective());

         const LoadType export_load_low = recip_branch.effDeficit()
            + surplus_per_eff_des * recip_branch.numProcsEffective();

         const LoadType export_load_high =
            tbox::MathUtilities<double>::Max(export_load_ideal,
               recip_branch.effMargin());

         if (d_print_steps) {
            tbox::plog << "Pushing to child " << ichild << " ("
                       << d_rank_tree->getChildRank(ichild) << ") " << export_load_ideal
                       << " [" << export_load_low << ", " << export_load_high << ']'
                       << std::endl;
         }

         t_local_load_moves->start();
         recip_branch.adjustOutboundLoad(
            unassigned,
            export_load_ideal,
            export_load_low,
            export_load_high);
         t_local_load_moves->stop();

         if (d_print_steps) {
            tbox::plog << "Packing data to child "
                       << ichild << ':' << d_rank_tree->getChildRank(ichild)
                       << std::endl;
         }
         tbox::MessageStream mstream;
         recip_branch.packDataToChild(mstream);
         child_sends[ichild].setSendTimer(t_child_send_wait);
         child_sends[ichild].setWaitTimer(t_child_send_wait);
         child_sends[ichild].beginSend(static_cast<const char *>(mstream.getBufferStart()),
            mstream.getCurrentSize());

      }

   }

   t_send_load_to_children->stop();

   if (d_print_steps) {
      tbox::plog << "After settling parent and children, unassigned is: ";
      unassigned.recursivePrint(tbox::plog, "  ", 1);
      tbox::plog << "Final local load (normalized) is "
                 << (balanced_work.getSumLoad() / group_avg_load)
                 << ", surplus = " << (balanced_work.getSumLoad() - group_avg_load)
                 << ", excess = "
                 << (balanced_work.getSumLoad() - (1 + d_flexible_load_tol) * group_avg_load)
                 << std::endl;
      tbox::plog << "Final branch:\n";
      my_branch.recursivePrint(tbox::plog, "  ");

      for (int ichild = 0; ichild < num_children; ++ichild) {
         tbox::plog << "Final child "
                    << ichild << ":" << d_rank_tree->getChildRank(ichild)
                    << " branch:\n";
         child_branches[ichild].recursivePrint(tbox::plog, "  ");
      }
   }

   /*
    * Finish messages before starting edge info exchange.
    * We have only sends to complete, so it should not take
    * long to advance them all to completion.
    */
   if (d_print_steps) {
      tbox::plog << d_object_name
                 << "::loadBalanceWithinRankGroup: waiting for sends to complete.\n";
   }

   t_finish_sends->start();
   parent_send_stage.advanceAll();
   parent_send_stage.clearCompletionQueue();
   child_send_stage.advanceAll();
   child_send_stage.clearCompletionQueue();
   t_finish_sends->stop();

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < num_children; ++i) {
      TBOX_ASSERT(child_sends[i].isDone());
      TBOX_ASSERT(child_recvs[i].isDone());
   }
   if (parent_send != 0) {
      TBOX_ASSERT(parent_send->isDone());
      TBOX_ASSERT(parent_recv->isDone());
   }
#endif
   if (d_print_steps) {
      tbox::plog << d_object_name << "::loadBalanceWithinRankGroup: completed sends.\n";
   }

   if (d_comm_graph_writer) {
      /*
       * To evaluate performance of the algorithm, record these edges:
       * - Two edges to parent: load shipped, and items shipped.
       * - Same two edges from parent, plus one more for message size.
       * - Two timer for edges from children, one from parent.
       * Record these nodes:
       * - Number of final boxes:
       */
      d_comm_graph_writer->addRecord(d_mpi, size_t(10), size_t(7));

      const int prank = (d_rank_tree->isRoot() ? -1 : d_rank_tree->getParentRank());

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(0),
         "load up",
         double(my_branch.getWantsWorkFromParent() ? 0 : my_branch.getShipmentLoad()),
         tbox::CommGraphWriter::TO,
         prank);

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(1),
         "items up",
         double(my_branch.getWantsWorkFromParent() ? 0 : my_branch.getShipmentPackageCount()),
         tbox::CommGraphWriter::TO,
         prank);

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(2),
         "origins up",
         double(my_branch.getWantsWorkFromParent() ? 0 : my_branch.getShipmentOriginatorCount()),
         tbox::CommGraphWriter::TO,
         prank);

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(3),
         "load down",
         double(my_branch.getWantsWorkFromParent() ? my_branch.getShipmentLoad() : 0),
         tbox::CommGraphWriter::FROM,
         (my_branch.getWantsWorkFromParent() ? prank : -1));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(4),
         "items down",
         double(my_branch.getWantsWorkFromParent() ? my_branch.getShipmentPackageCount() : 0),
         tbox::CommGraphWriter::FROM,
         (my_branch.getWantsWorkFromParent() ? prank : -1));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(5),
         "origins down",
         double(my_branch.getWantsWorkFromParent() ? my_branch.getShipmentOriginatorCount() : 0),
         tbox::CommGraphWriter::FROM,
         (my_branch.getWantsWorkFromParent() ? prank : -1));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(6),
         "bytes down",
         double(my_branch.getWantsWorkFromParent() ? parent_recv->getRecvSize() : int(0)),
         tbox::CommGraphWriter::FROM,
         (my_branch.getWantsWorkFromParent() ? prank : -1));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(7),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(0));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(8),
         "child wait",
         t_child_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         d_rank_tree->getChildRank(1));

      d_comm_graph_writer->setEdgeInCurrentRecord(
         size_t(9),
         "parent wait",
         t_parent_recv_wait->getTotalWallclockTime(),
         tbox::CommGraphWriter::FROM,
         (my_branch.getWantsWorkFromParent() ? prank : -1));

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(0),
         "initial box count",
         double(unbalanced_box_level.getLocalNumberOfBoxes()));

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(1),
         "initial load",
         double(unbalanced_box_level.getLocalNumberOfCells()) / group_avg_load);

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(2),
         "final box count",
         double(balanced_work.getNumberOfItems()));

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(3),
         "final load",
         double(balanced_work.getSumLoad()) / group_avg_load);

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(4),
         "final surplus",
         double(balanced_work.getSumLoad()) - group_avg_load);

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(5),
         "branch surplus",
         double(my_branch.surplus()));

      d_comm_graph_writer->setNodeValueInCurrentRecord(
         size_t(6),
         "unassigned highwater",
         double(unassigned_highwater));

   }

   destroyAsyncCommObjects(child_sends, parent_send);
   destroyAsyncCommObjects(child_recvs, parent_recv);

   t_distribute_load_across_rank_group->stop();
}

/*
 *************************************************************************
 * Compute the surplus per descendent in the effective tree rooted at
 * the local process.  This surplus, if any, is the difference between
 * the available load_for_descendents, spread over those descendents.
 *
 * Surplus per descendent will be zero if we don't need to push surplus
 * work to children.
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::computeSurplusPerEffectiveDescendent(
   const LoadType& unassigned_load,
   const LoadType& group_avg_load,
   const std::vector<BranchData>& child_branches,
   int first_child) const
{
   const LoadType load_for_me = group_avg_load * (1 + d_flexible_load_tol);

   const LoadType load_for_descendents = unassigned_load - load_for_me;

   LoadType ideal_export_to_children = 0.0;
   int num_effective_des = 0;
   for (size_t ichild = first_child; ichild < child_branches.size(); ++ichild) {
      if (child_branches[ichild].getWantsWorkFromParent()) {
         ideal_export_to_children += child_branches[ichild].effDeficit();
         num_effective_des += child_branches[ichild].numProcsEffective();
      }
   }

   const LoadType surplus_per_effective_descendent =
      (load_for_descendents - ideal_export_to_children) / num_effective_des;

   if (d_print_steps) {
      tbox::plog << "load_for_me = " << load_for_me
                 << ",  load_for_descendents = " << load_for_descendents
                 << ",  num_effective_des = " << num_effective_des
                 << ",  ideal_export_to_children = " << ideal_export_to_children
                 << ",  surplus_per_effective_descendent " << surplus_per_effective_descendent
                 << std::endl;
   }

   return surplus_per_effective_descendent;
}

/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void
TreeLoadBalancer::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& samrai_mpi)
{
   if (samrai_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      TBOX_ERROR(d_object_name << "::setSAMRAI_MPI error: Given\n"
                               << "communicator is invalid.");
   }

   if (d_mpi_is_dupe) {
      d_mpi.freeCommunicator();
   }

   // Enable private communicator.
   d_mpi.dupCommunicator(samrai_mpi);
   d_mpi_is_dupe = true;

   d_mca.setSAMRAI_MPI(d_mpi, true);
}

/*
 *************************************************************************
 * Set the MPI commuicator.
 *************************************************************************
 */
void
TreeLoadBalancer::freeMPICommunicator()
{
   if (d_mpi_is_dupe && d_mpi.getCommunicator() != MPI_COMM_NULL) {
      // Free the private communicator (if MPI has not been finalized).
      int flag;
      tbox::SAMRAI_MPI::Finalized(&flag);
      if (!flag) {
         d_mpi.freeCommunicator();
      }
   }
   d_mpi.setCommunicator(tbox::SAMRAI_MPI::commNull);
   d_mpi_is_dupe = false;
}

/*
 *************************************************************************
 * Create cycle-dependent groups to partition within.  Cycle fraction
 * is how far along we are in the number of cycles.  A value of 1 means
 * the last cycle.  Group size grows exponentially with cycle fraction
 * and includes all processes when cycle fraction is 1 (thus partitioning
 * is done across all prcesses).
 *
 * We chose this formula:
 * Group size = nprocs^cycle_fraction
 * thus number of groups = nprocs^(1-cycle_fraction)
 *
 * To avoid tiny groups that cannot effectively spread out severely
 * unbalanced load distribution, we round the number of groups down
 * and group size up.
 *************************************************************************
 */
void
TreeLoadBalancer::createBalanceRankGroupBasedOnCycles(
   tbox::RankGroup& rank_group,
   int& number_of_groups,
   int& group_num,
   double cycle_fraction) const
{
   number_of_groups =
      static_cast<int>(pow(static_cast<double>(d_mpi.getSize()),
                          1.0 - cycle_fraction));

   /*
    * All groups will have a base population count of
    * d_mpi.getSize()/number_of_groups.  The remainder from the
    * integer division is distributed to a subset of groups, starting
    * from group 0, so these groups will have one more than the base.
    */
   const int base_group_size = d_mpi.getSize() / number_of_groups;

   const int first_base_sized_group = d_mpi.getSize() % number_of_groups;

   const int first_rank_in_base_sized_group =
      first_base_sized_group * (1 + base_group_size);

   if (d_mpi.getRank() < first_rank_in_base_sized_group) {

      group_num = d_mpi.getRank() / (1 + base_group_size);

      const int group_first_rank = group_num * (1 + base_group_size);

      rank_group.setMinMax(group_first_rank,
         group_first_rank + base_group_size);
   } else {

      group_num = first_base_sized_group
         + (d_mpi.getRank() - first_rank_in_base_sized_group) / base_group_size;

      const int group_first_rank = first_rank_in_base_sized_group
         + (group_num - first_base_sized_group) * base_group_size;

      rank_group.setMinMax(group_first_rank,
         group_first_rank + base_group_size - 1);
   }
}

/*
 *************************************************************************
 * Set up the asynchronous communication objects for the process tree
 * containing ranks defined by the RankGroup.
 *
 * The process tree lay-out is defined by the BalancedDepthFirstTree
 * class, thus defining parent and children of the local process.
 * This method sets the AsyncCommPeer objects for communication with
 * children and parent.
 *************************************************************************
 */
void
TreeLoadBalancer::setupAsyncCommObjects(
   tbox::AsyncCommStage& child_stage,
   tbox::AsyncCommPeer<char> *& child_comms,
   tbox::AsyncCommStage& parent_stage,
   tbox::AsyncCommPeer<char> *& parent_comm,
   const tbox::RankGroup& rank_group) const
{

   child_comms = parent_comm = 0;

   const int num_children = d_rank_tree->getNumberOfChildren();

   if (num_children > 0) {

      child_comms = new tbox::AsyncCommPeer<char>[num_children];

      for (int child_num = 0; child_num < num_children; ++child_num) {

         const int child_rank_in_grp = d_rank_tree->getChildRank(child_num);
         const int child_true_rank = rank_group.getMappedRank(child_rank_in_grp);

         child_comms[child_num].initialize(&child_stage);
         child_comms[child_num].setPeerRank(child_true_rank);
         child_comms[child_num].setMPI(d_mpi);
         child_comms[child_num].setMPITag(TreeLoadBalancer_LOADTAG0,
            TreeLoadBalancer_LOADTAG1);
         child_comms[child_num].limitFirstDataLength(
            TreeLoadBalancer_FIRSTDATALEN);
      }
   }

   if (d_rank_tree->getParentRank() != tbox::RankTreeStrategy::getInvalidRank()) {

      const int parent_rank_in_grp = d_rank_tree->getParentRank();
      int parent_true_rank = rank_group.getMappedRank(parent_rank_in_grp);

      parent_comm = new tbox::AsyncCommPeer<char>;
      parent_comm->initialize(&parent_stage);
      parent_comm->setPeerRank(parent_true_rank);
      parent_comm->setMPI(d_mpi);
      parent_comm->setMPITag(TreeLoadBalancer_LOADTAG0,
         TreeLoadBalancer_LOADTAG1);
      parent_comm->limitFirstDataLength(
         TreeLoadBalancer_FIRSTDATALEN);

   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::destroyAsyncCommObjects(
   tbox::AsyncCommPeer<char> *& child_comms,
   tbox::AsyncCommPeer<char> *& parent_comm) const
{
   if (d_mpi.getSize() == 1) {
      TBOX_ASSERT(child_comms == 0);
      TBOX_ASSERT(parent_comm == 0);
   } else {
      if (child_comms != 0) {
         delete[] child_comms;
      }
      if (parent_comm != 0) {
         delete parent_comm;
      }
      child_comms = parent_comm = 0;
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::LoadType
TreeLoadBalancer::computeLocalLoad(
   const hier::BoxLevel& box_level) const
{
   t_compute_local_load->start();
   double load = 0.0;
   const hier::BoxContainer& boxes = box_level.getBoxes();
   for (hier::BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      double box_load = computeLoad(*ni);
      load += box_load;
   }
   t_compute_local_load->stop();
   return static_cast<LoadType>(load);
}

/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void
TreeLoadBalancer::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      d_print_steps = input_db->getBoolWithDefault("DEV_print_steps", false);
      d_check_connectivity =
         input_db->getBoolWithDefault("DEV_check_connectivity", d_check_connectivity);
      d_check_map =
         input_db->getBoolWithDefault("DEV_check_map", d_check_map);

      d_summarize_map = input_db->getBoolWithDefault("DEV_summarize_map",
            d_summarize_map);

      d_report_load_balance = input_db->getBoolWithDefault(
            "DEV_report_load_balance", d_report_load_balance);
      d_barrier_before = input_db->getBoolWithDefault("DEV_barrier_before",
            d_barrier_before);
      d_barrier_after = input_db->getBoolWithDefault("DEV_barrier_after",
            d_barrier_after);

      d_max_spread_procs =
         input_db->getIntegerWithDefault("max_spread_procs",
            d_max_spread_procs);

      d_flexible_load_tol =
         input_db->getDoubleWithDefault("flexible_load_tolerance",
            d_flexible_load_tol);

      d_allow_box_breaking =
         input_db->getBoolWithDefault("DEV_allow_box_breaking",
            d_allow_box_breaking);

      d_voucher_mode =
         input_db->getBoolWithDefault("DEV_voucher_mode",
            d_voucher_mode);

      if (input_db->isInteger("tile_size")) {
         input_db->getIntegerArray("tile_size", &d_tile_size[0], d_tile_size.getDim().getValue());
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (!(d_tile_size[i] >= 1)) {
               TBOX_ERROR("TreeLoadBalancer tile_size must be >= 1 in all directions.\n"
                  << "Input tile_size is " << d_tile_size);
            }
         }
      }

      if (input_db->isDouble("artificial_minimum_load")) {
         d_artificial_minimum =
            input_db->getDoubleVector("artificial_minimum_load");
      }
   }
}

/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void
TreeLoadBalancer::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      int flag;
      tbox::SAMRAI_MPI::Status mpi_status;
      int mpi_err = d_mpi.Iprobe(MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            &flag,
            &mpi_status);
      if (mpi_err != MPI_SUCCESS) {
         TBOX_ERROR("Error probing for possible lost messages.");
      }
      if (flag == true) {
         int count = -1;
         mpi_err = tbox::SAMRAI_MPI::Get_count(&mpi_status, MPI_INT, &count);
         TBOX_ERROR(
            "Library error!\n"
            << "TreeLoadBalancer detected before or\n"
            << "after using a private communicator that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n"
            << "current tags: "
            << ' ' << TreeLoadBalancer_LOADTAG0 << ' '
            << TreeLoadBalancer_LOADTAG1
            );
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancer::setTimers()
{
   /*
    * The first constructor gets timers from the TimerManager.
    * and sets up their deallocation.
    */
   if (!t_load_balance_box_level) {
      t_load_balance_box_level = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::loadBalanceBoxLevel()");

      t_get_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_map");
      t_use_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::use_map");

      t_constrain_size = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constrain_size");

      t_distribute_load_across_rank_group = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::distributeLoadAcrossRankGroup()");

      t_assign_to_local_and_populate_maps = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::assign_to_local_and_populate_maps");

      t_compute_local_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::computeLocalLoad");
      t_compute_global_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_global_load");
      t_compute_tree_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load");

      const int max_cycles_to_time = 5;
      t_compute_tree_load_for_cycle.resize(max_cycles_to_time, std::shared_ptr<tbox::Timer>());
      t_load_balance_for_cycle.resize(max_cycles_to_time, std::shared_ptr<tbox::Timer>());
      for (int i = 0; i < max_cycles_to_time; ++i) {
         t_compute_tree_load_for_cycle[i] = tbox::TimerManager::getManager()->
            getTimer(d_object_name + "::compute_tree_load_for_cycle["
               + tbox::Utilities::intToString(i) + "]");
         t_load_balance_for_cycle[i] = tbox::TimerManager::getManager()->
            getTimer(d_object_name + "::load_balance_for_cycle["
               + tbox::Utilities::intToString(i) + "]");
      }

      t_local_load_moves = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::local_load_moves");

      t_send_load_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_children");
      t_send_load_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_parent");
      t_get_load_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_children");
      t_get_load_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_parent");

      t_post_load_distribution_barrier = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::post_load_distribution_barrier");

      t_report_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::report_loads");

      t_finish_sends = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::finish_sends");

      t_child_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_send_wait");
      t_child_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_recv_wait");
      t_parent_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_send_wait");
      t_parent_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_recv_wait");

      t_barrier_before = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_before");
      t_barrier_after = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_after");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::BranchData::BranchData(
   const PartitioningParams& pparams,
   const TransitLoad& transit_load_prototype):
   d_num_procs(0),
   d_branch_load_current(0),
   d_branch_load_ideal(-1),
   d_branch_load_upperlimit(-1),
   d_eff_num_procs(0),
   d_eff_load_current(0),
   d_eff_load_ideal(-1),
   d_eff_load_upperlimit(-1),
   d_shipment(transit_load_prototype.clone()),
   d_wants_work_from_parent(false),
   d_pparams(&pparams),
   d_print_steps(false)
{
}

/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancer::BranchData::BranchData(
   const BranchData& other):
   d_num_procs(other.d_num_procs),
   d_branch_load_current(other.d_branch_load_current),
   d_branch_load_ideal(other.d_branch_load_ideal),
   d_branch_load_upperlimit(other.d_branch_load_upperlimit),
   d_eff_num_procs(other.d_eff_num_procs),
   d_eff_load_current(other.d_eff_load_current),
   d_eff_load_ideal(other.d_eff_load_ideal),
   d_eff_load_upperlimit(other.d_eff_load_upperlimit),
   d_shipment(other.d_shipment->clone()),
   d_wants_work_from_parent(other.d_wants_work_from_parent),
   d_pparams(other.d_pparams),
   t_pack_load(other.t_pack_load),
   t_unpack_load(other.t_unpack_load),
   d_print_steps(other.d_print_steps)
{
}

/*
 *************************************************************************
 * Set the starting load, which includes only the local processor's
 * contribution.  The parent's and children's contributions would be
 * added later.
 *************************************************************************
 */
void TreeLoadBalancer::BranchData::setStartingLoad(
   LoadType ideal,
   LoadType current,
   LoadType upperlimit)
{
   d_num_procs = 1;
   d_branch_load_ideal = ideal;
   d_branch_load_current = current;
   d_branch_load_upperlimit = upperlimit;

   d_eff_num_procs = d_num_procs;
   d_eff_load_ideal = d_branch_load_ideal;
   d_eff_load_current = d_branch_load_current;
   d_eff_load_upperlimit = d_branch_load_upperlimit;
}

/*
 *************************************************************************
 * Incorporate a child branch's data into this branch.
 *************************************************************************
 */
void
TreeLoadBalancer::BranchData::incorporateChild(
   const BranchData& child)
{
   d_num_procs += child.d_num_procs;
   d_branch_load_current += child.d_branch_load_current;
   d_branch_load_upperlimit += child.d_branch_load_upperlimit;
   d_branch_load_ideal += child.d_branch_load_ideal;

   if (child.d_wants_work_from_parent) {
      d_eff_num_procs += child.d_eff_num_procs;
      d_eff_load_current += child.d_eff_load_current;
      d_eff_load_upperlimit += child.d_eff_load_upperlimit;
      d_eff_load_ideal += child.d_eff_load_ideal;
   }

   d_branch_load_current += child.d_shipment->getSumLoad();
   d_eff_load_current += child.d_shipment->getSumLoad();
}

/*
 *************************************************************************
 * Incorporate a child branch's data into this branch.
 *************************************************************************
 */
TreeLoadBalancer::LoadType TreeLoadBalancer::BranchData::adjustOutboundLoad(
   TransitLoad& reserve,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load)
{
   LoadType actual_transfer = 0;

   if (low_load > d_pparams->getLoadComparisonTol()) {

      if (d_print_steps) {
         tbox::plog << "BranchData::adjustOutboundLoad adjusting shipment to "
                    << ideal_load << " [" << low_load << ", " << high_load << "]\n";
      }

      actual_transfer = d_shipment->getSumLoad();

      d_shipment->adjustLoad(
         reserve,
         ideal_load,
         low_load,
         high_load);

      actual_transfer = d_shipment->getSumLoad() - actual_transfer;

      d_branch_load_current -= d_shipment->getSumLoad();
      d_eff_load_current -= d_shipment->getSumLoad();

      if (d_print_steps) {
         tbox::plog << "BranchData::adjustOutboundLoad: Assigned to shipment ";
         d_shipment->recursivePrint(tbox::plog);
         tbox::plog << std::endl;
         tbox::plog << "Remaining in reserve: ";
         reserve.recursivePrint(tbox::plog, "  ", 0);
      }
   }

   return actual_transfer;
}

/*
 *************************************************************************
 * We could empty d_shipment at the end of this method, because it
 * is no longer essential.  We keep its contents only for diagnostics.
 *************************************************************************
 */
void TreeLoadBalancer::BranchData::moveInboundLoadToReserve(
   TransitLoad& reserve)
{
   reserve.insertAll(*d_shipment);
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::BranchData::packDataToParent(
   tbox::MessageStream& msg) const
{
   t_pack_load->start();
   msg << d_num_procs;
   msg << d_branch_load_current;
   msg << d_branch_load_ideal;
   msg << d_branch_load_upperlimit;
   msg << d_eff_num_procs;
   msg << d_eff_load_current;
   msg << d_eff_load_ideal;
   msg << d_eff_load_upperlimit;
   msg << d_wants_work_from_parent;
   d_shipment->putToMessageStream(msg);
   t_pack_load->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::BranchData::unpackDataFromChild(
   tbox::MessageStream& msg)
{
   t_unpack_load->start();
   msg >> d_num_procs;
   msg >> d_branch_load_current;
   msg >> d_branch_load_ideal;
   msg >> d_branch_load_upperlimit;
   msg >> d_eff_num_procs;
   msg >> d_eff_load_current;
   msg >> d_eff_load_ideal;
   msg >> d_eff_load_upperlimit;
   msg >> d_wants_work_from_parent;
   d_shipment->getFromMessageStream(msg);
   t_unpack_load->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::BranchData::packDataToChild(
   tbox::MessageStream& msg) const
{
   t_pack_load->start();
   d_shipment->putToMessageStream(msg);
   t_pack_load->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancer::BranchData::unpackDataFromParentAndIncorporate(
   tbox::MessageStream& msg)
{
   t_unpack_load->start();
   d_shipment->getFromMessageStream(msg);
   d_branch_load_current += d_shipment->getSumLoad();
   d_eff_load_current += d_shipment->getSumLoad();
   t_unpack_load->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancer::BranchData::setTimerPrefix(
   const std::string& timer_prefix)
{
   t_pack_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::pack_load");
   t_unpack_load = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::unpack_load");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancer::printStatistics(
   std::ostream& output_stream) const
{
   if (d_load_stat.empty()) {
      output_stream << "No statistics for TreeLoadBalancer.\n";
   } else {
      BalanceUtilities::reduceAndReportLoadBalance(
         d_load_stat,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         output_stream);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancer::BranchData::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int detail_depth) const
{
   os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
   os.precision(6);
   os << border
      << "Full nproc = " << d_num_procs
      << "   current = " << d_branch_load_current
      << "   ideal = " << d_branch_load_ideal
      << "   ratio = " << (d_branch_load_current / d_branch_load_ideal)
      << "   avg = " << (d_branch_load_current / d_num_procs)
      << "   upperlimit = " << d_branch_load_upperlimit
      << "   surplus = " << surplus()
      << "   excess =  " << excess()
      << '\n' << border
      << "Effective nproc = " << d_eff_num_procs
      << "   current = " << d_eff_load_current
      << "   ideal = " << d_eff_load_ideal
      << "   ratio = " << (d_eff_load_current / d_eff_load_ideal)
      << "   avg = " << (d_eff_load_current / d_eff_num_procs)
      << "   upperlimit = " << d_eff_load_upperlimit
      << "   surplus = " << effSurplus()
      << "   excess =  " << effExcess()
      << '\n' << border
      << "   wants work from parent = " << d_wants_work_from_parent
      << '\n' << border
      << "   shipment: ";
   d_shipment->recursivePrint(os, border + "   ", detail_depth - 1);
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
