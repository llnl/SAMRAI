/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using a "cascade" algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_CascadePartitioner_C
#define included_mesh_CascadePartitioner_C

#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/mesh/VoucherTransitLoad.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
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

const int CascadePartitioner::CascadePartitioner_LOADTAG0;
const int CascadePartitioner::CascadePartitioner_LOADTAG1;
const int CascadePartitioner::CascadePartitioner_FIRSTDATALEN;

const int CascadePartitioner::s_default_data_id = -1;

/*
 *************************************************************************
 * CascadePartitioner constructor.
 *************************************************************************
 */

CascadePartitioner::CascadePartitioner(
   const tbox::Dimension& dim,
   const std::string& name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(name),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mpi_is_dupe(false),
   d_workload_data_id(0),
   d_master_workload_data_id(s_default_data_id),
   d_tile_size(dim, 1),
   d_max_spread_procs(500),
   d_limit_supply_to_surplus(true),
   d_reset_obligations(true),
   d_flexible_load_tol(0.05),
   d_artificial_minimum(1,1.0),
   d_use_vouchers(false),
   d_mca(),
   // Shared data.
   d_workload_level(),
   d_balance_box_level(0),
   d_balance_to_reference(0),
   d_global_work_sum(-1),
   d_global_work_avg(-1),
   d_num_initial_owners(0),
   d_local_load(0),
   d_shipment(0),
   d_comm_stage(),
   // Performance evaluation and diagnostics.
   d_barrier_before(false),
   d_barrier_after(false),
   d_report_load_balance(false),
   d_summarize_map(false),
   d_print_steps(false),
   d_print_child_steps(false),
   d_check_connectivity(false),
   d_check_map(false)
{
   for (int i = 0; i < 4; ++i) d_comm_peer[i].initialize(&d_comm_stage);

   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);
   setTimers();
   d_comm_stage.setCommunicationWaitTimer(t_communication_wait);
   d_mca.setTimerPrefix(d_object_name);
}

/*
 *************************************************************************
 * CascadePartitioner constructor.
 *************************************************************************
 */

CascadePartitioner::~CascadePartitioner()
{
   freeMPICommunicator();
}

/*
 *************************************************************************
 * Accessory functions to get/set load balancing parameters.
 *************************************************************************
 */

bool
CascadePartitioner::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}

/*
 **************************************************************************
 **************************************************************************
 */
void
CascadePartitioner::setWorkloadPatchDataIndex(
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
         for (int i = asize; i < level_number - 1; i++) {
            d_workload_data_id[i] = d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = data_id;
      }
   } else {
      d_master_workload_data_id = data_id;
      for (int ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}

/*
 *************************************************************************
 * This method implements the abstract LoadBalanceStrategy interface.
 *************************************************************************
 */
void
CascadePartitioner::loadBalanceBoxLevel(
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
   NULL_USE(domain_box_level);
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
      if (static_cast<unsigned int>(level_number) < d_artificial_minimum.size()) {
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
         TBOX_ERROR("CascadePartitioner::loadBalanceBoxLevel:\n"
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
                 << balance_box_level.format("  ", 2)
                 << std::flush;
   }

   // Set effective_cut_factor to least common multiple of cut_factor and d_tile_size.
   const size_t nblocks = hierarchy->getGridGeometry()->getNumberBlocks();
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
         tbox::plog << d_object_name << "::loadBalanceBoxLevel"
                    << "  effective_cut_factor=" << effective_cut_factor
                    << "  d_tile_size=" << d_tile_size
                    << std::endl;
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

   d_workload_level.reset();
   t_load_balance_box_level->start();

   d_pparams = std::make_shared<PartitioningParams>(
         *balance_box_level.getGridGeometry(),
         balance_box_level.getRefinementRatio(),
         min_size, max_size, bad_interval, effective_cut_factor,
         minimum_cells,
         artificial_minimum,
         d_flexible_load_tol);

   LoadType local_load = computeLocalLoad(balance_box_level);

   globalWorkReduction(local_load,
                       (balance_box_level.getLocalNumberOfBoxes() != 0));

   d_global_work_avg = d_global_work_sum / rank_group.size();

   // Run the partitioning algorithm.
   partitionByCascade(
      balance_box_level,
      balance_to_reference,
      d_use_vouchers);

   t_load_balance_box_level->stop();

   int wrk_indx = getWorkloadDataId(level_number);

   /*
    * Do non-uniform load balance if a workload data id has been registered
    * and this is not a new finest level of the hierarchy.
    */
   if ((wrk_indx >= 0) && (hierarchy->getNumberOfLevels() > level_number)) {

      d_workload_level =
         std::make_shared<hier::PatchLevel>(balance_box_level,
                                              hierarchy->getGridGeometry(),
                                              hierarchy->getPatchDescriptor());

      d_workload_level->setLevelNumber(level_number);

      d_pparams->setWorkloadDataId(wrk_indx);
      d_pparams->setWorkloadPatchLevel(d_workload_level);

      /*
       * Set up workload_to_reference and reference_to_workload.  Since
       * d_workload_level is based on balance_box_level, the new Connectors
       * are effectively copies of balance_to_reference and its transpose.
       */
      std::shared_ptr<hier::Connector> workload_to_reference(
         std::make_shared<hier::Connector>(
            *d_workload_level->getBoxLevel(),
            balance_to_reference->getHead(),
            balance_to_reference->getConnectorWidth()));

      for (hier::Connector::ConstNeighborhoodIterator ei =
           balance_to_reference->begin();
           ei != balance_to_reference->end(); ++ei) {
         const hier::BoxId& box_id = *ei;
         for (hier::Connector::ConstNeighborIterator na =
              balance_to_reference->begin(ei);
              na != balance_to_reference->end(ei); ++na) {
            workload_to_reference->insertLocalNeighbor(*na, box_id);
         }
      }

      std::shared_ptr<hier::Connector> reference_to_workload(
         std::make_shared<hier::Connector>(
            balance_to_reference->getHead(),
            *d_workload_level->getBoxLevel(),
            balance_to_reference->getTranspose().getConnectorWidth()));

      for (hier::Connector::ConstNeighborhoodIterator ti =
           balance_to_reference->getTranspose().begin();
           ti != balance_to_reference->getTranspose().end(); ++ti) {
         const hier::BoxId& box_id = *ti;
         for (hier::Connector::ConstNeighborIterator ta =
              balance_to_reference->getTranspose().begin(ti);
              ta != balance_to_reference->getTranspose().end(ti); ++ta) {
            reference_to_workload->insertLocalNeighbor(*ta, box_id);
         }
      }

      /*
       * Cache the Connectors before calling setTranspose.
       */
      d_workload_level->cacheConnector(workload_to_reference);
      reference_to_workload->getBase().cacheConnector(reference_to_workload);
      reference_to_workload->setTranspose(workload_to_reference.get(), false);

      /*
       * Find the Connectors between the current level of the hierarchy and
       * the reference level.
       */
      std::shared_ptr<hier::PatchLevel> current_level(
         hierarchy->getPatchLevel(level_number));

      const hier::Connector& current_to_reference = 
         current_level->getBoxLevel()->findConnector(
            workload_to_reference->getHead(),
            hierarchy->getRequiredConnectorWidth(
               level_number,
               tbox::MathUtilities<int>::Max(level_number-1, 0)),
            hier::CONNECTOR_CREATE,
            true);

      const hier::Connector& reference_to_current =
         workload_to_reference->getHead().findConnector(
            *current_level->getBoxLevel(),
            hierarchy->getRequiredConnectorWidth(
               tbox::MathUtilities<int>::Max(level_number-1, 0),
               level_number),
            hier::CONNECTOR_CREATE,
            true);

      /*
       * All of the above Connector work was so that we can call these
       * bridge operations to connect the current and workload levels.
       */
      hier::OverlapConnectorAlgorithm oca;
      std::shared_ptr<hier::Connector> current_to_workload;
      oca.bridgeWithNesting(
         current_to_workload,
         current_to_reference,
         *reference_to_workload,
         hier::IntVector::getZero(d_dim),
         hier::IntVector::getZero(d_dim),
         hier::IntVector::getOne(d_dim),
         false);
      current_level->cacheConnector(current_to_workload);

      std::shared_ptr<hier::Connector> workload_to_current;
      oca.bridgeWithNesting(
         workload_to_current,
         *workload_to_reference,
         reference_to_current,                
         hier::IntVector::getZero(d_dim),
         hier::IntVector::getZero(d_dim),
         hier::IntVector::getOne(d_dim),
         false);
      d_workload_level->cacheConnector(workload_to_current);

      /*
       * Build and use a RefineSchedule to communicate workload data
       * from the current level to d_workload_level.
       */
      d_workload_level->allocatePatchData(wrk_indx);
   
      xfer::RefineAlgorithm fill_work_algorithm;

      std::shared_ptr<hier::RefineOperator> work_refine_op(
         std::make_shared<pdat::CellDoubleConstantRefine>());

      fill_work_algorithm.registerRefine(wrk_indx,
         wrk_indx,
         wrk_indx,
         work_refine_op);

      fill_work_algorithm.createSchedule(d_workload_level,
         current_level,
         level_number - 1,
         hierarchy)->fillData(0.0);

      t_load_balance_box_level->start();

      /*
       * Compute workloads for each box and run the partitioning algorithm
       */
      local_load =
         computeNonUniformWorkLoad(*d_workload_level);

      globalWorkReduction(local_load,
                          (balance_box_level.getLocalNumberOfBoxes() != 0));

      d_global_work_avg = d_global_work_sum / rank_group.size();

      /*
       * Run partitioning algorithm again, this time taking into account
       * the computed workloads.  This call always uses vouchers.
       */
      partitionByCascade(
         balance_box_level,
         balance_to_reference,
         true);

      d_workload_level.reset();
      t_load_balance_box_level->stop();

   }

   /*
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */

   hier::IntVector max_intvector(d_dim, tbox::MathUtilities<int>::getMax());
   if (max_size != max_intvector) {

      BalanceUtilities::constrainMaxBoxSizes(
         balance_box_level,
         balance_to_reference ? &balance_to_reference->getTranspose() : 0,
         *d_pparams);

      if (d_print_steps) {
         tbox::plog << " CascadePartitioner completed constraining box sizes."
                    << "\n";
      }

   }

   /*
    * Finished load balancing.  Clean up and wrap up.
    */

   d_pparams.reset();

   local_load = computeLocalLoad(balance_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_box_level.getBoxes().size()));

   if (d_print_steps) {
      tbox::plog << "Post balanced:\n" << balance_box_level.format("", 2)
                 << std::flush;
   }

   if (d_report_load_balance) {
      tbox::plog << d_object_name << "::loadBalanceBoxLevel results:" << std::endl;
      BalanceUtilities::reduceAndReportLoadBalance(
         std::vector<double>(1, local_load), balance_box_level.getMPI());
   }

   if (d_check_connectivity && balance_to_reference) {
      hier::Connector& reference_to_balance = balance_to_reference->getTranspose();
      tbox::plog << "CascadePartitioner checking balance-reference connectivity."
                 << std::endl;
      int errs = 0;
      if (reference_to_balance.checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in reference_to_balance!" << std::endl;
      }
      if (balance_to_reference->checkOverlapCorrectness(false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_reference!" << std::endl;
      }
      if (reference_to_balance.checkTransposeCorrectness(*balance_to_reference)) {
         ++errs;
         tbox::perr << "Error found in balance-reference transpose!" << std::endl;
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found.\n"
            << "reference_box_level:\n" << reference_to_balance.getBase().format("", 2)
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "reference_to_balance:\n" << reference_to_balance.format("", 2)
            << "balance_to_reference:\n" << balance_to_reference->format("", 2));
      }
      tbox::plog << "CascadePartitioner checked balance-reference connectivity."
                 << std::endl;
   }

   assertNoMessageForPrivateCommunicator();
}

/*
 *************************************************************************
 * This method implements the cascade partitioner algorithm.
 *
 * It calls distributeLoad to do distribute the load then assigns the
 * distributed the loads to their new owners and update Connectors.
 *************************************************************************
 */
void
CascadePartitioner::partitionByCascade(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_reference,
   bool use_vouchers) const
{
   if (d_print_steps) {
      tbox::plog << d_object_name << "::partitionByCascade: entered" << std::endl;
   }

   std::shared_ptr<TransitLoad> local_load;
   std::shared_ptr<TransitLoad> shipment;
   if (use_vouchers) {
      local_load = std::make_shared<VoucherTransitLoad>(*d_pparams);
      shipment = std::make_shared<VoucherTransitLoad>(*d_pparams);
      d_pparams->setUsingVouchers(true);
   } else {
      local_load = std::make_shared<BoxTransitSet>(*d_pparams);
      shipment = std::make_shared<BoxTransitSet>(*d_pparams);
   }
   local_load->setAllowBoxBreaking(true);
   local_load->setTimerPrefix(d_object_name);
   shipment->setTimerPrefix(d_object_name);

   const double ideal_box_width = pow(d_global_work_avg, 1.0 / d_dim.getValue());
   local_load->setThresholdWidth(ideal_box_width);
   shipment->setThresholdWidth(ideal_box_width);

   if (d_pparams->getArtificialMinimumLoad() >
       d_pparams->getMinBoxSizeProduct()) {
      local_load->insertAllWithArtificialMinimum(
         balance_box_level.getBoxes(),
         d_pparams->getArtificialMinimumLoad());
   } else {
      local_load->insertAll(balance_box_level.getBoxes());
   }
   if (d_workload_level) {
      local_load->setWorkload(*d_workload_level,
         getWorkloadDataId(d_workload_level->getLevelNumber()));
   }

   // Set up temporaries shared with the process groups.
   d_balance_box_level = &balance_box_level;
   d_balance_to_reference = balance_to_reference;
   d_local_load = &(*local_load);
   d_shipment = &(*shipment);

   CascadePartitionerTree groups(*this);
   groups.distributeLoad();

   d_balance_box_level = 0;
   d_balance_to_reference = 0;
   d_local_load = 0;
   d_shipment = 0;
   d_global_work_sum = -1;
   d_global_work_avg = -1;
   d_num_initial_owners = 0;

   if (d_print_steps) {
      tbox::plog << d_object_name << "::partitionByCascade: leaving" << std::endl;
   }
}

/*
 *************************************************************************
 * Update connectors according to the current boxes in d_local_load.
 * Only point-to-point communication should be used.  If there are not
 * power-of-two processes in d_mpi and d_max_spread_procs is small enough
 * to cause intermediate calls to updateConnectors, it is possible that
 * a subset of processes will be calling updateConnectors().  So don't
 * use collective MPI communications in this method!
 *************************************************************************
 */
void CascadePartitioner::updateConnectors() const
{
   t_update_connectors->start();

   if (d_print_steps) {
      tbox::plog
      << d_object_name << "::updateConnectors constructing unbalanced<==>balanced." << std::endl;
   }

   /*
    * Initialize empty balanced_box_level and mappings so they are
    * ready to be populated.
    */
   hier::IntVector zero_vector(hier::IntVector::getZero(d_dim));
   hier::BoxLevel balanced_box_level(
      d_balance_box_level->getRefinementRatio(),
      d_balance_box_level->getGridGeometry(),
      d_balance_box_level->getMPI());
   hier::MappingConnector balanced_to_unbalanced(balanced_box_level,
                                                 *d_balance_box_level,
                                                 hier::IntVector::getZero(d_dim));
   hier::MappingConnector unbalanced_to_balanced(*d_balance_box_level,
                                                 balanced_box_level,
                                                 hier::IntVector::getZero(d_dim));
   unbalanced_to_balanced.setTranspose(&balanced_to_unbalanced, false);

   t_assign_to_local_and_populate_maps->start();
   d_local_load->assignToLocalAndPopulateMaps(
      balanced_box_level,
      balanced_to_unbalanced,
      unbalanced_to_balanced,
      d_flexible_load_tol,
      d_mpi);
   t_assign_to_local_and_populate_maps->stop();

   if (d_summarize_map) {
      tbox::plog << d_object_name << "::updateConnectors unbalanced--->balanced map:" << std::endl
                 << unbalanced_to_balanced.format("\t", 0)
                 << "Map statistics:" << std::endl << unbalanced_to_balanced.formatStatistics("\t")
                 << d_object_name << "::updateConnectors balanced--->unbalanced map:" << std::endl
                 << balanced_to_unbalanced.format("\t", 0)
                 << "Map statistics:" << std::endl << balanced_to_unbalanced.formatStatistics("\t")
                 << std::endl;
   }

   if (d_check_map) {
      if (unbalanced_to_balanced.findMappingErrors() != 0) {
         TBOX_ERROR(
            d_object_name << "::updateConnectors Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            d_object_name << "::updateConnectors Transpose errors found!");
      }
   }

   if (d_summarize_map) {
      tbox::plog << d_object_name << "::updateConnectors: unbalanced--->balanced map:" << std::endl
                 << unbalanced_to_balanced.format("\t", 0)
                 << "Map statistics:" << std::endl << unbalanced_to_balanced.formatStatistics("\t")
                 << d_object_name << "::updateConnectors: balanced--->unbalanced map:" << std::endl
                 << balanced_to_unbalanced.format("\t", 0)
                 << "Map statistics:" << std::endl << balanced_to_unbalanced.formatStatistics("\t")
                 << std::endl;
   }

   if (d_balance_to_reference && d_balance_to_reference->hasTranspose()) {
      if (d_print_steps) {
         tbox::plog
         << d_object_name << "::updateConnectors applying unbalanced<==>balanced." << std::endl;
      }
      t_use_map->start();
      d_mca.modify(
         d_balance_to_reference->getTranspose(),
         unbalanced_to_balanced,
         d_balance_box_level,
         &balanced_box_level);
      t_use_map->stop();
   } else {
      hier::BoxLevel::swap(*d_balance_box_level, balanced_box_level);
   }

   d_local_load->insertAllWithExistingLoads(d_balance_box_level->getBoxes());

   if (d_print_steps) {
      tbox::plog
      << d_object_name << "::updateConnectors leaving." << std::endl;
   }

   t_update_connectors->stop();
}

/*
 *************************************************************************
 * Set d_global_work_sum, d_local_work_max, d_num_initial_owners.
 *************************************************************************
 */
void CascadePartitioner::globalWorkReduction(
   LoadType local_work,
   bool has_any_load) const
{
   t_global_work_reduction->start();

   d_global_work_sum = d_local_work_max = local_work;
   d_num_initial_owners = static_cast<size_t>(has_any_load);

   if (d_mpi.getSize() > 1) {
      double dtmp[2], dtmp_sum[2], dtmp_max[2];

      dtmp[0] = local_work;
      dtmp[1] = static_cast<double>(d_num_initial_owners);

      d_mpi.Allreduce(dtmp, dtmp_sum, 2, MPI_DOUBLE, MPI_SUM);
      d_global_work_sum = dtmp_sum[0];
      d_num_initial_owners = static_cast<size_t>(dtmp_sum[1]);

      d_mpi.Allreduce(dtmp, dtmp_max, 1, MPI_DOUBLE, MPI_MAX);
      d_local_work_max = dtmp_max[0];
   }

   if (d_print_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      tbox::plog.precision(6);
      tbox::plog << d_object_name << "::globalWorkReduction"
                 << " d_local_work_max=" << d_local_work_max
                 << " d_global_work_sum=" << d_global_work_sum
                 << " (initially born on "
                 << d_num_initial_owners << " procs) across all "
                 << d_mpi.getSize()
                 << " procs, averaging " << d_global_work_sum / d_mpi.getSize()
                 << " or " << pow(d_global_work_sum / d_mpi.getSize(), 1.0 / d_dim.getValue())
                 << "^" << d_dim << " per proc." << std::endl;
   }

   t_global_work_reduction->stop();
}

/*
 *************************************************************************
 * Compute log-base-2 of integer, rounded up.
 *************************************************************************
 */
int CascadePartitioner::lgInt(int s) {
   int lg_s = 0;
   while ((1 << lg_s) < s) {
      ++lg_s;
   }
   return lg_s;
}

/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void
CascadePartitioner::setSAMRAI_MPI(
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

   d_mca.setSAMRAI_MPI(d_mpi);
}

/*
 *************************************************************************
 * Set the MPI commuicator.
 *************************************************************************
 */
void
CascadePartitioner::freeMPICommunicator()
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
 *************************************************************************
 */
CascadePartitioner::LoadType
CascadePartitioner::computeLocalLoad(
   const hier::BoxLevel& box_level) const
{
   double artificial_load = 1.0;
   if (d_pparams) {
      artificial_load = d_pparams->getArtificialMinimumLoad();
   }
   double load = 0.0;
   const hier::BoxContainer& boxes = box_level.getBoxes();
   for (hier::BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      double box_load = static_cast<double>(ni->size());
      if (box_load < artificial_load) {
         box_load = artificial_load;
      }
      load += box_load;
   }
   return static_cast<LoadType>(load);


}

CascadePartitioner::LoadType
CascadePartitioner::computeNonUniformWorkLoad(
   const hier::PatchLevel& patch_level) const
{
   double load = 0.0;
   for (hier::PatchLevel::iterator ip(patch_level.begin());
        ip != patch_level.end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      double patch_work =
         BalanceUtilities::computeNonUniformWorkload(patch,
            getWorkloadDataId(patch_level.getLevelNumber()),
            patch->getBox());

      load += patch_work;
   }
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
CascadePartitioner::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      d_print_steps =
         input_db->getBoolWithDefault("DEV_print_steps", d_print_steps);
      d_print_child_steps =
         input_db->getBoolWithDefault("DEV_print_child_steps", d_print_child_steps);
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

      d_use_vouchers =
         input_db->getBoolWithDefault("use_vouchers", false);

      d_limit_supply_to_surplus =
         input_db->getBoolWithDefault("DEV_limit_supply_to_surplus",
            d_limit_supply_to_surplus);

      d_reset_obligations =
         input_db->getBoolWithDefault("DEV_reset_obligations",
            d_reset_obligations);

      d_flexible_load_tol =
         input_db->getDoubleWithDefault("flexible_load_tolerance",
            d_flexible_load_tol);

      if (input_db->isInteger("tile_size")) {
         input_db->getIntegerArray("tile_size", &d_tile_size[0], d_tile_size.getDim().getValue());
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (!(d_tile_size[i] >= 1)) {
               TBOX_ERROR("CascadePartitioner tile_size must be >= 1 in all directions.\n"
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
CascadePartitioner::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      tbox::SAMRAI_MPI::Status mpi_status;
      if (d_mpi.hasReceivableMessage(&mpi_status)) {
         int count = -1;
         tbox::SAMRAI_MPI::Get_count(&mpi_status, MPI_INT, &count);
         TBOX_ERROR(
            "Library error!\n"
            << "CascadePartitioner detected before or\n"
            << "after using a private communicator that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n"
            << "current tags: "
            << ' ' << CascadePartitioner_LOADTAG0 << ' '
            << CascadePartitioner_LOADTAG1
            );
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
CascadePartitioner::setTimers()
{
   /*
    * The first constructor gets timers from the TimerManager.
    * and sets up their deallocation.
    */
   if (!t_load_balance_box_level) {
      t_load_balance_box_level = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::loadBalanceBoxLevel()");

      t_use_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::use_map");

      t_update_connectors = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::updateConnectors()");

      t_global_work_reduction = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::globalWorkReduction()");

      t_assign_to_local_and_populate_maps = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::assign_to_local_and_populate_maps");

      t_communication_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::communication_wait");

      // These timers are shared by CascadePartitionerTree.
      t_distribute_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::distributeLoad()");
      t_combine_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::combineChildren()");
      t_balance_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::balanceChildren()");
      t_supply_work = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::supply_work");
      t_send_shipment = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::sendShipment()");
      t_receive_and_unpack_supplied_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::receiveAndUnpackSuppliedLoad()");

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
CascadePartitioner::printStatistics(
   std::ostream& output_stream) const
{
   if (d_load_stat.empty()) {
      output_stream << "No statistics for CascadePartitioner.\n";
   } else {
      BalanceUtilities::reduceAndReportLoadBalance(
         d_load_stat,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         output_stream);
   }
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
