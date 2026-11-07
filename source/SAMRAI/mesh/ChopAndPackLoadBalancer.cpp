/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Load balance routines for uniform and non-uniform workloads.
 *
 ************************************************************************/
#define ChopAndPackLoadBalancer_MARKLOADFORPOSTPROCESSING

#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstdlib>
#include <fstream>
#include <list>

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * Constructors and destructor for ChopAndPackLoadBalancer.
 *
 *************************************************************************
 */

ChopAndPackLoadBalancer::ChopAndPackLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(name),
   d_processor_layout_specified(false),
   d_processor_layout(d_dim),
   d_ignore_level_box_union_is_single_box(false),
   d_master_workload_data_id(-1),
   d_master_max_workload_factor(1.0),
   d_master_workload_tolerance(0.0),
   d_master_bin_pack_method("SPATIAL"),
   d_tile_size(dim, 1)
{
   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);
   setupTimers();
}

ChopAndPackLoadBalancer::ChopAndPackLoadBalancer(
   const tbox::Dimension& dim,
   const std::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name("ChopAndPackLoadBalancer"),
   d_processor_layout_specified(false),
   d_processor_layout(hier::IntVector::getZero(d_dim)),
   d_ignore_level_box_union_is_single_box(false),
   d_master_workload_data_id(-1),
   d_master_max_workload_factor(1.0),
   d_master_workload_tolerance(0.0),
   d_master_bin_pack_method("SPATIAL"),
   d_tile_size(dim, 1)

{
   getFromInput(input_db);
   setupTimers();
}

ChopAndPackLoadBalancer::~ChopAndPackLoadBalancer()
{
}

/*
 *************************************************************************
 *
 * Accessory functions to get/set load balancing parameters.
 *
 *************************************************************************
 */

bool
ChopAndPackLoadBalancer::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}

void
ChopAndPackLoadBalancer::setMaxWorkloadFactor(
   double factor,
   int level_number)
{
   TBOX_ASSERT(factor > 0.0);
   if (level_number >= 0) {
      int asize = static_cast<int>(d_max_workload_factor.size());
      if (asize < level_number + 1) {
         d_max_workload_factor.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; ++i) {
            d_max_workload_factor[i] = d_master_max_workload_factor;
         }
         d_max_workload_factor[level_number] = factor;
      }
   } else {
      d_master_max_workload_factor = factor;
      for (int ln = 0; ln < static_cast<int>(d_max_workload_factor.size()); ++ln) {
         d_max_workload_factor[ln] = d_master_max_workload_factor;
      }
   }
}

void
ChopAndPackLoadBalancer::setWorkloadTolerance(
   double tolerance,
   int level_number)
{
   TBOX_ASSERT(tolerance > 0.0);
   if (level_number >= 0) {
      int asize = static_cast<int>(d_workload_tolerance.size());
      if (asize < level_number + 1) {
         d_workload_tolerance.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; ++i) {
            d_workload_tolerance[i] = d_master_workload_tolerance;
         }
         d_workload_tolerance[level_number] = tolerance;
      }
   } else {
      d_master_workload_tolerance = tolerance;
      for (int ln = 0; ln < static_cast<int>(d_workload_tolerance.size()); ++ln) {
         d_workload_tolerance[ln] = d_master_workload_tolerance;
      }
   }
}

void
ChopAndPackLoadBalancer::setWorkloadPatchDataIndex(
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

void
ChopAndPackLoadBalancer::setUniformWorkload(
   int level_number)
{
   if (level_number >= 0) {
      int asize = static_cast<int>(d_workload_data_id.size());
      if (asize < level_number + 1) {
         d_workload_data_id.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; ++i) {
            d_workload_data_id[i] = d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = -1;
      }
   } else {
      d_master_workload_data_id = -1;
      for (int ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ++ln) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}

void
ChopAndPackLoadBalancer::setBinPackMethod(
   const std::string& method,
   int level_number)
{

   if (!(method == "GREEDY" ||
         method == "SPATIAL")) {
      TBOX_ERROR(
         d_object_name << " error: "
                       << "\n   String " << method
                       << " passed to setBinPackMethod()"
                       << " is not a valid method string identifier."
                       << std::endl);

   }

   if (level_number >= 0) {
      int asize = static_cast<int>(d_bin_pack_method.size());
      if (asize < level_number + 1) {
         d_bin_pack_method.resize(level_number + 1);
         for (int i = asize; i < level_number - 1; ++i) {
            d_bin_pack_method[i] = d_master_bin_pack_method;
         }
         d_bin_pack_method[level_number] = method;
      }
   } else {
      d_master_bin_pack_method = method;
      for (int ln = 0; ln < static_cast<int>(d_bin_pack_method.size()); ++ln) {
         d_bin_pack_method[ln] = d_master_bin_pack_method;
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   TBOX_ASSERT(!balance_to_anchor || balance_to_anchor->hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY6(d_dim,
      balance_box_level,
      min_size,
      max_size,
      domain_box_level,
      bad_interval,
      cut_factor);
   NULL_USE(rank_group);

   t_load_balance_box_level->start();

   hier::IntVector actual_max_size = max_size;
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (actual_max_size(d) < 0) {
         actual_max_size(d) = tbox::MathUtilities<int>::getMax();
      }
   }

   const size_t nblocks = balance_box_level.getGridGeometry()->getNumberBlocks();

   // Set effective_cut_factor to least common multiple of cut_factor and d_tile_size.
   hier::IntVector effective_cut_factor(cut_factor, nblocks);
   if (d_tile_size != hier::IntVector::getOne(d_dim)) {
      for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
         for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
            while (effective_cut_factor(b,d) / d_tile_size[d] * d_tile_size[d] !=
                   effective_cut_factor(b,d)) {
               effective_cut_factor(b,d) += cut_factor[d];
            }
         }
      }
   }

   t_get_global_boxes->barrierAndStart();
   hier::BoxLevel globalized_input_box_level(balance_box_level);
   globalized_input_box_level.setParallelState(hier::BoxLevel::GLOBALIZED);
   t_get_global_boxes->stop();

   hier::BoxContainer in_boxes;
   const hier::BoxContainer globalized_input_boxes(
      globalized_input_box_level.getGlobalBoxes());
   for (hier::RealBoxConstIterator bi(globalized_input_boxes.realBegin());
        bi != globalized_input_boxes.realEnd(); ++bi) {
      in_boxes.pushBack(*bi);
   }



   hier::BoxContainer physical_domain;
   domain_box_level.getGlobalBoxes(physical_domain);

   hier::BoxContainer out_boxes;
   hier::ProcessorMapping mapping;

   loadBalanceBoxes(
      out_boxes,
      mapping,
      in_boxes,
      hierarchy,
      level_number,
      physical_domain,
      balance_box_level.getRefinementRatio(),
      min_size,
      actual_max_size,
      effective_cut_factor,
      bad_interval);

   // Build up balance_box_level from old-style data.
   balance_box_level.initialize(
      hier::BoxContainer(),
      balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI(),
      hier::BoxLevel::GLOBALIZED);
   int i = 0;
   for (hier::BoxContainer::iterator itr = out_boxes.begin();
        itr != out_boxes.end(); ++itr, ++i) {
      hier::Box node(*itr, hier::LocalId(i),
                     mapping.getProcessorAssignment(i));
      balance_box_level.addBox(node);
   }
   // Reinitialize Connectors due to changed balance_box_level.
   if (balance_to_anchor) {
      hier::Connector& anchor_to_balance = balance_to_anchor->getTranspose();
      balance_to_anchor->clearNeighborhoods();
      balance_to_anchor->setBase(balance_box_level, true);
      anchor_to_balance.clearNeighborhoods();
      anchor_to_balance.setHead(balance_box_level, true);
      hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(*balance_to_anchor);
      oca.findOverlaps(anchor_to_balance, balance_box_level);
      balance_to_anchor->removePeriodicRelationships();
      anchor_to_balance.removePeriodicRelationships();
   }

   balance_box_level.setParallelState(hier::BoxLevel::DISTRIBUTED);

   t_load_balance_box_level->stop();
}

/*
 *************************************************************************
 *
 * This main load balance routine performs either uniform or
 * non-uniform load balancing operations on the given level depending
 * on user specifications.   In either case, the goal is to produce
 * a set of boxes and a mapping of those boxes to processors so that
 * the workload on each processor is close to the average workload.
 * The average workload is the total computational workload divided
 * by the number of processors.  In the uniform load balance case
 * (default), the workload is the number of cells in each box.  In the
 * non-uniform case, the workload is computed using weight data on the
 * grid hierarchy (i.e., a cell-centered double array on each patch.
 *
 * Typically, any box whose load is larger than the average is chopped.
 * A user can prescribe a parameter (the 'max workload factor') to alter
 * average load used in this computation. Chopping is done using the
 * BalanceUtilities::recursiveBisection()) method which is similar
 * to the Berger-Rigoutsos algorithm.
 *
 * Once the boxes are chopped into a collection os smaller boxes, they
 * are assigned to processors by a bin packing algorithm.
 *
 * The algorithm is summarized as follows:
 *
 * 1) Compute the estimated workload associated with each box.  In the
 *    uniform workload case, the load in each box region is the number
 *    of cells in the region.  Otherwise, the workload is computed using
 *    patch data defined by the d_workload_data_id array set by the user.
 *
 * 2) Compute the maximum workload allowed on any box.  This quantity is
 *    by default the total workload divided by the number of processors.
 *    The user may provide a maximum workload factor, either through the
 *    input file or through a member function, which can alter the
 *    average workload used in this computation.
 *
 * 3) Chop each box whose workload is more than the max allowed into a
 *    smaller set of boxes.
 *
 * 4) Check constraints placed on the boxes by the problem - i.e.
 *    verify boxes are within the maximum and minimum box size
 *    constraints and maintain a correct cut factor.
 *
 * 5) Sort boxes largest to smallest and form an array.  Also form an
 *    array of the workloads associated with each box.
 *
 * 6) Use a bin packing procedure to construct a processor mapping for
 *    the set of boxes.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::loadBalanceBoxes(
   hier::BoxContainer& out_boxes,
   hier::ProcessorMapping& mapping,
   const hier::BoxContainer& in_boxes,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::BoxContainer& physical_domain,
   const hier::IntVector& ratio_to_hierarchy_level_zero,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval) const
{
   t_load_balance_boxes->start();

   TBOX_ASSERT_OBJDIM_EQUALITY5(ratio_to_hierarchy_level_zero,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);
   TBOX_ASSERT(!physical_domain.empty());
   TBOX_ASSERT(min_size > hier::IntVector::getZero(d_dim));
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > hier::IntVector(d_dim,0));
   TBOX_ASSERT(bad_interval >= hier::IntVector::getZero(d_dim));

   /*
    * This method assumes in_boxes is not empty and will fail
    * if it is.  So shortcut it for empty in_boxes.
    */
   if (in_boxes.empty()) {
      out_boxes = hier::BoxContainer();
      return;
   }

   const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

   /*
    * If uniform load balancing is used and the level domain can be
    * expressed as a single box, we can construct an optimal box
    * layout across processors without more involved chopping operations.
    *
    * Otherwise, we chop each box individually to construct the new array
    * of boxes and associated array of workloads based on either uniform
    * or nonuniform workload estimates.
    */

   int wrk_indx = getWorkloadDataId(level_number);

   std::vector<double> workloads;

   if ((wrk_indx < 0) || (hierarchy->getNumberOfLevels() == 0)) {

      if (!d_ignore_level_box_union_is_single_box &&
          hierarchy->getGridGeometry()->getNumberBlocks() == 1) {
         hier::Box bbox = in_boxes.getBoundingBox();
         hier::BoxContainer difference(bbox);
         t_load_balance_boxes_remove_intersection->start();
         difference.removeIntersections(in_boxes);
         t_load_balance_boxes_remove_intersection->stop();

         if (difference.empty()) {

            t_chop_boxes->start();
            chopUniformSingleBox(out_boxes,
               workloads,
               bbox,
               min_size,
               max_size,
               cut_factor,
               bad_interval,
               physical_domain,
               mpi);
            t_chop_boxes->stop();

         } else {

            t_chop_boxes->start();
            chopBoxesWithUniformWorkload(out_boxes,
               workloads,
               in_boxes,
               hierarchy,
               level_number,
               min_size,
               max_size,
               cut_factor,
               bad_interval,
               physical_domain,
               mpi);
            t_chop_boxes->stop();

         }
      } else {
         t_chop_boxes->start();
         chopBoxesWithUniformWorkload(out_boxes,
            workloads,
            in_boxes,
            hierarchy,
            level_number,
            min_size,
            max_size,
            cut_factor,
            bad_interval,
            physical_domain,
            mpi);
         t_chop_boxes->stop();
      }

   } else {

      t_chop_boxes->start();
      chopBoxesWithNonuniformWorkload(out_boxes,
         workloads,
         in_boxes,
         hierarchy,
         level_number,
         ratio_to_hierarchy_level_zero,
         wrk_indx,
         min_size,
         max_size,
         cut_factor,
         bad_interval,
         physical_domain,
         mpi);
      t_chop_boxes->stop();

   }

#if 0
   /*
    * Assertion check for additional debugging - make sure new boxes
    * satisfy min size, cut_factor, bad_interval, and physical domain
    * constraints.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   const int nboxes = out_boxes.size();
   for (hier::BoxContainer::iterator ib = out_boxes.begin();
        ib != out_boxes.end(); ++ib) {
      hier::BoxUtilities::checkBoxConstraints(*ib,
         min_size,
         cut_factor,
         bad_interval,
         physical_domain);

   }
#endif
#endif

   /*
    * Generate mapping of boxes to processors using workload estimates.
    * Note boxes in array may be reordered during this operation.
    */

   binPackBoxes(out_boxes,
      mapping,
      workloads,
      getBinPackMethod(level_number));

   t_load_balance_boxes->stop();

   int my_rank = mpi.getRank();
   std::vector<int> mapping_vec = mapping.getProcessorMapping(); 
   size_t global_num_boxes = mapping_vec.size();
   double load = 0.0;
   TBOX_ASSERT(workloads.size() == global_num_boxes);
   for (size_t b = 0; b < global_num_boxes; ++b) {
      if (mapping_vec[b] == my_rank) {
         load += workloads[b];
      } else if (mapping_vec[b] > my_rank) {
         break;
      }
   }

   d_load_stat.push_back(load);

#if 0
   /*
    * For debugging, output load balance statistics
    * (assuming uniform load).
    */
   std::vector<double> procloads(tbox::SAMRAI_MPI::getNodes());
   for (int i = 0; i < static_cast<int>(procloads.size()); ++i) {
      procloads[i] = 0;
   }
   int itrCt = 0;
   for (hier::BoxContainer::iterator itr = out_boxes.begin();
        itr != out_boxes.end(); ++itr, ++itrCt) {
      int p = mapping.getProcessorAssignment(itrCt);
      procloads[p] += itr->size();
   }
   tbox::plog << "ChopAndPackLoadBalancer results (after):\n";
   reportLoadBalance(procloads);
#endif

#ifdef ChopAndPackLoadBalancer_MARKLOADFORPOSTPROCESSING
   // Performance: Output loads for global postprocessing.
   const std::vector<int>& local_indices = mapping.getLocalIndices();
   double local_load = 0;
   int local_indices_idx = 0;
   int idx = 0;
   int num_local_indices = static_cast<int>(local_indices.size());
   for (hier::BoxContainer::iterator itr = out_boxes.begin();
        itr != out_boxes.end() && local_indices_idx < num_local_indices;
        ++itr, ++idx) {
      if (local_indices[local_indices_idx] == idx) {
         local_load += static_cast<double>(itr->size());
         ++local_indices_idx;
      }
   }
   markLoadForPostprocessing(mpi.getSize(),
      local_load,
      num_local_indices);
#endif
}

/*
 *************************************************************************
 *
 * Private function that chops a single box into a set of boxes
 * that will map to the array of processors as uniformly as possible.
 * The routine assumes a spatially-uniform workload distribution.
 * Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopUniformSingleBox(
   hier::BoxContainer& out_boxes,
   std::vector<double>& out_workloads,
   const hier::Box& in_box,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{

   TBOX_ASSERT_OBJDIM_EQUALITY4(min_size, max_size, cut_factor, bad_interval);

   /*
    * Determine processor layout that corresponds to box size.
    */

   hier::IntVector processor_distribution(d_dim);
   if (d_processor_layout_specified) {
      processor_distribution = d_processor_layout;
   } else {
      BalanceUtilities::computeDomainDependentProcessorLayout(
         processor_distribution,
         mpi.getSize(),
         in_box);
   }

   /*
    * The ideal box size will be the size of the input box divided
    * by the number of processors in each direction.  Compute this
    * ideal size and then adjust as necessary to fit within min/max size
    * constraints.
    */

   hier::IntVector ideal_box_size(d_dim);
   for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
      ideal_box_size(i) = (int)ceil((double)in_box.numberCells(
               i) / (double)processor_distribution(i));

      ideal_box_size(i) = (ideal_box_size(i) > max_size(i) ?
                           max_size(i) : ideal_box_size(i));
      ideal_box_size(i) = (ideal_box_size(i) < min_size(i) ?
                           min_size(i) : ideal_box_size(i));
   }

   /*
    * Chop the single input box into a set of smaller boxes using the
    * ideal_box_size as the maximum size of each of the smaller boxes.
    */

   hier::BoxContainer tmp_box_list(in_box);

   hier::BoxUtilities::chopBoxes(tmp_box_list,
      ideal_box_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;

   const int nboxes = out_boxes.size();
   out_workloads.resize(nboxes);
   int ibCt = 0;
   for (hier::BoxContainer::iterator ib = out_boxes.begin();
        ib != out_boxes.end();
        ++ib, ++ibCt) {
      out_workloads[ibCt] = (double)(ib->size());
   }

}

/*
 *************************************************************************
 *
 * Private function that chops boxes on a list into another list of
 * boxes all of which have approximately or less than an average
 * workload estimate.  The routine assumes a spatially-uniform
 * workload distribution.   Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopBoxesWithUniformWorkload(
   hier::BoxContainer& out_boxes,
   std::vector<double>& out_workloads,
   const hier::BoxContainer& in_boxes,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{
   NULL_USE(hierarchy);
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY5(d_dim,
      *hierarchy,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxContainer tmp_in_boxes_list(in_boxes);

   /*
    * Chop any boxes in input box list that are larger than max box size
    * if possible.
    */

   hier::BoxUtilities::chopBoxes(tmp_in_boxes_list,
      max_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   double total_work = 0.0;
   for (hier::BoxContainer::iterator ib0 = tmp_in_boxes_list.begin();
        ib0 != tmp_in_boxes_list.end(); ++ib0) {
      total_work += static_cast<double>(ib0->size());
   }

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / mpi.getSize();

   hier::BoxContainer tmp_box_list;
   std::list<double> tmp_work_list;
   BalanceUtilities::recursiveBisectionUniform(tmp_box_list,
      tmp_work_list,
      tmp_in_boxes_list,
      average_work,
      getWorkloadTolerance(level_number),
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   if (tmp_box_list.size() != static_cast<int>(tmp_work_list.size())) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Number of boxes generated != number of workload values generated."
                       << std::endl);
   }

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;

   out_workloads.resize(out_boxes.size());
   int i = 0;
   for (std::list<double>::const_iterator il(tmp_work_list.begin());
        il != tmp_work_list.end(); ++il) {
      out_workloads[i] = *il;
      ++i;
   }

}

/*
 *************************************************************************
 *
 * Private function that chops boxes on a list into another list of
 * boxes all of which have approximately or less than an average
 * workload estimate.  The routine assumes a spatially-nonuniform
 * workload distribution.  Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopBoxesWithNonuniformWorkload(
   hier::BoxContainer& out_boxes,
   std::vector<double>& out_workloads,
   const hier::BoxContainer& in_boxes,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::IntVector& ratio_to_hierarchy_level_zero,
   int wrk_indx,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{

   TBOX_ASSERT_OBJDIM_EQUALITY5(ratio_to_hierarchy_level_zero,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxContainer tmp_in_boxes_list(in_boxes);

   hier::BoxUtilities::chopBoxes(tmp_in_boxes_list,
      max_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   /*
    * Create temporary patch level from in_boxes, distributed using simple
    * uniform workload estimate.  Then, fill the patch data on the level
    * corresponding to the non-uniform workload estimate.  Next, accumulate
    * the total work for the set of boxes.
    */

   hier::BoxContainer tmp_level_boxes(tmp_in_boxes_list);

   const int num_tmp_patches = tmp_level_boxes.size();
   std::vector<double> tmp_level_workloads(num_tmp_patches);
   int idx = 0;
   for (hier::BoxContainer::iterator i = tmp_level_boxes.begin();
        i != tmp_level_boxes.end();
        ++i, ++idx) {
      tmp_level_workloads[idx] = static_cast<double>(i->size());
   }

   hier::ProcessorMapping tmp_level_mapping;

   binPackBoxes(tmp_level_boxes,
      tmp_level_mapping,
      tmp_level_workloads,
      "GREEDY");

   std::shared_ptr<hier::BoxLevel> tmp_box_level(
      std::make_shared<hier::BoxLevel>(
         ratio_to_hierarchy_level_zero,
         hierarchy->getGridGeometry(),
         mpi,
         hier::BoxLevel::GLOBALIZED));
   idx = 0;
   for (hier::BoxContainer::iterator i = tmp_level_boxes.begin();
        i != tmp_level_boxes.end(); ++i, ++idx) {
      hier::Box node(*i, hier::LocalId(idx),
                     tmp_level_mapping.getProcessorAssignment(idx));
      tmp_box_level->addBox(node);
   }

   std::shared_ptr<hier::PatchLevel> tmp_level(
      std::make_shared<hier::PatchLevel>(*tmp_box_level,
                                           hierarchy->getGridGeometry(),
                                           hierarchy->getPatchDescriptor()));

   const hier::PatchLevel& hiercoarse =
      *hierarchy->getPatchLevel(level_number - 1);

   if (level_number != 0) {
      tmp_level->findConnectorWithTranspose(hiercoarse,
         hierarchy->getRequiredConnectorWidth(level_number, level_number - 1),
         hierarchy->getRequiredConnectorWidth(level_number - 1, level_number),
         hier::CONNECTOR_CREATE,
         true);
   }

   tmp_level->allocatePatchData(wrk_indx);

   xfer::RefineAlgorithm fill_work_algorithm;

   std::shared_ptr<hier::RefineOperator> work_refine_op(
      std::make_shared<pdat::CellDoubleConstantRefine>());

   fill_work_algorithm.registerRefine(wrk_indx,
      wrk_indx,
      wrk_indx,
      work_refine_op);

   std::shared_ptr<hier::PatchLevel> current_level;
   if (level_number <= hierarchy->getFinestLevelNumber()) {
      current_level = hierarchy->getPatchLevel(level_number);
   }

   fill_work_algorithm.createSchedule(tmp_level,
      current_level,
      level_number - 1,
      hierarchy)->fillData(0.0);

   double local_work = 0;
   for (hier::PatchLevel::iterator ip(tmp_level->begin());
        ip != tmp_level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      double patch_work =
         BalanceUtilities::computeNonUniformWorkload(patch,
            wrk_indx,
            patch->getBox());

      local_work += patch_work;
   }

   double total_work = local_work;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&total_work, 1, MPI_SUM);
   }

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / mpi.getSize();

   hier::BoxContainer tmp_box_list;
   std::list<double> tmp_work_list;
   BalanceUtilities::recursiveBisectionNonuniform(tmp_box_list,
      tmp_work_list,
      tmp_level,
      wrk_indx,
      average_work,
      getWorkloadTolerance(level_number),
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   tmp_level->deallocatePatchData(wrk_indx);
   tmp_level.reset();

   if (tmp_box_list.size() != static_cast<int>(tmp_work_list.size())) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Number of boxes generated != number of workload values generated."
                       << std::endl);
   }

   /*
    * Set local box array to list of chopped boxes and set local workload array.
    */

   hier::BoxContainer local_out_boxes(tmp_box_list);

   std::vector<double> local_out_workloads(local_out_boxes.size());

   int i = 0;
   for (std::list<double>::const_iterator il(tmp_work_list.begin());
        il != tmp_work_list.end(); ++il) {
      local_out_workloads[i] = *il;
      ++i;
   }

   /*
    * Gather local box and work arrays so that each processor has a copy.
    */
   exchangeBoxContainersAndWeightArrays(local_out_boxes,
      local_out_workloads,
      out_boxes,
      out_workloads,
      mpi);

}

/*
 *************************************************************************
 *
 * all-to-all exchange of box arrays and associated weights
 *
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::exchangeBoxContainersAndWeightArrays(
   const hier::BoxContainer& box_list_in,
   std::vector<double>& weights_in,
   hier::BoxContainer& box_list_out,
   std::vector<double>& weights_out,
   const tbox::SAMRAI_MPI& mpi) const
{
   TBOX_ASSERT(box_list_in.size() == static_cast<int>(weights_in.size()));

   /*
    * allocate send and receive buffers, and set array sizes
    * for the output arrays.
    */
   int size_in = box_list_in.size();
   int size_out = size_in;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&size_in, &size_out, 1, MPI_INT, MPI_SUM);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (size_out <= 0) {
      TBOX_ERROR("ChopAndPackLoadBalancer::exchangeBoxContainersAndWeightArrays() error"
         << "\n All box arrays have zero length!" << std::endl);
   }
#endif

   int buf_size_in = size_in * (d_dim.getValue() * 2 + 1);
   int buf_size_out = size_out * (d_dim.getValue() * 2 + 1);

   int curr_box_list_out_size = box_list_out.size();
   if (size_out > curr_box_list_out_size) {
      for (int i = curr_box_list_out_size; i < size_out; ++i) {
         box_list_out.emplaceBack(d_dim);
      }
   } else if (size_out < curr_box_list_out_size) {
      for (int i = size_out; i < curr_box_list_out_size; ++i) {
         box_list_out.popBack();
      }
   }
   weights_out.resize(size_out);

   std::vector<int> buf_in(buf_size_in);
   std::vector<int> buf_out(buf_size_out);

   int* buf_in_ptr = 0;
   int* buf_out_ptr = 0;
   double* wgts_in_ptr = 0;
   double* wgts_out_ptr = 0;

   if (size_in > 0) {
      buf_in_ptr = &buf_in[0];
      wgts_in_ptr = &weights_in[0];
   }
   if (size_out > 0) {
      wgts_out_ptr = &weights_out[0];
      buf_out_ptr = &buf_out[0];
   }

   /*
    * populate the buffers with data for sending
    */
   int offset = 0;
   for (hier::BoxContainer::const_iterator x = box_list_in.begin();
        x != box_list_in.end(); ++x) {
      for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
         buf_in_ptr[offset++] = x->lower(i);
         buf_in_ptr[offset++] = x->upper(i);
      }
      buf_in_ptr[offset++] = static_cast<int>(x->getBlockId().getBlockValue());
   }

   /*
    * exchange the data
    */
   std::vector<int> counts(mpi.getSize());
   mpi.Allgather(&size_in, 1, MPI_INT, &counts[0], 1, MPI_INT);
   std::vector<int> displs(mpi.getSize());
   displs[0] = 0;
   for (size_t i = 1; i < counts.size(); ++i) {
      displs[i] = displs[i - 1] + counts[i - 1];
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   size_t total_count = counts[0];
   for (size_t i = 1; i < counts.size(); ++i) {
      total_count += counts[i];
   }
   TBOX_ASSERT(weights_out.size() == total_count);
#endif

   mpi.Allgatherv(static_cast<void *>(wgts_in_ptr),
      size_in,
      MPI_DOUBLE,
      wgts_out_ptr,
      &counts[0],
      &displs[0],
      MPI_DOUBLE);

   const int ints_per_box = d_dim.getValue() * 2 + 1;
   for (size_t i = 0; i < displs.size(); ++i) {
      counts[i] *= ints_per_box;
      displs[i] *= ints_per_box;
   }
   mpi.Allgatherv(buf_in_ptr, buf_size_in, MPI_INT, buf_out_ptr, &counts[0], &displs[0], MPI_INT);

   /*
    * assemble the output array of boxes
    */
   offset = 0;
   for (hier::BoxContainer::iterator b = box_list_out.begin();
        b != box_list_out.end(); ++b) {
      for (tbox::Dimension::dir_t j = 0; j < d_dim.getValue(); ++j) {
         b->setLower(j, buf_out_ptr[offset++]);
         b->setUpper(j, buf_out_ptr[offset++]);
      }
      b->setBlockId(hier::BlockId(buf_out_ptr[offset++]));
   }

}

/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::printClassData(
   std::ostream& os) const
{
   os << "\nChopAndPackLoadBalancer::printClassData..." << std::endl;
   os << "\nChopAndPackLoadBalancer: this = "
      << (ChopAndPackLoadBalancer *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_processor_layout_specified = "
      << d_processor_layout_specified << std::endl;
   os << "d_processor_layout = "
      << d_processor_layout << std::endl;
   os << "d_master_workload_data_id = "
      << d_master_workload_data_id << std::endl;
   os << "d_master_max_workload_factor = "
      << d_master_max_workload_factor << std::endl;
   os << "d_workload_tolerance = "
      << d_master_workload_tolerance << std::endl;
   os << "d_master_bin_pack_method = "
      << d_master_bin_pack_method << std::endl;

   int ln;

   os << "d_workload_data_id..." << std::endl;
   for (ln = 0; ln < static_cast<int>(d_workload_data_id.size()); ++ln) {
      os << "    d_workload_data_id[" << ln << "] = "
         << d_workload_data_id[ln] << std::endl;
   }
   os << "d_max_workload_factor..." << std::endl;
   for (ln = 0; ln < static_cast<int>(d_max_workload_factor.size()); ++ln) {
      os << "    d_max_workload_factor[" << ln << "] = "
         << d_max_workload_factor[ln] << std::endl;
   }
   os << "d_workload_tolerance..." << std::endl;
   for (ln = 0; ln < static_cast<int>(d_workload_tolerance.size()); ++ln) {
      os << "    d_workload_tolerance[" << ln << "] = "
         << d_workload_tolerance[ln] << std::endl;
   }
   os << "d_bin_pack_method..." << std::endl;
   for (ln = 0; ln < static_cast<int>(d_bin_pack_method.size()); ++ln) {
      os << "    d_bin_pack_method[" << ln << "] = "
         << d_bin_pack_method[ln] << std::endl;
   }
   os << "d_ignore_level_box_union_is_single_box = "
      << d_ignore_level_box_union_is_single_box << std::endl;

}

/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      d_master_bin_pack_method =
         input_db->getStringWithDefault("bin_pack_method", "SPATIAL");
      if (!(d_master_bin_pack_method == "GREEDY" ||
            d_master_bin_pack_method == "SPATIAL")) {
         INPUT_VALUE_ERROR("bin_pack_method");
      }

      if (input_db->keyExists("max_workload_factor")) {
         d_max_workload_factor =
            input_db->getDoubleVector("max_workload_factor");
         for (int i = 0; i < static_cast<int>(d_max_workload_factor.size()); ++i) {
            if (!(d_max_workload_factor[i] >= 0)) {
               INPUT_RANGE_ERROR("max_workload_factor");
            }
         }

         // Use last entry in array as value for finer levels
         d_master_max_workload_factor =
            d_max_workload_factor[d_max_workload_factor.size() - 1];
      }

      if (input_db->keyExists("workload_tolerance")) {
         d_workload_tolerance =
            input_db->getDoubleVector("workload_tolerance");
         for (int i = 0; i < static_cast<int>(d_workload_tolerance.size()); ++i) {
            if (!(d_workload_tolerance[i] >= 0.0 &&
                  d_workload_tolerance[i] < 1.0)) {
               INPUT_RANGE_ERROR("workload_tolerance");
            }
         }

         // Use last entry in array as value for finer levels
         d_master_workload_tolerance =
            d_workload_tolerance[d_workload_tolerance.size() - 1];
      }

      d_ignore_level_box_union_is_single_box =
         input_db->getBoolWithDefault("ignore_level_box_union_is_single_box", false);

      d_processor_layout_specified = false;
      int temp_processor_layout[SAMRAI::MAX_DIM_VAL];
      if (input_db->keyExists("processor_layout")) {
         input_db->getIntegerArray("processor_layout",
            temp_processor_layout, d_dim.getValue());

         /* consistency check */
         int totprocs = 1;
         for (int n = 0; n < d_dim.getValue(); ++n) {
            totprocs *= temp_processor_layout[n];
         }

         if (totprocs != mpi.getSize()) {
            TBOX_WARNING(
               d_object_name << ": "
                             << "Input values for 'processor_layout' are inconsistent with"
                             << "\nnumber of processors.  Processor layout information will"
                             << "\nbe generated when needed."
                             << std::endl);
         } else {
            for (int n = 0; n < d_dim.getValue(); ++n) {
               d_processor_layout(n) = temp_processor_layout[n];
            }
            d_processor_layout_specified = true;
         }
      }

      if (input_db->isInteger("tile_size")) {
         input_db->getIntegerArray("tile_size", &d_tile_size[0], d_tile_size.getDim().getValue());
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (!(d_tile_size[i] >= 1)) {
               TBOX_ERROR("CascadePartitioner tile_size must be >= 1 in all directions.\n"
                  << "Input tile_size is " << d_tile_size);
            }
         }
      }

   }

}

/*
 *************************************************************************
 *
 * Private utility function to map boxes to processors.
 * Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::binPackBoxes(
   hier::BoxContainer& boxes,
   hier::ProcessorMapping& mapping,
   std::vector<double>& workloads,
   const std::string& bin_pack_method) const
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   t_bin_pack_boxes->start();
   /*
    * Sort boxes in order of highest to lowest workload and assign
    * boxes to processors by bin packing.
    */
   t_bin_pack_boxes_sort->start();
   BalanceUtilities::sortDescendingBoxWorkloads(boxes,
      workloads);
   t_bin_pack_boxes_sort->stop();

   /*
    * Finally, assign boxes to processors by bin packing.
    */

   int num_procs = mpi.getSize();

   t_bin_pack_boxes_pack->start();
   if (bin_pack_method == "SPATIAL") {

      BalanceUtilities::spatialBinPack(mapping,
         workloads,
         boxes,
         num_procs);

   } else if (bin_pack_method == "GREEDY") {

      BalanceUtilities::binPack(mapping,
         workloads,
         num_procs);

   } else {

      TBOX_ERROR(
         d_object_name << ": "
                       << "Unknown bin pack method "
                       << bin_pack_method << " found." << std::endl);

   }
   t_bin_pack_boxes_pack->stop();

   t_bin_pack_boxes->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::setupTimers()
{
   t_load_balance_box_level = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::loadBalanceBoxLevel()");

   t_load_balance_boxes = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::loadBalanceBoxes()");
   t_load_balance_boxes_remove_intersection = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::loadBalanceBoxes()_remove_intersection");
   t_get_global_boxes = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::get_global_boxes");
   t_bin_pack_boxes = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::binPackBoxes()");
   t_bin_pack_boxes_sort = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::binPackBoxes()_sort");
   t_bin_pack_boxes_pack = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::binPackBoxes()_pack");
   t_chop_boxes = tbox::TimerManager::getManager()->
      getTimer(d_object_name + "::chop_boxes");
}

void
ChopAndPackLoadBalancer::printStatistics(
   std::ostream& output_stream) const
{
   if (d_load_stat.empty()) {
      output_stream << "No statistics for ChopAndPackLoadBalancer.\n";
   } else {
      BalanceUtilities::reduceAndReportLoadBalance(
         d_load_stat,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         output_stream);
   }
}

}
}
