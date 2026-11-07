/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   utility routines useful for load balancing operations
 *
 ************************************************************************/

#ifndef included_mesh_BalanceUtilities
#define included_mesh_BalanceUtilities

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/mesh/PartitioningParams.h"
#include "SAMRAI/mesh/SpatialKey.h"
#include "SAMRAI/tbox/RankGroup.h"

#include <iostream>
#include <list>
#include <vector>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Utility class BalanceUtilities provides several functions
 * useful in various load balancing operations.  These utilities include
 * bin packing operations, box chopping by recursive bisection, and
 * computation of effective processor layouts for boxes.
 */

struct BalanceUtilities {
   /*!
    * Assign workloads to processors using a greedy algorithm that attempts
    * to distribute the sum of weights on each processor evenly across
    * the given number of processors.
    *
    * @return         double-valued estimate of the load balance efficiency
    *                 (ranges from zero to one hundred percent)
    *
    * @param mapping  Output processor mapping.
    * @param weights  std::vector of double-valued weights to distribute.
    * @param nproc    Integer number of processors, must be > 0.
    *
    * @pre nproc > 0
    * @pre for each memeber of weights, w, w >=0
    */
   static double
   binPack(
      hier::ProcessorMapping& mapping,
      std::vector<double>& weights,
      int nproc);

   /*!
    * Assign boxes to processors so that boxes spatially near each
    * other are likely to be assigned to processors near each other
    * (assuming that processor ordering is reflected in processor rank)
    * and so that the workload is approximately evenly distributed among
    * the processors.  The routine uses a Morton space-filling curve
    * algorithm.
    *
    * Note that this routine potentially reorders the boxes passed in to
    * achieve the first goal.
    *
    * @return         Double-valued estimate of the load balance efficiency
    *                 (ranges from zero to one hundred percent)
    *
    * @param mapping  Output processor mapping.
    * @param weights  std::vector of double-valued box weights to distribute.
    * @param boxes    hier::BoxContainer of boxes to distribute to processors.
    * @param nproc    Integer number of processors, must be > 0.
    *
    * @pre nproc > 0
    * @pre weights.size() == boxes.size()
    */
   static double
   spatialBinPack(
      hier::ProcessorMapping& mapping,
      std::vector<double>& weights,
      hier::BoxContainer& boxes,
      const int nproc);

   /*!
    * Recursively chop chops boxes in input boxlist until each box has
    * a workload less than the prescribed ideal workload or no more more
    * chopping is allowed by the given constraints.   A spatially-uniform
    * workload is assumed; i.e., all cells are weighted equally.  This routine
    * attempts to create as many boxes as possible with loads equal to or
    * slightly less than the ideal workload value so that they can be
    * mapped to processors effectively.
    *
    * @param out_boxes       Output box list.
    * @param out_workloads   Output list of box workloads.
    * @param in_boxes        Input boxlist for chopping.
    * @param ideal_workload  Input double ideal box workload, must be > 0.
    * @param workload_tolerance Input double workload tolerance, must be >= 0 and < 1.0
    * @param min_size        Input integer vector of minimum sizes for
    *                        output boxes. All entries must be > 0.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct sizes.  The box size in each
    *                        direction will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    * @param physical_domain hier::BoxContainer of boxes describing the
    *                        physical extent of the index space associated with
    *                        the in_boxes.  This box array cannot be empty.
    *
    * @pre (min_size.getDim() == cut_factor.getDim()) &&
    *      (min_size.getDim() == bad_interval.getDim())
    * @pre ideal_workload > 0
    * @pre (workload_tolerance >= 0) && (workload_tolerance < 1.0)
    * @pre min_size > hier::IntVector::getZero(min_size.getDim())
    * @pre cut_factor > hier::IntVector::getZero(min_size.getDim())
    * @pre bad_interval >= hier::IntVector::getZero(min_size.getDim())
    * @pre !physical_domain.empty()
    */
   static void
   recursiveBisectionUniform(
      hier::BoxContainer& out_boxes,
      std::list<double>& out_workloads,
      const hier::BoxContainer& in_boxes,
      double ideal_workload,
      const double workload_tolerance,
      const hier::IntVector& min_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain);

   /*!
    * Recursively chops boxes given by patches on input patch level until each
    * box has a workload less than the prescribed ideal workload or no more more
    * chopping is allowed by the given constraints.   A spatially-nonuniform
    * workload is assumed.  Cell weights must be given bydata defined by the
    * given patch data id on the given patch level.  This routine attempts to
    * create as many boxes as possible with loads equal to or slightly less
    * than the ideal workload value so that they can be mapped to processors
    * effectively.
    *
    * @param out_boxes       Output box list.
    * @param out_workloads   Output list of box workloads.
    * @param in_level        Input patch level whose patches describe input
    *                        box regions and whose patch data contain workload
    *                        estimate for each cell.
    * @param work_id         Input integer patch data id for cell-centered
    *                        double work estimate for each cell.
    * @param ideal_workload  Input double ideal box workload, must be > 0.
    * @param workload_tolerance Input double workload tolerance, must be >= 0 and < 1.0
    * @param min_size        Input integer vector of minimum sizes for
    *                        output boxes. All entries must be > 0.
    * @param cut_factor      Input integer vector used to create boxes with
    *                        correct sizes.  The box size in each
    *                        direction will be an integer multiple of the
    *                        corresponding cut factor vector entry.  All
    *                        vector entries must be > 0.  See hier::BoxUtilities
    *                        documentation for more details.
    * @param bad_interval    Input integer vector used to create boxes near
    *                        physical domain boundary with sufficient number
    *                        of cells.  No box face will be closer to the
    *                        boundary than the corresponding interval of cells
    *                        to the boundary (the corresponding value is given
    *                        by the normal direction of the box face) unless
    *                        the face coincides with the boundary itself.  The
    *                        point of this argument is to have no patch live
    *                        within a certain ghost cell width of the boundary
    *                        if its boundary does not coincide with that
    *                        boundary .  That is, all ghost cells along a face
    *                        will be either in the domain interior or outside
    *                        the domain.  All entries must be >= 0. See
    *                        hier::BoxUtilities documentation for more details.
    * @param physical_domain hier::BoxContainer of boxes describing the
    *                        physical extent of the index space associated with
    *                        the in_boxes.  This box array cannot be empty.
    *
    * @pre (min_size.getDim() == cut_factor.getDim()) &&
    *      (min_size.getDim() == bad_interval.getDim())
    * @pre ideal_workload > 0
    * @pre (workload_tolerance >= 0) && (workload_tolerance < 1.0)
    * @pre min_size > hier::IntVector::getZero(min_size.getDim())
    * @pre cut_factor > hier::IntVector::getZero(min_size.getDim())
    * @pre bad_interval >= hier::IntVector::getZero(min_size.getDim())
    * @pre !physical_domain.empty()
    */
   static void
   recursiveBisectionNonuniform(
      hier::BoxContainer& out_boxes,
      std::list<double>& out_workloads,
      const std::shared_ptr<hier::PatchLevel>& in_level,
      int work_id,
      double ideal_workload,
      const double workload_tolerance,
      const hier::IntVector& min_size,
      const hier::IntVector& cut_factor,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain);

   /*!
    * Compute factorization of processors corresponding to
    * size of given box.
    *
    * @param proc_dist  Output number of processors for each
    *                   coordinate direction.
    * @param num_procs  Input integer number of processors, must be > 0.
    * @param box        Input box to be distributed.
    *
    * @pre proc_dist.getDim() == box.getDim()
    * @pre num_procs > 0
    * @pre for each dimension, i, box.numberCells(i) > 0
    */
   static void
   computeDomainDependentProcessorLayout(
      hier::IntVector& proc_dist,
      int num_procs,
      const hier::Box& box);

   /*!
    * Compute a factorization of processors that does NOT necessarily
    * correspond to the dimensions of the supplied box.  For example, the
    * processor distribution in each direction may simply be a square root
    * (cube root in 3D) of the number of processors.  The box information
    * is used simply to determine a maximum number of processors in each
    * coordinate direction.
    *
    * @param proc_dist  Output number of processors for each
    *                   coordinate direction.
    * @param num_procs  Input integer number of processors, must be > 0.
    * @param box        Input box to be distributed.
    *
    * @pre proc_dist.getDim() == box.getDim()
    * @pre num_procs > 0
    * @pre for each dimension, i, box.numberCells(i) > 0
    */
   static void
   computeDomainIndependentProcessorLayout(
      hier::IntVector& proc_dist,
      int num_procs,
      const hier::Box& box);

   /*!
    * Sort box array in descending order of workload according to the
    * workload array.  Both the box array and the work array will be
    * sorted on return.
    *
    * Note that if you simply want to sort boxes based on their size,
    * see hier::BoxUtilities.
    *
    * @param boxes     Boxes to be sorted based on workload array.
    * @param workload  Workloads to use for sorting boxes.
    *
    * @pre boxes.size() == workload.size()
    */
   static void
   sortDescendingBoxWorkloads(
      hier::BoxContainer& boxes,
      std::vector<double>& workload);

   /*!
    * Compute total workload in region of argument box based on patch
    * data defined by given integer index.  The sum is computed on the
    * intersection of argument box and box over which data associated with
    * workload is defined.
    *
    * @return          Sum of workload values in box region.
    *
    * @param patch     Input patch on which workload data is defined.
    * @param wrk_indx  Input integer patch data identifier for work data.
    * @param box       Input box region
    *
    * Note that wrk_indx must refer to a valid cell-centered patch data entry.
    *
    * @pre patch
    * @pre patch->getDim() == box.getDim()
    */
   static double
   computeNonUniformWorkload(
      const std::shared_ptr<hier::Patch>& patch,
      int wrk_indx,
      const hier::Box& box);

   /*!
    * @brief Compute total workload in region of argument box based on patch
    * data defined by given integer index, while also associating weights with
    * corners of the box.
    *
    * The corner weights are fractional values that identify how much of the
    * workload is located in each quadrant or octant of the box.
    *
    * @return          Sum of workload values in the box.
    *
    * @param[out] corner_weights  Fraction of weights associated with
    *                             each corner of the patch. Vector must be
    *                             empty on input. 
    * @param[in]  patch     Patch on which workload data is defined.
    * @param[in]  wrk_indx  Patch data identifier for work data.
    * @param[in]  box       Box on which workload is computed
    *
    * @pre box.getBlockId() == patch->getBox().getBlockId()
    * @pre box.isSpatiallyEqual(box * patch->getBox())
    */
   static double
   computeNonUniformWorkloadOnCorners(
      std::vector<double>& corner_weights,
      const std::shared_ptr<hier::Patch>& patch,
      int wrk_indx,
      const hier::Box& box);

   /*!
    * @brief Find small boxes in a post-balance BoxLevel that are not
    * in a pre-balance BoxLevel.
    *
    * @param co Stream to report findings
    *
    * @param border Left border in report output
    *
    * @param [in] post_to_pre
    *
    * @param [in] min_width Report post-balance boxes smaller than
    * min_width in any direction.
    *
    * @param [in] min_vol Report post-balance boxes with fewer cells
    * than this.
    */
   static void
   findSmallBoxesInPostbalance(
      std::ostream& co,
      const std::string& border,
      const hier::MappingConnector& post_to_pre,
      const hier::IntVector& min_width,
      size_t min_vol);

   /*!
    * @brief Find small boxes in a post-balance BoxLevel that are not
    * in a pre-balance BoxLevel.
    *
    * This method does not scale.  It acquires and processes
    * globalized data.
    *
    * @param co Stream to report findings
    *
    * @param border Left border in report output
    *
    * @param [in] pre Pre-balance BoxLevel
    *
    * @param [in] post Post-balance BoxLevel
    *
    * @param [in] min_width Report post-balance boxes smaller than
    * min_width in any direction.
    *
    * @param [in] min_vol Report post-balance boxes with fewer cells
    * than this.
    */
   static void
   findSmallBoxesInPostbalance(
      std::ostream& co,
      const std::string& border,
      const hier::BoxLevel& pre,
      const hier::BoxLevel& post,
      const hier::IntVector& min_width,
      size_t min_vol);

   /*!
    * @brief Evaluate whether a new load is an improvement over a
    * current load based on their proximity to an ideal value or range
    * of acceptable values.
    *
    * There is a slight bias toward current load.  The new_load is better
    * only if it improves by at least pparams.getLoadComparisonTol().
    *
    * Return values in flags:
    * - [0]: -1, 0 or 1: degrades, leave-alone or improves in-range
    * - [1]: -1, 0 or 1: degrades, leave-alone or improves balance
    * - [2]: -1, 0 or 1: degrades, leave-alone or improves overall
    * - [3]: 0 or 1: whether new_load is within the range of [low, high]
    *
    * Return whether new_load is an improvement over current_load.
    */
   static bool
   compareLoads(
      int flags[],
      double current_load,
      double new_load,
      double ideal_load,
      double low_load,
      double high_load,
      const PartitioningParams &pparams);

   /*!
    * Compute and return load balance efficiency for a level.
    *
    * @return         Double-valued estimate of the load balance efficiency
    *                 (ranges from zero to one hundred percent)
    *
    * @param level            Input patch level to consider, can't be null.
    * @param os               Output stream for reporting load balance
    *                         details.
    * @param workload_data_id (Optional) Input integer id for workload
    *                         data on level.  If no value is given, the
    *                         calculation assumes spatially-uniform load.
    *
    * @pre level
    */
   static double
   computeLoadBalanceEfficiency(
      const std::shared_ptr<hier::PatchLevel>& level,
      std::ostream& os,
      int workload_data_id = -1);

   //@{

   //! @name Load balance reporting.

   /*!
    * @brief Globally reduce a sequence of workloads in an MPI group
    * and write out a summary of load balance efficiency.
    *
    * Each value in the sequence of workloads represent a certain load
    * the local process had over a sequence of load balancings.
    *
    * To be used for performance evaluation.  Not recommended for
    * general use.
    *
    * @param[in] local_loads Sequence of workloads of the local
    * process.  The size of @c local_loads is the number times load
    * balancing has been used.  It must be the same across all
    * processors in @c mpi.
    *
    * @param[in] mpi Represents all processes involved in the load balancing.
    *
    * @param[in] output_stream
    */
   static void
   reduceAndReportLoadBalance(
      const std::vector<double>& local_loads,
      const tbox::SAMRAI_MPI& mpi,
      std::ostream& output_stream = tbox::plog);

   //@}

   /*
    * Constrain maximum box sizes in the given BoxLevel and
    * update given Connectors to the changed BoxLevel.
    *
    * @pre !anchor_to_level || anchor_to_level->hasTranspose()
    */
   static void
   constrainMaxBoxSizes(
      hier::BoxLevel& box_level,
      hier::Connector* anchor_to_level,
      const PartitioningParams& pparams);

   static const int BalanceUtilities_PREBALANCE0 = 5;
   static const int BalanceUtilities_PREBALANCE1 = 6;

   /*!
    * Move Boxes in balance_box_level from ranks outside of
    * rank_group to ranks inside rank_group.  Modify the given connectors
    * to make them correct following this moving of boxes.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre !balance_to_anchor || (balance_to_anchor->getTranspose().checkTransposeCorrectness(*balance_to_anchor) == 0)
    * @pre !balance_to_anchor || (balance_to_anchor->checkTransposeCorrectness(balance_to_anchor->getTranspose()) == 0)
    */
   static void
   prebalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const tbox::RankGroup& rank_group);

private:
   struct RankAndLoad {
      int rank;
      double load;
   };

   static int
   qsortRankAndLoadCompareAscending(
      const void* v,
      const void* w);

   static int
   qsortRankAndLoadCompareDescending(
      const void* v,
      const void* w);

   static math::PatchCellDataNormOpsReal<double> s_norm_ops;

   static void
   privateHeapify(
      std::vector<int>& permutation,
      std::vector<double>& workload,
      const int index,
      const int heap_size);

   static void
   privateHeapify(
      std::vector<int>& permutation,
      std::vector<SpatialKey>& spatial_keys,
      const int index,
      const int heap_size);

   static void
   privateRecursiveProcAssign(
      const int wt_index_lo,
      const int wt_index_hi,
      std::vector<double>& weights,
      const int proc_index_lo,
      const int proc_index_hi,
      hier::ProcessorMapping& mapping,
      const double avg_weight);

   static void
   privatePrimeFactorization(
      const int N,
      std::vector<int>& p);

   static void
   privateResetPrimesArray(
      std::vector<int>& p);

   static bool
   privateBadCutPointsExist(
      const hier::BoxContainer& physical_domain);

   static void
   privateInitializeBadCutPointsForBox(
      std::vector<std::vector<bool> >& bad_cut_points,
      hier::Box& box,
      bool bad_domain_boundaries_exist,
      const hier::IntVector& bad_interval,
      const hier::BoxContainer& physical_domain);

   static bool
   privateFindBestCutDimension(
      tbox::Dimension::dir_t& cut_dim_out,
      const hier::Box& in_box,
      const hier::IntVector& min_size,
      const hier::IntVector& cut_factor,
      std::vector<std::vector<bool> >& bad_cut_points);

   static int
   privateFindCutPoint(
      double total_work,
      double ideal_workload,
      int mincut,
      int numcells,
      const std::vector<double>& work_in_slice,
      const std::vector<bool>& bad_cut_points);

   static void
   privateCutBoxesAndSetBadCutPoints(
      hier::Box& box_lo,
      std::vector<std::vector<bool> >& bad_cut_points_for_boxlo,
      hier::Box& box_hi,
      std::vector<std::vector<bool> >& bad_cut_points_for_boxhi,
      const hier::Box& in_box,
      tbox::Dimension::dir_t cutdim,
      int cut_index,
      const std::vector<std::vector<bool> >& bad_cut_points);

   static void
   privateRecursiveBisectionUniformSingleBox(
      hier::BoxContainer& out_boxes,
      std::list<double>& out_workloads,
      const hier::Box& in_box,
      double in_box_workload,
      double ideal_workload,
      const double workload_tolerance,
      const hier::IntVector& min_size,
      const hier::IntVector& cut_factor,
      std::vector<std::vector<bool> >& bad_cut_points);

   static void
   privateRecursiveBisectionNonuniformSingleBox(
      hier::BoxContainer& out_boxes,
      std::list<double>& out_workloads,
      const std::shared_ptr<hier::Patch>& patch,
      const hier::Box& in_box,
      double in_box_workload,
      int work_data_index,
      double ideal_workload,
      const double workload_tolerance,
      const hier::IntVector& min_size,
      const hier::IntVector& cut_factor,
      std::vector<std::vector<bool> >& bad_cut_points);

};

}
}
#endif
