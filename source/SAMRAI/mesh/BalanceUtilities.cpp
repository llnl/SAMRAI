/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   utility routines useful for load balancing operations
 *
 ************************************************************************/

#ifndef included_mesh_BalanceUtilities_C
#define included_mesh_BalanceUtilities_C

#include "SAMRAI/mesh/BalanceUtilities.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

math::PatchCellDataNormOpsReal<double> BalanceUtilities::s_norm_ops;

const int BalanceUtilities::BalanceUtilities_PREBALANCE0;
const int BalanceUtilities::BalanceUtilities_PREBALANCE1;

/*
 *************************************************************************
 *
 * Two internal functions for heap sorting an array of doubles
 * and an array of spatial key information.  See any standard sorting
 * references for more information about heapsort.
 *
 *************************************************************************
 */

void
BalanceUtilities::privateHeapify(
   std::vector<int>& permutation,
   std::vector<double>& workload,
   const int index,
   const int heap_size)
{
   const int l = 2 * index + 1;
   const int r = l + 1;
   int s = index;
   if ((l < heap_size) &&
       (workload[permutation[s]] > workload[permutation[l]])) {
      s = l;
   }
   if ((r < heap_size) &&
       (workload[permutation[s]] > workload[permutation[r]])) {
      s = r;
   }
   if (s != index) {
      const int tmp = permutation[s];
      permutation[s] = permutation[index];
      permutation[index] = tmp;
      privateHeapify(permutation, workload, s, heap_size);
   }
}

void
BalanceUtilities::privateHeapify(
   std::vector<int>& permutation,
   std::vector<SpatialKey>& spatial_keys,
   const int index,
   const int heap_size)
{
   const int l = 2 * index + 1;
   const int r = l + 1;
   int s = index;
   if ((l < heap_size) &&
       (spatial_keys[permutation[s]] < spatial_keys[permutation[l]])) {
      s = l;
   }
   if ((r < heap_size) &&
       (spatial_keys[permutation[s]] < spatial_keys[permutation[r]])) {
      s = r;
   }
   if (s != index) {
      const int tmp = permutation[s];
      permutation[s] = permutation[index];
      permutation[index] = tmp;
      privateHeapify(permutation, spatial_keys, s, heap_size);
   }
}

/*
 *************************************************************************
 *
 * Internal functions that recursively assigns weights to
 * processors (used in the spatial bin packing procedure).
 *
 *************************************************************************
 */
void
BalanceUtilities::privateRecursiveProcAssign(
   const int wt_index_lo,
   const int wt_index_hi,
   std::vector<double>& weights,
   const int proc_index_lo,
   const int proc_index_hi,
   hier::ProcessorMapping& mapping,
   const double avg_weight)
{

   TBOX_ASSERT(wt_index_hi >= wt_index_lo);
   TBOX_ASSERT(proc_index_hi >= proc_index_lo);
   TBOX_ASSERT((wt_index_hi - wt_index_lo) >= (proc_index_hi - proc_index_lo));

   int i;

   /*
    * if there is only one processor in range, then assign all boxes
    * in the weight index range to the processor
    */
   if (proc_index_hi == proc_index_lo) {
      for (i = wt_index_lo; i <= wt_index_hi; ++i) {
         mapping.setProcessorAssignment(i, proc_index_lo);
      }
   } else {  // otherwise recurse

      int nproc = proc_index_hi - proc_index_lo + 1;
      int left = nproc / 2;

      double cut_weight = left * avg_weight;

      int cut_index = wt_index_lo;
      double acc_weight = 0.0;

      /*
       * the loop ends with (cut_index-1) referring to the last
       * index that contributed to acc_weight (this is the convention
       * that is used for box cut points throughout the library)
       */
      while ((cut_index <= wt_index_hi) && (acc_weight < cut_weight)) {
         acc_weight += weights[cut_index];
         ++cut_index;
      }

      /*
       * move cut_index back if cumulative weight without last weight
       * added is closer to cut_weight, but only if cut_index is strictly
       * greater than (wt_index_lo + 1) so that we don't accidentally
       * get a cut_index of wt_index_lo which makes no sense.
       */
      double prev_weight = acc_weight - weights[cut_index - 1];
      if ((cut_index > wt_index_lo + 1) &&
          ((acc_weight - cut_weight) > (cut_weight - prev_weight))) {
         --cut_index;
      }

      /*
       * shift processors around to make sure that there are more procs
       * than weights for each of the smaller pieces
       */
      if (cut_index - wt_index_lo < left) {
         cut_index = left + wt_index_lo;
      } else if (wt_index_hi - cut_index + 1 < nproc - left) {
         cut_index = wt_index_hi + 1 - nproc + left;
      }

      /*
       * recurse on smaller segments of the processor mapping and weights
       * array.
       */
      privateRecursiveProcAssign(wt_index_lo, (cut_index - 1), weights,
         proc_index_lo, (proc_index_lo + left - 1), mapping, avg_weight);
      privateRecursiveProcAssign(cut_index, wt_index_hi, weights,
         (proc_index_lo + left), proc_index_hi, mapping, avg_weight);
   }

}

/*
 *************************************************************************
 *
 * Two internal functions for computing arrays of prime numbers.
 *
 * The first computes a prime factorization of N and stores the primes
 * in the array p.  The factorization algorithm uses trial division
 * described by Donald E. Knuth, The Art of Computer Programming,
 * 3rd edition, volume 2 Seminumerical Algorithms (Addison-Wesley,
 * 1998), section 4.5.4 (Factoring Into Primes), pp. 380-381.
 *
 * The second resets an array of primes by removing all instances of "1".
 *
 *************************************************************************
 */

void
BalanceUtilities::privatePrimeFactorization(
   const int N,
   std::vector<int>& p)
{
   /*
    * Input: N
    * Output: p[] - array of prime factors of N
    */

   // Step 1 - Initialization

   int k = 0;
   int t = 0;
   int n = N;
   int q;
   int r;
   if (p.size() < 1) p.resize(1);
   p[0] = 1;

   //  NOTE: d must hold the list of prime numbers up to sqrt(n).  We
   //  store prime numbers up to 101 which assures this method will work
   //  for up to 10201 processors.
   int d[] =
   { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
     67,
     71, 73, 79, 83, 89, 97, 101 };

   // Step 2 - once n is 1, terminate the algorithm

   while (n > 1) {

      // Step 3 - divide to form quotient, remainder:
      //          q = n/d[k],  r = n mod d[k]

      q = n / d[k];
      r = n % d[k];

      // Step 4 - zero remainder?

      if (r == 0) {

         // Step 5 - factor found. Increase t by 1, set p[t] = d[k], n = q.

         ++t;
         p.resize(t + 1);
         p[t] = d[k];
         n = q;

      } else {

         // Step 6 - low quotient? Increment k to try next prime.

         if (q > d[k]) {

            ++k;

         } else {

            // Step 7 - n is prime.  Increment t by 1, set p[t] = n, and terminate.

            ++t;
            p.resize(t + 1);
            p[t] = n;
            break;
         }

      }

   }

}

void
BalanceUtilities::privateResetPrimesArray(
   std::vector<int>& p)
{
   // keep a copy of the original p in array "temp"
   std::vector<int> temp;
   temp.resize(static_cast<int>(p.size()));
   int i;
   for (i = 0; i < static_cast<int>(p.size()); ++i) temp[i] = p[i];

   // resize p to only keep values > 1
   int newsize = 0;
   for (i = 0; i < static_cast<int>(p.size()); ++i) {
      if (p[i] > 1) ++newsize;
   }

   p.resize(newsize);
   newsize = 0;

   // set values in the new p array
   for (i = 0; i < static_cast<int>(temp.size()); ++i) {
      if (temp[i] > 1) {
         p[newsize] = temp[i];
         ++newsize;
      }
   }
}

/*
 *************************************************************************
 *
 * Internal function to determine whether bad cut points exist for
 * domain.  Note that no error checking is done.
 *
 *************************************************************************
 */

bool
BalanceUtilities::privateBadCutPointsExist(
   const hier::BoxContainer& physical_domain)
{
   bool bad_cuts_exist = false;

   std::map<hier::BlockId, hier::BoxContainer> domain_by_blocks;
   for (hier::BoxContainer::const_iterator itr = physical_domain.begin();
        itr != physical_domain.end(); ++itr) {
      const hier::BlockId& block_id = itr->getBlockId();
      domain_by_blocks[block_id].pushBack(*itr);
   }
   for (std::map<hier::BlockId, hier::BoxContainer>::iterator m_itr =
           domain_by_blocks.begin(); m_itr != domain_by_blocks.end(); ++m_itr) {
      hier::BoxContainer bounding_box(m_itr->second.getBoundingBox());
      bounding_box.removeIntersections(m_itr->second);
      if (!bounding_box.empty()) {
         bad_cuts_exist = true;
      }
   }

   return bad_cuts_exist;
}

/*
 *************************************************************************
 *
 * Internal function to initialize bad cut points along each coordinate
 * direction of given box based on relation to domain boundary.
 * Note that no error checking is done.
 *
 *************************************************************************
 */

void
BalanceUtilities::privateInitializeBadCutPointsForBox(
   std::vector<std::vector<bool> >& bad_cut_points,
   hier::Box& box,
   bool bad_domain_boundaries_exist,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(box, bad_interval);

   const tbox::Dimension dim(box.getDim());

   tbox::Dimension::dir_t ic, id;

   bool set_dummy_cut_points = true;

   if (bad_domain_boundaries_exist) {

      hier::IntVector tmp_max_gcw =
         hier::VariableDatabase::getDatabase()->
         getPatchDescriptor()->getMaxGhostWidth(dim);

      hier::BoxContainer bdry_list(box);
      bdry_list.grow(tmp_max_gcw);
      bdry_list.removeIntersections(physical_domain);
      if (!bdry_list.empty()) {
         set_dummy_cut_points = false;
      }

   }

   if (set_dummy_cut_points) {

      for (id = 0; id < dim.getValue(); ++id) {
         const int ncells = box.numberCells(id);
         bad_cut_points[id].resize(ncells);
         std::vector<bool>& arr_ref = bad_cut_points[id];
         for (ic = 0; ic < ncells; ++ic) {
            arr_ref[ic] = false;
         }
      }

   } else {

      for (id = 0; id < dim.getValue(); ++id) {
         hier::BoxUtilities::
         findBadCutPointsForDirection(id,
            bad_cut_points[id],
            box,
            physical_domain,
            bad_interval);
      }

   }

}

/*
 *************************************************************************
 *
 * Internal function to determine best cut direction for a box based
 * on constraints and adjust bad cut points as needed.  Return
 * value is true if some direction can be cut; false, otherwise.
 * If the box can be cut along some direction, then cut_dim_out is
 * set to the longest box direction that can be cut; otherwise,
 * cut_dim_out is set to the invalid value of SAMRAI::MAX_DIM_VAL.
 * Note no error checking is done.
 *
 *************************************************************************
 */

bool
BalanceUtilities::privateFindBestCutDimension(
   tbox::Dimension::dir_t& cut_dim_out,
   const hier::Box& in_box,
   const hier::IntVector& min_size,
   const hier::IntVector& cut_factor,
   std::vector<std::vector<bool> >& bad_cut_points)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(in_box, min_size, cut_factor);

   const tbox::Dimension& dim(in_box.getDim());

   bool can_cut_box = false;
   cut_dim_out = SAMRAI::MAX_DIM_VAL;

   hier::Box size_test_box(in_box);

   for (tbox::Dimension::dir_t id = 0; id < dim.getValue(); ++id) {
      int ncells = in_box.numberCells(id);
      if ((ncells < 2 * min_size(id)) ||
          (ncells % cut_factor(id))) {
         size_test_box.setLower(id, size_test_box.upper(id));
      }
   }

   if (size_test_box.size() > 1) {

      /*
       * Find good cut points along some box direction, starting with longest
       * direction, then trying next longest, etc., until good cut points found.
       */

      hier::Box test_box(size_test_box);
      tbox::Dimension::dir_t cutdim = test_box.longestDirection();
      int numcells = test_box.numberCells(cutdim);
      int cutfact = cut_factor(cutdim);
      int mincut = tbox::MathUtilities<int>::Max(min_size(cutdim), cutfact);

      int i;
      bool found_cut_point = false;
      while (!found_cut_point && (test_box.numberCells(cutdim) > 1)) {

         /*
          * Make potential cut points bad if they are sufficiently near the
          * box boundary or if they are are not divisible by cut factor.
          * Then, determine whether any good cut points exist along chosen
          * coordinate direction.
          */

         std::vector<bool>& bad_cuts_for_dir = bad_cut_points[cutdim];

         for (i = 0; i < mincut; ++i) {
            bad_cuts_for_dir[i] = true;
         }
         for (i = (numcells - mincut + 1); i < numcells; ++i) {
            bad_cuts_for_dir[i] = true;
         }
         for (i = 0; i < numcells; ++i) {
            if (i % cutfact) {
               bad_cuts_for_dir[i] = true;
            }
            found_cut_point = found_cut_point || !bad_cuts_for_dir[i];
         }

         if (!found_cut_point) {
            test_box.setLower(cutdim, test_box.upper(cutdim));
         }

         cutdim = test_box.longestDirection();
         numcells = test_box.numberCells(cutdim);
         cutfact = cut_factor(cutdim);
         mincut = tbox::MathUtilities<int>::Max(min_size(cutdim), cutfact);

      }

      if (found_cut_point) {
         can_cut_box = true;
         cut_dim_out = cutdim;
      }

   }

   return can_cut_box;

}

/*
 *************************************************************************
 *
 * Internal function to determine cut point for a single direction
 * given min cut, ideal workload, bad cut point constraints.
 * Note no error checking is done.
 *
 *************************************************************************
 */

int
BalanceUtilities::privateFindCutPoint(
   double total_work,
   double ideal_workload,
   int mincut,
   int numcells,
   const std::vector<double>& work_in_slice,
   const std::vector<bool>& bad_cut_points)
{

   int cut_index = 0;

   int half_num_pieces = (int)(total_work / ideal_workload + 1) / 2;
   double work_cutpt = half_num_pieces * ideal_workload;

   /*
    * Search for cutpoint closest to "work cutpoint"; i.e.,
    * point where work to "left" is closest to ideal workload.
    */

   double acc_work = 0.0;

   while (cut_index < mincut) {
      acc_work += work_in_slice[cut_index];
      ++cut_index;
   }

   int last = numcells - mincut;
   while ((acc_work < work_cutpt) && (cut_index < last)) {
      acc_work += work_in_slice[cut_index];
      ++cut_index;
   }

   /*
    * If estimated cutpoint is bad, search to left and right for a
    * good cut point.  Choose the one that will split box closer to the
    * work cut point.
    */

   if (bad_cut_points[cut_index]) {

      int l_index = cut_index;
      double l_work = acc_work;
      while ((bad_cut_points[l_index]) && (l_index > 2)) {
         l_work -= work_in_slice[l_index - 1];
         --l_index;
      }

      int r_index = cut_index;
      double r_work = acc_work;
      while ((bad_cut_points[r_index]) && (r_index < numcells - 1)) {
         r_work += work_in_slice[r_index];
         ++r_index;
      }

      if ((work_cutpt - l_work) < (r_work - work_cutpt)) {
         if (bad_cut_points[l_index]) {
            cut_index = r_index;
         } else {
            cut_index = l_index;
         }
      } else {
         if (bad_cut_points[r_index]) {
            cut_index = l_index;
         } else {
            cut_index = r_index;
         }
      }

   }

   return cut_index;

}

/*
 *************************************************************************
 *
 * Internal function to cut box in two at given cut point along given
 * direction.  box_lo, box_hi will be new disjoint boxes whose union
 * is the box to be cut (in_box).  bad_cut_points_for_boxlo, and
 * bad_cut_points_for_boxhi are associated arrays of bad cut points
 * defined by given bad cut point arrays for in_box.
 * Note no error checking is done.
 *
 *************************************************************************
 */

void
BalanceUtilities::privateCutBoxesAndSetBadCutPoints(
   hier::Box& box_lo,
   std::vector<std::vector<bool> >& bad_cut_points_for_boxlo,
   hier::Box& box_hi,
   std::vector<std::vector<bool> >& bad_cut_points_for_boxhi,
   const hier::Box& in_box,
   tbox::Dimension::dir_t cutdim,
   int cut_index,
   const std::vector<std::vector<bool> >& bad_cut_points)
{

   TBOX_ASSERT_OBJDIM_EQUALITY3(box_lo, box_hi, in_box);

   const tbox::Dimension& dim(box_lo.getDim());

   box_lo = in_box;
   box_lo.setUpper(cutdim, cut_index - 1);

   box_hi = in_box;
   box_hi.setLower(cutdim, cut_index);

   int i;
   for (tbox::Dimension::dir_t id = 0; id < dim.getValue(); ++id) {

      const std::vector<bool>& arr_ref_in = bad_cut_points[id];

      const int ncellslo = box_lo.numberCells(id);
      const int ncellshi = box_hi.numberCells(id);

      bad_cut_points_for_boxlo[id].resize(ncellslo);
      bad_cut_points_for_boxhi[id].resize(ncellshi);

      std::vector<bool>& arr_ref_cutlo = bad_cut_points_for_boxlo[id];
      for (i = 0; i < ncellslo; ++i) {
         arr_ref_cutlo[i] = arr_ref_in[i];
      }

      std::vector<bool>& arr_ref_cuthi = bad_cut_points_for_boxhi[id];

      if (id == cutdim) {
         int mark = box_lo.numberCells(cutdim);
         for (i = 0; i < ncellshi; ++i) {
            arr_ref_cuthi[i] = arr_ref_in[i + mark];
         }
      } else {
         for (i = 0; i < ncellshi; ++i) {
            arr_ref_cuthi[i] = arr_ref_in[i];
         }
      }

   }

}

/*
 *************************************************************************
 *
 * Internal recursive function to chop a single box into two
 * boxes, if possible, based on uniform workload estimates.  It is
 * assumed that the bad-cut-point arrays are set properly according
 * to the physical domain constraints.  Note no error checking is done.
 *
 *************************************************************************
 */

void
BalanceUtilities::privateRecursiveBisectionUniformSingleBox(
   hier::BoxContainer& out_boxes,
   std::list<double>& out_workloads,
   const hier::Box& in_box,
   double in_box_workload,
   double ideal_workload,
   const double workload_tolerance,
   const hier::IntVector& min_size,
   const hier::IntVector& cut_factor,
   std::vector<std::vector<bool> >& bad_cut_points)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(in_box, min_size, cut_factor);

   const tbox::Dimension dim(in_box.getDim());

   if (in_box_workload <= ((1. + workload_tolerance) * ideal_workload)) {

      out_boxes.pushFront(in_box);
      out_workloads.push_front(in_box_workload);

   } else {

      /*
       * Determine best direction to chop box.
       */
      tbox::Dimension::dir_t cut_dim;
      bool can_cut_box = privateFindBestCutDimension(
            cut_dim,
            in_box,
            min_size,
            cut_factor,
            bad_cut_points);

      if (can_cut_box) {

         int i;

         const int numcells = in_box.numberCells(cut_dim);
         int mincut =
            tbox::MathUtilities<int>::Max(min_size(cut_dim), cut_factor(cut_dim));

         /*
          * Search for chop point along chosen direction.
          */

         double work_in_single_slice = 1.0;
         for (tbox::Dimension::dir_t id = 0; id < dim.getValue(); ++id) {
            if (id != cut_dim) {
               work_in_single_slice *= (double)in_box.numberCells(id);
            }
         }

         std::vector<double> work_in_slices(numcells);
         for (i = 0; i < numcells; ++i) {
            work_in_slices[i] = work_in_single_slice;
         }

         int cut_index = privateFindCutPoint(in_box_workload,
               ideal_workload,
               mincut,
               numcells,
               work_in_slices,
               bad_cut_points[cut_dim]);
         cut_index += in_box.lower(cut_dim);

         /*
          * Create two new boxes based on cut index and define new
          * bad cut point arrays.  Then apply recursive bisection
          * to each new box.
          */

         hier::Box box_lo(dim);
         hier::Box box_hi(dim);

         std::vector<std::vector<bool> > bad_cut_points_for_boxlo(dim.getValue());
         std::vector<std::vector<bool> > bad_cut_points_for_boxhi(dim.getValue());

         privateCutBoxesAndSetBadCutPoints(box_lo,
            bad_cut_points_for_boxlo,
            box_hi,
            bad_cut_points_for_boxhi,
            in_box,
            cut_dim,
            cut_index,
            bad_cut_points);

         double box_lo_workload = (double)box_lo.size();
         privateRecursiveBisectionUniformSingleBox(out_boxes,
            out_workloads,
            box_lo,
            box_lo_workload,
            ideal_workload,
            workload_tolerance,
            min_size,
            cut_factor,
            bad_cut_points_for_boxlo);

         hier::BoxContainer boxes_hi;
         std::list<double> work_hi;

         double box_hi_workload = (double)box_hi.size();
         privateRecursiveBisectionUniformSingleBox(boxes_hi,
            work_hi,
            box_hi,
            box_hi_workload,
            ideal_workload,
            workload_tolerance,
            min_size,
            cut_factor,
            bad_cut_points_for_boxhi);

         out_boxes.spliceBack(boxes_hi);
         out_workloads.splice(out_workloads.end(), work_hi);

      } else {  // !can_cut_box

         out_boxes.pushFront(in_box);
         out_workloads.push_front(in_box_workload);

      }

   } // in_box_workload > ideal_workload

}

/*
 *************************************************************************
 *
 * Internal recursive function to chop a single box into two
 * boxes, if possible, based on nonuniform workload estimates.  It is
 * assumed that the bad-cut-point arrays are set properly according
 * to the physical domain constraints.  Note no error checking is done.
 *
 *************************************************************************
 */

void
BalanceUtilities::privateRecursiveBisectionNonuniformSingleBox(
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
   std::vector<std::vector<bool> >& bad_cut_points)
{
   TBOX_ASSERT(patch);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*patch, in_box, min_size, cut_factor);

   const tbox::Dimension dim(in_box.getDim());

   if (in_box_workload <= ((1. + workload_tolerance) * ideal_workload)) {

      out_boxes.pushFront(in_box);
      out_workloads.push_front(in_box_workload);

   } else {

      /*
       * Determine best direction to chop box.
       */
      tbox::Dimension::dir_t cut_dim;
      bool can_cut_box = privateFindBestCutDimension(
            cut_dim,
            in_box,
            min_size,
            cut_factor,
            bad_cut_points);

      if (can_cut_box) {

         int i;

         const int numcells = in_box.numberCells(cut_dim);
         int mincut =
            tbox::MathUtilities<int>::Max(min_size(cut_dim), cut_factor(cut_dim));

         /*
          * Search for chop point along chosen direction.
          */

         hier::Box slice_box = in_box;
         slice_box.setUpper(cut_dim, slice_box.lower(cut_dim));

         std::vector<double> work_in_slices(numcells);
         for (i = 0; i < numcells; ++i) {
            work_in_slices[i] =
               BalanceUtilities::computeNonUniformWorkload(patch,
                  work_data_index,
                  slice_box);
            slice_box.setLower(cut_dim, slice_box.lower(cut_dim) + 1);
            slice_box.setUpper(cut_dim, slice_box.lower(cut_dim));

         }

         int cut_index = privateFindCutPoint(in_box_workload,
               ideal_workload,
               mincut,
               numcells,
               work_in_slices,
               bad_cut_points[cut_dim]);
         cut_index += in_box.lower(cut_dim);

         /*
          * Create two new boxes based on cut index and define new
          * bad cut point arrays.  Then apply recursive bisection
          * to each new box.
          */

         hier::Box box_lo(dim);
         hier::Box box_hi(dim);

         std::vector<std::vector<bool> > bad_cut_points_for_boxlo(dim.getValue());
         std::vector<std::vector<bool> > bad_cut_points_for_boxhi(dim.getValue());

         privateCutBoxesAndSetBadCutPoints(box_lo,
            bad_cut_points_for_boxlo,
            box_hi,
            bad_cut_points_for_boxhi,
            in_box,
            cut_dim,
            cut_index,
            bad_cut_points);

         const int box_lo_ncells = box_lo.numberCells(cut_dim);
         double box_lo_workload = 0.0;
         for (i = 0; i < box_lo_ncells; ++i) {
            box_lo_workload += work_in_slices[i];
         }
         privateRecursiveBisectionNonuniformSingleBox(out_boxes,
            out_workloads,
            patch,
            box_lo,
            box_lo_workload,
            work_data_index,
            ideal_workload,
            workload_tolerance,
            min_size,
            cut_factor,
            bad_cut_points_for_boxlo);

         hier::BoxContainer boxes_hi;
         std::list<double> work_hi;

         double box_hi_workload = in_box_workload - box_lo_workload;
         privateRecursiveBisectionNonuniformSingleBox(boxes_hi,
            work_hi,
            patch,
            box_hi,
            box_hi_workload,
            work_data_index,
            ideal_workload,
            workload_tolerance,
            min_size,
            cut_factor,
            bad_cut_points_for_boxhi);

         out_boxes.spliceBack(boxes_hi);
         out_workloads.splice(out_workloads.end(), work_hi);

      } else {  // !can_cut_box

         out_boxes.pushFront(in_box);
         out_workloads.push_front(in_box_workload);

      }

   } // in_box_workload > ideal_workload

}

/*
 *************************************************************************
 *
 * Compute workload in box region of patch.
 *
 *************************************************************************
 */

double
BalanceUtilities::computeNonUniformWorkload(
   const std::shared_ptr<hier::Patch>& patch,
   int wrk_indx,
   const hier::Box& box)
{
   TBOX_ASSERT(patch);
   TBOX_ASSERT_OBJDIM_EQUALITY2(*patch, box);

   const std::shared_ptr<pdat::CellData<double> > work_data(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch->getPatchData(wrk_indx)));

   TBOX_ASSERT(work_data);

   double workload = s_norm_ops.L1Norm(work_data, box);

   return workload;
}

/*
 *************************************************************************
 *
 * Compute workload in patch while storing corner weights.
 *
 *************************************************************************
 */

double
BalanceUtilities::computeNonUniformWorkloadOnCorners(
   std::vector<double>& corner_weights,
   const std::shared_ptr<hier::Patch>& patch,
   int wrk_indx,
   const hier::Box& box)
{
   TBOX_ASSERT(patch);
   TBOX_ASSERT(box.getBlockId() == patch->getBox().getBlockId());
   TBOX_ASSERT(box.isSpatiallyEqual(patch->getBox()) ||
               box.isSpatiallyEqual(box * patch->getBox()));

   if (!corner_weights.empty()) {
      TBOX_ERROR("BalanceUtilities::computeNonUniformWorkloadOnCorners received a non-empty corner weights argument.");
   }

   double workload = 0.0;
   hier::Box work_box(box);
   hier::IntVector box_width(box.numberCells());
   hier::IntVector mid_point(box_width/2);
   const tbox::Dimension& dim = box.getDim();

   hier::IntVector corner_id(hier::IntVector::getZero(dim));
   do {
      // Set work_box to current corner.
      for (unsigned int d = 0; d < dim.getValue(); ++d) {
         if (corner_id[d] == 0) {
            work_box.setLower(d, box.lower()[d]);
            work_box.setUpper(d, box.lower()[d] + mid_point[d] - 1);
         } else {
            work_box.setLower(d, box.lower()[d] + mid_point[d]);
            work_box.setUpper(d, box.upper()[d]);
         }
      }

      // Increment n-dimensional corner id
      for (unsigned int d = 0; d < dim.getValue(); ++d) {
         if (corner_id[d] == 0) {
            corner_id[d] = 1;
            break;
         } else {
            corner_id[d] = 0;
         }
      }

      // Compute workload for current corner
      corner_weights.push_back(
         computeNonUniformWorkload(patch, wrk_indx, work_box));

      // End loop when all corners have been evaluated.
   } while (corner_id != hier::IntVector::getZero(dim));

   /*
    * Corner weights currently hold absolute workload.  The sum is the
    * the patch.
    */
   for (std::vector<double>::const_iterator itr = corner_weights.begin();
        itr != corner_weights.end(); ++itr) {
      workload += *itr;
   }

   /*
    * Change the corner weights to the fractions of the total workload.
    */
   for (std::vector<double>::iterator itr = corner_weights.begin();
        itr != corner_weights.end(); ++itr) {
      *itr /= workload;
   }

   return workload;
}


/*
 *************************************************************************
 *
 * Construct a processor mapping using the collection of weights and
 * the number of processors.  The return value provides an estimate of
 * the load balance efficiency from zero through one hundred.
 *
 *************************************************************************
 */

double
BalanceUtilities::binPack(
   hier::ProcessorMapping& mapping,
   std::vector<double>& weights,
   const int nproc)
{
   TBOX_ASSERT(nproc > 0);

   /*
    * Create the mapping array, find the average workload, and zero weights
    */

   const int nboxes = static_cast<int>(weights.size());
   mapping.setMappingSize(nboxes);

   double avg_work = 0.0;
   for (int w = 0; w < nboxes; ++w) {
      TBOX_ASSERT(weights[w] >= 0.0);
      avg_work += weights[w];
   }
   avg_work /= nproc;

   std::vector<double> work(nproc);
   for (int p = 0; p < nproc; ++p) {
      work[p] = 0.0;
   }

   /*
    * Assign each box to the processor with the lowest workload
    */
   for (int b = 0; b < nboxes; ++b) {
      const double weight = weights[b];

      int proc = 0;
      double diff = avg_work - (work[0] + weight);
      for (int p = 1; p < nproc; ++p) {
         const double d = avg_work - (work[p] + weight);
         if (((diff > 0.0) && (d >= 0.0) &&
              (d < diff)) || ((diff < 0.0) && (d > diff))) {
            diff = d;
            proc = p;
         }
      }

      work[proc] += weight;
      mapping.setProcessorAssignment(b, proc);
   }

   /*
    * Estimate load balance efficiency
    */

   double max_work = 0.0;
   for (int iw = 0; iw < nproc; ++iw) {
      if (work[iw] > max_work) max_work = work[iw];
   }

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
   return max_work == 0.0 ? 100.0 : 100.0 * avg_work / max_work;

}

/*
 *************************************************************************
 *
 * Construct a processor mapping using the collection of weights,
 * the number of processors, and the spatial distribution of the boxes.
 * The algorithm has two steps: (1) order the boxes based on the
 * location of the center of the box and (2) partition the boxes in
 * order so that the workload is fairly evenly distributed over the
 * processors.
 * The return value provides an estimate of the load balance efficiency
 * from zero through one hundred percent.
 *
 *************************************************************************
 */

double
BalanceUtilities::spatialBinPack(
   hier::ProcessorMapping& mapping,
   std::vector<double>& weights,
   hier::BoxContainer& boxes,
   const int nproc)
{
   TBOX_ASSERT(nproc > 0);
   TBOX_ASSERT(static_cast<int>(weights.size()) == boxes.size());

   const int nboxes = boxes.size();

   /*
    * compute offset which guarantees that the index space for all boxes
    * is positive.
    */
   hier::Index offset(boxes.front().lower());
   for (hier::BoxContainer::iterator itr = boxes.begin();
        itr != boxes.end(); ++itr) {
      offset.min(itr->lower());
   }

   /* construct array of spatialKeys */
   std::vector<SpatialKey> spatial_keys(nboxes);
   int i = 0;

   if (nboxes > 0) {
      const tbox::Dimension& dim = boxes.front().getDim();

      for (hier::BoxContainer::iterator itr = boxes.begin();
           itr != boxes.end(); ++itr) {

         /* compute center of box */
         hier::Index center = (itr->upper() + itr->lower()) / 2;

         if (dim == tbox::Dimension(1)) {
            spatial_keys[i].setKey(center(0) - offset(0));
         } else if (dim == tbox::Dimension(2)) {
            spatial_keys[i].setKey(center(0) - offset(0), center(1) - offset(1));
         } else if (dim == tbox::Dimension(3)) {
            spatial_keys[i].setKey(center(0) - offset(0), center(1) - offset(1),
               center(2) - offset(2));
         } else {
            TBOX_ERROR("BalanceUtilities::spatialBinPack error ..."
               << "\n not implemented for DIM>3" << std::endl);
         }
         ++i;
      }
   }

   /*
    * Sort boxes according to their spatial keys using a heapsort.
    */

   std::vector<int> permutation(nboxes);

   for (i = 0; i < nboxes; ++i) {
      permutation[i] = i;
   }

   for (i = nboxes / 2 - 1; i >= 0; --i) {
      privateHeapify(permutation, spatial_keys, i, nboxes);
   }
   for (i = nboxes - 1; i >= 1; --i) {
      const int tmp = permutation[0];
      permutation[0] = permutation[i];
      permutation[i] = tmp;
      privateHeapify(permutation, spatial_keys, 0, i);
   }

   /*
    * Copy unsorted data into temporaries and permute into sorted order
    */
   if (nboxes > 0) {
      const tbox::Dimension& dim = boxes.front().getDim();

      std::vector<hier::Box> unsorted_boxes(nboxes, hier::Box(dim));
      std::vector<double> unsorted_weights(nboxes);

      i = 0;
      for (hier::BoxContainer::iterator itr = boxes.begin();
           itr != boxes.end(); ++itr) {
         unsorted_boxes[i] = *itr;
         unsorted_weights[i] = weights[i];
         ++i;
      }

      i = 0;
      for (hier::BoxContainer::iterator itr = boxes.begin();
           itr != boxes.end(); ++itr) {
         *itr = unsorted_boxes[permutation[i]];
         weights[i] = unsorted_weights[permutation[i]];
         ++i;
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Verify that the spatial keys are sorted in non-decreasing order
       */

      std::vector<SpatialKey> unsorted_keys(nboxes);
      for (i = 0; i < nboxes; ++i) {
         unsorted_keys[i] = spatial_keys[i];
      }

      for (i = 0; i < nboxes; ++i) {
         spatial_keys[i] = unsorted_keys[permutation[i]];
      }

      for (i = 0; i < nboxes - 1; ++i) {
         TBOX_ASSERT(spatial_keys[i] <= spatial_keys[i + 1]);
      }
#endif

   }

   /* Find average workload */

   double avg_work = 0.0;
   for (i = 0; i < nboxes; ++i) {
      TBOX_ASSERT(weights[i] >= 0.0);
      avg_work += weights[i];
   }
   avg_work /= nproc;

   /*
    * Generate processor mapping.  (nboxes-1) as the maximum
    * processor number assignable if the number of processors
    * exceeds the number of boxes.
    */
   mapping.setMappingSize(nboxes);
   if (nproc <= nboxes) {
      privateRecursiveProcAssign(0, nboxes - 1, weights,
         0, nproc - 1, mapping, avg_work);
   } else {
      privateRecursiveProcAssign(0, nboxes - 1, weights,
         0, nboxes - 1, mapping, avg_work);
   }

   /* compute work load for each processor */
   std::vector<double> work(nproc);
   for (i = 0; i < nproc; ++i) {
      work[i] = 0.0;
   }
   for (i = 0; i < nboxes; ++i) {
      work[mapping.getProcessorAssignment(i)] += weights[i];
   }

   /*
    * Estimate load balance efficiency
    */

   double max_work = 0.0;
   for (i = 0; i < nproc; ++i) {
      if (work[i] > max_work) max_work = work[i];
   }

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
   return max_work == 0.0 ? 100.0 : 100.0 * avg_work / max_work;

}

/*
 **************************************************************************
 *
 * Recursively chops boxes in input list until all pieces have workload
 * less than the prescribed ideal workload or no more chopping is allowed
 * by the given constraints.   A spatially-uniform workload is assumed;
 * i.e., all cells are weighted equally.
 *
 **************************************************************************
 */

void
BalanceUtilities::recursiveBisectionUniform(
   hier::BoxContainer& out_boxes,
   std::list<double>& out_workloads,
   const hier::BoxContainer& in_boxes,
   const double ideal_workload,
   const double workload_tolerance,
   const hier::IntVector& min_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(min_size, cut_factor, bad_interval);

   const tbox::Dimension dim(min_size.getDim());

   TBOX_ASSERT(ideal_workload > 0);
   TBOX_ASSERT((workload_tolerance >= 0) && (workload_tolerance < 1.0));
   TBOX_ASSERT(min_size > hier::IntVector::getZero(dim));
   TBOX_ASSERT(cut_factor > hier::IntVector::getZero(dim));
   TBOX_ASSERT(bad_interval >= hier::IntVector::getZero(dim));
   TBOX_ASSERT(!physical_domain.empty());

   out_boxes.clear();
   out_workloads.clear();

   bool bad_domain_boundaries_exist =
      privateBadCutPointsExist(physical_domain);

   for (hier::BoxContainer::const_iterator ib = in_boxes.begin();
        ib != in_boxes.end(); ++ib) {

      hier::Box box2chop = *ib;

      TBOX_ASSERT(!box2chop.empty());

      double boxwork = static_cast<double>(box2chop.size());

      if (boxwork <= ((1.0 + workload_tolerance) * ideal_workload)) {

         out_boxes.pushFront(box2chop);
         out_workloads.push_front(boxwork);

      } else {

         std::vector<std::vector<bool> > bad_cut_points(dim.getValue());

         privateInitializeBadCutPointsForBox(bad_cut_points,
            box2chop,
            bad_domain_boundaries_exist,
            bad_interval,
            physical_domain);

         hier::IntVector box_cut_factor(cut_factor.getDim());
         if (cut_factor.getNumBlocks() == 1) {
            box_cut_factor = cut_factor;
         } else {
            box_cut_factor = cut_factor.getBlockVector(box2chop.getBlockId());
         }

         hier::BoxContainer tempboxes;
         std::list<double> temploads;
         privateRecursiveBisectionUniformSingleBox(
            tempboxes,
            temploads,
            box2chop,
            boxwork,
            ideal_workload,
            workload_tolerance,
            min_size,
            box_cut_factor,
            bad_cut_points);

         out_boxes.spliceBack(tempboxes);
         out_workloads.splice(out_workloads.end(), temploads);

      }
   }

}

/*
 **************************************************************************
 *
 * Recursively chops boxes described by patches in input patch level
 * until all pieces have workload less than the prescribed ideal workload
 * or no more chopping is allowed by the given constraints.  A spatially-
 * nonuniform workload is assumed; i.e., cell weights are given by the
 * patch data defined by the weight patch data id.
 *
 **************************************************************************
 */

void
BalanceUtilities::recursiveBisectionNonuniform(
   hier::BoxContainer& out_boxes,
   std::list<double>& out_workloads,
   const std::shared_ptr<hier::PatchLevel>& in_level,
   int work_id,
   double ideal_workload,
   const double workload_tolerance,
   const hier::IntVector& min_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(min_size, cut_factor, bad_interval);

   const tbox::Dimension dim(min_size.getDim());

   TBOX_ASSERT(ideal_workload > 0);
   TBOX_ASSERT((workload_tolerance >= 0) && (workload_tolerance < 1.0));
   TBOX_ASSERT(min_size > hier::IntVector::getZero(dim));
   TBOX_ASSERT(cut_factor > hier::IntVector::getZero(dim));
   TBOX_ASSERT(bad_interval >= hier::IntVector::getZero(dim));
   TBOX_ASSERT(!physical_domain.empty());

   out_boxes.clear();
   out_workloads.clear();

   bool bad_domain_boundaries_exist =
      privateBadCutPointsExist(physical_domain);

   for (hier::PatchLevel::iterator ip(in_level->begin());
        ip != in_level->end(); ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      hier::Box box2chop = patch->getBox();

      double boxwork = computeNonUniformWorkload(patch,
            work_id,
            box2chop);

      if (boxwork <= ((1. + workload_tolerance) * ideal_workload)) {

         out_boxes.pushFront(box2chop);
         out_workloads.push_front(boxwork);

      } else {

         std::vector<std::vector<bool> > bad_cut_points(dim.getValue());

         privateInitializeBadCutPointsForBox(bad_cut_points,
            box2chop,
            bad_domain_boundaries_exist,
            bad_interval,
            physical_domain);

         hier::IntVector box_cut_factor(cut_factor.getDim());
         if (cut_factor.getNumBlocks() == 1) {
            box_cut_factor = cut_factor;
         } else {
            box_cut_factor = cut_factor.getBlockVector(box2chop.getBlockId());
         }

         hier::BoxContainer tempboxes;
         std::list<double> temploads;
         privateRecursiveBisectionNonuniformSingleBox(
            tempboxes,
            temploads,
            patch,
            box2chop,
            boxwork,
            work_id,
            ideal_workload,
            workload_tolerance,
            min_size,
            box_cut_factor,
            bad_cut_points);

         out_boxes.spliceBack(tempboxes);
         out_workloads.splice(out_workloads.end(), temploads);
      }

   }

}

/*
 *************************************************************************
 *
 * Computes processor layout that corresponds, as closely as possible,
 * to the size of the supplied box.
 *
 * Inputs:
 *   num_procs - number of processors
 *   box - box which is to be layed out on processors
 * Output:
 *   proc_dist - vector describing processor layout
 *
 *************************************************************************
 */

void
BalanceUtilities::computeDomainDependentProcessorLayout(
   hier::IntVector& proc_dist,
   int num_procs,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(proc_dist, box);

   const tbox::Dimension& dim(proc_dist.getDim());

   tbox::Dimension::dir_t i;
   TBOX_ASSERT(num_procs > 0);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < dim.getValue(); ++i) {
      TBOX_ASSERT(box.numberCells(i) > 0);
   }
#endif

   /*
    * Initialize:
    *  - compute p[] - array of prime factors of tot_procs
    *  - set d[] initially set to box dims in each direction.
    *  - set proc_dist[] = 1 in each direction
    *  - set pnew[] (recomputed set of primes) initially to p
    */
   std::vector<int> p;
   privatePrimeFactorization(num_procs, p);

   hier::IntVector d = box.numberCells();
   for (i = 0; i < dim.getValue(); ++i) {
      proc_dist(i) = 1;
   }

   std::vector<int> pnew;
   pnew.resize(static_cast<int>(p.size()));
   for (i = 0; i < static_cast<int>(p.size()); ++i) {
      pnew[i] = p[i];
   }
   privateResetPrimesArray(pnew);

   /*
    *  main loop:  build up processor count with prime
    *              factors until # processors is reached
    *              or we have run out of prime factors.
    *  NOTE:  infinite loop conditions occur if no
    *         box directions can be divided by any of the prime factors
    *         of num_procs.  Adding a counter prevents this condition.
    */
   int counter = 0;
   while ((proc_dist.getProduct() < num_procs) &&
          (pnew.size() > 0) && (counter < num_procs)) {

      //  Loop over prime factors - largest to smallest
      for (int k = static_cast<int>(pnew.size()) - 1; k >= 0; --k) {

         //  determine i - direction in which d is largest
         i = 0;
         int nx = d[i];
         for (tbox::Dimension::dir_t j = 0; j < dim.getValue(); ++j) {
            if (d[j] > nx) i = j;
         }

         // Divide the length by the largest possible prime
         // factor and update processors accordingly. Remove the
         // chosen prime factor from the prime factors array.
         if (d[i] % pnew[k] == 0) {
            d[i] = d[i] / pnew[k];
            proc_dist[i] = proc_dist[i] * pnew[k];

            // Once a prime factor is used, remove it from the array of
            // primes by setting to one and calling the privateResetPrimesArray
            // (which removes any prime = 1).
            pnew[k] = 1;
            privateResetPrimesArray(pnew);
         }

         // Check if our iteration to build processor counts has
         // reach total_procs count.  If so, break out of loop.
         if (proc_dist.getProduct() == num_procs) break;

      } // loop over prime factors

      ++counter;
   } // while loop

   /*
    * This routine can fail under certain circumstances, such as
    * when no box direction exactly divides by any of the prime factors.
    * In this case, revert to the less stringent routine which simply
    * breaks up the domain into prime factors.
    */
   if (proc_dist.getProduct() != num_procs) {
      TBOX_WARNING("computeDomainDependentProcessorLayout(): could not \n"
         << "construct valid processor array - calling \n"
         << "computeDomainIndependentProcessorLayout() " << std::endl);
      computeDomainIndependentProcessorLayout(proc_dist, num_procs, box);
   }

}

/*
 *************************************************************************
 *
 * Computes processor layout that simply uses largest prime factors in
 * the decomposition.  The box is only used to determine the largest
 * size in each direction.  The processor decomposition will NOT
 * necessarily correspond to box dimensions.
 *
 *************************************************************************
 */

void
BalanceUtilities::computeDomainIndependentProcessorLayout(
   hier::IntVector& proc_dist,
   int num_procs,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(proc_dist, box);

   tbox::Dimension::dir_t i;
   const tbox::Dimension& dim(proc_dist.getDim());

   TBOX_ASSERT(num_procs > 0);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (i = 0; i < dim.getValue(); ++i) {
      TBOX_ASSERT(box.numberCells(i) > 0);
   }
#endif

   /*
    * Input:
    *   box - box which is to be layed out on processors
    * Returned:
    *   proc_dist - processors in each direction
    *
    *
    * Initialize:
    *  - compute p[] - array of prime factors of num_procs
    *  - set d[] initially set to box dims in each direction.
    *  - set proc_dist[] = 1 in each direction
    *  - set pnew[] (recomputed set of primes) initially to p
    */

   std::vector<int> p;
   privatePrimeFactorization(num_procs, p);

   hier::IntVector d = box.numberCells();
   for (i = 0; i < dim.getValue(); ++i) {
      proc_dist(i) = 1;
   }

   std::vector<int> pnew;
   pnew.resize(static_cast<int>(p.size()));
   for (i = 0; i < static_cast<int>(p.size()); ++i) pnew[i] = p[i];
   privateResetPrimesArray(pnew);

   /*
    *  main loop:  build up processor count with prime
    *              factors until # processors is reached
    *              or we have run out of prime factors.
    */
   while ((proc_dist.getProduct() < num_procs) && (pnew.size() > 0)) {

      //  determine i - direction in which d is largest
      i = 0;
      int nx = d[i];
      for (tbox::Dimension::dir_t j = 0; j < dim.getValue(); ++j) {
         if (d[j] > nx) i = j;
      }

      // Set proc_dist(i) to the largest prime factor
      int k = static_cast<int>(pnew.size()) - 1;
      d[i] = d[i] / pnew[k];
      proc_dist[i] = proc_dist[i] * pnew[k];

      // Once a prime factor is used, remove it from the array of
      // primes by setting to one and calling the privateResetPrimesArray
      // (which removes any prime = 1).
      pnew[k] = 1;
      privateResetPrimesArray(pnew);

      // Check if our iteration to build processor count has
      // reach total_procs count.  If so, break out of loop.
      if (proc_dist.getProduct() == num_procs) break;

   } // while loop

   // error check - in case any cases pop up where we have
   // run out of prime factors but have not yet built a valid
   // processor array.
   if (proc_dist.getProduct() != num_procs) {
      TBOX_ERROR(
         "BalanceUtilities::computeDomainIndependentProcessorLayout() error"
         << "\n  invalid processor array computed" << std::endl);
   }

}

/*
 *************************************************************************
 *
 * Sort box work loads in decreasing order using a heapsort.  Both
 * the box array and the work array will be returned in sorted order.
 *
 *************************************************************************
 */

void
BalanceUtilities::sortDescendingBoxWorkloads(
   hier::BoxContainer& boxes,
   std::vector<double>& workload)
{
   TBOX_ASSERT(boxes.size() == static_cast<int>(workload.size()));

   /*
    * Create the permutation array that represents indices in sorted order
    */

   const int nboxes = static_cast<int>(workload.size());
   std::vector<int> permutation(nboxes);

   for (int i = 0; i < nboxes; ++i) {
      permutation[i] = i;
   }

   /*
    * Execute the heapsort using static member function privateHeapify()
    */

   for (int j = nboxes / 2 - 1; j >= 0; --j) {
      privateHeapify(permutation, workload, j, nboxes);
   }
   for (int k = nboxes - 1; k >= 1; --k) {
      const int tmp = permutation[0];
      permutation[0] = permutation[k];
      permutation[k] = tmp;
      privateHeapify(permutation, workload, 0, k);
   }

   /*
    * Copy unsorted data into temporaries and permute into sorted order
    */
   if (nboxes > 0) {
      const tbox::Dimension& dim(boxes.front().getDim());

      std::vector<hier::Box> unsorted_boxes(nboxes, hier::Box(dim));
      std::vector<double> unsorted_workload(nboxes);

      int l = 0;
      for (hier::BoxContainer::iterator itr = boxes.begin();
           itr != boxes.end(); ++itr) {
         unsorted_boxes[l] = *itr;
         unsorted_workload[l] = workload[l];
         ++l;
      }

      int m = 0;
      for (hier::BoxContainer::iterator itr = boxes.begin();
           itr != boxes.end(); ++itr) {
         *itr = unsorted_boxes[permutation[m]];
         workload[m] = unsorted_workload[permutation[m]];
         ++m;
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Verify that the workload is sorted in nonincreasing order
       */

      for (int n = 0; n < nboxes - 1; ++n) {
         TBOX_ASSERT(workload[n] >= workload[n + 1]);
      }
#endif
   }
}

/*
 **************************************************************************
 *
 * Compute and return load balance efficiency for a level.
 *
 **************************************************************************
 */

double
BalanceUtilities::computeLoadBalanceEfficiency(
   const std::shared_ptr<hier::PatchLevel>& level,
   std::ostream& os,
   int workload_data_id)
{
   TBOX_ASSERT(level);

   NULL_USE(os);
   const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

   int i;

   const hier::ProcessorMapping& mapping = level->getProcessorMapping();

   const int nprocs = mpi.getSize();
   std::vector<double> work(nprocs);

   for (i = 0; i < nprocs; ++i) {
      work[i] = 0.0;
   }

   if (workload_data_id < 0) {

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         work[mapping.getProcessorAssignment(ip->getLocalId().getValue())] +=
            static_cast<double>((*ip)->getBox().size());
      }

   } else {

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         std::shared_ptr<pdat::CellData<double> > weight(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(workload_data_id)));

         TBOX_ASSERT(weight);

         work[mapping.getProcessorAssignment(ip->getLocalId().getValue())] +=
            s_norm_ops.L1Norm(weight, patch->getBox());
      }

   }

   if (mpi.getSize() > 1) {
      mpi.AllReduce(&work[0], 1, MPI_SUM);
   }

   double max_work = 0.0;
   double total_work = 0.0;
   for (i = 0; i < nprocs; ++i) {
      total_work += work[i];
      if (work[i] > max_work) max_work = work[i];
   }

   double efficiency = 100.0 * total_work / (max_work * nprocs);

   return efficiency;

}

/*
 *************************************************************************
 *************************************************************************
 */

void
BalanceUtilities::findSmallBoxesInPostbalance(
   std::ostream& co,
   const std::string& border,
   const hier::MappingConnector& post_to_pre,
   const hier::IntVector& min_width,
   size_t min_cells)
{
   const hier::BoxLevel& post = post_to_pre.getBase();
   const hier::BoxContainer& post_boxes = post.getBoxes();

   int local_new_min_count = 0;
   for (hier::BoxContainer::const_iterator bi = post_boxes.begin(); bi != post_boxes.end(); ++bi) {

      const hier::Box& post_box = *bi;

      if (post_box.numberCells() >= min_width && static_cast<size_t>(post_box.size()) >=
          min_cells) {
         continue;
      }

      if (post_to_pre.hasNeighborSet(post_box.getBoxId())) {
         hier::BoxContainer pre_neighbors;
         post_to_pre.getLocalNeighbors(pre_neighbors);

         bool small_width = true;
         bool small_cells = true;
         for (hier::BoxContainer::const_iterator na = pre_neighbors.begin();
              na != pre_neighbors.end(); ++na) {
            if (!(na->numberCells() >= min_width)) {
               small_width = false;
            }
            if (!(static_cast<size_t>(na->size()) >= min_cells)) {
               small_cells = false;
            }
         }
         if (small_width || small_cells) {
            ++local_new_min_count;
            co << border << "Post-box small_width=" << small_width
               << " small_cells=" << small_cells
               << ": " << post_box
               << post_box.numberCells() << '|' << post_box.size()
               << " from " << pre_neighbors.format(border, 2);
         }

      }

   }

   int global_new_min_count = local_new_min_count;
   post.getMPI().AllReduce(&global_new_min_count, 1, MPI_SUM);

   co << border
      << "  Total of " << local_new_min_count << " / "
      << global_new_min_count << " new minimums." << std::endl;
}

/*
 *************************************************************************
 *************************************************************************
 */

void
BalanceUtilities::findSmallBoxesInPostbalance(
   std::ostream& co,
   const std::string& border,
   const hier::BoxLevel& post,
   const hier::BoxLevel& pre,
   const hier::IntVector& min_width,
   size_t min_cells)
{
   hier::MappingConnector post_to_pre(post, pre, hier::IntVector::getZero(post.getDim()));
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps(post_to_pre);
   findSmallBoxesInPostbalance(co, border, post_to_pre, min_width, min_cells);
}

/*
 *************************************************************************
 *************************************************************************
 */

bool
BalanceUtilities::compareLoads(
   int flags[],
   double cur_load,
   double new_load,
   double ideal_load,
   double low_load,
   double high_load,
   const PartitioningParams& pparams)
{
   double cur_range_miss = cur_load >= high_load ? cur_load - high_load :
      (cur_load <= low_load ? low_load - cur_load : 0.0);
   double new_range_miss = new_load >= high_load ? new_load - high_load :
      (new_load <= low_load ? low_load - new_load : 0.0);
   flags[0] = new_range_miss < (cur_range_miss - pparams.getLoadComparisonTol()) ? 1 :
      (new_range_miss > cur_range_miss ? -1 : 0);

   double cur_diff = tbox::MathUtilities<double>::Abs(cur_load - ideal_load);
   double new_diff = tbox::MathUtilities<double>::Abs(new_load - ideal_load);

   flags[1] = new_diff < (cur_diff - pparams.getLoadComparisonTol()) ? 1 :
      (new_diff > cur_diff ? -1 : 0);

   flags[2] = flags[0] != 0 ? flags[0] : (flags[1] != 0 ? flags[1] : 0);

   flags[3] = (new_load <= high_load && new_load >= low_load);

   return flags[2] == 1;
}

/*
 *************************************************************************
 *************************************************************************
 */
void
BalanceUtilities::reduceAndReportLoadBalance(
   const std::vector<double>& loads,
   const tbox::SAMRAI_MPI& mpi,
   std::ostream& os)
{
   const int nproc = mpi.getSize();

   const double demarks[] = { 0.50,
                              0.70,
                              0.85,
                              0.92,
                              0.98,
                              1.02,
                              1.08,
                              1.15,
                              1.30,
                              1.50,
                              2.00 };
   const int ndemarks = 11;

   // Compute total, avg, min and max loads.

   std::vector<double> min_loads(loads);
   std::vector<int> min_ranks(loads.size());
   mpi.AllReduce(&min_loads[0], static_cast<int>(min_loads.size()), MPI_MINLOC, &min_ranks[0]);

   std::vector<double> max_loads(loads);
   std::vector<int> max_ranks(loads.size());
   mpi.AllReduce(&max_loads[0], static_cast<int>(max_loads.size()), MPI_MAXLOC, &max_ranks[0]);

   std::vector<double> total_loads(loads);
   mpi.AllReduce(&total_loads[0], static_cast<int>(total_loads.size()), MPI_SUM);

   const int n_population_zones = ndemarks + 1;
   std::vector<int> population(loads.size() * n_population_zones, 0);
   for (size_t iload = 0; iload < loads.size(); ++iload) {
      int izone;
      for (izone = 0; izone < ndemarks; ++izone) {
         if (loads[iload] / total_loads[iload] * mpi.getSize() < demarks[izone]) {
            break;
         }
      }
      population[iload * n_population_zones + izone] = 1;
   }
   mpi.AllReduce(&population[0], static_cast<int>(population.size()), MPI_SUM);

   for (size_t iload = 0; iload < loads.size(); ++iload) {

      const double total_load = total_loads[iload];
      const double avg_load = total_loads[iload] / mpi.getSize();
      const double min_load = min_loads[iload];
      const int r_min_load = min_ranks[iload];
      const double max_load = max_loads[iload];
      const int r_max_load = max_ranks[iload];

      os << "================ Sequence " << iload << " ===============\n";
      os << "total/avg loads: "
         << total_load << " / "
         << avg_load << "\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
      os << std::setprecision(6)
         << "min/max loads: "
         << min_load << " @ P" << r_min_load << " / "
         << max_load << " @ P" << r_max_load << "   "
         << "diffs: "
         << min_load - avg_load << " / "
         << max_load - avg_load << "   "
         << std::setprecision(4)
         << "normalized: "
         << (avg_load != 0 ? min_load / avg_load : 0.0) << " / "
         << (avg_load != 0 ? max_load / avg_load : 0.0) << "\n";

      const char bars[] = "----";
      const char space[] = "   ";
      os.setf(std::ios_base::fixed);
      os << bars;
      for (int izone = 0; izone < ndemarks; ++izone) {
         os << std::setw(4) << std::setprecision(2) << demarks[izone] << bars;
      }
      os << '\n';

#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
      for (int izone = 0; izone < n_population_zones; ++izone) {
         const int nrank = population[iload * n_population_zones + izone];
         const double percentage = 100.0 * nrank / nproc;
         os << std::setw(5) << percentage << space;
      }
      os << '\n';

   }

}

/*
 *************************************************************************
 * for use when sorting loads using the C-library qsort
 *************************************************************************
 */
int
BalanceUtilities::qsortRankAndLoadCompareDescending(
   const void* v,
   const void* w)
{
   const RankAndLoad* lv = (const RankAndLoad *)v;
   const RankAndLoad* lw = (const RankAndLoad *)w;
   if (lv->load > lw->load) {
      return -1;
   }

   if (lv->load < lw->load) {
      return 1;
   }

   return 0;
}

/*
 *************************************************************************
 * for use when sorting loads using the C-library qsort
 *************************************************************************
 */
int
BalanceUtilities::qsortRankAndLoadCompareAscending(
   const void* v,
   const void* w)
{
   const RankAndLoad* lv = (const RankAndLoad *)v;
   const RankAndLoad* lw = (const RankAndLoad *)w;
   if (lv->load < lw->load) {
      return -1;
   }

   if (lv->load > lw->load) {
      return 1;
   }

   return 0;
}

/*
 *************************************************************************
 * Constrain maximum box sizes in the given BoxLevel and
 * update given Connectors to the changed BoxLevel.
 *************************************************************************
 */
void
BalanceUtilities::constrainMaxBoxSizes(
   hier::BoxLevel& box_level,
   hier::Connector* anchor_to_level,
   const PartitioningParams& pparams)
{
   TBOX_ASSERT(!anchor_to_level || anchor_to_level->hasTranspose());

   const hier::IntVector& zero_vector(hier::IntVector::getZero(box_level.getDim()));

   hier::BoxLevel constrained(box_level.getRefinementRatio(),
                              box_level.getGridGeometry(),
                              box_level.getMPI());
   hier::MappingConnector unconstrained_to_constrained(box_level,
                                                       constrained,
                                                       zero_vector);

   const hier::BoxContainer& unconstrained_boxes = box_level.getBoxes();

   hier::LocalId next_available_index = box_level.getLastLocalId() + 1;

   for (hier::BoxContainer::const_iterator ni = unconstrained_boxes.begin();
        ni != unconstrained_boxes.end(); ++ni) {

      const hier::Box& box = *ni;

      const hier::IntVector box_size = box.numberCells();

      /*
       * If box already conform to max size constraint, keep it.
       * Else chop it up and keep the parts.
       */

      if (box_size <= pparams.getMaxBoxSize()) {

         constrained.addBox(box);

      } else {

         hier::BoxContainer chopped(box);
         hier::BoxUtilities::chopBoxes(
            chopped,
            pparams.getMaxBoxSize(),
            pparams.getMinBoxSize(),
            pparams.getCutFactor(),
            pparams.getBadInterval(),
            pparams.getDomainBoxes(box.getBlockId()));
         TBOX_ASSERT(!chopped.empty());

         if (chopped.size() != 1) {

            hier::Connector::NeighborhoodIterator base_box_itr =
               unconstrained_to_constrained.makeEmptyLocalNeighborhood(
                  box.getBoxId());

            for (hier::BoxContainer::iterator li = chopped.begin();
                 li != chopped.end(); ++li) {

               const hier::Box fragment = *li;

               const hier::Box new_box(fragment,
                                       next_available_index++,
                                       box_level.getMPI().getRank());
               TBOX_ASSERT(new_box.getBlockId() == ni->getBlockId());

               constrained.addBox(new_box);

               unconstrained_to_constrained.insertLocalNeighbor(
                  new_box,
                  base_box_itr);

            }

         } else {
            TBOX_ASSERT(box.isSpatiallyEqual(chopped.front()));
            constrained.addBox(box);
         }

      }

   }

   if (anchor_to_level && anchor_to_level->isFinalized()) {
      // Modify anchor<==>level Connectors and swap box_level with constrained.
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix("mesh::BalanceUtilities");
      mca.modify(*anchor_to_level,
         unconstrained_to_constrained,
         &box_level,
         &constrained);
   } else {
      // Swap box_level and constrained without touching anchor<==>level.
      hier::BoxLevel::swap(box_level, constrained);
   }

}

/*
 **************************************************************************
 * Move Boxes in balance_box_level from ranks outside of
 * rank_group to ranks inside rank_group.  Modify the given connectors
 * to make them correct following this moving of boxes.
 **************************************************************************
 */

void
BalanceUtilities::prebalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const tbox::RankGroup& rank_group)
{

   if (balance_to_anchor) {
      TBOX_ASSERT(balance_to_anchor->hasTranspose());
      TBOX_ASSERT(balance_to_anchor->getTranspose().checkTransposeCorrectness(
            *balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor->checkTransposeCorrectness(
            balance_to_anchor->getTranspose()) == 0);
   }

   /*
    * tmp_box_level will contain the same boxes as
    * balance_box_level, but all will live on the processors
    * specified in rank_group.
    */
   hier::BoxLevel tmp_box_level(balance_box_level.getRefinementRatio(),
                                balance_box_level.getGridGeometry(),
                                balance_box_level.getMPI());

   /*
    * If a rank is not in rank_group it is called a "sending" rank, as
    * it will send any Boxes it has to a rank in rank_group.
    */
   bool is_sending_rank = rank_group.isMember(balance_box_level.getMPI().getRank()) ? false : true;

   int output_nproc = rank_group.size();

   /*
    * the send and receive comm objects
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* box_send = 0;
   tbox::AsyncCommPeer<int>* box_recv = 0;
   tbox::AsyncCommPeer<int>* id_send = 0;
   tbox::AsyncCommPeer<int>* id_recv = 0;

   /*
    * A sending rank will send its Boxes to a receiving rank, and
    * that receiving processor will add it to its local set of Boxes.
    * When the box is added on the receiving processor, it will receive
    * a new LocalId.  This LocalId value needs to be sent back to
    * the sending processor, in order to construct the mapping connectors.
    *
    * Therefore the sending ranks construct comm objects for sending boxes
    * and receiving LocalIdes.
    *
    * Sending processors send to ranks in the rank_group determined by
    * a modulo heuristic.
    */
   if (is_sending_rank) {
      box_send = new tbox::AsyncCommPeer<int>;
      box_send->initialize(&comm_stage);
      box_send->setPeerRank(rank_group.getMappedRank(balance_box_level.getMPI().getRank()
            % output_nproc));
      box_send->setMPI(balance_box_level.getMPI());
      box_send->setMPITag(BalanceUtilities_PREBALANCE0 + 2 * balance_box_level.getMPI().getRank(),
         BalanceUtilities_PREBALANCE1 + 2 * balance_box_level.getMPI().getRank());

      id_recv = new tbox::AsyncCommPeer<int>;
      id_recv->initialize(&comm_stage);
      id_recv->setPeerRank(rank_group.getMappedRank(balance_box_level.getMPI().getRank()
            % output_nproc));
      id_recv->setMPI(balance_box_level.getMPI());
      id_recv->setMPITag(BalanceUtilities_PREBALANCE0 + 2 * balance_box_level.getMPI().getRank(),
         BalanceUtilities_PREBALANCE1 + 2 * balance_box_level.getMPI().getRank());
   }

   /*
    * The receiving ranks construct comm objects for receiving boxes
    * and sending LocalIdes.
    */
   int num_recvs = 0;
   if (rank_group.isMember(balance_box_level.getMPI().getRank())) {
      std::list<int> recv_ranks;
      for (int i = 0; i < balance_box_level.getMPI().getSize(); ++i) {
         if (!rank_group.isMember(i) &&
             rank_group.getMappedRank(i % output_nproc) == balance_box_level.getMPI().getRank()) {
            recv_ranks.push_back(i);
         }
      }
      num_recvs = static_cast<int>(recv_ranks.size());
      if (num_recvs > 0) {
         box_recv = new tbox::AsyncCommPeer<int>[num_recvs];
         id_send = new tbox::AsyncCommPeer<int>[num_recvs];
         int recv_count = 0;
         for (std::list<int>::const_iterator ri(recv_ranks.begin());
              ri != recv_ranks.end(); ++ri) {
            const int rank = *ri;
            box_recv[recv_count].initialize(&comm_stage);
            box_recv[recv_count].setPeerRank(rank);
            box_recv[recv_count].setMPI(balance_box_level.getMPI());
            box_recv[recv_count].setMPITag(BalanceUtilities_PREBALANCE0 + 2 * rank,
               BalanceUtilities_PREBALANCE1 + 2 * rank);

            id_send[recv_count].initialize(&comm_stage);
            id_send[recv_count].setPeerRank(rank);
            id_send[recv_count].setMPI(balance_box_level.getMPI());
            id_send[recv_count].setMPITag(BalanceUtilities_PREBALANCE0 + 2 * rank,
               BalanceUtilities_PREBALANCE1 + 2 * rank);

            ++recv_count;
         }
         TBOX_ASSERT(num_recvs == recv_count);
      }
   }

   /*
    * Construct the mapping Connectors which describe the mapping from the box
    * configuration of the given balance_box_level, to the new
    * configuration stored in tmp_box_level.  These mapping Connectors
    * are necessary to modify the two Connectors given in the argument list,
    * so that on return from this method, they will be correct for the new
    * balance_box_level.
    */
   const hier::IntVector& zero_vector(hier::IntVector::getZero(balance_box_level.getDim()));
   hier::MappingConnector balance_to_tmp(
      balance_box_level,
      tmp_box_level,
      zero_vector);

   hier::MappingConnector tmp_to_balance(
      tmp_box_level,
      balance_box_level,
      zero_vector);

   balance_to_tmp.setTranspose(&tmp_to_balance, false);

   /*
    * Where Boxes already exist on ranks in rank_group,
    * move them directly to tmp_box_level.
    */
   if (!is_sending_rank) {
      const hier::BoxContainer& unchanged_boxes =
         balance_box_level.getBoxes();

      for (hier::BoxContainer::const_iterator ni = unchanged_boxes.begin();
           ni != unchanged_boxes.end(); ++ni) {

         const hier::Box& box = *ni;
         tmp_box_level.addBox(box);
      }
   }

   const int buf_size = hier::Box::commBufferSize(balance_box_level.getDim());

   /*
    * On sending ranks, pack the Boxes into buffers and send.
    */
   if (is_sending_rank) {
      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      const int num_sending_boxes =
         static_cast<int>(sending_boxes.size());

      int* buffer = 0;
      if (num_sending_boxes > 0) {
         buffer = new int[buf_size * num_sending_boxes];
      }
      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         const hier::Box& box = *ni;

         box.putToIntBuffer(&buffer[box_count * buf_size]);
         ++box_count;
      }
      box_send->beginSend(buffer, static_cast<size_t>(buf_size * num_sending_boxes));

      if (buffer) {
         delete[] buffer;
      }
   }

   /*
    * On receiving ranks, complete the receives, add the boxes to local
    * tmp_box_level, insert boxes into tmp_to_balance, and then
    * send the new LocalIdes back to the sending processors.
    */
   if (!is_sending_rank && num_recvs > 0) {
      for (int i = 0; i < num_recvs; ++i) {
         box_recv[i].beginRecv();
      }
      int num_completed_recvs = 0;
      std::vector<bool> completed(num_recvs, false);
      while (num_completed_recvs < num_recvs) {
         for (int i = 0; i < num_recvs; ++i) {
            if (!completed[i] && box_recv[i].checkRecv()) {
               ++num_completed_recvs;
               completed[i] = true;
               const int num_boxes = box_recv[i].getRecvSize() / buf_size;
               const int* buffer = box_recv[i].getRecvData();
               int* id_buffer = 0;
               if (num_boxes > 0) {
                  id_buffer = new int[num_boxes];
               }

               for (int b = 0; b < num_boxes; ++b) {
                  hier::Box box(balance_box_level.getDim());

                  box.getFromIntBuffer(&buffer[b * buf_size]);

                  hier::BoxContainer::const_iterator tmp_iter =
                     tmp_box_level.addBox(box,
                        box.getBlockId());

                  hier::BoxId tmp_box_id = tmp_iter->getBoxId();

                  tmp_to_balance.insertLocalNeighbor(box, tmp_box_id);

                  id_buffer[b] = tmp_box_id.getLocalId().getValue();
               }
               id_send[i].beginSend(id_buffer, static_cast<size_t>(num_boxes));

               if (id_buffer) {
                  delete[] id_buffer;
               }
            }
         }
      }
      for (int i = 0; i < num_recvs; ++i) {
         if (!id_send[i].checkSend()) {
            id_send[i].completeCurrentOperation();
         }
      }
   }

   /*
    * On sending ranks, receive the LocalIds, and add the edges
    * to balance_to_tmp.
    */
   if (is_sending_rank) {
      if (!box_send->checkSend()) {
         box_send->completeCurrentOperation();
      }

      id_recv->beginRecv();

      if (!id_recv->checkRecv()) {
         id_recv->completeCurrentOperation();
      }
      const int* buffer = id_recv->getRecvData();

      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      TBOX_ASSERT(static_cast<int>(id_recv->getRecvSize()) == sending_boxes.size());

      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         hier::Box new_box(
            *ni,
            (hier::LocalId)buffer[box_count],
            rank_group.getMappedRank(balance_box_level.getMPI().getRank() % output_nproc));

         balance_to_tmp.insertLocalNeighbor(new_box, (*ni).getBoxId());
         ++box_count;
      }
   }

   if (balance_to_anchor && balance_to_anchor->hasTranspose()) {
      /*
       * This modify operation copies tmp_box_level to
       * balance_box_level, and changes balance_to_anchor and
       * its transpose such that they are correct for the new state
       * of balance_box_level.
       */
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix("mesh::BalanceUtilities");
      mca.modify(balance_to_anchor->getTranspose(),
         balance_to_tmp,
         &balance_box_level,
         &tmp_box_level);

      TBOX_ASSERT(balance_to_anchor->getTranspose().checkTransposeCorrectness(
            *balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor->checkTransposeCorrectness(
            balance_to_anchor->getTranspose()) == 0);
   } else {
      hier::BoxLevel::swap(balance_box_level, tmp_box_level);
   }

   /*
    * Clean up raw pointer allocation.
    */
   if (is_sending_rank) {
      delete box_send;
      delete id_recv;
   }
   if (num_recvs) {
      delete[] box_recv;
      delete[] id_send;
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
