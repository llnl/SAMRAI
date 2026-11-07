/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Binary tree of Boxes for overlap searches.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxTree.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

std::shared_ptr<tbox::Timer> BoxTree::t_build_tree[SAMRAI::MAX_DIM_VAL];
std::shared_ptr<tbox::Timer> BoxTree::t_search[SAMRAI::MAX_DIM_VAL];
unsigned int BoxTree::s_num_build[SAMRAI::MAX_DIM_VAL] =
{ 0 };
unsigned int BoxTree::s_num_generate[SAMRAI::MAX_DIM_VAL]
   =
   { 0 };
unsigned int BoxTree::s_num_duplicate[SAMRAI::MAX_DIM_VAL]
   =
   { 0 };
unsigned int BoxTree::s_num_search[SAMRAI::MAX_DIM_VAL] =
{ 0 };
unsigned int BoxTree::s_num_sorted_box[SAMRAI::MAX_DIM_VAL
] = { 0 };
unsigned int BoxTree::s_num_found_box[SAMRAI::MAX_DIM_VAL]
   =
   { 0 };
unsigned int BoxTree::s_max_sorted_box[SAMRAI::MAX_DIM_VAL
] = { 0 };
unsigned int BoxTree::s_max_found_box[SAMRAI::MAX_DIM_VAL]
   =
   { 0 };
unsigned int BoxTree::s_max_lin_search[SAMRAI::MAX_DIM_VAL
] = { 0 };

tbox::StartupShutdownManager::Handler
BoxTree::s_initialize_finalize_handler(
   BoxTree::initializeCallback,
   0,
   0,
   BoxTree::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 * Constructor for an empty BoxTree
 *************************************************************************
 */

BoxTree::BoxTree(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId()),
   d_partition_dir(0)
{
}

/*
 *************************************************************************
 * Constructor taking a BoxContainer
 *************************************************************************
 */
BoxTree::BoxTree(
   const tbox::Dimension& dim,
   const BoxContainer& boxes,
   int min_number):
   d_dim(dim),
   d_bounding_box(dim),
   d_block_id(BlockId::invalidId())
{
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] +=
      static_cast<int>(boxes.size());
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         static_cast<int>(boxes.size()));
#ifndef _OPENMP
   t_build_tree[d_dim.getValue() - 1]->start();
#endif
   min_number = (min_number < 1) ? 1 : min_number;

#ifdef DEBUG_CHECK_ASSERTIONS
   // Catch empty boxes so sorting logic does not have to.
   for (BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      TBOX_ASSERT(!ni->empty());
   }
#endif

   /*
    * Implementation note: We can simply copy boxes into
    * d_boxes and call privateGenerateTree using:
    *
    *   d_boxes.insert(boxes.begin(),
    *                         boxes.end());
    *   privateGenerateTree(d_boxes, min_number);
    *
    * However, this extra copy slows things down about 30%.
    * So we live with the repetitious code to do the same thing
    * that privateGenerateTree does.
    */

   /*
    * Compute the bounding box for the set of boxes.  Also get
    * BlockId from the given boxes.
    */
   if (!boxes.empty()) {
      TBOX_ASSERT(boxes.begin()->getBlockId() != BlockId::invalidId());
      d_block_id = boxes.begin()->getBlockId();
   }
   for (BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end(); ++ni) {
      d_bounding_box += (*ni);
      TBOX_ASSERT(ni->getBlockId() == d_block_id);
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (boxes.size() <= min_number) {
      for (BoxContainer::const_iterator ni = boxes.begin();
           ni != boxes.end(); ++ni) {
         d_boxes.push_back(&(*ni));
      }
   } else {

      /*
       * Partition the boxes into three sets, using the midpoint of
       * the longest direction of the bounding box:
       *
       * - those that belong to me (intersects the midpoint plane).  Put
       * these in d_boxes.
       *
       * - those that belong to my left child (lower than the midpoint
       * plane)
       *
       * - those that belong to my right child (higher than the midpoint
       * plane)
       */

      const IntVector bbsize = d_bounding_box.numberCells();
      d_partition_dir = 0;
      for (tbox::Dimension::dir_t d = 1; d < d_dim.getValue(); ++d) {
         if (bbsize(d_partition_dir) < bbsize(d)) {
            d_partition_dir = d;
         }
      }

      int midpoint =
         (d_bounding_box.lower(d_partition_dir)
          + d_bounding_box.upper(d_partition_dir)) / 2;

      std::list<const Box *> left_boxes, right_boxes;
      for (BoxContainer::const_iterator ni = boxes.begin();
           ni != boxes.end(); ++ni) {
         const Box& box = *ni;
         if (box.upper(d_partition_dir) <= midpoint) {
            left_boxes.push_back(&box);
         } else if (box.lower(d_partition_dir) > midpoint) {
            right_boxes.push_back(&box);
         } else {
            d_boxes.push_back(&box);
         }
      }

      setupChildren(min_number, left_boxes, right_boxes);

   }

   if (s_max_lin_search[d_dim.getValue() - 1] <
       static_cast<unsigned int>(d_boxes.size())) {
      s_max_lin_search[d_dim.getValue() - 1] =
         static_cast<unsigned int>(d_boxes.size());
   }

#ifndef _OPENMP
   t_build_tree[d_dim.getValue() - 1]->stop();
#endif
}

BoxTree::BoxTree(
   const std::list<const Box *> boxes,
   int min_number):
   d_dim((*(boxes.begin()))->getDim()),
   d_bounding_box(d_dim),
   d_boxes(boxes)
{
   ++s_num_build[d_dim.getValue() - 1];
   s_num_sorted_box[d_dim.getValue() - 1] +=
      static_cast<int>(boxes.size());
   s_max_sorted_box[d_dim.getValue() - 1] = tbox::MathUtilities<int>::Max(
         s_max_sorted_box[d_dim.getValue() - 1],
         static_cast<int>(boxes.size()));
#ifndef _OPENMP
   t_build_tree[d_dim.getValue() - 1]->start();
#endif
   min_number = (min_number < 1) ? 1 : min_number;

   privateGenerateTree(min_number);

#ifndef _OPENMP
   t_build_tree[d_dim.getValue() - 1]->stop();
#endif
}



/*
 *************************************************************************
 * Generate the tree from the boxes in d_boxes.
 *
 * Methods taking various input containers of Boxes could
 * simply copy the input Boxes into a BoxContainer, then call this
 * method.  However, we don't do that for efficiency reasons.  The
 * extra copy turns out to be significant.  Therefore, the
 * constructors have code similar to privateGenerateTree to split
 * the incoming Boxes into three groups.  These groups
 * are turned into child branches by setupChildren.
 *
 * This method is not timed using the Timers.  Only the public
 * itnerfaces are timed.  Isolating the recursive code in
 * privateGenerateTree also helps in timing the methods, because timer
 * starts/stops can be removed from the recursive codes.
 *************************************************************************
 */
void
BoxTree::privateGenerateTree(
   int min_number)
{
   ++s_num_generate[d_dim.getValue() - 1];

   if (d_boxes.size() > 0) {
      d_block_id = (**(d_boxes.begin())).getBlockId();
   }

   /*
    * Compute this tree's domain, which is the bounding box for the
    * constituent boxes.
    */
   for (std::list<const Box *>::const_iterator ni = d_boxes.begin();
        ni != d_boxes.end(); ++ni) {
      d_bounding_box += **ni;
   }

   /*
    * If the list of boxes is small enough, we won't
    * do any recursive stuff: we'll just let the boxes
    * live here.  In this case, there is no left child,
    * no right child, and no recursive d_center_child.
    */
   if (d_boxes.size() > static_cast<std::list<const Box *>::size_type>(min_number)) {
      /*
       * Partition the boxes into three sets, using the midpoint of
       * the longest direction of the bounding box:
       *
       * - those that belong to me (intersects the midpoint plane).  Put
       * these in d_boxes.
       *
       * - those that belong to my left child (lower than the midpoint
       * plane)
       *
       * - those that belong to my right child (higher than the midpoint
       * plane)
       */

      const IntVector bbsize = d_bounding_box.numberCells();
      d_partition_dir = 0;
      for (tbox::Dimension::dir_t d = 1; d < d_dim.getValue(); ++d) {
         if (bbsize(d_partition_dir) < bbsize(d)) {
            d_partition_dir = d;
         }
      }

      int midpoint =
         (d_bounding_box.lower(d_partition_dir)
          + d_bounding_box.upper(d_partition_dir)) / 2;

      std::list<const Box *> left_boxes, right_boxes;
      for (std::list<const Box *>::iterator ni = d_boxes.begin();
           ni != d_boxes.end(); ) {
         const Box* box = *ni;
         if (box->upper(d_partition_dir) <= midpoint) {
            left_boxes.push_back(box);
            std::list<const Box *>::iterator curr = ni;
            ++ni;
            d_boxes.erase(curr);
         } else if (box->lower(d_partition_dir) > midpoint) {
            right_boxes.push_back(box);
            std::list<const Box *>::iterator curr = ni;
            ++ni;
            d_boxes.erase(curr);
         } else {
            ++ni;
         }
      }

      setupChildren(min_number, left_boxes, right_boxes);
   }

   if (s_max_lin_search[d_dim.getValue() - 1] <
       static_cast<unsigned int>(d_boxes.size())) {
      s_max_lin_search[d_dim.getValue() - 1] =
         static_cast<unsigned int>(d_boxes.size());
   }
}

/*
 **************************************************************************
 * This method finishes the tree generation by setting up the child
 * branches.  It expects the Boxes to have been split into
 * left_boxes, right_boxes, and d_boxes.  It will
 * generate the d_left_child and d_right_child.  If d_boxes is
 * big enough, it will generate d_center_child.
 **************************************************************************
 */
void
BoxTree::setupChildren(
   const int min_number,
   std::list<const Box *>& left_boxes,
   std::list<const Box *>& right_boxes)
{
   const int total_size = static_cast<int>(
         left_boxes.size()
         + right_boxes.size()
         + d_boxes.size());

   /*
    * If all Boxes are in a single child, the child is just as
    * big as its parent, so there is no point recursing.  Put
    * everything into d_boxes so the check below will prevent
    * recursion.
    */
   if (left_boxes.size() == static_cast<std::list<const Box *>::size_type>(total_size)) {
      left_boxes.swap(d_boxes);
   } else if (right_boxes.size() == static_cast<std::list<const Box *>::size_type>(total_size)) {
      right_boxes.swap(d_boxes);
   }

#if 0
   tbox::plog << "Split " << d_boxes.size() << "  " << d_bounding_box
              << " across " << d_partition_dir << " at " << mid << " into "
              << ' ' << left_boxes.size()
              << ' ' << cent_boxes.size()
              << ' ' << right_boxes.size()
              << std::endl;
#endif
   /*
    * If d_boxes is big enough, generate a center child for it.
    */
   if (d_boxes.size() >
       static_cast<std::list<const Box *>::size_type>(min_number) /* recursion criterion */ &&
       d_boxes.size() <
       static_cast<std::list<const Box *>::size_type>(total_size) /* avoid infinite recursion */) {
      d_center_child.reset(new BoxTree(d_dim));
      d_boxes.swap(d_center_child->d_boxes);
      d_center_child->privateGenerateTree(min_number);
      d_boxes.clear();   // No longer needed for tree construction or search.
   }

   /*
    * Recurse to build this node's left and right children.
    */
   if (!left_boxes.empty()) {
      d_left_child.reset(new BoxTree(d_dim));
      left_boxes.swap(d_left_child->d_boxes);
      d_left_child->privateGenerateTree(min_number);
   }
   if (!right_boxes.empty()) {
      d_right_child.reset(new BoxTree(d_dim));
      right_boxes.swap(d_right_child->d_boxes);
      d_right_child->privateGenerateTree(min_number);
   }
}

/*
 **************************************************************************
 * Returns true if any Box in the tree intersects the argument.
 **************************************************************************
 */

bool
BoxTree::hasOverlap(
   const Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   bool has_overlap = false;
   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         has_overlap = d_center_child->hasOverlap(box);
      } else {
         for (std::list<const Box *>::const_iterator ni = d_boxes.begin();
              ni != d_boxes.end(); ++ni) {
            if (box.intersects(**ni)) {
               has_overlap = true;
               break;
            }
         }
      }

      if (!has_overlap && d_left_child) {
         has_overlap = d_left_child->hasOverlap(box);
      }

      if (!has_overlap && d_right_child) {
         has_overlap = d_right_child->hasOverlap(box);
      }
   }
   return has_overlap;
}

/*
 **************************************************************************
 * Fills the vector with pointers to Boxes that intersect the arguement
 **************************************************************************
 */
void
BoxTree::findOverlapBoxes(
   std::vector<const Box *>& overlap_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_boxes.size());
#ifndef _OPENMP
      t_search[d_dim.getValue() - 1]->start();
#endif
   }

   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT(box.getBlockId() == d_block_id);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_boxes, box, true);
      } else {
         for (std::list<const Box *>::const_iterator ni = d_boxes.begin();
              ni != d_boxes.end(); ++ni) {
            const Box* my_box = *ni;
            if (box.intersects(*my_box)) {
               overlap_boxes.push_back(my_box);
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_boxes, box, true);
      }
   }

   if (!recursive_call) {
#ifndef _OPENMP
      t_search[d_dim.getValue() - 1]->stop();
#endif
      num_found_box = static_cast<int>(overlap_boxes.size())
         - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}

/*
 **************************************************************************
 * Fills the vector with Boxes that intersect the arguement
 **************************************************************************
 */
void
BoxTree::findOverlapBoxes(
   BoxContainer& overlap_boxes,
   const Box& box,
   bool recursive_call) const
{
   int num_found_box = 0;
   if (!recursive_call) {
      ++s_num_search[d_dim.getValue() - 1];
      num_found_box = static_cast<int>(overlap_boxes.size());
#ifndef _OPENMP
      t_search[d_dim.getValue() - 1]->start();
#endif
   }

   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   TBOX_ASSERT(box.getBlockId() == d_block_id);

   if (box.intersects(d_bounding_box)) {

      if (d_center_child) {
         d_center_child->findOverlapBoxes(overlap_boxes, box, true);
      } else {
         if (overlap_boxes.isOrdered()) {
            for (std::list<const Box *>::const_iterator ni = d_boxes.begin();
                 ni != d_boxes.end(); ++ni) {
               const Box* this_box = *ni;
               if (box.intersects(*this_box)) {
                  overlap_boxes.insert(*this_box);
               }
            }
         } else {
            for (std::list<const Box *>::const_iterator ni = d_boxes.begin();
                 ni != d_boxes.end(); ++ni) {
               const Box* this_box = *ni;
               if (box.intersects(*this_box)) {
                  overlap_boxes.pushBack(*this_box);
               }
            }
         }
      }

      if (d_left_child) {
         d_left_child->findOverlapBoxes(overlap_boxes, box, true);
      }

      if (d_right_child) {
         d_right_child->findOverlapBoxes(overlap_boxes, box, true);
      }
   }

   if (!recursive_call) {
#ifndef _OPENMP
      t_search[d_dim.getValue() - 1]->stop();
#endif
      num_found_box = static_cast<int>(overlap_boxes.size()) - num_found_box;
      s_max_found_box[d_dim.getValue() - 1] =
         tbox::MathUtilities<int>::Max(s_max_found_box[d_dim.getValue() - 1],
            num_found_box);
      s_num_found_box[d_dim.getValue() - 1] += num_found_box;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTree::initializeCallback()
{
   for (int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      const std::string dim_str(tbox::Utilities::intToString(i + 1));
      t_build_tree[i] = tbox::TimerManager::getManager()->
         getTimer(std::string("hier::BoxTree::build_tree[") + dim_str + "]");
      t_search[i] = tbox::TimerManager::getManager()->
         getTimer(std::string("hier::BoxTree::search[") + dim_str + "]");
   }
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */
void
BoxTree::finalizeCallback()
{
   for (int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      t_build_tree[i].reset();
      t_search[i].reset();
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void
BoxTree::resetStatistics(
   const tbox::Dimension& dim)
{
   s_num_build[dim.getValue() - 1] = 0;
   s_num_generate[dim.getValue() - 1] = 0;
   s_num_duplicate[dim.getValue() - 1] = 0;
   s_num_search[dim.getValue() - 1] = 0;
   s_num_sorted_box[dim.getValue() - 1] = 0;
   s_num_found_box[dim.getValue() - 1] = 0;
   s_max_sorted_box[dim.getValue() - 1] = 0;
   s_max_found_box[dim.getValue() - 1] = 0;
   s_max_lin_search[dim.getValue() - 1] = 0;
}

/*
 ***************************************************************************
 ***************************************************************************
 */
void
BoxTree::printStatistics(
   const tbox::Dimension& dim)
{
   tbox::plog << "BoxTree local stats:"
              << "  build=" << s_num_build[dim.getValue() - 1]
              << "  generate=" << s_num_generate[dim.getValue() - 1]
              << "  duplicate=" << s_num_duplicate[dim.getValue() - 1]
              << "  search=" << s_num_search[dim.getValue() - 1]
              << "  sorted_box=" << s_num_sorted_box[dim.getValue() - 1]
              << "  found_box=" << s_num_found_box[dim.getValue() - 1]
              << "  max_sorted_box=" << s_max_sorted_box[dim.getValue() - 1]
              << "  max_found_box=" << s_max_found_box[dim.getValue() - 1]
              << "  max_lin_search=" << s_max_lin_search[dim.getValue() - 1]
              << std::endl;

   tbox::Statistician* st = tbox::Statistician::getStatistician();
   std::shared_ptr<tbox::Statistic> bdstat(st->getStatistic("num_build",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> gnstat(st->getStatistic("num_generate",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> dpstat(st->getStatistic("num_duplicate",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> srstat(st->getStatistic("num_search",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> sbstat(st->getStatistic("num_sorted_box",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> fbstat(st->getStatistic("num_found_box",
                                                "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> msbstat(st->getStatistic("max_sorted_box",
                                                 "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> mfbstat(st->getStatistic("max_found_box",
                                                 "PROC_STAT"));
   std::shared_ptr<tbox::Statistic> lsstat(st->getStatistic("max_lin_search",
                                                "PROC_STAT"));

   static int seq_num = 0;
   bdstat->recordProcStat(s_num_build[dim.getValue() - 1], seq_num);
   gnstat->recordProcStat(s_num_generate[dim.getValue() - 1], seq_num);
   dpstat->recordProcStat(s_num_duplicate[dim.getValue() - 1], seq_num);
   srstat->recordProcStat(s_num_search[dim.getValue() - 1], seq_num);
   sbstat->recordProcStat(s_num_sorted_box[dim.getValue() - 1], seq_num);
   fbstat->recordProcStat(s_num_found_box[dim.getValue() - 1], seq_num);
   lsstat->recordProcStat(s_max_lin_search[dim.getValue() - 1], seq_num);

   st->finalize(false);
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int nproc = mpi.getSize();

   double avg;
   double min, max;
   int rmin(0), rmax(0);

   int doublewidth = 6;
   int intwidth = 6;
   int namewidth = 20;

   min = max = avg = s_num_build[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << bdstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_generate[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << gnstat->getName()
              << "  " << std::setw(namewidth) << std::setw(doublewidth)
              << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_duplicate[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << dpstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_search[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << srstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_sorted_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << sbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_num_found_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << fbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_sorted_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << msbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_found_box[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << mfbstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   min = max = avg = s_max_lin_search[dim.getValue() - 1];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&min, 1, MPI_MIN, &rmin);
      mpi.AllReduce(&max, 1, MPI_MAX, &rmax);
      mpi.AllReduce(&avg, 1, MPI_SUM);
      avg /= nproc;
   }
   tbox::plog << std::setw(namewidth) << bdstat->getName()
              << "  " << std::setw(doublewidth) << std::setprecision(0) << avg
              << " [ " << std::setw(doublewidth) << std::setprecision(0)
              << min << " at " << std::setw(intwidth) << rmin
              << " -> " << std::setw(doublewidth) << std::setprecision(0)
              << max << " at " << std::setw(intwidth) << rmax
              << " ]"
              << std::endl;

   ++seq_num;
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
