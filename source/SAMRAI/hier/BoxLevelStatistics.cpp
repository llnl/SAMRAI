/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a BoxLevel.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include "SAMRAI/tbox/MathUtilities.h"

#include <cmath>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

std::string BoxLevelStatistics::s_quantity_names[NUMBER_OF_QUANTITIES];
int BoxLevelStatistics::s_longest_length;

tbox::StartupShutdownManager::Handler
BoxLevelStatistics::s_initialize_finalize_handler(
   BoxLevelStatistics::initializeCallback,
   0,
   0,
   BoxLevelStatistics::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ************************************************************************
 * Constructor.
 ************************************************************************
 */
BoxLevelStatistics::BoxLevelStatistics(
   const BoxLevel& box_level):
   d_mpi(box_level.getMPI()),
   d_dim(box_level.getDim())
{
   if (!box_level.isInitialized()) {
      TBOX_ERROR("BoxLevelStatistics requires an initialized BoxLevel.");
   }

   computeLocalBoxLevelStatistics(box_level);
   reduceStatistics();
}

/*
 ************************************************************************
 ************************************************************************
 */
BoxLevelStatistics::StatisticalQuantities::StatisticalQuantities()
{
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      d_values[i] = 0;
   }
}

/*
 ***********************************************************************
 * Compute statistics for local process.
 ***********************************************************************
 */

void BoxLevelStatistics::computeLocalBoxLevelStatistics(const BoxLevel& box_level)
{
   box_level.cacheGlobalReducedData();

   /*
    * Compute per-processor statistics.  Some quantities are readily
    * available while others are computed in the loop following.
    *
    * Aspect ratio uses a generalized formula that goes to 1 when box
    * has same length on all sides (regardless of dimension),
    * degenerates to the rectangular aspect ratio in 2D, and grows
    * appropriately for dimensions higher than 2.
    */

   d_sq.d_values[HAS_ANY_BOX] = (box_level.getLocalNumberOfBoxes() > 0);
   d_sq.d_values[NUMBER_OF_CELLS] =
      static_cast<double>(box_level.getLocalNumberOfCells());
   d_sq.d_values[NUMBER_OF_BOXES] =
      static_cast<double>(box_level.getLocalNumberOfBoxes());
   d_sq.d_values[MAX_BOX_VOL] = 0;
   d_sq.d_values[MIN_BOX_VOL] = tbox::MathUtilities<double>::getMax();
   d_sq.d_values[MAX_BOX_LEN] = 0;
   d_sq.d_values[MIN_BOX_LEN] = tbox::MathUtilities<double>::getMax();
   d_sq.d_values[MAX_ASPECT_RATIO] = 0;
   d_sq.d_values[SUM_ASPECT_RATIO] = 0;
   d_sq.d_values[SUM_SURFACE_AREA] = 0.;
   d_sq.d_values[SUM_NORM_SURFACE_AREA] = 0.;

   const BoxContainer& boxes = box_level.getBoxes();

   for (RealBoxConstIterator ni(boxes.realBegin());
        ni != boxes.realEnd(); ++ni) {

      const Box& box = *ni;
      const IntVector boxdims = box.numberCells();
      const double boxvol = static_cast<double>(boxdims.getProduct());
      const int longdim = boxdims.max();
      const int shortdim = boxdims.min();
      double aspect_ratio = 0.0;
      double surfarea = 0.;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         surfarea += 2 * double(boxvol) / boxdims(d);
         double tmp = static_cast<double>(boxdims(d)) / shortdim - 1.0;
         aspect_ratio += tmp * tmp;
      }
      aspect_ratio = 1.0 + sqrt(aspect_ratio);

      d_sq.d_values[MAX_BOX_VOL] =
         tbox::MathUtilities<double>::Max(d_sq.d_values[MAX_BOX_VOL],
            boxvol);
      d_sq.d_values[MIN_BOX_VOL] =
         tbox::MathUtilities<double>::Min(d_sq.d_values[MIN_BOX_VOL],
            boxvol);

      d_sq.d_values[MAX_BOX_LEN] =
         tbox::MathUtilities<double>::Max(d_sq.d_values[MAX_BOX_LEN],
            longdim);
      d_sq.d_values[MIN_BOX_LEN] =
         tbox::MathUtilities<double>::Min(d_sq.d_values[MIN_BOX_LEN],
            shortdim);

      d_sq.d_values[MAX_ASPECT_RATIO] =
         tbox::MathUtilities<double>::Max(d_sq.d_values[MAX_ASPECT_RATIO],
            aspect_ratio);

      d_sq.d_values[SUM_ASPECT_RATIO] += aspect_ratio;
      d_sq.d_values[SUM_SURFACE_AREA] += surfarea;

   }

   /*
    * Smallest surface area possible for the number of cells perfectly
    * distributed in d_mpi.
    */
   const double ideal_surfarea =
      2 * d_dim.getValue()
      * pow(double(box_level.getGlobalNumberOfCells()) / d_mpi.getSize(),
         double(d_dim.getValue() - 1) / d_dim.getValue());

   d_sq.d_values[SUM_NORM_SURFACE_AREA] =
      d_sq.d_values[SUM_SURFACE_AREA] / ideal_surfarea;

}

/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevelStatistics::reduceStatistics()
{
   d_sq_min = d_sq;
   d_sq_max = d_sq;
   d_sq_sum = d_sq;

   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce(
         d_sq_min.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MINLOC,
         d_rank_of_min);
      d_mpi.AllReduce(
         d_sq_max.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MAXLOC,
         d_rank_of_max);
      d_mpi.AllReduce(d_sq_sum.d_values, NUMBER_OF_QUANTITIES, MPI_SUM);
   } else {
      for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
         d_rank_of_min[i] = d_rank_of_max[i] = 0;
      }
   }
}

/*
 ***********************************************************************
 * Write out local and globally reduced statistics on the boxes.
 ***********************************************************************
 */

void BoxLevelStatistics::printBoxStats(
   std::ostream& co,
   const std::string& border) const
{
   co.unsetf(std::ios::fixed | std::ios::scientific);
   co.precision(3);

   /*
    * Smallest surface area possible for the number of cells perfectly
    * distributed in d_mpi.
    */
   const double ideal_width =
      pow(d_sq_sum.d_values[NUMBER_OF_CELLS] / d_mpi.getSize(),
         1.0 / d_dim.getValue());
   const double ideal_surfarea = 2 * d_dim.getValue()
      * pow(ideal_width, double(d_dim.getValue() - 1));

   co << border << "N = " << d_sq_sum.d_values[NUMBER_OF_BOXES]
      << " (global number of boxes)\n"
      << border << "P = " << d_mpi.getSize() << " (number of processes)\n"
      << border << "Ideal width (W) is " << ideal_width
      << ", surface area (A) is " << ideal_surfarea << " for "
      << (d_sq_sum.d_values[NUMBER_OF_CELLS] / d_mpi.getSize()) << " cells\n"
      << border << std::setw(s_longest_length) << std::string()
      << "    local        min               max             sum    sum/N    sum/P\n";

   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      co << border << std::setw(s_longest_length) << std::left
         << s_quantity_names[i]
         << ' ' << std::setw(8) << std::right << d_sq.d_values[i]
         << ' ' << std::setw(8) << std::right << d_sq_min.d_values[i] << " @ "
         << std::setw(6) << std::left << d_rank_of_min[i]
         << ' ' << std::setw(8) << std::right << d_sq_max.d_values[i] << " @ "
         << std::setw(6) << std::left << d_rank_of_max[i]
         << ' ' << std::setw(8) << std::right << d_sq_sum.d_values[i]
         << ' ' << std::setw(8)
         << std::right << d_sq_sum.d_values[i] / d_sq_sum.d_values[NUMBER_OF_BOXES]
         << ' ' << std::setw(8)
         << std::right << d_sq_sum.d_values[i] / d_mpi.getSize() << '\n';
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevelStatistics::initializeCallback()
{
   s_quantity_names[HAS_ANY_BOX] = "has any box";
   s_quantity_names[NUMBER_OF_CELLS] = "num cells";
   s_quantity_names[NUMBER_OF_BOXES] = "num boxes";
   s_quantity_names[MAX_BOX_VOL] = "max box vol";
   s_quantity_names[MIN_BOX_VOL] = "min box vol";
   s_quantity_names[MAX_BOX_LEN] = "max box len";
   s_quantity_names[MIN_BOX_LEN] = "min box len";
   s_quantity_names[MAX_ASPECT_RATIO] = "max aspect ratio";
   s_quantity_names[SUM_ASPECT_RATIO] = "sum aspect ratio";
   s_quantity_names[SUM_SURFACE_AREA] = "sum surf area";
   s_quantity_names[SUM_NORM_SURFACE_AREA] = "sum surf area/A";
   s_longest_length = 0;
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      s_longest_length = tbox::MathUtilities<int>::Max(
            s_longest_length, static_cast<int>(s_quantity_names[i].length()));
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */

void
BoxLevelStatistics::finalizeCallback()
{
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      s_quantity_names[i].clear();
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
