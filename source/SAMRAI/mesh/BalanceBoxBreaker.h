/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utilitiy for breaking boxes during partitioning.
 *
 ************************************************************************/

#ifndef included_mesh_BalanceBoxBreaker
#define included_mesh_BalanceBoxBreaker

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/mesh/PartitioningParams.h"

#include <vector>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Utilities for breaking up boxes during partitioning.
 */

class BalanceBoxBreaker
{

public:
   /*!
    * @brief Constructor.
    *
    * @param[in] pparams
    *
    * @param[in] print_break_steps Flag for debugging this class.
    */
   BalanceBoxBreaker(
      const PartitioningParams& pparams,
      bool print_break_steps = false);

   //! @brief Copy constructor.
   BalanceBoxBreaker(
      const BalanceBoxBreaker& other);

   /*!
    * @brief Break off a given load size from a given Box.
    *
    * Attempt to break off the ideal_load, or at least a load
    * inside the range [low_load, high_load].
    *
    * @param[out] breakoff Boxes broken off (usually just one).
    *
    * @param[out] leftover Remainder of box after breakoff is gone.
    *
    * @param[out] brk_load The load broken off.
    *
    * @param[in] box Box to break.
    *
    * @param[in] corner_weights Fraction of the total load associated with
    *                           each corner.
    *
    * @param[in] box_load Load of the box before breaking.
    *
    * @param[in] ideal_load Ideal load to break.
    *
    * @param[in] low_load
    *
    * @param[in] high_load
    *
    * @param[in] threshold_width Try to avoid making boxes thinner
    * than this width in any direction.
    *
    * @return whether a successful break was made.
    *
    * @pre ideal_load_to_break > 0
    */
   bool
   breakOffLoad(
      hier::BoxContainer& breakoff,
      hier::BoxContainer& leftover,
      double& brk_load,
      const hier::Box& box,
      double box_load,
      const std::vector<double>& corner_weights,
      double ideal_load,
      double low_load,
      double high_load,
      double threshold_width) const;

   /*!
    * @brief Set whether to print box breaking steps for debugging.
    */
   void setPrintBreakSteps(bool print_break_steps) {
      d_print_break_steps = print_break_steps;
   }

   /*!
    * @brief Break up box bursty against box solid and adds the pieces
    * in the given container.
    *
    * This version differs from that in BoxContainer in that it tries
    * to minimize slivers.
    */
   static void
   burstBox(
      hier::BoxContainer& boxes,
      const hier::Box& bursty,
      const hier::Box& solid);

   /*!
    * @brief Compute the penalty associated with an imbalance.
    */
   static double
   computeBalancePenalty(double imbalance)
   {
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   /*!
    * @brief Compute a score that is low for box widths smaller than
    * some threshold_width.
    */
   static double
   computeWidthScore(
      const hier::IntVector& box_size,
      double threshold_width);

   /*!
    * @brief Compute a combined width score for multiple boxes.
    */
   static double
   computeWidthScore(
      const hier::BoxContainer& boxes,
      double threshold_width);

private:
   //! @brief Information regarding a potential way to break a box.
   struct TrialBreak {
      TrialBreak(
         const PartitioningParams& pparams,
         double threshold_width,
         const hier::Box& whole_box,
         double whole_box_load,
         const std::vector<double>& corner_weights,
         const std::vector<std::vector<bool> >& bad_cuts,
         double ideal_load,
         double low_load,
         double high_load);
      //! @brief Construct TrialBreak where breakoff and leftover are reversed.
      TrialBreak(
         const TrialBreak& orig,
         bool make_reverse);
      //! @brief Compute data for breaking box from whole and store results.
      void
      computeBreakData(
         const hier::Box& box);

      //! @brief Compute load that would break off if the box is broken off.
      double computeBreakOffLoad(
         const hier::Box& box);

      void
      swapWithReversedTrial(
         TrialBreak& reversed);
      //! @brief Compute merits vs doing nothing and return improvement flag.
      bool
      computeMerits();
      //! @brief Whether this improves over another (or degrades or leaves alone).
      int
      improvesOver(
         const TrialBreak& other) const;

      //! @brief Swap this object with another.
      void
      swap(
         TrialBreak& other);

      double d_breakoff_load;
      hier::BoxContainer d_breakoff;
      hier::BoxContainer d_leftover;
      const double d_ideal_load;
      //! Part of the target break-off range [d_low_load,d_high_load]
      const double d_low_load;
      //! Part of the target break-off range [d_low_load,d_high_load]
      const double d_high_load;
      double d_width_score;
      double d_balance_penalty;
      //! @brief Flags from comparing this trial vs doing nothing.
      int d_flags[4];
      const PartitioningParams* d_pparams;
      const double d_threshold_width;
      const hier::Box& d_whole_box;
      double d_whole_box_load;
      const std::vector<double>& d_corner_weights;
      const std::vector<std::vector<bool> >& d_bad_cuts;
      std::vector<hier::Box> d_corner_box;
   };

private:
   bool
   breakOffLoad_planar(
      TrialBreak& trial) const;

   bool
   breakOffLoad_cubic(
      TrialBreak& trial) const;

   void
   setTimers();

   const PartitioningParams * d_pparams;

   //@{
   //! @name Debugging and diagnostic data

   bool d_print_break_steps;
   std::shared_ptr<tbox::Timer> t_break_off_load;
   std::shared_ptr<tbox::Timer> t_find_bad_cuts;

   //@}

public:
   friend std::ostream&
   operator << (
      std::ostream& os,
      const TrialBreak& tb);

};

}
}
#endif
