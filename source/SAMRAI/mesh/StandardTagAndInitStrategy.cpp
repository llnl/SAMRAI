/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for Richardson extrapolation error detection.
 *
 ************************************************************************/
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace mesh {

StandardTagAndInitStrategy::StandardTagAndInitStrategy()
{
}

StandardTagAndInitStrategy::~StandardTagAndInitStrategy()
{
}

/*
 *************************************************************************
 *
 * Default virtual function implementations.
 *
 *************************************************************************
 */

void
StandardTagAndInitStrategy::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(error_data_time);
   NULL_USE(tag_index);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("StandardTagAndInitStrategy::applyGradientDetector()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
StandardTagAndInitStrategy::coarsenDataForRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::PatchLevel>& coarser_level,
   const double coarsen_data_time,
   const bool before_advance)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(coarser_level);
   NULL_USE(coarsen_data_time);
   NULL_USE(before_advance);
   TBOX_ERROR("StandardTagAndInitStrategy::coarsenDataForRichardsonExtrapolation()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
StandardTagAndInitStrategy::applyRichardsonExtrapolation(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double error_data_time,
   const int tag_index,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_time,
   const bool uses_gradient_detector_too)
{
   NULL_USE(level);
   NULL_USE(error_data_time);
   NULL_USE(tag_index);
   NULL_USE(deltat);
   NULL_USE(error_coarsen_ratio);
   NULL_USE(initial_time);
   NULL_USE(uses_gradient_detector_too);
   TBOX_ERROR("StandardTagAndInitStrategy::applyRichardsonExtrapolation()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

double
StandardTagAndInitStrategy::getLevelDt(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double dt_time,
   const bool initial_time)
{
   NULL_USE(level);
   NULL_USE(dt_time);
   NULL_USE(initial_time);
   TBOX_ERROR("StandardTagAndInitStrategy::getLevelDt()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
   return 0.0;
}

void
StandardTagAndInitStrategy::resetTimeDependentData(
   const std::shared_ptr<hier::PatchLevel>& level,
   const double new_time,
   const bool can_be_refined)
{
   NULL_USE(level);
   NULL_USE(new_time);
   NULL_USE(can_be_refined);
   TBOX_ERROR("StandardTagAndInitStrategy::resetTimeDependentData()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

double
StandardTagAndInitStrategy::advanceLevel(
   const std::shared_ptr<hier::PatchLevel>& level,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const double current_time,
   const double new_time,
   const bool first_step,
   const bool last_step,
   const bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(hierarchy);
   NULL_USE(current_time);
   NULL_USE(new_time);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
   TBOX_ERROR("StandardTagAndInitStrategy::advanceLevel()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
   return 0.0;
}

void
StandardTagAndInitStrategy::resetDataToPreadvanceState(
   const std::shared_ptr<hier::PatchLevel>& level)
{
   NULL_USE(level);
   TBOX_ERROR("StandardTagAndInitStrategy::resetDataToPreadvanceState()"
      << "\nNo class derived supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
StandardTagAndInitStrategy::processHierarchyBeforeAddingNewLevel(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::BoxLevel>& new_box_level)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(new_box_level);
}

void
StandardTagAndInitStrategy::processLevelBeforeRemoval(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const std::shared_ptr<hier::PatchLevel>& old_level)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(old_level);
}

}
}
