/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic integration scheme.
 *
 ************************************************************************/
#include "SAMRAI/algs/HyperbolicPatchStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace algs {

HyperbolicPatchStrategy::HyperbolicPatchStrategy():
   xfer::RefinePatchStrategy(),
   xfer::CoarsenPatchStrategy(),
   d_data_context()
{
}

HyperbolicPatchStrategy::~HyperbolicPatchStrategy()
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
HyperbolicPatchStrategy::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("HyperbolicPatchStrategy::tagGradientDetectorCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
HyperbolicPatchStrategy::tagRichardsonExtrapolationCells(
   hier::Patch& patch,
   const int error_level_number,
   const std::shared_ptr<hier::VariableContext>& coarsened_fine,
   const std::shared_ptr<hier::VariableContext>& advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   NULL_USE(patch);
   NULL_USE(error_level_number);
   NULL_USE(coarsened_fine);
   NULL_USE(advanced_coarse);
   NULL_USE(regrid_time);
   NULL_USE(deltat);
   NULL_USE(error_coarsen_ratio);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_gradient_detector_too);
   TBOX_ERROR("HyperbolicPatchStrategy::tagRichardsonExtrapolationCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
HyperbolicPatchStrategy::setupLoadBalancer(
   HyperbolicLevelIntegrator* integrator,
   mesh::GriddingAlgorithm* gridding_algorithm)
{
   NULL_USE(integrator);
   NULL_USE(gridding_algorithm);
}

void
HyperbolicPatchStrategy::preprocessAdvanceLevelState(
   const std::shared_ptr<hier::PatchLevel>& level,
   double current_time,
   double dt,
   bool first_step,
   bool last_step,
   bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(current_time);
   NULL_USE(dt);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
}

void
HyperbolicPatchStrategy::postprocessAdvanceLevelState(
   const std::shared_ptr<hier::PatchLevel>& level,
   double current_time,
   double dt,
   bool first_step,
   bool last_step,
   bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(current_time);
   NULL_USE(dt);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
}

}
}
