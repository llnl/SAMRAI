/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parameters in load balancing.
 *
 ************************************************************************/

#ifndef included_mesh_PartitioningParams_C
#define included_mesh_PartitioningParams_C

#include "SAMRAI/mesh/PartitioningParams.h"

#include "SAMRAI/tbox/MathUtilities.h"

#include <cmath>

namespace SAMRAI {
namespace mesh {

PartitioningParams::LoadModel::LoadModel(
   Type type,
   double slope,
   double intercept,
   double artificial_minimum,
   const hier::IntVector& ghost_width):
   d_type(type),
   d_slope(slope),
   d_intercept(intercept),
   d_artificial_minimum(artificial_minimum),
   d_ghost_width(ghost_width)
{
}

PartitioningParams::LoadModel
PartitioningParams::LoadModel::cellCount(
   const tbox::Dimension& dim,
   double artificial_minimum)
{
   if (!std::isfinite(artificial_minimum) || artificial_minimum < 0.0) {
      TBOX_ERROR(
         "PartitioningParams::LoadModel artificial minimum must be "
         << "non-negative and finite.\n");
   }

   return LoadModel(
      CELL_COUNT,
      1.0,
      0.0,
      artificial_minimum,
      hier::IntVector(dim, 0));
}

PartitioningParams::LoadModel
PartitioningParams::LoadModel::linear(
   double slope,
   double intercept,
   const hier::IntVector& ghost_width,
   const hier::IntVector& minimum_box_size)
{
   if (ghost_width.getDim() != minimum_box_size.getDim()) {
      TBOX_ERROR(
         "PartitioningParams::LoadModel ghost width and minimum box size "
         << "must have the same dimension.\n");
   }
   if (!std::isfinite(slope) || !std::isfinite(intercept) || slope < 0.0) {
      TBOX_ERROR(
         "PartitioningParams::LoadModel linear coefficients must be finite "
         << "and the slope must be non-negative.\n");
   }

   for (int d = 0; d < ghost_width.getDim().getValue(); ++d) {
      if (ghost_width[d] < 0) {
         TBOX_ERROR(
            "PartitioningParams::LoadModel ghost widths must be "
            << "non-negative.\n");
      }
   }

   LoadModel model(LINEAR, slope, intercept, 0.0, ghost_width);
   for (hier::BlockId::block_t b = 0;
        b < minimum_box_size.getNumBlocks(); ++b) {
      double minimum_grown_volume = 1.0;
      for (int d = 0; d < minimum_box_size.getDim().getValue(); ++d) {
         minimum_grown_volume *=
            static_cast<double>(minimum_box_size(b, d)) +
            2.0 * static_cast<double>(ghost_width[d]);
      }
      const double minimum_load =
         slope * minimum_grown_volume + intercept;
      if (!std::isfinite(minimum_load) || minimum_load <= 0.0) {
         TBOX_ERROR(
            "PartitioningParams::LoadModel must produce a positive, finite "
            << "load for every minimum-size ghost-grown box.\n");
      }
   }
   return model;
}

double
PartitioningParams::LoadModel::computeBoxLoad(
   const hier::Box& box) const
{
   double load;
   if (d_type == LINEAR) {
      hier::Box grown_box(box);
      grown_box.grow(d_ghost_width);
      load = d_slope * static_cast<double>(grown_box.size()) + d_intercept;
   } else {
      load = tbox::MathUtilities<double>::Max(
         static_cast<double>(box.size()),
         d_artificial_minimum);
   }

   if (!std::isfinite(load) || load <= 0.0) {
      TBOX_ERROR(
         "PartitioningParams::LoadModel produced a non-positive or "
         << "non-finite box load.\n");
   }
   return load;
}

double
PartitioningParams::LoadModel::computeSplitWeight(
   const hier::Box& box) const
{
   return d_type == LINEAR ?
      computeBoxLoad(box) : static_cast<double>(box.size());
}

PartitioningParams::PartitioningParams(
   const hier::BaseGridGeometry& grid_geometry,
   const hier::IntVector& ratio_to_level_zero,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   size_t minimum_cells,
   double artificial_minimum_load,
   double flexible_load_tol):
   PartitioningParams(
      grid_geometry,
      ratio_to_level_zero,
      min_size,
      max_size,
      bad_interval,
      cut_factor,
      minimum_cells,
      LoadModel::cellCount(min_size.getDim(), artificial_minimum_load),
      flexible_load_tol)
{
}

PartitioningParams::PartitioningParams(
   const hier::BaseGridGeometry& grid_geometry,
   const hier::IntVector& ratio_to_level_zero,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   size_t minimum_cells,
   const LoadModel& load_model,
   double flexible_load_tol):
   d_min_size(min_size),
   d_max_size(max_size),
   d_bad_interval(bad_interval, grid_geometry.getNumberBlocks()),
   d_cut_factor(cut_factor),
   d_minimum_cells(minimum_cells),
   d_load_model(load_model),
   d_flexible_load_tol(flexible_load_tol),
   d_load_comparison_tol(1e-6),
   d_using_vouchers(false),
   d_work_data_id(-1)
{
   if (d_load_model.getGhostWidth().getDim() != min_size.getDim()) {
      TBOX_ERROR(
         "PartitioningParams load model has the wrong dimension.\n");
   }
   for (hier::BlockId::block_t bid(0); bid < grid_geometry.getNumberBlocks(); ++bid) {
      grid_geometry.computePhysicalDomain(
         d_block_domain_boxes[hier::BlockId(bid)], ratio_to_level_zero, hier::BlockId(bid));
   }
}

PartitioningParams::PartitioningParams(
   const PartitioningParams& other) = default;

double
PartitioningParams::computeBoxLoad(
   const hier::Box& box) const
{
   return d_load_model.computeBoxLoad(box);
}

double
PartitioningParams::computeSplitWeight(
   const hier::Box& box) const
{
   return d_load_model.computeSplitWeight(box);
}

void
PartitioningParams::setLoadModel(
   const LoadModel& load_model)
{
   if (load_model.getGhostWidth().getDim() != d_min_size.getDim()) {
      TBOX_ERROR(
         "PartitioningParams load model has the wrong dimension.\n");
   }
   d_load_model = load_model;
}

std::ostream& operator << (
   std::ostream& os,
   const PartitioningParams& pp)
{
   os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
   os.precision(6);
   os << "min_size=" << pp.d_min_size
   << "  max_size=" << pp.d_max_size
   << "  bad_interval=" << pp.d_bad_interval
   << "  cut_factor=" << pp.d_cut_factor
   << "  flexible_load_tol=" << pp.d_flexible_load_tol
   << "  using_linear_load=" << pp.d_load_model.isLinear()
   << "  load_slope=" << pp.d_load_model.getSlope()
   << "  load_intercept=" << pp.d_load_model.getIntercept()
   << "  load_comparison_tol=" << pp.d_load_comparison_tol
   << "  work_data_id=" << pp.d_work_data_id;
   for (std::map<hier::BlockId, hier::BoxContainer>::const_iterator mi =
           pp.d_block_domain_boxes.begin();
        mi != pp.d_block_domain_boxes.end(); ++mi) {
      os << ' ' << mi->first << ':' << mi->second.format();
   }
   return os;
}

}
}

#endif
