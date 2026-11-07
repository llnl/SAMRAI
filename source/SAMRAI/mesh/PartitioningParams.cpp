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

namespace SAMRAI {
namespace mesh {

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
   d_min_size(min_size),
   d_max_size(max_size),
   d_bad_interval(bad_interval, grid_geometry.getNumberBlocks()),
   d_cut_factor(cut_factor),
   d_minimum_cells(minimum_cells),
   d_artificial_minimum_load(artificial_minimum_load),
   d_flexible_load_tol(flexible_load_tol),
   d_load_comparison_tol(1e-6),
   d_using_vouchers(false),
   d_work_data_id(-1)
{
   for (hier::BlockId::block_t bid(0); bid < grid_geometry.getNumberBlocks(); ++bid) {
      grid_geometry.computePhysicalDomain(
         d_block_domain_boxes[hier::BlockId(bid)], ratio_to_level_zero, hier::BlockId(bid));
   }
}

PartitioningParams::PartitioningParams(
   const PartitioningParams& other):
   d_block_domain_boxes(other.d_block_domain_boxes),
   d_min_size(other.d_min_size),
   d_max_size(other.d_max_size),
   d_bad_interval(other.d_bad_interval),
   d_cut_factor(other.d_cut_factor),
   d_minimum_cells(other.d_minimum_cells),
   d_artificial_minimum_load(other.d_artificial_minimum_load),
   d_load_comparison_tol(other.d_load_comparison_tol),
   d_using_vouchers(other.d_using_vouchers),
   d_work_data_id(other.d_work_data_id)
{
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
