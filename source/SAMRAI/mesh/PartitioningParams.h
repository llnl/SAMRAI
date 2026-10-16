/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parameters in load balancing.
 *
 ************************************************************************/

#ifndef included_mesh_PartitioningParams
#define included_mesh_PartitioningParams

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Box.h"

#include <map>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Light weight class holding parameters generally used
 * in partitioning.
 */

class PartitioningParams
{
public:
   /*!
    * @brief Resolved load model for one partitioning operation.
    *
    * CascadePartitioner resolves hierarchy-level configuration and ghost
    * widths before constructing this value.  The resulting model is immutable
    * and valid for every box accepted by the corresponding partitioning
    * parameters.
    */
   struct LoadModel
   {
   public:
      enum Type {
         CELL_COUNT,
         LINEAR
      };

      static LoadModel
      cellCount(
         const tbox::Dimension& dim,
         double artificial_minimum);

      static LoadModel
      linear(
         double slope,
         double intercept,
         const hier::IntVector& ghost_width,
         const hier::IntVector& minimum_box_size);

      double
      computeBoxLoad(
         const hier::Box& box) const;

      double
      computeMinimumBoxLoad(
         const hier::IntVector& minimum_box_size) const;

      double
      computeSplitWeight(
         const hier::Box& box) const;

      bool isLinear() const {
         return d_type == LINEAR;
      }

      double getSlope() const {
         return d_slope;
      }

      double getIntercept() const {
         return d_intercept;
      }

      double getArtificialMinimum() const {
         return d_artificial_minimum;
      }

      const hier::IntVector& getGhostWidth() const {
         return d_ghost_width;
      }

   private:
      LoadModel(
         Type type,
         double slope,
         double intercept,
         double artificial_minimum,
         const hier::IntVector& ghost_width);

      double
      computeLoad(
         double modeled_cell_count) const;

      Type d_type;
      double d_slope;
      double d_intercept;
      double d_artificial_minimum;
      hier::IntVector d_ghost_width;
   };

   PartitioningParams(
      const hier::BaseGridGeometry& grid_geometry,
      const hier::IntVector& ratio_to_level_zero,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      size_t minimum_cells,
      double artificial_minimum_load,
      double flexible_load_tol);

   PartitioningParams(
      const hier::BaseGridGeometry& grid_geometry,
      const hier::IntVector& ratio_to_level_zero,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      size_t minimum_cells,
      const LoadModel& load_model,
      double flexible_load_tol);

   PartitioningParams(
      const PartitioningParams& other);

   double getMinBoxSizeProduct() const {
      return static_cast<double>(d_min_size.getProduct());
   }

   double getMinBoxLoad() const {
      return d_load_model.computeMinimumBoxLoad(d_min_size);
   }

   const hier::IntVector& getMinBoxSize() const {
      return d_min_size;
   }

   const hier::IntVector& getMaxBoxSize() const {
      return d_max_size;
   }

   const hier::BoxContainer& getDomainBoxes(const hier::BlockId& bid) const {
      return d_block_domain_boxes.find(bid)->second;
   }

   const hier::IntVector& getBadInterval() const {
      return d_bad_interval;
   }

   const hier::IntVector& getCutFactor() const {
      return d_cut_factor;
   }

   size_t getMinimumCellRequest() const {
      return d_minimum_cells;
   } 

   double getArtificialMinimumLoad() const {
      return d_load_model.getArtificialMinimum();
   }

   const tbox::Dimension& getDim() const {
      return d_min_size.getDim();
   }

   const double& getFlexibleLoadTol() const {
      return d_flexible_load_tol;
   }

   const double& getLoadComparisonTol() const {
      return d_load_comparison_tol;
   }

   bool usingLinearLoad() const {
      return d_load_model.isLinear();
   }

   /*!
    * @brief Compute the configured load for a box.
    *
    * Linear loading uses the volume after growing by the model's ghost width.
    * Cell-count loading applies the configured artificial minimum.
    */
   double
   computeBoxLoad(
      const hier::Box& box) const;

   /*!
    * @brief Compute a box's relative weight when apportioning split load.
    *
    * Unlike computeBoxLoad(), this does not apply the artificial minimum in
    * cell-count mode because split loads are apportioned by cell count.
    */
   double
   computeSplitWeight(
      const hier::Box& box) const;

   const LoadModel& getLoadModel() const {
      return d_load_model;
   }

   void
   setLoadModel(
      const LoadModel& load_model);

   const bool& usingVouchers() const {
      return d_using_vouchers;
   }

   void setUsingVouchers(bool using_vouchers) {
      d_using_vouchers = using_vouchers;
   }

   const int& getWorkloadDataId() const {
      return d_work_data_id;
   }

   void setWorkloadDataId(int work_data_id) {
      TBOX_ASSERT(work_data_id >= 0);
      d_work_data_id = work_data_id;
   }

   const hier::PatchLevel& getWorkloadPatchLevel() const {
      return *d_workload_level;
   }

   void setWorkloadPatchLevel(std::shared_ptr<hier::PatchLevel>& level) {
      TBOX_ASSERT(level.get());
      d_workload_level = level;
   }


   friend std::ostream&
   operator << (
      std::ostream& os,
      const PartitioningParams& pp);

private:
   std::map<hier::BlockId, hier::BoxContainer> d_block_domain_boxes;
   hier::IntVector d_min_size;
   hier::IntVector d_max_size;
   hier::IntVector d_bad_interval;
   hier::IntVector d_cut_factor;

   /*
    * @brief The requested minimum for number of cells in a patch
    */
   size_t d_minimum_cells;

   /*! @brief Validated load model for this partitioning operation. */
   LoadModel d_load_model;

   /*!
    * @brief Fraction of ideal load a process can accept over and
    * above the ideal.
    */
   double d_flexible_load_tol;

   /*!
    * @brief Tolerance for comparing floating point loads.
    *
    * Should be set to at least possible rounding errors.
    * Better if set between that and the greatest work value
    * that would be considered "no work".
    */
   double d_load_comparison_tol;

   /*!
    * @brief Flag for using voucher method or not.
    */
   bool d_using_vouchers;

   /*!
    * @brief Patch data id for nonuniform workload
    */
   int d_work_data_id;

   /*!
    * @brief Pointer to level holding nonuniform workload
    */
   std::shared_ptr<hier::PatchLevel> d_workload_level;
};

}
}

#endif
