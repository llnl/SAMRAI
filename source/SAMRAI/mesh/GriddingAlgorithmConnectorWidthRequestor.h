/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   GriddingAlgorihtm's implementation of PatchHierarchy
 *
 ************************************************************************/

#ifndef included_mesh_GriddingAlgorithmConnectorWidthRequestor
#define included_mesh_GriddingAlgorithmConnectorWidthRequestor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide GriddingAlgorithm needs Connectors
 * between hierarchy levels to be.
 */
class GriddingAlgorithmConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   GriddingAlgorithmConnectorWidthRequestor();

   /*!
    * @brief Compute Connector widths that this class requires in
    * order to work properly on a given hierarchy.
    *
    * Implements the pure virtual method
    * hier::PatchHierarchy::ConnectorWidthRequestorStrategy::computeRequiredConnectorWidths()
    *
    * @param[out] self_connector_widths Array of widths for Connectors
    * from a level to itself.
    *
    * @param[out] fine_connector_widths Array of widths for Connectors
    * from a level to the next finer level.
    *
    * @param[in]  patch_hierarchy
    */
   void
   computeRequiredConnectorWidths(
      std::vector<hier::IntVector>& self_connector_widths,
      std::vector<hier::IntVector>& fine_connector_widths,
      const hier::PatchHierarchy& patch_hierarchy) const;

   /*!
    * @brief Plan to specify enough Connector width to support the
    * given width between tag and cluster levels.
    *
    * This is only used by the GriddingAlgorithm that owns this object.
    */
   void
   setTagToClusterWidth(
      std::vector<hier::IntVector>& tag_to_cluster_width);

private:
   /*!
    * @brief Compute the Connector widths needed at a given level
    * number to support generating given Connector width at finer
    * levels.
    *
    * @param[out] coarse_to_fine_width Required width for the
    * coarse--->fine Connector.
    *
    * @param[out] coarse_to_coarse_width Required width for the
    * coarse--->coarse Connector.
    *
    * @param[in] fine_to_fine_width Fine--->fine width which the
    * coarse--->fine Connector must support.
    *
    * @param[in] fine_to_coarse_ratio The refinement ratio.
    *
    * @param[in] nesting_buffer_at_fine The nesting buffer between
    * fine and next finer level, measured in the fine index space.
    *
    * @param[in] max_stencil_width_at_coarse The maximum stencil width
    * to be used on the coarse level.
    *
    * @param[in] max_ghost_width_at_coarse The maximum ghost cell
    * width to be used on the coarse level.
    *
    * @pre (dim == coarse_to_fine_width.getDim()) &&
    *      (dim == coarse_to_coarse_width.getDim()) &&
    *      (dim == fine_to_fine_width.getDim()) &&
    *      (dim == fine_to_coarse_ratio.getDim()) &&
    *      (dim == nesting_buffer_at_fine.getDim()) &&
    *      (dim == max_stencil_width_at_coarse.getDim()) &&
    *      (dim == max_ghost_width_at_coarse.getDim())
    */
   void
   computeCoarserLevelConnectorWidthsFromFines(
      hier::IntVector& coarse_to_fine_width,
      hier::IntVector& coarse_to_coarse_width,
      const hier::IntVector& fine_to_fine_width,
      const hier::IntVector& fine_to_coarse_ratio,
      const hier::IntVector& nesting_buffer_at_fine,
      const hier::IntVector& max_stencil_width_at_coarse,
      const hier::IntVector& max_ghost_width_at_coarse) const;

   std::vector<hier::IntVector> d_tag_to_cluster_width;

};

}
}
#endif
