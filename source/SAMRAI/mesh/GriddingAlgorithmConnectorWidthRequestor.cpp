/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   GriddingAlgorihtm's implementation of PatchHierarchy
 *
 ************************************************************************/
#include "SAMRAI/mesh/GriddingAlgorithmConnectorWidthRequestor.h"

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

/*
 **************************************************************************
 **************************************************************************
 */
GriddingAlgorithmConnectorWidthRequestor::GriddingAlgorithmConnectorWidthRequestor()
{
}

/*
 **************************************************************************
 * Compute Connector widths that this class requires in order to work
 * properly on a given hierarchy.
 **************************************************************************
 */
void
GriddingAlgorithmConnectorWidthRequestor::computeRequiredConnectorWidths(
   std::vector<hier::IntVector>& self_connector_widths,
   std::vector<hier::IntVector>& fine_connector_widths,
   const hier::PatchHierarchy& patch_hierarchy) const
{
   if (!d_tag_to_cluster_width.empty()) {
      TBOX_ASSERT(static_cast<int>(d_tag_to_cluster_width.size()) >=
         patch_hierarchy.getMaxNumberOfLevels() - 1);
   }

   const tbox::Dimension& dim(patch_hierarchy.getDim());
   const int max_levels(patch_hierarchy.getMaxNumberOfLevels());

   const hier::IntVector max_ghost_width(
      patch_hierarchy.getPatchDescriptor()->getMaxGhostWidth(dim));

   const hier::IntVector max_stencil_width(
      patch_hierarchy.getGridGeometry()->getMaxTransferOpStencilWidth(dim));

   fine_connector_widths.resize(max_levels - 1,
      hier::IntVector::getZero(dim));
   self_connector_widths.resize(max_levels,
      hier::IntVector::getZero(dim));


   /*
    * Compute the Connector width needed to ensure all edges are found
    * during mesh generation.  It is safe to be conservative, but
    * carrying around a larger than necessary width requires more
    * memory and slows down Connector operations.
    *
    * We compute the finest level's requirement first because coarser
    * levels' requirements depend on finer levels'.  Connector widths
    * at coarser levels are computed in
    * computeCoarserLevelConnectorWidthsFromFines() once the finer
    * level's Connector widths are computed.
    */
   self_connector_widths[max_levels - 1].setAll(max_ghost_width);
   for (int ln = max_levels - 2; ln > -1; --ln) {

      computeCoarserLevelConnectorWidthsFromFines(
         fine_connector_widths[ln],
         self_connector_widths[ln],
         self_connector_widths[ln + 1],
         patch_hierarchy.getRatioToCoarserLevel(ln + 1),
         hier::IntVector(dim, patch_hierarchy.getProperNestingBuffer(ln + 1)),
         max_stencil_width,
         max_ghost_width);

      /*
       * Must be big enough for GriddingAlgorithm::computeProperNestingData().
       */
      self_connector_widths[ln].max(
         hier::IntVector(dim, patch_hierarchy.getProperNestingBuffer(ln)));

      /*
       * Must be big enough for GriddingAlgorithm to guarantee that
       * tag--->cluster width is at least d_tag_to_cluster_width[ln]
       * when bridging cluster<==>tag<==>tag.
       */
      if (!d_tag_to_cluster_width.empty()) {
         self_connector_widths[ln].max(d_tag_to_cluster_width[ln]);
      }
   }
}

/*
 *************************************************************************
 * During mesh generation, coarse-to-fine Connectors are commonly
 * used to compute fine-to-fine Connectors.  As a result, the
 * required coarse-to-fine width is a function of the required
 * fine-to-fine width, which is why the fine-level width data is
 * required for this method.
 *
 * Given the Connector width required to support mesh operations on a
 * fine level (and other relevant parameters), compute the width
 * required of Connectors at the next coarser level.  The width at the
 * coarer level must be big enough to support mesh generation
 * operations for itself and all finer levels.  We base the
 * computations on knowledge of how the GriddingAlgorithm generates
 * the mesh.  If those change, then this method should be updated.
 *
 * The formula for computing the width for coarse-level Connectors
 * reflects how SAMRAI generates the mesh.  This formula would have to
 * change if what it is based on changes.
 *************************************************************************
 */

void
GriddingAlgorithmConnectorWidthRequestor::computeCoarserLevelConnectorWidthsFromFines(
   hier::IntVector& coarse_to_fine_width,
   hier::IntVector& coarse_to_coarse_width,
   const hier::IntVector& fine_to_fine_width,
   const hier::IntVector& fine_to_coarse_ratio,
   const hier::IntVector& nesting_buffer_at_fine,
   const hier::IntVector& max_stencil_width_at_coarse,
   const hier::IntVector& max_ghost_width_at_coarse) const
{
   NULL_USE(max_stencil_width_at_coarse);

#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   const tbox::Dimension& dim(fine_to_fine_width.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY7(dim,
      coarse_to_fine_width,
      coarse_to_coarse_width,
      fine_to_fine_width,
      fine_to_coarse_ratio,
      nesting_buffer_at_fine,
      max_stencil_width_at_coarse,
      max_ghost_width_at_coarse);
#endif

   const size_t num_blocks = fine_to_coarse_ratio.getNumBlocks();

   /*
    * Coarse-to-fine width must be big enough to cover the width at the
    * finer level.  For example, the coarse level is used to bridge for
    * the fine-to-fine Connectors.  All requirements for width at the
    * finer level and above should be reflected in the fine-to-fine
    * width.
    */
   coarse_to_fine_width =
      hier::IntVector::ceilingDivide(fine_to_fine_width, fine_to_coarse_ratio);
   /*
    * Coarse-to-fine width must be big enough for the [ln] -> [ln+1]
    * Connector to see all the [ln+1] Boxes that are used to add
    * tags to [ln-1] for ensuring [ln] properly nests
    * [ln+1] when [ln] is being regenerated.
    *
    * The rationale for this adjustment is illustrated by the following
    * worst case example in 1D.  We are tagging on L0 to regenerate L1.
    * We must tag wherever [L2 grown by nesting buffer] intersects L0.
    *
    *          |--|         L2, 1-cell wide in this example.
    *
    *    |..|..|--|..|..|   L2, plus 2-cell nesting buffer
    *
    *          |--|         Current L1, to be replaced.
    *
    *    |--------------|   Minimal L1 to properly nest L2
    *
    *    |--|--|--|--|--|   Current L0, with 5 cells.  All cells must
    *                       be fully tagged to ensure new L1 meets it
    *                       minimal coverage.
    *
    * To determine where L0 must be tagged, we must bridge for
    * Connector L0--->L2.  And to make sure the left-most patch of L0
    * sees the patch on L2 via this Connector, Connector L0--->L1 must
    * be have at least a 2-cell GCW (in L2 index space).
    */
   coarse_to_fine_width.max(
      hier::IntVector::ceilingDivide(
         hier::IntVector(nesting_buffer_at_fine, num_blocks),
         fine_to_coarse_ratio));

   /*
    * Coarse-to-coarse Connectors must cover all data that the finer
    * level depends on.
    */
   coarse_to_coarse_width = coarse_to_fine_width;

   /*
    * Must cover coarse level's own ghost cells.
    */
   coarse_to_coarse_width.max(max_ghost_width_at_coarse);
}

/*
 **************************************************************************
 **************************************************************************
 */
void
GriddingAlgorithmConnectorWidthRequestor::setTagToClusterWidth(
   std::vector<hier::IntVector>& tag_to_cluster_width)
{
   d_tag_to_cluster_width = tag_to_cluster_width;
}

}
}
