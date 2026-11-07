/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   RefineSchedule's implementation of PatchHierarchy
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineScheduleConnectorWidthRequestor.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

tbox::StartupShutdownManager::Handler
RefineScheduleConnectorWidthRequestor::s_initialize_finalize_handler(
   RefineScheduleConnectorWidthRequestor::initializeCallback,
   0,
   0,
   RefineScheduleConnectorWidthRequestor::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

RefineScheduleConnectorWidthRequestor
RefineScheduleConnectorWidthRequestor::s_auto_registered_connector_width_requestor;

/*
 **************************************************************************
 **************************************************************************
 */
RefineScheduleConnectorWidthRequestor::RefineScheduleConnectorWidthRequestor(
   ):d_gcw_factor(1)
{
}

/*
 **************************************************************************
 **************************************************************************
 */
void
RefineScheduleConnectorWidthRequestor::setGhostCellWidthFactor(
   int gcw_factor)
{
   TBOX_ASSERT(gcw_factor >= 0);
   d_gcw_factor = gcw_factor;
}

/*
 **************************************************************************
 * Compute Connector widths that this class requires in order to work
 * properly on a given hierarchy.
 **************************************************************************
 */
void
RefineScheduleConnectorWidthRequestor::computeRequiredConnectorWidths(
   std::vector<hier::IntVector>& self_connector_widths,
   std::vector<hier::IntVector>& fine_connector_widths,
   const hier::PatchHierarchy& patch_hierarchy) const
{
   int max_levels = patch_hierarchy.getMaxNumberOfLevels();

   const tbox::Dimension& dim(patch_hierarchy.getDim());

   /*
    * Add one to max data ghost width to create overlaps of data
    * living on patch boundaries.
    */
   const hier::IntVector max_data_gcw(
      patch_hierarchy.getPatchDescriptor()->getMaxGhostWidth(dim) + 1);

   hier::IntVector max_stencil_width =
      patch_hierarchy.getGridGeometry()->getMaxTransferOpStencilWidth(dim);
   max_stencil_width.max(
      RefinePatchStrategy::getMaxRefineOpStencilWidth(dim));

   hier::IntVector zero_vector(hier::IntVector::getZero(dim),
                               patch_hierarchy.getNumberBlocks());

   /*
    * Compute the Connector width needed to ensure all edges are found
    * during mesh recursive refine schedule generation.  It is safe to
    * be conservative, but carrying around a larger than necessary
    * width requires more memory and slows down Connector operations.
    *
    * All Connectors to self need to be at least wide enough to
    * support the copy of data from the same level into ghost cells.
    * Thus, the width should be at least that of the max ghost data
    * width.  On the finest level, there is no other requirement.  For
    * other levels, we need enough width for:
    *
    * - refining the next finer level
    *
    * - refining recursively starting at each of the levels finer than
    *   it.
    */

   hier::IntVector self_width(max_data_gcw * d_gcw_factor,
                              patch_hierarchy.getNumberBlocks()); 
   self_connector_widths.clear();
   self_connector_widths.resize(max_levels, self_width);

   fine_connector_widths.clear();
   if (max_levels > 1) {
      fine_connector_widths.resize(max_levels - 1, zero_vector); // to be computed below.
   }

   /*
    * Note that the following loops go from fine to coarse.  This is
    * because Connector widths for coarse levels depend on those for
    * fine levels.
    */
   for (int ln = max_levels - 1; ln > -1; --ln) {
      computeRequiredFineConnectorWidthsForRecursiveRefinement(
         fine_connector_widths,
         max_data_gcw,
         max_stencil_width,
         patch_hierarchy,
         ln);
   }

}

/*
 **************************************************************************
 * Compute fine Connector width needed at each coarser level (lnc) for
 * recursive refinement starting with destination level ln.
 **************************************************************************
 */
void
RefineScheduleConnectorWidthRequestor::computeRequiredFineConnectorWidthsForRecursiveRefinement(
   std::vector<hier::IntVector>& fine_connector_widths,
   const hier::IntVector& data_gcw_on_initial_dst_ln,
   const hier::IntVector& max_stencil_width,
   const hier::PatchHierarchy& patch_hierarchy,
   int initial_dst_ln) const
{
   if (static_cast<int>(fine_connector_widths.size()) < initial_dst_ln) {
      fine_connector_widths.insert(
         fine_connector_widths.end(),
         initial_dst_ln - fine_connector_widths.size(),
         hier::IntVector(patch_hierarchy.getDim(), 0,
            patch_hierarchy.getGridGeometry()->getNumberBlocks()) );
   }

   const size_t nblocks = patch_hierarchy.getGridGeometry()->getNumberBlocks();

   hier::IntVector width_for_refining_recursively(
      data_gcw_on_initial_dst_ln * d_gcw_factor, nblocks);

   for (int lnc = initial_dst_ln - 1; lnc > -1; --lnc) {

      const hier::IntVector& ratio_to_coarser =
         patch_hierarchy.getRatioToCoarserLevel(lnc + 1);
      width_for_refining_recursively.ceilingDivide(ratio_to_coarser);

      /*
       * Data in the supplemental level in RefineSchedule may have ghost
       * cells as big as the stencil width.  Coarse_to_fine_width must be
       * big enough to allow coarse to bridge to fine's supplemental, and
       * the supplemental includes the stencil width at coarse.
       */
      width_for_refining_recursively += max_stencil_width;

      fine_connector_widths[lnc].max(width_for_refining_recursively);
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
