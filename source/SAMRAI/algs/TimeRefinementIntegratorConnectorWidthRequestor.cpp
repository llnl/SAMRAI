/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   TimeRefinementIntegrator's implementation of PatchHierarchy::ConnectorWidthRequestorStrategy
 *
 ************************************************************************/
#include "SAMRAI/algs/TimeRefinementIntegratorConnectorWidthRequestor.h"

#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace algs {

/*
 **************************************************************************
 **************************************************************************
 */
TimeRefinementIntegratorConnectorWidthRequestor::TimeRefinementIntegratorConnectorWidthRequestor()
{
}

/*
 **************************************************************************
 * Compute Connector widths that this class requires in order to work
 * properly on a given hierarchy.
 *
 * The only TimeRefinementIntegrator requirement is enough
 * self_connector_widths to fill the tag buffer that the integrator
 * passes to its GriddingAlgorithm.  For some reason, this ghost width
 * was not registered at the time the required Connector widths are
 * computed.  This appeared to be by design (see how it uses
 * GriddingAlgorithm::resetTagBufferingData), so I didn't change it,
 * but it probably should be redesigned.
 *
 * The fine connector widths are set to the self connector width divided
 * by the refinement ratio for each level (rounded up if there's a remainder).
 * This ensures that if the tag buffer causes a large self connector width,
 * the fine connector width is large enough to support bridge operations.
 **************************************************************************
 */
void
TimeRefinementIntegratorConnectorWidthRequestor::computeRequiredConnectorWidths(
   std::vector<hier::IntVector>& self_connector_widths,
   std::vector<hier::IntVector>& fine_connector_widths,
   const hier::PatchHierarchy& patch_hierarchy) const
{
   const tbox::Dimension& dim(patch_hierarchy.getDim());
   const int max_levels(patch_hierarchy.getMaxNumberOfLevels());

   fine_connector_widths.resize(max_levels - 1, hier::IntVector(dim, 0, patch_hierarchy.getNumberBlocks()));
   self_connector_widths.clear();
   self_connector_widths.reserve(max_levels);
   for (size_t ln = 0; ln < static_cast<size_t>(max_levels); ++ln) {
      hier::IntVector buffer(
         dim,
         d_tag_buffer.size() > ln ? d_tag_buffer[ln] : d_tag_buffer.back(),
         patch_hierarchy.getNumberBlocks());
      self_connector_widths.push_back(buffer);
      if (ln < static_cast<size_t>(max_levels-1)) {
         const hier::IntVector& ratio_to_coarse =
            patch_hierarchy.getRatioToCoarserLevel(ln+1);
        for (hier::BlockId::block_t b = 0;
              b < patch_hierarchy.getNumberBlocks(); ++b) {
            for (int d = 0; d < dim.getValue(); ++d) {
               int& fine_w = fine_connector_widths[ln](b,d);
               const int& self_w = self_connector_widths[ln](b,d);
               const int& ratio = ratio_to_coarse(b,d);
               if (fine_w * ratio < self_w) {
                  fine_w = self_w / ratio;
                  if (self_w % ratio) {
                     ++fine_w;
                  }
               }
            }
         }
      }
   }
}

/*
 **************************************************************************
 **************************************************************************
 */
void
TimeRefinementIntegratorConnectorWidthRequestor::setTagBuffer(
   const std::vector<int>& tag_buffer)
{
   d_tag_buffer = tag_buffer;
}

}
}
