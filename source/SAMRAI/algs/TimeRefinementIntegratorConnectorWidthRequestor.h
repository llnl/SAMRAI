/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   TimeRefinementIntegrator's implementation of PatchHierarchy::ConnectorWidthRequestorStrategy
 *
 ************************************************************************/

#ifndef included_mesh_TimeRefinementIntegratorConnectorWidthRequestor
#define included_mesh_TimeRefinementIntegratorConnectorWidthRequestor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

namespace SAMRAI {
namespace algs {

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide TimeRefinementIntegrator needs Connectors
 * between hierarchy levels to be.
 */
class TimeRefinementIntegratorConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   TimeRefinementIntegratorConnectorWidthRequestor();

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
    * given tag buffer.
    *
    * This is only used by the TimeRefinementIntegrator that owns this object.
    */
   void
   setTagBuffer(
      const std::vector<int>& tag_buffer);

private:
   std::vector<int> d_tag_buffer;

};

}
}
#endif
