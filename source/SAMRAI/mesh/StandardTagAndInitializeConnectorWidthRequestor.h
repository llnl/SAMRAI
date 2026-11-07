/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   StandardTagAndInitialize's implementation of PatchHierarchy
 *
 ************************************************************************/

#ifndef included_mesh_StandardTagAndInitializeConnectorWidthRequestor
#define included_mesh_StandardTagAndInitializeConnectorWidthRequestor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide StandardTagAndInitialize needs
 * Connectors between hierarchy levels to be.
 *
 * To do Richardson extrapolation, StandardTagAndInitialize will
 * coarsen a level and populate it with data.  A coarsened level has a
 * bigger ghost region because the coarse cells are bigger.  This
 * class is for telling the PatchHierarchy that
 * StandardTagAndInitialize will request Connectors based on the width
 * of the coarsened level.
 */
class StandardTagAndInitializeConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   StandardTagAndInitializeConnectorWidthRequestor();

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

private:
   /*!
    * @brief Return the coarsen ratio to be used with Richardson
    * extrapolation on the PatchHierarchy.
    *
    * @param[in] ratios_to_coarser Refinement ratios in a hierarchy.
    * @c ratios_to_coarser[ln] is the ratio between level ln and level
    * ln-1.
    *
    * @pre (ratios_to_coarser[1](0) % 2 == 0) ||
    *      (ratios_to_coarser[1](0) % 3 == 0)
    */
   int
   computeCoarsenRatio(
      const std::vector<hier::IntVector>& ratios_to_coarser) const;

};

}
}
#endif
