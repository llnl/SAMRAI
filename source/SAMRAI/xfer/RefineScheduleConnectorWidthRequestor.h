/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   RefineSchedule's implementation of PatchHierarchy
 *
 ************************************************************************/

#ifndef included_xfer_RefineScheduleConnectorWidthRequestor
#define included_xfer_RefineScheduleConnectorWidthRequestor

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide RefineSchedule needs Connectors
 * between hierarchy levels to be.
 */
class RefineScheduleConnectorWidthRequestor:
   public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   RefineScheduleConnectorWidthRequestor();

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
    * @brief Set the factor by which the ghost data widths in the
    * PatchHierarchy are multiplied to increase their effective values.
    * for the purpose of handling more ghost data than registered.
    *
    * We support values of ghost data width that are larger than those
    * registered by some factor.  An example of the need is when we
    * fill a PatchLevel that is a coarsened version of a current level
    * in the hierarchy.  Coarsening the level effectively increases
    * the size of the ghost data region because the coarse cells are
    * bigger.
    *
    * Note that the Connector widths are not linear functions of this
    * factor.  Increasing the factor is NOT the same as multiplying
    * the Connector widths by the same factor.
    *
    * @param[in] gcw_factor By default, @c gcw_factor=1.
    *
    * @pre gcw_factor >= 0
    */
   void
   setGhostCellWidthFactor(
      int gcw_factor);

   /*
    * @brief Compute fine Connector width needed at each coarser level
    * for recursive refinement starting with destination level ln.
    *
    * @param fine_connector_width [in,out] On input, contains width
    * requirement from other conditions, if any.  This method computes
    * the fine connector width for each level coarser than
    * initial_dst_ln, for recursive refinement with initial
    * destination level initial_dst_ln.  It sets fine_connector_width
    * to the max of its input value and the value computed.  It will
    * change the values of the first initial_dst_ln entries (or create
    * them if missing).
    *
    * @param data_gcw_on_initial_dst_ln Ghost data width to be filled
    * on the initial destination level.  Important: This input value
    * should include the 1 cell added for data on patch borders, if it
    * is needed.
    *
    * @param [out] max_stencil_width Max stencil width used by the
    * refinement operators.
    *
    * @param [in] patch_hierarchy
    *
    * @param [in] initial_dst_ln Level number corresponding to the
    * initial destination.
    */
   void
   computeRequiredFineConnectorWidthsForRecursiveRefinement(
      std::vector<hier::IntVector>& fine_connector_widths,
      const hier::IntVector& data_gcw_on_initial_dst_ln,
      const hier::IntVector& max_stencil_width,
      const hier::PatchHierarchy& patch_hierarchy,
      int initial_dst_ln) const;

private:
   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      hier::PatchHierarchy::registerAutoConnectorWidthRequestorStrategy(
         s_auto_registered_connector_width_requestor);
   }

   /*!
    * Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
   }

   /*!
    * @brief Static object that is auto-registered in PatchHierarchy
    * by default.
    */
   static RefineScheduleConnectorWidthRequestor
      s_auto_registered_connector_width_requestor;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

   /*
    * @brief The factor by which the ghost data widths in the
    * PatchHierarchy are multiplied to increase their effective
    * values.  for the purpose of handling more ghost data than
    * registered.
    */
   int d_gcw_factor;

};

}
}

#endif
