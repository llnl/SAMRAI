/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFullFillPattern
#define included_xfer_PatchLevelFullFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief PatchLevelFullFillPattern is a PatchLevelFillPattern that
 * fills the entire region the destination level, both interior and
 * ghost.
 *
 * For documentation on this interface see @ref PatchLevelFillPattern
 *
 * The fill boxes for this PatchLevelFillPattern will consist of
 * the entire region of the destination level that can be filled, both
 * interior and ghost regions.
 *
 * If a RefineSchedule is created using an
 * RefineAlgorithm::createSchedule which takes no
 * PatchLevelFillPattern argument, this class will be used as the
 * default PatchLevelFillPattern.
 *
 * @see RefineAlgorithm
 * @see RefineSchedule
 */

class PatchLevelFullFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelFullFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelFullFillPattern();

   /*!
    * @copydoc PatchLevelFillPattern::computeFillBoxesAndNeighborhoodSets()
    *
    * The computed fill_box_level for this fill pattern will be the
    * boxes of dst_box_level grown by the fill_ghost_width.
    *
    * @pre dst_box_level.getDim() == fill_ghost_width.getDim()
    */
   void
   computeFillBoxesAndNeighborhoodSets(
      std::shared_ptr<hier::BoxLevel>& fill_box_level,
      std::shared_ptr<hier::Connector>& dst_to_fill,
      const hier::BoxLevel& dst_box_level,
      const hier::IntVector& fill_ghost_width,
      bool data_on_patch_border);

   /*!
    * @copydoc PatchLevelFillPattern::needsToCommunicateDestinationFillBoxes()
    *
    * For this fill pattern, the source owner can compute fill boxes for
    * all of its destination neighbors using local data, so this method
    * returns false, allowing a communication step to be skipped.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @copydoc PatchLevelFillPattern::computeDestinationFillBoxesOnSourceProc()
    *
    * @pre dst_box_level.getDim() == fill_ghost_width.getDim()
    */
   void
   computeDestinationFillBoxesOnSourceProc(
      FillSet& dst_fill_boxes_on_src_proc,
      const hier::BoxLevel& dst_box_level,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @copydoc PatchLevelFillPattern::doesSourceLevelCommunicateToDestination()
    *
    * RefineSchedule should attempt to fill the destination level from
    * the source level on the same resolution to the extent possible.
    */
   bool
   doesSourceLevelCommunicateToDestination() const;

   /*!
    * @copydoc PatchLevelFillPattern::getMaxFillBoxes()
    */
   int
   getMaxFillBoxes() const;

   /*!
    * @copydoc PatchLevelFillPattern::fillingCoarseFineGhosts()
    */
   bool
   fillingCoarseFineGhosts() const;

   /*!
    * @copydoc PatchLevelFillPattern::fillingEnhancedConnectivityOnly()
    */
   bool
   fillingEnhancedConnectivityOnly() const;

private:
   PatchLevelFullFillPattern(
      const PatchLevelFullFillPattern&);             // not implemented
   PatchLevelFullFillPattern&
   operator = (
      const PatchLevelFullFillPattern&);             // not implemented

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;
};

}
}

#endif
