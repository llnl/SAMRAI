/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Level fill pattern for enhanced connectivity
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelEnhancedFillPattern
#define included_xfer_PatchLevelEnhancedFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/hier/Connector.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class PatchLevelEnhancedFillPattern is an implementation of the
 * abstract base class PatchLevelFillPattern.
 *
 * This class is used by the RefineSchedule to restrict filling to
 * only patches in ghost regions across an ehanced connectivity block boundary.
 * It is intended for users who wish to handle the filling of data at these
 * singularities separately from the filling of all other data.
 *
 * @see RefineSchedule
 */

class PatchLevelEnhancedFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelEnhancedFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelEnhancedFillPattern();

   /*!
    * @brief Compute the boxes representing the region that will
    *        be filled by a RefineSchedule.
    *
    * The fill boxes will exist only for destination patches that touch
    * a block boundary at enhanced connectivity.  The fill boxes will be
    * the portion of the ghost regions of those patches lying across the
    * enhanced connectivity boundary.
    *
    * @param[out] fill_box_level       Output BoxLevel to be filled
    * @param[out] dst_to_fill          Output Connector between
    *                                  dst_box_level and fill_box_level
    * @param[in] dst_box_level         destination level
    * @param[in] fill_ghost_width      Ghost width being filled by refine
    *                                  schedule
    * @param[in] data_on_patch_border  true if there is data living on patch
    *                                  borders
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
    * @brief Return true because source patch owner cannot compute fill
    * boxes across block boundaries.
    *
    * @return  true.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @brief Virtual method to compute the destination fill boxes.
    *
    * Since needsToCommunicateDestinationFillBoxes() returns true, this
    * method should never be called.  It is implemented here to satisfy
    * the pure virtual interface from the base class.  An error will result
    * if this is ever called.
    *
    * @pre needsToCommunicateDestinationFillBoxes()
    */
   void
   computeDestinationFillBoxesOnSourceProc(
      FillSet& dst_fill_boxes_on_src_proc,
      const hier::BoxLevel& dst_box_level,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @brief Tell RefineSchedule not to communicate data directly from source
    * to destination level.
    *
    * With this fill pattern, the RefineSchedule does not attempt to fill
    * as much of the destination level as possible from the source level at
    * the same level of resolution.  By returning 'false', this method
    * tells the RefineSchedule to skip that step.
    *
    * @return false
    */
   bool
   doesSourceLevelCommunicateToDestination() const;

   /*!
    * @brief Return the maximum number of fill boxes.
    *
    * This will not return a valid value until
    * computeFillBoxesAndNeighborhoodSets is fully implemented.  An
    * error will occur if this method is called
    */
   int
   getMaxFillBoxes() const;

   /*!
    * @brief Returns true because this fill pattern fills coarse fine ghost
    * data.
    */
   bool
   fillingCoarseFineGhosts() const;

   /*!
    * @brief Returns true because this fill pattern is specialized for
    * enhanced connectivity only.
    */
   bool
   fillingEnhancedConnectivityOnly() const;

private:
   PatchLevelEnhancedFillPattern(
      const PatchLevelEnhancedFillPattern&);           // not implemented

   PatchLevelEnhancedFillPattern&
   operator = (
      const PatchLevelEnhancedFillPattern&);           // not implemented

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;
};

}
}

#endif
