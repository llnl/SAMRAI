/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFillPattern
#define included_xfer_PatchLevelFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxNeighborhoodCollection.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Abstract base class for defining regions to fill on a PatchLevel.
 *
 * Class PatchLevelFillPattern is an abstract base class that provides
 * an interface used by the RefineSchedule to determine what spatial
 * regions will be filled on a given destination level.  For example,
 * a schedule may need to fill patch interiors only, ghost regions
 * only, or some combination thereof.  Concrete implementations of
 * this class will take a destination level and compute the desired
 * boxes to be filled as well as BoxNeighborhoodCollection information that
 * will be later used in the communications that fill those boxes.
 *
 * @see RefineSchedule
 * @see hier::BoxNeighborhoodCollection
 */

class PatchLevelFillPattern
{
public:
   typedef hier::BoxNeighborhoodCollection FillSet;

   /*!
    * @brief Default constructor
    */
   PatchLevelFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelFillPattern();

   /*!
    * @brief Compute the boxes to be filled and related communication data.
    *
    * This pure virtual method provides an interface to give the
    * RefineSchedule the information needed to fill particular spatial
    * regions on the destination level.  The specific regions (such as
    * patch interiors, ghost regions, or some combination thereof) will be
    * specified in the concrete implementations of this class.  Implementations
    * of this method should store the desired regions to be filled in the
    * BoxSet fill_box_level, and should compute a BoxNeighborhoodCollection
    * describing the relationship between dst_box_level and
    * fill_box_level.
    *
    * @param[out] fill_box_level       BoxLevel to be filled
    * @param[out] dst_to_fill          Connector between
    *                                  dst_box_level and fill_box_level
    * @param[in] dst_box_level         destination level
    * @param[in] fill_ghost_width      ghost width being filled by refine
    *                                  schedule
    * @param[in] data_on_patch_border  true if there is data living on patch
    *                                  borders
    */
   virtual void
   computeFillBoxesAndNeighborhoodSets(
      std::shared_ptr<hier::BoxLevel>& fill_box_level,
      std::shared_ptr<hier::Connector>& dst_to_fill,
      const hier::BoxLevel& dst_box_level,
      const hier::IntVector& fill_ghost_width,
      bool data_on_patch_border) = 0;

   /*!
    * @brief Return true if source owner can compute destination boxes on its
    * own.
    *
    * If the fill pattern is such that the source patch owner cannot compute
    * fill boxes for all of its destination neighbors using local data, then a
    * communication is needed to acquire this information, and this method
    * should return true.  If the source owner can compute the fill boxes
    * for its destination neighbors using its local data, then this method
    * should return false, which will prevent an unneeded communication step.
    *
    * If this method returns true,
    * computeDestinationFillBoxesOnSourceProc() will not be called.
    *
    * @return  true if communication needed, false otherwise.
    */
   virtual bool
   needsToCommunicateDestinationFillBoxes() const = 0;

   /*!
    * @brief Virtual method to compute the destination fill boxes.
    *
    * If needsToCommunicateDestinationFillBoxes() returns true, then this
    * method will not be called and the child classes do not need an
    * implementation.
    *
    * This method must be implemented if
    * needsToCommunicateDestinationFillBoxes() returns false, as that means
    * that the source patch owner will compute destination fill boxes using
    * its own local data with this method.
    *
    * @param[out] dst_fill_boxes_on_src_proc FillSet storing the destination
    *                                        neighbors of the source mapped
    *                                        boxes
    * @param[in] dst_box_level             destination level
    * @param[in] src_to_dst                Connector of source to destination
    * @param[in] fill_ghost_width          Ghost width being filled by refine
    *                                      schedule
    */
   virtual void
   computeDestinationFillBoxesOnSourceProc(
      FillSet& dst_fill_boxes_on_src_proc,
      const hier::BoxLevel& dst_box_level,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width) = 0;

   /*!
    * @brief Tell RefineSchedule whether it needs to communicate data
    * directly from source level to destination level.
    *
    * The return value of this method tells the RefineSchedule whether or
    * not to create an internal schedule to communicate as much data as
    * possible from the source level to the destination level on the same
    * level of resolution.  If it returns false, all of the fill boxes
    * on the destination level will be filled from a coarser level.  If
    * true, then any parts of the fill boxes that can be filled by
    * communicating data from the source level interior will be filled,
    * and whatever parts of the fill boxes are left unfilled by that
    * communication will be filled from a coarser level.
    *
    * Generally this method should return true for fill pattern implementations
    * that are used to fill the full interiors of patches and false
    * for patterns that fill only at or outside of patch boundaries.
    *
    * @return true if data should be communicated directly from source to
    * destination.
    */
   virtual bool
   doesSourceLevelCommunicateToDestination() const = 0;

   /*!
    * @brief Return the maximum number of fill boxes.
    *
    * Each destination patch has an associated set of fill boxes, and this
    * is the maximum size of those sets across all destination patches.
    *
    * The maximum number returned is for the patches on the
    * destination level provided to the most recent invocation of the
    * computeDestinationFillBoxesOnSourceProc method.
    *
    * @return maximum number of fill boxes.
    */
   virtual int
   getMaxFillBoxes() const = 0;

   /*!
    * @brief Return true if the fill pattern is intended to allow the filling
    * of ghost regions at coarse-fine boundaries.
    */
   virtual bool
   fillingCoarseFineGhosts() const = 0;

   /*!
    * @brief Return true if the fill pattern is specialized for filling
    * at enhanced connectivity block boundaries only, false for all other
    * cases.
    */
   virtual bool
   fillingEnhancedConnectivityOnly() const = 0;

private:
   PatchLevelFillPattern(
      const PatchLevelFillPattern&);                       // not implemented
   PatchLevelFillPattern&
   operator = (
      const PatchLevelFillPattern&);                     // not implemented

};

}
}
#endif
