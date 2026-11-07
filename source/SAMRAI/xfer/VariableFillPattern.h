/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_VariableFillPattern
#define included_xfer_VariableFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/PatchDataFactory.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class VariableFillPattern is an abstract base class that provides an
 * interface to create objects that can calculate overlaps which correspond
 * to a specific stencil.
 *
 * If an object of a concrete VariableFillPattern type is provided when
 * registering a refine operation with a RefineAlgorithm, then BoxOverlap
 * calculations associated with that refine operation will use the
 * VariableFillPattern implementation to restrict the overlaps to only a
 * desired subset of the the intersection between the source and destination
 * patches.  Thus only data within that subset will be filled when
 * RefineSchedule::fillData() is called.  For example, an implementation of
 * this class may be used to restrict the filling of data to locations on or
 * near a patch boundary.
 *
 * @see hier::BoxOverlap
 * @see RefineAlgorithm
 * @see RefineSchedule
 */

class VariableFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   VariableFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~VariableFillPattern();

   /*!
    * @brief This pure virtual method provides an interface to calculate
    * overlaps between the destination and source geometries for a copy
    * operation.
    *
    * Implementations of this method will restrict the calculated BoxOverlap
    * to a certain subset of the intersection between the destination and
    * source geometries.
    *
    * @param[in] dst_geometry    geometry object for destination box
    * @param[in] src_geometry    geometry object for source box
    * @param[in] dst_patch_box   box for the destination patch
    * @param[in] src_mask        the source mask, the box resulting from
    *                            transforming the source box
    * @param[in] fill_box        the box to be filled
    * @param[in] overwrite_interior  controls whether or not to include the
    *                                destination box interior in the overlap
    * @param[in] transformation  the transformation from source to
    *                            destination index space.
    *
    * @return                    std::shared_ptr to the calculated overlap
    *                            object
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   calculateOverlap(
      const hier::BoxGeometry& dst_geometry,
      const hier::BoxGeometry& src_geometry,
      const hier::Box& dst_patch_box,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation) const = 0;

   /*!
    * @brief This pure virtual method provides an interface for computing
    * overlaps which define the space to be filled by a refinement operation.
    *
    * Implementations of this method compute a BoxOverlap that covers a
    * desired subset of the space defined by fill_boxes.
    *
    * @param[in] fill_boxes  list representing all of the space on a patch
    *                        or its ghost region that may be filled by a
    *                        refine operator (cell-centered represtentation)
    * @param[in] node_fill_boxes node-centered representation of fill_boxes
    * @param[in] patch_box   box representing the patch where a refine operator
    *                        will fill data.  (cell-centered representation)
    * @param[in] data_box    box representing the full extent of the region
    *                        covered by a patch data object, including all
    *                        ghosts (cell-centered representation)
    * @param[in] patch_data_factory patch data factory for the data that is to
    *                               be filled
    *
    * @return                std::shared_ptr to the calculated overlap object
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   computeFillBoxesOverlap(
      const hier::BoxContainer& fill_boxes,
      const hier::BoxContainer& node_fill_boxes,
      const hier::Box& patch_box,
      const hier::Box& data_box,
      const hier::PatchDataFactory& patch_data_factory) const = 0;

   /*!
    * @brief Return the maximum ghost width of the stencil defined by the
    * VariableFillPattern implementation.
    */
   virtual const hier::IntVector&
   getStencilWidth() = 0;

   /*!
    * @brief Return a string name identifying the concrete subclass.
    */
   virtual const std::string&
   getPatternName() const = 0;

private:
   VariableFillPattern(
      const VariableFillPattern&);                     // not implemented
   VariableFillPattern&
   operator = (
      const VariableFillPattern&);                     // not implemented

};

}
}
#endif
