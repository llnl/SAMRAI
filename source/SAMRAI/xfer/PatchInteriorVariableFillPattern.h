/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class that fills PatchInteriors only
 *
 ************************************************************************/

#ifndef included_xfer_PatchInteriorVariableFillPattern
#define included_xfer_PatchInteriorVariableFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/VariableFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class PatchInteriorVariableFillPattern is an implementation of the
 * abstract base class VariableFillPattern.
 *
 * It is used to calculate overlaps hat consist of the intersections
 * between patches, including only patch interiors.  (The boundaries of a patch
 * are considered to be part of the interior.)
 */

class PatchInteriorVariableFillPattern:
   public VariableFillPattern
{
public:
   /*!
    * @brief Constructor
    *
    * @param[in] dim     Dimension
    */
   explicit PatchInteriorVariableFillPattern(
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   virtual ~PatchInteriorVariableFillPattern();

   /*!
    * @brief Calculate overlap between the destination and source geometries
    * using the geometries' own overlap calculation methods.
    *
    * The intersection between the given dst_geometry and src_geometry
    * will be calculated according to the properties of those geometries, but
    * will not include any ghost space outside of the patch represented by
    * dst_patch_box.
    *
    * If the argument overwrite_interior is set to false, then the returned
    * BoxOverlap will be empty, since this class only computes overlaps
    * representing patch interiors.
    *
    * @param[in] dst_geometry    geometry object for destination box
    * @param[in] src_geometry    geometry object for source box
    * @param[in] dst_patch_box   box for the destination patch
    * @param[in] src_mask        the source mask, the box resulting from
    *                            transforming the source box
    * @param[in] fill_box        the box to be filled
    * @param[in] overwrite_interior  controls whether or not to include the
    *                                destination box interior in the overlap.
    * @param[in] transformation  the transformation from source to
    *                            destination index space.
    *
    * @return                    std::shared_ptr to the calculated overlap
    *                            object
    *
    * @pre dst_patch_box.getDim() == src_mask.getDim()
    */
   std::shared_ptr<hier::BoxOverlap>
   calculateOverlap(
      const hier::BoxGeometry& dst_geometry,
      const hier::BoxGeometry& src_geometry,
      const hier::Box& dst_patch_box,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation) const;

   /*!
    * Computes a BoxOverlap object which defines the space to be filled by
    * a refinement operation.  For this implementation, that space is the
    * intersection between fill_boxes (computed by the RefineSchedule),
    * data_box, which specifies the extent of the destination data, and
    * patch_box, which specifies the patch interior of the destination.  The
    * patch data factory is used to compute the overlap with the appropriate
    * data centering, consistent with the centering of the data to be filled.
    *
    * @param[in] fill_boxes  list representing the all of the space on a patch
    *                        or its ghost region that may be filled by a
    *                        refine operator (cell-centered represtentation)
    * @param[in] unfilled_node_boxes node-centered representation of fill_boxes
    * @param[in] patch_box   box representing the patch where a refine operator
    *                        will fill data.  (cell-centered representation)
    * @param[in] data_box    box representing the full extent of the region
    *                        covered by a patch data object, including all
    *                        ghosts (cell-centered representation)
    * @param[in] pdf         patch data factory for the data that is to be
    *                        filled
    */
   std::shared_ptr<hier::BoxOverlap>
   computeFillBoxesOverlap(
      const hier::BoxContainer& fill_boxes,
      const hier::BoxContainer& unfilled_node_boxes,
      const hier::Box& patch_box,
      const hier::Box& data_box,
      const hier::PatchDataFactory& pdf) const;

   /*!
    * @brief Implementation of interface to get stencil width.
    *
    * @return    Zero IntVector since this fills patch interiors only.
    */
   const hier::IntVector&
   getStencilWidth();

   /*!
    * @brief Returns a string name identifier "PATCH_INTERIOR_FILL_PATTERN".
    */
   const std::string&
   getPatternName() const;

private:
   PatchInteriorVariableFillPattern(
      const PatchInteriorVariableFillPattern&);    // not implemented
   PatchInteriorVariableFillPattern&
   operator = (
      const PatchInteriorVariableFillPattern&);    // not implemented

   /*!
    * @brief The dimension of this object.
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Static string holding name identifier for this class.
    */
   static const std::string s_name_id;
};

}
}

#endif
