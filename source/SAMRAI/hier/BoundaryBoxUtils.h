/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Generic utilities for boundary box calculus.
 *
 ************************************************************************/

#ifndef included_hier_BoundaryBoxUtils
#define included_hier_BoundaryBoxUtils

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BoundaryBox.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Perform shifts, extensions, etc on a BoundaryBox using the box's
 * location index and type.
 *
 * @see BoundaryBox
 */

class BoundaryBoxUtils
{

public:
   /*!
    * @brief Construct with a boundary box.
    *
    * @param[in]  bbox  boundary box
    *
    * @see setBoundaryBox()
    */
   explicit BoundaryBoxUtils(
      const BoundaryBox& bbox);

   /*!
    * @brief Destructor.
    */
   ~BoundaryBoxUtils() = default;

   /*!
    * @brief Reset boundary box.
    *
    * All utility operations refer to this box.
    *
    * @param[in]  bbox  boundary box
    */
   void
   setBoundaryBox(
      const BoundaryBox& bbox)
   {
      d_bbox = bbox;
      computeOutwardShift();
   }

   /*!
    * @brief Get boundary box.
    *
    * @return The boundary box
    */
   const BoundaryBox&
   getBoundaryBox() const
   {
      return d_bbox;
   }

   /*!
    * @brief Get the outward direction in logical space.
    *
    * Each component of the returned outward direction will have
    * one of these values:
    *
    * <ul>
    *   <li> -1 if the outward direction is toward the lower indices
    *   <li> 0 for the direction is orthogonal to the outward direction.
    *   <li> 1 if the outward direction is toward the higher indices
    * </ul>
    *
    * The number of non-zero components should equal the boundary
    * type (codimension).
    *
    * @return IntVector containing the outward direction values
    */
   const IntVector&
   getOutwardShift() const
   {
      return d_outward;
   }

   /*!
    * @brief Stretch box outward by the given ghost cell width.
    *
    * The number of directions affected is the same as the
    * codimension of the boundary.
    *
    * Note that the BoundaryBox is defined to be one cell wide.  The
    * output of this method is the box held by the BoundaryBox, stretched to
    * cover the given ghost cell width.  This means that if ghost_cell_width
    * is one, the output is identical to the BoundaryBox.  If the ghost
    * width is zero in any direction, the output will shrink to nothing
    * in that direction.
    *
    * @param[out] box The stretched box.
    * @param[in] ghost_cell_width Ghost width to stretch the box.
    *
    * @pre getBoundaryBox().getDim() == box.getDim()
    * @pre ghost_cell_width >= IntVector::getZero(getBoundaryBox().getDim())
    */
   void
   stretchBoxToGhostWidth(
      Box& box,
      const IntVector& ghost_cell_width) const;

   /*!
    * @brief Extend box outward by the given amount.
    *
    * The output box is the box held by the BoundaryBox extended outward
    * according to the values in the extension IntVector.  If any
    * coordinate direction has a zero value in the IntVector returned by
    * getOutwardShift(), the box will not be extended in that direction.
    * The number of directions that will be extended is equal to the
    * codimension of the BoundaryBox.
    *
    * @param[out] box      The extended box
    * @param[in] extension IntVector telling how many cells to extend the box
    *                      in each coordinate direction.
    *
    * @pre getBoundaryBox().getDim() == box.getDim()
    */
   void
   extendBoxOutward(
      Box& box,
      const IntVector& extension) const;

   /*!
    * @brief Return the direction normal to the BoundaryBox.
    *
    * The normal direction is defined only for surface
    * boundaries (codimension 1).  A -1 is returned for
    * all other boundary types.
    *
    * @return The normal direction.
    */
   int
   normalDir() const
   {
      return getBoundaryBox().getLocationIndex() / 2;
   }

   /*!
    * @brief Trim a boundary box so that it does not stick out
    * past a limiting box in direction transverse to the boundary
    * normal.
    *
    * This method affects the only box directions parallel to
    * the boundary.  For methods affecting other box directions,
    * see stretchBoxToGhostWidth().
    *
    * The boundary type of the BoundaryBox that was given to the
    * BoundaryBoxUtils constructor must be less than dim.
    *
    * @param[in] limit_box Box to not stick out past
    *
    * @return New trimmed boundary box.
    *
    * @pre getBoundaryBox().getDim() == limit_box.getDim()
    * @pre getBoundaryBox().getBoundaryType() < getBoundaryBox().getDim().getValue())
    */
   BoundaryBox
   trimBoundaryBox(
      const Box& limit_box) const;

   /*!
    * @brief Return box describing the index space of the outer surface of
    * a boundary box.
    *
    * Define a box describing the indices of the surface of the
    * the input boundary box.  A surface is a face in 3D and an edge
    * in 2D.  These surfaces lie on the boundary itself.
    *
    * The input boundary_box must be of type 1
    * (see BoundaryBox::getBoundaryType()).
    *
    * This is a utility function for working with the surface
    * indices corresponding to a boundary box.
    *
    * @return a box to define the side indices corresponding to the
    * BoundaryBox
    *
    * @pre getBoundaryBox().getBoundaryType() == 1
    */
   Box
   getSurfaceBoxFromBoundaryBox() const;

private:
   // Unimplemented default constructor.
   BoundaryBoxUtils();

   /*!
    * @brief Compute the shift in the outward direction
    * (redundant data, function of boundary box) and store
    * in d_outward.
    *
    * @see getOutwardShift();
    */
   void
   computeOutwardShift();

   /*!
    * @brief Boundary box implicitly referred to by all methods.
    *
    * @see setBoundaryBox(), getBoundaryBox()
    */
   BoundaryBox d_bbox;

   /*!
    * @brief Vector pointing outward from patch.
    */
   IntVector d_outward;

};

}
}

#endif
