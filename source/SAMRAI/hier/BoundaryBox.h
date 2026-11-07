/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Box representing a portion of the AMR index space
 *
 ************************************************************************/

#ifndef included_hier_BoundaryBox
#define included_hier_BoundaryBox

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class BoundaryBox is is used to describe boundaries of a patch.
 *
 * Objects of this type are held by a PatchGeometry object.
 * The BoundaryBox consists of a Box, a boundary type (codimension), and a
 * location index.  The Box is one cell wide in at least one direction and is
 * located just outside of a patch boundary.  For example, a bondary box
 * along a patch face is one cell wide in the coordinate direction normal
 * to the patch face.  The boundary type identifies the type of patch bounadry:
 * face, edge, or corner (node).  The location index specifies the location
 * of the boundary box in relation to the patch boundary.
 * See the getBoundaryType() and getLocationIndex() methods for more
 * information.
 *
 * @see Box
 * @see PatchGeometry
 * @see BoundaryLookupTable
 */

class BoundaryBox
{
public:
   /*!
    * @brief Constructor that builds an undefined boundary box, with invalid
    * values other than dimension.
    *
    * @param[in] dim Dimension
    */
   explicit BoundaryBox(
      const tbox::Dimension& dim);

   /*!
    * @brief Construct a boundary box described by arguments.
    *
    * Assertions:  bdry_type must be between 1 and dim, inclusive.
    * location_index must be non-negative and less than the location index max,
    * which is the total number of boundaries of a particular codimension that
    * a single patch has.  This max value can be obtained from the class
    * BoundaryLookupTable using the getMaxLocationIndices method.
    *
    * @param[in] box
    * @param[in] bdry_type
    * @param[in] location_index
    *
    * @pre (bdry_type >= 1) && (bdry_type <= d_dim.getValue())
    * @pre location_index >= 0
    * @pre location_index < BoundaryLookupTable::getLookupTable(box.getDim())->getMaxLocationIndices()[bdry_type - 1]
    */
   BoundaryBox(
      const Box& box,
      const int bdry_type,
      const int location_index);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] boundary_box
    */
   BoundaryBox(
      const BoundaryBox& boundary_box) = default;

   BoundaryBox(
      BoundaryBox&& boundary_box) = default;

   /*!
    * @brief The destructor for BoundaryBox.
    */
   ~BoundaryBox() = default;

   /*!
    * @brief Return the Box member of the boundary box
    *
    * @return The Box
    */
   const Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the boundary type (codimension) of the boundary box.
    *
    * \verbatim
    * Convention:
    * ===========
    *
    * 1d
    * --
    * 1 = node
    *
    * 2d
    * --
    * 1 = edge
    * 2 = node
    *
    * 3d
    * --
    * 1 = face
    * 2 = edge
    * 3 = node
    * \endverbatim
    *
    * @return boundary type value
    */
   int
   getBoundaryType() const
   {
      return d_bdry_type;
   }

   /*!
    * @brief Return the location index for the boundary box.
    *
    * The location index is an integer which indicates the location of the
    * the boundary box in relation to the location of the associated patch.
    * The conventions for the location index depend on the dimension of
    * the problem and the boundary type (codimension) of the BoundaryBox.
    *
    * \verbatim
    * Conventions:
    * ============
    *
    * 1d
    * --
    * node (codimension 1):
    * x_lo : 0
    * x_hi : 1
    *
    * 2d
    * --
    * edge (codimension 1):
    * x_lo: 0
    * x_hi: 1
    * y_lo: 2
    * y_hi: 3
    *
    * node (codimension 2):
    * x_lo, y_lo: 0
    * x_hi, y_lo: 1
    * x_lo, y_hi: 2
    * x_hi, y_hi: 3
    *
    * 3d
    * --
    *
    * face (codimension 1):
    * x_lo: 0
    * x_hi: 1
    * y_lo: 2
    * y_hi: 3
    * z_lo: 4
    * z_hi: 5
    *
    * edge (codimension 2):
    * x_lo, y_lo: 0
    * x_hi, y_lo: 1
    * x_lo, y_hi: 2
    * x_hi, y_hi: 3
    * x_lo, z_lo: 4
    * x_hi, z_lo: 5
    * x_lo, z_hi: 6
    * x_hi, z_hi: 7
    * y_lo, z_lo: 8
    * y_hi, z_lo: 9
    * y_lo, z_hi: 10
    * y_hi, z_hi: 11
    *
    * node (codimension 3):
    * x_lo, y_lo, z_lo: 0
    * x_hi, y_lo, z_lo: 1
    * x_lo, y_hi, z_lo: 2
    * x_hi, y_hi, z_lo: 3
    * x_lo, y_lo, z_hi: 4
    * x_hi, y_lo, z_hi: 5
    * x_lo, y_hi, z_hi: 6
    * x_hi, y_hi, z_hi: 7
    *
    * \endverbatim
    *
    * @return The location index
    */
   int
   getLocationIndex() const
   {
      return d_location_index;
   }

   /*!
    * @brief Set the multiblock singularity flag to the argument value.
    *
    * In multiblock problems, the code setting up multiblock hierarchies and
    * patch levels will set this value to true when creating a BoundaryBox
    * that represents a patch boundary at a multiblock singularity.
    *
    * @param[in] is_mblk_singularity
    */
   void
   setIsMultiblockSingularity(
      bool is_mblk_singularity)
   {
      d_is_mblk_singularity = is_mblk_singularity;
   }

   /*!
    * @brief Get the value of the multiblock singularity flag.
    *
    * @return The singularity flag value.
    */
   bool
   getIsMultiblockSingularity() const
   {
      return d_is_mblk_singularity;
   }

   /*!
    * @brief The assignment operator copies all data components.
    *
    * @param[in] boundary_box
    */
   BoundaryBox&
   operator = (
      const BoundaryBox& boundary_box)
   {
      d_box = boundary_box.d_box;
      d_bdry_type = boundary_box.d_bdry_type;
      d_location_index = boundary_box.d_location_index;
      d_is_mblk_singularity = boundary_box.d_is_mblk_singularity;
      return *this;
   }

   /*!
    * @brief Enumerated type BoundaryOrientation is used to indicate where a
    * boundary box is located relative to a patch in a particular coordinate
    * direction.  MIDDLE means a boundary box is neither on the upper or lower
    * side of a patch in the given coordinate direction.  For example, an edge
    * boundary box on the right side of a patch in 2d is neither on the upper
    * or lower side of the patch in the J coordinate direction, so its
    * BoundaryOrientation value would be MIDDLE.  The same boundary box
    * would be UPPER in the I coordinate direction.
    */
   enum BoundaryOrientation {
      LOWER = -1,
      MIDDLE = 0,
      UPPER = 1
   };

   /*!
    * @brief Get which side of a patch the boundary box is on.
    *
    * Returns BoundaryOrientation value indicating whether the boundary
    * box is on the upper or lower side of the patch in the given coordinate
    * direction, or in the middle (neither upper nor lower).
    *
    * @return BoundaryOrientation value LOWER, MIDDLE, or UPPER
    *
    * @param[in] dir Coordinate direction on which to query
    *
    * @pre dir < getDim().getValue()
    */
   BoundaryOrientation
   getBoundaryOrientation(
      const int dir) const;

   /*!
    * @brief Return the dimension of this object.
    *
    * @return The dimension
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   /*!
    * @brief Dimension of the object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Box holding spatial location of the BoundaryBox.
    */
   Box d_box;

   /*!
    * @brief Codimension of the boundary.
    */
   int d_bdry_type;

   /*!
    * @brief Location index identifying relative location of the BoundaryBox
    */
   int d_location_index;

   /*!
    * @brief Flag telling whether the BoundaryBox is located at a
    * multiblock singularity.
    */
   bool d_is_mblk_singularity;
};

}
}

#endif
