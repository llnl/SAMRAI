/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianPatchGeometry
#define included_geom_CartesianPatchGeometry

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianPatchGeometry implements simple Cartesian mes
 * geometry management for a single patch in an AMR hierarchy.  The geometry is
 * limited to the mesh increments being specified by the DIM-tuple
 * (dx[0],...,dx[DIM-1]) associated with the patch, and the spatial
 * coordinates of the lower and upper corners of the patch within the
 * computational domain.  The grid data is set by CartesianGridGeometry
 * class.  This patch geometry class is derived from hier::PatchGeometry
 * base class.
 *
 * @see hier::BoundaryBox
 * @see hier::PatchGeometry
 * @see CartesianGridGeometry
 */

class CartesianPatchGeometry:
   public hier::PatchGeometry
{
public:
   typedef hier::PatchGeometry::TwoDimBool TwoDimBool;

   /**
    * Constructor for CartesianPatchGeometry class.  It simply passes
    * patch boundary information to hier::PatchGeometry base class constructor
    * and allocates storage for spatial coordinates on patch.
    *
    * @pre dx != 0
    * @pre x_lo != 0
    * @pre x_up != 0
    */
   CartesianPatchGeometry(
      const hier::IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const hier::BlockId& block_id,
      const double * dx,
      const double * x_lo,
      const double * x_hi);

   /**
    * Destructor for CartesianPatchGeometry deallocates the
    * storage for spatial coordinates on patch.
    */
   virtual ~CartesianPatchGeometry();

   /**
    * Return const pointer to dx array for patch.
    */
   const double *
   getDx() const
   {
      return d_dx;
   }

   /**
    * Return const pointer to lower spatial coordinate for patch.
    */
   const double *
   getXLower() const
   {
      return d_x_lo;
   }

   /**
    * Return const pointer to upper spatial coordinate for patch.
    */
   const double *
   getXUpper() const
   {
      return d_x_up;
   }

   /*!
    * @brief Put coordinates in a database in the Conduit blueprint format.
    *
    * This stores a description of this geometry as uniform coordinates
    * in the blueprint format.
    *
    * @param coords_db   Database to hold the coordinate information
    * @param box         Box corresponding to the patch that holds this
    *                    geometry.
    */ 
   void
   putBlueprintCoords(
      const std::shared_ptr<tbox::Database>& coords_db,
      const hier::Box& box) const;

   /**
    * Print CartesianPatchGeometry class data.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

private:
   // These are not implemented.
   CartesianPatchGeometry(
      const CartesianPatchGeometry&);
   CartesianPatchGeometry&
   operator = (
      const CartesianPatchGeometry&);

   double d_dx[SAMRAI::MAX_DIM_VAL];   // mesh increments for patch.
   double d_x_lo[SAMRAI::MAX_DIM_VAL]; // spatial coords of lower end of patch.
   double d_x_up[SAMRAI::MAX_DIM_VAL]; // spatial coords of upper end of patch.

};

}
}

#endif
