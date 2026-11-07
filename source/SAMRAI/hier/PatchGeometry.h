/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches
 *
 ************************************************************************/

#ifndef included_hier_PatchGeometry
#define included_hier_PatchGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/LocalId.h"
#include "SAMRAI/hier/PatchBoundaries.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <list>
#include <vector>

namespace SAMRAI {
namespace hier {

/**
 * Class PatchGeometry is the base class for geometry classes that
 * manage index space and mesh increment information on individual patches.
 * Patch geometry information is used for setting boundary conditions in
 * ghost cells and is used in the inter-level transfer operators for refining
 * or coarsening data between two patches associated with different index
 * spaces.  The boundary information for patches is actually computed by
 * the BaseGridGeometry class.
 *
 * @see BoundaryBox
 * @see BaseGridGeometry
 */

class PatchGeometry
{
public:
   /*!
    * @brief Array of 2*DIM booleans (with default constructor),
    * used to instantiate the sparse container map<LocalId,TwoDimBool>
    * (map<LocalId,bool[2*DIM]> does not work).
    */
   class TwoDimBool
   {
public:
      explicit TwoDimBool(
         const tbox::Dimension& dim);

      TwoDimBool(
         const tbox::Dimension& dim,
         bool v);

      void
      setAll(
         bool v)
      {
         for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
            d_data[i] = v;
         }
      }

      bool&
      operator () (
         int dim,
         int side)
      {
         TBOX_ASSERT(dim >= 0 && dim < d_dim.getValue());
         TBOX_ASSERT(side == 0 || side == 1);
         return d_data[2 * dim + side];
      }

      const bool&
      operator () (
         int dim,
         int side) const
      {
         TBOX_ASSERT(dim >= 0 && dim < d_dim.getValue());
         TBOX_ASSERT(side == 0 || side == 1);
         return d_data[2 * dim + side];
      }

      /**
       * Return the dimension of this object.
       */
      const tbox::Dimension&
      getDim() const
      {
         return d_dim;
      }

private:
      /*
       * Unimplemented default constructor.
       */
      TwoDimBool();

      const tbox::Dimension d_dim;
      bool d_data[2 * SAMRAI::MAX_DIM_VAL];
   };

   /**
    * The constructor for the patch geometry base class.
    *
    * @pre (ratio_to_level_zero.getDim() == touches_regular_bdry.getDim())
    * @pre all components of ratio_to_level_zero must be nonzero and all
    *      components of ratio_to_level_zero not equal to 1 must have the same
    *      sign
    */
   PatchGeometry(
      const IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const BlockId& block_id);

   /**
    * The virtual destructor for the patch geometry base class.
    */
   virtual ~PatchGeometry() = default;

   /**
    * Return const reference to patch boundary information.
    */
   const std::vector<std::vector<BoundaryBox> >&
   getPatchBoundaries() const
   {
      return d_patch_boundaries.getVectors();
   }

   /*!
    * @brief Set the boundary box vectors for this patch geometry.
    *
    * A vector of length DIM of std::vector<BoundaryBox> is passed
    * in to be stored as the boundary boxes for this patch geometry.
    *
    * @param bdry The vector of BoundaryBox vectors.
    */
   void
   setBoundaryBoxesOnPatch(
      const std::vector<std::vector<BoundaryBox> >& bdry);

   /**
    * Return const reference to ratio to level zero index space.
    */
   const IntVector&
   getRatio() const
   {
      return d_ratio_to_level_zero;
   }

   /**
    * Return a boolean value indicating whether the patch boundary
    * intersects the physical domain boundary in a non-periodic
    * direction.  In other words, the return value is true when the
    * patch has non-empty boundary boxes that lie outside the physical
    * domain.  Otherwise, the return value is false.  Note that when
    * a patch touches the "boundary" of the physical domain in a periodic
    * direction, there are no boundary boxes to fill; the data is filled
    * from the proper region of the domain interior in the periodic direction.
    */
   bool
   intersectsPhysicalBoundary() const
   {
      return d_has_regular_boundary;
   }

   /**
    * Return vector of boundary box components for patch each of which
    * intersects the patch at a single point (i.e., 0-dim intersection
    * between cells in patch and cells in boundary box).
    */
   const std::vector<BoundaryBox>&
   getNodeBoundaries() const
   {
      return d_patch_boundaries[getDim().getValue() - 1];
   }

   /**
    * Return vector of boundary box components for patch each of which
    * intersects the patch along a 1-dim edge (i.e., 1-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * @pre getDim().getValue() >= 2
    */
   const std::vector<BoundaryBox>&
   getEdgeBoundaries() const
   {
      if (getDim().getValue() < 2) {
         TBOX_ERROR("PatchGeometry error in getEdgeBoundary...\n"
            << "DIM < 2 not supported." << std::endl);
      }

      // The "funny" indexing prevents a warning when compiling for
      // DIM < 2.  This code is only reached if DIM >= 2 when
      // executing.
      return d_patch_boundaries[getDim().getValue() < 2 ? 0 : getDim().getValue() - 2];
   }

   /**
    * Return vector of boundary box components for patch each of which
    * intersects the patch along a 2-dim face (i.e., 2-dim intersection
    * between cells in patch and cells in boundary box).
    *
    * @pre getDim().getValue() >= 3
    */
   const std::vector<BoundaryBox>&
   getFaceBoundaries() const
   {
      if (getDim().getValue() < 3) {
         TBOX_ERROR("PatchGeometry error in getFaceBoundary...\n"
            << "DIM < 3 not supported." << std::endl);
      }

      // The "funny" indexing prevents a warning when compiling for
      // DIM < 3.  This code is only reached if DIM >= 3 when
      // executing.
      return d_patch_boundaries[getDim().getValue() < 3 ? 0 : getDim().getValue() - 3];
   }

   /**
    * Return vector of boundary box components for patch each of which
    * intersects the patch as a (DIM - codim)-dimensional object.
    * That is,
    *
    * if DIM == 1: (codim == 1) => same components as getNodeBoundaries.
    *
    * if DIM == 2, (codim == 1) => same components as getEdgeBoundaries.
    *              (codim == 2) => same components as getNodeBoundaries.
    *
    * if DIM == 3, (codim == 1) => same components as getFaceBoundaries.
    *              (codim == 2) => same components as getEdgeBoundaries.
    *              (codim == 3) => same components as getNodeBoundaries.
    *
    * @pre (codim > 0) && (codim <= getDim().getValue())
    * when codim < 0 or codim > DIM.
    */
   const std::vector<BoundaryBox>&
   getCodimensionBoundaries(
      const int codim) const
   {
      TBOX_ASSERT((codim > 0) && (codim <= getDim().getValue()));
      return d_patch_boundaries[codim - 1];
   }

   /**
    * Set the vector of boundary box components of the given codimension
    * for a patch.
    *
    * @pre (codim > 0) && (codim <= getDim().getValue())
    * @pre for each boundary_box in bdry_boxes, getBoundaryType() == codim
    */
   void
   setCodimensionBoundaries(
      const std::vector<BoundaryBox>& bdry_boxes,
      const int codim);

   /*!
    * @brief Compute a box outside a physical domain that needs to be filled.
    *
    * The patch box will be grown by the given ghost cell width and
    * then intersected with the boundary box.  The resulting intersection
    * will be grown to the needed ghost cell width in the direction
    * normal to the boundary.
    *
    * @param bbox BoundaryBox representing location and type of boundary
    * @param patch_box The box for the patch where data is being filled
    * @param gcw ghost cell width to fill
    *
    * @pre bbox.getDim() == patch_box.getDim()
    * @pre bbox.getDim() == gcw.getDim()
    * @pre all components of gcw are >= 0
    */
   Box
   getBoundaryFillBox(
      const BoundaryBox& bbox,
      const Box& patch_box,
      const IntVector& gcw) const;

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any non-periodic physical boundary
    */
   bool
   getTouchesRegularBoundary() const
   {
      return d_has_regular_boundary;
   }

   /*!
    * @brief Query whether patch touches a regular boundary
    *
    * Returns true if the Patch touches any periodic boundary
    */
   bool
   getTouchesPeriodicBoundary() const
   {
      return d_has_periodic_boundary;
   }

   /*!
    * @brief Query whether patch touches a specific regular boundary
    *
    * Returns true if the Patch touches a non-periodic physical boundary
    * on the side of the Patch specified in the argument list.  The side
    * is specified by an axis direction and a flag specified the upper or
    * lower side.
    *
    * @param axis       Axis direction normal to the side being checked
    * @param upperlower Flag should be 0 if checking the lower side in the
    *                   axis direction, or 1 if checking the upper side.
    *
    * @pre (axis >= 0) && (axis < getDim().getValue())
    * @pre (upperlower == 0) || (upperlower == 1)
    */
   bool
   getTouchesRegularBoundary(
      int axis,
      int upperlower) const
   {
      TBOX_ASSERT(axis >= 0 && axis < getDim().getValue());
      TBOX_ASSERT(upperlower == 0 || upperlower == 1);
      return d_touches_regular_bdry(axis, upperlower);
   }

   /**
    * Print object data to the specified output stream.
    */
   void
   printClassData(
      std::ostream& stream) const;

   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   const tbox::Dimension d_dim;

   bool d_has_regular_boundary;
   bool d_has_periodic_boundary;
   IntVector d_ratio_to_level_zero;
   PatchBoundaries d_patch_boundaries;

   TwoDimBool d_touches_regular_bdry;
   BlockId d_block_id;
};

}
}

#endif
