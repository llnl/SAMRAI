/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideGeometry
#define included_pdat_SideGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace pdat {

class SideIterator;

/*!
 * Class SideGeometry manages the mapping between the AMR index space
 * and the side-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between side-
 * centered box geometries for communication operations.
 *
 * See header file for SideData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see SideOverlap
 */

class SideGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef SideOverlap Overlap;

   /*!
    * @brief Convert an AMR index box space box into an side geometry box.
    * An side geometry box extends the given AMR index space box
    * by one at upper end for the side normal coordinate direction.
    *
    * @pre (side_normal >= 0) && (side_normal < box.getDim().getValue())
    */
   static hier::Box
   toSideBox(
      const hier::Box& box,
      tbox::Dimension::dir_t side_normal);

   /*!
    * @brief Transform a side-centered box.
    *
    * At input, the  Box is assumed to represent an side-centered box with
    * the given normal direction.  The Box will be transformed according to
    * the coordinate system transformation defined by the Transformation
    * object.  At output, normal_direction will represent the output box's
    * new normal direction.
    *
    * @param[in,out]  box
    * @param[in,out]  normal_direction
    * @param[in]      transformation
    */
   static void
   transform(
      hier::Box& box,
      int& normal_direction,
      const hier::Transformation& transformation);

   /*!
    * @brief Transform a SideIndex.
    *
    * This static method applies a coordinate system transformation to the
    * given SideIndex.
    *
    * @param[in,out]  index
    * @param[in]      transformation
    */
   static void
   transform(
      SideIndex& index,
      const hier::Transformation& transformation);

   static SideIterator
   begin(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   static SideIterator
   end(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   /*!
    * @brief Construct the side geometry object given an AMR index
    * space box, ghost cell width and directions vector indicating
    * which coordinate directions are allocated.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    * @pre directions.min() >= 0
    */
   SideGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts,
      const hier::IntVector& directions);

   /*!
    * @brief Construct the side geometry object given an AMR index
    * space box and ghost cell width.
    *
    * No directions vector is provided, so it is assumed that all
    * coordinate directions are allocated.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   SideGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~SideGeometry();

   /*!
    * @brief Compute the overlap in side-centered index space between
    * the source box geometry and the destination box geometry.
    *
    * @pre getBox().getDim() == src_mask.getDim()
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   calculateOverlap(
      const hier::BoxGeometry& dst_geometry,
      const hier::BoxGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const bool retry,
      const hier::BoxContainer& dst_restrict_boxes = hier::BoxContainer()) const;

   /*!
    * @brief Compute the side-centered destination boxes that represent
    * the overlap between the source box geometry and the destination
    * box geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset.getDim()
    */
   void
   computeDestinationBoxes(
      std::vector<hier::BoxContainer>& dst_boxes,
      const SideGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes = hier::BoxContainer()) const;

   /*!
    * @brief Set up a SideOverlap object based on the given boxes and the
    * transformation.
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this side centered box geometry
    * object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this side centered box
    * geometry object.
    */
   const hier::IntVector&
   getGhosts() const
   {
      return d_ghosts;
   }

   /*!
    * Return constant reference to vector describing which coordinate
    * directions managed by this side geometry object.
    *
    * A vector entry of zero indicates that this object will not perform
    * operations involving the corresponding coordinate direction.
    * A non-zero value indicates otherwise.
    */
   const hier::IntVector&
   getDirectionVector() const
   {
      return d_directions;
   }

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have side centered geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset().getDim()
    * @pre dst_geometry.getDirectionVector() == src_geometry.getDirectionVector()
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const SideGeometry& dst_geometry,
      const SideGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   static void
   rotateAboutAxis(
      SideIndex& index,
      const tbox::Dimension::dir_t axis,
      const int num_rotations);

   SideGeometry(
      const SideGeometry&);             // not implemented
   SideGeometry&
   operator = (
      const SideGeometry&);                     // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;
   hier::IntVector d_directions;

};

}
}

#endif
