/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_EdgeGeometry
#define included_pdat_EdgeGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace pdat {

class EdgeIterator;

/*!
 * Class EdgeGeometry manages the mapping between the AMR index space
 * and the edge-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between edge-
 * centered box geometries for communication operations.
 *
 * See header file for EdgeData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see EdgeOverlap
 */

class EdgeGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef EdgeOverlap Overlap;

   /*!
    * @brief Convert an AMR index box space box into an edge geometry box.
    * An edge geometry box extends the given AMR index box space box
    * by one in upper end for each coordinate direction not equal
    * to the axis direction.
    *
    * @pre (0 <= axis) && (axis < box.getDim().getValue())
    */
   static hier::Box
   toEdgeBox(
      const hier::Box& box,
      int axis_direction);

   /*!
    * @brief Transform an edge-centered box.
    *
    * At input, the  Box is assumed to represent an edge-centered box with
    * the given axis direction.  The Box will be transformed according to
    * the coordinate system transformation defined by the Transformation
    * object.  At output, axis_direction will represent the output box's
    * new axis direction.
    *
    * @param[in,out]  box
    * @param[in,out]  axis_direction
    * @param[in]      transformation
    */
   static void
   transform(
      hier::Box& box,
      int& axis_direction,
      const hier::Transformation& transformation);

   /*!
    * @brief Transform an EdgeIndex.
    *
    * This static method applies a coordinate system transformation to the
    * given EdgeIndex.
    *
    * @param[in,out]  index
    * @param[in]      transformation
    */
   static void
   transform(
      EdgeIndex& index,
      const hier::Transformation& transformation);

   static EdgeIterator
   begin(
      const hier::Box& box,
      int axis);

   static EdgeIterator
   end(
      const hier::Box& box,
      int axis);

   /*!
    * @brief Construct the edge geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   EdgeGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~EdgeGeometry();

   /*!
    * @brief Compute the overlap in edge-centered index space between
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
    * @brief Compute the edge-centered destination boxes that represent
    * the overlap between the source box geometry and the destination
    * box geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset.getDim()
    */
   void
   computeDestinationBoxes(
      std::vector<hier::BoxContainer>& dst_boxes,
      const EdgeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes = hier::BoxContainer()) const;

   /*!
    * @brief Set up a EdgeOverlap object based on the given boxes and the
    * transformation
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this edge centered box geometry
    * object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this edge centered box
    * geometry object.
    */
   const hier::IntVector&
   getGhosts() const
   {
      return d_ghosts;
   }

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have edge centered geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const EdgeGeometry& dst_geometry,
      const EdgeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   static void
   rotateAboutAxis(
      EdgeIndex& index,
      const int axis,
      const int num_rotations);

   EdgeGeometry(
      const EdgeGeometry&);             // not implemented
   EdgeGeometry&
   operator = (
      const EdgeGeometry&);                     // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
