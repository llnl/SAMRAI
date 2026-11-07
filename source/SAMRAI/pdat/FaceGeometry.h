/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceGeometry
#define included_pdat_FaceGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace pdat {

class FaceIterator;

/*!
 * Class FaceGeometry manages the mapping between the AMR index space
 * and the face-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between face-
 * centered box geometries for communication operations.
 *
 * See header file for FaceData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see FaceOverlap
 */

class FaceGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef FaceOverlap Overlap;

   /*!
    * @brief Convert an AMR index box space box into an face geometry box.
    * An face geometry box extends the given AMR index box space box
    * by one in upper end for the face normal coordinate direction.
    *
    * Recall that box indices are cyclically shifted such that the face normal
    * direction is the first coordinate index.  See SideData header file.
    *
    * @pre (face_normal >= 0) && (face_normal < box.getDim().getValue())
    */
   static hier::Box
   toFaceBox(
      const hier::Box& box,
      tbox::Dimension::dir_t face_normal);

   /*!
    * @brief Transform a face-centered box.
    *
    * At input, the  Box is assumed to represent an face-centered box with
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
    * @brief Transform a FaceIndex.
    *
    * This static method applies a coordinate system transformation to the
    * given FaceIndex.
    *
    * @param[in,out]  index
    * @param[in]      transformation
    */
   static void
   transform(
      FaceIndex& index,
      const hier::Transformation& transformation);

   static FaceIterator
   begin(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   static FaceIterator
   end(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   /*!
    * @brief Construct the face geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   FaceGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~FaceGeometry();

   /*!
    * @brief Compute the overlap in face-centered index space between
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
    * @brief Compute the face-centered destination boxes that represent
    * the overlap between the source box geometry and the destination
    * box geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset.getDim()
    */
   void
   computeDestinationBoxes(
      std::vector<hier::BoxContainer>& dst_boxes,
      const FaceGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes = hier::BoxContainer()) const;

   /*!
    * @brief Set up a FaceOverlap object based on the given boxes and the
    * transformation.
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this face centered box geometry
    * object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this face centered box
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
    * objects are guaranteed to have face centered geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const FaceGeometry& dst_geometry,
      const FaceGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   static void
   rotateAboutAxis(
      FaceIndex& index,
      const int axis,
      const int num_rotations);

   FaceGeometry(
      const FaceGeometry&);             // not implemented
   FaceGeometry&
   operator = (
      const FaceGeometry&);                     // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
