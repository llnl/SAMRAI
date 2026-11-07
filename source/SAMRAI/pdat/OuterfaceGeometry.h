/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceGeometry
#define included_pdat_OuterfaceGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

class FaceGeometry;

/*!
 * Class OuterfaceGeometry manages the mapping between the AMR index
 * and the outerface geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between outerface
 * box geometries and face or outerface box geometries for communication
 * operations.
 *
 * See header file for OuterfaceData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see FaceGeometry
 * @see FaceOverlap
 */

class OuterfaceGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef FaceOverlap Overlap;

   /*!
    * @brief Construct an outerface geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   OuterfaceGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OuterfaceGeometry();

   /*!
    * @brief Compute the overlap in face-centered index space on the
    * boundaries of the source box geometry and the destination box geometry.
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
    * @brief Set up a FaceOverlap object based on the given boxes and the
    * transformation
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this outerface box geometry object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this outerface box geometry object.
    */
   const hier::IntVector&
   getGhosts() const
   {
      return d_ghosts;
   }

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where the source
    * has outerface geometry and the destination face geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const FaceGeometry& dst_geometry,
      const OuterfaceGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const OuterfaceGeometry& dst_geometry,
      const OuterfaceGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   OuterfaceGeometry(
      const OuterfaceGeometry&);                // not implemented
   OuterfaceGeometry&
   operator = (
      const OuterfaceGeometry&);                    // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
