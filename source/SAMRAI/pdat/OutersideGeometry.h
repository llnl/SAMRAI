/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OutersideGeometry
#define included_pdat_OutersideGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

class SideGeometry;

/*!
 * Class OutersideGeometry manages the mapping between the AMR index
 * and the outerside geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between outerside
 * box geometries and side or outerside box geometries for communication
 * operations.
 *
 * See header file for OutersideData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see SideGeometry
 * @see SideOverlap
 */

class OutersideGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef SideOverlap Overlap;

   /*!
    * @brief Construct an outerside geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim
    * @pre ghosts.min() >= 0
    */
   OutersideGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OutersideGeometry();

   /*!
    * @brief Compute the overlap in side-centered index space on the
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
    * @brief Set up a SideOverlap object based on the given boxes and the
    * offset.
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this outerside box geometry object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this outerside box geometry object.
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
    * has outerside geometry and the destination side geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset().getDim()
    * @pre dst_geometry.getDirectionVector() == hier::IntVector::getOne(src_mask.getDim())
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const SideGeometry& dst_geometry,
      const OutersideGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const OutersideGeometry& dst_geometry,
      const OutersideGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   OutersideGeometry(
      const OutersideGeometry&);                // not implemented
   OutersideGeometry&
   operator = (
      const OutersideGeometry&);                    // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
