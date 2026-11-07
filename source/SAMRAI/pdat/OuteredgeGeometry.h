/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Box geometry information for edge centered objects
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeGeometry
#define included_pdat_OuteredgeGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

class EdgeGeometry;

/*!
 * Class OuteredgeGeometry manages the mapping between the AMR index
 * and the outeredge geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between outeredge
 * box geometries and edge or outeredge box geometries for communication
 * operations.
 *
 * See header file for OuteredgeData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see EdgeGeometry
 * @see EdgeOverlap
 */

class OuteredgeGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef EdgeOverlap Overlap;

   /*!
    * Convert a given box in the standard cell-centered AMR index space to an
    * outeredge geometry box for the specified axis, face normal, and
    * lower/upper side.   See OuteredgeData header file for a detailed
    * description of an outeredge box.
    *
    * @pre (0 <= axis) && (axis < dim.getValue())
    * @pre (0 <= face_normal) && (face_normal < dim.getValue())
    * @pre (face_normal != axis)
    * @pre (side == 0) || (side == 1)
    */
   static hier::Box
   toOuteredgeBox(
      const hier::Box& box,
      tbox::Dimension::dir_t axis,
      tbox::Dimension::dir_t face_normal,
      int side);

   /*!
    * @brief Construct an outeredge geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   OuteredgeGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OuteredgeGeometry();

   /*!
    * @brief Compute the overlap in edge-centered index space on the
    * boundaries of the source box geometry and the destination box geometry.
    *
    * @pre getBox().getDim() == src_mask.getDim()
    * @pre
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
    * @brief Set up a EdgeOverlap object based on the given boxes and the
    * transformation.
    */
   virtual std::shared_ptr<hier::BoxOverlap>
   setUpOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation) const;

   /*!
    * @brief Return the box for this outeredge box geometry object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this outeredge box geometry object.
    */
   const hier::IntVector&
   getGhosts() const
   {
      return d_ghosts;
   }

private:
   /*!
    * Compute overlap between a source outeredge geometry and a destination
    * edge geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const EdgeGeometry& dst_geometry,
      const OuteredgeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   /*!
    * Compute overlap between a source outeredge geometry and a destination
    * outeredge geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const OuteredgeGeometry& dst_geometry,
      const OuteredgeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes);

   OuteredgeGeometry(
      const OuteredgeGeometry&);               // not implemented
   OuteredgeGeometry&
   operator = (
      const OuteredgeGeometry&);                // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
