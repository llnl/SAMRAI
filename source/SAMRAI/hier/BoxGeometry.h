/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Box geometry description for overlap computations
 *
 ************************************************************************/

#ifndef included_hier_BoxGeometry
#define included_hier_BoxGeometry

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/**
 * Class BoxGeometry encapsulates the geometry information associated
 * with data defined over a box region, such as its ghost cell width and the
 * centering or geometry of the data index space.  The intersection (or overlap)
 * of two box geometries is generated via member function calculateOverlap().  The
 * form of this overlap depends on the particular geometry represented by the
 * subclass.
 *
 * Box geometry objects are created by the patch data factories since patch
 * data objects may not be available for patches that are distributed across
 * processor memories (patch data factories are always replicated).
 *
 * The concept of ``overlap'' or data dependency is more complex for generic
 * box geometry objects than for just cell-centered box indices in the abstract
 * AMR index space.  Problems arise in cases where data lies on the outside
 * corners, faces, or edges of a box.  For these data types, it is likely
 * that there will exist duplicate data values on different patches.
 *
 * The solution implemented here introduces the concept of ``priority''
 * between patches.  Data of patches with higher priority can overwrite
 * the interiors (face, node, or edge values associated with cells that
 * constitute the interior of the patch) of patches with lower priorities,
 * but lower priority patches can never overwrite the interiors of higher
 * priority patches.  This scheme introduces a total ordering of data and
 * therefore eliminates the duplicate information problem.
 *
 * In practice, this protocol means two things: (1) the communication
 * routines must always process copies from low priority sources to high
 * priority sources, and (2) patches must be given special permission to
 * overwrite their interior values during a write.  All destinations are
 * therefore represented by three quantities: (1) the box geometry of the
 * destination (which encodes the box, ghost cells, and geometry), (2) the
 * box geometry of the source, and (3) a flag indicating whether the source
 * has a higher priority than the destination (that is, whether the source
 * can overwrite the interior of the destination).  If the overwrite flag is
 * set, then data will be copied over the specified box domain and may write
 * into the interior of the destination.  If the overwrite flag is not set,
 * then data will be copied only into the ghost cell values and not the
 * interior values of the patch.
 *
 * @see BoxOverlap
 * @see PatchDataFactory
 * @see PatchData
 */

class BoxGeometry
{
public:
   /**
    * The default constructor for BoxGeometry does nothing interesting.
    */
   BoxGeometry() = default;

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~BoxGeometry() = default;

   /**
    * Calculate the overlap between two box geometry objects given the
    * source and destination (given by this) geometries, a source mask,
    * the priority overwrite flag, and a transformation from  source to
    * destination index spaces.  The box overlap description returned by
    * this function will be used in later copy and pack/unpack calls on
    * the patch data object.   The transformation is from the source space
    * into the destination index space.  That is, if p is in the source
    * index space, then p after being transformed  is the corresponding
    * point in the destination index space.  The overwrite flag is used to
    * represent priority between patches.  If it is set, then the copy is
    * allowed to modify the interior of the destination region.  Note that the
    * source and destination box geometries encode the geometry of the box
    * that they represent; thus, it is possible to calculate intersections
    * between different geometries.  This will be necessary when copying
    * data from flux sum counters into a face centered array in the AMR
    * flux synchronization algorithm.  The optional argument
    * dst_restrict_boxes can be used to add a further restriction on the
    * calculated overlap, so that none of the calculated overlap will lie
    * outside of the space covered by dst_restrict_boxes. If dst_restrict_boxes
    * is an empty BoxContainer, then it will have no effect on the overlap
    * calculation.
    */
   std::shared_ptr<BoxOverlap>
   calculateOverlap(
      const BoxGeometry& src_geometry,
      const Box& src_mask,
      const Box& fill_box,
      const bool overwrite_interior,
      const Transformation& transformation,
      const BoxContainer& dst_restrict_boxes = BoxContainer()) const
   {
      return this->calculateOverlap(
         *this,
         src_geometry,
         src_mask,
         fill_box,
         overwrite_interior,
         transformation,
         true,
         dst_restrict_boxes);
   }

   /**
    * Calculate the overlap between two box geometry objects given the
    * source and destination geometries.  This form calculateOverlap() is
    * redefined by the subclasses of BoxGeometry for the appropriate
    * intersection algorithms.  If calculateOverlap() cannot compute the
    * intersection between the two given geometries and retry is true, then
    * calculateOverlap() is called on the destination geometry object with
    * retry set to false (to avoid infinite recursion).  This protocol
    * makes it possible to add new box geometry types and still calculate
    * intersections with existing box geometry types.  The optional argument
    * dst_restrict_boxes can be used to add a further restriction on the
    * calculated overlap, so that none of the calculated overlap will lie
    * outside of the space covered by dst_restrict_boxes. If dst_restrict_boxes
    * is an empty BoxContainer, then it will have no effect on the overlap
    * calculation.
    */
   virtual std::shared_ptr<BoxOverlap>
   calculateOverlap(
      const BoxGeometry& dst_geometry,
      const BoxGeometry& src_geometry,
      const Box& src_mask,
      const Box& fill_box,
      const bool overwrite_interior,
      const Transformation& src_offset,
      const bool retry,
      const BoxContainer& dst_restrict_boxes = BoxContainer()) const = 0;

   /**
    * Set up a BoxOverlap object that consists simply of the given boxes
    * and the transformation.
    */
   virtual std::shared_ptr<BoxOverlap>
   setUpOverlap(
      const BoxContainer& boxes,
      const Transformation& offset) const = 0;

private:
   BoxGeometry(
      const BoxGeometry&) = delete;
   BoxGeometry&
   operator = (
      const BoxGeometry&) = delete;

};

}
}

#endif
