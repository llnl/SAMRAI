/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeGeometry
#define included_pdat_NodeGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

class NodeIterator;

/*!
 * Class NodeGeometry manages the mapping between the AMR index space
 * and the node-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between node-
 * centered box geometries for communication operations.
 *
 * See header file for NodeData<TYPE> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see NodeOverlap
 */

class NodeGeometry:public hier::BoxGeometry
{
public:
   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef NodeOverlap Overlap;

   /*!
    * @brief Convert an AMR index box space box into a node geometry box.
    * A node geometry box is extends the given AMR index box space box
    * by one at upper end for each coordinate direction.
    */
   static hier::Box
   toNodeBox(
      const hier::Box& box)
   {
      return box.empty() ?
             box :
             hier::Box(box.lower(), box.upper() + 1, box.getBlockId());
   }

   /*!
    * @brief Transform a node-centered box.
    *
    * At input, the Box is assumed to represent an node-centered box.
    * The Box will be transformed according to the coordinate system
    * transformation defined by the Transformation object.
    *
    * @param[in,out]  box
    * @param[in]      transformation
    */
   static void
   transform(
      hier::Box& box,
      const hier::Transformation& transformation);

   static NodeIterator
   begin(
      const hier::Box& box);

   static NodeIterator
   end(
      const hier::Box& box);

   /*!
    * @brief Transform a NodeIndex.
    *
    * This static method applies a coordinate system transformation to the
    * given NodeIndex.
    *
    * @param[in,out]  index
    * @param[in]      transformation
    */
   static void
   transform(
      NodeIndex& index,
      const hier::Transformation& transformation);

   /*!
    * @brief Construct the node geometry object given an AMR index
    * space box and ghost cell width.
    *
    * @pre box.getDim() == ghosts.getDim()
    * @pre ghosts.min() >= 0
    */
   NodeGeometry(
      const hier::Box& box,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~NodeGeometry();

   /*!
    * @brief Compute the overlap in node-centered index space between
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
    * @brief Compute the node-centered destination boxes that represent
    * the overlap between the source box geometry and the destination
    * box geometry.
    *
    * @pre src_mask.getDim() == transformation.getOffset.getDim()
    */
   void
   computeDestinationBoxes(
      hier::BoxContainer& dst_boxes,
      const NodeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
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
    * @brief Return the box for this node centered box geometry
    * object.
    */
   const hier::Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Return the ghost cell width for this node centered box
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
    * objects are guaranteed to have node centered geometry.
    */
   static std::shared_ptr<hier::BoxOverlap>
   doOverlap(
      const NodeGeometry& dst_geometry,
      const NodeGeometry& src_geometry,
      const hier::Box& src_mask,
      const hier::Box& fill_box,
      const bool overwrite_interior,
      const hier::Transformation& transformation,
      const hier::BoxContainer& dst_restrict_boxes)
   {
      hier::BoxContainer dst_boxes;
      dst_geometry.computeDestinationBoxes(dst_boxes,
         src_geometry,
         src_mask,
         fill_box,
         overwrite_interior,
         transformation,
         dst_restrict_boxes);

      // Create the node overlap data object using the boxes and source shift
      return std::make_shared<NodeOverlap>(dst_boxes, transformation);
   }

   static void
   rotateAboutAxis(
      NodeIndex& index,
      const int axis,
      const int num_rotations);

   NodeGeometry(
      const NodeGeometry&);             // not implemented
   NodeGeometry&
   operator = (
      const NodeGeometry&);                     // not implemented

   hier::Box d_box;
   hier::IntVector d_ghosts;

};

}
}

#endif
