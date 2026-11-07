/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating node data objects
 *
 ************************************************************************/

#ifndef included_pdat_NodeDataFactory
#define included_pdat_NodeDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/**
 * Class NodeDataFactory is a factory class used to allocate new
 * instances of NodeData objects.  It is a subclass of the patch
 * data factory class and node data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see NodeData
 * @see PatchDataFactory
 */

template<class TYPE>
class NodeDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The constructor for the node data factory class. The ghost cell width,
    * depth (number of components), and fine boundary representation arguments
    * give the defaults for all edge data objects created with this factory.
    * See the NodeVariable<TYPE> class header file for more information.
    *
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   NodeDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var);

   /**
    * Constructor for the node data factory class that takes a specific
    * tbox::ResourceAllocator.  The ghost cell width and depth (number of components)
    * arguments give the defaults for all node data objects created with this
    * factory.
    *
    * @pre depth > 0 @pre ghosts.min() >= 0
    */
   NodeDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var,
      tbox::ResourceAllocator allocator);

   /**
    * Virtual destructor for the node data factory class.
    */
   virtual ~NodeDataFactory();

   /**
    * @brief Abstract virtual function to clone a patch data factory.
    *
    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts default ghost cell width for concrete classes created from
    * the factory.
    *
    * @pre getDim() == ghosts.getDim()
    */
   virtual std::shared_ptr<hier::PatchDataFactory>
   cloneFactory(
      const hier::IntVector& ghosts);

   /**
    * Virtual factory function to allocate a concrete node data object.
    * The default information about the object (e.g., ghost cell width)
    * is taken from the factory.
    *
    * @pre getDim() == patch.getDim()
    */
   virtual std::shared_ptr<hier::PatchData>
   allocate(
      const hier::Patch& patch) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    *
    * @pre getDim() == box.getDim()
    */

   virtual std::shared_ptr<hier::BoxGeometry>
   getBoxGeometry(
      const hier::Box& box) const;

   /**
    * Get the depth (number of components).  This is the depth that
    * will be used in the instantiation of node data objects.
    */
   int
   getDepth() const;

   /**
    * Calculate the amount of memory needed to store the node data object,
    * including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean value indicating how data for the node quantity will be
    * treated on coarse-fine interfaces.  This value is passed into the
    * constructor.  See the NodeVariable<TYPE> class header file for more
    * information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the node data index space extends beyond the interior
    * of patches.  That is, node data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this NodeDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is NodeDataFactory or OuternodeDataFactory, false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   int d_depth;
   bool d_fine_boundary_represents_var;

   tbox::ResourceAllocator d_allocator;
   bool d_has_allocator;

};

}
}

#include "SAMRAI/pdat/NodeDataFactory.cpp"

#endif
