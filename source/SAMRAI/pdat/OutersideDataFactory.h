/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outerside data objects
 *
 ************************************************************************/

#ifndef included_pdat_OutersideDataFactory
#define included_pdat_OutersideDataFactory

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
 * Class OutersideDataFactory is a factory class used to allocate new
 * instances of OutersideData objects.  It is a subclass of the patch
 * data factory class and outerside data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see OutersideData
 * @see PatchDataFactory
 */

template<class TYPE>
class OutersideDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The constructor for the outerside data factory class.
    * The depth (number of components) sets the default for all of
    * the outerside data objects created with this factory.
    *
    * @pre depth > 0
    */
   OutersideDataFactory(
      const tbox::Dimension& dim,
      const int depth);

   /**
    * The constructor for the outerside data factory class.
    * The depth (number of components) sets the default for all of
    * the outerside data objects created with this factory.
    *
    * This constructor sets an Umpire allocator for the management of
    * data held in the outserside data objects.
    *
    * @pre depth > 0
    */
   OutersideDataFactory(
      const tbox::Dimension& dim,
      const int depth,
      tbox::ResourceAllocator allocator);

   /**
    * Virtual destructor for the outerside data factory class.
    */
   virtual ~OutersideDataFactory();

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
    * Virtual factory function to allocate a concrete outerside data object.
    * The default information about the object (e.g., depth) is taken from
    * the factory.
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
    * will be used in the instantiation of outerside data objects.
    */
   int
   getDepth() const;

   /**
    * Calculate the amount of memory needed to store the outerside data
    * object, including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean true value indicating that fine data for the outerside
    * quantity will take precedence on coarse-fine interfaces.  See the
    * OutersideVariable<TYPE> class header file for more information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the outerside data index space extends beyond the
    * interior of patches.  That is, outerside data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this OutersideDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is SideDataFactory or OutersideDataFactory, false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   int d_depth;
   hier::IntVector d_no_ghosts;
   tbox::ResourceAllocator d_allocator;
   bool d_has_allocator;
};

}
}

#include "SAMRAI/pdat/OutersideDataFactory.cpp"

#endif
