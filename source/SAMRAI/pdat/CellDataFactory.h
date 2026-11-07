/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating cell data objects
 *
 ************************************************************************/

#ifndef included_pdat_CellDataFactory
#define included_pdat_CellDataFactory

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
 * Class CellDataFactory is a factory class used to allocate new
 * instances of CellData objects.  It is a subclass of the patch
 * data factory class and cell data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see CellData
 * @see PatchDataFactory
 */

template<class TYPE>
class CellDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The default constructor for the cell data factory class.  The ghost
    * cell width and depth (number of components) arguments give the defaults
    * for all cell data objects created with this factory.
    *
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   CellDataFactory(
      int depth,
      const hier::IntVector& ghosts);

   /**
    * Constructor for the cell data factory class that takes a specific
    * tbox::ResourceAllocator.  The ghost cell width and depth (number of components)
    * arguments give the defaults for all cell data objects created with this
    * factory.
    *
    * @pre depth > 0 @pre ghosts.min() >= 0
    */
   CellDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      tbox::ResourceAllocator allocator);

   /**
    * Virtual destructor for the cell data factory class.
    */
   virtual ~CellDataFactory();

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
    * Virtual factory function to allocate a concrete cell data object.
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
    * will be used in the instantiation of cell data objects.
    */
   int
   getDepth() const;

   /**
    * Calculate the amount of memory needed to store the cell data object,
    * including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean true value indicating that the cell data quantities will
    * always be treated as though fine values represent them on coarse-fine
    * interfaces.
    * See the CellVariable<TYPE> class header file for more information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return false since the cell data index space matches the cell-centered
    * index space for AMR patches.  Thus, cell data does not live on patch
    * borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this CellDataFactory to the
    * supplied destination patch data factory. It will return true if
    * dst_pdf is a CellDataFactory, false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   int d_depth;

   tbox::ResourceAllocator d_allocator;
   bool d_has_allocator;
};

}
}

#include "SAMRAI/pdat/CellDataFactory.cpp"

#endif
