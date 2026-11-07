/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outernode data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeDataFactory
#define included_pdat_OuternodeDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/ResourceAllocator.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief
 * Class OuternodeDataFactory is a factory class used to allocate new
 * instances of OuternodeData objects.  It is a subclass of the patch
 * data factory class and outernode data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see OuternodeData
 * @see PatchDataFactory
 */

template<class TYPE>
class OuternodeDataFactory:public hier::PatchDataFactory
{
public:
   /*!
    * @brief The constructor for the outernode data factory class.
    *
    * The depth (number of components) gives the default for all of
    * the outernode data objects created with this factory.
    *
    * @pre depth > 0
    */
   OuternodeDataFactory(
      const tbox::Dimension& dim,
      int depth);

   /**
    * The constructor for the outernode data factory class.
    * The depth (number of components) sets the default for all of
    * the outernode data objects created with this factory.
    *
    * This constructor sets an Umpire allocator for the memory management
    * of the data held within outernode data objects.
    *
    * @pre depth > 0
    */
   OuternodeDataFactory(
      const tbox::Dimension& dim,
      int depth,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Virtual destructor for the outernode data factory class.
    */
   virtual ~OuternodeDataFactory();

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

   /*!
    * @brief Virtual factory function to allocate a concrete outernode data
    *        object.
    *
    * The default information about the object (e.g., depth) is taken from
    * the factory.
    *
    * @pre getDim() == patch.getDim()
    */
   virtual std::shared_ptr<hier::PatchData>
   allocate(
      const hier::Patch& patch) const;

   /*!
    * @brief Allocate the box geometry object associated with the patch data.
    *
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    *
    * @pre getDim() == box.getDim()
    */
   virtual std::shared_ptr<hier::BoxGeometry>
   getBoxGeometry(
      const hier::Box& box) const;

   /*!
    * @brief Get the depth (number of components).
    *
    * This is the depth that will be used in the instantiation of
    * outernode data objects.
    */
   int
   getDepth() const;

   /*!
    * @brief Calculate the amount of memory needed to store the outernode data
    *        object, including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean true value indicating that fine data for the outernode
    * quantity will take precedence on coarse-fine interfaces.  See the
    * OuternodeVariable class header file for more information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the outernode data index space extends beyond the
    * interior of patches.  That is, outernode data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this OuternodeDataFactory to the
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
   hier::IntVector d_no_ghosts;
   tbox::ResourceAllocator d_allocator;
   bool d_has_allocator;

};

}
}

#include "SAMRAI/pdat/OuternodeDataFactory.cpp"

#endif
