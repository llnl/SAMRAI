/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   IndexDataFactory
 *
 ************************************************************************/

#ifndef included_pdat_IndexDataFactory
#define included_pdat_IndexDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataFactory.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/**
 * Class IndexDataFactory<TYPE, BOX_GEOMETRY> is the factory class used to
 * allocate new instances of IndexData<TYPE, BOX_GEOMETRY> objects.  It is a
 * subclass of the hier::PatchDataFactory base class.
 *
 * @see IndexData
 * @see IndexVariable
 * @see hier::PatchDataFactory
 */

template<class TYPE, class BOX_GEOMETRY>
class IndexDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The default constructor for the IndexDataFactory<TYPE, BOX_GEOMETRY>
    * class.  The ghost cell width argument gives the default width for all
    * irregular data objects created with this factory.
    */
   explicit IndexDataFactory(
      const hier::IntVector& ghosts);

   /**
    * Virtual destructor for the irregular data factory class.
    */
   virtual ~IndexDataFactory();

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
    * Virtual factory function to allocate a concrete index data object.
    * The default information about the object (e.g., ghost cell width) is
    * taken from the factory.
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
    * Calculate the amount of memory needed to store the irregular data
    * object, including object data but not dynamically allocated data.
    * Because the irregular data list can grow and shrink, it would be
    * impossible to estimate the necessary amount of memory.  Instead,
    * dynamic data is allocated via the standard new/free mechanisms.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean true value indicating that the index data quantities
    * will always be treated as though fine values represent them on
    * coarse-fine interfaces. See the IndexVariable<TYPE, BOX_GEOMETRY> class
    * header file for more information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return false since the index data index space matches the cell-centered
    * index space for AMR patches.  Thus, index data does not live on patch
    * borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this IndexDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is an IndexDataFactory of the same type and dimension,
    * false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   IndexDataFactory(
      const IndexDataFactory&);
   IndexDataFactory&
   operator = (
      const IndexDataFactory&);

};

}
}

#include "SAMRAI/pdat/IndexDataFactory.cpp"

#endif
