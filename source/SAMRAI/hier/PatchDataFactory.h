/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory abstract base class for creating patch data objects
 *
 ************************************************************************/

#ifndef included_hier_PatchDataFactory
#define included_hier_PatchDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchData.h"

#include <memory>

namespace SAMRAI {
namespace hier {

class Patch;
class MultiblockDataTranslator;

/**
 * Class PatchDataFactory is an abstract base class used to allocate
 * new instances of patch data objects.  Recall that patch data objects (PDs)
 * are the data storage containers that exist within a patch.  PDs are
 * created using patch data factory (PDF) objects; this is an example of
 * the ``Abstract Factory'' method described in the Design Patterns book
 * by Gamma, et al.
 *
 * The separation of PDF from PD simplifies the creation of new concrete PD
 * classes since it separates the definition of the concrete class type from
 * the actual instantiation.  The actual concrete class associated with the
 * PDF is unknown to most of the framework; the PDF only defines enough
 * information to create the PD instance.  For example, to add a new type
 * of PD object MyPD (MyPatchData):
 * - Derive MyPDF from PDF and implement the abstract virtual
 *       function calls as appropriate for the concrete subclass;
 *       in particular, the allocate() function will return an instance
 *       of MyPD.
 * - Derive MyPD from PD and implement the abstract virtual
 *       function calls as appropriate for the concrete subclass.
 * - When defining the types of storage needed for a patch, add
 *       MyPDF to the patch descriptor list.
 * - Now whenever the PDF base class of MyPDF is asked to create
 *       a concrete class, it will create MyPD.
 * The creation of concrete PD objects is managed through the allocate()
 * interfaces in PDF.
 *
 * In addition to the generation of patch data, the patch data factory
 * also generates box geometry descriptions used to calculate the overlap
 * between two patch data objects.  The allocation of the box geometry
 * object is managed by the patch data factory instead of the patch data
 * object since patch data factories are guaranteed to exist on all of the
 * processors independent of the mapping of patches to processors.  Patch
 * data is guaranteed to exist only on those patches local to a processor.
 *
 * @see BoxGeometry
 * @see PatchData
 * @see PatchDescriptor
 */

class PatchDataFactory
{
public:
   /**
    * The constructor for the patch data factory class.
    *
    * @param ghosts ghost cell width for concrete classes created from
    * the factory.
    *
    * @pre ghosts.min() >= 0
    */
   explicit PatchDataFactory(
      const IntVector& ghosts);

   /**
    * @brief Virtual destructor for the patch data factory class.
    *
    */
   virtual ~PatchDataFactory() = default;

   /**
    * @brief Abstract virtual function to clone a patch data factory.
    *
    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts ghost cell width for concrete classes created from
    * the factory.
    */
   virtual std::shared_ptr<PatchDataFactory>
   cloneFactory(
      const IntVector& ghosts) = 0;

   /**
    * @brief Abstract virtual function to allocate a concrete patch data object.
    *
    */
   virtual std::shared_ptr<PatchData>
   allocate(
      const Patch& patch) const = 0;

   /**
    * @brief
    * Abstract virtual function to allocate a concrete box geometry
    * object.
    *
    * The box geometry object will be used in the calculation
    * of box intersections for the computation of data dependencies.
    */
   virtual std::shared_ptr<BoxGeometry>
   getBoxGeometry(
      const Box& box) const = 0;

   /**
    * @brief Get the ghost cell width.
    *
    * This is the ghost cell width that will be used in the
    * instantiation of concrete patch data instances.  The
    * ghost width is specified in the clone method.
    */
   const IntVector&
   getGhostCellWidth() const
   {
      return d_ghosts;
   }

   /**
    * @brief Abstract virtual function to compute the amount of memory needed to
    * allocate for object data and to represent the object itself.
    *
    * This includes any dynamic storage, such as arrays, needed by the
    * concrete patch data instance.  Although the patch data subclass
    * may choose not to allocate memory it must not use more memory
    * than requested here.
    */
   virtual size_t
   getSizeOfMemory(
      const Box& box) const = 0;

   /**
    * @brief Return true if the fine data values represent the data quantity
    * on coarse-fine interfaces if data lives on patch borders; false
    * otherwise.
    *
    * The boolean return value is supplied by the concrete
    * patch data factory subclass.
    */
   virtual bool
   fineBoundaryRepresentsVariable() const = 0;

   /**
    * @brief Return true if the variable data lives on patch borders; false otherwise.
    *
    * The boolean return value is supplied by the concrete patch data factory subclass.
    */
   virtual bool
   dataLivesOnPatchBorder() const = 0;

   /**
    * @brief Abstract virtual function that returns whether the
    * current PatchDataFactory can be copied to the supplied
    * destination PatchDataFactory.
    *
    * Mechanisms to check for valid types are implemented in the patch
    * data factory subclasses for particular datatypes.
    */
   virtual bool
   validCopyTo(
      const std::shared_ptr<PatchDataFactory>& dst_pdf) const = 0;

   virtual MultiblockDataTranslator *
   getMultiblockDataTranslator();

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_ghosts.getDim();
   }

protected:
   IntVector d_ghosts;

private:
   PatchDataFactory(
      const PatchDataFactory&) = delete;
   PatchDataFactory&
   operator = (
      const PatchDataFactory&) = delete;
   PatchDataFactory() = delete;  // must specify ghost width

};

}
}

#endif
