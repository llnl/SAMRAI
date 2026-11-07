/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating side data objects
 *
 ************************************************************************/

#ifndef included_pdat_SideDataFactory
#define included_pdat_SideDataFactory

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
 * Class SideDataFactory is a factory class used to allocate new
 * instances of SideData objects.  It is a subclass of the patch
 * data factory class and side data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * Note that it is possible to create a side data factory to allocate
 * and manage data for cell sides associated with a single coordinate
 * direction only.  See the constructor for more information.
 *
 * @see SideData
 * @see PatchDataFactory
 */

template<class TYPE>
class SideDataFactory:public hier::PatchDataFactory
{
public:
   /*!
    * @brief The constructor for the side data factory class.
    *
    * The ghost cell width, depth (number of components), and fine boundary
    * representation arguments give the defaults for all edge data objects
    * created with this factory.
    *
    * The directions vector describes the coordinate directions for which
    * data will be allocated on the sides of cells in the grid. A value of
    * 1 indicates that data will be allocated for that coordinate direction,
    * while a value of zero means data for that direction is not wanted.
    * To allocate data in all directions, provide an IntVector that is
    * 1 in all directions, or use the other constructor, which assumes by
    * default that all directions are desired.  See the
    * SideVariable<TYPE> class header file for more information.
    *
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    * @pre directions.min() >= 0
    */
   SideDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var,
      const hier::IntVector& directions);

   /**
    * Constructor for the cell data factory class that takes a specific
    * tbox::ResourceAllocator.  The ghost cell width and depth (number of components)
    * arguments give the defaults for all cell data objects created with this
    * factory.
    *
    * @pre depth > 0 @pre ghosts.min() >= 0
    */
   SideDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var,
      const hier::IntVector& directions,
      tbox::ResourceAllocator allocator);

   /*!
    * @brief Constructor for the side data factory class setting up allocation
    * of data in all coordinate directions.
    *
    * This constructor works the same as the other constructor, but
    * it takes no directions argument, meaning that all directions are going
    * to be allocated.
    *
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   SideDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var);

   /*!
    * @brief Constructor for the side data factory class setting up allocation
    * of data in all coordinate directions and uses an Umpire allocator
    *
    * This constructor works the same as the other constructor, but
    * it takes no directions argument, meaning that all directions are going
    * to be allocated.
    *
    * @pre depth > 0
    * @pre ghosts.min() >= 0
    */
   SideDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var,
      tbox::ResourceAllocator allocator);

   /**
    * Virtual destructor for the side data factory class.
    */
   virtual ~SideDataFactory();

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
    * Virtual factory function to allocate a concrete side data object.
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
    * will be used in the instantiation of side data objects.
    */
   int
   getDepth() const;

   /**
    * Return constant reference to vector describing which coordinate
    * directions have data associated with this side data object.
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction.  A non-zero
    * value indicates that a valid data array is maintained for that
    * coordinate direction.
    */
   const hier::IntVector&
   getDirectionVector() const;

   /**
    * Calculate the amount of memory needed to store the side data object,
    * including object data and dynamically allocated data.
    *
    * @pre getDim() == box.getDim()
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean value indicating how data for the side quantity will be
    * treated on coarse-fine interfaces.  This value is passed into the
    * constructor.  See the FaceVariable<TYPE> class header file for more
    * information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the side data index space extends beyond the interior
    * of patches.  That is, side data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this SideDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is SideDataFactory or OutersideDataFactory, false otherwise.
    *
    * @pre getDim() == dst_pdf->getDim()
    */
   bool
   validCopyTo(
      const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

   /*!
    * @brief Return true if this factory has an Umpire Allocator.
    */
   bool hasAllocator() const
   {
      return d_has_allocator;
   }

   /*!
    * @brief Get the Umpire Allocator.
    */
   tbox::ResourceAllocator getAllocator() const
   {
      return d_allocator;
   }

private:
   int d_depth;
   bool d_fine_boundary_represents_var;
   hier::IntVector d_directions;

   tbox::ResourceAllocator d_allocator;
   bool d_has_allocator;
};

} // Namespace pdat
} // Namespace SAMRAI

#include "SAMRAI/pdat/SideDataFactory.cpp"

#endif
