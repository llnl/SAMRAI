/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating edge data objects
 *
 ************************************************************************/

#ifndef included_pdat_EdgeDataFactory_C
#define included_pdat_EdgeDataFactory_C

#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/MemoryUtilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the default ghost cell width and depth.
 *
 *************************************************************************
 */

template<class TYPE>
EdgeDataFactory<TYPE>::EdgeDataFactory(
   int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_has_allocator(false)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
EdgeDataFactory<TYPE>::EdgeDataFactory(
   int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var,
   tbox::ResourceAllocator allocator):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_allocator(allocator),
   d_has_allocator(true)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
EdgeDataFactory<TYPE>::~EdgeDataFactory()
{
}

/*
 *************************************************************************
 *
 * Clone the factory and copy the default parameters to the new factory.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchDataFactory>
EdgeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);
   if (d_has_allocator) {   
      return std::make_shared<EdgeDataFactory<TYPE> >(
               d_depth,
               ghosts,
               d_fine_boundary_represents_var,
               d_allocator);
   } else {          
      return std::make_shared<EdgeDataFactory<TYPE> >(
             d_depth,
             ghosts,
             d_fine_boundary_represents_var);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete edge data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
EdgeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   if (d_has_allocator) {
      return std::make_shared<EdgeData<TYPE> >(
             patch.getBox(),
             d_depth,
             d_ghosts,
             d_allocator);
   } else {
      return std::make_shared<EdgeData<TYPE> >(
             patch.getBox(),
             d_depth,
             d_ghosts);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for edge data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
EdgeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   return std::make_shared<EdgeGeometry>(box, d_ghosts);
}

template<class TYPE>
int
EdgeDataFactory<TYPE>::getDepth() const
{
   return d_depth;
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the data object.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
EdgeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(EdgeData<TYPE>));
   const size_t data =
      EdgeData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from EdgeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
EdgeDataFactory<TYPE>::validCopyTo(
   const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are EdgeData and OuteredgeData.
    */
   if (!valid_copy) {
      std::shared_ptr<EdgeDataFactory<TYPE> > edf(
         std::dynamic_pointer_cast<EdgeDataFactory<TYPE>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (edf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      std::shared_ptr<OuteredgeDataFactory<TYPE> > oedf(
         std::dynamic_pointer_cast<OuteredgeDataFactory<TYPE>,
                                     hier::PatchDataFactory>(
            dst_pdf));
      if (oedf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean value indicating how data for the edge quantity will be
 * treated on coarse-fine interfaces.  This value is passed into the
 * constructor.
 *
 *************************************************************************
 */
template<class TYPE>
bool
EdgeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return d_fine_boundary_represents_var;
}

/**
 * Return true since the edge data index space extends beyond the interior of
 * patches.  That is, edge data lives on patch borders.
 */
template<class TYPE>
bool
EdgeDataFactory<TYPE>::dataLivesOnPatchBorder() const {
   return true;
}

}
}
#endif
