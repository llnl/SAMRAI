/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating face data objects
 *
 ************************************************************************/

#ifndef included_pdat_FaceDataFactory_C
#define included_pdat_FaceDataFactory_C

#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/hier/Patch.h"


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
FaceDataFactory<TYPE>::FaceDataFactory(
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

template <class TYPE>
FaceDataFactory<TYPE>::FaceDataFactory(int depth,
                                       const hier::IntVector& ghosts,
                                       bool fine_boundary_represents_var,
                                       tbox::ResourceAllocator allocator)
    : hier::PatchDataFactory(ghosts),
      d_depth(depth),
      d_fine_boundary_represents_var(fine_boundary_represents_var),
      d_allocator(allocator),
      d_has_allocator(true)
{
  TBOX_ASSERT(depth > 0);
  TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
FaceDataFactory<TYPE>::~FaceDataFactory()
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
FaceDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   if (d_has_allocator) {
      return std::make_shared<FaceDataFactory>(d_depth,
                                           ghosts,
                                           d_fine_boundary_represents_var,
                                           d_allocator);
   } else {
      return std::make_shared<FaceDataFactory>(
             d_depth,
             ghosts,
             d_fine_boundary_represents_var);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete face data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
FaceDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   if (d_has_allocator) {
      return std::make_shared<FaceData<TYPE> >(patch.getBox(), d_depth, d_ghosts, d_allocator);
   } else {
      return std::make_shared<FaceData<TYPE> >(
             patch.getBox(),
             d_depth,
             d_ghosts);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for face data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
FaceDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   return std::make_shared<FaceGeometry>(box, d_ghosts);
}

template<class TYPE>
int
FaceDataFactory<TYPE>::getDepth() const
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
FaceDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(FaceData<TYPE>));
   const size_t data =
      FaceData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from FaceData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
FaceDataFactory<TYPE>::validCopyTo(
   const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are FaceData and OuterfaceData.
    */
   if (!valid_copy) {
      std::shared_ptr<FaceDataFactory<TYPE> > fdf(
         std::dynamic_pointer_cast<FaceDataFactory<TYPE>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (fdf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      std::shared_ptr<OuterfaceDataFactory<TYPE> > ofdf(
         std::dynamic_pointer_cast<OuterfaceDataFactory<TYPE>,
                                     hier::PatchDataFactory>(
            dst_pdf));
      if (ofdf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 * Return a boolean value indicating how data for the face quantity will be
 * treated on coarse-fine interfaces.  This value is passed into the
 * constructor.  See the FaceVariable<TYPE> class header file for more
 * information.
 */
template<class TYPE>
bool
FaceDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return d_fine_boundary_represents_var;
}

/*
 * Return true since the face data index space extends beyond the interior of
 * patches.  That is, face data lives on patch borders.
 */
template<class TYPE>
bool
FaceDataFactory<TYPE>::dataLivesOnPatchBorder() const {
   return true;
}

}
}
#endif
