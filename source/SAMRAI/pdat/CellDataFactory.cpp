/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating cell data objects
 *
 ************************************************************************/

#ifndef included_pdat_CellDataFactory_C
#define included_pdat_CellDataFactory_C

#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/hier/Patch.h"


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

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
CellDataFactory<TYPE>::CellDataFactory(
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_has_allocator(false)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
CellDataFactory<TYPE>::CellDataFactory(
   int depth,
   const hier::IntVector& ghosts,
   tbox::ResourceAllocator allocator):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_allocator(allocator),
   d_has_allocator(true)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
CellDataFactory<TYPE>::~CellDataFactory()
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
CellDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   if (d_has_allocator) {
      return
         std::make_shared<CellDataFactory<TYPE> >(d_depth, ghosts, d_allocator);
   } else {
      return 
         std::make_shared<CellDataFactory<TYPE> >(d_depth, ghosts);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete cell data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
CellDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   if (d_has_allocator) {
      return std::make_shared<CellData<TYPE> >(
               patch.getBox(),
               d_depth,
               d_ghosts,
               d_allocator);
   } else {
      return std::make_shared<CellData<TYPE> >(
             patch.getBox(),
             d_depth,
             d_ghosts);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for cell data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
CellDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   return std::make_shared<CellGeometry>(box, d_ghosts);
}

template<class TYPE>
int
CellDataFactory<TYPE>::getDepth() const
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
CellDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(CellData<TYPE>));
   const size_t data =
      CellData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from CellData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
CellDataFactory<TYPE>::validCopyTo(
   const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Only valid option is CellData.
    */
   std::shared_ptr<CellDataFactory<TYPE> > cdf(
      std::dynamic_pointer_cast<CellDataFactory<TYPE>,
                                  hier::PatchDataFactory>(dst_pdf));
   if (cdf) {
      valid_copy = true;
   }
   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that the cell data quantities will
 * always be treated as though fine values represent them on coarse-fine
 * interfaces.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return false since the cell data index space matches the cell-centered
 * index space for AMR patches.  Thus, cell data does not live on patch
 * borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return false;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
