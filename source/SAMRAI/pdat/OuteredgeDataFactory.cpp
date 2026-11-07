/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outeredge data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeDataFactory_C
#define included_pdat_OuteredgeDataFactory_C

#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/pdat/OuteredgeGeometry.h"
#include "SAMRAI/hier/Patch.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the depth of the patch data.
 *
 *************************************************************************
 */

template<class TYPE>
OuteredgeDataFactory<TYPE>::OuteredgeDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim)),
   d_has_allocator(false)
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuteredgeDataFactory<TYPE>::OuteredgeDataFactory(
   const tbox::Dimension& dim,
   int depth,
   tbox::ResourceAllocator allocator):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim)),
   d_allocator(allocator),
   d_has_allocator(true)
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuteredgeDataFactory<TYPE>::~OuteredgeDataFactory()
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
OuteredgeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   if (d_has_allocator) {
      return std::make_shared<OuteredgeDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth,
             d_allocator);
   } else {
      return std::make_shared<OuteredgeDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete outeredge data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
OuteredgeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   if (d_has_allocator) {
      return std::make_shared<OuteredgeData<TYPE> >(patch.getBox(), d_depth, d_allocator);
   } else {
      return std::make_shared<OuteredgeData<TYPE> >(patch.getBox(), d_depth);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outeredge data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
OuteredgeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   return std::make_shared<OuteredgeGeometry>(box, zero_vector);
}

template<class TYPE>
int
OuteredgeDataFactory<TYPE>::getDepth() const
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
OuteredgeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuteredgeData<TYPE>));
   const size_t data = OuteredgeData<TYPE>::getSizeOfData(box,
         d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that fine data for the outeredge
 * quantity will take precedence on coarse-fine interfaces.  See the
 * OuteredgeVariable<TYPE> class header file for more information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuteredgeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return true;
}

/*
 *************************************************************************
 *
 * Return true since the outeredge data index space extends beyond the
 * interior of patches.  That is, outeredge data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuteredgeDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
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
OuteredgeDataFactory<TYPE>::validCopyTo(
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

}
}
#endif
