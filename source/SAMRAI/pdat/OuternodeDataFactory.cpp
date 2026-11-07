/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outernode data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeDataFactory_C
#define included_pdat_OuternodeDataFactory_C

#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/pdat/OuternodeData.h"
#include "SAMRAI/pdat/OuternodeGeometry.h"
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
OuternodeDataFactory<TYPE>::OuternodeDataFactory(
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
OuternodeDataFactory<TYPE>::OuternodeDataFactory(
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
OuternodeDataFactory<TYPE>::~OuternodeDataFactory()
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
OuternodeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   if (d_has_allocator) {
      return std::make_shared<OuternodeDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth,
             d_allocator);
   } else {
      return std::make_shared<OuternodeDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete outernode data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
OuternodeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   if (d_has_allocator) {
      return std::make_shared<OuternodeData<TYPE> >(patch.getBox(), d_depth, d_allocator);
   } else {
      return std::make_shared<OuternodeData<TYPE> >(patch.getBox(), d_depth);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outernode data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
OuternodeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   return std::make_shared<OuternodeGeometry>(box, zero_vector);
}

template<class TYPE>
int
OuternodeDataFactory<TYPE>::getDepth() const
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
OuternodeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuternodeData<TYPE>));
   const size_t data = OuternodeData<TYPE>::getSizeOfData(box,
         d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from OuternodeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
OuternodeDataFactory<TYPE>::validCopyTo(
   const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are NodeData and OuternodeData.
    */
   if (!valid_copy) {
      std::shared_ptr<NodeDataFactory<TYPE> > ndf(
         std::dynamic_pointer_cast<NodeDataFactory<TYPE>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (ndf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      std::shared_ptr<OuternodeDataFactory<TYPE> > ondf(
         std::dynamic_pointer_cast<OuternodeDataFactory<TYPE>,
                                     hier::PatchDataFactory>(
            dst_pdf));
      if (ondf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that fine data for the outernode
 * quantity will take precedence on coarse-fine interfaces.  See the
 * OuternodeVariable class header file for more information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuternodeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return true since the outernode data index space extends beyond the
 * interior of patches.  That is, outernode data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuternodeDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
}

}
}
#endif
