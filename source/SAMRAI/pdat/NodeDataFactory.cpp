/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating node data objects
 *
 ************************************************************************/

#ifndef included_pdat_NodeDataFactory_C
#define included_pdat_NodeDataFactory_C

#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
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
NodeDataFactory<TYPE>::NodeDataFactory(
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
NodeDataFactory<TYPE>::NodeDataFactory(
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
NodeDataFactory<TYPE>::~NodeDataFactory()
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
NodeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   if (d_has_allocator) {
      return
         std::make_shared<NodeDataFactory<TYPE> >(
            d_depth,
            ghosts,
            d_fine_boundary_represents_var,
            d_allocator);
   } else {
      return
         std::make_shared<NodeDataFactory<TYPE> >(
            d_depth,
            ghosts,
            d_fine_boundary_represents_var);
   }
}

/*
 *************************************************************************
 *
 * Allocate the concrete node data classes.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::PatchData>
NodeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);
   if (d_has_allocator) {
      return std::make_shared<NodeData<TYPE> >(
         patch.getBox(),
         d_depth,
         d_ghosts,
         d_allocator);
   } else {
      return std::make_shared<NodeData<TYPE> >(
         patch.getBox(),
         d_depth,
         d_ghosts);
   }
}

/*
 *************************************************************************
 *
 * Return the box geometry type for node data objects.
 *
 *************************************************************************
 */

template<class TYPE>
std::shared_ptr<hier::BoxGeometry>
NodeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   return std::make_shared<NodeGeometry>(box, d_ghosts);
}

template<class TYPE>
int
NodeDataFactory<TYPE>::getDepth() const
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
NodeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(NodeData<TYPE>));
   const size_t data =
      NodeData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from NodeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
NodeDataFactory<TYPE>::validCopyTo(
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
 * Return a boolean value indicating how data for the node quantity will be
 * treated on coarse-fine interfaces.  This value is passed into the
 * constructor.  See the NodeVariable<TYPE> class header file for more
 * information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
NodeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return d_fine_boundary_represents_var;
}

/*
 *************************************************************************
 *
 * Return true since the node data index space extends beyond the interior
 * of patches.  That is, node data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
NodeDataFactory<TYPE>::dataLivesOnPatchBorder() const {
   return true;
}

}
}
#endif
