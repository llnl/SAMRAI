/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   IndexDataFactory implementation
 *
 ************************************************************************/

#ifndef included_pdat_IndexDataFactory_C
#define included_pdat_IndexDataFactory_C

#include "SAMRAI/pdat/IndexDataFactory.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MemoryUtilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the default ghost cell width.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
IndexDataFactory<TYPE, BOX_GEOMETRY>::IndexDataFactory(
   const hier::IntVector& ghosts):
   hier::PatchDataFactory(ghosts)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataFactory<TYPE, BOX_GEOMETRY>::~IndexDataFactory()
{
}

/*
 *************************************************************************
 *
 * Clone the factory and copy the default parameters to the new factory.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
std::shared_ptr<hier::PatchDataFactory>
IndexDataFactory<TYPE, BOX_GEOMETRY>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   return std::make_shared<IndexDataFactory<TYPE, BOX_GEOMETRY> >(ghosts);
}

/*
 *************************************************************************
 *
 * Allocate the concrete irregular data class.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
std::shared_ptr<hier::PatchData>
IndexDataFactory<TYPE, BOX_GEOMETRY>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   return std::make_shared<IndexData<TYPE, BOX_GEOMETRY> >(
             patch.getBox(),
             d_ghosts);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for index data objects.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
std::shared_ptr<hier::BoxGeometry>
IndexDataFactory<TYPE, BOX_GEOMETRY>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   return std::make_shared<BOX_GEOMETRY>(box, d_ghosts);
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the object.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
size_t
IndexDataFactory<TYPE, BOX_GEOMETRY>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   NULL_USE(box);
   return tbox::MemoryUtilities::align(sizeof(IndexData<TYPE, BOX_GEOMETRY>));
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from IndexData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
bool
IndexDataFactory<TYPE, BOX_GEOMETRY>::validCopyTo(
   const std::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid option is another IndexData object of the same dimension
    * and type.
    */
   if (!valid_copy) {
      std::shared_ptr<IndexDataFactory<TYPE, BOX_GEOMETRY> > idf(
         std::dynamic_pointer_cast<IndexDataFactory<TYPE, BOX_GEOMETRY>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (idf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that the index data quantities
 * will always be treated as though fine values represent them on
 * coarse-fine interfaces. See the IndexVariable<TYPE, BOX_GEOMETRY> class
 * header file for more information.
 *
 *************************************************************************
 */
template<class TYPE, class BOX_GEOMETRY>
bool
IndexDataFactory<TYPE, BOX_GEOMETRY>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return false since the index data index space matches the cell-centered
 * index space for AMR patches.  Thus, index data does not live on patch
 * borders.
 *
 *************************************************************************
 */
template<class TYPE, class BOX_GEOMETRY>
bool
IndexDataFactory<TYPE, BOX_GEOMETRY>::dataLivesOnPatchBorder() const
{
   return false;
}

}
}
#endif
