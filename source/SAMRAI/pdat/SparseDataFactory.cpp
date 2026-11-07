/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation for SparseDataFactory
 *
 ************************************************************************/
#ifndef included_pdat_SparseDataFactory_C
#define included_pdat_SparseDataFactory_C

#include "SAMRAI/pdat/SparseDataFactory.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/SparseData.h"
#include "SAMRAI/tbox/MemoryUtilities.h"

#include <memory>

namespace SAMRAI {
namespace pdat {

/*
 * C'tors and d'tors.
 */
template<typename BOX_GEOMETRY>
SparseDataFactory<BOX_GEOMETRY>::SparseDataFactory(
   const hier::IntVector& ghosts,
   const std::vector<std::string>& dbl_attributes,
   const std::vector<std::string>& int_attributes):
   hier::PatchDataFactory(ghosts),
   d_dbl_attributes(dbl_attributes),
   d_int_attributes(int_attributes)
{
}

template<typename BOX_GEOMETRY>
SparseDataFactory<BOX_GEOMETRY>::~SparseDataFactory()
{
}

/*
 * Implementation of base class pure virtual functions
 */
template<typename BOX_GEOMETRY>
std::shared_ptr<hier::PatchDataFactory>
SparseDataFactory<BOX_GEOMETRY>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);
   return std::make_shared<SparseDataFactory<BOX_GEOMETRY> >(
             ghosts,
             d_dbl_attributes,
             d_int_attributes);
}

template<typename BOX_GEOMETRY>
std::shared_ptr<hier::PatchData>
SparseDataFactory<BOX_GEOMETRY>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);
   return std::make_shared<SparseData<BOX_GEOMETRY> >(
             patch.getBox(),
             d_ghosts,
             d_dbl_attributes,
             d_int_attributes);
}

template<typename BOX_GEOMETRY>
std::shared_ptr<hier::BoxGeometry>
SparseDataFactory<BOX_GEOMETRY>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   return std::make_shared<BOX_GEOMETRY>(box, d_ghosts);
}

template<typename BOX_GEOMETRY>
size_t
SparseDataFactory<BOX_GEOMETRY>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);
   NULL_USE(box);
   return tbox::MemoryUtilities::align(
      sizeof(SparseData<BOX_GEOMETRY>));
}

template<typename BOX_GEOMETRY>
bool
SparseDataFactory<BOX_GEOMETRY>::validCopyTo(
   const std::shared_ptr<PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);
   bool valid_copy = false;

   if (!valid_copy) {

      std::shared_ptr<SparseDataFactory<BOX_GEOMETRY> > idf(
         std::dynamic_pointer_cast<SparseDataFactory<BOX_GEOMETRY>,
                                     hier::PatchDataFactory>(dst_pdf));

      if (idf) {
         valid_copy = true;
      }
   }
   return valid_copy;
}

} // end namespace pdat
} // end namespace SAMRAI

#endif
