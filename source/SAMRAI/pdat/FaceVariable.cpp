/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceVariable_C
#define included_pdat_FaceVariable_C

#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for face variable objects
 *
 *************************************************************************
 */

template<class TYPE>
FaceVariable<TYPE>::FaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   const bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<FaceDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
FaceVariable<TYPE>::FaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth,
   const bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<FaceDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     allocator)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
FaceVariable<TYPE>::~FaceVariable()
{
}

template<class TYPE>
int FaceVariable<TYPE>::getDepth() const
{
   std::shared_ptr<FaceDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<FaceDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
