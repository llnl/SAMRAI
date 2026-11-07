/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_EdgeVariable_C
#define included_pdat_EdgeVariable_C

#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for edge variable objects
 *
 *************************************************************************
 */

template<class TYPE>
EdgeVariable<TYPE>::EdgeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   const bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<EdgeDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
EdgeVariable<TYPE>::~EdgeVariable()
{
}

template<class TYPE>
int EdgeVariable<TYPE>::getDepth() const
{
   std::shared_ptr<EdgeDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<EdgeDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
