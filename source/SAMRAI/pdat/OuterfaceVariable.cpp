/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceVariable_C
#define included_pdat_OuterfaceVariable_C

#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
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
OuterfaceVariable<TYPE>::OuterfaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  std::make_shared<OuterfaceDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuterfaceVariable<TYPE>::OuterfaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth):
   hier::Variable(name,
                  std::make_shared<OuterfaceDataFactory<TYPE> >(
                     dim, depth, allocator))
{
}

template<class TYPE>
OuterfaceVariable<TYPE>::~OuterfaceVariable()
{
}

template<class TYPE>
int OuterfaceVariable<TYPE>::getDepth() const
{
   std::shared_ptr<OuterfaceDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<OuterfaceDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
