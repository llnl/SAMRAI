/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OutersideVariable_C
#define included_pdat_OutersideVariable_C

#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/OutersideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for side variable objects
 *
 *************************************************************************
 */

template<class TYPE>
OutersideVariable<TYPE>::OutersideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  std::make_shared<OutersideDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OutersideVariable<TYPE>::OutersideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth):
   hier::Variable(name,
                  std::make_shared<OutersideDataFactory<TYPE> >(dim,
                                                                depth,
                                                                allocator))
{
}

template<class TYPE>
OutersideVariable<TYPE>::~OutersideVariable()
{
}

template<class TYPE>
int OutersideVariable<TYPE>::getDepth() const
{
   std::shared_ptr<OutersideDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<OutersideDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
