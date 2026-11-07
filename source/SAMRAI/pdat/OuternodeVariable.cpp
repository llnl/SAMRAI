/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outernode centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeVariable_C
#define included_pdat_OuternodeVariable_C

#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
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
OuternodeVariable<TYPE>::OuternodeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  std::make_shared<OuternodeDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuternodeVariable<TYPE>::~OuternodeVariable()
{
}

template<class TYPE>
int OuternodeVariable<TYPE>::getDepth() const
{
   std::shared_ptr<OuternodeDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<OuternodeDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
