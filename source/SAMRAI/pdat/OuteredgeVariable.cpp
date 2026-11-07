/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outeredge centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeVariable_C
#define included_pdat_OuteredgeVariable_C

#include "SAMRAI/pdat/OuteredgeVariable.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
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
OuteredgeVariable<TYPE>::OuteredgeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  std::make_shared<OuteredgeDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuteredgeVariable<TYPE>::~OuteredgeVariable()
{
}

template<class TYPE>
int OuteredgeVariable<TYPE>::getDepth() const
{
   std::shared_ptr<OuteredgeDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<OuteredgeDataFactory<TYPE>, hier::PatchData>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
