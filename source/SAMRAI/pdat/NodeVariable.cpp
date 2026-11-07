/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeVariable_C
#define included_pdat_NodeVariable_C

#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for node variable objects
 *
 *************************************************************************
 */

template<class TYPE>
NodeVariable<TYPE>::NodeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<NodeDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),

   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
NodeVariable<TYPE>::NodeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<NodeDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     allocator)),

   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
NodeVariable<TYPE>::~NodeVariable()
{
}

template<class TYPE>
int NodeVariable<TYPE>::getDepth() const
{
   std::shared_ptr<NodeDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<NodeDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
