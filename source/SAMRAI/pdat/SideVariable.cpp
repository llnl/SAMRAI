/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideVariable_C
#define included_pdat_SideVariable_C

#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/SideDataFactory.h"
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
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   const hier::IntVector& directions,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     directions)),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(directions)
{
   TBOX_ASSERT(directions.getDim() == getDim());
}

template<class TYPE>
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     hier::IntVector::getOne(dim))),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(hier::IntVector::getOne(dim))
{
}

template<class TYPE>
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   const hier::IntVector& directions,
   tbox::ResourceAllocator allocator,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     directions,
                     allocator)),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(directions)
{
   TBOX_ASSERT(directions.getDim() == getDim());
}

template<class TYPE>
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  std::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var,
                     hier::IntVector::getOne(dim),
                     allocator)),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(hier::IntVector::getOne(dim))
{
}


template<class TYPE>
SideVariable<TYPE>::~SideVariable()
{
}

template<class TYPE>
const hier::IntVector& SideVariable<TYPE>::getDirectionVector() const
{
   return d_directions;
}

template<class TYPE>
int SideVariable<TYPE>::getDepth() const
{
   std::shared_ptr<SideDataFactory<TYPE> > factory(
      SAMRAI_SHARED_PTR_CAST<SideDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
