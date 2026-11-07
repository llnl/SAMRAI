/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_CellVariable_C
#define included_pdat_CellVariable_C

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for cell variable objects
 *
 *************************************************************************
 */

template<class TYPE>
CellVariable<TYPE>::CellVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  std::make_shared<CellDataFactory<TYPE> >(depth,
                                                             hier::IntVector::getZero(dim))) // default zero ghost cells
{
}

template<class TYPE>
CellVariable<TYPE>::CellVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   tbox::ResourceAllocator allocator,
   int depth):
   hier::Variable(name,
                  std::make_shared<CellDataFactory<TYPE> >(depth,
                                                           hier::IntVector::getZero(dim),
                                                           allocator)) // default zero ghost cells
{
}

template<class TYPE>
CellVariable<TYPE>::~CellVariable()
{
}

template<class TYPE>
int
CellVariable<TYPE>::getDepth() const
{
   std::shared_ptr<CellDataFactory<TYPE> > cell_factory(
      SAMRAI_SHARED_PTR_CAST<CellDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(cell_factory);
   return cell_factory->getDepth();
}

/*
 *************************************************************************
 *
 * Return true indicating that cell data quantities will always
 * be treated as though fine values take precedence on coarse-fine
 * interfaces.  Note that this is really artificial since the cell
 * data index space matches the cell-centered index space for AMR
 * patches.  However, some value must be supplied for communication
 * operations.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellVariable<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return false indicating that cell data on a patch interior
 * does not exist on the patch boundary.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellVariable<TYPE>::dataLivesOnPatchBorder() const
{
   return false;
}

}
}
#endif
