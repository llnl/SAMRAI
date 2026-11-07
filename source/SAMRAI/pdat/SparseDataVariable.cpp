/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SparseDataVariable
 *
 ************************************************************************/

#ifndef included_pdat_SparseDataVariable_C
#define included_pdat_SparseDataVariable_C

#include "SAMRAI/pdat/SparseDataVariable.h"


#include "SAMRAI/pdat/SparseDataFactory.h"

namespace SAMRAI {
namespace pdat {

template<typename BOX_GEOMETRY>
SparseDataVariable<BOX_GEOMETRY>::SparseDataVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   const std::vector<std::string>& dbl_attributes,
   const std::vector<std::string>& int_attributes):
   hier::Variable(name,
                  std::make_shared<SparseDataFactory<BOX_GEOMETRY> >(
                     hier::IntVector::getZero(dim),
                     dbl_attributes, int_attributes))
{
}

template<typename BOX_GEOMETRY>
SparseDataVariable<BOX_GEOMETRY>::~SparseDataVariable()
{
}

} // end namespace pdat
} // end namespace SAMRAI
#endif
