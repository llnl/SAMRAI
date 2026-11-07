/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Copy operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_CopyOperation_C
#define included_pdat_CopyOperation_C

#include "SAMRAI/pdat/CopyOperation.h"

namespace SAMRAI {
namespace pdat {

/*
 * Member functions for CopyOperation
 */

template <class TYPE>
SAMRAI_INLINE
SAMRAI_HOST_DEVICE
void
CopyOperation<TYPE>::operator () (
   TYPE& vdst,
   const TYPE& vsrc) const
{
   vdst = vsrc;
}

}
}
#endif
