/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Sum operation on single array data elements templated on data type
 *
 ************************************************************************/

/*
 * Need this so includes script will think this is a templated class:
 * template <class TYPE>
 */

#ifndef included_pdat_SumOperation_C
#define included_pdat_SumOperation_C

#include "SAMRAI/pdat/SumOperation.h"

namespace SAMRAI {
namespace pdat {

/*
 * Member functions for SumOperation
 */

template<class TYPE>
SAMRAI_INLINE
SAMRAI_HOST_DEVICE
void
SumOperation<TYPE>::operator () (
   TYPE& vdst,
   const TYPE& vsrc) const
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   vdst += vsrc;
}

}
}
#endif
