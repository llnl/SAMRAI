/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Sum operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_SumOperation
#define included_pdat_SumOperation

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class SumOperation<TYPE> encapsulates a summation
 * operation into an object.
 */

template<class TYPE>
struct SumOperation
{
   /*!
    * The operator adds the source value to the destination.
    */
   SAMRAI_HOST_DEVICE
   void
   operator () (
      TYPE& vdst,
      const TYPE& vsrc) const;
};

}
}

#include "SAMRAI/pdat/SumOperation.cpp"

#endif
