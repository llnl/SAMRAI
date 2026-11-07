/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Copy operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_CopyOperation
#define included_pdat_CopyOperation

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class CopyOperation<TYPE> encapsulates a copy
 * operation into an object.
 */

template<class TYPE>
struct CopyOperation
{
   /*!
    * The operator copies the source value to the destination.
    */
   SAMRAI_HOST_DEVICE
   void
   operator () (
      TYPE& vdst,
      const TYPE& vsrc) const;
};

}
}

#include "SAMRAI/pdat/CopyOperation.cpp"

#endif
