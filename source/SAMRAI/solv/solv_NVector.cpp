/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to C++ vector implementation for Sundials package.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/solv_NVector.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"

using namespace SAMRAI;
using namespace solv;

#define SABSVEC_CAST(v) \
   (static_cast<SundialsAbstractVector *>(v \
                                          -> \
                                          content))

extern "C" {

void N_VPrint_SAMRAI(
   N_Vector v) {
   SABSVEC_CAST(v)->printVector();
}

}

#endif
