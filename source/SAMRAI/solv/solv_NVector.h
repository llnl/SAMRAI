/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   C interface to C++ vector implementation for Sundials package.
 *
 ************************************************************************/
#ifndef included_NVector_SAMRAI
#define included_NVector_SAMRAI

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_SUNDIALS

#include "sundials/sundials_nvector.h"

extern "C" {

/**
 * @brief Helper funtion for printing SAMRAI N_Vector.
 *
 */
void
N_VPrint_SAMRAI(
   N_Vector v);

}

#endif
#endif
