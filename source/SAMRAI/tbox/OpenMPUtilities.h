/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for using OpenMP.
 *
 ************************************************************************/

/*
 * @brief Macros defined for using OpenMP, with sensible definitions
 * when not using it.
 *
 * OpenMP symbols beginning with omp_ is prepended with TBOX_ to
 * indicate it is a SAMRAI toolbox macro.  Macros for OpenMP functions
 * must have argument list, even if it is empty.
 */

#ifndef included_tbox_OpenMPUtilities
#define included_tbox_OpenMPUtilities

#ifdef _OPENMP

#include "SAMRAI/SAMRAI_config.h"

#include <omp.h>

#define TBOX_omp_version _OPENMP

#define TBOX_omp_lock_t omp_lock_t

#define TBOX_omp_init_lock(LOCK_PTR) omp_init_lock(LOCK_PTR)
#define TBOX_omp_destroy_lock(LOCK_PTR) omp_destroy_lock(LOCK_PTR)

#define TBOX_omp_set_lock(LOCK_PTR) omp_set_lock(LOCK_PTR)
#define TBOX_omp_unset_lock(LOCK_PTR) omp_unset_lock(LOCK_PTR)

#define TBOX_omp_get_num_threads() omp_get_num_threads()
#define TBOX_omp_get_max_threads() omp_get_max_threads()

#define TBOX_IF_SINGLE_THREAD(CODE) \
   {   \
      if (omp_get_num_threads() == 1) {  \
         CODE      \
      }    \
   }

#define TBOX_IF_MULTI_THREAD(CODE) \
   {    \
      if (omp_get_num_threads() > 1) {   \
         CODE      \
      }    \
   }

#define TBOX_IF_IN_PARALLEL_REGION(CODE) \
   {      \
      if (omp_in_parallel()) {                \
         CODE      \
      }    \
   }

#define TBOX_IF_NOT_IN_PARALLEL_REGION(CODE) \
   {  \
      if (!omp_in_parallel()) {                \
         CODE      \
      }    \
   }

#define TBOX_IF_HAVE_OPENMP(CODE) { CODE }

#define TBOX_IF_NOT_HAVE_OPENMP(CODE)

#else

#define TBOX_omp_version 0

#define TBOX_omp_lock_t int

#define TBOX_omp_init_lock(LOCK_PTR)
#define TBOX_omp_destroy_lock(LOCK_PTR)

#define TBOX_omp_set_lock(LOCK_PTR)
#define TBOX_omp_unset_lock(LOCK_PTR)

#define TBOX_omp_get_num_threads() (1)
#define TBOX_omp_get_max_threads() (1)

#define TBOX_IF_SINGLE_THREAD(CODE) { CODE }

#define TBOX_IF_MULTI_THREAD(CODE)

#define TBOX_IF_IN_PARALLEL_REGION(CODE)

#define TBOX_IF_NOT_IN_PARALLEL_REGION(CODE) { CODE }

#define TBOX_IF_HAVE_OPENMP(CODE)

#define TBOX_IF_NOT_HAVE_OPENMP(CODE) { CODE }

#endif

#endif
