/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Accesses system times.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Clock.h"

#include <cstdlib>

namespace SAMRAI {
namespace tbox {

#ifdef SAMRAI_HAVE_SYS_TIMES_H
struct tms Clock::s_tms_buffer;
#endif
clock_t Clock::s_null_clock_t;

}
}
