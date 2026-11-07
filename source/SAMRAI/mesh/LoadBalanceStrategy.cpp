/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box load balancing routines.
 *
 ************************************************************************/
#include "SAMRAI/mesh/LoadBalanceStrategy.h"

#include <cstdlib>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

int LoadBalanceStrategy::s_sequence_number = 0;

/*
 *************************************************************************
 *
 * The constructor and destructor for LoadBalanceStrategy do
 * nothing that could be considered even remotely interesting.
 *
 *************************************************************************
 */

LoadBalanceStrategy::LoadBalanceStrategy()
{
}

LoadBalanceStrategy::~LoadBalanceStrategy()
{
}

/*
 *************************************************************************
 * Report the load balance on processor, primarily
 * for debugging and checking load balance quality.
 *************************************************************************
 */
void
LoadBalanceStrategy::markLoadForPostprocessing(
   int rank,
   double load,
   int nbox)
{
   tbox::plog << "Load mark " << s_sequence_number++
              << " proc " << rank
              << " load " << load
              << " nbox " << nbox
              << "\n";
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
