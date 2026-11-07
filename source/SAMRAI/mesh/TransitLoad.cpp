/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Container of loads for TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_TransitLoad_C
#define included_mesh_TransitLoad_C

#include "SAMRAI/mesh/TransitLoad.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 ***********************************************************************
 ***********************************************************************
 */
TransitLoad::TransitLoad():
   d_allow_box_breaking(true),
   d_threshold_width(1.0e-12)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
TransitLoad::TransitLoad(
   const TransitLoad& other):
   d_allow_box_breaking(other.d_allow_box_breaking),
   d_threshold_width(other.d_threshold_width)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TransitLoad::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   NULL_USE(detail_depth);
   co << border
      << getSumLoad() << " units in " << getNumberOfItems() << " items from "
      << getNumberOfOriginatingProcesses() << " processes.\n";
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

#endif
