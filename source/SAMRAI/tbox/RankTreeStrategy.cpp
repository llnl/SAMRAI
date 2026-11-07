/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#include "SAMRAI/tbox/RankTreeStrategy.h"

#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int RankTreeStrategy::s_invalid_rank(-1);
const unsigned int
RankTreeStrategy::s_invalid_child_number(
   MathUtilities<unsigned int>::getMax() / 2);

RankTreeStrategy::RankTreeStrategy()
{
}

RankTreeStrategy::~RankTreeStrategy()
{
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
