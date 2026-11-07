/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box generation routines.
 *
 ************************************************************************/
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * Default constructor and destructor for BoxGeneratorStrategy.
 *
 *************************************************************************
 */

BoxGeneratorStrategy::BoxGeneratorStrategy()
{
}

BoxGeneratorStrategy::~BoxGeneratorStrategy()
{
}

void
BoxGeneratorStrategy::setMinimumCellRequest(
   size_t num_cells)
{
   NULL_USE(num_cells);
}

}
}
