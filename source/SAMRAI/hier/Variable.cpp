/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for application-level variables
 *
 ************************************************************************/
#include "SAMRAI/hier/Variable.h"

namespace SAMRAI {
namespace hier {

int Variable::s_instance_counter = 0;

/*
 *************************************************************************
 *
 * The constructor copies the name of the variable, obtains a unique
 * instance number, and increments the number of global instances.
 * The destructor releases the name storage but does not decrease the
 * instance count, since instance numbers are never recycled.
 *
 *************************************************************************
 */

Variable::Variable(
   const std::string& name,
   const std::shared_ptr<PatchDataFactory>& factory):
   d_dim(factory->getDim()),
   d_name(name),
   d_factory(factory)
{
   d_instance = s_instance_counter++;
}

}
}
