/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple integer id and namestring variable context
 *
 ************************************************************************/
#include "SAMRAI/hier/VariableContext.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

int VariableContext::s_instance_counter = 0;

/*
 *************************************************************************
 *
 * The constructor copies the name of the variable context, obtains
 * a unique instance number, and increments the number of global
 * instances.  The destructor releases the name storage but does not
 * decrease the instance count; instance numbers are not recycled.
 *
 *************************************************************************
 */

VariableContext::VariableContext(
   const std::string& name)
{
   TBOX_ASSERT(!name.empty());

   d_index = s_instance_counter++;
   d_name = name;
}

}
}
