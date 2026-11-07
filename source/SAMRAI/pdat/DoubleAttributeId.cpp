/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and LICENSE.
*
* Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#include "SAMRAI/pdat/DoubleAttributeId.h"

namespace SAMRAI {
namespace pdat {

/**********************************************************************
 * c'tor
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId(
   int value):
   d_val(value)
{
}

/**********************************************************************
 * DoubleAttributeId c'tor (private)
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId():
   d_val(-1)
{
}

/**********************************************************************
 * copy c'tor
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId(
   const DoubleAttributeId& other):
   d_val(other.d_val)
{
}

/**********************************************************************
 * d'tor
 *********************************************************************/
DoubleAttributeId::~DoubleAttributeId()
{
}

}
}
