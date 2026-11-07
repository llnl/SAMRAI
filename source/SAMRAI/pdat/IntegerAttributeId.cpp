/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and LICENSE.
*
* Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#include "SAMRAI/pdat/IntegerAttributeId.h"

namespace SAMRAI {
namespace pdat {

/**********************************************************************
 * IntegerAttributesId ctor
 *********************************************************************/
IntegerAttributeId::IntegerAttributeId(
   int value):
   d_val(value)
{
}

/**********************************************************************
 * IntegerAttributesId copy ctor
 *********************************************************************/
IntegerAttributeId::IntegerAttributeId(
   const IntegerAttributeId& other):
   d_val(other.d_val)
{
}

/**********************************************************************
 * IntegerAttributesId c'tor (private)
 *********************************************************************/
IntegerAttributeId::IntegerAttributeId():
   d_val(-1)
{
}

/**********************************************************************
 * d'tor
 *********************************************************************/
IntegerAttributeId::~IntegerAttributeId()
{
}

}
}
