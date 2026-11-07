/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and LICENSE.
*
* Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#ifndef included_pdat_IntegerAttributeId_h
#define included_pdat_IntegerAttributeId_h

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Attribute identifying class for integer-valued attributes
 */
class IntegerAttributeId
{
public:
   explicit IntegerAttributeId(
      int value);
   IntegerAttributeId(
      const IntegerAttributeId& other);
   ~IntegerAttributeId();
   IntegerAttributeId&
   operator = (
      const IntegerAttributeId& rhs)
   {
      if (this != &rhs) {
         d_val = rhs.d_val;
      }
      return *this;
   }
   bool
   operator == (
      const IntegerAttributeId& other) const
   {
      return d_val == other.d_val;
   }
   bool
   operator != (
      const IntegerAttributeId& other) const
   {
      return !this->operator == (other);
   }
   int
   operator () () const
   {
      return d_val;
   }

   friend class Attributes;
private:
   IntegerAttributeId();

   int d_val;
}; // end class IntegerAttributeId.

}
}

#endif
