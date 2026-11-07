/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and LICENSE.
*
* Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#ifndef included_pdat_DoubleAttributeId_h
#define included_pdat_DoubleAttributeId_h

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Attribute identifying class for double-valued attributes
 */
class DoubleAttributeId
{
public:
   explicit DoubleAttributeId(
      int value);
   DoubleAttributeId(
      const DoubleAttributeId& other);
   ~DoubleAttributeId();
   DoubleAttributeId&
   operator = (
      const DoubleAttributeId& rhs)
   {
      if (this != &rhs) {
         d_val = rhs.d_val;
      }
      return *this;
   }
   bool
   operator == (
      const DoubleAttributeId& other) const
   {
      return d_val == other.d_val;
   }
   bool
   operator != (
      const DoubleAttributeId& other) const
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
   DoubleAttributeId();
   int d_val;
}; // end class DoubleAttributeId

}
}

#endif
