/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Generator of sequential LocalIds.
 *
 ************************************************************************/

#ifndef included_hier_SequentialLocalIdGenerator
#define included_hier_SequentialLocalIdGenerator

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MessageStream.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class for generating sequential LocalId.
 *
 * Objects of this class generate new LocalId by increasing the last
 * generated value by a constant increment (1 by default).  It is
 * possible to set both the last generated value and the increment.
 */
class SequentialLocalIdGenerator
{

public:
   /*!
    * @brief Default constructor.
    *
    * Construct an object generating sequential values starting with zero.
    */
   SequentialLocalIdGenerator():
      d_last_value(-1),
      d_increment(1) {
   }

   /*!
    * @brief Constructor.
    *
    * Construct an object with user-defined last value and increment.
    */
   SequentialLocalIdGenerator(
      const LocalId& last_value,
      const LocalId& increment = LocalId(1)):
      d_last_value(last_value),
      d_increment(increment) {
   }

   /*!
    * @brief Destructor.
    */
   ~SequentialLocalIdGenerator() = default;

   /*!
    * @brief Return a LocalId that is greater than the previous
    * value by the increment value.
    */
   LocalId nextValue() {
      d_last_value += d_increment;
      return d_last_value;
   }

   /*
    * @brief Reset the last value.
    *
    * @param [in] value
    */
   void setLastValue(const LocalId& last_value) {
      d_last_value = last_value;
   }

   /*
    * @brief Set the increment value
    *
    * @param [in] increment
    */
   void setIncrement(const LocalId& increment) {
      d_increment = increment;
   }

   //! @brief Pack into a MessageStream.
   friend tbox::MessageStream& operator << (
      tbox::MessageStream& msg,
      const SequentialLocalIdGenerator& id_gen)
   {
      msg << id_gen.d_last_value << id_gen.d_increment;
      return msg;
   }

   //! @brief Unpack from a MessageStream.
   friend tbox::MessageStream& operator >> (
      tbox::MessageStream& msg,
      SequentialLocalIdGenerator& id_gen)
   {
      msg >> id_gen.d_last_value >> id_gen.d_increment;
      return msg;
   }

private:
   LocalId d_last_value;
   LocalId d_increment;

};

}
}

#endif  // included_hier_SequentialLocalIdGenerator
