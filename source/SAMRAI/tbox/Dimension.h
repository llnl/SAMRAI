/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Dimension class for abstracting dimension
 *
 ************************************************************************/

#ifndef included_tbox_Dimension
#define included_tbox_Dimension

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <limits>

namespace SAMRAI {
namespace tbox {

class DatabaseBox;

/**
 * Class Dimension is used to represent the dimension of a SAMRAI
 * object.
 *
 * The maximum dimension is set at compile time using a flag to the
 * configure script.  This is used to allocate arrays in some lower
 * level classes such as IntVector.  If dynamic memory allocation is
 * used the performance impact is significant; a maximum dimension
 * allows for stack based memory allocation in performance critical
 * classes at the expense of wasting storage for objects with
 * dimension less than the maximum dimension.
 *
 * A class is used rather than a simple short or integer to provide
 * enhanced type safety.
 *
 */

class Dimension
{
public:
   //! @brief Primitive type for direction or dimension.
   using dir_t = unsigned short int; 

   /**
    * Constructor for Dimension, object is built using the specified dimension
    *
    * Note that the constructor is "explicit" thus making automatic
    * type conversions from integers impossible.  This is intentionally to
    * avoid unintended conversions.
    *
    * @pre (dim > 0) && (dim <= SAMRAI::MAX_DIM_VAL)
    */
   constexpr explicit Dimension(
      const unsigned short& dim) noexcept : d_dim(dim)
   {
      TBOX_CONSTEXPR_DIM_ASSERT(dim > 0 && dim <= SAMRAI::MAX_DIM_VAL);
   }

   /**
    * Construct a dimension equal to the argument.
    */
   constexpr Dimension(const Dimension& dimension) noexcept = default;

   /**
    * Equality operator.
    */
   constexpr bool
   operator == (
      const Dimension& rhs) const noexcept
   {
      return d_dim == rhs.d_dim;
   }

   /**
    * Inequality operator.
    */
   constexpr bool
   operator != (
      const Dimension& rhs) const noexcept
   {
      return d_dim != rhs.d_dim;
   }

   /**
    * Greater than operator.
    */
   constexpr bool
   operator > (
      const Dimension& rhs) const noexcept
   {
      return d_dim > rhs.d_dim;
   }

   /**
    * Greater than or equal operator.
    */
   constexpr bool
   operator >= (
      const Dimension& rhs) const noexcept
   {
      return d_dim >= rhs.d_dim;
   }

   /**
    * Less than operator.
    */
   constexpr bool
   operator < (
      const Dimension& rhs) const noexcept
   {
      return d_dim < rhs.d_dim;
   }

   /**
    * Less than or equal operator.
    */
   constexpr bool
   operator <= (
      const Dimension& rhs) const noexcept
   {
      return d_dim <= rhs.d_dim;
   }

   /**
    * Returns the dimension of the Dimension as an unsigned short.
    *
    * The method is provided to allow sizing of arrays based on the
    * dimension and for iteration.  In general this should not be
    * used for comparisons, the Dimension comparison operations are
    * better suited for that purpose.
    */
   constexpr unsigned short
   getValue() const noexcept
   {
      return d_dim;
   }

   /**
    * Returns the maximum dimension for the currently compiled library
    * as a Dimension object.
    *
    * When the SAMRAI library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    */
   static const Dimension&
   getMaxDimension() noexcept
   {
      static Dimension dim(SAMRAI::MAX_DIM_VAL);
      return dim;
   }

   /**
    * Output operator for debugging and error messages.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Dimension& rhs);

private:
   /*
    * Unimplemented default constructor.
    */
   Dimension() = delete;

   /**
    * Assignment operator is private to prevent dimensions
    * from being assigned.  This was done to improve type
    * safety.
    */
   Dimension&
   operator = (
      const Dimension& rhs) = delete;

   dir_t d_dim;
};

}
}

#endif
