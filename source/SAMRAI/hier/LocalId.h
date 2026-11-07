/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Generic identifier used on a single process.
 *
 ************************************************************************/

#ifndef included_hier_LocalId
#define included_hier_LocalId

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Generic identifier for identifying things on the local
 * process.
 *
 * The LocalId can be combined with a process rank, as is done in
 * GlobalId, to create global identifiers.
 *
 * Comparison operators are provided to define the sorted ordering of
 * objects.
 */
class LocalId
{

public:
   /*!
    * @brief Default constructor.
    */
   constexpr LocalId() noexcept : d_value(s_invalid_val)
   {
   }

   /*!
    * @brief Copy constructor.
    */
   constexpr LocalId(const LocalId& other) noexcept = default;

   /*!
    * @brief Construct from a numerical value.
    *
    * This method is explicit to prevent automatic conversion.
    */
   constexpr explicit LocalId(
      int value) noexcept : d_value(value)
   {
   }

   /*!
    * @brief Default destructor.
    */
   ~LocalId() noexcept = default;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   constexpr LocalId&
   operator = (
      const LocalId& rhs) noexcept = default;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   constexpr LocalId&
   operator = (
      int rhs) noexcept
   {
      d_value = rhs;
      return *this;
   }

   /*!
    * @brief Access the numerical value.
    */
   int&
   getValue() noexcept
   {
      return d_value;
   }

   /*!
    * @brief Access the numerical value.
    */
   constexpr const int&
   getValue() const noexcept
   {
    return d_value;
   }

   /*!
    * @brief Whether value is a valid one (not equal to getInvalidId()).
    */
   constexpr bool
   isValid() const noexcept
   {
      return d_value != s_invalid_val;
   }

   /*!
    * @brief Get the LocalId with a numerical value of zero.
    */
   static const LocalId&
   getZero() noexcept
   {
      return s_zero_id;
   }

   /*!
    * @brief Get the designated invalid value for this class.
    */
   static const LocalId&
   getInvalidId() noexcept
   {
      return s_invalid_id;
   }

   //@{

   //! @name Numerical operations.

   /*!
    * @brief Pre-increment iterator.
    *
    * Pre-increment increments the value and returns the incremented
    * state.
    */
   constexpr LocalId
   operator ++ () noexcept
   {
      ++d_value;
      return *this;
   }

   /*!
    * @brief Post-increment iterator.
    *
    * Post-increment saves the value, increment it and returns an
    * object with the saved value.
    */
   constexpr LocalId
   operator ++ (
      int) noexcept
   {
      int saved = d_value;
      ++d_value;
      return LocalId(saved);
   }

   /*!
    * @brief Addition.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator + (
      const LocalId& rhs) const noexcept
   {
      return LocalId(d_value + rhs.d_value);
   }

   /*!
    * @brief Subtraction.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator - (
      const LocalId& rhs) const noexcept
   {
      return LocalId(d_value - rhs.d_value);
   }

   /*!
    * @brief Multiplication.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator * (
      const LocalId& rhs) const noexcept
   {
      return LocalId(d_value * rhs.d_value);
   }

   /*!
    * @brief Division.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator / (
      const LocalId& rhs) const noexcept
   {
      return LocalId(d_value / rhs.d_value);
   }

   /*!
    * @brief Modulus.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator % (
      const LocalId& rhs) const noexcept
   {
      return LocalId(d_value % rhs.d_value);
   }

   /*!
    * @brief Addition and assignment.
    *
    * @param[in] rhs
    */
   constexpr LocalId&
   operator += (
      const LocalId& rhs) noexcept
   {
      d_value += rhs.d_value;
      return *this;
   }

   /*!
    * @brief Subtraction and assignment.
    *
    * @param[in] rhs
    */
   constexpr LocalId&
   operator -= (
      const LocalId& rhs) noexcept
   {
      d_value -= rhs.d_value;
      return *this;
   }

   /*!
    * @brief Integer addition.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator + (
      int rhs) const noexcept
   {
      return LocalId(d_value + rhs);
   }

   /*!
    * @brief Integer subtraction.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator - (
      int rhs) const noexcept
   {
      return LocalId(d_value - rhs);
   }

   /*!
    * @brief Integer multiplication.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator * (
      int rhs) const noexcept
   {
      return LocalId(d_value * rhs);
   }

   /*!
    * @brief Integer division.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator / (
      int rhs) const noexcept
   {
      return LocalId(d_value / rhs);
   }

   /*!
    * @brief Integer modulus.
    *
    * @param[in] rhs
    */
   constexpr LocalId
   operator % (
      int rhs) const noexcept
   {
      return LocalId(d_value % rhs);
   }

   /*!
    * @brief Integer addition and assignment.
    *
    * @param[in] rhs
    */
   constexpr LocalId&
   operator += (
      int rhs) noexcept
   {
      d_value += rhs;
      return *this;
   }

   /*!
    * @brief Integer subtraction and assignment.
    *
    * @param[in] rhs
    */
   constexpr LocalId&
   operator -= (
      int rhs) noexcept
   {
      d_value -= rhs;
      return *this;
   }

   //@}

   //@{

   //! @name Comparison with another LocalId.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   constexpr bool
   operator == (
      const LocalId& rhs) const noexcept
   {
      return d_value == rhs.d_value;
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator != (
      const LocalId& rhs) const noexcept
   {
      return d_value != rhs.d_value;
   }

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator < (
      const LocalId& rhs) const noexcept
   {
      return d_value < rhs.d_value;
   }

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator > (
      const LocalId& rhs) const noexcept
   {
      return d_value > rhs.d_value;
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator <= (
      const LocalId& rhs) const noexcept
   {
      return d_value <= rhs.d_value;
   }

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator >= (
      const LocalId& rhs) const noexcept
   {
      return d_value >= rhs.d_value;
   }

   //@}

   //@{

   //! @name Comparison with an integer.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   constexpr bool
   operator == (
      int rhs) const noexcept
   {
      return d_value == rhs;
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator != (
      int rhs) const noexcept
   {
      return d_value != rhs;
   }

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator < (
      int rhs) const noexcept
   {
      return d_value < rhs;
   }

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator > (
      int rhs) const noexcept
   {
      return d_value > rhs;
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator <= (
      int rhs) const noexcept
   {
      return d_value <= rhs;
   }

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const LocalId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator >= (
      int rhs) const noexcept
   {
      return d_value >= rhs;
   }

   //@}

   /*!
    * @brief Format and insert object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const LocalId& r);

private:
   /*!
    * @brief Numerical value of the identifier.
    */
   int d_value;

   inline static constexpr int s_zero_val = 0;
   inline static constexpr int s_invalid_val = tbox::MathUtilities<int>::getMax();

   /*!
    * @brief LocalId with a numerical value of zero.
    */
   static const LocalId s_zero_id;

   /*!
    * @brief Definition of invalid LocalId.
    */
   static const LocalId s_invalid_id;

};

}
}

#endif  // included_hier_LocalId
