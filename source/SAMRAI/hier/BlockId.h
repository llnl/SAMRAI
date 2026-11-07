/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/

#ifndef included_hier_BlockId
#define included_hier_BlockId

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Generic identifier for identifying the block id.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.
 */
class BlockId
{

public:

   //! @brief Primitive type for the number associated with a block.
   using block_t = unsigned int;

   /*!
    * @brief Default constructor sets the value to invalid.
    */
   constexpr BlockId() noexcept :
      d_value(s_invalid_val)
   {
   }

   /*!
    * @brief Copy constructor.
    */
   constexpr BlockId(
      const BlockId& other) noexcept = default;

   /*!
    * @brief Construct from an unsigned int. 
    *
    * This method is explicit to prevent automatic conversion.
    */
   constexpr explicit BlockId(
      const unsigned int& value) noexcept :
      d_value(value)
   {
   }

   /*!
    * @brief Construct from a signed int.
    *
    * This method is explicit to prevent automatic conversion.
    *
    * @pre value >= 0
    */
   constexpr explicit BlockId(
      const int& value) noexcept :
      d_value(static_cast<unsigned int>(value))
   {
      TBOX_CONSTEXPR_ASSERT(value >=0);
   }

   /*!
    * @brief Default destructor.
    */
   ~BlockId() noexcept = default;

   /*!
    * @brief Assignment operator.
    *
    * @param[in] rhs
    *
    * @return @c *this
    */
   constexpr BlockId&
   operator = (
      const BlockId& rhs) noexcept = default;

   /*!
    * @brief Set to an int value.
    *
    * @param[in] rhs
    */
   constexpr void
   setId(
      const int& rhs) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(rhs >= 0); 
      d_value = static_cast<block_t>(rhs);
   }

   constexpr void
   setId(
      const unsigned int& rhs) noexcept
   {
      d_value = rhs;
   }

   /*!
    * @brief Whether the value is valid.
    */
   constexpr bool
   isValid() const noexcept
   {
      return d_value != s_invalid_val;
   }

   /*!
    * @brief Access the numerical value.
    */
   constexpr const block_t&
   getBlockValue() const noexcept
   {
      return d_value;
   }

   /*!
    * @brief Get the BlockId with a numerical value of zero.
    */
   static constexpr const BlockId&
   zero() noexcept
   {
      return s_zero_id;
   }

   /*!
    * @brief Get the designated invalid value for this class.
    */
   static constexpr const BlockId&
   invalidId() noexcept
   {
      return s_invalid_id;
   }

   //@{

   //! @name Comparison with another BlockId.

   /*!
    * @brief Equality operator.
    *
    * All comparison operators compare the numerical value.
    *
    * @param[in] rhs
    */
   constexpr bool
   operator == (
      const BlockId& rhs) const noexcept
   {
      return d_value == rhs.d_value;
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator != (
      const BlockId& rhs) const noexcept
   {
      return d_value != rhs.d_value;
   }

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator < (
      const BlockId& rhs) const noexcept
   {
      return d_value < rhs.d_value;
   }

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator > (
      const BlockId& rhs) const noexcept
   {
      return d_value > rhs.d_value;
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator <= (
      const BlockId& rhs) const noexcept
   {
      return d_value <= rhs.d_value;
   }

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator >= (
      const BlockId& rhs) const noexcept
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
      const block_t& rhs) const noexcept
   {
      return d_value == rhs;
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator != (
      const block_t& rhs) const noexcept
   {
      return d_value != rhs;
   }

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator < (
      const block_t& rhs) const noexcept
   {
      return d_value < rhs;
   }

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator > (
      const block_t& rhs) const noexcept
   {
      return d_value > rhs;
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator <= (
      const block_t& rhs) const noexcept
   {
      return d_value <= rhs;
   }

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const BlockId&);
    *
    * @param[in] rhs
    */
   constexpr bool
   operator >= (
      const block_t& rhs) const noexcept
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
      const BlockId& r)
   {
      co << r.d_value;
      return co;
   }

private:
   /*!
    * @brief Numerical value of the identifier.
    */
   block_t d_value;

   static inline constexpr unsigned int s_zero_val = 0;
   static inline constexpr unsigned int s_invalid_val =
      tbox::MathUtilities<int>::getMax();

   /*!
    * @brief BlockId with a numerical value of zero.
    */
   static const BlockId s_zero_id;

   /*!
    * @brief Definition of invalid BlockId.
    */
   static const BlockId s_invalid_id;

};

}
}

#endif  // included_hier_BlockId
