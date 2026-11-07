/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple bit vector.
 *
 ************************************************************************/

#ifndef included_hier_ComponentSelector
#define included_hier_ComponentSelector

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <bitset>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Bit vector used to indicate on which patch data elements to
 * apply operations.
 *
 * Class ComponentSelector implements a bit vector of a fixed
 * length and is typically used to apply operations on subsets of entries
 * in the patch data array owned by a patch (e.g., allocate/deallocate).
 * All ComponentSelector objects have the same bit vector length that is
 * established by the SAMRAIManager utility. See the documentation
 * of the SAMRAIManager utility for information about changing this
 * maximum value.
 *
 * @see tbox::SAMRAIManager
 */
class ComponentSelector
{
public:
   /*!
    * @brief Create a ComponentSelector with all values initialized to
    * the specified flag value.
    *
    * @param[in]  flag @b Default: false
    */
   explicit ComponentSelector(
      const bool flag = false);

   /*!
    * @brief Copy constructor that create a component selector identical
    * to the argument.
    */
   ComponentSelector(
      const ComponentSelector& flags) = default;

   ComponentSelector(
      ComponentSelector&& flags) = default;

   /*!
    * @brief The destructor for a component selector does nothing interesting.
    */
   ~ComponentSelector() = default;

   /*!
    * @brief Get the size of the ComponentSelector.
    *
    * @return total number of flags (i.e., bits) in this component selector.
    */
   int
   getSize() const noexcept
   {
      return static_cast<int>(d_bit_vector.size()) * C_BITSET_SIZE;
   }

   /*!
    * @brief Copy assignment operator.
    */
   ComponentSelector&
   operator = (
      const ComponentSelector& flags) = default;

   ComponentSelector&
   operator = (
      ComponentSelector&& flags) = default;

   /*!
    * @brief Equality operator.  Two ComponentSelector objects are
    * equal when all their bits are the same.
    */
   bool
   operator == (
      const ComponentSelector& flags) const noexcept
   {
      return d_bit_vector == flags.d_bit_vector;
   }

   /*!
    * @brief Inequality operator.  Two ComponentSelector objects are
    * unequal if any of their bits are different.
    */
   bool
   operator != (
      const ComponentSelector& flags) const noexcept
   {
      return d_bit_vector != flags.d_bit_vector;
   }

   /*!
    * @brief Generate and return a component selector set to the bitwise
    * logical OR of this component selector and the argument component
    * selector.
    *
    * @param[in]  flags The ComponentSelector used to apply the logical OR
    *             operation
    *
    * @return The ComponentSelector which is the logical OR of this component.
    */
   ComponentSelector
   operator | (
      const ComponentSelector& flags) const;

   /*!
    * @brief Generate a component selector set to the bitwise logical AND
    * of this component selector and the argument component selector.
    *
    * @param[in]  flags The ComponentSelector used to apply the logical
    *             AND operation
    *
    * @return The ComponentSelector which is the logical AND of this component.
    */
   ComponentSelector
   operator& (
      const ComponentSelector& flags) const;

   /*!
    * @brief Generate a component selector set to the bitwise logical
    * negation of this component selector.
    *
    * @return The ComponentSelector which is the logical negation of this
    *         component.
    */
   ComponentSelector
   operator ! () const;

   /*!
    * @brief Set all bits in this component selector to the bitwise
    * logical OR of this component selector and the argument component
    * selector.
    *
    * This is a modifying operation on this ComponentSelector.
    *
    * @return This selector, modified.
    */
   ComponentSelector&
   operator |= (
      const ComponentSelector& flags);

   /*!
    * @brief Set all bits in this component selector to the bitwise
    * logical AND of this component selector and the argument component
    * selector.
    *
    * This is a modifying operation on this ComponentSelector.
    *
    * @return This selector, modified.
    */
   ComponentSelector&
   operator &= (
      const ComponentSelector& flags);

   /*!
    * @brief Check whether the specified bit vector position is true.
    *
    * @param[in]  i The position in the bit vector to check.
    *
    * @return True if the bit at position i is set to true.
    *
    * @pre (i >= 0)
    */
   bool
   isSet(
      const int i) const noexcept
   {
      TBOX_ASSERT(i >= 0);
      return i < getSize() && d_bit_vector[_index(i)].test(_element(i));
   }

   /*!
    * @brief Set the specified bit vector position to true.
    *
    * @param[in]  i The position in the bit vector to set to true.
    *
    * @pre (i >= 0)
    */
   void
   setFlag(
      const int i)
   {
      TBOX_ASSERT(i >= 0);
      if (i >= getSize()) {
         d_bit_vector.resize(d_bit_vector.size() + 1, d_bit_vector[0]);
         d_bit_vector[d_bit_vector.size()-1].reset();
      }
      d_bit_vector[_index(i)].set(_element(i));
   }

   /*!
    * @brief Set the specified bit vector position to false.
    *
    * @param[in]  i The position in the bit vector to set to false.
    *
    * @pre (i >= 0)
    */
   void
   clrFlag(
      const int i)
   {
      TBOX_ASSERT(i >= 0);
      if (i < getSize()) {
         d_bit_vector[_index(i)].reset(_element(i));
      } 
      d_max_bit_index = _findMaxIndex(d_bit_vector);
   }

   /*!
    * @brief Set all bit vector positions to true.
    */
   void
   setAllFlags()
   {
      for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
         (d_bit_vector[vi]).set();
      }
      d_max_bit_index =
         C_BITSET_SIZE * static_cast<int>(d_bit_vector.size()) - 1;
   }

   /*!
    * @brief Set all bit vector positions to false.
    */
   void
   clrAllFlags()
   {
      for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
         (d_bit_vector[vi]).reset();
      }
      d_max_bit_index = -1;
   }

   /*!
    * @brief Get the index of the highest set position.
    */
   int
   getMaxIndex() const noexcept
   {
      return d_max_bit_index;
   }

   /*!
    * @brief check if any bits in the vector are set to true.
    *
    * @return True if any bit in the vector is set true, otherwise false.
    */
   bool
   any() const noexcept;

   /*!
    * @brief check if no bits in the vector are set to true.
    *
    * This is essentially determining if all bits are false.
    *
    * @return True if no bits in the vector are set to true, otherwise false.
    */
   bool
   none() const noexcept
   {
      return !any();
   }

   /*!
    * @brief Print the bit vector data to the specified output stream.
    */
   void
   printClassData(
      std::ostream& os = tbox::plog) const;

private:
   /*
    *  Default length of std::bitset entries used in bit vector representation.
    */
   static constexpr int C_BITSET_SIZE = 1024;

   std::vector<std::bitset<C_BITSET_SIZE> > d_bit_vector;

   // the index of the highest bit set in the bit vector.
   // used for efficiency
   int d_max_bit_index;

   int
   _findMaxIndex(
      const std::vector<std::bitset<C_BITSET_SIZE> >& bits) const noexcept;
   // private function to return the index into the d_bit_vector
   constexpr int
   _index(
      const int i) const noexcept
   {
      return i / C_BITSET_SIZE;
   }

   // private function to return the element within the d_bit_vector[i]
   // bitset.
   constexpr int
   _element(
      const int i) const noexcept
   {
      return i % C_BITSET_SIZE;
   }
};

}
}

#endif
