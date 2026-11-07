/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Simple bit vector of a fixed length (128 bits)
 *
 ************************************************************************/
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace hier {

ComponentSelector::ComponentSelector(
   const bool flag):d_max_bit_index(-1)
{
   int pd_entries = tbox::SAMRAIManager::getMaxNumberPatchDataEntries();
   int num_bitset_elements = pd_entries / C_BITSET_SIZE;
   int num_mod = pd_entries % C_BITSET_SIZE;

   if (num_mod != 0) {
      ++num_bitset_elements;
   }

   std::bitset<C_BITSET_SIZE> l_bits;
   d_bit_vector.resize(num_bitset_elements, l_bits);

   if (flag) {
      // use the bitset "set" operation to set each vector element's
      // bitset values to "true".
      for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
         d_bit_vector[vi].set();
      }
      d_max_bit_index =
         (static_cast<int>(d_bit_vector.size()) * C_BITSET_SIZE) - 1;
   }
}

bool
ComponentSelector::any() const noexcept
{
   std::vector<std::bitset<C_BITSET_SIZE> >::const_iterator iter;
   bool set = false;
   for (iter = d_bit_vector.begin(); iter != d_bit_vector.end() && !set;
        ++iter) {
      set = iter->any();
   }
   return set;
}

int
ComponentSelector::_findMaxIndex(
   const std::vector<std::bitset<C_BITSET_SIZE> >& bits) const noexcept
{
   bool bits_set = false;
   int max_index = -1;
   for (size_t i = 0; i < bits.size() && !bits_set; ++i) {
      bits_set |= bits[i].any();
   }

   if (bits_set) {
      int j = C_BITSET_SIZE - 1;
      while (!bits[_index(j)].test(_element(j))) {
         --j;
      }
      max_index = j;
   }
   return max_index;
}

ComponentSelector
ComponentSelector::operator | (
   const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = d_bit_vector[vi] | flags.d_bit_vector[vi];
   }
   tmp.d_max_bit_index =
      tbox::MathUtilities<int>::Max(
         d_max_bit_index, flags.d_max_bit_index);
   return tmp;
}

ComponentSelector
ComponentSelector::operator & (
   const ComponentSelector& flags) const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = d_bit_vector[vi] & flags.d_bit_vector[vi];
   }
   tmp.d_max_bit_index = _findMaxIndex(tmp.d_bit_vector);
   return tmp;
}

ComponentSelector
ComponentSelector::operator ! () const
{
   ComponentSelector tmp;
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      tmp.d_bit_vector[vi] = ~(d_bit_vector[vi]);
   }
   tmp.d_max_bit_index = _findMaxIndex(tmp.d_bit_vector);
   return tmp;
}

ComponentSelector&
ComponentSelector::operator |= (
   const ComponentSelector& flags)
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      d_bit_vector[vi] |= flags.d_bit_vector[vi];
   }
   d_max_bit_index =
      tbox::MathUtilities<int>::Max(
         d_max_bit_index, flags.d_max_bit_index);
   return *this;
}

ComponentSelector&
ComponentSelector::operator &= (
   const ComponentSelector& flags)
{
   for (size_t vi = 0; vi < d_bit_vector.size(); ++vi) {
      d_bit_vector[vi] &= flags.d_bit_vector[vi];
   }
   d_max_bit_index = _findMaxIndex(d_bit_vector);
   return *this;
}

void
ComponentSelector::printClassData(
   std::ostream& os) const
{
   int i;
   const int number_of_bits = getSize();
   for (i = 0; i < number_of_bits; ++i) {
      os << " | Bit " << i << " = " << isSet(i);
   }
   os << "|\n";
}

}
}
