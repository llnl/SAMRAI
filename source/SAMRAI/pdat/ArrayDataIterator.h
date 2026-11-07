/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Iterator for array patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataIterator
#define included_pdat_ArrayDataIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class ArrayDataIterator is an iterator that provides methods for
 * stepping through the index space associated with a pdat_ArrayData object.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 *
 * \verbatim
 * hier::Box box;
 * ...
 * ArrayDataIterator cend(box, false);
 * for (ArrayDataIterator c(box, true); c != cend; c++) {
 *    // use index c of the box
 * }
 * \endverbatim
 *
 * Note that the array data iterator may not compile to efficient code,
 * depending on your compiler.  Many compilers are not smart enough to
 * optimize the looping constructs and indexing operations.
 *
 * @see ArrayData
 * @see hier::Index
 */

class ArrayDataIterator
{
public:
   /**
    * Constructor for the array data iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   ArrayDataIterator(
      const hier::Box& box,
      bool begin);

   /**
    * Copy constructor for the array data iterator
    */
   ArrayDataIterator(
      const ArrayDataIterator& iterator);

   /**
    * Assignment operator for the array data iterator.
    */
   ArrayDataIterator&
   operator = (
      const ArrayDataIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the array data iterator.
    */
   ~ArrayDataIterator();

   /**
    * Extract the index corresponding to the iterator position in the box.
    */
   const hier::Index&
   operator * () const
   {
      return d_index;
   }

   /**
    * Extract pointer to the index corresponding to the iterator position in
    * the box.
    */
   const hier::Index *
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   ArrayDataIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   ArrayDataIterator
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const ArrayDataIterator& iterator) const
   {
      return d_index == iterator.d_index;
   }

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const ArrayDataIterator& iterator) const
   {
      return d_index != iterator.d_index;
   }

private:
   hier::Index d_index;
   hier::Box d_box;
};

}
}

#endif
