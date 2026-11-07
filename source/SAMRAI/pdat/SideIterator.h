/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Iterator for side centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_SideIterator
#define included_pdat_SideIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class SideIterator is an iterator that provides methods for
 * stepping through the index space associated with a side centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 * \verbatim
 * hier::Box box;
 * ...
 * SideIterator cend(box, axis, false);
 * for (SideIterator c(box, axis, true); c != cend; ++c) {
 *    // use index c of the box
 * }
 * \endverbatim
 * Note that the side iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see SideData
 * @see SideGeometry
 * @see SideIndex
 */

class SideIterator
{
public:
   /**
    * Copy constructor for the side iterator
    */
   SideIterator(
      const SideIterator& iterator);

   /**
    * Assignment operator for the side iterator.
    */
   SideIterator&
   operator = (
      const SideIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the side iterator.
    */
   ~SideIterator();

   /**
    * Extract the side index corresponding to the iterator position in the box.
    */
   const SideIndex&
   operator * () const
   {
      return d_index;
   }

   /**
    * Extract a pointer to the side index corresponding to the iterator
    * position in the box.
    */
   const SideIndex *
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   SideIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   SideIterator
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const SideIterator& iterator) const
   {
      TBOX_ASSERT(d_box.isSpatiallyEqual(iterator.d_box));
      TBOX_ASSERT(d_box.isIdEqual(iterator.d_box));
      return d_index == iterator.d_index;
   }

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const SideIterator& iterator) const
   {
      TBOX_ASSERT(d_box.isSpatiallyEqual(iterator.d_box));
      TBOX_ASSERT(d_box.isIdEqual(iterator.d_box));
      return d_index != iterator.d_index;
   }

private:
   friend SideIterator
   SideGeometry::begin(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);
   friend SideIterator
   SideGeometry::end(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   /**
    * Constructor for the side iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   SideIterator(
      const hier::Box& box,
      const tbox::Dimension::dir_t axis,
      bool begin);

   // Unimplemented default constructor.
   SideIterator();

   SideIndex d_index;
   hier::Box d_box;
};

}
}

#endif
