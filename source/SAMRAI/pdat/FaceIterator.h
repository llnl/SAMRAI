/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Iterator for face centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_FaceIterator
#define included_pdat_FaceIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceIndex.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class FaceIterator is an iterator that provides methods for
 * stepping through the index space associated with a face centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 * \verbatim
 * hier::Box box;
 * ...
 * FaceIterator cend(box, axis, false);
 * for (FaceIterator c(box, axis, true); c != cend; ++c) {
 *    // use index c of the box
 * }
 * \endverbatim
 * Note that the face iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see FaceData
 * @see FaceGeometry
 * @see FaceIndex
 */

class FaceIterator
{
public:
   /**
    * Copy constructor for the face iterator
    */
   FaceIterator(
      const FaceIterator& iterator);

   /**
    * Assignment operator for the face iterator.
    */
   FaceIterator&
   operator = (
      const FaceIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the face iterator.
    */
   ~FaceIterator();

   /**
    * Extract the face index corresponding to the iterator position in the box.
    */
   const FaceIndex&
   operator * () const
   {
      return d_index;
   }

   /**
    * Extract a pointer to the face index corresponding to the iterator
    * position in the box.
    */
   const FaceIndex *
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   FaceIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   FaceIterator
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const FaceIterator& iterator) const
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
      const FaceIterator& iterator) const
   {
      TBOX_ASSERT(d_box.isSpatiallyEqual(iterator.d_box));
      TBOX_ASSERT(d_box.isIdEqual(iterator.d_box));
      return d_index != iterator.d_index;
   }

private:
   friend FaceIterator
   FaceGeometry::begin(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);
   friend FaceIterator
   FaceGeometry::end(
      const hier::Box& box,
      tbox::Dimension::dir_t axis);

   /**
    * Constructor for the face iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   FaceIterator(
      const hier::Box& box,
      const tbox::Dimension::dir_t axis,
      bool begin);

   // Unimplemented default constructor.
   FaceIterator();

   FaceIndex d_index;
   hier::Box d_box;
};

}
}

#endif
