/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A simple array template class
 *
 ************************************************************************/

#ifndef included_tbox_Array_C
#define included_tbox_Array_C

#include "SAMRAI/tbox/Array.h"

#include <new>
#include <cstdlib>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

template<class TYPE>
const typename Array<TYPE>::DoNotInitialize Array<TYPE>::UNINITIALIZED;

/*
 * Note that this class is specialized for the built-in types to avoid
 * invoking the default ctor for TYPE.   A simple assignment is
 * used for the built-in types.
 */

template<class TYPE>
Array<TYPE>::Array():
   d_objects(0),
   d_counter(0),
   d_elements(0)
{
}

template<class TYPE>
Array<TYPE>::Array(
   const Array<TYPE>& rhs):
   d_objects(rhs.d_objects),
   d_counter(rhs.d_counter),
   d_elements(rhs.d_elements)
{
   if (d_counter) {
      d_counter->addReference();
   }
}

template<class TYPE>
Array<TYPE>::Array(
   const int n,
   const TYPE& default_value)
{
   if (n > 0) {

      d_objects = reinterpret_cast<TYPE *>(malloc(sizeof(TYPE) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         void* p = &d_objects[i];
         (void)new (p)TYPE(default_value);
      }
   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<class TYPE>
Array<TYPE>::Array(
   const int n,
   const typename Array::DoNotInitialize& do_not_initialize_flag)
{
   NULL_USE(do_not_initialize_flag);

   if (n > 0) {
      d_objects = reinterpret_cast<TYPE *>(malloc(sizeof(TYPE) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;
   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<class TYPE>
Array<TYPE>::~Array()
{
   if (d_counter && d_counter->deleteReference()) {
      deleteObjects();
   }
}

template<class TYPE>
Array<TYPE>&
Array<TYPE>::operator = (
   const Array<TYPE>& rhs)
{
   if (this != &rhs) {
      if (d_counter && d_counter->deleteReference()) {
         deleteObjects();
      }
      d_objects = rhs.d_objects;
      d_counter = rhs.d_counter;
      d_elements = rhs.d_elements;
      if (d_counter) {
         d_counter->addReference();
      }
   }
   return *this;
}

template<class TYPE>
void
Array<TYPE>::resizeArray(
   const int n,
   const TYPE& default_value)
{
   if (n != d_elements) {
      Array<TYPE> array(n, default_value);
      const int s = (d_elements < n ? d_elements : n);
      for (int i = 0; i < s; ++i) {
         array.d_objects[i] = d_objects[i];
      }

      this->
      operator = (
         array);
   }
}

template<class TYPE>
void
Array<TYPE>::erase(
   const int position)
{
   TBOX_ASSERT((position >= 0) && (position < size()));

   if (d_elements > 1) {

      int new_d_elements(d_elements - 1);

      TYPE* new_d_objects = reinterpret_cast<TYPE *>(
            malloc(sizeof(TYPE) * new_d_elements));

      /* copy lower part of array */
      for (int j = 0; j < position; ++j) {
         void* p = &new_d_objects[j];
         (void)new (p)TYPE(d_objects[j]);
      }

      /* copy upper part of array */
      for (int j = position + 1; j < d_elements; ++j) {
         void* p = &new_d_objects[j - 1];
         (void)new (p)TYPE(d_objects[j]);
      }

      if (d_counter && d_counter->deleteReference()) {
         deleteObjects();
      }

      d_objects = new_d_objects;
      d_counter = new ReferenceCounter;
      d_elements = new_d_elements;

   } else {
      if (d_counter && d_counter->deleteReference()) {
         deleteObjects();
      }
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }

}

template<class TYPE>
void
Array<TYPE>::deleteObjects()
{
   if (d_objects) {
      for (int i = 0; i < d_elements; ++i) {
         d_objects[i].~TYPE();
      }
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<class TYPE>
TYPE&
Array<TYPE>::operator [] (
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < size()));

   return d_objects[i];
}

template<class TYPE>
const TYPE&
Array<TYPE>::operator [] (
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < size()));

   return d_objects[i];
}

template<class TYPE>
void
Array<TYPE>::setNull()
{
   if (d_counter && d_counter->deleteReference()) {
      deleteObjects();
   }
   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<class TYPE>
void
Array<TYPE>::clear()
{
   if (d_counter && d_counter->deleteReference()) {
      deleteObjects();
   }
   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<class TYPE>
bool
Array<TYPE>::isNull() const
{
   return !d_objects;
}

template<class TYPE>
bool
Array<TYPE>::empty() const
{
   return !d_objects;
}

template<class TYPE>
TYPE *
Array<TYPE>::getPointer(
   const int i)
{
   TBOX_ASSERT((i >= 0) && (i < size()));

   return &d_objects[i];
}

template<class TYPE>
const TYPE *
Array<TYPE>::getPointer(
   const int i) const
{
   TBOX_ASSERT((i >= 0) && (i < size()));

   return &d_objects[i];
}

template<class TYPE>
int
Array<TYPE>::getSize() const
{
   return d_elements;
}

template<class TYPE>
int
Array<TYPE>::size() const
{
   return d_elements;
}

template<class TYPE>
size_t
Array<TYPE>::align(
   const size_t bytes)
{
   size_t aligned = bytes + ALLOCATION_ALIGNMENT - 1;
   aligned -= aligned % ALLOCATION_ALIGNMENT;
   return aligned;
}

template<class TYPE>
void
Array<TYPE>::push_back(
   const TYPE& value)
{
   int i = d_elements;
   resizeArray(i + 1);
   d_objects[i] = value;
}

template<class TYPE>
const TYPE&
Array<TYPE>::back()
{
   TBOX_ASSERT(size() > 0);

   return d_objects[d_elements - 1];
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
