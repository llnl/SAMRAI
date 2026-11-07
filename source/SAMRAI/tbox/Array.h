/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A simple array template class
 *
 ************************************************************************/

#ifndef included_tbox_Array
#define included_tbox_Array

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/ReferenceCounter.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

/**
 * Class Array<TYPE> defines a smart pointer to an array of TYPE and
 * manages all reference counting and deallocation of the array (even
 * if the data was originally allocated from an arena).  When the
 * reference count on a Array<TYPE> object goes to zero, the array
 * objects are automatically deallocated.  The array class frees the
 * user from deleting and tracking aliases for object arrays.
 *
 * A block with references count and arena pointer is allocated for
 * all non-empty arrays.  These reference counted blocks are freed at
 * the end of the lifetime of the array.
 *
 * Class TYPE must define a copy constructor and an assignment
 * operator.
 *
 * @see ReferenceCounter
 */

template<class TYPE>
class Array
{
public:
   /*
    * This is a class used as a flag to ensure a different constructor
    * type signature for the uninitialized Array constructor.  Only
    * the Array::UNINITIALIZED value is every expected to be of this
    * type.
    */
   class DoNotInitialize
   {
public:
      DoNotInitialize() {
      }
   };

   /*
    * The flag value for use in the Array uninitialized constructor.
    */
   static const typename Array<TYPE>::DoNotInitialize UNINITIALIZED;

   /**
    * Create an array of zero elements.
    */
   Array();

   /**
    * Create an array of ``n'' elements.
    *
    * Elements will be initialized with "n" copies of default_value.
    * If not default_value is supplied the default_constructor is
    * invoked to create a default value.
    *
    */
   explicit Array(
      const int n,
      const TYPE& default_value = TYPE());

   /**
    * Create an array of ``n'' uninitialized elements.
    *
    * The Array::UNINITIALIZED value should be used for the
    * second argument to flag that the array is uninitialized.
    *
    * CAUTION: Invoking this constructor will potentially result in a
    * core dump as the element objects will not be initialized (the
    * default constructor is not invoked).  If TYPE is a builtin type
    * the arrays values should be assigned before use.  If TYPE is a
    * class, use the new placement operator on each array element
    * location to invoke a constructor.  This is shown in the
    * following example:
    *
    * \code
    * for(int i = 0; i < d_elements; ++i) {
    *       void *p = &d_objects[i];
    *       (void) new (p) TYPE(arg1, arg2);
    *    }
    *
    * \endcode
    *
    * This constructor may be used to optimize the construction of
    * Arrays when the elements are known to be assigned to some value
    * after construction and thus does not need to be initialized.  A
    * loop over the array elements doing an assiggment to the
    * default_value is avoided.
    *
    * @param n
    * @param do_not_initialize_flag
    */

   Array(
      const int n,
      const typename Array::DoNotInitialize& do_not_initialize_flag);

   /**
    * Copy constructor for the array.  This creates an alias to the
    * right hand side and increments the reference count.
    *
    * CAUTION: invoking resizeArray() forces a deep copy.
    * Upon return, two objects that formerly were aliases to the
    * same underlying data will point to separate data.  For this
    * reason, it is best to pass a Array by reference, instead
    * of by value.
    */
   Array(
      const Array& rhs);

   /**
    * Destructor for the array.  If the reference count for the array data
    * has gone to zero, then the array data is deallocated from the memory
    * arena from which it was allocated.
    */
   ~Array();

   /**
    * Array assignment.  The assignment operator copies a pointer to the
    * array data and increments the reference count.  Both array objects refer
    * to the same data, and changes to individual array entry values in one will
    * be reflected in the other array.  However, this assignment operation DOES NOT
    * involve a "deep copy" (see the resizeArray() routines below). Thus, changes
    * to one Array object container will not necessarily be reflected in the
    * other container.
    */
   Array&
   operator = (
      const Array& rhs);

   /**
    * Non-const array subscripting.  Return a reference the object at array
    * index ``i'' (between 0 and N-1, where N is the number of elements in
    * the array.
    *
    * @param i Array index of item whose reference is to be returned.
    *
    * @pre (i >= 0) && (i < size())
    */
   TYPE&
   operator [] (
      const int i);

   /**
    * Const array subscripting.  Return a const reference to the object
    * at array index ``i'' (between 0 and N-1, where N is the number of
    * elements in the array.
    *
    * @param i Array index of item whose reference is to be returned.
    *
    * @pre (i >= 0) && (i < size())
    */
   const TYPE&
   operator [] (
      const int i) const;

   /**
    * Test whether the array is NULL (has any elements).
    */
   bool
   isNull() const;

   /**
    * Test whether the array is empty (has no elements).
    *
    * Identical to isNull() but this method is common to several
    * container classes, including STL classes.
    */
   bool
   empty() const;

   /**
    * Set the length of the array to zero.  If the reference count for
    * the objects has dropped to zero, then the array data is deallocated.
    */
   void
   setNull();

   /**
    * Set the length of the array to zero.  If the reference count for
    * the objects has dropped to zero, then the array data is deallocated.
    *
    * Identical to setNull() but  this method is common to several
    * container classes, including STL classes.
    */
   void
   clear();

   /**
    * Return a non-const pointer to the i-th object.  The index must be
    * between 0 and N-1, where N is the number of elements in the array.
    *
    * @param i Array index of item whose reference is to be returned.
    *
    * @pre (i >= 0) && (i < size())
    */
   TYPE *
   getPointer(
      const int i = 0);

   /**
    * Return a const pointer to the i-th object.  The index must be
    * between 0 and N-1, where N is the number of elements in the array.
    *
    * @param i Array index of item whose reference is to be returned.
    *
    * @pre (i >= 0) && (i < size())
    */
   const TYPE *
   getPointer(
      const int i = 0) const;

   /**
    * Return the number of elements in the array.
    */
   int
   getSize() const;

   /**
    * Return the number of elements in the array.  Identical to getSize(),
    * but this method is common to several container classes.
    */
   int
   size() const;

   /**
    * Resize the array by allocating new array storage and copying from the
    * old array into the new; i.e., a "deep" copy.  Space for the new array
    * is allocated via the standard ``new'' operator.
    *
    * Elements added be initialized with copies of "default" using the
    * copy constructor for classes and assignment for internal types.
    */
   void
   resizeArray(
      const int n,
      const TYPE& default_value = TYPE());

   /**
    *
    * Adds a new element at the end of the array, after its current
    * last element. The content of this new element is initialized to
    * a copy of value.
    *
    */
   void
   push_back(
      const TYPE& value);

   /**
    * Returns a reference to the last element in the array container.
    *
    * @pre size() > 0
    */
   const TYPE&
   back();

   /**
    *
    * Removes from the array container a single element at position.
    *
    * @param position Array index of element to be removed.
    *
    * @pre (position >= 0) && (position < size())
    */
   void
   erase(
      const int position);

private:
   size_t
   align(
      const size_t bytes);

   static const size_t ALLOCATION_ALIGNMENT = 16;

   void
   deleteObjects();

   TYPE* d_objects;
   ReferenceCounter* d_counter;
   int d_elements;
};

template<>
Array<bool>::Array(
   const int n,
   const bool& default_value);
template<>
Array<char>::Array(
   const int n,
   const char& default_value);
template<>
Array<int>::Array(
   const int n,
   const int& default_value);
template<>
Array<float>::Array(
   const int n,
   const float& default_value);
template<>
Array<double>::Array(
   const int n,
   const double& default_value);

template<>
void
Array<bool>::deleteObjects();
template<>
void
Array<char>::deleteObjects();
template<>
void
Array<int>::deleteObjects();
template<>
void
Array<float>::deleteObjects();
template<>
void
Array<double>::deleteObjects();

}
}

/*
 * Default assume Array is not a standard type
 */

#include "SAMRAI/tbox/Array.cpp"

#endif
