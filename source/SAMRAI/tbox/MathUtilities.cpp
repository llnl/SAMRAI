/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utilities class to access common POSIX constants and math ops
 *
 ************************************************************************/

#ifndef included_tbox_MathUtilities_C
#define included_tbox_MathUtilities_C

#include <limits>

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Routines to initialize vectors and arrays to signaling NaNs.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setVectorToSignalingNaN(
   std::vector<TYPE>& vector)
{
   for (int i = 0; i < static_cast<int>(vector.size()); ++i) {
      vector[i] = getSignalingNaN();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToSignalingNaN(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; ++i) {
      array[i] = getSignalingNaN();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize vectors and arrays to max value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setVectorToMax(
   std::vector<TYPE>& vector)
{
   for (int i = 0; i < static_cast<int>(vector.size()); ++i) {
      vector[i] = getMax();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMax(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; ++i) {
      array[i] = getMax();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize vectors and arrays to min value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setVectorToMin(
   std::vector<TYPE>& vector)
{
   for (int i = 0; i < static_cast<int>(vector.size()); ++i) {
      vector[i] = getMin();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMin(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; ++i) {
      array[i] = getMin();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize vectors and arrays to epsilon value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setVectorToEpsilon(
   std::vector<TYPE>& vector)
{
   for (int i = 0; i < static_cast<int>(vector.size()); ++i) {
      vector[i] = getEpsilon();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToEpsilon(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; ++i) {
      array[i] = getEpsilon();
   }
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getZero()
{
   return static_cast< TYPE >( 0 );
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getOne()
{
   return static_cast< TYPE >( 1 );
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getSignalingNaN()
{
   return std::numeric_limits<TYPE>::signaling_NaN();
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getMin()
{
   return std::numeric_limits<TYPE>::min();
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getEpsilon()
{
   return std::numeric_limits<TYPE>::epsilon();
}

template<class TYPE>
bool
MathUtilities<TYPE>::isNaN(
   const TYPE& value)
{
   NULL_USE(value);
   return false;
}

template<class TYPE>
bool
MathUtilities<TYPE>::equalEps(
   const TYPE& a,
   const TYPE& b)
{
   return a == b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Min(
   TYPE a,
   TYPE b)
{
   return a < b ? a : b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Max(
   TYPE a,
   TYPE b)
{
   return a > b ? a : b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Abs(
   TYPE value)
{
   return value;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::round(
   TYPE x)
{
   return x;
}

}
}

#endif
