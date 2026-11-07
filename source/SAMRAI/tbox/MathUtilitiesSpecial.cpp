/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   MathUtilities routines to set up handlers and get
 *                signaling NaNs
 *
 ************************************************************************/

#include "SAMRAI/tbox/MathUtilities.h"

#ifdef HAVE_CMATH_ISNAN
#include <cmath>
#include <math.h>
#else
#include <math.h>
#endif

#include <float.h>
#include <limits.h>
#include <stdlib.h>

#include "SAMRAI/tbox/Complex.h"
#include "Utilities.h"

/*
 * The following lines setup assertion handling headers on the Sun.  If we
 * use Sun's native compiler, just pull in the <sunmath.h> include file.
 * If we are under solaris but use a different compiler (e.g. g++)
 * we have to explicitly define the functions that <sunmath.h> defines,
 * since we don't have access to this file.
 */
#ifdef __SUNPRO_CC
#include <sunmath.h>
#endif

#ifdef __INTEL_COMPILER
// Ignore Intel compiler remarks about non-pointer conversions.
#pragma warning (disable:2259)
#endif

namespace SAMRAI {
namespace tbox {

template<>
dcomplex MathUtilities<dcomplex>::getZero()
{
   return dcomplex(0.0,0.0);
}

template<>
dcomplex MathUtilities<dcomplex>::getOne()
{
   return dcomplex(1.0,0.0);
}

template<>
bool MathUtilities<float>::isNaN(
   const float& value)
{
   int i;
   /* This mess should be fixed when the next C++ standard comes out */
#if defined(HAVE_CMATH_ISNAN)
   i = std::isnan(value);
#elif defined(HAVE_ISNAN)
   i = isnan(value);
#elif defined(HAVE_ISNAND)
   i = __isnanf(value);
#elif defined(HAVE_INLINE_ISNAND)
   i = __inline_isnanf(value);
#else
   i = value != value;
#endif

   return (i != 0) ? true : false;
}

template<>
bool
MathUtilities<double>::isNaN(
   const double& value)
{
   int i;
   /* This mess should be fixed when the next C++ standard comes out */
#if defined(HAVE_CMATH_ISNAN)
   i = std::isnan(value);
#elif defined(HAVE_ISNAN)
   i = isnan(value);
#elif defined(HAVE_ISNAND)
   i = __isnand(value);
#elif defined(HAVE_INLINE_ISNAND)
   i = __inline_isnand(value);
#else
   i = value != value;
#endif

   return (i != 0) ? true : false;
}

template<>
bool
MathUtilities<dcomplex>::isNaN(
   const dcomplex& value)
{
   int i_re;
   int i_im;
#if defined(HAVE_CMATH_ISNAN)
   i_re = std::isnan(real(value));
   i_im = std::isnan(imag(value));
#elif defined(HAVE_ISNAN)
   i_re = isnan(real(value));
   i_im = isnan(imag(value));
#elif defined(HAVE_ISNAND)
   i_re = __isnand(real(value));
   i_im = __isnand(imag(value));
#elif defined(HAVE_INLINE_ISNAND)
   i_re = __inline_isnand(real(value));
   i_im = __inline_isnand(imag(value));
#else
   i_re = real(value) != real(value);
   i_im = imag(value) != imag(value);
#endif

   return ((i_re != 0) || (i_im != 0)) ? true : false;
}

template<>
dcomplex
MathUtilities<dcomplex>::getSignalingNaN()
{
   return dcomplex(std::numeric_limits<double>::signaling_NaN(),
      std::numeric_limits<double>::signaling_NaN());
}

template<>
bool
MathUtilities<float>::equalEps(
   const float& a,
   const float& b)
{
   float absmax = MathUtilities<float>::Max(
         MathUtilities<float>::Abs(a),
         MathUtilities<float>::Abs(b));
   float numerator = MathUtilities<float>::Abs(a - b);
   float denomenator =
      MathUtilities<float>::Max(absmax, MathUtilities<float>::getEpsilon());

   return numerator / denomenator < sqrt(MathUtilities<float>::getEpsilon());
}

template<>
bool
MathUtilities<double>::equalEps(
   const double& a,
   const double& b)
{
   double absmax = MathUtilities<double>::Max(
         MathUtilities<double>::Abs(a),
         MathUtilities<double>::Abs(b));
   double numerator = MathUtilities<double>::Abs(a - b);
   double denomenator =
      MathUtilities<double>::Max(absmax, MathUtilities<double>::getEpsilon());

   return numerator / denomenator < sqrt(MathUtilities<double>::getEpsilon());
}

template<>
bool
MathUtilities<dcomplex>::equalEps(
   const dcomplex& a,
   const dcomplex& b)
{
   double a_re = real(a);
   double a_im = imag(a);
   double b_re = real(b);
   double b_im = imag(b);

   return MathUtilities<double>::equalEps(a_re, b_re) &&
          MathUtilities<double>::equalEps(a_im, b_im);
}

template<>
dcomplex
MathUtilities<dcomplex>::Min(
   dcomplex a,
   dcomplex b)
{
   return norm(a) < norm(b) ? a : b;
}

template<>
dcomplex
MathUtilities<dcomplex>::Max(
   dcomplex a,
   dcomplex b)
{
   return norm(a) > norm(b) ? a : b;
}

template<>
int
MathUtilities<int>::Abs(
   int a)
{
   return a > 0 ? a : -a;
}

template<>
float
MathUtilities<float>::Abs(
   float a)
{
   return a > 0.0 ? a : -a;
}

template<>
double
MathUtilities<double>::Abs(
   double a)
{
   return a > 0.0 ? a : -a;
}

template<>
bool
MathUtilities<bool>::Rand(
   const bool& low,
   const bool& width)
{
   NULL_USE(low);
   NULL_USE(width);
   return mrand48() > 0 ? true : false;
}

template<>
char
MathUtilities<char>::Rand(
   const char& low,
   const char& width)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   return static_cast<char>(static_cast<double>(width) * drand48() + static_cast<double>(low));
}

template<>
int
MathUtilities<int>::Rand(
   const int& low,
   const int& width)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   return static_cast<int>(static_cast<double>(width) * drand48()) + low;
}

template<>
float
MathUtilities<float>::Rand(
   const float& low,
   const float& width)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif
   return static_cast<float>(static_cast<double>(width) * drand48()) + low;
}

template<>
double
MathUtilities<double>::Rand(
   const double& low,
   const double& width)
{
   return width * drand48() + low;
}

template<>
dcomplex
MathUtilities<dcomplex>::Rand(
   const dcomplex& low,
   const dcomplex& width)
{
   double real_part = real(width) * drand48() + real(low);
   double imag_part = imag(width) * drand48() + imag(low);
   return dcomplex(real_part, imag_part);
}

template<class TYPE>
TYPE
round_internal(
   TYPE x)
{
   /* algorithm used from Steven G. Kargl */
   double t;
   if (x >= 0.0) {
      t = ceil(x);
      if (t - x > 0.5)
         t -= 1.0;
      return static_cast<TYPE>(t);
   } else {
      t = ceil(-x);
      if (t + x > 0.5)
         t -= 1.0;
      return static_cast<TYPE>(-t);
   }
}

template float round_internal<float
                              >(
   float x);
template double round_internal<double
                               >(
   double x);

template<>
float
MathUtilities<float>::round(
   float x) {
   return round_internal<float>(x);
}

template<>
double
MathUtilities<double>::round(
   double x) {
   return round_internal<double>(x);
}

}
}
