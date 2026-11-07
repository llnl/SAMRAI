/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Array specializations
 *
 ************************************************************************/

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

template<>
Array<bool>::Array(
   const int n,
   const bool& default_value)
{
   if (n > 0) {
      d_objects = reinterpret_cast<bool *>(malloc(sizeof(bool) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         d_objects[i] = default_value;
      }

   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<>
Array<char>::Array(
   const int n,
   const char& default_value)
{
   if (n > 0) {
      d_objects = reinterpret_cast<char *>(malloc(sizeof(char) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         d_objects[i] = default_value;
      }

   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<>
Array<int>::Array(
   const int n,
   const int& default_value)
{
   if (n > 0) {
      d_objects = reinterpret_cast<int *>(malloc(sizeof(int) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         d_objects[i] = default_value;
      }

   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<>
Array<float>::Array(
   const int n,
   const float& default_value)
{
   if (n > 0) {
      d_objects = reinterpret_cast<float *>(malloc(sizeof(float) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         d_objects[i] = default_value;
      }

   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<>
Array<double>::Array(
   const int n,
   const double& default_value)
{
   if (n > 0) {
      d_objects = reinterpret_cast<double *>(malloc(sizeof(double) * n));
      d_counter = new ReferenceCounter;
      d_elements = n;

      for (int i = 0; i < d_elements; ++i) {
         d_objects[i] = default_value;
      }

   } else {
      d_objects = 0;
      d_counter = 0;
      d_elements = 0;
   }
}

template<>
void
Array<bool>::deleteObjects()
{
   if (d_objects) {
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<>
void
Array<char>::deleteObjects()
{
   if (d_objects) {
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<>
void
Array<int>::deleteObjects()
{
   if (d_objects) {
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<>
void
Array<float>::deleteObjects()
{
   if (d_objects) {
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

template<>
void
Array<double>::deleteObjects()
{
   if (d_objects) {
      free(reinterpret_cast<char *>(d_objects));
      delete d_counter;
   }

   d_objects = 0;
   d_counter = 0;
   d_elements = 0;
}

}
}
#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
