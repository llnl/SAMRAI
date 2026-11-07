/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   IEEE routines to set up handlers and get signaling NaNs
 *
 ************************************************************************/

#ifndef included_tbox_IEEE
#define included_tbox_IEEE

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * Class IEEE is a utility providing rotuines for managing IEEE trap
 * handlers and data set to signaling NaNs.  Signaling NaNs force
 * a trap if they are used in a numerical operation, so they are a
 * useful way to track uninitialized floating point data.  Signaling
 * NaN's may only be used for double and float data (and the real
 * ans imaginary parts of dcomplex data) and so operations are
 * provided here for those types only.
 *
 * IMPORTANT: To properly trap operations based on signaling NaN values,
 *            the routine IEEE::setupFloatingPointExceptionHandlers()
 *            must be called.  This is normally done in the
 *            SAMRAIManager::startup() routine.
 *
 * Note that all operations provided by this class (except for setting
 * up exception handling) are implemented in @see MathUtilities.
 * Operations are provided by this class since it is not templated on
 * data type and so calling the operations provided here may be easier
 * in some cases, such as in codes built based on earlier versions
 * of SAMRAI.  See the MathUtilities header file for details
 * about the routines.
 *
 * @see MathUtilities
 */

struct IEEE {
   /*!
    * Set up IEEE exception handlers so that normal IEEE exceptions will
    * cause a program abort.  This is useful for tracking down errors.
    * Note, however, that this may cause problems if your code relies on
    * IEEE underflow or overflow traps.
    */
   static void
   setupFloatingPointExceptionHandlers();

   /*!
    * Get the IEEE float signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static float
   getSignalingFloatNaN()
   {
      return MathUtilities<float>::getSignalingNaN();
   }

   /*!
    * Get the IEEE double signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static double
   getSignalingNaN()
   {
      return MathUtilities<double>::getSignalingNaN();
   }

   /*!
    * Get the dcomplex value with real and imaginary parts set to the
    * IEEE double signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static dcomplex
   getSignalingComplexNaN()
   {
      return MathUtilities<dcomplex>::getSignalingNaN();
   }

   /*!
    * Set supplied float value to the signaling NaN.
    */
   static void
   setNaN(
      float& f)
   {
      f = MathUtilities<float>::getSignalingNaN();
   }

   /*!
    * Set supplied double value to the signaling NaN.
    */
   static void
   setNaN(
      double& d)
   {
      d = MathUtilities<double>::getSignalingNaN();
   }

   /*!
    * Set real and imaginary parts of supplied dcomplex value to the
    * double signaling NaN.
    */
   static void
   setNaN(
      dcomplex& dc)
   {
      dc = MathUtilities<dcomplex>::getSignalingNaN();
   }

   /*!
    * Initialize a vector of floats to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeVectorToSignalingNaN(
      std::vector<float>& vector)
   {
      MathUtilities<float>::setVectorToSignalingNaN(vector);
   }

   /*!
    * Initialize a vector of doubles to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeVectorToSignalingNaN(
      std::vector<double>& vector)
   {
      MathUtilities<double>::setVectorToSignalingNaN(vector);
   }

   /*!
    * Initialize a vector of dcomplex to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeVectorToSignalingNaN(
      std::vector<dcomplex>& vector)
   {
      MathUtilities<dcomplex>::setVectorToSignalingNaN(vector);
   }

   /*!
    * Initialize an array of floats to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      float* array,
      int n = 1)
   {
      MathUtilities<float>::setArrayToSignalingNaN(array, n);
   }

   /*!
    * Initialize an array of doubles to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      double* array,
      int n = 1)
   {
      MathUtilities<double>::setArrayToSignalingNaN(array, n);
   }

   /*!
    * Initialize an array of dcomplex to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      dcomplex* array,
      int n = 1)
   {
      MathUtilities<dcomplex>::setArrayToSignalingNaN(array, n);
   }

   /*!
    * Return true if the supplied float value is NaN; else, false.
    */
   static bool
   isNaN(
      const float& f)
   {
      return MathUtilities<float>::isNaN(f);
   }

   /*!
    * Return true if the supplied double value is NaN; else, false.
    */
   static bool
   isNaN(
      const double& d)
   {
      return MathUtilities<double>::isNaN(d);
   }

   /*!
    * Return true if if either real and imaginary part of the supplied
    * dcomplex value is NaN; else, false.
    */
   static bool
   isNaN(
      const dcomplex& dc)
   {
      return MathUtilities<dcomplex>::isNaN(dc);
   }
};

}
}

#endif
