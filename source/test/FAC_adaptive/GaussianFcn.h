/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Gaussian function support for FAC solver tests.
 *
 ************************************************************************/
#ifndef included_GaussianFcn_h
#define included_GaussianFcn_h

#include "SAMRAI/tbox/Dimension.h"
#include <iostream>

using namespace SAMRAI;

/*!
 * @brief Gaussian function functor.
 *
 * Computes the function
 * @f[ e^{ \lambda |r-r_0|^2 } @f]
 *
 * lambda is generally, but not necessarily, negative.
 *
 * This functor exists for convenience, not efficiency.
 * If you use it heavily, you can improve your efficiency
 * by replacing its usage with lower level codes.
 */
class GaussianFcn
{

public:
   GaussianFcn();
   explicit GaussianFcn(
      const tbox::Dimension& dim);
   explicit GaussianFcn(
      const GaussianFcn& other);

   /*!
    * @brief Return the function value.
    */
   double
   operator () (
      double x) const;
   double
   operator () (
      double x,
      double y) const;
   double
   operator () (
      double x,
      double y,
      double z) const;

   /*!
    * @brief Set amplitude.
    */
   int
   setAmplitude(
      const double amp);

   /*!
    * @brief Set all wave numbers.
    *
    * Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
    */
   int
   setLambda(
      const double lambda);

   /*!
    * @brief Set all phase angles.
    *
    * Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
    */
   int
   setCenter(
      const double * center);

   /*!
    * @brief Get amplitude.
    */
   double
   getAmplitude() const;

   /*!
    * @brief Get lambda.
    */
   double
   getLambda() const;

   /*!
    * @brief Get center coordinates.
    */
   int
   getCenter(
      double * center) const;

   GaussianFcn&
   operator = (
      const GaussianFcn& r);

   //@{
   /*!
    * @name IO operators.
    */
   /*!
    * @brief Input from stream.
    *
    * Stream extract function reads input in the format used
    * by the stream insert function (see
    * opertor<< (std::ostream&, GaussianFcn &).
    * Except for allowing for missing centers (set to zero)
    * and lambda (set to 1), this function requires the correct
    * syntax or the result is undefined.
    *
    * @see GaussianFcn::ctype
    */
   friend std::istream&
   operator >> (
      std::istream& ci,
      GaussianFcn& cf);
   /*!
    * @brief Output to stream.
    *
    * Outputs sets of name=double values where the name is one
    * of nx, px, ny, nz or pz, and the double is the value of
    * the coefficient.
    *
    * Format of output with example values is
    * @verbatim
    * { lambda=1.0 cx=0 cy=0.0 cz=0.0 }
    * @endverbatim
    * where cx, cy and cz are the center of the Gaussian function.
    *
    * @see GaussianFcn::ctype
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const GaussianFcn& cf);
   //@}

private:
   const tbox::Dimension d_dim;

   //! Amplitude
   double d_amp;
   //! Center
   double d_center[SAMRAI::MAX_DIM_VAL];
   //! Lambda
   double d_lambda;

};

#endif
