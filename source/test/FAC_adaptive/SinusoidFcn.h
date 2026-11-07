/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Sinusoidal function functor in FAC solver test.
 *
 ************************************************************************/
#ifndef included_SinusoidFcn_h
#define included_SinusoidFcn_h

#include "SAMRAI/tbox/Dimension.h"

#include <iostream>

using namespace SAMRAI;

/*!
 * @brief Sinusoid function functor.
 *
 * This functor exists for convenience, not efficiency.
 * If you use it heavily, you can improve your efficiency
 * by replacing its usage with lower level codes.
 */
class SinusoidFcn
{

public:
   explicit SinusoidFcn(
      const tbox::Dimension& dim);

   SinusoidFcn(
      const SinusoidFcn& other);

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
   setWaveNumbers(
      const double * npi);

   /*!
    * @brief Set all phase angles.
    *
    * Wave numbers should be given in half-cycles, i.e., 1 -> @f$\pi@f$.
    */
   int
   setPhaseAngles(
      const double * ppi);

   /*!
    * @brief Get wave numbers.
    *
    * Wave numbers are be given in half-cycles, i.e., 1 -> @f$\pi@f$.
    */
   int
   getWaveNumbers(
      double * npi) const;

   /*!
    * @brief Get phase angles.
    *
    * Wave numbers are be given in half-cycles, i.e., 1 -> @f$\pi@f$.
    */
   int
   getPhaseAngles(
      double * ppi) const;

   //@{

   //! @name Differential operations.

   /*!
    * @brief Differentiate and return new function.
    *
    * @return Differentiated function.
    */
   SinusoidFcn
   differentiate(
      unsigned short int x
      ,
      unsigned short int y) const;
   SinusoidFcn
   differentiate(
      unsigned short int x
      ,
      unsigned short int y
      ,
      unsigned short int z) const;

   /*!
    * @brief Differentiate self and return self reference.
    *
    * @return Self reference
    */
   SinusoidFcn&
   differentiateSelf(
      unsigned short int x
      ,
      unsigned short int y);
   SinusoidFcn&
   differentiateSelf(
      unsigned short int x
      ,
      unsigned short int y
      ,
      unsigned short int z);

   //@}

   SinusoidFcn&
   operator = (
      const SinusoidFcn& r);

   //@{
   /*!
    * @name IO operators.
    */
   /*!
    * @brief Input from stream.
    *
    * Stream extract function reads input in the format used
    * by the stream insert function.  Except for allowing for
    * missing coefficients (set to zero), this function requires
    * the correct syntax or the result is undefined.
    *
    * @see SinusoidFcn::ctype
    */
   friend std::istream&
   operator >> (
      std::istream& ci,
      SinusoidFcn& cf);
   /*!
    * @brief Output to stream.
    *
    * Outputs sets of name=double values where the name is one
    * of nx, px, ny, nz or pz, and the double is the value of
    * the coefficient.
    *
    * @see SinusoidFcn::ctype
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const SinusoidFcn& cf);
   //@}

private:
   const tbox::Dimension d_dim;

   //! Amplitude
   double d_amp;
   //! Wave number in half-cycles
   double d_npi[SAMRAI::MAX_DIM_VAL];
   //! Phase shift in half-cycles
   double d_ppi[SAMRAI::MAX_DIM_VAL];

};

#endif
