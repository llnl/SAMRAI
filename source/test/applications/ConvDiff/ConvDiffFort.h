/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   F77 external declarations for SAMRAI Heat Equation example.
 *
 ************************************************************************/

#include <math.h>
#include <signal.h>

// Function argument list interfaces
extern "C" {

void SAMRAI_F77_FUNC(initsphere2d, INITSPHERE2D) (
   const double *, const double *, const double *,
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const int&);

void SAMRAI_F77_FUNC(initsphere3d, INITSPHERE3D) (
   const double *, const double *, const double *,
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const int&);

void SAMRAI_F77_FUNC(computerhs2d, COMPUTERHS2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const double *, // dx
   const double *, // d_convection_coeff
   const double&,  // d_diffusion_coeff
   const double&,  // d_source_coeff
   const double *, // prim_var_updated
   double *,       // function_eval
   const int&);    // NEQU

void SAMRAI_F77_FUNC(computerhs3d, COMPUTERHS3D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, // dx
   const double *, // d_convection_coeff
   const double&,  // d_diffusion_coeff
   const double&,  // d_source_coeff
   const double *, // prim_var_updated
   double *,       // function_eval
   const int&);    // NEQU

void SAMRAI_F77_FUNC(rkstep2d, RKSTEP2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const double&, const double&, const double&, const double&,
   const double *,
   const double&,
   const double&,
   double *,
   const double *,
   const double *,
   const int&);

void SAMRAI_F77_FUNC(rkstep3d, RKSTEP3D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double&, const double&, const double&, const double&,
   const double *,
   const double&,
   const double&,
   double *,
   const double *,
   const double *,
   const int&);

void SAMRAI_F77_FUNC(tagcells2d, TAGCELLS2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   int *,
   const double *,
   const int&,
   const double *,
   const int&);

void SAMRAI_F77_FUNC(tagcells3d, TAGCELLS3D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   int *,
   const double *,
   const int&,
   const double *,
   const int&);

}
