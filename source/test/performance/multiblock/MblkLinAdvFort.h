/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   F77 external declarations for SAMRAI linear advection example.
 *
 ************************************************************************/

#include <math.h>
#include <signal.h>

extern "C" {

// 2D:

void SAMRAI_F77_FUNC(linadvinit, LINADVINIT) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   double *,
   const int&,
   const double *, const double *);

void SAMRAI_F77_FUNC(linadvinitsine, LINADVINITSINE) (
   const int&, const double *, const double *,
   const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   double *,
   const int&,
   const double *, const double *,
   const double&, const double *);

void SAMRAI_F77_FUNC(initsphere, INITSPHERE) (
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   double *,
   double *,
   const double&, const double&,
   const double *, const double&);

void SAMRAI_F77_FUNC(stabledt, STABLEDT) (
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   double&);

void SAMRAI_F77_FUNC(inittraceflux, INITTRACEFLUX) (
   const int&, const int&,
   const int&, const int&,
   const double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(chartracing0, CHARTRACING0) (
   const double&, const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(chartracing1, CHARTRACING1) (
   const double&, const int&, const int&, const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(fluxcalculation, FLUXCALCULATION) (
   const double&, const int&, const int&,
   const double *,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxcorrec, FLUXCORREC) (
   const double&, const int&, const int&, const int&, const int&,
   const double *,
   const double *, const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(consdiff, CONSDIFF) (
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *, const double *,
   const double *,
   double *);

void SAMRAI_F77_FUNC(getbdry, GETBDRY) (
   const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&,
   const int&,
   const int&,
   const double *, const double&,
   double *,
   const double *, const double *, const int&);

void SAMRAI_F77_FUNC(detectgrad, DETECTGRAD) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void SAMRAI_F77_FUNC(detectshock, DETECTSHOCK) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&, const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void SAMRAI_F77_FUNC(stufprobc, STUFPROBC) (
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&);

// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *, double *);

// 3D:

void SAMRAI_F77_FUNC(linadvinit, LINADVINIT) (
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   double *,
   const int&,
   const double *, const double *);

void SAMRAI_F77_FUNC(linadvinitsine, LINADVINITSINE) (
   const int&, const double *, const double *,
   const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   double *,
   const int&,
   const double *, const double *,
   const double&, const double *);

void SAMRAI_F77_FUNC(initsphere, INITSPHERE) (
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   double *,
   double *,
   const double&, const double&,
   const double *, const double&);

void SAMRAI_F77_FUNC(stabledt, STABLEDT) (
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   double&);

void SAMRAI_F77_FUNC(inittraceflux, INITTRACEFLUX) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   double *, double *, double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(chartracing0, CHARTRACING0) (
   const double&, const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(chartracing1, CHARTRACING1) (
   const double&, const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(chartracing2, CHARTRACING2) (
   const double&, const int&, const int&, const int&, const int&,
   const int&, const int&,
   const int&, const double&, const double&, const int&,
   const double *,
   double *, double *,
   double *, double *,
   double *, double *);

void SAMRAI_F77_FUNC(fluxcalculation, FLUXCALCULATION) (
   const double&, const int&, const int&,
   const int&,
   const double *,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *,
   double *, double *, double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxcorrec2d, FLUXCORREC2D) (
   const double&, const int&, const int&, const int&, const int&,
   const int&, const int&,
   const double *, const double *, const int&,
   const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxcorrec3d, FLUXCORREC3D) (
   const double&, const int&, const int&, const int&, const int&,
   const int&, const int&,
   const double *, const double *,
   const double *,
   const double *, const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(consdiff, CONSDIFF) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const double *, const double *,
   const double *,
   const double *,
   double *);

void SAMRAI_F77_FUNC(getbdry, GETBDRY) (
   const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *, const double&,
   double *,
   const double *, const double *, const int&);

void SAMRAI_F77_FUNC(onethirdstate, ONETHIRDSTATE) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double *, const double *,
   const double *, const double *, const double *,
   double *);

void SAMRAI_F77_FUNC(fluxthird, FLUXTHIRD) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double *, const double *,
   const double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxcorrecjt, FLUXCORRECJT) (
   const double&, const double *, const int&,
   const int&, const int&, const int&, const int&, const int&, const int&,
   const double *, const double *,
   const double *, const double *, const double *,
   double *, double *, double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(detectgrad, DETECTGRAD) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void SAMRAI_F77_FUNC(detectshock, DETECTSHOCK) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double&, const double&,
   const int&, const int&,
   const double *,
   int *, int *);

void SAMRAI_F77_FUNC(stufprobc, STUFPROBC) (
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&);

// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
}
