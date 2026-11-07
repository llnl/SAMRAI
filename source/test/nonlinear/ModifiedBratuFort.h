/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fortran function declarations for modified-Bratu problem
 *
 ************************************************************************/

// Header file to define interfaces between Fortran and C++ code for
// the ModifiedBratuProblem class.

#ifndef included_modifiedBratuFort
#define included_modifiedBratuFort

extern "C"
{

// 1D

#define FORT_FILL1D fill1d_
#define FORT_EVALBRATU1D evalbratu1d_
#define FORT_EVALDIFFUSION1D evaldiffusioncoef1d_
#define FORT_EVALEXPONENTIAL1D evalexponential1d_
#define FORT_EVALFACEFLUXES1D evalfacefluxes1d_
#define FORT_EVALSOURCE1D evalsource1d_
#define FORT_EWBCFLUXFIX1D ewbcfluxfix1d_
#define FORT_NSBCFLUXFIX1D nsbcfluxfix1d_
#define FORT_TBBCFLUXFIX1D tbbcfluxfix1d_
#define FORT_EWFLUXCOPY1D ewfluxcopy1d_
#define FORT_NSFLUXCOPY1D nsfluxcopy1d_
#define FORT_TBFLUXCOPY1D tbfluxcopy1d_
#define FORT_BRATUJV1D bratujv1d_
#define FORT_SETBC1D setbc1d_
#define FORT_ERROR1D error1d_
#define FORT_EVALF1D evalf1d_
#define FORT_PROLONG1D prolong1d_

void
FORT_FILL1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const int&);

void
FORT_EVALBRATU1D(
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EVALDIFFUSION1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALEXPONENTIAL1D(
   const int&,
   const int&,
   const double *,
   const double&,
   const double *);

void
FORT_EVALFACEFLUXES1D(
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALSOURCE1D(
   const int&,
   const int&,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EWBCFLUXFIX1D(
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_EWFLUXCOPY1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_BRATUJV1D(
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_SETBC1D(
   const int&,
   const int&,
   const int&,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_ERROR1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double&,
   const double&);

void
FORT_EVALF1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_VAXPY1D(
   const int&,
   const int&,
   const double *,
   const double *,
   const double *);

void SAMRAI_F77_FUNC(compfacdiag1d, COMPFACDIAG1D) (const int&, const int&,
   const double&, const double&,
   const double *, const double *,
   double *);

void SAMRAI_F77_FUNC(compfacoffdiag1d, COMPFACOFFDIAG1D) (const int&, const int&,
   const double&, const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compdiffcoef1d, COMPDIFFCOEF1D) (const int&, const int&,
   const double *, const double *,
   const double&,
   double *);

void SAMRAI_F77_FUNC(compexpu1d, COMPEXPU1D) (const int&, const int&,
   const int&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrc1d, COMPSRC1D) (const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrcderv1d, COMPSRCDERV1D) (const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsideflux1d, COMPSIDEFLUX1D) (const int&, const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   double *);

void SAMRAI_F77_FUNC(fluxbdryfix1d, FLUXBDRYFIX1D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   double *);

void SAMRAI_F77_FUNC(fluxcopy01d, FLUXCOPY01D) (const int&, const int&,
   const int&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compresidual1d, COMPRESIDUAL1D) (const int&, const int&,
   const int&,
   const double *, const double&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   double *);

// Bonus function

void SAMRAI_F77_FUNC(adjcrsfineoffdiag1d, ADJCRSFINEOFFDIAG1D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   double *);

// 2D

#define FORT_FILL2D fill2d_
#define FORT_EVALBRATU2D evalbratu2d_
#define FORT_EVALDIFFUSION2D evaldiffusioncoef2d_
#define FORT_EVALEXPONENTIAL2D evalexponential2d_
#define FORT_EVALFACEFLUXES2D evalfacefluxes2d_
#define FORT_EVALSOURCE2D evalsource2d_
#define FORT_EWBCFLUXFIX2D ewbcfluxfix2d_
#define FORT_NSBCFLUXFIX2D nsbcfluxfix2d_
#define FORT_TBBCFLUXFIX2D tbbcfluxfix2d_
#define FORT_EWFLUXCOPY2D ewfluxcopy2d_
#define FORT_NSFLUXCOPY2D nsfluxcopy2d_
#define FORT_TBFLUXCOPY2D tbfluxcopy2d_
#define FORT_BRATUJV2D bratujv2d_
#define FORT_SETBC2D setbc2d_
#define FORT_ERROR2D error2d_
#define FORT_EVALF2D evalf2d_
#define FORT_PROLONG2D prolong2d_

void
FORT_FILL2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const int&);

void
FORT_EVALBRATU2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EVALDIFFUSION2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALEXPONENTIAL2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double&,
   const double *);

void
FORT_EVALFACEFLUXES2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALSOURCE2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EWBCFLUXFIX2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_NSBCFLUXFIX2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_EWFLUXCOPY2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_NSFLUXCOPY2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_BRATUJV2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_SETBC2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_ERROR2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double&,
   const double&);

void
FORT_EVALF2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_VAXPY2D(
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *);

/* These functions are in FACjacobian.m4 */

void SAMRAI_F77_FUNC(compjv2d, COMPJV2D) (
   const int& ifirst0, const int& ilast0,
   const int& ifirst1, const int& ilast1,
   const int& gwc,
   const double* diag,
   const double* flux0, const double* flux1,
   const double* v,
   const double* dx,
   const double& dt,
   double* jv);

void SAMRAI_F77_FUNC(compfacdiag2d, COMPFACDIAG2D) (const int&, const int&,
   const int&, const int&,
   const double&, const double&,
   const double *, const double *,
   double *);

void SAMRAI_F77_FUNC(compfacoffdiag2d, COMPFACOFFDIAG2D) (const int&, const int&,
   const int&, const int&,
   const double&, const double&,
   const double *, const double *,
   double *, double *);

void SAMRAI_F77_FUNC(compdiffcoef2d, COMPDIFFCOEF2D) (const int&, const int&,
   const int&, const int&,
   const double *, const double *,
   const double&,
   double *, double *);

void SAMRAI_F77_FUNC(compexpu2d, COMPEXPU2D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrc2d, COMPSRC2D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrcderv2d, COMPSRCDERV2D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsideflux2d, COMPSIDEFLUX2D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   const double *, const double *,
   const double *,
   double *, double *);

void SAMRAI_F77_FUNC(fluxbdryfix2d, FLUXBDRYFIX2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   double *, double *);

void SAMRAI_F77_FUNC(fluxcopy02d, FLUXCOPY02D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(fluxcopy12d, FLUXCOPY12D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compresidual2d, COMPRESIDUAL2D) (const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *, const double *,
   double *);

// 3D

#define FORT_FILL3D fill3d_
#define FORT_EVALBRATU3D evalbratu3d_
#define FORT_EVALDIFFUSION3D evaldiffusioncoef3d_
#define FORT_EVALEXPONENTIAL3D evalexponential3d_
#define FORT_EVALFACEFLUXES3D evalfacefluxes3d_
#define FORT_EVALSOURCE3D evalsource3d_
#define FORT_EWBCFLUXFIX3D ewbcfluxfix3d_
#define FORT_NSBCFLUXFIX3D nsbcfluxfix3d_
#define FORT_TBBCFLUXFIX3D tbbcfluxfix3d_
#define FORT_EWFLUXCOPY3D ewfluxcopy3d_
#define FORT_NSFLUXCOPY3D nsfluxcopy3d_
#define FORT_TBFLUXCOPY3D tbfluxcopy3d_
#define FORT_BRATUJV3D bratujv3d_
#define FORT_SETBC3D setbc3d_
#define FORT_ERROR3D error3d_
#define FORT_EVALF3D evalf3d_
#define FORT_PROLONG3D prolong3d_

void
FORT_FILL3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const int&);

void
FORT_EVALBRATU3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EVALDIFFUSION3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALEXPONENTIAL3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double&,
   const double *);

void
FORT_EVALFACEFLUXES3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_EVALSOURCE3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_EWBCFLUXFIX3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_NSBCFLUXFIX3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_TBBCFLUXFIX3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_EWFLUXCOPY3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_NSFLUXCOPY3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_TBFLUXCOPY3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const int&);

void
FORT_BRATUJV3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *,
   const double&,
   const double *);

void
FORT_SETBC3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const int *,
   const int *,
   const int&);

void
FORT_ERROR3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double&,
   const double *,
   const double *,
   const double *,
   const double&,
   const double&,
   const double&);

void
FORT_EVALF3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *,
   const double *);

void
FORT_VAXPY3D(
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const int&,
   const double *,
   const double *,
   const double *);

/* These functions are in FACjacobian.m4 */

void SAMRAI_F77_FUNC(compjv3d, COMPJV3D) (
   const int& ifirst0, const int& ilast0,
   const int& ifirst1, const int& ilast1,
   const int& ifirst2, const int& ilast2,
   const int& gwc,
   const double* diag,
   const double* flux0, const double* flux1, const double* flux2,
   const double* v,
   const double* dx,
   const double& dt,
   double* jv);

void SAMRAI_F77_FUNC(compfacdiag3d, COMPFACDIAG3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&, const double&,
   const double *, const double *,
   double *);
void SAMRAI_F77_FUNC(compfacoffdiag3d, COMPFACOFFDIAG3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&, const double&,
   const double *, const double *, const double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(compdiffcoef3d, COMPDIFFCOEF3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *, const double *,
   const double&,
   double *, double *, double *);

void SAMRAI_F77_FUNC(compexpu3d, COMPEXPU3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrc3d, COMPSRC3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsrcderv3d, COMPSRCDERV3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double *, const double&,
   const double&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compsideflux3d, COMPSIDEFLUX3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   const double *, const double *, const double *,
   const double *,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxbdryfix3d, FLUXBDRYFIX3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   const double&,
   double *, double *, double *);

void SAMRAI_F77_FUNC(fluxcopy03d, FLUXCOPY03D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   double *);
void SAMRAI_F77_FUNC(fluxcopy13d, FLUXCOPY13D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *,
   double *);

void SAMRAI_F77_FUNC(compresidual3d, COMPRESIDUAL3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const double *, const double&,
   const double *,
   const double *,
   const double *,
   const double *,
   const double *, const double *, const double *,
   double *);

}

#endif
