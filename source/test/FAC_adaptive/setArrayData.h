/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Misc array setting functions in FAC solver test.
 *
 ************************************************************************/
#ifndef include_setArrayData_h
#define include_setArrayData_h

#include "SAMRAI/pdat/MDA_Access.h"
#include "QuarticFcn.h"
#include "SinusoidFcn.h"

void
setArrayDataTo(
   int dim
   ,
   double* ptr
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   const double* coef = 0);

void
setArrayDataTo(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   const double* coef = 0);
void
setArrayDataTo(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   const double* coef = 0);

void
setArrayDataToSinusoidal(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   const double* npi,
   const double* ppi);
void
setArrayDataToSinusoidal(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   const double* npi,
   const double* ppi);

void
setArrayDataToSinusoidalGradient(
   int dim
   ,
   double** g_ptr
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h);

void
setArrayDataToConstant(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   double value);
void
setArrayDataToConstant(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   double value);

void
setArrayDataToLinear(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s,
   const int* lower,
   const int* upper,
   const double* xlo,
   const double* xhi,
   const double* h,
   double a0,
   double ax,
   double ay,
   double axy);
void
setArrayDataToLinear(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s,
   const int* lower,
   const int* upper,
   const double* xlo,
   const double* xhi,
   const double* h,
   double a0,
   double ax,
   double ay,
   double az,
   double axy,
   double axz,
   double ayz,
   double axyz);

void
setArrayDataToScaled(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower,
   const int* upper
   ,
   double factor);
void
setArrayDataToScaled(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower,
   const int* upper
   ,
   double factor);

void
setArrayDataToPerniceExact(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h);
void
setArrayDataToPerniceExact(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h);

void
setArrayDataToPerniceSource(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h);
void
setArrayDataToPerniceSource(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h);

void
setArrayDataToSinusoid(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* beg
   ,
   const int* end
   ,
   const int* ilo,
   const double* xlo,
   const double* h
   ,
   const SinusoidFcn& fcn);
void
setArrayDataToSinusoid(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* beg
   ,
   const int* end
   ,
   const int* ilo,
   const double* xlo,
   const double* h
   ,
   const SinusoidFcn& fcn);

void
setArrayDataToQuartic(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* beg
   ,
   const int* end
   ,
   const int* ilo,
   const double* xlo,
   const double* h
   ,
   const QuarticFcn& fcn);
void
setArrayDataToQuartic(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* beg
   ,
   const int* end
   ,
   const int* ilo,
   const double* xlo,
   const double* h
   ,
   const QuarticFcn& fcn);

#endif  // include_setArrayData_h
