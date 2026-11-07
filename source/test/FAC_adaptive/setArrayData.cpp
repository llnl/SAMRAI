/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Misc array setting functions in FAC solver test.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "setArrayData.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/tbox/Utilities.h"
#include <math.h>


/*!
 * \file
 * \brief Set array data when array is given as a pointer or
 * MultiDimArrayAccess object.
 *
 * There is no SAMRAI-specific code in this file.
 */

void setArrayDataToConstant(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   double value) {
   NULL_USE(xlo);
   NULL_USE(xhi);
   NULL_USE(h);

   for (int j = lower[1]; j <= upper[1]; ++j) {
      for (int i = lower[0]; i <= upper[0]; ++i) {
         s(i, j) = value;
      }
   }
}
void setArrayDataToConstant(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h
   ,
   double value) {
   NULL_USE(xlo);
   NULL_USE(xhi);
   NULL_USE(h);

   for (int k = lower[2]; k <= upper[2]; ++k) {
      for (int j = lower[1]; j <= upper[1]; ++j) {
         for (int i = lower[0]; i <= upper[0]; ++i) {
            s(i, j, k) = value;
         }
      }
   }
}

void setArrayDataTo(
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
   const double* coef) {
   if (dim == 2) {
      MDA_Access<double, 2, MDA_OrderColMajor<2> > s(ptr, lower, upper);
      setArrayDataTo(s, lower, upper, xlo, xhi, h, coef);
   } else if (dim == 3) {
      MDA_Access<double, 3, MDA_OrderColMajor<3> > s(ptr, lower, upper);
      setArrayDataTo(s, lower, upper, xlo, xhi, h, coef);
   }
}

void setArrayDataTo(
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
   const double* coef) {
   NULL_USE(xhi);

   const double ucoef[2] = { 1., 1. };
   if (coef == 0) coef = ucoef;
   for (int j = lower[1]; j <= upper[1]; ++j) {
      double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
      for (int i = lower[0]; i <= upper[0]; ++i) {
         double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
         s(i, j) = coef[0] * x + coef[1] * y;
      }
   }
}
void setArrayDataTo(
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
   const double* coef) {
   NULL_USE(xhi);

   const double ucoef[3] = { 1., 1., 1. };
   if (coef == 0) coef = ucoef;
   for (int k = lower[2]; k <= upper[2]; ++k) {
      double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            s(i, j, k) = coef[0] * x + coef[1] * y + coef[2] * z;
         }
      }
   }
}

void setArrayDataToSinusoidal(
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
   const double* ppi) {
   NULL_USE(xhi);

   double nx = npi[0], px = ppi[0];
   double ny = npi[1], py = ppi[1];
   for (int j = lower[1]; j <= upper[1]; ++j) {
      double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
      double siny = sin(M_PI * (ny * y + py));
      for (int i = lower[0]; i <= upper[0]; ++i) {
         double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
         double sinx = sin(M_PI * (nx * x + px));
         s(i, j) = sinx * siny;
      }
   }
}
void setArrayDataToSinusoidal(
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
   const double* ppi) {
   NULL_USE(xhi);

   double nx = npi[0], px = ppi[0];
   double ny = npi[1], py = ppi[1];
   double nz = npi[2], pz = ppi[2];
   for (int k = lower[2]; k <= upper[2]; ++k) {
      double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
      double sinz = sin(M_PI * (nz * z + pz));
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         double siny = sin(M_PI * (ny * y + py));
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            double sinx = sin(M_PI * (nx * x + px));
            s(i, j, k) = sinx * siny * sinz;
         }
      }
   }
}

void setArrayDataToSinusoidalGradient(
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
   const double* h) {
   NULL_USE(xhi);
   NULL_USE(h);
   if (dim == 2) {
      double* gx_ptr = g_ptr[0];
      MDA_Access<double, 2, MDA_OrderColMajor<2> > gx(gx_ptr, lower, upper);
      double* gy_ptr = g_ptr[1];
      MDA_Access<double, 2, MDA_OrderColMajor<2> > gy(gy_ptr, lower, upper);
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         double siny = sin(2 * M_PI * y);
         double cosy = cos(2 * M_PI * y);
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            double sinx = sin(2 * M_PI * x);
            double cosx = cos(2 * M_PI * x);
            gx(i, j) = 2 * M_PI * cosx * siny;
            gy(i, j) = sinx * 2 * M_PI * cosy;
         }
      }
   } else if (dim == 3) {
      double* gx_ptr = g_ptr[0];
      MDA_Access<double, 3, MDA_OrderColMajor<3> > gx(gx_ptr, lower, upper);
      double* gy_ptr = g_ptr[1];
      MDA_Access<double, 3, MDA_OrderColMajor<3> > gy(gy_ptr, lower, upper);
      double* gz_ptr = g_ptr[2];
      MDA_Access<double, 3, MDA_OrderColMajor<3> > gz(gz_ptr, lower, upper);
      for (int k = lower[2]; k <= upper[2]; ++k) {
         double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
         double sinz = sin(2 * M_PI * z);
         double cosz = cos(2 * M_PI * z);
         for (int j = lower[1]; j <= upper[1]; ++j) {
            double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
            double siny = sin(2 * M_PI * y);
            double cosy = cos(2 * M_PI * y);
            for (int i = lower[0]; i <= upper[0]; ++i) {
               double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
               double sinx = sin(2 * M_PI * x);
               double cosx = cos(2 * M_PI * x);
               gx(i, j, k) = 2 * M_PI * cosx * siny * sinz;
               gy(i, j, k) = sinx * 2 * M_PI * cosy * sinz;
               gz(i, j, k) = sinx * cosy * 2 * M_PI * cosz;
            }
         }
      }
   }
}

void setArrayDataToLinear(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s,
   const int* lower,
   const int* upper,
   const double* xlo,
   const double* xhi,
   const double* h,
   double a0,
   double ax,
   double ay,
   double axy) {
   NULL_USE(xhi);

   for (int j = lower[1]; j <= upper[1]; ++j) {
      double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
      for (int i = lower[0]; i <= upper[0]; ++i) {
         double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
         s(i, j) = a0 + ax * x + ay * y + axy * x * y;
      }
   }
}
void setArrayDataToLinear(
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
   double axyz) {
   NULL_USE(xhi);

   for (int k = lower[2]; k <= upper[2]; ++k) {
      double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            s(i, j, k) = a0 + ax * x + ay * y + az * z
               + axy * x * y + axz * x * z + ayz * y * z + axyz * x * y * z;
         }
      }
   }
}

void setArrayDataToScaled(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower,
   const int* upper
   ,
   double factor) {
   for (int j = lower[1]; j <= upper[1]; ++j) {
      for (int i = lower[0]; i <= upper[0]; ++i) {
         s(i, j) *= factor;
      }
   }
}
void setArrayDataToScaled(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower,
   const int* upper
   ,
   double factor) {
   for (int k = lower[2]; k <= upper[2]; ++k) {
      for (int j = lower[1]; j <= upper[1]; ++j) {
         for (int i = lower[0]; i <= upper[0]; ++i) {
            s(i, j, k) *= factor;
         }
      }
   }
}

/*!
 * \brief Set array to Michael's exact solution
 */
void setArrayDataToPerniceExact(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h) {
   NULL_USE(xhi);

   for (int j = lower[1]; j <= upper[1]; ++j) {
      double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
      for (int i = lower[0]; i <= upper[0]; ++i) {
         double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
         s(i, j) = x * (1 - x) * y * (1 - y);
      }
   }
}
void setArrayDataToPerniceExact(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h) {
   NULL_USE(xhi);

   for (int k = lower[2]; k <= upper[2]; ++k) {
      double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            s(i, j, k) = x * (1 - x) * y * (1 - y) * z * (1 - z);
         }
      }
   }
}

/*!
 * \brief Set array to Michael's source
 */
void setArrayDataToPerniceSource(
   MDA_Access<double, 2, MDA_OrderColMajor<2> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h) {
   NULL_USE(xhi);

   for (int j = lower[1]; j <= upper[1]; ++j) {
      double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
      for (int i = lower[0]; i <= upper[0]; ++i) {
         double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
         s(i, j) = -2 * (x * (1 - x) + y * (1 - y));
      }
   }
}
void setArrayDataToPerniceSource(
   MDA_Access<double, 3, MDA_OrderColMajor<3> >& s
   ,
   const int* lower
   ,
   const int* upper
   ,
   const double* xlo,
   const double* xhi,
   const double* h) {
   NULL_USE(xhi);

   for (int k = lower[2]; k <= upper[2]; ++k) {
      double z = xlo[2] + h[2] * (k - lower[2] + 0.5);
      for (int j = lower[1]; j <= upper[1]; ++j) {
         double y = xlo[1] + h[1] * (j - lower[1] + 0.5);
         for (int i = lower[0]; i <= upper[0]; ++i) {
            double x = xlo[0] + h[0] * (i - lower[0] + 0.5);
            s(i, j, k) = -2 * (x * (1 - x) + y * (1 - y) + z * (1 - z));
         }
      }
   }
}

void setArrayDataToSinusoid(
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
   const SinusoidFcn& fcn)
{
   for (int j = beg[1]; j <= end[1]; ++j) {
      double y = xlo[1] + h[1] * (j - ilo[1] + 0.5);
      for (int i = beg[0]; i <= end[0]; ++i) {
         double x = xlo[0] + h[0] * (i - ilo[0] + 0.5);
         s(i, j) = fcn(x, y);
      }
   }
}
void setArrayDataToSinusoid(
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
   const SinusoidFcn& fcn)
{
   for (int k = beg[2]; k <= end[2]; ++k) {
      double z = xlo[2] + h[2] * (k - ilo[2] + 0.5);
      for (int j = beg[1]; j <= end[1]; ++j) {
         double y = xlo[1] + h[1] * (j - ilo[1] + 0.5);
         for (int i = beg[0]; i <= end[0]; ++i) {
            double x = xlo[0] + h[0] * (i - ilo[0] + 0.5);
            s(i, j, k) = fcn(x, y, z);
         }
      }
   }
}

void setArrayDataToQuartic(
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
   const QuarticFcn& fcn)
{
   for (int j = beg[1]; j <= end[1]; ++j) {
      double y = xlo[1] + h[1] * (j - ilo[1] + 0.5);
      for (int i = beg[0]; i <= end[0]; ++i) {
         double x = xlo[0] + h[0] * (i - ilo[0] + 0.5);
         s(i, j) = fcn(x, y);
      }
   }
}
void setArrayDataToQuartic(
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
   const QuarticFcn& fcn)
{
   for (int k = beg[2]; k <= end[2]; ++k) {
      double z = xlo[2] + h[2] * (k - beg[2] + 0.5);
      for (int j = beg[1]; j <= end[1]; ++j) {
         double y = xlo[1] + h[1] * (j - ilo[1] + 0.5);
         for (int i = beg[0]; i <= end[0]; ++i) {
            double x = xlo[0] + h[0] * (i - ilo[0] + 0.5);
            s(i, j, k) = fcn(x, y, z);
         }
      }
   }
}
