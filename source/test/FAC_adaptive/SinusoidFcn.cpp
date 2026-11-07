/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Sinusoidal function functor in FAC solver test.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SinusoidFcn.h"
#include <math.h>
#include "SAMRAI/tbox/Utilities.h"

/*
 * Temporary fix for g++ lacking instantiations when --no-implicit-templates
 * is used (by SAMRAI)
 */
#define fill_n(p, n, v) { size_t _i; for (_i = 0; _i < n; ++_i) p[_i] = v; }
#define copy_n(s, n, d) { size_t _i; for (_i = 0; _i < n; ++_i) d[_i] = s[_i]; }

SinusoidFcn::SinusoidFcn(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_amp(1.0)
{
   fill_n(d_npi, d_dim.getValue(), 0.0)
   fill_n(d_ppi, d_dim.getValue(), 0.0)
}

SinusoidFcn::SinusoidFcn(
   const SinusoidFcn& other):
   d_dim(other.d_dim),
   d_amp(other.d_amp)
{
   copy_n(other.d_npi, d_dim.getValue(), d_npi)
   copy_n(other.d_ppi, d_dim.getValue(), d_ppi)
}

int SinusoidFcn::setAmplitude(
   const double amp) {
   d_amp = amp;
   return 0;
}

int SinusoidFcn::setWaveNumbers(
   const double* npi) {
   for (size_t i = 0; i < d_dim.getValue(); ++i) d_npi[i] = npi[i];
   return 0;
}

int SinusoidFcn::getWaveNumbers(
   double* npi) const {
   for (size_t i = 0; i < d_dim.getValue(); ++i) npi[i] = d_npi[i];
   return 0;
}

int SinusoidFcn::setPhaseAngles(
   const double* ppi) {
   for (size_t i = 0; i < d_dim.getValue(); ++i) d_ppi[i] = ppi[i];
   return 0;
}

int SinusoidFcn::getPhaseAngles(
   double* ppi) const {
   for (size_t i = 0; i < d_dim.getValue(); ++i) ppi[i] = d_ppi[i];
   return 0;
}

double SinusoidFcn::operator () (
   double x) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(1));
   double rval;
   rval = d_amp
      * sin(M_PI * (d_npi[0] * x + d_ppi[0]));
   return rval;
}
double SinusoidFcn::operator () (
   double x,
   double y) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(2));
   double rval;
   rval = d_amp
      * sin(M_PI * (d_npi[0] * x + d_ppi[0]))
      * sin(M_PI * (d_npi[1] * y + d_ppi[1]));
   return rval;
}
double SinusoidFcn::operator () (
   double x,
   double y,
   double z) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(3));
   double rval;
   rval = d_amp
      * sin(M_PI * (d_npi[0] * x + d_ppi[0]))
      * sin(M_PI * (d_npi[1] * y + d_ppi[1]))
      * sin(M_PI * (d_npi[2] * z + d_ppi[2]));
   return rval;
}

SinusoidFcn& SinusoidFcn::differentiateSelf(
   unsigned short int x
   ,
   unsigned short int y)
{
   /*
    * Since differentiation commutes,
    * simply differentiate one direction at a time.
    */
   // Differentiate in x direction.
   for ( ; x > 0; --x) {
      d_amp *= M_PI * d_npi[0];
      d_ppi[0] += 0.5;
   }
   while (d_ppi[0] > 2) d_ppi[0] -= 2.0;
   // Differentiate in y direction.
   for ( ; y > 0; --y) {
      d_amp *= M_PI * d_npi[1];
      d_ppi[1] += 0.5;
   }
   while (d_ppi[1] > 2) d_ppi[1] -= 2.0;
   return *this;
}

SinusoidFcn& SinusoidFcn::differentiateSelf(
   unsigned short int x
   ,
   unsigned short int y
   ,
   unsigned short int z)
{
   /*
    * Since differentiation commutes,
    * simply differentiate one direction at a time.
    */
   // Differentiate in x direction.
   for ( ; x > 0; --x) {
      d_amp *= M_PI * d_npi[0];
      d_ppi[0] += 0.5;
   }
   while (d_ppi[0] > 2) d_ppi[0] -= 2.0;
   // Differentiate in y direction.
   for ( ; y > 0; --y) {
      d_amp *= M_PI * d_npi[1];
      d_ppi[1] += 0.5;
   }
   while (d_ppi[1] > 2) d_ppi[1] -= 2.0;
   // Differentiate in z direction.
   for ( ; z > 0; --z) {
      d_amp *= M_PI * d_npi[2];
      d_ppi[2] += 0.5;
   }
   while (d_ppi[2] > 2) d_ppi[2] -= 2.0;
   return *this;
}

SinusoidFcn SinusoidFcn::differentiate(
   unsigned short int x
   ,
   unsigned short int y) const
{
   SinusoidFcn rval(*this);
   rval.differentiateSelf(x, y);
   return rval;
}

SinusoidFcn SinusoidFcn::differentiate(
   unsigned short int x
   ,
   unsigned short int y
   ,
   unsigned short int z) const
{
   SinusoidFcn rval(*this);
   rval.differentiateSelf(x, y, z);
   return rval;
}

#define EAT_WS(s) \
   { while (s.peek() == ' '                     \
            || s.peek() == '\t'                    \
            || s.peek() == '\n') { s.get(); } }

std::istream& operator >> (
   std::istream& ci,
   SinusoidFcn& sf) {
   fill_n(sf.d_npi, sf.d_dim.getValue(), 0.0)
   fill_n(sf.d_ppi, sf.d_dim.getValue(), 0.0)
   char dummy, name[2];
   EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
   ci >> dummy;
   TBOX_ASSERT(dummy == '{');
   EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
   while (ci.peek() != '}') {
      ci.read(name, 2);
      if (name[0] == 'a' && name[1] == 'm') {
         ci >> dummy;
         TBOX_ASSERT(dummy == 'p');
         ci >> dummy;
         TBOX_ASSERT(dummy == '=');
         ci >> sf.d_amp;
      } else {
         ci >> dummy;
         TBOX_ASSERT(dummy == '=');
         double * data(name[0] == 'n' ? sf.d_npi : sf.d_ppi);
         int dim(name[1] == 'x' ? 0 :
                 name[1] == 'y' ? 1 :
                 name[1] == 'z' ? 2 : 3);
         TBOX_ASSERT(dim < sf.d_dim.getValue());
         ci >> data[dim];
      }
      EAT_WS(ci) // ci >> std::skipws; // ci.ipfx(0);
   }
   return ci;
}

SinusoidFcn& SinusoidFcn::operator = (
   const SinusoidFcn& r) {
   TBOX_ASSERT(d_dim == r.d_dim);
   d_amp = r.d_amp;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_npi[i] = r.d_npi[i];
      d_ppi[i] = r.d_ppi[i];
   }
   return *this;
}

std::ostream& operator << (
   std::ostream& co,
   const SinusoidFcn& sf) {
   co << "{ amp=" << sf.d_amp;
   co << " nx=" << sf.d_npi[0] << " px=" << sf.d_npi[0];
   if (sf.d_dim >= tbox::Dimension(2)) {
      co << " ny=" << sf.d_npi[1] << " py=" << sf.d_npi[1];
   }
   if (sf.d_dim >= tbox::Dimension(3)) {
      co << " nz=" << sf.d_npi[2] << " pz=" << sf.d_npi[2];
   }
   co << " }";
   return co;
}
