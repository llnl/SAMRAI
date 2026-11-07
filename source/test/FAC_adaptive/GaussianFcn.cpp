/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Gaussian function support for FAC solver tests.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "GaussianFcn.h"
#include <math.h>
#include <stdlib.h>
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <string>
#include <string.h>

/*
 * Temporary fix for g++ lacking instantiations when --no-implicit-templates
 * is used (by SAMRAI)
 */
#define fill_n(p, n, v) { size_t _i; for (_i = 0; _i < n; ++_i) p[_i] = v; }
#define copy_n(s, n, d) { size_t _i; for (_i = 0; _i < n; ++_i) d[_i] = s[_i]; }

GaussianFcn::GaussianFcn():
   d_dim(2),
   d_amp(1.0),
   d_lambda(-1.0)
{
   fill_n(d_center, d_dim.getValue(), 0.0)
}

GaussianFcn::GaussianFcn(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_amp(1.0),
   d_lambda(-1.0)
{
   fill_n(d_center, d_dim.getValue(), 0.0)
}

GaussianFcn::GaussianFcn(
   const GaussianFcn& other):
   d_dim(other.d_dim),
   d_amp(other.d_amp),
   d_lambda(other.d_lambda)
{
   copy_n(other.d_center, d_dim.getValue(), d_center)
}

int GaussianFcn::setAmplitude(
   const double amp) {
   d_amp = amp;
   return 0;
}

int GaussianFcn::setLambda(
   const double lambda) {
   d_lambda = lambda;
   return 0;
}

int GaussianFcn::setCenter(
   const double* center) {
   for (size_t i = 0; i < d_dim.getValue(); ++i) d_center[i] = center[i];
   return 0;
}

double GaussianFcn::getAmplitude() const {
   return d_amp;
}

double GaussianFcn::getLambda() const {
   return d_lambda;
}

int GaussianFcn::getCenter(
   double* center) const {
   for (size_t i = 0; i < d_dim.getValue(); ++i) center[i] = d_center[i];
   return 0;
}

double GaussianFcn::operator () (
   double x) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(1));
   double rval;
   rval = (x - d_center[0]) * (x - d_center[0]);
   rval = exp(d_lambda * rval);
   return rval;
}
double GaussianFcn::operator () (
   double x,
   double y) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(2));
   double rval;
   rval =
      (x
       - d_center[0])
      * (x - d_center[0]) + (y - d_center[1]) * (y - d_center[1]);
   rval = exp(d_lambda * rval);
   return rval;
}
double GaussianFcn::operator () (
   double x,
   double y,
   double z) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(3));
   double rval;
   rval =
      (x - d_center[0]) * (x - d_center[0])
      + (y - d_center[1]) * (y - d_center[1])
      + (z - d_center[2]) * (z - d_center[2]);
   rval = exp(d_lambda * rval);
   return rval;
}

GaussianFcn& GaussianFcn::operator = (
   const GaussianFcn& r) {
   TBOX_ASSERT(d_dim == r.d_dim);
   d_amp = r.d_amp;
   d_lambda = r.d_lambda;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_center[i] = r.d_center[i];
   }
   return *this;
}

#define EAT_WS(s) \
   { while (s.peek() == ' '                     \
            || s.peek() == '\t'                    \
            || s.peek() == '\n') { s.get(); } }

std::istream& operator >> (
   std::istream& ci,
   GaussianFcn& gf) {
   fill_n(gf.d_center, gf.d_dim.getValue(), 0.0)
   gf.d_amp = 1.0;
   gf.d_lambda = -1.0;
   char dummy, name[6];
   EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
   ci >> dummy;
   TBOX_ASSERT(dummy == '{');
   EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
   while (ci.peek() != '}') {
      ci.read(name, 2);
      if (name[0] == 'l') {
         // Expect form lambda=<float>
         ci.read(name, 5);
         TBOX_ASSERT(!strncmp(name, "mbda=", 5));
         ci >> gf.d_lambda;
         EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
      } else if (name[0] == 'a') {
         // Expect form amp=<float>
         ci.read(name, 2);
         TBOX_ASSERT(!strncmp(name, "p=", 2));
         ci >> gf.d_amp;
         EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
      } else if (name[0] == 'c') {
         // Expect form c[xyz]=<float>
         int dim(name[1] == 'x' ? 0 :
                 name[1] == 'y' ? 1 :
                 name[1] == 'z' ? 2 : 3);
         TBOX_ASSERT(dim < gf.d_dim.getValue());
         ci >> dummy;
         TBOX_ASSERT(dummy == '=');
         ci >> gf.d_center[dim];
         EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
      } else {
         tbox::SAMRAI_MPI::abort();
      }
   }
   return ci;
}

std::ostream& operator << (
   std::ostream& co,
   const GaussianFcn& gf) {
   co << "{ amp=" << gf.d_amp << " lambda=" << gf.d_lambda
   << " cx=" << gf.d_center[0];
   if (gf.d_dim >= tbox::Dimension(2)) {
      co << " cy=" << gf.d_center[1];
   }
   if (gf.d_dim >= tbox::Dimension(3)) {
      co << " cz=" << gf.d_center[2];
   }
   co << " }";
   return co;
}
