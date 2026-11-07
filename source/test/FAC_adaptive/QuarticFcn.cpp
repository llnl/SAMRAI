/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Quartic function functor.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "QuarticFcn.h"
//#include IOSTREAM_HEADER_FILE

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#else
#ifndef included_strstream
#define included_strstream
#include <strstream.h>
#endif
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


#include "SAMRAI/tbox/Utilities.h"

/*
 * Temporary fix for g++ lacking instantiations when --no-implicit-templates
 * is used (by SAMRAI)
 */
#define fill_n(p, n, v) { size_t _i; for (_i = 0; _i < n; ++_i) p[_i] = v; }
#define copy_n(s, n, d) { size_t _i; for (_i = 0; _i < n; ++_i) d[_i] = s[_i]; }

QuarticFcn::QuarticFcn(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   fill_n(d_coefs, NUMBER_OF_3D_COEF, 0.0)
}

QuarticFcn::QuarticFcn(
   const QuarticFcn& r):
   d_dim(r.d_dim)
{
   copy_n(r.d_coefs, NUMBER_OF_3D_COEF, d_coefs)
}

int QuarticFcn::setPolynomialCoef(
   ctype type,
   double c) {
   d_coefs[type] = c;
   return 0;
}

double QuarticFcn::getPolynomialCoef(
   ctype type) const {
   return d_coefs[type];
}

double QuarticFcn::operator () (
   double x) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(1));
   double rval;
   rval = d_coefs[co_c]
      + d_coefs[co_x] * x
      + d_coefs[co_xx] * x * x
      + d_coefs[co_xxx] * x * x * x
      + d_coefs[co_xxxx] * x * x * x * x;
   return rval;
}
double QuarticFcn::operator () (
   double x,
   double y) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(2));
   double rval;
   rval = d_coefs[co_c]
      + d_coefs[co_x] * x + d_coefs[co_y] * y
      + d_coefs[co_xy] * x * y + d_coefs[co_xx] * x * x + d_coefs[co_yy] * y
      * y
      + d_coefs[co_xxx] * x * x * x + d_coefs[co_xxy] * x * x * y
      + d_coefs[co_xyy] * x * y * y + d_coefs[co_yyy] * y * y * y
      + d_coefs[co_xxxy] * x * x * x * y + d_coefs[co_xxyy] * x * x * y * y
      + d_coefs[co_xyyy] * x * y * y * y + d_coefs[co_yyyy] * y * y * y * y
   ;
   return rval;
}
double QuarticFcn::operator () (
   double x,
   double y,
   double z) const {
   TBOX_ASSERT(d_dim == tbox::Dimension(3));
   double rval;
   rval = d_coefs[co_c]
      + d_coefs[co_x] * x + d_coefs[co_y] * y + d_coefs[co_z] * z
      + d_coefs[co_xx] * x * x + d_coefs[co_yy] * y * y + d_coefs[co_zz] * z
      * z
      + d_coefs[co_xy] * x * y + d_coefs[co_xz] * x * z + d_coefs[co_yz] * y
      * z
      + d_coefs[co_zzz] * z * z * z + d_coefs[co_xxz] * x * x * z
      + d_coefs[co_xzz] * x * z * z
      + d_coefs[co_yyz] * y * y * z + d_coefs[co_yzz] * x * z * z
      + d_coefs[co_xyz] * x * y * z
      + d_coefs[co_xxxz] * x * x * x * z + d_coefs[co_xxyz] * x * x * y * z
      + d_coefs[co_xxzz] * x * x * z * x + d_coefs[co_xyyz] * x * y * y * z
      + d_coefs[co_xyzz] * x * y * z * z + d_coefs[co_xzzz] * x * z * z * z
      + d_coefs[co_yyyz] * y * y * y * z + d_coefs[co_yyzz] * y * y * z * z
      + d_coefs[co_yzzz] * y * z * z * z + d_coefs[co_zzzz] * z * z * z * z
   ;
   return rval;
}

QuarticFcn QuarticFcn::operator + (
   const QuarticFcn& r) const {
   QuarticFcn v(d_dim);
   size_t i;
   for (i = 0; i < NUMBER_OF_3D_COEF; ++i) {
      v.d_coefs[i] = d_coefs[i] + r.d_coefs[i];
   }
   return v;
}

QuarticFcn QuarticFcn::operator - (
   const QuarticFcn& r) const {
   QuarticFcn v(d_dim);
   size_t i;
   for (i = 0; i < NUMBER_OF_3D_COEF; ++i) {
      v.d_coefs[i] = d_coefs[i] - r.d_coefs[i];
   }
   return v;
}

QuarticFcn& QuarticFcn::operator = (
   const QuarticFcn& r) {
   size_t i;
   for (i = 0; i < NUMBER_OF_3D_COEF; ++i) {
      d_coefs[i] = r.d_coefs[i];
   }
   return *this;
}

QuarticFcn& QuarticFcn::differentiateSelf(
   unsigned short int x
   ,
   unsigned short int y)
{
   TBOX_ASSERT(d_dim == tbox::Dimension(2));
   /*
    * Since differentiation commutes,
    * simply differentiate one direction at a time.
    */
   // Differentiate in x direction.
   for ( ; x > 0; --x) {
      d_coefs[co_c] = d_coefs[co_x];
      d_coefs[co_x] = 2 * d_coefs[co_xx];
      d_coefs[co_xx] = 3 * d_coefs[co_xxx];
      d_coefs[co_xxx] = 4 * d_coefs[co_xxxx];
      d_coefs[co_xxxx] = 0;
      d_coefs[co_y] = d_coefs[co_xy];
      d_coefs[co_xy] = 2 * d_coefs[co_xxy];
      d_coefs[co_yy] = d_coefs[co_xyy];
      d_coefs[co_xxy] = 3 * d_coefs[co_xxxy];
      d_coefs[co_xyy] = 2 * d_coefs[co_xxyy];
      d_coefs[co_yyy] = d_coefs[co_xyyy];
      d_coefs[co_xxxy] = 0;
      d_coefs[co_xxyy] = 0;
      d_coefs[co_xyyy] = 0;
      d_coefs[co_yyyy] = 0;
   }
   // Differentiate in y direction.
   for ( ; y > 0; --y) {
      d_coefs[co_c] = d_coefs[co_y];
      d_coefs[co_x] = d_coefs[co_xy];
      d_coefs[co_y] = 2 * d_coefs[co_yy];
      d_coefs[co_xx] = d_coefs[co_xxy];
      d_coefs[co_xy] = 2 * d_coefs[co_xyy];
      d_coefs[co_yy] = 3 * d_coefs[co_yyy];
      d_coefs[co_xxx] = d_coefs[co_xxxy];
      d_coefs[co_xxy] = 2 * d_coefs[co_xxyy];
      d_coefs[co_xyy] = 3 * d_coefs[co_xyyy];
      d_coefs[co_yyy] = 4 * d_coefs[co_yyyy];
      d_coefs[co_xxxy] = 0;
      d_coefs[co_xxyy] = 0;
      d_coefs[co_xyyy] = 0;
      d_coefs[co_yyyy] = 0;
   }
   return *this;
}

QuarticFcn& QuarticFcn::differentiateSelf(
   unsigned short int x
   ,
   unsigned short int y
   ,
   unsigned short int z)
{
   TBOX_ASSERT(d_dim == tbox::Dimension(3));
   /*
    * Since differentiation commutes,
    * simply differentiate one direction at a time.
    */
   // Differentiate in x direction.
   for ( ; x > 0; --x) {
      d_coefs[co_c] = d_coefs[co_x];
      d_coefs[co_x] = 2 * d_coefs[co_xx];
      d_coefs[co_xx] = 3 * d_coefs[co_xxx];
      d_coefs[co_xxx] = 4 * d_coefs[co_xxxx];
      d_coefs[co_xxxx] = 0;
      d_coefs[co_y] = d_coefs[co_xy];
      d_coefs[co_xy] = 2 * d_coefs[co_xxy];
      d_coefs[co_yy] = d_coefs[co_xyy];
      d_coefs[co_xxy] = 3 * d_coefs[co_xxxy];
      d_coefs[co_xyy] = 2 * d_coefs[co_xxyy];
      d_coefs[co_yyy] = d_coefs[co_xyyy];
      d_coefs[co_xxxy] = 0;
      d_coefs[co_xxyy] = 0;
      d_coefs[co_xyyy] = 0;
      d_coefs[co_yyyy] = 0;
      d_coefs[co_z] = d_coefs[co_xz];
      d_coefs[co_xz] = 2 * d_coefs[co_xxz];
      d_coefs[co_yz] = d_coefs[co_xyz];
      d_coefs[co_zz] = d_coefs[co_xzz];
      d_coefs[co_xxz] = 3 * d_coefs[co_xxxz];
      d_coefs[co_xzz] = 2 * d_coefs[co_xxzz];
      d_coefs[co_yyz] = d_coefs[co_xyyz];
      d_coefs[co_yzz] = d_coefs[co_xyzz];
      d_coefs[co_xyz] = 2 * d_coefs[co_xxyz];
      d_coefs[co_zzz] = d_coefs[co_xzzz];
      d_coefs[co_xxxz] = 0;
      d_coefs[co_xxyz] = 0;
      d_coefs[co_xxzz] = 0;
      d_coefs[co_xyyz] = 0;
      d_coefs[co_xyzz] = 0;
      d_coefs[co_xzzz] = 0;
      d_coefs[co_yyyz] = 0;
      d_coefs[co_yyzz] = 0;
      d_coefs[co_yzzz] = 0;
      d_coefs[co_zzzz] = 0;
   }
   // Differentiate in y direction.
   for ( ; y > 0; --y) {
      d_coefs[co_c] = d_coefs[co_y];
      d_coefs[co_x] = d_coefs[co_xy];
      d_coefs[co_y] = 2 * d_coefs[co_yy];
      d_coefs[co_xx] = d_coefs[co_xxy];
      d_coefs[co_xy] = 2 * d_coefs[co_xyy];
      d_coefs[co_yy] = 3 * d_coefs[co_yyy];
      d_coefs[co_xxx] = d_coefs[co_xxxy];
      d_coefs[co_xxy] = 2 * d_coefs[co_xxyy];
      d_coefs[co_xyy] = 3 * d_coefs[co_xyyy];
      d_coefs[co_yyy] = 4 * d_coefs[co_yyyy];
      d_coefs[co_xxxy] = 0;
      d_coefs[co_xxyy] = 0;
      d_coefs[co_xyyy] = 0;
      d_coefs[co_yyyy] = 0;
      d_coefs[co_z] = d_coefs[co_yz];
      d_coefs[co_xz] = d_coefs[co_xyz];
      d_coefs[co_yz] = 2 * d_coefs[co_yyz];
      d_coefs[co_zz] = d_coefs[co_yzz];
      d_coefs[co_xxz] = d_coefs[co_xxyz];
      d_coefs[co_xzz] = d_coefs[co_xyzz];
      d_coefs[co_yyz] = 3 * d_coefs[co_yyyz];
      d_coefs[co_yzz] = 2 * d_coefs[co_yyzz];
      d_coefs[co_xyz] = 2 * d_coefs[co_xyyz];
      d_coefs[co_zzz] = d_coefs[co_yzzz];
      d_coefs[co_xxxz] = 0;
      d_coefs[co_xxyz] = 0;
      d_coefs[co_xxzz] = 0;
      d_coefs[co_xyyz] = 0;
      d_coefs[co_xyzz] = 0;
      d_coefs[co_xzzz] = 0;
      d_coefs[co_yyyz] = 0;
      d_coefs[co_yyzz] = 0;
      d_coefs[co_yzzz] = 0;
      d_coefs[co_zzzz] = 0;
   }
   // Differentiate in z direction.
   for ( ; z > 0; --z) {
      d_coefs[co_c] = d_coefs[co_z];
      d_coefs[co_x] = d_coefs[co_xz];
      d_coefs[co_y] = d_coefs[co_yz];
      d_coefs[co_z] = 2 * d_coefs[co_zz];
      d_coefs[co_xx] = d_coefs[co_xxz];
      d_coefs[co_xy] = d_coefs[co_xyz];
      d_coefs[co_yy] = d_coefs[co_yyz];
      d_coefs[co_xz] = 2 * d_coefs[co_xzz];
      d_coefs[co_yz] = 2 * d_coefs[co_yzz];
      d_coefs[co_zz] = 3 * d_coefs[co_zzz];
      d_coefs[co_xxx] = d_coefs[co_xxxz];
      d_coefs[co_xxy] = d_coefs[co_xxyz];
      d_coefs[co_xyy] = d_coefs[co_xyyz];
      d_coefs[co_yyy] = d_coefs[co_yyyz];
      d_coefs[co_xxz] = 2 * d_coefs[co_xxzz];
      d_coefs[co_xzz] = 3 * d_coefs[co_xzzz];
      d_coefs[co_yyz] = 2 * d_coefs[co_yyzz];
      d_coefs[co_yzz] = 3 * d_coefs[co_yzzz];
      d_coefs[co_xyz] = 2 * d_coefs[co_xyzz];
      d_coefs[co_zzz] = 4 * d_coefs[co_zzzz];
      d_coefs[co_xxxz] = 0;
      d_coefs[co_xxyz] = 0;
      d_coefs[co_xxzz] = 0;
      d_coefs[co_xyyz] = 0;
      d_coefs[co_xyzz] = 0;
      d_coefs[co_xzzz] = 0;
      d_coefs[co_yyyz] = 0;
      d_coefs[co_yyzz] = 0;
      d_coefs[co_yzzz] = 0;
      d_coefs[co_zzzz] = 0;
   }
   return *this;
}

QuarticFcn QuarticFcn::differentiate(
   unsigned short int x
   ,
   unsigned short int y) const
{
   TBOX_ASSERT(d_dim == tbox::Dimension(2));
   QuarticFcn rval(*this);
   rval.differentiateSelf(x, y);
   return rval;
}

QuarticFcn QuarticFcn::differentiate(
   unsigned short int x
   ,
   unsigned short int y
   ,
   unsigned short int z) const
{
   TBOX_ASSERT(d_dim == tbox::Dimension(3));
   QuarticFcn rval(*this);
   rval.differentiateSelf(x, y, z);
   return rval;
}

std::ostream& operator << (
   std::ostream& co,
   const QuarticFcn& qf) {
   unsigned short int i;
   co << '{' << 0 << '=' << qf.d_coefs[0];
   switch (qf.d_dim.getValue()) {
      case 1:
         for (i = 1; i < NUMBER_OF_1D_COEF; ++i) {
            co << ' ' << i << '=' << qf.d_coefs[i];
         }
         break;
      case 2:
         for (i = 1; i < NUMBER_OF_2D_COEF; ++i) {
            co << ' ' << i << '=' << qf.d_coefs[i];
         }
         break;
      case 3:
         for (i = 1; i < NUMBER_OF_3D_COEF; ++i) {
            co << ' ' << i << '=' << qf.d_coefs[i];
         }
         break;
   }
   co << '}';
   return co;
}

#define EAT_WS(s) \
   { while (s.peek() == ' '                     \
            || s.peek() == '\t'                    \
            || s.peek() == '\n') { s.get(); } }

std::istream& operator >> (
   std::istream& ci,
   QuarticFcn& qf) {
   fill_n(qf.d_coefs, NUMBER_OF_3D_COEF, 0.0);
   int number_of_coefs = 0;
   switch (qf.d_dim.getValue()) {
      case 1:
         number_of_coefs = NUMBER_OF_1D_COEF;
         break;
      case 2:
         number_of_coefs = NUMBER_OF_2D_COEF;
         break;
      case 3:
         number_of_coefs = NUMBER_OF_3D_COEF;
         break;
   }

   NULL_USE(number_of_coefs);

   unsigned short int j;
   int dummy;
   EAT_WS(ci) // ci >> std::noskipws; // ci.ipfx(0);
   dummy = ci.get();
   TBOX_ASSERT(dummy == '{');
   while (ci >> j) {
      TBOX_ASSERT(j < number_of_coefs);
      ci >> dummy >> qf.d_coefs[j];
   }
   return ci;
}
