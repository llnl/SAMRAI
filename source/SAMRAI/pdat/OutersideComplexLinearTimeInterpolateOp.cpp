/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for complex outerside data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/OutersideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/tbox/Complex.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"


/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in lintimint1d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx1d, LINTIMEINTOUTSIDECMPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx2d0,
                     LINTIMEINTOUTSIDECMPLX2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx2d1,
                     LINTIMEINTOUTSIDECMPLX2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
// in lintimint3d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d0,
                     LINTIMEINTOUTSIDECMPLX3D0) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d1,
                     LINTIMEINTOUTSIDECMPLX3D1) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
void SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d2,
                     LINTIMEINTOUTSIDECMPLX3D2) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
}

namespace SAMRAI {
namespace pdat {

OutersideComplexLinearTimeInterpolateOp::
OutersideComplexLinearTimeInterpolateOp():
   hier::TimeInterpolateOperator()
{
}

OutersideComplexLinearTimeInterpolateOp::~
OutersideComplexLinearTimeInterpolateOp()
{
}

void
OutersideComplexLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   NULL_USE(overlap);
   const tbox::Dimension& dim(where.getDim());

   const OutersideData<dcomplex>* old_dat =
      CPP_CAST<const OutersideData<dcomplex> *>(&src_data_old);
   const OutersideData<dcomplex>* new_dat =
      CPP_CAST<const OutersideData<dcomplex> *>(&src_data_new);
   OutersideData<dcomplex>* dst_dat =
      CPP_CAST<OutersideData<dcomplex> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT((where * old_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * new_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * dst_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

   const hier::Index& old_ilo = old_dat->getGhostBox().lower();
   const hier::Index& old_ihi = old_dat->getGhostBox().upper();
   const hier::Index& new_ilo = new_dat->getGhostBox().lower();
   const hier::Index& new_ihi = new_dat->getGhostBox().upper();

   const hier::Index& dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index& dst_ihi = dst_dat->getGhostBox().upper();

   const hier::Index& ifirst = where.lower();
   const hier::Index& ilast = where.upper();

   const double old_time = old_dat->getTime();
   const double new_time = new_dat->getTime();
   const double dst_time = dst_dat->getTime();

   TBOX_ASSERT((old_time < dst_time ||
                tbox::MathUtilities<double>::equalEps(old_time, dst_time)) &&
      (dst_time < new_time ||
       tbox::MathUtilities<double>::equalEps(dst_time, new_time)));

   double tfrac = dst_time - old_time;
   double denom = new_time - old_time;
   if (denom > tbox::MathUtilities<double>::getMin()) {
      tfrac /= denom;
   } else {
      tfrac = 0.0;
   }

   for (int d = 0; d < dst_dat->getDepth(); ++d) {
      // loop over lower and upper outerside arrays
      for (int i = 0; i < 2; ++i) {
         if (dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx1d,
               LINTIMEINTOUTSIDECMPLX1D) (ifirst(0), ilast(0),
               old_ilo(0), old_ihi(0),
               new_ilo(0), new_ihi(0),
               dst_ilo(0), dst_ihi(0),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
         } else if (dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx2d0,
               LINTIMEINTOUTSIDECMPLX2D0) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx2d1,
               LINTIMEINTOUTSIDECMPLX2D1) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(1, i, d),
               new_dat->getPointer(1, i, d),
               dst_dat->getPointer(1, i, d));
         } else if (dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d0,
               LINTIMEINTOUTSIDECMPLX3D0) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d1,
               LINTIMEINTOUTSIDECMPLX3D1) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(1, i, d),
               new_dat->getPointer(1, i, d),
               dst_dat->getPointer(1, i, d));
            SAMRAI_F77_FUNC(lintimeintoutsidecmplx3d2,
               LINTIMEINTOUTSIDECMPLX3D2) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(2, i, d),
               new_dat->getPointer(2, i, d),
               dst_dat->getPointer(2, i, d));
         } else {
            TBOX_ERROR(
               "OutersideComplexLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
               << std::endl);
         }
      }
   }
}

}
}
