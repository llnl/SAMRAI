/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for float outerside patch data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/OutersideFloatLinearTimeInterpolateOp.h"

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
void SAMRAI_F77_FUNC(lintimeintoutsidefloat1d, LINTIMEINTOUTSIDEFLOAT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidefloat2d0,
                     LINTIMEINTOUTSIDEFLOAT2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
void SAMRAI_F77_FUNC(lintimeintoutsidefloat2d1,
                     LINTIMEINTOUTSIDEFLOAT2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
// in lintimint3d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidefloat3d0,
                     LINTIMEINTOUTSIDEFLOAT3D0) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
void SAMRAI_F77_FUNC(lintimeintoutsidefloat3d1,
                     LINTIMEINTOUTSIDEFLOAT3D1) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
void SAMRAI_F77_FUNC(lintimeintoutsidefloat3d2,
                     LINTIMEINTOUTSIDEFLOAT3D2) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
}

namespace SAMRAI {
namespace pdat {

OutersideFloatLinearTimeInterpolateOp::OutersideFloatLinearTimeInterpolateOp():
   hier::TimeInterpolateOperator()
{
}

OutersideFloatLinearTimeInterpolateOp::~OutersideFloatLinearTimeInterpolateOp()
{
}

void
OutersideFloatLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   NULL_USE(overlap);
   const tbox::Dimension& dim(where.getDim());

   const OutersideData<float>* old_dat =
      CPP_CAST<const OutersideData<float> *>(&src_data_old);
   const OutersideData<float>* new_dat =
      CPP_CAST<const OutersideData<float> *>(&src_data_new);
   OutersideData<float>* dst_dat =
      CPP_CAST<OutersideData<float> *>(&dst_data);

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
            SAMRAI_F77_FUNC(lintimeintoutsidefloat1d,
               LINTIMEINTOUTSIDEFLOAT1D) (ifirst(0), ilast(0),
               old_ilo(0), old_ihi(0),
               new_ilo(0), new_ihi(0),
               dst_ilo(0), dst_ihi(0),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
         } else if (dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(lintimeintoutsidefloat2d0,
               LINTIMEINTOUTSIDEFLOAT2D0) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
            SAMRAI_F77_FUNC(lintimeintoutsidefloat2d1,
               LINTIMEINTOUTSIDEFLOAT2D1) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(1, i, d),
               new_dat->getPointer(1, i, d),
               dst_dat->getPointer(1, i, d));
         } else if (dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(lintimeintoutsidefloat3d0,
               LINTIMEINTOUTSIDEFLOAT3D0) (ifirst(0), ifirst(1), ifirst(2),
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
            SAMRAI_F77_FUNC(lintimeintoutsidefloat3d1,
               LINTIMEINTOUTSIDEFLOAT3D1) (ifirst(0), ifirst(1), ifirst(2),
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
            SAMRAI_F77_FUNC(lintimeintoutsidefloat3d2,
               LINTIMEINTOUTSIDEFLOAT3D2) (ifirst(0), ifirst(1), ifirst(2),
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
               "OutersideFloatLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
               << std::endl);
         }
      }
   }
}

}
}
