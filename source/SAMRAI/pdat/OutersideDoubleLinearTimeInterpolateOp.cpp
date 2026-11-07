/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for double outerside patch data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/OutersideDoubleLinearTimeInterpolateOp.h"

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
#pragma warning(disable : 1419)
#endif

// in lintimint1d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidedoub1d, LINTIMEINTOUTSIDEDOUB1D)(const int &,
                                                                       const int &,
                                                                       const int &, const int &,
                                                                       const int &, const int &,
                                                                       const int &, const int &,
                                                                       const double &,
                                                                       const double *, const double *,
                                                                       double *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidedoub2d0, LINTIMEINTOUTSIDEDOUB2D0)(const int &,
                                                                         const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const double &,
                                                                         const double *, const double *,
                                                                         double *);
void SAMRAI_F77_FUNC(lintimeintoutsidedoub2d1, LINTIMEINTOUTSIDEDOUB2D1)(const int &,
                                                                         const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const double &,
                                                                         const double *, const double *,
                                                                         double *);
// in lintimint3d.f:
void SAMRAI_F77_FUNC(lintimeintoutsidedoub3d0, LINTIMEINTOUTSIDEDOUB3D0)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const double &,
                                                                         const double *, const double *,
                                                                         double *);
void SAMRAI_F77_FUNC(lintimeintoutsidedoub3d1, LINTIMEINTOUTSIDEDOUB3D1)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const double &,
                                                                         const double *, const double *,
                                                                         double *);
void SAMRAI_F77_FUNC(lintimeintoutsidedoub3d2, LINTIMEINTOUTSIDEDOUB3D2)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const double &,
                                                                         const double *, const double *,
                                                                         double *);
}

namespace SAMRAI
{
namespace pdat
{

OutersideDoubleLinearTimeInterpolateOp::OutersideDoubleLinearTimeInterpolateOp() : hier::TimeInterpolateOperator()
{
}

OutersideDoubleLinearTimeInterpolateOp::~OutersideDoubleLinearTimeInterpolateOp()
{
}

void
OutersideDoubleLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   NULL_USE(overlap);
   const tbox::Dimension& dim(where.getDim());

   const OutersideData<double> *old_dat =
       CPP_CAST<const OutersideData<double> *>(&src_data_old);
   const OutersideData<double> *new_dat =
       CPP_CAST<const OutersideData<double> *>(&src_data_new);
   OutersideData<double> *dst_dat =
       CPP_CAST<OutersideData<double> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT((where * old_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * new_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * dst_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

   const hier::Index &old_ilo = old_dat->getGhostBox().lower();
   const hier::Index &old_ihi = old_dat->getGhostBox().upper();
   const hier::Index &new_ilo = new_dat->getGhostBox().lower();
   const hier::Index &new_ihi = new_dat->getGhostBox().upper();

   const hier::Index &dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index &dst_ihi = dst_dat->getGhostBox().upper();

   const hier::Index &ifirst = where.lower();
   const hier::Index &ilast = where.upper();

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
      for (int side = 0; side < 2; ++side) {
         if (dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(lintimeintoutsidedoub1d,
                            LINTIMEINTOUTSIDEDOUB1D)
            (ifirst(0), ilast(0),
             old_ilo(0), old_ihi(0),
             new_ilo(0), new_ihi(0),
             dst_ilo(0), dst_ihi(0),
             tfrac,
             old_dat->getPointer(0, side, d),
             new_dat->getPointer(0, side, d),
             dst_dat->getPointer(0, side, d));
         } else if (dim == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
            {
               SAMRAI::hier::Box d0_box = where;
               if (side == 0) {
                  d0_box.setLower(0, where.lower(0));
                  d0_box.setUpper(0, where.lower(0));
               } else if (side == 1) {
                  d0_box.setLower(0, where.upper(0));
                  d0_box.setUpper(0, where.upper(0));
               }
               auto old_array = old_dat->getConstView<2>(0, side, d);
               auto new_array = new_dat->getConstView<2>(0, side, d);
               auto dst_array = dst_dat->getView<2>(0, side, d);

               hier::parallel_for_all(d0_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            }
            {
               SAMRAI::hier::Box d1_box = where;
               if (side == 0) {
                  d1_box.setLower(1, where.lower(1));
                  d1_box.setUpper(1, where.lower(1));
               } else if (side == 1) {
                  d1_box.setLower(1, where.upper(1));
                  d1_box.setUpper(1, where.upper(1));
               }
               auto old_array = old_dat->getConstView<2>(1, side, d);
               auto new_array = new_dat->getConstView<2>(1, side, d);
               auto dst_array = dst_dat->getView<2>(1, side, d);

               hier::parallel_for_all(d1_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            }
#else
            SAMRAI_F77_FUNC(lintimeintoutsidedoub2d0,
                            LINTIMEINTOUTSIDEDOUB2D0)
            (ifirst(0), ifirst(1), ilast(0),
             ilast(1),
             old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
             new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
             dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
             tfrac,
             old_dat->getPointer(0, side, d),
             new_dat->getPointer(0, side, d),
             dst_dat->getPointer(0, side, d));
            SAMRAI_F77_FUNC(lintimeintoutsidedoub2d1,
                            LINTIMEINTOUTSIDEDOUB2D1)
            (ifirst(0), ifirst(1), ilast(0),
             ilast(1),
             old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
             new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
             dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
             tfrac,
             old_dat->getPointer(1, side, d),
             new_dat->getPointer(1, side, d),
             dst_dat->getPointer(1, side, d));
#endif  // test for RAJA
         } else if (dim == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
            {
               SAMRAI::hier::Box d0_box = where;
               if (side == 0) {
                  d0_box.setLower(0, where.lower(0));
                  d0_box.setUpper(0, where.lower(0));
               } else if (side == 1) {
                  d0_box.setLower(0, where.upper(0));
                  d0_box.setUpper(0, where.upper(0));
               }
               auto old_array = old_dat->getConstView<3>(0, side, d);
               auto new_array = new_dat->getConstView<3>(0, side, d);
               auto dst_array = dst_dat->getView<3>(0, side, d);

               hier::parallel_for_all(d0_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
            {
               SAMRAI::hier::Box d1_box = where;
               if (side == 0) {
                  d1_box.setLower(1, where.lower(1));
                  d1_box.setUpper(1, where.lower(1));
               } else if (side == 1) {
                  d1_box.setLower(1, where.upper(1));
                  d1_box.setUpper(1, where.upper(1));
               }

               auto old_array = old_dat->getConstView<3>(1, side, d);
               auto new_array = new_dat->getConstView<3>(1, side, d);
               auto dst_array = dst_dat->getView<3>(1, side, d);

               hier::parallel_for_all(d1_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
            {
               SAMRAI::hier::Box d2_box = where;
               if (side == 0) {
                  d2_box.setLower(2, where.lower(2));
                  d2_box.setUpper(2, where.lower(2));
               } else if (side == 1) {
                  d2_box.setLower(2, where.upper(2));
                  d2_box.setUpper(2, where.upper(2));
               }

               auto old_array = old_dat->getConstView<3>(2, side, d);
               auto new_array = new_dat->getConstView<3>(2, side, d);
               auto dst_array = dst_dat->getView<3>(2, side, d);

               hier::parallel_for_all(d2_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
#else
            SAMRAI_F77_FUNC(lintimeintoutsidedoub3d0,
                            LINTIMEINTOUTSIDEDOUB3D0)
            (ifirst(0), ifirst(1), ifirst(2),
             ilast(0), ilast(1), ilast(2),
             old_ilo(0), old_ilo(1), old_ilo(2),
             old_ihi(0), old_ihi(1), old_ihi(2),
             new_ilo(0), new_ilo(1), new_ilo(2),
             new_ihi(0), new_ihi(1), new_ihi(2),
             dst_ilo(0), dst_ilo(1), dst_ilo(2),
             dst_ihi(0), dst_ihi(1), dst_ihi(2),
             tfrac,
             old_dat->getPointer(0, side, d),
             new_dat->getPointer(0, side, d),
             dst_dat->getPointer(0, side, d));
            SAMRAI_F77_FUNC(lintimeintoutsidedoub3d1,
                            LINTIMEINTOUTSIDEDOUB3D1)
            (ifirst(0), ifirst(1), ifirst(2),
             ilast(0), ilast(1), ilast(2),
             old_ilo(0), old_ilo(1), old_ilo(2),
             old_ihi(0), old_ihi(1), old_ihi(2),
             new_ilo(0), new_ilo(1), new_ilo(2),
             new_ihi(0), new_ihi(1), new_ihi(2),
             dst_ilo(0), dst_ilo(1), dst_ilo(2),
             dst_ihi(0), dst_ihi(1), dst_ihi(2),
             tfrac,
             old_dat->getPointer(1, side, d),
             new_dat->getPointer(1, side, d),
             dst_dat->getPointer(1, side, d));
            SAMRAI_F77_FUNC(lintimeintoutsidedoub3d2,
                            LINTIMEINTOUTSIDEDOUB3D2)
            (ifirst(0), ifirst(1), ifirst(2),
             ilast(0), ilast(1), ilast(2),
             old_ilo(0), old_ilo(1), old_ilo(2),
             old_ihi(0), old_ihi(1), old_ihi(2),
             new_ilo(0), new_ilo(1), new_ilo(2),
             new_ihi(0), new_ihi(1), new_ihi(2),
             dst_ilo(0), dst_ilo(1), dst_ilo(2),
             dst_ihi(0), dst_ihi(1), dst_ihi(2),
             tfrac,
             old_dat->getPointer(2, side, d),
             new_dat->getPointer(2, side, d),
             dst_dat->getPointer(2, side, d));
#endif  // test for RAJA
         } else {
            TBOX_ERROR(
                "OutersideDoubleLienarTimeInterpolate::TimeInterpolate dim > 3 not supported"
                << std::endl);
         }
      }
   }
}

}  // namespace pdat
}  // namespace SAMRAI
