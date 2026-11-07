/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for side-centered double patch data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
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
void SAMRAI_F77_FUNC(lintimeintsidedoub1d, LINTIMEINTSIDEDOUB1D)(const int &,
                                                                 const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const double &,
                                                                 const double *, const double *,
                                                                 double *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintsidedoub2d0, LINTIMEINTSIDEDOUB2D0)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintsidedoub2d1, LINTIMEINTSIDEDOUB2D1)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintsidedoub3d0, LINTIMEINTSIDEDOUB3D0)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintsidedoub3d1, LINTIMEINTSIDEDOUB3D1)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintsidedoub3d2, LINTIMEINTSIDEDOUB3D2)(const int &,
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

SideDoubleLinearTimeInterpolateOp::SideDoubleLinearTimeInterpolateOp() : hier::TimeInterpolateOperator()
{
}

SideDoubleLinearTimeInterpolateOp::~SideDoubleLinearTimeInterpolateOp()
{
}

void
SideDoubleLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   const tbox::Dimension &dim(where.getDim());

   const SideData<double> *old_dat =
       CPP_CAST<const SideData<double> *>(&src_data_old);
   const SideData<double> *new_dat =
       CPP_CAST<const SideData<double> *>(&src_data_new);
   SideData<double> *dst_dat =
       CPP_CAST<SideData<double> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

   const hier::IntVector &directions = dst_dat->getDirectionVector();

   TBOX_ASSERT(directions ==
               hier::IntVector::min(directions, old_dat->getDirectionVector()));
   TBOX_ASSERT(directions ==
               hier::IntVector::min(directions, new_dat->getDirectionVector()));

   const hier::Index &old_ilo = old_dat->getGhostBox().lower();
   const hier::Index &old_ihi = old_dat->getGhostBox().upper();
   const hier::Index &new_ilo = new_dat->getGhostBox().lower();
   const hier::Index &new_ihi = new_dat->getGhostBox().upper();

   const hier::Index &dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index &dst_ihi = dst_dat->getGhostBox().upper();

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

   std::vector<hier::Box> side_where(dim.getValue(), hier::Box(dim));
   std::vector<hier::BoxContainer> ovlp_boxes(dim.getValue());

   const SideOverlap* side_ovlp = CPP_CAST<const SideOverlap*>(&overlap);
   for (int dir = 0; dir < dim.getValue(); ++dir) {
      if (directions(dir)) {
         side_where[dir] = SideGeometry::toSideBox(where, dir);
         side_ovlp->getSourceBoxContainer(ovlp_boxes[dir], dir);
      }
   }

   for (int d = 0; d < dst_dat->getDepth(); ++d) {
      if (dim == tbox::Dimension(1)) {
         if (directions(0)) {
            for (auto itr = ovlp_boxes[0].begin();
                 itr != ovlp_boxes[0].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[0]); 
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub1d, LINTIMEINTSIDEDOUB1D) (
                  ifirst(0), ilast(0),
                  old_ilo(0), old_ihi(0),
                  new_ilo(0), new_ihi(0),
                  dst_ilo(0), dst_ihi(0),
                  tfrac,
                  old_dat->getPointer(0, d),
                  new_dat->getPointer(0, d),
                  dst_dat->getPointer(0, d));
            }
         }
      } else if (dim == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
         if (directions(0)) {
            for (auto itr = ovlp_boxes[0].begin();
                 itr != ovlp_boxes[0].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[0]);
               auto old_array = old_dat->getConstView<2>(0, d);
               auto new_array = new_dat->getConstView<2>(0, d);
               auto dst_array = dst_dat->getView<2>(0, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            } // end iterate ovlp_boxes[0]
         } // directions(0)
         if (directions(1)) {
            for (auto itr = ovlp_boxes[1].begin();
                 itr != ovlp_boxes[1].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[1]);
               auto old_array = old_dat->getConstView<2>(1, d);
               auto new_array = new_dat->getConstView<2>(1, d);
               auto dst_array = dst_dat->getView<2>(1, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            } // end iterate ovlp_boxes[1]
         } // directions(1)
#else
         if (directions(0)) {
            for (auto itr = ovlp_boxes[0].begin();
                 itr != ovlp_boxes[0].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[0]);
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub2d0, LINTIMEINTSIDEDOUB2D0) (
                  ifirst(0), ifirst(1), ilast(0), ilast(1),
                  old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
                  new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
                  dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
                  tfrac,
                  old_dat->getPointer(0, d),
                  new_dat->getPointer(0, d),
                  dst_dat->getPointer(0, d));
            }
         }
         if (directions(1)) {
            for (auto itr = ovlp_boxes[1].begin();
                 itr != ovlp_boxes[1].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[1]);
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub2d1, LINTIMEINTSIDEDOUB2D1) (
                  ifirst(0), ifirst(1), ilast(0), ilast(1),
                  old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
                  new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
                  dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
                  tfrac,
                  old_dat->getPointer(1, d),
                  new_dat->getPointer(1, d),
                  dst_dat->getPointer(1, d));
            }
         }
#endif  // test for RAJA
      } else if (dim == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
         if (directions(0)) {
            for (auto itr = ovlp_boxes[0].begin();
                 itr != ovlp_boxes[0].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[0]);
               auto old_array = old_dat->getConstView<3>(0, d);
               auto new_array = new_dat->getConstView<3>(0, d);
               auto dst_array = dst_dat->getView<3>(0, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            } // end iterate ovlp_boxes[0]
         } // directions(0)
         if (directions(1)) {
            for (auto itr = ovlp_boxes[1].begin();
                 itr != ovlp_boxes[1].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[1]);
               auto old_array = old_dat->getConstView<3>(1, d);
               auto new_array = new_dat->getConstView<3>(1, d);
               auto dst_array = dst_dat->getView<3>(1, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            } // end iterate ovlp_boxes[1]
         } // directions(1)
         if (directions(2)) {
            for (auto itr = ovlp_boxes[2].begin();
                 itr != ovlp_boxes[2].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[2]);
               auto old_array = old_dat->getConstView<3>(2, d);
               auto new_array = new_dat->getConstView<3>(2, d);
               auto dst_array = dst_dat->getView<3>(2, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            } // end iterate over ovlp_boxes[2]
         } // directions(2)
#else
         if (directions(0)) {
            for (auto itr = ovlp_boxes[0].begin();
                 itr != ovlp_boxes[0].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[0]);
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub3d0, LINTIMEINTSIDEDOUB3D0) (
                  ifirst(0), ifirst(1), ifirst(2),
                  ilast(0), ilast(1), ilast(2),
                  old_ilo(0), old_ilo(1), old_ilo(2),
                  old_ihi(0), old_ihi(1), old_ihi(2),
                  new_ilo(0), new_ilo(1), new_ilo(2),
                  new_ihi(0), new_ihi(1), new_ihi(2),
                  dst_ilo(0), dst_ilo(1), dst_ilo(2),
                  dst_ihi(0), dst_ihi(1), dst_ihi(2),
                  tfrac,
                  old_dat->getPointer(0, d),
                  new_dat->getPointer(0, d),
                  dst_dat->getPointer(0, d));
            }
         }
         if (directions(1)) {
            for (auto itr = ovlp_boxes[1].begin();
                 itr != ovlp_boxes[1].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[1]);
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub3d1, LINTIMEINTSIDEDOUB3D1) (
                  ifirst(0), ifirst(1), ifirst(2),
                  ilast(0), ilast(1), ilast(2),
                  old_ilo(0), old_ilo(1), old_ilo(2),
                  old_ihi(0), old_ihi(1), old_ihi(2),
                  new_ilo(0), new_ilo(1), new_ilo(2),
                  new_ihi(0), new_ihi(1), new_ihi(2),
                  dst_ilo(0), dst_ilo(1), dst_ilo(2),
                  dst_ihi(0), dst_ihi(1), dst_ihi(2),
                  tfrac,
                  old_dat->getPointer(1, d),
                  new_dat->getPointer(1, d),
                  dst_dat->getPointer(1, d));
            }
         }
         if (directions(2)) {
            for (auto itr = ovlp_boxes[2].begin();
                 itr != ovlp_boxes[2].end(); ++itr) {
               hier::Box dest_box((*itr) * side_where[2]);
               const hier::Index& ifirst = dest_box.lower();
               const hier::Index& ilast = dest_box.upper();
               SAMRAI_F77_FUNC(lintimeintsidedoub3d2, LINTIMEINTSIDEDOUB3D2) (
                  ifirst(0), ifirst(1), ifirst(2),
                  ilast(0), ilast(1), ilast(2),
                  old_ilo(0), old_ilo(1), old_ilo(2),
                  old_ihi(0), old_ihi(1), old_ihi(2),
                  new_ilo(0), new_ilo(1), new_ilo(2),
                  new_ihi(0), new_ihi(1), new_ihi(2),
                  dst_ilo(0), dst_ilo(1), dst_ilo(2),
                  dst_ihi(0), dst_ihi(1), dst_ihi(2),
                  tfrac,
                  old_dat->getPointer(2, d),
                  new_dat->getPointer(2, d),
                  dst_dat->getPointer(2, d));
            }
         }
#endif  // test for RAJA
      } else {
         TBOX_ERROR(
             "SideDoubleLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
             << std::endl);
      }
   }
}

}  // namespace pdat
}  // namespace SAMRAI
