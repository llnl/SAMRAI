/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for edge-centered double patch data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"
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
void SAMRAI_F77_FUNC(lintimeintedgedoub1d, LINTIMEINTEDGEDOUB1D)(const int &,
                                                                 const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const double &,
                                                                 const double *, const double *,
                                                                 double *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintedgedoub2d0, LINTIMEINTEDGEDOUB2D0)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintedgedoub2d1, LINTIMEINTEDGEDOUB2D1)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintedgedoub3d0, LINTIMEINTEDGEDOUB3D0)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintedgedoub3d1, LINTIMEINTEDGEDOUB3D1)(const int &,
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
void SAMRAI_F77_FUNC(lintimeintedgedoub3d2, LINTIMEINTEDGEDOUB3D2)(const int &,
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

EdgeDoubleLinearTimeInterpolateOp::EdgeDoubleLinearTimeInterpolateOp() : hier::TimeInterpolateOperator()
{
}

EdgeDoubleLinearTimeInterpolateOp::~EdgeDoubleLinearTimeInterpolateOp()
{
}

void
EdgeDoubleLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   const tbox::Dimension &dim(where.getDim());

   const EdgeData<double> *old_dat =
       CPP_CAST<const EdgeData<double> *>(&src_data_old);
   const EdgeData<double> *new_dat =
       CPP_CAST<const EdgeData<double> *>(&src_data_new);
   EdgeData<double> *dst_dat =
       CPP_CAST<EdgeData<double> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

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

   std::vector<hier::Box> edge_where(dim.getValue(), hier::Box(dim));
   std::vector<hier::BoxContainer> ovlp_boxes(dim.getValue());

   const EdgeOverlap* edge_ovlp = CPP_CAST<const EdgeOverlap*>(&overlap);
   for (int dir = 0; dir < dim.getValue(); ++dir) {
      edge_where[dir] = EdgeGeometry::toEdgeBox(where, dir);
      edge_ovlp->getSourceBoxContainer(ovlp_boxes[dir], dir);
   }

   for (int d = 0; d < dst_dat->getDepth(); ++d) {
      if (dim == tbox::Dimension(1)) {
         for (auto itr = ovlp_boxes[0].begin();
              itr != ovlp_boxes[0].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[0]); 
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub1d, LINTIMEINTEDGEDOUB1D) (
               ifirst(0), ilast(0),
               old_ilo(0), old_ihi(0),
               new_ilo(0), new_ihi(0),
               dst_ilo(0), dst_ihi(0),
               tfrac,
               old_dat->getPointer(0, d),
               new_dat->getPointer(0, d),
               dst_dat->getPointer(0, d));
         }
      } else if (dim == tbox::Dimension(2)) {
         for (auto itr = ovlp_boxes[0].begin();
              itr != ovlp_boxes[0].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[0]);
#if defined(HAVE_RAJA)
            {
               auto old_array = old_dat->getConstView<2>(0, d);
               auto new_array = new_dat->getConstView<2>(0, d);
               auto dst_array = dst_dat->getView<2>(0, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            }
#else
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub2d0, LINTIMEINTEDGEDOUB2D0) (
               ifirst(0), ifirst(1), ilast(0), ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(0, d),
               new_dat->getPointer(0, d),
               dst_dat->getPointer(0, d));
#endif // HAVE_RAJA for 2D, 0 direction
         } // end iterate ovrer ovlp_boxes[0]

         for (auto itr = ovlp_boxes[1].begin();
              itr != ovlp_boxes[1].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[1]);

#if defined(HAVE_RAJA)
            {
               auto old_array = old_dat->getConstView<2>(1, d);
               auto new_array = new_dat->getConstView<2>(1, d);
               auto dst_array = dst_dat->getView<2>(1, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(j, k) = old_array(j, k) * oldfrac + new_array(j, k) * tfrac;
               });
            }
#else
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub2d1, LINTIMEINTEDGEDOUB2D1) (
               ifirst(0), ifirst(1), ilast(0), ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(1, d),
               new_dat->getPointer(1, d),
               dst_dat->getPointer(1, d));
#endif // HAVE_RAJA for 2D, 1 direction
         } // end iterate over ovlp_boxes[1]
      } else if (dim == tbox::Dimension(3)) {
         for (auto itr = ovlp_boxes[0].begin();
              itr != ovlp_boxes[0].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[0]);
#if defined(HAVE_RAJA)
            {
               auto old_array = old_dat->getConstView<3>(0, d);
               auto new_array = new_dat->getConstView<3>(0, d);
               auto dst_array = dst_dat->getView<3>(0, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
#else
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub3d0, LINTIMEINTEDGEDOUB3D0) (
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
         
#endif // HAVE RAJA for 3D, 0 direction
         } // end iterate ovlp_boxes[0]

         for (auto itr = ovlp_boxes[1].begin();
              itr != ovlp_boxes[1].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[1]);
#if defined(HAVE_RAJA)
            {
               auto old_array = old_dat->getConstView<3>(1, d);
               auto new_array = new_dat->getConstView<3>(1, d);
               auto dst_array = dst_dat->getView<3>(1, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
#else
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub3d1, LINTIMEINTEDGEDOUB3D1) (
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

#endif // HAVE RAJA for 3D, 1 direction
         } // end iterate ovlp_boxes[1]

         for (auto itr = ovlp_boxes[2].begin();
              itr != ovlp_boxes[2].end(); ++itr) {
            hier::Box dest_box((*itr) * edge_where[2]);
#if defined(HAVE_RAJA)
            {
               auto old_array = old_dat->getConstView<3>(2, d);
               auto new_array = new_dat->getConstView<3>(2, d);
               auto dst_array = dst_dat->getView<3>(2, d);

               hier::parallel_for_all(dest_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
                  const double oldfrac = 1.0 - tfrac;
                  dst_array(i, j, k) = old_array(i, j, k) * oldfrac + new_array(i, j, k) * tfrac;
               });
            }
#else
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintedgedoub3d2, LINTIMEINTEDGEDOUB3D2) (
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
#endif // HAVE RAJA for 3D, 2 direction
         } // end iterate ovlp_boxes[2]
      } else {
         TBOX_ERROR(
             "EdgeDoubleLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
             << std::endl);
      }
   }
}

}  // namespace pdat
}  // namespace SAMRAI
