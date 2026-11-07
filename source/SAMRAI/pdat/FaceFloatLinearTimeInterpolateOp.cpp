/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for face-centered float patch data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/FaceFloatLinearTimeInterpolateOp.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
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
void SAMRAI_F77_FUNC(lintimeintfacefloat1d, LINTIMEINTFACEFLOAT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const float *, const float *,
   float *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintfacefloat2d0, LINTIMEINTFACEFLOAT2D0) (const int&,
   const int&,
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
void SAMRAI_F77_FUNC(lintimeintfacefloat2d1, LINTIMEINTFACEFLOAT2D1) (const int&,
   const int&,
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
void SAMRAI_F77_FUNC(lintimeintfacefloat3d0, LINTIMEINTFACEFLOAT3D0) (const int&,
   const int&, const int&,
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
void SAMRAI_F77_FUNC(lintimeintfacefloat3d1, LINTIMEINTFACEFLOAT3D1) (const int&,
   const int&, const int&,
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
void SAMRAI_F77_FUNC(lintimeintfacefloat3d2, LINTIMEINTFACEFLOAT3D2) (const int&,
   const int&, const int&,
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

FaceFloatLinearTimeInterpolateOp::FaceFloatLinearTimeInterpolateOp():
   hier::TimeInterpolateOperator()
{
}

FaceFloatLinearTimeInterpolateOp::~FaceFloatLinearTimeInterpolateOp()
{
}

void
FaceFloatLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   const tbox::Dimension& dim(where.getDim());

   const FaceData<float>* old_dat =
      CPP_CAST<const FaceData<float> *>(&src_data_old);
   const FaceData<float>* new_dat =
      CPP_CAST<const FaceData<float> *>(&src_data_new);
   FaceData<float>* dst_dat =
      CPP_CAST<FaceData<float> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

   const hier::Index& old_ilo = old_dat->getGhostBox().lower();
   const hier::Index& old_ihi = old_dat->getGhostBox().upper();
   const hier::Index& new_ilo = new_dat->getGhostBox().lower();
   const hier::Index& new_ihi = new_dat->getGhostBox().upper();

   const hier::Index& dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index& dst_ihi = dst_dat->getGhostBox().upper();

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

   std::vector<hier::Box> face_where(dim.getValue(), hier::Box(dim));
   std::vector<hier::BoxContainer> ovlp_boxes(dim.getValue());

   const FaceOverlap* face_ovlp = CPP_CAST<const FaceOverlap*>(&overlap);
   for (int dir = 0; dir < dim.getValue(); ++dir) {
      face_where[dir] = FaceGeometry::toFaceBox(where, dir);
      face_ovlp->getSourceBoxContainer(ovlp_boxes[dir], dir);
   }

   for (int d = 0; d < dst_dat->getDepth(); ++d) {
      if (dim == tbox::Dimension(1)) {
         for (auto itr = ovlp_boxes[0].begin();
              itr != ovlp_boxes[0].end(); ++itr) {
            hier::Box dest_box((*itr) * face_where[0]); 
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat1d, LINTIMEINTFACEFLOAT1D) (
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
            hier::Box dest_box((*itr) * face_where[0]);
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat2d0, LINTIMEINTFACEFLOAT2D0) (
               ifirst(0), ifirst(1), ilast(0), ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(0, d),
               new_dat->getPointer(0, d),
               dst_dat->getPointer(0, d));
         }
         for (auto itr = ovlp_boxes[1].begin();
              itr != ovlp_boxes[1].end(); ++itr) {
            hier::Box dest_box((*itr) * face_where[1]);
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat2d1, LINTIMEINTFACEFLOAT2D1) (
               ifirst(0), ifirst(1), ilast(0), ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(1, d),
               new_dat->getPointer(1, d),
               dst_dat->getPointer(1, d));
         }
      } else if (dim == tbox::Dimension(3)) {
         for (auto itr = ovlp_boxes[0].begin();
              itr != ovlp_boxes[0].end(); ++itr) {
            hier::Box dest_box((*itr) * face_where[0]);
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat3d0, LINTIMEINTFACEFLOAT3D0) (
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
         for (auto itr = ovlp_boxes[1].begin();
              itr != ovlp_boxes[1].end(); ++itr) {
            hier::Box dest_box((*itr) * face_where[1]);
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat3d1, LINTIMEINTFACEFLOAT3D1) (
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
         for (auto itr = ovlp_boxes[2].begin();
              itr != ovlp_boxes[2].end(); ++itr) {
            hier::Box dest_box((*itr) * face_where[2]);
            const hier::Index& ifirst = dest_box.lower();
            const hier::Index& ilast = dest_box.upper();
            SAMRAI_F77_FUNC(lintimeintfacefloat3d2, LINTIMEINTFACEFLOAT3D2) (
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
      } else {
         TBOX_ERROR(
            "FaceFloatLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
            << std::endl);
      }
   }
}

}
}
