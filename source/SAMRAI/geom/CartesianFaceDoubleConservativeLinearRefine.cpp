/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for face-centered
 *                double data on a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include <cfloat>
#include <cmath>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Utilities.h"

#define SAMRAI_GEOM_MIN(a, b) (((b) < (a)) ? (b) : (a))

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

// in cartrefine1d.f:
void SAMRAI_F77_FUNC(cartclinreffacedoub1d, CARTCLINREFFACEDOUB1D)(const int &,
                                                                   const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *,
                                                                   double *, double *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinreffacedoub2d0, CARTCLINREFFACEDOUB2D0)(const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int *, const double *, const double *,
                                                                     const double *, double *,
                                                                     double *, double *, double *, double *);
void SAMRAI_F77_FUNC(cartclinreffacedoub2d1, CARTCLINREFFACEDOUB2D1)(const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int &, const int &, const int &, const int &,
                                                                     const int *, const double *, const double *,
                                                                     const double *, double *,
                                                                     double *, double *, double *, double *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinreffacedoub3d0, CARTCLINREFFACEDOUB3D0)(const int &,
                                                                     const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int *, const double *, const double *,
                                                                     const double *, double *,
                                                                     double *, double *, double *,
                                                                     double *, double *, double *);
void SAMRAI_F77_FUNC(cartclinreffacedoub3d1, CARTCLINREFFACEDOUB3D1)(const int &,
                                                                     const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int *, const double *, const double *,
                                                                     const double *, double *,
                                                                     double *, double *, double *,
                                                                     double *, double *, double *);
void SAMRAI_F77_FUNC(cartclinreffacedoub3d2, CARTCLINREFFACEDOUB3D2)(const int &,
                                                                     const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int &, const int &, const int &,
                                                                     const int *, const double *, const double *,
                                                                     const double *, double *,
                                                                     double *, double *, double *,
                                                                     double *, double *, double *);
}

namespace SAMRAI
{
namespace geom
{


CartesianFaceDoubleConservativeLinearRefine::
    CartesianFaceDoubleConservativeLinearRefine() : hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianFaceDoubleConservativeLinearRefine::~CartesianFaceDoubleConservativeLinearRefine()
{
}

int CartesianFaceDoubleConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianFaceDoubleConservativeLinearRefine::getStencilWidth(const tbox::Dimension &dim) const
{
   return hier::IntVector::getOne(dim);
}

void CartesianFaceDoubleConservativeLinearRefine::refine(
    hier::Patch &fine,
    const hier::Patch &coarse,
    const int dst_component,
    const int src_component,
    const hier::BoxOverlap &fine_overlap,
    const hier::IntVector &ratio) const
{
   const tbox::Dimension &dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim, coarse, ratio);

   std::shared_ptr<pdat::FaceData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
           coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::FaceData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
           fine.getPatchData(dst_component)));

   const pdat::FaceOverlap *t_overlap =
       CPP_CAST<const pdat::FaceOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index &cilo = cgbox.lower();
   const hier::Index &cihi = cgbox.upper();
   const hier::Index &filo = fdata->getGhostBox().lower();
   const hier::Index &fihi = fdata->getGhostBox().upper();

   const std::shared_ptr<CartesianPatchGeometry> cgeom(
       SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
           coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
       SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
           fine.getPatchGeometry()));

   for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {
      const hier::BoxContainer &boxes = t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {

         const hier::Box &face_box = *b;
         TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, face_box);

         hier::Box fine_box(dim);
         for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
            fine_box.setLower(
                static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()), face_box.lower(i));
            fine_box.setUpper(
                static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()), face_box.upper(i));
         }

         fine_box.setUpper(axis, fine_box.upper(axis) - 1);

         const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
         const hier::Index &ifirstc = coarse_box.lower();
         const hier::Index &ilastc = coarse_box.upper();
         const hier::Index &ifirstf = fine_box.lower();
         const hier::Index &ilastf = fine_box.upper();

         const hier::IntVector tmp_ghosts(dim, 0);
         std::vector<double> diff0_f(cgbox.numberCells(0) + 2);

         tbox::AllocatorDatabase *alloc_db = tbox::AllocatorDatabase::getDatabase();
         pdat::FaceData<double> slope0_f(cgbox, 1, tmp_ghosts, alloc_db->getTagAllocator());

         for (int d = 0; d < fdata->getDepth(); ++d) {
            if ((dim == tbox::Dimension(1))) {
               SAMRAI_F77_FUNC(cartclinreffacedoub1d, CARTCLINREFFACEDOUB1D)
               (
                   ifirstc(0), ilastc(0),
                   ifirstf(0), ilastf(0),
                   cilo(0), cihi(0),
                   filo(0), fihi(0),
                   &ratio[0],
                   cgeom->getDx(),
                   fgeom->getDx(),
                   cdata->getPointer(0, d),
                   fdata->getPointer(0, d),
                   &diff0_f[0], slope0_f.getPointer(0));
            } else if ((dim == tbox::Dimension(2))) {
#if defined(HAVE_RAJA)
               SAMRAI::hier::Box fine_box_plus = fine_box;
               SAMRAI::hier::Box diff_box = coarse_box;

               // Iteration space is slightly different between the directions

               if (axis == 1) {  // transpose boxes
                  fine_box_plus.setLower(0, fine_box.lower(1));
                  fine_box_plus.setLower(1, fine_box.lower(0));
                  fine_box_plus.setUpper(0, fine_box.upper(1));
                  fine_box_plus.setUpper(1, fine_box.upper(0));

                  diff_box.setLower(0, coarse_box.lower(1));
                  diff_box.setLower(1, coarse_box.lower(0));
                  diff_box.setUpper(0, coarse_box.upper(1));
                  diff_box.setUpper(1, coarse_box.upper(0));
               }
               SAMRAI::hier::Box slope_box = diff_box;

               fine_box_plus.growUpper(0, 1);
               diff_box.grow(0, 1);
               diff_box.growUpper(1, 1);

               if (axis == 0) {
                  slope_box.growUpper(0, 1);
               }

               // Array data setup for diff and slope uses dim as depth value and is not related to actual depth being processed
               // so that we have for 2d two components with the same box extents namely diff0, diff1
               // Lifetime of this Array data exits just for a given depth component being processed
               // We may want to redo this approach to avoid alloc/dealloc redundancies
               pdat::ArrayData<double> diff(diff_box, dim.getValue(), alloc_db->getDevicePool());
               pdat::ArrayData<double> slope(slope_box, dim.getValue(), alloc_db->getDevicePool());

               auto fine_array = fdata->getView<2>(axis, d);
               auto coarse_array = cdata->getConstView<2>(axis, d);

               auto diff0 = diff.getView<2>(0);
               auto diff1 = diff.getView<2>(1);

               auto slope0 = slope.getView<2>(0);
               auto slope1 = slope.getView<2>(1);

               const double *fdx = fgeom->getDx();
               const double *cdx = cgeom->getDx();
               const double fdx0 = fdx[0];
               const double fdx1 = fdx[1];
               const double cdx0 = cdx[0];
               const double cdx1 = cdx[1];

               const int r0 = ratio[0];
               const int r1 = ratio[1];
               if (axis == 0) {
                  hier::parallel_for_all(diff_box, [=] SAMRAI_HOST_DEVICE(int j /*fast*/, int k /*slow */) {
                     diff0(j, k) = coarse_array(j + 1, k) - coarse_array(j, k);
                     diff1(j, k) = coarse_array(j, k) - coarse_array(j, k - 1);
                  });

                  hier::parallel_for_all(slope_box, [=] SAMRAI_HOST_DEVICE(int j, int k) {
                     const double coef2j = 0.5 * (diff0(j - 1, k) + diff0(j, k));
                     const double boundj = 2.0 * SAMRAI_GEOM_MIN(fabs(diff0(j - 1, k)), fabs(diff0(j, k)));

                     if (diff0(j, k) * diff0(j - 1, k) > 0.0 && cdx0 != 0) {
                        slope0(j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2j), boundj), coef2j) / cdx0;
                     } else {
                        slope0(j, k) = 0.0;
                     }

                     const double coef2k = 0.5 * (diff1(j, k + 1) + diff1(j, k));
                     const double boundk = 2.0 * SAMRAI_GEOM_MIN(fabs(diff1(j, k + 1)), fabs(diff1(j, k)));

                     if (diff1(j, k) * diff1(j - 1, k) > 0.0 && cdx1 != 0) {
                        slope1(j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2k), boundk), coef2k) / cdx1;
                     } else {
                        slope1(j, k) = 0.0;
                     }
                  });

                  hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int j, int k) {
                     const int ic1 = (k < 0) ? (k + 1) / r1 - 1 : k / r1;
                     const int ic0 = (j < 0) ? (j + 1) / r0 - 1 : j / r0;

                     const int ir0 = j - ic0 * r0;
                     const int ir1 = k - ic1 * r1;
                     double deltax0, deltax1;

                     deltax1 = (static_cast<double>(ir1) + 0.5) * fdx1 - cdx1 * 0.5;
                     deltax0 = static_cast<double>(ir0) * fdx0;

                     fine_array(j, k) = coarse_array(ic0, ic1) + slope0(ic0, ic1) * deltax0 + slope1(ic0, ic1) * deltax1;
                  });

               }  // axis == 0
               else if (axis == 1) {
                  hier::parallel_for_all(diff_box, [=] SAMRAI_HOST_DEVICE(int j /*fast*/, int k /*slow */) {
                     diff0(j, k) = coarse_array(j, k) - coarse_array(j, k - 1);
                     diff1(j, k) = coarse_array(j + 1, k) - coarse_array(j, k);
                  });

                  hier::parallel_for_all(slope_box, [=] SAMRAI_HOST_DEVICE(int j, int k) {
                     const double coef2j = 0.5 * (diff0(j, k + 1) + diff0(j, k));
                     const double boundj = 2.0 * SAMRAI_GEOM_MIN(fabs(diff0(j, k + 1)), fabs(diff0(j, k)));

                     if (diff0(j, k) * diff0(j, k + 1) > 0.0 && cdx0 != 0) {
                        slope0(j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2j), boundj), coef2j) / cdx0;
                     } else {
                        slope0(j, k) = 0.0;
                     }

                     const double coef2k = 0.5 * (diff1(j - 1, k) + diff1(j, k));
                     const double boundk = 2.0 * SAMRAI_GEOM_MIN(fabs(diff1(j - 1, k)), fabs(diff1(j, k)));

                     if (diff1(j, k) * diff1(j - 1, k) > 0.0 && cdx1 != 0) {
                        slope1(j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2k), boundk), coef2k) / cdx1;
                     } else {
                        slope1(j, k) = 0.0;
                     }
                  });

                  hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int j, int k) {
                     const int ic0 = (k < 0) ? (k + 1) / r1 - 1 : k / r0;
                     const int ic1 = (j < 0) ? (j + 1) / r0 - 1 : j / r1;

                     const int ir0 = k - ic0 * r0;
                     const int ir1 = j - ic1 * r1;
                     double deltax0, deltax1;

                     deltax0 = (static_cast<double>(ir0) + 0.5) * fdx0 - cdx0 * 0.5;
                     deltax1 = static_cast<double>(ir1) * fdx1;

                     fine_array(j, k) = coarse_array(ic1, ic0) + slope1(ic1, ic0) * deltax1 + slope0(ic1, ic0) * deltax0;
                  });

               }  // axis == 1

#else  // Fortran Dimension 2
               std::vector<double> diff1_f(cgbox.numberCells(1) + 2);
               pdat::FaceData<double> slope1_f(cgbox, 1, tmp_ghosts, alloc_db->getTagAllocator());

               if (axis == 0) {
                  SAMRAI_F77_FUNC(cartclinreffacedoub2d0, CARTCLINREFFACEDOUB2D0)
                  (
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cgeom->getDx(),
                      fgeom->getDx(),
                      cdata->getPointer(0, d),
                      fdata->getPointer(0, d),
                      &diff0_f[0], slope0_f.getPointer(0),
                      &diff1_f[0], slope1_f.getPointer(0));
               } else if (axis == 1) {
                  SAMRAI_F77_FUNC(cartclinreffacedoub2d1, CARTCLINREFFACEDOUB2D1)
                  (
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cgeom->getDx(),
                      fgeom->getDx(),
                      cdata->getPointer(1, d),
                      fdata->getPointer(1, d),
                      &diff1_f[0], slope1_f.getPointer(1),
                      &diff0_f[0], slope0_f.getPointer(1));
               }
#endif
            } else if ((dim == tbox::Dimension(3))) {
#if defined(HAVE_RAJA)
               SAMRAI::hier::Box fine_box_plus = fine_box;
               SAMRAI::hier::Box diff_box = coarse_box;

               if (axis == 1) {  // transpose boxes <1,2,0>
                  fine_box_plus.setLower(0, fine_box.lower(1));
                  fine_box_plus.setLower(1, fine_box.lower(2));
                  fine_box_plus.setLower(2, fine_box.lower(0));

                  fine_box_plus.setUpper(0, fine_box.upper(1));
                  fine_box_plus.setUpper(1, fine_box.upper(2));
                  fine_box_plus.setUpper(2, fine_box.upper(0));

                  diff_box.setLower(0, coarse_box.lower(1));
                  diff_box.setLower(1, coarse_box.lower(2));
                  diff_box.setLower(2, coarse_box.lower(0));

                  diff_box.setUpper(0, coarse_box.upper(1));
                  diff_box.setUpper(1, coarse_box.upper(2));
                  diff_box.setUpper(2, coarse_box.upper(0));
               } else if (axis == 2) {  // <2,0,1>
                  fine_box_plus.setLower(0, fine_box.lower(2));
                  fine_box_plus.setLower(1, fine_box.lower(0));
                  fine_box_plus.setLower(2, fine_box.lower(1));

                  fine_box_plus.setUpper(0, fine_box.upper(2));
                  fine_box_plus.setUpper(1, fine_box.upper(0));
                  fine_box_plus.setUpper(2, fine_box.upper(1));

                  diff_box.setLower(0, coarse_box.lower(2));
                  diff_box.setLower(1, coarse_box.lower(0));
                  diff_box.setLower(2, coarse_box.lower(1));

                  diff_box.setUpper(0, coarse_box.upper(2));
                  diff_box.setUpper(1, coarse_box.upper(0));
                  diff_box.setUpper(2, coarse_box.upper(1));
               }

               SAMRAI::hier::Box slope_box = diff_box;
               slope_box.growUpper(0, 1);

               fine_box_plus.growUpper(0, 1);

               diff_box.grow(0, 1);
               diff_box.growUpper(1, 1);
               diff_box.growUpper(2, 1);

               pdat::ArrayData<double> diff(diff_box, dim.getValue(), alloc_db->getDevicePool());
               pdat::ArrayData<double> slope(slope_box, dim.getValue(), alloc_db->getDevicePool());

               auto fine_array = fdata->getView<3>(axis, d);
               auto coarse_array = cdata->getConstView<3>(axis, d);

               auto diff0 = diff.getView<3>(0);
               auto diff1 = diff.getView<3>(1);
               auto diff2 = diff.getView<3>(2);

               auto slope0 = slope.getView<3>(0);
               auto slope1 = slope.getView<3>(1);
               auto slope2 = slope.getView<3>(2);

               const double *fdx = fgeom->getDx();
               const double *cdx = cgeom->getDx();
               const double fdx0 = fdx[0];
               const double fdx1 = fdx[1];
               const double fdx2 = fdx[2];
               const double cdx0 = cdx[0];
               const double cdx1 = cdx[1];
               const double cdx2 = cdx[2];

               const int r0 = ratio[0];
               const int r1 = ratio[1];
               const int r2 = ratio[2];

               if (axis == 0) {

                  hier::parallel_for_all(diff_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k /*slow */) {
                     diff0(i, j, k) = coarse_array(i + 1, j, k) - coarse_array(i, j, k);
                     diff1(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j - 1, k);
                     diff2(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j, k - 1);
                  });

                  hier::parallel_for_all(slope_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     const double coef2i = 0.5 * (diff0(i - 1, j, k) + diff0(i, j, k));
                     const double boundi = 2.0 * SAMRAI_GEOM_MIN(fabs(diff0(i - 1, j, k)), fabs(diff0(i, j, k)));

                     if (diff0(i, j, k) * diff0(i - 1, j, k) > 0.0 && cdx0 != 0) {
                        slope0(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2i), boundi), coef2i) / cdx0;
                     } else {
                        slope0(i, j, k) = 0.0;
                     }

                     const double coef2j = 0.5 * (diff1(i, j + 1, k) + diff1(i, j, k));
                     const double boundj = 2.0 * SAMRAI_GEOM_MIN(fabs(diff1(i, j + 1, k)), fabs(diff1(i, j, k)));

                     if (diff1(i, j, k) * diff1(i, j + 1, k) > 0.0 && cdx1 != 0) {
                        slope1(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2j), boundj), coef2j) / cdx1;
                     } else {
                        slope1(i, j, k) = 0.0;
                     }

                     const double coef2k = 0.5 * (diff2(i, j, k + 1) + diff2(i, j, k));
                     const double boundk = 2.0 * SAMRAI_GEOM_MIN(fabs(diff2(i, j, k + 1)), fabs(diff2(i, j, k)));

                     if (diff2(i, j, k) * diff2(i, j, k + 1) > 0.0 && cdx2 != 0) {
                        slope2(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2k), boundk), coef2k) / cdx2;
                     } else {
                        slope2(i, j, k) = 0.0;
                     }
                  });

                  hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;
                     const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
                     const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;

                     const int ir0 = i - ic0 * r0;
                     const int ir1 = j - ic1 * r1;
                     const int ir2 = k - ic2 * r2;

                     double deltax0, deltax1, deltax2;
                     deltax0 = static_cast<double>(ir0) * fdx0;
                     deltax1 = (static_cast<double>(ir1) + 0.5) * fdx1 - cdx1 * 0.5;
                     deltax2 = (static_cast<double>(ir2) + 0.5) * fdx2 - cdx2 * 0.5;

                     fine_array(i, j, k) = coarse_array(ic0, ic1, ic2) + slope0(ic0, ic1, ic2) * deltax0 + slope1(ic0, ic1, ic2) * deltax1 + slope2(ic0, ic1, ic2) * deltax2;
                  });

               }                      // done axis == 0
               else if (axis == 1) {  //1,2,0
                  hier::parallel_for_all(diff_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k /*slow */) {
                     diff1(i, j, k) = coarse_array(i + 1, j, k) - coarse_array(i, j, k);
                     diff2(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j - 1, k);
                     diff0(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j, k - 1);
                  });

                  hier::parallel_for_all(slope_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     const double coef2i = 0.5 * (diff0(i, j, k + 1) + diff0(i, j, k));
                     const double boundi = 2.0 * SAMRAI_GEOM_MIN(fabs(diff0(i, j, k + 1)), fabs(diff0(i, j, k)));

                     if (diff0(i, j, k) * diff0(i, j, k + 1) > 0.0 && cdx0 != 0) {
                        slope0(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2i), boundi), coef2i) / cdx0;
                     } else {
                        slope0(i, j, k) = 0.0;
                     }

                     const double coef2j = 0.5 * (diff1(i - 1, j, k) + diff1(i, j, k));
                     const double boundj = 2.0 * SAMRAI_GEOM_MIN(fabs(diff1(i - 1, j, k)), fabs(diff1(i, j, k)));

                     if (diff1(i, j, k) * diff1(i - 1, j, k) > 0.0 && cdx1 != 0) {
                        slope1(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2j), boundj), coef2j) / cdx1;
                     } else {
                        slope1(i, j, k) = 0.0;
                     }

                     const double coef2k = 0.5 * (diff2(i, j + 1, k) + diff2(i, j, k));
                     const double boundk = 2.0 * SAMRAI_GEOM_MIN(fabs(diff2(i, j + 1, k)), fabs(diff2(i, j, k)));

                     if (diff2(i, j, k) * diff2(i, j + 1, k) > 0.0 && cdx2 != 0) {
                        slope2(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2k), boundk), coef2k) / cdx2;
                     } else {
                        slope2(i, j, k) = 0.0;
                     }
                  });

                  hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     // keep ic0 - ic2 consistent with i,j,k; just change ir0-ir2;  TODO redo 2dim case
                     const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;
                     const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
                     const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;

                     //1,2,0
                     const int ir1 = i - ic0 * r0;
                     const int ir2 = j - ic1 * r1;
                     const int ir0 = k - ic2 * r2;

                     double deltax0, deltax1, deltax2;
                     deltax0 = (static_cast<double>(ir0) + 0.5) * fdx0 - cdx0 * 0.5;
                     deltax1 = static_cast<double>(ir1) * fdx1;
                     deltax2 = (static_cast<double>(ir2) + 0.5) * fdx2 - cdx2 * 0.5;

                     fine_array(i, j, k) = coarse_array(ic0, ic1, ic2) + slope0(ic0, ic1, ic2) * deltax0 + slope1(ic0, ic1, ic2) * deltax1 + slope2(ic0, ic1, ic2) * deltax2;
                  });

               } else if (axis == 2) {  // 2,0,1

                  hier::parallel_for_all(diff_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k /*slow */) {
                     diff2(i, j, k) = coarse_array(i + 1, j, k) - coarse_array(i, j, k);
                     diff0(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j - 1, k);
                     diff1(i, j, k) = coarse_array(i, j, k) - coarse_array(i, j, k - 1);
                  });

                  hier::parallel_for_all(slope_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     const double coef2i = 0.5 * (diff0(i, j + 1, k) + diff0(i, j, k));
                     const double boundi = 2.0 * SAMRAI_GEOM_MIN(fabs(diff0(i, j + 1, k)), fabs(diff0(i, j, k)));

                     if (diff0(i, j, k) * diff0(i, j + 1, k) > 0.0 && cdx0 != 0) {
                        slope0(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2i), boundi), coef2i) / cdx0;
                     } else {
                        slope0(i, j, k) = 0.0;
                     }

                     const double coef2j = 0.5 * (diff1(i, j, k + 1) + diff1(i, j, k));
                     const double boundj = 2.0 * SAMRAI_GEOM_MIN(fabs(diff1(i, j, k + 1)), fabs(diff1(i, j, k)));

                     if (diff1(i, j, k) * diff1(i, j, k + 1) > 0.0 && cdx1 != 0) {
                        slope1(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2j), boundj), coef2j) / cdx1;
                     } else {
                        slope1(i, j, k) = 0.0;
                     }

                     const double coef2k = 0.5 * (diff2(i - 1, j, k) + diff2(i, j, k));
                     const double boundk = 2.0 * SAMRAI_GEOM_MIN(fabs(diff2(i - 1, j, k)), fabs(diff2(i, j, k)));

                     if (diff2(i, j, k) * diff2(i - 1, j, k) > 0.0 && cdx2 != 0) {
                        slope2(i, j, k) = copysign(SAMRAI_GEOM_MIN(fabs(coef2k), boundk), coef2k) / cdx2;
                     } else {
                        slope2(i, j, k) = 0.0;
                     }
                  });

                  hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                     // keep ic0 - ic2 consistent with i,j,k; just change ir0-ir2;  redo 2dim case
                     const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;
                     const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
                     const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;

                     //2,0,1
                     const int ir2 = i - ic0 * r0;
                     const int ir0 = j - ic1 * r1;
                     const int ir1 = k - ic2 * r2;

                     double deltax0, deltax1, deltax2;
                     deltax0 = (static_cast<double>(ir0) + 0.5) * fdx0 - cdx0 * 0.5;
                     deltax1 = (static_cast<double>(ir1) + 0.5) * fdx1 - cdx1 * 0.5;
                     deltax2 = static_cast<double>(ir2) * fdx2;

                     fine_array(i, j, k) = coarse_array(ic0, ic1, ic2) + slope0(ic0, ic1, ic2) * deltax0 + slope1(ic0, ic1, ic2) * deltax1 + slope2(ic0, ic1, ic2) * deltax2;
                  });
               }

#else
                      // Fortran Dimension3
                      // Iteration space is slightly different between the directions
               std::vector<double> diff1_f(cgbox.numberCells(1) + 2);
               std::vector<double> diff2_f(cgbox.numberCells(2) + 2);

               pdat::FaceData<double> slope1_f(cgbox, 1, tmp_ghosts, alloc_db->getTagAllocator());
               pdat::FaceData<double> slope2_f(cgbox, 1, tmp_ghosts, alloc_db->getTagAllocator());

               if (axis == 0) {
                  SAMRAI_F77_FUNC(cartclinreffacedoub3d0, CARTCLINREFFACEDOUB3D0)
                  (
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cgeom->getDx(),
                      fgeom->getDx(),
                      cdata->getPointer(0, d),
                      fdata->getPointer(0, d),
                      &diff0_f[0], slope0_f.getPointer(0),
                      &diff1_f[0], slope1_f.getPointer(0),
                      &diff2_f[0], slope2_f.getPointer(0));
               } else if (axis == 1) {
                  SAMRAI_F77_FUNC(cartclinreffacedoub3d1, CARTCLINREFFACEDOUB3D1)
                  (
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cgeom->getDx(),
                      fgeom->getDx(),
                      cdata->getPointer(1, d),
                      fdata->getPointer(1, d),
                      &diff1_f[0], slope1_f.getPointer(1),
                      &diff2_f[0], slope2_f.getPointer(1),
                      &diff0_f[0], slope0_f.getPointer(1));
               } else if (axis == 2) {
                  SAMRAI_F77_FUNC(cartclinreffacedoub3d2, CARTCLINREFFACEDOUB3D2)
                  (
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cgeom->getDx(),
                      fgeom->getDx(),
                      cdata->getPointer(2, d),
                      fdata->getPointer(2, d),
                      &diff2_f[0], slope2_f.getPointer(2),
                      &diff0_f[0], slope0_f.getPointer(2),
                      &diff1_f[0], slope1_f.getPointer(2));
               }
#endif                // test if defined RAJA
            } else {  // Not dimension 3
               TBOX_ERROR(
                   "CartesianFaceDoubleConservativeLinearRefine error...\n"
                   << "dim > 3 not supported." << std::endl);
            }
         }  // for depth
      }     // for boxes
   }        // for axis
}  // end refine call

}  // end namespace geom
}  // end namespace SAMRAI
