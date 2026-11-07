/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for face-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianFaceDoubleWeightedAverage.h"

#include <cfloat>
#include <cmath>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/Utilities.h"

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

// in cartcoarsen1d.f:
void SAMRAI_F77_FUNC(cartwgtavgfacedoub1d, CARTWGTAVGFACEDOUB1D)(const int &,
                                                                 const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const int *, const double *, const double *,
                                                                 const double *, double *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgfacedoub2d0, CARTWGTAVGFACEDOUB2D0)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);

void SAMRAI_F77_FUNC(cartwgtavgfacedoub2d1, CARTWGTAVGFACEDOUB2D1)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgfacedoub3d0, CARTWGTAVGFACEDOUB3D0)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgfacedoub3d1, CARTWGTAVGFACEDOUB3D1)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgfacedoub3d2, CARTWGTAVGFACEDOUB3D2)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
// in cartcoarsen4d.f:
void SAMRAI_F77_FUNC(cartwgtavgfacedoub4d0, CARTWGTAVGFACEDOUB4D0)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgfacedoub4d1, CARTWGTAVGFACEDOUB4D1)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgfacedoub4d2, CARTWGTAVGFACEDOUB4D2)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgfacedoub4d3, CARTWGTAVGFACEDOUB4D3)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
}

namespace SAMRAI
{
namespace geom
{


CartesianFaceDoubleWeightedAverage::CartesianFaceDoubleWeightedAverage() : hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianFaceDoubleWeightedAverage::~CartesianFaceDoubleWeightedAverage()
{
}

int CartesianFaceDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianFaceDoubleWeightedAverage::getStencilWidth(const tbox::Dimension &dim) const
{
   return hier::IntVector::getZero(dim);
}

void CartesianFaceDoubleWeightedAverage::coarsen(
    hier::Patch &coarse,
    const hier::Patch &fine,
    const int dst_component,
    const int src_component,
    const hier::Box &coarse_box,
    const hier::IntVector &ratio) const
{
   const tbox::Dimension &dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::FaceData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::FaceData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::FaceData<double>, hier::PatchData>(
           coarse.getPatchData(dst_component)));
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Index &filo = fdata->getGhostBox().lower();
   const hier::Index &fihi = fdata->getGhostBox().upper();
   const hier::Index &cilo = cdata->getGhostBox().lower();
   const hier::Index &cihi = cdata->getGhostBox().upper();

   const std::shared_ptr<CartesianPatchGeometry> fgeom(
       SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
           fine.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> cgeom(
       SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
           coarse.getPatchGeometry()));

   TBOX_ASSERT(fgeom);
   TBOX_ASSERT(cgeom);

   const hier::Index &ifirstc = coarse_box.lower();
   const hier::Index &ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartwgtavgfacedoub1d, CARTWGTAVGFACEDOUB1D)
         (ifirstc(0),
          ilastc(0),
          filo(0), fihi(0),
          cilo(0), cihi(0),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(0, d),
          cdata->getPointer(0, d));
      } else if ((dim == tbox::Dimension(2))) {
#if defined(HAVE_RAJA)

         const double *fdx = fgeom->getDx();
         const double *cdx = cgeom->getDx();

         const double fdx0 = fdx[0];
         const double fdx1 = fdx[1];
         const double cdx0 = cdx[0];
         const double cdx1 = cdx[1];

         const int r0 = ratio[0];
         const int r1 = ratio[1];

         SAMRAI::hier::Box coarse_box_plus = coarse_box;
         coarse_box_plus.growUpper(0, 1);

         auto fine_array = fdata->getConstView<2>(0, d);
         auto coarse_array = cdata->getView<2>(0, d);

         double lengthf = fdx1;
         double lengthc = cdx1;

         hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
            double spv = 0.0;
            int jf = j * r0;
            for (int ry = 0; ry < r1; ry++) {
               int kf = k * r1 + ry;
               spv += fine_array(jf, kf) * lengthf;
            }

            coarse_array(j, k) = spv / lengthc;
         });

         auto fine_array_t = fdata->getConstView<2>(1, d);
         auto coarse_array_t = cdata->getView<2>(1, d);

         SAMRAI::hier::Box coarse_box_transpose = coarse_box;
         coarse_box_transpose.setLower(0, coarse_box.lower(1));
         coarse_box_transpose.setLower(1, coarse_box.lower(0));
         coarse_box_transpose.setUpper(0, coarse_box.upper(1));
         coarse_box_transpose.setUpper(1, coarse_box.upper(0));
         coarse_box_transpose.growUpper(0, 1);

         lengthf = fdx0;
         lengthc = cdx0;

         hier::parallel_for_all(coarse_box_transpose, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
            double spv = 0.0;
            int jf = j * r1;  // careful here, ratios are also switched
            for (int ry = 0; ry < r0; ry++) {
               int kf = k * r0 + ry;
               spv += fine_array_t(jf, kf) * lengthf;
            }

            coarse_array_t(j, k) = spv / lengthc;
         });
#else  // Fortran Dim 2
         SAMRAI_F77_FUNC(cartwgtavgfacedoub2d0, CARTWGTAVGFACEDOUB2D0)
         (ifirstc(0),
          ifirstc(1), ilastc(0), ilastc(1),
          filo(0), filo(1), fihi(0), fihi(1),
          cilo(0), cilo(1), cihi(0), cihi(1),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(0, d),
          cdata->getPointer(0, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub2d1, CARTWGTAVGFACEDOUB2D1)
         (ifirstc(0),
          ifirstc(1), ilastc(0), ilastc(1),
          filo(0), filo(1), fihi(0), fihi(1),
          cilo(0), cilo(1), cihi(0), cihi(1),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(1, d),
          cdata->getPointer(1, d));
#endif
      } else if ((dim == tbox::Dimension(3))) {
#if defined(HAVE_RAJA)

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

         SAMRAI::hier::Box coarse_box_plus = coarse_box;
         coarse_box_plus.growUpper(0, 1);

         auto fine_array = fdata->getConstView<3>(0, d);
         auto coarse_array = cdata->getView<3>(0, d);

         double areaf = fdx1 * fdx2;
         double areac = cdx1 * cdx2;

         hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
            double spv = 0.0;
            int ii = i * r0;
            for (int rz = 0; rz < r2; rz++) {
               for (int ry = 0; ry < r1; ry++) {
                  int kk = k * r2 + rz;
                  int jj = j * r1 + ry;
                  spv += fine_array(ii, jj, kk) * areaf;
               }
            }

            coarse_array(i, j, k) = spv / areac;
         });

         //transpose to 1,2,0
         SAMRAI::hier::Box coarse_box_t1 = coarse_box;
         coarse_box_t1.setLower(0, coarse_box.lower(1));
         coarse_box_t1.setLower(1, coarse_box.lower(2));
         coarse_box_t1.setLower(2, coarse_box.lower(0));
         coarse_box_t1.setUpper(0, coarse_box.upper(1));
         coarse_box_t1.setUpper(1, coarse_box.upper(2));
         coarse_box_t1.setUpper(2, coarse_box.upper(0));
         coarse_box_t1.growUpper(0, 1);

         auto fine_array_t1 = fdata->getConstView<3>(1, d);
         auto coarse_array_t1 = cdata->getView<3>(1, d);

         areaf = fdx2 * fdx0;
         areac = cdx2 * cdx0;


         hier::parallel_for_all(coarse_box_t1, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
            double spv = 0.0;
            int ii = i * r1;
            for (int rz = 0; rz < r0; rz++) {
               for (int ry = 0; ry < r2; ry++) {
                  int kk = k * r0 + rz;
                  int jj = j * r2 + ry;
                  spv += fine_array_t1(ii, jj, kk) * areaf;
               }
            }

            coarse_array_t1(i, j, k) = spv / areac;
         });

         //transpose to 2,0,1
         SAMRAI::hier::Box coarse_box_t2 = coarse_box;
         coarse_box_t2.setLower(0, coarse_box.lower(2));
         coarse_box_t2.setLower(1, coarse_box.lower(0));
         coarse_box_t2.setLower(2, coarse_box.lower(1));
         coarse_box_t2.setUpper(0, coarse_box.upper(2));
         coarse_box_t2.setUpper(1, coarse_box.upper(0));
         coarse_box_t2.setUpper(2, coarse_box.upper(1));
         coarse_box_t2.growUpper(0, 1);

         auto fine_array_t2 = fdata->getConstView<3>(2, d);
         auto coarse_array_t2 = cdata->getView<3>(2, d);

         areaf = fdx0 * fdx1;
         areac = cdx0 * cdx1;

         hier::parallel_for_all(coarse_box_t2, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
            double spv = 0.0;
            int ii = i * r2;
            for (int rz = 0; rz < r1; rz++) {
               for (int ry = 0; ry < r0; ry++) {
                  int kk = k * r1 + rz;
                  int jj = j * r0 + ry;
                  spv += fine_array_t2(ii, jj, kk) * areaf;
               }
            }

            coarse_array_t2(i, j, k) = spv / areac;
         });

#else  // Fortran dim 3
         SAMRAI_F77_FUNC(cartwgtavgfacedoub3d0, CARTWGTAVGFACEDOUB3D0)
         (ifirstc(0),
          ifirstc(1), ifirstc(2),
          ilastc(0), ilastc(1), ilastc(2),
          filo(0), filo(1), filo(2),
          fihi(0), fihi(1), fihi(2),
          cilo(0), cilo(1), cilo(2),
          cihi(0), cihi(1), cihi(2),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(0, d),
          cdata->getPointer(0, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub3d1, CARTWGTAVGFACEDOUB3D1)
         (ifirstc(0),
          ifirstc(1), ifirstc(2),
          ilastc(0), ilastc(1), ilastc(2),
          filo(0), filo(1), filo(2),
          fihi(0), fihi(1), fihi(2),
          cilo(0), cilo(1), cilo(2),
          cihi(0), cihi(1), cihi(2),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(1, d),
          cdata->getPointer(1, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub3d2, CARTWGTAVGFACEDOUB3D2)
         (ifirstc(0),
          ifirstc(1), ifirstc(2),
          ilastc(0), ilastc(1), ilastc(2),
          filo(0), filo(1), filo(2),
          fihi(0), fihi(1), fihi(2),
          cilo(0), cilo(1), cilo(2),
          cihi(0), cihi(1), cihi(2),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(2, d),
          cdata->getPointer(2, d));
#endif
      } else if ((dim == tbox::Dimension(4))) {
         SAMRAI_F77_FUNC(cartwgtavgfacedoub4d0, CARTWGTAVGFACEDOUB4D0)
         (ifirstc(0),
          ifirstc(1), ifirstc(2), ifirstc(3),
          ilastc(0), ilastc(1), ilastc(2), ilastc(3),
          filo(0), filo(1), filo(2), filo(3),
          fihi(0), fihi(1), fihi(2), fihi(3),
          cilo(0), cilo(1), cilo(2), cilo(3),
          cihi(0), cihi(1), cihi(2), cihi(3),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(0, d),
          cdata->getPointer(0, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub4d1, CARTWGTAVGFACEDOUB4D1)
         (ifirstc(0),
          ifirstc(1), ifirstc(2), ifirstc(3),
          ilastc(0), ilastc(1), ilastc(2), ilastc(3),
          filo(0), filo(1), filo(2), filo(3),
          fihi(0), fihi(1), fihi(2), fihi(3),
          cilo(0), cilo(1), cilo(2), cilo(3),
          cihi(0), cihi(1), cihi(2), cihi(3),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(1, d),
          cdata->getPointer(1, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub4d2, CARTWGTAVGFACEDOUB4D2)
         (ifirstc(0),
          ifirstc(1), ifirstc(2), ifirstc(3),
          ilastc(0), ilastc(1), ilastc(2), ilastc(3),
          filo(0), filo(1), filo(2), filo(3),
          fihi(0), fihi(1), fihi(2), fihi(3),
          cilo(0), cilo(1), cilo(2), cilo(3),
          cihi(0), cihi(1), cihi(2), cihi(3),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(2, d),
          cdata->getPointer(2, d));
         SAMRAI_F77_FUNC(cartwgtavgfacedoub4d3, CARTWGTAVGFACEDOUB4D3)
         (ifirstc(0),
          ifirstc(1), ifirstc(2), ifirstc(3),
          ilastc(0), ilastc(1), ilastc(2), ilastc(3),
          filo(0), filo(1), filo(2), filo(3),
          fihi(0), fihi(1), fihi(2), fihi(3),
          cilo(0), cilo(1), cilo(2), cilo(3),
          cihi(0), cihi(1), cihi(2), cihi(3),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(3, d),
          cdata->getPointer(3, d));
      } else {
         TBOX_ERROR("CartesianFaceDoubleWeightedAverage error...\n"
                    << "dim > 4 not supported." << std::endl);
      }
   }
}

}  // namespace geom
}  // namespace SAMRAI
