/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for outerside double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
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
void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub1d, CARTWGTAVGOUTSIDEDOUB1D)(const int &,
                                                                       const int &,
                                                                       const int &, const int &,
                                                                       const int &, const int &,
                                                                       const int *, const double *, const double *,
                                                                       const double *, double *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub2d0, CARTWGTAVGOUTSIDEDOUB2D0)(const int &,
                                                                         const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int *, const double *, const double *,
                                                                         const double *, double *);

void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub2d1, CARTWGTAVGOUTSIDEDOUB2D1)(const int &,
                                                                         const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &,
                                                                         const int *, const double *, const double *,
                                                                         const double *, double *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d0, CARTWGTAVGOUTSIDEDOUB3D0)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int *, const double *, const double *,
                                                                         const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d1, CARTWGTAVGOUTSIDEDOUB3D1)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int *, const double *, const double *,
                                                                         const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d2, CARTWGTAVGOUTSIDEDOUB3D2)(const int &,
                                                                         const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int &, const int &, const int &,
                                                                         const int *, const double *, const double *,
                                                                         const double *, double *);
}

namespace SAMRAI
{
namespace geom
{


CartesianOutersideDoubleWeightedAverage::
    CartesianOutersideDoubleWeightedAverage() : hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianOutersideDoubleWeightedAverage::~CartesianOutersideDoubleWeightedAverage()
{
}

int CartesianOutersideDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianOutersideDoubleWeightedAverage::getStencilWidth(const tbox::Dimension &dim) const
{
   return hier::IntVector::getZero(dim);
}

void CartesianOutersideDoubleWeightedAverage::coarsen(
    hier::Patch &coarse,
    const hier::Patch &fine,
    const int dst_component,
    const int src_component,
    const hier::Box &coarse_box,
    const hier::IntVector &ratio) const
{
   const tbox::Dimension &dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::OutersideData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::OutersideData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
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
      // loop over lower and upper outerside arrays
      for (int side = 0; side < 2; ++side) {
         if ((dim == tbox::Dimension(1))) {
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub1d,
                            CARTWGTAVGOUTSIDEDOUB1D)
            (ifirstc(0), ilastc(0),
             filo(0), fihi(0),
             cilo(0), cihi(0),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(0, side, d),
             cdata->getPointer(0, side, d));
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

            int jf0bounds;  // setup capture variable representing fine lower/upper bounds depending on side for side-normal 0
            int kf1bounds;  // and similarly for side-normal 1

            // setup side-normal boxes
            SAMRAI::hier::Box coarse_box_sn0 = coarse_box;
            SAMRAI::hier::Box coarse_box_sn1 = coarse_box;

            if (side == 0) {
               coarse_box_sn0.setLower(0, coarse_box.lower(0));
               coarse_box_sn0.setUpper(0, coarse_box.lower(0));
               coarse_box_sn1.setLower(1, coarse_box.lower(1));
               coarse_box_sn1.setUpper(1, coarse_box.lower(1));
               jf0bounds = filo(0);  // for face-normal 0
               kf1bounds = filo(1);  // for face-normal 1
            } else { // side == 1
               coarse_box_sn0.setLower(0, coarse_box.upper(0));
               coarse_box_sn0.setUpper(0, coarse_box.upper(0));
               coarse_box_sn1.setLower(1, coarse_box.upper(1));
               coarse_box_sn1.setUpper(1, coarse_box.upper(1));
               jf0bounds = fihi(0);
               kf1bounds = fihi(1);
            }

            auto fine_array_0 = fdata->getConstView<2>(0, side, d);
            auto coarse_array_0 = cdata->getView<2>(0, side, d);

            double lengthf = fdx1;
            double lengthc = cdx1;

            hier::parallel_for_all(coarse_box_sn0, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
               double spv = 0.0;
               int jf = jf0bounds;
               for (int ry = 0; ry < r1; ry++) {
                  int kf = k * r1 + ry;
                  spv += fine_array_0(jf, kf) * lengthf;
               }

               coarse_array_0(j, k) = spv / lengthc;
            });

            auto fine_array_1 = fdata->getConstView<2>(1, side, d);
            auto coarse_array_1 = cdata->getView<2>(1, side, d);

            lengthf = fdx0;
            lengthc = cdx0;

            hier::parallel_for_all(coarse_box_sn1, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
               double spv = 0.0;
               int kf = kf1bounds;
               for (int rx = 0; rx < r0; rx++) {
                  int jf = j * r0 + rx;
                  spv += fine_array_1(jf, kf) * lengthf;
               }
               coarse_array_1(j, k) = spv / lengthc;
            });

#else  // Fortran Dim 2
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub2d0,
                            CARTWGTAVGOUTSIDEDOUB2D0)
            (ifirstc(0), ifirstc(1), ilastc(0),
             ilastc(1),
             filo(0), filo(1), fihi(0), fihi(1),
             cilo(0), cilo(1), cihi(0), cihi(1),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(0, side, d),
             cdata->getPointer(0, side, d));
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub2d1,
                            CARTWGTAVGOUTSIDEDOUB2D1)
            (ifirstc(0), ifirstc(1), ilastc(0),
             ilastc(1),
             filo(0), filo(1), fihi(0), fihi(1),
             cilo(0), cilo(1), cihi(0), cihi(1),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(1, side, d),
             cdata->getPointer(1, side, d));
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
            int if0bounds;  // setup capture variable representing fine lower/upper bounds depending on side for side-normal 0
            int jf1bounds;  // and similarly for side-normal 1
            int kf2bounds;  // and similarly for side-normal 2

            // setup side-normal boxes
            SAMRAI::hier::Box coarse_box_sn0 = coarse_box;
            SAMRAI::hier::Box coarse_box_sn1 = coarse_box;
            SAMRAI::hier::Box coarse_box_sn2 = coarse_box;

            if (side == 0) {
               coarse_box_sn0.setUpper(0, coarse_box.lower(0));
               if0bounds = filo(0);  // for side-normal 0
            } else { // side == 1
               coarse_box_sn0.setLower(0, coarse_box.upper(0));
               coarse_box_sn0.setUpper(0, coarse_box.upper(0));
               if0bounds = fihi(0);
            }

            auto fine_array_0 = fdata->getConstView<3>(0, side, d);
            auto coarse_array_0 = cdata->getView<3>(0, side, d);

            double areaf = fdx1 * fdx2;
            double areac = cdx1 * cdx2;

            hier::parallel_for_all(coarse_box_sn0, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
               double spv = 0.0;
               int ii = if0bounds;
               for (int rz = 0; rz < r2; rz++) {
                  for (int ry = 0; ry < r1; ry++) {
                     int kk = k * r2 + rz;
                     int jj = j * r1 + ry;
                     spv += fine_array_0(ii, jj, kk) * areaf;
                  }
               }

               coarse_array_0(i, j, k) = spv / areac;
            });

            if (side == 0) {
               coarse_box_sn1.setUpper(1, coarse_box.lower(1));
               jf1bounds = filo(1);  // for face-normal 1
            } else { // side == 1
               coarse_box_sn1.setLower(1, coarse_box.upper(1));
               jf1bounds = fihi(1);
            }

            auto fine_array_1 = fdata->getConstView<3>(1, side, d);
            auto coarse_array_1 = cdata->getView<3>(1, side, d);

            areaf = fdx2 * fdx0;
            areac = cdx2 * cdx0;

            hier::parallel_for_all(coarse_box_sn1, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
               double spv = 0.0;
               int jj = jf1bounds;
               for (int rx = 0; rx < r0; rx++) {
                  for (int rz = 0; rz < r2; rz++) {
                     int ii = i * r2 + rz;
                     int kk = k * r0 + rx;
                     spv += fine_array_1(ii, jj, kk) * areaf;
                  }
               }

               coarse_array_1(i, j, k) = spv / areac;
            });

            if (side == 0) {
               coarse_box_sn2.setUpper(2, coarse_box.lower(2));
               kf2bounds = filo(2);  // for side-normal 2
            } else { // side == 1
               coarse_box_sn2.setLower(2, coarse_box.upper(2));
               coarse_box_sn2.setUpper(2, coarse_box.upper(2));
               kf2bounds = fihi(2);
            }

            auto fine_array_2 = fdata->getConstView<3>(2, side, d);
            auto coarse_array_2 = cdata->getView<3>(2, side, d);

            areaf = fdx0 * fdx1;
            areac = cdx0 * cdx1;

            hier::parallel_for_all(coarse_box_sn2, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
               double spv = 0.0;
               int kk = kf2bounds;
               for (int ry = 0; ry < r1; ry++) {
                  for (int rx = 0; rx < r0; rx++) {
                     int ii = i * r0 + rx;
                     int jj = j * r1 + ry;
                     spv += fine_array_2(ii, jj, kk) * areaf;
                  }
               }

               coarse_array_2(i, j, k) = spv / areac;
            });

#else  // Fortran Dim 3
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d0,
                            CARTWGTAVGOUTSIDEDOUB3D0)
            (ifirstc(0), ifirstc(1), ifirstc(2),
             ilastc(0), ilastc(1), ilastc(2),
             filo(0), filo(1), filo(2),
             fihi(0), fihi(1), fihi(2),
             cilo(0), cilo(1), cilo(2),
             cihi(0), cihi(1), cihi(2),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(0, side, d),
             cdata->getPointer(0, side, d));
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d1,
                            CARTWGTAVGOUTSIDEDOUB3D1)
            (ifirstc(0), ifirstc(1), ifirstc(2),
             ilastc(0), ilastc(1), ilastc(2),
             filo(0), filo(1), filo(2),
             fihi(0), fihi(1), fihi(2),
             cilo(0), cilo(1), cilo(2),
             cihi(0), cihi(1), cihi(2),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(1, side, d),
             cdata->getPointer(1, side, d));
            SAMRAI_F77_FUNC(cartwgtavgoutsidedoub3d2,
                            CARTWGTAVGOUTSIDEDOUB3D2)
            (ifirstc(0), ifirstc(1), ifirstc(2),
             ilastc(0), ilastc(1), ilastc(2),
             filo(0), filo(1), filo(2),
             fihi(0), fihi(1), fihi(2),
             cilo(0), cilo(1), cilo(2),
             cihi(0), cihi(1), cihi(2),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(2, side, d),
             cdata->getPointer(2, side, d));
#endif
         } else {
            TBOX_ERROR("CartesianOutersideDoubleWeightedAverage error...\n"
                       << "dim > 3 not supported." << std::endl);
         }
      }
   }
}

}  // namespace geom
}  // namespace SAMRAI
