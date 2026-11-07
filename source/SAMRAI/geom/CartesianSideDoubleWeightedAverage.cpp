/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for side-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"

#include <cfloat>
#include <cmath>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
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
void SAMRAI_F77_FUNC(cartwgtavgsidedoub1d, CARTWGTAVGSIDEDOUB1D)(const int &,
                                                                 const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const int *, const double *, const double *,
                                                                 const double *, double *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgsidedoub2d0, CARTWGTAVGSIDEDOUB2D0)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);

void SAMRAI_F77_FUNC(cartwgtavgsidedoub2d1, CARTWGTAVGSIDEDOUB2D1)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgsidedoub3d0, CARTWGTAVGSIDEDOUB3D0)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgsidedoub3d1, CARTWGTAVGSIDEDOUB3D1)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *);
void SAMRAI_F77_FUNC(cartwgtavgsidedoub3d2, CARTWGTAVGSIDEDOUB3D2)(const int &,
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


CartesianSideDoubleWeightedAverage::CartesianSideDoubleWeightedAverage() : hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianSideDoubleWeightedAverage::~CartesianSideDoubleWeightedAverage()
{
}

int CartesianSideDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianSideDoubleWeightedAverage::getStencilWidth(const tbox::Dimension &dim) const
{
   return hier::IntVector::getZero(dim);
}

void CartesianSideDoubleWeightedAverage::coarsen(
    hier::Patch &coarse,
    const hier::Patch &fine,
    const int dst_component,
    const int src_component,
    const hier::Box &coarse_box,
    const hier::IntVector &ratio) const
{
   const tbox::Dimension &dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::SideData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::SideData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           coarse.getPatchData(dst_component)));

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::IntVector &directions(cdata->getDirectionVector());

   TBOX_ASSERT(directions ==
               hier::IntVector::min(directions, fdata->getDirectionVector()));

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
         if (directions(0)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub1d, CARTWGTAVGSIDEDOUB1D)
            (ifirstc(0),
             ilastc(0),
             filo(0), fihi(0),
             cilo(0), cihi(0),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(0, d),
             cdata->getPointer(0, d));
         }
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


         if (directions(0)) {
            SAMRAI::hier::Box coarse_box_plus = coarse_box;
            coarse_box_plus.growUpper(0, 1);

            auto fine_array = fdata->getConstView<2>(0, d);
            auto coarse_array = cdata->getView<2>(0, d);

            const double lengthf = fdx1;
            const double lengthc = cdx1;

            hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
               double spv = 0.0;
               int jf = j * r0;
               for (int ry = 0; ry < r1; ry++) {
                  int kf = k * r1 + ry;
                  spv += fine_array(jf, kf) * lengthf;
               }

               coarse_array(j, k) = spv / lengthc;
            });
         }
         if (directions(1)) {
            SAMRAI::hier::Box coarse_box_plus = coarse_box;
            coarse_box_plus.growUpper(1, 1);

            auto fine_array = fdata->getConstView<2>(1, d);
            auto coarse_array = cdata->getView<2>(1, d);

            const double lengthf = fdx0;
            const double lengthc = cdx0;

            hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
               double spv = 0.0;
               int kf = k * r1;
               for (int rx = 0; rx < r0; rx++) {
                  int jf = j * r0 + rx;
                  spv += fine_array(jf, kf) * lengthf;
               }

               coarse_array(j, k) = spv / lengthc;
            });
         }

#else
         if (directions(0)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub2d0, CARTWGTAVGSIDEDOUB2D0)
            (ifirstc(0),
             ifirstc(1), ilastc(0), ilastc(1),
             filo(0), filo(1), fihi(0), fihi(1),
             cilo(0), cilo(1), cihi(0), cihi(1),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(0, d),
             cdata->getPointer(0, d));
         }
         if (directions(1)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub2d1, CARTWGTAVGSIDEDOUB2D1)
            (ifirstc(0),
             ifirstc(1), ilastc(0), ilastc(1),
             filo(0), filo(1), fihi(0), fihi(1),
             cilo(0), cilo(1), cihi(0), cihi(1),
             &ratio[0],
             fgeom->getDx(),
             cgeom->getDx(),
             fdata->getPointer(1, d),
             cdata->getPointer(1, d));
         }
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

         if (directions(0)) {
            SAMRAI::hier::Box coarse_box_plus = coarse_box;
            coarse_box_plus.growUpper(0, 1);

            auto fine_array = fdata->getConstView<3>(0, d);
            auto coarse_array = cdata->getView<3>(0, d);

            const double areaf = fdx1 * fdx2;
            const double areac = cdx1 * cdx2;

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
         }
         if (directions(1)) {
            SAMRAI::hier::Box coarse_box_plus = coarse_box;
            coarse_box_plus.growUpper(1, 1);

            auto fine_array = fdata->getConstView<3>(1, d);
            auto coarse_array = cdata->getView<3>(1, d);

            const double areaf = fdx0 * fdx2;
            const double areac = cdx0 * cdx2;

            hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
               double spv = 0.0;
               int jj = j * r1;
               for (int rz = 0; rz < r2; rz++) {
                  for (int rx = 0; rx < r0; rx++) {
                     int ii = i * r0 + rx;
                     int kk = k * r2 + rz;
                     spv += fine_array(ii, jj, kk) * areaf;
                  }
               }

               coarse_array(i, j, k) = spv / areac;
            });
         }
         if (directions(2)) {
            SAMRAI::hier::Box coarse_box_plus = coarse_box;
            coarse_box_plus.growUpper(2, 1);

            auto fine_array = fdata->getConstView<3>(2, d);
            auto coarse_array = cdata->getView<3>(2, d);

            const double areaf = fdx0 * fdx1;
            const double areac = cdx0 * cdx1;

            hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
               double spv = 0.0;
               int kk = k * r2;
               for (int ry = 0; ry < r1; ry++) {
                  for (int rx = 0; rx < r0; rx++) {
                     int ii = i * r0 + rx;
                     int jj = j * r1 + ry;
                     spv += fine_array(ii, jj, kk) * areaf;
                  }
               }

               coarse_array(i, j, k) = spv / areac;
            });
         }

#else
         if (directions(0)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub3d0, CARTWGTAVGSIDEDOUB3D0)
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
         }
         if (directions(1)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub3d1, CARTWGTAVGSIDEDOUB3D1)
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
         }
         if (directions(2)) {
            SAMRAI_F77_FUNC(cartwgtavgsidedoub3d2, CARTWGTAVGSIDEDOUB3D2)
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
         }
#endif
      } else {
         TBOX_ERROR("CartesianCellSideDoubleWeightedAverage error...\n"
                    << "dim > 3 not supported." << std::endl);
      }
   }
}

}  // namespace geom
}  // namespace SAMRAI
