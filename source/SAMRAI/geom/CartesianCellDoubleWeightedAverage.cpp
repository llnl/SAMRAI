/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/ForAll.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

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
void SAMRAI_F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D)(const int &,
                                                                 const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &,
                                                                 const int *, const double *, const double *,
                                                                 const double *, double *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D)(const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int &, const int &, const int &, const int &,
                                                                 const int &, const int &, const int &, const int &,
                                                                 const int *, const double *, const double *,
                                                                 const double *, double *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D)(const int &,
                                                                 const int &, const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int &, const int &, const int &,
                                                                 const int *, const double *, const double *,
                                                                 const double *, double *);
// in cartcoarsen4d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub4d, CARTWGTAVGCELLDOUB4D)(const int &,
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

CartesianCellDoubleWeightedAverage::CartesianCellDoubleWeightedAverage() : hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianCellDoubleWeightedAverage::~CartesianCellDoubleWeightedAverage()
{
}

int CartesianCellDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleWeightedAverage::getStencilWidth(const tbox::Dimension &dim) const
{
   return hier::IntVector::getZero(dim);
}

void CartesianCellDoubleWeightedAverage::coarsen(
    hier::Patch &coarse,
    const hier::Patch &fine,
    const int dst_component,
    const int src_component,
    const hier::Box &coarse_box,
    const hier::IntVector &ratio) const
{
   RANGE_PUSH("WeightedAverage::coarsen", 4);

   const tbox::Dimension &dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::CellData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
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

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Index &ifirstc = coarse_box.lower();
   const hier::Index &ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D)
         (ifirstc(0),
          ilastc(0),
          filo(0), fihi(0),
          cilo(0), cihi(0),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(d),
          cdata->getPointer(d));
      } else if ((dim == tbox::Dimension(2))) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<2>(d);
         auto coarse_array = cdata->getView<2>(d);

         const double *fdx = fgeom->getDx();
         const double *cdx = cgeom->getDx();

         const double fdx0 = fdx[0];
         const double fdx1 = fdx[1];
         const double cdx0 = cdx[0];
         const double cdx1 = cdx[1];

         const int r0 = ratio[0];
         const int r1 = ratio[1];

         const double dVf = fdx0 * fdx1;
         const double dVc = cdx0 * cdx1;

         hier::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
            double spv = 0.0;

            for (int rx = 0; rx < r0; rx++) {
               for (int ry = 0; ry < r1; ry++) {
                  const int jf = j * r0 + rx;
                  const int kf = k * r1 + ry;
                  spv += fine_array(jf, kf) * dVf;
               }
            }

            coarse_array(j, k) = spv / dVc;
         });
#else
         SAMRAI_F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D)
         (ifirstc(0),
          ifirstc(1), ilastc(0), ilastc(1),
          filo(0), filo(1), fihi(0), fihi(1),
          cilo(0), cilo(1), cihi(0), cihi(1),
          &ratio[0],
          fgeom->getDx(),
          cgeom->getDx(),
          fdata->getPointer(d),
          cdata->getPointer(d));
#endif

      } else if ((dim == tbox::Dimension(3))) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<3>(d);
         auto coarse_array = cdata->getView<3>(d);

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

         const double dVf = fdx0 * fdx1 * fdx2;
         const double dVc = cdx0 * cdx1 * cdx2;

         hier::parallel_for_all(coarse_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
            double spv = 0.0;

            for (int rx = 0; rx < r0; rx++) {
               for (int ry = 0; ry < r1; ry++) {
                  for (int rz = 0; rz < r2; rz++) {
                     const int ii = i * r0 + rx;
                     const int jj = j * r1 + ry;
                     const int kk = k * r2 + rz;
                     spv += fine_array(ii, jj, kk) * dVf;
                  }
               }
            }

            coarse_array(i, j, k) = spv / dVc;
         });
#else
         SAMRAI_F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D)
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
          fdata->getPointer(d),
          cdata->getPointer(d));
#endif
      } else if ((dim == tbox::Dimension(4))) {
         SAMRAI_F77_FUNC(cartwgtavgcelldoub4d, CARTWGTAVGCELLDOUB4D)
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
          fdata->getPointer(d),
          cdata->getPointer(d));
      } else {
         TBOX_ERROR("CartesianCellDoubleWeightedAverage error...\n"
                    << "dim > 4 not supported." << std::endl);
      }
   }
   RANGE_POP
}

}  // namespace geom
}  // namespace SAMRAI
