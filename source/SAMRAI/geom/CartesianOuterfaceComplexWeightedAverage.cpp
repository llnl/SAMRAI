/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for outerface complex data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianOuterfaceComplexWeightedAverage.h"
#include "SAMRAI/tbox/Complex.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
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
#pragma warning (disable:1419)
#endif

// in cartcoarsen1d.f:
void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx1d, CARTWGTAVGOUTFACECPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx2d0, CARTWGTAVGOUTFACECPLX2D0) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);

void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx2d1, CARTWGTAVGOUTFACECPLX2D1) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d0, CARTWGTAVGOUTFACECPLX3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d1, CARTWGTAVGOUTFACECPLX3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d2, CARTWGTAVGOUTFACECPLX3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace geom {


CartesianOuterfaceComplexWeightedAverage::
CartesianOuterfaceComplexWeightedAverage():
   hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianOuterfaceComplexWeightedAverage::~
CartesianOuterfaceComplexWeightedAverage()
{
}

int
CartesianOuterfaceComplexWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianOuterfaceComplexWeightedAverage::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void
CartesianOuterfaceComplexWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::OuterfaceData<dcomplex> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<dcomplex>, hier::PatchData>(
         fine.getPatchData(src_component)));
   std::shared_ptr<pdat::OuterfaceData<dcomplex> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<dcomplex>, hier::PatchData>(
         coarse.getPatchData(dst_component)));
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();
   const hier::Index& cilo = cdata->getGhostBox().lower();
   const hier::Index& cihi = cdata->getGhostBox().upper();

   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));

   TBOX_ASSERT(fgeom);
   TBOX_ASSERT(cgeom);

   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); ++d) {
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; ++i) {
         if ((dim == tbox::Dimension(1))) {
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx1d,
               CARTWGTAVGOUTFACECPLX1D) (ifirstc(0), ilastc(0),
               filo(0), fihi(0),
               cilo(0), cihi(0),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
         } else if ((dim == tbox::Dimension(2))) {
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx2d0,
               CARTWGTAVGOUTFACECPLX2D0) (ifirstc(0), ifirstc(1), ilastc(0),
               ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx2d1,
               CARTWGTAVGOUTFACECPLX2D1) (ifirstc(0), ifirstc(1), ilastc(0),
               ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, i, d),
               cdata->getPointer(1, i, d));
         } else if ((dim == tbox::Dimension(3))) {
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d0,
               CARTWGTAVGOUTFACECPLX3D0) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d1,
               CARTWGTAVGOUTFACECPLX3D1) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, i, d),
               cdata->getPointer(1, i, d));
            SAMRAI_F77_FUNC(cartwgtavgoutfacecplx3d2,
               CARTWGTAVGOUTFACECPLX3D2) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(2, i, d),
               cdata->getPointer(2, i, d));
         } else {
            TBOX_ERROR("CartesianOuterfaceComplexWeightedAverage error...\n"
               << "dim > 3 not supported." << std::endl);
         }
      }
   }
}

}
}
