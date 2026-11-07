/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for edge-centered float data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianEdgeFloatWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"
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
void SAMRAI_F77_FUNC(cartwgtavgedgeflot1d, CARTWGTAVGEDGEFLOT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgedgeflot2d0, CARTWGTAVGEDGEFLOT2D0) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);

void SAMRAI_F77_FUNC(cartwgtavgedgeflot2d1, CARTWGTAVGEDGEFLOT2D1) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgedgeflot3d0, CARTWGTAVGEDGEFLOT3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);
void SAMRAI_F77_FUNC(cartwgtavgedgeflot3d1, CARTWGTAVGEDGEFLOT3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);
void SAMRAI_F77_FUNC(cartwgtavgedgeflot3d2, CARTWGTAVGEDGEFLOT3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *);
}

namespace SAMRAI {
namespace geom {


CartesianEdgeFloatWeightedAverage::CartesianEdgeFloatWeightedAverage():
   hier::CoarsenOperator("CONSERVATIVE_COARSEN")
{
}

CartesianEdgeFloatWeightedAverage::~CartesianEdgeFloatWeightedAverage()
{
}

int
CartesianEdgeFloatWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianEdgeFloatWeightedAverage::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void
CartesianEdgeFloatWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, coarse_box, ratio);

   std::shared_ptr<pdat::EdgeData<float> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<float>, hier::PatchData>(
         fine.getPatchData(src_component)));
   std::shared_ptr<pdat::EdgeData<float> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::EdgeData<float>, hier::PatchData>(
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
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartwgtavgedgeflot1d, CARTWGTAVGEDGEFLOT1D) (ifirstc(0),
            ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(0, d),
            cdata->getPointer(0, d));
      } else if ((dim == tbox::Dimension(2))) {
         SAMRAI_F77_FUNC(cartwgtavgedgeflot2d0, CARTWGTAVGEDGEFLOT2D0) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(0, d),
            cdata->getPointer(0, d));
         SAMRAI_F77_FUNC(cartwgtavgedgeflot2d1, CARTWGTAVGEDGEFLOT2D1) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(1, d),
            cdata->getPointer(1, d));
      } else if ((dim == tbox::Dimension(3))) {
         SAMRAI_F77_FUNC(cartwgtavgedgeflot3d0, CARTWGTAVGEDGEFLOT3D0) (ifirstc(0),
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
         SAMRAI_F77_FUNC(cartwgtavgedgeflot3d1, CARTWGTAVGEDGEFLOT3D1) (ifirstc(0),
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
         SAMRAI_F77_FUNC(cartwgtavgedgeflot3d2, CARTWGTAVGEDGEFLOT3D2) (ifirstc(0),
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
      } else {
         TBOX_ERROR("CartesianEdgeFloatWeightedAverage error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
