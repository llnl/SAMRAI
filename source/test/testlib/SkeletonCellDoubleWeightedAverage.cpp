/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Skeleton mesh.
 *
 ************************************************************************/

#include "SkeletonCellDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {
// in cartcoarsen1d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void SAMRAI_F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
}

using namespace SAMRAI;

SkeletonCellDoubleWeightedAverage::SkeletonCellDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator("SKELETON_CONSERVATIVE_COARSEN"),
   d_dim(dim)
{
}

SkeletonCellDoubleWeightedAverage::~SkeletonCellDoubleWeightedAverage()
{
}

int SkeletonCellDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
SkeletonCellDoubleWeightedAverage::getStencilWidth(const tbox::Dimension& dim) const {
   return hier::IntVector(dim, 0);
}

void SkeletonCellDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   std::shared_ptr<pdat::CellData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(dst_component)));
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   int flev_num = fine.getPatchLevelNumber();
   int clev_num = coarse.getPatchLevelNumber();

   // deal with levels not in hierarchy
   if (flev_num < 0) flev_num = clev_num + 1;
   if (clev_num < 0) clev_num = flev_num - 1;

   double cdx[SAMRAI::MAX_DIM_VAL];
   double fdx[SAMRAI::MAX_DIM_VAL];
   getDx(clev_num, cdx);
   getDx(flev_num, fdx);

   for (int d = 0; d < cdata->getDepth(); ++d) {
      if (fine.getDim() == tbox::Dimension(1)) {
         SAMRAI_F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (
            ifirstc(0), ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (
            ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (
            ifirstc(0), ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR("SkeletonCellDoubleWeightedAverage error...\n"
            << "dimension > 3 not supported." << std::endl);

      }
   }
}

void SkeletonCellDoubleWeightedAverage::setDx(
   const int level_number,
   const double* dx)
{
   if (level_number >= static_cast<int>(d_dx.size())) {
      d_dx.resize(level_number + 1);
      d_dx[level_number].resize(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); ++i) {
         d_dx[level_number][i] = dx[i];
      }
   }
}

void SkeletonCellDoubleWeightedAverage::getDx(
   const int level_number,
   double* dx) const
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      dx[i] = d_dx[level_number][i];
   }
}
