/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered
 *                double data on a Skeleton mesh.
 *
 ************************************************************************/

#include "SkeletonCellDoubleConservativeLinearRefine.h"

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
// in cartrefine1d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *, double *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
}

using namespace SAMRAI;

SkeletonCellDoubleConservativeLinearRefine::
SkeletonCellDoubleConservativeLinearRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator("SKELETON_CONSERVATIVE_LINEAR_REFINE"),
   d_dim(dim)
{
   const int max_levels = 10;
   d_dx.resize(max_levels);
   for (int n = 0; n < max_levels; ++n) {
      d_dx[n].resize(dim.getValue());
      for (int i = 0; i < dim.getValue(); ++i) {
         d_dx[n][i] = 1.;
      }
   }

}

SkeletonCellDoubleConservativeLinearRefine::~
SkeletonCellDoubleConservativeLinearRefine()
{
}

int
SkeletonCellDoubleConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
SkeletonCellDoubleConservativeLinearRefine::getStencilWidth(const tbox::Dimension& dim) const {
   return hier::IntVector(dim, 1);
}

void SkeletonCellDoubleConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      CPP_CAST<const pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      refine(fine,
         coarse,
         dst_component,
         src_component,
         *b,
         ratio);
   }
}

void SkeletonCellDoubleConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   std::shared_ptr<pdat::CellData<double> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         fine.getPatchData(dst_component)));

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   const hier::IntVector tmp_ghosts(fine.getDim(), 0);
   std::vector<double> diff0(cgbox.numberCells(0) + 1);
   pdat::CellData<double> slope0(cgbox, 1, tmp_ghosts);

   int flev_num = fine.getPatchLevelNumber();
   int clev_num = coarse.getPatchLevelNumber();

   // deal with levels not in hierarchy
   if (flev_num < 0) flev_num = clev_num + 1;
   if (clev_num < 0) clev_num = flev_num - 1;

   double cdx[SAMRAI::MAX_DIM_VAL];
   double fdx[SAMRAI::MAX_DIM_VAL];
   getDx(clev_num, cdx);
   getDx(flev_num, fdx);

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if (fine.getDim() == tbox::Dimension(1)) {
         SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D) (
            ifirstc(0), ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cdx,
            fdx,
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer());
      } else if (fine.getDim() == tbox::Dimension(2)) {

         std::vector<double> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<double> slope1(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D) (
            ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cdx,
            fdx,
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer(),
            &diff1[0], slope1.getPointer());
      } else if (fine.getDim() == tbox::Dimension(3)) {

         std::vector<double> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<double> slope1(cgbox, 1, tmp_ghosts);

         std::vector<double> diff2(cgbox.numberCells(2) + 1);
         pdat::CellData<double> slope2(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D) (
            ifirstc(0), ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            ifirstf(0), ifirstf(1), ifirstf(2),
            ilastf(0), ilastf(1), ilastf(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            &ratio[0],
            cdx,
            fdx,
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer(),
            &diff1[0], slope1.getPointer(),
            &diff2[0], slope2.getPointer());
      } else {
         TBOX_ERROR("SkeletonCellDoubleConservativeLinearRefine error...\n"
            << "dimension > 3 not supported." << std::endl);

      }
   }
}

void SkeletonCellDoubleConservativeLinearRefine::setDx(
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

void SkeletonCellDoubleConservativeLinearRefine::getDx(
   const int level_number,
   double* dx) const
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      dx[i] = d_dx[level_number][i];
   }
}
