/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered float
 *                data on a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianCellFloatConservativeLinearRefine.h"
#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
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

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in cartrefine1d.f:
void SAMRAI_F77_FUNC(cartclinrefcellflot1d, CARTCLINREFCELLFLOT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcellflot2d, CARTCLINREFCELLFLOT2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *, float *, float *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinrefcellflot3d, CARTCLINREFCELLFLOT3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *, float *,
   float *, float *, float *);
}

namespace SAMRAI {
namespace geom {


CartesianCellFloatConservativeLinearRefine::
CartesianCellFloatConservativeLinearRefine():
   hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianCellFloatConservativeLinearRefine::~
CartesianCellFloatConservativeLinearRefine()
{
}

int
CartesianCellFloatConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellFloatConservativeLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianCellFloatConservativeLinearRefine::refine(
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

void
CartesianCellFloatConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<pdat::CellData<float> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<float>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<float> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<float>, hier::PatchData>(
         fine.getPatchData(dst_component)));

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   const std::shared_ptr<CartesianPatchGeometry> cgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         coarse.getPatchGeometry()));
   const std::shared_ptr<CartesianPatchGeometry> fgeom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, hier::PatchGeometry>(
         fine.getPatchGeometry()));

   TBOX_ASSERT(cgeom);
   TBOX_ASSERT(fgeom);

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();
   const hier::Index& ifirstf = fine_box.lower();
   const hier::Index& ilastf = fine_box.upper();

   const hier::IntVector tmp_ghosts(dim, 0);
   std::vector<float> diff0(cgbox.numberCells(0) + 1);
   pdat::CellData<float> slope0(cgbox, 1, tmp_ghosts);

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartclinrefcellflot1d, CARTCLINREFCELLFLOT1D) (ifirstc(0),
            ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer());
      } else if ((dim == tbox::Dimension(2))) {
         std::vector<float> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<float> slope1(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcellflot2d, CARTCLINREFCELLFLOT2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer(),
            &diff1[0], slope1.getPointer());
      } else if ((dim == tbox::Dimension(3))) {
         std::vector<float> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<float> slope1(cgbox, 1, tmp_ghosts);

         std::vector<float> diff2(cgbox.numberCells(2) + 1);
         pdat::CellData<float> slope2(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcellflot3d, CARTCLINREFCELLFLOT3D) (ifirstc(0),
            ifirstc(1), ifirstc(2),
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
            cdata->getPointer(d),
            fdata->getPointer(d),
            &diff0[0], slope0.getPointer(),
            &diff1[0], slope1.getPointer(),
            &diff2[0], slope2.getPointer());
      } else {
         TBOX_ERROR("CartesianFloatConservativeLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);

      }
   }
}

}
}
