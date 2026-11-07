/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered
 *                omplex data on a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianCellComplexConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cfloat>
#include <cmath>

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
void SAMRAI_F77_FUNC(cartclinrefcellcplx1d, CARTCLINREFCELLCPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *,
   dcomplex *, dcomplex *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcellcplx2d, CARTCLINREFCELLCPLX2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *,
   dcomplex *, dcomplex *, dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(cartclinrefcellcplx3d, CARTCLINREFCELLCPLX3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *,
   dcomplex *, dcomplex *, dcomplex *,
   dcomplex *, dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace geom {


CartesianCellComplexConservativeLinearRefine::
CartesianCellComplexConservativeLinearRefine():
   hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianCellComplexConservativeLinearRefine::~
CartesianCellComplexConservativeLinearRefine()
{
}

int
CartesianCellComplexConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellComplexConservativeLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianCellComplexConservativeLinearRefine::refine(
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
CartesianCellComplexConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

   std::shared_ptr<pdat::CellData<dcomplex> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<dcomplex> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
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
   std::vector<dcomplex> diff0(cgbox.numberCells(0) + 1);
   pdat::CellData<dcomplex> slope0(cgbox, 1, tmp_ghosts);

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartclinrefcellcplx1d, CARTCLINREFCELLCPLX1D) (ifirstc(0),
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
         std::vector<dcomplex> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<dcomplex> slope1(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcellcplx2d, CARTCLINREFCELLCPLX2D) (ifirstc(0),
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
         std::vector<dcomplex> diff1(cgbox.numberCells(1) + 1);
         pdat::CellData<dcomplex> slope1(cgbox, 1, tmp_ghosts);

         std::vector<dcomplex> diff2(cgbox.numberCells(2) + 1);
         pdat::CellData<dcomplex> slope2(cgbox, 1, tmp_ghosts);

         SAMRAI_F77_FUNC(cartclinrefcellcplx3d, CARTCLINREFCELLCPLX3D) (ifirstc(0),
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
         TBOX_ERROR("CartesianCellComplexConservativeLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
