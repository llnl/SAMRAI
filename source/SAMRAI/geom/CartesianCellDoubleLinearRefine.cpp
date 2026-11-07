/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for cell-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"

#include <cfloat>
#include <cmath>
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
#pragma warning(disable : 1419)
#endif

// in cartrefine1d.f:
void SAMRAI_F77_FUNC(cartlinrefcelldoub1d, CARTLINREFCELLDOUB1D)(const int&,
                                                                 const int&,
                                                                 const int&, const int&,
                                                                 const int&, const int&,
                                                                 const int&, const int&,
                                                                 const int*, const double*, const double*,
                                                                 const double*, double*);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartlinrefcelldoub2d, CARTLINREFCELLDOUB2D)(const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&, const int&,
                                                                 const int&, const int&, const int&, const int&,
                                                                 const int&, const int&, const int&, const int&,
                                                                 const int*, const double*, const double*,
                                                                 const double*, double*);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartlinrefcelldoub3d, CARTLINREFCELLDOUB3D)(const int&,
                                                                 const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int&, const int&, const int&,
                                                                 const int*, const double*, const double*,
                                                                 const double*, double*);
}

namespace SAMRAI
{
namespace geom
{


CartesianCellDoubleLinearRefine::CartesianCellDoubleLinearRefine() : hier::RefineOperator("LINEAR_REFINE")
{
}

CartesianCellDoubleLinearRefine::~CartesianCellDoubleLinearRefine()
{
}

int CartesianCellDoubleLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellDoubleLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void CartesianCellDoubleLinearRefine::refine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const hier::BoxOverlap& fine_overlap,
    const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
       CPP_CAST<const pdat::CellOverlap*>(&fine_overlap);

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

void CartesianCellDoubleLinearRefine::refine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const hier::Box& fine_box,
    const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, coarse, fine_box, ratio);

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

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if ((dim == tbox::Dimension(1))) {
         SAMRAI_F77_FUNC(cartlinrefcelldoub1d, CARTLINREFCELLDOUB1D)
         (ifirstc(0),
          ilastc(0),
          ifirstf(0), ilastf(0),
          cilo(0), cihi(0),
          filo(0), fihi(0),
          &ratio[0],
          cgeom->getDx(),
          fgeom->getDx(),
          cdata->getPointer(d),
          fdata->getPointer(d));
      } else if ((dim == tbox::Dimension(2))) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<2>(d);
         auto coarse_array = cdata->getView<2>(d);
         const double* fdx = fgeom->getDx();
         const double* cdx = cgeom->getDx();
         const double fdx0 = fdx[0];
         const double fdx1 = fdx[1];
         const double cdx0 = cdx[0];
         const double cdx1 = cdx[1];

         const int r0 = ratio[0];
         const int r1 = ratio[1];

         hier::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE(int j /*fast*/, int k) {
            const int ic0 = (j < 0) ? (j + 1) / r0 - 1 : j / r0;
            const int ic1 = (k < 0) ? (k + 1) / r1 - 1 : k / r1;
            const int ir0 = j - ic0 * r0;
            const int ir1 = k - ic1 * r1;

            int jj = ic0;
            int kk = ic1;

            const double deltax0 = (static_cast<double>(ir0) + 0.5) * fdx0 - cdx0 * 0.5;
            const double deltax1 = (static_cast<double>(ir1) + 0.5) * fdx1 - cdx1 * 0.5;

            double x = deltax0 / cdx0;
            double y = deltax1 / cdx1;

            if (x < 0.0) {
               jj--;
               x += 1.0;
            }
            if (y < 0.0) {
               kk--;
               y += 1.0;
            }
            fine_array(j, k) = (coarse_array(jj, kk) + (coarse_array(jj + 1, kk) - coarse_array(jj, kk)) * x) * (1.0 - y) + (coarse_array(jj, kk + 1) + (coarse_array(jj + 1, kk + 1) - coarse_array(jj, kk + 1)) * x) * y;
         });
#else   // Fortran Dimension 2

         SAMRAI_F77_FUNC(cartlinrefcelldoub2d, CARTLINREFCELLDOUB2D)
         (ifirstc(0),
          ifirstc(1), ilastc(0), ilastc(1),
          ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
          cilo(0), cilo(1), cihi(0), cihi(1),
          filo(0), filo(1), fihi(0), fihi(1),
          &ratio[0],
          cgeom->getDx(),
          fgeom->getDx(),
          cdata->getPointer(d),
          fdata->getPointer(d));
#endif  // test for RAJA
      } else if ((dim == tbox::Dimension(3))) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<3>(d);
         auto coarse_array = cdata->getView<3>(d);
         const double* fdx = fgeom->getDx();
         const double* cdx = cgeom->getDx();
         const double fdx0 = fdx[0];
         const double fdx1 = fdx[1];
         const double fdx2 = fdx[2];
         const double cdx0 = cdx[0];
         const double cdx1 = cdx[1];
         const double cdx2 = cdx[2];

         const int r0 = ratio[0];
         const int r1 = ratio[1];
         const int r2 = ratio[2];

         hier::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE(int i /*fastest */, int j, int k) {
            const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;
            const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
            const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;

            const int ir0 = i - ic0 * r0;
            const int ir1 = j - ic1 * r1;
            const int ir2 = k - ic2 * r2;

            int ii = ic0;
            int jj = ic1;
            int kk = ic2;

            const double deltax0 = (static_cast<double>(ir0) + 0.5) * fdx0 - cdx0 * 0.5;
            const double deltax1 = (static_cast<double>(ir1) + 0.5) * fdx1 - cdx1 * 0.5;
            const double deltax2 = (static_cast<double>(ir2) + 0.5) * fdx2 - cdx2 * 0.5;

            double x = deltax0 / cdx0;
            double y = deltax1 / cdx1;
            double z = deltax2 / cdx2;

            if (x < 0.0) {
               ii--;
               x += 1.0;
            }
            if (y < 0.0) {
               jj--;
               y += 1.0;
            }
            if (z < 0.0) {
               kk--;
               z += 1.0;
            }

            fine_array(i, j, k) = ((coarse_array(ii, jj, kk) + (coarse_array(ii + 1, jj, kk) - coarse_array(ii, jj, kk)) * x) * (1.0 - y) + (coarse_array(ii, jj + 1, kk) + (coarse_array(ii + 1, jj + 1, kk) - coarse_array(ii, jj + 1, kk)) * x) * y) * (1.0 - z) + ((coarse_array(ii, jj, kk + 1) + (coarse_array(ii + 1, jj, kk + 1) - coarse_array(ii, jj, kk + 1)) * x) * (1.0 - y) + (coarse_array(ii, jj + 1, kk + 1) + (coarse_array(ii + 1, jj + 1, kk + 1) - coarse_array(ii, jj + 1, kk + 1)) * x) * y) * z;
         });

#else
         SAMRAI_F77_FUNC(cartlinrefcelldoub3d, CARTLINREFCELLDOUB3D)
         (ifirstc(0),
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
          fdata->getPointer(d));
#endif
      } else {
         TBOX_ERROR("CartesianCellDoubleLinearRefine error...\n"
                    << "dim > 3 not supported." << std::endl);
      }
   }
}

}  // namespace geom
}  // namespace SAMRAI
