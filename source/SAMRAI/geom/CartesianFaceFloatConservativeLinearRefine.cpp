/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for face-centered
 *                float data on a Cartesian mesh.
 *
 ************************************************************************/
#include "SAMRAI/geom/CartesianFaceFloatConservativeLinearRefine.h"
#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
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
void SAMRAI_F77_FUNC(cartclinreffaceflot1d, CARTCLINREFFACEFLOT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *);
// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinreffaceflot2d0, CARTCLINREFFACEFLOT2D0) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *, float *, float *);
void SAMRAI_F77_FUNC(cartclinreffaceflot2d1, CARTCLINREFFACEFLOT2D1) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *, float *, float *);
// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinreffaceflot3d0, CARTCLINREFFACEFLOT3D0) (const int&,
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
void SAMRAI_F77_FUNC(cartclinreffaceflot3d1, CARTCLINREFFACEFLOT3D1) (const int&,
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
void SAMRAI_F77_FUNC(cartclinreffaceflot3d2, CARTCLINREFFACEFLOT3D2) (const int&,
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

CartesianFaceFloatConservativeLinearRefine::
CartesianFaceFloatConservativeLinearRefine():
   hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianFaceFloatConservativeLinearRefine::~
CartesianFaceFloatConservativeLinearRefine()
{
}

int
CartesianFaceFloatConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianFaceFloatConservativeLinearRefine::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
CartesianFaceFloatConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(dim, coarse, ratio);

   std::shared_ptr<pdat::FaceData<float> > cdata(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<float>, hier::PatchData>(
         coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::FaceData<float> > fdata(
      SAMRAI_SHARED_PTR_CAST<pdat::FaceData<float>, hier::PatchData>(
         fine.getPatchData(dst_component)));

   const pdat::FaceOverlap* t_overlap =
      CPP_CAST<const pdat::FaceOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

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

   for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {

         const hier::Box& face_box = *b;
         TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, face_box);

         hier::Box fine_box(dim);
         for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
            fine_box.setLower(
               static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()), face_box.lower(i));
            fine_box.setUpper(
               static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()), face_box.upper(i));
         }

         fine_box.setUpper(axis, fine_box.upper(axis) - 1);

         const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
         const hier::Index& ifirstc = coarse_box.lower();
         const hier::Index& ilastc = coarse_box.upper();
         const hier::Index& ifirstf = fine_box.lower();
         const hier::Index& ilastf = fine_box.upper();

         const hier::IntVector tmp_ghosts(dim, 0);
         std::vector<float> diff0(cgbox.numberCells(0) + 2);
         pdat::FaceData<float> slope0(cgbox, 1, tmp_ghosts);

         for (int d = 0; d < fdata->getDepth(); ++d) {
            if ((dim == tbox::Dimension(1))) {
               SAMRAI_F77_FUNC(cartclinreffaceflot1d, CARTCLINREFFACEFLOT1D) (
                  ifirstc(0), ilastc(0),
                  ifirstf(0), ilastf(0),
                  cilo(0), cihi(0),
                  filo(0), fihi(0),
                  &ratio[0],
                  cgeom->getDx(),
                  fgeom->getDx(),
                  cdata->getPointer(0, d),
                  fdata->getPointer(0, d),
                  &diff0[0], slope0.getPointer(0));
            } else if ((dim == tbox::Dimension(2))) {
               std::vector<float> diff1(cgbox.numberCells(1) + 2);
               pdat::FaceData<float> slope1(cgbox, 1, tmp_ghosts);

               if (axis == 0) {
                  SAMRAI_F77_FUNC(cartclinreffaceflot2d0, CARTCLINREFFACEFLOT2D0) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cgeom->getDx(),
                     fgeom->getDx(),
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d),
                     &diff0[0], slope0.getPointer(0),
                     &diff1[0], slope1.getPointer(0));
               } else if (axis == 1) {
                  SAMRAI_F77_FUNC(cartclinreffaceflot2d1, CARTCLINREFFACEFLOT2D1) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cgeom->getDx(),
                     fgeom->getDx(),
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d),
                     &diff1[0], slope1.getPointer(1),
                     &diff0[0], slope0.getPointer(1));
               }
            } else if ((dim == tbox::Dimension(3))) {
               std::vector<float> diff1(cgbox.numberCells(1) + 2);
               pdat::FaceData<float> slope1(cgbox, 1, tmp_ghosts);

               std::vector<float> diff2(cgbox.numberCells(2) + 2);
               pdat::FaceData<float> slope2(cgbox, 1, tmp_ghosts);

               if (axis == 0) {
                  SAMRAI_F77_FUNC(cartclinreffaceflot3d0, CARTCLINREFFACEFLOT3D0) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d),
                     &diff0[0], slope0.getPointer(0),
                     &diff1[0], slope1.getPointer(0),
                     &diff2[0], slope2.getPointer(0));
               } else if (axis == 1) {
                  SAMRAI_F77_FUNC(cartclinreffaceflot3d1, CARTCLINREFFACEFLOT3D1) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d),
                     &diff1[0], slope1.getPointer(1),
                     &diff2[0], slope2.getPointer(1),
                     &diff0[0], slope0.getPointer(1));
               } else if (axis == 2) {
                  SAMRAI_F77_FUNC(cartclinreffaceflot3d2, CARTCLINREFFACEFLOT3D2) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(2, d),
                     fdata->getPointer(2, d),
                     &diff2[0], slope2.getPointer(2),
                     &diff0[0], slope0.getPointer(2),
                     &diff1[0], slope1.getPointer(2));
               }
            } else {
               TBOX_ERROR(
                  "CartesianFaceFloatConservativeLinearRefine error...\n"
                  << "dim > 3 not supported." << std::endl);
            }
         }
      }
   }
}

}
}
