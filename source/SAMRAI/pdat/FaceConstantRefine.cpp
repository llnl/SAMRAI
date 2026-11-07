/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for face-centered double data on
 *                a  mesh.
 *
 ************************************************************************/
#ifndef included_pdat_FaceConstantRefine_C
#define included_pdat_FaceConstantRefine_C

#include "SAMRAI/pdat/FaceConstantRefine.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"


namespace SAMRAI
{
namespace pdat
{

template <typename T>
void FaceConstantRefine<T>::refine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const hier::BoxOverlap& fine_overlap,
    const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());

   std::shared_ptr<FaceData<T> > cdata(
       SAMRAI_SHARED_PTR_CAST<FaceData<T>, hier::PatchData>(
           coarse.getPatchData(src_component)));
   std::shared_ptr<FaceData<T> > fdata(
       SAMRAI_SHARED_PTR_CAST<FaceData<T>, hier::PatchData>(
           fine.getPatchData(dst_component)));

   const FaceOverlap* t_overlap =
       CPP_CAST<const FaceOverlap*>(&fine_overlap);

   TBOX_ASSERT(t_overlap != 0);

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_ASSERT_OBJDIM_EQUALITY3(fine, coarse, ratio);

   const hier::Box& cgbox(cdata->getGhostBox());

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {

         const hier::Box& face_box = *b;
         TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, face_box);

         hier::Box fine_box(dim);
         for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
            fine_box.setLower(
                static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()),
                face_box.lower(i));
            fine_box.setUpper(
                static_cast<tbox::Dimension::dir_t>((axis + i) % dim.getValue()),
                face_box.upper(i));
         }

         fine_box.setUpper(axis, fine_box.upper(axis) - 1);

         const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
         const hier::Index& ifirstc = coarse_box.lower();
         const hier::Index& ilastc = coarse_box.upper();
         const hier::Index& ifirstf = fine_box.lower();
         const hier::Index& ilastf = fine_box.upper();

         for (int d = 0; d < fdata->getDepth(); ++d) {
            if (dim == tbox::Dimension(1)) {
               Call1dFortranFace(
                   ifirstc(0), ilastc(0),
                   ifirstf(0), ilastf(0),
                   cilo(0), cihi(0),
                   filo(0), fihi(0),
                   &ratio[0],
                   cdata->getPointer(0, d),
                   fdata->getPointer(0, d));
            } else if (dim == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
               SAMRAI::hier::Box fine_box_plus = fine_box;

               if (axis == 1) {  // transpose boxes
                  fine_box_plus.setLower(0, fine_box.lower(1));
                  fine_box_plus.setLower(1, fine_box.lower(0));
                  fine_box_plus.setUpper(0, fine_box.upper(1));
                  fine_box_plus.setUpper(1, fine_box.upper(0));
               }
               fine_box_plus.growUpper(0, 1);

               auto fine_array = fdata->template getView<2>(axis, d);
               auto coarse_array = cdata->template getConstView<2>(axis, d);

               const int r0 = ratio[0];
               const int r1 = ratio[1];

               hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int j, int k) {
                  const int ic1 = (k < 0) ? (k + 1) / r1 - 1 : k / r1;
                  const int ic0 = (j < 0) ? (j + 1) / r0 - 1 : j / r0;

                  fine_array(j, k) = coarse_array(ic0, ic1);
               });
#else   // Fortran Dimension 2
               if (axis == 0) {
                  Call2dFortranFace_d0(
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cdata->getPointer(0, d),
                      fdata->getPointer(0, d));
               } else if (axis == 1) {
                  Call2dFortranFace_d1(
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cdata->getPointer(1, d),
                      fdata->getPointer(1, d));
               }
#endif  // Test for RAJA
            } else if (dim == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
               SAMRAI::hier::Box fine_box_plus = fine_box;

               if (axis == 1) {  // transpose boxes <1,2,0>
                  fine_box_plus.setLower(0, fine_box.lower(1));
                  fine_box_plus.setLower(1, fine_box.lower(2));
                  fine_box_plus.setLower(2, fine_box.lower(0));

                  fine_box_plus.setUpper(0, fine_box.upper(1));
                  fine_box_plus.setUpper(1, fine_box.upper(2));
                  fine_box_plus.setUpper(2, fine_box.upper(0));
               } else if (axis == 2) {  // <2,0,1>
                  fine_box_plus.setLower(0, fine_box.lower(2));
                  fine_box_plus.setLower(1, fine_box.lower(0));
                  fine_box_plus.setLower(2, fine_box.lower(1));

                  fine_box_plus.setUpper(0, fine_box.upper(2));
                  fine_box_plus.setUpper(1, fine_box.upper(0));
                  fine_box_plus.setUpper(2, fine_box.upper(1));
               }
               fine_box_plus.growUpper(0, 1);

               auto fine_array = fdata->template getView<3>(axis, d);
               auto coarse_array = cdata->template getConstView<3>(axis, d);

               const int r0 = ratio[0];
               const int r1 = ratio[1];
               const int r2 = ratio[2];

               hier::parallel_for_all(fine_box_plus, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
                  const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;
                  const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
                  const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;

                  fine_array(i, j, k) = coarse_array(ic0, ic1, ic2);
               });

#else   // Fortran Dimension 3
               if (axis == 0) {
                  Call3dFortranFace_d0(
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cdata->getPointer(0, d),
                      fdata->getPointer(0, d));
               } else if (axis == 1) {
                  Call3dFortranFace_d1(
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cdata->getPointer(1, d),
                      fdata->getPointer(1, d));
               } else if (axis == 2) {
                  Call3dFortranFace_d2(
                      ifirstc(0), ifirstc(1), ifirstc(2),
                      ilastc(0), ilastc(1), ilastc(2),
                      ifirstf(0), ifirstf(1), ifirstf(2),
                      ilastf(0), ilastf(1), ilastf(2),
                      cilo(0), cilo(1), cilo(2),
                      cihi(0), cihi(1), cihi(2),
                      filo(0), filo(1), filo(2),
                      fihi(0), fihi(1), fihi(2),
                      &ratio[0],
                      cdata->getPointer(2, d),
                      fdata->getPointer(2, d));
               }
#endif  // Test for RAJA
            } else {
               TBOX_ERROR(
                   "FaceConstantRefine::refine() dimension > 3 not supported"
                   << std::endl);
            }
         }
      }
   }
}

}  // namespace pdat
}  // namespace SAMRAI


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
