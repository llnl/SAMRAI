/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for edge-centered double data on
 *                a  mesh.
 *
 ************************************************************************/
#ifndef included_pdat_EdgeConstantRefine_C
#define included_pdat_EdgeConstantRefine_C

#include "SAMRAI/pdat/EdgeConstantRefine.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */

namespace SAMRAI
{
namespace pdat
{

template <typename T>
int EdgeConstantRefine<T>::getOperatorPriority() const
{
   return 0;
}

template <typename T>
hier::IntVector
EdgeConstantRefine<T>::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

template <typename T>
void EdgeConstantRefine<T>::refine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const hier::BoxOverlap& fine_overlap,
    const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(fine.getDim());

   std::shared_ptr<EdgeData<T> > cdata(
       SAMRAI_SHARED_PTR_CAST<EdgeData<T>, hier::PatchData>(
           coarse.getPatchData(src_component)));
   std::shared_ptr<EdgeData<T> > fdata(
       SAMRAI_SHARED_PTR_CAST<EdgeData<T>, hier::PatchData>(
           fine.getPatchData(dst_component)));

   const EdgeOverlap* t_overlap = CPP_CAST<const EdgeOverlap*>(&fine_overlap);

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

   for (int axis = 0; axis < dim.getValue(); ++axis) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {

         hier::Box fine_box(*b);
         TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, fine_box);

         for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
            if (i != axis) {
               fine_box.setUpper(i, fine_box.upper(i) - 1);
            }
         }

         const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
         const hier::Index& ifirstc = coarse_box.lower();
         const hier::Index& ilastc = coarse_box.upper();
         const hier::Index& ifirstf = fine_box.lower();
         const hier::Index& ilastf = fine_box.upper();

         for (int d = 0; d < fdata->getDepth(); ++d) {
            if (dim == tbox::Dimension(1)) {
               Call1dFortranEdge<T>(
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

               if (axis == 0) {
                  fine_box_plus.growUpper(1, 1);
               } else if (axis == 1) {
                  fine_box_plus.growUpper(0, 1);
               }

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
                  Call2dFortranEdge_d0<T>(
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cdata->getPointer(0, d),
                      fdata->getPointer(0, d));
               } else if (axis == 1) {
                  Call2dFortranEdge_d1<T>(
                      ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                      ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                      cilo(0), cilo(1), cihi(0), cihi(1),
                      filo(0), filo(1), fihi(0), fihi(1),
                      &ratio[0],
                      cdata->getPointer(1, d),
                      fdata->getPointer(1, d));
               }
#endif  // test for RAJA
            } else if (dim == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
               SAMRAI::hier::Box fine_box_plus = fine_box;

               if (axis == 0) {
                  fine_box_plus.growUpper(1, 1);
                  fine_box_plus.growUpper(2, 1);
               } else if (axis == 1) {
                  fine_box_plus.growUpper(0, 1);
                  fine_box_plus.growUpper(2, 1);
               } else if (axis == 2) {
                  fine_box_plus.growUpper(0, 1);
                  fine_box_plus.growUpper(1, 1);
               }

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
                  Call3dFortranEdge_d0<T>(
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
                  Call3dFortranEdge_d1<T>(
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
                  Call3dFortranEdge_d2<T>(
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
                   "EdgeConstantRefine::refine dimension > 3 not supported"
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
