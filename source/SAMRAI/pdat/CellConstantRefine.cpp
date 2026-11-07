/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_CellConstantRefine_C
#define included_pdat_CellConstantRefine_C

#include "SAMRAI/pdat/CellConstantRefine.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>
#include <complex>


namespace SAMRAI
{
namespace pdat
{

template <typename T>
void CellConstantRefine<T>::refine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const hier::Box& fine_box,
    const hier::IntVector& ratio) const
{
   std::shared_ptr<CellData<T> > cdata(
       SAMRAI_SHARED_PTR_CAST<CellData<T>, hier::PatchData>(
           coarse.getPatchData(src_component)));
   std::shared_ptr<CellData<T> > fdata(
       SAMRAI_SHARED_PTR_CAST<CellData<T>, hier::PatchData>(
           fine.getPatchData(dst_component)));

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_ASSERT_OBJDIM_EQUALITY4(fine, coarse, fine_box, ratio);

   const hier::Box& cgbox(cdata->getGhostBox());

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();
   const hier::Index& ifirstf = fine_box.lower();
   const hier::Index& ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); ++d) {
      if (fine.getDim() == tbox::Dimension(1)) {
         Call1dFortranCell<T>(ifirstc(0), ilastc(0),
                              ifirstf(0), ilastf(0),
                              cilo(0), cihi(0),
                              filo(0), fihi(0),
                              &ratio[0],
                              cdata->getPointer(d),
                              fdata->getPointer(d));

      } else if (fine.getDim() == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->template getView<2>(d);
         auto coarse_array = cdata->template getView<2>(d);
         const int r0 = ratio[0];
         const int r1 = ratio[1];

         hier::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE(int j, int k) {
            const int ic1 = (k < 0) ? (k + 1) / r1 - 1 : k / r1;
            const int ic0 = (j < 0) ? (j + 1) / r0 - 1 : j / r0;

            fine_array(j, k) = coarse_array(ic0, ic1);
         });
#else  // Fortran Dimension 2
         Call2dFortranCell<T>(ifirstc(0), ifirstc(1),
                              ilastc(0), ilastc(1),
                              ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                              cilo(0), cilo(1), cihi(0), cihi(1),
                              filo(0), filo(1), fihi(0), fihi(1),
                              &ratio[0],
                              cdata->getPointer(d),
                              fdata->getPointer(d));

#endif  // test for RAJA
      } else if (fine.getDim() == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->template getView<3>(d);
         auto coarse_array = cdata->template getView<3>(d);
         const int r0 = ratio[0];
         const int r1 = ratio[1];
         const int r2 = ratio[2];

         hier::parallel_for_all(fine_box, [=] SAMRAI_HOST_DEVICE(int i, int j, int k) {
            const int ic2 = (k < 0) ? (k + 1) / r2 - 1 : k / r2;
            const int ic1 = (j < 0) ? (j + 1) / r1 - 1 : j / r1;
            const int ic0 = (i < 0) ? (i + 1) / r0 - 1 : i / r0;

            fine_array(i, j, k) = coarse_array(ic0, ic1, ic2);
         });
#else   // Fortran Dimension 3
         Call3dFortranCell<T>(ifirstc(0), ifirstc(1),
                              ifirstc(2),
                              ilastc(0), ilastc(1), ilastc(2),
                              ifirstf(0), ifirstf(1), ifirstf(2),
                              ilastf(0), ilastf(1), ilastf(2),
                              cilo(0), cilo(1), cilo(2),
                              cihi(0), cihi(1), cihi(2),
                              filo(0), filo(1), filo(2),
                              fihi(0), fihi(1), fihi(2),
                              &ratio[0],
                              cdata->getPointer(d),
                              fdata->getPointer(d));
#endif  // test for RAJA
      } else {
         TBOX_ERROR(
             "CellConstantRefine::refine dimension > 3 not supported"
             << std::endl);
      }  // dimension
   }     // for depth
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
