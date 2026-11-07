/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant averaging operator for node-centered double data on
 *                a  mesh.
 *
 ************************************************************************/
#ifndef included_pdat_NodeInjection_C
#define included_pdat_NodeInjection_C
#include "SAMRAI/pdat/NodeInjection.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"


namespace SAMRAI
{
namespace pdat
{

template <typename T>
void NodeInjection<T>::coarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio) const
{
   std::shared_ptr<NodeData<T> > fdata(
       SAMRAI_SHARED_PTR_CAST<NodeData<T>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<NodeData<T> > cdata(
       SAMRAI_SHARED_PTR_CAST<NodeData<T>, hier::PatchData>(
           coarse.getPatchData(dst_component)));

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_ASSERT_OBJDIM_EQUALITY4(coarse, fine, coarse_box, ratio);

   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();
   const hier::Index& cilo = cdata->getGhostBox().lower();
   const hier::Index& cihi = cdata->getGhostBox().upper();

   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); ++d) {
      if (fine.getDim() == tbox::Dimension(1)) {
         Call1dFortranNode(ifirstc(0), ilastc(0),
                           filo(0), fihi(0),
                           cilo(0), cihi(0),
                           &ratio[0],
                           fdata->getPointer(d),
                           cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->template getView<2>(d);
         auto coarse_array = cdata->template getView<2>(d);

         SAMRAI::hier::Box coarse_box_plus = coarse_box;
         coarse_box_plus.setUpper(0, cihi(0) + 1);
         coarse_box_plus.setUpper(1, cihi(1) + 1);

         const int r0 = ratio[0];
         const int r1 = ratio[1];

         hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int j /*fastest*/, int k) {
            int if0 = j * r0;
            int if1 = k * r1;
            coarse_array(j, k) = fine_array(if0, if1);
         });
#else
         Call2dFortranNode(ifirstc(0), ifirstc(1),
                           ilastc(0), ilastc(1),
                           filo(0), filo(1), fihi(0), fihi(1),
                           cilo(0), cilo(1), cihi(0), cihi(1),
                           &ratio[0],
                           fdata->getPointer(d),
                           cdata->getPointer(d));
#endif
      } else if (fine.getDim() == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->template getView<3>(d);
         auto coarse_array = cdata->template getView<3>(d);

         SAMRAI::hier::Box coarse_box_plus = coarse_box;
         coarse_box_plus.setUpper(0, cihi(0) + 1);
         coarse_box_plus.setUpper(1, cihi(1) + 1);
         coarse_box_plus.setUpper(2, cihi(2) + 1);

         const int r0 = ratio[0];
         const int r1 = ratio[1];
         const int r2 = ratio[2];

         hier::parallel_for_all(coarse_box_plus, [=] SAMRAI_HOST_DEVICE(int i /*fastest*/, int j, int k) {
            int if0 = i * r0;
            int if1 = j * r1;
            int if2 = k * r2;
            coarse_array(i, j, k) = fine_array(if0, if1, if2);
         });
#else
         Call3dFortranNode(ifirstc(0), ifirstc(1),
                           ifirstc(2),
                           ilastc(0), ilastc(1), ilastc(2),
                           filo(0), filo(1), filo(2),
                           fihi(0), fihi(1), fihi(2),
                           cilo(0), cilo(1), cilo(2),
                           cihi(0), cihi(1), cihi(2),
                           &ratio[0],
                           fdata->getPointer(d),
                           cdata->getPointer(d));
#endif
      } else {
         TBOX_ERROR(
             "NodeConstantRefine::coarsen dimension > 3 not supported"
             << std::endl);
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
