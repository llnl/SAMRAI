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
#include "SAMRAI/pdat/NodeDoubleInjection.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"

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

// in concoarsen1d.f:
void SAMRAI_F77_FUNC(conavgnodedoub1d, CONAVGNODEDOUB1D)(const int&, const int&,
                                                         const int&, const int&,
                                                         const int&, const int&,
                                                         const int*,
                                                         const double*, double*);
// in concoarsen2d.f:
void SAMRAI_F77_FUNC(conavgnodedoub2d, CONAVGNODEDOUB2D)(const int&, const int&,
                                                         const int&, const int&,
                                                         const int&, const int&, const int&, const int&,
                                                         const int&, const int&, const int&, const int&,
                                                         const int*,
                                                         const double*, double*);
// in concoarsen3d.f:
void SAMRAI_F77_FUNC(conavgnodedoub3d, CONAVGNODEDOUB3D)(const int&, const int&,
                                                         const int&,
                                                         const int&, const int&, const int&,
                                                         const int&, const int&, const int&,
                                                         const int&, const int&, const int&,
                                                         const int&, const int&, const int&,
                                                         const int&, const int&, const int&,
                                                         const int*,
                                                         const double*, double*);
}

namespace SAMRAI
{
namespace pdat
{

NodeDoubleInjection::NodeDoubleInjection() : hier::CoarsenOperator("CONSTANT_COARSEN")
{
}

NodeDoubleInjection::~NodeDoubleInjection()
{
}

int NodeDoubleInjection::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
NodeDoubleInjection::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void NodeDoubleInjection::coarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio) const
{
   std::shared_ptr<NodeData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<NodeData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<NodeData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<NodeData<double>, hier::PatchData>(
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
         SAMRAI_F77_FUNC(conavgnodedoub1d, CONAVGNODEDOUB1D)
         (ifirstc(0), ilastc(0),
          filo(0), fihi(0),
          cilo(0), cihi(0),
          &ratio[0],
          fdata->getPointer(d),
          cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(2)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<2>(d);
         auto coarse_array = cdata->getView<2>(d);

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
         SAMRAI_F77_FUNC(conavgnodedoub2d, CONAVGNODEDOUB2D)
         (ifirstc(0), ifirstc(1),
          ilastc(0), ilastc(1),
          filo(0), filo(1), fihi(0), fihi(1),
          cilo(0), cilo(1), cihi(0), cihi(1),
          &ratio[0],
          fdata->getPointer(d),
          cdata->getPointer(d));
#endif
      } else if (fine.getDim() == tbox::Dimension(3)) {
#if defined(HAVE_RAJA)
         auto fine_array = fdata->getView<3>(d);
         auto coarse_array = cdata->getView<3>(d);

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
         SAMRAI_F77_FUNC(conavgnodedoub3d, CONAVGNODEDOUB3D)
         (ifirstc(0), ifirstc(1),
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
             "NodeDoubleConstantRefine::coarsen dimension > 3 not supported"
             << std::endl);
      }
   }
}

}  // namespace pdat
}  // namespace SAMRAI
