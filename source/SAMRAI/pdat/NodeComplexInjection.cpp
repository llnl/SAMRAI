/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant averaging operator for node-centered complex data on
 *                a  mesh.
 *
 ************************************************************************/
#include "SAMRAI/pdat/NodeComplexInjection.h"
#include "SAMRAI/tbox/Complex.h"

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
#pragma warning (disable:1419)
#endif

// in concoarsen1d.f:
void SAMRAI_F77_FUNC(conavgnodecplx1d, CONAVGNODECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
// in concoarsen2d.f:
void SAMRAI_F77_FUNC(conavgnodecplx2d, CONAVGNODECPLX2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
// in concoarsen3d.f:
void SAMRAI_F77_FUNC(conavgnodecplx3d, CONAVGNODECPLX3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace pdat {

NodeComplexInjection::NodeComplexInjection():
   hier::CoarsenOperator("CONSTANT_COARSEN")
{
}

NodeComplexInjection::~NodeComplexInjection()
{
}

int
NodeComplexInjection::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
NodeComplexInjection::getStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getZero(dim);
}

void
NodeComplexInjection::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(coarse, fine, coarse_box, ratio);

   std::shared_ptr<NodeData<dcomplex> > fdata(
      SAMRAI_SHARED_PTR_CAST<NodeData<dcomplex>, hier::PatchData>(
         fine.getPatchData(src_component)));
   std::shared_ptr<NodeData<dcomplex> > cdata(
      SAMRAI_SHARED_PTR_CAST<NodeData<dcomplex>, hier::PatchData>(
         coarse.getPatchData(dst_component)));

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();
   const hier::Index& cilo = cdata->getGhostBox().lower();
   const hier::Index& cihi = cdata->getGhostBox().upper();

   const hier::Index& ifirstc = coarse_box.lower();
   const hier::Index& ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); ++d) {
      if (fine.getDim() == tbox::Dimension(1)) {
         SAMRAI_F77_FUNC(conavgnodecplx1d, CONAVGNODECPLX1D) (ifirstc(0), ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(2)) {
         SAMRAI_F77_FUNC(conavgnodecplx2d, CONAVGNODECPLX2D) (ifirstc(0), ifirstc(1),
            ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (fine.getDim() == tbox::Dimension(3)) {
         SAMRAI_F77_FUNC(conavgnodecplx3d, CONAVGNODECPLX3D) (ifirstc(0), ifirstc(1),
            ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR(
            "NodeComplexInjection::coarsen dimension > 3 not supported"
            << std::endl);
      }
   }
}

}
}
