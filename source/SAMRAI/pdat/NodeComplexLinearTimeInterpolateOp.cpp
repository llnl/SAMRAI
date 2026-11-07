/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for node-centered complex data.
 *
 ************************************************************************/
#include "SAMRAI/pdat/NodeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/tbox/Complex.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"


/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in lintimint1d.f:
void SAMRAI_F77_FUNC(lintimeintnodecmplx1d, LINTIMEINTNODECMPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
// in lintimint2d.f:
void SAMRAI_F77_FUNC(lintimeintnodecmplx2d, LINTIMEINTNODECMPLX2D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
// in lintimint3d.f:
void SAMRAI_F77_FUNC(lintimeintnodecmplx3d, LINTIMEINTNODECMPLX3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const dcomplex *, const dcomplex *,
   dcomplex *);
}

namespace SAMRAI {
namespace pdat {

NodeComplexLinearTimeInterpolateOp::NodeComplexLinearTimeInterpolateOp():
   hier::TimeInterpolateOperator()
{
}

NodeComplexLinearTimeInterpolateOp::~NodeComplexLinearTimeInterpolateOp()
{
}

void
NodeComplexLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::BoxOverlap& overlap,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   const tbox::Dimension& dim(where.getDim());

   const NodeData<dcomplex>* old_dat =
      CPP_CAST<const NodeData<dcomplex> *>(&src_data_old);
   const NodeData<dcomplex>* new_dat =
      CPP_CAST<const NodeData<dcomplex> *>(&src_data_new);
   NodeData<dcomplex>* dst_dat =
      CPP_CAST<NodeData<dcomplex> *>(&dst_data);

   TBOX_ASSERT(old_dat != 0);
   TBOX_ASSERT(new_dat != 0);
   TBOX_ASSERT(dst_dat != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(dst_data, where, src_data_old, src_data_new);

   hier::Box node_where = NodeGeometry::toNodeBox(where);

   const hier::Index& old_ilo = old_dat->getGhostBox().lower();
   const hier::Index& old_ihi = old_dat->getGhostBox().upper();
   const hier::Index& new_ilo = new_dat->getGhostBox().lower();
   const hier::Index& new_ihi = new_dat->getGhostBox().upper();

   const hier::Index& dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index& dst_ihi = dst_dat->getGhostBox().upper();

   const double old_time = old_dat->getTime();
   const double new_time = new_dat->getTime();
   const double dst_time = dst_dat->getTime();

   TBOX_ASSERT((old_time < dst_time ||
                tbox::MathUtilities<double>::equalEps(old_time, dst_time)) &&
      (dst_time < new_time ||
       tbox::MathUtilities<double>::equalEps(dst_time, new_time)));

   double tfrac = dst_time - old_time;
   double denom = new_time - old_time;
   if (denom > tbox::MathUtilities<double>::getMin()) {
      tfrac /= denom;
   } else {
      tfrac = 0.0;
   }

   const NodeOverlap* node_overlap = CPP_CAST<const NodeOverlap*>(&overlap);
   hier::BoxContainer ovlp_boxes;
   node_overlap->getSourceBoxContainer(ovlp_boxes);

   for (auto itr = ovlp_boxes.begin(); itr != ovlp_boxes.end(); ++itr) {
      hier::Box dest_box((*itr) * node_where);
      TBOX_ASSERT((dest_box * old_dat->getArrayData().getBox()).isSpatiallyEqual(dest_box));
      TBOX_ASSERT((dest_box * new_dat->getArrayData().getBox()).isSpatiallyEqual(dest_box));
      TBOX_ASSERT((dest_box * dst_dat->getArrayData().getBox()).isSpatiallyEqual(dest_box));

      const hier::Index& ifirst = dest_box.lower();
      const hier::Index& ilast = dest_box.upper();

      for (int d = 0; d < dst_dat->getDepth(); ++d) {
         if (dim == tbox::Dimension(1)) {
            SAMRAI_F77_FUNC(lintimeintnodecmplx1d, LINTIMEINTNODECMPLX1D) (
               ifirst(0), ilast(0),
               old_ilo(0), old_ihi(0),
               new_ilo(0), new_ihi(0),
               dst_ilo(0), dst_ihi(0),
               tfrac,
               old_dat->getPointer(d),
               new_dat->getPointer(d),
               dst_dat->getPointer(d));
         } else if (dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(lintimeintnodecmplx2d, LINTIMEINTNODECMPLX2D) (
               ifirst(0), ifirst(1), ilast(0), ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(d),
               new_dat->getPointer(d),
               dst_dat->getPointer(d));
         } else if (dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(lintimeintnodecmplx3d, LINTIMEINTNODECMPLX3D) (
               ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(d),
               new_dat->getPointer(d),
               dst_dat->getPointer(d));
         } else {
            TBOX_ERROR(
               "NodeComplexLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
               << std::endl);
         }
      }
   }
}

}
}
