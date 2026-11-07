/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered
 *                double data on a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianCellConservativeLinearRefine
#define included_geom_CartesianCellConservativeLinearRefine

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/InvokeOne.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>
#include <memory>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning(disable : 1419)
#endif

// in cartrefine1d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D)(const int &,
                                                                   const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *,
                                                                   double *, double *);
void SAMRAI_F77_FUNC(cartclinrefcellflot1d, CARTCLINREFCELLFLOT1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *);
void SAMRAI_F77_FUNC(cartclinrefcellcplx1d, CARTCLINREFCELLCPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *,
   dcomplex *, dcomplex *);

// in cartrefine2d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D)(const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int &, const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *,
                                                                   double *, double *, double *, double *);
void SAMRAI_F77_FUNC(cartclinrefcellflot2d, CARTCLINREFCELLFLOT2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const float *, float *,
   float *, float *, float *, float *);
void SAMRAI_F77_FUNC(cartclinrefcellcplx2d, CARTCLINREFCELLCPLX2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *,
   dcomplex *, dcomplex *, dcomplex *, dcomplex *);

// in cartrefine3d.f:
void SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D)(const int &,
                                                                   const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int &, const int &, const int &,
                                                                   const int *, const double *, const double *,
                                                                   const double *, double *,
                                                                   double *, double *, double *,
                                                                   double *, double *, double *);
void SAMRAI_F77_FUNC(cartclinrefcellflot3d, CARTCLINREFCELLFLOT3D) (const int&,
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

} // extern "C"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace geom {

template<typename T>
void Call1dFortranCellLinearRefine(const int& ifirstc, const int& ilastc,
   const int& ifirstf, const int& ilastf,
   const int& cilo, const int& cihi,
   const int& filo, const int& fihi,
   const int *ratio,
   const double *cdx, const double *fdx,
   const T *arrayc, T *arrayf,
   T *diff_f, T *slope_f)
{
   invokeOneOfThree(SAMRAI_F77_FUNC(cartclinrefcelldoub1d, CARTCLINREFCELLDOUB1D),
              SAMRAI_F77_FUNC(cartclinrefcellflot1d, CARTCLINREFCELLFLOT1D), 
              SAMRAI_F77_FUNC(cartclinrefcellcplx1d, CARTCLINREFCELLCPLX1D),
         ifirstc, ilastc, ifirstf, ilastf,
         cilo, cihi,filo, fihi,
         ratio,
         cdx,fdx,
         arrayc,arrayf,
         diff_f,slope_f);
}

template<typename T>
void Call2dFortranCellLinearRefine(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0,const int& ilastc1,   
   const int& ifirstf0, const int& ifirstf1,
   const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1,
   const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1,
   const int& fihi0, const int& fihi1,
   const int *ratio,
   const double *cdx, const double *fdx,
   const T *arrayc, T *arrayf,
   T *diff_f0, T *slope_f0,
   T *diff_f1, T *slope_f1)
{
   invokeOneOfThree(SAMRAI_F77_FUNC(cartclinrefcelldoub2d, CARTCLINREFCELLDOUB2D),
              SAMRAI_F77_FUNC(cartclinrefcellflot2d, CARTCLINREFCELLFLOT2D), 
              SAMRAI_F77_FUNC(cartclinrefcellcplx2d, CARTCLINREFCELLCPLX2D),
         ifirstc0, ifirstc1, ilastc0, ilastc1,
         ifirstf0, ifirstf1, ilastf0, ilastf1,
         cilo0, cilo1, cihi0, cihi1,filo0, filo1, fihi0, fihi1,
         ratio,
         cdx,fdx,
         arrayc,arrayf,
         diff_f0,slope_f0,
         diff_f1,slope_f1);

}


template<typename T>
void Call3dFortranCellLinearRefine(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
   const int& ilastc0,const int& ilastc1, const int& ilastc2,  
   const int& ifirstf0, const int& ifirstf1, const int& ifirstf2,
   const int& ilastf0, const int& ilastf1, const int& ilastf2,
   const int& cilo0, const int& cilo1, const int& cilo2,
   const int& cihi0, const int& cihi1, const int& cihi2,
   const int& filo0, const int& filo1, const int& filo2,
   const int& fihi0, const int& fihi1, const int& fihi2,
   const int *ratio,
   const double *cdx, const double *fdx,
   const T *arrayc, T *arrayf,
   T *diff_f0, T *slope_f0, 
   T *diff_f1, T *slope_f1, 
   T *diff_f2, T *slope_f2)
{
   invokeOneOfThree(SAMRAI_F77_FUNC(cartclinrefcelldoub3d, CARTCLINREFCELLDOUB3D),
              SAMRAI_F77_FUNC(cartclinrefcellflot3d, CARTCLINREFCELLFLOT3D), 
              SAMRAI_F77_FUNC(cartclinrefcellcplx3d, CARTCLINREFCELLCPLX3D),
         ifirstc0, ifirstc1, ifirstc2,
         ilastc0, ilastc1, ilastc2,
         ifirstf0, ifirstf1, ifirstf2,
         ilastf0, ilastf1, ilastf2,
         cilo0, cilo1, cilo2, cihi0, cihi1, cihi2,
         filo0, filo1, filo2, fihi0, fihi1, fihi2,
         ratio,
         cdx,fdx,
         arrayc,arrayf,
         diff_f0,slope_f0,
         diff_f1,slope_f1,
         diff_f2,slope_f2);

}

/**
 * Class CartesianCellDoubleConservativeLinearRefine implements
 * conservative linear interpolation for cell-centered double patch data
 * defined over a Cartesian mesh.  It is derived from the base class
 * hier::RefineOperator.  The numerical operations for the interpolation
 * use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */
template<typename T> // one of double,float,dcomplex
class CartesianCellConservativeLinearRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CartesianCellConservativeLinearRefine() : hier::RefineOperator("CONSERVATIVE_LINEAR_REFINE")
   {
   }

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CartesianCellConservativeLinearRefine()
   {
   }

   /**
    * The priority of cell-centered double conservative linear is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the conservative linear interpolation operator is
    * the vector of ones.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double conservative
    * linear interpolation operator.  Interpolation is performed on the
    * intersection of the destination patch and the boxes contained in
    * fine_overlap.  It is assumed that the coarse patch contains sufficient
    * data for the stencil width of the refinement operator.
    *
    * @pre dynamic_cast<const pdat::CellOverlap *>(&fine_overlap) != 0
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::BoxOverlap& fine_overlap,
      const hier::IntVector& ratio) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double conservative
    * linear interpolation operator.  Interpolation is performed on the
    * intersection of the destination patch and the fine box.  It is assumed
    * that the coarse patch contains sufficient data for the stencil width of
    * the refinement operator.  This differs from the above refine() method
    * only in that it operates on a single fine box instead of a BoxOverlap.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == fine_box.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    * @pre coarse.getPatchData(src_component) is actually a std::shared_ptr<pdat::CellData<double> >
    * @pre fine.getPatchData(dst_component) is actually a std::shared_ptr<pdat::CellData<double> >
    * @pre coarse.getPatchData(src_component)->getDepth() == fine.getPatchData(dst_component)->getDepth()
    * @pre (fine.getDim().getValue() == 1) ||
    *      (fine.getDim().getValue() == 2) || (fine.getDim().getValue() == 3)
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;

};

}
}
#endif
