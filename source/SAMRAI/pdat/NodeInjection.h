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

#ifndef included_pdat_NodeInjection
#define included_pdat_NodeInjection

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/InvokeOne.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in concoarsen1d.f:
void SAMRAI_F77_FUNC(conavgnodedoub1d, CONAVGNODEDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conavgnodeflot1d, CONAVGNODEFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conavgnodecplx1d, CONAVGNODECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conavgnodeintg1d, CONAVGNODEINTG1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);
// in concoarsen2d.f:
void SAMRAI_F77_FUNC(conavgnodedoub2d, CONAVGNODEDOUB2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conavgnodeflot2d, CONAVGNODEFLOT2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conavgnodecplx2d, CONAVGNODECPLX2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conavgnodeintg2d, CONAVGNODEINTG2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in concoarsen3d.f:
void SAMRAI_F77_FUNC(conavgnodedoub3d, CONAVGNODEDOUB3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conavgnodeflot3d, CONAVGNODEFLOT3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conavgnodecplx3d, CONAVGNODECPLX3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conavgnodeintg3d, CONAVGNODEINTG3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const int *, int *);
} // end extern "C"

namespace SAMRAI {
namespace pdat {

template<typename T>
void Call1dFortranNode(const int& ifirstc, const int& ilastc,
   const int& filo, const int& fihi,
   const int& cilo, const int& cihi,
   const int *ratio,
   const T *arrayf, T *arrayc)
{
   invokeOneOfFour(SAMRAI_F77_FUNC(conavgnodedoub1d, CONAVGNODEDOUB1D),
              SAMRAI_F77_FUNC(conavgnodeflot1d, CONAVGNODEFLOT1D), 
              SAMRAI_F77_FUNC(conavgnodecplx1d, CONAVGNODECPLX1D),
              SAMRAI_F77_FUNC(conavgnodeintg1d, CONAVGNODEINTG1D),
         ifirstc, ilastc,
         filo, fihi,cilo, cihi,
         ratio,
         arrayf,arrayc);

}

template<typename T>
void Call2dFortranNode(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int *ratio,
   const T *arrayf, T *arrayc)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conavgnodedoub2d, CONAVGNODEDOUB2D),
              SAMRAI_F77_FUNC(conavgnodeflot2d, CONAVGNODEFLOT2D), 
              SAMRAI_F77_FUNC(conavgnodecplx2d, CONAVGNODECPLX2D),
              SAMRAI_F77_FUNC(conavgnodeintg2d, CONAVGNODEINTG2D),
         ifirstc0,ifirstc1,
         ilastc0,ilastc1,
         filo0,filo1,fihi0,fihi1,
         cilo0,cilo1,cihi0,cihi1,
         ratio,
         arrayf,arrayc);
}

template<typename T>
void Call3dFortranNode(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
   const int& ilastc0, const int& ilastc1, const int& ilastc2,
   const int& filo0, const int& filo1, const int& filo2,
   const int& fihi0, const int& fihi1, const int& fihi2,
   const int& cilo0, const int& cilo1, const int& cilo2,
   const int& cihi0, const int& cihi1, const int& cihi2,
   const int *ratio,
   const T *arrayf, T *arrayc)
{

    invokeOneOfFour(SAMRAI_F77_FUNC(conavgnodedoub3d, CONAVGNODEDOUB3D),
              SAMRAI_F77_FUNC(conavgnodeflot3d, CONAVGNODEFLOT3D), 
              SAMRAI_F77_FUNC(conavgnodecplx3d, CONAVGNODECPLX3D),
              SAMRAI_F77_FUNC(conavgnodeintg3d, CONAVGNODEINTG3D),
         ifirstc0,ifirstc1,ifirstc2,
         ilastc0,ilastc1,ilastc2,
         filo0,filo1,filo2,fihi0,fihi1,fihi2,
         cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
         ratio,
         arrayf,arrayc);
}


/**
 * Class NodeDoubleInjection implements constant
 * averaging (i.e., injection) for node-centered double patch data defined
 * over a  mesh.  It is derived from the hier::CoarsenOperator base
 * class.  The numerical operations for theaveraging use FORTRAN numerical
 * routines.
 *
 * @see hier::CoarsenOperator
 */

template<typename T>  // one of double,float,dcomplex,int  
class NodeInjection:
   public hier::CoarsenOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   NodeInjection():
   hier::CoarsenOperator("CONSTANT_COARSEN")
   {
   }

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~NodeInjection()
   {
   }


   /**
    * The priority of node-centered constant averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const
   {
      return 0;
   }

   /**
    * The stencil width of the constant averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector::getZero(dim);
   }
   
   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the node-centered double constant
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

} // namespace pdat
} // namespace SAMRAI 
#include "SAMRAI/pdat/NodeInjection.cpp"
#endif
