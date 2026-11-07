/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to C++ vector implementation for Sundials package.
 *
 ************************************************************************/

#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cstdlib>

#ifdef HAVE_SUNDIALS

namespace SAMRAI {
namespace solv {

SundialsAbstractVector::SundialsAbstractVector()
{
   /* Create N vector */
   d_n_vector = (N_Vector)malloc(sizeof *d_n_vector);
   TBOX_ASSERT(d_n_vector != 0);

   /* Attach content and ops */
   d_n_vector->content = this;
   d_n_vector->ops = SundialsAbstractVector::createVectorOps();
}

SundialsAbstractVector::~SundialsAbstractVector()
{
   if (d_n_vector) {
      if (d_n_vector->ops) {
         free(d_n_vector->ops);
         d_n_vector->ops = 0;
      }
      free(d_n_vector);
      d_n_vector = 0;
   }
}

N_Vector_Ops
SundialsAbstractVector::createVectorOps()
{
   /* Create vector operation structure */

   N_Vector_Ops ops;

   ops = (N_Vector_Ops) std::calloc(1, sizeof(struct _generic_N_Vector_Ops));

// SGS TODO what about missing fns?
   ops->nvclone = SundialsAbstractVector::N_VClone_SAMRAI;
//   ops->nvcloneempty      = N_VCloneEmpty_SAMRAI;
   ops->nvdestroy = N_VDestroy_SAMRAI;
//   ops->nvspace           = N_VSpace_SAMRAI;
//   ops->nvgetarraypointer = N_VGetArrayPointer_SAMRAI;
//   ops->nvsetarraypointer = N_VSetArrayPointer_SAMRAI;
   ops->nvlinearsum = N_VLinearSum_SAMRAI;
   ops->nvconst = N_VConst_SAMRAI;
   ops->nvprod = N_VProd_SAMRAI;
   ops->nvdiv = N_VDiv_SAMRAI;
   ops->nvscale = N_VScale_SAMRAI;
   ops->nvabs = N_VAbs_SAMRAI;
   ops->nvinv = N_VInv_SAMRAI;
   ops->nvaddconst = N_VAddConst_SAMRAI;
   ops->nvdotprod = N_VDotProd_SAMRAI;
   ops->nvmaxnorm = N_VMaxNorm_SAMRAI;
//   ops->nvwrmsnormmask    = N_VWrmsNormMask_SAMRAI;
   ops->nvwrmsnorm = N_VWrmsNorm_SAMRAI;
   ops->nvmin = N_VMin_SAMRAI;
   ops->nvwl2norm = N_VWL2Norm_SAMRAI;
   ops->nvl1norm = N_VL1Norm_SAMRAI;
   ops->nvcompare = N_VCompare_SAMRAI;
   ops->nvinvtest = N_VInvTest_SAMRAI;
   ops->nvgetlength = N_VGetLength_SAMRAI;
//   ops->nvconstrmask      = N_VConstrMask_SAMRAI;
//   ops->nvminquotient     = N_VMinQuotient_SAMRAI;

   return ops;
}

}
}

#endif
