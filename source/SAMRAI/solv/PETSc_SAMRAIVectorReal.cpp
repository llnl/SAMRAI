/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between PETSc vector interface and SAMRAI vectors.
 *
 ************************************************************************/

#ifndef included_solv_PETSc_SAMRAIVectorReal_C
#define included_solv_PETSc_SAMRAIVectorReal_C

#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"

#ifdef HAVE_PETSC

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/PIO.h"

#include <cstdlib>

namespace SAMRAI {
namespace solv {

#define C_PSVEC_CAST(x) \
   (dynamic_cast<const PETSc_SAMRAIVectorReal< \
                    TYPE> *>(x))

/*
 *************************************************************************
 *
 * Static public member functions.
 *
 *************************************************************************
 */

template<class TYPE>
Vec
PETSc_SAMRAIVectorReal<TYPE>::createPETScVector(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& samrai_vec,
   MPI_Comm comm)
{
   TBOX_ASSERT(samrai_vec);

   static const bool vector_created_via_duplicate = false;

   PETSc_SAMRAIVectorReal<TYPE>* psv = new PETSc_SAMRAIVectorReal<TYPE>(
         samrai_vec, vector_created_via_duplicate, comm);

   return psv->getPETScVector();
}

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::destroyPETScVector(
   Vec petsc_vec)
{
   if (petsc_vec != 0) {
      PETSc_SAMRAIVectorReal<TYPE>* psv =
         static_cast<PETSc_SAMRAIVectorReal<TYPE> *>(petsc_vec->data);

      TBOX_ASSERT(psv != 0);

      delete psv;
   }
}

template<class TYPE>
std::shared_ptr<SAMRAIVectorReal<TYPE> >
PETSc_SAMRAIVectorReal<TYPE>::getSAMRAIVector(
   Vec petsc_vec)
{
   TBOX_ASSERT(petsc_vec != 0);

   PETSc_SAMRAIVectorReal<TYPE>* psv =
      static_cast<PETSc_SAMRAIVectorReal<TYPE> *>(petsc_vec->data);

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
   TBOX_ASSERT(psv != 0);
#endif

   return psv->d_samrai_vector;
}

/*
 *************************************************************************
 *
 * Protected constructor and destructor for PETSc_SAMRAIVectorReal.
 *
 *************************************************************************
 */

template<class TYPE>
PETSc_SAMRAIVectorReal<TYPE>::PETSc_SAMRAIVectorReal(
   const std::shared_ptr<SAMRAIVectorReal<TYPE> >& samrai_vector,
   bool vector_created_via_duplicate,
   MPI_Comm comm):
   PETScAbstractVectorReal<TYPE>(vector_created_via_duplicate, comm),
   d_samrai_vector(samrai_vector),
   d_vector_created_via_duplicate(vector_created_via_duplicate)
{
   // intentionally blank
}

template<class TYPE>
PETSc_SAMRAIVectorReal<TYPE>::~PETSc_SAMRAIVectorReal()
{
   // intentionally blank
}

/*
 *************************************************************************
 *
 * Other member functions
 *
 *************************************************************************
 */

template<class TYPE>
PETScAbstractVectorReal<TYPE> *
PETSc_SAMRAIVectorReal<TYPE>::makeNewVector()
{

   Vec petsc_vec = PETSc_SAMRAIVectorReal<TYPE>::getPETScVector();
   MPI_Comm comm;
   int ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(petsc_vec),
         &comm);
   PETSC_SAMRAI_ERROR(ierr);

   std::shared_ptr<SAMRAIVectorReal<TYPE> > sam_vec(
      d_samrai_vector->cloneVector(d_samrai_vector->getName()));
   sam_vec->allocateVectorData();
   const bool vector_created_via_duplicate = true;
   PETSc_SAMRAIVectorReal<TYPE>* out_vec =
      new PETSc_SAMRAIVectorReal<TYPE>(sam_vec,
                                       vector_created_via_duplicate,
                                       comm);
   return out_vec;
}

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::freeVector()
{

   if (d_vector_created_via_duplicate) {
      d_samrai_vector->freeVectorComponents();
      d_samrai_vector.reset();
      Vec petsc_vec = this->getPETScVector();

#ifdef DEBUG_CHECK_TBOX_ASSERTIONS
      TBOX_ASSERT(petsc_vec != 0);
#endif
      delete ((PETSc_SAMRAIVectorReal<TYPE> *)(petsc_vec->data));
   }
}

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::viewVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting PETSc_SAMRAIVectorReal..."
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::dotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
} // dotWith

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::TdotWith(
   const PETScAbstractVectorReal<TYPE>* y,
   bool local_only) const
{
   return d_samrai_vector->dot(C_PSVEC_CAST(y)->d_samrai_vector, local_only);
} // TdotWith

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::L1Norm(
   bool local_only) const
{
   return d_samrai_vector->L1Norm(local_only);
} // L1Norm

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::L2Norm(
   bool local_only) const
{
   return d_samrai_vector->L2Norm(local_only);
} // L2Norm

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::maxNorm(
   bool local_only) const
{
   return d_samrai_vector->maxNorm(local_only);
} // maxNorm

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::scaleVector(
   const TYPE alpha)
{
   d_samrai_vector->scale(alpha, d_samrai_vector);
} // scaleVector

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::copyVector(
   const PETScAbstractVectorReal<TYPE>* v_src)
{
   d_samrai_vector->copyVector(C_PSVEC_CAST(v_src)->d_samrai_vector);
} // copyVector

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::setToScalar(
   const TYPE alpha)
{
   d_samrai_vector->setToScalar(alpha);
} // setToScalar

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::swapWith(
   PETScAbstractVectorReal<TYPE>* v_other)
{
   d_samrai_vector->swapVectors(C_PSVEC_CAST(v_other)->d_samrai_vector);
} // swapWith

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::setAXPY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x)
{
   d_samrai_vector->axpy(alpha, C_PSVEC_CAST(
         x)->d_samrai_vector, d_samrai_vector);
} // setAXPY

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::setAXPBY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x,
   const TYPE beta)
{
   d_samrai_vector->linearSum(alpha, C_PSVEC_CAST(
         x)->d_samrai_vector, beta, d_samrai_vector);
} // setAXPBY

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::setWAXPY(
   const TYPE alpha,
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->axpy(alpha, C_PSVEC_CAST(x)->d_samrai_vector,
      C_PSVEC_CAST(y)->d_samrai_vector);
} // setWAXPY

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::pointwiseMultiply(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->multiply(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(
         y)->d_samrai_vector);
} // pointwiseMultiply

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::pointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* x,
   const PETScAbstractVectorReal<TYPE>* y)
{
   d_samrai_vector->divide(C_PSVEC_CAST(x)->d_samrai_vector, C_PSVEC_CAST(
         y)->d_samrai_vector);
} // pointwiseDivide

template<class TYPE>
double
PETSc_SAMRAIVectorReal<TYPE>::maxPointwiseDivide(
   const PETScAbstractVectorReal<TYPE>* y)
{
   return d_samrai_vector->maxPointwiseDivide(C_PSVEC_CAST(y)->d_samrai_vector);
} // maxPointwiseDivide

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::vecMax(
   int& i,
   TYPE& max) const
{
   static const bool interior_only = true;
   max = d_samrai_vector->max(interior_only);
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
} // vecMax

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::vecMin(
   int& i,
   TYPE& min) const
{
   static const bool interior_only = true;
   min = d_samrai_vector->min(interior_only);
   // Note: This is a bogus index value!
   //       Hopefully, PETSc doesn't use it for anything.
   i = 0;
} // vecMin

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::setRandomValues(
   const TYPE width,
   const TYPE low)
{
   d_samrai_vector->setRandomValues(width, low);
} // setRandomValues

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::getDataArray(
   TYPE** array)
{
   *array = 0;
} // getDataArray

template<class TYPE>
void
PETSc_SAMRAIVectorReal<TYPE>::restoreDataArray(
   TYPE** array)
{
   *array = 0;
} // restoreDataArray

template<class TYPE>
int
PETSc_SAMRAIVectorReal<TYPE>::getDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return 0;
} // getDataSize

template<class TYPE>
int
PETSc_SAMRAIVectorReal<TYPE>::getLocalDataSize() const
{
   // Note: This is a bogus value!
   //       But, PETSc requires some value to be returned.
   //       Hopefully, this will not cause problems.
   return 0;
} // getLocalDataSize

}
}
#endif
#endif
