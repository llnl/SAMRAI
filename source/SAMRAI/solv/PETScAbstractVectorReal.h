/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to C++ vector implementation for PETSc package.
 *
 ************************************************************************/

#ifndef included_solv_PETScAbstractVectorReal
#define included_solv_PETScAbstractVectorReal

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
 ************************************************************************
 */
#ifdef HAVE_PETSC

#ifdef REQUIRES_CMATH
#include <cmath>
#endif

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#ifndef included_petsc_vec
#define included_petsc_vec
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif
#ifdef OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif
#include "petscvec.h"
#endif
#ifndef included_petsc_vecimpl
#define included_petsc_vecimpl
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif
#ifdef OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif
#include "petsc/private/vecimpl.h"
#endif

namespace SAMRAI {
namespace solv {

/**
 * Class PETScAbstractVectorReal serves as an abstract base class for a
 * <TT>C++</TT> vector class that can be used with the PETSc solver framework.
 * Specifically, this class provides an interface for real-valued PETSc
 * vectors (i.e., where the data is either float or double).  PETSc allows
 * the use of user-defined vectors.  Thus, the intent of this base class is
 * that one may provide his/her own vector implementation in a subclass of
 * this base class that provides the the necessary vector data structures
 * and implements the pure virtual functions.  This class declares private
 * static member functions for linkage with the vector object function calls
 * understood by PETSc.  Each of the static member functions calls an
 * associated function in a subclass via the virtual function mechanism.
 * Note that the virtual members of this class are all protected.  They should
 * not be used outside of a subclass of this class.  The data member of the
 * PETSc vector object is set to an instantiation of the user-supplied vector
 * class, when an object of this class is constructed.
 *
 * PETSc was developed in the Mathematics and Computer Science Division at
 * Argonne National Laboratory (ANL).  For more information about PETSc,
 * see <TT>http://www-fp.mcs.anl.gov/petsc/</TT>.
 *
 * Important notes:
 *
 *
 *
 * - @b (1) The user-supplied vector subclass should only inherit from
 *            this base class. It MUST NOT employ multiple inheritance so that
 *            problems with casting from a base class pointer to a subclass
 *            pointer and the use of virtual functions works properly.
 * - @b (2) The user-supplied subclass that implements the vector data
 *            structures and operations is responsible for all parallelism,
 *            I/O, etc. associated with the vector objects.  PETSc only sees
 *            pointers to what it believes to be sequential vector objects
 *            associated with the local process only.  It has no knowledge
 *            of the structure of the vector data, nor the implementations
 *            of the individual vector routines.
 * - @b (3) Several of the operations defined in the PETSc <TT>_VecOps</TT>
 *            structure are left unimplemented in this class.  They will
 *            print an error message and throw an unrecoverable exeception
 *            if called which causes the program to abort.
 * - @b (4) By default, PETSc typdefs "Scalar" to <TT>double</TT>.  PETSc
 *            must be recompiled to use <TT>float</TT> data.  Also, PETSc
 *            support complex vector data.  A complex vector interface class
 *            similar to this class may be provided in the future if the
 *            need arises.
 *
 *
 *
 */

template<class TYPE>
class PETScAbstractVectorReal
{
protected:
   /**
    * Constructor PETScAbstractVectorReal class that provides a wrapper
    * for a SAMRAI vector so that it can be manipulated within PETSc.  The
    * constructor allocates the PETSc vector and sets its data structures and
    * member functions so that it can operate on the SAMRAI vector.
    */
   PETScAbstractVectorReal(
      bool vector_created_via_duplicate,
      MPI_Comm comm);

   /**
    * Destructor for PETScAbstractVectorReal class destroys the PETSc
    * vector created by the constructor.
    */
   virtual ~PETScAbstractVectorReal();

   /**
    * Return PETSc "Vec" object for this PETScAbstractVectorReal object.
    */
   Vec
   getPETScVector();

   /**
    * Clone the vector structure and allocate storage for the vector
    * data.  Then, return a pointer to the new vector instance.  This
    * function is distinct from the vector constructor since it is called
    * from within PETSc (via the duplicateVec(), duplicateVecs() functions)
    * to allocate new vector objects.
    */
   virtual PETScAbstractVectorReal<TYPE> *
   makeNewVector() = 0;

   /**
    * Destroy vector structure and its storage. This function is distinct
    * from the destructor since it will be called from within PETSc to
    * deallocate vectors (via the destroyVec(), destroyVecs() functions).
    */
   virtual void
   freeVector() = 0;

   /**
    * View all vector data.  Note that the user-supplied vector must
    * choose how to view the vector and its data; e.g., print to file,
    * dump to standard out, etc.
    */
   virtual void
   viewVector() const = 0;

   /**
    * Return @f$ (x,y) = \sum_i ( x_i * std::conj(y_i) ) @f$ , where @f$ x @f$  is this vector.
    * Note that for real vectors, this is the same as TdotWith().
    * If local_only is true, the operation is limited to parts owned by the
    * local process.
    */
   virtual TYPE
   dotWith(
      const PETScAbstractVectorReal<TYPE>* y,
      bool local_only = false) const = 0;

   /**
    * Limited to local data only,
    * return @f$ (x,y) = \sum_i ( x_i * (y_i) ) @f$ , where @f$ x @f$  is this vector.
    * Note that for real vectors, this is the same as dotWith().
    * If local_only is true, the operation is limited to parts owned by the
    * local process.
    */
   virtual TYPE
   TdotWith(
      const PETScAbstractVectorReal<TYPE>* y,
      bool local_only = false) const = 0;

   /**
    * Return @f$ L_1 @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual TYPE
   L1Norm(
      bool local_only = false) const = 0;

   /**
    * Return @f$ L_2 @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual TYPE
   L2Norm(
      bool local_only = false) const = 0;

   /**
    * Return @f$ L_{\infty} @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual TYPE
   maxNorm(
      bool local_only = false) const = 0;

   /**
    * Multiply each entry of this vector by given scalar.
    */
   virtual void
   scaleVector(
      const TYPE alpha) = 0;

   /**
    * Copy source vector data to this vector.
    */
   virtual void
   copyVector(
      const PETScAbstractVectorReal<TYPE>* v_src) = 0;

   /**
    * Set each entry of this vector to given scalar.
    */
   virtual void
   setToScalar(
      const TYPE alpha) = 0;

   /**
    * Swap data between this vector and argument vector.
    */
   virtual void
   swapWith(
      PETScAbstractVectorReal<TYPE>* v_other) = 0;

   /**
    * Set @f$ y = \alpha x + y @f$ , where @f$ y @f$  is this vector.
    */
   virtual void
   setAXPY(
      const TYPE alpha,
      const PETScAbstractVectorReal<TYPE>* x) = 0;

   /**
    * Set @f$ y = \alpha x + @beta y @f$ , where @f$ y @f$  is this vector.
    */
   virtual void
   setAXPBY(
      const TYPE alpha,
      const PETScAbstractVectorReal<TYPE>* x,
      const TYPE beta) = 0;

   /**
    * Set @f$ w = \alpha x + y @f$ , where @f$ w @f$  is this vector.
    */
   virtual void
   setWAXPY(
      const TYPE alpha,
      const PETScAbstractVectorReal<TYPE>* x,
      const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Set @f$ w_i = x_i y_i @f$ , where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   void
   pointwiseMultiply(
      const PETScAbstractVectorReal<TYPE>* x,
      const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Set @f$ w_i = x_i / y_i @f$ , where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   void
   pointwiseDivide(
      const PETScAbstractVectorReal<TYPE>* x,
      const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Compute @f$ max_i = abs(w_i / y_i) @f$ ,
    * where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   double
   maxPointwiseDivide(
      const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Find maximum vector entry and vector index at which maximum occurs.
    */
   virtual void
   vecMax(
      int& i,
      TYPE& max) const = 0;

   /**
    * Find minimum vector entry and vector index at which minimum occurs.
    */
   virtual void
   vecMin(
      int& i,
      TYPE& min) const = 0;

   /**
    * Set vector entries to random values.  Note that PETSc uses the
    * drand48() function and computes random value as width*drand48()+low.
    */
   virtual void
   setRandomValues(
      const TYPE width,
      const TYPE low) = 0;

   /**
    * Set argument to vector data in contiguous array (local to processor).
    */
   virtual void
   getDataArray(
      TYPE** array) = 0;

   /**
    * Return total length of vector.
    */
   virtual int
   getDataSize() const = 0;

   /**
    * Return length of vector (local to processor).
    */
   virtual int
   getLocalDataSize() const = 0;

   /*!
    * Restore pointer to vector data in contiguous array (local to
    * processor).
    */
   virtual void
   restoreDataArray(
      TYPE** array) = 0;

private:
   /*
    * PETSc vector object corresponding to this
    * PETScAbstractVectorReal object.
    */
   Vec d_petsc_vector;

   bool d_vector_created_via_duplicate;

   MPI_Comm d_comm;

   /*
    * Static member functions for linkage with PETSc solver package
    * routines.  Essentially, these functions match those in the
    * PETSc _VecOps structure.  Note that these operations are
    * actually implemented in a subclass of this base class using the
    * virtual function mechanism.
    */

   /*
    * Creates a new vector of the same type as an existing vector.
    */
   static PetscErrorCode
   vecDuplicate(
      Vec v,
      Vec* newv);

   /*
    * Creates an array of vectors of the same type as an existing vector.
    */
   static PetscErrorCode
   vecDuplicateVecs(
      Vec v,
      int n,
      Vec** varr_new);

   /*
    * Destroys an array of vectors.
    */
   static PetscErrorCode
   vecDestroyVecs(
      PetscInt n,
      Vec* v_arr);

   /*
    * Computes the vector dot product.
    */
   static PetscErrorCode
   vecDot(
      Vec x,
      Vec y,
      TYPE* val);

   /*
    * Computes vector multiple dot products.
    */
   static PetscErrorCode
   vecMDot(
      Vec x,
      PetscInt nv,
      const Vec* y,
      TYPE* val);

   /*
    * Computes the vector norm.
    *
    * Note that PETSc defines the following enumerated type (in
    * petscvec.h):
    *
    * typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
    *
    * @pre (type == NORM_1) || (type == NORM_2) || (type == NORM_INFINITY) ||
    *      (type == NORM_1_AND_2)
    */
   static PetscErrorCode
   vecNorm(
      Vec x,
      NormType type,
      TYPE* val);

   /*
    * Computes an indefinite vector dot product.  That is, this
    * routine does NOT use the complex conjugate.
    */
   static PetscErrorCode
   vecTDot(
      Vec x,
      Vec y,
      TYPE* val);

   /*
    * Computes indefinite vector multiple dot products.  That is, it
    * does NOT use the complex conjugate.
    */
   static PetscErrorCode
   vecMTDot(
      Vec x,
      PetscInt nv,
      const Vec* y,
      TYPE* val);

   /*
    * Scales a vector.
    */
   static PetscErrorCode
   vecScale(
      Vec x,
      TYPE alpha);

   /*
    * Copies a vector.
    */
   static PetscErrorCode
   vecCopy(
      Vec x,
      Vec y);

   /*
    *  Sets all components of a vector to a single scalar value.
    */
   static PetscErrorCode
   vecSet(
      Vec x,
      TYPE alpha);

   /*
    * Swaps the vectors x and y.
    */
   static PetscErrorCode
   vecSwap(
      Vec x,
      Vec y);

   /*
    * Computes y = alpha x + y.
    */
   static PetscErrorCode
   vecAXPY(
      Vec y,
      TYPE alpha,
      Vec x);

   /*
    * Computes y = alpha x + beta y.
    */
   static PetscErrorCode
   vecAXPBY(
      Vec y,
      TYPE alpha,
      TYPE beta,
      Vec x);

   /*
    * Computes y = y + sum alpha[j] x[j].
    */
   static PetscErrorCode
   vecMAXPY(
      Vec y,
      PetscInt nv,
      const TYPE* alpha,
      Vec* x);

   /*
    * Computes y = x + alpha y.
    */
   static PetscErrorCode
   vecAYPX(
      Vec y,
      TYPE alpha,
      Vec x);

   /*
    * Computes w = alpha x + y.
    */
   static PetscErrorCode
   vecWAXPY(
      Vec w,
      TYPE alpha,
      Vec x,
      Vec y);

   /*
    * Computes the component-wise multiplication w = x*y.
    */
   static PetscErrorCode
   vecPointwiseMult(
      Vec w,
      Vec x,
      Vec y);

   /*
    * Computes the component-wise division w = x/y.
    */
   static PetscErrorCode
   vecPointwiseDivide(
      Vec w,
      Vec x,
      Vec y);

   /*
    * Returns a pointer to a contiguous array that contains this
    * processor's portion of the vector data.
    */
   static PetscErrorCode
   vecGetArray(
      Vec x,
      TYPE** a);

   /*
    * Returns the global number of elements of the vector.
    */
   static PetscErrorCode
   vecGetSize(
      Vec x,
      PetscInt* size);

   /*
    * Returns the number of elements of the vector stored in local
    * memory.
    */
   static PetscErrorCode
   vecGetLocalSize(
      Vec x,
      PetscInt* size);

   /*
    * Restores a vector after VecGetArray() has been called.
    */
   static PetscErrorCode
   vecRestoreArray(
      Vec x,
      TYPE** a);

   /*
    * Determines the maximum vector component and its location.
    */
   static PetscErrorCode
   vecMax(
      Vec x,
      PetscInt* p,
      TYPE* val);

   /*
    * Determines the minimum vector component and its location.
    */
   static PetscErrorCode
   vecMin(
      Vec x,
      PetscInt* p,
      TYPE* val);

   /*
    * Sets all components of a vector to random numbers.
    */
   static PetscErrorCode
   vecSetRandom(
      Vec x,
      PetscRandom rctx);

   /*
    * Destroys a vector.
    */
   static PetscErrorCode
   vecDestroy(
      Vec v);

   /*
    * Views a vector object.
    */
   static PetscErrorCode
   vecView(
      Vec v,
      PetscViewer viewer);

   /*
    * Computes the vector dot product.
    */
   static PetscErrorCode
   vecDot_local(
      Vec x,
      Vec y,
      TYPE* val);

   /*
    * Computes an indefinite vector dot product.  That is, this
    * routine does NOT use the complex conjugate.
    */
   static PetscErrorCode
   vecTDot_local(
      Vec x,
      Vec y,
      TYPE* val);

   /*
    * Computes the vector norm.
    *
    * Note that PETSc defines the following enumerated type (in
    * petscvec.h):
    *
    * typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
    *
    * @pre (type == NORM_1) || (type == NORM_2) || (type == NORM_INFINITY) ||
    *      (type == NORM_1_AND_2)
    */
   static PetscErrorCode
   vecNorm_local(
      Vec x,
      NormType type,
      TYPE* val);

   /*
    * Computes vector multiple dot products.
    */
   static PetscErrorCode
   vecMDot_local(
      Vec x,
      PetscInt nv,
      const Vec* y,
      TYPE* val);

   /*
    * Computes indefinite vector multiple dot products.  That is, it
    * does NOT use the complex conjugate.
    */
   static PetscErrorCode
   vecMTDot_local(
      Vec x,
      PetscInt nv,
      const Vec* y,
      TYPE* val);

   /*
    * Computes the maximum of the component-wise division max = max_i
    * abs(x_i/y_i).
    */
   static PetscErrorCode
   vecMaxPointwiseDivide(
      Vec x,
      Vec y,
      TYPE* max);

   ///
   /// The remaining functions are not implemented and will result in
   /// an unrecoverable assertion being thrown and program abort if
   /// called.
   ///

   static PetscErrorCode
   vecSetValues(
      Vec x,
      PetscInt ni,
      const PetscInt* ix,
      const TYPE* y,
      InsertMode iora);

   static PetscErrorCode
   vecAssemblyBegin(
      Vec vec);

   static PetscErrorCode
   vecAssemblyEnd(
      Vec vec);

   static PetscErrorCode
   vecSetOption(
      Vec x,
      VecOption op,
      PetscBool result);

   static PetscErrorCode
   vecSetValuesBlocked(
      Vec x,
      PetscInt ni,
      const PetscInt* ix,
      const TYPE* y,
      InsertMode iora);

   static PetscErrorCode
   vecPlaceArray(
      Vec vec,
      const TYPE* array);

   static PetscErrorCode
   vecReplaceArray(
      Vec vec,
      const TYPE* array);

   static PetscErrorCode
   vecReciprocal(
      Vec vec);

   static PetscErrorCode
   vecViewNative(
      Vec v,
      PetscViewer viewer);

   static PetscErrorCode
   vecConjugate(
      Vec x);

   static PetscErrorCode
   vecSetLocalToGlobalMapping(
      Vec x,
      ISLocalToGlobalMapping mapping);

   static PetscErrorCode
   vecSetValuesLocal(
      Vec x,
      PetscInt ni,
      const PetscInt* ix,
      const TYPE* y,
      InsertMode iora);

   static PetscErrorCode
   vecResetArray(
      Vec vec);

   static PetscErrorCode
   vecSetFromOptions(
      Vec vec,
      PetscOptionItems* options);

   static PetscErrorCode
   vecLoad(
      Vec newvec,
      PetscViewer viewer);

   static PetscErrorCode
   vecPointwiseMax(
      Vec w,
      Vec x,
      Vec y);

   static PetscErrorCode
   vecPointwiseMaxAbs(
      Vec w,
      Vec x,
      Vec y);

   static PetscErrorCode
   vecPointwiseMin(
      Vec w,
      Vec x,
      Vec y);

   static PetscErrorCode
   vecGetValues(
      Vec x,
      PetscInt ni,
      const PetscInt* ix,
      PetscScalar* y);
};

}
}

#include "SAMRAI/solv/PETScAbstractVectorReal.cpp"

#endif
#endif
