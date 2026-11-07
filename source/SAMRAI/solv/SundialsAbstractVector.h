/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to C++ vector kernel operations for Sundials package.
 *
 ************************************************************************/

#ifndef included_solv_SundialsAbstractVector
#define included_solv_SundialsAbstractVector

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_SUNDIALS

#define SABSVEC_CAST(v) \
   (static_cast<SundialsAbstractVector *>(v \
                                          -> \
                                          content))

#ifndef included_sundials_nvector
#include "sundials/sundials_nvector.h"
#define included_sundials_nvector
#endif

namespace SAMRAI {
namespace solv {

/**
 * Class SundialsAbstractVector is an abstract base class that declares
 * operations provided by any <TT>C++</TT> class that may be used as the
 * vector kernel by the Sundials nonlinear solver package.  Sundials allows
 * arbitrarily defined vectors to be used within it as long as the proper
 * collection of operations are provided.  The intent of this base class
 * to provide the interface for one's own vector kernel.  One implements
 * a subclass that provides the functions declared herein as pure virtual
 * and which provides necessary vector data structures.  Note that the
 * virtual members of this class are all protected.  They should not be used
 * outside of a subclass of this class.
 *
 * Sundials was developed in the Center for Applied Scientific Computing (CASC)
 * at Lawrence Livermore National Laboratory (LLNL).  For more information
 * about Sundials, see A.G. Taylor and A.C. Hindmarsh, "User documentation for
 * Sundials, a nonlinear solver for sequential and parallel computers",
 * UCRL-ID-131185, Lawrence Livermore National Laboratory, 1998.
 *
 * Important notes:
 *
 *
 *
 * - \b (1) The user-supplied vector subclass should only inherit from
 *            this base class. It MUST NOT employ multiple inheritance so that
 *            problems with casting from a base class pointer to a subclass
 *            pointer and the use of virtual functions works properly.
 * - \b (2) The user-supplied subclass that implements the vector data
 *            structures and operations is responsible for all parallelism,
 *            I/O, etc. associated with the vector objects.  Sundials is
 *            implemented using a SPMD programming model.  It has no knowlege
 *            of the structure of the vector data, nor the implementations
 *            of the individual vector routines.
 * - \b (3) We assume the vector data is <TT>double</TT>, which is the
 *            default for Sundials.
 *
 *
 *
 *
 * @see SundialsSolver
 */

class SundialsAbstractVector
{
public:
   /**
    * Uninteresting constructor and destructor for SundialsAbstractVector.
    */
   SundialsAbstractVector();
   virtual ~SundialsAbstractVector();

   /**
    * Clone the vector structure and allocate storage for the vector
    * data.  Then, return a pointer to the new vector instance.  Note that
    * the new vector object must be distinct from the original.  This
    * function is distinct from the vector constructor since it will
    * be called from within Sundials to allocate vectors during the nonlinear
    * solution process.  The original solution vector must be setup by the
    * user's application code.
    */
   virtual SundialsAbstractVector *
   makeNewVector() = 0;

   /**
    * Destroy vector structure and its storage. This function is distinct
    * from the destructor since it will be called from within Sundials to
    * deallocate vectors during the nonlinear solution process.
    */
   virtual void
   freeVector() = 0;

   /**
    * Initialize all entries of this vector object to scalar \f$c\f$.
    */
   virtual void
   setToScalar(
      const double c) = 0;

   /**
    * Set this vector object to scalar \f$c x\f$, where \f$c\f$ is a scalar and
    * x is another vector.
    */
   virtual void
   scaleVector(
      const SundialsAbstractVector* x,
      const double c) = 0;

   /**
    * Set this vector object to \f$a x + b y\f$, where \f$a, b\f$ are scalars and
    * \f$x, y\f$ are vectors.
    */
   virtual void
   setLinearSum(
      const double a,
      const SundialsAbstractVector* x,
      const double b,
      const SundialsAbstractVector* y) = 0;

   /**
    * Set each entry of this vector: \f$v_i = x_i y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.
    */
   virtual void
   pointwiseMultiply(
      const SundialsAbstractVector* x,
      const SundialsAbstractVector* y) = 0;

   /**
    * Set each entry of this vector: \f$v_i = x_i / y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.  Based on the Sundials vector
    * implementation, it is not necessary to check for division by zero.
    */
   virtual void
   pointwiseDivide(
      const SundialsAbstractVector* x,
      const SundialsAbstractVector* y) = 0;

   /**
    * Set each entry of this vector to the absolute value of the
    * corresponding entry in vector \f$x\f$.
    */
   virtual void
   setAbs(
      const SundialsAbstractVector* x) = 0;

   /**
    * Set each entry of this vector: \f$v_i =  1 / x_i\f$, where \f$x_i\f$ is an entry
    * entry in vector \f$x\f$.  Based on the Sundials vector implementation,
    * it is not necessary to no check for division by zero.
    */
   virtual void
   pointwiseReciprocal(
      const SundialsAbstractVector* x) = 0;

   /**
    * Set each entry of this vector to the corresponding entry in vector \f$x\f$
    * plus the scalar \f$b\f$.
    */
   virtual void
   addScalar(
      const SundialsAbstractVector* x,
      const double b) = 0;

   /**
    * Return the dot product of this vector and the argument vector \f$x\f$.
    */
   virtual double
   dotWith(
      const SundialsAbstractVector* x) const = 0;

   /**
    * Return the max norm of this vector.
    */
   virtual double
   maxNorm() const = 0;

   /**
    * Return the \f$L_1\f$ norm of this vector.
    */
   virtual double
   L1Norm() const = 0;

   /**
    * Return the weighted-\f$L_2\f$ norm of this vector using the vector \f$x\f$
    * as the weighting vector.
    */
   virtual double
   weightedL2Norm(
      const SundialsAbstractVector* x) const = 0;

   /**
    * Return the weighted root mean squared norm of this vector using
    * the vector \f$x\f$ as the weighting vector.
    */
   virtual double
   weightedRMSNorm(
      const SundialsAbstractVector* x) const = 0;

   /**
    * Return the minimum entry of this vector.
    */
   virtual double
   vecMin() const = 0;

   /**
    * Return \f$0\f$ if \f$x_i \neq 0\f$ and \f$x_i z_i \leq 0\f$, for some \f$i\f$.
    * Here \f$z_i\f$ is an element of this vector. Otherwise, return \f$1\f$.
    */
   virtual int
   constrProdPos(
      const SundialsAbstractVector* x) const = 0;

   /**
    * Set each entry in this vector based on the vector \f$x\f$ as follows:
    * if \f$\mid x_i \mid \geq c\f$, then \f$v_i = 1\f$, else \f$v_i = 0\f$.
    */
   virtual void
   compareToScalar(
      const SundialsAbstractVector* x,
      const double c) = 0;

   /**
    * Set each entry of this vector: \f$v_i =  1 / x_i\f$, where \f$x_i\f$ is an
    * entry entry in the vector \f$x\f$, unless \f$x_i = 0\f$.  If \f$x_i = 0\f$,
    * then return \f$0\f$.  Otherwise, \f$1\f$ is returned.
    */
   virtual int
   testReciprocal(
      const SundialsAbstractVector* x) = 0;

   /*!
    * @brief Get the length of this vector.
    *
    * @return The length (number of elements in the underlying data)
    */
   virtual sunindextype getLength() const = 0;

   /**
    * Return the wrapped Sundials N_Vector.
    */
   N_Vector
   getNVector()
   {
      return d_n_vector;
   }

   /**
    * Print the vector data to the output stream used by the subclass
    * print routine.
    */
   virtual void
   printVector() const = 0;

private:
   N_Vector d_n_vector;

   /**
    * Create Sundials VectorOps structure
    *
    */
   static N_Vector_Ops
   createVectorOps();

   /**
    * The Sundials vector operations
    *
    */

//   static N_Vector N_VCloneEmpty_SAMRAI(N_Vector w);
   static N_Vector
   N_VClone_SAMRAI(
      N_Vector w)
   {
      /* Create content, which in this case is the SAMRAI
       * wrapper vector object */
      SundialsAbstractVector* v = SABSVEC_CAST(w)->makeNewVector();
      return v->getNVector();
   }

   static void
   N_VDestroy_SAMRAI(
      N_Vector v)
   {
      if (v) {
         SABSVEC_CAST(v)->freeVector();
      }
   }

// static void N_VSpace_SAMRAI(N_Vector v, long int *lrw, long int *liw);
// static realtype *N_VGetArrayPointer_SAMRAI(N_Vector v);
// static void N_VSetArrayPointer_SAMRAI(realtype *v_data, N_Vector v);
   static void
   N_VLinearSum_SAMRAI(
      realtype a,
      N_Vector x,
      realtype b,
      N_Vector y,
      N_Vector z)
   {
      SABSVEC_CAST(z)->setLinearSum(a, SABSVEC_CAST(x), b, SABSVEC_CAST(y));
   }

   static void
   N_VConst_SAMRAI(
      realtype c,
      N_Vector z)
   {
      SABSVEC_CAST(z)->setToScalar(c);
   }

   static void
   N_VProd_SAMRAI(
      N_Vector x,
      N_Vector y,
      N_Vector z)
   {
      SABSVEC_CAST(z)->pointwiseMultiply(SABSVEC_CAST(x), SABSVEC_CAST(y));
   }

   static void
   N_VDiv_SAMRAI(
      N_Vector x,
      N_Vector y,
      N_Vector z)
   {
      SABSVEC_CAST(z)->pointwiseDivide(SABSVEC_CAST(x), SABSVEC_CAST(y));
   }

   static void
   N_VScale_SAMRAI(
      realtype c,
      N_Vector x,
      N_Vector z)
   {
      SABSVEC_CAST(z)->scaleVector(SABSVEC_CAST(x), c);
   }

   static void
   N_VAbs_SAMRAI(
      N_Vector x,
      N_Vector z)
   {
      SABSVEC_CAST(z)->setAbs(SABSVEC_CAST(x));
   }

   static void
   N_VInv_SAMRAI(
      N_Vector x,
      N_Vector z)
   {
      SABSVEC_CAST(z)->pointwiseReciprocal(SABSVEC_CAST(x));
   }

   static void
   N_VAddConst_SAMRAI(
      N_Vector x,
      realtype b,
      N_Vector z)
   {
      SABSVEC_CAST(z)->addScalar(SABSVEC_CAST(x), b);
   }

   static realtype
   N_VDotProd_SAMRAI(
      N_Vector x,
      N_Vector y)
   {
      return SABSVEC_CAST(x)->dotWith(SABSVEC_CAST(y));
   }

   static realtype
   N_VMaxNorm_SAMRAI(
      N_Vector x)
   {
      return SABSVEC_CAST(x)->maxNorm();
   }

   static realtype
   N_VWrmsNorm_SAMRAI(
      N_Vector x,
      N_Vector w)
   {
      return SABSVEC_CAST(x)->weightedRMSNorm(SABSVEC_CAST(w));
   }

// static realtype N_VWrmsNormMask_SAMRAI(N_Vector x, N_Vector w, N_Vector id);
   static realtype
   N_VMin_SAMRAI(
      N_Vector x)
   {
      return SABSVEC_CAST(x)->vecMin();
   }

   static realtype
   N_VWL2Norm_SAMRAI(
      N_Vector x,
      N_Vector w)
   {
      return SABSVEC_CAST(x)->weightedL2Norm(SABSVEC_CAST(w));
   }

   static realtype
   N_VL1Norm_SAMRAI(
      N_Vector x)
   {
      return SABSVEC_CAST(x)->L1Norm();
   }

   static void
   N_VCompare_SAMRAI(
      realtype c,
      N_Vector x,
      N_Vector z)
   {
      SABSVEC_CAST(z)->compareToScalar(SABSVEC_CAST(x), c);
   }

   static booleantype
   N_VInvTest_SAMRAI(
      N_Vector x,
      N_Vector z)
   {
      return SABSVEC_CAST(z)->testReciprocal(SABSVEC_CAST(x));
   }

   static booleantype
   N_VConstrMask_SAMRAI(
      N_Vector c,
      N_Vector x,
      N_Vector m);

   static realtype
   N_VMinQuotient_SAMRAI(
      N_Vector num,
      N_Vector denom);

   static sunindextype
   N_VGetLength_SAMRAI(N_Vector x)
   {
      return SABSVEC_CAST(x)->getLength();
   }

};

}
}

#endif
#endif
