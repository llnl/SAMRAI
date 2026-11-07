/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between Sundials vector interface and SAMRAI vectors.
 *
 ************************************************************************/

#ifndef included_solv_Sundials_SAMRAIVector
#define included_solv_Sundials_SAMRAIVector

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT
 *  KINSOL -or- CVODE
 ************************************************************************
 */
#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"

#include <memory>

namespace SAMRAI {
namespace solv {

/**
 * Class Sundials_SAMRAIVector wraps a real-valued SAMRAI vector
 * (see SAMRAIVectorReal class) object so that it may be used with
 * the Sundials solver packages.  This class is derived from the
 * abstract base class SundialsAbstractVector, which defines a <TT>C++</TT>
 * interface for Sundials vectors.  It also maintains a pointer to a SAMRAI
 * vector object.  A SAMRAI vector is defined as a collection of patch data
 * components living on some subset of levels in a structured AMR mesh
 * hierarchy.
 *
 * Observe that there are only three public member functions in this class
 * They are used to create and destroy Sundials vector objects (i.e.,
 * "N_Vector"s), and to obtain the SAMRAI vector associated with the Sundials
 * vector.  In particular, note that the constructor and destructor of this
 * class are protected members.  The construction and destruction of instances
 * of this class may occur only through the static member functions that
 * create and destroy Sundials vector objects.
 *
 * Finally, we remark that this class provides vectors of type <TT>double</TT>,
 * which is the default for Sundials.
 *
 * @see SundialsAbstractVector
 * @see SAMRAIVectorReal
 */

class Sundials_SAMRAIVector:public SundialsAbstractVector
{
public:
   /**
    * Create and return a new SundialsAbstractVector vector object.  The
    * SAMRAI vector object is wrapped so that it may be manipulated
    * within Sundials as an N_Vector (which is typedef'd to
    * SundialsAbstractVector* in the abstract Sundials vector
    * interface).  It is important to note that this function does not
    * allocate storage for the vector data.  Data must be allocated
    * through the SAMRAI vector object directly.  For output of the
    * data through "N_VPrint" calls, the output stream to which the
    * SAMRAI vector object writes will be used.
    *
    * @pre samrai_vec
    */
   static SundialsAbstractVector *
   createSundialsVector(
      const std::shared_ptr<SAMRAIVectorReal<double> >& samrai_vec);

   /**
    * Destroy a given Sundials vector object. It is important to note that
    * this function does not deallocate storage for the vector data.
    * Vector data must be deallocated hrough the SAMRAI vector object.
    */
   static void
   destroySundialsVector(
      SundialsAbstractVector* sundials_vec);

   /**
    * Return pointer to the SAMRAI vector object associated with the
    * given Sundials wrapper vector.
    *
    * @pre sundials_vec != 0
    */
   static std::shared_ptr<SAMRAIVectorReal<double> >
   getSAMRAIVector(
      SundialsAbstractVector* sundials_vec);

   /**
    * Return pointer to the SAMRAI vector object associated with the
    * given Sundials vector.
    *
    * @pre sundials_vec != 0
    */
   static std::shared_ptr<SAMRAIVectorReal<double> >
   getSAMRAIVector(
      N_Vector sundials_vec);

   /*
    * Print the vector to the output stream used by the SAMRAI vector class.
    */
   void
   printVector() const;

protected:
   /*
    * Constructor for Sundials_SAMRAIVector.
    */
   explicit Sundials_SAMRAIVector(
      const std::shared_ptr<SAMRAIVectorReal<double> >& samrai_vector);

   /*
    * Virtual destructor for Sundials_SAMRAIVector.
    */
   virtual ~Sundials_SAMRAIVector();

private:
   /*
    * Return SAMRAI vector owned by this Sundials_SAMRAIVector object.
    */
   std::shared_ptr<SAMRAIVectorReal<double> >
   getSAMRAIVector()
   {
      return d_samrai_vector;
   }

   /*
    * The makeNewVector() function clones the vector structure and allocate
    * storage for the new vector data.  Then, a pointer to the new vector
    * instance is returned. This function is distinct from the constructor
    * since it will be called from within Sundials to allocate vectors used
    * during the nonlinear solution process.
    */
   SundialsAbstractVector *
   makeNewVector();

   /*
    * Destroy vector structure and its storage. This function is distinct
    * from the destructor since it will be called from within Sundials to
    * deallocate vectors during the nonlinear solution process.
    */
   virtual void
   freeVector();

   /*
    * Initialize all entries of this vector object to scalar \f$c\f$.
    */
   void
   setToScalar(
      const double c);

   /*
    * Set this vector object to scalar \f$c x\f$, where \f$c\f$ is a scalar and
    * \f$x\f$ is another vector.
    */
   void
   scaleVector(
      const SundialsAbstractVector* x,
      const double c);

   /*
    * Set this vector object to \f$a x + b y\f$, where \f$a, b\f$ are scalars and
    * \f$x, y\f$ are vectors.
    */
   void
   setLinearSum(
      const double a,
      const SundialsAbstractVector* x,
      const double b,
      const SundialsAbstractVector* y);

   /*
    * Set each entry of this vector: \f$v_i = x_i y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.
    */
   void
   pointwiseMultiply(
      const SundialsAbstractVector* x,
      const SundialsAbstractVector* y);

   /*
    * Set each entry of this vector: \f$v_i = \frac{x_i}{y_i}\f$, where
    * \f$x_i, y_i\f$ are entries in vectors \f$x\f$ and \f$y\f$.  Based on the Sundials
    * vector implementation, it is not necessary to no check for division by
    * zero.
    */
   void
   pointwiseDivide(
      const SundialsAbstractVector* x,
      const SundialsAbstractVector* y);

   /*
    * Set each entry of this vector to the absolute value of the
    * corresponding entry in vector \f$x\f$.
    */
   void
   setAbs(
      const SundialsAbstractVector* x);

   /*
    * Set each entry of this vector: \f$v_i = \frac{1}{x_i}\f$, where \f$x_i\f$ is
    * an entry in vector \f$x\f$.  Based on the Sundials vector implementation,
    * it is not necessary to no check for division by zero.
    */
   void
   pointwiseReciprocal(
      const SundialsAbstractVector* x);

   /*
    * Set each entry of this vector: \f$v_i = x_i + b\f$, where \f$x_i\f$ is an entry
    * in the vector \f$x\f$ and \f$b\f$ is a scalar.
    */
   void
   addScalar(
      const SundialsAbstractVector* x,
      const double b);

   /*
    * Return the dot product of this vector and the vector \f$x\f$.
    */
   double
   dotWith(
      const SundialsAbstractVector* x) const;

   /*
    * Return the max norm of this vector:
    * \f${\| v \|}_{\max} = \max_{i} (\mid v_i \mid)\f$.
    */
   double
   maxNorm() const;

   /*
    * Return the \f$L_1\f$ norm of this vector:
    * \f${\| v \|}_{L_1} = \sum_{i} (\mid v_i \mid)\f$ if no control volumes
    * are defined.  Otherwise,
    * \f${\| v \|}_{L_1} = \sum_{i} (\mid v_i \mid * cvol_i )\f$.
    */
   double
   L1Norm() const;

   /*
    * Return the weighted \f$L_2\f$ norm of this vector using the vector
    * \f$x\f$ as the weighting vector:
    * \f${\| v \|}_{WL2(x)} = \sqrt{ \sum_i( (x_i * v_i)^2 ) )}\f$ if no
    * control volumes are defined.  Otherwise,
    * \f${\| v \|}_{WL2(x)} = \sqrt{ \sum_i( (x_i * v_i)^2 cvol_i ) }\f$.
    */
   double
   weightedL2Norm(
      const SundialsAbstractVector* x) const;

   /*
    * Return the weighted root mean squared norm of this vector using
    * the vector \f$x\f$ as the weighting vector.   If control volumes are
    * not defined for the vector entries, the norm corresponds to the
    * weighted \f$L_2\f$-norm divided by the square root of the number of
    * vector entries.  Otherwise, the norm corresponds to the weighted
    * \f$L_2\f$-norm divided by the square root of the sum of the control volumes.
    */
   double
   weightedRMSNorm(
      const SundialsAbstractVector* x) const;

   /*
    * Return the minimum entry of this vector.
    */
   double
   vecMin() const;

   /*
    * Return \f$0\f$ if \f$x_i \neq 0\f$ and \f$x_i v_i \leq 0\f$, for some \f$i\f$.
    * Otherwise, return \f$1\f$.
    */
   int
   constrProdPos(
      const SundialsAbstractVector* x) const;

   /*
    * Set each entry in this vector based on the vector \f$x\f$ as follows:
    * if \f$\mid x_i \mid \geq c\f$, then \f$v_i = 1\f$, else \f$v_i = 0\f$.
    */
   void
   compareToScalar(
      const SundialsAbstractVector* x,
      const double c);

   /*
    * Set each entry of this vector: \f$v_i =  \frac{1}{x_i}\f$, where \f$x_i\f$
    * is an entry entry in the vector \f$x\f$, unless \f$x_i = 0\f$.  If \f$x_i = 0\f$,
    * then return \f$0\f$.  Otherwise, \f$1\f$ is returned.
    */
   int
   testReciprocal(
      const SundialsAbstractVector* x);

   /*!
    * @brief Get the length of this vector.
    *
    * @return The length (number of elements in the underlying data)
    */
   sunindextype getLength() const;

   /*
    * Vector data is maintained in SAMRAI vector structure.
    */
   std::shared_ptr<SAMRAIVectorReal<double> > d_samrai_vector;

};

}
}

#endif
#endif
