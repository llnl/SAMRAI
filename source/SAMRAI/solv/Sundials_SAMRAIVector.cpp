/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   "Glue code" between SAMRAI vector object and Sundials vector.
 *
 ************************************************************************/
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"

#ifdef HAVE_SUNDIALS

#define SKVEC_CAST(x) (((Sundials_SAMRAIVector *)x))

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Static public member functions.
 *
 *************************************************************************
 */

SundialsAbstractVector *
Sundials_SAMRAIVector::createSundialsVector(
   const std::shared_ptr<SAMRAIVectorReal<double> >& samrai_vec)
{
   TBOX_ASSERT(samrai_vec);
   SundialsAbstractVector* skv = new Sundials_SAMRAIVector(samrai_vec);

   return skv;
}

void
Sundials_SAMRAIVector::destroySundialsVector(
   SundialsAbstractVector* sundials_vec)
{
   if (sundials_vec) {
      delete (dynamic_cast<Sundials_SAMRAIVector *>(sundials_vec));
   }
}

std::shared_ptr<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector(
   SundialsAbstractVector* sundials_vec)
{
   TBOX_ASSERT(sundials_vec != 0);
   return (dynamic_cast<Sundials_SAMRAIVector *>(sundials_vec))->
          getSAMRAIVector();
}

std::shared_ptr<SAMRAIVectorReal<double> >
Sundials_SAMRAIVector::getSAMRAIVector(
   N_Vector sundials_vec)
{
   TBOX_ASSERT(sundials_vec != 0);
// sgs
   return static_cast<Sundials_SAMRAIVector *>(sundials_vec->content)->
          getSAMRAIVector();
}

/*
 *************************************************************************
 *
 * Constructor and destructor for Sundials_SAMRAIVector.
 *
 *************************************************************************
 */

Sundials_SAMRAIVector::Sundials_SAMRAIVector(
   const std::shared_ptr<SAMRAIVectorReal<double> >& samrai_vector):
   SundialsAbstractVector(),
   d_samrai_vector(samrai_vector)
{
}

Sundials_SAMRAIVector::~Sundials_SAMRAIVector()
{
}

/*
 *************************************************************************
 *
 * Other miscellaneous member functions
 *
 *************************************************************************
 */

SundialsAbstractVector *
Sundials_SAMRAIVector::makeNewVector()
{
   Sundials_SAMRAIVector* out_vec =
      new Sundials_SAMRAIVector(d_samrai_vector->cloneVector("out_vec"));
   out_vec->getSAMRAIVector()->allocateVectorData();
   return out_vec;
}

void
Sundials_SAMRAIVector::freeVector()
{
   d_samrai_vector->freeVectorComponents();
   d_samrai_vector.reset();
   delete this;
}

void
Sundials_SAMRAIVector::printVector() const
{
   std::ostream& s = d_samrai_vector->getOutputStream();
   s << "\nPrinting Sundials_SAMRAIVector..."
     << "\nthis = " << (Sundials_SAMRAIVector *)this
     << "\nSAMRAI vector structure and data: " << std::endl;
   d_samrai_vector->print(s);
   s << "\n" << std::endl;
}

void
Sundials_SAMRAIVector::setToScalar(
   const double c)
{
   d_samrai_vector->setToScalar(c);
}

void
Sundials_SAMRAIVector::scaleVector(
   const SundialsAbstractVector* x,
   const double c)
{
   d_samrai_vector->scale(c, SKVEC_CAST(x)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::setLinearSum(
   const double a,
   const SundialsAbstractVector* x,
   const double b,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->linearSum(a, SKVEC_CAST(x)->getSAMRAIVector(),
      b, SKVEC_CAST(y)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::pointwiseMultiply(
   const SundialsAbstractVector* x,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->multiply(SKVEC_CAST(x)->getSAMRAIVector(),
      SKVEC_CAST(y)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::pointwiseDivide(
   const SundialsAbstractVector* x,
   const SundialsAbstractVector* y)
{
   d_samrai_vector->divide(SKVEC_CAST(x)->getSAMRAIVector(),
      SKVEC_CAST(y)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::setAbs(
   const SundialsAbstractVector* x)
{
   d_samrai_vector->abs(SKVEC_CAST(x)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::pointwiseReciprocal(
   const SundialsAbstractVector* x)
{
   d_samrai_vector->reciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::addScalar(
   const SundialsAbstractVector* x,
   const double b)
{
   d_samrai_vector->addScalar(SKVEC_CAST(x)->getSAMRAIVector(), b);
}

double
Sundials_SAMRAIVector::dotWith(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->dot(SKVEC_CAST(x)->getSAMRAIVector());
}

double
Sundials_SAMRAIVector::maxNorm() const
{
   return d_samrai_vector->maxNorm();
}

double
Sundials_SAMRAIVector::L1Norm() const
{
   return d_samrai_vector->L1Norm();
}

double
Sundials_SAMRAIVector::weightedL2Norm(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->weightedL2Norm(SKVEC_CAST(x)->getSAMRAIVector());
}

double
Sundials_SAMRAIVector::weightedRMSNorm(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->weightedRMSNorm(SKVEC_CAST(x)->getSAMRAIVector());
}

double
Sundials_SAMRAIVector::vecMin() const
{
   return d_samrai_vector->min();
}

int
Sundials_SAMRAIVector::constrProdPos(
   const SundialsAbstractVector* x) const
{
   return d_samrai_vector->
          computeConstrProdPos(SKVEC_CAST(x)->getSAMRAIVector());
}

void
Sundials_SAMRAIVector::compareToScalar(
   const SundialsAbstractVector* x,
   const double c)
{
   d_samrai_vector->compareToScalar(SKVEC_CAST(x)->getSAMRAIVector(), c);
}

int
Sundials_SAMRAIVector::testReciprocal(
   const SundialsAbstractVector* x)
{
   return d_samrai_vector->testReciprocal(SKVEC_CAST(x)->getSAMRAIVector());
}

sunindextype
Sundials_SAMRAIVector::getLength() const
{
   return static_cast<sunindextype>(d_samrai_vector->getLength());
}

}
}

#endif
