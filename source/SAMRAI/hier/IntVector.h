/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A N-dimensional integer vector
 *
 ************************************************************************/

#ifndef included_hier_IntVector
#define included_hier_IntVector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <array>
#include <vector>
#include <iostream>
#include <optional>

namespace SAMRAI {

namespace hier {

class Index;

/*!
 * @brief Simple integer vector class with size based on a dimension value
 *
 * Class IntVector implements a vector of integers that, depending on the usage
 * context, has a length equal to the number of spatial dimensions being used
 * or its length is the number of dimension multiplied by a certain number of
 * blocks.
 *
 * The integers in IntVector are stored first in order of dimension, and then
 * in order of block number.  For a 3-dimensional IntVector on 4 blocks, the
 * order would be: { i0, j0, k0, i1, j1, k1, i2, j2, k2, i3, j3, k3 }.
 *
 * When used in the context of a single-block problem, the number of blocks
 * for an IntVector is always 1.  In a multiblock context, the number of
 * blocks for an IntVector may be either 1 or the number of blocks being used
 * in the problem.
 */

class IntVector
{
public:
   typedef tbox::Dimension::dir_t dir_t;

   /*!
    * @brief Creates an uninitialized IntVector for 1 block.
    *
    * @param dim
    */
   explicit IntVector(
      const tbox::Dimension& dim);

   /*!
    * @brief Creates an uninitialized IntVector of a given number of blocks.
    *
    * @pre num_blocks >=1
    *
    * @param num_blocks
    * @param dim
    */
   IntVector(
      size_t num_blocks,
      const tbox::Dimension& dim);

   /*!
    * @brief Construct an IntVector with all components equal to the
    * value argument.
    *
    * @pre num_blocks >=1
    *
    * @param dim
    * @param value
    * @param num_blocks
    */
   IntVector(
      const tbox::Dimension& dim,
      int value,
      size_t num_blocks = 1);

   /*!
    * @brief Construct an IntVector with the values provided by
    * an STL vector of ints.
    *
    * The dimension of the constructed IntVector will be the size of the
    * vec argument.  If num_blocks has a value greater than 1, then the
    * IntVector will be constructed with the values held by vec duplicated for
    * every block.
    *
    * @pre vec.size() >= 1
    *
    * @param vec Vector of integers with a size equal to the desired dimension
    * @param num_blocks
    */
   explicit IntVector(
      const std::vector<int>& vec,
      size_t num_blocks = 1);

   /*!
    * @brief Construct an IntVector with values provided by a raw array.
    *
    * This constructor assumes that the given array contains a number of values
    * equal to the dimension value.  As this constructor can do no error-
    * checking that the array is properly allocated and initialized, it is
    * up to the calling code to ensure that the array argument is valid.
    *
    * If num_blocks has a value greater than 1, then the IntVector
    * will be constructed with the values held by array duplicated for every
    * block.
    *
    * @param dim
    * @param array  Array of ints that should be allocated and initialized
    *               at a length equal to dim.getValue()
    * @param num_blocks
    */
   IntVector(
      const tbox::Dimension& dim,
      const int array[],
      size_t num_blocks = 1);

   /*!
    * @brief Copy constructor.
    */
   IntVector(const IntVector& rhs) = default;

   /*!
    * @brief Move copy constructor.
    */
   IntVector(IntVector&& rhs) = default;

   /*!
    * @brief Construct an IntVector from another IntVector.
    *
    * The main use case for this constructor is to use an IntVector sized for
    * one block to construct an IntVector sized for a larger number of blocks.
    * When used in this way, the constructed IntVector will duplicate the
    * contents of the argument for every block.
    *
    * If num_blocks is equal to the rhs argument's number of blocks, then this
    * constructor is equivalent to the copy constructor.
    *
    * @pre num_blocks >=1
    * @pre (rhs.getNumBlocks() == num_blocks || rhs.getNumBlocks() == 1)
    *
    * @param rhs
    * @param num_blocks
    */
   IntVector(
      const IntVector& rhs,
      size_t num_blocks);

   /*!
    * @brief Construct an IntVector from an Index.
    *
    * The constructed IntVector will have the same dimension value as the
    * Index.  If num_blocks is greater than 1, the values held by the Index
    * argument will be duplicated for every block.
    *
    * @param rhs
    * @param num_blocks
    *
    * @pre num_blocks >=1
    *
    * @param rhs
    * @param num_blocks
    */
   explicit IntVector(
      const Index& rhs,
      size_t num_blocks = 1);

   /*!
    * @brief The assignment operator sets the IntVector equal to the
    *        argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector&
   operator = (
      const IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      d_num_blocks = rhs.d_num_blocks;
      d_vector = rhs.d_vector;

      return *this;
   }

   /*!
    * @brief Move assignment operator
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector& operator=(IntVector&& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      d_num_blocks = rhs.d_num_blocks;
      d_vector = std::move(rhs.d_vector);
      return *this;
   }

   /*!
    * @brief Assignment operator assigning the values of an Index to an
    * IntVector.
    *
    * The assigned IntVector will be sized for one block.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector&
   operator = (
      const Index& rhs);

   /*!
    * @brief default destructor
    */
   ~IntVector() noexcept = default;

   /*!
    * @brief Return the number of blocks for this IntVector
    */
   size_t getNumBlocks() const noexcept
   {
      return d_num_blocks;
   }

   /*!
    * @brief Return an IntVector for one block extracted from a possibly
    * larger IntVector
    *
    * This constructs an IntVector sized for one block using the int values
    * associated with a given block.  The constructed IntVector is returned
    * by value.
    *
    * @pre block_id.getBlockValue() < getNumBlocks()
    *
    * @param block_id  BlockId indicates which block is associated with
    *                  the desired integer values.
    *
    * @return A constructed IntVector sized for one block
    */
   IntVector getBlockVector(const BlockId& block_id) const
   {
      TBOX_ASSERT(block_id.getBlockValue() < d_num_blocks);
      IntVector block_vec(d_dim,
                          &(d_vector[block_id.getBlockValue()*d_dim.getValue()]));

      return block_vec;
   }

   /*!
    * @brief Return reference to the specified component of the vector.  This
    * can only be used when the number of blocks is 1.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getNumBlocks() == 1
    */
   int&
   operator [] (
      const unsigned int i) noexcept
   {
      TBOX_ASSERT(i < d_dim.getValue());
      TBOX_ASSERT(d_num_blocks == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference.  This can only be used when the number of blocks is 1.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getNumBlocks() == 1
    */
   const int&
   operator [] (
      const unsigned int i) const noexcept
   {
      TBOX_ASSERT(i < d_dim.getValue());
      TBOX_ASSERT(d_num_blocks == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return reference to the specified component of the vector.  This
    * can only be used when the number of blocks is 1.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getNumBlocks() == 1
    */
   int&
   operator () (
      const unsigned int i) noexcept
   {
      TBOX_ASSERT(i < d_dim.getValue());
      TBOX_ASSERT(d_num_blocks == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference.  This can only be used when the number of blocks is 1.
    *
    * @pre (i >= 0) && (i < getDim().getValue())
    * @pre getNumBlocks() == 1
    */
   const int&
   operator () (
      const unsigned int i) const noexcept
   {
      TBOX_ASSERT(i < d_dim.getValue());
      TBOX_ASSERT(d_num_blocks == 1);
      return d_vector[i];
   }

   /*!
    * @brief Return the specified component of the vector.
    *
    * @pre (b >= 0) && (b < getNumBlocks())
    * @pre (i >= 0) && (i < getDim().getValue())
    *
    * The desired component is specified by the pair of the block number b
    * and the dimensional index i.
    */
   int&
   operator () (
      const BlockId::block_t b,
      const unsigned int i) noexcept
   {
      TBOX_ASSERT(b < d_num_blocks);
      TBOX_ASSERT(i < d_dim.getValue());
      return d_vector[b*d_dim.getValue() + i];
   }

   /*!
    * @brief Return the specified component of the vector as a const integer
    * reference
    *
    * @pre (b >= 0) && (b < getNumBlocks())
    * @pre (i >= 0) && (i < getDim().getValue())
    *
    * The desired component is specified by the pair of the block number b
    * and the dimensional index i.
    */
   const int&
   operator () (
      const BlockId::block_t b,
      const unsigned int i) const noexcept
   {
      TBOX_ASSERT(b < d_num_blocks);
      TBOX_ASSERT(i < d_dim.getValue());
      return d_vector[b*d_dim.getValue() + i];
   }

   /*!
    * @brief Plus-equals operator for two integer vectors.
    *
    * The following comments apply to most of the arithmetic operators and
    * some related methods:
    *
    * The arithmetic operators and related methods such as ceilingDivide
    * allow both full component-wise operations on IntVectors that have equal
    * numbers of blocks, as well as operations between multiblock IntVectors and
    * single-block IntVectors
    *
    * When an operator is invoked on a multiblock IntVector with a right-hand-
    * side single-block IntVector, the values within the right-hand-side are
    * applied to every block within the "this" IntVector.
    *
    * Assertion failures occur if operations are called using two multiblock
    * IntVectors that have different numbers of blocks, or if the
    * right-hand-side is a multiblock IntVector while "this" is single-block.
    *
    * There also are some operators that take a single integer as a right-
    * hand-side, in which case that integer value is applied in the operation
    * to every component of "this".
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   IntVector&
   operator += (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] += rhs.d_vector[i];
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            d_vector[i] += rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Plus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator + (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Plus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator += (
      const int rhs) noexcept
   {
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         d_vector[i] += rhs;
      }
      return *this;
   }

   /*!
    * @brief Plus operator for an integer vector and an integer.
    */
   IntVector
   operator + (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    *
    */
   IntVector&
   operator -= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] -= rhs.d_vector[i];
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            d_vector[i] -= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Minus operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator - (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator -= (
      const int rhs) noexcept
   {
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         d_vector[i] -= rhs;
      }
      return *this;
   }

   /*!
    * @brief Minus operator for an integer vector and an integer.
    */
   IntVector
   operator - (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Times-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    *
    */
   IntVector&
   operator *= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] *= rhs.d_vector[i];
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            d_vector[i] *= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Times operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator * (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp *= rhs;
      return tmp;
   }

   /*!
    * @brief Times-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator *= (
      const int rhs) noexcept
   {
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         d_vector[i] *= rhs;
      }
      return *this;
   }

   /*!
    * @brief Times operator for an integer vector and an integer.
    */
   IntVector
   operator * (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp *= rhs;
      return tmp;
   }

   /*!
    * @brief Quotient-equals operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   IntVector&
   operator /= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] /= rhs.d_vector[i];
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            d_vector[i] /= rhs.d_vector[i];
         }
      }
      return *this;
   }

   /*!
    * @brief Quotient operator for two integer vectors.
    *
    * @pre getDim() == rhs.getDim()
    */
   IntVector
   operator / (
      const IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      IntVector tmp(*this);
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Quotient-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator /= (
      const int rhs) noexcept
   {
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         d_vector[i] /= rhs;
      }
      return *this;
   }

   /*!
    * @brief Quotient operator for an integer vector and an integer.
    */
   IntVector
   operator / (
      const int rhs) const
   {
      IntVector tmp(*this);
      tmp /= rhs;
      return tmp;
   }


   /*!
    * @brief Component-wise ceilingDivide quotient (integer divide with
    *        rounding up).
    *
    * @pre getDim() == denominator.getDim()
    * @pre getNumBlocks() == denominator.getNumBlocks() || denominator.getNumBlocks() == 1
    */
   void
   ceilingDivide(
      const IntVector& denominator) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, denominator);
      TBOX_ASSERT(d_num_blocks == denominator.d_num_blocks ||
                  denominator.d_num_blocks == 1);

      /*
       * This is the formula for integer divide, rounding away from
       * zero.  It is meant as an extension of the ceilingDivide quotient of
       * 2 positive integers.
       *
       * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
       * the results of rounding.
       * - Add zero if there's no remainder to round.
       * - Round remainder to 1 if numerator and denominator has same sign.
       * - Round remainder to -1 if numerator and denominator has opposite sign.
       */
      if (denominator.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               d_vector[offset + i] = (d_vector[offset + i] / denominator[i]) +
               ((d_vector[offset + i] % denominator[i]) ?
                  ((d_vector[offset + i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            d_vector[i] = (d_vector[i] / denominator.d_vector[i]) +
               ((d_vector[i] % denominator.d_vector[i]) ?
               ((d_vector[i] > 0) == (denominator.d_vector[i] > 0) ? 1 : -1) : 0);
         }
      }
   }

   /*!
    * @brief Component-wise ceilingDivide quotient (integer divide with
    *        rounding up).
    *
    * @pre numerator.getDim() == denominator.getDim()
    * @pre numerator.getNumBlocks() == denominator.getNumBlocks() || denominator.getNumBlocks() == 1
    */
   static IntVector
   ceilingDivide(
      const IntVector& numerator,
      const IntVector& denominator)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(numerator, denominator);
      TBOX_ASSERT(numerator.d_num_blocks == denominator.d_num_blocks ||
                  denominator.d_num_blocks == 1);
      IntVector rval(numerator);
      rval.ceilingDivide(denominator);
      return rval;
   }

   /*!
    * @brief Unary minus to negate an integer vector.
    */
   IntVector
   operator - () const
   {
      IntVector tmp(*this);
      tmp *= -1;
      return tmp;
   }

   /*!
    * @brief Returns true if all components are equal to a given integer.
    */
   bool
   operator == (
      int rhs) const noexcept
   {
      bool result = true;
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; result && (i < length); ++i) {
         result = d_vector[i] == rhs;
      }
      return result;
   }

   /*!
    * @brief Returns true if some components are not equal to a given integer.
    */
   bool
   operator != (
      int rhs) const noexcept
   {
      return !(*this == rhs);
   }

   /*!
    * @brief Returns true if two vector objects are equal.  All components
    *        must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   bool
   operator == (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      bool result = true;
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; result && b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] == rhs.d_vector[i]);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] == rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if two vector objects are not equal.  Any of
    *        the components may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return !(*this == rhs);
   }

   /*!
    * @brief Returns true if each integer in vector is less than
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   bool
   operator < (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      bool result = true;
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; result && b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] < rhs.d_vector[i]);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] < rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is less or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   bool
   operator <= (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      bool result = true;
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; result && b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] <= rhs.d_vector[i]);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] <= rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is greater than
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   bool
   operator > (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      bool result = true;
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; result && b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] > rhs.d_vector[i]);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] > rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Returns true if each integer in vector is greater or equal to
    *        corresponding integer in comparison vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   bool
   operator >= (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      bool result = true;
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; result && b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; result && (i < d_dim.getValue()); ++i) {
               result = result && (d_vector[offset + i] >= rhs.d_vector[i]);
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; result && (i < length); ++i) {
            result = result && (d_vector[i] >= rhs.d_vector[i]);
         }
      }
      return result;
   }

   /*!
    * @brief Return the component-wise minimum of two integer vector objects.
    *
    * @pre getDim() == rhs.getDim()
    * @pre getNumBlocks() == rhs.getNumBlocks() || rhs.getNumBlocks() == 1
    */
   void
   min(
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               if (rhs.d_vector[i] < d_vector[offset + i]) {
                  d_vector[offset + i] = rhs.d_vector[i];
               }
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            if (rhs.d_vector[i] < d_vector[i]) {
               d_vector[i] = rhs.d_vector[i];
            }
         }
      }
   }

   /*!
    * @brief Return the minimum entry in an integer vector.
    */
   int
   min() const noexcept
   {
      int min = d_vector[0];
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         if (d_vector[i] < min) {
            min = d_vector[i];
         }
      }
      return min;
   }

   /*!
    * @brief Return the component-wise maximum of two integer vector objects.
    */
   void
   max(
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(d_num_blocks == rhs.d_num_blocks || rhs.d_num_blocks == 1);
      if (rhs.d_num_blocks == 1 && d_num_blocks != 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
               if (rhs.d_vector[i] > d_vector[offset + i]) {
                  d_vector[offset + i] = rhs.d_vector[i];
               }
            }
         }
      } else {
         size_t length = d_num_blocks * d_dim.getValue();
         for (unsigned int i = 0; i < length; ++i) {
            if (rhs.d_vector[i] > d_vector[i]) {
               d_vector[i] = rhs.d_vector[i];
            }
         }
      }
   }

   /*!
    * @brief Return the maximum entry in an integer vector.
    */
   int
   max() const noexcept
   {
      int max = d_vector[0];
      size_t length = d_num_blocks * d_dim.getValue();
      for (unsigned int i = 0; i < length; ++i) {
         if (d_vector[i] > max) {
            max = d_vector[i];
         }
      }
      return max;
   }

   /*!
    * @brief Utility function to take the minimum of two integer vector
    *        objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static IntVector
   min(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      IntVector tmp(a);
      tmp.min(b);
      return tmp;
   }

   /*!
    * @brief Utility function to take the maximum of two integer vector
    *        objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static IntVector
   max(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      IntVector tmp(a);
      tmp.max(b);
      return tmp;
   }

   /*!
    * @brief Set all block-wise components of an IntVector.
    *
    * If this IntVector and the argument IntVector are sized for the same
    * number of blocks, this is the equivalent of a copy operation.  If this
    * IntVector is multiblock while the argument is single-block, then
    * the values of the argument are copied to each block-wise component of
    * this IntVector.
    *
    * An error will occur if the numbers of blocks are unequal and the argument
    * IntVector is not single-block.
    *
    * @param vector  Input IntVector
    */
   void setAll(const IntVector& vector)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, vector);
      if (vector.d_num_blocks == d_num_blocks) {
         *this = vector;
      } else if (d_num_blocks > 1 && vector.d_num_blocks == 1) {
         for (BlockId::block_t b = 0; b < d_num_blocks; ++b) {
            unsigned int offset = b*d_dim.getValue();
            for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
               d_vector[offset + d] = vector[d];
            }
         }
      } else {
         TBOX_ERROR("IntVector::setAll() attempted with argument of non-compatible num_blocks." << std::endl);
      }
   }

   /*!
    * @brief Return the product of the entries in the integer vector
    *
    * If a BlockId argument is provided, the product of entries for that
    * block will be computed.  If no BlockId argument is provided, this
    * IntVector must be single-block.
    *
    * @param block_id  Optional block on which to compute the product.
    */
   long int
   getProduct(const BlockId& block_id = BlockId::invalidId()) const noexcept
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(block_id == BlockId::invalidId() ||
                  block_id.getBlockValue() < d_num_blocks);
      if (block_id == BlockId::invalidId()) {
         TBOX_ASSERT(d_num_blocks == 1);
      }
#endif
      BlockId::block_t b = block_id == BlockId::invalidId() ? 0 : block_id.getBlockValue();
      unsigned int offset = b * d_dim.getValue();
      long int prod = 1;
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         prod *= d_vector[offset + i];
      }
      return prod;
   }

   /*!
    * @brief Store the object state to the specified restart database
    *        with the provided name.
    *
    */
   void
   putToRestart(
      tbox::Database& restart_db,
      const std::string& name) const;

   /*!
    * @brief Restores the object giving it the provided name and getting its
    *        state from the specified restart database.
    *
    */
   void
   getFromRestart(
      tbox::Database& restart_db,
      const std::string& name);

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const noexcept
   {
      return d_dim;
   }

   /*!
    * @brief Return an IntVector of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.  The
    * returned IntVector is sized for one block.
    */
   static const IntVector&
   getZero(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_zeros[dim.getValue() - 1].has_value());
      return (s_zeros[dim.getValue() - 1].value());
   }

   /*!
    * @brief Return an IntVector of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.  The
    * returned IntVector is sized for one block.
    */
   static const IntVector&
   getOne(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_ones[dim.getValue() - 1].has_value());
      return (s_ones[dim.getValue() - 1].value());
   }

   /*!
    * @brief Set this IntVector to a sorted version of the given IntVector
    *
    * For an single-block IntVector, set the ith entry of this to the
    * position of the ith smallest value in the argument IntVector.
    *
    * If the IntVectors are multiblock, each section of the IntVector
    * associated with a block is sorted independently.
    *
    * Any values existing in this IntVector will be overwritten by this
    * method.  If this IntVector has a different block size than the input
    * vector, this IntVector will be resized to the input vector's size.
    *
    * @pre d_dim == values.getDim()
    *
    * @param values input IntVector to be sorted
    */
   void
   sortIntVector(
      const IntVector& values);

   /*!
    * @brief Read an integer vector from an input stream.  The format for
    *        the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      IntVector& rhs);

   /*!
    * @brief Write an integer vector into an output stream.  The format for
    *        the output is (i0,...,in) for an n-dimensional vector.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const IntVector& rhs);

private:

   /*
    * Unimplemented default constructor
    */
   IntVector() = delete;

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    *
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   tbox::Dimension d_dim;

   size_t d_num_blocks;

   std::vector<int> d_vector;

   static std::array<std::optional<IntVector>, SAMRAI::MAX_DIM_VAL> s_zeros;
   static std::array<std::optional<IntVector>, SAMRAI::MAX_DIM_VAL> s_ones;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
