/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index
#define included_hier_Index

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include <array>
#include <optional>

namespace SAMRAI {
namespace hier {

/**
 * Class Index implements a simple n-dimensional integer vector in the
 * AMR index space.  Index is used as lower and upper bounds when
 * creating a box and also when iterating over the cells in a box.  An
 * Index is essentially an integer vector but it carries along the
 * notion of indexing into AMR's abstract index space.
 *
 * @see Box
 * @see BoxIterator
 * @see IntVector
 */

class Index
{
public:

   typedef tbox::Dimension::dir_t dir_t;

   /**
    * @brief Creates an uninitialized Index.
    */
   explicit Index(
      const tbox::Dimension& dim) noexcept;

   /**
    * @brief Construct an Index with all components equal to the argument.
    */
   Index(
      const tbox::Dimension& dim,
      int value) noexcept;

   /**
    * @brief Construct a two-dimensional Index with the value (i,j).
    */
   Index(
      int i,
      int j) noexcept;

   /**
    * @brief Construct a three-dimensional Index with the value (i,j,k).
    */
   Index(
      int i,
      int j,
      int k) noexcept;

   /**
    * @brief Construct an n-dimensional Index with the values copied
    *        from the integer tbox::Array i of size n.
    *
    * The dimension of the constructed Index will be equal to the size of the
    * argument vector.
    *
    * @pre i.size() > 0
    */
   explicit Index(
      const std::vector<int>& i) noexcept;

   /**
    * @brief The copy constructor creates an Index equal to the argument.
    */
   Index(
      const Index& rhs) noexcept;

   /**
    * @brief Construct an Index equal to the argument IntVector.
    *
    * @pre rhs.getNumBlocks() == 1
    */
   explicit Index(
      const IntVector& rhs) noexcept;

   /**
    * @brief Construct an Index equal to the argument array.
    */
   Index(
      const tbox::Dimension& dim,
      const int array[]) noexcept;

   /**
    * @brief The assignment operator sets the Index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   Index&
   operator = (
      const Index& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] = rhs.d_index[i];
      }
      return *this;
   }

   /**
    * @brief The assignment operator sets the Index equal to the argument
    *        IntVector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   operator = (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] = rhs[i];
      }
      return *this;
   }

   /**
    * @brief The Index destructor does nothing interesting.
    */
   virtual ~Index() noexcept = default;

   /**
    * @brief Returns true if all components are equal.
    */
   constexpr bool
   operator == (
      const Index& rhs) const noexcept
   {
      bool result = true;
      for (unsigned int i = 0; result && (i < getDim().getValue()); ++i) {
         result = d_index[i] == rhs.d_index[i];
      }
      return result;
   }

   /**
    * @brief Returns true if one or more components are not equal.
    */
   constexpr bool
   operator != (
      const Index& rhs) const noexcept
   {
      return !(*this == rhs);
   }

   /**
    * @brief Plus-equals operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   operator += (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] += rhs[i];
      }
      return *this;
   }

   /**
    * @brief Plus operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index
   operator + (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * @brief Plus-equals operator for an Index
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr Index&
   operator += (
      const Index& rhs) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] += rhs.d_index[i];
      }
      return *this;
   }

   /**
    * @brief Plus operator for an Index and an another Index
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator + (
      const Index& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * @brief Plus-equals operator for an Index and an integer.
    */
   constexpr Index&
   operator += (
      int rhs) noexcept
   {
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] += rhs;
      }
      return *this;
   }

   /**
    * @brief Plus operator for an Index and an integer.
    */
   Index
   operator + (
      int rhs) const noexcept
   {
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * @brief Minus-equals operator for an Index
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr Index&
   operator -= (
      const Index& rhs) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] -= rhs.d_index[i];
      }
      return *this;
   }

   /**
    * @brief Minus operator for an Index
    *
    * @pre getDim() == rhs.getDim()
    */
   Index
   operator - (
      const Index& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Minus-equals operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   operator -= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] -= rhs[i];
      }
      return *this;
   }

   /**
    * @brief Minus operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index
   operator - (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Minus-equals operator for an Index and an integer.
    */
   constexpr Index&
   operator -= (
      int rhs) noexcept
   {
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         d_index[i] -= rhs;
      }
      return *this;
   }

   /**
    * @brief Minus operator for an Index and an integer.
    */
   Index
   operator - (
      int rhs) const noexcept
   {
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * @brief Times-equals operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   operator *= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         d_index[i] *= rhs[i];
      }
      return *this;
   }

   /**
    * @brief Times operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index
   operator * (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * @brief Times-equals operator for an Index and an integer.
    */
   constexpr Index&
   operator *= (
      int rhs) noexcept
   {
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         d_index[i] *= rhs;
      }
      return *this;
   }

   /**
    * @brief Times operator for an Index and an integer.
    */
   Index
   operator * (
      int rhs) const noexcept
   {
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * @brief Assign-quotient operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   operator /= (
      const IntVector& rhs) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         d_index[i] /= rhs[i];
      }
      return *this;
   }

   /**
    * @brief Quotient operator for an Index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index
   operator / (
      const IntVector& rhs) const noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      TBOX_ASSERT(rhs.getNumBlocks() == 1);
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * @brief Assign-quotient operator for an Index and an integer.
    */
   constexpr Index&
   operator /= (
      int rhs) noexcept
   {
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         d_index[i] /= rhs;
      }
      return *this;
   }

   /**
    * @brief Quotient operator for an Index and an integer.
    */
   Index
   operator / (
      int rhs) const noexcept
   {
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * @brief Return the specified component of the Index.
    *
    * @pre (i < getDim().getValue())
    */
   constexpr int&
   operator [] (
      const unsigned int i) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(i < getDim().getValue());
      return d_index[i];
   }

   /**
    * @brief Return the specified component of the vector as a const reference.
    *
    * @pre (i < getDim().getValue())
    */
   constexpr const int&
   operator [] (
      const unsigned int i) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(i < getDim().getValue());
      return d_index[i];
   }

   /**
    * @brief Return the specified component of the Index.
    *
    * @pre (i < getDim().getValue())
    */
   constexpr int&
   operator () (
      const unsigned int i) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(i < getDim().getValue());
      return d_index[i];
   }

   /**
    * @brief Return the specified component of the Index as a const reference.
    *
    * @pre (i < getDim().getValue())
    */
   constexpr const int&
   operator () (
      const unsigned int i) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(i < getDim().getValue());
      return d_index[i];
   }

   /**
    * @brief Returns true if each integer in Index is greater than
    *        corresponding integer in comparison Index.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr bool
   operator > (
      const Index& rhs) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      bool result = true;
      for (unsigned int i = 0; result && (i < getDim().getValue()); ++i) {
         result = result && (d_index[i] > rhs.d_index[i]);
      }
      return result;
   }

   /**
    * @brief Returns true if each integer in Index is greater or equal to
    *        corresponding integer in comparison Index.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr bool
   operator >= (
      const Index& rhs) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      bool result = true;
      for (unsigned int i = 0; result && (i < getDim().getValue()); ++i) {
         result = result && (d_index[i] >= rhs.d_index[i]);
      }
      return result;
   }

   /**
    * @brief Returns true if each integer in Index is less than
    *        corresponding integer in comparison Index.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr bool
   operator < (
      const Index& rhs) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      bool result = true;
      for (unsigned int i = 0; result && (i < getDim().getValue()); ++i) {
         result = result && (d_index[i] < rhs.d_index[i]);
      }
      return result;
   }

   /**
    * @brief Returns true if each integer in Index is less than or equal to
    *        corresponding integer in comparison Index.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr bool
   operator <= (
      const Index& rhs) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      bool result = true;
      for (unsigned int i = 0; result && (i < getDim().getValue()); ++i) {
         result = result && (d_index[i] <= rhs.d_index[i]);
      }
      return result;
   }

   /**
    * @brief Set Index the component-wise minimum of two Index objects.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr void
   min(
      const Index& rhs) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      for (dir_t i = 0; i < getDim().getValue(); ++i) {
         if (rhs.d_index[i] < d_index[i]) {
            d_index[i] = rhs.d_index[i];
         }
      }
   }

   /**
    * @brief Set Index the component-wise maximum of two Index objects.
    *
    * @pre getDim() == rhs.getDim()
    */
   constexpr void
   max(
      const Index& rhs) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(getDim() == rhs.getDim());
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         if (rhs.d_index[i] > d_index[i]) {
            d_index[i] = rhs.d_index[i];
         }
      }
   }

   /*!
    * @brief Coarsen the Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    *
    * @pre getDim() == ratio.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   Index&
   coarsen(
      const IntVector& ratio) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio);
      TBOX_ASSERT(ratio.getNumBlocks() == 1);
      for (unsigned int d = 0; d < getDim().getValue(); ++d) {
         (*this)(d) = coarsen((*this)(d), ratio(d));
      }
      return *this;
   }

   /*!
    * @brief Return an Index of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getZeroIndex(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_zeros[dim.getValue() - 1].has_value());
      return (s_zeros[dim.getValue() - 1].value());
   }

   /*!
    * @brief Return an Index of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getOneIndex(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_ones[dim.getValue() - 1].has_value());
      return (s_ones[dim.getValue() - 1].value());
   }

   /*!
    * @brief Return an Index with minimum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMinIndex(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_mins[dim.getValue() - 1].has_value());
      return (s_mins[dim.getValue() - 1].value());
   }

   /*!
    * @brief Return an Index with maximum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMaxIndex(
      const tbox::Dimension& dim)
   {
      TBOX_ASSERT(s_maxs[dim.getValue() - 1].has_value());
      return (s_maxs[dim.getValue() - 1].value());
   }

   /*!
    * @brief Coarsen an Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    *
    * @pre index.getDim() == ratio.getDim()
    * @pre rhs.getNumBlocks() == 1
    */
   static Index
   coarsen(
      const Index& index,
      const IntVector& ratio) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(index, ratio);
      TBOX_ASSERT(ratio.getNumBlocks() == 1);
      tbox::Dimension dim(index.getDim());
      Index tmp(dim);
      for (unsigned int d = 0; d < dim.getValue(); ++d) {
         tmp(d) = coarsen(index(d), ratio(d));
      }
      return tmp;
   }

   /*!
    * @brief Get the Dimension of the Index
    */
   constexpr const tbox::Dimension&
   getDim() const noexcept
   {
      return d_dim;
   }

   /**
    * @brief Read an input stream into an Index.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      Index& rhs);

   /**
    * @brief Write an integer index into an output stream.  The format for
    *        the output is (i0,...,in) for an n-dimensional index.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Index& rhs);


   /**
    * @brief Utility function to take the minimum of two Index objects.
    *
    * @pre a.getDim() == b.getDim()
    */
   static Index
   min(
      const Index& a,
      const Index& b) noexcept
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(a, b);
      Index tmp = a;
      tmp.min(b);
      return tmp;
   }

private:
   /*
    * Unimplemented default constructor
    */
   Index() = delete;

   static constexpr int
   coarsen(
      int index,
      int ratio) noexcept
   {
      return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
   }

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
    *
    */
   static void
   finalizeCallback();

   static std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> s_zeros;
   static std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> s_ones;
   static std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> s_maxs;
   static std::array<std::optional<Index>, SAMRAI::MAX_DIM_VAL> s_mins;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

   tbox::Dimension d_dim;

   std::array<int, SAMRAI::MAX_DIM_VAL> d_index;


};

}
}

#endif
