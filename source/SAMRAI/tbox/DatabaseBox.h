/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A box structure representing a portion of the AMR index space
 *
 ************************************************************************/

#ifndef included_tbox_DatabaseBox
#define included_tbox_DatabaseBox

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#include <array>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief POD data for class DatabaseBox
 *
 * The data in DatabaseBox need to reside in a POD class so that
 * HDF5's HOFFSET macro works.  (According to ANSI C++ standard,
 * it does not have to work with non-POD data.)
 */
struct DatabaseBox_POD {
   int d_dimension;
   std::array<int, SAMRAI::MAX_DIM_VAL> d_lo;
   std::array<int, SAMRAI::MAX_DIM_VAL> d_hi;
};

/**
 * Class DatabaseBox represents a box of up to MAX_DIM_VAL
 * dimensions in the AMR index space.  It is defined by lower and
 * upper bounds given by integer arrays.
 *
 * This box is an auxilliary data structure used by the database routines to
 * manipulate boxes.  This box type removes cyclic dependencies among the
 * database routines (which need a box) and the box (which needs the database
 * routines).  The box classes in the hierarchy package convert this box
 * structure into the standard SAMRAI box class used by the AMR algorithms.
 *
 * @internal This class should have @em NO data except for d_data.
 * See d_data for details.
 */

class DatabaseBox
{
public:
   /**
    * The default constructor creates a zero dimension empty box.
    */
   constexpr DatabaseBox()
   : d_data{}
   {
   }

   /**
    * Create a box of the specified dimension describing the index space
    * between lower and upper.
    */
   DatabaseBox(
      const Dimension& dim,
      const int * lower,
      const int * upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   constexpr DatabaseBox(
      const DatabaseBox& box) noexcept = default;

   /**
    * The assignment operator copies the index space of the argument box.
    */
   constexpr DatabaseBox&
   operator = (
      const DatabaseBox& box) noexcept = default;

   /**
    * The destructor does nothing interesting.
    */
   ~DatabaseBox() = default;

   /**
    * Return whether the box is empty.  A box is empty if it has dimension
    * zero or if any part of the upper index is less than its corresponding
    * part of the lower index.
    */
   constexpr bool empty() const noexcept
   {
      if (d_data.d_dimension == 0) {
         return true;
      }
      for (int i = 0; i < d_data.d_dimension; i++) {
         if (d_data.d_hi[i] < d_data.d_lo[i]) {
            return true;
         }
      }
      return false;
   }

   /**
    * Return the dimension of this object.
    */
   constexpr Dimension::dir_t
   getDimVal() const noexcept
   {
      return static_cast<Dimension::dir_t>(d_data.d_dimension);
   }

   void
   setDim(
      const Dimension& dim) noexcept
   {
      const int dim_val = dim.getValue();
      TBOX_ASSERT(dim_val >= 0 && dim_val <= SAMRAI::MAX_DIM_VAL);

      d_data.d_dimension = dim_val;
   }

   /**
    * Return the specified component (non-const) of the lower index of the box.
    *
    * @pre (i >= 0) && (i < getDimVal())
    */
   int&
   lower(
      const int i) noexcept
   {
      TBOX_ASSERT((i >= 0) && (i < getDimVal()));
      return d_data.d_lo[i];
   }

   /**
    * Return the specified component (non-const) of the upper index of the box.
    *
    * @pre (i >= 0) && (i < getDimVal())
    */
   int&
   upper(
      const int i) noexcept
   {
      TBOX_ASSERT((i >= 0) && (i < getDimVal()));
      return d_data.d_hi[i];
   }

   /**
    * Return the specified component (const) of the lower index of the box.
    *
    * @pre (i >= 0) && (i < getDimVal())
    */
   constexpr int
   lower(
      const int i) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT((i >= 0) && (i < getDimVal()));
      return d_data.d_lo[i];
   }

   /**
    * Return the specified component (const) of the upper index of the box.
    *
    * @pre (i >= 0) && (i < getDimVal())
    */
   constexpr int
   upper(
      const int i) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT((i >= 0) && (i < getDimVal()));
      return d_data.d_hi[i];
   }

   /**
    * Check whether two boxes represent the same portion of index space.
    */
   constexpr bool
   operator == (
      const DatabaseBox& box) const noexcept
   {
      const int dim = d_data.d_dimension;
      if (dim != box.d_data.d_dimension) {
         return false;
      }

      for (int i = 0; i < dim; ++i) {
         if (d_data.d_lo[i] != box.d_data.d_lo[i] ||
             d_data.d_hi[i] != box.d_data.d_hi[i]) {
            return false;
         }
      }
      return true;
   }

   /**
    * @brief All data members in a POD type.
    *
    * Due to the need to compute offsets for data members and that
    * offsets cannot be computed for non-POD data, we place all
    * data members in a POD struct and own an object of that
    * struct.
    *
    * Data members are public so that the HDFDatabase need not
    * mirror this structure in defining a compound type for HDF.
    */
   DatabaseBox_POD d_data;
};

}
}

#endif
