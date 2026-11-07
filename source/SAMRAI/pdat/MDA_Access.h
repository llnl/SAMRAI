/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Light-weight array class
 *
 ************************************************************************/
#ifndef included_MDA_Access_h
#define included_MDA_Access_h

#include "SAMRAI/SAMRAI_config.h"

#include <sys/types.h>
#include <cassert>
#include <iostream>

/* SGS TODO
 * This file is not written to SAMRAI coding style guidelines
 * And has horried commenting conventions (e.g. comments to left of code!)
 */

/*!
 * @file
 * @brief Provides classes supporting Fortran-style
 * multidimensional array accessing in C++.
 *
 * The classes are written for performance (or at least
 * to not degrade performance), so they are almost all
 * inlined with no run-time toggle-able error checking.
 * It is possible that this approach leads to long
 * compile times and large binaries if you are using
 * a not-so-smart compiler.  In theory though, because
 * these classes are not doing any extraneous computations,
 * it generates codes that are as optimizable as any
 * other code doing similar functions, including Fortran
 * codes.
 *
 * Five classes are defined in this file:
 * -# MDA_IndexRange: a class to define and
 *   manipulate index range objects.
 * -# MDA_OrderRowMajor: class with functions
 *   to compute order-dependent info from the index
 *   range.  This version is for row-major order.
 * -# MDA_OrderColMajor: the column-major
 *   counterpart of MDA_OrderRowMajor.
 * -# MDA_Access: class to allow access to individually
 *   indexed elements of a multidimensional array.
 * -# MDA_AccessConst: the const counterpart of
 *   MDA_Access, to allow read-only access to data.
 *
 * To give the compiler the maximum amount of information
 * with which to perform optimization, always use locally
 * scoped objects.  These classes are very light-weight,
 * so copying them is cheap.
 */

/*!
 * @brief Defines index ranges for multidimensional arrays
 *
 * Defines the abstract index range and methods for setting
 * and accessing it.
 *
 * Nothing is known about the ordering of the array.
 */
template<unsigned short MDA_DIM>
class MDA_IndexRange
{

public:
   /*!
    * @brief Type for the dimension counter.
    *
    * Type dir_t could be "unsigned short" instead of "short",
    * but using short causes the GNU compiler to issue unneeded
    * warnings about certain comparisons always being false
    * (when instantiating with MDA_DIM of 1).
    */
   typedef short dir_t;
   /*!
    * @brief Type for the index.
    */
   typedef int index_t;

protected:
   //! @brief Dimension number (to avoid compiler warnings, set to zero when unused).
   enum { D0 = 0,
          D1 = (MDA_DIM > 1 ? 1 : 0),
          D2 = (MDA_DIM > 2 ? 2 : 0),
          D3 = (MDA_DIM > 3 ? 3 : 0),
          D4 = (MDA_DIM > 4 ? 4 : 0) };

public:
   /*!
    * @brief Constructor for setting index data using size and
    * starting points.
    *
    * Any pointers that are NULL are not used.
    * The resulting default settings are:
    * - Array sizes are 0
    * - Array starting indices are 0
    * Since all arguments have default values,
    * this serves at the default constructor.
    *
    * There is another constructor which accepts the first and final
    * indices instead of the sizes and first indices.
    * @b NOTE: the place of the starting points is different than
    * it is for the constructor taking final indices instead of sizes.
    */
   MDA_IndexRange(
      /*! Array sizes                   */ const size_t* sz = ((size_t *)0),
      /*! Array starting indices        */ const index_t* st = ((index_t *)0))
   {
      dir_t i;
      if (st) {
         for (i = 0; i < MDA_DIM; ++i) {
            d_start[i] = st[i];
         }
      } else { for (i = 0; i < MDA_DIM; ++i) {
                  d_start[i] = 0;
               }
      }
      if (sz) {
         for (i = 0; i < MDA_DIM; ++i) {
            d_size[i] = sz[i];
         }
      } else { for (i = 0; i < MDA_DIM; ++i) {
                  d_size[i] = 0;
               }
      }
      setDependentData();
   }

   /*!
    * @brief Constructor for setting index data to a range.
    *
    * This version takes two @c index_t* arguments, for the initial
    * and final indices.  It does not support default arguments
    * until after the indices argument.  @b NOTE: the place of the
    * initial indices is different than it is for the constructor
    * taking sizes instead of final indices.
    *
    * If @c si is @c NULL, starting indices are set to 0.
    * If @c sf is @c NULL, sizes are set to zero.
    */
   MDA_IndexRange(
      /*! Array of initial indices      */ const index_t* si,
      /*! Array of final indices        */ const index_t* sf)
   {
      index_t i;
      if (si) {
         for (i = 0; i < MDA_DIM; ++i) d_start[i] = si[i];
      } else { for (i = 0; i < MDA_DIM; ++i) d_start[i] = 0;
      }
      if (sf) {
         for (i = 0; i < MDA_DIM; ++i) d_size[i] = 1 + sf[i] - d_start[i];
      } else { for (i = 0; i < MDA_DIM; ++i) d_size[i] = 0;
      }
      setDependentData();
   }

   /*!
    * @brief Virtual destructor to support inheritance.
    */
   virtual ~MDA_IndexRange()
   {
   }

   //@{ @name Functions to set indices

   /*!
    * Set size and starting indices.
    */
   void setSizeAndStart(
      /*! Array sizes (NULL for no change)      */ const size_t* sz = ((size_t *)0),
      /*! Starting indices (NULL for no change) */ const index_t* st =
         ((index_t *)0))
   {
      if (sz) for (dir_t i = 0; i < MDA_DIM; ++i) d_size[i] = sz[i];
      if (st) for (dir_t i = 0; i < MDA_DIM; ++i) d_start[i] = st[i];
      setDependentData();
   }

   /*!
    * Set first and final indices (inclusive).
    */
   void setInclusiveRange(
      /*! First valid indices (NULL for no change) */ const index_t first[
         MDA_DIM],
      /*! Final valid indices (NULL for no change) */ const index_t final[
         MDA_DIM])
   {
      if (first) for (dir_t i = 0; i < MDA_DIM; ++i) d_start[i] = first[i];
      if (final) for (dir_t i = 0; i < MDA_DIM; ++i) d_size[i] = final[i]
               - first[i] + 1;
      setDependentData();
   }

   /*!
    * @brief Adjust the array sizes
    *
    * Adjust the first and final indices.
    * Set the direction to adjust to >= MDA_DIM to adjust @em all dimensions.
    * @b Note: The third argument is the increment to the final index
    * and not the size.
    *
    * @b Note: No error checking is done, for example, to make sure that the
    * resulting size is non-negative.
    *
    * @return Adjusted MDA_IndexRange object
    */
   const MDA_IndexRange& adjustDim(
      /*! Direction to adjust      */ dir_t d,
      /*! Increment to first index */ index_t first,
      /*! Increment to final index */ index_t final)
   {
      if (d >= MDA_DIM) {
         dir_t i;
         for (i = 0; i < MDA_DIM; ++i) {
            d_start[i] += first;
            d_size[i] += final - first;
         }
      } else {
         d_start[d] += first;
         d_size[d] += final - first;
      }
      setDependentData();
      return *this;
   }

   //@}

   //@{ @name Comparison functions

   /*!
    * @brief Equivalence comparison
    */
   bool operator == (
      const MDA_IndexRange& r) const
   {
      dir_t d;
      for (d = 0; d < MDA_DIM; ++d) {
         if (d_start[d] != r.d_start[d]) return false;

         if (d_size[d] != r.d_size[d]) return false;
      }
      return true;
   }

   /*!
    * @brief Inequivalence comparison
    */
   bool operator != (
      const MDA_IndexRange& r) const
   {
      return !((*this) == r);
   }

   //@}

   //@{ @name IO functions
   /*!
    * @brief Output to ostream.
    */
   std::ostream& streamPut(
      std::ostream& os) const
   {
      os << MDA_DIM;
      for (dir_t i = 0; i < MDA_DIM; ++i)
         os << ' ' << d_start[i] << ' ' << d_size[i];
      return os;
   }
   /*!
    * @brief Input from istream.
    */
   std::istream& streamGet(
      std::istream& is)
   {
      dir_t dim;
      is >> dim;
      assert(dim == MDA_DIM);
      for (dir_t i = 0; i < MDA_DIM; ++i)
         is >> d_start[i] >> d_size[i];
      return is;
   }

   friend std::ostream& operator << (
      std::ostream& os,
      const MDA_IndexRange<MDA_DIM>& r) {
      return r.streamPut(os);
   }
   friend std::istream& operator >> (
      std::istream& is,
      MDA_IndexRange<MDA_DIM>& r) {
      return r.streamGet(is);
   }
   //@}

   //@{ @name Functions to facilitate looping

   /*!
    * @brief Give starting index of a given direction.
    */
   const index_t& beg(
      /*! index of direction */ size_t i) const
   {
      return d_start[i];
   }

   /*!
    * @brief Give ending index (one more than the last valid index)
    * of a given direction.
    */
   const index_t& end(
      /*! index of direction */ size_t i) const
   {
      return d_stop[i];
   }

   /*!
    * @brief Give size along a given direction.
    */
   const size_t& size(
      /*! index of direction */ size_t i) const
   {
      return d_size[i];
   }

   /*!
    * @brief Give size for all directions.
    */
   size_t totalSize() const
   {
      dir_t i = 0;
      size_t total_size = d_size[i];
      for (i = 1; i < MDA_DIM; ++i) total_size *= d_size[i];
      return total_size;
   }

   //@}

   //@{ @name Error checking functions

   //! Check if indices are in range.
   bool has(
      index_t i0) const
   {
      return (MDA_DIM == 1)
             && (i0 >= d_start[D0]) && (i0 < d_stop[D0]);
   }

   //! Check if indices are in range.
   bool has(
      index_t i0,
      index_t i1) const
   {
      return (MDA_DIM == 2)
             && (i0 >= d_start[D0]) && (i0 < d_stop[D0])
             && (i1 >= d_start[D1]) && (i1 < d_stop[D1]);
   }

   //! Check if indices are in range.
   bool has(
      index_t i0,
      index_t i1,
      index_t i2) const
   {
      return (MDA_DIM == 3)
             && (i0 >= d_start[D0]) && (i0 < d_stop[D0])
             && (i1 >= d_start[D1]) && (i1 < d_stop[D1])
             && (i2 >= d_start[D2]) && (i2 < d_stop[D2]);
   }

   //! Check if indices are in range.
   bool has(
      index_t i0,
      index_t i1,
      index_t i2,
      index_t i3) const
   {
      return (MDA_DIM == 4)
             && (i0 >= d_start[D0]) && (i0 < d_stop[D0])
             && (i1 >= d_start[D1]) && (i1 < d_stop[D1])
             && (i2 >= d_start[D2]) && (i2 < d_stop[D2])
             && (i3 >= d_start[D3]) && (i3 < d_stop[D3]);
   }

   //@}

private:
//! Set dependent data.
   void setDependentData() {
      index_t i;
      for (i = 0; i < MDA_DIM; ++i) {
         d_stop[i] = d_start[i] + static_cast<int>(d_size[i]);
      }
   }

protected:
   //! @brief Array of starting indices
   index_t d_start[MDA_DIM > 0 ? MDA_DIM : 1];

   //! @brief Array of stopping indices
   index_t d_stop[MDA_DIM > 0 ? MDA_DIM : 1];

   //! @brief Array of sizes
   size_t d_size[MDA_DIM > 0 ? MDA_DIM : 1];

};

/**********************************************************************/
/**********************************************************************/

/*!
 * @brief Performs computations based for row-major arrays.
 *
 * This class computes things that are dependent on
 * element order in memory, in this case, for the
 * row-major order.
 */
template<unsigned short MDA_DIM>
class MDA_OrderRowMajor:private MDA_IndexRange<MDA_DIM>
{
public:
   typedef int index_t;
   typedef MDA_IndexRange<MDA_DIM> range_t;
   typedef typename range_t::dir_t dir_t;
protected:
   enum { D0 = range_t::D0,
          D1 = range_t::D1,
          D2 = range_t::D2,
          D3 = range_t::D3 };
public:
   /*!
    * @brief Quantity of (MDA_DIM-1), used only if MDA_DIM > 1.
    * Otherwise defined to 1 to avoid out-of-range subscripts.
    */
   enum { MDA_Reduced_DIM = (MDA_DIM > 1 ? MDA_DIM - 1 : 1) };
   typedef MDA_OrderRowMajor<MDA_Reduced_DIM> reduced_order_t;
   //! @brief Similar to MDA_IndexRange constructor.
   MDA_OrderRowMajor(
      /*! Array sizes            */ const size_t* sz = ((size_t *)0),
      /*! Array starting indices */ const index_t* st = ((index_t *)0)):
      MDA_IndexRange<MDA_DIM>(sz, st)
   {
      computeSizeDependentData();
   }
   //! @brief Similar to MDA_IndexRange constructor.
   MDA_OrderRowMajor(
      /*! Array of initial indices */ const index_t* si,
      /*! Array of final indices   */ const index_t* sf):
      MDA_IndexRange<MDA_DIM>(si, sf)
   {
      computeSizeDependentData();
   }

   /*!
    * @brief Constructor for specifying index range object
    */
   MDA_OrderRowMajor(
      /*! Array index object */ const range_t& r):
      MDA_IndexRange<MDA_DIM>(r) {
      computeSizeDependentData();
   }

   //@{
   //! @name Access to index range
   /*!
    * @brief Const access to the index range object.
    *
    * The index range cannot be modified through this reference.
    * To modify the index range, use other member functions.
    */
   const range_t& range() const
   {
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::setSizeAndStart().
   const MDA_OrderRowMajor& setSizeAndStart(
      const size_t* sz = ((size_t *)0),
      const index_t* st = ((index_t *)0))
   {
      range_t::setSizeAndStart(sz, st);
      computeSizeDependentData();
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::setInclusiveRange().
   const MDA_OrderRowMajor& setInclusiveRange(
      const index_t first[MDA_DIM],
      const index_t final[MDA_DIM])
   {
      range_t::setInclusiveRange(first, final);
      computeSizeDependentData();
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::adjustDim().
   const MDA_OrderRowMajor& adjustDim(
      dir_t d,
      index_t first,
      index_t final)
   {
      range_t::adjustDim(d, first, final);
      computeSizeDependentData();
      return *this;
   }
//@}

//@{
//! @name Logical comparisons
/*!
 * @brief Equivalence comparison.
 *
 * Only independent data is compared, not dependent (redundant) data.
 */
   bool operator == (
      const MDA_OrderRowMajor& r) const
   {
      return range() == r.range();
   }
   /*!
    * @brief Inequivalence comparison.
    *
    * Only independent data is compared, not dependent (redundant) data.
    */
   bool operator != (
      const MDA_OrderRowMajor& r) const
   {
      return range() != r.range();
   }
//@}

//@{
//! @name Functions to compute offsets
/*!
 * @brief Compute offsets for arbitrary @c MDA_DIM
 *
 * This is flexible but not efficient!
 * You should use dimension-specific offset computations whenever possible.
 */
   index_t offset(
      const index_t i[MDA_DIM]) const
   {
      int d;
      size_t o = i[MDA_DIM > 0 ? MDA_DIM - 1 : 0] - this->beg(MDA_DIM - 1);
      for (d = MDA_DIM - 2; d >= 0; --d) o +=
            (i[d] - this->d_start[d]) * d_total_size[d + 1];
      return o;
   }
   index_t offset(
      index_t i0) const
   {
      return i0 - this->beg(0);
   }
   index_t offset(
      index_t i0,
      index_t i1) const
   {
      return (i0 - this->d_start[D0]) * d_total_size[D1]
             + (i1 - this->d_start[D1]);
   }
   index_t offset(
      index_t i0,
      index_t i1,
      index_t i2) const
   {
      return (i0 - this->d_start[D0]) * d_total_size[D1]
             + (i1 - this->d_start[D1]) * d_total_size[D2]
             + (i2 - this->d_start[D2]);
   }
   index_t offset(
      index_t i0,
      index_t i1,
      index_t i2,
      index_t i3) const
   {
      return (i0 - this->d_start[D0]) * d_total_size[D1]
             + (i1 - this->d_start[D1]) * d_total_size[D2]
             + (i2 - this->d_start[D2]) * d_total_size[D3]
             + (i3 - this->d_start[D3]);
   }
   index_t fixedOffset() const
   {
      return d_fixed_offset;
   }
   //@}
/*!
 * @brief Return the total size of subarray starting with direction d
 */
   size_t totalSize(
      unsigned short d) const
   {
      return d_total_size[d];
   }
/*!
 * @brief Computes the order object and offset for reducing the slowest
 * direction.
 *
 * A reduced array is the subarray resulting from fixing the slowest
 * (first) index.  The reduced array has one fewer dimension, a different
 * ordering object and its data starts at a different point in memory.
 * The change in starting point is the returned offset value, and the
 * new order object is returned in the referenced argument.
 *
 * @return Pointer offset (always positive) to the reduced array pointer.
 */
   size_t reduce(
      index_t i,
      reduced_order_t& new_order) const
   {
      new_order.setSizeAndStart(&this->size(1), &this->beg(1));
      return offset(i) * d_total_size[1];
   }
private:
/*!
 * @brief Recompute the total sizes array, which is dependent on sizes.
 */
   void computeSizeDependentData()
   {
      d_total_size[MDA_DIM - 1] = this->d_size[MDA_DIM - 1];
      d_fixed_offset = -this->d_start[MDA_DIM - 1];
      int i = MDA_DIM - 2;
      for ( ; i >= 0; --i) {
         d_total_size[i] = this->d_size[i] * d_total_size[i + 1];
         d_fixed_offset -= this->d_start[i] * d_total_size[i + 1];
      }
   }
/*!
 * @brief Total sizes of sub-dimensional arrays.
 *
 * @c d_total_size[i] is the size of the sub-matrix contained in the
 * last (fast) i directions of the array.  Incidentally, the stride
 * size of direction @c i is @c d_total_size[i+1] (and the stride size
 * for direction @c MDA_DIM-1 is 1.
 *
 * @c d_total_size[i] is really equal to
 * @f$ \prod_{j=i}^{MDA_DIM-1} size_j @f$
 *
 * This member simply caches size-dependent data.
 */
   size_t d_total_size[MDA_DIM > 0 ? MDA_DIM : 1];
/*!
 * @brief The fixed portions of offset calculations.
 *
 * Offsets can be separated into a fixed part (dependent only on range)
 * and a variable part (dependent on dereferencing indices).  To prevent
 * repeated computation of the fixed part, it is saved in this variable.
 * Note that a good optimizing compiler should already do this,
 * so doing it in the code may not really be needed.
 *
 * This member simply caches size-dependent data.
 */
   index_t d_fixed_offset;
};      // end MDA_OrderRowMajor

/**********************************************************************/
/**********************************************************************/

/*!
 * @brief Performs computations based for column-major arrays.
 *
 * This class computes things that are dependent on
 * element order in memory, in this case, for the
 * column-major order.
 */
template<unsigned short MDA_DIM>
class MDA_OrderColMajor:private MDA_IndexRange<MDA_DIM>
{
public:
   typedef int index_t;
   typedef MDA_IndexRange<MDA_DIM> range_t;
   typedef typename range_t::dir_t dir_t;
protected:
   enum { D0 = range_t::D0,
          D1 = range_t::D1,
          D2 = range_t::D2,
          D3 = range_t::D3 };
public:
   /*!
    * @brief Quantity of (MDA_DIM-1), used only if MDA_DIM > 1.
    * Otherwise defined to 1 to avoid out-of-range subscripts.
    */
   enum { MDA_Reduced_DIM = (MDA_DIM > 1 ? MDA_DIM - 1 : 1) };
   typedef MDA_OrderColMajor<MDA_Reduced_DIM> reduced_order_t;
   //! @brief Similar to MDA_IndexRange constructor.
   MDA_OrderColMajor(
      /*! Array sizes            */ const size_t* sz = ((size_t *)0),
      /*! Array starting indices */ const index_t* st = ((index_t *)0)):
      MDA_IndexRange<MDA_DIM>(sz, st)
   {
      computeSizeDependentData();
   }
   //! @brief Similar to MDA_IndexRange constructor.
   MDA_OrderColMajor(
      /*! Array of initial indices */ const index_t* si,
      /*! Array of final indices   */ const index_t* sf):
      MDA_IndexRange<MDA_DIM>(si, sf)
   {
      computeSizeDependentData();
   }

   /*!
    * @brief Constructor for specifying index range object
    */
   MDA_OrderColMajor(
      /*! Array index object */ const range_t& r):
      MDA_IndexRange<MDA_DIM>(r)
   {
      computeSizeDependentData();
   }

   //@{
   //! @name Access to index range (see MDA_IndexRange)
   /*!
    * @brief Const access to the index range object.
    *
    * The index range cannot be modified through this reference.
    * To modify the index range, use other member functions.
    */
   const range_t& range() const
   {
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::setSizeAndStart().
   const MDA_OrderColMajor& setSizeAndStart(
      const size_t* sz = ((size_t *)0),
      const index_t* st = ((index_t *)0))
   {
      range_t::setSizeAndStart(sz, st);
      computeSizeDependentData();
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::setInclusiveRange().
   const MDA_OrderColMajor& setInclusiveRange(
      const index_t first[MDA_DIM],
      const index_t final[MDA_DIM])
   {
      range_t::setInclusiveRange(first, final);
      computeSizeDependentData();
      return *this;
   }
   //! @brief Similar to MDA_IndexRange::adjustDim().
   const MDA_OrderColMajor& adjustDim(
      dir_t d,
      index_t first,
      index_t final)
   {
      range_t::adjustDim(d, first, final);
      computeSizeDependentData();
      return *this;
   }
//@}

//@{
//! @name Logical comparisons
/*!
 * @brief Equivalence comparison.
 *
 * Only independent data is compared, not dependent (redundant) data.
 */
   bool operator == (
      const MDA_OrderColMajor& r) const
   {
      return range() == r.range();
   }
   /*!
    * @brief Inequivalence comparison.
    *
    * Only independent data is compared, not dependent (redundant) data.
    */
   bool operator != (
      const MDA_OrderColMajor& r) const
   {
      return range() != r.range();
   }
//@}

//@{
//! @name Functions to compute offsets
/*!
 * @brief Compute offsets for arbitrary @c MDA_DIM
 *
 * This is flexible but not efficient!
 * You should use dimension-specific offset computations whenever possible.
 */
   index_t offset(
      const index_t i[MDA_DIM]) const
   {
      int d;
      size_t o = i[0] - this->beg(0);
      for (d = 1; d < MDA_DIM; ++d) o +=
            (i[d] - this->d_start[d]) * d_total_size[d - 1];
      return o;
   }
   index_t offset(
      index_t i0) const
   {
      return i0 - this->beg(0);
   }
   index_t offset(
      index_t i0,
      index_t i1) const
   {
      return (i0 - this->d_start[D0])
             + (i1 - this->d_start[D1])
             * static_cast<index_t>(d_total_size[D0]);
   }
   index_t offset(
      index_t i0,
      index_t i1,
      index_t i2) const
   {
      return (i0 - this->d_start[D0])
             + (i1 - this->d_start[D1])
             * static_cast<index_t>(d_total_size[D0])
             + (i2 - this->d_start[D2])
             * static_cast<index_t>(d_total_size[D1]);
   }
   index_t offset(
      index_t i0,
      index_t i1,
      index_t i2,
      index_t i3) const
   {
      return (i0 - this->d_start[D0])
             + (i1 - this->d_start[D1]) * d_total_size[D0]
             + (i2 - this->d_start[D2]) * d_total_size[D1]
             + (i3 - this->d_start[D3]) * d_total_size[D2];
   }
   index_t fixedOffset() const
   {
      return d_fixed_offset;
   }
//@}
/*!
 * @brief Return the total size of subarray starting with direction d
 */
   size_t totalSize(
      unsigned short d) const
   {
      return d_total_size[d];
   }
/*!
 * @brief Computes the order object and offset for reducing the slowest
 * direction.
 *
 * A reduced array is the subarray resulting from fixing the slowest
 * (last) index.  The reduced array has one fewer dimension, a different
 * ordering object and its data starts at a different point in memory.
 * The change in starting point is the returned offset value, and the
 * new order object is returned in the referenced argument.
 *
 * @return Pointer offset (always positive) to the reduced array pointer.
 */
   size_t reduce(
      index_t i,
      reduced_order_t& new_order) const
   {
      new_order.setSizeAndStart(&this->size(0), &this->beg(0));
      return (i
              - this->d_start[MDA_Reduced_DIM])
             * d_total_size[MDA_DIM > 1 ? MDA_DIM - 2 : 0];
      // return offset(i)*d_total_size[MDA_DIM-2];
   }
private:
/*!
 * @brief Recompute the total sizes array, which is dependent on sizes.
 */
   void computeSizeDependentData()
   {
      d_total_size[0] = this->d_size[0];
      d_fixed_offset = -this->d_start[0];
      int i = 1;
      for ( ; i < MDA_DIM; ++i) {
         d_total_size[i] = this->d_size[i] * d_total_size[i - 1];
         d_fixed_offset -= this->d_start[i]
            * static_cast<int>(d_total_size[i]);
      }
   }
/*!
 * @brief Total sizes of sub-dimensional arrays.
 *
 * @c d_total_size[i] is the size of the sub-matrix contained in the
 * first (fast) i+1 directions of the array.  Incidentally, the stride
 * size of direction @c i is @c d_total_size[i-1] (and the stride size
 * for direction @c 0 is 1.
 *
 * @c d_total_size[i] is really equal to
 * @f$ \prod_{j=0}^{i} size_j @f$
 *
 * This member simply caches size-dependent data.
 */
   size_t d_total_size[MDA_DIM > 0 ? MDA_DIM : 1];
/*!
 * @brief The fixed portions of offset calculations.
 *
 * Offsets can be separated into a fixed part (dependent only on range)
 * and a variable part (dependent on indices).  To prevent repeated
 * computation of the fixed part, it is saved in this variable.
 * Note that a good optimizing compiler should already do this,
 * so doing it in the code may not really be needed.
 *
 * This member simply caches size-dependent data.
 */
   index_t d_fixed_offset;
};      // end MDA_OrderColMajor

/**********************************************************************/
/**********************************************************************/

/*!
 * @brief Non-const multidimensional array access.
 *
 * This class @em never allocates or deallocates data.
 * It takes pointers to preallocated data
 * and provides an interface to that data.
 * Member functions are used to give that interface.
 *
 * This class provides functions for explicit index checking,
 * but it does @em NO implicit error checking on either the
 * dimensionality of the array or it size.
 * Such may be done through subclassing.
 *
 * The member functions should all be inlined for better
 * performance.
 *
 * This template class is set up to work with either
 * row-major or column-major data, depending on the
 * third template argument, which should be one of
 * -# @c MDA_OrderRowMajor (default if omitted)
 * -# @c MDA_OrderColMajor
 *
 * The reduce() function return a new array of smaller
 * dimensional that require less integer arithmetic
 * to access individual array members.  This should help
 * in optimizing code.  (My preliminary performance tests
 * using gcc and gprof on i686 Linux showed that the
 * MDA_Access functions run at half to slightly
 * faster than the speed of Fortran, depending on use of
 * reduced arrays.  However, note that gcc is not great at
 * optimizing Fortran.)
 */
template<class MDA_TYPE, unsigned short MDA_DIM, class OrderType =
            MDA_OrderRowMajor<MDA_DIM> >
class MDA_Access
{

public:
   /*!
    * @brief Type of data.
    */
   typedef MDA_TYPE value_t;
   typedef MDA_IndexRange<MDA_DIM> range_t;
   typedef typename range_t::dir_t dir_t;
   typedef typename range_t::index_t index_t;
   typedef OrderType order_t;
   typedef typename OrderType::reduced_order_t reduced_order_t;

   /*!
    * @brief Constructor for setting all data, with default values.
    *
    * Any pointers that are NULL are not used.  The resulting default
    * settings are:
    * - Data pointer is NULL
    * - Array sizes are 0
    * - Array starting indices are 0
    *
    * There is another constructor which accepts the first and final
    * indices instead of the sizes and first indices.
    * @b NOTE: the place of the initial indices is different than
    * it is for the constructor taking final indices instead of sizes.
    */
   MDA_Access(
      /*! Pointer to data        */ value_t* p = ((value_t *)0),
      /*! Array sizes            */ const size_t* sz = ((size_t *)0),
      /*! Array starting indices */ const index_t* st = ((index_t *)0)):
      d_ptr(p),
      d_order(sz, st)
   {
      setPtr1();
   }

   /*!
    * @brief Constructor for setting all data, with default values.
    *
    * Any pointers that are NULL are not used.
    * The resulting default settings are:
    * - Data pointer is NULL
    * - Array sizes are 0
    * - Array starting indices are 0
    *
    * This version takes two @c int* arguments, for the initial
    * and final indices.  It does not support default arguments
    * until after the indices argument.  @b NOTE: the place of the
    * initial indices is different than it is for the constructor
    * taking sizes instead of final indices.
    *
    * If @c si is @c NULL, starting indices are set to 0.
    * If @c sf is @c NULL, sizes are set to zero.
    */
   MDA_Access(
      /*! Pointer to data          */ value_t* p,
      /*! Array of initial indices */ const index_t* si,
      /*! Array of final indices   */ const index_t* sf):
      d_ptr(p),
      d_order(si, sf)
   {
      setPtr1();
   }

   /*!
    * @brief Constructor for specifying pointer and ordering object
    */
   MDA_Access(
      /*! Pointer to data    */ value_t* p,
      /*! Array index object */ const order_t& r):
      d_ptr(p),
      d_order(r)
   {
      setPtr1();
   }

   /*!
    * @brief Copy constructor
    */
   MDA_Access(
      /*! Copyee object */ const MDA_Access& r):
      d_ptr(r.d_ptr),
      d_order(r.d_order)
   {
      setPtr1();
   }

   /*!
    * @brief Virtual destructor to support inheritance.
    */
   virtual ~MDA_Access()
   {
   }

   /*!
    * @brief Assignment operator
    */
   MDA_Access&
   operator = (
      const MDA_Access& r)
   {
      d_ptr = r.d_ptr;
      d_ptr1 = r.d_ptr1;
      d_order = r.d_order;
      return *this;
   }

   /*!
    * @brief Conversion into boolean.
    *
    * @return true iff data pointer is not NULL.
    */
   operator bool () const
   {
      return d_ptr != (value_t *)0;
   }

   /*!
    * @brief Conversion into pointer.
    *
    * @return the data pointer.
    */
   operator value_t * () const
   {
      return d_ptr;
   }

   /*!
    * @brief Set the data pointer.
    */
   void setPointer(
      /*! Pointer value */ value_t* p)
   {
      d_ptr = p;
      setPtr1();
   }

   /*!
    * Set size and starting indices.
    *
    * @see MDA_IndexRange
    */
   void setSizeAndStart(
      /*! Array sizes (NULL for no change)      */ const size_t* sz = ((size_t *)0),
      /*! Starting indices (NULL for no change) */ const index_t* st =
         ((index_t *)0))
   {
      d_order.setSizeAndStart(sz, st);
      setPtr1();
   }

   /*!
    * Set first and final indices (inclusive).
    *
    * @see MDA_IndexRange
    */
   void setInclusiveRange(
      /*! First valid indices (NULL for no change) */ const index_t first[
         MDA_DIM],
      /*! Final valid indices (NULL for no change) */ const index_t final[
         MDA_DIM])
   {
      d_order.setInclusiveRange(first, final);
      setPtr1();
   }

   /*!
    * @brief Adjust the directions
    *
    * @see MDA_IndexRange::adjustDim.
    */
   const range_t& adjustDim(
      /*! Direction to adjust      */ dir_t d,
      /*! Increment to first index */ index_t first,
      /*! Increment to final index */ index_t final)
   {
      d_order.adjustDim(d, first, final);
      setPtr1();
      return d_order.range();
   }

//@{ @name Comparison functions

   /*!
    * @name Equivalence comparison
    */
   bool operator == (
      const MDA_Access& r) const
   {
      if (d_order != r.d_order) {
         return false;
      }
      if (d_ptr != r.d_ptr) {
         return false;
      }
      return true;
   }

   /*!
    * @name Inequivalence comparison
    */
   bool operator != (
      const MDA_Access& r) const
   {
      return !((*this) == r);
   }

//@}

//@{ @name Functions for accessing items

   const range_t& range() const
   {
      return d_order.range();
   }
   const index_t& beg(
      size_t i) const
   {
      return d_order.range().beg(i);
   }
   const index_t& end(
      size_t i) const
   {
      return d_order.range().end(i);
   }
   const size_t& size(
      size_t i) const
   {
      return d_order.range().size(i);
   }

   /*!
    * @brief Grant general access to item in an arbitrary dimensional array.
    *
    * This is flexible but not efficient!
    * You should use dimension-specific accesses whenever possible.
    */
   value_t& operator () (
      const index_t i[MDA_DIM]) const
   {
      return d_ptr[d_order.offset(i)];
   }

   /*!
    * @brief Grant general access to item in a 1D array.
    */
   value_t& operator () (
      index_t i0) const
   {
      return d_ptr[d_order.offset(i0)];
      /*
       * return d_ptr1[i0];
       */
   }

   /*!
    * @brief Grant general access to item in a 2D array.
    */
   value_t& operator () (
      index_t i0,
      index_t i1) const
   {
      return d_ptr[d_order.offset(i0, i1)];
   }

   /*!
    * @brief Grant general access to item in a 3D array.
    */
   value_t& operator () (
      index_t i0,
      index_t i1,
      index_t i2) const
   {
      return d_ptr[d_order.offset(i0, i1, i2)];
   }

   /*!
    * @brief Grant general access to item in a 4D array.
    */
   value_t& operator () (
      index_t i0,
      index_t i1,
      index_t i2,
      index_t i3) const
   {
      return d_ptr[d_order.offset(i0, i1, i2, i3)];
   }

   /*!
    * @brief Special case for 1D arrays, identical to @c operator(index_t),
    * using pre-added fixed offsets.
    *
    * This @em may be more efficient than @c (i) but it only works in 1D.
    * It is not guaranteed to work if the fixed offset is negative and
    * has greater value than the pointer address, since the addition of
    * the two gives a negative address, which the C standard leaves as
    * undefined behavior.
    */
   value_t& operator [] (
      index_t i0) const
   {
      return d_ptr1[i0];
      /*
       * return d_ptr[d_order.offset(i0)];
       */
   }

//@}

//@{ @name Functions to extract reduced-dimensional arrays.

/*!
 * @brief Fix the index of the slowest direction and return
 * the corresponding sub-array.
 *
 * This function is meant to facilitate optimization when using
 * this class.  In nested loops, the inner loops executes many
 * times with the indices corresponding to outer loops remaining
 * constants.  This leads to many many repeated integer arithmetics
 * operations that could be removed from the inner loop (but may
 * not be removed automatically by the compiler optimization).
 * To do this, reduce the array dimensionality one direction at a time,
 * by fixing index corresponding to the slowest varying direction.
 * (If you code is written to maximize cache data, this will NOT
 * be the index of the innermost loop.)
 *
 * To reduce multiple directions, string these calls together,
 * i.e. @c array.reduce(i).reduce(j).  However, since reduction
 * contains loops that cost O(MDA_DIM) and may be difficult for
 * compilers to optimize, you may want to save @c array.reduce(i)
 * and reuse it.
 *
 * @param i Index in slowest direction, which is the first
 * direction in a row-major array and the last direction in a
 * column-major array.
 *
 * @return The sub-array of dimension @c MDA_DIM-1, corresponding to
 * the index given.
 */
   MDA_Access<MDA_TYPE, OrderType::MDA_Reduced_DIM,
              typename OrderType::reduced_order_t> reduce(
      index_t i) const
   {
      typename OrderType::reduced_order_t new_order;
      int ptr_offset;
      ptr_offset = d_order.reduce(i, new_order);
      return MDA_Access<MDA_TYPE, OrderType::MDA_Reduced_DIM,
                        typename OrderType::reduced_order_t>(
                d_ptr + ptr_offset, new_order);
   }

//@}

//! Pointer to data.
private:
   value_t* d_ptr;
/*!
 * @brief Value of @c d_ptr-beg(0), used for optimizing 1D access.
 *
 * The use of precomputed @c d_ptr1=d_ptr-beg(0) speeds up 1D offset
 * computations by allowing us to compute  @c d_ptr1+i0 instead of
 * @c d_ptr+(i0-beg(0)), saving one integer subtraction for each
 * 1D access.  However, this could be a real problem if @c d_ptr<beg(0).
 * So far, that has not happened, and we are keeping our fingers
 * crossed.
 *
 * @see setPtr1()
 */
private:
   value_t* d_ptr1;
private:
   void setPtr1() {
      /*
       * If the following assert fails, our d_ptr1 optimization may
       * give undefined result.
       * assert( d_order.fixedOffset() > 0 ||
       *      (unsigned long)d_ptr > (unsigned long)(-d_order.fixedOffset()) );
       */
      d_ptr1 = d_ptr + d_order.fixedOffset();
   }
//! Offset computing object
private:
   order_t d_order;

};      // class MDA_Access

/**********************************************************************/
/**********************************************************************/

/*!
 * @brief Const data version of the multidimensional array access
 * template class MDA_Access.
 *
 * This class is almost exactly identical to its non-const
 * counterpart, MDA_Access.  It is used when the data
 * is const.
 *
 * This class differs only in that the value type is a const.
 * In fact, this class is trivial,
 * except for the public inheritance of MDA_Access
 * with the const type for the first template argument,
 * a constructor to build an object from a MDA_Access
 * object and an assignment operator to assign from a
 * MDA_Access object.
 * Other than that, see MDA_Access for documentations.
 *
 * The interfaces that are added by this class are trivial,
 * mirroring the interfaces defined in MDA_Access
 * with minor changes.
 *
 * @see MDA_Access
 */
template<class MDA_TYPE, unsigned short MDA_DIM, class OrderType =
            MDA_OrderRowMajor<MDA_DIM> >
class MDA_AccessConst:public MDA_Access<const MDA_TYPE, MDA_DIM, OrderType>
{
public:
   /*!
    * @brief Type of data.
    *
    * This declaration is redundant because it should already be inherited,
    * but the xlC compiler on ASCI Blue does not get it.
    */
   typedef const MDA_TYPE value_t;
   typedef MDA_IndexRange<MDA_DIM> range_t;
   typedef typename range_t::dir_t dir_t;
   typedef typename range_t::index_t index_t;
   typedef OrderType order_t;

   /*!
    * @brief See the MDA_Access version of this function.
    * @see MDA_Access::MDA_Access(value_t*,const size_t*,const index_t*)
    */
   MDA_AccessConst(
      /*! Pointer to data        */ value_t* p = ((MDA_TYPE *)0),
      /*! Array sizes            */ const size_t* sz = ((size_t *)0),
      /*! Array starting indices */ const index_t* st = ((index_t *)0)):
      MDA_Access<const MDA_TYPE, MDA_DIM, OrderType>(p, sz, st)
   {
   }
   /*!
    * @brief See the MDA_Access version of this function.
    * @see MDA_Access::MDA_Access(value_t*,const index_t*,const index_t*)
    */
   MDA_AccessConst(
      /*! Pointer to data          */ value_t* p,
      /*! Array of initial indices */ const index_t* si,
      /*! Array of final indices   */ const index_t* sf):
      MDA_Access<const MDA_TYPE, MDA_DIM, OrderType>(p, si, sf)
   {
   }
   /*!
    * @brief See the MDA_Access version of this function.
    * @see MDA_Access::MDA_Access(value_t*,const MDA_IndexRange<MDA_DIM>&)
    */
   MDA_AccessConst(
      /*! Pointer to data    */ value_t* p,
      /*! Array index object */ const MDA_IndexRange<MDA_DIM>& r):
      MDA_Access<const MDA_TYPE, MDA_DIM, OrderType>(p, r)
   {
   }
   /*!
    * @brief Construct from an object of the non-const version.
    * @see MDA_Access::MDA_Access(const MDA_Access<const MDA_TYPE,MDA_DIM>&)
    */
   MDA_AccessConst(
      const MDA_Access<MDA_TYPE, MDA_DIM, OrderType>& r):
      MDA_Access<const MDA_TYPE, MDA_DIM, OrderType>((const MDA_TYPE *)(r),
                                                     r.range())
   {
   }
   /*!
    * @brief Assign value from an object of the non-const version.
    */
   const MDA_AccessConst&
   operator = (
      const MDA_Access<MDA_TYPE, MDA_DIM, OrderType>& r) {
      (MDA_Access<MDA_TYPE, MDA_DIM, OrderType>&)(*this) = r;
      return *this;
   }
};

#endif  // included_MDA_Access_h
