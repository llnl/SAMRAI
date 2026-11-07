/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_CellIndex
#define included_pdat_CellIndex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <type_traits>

namespace SAMRAI {
namespace pdat {

/**
 * Class CellIndex implements a simple n-dimensional integer
 * vector for cell centered variables.  Cell indices contain an integer
 * index location in AMR index space and are identical to the AMR indices.
 *
 * @see hier::Index
 * @see CellData
 * @see CellGeometry
 * @see CellIterator
 */

class CellIndex:public hier::Index
{
public:
   /**
    * The default constructor for a cell index creates an uninitialized index.
    */
   explicit CellIndex(
      const tbox::Dimension& dim);

   /**
    * Construct a cell index from a regular AMR index.
    *
    */
   explicit CellIndex(
      const hier::Index& rhs);

   /**
    * The copy constructor creates a cell index equal to the argument.
    */
   CellIndex(
      const CellIndex& rhs);

   /**
    * The assignment operator sets the cell index equal to the argument.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator = (
      const CellIndex& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator = (rhs);
      return *this;
   }

   /**
    * The cell index destructor does nothing interesting.
    */
   ~CellIndex();

   /**
    * Plus-equals operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator += (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /*!
    * @brief Plus-equals operator for a cell index and hier::Index.
    *
    * The rhs could be hier::Index or CellIndex.  A deleted template overload
    * disallows any other derived Index types to be uses as the rhs.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator += (
      const hier::Index& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator += (rhs);
      return *this;
   }

   /*!
    * Deleted template overload to restrict rhs to Index or CellIndex.
    */
   template <class T,
      typename std::enable_if<
         std::is_base_of<hier::Index, T>::value &&
         !std::is_same<T, hier::Index>::value &&
         !std::is_same<T, CellIndex>::value,
         int>::type = 0>
   CellIndex&
   operator += (
      const T&) = delete;

   /**
    * Plus operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex
   operator + (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      CellIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /*!
    * @brief Plus operator for a cell index and hier::Index.
    *
    * The rhs could be hier::Index or CellIndex.  A deleted template overload
    * disallows any other derived Index types to be uses as the rhs.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex
   operator + (
      const hier::Index& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      CellIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /*!
    * Deleted template overload to restrict rhs to Index or CellIndex.
    */
   template <class T,
      typename std::enable_if<
         std::is_base_of<hier::Index, T>::value &&
         !std::is_same<T, hier::Index>::value &&
         !std::is_same<T, CellIndex>::value,
         int>::type = 0>
   CellIndex
   operator + (
      const T&) const = delete;

   /**
    * Plus-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator += (
      const int rhs)
   {
      hier::Index::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for a cell index and an integer.
    */
   CellIndex
   operator + (
      const int rhs) const
   {
      CellIndex tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator -= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex
   operator - (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      CellIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /*!
    * @brief Minus-equals operator for a cell index and hier::Index.
    *
    * The rhs could be hier::Index or CellIndex.  A deleted template overload
    * disallows any other derived Index types to be uses as the rhs.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator -= (
      const hier::Index& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator -= (rhs);
      return *this;
   }

   /*!
    * Deleted template overload to restrict rhs to Index or CellIndex.
    */
   template <class T,
      typename std::enable_if<
         std::is_base_of<hier::Index, T>::value &&
         !std::is_same<T, hier::Index>::value &&
         !std::is_same<T, CellIndex>::value,
         int>::type = 0>
   CellIndex&
   operator -= (
      const T&) = delete;

   /*!
    * @brief Minus-equals operator for a cell index and hier::Index.
    *
    * The rhs could be hier::Index or CellIndex.  A deleted template overload
    * disallows any other derived Index types to be uses as the rhs.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex
   operator - (
      const hier::Index& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      CellIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /*!
    * Deleted template overload to restrict rhs to Index or CellIndex.
    */
   template <class T,
      typename std::enable_if<
         std::is_base_of<hier::Index, T>::value &&
         !std::is_same<T, hier::Index>::value &&
         !std::is_same<T, CellIndex>::value,
         int>::type = 0>
   CellIndex
   operator - (
      const T&) const = delete;

   /**
    * Minus-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator -= (
      const int rhs)
   {
      hier::Index::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for a cell index and an integer.
    */
   CellIndex
   operator - (
      const int rhs) const
   {
      CellIndex tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex&
   operator *= (
      const hier::IntVector& rhs)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a cell index and an integer vector.
    *
    * @pre getDim() == rhs.getDim()
    */
   CellIndex
   operator * (
      const hier::IntVector& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      CellIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for a cell index and an integer.
    */
   CellIndex&
   operator *= (
      const int rhs)
   {
      hier::Index::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for a cell index and an integer.
    */
   CellIndex
   operator * (
      const int rhs) const
   {
      CellIndex tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Returns true if two cell index objects are equal.  All components
    * must be the same for equality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator == (
      const CellIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator == (rhs);
   }

   /**
    * Returns true if two cell index objects are not equal.  Any of
    * the components may be different for inequality.
    *
    * @pre getDim() == rhs.getDim()
    */
   bool
   operator != (
      const CellIndex& rhs) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);
      return ((hier::Index *)this)->operator != (rhs);
   }
};

}
}

#endif
