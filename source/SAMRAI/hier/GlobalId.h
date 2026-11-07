/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Globally unique identifier that can be locally determined.
 *
 ************************************************************************/
#ifndef included_hier_GlobalId
#define included_hier_GlobalId

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/LocalId.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A locally determined, globally unique identifier.
 *
 * GlobalIds are used as object identifiers in distributed sets of
 * objects, such as Patches and Boxes.  GlobalIds allow objects
 * to be created in parallel while maintaining globally-unique
 * identifiers for those objects.
 *
 * GlobalId includes an owner rank and a LocalId.  A GlobalId is
 * typically generated only by the owner, so global uniqueness is
 * guaranteed by uniqueness of the LocalId portion.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.  The owners and LocalIds are used for all comparisons.
 * The owner ranks are compared first, followed by the LocalIds.
 */
class GlobalId
{

public:
   /*!
    * @brief Default constructor is non-initializing.
    *
    * Sets rank to an invalid value provided by SAMRAI_MPI.
    *
    * The object can be initialized using the assignment operator or
    * the non-const versions of the getOwnerRank() and getLocalId() methods.
    */
   GlobalId() noexcept :
      d_owner_rank(tbox::SAMRAI_MPI::getInvalidRank()),
      d_local_id(LocalId::getInvalidId())
   {
   }

   /*!
    * @brief Initializing constructor.
    */
   constexpr GlobalId(
      const LocalId& local_id,
      const int owner_rank) noexcept :
      d_owner_rank(owner_rank),
      d_local_id(local_id)
   {
   }

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   constexpr GlobalId(
      const GlobalId& other) noexcept = default;

   /*!
    * @brief Destructor.
    */
   ~GlobalId() noexcept = default;

   /*!
    * @brief Access the owner rank.
    */
   int&
   getOwnerRank() noexcept
   {
      return d_owner_rank;
   }

   /*!
    * @brief Access the owner rank.
    */
   constexpr const int&
   getOwnerRank() const noexcept
   {
      return d_owner_rank;
   }

   /*!
    * @brief Access the LocalId.
    */
   LocalId&
   getLocalId() noexcept
   {
      return d_local_id;
   }

   /*!
    * @brief Access the LocalId.
    */
   constexpr const LocalId&
   getLocalId() const noexcept
   {
      return d_local_id;
   }

   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    * All comparison operators use the owner and rank.  Compare the
    * owner values first; if they compare equal, compare the LocalId
    * next.
    */
   constexpr GlobalId&
   operator = (
      const GlobalId& r) noexcept = default;

   constexpr bool
   operator == (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank == r.d_owner_rank) &&
             (d_local_id == r.d_local_id);
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const GlobalId&);
    */
   constexpr bool
   operator != (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank != r.d_owner_rank) ||
             (d_local_id != r.d_local_id);
   }

   /*!
    * @brief Less-than operator.
    *
    * See note on comparison for operator==(const GlobalId&);
    */
   constexpr bool
   operator < (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank < r.d_owner_rank) ||
             ((d_owner_rank == r.d_owner_rank) && (d_local_id < r.d_local_id));
   }

   /*!
    * @brief Greater-than operator.
    *
    * See note on comparison for operator==(const GlobalId&);
    */
   constexpr bool
   operator > (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank > r.d_owner_rank) ||
             ((d_owner_rank == r.d_owner_rank) && (d_local_id > r.d_local_id));
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    *
    * See note on comparison for operator==(const GlobalId&);
    */
   constexpr bool
   operator <= (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank < r.d_owner_rank) ||
             ((d_owner_rank == r.d_owner_rank) &&
              (d_local_id <= r.d_local_id));
   }

   /*!
    * @brief Greater-thanor-equal-to operator.
    *
    * See note on comparison for operator==(const GlobalId&);
    */
   constexpr bool
   operator >= (
      const GlobalId& r) const noexcept
   {
      return (d_owner_rank > r.d_owner_rank) ||
             ((d_owner_rank == r.d_owner_rank) &&
              (d_local_id >= r.d_local_id));
   }

   //@}

   /*!
    * @brief Format and insert rank and LocalId into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const GlobalId& r);

private:
   /*!
    * @brief Rank of owner of this id.
    */
   int d_owner_rank;

   /*!
    * @brief Local id on the owner process.
    */
   LocalId d_local_id;

};

}
}

#endif  // included_hier_GlobalId
