/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/

#ifndef included_hier_BoxId
#define included_hier_BoxId

#include <iostream>

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/GlobalId.h"
#include "SAMRAI/hier/PeriodicId.h"
#include "SAMRAI/tbox/MessageStream.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Identifier for a Box, consisting of a GlobalId and a PeriodicId
 *
 * Boxes are identified by their GlobalId and PeriodicId.
 * A Box and all its periodic images have the same GlobalId but
 * different PeriodicId.
 *
 * Comparison operators are provided to define sorted ordering of
 * objects.  The GlobalId and PeriodicId are used for all
 * comparisons.  The GlobalIds are compared first, followed by
 * the PeriodicIds.
 */
class BoxId
{

public:
   /*!
    * @brief Default constructor uses the default constructors for the
    * GlobalId and PeriodicId.
    *
    * The object can be changed using initialize() or by assignment.
    */
   BoxId() noexcept = default;

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id
    *
    * @pre periodic_id.isValid()
    */
   constexpr BoxId(
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero()) noexcept :
      d_global_id(local_id, owner_rank),
      d_periodic_id(periodic_id)
   {
      TBOX_CONSTEXPR_ASSERT(periodic_id.isValid());
   }

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] id
    *
    * @param[in] periodic_id
    *
    * @pre periodic_id.isValid()
    */
   constexpr explicit BoxId(
      const GlobalId& id,
      const PeriodicId& periodic_id = PeriodicId::zero()) noexcept :
      d_global_id(id),
      d_periodic_id(periodic_id)
   {
      TBOX_CONSTEXPR_ASSERT(periodic_id.isValid());
   }

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    *
    * @pre other.periodic_id.isValid()
    */
   constexpr BoxId(
      const BoxId& other) noexcept = default;

   /*!
    * @brief Assignment operator
    */
   constexpr BoxId&
   operator = (
      const BoxId& r) noexcept = default;

   /*!
    * @brief Destructor.
    */
   ~BoxId() noexcept = default;

   /*!
    * @brief Set all the attributes to given values.
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id
    */
   void
   initialize(
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero()) noexcept
   {
      TBOX_CONSTEXPR_ASSERT(periodic_id.isValid());
      d_global_id.getLocalId() = local_id;
      d_global_id.getOwnerRank() = owner_rank;
      d_periodic_id = periodic_id;
   }

   /*!
    * @brief Access the GlobalId.
    */
   constexpr const GlobalId&
   getGlobalId() const noexcept
   {
      return d_global_id;
   }

   /*!
    * @brief Access the owner rank.
    */
   constexpr int
   getOwnerRank() const noexcept
   {
      return d_global_id.getOwnerRank();
   }

   /*!
    * @brief Access the LocalId.
    */
   constexpr const LocalId&
   getLocalId() const noexcept
   {
      return d_global_id.getLocalId();
   }

   /*!
    * @brief Access the PeriodicId.
    */
   constexpr const PeriodicId&
   getPeriodicId() const noexcept
   {
      return d_periodic_id;
   }

   /*!
    * @brief Whether the PeriodicId refers to a periodic
    * image.
    */
   bool
   isPeriodicImage() const noexcept
   {
      return d_periodic_id != PeriodicId::zero();
   }

   /*!
    * @brief Whether the BoxId is valid--meaning it has a valid
    * GlobalId and PeriodicId.
    */
   constexpr bool
   isValid() const noexcept
   {
      return d_periodic_id.isValid() &&
             d_global_id.getLocalId() != LocalId::getInvalidId() &&
             d_global_id.getLocalId() >= 0 &&
             d_global_id.getOwnerRank() != tbox::SAMRAI_MPI::getInvalidRank();
   }

   //@{

   //! @name Comparison operators

   /*!
    * @brief Equality operator.
    *
    * All comparison operators use the GlobalId and PeriodicId.
    */

   constexpr bool
   operator == (
      const BoxId& r) const noexcept
   {
      bool rval = d_global_id == r.d_global_id &&
         d_periodic_id == r.d_periodic_id;
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      return rval;
   }

   /*!
    * @brief Inequality operator.
    *
    * See note on comparison for operator==(const BoxId&);
    */
   constexpr bool
   operator != (
      const BoxId& r) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      bool rval = d_global_id != r.d_global_id ||
         d_periodic_id != r.d_periodic_id;
      return rval;
   }

   /*!
    * @brief Less-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare the
    * LocalIds next; if they compare equal, compare the PeriodicIds.
    */
   constexpr bool
   operator < (
      const BoxId& r) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      return d_global_id.getOwnerRank() < r.d_global_id.getOwnerRank() ||
             (d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
              (d_global_id.getLocalId() < r.d_global_id.getLocalId() ||
               (d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
                (d_periodic_id < r.d_periodic_id))));
   }

   /*!
    * @brief Greater-than operator.
    *
    * Compare the owner ranks first; if they compare equal, compare the
    * LocalIds next; if they compare equal, compare the PeriodicIds.
    */
   constexpr bool
   operator > (
      const BoxId& r) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      return d_global_id.getOwnerRank() > r.d_global_id.getOwnerRank() ||
             (d_global_id.getOwnerRank() == r.d_global_id.getOwnerRank() &&
              (d_global_id.getLocalId() > r.d_global_id.getLocalId() ||
               (d_global_id.getLocalId() == r.d_global_id.getLocalId() &&
                (d_periodic_id > r.d_periodic_id))));
   }

   /*!
    * @brief Less-than-or-equal-to operator.
    */
   constexpr bool
   operator <= (
      const BoxId& r) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      return *this < r || *this == r;
   }

   /*!
    * @brief Greater-than-or-equal-to operator.
    */
   constexpr bool
   operator >= (
      const BoxId& r) const noexcept
   {
      TBOX_CONSTEXPR_ASSERT(d_periodic_id.isValid() && r.d_periodic_id.isValid());
      return *this > r || *this == r;
   }

   //@}

   //@{

   //! @name Support for message passing

   /*!
    * @brief Give number of ints required for putting a BoxId in
    * message passing buffer.
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   constexpr static int
   commBufferSize() noexcept
   {
      return 3;
   }

   /*!
    * @brief Put self into a int buffer.
    *
    * This is the opposite of getFromIntBuffer().  Number of ints
    * written is given by commBufferSize().
    */
   void
   putToIntBuffer(
      int* buffer) const noexcept
   {
      buffer[0] = getOwnerRank();
      buffer[1] = getLocalId().getValue();
      buffer[2] = getPeriodicId().getPeriodicValue();
   }

   /*!
    * @brief Set attributes according to data in int buffer.
    *
    * This is the opposite of putToIntBuffer().  Number of ints read
    * is given by commBufferSize().
    */
   void
   getFromIntBuffer(
      const int* buffer)
   {
      initialize(LocalId(buffer[1]),
         buffer[0],
         PeriodicId(buffer[2]));
   }

   /*!
    * @brief Put self into a MessageStream.
    *
    * This is the opposite of getFromMessageStream().  Number of ints
    * written is given by commBufferSize().
    */
   void
   putToMessageStream(
      tbox::MessageStream& msg) const
   {
      msg << getOwnerRank();
      msg << getLocalId().getValue();
      msg << getPeriodicId().getPeriodicValue();
   }

   /*!
    * @brief Set attributes according to data in MessageStream.
    *
    * This is the opposite of putToMessageStream().  Number of ints read
    * is given by commBufferSize().
    */
   void
   getFromMessageStream(
      tbox::MessageStream& msg)
   {
      int i1, i2, i3;
      msg >> i1;
      msg >> i2;
      msg >> i3;
      initialize(LocalId(i2), i1, PeriodicId(i3));
   }

   //@}

   /*!
    * @brief Format and insert the object into a stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BoxId& r);

private:
   GlobalId d_global_id;

   PeriodicId d_periodic_id;
};

}
}

#endif  // included_hier_BoxId
