/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility for cataloging periodic shift directions.
 *
 ************************************************************************/
#ifndef included_hier_PeriodicShiftCatalog
#define included_hier_PeriodicShiftCatalog

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PeriodicId.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Catalogs periodic shifts and assigns each one a unique PeriodicId
 *
 * In a periodic domain, a Box location (or rather its images)
 * can appear to be in many places.  We define the "shift distance" as
 * the IntVector distance in index space from the original grid
 * location to the image.
 *
 * This class catalogs all possible shifts for a periodic
 * configuration and assigns an integer index to each shift.  A shift
 * can be referenced by a single PeriodicId rather than its IntVector
 * distance.
 *
 * TODO: Periodic shift definition and management should probably be
 * be combined with multi-block relationships in the multi-block
 * re-design.
 */

class PeriodicShiftCatalog
{
public:

   /*!
    * @brief Constructs uninitialized object.
    *
    * After construction, the object should be populated with periodic
    * shifts, if any exist, using initializeShiftsByIndexDirections().
    */
   explicit PeriodicShiftCatalog(
      const tbox::Dimension& dim);

   ~PeriodicShiftCatalog() = default;

   /*!
    * @brief Return the shift distance corresponding to the given
    * shift number.
    *
    * @param[in] shift_number
    *
    * @pre shift_number.isValid()
    * @pre shift_number.getPeriodicValue() < getNumberOfShifts()
    */
   const IntVector&
   shiftNumberToShiftDistance(
      const PeriodicId& shift_number) const
   {
      TBOX_ASSERT(shift_number.isValid());
      TBOX_ASSERT(shift_number.getPeriodicValue() < (int)d_shifts.size());
      return d_shifts[shift_number.getPeriodicValue()];
   }

   /*!
    * @brief Return the shift number corresponding to the given shift
    * distance.
    *
    * The shift distance must correspond to an actual distance in the
    * catalog.  Otherwise, invalidShiftNumber() is returned.
    *
    * @param[in] shift_distance
    */
   PeriodicId
   shiftDistanceToShiftNumber(
      const IntVector& shift_distance) const
   {
      unsigned int s;
      for (s = 0; s < d_shifts.size(); ++s) {
         if (d_shifts[s] == shift_distance) {
            break;
         }
      }
      if (s == d_shifts.size()) {
         s = tbox::MathUtilities<int>::getMax();
      }
      return PeriodicId(s);
   }

   /*!
    * @brief Return the "invalid" shift number.
    */
   int
   getInvalidShiftNumber() const
   {
      return tbox::MathUtilities<int>::getMax();
   }

   /*!
    * @brief Return the shift number corresponding to the negative of
    * the shift distance of the given shift number.
    *
    * @param[in] shift_number
    *
    * @pre shift_number.isValid()
    * @pre shift_number.getPeriodicValue() < getNumberOfShifts()
    */
   const PeriodicId&
   getOppositeShiftNumber(
      const PeriodicId& shift_number) const
   {
      TBOX_ASSERT(shift_number.isValid());
      TBOX_ASSERT(shift_number.getPeriodicValue() < (int)d_shifts.size());
      return d_opposite_number[shift_number.getPeriodicValue()];
   }

   /*!
    * @brief Determine whether there is any periodicity.
    */
   bool
   isPeriodic() const
   {
      return d_shifts.size() > 1;
   }

   /*!
    * @brief Return the number of possible shifts.
    */
   int
   getNumberOfShifts() const
   {
      return (int)d_shifts.size();
   }

   /*!
    * @brief Return the shift number corresponding to no shift.
    */
   const PeriodicId&
   getZeroShiftNumber() const
   {
      return PeriodicId::zero();
   }

   /*!
    * @brief Compute all possible shifts for periodicity along the
    * index directions and initialize the shifts in the catalog.
    *
    * The shift_distances should be non-negative integers representing
    * the minimum distance (in reference index space) in each
    * direction that a domain has to be shifted before it repeats.  In
    * non-periodic directions, use the value of zero.
    *
    * Periodic shifts along index directions is the only periodic
    * shifting that SAMRAI fully supports.  There are 3 possible
    * shifts per direction: forward, backward and no-shift.  In
    * D-dimension, there can be up to 3^D images of the same location,
    * including diagonal shifts that results from combining multiple
    * shifts along different directions.
    *
    * This method must only be called once for each dimension
    * (indicated by the dimension of @c shift_distances).
    *
    * @param[in] shift_distance_along_index_directions The periodic
    * shift distance in each index direction.
    */
   void
   initializeShiftsByIndexDirections(
      const IntVector& shift_distance_along_index_directions);

private:

   /*!
    * @brief Set the shifts.
    *
    * Each IntVector in shifts must represent a minimum distance (in
    * reference index space) that the domain can be shifted for it to
    * overlay itself exactly.  For each shift given, its opposite (its
    * negative value) is automatically added.  Duplicate shifts are
    * disregarded.
    *
    * @param[in] dim
    *
    * @param[in] shifts All possible shifts.
    */
   void
   setShifts(
      const tbox::Dimension& dim,
      const std::vector<IntVector>& shifts);

   std::vector<IntVector> d_shifts;

   /*!
    * @brief Shift numbers corresponding to the shift in the opposite
    * direction.
    *
    * Always true: d_shifts[i] == -d_shifts[ d_opposite_number[i] ]
    */
   std::vector<PeriodicId> d_opposite_number;

   /*!
    * @brief The shift number corresponding to zero shift.
    */
   unsigned int d_zero_shift_number;

};

}
}

#endif
