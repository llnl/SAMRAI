/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Box representing a portion of the AMR index space
 *
 ************************************************************************/
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

std::array<std::optional<Box>, SAMRAI::MAX_DIM_VAL> Box::s_emptys{};
std::array<std::optional<Box>, SAMRAI::MAX_DIM_VAL> Box::s_universes{};

tbox::StartupShutdownManager::Handler
Box::s_initialize_finalize_handler(
   Box::initializeCallback,
   0,
   0,
   Box::finalizeCallback,
   tbox::StartupShutdownManager::priorityListElements);

#ifdef BOX_TELEMETRY
// These are to optionally track the cumulative number of Boxes constructed,
// the cumulative number of Box assignments, and the high water mark of
// Boxes in existance at any given time.
int Box::s_cumulative_constructed_ct = 0;

int Box::s_cumulative_assigned_ct = 0;

int Box::s_active_ct = 0;

int Box::s_high_water = 0;
#endif

Box::Box(
   const tbox::Dimension& dim):
   d_lo(dim, tbox::MathUtilities<int>::getMax()),
   d_hi(dim, tbox::MathUtilities<int>::getMin()),
   d_block_id(BlockId::invalidId()),
   d_id(GlobalId(), PeriodicId::zero()),
   d_id_locked(false),
   d_empty_flag(EmptyBoxState::BOX_EMPTY)
{
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and reset
   // the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const Index& lower,
   const Index& upper,
   const BlockId& block_id):
   d_lo(lower),
   d_hi(upper),
   d_block_id(block_id),
   d_id(GlobalId(), PeriodicId::zero()),
   d_id_locked(false),
   d_empty_flag(EmptyBoxState::BOX_UNKNOWN)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(lower, upper);
   TBOX_ASSERT(block_id != BlockId::invalidId());
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and reset
   // the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const Box& box):
   d_lo(box.d_lo),
   d_hi(box.d_hi),
   d_block_id(box.d_block_id),
   d_id(box.d_id),
   d_id_locked(false),
   d_empty_flag(box.d_empty_flag)
{
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and reset
   // the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const tbox::DatabaseBox& box):
   d_lo(tbox::Dimension(static_cast<unsigned short>(box.getDimVal())),
        tbox::MathUtilities<int>::getMax()),
   d_hi(tbox::Dimension(static_cast<unsigned short>(box.getDimVal())),
        tbox::MathUtilities<int>::getMin()),
   d_id(GlobalId(), PeriodicId::zero()),
   d_id_locked(false)
{
   set_Box_from_DatabaseBox(box);
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and reset
   // the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const Index& lower,
   const Index& upper,
   const BlockId& block_id,
   const LocalId& local_id,
   const int owner_rank,
   const PeriodicId& periodic_id):
   d_lo(lower),
   d_hi(upper),
   d_block_id(block_id),
   d_id(local_id, owner_rank, periodic_id),
   d_id_locked(false),
   d_empty_flag(EmptyBoxState::BOX_UNKNOWN)
{
   TBOX_ASSERT(periodic_id.isValid());
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and
   // reset the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const Box& box,
   const LocalId& local_id,
   const int owner,
   const PeriodicId& periodic_id):
   d_lo(box.d_lo),
   d_hi(box.d_hi),
   d_block_id(box.d_block_id),
   d_id(local_id, owner, periodic_id),
   d_id_locked(false),
   d_empty_flag(box.d_empty_flag)
{
   TBOX_ASSERT(periodic_id.isValid());
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and
   // reset the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const tbox::Dimension& dim,
   const GlobalId& global_id,
   const PeriodicId& periodic_id):
   d_lo(dim, tbox::MathUtilities<int>::getMax()),
   d_hi(dim, tbox::MathUtilities<int>::getMin()),
   d_block_id(BlockId::invalidId()),
   d_id(global_id, periodic_id),
   d_id_locked(false),
   d_empty_flag(EmptyBoxState::BOX_EMPTY)
{
   TBOX_ASSERT(periodic_id.isValid());
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and
   // reset the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

Box::Box(
   const tbox::Dimension& dim,
   const BoxId& box_id):
   d_lo(dim, tbox::MathUtilities<int>::getMax()),
   d_hi(dim, tbox::MathUtilities<int>::getMin()),
   d_id(box_id),
   d_id_locked(false),
   d_empty_flag(EmptyBoxState::BOX_EMPTY)
{
   TBOX_ASSERT(box_id.getPeriodicId().isValid());
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and
   // reset the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif
}

/*
 ******************************************************************************
 * Construct Box from the components of a reference Box
 * and possibly changing the periodic shift.
 ******************************************************************************
 */
Box::Box(
   const Box& other,
   const PeriodicId& periodic_id,
   const IntVector& refinement_ratio,
   const PeriodicShiftCatalog& shift_catalog):
   d_lo(other.d_lo),
   d_hi(other.d_hi),
   d_block_id(other.getBlockId()),
   d_id(other.getLocalId(), other.getOwnerRank(),
        periodic_id),
   d_id_locked(false),
   d_empty_flag(other.d_empty_flag)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, other, refinement_ratio);
#ifdef BOX_TELEMETRY
   // Increment the cumulative constructed count, active box count and
   // reset the high water mark of active boxes if necessary.
   ++s_cumulative_constructed_ct;
   ++s_active_ct;
   if (s_active_ct > s_high_water) {
      s_high_water = s_active_ct;
   }
#endif

   const tbox::Dimension& dim(d_lo.getDim());

   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < shift_catalog.getNumberOfShifts());

   if (refinement_ratio > IntVector::getZero(dim)) {

      if (other.getPeriodicId() != shift_catalog.getZeroShiftNumber()) {
         // Undo the shift that existed in other's Box.
         shift(-shift_catalog.shiftNumberToShiftDistance(other.
               getPeriodicId())
            * refinement_ratio);
      }

      if (periodic_id != shift_catalog.getZeroShiftNumber()) {
         // Apply the shift for this Box.
         shift(shift_catalog.shiftNumberToShiftDistance(periodic_id)
            * refinement_ratio);
      }

   } else if (refinement_ratio < IntVector::getZero(dim)) {

      if (other.getPeriodicId() != shift_catalog.getZeroShiftNumber()) {
         // Undo the shift that existed in other's Box.
         shift(shift_catalog.shiftNumberToShiftDistance(other.getPeriodicId())
            / refinement_ratio);
      }

      if (periodic_id != shift_catalog.getZeroShiftNumber()) {
         // Apply the shift for this Box.
         shift(-shift_catalog.shiftNumberToShiftDistance(periodic_id)
            / refinement_ratio);
      }

   } else {

      TBOX_ERROR(
         "Box::Box: Invalid refinement ratio "
         << refinement_ratio
         << "\nRefinement ratio must be completely positive or negative.");

   }
}

Box::~Box()
{
#ifdef BOX_TELEMETRY
   // There is one fewer box so decrement the active count.
   --s_active_ct;
#endif
}

Box&
Box::operator = (
   const Box& rhs)
{
   if (this != &rhs) {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, rhs);

      d_lo = rhs.d_lo;
      d_hi = rhs.d_hi;
      d_empty_flag = rhs.d_empty_flag;
      if (!d_id_locked) {
         d_block_id = rhs.d_block_id;
         d_id = rhs.d_id;
      } else if (d_block_id == rhs.d_block_id && d_id == rhs.d_id) {
         //No operation needed, the id objects were already equal. 
      } else {
         TBOX_ERROR("Attempted to change BoxId that is locked in an ordered BoxContainer.");
      }
#ifdef BOX_TELEMETRY
      // Increment the cumulative assigned count only.
      ++s_cumulative_assigned_ct;
#endif
   }
   return *this;
}

/*
 *******************************************************************************
 * Construct Box from a reference Box and possibly
 * changing the periodic shift.
 *
 * This method is not inlined because initializing a periodic-shifted
 * Box from another (possibly shifted) Box is more involved
 * and less frequently used.
 *
 * We inititalize d_id last so that we can support inititalizing an
 * object from a reference to itself.
 *******************************************************************************
 */
void
Box::initialize(
   const Box& other,
   const PeriodicId& periodic_id,
   const IntVector& refinement_ratio,
   const PeriodicShiftCatalog& shift_catalog)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, other, refinement_ratio);

   const tbox::Dimension& dim(d_lo.getDim());

   TBOX_ASSERT(periodic_id.isValid());
   TBOX_ASSERT(periodic_id.getPeriodicValue() < shift_catalog.getNumberOfShifts());

   d_lo = other.d_lo;
   d_hi = other.d_hi;
   d_empty_flag = other.d_empty_flag;

   if (refinement_ratio > IntVector::getZero(dim)) {

      if (other.getPeriodicId() != shift_catalog.getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         shift(-shift_catalog.shiftNumberToShiftDistance(other.
               getPeriodicId())
            * refinement_ratio);
      }

      if (periodic_id != shift_catalog.getZeroShiftNumber()) {
         // Apply the shift for this Box.
         shift(shift_catalog.shiftNumberToShiftDistance(periodic_id)
            * refinement_ratio);
      }

   } else if (refinement_ratio < IntVector::getZero(dim)) {

      if (other.getPeriodicId() != shift_catalog.getZeroShiftNumber()) {
         // Undo the shift that existed in r's Box.
         shift(shift_catalog.shiftNumberToShiftDistance(other.getPeriodicId())
            / refinement_ratio);
      }

      if (periodic_id != shift_catalog.getZeroShiftNumber()) {
         // Apply the shift for this Box.
         shift(-shift_catalog.shiftNumberToShiftDistance(periodic_id)
            / refinement_ratio);
      }

   } else {

      TBOX_ERROR(
         "Box::initialize: Invalid refinement ratio "
         << refinement_ratio
         << "\nRefinement ratio must be completely positive or negative.");

   }

   d_block_id = other.getBlockId();

   if (!d_id_locked) {
      d_id.initialize(
         other.getLocalId(), other.getOwnerRank(),
         periodic_id);
   } else {
      TBOX_ERROR("Attempted to change BoxId that is locked in an ordered BoxContainer.");
   }
}

Index
Box::index(
   const size_t offset) const
{
   TBOX_ASSERT(offset <= size());

   IntVector n(getDim());
   IntVector index(getDim());

   n = numberCells();

   size_t remainder = offset;

   for (int d = getDim().getValue() - 1; d > -1; --d) {
      /* Compute the stride for indexing */
      int stride = 1;
      for (int stride_dim = 0; stride_dim < d; ++stride_dim) {
         stride *= n[stride_dim];
      }

      /* Compute the local index */
      index[d] = static_cast<int>(remainder / stride);
      remainder -= index[d] * stride;

      /* Compute the global index */
      index[d] += lower(static_cast<tbox::Dimension::dir_t>(d));
   }

   return Index(index);
}

bool
Box::contains(
   const Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   if (box.empty()) {
      return true;
   }

   if (box.getBlockId() != d_block_id) {
      return false;
   }

   if (!contains(box.lower())) {
      return false;
   }

   if (!contains(box.upper())) {
      return false;
   }

   return true;
}

Box
Box::operator * (
   const Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   Box both(*this);

   if (d_block_id != box.d_block_id) {
      if (empty() || box.empty()) {
         both.setEmpty();
      } else {
         TBOX_ERROR("Attempted intersection of Boxes from different blocks.");
      }
   } else {
      both.d_lo.max(box.d_lo);
      both.d_hi.min(box.d_hi);
      both.d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
   }

   return both;
}

Box&
Box::operator *= (
   const Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   if (d_block_id != box.d_block_id) {
      if (empty() || box.empty()) {
         setEmpty();
      } else {
         TBOX_ERROR("Attempted intersection of Boxes from different blocks.");
      }
   } else {
      d_lo.max(box.d_lo);
      d_hi.min(box.d_hi);
      d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
   }

   return *this;
}

void
Box::intersect(
   const Box& other,
   Box& result) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, other);

   result = *this;
   if (result.d_block_id != other.d_block_id) {
      if (result.empty() || other.empty()) {
         result.setEmpty();
      } else {
         TBOX_ERROR("Attempted intersection of Boxes from different blocks.");
      }
   } else {
      result.d_lo.max(other.d_lo);
      result.d_hi.min(other.d_hi);
      result.d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
   }
}

bool
Box::intersects(
   const Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   if (d_block_id != box.d_block_id) {
      if (empty() || box.empty()) {
         return false;
      } else {
         TBOX_ERROR("Attempted intersection of Boxes from different blocks.");
      }
   }

   for (dir_t i = 0; i < getDim().getValue(); ++i) {
      if (tbox::MathUtilities<int>::Max(d_lo(i), box.d_lo(i)) >
          tbox::MathUtilities<int>::Min(d_hi(i), box.d_hi(i))) {
         return false;
      }
   }

   return true;
}

Box
Box::operator + (
   const Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   if (d_block_id != box.d_block_id) {
      if (!(empty()) || !(box.empty())) {
         TBOX_ERROR("Attempted bounding box of Boxes from different blocks.");
      }
   }

   Box bbox(*this);
   bbox += box;
   return bbox;
}

Box&
Box::operator += (
   const Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   if (!box.empty()) {
      if (empty()) {
         *this = box;
      } else if (d_block_id == box.d_block_id) {
         d_lo.min(box.d_lo);
         d_hi.max(box.d_hi);
         d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
      } else {
         TBOX_ERROR("Attempted bounding box of Boxes from different blocks.");
      }
   }
   return *this;
}

void
Box::lengthen(
   const dir_t direction,
   const int ghosts)
{
   TBOX_ASSERT((direction < getDim().getValue()));

   if (!empty()) {
      if (ghosts > 0) {
         d_hi(direction) += ghosts;
      } else {
         d_lo(direction) += ghosts;
      }
   }
}

void
Box::shorten(
   const dir_t direction,
   const int ghosts)
{
   TBOX_ASSERT((direction < getDim().getValue()));

   if (!empty()) {
      if (ghosts > 0) {
         d_hi(direction) -= ghosts;
      } else {
         d_lo(direction) -= ghosts;
      }
      d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
   }
}

void
Box::refine(
   const IntVector& ratio)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio);
   BlockId::block_t b = ratio.getNumBlocks() > 1 ? d_block_id.getBlockValue() : 0;
   TBOX_ASSERT(b < ratio.getNumBlocks());

   bool negative_ratio = false;
   for (unsigned int d = 0; d < getDim().getValue(); ++d) {
      if (ratio(b,d) < 0) {
         negative_ratio = true;
         break;
      }
   }

   if (!negative_ratio) {
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         d_lo(i) *= ratio(b,i);
         d_hi(i) = d_hi(i) * ratio(b,i) + (ratio(b,i) - 1);
      }
   } else {
      for (unsigned int i = 0; i < getDim().getValue(); ++i) {
         if (ratio(b,i) > 0) {
            d_lo(i) *= ratio(b,i);
            d_hi(i) = d_hi(i) * ratio(b,i) + (ratio(b,i) - 1);
         } else {
            d_lo(i) = coarsen(d_lo(i), -ratio(b,i));
            d_hi(i) = coarsen(d_hi(i), -ratio(b,i));
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Type Conversions
 *
 *************************************************************************
 */

tbox::DatabaseBox
Box::DatabaseBox_from_Box() const
{
   tbox::DatabaseBox new_Box;

   new_Box.setDim(getDim());

   for (dir_t i = 0; i < getDim().getValue(); ++i) {
      new_Box.lower(i) = d_lo(i);
      new_Box.upper(i) = d_hi(i);
   }

   return new_Box;
}

void
Box::set_Box_from_DatabaseBox(
   const tbox::DatabaseBox& box)
{
   for (dir_t i = 0; i < box.getDimVal(); ++i) {
      d_lo(i) = box.lower(i);
      d_hi(i) = box.upper(i);
   }
   d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
}

void
Box::putToIntBuffer(
   int* buffer) const
{
   buffer[0] = static_cast<int>(d_block_id.getBlockValue());
   ++buffer;

   d_id.putToIntBuffer(buffer);
   buffer += BoxId::commBufferSize();

   const dir_t dim(d_lo.getDim().getValue());
   for (dir_t d = 0; d < dim; ++d) {
      buffer[d] = d_lo(d);
      buffer[dim + d] = d_hi(d);
   }

}

void
Box::getFromIntBuffer(
   const int* buffer)
{
   d_block_id = BlockId(buffer[0]);
   ++buffer;

   d_id.getFromIntBuffer(buffer);
   buffer += BoxId::commBufferSize();

   const dir_t dim(d_lo.getDim().getValue());
   for (dir_t d = 0; d < dim; ++d) {
      d_lo(d) = buffer[d];
      d_hi(d) = buffer[dim + d];
   }
   d_empty_flag = EmptyBoxState::BOX_UNKNOWN;

}

void
Box::putToMessageStream(
   tbox::MessageStream& msg) const
{
   msg << d_block_id.getBlockValue();
   d_id.putToMessageStream(msg);
   msg.pack(&d_lo[0], d_lo.getDim().getValue());
   msg.pack(&d_hi[0], d_hi.getDim().getValue());
}

void
Box::getFromMessageStream(
   tbox::MessageStream& msg)
{
   int tmpi;
   msg >> tmpi;
   d_block_id = BlockId(tmpi);
   d_id.getFromMessageStream(msg);
   msg.unpack(&d_lo[0], d_lo.getDim().getValue());
   msg.unpack(&d_hi[0], d_hi.getDim().getValue());
   d_empty_flag = EmptyBoxState::BOX_UNKNOWN;
}

/*
 *************************************************************************
 *
 * Stream input/output operators: [(l0,...,ln),(u0,...,un)].8
 *
 *************************************************************************
 */

std::istream&
operator >> (
   std::istream& s,
   Box& box)
{
   while (s.get() != '[') ;
   Index tmp(box.getDim());
   s >> tmp;
   box.setLower(tmp);
   while (s.get() != ',') NULL_STATEMENT;
   s >> tmp;
   box.setUpper(tmp);
   while (s.get() != ']') NULL_STATEMENT;
   return s;
}

std::ostream&
operator << (
   std::ostream& s,
   const Box& box)
{
   if (box.empty()) {
      s << "[(),()]";
   } else {
      s << box.getBoxId() << ' ' << box.getBlockId()
      << '[' << box.lower() << ',' << box.upper() << ']';
   }
   return s;
}

/*
 *************************************************************************
 *
 * Static member function called from coalesceWith().  It attempts to
 * recursively coalesce intervals individual directions in index space.
 * If it is possible to coalesce two intervals (defined by a proper
 * overlap or adjacency relationship), the value true is returned.
 * If this is impossible, false is returned.
 *
 *************************************************************************
 */

bool
Box::coalesceIntervals(
   const int* lo1,
   const int* hi1,
   const int* lo2,
   const int* hi2,
   const int dim)
{
   bool retval = false;
   if (dim == 1) {
      // interval 1 to the right of interval 2.
      if ((lo1[0] <= hi2[0] + 1) && (hi2[0] <= hi1[0])) {
         retval = true;
         return retval;
      }
      // interval 1 to the left of interval 2.
      if ((lo1[0] <= lo2[0]) && (lo2[0] <= hi1[0] + 1)) {
         retval = true;
         return retval;
      }
   } else {
      for (dir_t id = 0; id < dim; ++id) {
         if ((lo1[id] == lo2[id]) && (hi1[id] == hi2[id])) {
            dir_t id2;
            int low1[SAMRAI::MAX_DIM_VAL];
            int high1[SAMRAI::MAX_DIM_VAL];
            int low2[SAMRAI::MAX_DIM_VAL];
            int high2[SAMRAI::MAX_DIM_VAL];
            for (id2 = 0; id2 < id; ++id2) {
               low1[id2] = lo1[id2];
               high1[id2] = hi1[id2];
               low2[id2] = lo2[id2];
               high2[id2] = hi2[id2];
            }
            for (id2 = static_cast<tbox::Dimension::dir_t>(id + 1); id2 < dim; ++id2) {
               dir_t id1 = static_cast<tbox::Dimension::dir_t>(id2 - 1);
               low1[id1] = lo1[id2];
               high1[id1] = hi1[id2];
               low2[id1] = lo2[id2];
               high2[id1] = hi2[id2];
            }
            if (coalesceIntervals(low1, high1, low2, high2, dim - 1)) {
               retval = true;
               return retval;
            }
         }
      }
   }

   return retval;
}

/*
 *************************************************************************
 *
 * Return true if this box can be coalesced with the argument box,
 * and set this box to the union of the boxes.  Otherwise, return false
 * and leave this box as is.  Two boxes may be coalesced if their union
 * is a box.  This routine attempts to coalesce the boxes along
 * each coordinate direction using the coalesceIntervals() function.
 *
 *************************************************************************
 */

bool
Box::coalesceWith(
   const Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   bool retval = false;

   if (empty() || box.empty()) {
      retval = true;
      *this += box;
   } else if (d_block_id == box.d_block_id) {
      dir_t id;
      const int* box_lo = &box.lower()[0];
      const int* box_hi = &box.upper()[0];
      int me_lo[SAMRAI::MAX_DIM_VAL];
      int me_hi[SAMRAI::MAX_DIM_VAL];
      for (id = 0; id < getDim().getValue(); ++id) {
         me_lo[id] = d_lo(id);
         me_hi[id] = d_hi(id);
      }
      if (coalesceIntervals(box_lo, box_hi, me_lo, me_hi, getDim().getValue())) {
         retval = true;
      } else { // test for one box containing the other...
         // test whether me contains box.
         retval = true;
         id = 0;
         while (retval && (id < getDim().getValue())) {
            retval = ((me_lo[id] <= box_lo[id]) && (me_hi[id] >= box_hi[id]));
            ++id;
         }
         if (!retval) { // me doesn't contain box; check other way around...
            retval = true;
            id = 0;
            while (retval && (id < getDim().getValue())) {
               retval = ((box_lo[id] <= me_lo[id])
                         && (box_hi[id] >= me_hi[id]));
               ++id;
            }
         }
      }
   } else { // BlockIds don't match, so don't coalesce.
      retval = false;
   }

   if (retval) *this += box;

   return retval;
}

/*
 *************************************************************************
 *                                                                       *
 * Rotates a 3-Dimensional box 90*num_rotations degrees around the given *
 * axis.                                                                 *
 *                                                                       *
 *************************************************************************
 */

void
Box::rotateAboutAxis(
   const dir_t axis,
   const int num_rotations)
{
   TBOX_ASSERT(axis < getDim().getValue());
   TBOX_ASSERT(getDim().getValue() == 3);

   const tbox::Dimension& dim(getDim());

   const int a = (axis + 1) % dim.getValue();
   const int b = (axis + 2) % dim.getValue();

   Index tmp_lo(dim);
   Index tmp_hi(dim);

   for (int j = 0; j < num_rotations; ++j) {
      tmp_lo = d_lo;
      tmp_hi = d_hi;
      d_lo(a) = tmp_lo(b);
      d_lo(b) = -tmp_hi(a) - 1;
      d_hi(a) = tmp_hi(b);
      d_hi(b) = -tmp_lo(a) - 1;
   }
}

/*
 *************************************************************************
 *                                                                       *
 * Rotate a box in the manner determined by the rotation number          *
 *                                                                       *
 *************************************************************************
 */

void
Box::rotate(
   const Transformation::RotationIdentifier rotation_ident)
{
   if (rotation_ident == Transformation::NO_ROTATE)
      return;

   TBOX_ASSERT(getDim().getValue() == 1 || getDim().getValue() == 2 ||
      getDim().getValue() == 3);

   if (getDim().getValue() == 1) {
      int rotation_number = static_cast<int>(rotation_ident);
      if (rotation_number > 1) {
         TBOX_ERROR("Box::rotate invalid 1D RotationIdentifier.");
      }
      if (rotation_number) {
         Index tmp_lo(d_lo);
         Index tmp_hi(d_hi);
         d_lo(0) = -tmp_hi(0) - 1;
         d_hi(0) = -tmp_lo(0) - 1;
      }
   } else if (getDim().getValue() == 2) {
      int rotation_number = static_cast<int>(rotation_ident);
      if (rotation_number > 3) {
         TBOX_ERROR("Box::rotate invalid 2D RotationIdentifier.");
      }
      for (int j = 0; j < rotation_number; ++j) {
         Index tmp_lo(d_lo);
         Index tmp_hi(d_hi);

         d_lo(0) = tmp_lo(1);
         d_lo(1) = -tmp_hi(0) - 1;
         d_hi(0) = tmp_hi(1);
         d_hi(1) = -tmp_lo(0) - 1;
      }
   } else {

      if (getDim().getValue() == 3) {
         if (rotation_ident == Transformation::NO_ROTATE) {
            return;
         } else if (rotation_ident == Transformation::KUP_IUP_JUP) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(2, 3);
         } else if (rotation_ident == Transformation::JUP_KUP_IUP) {
            rotateAboutAxis(1, 1);
            rotateAboutAxis(2, 1);
         } else if (rotation_ident == Transformation::IDOWN_KUP_JUP) {
            rotateAboutAxis(1, 2);
            rotateAboutAxis(0, 3);
         } else if (rotation_ident == Transformation::KUP_JUP_IDOWN) {
            rotateAboutAxis(1, 3);
         } else if (rotation_ident == Transformation::JUP_IDOWN_KUP) {
            rotateAboutAxis(2, 1);
         } else if (rotation_ident == Transformation::KDOWN_JUP_IUP) {
            rotateAboutAxis(1, 1);
         } else if (rotation_ident == Transformation::IUP_KDOWN_JUP) {
            rotateAboutAxis(0, 3);
         } else if (rotation_ident == Transformation::JUP_IUP_KDOWN) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(2, 3);
         } else if (rotation_ident == Transformation::KDOWN_IDOWN_JUP) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(2, 1);
         } else if (rotation_ident == Transformation::IDOWN_JUP_KDOWN) {
            rotateAboutAxis(1, 2);
         } else if (rotation_ident == Transformation::JUP_KDOWN_IDOWN) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(1, 3);
         } else if (rotation_ident == Transformation::JDOWN_IUP_KUP) {
            rotateAboutAxis(2, 3);
         } else if (rotation_ident == Transformation::IUP_KUP_JDOWN) {
            rotateAboutAxis(0, 1);
         } else if (rotation_ident == Transformation::KUP_JDOWN_IUP) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 1);
         } else if (rotation_ident == Transformation::JDOWN_KUP_IDOWN) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(1, 3);
         } else if (rotation_ident == Transformation::IDOWN_JDOWN_KUP) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 2);
         } else if (rotation_ident == Transformation::KUP_IDOWN_JDOWN) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(2, 1);
         } else if (rotation_ident == Transformation::JDOWN_KDOWN_IUP) {
            rotateAboutAxis(0, 3);
            rotateAboutAxis(1, 1);
         } else if (rotation_ident == Transformation::KDOWN_IUP_JDOWN) {
            rotateAboutAxis(0, 1);
            rotateAboutAxis(2, 3);
         } else if (rotation_ident == Transformation::IUP_JDOWN_KDOWN) {
            rotateAboutAxis(0, 2);
         } else if (rotation_ident == Transformation::JDOWN_IDOWN_KDOWN) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(2, 1);
         } else if (rotation_ident == Transformation::KDOWN_JDOWN_IDOWN) {
            rotateAboutAxis(0, 2);
            rotateAboutAxis(1, 3);
         } else if (rotation_ident == Transformation::IDOWN_KDOWN_JDOWN) {
            rotateAboutAxis(1, 2);
            rotateAboutAxis(0, 1);
         } else {
            TBOX_ERROR("Box::rotate invalid 3D RotationIdentifier.");
         }
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
Box::initializeCallback()
{
   for (unsigned short d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      tbox::Dimension dim(static_cast<unsigned short>(d + 1));
      s_emptys[d].emplace(dim);

      /*
       * Note we can't use Index getMin, getMax here as that
       * would create a dependency between static initializers
       */
      s_universes[d].emplace(
            Index(dim, tbox::MathUtilities<int>::getMin()),
            Index(dim, tbox::MathUtilities<int>::getMax()),
            BlockId(0));
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
Box::finalizeCallback()
{
   for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
      s_emptys[d].reset();
      s_universes[d].reset();
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
BoxIterator
Box::begin() const
{
   return iterator(*this, true);
}

/*
 *************************************************************************
 *************************************************************************
 */
BoxIterator
Box::end() const
{
   return iterator(*this, false);
}

BoxIterator::BoxIterator(
   const Box& box,
   bool begin):
   d_index(box.lower()),
   d_box(box)
{
   if (!d_box.empty() && !begin) {
      d_index(d_box.getDim().getValue() - 1) =
         d_box.upper(static_cast<tbox::Dimension::dir_t>(d_box.getDim().getValue() - 1)) + 1;
   }
}

}
}
