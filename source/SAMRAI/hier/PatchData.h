/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for patch data objects
 *
 ************************************************************************/

#ifndef included_hier_PatchData
#define included_hier_PatchData

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/**
 * Class PatchData is a pure virtual base class for the data storage
 * defined over a box.  Patch data objects are generally contained within
 * a patch.  Patch data defines the abstract virtual functions for data
 * management and communication that must be supplied by the subclasses.
 * Subclasses implement the virtual functions as appropriate for the derived
 * type; for example, cell centered objects will copy data differently than
 * face centered objects.
 *
 * Patch data objects are created by patch data factories and associated
 * subclasses.  This separation into abstract factory and concrete
 * implementation subclasses facilitates the creation of new patch data
 * subtypes.  See the Design Patterns book for more details about the
 * Abstract Factory pattern.
 *
 * The copy and pack/unpack functions in the patch data object take
 * box overlap descriptions that describe the index space over which
 * data is to be copied or packed/unpacked.  Box overlaps are computed
 * by the box geometry classes, which are accessed through the patch
 * data factories.  Box geometry classes are created by the patch data
 * factories instead of the patch data objects because patch data objects
 * are distributed across memory with patches and therefore may not exist
 * on a particular processor.  Patch data factories are guaranteed to
 * exist on all processors independent of the patch-to-processor mapping.
 *
 * @see BoxOverlap
 * @see BoxGeometry
 * @see Patch
 * @see PatchDataFactory
 * @see PatchDescriptor
 */

class PatchData
{
public:
   /**
    * The constructor for a patch data object.  Patch data objects will
    * manage the interior box over which they are defined and the associated
    * ghost cell width.
    *
    * @pre domain.getDim() == ghosts.getDim()
    */
   PatchData(
      const Box& domain,
      const IntVector& ghosts);

   /**
    * The virtual destructor for a patch data object.
    */
   virtual ~PatchData() = default;

   /**
    * Return the box over which this patch data object is defined.  All
    * objects in the same patch are defined over the same box, although
    * the patch data objects may interpret how to allocate storage for
    * that box in different ways.
    */
   const Box&
   getBox() const
   {
      return d_box;
   }

   /**
    * Return the ghost cell box.  The ghost cell box is defined to be
    * the interior box grown by the ghost cell width.
    */
   const Box&
   getGhostBox() const
   {
      return d_ghost_box;
   }

   /**
    * Get the ghost cell width associated with this patch data object.
    */
   const IntVector&
   getGhostCellWidth() const
   {
      return d_ghosts;
   }

   /**
    * Set the simulation time stamp for the patch data type.  The simulation
    * time is initialized to zero when the patch data type is created.
    */
   void
   setTime(
      const double timestamp)
   {
      d_timestamp = timestamp;
   }

   /**
    * Get the simulation time stamp for the patch data type.
    */
   double
   getTime() const
   {
      return d_timestamp;
   }

   /**
    * A fast copy between the source and destination.  Data is copied from
    * the source into the destination where there is overlap in the underlying
    * index space.  The copy is performed on the interior plus the ghost cell
    * width (for both the source and destination).  If this copy does not
    * understand how to copy data from the argument, then copy2() is called
    * on the source object.
    */
   virtual void
   copy(
      const PatchData& src) = 0;

   /**
    * A fast copy between the source and destination.  Data is copied from
    * the source into the destination where there is overlap in the underlying
    * index space.  The copy is performed on the interior plus the ghost cell
    * width (for both the source and destination).  If this copy does not
    * understand how to copy data from the destination, then it may throw
    * an assertion (aka dump core in a failed assertion).
    */
   virtual void
   copy2(
      PatchData& dst) const = 0;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects.  If this member function
    * cannot complete the copy from source (e.g., if it doesn't understand
    * the type of source), then copy2() is called on the source object.
    */
   virtual void
   copy(
      const PatchData& src,
      const BoxOverlap& overlap) = 0;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects.  The default implementation
    * of this method will call packStream without the fuser argument.
    */
   virtual void
   copyFuseable(
      const PatchData& src,
      const BoxOverlap& overlap);

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects If this member function
    * cannot complete the copy from the destination, then it may throw an
    * assertion (aka dump core in a failed assertion).
    */
   virtual void
   copy2(
      PatchData& dst,
      const BoxOverlap& overlap) const = 0;

   /**
    * Determines whether the patch data subclass can estimate the necessary
    * stream size using only index space information.  The return value will
    * most likely be true for data types that are fixed size (such as doubles)
    * but will be false for complex data types that allocate their own storage
    * (such as lists of particles).  This routine is used to estimate whether
    * a processor can estimate space for incoming messages or whether it needs
    * to receive a message size from the sending processor.
    */
   virtual bool
   canEstimateStreamSizeFromBox() const = 0;

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.  This estimate must be an upper
    * bound on the size of the data in the actual message stream.
    * The upper bound should be close, however, since buffer space
    * will be allocated according to these values, and excess buffer
    * space will waste memory resources.
    */
   virtual size_t
   getDataStreamSize(
      const BoxOverlap& overlap) const = 0;

   /**
    * Pack data lying on the specified index set into the output stream.
    * See the abstract stream virtual base class for more information about
    * the packing operators defined for streams.
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const BoxOverlap& overlap) const = 0;

   /**
    * Pack data lying on the specified index set into the output stream using
    * the given KernelFuser. The default implementation of this method will
    * call packStream without the fuser argument. See the abstract stream
    * virtual base class for more information about the packing operators
    * defined for streams.
    */
   virtual void
   packStreamFuseable(
      tbox::MessageStream& stream,
      const BoxOverlap& overlap) const;

   /**
    * Unpack data from the message stream into the specified index set.
    * See the abstract stream virtual base class for more information about
    * the packing operators defined for streams.
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const BoxOverlap& overlap) = 0;

   /**
    * Unpack data from the message stream into the specified index set using
    * the given KernelFuser. The default implementation of this method will
    * call unpackStream without the fuser argument. See the abstract stream
    * virtual base class for more information about the packing operators
    * defined for streams.
    */
   virtual void
   unpackStreamFuseable(
      tbox::MessageStream& stream,
      const BoxOverlap& overlap);

   /**
    * Checks that class version and restart file version are equal.  If so,
    * reads in the data members common to all patch data types from restart
    * database.
    *
    * @pre restart_db
    */
   virtual void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /**
    * Writes out the class version number to the restart database.  Then,
    * writes the data members common to all patch data types to database.
    *
    * @pre restart_db
    */
   virtual void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_box.getDim();
   }

protected:
   /**
    * This protected method is used by concrete patch data subclasses
    * to set the ghost box over which the patch data will be allocated.
    * Note that this allows the ghost box to be inconsistant with its
    * standard interpretation as the patch domain box grown by the ghost
    * cell width (as set in the constructor).
    *
    * This function is included to treat some special cases for concrete
    * patch data types and should be used with caution.
    *
    * @pre getDim() == ghost_box.getDim()
    * @pre (ghost_box * getBox()).isSpatiallyEqual(getBox())
    */
   void
   setGhostBox(
      const Box& ghost_box)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(d_box, ghost_box);
      TBOX_ASSERT((ghost_box * d_box).isSpatiallyEqual(d_box));
      d_ghost_box = ghost_box;
   }

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_PATCH_DATA_VERSION;

   PatchData(
      const PatchData&) = delete;
   PatchData&
   operator = (
      const PatchData&) = delete;

   Box d_box;                           // interior box description
   Box d_ghost_box;                     // interior box plus ghosts
   IntVector d_ghosts;                  // ghost cell width
   double d_timestamp;                          // timestamp for the data

};

}
}

#endif
