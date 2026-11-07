/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   spatial key used for generating space-filling curves.
 *
 ************************************************************************/

#ifndef included_mesh_SpatialKey
#define included_mesh_SpatialKey

#include "SAMRAI/SAMRAI_config.h"

#include <iostream>

namespace SAMRAI {
namespace mesh {

/**
 * Class SpatialKey provides a mapping from coordinates in
 * an abstract index space to a key that can be used to order the
 * the points with a Morton-type space filling curve (See Gutman,
 * Dr Dobb's Journal, July 1999, pg 115-21 for an introduction).
 * This curve orders a set of points in space so that points near
 * each other in space are likely, but not guaranteed, to be near
 * each other in the ordering.
 */

class SpatialKey
{
public:
   /**
    * The default constructor creates a spatial key with zero value.
    */
   SpatialKey();

   /**
    * Create a spatial key from given index space coordinates and
    * level number.
    */
   explicit SpatialKey(
      const unsigned int i,
      const unsigned int j = 0,
      const unsigned int k = 0,
      const unsigned int level_num = 0);

   /**
    * Copy constructor for spatial key.
    */
   SpatialKey(
      const SpatialKey& spatial_key);

   /**
    * The destructor for a spatial key does nothing interesting.
    */
   ~SpatialKey();

   /**
    * Assignment operator for spatial key.
    */
   SpatialKey&
   operator = (
      const SpatialKey& spatial_key)
   {
      for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; ++i) {
         d_key[i] = spatial_key.d_key[i];
      }
      return *this;
   }

   /**
    * Return true if argument key is equal to this key.  Otherwise,
    * return false.
    */
   bool
   operator == (
      const SpatialKey& spatial_key) const
   {
      bool are_equal = true;
      for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; ++i) {
         if (d_key[i] != spatial_key.d_key[i]) {
            are_equal = false;
            break;
         }
      }
      return are_equal;
   }

   /**
    * Return true if argument key is not equal to this key.  Otherwise,
    * return false.
    */
   bool
   operator != (
      const SpatialKey& spatial_key) const
   {
      return !((*this) == spatial_key);
   }

   /**
    * Return true if this key is less than argument key.  Otherwise,
    * return false.
    */
   bool
   operator < (
      const SpatialKey& spatial_key) const;

   /**
    * Return true if this key is less than or equal to argument key.
    * Otherwise, return false.
    */
   bool
   operator <= (
      const SpatialKey& spatial_key) const
   {
      return ((*this) < spatial_key) || ((*this) == spatial_key);
   }

   /**
    * Return true if this key is greater than argument key.  Otherwise,
    * return false.
    */
   bool
   operator > (
      const SpatialKey& spatial_key) const
   {
      return !((*this) < spatial_key) && ((*this) != spatial_key);
   }

   /**
    * Return true if this key is greater than or equal to argument key.
    * Otherwise, return false.
    */
   bool
   operator >= (
      const SpatialKey& spatial_key) const
   {
      return ((*this) > spatial_key) || ((*this) == spatial_key);
   }

   /**
    * Set this key to zero key.
    */
   void
   setToZero()
   {
      for (int i = 0; i < NUM_COORDS_MIXED_FOR_SPATIAL_KEY; ++i) {
         d_key[i] = 0;
      }
   }

   /**
    * Set this key from the index space coordinates and the level number.
    *
    * The values of i, j, and k default to 0 to handle cases
    * where the dimensions of the problem is not 3.
    * The default value of level_num is also 0.
    */
   void
   setKey(
      const unsigned int i = 0,
      const unsigned int j = 0,
      const unsigned int k = 0,
      const unsigned int level_num = 0);

   /**
    * Write a spatial key to an output stream.  The spatial key is
    * output in hex to avoid the binary to decimal conversion of the
    * key.
    *
    * Note that the proper functioning of this method depends on
    * having 32 bits per integer (more specifically, 8 hex characters
    * per integer).
    *
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const SpatialKey& spatial_key);

private:
   /*
    * Static integer constant.
    */
   static const int NUM_COORDS_MIXED_FOR_SPATIAL_KEY = 4;

   /*
    * Static integer constant.
    */
   static const int BITS_PER_BYTE;

   /*
    * Static integer constant.
    */
   static const int BITS_PER_HEX_CHAR;

   /*
    * Mix in one index space coordinate into the spatial
    * key.  coord is the value of the coordinate, and
    * coord_offset refers which coordinate is being blended
    * in.  coord_offset values of 3,2,1 indicate the
    * i,j,k coordinates respectively.  A coord_offset value
    * of 0 indicate the level number.
    */
   void
   blendOneCoord(
      const unsigned int coord,
      const int coord_offset);

   size_t d_bits_per_int;
   unsigned int d_key[NUM_COORDS_MIXED_FOR_SPATIAL_KEY];
};

}
}

#endif
