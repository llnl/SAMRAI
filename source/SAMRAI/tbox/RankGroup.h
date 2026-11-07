/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A class to manage a group of processor ranks
 *
 ************************************************************************/

#ifndef included_tbox_RankGroup
#define included_tbox_RankGroup

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * Class RankGroup is used to manage a set of ranks, namely a subset of
 * the available ranks in a parallel run.
 *
 * This class allows a group of ranks to be constructed as a generalized
 * subset of ranks, represented by an array passed to the constructor, or
 * as a contiguous ordered subset of ranks, the range from a minimum to a
 * maximum rank value.  It also has a default case where it represents the
 * full set of available ranks, from 0 to P-1, P being the number of
 * available processors.
 *
 * If N is the number of ranks in the represented subset, the class gives the
 * user a 1-to-1 mapping between the integer set [0,N-1] and the ranks
 * in the subset.
 */

class RankGroup
{
public:
   /*!
    * The default constructor constructs a RankGroup representing the full
    * set of available ranks.  It will use all of the ranks in the
    * MPI communicator used by the SAMRAI MPI object.  If a different
    * communicator is desired, then use a constructor with a communicator
    * argument.
    */
   RankGroup();

   /*!
    * This constructor constructs a RankGroup representing the full
    * set of available ranks in the given communicator.
    */
   explicit RankGroup(
      const SAMRAI_MPI& d_samrai_mpi);

   /*!
    * This constructor creates a RankGroup consisting of all ranks from min
    * to max, inclusive.  min must be >= 0 and max must be less than the
    * total number of available processors.
    *
    * @pre min >= 0
    * @pre min <= max
    */
   RankGroup(
      const int min,
      const int max,
      const SAMRAI_MPI& samrai_mpi =
         SAMRAI_MPI(SAMRAI_MPI::getSAMRAIWorld()));

   /*!
    * This constructor creates a RankGroup consisting of ranks corresponding
    * to the integers in the vector.  Each member of the vector must be >= 0,
    * less than the total number of available processors, and unique within
    * the vector.  Due to the use of an vector for storage, RankGroups
    * created with this constructor should be expected to be less efficient
    * than those created with the above min/max constructor.
    *
    * @pre !rank_group.empty()
    */
   explicit RankGroup(
      const std::vector<int>& rank_group,
      const SAMRAI_MPI& samrai_mpi =
         SAMRAI_MPI(SAMRAI_MPI::getSAMRAIWorld()));

   /*!
    * Copy constructor.
    */
   RankGroup(
      const RankGroup& other);

   /*!
    * Destructor
    */
   ~RankGroup();

   /*!
    * Assignment operator.
    */
   RankGroup&
   operator = (
      const RankGroup& rhs);

   /*!
    * Returns true if the RankGroup contains ranks for all available
    * processors.
    */
   bool
   containsAllRanks() const
   {
      return d_storage == USING_ALL;
   }

   /*!
    * Set the minimum and maximum ranks for the RankGroup.  All previous
    * state of this object will be eliminated.  The restrictions on the
    * arguments are the same as for the constructor that takes the min/max
    * arguments.
    */
   void
   setMinMax(
      const int min,
      const int max)
   {
      TBOX_ASSERT(min >= 0);
      TBOX_ASSERT(min <= max);
      d_storage = USING_MIN_MAX;
      d_ranks.resize(0);
      d_min = min;
      d_max = max;
   }

   /*!
    * Return true if the given rank is contained in the RankGroup.
    */
   bool
   isMember(
      const int rank) const;

   /*!
    * Return the size of the subset of ranks represented by this RankGroup.
    */
   int
   size() const;

   /*!
    * Given an integer identifier from the set [0,N-1], N being the size of
    * the RankGroup, return a unique rank from the RankGroup according to
    * a 1-to-1 mapping.
    */
   int
   getMappedRank(
      const int index) const;

   /*!
    * This is the inverse function of getMappedRank.  Given a rank that is
    * contained in the RankGroup, return a unique integer identifier in the
    * set [0,N-1], N being the size of the rank group, according to a 1-to-1
    * mapping.
    *
    * @pre rank >= 0
    */
   int
   getMapIndex(
      const int rank) const;

private:
   enum StorageType { USING_ALL,
                      USING_ARRAY,
                      USING_MIN_MAX };

   int d_min;
   int d_max;

   std::vector<int> d_ranks;

   StorageType d_storage;

   SAMRAI_MPI d_samrai_mpi;
};

}
}

#endif
