/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a Connector.
 *
 ************************************************************************/
#ifndef included_hier_ConnectorStatistics
#define included_hier_ConnectorStatistics

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief A utility for writing out various statistics of Connectors.
 */
class ConnectorStatistics
{

public:
   /*!
    * @brief Constructor.
    *
    * Compute and store statistics for the given Connector.  The
    * statistics reflects the current Connector state and can be
    * printed out afterwards.  All processes in the Connector's
    * SAMRAI_MPI must call this constructor because it requires
    * collective communication.
    *
    * @param[in] connector
    *
    * @pre connector.isFinalized()
    */
   explicit ConnectorStatistics(
      const Connector& connector);

   /*!
    * @brief Print out local and globally reduced statistics on the
    * Boxes.
    *
    * @param[in,out] os The output stream
    *
    * @param[in] border A string to print at the start of every line
    * in the output.
    */
   void
   printNeighborStats(
      std::ostream& os,
      const std::string& border) const;

private:
   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Indices for statistical quantites.
    *
    * Relationships (or edges) are from the perspective of the base
    * boxes, so it is possible for two base boxes to have
    * relationships to the same head box, in which case, we count two
    * relationships.  Neighbors are distinct head boxes, so in the
    * same example, there is only one neighbors even though there are
    * two relationships to that neighbor.
    *
    * The MIN_* and MAX_* quantities refer to the min/max over the
    * local base boxes.  When we compute the global min/max of these
    * quantities, that is exactly what they are: min of min, min of
    * max, max of min and max of max.
    */
   enum { NUMBER_OF_BASE_BOXES,
          NUMBER_OF_BASE_CELLS,

          HAS_ANY_NEIGHBOR_SETS,
          NUMBER_OF_NEIGHBOR_SETS,

          HAS_ANY_RELATIONSHIPS,
          NUMBER_OF_RELATIONSHIPS,
          MIN_NUMBER_OF_RELATIONSHIPS,
          MAX_NUMBER_OF_RELATIONSHIPS,

          NUMBER_OF_NEIGHBORS,
          NUMBER_OF_LOCAL_NEIGHBORS,
          NUMBER_OF_REMOTE_NEIGHBORS,
          NUMBER_OF_REMOTE_NEIGHBOR_OWNERS,

          OVERLAP_SIZE,
          LOCAL_OVERLAP_SIZE,
          REMOTE_OVERLAP_SIZE,

          NUMBER_OF_QUANTITIES };

   /*
    * @brief StatisticalQuantities to compute the min, avg and max for.
    *
    * These quantities will be computed locally on each process and
    * globally reduced.  Not all of these quantities are floating
    * points but all are represented as such.
    */
   struct StatisticalQuantities {
      StatisticalQuantities();
      double d_values[NUMBER_OF_QUANTITIES];
   };

   void
   computeLocalConnectorStatistics(
      const Connector& connector);

   void
   reduceStatistics();

   tbox::SAMRAI_MPI d_mpi;

   //! @brief Statistics of local process.
   StatisticalQuantities d_sq;
   //! @brief Global min of d_sq.
   StatisticalQuantities d_sq_min;
   //! @brief Global max of d_sq.
   StatisticalQuantities d_sq_max;
   //! @brief Global sum of d_sq.
   StatisticalQuantities d_sq_sum;

   int d_rank_of_min[NUMBER_OF_QUANTITIES];
   int d_rank_of_max[NUMBER_OF_QUANTITIES];

   /*!
    * @brief Names of the quantities in StatisticalQuantities.
    */
   static std::string s_quantity_names[NUMBER_OF_QUANTITIES];

   /*!
    * @brief Longest length in s_quantity_names.
    */
   static int s_longest_length;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif  // included_hier_ConnectorStatistics
