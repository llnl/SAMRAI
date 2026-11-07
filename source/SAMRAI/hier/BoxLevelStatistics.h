/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevelStatistics
#define included_hier_BoxLevelStatistics

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxLevel.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief A utility for writing out various statistics of Boxes.
 */
class BoxLevelStatistics
{

public:
   /*!
    * @brief Constructor.
    *
    * Compute and store statistics for the given BoxLevel.  The
    * statistics reflects the current BoxLevel state and can be
    * printed out afterwards.  All processes in the BoxLevel's
    * SAMRAI_MPI must call this constructor because it requires
    * collective communication.
    *
    * @param[in] box_level BoxLevel to compute statistics for.
    *
    * @pre box_level.isInitialized()
    */
   explicit BoxLevelStatistics(
      const BoxLevel& box_level);

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
   printBoxStats(
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
    */
   enum { HAS_ANY_BOX,
          NUMBER_OF_CELLS,
          NUMBER_OF_BOXES,
          MAX_BOX_VOL,
          MIN_BOX_VOL,
          MAX_BOX_LEN,
          MIN_BOX_LEN,
          MAX_ASPECT_RATIO,
          SUM_ASPECT_RATIO,
          SUM_SURFACE_AREA,
          SUM_NORM_SURFACE_AREA,
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
   computeLocalBoxLevelStatistics(
      const BoxLevel& box_level);

   void
   reduceStatistics();

   tbox::SAMRAI_MPI d_mpi;

   const tbox::Dimension d_dim;

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

   static tbox::StartupShutdownManager::Handler s_initialize_finalize_handler;

};

}
}

#endif  // included_hier_BoxLevelStatistics
