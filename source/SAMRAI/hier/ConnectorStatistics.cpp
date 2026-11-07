/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a Connector.
 *
 ************************************************************************/
#include "SAMRAI/hier/ConnectorStatistics.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

std::string ConnectorStatistics::s_quantity_names[NUMBER_OF_QUANTITIES];
int ConnectorStatistics::s_longest_length;

tbox::StartupShutdownManager::Handler
ConnectorStatistics::s_initialize_finalize_handler(
   ConnectorStatistics::initializeCallback,
   0,
   0,
   ConnectorStatistics::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ************************************************************************
 * Constructor.
 ************************************************************************
 */
ConnectorStatistics::ConnectorStatistics(
   const Connector& connector):
   d_mpi(connector.getMPI())
{
   if (!connector.isFinalized()) {
      TBOX_ERROR("ConnectorStatistics requires an finalized Connector.");
   }

   computeLocalConnectorStatistics(connector);
   reduceStatistics();
}

/*
 ************************************************************************
 ************************************************************************
 */
ConnectorStatistics::StatisticalQuantities::StatisticalQuantities()
{
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      d_values[i] = 0.;
   }
}

/*
 ***********************************************************************
 * Compute statistics for local process.
 ***********************************************************************
 */

void
ConnectorStatistics::computeLocalConnectorStatistics(const Connector& connector)
{
   if (!connector.isFinalized()) {
      TBOX_ERROR("ConnectorStatistics::computeLocalStatistics cannot compute\n"
         << "statistics for unfinalized Connector.");
      return;
   }

   connector.cacheGlobalReducedData();

   const BoxLevel& base(connector.getBase());
   const tbox::SAMRAI_MPI& mpi(base.getMPI());
   const tbox::Dimension& dim(base.getDim());

   /*
    * Whether head boxes should be refined or coarsened before box
    * calculus.
    */
   const bool refine_head = connector.getHeadCoarserFlag() &&
      (connector.getRatio() != IntVector::getOne(dim) ||
       !connector.ratioIsExact());
   const bool coarsen_head = !connector.getHeadCoarserFlag() &&
      (connector.getRatio() != IntVector::getOne(dim) ||
       !connector.ratioIsExact());

   /*
    * Compute per-processor statistics.  Some quantities are readily
    * available while others are computed in the loops following.
    */

   d_sq.d_values[NUMBER_OF_BASE_BOXES] =
      static_cast<double>(base.getLocalNumberOfBoxes());
   d_sq.d_values[NUMBER_OF_BASE_CELLS] =
      static_cast<double>(base.getLocalNumberOfCells());

   d_sq.d_values[HAS_ANY_NEIGHBOR_SETS] =
      static_cast<double>(connector.getLocalNumberOfNeighborSets() > 0);
   d_sq.d_values[NUMBER_OF_NEIGHBOR_SETS] =
      static_cast<double>(connector.getLocalNumberOfNeighborSets());

   d_sq.d_values[HAS_ANY_RELATIONSHIPS] =
      static_cast<double>(connector.getLocalNumberOfRelationships() > 0);
   d_sq.d_values[NUMBER_OF_RELATIONSHIPS] =
      static_cast<double>(connector.getLocalNumberOfRelationships());
   d_sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS] =
      tbox::MathUtilities<double>::getMax();
   d_sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS] = 0.;

   d_sq.d_values[NUMBER_OF_NEIGHBORS] = 0.;
   d_sq.d_values[NUMBER_OF_LOCAL_NEIGHBORS] = 0.;
   d_sq.d_values[NUMBER_OF_REMOTE_NEIGHBORS] = 0.;
   d_sq.d_values[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] = 0.;

   BoxContainer visible_neighbors(true); // All neighbors of local base boxes.

   for (Connector::ConstNeighborhoodIterator nbi = connector.begin();
        nbi != connector.end(); ++nbi) {

      const int num_relations = connector.numLocalNeighbors(*nbi);

      d_sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS] =
         tbox::MathUtilities<double>::Min(
            d_sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS],
            static_cast<double>(num_relations));

      d_sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS] =
         tbox::MathUtilities<double>::Max(
            d_sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS],
            static_cast<double>(num_relations));

      Box base_box = *base.getBoxStrict(*nbi);
      base_box.grow(connector.getConnectorWidth());

      for (Connector::ConstNeighborIterator ni = connector.begin(nbi);
           ni != connector.end(nbi); ++ni) {

         visible_neighbors.insert(*ni);
         Box neighbor = *ni;
         if (refine_head) {
            neighbor.refine(connector.getRatio());
         } else if (coarsen_head) {
            neighbor.coarsen(connector.getRatio());
         }
         if (neighbor.getBlockId() != base_box.getBlockId()) {
            base.getGridGeometry()->transformBox(neighbor,
               base.getRefinementRatio(),
               base_box.getBlockId(),
               neighbor.getBlockId());
         }
         neighbor *= base_box;
         const size_t size = neighbor.size();

         d_sq.d_values[OVERLAP_SIZE] += static_cast<double>(size);
         if (neighbor.getOwnerRank() == mpi.getRank()) {
            d_sq.d_values[LOCAL_OVERLAP_SIZE] += static_cast<double>(size);
         } else {
            d_sq.d_values[REMOTE_OVERLAP_SIZE] += static_cast<double>(size);
         }

      }

   }

   d_sq.d_values[NUMBER_OF_NEIGHBORS] =
      static_cast<double>(visible_neighbors.size());

   std::set<int> remote_neighbor_owners;
   for (BoxContainer::const_iterator bi = visible_neighbors.begin();
        bi != visible_neighbors.end(); ++bi) {
      const Box& neighbor = *bi;

      d_sq.d_values[NUMBER_OF_LOCAL_NEIGHBORS] +=
         neighbor.getOwnerRank() == mpi.getRank();

      d_sq.d_values[NUMBER_OF_REMOTE_NEIGHBORS] +=
         neighbor.getOwnerRank() != mpi.getRank();

      remote_neighbor_owners.insert(neighbor.getOwnerRank());

   }
   remote_neighbor_owners.erase(mpi.getRank());
   d_sq.d_values[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] =
      static_cast<double>(remote_neighbor_owners.size());
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void ConnectorStatistics::reduceStatistics()
{
   d_sq_min = d_sq;
   d_sq_max = d_sq;
   d_sq_sum = d_sq;

   if (d_mpi.getSize() > 1) {
      d_mpi.AllReduce(
         d_sq_min.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MINLOC,
         d_rank_of_min);
      d_mpi.AllReduce(
         d_sq_max.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MAXLOC,
         d_rank_of_max);
      d_mpi.AllReduce(d_sq_sum.d_values, NUMBER_OF_QUANTITIES, MPI_SUM);
   } else {
      for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
         d_rank_of_min[i] = d_rank_of_max[i] = 0;
      }
   }
}

/*
 ***********************************************************************
 * Write out local and globally reduced statistics on the relationships.
 ***********************************************************************
 */

void ConnectorStatistics::printNeighborStats(
   std::ostream& co,
   const std::string& border) const
{
   co.unsetf(std::ios::fixed | std::ios::scientific);
   co.precision(3);

   co << border << "N = " << d_sq_sum.d_values[NUMBER_OF_BASE_BOXES]
      << " (global number of boxes),  "
      << "P = " << d_mpi.getSize() << " (number of processes)\n"
      << border << std::setw(s_longest_length) << std::string()
      << "    local        min               max             sum    sum/N    sum/P\n";

   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      co << border << std::setw(s_longest_length) << std::left
         << s_quantity_names[i]
         << ' ' << std::setw(8) << std::right << d_sq.d_values[i]
         << ' ' << std::setw(8) << std::right << d_sq_min.d_values[i] << " @ "
         << std::setw(6) << std::left << d_rank_of_min[i]
         << ' ' << std::setw(8) << std::right << d_sq_max.d_values[i] << " @ "
         << std::setw(6) << std::left << d_rank_of_max[i]
         << ' ' << std::setw(8) << std::right << d_sq_sum.d_values[i]
         << ' ' << std::setw(8) << std::right
         << d_sq_sum.d_values[i] / d_sq_sum.d_values[NUMBER_OF_BASE_BOXES]
         << ' ' << std::setw(8) << std::right
         << d_sq_sum.d_values[i] / d_mpi.getSize() << '\n';
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
ConnectorStatistics::initializeCallback()
{

   s_quantity_names[NUMBER_OF_BASE_BOXES] = "num base boxes";
   s_quantity_names[NUMBER_OF_BASE_CELLS] = "num base cells";

   s_quantity_names[HAS_ANY_NEIGHBOR_SETS] = "has any neighbor sets";
   s_quantity_names[NUMBER_OF_NEIGHBOR_SETS] = "num neighbor sets";

   s_quantity_names[HAS_ANY_RELATIONSHIPS] = "has any relationships";
   s_quantity_names[NUMBER_OF_RELATIONSHIPS] = "num relationships";
   s_quantity_names[MIN_NUMBER_OF_RELATIONSHIPS] = "min num relationships";
   s_quantity_names[MAX_NUMBER_OF_RELATIONSHIPS] = "max num relationships";

   s_quantity_names[NUMBER_OF_NEIGHBORS] = "num neighbors";
   s_quantity_names[NUMBER_OF_LOCAL_NEIGHBORS] = "num local neighbors";
   s_quantity_names[NUMBER_OF_REMOTE_NEIGHBORS] = "num remote neighbors";
   s_quantity_names[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] = "num remote neighbor owners";

   s_quantity_names[OVERLAP_SIZE] = "overlap size";
   s_quantity_names[LOCAL_OVERLAP_SIZE] = "local overlap size";
   s_quantity_names[REMOTE_OVERLAP_SIZE] = "remote overlap size";

   s_longest_length = 0;
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      s_longest_length = tbox::MathUtilities<int>::Max(
            s_longest_length, static_cast<int>(s_quantity_names[i].length()));
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */

void
ConnectorStatistics::finalizeCallback()
{
   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      s_quantity_names[i].clear();
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
