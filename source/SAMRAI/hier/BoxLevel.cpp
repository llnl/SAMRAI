/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of boxes in a box_level of a distributed box graph.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxLevel.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"
#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int BoxLevel::HIER_BOX_LEVEL_VERSION = 0;
const int BoxLevel::BOX_LEVEL_NUMBER_OF_STATS = 20;

std::shared_ptr<tbox::Timer> BoxLevel::t_initialize_private;
std::shared_ptr<tbox::Timer> BoxLevel::t_acquire_remote_boxes;
std::shared_ptr<tbox::Timer> BoxLevel::t_cache_global_reduced_data;

const LocalId BoxLevel::s_negative_one_local_id(-1);

tbox::StartupShutdownManager::Handler
BoxLevel::s_initialize_finalize_handler(
   BoxLevel::initializeCallback,
   0,
   0,
   BoxLevel::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

BoxLevel::BoxLevel(
   const tbox::Dimension& dim,
   tbox::Database& restart_db,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom):

   d_mpi(tbox::SAMRAI_MPI::getSAMRAIWorld()),
   d_ratio(IntVector::getZero(dim)),

   d_local_number_of_cells(0),
   d_global_number_of_cells(0),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(0),

   d_max_number_of_boxes(0),
   d_min_number_of_boxes(0),
   d_max_number_of_cells(0),
   d_min_number_of_cells(0),

   d_local_max_box_size(),
   d_global_max_box_size(),
   d_local_min_box_size(),
   d_global_min_box_size(),

   d_local_bounding_box(),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(nullptr),
   d_persistent_overlap_connectors(nullptr),
   d_handle(),
   d_grid_geometry(),
   d_locked(false)
{
   getFromRestart(restart_db, grid_geom);
}

BoxLevel::BoxLevel(
   const BoxLevel& rhs):
   d_mpi(rhs.d_mpi),
   d_boxes(rhs.d_boxes),
   d_global_boxes(rhs.d_global_boxes),
   d_ratio(rhs.d_ratio),

   d_local_number_of_cells(rhs.d_local_number_of_cells),
   d_global_number_of_cells(rhs.d_global_number_of_cells),
   d_local_number_of_boxes(rhs.d_local_number_of_boxes),
   d_global_number_of_boxes(rhs.d_global_number_of_boxes),

   d_max_number_of_boxes(rhs.d_max_number_of_boxes),
   d_min_number_of_boxes(rhs.d_min_number_of_boxes),
   d_max_number_of_cells(rhs.d_max_number_of_cells),
   d_min_number_of_cells(rhs.d_min_number_of_cells),

   d_local_max_box_size(rhs.d_local_max_box_size),
   d_global_max_box_size(rhs.d_global_max_box_size),
   d_local_min_box_size(rhs.d_local_min_box_size),
   d_global_min_box_size(rhs.d_global_min_box_size),

   d_local_bounding_box(rhs.d_local_bounding_box),
   d_local_bounding_box_up_to_date(rhs.d_local_bounding_box_up_to_date),
   d_global_bounding_box(rhs.d_global_bounding_box),
   d_global_data_up_to_date(rhs.d_global_data_up_to_date),

   d_parallel_state(rhs.d_parallel_state),
   d_globalized_version(nullptr),
   d_persistent_overlap_connectors(nullptr),
   d_handle(),
   d_grid_geometry(rhs.d_grid_geometry),
   d_locked(false)
{
   // This cannot be the first constructor call, so no need to set timers.
}

BoxLevel::BoxLevel(
   const IntVector& ratio,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state):
   d_mpi(MPI_COMM_NULL),
   d_ratio(ratio),

   d_local_number_of_cells(0),
   d_global_number_of_cells(0),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(0),

   d_max_number_of_boxes(0),
   d_min_number_of_boxes(0),
   d_max_number_of_cells(0),
   d_min_number_of_cells(0),

   d_local_max_box_size(),
   d_global_max_box_size(),
   d_local_min_box_size(),
   d_global_min_box_size(),

   d_local_bounding_box(),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(nullptr),
   d_persistent_overlap_connectors(nullptr),
   d_handle(),
   d_grid_geometry(),
   d_locked(false)
{
   initialize(BoxContainer(), ratio, grid_geom, mpi, parallel_state);
}

BoxLevel::BoxLevel(
   const BoxContainer& boxes,
   const IntVector& ratio,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state):
   d_mpi(MPI_COMM_NULL),
   d_ratio(ratio),

   d_local_number_of_cells(0),
   d_global_number_of_cells(0),
   d_local_number_of_boxes(0),
   d_global_number_of_boxes(0),

   d_max_number_of_boxes(0),
   d_min_number_of_boxes(0),
   d_max_number_of_cells(0),
   d_min_number_of_cells(0),

   d_local_max_box_size(),
   d_global_max_box_size(),
   d_local_min_box_size(),
   d_global_min_box_size(),

   d_local_bounding_box(),
   d_local_bounding_box_up_to_date(false),
   d_global_bounding_box(),
   d_global_data_up_to_date(false),

   d_parallel_state(DISTRIBUTED),
   d_globalized_version(nullptr),
   d_persistent_overlap_connectors(nullptr),
   d_handle(),
   d_grid_geometry(),
   d_locked(false)
{
   initialize(boxes, ratio, grid_geom, mpi, parallel_state);
}

BoxLevel::~BoxLevel()
{
   d_locked = false;
   clear();
   d_persistent_overlap_connectors.reset();
}

/*
 ***********************************************************************
 * Assignment operator.
 ***********************************************************************
 */
BoxLevel&
BoxLevel::operator = (
   const BoxLevel& rhs)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::operator =: operating on locked BoxLevel."
         << std::endl);
   }
   if (&rhs != this) {
      /*
       * Protect this block from assignment to self because it is
       * inefficient and it removes d_boxes data before resetting it.
       */

      deallocateGlobalizedVersion();
      clearPersistentOverlapConnectors();
      detachMyHandle();

      d_parallel_state = rhs.d_parallel_state;
      d_mpi = rhs.d_mpi;
      d_ratio = rhs.d_ratio;

      d_local_number_of_cells = rhs.d_local_number_of_cells;
      d_local_number_of_boxes = rhs.d_local_number_of_boxes;
      d_global_number_of_cells = rhs.d_global_number_of_cells;
      d_global_number_of_boxes = rhs.d_global_number_of_boxes;

      d_local_max_box_size = rhs.d_local_max_box_size;
      d_global_max_box_size = rhs.d_global_max_box_size;
      d_local_min_box_size = rhs.d_local_min_box_size;
      d_global_min_box_size = rhs.d_global_min_box_size;

      d_local_bounding_box = rhs.d_local_bounding_box;
      d_local_bounding_box_up_to_date = rhs.d_local_bounding_box_up_to_date;
      d_global_bounding_box = rhs.d_global_bounding_box;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;

      d_boxes = rhs.d_boxes;
      d_global_boxes = rhs.d_global_boxes;
      d_grid_geometry = rhs.d_grid_geometry;
   }
   return *this;
}

void
BoxLevel::initialize(
   const BoxContainer& boxes,
   const IntVector& ratio,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::initialize(): operating on locked BoxLevel."
         << std::endl);
   }

   d_boxes = boxes;
   d_boxes.order();
   initializePrivate(
      ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void
BoxLevel::swapInitialize(
   BoxContainer& boxes,
   const IntVector& ratio,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::swapInitialize(): operating on locked BoxLevel."
         << std::endl);
   }
   TBOX_ASSERT(&boxes != &d_boxes);   // Library error if this fails.
   d_boxes.swap(boxes);
   d_boxes.order();
   initializePrivate(ratio,
      grid_geom,
      mpi,
      parallel_state);
}

void
BoxLevel::finalize()
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::finalize(): operating on locked BoxLevel."
         << std::endl);
   }

   // Erase non-local Boxes, if any, from d_boxes.
   for (BoxContainer::iterator mbi = d_boxes.begin();
        mbi != d_boxes.end(); /* incremented in loop */) {
      if (mbi->getOwnerRank() != d_mpi.getRank()) {
         d_boxes.erase(mbi++);
      } else {
         ++mbi;
      }
   }

   computeLocalRedundantData();
}

void
BoxLevel::initializePrivate(
   const IntVector& ratio,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom,
   const tbox::SAMRAI_MPI& mpi,
   const ParallelState parallel_state)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio);
   t_initialize_private->start();

   d_ratio = ratio;
   if (d_ratio.getNumBlocks() != grid_geom->getNumberBlocks()) {
      if (d_ratio.max() == d_ratio.min()) {
         size_t new_size = grid_geom->getNumberBlocks();
         d_ratio = IntVector(d_ratio, new_size);
      } else {
         TBOX_ERROR("BoxLevel::initializePrivate: anisotropic refinement\n"
            << "ratio " << ratio << " must be \n"
            << "defined for " << grid_geom->getNumberBlocks() << " blocks."
            << std::endl);
      }
   }

   clearForBoxChanges();

   d_mpi = mpi;
   d_grid_geometry = grid_geom;

   if (parallel_state == DISTRIBUTED) {
      d_global_boxes.clear();
   } else {
      d_global_boxes = d_boxes;
   }

   // Erase non-local Boxes, if any, from d_boxes.
   for (BoxContainer::iterator mbi = d_boxes.begin();
        mbi != d_boxes.end(); /* incremented in loop */) {
      if (mbi->getOwnerRank() != d_mpi.getRank()) {
         d_boxes.erase(mbi++);
      } else {
         ++mbi;
      }
   }

   d_parallel_state = parallel_state;
   d_global_number_of_cells = 0;
   d_global_number_of_boxes = 0;
   d_max_number_of_boxes = 0;
   d_min_number_of_boxes = 0;
   d_max_number_of_cells = 0;
   d_min_number_of_cells = 0;
   d_local_bounding_box_up_to_date = false;
   d_global_data_up_to_date = false;
   computeLocalRedundantData();

   t_initialize_private->stop();
}

bool
BoxLevel::operator == (
   const BoxLevel& r) const
{
   if (this == &r) {
      return true;
   }

   if (d_ratio != r.d_ratio) {
      return false;
   }

   if (d_mpi != r.d_mpi) {
      return false;
   }

   if (getBoxes() != r.getBoxes()) {
      return false;
   }

   return true;
}

bool
BoxLevel::operator != (
   const BoxLevel& r) const
{
   if (this == &r) {
      return false;
   }

   if (d_ratio != r.d_ratio) {
      return true;
   }

   if (d_mpi != r.d_mpi) {
      return true;
   }

   if (getBoxes() != r.getBoxes()) {
      return true;
   }

   return false;
}

/*
 ***********************************************************************
 * Clear data and reset them to unusuable values.
 *
 * Note: don't use IntVector::getOne here, because SAMRAI may have
 * already shut down.
 ***********************************************************************
 */
void
BoxLevel::removePeriodicImageBoxes()
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::removePeriodicImageBoxes(): operating on locked BoxLevel."
         << std::endl);
   }
   if (isInitialized()) {
      clearForBoxChanges();
      d_boxes.removePeriodicImageBoxes();
      if (d_parallel_state == GLOBALIZED) {
         d_global_boxes.removePeriodicImageBoxes();
      }
   }
}

/*
 ***********************************************************************
 * Clear data and reset them to unusuable values.
 *
 * Note: don't use IntVector::getOne here, because SAMRAI may have
 * already shut down.
 ***********************************************************************
 */
void
BoxLevel::clear()
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::clear(): operating on locked BoxLevel."
         << std::endl);
   }
   if (isInitialized()) {
      clearForBoxChanges();
      d_mpi = tbox::SAMRAI_MPI(MPI_COMM_NULL);
      d_boxes.clear();
      d_global_boxes.clear();
      d_ratio(0,0) = 0;
      d_local_number_of_cells = 0;
      d_global_number_of_cells = 0;
      d_local_number_of_boxes = 0;
      d_global_number_of_boxes = 0;
      d_max_number_of_boxes = 0;
      d_min_number_of_boxes = 0;
      d_max_number_of_cells = 0;
      d_min_number_of_cells = 0;
      d_local_bounding_box.clear();
      d_local_bounding_box_up_to_date = false;
      d_global_bounding_box.clear();
      d_global_data_up_to_date = false;
      d_local_max_box_size.clear();
      d_local_min_box_size.clear();
      d_global_max_box_size.clear();
      d_global_min_box_size.clear();
      d_parallel_state = DISTRIBUTED;
      d_grid_geometry.reset();
   }
}

void
BoxLevel::swap(
   BoxLevel& level_a,
   BoxLevel& level_b)
{
   if (level_a.locked() || level_b.locked()) {
      TBOX_ERROR("BoxLevel::initialize(): operating on locked BoxLevel."
         << std::endl);
   }

   if (&level_a != &level_b) {
      if (level_a.isInitialized() && level_b.isInitialized()) {
         TBOX_ASSERT_OBJDIM_EQUALITY2(level_a, level_b);
      }

      level_a.clearPersistentOverlapConnectors();
      level_b.clearPersistentOverlapConnectors();

      level_a.detachMyHandle();
      level_b.detachMyHandle();

      // Swap objects supporting swap operation.
      level_a.d_boxes.swap(level_b.d_boxes);
      level_a.d_global_boxes.swap(level_b.d_global_boxes);
      level_a.d_local_bounding_box.swap(level_b.d_local_bounding_box);
      level_a.d_local_min_box_size.swap(level_b.d_local_min_box_size);
      level_a.d_local_max_box_size.swap(level_b.d_local_max_box_size);
      level_a.d_global_bounding_box.swap(level_b.d_global_bounding_box);

      // Swap objects not supporting swap operation.

      int tmpint;
      bool tmpbool;
      Box tmpbox(level_a.getDim());
      ParallelState tmpstate;
      std::unique_ptr<const BoxLevel> tmpmbl;
      tbox::SAMRAI_MPI tmpmpi(MPI_COMM_NULL);
      std::shared_ptr<const BaseGridGeometry> tmpgridgeom(
         level_a.getGridGeometry());

      tmpstate = level_a.d_parallel_state;
      level_a.d_parallel_state = level_b.d_parallel_state;
      level_b.d_parallel_state = tmpstate;

      tmpmpi = level_a.d_mpi;
      level_a.d_mpi = level_b.d_mpi;
      level_b.d_mpi = tmpmpi;

      IntVector tmpvec = level_a.d_ratio;
      level_a.d_ratio = level_b.d_ratio;
      level_b.d_ratio = tmpvec;

      tmpint = static_cast<int>(level_a.d_local_number_of_cells);
      level_a.d_local_number_of_cells = level_b.d_local_number_of_cells;
      level_b.d_local_number_of_cells = tmpint;

      long int tmplongint = level_a.d_global_number_of_cells;
      level_a.d_global_number_of_cells = level_b.d_global_number_of_cells;
      level_b.d_global_number_of_cells = tmplongint;

      tmpint = static_cast<int>(level_a.d_local_number_of_boxes);
      level_a.d_local_number_of_boxes = level_b.d_local_number_of_boxes;
      level_b.d_local_number_of_boxes = tmpint;

      tmpint = level_a.d_global_number_of_boxes;
      level_a.d_global_number_of_boxes = level_b.d_global_number_of_boxes;
      level_b.d_global_number_of_boxes = tmpint;

      tmpbool = level_a.d_local_bounding_box_up_to_date;
      level_a.d_local_bounding_box_up_to_date = level_b.d_local_bounding_box_up_to_date;
      level_b.d_local_bounding_box_up_to_date = tmpbool;

      tmpbool = level_a.d_global_data_up_to_date;
      level_a.d_global_data_up_to_date = level_b.d_global_data_up_to_date;
      level_b.d_global_data_up_to_date = tmpbool;

      tmpmbl = std::move(level_a.d_globalized_version);
      level_a.d_globalized_version = std::move(level_b.d_globalized_version);
      level_b.d_globalized_version = std::move(tmpmbl);

      level_a.d_grid_geometry = level_b.d_grid_geometry;
      level_b.d_grid_geometry = tmpgridgeom;
   }
}

void
BoxLevel::computeLocalRedundantData()
{
   const IntVector max_vec(d_ratio.getDim(), tbox::MathUtilities<int>::getMax());
   const IntVector& zero_vec = IntVector::getZero(d_ratio.getDim());
   const size_t nblocks = d_grid_geometry->getNumberBlocks();

   d_local_number_of_boxes = 0;
   d_local_number_of_cells = 0;

   d_local_bounding_box.clear();
   d_local_min_box_size.clear();
   d_local_max_box_size.clear();
   d_local_bounding_box.resize(nblocks, Box(d_grid_geometry->getDim()));
   d_local_min_box_size.resize(nblocks, max_vec);
   d_local_max_box_size.resize(nblocks, zero_vec);

   for (RealBoxConstIterator ni(d_boxes.realBegin());
        ni != d_boxes.realEnd(); ++ni) {

      const BlockId::block_t& block_num = ni->getBlockId().getBlockValue();
      const IntVector boxdim(ni->numberCells());
      ++d_local_number_of_boxes;
      d_local_number_of_cells += boxdim.getProduct();
      d_local_bounding_box[block_num] += *ni;
      d_local_min_box_size[block_num].min(boxdim);
      d_local_max_box_size[block_num].max(boxdim);

   }

   d_local_bounding_box_up_to_date = true;
   d_global_data_up_to_date = false;
}

/*
 ****************************************************************************
 * Perform global reductions to get characteristics of the global data
 * without globalizing.  Data that can be reduced are combined into
 * arrays for reduction.  We do one sum reduction and one max reduction.
 * All the data we need fall into one or the other so it's all we need to
 * do.
 ****************************************************************************
 */
void
BoxLevel::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isInitialized());

   if (d_global_data_up_to_date) {
      return;
   }

   t_cache_global_reduced_data->barrierAndStart();

   const size_t nblocks = d_grid_geometry->getNumberBlocks();

   /*
    * Sum reduction is used to compute the global sums of box count
    * and cell count.
    */
   if (d_parallel_state == GLOBALIZED) {
      d_global_number_of_boxes = 0;
      d_global_number_of_cells = 0;
      for (RealBoxConstIterator ni(d_global_boxes.realBegin());
           ni != d_global_boxes.realEnd();
           ++ni) {
         ++d_global_number_of_boxes;
         d_global_number_of_cells += ni->size();
      }
   } else {
      if (d_mpi.getSize() > 1) {
         unsigned long int tmpa[2], tmpb[2];
         tmpa[0] = getLocalNumberOfBoxes();
         tmpa[1] = getLocalNumberOfCells();
         d_mpi.Allreduce(tmpa,
            tmpb,                        // Better to use MPI_IN_PLACE, but not some MPI's do not support.
            2,
            MPI_LONG,
            MPI_SUM);
         d_global_number_of_boxes = static_cast<int>(tmpb[0]);
         d_global_number_of_cells = static_cast<size_t>(tmpb[1]);
      } else {
         d_global_number_of_boxes = getLocalNumberOfBoxes();
         d_global_number_of_cells = getLocalNumberOfCells();
      }
   }

   if (d_global_bounding_box.size() != nblocks) {
      d_global_bounding_box.resize(nblocks, Box(getDim()));
      d_global_min_box_size.resize(nblocks, IntVector(getDim()));
      d_global_max_box_size.resize(nblocks, IntVector(getDim()));
   }

   /*
    * Max reduction is used to compute max/min box counts, max/min
    * cell counts, max/min box sizes, and bounding boxes.
    */
   if (d_mpi.getSize() == 1) {

      d_global_bounding_box = d_local_bounding_box;
      d_max_number_of_boxes = d_min_number_of_boxes =
            static_cast<int>(getLocalNumberOfBoxes());
      d_max_number_of_cells = d_min_number_of_cells =
            static_cast<int>(getLocalNumberOfCells());
      d_global_max_box_size = d_local_max_box_size;
      d_global_min_box_size = d_local_min_box_size;

   } else {

      if (d_mpi.getSize() > 1) {
         const tbox::Dimension& dim(getDim());

         std::vector<int> send_mesg;
         send_mesg.reserve(nblocks * 4 * dim.getValue() + 4);
         for (BlockId::block_t bn = 0; bn < nblocks; ++bn) {
            for (int i = 0; i < dim.getValue(); ++i) {
               send_mesg.push_back(-d_local_bounding_box[bn].lower()[i]);
               send_mesg.push_back(d_local_bounding_box[bn].upper()[i]);
               send_mesg.push_back(-d_local_min_box_size[bn][i]);
               send_mesg.push_back(d_local_max_box_size[bn][i]);
            }
         }
         send_mesg.push_back(static_cast<int>(getLocalNumberOfBoxes()));
         send_mesg.push_back(-static_cast<int>(getLocalNumberOfBoxes()));
         send_mesg.push_back(static_cast<int>(getLocalNumberOfCells()));
         send_mesg.push_back(-static_cast<int>(getLocalNumberOfCells()));

         std::vector<int> recv_mesg(send_mesg.size());
         d_mpi.Allreduce(
            &send_mesg[0],
            &recv_mesg[0],
            static_cast<int>(send_mesg.size()),
            MPI_INT,
            MPI_MAX);

         int tmpi = -1;
         for (BlockId::block_t bn = 0; bn < nblocks; ++bn) {
            for (int i = 0; i < dim.getValue(); ++i) {
               d_global_bounding_box[bn].setLower(static_cast<Box::dir_t>(i),
                  -recv_mesg[++tmpi]);
               d_global_bounding_box[bn].setUpper(static_cast<Box::dir_t>(i),
                  recv_mesg[++tmpi]);
               d_global_min_box_size[bn][i] = -recv_mesg[++tmpi];
               d_global_max_box_size[bn][i] = recv_mesg[++tmpi];
            }
            d_global_bounding_box[bn].setBlockId(BlockId(bn));
         }
         d_max_number_of_boxes = recv_mesg[++tmpi];
         d_min_number_of_boxes = -recv_mesg[++tmpi];
         d_max_number_of_cells = recv_mesg[++tmpi];
         d_min_number_of_cells = -recv_mesg[++tmpi];
         TBOX_ASSERT(tmpi == int(recv_mesg.size() - 1));

      } else {
         d_global_bounding_box = d_local_bounding_box;
      }
   }

   d_global_data_up_to_date = true;

   t_cache_global_reduced_data->stop();
}

int
BoxLevel::getLocalNumberOfBoxes(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_mpi.getRank()) {
      TBOX_ERROR(
         "Non-local boxes are not available in DISTRIBUTED mode."
         << std::endl);
   }
#endif
   TBOX_ASSERT(rank >= 0 && rank < d_mpi.getSize());

   if (rank == d_mpi.getRank()) {
      return d_local_number_of_boxes;
   } else {
      int count = 0;
      BoxContainerSingleOwnerIterator mbi(d_global_boxes.begin(rank));
      for ( ; mbi != d_global_boxes.end(rank); ++mbi) {
         if (!(*mbi).isPeriodicImage()) {
            ++count;
         }
      }

      return count;
   }
}

size_t
BoxLevel::getLocalNumberOfCells(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_mpi.getRank()) {
      TBOX_ERROR(
         "Non-local boxes are not available in DISTRIBUTED mode."
         << std::endl);
   }
#endif
   TBOX_ASSERT(rank >= 0 && rank < d_mpi.getSize());

   if (rank == d_mpi.getRank()) {
      return d_local_number_of_cells;
   } else {
      size_t count = 0;
      BoxContainerSingleOwnerIterator mbi(d_global_boxes.begin(rank));
      for ( ; mbi != d_global_boxes.end(rank); ++mbi) {
         if (!(*mbi).isPeriodicImage()) {
            count += (*mbi).size();
         }
      }

      return count;
   }
}

bool
BoxLevel::getSpatiallyEqualBox(
   const Box& box_to_match,
   const BlockId& block_id,
   Box& matching_box) const
{
   bool box_exists = false;
   for (BoxContainerSingleBlockIterator itr(d_boxes.begin(block_id));
        itr != d_boxes.end(block_id); ++itr) {
      if (box_to_match.isSpatiallyEqual(*itr)) {
         box_exists = true;
         matching_box = *itr;
         break;
      }
   }
   return box_exists;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::setParallelState(
   const ParallelState parallel_state)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::setParallelState(): operating on locked BoxLevel."
         << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isInitialized()) {
      TBOX_ERROR(
         "BoxLevel::setParallelState: Cannot change the parallel state of\n"
         << "an uninitialized BoxLevel.  See BoxLevel::initialize()"
         << std::endl);
   }
#endif
   if (d_parallel_state == DISTRIBUTED && parallel_state == GLOBALIZED) {
      acquireRemoteBoxes();
   } else if (d_parallel_state == GLOBALIZED &&
              parallel_state == DISTRIBUTED) {
      d_global_boxes.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::acquireRemoteBoxes()
{
   BoxLevel* object = this;
   acquireRemoteBoxes(1, &object);
}

/*
 ***********************************************************************
 * Acquire remote Boxes for multiple BoxLevels.
 * This method combines communication for the multiple
 * BoxLevels to increase message passing efficiency.
 *
 * Note: This method is stateless (could be static).
 ***********************************************************************
 */

void
BoxLevel::acquireRemoteBoxes(
   const int num_sets,
   BoxLevel* multiple_box_levels[])
{
   if (d_mpi.getSize() == 1) {
      // In single-proc mode, we already have all the Boxes already.
      for (int n = 0; n < num_sets; ++n) {
         multiple_box_levels[n]->d_global_boxes =
            multiple_box_levels[n]->d_boxes;
      }
      return;
   }

   t_acquire_remote_boxes->start();
   int n;

#ifdef DEBUG_CHECK_ASSERTIONS
   for (n = 0; n < num_sets; ++n) {
      if (multiple_box_levels[n]->getParallelState() !=
          DISTRIBUTED) {
         TBOX_ERROR("BoxLevel objects must be in distributed mode\n"
            << "when acquiring remote boxes.\n");
      }
   }
#endif

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack Boxes from all BoxLevels into a single message.
    */
   for (n = 0; n < num_sets; ++n) {
      const BoxLevel& box_level =
         *multiple_box_levels[n];
      box_level.acquireRemoteBoxes_pack(send_mesg);
   }
   int send_mesg_size = static_cast<int>(send_mesg.size());

   /*
    * Send and receive the data.
    */

   std::vector<int> recv_mesg_size(d_mpi.getSize());
   d_mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(d_mpi.getSize());
   int totl_size = 0;
   for (n = 0; n < d_mpi.getSize(); ++n) {
      proc_offset[n] = totl_size;
      totl_size += recv_mesg_size[n];
   }
   recv_mesg.resize(totl_size, BAD_INT);
   d_mpi.Allgatherv(&send_mesg[0],
      send_mesg_size,
      MPI_INT,
      &recv_mesg[0],
      &recv_mesg_size[0],
      &proc_offset[0],
      MPI_INT);

   /*
    * Extract Box info received from other processors.
    */
   for (n = 0; n < num_sets; ++n) {
      BoxLevel& box_level =
         *multiple_box_levels[n];
      box_level.acquireRemoteBoxes_unpack(recv_mesg,
         proc_offset);
   }

   t_acquire_remote_boxes->stop();

}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::acquireRemoteBoxes_pack(
   std::vector<int>& send_mesg) const
{
   const tbox::Dimension& dim(getDim());
   /*
    * Box acquisition occurs during globalization.  Thus, do not
    * rely on current value of d_parallel_state.
    */

   /*
    * Pack Box info from d_boxes into send_mesg,
    * starting at the offset location.
    */
   /*
    * Information to be packed:
    *   - Number of Boxes from self
    *   - Self Boxes
    */
   const int box_com_buf_size = Box::commBufferSize(dim);
   const int send_mesg_size = 1 + box_com_buf_size
      * static_cast<int>(d_boxes.size());
   const int old_size = static_cast<int>(send_mesg.size());
   send_mesg.resize(old_size + send_mesg_size, BAD_INT);

   int* ptr = &send_mesg[0] + old_size;
   *(ptr++) = static_cast<int>(d_boxes.size());

   for (BoxContainer::const_iterator i_boxes = d_boxes.begin();
        i_boxes != d_boxes.end();
        ++i_boxes) {
      (*i_boxes).putToIntBuffer(ptr);
      ptr += box_com_buf_size;
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::acquireRemoteBoxes_unpack(
   const std::vector<int>& recv_mesg,
   std::vector<int>& proc_offset)
{
   const tbox::Dimension& dim(getDim());
   /*
    * Unpack Box info from recv_mesg into d_global_boxes,
    * starting at the offset location.
    * Advance the proc_offset past the used data.
    */
   int n;
   int box_com_buf_size = Box::commBufferSize(dim);

   for (n = 0; n < d_mpi.getSize(); ++n) {
      if (n != d_mpi.getRank()) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int n_self_boxes = *(ptr++);
         proc_offset[d_mpi.getRank()] += (n_self_boxes) * box_com_buf_size;

         int i;
         Box box(dim);

         for (i = 0; i < n_self_boxes; ++i) {
            box.getFromIntBuffer(ptr);
            d_global_boxes.insert(
               d_global_boxes.end(), box);
            ptr += box_com_buf_size;
         }

      } else {
         d_global_boxes.insert(
            d_boxes.begin(), d_boxes.end());
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

BoxContainer::const_iterator
BoxLevel::addBox(
   const Box& box,
   const BlockId& block_id)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::addBox(): operating on locked BoxLevel."
         << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually adding Boxes is a local process\n"
         << "so it can only be performed in\n"
         << "DISTRIBUTED state." << std::endl);
   }
   if (box.getBlockId() != BlockId::invalidId()) {
      TBOX_ASSERT(box.getBlockId() == block_id);
   }
#endif

   clearForBoxChanges(false);

   const PeriodicShiftCatalog& shift_catalog =
      getGridGeometry()->getPeriodicShiftCatalog(); 

   BoxContainer::iterator new_iterator = d_boxes.begin();

   if (d_boxes.empty()) {
      Box new_box(
         box,
         LocalId::getZero(),
         d_mpi.getRank(),
         shift_catalog.getZeroShiftNumber());
      new_box.setBlockId(block_id);
      new_iterator = d_boxes.insert(d_boxes.end(), new_box);
   } else {
      // Set new_index to one more than the largest index used.
      BoxContainer::iterator ni = d_boxes.end();
      do {
         TBOX_ASSERT(ni != d_boxes.begin());   // There should not be all periodic images.
         --ni;
      } while (ni->isPeriodicImage());
      LocalId new_index = ni->getLocalId() + 1;
      Box new_box(
         box, new_index, d_mpi.getRank());
      new_box.setBlockId(block_id);
      new_iterator = d_boxes.insert(ni, new_box);
   }

   const IntVector box_size(new_iterator->numberCells());
   ++d_local_number_of_boxes;
   d_local_number_of_cells += box.size();
   d_local_bounding_box[block_id.getBlockValue()] += *new_iterator;
   d_local_max_box_size[block_id.getBlockValue()].max(box_size);
   d_local_min_box_size[block_id.getBlockValue()].min(box_size);
   d_global_data_up_to_date = false;

   return new_iterator;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addPeriodicBox(
   const Box& ref_box,
   const PeriodicId& shift_number)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::addPeriodicBox(): operating on locked BoxLevel."
         << std::endl);
   }

   const PeriodicShiftCatalog& shift_catalog =
      getGridGeometry()->getPeriodicShiftCatalog();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (shift_number == shift_catalog.getZeroShiftNumber()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox cannot be used to add regular box."
         << std::endl);
   }
#endif

   clearForBoxChanges(false);

   Box image_box(ref_box, shift_number, d_ratio, shift_catalog);

#ifdef DEBUG_CHECK_ASSERTIONS
   BoxContainer& boxes =
      d_parallel_state == DISTRIBUTED ? d_boxes : d_global_boxes;
   /*
    * Sanity checks:
    *
    * - Require that the real version of the reference Box exists
    *   before adding the periodic image Box.
    */
   Box real_box(getDim(),
                ref_box.getGlobalId(),
                shift_catalog.getZeroShiftNumber()); 
   if (boxes.find(real_box) == boxes.end()) {
      TBOX_ERROR(
         "BoxLevel::addPeriodicBox: cannot add periodic image Box "
         << image_box
         << "\nwithout the real Box (" << real_box
         << ") already in the BoxLevel.\n");
   }
#endif

   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(image_box);
   }
   if (image_box.getOwnerRank() == d_mpi.getRank()) {
      d_boxes.insert(image_box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::addBox(
   const Box& box)
{
   TBOX_ASSERT(box.getLocalId().isValid());
   if (locked()) {
      TBOX_ERROR("BoxLevel::addBox(): operating on locked BoxLevel."
         << std::endl);
   }
   clearForBoxChanges(false);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Sanity checks:
    * - Require that the real Box exists before adding the periodic image Box.
    */
   if (box.isPeriodicImage()) {
      Box real_box(getDim(),
                   box.getGlobalId(),
                   getGridGeometry()->getPeriodicShiftCatalog().
                      getZeroShiftNumber());
      BoxContainer& boxes = box.getOwnerRank() ==
         d_mpi.getRank() ? d_boxes : d_global_boxes;
      if (boxes.find(real_box) == boxes.end()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add periodic image Box "
            << box
            << "\nwithout the real Box (" << real_box
            << ") already in the BoxLevel.\n");
      }
      if (d_global_boxes.find(box) !=
          d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxLevel::addBox: cannot add Box "
            << box
            << "\nbecause it already exists ("
            << *boxes.find(box) << "\n");
      }
   }
#endif

   // Update counters.
   if (!box.isPeriodicImage()) {
      if (box.getOwnerRank() == d_mpi.getRank()) {
         const IntVector box_size(box.numberCells());
         ++d_local_number_of_boxes;
         d_local_number_of_cells += box.size();
         d_local_bounding_box[box.getBlockId().getBlockValue()] += box;
         d_local_max_box_size[box.getBlockId().getBlockValue()].max(box_size);
         d_local_min_box_size[box.getBlockId().getBlockValue()].min(box_size);
         d_global_data_up_to_date = false;
      }
      d_global_data_up_to_date = false;
      /*
       * TODO: bug: if some procs add a real Box and others do not,
       * their d_global_data_up_to_date flags will be inconsistent
       * resulting in incomplete participation in future collective
       * communication to compute that parameter.
       */
   }

   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(box);
   }
   if (box.getOwnerRank() == d_mpi.getRank()) {
      d_boxes.insert(box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   BoxContainer::iterator& ibox)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::eraseBox(): operating on locked BoxLevel."
         << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually erasing boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state." << std::endl);
   }
#endif

   clearForBoxChanges();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (ibox != d_boxes.find(*ibox)) {
      TBOX_ERROR("BoxLevel::eraseBox: Attempt to erase a\n"
         << "Box that does not belong to the BoxLevel\n"
         << "object.\n");
   }
#endif

   if (ibox->isPeriodicImage()) {
      d_boxes.erase(ibox++);
      // No need to update counters (they neglect periodic images).
   } else {
      /*
       * Update counters.  Bounding box cannot be updated (without
       * recomputing) because we don't know how the erased Box
       * affects the bounding box.
       */
      d_local_bounding_box_up_to_date = d_global_data_up_to_date = false;
      --d_local_number_of_boxes;
      d_local_number_of_cells -= ibox->size();
      // Erase real Box and its periodic images.
      const LocalId& local_id = ibox->getLocalId();
      do {
         d_boxes.erase(ibox++);
      } while (ibox != d_boxes.end() && ibox->getLocalId() ==
               local_id);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
BoxLevel::eraseBox(
   const Box& box)
{
   if (locked()) {
      TBOX_ERROR("BoxLevel::eraseBox(): operating on locked BoxLevel."
         << std::endl);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state != DISTRIBUTED) {
      TBOX_ERROR("Individually erasing Boxes is a local process\n"
         << "so it can only be performed in\n"
         << "distributed state." << std::endl);
   }
#endif

   clearForBoxChanges();

   d_local_bounding_box_up_to_date = d_global_data_up_to_date = false;

   BoxContainer::iterator ibox = d_boxes.find(box);
   if (ibox == d_boxes.end()) {
      TBOX_ERROR("BoxLevel::eraseBox: Box to be erased ("
         << box << ") is NOT a part of the BoxLevel.\n");
   }
   d_boxes.erase(ibox);
}

/*
 ****************************************************************************
 ****************************************************************************
 */
const BoxLevel&
BoxLevel::getGlobalizedVersion() const
{
   TBOX_ASSERT(isInitialized());

   if (d_parallel_state == GLOBALIZED) {
      return *this;
   }

   if (!d_globalized_version) {
      auto globalized_version = std::make_unique<BoxLevel>(*this);
      globalized_version->setParallelState(GLOBALIZED);
      TBOX_ASSERT(globalized_version->getParallelState() == GLOBALIZED);
      d_globalized_version = std::move(globalized_version);
   }

   TBOX_ASSERT(d_globalized_version->getParallelState() == GLOBALIZED);
   return *d_globalized_version;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
PersistentOverlapConnectors&
BoxLevel::getPersistentOverlapConnectors() const
{
   if (!d_persistent_overlap_connectors) {
      d_persistent_overlap_connectors.reset(new PersistentOverlapConnectors(*this));
   }
   return *d_persistent_overlap_connectors;
}

LocalId
BoxLevel::getFirstLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxContainer& boxes = getBoxes();
   if (boxes.empty()) {
      return s_negative_one_local_id;
   }
   BoxContainer::const_iterator ni = boxes.begin();
   while (ni->isPeriodicImage()) {
      TBOX_ASSERT(ni != boxes.end());   // There should be a real box!
      ++ni;
   }
   return ni->getLocalId();
}

LocalId
BoxLevel::getLastLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxContainer& boxes = getBoxes();
   if (boxes.empty()) {
      return s_negative_one_local_id;
   }
   LocalId last_local_id(0);
   for (BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end(); ++ni) {
      if (last_local_id < ni->getLocalId()) {
         last_local_id = ni->getLocalId();
      }
   }
   return last_local_id;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
BoxLevel::hasBox(
   const Box& box) const
{
   if (box.getOwnerRank() == d_mpi.getRank()) {

      BoxContainer::const_iterator ni = d_boxes.find(box);
      return ni != d_boxes.end();

   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state == DISTRIBUTED) {
         TBOX_ERROR("BoxLevel: Cannot check on remote Box "
            << box << " while in DISTRIBUTED mode.\n"
            << "See BoxLevel::setParallelState()." << std::endl);
      }
#endif
      BoxContainer::const_iterator ni = d_global_boxes.find(box);
      return ni != d_global_boxes.end();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
BoxContainer::const_iterator
BoxLevel::getBoxStrict(
   const Box& box) const
{
   if (box.getOwnerRank() == d_mpi.getRank()) {
      BoxContainer::const_iterator ni = d_boxes.find(box);
      if (ni == d_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level." << std::endl);
      }

      return ni;
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "BoxLevel::getBox: cannot get remote box "
            << box << " without being in globalized state." << std::endl);
      }
#endif
      BoxContainer::const_iterator ni = d_global_boxes.find(box);
      if (ni == d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level." << std::endl);
      }

      return ni;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
BoxContainer::const_iterator
BoxLevel::getBoxStrict(
   const BoxId& box_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (box_id.getOwnerRank() != d_mpi.getRank() && d_parallel_state != GLOBALIZED) {
      TBOX_ERROR(
         "BoxLevel::getBoxStrict: cannot get remote box " << box_id
                                                          <<
         " without being in globalized state." << std::endl);
   }
#endif

   Box box(getDim(), box_id);
   if (box.getOwnerRank() == d_mpi.getRank()) {
      BoxContainer::const_iterator ni = d_boxes.find(box);
      if (ni == d_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level." << std::endl);
      }
      return ni;
   } else {
      BoxContainer::const_iterator ni = d_global_boxes.find(box);
      if (ni == d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level." << std::endl);
      }
      return ni;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxLevel::getGlobalBoxes(BoxContainer& global_boxes) const
{
   for (BoxContainer::const_iterator itr = d_global_boxes.begin();
        itr != d_global_boxes.end(); ++itr) {
      global_boxes.pushBack(*itr);
   }
}

/*
 ***********************************************************************
 * Write the BoxLevel to a restart database.
 *
 * Write only local parts.
 ***********************************************************************
 */

void
BoxLevel::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   // This appears to be used in the RedistributedRestartUtility.
   restart_db->putBool("d_is_mapped_box_level", true);

   restart_db->putInteger(
      "HIER_MAPPED_BOX_LEVEL_VERSION", HIER_BOX_LEVEL_VERSION);
   restart_db->putInteger("d_nproc", d_mpi.getSize());
   restart_db->putInteger("d_rank", d_mpi.getRank());
   restart_db->putInteger("dim", d_ratio.getDim().getValue());
   const size_t nblocks = d_grid_geometry->getNumberBlocks();
   for (BlockId::block_t b =0; b < nblocks; ++b) {
      std::string ratio_name = "d_ratio_" +
         tbox::Utilities::intToString(static_cast<int>(b)); 
      std::vector<int> tmp_ratio(d_ratio.getDim().getValue());
      for (unsigned int d = 0; d < d_ratio.getDim().getValue(); ++d) {
         tmp_ratio[d] = d_ratio(b,d);
      }
      restart_db->putIntegerArray(ratio_name,
         &(tmp_ratio[0]),
         d_ratio.getDim().getValue());
   }
   getBoxes().putToRestart(restart_db->putDatabase("mapped_boxes"));
}

/*
 ***********************************************************************
 * Read the BoxLevel from a database.
 ***********************************************************************
 */

void
BoxLevel::getFromRestart(
   tbox::Database& restart_db,
   const std::shared_ptr<const BaseGridGeometry>& grid_geom)
{
   TBOX_ASSERT(restart_db.isInteger("dim"));
   const tbox::Dimension dim(static_cast<unsigned short>(
                                restart_db.getInteger("dim")));
   TBOX_ASSERT(getDim() == dim);

   const size_t nblocks = grid_geom->getNumberBlocks();

   IntVector ratio(nblocks, dim);
   for (BlockId::block_t b = 0; b < nblocks; ++b) {
      std::string ratio_name = "d_ratio_" +
         tbox::Utilities::intToString(static_cast<int>(b));
      std::vector<int> tmp_ratio(dim.getValue());
      restart_db.getIntegerArray(ratio_name, &tmp_ratio[0], dim.getValue());
      for (int d = 0; d < dim.getValue(); ++d) {
         ratio(b,d) = tmp_ratio[d];
      }
   }
 
#ifdef DEBUG_CHECK_ASSERTIONS
   const int version = restart_db.getInteger("HIER_MAPPED_BOX_LEVEL_VERSION");
   const int nproc = restart_db.getInteger("d_nproc");
   const int rank = restart_db.getInteger("d_rank");
#endif
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));
   TBOX_ASSERT(version <= HIER_BOX_LEVEL_VERSION);

   initialize(BoxContainer(), ratio, grid_geom);

   /*
    * Failing these asserts means that we don't have a compatible
    * restart database for the number of processors or we are reading another
    * processor's data.
    */
   TBOX_ASSERT(nproc == d_mpi.getSize());
   TBOX_ASSERT(rank == d_mpi.getRank());

   d_boxes.getFromRestart(*restart_db.getDatabase("mapped_boxes"));
   computeLocalRedundantData();

}

/*
 ***********************************************************************
 * Outputter copy constructor
 ***********************************************************************
 */

BoxLevel::Outputter::Outputter(
   const BoxLevel::Outputter& other):
   d_level(other.d_level),
   d_border(other.d_border),
   d_detail_depth(other.d_detail_depth),
   d_output_statistics(other.d_output_statistics)
{
}

/*
 ***********************************************************************
 * Construct a BoxLevel Outputter with formatting parameters.
 ***********************************************************************
 */

BoxLevel::Outputter::Outputter(
   const BoxLevel& box_level,
   const std::string& border,
   int detail_depth,
   bool output_statistics):
   d_level(box_level),
   d_border(border),
   d_detail_depth(detail_depth),
   d_output_statistics(output_statistics)
{
}

/*
 ***********************************************************************
 * Print out a BoxLevel according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream&
operator << (
   std::ostream& s,
   const BoxLevel::Outputter& format)
{
   if (format.d_output_statistics) {
      BoxLevelStatistics bls(format.d_level);
      bls.printBoxStats(s, format.d_border);
   } else {
      format.d_level.recursivePrint(s, format.d_border, format.d_detail_depth);
   }
   return s;
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxLevel to a stream.
 ***********************************************************************
 */

BoxLevel::Outputter
BoxLevel::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxLevel statistics to a stream.
 ***********************************************************************
 */

BoxLevel::Outputter
BoxLevel::formatStatistics(
   const std::string& border) const
{
   return Outputter(*this, border, 0, true);
}

/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 * Print out global bounding box only if it has been computed already.
 ***********************************************************************
 */

void
BoxLevel::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   if (detail_depth < 0) return;

   if (!isInitialized()) {
      co << border << "Uninitialized.\n";
      return;
   }
   co // << "Address        : " << (void*)this << '\n'
   << border << "Parallel state : "
   << (getParallelState() == DISTRIBUTED ? "DIST" : "GLOB") << '\n'
   << border << "Ratio          : " << getRefinementRatio() << '\n'
   << border << "Box count      : " << d_local_number_of_boxes << ", "
   << d_global_number_of_boxes << '\n'
   << border << "Cell count     : " << d_local_number_of_cells << ", "
   << d_global_number_of_cells << '\n'
   << border << "Bounding box   : " << getLocalBoundingBox(BlockId(0)) << ", "
   << (d_global_data_up_to_date ? getGlobalBoundingBox(BlockId(0)) : Box(getDim()))
   << '\n'
   << border << "Comm,rank,nproc: " << d_mpi.getCommunicator() << ", " << d_mpi.getRank()
   << ", " << d_mpi.getSize() << '\n'
   ;
   if (detail_depth > 0) {
      co << border << "Boxes:\n";
      if (getParallelState() == GLOBALIZED) {
         /*
          * Print boxes from all ranks.
          */
         for (BoxContainer::const_iterator bi = d_global_boxes.begin();
              bi != d_global_boxes.end();
              ++bi) {
            Box box = *bi;
            co << border << "    "
            << box << "   "
            << box.numberCells() << '\n';
         }
      } else {
         /*
          * Print local boxes only.
          */
         for (BoxContainer::const_iterator bi = d_boxes.begin();
              bi != d_boxes.end();
              ++bi) {
            Box box = *bi;
            co << border << "    "
            << box << "   "
            << box.numberCells() << '\n';
         }
      }
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
