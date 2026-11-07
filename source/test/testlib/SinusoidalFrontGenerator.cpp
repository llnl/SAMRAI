/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontGenerator class implementation
 *
 ************************************************************************/
#include "SinusoidalFrontGenerator.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

using namespace SAMRAI;

/*
 ***********************************************************************
 ***********************************************************************
 */
SinusoidalFrontGenerator::SinusoidalFrontGenerator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<tbox::Database>& database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(),
   d_time_shift(0.0),
   d_amplitude(0.2),
   d_buffer_distance(1, std::vector<double>(dim.getValue(), 0.0))
{
   TBOX_ASSERT(hier::VariableDatabase::getDatabase() != 0);

   std::vector<double> init_disp;
   std::vector<double> velocity;
   std::vector<double> period;

   // Parameters set by database, with defaults.
   std::shared_ptr<tbox::Database> sft_db; // SinusoidalFrontGenerator database.

   if (database) {

      sft_db = database->getDatabaseWithDefault("SinusoidalFrontGenerator", sft_db);

      if (database->isDouble("period")) {
         period = database->getDoubleVector("period");
         for (size_t i = 0; i < period.size(); ++i) {
            TBOX_ASSERT(period[i] > 0.0);
         }
      }
      if (database->isDouble("init_disp")) {
         init_disp = database->getDoubleVector("init_disp");
      }
      if (database->isDouble("velocity")) {
         velocity = database->getDoubleVector("velocity");
      }
      d_amplitude = database->getDoubleWithDefault("amplitude", d_amplitude);
      d_time_shift = database->getDoubleWithDefault("time_shift", d_time_shift);

      /*
       * Input parameters to determine whether to tag by buffering
       * fronts, and by how much.
       */
      const std::string bname("buffer_distance_");
      for (int ln = 0; ; ++ln) {
         const std::string lnstr(tbox::Utilities::intToString(ln));

         const std::string bnameln = bname + lnstr;

         std::vector<double> tmpa;

         if (database->isDouble(bnameln)) {
            tmpa = database->getDoubleVector(bnameln);
            if (static_cast<int>(tmpa.size()) != dim.getValue()) {
               TBOX_ERROR(bnameln << " input parameter must have " << dim << " values");
            }
         }

         if (!tmpa.empty()) {
            d_buffer_distance.resize(ln + 1);
            d_buffer_distance.back().swap(tmpa);
         } else {
            break;
         }

      }
   }

   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      d_init_disp[idim] =
         idim < static_cast<int>(init_disp.size()) ? init_disp[idim] : 0.0;
      d_velocity[idim] =
         idim < static_cast<int>(velocity.size()) ? velocity[idim] : 0.0;
      d_period[idim] =
         idim < static_cast<int>(period.size()) ? period[idim] : 1.0e20;
   }

   t_setup = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::setup");
   t_node_pos = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::node_pos");
   t_distance = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::distance");
   t_uval = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::uval");
   t_tag_cells = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::tag_cells");
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SinusoidalFrontGenerator::~SinusoidalFrontGenerator()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void SinusoidalFrontGenerator::applyGradientDetector(
   const std::shared_ptr<hier::PatchHierarchy>& base_hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
   NULL_USE(error_data_time);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation);
   TBOX_ASSERT(base_hierarchy_);
   std::shared_ptr<hier::PatchLevel> level_(
      base_hierarchy_->getPatchLevel(ln));
   TBOX_ASSERT(level_);

   hier::PatchLevel& level = *level_;

   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi) {
      hier::Patch& patch = **pi;

      std::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
      TBOX_ASSERT(tag_cell_data_);
      TBOX_ASSERT(tag_cell_data_->getTime() == error_data_time);

      // Compute tag data for patch.
      computePatchData(patch,
         0,
         tag_cell_data_.get(),
         patch.getBox());

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void SinusoidalFrontGenerator::setTags(
   bool& exact_tagging,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id)
{
   const std::shared_ptr<hier::PatchLevel>& tag_level(
      hierarchy->getPatchLevel(tag_ln));

   resetHierarchyConfiguration(hierarchy, 0, 1);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));

      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(tag_data_id)));
      TBOX_ASSERT(patch_geom);
      TBOX_ASSERT(tag_data);

      computeFrontsData(
         0 /* distance data */,
         0 /* uval data */,
         tag_data.get(),
         tag_data->getBox(),
         (static_cast<size_t>(tag_ln) < d_buffer_distance.size() ?
          d_buffer_distance[tag_ln] : d_buffer_distance.back()),
         patch_geom->getXLower(),
         patch_geom->getDx());

   }

   exact_tagging = false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void SinusoidalFrontGenerator::setDomain(
   hier::BoxContainer& domain,
   double xlo[],
   double xhi[],
   int autoscale_base_nprocs,
   const tbox::SAMRAI_MPI& mpi)
{
   TBOX_ASSERT(autoscale_base_nprocs <= mpi.getSize());
   TBOX_ASSERT(!domain.empty());

   hier::BoxContainer::const_iterator ii = domain.begin();
   ii->getDim();
   const tbox::Dimension& dim = domain.begin()->getDim();

   tbox::Dimension::dir_t doubling_dir = 1;
   while (autoscale_base_nprocs < mpi.getSize()) {
      for (hier::BoxContainer::iterator bi = domain.begin();
           bi != domain.end(); ++bi) {
         hier::Box& input_box = *bi;
         input_box.setUpper(doubling_dir,
            input_box.upper(doubling_dir) + input_box.numberCells(doubling_dir));
      }
      xhi[doubling_dir] += xhi[doubling_dir] - xlo[doubling_dir];
      doubling_dir = static_cast<tbox::Dimension::dir_t>((doubling_dir + 1) % dim.getValue());
      autoscale_base_nprocs *= 2;
      tbox::plog << "autoscale_base_nprocs = " << autoscale_base_nprocs << std::endl
                 << domain.format("IB: ", 2) << std::endl;
   }

   if (autoscale_base_nprocs != mpi.getSize()) {
      TBOX_ERROR("If autoscale_base_nprocs (" << autoscale_base_nprocs << ") is given,\n"
                                              << "number of processes (" << mpi.getSize()
                                              << ") must be\n"
                                              <<
         "a power-of-2 times the value of autoscale_base_nprocs.");
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */
void SinusoidalFrontGenerator::resetHierarchyConfiguration(
   /*! New hierarchy */ const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   /*! Coarsest level */ const int coarsest_level,
   /*! Finest level */ const int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   TBOX_ASSERT(new_hierarchy->getDim() == d_dim);
   d_hierarchy = new_hierarchy;
   TBOX_ASSERT(d_hierarchy);
}

/*
 * Compute the solution data for a patch.
 */
void SinusoidalFrontGenerator::computePatchData(
   const hier::Patch& patch,
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data,
   const hier::Box& fill_box) const
{
   t_setup->start();

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom.get() != 0);

   const double* xlo = patch_geom->getXLower();
   const double* dx = patch_geom->getDx();

   t_setup->stop();

   if (tag_data) {
      computeFrontsData(0, uval_data, tag_data,
         fill_box,
         (patch.getPatchLevelNumber() < static_cast<int>(d_buffer_distance.size()) ?
          d_buffer_distance[patch.getPatchLevelNumber()] : d_buffer_distance.back()),
         xlo, dx);
   } else {
      // Not computing tag => no tag buffer needed.
      computeFrontsData(0, uval_data, tag_data,
         fill_box,
         std::vector<double>(d_dim.getValue(), 0.0),
         xlo, dx);
   }
}

/*
 * Compute the various data due to the fronts.
 */
void SinusoidalFrontGenerator::computeFrontsData(
   pdat::NodeData<double>* dist_data,
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data,
   const hier::Box& fill_box,
   const std::vector<double>& buffer_distance,
   const double xlo[],
   const double dx[]) const
{
   t_setup->start();

   if (dist_data != 0 && tag_data != 0) {
      TBOX_ASSERT(dist_data->getBox().isSpatiallyEqual(tag_data->getBox()));
   }

   // Compute the buffer in terms of cells.
   hier::IntVector buffer_cells(d_dim);
   for (int i = 0; i < d_dim.getValue(); ++i) {
      buffer_cells(i) = static_cast<int>(0.5 + buffer_distance[i] / dx[i]);
   }

   hier::Box pbox(d_dim);
   double time = 0;
   if (tag_data != 0) {
      pbox = tag_data->getBox();
      time = tag_data->getTime() + d_time_shift;
   } else if (uval_data != 0) {
      pbox = uval_data->getBox();
      time = uval_data->getTime() + d_time_shift;
   } else {
      pbox = dist_data->getBox();
      time = dist_data->getTime() + d_time_shift;
   }

   double wave_number[SAMRAI::MAX_DIM_VAL];
   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      wave_number[idim] = 2 * 3.141592654 / d_period[idim];
   }

   t_setup->stop();

   /*
    * Initialize node x-distances from front.
    */

   t_node_pos->start();
   hier::Box front_box = fill_box;
   front_box.grow(buffer_cells);
   front_box.growUpper(hier::IntVector::getOne(d_dim));
   // Squash front_box to a single plane.
   front_box.setUpper(0, pbox.lower(0));
   front_box.setLower(0, pbox.lower(0));
   pdat::ArrayData<double> front_x(front_box, 1);
   pdat::ArrayData<int>::iterator aiend(front_x.getBox(), false);
   for (pdat::ArrayData<int>::iterator ai(front_x.getBox(), true);
        ai != aiend; ++ai) {

      const hier::Index& index = *ai;
      double y = 0.0, cosy = 0.0, z = 0.0, cosz = 1.0;

      y = xlo[1] + dx[1] * (index(1) - pbox.lower()[1]);
      cosy = cos(wave_number[1] * (y + d_init_disp[1] - d_velocity[1] * time));

      if (d_dim.getValue() > 2) {
         z = xlo[2] + dx[2] * (index(2) - pbox.lower()[2]);
         cosz = cos(wave_number[2] * (z + d_init_disp[2] - d_velocity[2] * time));
      }

      front_x(index, 0) = d_velocity[0] * time + d_init_disp[0]
         + d_amplitude * cosy * cosz;
      // tbox::plog << "index=" << index << "   y=" << y << "   cosy=" << cosy << "   z=" << z << "   cosz=" << cosz << "   front_x = " << front_x(index,0) << std::endl;
   }
   t_node_pos->stop();

   if (tag_data) {
      /*
       * Initialize tags to zero then tag specific cells.
       */
      tag_data->fill(0, fill_box);
      hier::Box buffered_box(fill_box);
      buffered_box.grow(buffer_cells);
      hier::BlockId blk0(0);
      pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(buffered_box));
      for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(buffered_box));
           ci != ciend; ++ci) {

         const pdat::CellIndex& cell_index = *ci;
         const hier::Box cell_box(cell_index, cell_index, blk0);

         double min_distance_to_front = tbox::MathUtilities<double>::getMax();
         double max_distance_to_front = -tbox::MathUtilities<double>::getMax();
         // tbox::plog << "initial distances to front: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;
         pdat::NodeIterator niend(pdat::NodeGeometry::end(cell_box));
         for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(cell_box)); ni != niend; ++ni) {

            const pdat::NodeIndex& node_index = *ni;
            hier::Index front_index = node_index;
            front_index(0) = pbox.lower(0);

            double node_x = xlo[0] + dx[0] * (node_index(0) - pbox.lower() (0));

            double distance_to_front = node_x - front_x(front_index, 0);
            min_distance_to_front = tbox::MathUtilities<double>::Min(
                  min_distance_to_front, distance_to_front);
            max_distance_to_front = tbox::MathUtilities<double>::Max(
                  max_distance_to_front, distance_to_front);
            // tbox::plog << "cell_index = " << cell_index << "   node_index = " << node_index << "   node_x = " << node_x << "   front_index = " << front_index << "   front_x = " << front_x(front_index,0) << "   distance to front(" << node_x << ") = " << distance_to_front << std::endl;

         }
         // tbox::plog << "distances to front: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;

         /*
          * Compute shifts needed to put distances in the range [ -.5*d_period[0], .5*d_period[0] ]
          * This makes the distances relative to the nearest front instead of front 0.
          */
#if 1
         const double cycles_up = min_distance_to_front > 0.5 * d_period[0] ?
            0.0 : static_cast<int>(0.5 - min_distance_to_front / d_period[0]);
         const double cycles_dn = max_distance_to_front < 0.5 * d_period[0] ?
            0.0 : static_cast<int>(0.5 + max_distance_to_front / d_period[0]);
         min_distance_to_front += (cycles_up - cycles_dn) * d_period[0];
         max_distance_to_front += (cycles_up - cycles_dn) * d_period[0];
#else
         // This is more readable than the above, but the short innner loop is too slow.
         while (min_distance_to_front < -0.5 * d_period[0]) {
            min_distance_to_front += d_period[0];
            max_distance_to_front += d_period[0];
         }
         while (max_distance_to_front > 0.5 * d_period[0]) {
            min_distance_to_front -= d_period[0];
            max_distance_to_front -= d_period[0];
         }
#endif
         // tbox::plog << "shifted ..........: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;
         if (min_distance_to_front <= 0 && max_distance_to_front >= 0) {
            // This cell has nodes on both sides of the front.  Tag it and the buffer_cells around it.
            hier::Box cell_and_buffer(cell_index, cell_index, blk0);
            cell_and_buffer.grow(buffer_cells);
            tag_data->fill(1, cell_and_buffer);
         }

      }
   }

   /*
    * Initialize U-value data.
    * The exact value of U increases by 1 across each front.
    */
   if (uval_data != 0) {
      t_uval->start();

      pdat::CellData<double>& uval(*uval_data);
      hier::Box uval_fill_box = uval.getGhostBox() * fill_box;
      uval.fill(0.0, uval_fill_box);
      const pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(uval_fill_box));
      for (pdat::CellData<double>::iterator ci = pdat::CellGeometry::begin(uval_fill_box);
           ci != ciend; ++ci) {
         const pdat::CellIndex& cindex = *ci;
         pdat::CellIndex squashed_cindex = cindex;
         squashed_cindex(0) = front_box.lower(0);
         double cellx = xlo[0] + dx[0] * (cindex(0) - pbox.lower(0) + 0.5);
         // Approximate cell's distance to front as average of its node distances.
         double dist_from_front = 0.0;
         if (d_dim == tbox::Dimension(2)) {
            dist_from_front = cellx - 0.5 * (
                  front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LowerLeft), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::UpperLeft), 0));
         } else if (d_dim == tbox::Dimension(3)) {
            dist_from_front = cellx - 0.25 * (
                  front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LLL), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LUL), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LLU), 0)
                  + front_x(pdat::NodeIndex(squashed_cindex, pdat::NodeIndex::LUU), 0));
         }
         const int k =
            (dist_from_front < 0 ? static_cast<int>(-dist_from_front / d_period[0] + 1) : 0);
         const int front_count = static_cast<int>(dist_from_front / d_period[0] + k) - k;
         uval(cindex) = front_count;
      }

      t_uval->stop();
   }

   /*
    * Initialize distance data.
    */
   if (dist_data != 0) {
      t_distance->start();

      pdat::NodeData<double>& dist_to_front(*dist_data);
      hier::Box dist_fill_box = dist_to_front.getGhostBox() * fill_box;
      pdat::NodeData<double>::iterator ni(pdat::NodeGeometry::begin(dist_fill_box));
      pdat::NodeData<double>::iterator niend(pdat::NodeGeometry::end(dist_fill_box));
      for ( ; ni != niend; ++ni) {
         const pdat::NodeIndex& index = *ni;
         pdat::NodeIndex front_index(index);
         front_index(0) = pbox.lower(0);
         dist_to_front(index) = xlo[0] + (index(0) - pbox.lower(0)) * dx[0]
            - front_x(front_index, 0);
      }

      t_distance->stop();
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */
#ifdef HAVE_HDF5
int SinusoidalFrontGenerator::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   /*
    * Register variables with plotter.
    */
   writer.registerDerivedPlotQuantity("Distance to front", "SCALAR", this, 1.0, "NODE");
   writer.registerDerivedPlotQuantity("U_Sinusoid", "SCALAR", this);
   writer.registerDerivedPlotQuantity("Tag value", "SCALAR", this);
   d_vis_owner_data.registerVariablesWithPlotter(writer);
   return 0;
}
#endif

/*
 ***********************************************************************
 ***********************************************************************
 */
bool SinusoidalFrontGenerator::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index,
   double simulation_time) const
{
   (void)depth_index;

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* xlo = patch_geom->getXLower();
   const double* dx = patch_geom->getDx();

   if (variable_name == "Distance to front") {
      pdat::NodeData<double> dist_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      dist_data.setTime(simulation_time);
      computeFrontsData(&dist_data, 0, 0, region,
         std::vector<double>(d_dim.getValue(), 0.0),
         xlo, dx);
      pdat::NodeData<double>::iterator ciend(pdat::NodeGeometry::end(patch.getBox()));
      for (pdat::NodeData<double>::iterator ci(pdat::NodeGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = dist_data(*ci);
      }
   } else if (variable_name == "U_Sinusoid") {
      pdat::CellData<double> u_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      u_data.setTime(simulation_time);
      computeFrontsData(0, &u_data, 0, region,
         std::vector<double>(d_dim.getValue(), 0.0),
         xlo, dx);
      pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(patch.getBox()));
      for (pdat::CellData<double>::iterator ci(pdat::CellGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = u_data(*ci);
      }
   } else if (variable_name == "Tag value") {
      pdat::CellData<int> tag_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      tag_data.setTime(simulation_time);
      computeFrontsData(0, 0, &tag_data, region,
         (static_cast<size_t>(patch.getPatchLevelNumber()) < d_buffer_distance.size() ?
          d_buffer_distance[patch.getPatchLevelNumber()] : d_buffer_distance.back()),
         xlo, dx);
      pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(patch.getBox()));
      for (pdat::CellData<double>::iterator ci(pdat::CellGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = tag_data(*ci);
      }
   } else {
      TBOX_ERROR("Unrecognized name " << variable_name);
   }
   return true;
}
