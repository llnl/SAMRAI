/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SphericalShellGenerator class implementation
 *
 ************************************************************************/
#include "SphericalShellGenerator.h"
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

SphericalShellGenerator::SphericalShellGenerator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<tbox::Database>& database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(),
   d_time_shift(0.0),
   d_radii(0),
   d_buffer_distance(1, std::vector<double>(dim.getValue(), 0.0))
{
   for (int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      d_init_center[i] = 0.0;
      d_velocity[i] = 0.0;
   }

   if (database) {

      if (database->isDouble("radii")) {

         d_radii = database->getDoubleVector("radii");

         if (static_cast<int>(d_radii.size()) % 2 != 0) {
            d_radii.push_back(tbox::MathUtilities<double>::getMax());
         }

         tbox::plog << "SphericalShellGenerator radii:\n";
         for (int i = 0; i < static_cast<int>(d_radii.size()); ++i) {
            tbox::plog << "\tradii[" << i << "] = " << d_radii[i] << '\n';
         }
      }

      if (database->isDouble("init_center")) {
         std::vector<double> tmpa = database->getDoubleVector("init_center");
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_init_center[d] = tmpa[d];
         }
      }
      if (database->isDouble("velocity")) {
         std::vector<double> tmpa = database->getDoubleVector("velocity");
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_velocity[d] = tmpa[d];
         }
      }

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

}

SphericalShellGenerator::~SphericalShellGenerator()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void SphericalShellGenerator::setTags(
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

      computeShellsData(0, tag_data.get(),
         tag_data->getBox(),
         (static_cast<size_t>(tag_ln) < d_buffer_distance.size() ?
          d_buffer_distance[tag_ln] : d_buffer_distance.back()),
         patch_geom->getXLower(),
         patch_geom->getDx());

   }

   exact_tagging = false;
}

void SphericalShellGenerator::setDomain(
   hier::BoxContainer& domain,
   double xlo[],
   double xhi[],
   int autoscale_base_nprocs,
   const tbox::SAMRAI_MPI& mpi)
{
   TBOX_ASSERT(!domain.empty());
   NULL_USE(xlo);
   NULL_USE(xhi);

   if (domain.size() != 1) {
      TBOX_ERROR("SphericalShellGenerator only supports single-box domains.");
   }

   hier::Box domain_box = domain.front();
   hier::IntVector tmp_intvec = domain_box.numberCells();
   const tbox::Dimension& dim = domain_box.getDim();

   double scale_factor = static_cast<double>(mpi.getSize()) / autoscale_base_nprocs;
   double linear_scale_factor = pow(scale_factor, 1.0 / dim.getValue());

   for (int d = 0; d < dim.getValue(); ++d) {
      // xhi[d] = xlo[d] + linear_scale_factor*(xhi[d]-xlo[d]);
      tmp_intvec(d) = static_cast<int>(0.5 + tmp_intvec(d) * linear_scale_factor);
   }
   tmp_intvec -= hier::IntVector::getOne(domain_box.getDim());
   tbox::plog << "SphericalShellGenerator::setDomain changing domain from "
              << domain_box << " to ";
   domain_box.setUpper(domain_box.lower() + tmp_intvec);
   tbox::plog << domain_box << '\n';

   domain.clear();
   domain.pushBack(domain_box);

}

void SphericalShellGenerator::resetHierarchyConfiguration(
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
void SphericalShellGenerator::computePatchData(
   const hier::Patch& patch,
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data,
   const hier::Box& fill_box) const
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom.get() != 0);

   const double* xlo = patch_geom->getXLower();
   const double* dx = patch_geom->getDx();

   if (tag_data) {
      computeShellsData(uval_data, tag_data,
         fill_box,
         (static_cast<size_t>(patch.getPatchLevelNumber()) < d_buffer_distance.size() ?
          d_buffer_distance[patch.getPatchLevelNumber()] : d_buffer_distance.back()),
         xlo, dx);
   } else {
      // Not computing tag => no tag buffer needed.
      computeShellsData(uval_data, tag_data,
         fill_box,
         std::vector<double>(d_dim.getValue(), 0.0),
         xlo, dx);
   }
}

/*
 * Compute the various data due to the shells.
 */
void SphericalShellGenerator::computeShellsData(
   pdat::CellData<double>* uval_data,
   pdat::CellData<int>* tag_data,
   const hier::Box& fill_box,
   const std::vector<double>& buffer_distance,
   const double xlo[],
   const double dx[]) const
{
   const int tag_val = 1;

   // Compute the buffer in terms of cells.
   hier::IntVector buffer_cells(d_dim);
   for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
      buffer_cells(d) = static_cast<int>(0.5 + buffer_distance[d] / dx[d]);
   }

   hier::Box pbox(d_dim), gbox(d_dim);
   double time = 0.0;
   if (tag_data != 0) {
      pbox = tag_data->getBox();
      gbox = tag_data->getGhostBox();
      time = tag_data->getTime() + d_time_shift;
   } else if (uval_data != 0) {
      pbox = uval_data->getBox();
      gbox = uval_data->getGhostBox();
      time = uval_data->getTime() + d_time_shift;
   }

   if (tag_data != 0) {
      /*
       * Compute radii of the nodes.  Tag cell if it has nodes farther than
       * the inner shell radius and nodes closer than the shell outer radius.
       */
      hier::Box radius_data_box = fill_box * gbox;
      pdat::NodeData<double> radius_data(radius_data_box, 1, buffer_cells);

      pdat::NodeData<int>::iterator niend(pdat::NodeGeometry::end(radius_data.getGhostBox()));
      for (pdat::NodeData<int>::iterator ni(pdat::NodeGeometry::begin(radius_data.getGhostBox()));
           ni != niend; ++ni) {
         const pdat::NodeIndex& idx = *ni;
         double r[SAMRAI::MAX_DIM_VAL];
         double rr = 0;
         for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
            r[d] = xlo[d]
               + dx[d] * (idx(d) - pbox.lower() (d))
               - (d_init_center[d] + time * d_velocity[d]);
            rr += r[d] * r[d];
         }
         rr = sqrt(rr);
         radius_data(idx) = rr;
      }

      hier::Box fbox = tag_data->getGhostBox() * fill_box;
      tag_data->getArrayData().fill(0, fill_box);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif


      pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(fbox));
      for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(fbox));
           ci != ciend; ++ci) {
         const pdat::CellIndex& cid = *ci;

         hier::Box check_box(cid, cid, pbox.getBlockId());
         check_box.grow(buffer_cells);
         pdat::NodeIterator node_itr_end(pdat::NodeGeometry::end(check_box));
         double min_node_radius = 1e20;
         double max_node_radius = 0;
         for (pdat::NodeIterator node_itr(pdat::NodeGeometry::begin(check_box));
              node_itr != node_itr_end; ++node_itr) {
            min_node_radius =
               tbox::MathUtilities<double>::Min(min_node_radius, radius_data(*node_itr));
            max_node_radius =
               tbox::MathUtilities<double>::Max(min_node_radius, radius_data(*node_itr));
         }
         for (int i = 0; i < static_cast<int>(d_radii.size()); i += 2) {
            if (d_radii[i] <= max_node_radius && min_node_radius < d_radii[i + 1]) {
               (*tag_data)(cid) = tag_val;
               break;
            }
         }
      }
   }

   if (uval_data != 0) {
      hier::Box fbox = uval_data->getGhostBox() * fill_box;
      uval_data->fill(static_cast<double>(d_radii.size()), fbox);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(fbox));
      for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(fbox));
           ci != ciend; ++ci) {
         const pdat::CellIndex& cid = *ci;

         double r[SAMRAI::MAX_DIM_VAL];
         double rr = 0;
         for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
            r[d] = xlo[d]
               + dx[d] * (cid(d) - pbox.lower() (d) + 0.5)
               - (d_init_center[d] + time * d_velocity[d]);
            rr += r[d] * r[d];
         }
         rr = sqrt(rr);

         for (int i = 0; i < static_cast<int>(d_radii.size()); ++i) {
            if (rr < d_radii[i]) {
               (*uval_data)(cid) = static_cast<double>(i);
               break;
            }
         }
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */
#ifdef HAVE_HDF5
int SphericalShellGenerator::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   writer.registerDerivedPlotQuantity("Tag value", "SCALAR", this);
   writer.registerDerivedPlotQuantity("U_Shells", "SCALAR", this);
   d_vis_owner_data.registerVariablesWithPlotter(writer);
   return 0;
}
#endif

/*
 ***********************************************************************
 ***********************************************************************
 */
bool SphericalShellGenerator::packDerivedDataIntoDoubleBuffer(
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

   if (variable_name == "U_Shells") {
      pdat::CellData<double> u_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      u_data.setTime(simulation_time);
      computeShellsData(&u_data, 0, region,
         std::vector<double>(d_dim.getValue(), 0.0),
         xlo, dx);
      pdat::CellData<double>::iterator ciend(pdat::CellGeometry::end(patch.getBox()));
      for (pdat::CellData<double>::iterator ci(pdat::CellGeometry::begin(patch.getBox()));
           ci != ciend; ++ci) {
         *(buffer++) = u_data(*ci);
      }
   } else if (variable_name == "Tag value") {

      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(patch_geom);

      pdat::CellData<int> tag_data(region, 1, hier::IntVector(d_dim, 0));
      tag_data.setTime(simulation_time);
      computeShellsData(
         0, &tag_data,
         region,
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
