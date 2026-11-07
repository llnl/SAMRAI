/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   ShrunkenLevelGenerator class implementation
 *
 ************************************************************************/
#include "ShrunkenLevelGenerator.h"
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

ShrunkenLevelGenerator::ShrunkenLevelGenerator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const std::shared_ptr<tbox::Database>& database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(),
   d_domain_scale_method('r'),
   d_shrink_distance(0)
{
   if (database) {

      /*
       * Input parameters to determine whether to tag by buffering
       * fronts or shrinking level, and by how much.
       */
      const std::string sname("shrink_distance_");
      for (int ln = 0; ; ++ln) {
         const std::string lnstr(tbox::Utilities::intToString(ln));

         // Look for buffer input first, then shrink input.
         const std::string snameln = sname + lnstr;

         std::vector<double> tmpa;

         if (database->isDouble(snameln)) {
            tmpa = database->getDoubleVector(snameln);
            if (static_cast<int>(tmpa.size()) != dim.getValue()) {
               TBOX_ERROR(snameln << " input parameter must have " << dim << " values");
            }
         }

         if (!tmpa.empty()) {
            d_shrink_distance.resize(d_shrink_distance.size() + 1);
            d_shrink_distance.back().insert(d_shrink_distance.back().end(),
               &tmpa[0],
               &tmpa[0] + static_cast<int>(tmpa.size()));
         } else {
            break;
         }

      }

      d_domain_scale_method =
         database->getCharWithDefault("domain_scale_method", d_domain_scale_method);

   }

}

ShrunkenLevelGenerator::~ShrunkenLevelGenerator()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void ShrunkenLevelGenerator::setTags(
   bool& exact_tagging,
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id)
{
   setTagsByShrinkingLevel(
      hierarchy,
      tag_ln,
      tag_data_id,
      hier::IntVector::getZero(d_dim),
      &d_shrink_distance[1][0]);
   exact_tagging = true;
}

void ShrunkenLevelGenerator::setDomain(
   hier::BoxContainer& domain,
   double xlo[],
   double xhi[],
   int autoscale_base_nprocs,
   const tbox::SAMRAI_MPI& mpi)
{
   TBOX_ASSERT(!domain.empty());
   NULL_USE(xlo);
   NULL_USE(xhi);

   if (d_domain_scale_method == 'r') {

      if (domain.size() != 1) {
         TBOX_ERROR("ShrunkenLevelGenerator resolution scaling only supports\n"
            << "single-box domains.");
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
      tbox::plog << "ShrunkenLevelGenerator::setDomain changing domain from "
                 << domain_box << " to ";
      domain_box.setUpper(domain_box.lower() + tmp_intvec);
      tbox::plog << domain_box << '\n';

      domain.clear();
      domain.pushBack(domain_box);

   } else {

      if (mpi.getSize() < autoscale_base_nprocs) {
         TBOX_ERROR("ShrunkenLevelGenerator::setDomain: When using\n"
            << "domain_scale_method = 't', autoscale_base_nprocs\n"
            << "cannot be smaller than number of processes.\n"
            << "Either set domain_scale_method = 'r', increase\n"
            << "autoscale_base_nprocs or run with mor processes.");
      }

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

}

void ShrunkenLevelGenerator::resetHierarchyConfiguration(
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

void ShrunkenLevelGenerator::setTagsByShrinkingLevel(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id,
   const hier::IntVector& shrink_cells,
   const double* shrink_distance)
{

   const tbox::Dimension dim(hierarchy->getDim());

   std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
         hierarchy->getGridGeometry()));
   TBOX_ASSERT(grid_geometry);

   const int tag_val = 1;

   const std::shared_ptr<hier::PatchLevel>& tag_level(
      hierarchy->getPatchLevel(tag_ln));

   const hier::BoxLevel& Ltag = *tag_level->getBoxLevel();

   /*
    * Compute shrinkage in terms of coarse cell count.  It should be
    * the largest of properly converted values for shrink_cells,
    * shrink_distance and nesting width.
    */
   const int nblocks =
      static_cast<int>(hierarchy->getGridGeometry()->getNumberBlocks());
   hier::IntVector shrink_width(dim, hierarchy->getProperNestingBuffer(tag_ln), nblocks);
   shrink_width.max(shrink_cells);

   const double *ref_dx = grid_geometry->getDx();
   for (int b = 0; b < nblocks; ++b) {
      const hier::IntVector& ref_ratio = Ltag.getRefinementRatio();
      for ( int i=0; i<dim.getValue(); ++i ) {
         double h = ref_dx[i]/ref_ratio(i);
         shrink_width(b,i) = tbox::MathUtilities<int>::Max(
            static_cast<int>(0.5 + shrink_distance[i]/h),
            shrink_width(b,i) );
      }
   }

   std::shared_ptr<hier::BoxLevel> tagfootprint;
   std::shared_ptr<hier::MappingConnector> Ltag_to_tagfootprint;
   const hier::Connector& Ltag_to_Ltag = Ltag.findConnector(Ltag,
         shrink_width,
         hier::CONNECTOR_CREATE);

   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts(tagfootprint,
      Ltag_to_tagfootprint,
      Ltag_to_Ltag,
      -shrink_width);
   tbox::plog << "Ltag_to_tagfootprint:\n" << Ltag_to_tagfootprint->format("Ltag->tagfootprint: ",
      2);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      std::shared_ptr<hier::Patch> patch = *pi;
      std::shared_ptr<pdat::CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch->getPatchData(tag_data_id)));
      TBOX_ASSERT(tag_data);

      tag_data->getArrayData().fillAll(0);

      if (!Ltag_to_tagfootprint->hasNeighborSet(patch->getBox().getBoxId())) {
         /*
          * Ltag_to_tagfootprint is a mapping Connector, so missing
          * neighbor set means the box has itself as its only
          * neighbor.
          */
         tag_data->getArrayData().fillAll(1);
      } else {
         hier::Connector::ConstNeighborhoodIterator ni =
            Ltag_to_tagfootprint->find(patch->getBox().getBoxId());

         for (hier::Connector::ConstNeighborIterator na = Ltag_to_tagfootprint->begin(ni);
              na != Ltag_to_tagfootprint->end(ni); ++na) {

            const hier::Box& tag_box = *na;
            tag_data->getArrayData().fillAll(tag_val, tag_box);

         }
      }

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
#ifdef HAVE_HDF5
int ShrunkenLevelGenerator::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   d_vis_owner_data.registerVariablesWithPlotter(writer);
   return 0;
}
#endif

/*
 ***********************************************************************
 ***********************************************************************
 */
bool ShrunkenLevelGenerator::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index,
   double simulation_time) const
{
   (void)buffer;
   (void)patch;
   (void)region;
   (void)variable_name;
   (void)depth_index;
   (void)simulation_time;
   return true;
}
