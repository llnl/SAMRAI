/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:  Manages data on stencils at coarse-fine boundaries
 *
 ************************************************************************/
#include "SAMRAI/xfer/CompositeBoundarySchedule.h"

#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/tbox/NVTXUtilities.h"


namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

CompositeBoundarySchedule::CompositeBoundarySchedule(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarse_level_num,
   int stencil_width,
   const std::set<int>& data_ids)
: d_hierarchy(hierarchy),
  d_stencil_width(stencil_width),
  d_data_ids(data_ids),
  d_coarse_level_num(coarse_level_num)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(coarse_level_num <= hierarchy->getFinestLevelNumber());
   TBOX_ASSERT(stencil_width >= 1);
   registerDataIds();
   createStencilForLevel();
}

/*
 *************************************************************************
 *
 * Destructor explicitly deallocates PatchData
 *
 *************************************************************************
 */

CompositeBoundarySchedule::~CompositeBoundarySchedule()
{
   if (d_stencil_level.get() != 0) { 
      for (std::set<int>::const_iterator ditr = d_data_ids.begin();
           ditr != d_data_ids.end(); ++ditr) {
      
         if (d_stencil_level->checkAllocated(*ditr)) {
            d_stencil_level->deallocatePatchData(*ditr);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Creates a stencil at the coarse-fine boundary of the given level
 *
 *************************************************************************
 */

void
CompositeBoundarySchedule::createStencilForLevel()
{
   int level_num = d_coarse_level_num;
   TBOX_ASSERT(level_num >= 0 &&
               level_num <= d_hierarchy->getFinestLevelNumber());

   const tbox::Dimension& dim = d_hierarchy->getDim();

   d_stencil_level.reset();

   if (level_num < d_hierarchy->getFinestLevelNumber()) {

      hier::IntVector stencil_vec(dim, d_stencil_width);

      hier::HierarchyNeighbors hier_nbrs(*d_hierarchy,
                                         level_num,
                                         level_num+1,
                                         false,
                                         d_stencil_width);

      const std::shared_ptr<hier::BaseGridGeometry>& grid_geom =
         d_hierarchy->getGridGeometry();

      const std::shared_ptr<hier::PatchLevel>& level =
         d_hierarchy->getPatchLevel(level_num);
      const std::shared_ptr<hier::PatchLevel>& finer_level =
         d_hierarchy->getPatchLevel(level_num+1);

      std::shared_ptr<hier::BoxLevel> stencil_box_level =
         std::make_shared<hier::BoxLevel>(
            finer_level->getRatioToLevelZero(),
            grid_geom);

      hier::IntVector ratio(d_hierarchy->getRatioToCoarserLevel(level_num+1)); 

      hier::IntVector s_to_c_width(
         d_hierarchy->getRequiredConnectorWidth(level_num+1, level_num));
      s_to_c_width.max(stencil_vec);

      unsigned int num_blocks = d_hierarchy->getNumberBlocks();
      for (hier::BlockId::block_t blk = 0; blk < num_blocks; ++blk) {
         for (unsigned int i = 0; i < dim.getValue(); ++i) {
            while (s_to_c_width(blk,i) % ratio(blk,i) != 0) {
               ++(s_to_c_width(blk,i));
            }             
         }
      }

      hier::IntVector c_to_s_width(s_to_c_width / ratio);

      hier::Connector coarse_to_stencil(
         *level->getBoxLevel(),
         *stencil_box_level,
         c_to_s_width); 

      hier::Connector stencil_to_coarse(
         *stencil_box_level,
         *level->getBoxLevel(),
         s_to_c_width);

      hier::LocalId stencil_local_id(-1);
      hier::LocalId coarse_local_id(-1);
      for (hier::PatchLevel::iterator itr = level->begin(); itr != level->end();
           ++itr) {
         const std::shared_ptr<hier::Patch>& patch(*itr);
         const hier::Box& patch_box = patch->getBox();
         const hier::BoxId& box_id = patch_box.getBoxId();

         const hier::BoxContainer& finer_nbrs =
            hier_nbrs.getFinerLevelNeighbors(patch_box, level_num);

         hier::BoxContainer finer_tree(finer_nbrs);
         finer_tree.makeTree(&(*grid_geom));

         const hier::IntVector& level_ratio(
            d_hierarchy->getRatioToCoarserLevel(level_num+1));
         hier::BoxContainer uncovered_fine(patch_box);
         uncovered_fine.refine(level_ratio);
         uncovered_fine.removeIntersections(finer_level->getRatioToLevelZero(),
                                        finer_tree,
                                        true);
         uncovered_fine.coalesce();

         hier::BoxContainer& uncovered_boxes =
            d_patch_to_uncovered[box_id];

         for (hier::BoxContainer::iterator uitr = uncovered_fine.begin();
              uitr != uncovered_fine.end(); ++uitr) { 
            ++coarse_local_id;
            hier::Box add_box(*uitr, coarse_local_id, box_id.getOwnerRank());
            add_box.coarsen(level_ratio);
            uncovered_boxes.insert(uncovered_boxes.end(), add_box);
         }

         for (hier::BoxContainer::iterator uitr = uncovered_boxes.begin();
              uitr != uncovered_boxes.end(); ++uitr) {
            hier::BoxContainer grow_boxes(*uitr);
            grow_boxes.refine(level_ratio);
            grow_boxes.grow(stencil_vec);
            grow_boxes.intersectBoxes(finer_level->getRatioToLevelZero(),
                                      finer_tree, true);
            grow_boxes.coalesce();

            for (hier::BoxContainer::iterator gitr = grow_boxes.begin();
                 gitr != grow_boxes.end(); ++gitr) {
               ++stencil_local_id;
               hier::Box stencil_box(*gitr, stencil_local_id, box_id.getOwnerRank());
               stencil_box_level->addBoxWithoutUpdate(stencil_box); 
               d_uncovered_to_stencil[uitr->getBoxId()].insert(
                  stencil_box.getBoxId()); 

               coarse_to_stencil.insertLocalNeighbor(stencil_box, box_id);
               stencil_to_coarse.insertLocalNeighbor(patch_box,
                                                     stencil_box.getBoxId());
            }
         }
      }

      stencil_box_level->finalize();

      coarse_to_stencil.setTranspose(&stencil_to_coarse, false);
      const hier::Connector& fine_to_coarse =
         finer_level->getBoxLevel()->findConnectorWithTranspose(
            *level->getBoxLevel(),
            s_to_c_width,
            c_to_s_width,
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      std::shared_ptr<hier::Connector> fine_to_stencil;
      hier::OverlapConnectorAlgorithm oca;
      oca.bridge(
         fine_to_stencil,
         fine_to_coarse,
         coarse_to_stencil,
         hier::IntVector::getZero(dim),
         true);

      const hier::Connector& fine_to_fine =
         finer_level->getBoxLevel()->findConnector(
            *finer_level->getBoxLevel(),
            d_hierarchy->getRequiredConnectorWidth(level_num+1, level_num+1),
            hier::CONNECTOR_IMPLICIT_CREATION_RULE,
            true);

      hier::BoxLevelConnectorUtils blcu;
      blcu.addPeriodicImagesAndRelationships(
         *stencil_box_level,
         fine_to_stencil->getTranspose(),
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         fine_to_fine); 

      finer_level->getBoxLevel()->cacheConnector(fine_to_stencil);

      d_stencil_level = std::make_shared<hier::PatchLevel>(
         stencil_box_level,
         grid_geom,
         level->getPatchDescriptor());
      d_stencil_level->setLevelNumber(level_num+1);

      d_stencil_created = true;

      setUpSchedule();

   } else {
      d_stencil_created = true;
      setUpSchedule();
   } 


}

void CompositeBoundarySchedule::registerDataIds()
{
   std::shared_ptr< xfer::VariableFillPattern >
      interior_fill(new PatchInteriorVariableFillPattern(d_hierarchy->getDim()));

   for (std::set<int>::const_iterator itr = d_data_ids.begin();
        itr != d_data_ids.end(); ++itr) {
      TBOX_ASSERT(*itr >= 0);
      int data_id = *itr;

      d_refine_algorithm.registerRefine(data_id, data_id, data_id,
         std::shared_ptr<hier::RefineOperator>(), interior_fill);
   }
}


/*
 *************************************************************************
 *
 * Private method sets up schedule for communicating data onto the stencil
 *
 *************************************************************************
 */

void
CompositeBoundarySchedule::setUpSchedule()
{
   TBOX_ASSERT(d_stencil_created);

   d_refine_schedule.reset();

   if (d_stencil_level) {
      std::shared_ptr< xfer::PatchLevelFillPattern >
         interior_level_fill(new PatchLevelInteriorFillPattern());
 
      d_refine_schedule =
         d_refine_algorithm.createSchedule(
            interior_level_fill,
            d_stencil_level,
            d_hierarchy->getPatchLevel(d_coarse_level_num+1));
   }
}

/*
 *************************************************************************
 *
 * Fill patch data on the stencil
 *
 *************************************************************************
 */
void
CompositeBoundarySchedule::fillData(double fill_time)
{
   if (!d_stencil_created) {
      TBOX_ERROR("CompositeBoundarySchedule::fillData error:  No stencil for boundary of Level " << d_coarse_level_num << " created." << std::endl);
   }
   if (d_refine_schedule == 0 && d_coarse_level_num != d_hierarchy->getFinestLevelNumber()) {
      TBOX_ERROR("CompositeBoundarySchedule::fillData error:  No schedule for filling stencil for boundary of Level " << d_coarse_level_num << " created." << std::endl);
   }

   if (d_stencil_level && d_refine_schedule) {
      for (std::set<int>::const_iterator ditr = d_data_ids.begin();
           ditr != d_data_ids.end(); ++ditr) {

         if (!d_stencil_level->checkAllocated(*ditr)) {
            d_stencil_level->allocatePatchData(*ditr);
         }
      }

      d_refine_schedule->fillData(fill_time);

   }

   const std::shared_ptr<hier::PatchLevel>& level =
      d_hierarchy->getPatchLevel(d_coarse_level_num);

   for (std::set<int>::const_iterator ditr = d_data_ids.begin();
        ditr != d_data_ids.end(); ++ditr) {

      const int& data_id = *ditr;
      d_data_map[data_id].clear();

      for (hier::PatchLevel::iterator itr = level->begin(); itr != level->end();
           ++itr) {

         const std::shared_ptr<hier::Patch>& patch(*itr);
         const hier::BoxId& box_id = patch->getBox().getBoxId();
         VectorPatchData& data_vec =
            d_data_map[data_id][box_id];
         data_vec.clear();

         std::map<hier::BoxId, std::shared_ptr<hier::PatchData> >&
            stencil_map = d_stencil_to_data[data_id];
         if (d_stencil_level) {
            const hier::BoxContainer& uncovered_boxes =
               d_patch_to_uncovered[box_id];

            for (hier::BoxContainer::const_iterator uitr = uncovered_boxes.begin();
                 uitr != uncovered_boxes.end(); ++uitr) {
               const std::set<hier::BoxId>& stencil_ids =
                  d_uncovered_to_stencil[uitr->getBoxId()];

               for (std::set<hier::BoxId>::const_iterator sitr = stencil_ids.begin();
                    sitr != stencil_ids.end(); ++sitr) {

                  data_vec.push_back(
                     d_stencil_level->getPatch(*sitr)->getPatchData(data_id)); 

                  stencil_map[*sitr] = data_vec.back();

               } 
            } 
         } 
      }
   }

   d_stencil_filled = true;  
}

/*
 *************************************************************************
 *
 * Get a vector of stencil PatchData
 *
 *************************************************************************
 */

const std::vector<std::shared_ptr<hier::PatchData> >&
CompositeBoundarySchedule::getBoundaryPatchData(
   const hier::Patch& patch,
   int data_id) const
{
   if (d_data_ids.find(data_id) == d_data_ids.end()) {
      TBOX_ERROR("CompositeBoundarySchedule::getBoundaryPatchData data_id " << data_id << " not registered with CompositeBoundarySchedule" << std::endl);
   }

   int level_num = d_coarse_level_num;
   if (patch.getPatchLevelNumber() != level_num) {
      TBOX_ERROR("CompositeBoundarySchedule::getBoundaryPatchData Patch is on level " << patch.getPatchLevelNumber() << " while CompositeBoundarySchedule was created for level " <<  level_num << "." << std::endl);
   }

   if (!d_stencil_filled) {
      TBOX_ERROR("CompositeBoundarySchedule::getBoundaryPatchData error:  stencil for boundary of Level " << level_num << " not filled." << std::endl);
   }

   std::map<int, MapBoxIdPatchData>::const_iterator data_map_itr = 
      d_data_map.find(data_id);

   if (data_map_itr == d_data_map.end()) {
      TBOX_ERROR("CompositeBoundarySchedule::getBoundaryPatchData error:  stencil data for patch data id " << data_id << " not found for level " << level_num << std::endl);
   }

   const MapBoxIdPatchData& inner_map = data_map_itr->second;  

   const hier::BoxId& box_id = patch.getBox().getBoxId();
   MapBoxIdPatchData::const_iterator inner_map_itr = inner_map.find(box_id);

   if (inner_map_itr == inner_map.end()) {
      TBOX_ERROR("CompositeBoundarySchedule::getBoundaryPatchData error:  stencil data for patch data id " << data_id << " not found for Patch " << box_id << " on level " << level_num << std::endl);
   }

   const std::vector<std::shared_ptr<hier::PatchData> >& data_vec =
      inner_map_itr->second;

   return data_vec;

}

std::vector<std::shared_ptr<hier::PatchData> >
CompositeBoundarySchedule::getDataForUncoveredBox(
   const hier::BoxId& uncovered_id,
   int data_id) const
{
   if (d_data_ids.find(data_id) == d_data_ids.end()) {
      TBOX_ERROR("CompositeBoundarySchedule:: data_id " << data_id << " not registered with CompositeBoundarySchedule" << std::endl);
   }

   if (!d_stencil_filled) {
      TBOX_ERROR("CompositeBoundarySchedule:: error:  stencil not filled." << std::endl);
   }

   std::map<int, std::map<hier::BoxId, std::shared_ptr<hier::PatchData> > >::const_iterator data_map_itr = d_stencil_to_data.find(data_id);
   if (data_map_itr == d_stencil_to_data.end()) {
      TBOX_ERROR("CompositeBoundarySchedule:: error:  stencil data for patch data id " << data_id << " not found." << std::endl);
   }

   const std::map<hier::BoxId, std::shared_ptr<hier::PatchData> >& data_map =
      data_map_itr->second;

   std::vector<std::shared_ptr<hier::PatchData> > data_vec;

   std::map<hier::BoxId, std::set<hier::BoxId> >::const_iterator st_itr =
      d_uncovered_to_stencil.find(uncovered_id);

   if (st_itr != d_uncovered_to_stencil.end()) {
      const std::set<hier::BoxId>& stencil_ids = st_itr->second;

      for (std::set<hier::BoxId>::const_iterator itr = stencil_ids.begin();
           itr != stencil_ids.end(); ++itr) {
         std::map<hier::BoxId, std::shared_ptr<hier::PatchData> >::const_iterator pd_itr = data_map.find(*itr);
         if (pd_itr == data_map.end()) {
            TBOX_ERROR("CompositeBoundarySchedule:: error:  No patch data for stencil box "  << *itr << " found." << std::endl); 
         }
         data_vec.push_back(pd_itr->second);
      }
   }

   return data_vec; 
}





}
}
