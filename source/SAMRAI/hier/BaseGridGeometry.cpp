/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/
#include "SAMRAI/hier/BaseGridGeometry.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/BoxTree.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/SingularityFinder.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <map>
#include <stdlib.h>
#include <vector>
#include <memory>

#define HIER_GRID_GEOMETRY_VERSION (3)

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

tbox::StartupShutdownManager::Handler
BaseGridGeometry::s_initialize_handler(
   BaseGridGeometry::initializeCallback,
   0,
   0,
   BaseGridGeometry::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<tbox::Timer> BaseGridGeometry::t_find_patches_touching_boundaries;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_touching_boundaries_init;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_touching_boundaries_loop;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_set_geometry_on_patches;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_set_boundary_boxes;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_set_geometry_data_on_patches;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_compute_boundary_boxes_on_level;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_get_boundary_boxes;
std::shared_ptr<tbox::Timer> BaseGridGeometry::t_adjust_multiblock_patch_level_boundaries;

/*
 *************************************************************************
 *
 * Constructors for BaseGridGeometry.  Both set up operator
 * handlers.  However, one initializes data members based on arguments.
 * The other initializes the object based on input database information.
 *
 *************************************************************************
 */
BaseGridGeometry::BaseGridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const std::shared_ptr<tbox::Database>& input_db,
   bool allow_multiblock):
   d_transfer_operator_registry(
      std::make_shared<TransferOperatorRegistry>(dim)),
   d_dim(dim),
   d_object_name(object_name),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_periodic_shift_catalog(d_dim),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_ratio_to_level_zero(1, IntVector::getOne(d_dim)),
   d_ratios_are_isotropic(true),
   d_has_enhanced_connectivity(false)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);

   tbox::RestartManager::getManager()->registerRestartItem(getObjectName(),
      this);

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart, allow_multiblock);

}

/*
 *************************************************************************
 *
 * Constructors for BaseGridGeometry.  Both set up operator
 * handlers.  However, one initializes data members based on arguments.
 * The other initializes the object based on input database information.
 *
 *************************************************************************
 */
BaseGridGeometry::BaseGridGeometry(
   const std::string& object_name,
   BoxContainer& domain):
   d_transfer_operator_registry(
      std::make_shared<TransferOperatorRegistry>(
         (*(domain.begin())).getDim())),
   d_dim((*(domain.begin())).getDim()),
   d_object_name(object_name),
   d_physical_domain(),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_periodic_shift_catalog(d_dim),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_number_of_block_singularities(0),
   d_ratio_to_level_zero(1, IntVector::getOne(d_dim)),
   d_ratios_are_isotropic(true),
   d_has_enhanced_connectivity(false)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!domain.empty());

   tbox::RestartManager::getManager()->
   registerRestartItem(getObjectName(), this);

   LocalId local_id(0);
   std::set<BlockId::block_t> block_numbers;
   for (BoxContainer::iterator itr = domain.begin(); itr != domain.end();
        ++itr) {
      block_numbers.insert(itr->getBlockId().getBlockValue());
      BoxId box_id(local_id++, 0);
      itr->setId(box_id);
   }
   d_number_blocks = block_numbers.size();
   if (d_ratio_to_level_zero[0].getNumBlocks() != d_number_blocks) {
      d_ratio_to_level_zero[0] =
         IntVector(IntVector::getOne(d_dim), d_number_blocks);
   }

   d_block_neighbors.resize(d_number_blocks);

   setPhysicalDomain(domain, d_number_blocks);

}

BaseGridGeometry::BaseGridGeometry(
   const std::string& object_name,
   BoxContainer& domain,
   const std::shared_ptr<TransferOperatorRegistry>& op_reg):
   d_transfer_operator_registry(op_reg),
   d_dim((*(domain.begin())).getDim()),
   d_object_name(object_name),
   d_physical_domain(),
   d_periodic_shift(IntVector::getZero(d_dim)),
   d_periodic_shift_catalog(d_dim),
   d_max_data_ghost_width(IntVector(d_dim, -1)),
   d_number_of_block_singularities(0),
   d_ratio_to_level_zero(1, IntVector::getOne(d_dim)),
   d_ratios_are_isotropic(true),
   d_has_enhanced_connectivity(false)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!domain.empty());

   tbox::RestartManager::getManager()->
   registerRestartItem(getObjectName(), this);

   LocalId local_id(0);
   std::set<BlockId::block_t> block_numbers;
   for (BoxContainer::iterator itr = domain.begin(); itr != domain.end();
        ++itr) {
      block_numbers.insert(itr->getBlockId().getBlockValue());
      BoxId box_id(local_id++, 0);
      itr->setId(box_id);
   }
   d_number_blocks = block_numbers.size();
   if (d_ratio_to_level_zero[0].getNumBlocks() != d_number_blocks) {
      d_ratio_to_level_zero[0] = 
         IntVector(IntVector::getOne(d_dim), d_number_blocks);
   }
   d_block_neighbors.resize(d_number_blocks);

   setPhysicalDomain(domain, d_number_blocks);
}

/*
 *************************************************************************
 *
 * Destructor.
 *
 *************************************************************************
 */

BaseGridGeometry::~BaseGridGeometry()
{
   tbox::RestartManager::getManager()->unregisterRestartItem(getObjectName());
}

/*
 *************************************************************************
 *
 * Compute boundary boxes for all patches in patch level.  The domain
 * array describes the interior of the level index space.  Note that
 * boundaries is assumed to be an array of DIM * #patches Arrays of
 * BoundaryBoxes.
 *
 *************************************************************************
 */

void
BaseGridGeometry::computeBoundaryBoxesOnLevel(
   std::map<BoxId, PatchBoundaries>& boundaries,
   const PatchLevel& level,
   const IntVector& periodic_shift,
   const IntVector& ghost_width,
   const std::vector<BoxContainer>& domain,
   bool do_all_patches) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(d_dim,
      level,
      periodic_shift,
      ghost_width);

   t_compute_boundary_boxes_on_level->start();

   TBOX_ASSERT(ghost_width >= IntVector::getZero(ghost_width.getDim()));
#ifdef DEBUG_CHECK_ASSERTIONS
   int num_per_dirs = 0;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      if (periodic_shift(i)) {
         ++num_per_dirs;
      }
   }
   if (num_per_dirs > 0) {
      TBOX_ASSERT(domain.size() == 1);
   }
#endif

   for (PatchLevel::iterator ip(level.begin()); ip != level.end(); ++ip) {
      const std::shared_ptr<Patch>& patch = *ip;
      const BoxId& patch_id = patch->getBox().getBoxId();
      const BlockId::block_t& block_num =
         patch->getBox().getBlockId().getBlockValue();

      if (patch->getPatchGeometry()->getTouchesRegularBoundary() ||
          do_all_patches) {

         const Box& box(patch->getBox());

         /*
          * patch_boundaries is an array of DIM BoxContainers for each patch.
          * patch_boundaries[DIM-1] will store boundary boxes of the
          * nodetype. If DIM > 1, patch_boundaries[DIM-2] will store
          * boundary boxes of the edge type, and if DIM > 2,
          * patch_boundaries[DIM-3] will store boundary boxes of the face
          * type.
          */

         /*
          * Create new map element if one does not exist.
          * Note can't use [] as this requires a default ctor which we do
          * not have for PatchBoundaries.
          */
         std::map<BoxId, PatchBoundaries>::iterator iter(
            boundaries.find(patch_id));
         if (iter == boundaries.end()) {
            std::pair<BoxId, PatchBoundaries> new_boundaries(patch_id,
                                                             PatchBoundaries(d_dim));
            iter = boundaries.insert(iter, new_boundaries);
         }
         getBoundaryBoxes((*iter).second, box, domain[block_num],
            ghost_width, periodic_shift);

#ifdef DEBUG_CHECK_ASSERTIONS
         for (int j = 0; j < d_dim.getValue(); ++j) {
            iter = (boundaries.find(patch_id));
            TBOX_ASSERT(iter != boundaries.end());
            for (int k = 0;
                 k < static_cast<int>(((*iter).second)[j].size()); ++k) {
               TBOX_ASSERT(checkBoundaryBox(((*iter).second)[j][k], *patch,
                     domain[block_num], num_per_dirs, ghost_width));
            }
         }
#endif
      }
   }
   t_compute_boundary_boxes_on_level->stop();
}

/*
 *************************************************************************
 *
 * For each patch in the level, use box intersection operation to
 * determine what kind of boundaries, if any the patch touches.  Call
 * Patch functions to set flags that store this information once it
 * is found.
 *
 *************************************************************************
 */

void
BaseGridGeometry::findPatchesTouchingBoundaries(
   std::map<BoxId, TwoDimBool>& touches_regular_bdry,
   std::map<BoxId, TwoDimBool>& touches_periodic_bdry,
   const PatchLevel& level) const
{
   t_find_patches_touching_boundaries->start();

   t_touching_boundaries_init->start();
   touches_regular_bdry.clear();
   touches_periodic_bdry.clear();
   t_touching_boundaries_init->stop();

   BoxContainer tmp_refined_periodic_domain_tree;
   if (level.getRatioToLevelZero() != IntVector::getOne(level.getDim())) {
      tmp_refined_periodic_domain_tree = d_domain_with_images;
      tmp_refined_periodic_domain_tree.refine(level.getRatioToLevelZero());
      tmp_refined_periodic_domain_tree.makeTree(this);
   }

   t_touching_boundaries_loop->start();
   for (PatchLevel::iterator ip(level.begin()); ip != level.end(); ++ip) {
      const std::shared_ptr<Patch>& patch = *ip;
      const Box& box(patch->getBox());

      std::map<BoxId, TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getBoxId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, TwoDimBool>(ip->getBox().getBoxId(), TwoDimBool(d_dim)));
      }

      std::map<BoxId, TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getBoxId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, TwoDimBool>(ip->getBox().getBoxId(), TwoDimBool(d_dim)));
      }

      computeBoxTouchingBoundaries(
         (*iter_touches_regular_bdry).second,
         (*iter_touches_periodic_bdry).second,
         box,
         level.getRatioToLevelZero(),
         tmp_refined_periodic_domain_tree.empty() ?
         d_physical_domain :
         //d_domain_search_tree :
         tmp_refined_periodic_domain_tree);
   }
   t_touching_boundaries_loop->stop();
   t_find_patches_touching_boundaries->stop();
}

void
BaseGridGeometry::computeBoxTouchingBoundaries(
   TwoDimBool& touches_regular_bdry,
   TwoDimBool& touches_periodic_bdry,
   const Box& box,
   const IntVector& refinement_ratio,
   const BoxContainer& refined_periodic_domain_tree) const
{

   /*
    * Create a list of boxes inside a layer of one cell outside the
    * patch.  Remove the intersections with the domain's interior, so that only
    * boxes outside the physical domain (if any) remain in the list.
    */
   BoxContainer bdry_list(box);
   bdry_list.grow(IntVector::getOne(d_dim));
   bdry_list.removeIntersections(refinement_ratio,
      refined_periodic_domain_tree);
   const bool touches_any_boundary = !bdry_list.empty();

   if (!touches_any_boundary) {
      for (int d = 0; d < d_dim.getValue(); ++d) {
         touches_regular_bdry(d, 0) = touches_periodic_bdry(d, 0) =
               touches_regular_bdry(d, 1) = touches_periodic_bdry(d, 1) = false;
      }
   } else {
      bool bdry_located = false;
      for (tbox::Dimension::dir_t nd = 0; nd < d_dim.getValue(); ++nd) {
         BoxContainer lower_list(bdry_list);
         BoxContainer upper_list(bdry_list);

         Box test_box(box);

         test_box.growLower(nd, 1);
         lower_list.intersectBoxes(test_box); // performance ok.  lower_list is short.

         test_box = box;
         test_box.growUpper(nd, 1);
         upper_list.intersectBoxes(test_box); // performance ok.  upper_list is short.

         if (!lower_list.empty()) {
            // Touches regular or periodic bdry on lower side.
            touches_periodic_bdry(nd, 0) = (d_periodic_shift(nd) != 0);
            touches_regular_bdry(nd, 0) = (d_periodic_shift(nd) == 0);
            bdry_located = true;
         }

         if (!upper_list.empty()) {
            // Touches regular or periodic bdry on upper side.
            touches_periodic_bdry(nd, 1) = (d_periodic_shift(nd) != 0);
            touches_regular_bdry(nd, 1) = (d_periodic_shift(nd) == 0);
            bdry_located = true;
         }
      }

      /*
       * By this point, bdry_located will have been set to true almost
       * every time whenever touches_any_boundary is true.  The only way
       * it will not be true is if the domain is not a parallelpiped, and
       * the patch touches the boundary only at a location such as the
       * concave corner of an L-shaped domain.
       */
      if (!bdry_located) {
         for (tbox::Dimension::dir_t nd = 0; nd < d_dim.getValue(); ++nd) {
            touches_periodic_bdry(nd, 0) = touches_periodic_bdry(nd, 1) = false;

            bool lower_side = false;
            bool upper_side = false;
            for (BoxContainer::iterator bl = bdry_list.begin();
                 bl != bdry_list.end(); ++bl) {
               if (bl->lower() (nd) < box.lower(nd)) {
                  lower_side = true;
               }
               if (bl->upper() (nd) > box.upper(nd)) {
                  upper_side = true;
               }
               if (lower_side && upper_side) {
                  break;
               }
            }
            touches_regular_bdry(nd, 0) = lower_side;
            touches_regular_bdry(nd, 1) = upper_side;
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Set geometry data for each patch on level.
 *
 *************************************************************************
 */

void
BaseGridGeometry::setGeometryOnPatches(
   PatchLevel& level,
   const IntVector& ratio_to_level_zero,
   const std::map<BoxId, TwoDimBool>& touches_regular_bdry,
   const bool defer_boundary_box_creation)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, level, ratio_to_level_zero);

   t_set_geometry_on_patches->start();

   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != 0);
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim.getValue() > 1) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         int i;
         for (i = 0; i < d_dim.getValue(); ++i) {
            bool pos0 = ratio_to_level_zero(b,i) > 0;
            bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
            TBOX_ASSERT(pos0 == pos1
               || (ratio_to_level_zero(b,i) == 1)
               || (ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
         }
      }
   }
#endif

   t_set_geometry_data_on_patches->start();
   for (PatchLevel::iterator ip(level.begin()); ip != level.end(); ++ip) {
      const std::shared_ptr<Patch>& patch = *ip;
      setGeometryDataOnPatch(*patch, ratio_to_level_zero,
         (*touches_regular_bdry.find(ip->getBox().getBoxId())).second);
   }
   t_set_geometry_data_on_patches->stop();

   if (!defer_boundary_box_creation) {
      setBoundaryBoxes(level);
   }

   t_set_geometry_on_patches->stop();
}

/*
 *************************************************************************
 *
 * Set boundary boxes for each patch on level.
 *
 *************************************************************************
 */

void
BaseGridGeometry::setBoundaryBoxes(
   PatchLevel& level)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, level);

   t_set_boundary_boxes->start();
   std::map<BoxId, PatchBoundaries> boundaries;

   const std::vector<BoxContainer>& domain(level.getPhysicalDomainArray());

   IntVector ghost_width(
      level.getPatchDescriptor()->getMaxGhostWidth(d_dim));

   if (d_max_data_ghost_width != -1 &&
       !(ghost_width <= d_max_data_ghost_width)) {

      TBOX_ERROR("Error in BaseGridGeometry object with name = "
         << d_object_name << ": in computeMaxGhostWidth():  "
         << "Cannot add variables and increase maximum ghost "
         << "width after creating the BaseGridGeometry!" << std::endl);
   }

   d_max_data_ghost_width = ghost_width;

   computeBoundaryBoxesOnLevel(
      boundaries,
      level,
      getPeriodicShift(IntVector::getOne(d_dim)),
      d_max_data_ghost_width,
      domain);

   for (std::map<BoxId, PatchBoundaries>::iterator mi = boundaries.begin();
        mi != boundaries.end(); ++mi) {
      std::shared_ptr<Patch> patch(level.getPatch((*mi).first));
      patch->getPatchGeometry()->setBoundaryBoxesOnPatch((*mi).second.getVectors());
   }

   t_set_boundary_boxes->stop();
}

/*
 *************************************************************************
 *
 * Create PatchGeometry geometry object, initializing its
 * boundary and assigning it to the given patch.
 *
 *************************************************************************
 */

void
BaseGridGeometry::setGeometryDataOnPatch(
   Patch& patch,
   const IntVector& ratio_to_level_zero,
   const PatchGeometry::TwoDimBool& touches_regular_bdry) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT_DIM_OBJDIM_EQUALITY3(dim, patch, ratio_to_level_zero,
      touches_regular_bdry);

   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != 0); 
   if (d_dim.getValue() > 1) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         int i;
         for (i = 0; i < d_dim.getValue(); ++i) {
            bool pos0 = ratio_to_level_zero(b,i) > 0;
            bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
            TBOX_ASSERT(pos0 == pos1
               || (ratio_to_level_zero(b,i) == 1)
               || (ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
         }
      }
   }
#endif

   std::shared_ptr<PatchGeometry> geometry(
      std::make_shared<PatchGeometry>(
         ratio_to_level_zero,
         touches_regular_bdry,
         patch.getBox().getBlockId()));

   patch.setPatchGeometry(geometry);

}

/*
 *************************************************************************
 * Checks to see if the version number for the class is the same as
 * as the version number of the restart file.
 * If they are equal, then the data from the database are read to local
 * variables and the setPhysicalDomain() method is called.
 *
 *************************************************************************
 */
void
BaseGridGeometry::getFromRestart()
{
   std::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(getObjectName())) {
      TBOX_ERROR("Restart database corresponding to "
         << getObjectName() << " not found in the restart file." << std::endl);
   }
   std::shared_ptr<tbox::Database> db(
      restart_db->getDatabase(getObjectName()));

   const tbox::Dimension dim(getDim());

   int ver = db->getInteger("HIER_GRID_GEOMETRY_VERSION");
   if (ver != HIER_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR(
         getObjectName() << ":  "
                         << "Restart file version is different than class version."
                         << std::endl);
   }

   d_number_blocks = static_cast<size_t>(db->getInteger("num_blocks"));
   if (d_ratio_to_level_zero[0].getNumBlocks() != d_number_blocks) {
      d_ratio_to_level_zero[0] = 
         IntVector(IntVector::getOne(d_dim), d_number_blocks);
   }

   d_singularity.resize(d_number_blocks);
   d_block_neighbors.resize(d_number_blocks);

   std::string domain_name;
   BoxContainer domain;
   LocalId local_id(0);

   for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
      std::string blk_string =
         tbox::Utilities::intToString(static_cast<int>(b));

      domain_name = "domain_boxes_" + blk_string;
      std::vector<tbox::DatabaseBox> db_box_vector =
         db->getDatabaseBoxVector(domain_name);
      BoxContainer block_domain_boxes(db_box_vector);

      for (BoxContainer::iterator itr = block_domain_boxes.begin();
           itr != block_domain_boxes.end(); ++itr) {
         Box box(*itr, local_id++, 0);
         box.setBlockId(BlockId(b));
         domain.pushBack(box);
      }

      if (d_number_blocks > 1) {

         std::string singularity_db_name =
            "Singularity_" + blk_string;
         std::shared_ptr<tbox::Database> singularity_db =
            db->getDatabase(singularity_db_name);
         d_singularity[b].getFromRestart(*singularity_db);

         std::string neighbors_db_name =
            "Neighbors_" + blk_string;
         std::shared_ptr<tbox::Database> neighbors_db =
            db->getDatabase(neighbors_db_name);
         int num_neighbors = neighbors_db->getInteger("num_neighbors");
         for (int count = 0; count < num_neighbors; ++count) {
            std::string neighbor_db_name =
               "neighbor_" + tbox::Utilities::intToString(count);
            std::shared_ptr<tbox::Database> neighbor_db =
               neighbors_db->getDatabase(neighbor_db_name);
            BlockId nbr_block_id(neighbor_db->getInteger("nbr_block_id"));
            BoxContainer nbr_transformed_domain;
            nbr_transformed_domain.getFromRestart(*neighbor_db);
            Transformation::RotationIdentifier nbr_rotation_ident =
               static_cast<Transformation::RotationIdentifier>(
                  neighbor_db->getInteger("rotation_identifier"));
            IntVector nbr_offset(dim);
            nbr_offset.getFromRestart(*neighbor_db, "offset");
            BlockId nbr_begin_block(neighbor_db->getInteger("begin_block"));
            BlockId nbr_end_block(neighbor_db->getInteger("end_block"));
            bool nbr_is_singularity = neighbor_db->getBool("d_is_singularity");
            Transformation nbr_transformation(nbr_rotation_ident,
               nbr_offset,
               nbr_begin_block,
               nbr_end_block);
            std::vector<Transformation> restart_transformation(
               1, nbr_transformation);
            Neighbor block_nbr(nbr_block_id,
               nbr_transformed_domain,
               restart_transformation);
            block_nbr.setSingularity(nbr_is_singularity);
            std::pair<BlockId, Neighbor> nbr_pair(nbr_block_id, block_nbr);
            d_block_neighbors[b].insert(nbr_pair);
         }
      }
   }
   setPhysicalDomain(domain, d_number_blocks);

   d_has_enhanced_connectivity = db->getInteger("d_has_enhanced_connectivity");

   d_number_of_block_singularities =
      db->getInteger("d_number_of_block_singularities");

   IntVector periodic_shift(dim);
   int* temp_shift = &periodic_shift[0];
   db->getIntegerArray("periodic_dimension", temp_shift, dim.getValue());
   initializePeriodicShift(periodic_shift);
}

/*
 *************************************************************************
 *
 * Data is read from input only if the simulation is not from restart.
 * Otherwise, all values specifed in the input database are ignored.
 * In this method data from the database are read to local
 * variables and the setPhysicalDomain() method is called.
 *
 *************************************************************************
 */

void
BaseGridGeometry::getFromInput(
   const std::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart,
   bool allow_multiblock)
{
   if (!is_from_restart && !input_db) {
      TBOX_ERROR(": BaseGridGeometry::getFromInput()\n"
         << "no input database supplied" << std::endl);
   }

   const tbox::Dimension dim(getDim());

   if (!is_from_restart) {

      d_number_blocks = 
         static_cast<BlockId::block_t>(
            input_db->getIntegerWithDefault("num_blocks", 1));
      if (!(d_number_blocks >= 1)) {
         INPUT_RANGE_ERROR("num_blocks");
      }

      if (d_number_blocks > 1 && !allow_multiblock) {
         TBOX_ERROR("BaseGridGeometry::getFromInput error...\n"
            << "num_blocks is >1 for an inherently single block grid geometry."
            << std::endl);
      }

      if (d_ratio_to_level_zero[0].getNumBlocks() != d_number_blocks) {
         d_ratio_to_level_zero[0] = 
            IntVector(IntVector::getOne(d_dim), d_number_blocks);
      }

      std::string domain_name;
      BoxContainer domain;
      LocalId local_id(0);

      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {

         domain_name = "domain_boxes_" +
            tbox::Utilities::intToString(static_cast<int>(b));

         BoxContainer block_domain_boxes;
         if (input_db->keyExists(domain_name)) {
            std::vector<tbox::DatabaseBox> db_box_vector =
               input_db->getDatabaseBoxVector(domain_name);
            block_domain_boxes = db_box_vector;
            if (block_domain_boxes.empty()) {
               TBOX_ERROR(
                  getObjectName() << ":  "
                                  << "No boxes for " << domain_name
                                  << " array found in input." << std::endl);
            }
         } else if (b == 0 && d_number_blocks == 1 &&
                    input_db->keyExists("domain_boxes")) {
            std::vector<tbox::DatabaseBox> db_box_vector =
               input_db->getDatabaseBoxVector("domain_boxes");
            block_domain_boxes = db_box_vector;
            if (block_domain_boxes.empty()) {
               TBOX_ERROR(
                  getObjectName() << ":  "
                                  << "No boxes for domain_boxes"
                                  << " array found in input." << std::endl);
            }
         } else {
            TBOX_ERROR(
               getObjectName() << ":  "
                               << "Key data '" << domain_name << "' not found in input."
                               << std::endl);
         }

         for (BoxContainer::iterator itr = block_domain_boxes.begin();
              itr != block_domain_boxes.end(); ++itr) {
            Box box(*itr, local_id++, 0);
            box.setBlockId(BlockId(b));
            domain.pushBack(box);
         }

      }

      int pbc[SAMRAI::MAX_DIM_VAL];
      IntVector per_bc(dim, 0);
      if (input_db->keyExists("periodic_dimension")) {
         input_db->getIntegerArray("periodic_dimension", pbc, dim.getValue());
         for (int i = 0; i < dim.getValue(); ++i) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }

      if (d_number_blocks > 1 && per_bc != IntVector::getZero(dim)) {
         TBOX_ERROR("BaseGridGeometry::getFromInput() error...\n"
            << "periodic boundaries are not currently supported for multiblock meshes."
            << std::endl);
      }

      initializePeriodicShift(per_bc);

      setPhysicalDomain(domain, d_number_blocks);

      readBlockDataFromInput(input_db);
   } else if (input_db) {
      bool read_on_restart =
         input_db->getBoolWithDefault("read_on_restart", false);
      int num_keys = static_cast<int>(input_db->getAllKeys().size());
      if (num_keys > 0 && read_on_restart) {
         TBOX_WARNING(
            "BaseGridGeometry::getFromInput() warning...\n"
            << "You want to override restart data with values from\n"
            << "an input database which is not allowed." << std::endl);
      }
   }
}

/*
 *************************************************************************
 *
 * Writes version number and data members for the class to restart database.
 *
 *************************************************************************
 */

void
BaseGridGeometry::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   const tbox::Dimension dim(getDim());

   restart_db->putInteger("HIER_GRID_GEOMETRY_VERSION",
      HIER_GRID_GEOMETRY_VERSION);

   restart_db->putInteger("num_blocks", static_cast<int>(d_number_blocks));

   for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {

      std::string blk_string =
         tbox::Utilities::intToString(static_cast<int>(b));
      std::string domain_name =
         "domain_boxes_" + blk_string;

      BoxContainer block_phys_domain(getPhysicalDomain(), BlockId(b));
      std::vector<tbox::DatabaseBox> temp_box_vector = block_phys_domain;

      restart_db->putDatabaseBoxVector(domain_name, temp_box_vector);

      if (d_number_blocks > 1) {

         std::string singularity_db_name =
            "Singularity_" + blk_string;
         std::shared_ptr<tbox::Database> singularity_db =
            restart_db->putDatabase(singularity_db_name);
         d_singularity[b].putToRestart(singularity_db);

         std::string neighbors_db_name =
            "Neighbors_" + blk_string;
         std::shared_ptr<tbox::Database> neighbors_db =
            restart_db->putDatabase(neighbors_db_name);
         neighbors_db->putInteger("num_neighbors",
            static_cast<int>(d_block_neighbors[b].size()));
         int count = 0;
         for (std::map<BlockId, Neighbor>::const_iterator ni = d_block_neighbors[b].begin();
              ni != d_block_neighbors[b].end(); ++ni) {
            const Neighbor& neighbor = ni->second;
            std::string neighbor_db_name =
               "neighbor_" +
               tbox::Utilities::intToString(static_cast<int>(count));
            std::shared_ptr<tbox::Database> neighbor_db =
               neighbors_db->putDatabase(neighbor_db_name);
            neighbor_db->putInteger("nbr_block_id",
               static_cast<int>(neighbor.getBlockId().getBlockValue()));
            neighbor.getTransformedDomain().putToRestart(neighbor_db);
            neighbor_db->putInteger("rotation_identifier",
               neighbor.getTransformation(0).getRotation());
            neighbor.getTransformation(0).getOffset().putToRestart(*neighbor_db,
               "offset");
            neighbor_db->putInteger("begin_block",
               static_cast<int>(
                  neighbor.getTransformation(0).getBeginBlock().getBlockValue()));
            neighbor_db->putInteger("end_block",
               static_cast<int>(
                  neighbor.getTransformation(0).getEndBlock().getBlockValue()));
            neighbor_db->putBool("d_is_singularity", neighbor.isSingularity());
            ++count;
         }
      }
   }

   restart_db->putInteger("d_has_enhanced_connectivity",
      d_has_enhanced_connectivity);

   restart_db->putInteger("d_number_of_block_singularities",
      d_number_of_block_singularities);

   IntVector level0_shift(getPeriodicShift(IntVector::getOne(dim)));
   int* temp_shift = &level0_shift[0];
   restart_db->putIntegerArray("periodic_dimension",
      temp_shift,
      dim.getValue());
}

/*
 * ************************************************************************
 *
 * Compute the valid periodic shifts for the given box.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::computeShiftsForBox(
   std::vector<IntVector>& shifts,
   const Box& box,
   const BoxContainer& domain_search_tree,
   const IntVector& periodic_shift) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, box, periodic_shift);

   shifts.clear();

   int num_periodic_dirs = 0;

   for (int i = 0; i < d_dim.getValue(); ++i) {
      if (periodic_shift(i) != 0) {
         ++num_periodic_dirs;
      }
   }

   if (num_periodic_dirs > 0) {

      shifts.reserve(d_periodic_shift_catalog.getNumberOfShifts());

      BoundaryLookupTable* blut =
         BoundaryLookupTable::getLookupTable(d_dim);

      const std::vector<int>& location_index_max =
         blut->getMaxLocationIndices();

      for (tbox::Dimension::dir_t d = 0; d < num_periodic_dirs; ++d) {

         const tbox::Dimension::dir_t codim = static_cast<tbox::Dimension::dir_t>(d + 1);

         for (int loc = 0; loc < location_index_max[d]; ++loc) {

            const std::vector<tbox::Dimension::dir_t>& dirs = blut->getDirections(loc, codim);

            bool need_to_test = true;
            for (int k = 0; k < static_cast<int>(dirs.size()); ++k) {
               if (periodic_shift(dirs[k]) == 0) {
                  need_to_test = false;
                  break;
               }
            }

            if (need_to_test) {

               Box border(box);
               IntVector border_shift(d_dim, 0);

               std::vector<bool> is_upper(codim);
               for (tbox::Dimension::dir_t j = 0; j < codim; ++j) {
                  if (blut->isUpper(loc, codim, j)) {
                     border.setLower(dirs[j], box.upper(dirs[j]));
                     border.setUpper(dirs[j], box.upper(dirs[j]));
                     border_shift(dirs[j]) = 1;
                     is_upper[j] = true;
                  } else {
                     border.setLower(dirs[j], box.lower(dirs[j]));
                     border.setUpper(dirs[j], box.lower(dirs[j]));
                     border_shift(dirs[j]) = -1;
                     is_upper[j] = false;
                  }
               }

               border.shift(border_shift);
               BoxContainer border_list(border);

               border_list.removeIntersections(domain_search_tree);

               if (!border_list.empty()) {

                  const Box& domain_bound_box =
                     domain_search_tree.getBoundingBox();

                  if (codim == 1) {

                     IntVector new_shift(d_dim, 0);
                     if (is_upper[0]) {
                        new_shift(dirs[0]) =
                           -domain_bound_box.numberCells(dirs[0]);
                     } else {
                        new_shift(dirs[0]) =
                           domain_bound_box.numberCells(dirs[0]);
                     }
                     // shifts.addItem(new_shift);
                     shifts.insert(shifts.end(), new_shift);

                  } else {

                     bool shift_to_add = true;
                     for (int c = 0; c < codim; ++c) {

                        if (is_upper[c]) {
                           if (border.upper(dirs[c]) <=
                               domain_bound_box.upper(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        } else {
                           if (border.lower(dirs[c]) >=
                               domain_bound_box.lower(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        }

                     }

                     if (shift_to_add) {
                        IntVector new_shift(d_dim, 0);
                        for (int b = 0; b < codim; ++b) {
                           if (is_upper[b]) {
                              new_shift(dirs[b]) =
                                 -domain_bound_box.numberCells(dirs[b]);
                           } else {
                              new_shift(dirs[b]) =
                                 domain_bound_box.numberCells(dirs[b]);
                           }
                        }
                        // shifts.addItem(new_shift);
                        shifts.insert(shifts.end(), new_shift);
                     }
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Decompose patch boundary region into pieces depending on spatial dim.
 *
 *************************************************************************
 */

void
BaseGridGeometry::getBoundaryBoxes(
   PatchBoundaries& patch_boundaries,
   const Box& box,
   const BoxContainer& domain_boxes,
   const IntVector& ghosts,
   const IntVector& periodic_shift) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(*this, box, ghosts, periodic_shift);

   t_get_boundary_boxes->start();

   const Index ifirst = box.lower();
   const Index ilast = box.upper();

   int num_per_dirs = 0;
   for (int d = 0; d < d_dim.getValue(); ++d) {
      num_per_dirs += (periodic_shift(d) ? 1 : 0);
   }

   if (num_per_dirs == d_dim.getValue()) {
      for (int k = 0; k < d_dim.getValue(); ++k) {
         patch_boundaries[k].clear();
      }
   } else {
      if (!domain_boxes.hasTree() && domain_boxes.size() > 10) {
         domain_boxes.makeTree(0);
      }

      BoxContainer per_domain_boxes;
      if (num_per_dirs != 0) {
         per_domain_boxes = domain_boxes;
         per_domain_boxes.grow(periodic_shift);
         if (per_domain_boxes.size() > 10) {
            per_domain_boxes.makeTree(0);
         }
      }

      BoundaryLookupTable* blut =
         BoundaryLookupTable::getLookupTable(d_dim);

      const std::vector<int>& location_index_max =
         blut->getMaxLocationIndices();
      std::vector<BoxContainer> codim_boxlist(d_dim.getValue());

      for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue() - num_per_dirs; ++d) {

         tbox::Dimension::dir_t codim = static_cast<tbox::Dimension::dir_t>(d + 1);

         for (int loc = 0; loc < location_index_max[d]; ++loc) {
            const std::vector<tbox::Dimension::dir_t>& dirs = blut->getDirections(loc, codim);

            bool all_is_per = true;
            for (tbox::Dimension::dir_t p = 0; p < codim; ++p) {
               if (periodic_shift(dirs[p]) == 0) {
                  all_is_per = false;
               }
            }

            if (!all_is_per) {
               Box border(box);
               IntVector border_shift(d_dim, 0);

               for (tbox::Dimension::dir_t i = 0; i < codim; ++i) {
                  if (blut->isUpper(loc, codim, i)) {
                     border.setLower(dirs[i], box.upper(dirs[i]));
                     border.setUpper(dirs[i], box.upper(dirs[i]));
                     border_shift(dirs[i]) = 1;
                  } else {
                     border.setLower(dirs[i], box.lower(dirs[i]));
                     border.setUpper(dirs[i], box.lower(dirs[i]));
                     border_shift(dirs[i]) = -1;
                  }
               }

               // grow in non-dirs directions
               for (tbox::Dimension::dir_t j = 0; j < d_dim.getValue(); ++j) {
                  bool dir_used = false;
                  for (tbox::Dimension::dir_t du = 0; du < codim; ++du) {
                     if (dirs[du] == j) {
                        dir_used = true;
                        break;
                     }
                  }
                  if (!dir_used) {
                     border.setUpper(j, ilast(j) + ghosts(j));
                     border.setLower(j, ifirst(j) - ghosts(j));
                  }
               }

               /*
                * Intersect border_list with domain, then shift so that
                * true boundary boxes are outside domain.  Then remove
                * intersections with the domain.
                */

               BoxContainer border_list(border);
               if (num_per_dirs != 0) {
                  border_list.intersectBoxes(per_domain_boxes);
               } else {
                  border_list.intersectBoxes(domain_boxes);
               }
               border_list.shift(border_shift);

               if (num_per_dirs != 0) {
                  border_list.removeIntersections(per_domain_boxes);
               } else {
                  border_list.removeIntersections(domain_boxes);
               }

               if (!border_list.empty()) {
                  for (int bd = 0; bd < d; ++bd) {
                     border_list.removeIntersections(codim_boxlist[bd]);

                     if (border_list.empty()) {
                        break;
                     }
                  }
               }

               if (!border_list.empty()) {
                  border_list.coalesce();
                  for (BoxContainer::iterator bl = border_list.begin();
                       bl != border_list.end(); ++bl) {

                     BoundaryBox boundary_box(*bl, codim, loc);

                     patch_boundaries[d].push_back(boundary_box);
                  }

                  codim_boxlist[d].spliceFront(border_list);
               }
            }

         }
      }
   }
   t_get_boundary_boxes->stop();
}

/*
 *************************************************************************
 *
 * Compute physical domain for index space related to reference domain
 * by specified ratio.  If any entry of ratio is negative, the reference
 * domain will be coarsened.  Otherwise, it will be refined.
 *
 *************************************************************************
 */

void
BaseGridGeometry::computePhysicalDomain(
   BoxContainer& domain_boxes,
   const IntVector& ratio_to_level_zero,
   const BlockId& block_id) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio_to_level_zero);
   const BlockId::block_t b = block_id.getBlockValue();

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   if (ratio_to_level_zero != 1) {
      int i;
      for (i = 0; i < d_dim.getValue(); ++i) {
         TBOX_ASSERT(ratio_to_level_zero(b,i) != 0);
      }
      if (d_dim.getValue() > 1) {
         for (i = 0; i < d_dim.getValue(); ++i) {
            bool pos0 = ratio_to_level_zero(b,i) > 0;
            bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
            TBOX_ASSERT(pos0 == pos1
               || (ratio_to_level_zero(b,i) == 1)
               || ((ratio_to_level_zero(b,i + 1) % d_dim.getValue()) == 1));
         }
      }
   }
#endif

   domain_boxes.clear();
   for (BoxContainer::const_iterator itr = d_physical_domain.begin();
        itr != d_physical_domain.end(); ++itr) {
      if (itr->getBlockId() == block_id) {
         domain_boxes.insert(*itr);
      }
   }

   if (ratio_to_level_zero != 1) {
      bool coarsen = false;
      IntVector tmp_rat = ratio_to_level_zero;
      for (int id = 0; id < d_dim.getValue(); ++id) {
         if (ratio_to_level_zero(b,id) < 0) coarsen = true;
         tmp_rat(b,id) = abs(ratio_to_level_zero(b,id));
      }
      if (coarsen) {
         domain_boxes.coarsen(tmp_rat);
      } else {
         domain_boxes.refine(tmp_rat);
      }
   }
}

/*
 *************************************************************************
 *
 * Compute physical domain for index space related to reference domain
 * by specified ratio.  If any entry of ratio is negative, the reference
 * domain will be coarsened.  Otherwise, it will be refined.
 *
 *************************************************************************
 */

void
BaseGridGeometry::computePhysicalDomain(
   BoxLevel& box_level,
   const IntVector& ratio_to_level_zero,
   const BlockId& block_id) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio_to_level_zero);
   const BlockId::block_t b = block_id.getBlockValue();

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   if (ratio_to_level_zero != 1) {
      int i;
      for (i = 0; i < d_dim.getValue(); ++i) {
         TBOX_ASSERT(ratio_to_level_zero(b,i) != 0);
      }
      if (d_dim.getValue() > 1) {
         for (i = 0; i < d_dim.getValue(); ++i) {
            bool pos0 = ratio_to_level_zero(b,i) > 0;
            bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
            TBOX_ASSERT(pos0 == pos1
               || (ratio_to_level_zero(b,i) == 1)
               || (ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
         }
      }
   }
#endif

   for (BoxContainer::const_iterator itr = d_domain_with_images.begin();
        itr != d_domain_with_images.end();
        ++itr) {
      if (itr->getBlockId() == block_id) {
         box_level.addBoxWithoutUpdate(*itr);
      }
   }

   if (ratio_to_level_zero != 1) {
      bool coarsen = false;
      IntVector tmp_rat = ratio_to_level_zero;
      for (int id = 0; id < d_dim.getValue(); ++id) {
         if (ratio_to_level_zero(b,id) < 0) {
            coarsen = true;
         }
         tmp_rat(id) = abs(ratio_to_level_zero(b,id));
      }
      if (coarsen) {
         box_level.coarsenBoxes(box_level, tmp_rat, IntVector::getOne(d_dim));
      } else {
         box_level.refineBoxes(box_level, tmp_rat, IntVector::getOne(d_dim));
      }
   }
}

/*
 *************************************************************************
 *
 * Compute physical domain for index space related to reference domain
 * by specified ratio.  If any entry of ratio is negative, the reference
 * domain will be coarsened.  Otherwise, it will be refined.
 *
 *************************************************************************
 */

void
BaseGridGeometry::computePhysicalDomain(
   BoxContainer& domain_boxes,
   const IntVector& ratio_to_level_zero) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   if (ratio_to_level_zero != 1) {
      /*
       * All components of ratio must be nonzero.  Additionally, all components
       * of ratio not equal to 1 must have the same sign.
       */
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         int i;
         for (i = 0; i < d_dim.getValue(); ++i) {
            TBOX_ASSERT(ratio_to_level_zero(b,i) != 0);
         }
         if (d_dim.getValue() > 1) {
            for (i = 0; i < d_dim.getValue(); ++i) {
               bool pos0 = ratio_to_level_zero(b,i) > 0;
               bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
               TBOX_ASSERT(pos0 == pos1
                  || (ratio_to_level_zero(b,i) == 1)
                  || (ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
            }
         }
      }
   } 
#endif

   domain_boxes = d_domain_with_images;

   if (ratio_to_level_zero != 1) {
      bool coarsen = false;
      IntVector tmp_rat = ratio_to_level_zero;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int id = 0; id < d_dim.getValue(); ++id) {
            if (ratio_to_level_zero(b,id) < 0) {
               coarsen = true;
            }
            tmp_rat(b,id) = abs(ratio_to_level_zero(b,id));
         }
      }
      if (coarsen) {
         domain_boxes.coarsen(tmp_rat);
      } else {
         domain_boxes.refine(tmp_rat);
      }
   }

}

void
BaseGridGeometry::computePhysicalDomain(
   BoxLevel& box_level,
   const IntVector& ratio_to_level_zero) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   if (ratio_to_level_zero != 1) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         int i;
         for (i = 0; i < d_dim.getValue(); ++i) {
            TBOX_ASSERT(ratio_to_level_zero(b,i) != 0);
         }
         if (d_dim.getValue() > 1) {
            for (i = 0; i < d_dim.getValue(); ++i) {
               bool pos0 = ratio_to_level_zero(b,i) > 0;
               bool pos1 = ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
               TBOX_ASSERT(pos0 == pos1
                  || (ratio_to_level_zero(b,i) == 1)
                  || (ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
            }
         }
      }
   }
#endif

   BoxContainer domain_boxes = d_domain_with_images;

   if (ratio_to_level_zero != 1) {
      bool coarsen = false;
      IntVector tmp_rat = ratio_to_level_zero;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int id = 0; id < d_dim.getValue(); ++id) {
            if (ratio_to_level_zero(b,id) < 0) {
               coarsen = true;
            }
            tmp_rat(b,id) = abs(ratio_to_level_zero(b,id));
         }
      }
      if (coarsen) {
         domain_boxes.coarsen(tmp_rat);
      } else {
         domain_boxes.refine(tmp_rat);
      }
   }

   for (BoxContainer::const_iterator bi = domain_boxes.begin();
        bi != domain_boxes.end(); ++bi) {

      box_level.addBoxWithoutUpdate(*bi);

   }
}

/*
 *************************************************************************
 *
 * Set physical domain data member from input box array and determine
 * whether domain is a single box.
 *
 *************************************************************************
 */

void
BaseGridGeometry::setPhysicalDomain(
   const BoxContainer& domain,
   const size_t number_blocks)
{
   TBOX_ASSERT(!domain.empty());
#ifdef DEBUG_CHECK_ASSERTIONS
   for (BoxContainer::const_iterator itr = domain.begin(); itr != domain.end();
        ++itr) {
      TBOX_ASSERT(itr->getBlockId().isValid());
      TBOX_ASSERT(itr->getBlockId().getBlockValue() < number_blocks);
   }
#endif

   d_physical_domain.clear();

   d_domain_is_single_box.resize(number_blocks);
   d_number_blocks = number_blocks;
   if (d_ratio_to_level_zero[0].getNumBlocks() != d_number_blocks) {
      d_ratio_to_level_zero[0] = 
         IntVector(IntVector::getOne(d_dim), d_number_blocks);
   }

   LocalId local_id(0);

   for (BlockId::block_t b = 0; b < number_blocks; ++b) {
      BlockId block_id(b);

      BoxContainer block_domain(domain, block_id);
      TBOX_ASSERT(!block_domain.empty());
      Box bounding_box(block_domain.getBoundingBox());
      BoxContainer bounding_cntnr(bounding_box);
      bounding_cntnr.removeIntersections(block_domain);
      if (bounding_cntnr.empty()) {
         d_domain_is_single_box[b] = true;
         Box box(bounding_box, local_id++, 0);
         box.setBlockId(block_id);
         d_physical_domain.pushBack(box);
      } else {
         d_domain_is_single_box[b] = false;
         for (BoxContainer::iterator itr = block_domain.begin();
              itr != block_domain.end(); ++itr) {
            d_physical_domain.pushBack(*itr);
         }
      }
   }

   if (d_physical_domain.size() == 1 &&
       d_periodic_shift != IntVector::getZero(d_dim)) {

      /*
       * If necessary, reset periodic shift amounts using the new
       * bounding box.
       */
      for (int id = 0; id < d_dim.getValue(); ++id) {
         d_periodic_shift(id) = ((d_periodic_shift(id) == 0) ? 0 : 1);
      }

      if (d_periodic_shift != IntVector::getZero(d_dim)) {
         /*
          * Check if the physical domain is valid for the specified
          * periodic conditions.  If so, compute the shift in each
          * direction based on the the number of cells.
          */
         if (checkPeriodicValidity(d_physical_domain)) {

            Box bounding_box(d_physical_domain.getBoundingBox());

            for (tbox::Dimension::dir_t id = 0; id < d_dim.getValue(); ++id) {
               d_periodic_shift(id) *= bounding_box.numberCells(id);
            }

         } else {
            TBOX_ERROR("Error in BaseGridGeometry object with name = "
               << d_object_name << ": in initializePeriodicShift():  "
               << "Domain is not periodic for one (or more) of the directions "
               << "specified in the geometry input file!" << std::endl);
         }
      }
   }

   resetDomainBoxContainer();

}

/*
 *************************************************************************
 *
 * Reset the domain BoxContainer based on current definition of
 * physical domain and periodic shift.
 *
 *************************************************************************
 */

void
BaseGridGeometry::resetDomainBoxContainer()
{
   d_physical_domain.makeTree(this);

   const bool is_periodic =
      d_periodic_shift != IntVector::getZero(d_periodic_shift.getDim());

   d_domain_with_images = d_physical_domain; // Images added next if is_periodic.

   if (is_periodic) {

      d_periodic_shift_catalog.initializeShiftsByIndexDirections(d_periodic_shift);

      const IntVector& one_vector(IntVector::getOne(d_dim));

      for (BoxContainer::const_iterator ni = d_physical_domain.begin();
           ni != d_physical_domain.end(); ++ni) {

         const Box& real_box = *ni;
         TBOX_ASSERT(real_box.getPeriodicId() == d_periodic_shift_catalog.getZeroShiftNumber());

         for (int ishift = 1; ishift < d_periodic_shift_catalog.getNumberOfShifts();
              ++ishift) {

            d_domain_with_images.emplaceBack(real_box,
                                             PeriodicId(ishift),
                                             one_vector,
                                             d_periodic_shift_catalog);
         }

      }
   }
   d_domain_with_images.makeTree(this);

}

/*
 *************************************************************************
 *
 * The argument is an IntVector of length DIM.  It is set to 1
 * for periodic directions and 0 for all other directions.  In the
 * periodic directions, the coarse-level shift is calculated and stored
 * in the IntVector d_periodic_shift. The shift is the number of cells
 * in each periodic direction and is zero in all other directions.
 *
 *************************************************************************
 */

void
BaseGridGeometry::initializePeriodicShift(
   const IntVector& directions)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, directions);

   d_periodic_shift = directions;

   if (d_physical_domain.size() == 1) {
      resetDomainBoxContainer();
   }
}

/*
 *************************************************************************
 *
 * This returns an IntVector of length d_dim that is set to the width of
 * the domain in periodic directions and 0 in all other directions.
 * the argument contains the refinement ratio relative to the coarsest
 * level, which is multiplied by d_periodic_shift to get the return
 * vector.
 *
 *************************************************************************
 */

IntVector
BaseGridGeometry::getPeriodicShift(
   const IntVector& ratio_to_level_zero) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ratio_to_level_zero);

   IntVector periodic_shift(d_dim);
   if (d_number_blocks > 1) {
      if (d_periodic_shift != IntVector::getZero(d_dim)) {
         TBOX_ERROR("BaseGridGeometry::getPeriodicShift() error...\n"
            << "A non-zero periodic shift cannont be used in a multiblock geometry."
            << std::endl);
      } else {
         periodic_shift = d_periodic_shift;
      }
   } else {
      /*
       * This is single-block, so only use zero BlockId.
       */
      const IntVector& block_ratio =
         ratio_to_level_zero;

#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * All components of ratio vector must be nonzero.  Additionally,
       * all components not equal to 1 must have the same sign.
       */
      int k;
      for (k = 0; k < d_dim.getValue(); ++k) {
         TBOX_ASSERT(block_ratio(k) != 0);
      }
      if (d_dim.getValue() > 1) {
         for (k = 0; k < d_dim.getValue(); ++k) {
            TBOX_ASSERT((block_ratio(k)
                         * block_ratio((k + 1) % d_dim.getValue()) > 0)
               || (block_ratio(k) == 1)
               || (block_ratio((k + 1) % d_dim.getValue()) == 1));
         }
      }
#endif

      for (int i = 0; i < d_dim.getValue(); ++i) {
         if (block_ratio(i) > 0) {
            periodic_shift(i) = d_periodic_shift(i) * block_ratio(i);
         } else {
            int abs_ratio = abs(block_ratio(i));
            periodic_shift(i) = d_periodic_shift(i) / abs_ratio;
         }
      }
   }
   return periodic_shift;
}

/*
 *************************************************************************
 *
 * This checks if the periodic directions given to the constructor are
 * valid for the domain.  Periodic directions are valid if the domain
 * has exactly two physical boundaries normal to the periodic direction.
 *
 *************************************************************************
 */

bool
BaseGridGeometry::checkPeriodicValidity(
   const BoxContainer& domain)
{
   bool is_valid = true;

   IntVector valid_direction(d_dim, 1);
   IntVector grow_direction(d_dim, 1);

   /*
    * Compute the bounding box of a "duplicate" domain + 1
    * cell and set the min and max indices of this grown box.
    */
   BoxContainer dup_domain(domain);

   Box domain_box = dup_domain.getBoundingBox();
   domain_box.grow(grow_direction);
   tbox::Dimension::dir_t i;
   Index min_index(d_dim, 0), max_index(d_dim, 0);
   for (i = 0; i < d_dim.getValue(); ++i) {
      //set min/max of the bounding box
      min_index(i) = domain_box.lower(i);
      max_index(i) = domain_box.upper(i);
   }

   /*
    * Next, for each direction, grow another "duplicate" domain
    * by 1.  Remove the intersections with the original domain,
    * and loop through the remaining box list, checking if the
    * upper index of the box matches the bounding box max or the
    * lower index of the box matches the bounding box min.  If
    * not, this direction is not a valid periodic direction.
    */
   for (i = 0; i < d_dim.getValue(); ++i) {
      BoxContainer dup_domain2(domain);
      IntVector grow_one(d_dim, 0);
      grow_one(i) = 1;
      dup_domain2.grow(grow_one);
      dup_domain2.unorder();
      dup_domain2.removeIntersections(domain);

      BoxContainer::iterator n = dup_domain2.begin();
      for ( ; n != dup_domain2.end(); ++n) {
         Box this_box = *n;
         Index box_lower = this_box.lower();
         Index box_upper = this_box.upper();
         if (d_periodic_shift(i) != 0) {
            if (!((box_lower(i) == min_index(i)) ||
                  (box_upper(i) == max_index(i)))) {
               valid_direction(i) = 0;
            }
         }
      }
   }

   for (i = 0; i < d_dim.getValue(); ++i) {
      if ((valid_direction(i) == 0) &&
          (d_periodic_shift(i) != 0)) {
         is_valid = false;
      }
   }

   return is_valid;
}

/*
 *************************************************************************
 *
 * Perform an error check on a recently-constructed boundary box to
 * make sure that it is the proper size, is adjacent to a patch, and is
 * outside the physical domain.
 *
 *************************************************************************
 */

bool
BaseGridGeometry::checkBoundaryBox(
   const BoundaryBox& boundary_box,
   const Patch& patch,
   const BoxContainer& domain,
   const int num_per_dirs,
   const IntVector& max_data_ghost_width) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(*this,
      boundary_box,
      patch,
      max_data_ghost_width);

   bool return_val = true;

   const Box& bbox = boundary_box.getBox();

   /*
    * Test to see that the box is of size 1 in at least 1 direction.
    */
   IntVector box_size(d_dim);

   for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
      box_size(i) = bbox.numberCells(i);
   }

   if (box_size.min() != 1) {
      return_val = false;
   }

   /*
    * Quick and dirty test to see that boundary box is adjacent to patch
    * boundary, or a patch boundary extended through the ghost region.
    */
   Box patch_box = patch.getBox();

   Box grow_patch_box(patch_box);

   grow_patch_box.grow(IntVector::getOne(d_dim));

   if (!grow_patch_box.isSpatiallyEqual((grow_patch_box + bbox))) {
      bool valid_box = false;
      grow_patch_box = patch_box;
      for (tbox::Dimension::dir_t j = 0; j < d_dim.getValue(); ++j) {
         if (num_per_dirs == 0) {

            for (tbox::Dimension::dir_t k = 1; k < d_dim.getValue(); ++k) {

               grow_patch_box.grow(static_cast<tbox::Dimension::dir_t>((j + k) % d_dim.getValue()),
                  max_data_ghost_width((j + k) % d_dim.getValue()));

            }

         } else {

            for (tbox::Dimension::dir_t k = 1; k < d_dim.getValue(); ++k) {

               grow_patch_box.grow(static_cast<tbox::Dimension::dir_t>((j + k) % d_dim.getValue()),
                  2 * max_data_ghost_width((j + k) % d_dim.getValue()));

            }

         }
         grow_patch_box.grow(j, 1);
         if (grow_patch_box.isSpatiallyEqual((grow_patch_box + bbox))) {
            valid_box = true;
         }
         grow_patch_box = patch_box;
      }
      if (!valid_box) {
         return_val = false;
      }
   }

   /*
    * check that the boundary box is outside the physical domain.
    */
   BoxContainer bbox_list(bbox);
   bbox_list.intersectBoxes(domain);

   if (!bbox_list.empty()) {
      return_val = false;
   }

   return return_val;
}

/*
 ***************************************************************************
 * Read multiblock metadata from input database
 ***************************************************************************
 */
void
BaseGridGeometry::readBlockDataFromInput(
   const std::shared_ptr<tbox::Database>& input_db)
{
   TBOX_ASSERT(input_db);

   d_singularity.resize(d_number_blocks);
   d_block_neighbors.resize(d_number_blocks);

   std::string sing_name;
   std::string neighbor_name;

   for (int bn = 0; true; ++bn) {
      neighbor_name = "BlockNeighbors" + tbox::Utilities::intToString(bn);

      if (!input_db->keyExists(neighbor_name)) {
         break;
      }
      std::shared_ptr<tbox::Database> pair_db(
         input_db->getDatabase(neighbor_name));

      BlockId block_a(pair_db->getInteger("block_a"));
      BlockId block_b(pair_db->getInteger("block_b"));
      Transformation::RotationIdentifier rotation_b_to_a;

      IntVector shift(d_dim, 0);
      if (d_dim.getValue() == 1) {
         rotation_b_to_a = Transformation::NO_ROTATE;
      } else {
         std::vector<std::string> rstr =
            pair_db->getStringVector("rotation_b_to_a");
         rotation_b_to_a = Transformation::getRotationIdentifier(rstr, d_dim);

         std::vector<int> b_array =
            pair_db->getIntegerVector("point_in_b_space");
         std::vector<int> a_array =
            pair_db->getIntegerVector("point_in_a_space");

         Index b_index(d_dim);
         Index a_index(d_dim);

         for (int p = 0; p < d_dim.getValue(); ++p) {
            b_index(p) = b_array[p];
            a_index(p) = a_array[p];
         }

         Box b_box(b_index, b_index, block_b);
         Box a_box(a_index, a_index, block_a);

         b_box.rotate(rotation_b_to_a);
         Index b_rotated_point(b_box.lower());
         Index a_point = (a_box.lower());

         shift = a_point - b_rotated_point;
      }

      registerNeighbors(block_a, block_b,
         rotation_b_to_a, shift);

   }

   /*
    * Each singularity exists where a certain set of blocks touch.
    * Each element of this vector represents a singularity, and the
    * set of integers are the block numbers for the blocks touching that
    * singularity.
    */
   std::set<std::set<BlockId> > singularity_blocks;
   if (d_number_blocks > 1) {
      findSingularities(singularity_blocks);
   }
   d_number_of_block_singularities =
      static_cast<int>(singularity_blocks.size());

   if (d_number_blocks == 1 && d_number_of_block_singularities > 0) {
      TBOX_ERROR("BaseGridGeometry::readBlockDataFromInput() error...\n"
         << "block singularities specified for single block problem."
         << std::endl);
   }

   if (d_number_of_block_singularities > 0) {

      /*
       * Process the singularites to determine if they are enhanced or
       * reduced connectivity, then compute and store needed internal data
       * for each case.
       */
      for (std::set<std::set<BlockId> >::iterator
           si = singularity_blocks.begin();
           si != singularity_blocks.end(); ++si) {

         for (std::set<BlockId>::iterator sbi = si->begin();
              sbi != si->end(); ++sbi) {

            const BlockId& cur_block_id = *sbi;
            const BlockId::block_t& cur_block = cur_block_id.getBlockValue();
            BoxContainer cur_grow(d_physical_domain, cur_block_id);
            cur_grow.unorder();
            cur_grow.grow(IntVector::getOne(d_dim));
            cur_grow.simplify();

            std::map<BlockId, Neighbor>& nbr_map =
               d_block_neighbors[cur_block];

            /*
             * nbr_ghost_buffer will contain buffers of width 1 covering the
             * space immediately across the block boundary with each
             * neighboring block.  The Boxes representing these buffers
             * are all transformed to the index spce of the current block.
             */
            std::map<BlockId, BoxContainer> nbr_ghost_buffer;

            for (std::map<BlockId, Neighbor>::iterator nei = nbr_map.begin();
                 nei != nbr_map.end(); ++nei) {

               const BlockId& nbr_blk = nei->second.getBlockId();

               if (si->find(nbr_blk) != si->end()) {

                  BoxContainer transformed_domain(
                     nei->second.getTransformedDomain());
                  transformed_domain.unorder();
                  transformed_domain.intersectBoxes(cur_grow);

                  nbr_ghost_buffer[nbr_blk] = transformed_domain;
               }
            }

            /*
             * Compare all the buffers in nbr_ghost_buffer to see if
             * any cover the same index space from the perspective of the
             * current block.  If any do, that means the current block has
             * more than one block neighbor that can fill the same ghost
             * region.  Those neighbors are enhanced connectivity neighbors
             * at this singularity.
             */
            std::set<BlockId::block_t> encon_nbrs;
            for (std::map<BlockId, BoxContainer>::iterator ng_itr =
                    nbr_ghost_buffer.begin();
                 ng_itr != nbr_ghost_buffer.end(); ++ng_itr) {

               const BoxContainer& ghost_buf = ng_itr->second;

               for (std::map<BlockId, BoxContainer>::iterator other = ng_itr;
                    other != nbr_ghost_buffer.end(); ++other) {

                  if (other != ng_itr) {

                     BoxContainer test_boxes(other->second);
                     test_boxes.intersectBoxes(ghost_buf);
                     if (!test_boxes.empty()) {
                        test_boxes.coalesce();
                        d_singularity[cur_block].spliceFront(test_boxes);
                        encon_nbrs.insert(ng_itr->first.getBlockValue());
                        encon_nbrs.insert(other->first.getBlockValue());
                     }
                  }
               }
            }

            /*
             * If neighboring blocks have been identified as enhanced
             * connectivity neighbors, set the flag in the Neighbor objects.
             */
            if (!encon_nbrs.empty()) {

               d_has_enhanced_connectivity = true;
               d_singularity[cur_block].coalesce();

               for (std::map<BlockId, Neighbor>::iterator nei = nbr_map.begin();
                    nei != nbr_map.end(); ++nei) {

                  const BlockId& nbr_blk = nei->second.getBlockId();
                  if (encon_nbrs.find(nbr_blk.getBlockValue()) !=
                      encon_nbrs.end()) {

                     nei->second.setSingularity(true);

                  }
               }

            } else {

               /*
                * If no enhanced connectivity neighbors have been found, we
                * must be at a reduced connectivity singularity.
                */

               /*
                * The location of the reduced connectivity singularity is
                * found by converting the block domain boxes to node-centered
                * representations and finding the intersection between
                * the blocks in singularity_blocks[si].
                */

               BoxContainer cur_domain(d_physical_domain, cur_block_id);
               cur_domain.coalesce();

               BoxContainer cur_domain_nodal(cur_domain);
               for (BoxContainer::iterator cdn = cur_domain_nodal.begin();
                    cdn != cur_domain_nodal.end(); ++cdn) {
                  cdn->setUpper(cdn->upper() + IntVector::getOne(d_dim));
               }

               for (std::map<BlockId, Neighbor>::iterator nei = nbr_map.begin();
                    nei != nbr_map.end(); ++nei) {

                  const BlockId& nbr_blk = nei->second.getBlockId();

                  if (si->find(nbr_blk) != si->end()) {

                     BoxContainer nbr_block_nodal(
                        nei->second.getTransformedDomain());
                     for (BoxContainer::iterator nbn = nbr_block_nodal.begin();
                          nbn != nbr_block_nodal.end(); ++nbn) {
                        nbn->setUpper(nbn->upper() + IntVector::getOne(d_dim));
                     }

                     cur_domain_nodal.intersectBoxes(nbr_block_nodal);
                  }
               }
               cur_domain_nodal.coalesce();
               TBOX_ASSERT(cur_domain_nodal.size() <= 1);

               /*
                * cur_domain_nodal now contains a node-centered box that
                * represents the location of the singularity.  Here we
                * convert that box to a cell-centered box located immediately
                * outside the current block's domain.
                */
               if (!cur_domain_nodal.empty()) {

                  const Box& sing_node_box = *(cur_domain_nodal.begin());
                  Box sing_box(d_dim);
                  sing_box.setBlockId(cur_block_id);

                  size_t sing_size = sing_node_box.size();

                  if (sing_size == 1) {

                     /*
                      * single point singularity.
                      */
                     const Index& sing_node = sing_node_box.lower();

                     std::unique_ptr<bool[]> found_corner =
                        std::make_unique<bool[]>(d_dim.getValue());

                     for (int d = 0; d < d_dim.getValue(); ++d) {
                        found_corner[d] = false;
                     }

                     bool use_box = true;

                     for (BoxContainer::iterator cd = cur_domain.begin();
                          cd != cur_domain.end(); ++cd) {

                        const Box& domain_box = *cd;

                        for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {

                           if (sing_node(d) == domain_box.lower(d)) {
                              sing_box.setLower(d, sing_node(d) - 1);
                              sing_box.setUpper(d, sing_node(d) - 1);
                              found_corner[d] = true;
                           } else if (sing_node(d) == domain_box.upper(d) + 1) {
                              sing_box.setLower(d, sing_node(d));
                              sing_box.setUpper(d, sing_node(d));
                              found_corner[d] = true;
                           }
                        }

                        for (int d = 0; d < d_dim.getValue(); ++d) {
                           if (!found_corner[d]) {
                              use_box = false;
                           }
                        }

                        if (use_box) {
                           d_singularity[cur_block].pushFront(sing_box);
                           break;
                        }
                     }
                  } else {
                     /*
                      * Line singularities
                      */
                     TBOX_ASSERT(d_dim.getValue() == 3);

                     IntVector width(sing_node_box.numberCells());
                     int num_width_one = 0;
                     int long_dir = -1;
                     for (int d = 0; d < d_dim.getValue(); ++d) {
                        if (width[d] == 1) {
                           ++num_width_one;
                        } else {
                           long_dir = d;
                        }
                     }
                     TBOX_ASSERT(long_dir >= 0);

                     if (num_width_one != 2) {
                        TBOX_ERROR("BaseGridGeometry::readBlockDataFromInput error...\n"
                           << "  object name = " << d_object_name
                           << " The computed singularity boundary for "
                           << sing_name << " is neither a point nor a line.  "
                           << " The user should verify that the input for "
                           << " the blocks at this singulary is correct."
                           << std::endl);
                        TBOX_ERROR("The computed singularity boundary between ");
                     }

                     const Index& sing_node = sing_node_box.lower();

                     std::unique_ptr<bool[]> found_corner =
                        std::make_unique<bool[]>(d_dim.getValue());

                     for (int d = 0; d < d_dim.getValue(); ++d) {
                        found_corner[d] = false;
                     }

                     bool use_box = true;

                     for (BoxContainer::iterator cd = cur_domain.begin();
                          cd != cur_domain.end(); ++cd) {

                        const Box& domain_box = *cd;

                        for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {

                           if (d != long_dir) {
                              if (sing_node(d) == domain_box.lower(d)) {
                                 sing_box.setLower(d, sing_node(d) - 1);
                                 sing_box.setUpper(d, sing_node(d) - 1);
                                 found_corner[d] = true;
                              } else if (sing_node(d) == domain_box.upper(d) + 1) {
                                 sing_box.setLower(d, sing_node(d));
                                 sing_box.setUpper(d, sing_node(d));
                                 found_corner[d] = true;
                              }
                           } else {
                              sing_box.setLower(d, sing_node_box.lower() (d));
                              sing_box.setUpper(d, sing_node_box.upper() (d) - 1);
                           }
                        }

                        for (int d = 0; d < d_dim.getValue(); ++d) {
                           if (d != long_dir && !found_corner[d]) {
                              use_box = false;
                           }
                        }
                        if (use_box) {
                           d_singularity[cur_block].pushFront(sing_box);
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

/*
 * ************************************************************************
 *
 * Get a BoxContainer representing all of the domain outside the given block.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::getDomainOutsideBlock(
   BoxContainer& domain_outside_block,
   const BlockId& block_id) const
{
   const std::map<BlockId, Neighbor>& nbr_map =
      d_block_neighbors[block_id.getBlockValue()];
   for (std::map<BlockId, Neighbor>::const_iterator nei = nbr_map.begin();
        nei != nbr_map.end(); ++nei) {
      BoxContainer transformed_domain(nei->second.getTransformedDomain());
      domain_outside_block.spliceFront(transformed_domain);
   }
}

/*
 * ************************************************************************
 *
 * Register a neighbor relationship between two blocks.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::registerNeighbors(
   const BlockId& block_a,
   const BlockId& block_b,
   const Transformation::RotationIdentifier rotation,
   const IntVector& shift_b_to_a)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, shift_b_to_a);

   const BlockId::block_t& a = block_a.getBlockValue();
   const BlockId::block_t& b = block_b.getBlockValue();
   BoxContainer b_domain_in_a_space(d_physical_domain, block_b);
   BoxContainer a_domain_in_b_space(d_physical_domain, block_a);
   b_domain_in_a_space.unorder();
   a_domain_in_b_space.unorder();

   Transformation::RotationIdentifier back_rotation =
      Transformation::getReverseRotationIdentifier(rotation, d_dim);
   IntVector back_shift(d_dim);

   if (d_dim.getValue() == 2 || d_dim.getValue() == 3) {
      Transformation::calculateReverseShift(back_shift,
         shift_b_to_a,
         rotation);
   } else {
      TBOX_ERROR("BaseGridGeometry::registerNeighbors error...\n"
         << "  object name = " << d_object_name
         << " Multiblock only works for 2D and 3D" << std::endl);
   }

   bool rotation_needed;
   if (rotation != 0) {
      rotation_needed = true;
   } else {
      rotation_needed = false;
   }

   if (rotation_needed) {
      b_domain_in_a_space.rotate(rotation);
      a_domain_in_b_space.rotate(back_rotation);
   }
   b_domain_in_a_space.shift(shift_b_to_a);
   a_domain_in_b_space.shift(back_shift);

   for (BoxContainer::iterator itr = b_domain_in_a_space.begin();
        itr != b_domain_in_a_space.end(); ++itr) {
      itr->setBlockId(block_a);
   }
   for (BoxContainer::iterator itr = a_domain_in_b_space.begin();
        itr != a_domain_in_b_space.end(); ++itr) {
      itr->setBlockId(block_b);
   }

   std::vector<Transformation> transformation;
   transformation.push_back(Transformation(rotation,
                                           shift_b_to_a,
                                           block_b,
                                           block_a));
   std::vector<Transformation> back_transformation;
   back_transformation.push_back(Transformation(back_rotation,
                                                back_shift,
                                                block_a,
                                                block_b));

   Neighbor neighbor_of_b(block_a, a_domain_in_b_space,
                          back_transformation);
   Neighbor neighbor_of_a(block_b, b_domain_in_a_space,
                          transformation);

   std::pair<BlockId, Neighbor> nbr_of_a_pair(block_b, neighbor_of_a);
   d_block_neighbors[a].insert(nbr_of_a_pair);
   std::pair<BlockId, Neighbor> nbr_of_b_pair(block_a, neighbor_of_b);
   d_block_neighbors[b].insert(nbr_of_b_pair);
}

/*
 * ************************************************************************
 *
 * Find singularities
 *
 * ************************************************************************
 */

void BaseGridGeometry::findSingularities(
   std::set<std::set<BlockId> >& singularity_blocks)
{
   TBOX_ASSERT(d_number_blocks > 1);
   TBOX_ASSERT(singularity_blocks.empty());

   BoxContainer chopped_domain;
   chopDomain(chopped_domain);

   chopped_domain.makeTree(this);

   std::map<BoxId, std::map<BoxId, int> > face_neighbors;

   for (BoxContainer::iterator b_itr = chopped_domain.begin();
        b_itr != chopped_domain.end(); ++b_itr) {

      const Box& base_box = *b_itr;
      const BlockId& base_block = base_box.getBlockId();
      const BoxId& base_id = base_box.getBoxId();

      Box grow_base(base_box);
      grow_base.grow(IntVector::getOne(d_dim));

      std::vector<const Box *> nbr_boxes;
      chopped_domain.findOverlapBoxes(nbr_boxes, grow_base, IntVector::getOne(d_dim));

      const std::map<BlockId, Neighbor>& nbrs_of_base =
         d_block_neighbors[base_block.getBlockValue()];

      for (std::vector<const Box *>::const_iterator n_itr = nbr_boxes.begin();
           n_itr != nbr_boxes.end(); ++n_itr) {

         const Box& nbr_box = **n_itr;
         const BoxId& nbr_id = nbr_box.getBoxId();

         if (nbr_id <= base_id) {
            continue;
         }

         const BlockId& nbr_block = nbr_box.getBlockId();
         if (base_block != nbr_block &&
             nbrs_of_base.find(nbr_block) == nbrs_of_base.end()) {
            continue;
         }

         if (face_neighbors[base_id].find(nbr_id) !=
             face_neighbors[base_id].end()) {
            TBOX_ASSERT(face_neighbors[nbr_id].find(base_id) !=
               face_neighbors[nbr_id].end());

            continue;
         }

         Box base_node_box(base_box);
         base_node_box.setUpper(
            base_node_box.upper() + IntVector::getOne(d_dim));

         Box transformed_nbr_box(nbr_box);
         if (nbr_block != base_block) {
            transformBox(transformed_nbr_box,
                         0,
                         base_block,
                         nbr_block);
         }

         Box& nbr_node_box = transformed_nbr_box;
         nbr_node_box.setUpper(
            nbr_node_box.upper() + IntVector::getOne(d_dim));

         Box face_box(base_node_box * nbr_node_box);

         bool is_face = false;
         int face_num = -1;

         if (!face_box.empty()) {
            IntVector box_size(face_box.numberCells());

            int num_width_one = 0;
            int normal_dir = -1;
            for (int d = 0; d < d_dim.getValue(); ++d) {
               if (box_size[d] == 1) {
                  ++num_width_one;
                  normal_dir = d;
               }
            }

            if (num_width_one == 1) {
               is_face = true;
               if (face_box.lower() (normal_dir) ==
                   base_node_box.lower() (normal_dir)) {
                  face_num = 2 * normal_dir;
               } else {
                  face_num = 2 * normal_dir + 1;
               }
            }
         }

         if (is_face) {
            TBOX_ASSERT(face_num >= 0);
            face_neighbors[base_id].insert(std::make_pair(nbr_id, face_num));

            Box transformed_base_box(base_box);
            if (nbr_block != base_block) {
               transformBox(transformed_base_box,
                            0,
                            nbr_block,
                            base_block);
            }

            transformed_base_box.setUpper(
               transformed_base_box.upper() + IntVector::getOne(d_dim));

            nbr_node_box = nbr_box;
            nbr_node_box.setUpper(
               nbr_node_box.upper() + IntVector::getOne(d_dim));

            face_box = transformed_base_box * nbr_node_box;

            if (!face_box.empty()) {
               IntVector box_size(face_box.numberCells());
               face_num = -1;
               int num_width_one = 0;
               int normal_dir = -1;

               for (int d = 0; d < d_dim.getValue(); ++d) {
                  if (box_size[d] == 1) {
                     ++num_width_one;
                     normal_dir = d;
                  }
               }

               if (num_width_one == 1) {
                  is_face = true;
                  if (face_box.lower() (normal_dir) ==
                      nbr_node_box.lower() (normal_dir)) {
                     face_num = 2 * normal_dir;
                  } else {
                     face_num = 2 * normal_dir + 1;
                  }
               } else {
                  TBOX_ERROR(
                     "BaseGridGeometry::findSingularities: Face found on one side of block boundary but not the other.");
               }

               face_neighbors[nbr_id].insert(std::make_pair(base_id, face_num));
            }
         }
      }
   }

   if (!face_neighbors.empty()) {
      if (d_singularity_finder.get() == nullptr) {
         d_singularity_finder.reset(new SingularityFinder(d_dim));
      }
      d_singularity_finder->findSingularities(singularity_blocks,
         chopped_domain,
         *this,
         face_neighbors);
   }
}

/*
 * ************************************************************************
 *
 * Chop domain to eliminate T-junctions.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::chopDomain(
   BoxContainer& chopped_domain)
{
   TBOX_ASSERT(chopped_domain.empty());

   chopped_domain = d_physical_domain;
   chopped_domain.order();

   std::map<BoxId, std::set<BoxId> > dont_chop;
   std::map<BoxId, std::set<BoxId> > neighbors;

   bool breaking_needed = true;
   bool chopped = false;
   while (breaking_needed) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         BlockId base_block(b);

         const std::map<BlockId, Neighbor>& nbrs_of_base = d_block_neighbors[b];

         BoxContainerSingleBlockIterator bi(chopped_domain.begin(base_block));

         chopped = false;
         for ( ; bi != chopped_domain.end(base_block); ++bi) {

            const Box& base_box = *bi;
            const BoxId& base_id = base_box.getBoxId();
            Box base_node_box(base_box);
            base_node_box.setUpper(
               base_node_box.upper() + IntVector::getOne(d_dim));
            IntVector base_node_size(base_node_box.numberCells());

            BoxContainerSingleBlockIterator si(
               chopped_domain.begin(base_block));

            for ( ; si != chopped_domain.end(base_block); ++si) {
               const BoxId& other_id = si->getBoxId();
               if (base_id != other_id) {
                  if (dont_chop[base_id].find(other_id) == dont_chop[base_id].end()) {
                     Box nbr_node_box(*si);
                     nbr_node_box.setUpper(
                        nbr_node_box.upper() + IntVector::getOne(d_dim));

                     Box intersect(base_node_box * nbr_node_box);
                     if (!intersect.empty()) {
                        neighbors[base_id].insert(other_id);
                        neighbors[other_id].insert(base_id);
                        IntVector intersect_size(intersect.numberCells());
                        for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
                           if (intersect_size[d] != 1) {
                              if (intersect_size[d] != base_node_size[d]) {
                                 bool chop_low;
                                 int chop;
                                 if (base_box.lower() (d) != si->lower() (d)) {
                                    chop = std::max<int>(base_box.lower(d),
                                                         si->lower(d));
                                    chop_low = true;
                                 } else {
                                    chop = std::min<int>(base_box.upper(d),
                                                         si->upper(d));
                                    chop_low = false;
                                 }

                                 BoxContainer::iterator box_itr =
                                    chopped_domain.find(base_box);

                                 LocalId local_id(chopped_domain.size());
                                 Box new_box(base_box, local_id, 0);
                                 if (chop_low == true) {
                                    new_box.setLower(d, chop);
                                    box_itr->setUpper(d, chop - 1);
                                 } else {
                                    new_box.setUpper(d, chop);
                                    box_itr->setLower(d, chop + 1);
                                 }
                                 chopped_domain.insert(chopped_domain.end(),
                                    new_box);
                                 chopped = true;

                                 dont_chop[base_id].clear();
                                 for (std::set<BoxId>::iterator nn_itr =
                                         neighbors[base_id].begin();
                                      nn_itr != neighbors[base_id].end();
                                      ++nn_itr) {

                                    dont_chop[*nn_itr].erase(base_id);
                                 }
                                 break;
                              }
                           }
                        }
                     }
                  }
                  if (chopped) {
                     break;
                  } else {
                     dont_chop[base_id].insert(other_id);
                  }
               }
            }
            if (chopped) break;

            for (std::map<BlockId, Neighbor>::const_iterator
                 itr = nbrs_of_base.begin();
                 itr != nbrs_of_base.end(); ++itr) {

               const BlockId& nbr_block = itr->second.getBlockId();

               BoxContainerSingleBlockIterator ni(
                  chopped_domain.begin(nbr_block));

               for ( ; ni != chopped_domain.end(nbr_block); ++ni) {

                  const BoxId& other_id = ni->getBoxId();
                  std::set<BoxId>& nochop = dont_chop[base_id];
                  if (nochop.find(other_id) == nochop.end()) {

                     Box nbr_box(*ni);
                     transformBox(nbr_box,
                                  0,
                                  base_block,
                                  nbr_block);
                     Box nbr_node_box(nbr_box);
                     nbr_node_box.setUpper(
                        nbr_node_box.upper() + IntVector::getOne(d_dim));

                     Box intersect(base_node_box * nbr_node_box);
                     if (!intersect.empty()) {
                        neighbors[base_id].insert(other_id);
                        neighbors[other_id].insert(base_id);
                        IntVector intersect_size(intersect.numberCells());
                        for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
                           if (intersect_size[d] != 1) {
                              if (intersect_size[d] != base_node_size[d]) {
                                 bool chop_low;
                                 int chop;
                                 if (base_box.lower() (d) != nbr_box.lower() (d)) {
                                    chop = std::max<int>(base_box.lower(d),
                                                         nbr_box.lower(d));
                                    chop_low = true;
                                 } else {
                                    chop = std::min<int>(base_box.upper(d),
                                                         nbr_box.upper(d));
                                    chop_low = false;
                                 }

                                 BoxContainer::iterator box_itr =
                                    chopped_domain.find(base_box);

                                 LocalId local_id(chopped_domain.size());
                                 Box new_box(base_box, local_id, 0);
                                 if (chop_low == true) {
                                    new_box.setLower(d, chop);
                                    box_itr->setUpper(d, chop - 1);
                                 } else {
                                    new_box.setUpper(d, chop);
                                    box_itr->setLower(d, chop + 1);
                                 }

                                 chopped_domain.insert(chopped_domain.end(),
                                    new_box);
                                 chopped = true;

                                 nochop.clear();
                                 for (std::set<BoxId>::iterator nn_itr =
                                         neighbors[base_id].begin();
                                      nn_itr != neighbors[base_id].end();
                                      ++nn_itr) {

                                    dont_chop[*nn_itr].erase(base_id);
                                 }
                                 break;
                              }
                           }
                        }
                     }
                  }

                  if (chopped) {
                     break;
                  } else {
                     nochop.insert(other_id);
                  }
               }
               if (chopped) break;
            }
            if (chopped) break;
         }
         if (chopped) break;
      }
      if (!chopped) {
         breaking_needed = false;
      }
   }
}

/*
 *************************************************************************
 * Rotate and shift a box according to the rotation and shift that is
 * used to transform the index space of input_block into the
 * index space of base_block.
 *************************************************************************
 */

bool
BaseGridGeometry::transformBox(
   Box& box,
   int level_number, 
   const BlockId& output_block,
   const BlockId& input_block) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const BlockId::block_t& out_blk = output_block.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[out_blk].find(input_block);
   if (itr != d_block_neighbors[out_blk].end()) {
      const Neighbor& neighbor = itr->second;
      const IntVector& refined_shift = neighbor.getShift(level_number);
      box.rotate(neighbor.getRotationIdentifier());
      box.shift(refined_shift);
      box.setBlockId(output_block);
      return true;
   }
   return false;
}

bool
BaseGridGeometry::transformBox(
   Box& box,
   const IntVector& ratio,
   const BlockId& output_block,
   const BlockId& input_block) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, box, ratio);

   return transformBox(box,
                       getEquivalentLevelNumber(ratio),
                       output_block,
                       input_block);
}

/*
 * ************************************************************************
 *
 * Rotate and shift the boxes in the given array according to the
 * rotation and shift that is used to transformed the index space of
 * input_block into the index space of output_block.
 *
 * ************************************************************************
 */

bool
BaseGridGeometry::transformBoxContainer(
   BoxContainer& boxes,
   const IntVector& ratio,
   const BlockId& output_block,
   const BlockId& input_block) const
{
   const BlockId::block_t& out_blk = output_block.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[out_blk].find(input_block);
   if (itr != d_block_neighbors[out_blk].end()) {
      const Neighbor& neighbor = itr->second;
      const IntVector& refined_shift = neighbor.getShift(getEquivalentLevelNumber(ratio));
      boxes.rotate(neighbor.getRotationIdentifier());
      boxes.shift(refined_shift);
      for (BoxContainer::iterator itr = boxes.begin(); itr != boxes.end();
           ++itr) {
         itr->setBlockId(output_block);
      }
      return true;
   }
   return false;
}

/*
 * ************************************************************************
 *
 * Set block to be the domain of transformed_block in the index space of
 * base_block.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::getTransformedBlock(
   BoxContainer& block,
   const BlockId& base_block,
   const BlockId& transformed_block)
{
   const BlockId::block_t& base_blk = base_block.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[base_blk].find(transformed_block);
   if (itr != d_block_neighbors[base_blk].end()) {
      block = itr->second.getTransformedDomain();
   }
}

/*
 * ************************************************************************
 *
 * Adjust all of the boundary boxes on the level so that they are
 * multiblock-aware.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::adjustMultiblockPatchLevelBoundaries(
   PatchLevel& patch_level)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch_level);
   TBOX_ASSERT(patch_level.getGridGeometry()->getNumberBlocks() == d_number_blocks);

   if (d_number_blocks > 1) {

      t_adjust_multiblock_patch_level_boundaries->start();

      const BoxContainer& d_boxes =
         patch_level.getBoxLevel()->getBoxes();

      IntVector gcw(patch_level.getPatchDescriptor()->getMaxGhostWidth(d_dim));

      for (BlockId::block_t nb = 0; nb < d_number_blocks; ++nb) {

         const BlockId block_id(nb);

         BoxContainer singularity(d_singularity[nb]);
         singularity.refine(patch_level.getRatioToLevelZero());

         BoxContainer pseudo_domain;

         std::map<BlockId, Neighbor>& nbr_map = d_block_neighbors[nb];
         for (std::map<BlockId, Neighbor>::iterator nei = nbr_map.begin();
              nei != nbr_map.end(); ++nei) {
            BoxContainer transformed_domain(nei->second.getTransformedDomain());
            pseudo_domain.spliceFront(transformed_domain);
         }

         pseudo_domain.refine(patch_level.getRatioToLevelZero());

         BoxContainer physical_domain(patch_level.getPhysicalDomain(block_id));
         BoxContainer sing_boxes(singularity);
         pseudo_domain.spliceFront(physical_domain);
         pseudo_domain.spliceFront(sing_boxes);
         pseudo_domain.coalesce();

         BoxContainerSingleBlockIterator mbi(d_boxes.begin(block_id));

         for ( ; mbi != d_boxes.end(block_id); ++mbi) {
            const BoxId& box_id = mbi->getBoxId();
            adjustBoundaryBoxesOnPatch(
               *patch_level.getPatch(box_id),
               pseudo_domain,
               gcw,
               singularity);
         }
      }

      t_adjust_multiblock_patch_level_boundaries->stop();
   }
}

/*
 * ************************************************************************
 *
 * Adjust all of the boundary boxes on the patch so that they are
 * multiblock-aware.
 *
 * ************************************************************************
 */

void
BaseGridGeometry::adjustBoundaryBoxesOnPatch(
   const Patch& patch,
   const BoxContainer& pseudo_domain,
   const IntVector& gcw,
   const BoxContainer& singularity)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, patch, gcw);

   /*
    * Avoid adjusting boundary boxes for the case where we just use
    * a single block, since this is equivalent to not using multiblocks
    * at all.
    */
   if (d_number_blocks > 1) {
      PatchBoundaries boundaries(d_dim);

      getBoundaryBoxes(boundaries,
         patch.getBox(),
         pseudo_domain,
         gcw,
         IntVector::getZero(d_dim));

      std::vector<BoundaryBox> codim_boundaries[SAMRAI::MAX_DIM_VAL];
      std::list<int> boundaries_in_sing[SAMRAI::MAX_DIM_VAL];
      for (int codim = 2; codim <= d_dim.getValue(); ++codim) {

         codim_boundaries[codim - 1] =
            patch.getPatchGeometry()->getCodimensionBoundaries(codim);

         int num_boxes = static_cast<int>(codim_boundaries[codim - 1].size());

         for (int n = 0; n < num_boxes; ++n) {
            Box border_box(codim_boundaries[codim - 1][n].getBox());
            BoxContainer sing_test_list(singularity);
            sing_test_list.intersectBoxes(border_box);
            if (!sing_test_list.empty()) {
               boundaries_in_sing[codim - 1].push_front(n);
            }
         }
      }

      for (int i = 0; i < d_dim.getValue(); ++i) {
         if (!boundaries_in_sing[i].empty()) {
            int old_size = static_cast<int>(boundaries[i].size());
            int nb = 0;
            for (std::list<int>::iterator b = boundaries_in_sing[i].begin();
                 b != boundaries_in_sing[i].end(); ++b) {
               boundaries[i].push_back(codim_boundaries[i][*b]);
               boundaries[i][old_size + nb].setIsMultiblockSingularity(true);
               ++nb;
            }
         }
         patch.getPatchGeometry()->setCodimensionBoundaries(boundaries[i],
            i + 1);
      }

   }

}

/*
 *************************************************************************
 *
 *************************************************************************
 */
Transformation::RotationIdentifier
BaseGridGeometry::getRotationIdentifier(
   const BlockId& dst,
   const BlockId& src) const
{
   TBOX_ASSERT(areNeighbors(dst, src));

   Transformation::RotationIdentifier rotate = Transformation::NO_ROTATE;
   const BlockId::block_t& dst_blk = dst.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[dst_blk].find(src);
   if (itr != d_block_neighbors[dst_blk].end()) {
      rotate = itr->second.getTransformation(0).getRotation();
   }

   return rotate;
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
const IntVector&
BaseGridGeometry::getOffset(
   const BlockId& dst,
   const BlockId& src,
   int level_num) const
{
   TBOX_ASSERT(areNeighbors(dst, src));

   const BlockId::block_t& dst_blk = dst.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[dst_blk].find(src);
   if (itr != d_block_neighbors[dst_blk].end()) {
      return itr->second.getTransformation(level_num).getOffset();
   }

   return IntVector::getOne(d_dim);
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
bool
BaseGridGeometry::areNeighbors(
   const BlockId& block_a,
   const BlockId& block_b) const
{
   bool are_neighbors = false;

   const BlockId::block_t& a_blk = block_a.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[a_blk].find(block_b);
   if (itr != d_block_neighbors[a_blk].end()) {
      are_neighbors = true;
   }

   return are_neighbors;
}

/*
 *************************************************************************
 *
 *************************************************************************
 */
bool
BaseGridGeometry::areSingularityNeighbors(
   const BlockId& block_a,
   const BlockId& block_b) const
{
   bool are_sing_neighbors = false;

   const BlockId::block_t& a_blk = block_a.getBlockValue();
   std::map<BlockId, Neighbor>::const_iterator itr =
      d_block_neighbors[a_blk].find(block_b);
   if (itr != d_block_neighbors[a_blk].end()) {
      if (itr->second.isSingularity()) {
         are_sing_neighbors = true;
      }
   }

   return are_sing_neighbors;
}

void
BaseGridGeometry::setUpRatios(
   const std::vector<IntVector>& ratio_to_coarser)
{
   int max_levels = static_cast<int>(ratio_to_coarser.size());
   TBOX_ASSERT(max_levels > 0);

   d_ratio_to_level_zero.resize(max_levels,
      IntVector(IntVector::getOne(d_dim), d_number_blocks));
   for (int ln = 1; ln < max_levels; ++ln) {
      d_ratio_to_level_zero[ln] = d_ratio_to_level_zero[ln-1] *
                                  ratio_to_coarser[ln];
   }

   if (d_number_blocks > 1) {
      for (int ln = 1; ln < max_levels; ++ln) {
         if (d_ratio_to_level_zero[ln].min() !=
             d_ratio_to_level_zero[ln].max()) {
            d_ratios_are_isotropic = false;
            break;
         }
      }
   }

   if (d_number_blocks > 1 && max_levels > 1) {
      setUpFineLevelTransformations();
   }
}

void
BaseGridGeometry::setUpFineLevelTransformations()
{
   int max_levels = static_cast<int>(d_ratio_to_level_zero.size());

   const IntVector& one_vector = IntVector::getOne(d_dim);

   for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {

      BlockId block_id(b);

      BoxContainer grow_domain;
      computePhysicalDomain(grow_domain, IntVector::getOne(d_dim), block_id);
      TBOX_ASSERT(d_number_blocks == 1 || grow_domain.size() == 1);
      Box domain_box(grow_domain.front());
      grow_domain.unorder();

      grow_domain.grow(one_vector);

      std::map<BlockId,Neighbor>& neighbors(d_block_neighbors[b]);

      for (int ln = 1; ln < max_levels; ++ln) {
         const IntVector& current_ratio = d_ratio_to_level_zero[ln]; 

         for (std::map<BlockId,Neighbor>::iterator ni = neighbors.begin();
              ni != neighbors.end(); ++ni) {

            const BlockId& nbr_block = ni->second.getBlockId();
            const BoxContainer& nbr_domain = ni->second.getTransformedDomain();
            TBOX_ASSERT(nbr_domain.size() == 1);
            const Box& nbr_box = nbr_domain.front();

            Box base_node_box(domain_box);
            base_node_box.setUpper(base_node_box.upper() + IntVector::getOne(d_dim)); 
            Box nbr_node_box(nbr_box);
            nbr_node_box.setUpper(nbr_node_box.upper() + IntVector::getOne(d_dim));

            Box shared_base_nodes(base_node_box*nbr_node_box);
            IntVector shared_size(shared_base_nodes.numberCells());

            IntVector base_bdry_dir(IntVector::getZero(d_dim));
            for (int d = 0; d < d_dim.getValue(); ++d) {
               if (shared_size[d] == 1) {
                  if (shared_base_nodes.upper(static_cast<Box::dir_t>(d)) ==
                      base_node_box.lower(static_cast<Box::dir_t>(d))) {
                     base_bdry_dir[d] = -1;
                  } else if (shared_base_nodes.lower(static_cast<Box::dir_t>(d)) ==
                             base_node_box.upper(static_cast<Box::dir_t>(d))) {
                     base_bdry_dir[d] = 1;
                  }
               }
            }

            Box transformed_base_box(domain_box);
            transformBox(transformed_base_box,
                         0,
                         nbr_block,
                         block_id);

            Box true_nbr_box(nbr_box);
            transformBox(true_nbr_box,
                         0,
                         nbr_block,
                         block_id);

            Box tran_base_node(transformed_base_box);
            tran_base_node.setUpper(tran_base_node.upper() + IntVector::getOne(d_dim));
            Box true_nbr_node(true_nbr_box);
            true_nbr_node.setUpper(true_nbr_node.upper() +  IntVector::getOne(d_dim));
            Box shared_nbr_nodes(tran_base_node * true_nbr_node);

            shared_size = shared_nbr_nodes.numberCells();

            IntVector nbr_bdry_dir(IntVector::getZero(d_dim));
            for (int d = 0; d < d_dim.getValue(); ++d) {
               if (shared_size[d] == 1) {
                  if (shared_nbr_nodes.upper(static_cast<Box::dir_t>(d)) ==
                      true_nbr_node.lower(static_cast<Box::dir_t>(d))) {
                     nbr_bdry_dir[d] = -1;
                  } else if (shared_nbr_nodes.lower(static_cast<Box::dir_t>(d)) ==
                             true_nbr_node.upper(static_cast<Box::dir_t>(d))) {
                     nbr_bdry_dir[d] = 1;
                  }
               }
            }

            Box refined_base(domain_box);
            refined_base.refine(current_ratio);
            Box refined_nbr(true_nbr_box);
            refined_nbr.refine(current_ratio);

            Box refined_base_node(refined_base);
            refined_base_node.setUpper(refined_base_node.upper() + IntVector::getOne(d_dim));
            Box refined_nbr_node(refined_nbr);
            refined_nbr_node.setUpper(refined_nbr_node.upper() + IntVector::getOne(d_dim));

            shared_base_nodes.refine(current_ratio);
            shared_nbr_nodes.refine(current_ratio);
            shared_base_nodes *= refined_base_node;
            shared_nbr_nodes *= refined_nbr_node;
             

            Index base_cell(refined_base.lower());
            Index nbr_cell(refined_nbr.lower());
            for (int d = 0; d < d_dim.getValue(); ++d) {
               if (base_bdry_dir[d] == -1) {
                  base_cell[d] =
                     refined_base.lower(static_cast<Box::dir_t>(d));
               } else if (base_bdry_dir[d] == 1) {
                  base_cell[d] =
                     refined_base.upper(static_cast<Box::dir_t>(d));
               } else {
                  base_cell[d] =
                     shared_base_nodes.lower(static_cast<Box::dir_t>(d)) +
                     ((shared_base_nodes.upper(static_cast<Box::dir_t>(d)) -
                       shared_base_nodes.lower(static_cast<Box::dir_t>(d)))/2);
               }
               if (nbr_bdry_dir[d] == -1) {
                  nbr_cell[d] = refined_nbr.lower(static_cast<Box::dir_t>(d));
               } else if (nbr_bdry_dir[d] == 1) {
                  nbr_cell[d] = refined_nbr.upper(static_cast<Box::dir_t>(d));
               } else {
                  nbr_cell[d] =
                     shared_nbr_nodes.lower(static_cast<Box::dir_t>(d)) +
                     ((shared_nbr_nodes.upper(static_cast<Box::dir_t>(d)) -
                       shared_nbr_nodes.lower(static_cast<Box::dir_t>(d)))/2);
               }
            }

            nbr_cell += nbr_bdry_dir;

            Transformation::RotationIdentifier rotation =
               ni->second.getRotationIdentifier();

            Box box_in_base(base_cell, base_cell, block_id);
            Box box_in_nbr(nbr_cell, nbr_cell, nbr_block);

            box_in_nbr.rotate(rotation);

            IntVector test_shift(box_in_base.lower() - box_in_nbr.lower());

            IntVector final_shift (d_dim);

            Box test_box(refined_nbr);
            test_box.rotate(rotation);
            test_box.shift(test_shift);
            test_box.setBlockId(block_id);
            Box test_box_node(test_box);
            test_box_node.setUpper(test_box_node.upper() + IntVector::getOne(d_dim));
            IntVector sh_num_cells(shared_nbr_nodes.numberCells());
            for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
               if (sh_num_cells[d] ==
                   current_ratio(shared_nbr_nodes.getBlockId().getBlockValue(),d)) {
                  shared_nbr_nodes.setLower(static_cast<Box::dir_t>(d),
                     shared_nbr_nodes.upper(static_cast<Box::dir_t>(d)));
               }
            }
            size_t sh_size = shared_nbr_nodes.size();

            if ((test_box_node * refined_base_node).size() == sh_size) {
               final_shift = test_shift;
            } else {
               IntVector adj_shift(d_dim);
               Box adj(Index(d_dim,-1),Index(d_dim,1),block_id);
               Box::iterator aend(adj.end());
               for (Box::iterator ai(adj.begin()); ai != aend; ++ai) {
                  adj_shift = test_shift + IntVector(*ai);
                  test_box_node = refined_nbr;
                  test_box_node.rotate(rotation);
                  test_box_node.shift(adj_shift);
                  test_box_node.setBlockId(block_id);
                  test_box_node.setUpper(test_box_node.upper() + IntVector::getOne(d_dim));
                   
                  if ((test_box_node * refined_base_node).size() == sh_size) {
                     final_shift = adj_shift;
                     break;
                  }
               }
               if (final_shift != adj_shift) {
                  const IntVector& coarse_shift = ni->second.getTransformation(ln-1).getOffset();;
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     if (coarse_shift[i] == 0) {
                        final_shift[i] = 0;
                     } else if (test_shift[i] % coarse_shift[i] == 0) {
                        final_shift[i] = test_shift[i];
                     } else if ((test_shift[i]-1) % coarse_shift[i] == 0) {
                        final_shift[i] = test_shift[i]-1; 
                     } else if ((test_shift[i]+1) % coarse_shift[i] == 0) {
                        final_shift[i] = test_shift[i]+1;
                     } else {
                        TBOX_ERROR("BaseGridGeometry error...\n"
                           << " Could not compute Level " << ln
                           << " shift between blocks "
                           << block_id.getBlockValue() 
                           << " and " << nbr_block.getBlockValue()
                           << std::endl);
                     }
                  }
               }
            }

            Transformation new_transformation(rotation,
                                              final_shift,
                                              nbr_block,
                                              block_id);

            ni->second.addTransformation(new_transformation, ln);
           
         }
      }
   }

}


/*
 *************************************************************************
 *
 * Print object data to the specified output stream.
 *
 *************************************************************************
 */

void
BaseGridGeometry::printClassData(
   std::ostream& stream) const
{

   stream << "\nBaseGridGeometry::printClassData..." << std::endl;
   stream << "BaseGridGeometry: this = "
          << (BaseGridGeometry *)this << std::endl;
   stream << "d_object_name = " << d_object_name << std::endl;

   const int n = d_physical_domain.size();
   stream << "Number of boxes describing physical domain = " << n << std::endl;
   stream << "Boxes describing physical domain..." << std::endl;
   d_physical_domain.print(stream);

   stream << "\nd_periodic_shift = " << d_periodic_shift << std::endl;

   stream << "d_max_data_ghost_width = " << d_max_data_ghost_width << std::endl;

   stream << "Block neighbor data:\n";

   for (BlockId::block_t bn = 0; bn < d_number_blocks; ++bn) {

      stream << "   Block " << bn << '\n';

      const BlockId block_id(bn);
      const BoxContainer& singularity_boxlist(
         getSingularityBoxContainer(block_id));

      for (ConstNeighborIterator li = begin(block_id);
           li != end(block_id); ++li) {
         const Neighbor& neighbor(*li);
         stream << "      neighbor block " << neighbor.getBlockId() << ':';
         stream << " singularity = " << neighbor.isSingularity() << '\n';
      }

      stream << "      singularity Boxes (" << singularity_boxlist.size() << ")\n";
      for (BoxContainer::const_iterator bi = singularity_boxlist.begin();
           bi != singularity_boxlist.end(); ++bi) {
         stream << "         " << *bi << '\n';
      }

   }

}

/*
 * ************************************************************************
 * ************************************************************************
 */

void
BaseGridGeometry::initializeCallback()
{
   t_find_patches_touching_boundaries = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::findPatchesTouchingBoundaries()");
   t_touching_boundaries_init = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::...TouchingBoundaries()_init");
   t_touching_boundaries_loop = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::...TouchingBoundaries()_loop");
   t_set_geometry_on_patches = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::setGeometryOnPatches()");
   t_set_boundary_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::setBoundaryBoxes()");
   t_set_geometry_data_on_patches = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::set_geometry_data_on_patches");
   t_compute_boundary_boxes_on_level = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::computeBoundaryBoxesOnLevel()");
   t_get_boundary_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::getBoundaryBoxes()");
   t_adjust_multiblock_patch_level_boundaries = tbox::TimerManager::getManager()->
      getTimer("hier::BaseGridGeometry::adjustMultiblockPatchLevelBoundaries()");
}

/*
 *************************************************************************
 *************************************************************************
 */

void
BaseGridGeometry::finalizeCallback()
{
   t_find_patches_touching_boundaries.reset();
   t_touching_boundaries_init.reset();
   t_touching_boundaries_loop.reset();
   t_set_geometry_on_patches.reset();
   t_set_boundary_boxes.reset();
   t_set_geometry_data_on_patches.reset();
   t_compute_boundary_boxes_on_level.reset();
   t_get_boundary_boxes.reset();
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::Neighbor::Neighbor(
   const BlockId& block_id,
   const BoxContainer& domain,
   const std::vector<Transformation>& transformation):
   d_block_id(block_id),
   d_transformed_domain(domain),
   d_transformation(transformation),
   d_is_singularity(false)
{
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::NeighborIterator::NeighborIterator(
   BaseGridGeometry* grid_geometry,
   const BlockId& block_id,
   bool from_start):
   d_grid_geom(grid_geometry),
   d_block_id(block_id)
{
   TBOX_ASSERT(grid_geometry != 0);

   if (from_start) {
      d_nbr_itr =
         grid_geometry->d_block_neighbors[block_id.getBlockValue()].begin();
   } else {
      d_nbr_itr =
         grid_geometry->d_block_neighbors[block_id.getBlockValue()].end();
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::NeighborIterator::NeighborIterator(
   BaseGridGeometry* grid_geometry,
   const BlockId& block_id,
   const BlockId& nbr_block_id):
   d_grid_geom(grid_geometry),
   d_block_id(block_id)
{
   TBOX_ASSERT(grid_geometry != 0);

   d_nbr_itr =
      grid_geometry->d_block_neighbors[block_id.getBlockValue()].find(nbr_block_id);
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::NeighborIterator::NeighborIterator(
   const NeighborIterator& other)
{
   d_grid_geom = other.d_grid_geom;
   d_block_id.setId(other.d_block_id.getBlockValue());
   d_nbr_itr = other.d_nbr_itr;
}

/*
 *************************************************************************
 *************************************************************************
 */

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::ConstNeighborIterator::ConstNeighborIterator(
   const BaseGridGeometry* grid_geometry,
   const BlockId& block_id,
   bool from_start):
   d_grid_geom(grid_geometry),
   d_block_id(block_id)
{
   TBOX_ASSERT(grid_geometry != 0);

   if (from_start) {
      d_nbr_itr =
         grid_geometry->d_block_neighbors[block_id.getBlockValue()].begin();
   } else {
      d_nbr_itr =
         grid_geometry->d_block_neighbors[block_id.getBlockValue()].end();
   }
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::ConstNeighborIterator::ConstNeighborIterator(
   const BaseGridGeometry* grid_geometry,
   const BlockId& block_id,
   const BlockId& nbr_block_id):
   d_grid_geom(grid_geometry),
   d_block_id(block_id)
{
   TBOX_ASSERT(grid_geometry != 0);

   d_nbr_itr =
      grid_geometry->d_block_neighbors[block_id.getBlockValue()].find(nbr_block_id);
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::ConstNeighborIterator::ConstNeighborIterator(
   const ConstNeighborIterator& other)
{
   d_grid_geom = other.d_grid_geom;
   d_block_id.setId(other.d_block_id.getBlockValue());
   d_nbr_itr = other.d_nbr_itr;
}

/*
 *************************************************************************
 *************************************************************************
 */

BaseGridGeometry::ConstNeighborIterator::ConstNeighborIterator(
   const NeighborIterator& other)
{
   d_grid_geom = other.d_grid_geom;
   d_block_id.setId(other.d_block_id.getBlockValue());
   d_nbr_itr = other.d_nbr_itr;
}

/*
 *************************************************************************
 *************************************************************************
 */

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
