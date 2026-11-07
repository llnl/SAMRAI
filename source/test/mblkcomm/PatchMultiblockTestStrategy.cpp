/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for patch data test operations.
 *
 ************************************************************************/

#include "PatchMultiblockTestStrategy.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

using namespace SAMRAI;

// These are used in the cell tagging routine.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/*
 *************************************************************************
 *
 * The constructor and destructor.
 *
 *************************************************************************
 */

PatchMultiblockTestStrategy::PatchMultiblockTestStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   d_variable_src_name.resize(0);
   d_variable_dst_name.resize(0);
   d_variable_depth.resize(0);
   d_variable_src_ghosts.resize(0, hier::IntVector(d_dim));
   d_variable_dst_ghosts.resize(0, hier::IntVector(d_dim));
   d_variable_refine_op.resize(0);
}

PatchMultiblockTestStrategy::~PatchMultiblockTestStrategy()
{
}

/*
 *************************************************************************
 *
 * Routines for reading variable and refinement data from input.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::readVariableInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   std::vector<std::string> var_keys = db->getAllKeys();
   int nkeys = static_cast<int>(var_keys.size());

   d_variable_src_name.resize(nkeys);
   d_variable_dst_name.resize(nkeys);
   d_variable_depth.resize(nkeys);
   d_variable_src_ghosts.resize(nkeys, hier::IntVector(d_dim, 0));
   d_variable_dst_ghosts.resize(nkeys, hier::IntVector(d_dim, 0));
   d_variable_refine_op.resize(nkeys);

   for (int i = 0; i < nkeys; ++i) {

      std::shared_ptr<tbox::Database> var_db(db->getDatabase(var_keys[i]));

      if (var_db->keyExists("src_name")) {
         d_variable_src_name[i] = var_db->getString("src_name");
      } else {
         TBOX_ERROR("Variable input error: No `src_name' string found for "
            << "key = " << var_keys[i] << std::endl);
      }

      if (var_db->keyExists("dst_name")) {
         d_variable_dst_name[i] = var_db->getString("dst_name");
      } else {
         TBOX_ERROR("Variable input error: No `dst_name' string found for "
            << "key = " << var_keys[i] << std::endl);
      }

      if (var_db->keyExists("depth")) {
         d_variable_depth[i] = var_db->getInteger("depth");
      } else {
         d_variable_depth[i] = 1;
      }

      if (var_db->keyExists("src_ghosts")) {
         int* tmp_ghosts = &d_variable_src_ghosts[i][0];
         var_db->getIntegerArray("src_ghosts", tmp_ghosts, d_dim.getValue());
      }

      if (var_db->keyExists("dst_ghosts")) {
         int* tmp_ghosts = &d_variable_dst_ghosts[i][0];
         var_db->getIntegerArray("dst_ghosts", tmp_ghosts, d_dim.getValue());
      }

      if (var_db->keyExists("refine_operator")) {
         d_variable_refine_op[i] = var_db->getString("refine_operator");
      } else {
         d_variable_refine_op[i] = "NO_REFINE";
      }

   }

}

void PatchMultiblockTestStrategy::readRefinementInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   std::vector<std::string> box_keys = db->getAllKeys();
   int nkeys = static_cast<int>(box_keys.size());

   d_refine_level_boxes.resize(nkeys);
   for (int i = 0; i < nkeys; ++i) {
      std::vector<tbox::DatabaseBox> db_box_vector =
         db->getDatabaseBoxVector(box_keys[i]);
      d_refine_level_boxes[i] = db_box_vector;
   }

}

/*
 *************************************************************************
 *
 * Tag cells on level specified in input box array for refinement.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::tagCellsInInputBoxes(
   hier::Patch& patch,
   int level_number,
   int tag_index)
{

   if (level_number < static_cast<int>(d_refine_level_boxes.size())) {

      std::shared_ptr<pdat::CellData<int> > tags(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
            patch.getPatchData(tag_index)));
      TBOX_ASSERT(tags);
      tags->fillAll(0);

      const hier::Box pbox = patch.getBox();

      for (hier::BoxContainer::iterator k =
              d_refine_level_boxes[level_number].begin();
           k != d_refine_level_boxes[level_number].end(); ++k) {
         tags->fill(1, *k * pbox, 0);
      }

   }

}

/*
 *************************************************************************
 *
 * Blank physical boundary and pre/postprocess
 * so tester isn't required to implement them when not needed.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(gcw);
}

void PatchMultiblockTestStrategy::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const std::shared_ptr<hier::VariableContext>& context,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(context);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PatchMultiblockTestStrategy::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const std::shared_ptr<hier::VariableContext>& context,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(context);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}
