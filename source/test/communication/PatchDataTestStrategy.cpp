/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for patch data test operations.
 *
 ************************************************************************/

#include "PatchDataTestStrategy.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

namespace SAMRAI {


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

PatchDataTestStrategy::PatchDataTestStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   d_variable_src_name.resize(0);
   d_variable_dst_name.resize(0);
   d_variable_depth.resize(0);
   d_variable_src_ghosts.resize(0, hier::IntVector(d_dim));
   d_variable_dst_ghosts.resize(0, hier::IntVector(d_dim));
   d_variable_coarsen_op.resize(0);
   d_variable_refine_op.resize(0);
}

PatchDataTestStrategy::~PatchDataTestStrategy()
{
}

/*
 *************************************************************************
 *
 * Routines for reading variable and refinement data from input.
 *
 *************************************************************************
 */

void PatchDataTestStrategy::readVariableInput(
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
   d_variable_coarsen_op.resize(nkeys);
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

      if (var_db->keyExists("coarsen_operator")) {
         d_variable_coarsen_op[i] = var_db->getString("coarsen_operator");
      } else {
         d_variable_coarsen_op[i] = "NO_COARSEN";
      }

      if (var_db->keyExists("refine_operator")) {
         d_variable_refine_op[i] = var_db->getString("refine_operator");
      } else {
         d_variable_refine_op[i] = "NO_REFINE";
      }

   }

}

/*
 *************************************************************************
 *
 * Blank physical boundary and pre/postprocess coarsen/refine operations
 * so tester isn't required to implement them when not needed.
 *
 *************************************************************************
 */

void PatchDataTestStrategy::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(gcw);
}

void PatchDataTestStrategy::preprocessRefine(
   const hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PatchDataTestStrategy::postprocessRefine(
   const hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PatchDataTestStrategy::preprocessCoarsen(
   const hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void PatchDataTestStrategy::postprocessCoarsen(
   const hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

}
