/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Lookup table to aid in BoundaryBox construction
 *
 ************************************************************************/
#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

BoundaryLookupTable *
BoundaryLookupTable::s_lookup_table_instance[SAMRAI::MAX_DIM_VAL] = { 0 };

tbox::StartupShutdownManager::Handler
BoundaryLookupTable::s_finalize_handler(
   0,
   0,
   0,
   BoundaryLookupTable::finalizeCallback,
   tbox::StartupShutdownManager::priorityBoundaryLookupTable);

/*
 *************************************************************************
 *
 * Lookup table constructor and destructor.
 *
 *************************************************************************
 */

BoundaryLookupTable::BoundaryLookupTable(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   if (d_table[0].empty()) {
      const tbox::Dimension::dir_t dim_val = d_dim.getValue();
      int factrl[SAMRAI::MAX_DIM_VAL + 1];
      factrl[0] = 1;
      for (int i = 1; i <= dim_val; ++i) {
         factrl[i] = i * factrl[i - 1];
      }

      d_ncomb.resize(dim_val);
      d_max_li.resize(dim_val);
      for (tbox::Dimension::dir_t codim = 1; codim <= dim_val; ++codim) {
         tbox::Dimension::dir_t cdm1 = static_cast<tbox::Dimension::dir_t>(codim - 1);
         d_ncomb[cdm1] = factrl[dim_val]
            / (factrl[codim] * factrl[dim_val - codim]);

         std::vector<int> work(codim * d_ncomb[cdm1]);

         int recursive_work[SAMRAI::MAX_DIM_VAL];
         int recursive_work_lvl = 0;
         int* recursive_work_ptr;
         buildTable(&work[0], recursive_work,
            recursive_work_lvl, recursive_work_ptr, codim, 1);

         d_table[cdm1].resize(d_ncomb[cdm1]);
         for (tbox::Dimension::dir_t j = 0; j < d_ncomb[cdm1]; ++j) {
            d_table[cdm1][j].resize(codim);
            for (tbox::Dimension::dir_t k = 0; k < codim; ++k) {
               d_table[cdm1][j][k] = static_cast<tbox::Dimension::dir_t>(work[j * codim + k] - 1);
            }
         }

         d_max_li[cdm1] = d_ncomb[cdm1] * (1 << codim);
      }
   }

   buildBoundaryDirectionVectors();
}

/*
 *************************************************************************
 *
 * Recursive function that computes the combinations in the lookup
 * table.
 *
 *************************************************************************
 */

void
BoundaryLookupTable::buildTable(
   int* table,
   int(&work)[SAMRAI::MAX_DIM_VAL],
   int& rec_level,
   int *& ptr,
   const int codim,
   const int ibeg)
{
   ++rec_level;
   if (rec_level == 1) {
      ptr = table;
   }
   int iend = d_dim.getValue() - codim + rec_level;
   for (int i = ibeg; i <= iend; ++i) {
      work[rec_level - 1] = i;
      if (rec_level != codim) {
         buildTable(ptr, work, rec_level, ptr, codim, i + 1);
      } else {
         for (int j = 0; j < codim; ++j) {
            *(ptr + j) = work[j];
         }
         ptr += codim;
      }
   }
   --rec_level;
}

/*
 *************************************************************************
 *
 * Build table of IntVectors indication locations of boundaries relative
 * to a patch.
 *
 *************************************************************************
 */

void
BoundaryLookupTable::buildBoundaryDirectionVectors()
{

   d_bdry_dirs.resize(d_dim.getValue());

   for (tbox::Dimension::dir_t i = 0; i < d_dim.getValue(); ++i) {
      d_bdry_dirs[i].resize(d_max_li[i], IntVector::getZero(d_dim));
      tbox::Dimension::dir_t codim = static_cast<tbox::Dimension::dir_t>(i + 1);

      for (int loc = 0; loc < d_max_li[i]; ++loc) {
         const std::vector<tbox::Dimension::dir_t>& dirs = getDirections(loc, codim);

         for (int d = 0; d < static_cast<int>(dirs.size()); ++d) {

            if (isUpper(loc, codim, d)) {

               d_bdry_dirs[i][loc](dirs[d]) = 1;

            } else {

               d_bdry_dirs[i][loc](dirs[d]) = -1;

            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Free statics.
 *
 *************************************************************************
 */
void
BoundaryLookupTable::finalizeCallback()
{
   for (int i = 0; i < SAMRAI::MAX_DIM_VAL; ++i) {
      if (s_lookup_table_instance[i]) {
         delete s_lookup_table_instance[i];
      }
      s_lookup_table_instance[i] = 0;
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
