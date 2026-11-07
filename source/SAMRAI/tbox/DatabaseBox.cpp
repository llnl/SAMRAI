/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A box structure representing a portion of the AMR index space
 *
 ************************************************************************/

#include "SAMRAI/tbox/DatabaseBox.h"

namespace SAMRAI {
namespace tbox {

DatabaseBox::DatabaseBox(
   const Dimension& dim,
   const int* lower,
   const int* upper)
{
   const int dim_val = dim.getValue();
   TBOX_ASSERT(dim_val >= 0);
   TBOX_ASSERT(dim_val <= SAMRAI::MAX_DIM_VAL);
   TBOX_ASSERT(dim_val == 0 || lower != nullptr);
   TBOX_ASSERT(dim_val == 0 || upper != nullptr);

#ifdef DEBUG_CHECK_ASSERTIONS
   if (dim_val > 0) {
      TBOX_ASSERT(lower != nullptr);
      TBOX_ASSERT(upper != nullptr);
   }
#endif

   d_data.d_dimension = dim_val;

   for (int i = 0; i < d_data.d_dimension; i++) {
      d_data.d_lo[i] = lower[i];
      d_data.d_hi[i] = upper[i];
   }
   for (int j = d_data.d_dimension; j < SAMRAI::MAX_DIM_VAL; j++) {
      d_data.d_lo[j] = 0;
      d_data.d_hi[j] = 0;
   }
}

}
}
