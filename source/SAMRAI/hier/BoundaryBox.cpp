/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   BoundaryBox representing a portion of the physical boundary
 *
 ************************************************************************/
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/BoundaryLookupTable.h"

namespace SAMRAI {
namespace hier {

BoundaryBox::BoundaryBox(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_box(dim),
   d_bdry_type(-1),
   d_location_index(-1)
{
}

BoundaryBox::BoundaryBox(
   const Box& box,
   const int bdry_type,
   const int location_index):
   d_dim(box.getDim()),
   d_box(box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   BoundaryLookupTable* blut = BoundaryLookupTable::getLookupTable(d_dim);
   const std::vector<int>& location_index_max = blut->getMaxLocationIndices();

   TBOX_ASSERT((bdry_type >= 1) && (bdry_type <= d_dim.getValue()));
   TBOX_ASSERT(location_index >= 0);
   TBOX_ASSERT(location_index < location_index_max[bdry_type - 1]);
#endif

   d_bdry_type = bdry_type;

   d_location_index = location_index;

   d_is_mblk_singularity = false;
}

BoundaryBox::BoundaryOrientation
BoundaryBox::getBoundaryOrientation(
   const int dir) const
{
   TBOX_ASSERT(dir < d_dim.getValue());

   BoundaryLookupTable* blut =
      BoundaryLookupTable::getLookupTable(d_dim);

   int bdry_dir =
      blut->getBoundaryDirections(d_bdry_type)[d_location_index](dir);

   TBOX_ASSERT(bdry_dir == -1 || bdry_dir == 0 || bdry_dir == 1);

   BoundaryOrientation retval;

   if (bdry_dir == -1) {
      retval = LOWER;
   } else if (bdry_dir == 0) {
      retval = MIDDLE;
   } else {
      retval = UPPER;
   }

   return retval;
}

}
}
