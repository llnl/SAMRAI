/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A flattened representation of a hierarchy
 *
 ************************************************************************/
#include "SAMRAI/hier/FlattenedHierarchy.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"


namespace SAMRAI {
namespace hier {

/*
 ***************************************************************************
 * Constructor evaluates the hierarchy and fills the visible boxes container
 ***************************************************************************
 */

FlattenedHierarchy::FlattenedHierarchy(
   const PatchHierarchy& hierarchy,
   int coarsest_level,
   int finest_level)
: d_coarsest_level(coarsest_level),
  d_finest_level(finest_level),
  d_patch_hierarchy(&hierarchy)
{
   int num_levels = hierarchy.getNumberOfLevels();
   TBOX_ASSERT(coarsest_level >= 0);
   TBOX_ASSERT(coarsest_level <= finest_level);
   TBOX_ASSERT(finest_level < num_levels);

   LocalId local_id(0);

   d_visible_boxes.resize(num_levels);
   std::vector<int> local_num_boxes(num_levels, 0);

   for (int ln = coarsest_level; ln <= finest_level; ++ln) {
      const std::shared_ptr<PatchLevel>& current_level =
         hierarchy.getPatchLevel(ln);

      int& num_boxes = local_num_boxes[ln];

      if (ln != finest_level) {

         const Connector& coarse_to_fine =
            current_level->findConnector(
               *(hierarchy.getPatchLevel(ln+1)),
               IntVector::getOne(hierarchy.getDim()),
               CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         const IntVector& connector_ratio = coarse_to_fine.getRatio();

         for (PatchLevel::iterator ip(current_level->begin());
              ip != current_level->end(); ++ip) {

            const std::shared_ptr<Patch>& patch = *ip;
            const Box& box = patch->getBox();
            const BlockId& block_id = box.getBlockId();
            const BoxId& box_id = box.getBoxId();
            BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];

            BoxContainer coarse_boxes(box);

            BoxContainer fine_nbr_boxes;
            if (coarse_to_fine.hasNeighborSet(box_id)) {
               coarse_to_fine.getNeighborBoxes(box_id, fine_nbr_boxes);
            }
            if (!fine_nbr_boxes.empty()) {
               BoxContainer fine_boxes;
               for (SAMRAI::hier::RealBoxConstIterator
                    nbr_itr = fine_nbr_boxes.realBegin();
                    nbr_itr != fine_nbr_boxes.realEnd(); ++nbr_itr) {
                  if (nbr_itr->getBlockId() == block_id) {
                     fine_boxes.pushBack(*nbr_itr);
                  }
               }

               fine_boxes.coarsen(connector_ratio);

               coarse_boxes.removeIntersections(fine_boxes);
               coarse_boxes.coalesce();
            }

            for (BoxContainer::iterator itr =
                 coarse_boxes.begin(); itr != coarse_boxes.end(); ++itr) {

               Box new_box(*itr, local_id, box_id.getOwnerRank());
               ++local_id;
               visible_boxes.insert(visible_boxes.end(), new_box);
               ++num_boxes;
            }
         }
      } else {
         for (PatchLevel::iterator ip(current_level->begin());
              ip != current_level->end(); ++ip) {
            const std::shared_ptr<Patch>& patch = *ip;
            const Box& box = patch->getBox();
            const BoxId& box_id = box.getBoxId();
            BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];

            Box new_box(box, local_id, box.getOwnerRank());
            ++local_id;
            visible_boxes.insert(visible_boxes.end(), new_box);
            ++num_boxes;
         }
      }
   }

   const tbox::SAMRAI_MPI& hier_mpi = hierarchy.getMPI();

   /*
    * Give each visible box a globally unique LocalId.
    */
   if (hier_mpi.getSize() > 1) {
      std::vector<int> global_num_boxes = local_num_boxes;
      hier_mpi.AllReduce(
         &global_num_boxes[0],
         num_levels,
         MPI_SUM);

      std::vector<int> scanned_num_boxes(num_levels);
      hier_mpi.Scan(
         &local_num_boxes[0],
         &scanned_num_boxes[0],
         num_levels, MPI_INT, MPI_SUM);


      for (int ln = coarsest_level; ln <= finest_level; ++ln) {
         auto& level_box_map = d_visible_boxes[ln];

         LocalId local_id(0);
         for (int i = 0; i < ln; ++i) {
            local_id += global_num_boxes[i];
         }
         local_id += (scanned_num_boxes[ln] - local_num_boxes[ln]);

         for (auto& boxes : level_box_map) {
            for (auto& box : boxes.second) {
               box.setLocalId(local_id);
               ++local_id;
            }
         }
      }
   }
}

/*
 **************************************************************************
 * Destructor
 **************************************************************************
 */

}
}
