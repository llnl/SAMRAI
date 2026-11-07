/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of edges incident from a box_level of a distributed
 *                box graph.
 *
 ************************************************************************/
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

static const std::string dbgbord;

namespace SAMRAI {
namespace hier {

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector::MappingConnector(
   const tbox::Dimension& dim):Connector(dim)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector::MappingConnector(
   const tbox::Dimension& dim,
   tbox::Database& restart_db):Connector(dim, restart_db)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector::MappingConnector(
   const MappingConnector& other):Connector(other)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector::MappingConnector(
   const BoxLevel& base_box_level,
   const BoxLevel& head_box_level,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state):
   Connector(base_box_level,
             head_box_level,
             base_width,
             parallel_state)
{
}


/*
 ***********************************************************************
 ***********************************************************************
 */
MappingConnector&
MappingConnector::operator = (
   const MappingConnector& rhs)
{
   Connector::operator = (rhs);
   return *this;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector *
MappingConnector::createLocalTranspose() const
{
   const IntVector transpose_gcw = convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         getConnectorWidth());

   MappingConnector* transpose = new MappingConnector(getHead(),
         getBase(),
         transpose_gcw);
   doLocalTransposeWork(transpose);
   return transpose;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnector *
MappingConnector::createTranspose() const
{
   MappingConnector* transpose =
      new MappingConnector(getHead(),
         getBase(),
         convertHeadWidthToBase(getBase().getRefinementRatio(),
            getHead().getRefinementRatio(),
            getConnectorWidth()));

   doTransposeWork(transpose);
   return transpose;
}

/*
 ***********************************************************************
 * Run findMappingErrors and assert that no errors are found.
 ***********************************************************************
 */

void
MappingConnector::assertMappingValidity(
   MappingType map_type) const
{
   size_t nerr = findMappingErrors(map_type);
   if (nerr != 0) {
      tbox::perr << "MappingConnector::assertMappingValidity found\n"
                 << nerr << " errors.\n"
                 << "mapping connector:\n" << format("MAP: ", 2)
                 << "pre-map:\n" << getBase().format("PRE: ", 2)
                 << "post-map:\n" << getHead().format("POST: ", 2)
                 << std::endl;
      TBOX_ERROR("MappingConnector::assertMappingValidity exiting due\n"
         << "to above errors." << std::endl);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t
MappingConnector::findMappingErrors(
   MappingType map_type) const
{
   const tbox::SAMRAI_MPI& mpi(getMPI());
   const int my_rank = mpi.getRank();

   // Need to know whether this is a local map.
   if (map_type == UNKNOWN) {
      if (mpi.getSize() > 1) {
         for (Connector::ConstNeighborhoodIterator ei = begin();
              ei != end(); ++ei) {
            for (Connector::ConstNeighborIterator ni = begin(ei);
                 ni != end(ei); ++ni) {
               if ((*ni).getOwnerRank() != my_rank) {
                  map_type = NOT_LOCAL;
                  break;
               }
            }
            if (map_type == NOT_LOCAL) {
               break;
            }
         }
         if (map_type == UNKNOWN) {
            map_type = LOCAL;
         }
         int tmpi = map_type == LOCAL ? 0 : 1;
         int tmpj; // For some reason, MPI_IN_PLACE is undeclared!
         mpi.Allreduce(&tmpi, &tmpj, 1, MPI_INT, MPI_MAX);
         if (tmpj > 0) {
            map_type = NOT_LOCAL;
         }
      } else {
         map_type = LOCAL;
      }
   }

   /*
    * If not a local map, we need a globalized copy of the head.
    */
   const BoxLevel& new_box_level =
      map_type == LOCAL ? getHead() :
      getHead().getGlobalizedVersion();

   int error_count = 0;

   /*
    * Find old Boxes that changed or disappeared on
    * the new BoxLevel.  There should be a mapping for each
    * Box that changed or disappeared.
    */
   const BoxContainer& old_boxes = getBase().getBoxes();
   for (RealBoxConstIterator ni(old_boxes.realBegin());
        ni != old_boxes.realEnd(); ++ni) {
      const Box& old_box = *ni;
      if (!new_box_level.hasBox(old_box)) {
         // old_box disappeared.  Require a mapping for old_box.
         if (!hasNeighborSet(old_box.getBoxId())) {
            ++error_count;
            tbox::perr << "MappingConnector::findMappingErrors("
                       << error_count
                       << "): old box " << old_box
                       << " disappeared without being mapped." << std::endl;
         }
      } else {
         const Box& new_box = *(new_box_level.getBoxStrict(old_box));
         if (!new_box.isSpatiallyEqual(old_box)) {
            // old_box has changed its box.  A mapping must exist for it.
            if (!hasNeighborSet(old_box.getBoxId())) {
               ++error_count;
               tbox::perr << "MappingConnector::findMappingErrors("
                          << error_count
                          << "): old box " << old_box
                          << " changed to " << new_box
                          << " without being mapped." << std::endl;
            }
         }
      }
   }

   /*
    * All mappings should point from a old Box to a new
    * set of Boxes.
    */
   for (Connector::ConstNeighborhoodIterator ei = begin();
        ei != end(); ++ei) {

      const BoxId& gid = *ei;

      if (!getBase().hasBox(gid)) {
         // Mapping does not go from a old box.
         ++error_count;
         tbox::perr << "MappingConnector::findMappingErrors("
                    << error_count
                    << "): mapping given for nonexistent index " << gid
                    << std::endl;
      } else {
         const Box& old_box = *(getBase().getBoxStrict(gid));

         for (Connector::ConstNeighborIterator ni = begin(ei);
              ni != end(ei); ++ni) {
            const Box& nabr = *ni;

            if (!new_box_level.hasBox(nabr)) {
               ++error_count;
               tbox::perr << "MappingConnector::findMappingErrors("
                          << error_count
                          << "): old box " << old_box
                          << " mapped to nonexistent new box "
                          << nabr << std::endl;
            } else {
               const Box& head_box =
                  *(new_box_level.getBoxStrict(nabr.getBoxId()));
               if (!nabr.isSpatiallyEqual(head_box)) {
                  ++error_count;
                  tbox::perr << "MappingConnector::findMappingErrors("
                             << error_count
                             << "): old box " << old_box
                             << " mapped to neighbor " << nabr
                             << ", which is inconsistent with head box "
                             << head_box << std::endl;
               }
            }

         }
      }
   }

   /*
    * After-boxes should nest in before-boxes grown by the mapping
    * width.
    */

   BoxLevelConnectorUtils blcu;
   std::shared_ptr<BoxLevel> bad_parts;
   std::shared_ptr<MappingConnector> pre_to_bad;
   const Connector* transpose = createTranspose();
   blcu.computeExternalParts(bad_parts,
      pre_to_bad,
      *transpose,
      getConnectorWidth());

   if (pre_to_bad->getLocalNumberOfRelationships() > 0) {
      tbox::perr << "MappingConnector::findMappingErrors() found bad nesting.\n"
                 << "Valid maps' head must nest in base, grown by\n"
                 << "the mapping Connector width.\n"
                 << "mapped boxes and their bad parts:\n"
                 << pre_to_bad->format()
                 << "mapping:\n" << format()
                 << "transpose mapping:\n" << transpose->format()
                 << std::endl;
      error_count += pre_to_bad->getLocalNumberOfRelationships();
   }

   delete transpose;

   return error_count;
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
