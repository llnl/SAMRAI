/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for patch data objects
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchData.h"

namespace SAMRAI {
namespace hier {

const int PatchData::HIER_PATCH_DATA_VERSION = 2;

PatchData::PatchData(
   const Box& domain,
   const IntVector& ghosts):
   d_box(domain),
   d_ghost_box(domain),
   d_ghosts(ghosts),
   d_timestamp(0.0)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(domain, ghosts);

   d_ghost_box.grow(ghosts);
}

void
PatchData::copyFuseable(
   const PatchData& src,
   const BoxOverlap& overlap)
{
   copy(src, overlap);
}

void
PatchData::packStreamFuseable(
   tbox::MessageStream& stream,
   const BoxOverlap& overlap) const
{
   packStream(stream, overlap);
}

void
PatchData::unpackStreamFuseable(
   tbox::MessageStream& stream,
   const BoxOverlap& overlap)
{
   unpackStream(stream, overlap);
}

/*
 *************************************************************************
 *
 * Checks that class and restart file version number are same.  If so,
 * reads in data members common to all patch data from restart database.
 *
 *************************************************************************
 */

void
PatchData::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   int ver = restart_db->getInteger("HIER_PATCH_DATA_VERSION");
   if (ver != HIER_PATCH_DATA_VERSION) {
      TBOX_ERROR("PatchData::getFromRestart() error...\n"
         << "  Restart file version different than class version" << std::endl);
   }

   d_box = restart_db->getDatabaseBox("d_box");
   d_ghost_box = restart_db->getDatabaseBox("d_ghost_box");
   d_timestamp = restart_db->getDouble("d_timestamp");
   restart_db->getIntegerArray("d_ghosts",
      &d_ghosts[0],
      d_ghosts.getDim().getValue());
}

/*
 *************************************************************************
 *
 * Write to restart database data members common to all patch data.
 *
 *************************************************************************
 */

void
PatchData::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("HIER_PATCH_DATA_VERSION", HIER_PATCH_DATA_VERSION);
   restart_db->putDatabaseBox("d_box", d_box);
   restart_db->putDatabaseBox("d_ghost_box", d_ghost_box);
   restart_db->putDouble("d_timestamp", d_timestamp);
   restart_db->putIntegerArray("d_ghosts",
      &d_ghosts[0],
      d_ghosts.getDim().getValue());
}

}
}
