/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Patch container class for patch data objects
 *
 ************************************************************************/
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include <typeinfo>
#include <string>

namespace SAMRAI {
namespace hier {

const int Patch::HIER_PATCH_VERSION = 2;

/*
 *************************************************************************
 *
 * Allocate a patch container but do not instantiate any components.
 *
 *************************************************************************
 */

Patch::Patch(
   const Box& box,
   const std::shared_ptr<PatchDescriptor>& descriptor):
   d_box(box),
   d_descriptor(descriptor),
   d_patch_data(d_descriptor->getMaxNumberRegisteredComponents()),
   d_patch_level_number(-1),
   d_patch_in_hierarchy(false)
{
   TBOX_ASSERT(box.getLocalId() >= 0);
}


/*
 *************************************************************************
 *
 * Calculate the amount of memory space required to allocate the
 * specified component(s).  This information can then be used by a
 * fixed-size memory allocator.
 *
 *************************************************************************
 */

size_t
Patch::getSizeOfPatchData(
   const ComponentSelector& components) const
{
   size_t size = 0;
   const int max_set_component = components.getMaxIndex();

   for (int i = 0; i < max_set_component && components.isSet(i); ++i) {
      size += d_descriptor->getPatchDataFactory(i)->getSizeOfMemory(
            d_box);
   }

   return size;
}

/*
 *************************************************************************
 *
 * Allocate the specified patch data object(s) on the patch.
 *
 *************************************************************************
 */

void
Patch::allocatePatchData(
   const int id,
   const double time)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();

   TBOX_ASSERT((id >= 0) && (id < ncomponents));

   if (ncomponents > static_cast<int>(d_patch_data.size())) {
      d_patch_data.resize(ncomponents);
   }

   if (!checkAllocated(id)) {
      d_patch_data[id] =
         d_descriptor->getPatchDataFactory(id)->allocate(*this);
   }
   d_patch_data[id]->setTime(time);
}

void
Patch::allocatePatchData(
   const ComponentSelector& components,
   const double time)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();
   if (ncomponents > static_cast<int>(d_patch_data.size())) {
      d_patch_data.resize(ncomponents);
   }

   for (int i = 0; i < ncomponents; ++i) {
      if (components.isSet(i)) {
         if (!checkAllocated(i)) {
            d_patch_data[i] =
               d_descriptor->getPatchDataFactory(i)->allocate(*this);
         }
         d_patch_data[i]->setTime(time);
      }
   }
}

/*
 *************************************************************************
 *
 * Deallocate (or set to null) the specified component(s).
 *
 *************************************************************************
 */

void
Patch::deallocatePatchData(
   const ComponentSelector& components)
{
   const int ncomponents = static_cast<int>(d_patch_data.size());
   for (int i = 0; i < ncomponents; ++i) {
      if (components.isSet(i)) {
         d_patch_data[i].reset();
      }
   }
}

/*
 *************************************************************************
 *
 * Set the time stamp for the specified components in the patch.
 *
 *************************************************************************
 */

void
Patch::setTime(
   const double timestamp,
   const ComponentSelector& components)
{
   const int ncomponents = static_cast<int>(d_patch_data.size());
   for (int i = 0; i < ncomponents; ++i) {
      if (components.isSet(i) && d_patch_data[i]) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

void
Patch::setTime(
   const double timestamp)
{
   const int ncomponents = static_cast<int>(d_patch_data.size());
   for (int i = 0; i < ncomponents; ++i) {
      if (d_patch_data[i]) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

/*
 *************************************************************************
 *
 * Checks that class and restart file version numbers are equal.  If so,
 * reads in data from database and have each patch_data item read
 * itself in from the database
 *
 *************************************************************************
 */

void
Patch::getFromRestart(
   const std::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   int ver = restart_db->getInteger("HIER_PATCH_VERSION");
   if (ver != HIER_PATCH_VERSION) {
      TBOX_ERROR("Patch::getFromRestart() error...\n"
         << "   Restart file version different than class version" << std::endl);
   }

   Box box(restart_db->getDatabaseBox("d_box"));
   const LocalId patch_local_id(restart_db->getInteger("d_patch_local_id"));
   int patch_owner = restart_db->getInteger("d_patch_owner");
   int block_id = restart_db->getInteger("d_block_id");
   box.setBlockId(BlockId(block_id));
   d_box.initialize(box,
      patch_local_id,
      patch_owner);

   d_patch_level_number = restart_db->getInteger("d_patch_level_number");
   d_patch_in_hierarchy = restart_db->getBool("d_patch_in_hierarchy");

   d_patch_data.resize(d_descriptor->getMaxNumberRegisteredComponents());

   int namelist_count = restart_db->getInteger("patch_data_namelist_count");
   std::vector<std::string> patch_data_namelist;
   if (namelist_count) {
      patch_data_namelist = restart_db->getStringVector("patch_data_namelist");
   }

   PatchDataRestartManager* pdrm = PatchDataRestartManager::getManager();
   ComponentSelector patch_data_read;

   for (int i = 0; i < static_cast<int>(patch_data_namelist.size()); ++i) {
      std::string& patch_data_name = patch_data_namelist[i];
      int patch_data_index;

      if (!restart_db->isDatabase(patch_data_name)) {
         TBOX_ERROR("Patch::getFromRestart() error...\n"
            << "   patch data" << patch_data_name
            << " not found in restart database" << std::endl);
      }
      std::shared_ptr<tbox::Database> patch_data_database(
         restart_db->getDatabase(patch_data_name));

      patch_data_index = d_descriptor->mapNameToIndex(patch_data_name);

      if ((patch_data_index >= 0) &&
          (pdrm->isPatchDataRegisteredForRestart(patch_data_index))) {
         std::shared_ptr<PatchDataFactory> patch_data_factory(
            d_descriptor->getPatchDataFactory(patch_data_index));
         d_patch_data[patch_data_index] = patch_data_factory->allocate(*this);
         d_patch_data[patch_data_index]->getFromRestart(patch_data_database);
         patch_data_read.setFlag(patch_data_index);
      }
   }

   if (!pdrm->registeredPatchDataMatches(patch_data_read)) {
      TBOX_WARNING("Patch::getFromRestart() warning...\n"
         << "   Some requested patch data components not "
         << "found in restart database" << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Write out the class version number to restart database.  Then,
 * writes out data to restart database and have each patch_data item write
 * itself out to the restart database.  The following data
 * members are written out: d_box, d_patch_number,
 * d_patch_level_number,
 * d_patch_in_hierarchy, d_patch_data[].
 * The database key for all data members is identical to the
 * name of the data member except for the d_patch_data.  These have
 * keys of the form "variable##context" which is the form that they
 * are stored by the patch descriptor.  In addition a list of the
 * patch_data names ("patch_data_namelist") and the number of patch data
 * items saved ("namelist_count") are also written to the database.
 * The PatchDataRestartManager determines which patchdata are written to
 * the database.
 *
 *************************************************************************
 */
void
Patch::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   int i;

   restart_db->putInteger("HIER_PATCH_VERSION", HIER_PATCH_VERSION);
   restart_db->putDatabaseBox("d_box", d_box);
   restart_db->putInteger("d_patch_local_id", d_box.getLocalId().getValue());
   restart_db->putInteger("d_patch_owner", d_box.getOwnerRank());
   restart_db->putInteger("d_block_id",
                          static_cast<int>(d_box.getBlockId().getBlockValue()));
   restart_db->putInteger("d_patch_level_number", d_patch_level_number);
   restart_db->putBool("d_patch_in_hierarchy", d_patch_in_hierarchy);

   int namelist_count = 0;
   const PatchDataRestartManager* pdrm = PatchDataRestartManager::getManager();
   for (i = 0; i < static_cast<int>(d_patch_data.size()); ++i) {
      if (pdrm->isPatchDataRegisteredForRestart(i) && checkAllocated(i)) {
         ++namelist_count;
      }
   }

   std::string patch_data_name;
   std::vector<std::string> patch_data_namelist(namelist_count);
   namelist_count = 0;
   for (i = 0; i < static_cast<int>(d_patch_data.size()); ++i) {
      if (pdrm->isPatchDataRegisteredForRestart(i) && checkAllocated(i)) {
         patch_data_namelist[namelist_count++] =
            patch_data_name = d_descriptor->mapIndexToName(i);
         std::shared_ptr<tbox::Database> patch_data_database(
            restart_db->putDatabase(patch_data_name));
         (d_patch_data[i])->putToRestart(patch_data_database);
      }
   }

   restart_db->putInteger("patch_data_namelist_count", namelist_count);
   if (namelist_count > 0) {
      restart_db->putStringVector("patch_data_namelist", patch_data_namelist);
   }
}

/*
 *************************************************************************
 *
 * Print information about the patch.
 *
 *************************************************************************
 */

int
Patch::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth) const
{
   NULL_USE(depth);

   const tbox::Dimension& dim(d_box.getDim());

   os << border
      << d_box
      << "\tdims: " << d_box.numberCells(0)
   ;
   for (tbox::Dimension::dir_t i = 1; i < dim.getValue(); ++i) {
      os << " X " << d_box.numberCells(i);
   }
   os << "\tsize: " << d_box.size()
      << "\n";
   return 0;
}

std::ostream&
operator << (
   std::ostream& s,
   const Patch& patch)
{
   const int ncomponents = static_cast<int>(patch.d_patch_data.size());
   s << "Patch::box = "
   << patch.d_box << std::endl << std::flush;
   s << "Patch::patch_level_number = " << patch.d_patch_level_number
   << std::endl << std::flush;
   s << "Patch::patch_in_hierarchy = " << patch.d_patch_in_hierarchy
   << std::endl << std::flush;
   s << "Patch::number_components = " << ncomponents
   << std::endl << std::flush;
   for (int i = 0; i < ncomponents; ++i) {
      s << "Component(" << i << ")=";
      if (!patch.d_patch_data[i]) {
         s << "NULL\n";
      } else {
         auto& p = *patch.d_patch_data[i];
         s << typeid(p).name()
         << " [GCW=" << patch.d_patch_data[i]->getGhostCellWidth() << "]\n";
      }
   }
   return s;
}

}
}
