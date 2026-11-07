/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Statistic.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int Statistic::ARRAY_INCREMENT = 100;
const int Statistic::TBOX_STATISTIC_VERSION = 1;
double Statistic::s_empty_seq_tag_entry = -99999999.;

/*
 *************************************************************************
 *
 * Statistic constructor and destructor.
 *
 *************************************************************************
 */

Statistic::Statistic(
   const std::string& name,
   const std::string& stat_type,
   int instance_id)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(!stat_type.empty());
   TBOX_ASSERT(instance_id > -1);

   if (!(stat_type == "PROC_STAT" || stat_type == "PATCH_STAT")) {
      TBOX_ERROR(
         "Invalid stat_type passed to Statistic constructor"
         << "   object name = " << name
         << ".\n   Valid stat types are `PROC_STAT' and `PATCH_STAT'"
         << std::endl);
   }

   d_object_name = name;
   d_instance_id = instance_id;

   if (stat_type == "PROC_STAT") {
      d_stat_type = PROC_STAT;
   } else if (stat_type == "PATCH_STAT") {
      d_stat_type = PATCH_STAT;
   }

   d_seq_counter = 0;
   d_total_patch_entries = 0;
   d_proc_stat_array_size = 0;
   d_patch_stat_array_size = 0;

}

Statistic::~Statistic()
{
   reset();
}

/*
 *************************************************************************
 *
 * Utility functions to record statistic record data.
 *
 *************************************************************************
 */

void
Statistic::recordProcStat(
   double value,
   int seq_num)
{
   if (getType() != "PROC_STAT") {
      TBOX_ERROR("Statistic::recordProcStat error...\n"
         << "    Statistic type is `PATCH_STAT'" << std::endl);
   }

   /*
    * Resize array of processor stats, if necessary.
    */
   checkArraySizes(seq_num);

   /*
    * The statistic maintains its own sequence counter but the user may
    * override the data maintained in the statistic by supplying a seq_num.
    * The seq_num argument for this method is -1 if the user does not
    * supply an alternative value.  The value is recorded using the following
    * logic:
    * 1) If seq_num < 0, increment counter and record a new value at seq_cnt.
    * 2) If 0 < seq_num < seq_cnt, overwrite the previous record at
    *    seq_num with the new value.
    * 3) If seq_cnt < seq_num, set seq_cnt = seq_num and add new record
    *    at seq_cnt.
    */
   if (seq_num < 0) {
      d_proc_array[d_seq_counter].value = value;
   } else if (seq_num < d_seq_counter) {
      d_proc_array[seq_num].value = value;
   } else {
      d_seq_counter = seq_num;
      d_proc_array[d_seq_counter].value = value;
   }
   ++d_seq_counter;
}

void
Statistic::recordPatchStat(
   int patch_num,
   double value,
   int seq_num)
{
   if (getType() != "PATCH_STAT") {
      TBOX_ERROR("Statistic::recordPatchStat error...\n"
         << "    Statistic type is `PROC_STAT'" << std::endl);
   }

   /*
    * Resize array of processor stats, if necessary.
    */
   checkArraySizes(seq_num);

   /*
    * The patch statistic differs from the processor statistic in that
    * each entry of the array of seq numbers contains a LIST of
    * PatchStatRecord entries, one for each patch on the processor.
    * The recording logic is thus slightly different than the Processor stat:
    *
    *   If seq_num < seq_counter {
    *     - Check the entries in the list of records at array index seq_num.
    *       Iterate through the list and see if the patch_id of any record
    *       matches the specified patch_num.
    *
    *       If patch_num entry exists {
    *          - overwrite existing entry
    *       } else {
    *          - create a new entry and append to end of list
    *
    *   }
    *
    *   If seq_num >= seq_counter
    *      - create new entry and append to end of list at the seq_num
    *        array index.
    *      - set seq_counter = seq_num
    *   }
    */
   if (seq_num < d_seq_counter) {
      std::list<Statistic::PatchStatRecord>& records =
         d_patch_array[seq_num].patch_records;
      bool found_patch_id = false;
      std::list<Statistic::PatchStatRecord>::iterator ir = records.begin();
      for ( ; ir != records.end(); ++ir) {
         if (ir->patch_id == patch_num) {
            ir->value = value;
            found_patch_id = true;
         }
      }
      if (!found_patch_id) {
         PatchStatRecord patchitem_record;
         patchitem_record.value = value;
         patchitem_record.patch_id = patch_num;
         d_patch_array[seq_num].patch_records.push_back(patchitem_record);
         ++d_total_patch_entries;
      }

   }

   if (seq_num >= d_seq_counter) {
      PatchStatRecord patchitem_record;
      patchitem_record.value = value;
      patchitem_record.patch_id = patch_num;
      d_patch_array[seq_num].patch_records.push_back(patchitem_record);
      ++d_total_patch_entries;
      d_seq_counter = seq_num + 1;
   }

}

/*
 *************************************************************************
 *
 * Utility function for communicating statistic data in parallel.
 *
 * Stream data size includes 4 ints (instance id, proc rank,
 *                                   stat type, seq length).
 *
 * Additionally, data stream size includes space needed for statistic
 * data values:
 *
 *    o for processor stat, this is 1 double (value) for each seq entry.
 *
 *    o for patch stat, this is 1 int (#patches) for each sequence
 *      entry + 1 int (patch_id) + 1 double (value) for each patch
 *      entry.
 *
 *************************************************************************
 */

size_t
Statistic::getDataStreamSize()
{
   size_t byte_size = MessageStream::getSizeof<int>(4);
   if (d_stat_type == PROC_STAT) {
      byte_size += MessageStream::getSizeof<double>(d_seq_counter);
   } else { // d_stat_type == PATCH_STAT
      byte_size += MessageStream::getSizeof<int>(d_seq_counter);
      byte_size += MessageStream::getSizeof<int>(d_total_patch_entries);
      byte_size += MessageStream::getSizeof<double>(d_total_patch_entries);
   }
   return byte_size;
}

void
Statistic::packStream(
   MessageStream& stream)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0) {
      TBOX_ERROR("Statistic::packStream error...\n"
         << "    Processor zero should not pack stat data" << std::endl);
   }

   int num_int = 4;
   int num_double = 0;
   if (d_stat_type == PROC_STAT) {
      num_double = d_seq_counter;
   }
   if (d_stat_type == PATCH_STAT) {
      num_int += d_seq_counter + d_total_patch_entries;
      num_double = d_total_patch_entries;
   }
   std::vector<int> idata(num_int);
   std::vector<double> ddata(num_double);

   idata[0] = mpi.getRank();
   idata[1] = d_instance_id;
   idata[2] = d_stat_type;
   idata[3] = d_seq_counter;

   int is = 0;
   if (d_stat_type == PROC_STAT) {

      for (is = 0; is < d_seq_counter; ++is) {
         ddata[is] = d_proc_array[is].value;
      }

   } else {  // d_stat_type == PATCH_STAT

      int mark = 4 + d_seq_counter;
      int isr = 0;

      for (is = 0; is < d_seq_counter; ++is) {
         std::list<Statistic::PatchStatRecord>& lrec =
            d_patch_array[is].patch_records;
         idata[4 + is] = static_cast<int>(lrec.size());

         std::list<Statistic::PatchStatRecord>::iterator ilr = lrec.begin();
         for ( ; ilr != lrec.end(); ++ilr) {
            idata[mark + isr] = ilr->patch_id;
            ddata[isr] = ilr->value;
            ++isr;
         }
      }
   }

   if (num_int > 0) {
      stream.pack(&idata[0], num_int);
   }
   if (num_double > 0) {
      stream.pack(&ddata[0], num_double);
   }

}

void
Statistic::unpackStream(
   MessageStream& stream)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() != 0) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Only processor zero should unpack stat data" << std::endl);
   }

   int src_rank, stat_id, stat_type, seq_len;

   stream >> src_rank;
   stream >> stat_id;
   stream >> stat_type;
   stream >> seq_len;

   if (src_rank == 0) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "     Processor zero should not send stat data" << std::endl);
   }
   if (stat_id != d_instance_id) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Incompatible statistic number ids" << std::endl);
   }
   if (stat_type != d_stat_type) {
      TBOX_ERROR("Statistic::unpackStream error...\n"
         << "    Incompatible statistic types" << std::endl);
   }

   int is;
   if (d_stat_type == PROC_STAT) {

      std::vector<double> ddata(seq_len);

      if (seq_len > 0) {
         stream.unpack(&ddata[0], seq_len);
         for (is = 0; is < seq_len; ++is) {
            recordProcStat(ddata[is], is);
         }
      }

   } else { // d_stat_type == PATCH_STAT

      if (seq_len > 0) {
         std::vector<int> inum_patches_data(seq_len);
         stream.unpack(&inum_patches_data[0], seq_len);

         int total_seq_items = 0;
         for (is = 0; is < seq_len; ++is) {
            total_seq_items += inum_patches_data[is];
         }

         std::vector<int> ipatch_num_data(total_seq_items);
         std::vector<double> ddata(total_seq_items);

         stream.unpack(&ipatch_num_data[0], total_seq_items);
         stream.unpack(&ddata[0], total_seq_items);

         int isr = 0;
         for (is = 0; is < seq_len; ++is) {
            for (int ipsr = 0; ipsr < inum_patches_data[is]; ++ipsr) {
               recordPatchStat(ipatch_num_data[isr], ddata[isr], is);
               ++isr;
            }
         }
      }
   }

}

void
Statistic::printClassData(
   std::ostream& stream,
   int precision) const
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   TBOX_ASSERT(precision > 0);

   stream.precision(precision);

   stream << "Local Data for " << getType() << " : " << getName() << std::endl;
   stream << "   Processor id = " << mpi.getRank() << std::endl;
   stream << "   Instance id = " << getInstanceId() << std::endl;
   stream << "   Sequence length = " << getStatSequenceLength() << std::endl;

   int is = 0;
   if (d_stat_type == PROC_STAT) {
      for (is = 0; is < d_seq_counter; ++is) {
         stream << "     sequence[" << is
                << "] : value = " << d_proc_array[is].value << std::endl;
      }
   } else {
      for (is = 0; is < d_seq_counter; ++is) {
         stream << "     sequence[" << is
                << "]" << std::endl;

         const std::list<Statistic::PatchStatRecord>& psrl =
            d_patch_array[is].patch_records;
         std::list<Statistic::PatchStatRecord>::const_iterator ilr =
            psrl.begin();
         for ( ; ilr != psrl.end(); ++ilr) {
            stream << "        patch # = " << ilr->patch_id
                   << " : value = " << ilr->value << std::endl;
         }
      }
   }

}

void
Statistic::checkArraySizes(
   int seq_num)
{
   /*
    * Verify that seq_num is less than array size.  If so, drop through.
    * If not, resize and initialize the array.
    */
   int high_mark = MathUtilities<int>::Max(seq_num, d_seq_counter);

   if (d_stat_type == PROC_STAT) {

      if (high_mark >= d_proc_stat_array_size) {
         int old_array_size = d_proc_stat_array_size;
         d_proc_stat_array_size += ARRAY_INCREMENT;
         d_proc_array.resize(d_proc_stat_array_size);
         for (int i = old_array_size; i < d_proc_stat_array_size; ++i) {
            d_proc_array[i].value = s_empty_seq_tag_entry;
         }

      }

   } else if (d_stat_type == PATCH_STAT) {

      if (high_mark >= d_patch_stat_array_size) {
         d_patch_stat_array_size += ARRAY_INCREMENT;
         d_patch_array.resize(d_patch_stat_array_size);
      }

   }

}

void
Statistic::putToRestart(
   const std::shared_ptr<Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("TBOX_STATISTIC_VERSION", TBOX_STATISTIC_VERSION);

   restart_db->putString("object_name", d_object_name);
   restart_db->putInteger("stat_type", d_stat_type);
   restart_db->putInteger("instance_id", d_instance_id);
   restart_db->putInteger("seq_counter", d_seq_counter);
   restart_db->putInteger("total_patch_entries", d_total_patch_entries);
   restart_db->putInteger("proc_stat_array_size", d_proc_stat_array_size);
   restart_db->putInteger("patch_stat_array_size", d_patch_stat_array_size);

   int i;

   if (d_stat_type == PROC_STAT) {
      std::vector<double> ddata(d_seq_counter);
      for (i = 0; i < d_seq_counter; ++i) {
         ddata[i] = d_proc_array[i].value;
      }

      if (d_seq_counter > 0) {
         restart_db->putDoubleVector("ddata", ddata);
      }

   }

   if (d_stat_type == PATCH_STAT) {
      std::vector<int> idata(d_seq_counter + d_total_patch_entries);
      std::vector<double> ddata(d_total_patch_entries);

      int il = 0;
      int mark = d_seq_counter;

      for (i = 0; i < d_seq_counter; ++i) {
         const std::list<Statistic::PatchStatRecord>& records =
            d_patch_array[i].patch_records;
         idata[i] = static_cast<int>(records.size());  // # patches at seq num
         std::list<Statistic::PatchStatRecord>::const_iterator ir =
            records.begin();
         for ( ; ir != records.end(); ++ir) {
            idata[mark + il] = ir->patch_id;
            ddata[il] = ir->value;
            ++il;
         }
      }

      if (d_seq_counter > 0) {
         restart_db->putIntegerVector("idata", idata);
         if (d_total_patch_entries > 0) {
            restart_db->putDoubleVector("ddata", ddata);
         }
      }
   }
}

void
Statistic::getFromRestart(
   const std::shared_ptr<Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   int ver = restart_db->getInteger("TBOX_STATISTIC_VERSION");
   if (ver != TBOX_STATISTIC_VERSION) {
      TBOX_ERROR("Restart file version different than class version.");
   }

   d_object_name = restart_db->getString("object_name");
   d_stat_type = restart_db->getInteger("stat_type");
   d_instance_id = restart_db->getInteger("instance_id");
   int seq_entries = restart_db->getInteger("seq_counter");
   int total_patches = restart_db->getInteger("total_patch_entries");
   d_proc_stat_array_size = restart_db->getInteger("proc_stat_array_size");
   d_patch_stat_array_size = restart_db->getInteger("patch_stat_array_size");

   d_proc_array.resize(d_proc_stat_array_size);
   d_patch_array.resize(d_patch_stat_array_size);

   int i;
   if (d_stat_type == PROC_STAT) {
      if (seq_entries > 0) {
         std::vector<double> ddata = restart_db->getDoubleVector("ddata");
         for (i = 0; i < seq_entries; ++i) {
            recordProcStat(ddata[i], i);
         }
      }
   }

   if (d_stat_type == PATCH_STAT) {
      if (seq_entries > 0) {
         std::vector<int> idata = restart_db->getIntegerVector("idata");

         std::vector<int> inum_patches(seq_entries);
         for (i = 0; i < seq_entries; ++i) {
            inum_patches[i] = idata[i];
         }

         if (total_patches > 0) {
            std::vector<double> ddata = restart_db->getDoubleVector("ddata");

            int il = 0;
            int mark = seq_entries;
            for (i = 0; i < seq_entries; ++i) {
               for (int ipsr = 0; ipsr < inum_patches[i]; ++ipsr) {
                  int patch_num = idata[mark + il];
                  double val = ddata[il];
                  recordPatchStat(patch_num, val, i);
                  ++il;
               }
            }
         }
      }
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
