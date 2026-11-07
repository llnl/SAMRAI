/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_tbox_Statistic
#define included_tbox_Statistic

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <list>

namespace SAMRAI {
namespace tbox {

class Statistician;

/**
 * Class Statistic defines a simple object that can be used to record
 * information generated during the course of a simulation for post-
 * processing later.  Each statistic object is created by the singleton
 * Statistian object and is defined by a name string
 * identifier and is characterized by the sort of information it may record.
 * Depending on how the object is created, it may record processor
 * information (i.e., a separate value for each processor), or patch
 * information (i.e., a separate value for each patch on each processor).
 * An example of the former may be the total number of cells on each
 * processor.  An example of the second may be the number of cells on each
 * patch.  Each recorded data item may be any numerical value, but it
 * will always be stored as a double for simplicity.  The string identifier
 * for a processor stat is "PROC_STAT" and the string identifier for a
 * patch stat is "PATCH_STAT".
 *
 * An example use of a Statistic to record the number of gridcells on each
 * processor is as follows:
 *
 *    std::shared_ptr<Statistic> stat_num_gridcells =
 *        Statistician::getStatistician()->
 *        getStatistic("NumberGridcells", "PROC_STAT");
 *    ...
 *    stat_num_gridcells->recordProcStat(num_cells_on_proc);
 *    ...
 *
 * The type of the statistic restricts the way in which the statistic
 * object may be used to record information.  For a "processor" stat
 * only the recordProcStat() functions can be used.  For a "patch"
 * stat only the recordPatchStat() functions can be used.
 *
 * Typically, the information is recorded to generate a time sequence of
 * values.  But this need not be the case always.  An optional time
 * stamp may be provided for each value as it is recorded.  In any case,
 * the sequence order of the values is determined by the recording order.
 *
 * Also, the Statistician class is used to manage Statistic
 * objects.  It provided a global point of access for creating and accessing
 * statistic objects and supports post-processing statistic information
 * in parallel.
 *
 * In some cases, it may be desirable to record information for each
 * level in a calculation; e.g., the number of cells on each processor
 * on level zero, level 1, etc.  In this case, one can cimply create a
 * separate statistic object for each level.
 *
 * @see Statistician
 */

class Statistic
{
   friend class Statistician;
public:
   /**
    * Virtual destructor destroys recorded object data.
    */
   ~Statistic();

   /**
    * Return string name identifier for statistic object.
    */
   const std::string&
   getName() const
   {
      return d_object_name;
   }

   /**
    * Return string statistic type identifier for statistic object.
    */
   std::string
   getType() const
   {
      return (d_stat_type == PROC_STAT) ? "PROC_STAT" : "PATCH_STAT";
   }

   /**
    * Return integer instance identifier for statistic object.
    */
   int
   getInstanceId() const
   {
      return d_instance_id;
   }

   /**
    * Return integer length of list of statistic sequence records.
    * This value is either the length of the processor statistic list
    * or the patch statistic list, whichever corresponds to the statistic
    * type.
    */
   int
   getStatSequenceLength() const
   {
      return d_seq_counter;
   }

   /**
    * Reset the state of the statistic information.
    */
   void
   reset()
   {
      d_proc_array.clear();
      d_patch_array.clear();
   }

   /**
    * Record double processor statistic value. The optional sequence number
    * argument identifies where in timestep sequence the value should be.
    * If the sequence number is not specified, an internal counter will
    * determine the appropriate sequence number.
    *
    * @pre getType() == "PROC_STAT"
    */
   void
   recordProcStat(
      double value,
      int seq_num = -1);

   /**
    * Record double patch statistic value.  The patch number refers to
    * the global patch number on a level.  The sequence number
    * argument identifies where in timestep sequence the value should be.
    * The sequence number MUST be explicitly specified because the number
    * of patches on each processor will generally be different at
    * each sequence step.
    *
    * @pre getType() == "PATCH_STAT"
    */
   void
   recordPatchStat(
      int patch_num,
      double value,
      int seq_num);

   /**
    * Return true if size of stream required to pack all statistic
    * data can be determined for all processors without exchanging
    * any details of structure of statistic data.  Otherwise, return false.
    */
   bool
   canEstimateDataStreamSize()
   {
      return false;
   }

   /**
    * Return number of bytes needed to stream the statistic data.
    * This is the amount needed by the stat transaction class.
    */
   size_t
   getDataStreamSize();

   /**
    * Pack contents of statistic data structure into message stream.
    *
    * @pre SAMRAI_MPI::getSAMRAIWorld().getRank() != 0
    */
   void
   packStream(
      MessageStream& stream);

   /**
    * Unpack contents of statistic data structure from message stream.
    *
    * @pre SAMRAI_MPI::getSAMRAIWorld().getRank() == 0
    */
   void
   unpackStream(
      MessageStream& stream);

   /**
    * Print statistic data to given output stream.  Floating point precision
    * can be specified (default is 12).
    *
    * @pre precision > 0
    */
   void
   printClassData(
      std::ostream& stream,
      int precision = 12) const;

   /**
    * Write statistic data members to restart database. The restart_db pointer
    * must be non-null.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<Database>& restart_db) const;

   /**
    * Read restarted times from restart database.
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<Database>& restart_db);

   /*
    * These structures are used to store statistic data entries.
    * They need to be declared public for the Sun CC compiler.
    */
   struct ProcStat {
      double value;        // stat record value
   };

   struct PatchStatRecord {
      int patch_id;         // global patch number
      double value;         // stat record value
   };

   struct PatchStat {
      std::list<Statistic::PatchStatRecord> patch_records; // stat record
   };

protected:
   /**
    * The constructor for the Statistic class sets the name string
    * and the statistic type for a statistic object.
    *
    * @pre !name.empty()
    * @pre !stat_type.empty()
    * @pre (stat_type == "PROC_STAT") || (stat_type = "PATCH_STAT")
    * @pre instance_id > -1
    */
   Statistic(
      const std::string& name,
      const std::string& stat_type,
      int instance_id);

   /**
    * Return const reference to list of processor records.
    */
   const std::vector<Statistic::ProcStat>&
   getProcStatSeqArray() const
   {
      return d_proc_array;
   }

   /**
    * Return const reference to list of patch records.
    */
   const std::vector<Statistic::PatchStat>&
   getPatchStatSeqArray() const
   {
      return d_patch_array;
   }

private:
   /*
    * Static double value used to indicate when a particular sequence entry
    * is skipped.
    */
   static double s_empty_seq_tag_entry;

   /*
    * Static integer constant describing this class's version number.
    */
   static const int TBOX_STATISTIC_VERSION;

   static const int ARRAY_INCREMENT;

   /*
    * Assess whether the processor or patch stat arrays need to be resized.
    */
   void
   checkArraySizes(
      int seq_num);

   // The following three members are not implemented
   Statistic();
   Statistic(
      const Statistic&);
   Statistic&
   operator = (
      const Statistic&);

   /*
    * The enumerated type maps statistic types to integers.
    */
   enum STATISTIC_RECORD_TYPE { PROC_STAT = 0, PATCH_STAT = 1 };

   /*
    * Name, instance id, and type identifier for this statistic object.
    */
   std::string d_object_name;
   int d_instance_id;
   int d_stat_type;            // see STATISTIC_RECORD_TYPE above.

   /*
    * Vectors of records.  Note that one of these will always be empty.
    * Integer sequence length refers to length of list corresponding
    * to stat type.
    */
   std::vector<Statistic::ProcStat> d_proc_array;
   std::vector<Statistic::PatchStat> d_patch_array;

   /*
    * Sequence and patch counters (NOTE: patch counter use for patch stats
    * only) and high-water-mark array sizes for proc and patch stats.
    */
   int d_seq_counter;
   int d_total_patch_entries;
   int d_proc_stat_array_size;
   int d_patch_stat_array_size;
};

}
}

#endif
