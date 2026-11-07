/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Singleton manager class for statistic objects.
 *
 ************************************************************************/

#ifndef included_tbox_Statistician
#define included_tbox_Statistician

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Statistic.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace tbox {

class StatisticRestartDatabase;

/**
 * Class Statistician is a Singleton class that manages a simple
 * database of Statistic objects.  This class provides a single point
 * of access to statistic objects so that any one of them may be updated
 * or recorded at any point in the code.  Access to the Singleton
 * statistician instance follows the standard SAMRAI implementation
 * found in other classes with similar Singleton behavior.  See static
 * member functions below for more information.
 *
 * Statistic objects can be to the database or accessed in code as follows:
 *
 *     std::shared_ptr<Statistic> stat =
 *           Statistician::getStatistician->
 *           getStatistic("name", "PROC_STAT");
 *
 * Here `name' is the name string identifier for the statistic object and
 * `PROC_STAT' is the type of statistic. See discussion for the getStatistic()
 * method below for more information.
 *
 * The statistic state is saved to restart files when restart file generation
 * is active.  This allows users to continue to accumulate timing
 * information when restarting a run.  If desired, all statistics can be
 * reset when restarting by calling the function resetAllStatistics().
 *
 * A variety of print options exist to dump statistics data.  Notably,
 * the printSpreadSheetOutput("print_dir") will write statistics data in
 * a tab-separated format to files in the supplied directory name.  The
 * naming convention for statistics data is "\<name\>-\<type\>.txt" where \<name\>
 * is the name of the statistic and \<type\> is either proc or patch stat.  The
 * files may be read in to a spreadsheet program such as MS Excel.
 *
 * For more information about data that can be recorded with statistics,
 * consult the header file for the Statistic class.
 *
 * @see Statistic
 */

class Statistician
{
   friend class StatisticRestartDatabase;
public:
   /**
    * Create the singleton instance of the statistic manager and return
    * a pointer to it.  This function is provided so that so that
    * users can control whether statistic information will be written
    * to/read from restart files.
    *
    * Statistics that exist in the restart file will be read from restart
    * when a run is restarted and the second argument is true.
    *
    * Generally, this routine should only be called once during program
    * execution.  If the statistician has been previously created (e.g.,
    * by an earlier call to this routine) this routine will do nothing
    * other than return the pointer to the existing singleton instance.
    */
   static Statistician *
   createStatistician(
      bool read_from_restart = true);

   /**
    * Return a pointer to the singleton statistician instance.
    * All access to the Statistician object is through the
    * getStatistician() function.  For example, to add a statistic
    * object with the name "my_stat" to the statistician, use the
    * following call:
    * Statistician::getStatistician()->addStatistic("my_stat").
    *
    * @pre s_statistician_instance
    */
   static Statistician *
   getStatistician();

   /**
    * Return pointer to statistic object with the given name string.
    * If a statistics with the given name already exists in the database
    * of statistics, the statistic with that name will be returned.
    * Otherwise, a new statistic will be created with that name.  The
    * stat_type string identifier is only used when a new statistic
    * object must be created.  The two avaible options are processor
    * statistics and patch statistics which are indicated by the strings
    * "PROC_STAT" and "PATCH_STAT", respectively.
    *
    * @pre !name.empty()
    * @pre !stat_type.empty()
    * @pre (stat_type == "PROC_STAT") || (stat_type == "PATCH_STAT")
    */
   std::shared_ptr<Statistic>
   getStatistic(
      const std::string& name,
      const std::string& stat_type);

   /**
    * Return true if a statistic whose name matches the argument string
    * exists in the database of statistics controlled by the statistician.
    * If a match is found, the statistic pointer in the argument list is set
    * to that statistic.  Otherwise, return false and return a null pointer.
    *
    * @pre !name.empty()
    */
   bool
   checkStatisticExists(
      std::shared_ptr<Statistic>& stat,
      const std::string& name) const;

   /**
    * Return integer number of local processor statistics maintained
    * by statistician.
    */
   int
   getNumberProcessorStats() const
   {
      return d_num_proc_stats;
   }

   /**
    * Return integer number of local patch statistics maintained
    * by statistician.
    */
   int
   getNumberPatchStats() const
   {
      return d_num_patch_stats;
   }

   /**
    * Reset all processor statistics to contain no information. The primary
    * intent of this function is to avoid using restarted statistic values
    * when performing a restarted run.  However, it can be called anytime.
    */
   void
   resetProcessorStatistics();

   /**
    * Reset all patch statistics to contain no information. The primary
    * intent of this function is to avoid using restarted statistic values
    * when performing a restarted run.  However, it can be called anytime.
    */
   void
   resetPatchStatistics();

   /**
    * Reset all patch and processor statistics to contain no information.
    */
   void
   resetStatistics();

   /**
    * Return integer instance identifier for processor statistic with given
    * name.  If this statistician object maintains no processor statistic
    * with that name, then a warning message results and the return value
    * will be the invalid instance identifier "-1".
    */
   int
   getProcStatId(
      const std::string& name) const;

   /**
    * Return number of sequence entries for processor statistic with given
    * integer identifier.   For convenience, the routine getProcStatId() is
    * provided to map the statistic string name to the proper integer
    * identifier.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_state &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_data.size())))
    */
   int
   getGlobalProcStatSequenceLength(
      int proc_stat_id);

   /**
    * Return statistic data value for processor statistic with given integer
    * identifier, sequence number, and processor number.   For convenience,
    * the routine getProcStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function
    * getGlobalProcStatSequenceLength() provides the sequence length for
    * a given processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_state &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_data[proc_stat_id].size())) &&
    *       (proc_num < SAMRAI_MPI::getSAMRAIWorld().getSize()))
    */
   double
   getGlobalProcStatValue(
      int proc_stat_id,
      int seq_num,
      int proc_num);

   /**
    * Return global sum of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_sum.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_sum[proc_stat_id].size())))
    */
   double
   getGlobalProcStatSum(
      int proc_stat_id,
      int seq_num);

   /**
    * Return global max of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_max.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_max[proc_stat_id].size())))
    */
   double
   getGlobalProcStatMax(
      int proc_stat_id,
      int seq_num);

   /**
    * Returns rank of processor holding global max for the processor
    * statistic specified by the given integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_imax.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_imax[proc_stat_id].size())))
    */
   int
   getGlobalProcStatMaxProcessorId(
      int proc_stat_id,
      int seq_num);

   /**
    * Return global min of processor statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getProcStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalProcStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_min.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_min[proc_stat_id].size())))
    */
   double
   getGlobalProcStatMin(
      int proc_stat_id,
      int seq_num);

   /**
    * Returns rank of processor holding global max for the processor
    * statistic specified by the given integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize &&
    *       (proc_stat_id >= 0) &&
    *       (proc_stat_id < static_cast<int>(d_global_proc_stat_imin.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_proc_stat_imin[proc_stat_id].size())))
    */
   int
   getGlobalProcStatMinProcessorId(
      int proc_stat_id,
      int seq_num);

   /**
    * Print global processor statistic data for a particular statistic
    * to given output stream.  Floating point precision may be specified
    * (default is 12).  Note that this method generates a general dump of
    * the data but does NOT generate it in tabulated form.  To generate
    * tabulated data, see the printGlobalProcStatDataFormatted() method.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (proc_stat_id >= 0) && (precision > 0))
    */
   void
   printGlobalProcStatData(
      int proc_stat_id,
      std::ostream& os,
      int precision = 12);

   /**
    * Print processor stat data in formatted output to given output
    * stream.  Floating point precision may be specified (default is 12).
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (proc_stat_id >= 0) && (precision > 0))
    */
   void
   printGlobalProcStatDataFormatted(
      int proc_stat_id,
      std::ostream& os,
      int precision = 12);

   /**
    * Print stat data for specified processor in formatted output to
    * given output stream.  Floating point precision may be specified
    * (default is 12).
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (proc_stat_id >= 0) && (proc_id >= 0) && (precision > 0))
    */
   void
   printGlobalProcStatDataFormatted(
      int proc_stat_id,
      int proc_id,
      std::ostream& os,
      int precision = 12);

   /**
    * Return integer instance identifier for patch statistic with given
    * name.  If this statistician object maintains no patch statistic
    * with that name, then a warning message results and the return value
    * will be the invalid instance identifier "-1".
    */
   int
   getPatchStatId(
      const std::string& name) const;

   /**
    * Return number of sequence entries for patch statistic with given
    * integer identifier.   For convenience, the routine getPatchStatId()
    * is provided to map the statistic string name to the proper integer
    * identifier.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())))
    */
   int
   getGlobalPatchStatSequenceLength(
      int patch_stat_id);

   /**
    * Return number of patch entries for patch statistic with given
    * integer identifier, and sequence number.   For convenience, the
    * routine getPatchStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatNumberPatches(
      int patch_stat_id,
      int seq_num);

   /**
    * Return global processor mapping for patch statistic with given integer
    * identifier, sequence number, and patch number.  For convenience,
    * the routine getPatchStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.  The function
    * getGlobalPatchStatNumberPatches() gives the number of patches
    * associated with a patch statistic and sequence number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_mapping.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_mapping[patch_stat_id].size())) &&
    *       (patch_num >= 0) &&
    *       (patch_num < static_cast<int>(d_global_patch_stat_mapping[patch_stat_id][seq_num].size())))
    */
   int
   getGlobalPatchStatPatchMapping(
      int patch_stat_id,
      int seq_num,
      int patch_num);

   /**
    * Return statistic data value for patch statistic with given integer
    * identifier, sequence number, and patch number.   For convenience,
    * the routine getPatchStatId() is provided to map the statistic string
    * name to the proper integer identifier.  The function
    * getGlobalPatchStatSequenceLength() provides the sequence length for
    * a given patch statistic.  The function
    * getGlobalPatchStatNumberPatches() gives the number of patches
    * associated with a patch statistic and sequence number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())) &&
    *       (patch_num >= 0) &&
    *       (patch_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id][seq_num].size())))
    */
   double
   getGlobalPatchStatValue(
      int patch_stat_id,
      int seq_num,
      int patch_num);

   /**
    * Return global sum of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatSum(
      int patch_stat_id,
      int seq_num);

   /**
    * Return global max of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatMax(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns ID of patch holding global max for the patch
    * statistic specified by the given integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMaxPatchId(
      int patch_stat_id,
      int seq_num);

   /**
    * Return global min of patch statistic with given integer
    * identifier and sequence number.   To identify the correct integer
    * identifier and valid sequence numbers, the method getPatchStatId() maps
    * the statistic string name to its integer identifier and the method
    * getGlobalPatchStatSequenceLength() returns the maximum sequence length
    * for the processor statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatMin(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns patch ID of patch holding global min for the patch
    * statistic specified by the given integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMinPatchId(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the sum of patch statistic information for a particular
    * processor.  The patch statistic is specified by its integer identifyer
    * and sequence number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (processor_id >= 0) &&
    *       (processor_id < SAMRAI_MPI::getSAMRAIWorld().getSize()) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatProcessorSum(
      int patch_stat_id,
      int processor_id,
      int seq_num);

   /**
    * Returns the maximum value of the patch statistic data summed
    * on each processor.  That is, patch statistic information is
    * summed on each processor, and this method returns the maximum
    * value, across all processors, of this summed data. The patch
    * statistic is specified by its integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatProcessorSumMax(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the processor ID which holds the maximum value of
    * summed patch statistic information across processors.  See
    * the discussion for the method getGlobalPatchStatProcessorSumMax()
    * for more information on the summed patch statistic information
    * on processors.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatProcessorSumMaxId(
      int patch_stat_id,
      int seq_num);
   /**
    * Returns the minimum value of the patch statistic data summed
    * on each processor.  That is, patch statistic information is
    * summed on each processor, and this method returns the minimum
    * value, across all processors, of this summed data. The patch
    * statistic is specified by its integer identifyer and sequence
    * number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   double
   getGlobalPatchStatProcessorSumMin(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the processor ID which holds the minimum value of
    * summed patch statistic information across processors.  See
    * the discussion for the method getGlobalPatchStatProcessorSumMin()
    * for more information on the summed patch statistic information
    * on processors.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatProcessorSumMinId(
      int patch_stat_id,
      int seq_num);

   /**
    * Return number of patches on the specified processor number for
    * patch statistic with given identifier, and sequence number.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())) &&
    *       (proc_id >= 0 && proc_id < SAMRAI_MPI::getSAMRAIWorld().getSize()))
    */
   int
   getGlobalPatchStatNumberPatchesOnProc(
      int patch_stat_id,
      int seq_num,
      int proc_id);

   /**
    * Returns the maximum number of patches per processor for the
    * specified patch statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMaxPatchesPerProc(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the processor ID holding the maximum number of patches
    * per processor for the specified patch statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMaxPatchesPerProcId(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the minimum number of patches per processor for the
    * specified patch statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())) &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMinPatchesPerProc(
      int patch_stat_id,
      int seq_num);

   /**
    * Returns the processor ID holding the minimum number of patches
    * per processor for the specified patch statistic.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats &&
    *       (patch_stat_id >= 0) &&
    *       (patch_stat_id < static_cast<int>(d_global_patch_stat_proc_data.size())) &&
    *       (seq_num >= 0) &&
    *       (seq_num < static_cast<int>(d_global_patch_stat_proc_data[patch_stat_id].size())))
    */
   int
   getGlobalPatchStatMinPatchesPerProcId(
      int patch_stat_id,
      int seq_num);
   /**
    * Print global processor statistic data for a particular statistic
    * to given output stream.  Floating point precision may be specified
    * (default is 12).  Note that this method generates a general dump of
    * the data but does NOT generate it in tabulated form.  To generate
    * tabulated data, see the printGlobalPatchStatDataFormatted() method.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (patch_stat_id >= 0 && precision > 0)
    */
   void
   printGlobalPatchStatData(
      int patch_stat_id,
      std::ostream& os,
      int precision = 12);

   /**
    * Print patch stat data in formatted output to given output
    * stream.  Floating point precision may be specified (default is 12).
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (patch_stat_id >= 0 && precision > 0)
    */
   void
   printGlobalPatchStatDataFormatted(
      int patch_stat_id,
      std::ostream& os,
      int precision = 12);

   /**
    * Typically, this routine is called after
    * a calculation has completed so that statistic data can be
    * retrieved, analized, printed to a file, etc.  It is not essential
    * that this routine be called, however, as each "get" and "print"
    * routine checks to see if statistic data has been finalized before
    * it peforms its function.
    *
    * If gather_individual_stats_on_proc_0 == true, the statistics are
    * gathered on proc 0 for further access.  If not, only
    * globally-reduced (min, max, sum) values for "PROC_STATS" are
    * available (on all processes).  Non-local "PATCH_STATS" are only available
    * if gather_individual_stats_on_proc_0 == true.
    */
   void
   finalize(
      bool gather_individual_stats_on_proc_0 = false);

   /**
    * Print data to given output stream for local statistics managed
    * by this statistician object.  Note that no fancy formatting is done.
    * Floating point precision can be specified (default is 12).
    */
   void
   printLocalStatData(
      std::ostream& os,
      int precision = 12) const;

   /**
    * Print global statistic data information to given output stream.
    * The data will NOT be in tabulated form.  Floating point precision
    * can be specified (default is 12).
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats)
    */
   void
   printAllGlobalStatData(
      std::ostream& os,
      int precision = 12);

   /**
    * Print sums of all global statistic data information to given
    * output stream. Floating point precision can be specified (default is 12).
    */
   void
   printAllSummedGlobalStatData(
      std::ostream& os,
      int precision = 12);

   /**
    * Print sums of all global statistic data information to specified
    * filename. Floating point precision can be specified (default is 12).
    */
   void
   printAllSummedGlobalStatData(
      const std::string& filename,
      int precision = 12);

   /**
    * Write all statistics data in tab-separated format to files in the
    * supplied directory name.  The naming convention used is "\<name\>-\<type\>.txt"
    * where \<name\> is the name of the statistic and \<type\> is either proc or
    * patch stat.  Floating point precision may be specified (default is 12).
    * The files may be read in to a spreadsheet program such as MS Excel.  If
    * no directory name is supplied, the files will be written to the directory
    * where the application is run.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats)
    */
   void
   printSpreadSheetOutput(
      const std::string& dirname = std::string(),
      int precision = 12);

   /**
    * Write tab-separated statistics data for specified processor.  This
    * method is identical to "printSpreadSheetOutput()" (above), but only
    * prints data for a single processor.  This may be useful for information
    * that is the same across all processors.  This method will only print
    * processor stats. Any patch stats will be ignored.
    *
    * @pre (SAMRAI_MPI::getSAMRAIWorld().getRank() != 0) ||
    *      (!d_must_call_finalize && d_has_gathered_stats)
    */
   void
   printSpreadSheetOutputForProcessor(
      const int proc_id,
      const std::string& dirname = std::string(),
      int precision = 12);

protected:
   /**
    * The constructor for Statistician is protected.  Consistent
    * with the definition of a Singleton class, only a statistician object
    * can have access to the constructor for the class.
    */
   Statistician();

   /**
    * Statistician is a Singleton class; its destructor is protected.
    */
   ~Statistician();

   /**
    * Initialize Singleton instance with instance of subclass.  This function
    * is used to make the singleton object unique when inheriting from this
    * base class.
    *
    * @pre !s_statistician_instance
    */
   void
   registerSingletonSubclassInstance(
      Statistician * subclass_instance);

   /**
    * During finalize() check statistic information on all processors
    * for consistency before generating vectors of data.
    *
    * @pre total_patches.size() == 0
    */
   void
   checkStatsForConsistency(
      std::vector<int>& total_patches);

   /**
    * Return true if a processor statistic whose name matches the
    * argument string exists in the database of statistics controlled
    * by the statistician.  If a match is found, the statistic pointer
    * in the argument list is set to that statistic.  Otherwise, return
    * false and return a null pointer.
    */
   bool
   checkProcStatExists(
      std::shared_ptr<Statistic>& stat,
      const std::string& name) const;

   /**
    * Return true if a patch statistic whose name matches the
    * argument string exists in the database of statistics controlled
    * by the statistician.  If a match is found, the statistic pointer
    * in the argument list is set to that statistic.  Otherwise, return
    * false and return a null pointer.
    */
   bool
   checkPatchStatExists(
      std::shared_ptr<Statistic>& stat,
      const std::string& name) const;

private:
   // Unimplemented copy constructor.
   Statistician(
      const Statistician& other);

   // Unimplemented assignment operator.
   Statistician&
   operator = (
      const Statistician& rhs);

   /*!
    * @brief Get global-reduction statistics without depending on an
    * MPI gather, which is slow and does not scale.
    */
   void
   reduceGlobalStatistics();

   /*
    * Gets the current maximum number of statistics.
    *
    * If trying to use more statistics than this value
    * the vectors should be resized.
    */
   int
   getMaximumNumberOfStatistics()
   {
      return static_cast<int>(d_proc_statistics.size());
   }

   /*
    * Set the maximum number of statistics.
    *
    * This will grow the internal vectors used to store values.
    */
   void
   setMaximumNumberOfStatistics(
      const int size)
   {
      if (size > static_cast<int>(d_proc_statistics.size())) {
         d_proc_statistics.resize(size);
         d_patch_statistics.resize(size);
      }
   }

   /**
    * Static data members to manage the singleton statistician instance.
    */
   static Statistician* s_statistician_instance;

   static void
   makeStatisticianInstance(
      bool read_from_restart = true);

   /*
    * Create and initialize state of restart database.
    */
   void
   initRestartDatabase(
      bool read_from_restart);

   /**
    * Allocate the Statistician instance.
    *
    * Automatically called by the StartupShutdownManager class.
    */
   static void
   initializeCallback();

   /**
    * Shutdown Statistician instance.
    *
    * Automatically called by the StartupShutdownManager class.
    */
   static void
   shutdownCallback();

   /**
    * Deallocate the Statistician instance.
    *
    * Automatically called by the StartupShutdownManager class.
    */
   static void
   finalizeCallback();

   /*
    * Internal database class for statistician restart capabilities.  See
    * class declaration below.
    */
   StatisticRestartDatabase* d_restart_database_instance;

   /*
    * Count of statistics registered with the statistician and vectors of
    * pointers to those statistics.
    */
   int d_num_proc_stats;
   std::vector<std::shared_ptr<Statistic> > d_proc_statistics;
   int d_num_patch_stats;
   std::vector<std::shared_ptr<Statistic> > d_patch_statistics;

   /*
    * Vectors of global statistic data assembled by the finalize() function.
    *
    * Global processor stat data is assembled as
    *    vector(stat id, seq id, proc id) = proc stat value.
    *
    * Global patch stat data is assembled as
    *    vector(stat_id, seq id, global patch id) = patch stat value.
    *
    * Global patch stat processor data is assembled as
    *    vector(stat_id, seq id, global proc id) = patch stats summed on
    *    different processors.
    *
    * The map of patches to processors is assembled as
    *    vector(stat_id, seq id, global patch id) = proc number.
    */
   bool d_must_call_finalize;
   bool d_has_gathered_stats;

   std::vector<std::vector<std::vector<double> > > d_global_proc_stat_data;

   std::vector<std::vector<std::vector<double> > > d_global_patch_stat_data;
   std::vector<std::vector<std::vector<double> > > d_global_patch_stat_proc_data;
   std::vector<std::vector<std::vector<int> > > d_global_patch_stat_mapping;

   /*!
    * @brief Vector of max-reduced processor stat data.
    *
    * d_global_proc_stat_max[i][j] is the max over all processors of
    * the stat id (i) and sequence id (j).
    */
   std::vector<std::vector<double> > d_global_proc_stat_max;

   /*!
    * @brief Processor owning the max value of processor stat data.
    *
    * d_global_proc_stat_imax[i][j] is the process corresponding to
    * d_global_proc_stat_max[i][j].
    */
   std::vector<std::vector<int> > d_global_proc_stat_imax;

   /*!
    * @brief Vector of min-reduced processor stat data.
    *
    * d_global_proc_stat_min[i][j] is the min over all processors of
    * the stat id (i) and sequence id (j).
    */
   std::vector<std::vector<double> > d_global_proc_stat_min;

   /*!
    * @brief Processor owning the min value of processor stat data.
    *
    * d_global_proc_stat_imin[i][j] is the process corresponding to
    * d_global_proc_stat_max[i][j].
    */
   std::vector<std::vector<int> > d_global_proc_stat_imin;

   /*!
    * @brief Vector of sum-reduced processor stat data.
    *
    * d_global_proc_stat_sum[i][j] is the sum over all processors of
    * the stat id (i) and sequence id (j).
    */
   std::vector<std::vector<double> > d_global_proc_stat_sum;

   /*
    * Internal value used to set and grow vectors for storing
    * statistics.
    */
   static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT;

   static StartupShutdownManager::Handler s_finalize_handler;
};

/*
 * Class StatisticRestartDatabase is a separate class used by the
 * statistician to provide restart capabilities. Each restartable
 * class must be derived from Serializable.  Since Statistician
 * is a Singleton, its destructor is protected.  StatisticRestartDatabase
 * has a publically accessible destructor.  To avoid improper use of this
 * class, it is privately derived from Serializable, and make it a
 * friend of Statistician.  In this way, its methods can only
 * be accessed via Serializable and Statistician.
 */

class StatisticRestartDatabase:private Serializable
{
   friend class Statistician;
public:
   /*
    * The StatisticRestartDatabase constructor caches a copy of the
    * database object name and registers the object with the restart
    * manager for subsequent restart files.  If the run is started from
    * a restart file and the boolean argument is true, we initialize
    * the statistics from restart.
    *
    * @pre !object_name.empty()
    */
   StatisticRestartDatabase(
      const std::string& object_name,
      bool read_from_restart);

   /*
    * The destructor for StatisticRestartDatabase unregisters
    * the database object with the restart manager.
    */
   virtual ~StatisticRestartDatabase();

   /*
    * Put all statistics and their state in the given restart database.
    * This function is inherited from Serializable.
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<Database>& restart_db) const;

   /*
    * Construct those statistics saved in the restart database.
    */
   void
   getFromRestart();

   /**
    * Return string name identifier for statistic object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   // Unimplemented default constructor.
   StatisticRestartDatabase();

   // Unimplemented default constructor.
   StatisticRestartDatabase(
      const StatisticRestartDatabase& other);

   // Unimplemented assignment operator.
   StatisticRestartDatabase&
   operator = (
      const StatisticRestartDatabase& rhs);

   std::string d_object_name;

   /*
    * Static integer constant describing this class's version number.
    */
   static const int TBOX_STATISTICRESTARTDATABASE_VERSION;

};

}
}

#endif
