/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#ifdef HAVE_HDF5

#include "SAMRAI/tbox/HDFDatabase.h"

#include <string>

using namespace SAMRAI;

/*
 * RedistributedRestartUtility is a utility class used to build the
 * restart-redistribute tool.  All methods are static and the class
 * contains no data.
 */

class RedistributedRestartUtility
{

public:
/*
 * Write redistributed restart files to a new restart directory.
 *
 * @param output_dirname      name of directory containing new restart files
 * @param input_dirname       name of directory containing files to be read
 * @param total_input_files   number of input files being read by the tool
 * @param total_output_files  number of output files being created by the tool
 * @param file_mapping        mapping between input and output files
 * @param restore_num         number identifying the restart dump being
 *                            processed
 */
   static void
   writeRedistributedRestartFiles(
      const std::string& output_dirname,
      const std::string& input_dirname,
      const int total_input_files,
      const int total_output_files,
      const std::vector<std::vector<int> >& file_mapping,
      const int restore_num);

private:
/*
 * A recursive function that does the reading and writing of data.  If the
 * key represents any type other than Database, data is read from input and
 * immediately written to output.  If it is a Database, this function is
 * called recursively on every key in that Database.  If the Database
 * represents a PatchLevel, readAndWritePatchLevelRestartData is called
 * instead, in order to properly distribute patches among the output files.
 *
 * The optional arguments are not used when this function is called for
 * reading and writing data from a patch, because the data for a patch
 * exists in only one input file and is written to only one output file.
 * Meaningful values should be given for all of the optional arguments in
 * all other cases.
 */
   static void
   readAndWriteRestartData(
      std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
      const std::vector<std::shared_ptr<tbox::Database> >& input_dbs,
      const std::string& key,
      const std::vector<std::vector<int> >* file_mapping = 0,
      int num_files_written = -1,
      int which_file_mapping = -1,
      int total_input_files = -1,
      int total_output_files = -1);

/*
 * Reads and writes data in a Database that represents a PatchLevel.
 * Global data are written to every output database, and data representing
 * patches are written to only one output database.
 */
   static void
   readAndWritePatchLevelRestartData(
      std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
      const std::vector<std::shared_ptr<tbox::Database> >& level_in_dbs,
      const std::string& key,
      const int num_files_written,
      const std::vector<int>& input_proc_nums,
      const int total_output_files);

/*
 * Reads and writes data in a Database that represents a BoxLevel.
 * Global data are written to every output database, and data representing
 * Boxes are written to only one output database.
 */
   static void
   readAndWriteBoxLevelRestartData(
      std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
      const std::vector<std::shared_ptr<tbox::Database> >& level_in_dbs,
      const std::string& key,
      const int num_files_written,
      const std::vector<int>& input_proc_nums,
      const int total_output_files);

/*
 * Reads patch restart data from input database and writes it to output
 * database.
 */
   static void
   readAndWritePatchRestartData(
      std::shared_ptr<tbox::Database>& patch_out_db,
      const std::shared_ptr<tbox::Database>& patch_in_db,
      const int output_proc);

};

#endif
