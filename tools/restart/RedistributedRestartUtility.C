/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include "RedistributedRestartUtility.h"

#ifdef HAVE_HDF5

#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cassert>
#include <list>

#define NAME_BUF_SIZE (32)

/*
 **************************************************************************
 * writeRedistributedRestartFiles
 **************************************************************************
 */

void RedistributedRestartUtility::writeRedistributedRestartFiles(
   const std::string& output_dirname,
   const std::string& input_dirname,
   const int total_input_files,
   const int total_output_files,
   const std::vector<std::vector<int> >& file_mapping,
   const int restore_num)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int nprocs = mpi.getSize();
   int rank = mpi.getRank();
   int output_files_per_proc = total_output_files / nprocs;
   int extra_output_files = total_output_files % nprocs;
   int num_files_written;
   if (total_input_files < total_output_files) {
      num_files_written = file_mapping[0][0];
   } else {
      num_files_written = output_files_per_proc * rank;
      if (extra_output_files) {
         if (rank >= extra_output_files) {
            num_files_written += extra_output_files;
         } else {
            num_files_written += rank;
         }
      }
   }
   int num_iterations = static_cast<int>(file_mapping.size());

   for (int icount = 0; icount < num_iterations; ++icount) {

      //We are writing to one file or reading only one file
      int num_files_to_read = (total_input_files < total_output_files) ?
         1 : static_cast<int>(file_mapping[icount].size());
      int num_files_to_write = (total_input_files < total_output_files) ?
         static_cast<int>(file_mapping[icount].size()) : 1;

      std::string restore_buf =
         "/restore." + tbox::Utilities::intToString(restore_num, 6);
      std::string nodes_buf =
         "/nodes." + tbox::Utilities::nodeToString(total_output_files);

      std::string restart_dirname = output_dirname + restore_buf + nodes_buf;

      //Make the subdirectories if this is the first iteration.
      if (icount == 0) {
         tbox::Utilities::recursiveMkdir(restart_dirname);
         tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
      }

      //Mount the output files on an array of output databases
      std::vector<std::shared_ptr<tbox::Database> >
      output_dbs(num_files_to_write);

      for (int i = 0; i < num_files_to_write; ++i) {

         int cur_out_file_id;
         if (total_input_files < total_output_files) {
            cur_out_file_id = file_mapping[icount][i];
         } else {
            cur_out_file_id = (output_files_per_proc * rank) + icount;
            if (extra_output_files) {
               if (rank >= extra_output_files) {
                  cur_out_file_id += extra_output_files;
               } else {
                  cur_out_file_id += rank;
               }
            }
         }

         std::string proc_buf =
            "/proc." + tbox::Utilities::processorToString(cur_out_file_id);

         std::string output_filename = restart_dirname + proc_buf;

         output_dbs[i].reset(new tbox::HDFDatabase(output_filename));

         int open_success = output_dbs[i]->create(output_filename);

         if (open_success < 0) {
            TBOX_ERROR(
               "Failed to open output file " << output_filename
                                             << "  HDF return code:  "
                                             << open_success);
         }

      }

      //Mount the input files on an array of input databases.
      std::vector<std::shared_ptr<tbox::Database> >
      input_dbs(num_files_to_read);

      nodes_buf = "/nodes." + tbox::Utilities::nodeToString(total_input_files);

      std::vector<std::string> input_keys(0);
      std::vector<std::string> test_keys(0);
      int num_keys = 0;

      int input_files_per_proc = total_input_files / nprocs;
      int extra_input_files = total_input_files % nprocs;
      for (int i = 0; i < num_files_to_read; ++i) {

         int cur_in_file_id;
         if (total_input_files < total_output_files) {
            cur_in_file_id = (input_files_per_proc * rank) + icount;
            if (extra_input_files) {
               if (rank >= extra_input_files) {
                  cur_in_file_id += extra_input_files;
               } else {
                  cur_in_file_id += rank;
               }
            }
         } else {
            cur_in_file_id = file_mapping[icount][i];
         }

         std::string proc_buf =
            "/proc." + tbox::Utilities::processorToString(cur_in_file_id);

         std::string restart_filename = input_dirname + restore_buf + nodes_buf
            + proc_buf;

         input_dbs[i].reset(new tbox::HDFDatabase(restart_filename));

         int open_success = input_dbs[i]->open(restart_filename);

         if (open_success < 0) {
            TBOX_ERROR(
               "Failed to open input file " << restart_filename
                                            << "  HDF return code:  "
                                            << open_success);
         }

         //Get the array of input keys.
         if (i == 0) {
            input_keys = input_dbs[i]->getAllKeys();
            num_keys = static_cast<int>(input_keys.size());
         } else {
            test_keys = input_dbs[i]->getAllKeys();
            if (static_cast<int>(test_keys.size()) != num_keys) {
               TBOX_ERROR("Input files contain differing number of keys");
            }
         }
      }

      //For every input key, call the recursive function that reads from the
      //input databases and writes to output databases.
      for (int i = 0; i < num_keys; ++i) {
         readAndWriteRestartData(output_dbs,
            input_dbs,
            input_keys[i],
            &file_mapping,
            num_files_written,
            icount,
            total_input_files,
            total_output_files);
      }

      //Unmount the databases.  This closes the files.
      for (int i = 0; i < num_files_to_read; ++i) {
         input_dbs[i]->close();
      }
      for (int i = 0; i < num_files_to_write; ++i) {
         output_dbs[i]->close();
      }

      num_files_written += num_files_to_write;
   }
}

/*
 **************************************************************************
 * readAndWriteRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWriteRestartData(
   std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
   const std::vector<std::shared_ptr<tbox::Database> >& input_dbs,
   const std::string& key,
   const std::vector<std::vector<int> >* file_mapping,  // = 0
   const int num_files_written, // = -1,
   const int which_file_mapping, // = -1
   const int total_input_files, // = -1,
   const int total_output_files) // = -1););
{
#ifdef DEBUG_CHECK_ASSERTIONS
   //One of the database arrays must be of size 1, and the other must be of
   //size >= 1.
   assert(output_dbs.size() >= 1);
   assert(input_dbs.size() >= 1);
   assert(input_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //This function works under the assumption that all of the input databases
   //contain the same keys, so we only need to check the type associated with
   //the key with one input database.

   //If the key is associated with any type other than Database, then the data
   //can be read from input_dbs[0] and written to every element of output_dbs.
   //The Database case is handled separately.

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int nprocs = mpi.getSize();
   int rank = mpi.getRank();

   if (input_dbs[0]->isDatabase(key)) {
      std::shared_ptr<tbox::Database> db = input_dbs[0]->getDatabase(key);

      if (db->keyExists("d_is_patch_level") &&
          db->getBool("d_is_patch_level")) {

         //Here we are handling the input database(s) for a PatchLevel.

         //Create array of level input databases.
         std::vector<std::shared_ptr<tbox::Database> > level_in_dbs(
            input_dbs.size());

         for (int i = 0; i < static_cast<int>(input_dbs.size()); ++i) {
            level_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         //input_proc_nums is an array that contains all of the processor
         //numbers that created the input databases that are currently
         //being processed.
         std::vector<int> input_proc_nums;
         if (total_input_files < total_output_files) {
            input_proc_nums.resize(1);
            int output_files_per_proc = total_output_files / nprocs;
            int extra_output_files = total_output_files % nprocs;
            int input_proc_num =
               (output_files_per_proc * rank) + which_file_mapping;
            if (extra_output_files) {
               if (rank >= extra_output_files) {
                  input_proc_num += extra_output_files;
               } else {
                  input_proc_num += rank;
               }
            }
            input_proc_nums[0] = input_proc_num;
         } else {
            input_proc_nums = (*file_mapping)[which_file_mapping];
         }

         //Call routine to write output according to the new processor mapping
         readAndWritePatchLevelRestartData(output_dbs,
            level_in_dbs,
            key,
            num_files_written,
            input_proc_nums,
            total_output_files);

      } else if (db->keyExists("d_is_mapped_box_level") &&
                 db->getBool("d_is_mapped_box_level")) {

         //Create array of level input databases.
         std::vector<std::shared_ptr<tbox::Database> > level_in_dbs(
            input_dbs.size());

         for (int i = 0; i < static_cast<int>(input_dbs.size()); ++i) {
            level_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         //input_proc_nums is an array that contains all of the processor
         //numbers that created the input databases that are currently
         //being processed.
         std::vector<int> input_proc_nums;
         if (total_input_files < total_output_files) {
            input_proc_nums.resize(1);
            int output_files_per_proc = total_output_files / nprocs;
            int extra_output_files = total_output_files % nprocs;
            int input_proc_num =
               (output_files_per_proc * rank) + which_file_mapping;
            if (extra_output_files) {
               if (rank >= extra_output_files) {
                  input_proc_num += extra_output_files;
               } else {
                  input_proc_num += rank;
               }
            }
            input_proc_nums[0] = input_proc_num;
         } else {
            input_proc_nums = (*file_mapping)[which_file_mapping];
         }

         readAndWriteBoxLevelRestartData(output_dbs,
            level_in_dbs,
            key,
            num_files_written,
            input_proc_nums,
            total_output_files);

      } else if (db->keyExists("d_is_edge_set") &&
                 db->getBool("d_is_edge_set")) {
         // don't write out edge sets.  They will be reconstructed from
         // PatchHierarchy during restarted run initialization.

         // no operation needed here.

      } else {
         //If this block is entered, then the key represents a database that
         //is not a patch level database.  We created child database arrays
         //for input_dbs and output_dbs, and then call readAndWriteRestartData
         //recursively.

         std::vector<std::shared_ptr<tbox::Database> > child_in_dbs(
            input_dbs.size());

         for (int i = 0; i < static_cast<int>(input_dbs.size()); ++i) {
            child_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         std::vector<std::shared_ptr<tbox::Database> > child_out_dbs(
            output_dbs.size());

         for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
            child_out_dbs[i] = output_dbs[i]->putDatabase(key);
         }

         std::vector<std::string> child_keys = db->getAllKeys();

         for (int j = 0; j < static_cast<int>(child_keys.size()); ++j) {
            readAndWriteRestartData(child_out_dbs,
               child_in_dbs,
               child_keys[j],
               file_mapping,
               num_files_written,
               which_file_mapping,
               total_input_files,
               total_output_files);
         }
      }
   } else if (input_dbs[0]->isInteger(key)) {

      std::vector<int> int_array = input_dbs[0]->getIntegerVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putIntegerVector(key, int_array);
      }

   } else if (input_dbs[0]->isDouble(key)) {

      std::vector<double> double_array = input_dbs[0]->getDoubleVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putDoubleVector(key, double_array);
      }

   } else if (input_dbs[0]->isBool(key)) {

      std::vector<bool> bool_array = input_dbs[0]->getBoolVector(key);

      if (key == "d_write_edges_for_restart") {
         for (int j = 0; j < static_cast<int>(bool_array.size()); ++j) {
            bool_array[j] = false;
         }
      }

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putBoolVector(key, bool_array);
      }

   } else if (input_dbs[0]->isDatabaseBox(key)) {

      std::vector<tbox::DatabaseBox> box_array =
         input_dbs[0]->getDatabaseBoxVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putDatabaseBoxVector(key, box_array);
      }

   } else if (input_dbs[0]->isString(key)) {

      std::vector<std::string> string_array = input_dbs[0]->getStringVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putStringVector(key, string_array);
      }

   } else if (input_dbs[0]->isComplex(key)) {

      std::vector<dcomplex> complex_array = input_dbs[0]->getComplexVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putComplexVector(key, complex_array);
      }

   } else if (input_dbs[0]->isChar(key)) {

      std::vector<char> char_array = input_dbs[0]->getCharVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putCharVector(key, char_array);
      }

   } else if (input_dbs[0]->isFloat(key)) {

      std::vector<float> float_array = input_dbs[0]->getFloatVector(key);

      for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
         output_dbs[i]->putFloatVector(key, float_array);
      }

   } else {

      TBOX_ERROR(
         "The key " << key
                    << " is invalid or not associated with a supported datatype.");

   }
}

/*
 **************************************************************************
 * readAndWritePatchLevelRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWritePatchLevelRestartData(
   std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
   const std::vector<std::shared_ptr<tbox::Database> >& level_in_dbs,
   const std::string& key,
   const int num_files_written,
   const std::vector<int>& input_proc_nums,
   const int total_output_files)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(output_dbs.size() >= 1);
   assert(level_in_dbs.size() >= 1);
   assert(level_in_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //Create an array of level output databases
   std::vector<std::shared_ptr<tbox::Database> > level_out_dbs(
      output_dbs.size());

   for (int i = 0; i < static_cast<int>(output_dbs.size()); ++i) {
      level_out_dbs[i] = output_dbs[i]->putDatabase(key);
   }

   //Read in data that is global to every processor
   bool is_patch_level = level_in_dbs[0]->getBool("d_is_patch_level");
   int version = level_in_dbs[0]->getInteger("HIER_PATCH_LEVEL_VERSION");
   std::vector<tbox::DatabaseBox> box_array(0);
   if (level_in_dbs[0]->keyExists("d_boxes")) {
      level_in_dbs[0]->getDatabaseBoxVector("d_boxes");
   }

   std::vector<int> ratio_to_zero =
      level_in_dbs[0]->getIntegerVector("d_ratio_to_level_zero");
   int number_blocks = level_in_dbs[0]->getInteger("d_number_blocks");
   std::vector<std::vector<tbox::DatabaseBox> > physical_domain(number_blocks);
   for (int nb = 0; nb < number_blocks; ++nb) {
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      physical_domain[nb] =
         level_in_dbs[0]->getDatabaseBoxVector(domain_name);
   }
   int level_number = level_in_dbs[0]->getInteger("d_level_number");
   int next_coarser_level =
      level_in_dbs[0]->getInteger("d_next_coarser_level_number");
   bool in_hierarchy = level_in_dbs[0]->getBool("d_in_hierarchy");
   std::vector<int> ratio_to_coarser =
      level_in_dbs[0]->getIntegerVector("d_ratio_to_coarser_level");

   const int out_size = static_cast<int>(level_out_dbs.size());

   //Write out global data.
   for (int i = 0; i < out_size; ++i) {
      level_out_dbs[i]->putBool("d_is_patch_level", is_patch_level);
      level_out_dbs[i]->putInteger("HIER_PATCH_LEVEL_VERSION", version);
      if (box_array.size() > 0) {
         level_out_dbs[i]->putDatabaseBoxVector("d_boxes", box_array);
      }
      level_out_dbs[i]->putIntegerVector("d_ratio_to_level_zero",
         ratio_to_zero);
      level_out_dbs[i]->putInteger("d_number_blocks", number_blocks);
      for (int nb = 0; nb < number_blocks; ++nb) {
         std::string domain_name = "d_physical_domain_"
            + tbox::Utilities::blockToString(nb);
         level_out_dbs[i]->putDatabaseBoxVector(domain_name,
            physical_domain[nb]);
      }
      level_out_dbs[i]->putInteger("d_level_number", level_number);
      level_out_dbs[i]->putInteger("d_next_coarser_level_number",
         next_coarser_level);
      level_out_dbs[i]->putBool("d_in_hierarchy", in_hierarchy);
      level_out_dbs[i]->putIntegerVector("d_ratio_to_coarser_level",
         ratio_to_coarser);

   }

   std::vector<std::shared_ptr<tbox::Database> >
   mapped_box_level_dbs_in(input_proc_nums.size());

   std::list<int> local_indices_used;
   int max_index_used = 0;

   //Each iteration of this loop processes the patches from one input
   //database.
   for (int i = 0; i < static_cast<int>(input_proc_nums.size()); ++i) {

      std::shared_ptr<tbox::Database> mbl_database =
         level_in_dbs[i]->getDatabase("mapped_box_level");

      mapped_box_level_dbs_in[i] = mbl_database;

      std::shared_ptr<tbox::Database> mapped_boxes_db =
         mbl_database->getDatabase("mapped_boxes");

      std::vector<int> local_indices(0);
      if (mapped_boxes_db->keyExists("local_indices")) {
         local_indices = mapped_boxes_db->getIntegerVector("local_indices");
      }
      std::vector<int> block_ids(0);
      if (mapped_boxes_db->keyExists("block_ids")) {
         block_ids = mapped_boxes_db->getIntegerVector("block_ids");
      }

      //This list will contain all of the patch numbers that came from a
      //single processor.
      int mbs_size = static_cast<int>(local_indices.size());
      std::list<int> input_local_patch_nums;
      std::list<int> input_local_block_ids;
      std::list<int> output_local_patch_nums;
      std::list<int> output_local_block_ids;
      int max_local_indices = 0;
      int min_local_indices = tbox::MathUtilities<int>::getMax();

      for (int j = 0; j < mbs_size; ++j) {
         bool new_patch_num = true;
         for (std::list<int>::iterator p(input_local_patch_nums.begin());
              p != input_local_patch_nums.end(); ++p) {
            if (*p == local_indices[j]) {
               new_patch_num = false;
               break;
            }
         }
         if (new_patch_num) {
            input_local_patch_nums.push_front(local_indices[j]);
            input_local_block_ids.push_front(block_ids[j]);
         }
      }

      if (out_size == 1) {
         bool recompute_local_patch_nums = false;
         if (local_indices_used.size() == 0) {
            for (std::list<int>::iterator ni(input_local_patch_nums.begin());
                 ni != input_local_patch_nums.end(); ++ni) {
               local_indices_used.push_front(*ni);
               max_index_used =
                  tbox::MathUtilities<int>::Max(max_index_used, *ni);
            }
         } else {
            for (std::list<int>::iterator ni(input_local_patch_nums.begin());
                 ni != input_local_patch_nums.end(); ++ni) {
               bool repeat_found = false;
               for (std::list<int>::iterator li(local_indices_used.begin());
                    li != local_indices_used.end(); ++li) {
                  if (*ni == *li) {
                     repeat_found = true;
                     break;
                  }
               }
               if (repeat_found) {
                  recompute_local_patch_nums = true;
                  int new_value = max_index_used + 1;
                  for (int a = 0; a < mbs_size; ++a) {
                     if (local_indices[a] == *ni) {
                        local_indices[a] = new_value;
                     }
                  }
                  local_indices_used.push_front(new_value);
                  max_index_used = new_value;
               } else {
                  local_indices_used.push_front(*ni);
                  max_index_used =
                     tbox::MathUtilities<int>::Max(max_index_used, *ni);
               }
            }
         }

         if (recompute_local_patch_nums) {
            for (int j = 0; j < mbs_size; ++j) {
               bool new_patch_num = true;
               for (std::list<int>::iterator p(output_local_patch_nums.begin());
                    p != output_local_patch_nums.end();
                    ++p) {
                  if (*p == local_indices[j]) {
                     new_patch_num = false;
                     break;
                  }
               }
               if (new_patch_num) {
                  output_local_patch_nums.push_front(local_indices[j]);
                  output_local_block_ids.push_front(block_ids[j]);
               }
            }
         } else {
            output_local_patch_nums = input_local_patch_nums;
            output_local_block_ids = input_local_block_ids;
         }
      } else {
         output_local_patch_nums = input_local_patch_nums;
         output_local_block_ids = input_local_block_ids;
      }

      for (int j = 0; j < mbs_size; ++j) {
         max_local_indices = tbox::MathUtilities<int>::Max(max_local_indices,
               local_indices[j]);
         min_local_indices = tbox::MathUtilities<int>::Min(min_local_indices,
               local_indices[j]);
      }

      int indices_range;
      if (mbs_size) {
         indices_range = max_local_indices - min_local_indices;
      } else {
         indices_range = 0;
         min_local_indices = 0;
      }

      std::vector<int> output_dist_cutoff(out_size);
      for (int j = 0; j < out_size; ++j) {
         output_dist_cutoff[j] = min_local_indices + j
            * (indices_range / out_size);
      }

      //For every patch number, get the patch database from input,
      //create a database for output, and call routine to read and
      //write patch database data.
      std::list<int>::iterator olp(output_local_patch_nums.begin());
      std::list<int>::iterator ilb(input_local_block_ids.begin());
      std::list<int>::iterator olb(output_local_block_ids.begin());
      for (std::list<int>::iterator ilp(input_local_patch_nums.begin());
           ilp != input_local_patch_nums.end(); ) {
         int output_id = 0;
         for (int a = 1; a < out_size; ++a) {
            if (*olp > output_dist_cutoff[a]) {
               output_id = a;
            }
         }

         std::string in_patch_name =
            "level_" + tbox::Utilities::levelToString(level_number)
            + "-patch_" + tbox::Utilities::patchToString(*ilp)
            + "-block_" + tbox::Utilities::blockToString(*ilb);
         std::string out_patch_name =
            "level_" + tbox::Utilities::levelToString(level_number)
            + "-patch_" + tbox::Utilities::patchToString(*olp)
            + "-block_" + tbox::Utilities::blockToString(*olb);

         std::shared_ptr<tbox::Database> patch_in_db =
            level_in_dbs[i]->getDatabase(in_patch_name);

         std::shared_ptr<tbox::Database> patch_out_db =
            level_out_dbs[output_id]->putDatabase(out_patch_name);

         int output_proc = num_files_written + output_id;
         readAndWritePatchRestartData(patch_out_db, patch_in_db, output_proc);

         ++ilp;
         ++olp;
         ++ilb;
         ++olb;
      }
   }

   readAndWriteBoxLevelRestartData(
      level_out_dbs, mapped_box_level_dbs_in,
      "mapped_box_level", num_files_written,
      input_proc_nums, total_output_files);

}

/*
 **************************************************************************
 * readAndWritePatchRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWritePatchRestartData(
   std::shared_ptr<tbox::Database>& patch_out_db,
   const std::shared_ptr<tbox::Database>& patch_in_db,
   const int output_proc)
{
   //Get the keys in the patch input database.
   std::vector<std::string> keys = patch_in_db->getAllKeys();

   //Place the database on arrays of length 1.
   std::vector<std::shared_ptr<tbox::Database> > in_db_array(1);
   std::vector<std::shared_ptr<tbox::Database> > out_db_array(1);

   in_db_array[0] = patch_in_db;
   out_db_array[0] = patch_out_db;

   //Call recursive function to read and write the data associated with each
   //key.
   for (int i = 0; i < static_cast<int>(keys.size()); ++i) {
      if (keys[i] == "d_patch_owner") {
         patch_out_db->putInteger(keys[i], output_proc);
      } else {
         readAndWriteRestartData(out_db_array, in_db_array, keys[i]);
      }
   }
}

/*
 **************************************************************************
 * readAndWriteBoxLevelRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWriteBoxLevelRestartData(
   std::vector<std::shared_ptr<tbox::Database> >& output_dbs,
   const std::vector<std::shared_ptr<tbox::Database> >& level_in_dbs,
   const std::string& key,
   const int num_files_written,
   const std::vector<int>& input_proc_nums,
   const int total_output_files)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(output_dbs.size() >= 1);
   assert(level_in_dbs.size() >= 1);
   assert(level_in_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   const int out_size = static_cast<int>(output_dbs.size());

   //Create an array of level output databases
   std::vector<std::shared_ptr<tbox::Database> > level_out_dbs(out_size);

   for (int i = 0; i < out_size; ++i) {
      level_out_dbs[i] = output_dbs[i]->putDatabase(key);
   }

   bool is_mapped_box_level = level_in_dbs[0]->getBool("d_is_mapped_box_level");
   int version = level_in_dbs[0]->getInteger("HIER_MAPPED_BOX_LEVEL_VERSION");
   int dim = level_in_dbs[0]->getInteger("dim");
   std::vector<std::vector<int> > ratio;
   int b = 0;
   for ( ; ; ++b ) {   
      std::string ratio_name = "d_ratio_" + tbox::Utilities::intToString(b);
      if (level_in_dbs[0]->isInteger(ratio_name)) {
         ratio.push_back(level_in_dbs[0]->getIntegerVector(ratio_name));
      } else {
         break;
      }
   }
   for (int i = 0; i < out_size; ++i) {
      level_out_dbs[i]->putBool("d_is_mapped_box_level", is_mapped_box_level);
      level_out_dbs[i]->putInteger("HIER_MAPPED_BOX_LEVEL_VERSION", version);
      level_out_dbs[i]->putInteger("dim", dim);
      for (int nb = 0; nb < static_cast<int>(ratio.size()); ++nb) {
         std::string ratio_name = "d_ratio_" + tbox::Utilities::intToString(nb);
         level_out_dbs[i]->putIntegerVector(ratio_name, ratio[nb]);
      }
      level_out_dbs[i]->putInteger("d_nproc", total_output_files);
      level_out_dbs[i]->putInteger("d_rank", num_files_written + i);
   }

   std::vector<int>* out_local_indices = new std::vector<int>[out_size];
   std::vector<int>* out_ranks = new std::vector<int>[out_size];
   std::vector<int>* out_periodic_ids =
      new std::vector<int>[out_size];
   std::vector<int>* out_block_ids = new std::vector<int>[out_size];

   std::vector<tbox::DatabaseBox>* out_box_array =
      new std::vector<tbox::DatabaseBox>[out_size];

   int out_vec_size = 0;
   std::vector<int> out_mbs_size(out_size, 0);

   std::list<int> local_indices_used;
   int max_index_used = 0;

   //Each iteration of this loop processes the patches from one input
   //database.
   version = level_in_dbs[0]->getDatabase("mapped_boxes")->getInteger(
         "HIER_BOX_CONTAINER_VERSION");
   for (int i = 0; i < static_cast<int>(input_proc_nums.size()); ++i) {

      std::shared_ptr<tbox::Database> mapped_boxes_in_db =
         level_in_dbs[i]->getDatabase("mapped_boxes");

      int mbs_size = mapped_boxes_in_db->getInteger("mapped_box_set_size");

      std::vector<int> local_indices;
      std::vector<int> ranks;
      std::vector<int> periodic_ids;
      std::vector<int> block_ids;
      std::vector<tbox::DatabaseBox> boxes;

      if (mapped_boxes_in_db->keyExists("local_indices")) {
         local_indices = mapped_boxes_in_db->getIntegerVector("local_indices");
      }
      if (mapped_boxes_in_db->keyExists("ranks")) {
         ranks = mapped_boxes_in_db->getIntegerVector("ranks");
      }
      if (mapped_boxes_in_db->keyExists("periodic_ids")) {
         periodic_ids =
            mapped_boxes_in_db->getIntegerVector("periodic_ids");
      }
      if (mapped_boxes_in_db->keyExists("block_ids")) {
         block_ids = mapped_boxes_in_db->getIntegerVector("block_ids");
      }
      if (mapped_boxes_in_db->keyExists("boxes")) {
         boxes = mapped_boxes_in_db->getDatabaseBoxVector("boxes");
      }

      if (out_size == 1) {
         std::list<int> new_indices;
         for (int k = 0; k < mbs_size; ++k) {
            bool is_new_index = true;
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ++ni) {
               if (local_indices[k] == *ni) {
                  is_new_index = false;
                  break;
               }
            }
            if (is_new_index) {
               new_indices.push_front(local_indices[k]);
            }
         }
         if (local_indices_used.size() == 0) {
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ++ni) {
               local_indices_used.push_front(*ni);
               max_index_used =
                  tbox::MathUtilities<int>::Max(max_index_used, *ni);
            }
         } else {
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ++ni) {
               bool repeat_found = false;
               for (std::list<int>::iterator li(local_indices_used.begin());
                    li != local_indices_used.end();
                    ++li) {
                  if (*ni == *li) {
                     repeat_found = true;
                     break;
                  }
               }
               if (repeat_found) {
                  int new_value = max_index_used + 1;
                  for (int a = 0; a < mbs_size; ++a) {
                     if (local_indices[a] == *ni) {
                        local_indices[a] = new_value;
                     }
                  }
                  local_indices_used.push_front(new_value);
                  max_index_used = new_value;
               } else {
                  local_indices_used.push_front(*ni);
                  max_index_used =
                     tbox::MathUtilities<int>::Max(max_index_used, *ni);
               }
            }
         }
      }

      int max_local_indices = 0;
      int min_local_indices = tbox::MathUtilities<int>::getMax();
      for (int j = 0; j < mbs_size; ++j) {
         max_local_indices = tbox::MathUtilities<int>::Max(max_local_indices,
               local_indices[j]);
         min_local_indices = tbox::MathUtilities<int>::Min(min_local_indices,
               local_indices[j]);
      }

      int indices_range = max_local_indices - min_local_indices;

      std::vector<int> output_dist_cutoff(out_size);
      for (int j = 0; j < out_size; ++j) {
         output_dist_cutoff[j] = min_local_indices + j
            * (indices_range / out_size);
      }

      if (mbs_size > 0) {
         out_vec_size += mbs_size;

         std::vector<int> output_ids(mbs_size);

         for (int k = 0; k < mbs_size; ++k) {
            output_ids[k] = 0;
            for (int a = 1; a < out_size; ++a) {
               if (local_indices[k] > output_dist_cutoff[a]) {
                  output_ids[k] = a;
               }
            }
         }

         for (int j = 0; j < out_size; ++j) {
            out_local_indices[j].reserve(out_vec_size);
            out_ranks[j].reserve(out_vec_size);
            out_periodic_ids[j].reserve(out_vec_size);
            out_block_ids[j].reserve(out_vec_size);

            out_box_array[j].reserve(out_vec_size);
            for (int k = 0; k < mbs_size; ++k) {
               if (output_ids[k] == j) {
                  int output_rank = num_files_written + output_ids[k];
                  out_local_indices[j].push_back(local_indices[k]);
                  out_ranks[j].push_back(output_rank);
                  out_periodic_ids[j].push_back(periodic_ids[k]);
                  out_block_ids[j].push_back(block_ids[k]);
                  out_box_array[j].push_back(boxes[k]);
                  ++out_mbs_size[j];
               }
            }
         }
      }
   }

   for (int j = 0; j < out_size; ++j) {
      std::shared_ptr<tbox::Database> mapped_boxes_out_db =
         level_out_dbs[j]->putDatabase("mapped_boxes");

      mapped_boxes_out_db->putInteger("HIER_BOX_CONTAINER_VERSION", version);
      mapped_boxes_out_db->putInteger("mapped_box_set_size", out_mbs_size[j]);

      if (out_mbs_size[j]) {
         mapped_boxes_out_db->putIntegerVector("local_indices",
            out_local_indices[j]);
         mapped_boxes_out_db->putIntegerVector("periodic_ids",
            out_periodic_ids[j]);
         mapped_boxes_out_db->putIntegerVector("block_ids",
            out_block_ids[j]);
         mapped_boxes_out_db->putIntegerVector("ranks",
            out_ranks[j]);
         mapped_boxes_out_db->putDatabaseBoxVector("boxes",
            out_box_array[j]);

      }
   }

   delete[] out_local_indices;
   delete[] out_ranks;
   delete[] out_periodic_ids;
   delete[] out_block_ids;
   delete[] out_box_array;
/*
 *    for (int j = 0; j < out_size; ++j) {
 *       std::vector<int> out_local_indices;
 *       std::vector<int> out_ranks;
 *       std::vector<int> out_periodic_ids;
 *       out_local_indices.reserve(mbs_size);
 *       out_ranks.reserve(mbs_size);
 *       out_periodic_ids.reserve(mbs_size);
 *
 *       std::vector<tbox::DatabaseBox> out_box_array(mbs_size);
 *
 *       int out_mbs_size = 0;
 *       for (int k = 0; k < mbs_size; ++k) {
 *          int output_id = 0;
 *          for (int a = 1; a < out_size; ++a) {
 *             if (local_indices[k] > output_dist_cutoff[a]) {
 *                output_id = a;
 *             }
 *          }
 *          if (output_id == j) {
 *             int output_rank = num_files_written + output_id;
 *             out_local_indices.push_back(local_indices[k]);
 *             out_ranks.push_back(output_rank);
 *             out_periodic_ids.push_back(periodic_ids[k]);
 *             out_box_array[out_mbs_size++] = boxes[k];
 *          }
 *       }
 *
 *       std::shared_ptr<tbox::Database> mapped_boxes_out_db =
 *          level_out_dbs[j]->putDatabase("mapped_boxes");
 *
 *       mapped_boxes_out_db->
 *          putInteger("HIER_BOX_CONTAINER_VERSION", version);
 *       mapped_boxes_out_db->
 *          putInteger("mapped_box_set_size", out_mbs_size);
 *
 *       if (out_mbs_size) {
 *          mapped_boxes_out_db->putIntegerVector("local_indices",
 *             out_local_indices);
 *          mapped_boxes_out_db->putIntegerVector("periodic_ids"
 *             out_periodic_ids);
 *          mapped_boxes_out_db->putIntegerVector("ranks", out_ranks);
 *          mapped_boxes_out_db->putDatabaseBoxVector("boxes", out_box_array);
 *
 *       }
 *    }
 */
}

#endif
