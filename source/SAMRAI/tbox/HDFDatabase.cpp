/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A database structure that stores HDF5 format data.
 *
 ************************************************************************/

#include "SAMRAI/tbox/HDFDatabase.h"

#ifdef HAVE_HDF5

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstring>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 * Macros starting with H5T_SAMRAI_ are for controlling the data
 * type that is actually written to the file.  As long as
 * these are not "native" types, the file should be portable.
 */

// Type used for writing simple (non-compound) data.
#define H5T_SAMRAI_INT H5T_STD_I32BE
#define H5T_SAMRAI_FLOAT H5T_IEEE_F32BE
#define H5T_SAMRAI_DOUBLE H5T_IEEE_F64BE
#define H5T_SAMRAI_BOOL H5T_STD_I8BE

// Type used for writing the data attribute key.
#define H5T_SAMRAI_ATTR H5T_STD_I8BE

/*
 *************************************************************************
 *
 * Macros to suppress the HDF5 messages sent to standard i/o; handle
 * errors explicity within this code.
 *
 *************************************************************************
 */

/*
 * SGS Note:  Can the new HDF5 stack stuff be a better solution to this?
 */
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
#define BEGIN_SUPPRESS_HDF5_WARNINGS                  \
   {                                                     \
      herr_t (* H5E_saved_efunc)( \
         hid_t, \
         void *) = 0;   \
      void* H5E_saved_edata = 0;                      \
      H5Eget_auto(H5E_DEFAULT, &H5E_saved_efunc, &H5E_saved_edata); \
      H5Eset_auto(H5E_DEFAULT, 0, 0);

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_DEFAULT, H5E_saved_efunc, H5E_saved_edata);  \
   }
#else
#define BEGIN_SUPPRESS_HDF5_WARNINGS                  \
   {                                                     \
      herr_t (* H5E_saved_efunc)( \
         void *) = 0;          \
      void* H5E_saved_edata = 0;                      \
      H5Eget_auto(&H5E_saved_efunc, &H5E_saved_edata);   \
      H5Eset_auto(0, 0);

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_saved_efunc, H5E_saved_edata);      \
   }
#endif

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Integer keys for identifying types in HDF5 database.  Negative
 * entries are used to distinguish arrays from scalars when printing
 * key information.
 *
 *************************************************************************
 */
const int HDFDatabase::KEY_DATABASE = 0;
const int HDFDatabase::KEY_BOOL_ARRAY = 1;
const int HDFDatabase::KEY_BOX_ARRAY = 2;
const int HDFDatabase::KEY_CHAR_ARRAY = 3;
const int HDFDatabase::KEY_COMPLEX_ARRAY = 4;
const int HDFDatabase::KEY_DOUBLE_ARRAY = 5;
const int HDFDatabase::KEY_FLOAT_ARRAY = 6;
const int HDFDatabase::KEY_INT_ARRAY = 7;
const int HDFDatabase::KEY_STRING_ARRAY = 8;
const int HDFDatabase::KEY_BOOL_SCALAR = -1;
const int HDFDatabase::KEY_BOX_SCALAR = -2;
const int HDFDatabase::KEY_CHAR_SCALAR = -3;
const int HDFDatabase::KEY_COMPLEX_SCALAR = -4;
const int HDFDatabase::KEY_DOUBLE_SCALAR = -5;
const int HDFDatabase::KEY_FLOAT_SCALAR = -6;
const int HDFDatabase::KEY_INT_SCALAR = -7;
const int HDFDatabase::KEY_STRING_SCALAR = -8;

/*
 *************************************************************************
 *
 * Static member function to iterate through the hdf5 data file and
 * assemble a list of desired (key, type) pairs.
 *
 *************************************************************************
 */

herr_t
HDFDatabase::iterateKeys(
   hid_t loc_id,
   const char* name,
   void* void_database)
{
   TBOX_ASSERT(name != 0);

   HDFDatabase* database = (HDFDatabase *)(void_database);

   if (database->d_still_searching) {

      H5G_stat_t statbuf;
      int type_key;
      herr_t errf;
      NULL_USE(errf);

      errf = H5Gget_objinfo(loc_id, name, 0, &statbuf);
      TBOX_ASSERT(errf >= 0);

      switch (statbuf.type) {
         case H5G_GROUP: {
            if (database->d_top_level_search_group == "/") {
               addKeyToList(name, KEY_DATABASE, void_database);
            } else if (!strcmp(name, database->d_group_to_search.c_str())) {
               hid_t grp;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
               grp = H5Gopen(loc_id, name, H5P_DEFAULT);
#else
               grp = H5Gopen(loc_id, name);
#endif

               TBOX_ASSERT(grp >= 0);

               database->d_found_group = true;
               database->d_still_searching =
                  H5Giterate(grp, ".", 0,
                     HDFDatabase::iterateKeys, void_database);
               TBOX_ASSERT(database->d_still_searching >= 0);

               database->d_found_group = false;
            } else {
               hid_t grp;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
               grp = H5Gopen(loc_id, name, H5P_DEFAULT);
#else
               grp = H5Gopen(loc_id, name);
#endif

               TBOX_ASSERT(grp >= 0);

               if (database->d_found_group) {
                  addKeyToList(name, KEY_DATABASE, void_database);
               } else {
                  errf = H5Giterate(grp, ".", 0,
                        HDFDatabase::iterateKeys, void_database);

                  TBOX_ASSERT(errf >= 0);

               }
            }
            break;
         }

         case H5G_DATASET: {
            if (database->d_still_searching && database->d_found_group) {
               hid_t this_set;
               BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
               this_set = H5Dopen(loc_id, name, H5P_DEFAULT);
#else
               this_set = H5Dopen(loc_id, name);
#endif
               END_SUPPRESS_HDF5_WARNINGS
               if (this_set > 0) {
                  hid_t attr = H5Aopen_name(this_set, "Type");
                  TBOX_ASSERT(attr >= 0);

                  errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
                  TBOX_ASSERT(errf >= 0);

                  hid_t this_space = H5Dget_space(this_set);
                  TBOX_ASSERT(this_space >= 0);

                  hsize_t nsel = H5Sget_select_npoints(this_space);
                  int array_size = int(nsel);
                  addKeyToList(name,
                     (array_size == 1 ? -type_key : type_key),
                     void_database);
                  errf = H5Sclose(this_space);
                  TBOX_ASSERT(errf >= 0);

                  errf = H5Aclose(attr);
                  TBOX_ASSERT(errf >= 0);

                  errf = H5Dclose(this_set);
                  TBOX_ASSERT(errf >= 0);
               }
            }
            break;
         }

         default: {
            TBOX_ERROR("HDFDatabase key search error....\n"
               << "   Unable to identify key = " << name
               << " as a known group or dataset" << std::endl);
         }
      }

   }
   return 0;
}

/*
 *************************************************************************
 *
 * Static member function to add key to list for database associated
 * with void* argument.
 *
 *************************************************************************
 */

void
HDFDatabase::addKeyToList(
   const char* name,
   int type,
   void* database)
{
   TBOX_ASSERT(name != 0);
   TBOX_ASSERT(database != 0);

   KeyData key_item;
   key_item.d_key = name;
   key_item.d_type = type;

   ((HDFDatabase *)database)->d_keydata.push_back(key_item);
}

/*
 *************************************************************************
 *
 * Public HDF database constructor creates an empty database with the
 * specified name.  It sets the group_ID to a default value of -1.
 * This data is used by member functions to track parent databases.
 *
 *************************************************************************
 */

HDFDatabase::HDFDatabase(
   const std::string& name):
   d_still_searching(0),
   d_found_group(0),
   d_is_file(false),
   d_file_id(-1),
   d_group_id(-1),
   d_database_name(name)
{

   TBOX_ASSERT(!name.empty());

   d_keydata.clear();
}

/*
 *************************************************************************
 *
 * Private HDF database constructor creates an empty database with the
 * specified name.  The group_ID is used privately within
 * the member functions to track parent databases.
 *
 *************************************************************************
 */

HDFDatabase::HDFDatabase(
   const std::string& name,
   hid_t group_ID):
   d_is_file(false),
   d_file_id(-1),
   d_group_id(group_ID),
   d_database_name(name)
{

   TBOX_ASSERT(!name.empty());

   d_keydata.clear();
}

/*
 *************************************************************************
 *
 * The database destructor closes the opened file or group.
 *
 *************************************************************************
 */

HDFDatabase::~HDFDatabase()
{
   herr_t errf;
   NULL_USE(errf);

   if (d_is_file) {
      close();
   }

   if (d_group_id != -1) {
      errf = H5Gclose(d_group_id);
      TBOX_ASSERT(errf >= 0);
   }

}

/*
 *************************************************************************
 *
 * Return true if the key exists within the database; false otherwise.
 *
 *************************************************************************
 */

bool
HDFDatabase::keyExists(
   const std::string& key)
{

   TBOX_ASSERT(!key.empty());

   bool key_exists = false;
   herr_t errf;
   NULL_USE(errf);

   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   this_set = H5Dopen(d_group_id, key.c_str());
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_set > 0) {
      key_exists = true;
      errf = H5Dclose(this_set);

      TBOX_ASSERT(errf >= 0);
   }
   if (!key_exists) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_group = H5Gopen(d_group_id, key.c_str());
#endif

      END_SUPPRESS_HDF5_WARNINGS
      if (this_group > 0) {
         key_exists = true;
         errf = H5Gclose(this_group);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return key_exists;
}

/*
 *************************************************************************
 *
 * Return all keys in the database.
 *
 *************************************************************************
 */

std::vector<std::string>
HDFDatabase::getAllKeys()
{
   performKeySearch();

   std::vector<std::string> tmp_keys(
      static_cast<std::vector<std::string>::size_type>(d_keydata.size()));

   size_t k = 0;
   for (std::list<KeyData>::iterator i = d_keydata.begin();
        i != d_keydata.end(); ++i) {
      tmp_keys[k] = i->d_key;
      ++k;
   }

   cleanupKeySearch();

   return tmp_keys;
}

/*
 *************************************************************************
 *
 * Get the type of the array entry associated with the specified key
 *
 *************************************************************************
 */
enum Database::DataType
HDFDatabase::getArrayType(
   const std::string& key) {

   enum Database::DataType type = Database::SAMRAI_INVALID;

   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {

      if (isDatabase(key)) {
         type = Database::SAMRAI_DATABASE;
      } else {

         hid_t this_set;
         BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
         this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
         this_set = H5Dopen(d_group_id, key.c_str());
#endif
         END_SUPPRESS_HDF5_WARNINGS;
         if (this_set > 0) {
            int type_key = readAttribute(this_set);

            switch (type_key) {
               case KEY_DATABASE:
                  type = Database::SAMRAI_DATABASE;
                  break;
               case KEY_BOOL_ARRAY:
                  type = Database::SAMRAI_BOOL;
                  break;
               case KEY_BOX_ARRAY:
                  type = Database::SAMRAI_BOX;
                  break;
               case KEY_CHAR_ARRAY:
                  type = Database::SAMRAI_CHAR;
                  break;
               case KEY_COMPLEX_ARRAY:
                  type = Database::SAMRAI_COMPLEX;
                  break;
               case KEY_DOUBLE_ARRAY:
                  type = Database::SAMRAI_DOUBLE;
                  break;
               case KEY_FLOAT_ARRAY:
                  type = Database::SAMRAI_FLOAT;
                  break;
               case KEY_INT_ARRAY:
                  type = Database::SAMRAI_INT;
                  break;
               case KEY_STRING_ARRAY:
                  type = Database::SAMRAI_STRING;
                  break;
               case KEY_BOOL_SCALAR:
                  type = Database::SAMRAI_BOOL;
                  break;
               case KEY_BOX_SCALAR:
                  type = Database::SAMRAI_BOX;
                  break;
               case KEY_CHAR_SCALAR:
                  type = Database::SAMRAI_CHAR;
                  break;
               case KEY_COMPLEX_SCALAR:
                  type = Database::SAMRAI_COMPLEX;
                  break;
               case KEY_DOUBLE_SCALAR:
                  type = Database::SAMRAI_DOUBLE;
                  break;
               case KEY_FLOAT_SCALAR:
                  type = Database::SAMRAI_FLOAT;
                  break;
               case KEY_INT_SCALAR:
                  type = Database::SAMRAI_INT;
                  break;
               case KEY_STRING_SCALAR:
                  type = Database::SAMRAI_STRING;
                  break;
            }

            errf = H5Dclose(this_set);
            TBOX_ASSERT(errf >= 0);
         }
      }
   }
   return type;
}

/*
 *************************************************************************
 *
 * Return the size of the array associated with the key.  If the key
 * does not exist, then zero is returned.
 * Array size is set based on the number of elements (points) within
 * the dataspace defined by the named dataset (or key).
 *
 *************************************************************************
 */

size_t
HDFDatabase::getArraySize(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   size_t array_size = 0;

   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   this_set = H5Dopen(d_group_id, key.c_str());
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_set > 0) {
      hid_t this_space = H5Dget_space(this_set);
      TBOX_ASSERT(this_space >= 0);

      hsize_t nsel;
      if (readAttribute(this_set) == KEY_CHAR_ARRAY) {
         hid_t dtype = H5Dget_type(this_set);
         TBOX_ASSERT(dtype >= 0);

         nsel = H5Tget_size(dtype);
      } else {
         nsel = H5Sget_select_npoints(this_space);
      }
      array_size = static_cast<size_t>(nsel);
      errf = H5Sclose(this_space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(this_set);
      TBOX_ASSERT(errf >= 0);

   }

   return array_size;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a database entry.  If the key does not exist, then false
 * is returned.  The key represents a database (or hdf group) if the
 * H5Gopen function on the key is successful.
 *
 *************************************************************************
 */

bool
HDFDatabase::isDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool is_database = false;
   herr_t errf;
   NULL_USE(errf);

   hid_t this_group;
   BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   this_group = H5Gopen(d_group_id, key.c_str());
#endif
   END_SUPPRESS_HDF5_WARNINGS;
   if (this_group > 0) {
      is_database = true;
      errf = H5Gclose(this_group);

      TBOX_ASSERT(errf >= 0);
   }

   return is_database;
}

/*
 *************************************************************************
 *
 * Create a new database with the specified key name.
 *
 *************************************************************************
 */

std::shared_ptr<Database>
HDFDatabase::putDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t this_group = H5Gcreate(d_group_id,
         key.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);
#else
   hid_t this_group = H5Gcreate(d_group_id, key.c_str(), 0);
#endif

   TBOX_ASSERT(this_group >= 0);

   std::shared_ptr<Database> new_database(
      std::make_shared<HDFDatabase>(key, this_group));

   return new_database;
}

/*
 ************************************************************************
 *
 * Get the database with the specified key name.
 *
 ************************************************************************
 */

std::shared_ptr<Database>
HDFDatabase::getDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isDatabase(key)) {
      TBOX_ERROR("HDFDatabase::getDatabase() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a database." << std::endl);
   }

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t this_group = H5Gopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   hid_t this_group = H5Gopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(this_group >= 0);

   std::shared_ptr<Database> database(
      std::make_shared<HDFDatabase>(key, this_group));

   return database;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a boolean entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isBool(
   const std::string& key)
{
   bool is_boolean = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_BOOL_ARRAY) {
            is_boolean = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_boolean;
}

/*
 *************************************************************************
 *
 * Create a boolean array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
HDFDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hsize_t dim[1] = { nelements };
      hid_t space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

      /*
       * We cannot be sure exactly what bool is because it is
       * represented differently on different platforms, and
       * it may have been redefined, i.e., by the Boolean
       * type.  We are unsure what the bool is so we convert it
       * to the native int type (H5T_NATIVE_INT) before giving
       * it to HDF.  When we write a bool, we write it the
       * shortest integer type we can find, the H5T_SAMRAI_BOOL
       * type.
       */
      std::vector<int> data1(nelements);
      for (size_t i = 0; i < nelements; ++i) data1[i] = data[i];

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_BOOL,
            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_BOOL,
            space, H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data1[0]);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_BOOL_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putBoolArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get boolean vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a boolean type.
 *
 ************************************************************************
 */

std::vector<bool>
HDFDatabase::getBoolVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isBool(key)) {
      TBOX_ERROR("HDFDatabase::getBoolVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a bool array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;
   herr_t errf;
   NULL_USE(errf);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   nsel = H5Sget_select_npoints(dspace);

   std::vector<bool> bool_array(
      static_cast<std::vector<bool>::size_type>(nsel));

   if (nsel > 0) {
      /*
       * We cannot be sure exactly what bool is because it is
       * represented differently on different platforms, and
       * it may have been redefined, i.e., by the Boolean
       * type.  So we read bools into native integer memory
       * then convert.
       */
      std::vector<int> data1(static_cast<std::vector<int>::size_type>(nsel));
      int* locPtr = &data1[0];
      errf = H5Dread(dset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            locPtr);
      TBOX_ASSERT(errf >= 0);

      // Convert what was just read in.
      for (unsigned int i = 0; i < nsel; ++i) bool_array[i] = data1[i];
   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return bool_array;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a box entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isDatabaseBox(
   const std::string& key)
{
   bool is_box = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_BOX_ARRAY) {
            is_box = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_box;
}

/*
 *************************************************************************
 *
 * Create a box array entry in the database with the specified key name.
 *
 *************************************************************************
 */

void
HDFDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   if (nelements > 0) {

      herr_t errf;
      NULL_USE(errf);

      // Memory type
      hid_t mtype = createCompoundDatabaseBox('n');
      // Storage type
      hid_t stype = createCompoundDatabaseBox('s');

      hsize_t length = nelements;
      hid_t space = H5Screate_simple(1, &length, 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset =
         H5Dcreate(d_group_id, key.c_str(), stype, space, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset =
         H5Dcreate(d_group_id, key.c_str(), stype, space, H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_BOX_ARRAY, dataset);

      errf = H5Tclose(mtype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(stype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putDatabaseBoxArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * A routine to get a box vector from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a box type.
 *
 ************************************************************************
 */

std::vector<DatabaseBox>
HDFDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isDatabaseBox(key)) {
      TBOX_ERROR("HDFDatabase::getDatabaseBoxVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a box array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;
   herr_t errf;
   NULL_USE(errf);

   // Memory type
   hid_t mtype = createCompoundDatabaseBox('n');

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif

   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   nsel = H5Sget_select_npoints(dspace);

   std::vector<DatabaseBox> boxVector(
      static_cast<std::vector<DatabaseBox>::size_type>(nsel));

   if (nsel > 0) {
      DatabaseBox* locPtr = &boxVector[0];
      errf = H5Dread(dset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);

      /*
       * This seems a little hacky but HDF and POD stuff is making this ugly.
       */
      for (size_t i = 0; i < static_cast<size_t>(nsel); ++i) {
         locPtr[i].setDim(Dimension(static_cast<unsigned short>(locPtr[i].d_data.d_dimension)));
      }
   }

   errf = H5Tclose(mtype);
   TBOX_ASSERT(errf >= 0);

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return boxVector;
}

hid_t
HDFDatabase::createCompoundDatabaseBox(
   char type_spec) const {

   herr_t errf;
   NULL_USE(errf);

   hid_t int_type_spec = H5T_SAMRAI_INT;
   switch (type_spec) {
      case 'n':
         // Use native type specs.
         int_type_spec = H5T_NATIVE_INT;
         break;
      case 's':
         // Use storage type specs.
         int_type_spec = H5T_SAMRAI_INT;
         break;
      default:
         TBOX_ERROR(
         "HDFDatabase::createCompundDatabaseBox() error in database "
         << d_database_name
         << "\n    Unknown type specifier found. " << std::endl);
   }
   hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(DatabaseBox));
   TBOX_ASSERT(type >= 0);

   errf = H5Tinsert(type, "dim", HOFFSET(DatabaseBox_POD, d_dimension),
         int_type_spec);
   TBOX_ASSERT(errf >= 0);

   const hsize_t box_dim = SAMRAI::MAX_DIM_VAL;
   insertArray(type, "lo", HOFFSET(DatabaseBox_POD, d_lo), 1, &box_dim,
      0, int_type_spec);
   insertArray(type, "hi", HOFFSET(DatabaseBox_POD, d_hi), 1, &box_dim,
      0, int_type_spec);
   return type;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a char entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isChar(
   const std::string& key)
{
   bool is_char = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_CHAR_ARRAY) {
            is_char = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_char;
}

/*
 *************************************************************************
 *
 * Create a char array entry in the database with the specified
 * key name. The charentry is defined by the hdf type H5T_C_S1.
 *
 *************************************************************************
 */

void
HDFDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hid_t atype, space, dataset;

      char* local_buf = new char[nelements];

      atype = H5Tcopy(H5T_C_S1);
      TBOX_ASSERT(atype >= 0);

      errf = H5Tset_size(atype, nelements);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
      TBOX_ASSERT(errf >= 0);

      for (size_t i = 0; i < nelements; ++i) {
         local_buf[i] = data[i];
      }

      space = H5Screate(H5S_SCALAR);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      dataset = H5Dcreate(d_group_id, key.c_str(), atype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      dataset = H5Dcreate(d_group_id, key.c_str(), atype, space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, local_buf);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_CHAR_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(atype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putCharArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get char vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a char type.
 *
 ************************************************************************
 */

std::vector<char>
HDFDatabase::getCharVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isChar(key)) {
      TBOX_ERROR("HDFDatabase::getCharVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a char array." << std::endl);
   }

   hid_t dset, dspace, dtype;
   size_t nsel = 0;
   herr_t errf;
   NULL_USE(errf);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   dtype = H5Dget_type(dset);
   TBOX_ASSERT(dtype >= 0);

   nsel = H5Tget_size(dtype);

   std::vector<char> charArray(
      static_cast<std::vector<char>::size_type>(nsel));

   if (nsel > 0) {
      char* locPtr = &charArray[0];
      errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);
   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Tclose(dtype);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return charArray;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a complex entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isComplex(
   const std::string& key)
{
   bool is_complex = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_COMPLEX_ARRAY) {
            is_complex = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_complex;
}

/*
 *************************************************************************
 *
 * Create a complex array entry in the database with the specified
 * key name.  The complex array is a compound type based on the hdf
 * type H5T_NATIVE_DOUBLE (for real and imag parts).
 *
 *************************************************************************
 */

void
HDFDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hid_t space, dataset;

      // Memory type
      hid_t mtype = createCompoundComplex('n');
      // Storage type
      hid_t stype = createCompoundComplex('s');

      hsize_t dim[] = { nelements };
      space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      dataset = H5Dcreate(d_group_id, key.c_str(), stype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      dataset = H5Dcreate(d_group_id, key.c_str(), stype, space,
            H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_COMPLEX_ARRAY, dataset);

      errf = H5Tclose(mtype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(stype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putComplexArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get complex vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a complex type.
 *
 ************************************************************************
 */

std::vector<dcomplex>
HDFDatabase::getComplexVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   if (!isComplex(key)) {
      TBOX_ERROR("HDFDatabase::getComplexVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a complex array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;

   // Memory type
   hid_t mtype = createCompoundComplex('n');

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   nsel = H5Sget_select_npoints(dspace);

   std::vector<dcomplex> complexArray(
      static_cast<std::vector<dcomplex>::size_type>(nsel));

   if (nsel > 0) {
      dcomplex* locPtr = &complexArray[0];
      errf = H5Dread(dset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);
   }

   errf = H5Tclose(mtype);
   TBOX_ASSERT(errf >= 0);

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return complexArray;
}

hid_t
HDFDatabase::createCompoundComplex(
   char type_spec) const {

   herr_t errf;
   NULL_USE(errf);

   hid_t double_type_spec = H5T_SAMRAI_DOUBLE;
   switch (type_spec) {
      case 'n':
         // Use native type specs.
         double_type_spec = H5T_NATIVE_DOUBLE;
         break;
      case 's':
         // Use storage type specs.
         double_type_spec = H5T_SAMRAI_DOUBLE;
         break;
      default:
         TBOX_ERROR("HDFDatabase::createCompundComplex() error in database "
         << d_database_name
         << "\n    Unknown type specifier found. " << std::endl);
   }
   hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(dcomplex));
   TBOX_ASSERT(type >= 0);

   errf = H5Tinsert(type, "real", HOFFSET(hdf_complex, re), double_type_spec);
   TBOX_ASSERT(errf >= 0);

   errf = H5Tinsert(type, "imag", HOFFSET(hdf_complex, im), double_type_spec);
   TBOX_ASSERT(errf >= 0);

   return type;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a double entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isDouble(
   const std::string& key)
{
   bool is_double = false;

   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_DOUBLE_ARRAY) {
            is_double = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_double;
}

/*
 *************************************************************************
 *
 * Create a double array entry in the database with the specified
 * key name.  The array type is based on the hdf type H5T_NATIVE_HDOUBLE.
 *
 *************************************************************************
 */

void
HDFDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hsize_t dim[] = { nelements };
      hid_t space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_DOUBLE,
            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_DOUBLE,
            space, H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_DOUBLE_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putDoubleArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get double vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a double type.
 *
 ************************************************************************
 */

std::vector<double>
HDFDatabase::getDoubleVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   if (!isDouble(key)) {
      TBOX_ERROR("HDFDatabase::getDoubleVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a double array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   nsel = H5Sget_select_npoints(dspace);

   std::vector<double> doubleArray(
      static_cast<std::vector<double>::size_type>(nsel));

   if (nsel > 0) {
      double* locPtr = &doubleArray[0];
      errf = H5Dread(dset, H5T_NATIVE_DOUBLE,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);
   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return doubleArray;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a float entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isFloat(
   const std::string& key)
{
   bool is_float = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_FLOAT_ARRAY) {
            is_float = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_float;
}

/*
 *************************************************************************
 *
 * Create a float array entry in the database with the specified
 * key name.  The array type is based on the hdf type H5T_NATIVE_HFLOAT.
 *
 *************************************************************************
 */

void
HDFDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hsize_t dim[] = { nelements };
      hid_t space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_FLOAT,
            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_FLOAT,
            space, H5P_DEFAULT);
#endif

      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_FLOAT_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putFloatArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get float vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a float type.
 *
 ************************************************************************
 */

std::vector<float>
HDFDatabase::getFloatVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   if (!isFloat(key)) {
      TBOX_ERROR("HDFDatabase::getFloatVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a float array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif

   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);
   nsel = H5Sget_select_npoints(dspace);

   std::vector<float> floatArray(
      static_cast<std::vector<float>::size_type>(nsel));

   if (nsel > 0) {
      float* locPtr = &floatArray[0];
      errf = H5Dread(dset, H5T_NATIVE_FLOAT,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);
   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return floatArray;

}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a integer entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isInteger(
   const std::string& key)
{
   bool is_int = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_INT_ARRAY) {
            is_int = true;
         }

         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_int;
}

/*
 *************************************************************************
 *
 * Create an integer array entry in the database with the specified
 * key name.  The array type is based on the hdf type H5T_NATIVE_HINT.
 *
 *************************************************************************
 */

void
HDFDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      hsize_t dim[] = { nelements };
      hid_t space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_INT,
            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_SAMRAI_INT,
            space, H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, data);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_INT_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

   } else {
      TBOX_ERROR("HDFDatabase::putIntegerArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get integer vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a integer type.
 *
 ************************************************************************
 */

std::vector<int>
HDFDatabase::getIntegerVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   if (!isInteger(key)) {
      TBOX_ERROR("HDFDatabase::getIntegerVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not an integer array." << std::endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   nsel = H5Sget_select_npoints(dspace);

   std::vector<int> intArray(static_cast<std::vector<int>::size_type>(nsel));

   if (nsel > 0) {
      int* locPtr = &intArray[0];
      errf = H5Dread(dset, H5T_NATIVE_INT,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
      TBOX_ASSERT(errf >= 0);

   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return intArray;
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a string entry.  If the key does not exist, then false
 * is returned.
 *
 *************************************************************************
 */

bool
HDFDatabase::isString(
   const std::string& key)
{
   bool is_string = false;
   herr_t errf;
   NULL_USE(errf);

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
      this_set = H5Dopen(d_group_id, key.c_str());
#endif
      END_SUPPRESS_HDF5_WARNINGS;
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_STRING_ARRAY) {
            is_string = true;
         }
         errf = H5Dclose(this_set);
         TBOX_ASSERT(errf >= 0);
      }
   }

   return is_string;
}

/*
 *************************************************************************
 *
 * Create a double array entry in the database with the specified
 * key name.  The array type is based on the hdf type H5T_C_S1.
 *
 *************************************************************************
 */

void
HDFDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   herr_t errf;
   NULL_USE(errf);

   if (nelements > 0) {

      int maxlen = 0;
      int current, data_size;
      size_t i;
      for (i = 0; i < nelements; ++i) {
         current = static_cast<int>(data[i].size());
         if (current > maxlen) maxlen = current;
      }

      char* local_buf = new char[nelements * (maxlen + 1)];
      for (i = 0; i < nelements; ++i) {
         strcpy(&local_buf[i * (maxlen + 1)], data[i].c_str());
         data_size = static_cast<int>(data[i].size());
         if (data_size < maxlen) {
            memset(&local_buf[i * (maxlen + 1)] + data_size + 1, 0,
               maxlen - data_size);
         }
      }

      hid_t atype = H5Tcopy(H5T_C_S1);
      TBOX_ASSERT(atype >= 0);

      errf = H5Tset_size(atype, maxlen + 1);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
      TBOX_ASSERT(errf >= 0);

      hsize_t dim[] = { nelements };
      hid_t space = H5Screate_simple(1, dim, 0);
      TBOX_ASSERT(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(),
            atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(),
            atype, space, H5P_DEFAULT);
#endif
      TBOX_ASSERT(dataset >= 0);

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, local_buf);
      TBOX_ASSERT(errf >= 0);

      // Write attribute so we know what kind of data this is.
      writeAttribute(KEY_STRING_ARRAY, dataset);

      errf = H5Sclose(space);
      TBOX_ASSERT(errf >= 0);

      errf = H5Tclose(atype);
      TBOX_ASSERT(errf >= 0);

      errf = H5Dclose(dataset);
      TBOX_ASSERT(errf >= 0);

      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putStringArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Two routines to get string arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a string type.
 *
 ************************************************************************
 */

std::vector<std::string>
HDFDatabase::getStringVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   herr_t errf;
   NULL_USE(errf);

   if (!isString(key)) {
      TBOX_ERROR("HDFDatabase::getStringVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a string array." << std::endl);
   }

   hsize_t nsel;
   size_t dsize;
   hid_t dset, dspace, dtype;

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
   dset = H5Dopen(d_group_id, key.c_str());
#endif
   TBOX_ASSERT(dset >= 0);

   dspace = H5Dget_space(dset);
   TBOX_ASSERT(dspace >= 0);

   dtype = H5Dget_type(dset);
   TBOX_ASSERT(dtype >= 0);

   dsize = H5Tget_size(dtype);
   nsel = H5Sget_select_npoints(dspace);

   std::vector<std::string> stringArray(
      static_cast<std::vector<std::string>::size_type>(nsel));

   if (nsel > 0) {
      char* local_buf = new char[nsel * dsize];

      errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_buf);
      TBOX_ASSERT(errf >= 0);

      for (std::vector<std::string>::size_type i = 0;
           i < static_cast<std::vector<std::string>::size_type>(nsel); ++i) {
         std::string* locPtr = &stringArray[i];
         *locPtr = &local_buf[i * dsize];
      }
      delete[] local_buf;
   }

   errf = H5Sclose(dspace);
   TBOX_ASSERT(errf >= 0);

   errf = H5Tclose(dtype);
   TBOX_ASSERT(errf >= 0);

   errf = H5Dclose(dset);
   TBOX_ASSERT(errf >= 0);

   return stringArray;
}

void
HDFDatabase::writeAttribute(
   int type_key,
   hid_t dataset_id)
{
   herr_t errf;
   NULL_USE(errf);

   hid_t attr_id = H5Screate(H5S_SCALAR);
   TBOX_ASSERT(attr_id >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_SAMRAI_ATTR,
         attr_id, H5P_DEFAULT, H5P_DEFAULT);
#else
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_SAMRAI_ATTR,
         attr_id, H5P_DEFAULT);
#endif
   TBOX_ASSERT(attr >= 0);

   errf = H5Awrite(attr, H5T_NATIVE_INT, &type_key);
   TBOX_ASSERT(errf >= 0);

   errf = H5Aclose(attr);
   TBOX_ASSERT(errf >= 0);

   errf = H5Sclose(attr_id);
   TBOX_ASSERT(errf >= 0);
}

int
HDFDatabase::readAttribute(
   hid_t dataset_id)
{
   herr_t errf;
   NULL_USE(errf);

   hid_t attr = H5Aopen_name(dataset_id, "Type");
   TBOX_ASSERT(attr >= 0);

   int type_key;
   errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
   TBOX_ASSERT(errf >= 0);

   errf = H5Aclose(attr);
   TBOX_ASSERT(errf >= 0);

   return type_key;
}

/*
 *************************************************************************
 *
 * Print contents of current database to the specified output stream.
 * Note that contents of subdatabases will not be printed.  This must
 * be done by iterating through all the subdatabases individually.
 *
 *************************************************************************
 */

void
HDFDatabase::printClassData(
   std::ostream& os)
{

   performKeySearch();

   if (d_keydata.size() == 0) {
      os << "Database named `" << d_database_name
         << "' has zero keys..." << std::endl;
   } else {
      os << "Printing contents of database named `"
         << d_database_name << "'..." << std::endl;
   }

   for (std::list<KeyData>::iterator i = d_keydata.begin();
        i != d_keydata.end(); ++i) {
      int t = i->d_type;
      switch (MathUtilities<int>::Abs(t)) {
         case KEY_DATABASE: {
            os << "   Data entry `" << i->d_key << "' is"
               << " a database" << std::endl;
            break;
         }
         case KEY_BOOL_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a boolean ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_BOX_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a box ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_CHAR_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a char ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_COMPLEX_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a complex ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_DOUBLE_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a double ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_FLOAT_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a float ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_INT_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " an integer ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         case KEY_STRING_ARRAY: {
            os << "   Data entry `" << i->d_key << "' is" << " a string ";
            os << ((t < 0) ? "scalar" : "array") << std::endl;
            break;
         }
         default: {
            TBOX_ERROR("HDFDatabase::printClassData error....\n"
               << "   Unable to identify key = " << i->d_key
               << " as a known group or dataset" << std::endl);
         }
      }
   }

   cleanupKeySearch();

}

/*
 *************************************************************************
 *
 * Create HDF data file specified by name.
 *
 *************************************************************************
 */

bool
HDFDatabase::create(
   const std::string& name)
{
   TBOX_ASSERT(!name.empty());

   bool status = false;

   hid_t file_id = 0;

   file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC,
         H5P_DEFAULT, H5P_DEFAULT);
   if (file_id < 0) {
      TBOX_ERROR("Unable to open HDF5 file " << name << "\n");
      status = false;
   } else {
      status = true;
      d_is_file = true;
      d_group_id = file_id;
      d_file_id = file_id;
   }

   return status;
}

/*
 *************************************************************************
 *
 * Open HDF data file specified by name
 *
 *************************************************************************
 */

bool
HDFDatabase::open(
   const std::string& name,
   const bool read_write_mode) {
   TBOX_ASSERT(!name.empty());

   bool status = false;

   hid_t file_id = 0;

   file_id = H5Fopen(name.c_str(),
         read_write_mode ? H5F_ACC_RDWR : H5F_ACC_RDONLY,
         H5P_DEFAULT);
   if (file_id < 0) {
      TBOX_ERROR("Unable to open HDF5 file " << name << "\n");
      status = false;
   } else {
      status = true;
      d_is_file = true;
      d_group_id = file_id;
      d_file_id = file_id;
   }

   return status;

}

/*
 *************************************************************************
 *
 * Close the open HDF data file specified by d_file_id.
 *
 *************************************************************************
 */

bool
HDFDatabase::close()
{
   herr_t errf = 0;
   NULL_USE(errf);

   if (d_is_file) {
      errf = H5Fclose(d_file_id);
      TBOX_ASSERT(errf >= 0);

      if (d_group_id == d_file_id) {
         d_group_id = -1;
      }
      d_file_id = -1;
      d_is_file = false;
   }

   if (errf >= 0) {
      return true;
   } else {
      return false;
   }
}

/*
 *************************************************************************
 *
 * Private helper function for writing arrays in HDF5.  This function
 * was deprecated in HDF5 1.4.  We replicate it here since it makes
 * arrays easier to use in this database class.
 *
 *************************************************************************
 */

void
HDFDatabase::insertArray(
   hid_t parent_id,
   const char* name,
   size_t offset,
   int ndims,
   const hsize_t dim[] /*ndims*/,
   const int* perm,
   hid_t member_id) const
{
   herr_t errf;
   NULL_USE(errf);
   NULL_USE(perm);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 2))

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
   hid_t array = H5Tarray_create(member_id, ndims, dim);
#else
   /*
    * Note that perm is NOT used by HDF, see HDF documentation.
    */
   hid_t array = H5Tarray_create(member_id, ndims, dim, perm);
#endif
   TBOX_ASSERT(array >= 0);

   errf = H5Tinsert(parent_id, name, offset, array);
   TBOX_ASSERT(errf >= 0);

   errf = H5Tclose(array);
   TBOX_ASSERT(errf >= 0);

#else
   size_t newdim[H5S_MAX_RANK];
   for (int i = 0; i < ndims; ++i) {
      newdim[i] = dim[i];
   }

   errf = H5Tinsert_array(parent_id,
         name,
         offset,
         ndims,
         newdim,
         perm,
         member_id);
   TBOX_ASSERT(errf >= 0);
#endif
}

/*
 *************************************************************************
 *
 * Private helper function for searching database keys.
 *
 *************************************************************************
 */

void
HDFDatabase::performKeySearch()
{
   herr_t errf;
   NULL_USE(errf);

   if (d_is_file) {
      d_group_to_search = "/";
      d_top_level_search_group = "/";
      d_found_group = 1;
   } else {
      d_group_to_search = d_database_name;
      d_top_level_search_group = std::string();
      d_found_group = 0;
   }

   d_still_searching = 1;

   errf = H5Giterate(d_group_id, "/", 0,
         HDFDatabase::iterateKeys, (void *)this);
   TBOX_ASSERT(errf >= 0);
}

void
HDFDatabase::cleanupKeySearch()
{
   d_top_level_search_group = std::string();
   d_group_to_search = std::string();
   d_still_searching = 0;
   d_found_group = 0;

   d_keydata.clear();
}

/*
 *************************************************************************
 * Attach to an already created HDF file.
 *************************************************************************
 */
bool
HDFDatabase::attachToFile(
   hid_t group_id)
{
   bool status = false;

   if (group_id > 0) {
      status = true;
      d_is_file = false;
      d_file_id = -1;
      d_group_id = group_id;
   } else {
      TBOX_ERROR("HDFDatabase: Invalid fileid supplied to attachToFile"
         << std::endl);
      status = false;
   }

   return status;
}

std::string
HDFDatabase::getName()
{
   return d_database_name;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
