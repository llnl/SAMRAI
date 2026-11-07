/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A database structure that stores Silo format data.
 *
 ************************************************************************/

#include "SAMRAI/tbox/SiloDatabase.h"

#ifdef HAVE_SILO

#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <vector>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 * Implementation notes.
 *
 * Silo allows only alphanumeric and "_" in directory names.  SAMRAI
 * key names are mangled to remove other characters from keynames when
 * put into the Silo database.
 *
 * Compound types in Silo have several listings in the directory
 * structure.  Rather than unscramble the different entries the
 * compound array types are stored in a subdirectory with that
 * key name and the actual compound type is a variable in
 * that subdirectory.
 *
 */

namespace SAMRAI {
namespace tbox {

const std::string SiloDatabase::DATABASE_BOX_NAME = "/__box__";
const std::string SiloDatabase::STRING_ARRAY_NAME = "/__string_array__";
const std::string SiloDatabase::COMPLEX_ARRAY_NAME = "/__complex_array__";

const std::string SiloDatabase::mangleID = "__U";

bool
SiloDatabase::IsValid(
   int i)
{
   return isalnum(i) || i == '_' || i == '/';
}

std::string
SiloDatabase::nameMangle(
   std::string name) {
   std::stringstream mangled_name;

   for (std::string::size_type i = 0; i < name.size(); ++i) {
      if (IsValid(name[i])) {
         mangled_name << name[i];
      } else {
         mangled_name << mangleID << std::hex << static_cast<int>(name[i]);
      }
   }

   return mangled_name.str();
}

std::string
SiloDatabase::nameDemangle(
   std::string name) {

   std::stringstream unmangled_name;
   std::string::size_type start = 0;
   std::string::size_type found_index = name.find(mangleID, start);

   unmangled_name << name.substr(start, (found_index - start));

   while (found_index != std::string::npos) {
      int character;

      std::stringstream hex_value(name.substr(
                                     found_index + mangleID.length(), 2));
      hex_value >> std::hex >> character;

      unmangled_name << (char)character;

      start = found_index + mangleID.length() + 2;
      found_index = name.find(mangleID, start);

      unmangled_name << name.substr(start, found_index - start);
   }

   return unmangled_name.str();
}

/*
 * Public Silo database constructor creates an empty database with the
 * specified name.  It sets the group_ID to a default value of -1.
 * This data is used by member functions to track parent databases.
 */
SiloDatabase::SiloDatabase(
   const std::string& name):
   d_is_file(false),
   d_file(0),
   d_directory("/"),
   d_database_name(name)
{
   TBOX_ASSERT(!name.empty());
}

/*
 * Private Silo database constructor creates an empty database with the
 * specified name.
 */
SiloDatabase::SiloDatabase(
   const std::string& name,
   DBfile* file,
   const std::string& directory,
   const bool create_in_file):
   d_is_file(false),
   d_file(file),
   d_directory(directory),
   d_database_name(name)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(!directory.empty());
   TBOX_ASSERT(file != 0);

   int err;

   if (create_in_file) {
      err = DBMkDir(d_file, nameMangle(d_directory).c_str());
      if (err < 0) {
         TBOX_ERROR("SiloDatabase MkDir failed " << d_directory << std::endl);
      }
   }
}

/*
 * The database destructor closes the opened file.
 */

SiloDatabase::~SiloDatabase()
{
   if (d_is_file) {
      close();
   }
}

/*
 *************************************************************************
 *
 * Create Silo data file specified by name.
 *
 *************************************************************************
 */

bool
SiloDatabase::create(
   const std::string& name)
{
   TBOX_ASSERT(!name.empty());

   bool status = false;

   if (d_file) {
      close();
   }

   d_file = DBCreate(name.c_str(), DB_CLOBBER, DB_LOCAL, 0, DB_PDB);

   if (d_file == 0) {

      TBOX_ERROR("Unable to open Silo file " << name << "\n");
      status = false;

   } else {

      d_is_file = true;
      d_directory = "";
      status = true;
   }

   return status;
}

/*
 *************************************************************************
 *
 * Open Silo data file specified by name
 *
 *************************************************************************
 */

bool
SiloDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   TBOX_ASSERT(!name.empty());

   bool status = false;

   if (d_file) {
      close();
   }

   d_file = DBOpen(name.c_str(),
         DB_UNKNOWN,
         read_write_mode ? DB_APPEND : DB_READ);

   if (d_file == 0) {

      TBOX_ERROR("Unable to open Silo file " << name << "\n");
      status = false;

   } else {

      d_is_file = true;
      d_directory = "/";
      status = true;
   }

   return status;
}

/*
 *************************************************************************
 *
 * Close the open Silo data file specified by d_file_id.
 *
 *************************************************************************
 */

bool
SiloDatabase::close()
{
   int err = 0;
   if (d_is_file) {
      err = DBClose(d_file);
#ifdef ASSERT_SILO_RETURN_VALUES
      TBOX_ASSERT(err >= 0);
#endif
      d_file = 0;
      d_is_file = false;
   }

   if (err >= 0) {
      return true;
   } else {
      return false;
   }
}

/*
 *************************************************************************
 * Attach to an already created Silo file.
 *************************************************************************
 */
bool
SiloDatabase::attachToFile(
   DBfile* file,
   const std::string& directory)
{
   bool status = false;

   if (file != 0) {
      status = true;
      d_is_file = false;
      d_file = file;
      d_directory = directory;

      std::string path = nameMangle(d_directory);
      if (DBInqVarType(d_file, path.c_str()) != DB_DIR) {
         int err = DBMkdir(d_file, path.c_str());
         if (err < 0) {
            TBOX_ERROR(
               "SiloDatabase: MkDir failed " << d_directory << std::endl);
         }
      }

   } else {
      TBOX_ERROR("SiloDatabase: Invalid file supplied to attachToFile"
         << std::endl);
      status = false;
   }

   return status;
}

/*
 *************************************************************************
 *
 * Return true if the key exists within the database; false otherwise.
 *
 *************************************************************************
 */

bool
SiloDatabase::keyExists(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   if (DBInqVarExists(d_file, path.c_str())) {
      return true;
   } else if (DBInqVarType(d_file, path.c_str()) == DB_DIR) {
      return true;
   } else {
      return false;
   }
}

/*
 *************************************************************************
 *
 * Return all keys in the database.
 *
 *************************************************************************
 */

std::vector<std::string>
SiloDatabase::getAllKeys()
{
   TBOX_ASSERT(hasDirectory());

   std::string path = nameMangle(d_directory);

   std::vector<std::string> tmp_keys;

   DBObjectType var_type = DBInqVarType(d_file, path.c_str());
   if (var_type == DB_DIR) {
      char original_path[256];

      DBGetDir(d_file, original_path);

      DBSetDir(d_file, path.c_str());

      DBtoc* toc = DBGetToc(d_file);

      tmp_keys.resize(toc->nvar + toc->ndir);

      for (int i = 0; i < toc->nvar; ++i) {
         tmp_keys[i] = toc->var_names[i];
      }

      for (int i = 0; i < toc->ndir; ++i) {
         tmp_keys[i + toc->nvar] = toc->dir_names[i];
      }

      DBSetDir(d_file, original_path);
   } else {
      TBOX_ERROR("Not a database " << d_directory << std::endl);
   }

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
SiloDatabase::getArrayType(
   const std::string& key) {
   TBOX_ASSERT(!key.empty());

   enum Database::DataType type = Database::SAMRAI_INVALID;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType var_type = DBInqVarType(d_file, path.c_str());

   if (var_type == DB_VARIABLE) {
      int obj_type = DBGetVarType(d_file, path.c_str());

      switch (obj_type) {
         case DB_INT:
            type = Database::SAMRAI_INT;
            break;
         case DB_SHORT:
            type = Database::SAMRAI_BOOL;
            break;
         case DB_FLOAT:
            type = Database::SAMRAI_FLOAT;
            break;
         case DB_DOUBLE:
            type = Database::SAMRAI_DOUBLE;
            break;
         case DB_CHAR:
            type = Database::SAMRAI_CHAR;
            break;
      }
   } else if (var_type == DB_DIR) {

      // Note that some of the types are stored in subdirectories
      // so check to see if this dir is one of the compound array types.
      std::string sub_path = path + COMPLEX_ARRAY_NAME;
      if (DBInqVarExists(d_file, sub_path.c_str())) {
         DBcompoundarray* ca = DBGetCompoundarray(d_file, sub_path.c_str());
         type = Database::SAMRAI_COMPLEX;
         DBFreeCompoundarray(ca);
      } else {
         sub_path = path + STRING_ARRAY_NAME;
         if (DBInqVarExists(d_file, sub_path.c_str())) {
            DBcompoundarray* ca = DBGetCompoundarray(d_file, sub_path.c_str());
            type = Database::SAMRAI_STRING;
            DBFreeCompoundarray(ca);
         } else {
            sub_path = path + DATABASE_BOX_NAME;
            if (DBInqVarExists(d_file, sub_path.c_str())) {
               DBcompoundarray* ca = DBGetCompoundarray(d_file,
                     sub_path.c_str());
               type = Database::SAMRAI_BOX;
               DBFreeCompoundarray(ca);
            } else {
               type = Database::SAMRAI_DATABASE;
            }
         }
      }
   } else {
      // Unrecognized type return INVALID
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
 *
 *************************************************************************
 */

size_t
SiloDatabase::getArraySize(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   int array_size = 0;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType var_type = DBInqVarType(d_file, path.c_str());

   if (var_type == DB_VARIABLE) {
      array_size = getSiloSimpleTypeLength(key);
   } else if (var_type == DB_DIR) {

      // Note that some of the types are stored in subdirectories
      // so check to see if this dir is one of the compound array types.

      std::string sub_path = path + COMPLEX_ARRAY_NAME;
      if (DBInqVarExists(d_file, sub_path.c_str())) {
         DBcompoundarray* ca = DBGetCompoundarray(d_file, sub_path.c_str());
         array_size = ca->elemlengths[0];
         DBFreeCompoundarray(ca);
      } else {
         sub_path = path + STRING_ARRAY_NAME;
         if (DBInqVarExists(d_file, sub_path.c_str())) {
            DBcompoundarray* ca = DBGetCompoundarray(d_file, sub_path.c_str());
            array_size = ca->nelems;
            DBFreeCompoundarray(ca);
         } else {
            sub_path = path + DATABASE_BOX_NAME;
            if (DBInqVarExists(d_file, sub_path.c_str())) {
               DBcompoundarray* ca = DBGetCompoundarray(d_file,
                     sub_path.c_str());
               array_size = ca->elemlengths[0];
               DBFreeCompoundarray(ca);
            } else {
               // Directory or some other unrecognized structure return 0
            }
         }
      }
   } else {
      // Unrecognized type return 0
   }

   return static_cast<size_t>(array_size);
}

/*
 *************************************************************************
 *
 * Return true or false depending on whether the specified key
 * represents a database entry.  If the key does not exist, then false
 * is returned.
 *************************************************************************
 */

bool
SiloDatabase::isDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   return isSiloType(key, DB_DIR);
}

/*
 *************************************************************************
 *
 * Create a new database with the specified key name.
 *
 *************************************************************************
 */

std::shared_ptr<Database>
SiloDatabase::putDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string path = d_directory + "/" + key;
   std::shared_ptr<Database> new_database(
      std::make_shared<SiloDatabase>(
         key,
         d_file,
         path,
         true));

   TBOX_ASSERT(new_database);

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
SiloDatabase::getDatabase(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isDatabase(key)) {
      TBOX_ERROR("SiloDatabase::getDatabase() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a database." << std::endl);
   }

   std::shared_ptr<Database> new_database(
      std::make_shared<SiloDatabase>(
         key,
         d_file,
         d_directory + "/" + key,
         false));

   return new_database;
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
SiloDatabase::isBool(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());
   return isSiloSimpleType(key, DB_SHORT);
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
SiloDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   short temp_array[nelements];

   for (size_t i = 0; i < nelements; ++i) {
      temp_array[i] = data[i];
   }

   putSiloSimpleType(key, temp_array, nelements, DB_SHORT);
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
SiloDatabase::getBoolVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isBool(key)) {
      TBOX_ERROR("SiloDatabase::getBoolVector() error in database "
         << d_database_name << std::endl
         << "    Key = " << key << " is not a bool array." << std::endl);
   }

   std::vector<bool> boolArray(getSiloSimpleTypeLength(key));

   short temp_array[getSiloSimpleTypeLength(key)];

   getSiloSimpleType(key, temp_array);

   for (int i = 0; i < getSiloSimpleTypeLength(key); ++i) {
      boolArray[i] = temp_array[i];
   }

   return boolArray;
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
SiloDatabase::isDatabaseBox(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool is_type = false;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType obj_type = DBInqVarType(d_file, path.c_str());
   if (obj_type == DB_DIR) {

      path += DATABASE_BOX_NAME;

      DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

      if (ca != 0) {
         if (ca->datatype == DB_INT) {
            is_type = true;
         }
      }

      DBFreeCompoundarray(ca);
   }

   return is_type;
}

/*
 *************************************************************************
 *
 * Create a box array entry in the database with the specified key name.
 *
 *************************************************************************
 */
void
SiloDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   const char* elemnames[3];
   int elemlengths[3];

   elemnames[0] = "dimension";
   elemnames[1] = "hi";
   elemnames[2] = "lo";

   int size = 0;
   for (size_t i = 0; i < nelements; ++i) {
      size += data[i].d_data.d_dimension;
   }

   elemlengths[0] = static_cast<int>(nelements);
   elemlengths[1] = static_cast<int>(size);
   elemlengths[2] = static_cast<int>(size);

   std::vector<int> values(nelements + size * 2);

   size_t offset = nelements;
   for (size_t i = 0; i < nelements; ++i) {
      values[i] = data[i].d_data.d_dimension;
      for (int d = 0; d < data[i].d_data.d_dimension; ++d) {
         values[offset] = data[i].d_data.d_lo[d];
         values[offset + size] = data[i].d_data.d_hi[d];
         ++offset;
      }
   }

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   int err = DBMkdir(d_file, path.c_str());
   if (err < 0) {
      TBOX_ERROR("SiloDatabase: MkDir failed " << d_directory << std::endl);
   }

   path = d_directory + "/" + key + DATABASE_BOX_NAME;
   path = nameMangle(path);

   err = DBPutCompoundarray(d_file, path.c_str(),
         const_cast<char **>(elemnames), elemlengths,
         3, &values[0], static_cast<int>(values.size()),
         DB_INT, 0);
   if (err < 0) {
      TBOX_ERROR(
         "SiloDatabase: DBPutCompoundarray failed " << d_directory
                                                    << std::endl);
   }

}

/*
 ************************************************************************
 *
 * A routine to get box vectors from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a box type.
 *
 ************************************************************************
 */

std::vector<DatabaseBox>
SiloDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string path = d_directory + "/" + key + DATABASE_BOX_NAME;
   path = nameMangle(path);

   DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

   std::vector<DatabaseBox> boxVector(ca->elemlengths[0]);

   int* values = static_cast<int *>(ca->values);
   int offset = ca->elemlengths[0];
   for (int i = 0; i < (ca->elemlengths[0]); ++i) {
      TBOX_ASSERT(values[i] > 0 && values[i] <= SAMRAI::MAX_DIM_VAL);
      boxVector[i].d_data.d_dimension = static_cast<Dimension::dir_t>(values[i]);
      /*
       * This preserves old behavior where boxes can be different dims but is
       * likely not supported anywhere else in the library.
       */
      boxVector[i].setDim(Dimension((unsigned short)values[i]));
      for (int d = 0; d < boxVector[i].d_data.d_dimension; ++d) {
         boxVector[i].d_data.d_lo[d] = values[offset];
         boxVector[i].d_data.d_hi[d] = values[offset + ca->elemlengths[1]];
         ++offset;
      }
   }

   DBFreeCompoundarray(ca);

   return boxVector;
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
SiloDatabase::isChar(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());
   return isSiloSimpleType(key, DB_CHAR);
}

/*
 *************************************************************************
 *
 * Create a char array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   putSiloSimpleType(key, data, nelements, DB_CHAR);
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
SiloDatabase::getCharVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isChar(key)) {
      TBOX_ERROR("SiloDatabase::getCharVector() error in database "
         << d_database_name << std::endl
         << "    Key = " << key << " is not a char array." << std::endl);
   }

   std::vector<char> charArray(getSiloSimpleTypeLength(key));

   getSiloSimpleType(key, &charArray[0]);

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
SiloDatabase::isComplex(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool is_type = false;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType obj_type = DBInqVarType(d_file, path.c_str());
   if (obj_type == DB_DIR) {

      path += COMPLEX_ARRAY_NAME;

      DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

      if (ca != 0) {
         if (ca->datatype == DB_DOUBLE) {
            is_type = true;
         }
      }

      DBFreeCompoundarray(ca);
   }

   return is_type;
}

/*
 *************************************************************************
 *
 * Create a complex array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   const char* elemnames[2];
   std::vector<double> values(nelements * 2);
   int elemlengths[2];

   elemnames[0] = "real";
   elemnames[1] = "imag";

   elemlengths[0] = static_cast<int>(nelements);
   elemlengths[1] = static_cast<int>(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      values[i] = data[i].real();
      values[i + nelements] = data[i].imag();
   }

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   int err = DBMkdir(d_file, path.c_str());
   if (err < 0) {
      TBOX_ERROR("SiloDatabase MkDir failed " << d_directory << std::endl);
   }

   path = d_directory + "/" + key + COMPLEX_ARRAY_NAME;
   path = nameMangle(path);

   err = DBPutCompoundarray(d_file, path.c_str(),
         const_cast<char **>(elemnames), elemlengths, 2,
         &values[0], static_cast<int>(values.size()),
         DB_DOUBLE, 0);
   if (err < 0) {
      TBOX_ERROR(
         "SiloDatabase DBPutCompoundarray failed " << d_directory
                                                   << std::endl);
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
SiloDatabase::getComplexVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isComplex(key)) {
      TBOX_ERROR("SiloDatabase::getComplexVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a complex array." << std::endl);
   }

   std::string path = d_directory + "/" + key + COMPLEX_ARRAY_NAME;
   path = nameMangle(path);

   DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

   std::vector<dcomplex> complexArray(ca->elemlengths[0]);

   for (int i = 0; i < ca->elemlengths[0]; ++i) {
      complexArray[i] = dcomplex(static_cast<double *>(ca->values)[i],
            static_cast<double *>(ca->values)[i + ca->elemlengths[0]]);
   }

   DBFreeCompoundarray(ca);

   return complexArray;
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
SiloDatabase::isDouble(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());
   return isSiloSimpleType(key, DB_DOUBLE);
}

/*
 *************************************************************************
 *
 * Create a double array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   putSiloSimpleType(key, data, nelements, DB_DOUBLE);
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
SiloDatabase::getDoubleVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isDouble(key)) {
      TBOX_ERROR("SiloDatabase::getDoubleVector() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a double array." << std::endl);
   }

   std::vector<double> doubleArray(getSiloSimpleTypeLength(key));

   getSiloSimpleType(key, &doubleArray[0]);

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
SiloDatabase::isFloat(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());
   return isSiloSimpleType(key, DB_FLOAT);
}

/*
 *************************************************************************
 *
 * Create a float array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   putSiloSimpleType(key, data, nelements, DB_FLOAT);
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
SiloDatabase::getFloatVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isFloat(key)) {
      TBOX_ERROR("SiloDatabase::getFloatVector() error in database "
         << d_database_name << std::endl
         << "    Key = " << key << " is not a float array." << std::endl);
   }

   std::vector<float> floatArray(getSiloSimpleTypeLength(key));

   getSiloSimpleType(key, &floatArray[0]);

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
SiloDatabase::isInteger(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());
   return isSiloSimpleType(key, DB_INT);
}

/*
 *************************************************************************
 *
 * Create an integer array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   putSiloSimpleType(key, data, nelements, DB_INT);
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
SiloDatabase::getIntegerVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isInteger(key)) {
      TBOX_ERROR("SiloDatabase::getIntegerVector() error in database "
         << d_database_name << std::endl
         << "    Key = " << key << " is not a integer array." << std::endl);
   }

   std::vector<int> intArray(getSiloSimpleTypeLength(key));

   getSiloSimpleType(key, &intArray[0]);

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
SiloDatabase::isString(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool is_type = false;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType obj_type = DBInqVarType(d_file, path.c_str());
   if (obj_type == DB_DIR) {

      path += STRING_ARRAY_NAME;

      DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

      if (ca != 0) {
         if (ca->datatype == DB_CHAR) {
            is_type = true;
         }
      }

      DBFreeCompoundarray(ca);
   }

   return is_type;
}

/*
 *************************************************************************
 *
 * Create a double array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
SiloDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   std::vector<std::string> strings(nelements);
   const char* elemnames[nelements];
   std::string values;
   int elemlengths[nelements];

   for (size_t i = 0; i < nelements; ++i) {
      strings[i] = Utilities::sizetToString(i);
      elemnames[i] = strings[i].c_str();
      elemlengths[i] = static_cast<int>(data[i].size());
      values.append(data[i]);
   }

   values.append("\0");

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   int err = DBMkdir(d_file, path.c_str());
   if (err < 0) {
      TBOX_ERROR("SiloDatabase DBMkdir failed " << d_directory << std::endl);
   }

   path = d_directory + "/" + key + STRING_ARRAY_NAME;
   path = nameMangle(path);

   DBPutCompoundarray(d_file, path.c_str(),
      const_cast<char **>(elemnames), elemlengths,
      static_cast<int>(nelements), const_cast<char *>(values.c_str()),
      static_cast<int>(values.size() + 1),
      DB_CHAR, 0);
   if (err < 0) {
      TBOX_ERROR(
         "SiloDatabase DBPutCompoundarray failed " << d_directory
                                                   << std::endl);
   }

}

/*
 ************************************************************************
 *
 * Two routines to get string vectors and arrays from the database with the
 * specified key name. In any case, an error message is printed and
 * the program exits if the specified key does not exist in the
 * database or is not associated with a string type.
 *
 ************************************************************************
 */

std::vector<std::string>
SiloDatabase::getStringVector(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   if (!isString(key)) {
      TBOX_ERROR("SiloDatabase::getStringVector() error in database "
         << d_database_name << std::endl
         << "    Key = " << key << " is not a string array." << std::endl);
   }

   std::string path = d_directory + "/" + key + STRING_ARRAY_NAME;
   path = nameMangle(path);

   DBcompoundarray* ca = DBGetCompoundarray(d_file, path.c_str());

   std::vector<std::string> stringArray(ca->nelems);

   std::string values = static_cast<char *>(ca->values);

   std::string::size_type start = 0;
   for (int i = 0; i < ca->nelems; ++i) {
      stringArray[i] = values.substr(start, ca->elemlengths[i]);
      start = start + ca->elemlengths[i];
   }

   DBFreeCompoundarray(ca);

   return stringArray;
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
SiloDatabase::printClassData(
   std::ostream& os)
{

   std::vector<std::string> keys = getAllKeys();

   if (keys.size() == 0) {
      os << "Database named `" << d_database_name
         << "' has zero keys..." << std::endl;
   } else {
      os << "Printing contents of database named `"
         << d_database_name << "'..." << std::endl;
   }

   for (int i = 0; i < static_cast<int>(keys.size()); ++i) {
      switch (getArrayType(keys[i])) {
         case Database::SAMRAI_INVALID: {
            os << "   Data entry `" << keys[i] << "' is"
               << " invalid" << std::endl;
            break;
         }
         case Database::SAMRAI_DATABASE: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a database" << std::endl;
            break;
         }
         case Database::SAMRAI_BOOL: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a boolean array"
               << std::endl;
            break;
         }
         case Database::SAMRAI_CHAR: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a char array" << std::endl;
            break;
         }
         case Database::SAMRAI_INT: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a integer array"
               << std::endl;
            break;
         }
         case Database::SAMRAI_COMPLEX: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a complex array"
               << std::endl;
            break;
         }
         case Database::SAMRAI_DOUBLE: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a complex array"
               << std::endl;
            break;
         }
         case Database::SAMRAI_FLOAT: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a float array" << std::endl;
            break;
         }
         case Database::SAMRAI_STRING: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a string array"
               << std::endl;
            break;
         }
         case Database::SAMRAI_BOX: {
            os << "   Data entry `" << keys[i] << "' is"
               << " a datbase box array"
               << std::endl;
            break;
         }

      }
   }

}

bool
SiloDatabase::isSiloType(
   const std::string& key,
   DBObjectType type)
{
   TBOX_ASSERT(!key.empty());

   bool is_type = false;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   DBObjectType obj_type = DBInqVarType(d_file, path.c_str());
   if (obj_type == type) {
      is_type = true;
   }

   return is_type;
}

bool
SiloDatabase::isSiloSimpleType(
   const std::string& key,
   const int simple_type)
{
   TBOX_ASSERT(!key.empty());

   bool is_type = false;

   // First check to see if it is a variable, if so then check to
   // see if it is the correct type of variable
   if (isSiloType(key, DB_VARIABLE)) {
      std::string path = d_directory + "/" + key;
      path = nameMangle(path);
      int obj_type = DBGetVarType(d_file, path.c_str());
      if (obj_type == simple_type) {
         is_type = true;
      }
   }

   return is_type;
}

bool
SiloDatabase::putSiloSimpleType(
   const std::string& key,
   const void* data,
   const size_t nelements,
   const int simple_type)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   int err;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   int dims[1];
   dims[0] = static_cast<int>(nelements);

   err = DBWrite(d_file,
         path.c_str(), const_cast<void *>(data), dims, 1, simple_type);
   if (err < 0) {
      TBOX_ERROR("SiloDatabase DBWrite failed " << key << std::endl);
   }

   return err < 0;
}

bool
SiloDatabase::getSiloSimpleType(
   const std::string& key,
   void* data)
{
   TBOX_ASSERT(!key.empty());
   TBOX_ASSERT(data != 0);

   int err;

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   err = DBReadVar(d_file, path.c_str(), data);
   if (err < 0) {
      TBOX_ERROR("SiloDatabase DBRead failed " << key << std::endl);
   }

   return err < 0;
}

int
SiloDatabase::getSiloSimpleTypeLength(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string path = d_directory + "/" + key;
   path = nameMangle(path);

   return DBGetVarLength(d_file, path.c_str());
}

std::string
SiloDatabase::getName(
   void)
{
   return d_database_name;
}

}
}

#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
