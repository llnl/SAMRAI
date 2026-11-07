/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An memory database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabase.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#include <stdlib.h>

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#define MEMORY_DB_ERROR(X) \
   do {                                         \
      pout << "MemoryDatabase: " << X << std::endl << std::flush;       \
      printClassData(pout);                                             \
      pout << "Program abort called..." << std::endl << std::flush;     \
      SAMRAI_MPI::abort();                                              \
   } while (0)

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 * o */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int MemoryDatabase::PRINT_DEFAULT = 1;
const int MemoryDatabase::PRINT_INPUT = 2;
const int MemoryDatabase::PRINT_UNUSED = 4;
const int MemoryDatabase::SSTREAM_BUFFER = 4096;

MemoryDatabase::MemoryDatabase(
   const std::string& name):
   d_database_name(name)
{
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

MemoryDatabase::~MemoryDatabase()
{
}

/*
 *************************************************************************
 *
 * Create memory data file specified by name.
 *
 *************************************************************************
 */

bool
MemoryDatabase::create(
   const std::string& name)
{
   d_database_name = name;
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Open memory data file specified by name
 *
 *************************************************************************
 */

bool
MemoryDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   if (read_write_mode == false) {
      TBOX_ERROR("MemoryDatabase::open: MemoryDatabase only supports\n"
         << "read-write mode.  The read_write_mode flag must be true.");

   }
   d_database_name = name;
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Close the open data file.
 *
 *************************************************************************
 */

bool
MemoryDatabase::close()
{
   d_database_name = "";
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Return whether the key exists in the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::keyExists(
   const std::string& key)
{
   return findKeyData(key) ? true : false;
}

/*
 *************************************************************************
 *
 * Return all of the keys in the database.
 *
 *************************************************************************
 */

std::vector<std::string>
MemoryDatabase::getAllKeys()
{
   std::vector<std::string> keys(d_keyvalues.size());

   std::vector<std::string>::size_type k = 0;
   for (std::map<std::string,KeyData>::iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); ++i) {
      keys[k++] = i->first;
   }

   return keys;
}

/*
 *************************************************************************
 *
 * Get the type of the array entry associated with the specified key
 *
 *************************************************************************
 */
enum Database::DataType
MemoryDatabase::getArrayType(
   const std::string& key) {
   KeyData* keydata = findKeyData(key);

   if (keydata) {
      return keydata->d_type;
   } else {
      return Database::SAMRAI_INVALID;
   }
}

/*
 *************************************************************************
 *
 * Get the size of the array entry associated with the specified key;
 * return 0 if the key does not exist.
 *
 *************************************************************************
 */

size_t
MemoryDatabase::getArraySize(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   if (keydata && keydata->d_type != Database::SAMRAI_DATABASE) {
      return keydata->d_array_size;
   } else {
      return 0;
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage the database values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDatabase(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_DATABASE : false;
}

std::shared_ptr<Database>
MemoryDatabase::putDatabase(
   const std::string& key)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_DATABASE;
   keydata.d_array_size = 1;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_database.reset(new MemoryDatabase(key));
   return keydata.d_database;
}

std::shared_ptr<Database>
MemoryDatabase::getDatabase(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_DATABASE) {
      MEMORY_DB_ERROR("Key=" << key << " is not a database...");
   }
   keydata->d_accessed = true;
   return keydata->d_database;
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isBool(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_BOOL : false;
}

void
MemoryDatabase::putBool(
   const std::string& key,
   const bool& data)
{
   putBoolArray(key, &data, 1);
}

void
MemoryDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_BOOL;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_boolean.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_boolean[i] = data[i];
   }

}

bool
MemoryDatabase::getBool(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_BOOL) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean scalar...");
   }
   keydata->d_accessed = true;
   return keydata->d_boolean[0];
}

bool
MemoryDatabase::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getBool(key);

   putBool(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<bool>
MemoryDatabase::getBoolVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_BOOL) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean...");
   }
   keydata->d_accessed = true;
   return keydata->d_boolean;
}

void
MemoryDatabase::getBoolArray(
   const std::string& key,
   bool* data,
   const size_t nelements)
{
   std::vector<bool> tmp = getBoolVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage box values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDatabaseBox(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_BOX : false;
}

void
MemoryDatabase::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   putDatabaseBoxArray(key, &data, 1);
}

void
MemoryDatabase::putDatabaseBoxVector(
   const std::string& key,
   const std::vector<DatabaseBox>& data)
{
   putDatabaseBoxArray(key, &data[0], data.size());
}

void
MemoryDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_BOX;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_box.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_box[i] = data[i];
   }
}

DatabaseBox
MemoryDatabase::getDatabaseBox(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_BOX) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single box...");
   }
   keydata->d_accessed = true;
   return keydata->d_box[0];
}

DatabaseBox
MemoryDatabase::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getDatabaseBox(key);

   putDatabaseBox(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<DatabaseBox>
MemoryDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_BOX) {
      MEMORY_DB_ERROR("Key=" << key << " is not a box...");
   }
   keydata->d_accessed = true;
   return keydata->d_box;
}

void
MemoryDatabase::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const size_t nelements)
{
   std::vector<DatabaseBox> tmp = getDatabaseBoxVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage character values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isChar(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_CHAR : false;
}

void
MemoryDatabase::putChar(
   const std::string& key,
   const char& data)
{
   putCharArray(key, &data, 1);
}

void
MemoryDatabase::putCharVector(
   const std::string& key,
   const std::vector<char>& data)
{
   putCharArray(key, &data[0], data.size());
}

void
MemoryDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_CHAR;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_char.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_char[i] = data[i];
   }
}

char
MemoryDatabase::getChar(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_CHAR) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single character...");
   }
   keydata->d_accessed = true;
   return keydata->d_char[0];
}

char
MemoryDatabase::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getChar(key);

   putChar(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<char>
MemoryDatabase::getCharVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_CHAR) {
      MEMORY_DB_ERROR("Key=" << key << " is not a character...");
   }
   keydata->d_accessed = true;
   return keydata->d_char;
}

void
MemoryDatabase::getCharArray(
   const std::string& key,
   char* data,
   const size_t nelements)
{
   std::vector<char> tmp = getCharVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage complex values within the database.
 * Note that complex numbers may be promoted from integers, floats,
 * and doubles.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isComplex(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_COMPLEX
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT
                              || keydata->d_type == Database::SAMRAI_DOUBLE);
}

void
MemoryDatabase::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   putComplexArray(key, &data, 1);
}

void
MemoryDatabase::putComplexVector(
   const std::string& key,
   const std::vector<dcomplex>& data)
{
   putComplexArray(key, &data[0], data.size());
}

void
MemoryDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_COMPLEX;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_complex.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_complex[i] = data[i];
   }
}

dcomplex
MemoryDatabase::getComplex(
   const std::string& key)
{
   dcomplex value(0.0, 0.0);
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = dcomplex((double)keydata->d_integer[0], 0.0);
         break;
      case Database::SAMRAI_FLOAT:
         value = dcomplex((double)keydata->d_float[0], 0.0);
         break;
      case Database::SAMRAI_DOUBLE:
         value = dcomplex(keydata->d_double[0], 0.0);
         break;
      case Database::SAMRAI_COMPLEX:
         value = keydata->d_complex[0];
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   keydata->d_accessed = true;
   return value;
}

dcomplex
MemoryDatabase::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getComplex(key);

   putComplex(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<dcomplex>
MemoryDatabase::getComplexVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   std::vector<dcomplex> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array.resize(keydata->d_integer.size());
         for (std::vector<dcomplex>::size_type i = 0;
              i < keydata->d_integer.size(); ++i) {
            array[i] = dcomplex((double)keydata->d_integer[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_FLOAT: {
         array.resize(keydata->d_float.size());
         for (std::vector<dcomplex>::size_type i = 0;
              i < keydata->d_float.size(); ++i) {
            array[i] = dcomplex((double)keydata->d_float[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_DOUBLE: {
         array.resize(keydata->d_double.size());
         for (std::vector<dcomplex>::size_type i = 0;
              i < keydata->d_double.size(); ++i) {
            array[i] = dcomplex(keydata->d_double[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_COMPLEX:
         array = keydata->d_complex;
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a complex...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const size_t nelements)
{
   std::vector<dcomplex> tmp = getComplexVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage double values within the database.
 * Note that doubles may be promoted from integers or floats.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDouble(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_DOUBLE
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT);
}

void
MemoryDatabase::putDouble(
   const std::string& key,
   const double& data)
{
   putDoubleArray(key, &data, 1);
}

void
MemoryDatabase::putDoubleVector(
   const std::string& key,
   const std::vector<double>& data)
{
   putDoubleArray(key, &data[0], data.size());
}

void
MemoryDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_DOUBLE;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_double.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_double[i] = data[i];
   }
}

double
MemoryDatabase::getDouble(
   const std::string& key)
{
   double value = 0.0;
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = (double)keydata->d_integer[0];
         break;
      case Database::SAMRAI_FLOAT:
         value = (double)keydata->d_float[0];
         break;
      case Database::SAMRAI_DOUBLE:
         value = keydata->d_double[0];
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   keydata->d_accessed = true;
   return value;
}

double
MemoryDatabase::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getDouble(key);

   putDouble(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<double>
MemoryDatabase::getDoubleVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   std::vector<double> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array.resize(keydata->d_integer.size());
         for (size_t i = 0; i < keydata->d_integer.size(); ++i) {
            array[i] = (double)keydata->d_integer[i];
         }
         break;
      }
      case Database::SAMRAI_FLOAT: {
         array.resize(keydata->d_float.size());
         for (size_t i = 0; i < keydata->d_float.size(); ++i) {
            array[i] = (double)keydata->d_float[i];
         }
         break;
      }
      case Database::SAMRAI_DOUBLE: {
         array = keydata->d_double;
         break;
      }
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a double...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getDoubleArray(
   const std::string& key,
   double* data,
   const size_t nelements)
{
   std::vector<double> tmp = getDoubleVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage float values within the database.
 * Note that floats may be promoted from integers or truncated from
 * doubles (without a warning).
 *
 *************************************************************************
 */

bool
MemoryDatabase::isFloat(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_DOUBLE
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT);
}

void
MemoryDatabase::putFloat(
   const std::string& key,
   const float& data)
{
   putFloatArray(key, &data, 1);
}

void
MemoryDatabase::putFloatVector(
   const std::string& key,
   const std::vector<float>& data)
{
   putFloatArray(key, &data[0], data.size());
}

void
MemoryDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_FLOAT;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_float.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_float[i] = data[i];
   }
}

float
MemoryDatabase::getFloat(
   const std::string& key)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   float value = 0.0;
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = static_cast<float>(keydata->d_integer[0]);
         break;
      case Database::SAMRAI_FLOAT:
         value = keydata->d_float[0];
         break;
      case Database::SAMRAI_DOUBLE:
         value = static_cast<float>(keydata->d_double[0]);
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   keydata->d_accessed = true;
   return value;
}

float
MemoryDatabase::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getFloat(key);

   putFloat(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<float>
MemoryDatabase::getFloatVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   std::vector<float> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array.resize(keydata->d_integer.size());
         for (size_t i = 0; i < keydata->d_integer.size(); ++i) {
            array[i] = static_cast<float>(keydata->d_integer[i]);
         }
         break;
      }
      case Database::SAMRAI_FLOAT:
         array = keydata->d_float;
         break;
      case Database::SAMRAI_DOUBLE: {
         array.resize(keydata->d_double.size());
         for (size_t i = 0; i < keydata->d_double.size(); ++i) {
            array[i] = static_cast<float>(keydata->d_double[i]);
         }
         break;
      }
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a float...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getFloatArray(
   const std::string& key,
   float* data,
   const size_t nelements)
{
   std::vector<float> tmp = getFloatVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage integer values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isInteger(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : keydata->d_type == Database::SAMRAI_INT;
}

void
MemoryDatabase::putInteger(
   const std::string& key,
   const int& data)
{
   putIntegerArray(key, &data, 1);
}

void
MemoryDatabase::putIntegerVector(
   const std::string& key,
   const std::vector<int>& data)
{
   putIntegerArray(key, &data[0], data.size());
}

void
MemoryDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_INT;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_integer.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_integer[i] = data[i];
   }
}

int
MemoryDatabase::getInteger(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_INT) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer scalar...");
   }
   keydata->d_accessed = true;
   return keydata->d_integer[0];
}

int
MemoryDatabase::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getInteger(key);

   putInteger(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<int>
MemoryDatabase::getIntegerVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_INT) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer...");
   }
   keydata->d_accessed = true;
   return keydata->d_integer;
}

void
MemoryDatabase::getIntegerArray(
   const std::string& key,
   int* data,
   const size_t nelements)
{
   std::vector<int> tmp = getIntegerVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage string values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isString(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : keydata->d_type == Database::SAMRAI_STRING;
}

void
MemoryDatabase::putString(
   const std::string& key,
   const std::string& data)
{
   putStringArray(key, &data, 1);
}

void
MemoryDatabase::putStringVector(
   const std::string& key,
   const std::vector<std::string>& data)
{
   putStringArray(key, &data[0], data.size());
}

void
MemoryDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   KeyData& keydata = d_keyvalues[key];
   keydata.d_type = Database::SAMRAI_STRING;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_string.resize(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      keydata.d_string[i] = data[i];
   }
}

std::string
MemoryDatabase::getString(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_STRING) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single string...");
   }
   keydata->d_accessed = true;
   return keydata->d_string[0];
}

std::string
MemoryDatabase::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getString(key);

   putString(key, defaultvalue);
   d_keyvalues[key].d_from_default = true;
   return defaultvalue;
}

std::vector<std::string>
MemoryDatabase::getStringVector(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_STRING) {
      MEMORY_DB_ERROR("Key=" << key << " is not a string...");
   }
   keydata->d_accessed = true;
   return keydata->d_string;
}

void
MemoryDatabase::getStringArray(
   const std::string& key,
   std::string* data,
   const size_t nelements)
{
   std::vector<std::string> tmp = getStringVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

std::string
MemoryDatabase::getName()
{
   return d_database_name;
}

std::string
MemoryDatabase::getName() const
{
   return d_database_name;
}

/*
 *************************************************************************
 *
 * Search the current database for a matching key.  If found, delete
 * that key and return true.  If the key does not exist, then return
 * false.
 *
 *************************************************************************
 */

bool MemoryDatabase::deleteKeyIfFound(
   const std::string& key)
{
   std::map<std::string,KeyData>::iterator itr = d_keyvalues.find(key);
   if (itr != d_keyvalues.end()) {
      d_keyvalues.erase(itr);
      return true;
   }
   return false;
}

/*
 *************************************************************************
 *
 * Find the key data associated with the specified key and return a
 * pointer to the record.  If no such key data exists, then return 0.
 *
 *************************************************************************
 */

MemoryDatabase::KeyData *
MemoryDatabase::findKeyData(
   const std::string& key)
{
   std::map<std::string,KeyData>::iterator itr = d_keyvalues.find(key);
   if (itr != d_keyvalues.end()) {
      return &(itr->second);
   }
   return 0;
}

/*
 *************************************************************************
 *
 * Find the key data associated with the specified key and return a
 * pointer to the record.  If no such key data exists, then exit with
 * an error message.
 *
 *************************************************************************
 */

MemoryDatabase::KeyData *
MemoryDatabase::findKeyDataOrExit(
   const std::string& key)
{
   std::map<std::string,KeyData>::iterator itr = d_keyvalues.find(key);
   if (itr != d_keyvalues.end()) {
      return &(itr->second);
   }
   MEMORY_DB_ERROR("Key ``" << key << "'' does not exist in the database...");
   return 0;
}

/*
 *************************************************************************
 *
 * Print the entire database to the specified output stream.
 *
 *************************************************************************
 */

void
MemoryDatabase::printClassData(
   std::ostream& os)
{
   printDatabase(os, 0, PRINT_DEFAULT | PRINT_INPUT | PRINT_UNUSED);
}

/*
 *************************************************************************
 *
 * Print database data to the specified output stream.
 *
 *************************************************************************
 */

void
MemoryDatabase::printDatabase(
   std::ostream& os,
   const int indent,
   const int toprint) const
{
   /*
    * Get the maximum key width in the output (excluding databases)
    */

   int width = 0;
   for (std::map<std::string,KeyData>::const_iterator k = d_keyvalues.begin();
        k != d_keyvalues.end(); ++k) {
      const KeyData& keydata = k->second;
      if (((keydata.d_from_default) && (toprint & PRINT_DEFAULT))
          || ((keydata.d_accessed) && (toprint & PRINT_INPUT))
          || (!(keydata.d_accessed) && (toprint & PRINT_UNUSED))) {
         if (keydata.d_type != Database::SAMRAI_DATABASE) {
            const int keywidth = static_cast<int>(k->first.length());
            if (keywidth > width) {
               width = keywidth;
            }
         }
      }
   }

   /*
    * Iterate over all non-database keys in the database and output key values
    */

   indentStream(os, indent);
   os << d_database_name << " {\n";
   for (std::map<std::string,KeyData>::const_iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); ++i) {
      const KeyData& keydata = i->second;
      if (((keydata.d_from_default) && (toprint & PRINT_DEFAULT))
          || ((keydata.d_accessed) && (toprint & PRINT_INPUT))
          || (!(keydata.d_accessed) && (toprint & PRINT_UNUSED))) {

#ifndef LACKS_SSTREAM
         std::ostringstream sstream;
#else
         char sstream_buffer[SSTREAM_BUFFER];
         std::ostrstream sstream(sstream_buffer, SSTREAM_BUFFER);
#endif

         switch (keydata.d_type) {

            case Database::SAMRAI_INVALID: {
               break;
            }

            case Database::SAMRAI_DATABASE: {
               break;
            }

            case Database::SAMRAI_BOOL: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<bool>::size_type n = keydata.d_boolean.size();
               for (std::vector<bool>::size_type j = 0; j < n; ++j) {
                  sstream << (keydata.d_boolean[j] ? "TRUE" : "FALSE");
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_BOX: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<DatabaseBox>::size_type n = keydata.d_box.size();
               for (std::vector<DatabaseBox>::size_type j = 0; j < n; ++j) {
                  const int m = keydata.d_box[j].getDimVal();
                  sstream << "[(";
                  for (int k = 0; k < m; ++k) {
                     sstream << keydata.d_box[j].lower(k);
                     if (k < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << "),(";
                  for (int l = 0; l < m; ++l) {
                     sstream << keydata.d_box[j].upper(l);
                     if (l < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << ")]";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_CHAR: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<char>::size_type n = keydata.d_char.size();
               for (std::vector<char>::size_type j = 0; j < n; ++j) {
                  sstream << "'" << keydata.d_char[j] << "'";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_COMPLEX: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<dcomplex>::size_type n = keydata.d_complex.size();
               for (std::vector<dcomplex>::size_type j = 0; j < n; ++j) {
                  sstream << keydata.d_complex[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_DOUBLE: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<double>::size_type n = keydata.d_double.size();
               for (std::vector<double>::size_type j = 0; j < n; ++j) {
                  sstream << keydata.d_double[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_FLOAT: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<float>::size_type n = keydata.d_float.size();
               for (std::vector<float>::size_type j = 0; j < n; ++j) {
                  sstream << keydata.d_float[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_INT: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<int>::size_type n = keydata.d_integer.size();
               for (std::vector<int>::size_type j = 0; j < n; ++j) {
                  sstream << keydata.d_integer[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_STRING: {
               indentStream(sstream, indent + 3);
               sstream << i->first;
               indentStream(sstream, width - static_cast<int>(i->first.length()));
               sstream << " = ";
               const std::vector<std::string>::size_type n = keydata.d_string.size();
               for (std::vector<std::string>::size_type j = 0; j < n; ++j) {
                  sstream << "\"" << keydata.d_string[j] << "\"";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }
         }

         /*
          * Output whether the key was used or default in column 60
          */

         if (keydata.d_type != Database::SAMRAI_DATABASE) {
#ifndef LACKS_SSTREAM
            const int tab = static_cast<int>(59 - sstream.str().length());
#else
            const int tab = static_cast<int>(59 - sstream.pcount());
#endif
            if (tab > 0) {
               indentStream(sstream, tab);
            }
            if (keydata.d_from_default) {
               sstream << " // from default";
            } else if (keydata.d_accessed) {
               sstream << " // input used";
            } else {
               sstream << " // input not used";
            }

//            sstream << std::endl << ends;
            sstream << std::endl;
            os << sstream.str();
         }
      }
   }

   /*
    * Finally, output all databases in the current key list
    */

   for (std::map<std::string,KeyData>::const_iterator j = d_keyvalues.begin();
        j != d_keyvalues.end(); ++j) {
      if (j->second.d_type == Database::SAMRAI_DATABASE) {
         std::shared_ptr<MemoryDatabase> db(
            SAMRAI_SHARED_PTR_CAST<MemoryDatabase, Database>(j->second.d_database));
         db->printDatabase(os, indent + 3, toprint);
      }
   }

   indentStream(os, indent);
   os << "}\n";
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
