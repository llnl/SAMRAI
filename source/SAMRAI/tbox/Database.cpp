/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for the SAMRAI database objects
 *
 ************************************************************************/

#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

#include <cstring>

namespace SAMRAI {
namespace tbox {

Database::Database()
{
}

Database::~Database()
{
}

/*
 ************************************************************************
 *
 * Get database entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a database type.
 *
 ************************************************************************
 */

std::shared_ptr<Database>
Database::getDatabaseWithDefault(
   const std::string& key,
   const std::shared_ptr<Database>& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      return getDatabase(key);
   } else {
      return defaultvalue;
   }

}

/*
 * Boolean
 */

/*
 * Create a boolean scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 */

void
Database::putBool(
   const std::string& key,
   const bool& data)
{
   TBOX_ASSERT(!key.empty());

   putBoolArray(key, &data, 1);
}

/*
 * Create a boolean array entry in the database with the specified
 * key
 */

void
Database::putBoolVector(
   const std::string& key,
   const std::vector<bool>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      int nbools = static_cast<int>(data.size());
      bool* bool_array = new bool[nbools];
      for (int i = 0; i < nbools; ++i) {
         bool_array[i] = data[i];
      }
      putBoolArray(key, bool_array, nbools);
      delete[] bool_array;
   } else {
      TBOX_ERROR("Database::putBoolVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool
Database::getBool(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return ret_val;
}

/*
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool
Database::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<bool> local_bool = getBoolVector(key);
      return local_bool.empty() ? defaultvalue : local_bool[0];
   } else {
      return defaultvalue;
   }
}

void
Database::getBoolArray(
   const std::string& key,
   bool* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<bool> tmp = getBoolVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getBoolVector() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }

}

/*
 *************************************************************************
 *
 * Create a box entry in the database with the specified
 * key name.  A box entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   TBOX_ASSERT(!key.empty());

   putDatabaseBoxArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a box vector entry in the database with the specified key name.
 *
 *************************************************************************
 */

void
Database::putDatabaseBoxVector(
   const std::string& key,
   const std::vector<DatabaseBox>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putDatabaseBoxArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putDatabaseBoxVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length vector with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get box scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a box type.
 *
 ************************************************************************
 */

DatabaseBox
Database::getDatabaseBox(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   DatabaseBox ret_val;
   getDatabaseBoxArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get box scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a box type.
 *
 ************************************************************************
 */

DatabaseBox
Database::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<DatabaseBox> local_box = getDatabaseBoxVector(key);
      return local_box.empty() ? defaultvalue : local_box[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<DatabaseBox> tmp = getDatabaseBoxVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDatabaseBoxArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Char
 */

/*
 *************************************************************************
 *
 * Create a char scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putChar(
   const std::string& key,
   const char& data)
{
   TBOX_ASSERT(!key.empty());

   putCharArray(key, &data, 1);

}

/*
 *************************************************************************
 *
 * Create a char vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putCharVector(
   const std::string& key,
   const std::vector<char>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putCharArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putCharVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get char scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a char type.
 *
 ************************************************************************
 */

char
Database::getChar(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   char ret_val;
   getCharArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get char scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a char type.
 *
 ************************************************************************
 */

char
Database::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<char> local_char = getCharVector(key);
      return local_char.empty() ? defaultvalue : local_char[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getCharArray(
   const std::string& key,
   char* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<char> tmp = getCharVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getCharArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Complex
 */

/*
 *************************************************************************
 *
 * Create a complex scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   TBOX_ASSERT(!key.empty());

   putComplexArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a complex vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putComplexVector(
   const std::string& key,
   const std::vector<dcomplex>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putComplexArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putComplexVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get complex scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a complex type.
 *
 ************************************************************************
 */

dcomplex
Database::getComplex(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   dcomplex ret_val;
   getComplexArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get complex scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a complex type.
 *
 ************************************************************************
 */

dcomplex
Database::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   dcomplex retval = defaultvalue;

   if (keyExists(key)) {
      std::vector<dcomplex> local_dcomplex = getComplexVector(key);
      if (!local_dcomplex.empty()) {
         retval = local_dcomplex[0];
      }
   }
   return retval;
}

void
Database::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<dcomplex> tmp = getComplexVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getComplexArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Float
 */

/*
 *************************************************************************
 *
 * Create a float scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putFloat(
   const std::string& key,
   const float& data)
{
   TBOX_ASSERT(!key.empty());

   putFloatArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a float vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putFloatVector(
   const std::string& key,
   const std::vector<float>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putFloatArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putFloatVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }

}

/*
 ************************************************************************
 *
 * Get float scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a float type.
 *
 ************************************************************************
 */

float
Database::getFloat(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   float ret_val;
   getFloatArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get float scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a float type.
 *
 ************************************************************************
 */

float
Database::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<float> local_float = getFloatVector(key);
      return local_float.empty() ? defaultvalue : local_float[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getFloatArray(
   const std::string& key,
   float* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<float> tmp = getFloatVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getFloatArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Double
 */

/*
 *************************************************************************
 *
 * Create a double scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putDouble(
   const std::string& key,
   const double& data)
{
   TBOX_ASSERT(!key.empty());

   putDoubleArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a double vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putDoubleVector(
   const std::string& key,
   const std::vector<double>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putDoubleArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putDoubleVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get double scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a double type.
 *
 ************************************************************************
 */

double
Database::getDouble(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   double ret_val;
   getDoubleArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get double scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a double type.
 *
 ************************************************************************
 */

double
Database::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<double> local_double = getDoubleVector(key);
      return local_double.empty() ? defaultvalue : local_double[0];
   } else {
      return defaultvalue;
   }
}

void
Database::getDoubleArray(
   const std::string& key,
   double* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<double> tmp = getDoubleVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDoubleArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Integer
 */

/*
 *************************************************************************
 *
 * Create a integer scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putInteger(
   const std::string& key,
   const int& data)
{
   TBOX_ASSERT(!key.empty());

   putIntegerArray(key, &data, 1);
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
Database::putIntegerVector(
   const std::string& key,
   const std::vector<int>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putIntegerArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putIntegerVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get integer scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a integer type.
 *
 ************************************************************************
 */

int
Database::getInteger(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   int ret_val;
   getIntegerArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get integer scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a integer type.
 *
 ************************************************************************
 */

int
Database::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<int> local_int = getIntegerVector(key);
      return local_int.empty() ? defaultvalue : local_int[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getIntegerArray(
   const std::string& key,
   int* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<int> tmp = getIntegerVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getIntegerArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * String
 */

/*
 *************************************************************************
 *
 * Create a string scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putString(
   const std::string& key,
   const std::string& data)
{
   TBOX_ASSERT(!key.empty());

   putStringArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a string vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putStringVector(
   const std::string& key,
   const std::vector<std::string>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putStringArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putStringVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get string scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a string type.
 *
 ************************************************************************
 */

std::string
Database::getString(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string ret_val;
   getStringArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get string scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a string type.
 *
 ************************************************************************
 */

std::string
Database::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<std::string> local_string = getStringVector(key);
      return local_string.empty() ? defaultvalue : local_string[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getStringArray(
   const std::string& key,
   std::string* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<std::string> tmp = getStringVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getStringArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }

}

bool
Database::isVector(
   const std::string& key)
{
   return isInteger(key + "_size");
}

void
Database::copyDatabase(const std::shared_ptr<Database>& database)
{
   std::vector<std::string> keys(database->getAllKeys());

   for (std::vector<std::string>::const_iterator k_itr = keys.begin();
        k_itr != keys.end(); ++k_itr) {

      const std::string& key = *k_itr;
      Database::DataType my_type = database->getArrayType(key);

      size_t size =  database->getArraySize(key);
      if (my_type == SAMRAI_DATABASE) {
         std::shared_ptr<Database> child_db = database->getDatabase(key);
         std::shared_ptr<Database> new_db = putDatabase(key);
         new_db->copyDatabase(child_db);
      } else if (my_type == SAMRAI_BOOL) {
         if (size == 1) {
            putBool(key, database->getBool(key));
         } else {
            std::vector<bool> bvec(database->getBoolVector(key));
            putBoolVector(key, bvec);
         }
      } else if (my_type == SAMRAI_CHAR) {
         if (size == 1) {
            putChar(key, database->getChar(key));
         } else {
            std::vector<char> bvec(database->getCharVector(key));
            putCharVector(key, bvec);
         }
      } else if (my_type == SAMRAI_INT) {
         if (size == 1) {
            putInteger(key, database->getInteger(key));
         } else {
            std::vector<int> bvec(database->getIntegerVector(key));
            putIntegerVector(key, bvec);
         }
      } else if (my_type == SAMRAI_COMPLEX) {
         if (size == 1) {
            putComplex(key, database->getComplex(key));
         } else {
            std::vector<dcomplex> bvec(database->getComplexVector(key));
            putComplexVector(key, bvec);
         }
      }  else if (my_type == SAMRAI_DOUBLE) {
         if (size == 1) {
            putDouble(key, database->getDouble(key));
         } else {
            std::vector<double> bvec(database->getDoubleVector(key));
            putDoubleVector(key, bvec);
         }
      } else if (my_type == SAMRAI_FLOAT) {
         if (size == 1) {
            putFloat(key, database->getFloat(key));
         } else {
            std::vector<float> bvec(database->getFloatVector(key));
            putFloatVector(key, bvec);
         }
      } else if (my_type == SAMRAI_STRING) {
         if (size == 1) {
            putString(key, database->getString(key));
         } else {
            std::vector<std::string> bvec(database->getStringVector(key));
            putStringVector(key, bvec);
         }
      } else if (my_type == SAMRAI_BOX) {
         if (size == 1) {
            putDatabaseBox(key, database->getDatabaseBox(key));
         } else {
            std::vector<DatabaseBox> bvec(database->getDatabaseBoxVector(key));
            putDatabaseBoxVector(key, bvec);
         }
      }
   }
}

#ifdef SAMRAI_HAVE_CONDUIT
void
Database::toConduitNode(conduit::Node& node)
{
   node.reset();
   std::vector<std::string> keys(getAllKeys());

   for (std::vector<std::string>::const_iterator k_itr = keys.begin();
        k_itr != keys.end(); ++k_itr) {

      const std::string& key = *k_itr;
      Database::DataType my_type = getArrayType(key);

      size_t size = getArraySize(key);
      if (my_type == SAMRAI_DATABASE) {
         std::shared_ptr<Database> child_db = getDatabase(key);
         child_db->toConduitNode(node[key]);
      } else if (my_type == SAMRAI_BOOL) {
         std::vector<bool> bool_vec(getBoolVector(key));
         node[key].set(conduit::DataType::uint8(size));
         conduit::uint8_array fill_array = node[key].as_uint8_array();
         for (size_t i = 0; i < size; ++i) {
            if (bool_vec[i]) {
               fill_array[i] = 1;
            } else {
               fill_array[i] = 0;
            }
         }
      } else if (my_type == SAMRAI_CHAR) {
         std::vector<char> char_vec(getCharVector(key));
         node[key].set(conduit::DataType::c_char(size));
         std::memcpy(node[key].data_ptr(), &(char_vec[0]), size*sizeof(char));
      } else if (my_type == SAMRAI_INT) {
         std::vector<int> int_vec(getIntegerVector(key));
         node[key].set(int_vec);
      } else if (my_type == SAMRAI_COMPLEX) {
         std::vector<dcomplex> cplx_vec(getComplexVector(key));
         node[key].set(conduit::DataType::c_double(2*size));
         conduit::double_array fill_array = node[key].as_double_array();
         for (size_t i = 0; i < size; ++i) {
            fill_array[i*2] = cplx_vec[i].real();
            fill_array[i*2+1] = cplx_vec[i].imag();
         }
      } else if (my_type == SAMRAI_DOUBLE) {
         node[key].set(getDoubleVector(key));
      } else if (my_type == SAMRAI_FLOAT) {
         node[key].set(getFloatVector(key));
      } else if (my_type == SAMRAI_STRING) {
         std::vector<std::string> str_vec(getStringVector(key));
         if (size == 1) {
            node[key].set(str_vec[0]);
         } else {
            size_t num_chars = 0;
            for (size_t i = 0; i < size; ++i) {
               num_chars += str_vec[i].size();
               ++num_chars;
            }
            node[key].set(conduit::DataType::c_char(num_chars));
            char* char_ptr = static_cast<char*>(node[key].data_ptr());
            for (size_t i = 0; i < size; ++i) {
               std::memcpy(char_ptr, str_vec[i].c_str(), (str_vec[i].size()+1)*sizeof(char));
               char_ptr += (str_vec[i].size()+1);
            }
         }
      } else if (my_type == SAMRAI_BOX) {
         std::vector<DatabaseBox> db_box(getDatabaseBoxVector(key));
         tbox::Dimension dim(db_box[0].getDimVal());
         size_t array_size = size * 2 * dim.getValue();
         node[key].set(conduit::DataType::int64(array_size));
         conduit::int64_array fill_array = node[key].as_int64_array();
         size_t i = 0;
         for (size_t b = 0; b < size; ++b) {
            for (unsigned short d = 0; d < dim.getValue(); ++d) {
               fill_array[i+d] = db_box[b].lower(d);
            }
            i += dim.getValue();
            for (unsigned short d = 0; d < dim.getValue(); ++d) {
               fill_array[i+d] = db_box[b].upper(d);
            }
            i += dim.getValue();
         }
      }
   }
}
#endif


}
}
