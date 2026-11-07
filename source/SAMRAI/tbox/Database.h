/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for the SAMRAI database objects
 *
 ************************************************************************/

#ifndef included_tbox_Database
#define included_tbox_Database

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <string>
#include <iostream>
#include <memory>
#ifdef SAMRAI_HAVE_CONDUIT
#include "conduit.hpp"
#endif

#define INPUT_RANGE_ERROR(param_name)                                     \
   TBOX_ERROR(getObjectName() << ": getFromInput() error\n" << param_name \
                              << " out of range.  Check documentation for valid input.\n")

#define INPUT_VALUE_ERROR(param_name)                                     \
   TBOX_ERROR(getObjectName() << ": getFromInput() error\n" << param_name \
                              << " invalid value.  Check documentation for valid input.\n")

namespace SAMRAI {
namespace tbox {

/**
 * @brief Class Database is an abstract base class for the input, restart,
 * and visualization databases.
 *
 * SAMRAI databases store (key,value) pairs in a hierarchical
 * database.  Each value may be another database or a boolean, box,
 * character, double complex, double, float, integer, or string.
 * DatabaseBoxes are stored using the toolbox box structure.
 *
 * Data is entered into the database through methods of the general form
 * putTYPE(key, TYPE) or putTYPEArray(key, TYPE array), where TYPE is the
 * type of value created.  If the specified key already exists in the
 * database, then the existing key is silently deleted.
 *
 * Data is extracted from the database through methods of the general form
 * TYPE = getTYPE(key), where TYPE is the type of value to be returned
 * from the database.  There are two general lookup methods.  In the first,
 * a default value is provided (for scalars only).  If the specified key is
 * not found in the database, then the specified default is returned.  In
 * the second form, no default is provided, and the database exists with
 * an error message and program exits if the key is not found.  The array
 * version of getTYPE() works in a similar fashion.
 */

class Database
{
public:
   /**
    * Enumerated type indicating what type of values is stored in
    * a database entry.  Returned from getType() method.
    *
    * Note: The SAMRAI_ prefix is needed since some poorly written
    *       packages do "#define CHAR" etc.
    */
   enum DataType { SAMRAI_INVALID,
                   SAMRAI_DATABASE,
                   SAMRAI_BOOL,
                   SAMRAI_CHAR,
                   SAMRAI_INT,
                   SAMRAI_COMPLEX,
                   SAMRAI_DOUBLE,
                   SAMRAI_FLOAT,
                   SAMRAI_STRING,
                   SAMRAI_BOX };

   /**
    * The constructor for the database base class does nothing interesting.
    */
   Database();

   /**
    * The virtual destructor for the database base class does nothing
    * interesting.
    */
   virtual ~Database();

   /**
    * Create a new database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    */
   virtual bool
   create(
      const std::string& name) = 0;

   /**
    * Open an existing database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    *
    * @param read_write_mode Open the database in read-write
    * mode instead of read-only mode.
    */
   virtual bool
   open(
      const std::string& name,
      const bool read_write_mode = false) = 0;

   /**
    * Close the database.
    *
    * Returns true if successful.
    *
    * If the database is currently open then close it.  This should
    * flush all data to the file (if the database is on disk).
    */
   virtual bool
   close() = 0;

   /**
    * Return true if the specified key exists in the database and false
    * otherwise.
    *
    * @param key Key name to lookup.
    */
   virtual bool
   keyExists(
      const std::string& key) = 0;

   /**
    * Return all keys in the database.
    */
   virtual std::vector<std::string>
   getAllKeys() = 0;

   /**
    * @brief Return the type of data associated with the key.
    *
    * If the key does not exist, then INVALID is returned
    *
    * @param key Key name in database.
    */
   virtual enum DataType
   getArrayType(
      const std::string& key) = 0;

   /**
    * @brief Return the size of the array associated with the key.
    *
    * If the key does not exist, then zero is returned.  If the key is
    * a database then zero is returned.
    *
    * @param key Key name in database.
    */
   virtual size_t
   getArraySize(
      const std::string& key) = 0;

   /**
    * Return whether the specified key represents a database entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isDatabase(
      const std::string& key) = 0;

   /**
    * Create a new database with the specified key name.  If the key already
    * exists in the database, then the old key record is deleted and the new
    * one is silently created in its place.
    *
    * @param key Key name in database.
    */
   virtual std::shared_ptr<Database>
   putDatabase(
      const std::string& key) = 0;

   /**
    * Get the database with the specified key name.  If the specified
    * key does not exist in the database or it is not a database, then
    * an error message is printed and the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::shared_ptr<Database>
   getDatabase(
      const std::string& key) = 0;

   /**
    * Get a database in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a database,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual std::shared_ptr<Database>
   getDatabaseWithDefault(
      const std::string& key,
      const std::shared_ptr<Database>& defaultvalue);

   /**
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isBool(
      const std::string& key) = 0;

   /**
    * Create a boolean scalar entry in the database with the specified
    * key name.  If thoe key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putBool(
      const std::string& key,
      const bool& data);

   /**
    * Create a boolean vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putBoolVector(
      const std::string& key,
      const std::vector<bool>& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const bool * const data,
      const size_t nelements) = 0;

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual bool
   getBool(
      const std::string& key);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a boolean scalar,
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    */
   virtual bool
   getBoolWithDefault(
      const std::string& key,
      const bool& defaultvalue);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<bool>
   getBoolVector(
      const std::string& key) = 0;

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getBoolArray(
      const std::string& key,
      bool* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a box entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isDatabaseBox(
      const std::string& key) = 0;

   /**
    * Create a box scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Data to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putDatabaseBox(
      const std::string& key,
      const DatabaseBox& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putDatabaseBoxVector(
      const std::string& key,
      const std::vector<DatabaseBox>& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const DatabaseBox * const data,
      const size_t nelements) = 0;

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * box scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual DatabaseBox
   getDatabaseBox(
      const std::string& key);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a box scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual DatabaseBox
   getDatabaseBoxWithDefault(
      const std::string& key,
      const DatabaseBox& defaultvalue);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<DatabaseBox>
   getDatabaseBoxVector(
      const std::string& key) = 0;

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getDatabaseBoxArray(
      const std::string& key,
      DatabaseBox* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isChar(
      const std::string& key) = 0;

   /**
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putChar(
      const std::string& key,
      const char& data);

   /**
    * Create a character vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putCharVector(
      const std::string& key,
      const std::vector<char>& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putCharArray(
      const std::string& key,
      const char * const data,
      const size_t nelements) = 0;

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual char
   getChar(
      const std::string& key);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a character scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual char
   getCharWithDefault(
      const std::string& key,
      const char& defaultvalue);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<char>
   getCharVector(
      const std::string& key) = 0;

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getCharArray(
      const std::string& key,
      char* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a complex entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isComplex(
      const std::string& key) = 0;

   /**
    * Create a complex scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putComplex(
      const std::string& key,
      const dcomplex& data);

   /**
    * Create a complex vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putComplexVector(
      const std::string& key,
      const std::vector<dcomplex>& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const dcomplex * const data,
      const size_t nelements) = 0;

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * complex scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual dcomplex
   getComplex(
      const std::string& key);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a complex scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual dcomplex
   getComplexWithDefault(
      const std::string& key,
      const dcomplex& defaultvalue);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<dcomplex>
   getComplexVector(
      const std::string& key) = 0;

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getComplexArray(
      const std::string& key,
      dcomplex* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isDouble(
      const std::string& key) = 0;

   /**
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putDouble(
      const std::string& key,
      const double& data);

   /**
    * Create a double vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putDoubleVector(
      const std::string& key,
      const std::vector<double>& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const double * const data,
      const size_t nelements) = 0;

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual double
   getDouble(
      const std::string& key);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a double scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual double
   getDoubleWithDefault(
      const std::string& key,
      const double& defaultvalue);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<double>
   getDoubleVector(
      const std::string& key) = 0;

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getDoubleArray(
      const std::string& key,
      double* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isFloat(
      const std::string& key) = 0;

   /**
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putFloat(
      const std::string& key,
      const float& data);

   /**
    * Create a float vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putFloatVector(
      const std::string& key,
      const std::vector<float>& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const float * const data,
      const size_t nelements) = 0;

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual float
   getFloat(
      const std::string& key);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a float scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual float
   getFloatWithDefault(
      const std::string& key,
      const float& defaultvalue);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<float>
   getFloatVector(
      const std::string& key) = 0;

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getFloatArray(
      const std::string& key,
      float* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents an integer entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isInteger(
      const std::string& key) = 0;

   /**
    * Create an integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putInteger(
      const std::string& key,
      const int& data);

   /**
    * Create an integer vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putIntegerVector(
      const std::string& key,
      const std::vector<int>& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const int * const data,
      const size_t nelements) = 0;

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual int
   getInteger(
      const std::string& key);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not an integer scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual int
   getIntegerWithDefault(
      const std::string& key,
      const int& defaultvalue);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<int>
   getIntegerVector(
      const std::string& key) = 0;

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getIntegerArray(
      const std::string& key,
      int* data,
      const size_t nelements);

   /**
    * Return whether the specified key represents a std::string entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isString(
      const std::string& key) = 0;

   /**
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    *
    * @pre !key.empty()
    */
   virtual void
   putString(
      const std::string& key,
      const std::string& data);

   /**
    * Create a string vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Vector with data to put into database.
    *
    * @pre !key.empty()
    * @pre data.size() > 0
    */
   virtual void
   putStringVector(
      const std::string& key,
      const std::vector<std::string>& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    */
   virtual void
   putStringArray(
      const std::string& key,
      const std::string * const data,
      const size_t nelements) = 0;

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    *
    * @pre !key.empty()
    */
   virtual std::string
   getString(
      const std::string& key);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a string scalar,
    * then an error message is printed and the program exits.
    *
    * @param key          Key name in database.
    * @param defaultvalue Default value to return if not found.
    *
    * @pre !key.empty()
    */
   virtual std::string
   getStringWithDefault(
      const std::string& key,
      const std::string& defaultvalue);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string vector, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    */
   virtual std::vector<std::string>
   getStringVector(
      const std::string& key) = 0;

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Array with data to put into database.
    * @param nelements Number of elements to write from array.
    *
    * @pre !key.empty()
    */
   virtual void
   getStringArray(
      const std::string& key,
      std::string* data,
      const size_t nelements);

   /**
    * Get a bool entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * bool scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      bool& scalar)
   {
      scalar = getBool(key);
   }

   /**
    * Get a bool entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a bool array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const bool scalar)
   {
      putBool(key, scalar);
   }

   /**
    * Get a bool entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a bool vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns vector that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<bool>& array)
   {
      array = getBoolVector(key);
   }

   /**
    * Create a bool vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<bool>& array)
   {
      putBoolVector(key, array);
   }

   /**
    * Get a char entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * char scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      char& scalar)
   {
      scalar = getChar(key);
   }

   /**
    * Get a char entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a char array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const char scalar)
   {
      putChar(key, scalar);
   }

   /**
    * Get a char entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a char vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<char>& array)
   {
      array = getCharVector(key);
   }

   /**
    * Create an char vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<char>& array)
   {
      putCharVector(key, array);
   }

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * complex scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      dcomplex& scalar)
   {
      scalar = getComplex(key);
   }

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const dcomplex scalar)
   {
      putComplex(key, scalar);
   }

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns array that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<dcomplex>& array)
   {
      array = getComplexVector(key);
   }

   /**
    * Create a complex vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<dcomplex>& array)
   {
      putComplexVector(key, array);
   }

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      float& scalar)
   {
      scalar = getFloat(key);
   }

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const float scalar)
   {
      putFloat(key, scalar);
   }

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns vector that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<float>& array)
   {
      array = getFloatVector(key);
   }

   /**
    * Create a float vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<float>& array)
   {
      putFloatVector(key, array);
   }

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      double& scalar)
   {
      scalar = getDouble(key);
   }

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const double scalar)
   {
      putDouble(key, scalar);
   }

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns vector that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<double>& array)
   {
      array = getDoubleVector(key);
   }

   /**
    * Create an double vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<double>& array)
   {
      putDoubleVector(key, array);
   }

   /**
    * Get a integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Returns scalar that was read.
    */
   void
   getScalar(
      const std::string& key,
      int& scalar)
   {
      scalar = getInteger(key);
   }

   /**
    * Get a integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void
   putScalar(
      const std::string& key,
      const int scalar)
   {
      putInteger(key, scalar);
   }

   /**
    * Get a integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer vector, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Returns vector that was read.
    */
   void
   getVector(
      const std::string& key,
      std::vector<int>& array)
   {
      array = getIntegerVector(key);
   }

   /**
    * Create an integer vector entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Vector to put into database.
    */
   void
   putVector(
      const std::string& key,
      const std::vector<int>& array)
   {
      putIntegerVector(key, array);
   }

   /**
    * Return whether the specified key represents a vector entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool
   isVector(
      const std::string& key);

   /**
    * Get a std::vector<TYPE> of a collection of objects (as opposed to
    * primitives) from the database with the specified key name. If the
    * specified key does not exist in the database or is not a vector,
    * then an error message is printed and the program exits.
    *
    * TYPE must implement the Database::Serializable interface.
    *
    * @param key     Key name in database.
    * @param vector  Returns the filled in vector.
    */
   template<class TYPE>
   void
   getObjectVector(
      const std::string& key,
      std::vector<TYPE>& vector)
   {
      size_t size = getInteger(key + "_size");
      for (unsigned int i = 0; i < size; ++i) {
         const std::string index_str = Utilities::intToString(i);
         vector[i].getFromRestart(*this, key + "_" + index_str);
      }
   }

   /**
    * Create a vector entry of a collection of objects (as opposed to
    * primitives) in the database with the specified key name.  If the
    * key already exists in the database, then the old key record is
    * deleted and the new one is silently created in its place.
    *
    * TYPE must implement the Database::Serializable interface.
    *
    * @param key    Key name in database.
    * @param vector Vector to put into database.
    */
   template<class TYPE>
   void
   putObjectVector(
      const std::string& key,
      const std::vector<TYPE>& vector)
   {
      unsigned int size = static_cast<int>(vector.size());
      putInteger(key + "_size", size);
      for (unsigned int i = 0; i < size; ++i) {
         const std::string index_str = Utilities::intToString(i);
         vector[i].putToRestart(*this, key + "_" + index_str);
      }
   }

   /**
    * @brief Returns the name of this database.
    *
    * The name for the root of the database is the name supplied when creating
    * it.  Names for nested databases are the keyname of the database.
    *
    */
   virtual std::string
   getName() = 0;

   /*!
    * @brief Full copy of a database
    *
    * @param database  Database to be copied
    */
   virtual void copyDatabase(const std::shared_ptr<Database>& database);

#ifdef SAMRAI_HAVE_CONDUIT
   /*!
    * @brief Write data held in this database to a Conduit Node
    *
    * The hierarchical structure of a SAMRAI database will be replicated in
    * Conduit's hierarchical format.
    *
    * @param node  Output node
    */
   virtual void toConduitNode(conduit::Node& node);
#endif

   /**
    * Print the current database to the specified output stream.  If
    * no output stream is specified, then data is written to stream pout.
    *
    * @param os Output stream.
    */
   virtual void
   printClassData(
      std::ostream& os = pout) = 0;

private:
   // Unimplemented copy constructor.
   Database(
      const Database& other);

   // Unimplemented assignment operator.
   Database&
   operator = (
      const Database& rhs);

};

}
}

#endif
