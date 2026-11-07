/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A null database that does nothing for all database methods.
 *
 ************************************************************************/

#ifndef included_tbox_NullDatabase
#define included_tbox_NullDatabase

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

#include <memory>

namespace SAMRAI {
namespace tbox {

/**
 * The NullDatabase provides an implementation of the Database
 * interface with empty methods for the purpose of reducing the
 * the number of guards necessary in methods from other classes that
 * use databases.
 *
 * See the Database class documentation for a description of the
 * generic database interface.
 *
 */

class NullDatabase:public Database
{
public:
   /**
    * The null database constructor creates an empty database with
    * the name "null".
    */
   NullDatabase();

   /**
    * The input database destructor deallocates the data in the database.
    */
   virtual ~NullDatabase();

   /**
    * Create a new database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    */
   virtual bool
   create(
      const std::string& name);

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
      const bool read_write_mode = false);

   /**
    * Close the database.
    *
    * Returns true if successful.
    *
    * If the database is currently open then close it.  This should
    * flush all data to the file (if the database is on disk).
    */
   virtual bool
   close();

   /**
    * Always returns true.
    */
   virtual bool
   keyExists(
      const std::string& key);

   /**
    * Return an empty std::vector<string>.
    */
   virtual std::vector<std::string>
   getAllKeys();

   /**
    * Return INVALID.
    */
   virtual enum DataType
   getArrayType(
      const std::string& key);

   /**
    * Always returns 0.
    */
   virtual size_t
   getArraySize(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isDatabase(
      const std::string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual std::shared_ptr<Database>
   putDatabase(
      const std::string& key);

   /**
    * Returns a pointer to the null database.
    */
   virtual std::shared_ptr<Database>
   getDatabase(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isBool(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const bool * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<bool>.
    */
   virtual std::vector<bool>
   getBoolVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isDatabaseBox(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const DatabaseBox * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<box>.
    */
   virtual std::vector<DatabaseBox>
   getDatabaseBoxVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isChar(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putCharArray(
      const std::string& key,
      const char * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<char>.
    */
   virtual std::vector<char>
   getCharVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isComplex(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const dcomplex * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<dcomplex>.
    */
   virtual std::vector<dcomplex>
   getComplexVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isDouble(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const double * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<double>.
    */
   virtual std::vector<double>
   getDoubleVector(
      const std::string& key);

   /**
    * Always return true.
    */
   virtual bool
   isFloat(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const float * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<float>.
    */
   virtual std::vector<float>
   getFloatVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isInteger(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const int * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<int>.
    */
   virtual std::vector<int>
   getIntegerVector(
      const std::string& key);

   /**
    * Always returns true.
    */
   virtual bool
   isString(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual void
   putStringArray(
      const std::string& key,
      const std::string * const data,
      const size_t nelements);

   /**
    * Returns an empty std::vector<std::string>.
    */
   virtual std::vector<std::string>
   getStringVector(
      const std::string& key);

   /**
    * Does nothing.
    */
   virtual std::string
   getName();

   /**
    * Does nothing.
    */
   virtual void
   printClassData(
      std::ostream& os = pout);

   using Database::putBoolArray;
   using Database::getBoolArray;
   using Database::putDatabaseBoxArray;
   using Database::getDatabaseBoxArray;
   using Database::getDatabaseBoxVector;
   using Database::putCharArray;
   using Database::getCharArray;
   using Database::putComplexArray;
   using Database::getComplexArray;
   using Database::putFloatArray;
   using Database::getFloatArray;
   using Database::putDoubleArray;
   using Database::getDoubleArray;
   using Database::putIntegerArray;
   using Database::getIntegerArray;
   using Database::putStringArray;
   using Database::getStringArray;

private:
   NullDatabase(
      const NullDatabase&);             // not implemented
   NullDatabase&
   operator = (
      const NullDatabase&);             // not implemented

};

}
}

#endif
