/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A database structure that stores HDF5 format data.
 *
 ************************************************************************/

#ifndef included_tbox_HDFDatabase
#define included_tbox_HDFDatabase

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT HDF5
 ************************************************************************
 */
#ifdef HAVE_HDF5

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/PIO.h"

#ifdef RCSID
#undef RCSID
#endif
#include "hdf5.h"

#include <string>
#include <list>
#include <memory>

namespace SAMRAI {
namespace tbox {

/**
 * Class HDFDatabase implements the interface of the Database
 * class to store data in the HDF5 (Hierarchical Data Format) data format.
 *
 * It is assumed that all processors will access the database in the same
 * manner.  Error reporting is done using the SAMRAI error reporting macros.
 *
 * @see Database
 */

class HDFDatabase:public Database
{
public:
   /**
    * The HDF database constructor creates an empty database with the
    * specified name.  By default the database will not be associated
    * with a file until it is mounted, using the mount() member function.
    *
    * The name string is *NOT* the name of the HDF file.
    *
    * @pre !name.empty()
    */
   explicit HDFDatabase(
      const std::string& name);

   /**
    * Constructor used to create sub-databases.
    *
    * @pre !name.empty()
    */
   HDFDatabase(
      const std::string& name,
      hid_t group_ID);

   /**
    * The database destructor closes the HDF5 group or data file.
    */
   virtual ~HDFDatabase();

   /**
    * Return true if the specified key exists in the database
    * and false otherwise.
    *
    * @pre !key.empty()
    */
   virtual bool
   keyExists(
      const std::string& key);

   /**
    * Return a vector of all keys in the current database.  Note that
    * no keys from subdatabases contained in the current database
    * will appear in the array.  To get the keys of any other
    * database, you must call this routine for that database.
    */
   virtual std::vector<std::string>
   getAllKeys();

   /**
    * @brief Return the type of data associated with the key.
    *
    * If the key does not exist, then INVALID is returned
    *
    * @param key Key name in database.
    */
   virtual enum DataType
   getArrayType(
      const std::string& key);

   /**
    * Return the size of the array associated with the key.  If the key
    * does not exist, then zero is returned.
    *
    * @pre !key.empty()
    */
   virtual size_t
   getArraySize(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a database entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    *
    * @pre !key.empty()
    */
   virtual bool
   isDatabase(
      const std::string& key);

   /**
    * Create a new database with the specified key name and return a
    * pointer to it.  A new group with the key name is also created
    * in the data file.
    *
    * @pre !key.empty()
    */
   virtual std::shared_ptr<Database>
   putDatabase(
      const std::string& key);

   /**
    * Get the database with the specified key name.  If the specified
    * key does not represent a database entry in the database, then
    * an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isDatabase(key)
    */
   virtual std::shared_ptr<Database>
   getDatabase(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a boolean entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isBool(
      const std::string& key);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const bool * const data,
      const size_t nelements);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database,
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isBool(key)
    */
   virtual std::vector<bool>
   getBoolVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a box entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isDatabaseBox(
      const std::string& key);

   /**
    * Create a box array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const DatabaseBox * const data,
      const size_t nelements);

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
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a char entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isChar(
      const std::string& key);

   /**
    * Create a character array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putCharArray(
      const std::string& key,
      const char * const data,
      const size_t nelements);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database,
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isChar(key)
    */
   virtual std::vector<char>
   getCharVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a complex entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isComplex(
      const std::string& key);

   /**
    * Create a complex array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const dcomplex * const data,
      const size_t nelements);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isComplex(key)
    */
   virtual std::vector<dcomplex>
   getComplexVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a double entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isDouble(
      const std::string& key);

   /**
    * Create a double array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const double * const data,
      const size_t nelements);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isDouble(key)
    */
   virtual std::vector<double>
   getDoubleVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a float entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isFloat(
      const std::string& key);

   /**
    * Create a float array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const float * const data,
      const size_t nelements);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isFloat(key)
    */
   virtual std::vector<float>
   getFloatVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents an integer entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isInteger(
      const std::string& key);

   /**
    * Create an integer array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const int * const data,
      const size_t nelements);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isInteger(key)
    */
   virtual std::vector<int>
   getIntegerVector(
      const std::string& key);

   /**
    * Return true or false depending on whether the specified key
    * represents a string entry.  If the key does not exist or if
    * the string is empty, then false is returned.
    */
   virtual bool
   isString(
      const std::string& key);

   /**
    * Create a string array entry in the database with the specified
    * key name.
    *
    * @pre !key.empty()
    * @pre data != 0
    */
   virtual void
   putStringArray(
      const std::string& key,
      const std::string * const data,
      const size_t nelements);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database
    * then an error message is printed and the program exits.
    *
    * @pre !key.empty()
    * @pre isString(key)
    */
   virtual std::vector<std::string>
   getStringVector(
      const std::string& key);

   /**
    * Print contents of current database to the specified output stream.
    * If no output stream is specified, then data is written to stream pout.
    * Note that none of the subdatabases contained in the current database
    * will have their contents printed.  To view the contents of any other
    * database, you must call this print routine for that database.
    */
   virtual void
   printClassData(
      std::ostream& os = pout);

   /**
    * Create a new database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    *
    * @pre !name.empty()
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
    *
    * @pre !name.empty()
    */
   virtual bool
   open(
      const std::string& name,
      bool read_write_mode = false);

   /**
    * @brief Attach the Database to an existing HDF file.
    *
    * If an application has an existing HDF file used for restart this
    * method can be used to write SAMRAI restart information to the
    * existing file instead of SAMRAI creating a distinct file.
    *
    * The group_id should be a hid returned from a H5Gcreate call.
    * SAMRAI data will be written within this group.
    *
    * Returns true if attach was successful.
    *
    * @pre group_id > 0
    */
   virtual bool
   attachToFile(
      hid_t group_id);

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
    * @brief Returns the name of this database.
    *
    * The name for the root of the database is the name supplied when creating
    * it.  Names for nested databases are the keyname of the database.
    *
    */
   virtual std::string
   getName();

   /**
    * Return the group_id so VisIt can access an object's HDF database.
    */
   hid_t
   getGroupId()
   {
      return d_group_id;
   }

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
   HDFDatabase();                     // not implemented
   HDFDatabase(
      const HDFDatabase&);            // not implemented
   HDFDatabase&
   operator = (
      const HDFDatabase&);            // not implemented

   /*
    * Static function passed HDF5 iterator routine to look up database keys.
    */
   static herr_t
   iterateKeys(
      hid_t loc_id,
      const char* name,
      void* database);

   /*
    * Static member used to construct list of keys when searching for
    * database keys via the HDF5 iterator routine.
    */
   static void
   addKeyToList(
      const char* name,
      int type,
      void* database);

   /*
    * Private utility routine for inserting array data in the database
    */
   void
   insertArray(
      hid_t parent_id,
      const char* name,
      size_t offset,
      int ndims,
      const hsize_t dim[] /*ndims*/,
      const int* perm,
      hid_t member_id) const;

   /*!
    * @brief Create an HDF compound type for box.
    *
    * When finished, the type should be closed using H5Tclose(hid_t).
    *
    * @param type_spec 'n' mean use native types; 'f' = means use
    *        types for file.
    * @return the compound type id.
    * @internal We currently create a new compound type every
    * time we write a compound.  It would be more efficient to
    * cache the type id for the file.
    *
    * @pre (type_spec == 'n') || (type_spec == 's')
    */
   hid_t
   createCompoundDatabaseBox(
      char type_spec) const;

   /*!
    * @brief Create an HDF compound type for complex.
    *
    * When finished, the type should be closed using H5Tclose(hid_t).
    *
    * @param type_spec 'n' mean use native types; 'f' = means use
    *        types for file.
    * @return the compound type id.
    * @internal We currently create a new compound type every
    * time we write a compound.  It would be more efficient to
    * cache the type id for the file.
    */
   hid_t
   createCompoundComplex(
      char type_spec) const;

   /*
    * Private utility routines for searching keys in database;
    */
   void
   performKeySearch();
   void
   cleanupKeySearch();

   /*!
    * @brief Write attribute for a given dataset.
    *
    * Currently only one attribute is kept for each dataset: its type.
    * The type attribute is used to determine what kind of data the
    * dataset represent.
    *
    * @param type_key Type identifier for the dataset
    * @param dataset_id The HDF dataset id
    */
   void
   writeAttribute(
      int type_key,
      hid_t dataset_id);

   /*!
    * @brief Read attribute for a given dataset.
    *
    * Currently only one attribute is kept for each dataset: its type.
    * The type attribute is returned.
    *
    * @param dataset_id The HDF dataset id
    * @return type attribute
    */
   int
   readAttribute(
      hid_t dataset_id);

   struct hdf_complex {
      double re;
      double im;
   };

   /*
    * The following structure is used to store (key,type) pairs when
    * searching for keys in the database.
    */
   struct KeyData {
      std::string d_key;   // group or dataset name
      int d_type;     // type of entry
   };

   std::string d_top_level_search_group;
   std::string d_group_to_search;
   int d_still_searching;
   int d_found_group;

   /*
    * HDF5 file and group id, boolean flag indicating whether database
    * is associated with a mounted file, and name of this database object.
    */
   /*!
    * @brief Whether database is mounted to a file
    */
   bool d_is_file;
   /*!
    * @brief ID of file attached to database
    *
    * Is either -1 (not mounted to file) or set to the return value from
    * opening a file.
    * Set to -1 on unmounting the file.
    */
   hid_t d_file_id;
   /*!
    * @brief ID of group attached to database
    *
    * A database object is always attached to a group.
    * The group id is set in the constructor when constructing from a group.
    * If the object mounts a file, the group id is the file id.
    */
   hid_t d_group_id;

   /*
    * Name of this database object (passed into constructor)
    */
   const std::string d_database_name;

   /*
    * List of (key,type) pairs assembled when searching for keys.
    */
   std::list<KeyData> d_keydata;

   /*
    *************************************************************************
    *
    * Integer keys for identifying types in HDF5 database.  Negative
    * entries are used to distinguish arrays from scalars when printing
    * key information.
    *
    *************************************************************************
    */
   static const int KEY_DATABASE;
   static const int KEY_BOOL_ARRAY;
   static const int KEY_BOX_ARRAY;
   static const int KEY_CHAR_ARRAY;
   static const int KEY_COMPLEX_ARRAY;
   static const int KEY_DOUBLE_ARRAY;
   static const int KEY_FLOAT_ARRAY;
   static const int KEY_INT_ARRAY;
   static const int KEY_STRING_ARRAY;
   static const int KEY_BOOL_SCALAR;
   static const int KEY_BOX_SCALAR;
   static const int KEY_CHAR_SCALAR;
   static const int KEY_COMPLEX_SCALAR;
   static const int KEY_DOUBLE_SCALAR;
   static const int KEY_FLOAT_SCALAR;
   static const int KEY_INT_SCALAR;
   static const int KEY_STRING_SCALAR;

};

}
}

#endif
#endif
