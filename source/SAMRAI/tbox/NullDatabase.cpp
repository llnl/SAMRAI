/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An input database structure that stores (key,value) pairs
 *
 ************************************************************************/

#include "SAMRAI/tbox/NullDatabase.h"

namespace SAMRAI {
namespace tbox {

NullDatabase::NullDatabase()
{
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

NullDatabase::~NullDatabase()
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
NullDatabase::create(
   const std::string& name)
{
   NULL_USE(name);
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
NullDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   NULL_USE(name);
   NULL_USE(read_write_mode);
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
NullDatabase::close()
{
   return true;
}

/*
 *************************************************************************
 *
 * Always returns true.
 *
 *************************************************************************
 */

bool
NullDatabase::keyExists(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

/*
 *************************************************************************
 *
 * Return an empty std::vector<std::string>.
 *
 *************************************************************************
 */

std::vector<std::string>
NullDatabase::getAllKeys()
{
   std::vector<std::string> keys(0);
   return keys;
}

/*
 *************************************************************************
 *
 * Always returns INVALID.
 *
 *************************************************************************
 */

Database::DataType
NullDatabase::getArrayType(
   const std::string& key)
{
   NULL_USE(key);
   return Database::SAMRAI_INVALID;
}

/*
 *************************************************************************
 *
 * Always returns 0.
 *
 *************************************************************************
 */

size_t
NullDatabase::getArraySize(
   const std::string& key)
{
   NULL_USE(key);
   return 0;
}

/*
 *************************************************************************
 *
 * Member functions that manage the database values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

std::shared_ptr<Database>
NullDatabase::putDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return std::shared_ptr<Database>(this);
}

std::shared_ptr<Database>
NullDatabase::getDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return std::make_shared<NullDatabase>();
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isBool(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<bool>
NullDatabase::getBoolVector(
   const std::string& key)
{
   NULL_USE(key);
   std::vector<bool> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage box values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isDatabaseBox(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<DatabaseBox>
NullDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   NULL_USE(key);

   std::vector<DatabaseBox> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage character values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isChar(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<char>
NullDatabase::getCharVector(
   const std::string& key)
{
   NULL_USE(key);

   std::vector<char> empty(0);
   return empty;
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
NullDatabase::isComplex(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<dcomplex>
NullDatabase::getComplexVector(
   const std::string& key)
{
   NULL_USE(key);

   std::vector<dcomplex> empty(0);
   return empty;
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
NullDatabase::isDouble(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<double>
NullDatabase::getDoubleVector(
   const std::string& key)
{
   NULL_USE(key);
   std::vector<double> empty(0);
   return empty;
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
NullDatabase::isFloat(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<float>
NullDatabase::getFloatVector(
   const std::string& key)
{
   NULL_USE(key);

   std::vector<float> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage integer values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isInteger(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<int>
NullDatabase::getIntegerVector(
   const std::string& key)
{
   NULL_USE(key);

   std::vector<int> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage string values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isString(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

std::vector<std::string>
NullDatabase::getStringVector(
   const std::string& key)
{
   NULL_USE(key);
   std::vector<std::string> empty(0);
   return empty;
}

std::string
NullDatabase::getName(
   void)
{
   return std::string();
}

/*
 *************************************************************************
 *
 * Does nothing.
 *
 *************************************************************************
 */

void
NullDatabase::printClassData(
   std::ostream& os)
{
   NULL_USE(os);
}

}
}
