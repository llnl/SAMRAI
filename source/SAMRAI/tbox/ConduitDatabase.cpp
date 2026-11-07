/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An memory database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#ifdef SAMRAI_HAVE_CONDUIT

#include "SAMRAI/tbox/ConduitDatabase.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#include <cstring>

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#define TBOX_CONDUIT_DB_ERROR(X) \
   do {                                         \
      pout << "ConduitDatabase: " << X << std::endl << std::flush;       \
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

const int ConduitDatabase::PRINT_DEFAULT = 1;
const int ConduitDatabase::PRINT_INPUT = 2;
const int ConduitDatabase::PRINT_UNUSED = 4;
const int ConduitDatabase::SSTREAM_BUFFER = 4096;

ConduitDatabase::ConduitDatabase(
   const std::string& name):
   d_database_name(name)
{
   d_node = new conduit::Node();
   d_node->reset();
   d_node->set_dtype(conduit::DataType::object());
}

ConduitDatabase::ConduitDatabase(
   const std::string& name,
   conduit::Node* node):
   d_database_name(name),
   d_node(node)
{
   TBOX_ASSERT(d_node);
   TBOX_ASSERT(d_node->dtype().is_object());
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

ConduitDatabase::~ConduitDatabase()
{
   d_child_dbs.clear();
   if (!(d_node->parent())) {
      delete d_node;
   }
}

/*
 *************************************************************************
 *
 * create not implemented
 *
 *************************************************************************
 */

bool
ConduitDatabase::create(
   const std::string& name)
{
   NULL_USE(name);
   TBOX_ERROR("ConduitDatabase::create: ConduitDatabase::create() not\n"
      << "implemented.");
   return false;
}

/*
 *************************************************************************
 *
 * open not implemented
 *
 *************************************************************************
 */

bool
ConduitDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   NULL_USE(name);
   NULL_USE(read_write_mode);
   TBOX_ERROR("ConduitDatabase::open: ConduitDatabase::open() not\n"
      << "implemented.");
   return false;
}

/*
 *************************************************************************
 *
 * close not implemented
 *
 *************************************************************************
 */

bool
ConduitDatabase::close()
{
   TBOX_ERROR("ConduitDatabase::close: ConduitDatabase::close() not\n"
      << "implemented.");
   return false;
}

/*
 *************************************************************************
 *
 * Return whether the key exists in the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::keyExists(
   const std::string& key)
{
   return (d_node->has_child(key));
}

/*
 *************************************************************************
 *
 * Return all of the keys in the database.
 *
 *************************************************************************
 */

std::vector<std::string>
ConduitDatabase::getAllKeys()
{
   return d_node->child_names();
}

/*
 *************************************************************************
 *
 * Get the type of the array entry associated with the specified key
 *
 *************************************************************************
 */
enum Database::DataType
ConduitDatabase::getArrayType(
   const std::string& key)
{
   if (isBool(key)) {
      return SAMRAI_BOOL;
   } else if (isChar(key)) {
      return SAMRAI_CHAR;
   } else if (isComplex(key)) {
      return SAMRAI_COMPLEX;
   } else if (isDatabase(key)) {
      return SAMRAI_DATABASE;
   } else if (isDatabaseBox(key)) {
      return SAMRAI_BOX;
   } else if (isDouble(key)) {
      return SAMRAI_DOUBLE;
   } else if (isFloat(key)) {
      return SAMRAI_FLOAT;
   } else if (isInteger(key)) {
      return SAMRAI_INT;
   } else if (isString(key)) {
      return SAMRAI_STRING;
   } else {
      return SAMRAI_INVALID;
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
ConduitDatabase::getArraySize(
   const std::string& key)
{
   if (!(*d_node).has_child(key) || isDatabase(key)) {
      return 0;
   } else {
      if (isBool(key)) {
         return (*d_node)[key]["data"].dtype().number_of_elements();
      } else if (isDatabaseBox(key)) {
         return (*d_node)[key]["dimension"].dtype().number_of_elements();
      } else if (isString(key)) {
         if ((*d_node)[key].dtype().is_string()) {
            return 1;
         } else {
            return (*d_node)[key].number_of_children();
         }
      } else if (isComplex(key)) {
         return ((*d_node)[key].dtype().number_of_elements() / 2);
      } else {
         return (*d_node)[key].dtype().number_of_elements();
      }
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
ConduitDatabase::isDatabase(
   const std::string& key)
{
   bool is_database = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_DATABASE) {
      is_database = true;
   }
   return is_database;
}

std::shared_ptr<Database>
ConduitDatabase::putDatabase(
   const std::string& key)
{
   deleteKeyIfFound(key);

   (*d_node)[key].set_dtype(conduit::DataType::object());
   d_child_dbs[key] = std::make_shared<ConduitDatabase>(key, &((*d_node)[key]));

   d_types[key] = SAMRAI_DATABASE;

   return d_child_dbs[key];
}

std::shared_ptr<Database>
ConduitDatabase::getDatabase(
   const std::string& key)
{
   if (!isDatabase(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a database...");
   }
   return d_child_dbs[key];
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
ConduitDatabase::isBool(
   const std::string& key)
{
   bool is_bool = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_BOOL) {
      is_bool = true;
   }
   return is_bool;
}

void
ConduitDatabase::putBool(
   const std::string& key,
   const bool& data)
{
   putBoolArray(key, &data, 1);
}

void
ConduitDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   std::vector<conduit::uint8> uint8_vec(nelements, 0);
   for (unsigned int i = 0; i < nelements; ++i) {
      if (data[i]) {
         uint8_vec[i] = 1;
      }
   }
   (*d_node)[key]["data"].set(&(uint8_vec[0]), nelements);
   (*d_node)[key]["bool"] = "true";
   d_types[key] = SAMRAI_BOOL;
}

bool
ConduitDatabase::getBool(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isBool(key) ||
       child["data"].dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a boolean scalar...");
   }
   return static_cast<bool>(child["data"].as_uint8());
}

bool
ConduitDatabase::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   if (d_node->has_child(key)) return getBool(key);

   putBool(key, defaultvalue);
   return defaultvalue;
}

std::vector<bool>
ConduitDatabase::getBoolVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isBool(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a boolean...");
   }
   conduit::uint8_array int_vals = child["data"].as_uint8_array();
   unsigned int vec_size = child["data"].dtype().number_of_elements();
   std::vector<bool> bool_vec(vec_size, false);
   for (unsigned int i = 0; i < vec_size; ++i) {
      if (int_vals[i] != 0) {
         bool_vec[i] = true;
      }
   }
   return bool_vec;
}

void
ConduitDatabase::getBoolArray(
   const std::string& key,
   bool* data,
   const size_t nelements)
{
   std::vector<bool> tmp = getBoolVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
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
ConduitDatabase::isDatabaseBox(
   const std::string& key)
{
   bool is_box = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_BOX) {
      is_box = true;
   }
   return is_box;
}

void
ConduitDatabase::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   putDatabaseBoxArray(key, &data, 1);
}

void
ConduitDatabase::putDatabaseBoxVector(
   const std::string& key,
   const std::vector<DatabaseBox>& data)
{
   putDatabaseBoxArray(key, &data[0], data.size());
}

void
ConduitDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   std::vector<conduit::uint8> dim_vec(nelements);
   std::vector<int> lo_vec(nelements * SAMRAI::MAX_DIM_VAL);
   std::vector<int> hi_vec(nelements * SAMRAI::MAX_DIM_VAL);
   for (unsigned int i = 0; i < nelements; ++i) {
      dim_vec[i] = data[i].getDimVal();

      for (int d = 0; d < dim_vec[i]; ++d) {
         lo_vec[i*SAMRAI::MAX_DIM_VAL + d] = data[i].lower(d);
         hi_vec[i*SAMRAI::MAX_DIM_VAL + d] = data[i].upper(d);
      }
      for (int d = dim_vec[i]; d < SAMRAI::MAX_DIM_VAL; ++d) {
         lo_vec[i*SAMRAI::MAX_DIM_VAL + d] = 0;
         hi_vec[i*SAMRAI::MAX_DIM_VAL + d] = 0;
      }
   }

   (*d_node)[key]["box"] = "true";
   (*d_node)[key]["dimension"].set(&(dim_vec[0]), nelements);
   (*d_node)[key]["lower"].set(&(lo_vec[0]), nelements*SAMRAI::MAX_DIM_VAL);
   (*d_node)[key]["upper"].set(&(hi_vec[0]), nelements*SAMRAI::MAX_DIM_VAL);
   d_types[key] = SAMRAI_BOX;
}

DatabaseBox
ConduitDatabase::getDatabaseBox(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isDatabaseBox(key) ||
       child["dimension"].dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single box...");
   }

   tbox::Dimension dim(child["dimension"].as_uint8());
   int* lower = child["lower"].as_int_ptr(); 
   int* upper = child["upper"].as_int_ptr(); 

   DatabaseBox db_box(dim, lower, upper);

   return db_box;
}

DatabaseBox
ConduitDatabase::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   if (d_node->has_child(key)) return getDatabaseBox(key);

   putDatabaseBox(key, defaultvalue);
   return defaultvalue;
}

std::vector<DatabaseBox>
ConduitDatabase::getDatabaseBoxVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isDatabaseBox(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a DatabaseBox...");
   }
   std::vector<DatabaseBox> box_vec;
   size_t vec_size = child["dimension"].dtype().number_of_elements();
   conduit::uint8_array dim_vals = child["dimension"].as_uint8_array();
   conduit::int_array lo_vals = child["lower"].as_int_array();
   conduit::int_array hi_vals = child["upper"].as_int_array();
   for (size_t i = 0; i < vec_size; ++i) {
      tbox::Dimension dim(dim_vals[i]);
      int * lower = &(lo_vals[i*SAMRAI::MAX_DIM_VAL]);
      int * upper = &(hi_vals[i*SAMRAI::MAX_DIM_VAL]);
      DatabaseBox db_box(dim, lower, upper);
      box_vec.push_back(db_box); 
   }
   return box_vec;
}

void
ConduitDatabase::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const size_t nelements)
{
   std::vector<DatabaseBox> tmp = getDatabaseBoxVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
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
ConduitDatabase::isChar(
   const std::string& key)
{
   bool is_char = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_CHAR) {
      is_char = true;
   }
   return is_char;
}

void
ConduitDatabase::putChar(
   const std::string& key,
   const char& data)
{
   putCharArray(key, &data, 1);
}

void
ConduitDatabase::putCharVector(
   const std::string& key,
   const std::vector<char>& data)
{
   deleteKeyIfFound(key);
   conduit::DataType char_type = conduit::DataType::c_char(data.size());

   char* char_data = new char[data.size()];
   std::memcpy(char_data, &(data[0]), data.size()*sizeof(char));
   (*d_node)[key].set_data_using_dtype(char_type,
                                       static_cast<void*>(char_data));
   delete[] char_data;

   d_types[key] = SAMRAI_CHAR;
}

void
ConduitDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   conduit::DataType char_type = conduit::DataType::c_char(nelements);

   char* char_data = new char[nelements];
   std::memcpy(char_data, data, nelements*sizeof(char));
   (*d_node)[key].set_data_using_dtype(char_type, 
                                       static_cast<void*>(char_data));
   delete[] char_data;

   d_types[key] = SAMRAI_CHAR;
}

char
ConduitDatabase::getChar(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isChar(key) ||
       !child.dtype().is_char() ||
       child.dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single char...");
   }
   return child.as_char();
}

char
ConduitDatabase::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   if (d_node->has_child(key)) return getChar(key);

   putChar(key, defaultvalue);
   return defaultvalue;
}

std::vector<char>
ConduitDatabase::getCharVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isChar(key) || !child.dtype().is_char()) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a char...");
   }
   size_t vec_size = child.dtype().number_of_elements();
   std::vector<char> char_vec(vec_size);
   const char* char_array = static_cast<char *>(child.data_ptr());
   for (size_t i = 0; i < vec_size; ++i) {
      char_vec[i] = char_array[i];
   }
   return char_vec;
}

void
ConduitDatabase::getCharArray(
   const std::string& key,
   char* data,
   const size_t nelements)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isChar(key) || !child.dtype().is_char()) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a char...");
   }
   size_t array_size = child.dtype().number_of_elements();

   if (array_size != nelements) {
      TBOX_CONDUIT_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << array_size << "...");
   }

   const char* char_array = static_cast<char *>(child.data_ptr());
   for (size_t i = 0; i < array_size; ++i) {
      data[i] = char_array[i];
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
ConduitDatabase::isComplex(
   const std::string& key)
{
   bool is_complex = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_COMPLEX) {
      is_complex = true;
   }
   return is_complex;
}

void
ConduitDatabase::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   putComplexArray(key, &data, 1);
}

void
ConduitDatabase::putComplexVector(
   const std::string& key,
   const std::vector<dcomplex>& data)
{
   deleteKeyIfFound(key);
   conduit::DataType cplx_type = conduit::DataType::c_double(2*data.size());

   double* dbl_data = new double[2*data.size()];
   for (size_t i = 0; i < data.size(); ++i) {
      dbl_data[i*2] = data[i].real();
      dbl_data[i*2+1] = data[i].imag();
   } 
   (*d_node)[key].set_data_using_dtype(cplx_type,
                                       static_cast<void*>(dbl_data));
   delete[] dbl_data;

   d_types[key] = SAMRAI_COMPLEX;
}

void
ConduitDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const size_t nelements)
{
   deleteKeyIfFound(key);
   conduit::DataType cplx_type = conduit::DataType::c_double(2*nelements);

   double* dbl_data = new double[2*nelements];
   for (size_t i = 0; i < nelements; ++i) {
      dbl_data[i*2] = data[i].real();
      dbl_data[i*2+1] = data[i].imag();
   }
   (*d_node)[key].set_data_using_dtype(cplx_type,
                                       static_cast<void*>(dbl_data));
   delete[] dbl_data;

   d_types[key] = SAMRAI_COMPLEX;
}

dcomplex
ConduitDatabase::getComplex(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isComplex(key) ||
       !child.dtype().is_double() ||
       child.dtype().number_of_elements() != 2) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single dcomplex...");
   }
   conduit::double_array dbl_vals = child.as_double_array();
   dcomplex rval(dbl_vals[0], dbl_vals[1]);
   return rval;
}

dcomplex
ConduitDatabase::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   if (d_node->has_child(key)) return getComplex(key);

   putComplex(key, defaultvalue);
   return defaultvalue;
}

std::vector<dcomplex>
ConduitDatabase::getComplexVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isComplex(key) || !child.dtype().is_double()) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a dcomplex...");
   }

   size_t data_size = child.dtype().number_of_elements();
   if (data_size % 2 != 0) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " does not have real and imag component for all dcomplex values.");
   }

   std::vector<dcomplex> cplx_vec(data_size/2);
   const double* dbl_array = static_cast<double *>(child.data_ptr());
   for (size_t i = 0; i < data_size/2; ++i) {
      cplx_vec[i] = {dbl_array[i*2], dbl_array[i*2+1]};
   }
   return cplx_vec;
}

void
ConduitDatabase::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const size_t nelements)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isComplex(key) || !child.dtype().is_double()) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a dcomplex...");
   }
   size_t data_size = child.dtype().number_of_elements();
   if (data_size % 2 != 0) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " does not have real and imag component for all dcomplex values.");
   }
   if (data_size != 2*nelements) {
      TBOX_CONDUIT_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << data_size/2 << "...");
   }

   const double* dbl_array = static_cast<double *>(child.data_ptr());
   for (size_t i = 0; i < data_size/2; ++i) {
      data[i] = {dbl_array[i*2], dbl_array[i*2+1]};
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
ConduitDatabase::isDouble(
   const std::string& key)
{
   bool is_double = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_DOUBLE) {
      is_double = true;
   }
   return is_double;
}

void
ConduitDatabase::putDouble(
   const std::string& key,
   const double& data)
{
   putDoubleArray(key, &data, 1);
}

void
ConduitDatabase::putDoubleVector(
   const std::string& key,
   const std::vector<double>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_DOUBLE;
}

void
ConduitDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const size_t nelements)
{
   std::vector<double> dbl_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      dbl_vec[i] = data[i];
   }
   putDoubleVector(key, dbl_vec);
}

double
ConduitDatabase::getDouble(
   const std::string& key)
{
   double value = 0.0;
   conduit::Node& child = getChildNodeOrExit(key);

   if (!isDouble(key) || child.dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single double...");
   }

   const conduit::DataType& datatype = child.dtype();
   if (datatype.is_double()) {
      value = child.as_double();
   } else if (datatype.is_integer()) {
      value = static_cast<double>(child.as_int());
   } else if (datatype.is_float()) {
      value = static_cast<double>(child.as_float());
   } else {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single double...");
   }

   return value;
}

double
ConduitDatabase::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   if (d_node->has_child(key)) return getDouble(key);

   putDouble(key, defaultvalue);
   return defaultvalue;
}

std::vector<double>
ConduitDatabase::getDoubleVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);

   if (!isDouble(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a double...");
   }

   std::vector<double> dbl_vec;

   const conduit::DataType& datatype = child.dtype();
   size_t vec_size = datatype.number_of_elements();
   dbl_vec.resize(vec_size);
   if (datatype.is_double()) {
      conduit::double_array dbl_array = child.as_double_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = dbl_array[i];
      }
   } else if (datatype.is_integer()) {
      conduit::int_array int_vals = child.as_int_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = static_cast<double>(int_vals[i]);
      }
   } else if (datatype.is_float()) {
      conduit::float_array float_array = child.as_float_array();
      for (size_t i = 0; i < vec_size; ++i) {
         dbl_vec[i] = static_cast<double>(float_array[i]);
      }
   } else {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single double...");
   }

   return dbl_vec;
}

void
ConduitDatabase::getDoubleArray(
   const std::string& key,
   double* data,
   const size_t nelements)
{
   std::vector<double> tmp = getDoubleVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
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
ConduitDatabase::isFloat(
   const std::string& key)
{
   bool is_float = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_FLOAT) {
      is_float = true;
   }
   return is_float;
}

void
ConduitDatabase::putFloat(
   const std::string& key,
   const float& data)
{
   putFloatArray(key, &data, 1);
}

void
ConduitDatabase::putFloatVector(
   const std::string& key,
   const std::vector<float>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_FLOAT;
}

void
ConduitDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const size_t nelements)
{
   std::vector<float> flt_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      flt_vec[i] = data[i];
   }

   putFloatVector(key, flt_vec);
}

float
ConduitDatabase::getFloat(
   const std::string& key)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   float value = 0.0;
   conduit::Node& child = getChildNodeOrExit(key);

   if (!isFloat(key) || child.dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single float...");
   }

   const conduit::DataType& datatype = child.dtype();
   if (datatype.is_float()) {
      value = child.as_float();
   } else if (datatype.is_integer()) {
      value = static_cast<float>(child.as_int());
   } else if (datatype.is_double()) {
      value = static_cast<float>(child.as_double());
   } else {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single float...");
   }

   return value;
}

float
ConduitDatabase::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   if (d_node->has_child(key)) return getFloat(key);

   putFloat(key, defaultvalue);
   return defaultvalue;
}

std::vector<float>
ConduitDatabase::getFloatVector(
   const std::string& key)
{
// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   if (!isFloat(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a float...");
   }

   conduit::Node& child = getChildNodeOrExit(key);
   std::vector<float> flt_vec;

   const conduit::DataType& datatype = child.dtype();
   size_t vec_size = datatype.number_of_elements();
   flt_vec.resize(vec_size);
   if (datatype.is_float()) {
      conduit::float_array flt_array = child.as_float_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = flt_array[i];
      }
   } else if (datatype.is_integer()) {
      conduit::int_array int_vals = child.as_int_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = static_cast<float>(int_vals[i]);
      }
   } else if (datatype.is_double()) {
      conduit::double_array dbl_array = child.as_double_array();
      for (size_t i = 0; i < vec_size; ++i) {
         flt_vec[i] = static_cast<float>(dbl_array[i]);
      }
   } else {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single float...");
   }

   return flt_vec;

}

void
ConduitDatabase::getFloatArray(
   const std::string& key,
   float* data,
   const size_t nelements)
{
   std::vector<float> tmp = getFloatVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
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
ConduitDatabase::isInteger(
   const std::string& key)
{
   bool is_int = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_INT) {
      is_int = true;
   }
   return is_int;
}

void
ConduitDatabase::putInteger(
   const std::string& key,
   const int& data)
{
   putIntegerArray(key, &data, 1);
}

void
ConduitDatabase::putIntegerVector(
   const std::string& key,
   const std::vector<int>& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set(data);
   d_types[key] = SAMRAI_INT;
}

void
ConduitDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const size_t nelements)
{
   std::vector<int> int_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      int_vec[i] = data[i];
   }
   putIntegerVector(key, int_vec);
}

int
ConduitDatabase::getInteger(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isInteger(key) ||
       !child.dtype().is_int() ||
       child.dtype().number_of_elements() != 1) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not an integer scalar...");
   }
   return child.as_int();
}

int
ConduitDatabase::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   if (d_node->has_child(key)) return getInteger(key);

   putInteger(key, defaultvalue);
   return defaultvalue;
}

std::vector<int>
ConduitDatabase::getIntegerVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isInteger(key) || !child.dtype().is_int()) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not an integer...");
   }
   size_t vec_size = child.dtype().number_of_elements();
   std::vector<int> int_vec(vec_size);
   conduit::int_array int_vals = child.as_int_array();
   for (size_t i = 0; i < vec_size; ++i) {
      int_vec[i] = int_vals[i];
   }
   return int_vec;
}

void
ConduitDatabase::getIntegerArray(
   const std::string& key,
   int* data,
   const size_t nelements)
{
   std::vector<int> tmp = getIntegerVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
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
ConduitDatabase::isString(
   const std::string& key)
{
   bool is_string = false;
   if (d_node->has_child(key) && d_types[key] == SAMRAI_STRING) {
      is_string = true;
   }
   return is_string;
}

void
ConduitDatabase::putString(
   const std::string& key,
   const std::string& data)
{
   deleteKeyIfFound(key);
   (*d_node)[key].set_string(data);
   d_types[key] = SAMRAI_STRING;
}

void
ConduitDatabase::putStringVector(
   const std::string& key,
   const std::vector<std::string>& data)
{
   deleteKeyIfFound(key);
   int i = 0;
   for (std::vector<std::string>::const_iterator itr = data.begin();
        itr != data.end(); ++itr) {
      std::stringstream ss;
      ss << i;
      std::string id = "str" + ss.str();
      (*d_node)[key][id].set_string(*itr);
      ++i; 
   }
   d_types[key] = SAMRAI_STRING;
}

void
ConduitDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const size_t nelements)
{
   std::vector<std::string> str_vec(nelements);

   for (size_t i = 0; i < nelements; ++i) {
      str_vec[i] = data[i];
   }
   putStringVector(key, str_vec);
}

std::string
ConduitDatabase::getString(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isString(key) ||
       child.has_child("str0")) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a single string ...");
   }
   return child.as_string();
}

std::string
ConduitDatabase::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   if (d_node->has_child(key)) return getString(key);

   putString(key, defaultvalue);
   return defaultvalue;
}

std::vector<std::string>
ConduitDatabase::getStringVector(
   const std::string& key)
{
   conduit::Node& child = getChildNodeOrExit(key);
   if (!isString(key)) {
      TBOX_CONDUIT_DB_ERROR("Key=" << key << " is not a string...");
   }

   std::vector<std::string> str_vec;

   if (child.dtype().is_string()) {
      str_vec.push_back(child.as_string());
   } else {

      size_t nelements = child.number_of_children();

      for (size_t i = 0; i < nelements; ++i) {
         std::stringstream ss;
         ss << i;
         std::string id = "str" + ss.str();
         str_vec.push_back(child[id].as_string());
      } 
   }

   return str_vec;
}

void
ConduitDatabase::getStringArray(
   const std::string& key,
   std::string* data,
   const size_t nelements)
{
   std::vector<std::string> tmp = getStringVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_CONDUIT_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

std::string
ConduitDatabase::getName()
{
   return d_database_name;
}

std::string
ConduitDatabase::getName() const
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

bool ConduitDatabase::deleteKeyIfFound(
   const std::string& key)
{
   if (d_node->has_child(key)) {
      if (isDatabase(key)) {
         d_child_dbs.erase(key);
      }
      d_node->remove(key);
      d_types.erase(key);
      return true;
   } else {
      return false;
   }
}

/*
 *************************************************************************
 *
 * Find the child node associated with the specified key, exit with error
 * if not found.
 *
 *************************************************************************
 */


conduit::Node&
ConduitDatabase::getChildNodeOrExit(
   const std::string& key)
{
   if (!d_node->has_child(key)) {
      TBOX_CONDUIT_DB_ERROR("Key ``" << key << "'' does not exist in the database...");
   }
   return (*d_node)[key];
}

/*
 *************************************************************************
 *
 * Print the entire database to the specified output stream.
 *
 *************************************************************************
 */

void
ConduitDatabase::printClassData(
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
ConduitDatabase::printDatabase(
   std::ostream& os,
   const int indent,
   const int toprint) const
{
   NULL_USE(os);
   NULL_USE(indent);
   NULL_USE(toprint);
// need to implement print
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
