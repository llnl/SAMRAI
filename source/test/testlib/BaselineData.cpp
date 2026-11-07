/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 ************************************************************************/

#include "BaselineData.h"

#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/Parser.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string_view>

namespace {

// Parse a rank-local text baseline into a memory database.
std::shared_ptr<SAMRAI::tbox::Database> parseTextDatabase(
   const std::string& filename)
{
   std::unique_ptr<FILE, decltype(&fclose)> input(
      fopen(filename.c_str(), "r"),
      fclose);
   if (!input) {
      return std::shared_ptr<SAMRAI::tbox::Database>();
   }

   std::shared_ptr<SAMRAI::tbox::MemoryDatabase> database(
      std::make_shared<SAMRAI::tbox::MemoryDatabase>("BaselineData input"));
   SAMRAI::tbox::Parser parser;
   const int errors = parser.parse(filename, input.get(), database, true);
   if (errors != 0) {
      return std::shared_ptr<SAMRAI::tbox::Database>();
   }

   return database;
}

// Write comma-separated values without a trailing separator.
template<class T>
void writeValues(std::ostream& os, const std::vector<T>& values)
{
   for (size_t i = 0; i != values.size(); ++i) {
      if (i != 0) {
         os << ", ";
      }
      os << values[i];
   }
}

}

namespace testlib {

BaselineData::BaselineData(const std::string& name):
   d_name(name),
   d_writing(false)
{
}

BaselineData::~BaselineData()
{
   close();
}

std::string
BaselineData::addBaselineSuffix(const std::string& filename)
{
   constexpr std::string_view suffix(".baseline");
   return filename.ends_with(suffix) ? filename : filename + suffix.data();
}

bool
BaselineData::create(const std::string& filename)
{
   d_filename = addBaselineSuffix(filename);
   d_writing = true;
   d_database = std::make_shared<SAMRAI::tbox::MemoryDatabase>(d_name);
   return true;
}

bool
BaselineData::open(const std::string& filename)
{
   d_writing = false;
   d_filename.clear();
   const std::string text_filename(addBaselineSuffix(filename));
   d_database = parseTextDatabase(text_filename);
   return d_database != nullptr;
}

bool
BaselineData::close()
{
   if (!d_writing || !d_database) {
      return true;
   }
   std::ofstream output(d_filename);
   if (!output.good()) {
      return false;
   }
   output << "// Generated SAMRAI regression baseline.\n";
   writeDatabase(output, d_database, 0);
   output.close();
   d_writing = false;
   return output.good();
}

std::shared_ptr<SAMRAI::tbox::Database>
BaselineData::getDatabase() const
{
   return d_database;
}

void
BaselineData::writeDatabase(
   std::ostream& os,
   const std::shared_ptr<SAMRAI::tbox::Database>& database,
   int indent) const
{
   std::vector<std::string> keys(database->getAllKeys());
   std::sort(keys.begin(), keys.end());
   const std::string pad(static_cast<size_t>(indent), ' ');
   os << std::setprecision(std::numeric_limits<double>::max_digits10);
   for (const std::string& key : keys) {
      const SAMRAI::tbox::Database::DataType type(database->getArrayType(key));
      if (type == SAMRAI::tbox::Database::SAMRAI_DATABASE) {
         os << pad << key << " {\n";
         writeDatabase(os, database->getDatabase(key), indent + 3);
         os << pad << "}\n";
         continue;
      }
      os << pad << key << " = ";
      switch (type) {
         case SAMRAI::tbox::Database::SAMRAI_BOOL:
         {
            const std::vector<bool> values(database->getBoolVector(key));
            for (size_t i = 0; i != values.size(); ++i) {
               os << (i ? ", " : "") << (values[i] ? "TRUE" : "FALSE");
            }
            break;
         }
         case SAMRAI::tbox::Database::SAMRAI_INT:
            writeValues(os, database->getIntegerVector(key));
            break;
         case SAMRAI::tbox::Database::SAMRAI_FLOAT:
            writeValues(os, database->getFloatVector(key));
            break;
         case SAMRAI::tbox::Database::SAMRAI_DOUBLE:
            writeValues(os, database->getDoubleVector(key));
            break;
         case SAMRAI::tbox::Database::SAMRAI_COMPLEX:
            writeValues(os, database->getComplexVector(key));
            break;
         case SAMRAI::tbox::Database::SAMRAI_BOX:
         {
            const std::vector<SAMRAI::tbox::DatabaseBox> values(database->getDatabaseBoxVector(key));
            for (size_t i = 0; i != values.size(); ++i) {
               if (i) {
                  os << ", ";
               }
               os << "[(";
               for (int d = 0; d != values[i].getDimVal(); ++d) {
                  os << (d ? "," : "") << values[i].lower(d);
               }
               os << "),(";
               for (int d = 0; d != values[i].getDimVal(); ++d) {
                  os << (d ? "," : "") << values[i].upper(d);
               }
               os << ")]";
            }
            break;
         }
         default:
            TBOX_ERROR("BaselineData: unsupported database type for key " << key << "\n");
      }
      os << "\n";
   }
}

}
