/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 ************************************************************************/

#ifndef included_BaselineData
#define included_BaselineData

#include "SAMRAI/tbox/Database.h"

#include <memory>
#include <string>

namespace testlib {

/**
 * @brief Text-backed regression baseline database for BoxLevels and Connectors.
 *
 * This class is intended for testing baselines of BoxLevels and Connectors
 * for various executables in this test suite.  It is not a
 * general-purpose serialization format and does not have the full features
 * of the Database classes in the tbox directory.
 *
 * The on-disk representation is named @p filename with a `.baseline` suffix.
 */
class BaselineData
{
public:
   explicit BaselineData(const std::string& name = "BaselineData");
   ~BaselineData();

   /**
    * Create an empty text baseline database.
    *
    * The database is written to @p filename with a `.baseline` suffix when
    * close() is called.
    *
    * @param filename Base filename for the baseline.
    * @return True.
    */
   bool create(const std::string& filename);

   /**
    * Open a text baseline database.
    *
    * @param filename Base filename for the baseline.
    * @return True if @p filename with a `.baseline` suffix was parsed
    *         successfully.
    */
   bool open(const std::string& filename);

   /**
    * Write a baseline created with create().
    *
    * @return True if no write is needed or if the baseline was written
    *         successfully.
    */
   bool close();

   /**
    * Return the database being read or written.
    *
    * @return The baseline database, or a null pointer if open() failed.
    */
   std::shared_ptr<SAMRAI::tbox::Database> getDatabase() const;

private:
   static std::string addBaselineSuffix(const std::string& filename);
   void writeDatabase(std::ostream& os,
      const std::shared_ptr<SAMRAI::tbox::Database>& database,
      int indent) const;

   std::string d_name;
   std::string d_filename;
   bool d_writing;
   std::shared_ptr<SAMRAI::tbox::Database> d_database;
};

}

#endif
