/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Some simple generic database test functions
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"

#include <string>
#include <memory>

using namespace SAMRAI;

// Number of (non-abortive) failures.
extern int number_of_failures;

/**
 * Write database and test contents.
 */
void
setupTestData(
   void);

/**
 * Write database and test contents.
 */
void
writeTestData(
   std::shared_ptr<tbox::Database> db);

/**
 * Read database and test contents.
 */
void
readTestData(
   std::shared_ptr<tbox::Database> db);

/**
 * Test contents of database.
 */
void
testDatabaseContents(
   std::shared_ptr<tbox::Database> db,
   const std::string& tag);
