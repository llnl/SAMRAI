/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for a HDFDatabaseFactory
 *
 ************************************************************************/

#include "SAMRAI/tbox/HDFDatabaseFactory.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace tbox {

HDFDatabaseFactory::HDFDatabaseFactory()
{
}

HDFDatabaseFactory::~HDFDatabaseFactory()
{
}

HDFDatabaseFactory::HDFDatabaseFactory(
   const HDFDatabaseFactory& other):
   DatabaseFactory()
{
   NULL_USE(other);
}

HDFDatabaseFactory&
HDFDatabaseFactory::operator = (
   const HDFDatabaseFactory& rhs)
{
   NULL_USE(rhs);
   return *this;
}

/**
 * Build a new Database object.
 */
std::shared_ptr<Database>
HDFDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_HDF5
   std::shared_ptr<HDFDatabase> database(
      std::make_shared<HDFDatabase>(name));
   return database;

#else
   NULL_USE(name);
   TBOX_WARNING("HDF5DatabaseFactory: Cannot allocate an HDFDatabase.\n"
      << "SAMRAI was not configured with HDF." << std::endl);
   return std::shared_ptr<Database>();

#endif
}

}
}
