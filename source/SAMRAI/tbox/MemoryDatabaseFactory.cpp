/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A factory for building MemoryDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabaseFactory.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/Utilities.h"


namespace SAMRAI {
namespace tbox {

MemoryDatabaseFactory::MemoryDatabaseFactory()
{
}

MemoryDatabaseFactory::~MemoryDatabaseFactory()
{
}

MemoryDatabaseFactory::MemoryDatabaseFactory(
   const MemoryDatabaseFactory& other):
   DatabaseFactory()
{
   NULL_USE(other);
}

MemoryDatabaseFactory&
MemoryDatabaseFactory::operator = (
   const MemoryDatabaseFactory& rhs)
{
   NULL_USE(rhs);
   return *this;
}

/**
 * Build a new MemoryDatabase object.
 */
std::shared_ptr<Database>
MemoryDatabaseFactory::allocate(
   const std::string& name) {
   return std::make_shared<MemoryDatabase>(name);
}

}
}
