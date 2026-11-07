/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for a DatabaseFactory
 *
 ************************************************************************/

#ifndef included_tbox_DatabaseFactory
#define included_tbox_DatabaseFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief Abstract base class factory used to build Database objects.
 *
 * Used to build database objects.  For example, RestartManager
 * may use a DatabaseFactory to build databases when creating
 * a restart file.
 */
class DatabaseFactory
{
public:
   /**
    * Default constructor
    */
   DatabaseFactory();

   /**
    * Destructor
    */
   virtual ~DatabaseFactory();

   /**
    * Build a new Database instance.
    */
   virtual std::shared_ptr<Database>
   allocate(
      const std::string& name) = 0;

private:
   // Unimplemented copy constructor.
   DatabaseFactory(
      const DatabaseFactory& other);

   // Unimplemented assignment operator.
   DatabaseFactory&
   operator = (
      const DatabaseFactory& rhs);
};

}
}

#endif
