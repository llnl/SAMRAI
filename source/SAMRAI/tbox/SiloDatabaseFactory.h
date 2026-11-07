/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A factory for building SiloDatabases
 *
 ************************************************************************/

#ifndef included_tbox_SiloDatabaseFactory
#define included_tbox_SiloDatabaseFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/DatabaseFactory.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief SiloDatabase factory.
 *
 * Builds a new SiloDatabase.
 */
class SiloDatabaseFactory:public DatabaseFactory
{
public:
   /**
    * Default constructor.
    */
   SiloDatabaseFactory();

   /**
    * Copy constructor.
    */
   SiloDatabaseFactory(
      const SiloDatabaseFactory& other);

   /**
    * Assignment operator.
    */
   SiloDatabaseFactory&
   operator = (
      const SiloDatabaseFactory& rhs);

   /**
    * Destructor.
    */
   ~SiloDatabaseFactory();

   /**
    * Build a new Database object.
    */
   virtual std::shared_ptr<Database>
   allocate(
      const std::string& name);
};

}
}

#endif
