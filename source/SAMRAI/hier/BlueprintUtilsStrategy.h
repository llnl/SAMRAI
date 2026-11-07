/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   BlueprintUtilsStrategy
 *
 ************************************************************************/
#ifndef included_hier_BlueprintUtilsStrategy
#define included_hier_BlueprintUtilsStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>
#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Strategy class for use by BlueprintUtils
 *
 * This abstract base class can be used by BlueprintUtils for call-backs
 * to application code.
 *
 * @see BlueprintUtils 
 */
class BlueprintUtilsStrategy
{

public:

   /*!
    * @brief Constructor
    */
   BlueprintUtilsStrategy() = default;

   /*!
    * @brief Destructor
    */
   virtual ~BlueprintUtilsStrategy() = default;

   /*!
    * @brief Put blueprint coordinate information into a database
    *
    * This pure virtual function provides an interface to call into application
    * code to put coordinate information into a database that is part of
    * a blueprint description of the mesh.  The coordinates for a single patch
    * should be put into the database, using one of blueprint's recognized
    * coordinate types:  uniform, rectilinear, or explicit.
    *
    * @param coords_db   Database to hold coordinate information
    * @param patch       Patch for which coordinates will be described
    * @param box         Coordinate information should be described for
    *                    intersection of this box and the patch's box
    */
   virtual void putCoordinatesToDatabase(
      std::shared_ptr<tbox::Database>& coords_db,
      const Patch& patch,
      const Box& box) = 0;

   /*!
    * @brief Put blueprint field data into a domain database
    *
    * This virtual function provides an interface to call into application
    * code to add field data to a domain that is part of a blueprint
    * mesh, according to the blueprints format for fields.
    *
    * The default no-op implementation allows for meshes to be created with
    * no fields.
    *
    * @param domain_db   Database for the blueprint mesh domain.  The domain
    *                    should represent the spatial intersection of the
    *                    patch and box arguments
    * @param patch       Patch holding the data.
    * @param box         Only data existing on the intersection of this
    *                    box and the patch's box should be added to the
    *                    domain database
    * @param topology_name   Name of the topology to be associated with the
    *                        fields
    *
    */
   virtual void putFieldsToDomainDatabase(
      std::shared_ptr<tbox::Database>& domain_db,
      const Patch& patch,
      const Box& box,
      const std::string& topology_name)
   {
      NULL_USE(domain_db);
      NULL_USE(patch);
      NULL_USE(box);
      NULL_USE(topology_name);
   }

private:


};

}
}

#endif  // included_hier_BlueprintUtils
