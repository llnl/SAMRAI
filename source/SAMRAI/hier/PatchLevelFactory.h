/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/

#ifndef included_hier_PatchLevelFactory
#define included_hier_PatchLevelFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Factory used to create new patch levels.
 *
 * New types of patch level objects can be introduced into SAMRAI by deriving
 * from PatchLevelFactory and re-defining allocate.
 *
 * @see PatchLevel
 */
class PatchLevelFactory
{
public:
   /*!
    * @brief Construct a patch level factory object.
    */
   PatchLevelFactory() = default;

   /*!
    * @brief Virtual destructor for patch level factory objects.
    */
   virtual ~PatchLevelFactory() = default;

   /*!
    * @brief Allocate a patch level with the specified boxes and processor
    * mappings.
    *
    * This method results in the allocated PatchLevel making a COPY of the
    * supplied BoxLevel.  If the caller intends to modify the supplied BoxLevel
    * for other purposes after allocating the new PatchLevel, then this method
    * must be used rather than the method taking a std::shared_ptr<BoxLevel>.
    *
    * Redefine this function to change the method for creating patch levels.
    *
    * @return A std::shared_ptr to the newly created PatchLevel.
    *
    * @param[in]  box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  factory @b Default: a std::shared_ptr to the standard
    *             PatchFactory
    *
    * @pre box_level.getDim() == grid_geometry->getDim()
    */
   virtual std::shared_ptr<PatchLevel>
   allocate(
      const BoxLevel& box_level,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory =
         std::shared_ptr<PatchFactory>()) const;

   /*!
    * @brief Allocate a patch level with the specified boxes and processor
    * mappings.
    *
    * This method results in the allocated PatchLevel ACQUIRING the supplied
    * BoxLevel.  If the caller will not modify the supplied BoxLevel for other
    * purposes after allocating the new PatchLevel, then this method may be
    * used rather than the method taking a BoxLevel&.  Use of this method where
    * permitted is more efficient as it avoids copying an entire BoxLevel.
    * Note that this method results in the supplied BoxLevel being locked so
    * that any attempt by the caller to modify it after calling this method
    * will result in an unrecoverable error.
    *
    * Redefine this function to change the method for creating patch levels.
    *
    * @return A std::shared_ptr to the newly created PatchLevel.
    *
    * @param[in]  box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  factory @b Default: a std::shared_ptr to the standard
    *             PatchFactory
    *
    * @pre box_level.getDim() == grid_geometry->getDim()
    */
   virtual std::shared_ptr<PatchLevel>
   allocate(
      const std::shared_ptr<BoxLevel> box_level,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory =
         std::shared_ptr<PatchFactory>()) const;

   /*!
    * @brief Allocate a patch level using the data from the database to
    * initialize it.
    *
    * Redefine this function to change the method for creating
    * patch levels from a database.
    *
    * @return A std::shared_ptr to the newly created PatchLevel.
    *
    * @param[in]  database
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  factory @b Default: a std::shared_ptr to the standard
    *             PatchFactory
    * @param[in]  defer_boundary_box_creation @b Default: false
    */
   virtual std::shared_ptr<PatchLevel>
   allocate(
      const std::shared_ptr<tbox::Database>& database,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory =
         std::shared_ptr<PatchFactory>(),
      const bool defer_boundary_box_creation = false) const;

private:
   /*
    * Copy constructor and assignment are not implemented.
    */
   PatchLevelFactory(
      const PatchLevelFactory&) = delete;
   PatchLevelFactory&
   operator = (
      const PatchLevelFactory&) = delete;

};

}
}

#endif
