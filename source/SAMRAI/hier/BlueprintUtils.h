/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Blueprint utilities.
 *
 ************************************************************************/
#ifndef included_hier_BlueprintUtils
#define included_hier_BlueprintUtils

#include "SAMRAI/SAMRAI_config.h"

#ifdef SAMRAI_HAVE_CONDUIT
#include "SAMRAI/hier/BlueprintUtilsStrategy.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for performing common tasks for describing the mesh
 * in the Conduit blueprint format and writing to files.
 *
 * See https://llnl-conduit.readthedocs.io for documentation of the Conduit library.
 *
 * This class can use a BlueprintUtilsStrategy to allow a call-back to
 * application code to fill a blueprint with problem-specific coordinate
 * information.
 *
 * Combined with a BlueprintUtilsStrategy, this loops over a hierarchy,
 * fills Blueprint topology and coordinate entries, calls back to user code
 * for specific coord choices:  uniform, rectilinear, explicit
 */

class PatchHierarchy;
class FlattenedHierarchy;

class BlueprintUtils
{

public:

   /*!
    * @brief Constructor
    *
    * @param strategy    Strategy pointer for callback to application code.
    */
   BlueprintUtils(BlueprintUtilsStrategy* strategy);

   /*!
    * @brief Destructor
    */
   virtual ~BlueprintUtils() = default;

   /*!
    * @brief Put topology and coordinates to the database
    *
    * Using the BlueprintUtilsStrategy given to the constructor of this
    * object, this loops over a hierarchy, fills blueprint topology and
    * and coordinate entries, and calls back to user code to specify the
    * coordinates using the types of coordinates recognized by the
    * blueprint:  uniform, rectilinear, or explicit.
    *
    * @param blueprint_db  Top-level blueprint database holding all local
    *                      domain information
    * @param hierarchy     The hierarchy being described
    * @param topology_name Name of the topology
    */
   void putTopologyAndCoordinatesToDatabase(
      const std::shared_ptr<tbox::Database>& blueprint_db,
      const PatchHierarchy& hierarchy,
      const std::string& topology_name) const;

   /*!
    * @brief Put topology and coordinates to the database
    *
    * Using the BlueprintUtilsStrategy given to the constructor of this
    * object, this loops over a hierarchy, fills blueprint topology and
    * and coordinate entries, and calls back to user code to specify the
    * coordinates using the types of coordinates recognized by the
    * blueprint:  uniform, rectilinear, or explicit.
    *
    * This overloaded version of the method includes a FlattenedHierarchy
    * argument to restrict the filling of coordinates to the finest available
    * level of resolution as represented by the flattened version of the
    * hierarchy.
    *
    * @param blueprint_db    Top-level blueprint database holding all local
    *                        domain information
    * @param hierarchy       The full AMR hierarchy being described
    * @param flat_hierarchy  The flattened version of the AMR hierarchy.  
    * @param topology_name   Name of the topology
    */
   void putTopologyAndCoordinatesToDatabase(
      const std::shared_ptr<tbox::Database>& blueprint_db,
      const PatchHierarchy& hierarchy,
      const FlattenedHierarchy& flat_hierarchy,
      const std::string& topology_name) const;

   /*!
    * @brief Put fields to the database
    *
    * Using the BlueprintUtilsStrategy given to the constructor of this
    * object, this loops over a hierarchy and calls back to user code to
    * add fields data to the mesh in the blueprint_db, according to
    * Blueprint schema defined by conduit.
    *
    * The fields will be associated with the topology identified by the
    * topology_name string.
    *
    * The FlattenedHierarchy argument restricts the filling of field data
    * to the finest level of resolution as represented by the flattened
    * version of the hierarchy.
    *
    * @param blueprint_db    Top-level blueprint database holding all local
    *                        domain information
    * @param hierarchy       The full AMR hierarchy being described
    * @param flat_hierarchy  The flattened version of the AMR hierarchy.  
    * @param topology_name   Name of the topology
    */
   void putFieldsToDatabase(
      const std::shared_ptr<tbox::Database>& blueprint_db,
      const PatchHierarchy& hierarchy,
      const FlattenedHierarchy& flat_hierarchy,
      const std::string& topology_name) const;

   /*!
    * @brief Write blueprint to files.
    *
    * Given a conduit node that holds a blueprint, write it in parallel files.
    * There will be a root file holding the blueprint index and a data files
    * for each rank holding all of the parallel domain data from each rank's
    * local patches.
    *
    * The output files will be {rootfile_name}.root and
    * {data_name}NNNNNN.{io_protocol}, with the Ns representing a six-digit
    * rank number.
    *
    * @param blueprint  Node holding the local blueprint on all ranks
    * @param samrai_mpi SAMRAI_MPI object containing all ranks in use
    * @param num_global_domains  Global number of domains (patches) in the
    *                            hierarchy described by the blueprint
    * @param mesh_name  A name for the mesh
    * @param data_name  File name stub for the files holding the domain data
    * @param rootfile_name  File name sub for the root file.
    * @param io_protocol    I/O protocol string identifier--must be a valid
    *                       protocol as described in Conduit documentation
    */
   void writeBlueprintMesh(
      const conduit::Node& blueprint,
      const tbox::SAMRAI_MPI& samrai_mpi,
      const int num_global_domains,
      const std::string& mesh_name,
      const std::string& data_name,
      const std::string& rootfile_name,
      const std::string& io_protocol) const;

   void setTopologyName(const std::string& topo_name)
   {
      m_topology_name = topo_name;
   }

   const std::string& getTopologyName() const
   {
      return m_topology_name;
   }

private:

   BlueprintUtilsStrategy* d_strategy;

   std::string m_topology_name = "mesh";

};

}
}

#endif // SAMRAI_HAVE_CONDUIT

#endif  // included_hier_BlueprintUtils
