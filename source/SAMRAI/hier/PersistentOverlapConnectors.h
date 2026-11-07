/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Manager of Connectors incident from a common BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_PersistentOverlapConnectors
#define included_hier_PersistentOverlapConnectors

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Action to take when Connectors between BoxLevels are not found.
 */
enum ConnectorNotFoundAction {
   CONNECTOR_ERROR,                 // If the Connector is not found then error
   CONNECTOR_CREATE,                // If the Connector is not found silently
                                    // create it
   CONNECTOR_IMPLICIT_CREATION_RULE // If the Connector is not found take
                                    // action specified by
                                    // s_implicit_connector_creation_rule
};

class Connector;
class BoxLevel;

/*!
 * @brief A managager of overlap Connectors incident from a
 * BoxLevel, used to store and, if needed, generate overlap
 * Connectors in BoxLevel.
 *
 * PersistantOverlapConnectors provides a mechanism for objects using
 * the BoxLevel to look up overlap Connectors for the
 * BoxLevel.  Connectors created or returned by this class are
 * complete overlap Connectors that always contain the correct overlap
 * neighbor data.
 *
 * For improved scalability, Connectors can be constructed externally
 * and copied into the collection.  Connectors can also be
 * automatically computed using a non-scalable global search.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b>
 *
 *    - \b implicit_connector_creation_rule
 *      How to proceed when findConnector() cannot find any suitable overlap
 *      Connector.  Values can be "ERROR", "WARN" (default) or "SILENT".  If
 *      "SILENT", silently get a globalized version of the head BoxLevel and
 *      look for overlaps.  If "WARN", do the same thing but write a warning to
 *      the log.  If "ERROR", exit with an error.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>implicit_connector_creation_rule</td>
 *     <td>string</td>
 *     <td>"WARN"</td>
 *     <td>"ERROR", "WARN", "SILENT"</td>
 *     <td>opt</td>
 *     <td>Not read from restart</td>
 *   </tr>
 * </table>
 *
 * @note Creating overlap Connectors by global search is not scalable
 * Nevertheless, the default for implicit_connector_creation_rule is "WARN",
 * so that application development need not worry about missing overlap
 * Connectors.  To selectively enable automatic Connector generation, set this
 * input paramter to "ERROR" and call findConnector() with "CREATE" where you
 * are unsure if the Connector has been created.
 *
 * @see findConnector()
 * @see Connector
 */
class PersistentOverlapConnectors
{

public:
   /*!
    * @brief Set whether to create empty neighbor containers when a
    * base Box has no neighbor.
    *
    * Setting a value of true means that all Connectors returned will
    * have a neighbor container for each base box.  False means that
    * base boxes with no neighbors will not have a neighbor container.
    * The default is false.
    */
   static void
   setCreateEmptyNeighborContainers(
      bool create_empty_neighbor_containers);

   /*!
    * @brief Deletes all Connectors to and from this object
    */
   ~PersistentOverlapConnectors();

private:

   /*!
    * @brief Create an overlap Connector, computing relationships by
    * globalizing data.
    *
    * The base will be the BoxLevel that owns this object.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnector(
      const BoxLevel& head,
      const IntVector& connector_width);

   /*!
    * @brief Create an overlap Connector with its transpose, computing
    * relationships by globalizing data.
    *
    * The base will be the BoxLevel that owns this object.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head
    * @param[in] connector_width
    * @param[in] transpose_connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   createConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& connector_width,
      const IntVector& transpose_connector_width);

   /*!
    * @brief Cache the supplied overlap Connector and its transpose
    * if it exists.
    *
    * @param[in] connector
    *
    * @pre connector
    * @pre myBoxLevel().isInitialized()
    * @pre myBoxLevel() == connector->getBase()
    */
   void
   cacheConnector(
      std::shared_ptr<Connector>& connector);

   /*!
    * @brief Find an overlap Connector with the given head and minimum
    * Connector width.  If the specified Connector is not found, take the
    * specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * @par Assertions
    * If no Connector fits the criteria and not_found_action == ERROR, an
    * unrecoverable error will be generated.  If not_found_action == CREATE,
    * the Connector will be generated using an unscalable algorithm.  If
    * not_found_action == IMPLICIT_CREATION_RULE, the behavior will be
    * determined by the @c implicit_connector_creation_rule input parameter.
    * If it is "ERROR", an unrecoverable error will be generated.  If it is
    * "WARN" or "SILENT" the Connector will be generated using an unscalable
    * algorithm and either a warning will be generated or not.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true);

   /*!
    * @brief Find an overlap Connector with its transpose with the given head
    * and minimum Connector widths.  If the specified Connector is not found,
    * take the specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * TODO: The criterion for selecting a single Connector is
    * arbitrary and should be re-examined.
    *
    * @par Assertions
    * If no Connector fits the criteria and not_found_action == ERROR, an
    * unrecoverable error will be generated.  If not_found_action == CREATE,
    * the Connector will be generated using an unscalable algorithm.  If
    * not_found_action == IMPLICIT_CREATION_RULE, the behavior will be
    * determined by the @c implicit_connector_creation_rule input parameter.
    * If it is "ERROR", an unrecoverable error will be generated.  If it is
    * "WARN" or "SILENT" the Connector will be generated using an unscalable
    * algorithm and either a warning will be generated or not.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] transpose_min_connector_width Find the transpose overlap
    *      Connector satisfying this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre myBoxLevel().isInitialized()
    * @pre head.isInitialized()
    */
   const Connector&
   findConnectorWithTranspose(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      const IntVector& transpose_min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true);

   /*!
    * @brief Returns whether the object has overlap
    * Connectors with the given head and minimum Connector
    * width.
    *
    * TODO:  does the following comment mean that this must be called
    * before the call to findConnector?
    *
    * If this returns true, the Connector fitting the specification
    * exists and findConnector() will not throw an assertion.
    *
    * @param[in] head Find the overlap Connector with this specified head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    *
    * @return True if a Connector is found, otherwise false.
    */
   bool
   hasConnector(
      const BoxLevel& head,
      const IntVector& min_connector_width) const;

   /*!
    * @brief Delete stored Connectors.
    */
   void
   clear();

   const BoxLevel&
   myBoxLevel()
   {
      return d_my_box_level;
   }

   // Unimplemented default constructor.
   PersistentOverlapConnectors();

   //@{ @name Methods meant only for BoxLevel to use.

   /*!
    * @brief Constructor, to be called from the BoxLevel
    * allocating the object.
    *
    * @param my_box_level The BoxLevel served by this
    * object.
    */
   explicit PersistentOverlapConnectors(
      const BoxLevel& my_box_level);

   /*
    * Read from the input database.
    */
   void
   getFromInput();

   /*
    * Method which does work of findConnector and findConnectorWithTranspose.
    */
   std::shared_ptr<Connector>
   doFindConnectorWork(
      const BoxLevel& head,
      const IntVector& min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only);

   /*
    * Method which does work of cacheConnector.
    */
   void
   doCacheConnectorWork(
      const BoxLevel& head,
      std::shared_ptr<Connector>& connector);

   /*
    * @brief Make sure all base boxes have a neighbor set or remove
    * empty neighbor sets, depending on
    * s_create_empty_neighbor_containers.
    */
   void
   postprocessForEmptyNeighborContainers(
      Connector& connector);

   //@}

   //@{
   /*!
    * @brief Only BoxLevels are allowed to construct a
    * PersistentOverlapConnectors.
    */
   friend class BoxLevel;
   //@}

   typedef std::vector<std::shared_ptr<Connector> > ConVect;

   /*!
    * @brief Persistent overlap Connectors incident from me.
    */
   ConVect d_cons_from_me;

   /*!
    * @brief Persistent overlap Connectors incident to me.
    */
   ConVect d_cons_to_me;

   /*!
    * @brief Reference to the BoxLevel served by this object.
    */
   const BoxLevel& d_my_box_level;

   /*!
    * @brief Whether to check overlap Connectors when they are created.
    */
   static char s_check_created_connectors;

   /*!
    * @brief Whether to check overlap Connectors when they are accessed.
    */
   static char s_check_accessed_connectors;

   /*!
    * @brief Whether to create empty neighbor containers when a base
    * Box has no neighbor.
    */
   static bool s_create_empty_neighbor_containers;

   /*!
    * @brief Whether to force Connector finding functions to create
    * connectors that are missing.
    *
    * See input parameter implicit_connector_creation_rule.
    */
   static char s_implicit_connector_creation_rule;

   /*
    * @brief Count of how many times we have done implicit global searches.
    */
   static size_t s_num_implicit_global_searches;

};

}
}

#endif // included_hier_PersistentOverlapConnectors
