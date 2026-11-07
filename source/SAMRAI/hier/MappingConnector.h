/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Set of distributed box-graph relationships from one BoxLevel
 *                to another describing Box mappings.
 *
 ************************************************************************/
#ifndef included_hier_MappingConnector
#define included_hier_MappingConnector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief A container which holds overlap relationship connections between two
 * BoxLevels.
 */

class MappingConnector:public Connector
{
public:
   /*!
    * @brief Creates an uninitialized MappingConnector object in the
    * distributed state.
    *
    * @param dim The dimension of the head and base BoxLevels that
    * this object will eventually connect.
    *
    * @see setBase()
    * @see setHead()
    * @see setWidth()
    */
   explicit MappingConnector(
      const tbox::Dimension& dim);

   /*!
    * @brief Creates an MappingConnector which is initialized from a restart
    * database.
    *
    * @param dim The dimension of the head and base BoxLevels that
    * this object will eventually connect.
    *
    * @param restart_db Restart Database written by a Connector.
    *
    * @see setBase()
    * @see setHead()
    * @see setWidth()
    */
   MappingConnector(
      const tbox::Dimension& dim,
      tbox::Database& restart_db);

   /*!
    * @brief Copy constructor.
    *
    * @param[in] other
    */
   MappingConnector(
      const MappingConnector& other);

   /*!
    * @brief Initialize a MappingConnector with no defined relationships.
    *
    * The MappingConnector's relationships are initialized to a dummy state.
    *
    * @param[in] base_box_level
    * @param[in] head_box_level
    * @param[in] base_width
    * @param[in] parallel_state
    *
    * @pre (base_box_level.getDim() == head_box_level.getDim()) &&
    *      (base_box_level.getDim() == base_width.getDim())
    */
   MappingConnector(
      const BoxLevel& base_box_level,
      const BoxLevel& head_box_level,
      const IntVector& base_width,
      const BoxLevel::ParallelState parallel_state = BoxLevel::DISTRIBUTED);

   /*!
    * @brief Destructor.
    */
   ~MappingConnector() = default;

   /*!
    * @brief Assignment operator
    */
   MappingConnector&
   operator = (
      const MappingConnector& rhs);

   /*!
    * @brief Create and return this MappingConnector's transpose, assuming that
    * all relationships are local (no remote neighbors).
    *
    * If any remote neighbor is found an unrecoverable assertion is
    * thrown.
    *
    * Non-periodic relationships in are simply reversed to get the transpose
    * relationship.  For each periodic relationship we create a periodic
    * relationship incident from the unshifted head neighbor to the shifted
    * base neighbor.  This is because all relationships must be incident from a
    * real (unshifted) Box.
    */
   MappingConnector *
   createLocalTranspose() const;

   /*!
    * @brief Create and return this MappingConnector's transpose.
    *
    * Similar to createLocalTranspose(), but this method allows
    * non-local edges.  Global data is required, so this method
    * is not scalable.
    */
   virtual MappingConnector *
   createTranspose() const;

   /*!
    * @brief Types of mappings for use in findMappingErrors() and
    *        assertMappingValidity().
    */
   enum MappingType { LOCAL, NOT_LOCAL, UNKNOWN };

   /*!
    * @brief Check if the MappingConnector has a valid mapping.
    *
    * This function can be called prior to calling modify().  It
    * is intended to check the MappingConnector(s) argument(s)
    * in modify() to determine if it has a valid mapping.
    * In other words, it checks to see if the mapping can be used
    * in modify() without logic errors.  It does no other checks.
    *
    * @param[in] map_type LOCAL means assume the mapping is local.  NOT_LOCAL
    * means the mapping is not local.  UNKNOWN means find out whether the
    * map is local or not (communication required) and act
    * accordingly.
    *
    * @return number of errors found.
    */
   size_t
   findMappingErrors(
      MappingType map_type = UNKNOWN) const;

   /*!
    * @brief Run findMappingErrors and abort if any errors are found.
    *
    * @param[in] map_type LOCAL means assume the mapping is local.  NOT_LOCAL
    * means the mapping is not local.  UNKNOWN means find out whether the
    * map is local or not (communication required) and act
    * accordingly.
    */
   void
   assertMappingValidity(
      MappingType map_type = UNKNOWN) const;

};

}
}

#endif // included_hier_MappingConnector
