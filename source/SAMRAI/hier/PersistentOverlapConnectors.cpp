/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Registry of PersistentOverlapConnectorss incident from a common BoxLevel.
 *
 ************************************************************************/
#include "SAMRAI/hier/PersistentOverlapConnectors.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/InputManager.h"

#include <ctype.h>

namespace SAMRAI {
namespace hier {

char PersistentOverlapConnectors::s_check_created_connectors('\0');
char PersistentOverlapConnectors::s_check_accessed_connectors('\0');
bool PersistentOverlapConnectors::s_create_empty_neighbor_containers(false);
char PersistentOverlapConnectors::s_implicit_connector_creation_rule('w');
size_t PersistentOverlapConnectors::s_num_implicit_global_searches(0);

/*
 ************************************************************************
 * This private constructor can only be used by the friend
 * class BoxLevel.
 ************************************************************************
 */
PersistentOverlapConnectors::PersistentOverlapConnectors(
   const BoxLevel& my_box_level):
   d_my_box_level(my_box_level)
{
   getFromInput();
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
PersistentOverlapConnectors::~PersistentOverlapConnectors()
{
   clear();
}

/*
 ************************************************************************
 * Read input parameters.
 ************************************************************************
 */
void
PersistentOverlapConnectors::getFromInput()
{
   if (s_check_created_connectors == '\0') {
      s_check_created_connectors = 'n';
      s_check_accessed_connectors = 'n';
      if (tbox::InputManager::inputDatabaseExists()) {
         std::shared_ptr<tbox::Database> input_db(
            tbox::InputManager::getInputDatabase());
         if (input_db->isDatabase("PersistentOverlapConnectors")) {
            std::shared_ptr<tbox::Database> pocdb(
               input_db->getDatabase("PersistentOverlapConnectors"));

            const bool check_created_connectors(
               pocdb->getBoolWithDefault("DEV_check_created_connectors", false));
            s_check_created_connectors = check_created_connectors ? 'y' : 'n';

            const bool check_accessed_connectors(
               pocdb->getBoolWithDefault("DEV_check_accessed_connectors", false));
            s_check_accessed_connectors =
               check_accessed_connectors ? 'y' : 'n';

            if (pocdb->isString("implicit_connector_creation_rule")) {

               std::string implicit_connector_creation_rule =
                  pocdb->getString("implicit_connector_creation_rule");

               if (implicit_connector_creation_rule != "ERROR" &&
                   implicit_connector_creation_rule != "WARN" &&
                   implicit_connector_creation_rule != "SILENT") {
                  TBOX_ERROR("PersistentOverlapConnectors::getFromInput error:\n"
                     << "implicit_connector_creation_rule must be set to\n"
                     << "\"ERROR\", \"WARN\" or \"SILENT\".\n");
               }

               s_implicit_connector_creation_rule =
                  char(tolower(implicit_connector_creation_rule[0]));
            }
         }
      }
   }
}

/*
 ************************************************************************
 * Create Connector using global search for edges.
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::createConnector(
   const BoxLevel& head,
   const IntVector& connector_width)
{
   TBOX_ASSERT(d_my_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   const size_t num_blocks = head.getRefinementRatio().getNumBlocks();
   IntVector width(connector_width);
   if (width.getNumBlocks() == 1 && num_blocks != 1) {
      if (width.max() == width.min()) {
         width = IntVector(width, num_blocks);
      } else {
         TBOX_ERROR("Anisotropic head width argument for PersistentOverlapConnectors::createConnector must be of size equal to the number of blocks." << std::endl);
      }
   }

   for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {
      if (&d_cons_from_me[i]->getHead() == &head &&
          d_cons_from_me[i]->getConnectorWidth() == width) {
         TBOX_ERROR(
            "PersistentOverlapConnectors::createConnector:\n"
            << "Cannot create duplicate Connectors.");
      }
   }

   std::shared_ptr<Connector> new_connector;
   OverlapConnectorAlgorithm oca;
   oca.findOverlaps(new_connector,
      d_my_box_level,
      head,
      width);

   postprocessForEmptyNeighborContainers(*new_connector);

   d_cons_from_me.push_back(new_connector);
   head.getPersistentOverlapConnectors().d_cons_to_me.push_back(new_connector);

   return *d_cons_from_me.back();
}

/*
 ************************************************************************
 * Create Connector and with transpose using global search for edges.
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::createConnectorWithTranspose(
   const BoxLevel& head,
   const IntVector& connector_width,
   const IntVector& transpose_connector_width)
{
   TBOX_ASSERT(d_my_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   const Connector& forward = createConnector(head, connector_width);
   if (&d_my_box_level != &head) {
      head.createConnector(d_my_box_level, transpose_connector_width);
      d_cons_from_me.back()->setTranspose(
         head.getPersistentOverlapConnectors().d_cons_from_me.back().get(),
         false);
   }

   return forward;
}

/*
 ************************************************************************
 * Cache the user-provided Connector and its transpose if it exists.
 ************************************************************************
 */
void
PersistentOverlapConnectors::cacheConnector(
   std::shared_ptr<Connector>& connector)
{
   TBOX_ASSERT(connector);
   TBOX_ASSERT(d_my_box_level.isInitialized());
   TBOX_ASSERT(d_my_box_level == connector->getBase());

   const BoxLevel& head = connector->getHead();
   doCacheConnectorWork(head, connector);
   if (connector->hasTranspose()) {
      if (connector.get() != &connector->getTranspose()) {
         std::shared_ptr<Connector> transpose(&connector->getTranspose());
         head.getPersistentOverlapConnectors().doCacheConnectorWork(
            d_my_box_level,
            transpose);
      }
      connector->setTranspose(&connector->getTranspose(), false);
   } else {
      connector->setTranspose(0, false);
   }
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::findConnector(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   ConnectorNotFoundAction not_found_action,
   bool exact_width_only)
{
   TBOX_ASSERT(d_my_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   std::shared_ptr<Connector> found = doFindConnectorWork(head,
         min_connector_width,
         not_found_action,
         exact_width_only);
   if (found->hasTranspose()) {
      found->setTranspose(&found->getTranspose(), false);
   } else {
      if (&d_my_box_level == &head) {
         found->setTranspose(found.get(), false);
      } else {
         found->setTranspose(0, false);
      }
   }
   return *found;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
const Connector&
PersistentOverlapConnectors::findConnectorWithTranspose(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   const IntVector& transpose_min_connector_width,
   ConnectorNotFoundAction not_found_action,
   bool exact_width_only)
{
   TBOX_ASSERT(d_my_box_level.isInitialized());
   TBOX_ASSERT(head.isInitialized());

   std::shared_ptr<Connector> forward = doFindConnectorWork(head,
         min_connector_width,
         not_found_action,
         exact_width_only);
   if (&d_my_box_level != &head) {
      std::shared_ptr<Connector> transpose =
         head.getPersistentOverlapConnectors().doFindConnectorWork(
            d_my_box_level,
            transpose_min_connector_width,
            not_found_action,
            exact_width_only);
      forward->setTranspose(transpose.get(), false);
   } else {
      forward->setTranspose(forward.get(), false);
   }

   return *forward;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
bool
PersistentOverlapConnectors::hasConnector(
   const BoxLevel& head,
   const IntVector& min_connector_width) const
{
   const size_t num_blocks = head.getRefinementRatio().getNumBlocks();
   IntVector min_width(min_connector_width);
   if (min_width.getNumBlocks() == 1 && num_blocks != 1) {
      if (min_width.max() == min_width.min()) {
         min_width = IntVector(min_width, num_blocks);
      } else {
         TBOX_ERROR("Anisotropic head width argument for PersistentOverlapConnectors::doFindConnectorWork must be of size equal to the number of blocks." << std::endl);
      }
   }
   for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {
      if (&d_cons_from_me[i]->getHead() == &head &&
          d_cons_from_me[i]->getConnectorWidth() >= min_width) {
         return true;
      }
   }
   return false;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
void
PersistentOverlapConnectors::clear()
{
   if (d_cons_from_me.empty() && d_cons_to_me.empty()) {
      return;
   }

   /*
    * Delete Connectors from me.
    */
   for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {

      const Connector* delete_me = d_cons_from_me[i].get();

      ConVect& cons_at_head =
         delete_me->getHead().getPersistentOverlapConnectors().d_cons_to_me;

      for (ConVect::iterator j = cons_at_head.begin();
           j != cons_at_head.end(); ++j) {
         if (j->get() == delete_me) {
            j->reset();
            cons_at_head.erase(j);
            break;
         }
      }

#ifdef DEBUG_CHECK_ASSERTIONS

      for (int j = 0; j < static_cast<int>(cons_at_head.size()); ++j) {
         TBOX_ASSERT(cons_at_head[j].get() != delete_me);
      }
#endif

      d_cons_from_me[i].reset();
   }
   d_cons_from_me.clear();

   /*
    * Delete Connectors to me.
    */
   for (int i = 0; i < static_cast<int>(d_cons_to_me.size()); ++i) {

      const Connector* delete_me = d_cons_to_me[i].get();

      // Remove reference held by other end of Connector.
      ConVect& cons_at_base =
         delete_me->getBase().getPersistentOverlapConnectors().d_cons_from_me;

      for (ConVect::iterator j = cons_at_base.begin();
           j != cons_at_base.end(); ++j) {
         if (j->get() == delete_me) {
            j->reset();
            cons_at_base.erase(j);
            break;
         }
      }

#ifdef DEBUG_CHECK_ASSERTIONS

      for (int j = 0; j < static_cast<int>(cons_at_base.size()); ++j) {
         TBOX_ASSERT(cons_at_base[j].get() != delete_me);
      }

#endif

      d_cons_to_me[i].reset();

   }
   d_cons_to_me.clear();
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
std::shared_ptr<Connector>
PersistentOverlapConnectors::doFindConnectorWork(
   const BoxLevel& head,
   const IntVector& min_connector_width,
   ConnectorNotFoundAction not_found_action,
   bool exact_width_only)
{
   const size_t num_blocks = head.getRefinementRatio().getNumBlocks();
   IntVector min_width(min_connector_width);
   if (min_width.getNumBlocks() == 1 && num_blocks != 1) {
      if (min_width.max() == min_width.min()) {
         min_width = IntVector(min_width, num_blocks);
      } else {
         TBOX_ERROR("Anisotropic head width argument for PersistentOverlapConnectors::doFindConnectorWork must be of size equal to the number of blocks." << std::endl);
      }
   }

   std::shared_ptr<Connector> found;
   for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());

      if (&(d_cons_from_me[i]->getHead()) == &head) {
         if (d_cons_from_me[i]->getConnectorWidth() >= min_width) {
            if (!found) {
               found = d_cons_from_me[i];
            } else {
               IntVector vdiff =
                  d_cons_from_me[i]->getConnectorWidth()
                  - found->getConnectorWidth();

               TBOX_ASSERT(vdiff != 0);

               size_t nblocks = head.getGridGeometry()->getNumberBlocks();
               int diff = 0;
               for (BlockId::block_t b = 0; b < nblocks; ++b) {
                  for (unsigned int j = 0; j < vdiff.getDim().getValue(); ++j) {
                     diff += vdiff(b,j);
                  }
               }
               if (diff < 0) {
                  found = d_cons_from_me[i];
               }
            }
            if (found->getConnectorWidth() == min_width) {
               break;
            }
         }
      }
   }

   bool fail = false;
   bool warn = false;
   bool create = false;
   if (not_found_action == CONNECTOR_ERROR) {
      fail = true;
   } else if (not_found_action == CONNECTOR_CREATE) {
      create = true;
   } else if (s_implicit_connector_creation_rule == 'e') {
      fail = true;
   } else if (s_implicit_connector_creation_rule == 'w') {
      warn = true;
      create = true;
   } else if (s_implicit_connector_creation_rule == 's') {
      create = true;
   }
   if (!found) {
      if (fail) {
         tbox::perr
         << "PersistentOverlapConnectors::findConnector: Failed to find Connector\n"
         << &d_my_box_level << "--->" << &head
         << " with " << (exact_width_only ? "exact" : "min")
         << " width of " << min_width << ".\n"
         << "base:\n" << d_my_box_level.format("B: ")
         << "head:\n" << head.format("H: ")
         << "\nThe available Connectors have these widths:\n";
         for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {
            if (&(d_cons_from_me[i]->getHead()) == &head) {
               tbox::perr << "\t" << d_cons_from_me[i]->getConnectorWidth()
                          << '\n';
            }
         }
         TBOX_ERROR("To automatically create the missing "
            << "connector, call findConnector with CREATE or call\n"
            << "with IMPLICIT_CREATION_RULE and\n"
            << "implicit_connector_creation_rule = \"WARN\" in the\n"
            << "PersistentOverlapConnectors input database.");
      } else if (create) {
         ++s_num_implicit_global_searches;
         if (warn) {
            TBOX_WARNING("PersistentOverlapConnectors::findConnector is resorting\n"
               << "to a global search to find overlaps between "
               << &d_my_box_level << " and " << &head << ".\n"
               << "This relies on unscalable data or triggers unscalable operations.\n"
               << "Number of implicit global searches: " << s_num_implicit_global_searches << '\n');
         }

         createConnector( head, min_width );
         found = d_cons_from_me.back();
      }
   } else if (exact_width_only &&
              found->getConnectorWidth() != min_width) {

      /*
       * Found a sufficient Connector, but it is too wide.  Extract
       * relevant neighbors from it to make a Connector with the exact
       * width.  This is scalable!
       */

      OverlapConnectorAlgorithm oca;
      std::shared_ptr<Connector> new_connector(std::make_shared<Connector>(
         d_my_box_level,
         head,
         min_width));
      oca.extractNeighbors(*new_connector, *found, min_width);

      postprocessForEmptyNeighborContainers(*new_connector);

      d_cons_from_me.push_back(new_connector);
      head.getPersistentOverlapConnectors().d_cons_to_me.push_back(
         new_connector);

      found = new_connector;

   }

   if (s_check_accessed_connectors == 'y') {
      if (found->checkOverlapCorrectness() != 0) {
         TBOX_ERROR("PersistentOverlapConnectors::findConnector errror:\n"
            << "Bad overlap Connector found.");
      }
   }
   return found;
}

/*
 ************************************************************************
 *
 ************************************************************************
 */
void
PersistentOverlapConnectors::doCacheConnectorWork(
   const BoxLevel& head,
   std::shared_ptr<Connector>& connector)
{
   for (int i = 0; i < static_cast<int>(d_cons_from_me.size()); ++i) {
      TBOX_ASSERT(d_cons_from_me[i]->isFinalized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().isInitialized());
      TBOX_ASSERT(d_cons_from_me[i]->getBase().getBoxLevelHandle());
      TBOX_ASSERT(d_cons_from_me[i]->getHead().getBoxLevelHandle());
      TBOX_ASSERT(&d_cons_from_me[i]->getBase() ==
         &d_cons_from_me[i]->getBase().getBoxLevelHandle()->
         getBoxLevel());
      TBOX_ASSERT(&d_cons_from_me[i]->getHead() ==
         &d_cons_from_me[i]->getHead().getBoxLevelHandle()->
         getBoxLevel());
      if (&(d_cons_from_me[i]->getHead()) == &head &&
          d_cons_from_me[i]->getConnectorWidth() == connector->getConnectorWidth()) {
         TBOX_ERROR(
            "PersistentOverlapConnectors::createConnector:\n"
            << "Cannot create duplicate Connectors.");
      }
   }

   connector->setBase(d_my_box_level);
   connector->setHead(head, true);

   postprocessForEmptyNeighborContainers(*connector);

   if (s_check_created_connectors == 'y') {
      if (connector->checkOverlapCorrectness() != 0) {
         TBOX_ERROR("PersistentOverlapConnectors::cacheConnector errror:\n"
            << "Bad overlap Connector found.");
      }
   }

   d_cons_from_me.push_back(connector);
   head.getPersistentOverlapConnectors().d_cons_to_me.push_back(connector);
}

/*
 ************************************************************************
 ************************************************************************
 */
void
PersistentOverlapConnectors::setCreateEmptyNeighborContainers(
   bool create_empty_neighbor_containers)
{
   s_create_empty_neighbor_containers = create_empty_neighbor_containers;
}

/*
 ************************************************************************
 ************************************************************************
 */
void
PersistentOverlapConnectors::postprocessForEmptyNeighborContainers(
   Connector& connector)
{
   if (s_create_empty_neighbor_containers) {
      const BoxContainer& base_boxes = connector.getBase().getBoxes();
      for (RealBoxConstIterator bi(base_boxes.realBegin());
           bi != base_boxes.realEnd(); ++bi) {
         connector.makeEmptyLocalNeighborhood(bi->getBoxId());
      }
   } else {
      /*
       * Remove empty neighborhood sets.  They are not essential to an
       * overlap Connector.
       */
      connector.eraseEmptyNeighborSets();
   }
}

}
}
