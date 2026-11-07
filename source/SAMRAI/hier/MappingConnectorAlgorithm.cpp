/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with MappingConnectors.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxContainerUtils.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>

namespace SAMRAI {
namespace hier {

const std::string MappingConnectorAlgorithm::s_default_timer_prefix(
   "hier::MappingConnectorAlgorithm");
std::map<std::string,
         MappingConnectorAlgorithm::TimerStruct> MappingConnectorAlgorithm::s_static_timers;
char MappingConnectorAlgorithm::s_ignore_external_timer_prefix('n');

char MappingConnectorAlgorithm::s_print_steps = '\0';

const std::string MappingConnectorAlgorithm::s_dbgbord;

int MappingConnectorAlgorithm::s_operation_mpi_tag = 0;
/*
 * Do we even need to use different tags each time we modify???
 * Unique tags were used to help debug, but the methods may work
 * with reused tags anyway.
 */

tbox::StartupShutdownManager::Handler
MappingConnectorAlgorithm::s_initialize_finalize_handler(
   MappingConnectorAlgorithm::initializeCallback,
   0,
   0,
   MappingConnectorAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnectorAlgorithm::MappingConnectorAlgorithm():
   d_mpi(MPI_COMM_NULL),
   d_mpi_is_exclusive(false),
   d_object_timers(NULL),
   d_barrier_before_communication(false),
   d_sanity_check_inputs(false),
   d_sanity_check_outputs(false)
{
   getFromInput();
   setTimerPrefix(s_default_timer_prefix);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

MappingConnectorAlgorithm::~MappingConnectorAlgorithm()
{
   if (d_mpi_is_exclusive) {
      d_mpi.freeCommunicator();
      d_mpi_is_exclusive = false;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void MappingConnectorAlgorithm::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& mpi,
   bool make_duplicate)
{
   if (d_mpi_is_exclusive) {
      d_mpi.freeCommunicator();
      d_mpi_is_exclusive = false;
   }
   if (make_duplicate) {
      d_mpi.dupCommunicator(mpi);
      d_mpi_is_exclusive = true;
   } else {
      d_mpi = mpi;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
MappingConnectorAlgorithm::getFromInput()
{
   if (s_print_steps == '\0') {
      s_print_steps = 'n';
      if (tbox::InputManager::inputDatabaseExists()) {
         std::shared_ptr<tbox::Database> idb(
            tbox::InputManager::getInputDatabase());
         if (idb->isDatabase("MappingConnectorAlgorithm")) {
            std::shared_ptr<tbox::Database> mca_db(
               idb->getDatabase("MappingConnectorAlgorithm"));
            s_print_steps =
               mca_db->getCharWithDefault("DEV_print_modify_steps", 'n');
            if (!(s_print_steps == 'n' || s_print_steps == 'y')) {
               INPUT_VALUE_ERROR("DEV_print_modify_steps");
            }
            s_ignore_external_timer_prefix =
               mca_db->getCharWithDefault("DEV_ignore_external_timer_prefix",
                  'n');
            if (!(s_ignore_external_timer_prefix == 'n' ||
                  s_ignore_external_timer_prefix == 'y')) {
               INPUT_VALUE_ERROR("DEV_ignore_external_timer_prefix");
            }
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
MappingConnectorAlgorithm::modify(
   Connector& anchor_to_mapped,
   const MappingConnector& old_to_new,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{
   d_object_timers->t_modify_public->start();

   Connector* old_to_anchor = 0;
   if (anchor_to_mapped.hasTranspose()) {
      old_to_anchor = &anchor_to_mapped.getTranspose();
   } else {
      old_to_anchor = anchor_to_mapped.createLocalTranspose();
   }
   Connector& mapped_to_anchor = *old_to_anchor;

   const MappingConnector* new_to_old = 0;
   if (old_to_new.hasTranspose()) {
      new_to_old =
         static_cast<MappingConnector *>(&old_to_new.getTranspose());
   }

   /*
    * Ensure that head and base BoxLevels in argument agree with each
    * other and that transpose Connectors are really transposes.
    */
   const Connector& anchor_to_old = anchor_to_mapped;

   const BoxLevel* old = &old_to_new.getBase();

   if (d_sanity_check_inputs) {
      if (!d_mpi.hasNullCommunicator() && !d_mpi.isCongruentWith(old->getMPI())) {
         TBOX_ERROR("MappingConnectorAlgorithm::modify input error: Input BoxLevel\n"
            << "has SAMRAI_MPI that is incongruent with MappingConnectorAlgorithm's.\n"
            << "See MappingConnectorAlgorithm::setSAMRAI_MPI.\n");
      }
   }

   if ((new_to_old && (old != &new_to_old->getHead())) ||
       old != &anchor_to_mapped.getHead() ||
       old != &mapped_to_anchor.getBase()) {
      if (new_to_old) {
         TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
            << "Given Connectors to base and head of modify are not incident\n"
            << "from the same old in MappingConnectorAlgorithm::modify:\n"
            << "anchor_to_old is  TO  " << &anchor_to_old.getHead() << "\n"
            << "old_to_new is FROM " << &old_to_new.getBase()
            << "\n"
            << "new_to_old is  TO  " << &new_to_old->getHead()
            << "\n"
            << "old_to_anchor is FROM " << &old_to_anchor->getBase() << "\n");
      } else {
         TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
            << "Given Connectors to base and head of modify are not incident\n"
            << "from the same old in MappingConnectorAlgorithm::modify:\n"
            << "anchor_to_old is  TO  " << &anchor_to_old.getHead() << "\n"
            << "old_to_new is FROM " << &old_to_new.getBase()
            << "\n"
            << "old_to_anchor is FROM " << &old_to_anchor->getBase() << "\n");
      }
   }
   if (&anchor_to_old.getBase() != &old_to_anchor->getHead()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from anchor of modify do not refer\n"
         << "to the same BoxLevel:\n"
         << "anchor_to_old is FROM " << &anchor_to_old.getBase() << "\n"
         << "old_to_anchor is  TO  " << &old_to_anchor->getHead() << "\n");
   }
   if (new_to_old && (&old_to_new.getHead() != &new_to_old->getBase())) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from new of modify do not refer\n"
         << "to the same BoxLevel:\n"
         << "new_to_old is FROM " << &new_to_old->getBase() << "\n"
         << "old_to_new is  TO  " << &old_to_new.getHead() << "\n");
   }
   if (!anchor_to_old.isTransposeOf(*old_to_anchor)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between anchor and old of modify\n"
         << "are not transposes of each other.\n"
         << "See MappingConnectorAlgorithm::isTransposeOf().\n");
   }
   if (new_to_old && !new_to_old->isTransposeOf(old_to_new)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between new and old of modify\n"
         << "are not transposes of each other.\n"
         << "See MappingConnectorAlgorithm::isTransposeOf().\n");
   }
   if (anchor_to_old.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

   if (s_print_steps == 'y') {
      tbox::plog
      << "MappingConnectorAlgorithm::modify: old box_level:\n"
      << old_to_new.getBase().format(s_dbgbord, 2);
      if (new_to_old) {
         tbox::plog
         << "MappingConnectorAlgorithm::modify: new box_level:\n"
         << new_to_old->getBase().format(s_dbgbord, 2);
      }
      tbox::plog
      << "MappingConnectorAlgorithm::modify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 2);
      if (new_to_old) {
         tbox::plog
         << "MappingConnectorAlgorithm::modify: new_to_old:\n"
         << new_to_old->format(s_dbgbord, 2);
      }
   }

   if (0) {
      // Expensive input checking.
      const BoxLevel& anchor_box_level = anchor_to_old.getBase();
      const BoxLevel& old_box_level = old_to_new.getBase();

      tbox::plog
      << "anchor box_level:\n" << anchor_box_level.format(s_dbgbord, 2)
      << "anchor_to_old:\n" << anchor_to_old.format(s_dbgbord, 2)
      << "old box_level:\n" << old_box_level.format(s_dbgbord, 2)
      << "old_to_new:\n" << old_to_new.format(s_dbgbord, 2);
      if (new_to_old) {
         tbox::plog
         << "new_to_old:\n" << new_to_old->format(s_dbgbord, 2);
      }

      TBOX_ASSERT(anchor_to_old.checkOverlapCorrectness() == 0);
      TBOX_ASSERT(old_to_anchor->checkOverlapCorrectness() == 0);
      TBOX_ASSERT(old_to_anchor->checkTransposeCorrectness(anchor_to_old,
            true) == 0);
      TBOX_ASSERT(old_to_new.checkOverlapCorrectness(true, false) == 0);
      TBOX_ASSERT(!new_to_old ||
                  new_to_old->checkOverlapCorrectness(true, false) == 0);
      TBOX_ASSERT(!new_to_old ||
                  old_to_new.checkTransposeCorrectness(*new_to_old, true) == 0);
   }

   d_object_timers->t_modify_public->stop();
   privateModify(anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      new_to_old,
      mutable_new,
      mutable_old);
   d_object_timers->t_modify_public->start();

   if (d_sanity_check_outputs) {
      anchor_to_mapped.assertTransposeCorrectness(mapped_to_anchor);
      mapped_to_anchor.assertTransposeCorrectness(anchor_to_mapped);
   }

   if (!anchor_to_mapped.hasTranspose()) {
      delete old_to_anchor;
   }

   d_object_timers->t_modify_public->stop();
}

/*
 ***********************************************************************
 * This method modifies overlap Connectors based on the described
 * changes to their base or head BoxLevels.  The change is
 * described as a mapping from the old state to the new.
 *
 * Essential nomenclature:
 * - mapped: The BoxLevel that is being changed.
 * - anchor: A BoxLevel that is NOT being changed.
 *   This method modifies the overlap Connectors between anchor
 *   and mapped.
 * - old: The state of mapped before the change.
 * - new: the state of mapped after the change.
 * - old_to_new: Desription of the change.  The NeighborSet of
 *   an old Box is what the old Box will become.
 *   By convention, the NeighborSet of a un-changing Box
 *   is not required.  However, an empty NeighborSet means that
 *   the old Box will disappear.
 *
 * While modify adjusts the Connector to reflect changes in the mapped
 * BoxLevel, it does NOT manipulate any BoxLevel objects
 * (other than initializing one with another; see in-place changes).
 *
 * The Connector width of the mapping must be at least equal to the
 * ammount that new boxes protrude from their old boxes.  The protusion
 * may generate undetected overlaps between anchor and mapped.  To
 * avoid generating undetected overlaps, the width of anchor<==>mapped
 * are shrunken by the width of the mapping.
 *
 * At this time, new_to_old is only used to determine the remote
 * owners that must be notified of a local box being mapped.
 * The Connector object contains more info than required, so it can be
 * replaced by something more concise.  However, requiring a transpose
 * Connector lets us use some canned sanity checks.  The use of the
 * transpose Connector is fairly efficient.
 ***********************************************************************
 */

void
MappingConnectorAlgorithm::privateModify(
   Connector& anchor_to_mapped,
   Connector& mapped_to_anchor,
   const MappingConnector& old_to_new,
   const MappingConnector* new_to_old,
   BoxLevel* mutable_new,
   BoxLevel* mutable_old) const
{
   const tbox::SAMRAI_MPI& mpi = d_mpi.hasNullCommunicator() ?
      old_to_new.getBase().getMPI() : d_mpi;

   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Errant message detected.");
   }

   if (d_barrier_before_communication) {
      mpi.Barrier();
   }
   d_object_timers->t_modify->start();
   d_object_timers->t_modify_misc->start();

   if (s_print_steps == 'y') {
      tbox::plog
      << "MappingConnectorAlgorithm::privateModify: anchor_to_old:\n"
      << anchor_to_mapped.format(s_dbgbord, 3)
      << "MappingConnectorAlgorithm::privateModify: mapped_to_anchor:\n"
      << mapped_to_anchor.format(s_dbgbord, 3)
      << "MappingConnectorAlgorithm::privateModify: old_to_new:\n"
      << old_to_new.format(s_dbgbord, 3);
      if (new_to_old) {
         tbox::plog
         << "MappingConnectorAlgorithm::privateModify: new_to_old:\n"
         << new_to_old->format(s_dbgbord, 3);
      }
   }

   privateModify_checkParameters(
      anchor_to_mapped,
      mapped_to_anchor,
      old_to_new,
      new_to_old);

   /*
    * anchor<==>mapped start out as Connectors between
    * old and anchor.  Make copies of Connectors to and from old
    * so we can modify anchor_to_mapped/mapped_to_anchor without
    * losing needed data.
    */
   const Connector anchor_to_old = anchor_to_mapped;
   const Connector old_to_anchor = mapped_to_anchor;

   /*
    * We will modify the mapped BoxLevel to make it the new
    * BoxLevel.  "mapped" is the name for the changing
    * BoxLevel.  Because this BoxLevel will be
    * modified for output, it is synonymous with "new".
    */
   Connector& anchor_to_new = anchor_to_mapped;
   Connector& new_to_anchor = mapped_to_anchor;

   /*
    * Shorthand for the three BoxLevels and their Refinement
    * ratios.
    */
   const BoxLevel& old = old_to_anchor.getBase();
   const BoxLevel& anchor = anchor_to_old.getBase();
   const BoxLevel& new_level = old_to_new.getHead();
   const IntVector& old_ratio = old.getRefinementRatio();
   const IntVector& anchor_ratio = anchor.getRefinementRatio();
   const IntVector& new_ratio = new_level.getRefinementRatio();

   const tbox::Dimension dim(old.getDim());

   /*
    * The width of old-->new indicates the maximum amount of box
    * growth caused by the change.  A value of zero means no growth.
    *
    * Growing old Boxes may generate new overlaps that we cannot
    * detect because they lie outside the current anchor<==>mapped
    * widths.  To reflect that we do not see any farther just because
    * the Boxes have grown, we shrink the widths by (nominally)
    * the amount of growth.  To ensure the shrinkage is consistent
    * between transpose pairs of Connectors, it is converted to the
    * coarser index space then converted into the final index space.
    * Calls to Connector::convertHeadWidthToBase() perform the
    * conversions.
    */
   const IntVector shrinkage_in_new_index_space =
      Connector::convertHeadWidthToBase(
         old_ratio,
         new_ratio,
         old_to_new.getConnectorWidth());
   const IntVector shrinkage_in_anchor_index_space =
      Connector::convertHeadWidthToBase(
         anchor_ratio,
         new_ratio,
         shrinkage_in_new_index_space);
   const IntVector anchor_to_new_width =
      anchor_to_old.getConnectorWidth() - shrinkage_in_anchor_index_space;
   if (!(anchor_to_new_width >= IntVector::getZero(dim))) {
      TBOX_ERROR(
         "MappingConnectorAlgorithm::privateModify error:\n"
         << "Mapping connector allows mapped BoxLevel to grow\n"
         << "too much.  The growth may generate new overlaps\n"
         << "that cannot be detected by mapping algorithm, thus\n"
         << "causing output overlap Connectors to be incomplete.\n"
         << "MappingConnector:\n" << old_to_new.format("", 0)
         << "anchor--->mapped:\n" << anchor_to_new.format("", 0)
         << "Connector width of anchor--->mapped will shrink\n"
         << "by " << shrinkage_in_anchor_index_space << " which\n"
         << "will result in a non-positive width." << std::endl);
   }
   const IntVector new_to_anchor_width =
      Connector::convertHeadWidthToBase(
         new_ratio,
         anchor.getRefinementRatio(),
         anchor_to_new_width);

   anchor_to_new.shrinkWidth(anchor_to_new_width);
   new_to_anchor.shrinkWidth(new_to_anchor_width);

   d_object_timers->t_modify_misc->stop();

   /*
    * The essential modify algorithm is in this block.
    */

   d_object_timers->t_modify_misc->start();

   /*
    * Initialize the output connectors with the correct new
    * BoxLevel.  (As inputs, they were referencing the old
    * BoxLevel.)
    */

   /*
    * Determine which ranks we have to communicate with.  They are
    * the ones who owns Boxes that the local Boxes will
    * be mapped to.
    */
   std::set<int> incoming_ranks, outgoing_ranks;
   old_to_anchor.getLocalOwners(outgoing_ranks);
   anchor_to_old.getLocalOwners(incoming_ranks);
   old_to_new.getLocalOwners(outgoing_ranks);
   if (new_to_old) {
      new_to_old->getLocalOwners(incoming_ranks);
   }

   // We don't need to communicate locally.
   incoming_ranks.erase(mpi.getRank());
   outgoing_ranks.erase(mpi.getRank());

   /*
    * visible_anchor_nabrs, visible_new_nabrs are the neighbors that
    * are seen by the local process.  For communication efficiency, we
    * use a looping construct that assumes that they are sorted by
    * owners first.  Note the comparator BoxOwnerFirst used to
    * achieve this ordering.
    */
   bool ordered = true;
   BoxContainer visible_anchor_nabrs(ordered), visible_new_nabrs(ordered);
   InvertedNeighborhoodSet anchor_eto_old, new_eto_old;
   for (Connector::ConstNeighborhoodIterator ei = old_to_anchor.begin();
        ei != old_to_anchor.end(); ++ei) {
      const BoxId& old_gid = *ei;
      for (Connector::ConstNeighborIterator na = old_to_anchor.begin(ei);
           na != old_to_anchor.end(ei); ++na) {
         visible_anchor_nabrs.insert(*na);
         if (old_to_new.hasNeighborSet(old_gid)) {
            /*
             * anchor_eto_old is an InvertedNeighborhoodSet mapping visible anchor
             * Boxes to local old Boxes that are changing (excludes
             * old Boxes that do not change).
             */
            anchor_eto_old[*na].insert(old_gid);
         }
      }
   }
   for (Connector::ConstNeighborhoodIterator ei = old_to_new.begin();
        ei != old_to_new.end(); ++ei) {
      const BoxId& old_gid = *ei;
      for (Connector::ConstNeighborIterator na = old_to_new.begin(ei);
           na != old_to_new.end(ei); ++na) {
         visible_new_nabrs.insert(visible_new_nabrs.end(), *na);
         new_eto_old[*na].insert(old_gid);
      }
   }

   /*
    * Object for communicating relationship changes.
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int> * all_comms(0);

   d_object_timers->t_modify_misc->stop();

   /*
    * Set up communication mechanism (and post receives).
    */
   d_object_timers->t_modify_setup_comm->start();

   s_operation_mpi_tag = 0;

   setupCommunication(
      all_comms,
      comm_stage,
      mpi,
      incoming_ranks,
      outgoing_ranks,
      d_object_timers->t_modify_MPI_wait,
      s_operation_mpi_tag,
      s_print_steps == 'y');

   d_object_timers->t_modify_setup_comm->stop();

   /*
    * There are three major parts to computing the new neighbors:
    * (1) remove relationships for Boxes being mapped into
    * something else, (2) discover new relationships for the new
    * Boxes and (3) share information with other processes.
    *
    * In steps 1 and 2, the owners of the Boxes being
    * modified determine which relationships to remove and which to
    * add.  Some of this information is kept locally while the rest
    * are to be sent to the owners of the affected anchor and new
    * Boxes.
    *
    * The three parts are done in the 3 steps following.  Note that
    * these steps do not correspond to the parts.
    */

   /*
    * Messages for other processors describing removed and added relationships.
    */
   std::map<int, std::vector<int> > send_mesgs;
   for (std::set<int>::const_iterator itr(outgoing_ranks.begin());
        itr != outgoing_ranks.end(); ++itr) {
      send_mesgs[*itr];
   }

   /*
    * First step: Remove neighbor data for Boxes that are
    * going away and cache information to be sent out.
    */
   privateModify_removeAndCache(
      send_mesgs,
      anchor_to_new,
      &new_to_anchor,
      old_to_new);

   /*
    * Second step: Discover overlaps for new Boxes.  Send
    * messages with information about what relationships to remove
    * or create.
    */
   privateModify_discoverAndSend(
      send_mesgs,
      anchor_to_new,
      &new_to_anchor,
      incoming_ranks,
      outgoing_ranks,
      all_comms,
      visible_new_nabrs,
      visible_anchor_nabrs,
      anchor_eto_old,
      new_eto_old,
      old_to_anchor,
      anchor_to_old,
      old_to_new);

   /*
    * Third step: Receive and unpack messages from incoming_ranks.
    * These message contain information about what relationships to
    * remove or add.
    */
   receiveAndUnpack(
      anchor_to_new,
      &new_to_anchor,
      incoming_ranks,
      all_comms,
      comm_stage,
      d_object_timers->t_modify_receive_and_unpack,
      s_print_steps == 'y');

   d_object_timers->t_modify_misc->start();

   if (all_comms) {
      delete[] all_comms;
   }

   /*
    * Now we have set up the NeighborhoodSets for anchor_to_new and
    * new_to_anchor so we can initialize these Connectors.
    */
   anchor_to_new.setHead(new_level);
   anchor_to_new.setWidth(anchor_to_new_width, true);
   new_to_anchor.setBase(new_level);
   new_to_anchor.setHead(anchor);
   new_to_anchor.setWidth(new_to_anchor_width, true);

   if (!anchor_to_new.ratioIsExact()) {
      TBOX_WARNING("MappingConnectorAlgorithm::privateModify: generated\n"
         << "overlap Connectors with non-integer ratio between\n"
         << "the base and head.  The results are not guaranteed\n"
         << "to be complete overlap Connectors." << std::endl);
   }
   d_object_timers->t_modify_misc->stop();

   /*
    * Optional in-place changes:
    *
    * Note that the old and new BoxLevels gotten from Connectors
    * are const, so this method cannot modify them.  Only by specifying
    * the mutable BoxLevels, can this method modify them.
    *
    * If users provide mutable object to initialize to the new
    * BoxLevel, this method initializes it to the new
    * BoxLevel and uses it in the output Connectors.
    */
   d_object_timers->t_modify_misc->start();
   if (mutable_new == &old_to_new.getBase() &&
       mutable_old == &old_to_new.getHead()) {
      /*
       * Since mutable_new is old and mutable_old is new, shortcut
       * two assignments by swapping.
       */
      BoxLevel::swap(*mutable_new, *mutable_old);
      new_to_anchor.setBase(*mutable_new, true);
      anchor_to_new.setHead(*mutable_new, true);
   } else {
      if (mutable_new != 0) {
         *mutable_new = old_to_new.getHead();
         new_to_anchor.setBase(*mutable_new, true);
         anchor_to_new.setHead(*mutable_new, true);
      }
      if (mutable_old != 0) {
         *mutable_old = old_to_new.getBase();
      }
   }
   d_object_timers->t_modify_misc->stop();

   d_object_timers->t_modify->stop();

   if (mpi.hasReceivableMessage(0, MPI_ANY_SOURCE, MPI_ANY_TAG)) {
      TBOX_ERROR("Errant message detected.");
   }
}

/*
 ***********************************************************************
 * Do some standard (not expensive) checks on the arguments of modify.
 ***********************************************************************
 */

void
MappingConnectorAlgorithm::privateModify_checkParameters(
   const Connector& anchor_to_mapped,
   const Connector& mapped_to_anchor,
   const MappingConnector& old_to_new,
   const MappingConnector* new_to_old) const
{
   const BoxLevel& old = mapped_to_anchor.getBase();

   /*
    * Ensure that Connectors incident to and from the old agree on
    * what the old is.
    */
   if ((&old != &old_to_new.getBase()) ||
       (new_to_old && new_to_old->isFinalized() &&
        (&old != &new_to_old->getHead())) ||
       (&old != &anchor_to_mapped.getHead())) {
      if (new_to_old) {
         TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
            << "Given Connectors to anchor and new of modify are not incident\n"
            << "from the same old in MappingConnectorAlgorithm::modify:\n"
            << "anchor_to_mapped is  TO  " << &anchor_to_mapped.getHead() << "\n"
            << "old_to_new is FROM " << &old_to_new.getBase()
            << "\n"
            << "new_to_old is  TO  " << &new_to_old->getHead()
            << "\n"
            << "mapped_to_anchor is FROM " << &mapped_to_anchor.getBase() << "\n");
      } else {
         TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
            << "Given Connectors to anchor and new of modify are not incident\n"
            << "from the same old in MappingConnectorAlgorithm::modify:\n"
            << "anchor_to_mapped is  TO  " << &anchor_to_mapped.getHead() << "\n"
            << "old_to_new is FROM " << &old_to_new.getBase()
            << "\n"
            << "mapped_to_anchor is FROM " << &mapped_to_anchor.getBase() << "\n");
      }
   }
   /*
    * Ensure that new and anchor box_levels in argument agree with
    * new and anchor in the object.
    */
   if (&mapped_to_anchor.getHead() != &anchor_to_mapped.getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from anchor of modify do not refer\n"
         << "to the anchor of the modify in MappingConnectorAlgorithm::modify:\n"
         << "anchor_to_mapped is FROM " << &anchor_to_mapped.getBase() << "\n"
         << "mapped_to_anchor is  TO  " << &mapped_to_anchor.getHead() << "\n"
         << "anchor of modify is    " << &anchor_to_mapped.getBase() << "\n");
   }
   if (new_to_old && new_to_old->isFinalized() &&
       &old_to_new.getHead() != &new_to_old->getBase()) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors to and from new of modify do not refer\n"
         << "to the new of the modify in MappingConnectorAlgorithm::modify:\n"
         << "new_to_old is FROM " << &new_to_old->getBase()
         << "\n"
         << "old_to_new is  TO  " << &old_to_new.getHead()
         << "\n"
         << "new of modify is    " << &anchor_to_mapped.getHead() << "\n");
   }
   if (!anchor_to_mapped.isTransposeOf(mapped_to_anchor)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between anchor and mapped of modify\n"
         << "are not transposes of each other.\n"
         << "See Connector::isTransposeOf().\n");
   }
   if (new_to_old && new_to_old->isFinalized() &&
       !new_to_old->isTransposeOf(old_to_new)) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "Given Connectors between new and old of modify\n"
         << "are not transposes of each other.\n"
         << "See Connector::isTransposeOf().\n");
   }
   if (anchor_to_mapped.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }
   if (mapped_to_anchor.getParallelState() != BoxLevel::DISTRIBUTED) {
      TBOX_ERROR("Bad input for MappingConnectorAlgorithm::modify:\n"
         << "bridging is currently set up for DISTRIBUTED\n"
         << "mode only.\n");
   }

   // Expensive sanity checks:
   if (d_sanity_check_inputs) {
      anchor_to_mapped.assertTransposeCorrectness(mapped_to_anchor);
      mapped_to_anchor.assertTransposeCorrectness(anchor_to_mapped);
      if (new_to_old && new_to_old->isFinalized()) {
         /*
          * Not sure if the following are valid checks for modify operation.
          * Modify *may* have different restrictions on the mapping Connector.
          */
         new_to_old->assertTransposeCorrectness(old_to_new);
         old_to_new.assertTransposeCorrectness(*new_to_old);
      }
      size_t nerrs = old_to_new.findMappingErrors();
      if (nerrs != 0) {
         TBOX_ERROR("MappingConnectorUtil::privateModify: found errors in\n"
            << "mapping Connector.\n"
            << "old:\n" << old_to_new.getBase().format("OLD: ")
            << "new:\n" << old_to_new.getHead().format("NEW: ")
            << "old_to_new:\n" << old_to_new.format("O->N: ")
            << std::endl);
      }
   }
}

/*
 ***********************************************************************
 * Remove relationships made obsolete by mapping.  Cache outgoing
 * information in message buffers.
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::privateModify_removeAndCache(
   std::map<int, std::vector<int> >& send_mesgs,
   Connector& anchor_to_new,
   Connector* new_to_anchor,
   const MappingConnector& old_to_new) const
{
   d_object_timers->t_modify_remove_and_cache->start();

   const tbox::Dimension& dim(old_to_new.getBase().getDim());
   const tbox::SAMRAI_MPI& mpi = d_mpi.getCommunicator() == MPI_COMM_NULL ?
      old_to_new.getBase().getMPI() : d_mpi;
   const int rank(mpi.getRank());

   /*
    * Remove relationships with old boxes (because
    * they are going away). These are Boxes mapped by
    * old_to_new.
    *
    * Erase local old Boxes from new_to_anchor.
    *
    * If the old boxes have neighbors in the anchor
    * BoxLevel, some relationships from anchor_eto_old should be
    * erased also.  For each neighbor from a remote anchor box to a
    * local old box, add data to mesg_to_owners saying what
    * Box is disappearing and what anchor Box should
    * no longer reference it.
    */
   for (Connector::ConstNeighborhoodIterator iold = old_to_new.begin();
        iold != old_to_new.end(); ++iold) {

      const BoxId& old_gid_gone = *iold;
      const Box old_box_gone(dim, old_gid_gone);

      if (new_to_anchor->hasNeighborSet(old_gid_gone)) {
         // old_gid_gone exists in new_to_anchor.  Remove it.

         Connector::ConstNeighborhoodIterator affected_anchor_nbrhd =
            new_to_anchor->find(old_gid_gone);

         if (s_print_steps == 'y') {
            tbox::plog << "Box " << old_box_gone
                       << " is gone." << std::endl;
         }

         for (Connector::ConstNeighborIterator ianchor =
                 new_to_anchor->begin(affected_anchor_nbrhd);
              ianchor != new_to_anchor->end(affected_anchor_nbrhd); /* incremented in loop */) {

            if (s_print_steps == 'y') {
               tbox::plog << "  Box " << *ianchor
                          << " is affected." << std::endl;
            }

            const int anchor_nabr_owner = ianchor->getOwnerRank();
            if (anchor_nabr_owner == rank) {
               // Erase local relationship from anchor to old_gid_gone.
               do {

                  if (s_print_steps == 'y') {
                     tbox::plog << "  Fixing affected box " << *ianchor
                                << std::endl;
                  }

                  TBOX_ASSERT(anchor_to_new.hasNeighborSet(ianchor->getBoxId()));

                  if (s_print_steps == 'y') {
                     anchor_to_new.writeNeighborhoodToStream(
                        tbox::plog,
                        ianchor->getBoxId());
                     tbox::plog << std::endl;
                  }
                  if (anchor_to_new.hasLocalNeighbor(ianchor->getBoxId(),
                                                     old_box_gone)) {
                     if (s_print_steps == 'y') {
                        tbox::plog << "    Removing neighbor " << old_box_gone
                                   << " from list for " << *ianchor << std::endl;
                     }
                     anchor_to_new.eraseNeighbor(old_box_gone,
                                                 ianchor->getBoxId());
                  }

                  ++ianchor;

                  // Skip past periodic image Boxes.
                  while (ianchor != new_to_anchor->end(affected_anchor_nbrhd) &&
                         ianchor->isPeriodicImage()) {
                     ++ianchor;
                  }

               } while (ianchor != new_to_anchor->end(affected_anchor_nbrhd) &&
                        ianchor->getOwnerRank() == rank);
            } else {
               // Tell owner of nabr to erase references to old_gid_gone.
               std::vector<int>& mesg = send_mesgs[anchor_nabr_owner];
               // mesg[0] is the counter for how many boxes are removed.
               if (mesg.empty()) {
                  mesg.insert(mesg.end(), 1);
               } else {
                  ++mesg[0];
               }
               mesg.insert(mesg.end(), old_gid_gone.getLocalId().getValue());
               mesg.insert(mesg.end(), -1);
               int i_count = static_cast<int>(mesg.size());
               mesg.insert(mesg.end(), 0);
               do {
                  mesg.insert(mesg.end(), ianchor->getLocalId().getValue());
                  mesg.insert(mesg.end(), 
                     static_cast<int>(ianchor->getBlockId().getBlockValue()));
                  ++mesg[i_count];
                  if (s_print_steps == 'y') tbox::plog
                     << "    Request change " << mesg[i_count]
                     << " to neighbors fo " << *ianchor << std::endl;
                  ++ianchor;
               } while (ianchor != new_to_anchor->end(affected_anchor_nbrhd) &&
                        ianchor->getOwnerRank() == anchor_nabr_owner);
            }
         }

         /*
          * Erase relationships from old_box_gone to anchor
          * box_level.
          */
         new_to_anchor->eraseLocalNeighborhood(old_gid_gone);

      }

   }

   d_object_timers->t_modify_remove_and_cache->stop();
}

/*
 ***********************************************************************
 * Discover overlaps with new Boxes and send outgoing messages.  We thread
 * the discovery process for non-local overlaps and then send the overlap
 * results.  Then we discover the local overlaps and set them in the
 * Connector(s).
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::privateModify_discoverAndSend(
   std::map<int, std::vector<int> >& send_mesgs,
   Connector& anchor_to_new,
   Connector* new_to_anchor,
   const std::set<int>& incoming_ranks,
   const std::set<int>& outgoing_ranks,
   tbox::AsyncCommPeer<int>* all_comms,
   BoxContainer& visible_new_nabrs,
   BoxContainer& visible_anchor_nabrs,
   const InvertedNeighborhoodSet& anchor_eto_old,
   const InvertedNeighborhoodSet& new_eto_old,
   const Connector& old_to_anchor,
   const Connector& anchor_to_old,
   const MappingConnector& old_to_new) const
{
   if (visible_anchor_nabrs.empty() && visible_new_nabrs.empty()) {
      return;
   }

   /*
    * Discover overlaps.  Overlaps are either locally stored or
    * packed into a message for sending.
    */

   d_object_timers->t_modify_discover_and_send->start();

   const BoxLevel& old(old_to_new.getBase());

   const tbox::Dimension& dim(old.getDim());
   const tbox::SAMRAI_MPI& mpi = d_mpi.getCommunicator() == MPI_COMM_NULL ?
      old_to_new.getBase().getMPI() : d_mpi;
   const int rank = mpi.getRank();

   /*
    * Local process can find some neighbors for the (local and
    * remote) Boxes in visible_anchor_nabrs and visible_new_nabrs.
    * Separate this into 2 parts: discovery of remote Boxes which
    * may be threaded, and discovery of local Boxes can not be.
    * In either case we loop through the visible_anchor_nabrs and
    * compare each to visible_new_nabrs, looking for overlaps.
    * Then vice versa.  Since each of these NeighborSets is
    * ordered by processor owner first and we know each non-local
    * processor we can construct each non-local message in a
    * separate thread and then find and set all the local overlaps.
    *
    * To do this we first separate visible_anchor_nabrs into 2 groups
    * non-local and local neighbors.  Also do the same for
    * visible_new_nabrs.
    */
   bool ordered = true;
   BoxContainer visible_local_anchor_nabrs(ordered);
   BoxContainer visible_local_new_nabrs(ordered);
   const Box this_proc_start(dim, GlobalId(LocalId::getZero(), rank));
   BoxContainer::iterator anchor_ni =
      visible_anchor_nabrs.lowerBound(this_proc_start);
   BoxContainer::iterator new_ni =
      visible_new_nabrs.lowerBound(this_proc_start);
   while (anchor_ni != visible_anchor_nabrs.end() &&
          anchor_ni->getOwnerRank() == rank) {
      visible_local_anchor_nabrs.insert(*anchor_ni);
      visible_anchor_nabrs.erase(anchor_ni++);
   }
   while (new_ni != visible_new_nabrs.end() &&
          new_ni->getOwnerRank() == rank) {
      visible_local_new_nabrs.insert(*new_ni);
      visible_new_nabrs.erase(new_ni++);
   }

   // Discover all non-local overlaps.
   int i = 0;
   int imax = static_cast<int>(outgoing_ranks.size());
   std::vector<int> another_outgoing_ranks(outgoing_ranks.size());
   for (std::set<int>::const_iterator outgoing_ranks_itr(outgoing_ranks.begin());
        outgoing_ranks_itr != outgoing_ranks.end(); ++outgoing_ranks_itr) {
      another_outgoing_ranks[i++] = *outgoing_ranks_itr;
   }
#ifdef HAVE_OPENMP
#pragma omp parallel private(i) num_threads(4)
   {
#pragma omp for schedule(dynamic) nowait
#endif
   for (i = 0; i < imax; ++i) {
      BoxId outgoing_proc_start_id(
         LocalId::getZero(),
         another_outgoing_ranks[i]);
      Box outgoing_proc_start(dim, outgoing_proc_start_id);
      BoxContainer::const_iterator thread_anchor_ni =
         visible_anchor_nabrs.lowerBound(outgoing_proc_start);
      BoxContainer::const_iterator thread_new_ni =
         visible_new_nabrs.lowerBound(outgoing_proc_start);
      privateModify_discover(
         send_mesgs[another_outgoing_ranks[i]],
         anchor_to_new,
         new_to_anchor,
         visible_anchor_nabrs,
         visible_new_nabrs,
         thread_anchor_ni,
         thread_new_ni,
         another_outgoing_ranks[i],
         dim,
         rank,
         anchor_eto_old,
         new_eto_old,
         old_to_anchor,
         anchor_to_old,
         old_to_new);
   }
#ifdef HAVE_OPENMP
}
#endif

   /*
    * Send all non-local overlap messages.
    * As an optimization, send to the next higher ranked process first followed
    * by successively higher processes and finally looping around to process
    * 0 through the next lower ranked process.  This spreads out the sends more
    * evenly and prevents everyone from sending to the same processor (like
    * process 0) at the same time.
    */
   int num_outgoing_ranks = static_cast<int>(outgoing_ranks.size());
   int num_incoming_ranks = static_cast<int>(incoming_ranks.size());
   int num_comms = num_outgoing_ranks + num_incoming_ranks;
   std::set<int>::const_iterator outgoing_ranks_itr(
      outgoing_ranks.lower_bound(rank + 1));
   if (outgoing_ranks_itr == outgoing_ranks.end()) {
      outgoing_ranks_itr = outgoing_ranks.begin();
   }
   int comm_offset = num_incoming_ranks;
   for ( ; comm_offset < num_comms; ++comm_offset) {
      if (all_comms[comm_offset].getPeerRank() == *outgoing_ranks_itr) {
         break;
      }
   }
   TBOX_ASSERT(num_outgoing_ranks == 0 || comm_offset < num_comms);
   for (int outgoing_ranks_ctr = 0;
        outgoing_ranks_ctr < num_outgoing_ranks; ++outgoing_ranks_ctr) {
      std::vector<int>& send_mesg = send_mesgs[*outgoing_ranks_itr];
      tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[comm_offset];
      outgoing_comm.beginSend(
         &send_mesg[0],
         send_mesg.size());
      ++comm_offset;
      ++outgoing_ranks_itr;
      TBOX_ASSERT((outgoing_ranks_itr == outgoing_ranks.end()) ==
         (comm_offset == num_comms));
      if (outgoing_ranks_itr == outgoing_ranks.end()) {
         outgoing_ranks_itr = outgoing_ranks.begin();
      }
      if (comm_offset == num_comms) {
         comm_offset = num_incoming_ranks;
      }
      if (s_print_steps == 'y') {
         tbox::plog << "Sent to " << outgoing_comm.getPeerRank() << std::endl;
      }
   }

   // Discover all local overlaps and store them in the Connector(s).
   BoxContainer::const_iterator anchor_local_ni =
      visible_local_anchor_nabrs.lowerBound(this_proc_start);
   BoxContainer::const_iterator new_local_ni =
      visible_local_new_nabrs.lowerBound(this_proc_start);
   privateModify_discover(
      send_mesgs[rank],
      anchor_to_new,
      new_to_anchor,
      visible_local_anchor_nabrs,
      visible_local_new_nabrs,
      anchor_local_ni,
      new_local_ni,
      rank,
      dim,
      rank,
      anchor_eto_old,
      new_eto_old,
      old_to_anchor,
      anchor_to_old,
      old_to_new);

   d_object_timers->t_modify_discover_and_send->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::privateModify_discover(
   std::vector<int>& send_mesg,
   Connector& anchor_to_new,
   Connector* new_to_anchor,
   const BoxContainer& visible_anchor_nabrs,
   const BoxContainer& visible_new_nabrs,
   BoxContainer::const_iterator& anchor_ni,
   BoxContainer::const_iterator& new_ni,
   int curr_owner,
   const tbox::Dimension& dim,
   int rank,
   const InvertedNeighborhoodSet& anchor_eto_old,
   const InvertedNeighborhoodSet& new_eto_old,
   const Connector& old_to_anchor,
   const Connector& anchor_to_old,
   const MappingConnector& old_to_new) const
{
   if (s_print_steps == 'y') {
      tbox::plog << "visible_anchor_nabrs:" << std::endl;
      for (BoxContainer::const_iterator na = visible_anchor_nabrs.begin();
           na != visible_anchor_nabrs.end(); ++na) {
         tbox::plog << "  " << *na << std::endl;
      }
      tbox::plog << "visible_new_nabrs:" << std::endl;
      for (BoxContainer::const_iterator na = visible_new_nabrs.begin();
           na != visible_new_nabrs.end(); ++na) {
         tbox::plog << "  " << *na << std::endl;
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   // Owners that we have sent messages to.  Used for debugging.
   std::set<int> owners_sent_to;
#endif

   TBOX_ASSERT(owners_sent_to.find(curr_owner) == owners_sent_to.end());

   /*
    * Set up send_message to contain info discovered
    * locally and should be sent to curr_owner.
    *
    * Content of send_mesg:
    * - neighbor-removal section cached in send_mesg.
    * - offset to the reference section (see below)
    * - number of anchor boxes for which neighbors are found
    * - number of new boxes for which neighbors are found
    *   - id of base/head box
    *   - number of neighbors found for base/head box.
    *     - BoxId of neighbors found.
    *       Boxes of found neighbors are given in the
    *       reference section of the message.
    * - reference section: all the boxes referenced as
    *   neighbors (accumulated in referenced_anchor_nabrs
    *   and referenced_new_nabrs).
    *   - number of referenced base neighbors
    *   - number of referenced head neighbors
    *   - referenced base neighbors
    *   - referenced head neighbors
    *
    * The purpose of factoring out info on the neighbors referenced
    * is to reduce redundant data that can eat up lots of memory
    * and message passing bandwidth when there are lots of Boxes
    * with the same neighbors.
    */

   /*
    * The first section of the send_mesg is the remote neighbor-removal
    * section (computed above).
    */
   if (curr_owner != rank && send_mesg.empty()) {
      // No neighbor-removal data found for curr_owner.
      send_mesg.insert(send_mesg.end(), 0);
   }

   // Indices of certain positions in send_mesg.
   const int idx_offset_to_ref = static_cast<int>(send_mesg.size());
   const int idx_num_anchor_boxes = idx_offset_to_ref + 1;
   const int idx_num_new_boxes = idx_offset_to_ref + 2;
   send_mesg.insert(send_mesg.end(), 3, 0);

   // Boxes referenced in the message, used when adding ref section.
   BoxContainer referenced_anchor_nabrs;
   BoxContainer referenced_new_nabrs;

   /*
    * Find locally visible new neighbors for all anchor Boxes owned by
    * curr_owner.
    */
   privateModify_findOverlapsForOneProcess(
      curr_owner,
      visible_anchor_nabrs,
      anchor_ni,
      send_mesg,
      idx_num_anchor_boxes,
      anchor_to_new,
      referenced_new_nabrs,
      anchor_to_old,
      old_to_anchor,
      old_to_new,
      anchor_eto_old,
      old_to_new.getHead().getRefinementRatio());

   /*
    * Find locally visible anchor neighbors for all new Boxes owned by
    * curr_owner.
    */
   privateModify_findOverlapsForOneProcess(
      curr_owner,
      visible_new_nabrs,
      new_ni,
      send_mesg,
      idx_num_new_boxes,
      *new_to_anchor,
      referenced_anchor_nabrs,
      *new_to_anchor,
      anchor_to_new,
      old_to_anchor,
      new_eto_old,
      old_to_anchor.getHead().getRefinementRatio());

   if (curr_owner != rank) {
      /*
       * If this discovery is off processor then the send message must be
       * filled with the referenced neighbors.
       */

      packReferencedNeighbors(
         send_mesg,
         idx_offset_to_ref,
         referenced_new_nabrs,
         referenced_anchor_nabrs,
         dim,
         s_print_steps == 'y');

#ifdef DEBUG_CHECK_ASSERTIONS
      owners_sent_to.insert(curr_owner);
#endif

   } // Block to send discoveries to curr_owner.
}

/*
 ***********************************************************************
 *
 * Find overlaps from visible_base_nabrs to head_rbbt.  Find only
 * overlaps for Boxes owned by owner_rank.
 *
 * On entry, base_ni points to the first Box in visible_base_nabrs
 * owned by owner_rank.  Increment base_ni past those Boxes
 * processed and remove them from visible_base_nabrs.
 *
 * Save local and semilocal overlaps in bridging_connector.  For
 * remote overlaps, pack in send_mesg, add head Box to
 * referenced_head_nabrs and increment
 * send_mesg[remote_box_counter_index].
 *
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::privateModify_findOverlapsForOneProcess(
   const int owner_rank,
   const BoxContainer& visible_base_nabrs,
   BoxContainer::const_iterator& base_ni,
   std::vector<int>& send_mesg,
   int remote_box_counter_index,
   Connector& mapped_connector,
   BoxContainer& referenced_head_nabrs,
   const Connector& unmapped_connector,
   const Connector& unmapped_connector_transpose,
   const Connector& mapping_connector,
   const InvertedNeighborhoodSet& inverted_nbrhd,
   const IntVector& head_refinement_ratio) const
{
#ifndef HAVE_OPENMP
   d_object_timers->t_modify_find_overlaps_for_one_process->start();
#endif

   const BoxLevel& old = mapping_connector.getBase();
   const std::shared_ptr<const BaseGridGeometry>& grid_geometry(
      old.getGridGeometry());
   const tbox::SAMRAI_MPI& mpi = d_mpi.getCommunicator() == MPI_COMM_NULL ? old.getMPI() : d_mpi;
   const int rank = mpi.getRank();

   while (base_ni != visible_base_nabrs.end() &&
          base_ni->getOwnerRank() == owner_rank) {
      const Box& base_box = *base_ni;
      if (s_print_steps == 'y') {
         tbox::plog << "Finding neighbors for base_box "
                    << base_box << std::endl;
      }
      Box compare_box = base_box;
      BoxContainer compare_boxes;

      if (grid_geometry->getNumberBlocks() == 1 ||
          grid_geometry->hasIsotropicRatios()) {
         compare_box.grow(mapped_connector.getConnectorWidth());
         if (unmapped_connector.getHeadCoarserFlag()) {
            compare_box.coarsen(unmapped_connector.getRatio());
         }
         else if (unmapped_connector_transpose.getHeadCoarserFlag()) {
            compare_box.refine(unmapped_connector_transpose.getRatio());
         }
         compare_boxes.push_back(compare_box);
      } else {
         TBOX_ASSERT(unmapped_connector.getRatio() ==
                     unmapped_connector_transpose.getRatio());
         BoxUtilities::growAndAdjustAcrossBlockBoundary(
            compare_boxes,
            compare_box,
            grid_geometry,
            mapped_connector.getBase().getRefinementRatio(),
            unmapped_connector.getRatio(),
            mapped_connector.getConnectorWidth(), 
            unmapped_connector_transpose.getHeadCoarserFlag(),
            unmapped_connector.getHeadCoarserFlag());
      }

      std::vector<Box> found_nabrs;
      for (BoxContainer::iterator c_itr = compare_boxes.begin();
           c_itr != compare_boxes.end(); ++c_itr) {
         const Box& comp_box = *c_itr;
         BlockId compare_box_block_id(comp_box.getBlockId());
         Box transformed_compare_box(comp_box);

         InvertedNeighborhoodSet::const_iterator ini =
            inverted_nbrhd.find(base_box);
         if (ini != inverted_nbrhd.end()) {
            const BoxIdSet& old_indices = ini->second;

            for (BoxIdSet::const_iterator na = old_indices.begin();
                 na != old_indices.end(); ++na) {
               Connector::ConstNeighborhoodIterator nbrhd =
                  mapping_connector.findLocal(*na);
               if (nbrhd != mapping_connector.end()) {
                  /*
                   * There are anchor Boxes with relationships to
                   * the old Box identified by *na.
                   */
                  for (Connector::ConstNeighborIterator naa =
                       mapping_connector.begin(nbrhd);
                       naa != mapping_connector.end(nbrhd); ++naa) {
                     const Box& new_nabr(*naa);
                     transformed_compare_box = comp_box;
                     bool do_intersect = true;
                     if (compare_box_block_id != new_nabr.getBlockId()) {
                        // Re-transform compare_box and note its new BlockId.
                        do_intersect = 
                           grid_geometry->transformBox(
                              transformed_compare_box,
                              head_refinement_ratio,
                              new_nabr.getBlockId(),
                              compare_box_block_id);
                     }
                     if (do_intersect) {
                        if (transformed_compare_box.intersects(new_nabr)) {
                           found_nabrs.insert(found_nabrs.end(), *naa);
                        }
                     }
                  }
               }
            }
         }
      }
      if (s_print_steps == 'y') {
         tbox::plog << "Found " << found_nabrs.size() << " neighbors :";
         BoxContainerUtils::recursivePrintBoxVector(
            found_nabrs,
            tbox::plog,
            "\n");
         tbox::plog << std::endl;
      }
      if (!found_nabrs.empty()) {
         if (base_box.getOwnerRank() != rank) {
            // Pack up info for sending.
            ++send_mesg[remote_box_counter_index];
            const int subsize = 3 +
               BoxId::commBufferSize() * static_cast<int>(found_nabrs.size());
            send_mesg.insert(send_mesg.end(), subsize, -1);
            int* submesg = &send_mesg[send_mesg.size() - subsize];
            *(submesg++) = base_box.getLocalId().getValue();
            *(submesg++) = static_cast<int>(
               base_box.getBlockId().getBlockValue());
            *(submesg++) = static_cast<int>(found_nabrs.size());
            for (std::vector<Box>::const_iterator na = found_nabrs.begin();
                 na != found_nabrs.end(); ++na) {
               const Box& nabr = *na;
               referenced_head_nabrs.insert(nabr);
               nabr.getBoxId().putToIntBuffer(submesg);
               submesg += BoxId::commBufferSize();
            }
         } else {
            /*
             * Save neighbor info locally.
             *
             * To improve communication time, we should really send
             * the head neighbors before doing anything locally.
             */
            if (!found_nabrs.empty()) {
               Connector::NeighborhoodIterator base_box_itr =
                  mapped_connector.makeEmptyLocalNeighborhood(base_box.getBoxId());
               for (std::vector<Box>::const_iterator na = found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  mapped_connector.insertLocalNeighbor(*na, base_box_itr);
               }
            }
         }
      }
      if (s_print_steps == 'y') {
         tbox::plog << "Erasing visible base nabr " << (*base_ni) << std::endl;
      }
      ++base_ni;
      if (s_print_steps == 'y') {
         if (base_ni == visible_base_nabrs.end()) {
            tbox::plog << "Next base nabr: end" << std::endl;
         } else {
            tbox::plog << "Next base nabr: " << *base_ni << std::endl;
         }
      }
   }

#ifndef HAVE_OPENMP
   d_object_timers->t_modify_find_overlaps_for_one_process->stop();
#endif
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
MappingConnectorAlgorithm::initializeCallback()
{
   // Initialize timers with default prefix.
   getAllTimers(s_default_timer_prefix,
      s_static_timers[s_default_timer_prefix]);

}

/*
 ***************************************************************************
 * Free statics.  To be called by shutdown registry to make sure
 * memory for statics do not leak.
 ***************************************************************************
 */

void
MappingConnectorAlgorithm::finalizeCallback()
{
   s_static_timers.clear();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::string timer_prefix_used;
   if (s_ignore_external_timer_prefix == 'y') {
      timer_prefix_used = s_default_timer_prefix;
   } else {
      timer_prefix_used = timer_prefix;
   }
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix_used));
   if (ti == s_static_timers.end()) {
      d_object_timers = &s_static_timers[timer_prefix_used];
      getAllTimers(timer_prefix_used, *d_object_timers);
   } else {
      d_object_timers = &(ti->second);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
MappingConnectorAlgorithm::getAllTimers(
   const std::string& timer_prefix,
   TimerStruct& timers)
{
   timers.t_modify_public = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::modify()_public");
   timers.t_modify = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify()");
   timers.t_modify_setup_comm = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::setupCommunication()");
   timers.t_modify_remove_and_cache = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify_removeAndCache()");
   timers.t_modify_discover_and_send = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify_discoverAndSend()");
   timers.t_modify_find_overlaps_for_one_process = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify_findOverlapsForOneProcess()");
   timers.t_modify_receive_and_unpack = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::receiveAndUnpack()");
   timers.t_modify_MPI_wait = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify()_MPI_wait");
   timers.t_modify_misc = tbox::TimerManager::getManager()->
      getTimer(timer_prefix + "::privateModify()_misc");
}

}
}
