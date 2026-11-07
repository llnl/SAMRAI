/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A class describing the adjacency of Boxes.
 *
 ************************************************************************/
#include "SAMRAI/hier/BoxNeighborhoodCollection.h"
#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace hier {

const int BoxNeighborhoodCollection::HIER_BOX_NBRHD_COLLECTION_VERSION = 0;

BoxNeighborhoodCollection::BoxNeighborhoodCollection(
   const BoxContainer& base_boxes)
{
   // For each base Box in base_boxes create an empty neighborhood.
   for (BoxContainer::const_iterator itr = base_boxes.begin();
        itr != base_boxes.end(); ++itr) {
      insert(itr->getBoxId());
   }
}

BoxNeighborhoodCollection::BoxNeighborhoodCollection(
   const BoxNeighborhoodCollection& other)
{
   // Iterate through the other collection and create in this the same
   // neighborhoods that the other contains.
   for (ConstIterator base_boxes_itr(other.begin());
        base_boxes_itr != other.end(); ++base_boxes_itr) {
      Iterator new_base_box = insert(*base_boxes_itr).first;
      for (ConstNeighborIterator nbr_itr(other.begin(base_boxes_itr));
           nbr_itr != other.end(base_boxes_itr); ++nbr_itr) {
         insert(new_base_box, *nbr_itr);
      }
   }
}

BoxNeighborhoodCollection&
BoxNeighborhoodCollection::operator = (
   const BoxNeighborhoodCollection& rhs)
{
   // Empty this container then iterate through the other collection and
   // create in this the same neighborhoods that the other contains.
   clear();
   for (ConstIterator base_boxes_itr(rhs.begin());
        base_boxes_itr != rhs.end(); ++base_boxes_itr) {
      Iterator new_base_box = insert(*base_boxes_itr).first;
      for (ConstNeighborIterator nbr_itr(rhs.begin(base_boxes_itr));
           nbr_itr != rhs.end(base_boxes_itr); ++nbr_itr) {
         insert(new_base_box, *nbr_itr);
      }
   }
   return *this;
}

bool
BoxNeighborhoodCollection::operator == (
   const BoxNeighborhoodCollection& rhs) const
{
   // Both collections must belong to the same process and have the same
   // number of box neighborhoods.
   bool result = numBoxNeighborhoods() == rhs.numBoxNeighborhoods();

   if (result) {
      // Both collections have the same number of box neighborhoods so now
      // check that the neighborhoods are the same.
      ConstIterator base_boxes_itr(begin());
      ConstIterator rhs_base_boxes_itr(rhs.begin());
      for ( ; base_boxes_itr != end(); ++base_boxes_itr, ++rhs_base_boxes_itr) {
         // Check that this base Box is the same in each collection and that
         // the number of neighbors of this base Box is the same in each
         // collection.
         if (*base_boxes_itr != *rhs_base_boxes_itr) {
            result = false;
            break;
         }
         if (numNeighbors(base_boxes_itr) !=
             rhs.numNeighbors(rhs_base_boxes_itr)) {
            result = false;
            break;
         }

         // This base Box and its number of neighbors match so check that each
         // neighbor from each collection matches.
         ConstNeighborIterator nbr_itr(begin(base_boxes_itr));
         ConstNeighborIterator rhs_nbr_itr(rhs.begin(rhs_base_boxes_itr));
         for ( ; nbr_itr != end(base_boxes_itr); ++nbr_itr, ++rhs_nbr_itr) {
            if (!nbr_itr->isSpatiallyEqual(*rhs_nbr_itr) ||
                !nbr_itr->isIdEqual(*rhs_nbr_itr)) {
               result = false;
               break;
            }
         }

         // If one of the neighbors of this base Box didn't match bail out.
         if (!result) {
            break;
         }
      }
   }

   return result;
}

bool
BoxNeighborhoodCollection::operator != (
   const BoxNeighborhoodCollection& rhs) const
{
   return !(*this == rhs);
}

int
BoxNeighborhoodCollection::sumNumNeighbors() const
{
   // Count the neighbors in each base Box.
   int ct = 0;
   for (ConstIterator base_boxes_itr(begin());
        base_boxes_itr != end(); ++base_boxes_itr) {
      ct += numNeighbors(base_boxes_itr);
   }
   return ct;
}

bool
BoxNeighborhoodCollection::hasNeighbor(
   const ConstIterator& base_box_itr,
   const Box& nbr) const
{
   if (base_box_itr == end()) {
      return false;
   } else {
      HeadBoxPool::const_iterator nbrs_itr = d_nbrs.find(nbr);
      if (nbrs_itr == d_nbrs.end()) {
         return false;
      } else {
         return base_box_itr.d_itr->second.find(&(*nbrs_itr)) !=
                base_box_itr.d_itr->second.end();
      }
   }
}

bool
BoxNeighborhoodCollection::neighborhoodEqual(
   const BoxId& base_box_id,
   const BoxNeighborhoodCollection& other) const
{
   if (numNeighbors(base_box_id) != other.numNeighbors(base_box_id)) {
      return false;
   }
   bool result = true;
   ConstNeighborIterator this_ni = begin(base_box_id);
   ConstNeighborIterator other_ni = other.begin(base_box_id);
   for ( ; this_ni != end(base_box_id); ++this_ni, ++other_ni) {
      if (!this_ni->isIdEqual(*other_ni)) {
         result = false;
         break;
      }
   }
   return result;
}

bool
BoxNeighborhoodCollection::isLocal(
   int rank) const
{
   // See if all the neighbors of all base Boxes belong to the same processor
   // as this object.
   bool is_local = true;
   for (ConstIterator base_box_itr(begin()); is_local && base_box_itr != end();
        ++base_box_itr) {
      for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
           nbr_itr != end(base_box_itr); ++nbr_itr) {
         if (nbr_itr->getOwnerRank() != rank) {
            is_local = false;
            break;
         }
      }
   }
   return is_local;
}

void
BoxNeighborhoodCollection::getOwners(
   std::set<int>& owners) const
{
   // Put the processor id of all neigbhors of all base Boxes into owners.
   for (ConstIterator base_box_itr(begin());
        base_box_itr != end(); ++base_box_itr) {
      getOwners(base_box_itr, owners);
   }
}

void
BoxNeighborhoodCollection::getOwners(
   ConstIterator& base_box_itr,
   std::set<int>& owners) const
{
   // Put the processor id of all neigbhors of the base Box pointed to by
   // base_box_itr into owners.
   for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
        nbr_itr != end(base_box_itr); ++nbr_itr) {
      owners.insert(nbr_itr->getOwnerRank());
   }
}

void
BoxNeighborhoodCollection::insert(
   Iterator& base_box_itr,
   const Box& new_nbr)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   // First add the new_nbr to the collection of neighbors if it is not there.
   HeadBoxPool::iterator nbr_itr = d_nbrs.find(new_nbr);
   if (nbr_itr == d_nbrs.end()) {
      nbr_itr = d_nbrs.insert(new_nbr).first;
      const Box& tmp = *nbr_itr;
      d_nbr_link_ct[&tmp] = 0;
   }
   const Box& nbr_in_d_nbrs = *nbr_itr;

   // Now add the neighbor to the base Box's neighborhood.
   bool new_nbr_ref = base_box_itr.d_itr->second.insert(&nbr_in_d_nbrs).second;

   if (new_nbr_ref) {
      ++(d_nbr_link_ct.find(&nbr_in_d_nbrs)->second);
   }
}

void
BoxNeighborhoodCollection::insert(
   Iterator& base_box_itr,
   const BoxContainer& new_nbrs)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   // Add each neighbor in the container to the base Box.
   for (BoxContainer::const_iterator new_nbr_itr = new_nbrs.begin();
        new_nbr_itr != new_nbrs.end(); ++new_nbr_itr) {

      // First add this new neighbor to the collection of neighbors if it is
      // not there.
      const Box& new_nbr = *new_nbr_itr;
      HeadBoxPool::iterator nbr_itr = d_nbrs.find(new_nbr);
      if (nbr_itr == d_nbrs.end()) {
         nbr_itr = d_nbrs.insert(new_nbr).first;
         const Box& tmp = *nbr_itr;
         d_nbr_link_ct[&tmp] = 0;
      }
      const Box& nbr_in_d_nbrs = *nbr_itr;

      // Now add this new neighbor to the base Box's neighborhood.
      bool new_nbr_ref =
         base_box_itr.d_itr->second.insert(&nbr_in_d_nbrs).second;

      if (new_nbr_ref) {
         ++(d_nbr_link_ct.find(&nbr_in_d_nbrs)->second);
      }
   }
}

void
BoxNeighborhoodCollection::erase(
   Iterator& base_box_itr,
   const Box& nbr)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   HeadBoxPool::iterator nbr_itr = d_nbrs.find(nbr);
   TBOX_ASSERT(nbr_itr != d_nbrs.end());

   const Box& nbr_in_d_nbrs = *nbr_itr;
   Neighborhood::size_type nbr_erased =
      base_box_itr.d_itr->second.erase(&nbr_in_d_nbrs);
   if (nbr_erased != 0) {
      int& link_ct = d_nbr_link_ct[&nbr_in_d_nbrs];
      if (link_ct == 1) {
         d_nbr_link_ct.erase(&nbr_in_d_nbrs);
         d_nbrs.erase(nbr_itr);
      } else {
         --link_ct;
      }
   }
}

void
BoxNeighborhoodCollection::erase(
   Iterator& base_box_itr,
   const BoxContainer& nbrs)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   // Remove each neighbor in the container from the base Box.
   for (BoxContainer::const_iterator old_nbr_itr = nbrs.begin();
        old_nbr_itr != nbrs.end(); ++old_nbr_itr) {

      HeadBoxPool::iterator nbr_itr = d_nbrs.find(*old_nbr_itr);
      TBOX_ASSERT(nbr_itr != d_nbrs.end());

      const Box& nbr_in_d_nbrs = *nbr_itr;
      Neighborhood::size_type nbr_erased =
         base_box_itr.d_itr->second.erase(&nbr_in_d_nbrs);
      if (nbr_erased != 0) {
         int& link_ct = d_nbr_link_ct[&nbr_in_d_nbrs];
         if (link_ct == 1) {
            d_nbr_link_ct.erase(&nbr_in_d_nbrs);
            d_nbrs.erase(nbr_itr);
         } else {
            --link_ct;
         }
      }
   }
}

BoxNeighborhoodCollection::InsertRetType
BoxNeighborhoodCollection::insert(
   const BoxId& new_base_box)
{
   // First, add the base Box to the pool of base Boxes.  If it's already there
   // this is a no-op.
   std::pair<BaseBoxPoolItr, bool> base_box_insert_info =
      d_base_boxes.insert(new_base_box);

   // Second, add the empty neighborhood for the base Box to the adjacency list
   // if it is not there.
   if (base_box_insert_info.second == true) {
      AdjListItr base_box_itr = d_adj_list.insert(
            d_adj_list.end(),
            std::make_pair(&(*(base_box_insert_info.first)), Neighborhood()));
      return std::make_pair(Iterator(*this, base_box_itr), true);
   } else {
      AdjListItr base_box_itr(d_adj_list.find(&(*(base_box_insert_info.first))));
      return std::make_pair(Iterator(*this, base_box_itr), false);
   }
}

void
BoxNeighborhoodCollection::erase(
   Iterator& base_box_itr)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   // Erasing base Boxes so clobber entire d_adj_list entry and d_base_boxes
   // entry.
   d_base_boxes.erase(base_box_itr.d_base_boxes_itr);
   d_adj_list.erase(base_box_itr.d_itr);
}

void
BoxNeighborhoodCollection::erase(
   Iterator& first_base_box_itr,
   Iterator& last_base_box_itr)
{
   // For each base Box in the range erase it.
   for (Iterator base_box_itr(first_base_box_itr);
        base_box_itr != last_base_box_itr; ) {
      Iterator current_base_box(base_box_itr);
      ++base_box_itr;
      erase(current_base_box);
   }
}

void
BoxNeighborhoodCollection::eraseNonLocalNeighborhoods(
   int rank)
{
   // Find all base Boxes which do not belong to the same processor as this
   // object and remove them entirely.
   for (Iterator base_box_itr(begin()); base_box_itr != end(); ) {
      if (base_box_itr->getOwnerRank() != rank) {
         Iterator current_base_box(base_box_itr);
         ++base_box_itr;
         erase(current_base_box);
      } else {
         ++base_box_itr;
      }
   }
}

void
BoxNeighborhoodCollection::eraseEmptyNeighborhoods()
{
   // Find all base Boxes which have no neighbors and remove them entirely.
   for (Iterator base_box_itr(begin()); base_box_itr != end(); ) {
      if (base_box_itr.d_itr->second.empty()) {
         Iterator current_base_box(base_box_itr);
         ++base_box_itr;
         erase(current_base_box);
      } else {
         ++base_box_itr;
      }
   }
}

void
BoxNeighborhoodCollection::erasePeriodicNeighbors()
{
   for (Iterator base_box_itr(begin()); base_box_itr != end(); ++base_box_itr) {
      for (NeighborIterator nbr(begin(base_box_itr)); nbr != end(base_box_itr); ) {
         const Box& nbr_box = *nbr;
         ++nbr;
         if (nbr_box.isPeriodicImage()) {
            erase(base_box_itr, nbr_box);
         }
      }
   }
}

void
BoxNeighborhoodCollection::clear()
{
   // Entirely remove all base Boxes and their neighborhoods.
   d_adj_list.clear();
   d_base_boxes.clear();
   d_nbr_link_ct.clear();
   d_nbrs.clear();
}

void
BoxNeighborhoodCollection::coarsenNeighbors(
   const IntVector& ratio)
{
   for (HeadBoxPool::iterator nbr_itr(d_nbrs.begin());
        nbr_itr != d_nbrs.end(); ++nbr_itr) {
      Box& box_to_coarsen = const_cast<Box&>(*nbr_itr);
      box_to_coarsen.coarsen(ratio);
   }
}

void
BoxNeighborhoodCollection::refineNeighbors(
   const IntVector& ratio)
{
   for (HeadBoxPool::iterator nbr_itr(d_nbrs.begin());
        nbr_itr != d_nbrs.end(); ++nbr_itr) {
      Box& box_to_refine = const_cast<Box&>(*nbr_itr);
      box_to_refine.refine(ratio);
   }
}

void
BoxNeighborhoodCollection::growNeighbors(
   const IntVector& growth)
{
   for (HeadBoxPool::iterator nbr_itr(d_nbrs.begin());
        nbr_itr != d_nbrs.end(); ++nbr_itr) {
      Box& box_to_grow = const_cast<Box&>(*nbr_itr);
      box_to_grow.grow(growth);
   }
}

void
BoxNeighborhoodCollection::getNeighbors(
   BoxContainer& neighbors) const
{
   // Iterate through the neighbors of each neighborhood and dump them into the
   // BoxContainer.
   if (neighbors.isOrdered()) {
      for (ConstIterator base_box_itr(begin());
           base_box_itr != end(); ++base_box_itr) {
         for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
              nbr_itr != end(base_box_itr); ++nbr_itr) {
            neighbors.insert(*nbr_itr);
         }
      }
   } else {
      for (ConstIterator base_box_itr(begin());
           base_box_itr != end(); ++base_box_itr) {
         for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
              nbr_itr != end(base_box_itr); ++nbr_itr) {
            neighbors.pushBack(*nbr_itr);
         }
      }
   }
}

void
BoxNeighborhoodCollection::getNeighbors(
   BoxContainer& neighbors,
   const BlockId& block_id) const
{
   // Iterate through the neighbors of each neighborhood and dump those with
   // the requested BlockId into the BoxContainer.
   for (ConstIterator base_box_itr(begin());
        base_box_itr != end(); ++base_box_itr) {
      for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
           nbr_itr != end(base_box_itr); ++nbr_itr) {
         const Box& nbr = *nbr_itr;
         if (nbr.getBlockId() == block_id) {
            neighbors.insert(nbr);
         }
      }
   }
}

void
BoxNeighborhoodCollection::getNeighbors(
   std::map<BlockId, BoxContainer>& neighbors) const
{
   for (ConstIterator base_box_itr(begin());
        base_box_itr != end(); ++base_box_itr) {
      for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
           nbr_itr != end(base_box_itr); ++nbr_itr) {
         const Box& nbr = *nbr_itr;
         neighbors[nbr.getBlockId()].insert(nbr);
      }
   }
}

void
BoxNeighborhoodCollection::getNeighbors(
   const BoxId& base_box_id,
   BoxContainer& neighbors) const
{
   ConstIterator nbrhd = find(base_box_id);
   for (ConstNeighborIterator ni = begin(nbrhd); ni != end(nbrhd); ++ni) {
      neighbors.insert(*ni);
   }
}

void
BoxNeighborhoodCollection::getPeriodicNeighbors(
   BoxContainer& result) const
{
   // Iterate through each neighbor of each base Box and place each neighbor
   // which is a periodic image into the supplied BoxContainer.
   for (ConstIterator base_box_itr(begin());
        base_box_itr != end(); ++base_box_itr) {
      for (ConstNeighborIterator nbr_itr(begin(base_box_itr));
           nbr_itr != end(base_box_itr); ++nbr_itr) {
         const Box& nbr = *nbr_itr;
         if (nbr.isPeriodicImage()) {
            result.insert(nbr);
         }
      }
   }
}

void
BoxNeighborhoodCollection::putToIntBuffer(
   std::vector<int>& send_mesg,
   const tbox::Dimension& dim,
   int buff_init) const
{
   /*
    * Pack neighborhood info into send_mesg.
    */
   /*
    * Information to be packed:
    *   - Number of local boxes with neighbors (no info to send
    *     for those boxes without neighbors)
    *   - For each local box,
    *     - box's local index
    *     - number of neighbors
    *     - neighbors
    */
   int box_id_com_buf_size = BoxId::commBufferSize();
   int box_com_buf_size = Box::commBufferSize(dim);
   const int num_box_neighborhoods = numBoxNeighborhoods();
   int num_nabrs = sumNumNeighbors();

   /* Message size is 1 for number of local base Boxes with neighbor lists plus
    * the number of neighbors for each base Box plus the amount of box data
    * written for each base Box and each base Box's neighbors. */
   const int mesg_size = 1 + num_box_neighborhoods
      + (num_box_neighborhoods * box_id_com_buf_size)
      + (num_nabrs * box_com_buf_size);

   send_mesg.resize(mesg_size, buff_init);
   send_mesg[0] = num_box_neighborhoods;
   int imesg = 1;

   for (ConstIterator ci(begin()); ci != end(); ++ci) {

      /* Info for each base Box. */
      ci->putToIntBuffer(&send_mesg[imesg]);
      imesg += box_id_com_buf_size;

      /* Number of neighbors for this base Box. */
      send_mesg[imesg++] = numNeighbors(ci);

      for (ConstNeighborIterator ni(begin(ci)); ni != end(ci); ++ni) {
         /* Info for each neighbor. */
         ni->putToIntBuffer(&send_mesg[imesg]);
         imesg += box_com_buf_size;
      }

   }

   TBOX_ASSERT(imesg == mesg_size);
}

void
BoxNeighborhoodCollection::getFromIntBuffer(
   const std::vector<int>& recv_mesg,
   const std::vector<int>& proc_offset,
   const tbox::Dimension& dim,
   int num_proc,
   int rank)
{
   /*
    * Unpack neighborhood info from recv_mesg.
    */
   int box_id_com_buf_size = BoxId::commBufferSize();
   int box_com_buf_size = Box::commBufferSize(dim);

   BoxId base_box;
   Box nabr(dim);
   for (int n = 0; n < num_proc; ++n) {

      if (n != rank) {

         const int* ptr = &recv_mesg[0] + proc_offset[n];
         const int num_nabr_lists = *(ptr++);

         for (int i = 0; i < num_nabr_lists; ++i) {

            base_box.getFromIntBuffer(ptr);
            ptr += box_id_com_buf_size;

            const int num_nabrs = (*ptr++);

            Iterator base_box_loc = insert(base_box).first;

            for (int nn = 0; nn < num_nabrs; ++nn) {
               nabr.getFromIntBuffer(ptr);
               ptr += box_com_buf_size;
               insert(base_box_loc, nabr);
            }

         }

      }
   }
}

void
BoxNeighborhoodCollection::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   // This appears to be used in the RedistributedRestartUtility.
   restart_db->putBool("d_is_edge_set", true);

   restart_db->putInteger(
      "HIER_BOX_NBRHD_COLLECTION_VERSION",
      HIER_BOX_NBRHD_COLLECTION_VERSION);
   const int num_neighborhoods = numBoxNeighborhoods();
   restart_db->putInteger("number_of_sets", num_neighborhoods);

   if (num_neighborhoods > 0) {

      std::vector<int> owners(num_neighborhoods);
      std::vector<int> local_indices(num_neighborhoods);
      std::vector<int> periodic_ids(num_neighborhoods);
      int counter = 0;
      for (ConstIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& base_box_id = *ei;
         owners[counter] = base_box_id.getOwnerRank();
         local_indices[counter] = base_box_id.getLocalId().getValue();
         periodic_ids[counter] = base_box_id.getPeriodicId().getPeriodicValue();
         ++counter;
      }

      restart_db->putIntegerVector("owners", owners);
      restart_db->putIntegerVector("local_indices", local_indices);
      restart_db->putIntegerVector("periodic_ids", periodic_ids);
      const std::string set_db_string("set_for_local_id_");

      for (ConstIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& mbid = *ei;
         const std::string set_name =
            set_db_string
            + tbox::Utilities::processorToString(mbid.getOwnerRank())
            + tbox::Utilities::patchToString(mbid.getLocalId().getValue())
            + tbox::Utilities::intToString(mbid.getPeriodicId().getPeriodicValue());
         std::shared_ptr<tbox::Database> nbr_db =
            restart_db->putDatabase(set_name);

         const int mbs_size = numNeighbors(ei);
         nbr_db->putInteger("mapped_box_set_size", mbs_size);

         if (mbs_size > 0) {

            std::vector<int> local_ids(mbs_size);
            std::vector<int> ranks(mbs_size);
            std::vector<int> block_ids(mbs_size);
            std::vector<int> periodic_ids(mbs_size);

            std::vector<tbox::DatabaseBox> db_box_array(mbs_size);

            counter = 0;

            for (ConstNeighborIterator ni = begin(ei); ni != end(ei); ++ni) {
               const Box& nbr = *ni;
               local_ids[counter] = nbr.getLocalId().getValue();
               ranks[counter] = nbr.getOwnerRank();
               block_ids[counter] =
                  static_cast<int>(nbr.getBlockId().getBlockValue());
               periodic_ids[counter] = nbr.getPeriodicId().getPeriodicValue();
               db_box_array[counter] = nbr;
               ++counter;
            }

            nbr_db->putIntegerVector("local_indices", local_ids);
            nbr_db->putIntegerVector("ranks", ranks);
            nbr_db->putIntegerVector("block_ids", block_ids);
            nbr_db->putIntegerVector("periodic_ids", periodic_ids);
            nbr_db->putDatabaseBoxVector("boxes", db_box_array);
         }
      }
   }
}

void
BoxNeighborhoodCollection::getFromRestart(
   tbox::Database& restart_db)
{
   int version = restart_db.getInteger("HIER_BOX_NBRHD_COLLECTION_VERSION");
   if (version != HIER_BOX_NBRHD_COLLECTION_VERSION) {
      TBOX_ERROR("BoxNeighborhoodCollection::getFromRestart() error...\n"
         << "   Restart file version different than class version.");
   }

   const unsigned int number_of_sets = restart_db.getInteger("number_of_sets");
   if (number_of_sets > 0) {

      std::vector<int> owners =
         restart_db.getIntegerVector("owners");
      std::vector<int> local_indices =
         restart_db.getIntegerVector("local_indices");
      std::vector<int> periodic_ids =
         restart_db.getIntegerVector("periodic_ids");

      const std::string set_db_string("set_for_local_id_");

      for (size_t i = 0; i < number_of_sets; ++i) {
         BoxId box_id(
            LocalId(local_indices[i]),
            owners[i],
            PeriodicId(periodic_ids[i]));
         const std::string set_name =
            set_db_string
            + tbox::Utilities::processorToString(box_id.getOwnerRank())
            + tbox::Utilities::patchToString(box_id.getLocalId().getValue())
            + tbox::Utilities::intToString(box_id.getPeriodicId().getPeriodicValue());
         std::shared_ptr<tbox::Database> nbr_db(
            restart_db.getDatabase(set_name));
         const unsigned int mbs_size =
            nbr_db->getInteger("mapped_box_set_size");
         Iterator base_box_loc = insert(box_id).first;
         if (mbs_size > 0) {

            std::vector<int> local_ids =
               nbr_db->getIntegerVector("local_indices");
            std::vector<int> ranks =
               nbr_db->getIntegerVector("ranks");
            std::vector<int> block_ids =
               nbr_db->getIntegerVector("block_ids");
            std::vector<int> periodic_ids =
               nbr_db->getIntegerVector("periodic_ids");
            std::vector<tbox::DatabaseBox> db_box_array =
               nbr_db->getDatabaseBoxVector("boxes");

            for (unsigned int i = 0; i < mbs_size; ++i) {
               Box nbr(db_box_array[i]);
               nbr.setBlockId(BlockId(block_ids[i]));
               BoxId box_id(LocalId(local_ids[i]),
                            ranks[i],
                            PeriodicId(periodic_ids[i]));
               nbr.setId(box_id);
               insert(base_box_loc, nbr);
            }
         }
      }
   }
}

BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   bool from_start):
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
         nbrhds.d_adj_list.end()),
   d_base_boxes_itr(from_start ? nbrhds.d_base_boxes.begin() :
                    nbrhds.d_base_boxes.end())
{
}

BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   AdjListItr itr):
   d_collection(&nbrhds),
   d_itr(itr),
   d_base_boxes_itr(nbrhds.d_base_boxes.find(*(itr->first)))
{
}

BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   bool from_start):
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
         nbrhds.d_adj_list.end()),
   d_base_boxes_itr(from_start ? nbrhds.d_base_boxes.begin() :
                    nbrhds.d_base_boxes.end())
{
}

BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   AdjListConstItr itr):
   d_collection(&nbrhds),
   d_itr(itr),
   d_base_boxes_itr(nbrhds.d_base_boxes.find(*(itr->first)))
{
}

BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const Iterator& other):
   d_collection(other.d_collection),
   d_itr(other.d_itr),
   d_base_boxes_itr(other.d_base_boxes_itr)
{
}

BoxNeighborhoodCollection::NeighborIterator::NeighborIterator(
   Iterator& base_box_itr,
   bool from_start):
   d_collection(base_box_itr.d_collection),
   d_base_box(base_box_itr.d_itr->first),
   d_itr(from_start ? base_box_itr.d_itr->second.begin() :
         base_box_itr.d_itr->second.end())
{
}

BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const ConstIterator& base_box_itr,
   bool from_start):
   d_collection(base_box_itr.d_collection),
   d_base_box(base_box_itr.d_itr->first),
   d_itr(from_start ? base_box_itr.d_itr->second.begin() :
         base_box_itr.d_itr->second.end())
{
}

BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const NeighborIterator& other):
   d_collection(other.d_collection),
   d_base_box(other.d_base_box),
   d_itr(other.d_itr)
{
}

}
}
