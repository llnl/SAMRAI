/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/
#include "SAMRAI/hier/UncoveredBoxIterator.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/FlattenedHierarchy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

UncoveredBoxIterator::UncoveredBoxIterator(
   const PatchHierarchy* hierarchy,
   bool begin):
   d_hierarchy(hierarchy),
   d_uncovered_boxes_itr(BoxContainer().begin()),
   d_uncovered_boxes_itr_end(BoxContainer().end()),
   d_item(nullptr)
{
   TBOX_ASSERT(hierarchy);

   d_finest_level_num = d_hierarchy->getFinestLevelNumber();
   if (begin) {
      d_level_num = -1;
      d_flattened_hierarchy = new FlattenedHierarchy(*d_hierarchy, 0, d_finest_level_num);
      d_allocated_flattened_hierarchy = true;
      findFirstUncoveredBox();
   } else {
      d_level_num = d_finest_level_num + 1;
      d_flattened_hierarchy = nullptr;
      d_allocated_flattened_hierarchy = false;
   }
}


UncoveredBoxIterator::UncoveredBoxIterator(
   const FlattenedHierarchy* flattened_hierarchy,
   bool begin):
   d_hierarchy(&(flattened_hierarchy->getPatchHierarchy())),
   d_uncovered_boxes_itr(BoxContainer().begin()),
   d_uncovered_boxes_itr_end(BoxContainer().end()),
   d_item(nullptr)
{
   TBOX_ASSERT(flattened_hierarchy);

   d_finest_level_num = d_hierarchy->getFinestLevelNumber();
   if (begin) {
      d_level_num = -1;
      d_flattened_hierarchy = flattened_hierarchy;
      d_allocated_flattened_hierarchy = false;
      findFirstUncoveredBox();
   } else {
      d_level_num = d_finest_level_num + 1;
      d_flattened_hierarchy = nullptr;
      d_allocated_flattened_hierarchy = false;
   }
}

UncoveredBoxIterator::UncoveredBoxIterator(
   const UncoveredBoxIterator& other):
   d_hierarchy(other.d_hierarchy),
   d_flattened_hierarchy(nullptr),
   d_allocated_flattened_hierarchy(false),
   d_level_num(other.d_level_num),
   d_current_patch_id(other.d_current_patch_id),
   d_uncovered_boxes_itr(other.d_uncovered_boxes_itr),
   d_uncovered_boxes_itr_end(other.d_uncovered_boxes_itr_end),
   d_item(nullptr),
   d_finest_level_num(other.d_finest_level_num)
{
   if (other.d_item) {
      d_item = new std::pair<std::shared_ptr<Patch>, Box>(*other.d_item);
   }
   if (other.d_flattened_hierarchy) {
      d_flattened_hierarchy =
         new FlattenedHierarchy(*other.d_flattened_hierarchy);
      d_allocated_flattened_hierarchy = true;
      if (d_level_num <= d_finest_level_num) {
         TBOX_ASSERT(d_item);
         const Box& patch_box = d_item->first->getBox();
         const BoxContainer& visible_boxes =
            d_flattened_hierarchy->getVisibleBoxes(patch_box, d_level_num);

         BoxContainer::const_iterator itr = visible_boxes.begin(); 
            for( ; itr != visible_boxes.end(); ++itr) {

            if (itr->getBoxId() == d_item->second.getBoxId() &&
                itr->isSpatiallyEqual(d_item->second)) {

               d_uncovered_boxes_itr = itr;
               d_uncovered_boxes_itr_end = visible_boxes.end();
               break;
            }
         }

         if (itr == visible_boxes.end()) {
            d_uncovered_boxes_itr = itr;
            d_uncovered_boxes_itr_end = itr;
         }
      }
   }
}

UncoveredBoxIterator::~UncoveredBoxIterator()
{
   if (d_item) {
      delete d_item;
   }
   if (d_flattened_hierarchy && d_allocated_flattened_hierarchy) {
      delete d_flattened_hierarchy;
   } 
}

UncoveredBoxIterator&
UncoveredBoxIterator::operator = (
   const UncoveredBoxIterator& rhs)
{
   if (this != &rhs) {
      d_hierarchy = rhs.d_hierarchy;
      d_level_num = rhs.d_level_num;
      d_uncovered_boxes_itr = rhs.d_uncovered_boxes_itr;
      d_uncovered_boxes_itr_end = rhs.d_uncovered_boxes_itr_end;
      d_current_patch_id = rhs.d_current_patch_id; 
      if (d_item) {
         delete d_item;
      }
      if (rhs.d_item) {
         d_item = new std::pair<std::shared_ptr<Patch>, Box>(*rhs.d_item);
         d_item->first = rhs.d_item->first;
         d_item->second = rhs.d_item->second;
      } else {
         d_item = nullptr;
      } 
      d_finest_level_num = rhs.d_finest_level_num;
      if (d_flattened_hierarchy && d_allocated_flattened_hierarchy) {
         delete d_flattened_hierarchy;
      }
      d_flattened_hierarchy = nullptr;
      d_allocated_flattened_hierarchy = false;
      if (rhs.d_flattened_hierarchy) {
         d_flattened_hierarchy =
            new FlattenedHierarchy(*rhs.d_flattened_hierarchy);
         d_allocated_flattened_hierarchy = true;
         if (d_level_num <= d_finest_level_num) {
            TBOX_ASSERT(d_item);
            const Box& patch_box = d_item->first->getBox();
            const BoxContainer& visible_boxes =
            d_flattened_hierarchy->getVisibleBoxes(patch_box, d_level_num);

            BoxContainer::const_iterator itr = visible_boxes.begin();
               for( ; itr != visible_boxes.end(); ++itr) {

               if (itr->getBoxId() == d_item->second.getBoxId() &&
                   itr->isSpatiallyEqual(d_item->second)) {

                  d_uncovered_boxes_itr = itr;
                  d_uncovered_boxes_itr_end = visible_boxes.end();
                  break;
               }
            }

            if (itr == visible_boxes.end()) {
               d_uncovered_boxes_itr = itr;
               d_uncovered_boxes_itr_end = itr;
            }
         }
      }
   }

   return *this;
}

const std::pair<std::shared_ptr<Patch>, Box>&
UncoveredBoxIterator::operator * () const
{
   return *d_item;
}

const std::pair<std::shared_ptr<Patch>, Box> *
UncoveredBoxIterator::operator -> () const
{
   return d_item;
}

bool
UncoveredBoxIterator::operator == (
   const UncoveredBoxIterator& rhs) const
{
   // Frist check and see if the iterators are working on the same hierarchies
   // and levels.  If not then they are not equal.
   bool result = d_hierarchy == rhs.d_hierarchy &&
      d_level_num == rhs.d_level_num;

   if (d_flattened_hierarchy == nullptr && rhs.d_flattened_hierarchy != nullptr) {
      result = false;
   }
   if (d_flattened_hierarchy != nullptr && rhs.d_flattened_hierarchy == nullptr) {
      result = false;
   }
   if (d_item == nullptr && rhs.d_item != nullptr) {
      result = false;
   }
   if (d_item != nullptr && rhs.d_item == nullptr) {
      result = false;
   }
   if (result) {
      if (d_item == nullptr && rhs.d_item == nullptr) {
         result = true;
      }
      if (d_item && rhs.d_item) {
         if (d_item->second.isIdEqual(rhs.d_item->second) &&
             d_item->second.isSpatiallyEqual(rhs.d_item->second) &&
             d_item->first->getBox().isIdEqual(rhs.d_item->first->getBox()) &&
             d_item->first->getBox().isSpatiallyEqual(rhs.d_item->first->getBox())) {
            result = true;
         } else {
            result = false;
         }
      }
   }

   return result;
}

bool
UncoveredBoxIterator::operator != (
   const UncoveredBoxIterator& rhs) const
{
   return !(*this == rhs);
}

UncoveredBoxIterator&
UncoveredBoxIterator::operator ++ ()
{
   incrementIterator();
   return *this;
}

UncoveredBoxIterator
UncoveredBoxIterator::operator ++ (
   int)
{
   // Save the state of the iterator.
   UncoveredBoxIterator tmp(*this);

   incrementIterator();

   // Return iterator in original state.
   return tmp;
}

void
UncoveredBoxIterator::incrementIterator()
{
   /*
    * Increment the iterator over the uncovered boxes for the current patch.
    * If we reach the end of the uncovered boxes for the current patch,
    * look for the next patch with uncovered boxes, moving to finer levels if
    * necessary
    */
   ++d_uncovered_boxes_itr;
   if (d_uncovered_boxes_itr != d_uncovered_boxes_itr_end) {
      d_current_patch_id = d_item->first->getBox().getBoxId();
      setIteratorItem();
   } else {

      bool id_found = false;

      bool new_level = false;
      while (d_level_num <= d_finest_level_num) {

         std::shared_ptr<PatchLevel> this_level =
            d_hierarchy->getPatchLevel(d_level_num);

         const BoxContainer& this_level_boxes =
            this_level->getBoxLevel()->getBoxes();

         for (RealBoxConstIterator this_itr = this_level_boxes.realBegin();
              this_itr != this_level_boxes.realEnd(); ++this_itr) {

            if (!new_level &&
                this_itr->getBoxId() <= d_item->first->getBox().getBoxId()) {
               continue;
            }
            const BoxContainer& uncovered_boxes =
               d_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);

            if (!uncovered_boxes.empty()) {
               d_uncovered_boxes_itr = uncovered_boxes.begin();
               d_uncovered_boxes_itr_end = uncovered_boxes.end();
               d_current_patch_id = this_itr->getBoxId();
               id_found = true;
               break;
            }
         }
         if (id_found) {
            break;
         } else {
            ++d_level_num;
            new_level = true;
         }
      }

      if (id_found) {
         setIteratorItem();
      } else {
         if (d_flattened_hierarchy && d_allocated_flattened_hierarchy) {
            delete d_flattened_hierarchy;
         }
         d_flattened_hierarchy = nullptr;
         if (d_item) {
            delete d_item;
            d_item = nullptr; 
         }
      }
   }
}

void
UncoveredBoxIterator::findFirstUncoveredBox()
{
   ++d_level_num;
   std::shared_ptr<PatchLevel> this_level =
      d_hierarchy->getPatchLevel(d_level_num);

   bool id_found = false;

   while (d_level_num <= d_finest_level_num) { 

      std::shared_ptr<PatchLevel> this_level =
         d_hierarchy->getPatchLevel(d_level_num);

      const BoxContainer& this_level_boxes =
         this_level->getBoxLevel()->getBoxes();

      for (RealBoxConstIterator this_itr = this_level_boxes.realBegin();
           this_itr != this_level_boxes.realEnd(); ++this_itr) {
         const BoxContainer& uncovered_boxes =
            d_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);

         if (!uncovered_boxes.empty()) {
            d_uncovered_boxes_itr = uncovered_boxes.begin();
            d_uncovered_boxes_itr_end = uncovered_boxes.end();
            d_current_patch_id = this_itr->getBoxId(); 
            id_found = true;
            break;
         }
      }
      if (id_found) {
         break;
      } else {
         ++d_level_num;
      }
   }

   if (id_found) {
      setIteratorItem();
   } else {
      if (d_item) {
         delete d_item;
         d_item = nullptr;
      }
      if (d_flattened_hierarchy && d_allocated_flattened_hierarchy) {
         delete d_flattened_hierarchy;
      }
      d_flattened_hierarchy = nullptr;
      if (d_item) {
         delete d_item;
         d_item = nullptr;
      }
   }
}

void
UncoveredBoxIterator::setIteratorItem()
{
   // Get the current uncovered box.
   const Box& cur_box = *d_uncovered_boxes_itr;
   std::shared_ptr<PatchLevel> this_level =
      d_hierarchy->getPatchLevel(d_level_num);

   // Update item with the current originating patch and the current box.
   if (d_item) {
      d_item->first =
         this_level->getPatch(d_current_patch_id);
      d_item->second = cur_box;
   } else {
      d_item =
         new std::pair<std::shared_ptr<Patch>, Box>(
            this_level->getPatch(d_current_patch_id),
            cur_box);
   }
}

}
}
