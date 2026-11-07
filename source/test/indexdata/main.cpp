/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include <cassert>
#include <cstdlib>

#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/pdat/IndexVariable.cpp"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/pdat/CellData.cpp"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexData.cpp"
#include "SAMRAI/pdat/IndexDataFactory.h"
#include "SAMRAI/pdat/IndexDataFactory.cpp"

#include <list>
#include <memory>

using namespace SAMRAI;
using namespace hier;
using namespace pdat;
using namespace tbox;


#define NN 10

class Item
{

public:
   Item()
   {
   }

   ~Item()
   {
   }

   void copySourceItem(
      const hier::Index& idx,
      const hier::IntVector& src_offset,
      const Item& src_item)
   {
      NULL_USE(idx);
      NULL_USE(src_offset);
      for (int n = 0; n < NN; ++n) {
         x[n] = src_item.x[n];
      }
   }

   size_t getDataStreamSize()
   {
      return 0;
   }

   void packStream(
      MessageStream& stream)
   {
      NULL_USE(stream);
   }

   void unpackStream(
      MessageStream& stream,
      const hier::IntVector offset)
   {
      NULL_USE(stream);
      NULL_USE(offset);
   }

   void putToRestart(
      std::shared_ptr<tbox::Database> dbase)
   {
      NULL_USE(dbase);
   }
   void getFromRestart(
      std::shared_ptr<tbox::Database> dbase)
   {
      NULL_USE(dbase);
   }

   double x[NN];
};

int main(
   int argc,
   char* argv[])
{
   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();

   {
      tbox::Dimension dim(2);

      Index box_lo = Index(dim, 0);
      Index box_hi = Index(dim, 100);
      Box box(box_lo, box_hi, BlockId(0));

      srand(1);

      hier::IntVector v(dim, 0);
      hier::IntVector ghosts(dim, 0);

      /******************************************************************************
      * InedxData interface tests.
      ******************************************************************************/
      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         Item* item = new Item;

         v[0] = 0;
         v[1] = 0;
         Index idx(v);
         idx_data.addItemPointer(idx, item);

         // isElement()
         assert(idx_data.isElement(idx));
         v[0] = 1;
         v[1] = 0;
         Index idx2(v);
         assert(!idx_data.isElement(idx2));
         v[0] = 0;
         v[1] = 1;
         Index idx3(v);
         assert(!idx_data.isElement(idx3));

         // addItem()/getItem()
         assert(idx_data.getItem(idx) == item);

         assert(idx_data.getNumberOfItems() == 1);

         // removeItem()
         idx_data.removeItem(idx);
         assert(!idx_data.isElement(idx));

         assert(idx_data.getNumberOfItems() == 0);
      }

      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         Item* item = new Item;

         v[0] = 0;
         v[1] = 0;
         Index idx(v);
         idx_data.addItem(idx, *item);
         delete item;

         // isElement()
         assert(idx_data.isElement(idx));
         v[0] = 1;
         v[1] = 0;
         Index idx2(v);
         assert(!idx_data.isElement(idx2));
         v[0] = 0;
         v[1] = 1;
         Index idx3(v);
         assert(!idx_data.isElement(idx3));

         // addItem()/getItem()
         assert(idx_data.getNumberOfItems() == 1);

         // removeItem()
         idx_data.removeItem(idx);
         assert(!idx_data.isElement(idx));

         assert(idx_data.getNumberOfItems() == 0);
      }

      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         Item* item = new Item;
         v[0] = 0;
         v[1] = 0;
         Index idx(v);
         idx_data.replaceAddItem(idx, *item);
         delete item;

         // isElement()
         assert(idx_data.isElement(idx));
         v[0] = 1;
         v[1] = 0;
         Index idx2(v);
         assert(!idx_data.isElement(idx2));
         v[0] = 0;
         v[1] = 1;
         Index idx3(v);
         assert(!idx_data.isElement(idx3));

         // addItem()/getItem()
         assert(idx_data.getNumberOfItems() == 1);

         item = new Item;
         idx_data.replaceAddItem(idx, *item);
         delete item;

         assert(idx_data.getNumberOfItems() == 1);

         // removeItem()
         idx_data.removeItem(idx);
         assert(!idx_data.isElement(idx));
      }

      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         // getNumberItems()

         v[0] = 0;
         v[1] = 0;
         Index idx1(v);
         idx_data.addItemPointer(idx1, new Item);

         v[0] = 1;
         v[1] = 0;
         Index idx2(v);
         idx_data.addItemPointer(idx2, new Item);

         assert(idx_data.getNumberOfItems() == 2);

         // remove 1
         idx_data.removeItem(idx1);
         assert(idx_data.getNumberOfItems() == 1);

         // replace 1 at same index, no change
         idx_data.addItemPointer(idx2, new Item);
         assert(idx_data.getNumberOfItems() == 1);
      }

      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         // removeInsideBox()
         v[0] = 2;
         v[1] = 2;
         Index lo(v);

         v[0] = 3;
         v[1] = 5;
         Index hi(v);

         Box box1(lo, hi, BlockId(0));
         hier::Box::iterator biend(box1.end());
         for (Box::iterator bi(box1.begin()); bi != biend; ++bi) {

            Index idx = *bi;

            idx_data.addItemPointer(idx, new Item);

         }

         assert(static_cast<size_t>(idx_data.getNumberOfItems()) == box1.size());

         idx_data.removeInsideBox(box1);

         assert(idx_data.getNumberOfItems() == 0);
      }
      {
         IndexData<Item, pdat::CellGeometry> idx_data(box, ghosts);

         // removeAllItems()
         v[0] = 0;
         v[1] = 0;
         Index lo(v);

         v[0] = 1;
         v[1] = 1;
         Index hi(v);

         Box box1(lo, hi, BlockId(0));
         hier::Box::iterator biend(box1.end());
         for (Box::iterator bi(box1.begin()); bi != biend; ++bi) {

            Index idx = *bi;

            idx_data.addItemPointer(idx, new Item);

         }

         assert(static_cast<size_t>(idx_data.getNumberOfItems()) == box1.size());

         idx_data.removeAllItems();

         assert(idx_data.getNumberOfItems() == 0);
      }

      {
         // copy() where src and dst are same box

         IndexData<Item, pdat::CellGeometry> src(box, ghosts);
         IndexData<Item, pdat::CellGeometry> dst(box, ghosts);

         v[0] = 0;
         v[1] = 0;
         Index lo(v);

         v[0] = 1;
         v[1] = 1;
         Index hi(v);

         Box box1(lo, hi, BlockId(0));
         hier::Box::iterator biend(box1.end());
         for (Box::iterator bi(box1.begin()); bi != biend; ++bi) {
            src.addItemPointer(*bi, new Item);
         }

         assert(static_cast<size_t>(src.getNumberOfItems()) == box1.size());
         assert(static_cast<size_t>(dst.getNumberOfItems()) == 0);

         dst.copy(src);

         assert(dst.getNumberOfItems() == src.getNumberOfItems());
      }

      {

         // copy() where src and dst partially overlap, and only
         // some of src's items are contained in overlap.

         v[0] = 0;
         v[1] = 0;
         Index lo_src(v);

         v[0] = 2;
         v[1] = 2;
         Index hi_src(v);

         Box box_src(lo_src, hi_src, BlockId(0));
         IndexData<Item, pdat::CellGeometry> src(box_src, ghosts);

         // Two of these three items should end up in dst
         v[0] = 0;
         v[1] = 0;
         Index idx_item1(v);
         src.addItemPointer(idx_item1, new Item);

         v[0] = 1;
         v[1] = 1;
         Index idx_item2(v);
         src.addItemPointer(idx_item2, new Item);

         v[0] = 2;
         v[1] = 2;
         Index idx_item3(v);
         src.addItemPointer(idx_item3, new Item);

         v[0] = 1;
         v[1] = 1;
         Index lo_dst(v);
         v[0] = 3;
         v[1] = 3;
         Index hi_dst(v);
         Box box_dst(lo_dst, hi_dst, BlockId(0));

         IndexData<Item, pdat::CellGeometry> dst(box_dst, ghosts);

         assert(src.getNumberOfItems() == 3);
         assert(dst.getNumberOfItems() == 0);

         dst.copy(src);

         assert(dst.getNumberOfItems() == 2);
      }

      {
         // copy() with overlap argument

         // src
         // x . . . . . x "
         // .   .   . 3 . " // 2
         // . . . . . . . "
         // .   . 2 .   . " // 1
         // . . . . . . . "
         // . 1 .   .   . " // 0
         // x . . . . . x "

         // dst orig
         // . . x . . . x "
         // .   . 4 .   . " // 2
         // . . . . . . . "
         // .   .   .   . " // 1
         // . . x . . . x "
         // .   .   .   . " // 0
         // . . . . . . . "

         // dst expected
         // . . x . . . x "
         // .   .   . 3 . "
         // . . . . . . . "
         // .   . 2 .   . "
         // . . x . . . x "
         // .   .   .   . "
         // . . . . . . . "

         v[0] = 0;
         v[1] = 0;
         Index lo_src(v);
         v[0] = 2;
         v[1] = 2;
         Index hi_src(v);
         Box box_src(lo_src, hi_src, BlockId(0));
         IndexData<Item, pdat::CellGeometry> src(box_src, ghosts);

         // Two of these three items should end up in dst
         v[0] = 0;
         v[1] = 0;
         Index idx_item1(v);
         src.addItemPointer(idx_item1, new Item);

         v[0] = 1;
         v[1] = 1;
         Index idx_item2(v);
         src.addItemPointer(idx_item2, new Item);

         v[0] = 2;
         v[1] = 2;
         Index idx_item3(v);
         src.addItemPointer(idx_item3, new Item);

         v[0] = 1;
         v[1] = 1;
         Index lo_dst(v);

         v[0] = 2;
         v[1] = 2;
         Index hi_dst(v);
         Box box_dst(lo_dst, hi_dst, BlockId(0));

         IndexData<Item, pdat::CellGeometry> dst(box_dst, ghosts);

         // This item should be removed
         v[0] = 1;
         v[1] = 2;
         Index idx_item4(v);
         dst.addItemPointer(idx_item4, new Item);

         assert(src.getNumberOfItems() == 3);
         assert(dst.getNumberOfItems() == 1);

         IntVector src_offset(dim, 0);
         BoxContainer boxes(box_src);
         boxes.pushFront(box_dst);
         BoxContainer intersection(box_src * box_dst);
         CellOverlap overlap(intersection, hier::Transformation(src_offset));

         BoxContainer dst_boxlist(overlap.getDestinationBoxContainer());

         dst.copy(src, overlap);

         assert(dst.getNumberOfItems() == 2);
      }

      {
         // copy(): Same as test7 using copy2 which reverses src and dst

         v[0] = 0;
         v[1] = 0;
         Index lo_src(v);

         v[0] = 2;
         v[1] = 2;
         Index hi_src(v);
         Box box_src(lo_src, hi_src, BlockId(0));
         IndexData<Item, pdat::CellGeometry> src(box_src, ghosts);

         // Two of these three items should end up in dst
         v[0] = 0;
         v[1] = 0;
         Index idx_item1(v);
         src.addItemPointer(idx_item1, new Item);

         v[0] = 1;
         v[1] = 1;
         Index idx_item2(v);
         src.addItemPointer(idx_item2, new Item);

         v[0] = 2;
         v[1] = 2;
         Index idx_item3(v);
         src.addItemPointer(idx_item3, new Item);

         v[0] = 1;
         v[1] = 1;
         Index lo_dst(v);

         v[0] = 3;
         v[1] = 3;
         Index hi_dst(v);
         Box box_dst(lo_dst, hi_dst, BlockId(0));

         IndexData<Item, pdat::CellGeometry> dst(box_dst, ghosts);

         assert(src.getNumberOfItems() == 3);
         assert(dst.getNumberOfItems() == 0);

         src.copy2(dst);

         assert(dst.getNumberOfItems() == 2);
      }

      {
         v[0] = 0;
         v[1] = 0;
         Index lo(v);

         v[0] = 2;
         v[1] = 2;
         Index hi(v);
         Box data_box(lo, hi, BlockId(0));
         IndexData<Item, pdat::CellGeometry> data(data_box, ghosts);

         // Add three items
         v[0] = 0;
         v[1] = 0;
         Index idx_item1(v);
         data.addItemPointer(idx_item1, new Item);

         v[0] = 0;
         v[1] = 1;
         Index idx_item2(v);
         data.addItemPointer(idx_item2, new Item);

         v[0] = 2;
         v[1] = 1;
         Index idx_item3(v);
         data.addItemPointer(idx_item3, new Item);

         int count = 0;
         IndexIterator<Item, pdat::CellGeometry> itend(data, false);
         for (IndexIterator<Item, pdat::CellGeometry> it(data, true);
              it != itend; ++it) {
            ++count;
         }
         assert(3 == count);
      }

      int size = 100;
      {
         std::shared_ptr<tbox::Timer> timer(
            tbox::TimerManager::getManager()->
            getTimer("IndexDataAppendItemSequential", true));

         tbox::plog << "Begin Timing" << std::endl;

         Index lo = Index(dim, 0);
         Index hi = Index(dim, size);
         Box data_box(lo, hi, BlockId(0));

         IndexData<Item, pdat::CellGeometry> idx_data(data_box, ghosts);

         timer->start();

         for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
               v[0] = i;
               v[1] = j;
               Index idx(v);

               Item new_item;
               idx_data.appendItem(idx, new_item);
            }
         }

         size_t numberOfItems = idx_data.getNumberOfItems();
         timer->stop();

         tbox::plog << numberOfItems << std::endl;

         tbox::plog.precision(16);

         tbox::plog << "IndexData appendItem Sequential insert time : "
                    << timer->getTotalWallclockTime() << std::endl;

         tbox::plog << "End Timing" << std::endl;
      }

      {
         std::shared_ptr<tbox::Timer> timer(
            tbox::TimerManager::getManager()->
            getTimer("IndexDataAppendItemPointerSequential", true));

         tbox::plog << "Begin Timing" << std::endl;

         Index lo = Index(dim, 0);
         Index hi = Index(dim, size);
         Box data_box(lo, hi, BlockId(0));

         IndexData<Item, pdat::CellGeometry> idx_data(data_box, ghosts);

         timer->start();

         for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
               v[0] = i;
               v[1] = j;
               Index idx(v);

               Item* new_item = new Item();
               idx_data.appendItemPointer(idx, new_item);
            }
         }

         timer->stop();

         tbox::plog.precision(16);

         tbox::plog << "IndexData appendItemPointer sequential insert time : "
                    << timer->getTotalWallclockTime() << std::endl;

         tbox::plog << "End Timing" << std::endl;
      }

      size = 100;
      int num_inserts = 100000;

      {
         std::shared_ptr<tbox::Timer> timer(
            tbox::TimerManager::getManager()->
            getTimer("IndexDataAppendItemRandom", true));

         tbox::plog << "Begin Timing" << std::endl;

         Index lo = Index(dim, 0);
         Index hi = Index(dim, size);
         Box data_box(lo, hi, BlockId(0));

         IndexData<Item, pdat::CellGeometry> idx_data(data_box, ghosts);

         timer->start();

         for (int n = 0; n < num_inserts; ++n) {
            int i = rand() % size;
            int j = rand() % size;

            v[0] = i;
            v[1] = j;
            Index idx(v);

            Item new_item;
            idx_data.appendItem(idx, new_item);
         }

         size_t numberOfItems = idx_data.getNumberOfItems();
         timer->stop();

         tbox::plog << numberOfItems << std::endl;

         tbox::plog.precision(16);

         tbox::plog << "IndexData appendItem random insert time : "
                    << timer->getTotalWallclockTime() << std::endl;

         tbox::plog << "End Timing" << std::endl;
      }

      {
         std::shared_ptr<tbox::Timer> timer(
            tbox::TimerManager::getManager()->
            getTimer("IndexDataAppendItemPointerRandom", true));

         tbox::plog << "Begin Timing" << std::endl;

         Index lo = Index(dim, 0);
         Index hi = Index(dim, size);
         Box data_box(lo, hi, BlockId(0));

         IndexData<Item, pdat::CellGeometry> idx_data(data_box, ghosts);

         timer->start();

         for (int n = 0; n < num_inserts; ++n) {
            int i = rand() % size;
            int j = rand() % size;

            v[0] = i;
            v[1] = j;
            Index idx(v);

            Item* new_item = new Item();
            idx_data.appendItemPointer(idx, new_item);
         }

         timer->stop();

         tbox::plog.precision(16);

         tbox::plog << "IndexData appendItemPointer random insert time : "
                    << timer->getTotalWallclockTime() << std::endl;

         tbox::plog << "End Timing" << std::endl;
      }

      size = 100;

      {
         std::shared_ptr<tbox::Timer> timer(
            tbox::TimerManager::getManager()->
            getTimer("IndexDataReplace", true));

         tbox::plog << "Begin Timing" << std::endl;

         Index lo = Index(dim, 0);
         Index hi = Index(dim, size);
         Box data_box(lo, hi, BlockId(0));

         IndexData<Item, pdat::CellGeometry> idx_data(data_box, ghosts);

         timer->start();

         for (int n = 0; n < num_inserts; ++n) {
            int i = rand() % size;
            int j = rand() % size;

            v[0] = i;
            v[1] = j;
            Index idx(v);

            Item* new_item = new Item();
            idx_data.replaceAddItemPointer(idx, new_item);
         }

         timer->stop();

         tbox::plog.precision(16);

         tbox::plog << "IndexData replaceAddItemPointer random insert time : "
                    << timer->getTotalWallclockTime() << std::endl;

         tbox::plog << "End Timing" << std::endl;
      }
   }

   tbox::pout << "PASSED" << std::endl;

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();

   exit(0);
}
