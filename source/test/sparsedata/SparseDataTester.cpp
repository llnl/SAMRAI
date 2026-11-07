/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test class for SparseData (implementation).
 *
 ************************************************************************/
#include "SparseDataTester.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputDatabase.h"

#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

using namespace SAMRAI;

namespace sam_test {

SparseDataTester::SparseDataTester(
   const tbox::Dimension& dim):
   d_initialized(false),
   d_dim(dim)
{
}

SparseDataTester::~SparseDataTester() {
   d_sparse_data->clear();
}

bool
SparseDataTester::testConstruction()
{
   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));

   hier::IntVector vec(d_dim, 0);
   hier::IntVector ghosts(d_dim, 0);
   std::vector<std::string> dkeys;
   _getDblKeys(dkeys);
   std::vector<std::string> ikeys;
   _getIntKeys(ikeys);

   std::shared_ptr<SparseDataType> sparse(
      new SparseDataType(box, ghosts, dkeys, ikeys));
   d_sparse_data = sparse;

   d_sparse_data->printNames(tbox::plog);

   bool passed = true;
   if (!d_sparse_data->empty()) {
      passed = false;
      tbox::perr << "Sparse data should be empty and is not" << std::endl;
   } else {
      d_initialized = true;
   }

   SparseDataType::iterator iter = d_sparse_data->registerIndex(hi);

   SparseDataType::AttributeIterator index_iter(d_sparse_data->begin(hi)),
   index_iterend(d_sparse_data->end(hi));

   if (index_iter != index_iterend) {
      passed = false;
      tbox::perr << "something is wrong with index2's list of attributes"
                 << std::endl;
   }

   d_sparse_data->clear();
   return passed;

}

bool
SparseDataTester::testCopy()
{

   //_fillObject(d_sparse_data);

   // ensure d_sparse_data is empty before we start
   d_sparse_data->clear();
   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());
   bool success = _testCopy(d_sparse_data, sample);
   // clean up d_sparse_data
   d_sparse_data->clear();
   sample->clear();
   return success;

}

bool
SparseDataTester::testCopy2()
{
   // ensure the tester's copy of d_sparse_data is empty before
   // we start
   d_sparse_data->clear();
   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());
   _fillObject(sample);

   bool success = _testCopy(sample, d_sparse_data);
   // clean up d_sparse_data
   d_sparse_data->clear();
   sample->clear();
   return success;

}

bool
SparseDataTester::testAdd()
{
   bool success = true;
   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());

   success = success && (sample->empty() ? true : false);
   if (success)
      _fillObject(sample);

   success = success && (!sample->empty() ? true : false);
   if (success)
      sample->clear();

   success = success && (sample->empty() ? true : false);
   return success;
}

bool
SparseDataTester::testRemove()
{
   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());
   _fillObject(sample);
   bool success = (!sample->empty() ? true : false);

   hier::Index idx = _getRandomIndex();

   int size = d_sparse_data->size();
   SparseDataType::iterator iter = d_sparse_data->begin();
   while (iter != d_sparse_data->end()) {
      if (iter.getIndex() == idx) {
         sample->remove(iter);
      } else {
         ++iter;
      }
   }

   int newsize = sample->size();
   success = success && (newsize = size - 1) ? true : false;

   if (success) {
      sample->clear();
      success = success && (sample->empty() ? true : false);
   }
 
   return success;

}

bool
SparseDataTester::testPackStream(int num_indices)
{
   bool success = true;

   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());
   _fillObject(sample, num_indices);

   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));
   hier::BoxContainer blist(box);
   hier::Transformation trans(hier::IntVector::getZero(d_dim));

   pdat::CellOverlap overlap(blist, trans);

   size_t strsize = sample->getDataStreamSize(overlap);

   SparseDataType::iterator iter(sample.get());

   tbox::MessageStream str(strsize, tbox::MessageStream::Write);
   sample->packStream(str, overlap);
   tbox::plog << "Printing sample1" << std::endl;
   sample->printAttributes(tbox::plog);

   std::shared_ptr<SparseDataType> sample2(_createEmptySparseData());
   tbox::MessageStream upStr(strsize, tbox::MessageStream::Read,
                             str.getBufferStart());

   sample2->unpackStream(upStr, overlap);

//   SparseDataType::iterator iter2(sample2.get());

   sample2->printNames(tbox::plog);
   tbox::plog << "Printing sample2" << std::endl;
   sample2->printAttributes(tbox::plog);
/*
 * for ( ; iter != sample->end() && iter2 != sample2->end(); ++iter, ++iter2) {
 *    tbox::plog << "iter1 node: " << std::endl;
 *    tbox::plog << iter;
 *    tbox::plog << "iter2 node: " << std::endl;
 *    tbox::plog << iter2;
 *    if (!iter.equals(iter2)) {
 *       success = false;
 *    }
 *    tbox::plog << std::endl;
 * }
 */
   if (num_indices == 0) {
      if (!sample2->empty()) {
         success = false;
      }
   } else {
      for ( ; iter != sample->end(); ++iter) {
         tbox::plog << "iter1 node: " << std::endl;
         tbox::plog << iter;
      
         bool found = false;
         SparseDataType::iterator iter2(sample2.get());
         for ( ; iter2 != sample2->end(); ++iter2) {
            if (iter.equals(iter2)) {
               found = true;
               tbox::plog << "iter2 node: " << std::endl;
               tbox::plog << iter2;
               break;
            }
         }
   
         if (!found) {
            success = false;
         }
         tbox::plog << std::endl;
      }

      const pdat::DoubleAttributeId did =
         sample2->getDblAttributeId("double_key_0");

      if (!sample2->isValid(did)) {
         success = false;
      }
   }
   sample->clear();
   sample2->clear();
   return success;
}

bool
SparseDataTester::testDatabaseInterface()
{
   bool success = true;

   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());
   _fillObject(sample);
   std::shared_ptr<tbox::Database> input_db(
      new tbox::InputDatabase("input_db"));
   sample->putToRestart(input_db);

   std::shared_ptr<SparseDataType> sample2(_createEmptySparseData());
   sample2->getFromRestart(input_db);

   SparseDataType::iterator iter1(sample.get());

   for ( ; iter1 != sample->end(); ++iter1) {
      bool found = false;
      SparseDataType::iterator iter2(sample2.get());
      for ( ; iter2 != sample2->end(); ++iter2) {
         if (iter1.equals(iter2)) {
            found = true;
            break;
         }
      }

      if (!found) {
         success = false;
      }
   }

   sample->clear();
   sample2->clear();
   return success;
}

void
SparseDataTester::testTiming()
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("SparseDataAddItem", true));

   hier::IntVector v(d_dim, 0);
   tbox::plog << "Begin Timing" << std::endl;
   std::shared_ptr<SparseDataType> sample(_createEmptySparseData());

   timer->start();
   SparseDataType::iterator iter;
   for (int i = 0; i < 100; ++i) {
      for (int j = 0; j < 100; ++j) {
         v[0] = i;
         v[1] = j;
         hier::Index idx(v);
         double dvalues[DSIZE];
         _getDblValues(dvalues);
         int ivalues[ISIZE];
         _getIntValues(ivalues);

         iter = sample->registerIndex(idx);
         for (int m = i; m < j + 1; ++m) {
            iter.insert(dvalues, ivalues);
         }
      }
   }

   timer->stop();
   int numItems = sample->size();
   tbox::plog << numItems << std::endl;
   tbox::plog.precision(16);
   tbox::plog << "SparseData addItem insert time : "
              << timer->getTotalWallclockTime() << std::endl;
   tbox::plog << "End Timing" << std::endl;
   sample->clear();
}


bool
SparseDataTester::_testCopy(
   std::shared_ptr<SparseDataType> src,
   std::shared_ptr<SparseDataType> dst)
{
   bool success = true;
   src->copy(*dst);
   TBOX_ASSERT(src->size() == dst->size());
   pdat::SparseData<pdat::CellGeometry>::iterator me(src.get());
   pdat::SparseData<pdat::CellGeometry>::iterator me_end = src->end();
   pdat::SparseData<pdat::CellGeometry>::iterator other(dst.get());

   for ( ; me != me_end && success != false; ++me, ++other) {
      if (me != other) {
         success = false;
      }
   }
   return success;
}

hier::Index
SparseDataTester::_getRandomIndex()
{
   // return a random index in the range that we created
   // with _fillObject
   int val = (rand() % d_num_indices);
   hier::IntVector v(d_dim, 0);
   v[0] = val;
   v[1] = val;
   return hier::Index(v);
}

void
SparseDataTester::_fillObject(
   std::shared_ptr<SparseDataType> sparse_data,
   int num_indices)
{
   d_num_indices = num_indices;

   hier::IntVector v(d_dim, 0);
   double* dvalues = new double[DSIZE];
   _getDblValues(dvalues);
   int* ivalues = new int[ISIZE];
   _getIntValues(ivalues);

   SparseDataType::iterator iter;
   for (int i = 1; i < d_num_indices; ++i) {
      v[0] = i;
      v[1] = i;
      hier::Index idx(v);
      iter = sparse_data->registerIndex(idx);
      for (int k = 0; k < i; ++k) {
         iter.insert(dvalues, ivalues);
      }
   }

   delete[] dvalues;
   delete[] ivalues;
}

std::shared_ptr<pdat::SparseData<pdat::CellGeometry> >
SparseDataTester::_createEmptySparseData()
{
   hier::Index lo = hier::Index(d_dim, 0);
   hier::Index hi = hier::Index(d_dim, 100);
   hier::Box box(lo, hi, hier::BlockId(0));
   hier::IntVector ghosts(d_dim, 0);

   std::vector<std::string> dkeys;
   _getDblKeys(dkeys);
   std::vector<std::string> ikeys;
   _getIntKeys(ikeys);
   return std::shared_ptr<SparseDataType>(
             new SparseDataType(box, ghosts, dkeys, ikeys));
}

void
SparseDataTester::_getDblKeys(std::vector<std::string>& keys)
{
   for (int i = 0; i < DSIZE; ++i) {
      std::stringstream key_name;
      key_name << "DOUBLE_KEY_" << i;
      keys.push_back(key_name.str());
   }
}

void
SparseDataTester::_getIntKeys(std::vector<std::string>& keys)
{
   for (int i = 0; i < ISIZE; ++i) {
      std::stringstream key_name;
      key_name << "INTEGER_KEY_" << i;
      keys.push_back(key_name.str());
   }
}

void
SparseDataTester::_getDblValues(double* values)
{
   for (int i = 0; i < DSIZE; ++i) {
      values[i] = (double)i;
   }
}

void
SparseDataTester::_getIntValues(int* values)
{
   for (int i = 0; i < ISIZE; ++i) {
      values[i] = i;
   }
}


} // end namespace sam_test
