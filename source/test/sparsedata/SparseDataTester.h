/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test class for SparseData.
 *
 ************************************************************************/
#ifndef included_SparseDataTester_h
#define included_SparseDataTester_h

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/pdat/SparseData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/hier/Index.h"

#include <string>
#include <vector>
#include <memory>

#define SPARSE_NUM_INDICES 5
namespace sam_test {

using namespace SAMRAI;

class SparseDataTester
{
public:
   SparseDataTester(
      const tbox::Dimension& dim);
   ~SparseDataTester();

   bool
   testConstruction();
   bool
   testAdd();
   bool
   testRemove();
   bool
   testCopy();
   bool
   testCopy2();
   void
   testTiming();
   bool
   testPackStream(int num_indices = SPARSE_NUM_INDICES);
   bool
   testDatabaseInterface();

private:
   static const int DSIZE = 7;
   static const int ISIZE = 3;

   typedef pdat::SparseData<pdat::CellGeometry> SparseDataType;

   void
   _fillObject(
      std::shared_ptr<SparseDataType> sparse_data,
      int num_indices = SPARSE_NUM_INDICES);
   void
   _getDblKeys(
      std::vector<std::string>& keys);
   void
   _getIntKeys(
      std::vector<std::string>& keys);
   void
   _getDblValues(
      double* values);
   void
   _getIntValues(
      int* values);
   bool
   _testCopy(
      std::shared_ptr<SparseDataType> src,
      std::shared_ptr<SparseDataType> dst);
   std::shared_ptr<SparseDataType>
   _createEmptySparseData();
   hier::Index
   _getRandomIndex();

   std::shared_ptr<SparseDataType> d_sparse_data;

   bool d_initialized;
   tbox::Dimension d_dim;

   int d_num_indices = 0;
};

} // end namespace sam_test

#endif
