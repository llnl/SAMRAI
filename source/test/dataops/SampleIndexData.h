/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Example user defined index data type used in indx_dataops
 *                test.
 *
 ************************************************************************/

#ifndef included_SampleIndexDataXD
#define included_SampleIndexDataXD

#include "SAMRAI/SAMRAI_config.h"

//#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IOStream.h"

#include <memory>

/**
 * The SampleClass struct holds some dummy data and methods.  It's intent
 * is to indicate how a user could construct their own index data type.
 */

using namespace SAMRAI;

class SampleIndexData
{
public:
   SampleIndexData();

   /**
    * The destructor for SampleIndexData.
    */
   ~SampleIndexData();

   /**
    * Sets a dummy integer in this class.
    */
   void
   setInt(
      const int dummy);

   /**
    * Returns a dummy integer in this class.
    */
   int
   getInt() const;

   /**
    * The copySourceItem() method allows SampleIndexData to be a templated
    * data type for IndexData - i.e. IndexData<SampleIndexData>.  In
    * addition to this method, the other methods that must be defined are
    * getDataStreamSize(), packStream(), unpackStream() for communication,
    * putToRestart(), getFromRestart() for restart.  These are
    * described below.
    */
   void
   copySourceItem(
      const hier::Index& index,
      const hier::IntVector& src_offset,
      const SampleIndexData& src_item);

   /**
    * The following functions enable parallel communication with
    * SampleIndexDatas.  They are used in SAMRAI communication infrastructure
    * to specify the number of bytes of data stored in each SampleIndexData
    * object, and to pack and unpack the data to the specified stream.
    */
   size_t
   getDataStreamSize();
   void
   packStream(
      tbox::MessageStream& stream);
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::IntVector& offset);

   /**
    * These functions are used to read/write SampleIndexData data to/from
    * restart.
    */
   void
   getFromRestart(
      std::shared_ptr<tbox::Database>& restart_db);
   void
   putToRestart(
      std::shared_ptr<tbox::Database>& restart_db) const;

private:
   /*
    * Dummy int data
    */
   int d_dummy_int;

   // ADD ANY OTHER DATA HERE
};
#endif
