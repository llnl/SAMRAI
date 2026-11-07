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

#include "SampleIndexData.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Dimension.h"

#include <iostream>

/*
 *************************************************************************
 *
 * Constructor providing cell index.
 *
 *************************************************************************
 */

using namespace SAMRAI;

SampleIndexData::SampleIndexData():
   d_dummy_int(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

SampleIndexData::~SampleIndexData()
{
}

/*
 *************************************************************************
 *
 * Set dummy int data
 *
 *************************************************************************
 */
void SampleIndexData::setInt(
   const int dummy)
{
   d_dummy_int = dummy;
}

/*
 *************************************************************************
 *
 *  Return dummy int data
 *
 *************************************************************************
 */
int SampleIndexData::getInt() const
{
   return d_dummy_int;
}

/*
 *************************************************************************
 *
 * The copySourceItem() method allows SampleIndexData to be a templated
 * data type for IndexData - i.e. IndexData<SampleIndexData>.
 *
 *************************************************************************
 */
void SampleIndexData::copySourceItem(
   const hier::Index& index,
   const hier::IntVector& src_offset,
   const SampleIndexData& src_item)
{
   NULL_USE(index);
   NULL_USE(src_offset);
   d_dummy_int = src_item.d_dummy_int;
}

/*
 *************************************************************************
 *
 * The getDataStreamSize(), packStream(), and unpackStream() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to communicate SampleIndexData,
 * specifying how many bytes will be packed during the "packStream()"
 * method.
 *
 *************************************************************************
 */

size_t SampleIndexData::getDataStreamSize()
{
   return 0;
}

void SampleIndexData::packStream(
   tbox::MessageStream& stream)
{
   NULL_USE(stream);
}

void SampleIndexData::unpackStream(
   tbox::MessageStream& stream,
   const hier::IntVector& offset)
{
   NULL_USE(stream);
   NULL_USE(offset);
}

/*
 *************************************************************************
 *
 * The putToRestart() and getFromRestart() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to write/read SampleIndexData,
 * data to/from the restart database.
 *
 *************************************************************************
 */

void SampleIndexData::putToRestart(
   std::shared_ptr<tbox::Database>& restart_db) const
{
   NULL_USE(restart_db);
}

void SampleIndexData::getFromRestart(
   std::shared_ptr<tbox::Database>& restart_db)
{
   NULL_USE(restart_db);
}

/*
 *****************************************************************
 *
 *  Templates used for SampleIndexData
 *
 *****************************************************************
 */

//#include "SampleIndexData.h"
//#include "SAMRAI/pdat/IndexData.cpp"
//#include "SAMRAI/pdat/IndexDataFactory.cpp"
//#include "SAMRAI/pdat/IndexVariable.cpp"
//#include "SAMRAI/pdat/CellGeometry.h"
//
//namespace SAMRAI {
//
//template class pdat::SparseData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::SparseDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataNode<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexIterator<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexVariable<SampleIndexData, pdat::CellGeometry>;
//template class std::vector<SampleIndexData>;
//template class std::vector<pdat::IndexDataNode<SampleIndexData,
//                                               pdat::CellGeometry> >;
//template class std::shared_ptr<pdat::IndexData<SampleIndexData,
//                                                 pdat::CellGeometry> >;
//template class std::shared_ptr<pdat::IndexVariable<SampleIndexData,
//                                                     pdat::CellGeometry> >;
//template class std::shared_ptr<pdat::IndexDataFactory<SampleIndexData,
//                                                        pdat::CellGeometry> >;
//
//}
