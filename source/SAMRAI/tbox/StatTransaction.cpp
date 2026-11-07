/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction structure for statistic data copies
 *
 ************************************************************************/

#include "SAMRAI/tbox/StatTransaction.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

StatTransaction::StatTransaction(
   const std::shared_ptr<Statistic>& stat,
   int src_proc_id,
   int dst_proc_id):
   d_stat(stat),
   d_src_id(src_proc_id),
   d_dst_id(dst_proc_id)
{
}

StatTransaction::~StatTransaction()
{
}

bool
StatTransaction::canEstimateIncomingMessageSize()
{
   return d_stat->canEstimateDataStreamSize();
}

size_t
StatTransaction::computeIncomingMessageSize()
{
   return d_stat->getDataStreamSize();
}

size_t
StatTransaction::computeOutgoingMessageSize()
{
   return d_stat->getDataStreamSize();
}

int
StatTransaction::getSourceProcessor()
{
   return d_src_id;
}

int
StatTransaction::getDestinationProcessor()
{
   return d_dst_id;
}

void
StatTransaction::packStream(
   MessageStream& stream)
{
   d_stat->packStream(stream);
}

void
StatTransaction::unpackStream(
   MessageStream& stream)
{
   d_stat->unpackStream(stream);
}

void
StatTransaction::copyLocalData()
{
   // Nothing to do here!
}

void
StatTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Stat Transaction" << std::endl;
   stream << "   source processor:   " << d_src_id << std::endl;
   stream << "   destination processor:   " << d_dst_id << std::endl;
   stream << "   stat name:   " << d_stat->getName() << std::endl;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
