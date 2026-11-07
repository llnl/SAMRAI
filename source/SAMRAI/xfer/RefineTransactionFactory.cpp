/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface for factory objects that create transactions for
 *                refine schedules.
 *
 ************************************************************************/
#include "SAMRAI/xfer/RefineTransactionFactory.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

RefineTransactionFactory::RefineTransactionFactory()
{
}

RefineTransactionFactory::~RefineTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Default no-op implementations of optional virtual functions.
 *
 *************************************************************************
 */

void
RefineTransactionFactory::setTransactionTime(
   double fill_time)
{
   NULL_USE(fill_time);
}

}
}
