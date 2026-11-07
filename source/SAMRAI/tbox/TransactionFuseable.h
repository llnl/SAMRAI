/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   base class for fuseable schedule transactions
 *
 ************************************************************************/

#ifndef included_tbox_TransactionFuseable
#define included_tbox_TransactionFuseable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"

#include <iostream>

/*!
 * @brief Class TransactionFuseable inherits from tbox::Transaction,
 * providing a type of transaction that can use a kernel fuser object.
 *
 * This class optionally holds a pointer to a StagedKernelFusers object
 * that can be used by its child classes.
 */


namespace SAMRAI {
namespace tbox {

class TransactionFuseable :
   public Transaction
{
public:

   /*!
    * @brief Default constructor
    */
   TransactionFuseable();

   /*!
    * @brief Destructor
    */
   virtual ~TransactionFuseable();


   /*!
    * @brief Provide a pointer to a StagedKernelFusers object.
    */
   void setKernelFuser(StagedKernelFusers* fuser);

   /*!
    * @brief Get a pointer to the StagedKernelFusers object.
    */
   StagedKernelFusers* getKernelFuser();

private:
   StagedKernelFusers* d_fuser{nullptr};

};

}
}

#endif
