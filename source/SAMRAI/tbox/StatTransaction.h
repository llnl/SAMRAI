/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction structure for statistic data copies
 *
 ************************************************************************/

#ifndef included_tbox_StatTransaction
#define included_tbox_StatTransaction

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Transaction.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace tbox {

/**
 * A stattistic transaction represents a simple copy communication
 * transaction between two processors for sending and gathering statistic
 * information generated on different processors.
 *
 * @see Schedule
 * @see Transaction
 */

class StatTransaction:public Transaction
{
public:
   /**
    * Create a transaction for communicating local statistic information.
    * This transaction will be responsible for either: (1) packing a
    * message stream with statistic information if this processor number
    * is the same as the given source processor id, or (2) unpacking
    * tatistic information from the message stream if this processor number
    * is the same as the given destination processor.  The statistic pointer
    * passed in through the argument list must be non-null and will be used
    * as either the source or destination statistic depending on whether
    * case (1) or (2) applies.
    *
    * Note that generally this transaction class is used to pass information
    * between two different processors and unexpected behavior may result
    * if the source and destination processors are the same.  Also, note
    * that the copyLocalData() routine has an empty implementation.
    */
   StatTransaction(
      const std::shared_ptr<Statistic>& stat,
      int src_proc_id,
      int dst_proc_id);

   /**
    * The virtual destructor for the copy transaction releases all
    * memory associated with the transaction.
    */
   virtual ~StatTransaction();

   /**
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communications protocol kicks in and the message size is transmitted
    * between nodes.
    */
   virtual bool
   canEstimateIncomingMessageSize();

   /**
    * Return the amount of buffer space needed for the incoming message.
    * This routine is only called if the transaction can estimate the
    * size of the incoming message.
    */
   virtual size_t
   computeIncomingMessageSize();

   /**
    * Return the buffer space needed for the outgoing message.
    */
   virtual size_t
   computeOutgoingMessageSize();

   /**
    * Return the sending processor for the communications transaction.
    */
   virtual int
   getSourceProcessor();

   /**
    * Return the receiving processor for the communications transaction.
    */
   virtual int
   getDestinationProcessor();

   /**
    * Pack the transaction data into the message stream.
    */
   virtual void
   packStream(
      MessageStream& stream);

   /**
    * Unpack the transaction data from the message stream.
    */
   virtual void
   unpackStream(
      MessageStream& stream);

   /**
    * Perform the local data copy for the transaction.  This function
    * drops through as it is not needed.
    */
   virtual void
   copyLocalData();

   /**
    * Print transaction information to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& stream) const;

private:
   StatTransaction();                           // not implemented
   StatTransaction(
      const StatTransaction&);                  // not implemented
   StatTransaction&
   operator = (
      const StatTransaction&);                  // not implemented

   std::shared_ptr<Statistic> d_stat;
   int d_src_id;
   int d_dst_id;

};

}
}
#endif
