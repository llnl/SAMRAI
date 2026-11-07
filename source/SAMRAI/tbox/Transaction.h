/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for all schedule transactions
 *
 ************************************************************************/

#ifndef included_tbox_Transaction
#define included_tbox_Transaction

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MessageStream.h"

#include <iostream>

namespace SAMRAI {
namespace tbox {

/**
 * Class Transaction describes a single communication between two
 * processors or a local data copy.  It is an abstract base class for each
 * data transaction in a communication schedule.
 */

class Transaction
{
public:
   /**
    * The constructor for transaction does nothing interesting.
    */
   Transaction();

   /**
    * The virtual destructor for transaction does nothing interesting.
    */
   virtual ~Transaction();

   /**
    * Return a boolean indicating whether this transaction can estimate
    * the size of an incoming message.  If this is false, then a different
    * communications protocol kicks in and the message size is transmitted
    * between nodes.
    *
    * Note that the message size estimate may be an overestimate but
    * should not be an underestimate.  Also, the receiver should never
    * estimate lower than the sender, because that may lead to
    * fatal MPI errors due to truncated messages.
    */
   virtual bool
   canEstimateIncomingMessageSize() = 0;

   /**
    * Return the amount of buffer space needed for the incoming message.
    * This routine is only called if the transaction can estimate the
    * size of the incoming message.
    *
    * @see canEstimateIncomingMessageSize().
    */
   virtual size_t
   computeIncomingMessageSize() = 0;

   /**
    * Return the buffer space needed for the outgoing message.
    *
    * @see canEstimateIncomingMessageSize().
    */
   virtual size_t
   computeOutgoingMessageSize() = 0;

   /**
    * Return the sending processor for the communications transaction.
    */
   virtual int
   getSourceProcessor() = 0;

   /**
    * Return the receiving processor for the communications transaction.
    */
   virtual int
   getDestinationProcessor() = 0;

   /**
    * Pack the transaction data into the message stream.
    */
   virtual void
   packStream(
      MessageStream& stream) = 0;

   /**
    * Unpack the transaction data from the message stream.
    */
   virtual void
   unpackStream(
      MessageStream& stream) = 0;

   /**
    * Perform the local data copy for the transaction.
    */
   virtual void
   copyLocalData() = 0;

   /**
    * Print out transaction information.
    */
   virtual void
   printClassData(
      std::ostream& stream) const = 0;

private:
   Transaction(
      const Transaction&);              // not implemented
   Transaction&
   operator = (
      const Transaction&);              // not implemented

};

}
}

#endif
