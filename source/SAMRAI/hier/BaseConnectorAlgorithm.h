/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Algorithms to work with maping Connectors.
 *
 ************************************************************************/
#ifndef included_hier_BaseConnectorAlgorithm
#define included_hier_BaseConnectorAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"

namespace SAMRAI {
namespace hier {

class BaseConnectorAlgorithm
{
protected:
   /*!
    * @brief Constructor
    */
   BaseConnectorAlgorithm() = default;

   /*!
    * @brief Destructor.
    */
   virtual ~BaseConnectorAlgorithm() = default;

   /*!
    * @brief Set up communication objects for use in privateBridge/Modify.
    */
   void
   setupCommunication(
      tbox::AsyncCommPeer<int> *& all_comms,
      tbox::AsyncCommStage& comm_stage,
      const tbox::SAMRAI_MPI& mpi,
      const std::set<int>& incoming_ranks,
      const std::set<int>& outgoing_ranks,
      const std::shared_ptr<tbox::Timer>& mpi_wait_timer,
      int& operation_mpi_tag,
      bool print_steps) const;

   /*!
    * @brief Pack referenced neighbors discovered during privateBridge/Modify
    * into message for one processor.
    */
   void
   packReferencedNeighbors(
      std::vector<int>& send_mesg,
      int idx_offset_to_ref,
      const BoxContainer& referenced_new_head_nabrs,
      const BoxContainer& referenced_new_base_nabrs,
      const tbox::Dimension& dim,
      bool print_steps) const;

   /*!
    * @brief Receive messages and unpack info sent from other processes.
    */
   void
   receiveAndUnpack(
      Connector& new_base_to_new_head,
      Connector* new_head_to_new_base,
      const std::set<int>& incoming_ranks,
      tbox::AsyncCommPeer<int>* all_comms,
      tbox::AsyncCommStage& comm_stage,
      const std::shared_ptr<tbox::Timer>& receive_and_unpack_timer,
      bool print_steps) const;

private:
   /*
    * Data length limit on first message of a communication.
    */
   static const int BASE_CONNECTOR_ALGORITHM_FIRST_DATA_LENGTH;

   //! @brief Unpack message sent by sendDiscoverytoOneProcess().
   void
   unpackDiscoveryMessage(
      const tbox::AsyncCommPeer<int>* incoming_comm,
      Connector& west_to_east,
      Connector* east_to_west,
      bool print_steps) const;
};

}
}

#endif // included_hier_BaseConnectorAlgorithm
