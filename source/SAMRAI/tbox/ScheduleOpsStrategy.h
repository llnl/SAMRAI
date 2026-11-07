/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for Schedule operations
 *
 ************************************************************************/

#ifndef included_tbox_ScheduleOpsStrategy
#define included_tbox_ScheduleOpsStrategy

#include "SAMRAI/SAMRAI_config.h"


namespace SAMRAI {
namespace tbox {

/*!
 * ScheduleOpsStrategy is a virtual base class that provides pure virtual
 * methods that can serve as callback hooks into user code before or after
 * certain operations inside tbox::Schedule.
 *
 * tbox::Schedule contains a method where a pointer to a concrete instance
 * of this class can be passed into the tbox::Schedule object.
 *
 * @see Schedule 
 */

class ScheduleOpsStrategy
{
public:
   /*!
    * @brief Default constructor
    */
   ScheduleOpsStrategy();

   /*!
    * @brief Virtual destructor
    */
   virtual ~ScheduleOpsStrategy();

   /*!
    * @brief Callback for user code operations before the communication
    * operations in Schedule::communicate
    */
   virtual void preCommunicate() = 0;

   /*!
    * @brief Callback for user code operations before local data copies
    */
   virtual void preCopy() = 0;

   /*!
    * @brief Callback for user code operations before the packing of
    * message streams
    */
   virtual void prePack() = 0;

   /*!
    * @brief Callback for user code operations before the packing of
    * message streams
    */
   virtual void preUnpack() = 0;

   /*!
    * @brief Callback for user code operations after the communication
    * operations in Schedule::communicate
    */
   virtual void postCommunicate() = 0;

   /*!
    * @brief Callback for user code operations after local data copies
    */
   virtual void postCopy() = 0;

   /*!
    * @brief Callback for user code operations after the packing of
    * message streams
    */
   virtual void postPack() = 0;

   /*!
    * @brief Callback for user code operations after the packing of
    * message streams
    */
   virtual void postUnpack() = 0;

   /*!
    * @brief Callback to direct tbox::Schedule to defer all MPI sends until
    * all MessageStreams are packed.
    *
    * If this returns true, all subsequent packed streams will not be sent
    * until after postPack is called.
    */
   virtual bool deferMessageSend() = 0;

   /*!
    * @brief Tell whether a host-device synchronize is needed.
    *
    * This virtual method is provided so that the child implementations
    * can be able to indicate whether their operations have been left
    * in a state that requires tbox::Schedule to call a host-device
    * synchronize.  For example if child implmentation of this class
    * has just called a synchronize, it can return false, which tells
    * tbox::Schedule to skip the next synchronize call in its own code.
    *
    * The default implementation returns true as the safest default behavior,
    * as it will make tbox::Schedule call synchronize every time there
    * is any uncertainty about whether it is required.
    */
   virtual bool needSynchronize()
   {
      return true;
   }

private:
   ScheduleOpsStrategy(
      const ScheduleOpsStrategy&);              // not implemented
   ScheduleOpsStrategy&
   operator = (
      const ScheduleOpsStrategy&);              // not implemented

};

}
}

#endif
