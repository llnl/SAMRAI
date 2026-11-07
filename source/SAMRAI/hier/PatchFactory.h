/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch classes
 *
 ************************************************************************/

#ifndef included_hier_PatchFactory
#define included_hier_PatchFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/tbox/Database.h"

#include <memory>

namespace SAMRAI {
namespace hier {

/**
 * Class PatchFactory is a factory object used to create new patches.
 * New types of patch objects can be introduced into the hierarchy through
 * derivation and re-defining the allocate member function.  There should
 * be no direct calls to the patch constructor (other than through the
 * patch factory).
 *
 * @see Patch
 */

class PatchFactory
{
public:
   /**
    * Construct a patch factory object.
    */
   PatchFactory() = default;

   /**
    * Virtual destructor for patch factory objects.
    */
   virtual ~PatchFactory() = default;

   /**
    * Allocate a patch with the specified domain and patch descriptor.
    */
   virtual std::shared_ptr<Patch>
   allocate(
      const Box& box_level_box,
      const std::shared_ptr<PatchDescriptor>& descriptor) const;

private:
   PatchFactory(
      const PatchFactory&) = delete;
   PatchFactory&
   operator = (
      const PatchFactory&) = delete;

};

}
}

#endif
