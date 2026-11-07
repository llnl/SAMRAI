/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_SingularityPatchStrategy
#define included_xfer_SingularityPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include <memory>

namespace SAMRAI {

namespace hier {
class Box;
class Connector;
class BaseGridGeometry;
class BoundaryBox;
class Patch;
class PatchLevel;
}

namespace xfer {

/*!
 * @brief Abstract base class for setting ghost data when refining at
 * a multiblock singularity.
 *
 * This class would be needed where ever you need a
 * RefinePatchStrategy and have a multiblock mesh with singularities.
 * Otherwise, you can ignore it.
 *
 * To use this base class, inherit it in the concrete class where you
 * inherit RefinePatchStrategy and implement its pure virtual methods.
 *
 * SingularityPatchStrategy provides an interface for a user to supply
 * methods for application-specific refining of data at block
 * singularities in the mesh.  A concrete subclass must define
 * the method fillSingularityBoundaryConditions.
 */

class SingularityPatchStrategy
{
public:
   /*!
    * @brief Constructor.
    */
   SingularityPatchStrategy();

   /*!
    * @brief Destructor
    */
   virtual ~SingularityPatchStrategy();

   /*!
    * @brief Set the ghost data at a multiblock singularity.
    *
    * This virtual method allows for a user-defined implemenation to
    * fill ghost data at ghost regions located at reduced or enhanced
    * connectivity multiblock singularities.  The encon_level and
    * dst_to_encon arguments may be ignored if the patch touches no
    * enhanced connectivity singularities.
    *
    * The patches in encon level are in the coordinate system of the blocks
    * where they originated, not in that of the destination patch, so the
    * filling operation must take into account the transformation between
    * blocks.
    *
    * @param patch The patch containing the data to be filled
    * @param encon_level  Level representing enhanced connectivity ghost
    *                     regions
    * @param dst_to_encon  Connector from destination level to encon_level
    * @param fill_box Box covering maximum amount of ghost cells to be filled
    * @param boundary_box BoundaryBox describing location of singularity in
    *                     relation to patch
    * @param[in] grid_geometry
    */
   virtual void
   fillSingularityBoundaryConditions(
      hier::Patch& patch,
      const hier::PatchLevel& encon_level,
      std::shared_ptr<const hier::Connector> dst_to_encon,
      const hier::Box& fill_box,
      const hier::BoundaryBox& boundary_box,
      const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry) = 0;

   bool
   needSynchronize()
   {
      bool flag = d_need_synchronize;
      d_need_synchronize = true;
      return flag;
   }

protected:

   void
   setNeedSingularitySynchronize(bool flag)
   {
      d_need_synchronize = flag;
   }

private:

   bool d_need_synchronize = true;


};

}
}

#endif
