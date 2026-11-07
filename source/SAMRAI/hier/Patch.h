/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Patch container class for patch data objects
 *
 ************************************************************************/

#ifndef included_hier_Patch
#define included_hier_Patch

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A container for patch data objects defined over a box.
 *
 * Class Patch is a container for patch data objects defined over
 * some box region.  Since the SAMRAI convention is that the structured
 * mesh index space is always defined as cell-centerd, the box of a patch
 * may not match the box describing the patch data (in fact, it is not in
 * general).  Users must be clear about how they relate the patch box to
 * the data living on the patch.
 *
 * Each of the patch data objects that exist on the patch is a subclass
 * of the patch data pure virtual base class.  The construction
 * of the patch data objects on a patch are managed via the factory classes
 * maintained by the patch descriptor which is shared by all patches.
 *
 * By default, the patch constructor does not allocate space for the various
 * patch data components that live on the patch.  Individual components or sets
 * of components can be created or destroyed via patch member functions.
 *
 * @see Box
 * @see PatchDescriptor
 * @see PatchData
 * @see PatchDataFactory
 * @see PatchGeometry
 */
class Patch
{
public:
   /*!
    * @brief  Allocate a patch container over the box.
    *
    * @note
    * Patch data components are not allocated/instantiated.
    *
    * @param[in]  box
    * @param[in]  descriptor
    *
    * @pre box.getLocalId() >= 0
    */
   Patch(
      const Box& box,
      const std::shared_ptr<PatchDescriptor>& descriptor);

   /*!
    * @brief Virtual destructor for patch objects.
    */
   virtual ~Patch() = default;

   /*!
    * @brief Get the box over which the patch is defined.
    *
    * Note that the patch data objects on a patch are free to interpret
    * this box as appropriate for the "geometry" of their data
    * (e.g., face-centered, node-centered, etc.)
    *
    * @return The box over which this patch is defined.
    */
   const Box&
   getBox() const
   {
      return d_box;
   }

   /*!
    * @brief Get the GlobalId for this patch.
    *
    * The GlobalId is the patch's unique identifier within the
    * PatchLevel.  It is identical to the GlobalId of the Patch's
    * Box.
    *
    * @return GlobalId for this patch.
    */
   const GlobalId&
   getGlobalId() const
   {
      return d_box.getGlobalId();
   }

   /*!
    * @brief Get the patch's LocalId.
    *
    * The LocalId is the same as that of the Box used to
    * construct the Patch.
    *
    * @return The LocalId of this patch.
    */
   const LocalId&
   getLocalId() const
   {
      return d_box.getLocalId();
   }

   /*!
    * @brief Get the patch descriptor.
    *
    * The patch descriptor describes the patch data types that may exist on
    * the patch.
    *
    * @return the patch descriptor for this patch.
    */
   std::shared_ptr<PatchDescriptor>
   getPatchDescriptor() const
   {
      return d_descriptor;
   }

   /*!
    * @brief Get the patch data identified by the specified id.
    *
    * Typically, this function should only be called to access
    * patch data objects that have already been explicitly allocated through
    * a call to one of the patch allocation routines.  A NULL pointer will
    * be returned if the patch data object is not allocated.
    *
    * @return a pointer to the patch data object associated with the specified
    * identifier.
    *
    * @param[in]  id
    *
    * @pre (id >= 0) && (id < numPatchData())
    */
   std::shared_ptr<PatchData>
   getPatchData(
      const int id) const
   {
      TBOX_ASSERT((id >= 0) && (id < numPatchData()));
      return d_patch_data[id];
   }

   /*!
    * @brief Get the patch data associated with the specified variable and
    * context.
    *
    * Typically, this function should only be called
    * to access patch data objects that have already been explicitly allocated
    * through a call to one of the patch allocation routines.  A NULL pointer
    * will be returned if the patch data object is not allocated.
    *
    * @return a pointer to the patch data object associated with the specified
    * variable and context.
    *
    * @param[in]  variable
    * @param[in]  context
    *
    * @pre getDim() == variable->getDim()
    * @pre (id >= 0) && (id < numPatchData())
    */
   std::shared_ptr<PatchData>
   getPatchData(
      const std::shared_ptr<Variable>& variable,
      const std::shared_ptr<VariableContext>& context) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *variable);
      int id = VariableDatabase::getDatabase()->
         mapVariableAndContextToIndex(variable, context);
      TBOX_ASSERT((id >= 0) && (id < numPatchData()));
      return d_patch_data[id];
   }

   /*!
    * @brief Set the patch data Pointer associated with the specified
    * identifier.
    *
    * @note
    * This member function must be used with caution.  It can only
    * be called for patch data indices that were previously allocated through
    * one of the patch routines.  This member function does not check to see
    * whether the patch data types are already allocated, consistent, or
    * whether they have the same factory.  So, for example, a face centered
    * data type could be assigned to a location reserved for a cell centered
    * data type (with potentially different box and ghost cell width).
    *
    * @param[in]  id
    * @param[out] data
    *
    * @pre getDim() == data->getDim();
    * @pre (id >= 0) && (id < numPatchData())
    */
   void
   setPatchData(
      const int id,
      const std::shared_ptr<PatchData>& data)
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *data);
      TBOX_ASSERT((id >= 0) && (id < numPatchData()));
      d_patch_data[id] = data;
   }

   /*!
    * @brief Check whether the specified component has been allocated.
    *
    * @return True if the component has been allocated.
    *
    * @param[in]  id
    */
   bool
   checkAllocated(
      const int id) const
   {
      return (id < numPatchData()) && d_patch_data[id];
   }

   /*!
    * @brief Get the size of the patch data
    *
    * @return the amount of memory needed to allocate the specified component.
    *
    * @param[in]  id
    */
   size_t
   getSizeOfPatchData(
      const int id) const
   {
      return getPatchDescriptor()->getPatchDataFactory(id)->getSizeOfMemory(d_box);
   }

   /*!
    * @brief Get the size of the patch data for all components specified
    *
    * @return the amount of memory needed to allocate the specified components.
    *
    * @param[in]  components
    */
   size_t
   getSizeOfPatchData(
      const ComponentSelector& components) const;

   /*!
    * @brief Allocate the specified component on the patch.
    * @par Assertions
    * An assertion will result if the component is already allocated.
    * This provides a key bit of debugging information that may be useful
    * to application developers.
    *
    * @param[in]  id
    * @param[in]  time
    *
    * @pre (id >= 0) && (id < getPatchDescriptor()->getMaxNumberRegisteredComponents())
    */
   void
   allocatePatchData(
      const int id,
      const double time = 0.0);

   /*!
    * @brief Allocate the specified components on the patch.
    * @par Assertions
    * An assertion will result if any requested component is
    * already allocated.  This provides a key bit of debugging information
    * that may be useful to application developers.
    *
    * @param[in]  components
    * @param[in]  time
    */
   void
   allocatePatchData(
      const ComponentSelector& components,
      const double time = 0.0);

   /*!
    * @brief Deallocate the specified component.
    *
    * This component will need to be reallocated before its next use.
    *
    * @param[in]  id
    *
    * @pre (id >= 0) && (id < getPatchDescriptor()->getMaxNumberRegisteredComponents())
    */
   void
   deallocatePatchData(
      const int id)
   {
      TBOX_ASSERT((id >= 0) &&
         (id < getPatchDescriptor()->getMaxNumberRegisteredComponents()));
      if (id < numPatchData()) {
         d_patch_data[id].reset();
      }
   }

   /*!
    * @brief Deallocate the specified components.
    *
    * These components will need to be reallocated before their next use.
    *
    * @param[in]  components
    */
   void
   deallocatePatchData(
      const ComponentSelector& components);

   /*!
    * @brief Set the geometry specification for the patch.
    *
    * This includes patch boundary information and grid data.
    *
    * @param[in]  geometry
    */
   void
   setPatchGeometry(
      const std::shared_ptr<PatchGeometry>& geometry)
   {
      d_patch_geometry = geometry;
   }

   /*!
    * @brief Get the patch geometry
    *
    * @return pointer to patch geometry object.
    */
   std::shared_ptr<PatchGeometry>
   getPatchGeometry() const
   {
      return d_patch_geometry;
   }

   /*!
    * @brief Set the timestamp value for the specified patch component.
    *
    * @param[in]  timestamp
    * @param[in]  id
    *
    * @pre (id >= 0) && (id < numPatchData())
    * @pre getPatchData(id)
    */
   void
   setTime(
      const double timestamp,
      const int id)
   {
      TBOX_ASSERT((id >= 0) && (id < numPatchData()));
      TBOX_ASSERT(d_patch_data[id]);
      d_patch_data[id]->setTime(timestamp);
   }

   /*!
    * @brief Set the timestamp value for the specified patch components.
    *
    * @param[in]  timestamp
    * @param[in]  components
    */
   void
   setTime(
      const double timestamp,
      const ComponentSelector& components);

   /*!
    * @brief Set the timestamp value for all allocated patch components.
    *
    * @param[in]  timestamp
    */
   void
   setTime(
      const double timestamp);

   /*!
    * @brief Get the level number of the patch level in the patch hierarchy
    * where this patch resides.
    *
    * @note
    * This value can be a valid level number (i.e., >=0) even
    * when the level is not in a hierarchy.  In this case the level
    * number will be the number of the hierarchy level matching the
    * index space of the patch level holding this patch.  If the patch
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.
    *
    * @return The valid level number >= 0 or undefined (-1) if the patch
    * level does not align with the index space of a level in the hierarchy.
    *
    * @see inHierarchy()
    */
   int
   getPatchLevelNumber() const
   {
      return d_patch_level_number;
   }

   /*!
    * @brief Set the patch level number for this patch.
    *
    * When the index space of the level owning this patch aligns with the
    * index space of some valid hierarchy level, sets the patch level number
    * The default level number is -1.
    *
    * @param[in]  level_number
    */
   void
   setPatchLevelNumber(
      const int level_number)
   {
      d_patch_level_number = level_number;
   }
   /*!
    * @brief Determine if the level holding this patch resides in a hierarchy.
    *
    * @return True if the level holding this patch resides in a hierarchy,
    * otherwise false.
    */
   bool
   inHierarchy() const
   {
      return d_patch_in_hierarchy;
   }

   /*!
    * @brief Set a flag determining if the patch resides in a hierarchy.
    *
    * Set to true if the level holding this patch resides in a hierarchy;
    * false otherwise.  The default setting is false.
    *
    * @param[in]  in_hierarchy The flag that indicates whether the patch
    *             resides in a hierarchy.
    */
   void
   setPatchInHierarchy(
      bool in_hierarchy)
   {
      d_patch_in_hierarchy = in_hierarchy;
   }

   /*!
    * @brief Get the patch data items from the restart database.
    *
    * Patch state is read in from the database and all patch
    * data objects specified in the PatchDataRestartManager are created.
    *
    * The class version and restart file version must be equal.
    *
    * @par Assertions
    * Checks that data retrieved from the database are of the type
    * expected, and that the patch_number read in from the database
    * matches the patch number assigned to this Patch.
    * @note
    * A warning will be printed to the log file if some patch data components
    * that were requested through the PatchDataRestartManager are not found in
    * the database.
    *
    * @param[in]  restart_db
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * @brief Write patch data and other patch information to the restart
    * database.
    *
    * Class version number and the state of the patch object are written.
    * Patch data objects specified in the PatchDataRestartManager are also
    * written.
    *
    * @param[in]  restart_db
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * @brief Print a patch (for debugging).
    *
    * Depth is kept for consistency with other recursivePrint methods,
    * but is not used because this is the lowest level of recursion
    * currently supported.
    *
    * @param[in,out]  os
    * @param[in]  border
    * @param[in]  depth
    */
   int
   recursivePrint(
      std::ostream& os,
      const std::string& border = std::string(),
      int depth = 0) const;

   /*!
    * @brief Get the dimension of this object.
    *
    * @return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_box.getDim();
   }

   /*!
    * @brief Output patch information (box and number of components).
    *
    * @param[in,out]  s The output stream
    * @param[in]      patch
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Patch& patch);

   /*!
    * @brief Number of PatchData objects defined on this Patch.
    */
   int
   numPatchData() const
   {
      return static_cast<int>(d_patch_data.size());
   }

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_PATCH_VERSION;

   /*
    * Copy constructor and assignment operator are not implemented.
    */
   Patch(
      const Patch&);
   Patch&
   operator = (
      const Patch&);

   /*
    * The box defining the extent of this patch.
    */
   Box d_box;

   std::shared_ptr<PatchDescriptor> d_descriptor;

   std::shared_ptr<PatchGeometry> d_patch_geometry;

   std::vector<std::shared_ptr<PatchData> > d_patch_data;

   int d_patch_level_number;

   bool d_patch_in_hierarchy;

};

}
}

#endif
