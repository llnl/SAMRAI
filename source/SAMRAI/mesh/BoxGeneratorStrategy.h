/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box generation routines.
 *
 ************************************************************************/

#ifndef included_mesh_BoxGeneratorStrategy
#define included_mesh_BoxGeneratorStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"

namespace SAMRAI {
namespace mesh {

/**
 * Class BoxGeneratorStrategy is an abstract base class that defines
 * a Strategy pattern interface for operations to build boxes that cover a
 * collection of tagged cells on a single AMR patch hierarchy level.
 *
 * @see hier::PatchLevel
 */

class BoxGeneratorStrategy
{
public:
   /**
    * Default constructor.
    */
   BoxGeneratorStrategy();

   /**
    * Virtual destructor.
    */
   virtual ~BoxGeneratorStrategy();

   /*!
    * @brief Cluster tags using the DLBG interfaces.
    *
    * @param tag_to_new_width [in] Width that tag_to_new should have.
    * If implementation does not provide this width for tag_to_new,
    * then it should set the width to zero.
    *
    * @param[out] new_box_level BoxLevel containing Boxes of clustered tagged
    * cells.
    * @param[out] tag_to_new Connector from the tagged to the new BoxLevels.
    * @param[in] tag_level Tagged PatchLevel.
    * @param[in] tag_data_index Index of PatchData used to denote tagging.
    * @param[in] tag_val Value of PatchData indicating a tagged cell.
    * @param[in] bound_boxes Collection of Boxes describing the bounding box
    * of each block in the tag level.
    * @param[in] min_box Smallest box size resulting from clustering.
    * @param[in] tag_to_new_width Width of tag_to_new Connector.
    */
   virtual void
   findBoxesContainingTags(
      std::shared_ptr<hier::BoxLevel>& new_box_level,
      std::shared_ptr<hier::Connector>& tag_to_new,
      const std::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::BoxContainer& bound_boxes,
      const hier::IntVector& min_box,
      const hier::IntVector& tag_to_new_width) = 0;

   /*!
    * @brief Set a minimum cell request value.
    *
    * This serves as a virtual interface for child classes to receive a
    * value that requests a minimum number of cells for the boxes produced
    * by a box generation class.  A default no-op implementation is
    * provided, as not all child classes may use such a feature.
    */
   virtual void
   setMinimumCellRequest(
      size_t num_cells);

   /*!
    * @brief Set the ratio to the new level that will be produced.
    */
   virtual void
   setRatioToNewLevel( 
      const hier::IntVector& ratio) = 0;


private:
   // The following are not implemented:
   BoxGeneratorStrategy(
      const BoxGeneratorStrategy&);
   BoxGeneratorStrategy&
   operator = (
      const BoxGeneratorStrategy&);

};

}
}
#endif
