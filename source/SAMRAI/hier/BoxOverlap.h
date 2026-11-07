/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class that describes intersections between AMR boxes
 *
 ************************************************************************/

#ifndef included_hier_BoxOverlap
#define included_hier_BoxOverlap

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Transformation.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class BoxOverlap is an abstract base class used to represent a region
 * where data will be communicated between two AMR patches.
 *
 * The exact form of the overlap will be determined by the subclass, which
 * will be implemented to handle a particular patch data class.  For example,
 * the rules for data intersection for face-centered data are different
 * from that for cell-centered data.
 *
 * The BoxOverlap class provides three functions.  First, it serves
 * as a base class that can answer the question whether an intersection
 * is empty, and is therefore useful for determining communication
 * dependencies.  Second, it holds information about any transformation (such
 * as a periodic shift) between a source and destination patch.  Third, it
 * is a storage location for the exact form of the intersection of the data
 * residing on two patches, which can be quite complicated (for example,
 * for face centered data).  To access the information about the intersection,
 * type-safe type casting should be used to access the subclass and its
 * member functions.
 *
 * @see BoxGeometry
 */

class BoxOverlap
{
public:
   /*!
    * @brief The default constructor for BoxOverlap.
    */
   BoxOverlap() = default;

   /*!
    * @brief The virtual destructor.
    */
   virtual ~BoxOverlap() = default;

   /*!
    * @brief Return true if overlap object represents an empty data
    * overlap; i.e., there is no data intersection between patches.
    *
    * Note that two patches may communicate even if they do not intersect in
    * the underlying AMR index space (e.g., if data values exist at the
    * a faces or corners of cells).
    *
    * @return  Returns true for empty overlaps and false otherwise.
    */
   virtual bool
   isOverlapEmpty() const = 0;

   /*!
    * @brief Return the offset between the destination and source index spaces.
    *
    * The destination index space is the source index space shifted
    * by this amount.
    *
    * @return  The offset.
    */
   virtual const IntVector&
   getSourceOffset() const = 0;

   /*!
    * @brief Refturn the transformation between the destination and source
    * index spaces.
    *
    * The destination index space is the source index space rotated and shifted
    * by the returned transformation.
    *
    * @return  The transformation.
    */
   virtual const Transformation&
   getTransformation() const = 0;

   /*!
    * @brief Print BoxOverlap object data.
    *
    * @param[in]     os the std::ostream to which to print
    */
   virtual void
   print(
      std::ostream& os) const;

private:
   BoxOverlap(
      const BoxOverlap&) = delete;
   BoxOverlap&
   operator = (
      const BoxOverlap&) = delete;

};

}
}

#endif
