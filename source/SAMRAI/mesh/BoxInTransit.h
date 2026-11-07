/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer.
 *
 ************************************************************************/

#ifndef included_mesh_BoxInTransit
#define included_mesh_BoxInTransit

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/mesh/TransitLoad.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MessageStream.h"

namespace SAMRAI {
namespace mesh {

/*!
 * @brief A Box moving around during load balancing.
 *
 * The purpose of the BoxInTransit is to associate extra data with
 * a Box as it is broken up and passed from process to process.  A
 * BoxInTransit is a Box going through these changes.  It has a
 * current work load and an orginating Box.
 */
struct BoxInTransit {

   /*!
    * @brief Constructor
    *
    * @param[in] dim
    */
   explicit BoxInTransit(
      const tbox::Dimension& dim):
      d_box(dim),
      d_orig_box(dim) {
   }

   /*!
    * @brief Construct a new BoxInTransit from an originating box.
    *
    * @param[in] origin
    */
   BoxInTransit(
      const hier::Box& origin):
      d_box(origin),
      d_orig_box(origin),
      d_boxload(static_cast<double>(origin.size()))
   {
      d_boxsize = d_boxload;
   }

   /*!
    * @brief Construct new object like an existing object but with a new ID.
    *
    * @param[in] other
    * @param[in] box
    * @param[in] rank
    * @param[in] local_id
    *
    * @param[in] load Box load.  If omitted, set to box's volume.
    */
   BoxInTransit(
      const BoxInTransit& other,
      const hier::Box& box,
      int rank,
      hier::LocalId local_id,
      double load = -1.):
      d_box(box, local_id, rank),
      d_orig_box(other.d_orig_box),
      d_boxsize(static_cast<double>(d_box.size()))
   {
      d_boxload = load >= 0 ? load : other.d_boxload;
      if (d_boxload > 0) {
         d_corner_weights = other.d_corner_weights;
      }
   }

   /*!
    * @brief Copy constructor
    *
    * @param[in] other   object to be copied
    */
   BoxInTransit(
      const BoxInTransit& other):
      d_box(other.d_box),
      d_orig_box(other.d_orig_box),
      d_boxload(other.d_boxload),
      d_boxsize(other.d_boxsize),
      d_corner_weights(other.d_corner_weights)
   {
   }

   /*!
    * @brief Assignment operator
    *
    * @param[in] other
    */
   BoxInTransit& operator = (const BoxInTransit& other) {
      d_box = other.d_box;
      d_orig_box = other.d_orig_box;
      d_boxload = other.d_boxload;
      d_boxsize = other.d_boxsize;
      d_corner_weights = other.d_corner_weights;
      return *this;
   }

   //! @brief Return the owner rank.
   int getOwnerRank() const {
      return d_box.getOwnerRank();
   }

   //! @brief Return the LocalId.
   hier::LocalId getLocalId() const {
      return d_box.getLocalId();
   }

   //! @brief Return the Box.
   hier::Box& getBox() {
      return d_box;
   }

   //! @brief Return the Box.
   const hier::Box& getBox() const {
      return d_box;
   }

   //! @brief Return the original Box.
   hier::Box& getOrigBox() {
      return d_orig_box;
   }

   //! @brief Return the original Box.
   const hier::Box& getOrigBox() const {
      return d_orig_box;
   }

   //! @brief Return the Box load.
   double getLoad() const {
      return d_boxload;
   }

   //! @brief Set the Box load.
   void setLoad(double boxload) {
      d_boxload = boxload;
   }

   //! @brief Set the corner weights vector
   void setCornerWeights(const std::vector<double>& corner_weights) {
      d_corner_weights = corner_weights;
   }

   //! @brief Get the box size
   double getSize() const {
      return d_boxsize;
   }

   //! @brief Get reference to the corner weights vector
   const std::vector<double>& getCornerWeights() const {
      return d_corner_weights;
   }

   //! @brief Return whether this is an original box.
   bool isOriginal() const {
      return d_box.isIdEqual(d_orig_box);
   }

   //! @brief Put self into a MessageStream.
   void putToMessageStream(tbox::MessageStream& mstream) const {
      d_box.putToMessageStream(mstream);
      d_orig_box.putToMessageStream(mstream);
      mstream << d_boxload;
      mstream << d_boxsize;
      int num_corners = static_cast<int>(d_corner_weights.size());
      mstream << num_corners;
      for (int nc = 0; nc < num_corners; ++nc) {
         mstream << d_corner_weights[nc]; 
      }
   }

   //! @brief Set attributes according to data in a MessageStream.
   void getFromMessageStream(tbox::MessageStream& mstream) {
      d_box.getFromMessageStream(mstream);
      d_orig_box.getFromMessageStream(mstream);
      mstream >> d_boxload;
      mstream >> d_boxsize;
      int num_corners = 0;
      mstream >> num_corners;
      d_corner_weights.resize(num_corners);
      for (int nc = 0; nc < num_corners; ++nc) {
         mstream >> d_corner_weights[nc]; 
      }
   }

   /*!
    * @brief Insert BoxInTransit into an output stream.
    */
   friend std::ostream& operator << (std::ostream& co,
                                     const BoxInTransit& r) {
      co << r.d_box
      << r.d_box.numberCells() << '|' << r.d_box.size() << '-'
      << r.d_orig_box
      << r.d_orig_box.numberCells() << '|' << r.d_orig_box.size();
      return co;
   }

private:

   //! @brief Current Box
   hier::Box d_box;

   //! @brief Originating Box (the oldest one leading to this one).
   hier::Box d_orig_box;

   //! @brief Work load in this box.
   double d_boxload;

   //! @brief Size of this box
   double d_boxsize;

   //! @brief Vector holding fractional workload weighting for each corner
   std::vector<double> d_corner_weights;
};

}
}

#endif
