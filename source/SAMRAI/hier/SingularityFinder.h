/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class for finding multiblock singularities
 *
 ************************************************************************/

#ifndef included_hier_SingularityFinder
#define included_hier_SingularityFinder

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"


namespace SAMRAI {
namespace hier {

/*!
 * @brief Class SingularityFinder implements an algorithm to find
 * multiblock singularities based on the relationships between blocks that
 * touch on an entire face.
 */

class SingularityFinder
{
public:
   SingularityFinder(
      const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   ~SingularityFinder() = default;

   /*!
    * @brief Find the singularities based on the relationship between
    * face neighbors.
    *
    * The singularity information is computed based on the face_neighbors
    * input. In the nested map of maps face_neighbors, the key for the outer
    * map is a BoxId of a Box in the domain_boxes input (the reference Box).
    * For each reference Box,the inner map contains information about all its
    * face neighbors.
    *
    * Each face neighbor has a BoxId as its identifying key and
    * an integer to indicate which face of the reference Box the neighbor
    * touches. The formula for the integer identifier is:  For each normal
    * direction D = 0,...,NDIM-1, the lower face identifier is 2*D and the
    * upper face identifier is 2*D+1.
    *
    * The domain_boxes input is a representation of the physical domain that
    * must be set up such that all face neighbor relationships between
    * Boxes touch each other across the entirety of each Box face.
    *
    * @param[out] singularity_blocks Output for the computed singularities.
    *                                The outer set is the container of all
    *                                singularities, while each inner set
    *                                contains all of the BlockIds for the
    *                                blocks that touch one singularity.
    * @param[in] domain_boxes   Representation of the physical domain.
    * @param[in] grid_geometry  The grid geometry describing the multiblock
    *                           layout.
    * @param[in] face_neighbors Input information about face neighbor
    *                           relationships, described above.
    */
   void
   findSingularities(
      std::set<std::set<BlockId> >& singularity_blocks,
      const BoxContainer& domain_boxes,
      const BaseGridGeometry& grid_geometry,
      const std::map<BoxId, std::map<BoxId, int> >& face_neighbors);

private:
   struct Block;
   struct Face;
   struct Edge;
   struct Point;

   /*!
    * @brief Private struct Block holds inforation about all
    * the faces, edges, and corner points of a block.
    */
   struct Block {
      Block(
         const tbox::Dimension& dim) {

         if (dim.getValue() == 1) {
            d_nfaces = 2;
            d_nedges = 0;
            d_npoints = 0;
         } else if (dim.getValue() == 2) {
            d_nfaces = 4;
            d_nedges = 0;
            d_npoints = 4;
         } else if (dim.getValue() == 3) {
            d_nfaces = 6;
            d_nedges = 12;
            d_npoints = 8;
         } else {
            d_nfaces = 0;
            d_nedges = 0;
            d_npoints = 0;
         }

         if (d_nfaces) {
            d_face.resize(d_nfaces);
         }
         if (d_nedges) {
            d_edge.resize(d_nedges);
         }
         if (d_npoints) {
            d_point.resize(d_npoints);
         }

      }

      int d_nfaces;
      int d_nedges;
      int d_npoints;
      std::vector<std::shared_ptr<Face> > d_face;
      std::vector<std::shared_ptr<Edge> > d_edge;
      std::vector<std::shared_ptr<Point> > d_point;
   };

   /*!
    * @brief Private struct Face
    *
    * Each block face is represented by this struct.  When a face is shared by
    * two blocks, a single instance of this struct is shared by both blocks.
    */
   struct Face {

      Face() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int, int> d_block_to_face;
   };

   /*!
    * @brief Private struct Edge
    *
    * In 3D, each block edge is represented by this struct.  When an edge is
    * shared by multiple blocks, a single instance of this struct is shared by
    * all of the blocks.
    */
   struct Edge {

      Edge() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int, int> d_block_to_edge;
   };

   /*!
    * @brief Private struct Point
    *
    * Each corner point of a block is represented by this struct.  When a
    * point is shared by multiple blocks, a single instance of this struct is
    * shared by all of the blocks.
    */
   struct Point {

      Point() {
         d_bdry = false;
      }

      bool d_bdry;
      std::set<int> d_blocks;
      std::map<int, int> d_block_to_point;
   };

   /*
    * Unimplemented default constructor.
    */
   SingularityFinder();

   /*!
    * @brief Find all shared information between two Boxes.
    *
    * Given two BoxIds representing two Boxes that share a face,
    * compute all the information about shared points and edges.
    *
    * @param[in] id_a           ID of Box A
    * @param[in] id_b           ID of Box B
    * @param[in] face_a         Face of Box A that touches Box B
    * @param[in] face_b         Face of Box B that touches Box A
    * @param[in] domain_boxes   Boxes for the physical domain
    * @param[in] grid_geometry
    */
   void
   connect(
      const BoxId& id_a,
      const BoxId& id_b,
      int face_a,
      int face_b,
      const BoxContainer& domain_boxes,
      const BaseGridGeometry& grid_geometry);

   /*!
    * @brief Determine which block faces are not shared by two blocks.  These
    * are boundary faces and must be excluded from the singularity computation.
    */
   void
   findBoundaryFaces();

   /*!
    * @brief Determine which edges are shared between two boxes.
    */
   void
   findCoincidentEdges(
      std::map<int, int>& map_of_edges,
      const BoxId& id_a,
      const BoxId& id_b,
      int facea,
      int faceb,
      const BoxContainer& domain_boxes,
      const BaseGridGeometry& grid_geometry);

   /*!
    * @brief Determine which corner points are shared between two boxes.
    */
   void
   findCoincidentPoints(
      std::map<int, int>& map_of_points,
      const BoxId& id_a,
      const BoxId& id_b,
      int facea,
      const BoxContainer& domain_boxes,
      const BaseGridGeometry& grid_geometry);

   tbox::Dimension d_dim;

   std::vector<std::shared_ptr<Block> > d_blocks;
   std::vector<std::shared_ptr<Face> > d_faces;
   std::list<std::shared_ptr<Edge> > d_edges;
   std::list<std::shared_ptr<Point> > d_points;

   /*!
    * @brief Static vectors that serve as a table to map face identifiers
    * to edge and node integer identifiers.
    */
   static std::vector<std::vector<int> > s_face_edges;
   static std::vector<std::vector<int> > s_face_nodes;

};

}
}

#endif
