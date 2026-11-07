/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Class for finding multiblockSingularities
 *
 ************************************************************************/
#include "SAMRAI/hier/SingularityFinder.h"


namespace SAMRAI {
namespace hier {

std::vector<std::vector<int> > SingularityFinder::s_face_edges;
std::vector<std::vector<int> > SingularityFinder::s_face_nodes;

/*
 * ************************************************************************
 *
 * Constructor
 *
 * ************************************************************************
 */

SingularityFinder::SingularityFinder(
   const tbox::Dimension& dim):
   d_dim(dim)
{

   if (d_dim.getValue() == 3 && s_face_edges.empty()) {
      std::vector<int> edges;
      edges.resize(4);
      edges[0] = 0;
      edges[1] = 2;
      edges[2] = 4;
      edges[3] = 6;
      s_face_edges.push_back(edges);
      edges[0] = 1;
      edges[1] = 3;
      edges[2] = 5;
      edges[3] = 7;
      s_face_edges.push_back(edges);
      edges[0] = 0;
      edges[1] = 1;
      edges[2] = 8;
      edges[3] = 10;
      s_face_edges.push_back(edges);
      edges[0] = 2;
      edges[1] = 3;
      edges[2] = 9;
      edges[3] = 11;
      s_face_edges.push_back(edges);
      edges[0] = 4;
      edges[1] = 5;
      edges[2] = 8;
      edges[3] = 9;
      s_face_edges.push_back(edges);
      edges[0] = 6;
      edges[1] = 7;
      edges[2] = 10;
      edges[3] = 11;
      s_face_edges.push_back(edges);
   }

   if (s_face_nodes.empty()) {
      std::vector<int> nodes;
      if (d_dim.getValue() == 3) {
         nodes.resize(4);
         nodes[0] = 0;
         nodes[1] = 2;
         nodes[2] = 4;
         nodes[3] = 6;
         s_face_nodes.push_back(nodes);
         nodes[0] = 1;
         nodes[1] = 3;
         nodes[2] = 5;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         nodes[2] = 4;
         nodes[3] = 5;
         s_face_nodes.push_back(nodes);
         nodes[0] = 2;
         nodes[1] = 3;
         nodes[2] = 6;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         nodes[2] = 2;
         nodes[3] = 3;
         s_face_nodes.push_back(nodes);
         nodes[0] = 4;
         nodes[1] = 5;
         nodes[2] = 6;
         nodes[3] = 7;
         s_face_nodes.push_back(nodes);
      } else if (d_dim.getValue() == 2) {
         nodes.resize(2);
         nodes[0] = 0;
         nodes[1] = 2;
         s_face_nodes.push_back(nodes);
         nodes[0] = 1;
         nodes[1] = 3;
         s_face_nodes.push_back(nodes);
         nodes[0] = 0;
         nodes[1] = 1;
         s_face_nodes.push_back(nodes);
         nodes[0] = 2;
         nodes[1] = 3;
         s_face_nodes.push_back(nodes);
      }
   }
}

/*
 * ************************************************************************
 *
 * Destructor
 *
 * ************************************************************************
 */

/*
 * ************************************************************************
 *
 * Find the singularities.
 *
 * ************************************************************************
 */

void
SingularityFinder::findSingularities(
   std::set<std::set<BlockId> >& singularity_blocks,
   const BoxContainer& domain_boxes,
   const BaseGridGeometry& grid_geometry,
   const std::map<BoxId, std::map<BoxId, int> >& face_neighbors)
{
   if (face_neighbors.empty()) return;

   if (grid_geometry.getNumberBlocks() == 1) return;

   TBOX_ASSERT(singularity_blocks.empty());

   std::map<BoxId, std::map<BoxId, int> > unprocessed = face_neighbors;
   std::list<BoxId> to_be_processed;
   to_be_processed.push_back(face_neighbors.begin()->first);

   do {
      BoxId base_id = *(to_be_processed.begin());
      to_be_processed.pop_front();

      if (unprocessed.find(base_id) == unprocessed.end()) continue;

      TBOX_ASSERT(face_neighbors.find(base_id) != face_neighbors.end());
      const std::map<BoxId, int>& nbr_ids = face_neighbors.find(base_id)->second;
      for (std::map<BoxId, int>::const_iterator nbr_itr = nbr_ids.begin();
           nbr_itr != nbr_ids.end(); ++nbr_itr) {

         const std::pair<BoxId, int>& nbr_face = *nbr_itr;
         const BoxId& nbr_id = nbr_face.first;
         int facea = nbr_face.second;

         TBOX_ASSERT(face_neighbors.find(nbr_id) != face_neighbors.end());
         const std::map<BoxId, int>& nbr_of_nbr = face_neighbors.find(nbr_id)->second;
         TBOX_ASSERT(nbr_of_nbr.find(base_id) != nbr_of_nbr.end());
         int faceb = nbr_of_nbr.find(base_id)->second;

         if (unprocessed[base_id].find(nbr_id) != unprocessed[base_id].end()) {
            connect(base_id, nbr_id, facea, faceb, domain_boxes, grid_geometry);

            to_be_processed.push_back(nbr_id);

            unprocessed[nbr_id].erase(base_id);
            if (unprocessed[nbr_id].empty()) {
               unprocessed.erase(nbr_id);
            }
            unprocessed[base_id].erase(nbr_id);
            if (unprocessed[base_id].empty()) {
               unprocessed.erase(base_id);
               break;
            }
         }
      }
   } while (!unprocessed.empty());

   findBoundaryFaces();

   std::set<std::set<BlockId> > sing_set;

   for (std::list<std::shared_ptr<Edge> >::iterator e_itr = d_edges.begin();
        e_itr != d_edges.end(); ++e_itr) {
      if ((**e_itr).d_bdry) continue;
      int nblocks = static_cast<int>((**e_itr).d_blocks.size());
      bool enhanced = false;
      bool reduced = false;
      if (nblocks < 4) {
         reduced = true;
      } else if (nblocks > 4) {
         enhanced = true;
      }

      if (reduced || enhanced) {
         std::set<BlockId> sing_set;
         for (std::set<int>::iterator s_itr = (**e_itr).d_blocks.begin();
              s_itr != (**e_itr).d_blocks.end(); ++s_itr) {
            BoxId box_id(LocalId(*s_itr), 0);
            const BlockId& block_id =
               domain_boxes.find(Box(d_dim, box_id))->getBlockId();
            sing_set.insert(block_id);
         }
         singularity_blocks.insert(sing_set);
      }
   }

   for (std::list<std::shared_ptr<Point> >::iterator p_itr = d_points.begin();
        p_itr != d_points.end(); ++p_itr) {
      if ((**p_itr).d_bdry) continue;
      int nblocks = static_cast<int>((**p_itr).d_blocks.size());
      bool enhanced = false;
      bool reduced = false;
      if (nblocks < (1 << d_dim.getValue())) {
         reduced = true;
      } else if (nblocks > (1 << d_dim.getValue())) {
         enhanced = true;
      }

      if (reduced || enhanced) {
         std::set<BlockId> sing_set;
         for (std::set<int>::iterator s_itr = (**p_itr).d_blocks.begin();
              s_itr != (**p_itr).d_blocks.end(); ++s_itr) {
            BoxId box_id(LocalId(*s_itr), 0);
            const BlockId& block_id =
               domain_boxes.find(Box(d_dim, box_id))->getBlockId();
            sing_set.insert(block_id);
         }
         singularity_blocks.insert(sing_set);
      }
   }

}

/*
 * ************************************************************************
 *
 * Figure out the connection between Box A and Box B
 *
 * ************************************************************************
 */

void
SingularityFinder::connect(const BoxId& id_a,
                           const BoxId& id_b,
                           int facea,
                           int faceb,
                           const BoxContainer& domain_boxes,
                           const BaseGridGeometry& grid_geometry)
{

   std::shared_ptr<Face> face = std::make_shared<Face>();
   d_faces.push_back(face);

   if (d_blocks.empty()) {
      d_blocks.resize(domain_boxes.size());
   }

   int a = id_a.getLocalId().getValue();
   int b = id_b.getLocalId().getValue();
   if (d_blocks[a].get() == nullptr) {
      d_blocks[a] = std::make_shared<Block>(d_dim);
   }
   if (d_blocks[b].get() == nullptr) {
      d_blocks[b] = std::make_shared<Block>(d_dim);
   }

   d_blocks[a]->d_face[facea] = face;
   d_blocks[b]->d_face[faceb] = face;
   face->d_block_to_face[a] = facea;
   face->d_block_to_face[b] = faceb;

   /*
    * Map from edge on 'a' to edge on 'b'
    */
   std::map<int, int> map_of_edges;
   if (d_dim.getValue() == 3) {
      findCoincidentEdges(map_of_edges,
         id_a,
         id_b,
         facea,
         faceb,
         domain_boxes,
         grid_geometry);
   }

   for (std::map<int, int>::const_iterator e_itr = map_of_edges.begin();
        e_itr != map_of_edges.end(); ++e_itr) {

      std::shared_ptr<Edge>& edgea = d_blocks[a]->d_edge[e_itr->first];
      std::shared_ptr<Edge>& edgeb = d_blocks[b]->d_edge[e_itr->second];

      if (edgea.get() == nullptr && edgeb.get() == nullptr) {
         edgea.reset(new Edge());
         edgeb = edgea;
         d_edges.push_back(edgea);
      } else if (edgea.get() != nullptr && edgeb.get() == nullptr) {
         edgeb = edgea;
      } else if (edgeb.get() != nullptr && edgea.get() == nullptr) {
         edgea = edgeb;
      } else if (edgea.get() == edgeb.get()) {
         // nothing needed
      } else {
         TBOX_ASSERT(edgea.get() != nullptr && edgeb.get() != nullptr);
         for (std::set<int>::iterator b_itr = edgeb->d_blocks.begin();
              b_itr != edgeb->d_blocks.end(); ++b_itr) {
            edgea->d_blocks.insert(*b_itr);
            edgea->d_block_to_edge[*b_itr] = edgeb->d_block_to_edge[*b_itr];
            if (*b_itr != a && *b_itr != b) {
               d_blocks[*b_itr]->d_edge[edgea->d_block_to_edge[*b_itr]] = edgea;
            }
         }
         for (std::list<std::shared_ptr<Edge> >::iterator e_itr =
                 d_edges.begin(); e_itr != d_edges.end(); ++e_itr) {
            if (e_itr->get() == edgeb.get()) {
               d_edges.erase(e_itr);
               break;
            }
         }
         edgeb = edgea;
      }

      edgea->d_blocks.insert(a);
      edgea->d_blocks.insert(b);
      edgea->d_block_to_edge[a] = e_itr->first;
      edgea->d_block_to_edge[b] = e_itr->second;
   }

   std::map<int, int> map_of_points;
   findCoincidentPoints(map_of_points,
      id_a,
      id_b,
      facea,
      domain_boxes,
      grid_geometry);

   for (std::map<int, int>::const_iterator e_itr = map_of_points.begin();
        e_itr != map_of_points.end(); ++e_itr) {

      std::shared_ptr<Point>& pointa = d_blocks[a]->d_point[e_itr->first];
      std::shared_ptr<Point>& pointb = d_blocks[b]->d_point[e_itr->second];

      if (pointa.get() == nullptr && pointb.get() == nullptr) {
         pointa.reset(new Point());
         pointb = pointa;
         d_points.push_back(pointa);
      } else if (pointa.get() != nullptr && pointb.get() == nullptr) {
         pointb = pointa;
      } else if (pointb.get() != nullptr && pointa.get() == nullptr) {
         pointa = pointb;
      } else if (pointa.get() == pointb.get()) {
         // nothing needed
      } else {
         TBOX_ASSERT(pointa.get() != nullptr && pointb.get() != nullptr);
         for (std::set<int>::iterator b_itr = pointb->d_blocks.begin();
              b_itr != pointb->d_blocks.end(); ++b_itr) {
            pointa->d_blocks.insert(*b_itr);
            pointa->d_block_to_point[*b_itr] =
               pointb->d_block_to_point[*b_itr];
            if (*b_itr != a && *b_itr != b) {
               d_blocks[*b_itr]->d_point[pointa->d_block_to_point[*b_itr]] =
                  pointa;
            }
         }
         for (std::list<std::shared_ptr<Point> >::iterator e_itr =
                 d_points.begin(); e_itr != d_points.end(); ++e_itr) {
            if (e_itr->get() == pointb.get()) {
               d_points.erase(e_itr);
               break;
            }
         }
         pointb = pointa;
      }

      pointa->d_blocks.insert(a);
      pointa->d_blocks.insert(b);
      pointa->d_block_to_point[a] = e_itr->first;
      pointa->d_block_to_point[b] = e_itr->second;
   }

}

/*
 * ************************************************************************
 *
 * Figure the boundary faces, those not shared between blocks or boxes
 *
 * ************************************************************************
 */

void
SingularityFinder::findBoundaryFaces()
{
   for (int iblock = 0; iblock < static_cast<int>(d_blocks.size()); ++iblock) {

      std::shared_ptr<Block>& block = d_blocks[iblock];

      for (int iface = 0; iface < 2 * d_dim.getValue(); ++iface) {

         std::shared_ptr<Face>& face = block->d_face[iface];

         if (face.get() == nullptr) {

            face.reset(new Face());
            d_faces.push_back(face);
            face->d_bdry = true;
            face->d_blocks.insert(iblock);

            if (d_dim.getValue() > 2) {
               for (int iedge = 0; iedge < 4; ++iedge) {
                  int edge_idx = s_face_edges[iface][iedge];
                  if (!block->d_edge[edge_idx]) {
                     std::shared_ptr<Edge> edge = std::make_shared<Edge>();
                     d_edges.push_back(edge);
                     edge->d_blocks.insert(iblock);
                     block->d_edge[edge_idx] = edge;
                  }
                  block->d_edge[edge_idx]->d_bdry = true;
               }
            }

            int num_pts = 1 << (d_dim.getValue() - 1);
            for (int ipoint = 0; ipoint < num_pts; ++ipoint) {
               int point_idx = s_face_nodes[iface][ipoint];
               if (!block->d_point[point_idx]) {
                  std::shared_ptr<Point> point = std::make_shared<Point>();
                  point->d_blocks.insert(iblock);
                  block->d_point[point_idx] = point;
               }
               block->d_point[point_idx]->d_bdry = true;
            }
         }
      }
   }
}

/*
 * ************************************************************************
 *
 * Given two boxes that are face neighbors, find their shared edges
 *
 * ************************************************************************
 */

void
SingularityFinder::findCoincidentEdges(
   std::map<int, int>& map_of_edges,
   const BoxId& id_a,
   const BoxId& id_b,
   int facea,
   int faceb,
   const BoxContainer& domain_boxes,
   const BaseGridGeometry& grid_geometry)
{
   if (d_dim.getValue() != 3) return;

   TBOX_ASSERT(map_of_edges.empty());

   Box a_box = *(domain_boxes.find(Box(d_dim, id_a)));
   Box b_box = *(domain_boxes.find(Box(d_dim, id_b)));

   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (a_box.lower() (d) == a_box.upper() (d)) {
         a_box.setLower(static_cast<Box::dir_t>(d),
            a_box.lower(static_cast<Box::dir_t>(d)) - 1);
         a_box.setUpper(static_cast<Box::dir_t>(d),
            a_box.upper(static_cast<Box::dir_t>(d)) + 1);
      }
      if (b_box.lower() (d) == b_box.upper() (d)) {
         b_box.setLower(static_cast<Box::dir_t>(d),
            b_box.lower(static_cast<Box::dir_t>(d)) - 1);
         b_box.setUpper(static_cast<Box::dir_t>(d),
            b_box.upper(static_cast<Box::dir_t>(d)) + 1);
      }
   }

   Box b_node_box(b_box);
   b_node_box.setUpper(b_node_box.upper() + IntVector::getOne(d_dim));
   IntVector b_box_size(b_node_box.numberCells());
   BoxContainer b_edge_boxes;

   int nedges_per_face = 4;

   for (int i = 0; i < nedges_per_face; ++i) {

      int edgea_idx = s_face_edges[facea][i];
      int edgeb_idx = -1;

      Box edge_box(a_box);

      switch (edgea_idx) {

         case 0:
            edge_box.setUpper(0, edge_box.lower(0));
            edge_box.setUpper(1, edge_box.lower(1));
            break;
         case 1:
            edge_box.setLower(0, edge_box.upper(0));
            edge_box.setUpper(1, edge_box.lower(1));
            break;
         case 2:
            edge_box.setUpper(0, edge_box.lower(0));
            edge_box.setLower(1, edge_box.upper(1));
            break;
         case 3:
            edge_box.setLower(0, edge_box.upper(0));
            edge_box.setLower(1, edge_box.upper(1));
            break;
         case 4:
            edge_box.setUpper(0, edge_box.lower(0));
            edge_box.setUpper(2, edge_box.lower(2));
            break;
         case 5:
            edge_box.setLower(0, edge_box.upper(0));
            edge_box.setUpper(2, edge_box.lower(2));
            break;
         case 6:
            edge_box.setUpper(0, edge_box.lower(0));
            edge_box.setLower(2, edge_box.upper(2));
            break;
         case 7:
            edge_box.setLower(0, edge_box.upper(0));
            edge_box.setLower(2, edge_box.upper(2));
            break;
         case 8:
            edge_box.setUpper(1, edge_box.lower(1));
            edge_box.setUpper(2, edge_box.lower(2));
            break;
         case 9:
            edge_box.setLower(1, edge_box.upper(1));
            edge_box.setUpper(2, edge_box.lower(2));
            break;
         case 10:
            edge_box.setUpper(1, edge_box.lower(1));
            edge_box.setLower(2, edge_box.upper(2));
            break;
         case 11:
            edge_box.setLower(1, edge_box.upper(1));
            edge_box.setLower(2, edge_box.upper(2));
            break;
         default:
            break;
      }

      if (a_box.getBlockId() != b_box.getBlockId()) {
         bool transformed = grid_geometry.transformBox(edge_box,
                                                       0,
                                                       b_box.getBlockId(),
                                                       a_box.getBlockId());
#ifndef DEBUG_CHECK_ASSERTIONS
         NULL_USE(transformed);
#endif
         TBOX_ASSERT(transformed);
      }
      edge_box.setUpper(edge_box.upper() + IntVector::getOne(d_dim));
      Box b_edge(edge_box * b_node_box);

      IntVector b_edge_dirs(b_edge.numberCells());
      int num_zero_dirs = 0;
      for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
         TBOX_ASSERT(b_edge_dirs[d] >= 1);
         if (b_edge_dirs[d] == b_box_size[d]) {
            b_edge_dirs[d] = 0;
            ++num_zero_dirs;
         } else if (b_edge.lower() (d) == b_node_box.lower(d)) {
            b_edge_dirs[d] = -1;
         } else if (b_edge.upper() (d) == b_node_box.upper(d)) {
            b_edge_dirs[d] = 1;
         } else {
            TBOX_ERROR(
               "SingularityFinder::findCoincidentEdges: Transformed box is not at the edge of the reference box.");
         }
      }

      if (num_zero_dirs != 1) {
         if (b_edge_boxes.empty()) {
            for (int e = 0; e < nedges_per_face; ++e) {

               int edge_idx = s_face_edges[faceb][e];

               Box add_box(b_box);
               add_box.setUpper(add_box.upper() + IntVector::getOne(d_dim));

               switch (edge_idx) {

                  case 0:
                     add_box.setUpper(0, add_box.lower(0));
                     add_box.setUpper(1, add_box.lower(1));
                     break;
                  case 1:
                     add_box.setLower(0, add_box.upper(0));
                     add_box.setUpper(1, add_box.lower(1));
                     break;
                  case 2:
                     add_box.setUpper(0, add_box.lower(0));
                     add_box.setLower(1, add_box.upper(1));
                     break;
                  case 3:
                     add_box.setLower(0, add_box.upper(0));
                     add_box.setLower(1, add_box.upper(1));
                     break;
                  case 4:
                     add_box.setUpper(0, add_box.lower(0));
                     add_box.setUpper(2, add_box.lower(2));
                     break;
                  case 5:
                     add_box.setLower(0, add_box.upper(0));
                     add_box.setUpper(2, add_box.lower(2));
                     break;
                  case 6:
                     add_box.setUpper(0, add_box.lower(0));
                     add_box.setLower(2, add_box.upper(2));
                     break;
                  case 7:
                     add_box.setLower(0, add_box.upper(0));
                     add_box.setLower(2, add_box.upper(2));
                     break;
                  case 8:
                     add_box.setUpper(1, add_box.lower(1));
                     add_box.setUpper(2, add_box.lower(2));
                     break;
                  case 9:
                     add_box.setLower(1, add_box.upper(1));
                     add_box.setUpper(2, add_box.lower(2));
                     break;
                  case 10:
                     add_box.setUpper(1, add_box.lower(1));
                     add_box.setLower(2, add_box.upper(2));
                     break;
                  case 11:
                     add_box.setLower(1, add_box.upper(1));
                     add_box.setLower(2, add_box.upper(2));
                     break;
                  default:
                     break;

               }

               b_edge_boxes.pushBack(add_box);
            }
         }

         BoxContainer b_edge_cntnr(b_edge);
         b_edge_cntnr.intersectBoxes(b_edge_boxes);
         b_edge_cntnr.coalesce();
         TBOX_ASSERT(b_edge_cntnr.size() == 1);

         b_edge = *(b_edge_cntnr.begin());
         b_edge_dirs = b_edge.numberCells();
         num_zero_dirs = 0;
         for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
            TBOX_ASSERT(b_edge_dirs[d] >= 1);
            if (b_edge_dirs[d] == b_box_size[d]) {
               b_edge_dirs[d] = 0;
               ++num_zero_dirs;
            } else if (b_edge.lower() (d) == b_node_box.lower(d)) {
               b_edge_dirs[d] = -1;
            } else if (b_edge.upper() (d) == b_node_box.upper(d)) {
               b_edge_dirs[d] = 1;
            } else {
               TBOX_ERROR(
                  "SingularityFinder::findCoincidentEdges: Transformed box is not at the edge of the reference box.");
            }
         }

         TBOX_ASSERT(num_zero_dirs == 1);
      }

      if (b_edge_dirs[0] == 0) {
         TBOX_ASSERT(b_edge_dirs[1] != 0 && b_edge_dirs[2] != 0);
         if (b_edge_dirs[1] == -1) {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 8;
            } else {
               edgeb_idx = 10;
            }
         } else {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 9;
            } else {
               edgeb_idx = 11;
            }
         }
      } else if (b_edge_dirs[1] == 0) {
         TBOX_ASSERT(b_edge_dirs[0] != 0 && b_edge_dirs[2] != 0);
         if (b_edge_dirs[0] == -1) {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 4;
            } else {
               edgeb_idx = 6;
            }
         } else {
            if (b_edge_dirs[2] == -1) {
               edgeb_idx = 5;
            } else {
               edgeb_idx = 7;
            }
         }
      } else if (b_edge_dirs[2] == 0) {
         TBOX_ASSERT(b_edge_dirs[0] != 0 && b_edge_dirs[1] != 0);
         if (b_edge_dirs[0] == -1) {
            if (b_edge_dirs[1] == -1) {
               edgeb_idx = 0;
            } else {
               edgeb_idx = 2;
            }
         } else {
            if (b_edge_dirs[1] == -1) {
               edgeb_idx = 1;
            } else {
               edgeb_idx = 3;
            }
         }
      } else {
         TBOX_ERROR("SingularityFinder::findCoincidentEdges failed to find location of edges.");
      }

      TBOX_ASSERT(edgeb_idx >= 0);

      map_of_edges[edgea_idx] = edgeb_idx;

   }

}

/*
 * ************************************************************************
 *
 * Given two boxes that are face neighbors, find their shared corner points
 *
 * ************************************************************************
 */

void
SingularityFinder::findCoincidentPoints(
   std::map<int, int>& map_of_points,
   const BoxId& id_a,
   const BoxId& id_b,
   int facea,
   const BoxContainer& domain_boxes,
   const BaseGridGeometry& grid_geometry)
{

   Box a_box = *(domain_boxes.find(Box(d_dim, id_a)));
   Box b_box = *(domain_boxes.find(Box(d_dim, id_b)));

   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (a_box.lower() (d) == a_box.upper() (d)) {
         a_box.setLower(static_cast<Box::dir_t>(d),
            a_box.lower(static_cast<Box::dir_t>(d)) - 1);
         a_box.setUpper(static_cast<Box::dir_t>(d),
            a_box.upper(static_cast<Box::dir_t>(d)) + 1);
      }
      if (b_box.lower() (d) == b_box.upper() (d)) {
         b_box.setLower(static_cast<Box::dir_t>(d),
            b_box.lower(static_cast<Box::dir_t>(d)) - 1);
         b_box.setUpper(static_cast<Box::dir_t>(d),
            b_box.upper(static_cast<Box::dir_t>(d)) + 1);
      }
   }

   Box b_node_box(b_box);
   b_node_box.setUpper(b_node_box.upper() + IntVector::getOne(d_dim));
   IntVector b_box_size(b_node_box.numberCells());

   int npoints_per_face = 1 << (d_dim.getValue() - 1);

   for (int i = 0; i < npoints_per_face; ++i) {

      int pointa_idx = s_face_nodes[facea][i];

      Box point_box(a_box);

      IntVector corner_dirs(d_dim, 0);
      corner_dirs[0] = (pointa_idx % 2 == 0) ? -1 : 1;
      if (d_dim.getValue() > 1) {
         corner_dirs[1] = ((pointa_idx / 2) % 2 == 0) ? -1 : 1;
      }
      if (d_dim.getValue() > 2) {
         corner_dirs[2] = ((pointa_idx / 4) % 2 == 0) ? -1 : 1;
      }

      for (tbox::Dimension::dir_t d = 0; d < d_dim.getValue(); ++d) {
         TBOX_ASSERT(corner_dirs[d] == -1 || corner_dirs[d] == 1);
         if (corner_dirs[d] == -1) {
            point_box.setUpper(d, point_box.lower(d));
         } else {
            point_box.setLower(d, point_box.upper(d));
         }
      }

      if (a_box.getBlockId() != b_box.getBlockId()) {
         bool transformed = grid_geometry.transformBox(point_box,
                                                       0,
                                                       b_box.getBlockId(),
                                                       a_box.getBlockId());
#ifndef DEBUG_CHECK_ASSERTIONS
         NULL_USE(transformed);
#endif
         TBOX_ASSERT(transformed);
      }

      point_box.setUpper(point_box.upper() + IntVector::getOne(d_dim));
      Box b_point(point_box * b_node_box);

      IntVector b_point_dirs(d_dim, 0);
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (b_point.lower() (d) == b_node_box.lower() (d)) {
            b_point_dirs[d] = -1;
         } else if (b_point.upper() (d) == b_node_box.upper() (d)) {
            b_point_dirs[d] = 1;
         } else {
            TBOX_ERROR(
               "SingularityFinder::findCoincidentPoints: Transformed box is not located at the corner of the reference box");
         }
      }

      int pointb_idx = 0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (b_point_dirs[d] == 1) {
            pointb_idx += (1 << d);
         }
      }

      TBOX_ASSERT(pointb_idx >= 0);

      map_of_points[pointa_idx] = pointb_idx;

   }

}

}
}
