/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for time interpolation
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   /*
    * Initialize tbox::MPI and SAMRAI, enable logging, and process command line.
    * Note this example is set up to run in serial only.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);

   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {

      if (argc != 2) {
         TBOX_ERROR(
            "USAGE:  " << argv[0] << " <input filename> "
                       << "<restart dir> <restore number> [options]\n"
                       << "  options:\n"
                       << "  none at this time"
                       << std::endl);
         return -1;
      }

      std::string input_filename = argv[1];

      /*
       * Create input database and parse all data in input file.
       */

      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Read "Main" input data.
       */

      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      int data_depth = main_db->getIntegerWithDefault("depth", 1);
      int ghosts = main_db->getIntegerWithDefault("ghosts", 4);
      double frac = main_db->getDoubleWithDefault(
        "interp_fraction", 0.31); 
      TBOX_ASSERT(data_depth >= 1);
      TBOX_ASSERT(ghosts >= 0);
      TBOX_ASSERT(frac >= 0.0 && frac <= 1.0);

      bool bdry_only = main_db->getBoolWithDefault("bdry_only", false);

      std::string log_file_name = "time_interp.log";
      if (main_db->keyExists("log_file_name")) {
         log_file_name = main_db->getString("log_file_name");
      }
      tbox::PIO::logOnlyNodeZero(log_file_name);

      double test_eps = sqrt(tbox::MathUtilities<double>::getEpsilon());

      hier::Transformation zero_transformation(
         hier::Transformation::NO_ROTATE,
         hier::IntVector::getZero(dim),
         hier::BlockId(0),
         hier::BlockId(0));

      hier::Box box(main_db->getDatabaseBox("box"));
      box.setBlockId(hier::BlockId(0));

      hier::IntVector ghost_vec(dim, ghosts);
      pdat::CellData<double> cell_old(box, data_depth, ghost_vec); 
      pdat::CellData<double> cell_new(box, data_depth, ghost_vec); 
      pdat::CellData<double> cell_dst(box, data_depth, ghost_vec); 
      pdat::CellData<double> cell_expected(box, data_depth, ghost_vec); 

      cell_old.setTime(0.0);
      cell_new.setTime(1.0);
      cell_dst.setTime(frac);
      cell_expected.setTime(frac);

      const hier::Box& ghost_box(cell_old.getGhostBox());
      auto cell_end(pdat::CellGeometry::end(ghost_box));
      for (int dd = 0; dd < data_depth; ++dd) { 
         for (auto cell_itr(pdat::CellGeometry::begin(ghost_box));
              cell_itr != cell_end; ++cell_itr) {
            double old_val = 3.7 * static_cast<double>((*cell_itr)[0]) + dd;
            double new_val = -1.5 * static_cast<double>((*cell_itr)[0]) + dd; 
            if (dim.getValue() > 1) {
               old_val += -4.4 * static_cast<double>((*cell_itr)[1]) + 2*dd;
               new_val += 10.3 * static_cast<double>((*cell_itr)[1]) + 2*dd;
            }
            if (dim.getValue() > 2) {
               old_val += -2.7 * static_cast<double>((*cell_itr)[2]) + 3*dd;
               new_val += 2.1 * static_cast<double>((*cell_itr)[2]) + 3*dd;
            }
            cell_old(*cell_itr, dd) = old_val;
            cell_new(*cell_itr, dd) = new_val;
            cell_expected(*cell_itr, dd) = old_val + frac*(new_val-old_val);
         }
      }

      hier::BoxContainer ghost_cntnr(ghost_box);
      pdat::CellOverlap cell_ovlp(ghost_cntnr, zero_transformation);
      pdat::CellDoubleLinearTimeInterpolateOp cell_op;
      cell_op.timeInterpolate(cell_dst, ghost_box, cell_ovlp, cell_old, cell_new);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      for (int dd = 0; dd < data_depth; ++dd) {
         for (auto cell_itr(pdat::CellGeometry::begin(ghost_box));
              cell_itr != cell_end; ++cell_itr) {
            double correct = cell_expected(*cell_itr, dd);
            double result = cell_dst(*cell_itr, dd);
            if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
               if (tbox::MathUtilities<double>::Abs(correct) < test_eps && 
                   tbox::MathUtilities<double>::Abs(result) < test_eps) {
                  continue;
               }
               tbox::perr << "Cell time interp test FAILED: ...."
                          << " : cell index = " << *cell_itr
                          << " of box"
                          << " " << ghost_box << std::endl;
               tbox::perr << "    result = " << result
                          << " : correct = " << correct << std::endl;
               ++fail_count; 
            }
         }
      }
#if 1
      pdat::NodeData<double> node_old(box, data_depth, ghost_vec);
      pdat::NodeData<double> node_new(box, data_depth, ghost_vec);
      pdat::NodeData<double> node_dst(box, data_depth, ghost_vec);
      pdat::NodeData<double> node_expected(box, data_depth, ghost_vec);

      node_old.setTime(0.0);
      node_new.setTime(1.0);
      node_dst.setTime(frac);
      node_expected.setTime(frac);

      auto node_end(pdat::NodeGeometry::end(ghost_box));
      for (int dd = 0; dd < data_depth; ++dd) {
         for (auto node_itr(pdat::NodeGeometry::begin(ghost_box));
              node_itr != node_end; ++node_itr) {
            double old_val = -7.5 * static_cast<double>((*node_itr)[0]) + dd;
            double new_val = -9.1 * static_cast<double>((*node_itr)[0]) + dd;
            if (dim.getValue() > 1) {
               old_val += -9.9 * static_cast<double>((*node_itr)[1]) + 2*dd;
               new_val += 8.5 * static_cast<double>((*node_itr)[1]) + 2*dd;
            }
            if (dim.getValue() > 2) {
               old_val += -6.4 * static_cast<double>((*node_itr)[2]) + 3*dd;
               new_val += -3.8 * static_cast<double>((*node_itr)[2]) + 3*dd;
            }
            node_old(*node_itr, dd) = old_val;
            node_new(*node_itr, dd) = new_val;
            node_expected(*node_itr, dd) = old_val + frac*(new_val-old_val);
         }
      }

      ghost_cntnr.clear();
      if (!bdry_only) {
         ghost_cntnr.push_back(pdat::NodeGeometry::toNodeBox(ghost_box));
      } else {
         hier::Box flat_box(pdat::NodeGeometry::toNodeBox(box));
         flat_box.setUpper(0, flat_box.lower(0));
         ghost_cntnr.push_back(flat_box);
      }
      pdat::NodeOverlap node_ovlp(ghost_cntnr, zero_transformation);
      pdat::NodeDoubleLinearTimeInterpolateOp node_op;

      node_op.timeInterpolate(node_dst, ghost_box, node_ovlp, node_old, node_new);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      hier::BoxContainer ovlp_boxes;
      node_ovlp.getSourceBoxContainer(ovlp_boxes);
      for (auto bitr = ovlp_boxes.begin();
           bitr != ovlp_boxes.end(); ++bitr) {
         for (int dd = 0; dd < data_depth; ++dd) {
            for (auto node_itr(pdat::NodeGeometry::begin(ghost_box));
                 node_itr != node_end; ++node_itr) {
               if ((*bitr).contains(*node_itr)) {
                  double correct = node_expected(*node_itr, dd);
                  double result = node_dst(*node_itr, dd);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                         tbox::MathUtilities<double>::Abs(result) < test_eps) {
                        continue;
                     }
                     tbox::perr << "Node time interp test FAILED: ...."
                                << " : node index = " << *node_itr
                                << " of box"
                                << " " << ghost_box << std::endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << std::endl;
                     ++fail_count;
                  }
               }
            }
         }
      }
#endif
      pdat::FaceData<double> face_old(box, data_depth, ghost_vec);
      pdat::FaceData<double> face_new(box, data_depth, ghost_vec);
      pdat::FaceData<double> face_dst(box, data_depth, ghost_vec);
      pdat::FaceData<double> face_expected(box, data_depth, ghost_vec);

      face_old.setTime(0.0);
      face_new.setTime(1.0);
      face_dst.setTime(frac);
      face_expected.setTime(frac);

      std::vector<hier::BoxContainer> ghost_bxs(dim.getValue());
      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         auto face_end(pdat::FaceGeometry::end(ghost_box, axis));
         for (int dd = 0; dd < data_depth; ++dd) {
            for (auto face_itr(pdat::FaceGeometry::begin(ghost_box, axis));
                 face_itr != face_end; ++face_itr) {
               double old_val = -8.1 * static_cast<double>((*face_itr)[0]) + dd;
               double new_val = -7.3 * static_cast<double>((*face_itr)[0]) + dd;
               if (dim.getValue() > 1) {
                  old_val += 6.9 * static_cast<double>((*face_itr)[1]) + 2*dd;
                  new_val += 8.5 * static_cast<double>((*face_itr)[1]) + 2*dd;
               }
               if (dim.getValue() > 2) {
                  old_val += 2.4 * static_cast<double>((*face_itr)[2]) + 3*dd;
                  new_val += -4.7 * static_cast<double>((*face_itr)[2]) + 3*dd;
               }
               old_val += static_cast<double>(axis);
               new_val -= static_cast<double>(axis);
               face_old(*face_itr, dd) = old_val;
               face_new(*face_itr, dd) = new_val;
               face_expected(*face_itr, dd) = old_val + frac*(new_val-old_val);
            }
         }
         if (!bdry_only) {
            ghost_bxs[axis].push_back(
               pdat::FaceGeometry::toFaceBox(ghost_box, axis));
         } else {
            hier::Box flat_box(pdat::FaceGeometry::toFaceBox(box, axis));
            flat_box.setUpper(axis, flat_box.lower(axis));
            ghost_bxs[axis].push_back(flat_box); 
         }
      }

      pdat::FaceOverlap face_ovlp(ghost_bxs, zero_transformation);
      pdat::FaceDoubleLinearTimeInterpolateOp face_op;
      face_op.timeInterpolate(face_dst, ghost_box, face_ovlp, face_old, face_new);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         hier::BoxContainer ovlp_boxes;
         int normal = axis;
         face_ovlp.getSourceBoxContainer(ovlp_boxes, normal);
         for (auto bitr = ovlp_boxes.begin();
              bitr != ovlp_boxes.end(); ++bitr) {
            auto face_end(pdat::FaceGeometry::end(ghost_box, axis));
            for (int dd = 0; dd < data_depth; ++dd) {
               for (auto face_itr(pdat::FaceGeometry::begin(ghost_box, axis));
                    face_itr != face_end; ++face_itr) {
                  if ((*bitr).contains(*face_itr)) { 
                     double correct = face_expected(*face_itr, dd);
                     double result = face_dst(*face_itr, dd);
                     if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                        if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                            tbox::MathUtilities<double>::Abs(result) < test_eps) {
                           continue;
                        }
                        tbox::perr << "Face time interp test FAILED: ...."
                                   << " : face index = " << *face_itr
                                   << " of box"
                                   << " " << ghost_box << std::endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << std::endl;
                        ++fail_count;
                     }
                  }
               }
            }
         }
      }

      pdat::OuterfaceData<double> oface_old(box, data_depth);
      pdat::OuterfaceData<double> oface_new(box, data_depth);
      pdat::OuterfaceData<double> oface_dst(box, data_depth);

      oface_old.setTime(0.0);
      oface_new.setTime(1.0);
      oface_dst.setTime(frac);

      oface_old.copy(face_old);
      oface_new.copy(face_new);

      pdat::OuterfaceDoubleLinearTimeInterpolateOp oface_op;

      oface_op.timeInterpolate(oface_dst, box, face_ovlp, oface_old, oface_new);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         for (int face = 0; face < 2; ++face) {
            const hier::Box& databox =
               oface_dst.getArrayData(axis, face).getBox();
            auto face_end(pdat::FaceGeometry::end(box, axis));
            for (int dd = 0; dd < data_depth; ++dd) {
               for (auto face_itr(pdat::FaceGeometry::begin(box, axis));
                    face_itr != face_end; ++face_itr) {
                  if (!databox.contains(*face_itr)) {
                     continue;
                  }
                  double correct = face_expected(*face_itr, dd);
                  double result = oface_dst(*face_itr, face, dd);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                         tbox::MathUtilities<double>::Abs(result) < test_eps) {
                        continue;
                     }
                     tbox::perr << "Outerface time interp test FAILED: ...."
                                << " : face index = " << *face_itr
                                << " of box"
                                << " " << ghost_box << std::endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << std::endl;
                     ++fail_count;
                  }
               }
            }  
         }  
      }

      pdat::SideData<double> side_old(box, data_depth, ghost_vec);
      pdat::SideData<double> side_new(box, data_depth, ghost_vec);
      pdat::SideData<double> side_dst(box, data_depth, ghost_vec);
      pdat::SideData<double> side_expected(box, data_depth, ghost_vec);

      side_old.setTime(0.0);
      side_new.setTime(1.0);
      side_dst.setTime(frac);
      side_expected.setTime(frac);

      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         auto side_end(pdat::SideGeometry::end(ghost_box, axis));
         for (int dd = 0; dd < data_depth; ++dd) {
            for (auto side_itr(pdat::SideGeometry::begin(ghost_box, axis));
                 side_itr != side_end; ++side_itr) {
               double old_val = -8.3 * static_cast<double>((*side_itr)[0]) + dd;
               double new_val = -9.2 * static_cast<double>((*side_itr)[0]) + dd;
               if (dim.getValue() > 1) {
                  old_val += 7.6 * static_cast<double>((*side_itr)[1]) + 2*dd;
                  new_val += 4.7 * static_cast<double>((*side_itr)[1]) + 2*dd;
               }
               if (dim.getValue() > 2) {
                  old_val += 5.8 * static_cast<double>((*side_itr)[2]) + 3*dd;
                  new_val += -5.4 * static_cast<double>((*side_itr)[2]) + 3*dd;
               }
               old_val += static_cast<double>(axis);
               new_val -= static_cast<double>(axis);
               side_old(*side_itr, dd) = old_val;
               side_new(*side_itr, dd) = new_val;
               side_expected(*side_itr, dd) = old_val + frac*(new_val-old_val);
            }
         }
         ghost_bxs[axis].clear();
         if (!bdry_only) {
            ghost_bxs[axis].push_back(
               pdat::SideGeometry::toSideBox(ghost_box, axis));
         } else {
            hier::Box flat_box(pdat::SideGeometry::toSideBox(box, axis));
            flat_box.setUpper(axis, flat_box.lower(axis));
            ghost_bxs[axis].push_back(flat_box);
         }
      }

      pdat::SideOverlap side_ovlp(ghost_bxs, zero_transformation);
      pdat::SideDoubleLinearTimeInterpolateOp side_op;

      side_op.timeInterpolate(side_dst, ghost_box, side_ovlp, side_old, side_new);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         hier::BoxContainer ovlp_boxes;
         int normal = axis;
         side_ovlp.getSourceBoxContainer(ovlp_boxes, normal);
         for (auto bitr = ovlp_boxes.begin();
              bitr != ovlp_boxes.end(); ++bitr) {
            auto side_end(pdat::SideGeometry::end(ghost_box, axis));
            for (int dd = 0; dd < data_depth; ++dd) {
               for (auto side_itr(pdat::SideGeometry::begin(ghost_box, axis));
                    side_itr != side_end; ++side_itr) {
                  if ((*bitr).contains(*side_itr)) {
                     double correct = side_expected(*side_itr, dd);
                     double result = side_dst(*side_itr, dd);
                     if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                        if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                            tbox::MathUtilities<double>::Abs(result) < test_eps) {
                           continue;
                        }
                        tbox::perr << "Side time interp test FAILED: ...."
                                 << " : side index = " << *side_itr
                                 << " of box"
                                 << " " << ghost_box << std::endl;
                        tbox::perr << "    result = " << result
                                 << " : correct = " << correct << std::endl;
                        ++fail_count;
                     }
                  }
               }
            }
         }
      }

      pdat::OutersideData<double> oside_old(box, data_depth);
      pdat::OutersideData<double> oside_new(box, data_depth);
      pdat::OutersideData<double> oside_dst(box, data_depth);

      oside_old.setTime(0.0);
      oside_new.setTime(1.0);
      oside_dst.setTime(frac);

      oside_old.copy(side_old);
      oside_new.copy(side_new);

      pdat::OutersideDoubleLinearTimeInterpolateOp oside_op;

      oside_op.timeInterpolate(oside_dst, box, side_ovlp, oside_old, oside_new);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         for (int side = 0; side < 2; ++side) {
            const hier::Box& databox =
               oside_dst.getArrayData(axis, side).getBox();
            auto side_end(pdat::SideGeometry::end(box, axis));
            for (int dd = 0; dd < data_depth; ++dd) {
               for (auto side_itr(pdat::SideGeometry::begin(box, axis));
                    side_itr != side_end; ++side_itr) {
                  if (!databox.contains(*side_itr)) {
                     continue;
                  }
                  double correct = side_expected(*side_itr, dd);
                  double result = oside_dst(*side_itr, side, dd);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                         tbox::MathUtilities<double>::Abs(result) < test_eps) {
                        continue;
                     }  
                     tbox::perr << "Outerside time interp test FAILED: ...."
                                << " : side index = " << *side_itr
                                << " of box"
                                << " " << ghost_box << std::endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << std::endl;
                     ++fail_count;
                  }
               }
            }
         }
      }

      pdat::EdgeData<double> edge_old(box, data_depth, ghost_vec);
      pdat::EdgeData<double> edge_new(box, data_depth, ghost_vec);
      pdat::EdgeData<double> edge_dst(box, data_depth, ghost_vec);
      pdat::EdgeData<double> edge_expected(box, data_depth, ghost_vec);

      edge_old.setTime(0.0);
      edge_new.setTime(1.0);
      edge_dst.setTime(frac);
      edge_expected.setTime(frac);

      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         auto edge_end(pdat::EdgeGeometry::end(ghost_box, axis));
         for (int dd = 0; dd < data_depth; ++dd) {
            for (auto edge_itr(pdat::EdgeGeometry::begin(ghost_box, axis));
                 edge_itr != edge_end; ++edge_itr) {
               double old_val = -6.0 * static_cast<double>((*edge_itr)[0]) + dd;
               double new_val = 3.5 * static_cast<double>((*edge_itr)[0]) + dd;
               if (dim.getValue() > 1) {
                  old_val += -8.9 * static_cast<double>((*edge_itr)[1]) + 2*dd;
                  new_val += -9.5 * static_cast<double>((*edge_itr)[1]) + 2*dd;
               }
               if (dim.getValue() > 2) {
                  old_val += -3.3 * static_cast<double>((*edge_itr)[2]) + 3*dd;
                  new_val += 6.6 * static_cast<double>((*edge_itr)[2]) + 3*dd;
               }
               old_val += static_cast<double>(axis);
               new_val -= static_cast<double>(axis);
               edge_old(*edge_itr, dd) = old_val;
               edge_new(*edge_itr, dd) = new_val;
               edge_expected(*edge_itr, dd) = old_val + frac*(new_val-old_val);
            }
         }
         ghost_bxs[axis].clear();
         if (!bdry_only) {
            ghost_bxs[axis].push_back(
               pdat::EdgeGeometry::toEdgeBox(ghost_box, axis));
         } else {
            hier::Box flat_box(pdat::EdgeGeometry::toEdgeBox(box, axis));
            flat_box.setUpper(axis, flat_box.lower(axis));
            ghost_bxs[axis].push_back(flat_box);
         }
      }

      pdat::EdgeOverlap edge_ovlp(ghost_bxs, zero_transformation);
      pdat::EdgeDoubleLinearTimeInterpolateOp edge_op;

      edge_op.timeInterpolate(edge_dst, ghost_box, edge_ovlp, edge_old, edge_new);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      for (unsigned short axis = 0; axis < dim.getValue(); ++axis) {
         hier::BoxContainer ovlp_boxes;
         int normal = axis;
         edge_ovlp.getSourceBoxContainer(ovlp_boxes, normal);
         for (auto bitr = ovlp_boxes.begin();
              bitr != ovlp_boxes.end(); ++bitr) {
            auto edge_end(pdat::EdgeGeometry::end(ghost_box, axis));
            for (int dd = 0; dd < data_depth; ++dd) {
               for (auto edge_itr(pdat::EdgeGeometry::begin(ghost_box, axis));
                    edge_itr != edge_end; ++edge_itr) {
                  if ((*bitr).contains(*edge_itr)) {
                     double correct = edge_expected(*edge_itr, dd);
                     double result = edge_dst(*edge_itr, dd);
                     if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                        if (tbox::MathUtilities<double>::Abs(correct) < test_eps &&
                            tbox::MathUtilities<double>::Abs(result) < test_eps) {
                           continue;
                        }
                        tbox::perr << "Edge time interp test FAILED: ...."
                                   << " : edge index = " << *edge_itr
                                   << " of box"
                                   << " " << ghost_box << std::endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << std::endl;
                        ++fail_count;
                     }
                  }
               }
            }
         }
      }
      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  time interp test" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}
