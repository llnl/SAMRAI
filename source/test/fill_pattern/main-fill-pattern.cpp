/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/FirstLayerCellVariableFillPattern.h"
#include "SAMRAI/pdat/FirstLayerCellNoCornersVariableFillPattern.h"
#include "SAMRAI/pdat/SecondLayerNodeVariableFillPattern.h"
#include "SAMRAI/pdat/SecondLayerNodeNoCornersVariableFillPattern.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <cstring>
#include <stdlib.h>

using namespace SAMRAI;

/*
 *
 * Test program to test VariableFillPattern implementations
 *
 */

void txt2boxes(
   const char* txt,
   hier::BoxContainer& boxes)
{
   // algorithm:
   // find width
   // find height
   // find x locations in i,j
   // foreach x1 in x:
   //   find x' where i>i1, j>j1
   //     foreach x2 in x':
   //       if 4 corners && interior cells blank, (x1,x2) is a box
   // translate coordinates into cell-centered, lower-left origin based

   int width = -1;
   for (unsigned int idx = 0; idx < strlen(txt) - 1; ++idx) {
      if (('x' == txt[idx] || '.' == txt[idx]) &&
          ('.' == txt[idx + 1] || '|' == txt[idx + 1])) {
         width = idx + 1;
         break;
      }
   }
   if (-1 == width) {
      std::cout << "error in box txt" << std::endl;
      exit(1);
   }

   int height = static_cast<int>(strlen(txt)) / width;

   // Find cell height
   int cell_height = (height - 1) / 2;
   int cell_max = cell_height - 1;

   // make vector of x locations
   std::vector<std::pair<int, int> > ix;
   for (unsigned int idx = 0; idx < strlen(txt); ++idx) {
      if ('x' == txt[idx]) {
         int j = idx / width;
         int i = idx - j * width;
         ix.push_back(std::pair<int, int>(i, j));
      }
   }

   // foreach x1 in x
   std::vector<std::pair<int, int> >::iterator it;
   for (it = ix.begin(); it != ix.end(); ++it) {

      std::vector<std::pair<int, int> >::iterator it2;

      // We need to gather all potential boxes rooted here, and then
      // only take the smallest one.

      std::vector<hier::Box> boxes_here;
      boxes_here.clear();

      for (it2 = ix.begin(); it2 != ix.end(); ++it2) {

         if (it2->first > it->first &&
             it2->second > it->second) {

            bool isbox = true;

            // If the two other corners exist, and...
            int i1 = it->first;
            int j1 = it2->second;
            int idx1 = j1 * width + i1;
            if (txt[idx1] != 'x') isbox = false;

            int i2 = it2->first;
            int j2 = it->second;
            int idx2 = j2 * width + i2;
            if (txt[idx2] != 'x') isbox = false;

            // ...interior cells contain no corners
            for (int i = it->first + 1; i < it2->first; ++i) {
               for (int j = it->second + 1; j < it2->second; ++j) {
                  int idx = j * width + i;
                  if ('x' == txt[idx]) isbox = false;
                  if ('-' == txt[idx]) isbox = false;
                  if ('|' == txt[idx]) isbox = false;
               }
            }

            if (isbox) {

               // Translate indices into node centered coords
               int i0 = it->first / 4;
               i1 = it2->first / 4;
               int j0 = it->second / 2;
               j1 = it2->second / 2;

               --i1;
               --j1;

               // Flip coordinates vertically.
               j0 = cell_max - j0;
               j1 = cell_max - j1;

               // Lower left uses j1, upper right j0
               int tmp = j1;
               j1 = j0;
               j0 = tmp;

               hier::Box abox(hier::Index(i0, j0),
                              hier::Index(i1, j1),
                              hier::BlockId(0));
               boxes_here.push_back(abox);
            }
         }
      }

      // Find smallest box at this 'x'
      if (boxes_here.size()) {

         hier::Box smallest_box(boxes_here[0]);

         for (std::vector<hier::Box>::iterator itb = boxes_here.begin();
              itb != boxes_here.end(); ++itb) {
            if ((*itb).numberCells() < smallest_box.numberCells()) {
               smallest_box = *itb;
            }
         }

         boxes.pushBack(smallest_box);
      }

   }

   // Shift all boxes into SAMRAI coordinates
   for (hier::BoxContainer::iterator itr = boxes.begin();
        itr != boxes.end(); ++itr) {
      itr->shift(-hier::IntVector(tbox::Dimension(2), 2));
   }
}

int txt_width(
   const char* txt)
{
   int width = -1;
   for (unsigned int idx = 0; idx < strlen(txt) - 1; ++idx) {
      if (('x' == txt[idx] || '.' == txt[idx]) &&
          ('.' == txt[idx + 1] || '|' == txt[idx + 1])) {
         width = idx + 1;
         break;
      }
   }
   if (-1 == width) {
      std::cout << "error in box txt" << std::endl;
      exit(1);
   }
   return width;
}

bool txt_next_val(
   const char* txt,
   int& idx,
   const hier::PatchData& data,
   int* datapt,
   bool is_node)
{
   // Find text size
   int txt_w = txt_width(txt);
   int txt_h = static_cast<int>(strlen(txt)) / txt_w;

   // Find grid size
   int grid_height = (txt_h - 1) / 2;
   int grid_width = (txt_w - 1) / 4;
   int grid_max = grid_height - 1;

   int cnt_max = 10000; // limit infinite loop possibility
   int cnt = 0;
   do {

      //
      // Translate domain local idx into grid idx
      //

      //const hier::Box& ghost_box = data.getGhostBox();
      hier::Box ghost_box(data.getGhostBox().getDim());
      if (is_node) {
         ghost_box = pdat::NodeGeometry::toNodeBox(data.getGhostBox());
      } else {
         ghost_box = data.getGhostBox();
      }
      // Translate domain idx to domain coordinates
      int domain_i = idx % ghost_box.numberCells(0);
      int domain_j = idx / ghost_box.numberCells(0);

      tbox::Dimension dim(ghost_box.getDim());
      hier::Box shifted_box(hier::Box::shift(ghost_box,
                               hier::IntVector(dim, 2)));
      // Translate domain coordinates into grid zone coordintes
      int di = shifted_box.lower() (0);
      int dj = shifted_box.lower() (1);

      int grid_i = domain_i + di;
      int grid_j = domain_j + dj;

      // If we outside the grid, there cannot be a value here
      if (grid_i < 0 || grid_j < 0) {
         ++idx;
         continue;
      }
      if (grid_i > grid_width || grid_j > grid_height) {
         ++idx;
         continue;
      }
      // Translate grid coords to text coordinates.  Text coordinates
      // have j increasing downwards.
      int txt_zone_i = grid_i * 4 + 2;
      int txt_zone_j = (grid_max - grid_j) * 2 + 1;

      int txt_node_i = grid_i * 4;
      int txt_node_j = (grid_max - grid_j) * 2 + 2;

      // Translate text coordinates to txt idx
      unsigned int txt_zone_idx = txt_zone_i + txt_zone_j * txt_w;
      int txt_node_idx = txt_node_i + txt_node_j * txt_w;

      // If we're past the end of the txt, return false
      if (txt_zone_idx > strlen(txt)) {
         return false;
      }

      // Check for non-zero zone data
      if (' ' != txt[txt_zone_idx]) {

         std::istringstream valstr(&txt[txt_zone_idx]);
         valstr >> *datapt;
         return true;
      }

      // Check for numeric node data
      if ('0' == txt[txt_node_idx] ||
          '1' == txt[txt_node_idx] ||
          '2' == txt[txt_node_idx] ||
          '3' == txt[txt_node_idx] ||
          '4' == txt[txt_node_idx] ||
          '5' == txt[txt_node_idx] ||
          '6' == txt[txt_node_idx] ||
          '7' == txt[txt_node_idx] ||
          '8' == txt[txt_node_idx] ||
          '9' == txt[txt_node_idx]) {

         std::istringstream valstr(&txt[txt_node_idx]);
         valstr >> *datapt;
         return true;
      }

      ++idx; // advance to next domain idx

   } while (cnt++ < cnt_max);

   std::cout << "Data reading loop exceeded maximum iterations"
        << __LINE__ << " in "
        << __FILE__ << std::endl;

   exit(1);
   return false;
}

void txt2data(
   const char* txt,
   const hier::PatchData& data,
   int* datptr,
   bool zero_out,
   bool is_node)
{
   if (zero_out) memset(datptr, 0, data.getGhostBox().size() * sizeof(int));

   int idx = 0;
   int datapt;

   while (txt_next_val(txt, idx, data, &datapt, is_node)) {
      datptr[idx++] = datapt;
   }
}

/*
 * Acceptance test cases.  First iteration at an executable
 * specification of the desired change.
 */

bool SingleLevelTestCase(
   const char* levelboxes_txt,
   const char* initialdata_txt[],
   const char* finaldata_txt[],
   std::shared_ptr<hier::Variable> variable,
   std::shared_ptr<xfer::VariableFillPattern> fill_pattern,
   tbox::Dimension& dim)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const std::string& pattern_name = fill_pattern->getPatternName();

   hier::BoxContainer level_boxes;
   txt2boxes(levelboxes_txt, level_boxes);

   std::shared_ptr<geom::GridGeometry> geom(
      new geom::GridGeometry(
         "GridGeometry",
         level_boxes));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
      new hier::PatchHierarchy("hier", geom));

   std::shared_ptr<hier::BoxLevel> mblevel(
      std::make_shared<hier::BoxLevel>(hier::IntVector(dim, 1), geom));

   const int num_nodes = mpi.getSize();
   const int num_boxes = level_boxes.size();
   hier::LocalId local_id(0);
   hier::BoxContainer::iterator level_boxes_itr = level_boxes.begin();
   for (int i = 0; i < num_boxes; ++i, ++level_boxes_itr) {

      int proc;
      if (i < num_boxes / num_nodes) {
         proc = 0;
      } else {
         proc = 1;
      }

      if (proc == mpi.getRank()) {
         mblevel->addBox(hier::Box(*level_boxes_itr, local_id, proc));
         ++local_id;
      }

   }

   int level_no = 0;
   hierarchy->makeNewPatchLevel(level_no, mblevel);

   std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(0));

   // There is one variable-context pair with a gcw of 2

   xfer::RefineAlgorithm refine_alg;

   std::shared_ptr<hier::VariableContext> context(
      hier::VariableDatabase::getDatabase()->getContext("CONTEXT"));

   hier::IntVector ghost_cell_width(dim, 2);

   int data_id =
      hier::VariableDatabase::getDatabase()->registerVariableAndContext(
         variable, context, ghost_cell_width);

   refine_alg.registerRefine(data_id, data_id, data_id,
      std::shared_ptr<hier::RefineOperator>(),
      fill_pattern);

   level->allocatePatchData(data_id);

   if (pattern_name == "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN" ||
       pattern_name == "FIRST_LAYER_CELL_FILL_PATTERN") {
      // Loop over each patch and initialize data
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch(*p);
         std::shared_ptr<pdat::CellData<int> > cdata(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
               patch->getPatchData(data_id)));
         TBOX_ASSERT(cdata);

         int data_txt_id = patch->getBox().getLocalId().getValue();
         if (mpi.getRank() == 1) {
            data_txt_id += (num_boxes / num_nodes);
         }

         txt2data(initialdata_txt[data_txt_id], *cdata,
            cdata->getPointer(), false, false);
      }
   } else if (pattern_name == "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN" ||
              pattern_name == "SECOND_LAYER_NODE_FILL_PATTERN") {
      // Loop over each patch and initialize data
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch(*p);
         std::shared_ptr<pdat::NodeData<int> > ndata(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<int>, hier::PatchData>(
               patch->getPatchData(data_id)));
         TBOX_ASSERT(ndata);

         int data_txt_id = patch->getBox().getLocalId().getValue();
         if (mpi.getRank() == 1) {
            data_txt_id += (num_boxes / num_nodes);
         }

         txt2data(initialdata_txt[data_txt_id], *ndata,
            ndata->getPointer(), false, true);
      }
   }

   // Cache Connector required for the schedule generation.
   level->findConnector(*level,
      hier::IntVector(dim, 2),
      hier::CONNECTOR_CREATE);

   // Create and run comm schedule
   refine_alg.createSchedule(level)->fillData(0.0, false);

   // Check for expected data
   bool failed = false;

   if (pattern_name == "FIRST_LAYER_CELL_NO_CORNERS_FILL_PATTERN" ||
       pattern_name == "FIRST_LAYER_CELL_FILL_PATTERN") {
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch(*p);
         std::shared_ptr<pdat::CellData<int> > cdata(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
               patch->getPatchData(data_id)));
         TBOX_ASSERT(cdata);

         pdat::CellData<int> expected(cdata->getBox(),
                                      cdata->getDepth(),
                                      ghost_cell_width);

         int data_txt_id = patch->getBox().getLocalId().getValue();
         if (mpi.getRank() == 1) {
            data_txt_id += (num_boxes / num_nodes);
         }

         txt2data(finaldata_txt[data_txt_id],
            expected, expected.getPointer(), false, false);

         pdat::CellData<int>::iterator ciend(pdat::CellGeometry::end(cdata->getGhostBox()));
         for (pdat::CellData<int>::iterator ci(pdat::CellGeometry::begin(cdata->getGhostBox()));
              ci != ciend; ++ci) {
            if ((*cdata)(*ci) != expected(*ci)) {
               failed = true;
            }
         }

      }
   } else if (pattern_name == "SECOND_LAYER_NODE_NO_CORNERS_FILL_PATTERN" ||
              pattern_name == "SECOND_LAYER_NODE_FILL_PATTERN") {
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch(*p);
         std::shared_ptr<pdat::NodeData<int> > ndata(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<int>, hier::PatchData>(
               patch->getPatchData(data_id)));
         TBOX_ASSERT(ndata);

         pdat::NodeData<int> expected(ndata->getBox(),
                                      ndata->getDepth(),
                                      ghost_cell_width);

         int data_txt_id = patch->getBox().getLocalId().getValue();
         if (mpi.getRank() == 1) {
            data_txt_id += (num_boxes / num_nodes);
         }

         txt2data(finaldata_txt[data_txt_id],
            expected, expected.getPointer(), false, true);

         pdat::NodeData<int>::iterator niend(pdat::NodeGeometry::end(ndata->getGhostBox()));
         for (pdat::NodeData<int>::iterator ni(pdat::NodeGeometry::begin(ndata->getGhostBox()));
              ni != niend; ++ni) {
            if ((*ndata)(*ni) != expected(*ni)) {
               failed = true;
            }
         }

      }
   }

   if (failed) {
      tbox::perr << "FAILED: - Test of " << pattern_name << std::endl;
   }

   return failed;
}

/*
 * This tests FirstLayerCellNoCornersVariableFillPattern ..
 */

bool Test_FirstLayerCellNoCornersVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 4
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 2
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 2 . 2 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 3 . 3 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 0 data after comm

   const char* final0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 0 . 0 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 1 . 1 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   const char* initial_txt[4] = { initial0_txt, initial1_txt,
                                  initial2_txt, initial3_txt };
   const char* final_txt[4] = { final0_txt, final1_txt,
                                final2_txt, final3_txt };
   tbox::Dimension dim(2);

   std::shared_ptr<pdat::CellVariable<int> > var(
      new pdat::CellVariable<int>(dim, "1cellnocorners"));

   std::shared_ptr<pdat::FirstLayerCellNoCornersVariableFillPattern> fill_pattern(
      new pdat::FirstLayerCellNoCornersVariableFillPattern(dim));

   return SingleLevelTestCase(levelboxes_txt,
      initial_txt,
      final_txt,
      var,
      fill_pattern,
      dim);
}

/*
 * This tests FirstLayerCellVariableFillPattern
 */

bool Test_FirstLayerCellVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 4
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 2
      ". . . . x . . . x . . . x . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   const char* final0_txt =
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 6
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 1 . 0 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 0 . 0 . 2 . 2 . 3 . 0 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 0 . 0 . 0 . 0 . 0 . 0 .   .   ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 6
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 0 . 1 . 1 . 1 . 1 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 1 . 2 . 3 . 3 . 1 . 1 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 1 . 1 . 1 . 1 . 1 . 1 ." // 2
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 5
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 0 . 0 . 1 . 2 .   .   ." // 4
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 3
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 3 . 2 .   .   ." // 2
      ". . . . x . . . x . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 1
      ". . . . . . . . . . . . . . . . ."
      ". 2 . 2 . 2 . 2 . 2 . 2 .   .   ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 7
      ". . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 5
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 0 . 1 . 1 . 3 . 3 ." // 4
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 3
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 2 . 3 . 3 . 3 . 3 ." // 2
      ". . . . . . . . x . . . x . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 1
      ". . . . . . . . . . . . . . . . ."
      ".   .   . 3 . 3 . 3 . 3 . 3 . 3 ." // 0
      ". . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7
   ;

   const char* initial_txt[4] = { initial0_txt, initial1_txt,
                                  initial2_txt, initial3_txt };
   const char* final_txt[4] = { final0_txt, final1_txt,
                                final2_txt, final3_txt };

   tbox::Dimension dim(2);

   std::shared_ptr<pdat::CellVariable<int> > var(
      new pdat::CellVariable<int>(dim, "1cell"));

   std::shared_ptr<pdat::FirstLayerCellVariableFillPattern> fill_pattern(
      new pdat::FirstLayerCellVariableFillPattern(dim));

   return SingleLevelTestCase(levelboxes_txt,
      initial_txt,
      final_txt,
      var,
      fill_pattern,
      dim);
}

bool Test_SecondLayerNodeNoCornersVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 8 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 8 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 1 . 8 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 1 . 1 . 8 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 8 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 8 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 3 . 3 . 8 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 3 . 8 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   const char* final0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 8 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 8 . 0 . 8 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 2 . 2 . 8 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 0 . 8 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 8 . 1 . 8 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 8 . 3 . 3 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 0 . 0 . 8 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 8 . 2 . 8 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 8 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 8 . 1 . 1 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 8 . 3 . 8 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 2 . 8 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   const char* initial_txt[4] = { initial0_txt, initial1_txt,
                                  initial2_txt, initial3_txt };
   const char* final_txt[4] = { final0_txt, final1_txt,
                                final2_txt, final3_txt };

   tbox::Dimension dim(2);

   std::shared_ptr<pdat::NodeVariable<int> > var(
      new pdat::NodeVariable<int>(
         dim,
         "secondnodenocorners"));

   std::shared_ptr<pdat::SecondLayerNodeNoCornersVariableFillPattern>
   fill_pattern(
      new pdat::SecondLayerNodeNoCornersVariableFillPattern(dim));

   return SingleLevelTestCase(levelboxes_txt,
      initial_txt,
      final_txt,
      var,
      fill_pattern,
      dim);
}

bool Test_SecondLayerNodeVariableFillPattern()
{
   const char* levelboxes_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . x . . . x . . . x . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 0 data before comm

   const char* initial0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 8 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 8 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 1 data before comm

   const char* initial1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 1 . 8 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 1 . 1 . 8 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 2 data before comm

   const char* initial2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 8 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 8 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // patch 3 data before comm

   const char* initial3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 3 . 3 . 8 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 3 . 8 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   const char* final0_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 0 . 0 . 0 . 0 . 0 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 0 . 0 . 0 . 0 . 8 . 1 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 0 . 0 . 0 . 8 . 0 . 8 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 0 . 0 . 2 . 2 . 8 . 3 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 0 . 0 . 0 . 0 . 0 . 0 . 0 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 1 data after comm

   const char* final1_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 1 . 0 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 1 . 0 . 8 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 1 . 8 . 1 . 8 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 1 . 2 . 8 . 3 . 3 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 1 . 1 . 1 . 1 . 1 . 1 . 1 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 2 data after comm

   const char* final2_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . 2 . 2 . 0 . 0 . 8 . 1 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . 2 . 2 . 2 . 8 . 2 . 8 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . 2 . 2 . 2 . 2 . 8 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . 2 . 2 . 2 . 2 . 2 . 3 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . 2 . 2 . 2 . 2 . 2 . 2 . 2 . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   // expected patch 3 data after comm

   const char* final3_txt =
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 9
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 8
      ". . . . . . . . . . . . . . . . . . . . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 7
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 6
      ". . . . . . 3 . 0 . 8 . 1 . 1 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 5
      ". . . . . . 3 . 8 . 3 . 8 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 4
      ". . . . . . 3 . 2 . 8 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 3
      ". . . . . . 3 . 2 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 2
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 1
      ". . . . . . 3 . 3 . 3 . 3 . 3 . 3 . 3 . ."
      ".   .   .   .   .   .   .   .   .   .   ." // 0
      ". . . . . . . . . . . . . . . . . . . . ."

      // 0   1   2   3   4   5   6   7   8   9
   ;

   const char* initial_txt[4] = { initial0_txt, initial1_txt,
                                  initial2_txt, initial3_txt };
   const char* final_txt[4] = { final0_txt, final1_txt,
                                final2_txt, final3_txt };

   tbox::Dimension dim(2);

   std::shared_ptr<pdat::NodeVariable<int> > var(
      new pdat::NodeVariable<int>(dim, "secondnode"));

   std::shared_ptr<pdat::SecondLayerNodeVariableFillPattern> fill_pattern(
      new pdat::SecondLayerNodeVariableFillPattern(dim));

   return SingleLevelTestCase(levelboxes_txt,
      initial_txt,
      final_txt,
      var,
      fill_pattern,
      dim);
}

int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int failures = 0;

   failures += Test_FirstLayerCellNoCornersVariableFillPattern();
   failures += Test_FirstLayerCellVariableFillPattern();
   failures += Test_SecondLayerNodeNoCornersVariableFillPattern();
   failures += Test_SecondLayerNodeVariableFillPattern();

   if (failures == 0) {
      tbox::pout << "\nPASSED:  fill_pattern" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return failures;
}
