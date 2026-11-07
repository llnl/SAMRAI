/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for testing SAMRAI data access
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <stdlib.h>

#include <concepts>

using namespace SAMRAI;

namespace {

template <class L, class R>
concept HasPlus = requires(const L& lhs, const R& rhs) {
   lhs + rhs;
};

template <class L, class R>
concept HasPlusEquals = requires(L& lhs, const R& rhs) {
   lhs += rhs;
};

template <class L, class R>
concept HasMinus = requires(const L& lhs, const R& rhs) {
   lhs - rhs;
};

template <class L, class R>
concept HasMinusEquals = requires(L& lhs, const R& rhs) {
   lhs -= rhs;
};

static_assert(HasPlus<pdat::CellIndex, hier::Index>);
static_assert(HasPlus<pdat::CellIndex, pdat::CellIndex>);
static_assert(!HasPlus<pdat::CellIndex, pdat::NodeIndex>);
static_assert(!HasPlus<pdat::CellIndex, pdat::FaceIndex>);
static_assert(!HasPlus<pdat::CellIndex, pdat::SideIndex>);
static_assert(!HasPlus<pdat::CellIndex, pdat::EdgeIndex>);
static_assert(HasPlusEquals<pdat::CellIndex, hier::Index>);
static_assert(HasPlusEquals<pdat::CellIndex, pdat::CellIndex>);
static_assert(HasMinus<pdat::CellIndex, hier::Index>);
static_assert(HasMinus<pdat::CellIndex, pdat::CellIndex>);
static_assert(!HasMinus<pdat::CellIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::CellIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::CellIndex, pdat::SideIndex>);
static_assert(!HasMinus<pdat::CellIndex, pdat::EdgeIndex>);
static_assert(HasMinusEquals<pdat::CellIndex, hier::Index>);
static_assert(HasMinusEquals<pdat::CellIndex, pdat::CellIndex>);

static_assert(!HasPlus<pdat::NodeIndex, hier::Index>);
static_assert(HasPlus<pdat::NodeIndex, pdat::NodeIndex>);
static_assert(!HasPlus<pdat::NodeIndex, pdat::CellIndex>);
static_assert(!HasPlus<pdat::NodeIndex, pdat::FaceIndex>);
static_assert(!HasPlus<pdat::NodeIndex, pdat::SideIndex>);
static_assert(!HasPlus<pdat::NodeIndex, pdat::EdgeIndex>);
static_assert(!HasPlusEquals<pdat::NodeIndex, hier::Index>);
static_assert(HasPlusEquals<pdat::NodeIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::NodeIndex, hier::Index>);
static_assert(HasMinus<pdat::NodeIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::NodeIndex, pdat::CellIndex>);
static_assert(!HasMinus<pdat::NodeIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::NodeIndex, pdat::SideIndex>);
static_assert(!HasMinus<pdat::NodeIndex, pdat::EdgeIndex>);
static_assert(!HasMinusEquals<pdat::NodeIndex, hier::Index>);
static_assert(HasMinusEquals<pdat::NodeIndex, pdat::NodeIndex>);

static_assert(!HasPlus<pdat::FaceIndex, hier::Index>);
static_assert(HasPlus<pdat::FaceIndex, pdat::FaceIndex>);
static_assert(!HasPlus<pdat::FaceIndex, pdat::CellIndex>);
static_assert(!HasPlus<pdat::FaceIndex, pdat::NodeIndex>);
static_assert(!HasPlus<pdat::FaceIndex, pdat::SideIndex>);
static_assert(!HasPlus<pdat::FaceIndex, pdat::EdgeIndex>);
static_assert(!HasPlusEquals<pdat::FaceIndex, hier::Index>);
static_assert(HasPlusEquals<pdat::FaceIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::FaceIndex, hier::Index>);
static_assert(HasMinus<pdat::FaceIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::FaceIndex, pdat::CellIndex>);
static_assert(!HasMinus<pdat::FaceIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::FaceIndex, pdat::SideIndex>);
static_assert(!HasMinus<pdat::FaceIndex, pdat::EdgeIndex>);
static_assert(!HasMinusEquals<pdat::FaceIndex, hier::Index>);
static_assert(HasMinusEquals<pdat::FaceIndex, pdat::FaceIndex>);

static_assert(!HasPlus<pdat::SideIndex, hier::Index>);
static_assert(HasPlus<pdat::SideIndex, pdat::SideIndex>);
static_assert(!HasPlus<pdat::SideIndex, pdat::CellIndex>);
static_assert(!HasPlus<pdat::SideIndex, pdat::NodeIndex>);
static_assert(!HasPlus<pdat::SideIndex, pdat::FaceIndex>);
static_assert(!HasPlus<pdat::SideIndex, pdat::EdgeIndex>);
static_assert(!HasPlusEquals<pdat::SideIndex, hier::Index>);
static_assert(HasPlusEquals<pdat::SideIndex, pdat::SideIndex>);
static_assert(!HasMinus<pdat::SideIndex, hier::Index>);
static_assert(HasMinus<pdat::SideIndex, pdat::SideIndex>);
static_assert(!HasMinus<pdat::SideIndex, pdat::CellIndex>);
static_assert(!HasMinus<pdat::SideIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::SideIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::SideIndex, pdat::EdgeIndex>);
static_assert(!HasMinusEquals<pdat::SideIndex, hier::Index>);
static_assert(HasMinusEquals<pdat::SideIndex, pdat::SideIndex>);

static_assert(!HasPlus<pdat::EdgeIndex, hier::Index>);
static_assert(HasPlus<pdat::EdgeIndex, pdat::EdgeIndex>);
static_assert(!HasPlus<pdat::EdgeIndex, pdat::CellIndex>);
static_assert(!HasPlus<pdat::EdgeIndex, pdat::NodeIndex>);
static_assert(!HasPlus<pdat::EdgeIndex, pdat::FaceIndex>);
static_assert(!HasPlus<pdat::EdgeIndex, pdat::SideIndex>);
static_assert(!HasPlusEquals<pdat::EdgeIndex, hier::Index>);
static_assert(HasPlusEquals<pdat::EdgeIndex, pdat::EdgeIndex>);
static_assert(!HasMinus<pdat::EdgeIndex, hier::Index>);
static_assert(HasMinus<pdat::EdgeIndex, pdat::EdgeIndex>);
static_assert(!HasMinus<pdat::EdgeIndex, pdat::CellIndex>);
static_assert(!HasMinus<pdat::EdgeIndex, pdat::NodeIndex>);
static_assert(!HasMinus<pdat::EdgeIndex, pdat::FaceIndex>);
static_assert(!HasMinus<pdat::EdgeIndex, pdat::SideIndex>);
static_assert(!HasMinusEquals<pdat::EdgeIndex, hier::Index>);
static_assert(HasMinusEquals<pdat::EdgeIndex, pdat::EdgeIndex>);

}  // namespace

int main(
   int argc,
   char* argv[])
{
   int error_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      /*
       * Process command line arguments.
       */
      std::string input_filename;

      if (argc != 2) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << std::endl;
         exit(-1);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */
      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /**************************************************************************
      * Read input data and setup objects.
      **************************************************************************/

      /*
       * Retreive "Main" section of input db.
       */
      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim((unsigned short)main_db->getInteger("dim"));

      const std::string log_fn =
         "iteratortest" + tbox::Utilities::intToString(dim.getValue()) + "d.log";
      tbox::PIO::logAllNodes(log_fn);

      /*
       * Regular pointer tests.
       */

      hier::Index box_lower(dim, 0);
      hier::Index box_upper(dim);

      for (int d = 0; d < dim.getValue(); ++d) {
         box_upper(d) = (d + 4) * 3;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));

      pdat::CellData<double> cell_data(box, 1, hier::IntVector(dim, 0));
      pdat::FaceData<double> face_data(box, 1, hier::IntVector(dim, 0));
      pdat::NodeData<double> node_data(box, 1, hier::IntVector(dim, 0));
      pdat::EdgeData<double> edge_data(box, 1, hier::IntVector(dim, 0));
      pdat::SideData<double> side_data(box, 1, hier::IntVector(dim, 0));

      /*
       * The tests of the iterators first fill the patch data by directly
       * accessing the pointer to the double data.  The iterators should access
       * the data in exactly the same order that it was filled.
       */

      /*
       * The tests of the index classes loop over the box with a CellIterator.
       * The CellIterator gives each cell-centered index within the box.
       * For each cell-centered index, the tests loop over every possible
       * set of parameters for the datatype-specific index class (EdgeIndex,
       * SideIndex, etc.).  Then the test checks that accessing the data using
       * the index class retrieves data from the correct point in the data array.
       */

      /*
       * Test CellIterator
       */

      double* cell_ptr = cell_data.getPointer();

      size_t cell_data_size = box.size();

      for (size_t i = 0; i < cell_data_size; ++i) {
         cell_ptr[i] = (double)i;
      }

      int j = 0;
      pdat::CellIterator ciend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator ci(pdat::CellGeometry::begin(box));
           ci != ciend; ++ci) {
         if (!tbox::MathUtilities<double>::equalEps(cell_data(*ci),
                cell_ptr[j])) {
            tbox::perr << "FAILED: - CellIterator test" << std::endl;
            ++error_count;
         }
         ++j;
      }

      /*
       * Test FaceIterator
       */

      double* face_ptr[SAMRAI::MAX_DIM_VAL];

      for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {

         face_ptr[axis] = face_data.getPointer(axis);

         hier::Box face_box = pdat::FaceGeometry::toFaceBox(box,
               axis);

         size_t face_data_size = face_box.size();

         for (size_t i = 0; i < face_data_size; ++i) {
            face_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::FaceIterator fiend(pdat::FaceGeometry::end(box, axis));
         for (pdat::FaceIterator fi(pdat::FaceGeometry::begin(box, axis));
              fi != fiend; ++fi) {
            if (!tbox::MathUtilities<double>::equalEps(face_data(*fi),
                   face_ptr[axis][j])) {
               tbox::perr << "FAILED: - FaceIterator test" << std::endl;
               ++error_count;
            }
            ++j;
         }
      }

      /*
       * Test FaceIndex
       */

      pdat::CellIterator ifcend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator ifc(pdat::CellGeometry::begin(box));
           ifc != ifcend; ++ifc) {

         for (tbox::Dimension::dir_t a = 0; a < dim.getValue(); ++a) {

            hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, a);
            hier::Index flo = face_box.lower();
            hier::Index fhi = face_box.upper();
            for (int f = 0; f < 2; ++f) {

               pdat::FaceIndex findx(*ifc, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; --i) {
                  offset = (fhi(i - 1) - flo(i - 1) + 1)
                     * (findx(i) - flo(i) + offset);
               }
               offset += findx(0) - flo(0);

               if (!tbox::MathUtilities<double>::equalEps(face_data(findx),
                      face_ptr[a][offset])) {
                  tbox::perr << "FAILED: - FaceIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      /*
       * Test NodeIterator
       */

      double* node_ptr = node_data.getPointer();

      hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);

      hier::Index nlo = node_box.lower();
      hier::Index nhi = node_box.upper();

      size_t node_data_size = node_box.size();

      for (size_t i = 0; i < node_data_size; ++i) {
         node_ptr[i] = (double)i;
      }

      j = 0;
      pdat::NodeIterator niend(pdat::NodeGeometry::end(box));
      for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(box));
           ni != niend; ++ni) {
         if (!tbox::MathUtilities<double>::equalEps(node_data(*ni),
                node_ptr[j])) {
            tbox::perr << "FAILED: - NodeIterator test" << std::endl;
            ++error_count;
         }
         ++j;
      }

      /*
       * Test NodeIndex
       */

      pdat::CellIterator incend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator inc(pdat::CellGeometry::begin(box));
           inc != incend; ++inc) {

         hier::IntVector corner(dim, 0);

         bool all_corners_complete = false;

         while (!all_corners_complete) {

            pdat::NodeIndex nindx(*inc, corner);

            int offset = 0;
            for (int i = dim.getValue() - 1; i > 0; --i) {
               offset = (nhi(i - 1) - nlo(i - 1) + 1)
                  * (nindx(i) - nlo(i) + offset);
            }
            offset += nindx(0) - nlo(0);

            if (!tbox::MathUtilities<double>::equalEps(node_data(nindx),
                   node_ptr[offset])) {
               tbox::perr << "FAILED: - NodeIndex test" << std::endl;
               ++error_count;
            }

            int u;
            for (u = 0; u < dim.getValue(); ++u) {
               if (corner(u) == 1) {
                  corner(u) = 0;
               } else {
                  ++corner(u);
                  break;
               }
            }
            if (u == dim.getValue()) {
               all_corners_complete = true;
            }
         }
      }

      /*
       * Test EdgeIterator
       */

      double* edge_ptr[SAMRAI::MAX_DIM_VAL];

      for (int axis = 0; axis < dim.getValue(); ++axis) {

         edge_ptr[axis] = edge_data.getPointer(axis);

         hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box,
               axis);
         size_t edge_data_size = edge_box.size();

         for (size_t i = 0; i < edge_data_size; ++i) {
            edge_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::EdgeIterator eiend(pdat::EdgeGeometry::end(box, axis));
         for (pdat::EdgeIterator ei(pdat::EdgeGeometry::begin(box, axis));
              ei != eiend; ++ei) {
            if (!tbox::MathUtilities<double>::equalEps(edge_data(*ei),
                   edge_ptr[axis][j])) {
               tbox::perr << "FAILED: - EdgeIterator test" << std::endl;
               ++error_count;
            }
            ++j;
         }
      }

      /*
       * Test EdgeIndex
       */

      pdat::CellIterator iecend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator iec(pdat::CellGeometry::begin(box));
           iec != iecend; ++iec) {
         for (int a = 0; a < dim.getValue(); ++a) {

            hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box,
                  a);
            hier::Index elo = edge_box.lower();
            hier::Index ehi = edge_box.upper();

            for (int f = 0; f < (1 << (dim.getValue() - 1)); ++f) {
               pdat::EdgeIndex eindx(*iec, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; --i) {
                  offset = (ehi(i - 1) - elo(i - 1) + 1)
                     * (eindx(i) - elo(i) + offset);
               }
               offset += eindx(0) - elo(0);

               if (!tbox::MathUtilities<double>::equalEps(edge_data(eindx),
                      edge_ptr[a][offset])) {
                  tbox::perr << "FAILED: - EdgeIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      /*
       * Test SideIterator
       */

      double* side_ptr[SAMRAI::MAX_DIM_VAL];

      for (tbox::Dimension::dir_t axis = 0; axis < dim.getValue(); ++axis) {

         side_ptr[axis] = side_data.getPointer(axis);

         hier::Box side_box = pdat::SideGeometry::toSideBox(box, axis);
         size_t side_data_size = side_box.size();

         for (size_t i = 0; i < side_data_size; ++i) {
            side_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::SideIterator siend(pdat::SideGeometry::end(box, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(box, axis));
              si != siend; ++si) {
            if (!tbox::MathUtilities<double>::equalEps(side_data(*si),
                   side_ptr[axis][j])) {
               tbox::perr << "FAILED: - SideIterator test" << std::endl;
               ++error_count;
            }
            ++j;
         }
      }

      /*
       * Test SideIndex
       */

      pdat::CellIterator iscend(pdat::CellGeometry::end(box));
      for (pdat::CellIterator isc(pdat::CellGeometry::begin(box));
           isc != iscend; ++isc) {
         for (tbox::Dimension::dir_t a = 0; a < dim.getValue(); ++a) {

            hier::Box side_box = pdat::SideGeometry::toSideBox(box, a);
            hier::Index slo = side_box.lower();
            hier::Index shi = side_box.upper();

            for (int f = 0; f < 2; ++f) {
               pdat::SideIndex sindx(*isc, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; --i) {
                  offset = (shi(i - 1) - slo(i - 1) + 1)
                     * (sindx(i) - slo(i) + offset);
               }
               offset += sindx(0) - slo(0);

               if (!tbox::MathUtilities<double>::equalEps(side_data(sindx),
                      side_ptr[a][offset])) {
                  tbox::perr << "FAILED: - SideIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      if (error_count == 0) {
         tbox::pout << "\nPASSED:  dataaccess" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return error_count;
}
