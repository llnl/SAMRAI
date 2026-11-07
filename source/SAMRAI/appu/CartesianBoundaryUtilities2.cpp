/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating 2D Cartesian boundary data
 *
 ************************************************************************/

#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"

#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in
 * boundary condition implementation.
 *
 *************************************************************************
 */

extern "C" {

// Disable Intel warnings about extern declarations in C file
#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

void SAMRAI_F77_FUNC(stufcartbdryloc2d, STUFCARTBDRYLOC2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&);

void SAMRAI_F77_FUNC(stufcartbdrycond2d, STUFCARTBDRYCOND2D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&);

void SAMRAI_F77_FUNC(getcartedgebdry2d, GETCARTEDGEBDRY2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const int&,
   const int&,
   const double *,
   double *,
   const int&);

void SAMRAI_F77_FUNC(getcartnodebdry2d, GETCARTNODEBDRY2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   const int&,
   const int&,
   const double *,
   double *,
   const int&);

}

namespace SAMRAI {
namespace appu {

bool CartesianBoundaryUtilities2::s_fortran_constants_stuffed = false;

/*
 * This function reads 2D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control
 * is passed to the BoundaryStrategy to read the boundary state data
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic
 * domain information is used to determine which boundary edges or
 * node entries are not required from input.  Error checking requires
 * that node boundary conditions are consistent with those
 * specified for the edges.
 *
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN data
 *    bdry_db .......... input database containing all boundary data
 *    edge_conds ....... array into which integer boundary conditions
 *                       for edges are read
 *    node_conds ....... array into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry2::getPeriodicShift())
 */

void
CartesianBoundaryUtilities2::getFromInput(
   BoundaryUtilityStrategy* bdry_strategy,
   const std::shared_ptr<tbox::Database>& input_db,
   std::vector<int>& edge_conds,
   std::vector<int>& node_conds,
   const hier::IntVector& periodic)

{
   TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
   TBOX_ASSERT(bdry_strategy != 0);
   TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
   TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);

   if (!input_db) {
      TBOX_ERROR(": CartesianBoundaryUtility2::getFromInput()\n"
         << "no input database supplied" << std::endl);
   }

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   read2dBdryEdges(bdry_strategy,
      input_db,
      edge_conds,
      periodic);

   read2dBdryNodes(input_db,
      edge_conds,
      node_conds,
      periodic);

}

/*
 * Function to fill edge boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_edge_conds ...... array of boundary conditions for patch edges
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void
CartesianBoundaryUtilities2::fillEdgeBoundaryData(
   const std::string& varname,
   const std::shared_ptr<pdat::CellData<double> >& vardata,
   const hier::Patch& patch,
   const hier::IntVector& ghost_fill_width,
   const std::vector<int>& bdry_edge_conds,
   const std::vector<double>& bdry_edge_values)
{
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(vardata);
   TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_2D_EDGES);
   TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) ==
      NUM_2D_EDGES * (vardata->getDepth()));

   TBOX_DIM_ASSERT(ghost_fill_width.getDim() == tbox::Dimension(2));
   TBOX_ASSERT_OBJDIM_EQUALITY3(*vardata, patch, ghost_fill_width);

   NULL_USE(varname);

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const double* dx = pgeom->getDx();

   const hier::Box& interior(patch.getBox());
   const hier::Index& ifirst(interior.lower());
   const hier::Index& ilast(interior.upper());

   const hier::IntVector& ghost_cells(vardata->getGhostCellWidth());

   hier::IntVector gcw_to_fill(hier::IntVector::min(ghost_cells,
                                  ghost_fill_width));

   const std::vector<hier::BoundaryBox>& edge_bdry =
      pgeom->getCodimensionBoundaries(Bdry::EDGE2D);
   for (int i = 0; i < static_cast<int>(edge_bdry.size()); ++i) {

      TBOX_ASSERT(edge_bdry[i].getBoundaryType() == Bdry::EDGE2D);

      int bedge_loc = edge_bdry[i].getLocationIndex();

      hier::Box fill_box(pgeom->getBoundaryFillBox(edge_bdry[i],
                            interior,
                            gcw_to_fill));
      const hier::Index& ibeg(fill_box.lower());
      const hier::Index& iend(fill_box.upper());

      SAMRAI_F77_FUNC(getcartedgebdry2d, GETCARTEDGEBDRY2D) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ibeg(0), iend(0),
         ibeg(1), iend(1),
         ghost_cells(0), ghost_cells(1),
         dx,
         bedge_loc,
         bdry_edge_conds[bedge_loc],
         &bdry_edge_values[0],
         vardata->getPointer(),
         vardata->getDepth());

   }

}

/*
 * Function to fill node boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 *    bdry_edge_values ..... array of boundary values for edges
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void
CartesianBoundaryUtilities2::fillNodeBoundaryData(
   const std::string& varname,
   const std::shared_ptr<pdat::CellData<double> >& vardata,
   const hier::Patch& patch,
   const hier::IntVector& ghost_fill_width,
   const std::vector<int>& bdry_node_conds,
   const std::vector<double>& bdry_edge_values)
{
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(vardata);
   TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_2D_NODES);
   TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) ==
      NUM_2D_EDGES * (vardata->getDepth()));

   TBOX_DIM_ASSERT(ghost_fill_width.getDim() == tbox::Dimension(2));
   TBOX_ASSERT_OBJDIM_EQUALITY3(*vardata, patch, ghost_fill_width);

   NULL_USE(varname);

   if (!s_fortran_constants_stuffed) {
      stuff2dBdryFortConst();
   }

   const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const double* dx = pgeom->getDx();

   const hier::Box& interior(patch.getBox());
   const hier::Index& ifirst(interior.lower());
   const hier::Index& ilast(interior.upper());

   const hier::IntVector& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_fill_width);

   const std::vector<hier::BoundaryBox>& node_bdry =
      pgeom->getCodimensionBoundaries(Bdry::NODE2D);
   for (int i = 0; i < static_cast<int>(node_bdry.size()); ++i) {

      TBOX_ASSERT(node_bdry[i].getBoundaryType() == Bdry::NODE2D);

      int bnode_loc = node_bdry[i].getLocationIndex();

      hier::Box fill_box(pgeom->getBoundaryFillBox(node_bdry[i],
                            interior,
                            gcw_to_fill));
      const hier::Index& ibeg(fill_box.lower());
      const hier::Index& iend(fill_box.upper());

      SAMRAI_F77_FUNC(getcartnodebdry2d, GETCARTNODEBDRY2D) (ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ibeg(0), iend(0),
         ibeg(1), iend(1),
         ghost_cells(0), ghost_cells(1),
         dx,
         bnode_loc,
         bdry_node_conds[bnode_loc],
         &bdry_edge_values[0],
         vardata->getPointer(),
         vardata->getDepth());

   }

}

/*
 * Function that returns the integer edge boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown,
 * or the boundary condition type is inconsistant with the node location
 * an error results.
 */

int
CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
   int node_loc,
   int node_btype)
{

   int ret_edge = -1;

   switch (node_btype) {
      case BdryCond::XFLOW:
      case BdryCond::XREFLECT:
      case BdryCond::XDIRICHLET:
      case BdryCond::XNEUMANN:
      {
         if (node_loc == NodeBdyLoc2D::XLO_YLO ||
             node_loc == NodeBdyLoc2D::XLO_YHI) {
            ret_edge = BdryLoc::XLO;
         } else {
            ret_edge = BdryLoc::XHI;
         }
         break;
      }
      case BdryCond::YFLOW:
      case BdryCond::YREFLECT:
      case BdryCond::YDIRICHLET:
      case BdryCond::YNEUMANN:
      {
         if (node_loc == NodeBdyLoc2D::XLO_YLO ||
             node_loc == NodeBdyLoc2D::XHI_YLO) {
            ret_edge = BdryLoc::YLO;
         } else {
            ret_edge = BdryLoc::YHI;
         }
         break;
      }
      default: {
         TBOX_ERROR("Unknown node boundary condition type = "
            << node_btype << " passed to \n"
            << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
            << std::endl);
      }
   }

   if (ret_edge == -1) {
      TBOX_ERROR("Node boundary condition type = "
         << node_btype << " and node location = " << node_loc
         << "\n passed to "
         << "CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry"
         << " are inconsistant." << std::endl);
   }

   return ret_edge;

}

/*
 * Function to check 2D boundary data filling.  Arguments are:
 *
 *    varname ..... name of variable (for error reporting)
 *    patch ....... patch on which boundary data to check lives
 *    data_id ..... patch data index on patch
 *    depth ....... depth index of data to check
 *    gcw_to_check. boundary ghost width to fill
 *    bbox ........ boundary box to check
 *    bcase ....... boundary condition case for edge or a node to check
 *    bstate ...... boundary state that applies when such a value is
 *                  required, such as when using Dirichlet conditions
 */

int
CartesianBoundaryUtilities2::checkBdryData(
   const std::string& varname,
   const hier::Patch& patch,
   int data_id,
   int depth,
   const hier::IntVector& gcw_to_check,
   const hier::BoundaryBox& bbox,
   int bcase,
   const double& bstate)
{
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(data_id >= 0);
   TBOX_ASSERT(depth >= 0);

   TBOX_DIM_ASSERT(gcw_to_check.getDim() == tbox::Dimension(2));
   TBOX_ASSERT_OBJDIM_EQUALITY3(patch, gcw_to_check, bbox);

   int num_bad_values = 0;

   int btype = bbox.getBoundaryType();
   int bloc = bbox.getLocationIndex();

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);

   std::shared_ptr<pdat::CellData<double> > vardata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(data_id)));
   TBOX_ASSERT(vardata);

   std::string bdry_type_str;
   if (btype == Bdry::EDGE2D) {
      bdry_type_str = "EDGE";
   } else if (btype == Bdry::NODE2D) {
      bdry_type_str = "NODE";
   } else {
      TBOX_ERROR(
         "Unknown btype " << btype
                          << " passed to CartesianBoundaryUtilities2::checkBdryData()! "
                          << std::endl);
   }

   tbox::plog << "\n\nCHECKING 2D " << bdry_type_str << " BDRY DATA..."
              << std::endl;
   tbox::plog << "varname = " << varname << " : depth = " << depth << std::endl;
   tbox::plog << "bbox = " << bbox.getBox() << std::endl;
   tbox::plog << "btype, bloc, bcase = "
              << btype << ", = " << bloc << ", = " << bcase << std::endl;

   tbox::Dimension::dir_t idir;
   double valfact = 0.0, constval = 0.0, dxfact = 0.0;
   int offsign;

   get2dBdryDirectionCheckValues(idir, offsign,
      btype, bloc, bcase);

   if (btype == Bdry::EDGE2D) {

      if (bcase == BdryCond::FLOW) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::REFLECT) {
         valfact = -1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::DIRICHLET) {
         valfact = 0.0;
         constval = bstate;
         dxfact = 0.0;
      } else if (bcase == BdryCond::NEUMANN) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR(
            "Unknown bcase " << bcase
                             << " passed to CartesianBoundaryUtilities2::checkBdryData()"
                             << "\n for " << bdry_type_str
                             << " at location " << bloc << std::endl);
      }

   } else if (btype == Bdry::NODE2D) {

      if (bcase == BdryCond::XFLOW || bcase == BdryCond::YFLOW) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XREFLECT || bcase == BdryCond::YREFLECT) {
         valfact = -1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XDIRICHLET ||
                 bcase == BdryCond::YDIRICHLET) {
         valfact = 0.0;
         constval = bstate;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XNEUMANN || bcase == BdryCond::YNEUMANN) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = bstate * pgeom->getDx()[idir];
      } else {
         TBOX_ERROR(
            "Unknown bcase " << bcase
                             << " passed to CartesianBoundaryUtilities2::checkBdryData()"
                             << "\n for " << bdry_type_str
                             << " at location " << bloc << std::endl);
      }

   }

   hier::Box gbox_to_check =
      vardata->getGhostBox() * pgeom->getBoundaryFillBox(bbox,
         patch.getBox(),
         gcw_to_check);

   hier::Box cbox(gbox_to_check);
   hier::Box dbox(gbox_to_check);
   hier::Index ifirst(vardata->getBox().lower());
   hier::Index ilast(vardata->getBox().upper());

   if (offsign == -1) {
      cbox.setLower(idir, ifirst(idir) - 1);
      cbox.setUpper(idir, ifirst(idir) - 1);
      dbox.setLower(idir, ifirst(idir));
      dbox.setUpper(idir, ifirst(idir));
   } else {
      cbox.setLower(idir, ilast(idir) + 1);
      cbox.setUpper(idir, ilast(idir) + 1);
      dbox.setLower(idir, ilast(idir));
      dbox.setUpper(idir, ilast(idir));
   }

   pdat::CellIterator id(pdat::CellGeometry::begin(dbox));
   pdat::CellIterator icend(pdat::CellGeometry::end(cbox));
   for (pdat::CellIterator ic(pdat::CellGeometry::begin(cbox));
        ic != icend; ++ic) {
      double checkval = valfact * (*vardata)(*id, depth) + constval;
      pdat::CellIndex check = *ic;
      for (int p = 0; p < gbox_to_check.numberCells(idir); ++p) {
         double offcheckval = checkval + dxfact * (p + 1);

#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
         if ((*vardata)(check, depth) != offcheckval) {
            ++num_bad_values;
            TBOX_WARNING("Bad " << bdry_type_str
                                << " boundary value for " << varname
                                << " found in cell " << check
                                << "\n   found = " << (*vardata)(check, depth)
                                << " : correct = " << offcheckval << std::endl);
         }
         check(idir) += offsign;
      }
      ++id;
   }

   return num_bad_values;

}

/*
 * Private function to read 2D edge boundary data from input database.
 */

void
CartesianBoundaryUtilities2::read2dBdryEdges(
   BoundaryUtilityStrategy* bdry_strategy,
   const std::shared_ptr<tbox::Database>& input_db,
   std::vector<int>& edge_conds,
   const hier::IntVector& periodic)
{
   TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));

   TBOX_ASSERT(bdry_strategy != 0);
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);

   int num_per_dirs = 0;
   for (int id = 0; id < 2; ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < 2) { // face boundary input required

      for (int s = 0; s < NUM_2D_EDGES; ++s) {

         std::string bdry_loc_str;
         switch (s) {
            case BdryLoc::XLO: { bdry_loc_str = "boundary_edge_xlo";
                                 break;
            }
            case BdryLoc::XHI: { bdry_loc_str = "boundary_edge_xhi";
                                 break;
            }
            case BdryLoc::YLO: { bdry_loc_str = "boundary_edge_ylo";
                                 break;
            }
            case BdryLoc::YHI: { bdry_loc_str = "boundary_edge_yhi";
                                 break;
            }
            default: NULL_STATEMENT;
         }

         bool need_data_read = true;
         if (num_per_dirs > 0) {
            if (periodic(0) && (s == BdryLoc::XLO || s == BdryLoc::XHI)) {
               need_data_read = false;
            } else if (periodic(1) && (s == BdryLoc::YLO || s == BdryLoc::YHI)) {
               need_data_read = false;
            }
         }

         if (need_data_read) {
            std::shared_ptr<tbox::Database> bdry_loc_db(
               input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
               bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "FLOW") {
               edge_conds[s] = BdryCond::FLOW;
            } else if (bdry_cond_str == "REFLECT") {
               edge_conds[s] = BdryCond::REFLECT;
            } else if (bdry_cond_str == "DIRICHLET") {
               edge_conds[s] = BdryCond::DIRICHLET;
               bdry_strategy->readDirichletBoundaryDataEntry(bdry_loc_db,
                  bdry_loc_str,
                  s);
            } else if (bdry_cond_str == "NEUMANN") {
               edge_conds[s] = BdryCond::NEUMANN;
               bdry_strategy->readNeumannBoundaryDataEntry(bdry_loc_db,
                  bdry_loc_str,
                  s);
            } else {
               TBOX_ERROR("Unknown edge boundary string = "
                  << bdry_cond_str << " found in input." << std::endl);
            }
         } // if (need_data_read)

      } // for (int s = 0 ...

   } // if (num_per_dirs < 2)

}

/*
 * Private function to read 2D node boundary data from input database.
 */

void
CartesianBoundaryUtilities2::read2dBdryNodes(
   const std::shared_ptr<tbox::Database>& input_db,
   const std::vector<int>& edge_conds,
   std::vector<int>& node_conds,
   const hier::IntVector& periodic)
{
   TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));

   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
   TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);

   int num_per_dirs = 0;
   for (int id = 0; id < 2; ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < 1) { // node boundary data required

      for (int s = 0; s < NUM_2D_NODES; ++s) {

         std::string bdry_loc_str;
         switch (s) {
            case NodeBdyLoc2D::XLO_YLO: {
               bdry_loc_str = "boundary_node_xlo_ylo";
               break;
            }
            case NodeBdyLoc2D::XHI_YLO: {
               bdry_loc_str = "boundary_node_xhi_ylo";
               break;
            }
            case NodeBdyLoc2D::XLO_YHI: {
               bdry_loc_str = "boundary_node_xlo_yhi";
               break;
            }
            case NodeBdyLoc2D::XHI_YHI: {
               bdry_loc_str = "boundary_node_xhi_yhi";
               break;
            }
            default: NULL_STATEMENT;
         }

         std::shared_ptr<tbox::Database> bdry_loc_db(
            input_db->getDatabase(bdry_loc_str));
         std::string bdry_cond_str =
            bdry_loc_db->getString("boundary_condition");
         if (bdry_cond_str == "XFLOW") {
            node_conds[s] = BdryCond::XFLOW;
         } else if (bdry_cond_str == "YFLOW") {
            node_conds[s] = BdryCond::YFLOW;
         } else if (bdry_cond_str == "XREFLECT") {
            node_conds[s] = BdryCond::XREFLECT;
         } else if (bdry_cond_str == "YREFLECT") {
            node_conds[s] = BdryCond::YREFLECT;
         } else if (bdry_cond_str == "XDIRICHLET") {
            node_conds[s] = BdryCond::XDIRICHLET;
         } else if (bdry_cond_str == "YDIRICHLET") {
            node_conds[s] = BdryCond::YDIRICHLET;
         } else if (bdry_cond_str == "XNEUMANN") {
            node_conds[s] = BdryCond::XNEUMANN;
         } else if (bdry_cond_str == "YNEUMANN") {
            node_conds[s] = BdryCond::YNEUMANN;
         } else {
            TBOX_ERROR("Unknown node boundary string = "
               << bdry_cond_str << " found in input." << std::endl);
         }

         std::string proper_edge;
         std::string proper_edge_data;
         bool no_edge_data_found = false;
         if (bdry_cond_str == "XFLOW" ||
             bdry_cond_str == "XDIRICHLET" ||
             bdry_cond_str == "XNEUMANN" ||
             bdry_cond_str == "XREFLECT") {
            if (s == NodeBdyLoc2D::XLO_YLO ||
                s == NodeBdyLoc2D::XLO_YHI) {
               proper_edge = "XLO";
               if (bdry_cond_str == "XFLOW" &&
                   edge_conds[BdryLoc::XLO] != BdryCond::FLOW) {
                  no_edge_data_found = true;
                  proper_edge_data = "FLOW";
               }
               if (bdry_cond_str == "XDIRICHLET" &&
                   edge_conds[BdryLoc::XLO] != BdryCond::DIRICHLET) {
                  no_edge_data_found = true;
                  proper_edge_data = "DIRICHLET";
               }
               if (bdry_cond_str == "XNEUMANN" &&
                   edge_conds[BdryLoc::XLO] != BdryCond::NEUMANN) {
                  no_edge_data_found = true;
                  proper_edge_data = "NEUMANN";
               }
               if (bdry_cond_str == "XREFLECT" &&
                   edge_conds[BdryLoc::XLO] != BdryCond::REFLECT) {
                  no_edge_data_found = true;
                  proper_edge_data = "REFLECT";
               }
            } else {
               proper_edge = "XHI";
               if (bdry_cond_str == "XFLOW" &&
                   edge_conds[BdryLoc::XHI] != BdryCond::FLOW) {
                  no_edge_data_found = true;
                  proper_edge_data = "FLOW";
               }
               if (bdry_cond_str == "XDIRICHLET" &&
                   edge_conds[BdryLoc::XHI] != BdryCond::DIRICHLET) {
                  no_edge_data_found = true;
                  proper_edge_data = "DIRICHLET";
               }
               if (bdry_cond_str == "XNEUMANN" &&
                   edge_conds[BdryLoc::XHI] != BdryCond::NEUMANN) {
                  no_edge_data_found = true;
                  proper_edge_data = "NEUMANN";
               }
               if (bdry_cond_str == "XREFLECT" &&
                   edge_conds[BdryLoc::XHI] != BdryCond::REFLECT) {
                  no_edge_data_found = true;
                  proper_edge_data = "REFLECT";
               }
            }
         } else if (bdry_cond_str == "YFLOW" ||
                    bdry_cond_str == "YDIRICHLET" ||
                    bdry_cond_str == "YNEUMANN" ||
                    bdry_cond_str == "YREFLECT") {
            if (s == NodeBdyLoc2D::XLO_YLO ||
                s == NodeBdyLoc2D::XHI_YLO) {
               proper_edge = "YLO";
               if (bdry_cond_str == "YFLOW" &&
                   edge_conds[BdryLoc::YLO] != BdryCond::FLOW) {
                  no_edge_data_found = true;
                  proper_edge_data = "FLOW";
               }
               if (bdry_cond_str == "YDIRICHLET" &&
                   edge_conds[BdryLoc::YLO] != BdryCond::DIRICHLET) {
                  no_edge_data_found = true;
                  proper_edge_data = "DIRICHLET";
               }
               if (bdry_cond_str == "YNEUMANN" &&
                   edge_conds[BdryLoc::YLO] != BdryCond::NEUMANN) {
                  no_edge_data_found = true;
                  proper_edge_data = "NEUMANN";
               }
               if (bdry_cond_str == "YREFLECT" &&
                   edge_conds[BdryLoc::YLO] != BdryCond::REFLECT) {
                  no_edge_data_found = true;
                  proper_edge_data = "REFLECT";
               }
            } else {
               proper_edge = "YHI";
               if (bdry_cond_str == "YFLOW" &&
                   edge_conds[BdryLoc::YHI] != BdryCond::FLOW) {
                  no_edge_data_found = true;
                  proper_edge_data = "FLOW";
               }
               if (bdry_cond_str == "YDIRICHLET" &&
                   edge_conds[BdryLoc::YHI] != BdryCond::DIRICHLET) {
                  no_edge_data_found = true;
                  proper_edge_data = "DIRICHLET";
               }
               if (bdry_cond_str == "YNEUMANN" &&
                   edge_conds[BdryLoc::YHI] != BdryCond::NEUMANN) {
                  no_edge_data_found = true;
                  proper_edge_data = "NEUMANN";
               }
               if (bdry_cond_str == "YREFLECT" &&
                   edge_conds[BdryLoc::YHI] != BdryCond::REFLECT) {
                  no_edge_data_found = true;
                  proper_edge_data = "REFLECT";
               }
            }
         }
         if (no_edge_data_found) {
            TBOX_ERROR(
               "Bdry condition " << bdry_cond_str
                                 << " found for "
                                 << bdry_loc_str
                                 << "\n but no "
                                 << proper_edge_data
                                 << " data found for edge "
                                 << proper_edge << std::endl);
         }

      } // for (int s = 0 ...

   } // if (num_per_dirs < 1)

}

/*
 * Private function to get boundary orientation information for
 * 2D boundary condition checking.  Called from checkBdryData().
 */

void
CartesianBoundaryUtilities2::get2dBdryDirectionCheckValues(
   tbox::Dimension::dir_t& idir,
   int& offsign,
   int btype,
   int bloc,
   int bcase)
{

   std::string bdry_type_str;

   if (btype == Bdry::EDGE2D) {

      bdry_type_str = "NODE";

      if (bloc == BdryLoc::XLO || bloc == BdryLoc::XHI) {
         idir = 0;
         if (bloc == BdryLoc::XLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bloc == BdryLoc::YLO || bloc == BdryLoc::YHI) {
         idir = 1;
         if (bloc == BdryLoc::YLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else {
         TBOX_ERROR(
            "Unknown boundary location " << bloc
                                         <<
            " passed to CartesianBoundaryUtilities2::checkBdryData()"
                                         << "\n for "
                                         << bdry_type_str << " boundary " << std::endl);
      }

   } else if (btype == Bdry::NODE2D) {

      bdry_type_str = "NODE";

      if (bcase == BdryCond::XFLOW || bcase == BdryCond::XREFLECT ||
          bcase == BdryCond::XDIRICHLET || bcase == BdryCond::XNEUMANN) {
         idir = 0;
         if (bloc == NodeBdyLoc2D::XLO_YLO || bloc == NodeBdyLoc2D::XLO_YHI) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == BdryCond::YFLOW || bcase == BdryCond::YREFLECT ||
                 bcase == BdryCond::YDIRICHLET ||
                 bcase == BdryCond::YNEUMANN) {
         idir = 1;
         if (bloc == NodeBdyLoc2D::XLO_YLO || bloc == NodeBdyLoc2D::XHI_YLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      }

   } else {
      TBOX_ERROR(
         "Unknown boundary type " << btype
                                  << " passed to CartesianBoundaryUtilities2::checkBdryData()"
                                  << "\n for " << bdry_type_str
                                  << " at location " << bloc
                                  << std::endl);
   }

}

/*
 * Private function to stuff 2D boundary contants into Fortran common blocks
 */

void
CartesianBoundaryUtilities2::stuff2dBdryFortConst()
{
   SAMRAI_F77_FUNC(stufcartbdryloc2d, STUFCARTBDRYLOC2D) (BdryLoc::XLO, BdryLoc::XHI,
      BdryLoc::YLO, BdryLoc::YHI, NodeBdyLoc2D::XLO_YLO, NodeBdyLoc2D::XHI_YLO,
      NodeBdyLoc2D::XLO_YHI, NodeBdyLoc2D::XHI_YHI);
   SAMRAI_F77_FUNC(stufcartbdrycond2d, STUFCARTBDRYCOND2D) (BdryCond::FLOW,
      BdryCond::XFLOW, BdryCond::YFLOW,
      BdryCond::REFLECT,
      BdryCond::XREFLECT, BdryCond::YREFLECT,
      BdryCond::DIRICHLET,
      BdryCond::XDIRICHLET, BdryCond::YDIRICHLET,
      BdryCond::NEUMANN,
      BdryCond::XNEUMANN, BdryCond::YNEUMANN);
   s_fortran_constants_stuffed = true;
}

}
}
