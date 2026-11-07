/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating 3D Skeleton boundary data
 *
 ************************************************************************/

#include "SkeletonBoundaryUtilities3.h"

#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in
 * boundary condition implementation.
 *
 *************************************************************************
 */

extern "C" {

void SAMRAI_F77_FUNC(stufcartbdryloc3d, STUFCARTBDRYLOC3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&);

void SAMRAI_F77_FUNC(stufcartbdrycond3d, STUFCARTBDRYCOND3D) (
   const int&,
   const int&, const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&,
   const int&, const int&, const int&);

void SAMRAI_F77_FUNC(getcartfacebdry3d, GETCARTFACEBDRY3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const int&,
   const int&,
   const double *,
   double *,
   const int&);

void SAMRAI_F77_FUNC(getcartedgebdry3d, GETCARTEDGEBDRY3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const int&,
   const int&,
   const double *,
   double *,
   const int&);

void SAMRAI_F77_FUNC(getcartnodebdry3d, GETCARTNODEBDRY3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const int&,
   const int&,
   const double *,
   double *,
   const int&);

}
using namespace SAMRAI;
using namespace appu;

bool SkeletonBoundaryUtilities3::s_fortran_constants_stuffed = false;

/*
 * This function reads 3D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control
 * is passed to the BoundaryStrategy to read the boundary state data
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic
 * domain information is used to determine which boundary face, edge, or
 * node entries are not required from input.  Error checking requires
 * that node and edge boundary conditions are consistent with those
 * specified for the faces.
 *
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN conditions
 *    input_db ......... input database containing all boundary data
 *    face_conds ....... vector into which integer boundary conditions
 *                       for faces are read
 *    edge_conds ....... vector into which integer boundary conditions
 *                       for edges are read
 *    node_conds ....... vector into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry3::getPeriodicShift())
 */

void SkeletonBoundaryUtilities3::getFromInput(
   BoundaryUtilityStrategy* bdry_strategy,
   const std::shared_ptr<tbox::Database>& input_db,
   std::vector<int>& face_conds,
   std::vector<int>& edge_conds,
   std::vector<int>& node_conds,
   const hier::IntVector& periodic)
{
   TBOX_ASSERT(bdry_strategy != 0);
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
   TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);
   TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   read3dBdryFaces(bdry_strategy,
      input_db,
      face_conds,
      periodic);

   read3dBdryEdges(input_db,
      face_conds,
      edge_conds,
      periodic);

   read3dBdryNodes(input_db,
      face_conds,
      node_conds,
      periodic);

}

/*
 * Function to fill face boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_face_conds ...... vector of boundary conditions for patch faces
 *    bdry_face_values ..... vector of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void SkeletonBoundaryUtilities3::fillFaceBoundaryData(
   const std::string& varname,
   std::shared_ptr<pdat::CellData<double> >& vardata,
   const hier::Patch& patch,
   const hier::IntVector& ghost_fill_width,
   const std::vector<int>& bdry_face_conds,
   const std::vector<double>& bdry_face_values)
{
   NULL_USE(varname);
   TBOX_ASSERT(vardata);
   TBOX_ASSERT(static_cast<int>(bdry_face_conds.size()) == NUM_3D_FACES);
//   TBOX_ASSERT(static_cast<int>(static_cast<int>(bdry_face_values.size())) ==
//               NUM_3D_FACES*(vardata->getDepth()));

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());
   //const double* dx = pgeom->getDx();
   const double dx[3] = { 0., 0., 0. };

   const hier::Box& interior(patch.getBox());
   const hier::Index& ifirst(interior.lower());
   const hier::Index& ilast(interior.upper());

   const hier::IntVector& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_fill_width);
   const std::vector<hier::BoundaryBox>& face_bdry =
      pgeom->getCodimensionBoundaries(Bdry::FACE3D);
   for (int i = 0; i < static_cast<int>(face_bdry.size()); ++i) {
      TBOX_ASSERT(face_bdry[i].getBoundaryType() == Bdry::FACE3D);

      int bface_loc = face_bdry[i].getLocationIndex();

      hier::Box fill_box(pgeom->getBoundaryFillBox(face_bdry[i],
                            interior,
                            gcw_to_fill));
      const hier::Index& ibeg(fill_box.lower());
      const hier::Index& iend(fill_box.upper());

      SAMRAI_F77_FUNC(getcartfacebdry3d, GETCARTFACEBDRY3D) (
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         ibeg(0), iend(0),
         ibeg(1), iend(1),
         ibeg(2), iend(2),
         ghost_cells(0), ghost_cells(1), ghost_cells(2),
         dx,
         bface_loc,
         bdry_face_conds[bface_loc],
         &bdry_face_values[0],
         vardata->getPointer(),
         vardata->getDepth());

   }

}

/*
 * Function to fill edge boundary values.
 *
 * Arguments are:
 *    varname .............. name of variable (for error reporting)
 *    vardata .............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    ghost_width_to_fill .. width of ghost region to fill
 *    bdry_edge_conds ...... vector of boundary conditions for patch edges
 *    bdry_face_values ..... vector of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void SkeletonBoundaryUtilities3::fillEdgeBoundaryData(
   const std::string& varname,
   std::shared_ptr<pdat::CellData<double> >& vardata,
   const hier::Patch& patch,
   const hier::IntVector& ghost_fill_width,
   const std::vector<int>& bdry_edge_conds,
   const std::vector<double>& bdry_face_values)
{
   NULL_USE(varname);

   TBOX_ASSERT(vardata);
   TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_3D_EDGES);
   TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) ==
      NUM_3D_FACES * (vardata->getDepth()));

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());
   //const double* dx = pgeom->getDx();
   const double dx[3] = { 0., 0., 0. };

   const hier::Box& interior(patch.getBox());
   const hier::Index& ifirst(interior.lower());
   const hier::Index& ilast(interior.upper());

   const hier::IntVector& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_fill_width);

   const std::vector<hier::BoundaryBox>& edge_bdry =
      pgeom->getCodimensionBoundaries(Bdry::EDGE3D);
   for (int i = 0; i < static_cast<int>(edge_bdry.size()); ++i) {
      TBOX_ASSERT(edge_bdry[i].getBoundaryType() == Bdry::EDGE3D);

      int bedge_loc = edge_bdry[i].getLocationIndex();

      hier::Box fill_box(pgeom->getBoundaryFillBox(edge_bdry[i],
                            interior,
                            gcw_to_fill));
      const hier::Index& ibeg(fill_box.lower());
      const hier::Index& iend(fill_box.upper());

      SAMRAI_F77_FUNC(getcartedgebdry3d, GETCARTEDGEBDRY3D) (
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         ibeg(0), iend(0),
         ibeg(1), iend(1),
         ibeg(2), iend(2),
         ghost_cells(0), ghost_cells(1), ghost_cells(2),
         dx,
         bedge_loc,
         bdry_edge_conds[bedge_loc],
         &bdry_face_values[0],
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
 *    bdry_node_conds ...... vector of boundary conditions for patch nodes
 *    bdry_face_values ..... vector of boundary values for faces
 *                           (this must be consistent with boundary
 *                           condition types)
 */

void SkeletonBoundaryUtilities3::fillNodeBoundaryData(
   const std::string& varname,
   std::shared_ptr<pdat::CellData<double> >& vardata,
   const hier::Patch& patch,
   const hier::IntVector& ghost_fill_width,
   const std::vector<int>& bdry_node_conds,
   const std::vector<double>& bdry_face_values)
{
   NULL_USE(varname);

   TBOX_ASSERT(vardata);
   TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_3D_NODES);
   TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) ==
      NUM_3D_FACES * (vardata->getDepth()));

   if (!s_fortran_constants_stuffed) {
      stuff3dBdryFortConst();
   }

   const std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());
   //const double* dx = pgeom->getDx();
   const double dx[3] = { 0., 0., 0. };

   const hier::Box& interior(patch.getBox());
   const hier::Index& ifirst(interior.lower());
   const hier::Index& ilast(interior.upper());

   const hier::IntVector& ghost_cells = vardata->getGhostCellWidth();

   hier::IntVector gcw_to_fill = hier::IntVector::min(ghost_cells,
         ghost_fill_width);

   const std::vector<hier::BoundaryBox>& node_bdry =
      pgeom->getCodimensionBoundaries(Bdry::NODE3D);
   for (int i = 0; i < static_cast<int>(node_bdry.size()); ++i) {
      TBOX_ASSERT(node_bdry[i].getBoundaryType() == Bdry::NODE3D);

      int bnode_loc = node_bdry[i].getLocationIndex();

      hier::Box fill_box(pgeom->getBoundaryFillBox(node_bdry[i],
                            interior,
                            gcw_to_fill));
      const hier::Index& ibeg(fill_box.lower());
      const hier::Index& iend(fill_box.upper());

      SAMRAI_F77_FUNC(getcartnodebdry3d, GETCARTNODEBDRY3D) (
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
         ifirst(2), ilast(2),
         ibeg(0), iend(0),
         ibeg(1), iend(1),
         ibeg(2), iend(2),
         ghost_cells(0), ghost_cells(1), ghost_cells(2),
         dx,
         bnode_loc,
         bdry_node_conds[bnode_loc],
         &bdry_face_values[0],
         vardata->getPointer(),
         vardata->getDepth());

   }

}

/*
 * Function that returns the integer face boundary location
 * corresponding to the given edge location and edge boundary
 * condition.
 *
 * If the edge boundary condition type or edge location are unknown,
 * or the boundary condition type is inconsistant with the edge location
 * an error results.
 */

int SkeletonBoundaryUtilities3::getFaceLocationForEdgeBdry(
   int edge_loc,
   int edge_btype)
{

   int ret_face = -1;

   switch (edge_btype) {
      case BdryCond::XFLOW:
      case BdryCond::XREFLECT:
      case BdryCond::XDIRICHLET:
      case BdryCond::XNEUMANN:
      {
         if (edge_loc == EdgeBdyLoc3D::XLO_ZLO || edge_loc == EdgeBdyLoc3D::XLO_ZHI ||
             edge_loc == EdgeBdyLoc3D::XLO_YLO || edge_loc == EdgeBdyLoc3D::XLO_YHI) {
            ret_face = BdryLoc::XLO;
         } else if (edge_loc == EdgeBdyLoc3D::XHI_ZLO ||
                    edge_loc == EdgeBdyLoc3D::XHI_ZHI ||
                    edge_loc == EdgeBdyLoc3D::XHI_YLO ||
                    edge_loc == EdgeBdyLoc3D::XHI_YHI) {
            ret_face = BdryLoc::XHI;
         }
         break;
      }
      case BdryCond::YFLOW:
      case BdryCond::YREFLECT:
      case BdryCond::YDIRICHLET:
      case BdryCond::YNEUMANN:
      {
         if (edge_loc == EdgeBdyLoc3D::YLO_ZLO || edge_loc == EdgeBdyLoc3D::YLO_ZHI ||
             edge_loc == EdgeBdyLoc3D::XLO_YLO || edge_loc == EdgeBdyLoc3D::XHI_YLO) {
            ret_face = BdryLoc::YLO;
         } else if (edge_loc == EdgeBdyLoc3D::YHI_ZLO ||
                    edge_loc == EdgeBdyLoc3D::YHI_ZHI ||
                    edge_loc == EdgeBdyLoc3D::XLO_YHI ||
                    edge_loc == EdgeBdyLoc3D::XHI_YHI) {
            ret_face = BdryLoc::YHI;
         }
         break;
      }
      case BdryCond::ZFLOW:
      case BdryCond::ZREFLECT:
      case BdryCond::ZDIRICHLET:
      case BdryCond::ZNEUMANN:
      {
         if (edge_loc == EdgeBdyLoc3D::YLO_ZLO || edge_loc == EdgeBdyLoc3D::YHI_ZLO ||
             edge_loc == EdgeBdyLoc3D::XLO_ZLO || edge_loc == EdgeBdyLoc3D::XHI_ZLO) {
            ret_face = BdryLoc::ZLO;
         } else if (edge_loc == EdgeBdyLoc3D::YLO_ZHI ||
                    edge_loc == EdgeBdyLoc3D::YHI_ZHI ||
                    edge_loc == EdgeBdyLoc3D::XLO_ZHI ||
                    edge_loc == EdgeBdyLoc3D::XHI_ZHI) {
            ret_face = BdryLoc::ZHI;
         }
         break;
      }
      default: {
         TBOX_ERROR("Unknown edge boundary condition type = "
            << edge_btype << " passed to \n"
            << "SkeletonBoundaryUtilities3::getFaceLocationForEdgeBdry"
            << std::endl);
      }
   }

   if (ret_face == -1) {
      TBOX_ERROR("Edge boundary condition type = "
         << edge_btype << " and edge location = " << edge_loc
         << "\n passed to "
         << "SkeletonBoundaryUtilities3::getFaceLocationForEdgeBdry"
         << " are inconsistant." << std::endl);
   }

   return ret_face;

}

/*
 * Function that returns the integer face boundary location
 * corresponding to the given node location and node boundary
 * condition.
 *
 * If the node boundary condition type or node location are unknown,
 * or the boundary condition type is inconsistant with the node location
 * an error results.
 */

int SkeletonBoundaryUtilities3::getFaceLocationForNodeBdry(
   int node_loc,
   int node_btype)
{

   int ret_face = -1;

   switch (node_btype) {
      case BdryCond::XFLOW:
      case BdryCond::XREFLECT:
      case BdryCond::XDIRICHLET:
      case BdryCond::XNEUMANN:
      {
         if (node_loc == NodeBdyLoc3D::XLO_YLO_ZLO ||
             node_loc == NodeBdyLoc3D::XLO_YHI_ZLO ||
             node_loc == NodeBdyLoc3D::XLO_YLO_ZHI ||
             node_loc == NodeBdyLoc3D::XLO_YHI_ZHI) {
            ret_face = BdryLoc::XLO;
         } else {
            ret_face = BdryLoc::XHI;
         }
         break;
      }
      case BdryCond::YFLOW:
      case BdryCond::YREFLECT:
      case BdryCond::YDIRICHLET:
      case BdryCond::YNEUMANN:
      {
         if (node_loc == NodeBdyLoc3D::XLO_YLO_ZLO ||
             node_loc == NodeBdyLoc3D::XHI_YLO_ZLO ||
             node_loc == NodeBdyLoc3D::XLO_YLO_ZHI ||
             node_loc == NodeBdyLoc3D::XHI_YLO_ZHI) {
            ret_face = BdryLoc::YLO;
         } else {
            ret_face = BdryLoc::YHI;
         }
         break;
      }
      case BdryCond::ZFLOW:
      case BdryCond::ZREFLECT:
      case BdryCond::ZDIRICHLET:
      case BdryCond::ZNEUMANN:
      {
         if (node_loc == NodeBdyLoc3D::XLO_YLO_ZLO ||
             node_loc == NodeBdyLoc3D::XHI_YLO_ZLO ||
             node_loc == NodeBdyLoc3D::XLO_YHI_ZLO ||
             node_loc == NodeBdyLoc3D::XHI_YHI_ZLO) {
            ret_face = BdryLoc::ZLO;
         } else {
            ret_face = BdryLoc::ZHI;
         }
         break;
      }
      default: {
         TBOX_ERROR("Unknown node boundary condition type = "
            << node_btype << " passed to \n"
            << "SkeletonBoundaryUtilities3::getFaceLocationForNodeBdry"
            << std::endl);
      }
   }

   if (ret_face == -1) {
      TBOX_ERROR("Node boundary condition type = "
         << node_btype << " and node location = " << node_loc
         << "\n passed to "
         << "SkeletonBoundaryUtilities3::getFaceLocationForNodeBdry"
         << " are inconsistant." << std::endl);
   }

   return ret_face;

}

/*
 * Function to check 3D boundary data filling.  Arguments are:
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

int SkeletonBoundaryUtilities3::checkBdryData(
   const std::string& varname,
   const hier::Patch& patch,
   int data_id,
   int depth,
   const hier::IntVector& gcw_to_check,
   const hier::BoundaryBox& bbox,
   int bcase,
   double bstate)
{
   TBOX_ASSERT(!varname.empty());
   TBOX_ASSERT(data_id >= 0);
   TBOX_ASSERT(depth >= 0);

   int num_bad_values = 0;

   int btype = bbox.getBoundaryType();
   int bloc = bbox.getLocationIndex();

   std::shared_ptr<hier::PatchGeometry> pgeom(
      patch.getPatchGeometry());

   std::shared_ptr<pdat::CellData<double> > vardata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(data_id)));
   TBOX_ASSERT(vardata);

   std::string bdry_type_str;
   if (btype == Bdry::FACE3D) {
      bdry_type_str = "FACE";
   } else if (btype == Bdry::EDGE3D) {
      bdry_type_str = "EDGE";
   } else if (btype == Bdry::NODE3D) {
      bdry_type_str = "NODE";
   } else {
      TBOX_ERROR(
         "Unknown btype " << btype
                          << " passed to SkeletonBoundaryUtilities3::checkBdryData()! "
                          << std::endl);
   }

   tbox::plog << "\n\nCHECKING 3D " << bdry_type_str << " BDRY DATA..." << std::endl;
   tbox::plog << "varname = " << varname << " : depth = " << depth << std::endl;
   tbox::plog << "bbox = " << bbox.getBox() << std::endl;
   tbox::plog << "btype, bloc, bcase = "
              << btype << ", = " << bloc << ", = " << bcase << std::endl;

   tbox::Dimension::dir_t idir;
   double valfact = 0.0, constval = 0.0, dxfact = 0.0;
   int offsign;

   get3dBdryDirectionCheckValues(idir, offsign,
      btype, bloc, bcase);

   if (btype == Bdry::FACE3D) {

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
      } else {
         TBOX_ERROR(
            "Unknown bcase " << bcase
                             << " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                             << "\n for " << bdry_type_str
                             << " at location " << bloc << std::endl);
      }

   } else if (btype == Bdry::EDGE3D) {

      if (bcase == BdryCond::XFLOW || bcase == BdryCond::YFLOW ||
          bcase == BdryCond::ZFLOW) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XREFLECT || bcase == BdryCond::YREFLECT ||
                 bcase == BdryCond::ZREFLECT) {
         valfact = -1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XDIRICHLET || bcase == BdryCond::YDIRICHLET ||
                 bcase == BdryCond::ZDIRICHLET) {
         valfact = 0.0;
         constval = bstate;
         dxfact = 0.0;
      } else {
         TBOX_ERROR(
            "Unknown bcase " << bcase
                             << " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                             << "\n for " << bdry_type_str
                             << " at location " << bloc << std::endl);
      }

   } else if (btype == Bdry::NODE3D) {

      if (bcase == BdryCond::XFLOW || bcase == BdryCond::YFLOW ||
          bcase == BdryCond::ZFLOW) {
         valfact = 1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XREFLECT || bcase == BdryCond::YREFLECT ||
                 bcase == BdryCond::ZREFLECT) {
         valfact = -1.0;
         constval = 0.0;
         dxfact = 0.0;
      } else if (bcase == BdryCond::XDIRICHLET || bcase == BdryCond::YDIRICHLET ||
                 bcase == BdryCond::ZDIRICHLET) {
         valfact = 0.0;
         constval = bstate;
         dxfact = 0.0;
      } else {
         TBOX_ERROR(
            "Unknown bcase " << bcase
                             << " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                             << "\n for " << bdry_type_str
                             << " at location " << bloc << std::endl);
      }

   }

   hier::Box gbox_to_check(vardata->getGhostBox() * pgeom->getBoundaryFillBox(
                              bbox,
                              patch.getBox(),
                              gcw_to_check));

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
 * Private function to read 3D face boundary data from input database.
 */

void SkeletonBoundaryUtilities3::read3dBdryFaces(
   BoundaryUtilityStrategy* bdry_strategy,
   std::shared_ptr<tbox::Database> input_db,
   std::vector<int>& face_conds,
   const hier::IntVector& periodic)
{
   TBOX_ASSERT(bdry_strategy != 0);
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);

   int num_per_dirs = 0;
   for (int id = 0; id < 3; ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < 3) { // face boundary input required

      for (int s = 0; s < NUM_3D_FACES; ++s) {

         std::string bdry_loc_str;
         switch (s) {
            case BdryLoc::XLO: { bdry_loc_str = "boundary_face_xlo";
                                 break;
            }
            case BdryLoc::XHI: { bdry_loc_str = "boundary_face_xhi";
                                 break;
            }
            case BdryLoc::YLO: { bdry_loc_str = "boundary_face_ylo";
                                 break;
            }
            case BdryLoc::YHI: { bdry_loc_str = "boundary_face_yhi";
                                 break;
            }
            case BdryLoc::ZLO: { bdry_loc_str = "boundary_face_zlo";
                                 break;
            }
            case BdryLoc::ZHI: { bdry_loc_str = "boundary_face_zhi";
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
            } else if (periodic(2) && (s == BdryLoc::ZLO || s == BdryLoc::ZHI)) {
               need_data_read = false;
            }
         }

         if (need_data_read) {
            if (input_db->keyExists(bdry_loc_str)) {
               std::shared_ptr<tbox::Database> bdry_loc_db(
                  input_db->getDatabase(bdry_loc_str));
               if (bdry_loc_db) {
                  if (bdry_loc_db->keyExists("boundary_condition")) {
                     std::string bdry_cond_str =
                        bdry_loc_db->getString("boundary_condition");
                     if (bdry_cond_str == "FLOW") {
                        face_conds[s] = BdryCond::FLOW;
                     } else if (bdry_cond_str == "REFLECT") {
                        face_conds[s] = BdryCond::REFLECT;
                     } else if (bdry_cond_str == "DIRICHLET") {
                        face_conds[s] = BdryCond::DIRICHLET;
                        bdry_strategy->
                        readDirichletBoundaryDataEntry(bdry_loc_db,
                           bdry_loc_str,
                           s);
                     } else if (bdry_cond_str == "NEUMANN") {
                        face_conds[s] = BdryCond::NEUMANN;
                        bdry_strategy->
                        readNeumannBoundaryDataEntry(bdry_loc_db,
                           bdry_loc_str,
                           s);
                     } else {
                        TBOX_ERROR("Unknown face boundary string = "
                           << bdry_cond_str << " found in input." << std::endl);
                     }
                  } else {
                     TBOX_ERROR("'boundary_condition' entry missing from "
                        << bdry_loc_str << " input database." << std::endl);
                  }
               }
            } else {
               TBOX_ERROR(bdry_loc_str
                  << " database entry not found in input." << std::endl);
            }
         } // if (need_data_read)

      } // for (int s = 0 ...

   } // if (num_per_dirs < 3)

}

/*
 * Private function to read 3D edge boundary data from input database.
 */

void SkeletonBoundaryUtilities3::read3dBdryEdges(
   std::shared_ptr<tbox::Database> input_db,
   const std::vector<int>& face_conds,
   std::vector<int>& edge_conds,
   const hier::IntVector& periodic)
{
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
   TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);

   int num_per_dirs = 0;
   for (int id = 0; id < 3; ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < 2) {  // edge boundary input required

      for (int s = 0; s < NUM_3D_EDGES; ++s) {

         std::string bdry_loc_str;
         switch (s) {
            case EdgeBdyLoc3D::YLO_ZLO: {
               bdry_loc_str = "boundary_edge_ylo_zlo";
               break;
            }
            case EdgeBdyLoc3D::YHI_ZLO: {
               bdry_loc_str = "boundary_edge_yhi_zlo";
               break;
            }
            case EdgeBdyLoc3D::YLO_ZHI: {
               bdry_loc_str = "boundary_edge_ylo_zhi";
               break;
            }
            case EdgeBdyLoc3D::YHI_ZHI: {
               bdry_loc_str = "boundary_edge_yhi_zhi";
               break;
            }
            case EdgeBdyLoc3D::XLO_ZLO: {
               bdry_loc_str = "boundary_edge_xlo_zlo";
               break;
            }
            case EdgeBdyLoc3D::XLO_ZHI: {
               bdry_loc_str = "boundary_edge_xlo_zhi";
               break;
            }
            case EdgeBdyLoc3D::XHI_ZLO: {
               bdry_loc_str = "boundary_edge_xhi_zlo";
               break;
            }
            case EdgeBdyLoc3D::XHI_ZHI: {
               bdry_loc_str = "boundary_edge_xhi_zhi";
               break;
            }
            case EdgeBdyLoc3D::XLO_YLO: {
               bdry_loc_str = "boundary_edge_xlo_ylo";
               break;
            }
            case EdgeBdyLoc3D::XHI_YLO: {
               bdry_loc_str = "boundary_edge_xhi_ylo";
               break;
            }
            case EdgeBdyLoc3D::XLO_YHI: {
               bdry_loc_str = "boundary_edge_xlo_yhi";
               break;
            }
            case EdgeBdyLoc3D::XHI_YHI: {
               bdry_loc_str = "boundary_edge_xhi_yhi";
               break;
            }
            default: NULL_STATEMENT;
         }

         bool need_data_read = false;
         if (num_per_dirs == 0) {
            need_data_read = true;
         } else if (periodic(0) &&
                    (s == EdgeBdyLoc3D::YLO_ZLO || s == EdgeBdyLoc3D::YHI_ZLO ||
                     s == EdgeBdyLoc3D::YLO_ZHI || s == EdgeBdyLoc3D::YHI_ZHI)) {
            need_data_read = true;
         } else if (periodic(1) &&
                    (s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::XLO_ZHI ||
                     s == EdgeBdyLoc3D::XHI_ZLO || s == EdgeBdyLoc3D::XHI_ZHI)) {
            need_data_read = true;
         } else if (periodic(2) &&
                    (s == EdgeBdyLoc3D::XLO_YLO || s == EdgeBdyLoc3D::XHI_YLO ||
                     s == EdgeBdyLoc3D::XLO_YHI || s == EdgeBdyLoc3D::XHI_YHI)) {
            need_data_read = true;
         }

         if (need_data_read) {
            if (input_db->keyExists(bdry_loc_str)) {
               std::shared_ptr<tbox::Database> bdry_loc_db(
                  input_db->getDatabase(bdry_loc_str));
               if (bdry_loc_db) {
                  if (bdry_loc_db->keyExists("boundary_condition")) {
                     std::string bdry_cond_str =
                        bdry_loc_db->getString("boundary_condition");
                     if (bdry_cond_str == "XFLOW") {
                        edge_conds[s] = BdryCond::XFLOW;
                     } else if (bdry_cond_str == "YFLOW") {
                        edge_conds[s] = BdryCond::YFLOW;
                     } else if (bdry_cond_str == "ZFLOW") {
                        edge_conds[s] = BdryCond::ZFLOW;
                     } else if (bdry_cond_str == "XREFLECT") {
                        edge_conds[s] = BdryCond::XREFLECT;
                     } else if (bdry_cond_str == "YREFLECT") {
                        edge_conds[s] = BdryCond::YREFLECT;
                     } else if (bdry_cond_str == "ZREFLECT") {
                        edge_conds[s] = BdryCond::ZREFLECT;
                     } else if (bdry_cond_str == "XDIRICHLET") {
                        edge_conds[s] = BdryCond::XDIRICHLET;
                     } else if (bdry_cond_str == "YDIRICHLET") {
                        edge_conds[s] = BdryCond::YDIRICHLET;
                     } else if (bdry_cond_str == "ZDIRICHLET") {
                        edge_conds[s] = BdryCond::ZDIRICHLET;
                     } else if (bdry_cond_str == "XNEUMANN") {
                        edge_conds[s] = BdryCond::XNEUMANN;
                     } else if (bdry_cond_str == "YNEUMANN") {
                        edge_conds[s] = BdryCond::YNEUMANN;
                     } else if (bdry_cond_str == "ZNEUMANN") {
                        edge_conds[s] = BdryCond::ZNEUMANN;
                     } else {
                        TBOX_ERROR("Unknown edge boundary string = "
                           << bdry_cond_str << " found in input." << std::endl);
                     }

                     bool ambiguous_type = false;
                     if (bdry_cond_str == "XFLOW" ||
                         bdry_cond_str == "XREFLECT" ||
                         bdry_cond_str == "XDIRICHLET" ||
                         bdry_cond_str == "XNEUMANN") {
                        if (s == EdgeBdyLoc3D::YLO_ZLO || s == EdgeBdyLoc3D::YHI_ZLO ||
                            s == EdgeBdyLoc3D::YLO_ZHI || s == EdgeBdyLoc3D::YHI_ZHI) {
                           ambiguous_type = true;
                        }
                     } else if (bdry_cond_str == "YFLOW" ||
                                bdry_cond_str == "YREFLECT" ||
                                bdry_cond_str == "YDIRICHLET" ||
                                bdry_cond_str == "YNEUMANN") {
                        if (s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::XLO_ZHI ||
                            s == EdgeBdyLoc3D::XHI_ZLO || s == EdgeBdyLoc3D::XHI_ZHI) {
                           ambiguous_type = true;
                        }
                     } else if (bdry_cond_str == "ZFLOW" ||
                                bdry_cond_str == "ZREFLECT" ||
                                bdry_cond_str == "ZDIRICHLET" ||
                                bdry_cond_str == "ZNEUMANN") {
                        if (s == EdgeBdyLoc3D::XLO_YLO || s == EdgeBdyLoc3D::XHI_YLO ||
                            s == EdgeBdyLoc3D::XLO_YHI || s == EdgeBdyLoc3D::XHI_YHI) {
                           ambiguous_type = true;
                        }
                     }
                     if (ambiguous_type) {
                        TBOX_ERROR("Ambiguous bdry condition "
                           << bdry_cond_str
                           << " found for " << bdry_loc_str << std::endl);
                     }

                     std::string proper_face;
                     std::string proper_face_data;
                     bool no_face_data_found = false;
                     if (bdry_cond_str == "XFLOW" ||
                         bdry_cond_str == "XDIRICHLET" ||
                         bdry_cond_str == "XNEUMANN" ||
                         bdry_cond_str == "XREFLECT") {
                        if (s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::XLO_ZHI ||
                            s == EdgeBdyLoc3D::XLO_YLO || s == EdgeBdyLoc3D::XLO_YHI) {
                           proper_face = "XLO";
                           if (bdry_cond_str == "XFLOW" &&
                               face_conds[BdryLoc::XLO] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "XDIRICHLET" &&
                               face_conds[BdryLoc::XLO] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "XNEUMANN" &&
                               face_conds[BdryLoc::XLO] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "XREFLECT" &&
                               face_conds[BdryLoc::XLO] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "XHI";
                           if (bdry_cond_str == "XFLOW" &&
                               face_conds[BdryLoc::XHI] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "XDIRICHLET" &&
                               face_conds[BdryLoc::XHI] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "XNEUMANN" &&
                               face_conds[BdryLoc::XHI] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "XREFLECT" &&
                               face_conds[BdryLoc::XHI] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     } else if (bdry_cond_str == "YFLOW" ||
                                bdry_cond_str == "YDIRICHLET" ||
                                bdry_cond_str == "YNEUMANN" ||
                                bdry_cond_str == "YREFLECT") {
                        if (s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::YLO_ZHI ||
                            s == EdgeBdyLoc3D::XLO_YLO || s == EdgeBdyLoc3D::XHI_YLO) {
                           proper_face = "YLO";
                           if (bdry_cond_str == "YFLOW" &&
                               face_conds[BdryLoc::YLO] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "YDIRICHLET" &&
                               face_conds[BdryLoc::YLO] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "YNEUMANN" &&
                               face_conds[BdryLoc::YLO] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "YREFLECT" &&
                               face_conds[BdryLoc::YLO] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "YHI";
                           if (bdry_cond_str == "YFLOW" &&
                               face_conds[BdryLoc::YHI] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "YDIRICHLET" &&
                               face_conds[BdryLoc::YHI] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "YNEUMANN" &&
                               face_conds[BdryLoc::YHI] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "YREFLECT" &&
                               face_conds[BdryLoc::YHI] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     } else if (bdry_cond_str == "ZFLOW" ||
                                bdry_cond_str == "ZDIRICHLET" ||
                                bdry_cond_str == "ZNEUMANN" ||
                                bdry_cond_str == "ZREFLECT") {
                        if (s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::YHI_ZLO ||
                            s == EdgeBdyLoc3D::XLO_ZLO || s == EdgeBdyLoc3D::XHI_ZLO) {
                           proper_face = "ZLO";
                           if (bdry_cond_str == "ZFLOW" &&
                               face_conds[BdryLoc::ZLO] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "ZDIRICHLET" &&
                               face_conds[BdryLoc::ZLO] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "ZNEUMANN" &&
                               face_conds[BdryLoc::ZLO] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "ZREFLECT" &&
                               face_conds[BdryLoc::ZLO] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        } else {
                           proper_face = "ZHI";
                           if (bdry_cond_str == "ZFLOW" &&
                               face_conds[BdryLoc::ZHI] != BdryCond::FLOW) {
                              no_face_data_found = true;
                              proper_face_data = "FLOW";
                           }
                           if (bdry_cond_str == "ZDIRICHLET" &&
                               face_conds[BdryLoc::ZHI] != BdryCond::DIRICHLET) {
                              no_face_data_found = true;
                              proper_face_data = "DIRICHLET";
                           }
                           if (bdry_cond_str == "ZNEUMANN" &&
                               face_conds[BdryLoc::ZHI] != BdryCond::NEUMANN) {
                              no_face_data_found = true;
                              proper_face_data = "NEUMANN";
                           }
                           if (bdry_cond_str == "ZREFLECT" &&
                               face_conds[BdryLoc::ZHI] != BdryCond::REFLECT) {
                              no_face_data_found = true;
                              proper_face_data = "REFLECT";
                           }
                        }
                     }
                     if (no_face_data_found) {
                        TBOX_ERROR(
                           "Bdry condition " << bdry_cond_str
                                             << " found for "
                                             << bdry_loc_str
                                             << "\n but no "
                                             << proper_face_data
                                             << " data found for face "
                                             << proper_face << std::endl);
                     }
                  } else {
                     TBOX_ERROR("'boundary_condition' entry missing from "
                        << bdry_loc_str << " input database." << std::endl);
                  }
               }
            } else {
               TBOX_ERROR(bdry_loc_str
                  << " database entry not found in input." << std::endl);
            }

         } // if (need_data_read)

      } // for (int s = 0 ...

   } // if (num_per_dirs < 2)

}

/*
 * Private function to read 3D node boundary data from input database.
 */

void SkeletonBoundaryUtilities3::read3dBdryNodes(
   std::shared_ptr<tbox::Database> input_db,
   const std::vector<int>& face_conds,
   std::vector<int>& node_conds,
   const hier::IntVector& periodic)
{
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
   TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);

   int num_per_dirs = 0;
   for (int id = 0; id < 3; ++id) {
      if (periodic(id)) ++num_per_dirs;
   }

   if (num_per_dirs < 1) { // node boundary data required

      for (int s = 0; s < NUM_3D_NODES; ++s) {

         std::string bdry_loc_str;
         switch (s) {
            case NodeBdyLoc3D::XLO_YLO_ZLO: {
               bdry_loc_str = "boundary_node_xlo_ylo_zlo";
               break;
            }
            case NodeBdyLoc3D::XHI_YLO_ZLO: {
               bdry_loc_str = "boundary_node_xhi_ylo_zlo";
               break;
            }
            case NodeBdyLoc3D::XLO_YHI_ZLO: {
               bdry_loc_str = "boundary_node_xlo_yhi_zlo";
               break;
            }
            case NodeBdyLoc3D::XHI_YHI_ZLO: {
               bdry_loc_str = "boundary_node_xhi_yhi_zlo";
               break;
            }
            case NodeBdyLoc3D::XLO_YLO_ZHI: {
               bdry_loc_str = "boundary_node_xlo_ylo_zhi";
               break;
            }
            case NodeBdyLoc3D::XHI_YLO_ZHI: {
               bdry_loc_str = "boundary_node_xhi_ylo_zhi";
               break;
            }
            case NodeBdyLoc3D::XLO_YHI_ZHI: {
               bdry_loc_str = "boundary_node_xlo_yhi_zhi";
               break;
            }
            case NodeBdyLoc3D::XHI_YHI_ZHI: {
               bdry_loc_str = "boundary_node_xhi_yhi_zhi";
               break;
            }
            default: NULL_STATEMENT;
         }

         if (input_db->keyExists(bdry_loc_str)) {
            std::shared_ptr<tbox::Database> bdry_loc_db(
               input_db->getDatabase(bdry_loc_str));
            if (bdry_loc_db) {
               if (bdry_loc_db->keyExists("boundary_condition")) {
                  std::string bdry_cond_str =
                     bdry_loc_db->getString("boundary_condition");
                  if (bdry_cond_str == "XFLOW") {
                     node_conds[s] = BdryCond::XFLOW;
                  } else if (bdry_cond_str == "YFLOW") {
                     node_conds[s] = BdryCond::YFLOW;
                  } else if (bdry_cond_str == "ZFLOW") {
                     node_conds[s] = BdryCond::ZFLOW;
                  } else if (bdry_cond_str == "XREFLECT") {
                     node_conds[s] = BdryCond::XREFLECT;
                  } else if (bdry_cond_str == "YREFLECT") {
                     node_conds[s] = BdryCond::YREFLECT;
                  } else if (bdry_cond_str == "ZREFLECT") {
                     node_conds[s] = BdryCond::ZREFLECT;
                  } else if (bdry_cond_str == "XDIRICHLET") {
                     node_conds[s] = BdryCond::XDIRICHLET;
                  } else if (bdry_cond_str == "YDIRICHLET") {
                     node_conds[s] = BdryCond::YDIRICHLET;
                  } else if (bdry_cond_str == "ZDIRICHLET") {
                     node_conds[s] = BdryCond::ZDIRICHLET;
                  } else if (bdry_cond_str == "XNEUMANN") {
                     node_conds[s] = BdryCond::XNEUMANN;
                  } else if (bdry_cond_str == "YNEUMANN") {
                     node_conds[s] = BdryCond::YNEUMANN;
                  } else if (bdry_cond_str == "ZNEUMANN") {
                     node_conds[s] = BdryCond::ZNEUMANN;
                  } else {
                     TBOX_ERROR("Unknown node boundary string = "
                        << bdry_cond_str << " found in input." << std::endl);
                  }

                  std::string proper_face;
                  std::string proper_face_data;
                  bool no_face_data_found = false;
                  if (bdry_cond_str == "XFLOW" ||
                      bdry_cond_str == "XDIRICHLET" ||
                      bdry_cond_str == "XNEUMANN" ||
                      bdry_cond_str == "XREFLECT") {
                     if (s == NodeBdyLoc3D::XLO_YLO_ZLO ||
                         s == NodeBdyLoc3D::XLO_YHI_ZLO ||
                         s == NodeBdyLoc3D::XLO_YLO_ZHI ||
                         s == NodeBdyLoc3D::XLO_YHI_ZHI) {
                        proper_face = "XLO";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[BdryLoc::XLO] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[BdryLoc::XLO] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[BdryLoc::XLO] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[BdryLoc::XLO] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "XHI";
                        if (bdry_cond_str == "XFLOW" &&
                            face_conds[BdryLoc::XHI] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "XDIRICHLET" &&
                            face_conds[BdryLoc::XHI] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "XNEUMANN" &&
                            face_conds[BdryLoc::XHI] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "XREFLECT" &&
                            face_conds[BdryLoc::XHI] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  } else if (bdry_cond_str == "YFLOW" ||
                             bdry_cond_str == "YDIRICHLET" ||
                             bdry_cond_str == "YNEUMANN" ||
                             bdry_cond_str == "YREFLECT") {
                     if (s == NodeBdyLoc3D::XLO_YLO_ZLO ||
                         s == NodeBdyLoc3D::XHI_YLO_ZLO ||
                         s == NodeBdyLoc3D::XLO_YLO_ZHI ||
                         s == NodeBdyLoc3D::XHI_YLO_ZHI) {
                        proper_face = "YLO";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[BdryLoc::YLO] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[BdryLoc::YLO] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[BdryLoc::YLO] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[BdryLoc::YLO] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "YHI";
                        if (bdry_cond_str == "YFLOW" &&
                            face_conds[BdryLoc::YHI] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "YDIRICHLET" &&
                            face_conds[BdryLoc::YHI] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "YNEUMANN" &&
                            face_conds[BdryLoc::YHI] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "YREFLECT" &&
                            face_conds[BdryLoc::YHI] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  } else if (bdry_cond_str == "ZFLOW" ||
                             bdry_cond_str == "ZDIRICHLET" ||
                             bdry_cond_str == "ZNEUMANN" ||
                             bdry_cond_str == "ZREFLECT") {
                     if (s == NodeBdyLoc3D::XLO_YLO_ZLO ||
                         s == NodeBdyLoc3D::XHI_YLO_ZLO ||
                         s == NodeBdyLoc3D::XLO_YHI_ZLO ||
                         s == NodeBdyLoc3D::XHI_YHI_ZLO) {
                        proper_face = "ZLO";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[BdryLoc::ZLO] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[BdryLoc::ZLO] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[BdryLoc::ZLO] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[BdryLoc::ZLO] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     } else {
                        proper_face = "ZHI";
                        if (bdry_cond_str == "ZFLOW" &&
                            face_conds[BdryLoc::ZHI] != BdryCond::FLOW) {
                           no_face_data_found = true;
                           proper_face_data = "FLOW";
                        }
                        if (bdry_cond_str == "ZDIRICHLET" &&
                            face_conds[BdryLoc::ZHI] != BdryCond::DIRICHLET) {
                           no_face_data_found = true;
                           proper_face_data = "DIRICHLET";
                        }
                        if (bdry_cond_str == "ZNEUMANN" &&
                            face_conds[BdryLoc::ZHI] != BdryCond::NEUMANN) {
                           no_face_data_found = true;
                           proper_face_data = "NEUMANN";
                        }
                        if (bdry_cond_str == "ZREFLECT" &&
                            face_conds[BdryLoc::ZHI] != BdryCond::REFLECT) {
                           no_face_data_found = true;
                           proper_face_data = "REFLECT";
                        }
                     }
                  }
                  if (no_face_data_found) {
                     TBOX_ERROR(
                        "Bdry condition " << bdry_cond_str
                                          << " found for "
                                          << bdry_loc_str
                                          << "\n but no "
                                          << proper_face_data
                                          << " data found for face "
                                          << proper_face << std::endl);
                  }

               } else {
                  TBOX_ERROR("'boundary_condition' entry missing from "
                     << bdry_loc_str << " input database." << std::endl);
               }
            }
         } else {
            TBOX_ERROR(bdry_loc_str
               << " database entry not found in input." << std::endl);
         }

      } // for (int s = 0 ...

   } // if (num_per_dirs < 1)

}

/*
 * Private function to get boundary orientation information for
 * 3D boundary condition checking.  Called from checkBdryData().
 */

void SkeletonBoundaryUtilities3::get3dBdryDirectionCheckValues(
   tbox::Dimension::dir_t& idir,
   int& offsign,
   int btype,
   int bloc,
   int bcase)
{

   std::string bdry_type_str;

   if (btype == Bdry::FACE3D) {

      bdry_type_str = "FACE";

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
      } else if (bloc == BdryLoc::ZLO || bloc == BdryLoc::ZHI) {
         idir = 2;
         if (bloc == BdryLoc::ZLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else {
         TBOX_ERROR(
            "Unknown boundary location " << bloc
                                         <<
            " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                                         << "\n for "
                                         << bdry_type_str << " boundary " << std::endl);
      }

   } else if (btype == Bdry::EDGE3D) {

      bdry_type_str = "EDGE";

      bool bad_case = false;
      if (bcase == BdryCond::XFLOW || bcase == BdryCond::XREFLECT ||
          bcase == BdryCond::XDIRICHLET || bcase == BdryCond::XNEUMANN) {
         idir = 0;
         if (bloc == EdgeBdyLoc3D::XLO_ZLO || bloc == EdgeBdyLoc3D::XLO_ZHI ||
             bloc == EdgeBdyLoc3D::XLO_YLO || bloc == EdgeBdyLoc3D::XLO_YHI) {
            offsign = -1;
         } else if (bloc == EdgeBdyLoc3D::XHI_ZLO || bloc == EdgeBdyLoc3D::XHI_ZHI ||
                    bloc == EdgeBdyLoc3D::XHI_YLO || bloc == EdgeBdyLoc3D::XHI_YHI) {
            offsign = 1;
         } else {
            bad_case = true;
         }
      } else if (bcase == BdryCond::YFLOW || bcase == BdryCond::YREFLECT ||
                 bcase == BdryCond::YDIRICHLET || bcase == BdryCond::YNEUMANN) {
         idir = 1;
         if (bloc == EdgeBdyLoc3D::YLO_ZLO || bloc == EdgeBdyLoc3D::YLO_ZHI ||
             bloc == EdgeBdyLoc3D::XLO_YLO || bloc == EdgeBdyLoc3D::XHI_YLO) {
            offsign = -1;
         } else if (bloc == EdgeBdyLoc3D::YHI_ZLO || bloc == EdgeBdyLoc3D::YHI_ZHI ||
                    bloc == EdgeBdyLoc3D::XLO_YHI || bloc == EdgeBdyLoc3D::XHI_YHI) {
            offsign = 1;
         } else {
            bad_case = true;
         }
      } else if (bcase == BdryCond::ZFLOW || bcase == BdryCond::ZREFLECT ||
                 bcase == BdryCond::ZDIRICHLET || bcase == BdryCond::ZNEUMANN) {
         idir = 2;
         if (bloc == EdgeBdyLoc3D::YLO_ZLO || bloc == EdgeBdyLoc3D::YHI_ZLO ||
             bloc == EdgeBdyLoc3D::XLO_ZLO || bloc == EdgeBdyLoc3D::XHI_ZLO) {
            offsign = -1;
         } else if (bloc == EdgeBdyLoc3D::YLO_ZHI || bloc == EdgeBdyLoc3D::YHI_ZHI ||
                    bloc == EdgeBdyLoc3D::XLO_ZHI || bloc == EdgeBdyLoc3D::XHI_ZHI) {
            offsign = 1;
         } else {
            bad_case = true;
         }
      }
      if (bad_case) {
         TBOX_ERROR(
            "Unknown or ambigous bcase " << bcase
                                         <<
            " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                                         << "\n for " << bdry_type_str
                                         << " at location " << bloc
                                         << std::endl);
      }

   } else if (btype == Bdry::NODE3D) {

      bdry_type_str = "NODE";

      if (bcase == BdryCond::XFLOW || bcase == BdryCond::XREFLECT ||
          bcase == BdryCond::XDIRICHLET || bcase == BdryCond::XNEUMANN) {
         idir = 0;
         if (bloc == NodeBdyLoc3D::XLO_YLO_ZLO || bloc == NodeBdyLoc3D::XLO_YHI_ZLO ||
             bloc == NodeBdyLoc3D::XLO_YLO_ZHI || bloc == NodeBdyLoc3D::XLO_YHI_ZHI) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == BdryCond::YFLOW || bcase == BdryCond::YREFLECT ||
                 bcase == BdryCond::YDIRICHLET || bcase == BdryCond::YNEUMANN) {
         idir = 1;
         if (bloc == NodeBdyLoc3D::XLO_YLO_ZLO || bloc == NodeBdyLoc3D::XHI_YLO_ZLO ||
             bloc == NodeBdyLoc3D::XLO_YLO_ZHI || bloc == NodeBdyLoc3D::XHI_YLO_ZHI) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      } else if (bcase == BdryCond::ZFLOW || bcase == BdryCond::ZREFLECT ||
                 bcase == BdryCond::ZDIRICHLET || bcase == BdryCond::ZNEUMANN) {
         idir = 2;
         if (bloc == NodeBdyLoc3D::XLO_YLO_ZLO || bloc == NodeBdyLoc3D::XHI_YLO_ZLO ||
             bloc == NodeBdyLoc3D::XLO_YHI_ZLO || bloc == NodeBdyLoc3D::XHI_YHI_ZLO) {
            offsign = -1;
         } else {
            offsign = 1;
         }
      }

   } else {
      TBOX_ERROR(
         "Unknown boundary type " << btype
                                  << " passed to SkeletonBoundaryUtilities3::checkBdryData()"
                                  << "\n for " << bdry_type_str
                                  << " at location " << bloc << std::endl);
   }

}

/*
 * Private function to stuff 3D boundary contants into Fortran common blocks
 */

void SkeletonBoundaryUtilities3::stuff3dBdryFortConst()
{
   SAMRAI_F77_FUNC(stufcartbdryloc3d, STUFCARTBDRYLOC3D) (
      BdryLoc::XLO, BdryLoc::XHI, BdryLoc::YLO, BdryLoc::YHI, BdryLoc::ZLO, BdryLoc::ZHI,
      EdgeBdyLoc3D::YLO_ZLO, EdgeBdyLoc3D::YHI_ZLO, EdgeBdyLoc3D::YLO_ZHI,
      EdgeBdyLoc3D::YHI_ZHI, EdgeBdyLoc3D::XLO_ZLO, EdgeBdyLoc3D::XLO_ZHI,
      EdgeBdyLoc3D::XHI_ZLO, EdgeBdyLoc3D::XHI_ZHI, EdgeBdyLoc3D::XLO_YLO,
      EdgeBdyLoc3D::XHI_YLO, EdgeBdyLoc3D::XLO_YHI, EdgeBdyLoc3D::XHI_YHI,
      NodeBdyLoc3D::XLO_YLO_ZLO, NodeBdyLoc3D::XHI_YLO_ZLO, NodeBdyLoc3D::XLO_YHI_ZLO,
      NodeBdyLoc3D::XHI_YHI_ZLO, NodeBdyLoc3D::XLO_YLO_ZHI, NodeBdyLoc3D::XHI_YLO_ZHI,
      NodeBdyLoc3D::XLO_YHI_ZHI, NodeBdyLoc3D::XHI_YHI_ZHI);
   SAMRAI_F77_FUNC(stufcartbdrycond3d, STUFCARTBDRYCOND3D) (
      BdryCond::FLOW,
      BdryCond::XFLOW, BdryCond::YFLOW, BdryCond::ZFLOW,
      BdryCond::REFLECT,
      BdryCond::XREFLECT, BdryCond::YREFLECT, BdryCond::ZREFLECT,
      BdryCond::DIRICHLET,
      BdryCond::XDIRICHLET, BdryCond::YDIRICHLET, BdryCond::ZDIRICHLET,
      BdryCond::NEUMANN,
      BdryCond::XNEUMANN, BdryCond::YNEUMANN, BdryCond::ZNEUMANN);
   s_fortran_constants_stuffed = true;
}
