/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility routines for manipulating Skeleton 3d boundary data
 *
 ************************************************************************/

#ifndef included_appu_SkeletonBoundaryUtilities3
#define included_appu_SkeletonBoundaryUtilities3

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <vector>
#include <memory>

using namespace SAMRAI;
using namespace appu;

/*!
 * @brief Class SkeletonBoundaryUtilities2 is a utility class that
 * simplifies the implementation of simple physical boundary data in
 * 3 spatial dimensions.  It contains routines for reading boundary data
 * information from input files, applying those boundary conditions,
 * and error checking boundary data.  These routines apply to the
 * case of cell-centered double data only.  One may use all of these
 * capabilities, or use the input reading, boundary setting, and error
 * checking routines independently.
 *
 * To use the boundary condition input reading capabilities, the format
 * of the input file section containing the boundary information must
 * be as described next.  Boundary face, node, and edge entries are only
 * required for those that are not filled automatically when periodic
 * conditions apply.
 *
 * The boundary condition for face "*" is provided in a section as follows:
 *
 * \verbatim
 *
 *    boundary_face_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *       // Any problem-specific boundary data read by user routines
 *       // is placed here...
 *    }
 *
 * Allowable face identifiers (i.e., values for "*") are:
 *       xlo, xhi, ylo, yhi, zlo, zhi
 * Supported face boundary_condition string values are:
 *       "FLOW", "REFLECT", "DIRICHLET", "NEUMANN"
 *
 * \endverbatim
 *
 * The boundary condition for edge "*" is provided in a section as follows:
 *
 * \verbatim
 *
 *    boundary_edge_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 *
 * Allowable edge identifiers (i.e., values for "*") are:
 *       ylo_zlo, yhi_zlo, ylo_zhi, yhi_zhi,
 *       xlo_zlo, xlo_zhi, xhi_zlo, xhi_zhi,
 *       xlo_ylo, xhi_ylo, xlo_yhi, xhi_yhi
 * Supported edge boundary_condition string values are:
 *       "XFLOW", "YFLOW", "ZFLOW",
 *       "XREFLECT", "YREFLECT", "ZREFLECT",
 *       "XDIRICHLET", "YDIRICHLET", "ZDIRICHLET"
 *       "XNEUMANN", "YNEUMANN", "ZNEUMANN"
 *
 * \endverbatim
 *
 * Note that edge conditions must be consistent with adjacent face conditions.
 *
 * The boundary condition for node "*" is provided in a section as follows:
 *
 * \verbatim
 *
 *    boundary_node_* {
 *       boundary_condition  = ...  // boundary condition string identifier
 *    }
 *
 * Allowable node identifiers (i.e., values for "*") are:
 *       ylo_zlo, yhi_zlo, ylo_zhi, yhi_zhi,
 *       xlo_zlo, xlo_zhi, xhi_zlo, xhi_zhi,
 *       xlo_ylo, xhi_ylo, xlo_yhi, xhi_yhi
 * Supported node boundary_condition values are:
 *       "XFLOW", "YFLOW", "ZFLOW",
 *       "XREFLECT", "YREFLECT", "ZREFLECT",
 *       "XDIRICHLET", "YDIRICHLET", "ZDIRICHLET"
 *       "XNEUMANN", "YNEUMANN", "ZNEUMANN"
 *
 * \endverbatim
 *
 * Note that node conditions must be consistent with adjacent face conditions.
 *
 * See the include file SkeletonBoundaryDefines.h for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.  If you choose to use the input reading
 * capabilities only and write your own boundary condition routines in
 * FORTRAN, you should note that the integer constants for the various
 * boundary condition types and locations are automatically "stuffed" into
 * FORTRAN common blocks.  This avoids potential problems with
 * inconsistencies between C++ and FORTRAN usage.  Please see the
 * FORTRAN include file cartbdryparams3d.i for details.
 *
 * @see appu::BoundaryUtilityStrategy3
 */

struct SkeletonBoundaryUtilities3 {
public:
   /*!
    * Function to read 3d boundary data from input database.
    * The integer boundary condition types are placed in the integer
    * arrays supplied by the caller (typically, the concrete
    * BoundaryUtilityStrategy object provided).  When DIRICHLET or
    * NEUMANN conditions are specified, control is passed to the
    * BoundaryUtilityStrategy to read the boundary state data specific to
    * the problem.
    *
    * Errors will be reported and the program will abort whenever necessary
    * boundary condition information is missing in the input database, or
    * when the data read in is either unknown or inconsistent.  The periodic
    * domain information is used to determine which boundary face, edge, or
    * node entries are not required from input.  Error checking
    * requires that node and edge boundary conditions are consistent
    * with those specified for the faces.
    *
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param bdry_strategy user-defined object that reads DIRICHLET or NEUMANN
    *                      conditions
    * @param input_db      input database containing all boundary data
    * @param face_conds    vector into which integer face boundary condition
    *                      types are read
    * @param edge_conds    vector into which integer edge boundary condition
    *                      types are read
    * @param node_conds    vector into which integer node boundary condition
    *                      types are read
    * @param periodic      integer vector specifying which coordinate
    *                      directions are periodic (e.g., value returned from
    *                      GridGeometry2::getPeriodicShift())
    */
   static void
   getFromInput(
      appu::BoundaryUtilityStrategy* bdry_strategy,
      const std::shared_ptr<tbox::Database>& input_db,
      std::vector<int>& face_conds,
      std::vector<int>& edge_conds,
      std::vector<int>& node_conds,
      const hier::IntVector& periodic);

   /*!
    * Function to fill 3d face boundary values for a patch.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param varname             String name of variable (for error reporting).
    * @param vardata             Cell-centered patch data object to fill.
    * @param patch               hier::Patch on which data object lives.
    * @param ghost_width_to_fill Width of ghost region to fill.
    * @param bdry_face_conds     std::vector of boundary condition types for
    *                            patch faces.
    * @param bdry_face_values    std::vector of boundary values for patch faces.
    */
   static void
   fillFaceBoundaryData(
      const std::string& varname,
      std::shared_ptr<pdat::CellData<double> >& vardata,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const std::vector<int>& bdry_face_conds,
      const std::vector<double>& bdry_face_values);

   /*!
    * Function to fill 3d edge boundary values for a patch.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param varname             String name of variable (for error reporting).
    * @param vardata             Cell-centered patch data object to fill.
    * @param patch               hier::Patch on which data object lives.
    * @param ghost_width_to_fill Width of ghost region to fill.
    * @param bdry_edge_conds     std::vector of boundary condition types for
    *                            patch edges.
    * @param bdry_face_values    std::vector of boundary values for patch faces.
    */
   static void
   fillEdgeBoundaryData(
      const std::string& varname,
      std::shared_ptr<pdat::CellData<double> >& vardata,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const std::vector<int>& bdry_edge_conds,
      const std::vector<double>& bdry_face_values);

   /*!
    * Function to fill 3d node boundary values for a patch.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @param varname             String name of variable (for error reporting).
    * @param vardata             Cell-centered patch data object to fill.
    * @param patch               hier::Patch on which data object lives.
    * @param ghost_width_to_fill Width of ghost region to fill.
    * @param bdry_node_conds     std::vector of boundary condition types for
    *                            patch nodes.
    * @param bdry_face_values    std::vector of boundary values for patch faces.
    */
   static void
   fillNodeBoundaryData(
      const std::string& varname,
      std::shared_ptr<pdat::CellData<double> >& vardata,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_fill,
      const std::vector<int>& bdry_node_conds,
      const std::vector<double>& bdry_face_values);

   /*!
    * Function that returns the integer face boundary location
    * corresponding to the given edge location and edge boundary
    * condition.
    *
    * If the edge boundary condition type or edge location are unknown,
    * or the boundary condition type is inconsistant with the edge location
    * an error results.
    *
    * @return Integer face location for edge location and boundary condition type.
    *
    * @param edge_loc   Integer location for edge.
    * @param edge_btype Integer boundary condition type for edge.
    */
   static int
   getFaceLocationForEdgeBdry(
      int edge_loc,
      int edge_btype);

   /*!
    * Function that returns the integer face boundary location
    * corresponding to the given node location and node boundary
    * condition.
    *
    * If the node boundary condition type or node location are unknown,
    * or the boundary condition type is inconsistant with the node location
    * an error results.
    *
    * @return Integer face location for node location and boundary condition type.
    *
    * @param node_loc   Integer location for node.
    * @param node_btype Integer boundary condition type for node.
    */
   static int
   getFaceLocationForNodeBdry(
      int node_loc,
      int node_btype);

   /*!
    * Function to check 3d boundary data for a patch data quantity on
    * a patch after it is set.  A warning message will be sent to log
    * file for each bad boundary value that is found.
    *
    * When assertion checking is active, assertions will result when any
    * of the pointer arguments is null, or an array is passed in with the
    * the wrong size.
    *
    * @return Integer number of bad boundary values found.
    *
    * @param varname       String name of variable (for error reporting).
    * @param patch         hier::Patch on which data object lives.
    * @param data_id       hier::Patch data index for data on patch.
    * @param depth         Depth index of patch data to check.
    * @param gcw_to_check  Width of ghost region to check.
    * @param bbox          Boundary box to check.
    * @param bcase         Boundary condition type for given edge or node.
    * @param bstate        Boundary value that applies in DIRICHLET or NEUMANN case.
    */
   static int
   checkBdryData(
      const std::string& varname,
      const hier::Patch& patch,
      int data_id,
      int depth,
      const hier::IntVector& gcw_to_check,
      const hier::BoundaryBox& bbox,
      int bcase,
      double bstate);

private:
   static bool s_fortran_constants_stuffed;

   static void
   read3dBdryFaces(
      appu::BoundaryUtilityStrategy* bdry_strategy,
      std::shared_ptr<tbox::Database> input_db,
      std::vector<int>& face_conds,
      const hier::IntVector& periodic);

   static void
   read3dBdryEdges(
      std::shared_ptr<tbox::Database> input_db,
      const std::vector<int>& face_conds,
      std::vector<int>& edge_conds,
      const hier::IntVector& periodic);

   static void
   read3dBdryNodes(
      std::shared_ptr<tbox::Database> input_db,
      const std::vector<int>& face_conds,
      std::vector<int>& node_conds,
      const hier::IntVector& periodic);

   static void
   get3dBdryDirectionCheckValues(
      tbox::Dimension::dir_t& idir,
      int& offsign,
      int btype,
      int bloc,
      int bcase);

   static void
   stuff3dBdryFortConst();
};

#endif
