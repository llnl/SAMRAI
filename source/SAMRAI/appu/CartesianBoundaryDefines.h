/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Defines for boundary condition integer constants
 *
 ************************************************************************/

#ifndef included_appu_CartesianBoundaryDefines
#define included_appu_CartesianBoundaryDefines

/*
 * Definitions for boundary types in 1d, 2d, and 3d:
 */

//@{
//! @name Definitions for boundary types in 1d, 2d, and 3d:
namespace Bdry {
enum Type {
   UNDEFINED = -1,

   FACE3D = 1,
   EDGE3D = 2,
   NODE3D = 3,

   EDGE2D = 1,
   NODE2D = 2,

   NODE1D = 1
};
}
//@}

/*
 * Definitions for boundary array sizes in 1d, 2d, or 3d:
 */

//@{
//! @name Definitions for boundary array sizes in 1d, 2d, or 3d:
const int NUM_1D_NODES = 2;

const int NUM_2D_EDGES = 4;
const int NUM_2D_NODES = 4;

const int NUM_3D_FACES = 6;
const int NUM_3D_EDGES = 12;
const int NUM_3D_NODES = 8;
//@}

/*
 * Definitions for Face, Edge, and Node boundary locations:
 *
 * Note that these definitions are used only for:
 * - Node boundary locations in 1d (XLO, XHI only), or
 * - Edge boundary locations in 2d (XLO, XHI, YLO, YHI only), or
 * - Face boundary locations in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary locations (see source code for more information):
namespace BdryLoc {
enum Type {
   XLO = 0,
   XHI = 1,
   YLO = 2,
   YHI = 3,
   ZLO = 4,
   ZHI = 5
};
}
//@}

/*
 * Definitions for Node boundary locations in 2d:
 */

//@{
//! @name Definitions for Node boundary locations in 2d:
namespace NodeBdyLoc2D {
enum Type {
   XLO_YLO = 0,
   XHI_YLO = 1,
   XLO_YHI = 2,
   XHI_YHI = 3
};
}
//@}

/*
 * Definitions for Edge boundary locations in 3d:
 */

//@{
//! @name Definitions for Edge boundary locations in 3d:
namespace EdgeBdyLoc3D {
enum Type {
   XLO_YLO = 0,
   XHI_YLO = 1,
   XLO_YHI = 2,
   XHI_YHI = 3,
   XLO_ZLO = 4,
   XHI_ZLO = 5,
   XLO_ZHI = 6,
   XHI_ZHI = 7,
   YLO_ZLO = 8,
   YHI_ZLO = 9,
   YLO_ZHI = 10,
   YHI_ZHI = 11
};
}
//@}

/*
 * Definitions for Node boundary locations in 3d:
 */

//@{
//! @name Definitions for Node boundary locations in 3d:
namespace NodeBdyLoc3D {
enum Type {
   XLO_YLO_ZLO = 0,
   XHI_YLO_ZLO = 1,
   XLO_YHI_ZLO = 2,
   XHI_YHI_ZLO = 3,
   XLO_YLO_ZHI = 4,
   XHI_YLO_ZHI = 5,
   XLO_YHI_ZHI = 6,
   XHI_YHI_ZHI = 7
};
}
//@}

/*
 * Definitions for Face, Edge, and Node boundary conditions:
 *
 * Note that FLOW, REFLECT, DIRICHLET and NEUMANN are used only for:
 * - Node boundary conditions in 1d, or
 * - Edge boundary conditions in 2d, or
 * - Face boundary conditions in 3d.
 *
 * Note that [X, Y, Z]FLOW, [X, Y, Z]REFLECT, [X, Y, Z]DIRICHLET, and
 * [X, Y, Z]NEUMANN are used only for:
 * - Node boundary conditions in 2d (X and Y cases only), or
 * - Edge and Node boundary conditions in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary conditions (see source code for more information):
namespace BdryCond {
enum Type {
   FLOW = 90,
   REFLECT = 91,
   DIRICHLET = 92,
   NEUMANN = 93,
   XFLOW = 900,
   YFLOW = 901,
   ZFLOW = 902,
   XREFLECT = 910,
   YREFLECT = 911,
   ZREFLECT = 912,
   XDIRICHLET = 920,
   YDIRICHLET = 921,
   ZDIRICHLET = 922,
   XNEUMANN = 930,
   YNEUMANN = 931,
   ZNEUMANN = 932
};
}
//@}

#endif
