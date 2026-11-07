c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for debugging 1d patch data types.
c
include(PDAT_FORTDIR/pdat_dbugstuff1d.i)dnl
c
c***********************************************************************
c Debugging routines for 1d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub1d(
pdat_debug_cell_1d(`double precision')dnl
c
      subroutine dbugcellflot1d(
pdat_debug_cell_1d(`real')dnl
c
      subroutine dbugcellcplx1d(
pdat_debug_cell_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub1d(
pdat_debug_face_1d(`double precision')dnl
c
      subroutine dbugfaceflot1d(
pdat_debug_face_1d(`real')dnl
c
      subroutine dbugfacecplx1d(
pdat_debug_face_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub1d(
pdat_debug_node_1d(`double precision')dnl
c
      subroutine dbugnodeflot1d(
pdat_debug_node_1d(`real')dnl
c
      subroutine dbugnodecplx1d(
pdat_debug_node_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub1d(
pdat_debug_outerface_1d(`double precision')dnl
c
      subroutine dbugoutfaceflot1d(
pdat_debug_outerface_1d(`real')dnl
c
      subroutine dbugoutfacecplx1d(
pdat_debug_outerface_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub1d(
pdat_debug_outerside_1d(`double precision')dnl
c
      subroutine dbugoutsideflot1d(
pdat_debug_outerside_1d(`real')dnl
c
      subroutine dbugoutsidecplx1d(
pdat_debug_outerside_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub1d(
pdat_debug_side_1d(`double precision')dnl
c
      subroutine dbugsideflot1d(
pdat_debug_side_1d(`real')dnl
c
      subroutine dbugsidecplx1d(
pdat_debug_side_1d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 1d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub1d(
pdat_debug_edge_1d(`double precision')dnl
c
      subroutine dbugedgeflot1d(
pdat_debug_edge_1d(`real')dnl
c
      subroutine dbugedgecplx1d(
pdat_debug_edge_1d(`double complex')dnl
c
