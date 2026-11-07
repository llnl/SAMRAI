c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for debugging 3d patch data types.
c
include(PDAT_FORTDIR/pdat_dbugstuff3d.i)dnl
c
c***********************************************************************
c Debugging routines for 3d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub3d(
pdat_debug_cell_3d(`double precision')dnl
c
      subroutine dbugcellflot3d(
pdat_debug_cell_3d(`real')dnl
c
      subroutine dbugcellcplx3d(
pdat_debug_cell_3d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 3d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub3d0(
pdat_debug_face_3d(`double precision',0,1,2)dnl
c
      subroutine dbugfacedoub3d1(
pdat_debug_face_3d(`double precision',1,2,0)dnl
c
      subroutine dbugfacedoub3d2(
pdat_debug_face_3d(`double precision',2,0,1)dnl
c
      subroutine dbugfaceflot3d0(
pdat_debug_face_3d(`real',0,1,2)dnl
c
      subroutine dbugfaceflot3d1(
pdat_debug_face_3d(`real',1,2,0)dnl
c
      subroutine dbugfaceflot3d2(
pdat_debug_face_3d(`real',2,0,1)dnl
c
      subroutine dbugfacecplx3d0(
pdat_debug_face_3d(`double complex',0,1,2)dnl
c
      subroutine dbugfacecplx3d1(
pdat_debug_face_3d(`double complex',1,2,0)dnl
c
      subroutine dbugfacecplx3d2(
pdat_debug_face_3d(`double complex',2,0,1)dnl
c
c
c***********************************************************************
c Debugging routines for 3d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub3d(
pdat_debug_node_3d(`double precision')dnl
c
      subroutine dbugnodeflot3d(
pdat_debug_node_3d(`real')dnl
c
      subroutine dbugnodecplx3d(
pdat_debug_node_3d(`double complex')dnl
c
c
c***********************************************************************
c Debugging routines for 3d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub3d0(
pdat_debug_outerface_3d(`double precision',0,1,2)dnl
c
      subroutine dbugoutfacedoub3d1(
pdat_debug_outerface_3d(`double precision',1,2,0)dnl
c
      subroutine dbugoutfacedoub3d2(
pdat_debug_outerface_3d(`double precision',2,0,1)dnl
c
      subroutine dbugoutfaceflot3d0(
pdat_debug_outerface_3d(`real',0,1,2)dnl
c
      subroutine dbugoutfaceflot3d1(
pdat_debug_outerface_3d(`real',1,2,0)dnl
c
      subroutine dbugoutfaceflot3d2(
pdat_debug_outerface_3d(`real',2,0,1)dnl
c
      subroutine dbugoutfacecplx3d0(
pdat_debug_outerface_3d(`double complex',0,1,2)dnl
c
      subroutine dbugoutfacecplx3d1(
pdat_debug_outerface_3d(`double complex',1,2,0)dnl
c
      subroutine dbugoutfacecplx3d2(
pdat_debug_outerface_3d(`double complex',2,0,1)dnl
c
c
c***********************************************************************
c Debugging routines for 3d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub3d0(
pdat_debug_outerside_3d(`double precision',0,1,2)dnl
c
      subroutine dbugoutsidedoub3d1(
pdat_debug_outerside_3d(`double precision',1,2,0)dnl
c
      subroutine dbugoutsidedoub3d2(
pdat_debug_outerside_3d(`double precision',2,0,1)dnl
c
      subroutine dbugoutsideflot3d0(
pdat_debug_outerside_3d(`real',0,1,2)dnl
c
      subroutine dbugoutsideflot3d1(
pdat_debug_outerside_3d(`real',1,2,0)dnl
c
      subroutine dbugoutsideflot3d2(
pdat_debug_outerside_3d(`real',2,0,1)dnl
c
      subroutine dbugoutsidecplx3d0(
pdat_debug_outerside_3d(`double complex',0,1,2)dnl
c
      subroutine dbugoutsidecplx3d1(
pdat_debug_outerside_3d(`double complex',1,2,0)dnl
c
      subroutine dbugoutsidecplx3d2(
pdat_debug_outerside_3d(`double complex',2,0,1)dnl
c
c
c***********************************************************************
c Debugging routines for 3d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub3d0(
pdat_debug_side_3d(`double precision',0,1,2)dnl
c
      subroutine dbugsidedoub3d1(
pdat_debug_side_3d(`double precision',1,2,0)dnl
c
      subroutine dbugsidedoub3d2(
pdat_debug_side_3d(`double precision',2,0,1)dnl
c
      subroutine dbugsideflot3d0(
pdat_debug_side_3d(`real',0,1,2)dnl
c
      subroutine dbugsideflot3d1(
pdat_debug_side_3d(`real',1,2,0)dnl
c
      subroutine dbugsideflot3d2(
pdat_debug_side_3d(`real',2,0,1)dnl
c
      subroutine dbugsidecplx3d0(
pdat_debug_side_3d(`double complex',0,1,2)dnl
c
      subroutine dbugsidecplx3d1(
pdat_debug_side_3d(`double complex',1,2,0)dnl
c
      subroutine dbugsidecplx3d2(
pdat_debug_side_3d(`double complex',2,0,1)dnl
c
c
c***********************************************************************
c Debugging routines for 3d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub3d0(
pdat_debug_edge_3d(`double precision',0,1,2)dnl
c
      subroutine dbugedgedoub3d1(
pdat_debug_edge_3d(`double precision',1,2,0)dnl
c
      subroutine dbugedgedoub3d2(
pdat_debug_edge_3d(`double precision',2,0,1)dnl
c
      subroutine dbugedgeflot3d0(
pdat_debug_edge_3d(`real',0,1,2)dnl
c
      subroutine dbugedgeflot3d1(
pdat_debug_edge_3d(`real',1,2,0)dnl
c
      subroutine dbugedgeflot3d2(
pdat_debug_edge_3d(`real',2,0,1)dnl
c
      subroutine dbugedgecplx3d0(
pdat_debug_edge_3d(`double complex',0,1,2)dnl
c
      subroutine dbugedgecplx3d1(
pdat_debug_edge_3d(`double complex',1,2,0)dnl
c
      subroutine dbugedgecplx3d2(
pdat_debug_edge_3d(`double complex',2,0,1)dnl
c
