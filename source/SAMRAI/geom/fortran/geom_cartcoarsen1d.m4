c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial coarsening of 1d patch data
c                on a regular Cartesian mesh.
c
include(FORTDIR/geom_m4cartcoarsenops1d.i)dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub1d(
cart_wgtavg_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot1d(
cart_wgtavg_op_cell_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx1d(
cart_wgtavg_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub1d(
cart_wgtavg_op_edge_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot1d(
cart_wgtavg_op_edge_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx1d(
cart_wgtavg_op_edge_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub1d(
cart_wgtavg_op_face_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot1d(
cart_wgtavg_op_face_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx1d(
cart_wgtavg_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub1d(
cart_wgtavg_op_outerface_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot1d(
cart_wgtavg_op_outerface_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx1d(
cart_wgtavg_op_outerface_1d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 1d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub1d(
cart_wgtavg_op_outerside_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub1d(
cart_wgtavg_op_side_1d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot1d(
cart_wgtavg_op_side_1d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 1d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx1d(
cart_wgtavg_op_side_1d(`double complex')dnl
c
