c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial refining of 1d patch data
c                on a regular Cartesian mesh.
c
include(PDAT_FORTDIR/pdat_m4conrefineops1d.i)dnl
c
c***********************************************************************
c Constant interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub1d(
conref_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot1d(
conref_op_cell_1d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx1d(
conref_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 1d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg1d(
conref_op_cell_1d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub1d(
conref_op_edge_1d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot1d(
conref_op_edge_1d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 1d edge-centered complex data
c***********************************************************************
c
      subroutine conrefedgecplx1d(
conref_op_edge_1d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 1d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg1d(
conref_op_edge_1d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub1d(
conref_op_face_1d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot1d(
conref_op_face_1d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 1d face-centered complex data
c***********************************************************************
c
      subroutine conreffacecplx1d(
conref_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 1d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg1d(
conref_op_face_1d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub1d(
conref_op_outerface_1d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot1d(
conref_op_outerface_1d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 1d outerface complex data
c***********************************************************************
c
      subroutine conrefoutfacecplx1d(
conref_op_outerface_1d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 1d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg1d(
conref_op_outerface_1d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub1d(
conref_op_side_1d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot1d(
conref_op_side_1d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 1d side-centered complex data
c***********************************************************************
c
      subroutine conrefsidecplx1d(
conref_op_side_1d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 1d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg1d(
conref_op_side_1d(`integer')dnl
c
