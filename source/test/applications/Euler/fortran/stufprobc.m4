c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine.
c
      subroutine stufprobc(
     &  APPROX_RIEM_SOLVEin,EXACT_RIEM_SOLVEin,HLLC_RIEM_SOLVEin,
     &  PIECEWISE_CONSTANT_Xin,PIECEWISE_CONSTANT_Yin,
     &  PIECEWISE_CONSTANT_Zin,
     &  SPHEREin,STEPin,
     &  CELLGin,FACEGin,FLUXGin)
      implicit none
      integer 
     &  APPROX_RIEM_SOLVEin,EXACT_RIEM_SOLVEin,HLLC_RIEM_SOLVEin,
     &  PIECEWISE_CONSTANT_Xin,PIECEWISE_CONSTANT_Yin,
     &  PIECEWISE_CONSTANT_Zin,
     &  SPHEREin,STEPin,
     &  CELLGin,FACEGin,FLUXGin
include(FORTDIR/probparams.i)dnl

      APPROX_RIEM_SOLVE=APPROX_RIEM_SOLVEin
      EXACT_RIEM_SOLVE=EXACT_RIEM_SOLVEin
      HLLC_RIEM_SOLVE=HLLC_RIEM_SOLVEin
      PIECEWISE_CONSTANT_X=PIECEWISE_CONSTANT_Xin
      PIECEWISE_CONSTANT_Y=PIECEWISE_CONSTANT_Yin
      PIECEWISE_CONSTANT_Z=PIECEWISE_CONSTANT_Zin
      SPHERE=SPHEREin
      STEP=STEPin
      CELLG=CELLGin
      FACEG=FACEGin
      FLUXG=FLUXGin

      return
      end
