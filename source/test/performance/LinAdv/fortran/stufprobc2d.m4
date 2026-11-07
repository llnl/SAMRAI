c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine in 2d.
c
      subroutine stufprobc2d(
     &  PIECEWISE_CONSTANT_Xin,PIECEWISE_CONSTANT_Yin,
     &  PIECEWISE_CONSTANT_Zin,
     &  SINE_CONSTANT_Xin,SINE_CONSTANT_Yin,SINE_CONSTANT_Zin,
     &  SPHEREin,
     &  CELLGin,FACEGin,FLUXGin)
      implicit none
      integer 
     &  PIECEWISE_CONSTANT_Xin,PIECEWISE_CONSTANT_Yin,
     &  PIECEWISE_CONSTANT_Zin,
     &  SINE_CONSTANT_Xin,SINE_CONSTANT_Yin,SINE_CONSTANT_Zin,
     &  SPHEREin,
     &  CELLGin,FACEGin,FLUXGin
include(FORTDIR/probparams.i)dnl

      PIECEWISE_CONSTANT_X=PIECEWISE_CONSTANT_Xin
      PIECEWISE_CONSTANT_Y=PIECEWISE_CONSTANT_Yin
      PIECEWISE_CONSTANT_Z=PIECEWISE_CONSTANT_Zin
      SINE_CONSTANT_X=SINE_CONSTANT_Xin
      SINE_CONSTANT_Y=SINE_CONSTANT_Yin  
      SINE_CONSTANT_Z=SINE_CONSTANT_Zin
      SPHERE=SPHEREin
      CELLG=CELLGin
      FACEG=FACEGin
      FLUXG=FLUXGin

      return
      end
