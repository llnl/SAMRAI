c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for gradient computation in 2d.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine detectgrad2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  dx,
     &  gradtol,
     &  dotag,donttag,
     &  var,
     &  tags,temptags)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  dotag,donttag,
     &  vghost0,vghost1,
     &  tagghost0,tagghost1,
     &  ttagghost0,ttagghost1
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
c variables indexed as 2dimensional
      REAL
     &  var(CELL2dVECG(ifirst,ilast,vghost))
      integer
     &  tags(CELL2dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL2dVECG(ifirst,ilast,ttagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1,diag01
      logical tagcell
      integer ic0,ic1
c
c***********************************************************************
c
      tol = gradtol
      diag01 = sqrt(dx(0)**2+dx(1)**2)

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0

          if (tags(ic0,ic1) .ne. 0) then
            loctol = 0.125*tol
          else 
            loctol = tol
          endif
   
          tagcell = .false.

          presm1 = var(ic0-1,ic1)
          presp1 = var(ic0+1,ic1)
          facejump = abs(var(ic0,ic1)-presm1)
          facejump = max(facejump,abs(var(ic0,ic1)-presp1))
          tagcell = ((facejump).gt.(loctol*dx(0)))
          if (.not.tagcell) then
            presm1 = var(ic0,ic1-1)
            presp1 = var(ic0,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*dx(1)))
          endif

          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1-1)
            presp1 = var(ic0+1,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif
          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1+1)
            presp1 = var(ic0+1,ic1-1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif

          if ( tagcell ) then
            temptags(ic0,ic1) = dotag
          endif
        enddo
      enddo
      return
      end

      subroutine detectshock2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  dx,
     &  gradtol,gradonset,
     &  dotag,donttag,
     &  var,
     &  tags,temptags)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  dotag,donttag,
     &  vghost0,vghost1,
     &  tagghost0,tagghost1,
     &  ttagghost0,ttagghost1
      REAL
     &  dx(0:NDIM-1),
     &  gradtol,gradonset
c variables indexed as 2dimensional
      REAL
     &  var(CELL2dVECG(ifirst,ilast,vghost))
      integer
     &  tags(CELL2dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL2dVECG(ifirst,ilast,ttagghost))
c
      REAL tol,onset
      REAL jump1, jump2, facejump, loctol,locon
      REAL presm1,presm2,presp1,presp2
      REAL diag01 
      logical tagcell
      integer ic0,ic1
c
c***********************************************************************
c
      tol = gradtol
      onset = gradonset
      diag01 = sqrt(dx(0)**2+dx(1)**2)

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0

          if (tags(ic0,ic1) .ne. 0) then
            loctol = 0.125*tol
            locon = 0.66*onset
          else 
            loctol = tol
            locon = onset
          endif
   
          tagcell = .false.

          presm1 = var(ic0-1,ic1)
          presm2 = var(ic0-2,ic1)
          presp1 = var(ic0+1,ic1)
          presp2 = var(ic0+2,ic1)
          jump2 = presp2-presm2
          jump1 = presp1-presm1
          facejump = abs(var(ic0,ic1)-presm1)
          facejump = max(facejump,abs(var(ic0,ic1)-presp1))
          tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                ((jump1*jump2).lt.zero)).and.
     &               ((facejump).gt.(loctol*dx(0))))
          if (.not.tagcell) then
            presm1 = var(ic0,ic1-1)
            presm2 = var(ic0,ic1-2)
            presp1 = var(ic0,ic1+1)
            presp2 = var(ic0,ic1+2)
            jump2 = presp2-presm2
            jump1 = presp1-presm1
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                ((jump1*jump2).lt.zero)).and.
     &               ((facejump).gt.(loctol*dx(1)))) 
          endif

          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1-1)
            presp1 = var(ic0+1,ic1+1)
            presm2 = var(ic0-2,ic1-2)
            presp2 = var(ic0+2,ic1+2)
            jump1 = presp1-presm1
            jump2 = presp2-presm2
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag01)))
          endif
          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1+1)
            presp1 = var(ic0+1,ic1-1)
            presm2 = var(ic0-2,ic1+2)
            presp2 = var(ic0+2,ic1-2)
            jump1 = presp1-presm1
            jump2 = presp2-presm2
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag01)))
          endif

          if ( tagcell ) then
            temptags(ic0,ic1) = dotag
          endif
        enddo
      enddo
      return
      end


