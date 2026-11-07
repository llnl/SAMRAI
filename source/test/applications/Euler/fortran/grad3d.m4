c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to compute gradients for 3d euler equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine detectgrad3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  vghost2,tagghost2,ttagghost2,
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
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dotag,donttag,
     &  vghost0,vghost1,vghost2,
     &  tagghost0,tagghost1,tagghost2,
     &  ttagghost0,ttagghost1,ttagghost2
      REAL
     &  dx(0:NDIM-1),
     &  gradtol
c variables indexed as 3dimensional
      REAL
     &  var(CELL3dVECG(ifirst,ilast,vghost))
      integer
     &  tags(CELL3dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL3dVECG(ifirst,ilast,ttagghost))
c
      REAL tol
      REAL facejump, loctol
      REAL presm1,presp1
      REAL diag(0:NDIM-1),diag012 
      logical tagcell
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      tol = gradtol
      diag(0) = sqrt(dx(2)**2+dx(1)**2)
      diag(1) = sqrt(dx(0)**2+dx(2)**2)
      diag(2) = sqrt(dx(0)**2+dx(1)**2)
      diag012= sqrt(dx(0)**2+dx(1)**2+dx(2)**2)

      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0

            if (tags(ic0,ic1,ic2) .ne. 0) then
              loctol = 0.125*tol
            else 
              loctol = tol
            endif
     
            tagcell = .false.
  
            presm1 = var(ic0-1,ic1,ic2)
            presp1 = var(ic0+1,ic1,ic2)
            facejump = abs(var(ic0,ic1,ic2)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
            tagcell = ((facejump).gt.(loctol*dx(0)))
            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2)
              presp1 = var(ic0,ic1+1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*dx(1)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1,ic2-1)
              presp1 = var(ic0,ic1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*dx(2)))
            endif


c   2Dimensional diagonals

            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2-1)
              presp1 = var(ic0,ic1+1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(0)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1+1,ic2-1)
              presp1 = var(ic0,ic1-1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(0)))
            endif
  
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2-1)
              presp1 = var(ic0+1,ic1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(1)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2+1)
              presp1 = var(ic0+1,ic1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(1)))
            endif

            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2)
              presp1 = var(ic0+1,ic1+1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(2)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2)
              presp1 = var(ic0+1,ic1-1,ic2)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag(2)))
            endif

c   End 2Dimensional diagonals
c   3Dimensional diagonals
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2-1)
              presp1 = var(ic0+1,ic1+1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2+1)
              presp1 = var(ic0+1,ic1+1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2-1)
              presp1 = var(ic0+1,ic1-1,ic2+1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2+1)
              presp1 = var(ic0+1,ic1-1,ic2-1)
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((facejump).gt.(loctol*diag012))
            endif

c   End 3Dimensional diagonals

            if ( tagcell ) then
              temptags(ic0,ic1,ic2) = dotag
            endif
          enddo
        enddo
      enddo
      return
      end

      subroutine detectshock3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  vghost2,tagghost2,ttagghost2,
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
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dotag,donttag,
     &  vghost0,vghost1,vghost2,
     &  tagghost0,tagghost1,tagghost2,
     &  ttagghost0,ttagghost1,ttagghost2
      REAL
     &  dx(0:NDIM-1),
     &  gradtol,gradonset
c variables indexed as 3dimensional
      REAL
     &  var(CELL3dVECG(ifirst,ilast,vghost))
      integer
     &  tags(CELL3dVECG(ifirst,ilast,tagghost)),
     &  temptags(CELL3dVECG(ifirst,ilast,ttagghost))
c
      REAL tol,onset
      REAL jump1, jump2, facejump, loctol,locon
      REAL presm1,presm2,presp1,presp2
      REAL diag(0:NDIM-1),diag012 
      logical tagcell
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      tol = gradtol
      onset = gradonset
      diag(0) = sqrt(dx(2)**2+dx(1)**2)
      diag(1) = sqrt(dx(0)**2+dx(2)**2)
      diag(2) = sqrt(dx(0)**2+dx(1)**2)
      diag012= sqrt(dx(0)**2+dx(1)**2+dx(2)**2)

      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0

            if (tags(ic0,ic1,ic2) .ne. 0) then
              loctol = 0.125*tol
              locon = 0.66*onset
            else 
              loctol = tol
              locon = onset
            endif
     
            tagcell = .false.
  
            presm1 = var(ic0-1,ic1,ic2)
            presm2 = var(ic0-2,ic1,ic2)
            presp1 = var(ic0+1,ic1,ic2)
            presp2 = var(ic0+2,ic1,ic2)
            jump2 = presp2-presm2
            jump1 = presp1-presm1
            facejump = abs(var(ic0,ic1,ic2)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
            tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                ((jump1*jump2).lt.zero)).and.
     &               ((facejump).gt.(loctol*dx(0))))
            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2)
              presm2 = var(ic0,ic1-2,ic2)
              presp1 = var(ic0,ic1+1,ic2)
              presp2 = var(ic0,ic1+2,ic2)
              jump2 = presp2-presm2
              jump1 = presp1-presm1
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                    ((jump1*jump2).lt.zero)).and.
     &                   ((facejump).gt.(loctol*dx(1)))) 
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1,ic2-1)
              presm2 = var(ic0,ic1,ic2-2)
              presp1 = var(ic0,ic1,ic2+1)
              presp2 = var(ic0,ic1,ic2+2)
              jump2 = presp2-presm2
              jump1 = presp1-presm1
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*dx(2)))) 
            endif


c   2Dimensional diagonals

            if (.not.tagcell) then
              presm1 = var(ic0,ic1-1,ic2-1)
              presp1 = var(ic0,ic1+1,ic2+1)
              presm2 = var(ic0,ic1-2,ic2-2)
              presp2 = var(ic0,ic1+2,ic2+2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                    ((jump1*jump2).lt.zero)).and.
     &                   ((facejump).gt.(loctol*diag(0))))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0,ic1+1,ic2-1)
              presp1 = var(ic0,ic1-1,ic2+1)
              presm2 = var(ic0,ic1+2,ic2-2)
              presp2 = var(ic0,ic1-2,ic2+2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                    ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag(0))))
            endif
  
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2-1)
              presp1 = var(ic0+1,ic1,ic2+1)
              presm2 = var(ic0-2,ic1,ic2-2)
              presp2 = var(ic0+2,ic1,ic2+2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag(1))))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1,ic2+1)
              presp1 = var(ic0+1,ic1,ic2-1)
              presm2 = var(ic0-2,ic1,ic2+2)
              presp2 = var(ic0+2,ic1,ic2-2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                    ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag(1))))
            endif

            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2)
              presp1 = var(ic0+1,ic1+1,ic2)
              presm2 = var(ic0-2,ic1-2,ic2)
              presp2 = var(ic0+2,ic1+2,ic2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                    ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag(2))))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2)
              presp1 = var(ic0+1,ic1-1,ic2)
              presm2 = var(ic0-2,ic1+2,ic2)
              presp2 = var(ic0+2,ic1-2,ic2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag(2))))
            endif

c   End 2Dimensional diagonals
c   3Dimensional diagonals
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2-1)
              presp1 = var(ic0+1,ic1+1,ic2+1)
              presm2 = var(ic0-2,ic1-2,ic2-2)
              presp2 = var(ic0+2,ic1+2,ic2+2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                   ((facejump).gt.(loctol*diag012)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1-1,ic2+1)
              presp1 = var(ic0+1,ic1+1,ic2-1)
              presm2 = var(ic0-2,ic1-2,ic2+2)
              presp2 = var(ic0+2,ic1+2,ic2-2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag012)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2-1)
              presp1 = var(ic0+1,ic1-1,ic2+1)
              presm2 = var(ic0-2,ic1+2,ic2-2)
              presp2 = var(ic0+2,ic1-2,ic2+2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                   ((facejump).gt.(loctol*diag012)))
            endif
            if (.not.tagcell) then
              presm1 = var(ic0-1,ic1+1,ic2+1)
              presp1 = var(ic0+1,ic1-1,ic2-1)
              presm2 = var(ic0-2,ic1+2,ic2+2)
              presp2 = var(ic0+2,ic1-2,ic2-2)
              jump1 = presp1-presm1
              jump2 = presp2-presm2
              facejump = abs(var(ic0,ic1,ic2)-presm1)
              facejump = max(facejump,abs(var(ic0,ic1,ic2)-presp1))
              tagcell = ((((abs(jump2)*locon).le.abs(jump1)).or.
     &                  ((jump1*jump2).lt.zero)).and.
     &                 ((facejump).gt.(loctol*diag012)))
            endif

c   End 3Dimensional diagonals

            if ( tagcell ) then
              temptags(ic0,ic1,ic2) = dotag
            endif
          enddo
        enddo
      enddo
      return
      end

