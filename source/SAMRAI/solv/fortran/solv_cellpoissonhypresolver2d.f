c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for scalar Poisson Hypre solver.
c
c***********************************************************************
c***********************************************************************
      subroutine compdiagvariablec2d(
     &  diag, c,
     &  offdiagi, offdiagj,
     &  ifirst, ilast, jfirst, jlast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast
      double precision diag(ifirst:ilast,jfirst:jlast)
      double precision c(ifirst:ilast,jfirst:jlast)
      double precision offdiagi(ifirst:ilast+1,jfirst:jlast)
      double precision offdiagj(ifirst:ilast,jfirst:jlast+1)
      double precision cscale, dscale
      integer i, j
c     Assume g value of zero
      do i=ifirst,ilast
         do j=jfirst,jlast
            diag(i,j) = cscale*c(i,j) - 
     &        ( offdiagi(i,j) + offdiagi(i+1,j) +
     &          offdiagj(i,j) + offdiagj(i,j+1) )
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compdiagscalarc2d(
     &  diag, c,
     &  offdiagi, offdiagj,
     &  ifirst, ilast, jfirst, jlast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast
      double precision diag(ifirst:ilast,jfirst:jlast)
      double precision c
      double precision offdiagi(ifirst:ilast+1,jfirst:jlast)
      double precision offdiagj(ifirst:ilast,jfirst:jlast+1)
      double precision cscale, dscale
      integer i, j
c     Assume g value of zero
      do i=ifirst,ilast
         do j=jfirst,jlast
            diag(i,j) = cscale*c - 
     &        ( offdiagi(i,j) + offdiagi(i+1,j) +
     &          offdiagj(i,j) + offdiagj(i,j+1) )
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compdiagzeroc2d(
     &  diag,
     &  offdiagi, offdiagj,
     &  ifirst, ilast, jfirst, jlast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast
      double precision diag(ifirst:ilast,jfirst:jlast)
      double precision offdiagi(ifirst:ilast+1,jfirst:jlast)
      double precision offdiagj(ifirst:ilast,jfirst:jlast+1)
      double precision cscale, dscale
      integer i, j
c     Assume g value of zero
      do i=ifirst,ilast
         do j=jfirst,jlast
            diag(i,j) = -( offdiagi(i,j) + offdiagi(i+1,j) +
     &                     offdiagj(i,j) + offdiagj(i,j+1) )
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjbdry2d(
     &  diag,
     &  offdiagi, offdiagj,
     &  pifirst, pilast, pjfirst, pjlast,
     &  acoef,
     &  bcoef,
     &  aifirst, ailast, ajfirst, ajlast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast,
     &  lower, upper,
     &  location, h )
c***********************************************************************
      implicit none
      integer pifirst, pilast, pjfirst, pjlast
      double precision diag(pifirst:pilast,pjfirst:pjlast)
      double precision offdiagi(pifirst:pilast+1,pjfirst:pjlast)
      double precision offdiagj(pifirst:pilast,pjfirst:pjlast+1)
      integer aifirst, ailast, ajfirst, ajlast
      double precision acoef(aifirst:ailast,ajfirst:ajlast)
      double precision bcoef(aifirst:ailast,ajfirst:ajlast)
      integer kifirst, kilast, kjfirst, kjlast
      double precision auk0(kifirst:kilast,kjfirst:kjlast)
      integer lower(0:1), upper(0:1)
      integer location
      double precision h(0:1), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      double precision uk0, k1
c     Nomenclature for indices: gho=ghost, int=interior,
c     fac=surface, beg=beginning, end=ending.
      integer i, j
      hh = h(location/2)
      if ( location .eq. 0 ) then
c        min i edge
         igho = upper(0)
         iint = igho + 1
         ifac = igho + 1
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j) + k1*offdiagi(ifac,j)
            auk0(ifac,j) = uk0*offdiagi(ifac,j)
            offdiagi(ifac,j) = 0.0
         enddo
      elseif ( location .eq. 1 ) then
c        min i edge
         igho = lower(0)
         iint = igho - 1
         ifac = igho
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j) + k1*offdiagi(ifac,j)
            auk0(ifac,j) = uk0*offdiagi(ifac,j)
            offdiagi(ifac,j) = 0.0
         enddo
      elseif ( location .eq. 2 ) then
c        min i edge
         jgho = upper(1)
         jint = jgho + 1
         jfac = jgho + 1
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint) + k1*offdiagj(i,jfac)
            auk0(i,jfac) = uk0*offdiagj(i,jfac)
            offdiagj(i,jfac) = 0.0
         enddo
      elseif ( location .eq. 3 ) then
c        min i edge
         jgho = lower(1)
         jint = jgho - 1
         jfac = jgho
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint) + k1*offdiagj(i,jfac)
            auk0(i,jfac) = uk0*offdiagj(i,jfac)
            offdiagj(i,jfac) = 0.0
         enddo
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjbdryconstoffdiags2d(
     &  diag,
     &  offdiag,
     &  pifirst, pilast, pjfirst, pjlast,
     &  acoef,
     &  aifirst, ailast, ajfirst, ajlast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast,
     &  lower, upper,
     &  location, h )
c***********************************************************************
      implicit none
      integer pifirst, pilast, pjfirst, pjlast
      double precision diag(pifirst:pilast,pjfirst:pjlast)
      double precision offdiag(0:1)
      integer aifirst, ailast, ajfirst, ajlast
      double precision acoef(aifirst:ailast,ajfirst:ajlast)
      integer kifirst, kilast, kjfirst, kjlast
      double precision auk0(kifirst:kilast,kjfirst:kjlast)
      integer lower(0:1), upper(0:1)
      integer location
      double precision h(0:1), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      double precision uk0, k1
c     Nomenclature for indices: gho=ghost, int=interior,
c     fac=surface, beg=beginning, end=ending.
      integer i, j
      hh = h(location/2)
      if ( location .eq. 0 ) then
c        min i edge
         igho = upper(0)
         iint = igho + 1
         ifac = igho + 1
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j) + k1*offdiag(0)
            auk0(ifac,j) = uk0*offdiag(0)
         enddo
      elseif ( location .eq. 1 ) then
c        min i edge
         igho = lower(0)
         iint = igho - 1
         ifac = igho
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j) + k1*offdiag(0)
            auk0(ifac,j) = uk0*offdiag(0)
         enddo
      elseif ( location .eq. 2 ) then
c        min i edge
         jgho = upper(1)
         jint = jgho + 1
         jfac = jgho + 1
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint) + k1*offdiag(1)
            auk0(i,jfac) = uk0*offdiag(1)
         enddo
      elseif ( location .eq. 3 ) then
c        min i edge
         jgho = lower(1)
         jint = jgho - 1
         jfac = jgho
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint) + k1*offdiag(1)
            auk0(i,jfac) = uk0*offdiag(1)
         enddo
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjustrhs2d(
     &  rhs,
     &  rifirst, rilast, rjfirst, rjlast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast,
     &  gcoef,
     &  aifirst, ailast, ajfirst, ajlast,
     &  lower, upper,
     &  location )
c***********************************************************************
      implicit none
      integer rifirst, rilast, rjfirst, rjlast
      double precision rhs(rifirst:rilast,rjfirst:rjlast)
      integer kifirst, kilast, kjfirst, kjlast
      double precision auk0(kifirst:kilast,kjfirst:kjlast)
      integer aifirst, ailast, ajfirst, ajlast
      double precision gcoef(aifirst:ailast,ajfirst:ajlast)
      integer lower(0:1), upper(0:1)
      integer location
      integer igho, iint, ifac
      integer jgho, jint, jfac
c     Nomenclature for indices: cel=first-cell, gho=ghost,
c     beg=beginning, end=ending.
      integer i, j
      if ( location .eq. 0 ) then
c        min i edge
         igho = upper(0)
         ifac = igho + 1
         iint = igho + 1
         do j=lower(1),upper(1)
            rhs(iint,j) = rhs(iint,j) - auk0(ifac,j)*gcoef(ifac,j)
         enddo
      elseif ( location .eq. 1 ) then
c        max i edge
         igho = lower(0)
         ifac = igho
         iint = igho - 1
         do j=lower(1),upper(1)
            rhs(iint,j) = rhs(iint,j) - auk0(ifac,j)*gcoef(ifac,j)
         enddo
      elseif ( location .eq. 2 ) then
c        min j edge
         jgho = upper(1)
         jfac = jgho + 1
         jint = jgho + 1
         do i=lower(0),upper(0)
            rhs(i,jint) = rhs(i,jint) - auk0(i,jfac)*gcoef(i,jfac)
         enddo
      elseif ( location .eq. 3 ) then
c        max j edge
         jgho = lower(1)
         jfac = jgho
         jint = jgho - 1
         do i=lower(0),upper(0)
            rhs(i,jint) = rhs(i,jint) - auk0(i,jfac)*gcoef(i,jfac)
         enddo
      endif
      return
      end
c***********************************************************************
