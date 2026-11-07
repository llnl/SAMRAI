c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for scalar Poisson Hypre solver.
c
c***********************************************************************
c***********************************************************************
      subroutine compdiagvariablec3d(
     &  diag, c,
     &  offdiagi, offdiagj, offdiagk,
     &  ifirst, ilast, jfirst, jlast, kfirst, klast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      double precision diag(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision c(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision 
     &  offdiagi(ifirst:ilast+1,jfirst:jlast,kfirst:klast),
     &  offdiagj(ifirst:ilast,jfirst:jlast+1,kfirst:klast),
     &  offdiagk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)
      double precision cscale, dscale
      integer i, j, k
c     Assume g value of zero
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               diag(i,j,k) = cscale*c(i,j,k) - 
     &              ( offdiagi(i,j,k) + offdiagi(i+1,j,k) +
     &                offdiagj(i,j,k) + offdiagj(i,j+1,k) +
     &                offdiagk(i,j,k) + offdiagk(i,j,k+1) )
            enddo
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compdiagscalarc3d(
     &  diag, c,
     &  offdiagi, offdiagj, offdiagk,
     &  ifirst, ilast, jfirst, jlast, kfirst, klast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      double precision diag(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision c
      double precision 
     &  offdiagi(ifirst:ilast+1,jfirst:jlast,kfirst:klast),
     &  offdiagj(ifirst:ilast,jfirst:jlast+1,kfirst:klast),
     &  offdiagk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)
      double precision cscale, dscale
      integer i, j, k
c     Assume g value of zero
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               diag(i,j,k) = cscale*c - 
     &              ( offdiagi(i,j,k) + offdiagi(i+1,j,k) +
     &                offdiagj(i,j,k) + offdiagj(i,j+1,k) +
     &                offdiagk(i,j,k) + offdiagk(i,j,k+1) )
            enddo
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compdiagzeroc3d(
     &  diag,
     &  offdiagi, offdiagj, offdiagk,
     &  ifirst, ilast, jfirst, jlast, kfirst, klast,
     &  cscale, dscale )
c***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      double precision diag(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision 
     &  offdiagi(ifirst:ilast+1,jfirst:jlast,kfirst:klast),
     &  offdiagj(ifirst:ilast,jfirst:jlast+1,kfirst:klast),
     &  offdiagk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)
      double precision cscale, dscale
      integer i, j, k
c     Assume g value of zero
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               diag(i,j,k) = - ( offdiagi(i,j,k) + offdiagi(i+1,j,k) +
     &                           offdiagj(i,j,k) + offdiagj(i,j+1,k) +
     &                           offdiagk(i,j,k) + offdiagk(i,j,k+1) )
            enddo
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjbdry3d(
     &  diag,
     &  offdiagi, offdiagj, offdiagk,
     &  pifirst, pilast, pjfirst, pjlast, pkfirst, pklast,
     &  acoef,
     &  bcoef,
     &  aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &  lower, upper,
     &  location, h )
c***********************************************************************
      implicit none
      integer pifirst, pilast, pjfirst, pjlast, pkfirst, pklast
      double precision 
     &  diag(pifirst:pilast,pjfirst:pjlast,pkfirst:pklast)
      double precision
     &  offdiagi(pifirst:pilast+1,pjfirst:pjlast,pkfirst:pklast),
     &  offdiagj(pifirst:pilast,pjfirst:pjlast+1,pkfirst:pklast),
     &  offdiagk(pifirst:pilast,pjfirst:pjlast,pkfirst:pklast+1)
      integer aifirst, ailast, ajfirst, ajlast, akfirst, aklast
      double precision 
     &  acoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      double precision 
     &  bcoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      integer kifirst, kilast, kjfirst, kjlast, kkfirst, kklast
      double precision
     &  auk0(kifirst:kilast,kjfirst:kjlast,kkfirst:kklast)
      integer lower(0:2), upper(0:2)
      integer location
      double precision h(0:2), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      integer kgho, kint, kfac
      double precision uk0, k1
c     Nomenclature for indices: gho=ghost, int=interior,
c     fac=surface, beg=beginning, end=ending.
      integer i, j, k
      hh = h(location/2)
      if ( location .eq. 0 ) then
c     min i side
         igho = upper(0)
         ifac = igho + 1
         iint = igho + 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiagi(ifac,j,k)
               auk0(ifac,j,k) = uk0*offdiagi(ifac,j,k)
               offdiagi(ifac,j,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 1 ) then
c     max i side
         igho = lower(0)
         ifac = igho
         iint = igho - 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiagi(ifac,j,k)
               auk0(ifac,j,k) = uk0*offdiagi(ifac,j,k)
               offdiagi(ifac,j,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 2 ) then
c     min j side
         jgho = upper(1)
         jfac = jgho + 1
         jint = jgho + 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiagj(i,jfac,k)
               auk0(i,jfac,k) = uk0*offdiagj(i,jfac,k)
               offdiagj(i,jfac,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 3 ) then
c     max j side
         jgho = lower(1)
         jfac = jgho
         jint = jgho - 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiagj(i,jfac,k)
               auk0(i,jfac,k) = uk0*offdiagj(i,jfac,k)
               offdiagj(i,jfac,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 4 ) then
c     min k side
         kgho = upper(2)
         kfac = kgho + 1
         kint = kgho + 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiagk(i,j,kfac)
               auk0(i,j,kfac) = uk0*offdiagk(i,j,kfac)
               offdiagk(i,j,kfac) = 0.0
            enddo
         enddo
      elseif ( location .eq. 5 ) then
c     max k side
         kgho = lower(2)
         kfac = kgho
         kint = kgho - 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiagk(i,j,kfac)
               auk0(i,j,kfac) = uk0*offdiagk(i,j,kfac)
               offdiagk(i,j,kfac) = 0.0
            enddo
         enddo
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjbdryconstoffdiags3d(
     &  diag,
     &  offdiag,
     &  pifirst, pilast, pjfirst, pjlast, pkfirst, pklast,
     &  acoef,
     &  aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &  lower, upper,
     &  location, h )
c***********************************************************************
      implicit none
      integer pifirst, pilast, pjfirst, pjlast, pkfirst, pklast
      double precision 
     &  diag(pifirst:pilast,pjfirst:pjlast,pkfirst:pklast)
      double precision
     &  offdiag(0:2)
      integer aifirst, ailast, ajfirst, ajlast, akfirst, aklast
      double precision 
     &  acoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      integer kifirst, kilast, kjfirst, kjlast, kkfirst, kklast
      double precision
     &  auk0(kifirst:kilast,kjfirst:kjlast,kkfirst:kklast)
      integer lower(0:2), upper(0:2)
      integer location
      double precision h(0:2), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      integer kgho, kint, kfac
      double precision uk0, k1
c     Nomenclature for indices: gho=ghost, int=interior,
c     fac=surface, beg=beginning, end=ending.
      integer i, j, k
      hh = h(location/2)
      if ( location .eq. 0 ) then
c     min i side
         igho = upper(0)
         ifac = igho + 1
         iint = igho + 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiag(0)
               auk0(ifac,j,k) = uk0*offdiag(0)
            enddo
         enddo
      elseif ( location .eq. 1 ) then
c     max i side
         igho = lower(0)
         ifac = igho
         iint = igho - 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiag(0)
               auk0(ifac,j,k) = uk0*offdiag(0)
            enddo
         enddo
      elseif ( location .eq. 2 ) then
c     min j side
         jgho = upper(1)
         jfac = jgho + 1
         jint = jgho + 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiag(1)
               auk0(i,jfac,k) = uk0*offdiag(1)
            enddo
         enddo
      elseif ( location .eq. 3 ) then
c     max j side
         jgho = lower(1)
         jfac = jgho
         jint = jgho - 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiag(1)
               auk0(i,jfac,k) = uk0*offdiag(1)
            enddo
         enddo
      elseif ( location .eq. 4 ) then
c     min k side
         kgho = upper(2)
         kfac = kgho + 1
         kint = kgho + 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiag(2)
               auk0(i,j,kfac) = uk0*offdiag(2)
            enddo
         enddo
      elseif ( location .eq. 5 ) then
c     max k side
         kgho = lower(2)
         kfac = kgho
         kint = kgho - 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiag(2)
               auk0(i,j,kfac) = uk0*offdiag(2)
            enddo
         enddo
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine adjustrhs3d(
     &  rhs,
     &  rifirst, rilast, rjfirst, rjlast, rkfirst, rklast,
     &  auk0,
     &  kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &  gcoef,
     &  aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &  lower, upper,
     &  location )
c***********************************************************************
      implicit none
      integer rifirst, rilast, rjfirst, rjlast, rkfirst, rklast
      double precision rhs(rifirst:rilast,rjfirst:rjlast,rkfirst:rklast)
      integer kifirst, kilast, kjfirst, kjlast, kkfirst, kklast
      double precision
     & auk0(kifirst:kilast,kjfirst:kjlast,kkfirst:kklast)
      integer aifirst, ailast, ajfirst, ajlast, akfirst, aklast
      double precision
     & gcoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      integer lower(0:2), upper(0:2)
      integer location
      integer igho, icel, ifac
      integer jgho, jcel, jfac
      integer kgho, kcel, kfac
c     Nomenclature for indices: cel=first-cell, gho=ghost,
c     beg=beginning, end=ending.
      integer i, j, k
      if ( location .eq. 0 ) then
c     min i side
         igho = upper(0)
         ifac = igho + 1
         icel = igho + 1
         do j=lower(1),upper(1)
            do k=lower(2),upper(2)
               rhs(icel,j,k) = rhs(icel,j,k)
     &                       - auk0(ifac,j,k)*gcoef(ifac,j,k)
            enddo
         enddo
      elseif ( location .eq. 1 ) then
c     max i side
         igho = lower(0)
         ifac = igho
         icel = igho - 1
         do j=lower(1),upper(1)
            do k=lower(2),upper(2)
               rhs(icel,j,k) = rhs(icel,j,k)
     &                       - auk0(ifac,j,k)*gcoef(ifac,j,k)
            enddo
         enddo
      elseif ( location .eq. 2 ) then
c     min j side
         jgho = upper(1)
         jfac = jgho + 1
         jcel = jgho + 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               rhs(i,jcel,k) = rhs(i,jcel,k)
     &                       - auk0(i,jfac,k)*gcoef(i,jfac,k)
            enddo
         enddo
      elseif ( location .eq. 3 ) then
c     max j side
         jgho = lower(1)
         jfac = jgho
         jcel = jgho - 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               rhs(i,jcel,k) = rhs(i,jcel,k)
     &                       - auk0(i,jfac,k)*gcoef(i,jfac,k)
            enddo
         enddo
      elseif ( location .eq. 4 ) then
c     min k side
         kgho = upper(2)
         kfac = kgho + 1
         kcel = kgho + 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               rhs(i,j,kcel) = rhs(i,j,kcel)
     &                       - auk0(i,j,kfac)*gcoef(i,j,kfac)
            enddo
         enddo
      elseif ( location .eq. 5 ) then
c     max k side
         kgho = lower(2)
         kfac = kgho
         kcel = kgho - 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               rhs(i,j,kcel) = rhs(i,j,kcel)
     &                       - auk0(i,j,kfac)*gcoef(i,j,kfac)
            enddo
         enddo
      endif
      return
      end
c***********************************************************************
