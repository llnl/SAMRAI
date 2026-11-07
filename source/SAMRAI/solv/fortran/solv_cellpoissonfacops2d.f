c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for scalar Poisson FAC operator
c
c***********************************************************************
c***********************************************************************
      subroutine compfluxvardc2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &xdc , ydc , dcgi, dcgj ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer dcgi, dcgj, fluxgi, fluxgj,
     &        solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision dx(0:1)

      double precision dxi, dyi
      integer i, j

      dxi = 1./dx(0)
      dyi = 1./dx(1)

      do j=jfirst,jlast
      do i=ifirst,ilast+1
         xflux(i,j) = dxi*xdc(i,j)*( soln(i,j) - soln(i-1,j) )
      enddo
      enddo
      do j=jfirst,jlast+1
      do i=ifirst,ilast
         yflux(i,j) = dyi*ydc(i,j)*( soln(i,j) - soln(i,j-1) )
      enddo
      enddo

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compfluxcondc2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &dc ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj,
     &        solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision dx(0:1)

      double precision dxi, dyi
      integer i, j

      dxi = 1./dx(0)
      dyi = 1./dx(1)

      do j=jfirst,jlast
      do i=ifirst,ilast+1
         xflux(i,j) = dxi*dc*( soln(i,j) - soln(i-1,j) )
      enddo
      enddo
      do j=jfirst,jlast+1
      do i=ifirst,ilast
         yflux(i,j) = dyi*dc*( soln(i,j) - soln(i,j-1) )
      enddo
      enddo

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxvardcvarsf2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &xdc , ydc , dcgi, dcgj ,
     &rhs , rhsgi, rhsgj ,
     &scalar_field , sfgi, sfgj ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer dcgi, dcgj, fluxgi, fluxgj, rhsgi, rhsgj,
     &        solngi, solngj, sfgi, sfgj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj)
      double precision dx(0:1)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi
      double precision dudr
      double precision rcoef
      integer i, j
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      rcoef = 1.0

      maxres = 0.0

      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (jfirst+j)-((jfirst+j)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
            residual
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field(i,j)*soln(i,j)
            dudr = 1./( ( dxi*dxi*( xdc(i+1,j) + xdc(i,j) )
     &                  + dyi*dyi*( ydc(i,j+1) + ydc(i,j) ) )
     &                  - scalar_field(i,j) )
            du = -residual*dudr
            soln(i,j) = soln(i,j) + du*rcoef
            if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxcondcvarsf2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &dc ,
     &rhs , rhsgi, rhsgj ,
     &scalar_field , sfgi, sfgj ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj, rhsgi, rhsgj,
     &        solngi, solngj, sfgi, sfgj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj)
      double precision dx(0:1)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi
      double precision dudr
      double precision rcoef
      integer i, j
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      rcoef = 1.0

      maxres = 0.0

      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (jfirst+j)-((jfirst+j)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
            residual
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field(i,j)*soln(i,j)
            dudr = 1./( ( dxi*dxi*( dc + dc )
     &                  + dyi*dyi*( dc + dc ) )
     &                  - scalar_field(i,j) )
            du = -residual*dudr
            soln(i,j) = soln(i,j) + du*rcoef
            if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxvardcconsf2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &xdc , ydc , dcgi, dcgj ,
     &rhs , rhsgi, rhsgj ,
     &scalar_field ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer dcgi, dcgj, fluxgi, fluxgj, rhsgi, rhsgj,
     &        solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field
      double precision dx(0:1)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi
      double precision dudr
      double precision rcoef
      integer i, j
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      rcoef = 1.0

      maxres = 0.0

      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (jfirst+j)-((jfirst+j)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
            residual
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field*soln(i,j)
            dudr = 1./( ( dxi*dxi*( xdc(i+1,j) + xdc(i,j) )
     &                  + dyi*dyi*( ydc(i,j+1) + ydc(i,j) ) )
     &                  - scalar_field )
            du = -residual*dudr
            soln(i,j) = soln(i,j) + du*rcoef
            if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxcondcconsf2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &dc ,
     &rhs , rhsgi, rhsgj ,
     &scalar_field ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj, rhsgi, rhsgj,
     &        solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field
      double precision dx(0:1)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi
      double precision dudr
      double precision rcoef
      integer i, j
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      rcoef = 1.0

      maxres = 0.0

      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (jfirst+j)-((jfirst+j)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
            residual
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field*soln(i,j)
            dudr = 1./( ( dxi*dxi*( dc + dc )
     &                  + dyi*dyi*( dc + dc ) )
     &                  - scalar_field )
            du = -residual*dudr
            soln(i,j) = soln(i,j) + du*rcoef
            if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      return
      end
c***********************************************************************

c***********************************************************************
      subroutine compresvarsca2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &rhs , rhsgi, rhsgj ,
     &residual , residualgi, residualgj ,
     &scalar_field , sfgi, sfgj ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj, rhsgi, rhsgj,
     &        residualgi, residualgj, solngi, solngj, sfgi, sfgj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision residual(ifirst-rhsgi:ilast+rhsgi,
     &                          jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj)
      double precision dx(0:1)

      double precision dxi, dyi
      integer i, j

      dxi = 1./dx(0)
      dyi = 1./dx(1)

      do j=jfirst,jlast
         do i=ifirst,ilast
            residual(i,j)
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field(i,j)*soln(i,j)
         enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compresconsca2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &rhs , rhsgi, rhsgj ,
     &residual , residualgi, residualgj ,
     &scalar_field ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj, rhsgi, rhsgj,
     &        residualgi, residualgj, solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj)
      double precision residual(ifirst-rhsgi:ilast+rhsgi,
     &                          jfirst-rhsgj:jlast+rhsgj)
      double precision scalar_field
      double precision dx(0:1)

      double precision dxi, dyi
      integer i, j

      dxi = 1./dx(0)
      dyi = 1./dx(1)

      do j=jfirst,jlast
         do i=ifirst,ilast
            residual(i,j)
     &         = rhs(i,j)
     &         - ( dxi*( xflux(i+1,j) - xflux(i,j) )
     &           + dyi*( yflux(i,j+1) - yflux(i,j) ) )
     &         - scalar_field*soln(i,j)
         enddo
      enddo
      return
      end
c***********************************************************************

c***********************************************************************
      subroutine ewingfixfluxvardc2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &xdc , ydc , dcgi, dcgj ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &location_index ,
     &ratio_to_coarser ,
     &blower, bupper,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer dcgi, dcgj, fluxgi, fluxgj,
     &        solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision dx(0:1)
      integer location_index
      integer ratio_to_coarser(0:1)
c     Lower and upper corners of boundary box
      integer blower(0:1), bupper(0:1)

      double precision h
      integer i, ibeg, iend, igho, j, jbeg, jend, jgho
c     Fine grid indices inside one coarse grid.
      integer ip, jp
c     Fine grid indices for point diametrically opposite from (ip,jp).
      integer iq, jq
c     Weights associated with longtitudinal and transverse
c     (with respect to boundary normal) gradients.
      double precision tranwt, longwt

      if ( location_index .eq. 0 ) then
c        min i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = bupper(0)
         i = igho+1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do j=jbeg,jend,ratio_to_coarser(1)
            do jp=0,ratio_to_coarser(1)-1
               jq = ratio_to_coarser(1) - jp - 1
               xflux(i,j+jp)
     &            = longwt*xflux(i,j+jp)
     &            + tranwt*xdc(i,j+jp)*(
     &               soln(i,j+jq) - soln(i,j+jp) )/h
            enddo
         enddo
      elseif ( location_index .eq. 1 ) then
c        max i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = blower(0)
         i = igho-1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do j=jbeg,jend,ratio_to_coarser(1)
            do jp=0,ratio_to_coarser(1)-1
               jq = ratio_to_coarser(1) - jp - 1
               xflux(igho,j+jp)
     &            = longwt*xflux(igho,j+jp)
     &            - tranwt*xdc(igho,j+jp)*(
     &               soln(i,j+jq) - soln(i,j+jp) )/h
            enddo
         enddo
      elseif ( location_index .eq. 2 ) then
c        min j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = bupper(1)
         j = jgho+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do ip=0,ratio_to_coarser(0)-1
               iq = ratio_to_coarser(0) - ip - 1
               yflux(i+ip,j)
     &            = longwt*yflux(i+ip,j)
     &            + tranwt*ydc(i+ip,j)*(
     &               soln(i+iq,j) - soln(i+ip,j) )/h
            enddo
         enddo
      elseif ( location_index .eq. 3 ) then
c        max j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = blower(1)
         j = jgho-1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do ip=0,ratio_to_coarser(0)-1
               iq = ratio_to_coarser(0) - ip - 1
               yflux(i+ip,jgho)
     &            = longwt*yflux(i+ip,jgho)
     &            - tranwt*ydc(i+ip,jgho)*(
     &               soln(i+iq,j) - soln(i+ip,j) )/h
            enddo
         enddo
      endif

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine ewingfixfluxcondc2d(
     &xflux , yflux , fluxgi, fluxgj ,
     &dc ,
     &soln , solngi, solngj ,
     &ifirst, ilast, jfirst, jlast ,
     &location_index ,
     &ratio_to_coarser ,
     &blower, bupper,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast
      integer fluxgi, fluxgj, solngi, solngj
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj)
      double precision dx(0:1)
      integer location_index
      integer ratio_to_coarser(0:1)
c     Lower and upper corners of boundary box
      integer blower(0:1), bupper(0:1)

      double precision h
      integer i, ibeg, iend, igho, j, jbeg, jend, jgho
c     Fine grid indices inside one coarse grid.
      integer ip, jp
c     Fine grid indices for point diametrically opposite from (ip,jp).
      integer iq, jq
c     Weights associated with longtitudinal and transverse
c     (with respect to boundary normal) gradients.
      double precision tranwt, longwt

      if ( location_index .eq. 0 ) then
c        min i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = bupper(0)
         i = igho+1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do j=jbeg,jend,ratio_to_coarser(1)
            do jp=0,ratio_to_coarser(1)-1
               jq = ratio_to_coarser(1) - jp - 1
c     write(*,*) i, j, jq, jp
               xflux(i,j+jp)
     &            = longwt*xflux(i,j+jp)
     &            + tranwt*dc*( soln(i,j+jq) - soln(i,j+jp) )/h
            enddo
         enddo
      elseif ( location_index .eq. 1 ) then
c        max i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = blower(0)
         i = igho-1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do j=jbeg,jend,ratio_to_coarser(1)
            do jp=0,ratio_to_coarser(1)-1
               jq = ratio_to_coarser(1) - jp - 1
               xflux(igho,j+jp)
     &            = longwt*xflux(igho,j+jp)
     &            - tranwt*dc*( soln(i,j+jq) - soln(i,j+jp) )/h
            enddo
         enddo
      elseif ( location_index .eq. 2 ) then
c        min j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = bupper(1)
         j = jgho+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do ip=0,ratio_to_coarser(0)-1
               iq = ratio_to_coarser(0) - ip - 1
               yflux(i+ip,j)
     &            = longwt*yflux(i+ip,j)
     &            + tranwt*dc*( soln(i+iq,j) - soln(i+ip,j) )/h
            enddo
         enddo
      elseif ( location_index .eq. 3 ) then
c        max j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = blower(1)
         j = jgho-1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do ip=0,ratio_to_coarser(0)-1
               iq = ratio_to_coarser(0) - ip - 1
               yflux(i+ip,jgho)
     &            = longwt*yflux(i+ip,jgho)
     &            - tranwt*dc*( soln(i+iq,j) - soln(i+ip,j) )/h
            enddo
         enddo
      endif

      return
      end
c***********************************************************************
