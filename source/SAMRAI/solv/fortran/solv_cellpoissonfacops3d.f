c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for scalar Poisson FAC operator.
c
c***********************************************************************
c***********************************************************************
      subroutine compfluxvardc3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &xdc , ydc , zdc, dcgi, dcgj , dcgk ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer dcgi, dcgj, dcgk, fluxgi, fluxgj, fluxgk, 
     &        solngi, solngj, solngk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision zdc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+1+dcgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision dx(0:2)

      double precision dxi, dyi, dzi
      integer i, j, k

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)

      do k=kfirst,klast
      do j=jfirst,jlast
      do i=ifirst,ilast+1
         xflux(i,j,k) = dxi*xdc(i,j,k)*( soln(i,j,k) - soln(i-1,j,k) )
      enddo
      enddo
      enddo
      do k=kfirst,klast
      do j=jfirst,jlast+1
      do i=ifirst,ilast
         yflux(i,j,k) = dyi*ydc(i,j,k)*( soln(i,j,k) - soln(i,j-1,k) )
      enddo
      enddo
      enddo
      do k=kfirst,klast+1
      do j=jfirst,jlast
      do i=ifirst,ilast
         zflux(i,j,k) = dzi*zdc(i,j,k)*( soln(i,j,k) - soln(i,j,k-1) )
      enddo
      enddo
      enddo

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compfluxcondc3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &dc ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer fluxgi, fluxgj, fluxgk, 
     &        solngi, solngj, solngk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision dx(0:2)

      double precision dxi, dyi, dzi
      integer i, j, k

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)

      do k=kfirst,klast
      do j=jfirst,jlast
      do i=ifirst,ilast+1
         xflux(i,j,k) = dxi*dc*( soln(i,j,k) - soln(i-1,j,k) )
      enddo
      enddo
      enddo
      do k=kfirst,klast
      do j=jfirst,jlast+1
      do i=ifirst,ilast
         yflux(i,j,k) = dyi*dc*( soln(i,j,k) - soln(i,j-1,k) )
      enddo
      enddo
      enddo
      do k=kfirst,klast+1
      do j=jfirst,jlast
      do i=ifirst,ilast
         zflux(i,j,k) = dzi*dc*( soln(i,j,k) - soln(i,j,k-1) )
      enddo
      enddo
      enddo

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxvardcvarsf3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &xdc , ydc , zdc, dcgi, dcgj , dcgk ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &scalar_field , sfgi, sfgj , sfgk ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer dcgi, dcgj, dcgk, fluxgi, fluxgj, fluxgk, 
     &        rhsgi, rhsgj, rhsgk, solngi, solngj, solngk,
     &        sfgi, sfgj, sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision zdc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+1+dcgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj,
     &                              kfirst-sfgk:klast+sfgk)
      double precision dx(0:2)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi, dzi, dxi2, dyi2, dzi2
      double precision dudr
      double precision rcoef
      integer i, j, k
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)
      dxi2 = dxi*dxi
      dyi2 = dyi*dyi
      dzi2 = dzi*dzi
      rcoef = 1.0

      maxres = 0.0

      do k=kfirst,klast
      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (ifirst+j+k)-((ifirst+j+k)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
             residual
     &          = rhs(i,j,k)
     &          - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &            + dyi*( yflux(i,j+1,k) - yflux(i,j,k) )
     &            + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &          - scalar_field(i,j,k)*soln(i,j,k)
             dudr = 1./( ( dxi2*( xdc(i+1,j,k)+xdc(i,j,k) )
     &                   + dyi2*( ydc(i,j+1,k)+ydc(i,j,k) )
     &                   + dzi2*( zdc(i,j,k+1)+zdc(i,j,k) ) )
     &                   - scalar_field(i,j,k) )
             du = -residual*dudr
             soln(i,j,k) = soln(i,j,k) + du*rcoef
             if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxcondcvarsf3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &dc ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &scalar_field , sfgi, sfgj , sfgk ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer fluxgi, fluxgj, fluxgk, 
     &        rhsgi, rhsgj, rhsgk, solngi, solngj, solngk,
     &        sfgi, sfgj, sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj,
     &                              kfirst-sfgk:klast+sfgk)
      double precision dx(0:2)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi, dzi, dxi2, dyi2, dzi2
      double precision dudr
      double precision rcoef
      integer i, j, k
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)
      dxi2 = dxi*dxi
      dyi2 = dyi*dyi
      dzi2 = dzi*dzi
      rcoef = 1.0

      maxres = 0.0

      do k=kfirst,klast
      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (ifirst+j+k)-((ifirst+j+k)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
             residual
     &          = rhs(i,j,k)
     &          - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &            + dyi*( yflux(i,j+1,k) - yflux(i,j,k) )
     &            + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &          - scalar_field(i,j,k)*soln(i,j,k)
             dudr = 1./( ( dxi2*( dc + dc )
     &                   + dyi2*( dc + dc )
     &                   + dzi2*( dc + dc ) )
     &                   - scalar_field(i,j,k) )
             du = -residual*dudr
             soln(i,j,k) = soln(i,j,k) + du*rcoef
             if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxvardcconsf3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &xdc , ydc , zdc, dcgi, dcgj , dcgk ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &scalar_field ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer dcgi, dcgj, dcgk, fluxgi, fluxgj, fluxgk, 
     &        rhsgi, rhsgj, rhsgk, solngi, solngj, solngk,
     &        sfgi, sfgj, sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision zdc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+1+dcgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field
      double precision dx(0:2)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi, dzi, dxi2, dyi2, dzi2
      double precision dudr
      double precision rcoef
      integer i, j, k
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)
      dxi2 = dxi*dxi
      dyi2 = dyi*dyi
      dzi2 = dzi*dzi
      rcoef = 1.0

      maxres = 0.0

      do k=kfirst,klast
      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (ifirst+j+k)-((ifirst+j+k)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
             residual
     &          = rhs(i,j,k)
     &          - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &            + dyi*( yflux(i,j+1,k) - yflux(i,j,k) )
     &            + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &          - scalar_field*soln(i,j,k)
             dudr = 1./( ( dxi2*( xdc(i+1,j,k)+xdc(i,j,k) )
     &                   + dyi2*( ydc(i,j+1,k)+ydc(i,j,k) )
     &                   + dzi2*( zdc(i,j,k+1)+zdc(i,j,k) ) )
     &                   - scalar_field )
             du = -residual*dudr
             soln(i,j,k) = soln(i,j,k) + du*rcoef
             if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine rbgswithfluxmaxcondcconsf3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &dc ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &scalar_field ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &dx ,
     &offset, maxres )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer fluxgi, fluxgj, fluxgk, 
     &        rhsgi, rhsgj, rhsgk, solngi, solngj, solngk,
     &        sfgi, sfgj, sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field
      double precision dx(0:2)
      integer offset
      double precision maxres

      double precision residual, du
      double precision dxi, dyi, dzi, dxi2, dyi2, dzi2
      double precision dudr
      double precision rcoef
      integer i, j, k
      integer ioffset

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)
      dxi2 = dxi*dxi
      dyi2 = dyi*dyi
      dzi2 = dzi*dzi
      rcoef = 1.0

      maxres = 0.0

      do k=kfirst,klast
      do j=jfirst,jlast
c        offset must be 0 (red) or 1 (black)
         if ( (ifirst+j+k)-((ifirst+j+k)/2*2) .ne. offset ) then
            ioffset = 1
         else
            ioffset = 0
         endif
         do i=ifirst+ioffset,ilast,2
             residual
     &          = rhs(i,j,k)
     &          - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &            + dyi*( yflux(i,j+1,k) - yflux(i,j,k) )
     &            + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &          - scalar_field*soln(i,j,k)
             dudr = 1./( ( dxi2*( dc + dc )
     &                   + dyi2*( dc + dc )
     &                   + dzi2*( dc + dc ) )
     &                   - scalar_field )
             du = -residual*dudr
             soln(i,j,k) = soln(i,j,k) + du*rcoef
             if ( maxres .lt. abs(residual) ) maxres = abs(residual)
         enddo
      enddo
      enddo
      return
      end
c***********************************************************************

c***********************************************************************
      subroutine compresvarsca3d(
     &xflux , yflux , zflux , fluxgi, fluxgj , fluxgk ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &residual , residualgi, residualgj , residualgk ,
     &scalar_field , sfgi, sfgj , sfgk ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst , klast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast , kfirst , klast
      integer fluxgi, fluxgj, fluxgk , rhsgi, rhsgj, rhsgk ,
     &        residualgi, residualgj, residualgk ,
     &        solngi, solngj, solngk , sfgi, sfgj , sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision residual(ifirst-rhsgi:ilast+rhsgi,
     &                          jfirst-rhsgj:jlast+rhsgj,
     &                          kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field(ifirst-sfgi:ilast+sfgi,
     &                              jfirst-sfgj:jlast+sfgj,
     &                              kfirst-sfgk:klast+sfgk)
      double precision dx(0:2)

      double precision dxi, dyi , dzi
      integer i, j, k

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)

      do k=kfirst,klast
      do j=jfirst,jlast
      do i=ifirst,ilast
         residual(i,j,k)
     &      = rhs(i,j,k)
     &      - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &        + dyi*( yflux(i,j+1,k) - yflux(i,j,k) ) 
     &        + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &      - scalar_field(i,j,k)*soln(i,j,k)
      enddo
      enddo
      enddo
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine compresconsca3d(
     &xflux , yflux , zflux , fluxgi, fluxgj , fluxgk ,
     &rhs , rhsgi, rhsgj , rhsgk ,
     &residual , residualgi, residualgj , residualgk ,
     &scalar_field ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst , klast ,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast , kfirst , klast
      integer fluxgi, fluxgj, fluxgk , rhsgi, rhsgj, rhsgk ,
     &        residualgi, residualgj, residualgk ,
     &        solngi, solngj, solngk , sfgi, sfgj , sfgk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision rhs(ifirst-rhsgi:ilast+rhsgi,
     &                     jfirst-rhsgj:jlast+rhsgj,
     &                     kfirst-rhsgk:klast+rhsgk)
      double precision residual(ifirst-rhsgi:ilast+rhsgi,
     &                          jfirst-rhsgj:jlast+rhsgj,
     &                          kfirst-rhsgk:klast+rhsgk)
      double precision scalar_field
      double precision dx(0:2)

      double precision dxi, dyi , dzi
      integer i, j, k

      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)

      do k=kfirst,klast
      do j=jfirst,jlast
      do i=ifirst,ilast
         residual(i,j,k)
     &      = rhs(i,j,k)
     &      - ( dxi*( xflux(i+1,j,k) - xflux(i,j,k) )
     &        + dyi*( yflux(i,j+1,k) - yflux(i,j,k) ) 
     &        + dzi*( zflux(i,j,k+1) - zflux(i,j,k) ) )
     &      - scalar_field*soln(i,j,k)
      enddo
      enddo
      enddo
      return
      end
c***********************************************************************

c***********************************************************************
      subroutine ewingfixfluxvardc3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &xdc , ydc , zdc, dcgi, dcgj , dcgk ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &location_index ,
     &ratio_to_coarser ,
     &blower, bupper,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer dcgi, dcgj, dcgk, fluxgi, fluxgj, fluxgk, 
     &        solngi, solngj, solngk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision xdc(ifirst-dcgi:ilast+1+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision ydc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+1+dcgj,
     &                     kfirst-dcgk:klast+dcgk)
      double precision zdc(ifirst-dcgi:ilast+dcgi,
     &                     jfirst-dcgj:jlast+dcgj,
     &                     kfirst-dcgk:klast+1+dcgk)
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision dx(0:2)
      integer location_index
      integer ratio_to_coarser(0:2)
c     Lower and upper corners of boundary box
      integer blower(0:2), bupper(0:2)

      double precision h
      integer i, ibeg, iend, ibnd, igho,
     &        j, jbeg, jend, jbnd, jgho,
     &        k, kbeg, kend, kbnd, kgho
c     Fine grid indices inside one coarse grid.
      integer ip, jp, kp
c     Fine grid indices for point diametrically opposite from (ip,jp).
      integer iq, jq, kq
c     Weights associated with longtitudinal and transverse
c     (with respect to boundary normal) gradients.
      double precision tranwt, longwt

      if ( location_index .eq. 0 ) then
c        min i face
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         i = bupper(0)+1
         ibnd = bupper(0)+1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do k=kbeg,kend,ratio_to_coarser(2)
            do j=jbeg,jend,ratio_to_coarser(1)
               do kp=0,ratio_to_coarser(2)-1
                  kq = ratio_to_coarser(2) - kp - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     xflux(ibnd,j+jp,k+kp)
     &                  = longwt*xflux(ibnd,j+jp,k+kp)
     &                  + tranwt*xdc(ibnd,j+jp,k+kp)*( 
     &                    soln(i,j+jq,k+kq) - soln(i,j+jp,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 1 ) then
c        max i face
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         i = blower(0)-1
         ibnd = blower(0)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do k=kbeg,kend,ratio_to_coarser(2)
            do j=jbeg,jend,ratio_to_coarser(1)
               do kp=0,ratio_to_coarser(2)-1
                  kq = ratio_to_coarser(2) - kp - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     xflux(ibnd,j+jp,k+kp)
     &                  = longwt*xflux(ibnd,j+jp,k+kp)
     &                  - tranwt*xdc(ibnd,j+jp,k+kp)*( 
     &                    soln(i,j+jq,k+kq) - soln(i,j+jp,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 2 ) then
c        min j face
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         j = bupper(1)+1
         jbnd = bupper(1)+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do k=kbeg,kend,ratio_to_coarser(2)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do kp=0,ratio_to_coarser(2)-1
                     kq = ratio_to_coarser(2) - kp - 1
                     yflux(i+ip,jbnd,k+kp)
     &                  = longwt*yflux(i+ip,jbnd,k+kp)
     &                  + tranwt*ydc(i+ip,jbnd,k+kp)*( 
     &                    soln(i+iq,j,k+kq) - soln(i+ip,j,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 3 ) then
c        max j face
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         j = blower(1)-1
         jbnd = blower(1)
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do k=kbeg,kend,ratio_to_coarser(2)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do kp=0,ratio_to_coarser(2)-1
                     kq = ratio_to_coarser(2) - kp - 1
                     yflux(i+ip,jbnd,k+kp)
     &                  = longwt*yflux(i+ip,jbnd,k+kp)
     &                  - tranwt*ydc(i+ip,jbnd,k+kp)*( 
     &                    soln(i+iq,j,k+kq) - soln(i+ip,j,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 4 ) then
c        min k face
         tranwt = 1.0/(1+ratio_to_coarser(2))
         longwt = 2*tranwt
         h = dx(2)
         k = bupper(2)+1
         kbnd = bupper(2)+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do j=jbeg,jend,ratio_to_coarser(1)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     zflux(i+ip,j+jp,kbnd)
     &                  = longwt*zflux(i+ip,j+jp,kbnd)
     &                  + tranwt*zdc(i+ip,j+jp,kbnd)*( 
     &                    soln(i+iq,j+jq,k) - soln(i+ip,j+jp,k) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 5 ) then
c        max k face
         tranwt = 1.0/(1+ratio_to_coarser(2))
         longwt = 2*tranwt
         h = dx(2)
         k = blower(2)-1
         kbnd = blower(2)
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do j=jbeg,jend,ratio_to_coarser(1)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     zflux(i+ip,j+jp,kbnd)
     &                  = longwt*zflux(i+ip,j+jp,kbnd)
     &                  - tranwt*zdc(i+ip,j+jp,kbnd)*( 
     &                    soln(i+iq,j+jq,k) - soln(i+ip,j+jp,k) )/h
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine ewingfixfluxcondc3d(
     &xflux , yflux , zflux, fluxgi, fluxgj , fluxgk,
     &dc ,
     &soln , solngi, solngj , solngk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast ,
     &location_index ,
     &ratio_to_coarser ,
     &blower, bupper,
     &dx )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer fluxgi, fluxgj, fluxgk, 
     &        solngi, solngj, solngk
      double precision xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision yflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+1+fluxgj,
     &                       kfirst-fluxgk:klast+fluxgk)
      double precision zflux(ifirst-fluxgi:ilast+fluxgi,
     &                       jfirst-fluxgj:jlast+fluxgj,
     &                       kfirst-fluxgk:klast+1+fluxgk)
      double precision dc
      double precision soln(ifirst-solngi:ilast+solngi,
     &                      jfirst-solngj:jlast+solngj,
     &                      kfirst-solngk:klast+solngk)
      double precision dx(0:2)
      integer location_index
      integer ratio_to_coarser(0:2)
c     Lower and upper corners of boundary box
      integer blower(0:2), bupper(0:2)

      double precision h
      integer i, ibeg, iend, ibnd, igho,
     &        j, jbeg, jend, jbnd, jgho,
     &        k, kbeg, kend, kbnd, kgho
c     Fine grid indices inside one coarse grid.
      integer ip, jp, kp
c     Fine grid indices for point diametrically opposite from (ip,jp).
      integer iq, jq, kq
c     Weights associated with longtitudinal and transverse
c     (with respect to boundary normal) gradients.
      double precision tranwt, longwt

      if ( location_index .eq. 0 ) then
c        min i face
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         i = bupper(0)+1
         ibnd = bupper(0)+1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do k=kbeg,kend,ratio_to_coarser(2)
            do j=jbeg,jend,ratio_to_coarser(1)
               do kp=0,ratio_to_coarser(2)-1
                  kq = ratio_to_coarser(2) - kp - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     xflux(ibnd,j+jp,k+kp)
     &                  = longwt*xflux(ibnd,j+jp,k+kp)
     &                  + tranwt*dc*( 
     &                    soln(i,j+jq,k+kq) - soln(i,j+jp,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 1 ) then
c        max i face
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         i = blower(0)-1
         ibnd = blower(0)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do k=kbeg,kend,ratio_to_coarser(2)
            do j=jbeg,jend,ratio_to_coarser(1)
               do kp=0,ratio_to_coarser(2)-1
                  kq = ratio_to_coarser(2) - kp - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     xflux(ibnd,j+jp,k+kp)
     &                  = longwt*xflux(ibnd,j+jp,k+kp)
     &                  - tranwt*dc*( 
     &                    soln(i,j+jq,k+kq) - soln(i,j+jp,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 2 ) then
c        min j face
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         j = bupper(1)+1
         jbnd = bupper(1)+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do k=kbeg,kend,ratio_to_coarser(2)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do kp=0,ratio_to_coarser(2)-1
                     kq = ratio_to_coarser(2) - kp - 1
                     yflux(i+ip,jbnd,k+kp)
     &                  = longwt*yflux(i+ip,jbnd,k+kp)
     &                  + tranwt*dc*( 
     &                    soln(i+iq,j,k+kq) - soln(i+ip,j,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 3 ) then
c        max j face
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         j = blower(1)-1
         jbnd = blower(1)
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         kbeg = max(blower(2),kfirst)
         kend = min(bupper(2),klast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do k=kbeg,kend,ratio_to_coarser(2)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do kp=0,ratio_to_coarser(2)-1
                     kq = ratio_to_coarser(2) - kp - 1
                     yflux(i+ip,jbnd,k+kp)
     &                  = longwt*yflux(i+ip,jbnd,k+kp)
     &                  - tranwt*dc*( 
     &                    soln(i+iq,j,k+kq) - soln(i+ip,j,k+kp) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 4 ) then
c        min k face
         tranwt = 1.0/(1+ratio_to_coarser(2))
         longwt = 2*tranwt
         h = dx(2)
         k = bupper(2)+1
         kbnd = bupper(2)+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do j=jbeg,jend,ratio_to_coarser(1)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     zflux(i+ip,j+jp,kbnd)
     &                  = longwt*zflux(i+ip,j+jp,kbnd)
     &                  + tranwt*dc*( 
     &                    soln(i+iq,j+jq,k) - soln(i+ip,j+jp,k) )/h
                  enddo
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 5 ) then
c        max k face
         tranwt = 1.0/(1+ratio_to_coarser(2))
         longwt = 2*tranwt
         h = dx(2)
         k = blower(2)-1
         kbnd = blower(2)
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do i=ibeg,iend,ratio_to_coarser(0)
            do j=jbeg,jend,ratio_to_coarser(1)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  do jp=0,ratio_to_coarser(1)-1
                     jq = ratio_to_coarser(1) - jp - 1
                     zflux(i+ip,j+jp,kbnd)
     &                  = longwt*zflux(i+ip,j+jp,kbnd)
     &                  - tranwt*dc*( 
     &                    soln(i+iq,j+jq,k) - soln(i+ip,j+jp,k) )/h
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
c***********************************************************************
