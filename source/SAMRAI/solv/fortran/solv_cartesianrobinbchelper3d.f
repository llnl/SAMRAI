c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for Cartesian 3d Robin boundary conditions.
c
c***********************************************************************
c***********************************************************************
      subroutine settype1cells3d(
     &  data, 
     &  difirst, dilast, djfirst, djlast, dkfirst, dklast,
     &  a, b, g,
     &  ifirst, ilast, jfirst, jlast, kfirst, klast,
     &  ibeg, iend, jbeg, jend, kbeg, kend,
     &  face, ghos, inte, location,
     &  h, zerog )
c***********************************************************************
      implicit none
      integer difirst, dilast, djfirst, djlast, dkfirst, dklast
      double precision 
     &  data(difirst:dilast,djfirst:djlast,dkfirst:dklast)
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      double precision a(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision b(ifirst:ilast,jfirst:jlast,kfirst:klast)
      double precision g(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer ibeg, iend
      integer jbeg, jend
      integer kbeg, kend
      integer face, ghos, inte, location, zerog
      double precision h
      integer i, j, k
      if ( zerog .eq. 1 ) then
c        Assume g value of zero
         if ( location/2 .eq. 0 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(ghos,j,k) 
     &         = ( 0
     &         + data(inte,j,k)*( b(face,j,k)-0.5*h*a(face,j,k) ) )
     &         / ( b(face,j,k)+0.5*h*a(face,j,k) )
            enddo
            enddo
            enddo
         elseif ( location/2 .eq. 1 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(i,ghos,k) 
     &         = ( 0
     &         + data(i,inte,k)*( b(i,face,k)-0.5*h*a(i,face,k) ) )
     &         / ( b(i,face,k)+0.5*h*a(i,face,k) )
            enddo
            enddo
            enddo
         elseif ( location/2 .eq. 2 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(i,j,ghos) 
     &         = ( 0
     &         + data(i,j,inte)*( b(i,j,face)-0.5*h*a(i,j,face) ) )
     &         / ( b(i,j,face)+0.5*h*a(i,j,face) )
            enddo
            enddo
            enddo
         endif
      else
c        Assume finite g
         if ( location/2 .eq. 0 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(ghos,j,k) 
     &         = ( h*g(face,j,k)
     &         + data(inte,j,k)*( b(face,j,k)-0.5*h*a(face,j,k) ) )
     &         / ( b(face,j,k)+0.5*h*a(face,j,k) )
            enddo
            enddo
            enddo
         elseif ( location/2 .eq. 1 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(i,ghos,k) 
     &         = ( h*g(i,face,k)
     &         + data(i,inte,k)*( b(i,face,k)-0.5*h*a(i,face,k) ) )
     &         / ( b(i,face,k)+0.5*h*a(i,face,k) )
            enddo
            enddo
            enddo
         elseif ( location/2 .eq. 2 ) then
            do i=ibeg,iend
            do j=jbeg,jend
            do k=kbeg,kend
               data(i,j,ghos) 
     &         = ( h*g(i,j,face)
     &         + data(i,j,inte)*( b(i,j,face)-0.5*h*a(i,j,face) ) )
     &         / ( b(i,j,face)+0.5*h*a(i,j,face) )
            enddo
            enddo
            enddo
         endif
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine settype2cells3d(
     &  data, 
     &  difirst, dilast, djfirst, djlast, dkfirst, dklast,
     &  lower, upper, location )
c***********************************************************************
      implicit none
      integer difirst, dilast, djfirst, djlast, dkfirst, dklast
      double precision 
     &  data(difirst:dilast,djfirst:djlast,dkfirst:dklast)
      integer lower(0:2), upper(0:2), location
      integer i, j, k
      integer ibeg, iend, jbeg, jend, kbeg, kend
      if     ( location .eq. 0 ) then
c        min i min j edge, running along k
         i = lower(0)
         j = lower(1)
         kbeg = lower(2)
         kend = upper(2)
         do k=kbeg,kend
            data(i,j,k) = -data(i+1,j+1,k)
     &         + ( data(i+1,j,k) + data(i,j+1,k) )
         enddo
      elseif ( location .eq. 1 ) then
c        max i min j edge, running along k
         i = lower(0)
         j = lower(1)
         kbeg = lower(2)
         kend = upper(2)
         do k=kbeg,kend
            data(i,j,k) = -data(i-1,j+1,k)
     &         + ( data(i-1,j,k) + data(i,j+1,k) )
         enddo
      elseif ( location .eq. 2 ) then
c        min i max j edge, running along k
         i = lower(0)
         j = lower(1)
         kbeg = lower(2)
         kend = upper(2)
         do k=kbeg,kend
            data(i,j,k) = -data(i+1,j-1,k)
     &         + ( data(i+1,j,k) + data(i,j-1,k) )
         enddo
      elseif ( location .eq. 3 ) then
c        max i max j edge, running along k
         i = lower(0)
         j = lower(1)
         kbeg = lower(2)
         kend = upper(2)
         do k=kbeg,kend
            data(i,j,k) = -data(i-1,j-1,k)
     &         + ( data(i,j-1,k) + data(i-1,j,k) )
         enddo
      elseif ( location .eq. 4 ) then
c        min i min k edge, running along j
         i = lower(0)
         k = lower(2)
         jbeg = lower(1)
         jend = upper(1)
         do j=jbeg,jend
            data(i,j,k) = -data(i+1,j,k+1)
     &         + ( data(i+1,j,k) + data(i,j,k+1) )
      enddo
      elseif ( location .eq. 5 ) then
c        max i min k edge, running along j
         i = lower(0)
         k = lower(2)
         jbeg = lower(1)
         jend = upper(1)
         do j=jbeg,jend
            data(i,j,k) = -data(i-1,j,k+1)
     &         + ( data(i-1,j,k) + data(i,j,k+1) )
         enddo
      elseif ( location .eq. 6 ) then
c        min i max k edge, running along j
         i = lower(0)
         k = lower(2)
         jbeg = lower(1)
         jend = upper(1)
         do j=jbeg,jend
            data(i,j,k) = -data(i+1,j,k-1)
     &         + ( data(i+1,j,k) + data(i,j,k-1) )
         enddo
      elseif ( location .eq. 7 ) then
c        max i max k edge, running along j
         i = lower(0)
         k = lower(2)
         jbeg = lower(1)
         jend = upper(1)
         do j=jbeg,jend
            data(i,j,k) = -data(i-1,j,k-1)
     &         + ( data(i,j,k-1) + data(i-1,j,k) )
         enddo
      elseif ( location .eq. 8 ) then
c        min j min k edge, running along i
         j = lower(1)
         k = lower(2)
         ibeg = lower(0)
         iend = upper(0)
         do i=ibeg,iend
            data(i,j,k) = -data(i,j+1,k+1)
     &         + ( data(i,j,k+1) + data(i,j+1,k) )
         enddo
      elseif ( location .eq. 9 ) then
c        max j min k edge, running along i
         j = lower(1)
         k = lower(2)
         ibeg = lower(0)
         iend = upper(0)
         do i=ibeg,iend
            data(i,j,k) = -data(i,j-1,k+1)
     &         + ( data(i,j,k+1) + data(i,j-1,k) )
         enddo
      elseif ( location .eq. 10 ) then
c        min j max k edge, running along i
         j = lower(1)
         k = lower(2)
         ibeg = lower(0)
         iend = upper(0)
         do i=ibeg,iend
            data(i,j,k) = -data(i,j+1,k-1)
     &         + ( data(i,j,k-1) + data(i,j+1,k) )
         enddo
      elseif ( location .eq. 11 ) then
c        max j max k edge, running along i
         j = lower(1)
         k = lower(2)
         ibeg = lower(0)
         iend = upper(0)
         do i=ibeg,iend
            data(i,j,k) = -data(i,j-1,k-1)
     &         + ( data(i,j-1,k) + data(i,j,k-1) )
         enddo
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine settype3cells3d(
     &  data, 
     &  difirst, dilast, djfirst, djlast, dkfirst, dklast,
     &  lower, upper, location )
c***********************************************************************
      implicit none
      integer difirst, dilast, djfirst, djlast, dkfirst, dklast
      double precision 
     &  data(difirst:dilast,djfirst:djlast,dkfirst:dklast)
      integer lower(0:2), upper(0:2), location
      integer i, j, k
      if     ( location .eq. 0 ) then
c        min i min j min k node
         i = lower(0)
         j = lower(1)
         k = lower(2)
         data(i,j,k) = -data(i+1,j+1,k+1)
     &               + 2./3.*( data(i,j+1,k+1) 
     &                       + data(i+1,j,k+1) 
     &                       + data(i+1,j+1,k) )
      elseif ( location .eq. 1 ) then
c        max i min j min k node
         i = upper(0)
         j = lower(1)
         k = lower(2)
         data(i,j,k) = -data(i-1,j+1,k+1)
     &               + 2./3.*( data(i,j+1,k+1) 
     &                       + data(i-1,j,k+1) 
     &                       + data(i-1,j+1,k) )
      elseif ( location .eq. 2 ) then
c        min i max j min k node
         i = lower(0)
         j = upper(1)
         k = lower(2)
         data(i,j,k) = -data(i+1,j-1,k+1)
     &               + 2./3.*( data(i,j-1,k+1) 
     &                       + data(i+1,j,k+1) 
     &                       + data(i+1,j-1,k) )
      elseif ( location .eq. 3 ) then
c        max i max j min k node
         i = upper(0)
         j = upper(1)
         k = lower(2)
         data(i,j,k) = -data(i-1,j-1,k+1)
     &               + 2./3.*( data(i,j-1,k+1) 
     &                       + data(i-1,j,k+1) 
     &                       + data(i-1,j-1,k) )
      elseif ( location .eq. 4 ) then
c        min i min j max k node
         i = lower(0)
         j = lower(1)
         k = upper(2)
         data(i,j,k) = -data(i+1,j+1,k-1)
     &               + 2./3.*( data(i,j+1,k-1) 
     &                       + data(i+1,j,k-1) 
     &                       + data(i+1,j+1,k) )
      elseif ( location .eq. 5 ) then
c        max i min j max k node
         i = upper(0)
         j = lower(1)
         k = upper(2)
         data(i,j,k) = -data(i-1,j+1,k-1)
     &               + 2./3.*( data(i,j+1,k-1) 
     &                       + data(i-1,j,k-1) 
     &                       + data(i-1,j+1,k) )
      elseif ( location .eq. 6 ) then
c        min i max j max k node
         i = lower(0)
         j = upper(1)
         k = upper(2)
         data(i,j,k) = -data(i+1,j-1,k-1)
     &               + 2./3.*( data(i,j-1,k-1) 
     &                       + data(i+1,j,k-1) 
     &                       + data(i+1,j-1,k) )
      elseif ( location .eq. 7 ) then
c        max i max j max k node
         i = upper(0)
         j = upper(1)
         k = upper(2)
         data(i,j,k) = -data(i-1,j-1,k-1)
     &               + 2./3.*( data(i,j-1,k-1) 
     &                       + data(i-1,j,k-1) 
     &                       + data(i-1,j-1,k) )
      endif
      return
      end
c***********************************************************************
