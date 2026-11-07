c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   This code is post-processed fortran that has been cleaned up.
c

c     ----------------------------------------------------------------------

      subroutine bcmultiblock(
     *    x_ilo, x_ihi, x_jlo, x_jhi, x_klo, x_khi, x, 
     *    d_ilo, d_ihi, d_jlo, d_jhi, d_klo, d_khi, d,  nd,
     *    dlo, dhi,  
     *    glo, ghi,  
     *    lo,  hi,
     *    dir, side )
c
c     Extrapolate boundary conditions
c
c     x = nodal positions
c     d = state array
c     dlo, dhi = domain bounds
c     glo, ghi = ghost box bounds
c     lo,  hi  = patch bounds
c
      implicit none

c     ... input
      integer x_ilo,x_ihi,x_jlo,x_jhi,x_klo,x_khi
      real*8 x(x_ilo:x_ihi+1,x_jlo:x_jhi+1,x_klo:x_khi+1,3)

      integer nd
      integer d_ilo,d_ihi,d_jlo,d_jhi,d_klo,d_khi
      real*8 d(d_ilo:d_ihi,d_jlo:d_jhi,d_klo:d_khi,nd)

      integer dlo(3), dhi(3)  ! the index space high and low
      integer glo(3), ghi(3)  ! the zonal grown high and low
      integer lo(3),  hi(3)   ! the interior high and low
      integer dir, side       ! the direction we are filling in and low(0) or high(1)

c     ... local
      integer i, j, k, n
      integer ii, jj, kk

c     -------------------------- begin here ---------------------------

c
c first we sweep out i-1 ghosts then i+1 ghosts, note to get at the edges
c and corners, we make sure that we populate the ghosts used in the extrapolation
c with valid data, hence the changing expanding bounds when we then go to the j, and k
c ghosts
c
c     for zones we just extrapolate piece-wise constant (even bc)
c     for nodal positions we extrapolate linearly (preserves mesh size and spacing)
c

c     ilo conditions
      if (dir.eq.0 .and. side.eq.0) then
         
         do k = max(glo(3),dlo(3)), min(ghi(3),dhi(3)) 
            do j = max(glo(2),dlo(2)), min(ghi(2),dhi(2)) 
               do i = glo(1), dlo(1)-1
                  do n = 1, nd
                     d(i,j,k,n) = d(dlo(1),j,k,n)
                  end do
               end do
            end do 
         end do 
         
         do k = max(glo(3),dlo(3)), min(ghi(3)+1,dhi(3)+1) 
            do j = max(glo(2),dlo(2)), min(ghi(2)+1,dhi(2)+1) 
               do i = glo(1), dlo(1)-1
                  ii = 2*(dlo(1))-i
                  x(i,j,k,1) = 2.*x(dlo(1),j,k,1) - x(ii,j,k,1)
                  x(i,j,k,2) = 2.*x(dlo(1),j,k,2) - x(ii,j,k,2) 
                  x(i,j,k,3) = 2.*x(dlo(1),j,k,3) - x(ii,j,k,3) 
               end do
            end do 
         end do 
         
c     ihi conditions
      else if (dir.eq.0 .and. side.eq.1) then
         
         do k = max(glo(3),dlo(3)), min(ghi(3),dhi(3)) 
            do j = max(glo(2),dlo(2)), min(ghi(2),dhi(2)) 
               do i = dhi(1)+1, ghi(1)
                  do n = 1, nd
                     d(i,j,k,n) = d(dhi(1),j,k,n)
                  end do
               end do
            end do 
         end do 
         
         do k = max(glo(3),dlo(3)), min(ghi(3)+1,dhi(3)+1) 
            do j = max(glo(2),dlo(2)), min(ghi(2)+1,dhi(2)+1) 
               do i = dhi(1)+2,ghi(1)+1
                  ii = 2*(dhi(1)+1)-i
                  x(i,j,k,1) = 2.*x(dhi(1)+1,j,k,1) - x(ii,j,k,1)
                  x(i,j,k,2) = 2.*x(dhi(1)+1,j,k,2) - x(ii,j,k,2) 
                  x(i,j,k,3) = 2.*x(dhi(1)+1,j,k,3) - x(ii,j,k,3) 
               end do
            end do 
         end do 
         
      
c     jlo conditions
      else if (dir.eq.1 .and. side.eq.0) then
         do k = max(glo(3),dlo(3)), min(ghi(3),dhi(3)) 
            do j = glo(2), dlo(2)-1
               do i = glo(1), ghi(1)
                  do n = 1, nd
                     d(i,j,k,n) = d(i,dlo(2),k,n)
                  end do
               end do
            end do
         end do 
         
         do k = max(glo(3),dlo(3)), min(ghi(3)+1,dhi(3)+1) 
            do j = glo(2), dlo(2)-1
               do i = glo(1), ghi(1)+1
                  jj = 2*(dlo(2))-j
                  x(i,j,k,1) = 2.*x(i,dlo(2),k,1) - x(i,jj,k,1)
                  x(i,j,k,2) = 2.*x(i,dlo(2),k,2) - x(i,jj,k,2)
                  x(i,j,k,3) = 2.*x(i,dlo(2),k,3) - x(i,jj,k,3)
               end do
            end do
         end do 
         
c     jhi conditions
      else if (dir.eq.1 .and. side.eq.1) then
         
         do k = max(glo(3),dlo(3)), min(ghi(3),dhi(3)) 
            do j = dhi(2)+1, ghi(2)
               do i = glo(1), ghi(1)
                  do n = 1, nd
                     d(i,j,k,n) = d(i,dhi(2),k,n)
                  end do
               end do
            end do
         end do 
         
         do k = max(glo(3),dlo(3)), min(ghi(3)+1,dhi(3)+1) 
            do j = dhi(2)+2, ghi(2)+1
               do i = glo(1), ghi(1)+1
                  jj = 2*(dhi(2)+1)-j
                  x(i,j,k,1) = 2.*x(i,dhi(2)+1,k,1) - x(i,jj,k,1)
                  x(i,j,k,2) = 2.*x(i,dhi(2)+1,k,2) - x(i,jj,k,2)
                  x(i,j,k,3) = 2.*x(i,dhi(2)+1,k,3) - x(i,jj,k,3)
               end do
            end do
         end do 
         
c     klo conditions
      else if (dir.eq.2 .and. side.eq.0) then
         
         do k = glo(3), dlo(3)-1
            do j = glo(2), ghi(2)
               do i = glo(1), ghi(1)
                  do n = 1, nd
                     d(i,j,k,n) = d(i,j,dlo(3),n)
                  end do
               end do
            end do
         end do
         
         do k = glo(3), dlo(3)-1
            do j = glo(2), ghi(2)+1
               do i = glo(1), ghi(1)+1
                  kk = 2*(dlo(3))-k
                  x(i,j,k,1) = 2.*x(i,j,dlo(3),1) - x(i,j,kk,1)
                  x(i,j,k,2) = 2.*x(i,j,dlo(3),2) - x(i,j,kk,2)
                  x(i,j,k,3) = 2.*x(i,j,dlo(3),3) - x(i,j,kk,3)
               end do
            end do
         end do
         
c     khi conditions
      else if (dir.eq.2 .and. side.eq.1) then
         
         do k = dhi(3)+1, ghi(3)
            do j = glo(2), ghi(2)
               do i = glo(1), ghi(1)
                  do n = 1, nd
                     d(i,j,k,n) = d(i,j,dhi(3),n)
                  end do
               end do
            end do
         end do
         
         do k = dhi(3)+2, ghi(3)+1
            do j = glo(2), ghi(2)+1
               do i = glo(1), ghi(1)+1
                  kk = 2*(dhi(3)+1)-k
                  x(i,j,k,1) = 2.*x(i,j,dhi(3)+1,1) - x(i,j,kk,1)
                  x(i,j,k,2) = 2.*x(i,j,dhi(3)+1,2) - x(i,j,kk,2)
                  x(i,j,k,3) = 2.*x(i,j,dhi(3)+1,3) - x(i,j,kk,3)
               end do
            end do
         end do
         
      end if
c
c     end of subroutine
c
      return
      end

c---------------------------------------------------------------

      subroutine getmyfacebdry(
     &     d_ilo, d_ihi, d_jlo, d_jhi, d_klo, d_khi, d,  nd,
     &     ifirst0, ilast0,
     $     ifirst1, ilast1,
     $     ifirst2, ilast2,      ! interior bounds
     &     ibeg0, iend0,
     $     ibeg1, iend1,
     $     ibeg2, iend2,         ! bounds to fill
     &     bdry_loc )
c
c     fill in the faces
c
c     input
c
      implicit none
      integer nd
      integer d_ilo,d_ihi,d_jlo,d_jhi,d_klo,d_khi
      real*8 d(d_ilo:d_ihi,d_jlo:d_jhi,d_klo:d_khi,nd)

      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer bdry_loc
c
c     locals
c
      include 'bc_common.f'
      integer ic2, ic1, ic0, ict0, ict1, ict2, k
      integer ipivot0, ipivot1, ipivot2
c
c     -------------------
c
      if ( (bdry_loc .eq. XLO).or.(bdry_loc .eq. XHI) ) then

         if ( bdry_loc .eq. XLO ) then
            ipivot0 = ifirst0
         else
            ipivot0 = ilast0
         endif

         do ic2=ibeg2,iend2
            ict2 = ic2
            do ic1=ibeg1,iend1
               ict1 = ic1
               do ic0=ibeg0,iend0
                  ict0 = ipivot0
                  do k=0,nd-1
                     d(ic0,ic1,ic2,k) = d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo

      else if ( (bdry_loc .eq. YLO).or.(bdry_loc .eq. YHI) ) then

         if ( bdry_loc .eq. YLO ) then
            ipivot1 = ifirst1
         else
            ipivot1 = ilast1
         endif

         do ic2 = ibeg2,iend2
            ict2 = ic2
            do ic1 = ibeg1,iend1
               ict1 = ipivot1
               do ic0 = ibeg0,iend0
                  ict0 = ic0
                  do k = 1, nd
                     d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo

      else if ( (bdry_loc .eq. ZLO).or.(bdry_loc .eq. ZHI) ) then

         if ( bdry_loc .eq. ZLO ) then
            ipivot2 = ifirst2
         else
            ipivot2 = ilast2
         endif

         do ic2 = ibeg2,iend2
            ict2 = ipivot2
            do ic1 = ibeg1,iend1
               ict1 = ic1
               do ic0 = ibeg0,iend0
                  ict0 = ic0
                  do k = 1, nd
                     d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo

      endif
c
c     end of subroutine
c
      return
      end

c---------------------------------------------------------------

      subroutine getmyedgebdry(
     &     d_ilo, d_ihi, d_jlo, d_jhi, d_klo, d_khi, d,  nd,
     &     ifirst0, ilast0,
     $     ifirst1, ilast1,
     $     ifirst2, ilast2,      ! interior bounds
     &     ibeg0, iend0,
     $     ibeg1, iend1,
     $     ibeg2, iend2,         ! bounds to fill
     &     bdry_loc )
c
c     fill in the edges
c
c     input
c
      implicit none
      integer nd
      integer d_ilo,d_ihi,d_jlo,d_jhi,d_klo,d_khi
      real*8 d(d_ilo:d_ihi,d_jlo:d_jhi,d_klo:d_khi,nd)

      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer bdry_loc
c
c     locals
c
      include 'bc_common.f'
      integer ic2, ic1, ic0, ict0, ict1, ict2
      integer k
      integer ipivot
c
c     ------------------- the i edges
c
      if (  (bdry_loc .eq. Y0Z0).or.
     $     (bdry_loc .eq. Y1Z0).or.
     &     (bdry_loc .eq. Y0Z1).or.
     $     (bdry_loc .eq. Y1Z1) ) then

         if ( bdry_loc .eq. Y0Z0 ) then
            ipivot = ifirst1
         else if ( bdry_loc .eq. Y1Z0 ) then
            ipivot = ilast1
         else if ( bdry_loc .eq. Y0Z1 ) then
            ipivot = ifirst1
         else
            ipivot = ilast1
         endif

         do ic2 = ibeg2,iend2
            ict2 = ic2
            do ic1 = ibeg1,iend1
               ict1 = ipivot
               do ic0 = ibeg0,iend0
                  ict0 = ic0
                  do k = 1, nd
                     d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo
c
c     ------------------- the j edges
c
      else if ( 
     $        (bdry_loc .eq. X0Z0).or.
     $        (bdry_loc .eq. X0Z1).or.
     &        (bdry_loc .eq. X1Z0).or.
     $        (bdry_loc .eq. X1Z1) ) then

         if ( bdry_loc .eq. X0Z0 ) then
            ipivot = ifirst0
         else if ( bdry_loc .eq. X0Z1 ) then
            ipivot = ifirst0
         else if ( bdry_loc .eq. X1Z0 ) then
            ipivot = ilast0
         else
            ipivot = ilast0
         endif

         do ic2 = ibeg2,iend2
            ict2 = ic2
            do ic1 = ibeg1,iend1
               ict1 = ic1
               do ic0 = ibeg0,iend0
                  ict0 = ipivot
                  do k =1, nd
                     d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo

c
c     ------------------- the k edges
c
      else if ( 
     $        (bdry_loc .eq. X0Y0).or.
     $        (bdry_loc .eq. X1Y0).or.
     &        (bdry_loc .eq. X0Y1).or.
     $        (bdry_loc .eq. X1Y1) ) then

         if ( bdry_loc .eq. X0Y0 ) then
            ipivot = ifirst0
         else if ( bdry_loc .eq. X1Y0 ) then
            ipivot = ilast0
         else if ( bdry_loc .eq. X0Y1 ) then
            ipivot = ifirst0
         else
            ipivot = ilast0
         endif

         do ic2 = ibeg2,iend2
            ict2 = ic2
            do ic1 = ibeg1,iend1
               ict1 = ic1
               do ic0 = ibeg0,iend0
                  ict0 = ipivot
                  do k= 1, nd
                     d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
                  enddo
               enddo
            enddo
         enddo
      
      endif
c
c     end of subroutine
c
      return
      end


c-------------------------------------------------------------------------


      subroutine getmynodebdry(
     &     d_ilo, d_ihi, d_jlo, d_jhi, d_klo, d_khi, d,  nd,
     &     ifirst0, ilast0,
     $     ifirst1, ilast1,
     $     ifirst2, ilast2,      ! interior bounds
     &     ibeg0, iend0,
     $     ibeg1, iend1,
     $     ibeg2, iend2,         ! bounds to fill
     &     bdry_loc )
c
c     fill in the nodes (corners)
c
      implicit none
      integer nd
      integer d_ilo,d_ihi,d_jlo,d_jhi,d_klo,d_khi
      real*8 d(d_ilo:d_ihi,d_jlo:d_jhi,d_klo:d_khi,nd)

      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer bdry_loc
c
c     locals
c
      include 'bc_common.f'
      integer ic2, ic1, ic0, ict0, ict1, ict2, k
      integer ipivot
c
c     -------------------
c
      if ( bdry_loc .eq. X0Y0Z0 ) then
         ipivot = ifirst0

      else if ( bdry_loc .eq. X1Y0Z0 ) then
         ipivot = ilast0

      else if ( bdry_loc .eq. X0Y1Z0 ) then
         ipivot = ifirst0

      else if ( bdry_loc .eq. X1Y1Z0 ) then
         ipivot = ilast0

      else if ( bdry_loc .eq. X0Y0Z1 ) then
         ipivot = ifirst0

      else if ( bdry_loc .eq. X1Y0Z1 ) then
         ipivot = ilast0

      else if ( bdry_loc .eq. X0Y1Z1 ) then
         ipivot = ifirst0

      else if ( bdry_loc .eq. X1Y1Z1 ) then
         ipivot = ilast0

      endif

      do ic2 = ibeg2,iend2
         ict2 = ic2
         do ic1 = ibeg1, iend1
            ict1 = ic1
            do ic0 = ibeg0, iend0
               ict0 = ipivot
               do k= 1, nd
                  d(ic0,ic1,ic2,k)=d(ict0,ict1,ict2,k)
               enddo
            enddo
         enddo
      enddo
c
c     end of subroutine
c
      return
      end
