!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: hcsevp.f

      subroutine hcsevp (h,s,c,n,e,sigma,np)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine solves the surface eigenvalue problem
!     (H-E*S)*C=0 by canonical orthogonalisation of the basis.
!     -----------------------------------------------------------------
!
      parameter (CUTOFF=1.d-04)
      dimension h(n,n),s(n,n),c(n,n)
      dimension e(n),sigma(n)
!
      stest = 0.d0
      do j = 1,n
         stestj = 0.d0
         do i = 1,n
            stestj = stestj+abs(s(i,j))
         enddo
         stestj = stestj+abs(1.d0-s(j,j))-abs(s(j,j))
         stest = max(stest,stestj)
      enddo
      if (stest .lt. 1.d-08) then
         np = n
         do j = 1,n
            do i = 1,n
               c(i,j) = h(i,j)
            enddo
         enddo
         call symevp (c,n,n,e,ierr)
         if (ierr .ne. 0) stop 'hcsevp 1'
      else
         do j = 1,n
            do i = 1,n
               s(i,j) = -s(i,j)
            enddo
         enddo
         call symevp (s,n,n,sigma,ierr)
         if (ierr .ne. 0) stop 'hcsevp 2'
         np = 0
         do j = 1,n
            sigma(j) = -sigma(j)
            if (sigma(j) .gt. CUTOFF) then
               scale = 1.d0/sqrt(sigma(j))
               do i = 1,n
                  s(i,j) = scale*s(i,j)
               enddo
               np = j
            endif
         enddo
         if (np .gt. 0) then
            call dgemm ('n','n',n,np,n,1.d0,h,n,s,n,0.d0,c,n)
            call dgemm ('t','n',np,np,n,1.d0,s,n,c,n,0.d0,h,n)
            call symevp (h,n,np,e,ierr)
            if (ierr .ne. 0) stop 'hcsevp 3'
            call dgemm ('n','n',n,np,np,1.0d0,s,n,h,n,0.0d0,c,n)
         endif
      endif
      return
      end
