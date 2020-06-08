!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rstevp.f

      subroutine rstevp (d,e,n,v,ldv,ierr)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses LAPACK DSTEV to diagonalise
!     a real symmetric tridiagonal matrix.
!     -----------------------------------------------------------------
!
      dimension d(n),e(n),v(ldv,n)

!     dimension work(2*n)
      allocatable work(:)

      allocate (work(2*n))

!
      do j = 1,n-1
         e(j) = e(j+1)
      enddo
      e(n) = 0.d0
      call dstev ('V',n,d,e,v,ldv,work,ierr)
      return
      end
