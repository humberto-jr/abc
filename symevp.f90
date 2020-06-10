!     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: symevp.f

      subroutine symevp (a,lda,n,d,ierr)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses LAPACK DSYEV to
!     diagonalise a real symmetric matrix.
!     -----------------------------------------------------------------
!
      dimension a(lda,n),d(n)

!     dimension work(34*n)
      allocatable work(:)

      allocate (work(34*n))
!
      lwork = 34*n
      call dsyev ('v','l',n,a,lda,d,work,lwork,ierr)
      return
      end
