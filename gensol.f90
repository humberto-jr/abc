!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: gensol.f

      subroutine gensol (a,lda,n,b,ldb,m,ierr)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses LAPACK DGETRF and DGETRS to
!     solve the linear equations A*X = B.
!     -----------------------------------------------------------------
!
      dimension a(lda,n),b(ldb,m)

!     dimension ipiv(n)
      allocatable ipiv(:)

      allocate (ipiv(n))

!
      call dgetrf (n,n,a,lda,ipiv,ierr)
      if (ierr .ne. 0) return
      call dgetrs ('N',n,m,a,lda,ipiv,b,ldb,ierr)
      return
      end
