c     lin.f
c     ----------------------------------------------------------------- 
c     Linear algebra routines for the CCP6 reactive scattering program.
c     This version with LAPACK calls dated 31 March 2000.
c     ----------------------------------------------------------------- 
c
      subroutine gensol (a,lda,n,b,ldb,m,ierr)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses LAPACK DGETRF and DGETRS to 
c     solve the linear equations A*X = B.
c     ----------------------------------------------------------------- 
c
      dimension a(lda,n),b(ldb,m)
      dimension ipiv(n)
c
      call dgetrf (n,n,a,lda,ipiv,ierr)
      if (ierr .ne. 0) return
      call dgetrs ('N',n,m,a,lda,ipiv,b,ldb,ierr)
      return
      end 
 
      subroutine syminv (a,lda,n,ierr)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses LAPACK DSYTRF and DSYTRI 
c     to invert a real symmetric indefinite matrix. 
c     ----------------------------------------------------------------- 
c
      dimension a(lda,n)
      dimension ipiv(n),work(32*n)
c
      lwork = 32*n
      call dsytrf ('U',n,a,lda,ipiv,work,lwork,ierr)
      if (ierr .ne. 0) return
      call dsytri ('U',n,a,lda,ipiv,work,ierr)
      do j = 2,n
         do i = 1,j-1
            a(j,i) = a(i,j)
         enddo
      enddo
      return
      end

      subroutine symevp (a,lda,n,d,ierr)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses LAPACK DSYEV to
c     diagonalise a real symmetric matrix.
c     ----------------------------------------------------------------- 
c
      dimension a(lda,n),d(n)
      dimension work(34*n)
c
      lwork = 34*n
      call dsyev ('V','L',n,a,lda,d,work,lwork,ierr)
      return
      end

      subroutine rstevp (d,e,n,v,ldv,ierr)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses LAPACK DSTEV to diagonalise
c     a real symmetric tridiagonal matrix.
c     ----------------------------------------------------------------- 
c
      dimension d(n),e(n),v(ldv,n)
      dimension work(2*n)
c
      do j = 1,n-1
         e(j) = e(j+1)
      enddo
      e(n) = 0.d0
      call dstev ('V',n,d,e,v,ldv,work,ierr)
      return
      end
