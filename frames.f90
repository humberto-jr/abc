!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: frames.f

      subroutine frames (x,c,cro,jlev,klev,llev,mode)
      implicit double precision (a-h,o-z)
      double precision llev
!
!     -----------------------------------------------------------------
!     This subroutine transforms the symmetric matrix x from a
!     (possibly truncated) BF basis set to a (possibly truncated)
!     SF basis set if mode = +1, or vice versa if mode = -1.
!     The (possibly fractional) orbital angular momentum quantum
!     numbers of the SF basis are returned in llev.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
!
!     input arrays
!
      dimension x(n,n),c(n,n)
      dimension cro(3,n)
      dimension jlev(n),klev(n),llev(n)
!
!     local arrays
!
!     dimension d(n),e(n)
      allocatable d(:),e(:)

      allocate (d(n),e(n))

!
!     loop over blocks of channels with the same a,v, and j,
!     assuming an a,v,j,k (or a,v,j,l) channel ordering
!
      nlo = 1
   1  continue
         j = jlev(nlo)
         kmaxj = min(j,jtot,kmax)
         ndo = kmaxj-kmin+1
         nhi = nlo+ndo-1
         if (nhi .gt. n) stop 'frames 1'
!
!        form the BF to SF transformation matrix by diagonalising
!        the operator L^2 stored in cro, making sure to apply a
!        consistent phase convention
!
         do i = nlo,nhi
            d(i) = cro(2,i)
            e(i) = cro(1,i)
         enddo
         call rstevp (d(nlo),e(nlo),ndo,c(nlo,nlo),n,ierr)
         if (ierr .ne. 0) stop 'frames 2'
         phase = (-1)**(j+kmaxj)
         do k = nlo,nhi
            if (phase*c(nhi,k) .lt. 0.d0) then
               do i = nlo,nhi
                  c(i,k) = -c(i,k)
               enddo
            endif
            llev(k) = sqrt(d(k)+0.25d0)-0.5d0
            if (kmaxj .eq. min(j,jtot)) llev(k) = nint(llev(k))
         enddo
!
!        transform the present block of the matrix x, applying
!        the inverse transformation if mode = -1
!
         if (ndo .gt. 1) then
            if (mode .lt. 0) then
               do j = nlo+1,nhi
                  do i = nlo,j-1
                     cswap = c(i,j)
                     c(i,j) = c(j,i)
                     c(j,i) = cswap
                  enddo
               enddo
            endif
            do i = 1,n
               do j = nlo,nhi
                  d(j) = 0.d0
                  do k = nlo,nhi
                     d(j) = d(j)+x(i,k)*c(k,j)
                  enddo
               enddo
               do j = nlo,nhi
                  x(i,j) = d(j)
               enddo
            enddo
            do j = 1,n
               do i = nlo,nhi
                  d(i) = 0.d0
                  do k = nlo,nhi
                     d(i) = d(i)+c(k,i)*x(k,j)
                  enddo
               enddo
               do i = nlo,nhi
                  x(i,j) = d(i)
               enddo
            enddo
         endif
!
!        and increment nlo for the next block
!
         nlo = nhi+1
      if (nlo .le. n) go to 1
      return
      end
