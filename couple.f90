!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: couple.f

      subroutine couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs the Hamiltonian and overlap matrices
!     h and s that are needed to solve the hyperspherical coordinate
!     surface eigenvalue problem at the given value of rho.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension eint(n),ilev(n),jlev(n),klev(n)
      dimension h(n,n),s(n,n)
!
!     local arrays
!
!     dimension v(mro,mvi)
      allocatable v(:,:)

      allocate (v(mro,mvi))

!
      do j = 1,n
         do i = 1,n
            h(i,j) = 0.d0
            s(i,j) = 0.d0
         enddo
         s(j,j) = 1.d0
      enddo
      na = 3-iabs(jpar)
      do ia = 1,na
         call direct (rho,cvi,cro,ilev,jlev,klev,h,v,ia)
         do ib = 1,na
            if (ib .ne. ia) then
               call exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
            endif
         enddo
      enddo
      if (na .eq. 2) then
         ia = 2
         ib = 3
         call exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
      endif
      do j = 1,n
         do i = 1,n
            h(i,j) = h(i,j)+eint(i)*s(i,j)
         enddo
      enddo
      heps = 0.d0
      seps = 0.d0
      do j = 2,n
         do i = 1,j-1
            heps = max(heps,abs(h(i,j)-h(j,i)))
            h(i,j) = 0.5d0*(h(i,j)+h(j,i))
            h(j,i) = h(i,j)
            seps = max(seps,abs(s(i,j)-s(j,i)))
            s(i,j) = 0.5d0*(s(i,j)+s(j,i))
            s(j,i) = s(i,j)
         enddo
      enddo
!
!     couple summary
!
      write (6,61) rho,heps/rmlmda,seps
  61  format(/1x,'COUPLE:', &
      '   rho = ',f8.3, &
      '   heps = ',f8.3, &
      '  seps = ',1p,e8.1)
      return
      end
