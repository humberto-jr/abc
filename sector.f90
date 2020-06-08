!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: sector.f

      subroutine sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine sets up and solves the hyperspherical
!     coordinate surface eigenvalue problem at hyperradius rho.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /scales/ scale(3),rmlmda
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension h(n,n),s(n,n),c(n,n),e(n)
!
!     local arrays
!
!     dimension eint(n),sigma(n)
      allocatable eint(:),sigma(:)

      allocate (eint(n),sigma(n))

!
!     surface basis set
!
      call sbasis (rho,cvi,eint,ilev,jlev,klev,nlev,0)
!
!     hamiltonian and overlap matrices
!
      call couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
!
!     surface eigenvalue problem
!
      call hcsevp (h,s,c,n,e,sigma,np)
!
!     sector summary
!
      write (6,61) e(1)/rmlmda,e(np)/rmlmda,np

  61  format(1x,'SECTOR:', &
      '  e(1) = ',f8.3, &
      '  e(np) = ',f8.3, &
      '    np = ',i8)
      return
      end
