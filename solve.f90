!     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: solve.f

      subroutine solve (cvi,cro,ilev,jlev,klev,nlev,ered,nnrg)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine solves the hyperspherical coordinate
!     coupled-channel equations between rho = rmin and rmax.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /inputs/ emax,mtr
      common /ranges/ rmin,rmax,smax
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension ered(nnrg)
!
!     local arrays
!
!     dimension cvp(nvi,n)
      allocatable cvp(:,:)

!     dimension eint(n),e(n)
      allocatable eint(:),e(:)

!     dimension h(n,n),s(n,n),c(n,n)
      allocatable h(:,:),s(:,:),c(:,:)

      allocate (cvp(nvi,n))
      allocate (eint(n),e(n))
      allocate (h(n,n),s(n,n),c(n,n))

!
!     loop over sectors
!
      npp = 0
      drho = (rmax-rmin)/mtr
      rhop = rmin-0.5d0*drho
      do ktr = 1,mtr
!
!        surface eigenvalue problem
!
         rho = rhop+drho
         call sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
!
!        sector-to-sector transformation matrix
!
         if (npp .eq. 0) then
            call putrec (c,n,0)
         else
            call metric (rhop,cvp,npp,rho,cvi,np,ilev,jlev,klev,h,s,c)
         endif
!
!        log derivative propagation
!
         call logder (drho,ered,nnrg,h,s,c,e,n,np,npp)
!
!        setup for the next sector
!
         do j = 1,n
            do k = 1,nvi
               cvp(k,j) = cvi(k,j)
            enddo
         enddo
         npp = np
         rhop = rho
      enddo
!
!     transformation to final coupled-channel basis set at rmax
!
      npp = n+1
      rho = rmax
      call sbasis (rho,cvi,e,ilev,jlev,klev,nlev,1)
      call overlp (rhop,cvp,rho,cvi,ilev,jlev,klev,h)
!     call overlp (rho,cvi,rhop,cvp,ilev,jlev,klev,s)
!     do j = 1,n
!     do i = 1,n
!     h(i,j) = 0.5d0*(h(i,j)+s(j,i))
!     enddo
!     enddo
      call getrec (s,n,0)
      call dgemm  ('t','n',n,np,n,1.d0,h,n,s,n,0.d0,c,n)
      call logder (drho,ered,nnrg,h,s,c,e,n,np,npp)
      return
      end
