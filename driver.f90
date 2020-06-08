!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: driver.f

      subroutine driver (nmax)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Main driving subroutine for hyperspherical coordinate calculation
!     -----------------------------------------------------------------
!
!     common blocks
!
      double precision mass,mtot,mred
      common /arrays/ mro,mvi,nvi,n
      common /energy/ enrg,dnrg,nnrg
      common /inputs/ emax,mtr
      common /masses/ mass(3),mtot,mred
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
!
!     global arrays
!
!     parameter (nmax = 1500)

!     dimension ered(nmax)
      allocatable ered(:)

!     dimension cvi(nvi,nmax),cro(3,nmax),eint(nmax)
      allocatable cvi(:,:),cro(:,:),eint(:)

!     dimension ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax)
      allocatable ilev(:),jlev(:),klev(:),nlev(:)

!     dimension cvr(nvi,nmax)
      allocatable cvr(:,:)

      allocate (ered(nmax))
      allocate (cvi(nvi,nmax),cro(3,nmax),eint(nmax))
      allocate (ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax))
      allocate (cvr(nvi,nmax))

!
!     rovibrational basis set
!
      call basis (cvi,cro,eint,ilev,jlev,klev,nlev,nmax)
!
!     scattering energies
!
      eminr = 1.d+30
      eminp = 1.d+30
      do j = 1,n
         if (ilev(j) .eq. 1) then
            eminr = min(eminr,eint(j))
         else
            eminp = min(eminp,eint(j))
         endif
      enddo
      eminr = eminr/rmlmda
      eminp = eminp/rmlmda
      emin = max(eminr,eminp)
      nrg = 0
      do inrg = 1,nnrg
         if (enrg .gt. emin) then
            nrg = nrg+1
            if (nrg .le. nmax) then
               ered(nrg) = rmlmda*enrg
            endif
         else
            write (6,61) enrg
  61        format(1x,'*** No open reaction channels at ' &
            ,f7.4,' eV. Energy skipped.')
         endif
         enrg = enrg+dnrg
      enddo
      nnrg = nrg
      if (nnrg .gt. nmax) then
         write (6,99) nnrg,nmax
  99     format(/1x,'Error: nnrg = ',i6,' > nmax = ',i6/1x, &
                    'Run halted in subroutine driver')
         stop
      endif
      if (nnrg .eq. 0) return
!
!     coupled channel propagation
!
      call setrec (n,nnrg)
      call solve (cvr,cro,ilev,jlev,klev,nlev,ered,nnrg)
!
!     asymptotic matching and output
!
      call match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
      call endrec
      return
      end
