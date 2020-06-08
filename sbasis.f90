!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: sbasis.f

      subroutine sbasis (rho,cvi,eint,ilev,jlev,klev,nlev,mode)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs a multiple-arrangement
!     hyperspherical basis set at hyperradius rho.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
!
!     input arrays
!
      dimension cvi(nvi,n)
      dimension eint(n),ilev(n),jlev(n),klev(n),nlev(n)
!
!     local arrays
!
!     dimension h(nvi,nvi),t(nvi),e(nvi)
      allocatable h(:,:),t(:),e(:)

!     dimension c(0:2*nvi),s(0:2*nvi),v(0:2*nvi)
      allocatable c(:),s(:),v(:)

!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

      allocate (h(nvi,nvi),t(nvi),e(nvi))
      allocate (c(0:2*nvi),s(0:2*nvi),v(0:2*nvi))
      allocate (wvi(mvi),xvi(mvi))

!
!     vibrational quadrature rule
!
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
!
!     kinetic energy
!
      pi = acos(-1.d0)
      tscale = pi/tmax
      do jvi = 1,nvi
         t(jvi) = ((jvi*tscale)**2-0.25d0)/rho**2
      enddo
      do jvi = 0,2*nvi
         c(jvi) = 0.d0
         s(jvi) = 0.d0
      enddo
      do kvi = 1,mvi
         theta = xvi(kvi)
         weight = wvi(kvi)/tmax
         oncos2 = weight/cos(theta)**2
         onsin2 = weight/sin(theta)**2
         arg = tscale*theta
         do jvi = 0,2*nvi
            cosine = cos(jvi*arg)
            c(jvi) = c(jvi)+cosine*oncos2
            s(jvi) = s(jvi)+cosine*onsin2
         enddo
      enddo
!
!     loop over arrangements
!
      na = 3-iabs(jpar)
      do ia = 1,na
!
!        reference potential
!
         do jvi = 0,2*nvi
            v(jvi) = 0.d0
         enddo
         do kvi = 1,mvi
            theta = xvi(kvi)
            weight = wvi(kvi)/tmax
            sa = rho*sin(theta)
            call potenl (100.d0,sa,0.d0,va,ia)
            vtheta = weight*va
            arg = tscale*theta
            do jvi = 0,2*nvi
               cosine = cos(jvi*arg)
               v(jvi) = v(jvi)+cosine*vtheta
            enddo
         enddo
!
!        loop over rotational quantum numbers
!
         if (ia .eq. 1) then
            jmin = (1-jpar)/2
            jinc = 1+iabs(jpar)
         else
            jmin = 0
            jinc = 1
         endif
         do j = jmin,jmax,jinc
            kmaxj = min(j,kmax)
            do k = kmin,kmaxj
!
!              vibrational eigenvalue problem
!
               if (mode .eq. 0) then
!
!                 k-dependent vibrational functions
!                 in the interaction region
!
                  ctot = (jtot*(jtot+1)+j*(j+1)-2*k*k)/rho**2
               else
!
!                 k-independent vibrational functions
!                 in the asymptotic region
!
                  ctot = 0.d0
               endif
               cent = j*(j+1)/rho**2
               do jvi = 1,nvi
                  do ivi = 1,jvi
                     h(ivi,jvi) = ctot*(c(jvi-ivi)-c(jvi+ivi)) &
                                + cent*(s(jvi-ivi)-s(jvi+ivi)) &
                                +      (v(jvi-ivi)-v(jvi+ivi))
                     h(jvi,ivi) = h(ivi,jvi)
                  enddo
                  h(jvi,jvi) = h(jvi,jvi)+t(jvi)
               enddo
               call symevp (h,nvi,nvi,e,ierr)
               if (ierr .ne. 0) stop 'sbasis 1'
!
!              specified surface basis functions
!
               do kvi = 1,nvi
                  do 1 i = 1,n
                     if (ilev(i) .ne. ia) go to 1
                     if (jlev(i) .ne. j) go to 1
                     if (klev(i) .ne. k) go to 1
                     if (nlev(i) .ne. kvi-1) go to 1
                     do jvi = 1,nvi
                        cvi(jvi,i) = h(jvi,kvi)
                     enddo
                     eint(i) = e(kvi)
   1              continue
               enddo
            enddo
         enddo
      enddo
      return
      end
