!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: direct.f

      subroutine direct (rho,cvi,cro,ilev,jlev,klev,h,v,ia)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine fills in the direct matrix elements of h
!     within arrangement ia at hyperradius rho. Last modified 9/6/99.
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
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension h(n,n),v(mro,mvi)
!
!     local arrays
!
!     dimension wro(mro),xro(mro)
      allocatable wro(:),xro(:)

!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

!     dimension pro(0:jmax,0:kmax,mro)
      allocatable pro(:,:,:)

!     dimension pvi(mvi,n)
      allocatable pvi(:,:)

!     dimension hro(0:jmax,0:jmax,0:kmax)
      allocatable hro(:,:,:)

      allocate (wro(mro),xro(mro))
      allocate (wvi(mvi),xvi(mvi))
      allocate (pro(0:jmax,0:kmax,mro))
      allocate (pvi(mvi,n))
      allocate (hro(0:jmax,0:jmax,0:kmax))

!
!     arrangement indices
!
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      if (na .lt. 1) return
!
!     quadrature rules
!
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
!
!     rotational basis functions
!
      do kro = 1,mro
         cosa = xro(kro)
         call prot (pro(0,0,kro),jmax,kmax,cosa)
      enddo
!
!     vibrational basis functions
!
      call pvibv (tmin,xvi,tmax,cvi(1,nla),nvi,na,pvi(1,nla),mvi,0)
!
!     loop over quadrature points
!
      do kvi = 1,mvi
         ta = xvi(kvi)
         ra = rho*cos(ta)
         sa = rho*sin(ta)
         call potenl (100.d0,sa,0.d0,va,ia)
!
!        rotational integrals
!
         do k = kmin,kmax
            do jb = k,jmax
               do ja = k,jb
                  hro(ja,jb,k) = 0.d0
               enddo
            enddo
         enddo
         do kro = 1,mro
            cosa = xro(kro)
            call potenl (ra,sa,cosa,vabc,ia)
            v(kro,kvi) = vabc-va
            sab = wro(kro)*wvi(kvi)
            vab = sab*v(kro,kvi)
            nj = jmax+1
            do k = kmin,kmax
               nk = nj-k
               call dsyr ('u',nk,vab,pro(k,k,kro),1,hro(k,k,k),nj)
            enddo
         enddo
         do k = kmin,kmax
            do jb = k,jmax
               do ja = k,jb
                  hro(jb,ja,k) = hro(ja,jb,k)
               enddo
            enddo
         enddo
!
!        vibrational integrals and Coriolis coupling
!
         cent = wvi(kvi)/ra**2
         do lb = nla,nha
            jb = jlev(lb)
            k = klev(lb)
            cmb = cent*cro(1,lb)*pvi(kvi,lb)
            cpb = cent*cro(3,lb)*pvi(kvi,lb)
            do la = nla,lb
               ja = jlev(la)
               ka = klev(la)
               if (ka .eq. k) then
                  h(la,lb) = h(la,lb)+ &
                             pvi(kvi,la)*hro(ja,jb,k)*pvi(kvi,lb)
               else if (ja .eq. jb) then
                  if (ka .eq. k-1) then
                     h(la,lb) = h(la,lb)+pvi(kvi,la)*cmb
                  else if (ka .eq. k+1) then
                     h(la,lb) = h(la,lb)+pvi(kvi,la)*cpb
                  endif
               endif
            enddo
         enddo
      enddo
      do lb = nla,nha
         do la = nla,lb
            h(lb,la) = h(la,lb)
         enddo
      enddo
      return
      end
