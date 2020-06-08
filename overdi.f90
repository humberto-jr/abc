!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: overdi.f

      subroutine overdi (rhoa,cva,rhob,cvb,ilev,jlev,klev,s,ia)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine fills in the matrix elements of the setor-
!     to-sector overlap matrix s within arrangement ia.
!     Compare with subroutine direct.
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
      dimension cva(nvi,n),cvb(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension s(n,n)
!
!     local arrays
!
!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

!     dimension pva(mvi,n),pvb(mvi,n)
      allocatable pva(:,:),pvb(:,:)

      allocate (wvi(mvi),xvi(mvi))
      allocate (pva(mvi,n),pvb(mvi,n))

!
!     arrangement indices
!
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      if (na .lt. 1) return
!
!     vibrational quadrature rule
!
      tmin = 0.d0
      tamax = asin(min(1.d0,smax/rhoa))
      tbmax = asin(min(1.d0,smax/rhob))
      call qvib (tmin,tamax,mvi,wvi,xvi)
!
!     vibrational basis functions
!
      call pvibv (tmin,xvi,tamax,cva(1,nla),nvi,na,pva(1,nla),mvi,0)
      call pvibv (tmin,xvi,tbmax,cvb(1,nla),nvi,na,pvb(1,nla),mvi,0)
!
!     vibrational integrals
!
      do kvi = 1,mvi
         tb = xvi(kvi)
         sab = wvi(kvi)
         if (tb .lt. tbmax) then
            do j = nla,nha
               svb = sab*pvb(kvi,j)
               do i = nla,nha
                  if (jlev(i).eq.jlev(j) .and. klev(i).eq.klev(j)) then
                     s(i,j) = s(i,j)+pva(kvi,i)*svb
                  endif
               enddo
            enddo
         endif
      enddo
      return
      end
