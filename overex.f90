!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: overex.f

      subroutine overex (rhoa,cva,rhob,cvb,ilev,jlev,klev,s,ia,ib)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine fills in the matrix elements of the sector-
!     to-sector overlap matrix s between arrangements ia and ib.
!     Compare with subroutine exchng.
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
      parameter (nblock = 32)

!     dimension wro(mro),xro(mro)
      allocatable wro(:),xro(:)

!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

!     dimension ab2(0:jmax),pvi(n)
      allocatable ab2(:),pvi(:)

!     dimension pro(0:jmax,0:kmax)
      allocatable pro(:,:)

!     dimension dro(0:kmax,0:kmax)
      allocatable dro(:,:)

!     dimension prv(n,nblock,0:kmax)
      allocatable prv(:,:,:)

!     dimension srv(nblock,0:jmax,0:kmax)
      allocatable srv(:,:,:)

!     dimension sro(n,0:jmax,0:kmax)
      allocatable sro(:,:,:)

      allocate (wro(mro),xro(mro))
      allocate (wvi(mvi),xvi(mvi))
      allocate (ab2(0:jmax),pvi(n))
      allocate (pro(0:jmax,0:kmax))
      allocate (dro(0:kmax,0:kmax))
      allocate (prv(n,nblock,0:kmax))
      allocate (srv(nblock,0:jmax,0:kmax))
      allocate (sro(n,0:jmax,0:kmax))

!
!     arrangement indices
!
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      call arrang (ilev,n,jpar,ib,nlb,nhb,nb)
      if (na.lt.1 .or. nb.lt.1) return
!
!     A+B2 permutation symmetry
!
      do j = 0,jmax
         ab2(j) = 1.d0
         if (jpar .ne. 0) then
            if (ia+ib .eq. 3) ab2(j) = sqrt(2.d0)
            if (ia+ib .eq. 5) ab2(j) = jpar*(-1)**j
         endif
      enddo
!
!     quadrature rules
!
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tamax = asin(min(1.d0,smax/rhoa))
      tbmax = asin(min(1.d0,smax/rhob))
      call qvib (tmin,tamax,mvi,wvi,xvi)
!
!     loop over quadrature points
!
      do kvi = 1,mvi
         inner = 0
         iblock = 0
         do kro = 1,mro
            ta = xvi(kvi)
            cosa = xro(kro)
!
!           coordinate transformation
!
            call coords (ta,cosa,ia,tb,cosb,ib,xab)
            if (tb .lt. tbmax) then
               inner = inner+1
               if (inner .eq. 1) then
                  do ka = kmin,kmax
                     do ja = ka,jmax
                        do lb = nlb,nhb
                           sro(lb,ja,ka) = 0.d0
                        enddo
                     enddo
                  enddo
               endif
               iblock = iblock+1
!
!              jacobian factor
!
               weight = wro(kro)*wvi(kvi)
               sab = weight*sin(2.d0*ta)/sin(2.d0*tb)
!
!              parity-adapted reduced rotation matrix elements
!
               do kb = kmin,kmax
                  do ka = kmin,kmax
                     dro(ka,kb) = rotmel(jtot,ipar,ka,kb,xab)
                  enddo
               enddo
!
!              rovibrational functions in arrangement ib
!
               call prot (pro,jmax,kmax,cosb)
               call pvib (tmin,tb,tbmax,cvb(1,nlb),nvi,nb,pvi(nlb),0)
               do lb = nlb,nhb
                  jb = jlev(lb)
                  kb = klev(lb)
                  pb = ab2(jb)*pro(jb,kb)*pvi(lb)
                  do ka = kmin,kmax
                     prv(lb,iblock,ka) = dro(ka,kb)*pb
                  enddo
               enddo
!
!              rotational functions in arrangement ia
!
               call prot (pro,jmax,kmax,cosa)
               do ka = kmin,kmax
                  do ja = ka,jmax
                     srv(iblock,ja,ka) = sab*pro(ja,ka)
                  enddo
               enddo
!
!              partial overlap matrix elements
!
               if (iblock .eq. nblock) then
                  do ka = kmin,kmax
                     nk = jmax-ka+1
                     call dgemm ('n','n',nb,nk,nblock,1.d0, &
                                 prv(nlb,1,ka),n,srv(1,ka,ka),nblock, &
                                 1.d0,sro(nlb,ka,ka),n)
                  enddo
                  iblock = 0
               endif
            endif
         enddo
         if (iblock .gt. 0) then
            do ka = kmin,kmax
               nk = jmax-ka+1
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n, &
                           srv(1,ka,ka),nblock,1.d0,sro(nlb,ka,ka),n)
            enddo
         endif
         if (inner .gt. 0) then
!
!           vibrational functions in arrangement ia
!
            call pvib (tmin,ta,tamax,cva(1,nla),nvi,na,pvi(nla),0)
!
!           full overlap matrix elements
!
            do la = nla,nha
               ja = jlev(la)
               ka = klev(la)
               do lb = nlb,nhb
                  s(la,lb) = s(la,lb)+pvi(la)*sro(lb,ja,ka)
               enddo
            enddo
         endif
      enddo
      return
      end
