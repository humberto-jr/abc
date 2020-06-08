!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: exchng.f

      subroutine exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine fills in the exchange matrix elements of h
!     and s between arrangements ia and ib at hyperradius rho.
!     Last modified 9/6/99.
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
      dimension h(n,n),s(n,n),v(mro,mvi)
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

!     dimension pro(0:jmax,0:kmax+1)
      allocatable pro(:,:)

!     dimension dro(0:kmax+1,0:kmax+1)
      allocatable dro(:,:)

!     dimension prv(n,nblock,0:kmax+1)
      allocatable prv(:,:,:)

!     dimension hrv(nblock,0:jmax,0:kmax+1)
      allocatable hrv(:,:,:)

!     dimension srv(nblock,0:jmax,0:kmax+1)
      allocatable srv(:,:,:)

!     dimension hro(n,0:jmax,0:kmax+1)
      allocatable hro(:,:,:)

!     dimension sro(n,0:jmax,0:kmax+1)
      allocatable sro(:,:,:)

      allocate (wro(mro),xro(mro))
      allocate (wvi(mvi),xvi(mvi))
      allocate (ab2(0:jmax),pvi(n))
      allocate (pro(0:jmax,0:kmax+1))
      allocate (dro(0:kmax+1,0:kmax+1))
      allocate (prv(n,nblock,0:kmax+1))
      allocate (hrv(nblock,0:jmax,0:kmax+1))
      allocate (srv(nblock,0:jmax,0:kmax+1))
      allocate (hro(n,0:jmax,0:kmax+1))
      allocate (sro(n,0:jmax,0:kmax+1))

!
!     arrangement indices
!
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      call arrang (ilev,n,jpar,ib,nlb,nhb,nb)
      if (na.lt.1 .or. nb.lt.1) return
!
!     projection ranges (allowing for Coriolis terms)
!
      kamax = kmax
      if (kmax .lt. min(jtot,jmax)) kamax = kmax+1
      kbmax = kmax
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
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
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
            if (tb .lt. tmax) then
               inner = inner+1
               if (inner .eq. 1) then
                  do ka = kmin,kamax
                     do ja = ka,jmax
                        do lb = nlb,nhb
                           hro(lb,ja,ka) = 0.d0
                           sro(lb,ja,ka) = 0.d0
                        enddo
                     enddo
                  enddo
               endif
               iblock = iblock+1
!
!              jacobian factor and interaction potential
!
               weight = wro(kro)*wvi(kvi)
               sab = weight*sin(2.d0*ta)/sin(2.d0*tb)
               vab = sab*v(kro,kvi)
!
!              parity-adapted reduced rotation matrix elements
!
               do kb = kmin,kbmax
                  do ka = kmin,kamax
                     dro(ka,kb) = rotmel(jtot,ipar,ka,kb,xab)
                  enddo
               enddo
!
!              rovibrational functions in arrangement ib
!
               call prot (pro,jmax,kbmax,cosb)
               call pvib (tmin,tb,tmax,cvi(1,nlb),nvi,nb,pvi(nlb),0)
               do lb = nlb,nhb
                  jb = jlev(lb)
                  kb = klev(lb)
                  pb = ab2(jb)*pro(jb,kb)*pvi(lb)
                  do ka = kmin,kamax
                     prv(lb,iblock,ka) = dro(ka,kb)*pb
                  enddo
               enddo
!
!              rotational functions in arrangement ia
!
               call prot (pro,jmax,kamax,cosa)
               do ka = kmin,kamax
                  do ja = ka,jmax
                     hrv(iblock,ja,ka) = vab*pro(ja,ka)
                     srv(iblock,ja,ka) = sab*pro(ja,ka)
                  enddo
               enddo
!
!              partial potential and overlap matrix elements
!
               if (iblock .eq. nblock) then
                  do ka = kmin,kamax
                     nk = jmax-ka+1
                     call dgemm ('n','n',nb,nk,nblock,1.d0, &
                                 prv(nlb,1,ka),n,hrv(1,ka,ka),nblock, &
                                 1.d0,hro(nlb,ka,ka),n)
                     call dgemm ('n','n',nb,nk,nblock,1.d0, &
                                 prv(nlb,1,ka),n,srv(1,ka,ka),nblock, &
                                 1.d0,sro(nlb,ka,ka),n)
                  enddo
                  iblock = 0
               endif
            endif
         enddo
         if (iblock .gt. 0) then
            do ka = kmin,kamax
               nk = jmax-ka+1
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n, &
                           hrv(1,ka,ka),nblock,1.d0,hro(nlb,ka,ka),n)
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n, &
                           srv(1,ka,ka),nblock,1.d0,sro(nlb,ka,ka),n)
            enddo
         endif
         if (inner .gt. 0) then
!
!           vibrational functions in arrangement ia
!
            call pvib (tmin,ta,tmax,cvi(1,nla),nvi,na,pvi(nla),0)
!
!           full potential, overlap, and Coriolis matrix elements
!
            cent = 1.d0/(rho*cos(ta))**2
            do la = nla,nha
               ja = jlev(la)
               ka = klev(la)
               do lb = nlb,nhb
                  h(la,lb) = h(la,lb)+pvi(la)*hro(lb,ja,ka)
                  s(la,lb) = s(la,lb)+pvi(la)*sro(lb,ja,ka)
               enddo
               if (ka .gt. kmin) then
                  crv = cro(1,la)*cent*pvi(la)
                  do lb = nlb,nhb
                     h(la,lb) = h(la,lb)+crv*sro(lb,ja,ka-1)
                  enddo
               endif
               if (ka.lt.ja .and. ka.lt.kamax) then
                  crv = cro(3,la)*cent*pvi(la)
                  do lb = nlb,nhb
                     h(la,lb) = h(la,lb)+crv*sro(lb,ja,ka+1)
                  enddo
               endif
            enddo
         endif
      enddo
      return
      end
