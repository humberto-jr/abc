!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: bessel.f

      subroutine bessel (x,y,z,cvr,cvi,eint, &
                         ilev,jlev,klev,llev,nlev,ered)
      implicit double precision (a-h,o-z)
      double precision llev
!
!     -----------------------------------------------------------------
!     This subroutine forms the asymptotic solution and derivative
!     matrices a,b,c, and d needed for matching the asymptotic Delves
!     coordinate log derivative matrix y onto a reactance matrix k,
!     as in eqs. (116) to (121) of Pack and Parker, and then proceeds
!     to calculate the reactance matrix in the array x.
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
      dimension x(n,n),y(n,n),z(n,n)
      dimension cvr(nvi,n),cvi(nvi,n),eint(n)
      dimension ilev(n),jlev(n),klev(n),llev(n),nlev(n)
!
!     local arrays
!
!     dimension wvi(mvi),xvi(mvi),pvi(nvi,2)
      allocatable wvi(:),xvi(:),pvi(:,:)

!     dimension fr(n),fvi(n),gvi(n)
      allocatable fr(:),fvi(:),gvi(:)

!     dimension fa(n),fb(n),fc(n),fd(n)
      allocatable fa(:),fb(:),fc(:),fd(:)

!     dimension a(n,0:nvi-1),b(n,0:nvi-1)
      allocatable a(:,:),b(:,:)

!     dimension erow(n),ecol(n),srow(n),scol(n)
      allocatable erow(:),ecol(:),srow(:),scol(:)

      allocate (wvi(mvi),xvi(mvi),pvi(nvi,2))
      allocate (fr(n),fvi(n),gvi(n))
      allocate (fa(n),fb(n),fc(n),fd(n))
      allocate (a(n,0:nvi-1),b(n,0:nvi-1))
      allocate (erow(n),ecol(n),srow(n),scol(n))

!
!     initialisation
!
      do nj = 0,nvi-1
         do i = 1,n
            a(i,nj) = 0.d0
            b(i,nj) = 0.d0
         enddo
      enddo
      do j = 1,n
         do i = 1,n
            x(i,j) = 0.d0
            z(i,j) = 0.d0
         enddo
      enddo
!
!     vibrational quadrature rule
!
      pi = acos(-1.d0)
      piby2 = 0.5d0*pi
      rtrho = sqrt(rmax)
      tworho = 2.d0*rmax
      smin = 0.d0
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rmax))
      call qvib (tmin,tmax,mvi,wvi,xvi)
!
!     Delves vibrational quadrature
!
      do kvi = 1,mvi
         weight = rtrho*wvi(kvi)
         ta = xvi(kvi)
         cta = cos(ta)
         sta = sin(ta)
         ra = rmax*cta
         sa = rmax*sta
         if (sa .lt. smax) then
!
!           Delves vibrational functions
!
            call pvib (tmin,ta,tmax,cvr,nvi,n,fr,0)
            do i = 1,n
               fr(i) = weight*fr(i)
            enddo
!
!           Jacobi vibrational functions
!
            call pvib (smin,sa,smax,cvi,nvi,n,fvi,0)
            call pvib (smin,sa,smax,cvi,nvi,n,gvi,1)
!
!           Jacobi translational (Riccati-Bessel) functions
!
            do i = 1,n
               ell = llev(i)
               psq = ered-eint(i)
               if (psq .gt. 0.d0) then
                  p = sqrt(psq)
                  arg = p*ra
                  call rbesjy (ell,arg,cj,dj,ej,cy,dy,ey)
!
!                 exponential scaling of j and y
!                 (for stability near channel thresholds)
!
                  if (kvi .eq. 1) then
                     ecol(i) = ej
                     scol(i) = exp(ej)
                     erow(i) = ey
                     srow(i) = exp(-ey)
                  else
                     sj = exp(ej-ecol(i))
                     cj = sj*cj
                     dj = sj*dj
                     sy = exp(ey-erow(i))
                     cy = sy*cy
                     dy = sy*dy
                  endif
                  rtp = sqrt(p)
                  atr = cj/rtp
                  btr = cy/rtp
                  ctr = dj*rtp
                  dtr = dy*rtp
               else
                  p = sqrt(-psq)
                  arg = p*ra
                  call rbessk (ell,arg,ck,dk,ek)
!
!                 exponential scaling of k
!                 (for the same reason)
!
                  if (kvi .eq. 1) then
                     ecol(i) = 0.d0
                     scol(i) = 0.d0
                     erow(i) = ek
                     srow(i) = 0.d0
                  else
                     sk = exp(ek-erow(i))
                     ck = sk*ck
                     dk = sk*dk
                  endif
                  rtp = sqrt(p)
                  atr = 0.d0
                  btr = ck/rtp
                  ctr = 0.d0
                  dtr = dk*rtp
               endif
!
!              Pack and Parker eqs. (116) to (121)
!
               fa(i) = atr*fvi(i)
               fb(i) = btr*fvi(i)
               fc(i) = cta*ctr*fvi(i)+sta*atr*gvi(i)
               fd(i) = cta*dtr*fvi(i)+sta*btr*gvi(i)
               fc(i) = fa(i)/tworho+fc(i)
               fd(i) = fb(i)/tworho+fd(i)
            enddo
!
!           integral accumulation
!
            do j = 1,n
               ij = ilev(j)
               jj = jlev(j)
               kj = klev(j)
               nj = nlev(j)
               do i = 1,n
                  ii = ilev(i)
                  ji = jlev(i)
                  ki = klev(i)
                  if (ii.eq.ij .and. ji.eq.jj .and. ki.eq.kj) then
                     a(i,nj) = a(i,nj)+fr(i)*fa(j)
                     b(i,nj) = b(i,nj)+fr(i)*fb(j)
                     x(i,j)  = x(i,j) -fr(i)*fc(j)
                     z(i,j)  = z(i,j) -fr(i)*fd(j)
                  endif
               enddo
            enddo
         endif
      enddo
!
!     y to k
!
      do j = 1,n
         ij = ilev(j)
         jj = jlev(j)
         kj = klev(j)
         nj = nlev(j)
         do k = 1,n
            ik = ilev(k)
            jk = jlev(k)
            kk = klev(k)
            if (ik.eq.ij .and. jk.eq.jj .and. kk.eq.kj) then
               do i = 1,n
                  x(i,j) = x(i,j)+y(i,k)*a(k,nj)
                  z(i,j) = z(i,j)+y(i,k)*b(k,nj)
               enddo
            endif
         enddo
      enddo
      call gensol (z,n,n,x,n,n,ierr)
      if (ierr .ne. 0) stop 'bessel 1'
!
!     elimination of exponential scaling factors
!     from the reactance matrix
!
      do j = 1,n
         do i = 1,n
            x(i,j) = srow(i)*x(i,j)*scol(j)
         enddo
      enddo
      return
      end
