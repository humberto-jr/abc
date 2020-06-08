      program abc
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     CCP6 hyperspherical coordinate reactive scattering program ABC.
c     This version dated 31 March 2000.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      double precision mass,mtot,mred
      common /arrays/ mro,mvi,nvi,n
      common /energy/ enrg,dnrg,nnrg
      common /inputs/ emax,mtr
      common /masses/ mass(3),mtot,mred
      common /output/ nout,jout
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     input parameters
c
      namelist /input/ mass,jtot,ipar,jpar,jmax,kmax,rmax,mtr,emax,
     +                 enrg,dnrg,nnrg,nout,jout
      read (5,input)
c
c     universal constants
c     (energies in eV, lengths in bohr, masses in amu)
c
      do ia = 1,3
         imass = int(mass(ia))
         if (imass .eq.  1) then
            mass(ia) = 1.007825d0
         else if (imass .eq.  2) then
            mass(ia) = 2.014000d0
         else if (imass .eq. 19) then
            mass(ia) = 18.99840d0
         else if (imass .eq. 35) then
            mass(ia) = 34.96885d0
         else if (imass .eq. 37) then
            mass(ia) = 36.96590d0
         else
            stop 'abc 1'
         endif
      enddo
      hbarsq = 0.014927625d0
c
c     scale factors for MSJ coordinates
c
      mtot = mass(1)+mass(2)+mass(3)
      mred = sqrt(mass(1)*mass(2)*mass(3)/mtot)
      scale(1) = sqrt((mass(1)/mred)*(1.d0-mass(1)/mtot))
      scale(2) = sqrt((mass(2)/mred)*(1.d0-mass(2)/mtot))
      scale(3) = sqrt((mass(3)/mred)*(1.d0-mass(3)/mtot))
      rmlmda = 2.0d0*mred/hbarsq
c
c     angular momentum projections
c
      if ((-1)**jtot .eq. ipar) then
         kmin = 0
      else
         kmin = 1
      endif
      kmax = min(jtot,jmax,kmax)
      if (mass(2) .ne. mass(3)) jpar = 0
      nout = max(nout,0)
      if (jpar .eq. -1) then
         jout = max(jout,1)
      else
         jout = max(jout,0)
      endif
c
c     input summary
c
      write (6,61) mass,jtot,ipar,jpar,jmax,kmin,kmax,rmax,mtr,emax,
     +             enrg,dnrg,nnrg,nout,jout
  61  format(/1x,'INPUT:'/1x,70('-')/1x,
     +'  mass = ',3f16.3/1x,
     +'  jtot = ',i32/1x,
     +'  ipar = ',i32/1x,
     +'  jpar = ',i32/1x,
     +'  jmax = ',i32/1x,
     +'  kmin = ',i32/1x,
     +'  kmax = ',i32/1x,
     +'  rmax = ',f32.3/1x,
     +'   mtr = ',i32/1x,
     +'  emax = ',f32.3/1x,
     +'  enrg = ',f32.3/1x,
     +'  dnrg = ',f32.3/1x,
     +'  nnrg = ',i32/1x,
     +'  nout = ',i32/1x,
     +'  jout = ',i32/1x,70('-'))
c
c     calculation
c
      if (kmin .gt. kmax) stop 'abc 2'
      call parset
      call driver 
      stop
      end

      subroutine parset 
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     This subroutine uses the input parameters emax and jmax 
c     to determine rmin, smax, mvi, mro, and nvi.
c     -----------------------------------------------------------------   
c
c     common blocks
c
      double precision mass,mtot,mred
      common /arrays/ mro,mvi,nvi,n
      common /masses/ mass(3),mtot,mred
      common /inputs/ emax,mtr
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     local arrays
c
      parameter (m = 1000)
      dimension v(m)
      dimension smin(3),vmid(3)
c
c     smax
c
      ds = 0.01d0
      smax = 0.d0
      scut = 5.d0
      ered = rmlmda*emax
      do ia = 1,3
         lmid = 1
         do l = 1,m
            sa = l*ds
            call potenl (100.d0,sa,0.d0,v(l),ia)
            if (v(l) .lt. v(lmid)) lmid = l
         enddo
         vmid(ia) = v(lmid)
         ll = 0
         do l = 1,lmid
            if (v(l) .lt. ered) go to 1
            ll = l
         enddo
   1     sum = 0.d0
         lmin = 0
         do l = ll,1,-1
            lmin = l
            sum = sum+ds*sqrt(v(l)-ered)
            if (sum .gt. scut) go to 2
         enddo
   2     smin(ia) = lmin*ds
         lr = m+1
         do l = m,lmid,-1
            if (v(l) .lt. ered) go to 3
            lr = l
         enddo
   3     sum = 0.0d0
         lmax = m+1
         do l = lr,m
            lmax = l
            sum = sum+ds*sqrt(v(l)-ered)
            if (sum .gt. scut) go to 4
         enddo
   4     smax = max(smax,lmax*ds)
      enddo
c
c     rmin
c
      pi = acos(-1.d0)
      rm2 = 0.d0
      rm2 = rm2+(1.d0-mass(1)/mtot)*smin(1)**2
      rm2 = rm2+(1.d0-mass(2)/mtot)*smin(2)**2
      rm2 = rm2+(1.d0-mass(3)/mtot)*smin(3)**2
      rmin = sqrt(rm2)
c
c     nvi
c
      nvi = 0            
      range = 2.d0*smax/pi
      do ia = 1,3
         if (ered .gt. vmid(ia)) then
            wave = sqrt(ered-vmid(ia))
            nvia = range*wave
            nvi = max(nvi,nvia)    
         endif
      enddo
c
c     mro and mvi
c
      mro = 3*(2*(jmax+1)+(nvi+1))/4
      mvi = 3*(2*(nvi+1)+(jmax+1))/4
c
c     parameter summary
c
      write (6,61) rmin,mro,mvi,nvi,smax
      if (rmin .ge. rmax) stop 'parset 1'
      return
  61  format(/1x,'PARSET:'/1x,70('-')/1x,
     +'  rmin = ',f32.2/1x, 
     +'   mro = ',i32/1x,
     +'   mvi = ',i32/1x,
     +'   nvi = ',i32/1x,
     +'  smax = ',f32.2/1x,70('-'))
      end     

      subroutine driver 
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Main driving subroutine for hyperspherical coordinate calculation
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      double precision mass,mtot,mred
      common /arrays/ mro,mvi,nvi,n
      common /energy/ enrg,dnrg,nnrg
      common /inputs/ emax,mtr
      common /masses/ mass(3),mtot,mred
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     global arrays
c
      parameter (nmax = 1500)
      dimension ered(nmax)
      dimension cvi(nvi,nmax),cro(3,nmax),eint(nmax)
      dimension ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax)
      dimension cvr(nvi,nmax)
c
c     rovibrational basis set
c
      call basis (cvi,cro,eint,ilev,jlev,klev,nlev,nmax)
c
c     scattering energies
c
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
  61        format(1x,'*** No open reaction channels at '
     +      ,f7.4,' eV. Energy skipped.')
         endif 
         enrg = enrg+dnrg
      enddo
      nnrg = nrg
      if (nnrg .gt. nmax) then
         write (6,99) nnrg,nmax
  99     format(/1x,'Error: nnrg = ',i6,' > nmax = ',i6/1x,
     +              'Run halted in subroutine driver')
         stop
      endif 
      if (nnrg .eq. 0) return
c
c     coupled channel propagation
c
      call setrec (n,nnrg)
      call solve (cvr,cro,ilev,jlev,klev,nlev,ered,nnrg)
c
c     asymptotic matching and output
c
      call match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
      call endrec
      return
      end

      subroutine basis (cvi,cro,eint,ilev,jlev,klev,nlev,nmax)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs an asymptotic multiple-arrangement 
c     rovibrational basis set containing all channels with j.le.jmax
c     and eint.le.emax.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /inputs/ emax,mtr
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     input arrays
c
      dimension cvi(nvi,nmax),cro(3,nmax),eint(nmax)
      dimension ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax)
c
c     local arrays
c
      dimension h(nvi,nvi),t(nvi),e(nvi)
      dimension s(0:2*nvi),v(0:2*nvi)
      dimension wvi(mvi),xvi(mvi)
c
c     vibrational quadrature rule
c
      smin = 0.d0
      call qvib (smin,smax,mvi,wvi,xvi)
c
c     kinetic energy 
c
      pi = acos(-1.d0)
      tscale = pi/smax 
      do jvi = 1,nvi
         t(jvi) = (jvi*tscale)**2
      enddo
      do jvi = 0,2*nvi
         s(jvi) = 0.d0
      enddo
      do kvi = 1,mvi
         sa = xvi(kvi)
         weight = wvi(kvi)/smax
         onsa2 = weight/sa**2
         arg = tscale*sa
         do jvi = 0,2*nvi 
            cosine = cos(jvi*arg) 
            s(jvi) = s(jvi)+cosine*onsa2
         enddo
      enddo
c
c     loop over arrangements
c
      n = 0
      ered = rmlmda*emax
      na = 3-iabs(jpar)
      do ia = 1,na
c
c        reference potential 
c
         do jvi = 0,2*nvi
            v(jvi) = 0.d0
         enddo
         do kvi = 1,mvi
            sa = xvi(kvi)
            weight = wvi(kvi)/smax
            call potenl (100.d0,sa,0.d0,va,ia)
            va = weight*va
            arg = tscale*sa
            do jvi = 0,2*nvi 
               cosine = cos(jvi*arg) 
               v(jvi) = v(jvi)+cosine*va
            enddo
         enddo
c
c        loop over rotational quantum numbers
c
         if (ia .eq. 1) then
            jmin = (1-jpar)/2
            jinc = 1+iabs(jpar)
         else
            jmin = 0
            jinc = 1
         endif
         do j = jmin,jmax,jinc
            kmaxt = min(j,jtot)
            kmaxj = min(j,kmax)
c
c           vibrational eigenvalue problem
c
            ctot = jtot*(jtot+1)
            cent = j*(j+1)
            do jvi = 1,nvi
               do ivi = 1,jvi
                  h(ivi,jvi) = cent*(s(jvi-ivi)-s(jvi+ivi))
     +                       +      (v(jvi-ivi)-v(jvi+ivi))
                  h(jvi,ivi) = h(ivi,jvi)
               enddo
               h(jvi,jvi) = h(jvi,jvi)+t(jvi)
            enddo
            call symevp (h,nvi,nvi,e,ierr)
            if (ierr .ne. 0) stop 'basis 1'
c
c           basis functions with eint .le. ered
c
            do kvi = 1,nvi
               if (e(kvi) .le. ered) then
                  do k = kmin,kmaxj
                     n = n+1
                     if (n .le. nmax) then
c
c                       vibrational expansion coefficients
c
                        do jvi = 1,nvi
                           cvi(jvi,n) = h(jvi,kvi)
                        enddo
c
c                       Coriolis coupling matrix elements
c
                        if (k .eq. kmin) then
                           cro(1,n) = 0.d0
                        else if (k .eq. 1) then
                           cro(1,n) = -sqrt(2.d0*ctot*cent)
                        else
                           ckm = k*(k-1)
                           cro(1,n) = -sqrt((ctot-ckm)*(cent-ckm))
                        endif
                        cro(2,n) = ctot+cent-2*k*k
                        if (k .eq. kmaxt) then
                           cro(3,n) = 0.d0
                        else if (k .eq. 0) then
                           cro(3,n) = -sqrt(2.d0*ctot*cent)
                        else 
                           ckp = k*(k+1)
                           cro(3,n) = -sqrt((ctot-ckp)*(cent-ckp))
                        endif
c
c                       quantum numbers
c
                        eint(n) = e(kvi)
                        ilev(n) = ia
                        jlev(n) = j
                        klev(n) = k
                        nlev(n) = kvi-1
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
      if (n .gt. nmax) then
         write (6,99) n,nmax
  99     format(/1x,'Error: n = ',i6,' > nmax = ',i6/1x,
     +              'Run halted in subroutine basis')
         stop
      endif 
c
c     channel ordering
c
      call basort (cvi,nvi,cro,eint,ilev,jlev,klev,nlev,n)
c
c     phase convention
c
      call bclean (cvi,nvi,n,smin,smax,mvi)
c
c     basis summary
c
      write (6,61) 
  61  format(/1x,'BASIS:'/1x,70('-')/1x,
     + ' Channel       a       v       j       k       Energy (eV) '/1x,
     + 70('-'))
      do i = 1,n
         temp = eint(i)/rmlmda
         write (6,62) i,ilev(i),nlev(i),jlev(i),klev(i),temp
  62     format(1x,5i8,f16.5)
      enddo
      write (6,63)
  63  format(1x,70('-'))
      return
      end

      subroutine basort (cvi,nvi,cro,eint,ilev,jlev,klev,nlev,n)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine sorts the channels in the input  
c     basis set in order of increasing a,v,j,k.
c     ----------------------------------------------------------------- 
c
      dimension cvi(nvi,n),cro(3,n)
      dimension eint(n),ilev(n),jlev(n),klev(n),nlev(n)
c
      do j = 1,n-1
         k = j
         do i = j+1,n
            if (ilev(i) .lt. ilev(k)) k = i
            if (ilev(i) .gt. ilev(k)) go to 1
            if (nlev(i) .lt. nlev(k)) k = i
            if (nlev(i) .gt. nlev(k)) go to 1
            if (jlev(i) .lt. jlev(k)) k = i
            if (jlev(i) .gt. jlev(k)) go to 1
            if (klev(i) .lt. klev(k)) k = i
   1        continue
         enddo
         if (k .ne. j) then
            do i = 1,nvi
               cswap = cvi(i,j)
               cvi(i,j) = cvi(i,k)
               cvi(i,k) = cswap
            enddo
            do i = 1,3
               cswap = cro(i,j)
               cro(i,j) = cro(i,k)
               cro(i,k) = cswap
            enddo 
            eswap = eint(j)
            eint(j) = eint(k)
            eint(k) = eswap
            iswap = ilev(j)
            ilev(j) = ilev(k)
            ilev(k) = iswap
            jswap = jlev(j)
            jlev(j) = jlev(k)
            jlev(k) = jswap
            kswap = klev(j)
            klev(j) = klev(k)
            klev(k) = kswap
            nswap = nlev(j)
            nlev(j) = nlev(k)
            nlev(k) = nswap
         endif
      enddo 
      return
      end

      subroutine bclean (cvi,nvi,n,smin,smax,mvi)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine removes the phase ambiguity from the vibrational
c     expansion coefficient array cvi by ensuring that each vibrational
c     eigenfunction is positive at its inner turning point.
c     -----------------------------------------------------------------  
c     
      dimension cvi(nvi,n) 
      dimension pvi(n,mvi)
c
      ds = (smax-smin)/mvi
      do k = 1,mvi
         sa = smin+(k-0.5d0)*ds
         call pvib (smin,sa,smax,cvi,nvi,n,pvi(1,k),0)
      enddo
      do 1 j = 1,n
         pmax = 0.d0
         do k = 1,mvi
            pmax = max(pmax,abs(pvi(j,k)))
         enddo
         ptest = 0.1d0*pmax
         do k = 1,mvi
            if (abs(pvi(j,k)) .gt. ptest) then
               if (pvi(j,k) .lt. 0.d0) then
                  do i = 1,nvi
                     cvi(i,j) = -cvi(i,j)
                  enddo
               endif
               go to 1
            endif
         enddo
   1  continue
      return
      end

      subroutine solve (cvi,cro,ilev,jlev,klev,nlev,ered,nnrg)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine solves the hyperspherical coordinate
c     coupled-channel equations between rho = rmin and rmax.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /inputs/ emax,mtr
      common /ranges/ rmin,rmax,smax
c
c     input arrays
c
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension ered(nnrg)
c
c     local arrays
c
      dimension cvp(nvi,n)
      dimension eint(n),e(n)
      dimension h(n,n),s(n,n),c(n,n)
c
c     loop over sectors
c
      npp = 0 
      drho = (rmax-rmin)/mtr
      rhop = rmin-0.5d0*drho
      do ktr = 1,mtr
c
c        surface eigenvalue problem
c
         rho = rhop+drho
         call sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
c
c        sector-to-sector transformation matrix
c
         if (npp .eq. 0) then
            call putrec (c,n,0)
         else
            call metric (rhop,cvp,npp,rho,cvi,np,ilev,jlev,klev,h,s,c)
         endif
c
c        log derivative propagation
c
         call logder (drho,ered,nnrg,h,s,c,e,n,np,npp)
c
c        setup for the next sector
c
         do j = 1,n
            do k = 1,nvi
               cvp(k,j) = cvi(k,j)
            enddo
         enddo
         npp = np
         rhop = rho
      enddo
c
c     transformation to final coupled-channel basis set at rmax
c
      npp = n+1
      rho = rmax
      call sbasis (rho,cvi,e,ilev,jlev,klev,nlev,1)
      call overlp (rhop,cvp,rho,cvi,ilev,jlev,klev,h)
c     call overlp (rho,cvi,rhop,cvp,ilev,jlev,klev,s)
c     do j = 1,n
c        do i = 1,n
c           h(i,j) = 0.5d0*(h(i,j)+s(j,i))
c        enddo
c     enddo
      call getrec (s,n,0)
      call dgemm  ('t','n',n,np,n,1.d0,h,n,s,n,0.d0,c,n)
      call logder (drho,ered,nnrg,h,s,c,e,n,np,npp)
      return
      end

      subroutine logder (drho,ered,nnrg,y,z,t,e,n,np,npp)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Constant reference potential log derivative propagator, 
c     with nnrg energies propagated simultaneously.
c     ----------------------------------------------------------------- 
c
c     input arrays
c
      dimension ered(nnrg)
      dimension y(n,n),z(n,n),t(n,n),e(n)
c
c     local arrays
c
      dimension y1(n),y2(n),y3(n),y4(n)
c
c     loop over energies
c
      do nrg = 1,nnrg
c
c        constant reference potential propagators
c
         if (npp .le. n) then
            h = drho
            eps = 1.d-16/(h*h)
            do i = 1,np
               psq = e(i)-ered(nrg)
               if (psq .gt. eps) then
                  p = sqrt(psq)
                  y1(i) = p/tanh(p*h)
                  y2(i) = p/sinh(p*h)
               else if (psq .lt. -eps) then
                  p = sqrt(-psq)
                  y1(i) = p/tan(p*h)
                  y2(i) = p/sin(p*h)
               else
                  y1(i) = 1.d0/h
                  y2(i) = 1.d0/h 
               endif
               y3(i) = y2(i)
               y4(i) = y1(i)
            enddo
         endif
c
c        log derivative propagation
c
         if (npp .le. 0) then
            do j = 1,np
               do i = 1,np
                  y(i,j) = 0.d0
               enddo
               y(j,j) = y4(j)
            enddo
         else if (npp .le. n) then
            call getrec (y,n,nrg)
            call dgemm ('n','n',npp,np,npp,1.d0,y,n,t,n,0.d0,z,n)
            call dgemm ('t','n',np, np,npp,1.d0,t,n,z,n,0.d0,y,n)
            do j = 1,np
               do i = 1,j
                  y(i,j) = 0.5d0*(y(i,j)+y(j,i))
               enddo
               y(j,j) = y(j,j)+y1(j)
            enddo
            call syminv (y,n,np,ierr)
            if (ierr .ne. 0) stop 'logder 1'
            do j = 1,np
               do i = 1,j
                  y(i,j) = -y3(i)*y(i,j)*y2(j)
                  y(j,i) = y(i,j)
               enddo
               y(j,j) = y(j,j)+y4(j)
            enddo
         else if (npp .gt. n) then
            call getrec (y,n,nrg)
            call dgemm ('n','t',np,n,np,1.d0,y,n,t,n,0.d0,z,n)
            call dgemm ('n','n',n, n,np,1.d0,t,n,z,n,0.d0,y,n)
            do j = 1,np
               do i = 1,j
                  y(i,j) = 0.5d0*(y(i,j)+y(j,i))
                  y(j,i) = y(i,j)
               enddo
            enddo
         endif
         call putrec (y,n,nrg)
      enddo 
      return
      end

      subroutine setrec (n,nnrg)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------   
c     This subroutine opens a direct access scratch file on unit 10.
c     ----------------------------------------------------------------- 
c
      diskmb = 8.0d-6*n*n*(nnrg+1)
      if (diskmb .lt. 1000.d0) then
         open (unit=10,status='scratch',form='unformatted',
     +         access='direct',recl=8*n*n)
         write (6,61) diskmb
      else
         nnrgmx = (1000.d0/diskmb)*(nnrg+1)-1
         write (6,62) nnrgmx
         stop
      endif
      return
  61  format(/1x,'SETREC:'/1x,70('-')/1x,
     + '*** This run is using ',f8.3,'  Mb of scratch disk space'/1x,
     + 70('-'))
  62  format(/1x,'SETREC:'/1x,70('-')/1x,
     + '*** This run needs more than 1 Gb of scratch disk space'/1x,
     + '*** Reduce nnrg to ',i4,' and try again'/1x,70('-'))
      end

      subroutine putrec (a,n,nrg)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine writes a matrix to the scratch file on unit 10.
c     ----------------------------------------------------------------- 
c
      dimension a(n,n)
      irec = nrg+1
      write (unit=10,rec=irec) a
      return
      end

      subroutine getrec (a,n,nrg)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine reads a matrix from the scratch file on unit 10. 
c     ----------------------------------------------------------------- 
c
      dimension a(n,n)
      irec = nrg+1
      read (unit=10,rec=irec) a
      return
      end

      subroutine endrec 
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------  
c     This subroutine closes and deletes the scratch file on unit 10.
c     ----------------------------------------------------------------- 
c
      close (unit=10,status='delete')
      return
      end

      subroutine sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine sets up and solves the hyperspherical 
c     coordinate surface eigenvalue problem at hyperradius rho. 
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /scales/ scale(3),rmlmda
c
c     input arrays
c
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension h(n,n),s(n,n),c(n,n),e(n)
c
c     local arrays
c
      dimension eint(n),sigma(n)
c
c     surface basis set
c
      call sbasis (rho,cvi,eint,ilev,jlev,klev,nlev,0)
c
c     hamiltonian and overlap matrices
c
      call couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
c
c     surface eigenvalue problem
c
      call hcsevp (h,s,c,n,e,sigma,np) 
c
c     sector summary
c
      write (6,61) e(1)/rmlmda,e(np)/rmlmda,np
  61  format(1x,'SECTOR:',
     +'  e(1) = ',f8.3,
     +'  e(np) = ',f8.3,
     +'    np = ',i8)
      return
      end 

      subroutine sbasis (rho,cvi,eint,ilev,jlev,klev,nlev,mode)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs a multiple-arrangement 
c     hyperspherical basis set at hyperradius rho.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c
c     input arrays
c
      dimension cvi(nvi,n)
      dimension eint(n),ilev(n),jlev(n),klev(n),nlev(n)
c
c     local arrays
c
      dimension h(nvi,nvi),t(nvi),e(nvi)
      dimension c(0:2*nvi),s(0:2*nvi),v(0:2*nvi)
      dimension wvi(mvi),xvi(mvi)
c
c     vibrational quadrature rule
c
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
c
c     kinetic energy 
c
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
c
c     loop over arrangements
c
      na = 3-iabs(jpar)
      do ia = 1,na
c
c        reference potential 
c
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
c
c        loop over rotational quantum numbers
c
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
c
c              vibrational eigenvalue problem
c
               if (mode .eq. 0) then
c
c                 k-dependent vibrational functions 
c                 in the interaction region
c
                  ctot = (jtot*(jtot+1)+j*(j+1)-2*k*k)/rho**2 
               else
c
c                 k-independent vibrational functions
c                 in the asymptotic region 
c
                  ctot = 0.d0
               endif
               cent = j*(j+1)/rho**2
               do jvi = 1,nvi
                  do ivi = 1,jvi
                     h(ivi,jvi) = ctot*(c(jvi-ivi)-c(jvi+ivi))
     +                          + cent*(s(jvi-ivi)-s(jvi+ivi))
     +                          +      (v(jvi-ivi)-v(jvi+ivi))
                     h(jvi,ivi) = h(ivi,jvi)
                  enddo
                  h(jvi,jvi) = h(jvi,jvi)+t(jvi)
               enddo
               call symevp (h,nvi,nvi,e,ierr)
               if (ierr .ne. 0) stop 'sbasis 1'
c
c              specified surface basis functions
c
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

      subroutine couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs the Hamiltonian and overlap matrices
c     h and s that are needed to solve the hyperspherical coordinate
c     surface eigenvalue problem at the given value of rho.
c     ----------------------------------------------------------------- 
c     
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     input arrays
c
      dimension cvi(nvi,n),cro(3,n)
      dimension eint(n),ilev(n),jlev(n),klev(n)
      dimension h(n,n),s(n,n)
c
c     local arrays
c
      dimension v(mro,mvi)
c
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
c
c     couple summary
c
      write (6,61) rho,heps/rmlmda,seps
  61  format(/1x,'COUPLE:',
     +'   rho = ',f8.3,
     +'   heps = ',f8.3,
     +'  seps = ',1p,e8.1)
      return
      end 

      subroutine direct (rho,cvi,cro,ilev,jlev,klev,h,v,ia)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine fills in the direct matrix elements of h 
c     within arrangement ia at hyperradius rho. Last modified 9/6/99.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c      
c     input arrays
c
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension h(n,n),v(mro,mvi)
c
c     local arrays
c
      dimension wro(mro),xro(mro)
      dimension wvi(mvi),xvi(mvi)
      dimension pro(0:jmax,0:kmax,mro)
      dimension pvi(mvi,n)
      dimension hro(0:jmax,0:jmax,0:kmax)
c
c     arrangement indices
c
      call arrang (ilev,n,jpar,ia,nla,nha,na) 
      if (na .lt. 1) return
c
c     quadrature rules
c
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
c
c     rotational basis functions
c
      do kro = 1,mro
         cosa = xro(kro)
         call prot (pro(0,0,kro),jmax,kmax,cosa)
      enddo
c
c     vibrational basis functions
c
      call pvibv (tmin,xvi,tmax,cvi(1,nla),nvi,na,pvi(1,nla),mvi,0)
c
c     loop over quadrature points
c
      do kvi = 1,mvi
         ta = xvi(kvi)
         ra = rho*cos(ta)
         sa = rho*sin(ta)
         call potenl (100.d0,sa,0.d0,va,ia)
c
c        rotational integrals
c
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
c
c        vibrational integrals and Coriolis coupling 
c
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
                  h(la,lb) = h(la,lb)+
     +                       pvi(kvi,la)*hro(ja,jb,k)*pvi(kvi,lb)
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

      subroutine exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine fills in the exchange matrix elements of h 
c     and s between arrangements ia and ib at hyperradius rho.
c     Last modified 9/6/99.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c      
c     input arrays
c
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension h(n,n),s(n,n),v(mro,mvi)
c
c     local arrays
c
      parameter (nblock = 32)
      dimension wro(mro),xro(mro)
      dimension wvi(mvi),xvi(mvi)
      dimension ab2(0:jmax),pvi(n)
      dimension pro(0:jmax,0:kmax+1)
      dimension dro(0:kmax+1,0:kmax+1)
      dimension prv(n,nblock,0:kmax+1)
      dimension hrv(nblock,0:jmax,0:kmax+1)
      dimension srv(nblock,0:jmax,0:kmax+1)
      dimension hro(n,0:jmax,0:kmax+1)
      dimension sro(n,0:jmax,0:kmax+1)
c
c     arrangement indices
c
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      call arrang (ilev,n,jpar,ib,nlb,nhb,nb)
      if (na.lt.1 .or. nb.lt.1) return
c
c     projection ranges (allowing for Coriolis terms)
c
      kamax = kmax
      if (kmax .lt. min(jtot,jmax)) kamax = kmax+1
      kbmax = kmax
c
c     A+B2 permutation symmetry
c
      do j = 0,jmax
         ab2(j) = 1.d0
         if (jpar .ne. 0) then
            if (ia+ib .eq. 3) ab2(j) = sqrt(2.d0)
            if (ia+ib .eq. 5) ab2(j) = jpar*(-1)**j
         endif
      enddo
c
c     quadrature rules
c
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
c
c     loop over quadrature points
c
      do kvi = 1,mvi
         inner = 0
         iblock = 0
         do kro = 1,mro
            ta = xvi(kvi)
            cosa = xro(kro)
c
c           coordinate transformation
c
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
c
c              jacobian factor and interaction potential
c
               weight = wro(kro)*wvi(kvi)
               sab = weight*sin(2.d0*ta)/sin(2.d0*tb)
               vab = sab*v(kro,kvi)
c
c              parity-adapted reduced rotation matrix elements
c
               do kb = kmin,kbmax
                  do ka = kmin,kamax
                     dro(ka,kb) = rotmel(jtot,ipar,ka,kb,xab)
                  enddo
               enddo
c
c              rovibrational functions in arrangement ib
c
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
c
c              rotational functions in arrangement ia
c
               call prot (pro,jmax,kamax,cosa)
               do ka = kmin,kamax
                  do ja = ka,jmax
                     hrv(iblock,ja,ka) = vab*pro(ja,ka)
                     srv(iblock,ja,ka) = sab*pro(ja,ka)
                  enddo
               enddo
c
c              partial potential and overlap matrix elements
c
               if (iblock .eq. nblock) then
                  do ka = kmin,kamax
                     nk = jmax-ka+1
                     call dgemm ('n','n',nb,nk,nblock,1.d0,
     +                           prv(nlb,1,ka),n,hrv(1,ka,ka),nblock,
     +                           1.d0,hro(nlb,ka,ka),n)
                     call dgemm ('n','n',nb,nk,nblock,1.d0,
     +                           prv(nlb,1,ka),n,srv(1,ka,ka),nblock,
     +                           1.d0,sro(nlb,ka,ka),n)
                  enddo
                  iblock = 0
               endif
            endif
         enddo
         if (iblock .gt. 0) then
            do ka = kmin,kamax
               nk = jmax-ka+1
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n,
     +                     hrv(1,ka,ka),nblock,1.d0,hro(nlb,ka,ka),n)
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n,
     +                     srv(1,ka,ka),nblock,1.d0,sro(nlb,ka,ka),n)
            enddo
         endif
         if (inner .gt. 0) then
c
c           vibrational functions in arrangement ia
c
            call pvib (tmin,ta,tmax,cvi(1,nla),nvi,na,pvi(nla),0)
c
c           full potential, overlap, and Coriolis matrix elements
c
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

      subroutine hcsevp (h,s,c,n,e,sigma,np) 
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine solves the surface eigenvalue problem 
c     (H-E*S)*C=0 by canonical orthogonalisation of the basis.
c     ----------------------------------------------------------------- 
c
      parameter (CUTOFF=1.d-04)
      dimension h(n,n),s(n,n),c(n,n)
      dimension e(n),sigma(n)
c
      stest = 0.d0
      do j = 1,n
         stestj = 0.d0
         do i = 1,n
            stestj = stestj+abs(s(i,j))
         enddo
         stestj = stestj+abs(1.d0-s(j,j))-abs(s(j,j))
         stest = max(stest,stestj)
      enddo
      if (stest .lt. 1.d-08) then
         np = n
         do j = 1,n
            do i = 1,n
               c(i,j) = h(i,j)
            enddo
         enddo
         call symevp (c,n,n,e,ierr)
         if (ierr .ne. 0) stop 'hcsevp 1'
      else
         do j = 1,n
            do i = 1,n
               s(i,j) = -s(i,j)
            enddo
         enddo
         call symevp (s,n,n,sigma,ierr)
         if (ierr .ne. 0) stop 'hcsevp 2' 
         np = 0
         do j = 1,n
            sigma(j) = -sigma(j)
            if (sigma(j) .gt. CUTOFF) then
               scale = 1.d0/sqrt(sigma(j))         
               do i = 1,n
                  s(i,j) = scale*s(i,j)      
               enddo
               np = j
            endif
         enddo
         if (np .gt. 0) then
            call dgemm ('n','n',n,np,n,1.d0,h,n,s,n,0.d0,c,n)
            call dgemm ('t','n',np,np,n,1.d0,s,n,c,n,0.d0,h,n)
            call symevp (h,n,np,e,ierr)
            if (ierr .ne. 0) stop 'hcsevp 3'
            call dgemm ('n','n',n,np,np,1.0d0,s,n,h,n,0.0d0,c,n)
         endif
      endif
      return
      end

      subroutine metric (rhop,cvp,npp,rho,cvi,np,ilev,jlev,klev,h,s,c)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs the sector-to-sector transformation
c     matrix t between sectors centred at hyperradii rhop and rho.
c     Compare with subroutine sector. 
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
c
c     input arrays
c
      dimension cvp(nvi,n),cvi(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension h(n,n),s(n,n),c(n,n)
c
c     primitive overlap integrals
c
      call overlp (rhop,cvp,rho,cvi,ilev,jlev,klev,h)
c     call overlp (rho,cvi,rhop,cvp,ilev,jlev,klev,s)
c     seps = 0.d0
c     do j = 1,n
c        do i = 1,n
c           seps = max(seps,abs(h(i,j)-s(j,i)))
c           h(i,j) = 0.5d0*(h(i,j)+s(j,i))
c        enddo
c     enddo
c
c     sector-to-sector transformation matrix
c
      call dgemm ('n','n',n,np,n,1.d0,h,n,c,n,0.d0,s,n)
      call getrec (h,n,0)
      call putrec (c,n,0)
      call dgemm ('t','n',npp,np,n,1.d0,h,n,s,n,0.d0,c,n)
c
c     metric summary
c
c     write (6,61) rhop,rho,seps
c 61  format(1x,'METRIC:',
c    +'  rhop = ',f8.3,
c    +'    rho = ',f8.3,
c    +'  seps = ',1p,e8.1)
      return
      end 

      subroutine overlp (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs the sector-to-sector overlap matrix
c     s between sectors centred at hyperradii rhoa and rhob.
c     Compare with subroutine couple.
c     ----------------------------------------------------------------- 
c     
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c
c     input arrays
c
      dimension cvia(nvi,n),cvib(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension s(n,n)
c
      do j = 1,n
         do i = 1,n
            s(i,j) = 0.d0
         enddo
      enddo
      na = 3-iabs(jpar) 
      do ia = 1,na
         call overdi (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia)
         do ib = 1,na
            if (ib .ne. ia) then
             call overex (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia,ib)
            endif
         enddo
      enddo
      if (na .eq. 2) then
         ia = 2
         ib = 3
         call overex (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia,ib)
      endif
      return
      end 

      subroutine overdi (rhoa,cva,rhob,cvb,ilev,jlev,klev,s,ia)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine fills in the matrix elements of the setor-
c     to-sector overlap matrix s within arrangement ia.
c     Compare with subroutine direct.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c      
c     input arrays
c
      dimension cva(nvi,n),cvb(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension s(n,n)
c
c     local arrays
c
      dimension wvi(mvi),xvi(mvi)
      dimension pva(mvi,n),pvb(mvi,n)
c
c     arrangement indices
c
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      if (na .lt. 1) return
c
c     vibrational quadrature rule
c
      tmin = 0.d0
      tamax = asin(min(1.d0,smax/rhoa))
      tbmax = asin(min(1.d0,smax/rhob))
      call qvib (tmin,tamax,mvi,wvi,xvi)
c
c     vibrational basis functions
c
      call pvibv (tmin,xvi,tamax,cva(1,nla),nvi,na,pva(1,nla),mvi,0)
      call pvibv (tmin,xvi,tbmax,cvb(1,nla),nvi,na,pvb(1,nla),mvi,0)
c
c     vibrational integrals
c
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

      subroutine overex (rhoa,cva,rhob,cvb,ilev,jlev,klev,s,ia,ib)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine fills in the matrix elements of the sector-
c     to-sector overlap matrix s between arrangements ia and ib.
c     Compare with subroutine exchng.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c      
c     input arrays
c
      dimension cva(nvi,n),cvb(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension s(n,n)
c
c     local arrays
c
      parameter (nblock = 32)
      dimension wro(mro),xro(mro)
      dimension wvi(mvi),xvi(mvi)
      dimension ab2(0:jmax),pvi(n)
      dimension pro(0:jmax,0:kmax)
      dimension dro(0:kmax,0:kmax)
      dimension prv(n,nblock,0:kmax)
      dimension srv(nblock,0:jmax,0:kmax)
      dimension sro(n,0:jmax,0:kmax)
c
c     arrangement indices
c
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      call arrang (ilev,n,jpar,ib,nlb,nhb,nb)
      if (na.lt.1 .or. nb.lt.1) return
c
c     A+B2 permutation symmetry
c
      do j = 0,jmax
         ab2(j) = 1.d0
         if (jpar .ne. 0) then
            if (ia+ib .eq. 3) ab2(j) = sqrt(2.d0)
            if (ia+ib .eq. 5) ab2(j) = jpar*(-1)**j
         endif
      enddo
c
c     quadrature rules
c
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tamax = asin(min(1.d0,smax/rhoa))
      tbmax = asin(min(1.d0,smax/rhob))
      call qvib (tmin,tamax,mvi,wvi,xvi)
c
c     loop over quadrature points
c
      do kvi = 1,mvi
         inner = 0
         iblock = 0
         do kro = 1,mro
            ta = xvi(kvi)
            cosa = xro(kro)
c
c           coordinate transformation
c
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
c
c              jacobian factor 
c
               weight = wro(kro)*wvi(kvi)
               sab = weight*sin(2.d0*ta)/sin(2.d0*tb)
c
c              parity-adapted reduced rotation matrix elements
c
               do kb = kmin,kmax
                  do ka = kmin,kmax
                     dro(ka,kb) = rotmel(jtot,ipar,ka,kb,xab)
                  enddo
               enddo
c
c              rovibrational functions in arrangement ib
c
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
c
c              rotational functions in arrangement ia
c
               call prot (pro,jmax,kmax,cosa)
               do ka = kmin,kmax
                  do ja = ka,jmax
                     srv(iblock,ja,ka) = sab*pro(ja,ka)
                  enddo
               enddo
c
c              partial overlap matrix elements
c
               if (iblock .eq. nblock) then
                  do ka = kmin,kmax
                     nk = jmax-ka+1
                     call dgemm ('n','n',nb,nk,nblock,1.d0,
     +                           prv(nlb,1,ka),n,srv(1,ka,ka),nblock,
     +                           1.d0,sro(nlb,ka,ka),n)
                  enddo
                  iblock = 0
               endif
            endif
         enddo
         if (iblock .gt. 0) then
            do ka = kmin,kmax
               nk = jmax-ka+1
               call dgemm ('n','n',nb,nk,iblock,1.d0,prv(nlb,1,ka),n,
     +                     srv(1,ka,ka),nblock,1.d0,sro(nlb,ka,ka),n)
            enddo
         endif
         if (inner .gt. 0) then
c
c           vibrational functions in arrangement ia
c
            call pvib (tmin,ta,tamax,cva(1,nla),nvi,na,pvi(nla),0)
c
c           full overlap matrix elements
c
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

      subroutine arrang (ilev,n,jpar,ia,nla,nha,na)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine calculates arrangement channel indices nla, nha
c     and na such that the na = nha-nla+1 channels in arrangement ia 
c     occur between nla and nha (assuming that the channels have been
c     sorted in a,v,j order by subroutine basort).
c     ----------------------------------------------------------------- 
c
      dimension ilev(n)
      dimension nlo(3),nhi(3)
c
      do i = 1,3
         nlo(i) = n+1
         nhi(i) = 0
      enddo
      do j = 1,n
         nlo(ilev(j)) = min(nlo(ilev(j)),j) 
         nhi(ilev(j)) = max(nhi(ilev(j)),j)
      enddo
      if (jpar .ne. 0) then
         nlo(3) = nlo(2)
         nhi(3) = nhi(2)
      endif
      nla = nlo(ia)
      nha = nhi(ia)
      na = nha-nla+1
      return
      end

      subroutine coords (ta,cosa,ia,tb,cosb,ib,xab)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------  
c     This subroutine calculates ib-arrangement from ia-arrangement 
c     Delves hyperangular coordinates, along with the angle xab 
c     between the body-frame z-axes of the two arrangements. 
c     -----------------------------------------------------------------  
c
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
      common /scales/ scale(3),rmlmda
c
      if (ib .eq. ia) then
         tb = ta
         cosb = cosa
         xab = 0.d0
      else
         sab = 1.d0/(scale(ia)*scale(ib))
         cab = -mred*sab/mass(6-ia-ib)
         if (ia-ib.eq.1 .or. ib-ia.eq.2) sab = -sab
         ra = cos(ta)
         sa = sin(ta)
         rsca = ra*sa*cosa
         rb = sqrt((cab*ra)**2-2.d0*cab*sab*rsca+(sab*sa)**2)
         sb = sqrt((sab*ra)**2+2.d0*cab*sab*rsca+(cab*sa)**2)
         rscb = cab*sab*(ra-sa)*(ra+sa)+rsca*(cab-sab)*(cab+sab)
         cosb = rscb/(rb*sb)
         cosb = max(cosb,-1.d0)
         cosb = min(cosb,+1.d0)
         tb = atan(sb/rb)
         cxab = (rb*cab+sb*cosb*sab)/ra
         cxab = max(cxab,-1.d0)
         cxab = min(cxab,+1.d0)
         xab = acos(cxab)
         if (ia-ib.eq.1 .or. ib-ia.eq.2) xab = -xab 
      endif
      return
      end

      subroutine potenl (ra,sa,cosa,va,ia)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     This subroutine returns the potential energy surface va as a 
c     function of mass-scaled Jacobi coordinates in arrangement ia.
c     ----------------------------------------------------------------- 
c
      dimension r(3)
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
      common /scales/ scale(3),rmlmda
c
      ib = ia+1
      if (ib .gt. 3) ib = 1
      ic = 6-ia-ib
      rap = ra/scale(ia)
      cap = 2.0d0*cosa*rap
      sap = scale(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap 
      r(ia) = sap
      r(ib) = sqrt(rap**2-cap*sbp+sbp**2)
      r(ic) = sqrt(rap**2+cap*scp+scp**2)
      call potsub (r,vev)
      va = rmlmda*vev
      return
      end

      subroutine pvib (a,x,b,c,n,m,p,mode)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine calculates m normalised vibrational 
c     eigenfunctions in a < x < b (if mode = 0), or their  
c     first derivatives with respect to x (if mode = 1).
c     ----------------------------------------------------------------- 
c
      dimension c(n,m),p(m)
      dimension q(0:n)
c
c     primitive sine basis functions
c
      pi = acos(-1.d0)
      scale = sqrt(2.d0/(b-a))
      shift = pi/(b-a)
      y = shift*(x-a)
      cosy = cos(y)
      recur = 2*cosy
      if (mode .eq. 0) then
         q(0) = 0.d0
         q(1) = scale*sin(y)
         do j = 2,n
            q(j) = recur*q(j-1)-q(j-2)
         enddo
      else 
         scale = shift*scale
         q(0) = scale
         q(1) = scale*cosy
         do j = 2,n
            q(j) = recur*q(j-1)-q(j-2)
         enddo
         do j = 2,n
            q(j) = j*q(j)
         enddo
      endif
c
c     contracted vibrational eigenfunctions
c
      call dgemv ('t',n,m,1.d0,c,n,q(1),1,0.d0,p,1)
      return
      end 

      subroutine pvibv (a,x,b,c,n,m,p,l,mode)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Vector version of subroutine pvib for l different values of x.
c     ----------------------------------------------------------------- 
c
      dimension x(l),c(n,m),p(l,m)
      dimension q(l,0:n)
c
c     primitive sine basis functions
c
      pi = acos(-1.d0)
      shift = pi/(b-a)
      scale = sqrt(2.d0/(b-a))
      if (mode .eq. 0) then
         do k = 1,l
            y = shift*(x(k)-a)
            recur = 2*cos(y)
            q(k,0) = 0.d0
            q(k,1) = scale*sin(y)
            do j = 2,n
               q(k,j) = recur*q(k,j-1)-q(k,j-2)
            enddo
         enddo
      else 
         scale = shift*scale
         do k = 1,l
            y = shift*(x(k)-a)
            cosy = cos(y)
            recur = 2*cosy
            q(k,0) = scale
            q(k,1) = scale*cosy
            do j = 2,n
               q(k,j) = recur*q(k,j-1)-q(k,j-2)
            enddo
            do j = 2,n
               q(k,j) = j*q(k,j)
            enddo
         enddo
      endif
c
c     contracted vibrational eigenfunctions
c
      call dgemm ('n','n',l,m,n,1.d0,q(1,1),l,c,n,0.d0,p,l)
      return
      end 

      subroutine qvib (a,b,m,w,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs an m point particle-in-a-box
c     (midpoint) quadrature rule in the interval a < x < b.
c     ----------------------------------------------------------------- 
c     
      dimension w(m),x(m)
c
      dx = (b-a)/m
      do k = 1,m
         w(k) = dx
         x(k) = a+(k-0.5d0)*dx
      enddo
      return
      end

      subroutine prot (p,jmax,kmax,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses the algorithm in Section 6.8 of Numerical
c     Recipes to calculate an array of normalised associated Legendre 
c     polynomials P(j,k,x) = sqrt(2*pi)*Y(j,k;acos(x),0).
c     ----------------------------------------------------------------- 
c
      parameter (jtop = 100)
      dimension p(0:jmax,0:kmax) 
      dimension a(0:jtop,0:jtop),b(0:jtop,0:jtop),c(0:jtop,0:jtop)
      save init,a,b,c
      data init /0/
c      
      if (init .eq. 0) then
         do j = 0,jtop
            do k = 0,j-2
               a(j,k) = (j+j-1.d0)/(j-k)
               b(j,k) = (j+k-1.d0)/(j-k)
            enddo
            fac = j+0.5d0
            c(j,0) = sqrt(fac)
            do k = 1,j
               fac = fac/((j+k)*(j-k+1))
               c(j,k) = sqrt(fac)
            enddo
         enddo
         init = 1 
      endif
c
      if (jmax .gt. jtop) stop 'prot 0'
      if (abs(x) .gt. 1.d0) stop 'prot 1'
c
      do j = 0,jmax
         do k = j+1,kmax
            p(j,k) = 0.d0
         enddo
      enddo
      p(0,0) = 1.d0
      if (kmax .gt. 0) then
         sx = -sqrt((1.d0-x)*(1.d0+x))
         do k = 1,min(jmax,kmax)
            p(k,k) = (2*k-1)*sx*p(k-1,k-1)
         enddo
      endif
      if (jmax .gt. 0) then
         p(1,0) = x
         do k = 1,min(jmax-1,kmax)
            p(k+1,k) = (2*k+1)*x*p(k,k)
         enddo
      endif
      do k = 0,kmax
         do j = k+2,jmax
            p(j,k) = a(j,k)*x*p(j-1,k)-b(j,k)*p(j-2,k)
         enddo
      enddo
      do k = 0,kmax
         do j = k,jmax
            p(j,k) = c(j,k)*p(j,k)
         enddo
      enddo
      return
      end

      subroutine qrot (m,w,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine constructs an m point Gauss-Legendre 
c     quadrature rule in the interval -1 < x < 1.
c     ----------------------------------------------------------------- 
c
      dimension w(m),x(m)
      dimension a2(m),a3(m)
c
      do j = 1,m
         a2(j) = (j+j-1.d0)/j
         a3(j) = (j-1.d0)/j
      enddo
      n = (m+1)/2
      pi = acos(-1.d0)
      do k = 1,n
         z = cos(pi*(k-0.25d0)/(m+0.5d0))
         do i = 1,6
            p2 = 0.d0
            p1 = 1.d0
            do j = 1,m
               p3 = p2
               p2 = p1
               p1 = a2(j)*z*p2-a3(j)*p3
            enddo
            p2 = m*(p2-z*p1)/(1.d0-z*z)
            z = z-p1/p2
         enddo
         x(k) = -z
         x(m+1-k) = +z
         w(k) = 2.d0/((1.d0-z*z)*p2*p2)
         w(m+1-k) = w(k)
      enddo
      return
      end

      subroutine match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
      implicit double precision (a-h,o-z)
      double precision llev
c
c     -----------------------------------------------------------------  
c     This subroutine uses the asymptotic log derivative matrix
c     at each energy to calculate and output the scattering matrix.
c     ----------------------------------------------------------------- 
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /output/ nout,jout
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
c
c     input arrays
c
      dimension cvr(nvi,n),cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension eint(n),ered(nnrg)
c
c     local arrays
c
      dimension x(n,n),y(n,n),z(n,n)
      dimension llev(n),crp(3,nnrg)
c
c     loop over energies
c
      do 3 nrg = 1,nnrg
c
c        BF to SF transformation of Y
c
         call getrec (y,n,nrg)
         call frames (y,z,cro,jlev,klev,llev,+1)
c
c        asymptotic matching
c
         enrg = ered(nrg)
         call bessel (x,y,z,cvr,cvi,eint,ilev,jlev,klev,llev,nlev,enrg)
         call smatrx (x,y,z,eint,llev,n,seps,ierr,enrg)
         if (ierr .ne. 0) stop 'match 1'
         enrg = enrg/rmlmda
c
c        SF to BF transformation of S
c
         call frames (x,z,cro,jlev,klev,llev,-1)
         call frames (y,z,cro,jlev,klev,llev,-1)
c
c        phase convention for parity-adapted BF S-matrix
c
         if (ipar .eq. -1) then
            do j = 1,n
               do i = 1,n
                  x(i,j) = -x(i,j)
                  y(i,j) = -y(i,j)
               enddo
            enddo
         endif
c
c        cumulative reaction probabilities
c
         do ia = 1,3
            crp(ia,nrg) = 0.d0
         enddo
         do i = 1,n
            if (ered(nrg).gt.eint(i) .and. ilev(i).eq.1) then
               do j = 1,n
                  if (ered(nrg).gt.eint(j) .and. ilev(j).ne.1) then
                     tp = x(i,j)**2+y(i,j)**2
                     crp(ilev(j),nrg) = crp(ilev(j),nrg)+tp
                  endif
               enddo
            endif
         enddo 
c
c        OUTPUT:
c        This outputs S-matrix elements and transition probabilities 
c        from all open A+BC(v,j) reactant channels with v.le.nout
c        and j.le.jout to all open A+BC, B+CA, C+AB product channels: 
c
         unit = 0.d0
         cols = 0.d0
         write (6,61) enrg,seps
         do 2 i = 1,n
            if (ered(nrg) .le. eint(i)) go to 2
            ii = ilev(i)
            ji = jlev(i)
            ki = klev(i)
            ni = nlev(i)
            if (ii .ne. 1) go to 2
            if (ni .gt. nout) go to 2
            if (ji .gt. jout) go to 2
            jp = -1
            sumj = 0.d0
            cols = cols+1.d0
            do 1 j = 1,n
               if (ered(nrg) .le. eint(j)) go to 1
               if = ilev(j)
               jf = jlev(j)
               kf = klev(j)
               nf = nlev(j)
               tp = x(i,j)**2+y(i,j)**2 
               unit = unit+tp
               if (jf .lt. jp) then
                  write (6,63) ii,ni,ji,ki,ip,np,sumj
                  sumj = 0.d0
               endif 
               ip = if
               np = nf
               jp = jf
               sumj = sumj+tp
               write (6,62) ii,ni,ji,ki,if,nf,jf,kf,x(i,j),y(i,j),tp
   1        continue
            write (6,63) ii,ni,ji,ki,ip,np,sumj
   2     continue 
         if (cols .gt. 0.d0) then
            unit = unit/cols
            write (6,64) unit
         endif
   3  continue 
c
c     And this outputs cumulative reaction probabilities:
c
      if (jpar .eq. 0) then
         write (6,65)
         do nrg = 1,nnrg
            enrg = ered(nrg)/rmlmda
            write (6,66) enrg,crp(2,nrg),crp(3,nrg)
         enddo
         write (6,67)
      else
         write (6,68)
         do nrg = 1,nnrg
            enrg = ered(nrg)/rmlmda
            write (6,69) enrg,crp(2,nrg)
         enddo
         write (6,67)
      endif
      return
  61  format(/1x,'MATCH:'/1x,70('-')/1x,
     + ' E(eV) = ',f8.5,27x,'Seps = ',1p,e8.1,0p/1x,70('-')/1x,
     + '  a    v    j    k    a''   v''   j''   k''     Re(S)',
     + '     Im(S)     |S|^2 '/1x,70('-'))
  62  format(1x,i3,7i5,1x,3f10.5)
  63  format(1x,i3,5i5,5x,'  all ',f30.5)
  64  format(1x,70('-')/1x,29x,'Unitarity',f31.5/1x,70('-'))
  65  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x,
     + '       E(eV)          A+BC --> B+CA          A+BC --> C+AB '
     + /1x,70('-'))
  66  format(1x,f13.5,f19.5,f23.5)
  67  format(1x,70('-'))
  68  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x,
     +'       E(eV)                                     A+B2 --> B+AB '
     + /1x,70('-'))
  69  format(1x,f13.5,f46.5)
      end

      subroutine frames (x,c,cro,jlev,klev,llev,mode)
      implicit double precision (a-h,o-z)
      double precision llev
c
c     ----------------------------------------------------------------- 
c     This subroutine transforms the symmetric matrix x from a 
c     (possibly truncated) BF basis set to a (possibly truncated) 
c     SF basis set if mode = +1, or vice versa if mode = -1.
c     The (possibly fractional) orbital angular momentum quantum
c     numbers of the SF basis are returned in llev. 
c     ----------------------------------------------------------------- 
c     
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c
c     input arrays
c
      dimension x(n,n),c(n,n)
      dimension cro(3,n)
      dimension jlev(n),klev(n),llev(n)
c
c     local arrays
c
      dimension d(n),e(n)
c
c     loop over blocks of channels with the same a,v, and j,
c     assuming an a,v,j,k (or a,v,j,l) channel ordering
c
      nlo = 1
   1  continue
         j = jlev(nlo)
         kmaxj = min(j,jtot,kmax)
         ndo = kmaxj-kmin+1
         nhi = nlo+ndo-1
         if (nhi .gt. n) stop 'frames 1'
c
c        form the BF to SF transformation matrix by diagonalising
c        the operator L^2 stored in cro, making sure to apply a
c        consistent phase convention
c
         do i = nlo,nhi
            d(i) = cro(2,i)
            e(i) = cro(1,i)
         enddo
         call rstevp (d(nlo),e(nlo),ndo,c(nlo,nlo),n,ierr)
         if (ierr .ne. 0) stop 'frames 2'
         phase = (-1)**(j+kmaxj)
         do k = nlo,nhi
            if (phase*c(nhi,k) .lt. 0.d0) then
               do i = nlo,nhi
                  c(i,k) = -c(i,k)
               enddo
            endif
            llev(k) = sqrt(d(k)+0.25d0)-0.5d0
            if (kmaxj .eq. min(j,jtot)) llev(k) = nint(llev(k))
         enddo
c
c        transform the present block of the matrix x, applying
c        the inverse transformation if mode = -1
c
         if (ndo .gt. 1) then 
            if (mode .lt. 0) then
               do j = nlo+1,nhi
                  do i = nlo,j-1
                     cswap = c(i,j)
                     c(i,j) = c(j,i)
                     c(j,i) = cswap
                  enddo
               enddo
            endif
            do i = 1,n
               do j = nlo,nhi
                  d(j) = 0.d0
                  do k = nlo,nhi
                     d(j) = d(j)+x(i,k)*c(k,j)
                  enddo
               enddo
               do j = nlo,nhi
                  x(i,j) = d(j)
               enddo
            enddo
            do j = 1,n
               do i = nlo,nhi
                  d(i) = 0.d0
                  do k = nlo,nhi
                     d(i) = d(i)+c(k,i)*x(k,j)
                  enddo
               enddo
               do i = nlo,nhi
                  x(i,j) = d(i)
               enddo
            enddo
         endif
c
c        and increment nlo for the next block
c
         nlo = nhi+1
      if (nlo .le. n) go to 1
      return
      end 

      subroutine bessel (x,y,z,cvr,cvi,eint,
     +                   ilev,jlev,klev,llev,nlev,ered)
      implicit double precision (a-h,o-z)
      double precision llev
c
c     -----------------------------------------------------------------  
c     This subroutine forms the asymptotic solution and derivative
c     matrices a,b,c, and d needed for matching the asymptotic Delves
c     coordinate log derivative matrix y onto a reactance matrix k,
c     as in eqs. (116) to (121) of Pack and Parker, and then proceeds
c     to calculate the reactance matrix in the array x.
c     -----------------------------------------------------------------  
c
c     common blocks
c
      common /arrays/ mro,mvi,nvi,n
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
c
c     input arrays
c
      dimension x(n,n),y(n,n),z(n,n)
      dimension cvr(nvi,n),cvi(nvi,n),eint(n)
      dimension ilev(n),jlev(n),klev(n),llev(n),nlev(n)
c
c     local arrays
c
      dimension wvi(mvi),xvi(mvi),pvi(nvi,2)
      dimension fr(n),fvi(n),gvi(n)
      dimension fa(n),fb(n),fc(n),fd(n)
      dimension a(n,0:nvi-1),b(n,0:nvi-1)
      dimension erow(n),ecol(n),srow(n),scol(n)
c
c     initialisation
c
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
c
c     vibrational quadrature rule
c
      pi = acos(-1.d0)
      piby2 = 0.5d0*pi
      rtrho = sqrt(rmax)
      tworho = 2.d0*rmax
      smin = 0.d0
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rmax))
      call qvib (tmin,tmax,mvi,wvi,xvi)
c
c     Delves vibrational quadrature
c
      do kvi = 1,mvi
         weight = rtrho*wvi(kvi)
         ta = xvi(kvi)
         cta = cos(ta)
         sta = sin(ta)
         ra = rmax*cta
         sa = rmax*sta
         if (sa .lt. smax) then
c
c           Delves vibrational functions
c
            call pvib (tmin,ta,tmax,cvr,nvi,n,fr,0)
            do i = 1,n
               fr(i) = weight*fr(i)
            enddo
c
c           Jacobi vibrational functions
c        
            call pvib (smin,sa,smax,cvi,nvi,n,fvi,0)
            call pvib (smin,sa,smax,cvi,nvi,n,gvi,1)
c
c           Jacobi translational (Riccati-Bessel) functions
c
            do i = 1,n
               ell = llev(i)
               psq = ered-eint(i)
               if (psq .gt. 0.d0) then
                  p = sqrt(psq)
                  arg = p*ra
                  call rbesjy (ell,arg,cj,dj,ej,cy,dy,ey)
c
c                 exponential scaling of j and y  
c                 (for stability near channel thresholds)
c 
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
c
c                 exponential scaling of k 
c                 (for the same reason)
c
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
c
c              Pack and Parker eqs. (116) to (121)
c
               fa(i) = atr*fvi(i)
               fb(i) = btr*fvi(i)
               fc(i) = cta*ctr*fvi(i)+sta*atr*gvi(i)
               fd(i) = cta*dtr*fvi(i)+sta*btr*gvi(i)
               fc(i) = fa(i)/tworho+fc(i)
               fd(i) = fb(i)/tworho+fd(i)
            enddo
c
c           integral accumulation
c
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
c
c     y to k
c
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
c
c     elimination of exponential scaling factors
c     from the reactance matrix
c
      do j = 1,n
         do i = 1,n
            x(i,j) = srow(i)*x(i,j)*scol(j)
         enddo
      enddo
      return
      end

      subroutine smatrx (x,y,z,eint,llev,n,seps,ierr,ered)
      implicit double precision (a-h,o-z)
      double precision llev
c
c     ----------------------------------------------------------------- 
c     This subroutine uses the reactance matrix k returned by 
c     subroutine bessel to calculate the scattering matrix s, the 
c     real and imaginary parts of which are returned in x and y.
c     ----------------------------------------------------------------- 
c
      dimension x(n,n),y(n,n),z(n,n)
      dimension eint(n),llev(n),t(n)
      dimension arg(n)
c
c     k to k(open,open)
c
      j = 0
      sum = 1.d0
      dif = 0.d0
      do jc = 1,n
         if (ered .gt. eint(jc)) then
            j = j+1
            i = 0
            do ic = 1,jc
               if (ered .gt. eint(ic)) then
                  i = i+1
                  sum = sum+abs(x(ic,jc)+x(jc,ic))
                  dif = dif+abs(x(ic,jc)-x(jc,ic))
                  z(i,j) = 0.5d0*(x(ic,jc)+x(jc,ic))
                  z(j,i) = z(i,j)
               endif
            enddo
         endif
      enddo
      seps = dif/sum
      ierr = 0
      nopen = j
      if (nopen .le. 0) return
c
c     k(open,open) to s
c
      call symevp (z,n,nopen,t,ierr)
      if (ierr .ne. 0) return
      do jc = 1,n
         do ic = 1,jc
            x(ic,jc) = 0.d0
            y(ic,jc) = 0.d0
         enddo
      enddo
      do k = 1,nopen
         denom = 1.d0+t(k)**2
         sr = (1.d0-t(k)**2)/denom
         si = 2.d0*t(k)/denom
         j = 0
         do jc = 1,n
            if (ered .gt. eint(jc)) then
               j = j+1
               zr = z(j,k)*sr
               zi = z(j,k)*si
               i = 0
               do ic = 1,jc
                  if (ered .gt. eint(ic)) then
                     i = i+1
                     x(ic,jc) = x(ic,jc)+z(i,k)*zr
                     y(ic,jc) = y(ic,jc)+z(i,k)*zi
                  endif
               enddo
            endif
         enddo
      enddo
      do jc = 2,n
         do ic = 1,jc-1
            x(jc,ic) = x(ic,jc)
            y(jc,ic) = y(ic,jc)
         enddo
      enddo
c
c     phase factor of (-i)**(l+l') for SF to BF transformation:
c     (see Eq. (27) of my Asymptotic Matching I notes)
c
      pi = acos(-1.d0)
      piby2 = 0.5d0*pi
      do j = 1,n
         if (ered .gt. eint(j)) then
            arg(j) = -llev(j)*piby2
            do i = 1,j
               if (ered .gt. eint(i)) then
                  aij = arg(i)+arg(j)
                  cij = cos(aij)
                  sij = sin(aij)
                  xij = cij*x(i,j)-sij*y(i,j)
                  yij = cij*y(i,j)+sij*x(i,j)
                  x(i,j) = xij
                  y(i,j) = yij
                  x(j,i) = xij
                  y(j,i) = yij
               endif
            enddo
         endif
      enddo
      return
      end
