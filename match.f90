!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: match.f

      subroutine match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
      implicit double precision (a-h,o-z)
      double precision llev
!
!     -----------------------------------------------------------------
!     This subroutine uses the asymptotic log derivative matrix
!     at each energy to calculate and output the scattering matrix.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /output/ nout,jout
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
!
!     input arrays
!
      dimension cvr(nvi,n),cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension eint(n),ered(nnrg)
!
!     local arrays
!
!     dimension x(n,n),y(n,n),z(n,n)
      allocatable x(:,:),y(:,:),z(:,:)

!     dimension llev(n),crp(3,nnrg)
      allocatable llev(:),crp(:,:)

      allocate (x(n,n),y(n,n),z(n,n))
      allocate (llev(n),crp(3,nnrg))

!
!     loop over energies
!
      do 3 nrg = 1,nnrg
!
!        BF to SF transformation of Y
!
         call getrec (y,n,nrg)
         call frames (y,z,cro,jlev,klev,llev,+1)
!
!        asymptotic matching
!
         enrg = ered(nrg)
         call bessel (x,y,z,cvr,cvi,eint,ilev,jlev,klev,llev,nlev,enrg)
         call smatrx (x,y,z,eint,llev,n,seps,ierr,enrg)
         if (ierr .ne. 0) stop 'match 1'
         enrg = enrg/rmlmda
!
!        SF to BF transformation of S
!
         call frames (x,z,cro,jlev,klev,llev,-1)
         call frames (y,z,cro,jlev,klev,llev,-1)
!
!        phase convention for parity-adapted BF S-matrix
!
         if (ipar .eq. -1) then
            do j = 1,n
               do i = 1,n
                  x(i,j) = -x(i,j)
                  y(i,j) = -y(i,j)
               enddo
            enddo
         endif
!
!        cumulative reaction probabilities
!
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
!
!        OUTPUT:
!        This outputs S-matrix elements and transition probabilities
!        from all open A+BC(v,j) reactant channels with v.le.nout
!        and j.le.jout to all open A+BC, B+CA, C+AB product channels:
!
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
!
!     And this outputs cumulative reaction probabilities:
!
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
  61  format(/1x,'MATCH:'/1x,70('-')/1x, &
       ' E(eV) = ',es20.12,24x,'Seps = ',1p,e8.1,0p/1x,70('-')/1x, &
       '  a    v    j    k    a''   v''   j''   k''     Re(S)', &
       '     Im(S)     |S|^2 '/1x,70('-'))
  62  format(1x,i3,7i5,1x,3f10.5)
  63  format(1x,i3,5i5,5x,'  all ',f30.5)
  64  format(1x,70('-')/1x,29x,'Unitarity',f31.5/1x,70('-'))
  65  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x, &
       '       E(eV)          A+BC --> B+CA          A+BC --> C+AB ' &
       /1x,70('-'))
  66  format(1x,f13.5,f19.5,f23.5)
  67  format(1x,70('-'))
  68  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x, &
      '       E(eV)                                     A+B2 --> B+AB ' &
       /1x,70('-'))
  69  format(1x,f13.5,f46.5)
      end
