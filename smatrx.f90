!     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: smatrx.f

      subroutine smatrx (x,y,z,eint,llev,n,seps,ierr,ered)
      implicit double precision (a-h,o-z)
      double precision llev
!
!     -----------------------------------------------------------------
!     This subroutine uses the reactance matrix k returned by
!     subroutine bessel to calculate the scattering matrix s, the
!     real and imaginary parts of which are returned in x and y.
!     -----------------------------------------------------------------
!
      dimension x(n,n),y(n,n),z(n,n)
      dimension eint(n),llev(n),t(n)

!     dimension arg(n)
      allocatable arg(:)

      allocate (arg(n))

!
!     k to k(open,open)
!
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
!
!     k(open,open) to s
!
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
!
!     phase factor of (-i)**(l+l') for SF to BF transformation:
!     (see Eq. (27) of my Asymptotic Matching I notes)
!
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
