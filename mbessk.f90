!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: mbessk.f

      subroutine mbessk (v,x,ck,dk,ek)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses Temme's method [ N.M.Temme, J Comput Phys
!     19 (1975) 324-337 ] to calculate the Modified Bessel function
!
!     K(v,x) = ck BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
!
!     and its first derivative with respect to x
!
!     d/dx K(v,x) = dk BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
!
!     for a given real order v >= 0 and real argument x > 0.
!     Note the exponential scaling, which is used to avoid
!     overflow of K(v,x) for v >> x and underflow for v << x.
!     -----------------------------------------------------------------
!
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
!
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'mbessk 0'
      pi = acos(-1.d0)
      xmin = 1.d0
!
!     begin by calculating K(a,x) and K(a+1,x) for |a| <= 1/2
!
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
!
!        using Temme's series for small x
!
         b = x/2.d0
         d = -dlog(b)
         e = a*d
         c = a*pi
         if (abs(c) .lt. eps) then
            c = 1.d0
         else
            c = c/sin(c)
         endif
         if (abs(e) .lt. eps) then
            s = 1.d0
         else
            s = sinh(e)/e
         endif
         e = exp(e)
         g = e*rgamma(a,p,q)
         e = (e+1.d0/e)/2.d0
         f = c*(p*e+q*s*d)
         e = a*a
         p = 0.5d0*g*c
         q = 0.5d0/g
         c = 1.d0
         d = b*b
         ak = f
         ak1 = p
         do n = 1,maxit
            f = (f*n+p+q)/(n*n-e)
            c = c*d/n
            p = p/(n-a)
            q = q/(n+a)
            g = c*(p-n*f)
            h = c*f
            ak = ak+h
            ak1 = ak1+g
            if (h/ak+abs(g)/ak1 .lt. eps) go to 1
         enddo
         stop 'mbessk 1'
   1     f = ak
         g = ak1/b
         ex = 0.d0
      else if (x .ge. xmin) then
!
!        and Temme's PQ method for large x
!
         c = 0.25d0-a*a
         g = 1.d0
         f = 0.d0
         e = x*cos(a*pi)/pi/eps
         do n = 1,maxit
            h = (2*(n+x)*g-(n-1+c/n)*f)/(n+1)
            f = g
            g = h
            if (h*n .gt. e) go to 2
         enddo
         stop 'mbessk 2'
   2     p = f/g
         q = p
         b = x+x
         e = b-2.d0
         do m = n,1,-1
            p = (m-1+c/m)/(e+(m+1)*(2.d0-p))
            q = p*(q+1.d0)
         enddo
         f = sqrt(pi/b)/(1.d0+q)
         g = f*(a+x+0.5d0-p)/x
         ex = x
      endif
!
!     now recur upwards from K(a,x) to K(v,x),
!     scaling to avoid overflow along the way
!
      p = 0.d0
      if (na .gt. 0) then
         y = 2.d0/x
         do n = 1,na
            h = y*(a+n)*g+f
            f = g
            g = h
   3        if (abs(f) .gt. 4.d0) then
               p = p+1.d0
               f = 0.0625d0*f
               g = 0.0625d0*g
               go to 3
            endif
         enddo
      endif
      ck = f
      dk = (v/x)*f-g
      sk = sqrt(ck*ck+dk*dk)
      ck = ck/sk
      dk = dk/sk
      ek = dlog(sk)+p*dlog(16.d0)-ex
      return
      end
