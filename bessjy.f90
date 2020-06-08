!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: bessjy.f

      subroutine bessjy (v,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses a combination of methods (mostly due
!     to Temme) to calculate the Ordinary Bessel functions
!
!     J(v,x) = cj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
!     Y(v,x) = cy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
!
!     and their first derivatives with respect to x
!
!     d/dx J(v,x) = dj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
!     d/dx Y(v,x) = dy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
!
!     for a given real order v >= 0 and real argument x > 0.
!     Note the exponential scaling, which is used to avoid
!     overflow of Y(v,x) and underflow of J(v,x) for v >> x.
!     -----------------------------------------------------------------
!
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
!
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'bessjy 0'
      pi = acos(-1.d0)
      xmin = 3.d0
      xmax = 5.d0-dlog10(eps)
!
!     begin by calculating Y(a,x) and Y(a+1,x) for |a| <= 1/2
!
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
!
!        using Temme's series (bessya) for small x
!        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
!
         b = x/2.d0
         d = -dlog(b)
         e = a*d
         if (abs(a) .lt. eps) then
            c = 1.d0/pi
         else
            c = a/sin(a*pi)
         endif
         if (abs(e) .lt. eps) then
            s = 1.d0
         else
            s = sinh(e)/e
         endif
         e = exp(e)
         g = e*rgamma(a,p,q)
         e = (e+1.d0/e)/2.d0
         f = 2*c*(p*e+q*s*d)
         e = a*a
         p = g*c
         q = 1.d0/g/pi
         c = a*pi/2.d0
         if (abs(c) .lt. eps) then
            r = 1.d0
         else
            r = sin(c)/c
         endif
         r = pi*c*r*r
         c = 1.d0
         d = -b*b
         ya = f+r*q
         ya1 = p
         do n = 1,maxit
            f = (f*n+p+q)/(n*n-e)
            c = c*d/n
            p = p/(n-a)
            q = q/(n+a)
            g = c*(f+r*q)
            h = c*p-n*g
            ya = ya+g
            ya1 = ya1+h
            del = abs(g)/(1.d0+abs(ya))
            del1 = abs(h)/(1.d0+abs(ya1))
            if (del+del1 .lt. eps) go to 1
         enddo
         stop 'bessjy 1'
   1     f = -ya
         g = -ya1/b
      else if (x.ge.xmin .and. x.lt.xmax) then
!
!        Temme's PQ method (besspqa) for intermediate x
!        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
!
         c = 0.25d0-a*a
         b = x+x
         p = pi
         e = (x*cos(a*pi)/pi/eps)**2
         p = 1.d0
         q = -x
         r = 1.d0+x*x
         s = r
         do n = 2,maxit
            d = (n-1+c/n)/s
            p = (2*n-p*d)/(n+1)
            q = (-b+q*d)/(n+1)
            s = p*p+q*q
            r = r*s
            if (r*n*n .gt. e) go to 2
         enddo
         stop 'bessjy 2'
   2     p = p/s
         f = p
         q = -q/s
         g = q
         do m = n,1,-1
            r = (m+1)*(2.d0-p)-2.d0
            s = b+(m+1)*q
            d = (m-1+c/m)/(r*r+s*s)
            p = d*r
            q = d*s
            e = f+1.d0
            f = p*e-g*q
            g = q*e+p*g
         enddo
         f = 1.d0+f
         d = f*f+g*g
         pa = f/d
         qa = -g/d
         d = a+0.5d0-p
         q = q+x
         pa1 = (pa*q-qa*d)/x
         qa1 = (qa*q+pa*d)/x
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      else if (x .ge. xmax) then
!
!        and Hankel's asymptotic expansions for large x
!        [ Abramowitz and Stegun, Section 9.2 ]
!
         p = 0.d0
         q = 0.d0
         do ia = 0,1
            pa = p
            qa = q
            y = 4.d0*(a+ia)**2
            z = 8.d0*x
            d = 0.d0
            w = -1.d0
            p = 1.d0
            q = 0.d0
            tp = 1.d0
            do k = 1,maxit
               d = d+z
               w = w+2.d0
               tq = +tp*(y-w*w)/d
               q = q+tq
               d = d+z
               w = w+2.d0
               tp = -tq*(y-w*w)/d
               p = p+tp
               if (abs(tp)+abs(tq) .lt. eps) go to 3
            enddo
            stop 'bessjy 3'
   3        p = p-0.5d0*tp
            q = q-0.5d0*tq
         enddo
         pa1 = p
         qa1 = q
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      endif
!
!     now recur upwards from Y(a,x) to Y(v,x),
!     scaling to avoid overflow along the way
!
      p = 0.d0
      if (na .gt. 0) then
         y = 2.d0/x
         do n = 1,na
            h = y*(a+n)*g-f
            f = g
            g = h
   4        if (abs(f) .gt. 4.d0) then
               p = p+1.d0
               f = 0.0625d0*f
               g = 0.0625d0*g
               go to 4
            endif
         enddo
      endif
      cy = f
      dy = (v/x)*f-g
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = dlog(sy)+p*dlog(16.d0)
!
!     finally, calculate J(v,x) and dJ(v,x)/dx
!
      vv = max(xmin,v)
      if (x .ge. vv) then
!
!        using upward recursion in the classically allowed region
!
         f = d*(pa*c-qa*s)
         g = d*(qa1*c+pa1*s)
         if (na .gt. 0) then
            y = 2.d0/x
            do n = 1,na
               h = y*(a+n)*g-f
               f = g
               g = h
            enddo
         endif
         cj = f
         dj = (v/x)*f-g
         sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)
      else
!
!        and CF1 in the classically forbidden region
!        [ Numerical Recipes, 2nd Edition, Section 6.7 ]
!
         ap = 1.d0
         a = v/x
         bp = 0.d0
         b = 1.d0
         f = 0.d0
         g = 0.d0
         y = 2.d0/x
         w = y/pi
         do n = 1,maxit
            an = y*(v+n)*a-ap
            ap = a
            a = an
            bn = y*(v+n)*b-bp
            bp = b
            b = bn
            if (abs(b) .gt. abs(a)) then
               ap = ap/b
               a = a/b
               bp = bp/b
               b = 1.d0
               if (abs(a-f) .lt. eps*abs(f)) then
                  cj = w/(dy-cy*a)
                  dj = a*cj
                  go to 5
               endif
               f = a
            else
               bp = bp/a
               b = b/a
               ap = ap/a
               a = 1.d0
               if (abs(b-g) .lt. eps*abs(g)) then
                  dj = w/(dy*b-cy)
                  cj = b*dj
                  go to 5
               endif
               g = b
            endif
         enddo
         stop 'bessjy 4'
   5     sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)-ey
      endif
      return
      end
