c     fun.f
c     ----------------------------------------------------------------- 
c     Special function subroutines for CCP6 hyperspherical 
c     coordinate reactive scattering program ABC.
c     This version with translations of Temme's Bessel 
c     function routines dated 31 March 2000.
c     ----------------------------------------------------------------- 
c
      function rotmel (jtot,ipar,k,m,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses redrot to calculate a parity adapted reduced 
c     rotation matrix element with angular momentum jtot, parity ipar,
c     angular momentum projections k and m, and argument beta.
c     ----------------------------------------------------------------- 
c     
      rotmel = 0.d0
      if (k.gt.jtot .or. k.lt.0) return
      if (m.gt.jtot .or. m.lt.0) return
      sign = ipar*(-1)**(jtot+k)
      rj = jtot
      rk = k
      rm = m
      rotmel = redrot(rj,rk,rm,beta)+sign*redrot(rj,-rk,rm,beta)
      if (k .eq. 0) rotmel = sqrt(0.5d0)*rotmel
      if (m .eq. 0) rotmel = sqrt(0.5d0)*rotmel
      return
      end

      function redrot (rj,rk,rm,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses eq. (4.1.23) of Edmonds 
c     to calculate the reduced rotation matrix element 
c     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
c     ----------------------------------------------------------------- 
c     
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
c     half integer angular momenta
c
      sj = half*nint(two*rj)
      sk = half*nint(two*rk)
      sm = half*nint(two*rm)
c
c     projection ranges
c
      redrot = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return      
c
c     reflection symmetries
c      
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
c
c     evaluation
c
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      redrot = sign(d4,d3)
      if (mod(isign,2) .ne. 0) redrot = -redrot
      return
      end

      function pjacob (n,a,b,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Jacobi polynomial p(n,a,b;x)
c     Abramowitz and Stegun eq. (22.7.1)
c     ----------------------------------------------------------------- 
c
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two 
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob = fp
      return
      end

      subroutine rbesjy (ell,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Riccati-Bessel functions of fractional order  
c     and their first derivatives with respect to x:
c
c     j(ell,x) = cj * exp(ej)
c     y(ell,x) = cy * exp(ey)
c     d/dx j(ell,x) = dj * exp(ej)
c     d/dx y(ell,x) = dy * exp(ey)
c     ----------------------------------------------------------------- 
c
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbesjy 0'       
      v = ell+0.5d0
      call bessjy (v,x,cj,dj,ej,cy,dy,ey)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dj = dj+cj/(2.d0*x)
      sj = sqrt(cj*cj+dj*dj)
      cj = cj/sj
      dj = dj/sj
      ej = ej+dlog(sj)+ex
      dy = dy+cy/(2.d0*x)
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = ey+dlog(sy)+ex
      return
      end

      subroutine rbessk (ell,x,ck,dk,ek)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Modified Riccati-Bessel function of the third kind
c     and its first derivative with respect to x:
c
c     k(ell,x) = ck * exp(ek)
c     d/dx k(ell,x) = dk * exp(ek)
c     ----------------------------------------------------------------- 
c
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbessk 0'       
      v = ell+0.5d0
      call mbessk (v,x,ck,dk,ek)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dk = dk+ck/(2.d0*x)
      sk = sqrt(ck*ck+dk*dk)
      ck = ck/sk
      dk = dk/sk
      ek = ek+dlog(sk)+ex
      return
      end

      subroutine bessjy (v,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses a combination of methods (mostly due
c     to Temme) to calculate the Ordinary Bessel functions
c
c     J(v,x) = cj * exp(ej)
c     Y(v,x) = cy * exp(ey)
c
c     and their first derivatives with respect to x
c
c     d/dx J(v,x) = dj * exp(ej)
c     d/dx Y(v,x) = dy * exp(ey)
c
c     for a given real order v >= 0 and real argument x > 0.  
c     Note the exponential scaling, which is used to avoid
c     overflow of Y(v,x) and underflow of J(v,x) for v >> x.
c     ----------------------------------------------------------------- 
c
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
c
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'bessjy 0'
      pi = acos(-1.d0)
      xmin = 3.d0
      xmax = 5.d0-dlog10(eps)
c
c     begin by calculating Y(a,x) and Y(a+1,x) for |a| <= 1/2
c
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
c
c        using Temme's series (bessya) for small x 
c        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
c
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
c
c        Temme's PQ method (besspqa) for intermediate x  
c        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
c
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
c 
c        and Hankel's asymptotic expansions for large x           
c        [ Abramowitz and Stegun, Section 9.2 ]
c
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
c
c     now recur upwards from Y(a,x) to Y(v,x),
c     scaling to avoid overflow along the way
c
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
c
c     finally, calculate J(v,x) and dJ(v,x)/dx
c
      vv = max(xmin,v)
      if (x .ge. vv) then
c
c        using upward recursion in the classically allowed region
c
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
c
c        and CF1 in the classically forbidden region
c        [ Numerical Recipes, 2nd Edition, Section 6.7 ]
c
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

      subroutine mbessk (v,x,ck,dk,ek)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses Temme's method [ N.M.Temme, J Comput Phys
c     19 (1975) 324-337 ] to calculate the Modified Bessel function
c
c     K(v,x) = ck * exp(ek)
c
c     and its first derivative with respect to x
c
c     d/dx K(v,x) = dk * exp(ek)
c
c     for a given real order v >= 0 and real argument x > 0.  
c     Note the exponential scaling, which is used to avoid
c     overflow of K(v,x) for v >> x and underflow for v << x.
c     ----------------------------------------------------------------- 
c
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
c
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'mbessk 0'
      pi = acos(-1.d0)
      xmin = 1.d0
c
c     begin by calculating K(a,x) and K(a+1,x) for |a| <= 1/2
c
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
c
c        using Temme's series for small x 
c
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
c
c        and Temme's PQ method for large x  
c
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
c
c     now recur upwards from K(a,x) to K(v,x),
c     scaling to avoid overflow along the way
c
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

      function rgamma(x,odd,even)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Direct fortran translation of Temme's algol routine for computing
c     rgamma = 1/Gamma(1-x), along with its odd and even parts, for
c     abs(x) .le. 0.5. [ N.M.Temme, J Comput Phys 19 (1975) 324-337 ]
c     ----------------------------------------------------------------- 
c
      dimension b(12)
      data b / -0.283876542276024d0, -0.076852840844786d0,
     *         +0.001706305071096d0, +0.001271927136655d0,
     *         +0.000076309597586d0, -0.000004971736704d0,
     *         -0.000000865920800d0, -0.000000033126120d0,
     *         +0.000000001745136d0, +0.000000000242310d0,
     *         +0.000000000009161d0, -0.000000000000170d0 / 
      save b
c
      x2 = x*x*8.d0
      alfa = -0.000000000000001d0 
      beta = 0.d0
      do i = 12,2,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      even = (beta/2.d0+alfa)*x2-alfa+0.921870293650453d0
      alfa = -0.000000000000034d0
      beta = 0.d0
      do i = 11,1,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      odd = 2*(alfa+beta)
      rgamma = odd*x+even
      return
      end      
