c     pot.f
c     ----------------------------------------------------------------- 
c     Potential energy surface subroutines for CCP6 hyperspherical
c     coordinate reactive scattering program ABC. This version with 
c     tried and tested H+H2, F+H2 and Cl+H2 surfaces dated 31/03/2000.
c     ----------------------------------------------------------------- 
 
      subroutine potsub (r,vev)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine chooses which potential to use on the basis 
c     of the atomic masses in common /masses/, and also ensures
c     that the potential is called with the bond lengths in the
c     correct order.
c     ----------------------------------------------------------------- 
c
      dimension r(3),s(3),m(3)
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred 
c
      imax = 1
      imin = 1
      do i = 1,3
         m(i) = nint(mass(i))
         if (m(i) .gt. m(imax)) then
            imax = i
         else if (m(i) .lt. m(imin)) then
            imin = i
         endif 
      enddo
      if (imax .eq. imin) then
         imin = 1
         imid = 2
         imax = 3
      else
         imid = 6-imax-imin
      endif
      s(1) = r(imax)
      s(2) = r(imid)
      s(3) = r(imin)
      mmax = m(imax)
      mmid = m(imid)
      mmin = m(imin)
      if (mmin .lt. 1) stop 'potsub 1'
      if (mmid .gt. 2) stop 'potsub 2'
      if (mmax .le. 2) then
         call hh2pot (s,vev)
      else if (mmax .eq. 19) then
         call fh2pot (s,vev)
      else if (mmax .eq. 35) then
         call clh2pt (s,vev)
      else if (mmax .eq. 37) then
         call clh2pt (s,vev)
      else
         stop 'potsub 3'
      endif
      return
      end

      subroutine hh2pot (r,vev)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     LSTH H+H2 potential energy surface.
c
c     r(1) = H-H' distance in bohr
c     r(2) = H'-H" distance in bohr
c     r(3) = H"-H distance in bohr
c
c     vev = potential in eV from bottom of asymptotic H+H2 valley.
c     ----------------------------------------------------------------- 
c
      dimension r(3)
c
      call lsth (r,vev)
      return
      end
C
C     -----------------------------------------------------------------
C     LSTH H3 POTENTIAL ENERGY SURFACE
C     DIST MUST BE IN BOHR.
C     POTNTL IS THE POTENTIAL RELATIVE TO THE H2 WELL MINIMUM IN EV.
C     -----------------------------------------------------------------
C
      SUBROUTINE LSTH (DIST,POTNTL)
      implicit double precision (a-h,o-z)
      COMMON / LEPS / DE(3), BETA(3), R0(3), DELTA(3), VBASE
      COMMON / POTCOM / C6, C8, RKW(87), EKW(87), WKW(87)
      COMMON / UVCHI / U, V, CHI, R, D, RS, DS, UVX(16)
      DIMENSION DIST(3)
      EXTERNAL LSTHDATA
      DATA CONV1,CONV2,CONV3 / .52917706D-8,27.211608,1.6021892D-12 /
      SAVE NCALL
      DATA NCALL / 0 /
      IF (NCALL .EQ. 0) THEN
         write (6,61) 
         DO I = 1, 87
            RKW(I) = RKW(I)/CONV1
            EKW(I) = EKW(I)/(CONV2*CONV3)
            WKW(I) = WKW(I)/(CONV2*CONV3)*CONV1*CONV1
         ENDDO
         C6 = C6/(CONV2*CONV3*CONV1**6)
         C8 = C8/(CONV2*CONV3*CONV1**8)
         NCALL = 1
      ENDIF 
      CALL LSTHPOT ( DIST, V )
      POTNTL = (V+0.174475D0)*CONV2
      RETURN
  61  format(/1x,
     +'This calculation is using the LSTH H+H2 PES'/1x,
     +'Please cite the following papers:'/1x,
     +'1. P.Siegbahn and B.Liu,  J. Chem. Phys. 68, 2457 (1978)'/1x,
     +'2. D.G.Truhlar and C.J.Horowitz,   ibid. 68, 2466 (1978)'/1x,
     +'3. D.G.Truhlar and C.J.Horowitz,   ibid. 71, 1514 (1979)')
      END

      SUBROUTINE LSTHPOT(X,E)
      implicit double precision (a-h,o-z)
      COMMON /ECOM/ EB2,EB
      COMMON /POTCOM/ C6,C8,RKW(87),EKW(87),WKW(87)
      DIMENSION X(3),TAB(3),XQ(3),XJ(3)
      DATA XXN1,XXN2,XXN3 /.0012646477,-.0001585792,.0000079707/
      DATA XXN4,FN /-.0000001151,.0035/
      DATA A,A1,C,F /-1.514663474,-1.46,-1.2148730613,2.088442/
      DATA FB,B1,B2,B3 /.052,1.7732141742,-2.0979468223,-3.9788502170/
      DATA FB2,C1,C2 /.52,3.0231771503,-1.08935219/
      DATA D1,D2,D3 /.4908116374,-.8718696387,.1612118092/
      DATA D4,FL,XL1,XL2 /-.1273731045,.79,-13.3599568553,.9877930913/
      WN=ABS((X(1)-X(2))*(X(2)-X(3))*(X(1)-X(3)))
      COS1=(X(1)*X(1)-X(2)*X(2)-X(3)*X(3))/X(2)/X(3)
      COS2=(X(2)*X(2)-X(1)*X(1)-X(3)*X(3))/X(1)/X(3)
      COS3=(X(3)*X(3)-X(2)*X(2)-X(1)*X(1))/X(1)/X(2)
      WB=(COS1+COS2+COS3)/2.0+1.0
      DO I=1,3
         T=C*(A+X(I)*(1.+X(I)*A1))*EXP(-F*X(I))
         TAB(1) = 0.0
         IF ( X(I) .LT. 10.0 ) THEN
            CALL LSTHSPL(87,RKW,EKW,WKW,1,X(I),TAB)
         ENDIF
         XQ(I)=(TAB(1)+T)/2.0
         XJ(I)=(TAB(1)-T)/2.0
      ENDDO
      XR=X(1)+X(2)+X(3)
      XX = EXP(-FB*XR*XR)
      XX1= EXP(-FB2*XR)
      XRI=1.0/X(1)+1.0/X(2)+1.0/X(3)
      EB2=WB*(C1+XR*C2+D1*XRI)*XX1
      EQ=(X(1)-X(2))**2+(X(1)-X(3))**2+(X(3)-X(2))**2
      EB3=EQ*WB*(D4*XX+D3*XX1)
      EB=(WB**2)*(B1+WB*(B2+WB*(B3))+D2*XRI)*XX
      ELONG=WB*(XL1+XL2*XR*XR)*EXP(-FL*XR)
      XLSP=(WN**2)*(XXN1+WN*(XXN2+WN*(XXN3+WN*XXN4)))
      EN = XLSP*EXP(-FN*XR**3)
      ELP=SQRT(((XJ(1)-XJ(2))**2+(XJ(1)-XJ(3))**2+(XJ(2)-XJ(3))**2)/2.)
      ELP=XQ(1)+XQ(2)+XQ(3)-ELP
      E=(ELP+EN+EB+EB2+EB3+ELONG)
      RETURN
      END

      SUBROUTINE LSTHSPL(N,X,F,W,IJ,Y,TAB)
      implicit double precision (a-h,o-z)
      DIMENSION X(N),F(N),W(N),TAB(3)
      DATA R6 / .166666666666667D0 /
      IF (Y .LE. X(1)) THEN
         I = 1
      ELSE IF (Y .GE. X(N)) THEN
         I = N-1
      ELSE
         I = 0
         DO K = 1,N
            IF (X(K) .GT. Y) GO TO 1
            I = I+1
         ENDDO
   1     CONTINUE
      ENDIF 
      MI =(I-1)*IJ+1
      KI=MI+IJ
      FLK=X(I+1)-X(I)
      RFLK = 1.0D0/FLK
      A=(W(MI)*(X(I+1)-Y)**3+W(KI)*(Y-X(I))**3)*R6*RFLK
      B=(F(KI)*RFLK-W(KI)*FLK*R6)*(Y-X(I))
      C=(F(MI)*RFLK-FLK*W(MI)*R6)*(X(I+1)-Y)
      TAB(1)=A+B+C
      RETURN
      END

      BLOCK DATA LSTHDATA
      implicit double precision (a-h,o-z)
      COMMON /POTCOM/ C6,C8,RKW(87),EKW(87),WKW(87)
      COMMON / LEPS / DE(3), BETA(3), R0(3), DELTA(3), VBASE
      DATA DE, BETA, R0, DELTA, VBASE
     > / 3*4.74764, 3*1.04435, 3*1.4011, 4*0.D0 /
      DATA C6,C8/ .66047118D-59,  .5896753D-74/
      DATA (RKW(I),I=1,44)/
     1  .21166680D-08,  .23812515D-08,  .26458350D-08,  .29104185D-08,
     2  .31750020D-08,  .34395855D-08,  .37041690D-08,  .39687525D-08,
     3  .42333360D-08,  .47625030D-08,  .52916700D-08,  .58208370D-08,
     4  .63500040D-08,  .68791710D-08,  .71437545D-08,  .73554213D-08,
     5  .74083380D-08,  .74136297D-08,  .74141588D-08,  .74612547D-08,
     6  .76729215D-08,  .79375050D-08,  .84666720D-08,  .89958390D-08,
     7  .95250060D-08,  .10054173D-07,  .10583340D-07,  .11112507D-07,
     8  .11641674D-07,  .12170841D-07,  .12700008D-07,  .13229175D-07,
     9  .13758342D-07,  .14287509D-07,  .14816676D-07,  .15345843D-07,
     9  .15875010D-07,  .16404177D-07,  .16933344D-07,  .17462511D-07,
     9  .17991678D-07,  .18520845D-07,  .19050012D-07,  .19579179D-07/
      DATA (RKW(I),I=45,87)/
     1  .20108346D-07,  .20637513D-07,  .21166680D-07,  .21695847D-07,
     2  .22225014D-07,  .22754181D-07,  .23283348D-07,  .23812515D-07,
     3  .24341682D-07,  .24870849D-07,  .25400016D-07,  .25929183D-07,
     4  .26458350D-07,  .26987517D-07,  .27516684D-07,  .28045851D-07,
     5  .28575018D-07,  .29104185D-07,  .29633352D-07,  .30162519D-07,
     6  .30691686D-07,  .31220853D-07,  .31750020D-07,  .32279187D-07,
     7  .32808354D-07,  .33337521D-07,  .33866688D-07,  .34395855D-07,
     8  .34925022D-07,  .35454189D-07,  .35983356D-07,  .36512523D-07,
     9  .37041690D-07,  .38100024D-07,  .39158358D-07,  .40216692D-07,
     9  .41275026D-07,  .42333360D-07,  .43656277D-07,  .44979195D-07,
     9  .47625030D-07,  .50270865D-07,  .52916700D-07/
      DATA (EKW(I),I=1,44)/
     1  .38355782D-10,  .28297040D-10,  .20637246D-10,  .14701904D-10,
     2  .10043069D-10,  .63492841D-11,  .33993629D-11,  .10316727D-11,
     3 -.87434721D-12, -.36464790D-11, -.54294007D-11, -.65418745D-11,
     4 -.71904983D-11, -.75136199D-11, -.75841062D-11, -.76054247D-11,
     5 -.76064144D-11, -.76064231D-11, -.76064231D-11, -.76057822D-11,
     6 -.75881650D-11, -.75357581D-11, -.73494368D-11, -.70825018D-11,
     7 -.67603262D-11, -.64020791D-11, -.60219903D-11, -.56307147D-11,
     8 -.52369150D-11, -.48466944D-11, -.44648008D-11, -.40948699D-11,
     9 -.37397225D-11, -.34012157D-11, -.30809408D-11, -.27798701D-11,
     9 -.24985746D-11, -.22372897D-11, -.19959718D-11, -.17743769D-11,
     9 -.15719817D-11, -.13881105D-11, -.12218784D-11, -.10723697D-11/
      DATA (EKW(I),I=45,87)/
     1 -.93861231D-12, -.81946400D-12, -.71362123D-12, -.62014234D-12,
     2 -.53785022D-12, -.46565062D-12, -.40240566D-12, -.34738295D-12,
     3 -.29951872D-12, -.25799337D-12, -.22200474D-12, -.19087272D-12,
     4 -.16403492D-12, -.14085484D-12, -.12093576D-12, -.10375887D-12,
     5 -.89036444D-13, -.76384837D-13, -.65525033D-13, -.56234691D-13,
     6 -.48256593D-13, -.41407636D-13, -.35530873D-13, -.30526033D-13,
     7 -.26288486D-13, -.22504339D-13, -.19470046D-13, -.16845557D-13,
     8 -.14508803D-13, -.12669045D-13, -.10750814D-13, -.93906135D-14,
     9 -.82353152D-14, -.62516898D-14, -.47345433D-14, -.37841469D-14,
     9 -.29732583D-14, -.23018774D-14, -.17612850D-14, -.13689195D-14,
     9 -.80652901D-15, -.52751357D-15, -.39672508D-15/
      DATA (WKW(I),I=1,44)/
     1  .47955853D+08,  .33383238D+08,  .24122353D+08,  .17928049D+08,
     2  .13573565D+08,  .10490948D+08,  .82183730D+07,  .65379276D+07,
     3  .51991523D+07,  .34159125D+07,  .23332144D+07,  .16170628D+07,
     4  .11375504D+07,  .80733801D+06,  .68676860D+06,  .59701325D+06,
     5  .58220055D+06,  .55824440D+06,  .57950581D+06,  .55530308D+06,
     6  .48691419D+06,  .40716760D+06,  .28151647D+06,  .19409292D+06,
     7  .12576620D+06,  .75753358D+05,  .39226713D+05,  .70412432D+04,
     8 -.13304725D+05, -.30515597D+05, -.43054465D+05, -.53595765D+05,
     9 -.59330838D+05, -.65643799D+05, -.68753151D+05, -.70834202D+05,
     9 -.71637937D+05, -.71386325D+05, -.70654894D+05, -.68601021D+05,
     9 -.66338213D+05, -.62964080D+05, -.59760288D+05, -.56332558D+05/
      DATA (WKW(I),I=45,87)/
     1 -.52415850D+05, -.47035830D+05, -.44541699D+05, -.39720725D+05,
     2 -.36276854D+05, -.31426288D+05, -.29891253D+05, -.25188332D+05,
     3 -.22741928D+05, -.19668552D+05, -.17220204D+05, -.15514319D+05,
     4 -.12735742D+05, -.11917426D+05, -.94685547D+04, -.89660376D+04,
     5 -.72596243D+04, -.63673240D+04, -.56644142D+04, -.46042174D+04,
     6 -.40364625D+04, -.34442724D+04, -.30178681D+04, -.31671433D+04,
     7 -.75450027D+03, -.35299574D+04, -.11929538D+04, -.47918458D+03,
     8 -.30556609D+04,  .20525819D+04, -.34732069D+04, -.11680244D+03,
     9 -.45006195D+03, -.34293366D+03, -.67703965D+03,  .15122975D+02,
     9 -.13076777D+03, -.23936742D+03, -.57485345D+02, -.38865755D+02,
     9 -.45225520D+02, -.23109707D+02,  .10620709D+02/
      END

      subroutine fh2pot (r,vev)
      implicit double precision (a-h,o-z)
      integer b,c,d
c
c     ----------------------------------------------------------------- 
c     Stark-Werner F+H2 potential energy surface
c
c     r(1) = H-H distance in bohr
c     r(2) = F-H distance in bohr
c     r(3) = H-F distance in bohr
c
c     vev = potential in eV from bottom of asymptotic F+H2 valley
c     ----------------------------------------------------------------- 
c
      dimension r(3)
      parameter (nmx=200)
      dimension a(nmx),b(nmx),c(nmx),d(nmx),p(nmx)
      save icall,a,b,c,d,p,ma,np
      data icall/0/
c
      if (icall .eq. 0) then
         write (6,61)
         open (unit=1,file='SW.3p',status='old')
         i=1
   2     read (1,*,end=3) nparm,b(i),c(i),d(i),a(i)
            i=i+1
            goto 2
   3     ma=i-1
         close (unit=1)
         open (unit=1,file='SW.2p',status='old')
         i=1
   4     read (1,*,end=5) p(i)
            i=i+1
            goto 4
   5     np=i-1
         close (unit=1)
         icall = 1
      endif
c
c     (Modified) Aguado-Paniagua type functions:
c
      x = min(r(2),r(3))
      y = r(1)
      z = max(r(2),r(3))
      b1 = a(ma-5)
      b2 = a(ma-4)
      b3 = a(ma-3)
      x0 = a(ma-2)
      y0 = a(ma-1)
      z0 = a(ma)
      fit = 0.0d0
      do i = 1,ma-6
         expon = b(i)*b1*(x-x0)+c(i)*b2*(y-y0)+d(i)*b3*(z-z0)
         fex = exp(-expon)
         fxy = ((x**b(i))*(y**c(i))*(z**d(i)))*fex
         fit = fit+a(i)*fxy
      enddo
      xr=x-p(3)
      yr=y-p(9)
      zr=z-p(3)
      fx=dexp(-p(2)*xr)
      fy=dexp(-p(8)*yr)
      fz=dexp(-p(2)*zr)
      xval=-p(1)*(1.0d0+p(2)*xr+p(4)*xr**2+p(5)*xr**3)*fx+p(6)
      yval=-p(7)*(1.0d0+p(8)*yr+p(10)*yr**2+p(11)*yr**3)*fy+p(12)
      zval=-p(1)*(1.0d0+p(2)*zr+p(4)*zr**2+p(5)*zr**3)*fz+p(6)
      vagpan=fit+xval+yval+zval
      vev = 27.2113961d0*vagpan
      return
  61  format(/1x,
     +'This calculation is using the SW F+H2 PES'/1x,
     +'Please cite: ',
     +'K.Stark and H-J.Werner, J. Chem. Phys. 104, 6515 (1996)')
      end

      subroutine clh2pt (r,vev)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Bian-Werner Cl+HD potential energy surface
c
c     r(1) = H-D distance in bohr
c     r(2) = D-Cl distance in bohr
c     r(3) = H-Cl distance in bohr
c
c     vev = potential in eV from bottom of Cl+HD valley
c     ----------------------------------------------------------------- 
c
      parameter (VMAX = 50.d0)
      dimension r(3)
c
      x = r(2)
      y = r(1)
      z = r(3)
      call bw4pot (x,y,z,vau)
      vev = 27.21139610d0*vau
      vev = min(vev,VMAX)
      return
      end

      subroutine bw4pot (x,y,z,v)
      implicit double precision (a-h,o-z)
c
c-------------------------------------------------------------------
c
c     System:   ClH2
c     Name:     BW4
c     Author:  Wensheng Bian and Joachim Werner 
c     Functional form: Aguado-Paniagua 
c     Energy Zero Point: the asymptote Cl+H2(re) in a.u.
c
c     This subroutine calculates the potential energy for the
c     NON spin-orbit corrected MRCI+Q surface for the system
c     **Scale fact=.948d0**
c
c            Cl + H2
c
c     Input are the three distances x,y,z 
c         Cl    H1    H2
c         |_____|
c            x 
c               |_____|
c                  y
c         |___________|
c               z 
c
c-------------------------------------------------------------------
c
      parameter(nmx=600,np=27,mma=20)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      dimension xex(0:mma),xey(0:mma),xez(0:mma)
      dimension p(np) 
      save ifirst,p
      data ifirst/-1/
      data p/14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021,
     1       1.079946167,-0.02206944094,-7.109456997,
     1        36.79845478,-109.3716794,176.4925683,
     1        -120.4407534,2.351569228,1.082968302,
     1       14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021/
c
c on first call of this subroutine, read in three-body parameter
c
      if (ifirst.eq.-1) then
         write (6,61)
         call bw4ini
         ifirst=0
      end if
c 
c.... Three-Body-Potential : Aguado-Paniagua
c
c.... initialize the non-linear parameters
c
      b1 = a(ma-2)
      b2 = a(ma-1)
      b3 = a(ma)
      xexpon = b1*x
      yexpon = b2*y
      zexpon = b3*z
      exponx=dexp(-xexpon)
      expony=dexp(-yexpon)
      exponz=dexp(-zexpon)
      fex = x*exponx
      fey = y*expony
      fez = z*exponz
      xex(0)=1
      xey(0)=1
      xez(0)=1
      xex(1)=fex
      xey(1)=fey
      xez(1)=fez
      do m=2,mmax-1
         xex(m)=xex(m-1)*fex
         xey(m)=xey(m-1)*fey
         xez(m)=xez(m-1)*fez
      enddo
      fit = 0.0d0
      do i=1,ma-3
         fit=fit+xex(ix(i))*xey(iy(i))*xez(iz(i))*a(i)
      enddo
c
c.... Two-Body-Potential : Aguado-Paniagua
c
c       c0      c1      c2      c3      c4     c5    c6     alpha  beta
c  ClH  p(1)    p(2)    p(3)    p(4)    p(5)   p(6)  p(7)   p(8)   p(9)
c  HH   p(10)   p(11)   p(12)   p(13)   p(14)  p(15) p(16)  p(17)  p(18)
c  ClH  p(19)   p(20)   p(21)   p(22)   p(23)  p(24) p(25)  p(26)  p(27)
c
      rhox=x*dexp(-p(9)*x)
      rhoy=y*dexp(-p(18)*y)
      rhoz=z*dexp(-p(27)*z)
      xval=p(1)*dexp(-p(8)*x)/x
      yval=p(10)*dexp(-p(17)*y)/y
      zval=p(19)*dexp(-p(26)*z)/z
      do i=1,6
         xval=xval+p(i+1)*rhox**i
         yval=yval+p(i+10)*rhoy**i
         zval=zval+p(i+19)*rhoz**i
      enddo
c
c.... Total Potential in atomic units
c
      v=fit+xval+yval+zval
c
c.... Relative to Cl+H2(re)
c
      v = v+0.1747737310d0
c
c     if (x.lt.1.75d0.or.z.lt.1.75d0.or.y.lt.0.8d0) then
c        v=1.0d0
c     end if
c
      return
  61  format(/1x,
     +'This calculation is using the BW Cl+H2 PES'/1x,
     +'Please cite: ',
     +'W.Bian and H-J.Werner, J. Chem. Phys. 112, 220 (2000)')
      end

c------------------------------------------------------------------------
      subroutine bw4ini 
      implicit double precision (a-h,o-z)
      parameter(nmx=600)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      open(1,file='BW.3p',status='old')
      i=1
      rewind 1
10    read(1,*,end=100) nparm,ix(i),iy(i),iz(i),a(i)
      i=i+1
          goto 10
100   ma=i-1
      close(1)
      m=ix(ma-3)
      mmax=m+1
      return
      end
