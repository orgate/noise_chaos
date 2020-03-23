      SUBROUTINE RK_CUST(y,dydx,n,h,yout,a,b)
! Y=Y_i, DYDX= Y'_i, n-dim,  h -step, YOUT=Y_i+1, derivs calculates derivs
! FOR AUTONOMOUS SYSTEM WITHOUT TIME DEPENDENCE
      INTEGER n,NMAX
      DOUBLE PRECISION h,dydx(n),y(n),yout(n),a(n,n,n),b(n,n)
      EXTERNAL derivs
      PARAMETER (NMAX=200)
      INTEGER i
      DOUBLE PRECISION h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5
      h6=h/6.
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(yt,dyt,n,a,b)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(yt,dym,n,a,b)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(yt,dyt,n,a,b)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))

14    continue
      return
      END SUBROUTINE RK_CUST
