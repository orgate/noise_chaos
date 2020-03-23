  PROGRAM chaos !simulation of chaotic attractors
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: Y(:),DY(:),YN(:) 
  DOUBLE PRECISION, ALLOCATABLE :: X(:),DX(:),XN(:) 
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:,:),B(:,:),AL(:,:),LEXPMAV(:),DME(:) 
  DOUBLE PRECISION, ALLOCATABLE :: AS(:,:,:),BS(:,:),XS(:) 
  DOUBLE PRECISION, ALLOCATABLE :: X1(:),X2(:),Y1(:),Y2(:),WR(:),WI(:),Z1(:),Z2(:) 
  DOUBLE PRECISION, ALLOCATABLE :: YM1(:),YM2(:),ZM1(:),ZM2(:),C1(:),C2(:),H(:), HX(:) 
  INTEGER N,K,L,IW,K1,K2,K3,CC,NP,MCH,MP,MF,MC,K4,MMCH,MMP,MMF,MATKEY,IC,KC,NH,NX
  INTEGER, ALLOCATABLE :: indx(:)
  DOUBLE PRECISION T, DT, TFIN, TMEAS, TM, DIF,ran3,CH,LX,LT,LEXP,TEXP,EXPAV, NORM,LH,LHMAX,XMAX
  DOUBLE PRECISION GS,XG,GW,XP,DEXPAV,NLL,ML,DL,LL,WF,WV,AU,TR,TRR,EMA,EMM,EMAX,EMIN,WX
  GS(XP,GW)=DEXP(-XP*XP/(GW*GW*2))! for picking random coefficents for each dimension
  N=20   !SPACE DIMENSIONALITY
  TFIN=25.D0/N/N !d-4 for n=200 RELAXATION TIME 
  TEXP=2000.D0/N/N ! 400 for all stat 5000 for nice pitct! L-exp meas. time
  DT=1.d-1/N/N ! 2.d-6 for n=200 ! INTEGRATION TIME STEP
  LT=1.D-3 ! TIME EVOLUTION BEFORE RESCALING FOR LYAPUNOV EXPONENT
  LX=1.D-3 ! DIFFERENCE IN INITIAL CONDITIONS FOR LYAPUNOV EXPONENT
  GW=1.D0/DSQRT(1.D0) ! WIDTH OF A GAUSSIAN DISTRIBUTION 
  NP=1000000 ! MAXIMUM NUMBER OF POINTS IN THE OUTPUT.
  CH=1.D-1 ! CHAOS CONSTANT, IF IN ONE UNIT TRAJECTORIES DIVERGE MORE --> CHAOS
  TR=TFIN*0.9/NP ! INTERVAL BETWEEN RECORDING POINTS
  MC=200! mc/ic -- NUMBER OF REPLICAS OF COEFFICIENTS
  IC=4! NUMBER OF REPLICAS IN INITIAL CONDITIONS
  NH=15 ! NUMBER OF POINTS IN LYAPUNOV EXPONENT HISTOGRAM
  LHMAX=0.3*N*N ! MAXIMUM OF HISTOGRAM FOR LH
  MATKEY=1 ! MATKEY =1 -- NO MATRIX DIAGONALIZATION, ANY ELSE, WITH
   NX=150 ! NUMBER OF POINTS IN  X HISTOGRAM
  XMAX=1.5D0*N ! MAXIMUM OF HISTOGRAM 
  ALLOCATE(Y(N),DY(N),YN(N),LEXPMAV(N),DME(N))
  ALLOCATE(X(N),DX(N),XN(N),WR(N),WI(N),H(-NH:NH),HX(-NX:NX))
  ALLOCATE(A(N,N,N),B(N,N),AL(N,N),indx(N))
  ALLOCATE(AS(N,N,N),BS(N,N),XS(N))
  ALLOCATE(X1(NP),X2(NP),Y1(NP),Y2(NP),Z1(NP),Z2(NP))
  ALLOCATE(YM1(NP),YM2(NP),ZM1(NP),ZM2(NP),C1(NP),C2(NP))
  open(unit=12,file='rand.dat.txt')
  read(12,*)iw
  close(12)
  MCH=0 ! NUMBER OF CHAOTIC TRAJECTORIES
  MP=0  ! NUMBER OF PERIODIC TRAJECTORIES
  MF=0  ! NUMBER OF FIXED POINT TRAJECTORIES
  MMCH=0 ! NUMBER OF CHAOTIC TRAJECTORIES FOR EIGEN VALUE ignore this fornow
  MMP=0  ! NUMBER OF PERIODIC TRAJECTORIES FOR EIGEN VALUE ignore this for now
  MMF=0  ! NUMBER OF FIXED POINT TRAJECTORIES FOR EIGEN VALUE ignore this for now...
  NLL=0  ! NUMBER OF POSITIVE EIGENVALUES
  EXPAV=0.D0
  DEXPAV=0.D0
  ML=0.D0
  DL=0.D0
  EMA=0.D0 ! MAXIMUM EXPONENT
  EMM=1.D+2 ! MINIMUM EXPONENT
  EMAX=0.D0 ! MAXIMUM DIRECT EXPONENT
  EMIN=1.D+4! MINIMUM DIRECT EXPONENT
       DO CC=1,MC/IC ! CYCLE ON REPLICAS NOW EACH REPLICA HAS A NUMBER OF IC
  A=0.d0
  B=0.d0
!  iw=-5462456 ! for random number generator
  DO K1=1,N ! RANDOM MATRIX OF COEFFICIENTS
 !    XP=RAN3(IW)*2-1.D0 ! UNIFORM ASSIGNMENT OF INITIAL CONDITIONS
     DO K2=1,N
        DO K3=1,N
!           A(K1,K2,K3)=RAN3(IW)*2-1.D0
54     XP=(RAN3(IW)*2-1.D0)*4*GW ! GAUSSIAN ASSIGNMENT
       IF(GS(XP,GW).LE.RAN3(IW)) GO TO 54
         A(K1,K2,K3)=XP
        END DO
!        B(K1,K2)=RAN3(IW)*2-1.D0
64     XP=(RAN3(IW)*2-1.D0)*4*GW ! GAUSSIAN ASSIGNMENT
       IF(GS(XP,GW).LE.RAN3(IW)) GO TO 64
       B(K1,K2)=XP   
    END DO
  END DO
  WRITE(*,*) 'NEW COEFFICIENTS'
!  b=0.d0 !!!!!!!!!!!!!!!COMMENT THIS OUT!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO KC=1,IC
  DO K1=1,N ! RANDOM INITIAL CONDITIONS
 !    XP=RAN3(IW)*2-1.D0 ! UNIFORM ASSIGNMENT OF INITIAL CONDITIONS
44     XP=(RAN3(IW)*2-1.D0)*4*GW ! GAUSSIAN ASSIGNMENT
       IF(GS(XP,GW).LE.RAN3(IW)) GO TO 44
       X(K1)=XP*N ! THE SIZE OF ATTRACTOR SCALES WITH N
       END DO
  K4=0
  T=0.D0
  H=0.D0
  DO WHILE (T.LE.TFIN)!.AND.DIF.LE.EPS) !evolving the trajectory till initial transients are dealth with and traj settles on an attractor
     CALL DERIVS (X,DX,N,A,B)
     CALL RK_CUST(X,DX,N,DT,XN,A,B)
     T=T+DT
     X=XN 
  END DO

  Y=X + LX/DSQRT(1.D0*N)! i guess the displaced trajectory (Y, displaced by a small amount from X) in order to calculate lyapunov exponent
  T=0.D0
  LEXP=0.D0
  LEXPMAV=0.D0
  DME=0.D0
  WF=0.D0
  WV=0.D0
  WX=0.D0
  TMEAS=LT
  DO WHILE (T.LE.TEXP) ! LOOP ON MEASURING L-EXPONENT 
  DO WHILE (T.LE.TMEAS)! LOOP ON A SINGLE STRETCHING
     CALL DERIVS (X,DX,N,A,B)
     CALL RK_CUST(X,DX,N,DT,XN,A,B)
     CALL DERIVS (Y,DY,N,A,B)
     CALL RK_CUST(Y,DY,N,DT,YN,A,B)
     T=T+DT
     X=XN
     Y=YN
     IF(K4.LE.NP-1) THEN ! KEEPINT THE RECORD OF THE FIRST TWO POINTS
     K4=K4+1
        X1(K4)=X(1)! for plotting the trajectory
        X2(K4)=X(2) ! for plotting the trajectory..both X1 and X2 help us in plotting the 2d cross section of the trajecotry.
      END IF
   END DO

!!$   IF(MATKEY.EQ.1) GO TO 77
!!$   AL=0.D0 ! LINEARIZATION AND MATRIX
!!$   DO K1=1,N
!!$      AL(K1,K1)=-3*X(K1)*X(K1)
!!$      DO K2=1,N
!!$         AL(K1,K2)=AL(K1,K2)+B(K1,K2)
!!$         DO K3=1,N
!!$            AL(K1,K2)=AL(K1,K2)+(A(K1,K2,K3)+A(K1,K3,K2))*X(K3)
!!$         END DO
!!$      END DO
!!$   END DO
!!$  CALL balanc(al,n,n)
!!$  CALL elmhes(al,n,n)
!!$  CALL hqr(al,n,n,wr,wi)
!!$  CALL indexx(n,wr,indx)
!!$     DO K=1,N
!!$        AU=WR(INDX(N-K+1))
!!$        LEXPMAV(K)=LEXPMAV(K)+AU
!!$        DME(K)=DME(K)+AU*AU
!!$     END DO

! END MARTIX PART  

77   DIF=0.D0
     DO K=1,N
        DIF=DIF+(Y(K)-X(K))**2
        WX=WX+X(K)**2
     END DO
     DIF=DSQRT(DIF)/LX
     AU=DLOG(DIF)
     LH=AU/LT
     WV=WV+AU*AU
     LEXP=LEXP + AU
     K=NH*LH/LHMAX ! HISTOGRAM OF EXPONENT
     IF(dabs(1.d0*K).LE.NH) H(K)=H(K)+1   ! some lyapunov histogram

       DO K=1,N
        K2=(NX*X(K))/XMAX
        IF (IABS(K2).LE.NX) HX(K2)=HX(K2)+1 !HISTOGRAM OF THE X. 
     END DO


     DO K=1,N
        Y(K)=X(K)+(Y(K)-X(K))/DIF
     END DO
     TMEAS=TMEAS+LT ! this is becuase at first iterations T goes till TMEAS. at the seocnd iteration (loop on a single stretching) T is already TMEAS. Therefor,
     !we need T going from TMEAS (old value) to TMEAS +LT....the second stretch.
  END DO ! end of the loop on sampling points
  LEXP=LEXP/T          ! AVERAGE EXPONENT FOR THE DIRECT METHOD
  WV=WV/T/LT           ! DISPERSION FOR DIRECT METHOD
  WV=DSQRT(DABS(WV-LEXP*LEXP))
  WX=DSQRT(WX/N/T*LT)
!!$  LEXPMAV=LEXPMAV/T*LT ! AVERAGE EXPONENT FOR THE MATRIX METHOD
!!$  DME=DME/T*LT           ! DISPERSION FOR THE MATRIX METHOD
!!$  LL=LEXPMAV(1) ! LARGEST MATRIX LYAPUNOV EXPONENT 
!!$  DO K=1,N
!!$  DME(K)=DSQRT(DME(K)-LEXPMAV(K)**2)
!!$  END DO
!!$  WF=DME(1)
  
     IF(LEXP.GE.CH) THEN ! FOR TRADITIONAL EXPONENTS
        MCH=MCH+1
        EXPAV=EXPAV+LEXP
        DEXPAV=DEXPAV+LEXP*LEXP
        IF(LEXP.GE.EMAX) THEN 
           EMAX=LEXP
        y1=x1  
        y2=x2
        AS=A ! WRITING THE COEFFICIENTS AND INITIAL CONDITIONS  
        BS=B ! FOR THE LARGEST CHAOS
        XS=X
       END IF

!!$        K=1
!!$        DO WHILE (LEXPMAV(K).GE.CH) ! WE COUNT THE NUMBER OF POSITIVE MATRIX EXPONENTS
!!$        K=K+1  ! ONLY FOR DIRECTLY DETERMINED CHAOS
!!$        NLL=NLL+1.D0
!!$        END DO
     ELSE IF(LEXP.LT.CH.AND.LEXP.GT.-CH) THEN
        MP=MP+1
           c1=x1
           c2=x2
     ELSE IF(LEXP.LE.-CH) THEN
      IF(LEXP.LE.EMIN) THEN
           EMIN=LEXP
           Z1=X1
           Z2=X2
        END IF
        MF=MF+1
     END IF
 OPEN(UNIT=90, FILE='hist.dat.txt') ! writing the current trajectory
  DO K1=-NH,NH
     WRITE(90,*)K1*LHMAX/NH,H(K1)
  END DO
  CLOSE(90)
  OPEN(UNIT=110, FILE='histx3.dat.txt') ! writing the current trajectory
  DO K1=-NX,NX
     WRITE(110,*)K1*XMAX/NX,HX(K1)
  END DO
  CLOSE(110)


  OPEN(UNIT=90, FILE='xy.dat.txt') ! writing the current trajectory
  DO K1=1,K4
     WRITE(90,*)x1(K1),x2(K1)
  END DO
  CLOSE(90)
  OPEN(UNIT=10, FILE='xymax_st.dat.txt')
  DO K1=1,K4
     WRITE(10,*)Y1(K1),Y2(K1)
  END DO
  CLOSE(10)
  open(unit=10,file='coeff.dat.txt') ! writing largest chaos  coefficients
  do k1=1,n
     do k2=1,n
        do k3=1,n
           write (10,*) as(k1,k2,k3)
        end do
        write (10,*) bs(k1,k2)
     end do
     write (10,*) xs(k1)
  end do 
  close(10)

 OPEN(UNIT=15, FILE='xymin_st.dat.txt') ! printing out the fixed point as Lyapunmov expoenent is less than -1
  DO K1=1,K4
     WRITE(15,*)Z1(K1),Z2(K1)
  END DO
  CLOSE(15)
 OPEN(UNIT=13, FILE='xycycle_st.dat.txt') ! printing out the cycle or quasiperiodic/periodic  orbits as lyapunov exponent is between -1 and 1
  DO K1=1,K4
     WRITE(13,*)C1(K1),C2(K1)
  END DO
  CLOSE(13)

!!$    IF(LL.GE.CH) THEN ! FOR MATRIX EXPONENTS
!!$        MMCH=MMCH+1
!!$        ML=ML+LL
!!$        DL=DL+LL*LL
!!$        IF(LL.GE.EMA) THEN 
!!$           EMA=LL
!!$        ym1=x1
!!$        ym2=x2
!!$       END IF
!!$      IF(LL.LE.EMM) THEN
!!$           EMM=LL
!!$           ZM1=X1
!!$           ZM2=X2
!!$        END IF
  OPEN(UNIT=25,FILE='statexp_M.dat.txt',position='append')
  write(25,*) N,LEXP,WV,WX,MCH,MP,MF,CC 
!!$  write(25,*) LL,WF,MMCH,MMP,MMF,CC
!!$  WRITE(25,22) (LEXPMAV(K),DME(K),K=1,MIN0(N,8))
  CLOSE(25)
!!$     ELSE IF(LL.LT.CH.AND.LL.GT.-CH) THEN
!!$        MMP=MMP+1
!!$     ELSE IF(LL.LE.-1.D-1) THEN
!!$        MMF=MMF+1
!!$     END IF
     WRITE(*,*)LEXP,WV,WX,MCH,MP,MF,CC ! DATA FROM A PARTICULAR REALIZATION
!!$     WRITE(*,*) LL,WF,MMCH,MMP,MMF,CC
!!$     WRITE(*,*)(LEXPMAV(K),K=1,4)
     write(*,*)'lexpmax, lexpmin', emax, emin
  WRITE(*,*)

22 FORMAT(20E12.4)
END DO ! END OF INITICAL CONDITION LOOP
  END DO ! END OF CONFIGURATION LOOP


  EXPAV=EXPAV/MCH  !FINAL OUTPUT
  DEXPAV=DEXPAV/MCH
  DEXPAV=DSQRT(DEXPAV-EXPAV*EXPAV)

  ML=ML/MMCH
  DL=DL/MMCH
  DL=DSQRT(DL-ML*ML)
  NLL=NLL/MCH ! AS IT IS COUNTED ONLY FOR THE DIRECTLY DETERMINED CHAOS

  OPEN(UNIT=30, FILE='statd.dat.txt',position='append')
!!$  OPEN(UNIT=31, FILE='statm.dat.txt',position='append')
  WRITE (30,*)N,EXPAV,DEXPAV,MCH,MP,MF
!!$  WRITE (31,*)N,ML,DL,NLL,MMCH,MMP,MMF
  CLOSE(30)
!!$  CLOSE(31)

!!$  OPEN(UNIT=10, FILE='xymax_st.dat.txt')
!!$  DO K1=1,K4
!!$     WRITE(10,*)Y1(K1),Y2(K1)
!!$  END DO
!!$  CLOSE(10)
!!$ OPEN(UNIT=15, FILE='xymin_st.dat.txt')
!!$  DO K1=1,K4
!!$     WRITE(15,*)Z1(K1),Z2(K1)
!!$  END DO
!!$  CLOSE(15)
!!$ OPEN(UNIT=13, FILE='xycycle_st.dat.txt')
!!$  DO K1=1,K4
!!$     WRITE(13,*)C1(K1),C2(K1)
!!$  END DO
!!$  CLOSE(13)


!!$  OPEN(UNIT=25, FILE='xymax_ma.dat.txt')
!!$  DO K1=1,K4
!!$     WRITE(25,*)YM1(K1),YM2(K1)
!!$  END DO
!!$  CLOSE(25)
!!$ OPEN(UNIT=35, FILE='xymin_ma.dat.txt')
!!$  DO K1=1,K4
!!$     WRITE(35,*)ZM1(K1),ZM2(K1)
!!$  END DO
!!$  CLOSE(35)
  open(unit=12, file='rand.dat.txt')
  write(12,*)nint(-ran3(iw)*1.d+6)
  close(12)
END PROGRAM CHAOS
!==========================================================================
!  ran3 from Numerical Recipes, 2nd edition
!--------------------------------------------------------------------------
      function ran3(idum)
!==========================================================================

!         implicit real*4(m)

      integer idum
      integer mbig,seed,mz
!      real mbig,seed,mz
      double precision ran3,fac

!         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
!      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
 
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
!      real  mj,mk,ma(55)

      save iff,inext,inextp,ma

      data iff /0/

1      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if(ran3.le.0.or.ran3.ge.1) goto 1
!!      if(ran3.le.0.or.ran3.ge.1) Then
      !! write(6,*)'RAN3 failed: ', ran3
      ! stop
!!      endif
      return
      end function ran3



       SUBROUTINE derivs(y,dy,n,a,b)
! Y=Y_i, DY= Y'_i, n-dim
      INTEGER n,k1,k2,k3
      DOUBLE PRECISION y(n),dy(n),a(n,n,n),b(n,n)
      DOUBLE PRECISION S,R,BE,AA
      AA=1.D0
      dy=0.d0
      do k1=1,n
         dy(k1)=-y(k1)**3
!         dy(k1)=-y(k1)*dabs(y(k1))*n*n
          do k2=1,n
            dy(k1)=dy(k1) + b(k1,k2)*y(k2) 
            do k3=1,n
            dy(k1)=dy(k1) + a(k1,k2,k3)*y(k2)*y(k3)  
            end do   
            end do
            end do


      RETURN
       END SUBROUTINE DERIVS
         
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




       
