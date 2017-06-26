        FUNCTION ran1(IDUM)

       implicit none

!  RAN1 returns a unifom random deviate on the interval [0,1]
! __________________________________________________________________
!
      INTEGER :: IDUM
      REAL*8 :: RAN2,ran1
      integer,parameter :: IM1=2147483563,IM2=2147483399
      integer,parameter :: IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, NTAB=32
        integer,parameter :: NDIV=1+IMM1/NTAB
         real*8,parameter :: EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1
        INTEGER :: IDUM2,J,K,IV(NTAB),IY
        DATA IDUM2/123456789/, iv/NTAB*0/, iy/0/
        IF (IDUM.LE.0) THEN
          IDUM=MAX(-IDUM,1)
          IDUM2=IDUM
          DO  J=NTAB+8,1,-1
             K=IDUM/IQ1
             IDUM=IA1*(IDUM-K*IQ1)-K*IR1
             IF (IDUM.LT.0) IDUM=IDUM+IM1
             IF (J.LE.NTAB) IV(J)=IDUM
          end do
          IY=IV(1)

        ENDIF
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF (IDUM.LT.0) IDUM=IDUM+IM1
        K=IDUM2/IQ2
        IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
        IF (IDUM2.LT.0) IDUM2=IDUM2+IM2

         J=1+IY/NDIV
        IY=IV(J)-IDUM2
        IV(J)=IDUM
        IF(IY.LT.1)IY=IY+IMM1
        RAN1=MIN(AM*IY,RNMX)
        END function ran1

