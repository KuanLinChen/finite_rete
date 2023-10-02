      SUBROUTINE HEAT2A(DENS,TEM,PRE,RMXBAR,GAMBAR,CPBAR,HHBAR,VOFP)
C*START*CHEMISTRY PACKAGE*********************************************
C---------------------------------------------------------------------
C     FINITE-RATE CHEMISTRY DRIVER
C---------------------------------------------------------------------
cnvx*CALL fdns01
      include 'fdns01'
cnvx*CALL fdns02
      include 'fdns02'
cnvx*CALL fdns03
      include 'fdns03'
cnvx*CALL fdns04
      include 'fdns04'
      include 'fdns10'
      include 'fdns11'
      COMMON /THRMPT/CP1(NSPM),HP1(NSPM),ST1(NSPM),GF1(NSPM)
C-GC  DIMENSION XN1(NSPM),XN2(NSPM),XN3(NSPM)
C-GC &         ,XN0(NSPM)
      DIMENSION XN1(50),XN2(50),XN3(50),XN0(50)
      DIMENSION WDT(NSPM),F(NSPM),FM3(50),XKN(50),XKF(50)
C-GC
      DIMENSION A(50,50),B(50),BB(50),XLAM(50,50),GMOLCC(50),A2(50,50)
      DIMENSION CUBA(4),CUBR(3)
C-GC
      ISPROD=1+INT(INSO(12)/10)
      IF(INSO(12) .GE. 100) ISPROD=ISPROD-10
C
C     ISPROD = 1: EXPLICIT PRODUCTION SOURCE TERM
C              2: IMPLICIT PRODUCTION SOURCE TERM WITH PSEUDO-TIME STEP
C              3: IMPLICIT PRODUCTION SOURCE TERM WITH REAL TIME STEP
C              4: IMPLICIT PRODUCTION SOURCE TERM WITH REAL TIME STEP,
C                 AND SUB-ITERATION, CONSTANT TEMPERATURE
C              5: IMPLICIT PRODUCTION SOURCE TERM WITH REAL TIME STEP,
C                 AND SUB-ITERATION, CONSTANT ENTHALPY
C-GC
C
C-----INITIALIZE PARAMETERS
      RGG   = 8314.6
      PQ1   = PRE
      TQ1   = AMAX1(1.000,TEM)
C     RTQ1  = 1.0/TQ1
C     XTQ1  = ALOG(TQ1)
C-----START CALCULATIONS
      DO 100 J=1,NGAS
        IF(ALPHA(J).LE.1.0E-29) ALPHA(J) = 0.0
        WDT(J) = 0.0
        XN2(J) = 0.0
 100  CONTINUE
C-----BASIC THERMODYNAMICS PROPERTIES
      CALL CPHG(TQ1)
C-GC
      IF(INSO(7) .EQ. 0) THEN
        DUM   = 0.0
        RMXBAR= 0.0
        HHBAR = 0.0
        CPBAR = 0.0
        DO J=1,NGAS
          XN3(J)= ALPHA(J)*WTMOLE(J)
          DUM   = DUM+XN3(J)
        ENDDO
        DO J=1,NGAS
          XN3(J)= XN3(J)/DUM
          RMXBAR= RMXBAR+ALPHA(J)*RG(J)
          HHBAR = HHBAR+XN3(J)*HP1(J)
          CPBAR = CPBAR+XN3(J)*CP1(J)
        ENDDO
        GAMBAR= CPBAR/(CPBAR-1.)
        CPBAR = CPBAR*RMXBAR
        HHBAR = HHBAR*RMXBAR*TQ1
      ENDIF
C
      DUM   = 0.0
      DO J=1,NGAS
        XN3(J) = ALPHA(J)*WTMOLE(J)
        DUM    = DUM+XN3(J)
        XN1(J) = ALPHA(J)
C       IF(INSO(7) .EQ. J) ALPHA(J) = ALPHA(J)*VOFP
      END DO
      DO J=1,NGAS
        XN3(J)= XN3(J)/DUM
      END DO
C-GC
C-----FOR FINITE-RATE CHEMICAL REACTION WHEN NREACT > 0
C     IF(NREACT.GT.0.AND.TQ1.GT.300.0.AND.DENS.GT.0.0) THEN
      IF(NREACT.GT.0.0.AND.DENS.GT.0.0) THEN
        IF(ISPROD .EQ. 1) THEN
C-GC
C
C-----GET CHEMICAL POTENTIALS (GIBBS FREE ENERGY)
          FACT  = 1.0E-03
          FACTR = 1.0E+03
          SUM3  = 0.0
          DO 1100 J=1,NGAS
            GMOLCC(J)=DENS*ALPHA(J)*WTMOLE(J)*FACT
            SUM3   = SUM3+GMOLCC(J)
            F(J)   = HP1(J)-ST1(J)
 1100     CONTINUE
C-------EQUILIBRIUM AND FORWARD-RATE CONSTANTS
          RGGCGS= 82.047163
c         TQZ   = AMAX1(1200.,TQ1)
          TQZ   = TQ1
          TQX   = 1.0/(RGGCGS*TQZ)
          TQY   = 1.0/TQZ
          DO 1225 I=1,NREACT
            RNT = 0.0
            RRF = 0.0
            FM3(I) = 1.0
            IF(ITHIRD(I).GE.1.AND.ITHIRD(I).LE.998) THEN
              FM3(I) = ALPHA(ITHIRD(I))*WTMOLE(ITHIRD(I))*DENS*FACT
            ELSE IF(ITHIRD(I).GE.999) THEN
              FM3(I) = SUM3
            ENDIF
            DO 1226 J=1,NGAS
              RNT = RNT+STCOEF(I,J)
              RRF = RRF+STCOEF(I,J)*F(J)
 1226       CONTINUE
            XKN(I) = TQX**RNT*EXP(AMAX1(-65.0,AMIN1(65.0,-RRF)))
C-GC        XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
C-GC        ******  SOOT OXIDATION & FORMATION ******
            IF(ARRHA(I) .LE. 1.E-30) THEN
              WMSOOT=1./WTMOLE(NGAS)
              FNO2  =XN3(2)
              DNSOOT=1.86
              DMSOOT=ARRHB(I)
              CALL SOOTOX(PRE,TQ1,DNSOOT,DMSOOT,WMSOOT,FNO2,XKFS)
              XKF(I) = XKFS
            ELSE IF(ARRHA(I) .LE. 1.0 .AND. ARRHA(I) .GT. 1.E-30) THEN
              XKF(I) = ARRHA(I)*EXP(-ARRHB(I)*(1.-TQ1/1700.)**2)
            ELSE
              XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
            ENDIF
C-GC        ******************************
 1225     CONTINUE
C-------RATE CALCULATION
          DO 1300 I=1,NREACT
            RATEF = 1.0
            RATEB = 1.0
C--REGULAR
            IF(IGLOB(I).EQ.0) THEN
              DO 1310 J=1,NGAS
C-GC            DENAW = DENS*ALPHA(J)*WTMOLE(J)*FACT
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(ALPHA(J).GT.1.0E-30) THEN
C-GC                RATEF = RATEF*DENAW**(-STCOEF(I,J))
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(ALPHA(J).GT.1.0E-30) THEN
C-GC                RATEB = RATEB*DENAW**( STCOEF(I,J))
                    RATEB = RATEB*GMOLCC(J)**( STCOEF(I,J))
                  ELSE
                    RATEB = 0.0
                  ENDIF
                ENDIF
 1310         CONTINUE
C--GLOBAL-1
            ELSE IF(IGLOB(I).EQ.1) THEN
              RATEB = 0.0
              DO 1311 J=1,NGAS
C-GC            DENAW = DENS*ALPHA(J)*WTMOLE(J)*FACT
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(ALPHA(J).GT.1.0E-30) THEN
C-GC                RATEF = RATEF*DENAW**(-STCOEF(I,J))
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ENDIF
 1311         CONTINUE
C--GLOBAL-2
            ELSE IF(IGLOB(I).EQ.2) THEN
              RATEB = 0.0
              DO 1312 J=1,NGAS
C-GC            DENAW = DENS*ALPHA(J)*WTMOLE(J)*FACT
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(ALPHA(J).GT.1.0E-30) THEN
C-GC                RATEF = RATEF*DENAW**(-STCOEG(I,J))
                    RATEF = RATEF*GMOLCC(J)**(-STCOEG(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  RATEF = RATEF*GMOLXX**(-STCOEG(I,J))
C-GC     *************************************
                ENDIF
 1312         CONTINUE
            ENDIF
C--FINAL ASSEMBLY
            DO 1320 J=1,NGAS
              WDT(J) = WDT(J)+FM3(I)*FACTR*(STCOEF(I,J)/WTMOLE(J))
     &                *(RATEF-RATEB/AMAX1(1.E-30,XKN(I)))*XKF(I)
              ALPHA(J)=XN1(J)
 1320       CONTINUE
 1300     CONTINUE
C-GC
        ELSE IF(ISPROD .EQ. 2) THEN
          NORDER=MIN(2,MAX(1,MOD(INSO(12),10)))
          PORDER=1./FLOAT(NORDER)
          DTQ=1./AMIN1(1.00,ABS(DTT))
          DTDIM=DTQ
          FACT=1.0E-03
          FACTR=1.0E+03
C-GC      DO 2100 I=1,NREACT
          DO I=1,NREACT
          DO 2100 J=1,NGAS  
            XLAM(I,J)=0.0
 2100     CONTINUE
          ENDDO
          SUM3=0.0
          DO 2150 J=1,NGAS
            GMOLCC(J)=DENS*ALPHA(J)*WTMOLE(J)*FACT
            SUM3=SUM3+GMOLCC(J)
            F(J)=HP1(J)-ST1(J)
 2150     CONTINUE
C
C-------EQUILIBRIUM AND FORWARD-RATE CONSTANTS
C
          RGGCGS= 82.047163
          TQZ   = TQ1
C         TQZ   = AMAX1(1500.,TQ1)
          TQX   = 1.0/(RGGCGS*TQZ)
          TQY   = 1.0/TQZ
          DO 2200 I=1,NREACT
            RNT = 0.0
            RRF = 0.0
            FM3(I) = 1.0
            IF(ITHIRD(I).GE.1.AND.ITHIRD(I).LE.998) THEN
              FM3(I) = ALPHA(ITHIRD(I))*WTMOLE(ITHIRD(I))*DENS*FACT
            ELSE IF(ITHIRD(I).GE.999) THEN
              FM3(I) = SUM3
            ENDIF
            DO 2250 J=1,NGAS
              RNT = RNT+STCOEF(I,J)
              RRF = RRF+STCOEF(I,J)*F(J)
 2250       CONTINUE
            XKN(I) = TQX**RNT*EXP(AMAX1(-65.0,AMIN1(65.0,-RRF)))
            XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
 2200     CONTINUE
C
C     CALCULATE RATES FOR EACH REACTION
C
          DO 2300 I=1,NREACT
            RATEF = FM3(I)
            RATEB = FM3(I)
C--REGULAR
            IF(IGLOB(I).EQ.0) THEN
              DO 2400 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEB = RATEB*GMOLCC(J)**( STCOEF(I,J))
                  ELSE
                    RATEB = 0.0
                  ENDIF
                ENDIF
 2400         CONTINUE
              RATEF=XKF(I)*RATEF
              RATEB=XKF(I)*RATEB/XKN(I)
              DO 2450 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEB/GMOLCC(J)
                  ENDIF
                ENDIF
 2450         CONTINUE
              BB(I)=RATEF-RATEB
C--GLOBAL-1
            ELSE IF(IGLOB(I).EQ.1) THEN
              RATEB = 0.0
              DO 2500 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ENDIF
 2500         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 2550 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ENDIF
 2550         CONTINUE
              BB(I)=RATEF
C--GLOBAL-2
            ELSE IF(IGLOB(I).EQ.2) THEN
              RATEB = 0.0
              DO 2600 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEG(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  RATEF = RATEF*GMOLXX**(-STCOEG(I,J))
C-GC
                ENDIF
 2600         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 2650 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLCC(J)
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLXX
C-GC
                ENDIF
 2650         CONTINUE
              BB(I)=RATEF
            ENDIF
 2300     CONTINUE 
c
c     Fill arrays for Gauss elimination
c
          DO 2700 K=1,NGAS
            B(K)=0.0
            DO 2750 I=1,NREACT
              B(K)=B(K)+STCOEF(I,K)*BB(I)
 2750       CONTINUE
            DO 2800 J=1,NGAS
              A(K,J)=0.0
              IF(J .EQ. K) A(K,J)=DTDIM
              DO 2850 I=1,NREACT
                A(K,J)=A(K,J)-PORDER*STCOEF(I,K)*XLAM(I,J)
 2850         CONTINUE
 2800       CONTINUE
 2700     CONTINUE
C
C
          CALL GAUSS(A,B,GMOLCC,NGAS)
c
c     chemistry rate production terms
c
          DO 2900 J=1,NGAS
            WDT(J)=GMOLCC(J)*DTDIM*FACTR/WTMOLE(J)
            ALPHA(J)=XN1(J)
 2900     CONTINUE
        ELSE IF(ISPROD .EQ. 3) THEN
          NORDER=MIN(2,MAX(1,MOD(INSO(12),10)))
          PORDER=1./FLOAT(NORDER)
          DTQ=ABS(DTT)
          DTDIM=DTQ*XREF/UREF
          FACT=1.0E-03
          FACTR=1.0E+03
C-GC      DO 3100 I=1,NREACT
          DO I=1,NREACT
          DO 3100 J=1,NGAS  
            XLAM(I,J)=0.0
 3100     CONTINUE
          ENDDO
          SUM3=0.0
          DO 3150 J=1,NGAS
            GMOLCC(J)=DENS*ALPHA(J)*WTMOLE(J)*FACT
            SUM3=SUM3+GMOLCC(J)
            F(J)=HP1(J)-ST1(J)
 3150     CONTINUE
C
C-------EQUILIBRIUM AND FORWARD-RATE CONSTANTS
C
          RGGCGS= 82.047163
          TQZ   = TQ1
C         TQZ   = AMAX1(1500.,TQ1)
          TQX   = 1.0/(RGGCGS*TQZ)
          TQY   = 1.0/TQZ
          DO 3200 I=1,NREACT
            RNT = 0.0
            RRF = 0.0
            FM3(I) = 1.0
            IF(ITHIRD(I).GE.1.AND.ITHIRD(I).LE.998) THEN
              FM3(I) = ALPHA(ITHIRD(I))*WTMOLE(ITHIRD(I))*DENS*FACT
            ELSE IF(ITHIRD(I).GE.999) THEN
              FM3(I) = SUM3
            ENDIF
            DO 3250 J=1,NGAS
              RNT = RNT+STCOEF(I,J)
              RRF = RRF+STCOEF(I,J)*F(J)
 3250       CONTINUE
            XKN(I) = TQX**RNT*EXP(AMAX1(-65.0,AMIN1(65.0,-RRF)))
            XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
 3200     CONTINUE
C
C     CALCULATE RATES FOR EACH REACTION
C
          DO 3300 I=1,NREACT
            RATEF = FM3(I)
            RATEB = FM3(I)
C--REGULAR
            IF(IGLOB(I).EQ.0) THEN
              DO 3400 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEB = RATEB*GMOLCC(J)**( STCOEF(I,J))
                  ELSE
                    RATEB = 0.0
                  ENDIF
                ENDIF
 3400         CONTINUE
              RATEF=XKF(I)*RATEF
              RATEB=XKF(I)*RATEB/XKN(I)
              DO 3450 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEB/GMOLCC(J)
                  ENDIF
                ENDIF
 3450         CONTINUE
              BB(I)=RATEF-RATEB
C--GLOBAL-1
            ELSE IF(IGLOB(I).EQ.1) THEN
              RATEB = 0.0
              DO 3500 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ENDIF
 3500         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 3550 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ENDIF
 3550         CONTINUE
              BB(I)=RATEF
C--GLOBAL-2
            ELSE IF(IGLOB(I).EQ.2) THEN
              RATEB = 0.0
              DO 3600 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEG(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  RATEF = RATEF*GMOLXX**(-STCOEG(I,J))
C-GC
                ENDIF
 3600         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 3650 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLCC(J)
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLXX
C-GC
                ENDIF
 3650         CONTINUE
              BB(I)=RATEF
            ENDIF
 3300     CONTINUE 
c
c     Fill arrays for Gauss elimination
c
          DO 3700 K=1,NGAS
            B(K)=0.0
            DO 3750 I=1,NREACT
              B(K)=B(K)+STCOEF(I,K)*BB(I)
 3750       CONTINUE
            DO 3800 J=1,NGAS
              A(K,J)=0.0
              IF(J .EQ. K) A(K,J)=1.
              DO 3850 I=1,NREACT
                A(K,J)=A(K,J)-PORDER*STCOEF(I,K)*XLAM(I,J)*DTDIM
 3850         CONTINUE
 3800       CONTINUE
 3700     CONTINUE
C
          CALL GAUSS(A,B,GMOLCC,NGAS)
C
          DO 3900 J=1,NGAS
            WDT(J)=GMOLCC(J)*FACTR/WTMOLE(J)
            ALPHA(J)=XN1(J)
 3900     CONTINUE
        ELSE IF(ISPROD .GE. 4 .AND. ISPROD .LE. 6) THEN
          NORDER=MIN(3,MOD(INSO(12),10))
          DTQ=ABS(DTT)*XREF/UREF
          DTDIM=1./DTQ
C-----  Time step control
c         XDT=1.0E-8
c         DTMIN=1.E-15
c         ATOL=1.E-8
c         ESP=0.05
c         NSTEP=2000
C-----GROUP 1 ACCURACY CONTROL 
          ESP   = 0.20     ! time step increment control 
          DTMIN = 1.0E-10  ! minimum time step size
          ATOL  = 1.0E-8   ! species sensitivity
          NSTEP = 500      ! maximum integration steps
C
C-----GROUP 2 ACCURACY CONTROL 
c         ESP   = 0.02     ! time step increment control 
c         DTMIN = 1.0E-12  ! minimum time step size
c         ATOL  = 1.0E-10  ! species sensitivity
c         NSTEP = 5000     ! maximum integration steps
C
C-----GROUP 3 ACCURACY CONTROL 
c         ESP   = 0.01     ! time step increment control 
c         DTMIN = 1.0E-15  ! minimum time step size
c         ATOL  = 1.0E-12  ! species sensitivity
c         NSTEP = 50000    ! maximum integration steps
C
          XDT = DTMIN 
          DTMAX=0.5*DTQ
          DTINC=2.
          TOTIME=0.
          ISTEP=0
          FACT=1.0E-03
          FACTR=1.0E+03
          TMOLN=0.
C------ 1st order implicit
          C1JAC=-1.
          C2JAC=0.
C------ 2nd order parasol 
          IF(NORDER .EQ. 2) THEN
            C1JAC=-0.5
C------ 4th order parasol 
          ELSE IF(NORDER .EQ. 3) THEN
            C1JAC=-0.5
            C2JAC=1./12
          END IF
          DO 4010 J=1,NGAS
            XN2(J)=ALPHA(J)
            XN3(J)=XN1(J)*WTMOLE(J)
C           XN3(J)=ALPHA(J)*WTMOLE(J)
            TMOLN=TMOLN+XN3(J)
 4010     CONTINUE
          DO 4020 J=1,NGAS
            XN3(J)=XN3(J)/TMOLN
 4020     CONTINUE
C------- sub-iteration loop for implicit chemistry
C
 4001     CONTINUE
          IF(TOTIME .GE. 0.999*DTQ .OR. ISTEP .GE. NSTEP) GO TO 4999
          ISTEP = ISTEP + 1
          IF(ISTEP .GT. 1 .AND. ISPROD .EQ. 5) THEN
            HHTMP=HHBAR
            CALL HEAT2B(DENS,TQ1,PRE,HHTMP,TMTMP)
            IF(ABS(TQ1-TMTMP) .GE. 1.) CALL CPHG(TMTMP)
C           IF(ABS(TQ1-TMTMP) .GE. 1. .AND. TMTP .GT. 1000.) THEN
C             CALL CPHG(TMTMP)
              TQ1=TMTMP
C           ENDIF
          END IF
C-GC      DO 4100 I=1,NREACT
          DO I=1,NREACT
          DO 4100 J=1,NGAS  
            XLAM(I,J)=0.0
 4100     CONTINUE
          ENDDO
          SUM3=0.0
          DO 4110 J=1,NGAS
            GMOLCC(J)=DENS*ALPHA(J)*WTMOLE(J)*FACT
            SUM3=SUM3+GMOLCC(J)
            XN0(J)=0.
            F(J)=HP1(J)-ST1(J)
 4110     CONTINUE
C
C-------EQUILIBRIUM AND FORWARD-RATE CONSTANTS
C
          RGGCGS= 82.047163
          TQZ   = TQ1
c         TQZ   = AMAX1(1200.,TQ1)
          TQX   = 1.0/(RGGCGS*TQZ)
          TQY   = 1.0/TQZ
          DO 4200 I=1,NREACT
            RNT = 0.0
            RRF = 0.0
            FM3(I) = 1.0
            IF(ITHIRD(I).GE.1.AND.ITHIRD(I).LE.998) THEN
              FM3(I) = ALPHA(ITHIRD(I))*WTMOLE(ITHIRD(I))*DENS*FACT
            ELSE IF(ITHIRD(I).GE.999) THEN
              FM3(I) = SUM3
            ENDIF
            DO 4250 J=1,NGAS
              RNT = RNT+STCOEF(I,J)
              RRF = RRF+STCOEF(I,J)*F(J)
 4250       CONTINUE
            XKN(I) = TQX**RNT*EXP(AMAX1(-65.0,AMIN1(65.0,-RRF)))
C-GC        XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
C-GC        ******  SOOT OXIDATION & FORMATION ******
            IF(ARRHA(I) .LE. 1.E-30) THEN
              WMSOOT=1./WTMOLE(NGAS)
              FNO2  =XN3(2)
              DNSOOT=1.86
              DMSOOT=ARRHB(I)
              CALL SOOTOX(PRE,TQ1,DNSOOT,DMSOOT,WMSOOT,FNO2,XKFS)
              XKF(I) = XKFS
            ELSE IF(ARRHA(I) .LE. 1.0 .AND. ARRHA(I) .GT. 1.E-30) THEN
              XKF(I) = ARRHA(I)*EXP(-ARRHB(I)*(1.-TQ1/1700.)**2)
            ELSE
              XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
            ENDIF
C-GC        ******************************
 4200     CONTINUE
C
C     CALCULATE RATES FOR EACH REACTION
C
          DO 4300 I=1,NREACT
            RATEF = FM3(I)
            RATEB = FM3(I)
C--REGULAR
            IF(IGLOB(I).EQ.0) THEN
              DO 4400 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEB = RATEB*GMOLCC(J)**( STCOEF(I,J))
                  ELSE
                    RATEB = 0.0
                  ENDIF
                ENDIF
 4400         CONTINUE
              RATEF=XKF(I)*RATEF
              RATEB=XKF(I)*RATEB/XKN(I)
              DO 4450 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ELSE IF(STCOEF(I,J).GT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEB/GMOLCC(J)
                  ENDIF
                ENDIF
 4450         CONTINUE
              BB(I)=RATEF-RATEB
C--GLOBAL-1
            ELSE IF(IGLOB(I).EQ.1) THEN
              RATEB = 0.0
              DO 4500 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEF(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
                ENDIF
 4500         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 4550 J=1,NGAS
                IF(STCOEF(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEF(I,J)*RATEF/GMOLCC(J)
                  ENDIF
                ENDIF
 4550         CONTINUE
              BB(I)=RATEF
C--GLOBAL-2
            ELSE IF(IGLOB(I).EQ.2) THEN
              RATEB = 0.0
              DO 4600 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    RATEF = RATEF*GMOLCC(J)**(-STCOEG(I,J))
                  ELSE
                    RATEF = 0.0
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  RATEF = RATEF*GMOLXX**(-STCOEG(I,J))
C-GC     *************************************
                ENDIF
 4600         CONTINUE
              RATEF=XKF(I)*RATEF
              DO 4650 J=1,NGAS
                IF(STCOEG(I,J).LT.0.0) THEN
                  IF(GMOLCC(J).GT.1.0E-30) THEN
                    XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLCC(J)
                  ENDIF
C-GC     *****  for Soot formation  **********
                ELSE IF(STCOEG(I,J) .GT. 0.0) THEN
                  GMOLXX= DENS*(ALPHA(J)+0.001)*WTMOLE(J)*FACT 
                  XLAM(I,J)=-STCOEG(I,J)*RATEF/GMOLXX 
C-GC     *************************************
                ENDIF
 4650         CONTINUE
              BB(I)=RATEF
            ENDIF
 4300     CONTINUE 
c
c     Fill arrays for Gauss elimination
c
          DO 4700 K=1,NGAS
            B(K)=0.0
            DO 4705 I=1,NREACT
C-GC          B(K)=B(K)+STCOEF(I,K)*BB(I)
C-GC        ***********  SOOT OXIDATION  ***********
              IF(ARRHA(I) .LE. 1.E-30) THEN
                B(K)=B(K)+STCOEF(I,K)*BB(I)*WTMOLE(K) 
              ELSE
                B(K)=B(K)+STCOEF(I,K)*BB(I)
              ENDIF
C-GC        ****************************************
 4705       CONTINUE
            DO 4710 J=1,NGAS
              A(K,J)=0.0
              A2(K,J)=0.0
              DO 4720 I=1,NREACT
C-GC            A(K,J)=A(K,J)+STCOEF(I,K)*XLAM(I,J)
C-GC            ***********  SOOT OXIDATION  ***********
                IF(ARRHA(I) .LE. 1.E-30) THEN
                  A(K,J)=A(K,J)+STCOEF(I,K)*XLAM(I,J)*WTMOLE(K)
                ELSE
                  A(K,J)=A(K,J)+STCOEF(I,K)*XLAM(I,J)
                ENDIF
C-GC            ****************************************
 4720         CONTINUE
 4710       CONTINUE
 4700     CONTINUE
          IF(ISTEP .EQ. 1) THEN
            XDT=DTQ
            DO 4730 J=1,NGAS
              B1=B(J)*FACTR/DENS/TMOLN
              XDT=AMIN1(XDT,0.1*(XN3(J)+ATOL)/ABS(B1+ATOL))
 4730       CONTINUE
            XDT=AMIN1(DTMAX,AMAX1(XDT,DTMIN))
            XDT1=XDT
          END IF
          IF(NORDER .EQ. 3) THEN 
C-GC        DO 4740 I=1,NGAS
            DO I=1,NGAS
            DO 4740 J=1,NGAS
              SQJAC=0.0
              DO 4750 K=1,NGAS
                SQJAC=SQJAC+A(I,K)*A(K,J)
 4750         CONTINUE
              A2(I,J)=SQJAC
 4740       CONTINUE
            ENDDO
          END IF
          DO 4760 I=1,NGAS
          DO 4770 J=1,NGAS
C-GC        A(I,J)=C1JAC*A(I,J)+C2JAC*A2(I,J)*XDT
            A(I,J)=C1JAC*A(I,J)*XDT+C2JAC*A2(I,J)*XDT*XDT
 4770       CONTINUE
C-GC        A(I,I)=A(I,I)+1./XDT
            A(I,I)=A(I,I)+1.
            B(I)=B(I)*XDT
 4760     CONTINUE
C
          CALL GAUSS(A,B,GMOLCC,NGAS)
C------  Time Step Control
          DALM=0.
          DO 4800 J=1,NGAS
            XK1=GMOLCC(J)*FACTR/DENS/TMOLN
            IF(XN3(J)+XK1 .LT. 0.) XK1=0.
            ALFA1=XN3(J)+XK1
            AUX1=XK1/AMIN1(XN3(J)+ATOL,ALFA1+ATOL)
            AUX2=XN0(J)/(XN3(J)+ATOL)
            AUX3=AMIN1(ABS(AUX1-AUX2),2.*ABS(AUX1))
            DALM=AMAX1(DALM,AUX3)
 4800     CONTINUE
          IF(ABS(DALM) .LT. 1.E-15) GO TO 4999
          IF(XDT .LE. 1.1*DTMIN) DTINC=1.
          DTINC1=DTINC
          DTINC=SQRT(ESP/DALM)
          DTINC=AMIN1(10.,AMAX1(0.3,DTINC,0.2*DTINC1),2.+DTINC1)
          XDT=DTINC*XDT
          XDTEND=DTQ-TOTIME-XDT1
          IF(XDT .GT. 0.8*XDTEND) XDT=XDTEND
          XDT=AMAX1(XDT,DTMIN)
          TOTIME=TOTIME+XDT1
          XDT1=XDT
          DO 4810 J=1,NGAS
            XK1=GMOLCC(J)*FACTR/DENS/TMOLN
            IF(XN3(J)+XK1 .LT. 0.) XK1=0.
            XN0(J)=XK1
            XN3(J)=XN3(J)+XK1
            ALPHA(J)=ALPHA(J)+XK1/WTMOLE(J)*TMOLN
 4810     CONTINUE
          GO TO 4001
 4999     CONTINUE
          IF(ISTEP .GE. NSTEP) THEN
            WRITE(36,*) '@@@ SPECIES RATE EXCEED NSTEP AT IJK = ',
     &        IJKQC,', DT =',TOTIME
            DTDIM=1./TOTIME
          ENDIF
c
c     chemistry rate production terms
c
          DO 4950 J=1,NGAS
            WDT(J)=(ALPHA(J)-XN2(J))*DENS*DTDIM
            ALPHA(J)=XN1(J)
            XN2(J)=0.
 4950     CONTINUE
          TFMIN(IJKQC) = DTDIM*XREF/UREF
        END IF
C-GC
      END IF
C-END-FINITE RATE CHEMISTRY-------------------------------------------
C-GC
      IF(DENS.GT.0.0) THEN
        DO 305 J=1,NGAS
          FMDOT(IJKQC,J) = WDT(J)
 305    CONTINUE
      ENDIF
C-END-HEAT2A----------------------------------------------------------
      RETURN
      END
C-GC
      SUBROUTINE SOOTOX(PRE,TEM,DNSOOT,DMSOOT,WMSOOT,FNO2,XKF)
C
      DIMENSION XRA(4),XRER(4),XRK(4)
      DATA XRA/2.E+1,4.46E-3,1.51E+5,2.13E+1/
      DATA XRER/1.509E+4,7.649E+3,4.882E+4,-2.063E+3/
C
      RATM=82.06
      DO 100 I=1,4
        XRK(I)=XRA(I)*EXP(-XRER(I)/TEM)
  100 CONTINUE
      PREO2=FNO2*PRE*1.0E-5
      DUM=PREO2*XRK(3)/XRK(2)
      XXX=1./(1.+DUM)
      P1=XRK(1)*XXX/(1.+XRK(4)*PREO2)
      P2=(1.-XXX)*XRK(2)
      XKF=6.*WMSOOT*RATM*TEM*(P1+P2)/(DNSOOT*DMSOOT)
      RETURN
      END
C
      Subroutine GAUSS(A,B,X,N)
c  
c     Gaussian Elimination Using Pivoting for Ax=b
c  
      dimension A(50,50),B(50),X(50),S(50)
      dimension L(50)
c  
c     Scale vector S equals reciprocal of largest element in row i
c  
      do 200 i=1,n
      L(i)=i
      smax=0.0
      do 100 j=1,n
      smax=amax1(smax,abs(A(i,j)))
  100 continue
      S(i)=1.0/smax
  200 continue
      do 600 j=1,n-1 
c  
c     Index vector L orders rows by largest |A(i,j)|/S(i)
c  
      rmax=0.0
      do 300 i=j,n
      irow=L(i)
      ratio=abs(A(irow,j))*S(irow)
      if(ratio.gt.rmax) then
        inew=i
        rmax=ratio
      endif
  300 continue
      iold=L(inew)
      L(inew)=L(j)
      L(j)=iold
c
c     
      do 500 i=j+1,n 
      irow=L(i)
      Arat=A(irow,j)/A(iold,j)
      do 400 jj=j+1,n
      A(irow,jj)=A(irow,jj)-Arat*a(iold,jj)
  400 continue
      B(irow)=B(irow)-Arat*B(iold)
  500 continue
  600 continue
c
c     Back Substitution
c
      irow=L(n)
      X(n)=B(irow)/A(irow,n)
      do 1000 i=n-1,1,-1
      irow=L(i)
      sum=B(irow)
      do 900 j=i+1,n
      sum=sum-A(irow,j)*X(j)
  900 continue
      X(i)=sum/A(irow,i)
 1000 continue
      return
      end
C
      SUBROUTINE CPHG(TEM)
      include 'fdns01'
      include 'fdns02'
      COMMON /THRMPT/CP1(NSPM),HP1(NSPM),ST1(NSPM),GF1(NSPM)
C-----INITIALIZE PARAMETERS
      TQ1   = AMAX1(1.000,TEM)
      RTQ1  = 1.0/TQ1
      XTQ1  = ALOG(TQ1)
C-----BASIC THERMODYNAMICS PROPERTIES
      K      = 2
      IF(TQ1.GT.1000.0) K = 1
      IF(TQ1.GE.300.0.AND.TQ1.LE.5000.0) THEN
C-------FOR 300 <= TQ1 <= 5000
        DO 150 J=1,NGAS
          CP1(J) = HF(1,K,J)+(HF(2,K,J)+(HF(3,K,J)+(HF(4,K,J)
     &    +HF(5,K,J)*TQ1)*TQ1)*TQ1)*TQ1
          HP1(J) = HF(1,K,J)+(.5*HF(2,K,J)+(HF(3,K,J)/3.+(.25
     &    *HF(4,K,J)+.20*HF(5,K,J)*TQ1)*TQ1)*TQ1)*TQ1+HF(6,K,J)*RTQ1
          ST1(J) = HF(1,K,J)*XTQ1+(HF(2,K,J)+(.5*HF(3,K,J)+(.3333333
     &    *HF(4,K,J)+.25*HF(5,K,J)*TQ1)*TQ1)*TQ1)*TQ1+HF(7,K,J)
          GF1(J) = HP1(J)-ST1(J)
 150    CONTINUE
      ELSE
C-------FOR 5000 < TQ1 OR 300 > TQ1
        TQQ = 5000.0
        IF(TQ1.LT.300.0) TQQ = 300.0
        RTQQ = 1.0/TQQ
        XTQQ = ALOG(TQQ)
        DO 151 J=1,NGAS
          CP1(J) = HF(1,K,J)+(HF(2,K,J)+(HF(3,K,J)+(HF(4,K,J)
     &    +HF(5,K,J)*TQQ)*TQQ)*TQQ)*TQQ
          HP1(J) = ((HF(1,K,J)+(.5*HF(2,K,J)+(HF(3,K,J)/3.+(.25
     &    *HF(4,K,J)+.20*HF(5,K,J)*TQQ)*TQQ)*TQQ)*TQQ)*TQQ+HF(6,K,J)
     &    +CP1(J)*(TQ1-TQQ))*RTQ1
          ST1(J) = HF(1,K,J)*XTQQ+(HF(2,K,J)+(.5*HF(3,K,J)+(.3333333
     &    *HF(4,K,J)+.25*HF(5,K,J)*TQQ)*TQQ)*TQQ)*TQQ+HF(7,K,J)
     &    +CP1(J)*(XTQ1-XTQQ)
          GF1(J) = HP1(J)-ST1(J)
 151    CONTINUE
      ENDIF
      RETURN
      END
C-GC
      SUBROUTINE HEAT2B(DENS,TEM,PRE,HHBAR,TMTMP)
C---------------------------------------------------------------------
C     GEN TEMPERATURE FROM KNOWN ENTHALPY
C---------------------------------------------------------------------
cnvx*CALL fdns01
      include 'fdns01'
cnvx*CALL fdns02
      include 'fdns02'
cnvx*CALL fdns11
      include 'fdns11'
C
C-----INITIALIZE PARAMETERS
        RGG   = 8314.6
        RMXBAR= 0.0
        DUM   = 0.0
        TQ1   = AMAX1(1.000,TEM)
        TMTMP = TQ1
C-----START CALCULATIONS
        DO 100 J=1,NGAS
         IF(ALPHA(J).LE.1.0E-29) ALPHA(J) = 0.0
         DUM    = DUM+ALPHA(J)*WTMOLE(J)
         RMXBAR = RMXBAR+ALPHA(J)*RG(J)
 100    CONTINUE
        HHBAR = HHBAR/RMXBAR
        WTMLT = 1.0/DUM
C-----NEWTON ITERATION
        DO 200 II=1,20
          HPTMP = 0.0
          DHDT  = 0.0
          K     = 2
          IF(TQ1.GT.1000.0) K = 1
          IF(TQ1.GE.300.0.AND.TQ1.LE.5000.0) THEN
C-------FOR 300 <= TQ1 <= 5000
          DO 220 J=1,NGAS
            XN1 = ALPHA(J)*WTMOLE(J)*WTMLT
            CP1 = HF(1,K,J)+(HF(2,K,J)+(HF(3,K,J)+(HF(4,K,J)
     &      +HF(5,K,J)*TQ1)*TQ1)*TQ1)*TQ1
            HP1 =(HF(1,K,J)+(.5*HF(2,K,J)+(HF(3,K,J)/3.+(.25
     &      *HF(4,K,J)+.20*HF(5,K,J)*TQ1)*TQ1)*TQ1)*TQ1)*TQ1+HF(6,K,J)
            HPTMP = HPTMP+XN1*HP1
            DHDT  = DHDT +XN1*CP1
 220      CONTINUE
          ELSE
C-------FOR 5000 < TQ1 OR 300 > TQ1
          TQQ = 5000.0
          IF(TQ1.LT.300.0) TQQ = 300.0
          DO 221 J=1,NGAS
            XN1 = ALPHA(J)*WTMOLE(J)*WTMLT
            CP1 = HF(1,K,J)+(HF(2,K,J)+(HF(3,K,J)+(HF(4,K,J)
     &      +HF(5,K,J)*TQQ)*TQQ)*TQQ)*TQQ
            HP1 = (HF(1,K,J)+(.5*HF(2,K,J)+(HF(3,K,J)/3.+(.25
     &      *HF(4,K,J)+.20*HF(5,K,J)*TQQ)*TQQ)*TQQ)*TQQ)*TQQ+HF(6,K,J)
     &      +CP1*(TQ1-TQQ)
            HPTMP = HPTMP+XN1*HP1
            DHDT  = DHDT +XN1*CP1
 221      CONTINUE
          ENDIF
C-----NEWTON EXTRAPOLATION
          HHTMP = HPTMP-HHBAR
C-GC1
c         DHHR = HHTMP/AMAX1(1.E-15,ABS(HHBAR))
c         IF(ABS(DHHR) .LE. 1.E-5) GO TO 299
C-GC2
          IF(ABS(HHTMP).LE.1.0E-10) GO TO 299
          IF(ABS(DHDT).LE.1.0E-20) GO TO 299
          TQ1   = TQ1-HHTMP/DHDT
 200    CONTINUE
 299    TMTMP = AMAX1(10.00,TQ1)
C-END-HEAT2B----------------------------------------------------------
      RETURN
      END
