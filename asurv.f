C
C      **********************************************************************
C      *********************** SUBROUTINE AARRAY ****************************
C      **********************************************************************
C
       SUBROUTINE AARRAY(Z,IND,ISTA,IS,NTOT,NDAT,MVAR,NG1,NG2,XY,
     +            ID1,ID2,ISAVE)
C
C*******           ISIGN,IFULL IS ADDED ON "COMMON' STATEMENT               *
C      *                                                                    *
C      *     INPUT       Z(I,J)  : DATA TO BE TESTED                        *
C      *                 IND(I,J): INDICATOR OF CENSORING                   *
C      *                 ISTA(I) : INDICATOR OF GROUP                       *
C      *                   IS    : IS-TH SUB-DATA SET                       *
C      *                   NG1   : INDICATOR OF THE FIRST GROUP             *
C      *                   NG2   : INDICATOR OF THE SECOND GROUP            *
C      *                   NTOT  : TOTAL NUMBER OF DATA POINTS              *
C      *                   LL    : INDICATOR OF OUTPUT FILE                 *
C      *                  IPR    : INDICATOR FOR PRINTING                   *
C      *                                                                    *
C      *     OUTPUT         N    : NTOT                                     *
C      *                    N1   : NUMBER OF DATA POINTS IN GROUP 1         *
C      *                    N2   : NUMBER OF DATA POINTS IN GROUP 2         *
C      *                   NCEN  : NUMBER OF CENSORED DATA POINTS           *
C      *                   ISIGN : INDICATOR OF LOWER/UPPER LIMITS          *
C      *                                                                    *
C      *     PUT ALL OBS. IN ARRAY XY AND FORM ARRAYS ID1 AND ID2           *
C      *           ID1(I)=0  : ITH OBS. IS UNCENSORED                       *
C      *                  1  : ITH OBS. IS CENSORED                         *
C      *           ID2(I)=J  : ITH OBS. IS FROM ITH SAMPLE, J=1,2           *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *           SORT2                                                    *
C
C*******     ALTHOUGH THIS SUBROUTINE HAS THE SAME NAME AS A PROGRAM FROM   *
C*******     "STATISTICAL METHODS FOR SURVIVAL DATA ANALYSIS" BY ELISA T.   *
C*******     LEE, 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA),        *
C*******     IT IS DIFFERENT EXCEPT IN THE GENERAL PURPOSE.                 *
C*******     ID1(I) IS ASSIGNED IN THE OPPOSITE WAY SO THAT THE PPROGRAM    *
C*******     CAN USE THE DATA SETS WHICH ARE MADE FOR OTHER PROGRAMS.       *
C

        IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

        DIMENSION Z(MVAR,NDAT),IND(MVAR,NDAT),ISTA(NTOT)
        DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT),ISAVE(NTOT)
        COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
C
        IU=0
        NCEN=0
        NCOMP=0 
        N1=0
        N2=0
        ISIGN=1
C
C      *    FIND THE CENSORSHIP OF THE DATA SET. -1 FOR UPPER LIMITS       *
C      *    AND 1 FOR LOWER LIMITS                                         *
C
       DO 100 I=1,NTOT
          ISAVE(I) = 0
          IF(IND(IS,I) .EQ. 0) GOTO 100
            ISIGN=IND(IS,I)/IABS(IND(IS,I))
            ISAVE(I) = ISIGN
  100  CONTINUE

C      * CHECK WHETHER THE UPPER AND LOWER LIMITS ARE MIXED IN THE SAME   *
C      * VARIABLE. IF SO, THE PROGRAM IS TERMINATED.                      *
C      * THIS TEST WAS ADDED.                                             *
C
       DO 110 I = 1, NTOT
          IF(ISAVE(I) .EQ. 0) GOTO 110
          IF(ISAVE(I) .NE. ISIGN) THEN
          PRINT *
             PRINT *,'YOU CANNOT HAVE BOTH UPPER AND LOWER LIMITS'
             PRINT *,'IN ONE VARIABLE AT THE SAME TIME.'
             PRINT *,'PLEASE CHECK YOUR DATA.'
             PRINT *,'THE PROGRAM HAS BEEN TERMINATED.'
             PRINT *
             STOP
          ENDIF
  110  CONTINUE
C     
C      *        COUNT NUMBER OF DATA POINTS IN THE TWO SUBSAMPLES        *
C

        DO 400 I = 1, NTOT
        IF((ISTA(I) .EQ. NG1) .OR. (ISTA(I) .EQ. NG2)) THEN
           NCOMP = NCOMP + 1
           XY(NCOMP) = ISIGN*Z(IS,I)
           IF(ISTA(I) .EQ. NG1) ID2(NCOMP) = 1
           IF(ISTA(I) .EQ. NG2) ID2(NCOMP) = 2
           IF(IABS(IND(IS,I)) .NE. 1) THEN
              ID1(NCOMP) = 0
              IU = IU + 1
              IF(ID2(NCOMP) .EQ. 1) N1 = N1 + 1
              IF(ID2(NCOMP) .EQ. 2) N2 = N2 + 1
           ELSE
              ID1(NCOMP) = 1
              NCEN = NCEN + 1
              IF(ID2(NCOMP) .EQ. 1) N1 = N1 + 1
              IF(ID2(NCOMP) .EQ. 2) N2 = N2 + 1
        ENDIF
        ENDIF
  400   CONTINUE

        CALL SORT2(XY, ID1, ID2, NTOT)

        RETURN
        END


C
C      **********************************************************************
C      ********************* FUNCTION AGAUSS  *******************************
C      **********************************************************************
C
       FUNCTION AGAUSS(Z)
C
C      *        EVALUATES THE INTEGRAL OF THE GAUSSIAN PROBABILITY FUNCTION *
C      *        OBTAINED FROM PROGRAM 3-5 ON P. 35 OF "DATA REDUCTION AND   *
C      *        ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", P. R. BEVINGTON, *
C      *        1969, McGRAW HILL, (NY:NY).                                 *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C     
C
       Z=DABS(Z)
       AGAUSS=1.0
C
C      *         IF Z>5.0, USE APPROXIMATION FOR PROB TO AVOID ERROR        *
C
       IF(Z.LE.5.0) THEN
          DENOM=1.0
          IF(Z .GT. 0.0) THEN
             TERM=0.7071067812D00*Z
             SUM=TERM
             Y2=(Z**2)/2.0
   31        DENOM=DENOM+2.0
             TERM=TERM*(Y2*2.0/DENOM)
             SUM=SUM+TERM
             IF(TERM/SUM-1.0E-10 .GT. 0.0) THEN
                GOTO 31
             ELSE
                AGAUSS=1.128379167D00*SUM*DEXP(-Y2)
             ENDIF
          ELSE
             AGAUSS = 0.0
          ENDIF
       ENDIF

       RETURN
       END

C
C*************************************************************************
C********************* SUBROUTINE AKRANK *********************************
C*************************************************************************
C
C
       SUBROUTINE AKRANK(IND, X, NTOT, IP, R, MVAR,ZU, ZC,
     +                 PL, F, V, FMASS, ITEMP, PTEMP,Z1,
     +                 WRK1,WRK2,WRK3,DWRK1,IWRK1,SWRK1) 
C
C     *   THIS SUBROUTINE COMPUTES AKRITAS' RANK                         *
C     *                                                                  *
C     *   REFERENCE                                                      *
C     *         PENN STATE UNIVERSITY, DEPARTMENT OF STATISTICS,         *
C     *         TECHNICAL REPORTS AND PREPRINTS SERIES, NUMBER 87,       *
C     *         "ALIGNED RANK TESTS FOR REGRESSION WITH CENSORED DATA",  *
C     *         MICHAEL G. AKRITAS, SEPTEMBER 1989                       *
C     *   INPUT                                                          *
C     *         IND      : INDICATOR OF CENSORSHIP                       *
C     *         X        : VARIABLE                                      *
C     *         NTOT     : TOTAL NUMBER OF DATA POINTS                   *
C     *                                                                  *
C     *   OUTPUT                                                         *
C     *         R        : RANK                                          *
C     *         PL       : PL ESTIMATOR                                  *
C     *         F        : 1.0 - PL  (DISTRIBUTION FUNCTION)             *
C     *                                                                  *
C     *   OTHER VARIABLES                                                *
C     *         IP       : INDEX OF VARIABLE BEING RANKED                *
C     *         MVAR     : NUMBER OF VARIABLES                           *
C     *         ZU       : DETECTED DATA                                 *
C     *         ZC       : CENSORED DATA                                 *
C     *         FMASS    : JUMPS IN PL ESTIMATOR                         *
C     *         IU       : NUMBER OF DETECTIONS                          *
C     *         IC       : NUMBER OF CENSORED DATA POINTS                *
C     *         PTEMP    : TEMPORARY STORAGE OF PL ESTIMATOR             *
C     *                                                                  *
C     *   SUBROUTINES                                                    *
C     *         XVAR, PLESTM, SORT1                                      *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(MVAR, NTOT), X(MVAR, NTOT), ZU(NTOT), ZC(NTOT)
       DIMENSION PL(NTOT), R(MVAR, NTOT), F(NTOT), V(NTOT), FMASS(NTOT)
       DIMENSION Z1(MVAR, NTOT), ITEMP(NTOT), PTEMP(NTOT)
       DIMENSION IWRK1(NTOT),DWRK1(MVAR,NTOT),SWRK1(MVAR)
       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT)

C
C     *     CALL SUBROUTINE XVAR : DISTINGUISH DETECTIONS AND CENSORED   *
C     *     DATA POINTS.                                                 *
C

       CALL XVAR(IND,X,IP,NTOT,ISIGN,ZU,ZC,IU,IC,IWRK1,WRK1,WRK2,WRK3,
     +           DWRK1,SWRK1,LTOT,MVAR,ITEMP)

C
C     *     CALL PLESTM : PL ESTIMATOR COMPUTATION                       *
C
C
       DO 5 I = 1, NTOT
          ITEMP(I) = 0
          Z1(1, I) = 0.0
          IWRK1(I) = IND(IP,I)
    5  CONTINUE
    
       IF(IU .EQ. 0) THEN
          WRITE(6,3)
    3     FORMAT('NO DETECTIONS: PROGRAM IS TERMINATED')
          STOP
       ENDIF
       
       CALL SORT1(IWRK1,Z1,ZU,IU,1,ITEMP,SWRK1,MVAR)
C
       IF(IC .NE. 0) CALL SORT1(IWRK1,Z1,ZC,IC,1,ITEMP,SWRK1,MVAR)

C
C  >>>>  Bug fixed Sept. 1996.  NCH was missing from following line  <<<<

       CALL PLESTM(ZU, ZC, IU, IC, PL, V, NTOT,SMEAN,SIGM,ICH,NCH,IWRK1)
       
C
C     *  IF THE DATA CONTAINS CENSORED DATA POINTS, THE PRODUCT LIMIT    *
C     *  ESTIMATOR MUST BE ADJUSTED TO INCLUDE CENSORED DATA POINTS.     *
C
       IF(IC .NE. 0) THEN
       
C     *   IF THE DATA HAS UPPER LIMITS, FIRST THE PRODUCT LIMIT ESTIMATOR*
C     *   MUST BE ADJUSTED.                                              *

          IF(ISIGN .LT. 0) THEN
             FMASS(1) = 1.0 - PL(1)
             DO 10 I = 2, IU
                FMASS(I) = PL(I-1)-PL(I)
   10        CONTINUE

             J = IU/2
             DO 20 I = 1, J
                FTEMP=FMASS(I)
                FMASS(I)=FMASS(IU-I+1)
                FMASS(IU-I+1)=FTEMP
   20        CONTINUE

             DO 40 I = 1, IU
                PTEMP(I) = 1.0
                DO 30 J = 1, I
                   PTEMP(I)=PTEMP(I)-FMASS(J)
   30           CONTINUE
   40        CONTINUE
          ELSE
             DO 50 I = 1, IU
                PTEMP(I)=PL(I)
   50        CONTINUE
          ENDIF

C      *  NOW, PRODUCT LIMIT ESTIMATOR VALUES ARE ASSIGNED TO CENSORED  *
C      *  DATA POINTS.                                                  *

          IF(IND(IP,1) .EQ. 0) THEN
             PL(1) = PTEMP(1)
             J = 1
          ELSE
             PL(1) = 1.0
             J = 0
          ENDIF
       
          DO 60 I = 2, NTOT
             IF(IND(IP, I) .EQ. 0) THEN
                J = J + 1
                PL(I) = PTEMP(J)
             ELSE 
                PL(I) = PL(I-1)
             ENDIF
   60     CONTINUE
   
       ENDIF

C     *  THE PRODUCT LIMIT ESTIMATE IS NOW USED TO ESTIMATE THE       * 
C     *  DISTRIBUTION FUNCTION (F) AT ALL POINTS.                     *

       DO 65 I = 1, NTOT
          F(I) = 1.0 - PL(I)
   65  CONTINUE

C
C     *   COMPUTE HERE AKRITAS' RANK USING F-VALUES                 *
C
       DO 90 I = 1, NTOT
         IF(IND(IP, I) .EQ. 0) THEN
            R(IP, I) = REAL(NTOT)*F(I)
         ELSEIF(IND(IP, I) .GT. 0) THEN
            R(IP, I) = REAL(NTOT)*(0.5 + 0.5*F(I))
         ELSE
            R(IP, I) = NTOT*(0.5*F(I))
         ENDIF
   90  CONTINUE

       RETURN
       END

C
C      **********************************************************************
C      ******************** SUBROUTINE ARISK  *******************************
C      **********************************************************************
C
       SUBROUTINE ARISK(R,XM,X,E1,NG,H,XY,ID1,NTOT)
C
C  
C      *       THIS SUBROUTINE COMPUTES THE FOLLOWING FOUR                  *
C      *       ARRAYS FOR SUBROUTINE COX, LRANK, AND PWLCXN.                *
C      *         R(I) : NO. OF OBSERVATIONS IN RISK SET AT THE              *
C      *                I-TH DISTINCT FAILURE TIME.                         *
C      *        XM(I) : MULTIPLICITY OF THE I-TH DISTINCT                   *
C      *                FAILURE TIME.                                       *
C      *        E1(I) : XM(I)/R(I)                                          *
C      *         H(I) : KAPLAN AND MEIER'S ESTIMATES OF THE                 *
C      *                SURVIVOR FUNCTION                                   *
C      *                                                                    *
C      *         X(I) : THE ARRAY OF DISTINCT FAILURE TIMES                 *
C      *         NG   : NO OF X                                             *
C      *  THIS SUBROUTINE IS OBTAINED FROM ELISA T. LEE, "STATISTICAL       *
C      *  METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING      *
C      *  PUBLICATIONS (BELMONT:CA); BUT HAS BEEN SIGNIFICANTLY MODIFIED.   *
C      *                                                                    *
C

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION R(NTOT),XM(NTOT),X(NTOT),H(NTOT),E1(NTOT)
       DIMENSION XY(NTOT),ID1(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
C
       L=1
       I=1
       R(L)=REAL(NCOMP)
C
C      *             COMPUTE RISK SETS, AND OTHER QUANTITIES                *
C
  24   IF(ID1(I).NE.0) THEN
          R(L)=R(L)-1.0
          I=I+1
          GOTO 24
       ENDIF
  25   XM(L)=1.0
       XNC=0.0
       TEMP=XY(I)
       X(L)=TEMP

  21   IF(I.NE.NCOMP) THEN
          I=I+1
C
          IF(ID1(I).NE.1) THEN
             IF(TEMP.NE.XY(I)) GOTO 20
             XM(L)=XM(L)+1.0
             GOTO 21
          ENDIF

  26      XNC=XNC+1.0
          X(L)=TEMP
          GOTO 21

  20      L=L+1
          R(L)=R(L-1)-XM(L-1)-XNC
          GOTO 25
       ENDIF

  23   X(L)=TEMP
       NG=L
C    
C      *          COMPUTE KM ESTIMATOR                                      *

       DO 30 I=1,NG
          E1(I)=XM(I)/R(I)
  30   CONTINUE

       H(1)=1.0
       NG1=NG+1

       DO 31 I=2,NG1
          H(I)=H(I-1)*(1.0-E1(I-1))
  31   CONTINUE

       RETURN
       END


C
C
C      *  ASURV:   SURVIVAL ANALYSIS PACKAGE FOR ASTRONOMERS                *
C      *                                                                    *
C      *  DEVELOPED BY:        TAKASHI ISOBE                                *
C      *                 CENTER FOR SPACE RESEARCH                          * 
C      *           THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY                *
C      *                                                                    *
C      *                      MICHAEL LAVALLEY                              *
C      *                  DEPARTMENT OF STATISTICS                          *
C      *               THE PENSYLVANIA STATE UNIVERSITY                     *
C      *       330A CLASSROOM BUILDING, UNIVERSITY PARK PA 16802            *
C      *                 INTERNET: MLV@STAT.PSU.EDU                         *
C      *                                                                    *
C      *                      ERIC FEIGELSON                                *
C      *             DEPARTMENT OF ASTRONOMY AND ASTROPHYSICS               * 
C      *                THE PENSYLVANIA STATE UNIVERSITY                    *
C      *              525 DAVEY LAB. UNIVERSITY PARK PA 16802               *
C      *                                                                    *
C      *  REV. 1.2  SECOND UPDATE   SUMMER 1992                             *
C      *                                                                    *
C      *     THIS PACKAGE IS WRITTEN TO PROVIDE SEVERAL                     *
C      *  SURVIVAL ANALYSIS METHODS WHICH ARE USEFUL IN ANALYZING           *
C      *  ASTRONOMICAL DATA. SURVIVAL ANALYSIS IS A GROUP OF STATISTICAL    *
C      *  METHODS WHICH TREAT PROBLEMS WITH CENSORED DATA (UPPER OR LOWER   *
C      *  LIMITS). THIS PACKAGE INCLUDES SOME TECHNIQUES DEVELOPED IN       *
C      *  IN OTHER FIELDS (E.G. ECONOMICS, ACTUARIAL SCIENCE, RELIABILITY   *
C      *  MATHEMATICS), AND A FEW METHODS DEVELOPED BY ASTRONOMERS.         *
C      *                                                                    *
C      *   THE METHODS PROVIDED IN THIS PACKAGE ARE :                       *
C      *                                                                    *
C      *   UNIVARIATE DISTRIBUTION :  KAPLAN-MEIER ESTIMATOR                *
C      *   TWO-SAMPLE TESTS        :  GEHAN TEST                            *
C      *                              LOGRANK TEST                          *
C      *                              PETO AND PETO TEST                    *
C      *                              PETO AND PRENTICE TEST                *
C      *   CORRELATION TESTS       :  COX PROPORTIONAL HAZARDS MODEL        *
C      *                              GENERALIZED KENDALL'S TAU (BHK METHOD)*
C      *                              GENERALIZED SPEARMAN'S RHO            *
C      *                                   (AKRITAS' METHOD)                *
C      *   LINEAR REGRESSIONS      :  EM ALGORITHM WITH NORMAL DISTRIBUTION *
C      *                              BUCKLEY-JAMES METHOD                  *
C      *                              TWO-DIMENSIONAL KAPLAN-MEIER          *
C      *                                  REGRESSION FOR DUAL-CENSORED DATA *
C      *                                                                    *
C      *                                                                    *
C      *   INPUTS                                                           *
C      *                                                                    *
C      *       IS0     :  IF 1 : UNIVARIATE PROBLEM                         *
C      *                     2 : CORRELATION/REGRESSION PROBLEM             *
C      *                     3 : EXIT                                       *
C      *                                                                    *
C      *   SUBROUTINES DATA1, UNIVAR, BIVAR                                 *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 BLANK
C       OPEN(6,CARRIAGECONTROL='LIST',STATUS='OLD')
C
C
C
       PRINT *
       PRINT *,'    ***************************************************'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                WELCOME TO ASURV                 *'
       PRINT *,'    *           SURVIVAL ANALYSIS PACKAGE             *' 
       PRINT *,'    *                FOR ASTRONOMERS                  *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                  DEVELOPED BY:                  *'
       PRINT *,'    *                  TAKASHI ISOBE                  *' 
       PRINT *,'    *         (CENTER FOR SPACE RESEARCH, MIT)        *'
       PRINT *,'    *                 MICHAEL LAVALLEY                *'
       PRINT *,'    *         (DEPT. OF STATISTICS, PENN STATE)       *'
       PRINT *,'    *                  ERIC FEIGELSON                 *'
       PRINT *,'    * (DEPT. OF ASTRONOMY & ASTROPHYSICS, PENN STATE) *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *              REV 1.2  SUMMER 1992               *'
       PRINT *,'    ***************************************************'
       PRINT *
       PRINT *
       PRINT *
       PRINT *
       PRINT *,'              (CARRIAGE RETURN TO CONTINUE) '
       READ(5,50) BLANK
   50  FORMAT(A1)
       PRINT *
C
C      *           START CONVERSATION WITH THE USER                         *
C
       PRINT *
       PRINT *
       PRINT *
  100  PRINT *,'                          MENU  '
       PRINT *
       PRINT *
       PRINT *,'       UNIVARIATE DATA           BIVARIATE DATA '
       PRINT *
       PRINT *
       PRINT *,'     DISTRIBUTION FUNCTION       CORRELATION '
       PRINT *,'   1 KAPLAN-MEIER ESTIMATOR    1 COX REGRESSION '
       PRINT *,'                               2 GEN. KENDALL TAU'
       PRINT *,'                               3 GEN. SPEARMAN RHO'
       PRINT *
       PRINT *
       PRINT *,'     TWO-SAMPLE TESTS            LINEAR REGRESSION '
       PRINT *,'   1 GEHAN TESTS               1 EM ALGORITHM WITH  '
       PRINT *,'   2 LOGRANK TEST                 GAUSSIAN RESIDUALS ' 
       PRINT *,'   3 PETO AND PETO TEST        2 BUCKLEY-JAMES METHOD ' 
       PRINT *,'   4 PETO AND PRENTICE TEST       WITH KM RESIDUALS '
       PRINT *,'                               3 SCHMITT METHOD FOR '
       PRINT *,'                                  DUAL CENSORED DATA ' 
       PRINT *
       PRINT *
       PRINT *
       PRINT *,'            (CARRIAGE RETURN TO CONTINUE) '
       READ(5,50) BLANK
C
       PRINT *
C
C      *  CHOICE : UNIVARIATE PROBLEM OR CORRELATION/REGRESSION PROBLEM     *
C
       PRINT *
       PRINT *,'    SELECT DATA TYPE: ' 
       PRINT *,'     1 UNIVARIATE DATA '
       PRINT *,'     2 BIVARIATE DATA ' 
       PRINT *,'     3 EXIT '
  200  WRITE(6,210)
  210  FORMAT(' CHOICE ? ')
C  210  FORMAT('          CHOICE ? ',$)
C
       CALL DATA1(IS0)
C
       IF((IS0.EQ.1).OR.(IS0.EQ.2).OR.(IS0.EQ.3)) GOTO 300
       PRINT *,'PLEASE TYPE ONCE MORE'
       GOTO 200
C
  300  IBACK=0
       IF(IS0.EQ.1) CALL UNIVAR(IBACK)
       IF(IS0.EQ.2) CALL BIVAR(IBACK)
       IF(IS0.EQ.3) STOP
C
       IF(IBACK.EQ.1) GOTO 100
       STOP
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE BHK  *******************************
C      **********************************************************************
C
       SUBROUTINE BHK(IND,XX,YY,NTOT,OUTPUT,X,Y,IAA,IBB,IP,MVAR)
C
C      *             GENERALIZED KENDALL'S TAU CORRELATION COEFFICIENT      *
C      *                         FOR CENSORED DATA                          *
C
C      *      THIS PROGRAM COMPUTES KENDALL'S TAU FOR BIVARIATE DATA        *
C      *      SETS. THE DATA SETS CAN CONTAIN CENSORED POINTS IN THE        *
C      *      INDEPENDENT VARIABLE AND/OR THE DEPENDENT VARIABLE.           *
C      *      ALTHOUGH THIS PROGRAM GIVES ANSWERS FOR DATA SETS WHICH       *
C      *      CONTAIN TIES, IT MAY NOT BE ACCURATE.                         *
C      *    PARAMETERS :                                                    *
C      *     INPUT                                                          *
C      *     NTOT          : NUMBER OF OBSERVATIONS                         *
C      *     XX(1,I)       : INDEPENDENT PARAMETER OF I-TH OBSERVATION      *
C      *     YY(I)         : DEPENDENT PARAMETER OF I-TH OBSERVATION        *
C      *     IND(I)        : INDICATOR OF CENSORED STATUS                   *
C      *      EACH POINT MUST BE SPECIFIED ITS CENSORED STATUS :            *
C      *        FOR THE LOWER LIMITS                                        *
C      *                0   :   DETECTED POINT                              *
C      *                1   :   ONLY DEPENDENT VARIABLE IS LOWER LIMIT      *
C      *                2   :   ONLY INDEPENDENT VARIABLE IS LOWER LIMIT    *
C      *                3   :   BOTH VARIABLES ARE LOWER LIMIT              *
C      *                4   :   INDEPENDENT VARIABLE IS LOWER LIMIT AND     *
C      *                        DEPENDENT VARIABLE IS UPPER LIMIT           *
C      *      FOR THE UPPER LIMITS, CHANGE THE SIGN OF ABOVE INDICATORS.    *
C      *                                                                    *
C      *     WORK                                                           *
C      *     X(I)           : =XX(1,I)                                      *
C      *     Y(I)           : =YY(I)                                        *
C      *     IP(I)          : =IND(I)                                       *
C      *     IAA(I)         : CONCORDANCE INFORMATION FOR X                 *
C      *     IBB(I)         : CONCORDANCE INFORMATION FOR Y                 *
C      *     OUTPUT                                                         *
C      *     PROB          : SIGNIFICANCE LEVEL FOR THE HYPOTHESIS THAT     *
C      *                     X AND Y ARE NOT CORRELATED UNDER THE           *
C      *                     GAUSSIAN DISTRIBUTION                          *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *               CENS, COEFF                                          *
C      *                                                                    *
C      *      REF. BROWN, HOLLANDER, AND KORWAR 1974, IN RELIABILITY        *
C      *           AND BIOMETRY P.327, EQNS 1 TO 8, PROSCHAN AND            *
C      *           SERFLING EDS (SIAM)                                      *
C
C      *     NOTE:  THIS PROGRAM IS QUITE CPU INTENSIVE FOR LARGE DATA      *
C      *            SETS (MORE THAN A FEW HUNDRED POINTS).                  *
C      *                                                                    *
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION XX(MVAR,NTOT),YY(NTOT),IND(NTOT),X(NTOT),Y(NTOT)
       DIMENSION IP(NTOT),IAA(NTOT),IBB(NTOT)
       CHARACTER*9 OUTPUT 
C
       SIS =0.0
       ASUM =0.0
       BSUM =0.0
       AASUM=0.0
       BBSUM=0.0
C
C      *       SUBSTITUE XX AND YY TO X AND Y SO THAT THE ORIGINAL VALUES   *
C      *       WON'T BE CHANGED.                                            *
C
       DO 90 I=1,NTOT
          X(I)  = XX(1,I)
          Y(I)  = YY(I)
          IP(I) = IND(I)
   90  CONTINUE
C
C
C      *      THE SUBROUTINE CENS ADDS OR SUBTRACTS A SMALL NUMBER          *
C      *      FROM EACH CENSORED POINT SO THAT NO TIES WITH DETECTED        *
C      *      POINTS OCCUR.                                                 *
C
C
       CALL CENS(X,Y,IP,NTOT)
C
C
C      *      START MAKING INFORMATION FOR CONCORDANCE                      *
C
C
       DO 1900 I=1,NTOT
C
C      *      INFORMATION OF CONCORDANCE  FOR THE INDEPENDENT VAR.          *
C
          IA=2
          IB=3
          IC=4
          ID=-2
          IE=-3
          IG=-4
          IH=1
          IJ=-1
C
C      *       SUBROUTINE WHICH FINDS CONCORDANCE INFORMATION               *
C
          CALL COEFF(I,X,IP,NTOT,IAA,IA,IB,IC,ID,IE,IG,IH,IJ)
C
C      *       INFORMATION OF CONCORDANCE FOR THE DEPENDENT VAR.            *
C
          IA=1
          IB=3
          IC=-4
          ID=-1
          IE=-3
          IG=4
          IH=2
          IJ=-2

          CALL COEFF(I,Y,IP,NTOT,IBB,IA,IB,IC,ID,IE,IG,IH,IJ)
C
C      *        START COMPUTING QUANTITIES IS, IASUM, IBSUM,                *
C      *        IAASUM, AND IBBSUM.                                         * 
C
          DO 1800 J=1,NTOT
             IF((IAA(J).EQ.0).AND.(IBB(J).EQ.0)) GOTO 1800
             SIS=SIS+IAA(J)*IBB(J)
             ASUM=ASUM+IAA(J)**2
             BSUM=BSUM+IBB(J)**2

 1650        DO 1700 K=1,NTOT
                IF(IAA(J).NE.0) THEN
                   IF(IAA(K).NE.0) THEN
                      AASUM=AASUM+IAA(J)*IAA(K)
                   ENDIF
                ENDIF
 1670           IF(IBB(J).NE.0) THEN
                   IF(IBB(K).NE.0) THEN
                      BBSUM=BBSUM+IBB(J)*IBB(K)
                   ENDIF
                ENDIF
 1700        CONTINUE
 1800     CONTINUE
 1900  CONTINUE
C
C      *    NOW COMPUTE THE STATISTIC AND THE PROBABILITY                   *
C 
       D1=REAL(NTOT*(NTOT-1))
       D2=REAL(D1*(NTOT-2))
       ALP=2.0*(ASUM*BSUM)/D1
       GAM=4.0*((AASUM-ASUM)*(BBSUM-BSUM))/D2
       VAR=ALP+GAM
       SIGMA=DSQRT(VAR)
       Z=SIS/SIGMA
       PROB=1.0-AGAUSS(Z)
C
       IF(OUTPUT.EQ.'         ') THEN
          WRITE(6,2030) 
          WRITE(6,2003)
          WRITE(6,2030)
          WRITE(6,2005) Z
          WRITE(6,2007) PROB
          WRITE(6,2030) 
       ELSE
          WRITE(60,2030) 
          WRITE(60,2003)
          WRITE(60,2030)
          WRITE(60,2005) Z
          WRITE(60,2007) PROB
          WRITE(60,2030) 
       ENDIF
 2003  FORMAT(5X,'CORRELATION TEST BY GENERALIZED KENDALL`S TAU')
 2005  FORMAT(7X,'Z-VALUE      =',F12.3)
 2007  FORMAT(7X,'PROBABILITY  =',F13.4,/,
     +    ' (PROBABILITY THAT A CORRELATION IS NOT PRESENT)')
 2030  FORMAT('      ')

       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE BIN  *******************************
C      **********************************************************************
C
       SUBROUTINE BIN(NTOT,MX,MY,ISKIP,ICENS,DELX,DELY,XORG,YORG,MM,
     +                M1,M2,M3,M4,M5,M6,M7,M8,INDEX,LP,XT,YT,Z,SWRK1,
     +                X,Y,NP,XB,YB,F,N,N1,N2,N3,N4,N5,N6,N7,N8,IB,MVAR)

C
C
C      *                                                                    *
C      *    THIS SUBROUTINE DOES BINNING AND CHANGES CENSORED POINTS        *
C      *    WHICH DO NOT HAVE DETECTED POINTS ABOVE (OR BELOW)              *
C      *    TO DETECTED POINTS.                                             *
C      *                                                                    *
C      *             WARNING   WARNING   WARNING   WARNING                  *
C      *                                                                    *
C      *    THE USER SHOULD BE WARNED THAT THIS SUBROUTINE ACTUALLY         *
C      *    CHANGES THE DATA!!  FIRST, IT REDEFINES SOME LIMITS TO          *
C      *    DETECTIONS.  IF THE BINS ARE CHOSEN TO BE TOO NARROW, THEN      *
C      *    VIRTUALLY ALL LIMITS COULD BE CHANGED.  SECOND, IT PUSHES       *
C      *    EACH LIMIT INTO THE ADJACENT BIN.  IF THE BINS ARE CHOSEN TO    *
C      *    TO BE TOO WIDE, THIS SUBSTANTIALLY ALTERS THE MEASURED VALUES.  *
C      *    THUS, THE USER MUST TREAD A FINE LINE IN CHOSING BIN SIZES.     *
C      *                                                                    * 
C      *                                                                    *
C      *    INPUT                                                           *
C      *           X(I)  : INDEPENDENT VARIABLE                             *
C      *           Y(I)  : DEPENDENT VARIABLE                               *
C      *          NP(I)  : INDICATOR OF CENSORING                           *
C      *          NTOT   : TOTAL NUMBER OF DATA                             *
C      *            MX   : NUMBER OF BINS IN X                              *
C      *            MY   : NUMBER OF BINS IN Y                              *
C      *          ISKIP  : INDICATOR OF BINNING PROCESS                     *
C      *          ICENS  : CENSORING STATUS OF THE DATA SET                 *
C      *      IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :                *
C      *          DELX   : BIN SIZE OF X AXIS                               *
C      *          DELY   : BIN SIZE OF Y AXIS                               *
C      *          XORIG  : ORIGIN OF X                                      *
C      *          YORIG  : ORIGIN OF Y                                      *
C      *                                                                    *
C      *     WORK                                                           *
C      *           YT(I) : COPY OF Y(I) FOR SORTING PROGRAM.                *
C      *            M1   : # OF Y LOWER LIMITS CHANGED TO DETECTIONS        *
C      *            M2   : # OF X LOWER LIMITS CHANGED TO DETECTIONS        *
C      *            M3   : # OF DOUBLE LOWER LIMITS CHANGED TO              *
C      *                   DETECTIONS                                       *
C      *            M4   : # OF Y LOWER , X UPPER LIMITS CHANGED TO         *
C      *                   DETECTIONS                                       *
C      *            M5   : # OF Y UPPER LIMITS CHANGED TO DETECTIONS        *
C      *            M6   : # OF X LOWER LIMITS CHANGED TO DETECTIONS        *
C      *            M7   : # OF DOUBLE UPPER LIMITS CHANGED TO              *
C      *                   DETECTIONS                                       *
C      *            M8   : # OF Y UPPER , X LOWER LIMITS CHANGED TO         *
C      *                   DETECTIONS                                       *
C      *            NC1, NC2,...,NC8 : # OF CENSORED POINTS. SEE THE        *
C      *                  MAIN PROGRAM FOR THE DEFINITIONS                  *
C      *            IB   : DIMENSION SIZE OF BINS                           *
C      *                                                                    *
C      *    OUTPUT                                                          *
C      *           F(I,J): INITIAL GUESS OF THE PROBABILITY OF THE          *
C      *                   BIN(I,J)                                         *
C      *           N(I,J): NUMBER OF DETECTED POINTS IN THE BIN(I,J)        *
C      *          N1(I,J): NUMBER OF Y LOWER LIMITS IN THE BIN(I,J)         *
C      *          N2(I,J): NUMBER OF X LOWER LIMITS IN THE BIN(I,J)         *
C      *          N3(I,J): NUMBER OF DOUBLE LOWER LIMITS IN THE BIN(I,J)    *
C      *          N4(I,J): NUMBER OF Y LOWER, X UPPER LIMITS IN THE         *
C      *                   BIN(I,J)                                         *
C      *          N5(I,J): NUMBER OF Y UPPER LIMITS IN THE BIN(I,J)         *
C      *          N6(I,J): NUMBER OF X UPPER LIMITS IN THE BIN(I,J)         *
C      *          N7(I,J): NUMBER OF DOUBLE UPPER LIMITS IN THE BIN(I,J)    *
C      *          N8(I,J): NUMBER OF Y UPPER, X LOWER LIMITS IN THE         *
C      *                   BIN(I,J)                                         *
C      *           XB(I) : COORDINATE OF CENTER OF THE BIN IN X             *
C      *           YB(I) : COORDINATE OF CENTER OF THE BINS IN Y            *
C      *     IF ISKIP=0, THE NEXT VALUES ARE OUTPUTS  :                     *
C      *          DELX   : BIN SIZE OF X AXIS                               *
C      *          DELY   : BIN SIZE OF Y AXIS                               *
C      *          XORIG  : ORIGIN OF X                                      *
C      *          YORIG  : ORIGIN OF Y                                      *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *                   SORT1                                            *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION INDEX(NTOT),LP(NTOT),XT(NTOT),YT(NTOT),Z(MVAR,NTOT)
       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB),SWRK1(MVAR)
       DIMENSION F(IB,IB),N(IB,IB),N1(IB,IB),N2(IB,IB),N3(IB,IB)
       DIMENSION N4(IB,IB),N5(IB,IB),N6(IB,IB),N7(IB,IB),N8(IB,IB)
       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
C
C
C      *  SUBSTITUE NP, X, AND Y TO LP, XT, AND YT SO THAT THE ORIGINAL DATA*
C      *  WON'T BE CHANGED.                                                 *
C
       DO 100 J=1,NTOT
          LP(J)=NP(J)
          XT(J)=X(J)
          YT(J)=Y(J)
          Z(1,J)=1.0
  100  CONTINUE
C
C      *    CALL THE SUBROUTINE SORT1, AND FIND MIN. AND MAX. OF X AND Y.   *
C      *    IF ISKIP=0, THE ORIGIN AND BIN SIZES ARE ALREADY GIVEN          *
C
       IF(ISKIP.EQ.0) THEN
C
C      *                 SORTING X                                          *
C
          CALL SORT1(LP,Z,XT,NTOT,1,INDEX,SWRK1,MVAR)
C
C      *                 SORTING Y                                          *
C
          CALL SORT1(LP,Z,YT,NTOT,1,INDEX,SWRK1,MVAR)
C
C      *              FIND THE SIZES OF BINS                                *
C
          DELX=XT(NTOT)-XT(1)
          DELY=YT(NTOT)-YT(1)
          DELX=DELX/FLOAT(MX-2)
          DELY=DELY/FLOAT(MY-2)
C
C      *             FIND THE ORIGIN OF THE GRID                            *
C
          XORG=XT(1)-1.5*DELX
          YORG=YT(1)-1.5*DELY
       ENDIF
C
C
C      *           INITIALIZE N, N1,....,N8, AND F                          *
C
       DO 300 I=1,MX
          DO 200 J=1,MY
             N(I,J) =0
             N1(I,J)=0
             N2(I,J)=0
             N3(I,J)=0
             N4(I,J)=0
             N5(I,J)=0
             N6(I,J)=0
             N7(I,J)=0
             N8(I,J)=0
             F(I,J)=0.0
  200     CONTINUE
  300  CONTINUE
C
       DO 390 I=1,NTOT
C
C      *    FIND POSITION OF I-TH DATA POINT IN THE GRID AND COUNT          *
C      *    NUMBERS OF N,N1,N2,.....,N8.                                    *
C
          IP=INT((X(I)-XORG)/DELX)+1
          JP=INT((Y(I)-YORG)/DELY)+1

C
C      *      FOR CONVENIENCE CENSORED POINTS ARE ASSIGNED TO THE NEXT BIN  *
C
C
C      *              DETECTIONS                                            *
C
          IF(NP(I).EQ.0) THEN
             N(IP,JP)=N(IP,JP)+1
      
C
C      *              Y LOWER LIMITS                                        *
C
          ELSEIF(NP(I).EQ.1) THEN
             N1(IP,JP+1)=N1(IP,JP+1)+1
          
C
C      *              X LOWER LIMITS                                        *
C
          ELSEIF(NP(I).EQ.2) THEN
             N2(IP+1,JP)=N2(IP+1,JP)+1
          
C
C      *              DOUBLE LOWER LIMITS                                   *
C
          ELSEIF(NP(I).EQ.3) THEN
             N3(IP+1,JP+1)=N3(IP+1,JP+1)+1
          
C
C      *              Y LOWER LIMITS, X UPPER LIMITS                        *
C
          ELSEIF(NP(I).EQ.4) THEN
             N4(IP+1,JP-1)=N4(IP+1,JP-1)+1
          
C
C      *              Y UPPER LIMITS                                        *
C
          ELSEIF(NP(I).EQ.-1) THEN
             N5(IP,JP-1)=N5(IP,JP-1)+1
          
C
C      *              X UPPER LIMITS                                        *
C
          ELSEIF(NP(I).EQ.-2) THEN
             N6(IP-1,JP)=N6(IP-1,JP)+1
          
C
C      *              DOUBLE  UPPER LIMITS                                  *
C
          ELSEIF(NP(I).EQ.-3) THEN
             N7(IP-1,JP-1)=N7(IP-1,JP-1)+1
          
C
C      *              Y UPPER LIMITS, X LOWER LIMITS                        *
C
          ELSEIF(NP(I).EQ.-4) THEN
             N8(IP-1,JP+1)=N8(IP-1,JP+1)+1

       ELSE
          PRINT *,' THE CENSORSHIP INDICATOR IS NOT RECOGNIZED'
          RETURN
       ENDIF
  390  CONTINUE
C
C      *    SET THE COORDINATES OF THE EACH BIN                             *
C
       DO 410 I=1,MX
          XB(I)=XORG+DELX/2.0+DELX*(I-1)
  410  CONTINUE

       DO 420 I=1,MY
          YB(I)=YORG+DELY/2.0+DELY*(I-1)
  420  CONTINUE
C
C      *    START CHECKING THE RELATION BETWEEN CENSORED POINTS AND         *
C      *    DETECTED POINTS. IF THE CENSORED POINTS ARE LOCATED  SO         *
C      *    THAT THEY CANNOT GIVE WEIGHT TO DETECTED POINTS, THE            *
C      *    CENSORED POINTS ARE CHANGED TO DETECTIONS.                      *
C
       M1=0
       M2=0
       M3=0
       M4=0
       M5=0
       M6=0
       M7=0
       M8=0
C
C
C      *              Y LOWER LIMITS                                        *
C

       IF(NC1.NE.0) THEN
          DO 600 I=1,MX
             DO 500 J=1,MY
                JJ=MY-J+1
                IF(N1(I,JJ).NE.0) THEN
                   K=JJ
  450              IF(N(I,K).EQ.0) THEN
                      K=K+1
                      IF(K.LE.MY) GOTO 450
                      M1=M1+N1(I,JJ)
                      N(I,JJ)=N(I,JJ)+N1(I,JJ)
                      N1(I,JJ)=0
                   ENDIF
                ENDIF
  500        CONTINUE
  600     CONTINUE
       ENDIF

C
C
C      *              X LOWER LIMITS                                        *
C
       IF(NC2.NE.0) THEN
          DO 800 J=1,MY
             DO 700 I=1,MX
                II=MX-I+1
                IF(N2(II,J).NE.0) THEN
                   L=II
  650              IF(N(L,J).EQ.0) THEN
                      L=L+1
                      IF(L.LE.MX) GOTO 650
                      M2=M2+N2(II,J)
                      N(II,J)=N(II,J)+N2(II,J)
                      N2(II,J)=0
                   ENDIF
                ENDIF
  700        CONTINUE
  800     CONTINUE
       ENDIF

C
C
C      *              DOUBLE LOWER LIMITS                                   *
C
       IF(NC3.NE.0) THEN
          DO 1000 I=1,MX
             II=MX-I+1
             DO 950  J=1,MY
                JJ=MY-J+1
                IF(N3(II,JJ).NE.0) THEN
                   L=II
  850              K=JJ
  900              IF(N(II,JJ).EQ.0) THEN
                      K=K+1
                      IF(K.LE.MY) GOTO 900
                      L=L+1
                      IF(L.LE.MX) GOTO 850
                      M3=M3+N3(II,JJ)
                      N(II,JJ)=N(II,JJ)+N3(II,JJ)
                      N3(II,JJ)=0
                   ENDIF
                ENDIF
  950        CONTINUE
 1000     CONTINUE
       ENDIF

C
C
C      *              Y LOWER LIMITS, X UPPER LIMITS                        *
C
       IF(NC4.NE.0) THEN
          DO 1300 I=1,MX
             II=MX-I+1
             DO 1200 J=1,MY
                IF(N4(II,J).NE.0) THEN
                   L=II
 1050              K=J
 1100              IF(N(L,K).EQ.0) THEN
                      K=K-1
                      IF(K.GE.1) GOTO 1100
                      L=L+1
                      IF(L.LE.MX) GOTO 1050
                      M4=M4+N4(II,J)
                      N(II,J)=N(II,J)+N4(II,J)
                      N4(II,J)=0
                   ENDIF
                ENDIF
 1200        CONTINUE
 1300     CONTINUE
       ENDIF

C
C
C      *              Y UPPER LIMITS                                        *
C
       IF(NC5.NE.0) THEN
          DO 1600 I=1,MX
             DO 1500 J=1,MY
                IF(N5(I,J).NE.0) THEN
                   K=J
 1450              IF(N(I,K).EQ.0) THEN
                      K=K-1
                      IF(K.GE.1) GOTO 1450
                      M5=M5+N5(I,J)
                      N(I,J) = N(I,J) + N5(I,J)
                      N5(I,J)=0
                   ENDIF
                ENDIF
 1500        CONTINUE
 1600     CONTINUE
       ENDIF

C
C
C      *              X UPPER LIMITS                                        *
C
       IF(NC6.NE.0) THEN
          DO 1800 J=1,MY
             DO 1700 I=1,MX
                IF(N6(I,J).NE.0) THEN
                   L=I
 1650              IF(N(L,J).EQ.0) THEN
                      L=L-1
                      IF(L.GE.1) GOTO 1650
                      M6=M6+N6(I,J)
                      N(I,J)=N(I,J)+N6(I,J)
                      N6(I,J)=0
                   ENDIF
                ENDIF
 1700        CONTINUE
 1800     CONTINUE
       ENDIF

C
C
C      *              DOUBLE UPPER LIMITS                                   *
C
       IF(NC7.NE.0) THEN
          DO 2000 I=1,MX
             DO 1950 J=1,MY
                IF(N7(I,J).NE.0) THEN
                   L=I
 1850              K=J
 1900              IF(N(L,K).EQ.0) THEN
                      K=K-1
                      IF(K.GE.1) GOTO 1900
                      L=L-1
                      IF(L.GE.1) GOTO 1850
                      M7=M7+N7(I,J)
                      N(I,J)=N(I,J)+N7(I,J)
                      N7(I,J)=0
                   ENDIF
                ENDIF
 1950        CONTINUE
 2000     CONTINUE
       ENDIF

C
C
C      *              Y UPPER LIMITS, X LOWER LIMITS                        *
C
       IF(NC8.NE.0) THEN
          DO 2300 I=1,MX
             DO 2200 J=1,MY
                JJ=MY-J+1
                IF(N8(I,JJ).NE.0) THEN
                   L=I
 2050              K=JJ
 2100              IF(N(L,K).EQ.0) THEN
                      K=K+1
                      IF(K.LE.MY) GOTO 2100
                      L=L-1
                      IF(L.GE.1) GOTO 2050
                      M8=M8+N8(I,JJ)
                      N(I,JJ) = N(I,JJ)+N8(I,JJ)
                      N8(I,JJ)=0
                   ENDIF
                ENDIF
 2200        CONTINUE
 2300     CONTINUE
       ENDIF


C
       MM=M1+M2+M3+M4
C
C      *             INITIAL GUESS OF F                                     *
C
       SNT=NTOT
       DO 2440 I=1,MX
          DO 2430 J=1,MY
             IF(N(I,J).NE.0) F(I,J)=FLOAT(N(I,J))/SNT
 2430     CONTINUE
 2440  CONTINUE

       RETURN
       END



C
C      **********************************************************************
C      ******************** SUBROUTINE BJ   *********************************
C      **********************************************************************
C
       SUBROUTINE BJ(IND,X,Y,NTOT,TOL,MAX,NVAR,ND,NC,ICENS,OUTPUT,
     +               ALPHA,SIGMAA,
     +               IWRK1,IWRK2,IWRK4,IWRK5,IWRK6,IWRK7,
     +               IWRK8,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
     +               SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
     +               SWRK8,SWRK9,DWRK1,EWRK1,MVAR)
C
C
C      *    LINEAR REGRESSION WITH CENSORED DATA : BUCKLEY-JAMES METHOD     *
C      *                                                                    *
C      *        USING A NONPARAMETRIC METHOD, THIS PROGRAM CALCULATES       *
C      *    LINEAR REGRESSION COEFFICIENTS FOR DATA WHICH CONTAINS SOME     *
C      *    CENSORED OBSERVATIONS.                                          *
C      *                                                                    *
C      *    PARAMETERS :                                                    *
C      *     INPUT                                                          *
C      *     NTOT          : NUMBER OF OBSERVATIONS                         *
C      *     NVAR          : NUMBER OF INDEPENDENT VARIABLE                 *
C      *     ALPHA         : INITIAL REGRESSION COEFFICIENT ESTIMATES       *
C      *                       (PLEASE ALWAYS USE 0.0 IN THIS PROGRAM).     *
C      *     TOL           : TOLERANCE FOR CONVERGENCE  (E.G. 1.0E-05)      *
C      *     MAX           : MAXIMUM ITERATION (E.G. 20)                    *
C      *      X(J,I)       : THE MATRIX CONTAINS THE COEFF.OF J-TH          *
C      *                     LOCATION PARAMETER AND I-TH OBSERVATION        *
C      *      Y(I)         : DEPENDENT PARAMETER OF I-TH OBSERVATION        *
C      *     IND(I)        : INDICATOR OF CENSORED DATA ...                 *
C      *                     IF IND(I)= 1  : LOWER LIMIT                    *
C      *                              = 0  : UNCENSORED POINT               *
C      *                              =-1  : UPPER LIMIT                    *
C      *      OUTPUT                                                        *
C      *     ALPHA(1)      : INTERCEPT COEFFICIENT                          *
C      *     ALPHA(J)      : J>1, J-TH SLOPE COEFFICIENTS                   *
C      *     ALPHA(MPLONE) : STANDARD DEVIATION                             *
C      *     SIGMAA(J)     : STANDARD DEVIATION OF J-TH COEFFICIENT         *
C      *                                                                    *
C      *    SUBROUTINES                                                     *
C      *                     BUCKLY                                         *
C
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(NTOT),X(MVAR,NTOT),Y(NTOT),ALPHA(MVAR),SIGMAA(MVAR)
       DIMENSION IWRK1(NTOT),IWRK2(NTOT),IWRK4(NTOT)
       DIMENSION IWRK5(NTOT),IWRK6(NTOT),IWRK7(NTOT),IWRK8(NTOT)
       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT)
       DIMENSION WRK5(NTOT),WRK6(NTOT),WRK7(NTOT),WRK8(NTOT)
       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
       DIMENSION SWRK5(MVAR),SWRK6(MVAR),SWRK7(MVAR),SWRK8(MVAR)
       DIMENSION SWRK9(MVAR),DWRK1(MVAR,NTOT),EWRK1(MVAR,MVAR)
       CHARACTER*9 OUTPUT
C
       NVAR1=NVAR+1
C
C      *      IF CENSORING IS DUE TO UPPER LIMITS, CHANGE THE SIGNS OF DATA *
C      *      X(I) AND Y(I) BECAUSE B-J METHOD ASSUMES LOWER LIMITS.        *
C
       IF(ICENS .LT. 0) THEN
          DO 1222 I=1,NTOT
             DO 1111 J=1,NVAR
                X(J,I)=-X(J,I)
 1111        CONTINUE
             Y(I)=-Y(I)
 1222     CONTINUE
       ENDIF
C
C
C      *      BUCKLY   : THE SUBROUTINE WHICH PERFORMS THE BUCKLEY AND      *
C      *                 JAMES METHOD.                                      *
C
C
       CALL BUCKLY(X,Y,ALPHA,IND,TOL,SIGMAA,NTOT,NVAR,ND,NC,MAX,ITE,
     +                   IWRK1,DWRK1,IWRK2,IWRK4,WRK1,WRK2,IWRK5,
     +                   WRK3,WRK4,WRK5,WRK6,WRK7,SWRK1,SWRK2,SWRK3,
     +                   IWRK6,IWRK7,IWRK8,SWRK4,SWRK5,SWRK6,SWRK7,
     +                   SWRK8,SWRK9,EWRK1,WRK8,MVAR)
C
C
C       *    CORRECT THE SIGNS OF THE DATA TO THE ORIGINAL ONES, IF THE     *
C       *    CENSORING IS UPPER LIMIT.                                      *
C
       IF(ICENS.LT.0) THEN
          DO 1223 I=1,NTOT

             DO 1113 J=1,NVAR
                X(J,I)=-X(J,I)
 1113        CONTINUE
             Y(I)=-Y(I)
 1223     CONTINUE
C
          ALPHA(1)=-ALPHA(1)
C
       ENDIF

  320  IF(OUTPUT.EQ.'         ') THEN
          PRINT 1050
          PRINT 1020
          PRINT 1050
C
          PRINT 1200,ALPHA(1)

          DO 452 J=2,NVAR1
             JI=J-1
             PRINT 1250,JI,ALPHA(J),SIGMAA(J)
  452     CONTINUE

          PRINT 1300,ALPHA(NVAR+2)
          PRINT 1350,ITE
          PRINT 1050
       ELSE
          WRITE(60,1050)
          WRITE(60,1020)
          WRITE(60,1050)
C
          WRITE(60,1200) ALPHA(1)

          DO 450 J=2,NVAR1
             JI=J-1
             WRITE(60,1250) JI,ALPHA(J),SIGMAA(J)
  450     CONTINUE

          WRITE(60,1300) ALPHA(NVAR+2)
          WRITE(60,1350) ITE
          WRITE(60,1050)
       ENDIF
C
C
 1020  FORMAT(T5,'LINEAR REGRESSION BY BUCKLEY-JAMES METHOD' )
 1050  FORMAT(T5,'      ')
 1100  FORMAT(T8,'DATA TITLE :',T25,60A1)
 1200  FORMAT(T8,'INTERCEPT COEFF    :',F8.4)
 1250  FORMAT(T8,'SLOPE COEFF ',I1,'      :',F8.4,T38,'+/-',T41,
     +         F8.4)
 1300  FORMAT(T8,'STANDARD DEVIATION :',F8.4)
 1350  FORMAT(T8,'ITERATIONS         :',I3)
 4000  RETURN
       END

C
C      **********************************************************************
C      *********************** SUBROUTINE BUCKLY  ***************************
C      **********************************************************************
C
       SUBROUTINE BUCKLY(X,Y,ALPHA,IND,TOL,SIGMAA,NTOT,
     +                   NVAR,NU,NC,MAX,ITE,IND2,XX,IPT,IR,ND,TY,
     +                   T,NO,Z,W,WX,ZY,V,TEST,TEST2,BU,
     +                   IWRK1,IWRK2,SWRK1,SWRK2,SWRK3,SWRK4,
     +                   SWRK5,SWRK6,EWRK1,WRK1,MVAR)
C
C
C      *     THIS IS A SUBPROGRAM WHICH PERFORMS THE BUCKLEY-JAMES          *
C      *     METHOD. THIS SUBROUTINE WAS ADAPTED FROM CODE BY J. HALPERN    *
C      *     (STANFORD UNIVERSITY SCHOOL OF MEDICINE, DEPARTMENT            *
C      *        OF FAMILY, COMMUNITY AND PREVENTIVE MEDICINE.)              *
C      *                                                                    *
C      *  INPUT                                                             *
C      *         X(J,I)   : INDEPENDENT VARIABLES                           *
C      *         Y(I)     : DEPENDENT VARIABLE                              *
C      *         IND(I)   : INDICATOR OF CENSORING                          *
C      *         TOL      : TOLERANCE LEVEL                                 *
C      *         NTOT     : NUMBER OF DATA POINTS                           *
C      *         NVAR     : NUMBER OF INDEPENDENT VARIABLES                 *
C      *         NU       : NUMBER OF DETECTED POINTS                       *
C      *         NC       : NUMBER OF CENSORED POINTS                       *
C      *         MAX      : MAXIMUM ITERATION                               *
C      *                                                                    *
C      *   WORK                                                             *
C      *          V(J)    : AVERAGE OF J-TH DETECTED INDEPENDENT VARIABLE   *
C      *          BU(J)   : VARIANCE OF J-TH DETECTED INDEPENDENT VARIABLE  *
C      *         TEST(J)  : STORE OF THE PREVIOUS STEP ESTIMATIONS OF       *
C      *                    ALPHA(J)                                        *
C      *          IR(I)   : SORTING ORDER                                   *
C      *          Z(I)    : RESIDUALS                                       *
C      *          W(I)    : KM ESTIMATOR                                    *
C      *          WX(I)   : WEIGHT                                          *
C      *                                                                    *
C      *  OUTPUT                                                            *
C      *       ALPHA(J)   : REGRESSION COEFFICIENTS                         *
C      *       SIGMA(J)   : ERROR                                           *
C      *       ITE        : ITERATION NUMBER                                *
C      *                                                                    *
C      *    SUBROUTINES                                                     *
C      *                    SORT1, REGRES                                   *
C
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)

       DIMENSION IND(NTOT),IND2(NTOT),XX(MVAR,NTOT),IPT(NTOT)
       DIMENSION X(MVAR,NTOT),Y(NTOT),IR(NTOT),ND(NTOT),TY(NTOT)
       DIMENSION T(NTOT),NO(NTOT),Z(NTOT),W(NTOT),WX(NTOT),ZY(NTOT)
       DIMENSION V(MVAR),ALPHA(MVAR),TEST(MVAR),BU(MVAR),SIGMAA(MVAR)
       DIMENSION TEST2(MVAR),IWRK1(NTOT),IWRK2(NTOT)
       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
       DIMENSION SWRK5(MVAR),SWRK6(MVAR),EWRK1(MVAR,MVAR),WRK1(NTOT)
C
C      *               INITIALIZATION                                       *
C
       ITE=0
       NB=NVAR+1
       DO 343 J=1,NVAR
          V(J) =0.0
          BU(J)=0.0
  343  CONTINUE

       DO 392 IN=1,NB
          TEST(IN)=0.0
       TEST2(IN) =0.0
  392  CONTINUE
C
C      *       CALCULATE SOME VALUES FOR THE STANDARD DEVIATION             *
C
       DO 5 I=1,NTOT
          IF(IND(I).EQ.0) THEN
             DO 63 J=1,NVAR
                V(J)=V(J)+X(J,I)
   63        CONTINUE
       ENDIF
    5  CONTINUE

       DO 68 J=1,NVAR
          V(J)=V(J)/REAL(NU)
   68  CONTINUE

       DO 51 I=1,NTOT
          IF(IND(I).EQ.0) THEN
             DO 53 J=1,NVAR
                BU(J)=BU(J)+(X(J,I)-V(J))**2
   53        CONTINUE
       ENDIF
   51  CONTINUE
C
C      *      REGRES : SUBPROGRAM FOR LINEAR REGRESSION WITHOUT             *
C      *               CONSIDERING CENSORING STATUS                         *
C
       CALL REGRES(X,Y,NTOT,NVAR,ALPHA,RMUL,SWRK1,SWRK2,WRK1,
     +             IWRK1,IWRK2,SWRK3,SWRK4,SWRK5,EWRK1,SWRK6,MVAR)
C
C      *              GET RESIDUALS Z(I)                                    *
C
C      *       START ITERATION : 2000 LOOP.                                 *
C
C
 2000  DO 31 I=1,NTOT
          T(I)=-400.0
          IND2(I)=IND(I)
C       
          ZS=0.0
          DO 61 J=1,NVAR
             JJ=J+1
             ZS=ZS+ALPHA(JJ)*X(J,I)
             XX(J,I)=X(J,I)
   61     CONTINUE
C
          Z(I)=Y(I)-ZS
   31  CONTINUE
C
C      *            SORTING .... INCREASING ORDER                           *
C
       CALL SORT1(IND2,XX,Z,NTOT,NVAR,IR,SWRK1,MVAR)
C
       DO 311 I=1,NTOT
          TY(I)=Y(IR(I))
          ZY(I)=Z(I)
  311  CONTINUE
C
C      *       THE LARGEST RESIDUAL MUST BE UNCENSORED.                     *
C
       IND2(NTOT)=0
C
C      *      ESTIMATE  VALUES FOR CENSORED DATA                            *
C      *                                                                    *
C      *    TY(I)=YY(I)*DEL+((ALPHA*X+SUM(WXX(K)*Z(K))/(1-W(I)))*(1-DEL)    *
C      *     WHERE                                                          *
C      *          TY   : ESTIMATED DEPENDENT VALUE                          *
C      *          DEL  : IF THE DATA IS UNCENSORED :DEL=1.0                 *
C      *                 IF THE DATA IS CENSORED  :DEL=0.0                  *
C      *          SUM  : SUM OVER UNCENSORED DATA Z(K)<Z(I)                 *
C      *          WX   : WEIGHT ... W(I-1)-W(I)                             *
C      *          W    : KAPLAN-MEIER PRODUCT LIMIT ESTIMATOR               *
C
C
       K=0
       DO 21 I=1,NTOT
          IF((IND2(I).NE.0).OR.(K.NE.0)) THEN
             IF(IND2(I).NE.0) THEN
                IPT(I)=K
                GOTO 21
             ENDIF
             IF((IND2(I).EQ.0).AND.(ZY(I).EQ.T(K))) THEN
                ND(K)=ND(K)+1
                IPT(I)=K
                GOTO 21
             ENDIF
          ENDIF
          ND(K+1)=1
          T(K+1)=ZY(I)
          IPT(I)=K+1
          K=K+1
   21  CONTINUE
C
       NI=K
       DO 28 I=1,NI
          NO(I)=0
          IDZ=0
   28  CONTINUE
C
       DO 29 I=1,NTOT
C
C      *      IF THE FIRST POINT IS A CENSORED VALUE, DROP THE POINT.       *
C      *      PT RUNS FROM 1 TO NU.                                         *
C
          IF(IPT(I).EQ.0) IDZ=IDZ+1
          IF(IPT(I).GT.0) NO(IPT(I))=NO(IPT(I))+1
   29  CONTINUE

       DENOM=REAL(NTOT-IDZ)
       W(1)=1.0-ND(1)/DENOM
       WJ=NTOT-NO(1)-IDZ
C
       DO 30 I=2,NI
          W(I)=W(I-1)*(1.0-ND(I)/WJ)
          WJ=WJ-NO(I)
   30  CONTINUE

       WX(1)=1.0-W(1)
C
       DO 41 I=2,NI
          WX(I)=W(I-1)-W(I)
   41  CONTINUE
C
       DO 83 I=1,NI
          WX(I)=WX(I)/REAL(ND(I))
   83  CONTINUE
C
       Z(NTOT)=0.0
       DO 36 JJ=2,NTOT
          I=NTOT-JJ+1
          Z(I)=Z(I+1)
          IF((ZY(I+1).GT.ZY(I)).AND.(IND2(I+1).EQ.0)) 
     +                           Z(I)=Z(I+1)+WX(IPT(I+1))*ZY(I+1)
   36  CONTINUE
C
       III=NTOT-1
       DO 32 I=1,III
          IF(IPT(I).GT.0) Z(I)=Z(I)/W(IPT(I))
          IF(IND2(I).NE.0) THEN
             ZS=0.0
             DO 62 J=1,NVAR
                JJ=J+1
                ZS=ZS+ALPHA(JJ)*XX(J,I)
   62        CONTINUE
             TY(I)=Z(I)+ZS
          ENDIF
   32  CONTINUE
C
       CALL REGRES(XX,TY,NTOT,NVAR,ALPHA,RMUL,SWRK1,SWRK2,WRK1,
     +             IWRK1,IWRK2,SWRK3,SWRK4,SWRK5,EWRK1,SWRK6,MVAR)
       ITE=ITE+1
C
C       *        TEST FOR CONVERGENCE                                       *
C       *     TEST(I) CONTAINS A PREVIOUS VALUE OF ALPHA(I).                *
C       *     IF NUMBER OF ITERATION EXCEEDS MAXITS, THE PROGRAM STOPS      *
C       *     TO CONTINUE ITERATION, EVEN IF TOLERANCE IS LARGER THAN       *
C       *     THE ASSIGNED VALUE. SINCE THE FINAL VALUES ARE OFTEN          *
C       *     TRAPPED IN OSCILLATION, THE PROGRAM TAKES AVERAGE OF THE      *
C       *     LAST TWO VALUES FOR THE FINAL OUTPUT.                         *
C
       IF(ITE.LE.MAX) THEN

          SUM=0.0
          DO 154 K=1,NB
             SUM=SUM+(ALPHA(K)-TEST2(K))**2
             TEST2(K)=TEST(K)
             TEST(K)=ALPHA(K)
  154     CONTINUE
C
C       *  IF THE DIFFERENCE IS LARGER THAN THE TOLERANCE, REPEAT ITERATION *
C
  234     IF(DSQRT(SUM).GT.TOL) GOTO 2000
       ENDIF
C
C       *  THE CONVERSION IS OBTAINED OR THE ITERATION REACHED THE MAX.     *
C
       DO 270 K=1,NB
          ALPHA(K)=(ALPHA(K)+TEST2(K))/2.0
  270  CONTINUE
C
C      *             CALCULATION OF VARIANCE ETC.                           *
C
       STD=0.0
       EM =0.0
       K=NTOT
       DO 35 I=1,NTOT
          IF(IND(I).EQ.0) THEN

             ZS=0.0
             DO 84 J=1,NVAR
                JJ=J+1
                ZS=ZS+ALPHA(JJ)*X(J,I)
   84        CONTINUE

             Z(I)=Y(I)-ZS
             EM=EM+Z(I)
          ENDIF
   35  CONTINUE

       EM=EM/REAL(NU)
       DO 37 I = 1,NTOT
          IF(IND(I) .EQ. 0) THEN
             EM2 = Z(I) - EM
             STD = STD+EM2*EM2
          ENDIF
   37  CONTINUE

       STD=STD/FLOAT(NU-NB)

       DO 76 I=1,NVAR
          SIGMAA(I+1)=DSQRT(STD/BU(I))
   76  CONTINUE

       ALPHA(NB+1)=DSQRT(STD)

       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE CENS  ******************************
C      **********************************************************************
C
        SUBROUTINE CENS(X,Y,IP,NTOT)
C
C      *         THIS SUBROUTINE ADDS OR SUBTRACTS 0.00001 TIMES THE        *
C      *         VALUE FROM CENSORED POINTS SO THAT IF THERE ARE TIES WITH  *
C      *         DETECTED POINTS, THE CENSORED POINTS CAN BE                *
C      *         DISTINGUISHED FROM THE DETECTIONS.                         *
C      *         IF THE SMALLEST DIGIT IS LESS THAN OR EQUAL TO             *
C      *         0.0001, THEN 'CONST' SHOULD BE CHANGED TO A                *
C      *         SMALLER VALUE.                                             *
C      *                                                                    *
C      *         INPUT AND OUTPUT:                                          *
C      *                   X(I)    : INDEPENDENT VARIABLE                   *
C      *                   Y(I)    : DEPENDENT VARIABLE                     *
C      *                  IP(I)    : CENSORED STATUS INDICATOR              *
C
        IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
        DIMENSION X(NTOT),Y(NTOT),IP(NTOT)
C
C
        CONST=1.000001
        DO 500 I=1,NTOT
C
C      *                 UPPER LIMITS CASES                                 *
C
        IF(IP(I).EQ.-1) THEN 
              Y(I)=Y(I)/CONST

          ELSEIF(IP(I).EQ.-2) THEN
              X(I)=X(I)/CONST

          ELSEIF(IP(I).EQ.-3) THEN
              X(I)=X(I)/CONST
              Y(I)=Y(I)/CONST

          ELSEIF(IP(I).EQ.-4) THEN
              X(I)=X(I)/CONST
              Y(I)=Y(I)*CONST

C
C      *                  LOWER LIMIT CASES                                 *
C
          ELSEIF(IP(I).EQ.1)  THEN
              Y(I)=Y(I)*CONST

          ELSEIF(IP(I).EQ.2)  THEN
              X(I)=X(I)*CONST

          ELSEIF(IP(I).EQ.3)  THEN
              X(I)=X(I)*CONST 
              Y(I)=Y(I)*CONST

          ELSEIF(IP(I).EQ.4)  THEN
              X(I)=X(I)*CONST
              Y(I)=Y(I)/CONST

        ENDIF
  500   CONTINUE
        RETURN
        END
C
C      **********************************************************************
C      ********************** SUBROUTINE BIVAR ******************************
C      **********************************************************************
C
       SUBROUTINE BIVAR(IBACK)
          
C
C      *   CORRELATION AND REGRESSION                                       *
C      *   PARAMETERS                                                       *
C      *      MVAR     :  THE MAXIMUM NUMBER OF VARIABLES (COLUMNS)         *
C      *                   ALLOWED IN A DATA SET.                           *
C      *      NDAT     :  THE MAXIMUM NUMBER OF DATA POINTS (ROWS)          *
C      *                   ALLOWED IN A DATA SET.                           *
C      *      IBIN     :  THE DIMENSION SIZE FOR BINS - USED IN SCHMITT'S   *
C      *                   PROCEDURE COMPUTATIONS.                          *
C      *      LENG     :  (MVAR+1)+NTOT - USED IN EM ALGORITHM COMPUTATIONS *
C      *      LEGWRK   :  (MVAR+1)*NTOT - USED IN EM ALGORITHM COMPUTATIONS *
C      *   INPUT                                                            *
C      *      FILE     :  NAME OF DATA FILE (9 LETTERS)                     *
C      *      TITLE    :  TITLE OF THE PROBLEM (80 LETTERS)                 *
C      *      NVAR     :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
C      *      NTOT     :  THE ACTUAL NUMBER OF DATA POINTS IN THE DATA SET  *
C      *      ICOL     :  INDICATOR OF VARIABLE (<=NVAR)                    *
C      *                  IF A MULTIVARIATE PROBLEM IS NEEDED, SET ICOL=0   *
C      *      COLM     :  NAME OF THE INDEPENDENT VARIABLE                  *
C      *      YNAME    :  NAME OF THE DEPENDENT VARIABLE                    *
C      *      COMMAND  :  NAME OF THE "COMMAND" FILE                        *
C      *      OUTPUT   :  NAME OF THE OUTPUT FILE                           *
C      *      IND(1,I) :  INDICATOR OF CENSORING                            *
C      *                    IF =0,  DETECTED                                *
C      *                       =1,  Y LOWER LIMIT                           *
C      *                       =2,  X LOWER LIMIT                           *
C      *                       =3,  DOUBLE LOWER LIMIT                      *
C      *                       =4,  X UPPER LIMIT AND Y LOWER LIMIT         *
C      *                  FOR THE UPPER LIMITS, CHANGE THE SIGN             *
C      *                  2, 3, AND 4 CAN BE USED ONLY IN GEN. KENDALL'S    *
C      *                  TAU, GEN. SPEARMAN'S RHO, AND SCHMITT'S METHOD    *
C      *      X(J,I)   :  INDEPENDENT VARIABLES                             *
C      *      Y(I)     :  DEPENDENT VARIABLE                                *
C      *     IPROG(I)  :  INDICATOR OF METHODS                              *
C      *     NOTEST    :  NUMBERS OF TEST                                   *
C      *  INPUT FOR EM ALGORITHM                                            *
C      *      TOL      :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX      :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      IBET     :  IF 0, NO DEPENDENT VARIABLE IS CONFINED BETWEEN   *
C      *                        TWO VALUES                                  *
C      *                     1, THERE ARE SOME DEPENDENT VARIABLE WHICH     *
C      *                        ARE CONFINED BETWEEN TWO VALUES             *
C      *    ALPHA(K)   :  INITIAL ESTIMATE OF REGRESSION COEFFICIENTS       *
C      *  INPUTS FOR BUCKLEY-JAMES METHOD                                   *
C      *      TOL1     :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX1     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *  INPUTS FOR SCHMITT'S BINNING METHOD                               *
C      *      MX       :  BIN NUMBER OF X AXES                              *
C      *      MY       :  BIN NUMBER OF Y AXES                              *
C      *      TOL3     :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                  *
C      *      MAX3     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      XBIN     :  BIN SIZE FOR X AXES                               *
C      *      YBIN     :  BIN SIZE FOR Y AXES                               *
C      *      XORG     :  ORIGIN OF X AXES                                  *
C      *      YORG     :  ORIGIN OF Y AXES                                  *
C      *      ISKIP    :  IF 0, THE PROGRAM WILL PROVIDE XBIN, YBIN, XORG,  *
C      *                        AND YORG.                                   *
C      *                    >0, THESE VALUES MUST BE PROVIDED BY THE USER   *
C      *      IPIRNT   :  IF 0, NO TWO DIMENSIONAL K-M ESTIMATOR WILL BE    *
C      *                        PRINTED                                     *
C      *                    >0, TWO DIMENSIONAL K-M ESTIMATOR WILL BE       *
C      *                        PRINTED                                     *
C      *                                                                    *
C      *    WORKING VARIABLES AND ARRAYS:                                   *
C      *      NTOT     :  NUMBER OF DATA POINTS                             *
C      *      ND       :  NUMBER OF DETECTED POINTS                         *
C      *      NC1      :  NUMBER OF Y LOWER LIMITS                          *
C      *      NC2      :  NUMBER OF X LOWER LIMITS                          * 
C      *      NC3      :  NUMBER OF DOUBLE LOWER LIMITS                     * 
C      *      NC4      :  NUMBER OF Y LOWER AND X UPPER LIMITS              *
C      *      NC5      :  NUMBER OF Y UPPER LIMITS                          *
C      *      NC6      :  NUMBER OF X UPPER LIMITS                          *
C      *      NC7      :  NUMBER OF DOUBLE UPPER LIMITS                     *
C      *      NC8      :  NUMBER OF Y UPPER AND X LOWER LIMITS              *
C      *      ICENS    :  IF 0, CENSORING IS MIXED                          *
C      *                     1, CENSORING IS Y LOWER LIMITS ONLY            *
C      *                    -1, CENSORING IS Y UPPER LIMITS ONLY            *
C      *      NYC      :  NC1+NC2                                           *
C      *      NXC      :  NC3+NC4                                           *
C      *      NBC      :  NC5+NC6+NC7+NC8                                   *
C      *      IDO      :  NXC+NBC                                           *
C      *      IMUL     :  INDICATOR OF MULTIVARIATE PROBLEM                 *
C      *      XX(J,I)  :  =X(ICOL,I), EXCEPT FOR MULTI INDEPENDENT VARIABLE *
C      *                  CASE (J=1,NVAR).                                  *
C      *      IND2(I)  :  =IND(1,I)                                         *
C      *                                                                    *
C      *  OUTPUT                                                            *
C      *     COXREG                                                         *
C      *      CHI      : GLOBAL CHI-SQUARE                                  *
C      *      PROB     : PROBABILITY FOR NULL                               *
C      *     BHK (GNERALIZED KENDALL'S TAU)                                 *
C      *       Z       : DEVIATION                                          *
C      *      PROB     : PROBABILITY FOR NULL                               *
C      *     EM ALGORITHM                                                   *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS  (K=1,NVAR+1)       *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     ITE       : NUMBER OF ITERATION                                *
C      *     BUCKLEY-JAMES                                                  *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS (K=1,NVAR+1)        *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     SCHMITT                                                        *
C      *     ALPHA     : INTERCEPT COEFFICIENT                              *
C      *     BETA      : SLOPE COEFFICIENT                                  *
C      *   *****  ALL OUTPUTS ARE INSIDE OF EACH SUBROUTINE                 *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *     DATA1, DATREG, DATA2, MULVAR                                   *
C      *                                                                    *

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

C      *   THIS PARAMETER STATEMENT AND THE ONE IN UNIVAR.F ARE THE ONLY    *
C      *   STATEMENTS THAT NEED TO BE ADJUSTED IF THE USER WISHES TO        *
C      *   ANALYZE DATA SETS OF MORE THAN 500 OBSERVATIONS OR MORE THAN     *
C      *   VARIABLES.                                                       *

C  **************************************************************************
       PARAMETER(MVAR=4, NDAT=500, IBIN=50)
C  **************************************************************************

       CHARACTER*1 CHECK,CHAR(4,10)
       CHARACTER*7 BB(10),YY
       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM
       CHARACTER*9 YNAME,DUMMY1
       CHARACTER*80 TITLE 
       
       DIMENSION IND(MVAR,NDAT),X(MVAR,NDAT),Y(NDAT)
       DIMENSION IPROG(6),IIND(NDAT)

       DIMENSION DWRK1(MVAR,NDAT),DWRK2(MVAR,NDAT)
       DIMENSION DWRK3(MVAR,NDAT),DWRK4(MVAR,NDAT)
       DIMENSION DWRK5(MVAR,NDAT),DWRK6(MVAR,NDAT)
       DIMENSION DWRK8(MVAR,NDAT)
       
       DIMENSION EWRK1(4,4),RWRK1(NDAT,MVAR)

       DIMENSION AWRK(5,IBIN)
       DIMENSION WWRK1((MVAR+1)+NDAT)
       DIMENSION WWRK2((MVAR+1)+NDAT)
       DIMENSION VWRK1((MVAR+1)*NDAT)
       DIMENSION VWRK2((MVAR+1)*NDAT)

       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
       DIMENSION WRK12(NDAT)

       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR)
       DIMENSION SWRK4(MVAR),SWRK5(MVAR),SWRK6(MVAR)
       DIMENSION SWRK7(MVAR),SWRK8(MVAR),SWRK9(MVAR)
       DIMENSION SWRK10(MVAR),SWRK11(MVAR),SWRK17(MVAR)

       DIMENSION LWRK1(MVAR,NDAT), LWRK2(MVAR,NDAT)
       DIMENSION LWRK3(MVAR,NDAT)
       DIMENSION IWRK1(NDAT),IWRK2(NDAT),IWRK3(NDAT)
       DIMENSION IWRK4(NDAT),IWRK5(NDAT),IWRK6(NDAT)
       DIMENSION IWRK7(NDAT),IWRK8(NDAT)

       DIMENSION IWRK9(NDAT),CWRK1(IBIN),CWRK2(IBIN)

       DIMENSION IBWRK1(IBIN,IBIN),IBWRK2(IBIN,IBIN)
       DIMENSION IBWRK3(IBIN,IBIN),IBWRK4(IBIN,IBIN)
       DIMENSION IBWRK5(IBIN,IBIN),IBWRK6(IBIN,IBIN)
       DIMENSION IBWRK7(IBIN,IBIN),IBWRK8(IBIN,IBIN)
       DIMENSION IBWRK9(IBIN,IBIN)
       DIMENSION BWRK1(IBIN,IBIN),BWRK2(IBIN,IBIN)

       LENG = (MVAR+1)+NDAT
       LEGWRK = (MVAR+1)*NDAT

       DO 5000 K=1,10
       BB(K)='     X '
 5000  CONTINUE
       YY='     Y '
C
 6000  PRINT *
       PRINT *
       PRINT *,'      CORRELATION AND REGRESSION CALCULATIONS'
       PRINT *
       PRINT *,' CORRELATION OPTIONS        LINEAR REGRESSION OPTIONS'
       PRINT *,' 1. COX HAZARD MODEL        4. EM ALGORITHM WITH '
       PRINT *,'                               NORMAL DISTRIBUTION'
       PRINT *,' 2. GEN. KENDALL`S TAU      5. BUCKLEY-JAMES METHOD'
       PRINT *,' 3. GEN. SPEARMAN`S RHO     6. SCHMITT`S BINNING METHOD'
       PRINT *
       PRINT *,' DATA SETS WITH CENSORING IN ONLY ONE DIRECTION OF THE'
       PRINT *,' DEPENDENT VARIABLE CAN USE ALL METHODS.'
       PRINT *
       PRINT *,' DATA SETS WITH SEVERAL INDEPENDENT AND ONE DEPENDENT'
       PRINT *,' VARIABLE CAN USE ONLY THE COX PROPORTIONAL HAZARD'
       PRINT *,' MODEL,EM ALGORITHM, OR BUCKLEY-JAMES METHOD.  ONLY'
       PRINT *,' ONE TYPE OF CENSORING IN THE DEPENDENT VARIABLE IS'
       PRINT *,' ALLOWED.'
       PRINT *
       PRINT *,' IF YOUR DATA SET HAS CENSORED POINTS IN THE '
       PRINT *,' INDEPENDENT VARIABLE AND/OR DUAL CENSORED POINTS,'
       PRINT *,' YOU CAN USE ONLY THE GEN. KENDALL`S TAU OR GEN.'
       PRINT *,' SPEARMAN`S RHO CORRELATION COEFFICIENT, OR'
       PRINT *,' SCHMITT`S BINNED LINEAR REGRESSION.'
       PRINT *
 6010  PRINT *
C
C      *  CHECK WHETHER THE USER WANTS TO USE COMMAND FILE INPUTS. IF SO,   *
C      *  GO TO 6660                                                        *
C
   50  FORMAT(A1)
 1380  FORMAT(A9)
C
       OUTPUT='         '
       ICOMM=0
       ICOL=1
       PRINT *,'DO YOU WANT TO READ ALL INFORMATION'
       WRITE(6,6020)
 6020  FORMAT('      FROM A COMMAND FILE (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6660
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6025
       GOTO 6010
C
C      *          READ FROM THE TERMINAL                                    *
C
 6025  PRINT *
       PRINT *,'          OK, LET US READ FROM THE TERMINAL '
C
C      *           READ TITLE                                               *
C
 6030  PRINT *
       WRITE(6,6040)
 6040  FORMAT('WHAT IS THE TITLE OF THE PROBLEM ? ')
       READ(5,6050) TITLE
 6050  FORMAT(A80)
C
C      *           READ DATA FILE NAME                                      *
C
 6051  PRINT *
       WRITE(6,6052)
 6052  FORMAT('WHAT IS THE DATA FILE NAME ? ')
       READ(5,1380) FILE
C
C      *           READ NUMBER OF INDEPENDENT VARIABLES                     *
C
 6060  PRINT *
       WRITE(6,6070)
 6070  FORMAT('HOW MANY INDEPENDENT VARIABLES DO YOU HAVE ? ')
       CALL DATA1(NVAR)
       IF((NVAR.GE.1).AND.(NVAR.LE.MVAR-2)) GOTO 6080
       PRINT *
       PRINT *,'    YOUR CHOICE IS NOT ACCEPTABLE. PLEASE TYPE IT AGAIN'
       GOTO 6060
C
C      *    CALL SUBROUTINE "DATREG" TO READ DATA                           *
C
 6080  CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
       DO 6090 I = 1, NTOT
            IIND(I) = IND(1,I)    
 6090  CONTINUE
C
C      *        CHECK WHICH METHODS THE USER CAN USE                        *
C
       IDC=NXC+NBC
       IF((NVAR.EQ.1).AND.(IDC.EQ.0)) GOTO 6530
       IF((NVAR.NE.1).AND.(IDC.NE.0)) GOTO 6340
       IF((NVAR.EQ.1).AND.(IDC.NE.0)) GOTO 6400
 6170  PRINT *
       WRITE(6,6180)
 6180  FORMAT('IS THIS A MULTIVARIATE PROBLEM  (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6220
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6340
       GOTO 6170
C
C      *      DATA SET WITH MORE THAN ONE INDEPENDENT VARIABLES             *
C
 6220  PRINT *
       PRINT *,'           YOU CAN USE THE NEXT METHODS   '
       PRINT *
       PRINT *,'         1. COX HAZARD METHOD'
       PRINT *,'         4. EM ALGORITHM WITH NORMAL DISTRIBUTION'
       PRINT *,'         5. BUCKLEY-JAMES METHOD'
C
       ICOL=0
       J=1
 6230  PRINT *
       WRITE(6,6240)
 6240  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
       CALL DATA1(IPROG(J))
       IF((IPROG(J).EQ.1).OR.(IPROG(J).EQ.4).OR.(IPROG(J).EQ.5))
     +                                                    GOTO 6245
       GOTO 6230
 6245  IF(J.EQ.1) GOTO 6260
       J1=J-1
       DO 6250 K=1,J1
       IF(IPROG(K).NE.IPROG(J)) GOTO 6250
       PRINT *
       PRINT *,' YOU ALREADY CHOSE THAT METHOD.'
       PRINT *,'  PLEASE CHOOSE ANOTHER ONE'
       GOTO 6230
 6250  CONTINUE
 6260  IF(J.GE.3) GOTO 6280
 6265  PRINT *
       WRITE(6,6270)
 6270  FORMAT('DO YOU WANT TO USE ANY OTHER METHOD (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6230
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6280
       GOTO 6265
C
 6280  NOTEST=J
       GOTO 6604
C
C      *    FOR THE CASE  THAT THE DATA SET CONTAINS MIXED CENSORING        *
C      *    (THAT IS, UPPER AND LOWER LIMITS SIMULTANEOUSLY AND/OR          *
C      *    CENSORING IN BOTH VARIABLES).                                   *
C
 6340  PRINT *
       WRITE(6,6350)
 6350  FORMAT('WHICH INDEPENDENT VARIABLE DO YOU WANT TO USE ? ')
       CALL DATA1(ICOL)
       IF(ICOL.GT.NVAR) GOTO 6340
       IF(ICOL.LE.0) GOTO 6340
C
 6400  IF(NBC.EQ.0) GOTO 6530
       J=1
       PRINT *
       PRINT *,'          YOU CAN USE THE FOLLOWING METHODS'
       PRINT *,'          2. GEN. KENDALL`S TAU METHOD'
       PRINT *,'          3. GEN. SPEARMAN`S RHO METHOD'
       PRINT *,'          6. SCHMITT`S BINNING METHOD'
 6410  PRINT *
       WRITE(6,6420)
 6420  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
       CALL DATA1(IPROG(J))
       IF((IPROG(J).EQ.2).OR.(IPROG(J).EQ.3).OR.(IPROG(J).EQ.6)) 
     +                                                     GOTO 6425
       GOTO 6410
 6425  IF(J.EQ.1) GOTO 6440
       J1=J-1
       DO 6430 K=1,J1
       IF(IPROG(K).NE.IPROG(J)) GOTO 6430
       PRINT *
       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
       PRINT *,'    PLEASE CHOOSE THE OTHER ONE'
       GOTO 6410
 6430  CONTINUE
 6440  IF(J.EQ.3) GOTO 6600
 6450  PRINT *
       WRITE(6,6460)
 6460  FORMAT('DO YOU WANT TO USE THE OTHER METHOD, TOO (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6410
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
       GOTO 6450
C
C      *  FOR THE CASE THAT THE DATA SET CONTAINS ONE INDEPENDENT AND ONE   *
C      *  DEPENDENT VARIABLES AND ONE KIND OF CENSORING IN THE DEPENDENT    *
C      *  VARIABLE.                                                         *
C
 6530  PRINT *
       PRINT *,'     YOU CAN USE THE FOLLOWING METHODS'
       PRINT *
       PRINT *,'     1. COX HAZARD MODEL    4. EM ALGORITHM WITH'
       PRINT *,'                               NORMAL DISTRIBUTION'
       PRINT *,'     2. KENDALL`S TAU       5. BUCKLEY-JAMES REGRESSION'
       PRINT *,'     3. SPEARMAN`S RHO',
     +         '      6. SCHMITT`S BINNED REGRESSION'
       PRINT *
       J=1
 6540  PRINT *
       WRITE(6,6550)
 6550  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
       CALL DATA1(IPROG(J))
       IF(IPROG(J).LT.1) GOTO 6540
       IF(IPROG(J).GT.6) GOTO 6570
       IF(J.EQ.1) GOTO 6580
       J1=J-1
       DO 6560 K=1,J1
       IF(IPROG(K).NE.IPROG(J)) GOTO 6560
       PRINT *
       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
       PRINT *,'    PLEASE CHOOSE ANOTHER ONE.'
       GOTO 6540
 6560  CONTINUE
 6570  IF(J.GE.6) GOTO 6600
 6580  PRINT *
       WRITE(6,6590)
 6590  FORMAT('DO YOU WANT TO USE ANOTHER METHOD (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6540
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
       GOTO 6580
C
 6600  NOTEST=J
C
C      *        NAME THE VARIABLES                                          *
C
 6601  PRINT *
       PRINT *,'          PLEASE NAME THE VARIABLES : '
       PRINT *
       WRITE(6,6602)
 6602  FORMAT('WHAT IS THE NAME OF THE INDEPENDENT VARIABLE ? ')
       READ(5,1380) COLM
C
       PRINT *
       WRITE(6,6603)
 6603  FORMAT('WHAT IS THE NAME OF THE DEPENDENT VARIABLE ? ')
       READ(5,1380) YNAME
C
 6604  PRINT *
       WRITE(6,6605)
 6605  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA  (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 6609
       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 6609
       GOTO 6604
C
C      *   CHECK WHETHER THE USER WANT TO SAVE THE RESULT IN AN OUTPUT FILE *
C
 6609  PRINT *
       PRINT *,'DO YOU WANT TO SAVE THE RESULT '
       WRITE(6,6610)
 6610  FORMAT('     IN AN OUTPUT FILE (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6620
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 7116
       GOTO 6609
 6620  PRINT *
       WRITE(6,6630)
 6630  FORMAT('WHAT IS THE NAME OF THE OUTPUT FILE ? ')
       READ(5,1380) OUTPUT
       GOTO 7116
C
C
C
C      *    USE "COMMAND" FILE FOR INPUTS                                   *
C
 6660  PRINT *
       WRITE(6,6670)
 6670  FORMAT('WHAT IS THE NAME OF THE COMMAND FILE ? ')
       READ(5,1380) COMMAND
C
 6700  OPEN(UNIT=50, FILE=COMMAND, STATUS='OLD', FORM='FORMATTED')
       ICOMM=1
C
C      *   READ TITLE OF THE PROBLEM ; NAME OF THE DATA FILE                *
C
       READ(50,6710) TITLE
 6710  FORMAT(A80)
       READ(50,1380) FILE
C
C      * READ NUMBER OF VARIABLES; WHICH VARIABLE WILL BE USED; AND HOW     *
C      * MANY METHODS THE USER WANTS TO USE.                                *
C
       READ(50,6720) ((CHAR(I,J),I=1,4),J=1,3)
 6720  FORMAT(20A1)
       CALL DATA2(CHAR,1,3,NVAR,LIND)
       IF(LIND.EQ.0) GOTO 6750
 6730  PRINT *
       PRINT *,'   TOTAL NUMBER OF INDEPENDENT VARIABLES IS NOT CLEAR.'
       STOP
 6750  IF(NVAR.LT.1) GOTO 6730
       IF(NVAR.NE.1) GOTO 6760
       ICOL=1
       GOTO 6915
C 
 6760  CALL DATA2(CHAR,50,3,ICOL,LIND)
       IF(LIND.EQ.0) GOTO 6900
 6860  PRINT *
       PRINT *,'     THE CHOICE OF THE VARIABLE IS NOT CLEAR'
       STOP
 6900  IF(ICOL.LE.0) GOTO 6860
       IF(ICOL.GT.NVAR) GOTO 6860
C
C      *         CHOICE OF THE METHODS                                      *
C
 6915  CALL DATA2(CHAR,3,3,NOTEST,LIND)
       IF(LIND.EQ.0) GOTO 6950
 6930  PRINT *
       PRINT *,'    IT IS NOT CLEAR HOW MANY METHODS YOU WANT TO USE '
       STOP
 6950  IF(NOTEST.LE.0) GOTO 6930
       IF(NOTEST.GT.6) GOTO 6930
C
       READ(50,6960) ((CHAR(I,J),I=1,4),J=1,NOTEST)
 6960  FORMAT(30A1)
       DO 7020 I=1,NOTEST
       CALL DATA2(CHAR,I,NOTEST,IPROG(I),LIND)
       IF(LIND.EQ.0) GOTO 7010
 6970  PRINT *
       IF(I.EQ.1) PRINT *,'     FIRST PROGRAM NUMBER IS NOT CLEAR'
       IF(I.EQ.2) PRINT *,'     SECOND PROGRAM NUMBER IS NOT CLEAR'
       IF(I.GE.3) WRITE(6,6780) I
 6780  FORMAT(5X,I4,'-TH PROGRAM NUMBER IS NOT CLEAR')
       STOP
 7010  IF(IPROG(I).LE.0) GOTO 6970
       IF(IPROG(I).GT.6) GOTO 6970
 7020  CONTINUE
C
C      *  READ NAMES OF THE INDEPENDENT AND DEPENDENT VARIABLES. IF THE     *
C      *  PROBLEM HAS MULTI-INDEPENDENT VARIABLES, THESE NAMES WIL BE       *
C      *  IGNORED.                                                          *
C
       READ(50,7022) COLM,YNAME
 7022  FORMAT(2A9)
C
       CLOSE(UNIT=50,STATUS='KEEP')
C
C      *      CALL SUBROUTINE "DATREG" TO READ DATA                         *
C
       CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
C
       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
C
C      * THE NEXT SEVERAL LINES READ IN DUMMY VALUES TO PREVENT READING     *
C      * THE COMMANDS A SECOND TIME.                                        *
C
       READ(50,1380) DUMMY1
       READ(50,1380) DUMMY1
       READ(50,7029) IDUMMY
       READ(50,7029) IDUMMY
       READ(50,1380) DUMMY1
 7029  FORMAT(I4)
C
C      *  CHECK WHETHER THE ASSIGNED METHODS CAN BE USED FOR THE DATA       *
C
       IF(NVAR.GE.2) GOTO 7070
       IF((NXC.EQ.0).AND.(NBC.EQ.0)) GOTO 7110
C
C      *   THE CASE WITH MIXED CENSORING IN DATA                            *
C
       I=1
 7030  IF(IPROG(I).NE.1) GOTO 7040
       PRINT *
       PRINT *,'      YOU CANNOT USE COX HAZARD MODEL FOR THIS DATA SET'
       PRINT *,'      THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7040  IF(IPROG(I).NE.4) GOTO 7050
       PRINT *
       PRINT *,'       YOU CANNOT USE EM ALGORITHM FOR THIS DATA SET'
       PRINT *,'       THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7050  IF(IPROG(I).NE.5) GOTO 7060
       PRINT *
       PRINT *,'       YOU CANNOT USE BUCKLEY-JAMES METHOD FOR THIS'
       PRINT *,'        DATA SET. THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7060  IF(I.GE.NOTEST) GOTO 7110
       I=I+1
       GOTO 7030
C
C      *     THE CASE WITH MORE THAN ONE INDEPENDENT VARIABLES              *
C
 7070  I=1
 7080  IF(IPROG(I).NE.2) GOTO 7085
       PRINT *
       PRINT *,'         YOU CANNOT USE THE KENDALL`S TAU METHOD FOR'
       PRINT *,'         THIS DATA SET'
       PRINT *,'         THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7085  IF(IPROG(I).NE.3) GOTO 7090
       PRINT *
       PRINT *,'       YOU CANNOT USE SPEARMAN`S RHO FOR THIS DATA SET'
       PRINT *,'       THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7090  IF(IPROG(I).NE.6) GOTO 7100
       PRINT *
       PRINT *,'         YOU CANNOT USE SCHMITT`S BINNED REGRESSION'
       PRINT *,'         THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7100  IF(I.EQ.NOTEST) GOTO 7110
       I=I+1
       GOTO 7080
C
C      *        READ PRINT OUT INDICATOR FOR THE DATA                       *
C
 7110  READ(50,6960) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 7114
 7112  PRINT *
       PRINT *,'     THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       STOP
 7114  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 7115
       GOTO 7112
C
C      *       READ OUTPUT FILE NAME                                        *
C
 7115  READ(50,1380) OUTPUT
C
C      * CALL SUBROUTINE "MULVAR" TO COMPUTE CORRELATION/REGRESSION PROBLEMS*
C

 7116  IF(OUTPUT .NE. '         ') OPEN(UNIT=60,FILE=OUTPUT,
     +                                 STATUS='NEW'
C  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
C  VAX/VMS MACHINES.
C     +                                 ,CARRIAGECONTROL='LIST'
     +                                 )


       CALL  MULVAR(X,Y,IND,NTOT,ICOL,NVAR,NOTEST,IPROG,ICOMM,
     +              OUTPUT,COLM,FILE,YNAME,TITLE,ND,NYC,ICENS,
     +              NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,MVAR,
     +              LENG,LEGWRK,IBIN,DWRK1,IWRK9,SWRK17,DWRK2,
     +              DWRK3,DWRK4,DWRK5,DWRK6,DWRK8,RWRK1,
     +              EWRK1,AWRK,WWRK1,WWRK2,
     +              VWRK1,VWRK2,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
     +              WRK7,WRK8,WRK9,WRK10,WRK11,WRK12,
     +              SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
     +              SWRK8,SWRK9,SWRK10,SWRK11,LWRK1,LWRK2,LWRK3,
     +              IWRK1,IWRK2,IWRK3,IWRK4,IWRK5,IWRK6,IWRK7,
     +              IWRK8,CWRK1,CWRK2,IBWRK1,IBWRK2,IBWRK3,
     +              IBWRK4,IBWRK5,IBWRK6,IBWRK7,IBWRK8,IBWRK9,
     +              BWRK1,BWRK2)
C
       IF(IDATA.EQ.0) GOTO 7219
       IF(OUTPUT.NE.'         ') WRITE(60,7140)
       IF(OUTPUT.NE.'         ') WRITE(60,7117) FILE
       IF(OUTPUT.EQ.'         ') PRINT 7140
       IF(OUTPUT.EQ.'         ') PRINT 7117, FILE
 7117  FORMAT(5X,'INPUT DATA FILE : ',A9)
       IF(ICOL.NE.0) GOTO 7130
       IF(OUTPUT.NE.'         ') WRITE(60,7118) (BB(K),K,K=1,NVAR),YY
       IF(OUTPUT.EQ.'         ') PRINT 7118,(BB(K),K,K=1,NVAR),YY
 7118  FORMAT(4X,'CENSORSHIP',12(A7,I2,1X))
       DO 7119 I=1,NTOT
       IF(OUTPUT.NE.'         ') WRITE(60,7120) IIND(I),
     +                                      (X(J,I),J=1,NVAR),Y(I)
       IF(OUTPUT.EQ.'         ') PRINT 7120,IIND(I),
     +                                      (X(J,I),J=1,NVAR),Y(I)
 7119  CONTINUE
 7120  FORMAT(7X,I4,3X,10F10.3)
       GOTO 7219
 7130  IF(OUTPUT.NE.'         ') WRITE(60,7133)
       IF(OUTPUT.EQ.'         ') PRINT 7133
 7133  FORMAT(5X,' CENSORSHIP        X        Y')
       DO 7134 I=1,NTOT
       IF(OUTPUT.NE.'         ') WRITE(60,7135) IIND(I),X(ICOL,I),Y(I)
       IF(OUTPUT.EQ.'         ') PRINT 7135,IIND(I),X(ICOL,I),Y(I)
 7134  CONTINUE
 7135  FORMAT(8X,I4,5X,2F10.3)
 7140  FORMAT('     ')
C
 7219  IF(OUTPUT .NE. '         ') CLOSE(UNIT=60)
       PRINT *
       PRINT *
       PRINT *,'    COMPUTATIONS FOR CORRELATION/REGRESSION'
       PRINT *,'                          PROBLEMS ARE FINISHED'
 7220  PRINT *
       WRITE(6,7230)
 7230  FORMAT('DO YOU WANT TO DO ANY OTHER ANALYSIS (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1 
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
       GOTO 7220
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE CHOL  *******************************
C      **********************************************************************
C
       SUBROUTINE CHOL(A,N,U,NULLTY,NA,NU,IFAULT)
C
C      *      ALGORITHM AS 6 J.R.STATIST.SOC.C.(1968) VOL.17, NO.2          *
C      *                                                                    *
C      *      GIVEN A SYMMETRIC MATRIX OF ORDER N AS A LOWER TRIANGLE       *
C      *      IN A( ),  CALCULATE AN UPPER TRIANGLE, U( ), SUCH THAT        *
C      *      UPRIME*U=A. U( ) MAY COINCIDE WITH A( ). A( ) MUST BE         *
C      *      POSITIVE SEMIDEFINITE.                                        *
C      *      ETA IS SET TO MULTIPLYING FACTOR DETERMINING THE              *
C      *      EFFECTIVE  ZERO FOR PIVOT.                                    *
C      *      NULLTY IS RETURNED AS NO. OF EFFECTIVE ZERO PIVOTS.           *
C      *      IFAULT IS RETURNED AS 1,IF N.LE.0, 2,IF A( ) IS NOT           *
C      *      POSITIVE SEMI-DEFINITE WITHIN THE TOLERANCE BY ETA.           *
C      *      OTHERWISE ZERO.                                               *
C
C      *        NOTE : VARIABLES NA,NU, HAVE BEEN ADDED TO THE              *
C      *               ARGUMENT LIST AND USED TO DIMENSION TO ARRAYS        *
C      *               A AND U, RESPECTIVELY. (BY WOLYNETZ (1979))          *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION A(NA),U(NU)
C
       DATA ETA /1.0E-9/
C
C      *       THE VALUE OF ETA WILL DEPEND ON THE WORD LENGTH OF           *
C      *       THE COMPUTER BEING USED.                                     *
C
       IFAULT=1
       IF(N.GT.0) THEN
          IFAULT=2
          NULLTY=0
          J=1
          K=0

          DO 10 ICOL=1,N
             L=0

             DO 11 IROW=1,ICOL
                K=K+1
                W=A(K)
                M=J

                DO 12 I=1,IROW
                   L=L+1
                   IF(I.EQ.IROW) GOTO 13
                   W=W-U(L)*U(M)
                   M=M+1
   12           CONTINUE

   13           IF(IROW.EQ.ICOL) GOTO 14
                IF(U(L).EQ.0.0) THEN
                   U(K) = 0.0
                ELSE
                   U(K)=W/U(L)
                ENDIF
   11        CONTINUE

   14        IF(DABS(W).GE.DABS(ETA*A(K))) THEN
                IF(W.LT.0.0) GOTO 100
                U(K)=DSQRT(W)
             ELSE
                U(K)=0.0
                NULLTY=NULLTY+1
             ENDIF
             J=J+ICOL
   10     CONTINUE

          IFAULT=0

       ENDIF
  100  RETURN
       END

C
C      **********************************************************************
C      ************************ SUBROUTINE COEFF  ***************************
C      **********************************************************************
C
       SUBROUTINE COEFF(I,X,IP,NTOT,ICOEFF,IA,IB,IC,ID,IE,IG,IH,IJ)
C
C      *        SUBROUTINE WHICH FINDS CONCORDANCE INFORMATION  OF          *
C      *        THE QUANTITY X(I).                                          *
C      *                                                                    *
C      *      INPUT  :   X(I)      : THE QUANTITIY TO BE EXAMINED           *
C      *                IP(I)      : CENSORED STATUS OF X(I)                *
C      *                NTOT       : NUMBER OF DATA                         *
C      *      OUTPUT : ICOEFF(I)   : CONCORDANCE INFORMATION:               *
C      *                             FOR X(I) AND X(J) WITH I<J,            *
C      *                             IF X(I)<X(J), ICOEFF= 1                *
C      *                             IF X(I)>X(J), ICOEFF=-1                *
C      *                             OTHERWISE,    ICOEFF= 0                *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION X(NTOT),IP(NTOT),ICOEFF(NTOT)
C
       DO 100 J=1,NTOT
          ICOEFF(J)=0
          IF(X(I).LT.X(J)) THEN
             IF(IP(I).EQ.IA) GOTO 100
             IF(IP(J).EQ.ID) GOTO 100
             IF(IP(I).EQ.IB) GOTO 100
             IF(IP(J).EQ.IE) GOTO 100
             IF(IP(I).EQ.IC) GOTO 100
             IF(IP(J).EQ.IG) GOTO 100

             ICOEFF(J)=1

          ELSEIF(X(I).GT.X(J)) THEN
   50        IF(IP(I).EQ.ID) GOTO 100
             IF(IP(J).EQ.IA) GOTO 100
             IF(IP(I).EQ.IE) GOTO 100
             IF(IP(J).EQ.IB) GOTO 100
             IF(IP(I).EQ.IG) GOTO 100
             IF(IP(J).EQ.IC) GOTO 100

             ICOEFF(J)=-1
          ENDIF
  100  CONTINUE
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE COXREG  *****************************
C      **********************************************************************
C
       SUBROUTINE COXREG(IND,XX,YY,NTOT,NVAR,OUTPUT,ICENS,
     +                      RINFO,SCORE,FINFO,IL,IM,IP,Y,X,
     +                      SWRK1,IWRK1,IWRK2,MVAR)
C
C      *      THIS PROGRAM COMPUTES A CORRELATION PROBABILITY ACCORDING     *
C      *      TO COX'S (1972) PROPORTIONAL HAZARDS MODEL.                   *
C      *      ONLY ONE TYPE OF CENSORING (I.E. LOWER OR UPPER)              *
C      *      IS ALLOWED IN Y, BUT UP TO NVAR INDEPENDENT VARIABLES CAN     *
C      *      BE USED. THE HYPOTHESIS TESTED IS THE ABSENCE OF CORRELATION  *
C      *      BETWEEN THE DEPENDENT VARIABLE AND INDEPENDENT VARIABLES.     *
C      *      THEREFORE, THE  REGRESSION COEFFICIENT IN COX MODEL BETA      *
C      *      IS SET TO ZERO.                                               *
C
C      * NOTE NOTE NOTE:   THE PROBABILITY CALCULATED MAY NOT BE ACCURATE   *
C      *      WHEN THERE ARE A LARGE NUMBER OF TIED OBSERVATIONS (CF.       *
C      *      R. G. MILLER, SURVIVAL ANALYSIS, 1981, PP. 136-7).            *
C
C      *                                                                    *
C      *      INPUT    IND(I)  :  INDICATOR OF CENSORING                    *
C      *                           0 : UNCENSORED DATA POINT                *
C      *                           1 : Y(I) IS LOWER LIMIT                  *
C      *                          -1 : Y(I) IS UPPER LIMIT                  *
C      *              XX(J,I)  :  INDEPENDENT VARIABLES (J=1,..NVAR)        *
C      *              YY(I)    :  DEPENDENT VARIABLE                        *
C      *               NTOT    :  TOTAL NO. OF OBSERVATIONS                 *
C      *               NVAR    :  NO. OF INDEPENDENT VARIABLES              *
C      *                                                                    *
C      *       WORK      DF    :  DEGREE OF FREEDOM                         *
C      *                X(J,I) :  =XX(J,I)                                  *
C      *                Y(I)   :  =YY(I)                                    *
C      *                IP(I)  :  =IND(I)                                   *
C      *                IL(I)  :  INDICATOR OF TIES (#  OF TIES)            *
C      *                IM(I)  :  INDICATOR OF TIES (POSITION)              *
C      *               RINFO(I):  INFORMATION MATRIX AND ITS INVERSE        *
C      *                          MATRIX AFTER CALLING SUBROUTINE           *
C      *                          MATINV.                                   *
C      *               SCORE(I):  SCORE VECTOR                              *
C      *                                                                    *
C      *     OUTPUT       CHI  :  GLOBAL CHI-SQUARE                         *
C      *                 PROB  :  PROBABILITY OF CORRELATION                *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *                          SORT1, TIE, MATINV, PCHISQ                *
C      *                                                                    *
C      *     REFERENCE:    RUPERT G. MILLER JR., "SURVIVAL ANALYSIS", 1981, *
C      *                         JOHN WILEY & SONS (NY:NY)                  *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION RINFO(MVAR,MVAR),SCORE(MVAR),FINFO(MVAR)
       DIMENSION IND(NTOT),IL(NTOT),IM(NTOT),IP(NTOT),Y(NTOT),YY(NTOT)
       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT)
       DIMENSION SWRK1(MVAR), IWRK1(MVAR),IWRK2(MVAR)
       CHARACTER*9 OUTPUT
C
       DF=NVAR
C
C      *   SUBSTITUTE XX,YY,AND IND TO X, Y, IP TO AVOID ALTERATION OF      *
C      *   THE ORIGINAL DATA                                                *
C
       DO 20 I=1,NTOT
          IP(I)=IND(I)
          Y(I)=YY(I)
          IF(ICENS.EQ.-1) Y(I)=-YY(I)
C
C      *  IF THE OBSERVATION IS CENSORED, ADD A SMALL NUMBER TO AVOID TIES  *
C      *  WITH DETECTED VALUE.                                              *
C
          IF(IP(I).NE.0) Y(I)=Y(I)*(1.0+FLOAT(ICENS)*0.0000001)
C
          DO 10 J=1,NVAR
             X(J,I)=XX(J,I)
             IF(ICENS.EQ.-1) X(J,I)=-X(J,I)
   10     CONTINUE
   20  CONTINUE
C
C      *           SORT Y IN ASCENDING ORDER                                *
C
       CALL SORT1(IP,X,Y,NTOT,NVAR,IL,SWRK1,MVAR)
C 
C      *          CHECK TIED DATA POINTS AND GIVE THEM A SPECIAL FLAG.      *
C
       CALL TIE(IP,X,Y,NTOT,NVAR,IL,IM,SWRK1,MVAR)
C
C      *         COMPUTE SCORE VECTOR. DIMENSION IS NVAR                    *
C
       DO 600 I=1,NVAR
          SCORE(I)=0.0
C
          DO 500 J=1,NTOT 
             IF(IP(J).EQ.0) THEN
                IF(IL(J).EQ.1) THEN
                   SUM=0.0
C
                   DO 400 K=J,NTOT
                      SUM=SUM+X(I,K)
  400              CONTINUE
C
                   JJ=J+IM(J)-1
                   XSUM=0.0
                   DO 450 KL=J,JJ
                      XSUM=XSUM+X(I,KL)
  450              CONTINUE
C
                   DEN=REAL(NTOT+1-J)
                   SCORE(I)=SCORE(I)+XSUM-IM(J)*SUM/DEN
                ENDIF
             ENDIF
  500     CONTINUE
  600  CONTINUE
C
C      *    COMPUTE THE INFORMATION MATRIX. DIMENSION IS NVAR BY NVAR       *
C
       DO 1000 I=1,NVAR
          DO  900 J=I,NVAR
             RINFO(I,J)=0.0
C
             DO 800 K=1,NTOT 
                IF(IP(K).EQ.0) THEN
                   IF(IL(K).EQ.1) THEN
                      SUM1=0.0
                      SUM2=0.0
                      SUM3=0.0
C
                      DO 700 L=K,NTOT
                         SUM1=SUM1+X(I,L)
                         SUM2=SUM2+X(J,L)
                         SUM3=SUM3+X(I,L)*X(J,L)
  700                 CONTINUE
                      DEN=NTOT+1-K
                      RINFO(I,J)=RINFO(I,J)-REAL(IM(K))
     +                     *(SUM1*SUM2/DEN**2-SUM3/DEN)
                   ENDIF
                ENDIF
  800        CONTINUE
             RINFO(J,I)=RINFO(I,J)
  900     CONTINUE
 1000  CONTINUE
C
C      *     INVERT INFORMATION MATRX RINFO(I,J). THE INVERTED MATRIX       *
C      *     IS STORED IN RINFO(I,J).                                       *
C
       CALL MATINV(RINFO,NVAR,DET,IWRK1,IWRK2,MVAR)
C
C      *      COMPUTE GLOBAL CHI-SQUARE:                                    *
C      *                                                                    *
C      *       CHI = [SCORE(I)**T] X [RINFO(I,J)**-1] X [SCORE(J)]          *
C      *         WHERE T IS TRANSVERSE.                                     *
C
       DO 1200 I=1,NVAR
          FINFO(I)=0.0
          DO 1100 K=1,NVAR
             FINFO(I)=FINFO(I)+RINFO(I,K)*SCORE(K)
 1100     CONTINUE
 1200  CONTINUE
       CHI=0.0
C
       DO 1300 L=1,NVAR
          CHI=CHI+FINFO(L)*SCORE(L)
 1300  CONTINUE

C
C      *              GET THE REDUCED CHI-SQUARE                            *
C
       RCHI=CHI/DF
C
C      *           COMPUTE CHI-SQUARE PROBABILITY                           *
C
       PROB=PCHISQ(RCHI,NVAR)
C
       IF(OUTPUT.EQ.'         ') THEN
          PRINT *
          PRINT 1450
          PRINT 1400
          PRINT 1600,CHI
          PRINT 1651,NVAR
          PRINT 1650,PROB
          PRINT *
       ELSE
          WRITE(60,1400)
          WRITE(60,1450)
          WRITE(60,1400)
          WRITE(60,1600) CHI
          WRITE(60,1651) NVAR
          WRITE(60,1650) PROB
          WRITE(60,1400)
       ENDIF
 1400  FORMAT('    ')
 1450  FORMAT(5X,'CORRELATION TEST BY COX PROPORTIONAL HAZARD MODEL')
 1600  FORMAT(6X,' GLOBAL CHI SQUARE =',F9.3,' WITH ')
 1651  FORMAT(11X,I5,' DEGREES OF FREEDOM')
 1650  FORMAT(6X,' PROBABILITY    =',F10.4)

       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE DATA1 ******************************
C      **********************************************************************
C
C      *   THIS SUBROUTINE'S PURPOSE IS TO ENABLE EASY, FOOL-PROOF          *
C      *   KEYBOARD ENTRY OF INTEGER INPUT DATA.                            *
C
       SUBROUTINE DATA1(INTEG)
C
C      *               VARIABLE DECLARATIONS                                *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 B(4)
       INTEGER*4 CUTOFF,I1,I2,INTEG,TOTAL
       INTEGER*4 N(4)
C
C
C      *    READ IN NUMBER AFTER WRITING PROMPT. FIELD SIZE = 4             *
C

 3     READ(5,1) (B(I1),I1=1,4)
 1      FORMAT(4(A1))

C
C      *             ANALYZE DIGITS OF NUMBER                               *
C
       DO 2 I1=1,4
       IF(B(I1).EQ.'0') N(I1)=0
       IF(B(I1).EQ.'1') N(I1)=1
       IF(B(I1).EQ.'2') N(I1)=2
       IF(B(I1).EQ.'3') N(I1)=3
       IF(B(I1).EQ.'4') N(I1)=4
       IF(B(I1).EQ.'5') N(I1)=5
       IF(B(I1).EQ.'6') N(I1)=6
       IF(B(I1).EQ.'7') N(I1)=7
       IF(B(I1).EQ.'8') N(I1)=8
       IF(B(I1).EQ.'9') N(I1)=9
C      IF((B(I1).EQ.' ').AND.(I1.EQ.1)) PRINT *,'ENTER NUMBER.'
       IF((B(I1).EQ.' ').AND.(I1.EQ.1)) GOTO 3
C
       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
     +  (B(I1).NE.'9').AND.(B(I1).EQ.' ')) CUTOFF=I1
C
       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
     +  (B(I1).NE.'9').AND.(B(I1).EQ.' ')) GOTO 4
C
       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
     +  (B(I1).NE.'9').AND.(B(I1).NE.' '))
     +        PRINT *,'POSITIVE INTEGER ONLY'
C
       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
     +  (B(I1).NE.'9').AND.(B(I1).NE.' ')) GOTO 3
 2      CONTINUE
C
C      *                   TOTALS UP THE NUMBER                             *
C
       CUTOFF=5
 4     TOTAL=0
C
       LAST=CUTOFF-1
       DO 5 I2=1,LAST
       TOTAL=TOTAL+N(I2)*(10**(CUTOFF-I2-1))
 5     CONTINUE
C
       INTEG=TOTAL
       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE DATA2  *****************************
C      **********************************************************************
C
       SUBROUTINE DATA2(B,J,IT,INTEG,LIND)
C
C      *   PURPOSE IS TO ENABLE EASY, FOOL-PROOF KEYBOARD ENTRY OF INTEGER  *
C      *   INPUT DATA.                                                      *
C      *  INPUT                                                             *
C      *          B(4,IT)  : INPUT IN CHARACTER FORM. THIS WILL BE TESTED   *
C      *                     AND CHANGED TO INTEGER                         *
C      *          J        : J-TH INPUT, <=IT.                              *
C      *          IT       : DIMENSION OF B                                 *
C      *  OUTPUT                                                            *
C      *          INTEG    : INTEGER OUTPUT                                 *
C      *          LIND     : INDICATOR OF READABILITY OF INPUT              *
C      *                     IF LIND=0, B IS SUCCESSFULLY CHANGED TO INTEG  *
C      *                     IF LIND=1, B CANNOT BE CONVERTED TO  INTEGER   *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 B(4,IT)
       INTEGER CUTOFF,TOTAL
       INTEGER N(4)
C
C      *             ANALYZE DIGITS OF NUMBER                               *
C
        LIND=0
       DO 2 I1=1,4
       IF(B(I1,J).EQ.'0') GOTO 3 
       IF(B(I1,J).EQ.'1') GOTO 3 
       IF(B(I1,J).EQ.'2') GOTO 3 
       IF(B(I1,J).EQ.'3') GOTO 3 
       IF(B(I1,J).EQ.'4') GOTO 3 
       IF(B(I1,J).EQ.'5') GOTO 3
       IF(B(I1,J).EQ.'6') GOTO 3 
       IF(B(I1,J).EQ.'7') GOTO 3 
       IF(B(I1,J).EQ.'8') GOTO 3 
       IF(B(I1,J).EQ.'9') GOTO 3 
           IF(B(I1,J).EQ.' ') GOTO 3
           LIND=1
           GOTO 6
C
   3   IF(B(I1,J).EQ.'0') N(I1)=0
       IF(B(I1,J).EQ.'1') N(I1)=1
       IF(B(I1,J).EQ.'2') N(I1)=2
       IF(B(I1,J).EQ.'3') N(I1)=3
       IF(B(I1,J).EQ.'4') N(I1)=4
       IF(B(I1,J).EQ.'5') N(I1)=5
       IF(B(I1,J).EQ.'6') N(I1)=6
       IF(B(I1,J).EQ.'7') N(I1)=7
       IF(B(I1,J).EQ.'8') N(I1)=8
       IF(B(I1,J).EQ.'9') N(I1)=9
C
       IF((B(I1,J).EQ.' ').AND.(I1.EQ.1)) LIND=1
           IF((B(I1,J).EQ.' ').AND.(I1.EQ.1)) GOTO 6
C
        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.
     +  (B(I1,J).NE.'5').AND.   (B(I1,J).NE.'6').AND.(B(I1,J).NE.'7')
     +  .AND.(B(I1,J).NE.'8').AND.(B(I1,J).NE.'9').AND.(B(I1,J).EQ.' '))
     +            CUTOFF=I1
C
        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.
     +   (B(I1,J).NE.'5').AND.(B(I1,J).NE.'6').AND.(B(I1,J).NE.'7')
     +   .AND.(B(I1,J).NE.'8').AND.(B(I1,J).NE.'9')
     +   .AND.(B(I1,J).EQ.' ')) GOTO 4
C     
        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.(B(I1,J).NE.'5')
     +   .AND.(B(I1,J).NE.'6').AND.(B(I1,J).NE.'7').AND.(B(I1,J).NE.'8')
     +  .AND.(B(I1,J).NE.'9').AND.(B(I1,J).NE.' ')) LIND=1   
 2      CONTINUE
C
C      *                  TOTAL UP NUMBER                                   *
C
       CUTOFF=5
 4      TOTAL=0
C
       LAST=CUTOFF-1
       DO 5 I2=1,LAST
       TOTAL=TOTAL+N(I2)*(10**(CUTOFF-I2-1))
 5     CONTINUE
       INTEG=TOTAL
6      RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE DATREG  ****************************
C      **********************************************************************
C
       SUBROUTINE DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,
     +                   NC5,NC6,NC7,NC8,ICENS,NYC,NXC,NBC,NDIM,NDATA)
C
C      * THIS SUBROUTINE READS DATA FROM THE FILE FOR CORRELATION/REGRESSION*
C      * PROBLEMS                                                           *
C      *                                                                    *
C      * INPUT        NVAR     : NUMBER OF THE INDEPENDENT VARIABLE         *
C      *              FILE     : NAME OF THE DATA FILE                      *
C      *                                                                    *
C      * OUTPUT       IND(1,I) : INDICATOR OF CENSORSHIP                    *
C      *              X(J,I)   : INDEPENDENT VARIABLES                      *
C      *              Y(I)     : DEPENDENT VARIABLES                        *
C      *              NTOT     : NUMBER OF DATA POINTS                      *
C      *              ND       : NUMBER OF DETECTED POINTS                  *
C      *              NC1      : NUMBER OF Y LOWER LIMITS                   *
C      *              NC2      : NUMBER OF X LOWER LIMITS                   *
C      *              NC3      : NUMBER OF DOUBLE LOWER LIMITS              *
C      *              NC4      : NUMBER OF Y LOWER, X UPPER LIMITS          *
C      *              NC5      : NUMBER OF Y UPPER LIMITS                   *
C      *              NC6      : NUMBER OF X UPPER LIMITS                   *
C      *              NC7      : NUMBER OF DOUBLE UPPER LIMITS              *
C      *              NC8      : NUMBER OF Y UPPER, X LOWER LIMITS          *
C      *              ICENS    : INDICATOR OF CENSORING                     *
C      *              NYC      : NC1+NC5                                    *
C      *              NXC      : NC2+NC6                                    *
C      *              NBC      : NC3+NC4+NC7+NC8                            *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(NDIM, NDATA),X(NDIM, NDATA),Y(NDATA)
       CHARACTER*9 FILE
C
       OPEN(UNIT=40,FILE=FILE, STATUS='OLD', FORM='FORMATTED')
C
       NTOT=1
       ND  =0
       NC1 =0
       NC2 =0
       NC3 =0
       NC4 =0
       NC5 =0
       NC6 =0
       NC7 =0
       NC8 =0
C
   10  READ(40,20,END=30) IND(1,NTOT),(X(J,NTOT),J=1,NVAR),Y(NTOT)
   20  FORMAT(I4,11F10.3)
C
       IF(IND(1,NTOT).EQ.0)  ND =ND+1
       IF(IND(1,NTOT).EQ.1)  NC1=NC1+1
       IF(IND(1,NTOT).EQ.2)  NC2=NC2+1
       IF(IND(1,NTOT).EQ.3)  NC3=NC3+1
       IF(IND(1,NTOT).EQ.4)  NC4=NC4+1
       IF(IND(1,NTOT).EQ.-1) NC5=NC5+1
       IF(IND(1,NTOT).EQ.-2) NC6=NC6+1
       IF(IND(1,NTOT).EQ.-3) NC7=NC7+1
       IF(IND(1,NTOT).EQ.-4) NC8=NC8+1
       NTOT=NTOT+1

       IF(NTOT.GT.NDATA) THEN
          WRITE(6,12) NDATA
          WRITE(6,14)
   12     FORMAT(' ****WARNING!!  THERE ARE MORE THAN',I5,' POINTS. ')
   14     FORMAT(' ARRAYS WILL OVERFLOW UNLESS DIMENSIONS ARE CHANGED')
          STOP       
       ENDIF
C
C      * THE NEXT BLOCK OF CODE CHECKS FOR LINES WHICH HAVE ALL INPUTS     *
C      * EQUAL TO ZERO.  THIS WOULD BE THE EFFECT OF A BLANK LINE IN THE   *
C      * DATA FILE, AND COULD RESULT IN INCORRECT VALUES COMPUTED.         *
C
       K = 0
       IF(IND(1,NTOT-1) .NE. 0) K =1
       IF(Y(NTOT-1) .NE. 0.0) K = 1
       DO 21 J = 1, NVAR
          IF(X(J,NTOT-1) .NE. 0.0) K = 1
 21    CONTINUE
       IF(K .EQ. 0) THEN
          WRITE(6,22)
          WRITE(6,24) NTOT
       ENDIF

 22    FORMAT('         ')
 24    FORMAT('WARNING: LINE ',I4,' IN THE DATA FILE MAY BE BLANK')

       GOTO 10
C
   30  NTOT=NTOT-1
       NYC=NC1+NC5
       NXC=NC2+NC6
       NBC=NC3+NC4+NC7+NC8
       NLC=NC1+NC2+NC3+NC4
       NUC=NC5+NC6+NC7+NC8
       ICENS=0
       IF(NLC.EQ.0) ICENS=-1
       IF(NUC.EQ.0) ICENS=1
C
       CLOSE(UNIT=40)
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE EM  *********************************
C      **********************************************************************
C
       SUBROUTINE EM(IND,X,Y,NTOT,TOL,MAXITS,NVAR,IBET,ND,NC,OUTPUT,
     +               FILE,ALPHA,XX,Y2,W,WCEN,VCOV,WORK,SIGMAA,TOLA,
     +               LENG,LEGWRK,MVAR)
C
C      *           MAXIMUM LIKELIHOOD ESTIMATION IN A LINEAR MODEL          *
C      *           FROM CONFINED AND CENSORED NORMAL DATA                   *
C      *                                                                    *
C      *       FORMAL PARAMETERS :                                          *
C      *                                                                    *
C      *       NTOT    INPUT : THE NUMBER OF OBSERVATIONS                   *
C      *      MPLONE   INPUT : THE TOTAL NUMBER OF PARAMETERS TO BE         *
C      *                       ESTIMATED (I.E. NB+1)                        *
C      *      MAXITS   INPUT : THE MAXIMUM NUMBER OF ITERATIONS ALLOWED     *
C      *      IBET     INPUT : IF THERE IS DATA WHICH IS CONFINED           *
C      *                       BETWEEN TWO VALUES, IBET=1 OTHERWISE 0       *
C      *      ALPHA    INPUT : IF ALPHA(NVAR+2).LE.0.0, THE SUBROUTINE      *
C      *                       CALCULATES INITIAL PARAMETER ESTIMATES.      *
C      *                       IF >0.0, IT CONTAINS THE INITIAL             *
C      *                       ESTIMATE OF THE J-TH LOCATION PARAMETER      *
C      *                       FOR J=1,2,.....,NB.                          *
C      *               OUTPUT: THE MOST RECENT PARAMETER ESTIMATES          *
C      *                       BEFORE EXIT FROM THE SUBROUTINE.             *
C      *      TOL      INPUT : CONVERGENCE TO THE MAXIMUM LIKELIHOOD        *
C      *                       PARAMETER ESTIMATE IS REACHED IF THE         *
C      *                       DIFFERENCE BETWEEN CONSECUTIVE ESTIMATES     *
C      *                       OF THE J-TH PARAMETER IS LESS THAN TOL(J)    *
C      *                       FOR J=1,2,........,M.                        *
C      *       IND     INPUT : INDICATOR OF CENSORED DATA                   *
C      *       XX      INPUT : THE DESIGN MATRIX XX(I,J) CONTAINS THE       *
C      *                       COEFFICIENT OF THE J-TH LOCATION PARA-       *
C      *                       METER FOR I-TH OBSERVATION.                  *
C      *       Y       INPUT : IF PP(I)=0, THE I-TH OBSERVATION IS          *
C      *                       COMPLETELY SPECIFIED IN Y(I); IF             *
C      *                       PP(I)=-1, Y(I)  IS LEFT-CENSORED: IF         *
C      *                       PP(I)=1, Y(I) IS RIGHT-CENSORED: IF          *
C      *                       PP(I)=5, Y(I) IS CONFINED BETWEEN TWO        *
C      *                       VALUES.                                      *
C      *      ROWX     WORK  : THE NUMBER OF ROWS OF X (=NTOT)              *
C      *      COLX     WORK  : THE NUMBER OF COLUMNS OF X (=NB)             *
C      *       W       WORK  :                                              *
C      *      LENW     WORK  :   = NB+NTOT                                  *
C      *      WORK     WORK  :                                              *
C      *      LENWRK   WORK  :   = NB*NTOT                                  *
C      *      LEG            : DIMENSION SIZE >= LENW                       *
C      *      LEGWRK         : DIMENSION SIZE >= LWNWRK                     *
C      *      MVAR           : DIMENSION SIZE >= MPLONE                     *
C      *      VCOV     OUTPUT: IF THE PROCEDURE CONVERGED TO THE MAX-       *
C      *                       IMUM LIKELIHOOD ESTIMATES, THE FIRST         *
C      *                       (NB+1)*(NB+1) POSITIONS CONTAIN AN ESTI-     *
C      *                       MATE OF THE VARIANCE-COVARIANCE MATRIX       *
C      *                       OF THESE ESTIMATES.                          *
C      *      SIGMAA   OUTPUT: (J,J) COMPONENT CONTAINS THE STANDARD        *
C      *                       DEVIATION OF THE J-TH PARAMETER.             *
C      *      IFAULT   OUTPUT: FAILURE INDICATOR                            *
C      *      ICHECK   OUTPUT: IF ERROR ANALYSIS IS NOT COMPLETED,          *
C      *                       ICHECK HAS A VALUE : 1 ,OTHERWISE 0 .        *
C      *                                                                    *
C      *       SUBROUTINES                                                  *
C      *                        EMALGO                                      *
C      *                                                                    *
C      *       REF : M.S.WOLYNETZ AS 139 APL.STATIST.VOL.28 195 (1979)      *
C      *             PLUS CORRECTIONS IN LATER ISSUES OF APPLIED STATISTICS.*
C      *                                                                    *
C      *       SUBROUTINE EMALGO IS COPIED DIRECTLY FROM WOLYNETZ, EXCEPT   *
C      *       FOR A FEW CHANGES TO PRINT AND COMMENT STATEMENTS.           *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
       INTEGER   ROWX,COLX

       DIMENSION X(MVAR,NTOT),Y(NTOT),IND(NTOT),XX(NTOT,MVAR),Y2(NTOT)
       DIMENSION W(LENG),WCEN(LENG),VCOV(LEGWRK),WORK(LEGWRK)
       DIMENSION ALPHA(MVAR),SIGMAA(MVAR,MVAR),TOLA(MVAR)
       CHARACTER*9 FILE,OUTPUT
C
C
C      *                  INITIALIZATION                                    *
C
       MPLONE=NVAR+2
       NB=NVAR+1
       ROWX=NTOT
       COLX=MPLONE
       LENW=NTOT+MPLONE
       LENWRK=NTOT*MPLONE
C
       DO 10 I=1,MPLONE
          TOLA(I)=TOL
   10  CONTINUE
C
C      *                                                                    *
C      *   IF THERE ARE DATA POINTS CONFINED BETWEEN TWO VALUES, READ THE   *
C      *   INPUT DATA AGAIN.                                                *
C      *                                                                    *
       IF(IBET.NE.0) THEN

          OPEN(UNIT=40,FILE=FILE,STATUS='OLD',FORM='FORMATTED')

          DO 40 I=1,NTOT
             READ(40,30) IND(I),(XX(I,J),J=2,NB),Y(I),Y2(I)
   30        FORMAT(I3,20F10.3)
             XX(I,1) = 1.0
   40     CONTINUE
       CLOSE(UNIT=40)
C
       ELSE
          DO 100 I=1,NTOT
             XX(I,1)=1.0
             DO 90 J=1,NVAR
                JJ=J+1
                XX(I,JJ)=X(J,I)
   90        CONTINUE
  100     CONTINUE
       ENDIF
C
       ICHECK=0
C
C      *   CALL THE SUBPROGRAM FOR THE MAXIMUM LIKELIHOOD ESTIMATES         *
C
       CALL EMALGO(NTOT,Y,Y2,IND,MPLONE,XX,ROWX,COLX,W,WCEN,LENW,
     +      VCOV,WORK,LENWRK,ALPHA,TOLA,MAXITS,IFAULT,ICHECK,NC)
C
C
C      *    PRINT OUT THE RESULTS                                           *
C
C
       IF(IFAULT.GT.0) THEN
          IF(ICHECK.EQ.0) THEN
C
C      *        CHANGE VCOV ARRAY FROM ONE DIM. TO TWO DIM.                 *
C
             K=1
             DO 300 I=1,MPLONE
                DO 200 J=1,MPLONE
                   SIGMAA(I,J)=DSQRT(DABS(VCOV(K)))
                   K=K+1
  200           CONTINUE
  300        CONTINUE
          ENDIF

  360     IF(OUTPUT.EQ.'         ') THEN
C
C      PRINT FINAL REGRESSION COEFFICIENTS AT THE TERMINAL. 
C
             PRINT 1050
             PRINT 1020         
             PRINT 1050
             PRINT 1200,ALPHA(1),SIGMAA(1,1)

             DO 450 J=2,NB
                JI=J-1
                PRINT 1250,JI,ALPHA(J),SIGMAA(J,J)
  450        CONTINUE

             PRINT 1300,ALPHA(MPLONE)
             IF(ICHECK.EQ.0) THEN
                ITE=IFAULT
             ELSE
                ITE=ICHECK
             ENDIF
             PRINT 1350,ITE
             PRINT 1050
          ELSE
C
C      PRINT FINAL REGRESSION COEFFICIENTS IN THE OUTPUT FILE.
C
             WRITE(60,1050)
             WRITE(60,1020)        
             WRITE(60,1050)
             WRITE(60,1200) ALPHA(1),SIGMAA(1,1)
   
             DO 455 J=2,NB
                JI=J-1
                WRITE(60,1250) JI,ALPHA(J),SIGMAA(J,J)
  455        CONTINUE
             WRITE(60,1300) ALPHA(MPLONE)
             IF(ICHECK.EQ.0) THEN
                ITE=IFAULT
             ELSE
                ITE=ICHECK
             ENDIF
             WRITE(60,1350) ITE
             WRITE(60,1050)
          ENDIF
C
C      IF AN ERROR OCCURED DURING THE COMPUTATION, PRINT OUT AN ERROR
C      MESSAGE.
C
C      IN THE FOLLOWING, WE HAVE COLLECTED ALL OF WOLYNETZ'S ERROR FLAGS
C      AND PRINT OUT THE APPROPRIATE ERROR MESSAGE. 
C

       ELSE

C
          IF(OUTPUT.EQ.'         ') THEN
C
C      PRINT FINAL REGRESSION COEFFICIENTS ON THE TERMINAL. 
C
             PRINT 1000
             PRINT 1020
             PRINT 1050
             PRINT 2001
             IF(IFAULT.EQ.-1) THEN
                PRINT 2002
                PRINT 2003
             ELSEIF(IFAULT.EQ.-4) THEN
                PRINT 2012
                PRINT 2014
             ELSEIF(IFAULT.EQ.-5) THEN
                PRINT 2021
                PRINT 2022
                PRINT 2023
                PRINT 2024
                PRINT 2025
                PRINT 2026
                PRINT 2027
             ELSEIF(IFAULT.EQ.-6) THEN
                PRINT 2031
                PRINT 2032
                PRINT 2033
                PRINT 2034
                PRINT 2035
                PRINT 2036
                PRINT 2037
             ELSEIF(IFAULT.EQ.-7) THEN
                PRINT 2045
             ELSEIF(IFAULT.EQ.-8) THEN
                PRINT 2055
             ELSE
                PRINT 2070
             ENDIF
C
       ELSE
C
C      PRINT FINAL REGRESSION COEFFICIENTS IN THE OUTPUT FILE.
C
             WRITE(60,1000)
             WRITE(60,1020)
             WRITE(60,1050)
             WRITE(60,2001)
             IF(IFAULT.EQ.-1) THEN
                WRITE(60,2002)
                WRITE(60,2003)
             ELSEIF(IFAULT.EQ.-4) THEN
                WRITE(60,2012)
                WRITE(60,2014)
             ELSEIF(IFAULT.EQ.-5) THEN
                WRITE(60,2021)
                WRITE(60,2022)
                WRITE(60,2023)
                WRITE(60,2024)
                WRITE(60,2025)
                WRITE(60,2026)
                WRITE(60,2027)
             ELSEIF(IFAULT.EQ.-6) THEN
                WRITE(60,2031)
                WRITE(60,2032)
                WRITE(60,2033)
                WRITE(60,2034)
                WRITE(60,2035)
                WRITE(60,2036)
                WRITE(60,2037)
             ELSEIF(IFAULT.EQ.-7) THEN
                WRITE(60,2045)
             ELSEIF(IFAULT.EQ.-8) THEN
                WRITE(60,2055)
             ELSE
                WRITE(60,2070)
             ENDIF
          ENDIF
C
       ENDIF
C
 1000  FORMAT(T10,'REGRESSION ANALYSIS WITH CENSORED DATA')
 1020  FORMAT(T5,'LINEAR REGRESSION BY PARAMETRIC EM ALGORITHM')
 1050  FORMAT(T5,'      ')
 1100  FORMAT(T8,'DATA TITLE :',T25,60A1)
 1150  FORMAT(T8,'TOTAL # OF DATA :',T26,I3,T33,'CENSORED DATA :'
     +           ,T48,I3)
 1200  FORMAT(T8,'INTERCEPT COEFF    :',F8.4,T38,'+/-',T41,F8.4)
 1250  FORMAT(T8,'SLOPE COEFF ',I1,'      :',F8.4,T38,'+/-',T41,
     +         F8.4,5X)
 1300  FORMAT(T8,'STANDARD DEVIATION :',F8.4)
 1350  FORMAT(T8,'ITERATIONS         :',I3)
 2001  FORMAT(T5,'NOTICE :')
 2002  FORMAT(T5,'MAXIMUM NUMBER OF ITERATION REACHED AND')
 2003  FORMAT(T5,'CONVERGENCE HAS NOT BEEN OBTAINED.')
 2012  FORMAT(T5,'NUMBER OF COMPLETELY SPECIFIED DATA IS LESS')
 2014  FORMAT(T5,'THAN NB+1')
 2021  FORMAT(T5,'THE MATRIX IS NOT POSITIVE DEFINITE,AS')
 2022  FORMAT(T5,'DETERMINED BY SUBROUTINE "SYMINV",A MATRIX')
 2023  FORMAT(T5,'INVERSION PROCEDURE(HEALY,1968); THE VALUE')
 2024  FORMAT(T5,'OF "NULLTY" AND "IFAULT" RETURNED BY SYMINV')
 2025  FORMAT(T5,'ARE PLACED IN THE FIRST TWO POSITIONS OF ')
 2026  FORMAT(T5,'THE ARRAY "VCOV" BEFORE RETURNING TO THE ')
 2027  FORMAT(T5,'CALLING PROGRAM.')
 2031  FORMAT(T5,'THE ESTIMATE OF THE VARIANCE- COVARIANCE')
 2032  FORMAT(T5,'MATRIX IS NOT POSITIVE DEFINITE, AS DETER-')
 2033  FORMAT(T5,'MINED BY SUBROUTINE "SYMINV"(HEALY,1968):')
 2034  FORMAT(T5,'THE VALUES OF "NULLTY" AND "IFAULT", RETURNED')
 2035  FORMAT(T5,'BY SYMINV ,ARE PLACED IN THE FIRST TWO ')
 2036  FORMAT(T5,'POSITIONS OF THE ARRY "VCOV" BEFORE RETURNING')
 2037  FORMAT(T5,'TO THE CALLING PROGRAM')
 2045  FORMAT(T5,'"ROWX" IS LESS THAN NTOT')
 2055  FORMAT(T5,'"COLX" IS LESS THAN NB')
 2070  FORMAT(T5,'THE PROGRAM IS STOPPED FOR UNKNOWN REASONS')

       RETURN
       END

C
C      **********************************************************************
C      *********************** FUNCTION FACTOR  *****************************
C      **********************************************************************
C
       FUNCTION FACTOR(N)
C
C      *     COMPUTES THE FACTORIAL FUNCTION FOR INTEGERS.                  *
C      *     THIS FUNCTION IS BASED ON PROG. 3-2 ON P. 32 OF "DATA REDUCTION*
C      *     AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", PHILIP R.       *
C      *     BEVINGTON, 1969, McGRAW HILL (NY:NY)                           *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
       FACTOR=1.0
       IF(N-1 .GT. 0) THEN
   13     IF(N-10 .LE. 0) THEN
C
C      *                N LESS THAN 11                                      *
C
   21        DO 23 I=2,N
                FI=I
                FACTOR=FACTOR*FI
   23        CONTINUE
C
C      *               N GREATER THAN 10                                    *
C
          ELSE
   31        SUM=0.0
             DO 34 I=11,N
                FI=I
                SUM=SUM+DLOG(FI)
   34        CONTINUE
             FACTOR=3628800.0D00*DEXP(SUM)
          ENDIF
       ENDIF
   40  RETURN
       END
C
C      **********************************************************************
C      ********************** SUBROUTINE DATAIN  ****************************
C      **********************************************************************
C
       SUBROUTINE DATAIN(IUNI,FILE,NVAR,ISTA,IND,X,NTOT,NDATA,MVAR)
C
C      *  THIS SUBROUTINE READS DATA FOR UNIVARIATE AND TWO SAMPLE PROBLEMS *
C      *                                                                    *
C      *  INPUT                                                             *
C      *          INUI  : INDICATOR OF WHICH PROBLEM (KM OR TWO SAMPLE) IS  *
C      *                  NEEDED                                            *
C      *          FILE  : DATA FILE NAME                                    *
C      *          NVAR  : NUMBER OF VARIABLES                               *
C      *          NDATA : DIMENSION FOR THE VARIABLES                       *
C      *  OUTPUT                                                            *
C      *         ISTA(I): INDICATOR OF GROUPS                               *
C      *        IND(J,I): INDICATOR OF CENSORSHIP                           *
C      *          X(J,I): VARIABLES                                         *
C      *         NTOT   : NUMBER OF DATA POINTS                             *
C      *                                                                    *
C      *    THIS FILE WAS MODIFIED ON 4/13/90                               *
C      *       NDATA WAS ADDED FOR THE DIMENSION DECLARATION.              *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(MVAR, NDATA), X(MVAR,NDATA), ISTA(NDATA)
       CHARACTER*9 FILE
C
C      *        OPEN THE DATA FILE AND READ DATA                            *
C
       OPEN(UNIT=40, FILE=FILE, STATUS='OLD', FORM='FORMATTED')
C
       J=0
       IF(IUNI .EQ. 1) THEN
C
C      *            K-M ESTIMATOR FORMAT                                    *
C
   10     J=J+1
          READ(40,30,END=50) (IND(I,J),X(I,J),I=1,NVAR)
          K=0
          DO 15 I = 1,NVAR
             IF(IND(I,J).NE.0.OR.X(I,J).NE.0.0) K = 1
 15       CONTINUE
          IF(K .EQ. 0) THEN
             WRITE(6,44)
             WRITE(6,45) J
          ENDIF
          GOTO 10
C
C      *            TWO SAMPLE TEST FORMAT                                  *
C
       ELSE
   20     J=J+1
          READ(40,40,END=50) ISTA(J),(IND(I,J),X(I,J),I=1,NVAR)
          K = 0
          IF(ISTA(J) .NE. 0) K=1
          DO 25 I = 1, NVAR
             IF(IND(I,J) .NE. 0 .OR. X(I,J) .NE. 0.0) K = 1
 25       CONTINUE
          IF(K .EQ. 0) THEN
             WRITE(6,44)
             WRITE(6,45) J
          ENDIF
          GOTO 20
       ENDIF
C
   30  FORMAT(10(I4,F10.3))
   40  FORMAT(I4,10(I4,F10.3))
   44  FORMAT('         ')  
   45  FORMAT('WARNING: LINE ',I4,
     +         ' IN THE DATA FILE MAY BE BLANK') 
   50  NTOT=J-1
       CLOSE(UNIT=40)
     
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE EMALGO  *****************************
C      **********************************************************************
C
       SUBROUTINE EMALGO(NTOT,Y1,Y2,P,MPLONE,X,ROWX,COLX,W,WCEN,LENW,
     +          VCOV,WORK,LENWRK,ALPHA,TOL,MAXITS,IFAULT,ICHECK,NC)
C
C
C      *       ALGORITHM AS 139 APPL.STATIST. (1979) VOL.28., NO2           *
C      *                                                                    *
C      *       COMPUTES MAXIMUM LIKELIHOOD ESTIMATES                        *
C      *       FROM A LINEAR MODEL WITH NORMAL HETEROGENEOUS                *
C      *       VARIANCE. THE DESIGN MATRIX MUST BE NON-SINGULAR.            *
C      *       THE DEPENDENT VARIABLE MAY INCLUDE OBSERVATIONS              *
C      *       CENSORED IN EITHER TAIL AND/OR OBSERVATIONS CONFINED         *
C      *       BETWEEN FINITE LIMITS.                                       *
C      *                                                                    *
C      *  SUBROUTINES                                                       *
C      *             SYMINV, UNPACK, RMILLS                                 *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       INTEGER ROWX,COLX,P

       DIMENSION X(NTOT,MPLONE),TOL(MPLONE),ALPHA(MPLONE)
       DIMENSION Y1(NTOT),Y2(NTOT),P(NTOT)
       DIMENSION VCOV(LENWRK),WORK(LENWRK),W(LENW),WCEN(LENW)
       DATA QLIMIT /0.00001/, RLIMIT /0.00001/
       DATA C /0.39894228/
C
C      *                    INITIALIZATION                                  *
C      *       THE NEXT LINE IS LOCATED IN A DIFFERENT PLACE IN THE         *
C      *       ORIGINAL PROGRAM BY WOLYNETZ                                 *
C
       M=MPLONE-1
C
C      *                 CHECK ARRAY SIZES, ETC                             *
C
       IFAULT=-7
       IF(ROWX.LT.NTOT) RETURN
       IFAULT=-8
       IF(COLX.LT.M) RETURN
       IFAULT=-9
C
C
C      *       THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 31       *
C      *       APPL.STAT. VOL.29 P.228 , (1980).                            *
C
C
       IF(LENW.LT.(MPLONE+NTOT)) RETURN
       IFAULT=-10
       IF(LENWRK.LT.(M*NTOT)) RETURN
C
C      *          COMPUTE X'X IN LOWER TRIANGULAR FORM                      *
C
       II=0
       DO 53 I=1,M
          DO 50 J=1,I
             TEMP=0.0
             DO 40 K=1,NTOT
                TEMP=TEMP+X(K,I)*X(K,J)
   40        CONTINUE
             II=II+1
             VCOV(II)=TEMP
   50     CONTINUE
   53  CONTINUE

       CALL SYMINV(VCOV,M,WORK,W,NUL,LENWRK,LENWRK,LENW,IFAULT)

       IF((IFAULT.NE.0) .OR. (NUL.NE.0)) THEN
          VCOV(2)=REAL(IFAULT)
          VCOV(1)=REAL(NUL)
          IFAULT=-5
          RETURN
       ENDIF

C
C      *      THE MATRIX IS NON-SINGULAR AND WE HAVE ITS INVERSE.  NOW WE   *
C      *      COMPUTE (X'X) INVERSE*X.                                      *
C      *      THE FOLLOWING SCHEME IS USED TO REDUCE THE NUMBER OF STORAGE  *
C      *      ARRAYS NEEDED, BY EXPANDING FROM THE TRIANGULAR TO A SQUARE   *
C      *      MATRIX.                                                       *
C
C      *         THE NEXT LINE IS NOT IN WOLYNETZ.                          *
C
       IF(M.NE.1) THEN
          CALL UNPACK(WORK,M,LENWRK)
       ENDIF
C
C      *      DO MULTIPLICATION--ONE ROW AT A TIME--STARTING WITH           *
C      *      THE LAST ROW                                                  *
C
       JJ=NTOT*M
       II=M*M
       DO 220 I=1,M
          II=II-M

          DO 200 J=1,NTOT
             TEMP=0.0

             DO 170 K=1,M
                IIK=II+K
                TEMP=TEMP+WORK(IIK)*X(J,K)
  170        CONTINUE
             W(J)=TEMP
  200     CONTINUE

          DO 210 J=1,NTOT
             IJ=NTOT+1-J
             WORK(JJ)=W(IJ)
             JJ=JJ-1
  210     CONTINUE
  220  CONTINUE
C
       XSIG=ALPHA(MPLONE)
       IF(XSIG.LE.0.0) THEN
C
C      *       NO ACCEPTABLE INITIAL VALUE FOR SIGMA HAS BEEN INPUT,        *
C      *       OBTAIN INITIAL ESTIMATES FROM EXACTLY SPECIFIED              *
C      *       OBSERVATIONS ONLY (ALTHOUGH THE MATRIX IS BASED ON ALL       *
C      *       OBSERVATIONS) AND CONFINED OBSERVATIONS.                     *
C
          II=-NTOT
          DO 300 I=1,M
             II=II+NTOT
             TEMP=0.0
             DO 280 J=1,NTOT
                IIJ=II+J
                IPT=P(J)
                IF(IPT.EQ.0) THEN
                   TEMP=TEMP+WORK(IIJ)*Y1(J)
                ELSEIF(IPT.EQ.5) THEN 
                   TEMP=TEMP+WORK(IIJ)*(Y1(J)+Y2(J))*0.5
                ENDIF
  280        CONTINUE
             ALPHA(I)=TEMP
  300     CONTINUE
C
C      *           CALCULATE THE INITIAL ESTIMATE OF SIGMA                  *
C
          SUM2=0.0
          TEMP=0.0
          DO 350 I=1,NTOT
             IPT=P(I)
             IF(IABS(IPT).NE.1) THEN
                DEMP=Y1(I)
                IF(IPT.EQ.5) DEMP=(DEMP+Y2(I))*0.5

                DO 320 J=1,M
                   DEMP=DEMP-ALPHA(J)*X(I,J)
  320           CONTINUE

                SUM2=SUM2+DEMP**2
                TEMP=TEMP+1.0
             ENDIF
  350     CONTINUE

          XSIG=DSQRT(SUM2/TEMP)
       ENDIF
C
C      *         COMPUTE SOME CONSTANTS NEEDED THROUGHOUT THE SUBROUTINE    *
C
       R=0.0
       R2=0.0
       IFAULT=-2
       DO 600 I=1,NTOT
          IPT=P(I)
          IF(IPT.EQ.0)  THEN
             R=R+1.0
             W(I)=Y1(I)
C
C      *       THE NEXT LINE IS LOCATED IN A DIFFERENT PLACE IN THE         *
C      *       ORIGINAL PROGRAM BY WOLYNETZ                                 *
C
          ELSEIF(IPT.EQ.5) THEN
             IF(DABS(Y1(I)-Y2(I)) .GT. QLIMIT*DABS(Y1(I))) THEN
                R2=R2+1.0
                IF(Y1(I).LT.Y2(I)) GOTO 600
                RETURN
             ENDIF
             P(I)=0
             R=R+1.0
             W(I)=Y1(I)
          ENDIF
  600  CONTINUE

       I=INT(R+R2+0.01)
       IFAULT=-4
       IF(I.LT.MPLONE) RETURN
       IFAULT=0
C
C      *             START OF THE ITERATION PROCEDURE                       *
C
  620  TD=R
       SUM2=0.0
C
C      *             COMPLETE W-VECTOR                                      *
C
       DO 1000 I=1,NTOT
          IPT=P(I)
          YMEAN=0.0

          DO 650 J=1,M
             YMEAN=YMEAN+ALPHA(J)*X(I,J)
  650     CONTINUE
C
C      *       OBSERVATION IS NOT EXACTLY SPECIFIED: START FROM LINE 990    *
C
          IF(IPT.NE.0) THEN
             TEMP=(Y1(I)-YMEAN)/XSIG
C
C      *        OBSERVATION CENSORED FROM ABOVE:  LOWER BOUND IS KNOWN      *
C
             IF((IPT-1) .EQ. 0) THEN
  700           CALL RMILLS(TEMP,F,RLIMIT)
                W(I)=YMEAN+XSIG*F
                TD=TD+F*(F-TEMP)
C
C      *         OBSERVATON CENSORED FROM BELOW:  UPPER BOUND IS KNOWN      *
C
             ELSEIF((IPT-1) .LT. 0) THEN
  750           CALL RMILLS(-TEMP,F,RLIMIT)
                W(I)=YMEAN-XSIG*F
                TD=TD+F*(F+TEMP)
C
C      *       OBSERVATION CONFINED TO LIE BETWEEN TWO FINITE LIMITS        *
C
             ELSEIF((IPT-1) .GT. 0) THEN
  800           YN=DEXP(-0.5*TEMP**2)*C
                CALL RMILLS(TEMP,F,RLIMIT)
                YQ=YN/F
                TMPU=(Y2(I)-YMEAN)/XSIG
                YNU=DEXP(-0.5*TMPU**2)*C
                CALL RMILLS(TMPU,FU,RLIMIT)
                YQU=YNU/FU
                TINT=YQ-YQU

                IF(TINT.LT.QLIMIT) THEN
                   IFAULT=-3
                   RETURN
                ENDIF
C
C      *       AFTER STANDARDIZING, UPPER AND LOWER LIMITS RESULT IN        *
C      *       THE SAME PROBABILITY INTEGRAL                                *
C
  820           A=(YN-YNU)/TINT
                W(I)=YMEAN+XSIG*A
                TD=TD+(A**2+(TMPU*YNU-TEMP*YN)/TINT)
             ENDIF
          ENDIF
C
C      *        CALCULATE RESIDUAL SUM OF SQUARES                           *
C
  990     SUM2=SUM2+(W(I)-YMEAN)**2
 1000  CONTINUE
C
C      *    UPDATE PARAMETER ESTIMATES-STORE IN THE END OF THE W-VECTOR     *
C
       JJ=-NTOT
       DO 1200 J=1,M
          JJ=JJ+NTOT
          TEMP=0.0

          DO 1100 I=1,NTOT
             JJI=JJ+I
             TEMP=TEMP+WORK(JJI)*W(I)
 1100     CONTINUE
          NJ=NTOT+J
          W(NJ)=TEMP
 1200  CONTINUE

       NJ=NTOT+MPLONE
       W(NJ)=DSQRT(SUM2/TD)
C
C      *       STORE THE ESTIMATES FOR THE CENSORED POINTS                  *
C
       KC=0
       DO 1250 I=1,NTOT
          IF(P(I).EQ.0) GOTO 1250
          KC=KC+1
          WCEN(KC)=W(I)
 1250  CONTINUE
C
C      *             TEST FOR CONVERGENCE                                   *
C
       DO 1300 J=1,MPLONE
          NJ=NTOT+J
          IF(DABS(ALPHA(J)-W(NJ)).GE.TOL(J)) GOTO 1400
 1300  CONTINUE
C
C      *          IF WE REACH HERE, CONVERGENCE IS OBTAINED                 *
C
       IJ=IFAULT
       IFAULT=-1
C
C      *                UPDATE VALUES                                       *
C
 1400  DO 1450 J=1,MPLONE
          NJ=NTOT+J
          ALPHA(J)=W(NJ)
 1450  CONTINUE

       XSIG=ALPHA(MPLONE)
       IFAULT=IFAULT+1
       IF(IFAULT.NE.0) THEN
C
C      *      IF THE NUMBER OF ITERATIONS HAS NOT REACHED THE MAX., TRY     *
C      *      ANOTHER ITERATION.                                            *
C
          IF(IFAULT.LE.MAXITS) GOTO 620
          IFAULT=-1
          RETURN
       ENDIF
C
C      *        CONVERGENCE OBTAINED.  COMPUTE VARIANCE--COVARIANCE         *
C      *        MATRIX, AND INITIALIZE THE WORK ARRAY                       *
C
 1600  II=MPLONE*(MPLONE+1)/2
       DO 1650 I=1,II
          WORK(I)=0.0
 1650  CONTINUE

       DO 2500 I=1,NTOT
          IPT=P(I)
          YS=Y1(I)

          DO 1680 J=1,M
             YS=YS-ALPHA(J)*X(I,J)
 1680     CONTINUE

          YS=YS/XSIG
          JJ=0
          IF(IPT.EQ.0) THEN

C
C      *             EXACTLY SPECIFIED OBSERVATION                          *
C

             DO 1750 K=1,M
                DO 1720 J=1,K
                   JJ=JJ+1
                   WORK(JJ)=WORK(JJ)+X(I,K)*X(I,J)
 1720           CONTINUE
C
C
C      *      THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32        *
C      *      APPL.STAT. VOL 30, P. 105 (1981).                             *
C
C
                KK=II-MPLONE+K
                WORK(KK)=WORK(KK)+YS*X(I,K)
 1750        CONTINUE
             WORK(II)=WORK(II)+1.0+YS**2

C
C      *      OBSERVATION CENSORED FROM ABOVE:  LOWER BOUND IS KNOWN        *
C

          ELSEIF((IPT-1) .LE. 0) THEN
             IF((IPT-1) .EQ. 0) THEN
                CALL RMILLS(YS,F,RLIMIT)
                TEMP=F*(F-YS)

C
C      *      OBSERVATION CENSORED FROM BELOW:  UPPER BOUND IS KNOWN        *
C

             ELSE
                CALL RMILLS(-YS,F,RLIMIT)
                TEMP=F*(F+YS)
             ENDIF
C
C      *         ROUTINE FOR CENSORED OBSERVATIONS                          *
C
             DO 2190 K=1,M
                DO 2170 J=1,K
                   JJ=JJ+1
                   WORK(JJ)=WORK(JJ)+X(I,J)*X(I,K)*TEMP
 2170           CONTINUE
C
C      *      THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32        *
C      *      APPL.STAT. VOL 30, P. 105 (1981).                             *
C
                KK=II-MPLONE+K
                WORK(KK)=WORK(KK)+YS*X(I,K)*TEMP
 2190        CONTINUE
             WORK(II)=WORK(II)+YS**2*TEMP

C
C      *       OBSERVATION CONFINED BETWEEN TWO FINITE LIMITS               *
C

          ELSEIF((IPT-1) .GT. 0) THEN
             YN=DEXP(-0.5*YS**2)*C
             CALL RMILLS(YS,F,RLIMIT)
             YQ=YN/F
             YSU=YS+(Y2(I)-Y1(I))/XSIG
             CALL RMILLS(YSU,FU,RLIMIT)
             YNU=DEXP(-0.5*YSU**2)*C
             YQU=YNU/FU
             TINT=YQ-YQU
             A=(YN-YNU)/TINT
             B=(YNU*YSU-YN*YS)/TINT
             TEMP=A**2+B

             DO 2350 K=1,M

                DO 2330 J=1,K
                   JJ=JJ+1
                   WORK(JJ)=WORK(JJ)+X(I,J)*X(I,K)*TEMP
 2330           CONTINUE
                TEMP=(YS**2*YN-YSU**2*YNU)/TINT
C
C
C      *     THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32         *
C      *     APPL.STAT. VOL 30, P. 105 (1981).                              *
C
C
                KK=II-MPLONE+K
                WORK(KK)=WORK(KK)-(TEMP+A*B)*X(I,K)
 2350        CONTINUE

             TEMP=(YS**3*YN-YSU**3*YNU)/TINT
             WORK(II)=WORK(II)-TEMP+B**2
           ENDIF
 2500  CONTINUE
C
C      *                   INVERT THE MATRIX                                *
C
       CALL SYMINV(WORK,MPLONE,VCOV,W,NUL,LENWRK,LENWRK,LENW,IFAULT)

       IF((IFAULT.NE.0).OR.(NUL.NE.0)) THEN
          VCOV(2)=REAL(IFAULT)
          VCOV(1)=REAL(NUL)
          IFAULT=-6
          ICHECK=IJ
          RETURN
       ENDIF
C
C      *               RESTORE THE ITERATION COUNTER                         *
C
       IFAULT=IJ
C
C      *               MULTIPLY BY SIGMA-SQUARED                            *
C
       TEMP=XSIG**2
       DO 2580 I=1,II
          VCOV(I)=VCOV(I)*TEMP
 2580  CONTINUE

C
C      *                UNPACK THE MATRIX                                   *
C
       CALL UNPACK(VCOV,MPLONE,LENWRK)

       RETURN
       END

C
C      **********************************************************************
C      ********************* FUNCTION GAMMA  ********************************
C      **********************************************************************
C
       FUNCTION GAMMA(X)
C
C      *     THIS FUNCTION IS OBTAINED FROM PHILIP R. BEVINGTON, "DATA      *
C      *     REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", 1969, *
C      *     McGRAW HILL (NY:NY), PROGRAM 7-2 P. 126                        *
C      *     THIS COMPUTES THE GAMMA FUNCTION FOR INTEGERS AND HALF-INTEGERS*
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
C      *                INTEGERIZE ARGUMENT                                 *
C
       N = INT(X - 0.25)
       XN=REAL(N)
   13  IF((X-XN-0.75) .GT. 0.0) THEN 
C
C      *                ARGUMENT IS INTEGER                                 *
C
          GAMMA=FACTOR(N)
C
C      *                ARGUMENT IS HALF-INTEGER                            *
C
       ELSE
   31     PROD=1.77245385D00
          IF(N .LE. 0) THEN
             GAMMA=PROD
             GOTO 56
          ENDIF
          IF(N-10 .LE. 0) THEN
             DO 43 I=1,N
                FI=REAL(I)
                PROD=PROD*(FI-0.5)
   43        CONTINUE
             GAMMA=PROD
          ELSE
   51        SUM=0.0
             DO 54 I=11,N
                FI=I
                SUM=SUM+DLOG(FI-0.5)
   54        CONTINUE
   55        GAMMA=PROD*639383.8623D00*DEXP(SUM)
          ENDIF
   56  ENDIF
   60  RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE GEHAN  ******************************
C      **********************************************************************
C
       SUBROUTINE GEHAN(R1,R2,TEST,PROB,XY,ID1,ID2,NTOT)
C
C 
C      * THIS SUBROUTINE COMPUTES GEHAN'S GENERALIZED WILCOXON TEST         *
C      * STATISTIC.  THE COMPUTATIONAL FORM IS FROM E.T. LEE , STATISTICAL  *
C      * METHODS FOR SURVIVAL DATA ANALYSIS, 1980, LIFETIME LEARNING        *
C      * PUBLICATIONS, BELM0NT, CA. THE FORM USED IS THE MANTEL METHOD OF   *
C      * ASSIGNING A SCORE TO EACH OBSERVATION BASED ON ITS RELATIVE RANK,  *
C      * FROM EQUATION 5.5 AND EXAMPLE 5.1                                  *
C      *                                                                    *
C      *         SUBROUTINES                                                *
C      *                   STAT                                             *
C

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION R1(NTOT),R2(NTOT)
       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
C
C      *        COMPUTATION OF R1                                           *
C      *        STEP 1 AND 2 : RANK FROM LEFT TO RIGHT, OMITTING            *
C      *        RIGHT CENSORED VALUES. ASSIGN NEXT HIGHER RANK              *
C      *        TO RIGHT CENSORED VALUES                                    *
C
C
       IRANK=0
       DO 90 I=1,NCOMP
          IF(ID1(I).EQ.1) THEN
             R1(I)=IRANK+1
          ELSE
             IRANK=IRANK+1
             R1(I)=REAL(IRANK)
          ENDIF
  90   CONTINUE
C
C      *        STEP3 : REDUCE THE RANK OF TIED OBSERVATIONS TO             *
C      *        THE LOWEST RANK FOR THE VALUE                               *
C
       K1=NCOMP-1
       L1=1
  12   IF(XY(L1).EQ.XY(L1+1)) THEN
C
          JEMP=IABS(ID1(L1)-1)*IABS(ID1(L1+1)-1)
          IF(JEMP.NE.0) THEN
             R1(L1+1)=R1(L1)
             IF(L1.EQ.K1) GOTO 13
             L1=L1+1
             GOTO 12
          ENDIF
       ENDIF
       IF(L1.NE.K1) THEN
          L1=L1+1
          GOTO 12
       ENDIF
C
C      *             COMPUTATION OF R2                                      *
C      *             STEP 1 : RANK FROM RIGHT TO LEFT                       *
C
  13   DO 14 I=1,NCOMP
          R2(I)=REAL(NCOMP-I+1)
  14   CONTINUE
C
C      *          STEP2 : REDUCE THE RANK OF TIED OBSERVATIONS              *
C      *          TO THE LOWEST FOR THE VALUE                               *
C
       L1=NCOMP
  22   IF(XY(L1).EQ.XY(L1-1)) THEN
C
          JEMP=IABS(ID1(L1)-1)*IABS(ID1(L1-1)-1)
          IF(JEMP.NE.0) THEN
             R2(L1-1)=R2(L1)
             IF(L1.EQ.2) GOTO 23
             L1=L1-1
             GOTO 22
          ENDIF
       ENDIF
       IF(L1.NE.2) THEN
          L1=L1-1
          GOTO 22
       ENDIF

  23   IF(NCEN.NE.0) THEN
C
C      *         STEP 3 : REDUCE THE RANK OF RIGHT CENSORED                 *
C      *         OBSERVATION TO UNITY                                       *
C
          DO 24 I=1,NCOMP
             IF(ID1(I).NE.0)  R2(I)=1.0
  24      CONTINUE
       ENDIF
C
C      *               COMPUTE FINAL SCORES - R1(I)                         *
       DO 25 I=1,NCOMP
          R1(I)=R1(I)-R2(I)
  25   CONTINUE

       CALL STATS(R1,TEST,XY,ID1,ID2,NTOT)
       PROB=1.0-AGAUSS(TEST)

       RETURN
       END



C      **********************************************************************
C      ******************** SUBROUTINE GRDPROB  *****************************
C      **********************************************************************
C
       SUBROUTINE GRDPRB(NTOT,MX,MY,SUM,ISKIP,ICENS,DELX,
     +                   DELY,XORG,YORG,TOL,MAX,MM,M1,M2,M3,M4,M5,
     +                   M6,M7,M8,X,Y,NP,XB,YB,F,FC,N,N1,N2,N3,
     +                   N4,N5,N6,N7,N8,IWRK1,IWRK2,WRK1,WRK2,
     +                   SWRK1,DWRK1,IB,MVAR)
C
C
C      *                                                                    *
C      *     THIS SUBPROGRAM COMPUTES THE PROBABILITY OF BIN(I,J)           *
C      *     IN WHICH DETECTED POINTS EXIST. ONLY BINS WHICH HAVE           *
C      *     DETECTED POINTS CAN HAVE NON-ZERO PROBABILITY.                 *
C      *                                                                    *
C      *     INPUT                                                          *
C      *              X(I)  : INDEPENDENT VARIABLE                          *
C      *              Y(I)  : DEPENDENT VARIABLE                            *
C      *             NP(I)  : INDICATOR OF CENSORED STATUS                  *
C      *              NT    : TOTAL NUMBER OF DATA                          *
C      *              MX    : BIN NUMBER OF X                               *
C      *              MY    : BIN NUMBER OF Y                               *
C      *             ISKIP  : INDICATOR OF BINNING                          *
C      *             ICENS  : CENSORING STATUS OF DATA SET                  *
C      *             TOL    : TOLERANCE FOR COMPUTAION F(I,J)               *
C      *             MAX    : MAXIMUM ITERATION                             *
C      *        IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED  :             *
C      *             DELX   : BIN SIZE OF X AXIS                            *
C      *             DELY   : BIN SIZE OF Y AXIS                            *
C      *             XORG   : ORIGIN OF X                                   *
C      *             YORG   : ORIGIN OF Y                                   *
C      *      WORK                                                          *
C      *             FC(I,J): COPY OF F(I,J) (TO CHECK THE CONVERGENCE)     *
C      *             N1(I,J): NUMBER OF Y LOWER LIMITS IN THE BIN (I,J)     *
C      *             N2(I,J): NUMBER OF X LOWER LIMITS IN THE BIN (I,J)     *
C      *             N3(I,J): NUMBER OF DOUBLE LOWER LIMITS IN THE          *
C      *                      BIN(I,J)                                      *
C      *             N4(I,J): NUMBER OF Y LOWER X UPPER LIMTS IN THE        *
C      *                      BIN(I,J)                                      *
C      *             N5(I,J): NUMBER OF Y UPPER LIMITS IN THE BIN (I,J)     *
C      *             N6(I,J): NUMBER OF X UPPER LIMITS IN THE BIN (I,J)     *
C      *             N7(I,J): NUMBER OF DOUBLE UPPER LIMITS IN THE          *
C      *                      BIN(I,J)                                      *
C      *             N8(I,J): NUMBER OF Y UPPER X LOWER LIMTS IN THE        *
C      *                      BIN(I,J)                                      *
C      *              SUM1  : WEIGHT ON BIN(I,J) FROM Y LOWER LIMIT         *
C      *              SUM2  : WEIGHT ON BIN(I,J) FROM X LOWER LIMITS        *
C      *              SUM3  : WEIGHT ON BIN(I,J) FROM DOUBLE LOWER          *
C      *                      LIMITS                                        *
C      *              SUM4  : WEIGHT ON BIN(I,J) FROM Y LOWER, X UPPER      *
C      *                      LIMITS                                        *
C      *              SUM5  : WEIGHT ON BIN(I,J) FROM Y UPPER LIMITS        *
C      *              SUM6  : WEIGHT ON BIN(I,J) FROM X UPPER LIMITS        *
C      *              SUM7  : WEIGHT ON BIN(I,J) FROM DOUBLE UPPER          *
C      *                      LIMITS                                        *
C      *              SUM8  : WEIGHT ON BIN(I,J) FROM Y UPPER, X LOWER      *
C      *                      LIMITS                                        *
C      *              ITE   : NUMBER OF ITERATIONS                          *
C      *              DEL   : TOLERANCE [SUM (F(I,J)-FC(I,J))**2]           *
C      *     OUTPUT                                                         *
C      *              F(I,J): NUMBER OF DATA POINTS IN BIN(I,J)             *
C      *              XB(I) : POSITION OF  X BIN                            *
C      *              YB(I) : POSITION OF  Y BIN                            *
C      *                                                                    *
C      *     SUBROUTINE : BIN                                               *
C      *                                                                    *
C
       IMPLICIT REAL*8(A-H,O-Z), INTEGER (I-N)
       DIMENSION FC(IB,IB)
       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB)
       DIMENSION F(IB,IB),N(IB,IB),N1(IB,IB),N2(IB,IB),N3(IB,IB)
       DIMENSION N4(IB,IB),N5(IB,IB),N6(IB,IB),N7(IB,IB),N8(IB,IB)

       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
       DIMENSION DWRK1(MVAR,NTOT),SWRK1(MVAR)

       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
C
C      *              CALL SUBROUTINE BIN                                   *
C
       CALL BIN(NTOT,MX,MY,ISKIP,ICENS,DELX,DELY,XORG,YORG,MM,M1,M2,
     +          M3,M4,M5,M6,M7,M8,IWRK1,IWRK2,WRK1,WRK2,DWRK1,SWRK1,
     +          X,Y,NP,XB,YB,F,N,N1,N2,N3,N4,N5,N6,N7,N8,IB,MVAR)

C
C      *             INITIAL SETTING OF FC(I,J)                             *
C
       DO 300 I=1,MX
          DO 200 J=1,MY
             FC(I,J)=0.0
  200     CONTINUE
  300  CONTINUE
C
C
C      *             START ITERATIONS TO FIND F(I,J) LOOP 500               *
C
       SNT=NTOT
       ITE=1

  500  DEL=0.0
       DO 1600 I=1,MX
          DO 1500 J=1,MY
             IF(F(I,J).NE.0.0) THEN
                SUM1=0
                SUM2=0
                SUM3=0
                SUM4=0
                SUM5=0
                SUM6=0
                SUM7=0
                SUM8=0
C
C      *      COMPUTE THE INFLUENCE OF CENSORED DATA POINTS ON THE DETECTED *
C      *      POINT AT I,J.                                                 *
C
C
C   
C      *              Y LOWER LIMITS                                        *
C
                IF(NC1.NE.0) THEN
                   DO 600 L=1,J
                      IF(N1(I,L).NE.0) THEN
                         SUMF1=0.0
                         DO 550 L1=L,MY
                            SUMF1=SUMF1+F(I,L1)
  550                    CONTINUE

               SUM1=SUM1+(FLOAT(N1(I,L))/SNT)*(F(I,J)/SUMF1)

                      ENDIF
  600              CONTINUE
                ENDIF
C
C
C      *             X LOWER LIMITS                                         *
C
                IF(NC2.NE.0) THEN
                   DO 700 K=1,I
                      IF(N2(K,J).NE.0) THEN
                         SUMF2=0.0
                         DO 650 K1=K,MX
                            SUMF2=SUMF2+F(K1,J)
  650                    CONTINUE

                SUM2=SUM2+(FLOAT(N2(K,J))/SNT)*(F(I,J)/SUMF2)

                      ENDIF
  700              CONTINUE
                ENDIF
C
C
C      *            DOUBLE LOWER LIMITS                                     *
C
                IF(NC3.NE.0) THEN
                   DO 800 K=1,I
                      DO 790 L=1,J
                         IF(N3(K,L).NE.0) THEN
                            SUMF3=0.0
                            DO 780 K1=K,MX
                               DO 770 L1=L,MY
                                  SUMF3=SUMF3+F(K1,L1)
  770                          CONTINUE
  780                       CONTINUE

                SUM3=SUM3+(FLOAT(N3(K,L))/SNT)*(F(I,J)/SUMF3)

                         ENDIF
  790                 CONTINUE
  800              CONTINUE
                ENDIF
C
C
C      *            Y LOWER, X UPPER LIMITS                                 *
C
                IF(NC4.NE.0) THEN
                   DO 900 K=1,MX
                      KK=MX-K+1
                      IF(KK.LT.I) GOTO 910
                      DO 890 L=1,J
                         IF(N4(KK,L).NE.0) THEN
                            SUMF4=0.0
                            DO 880 K1=1,KK
                               DO 870 L1=L,MY
                                  SUMF4=SUMF4+F(K1,L1)
  870                          CONTINUE
  880                       CONTINUE

                SUM4=SUM4+(FLOAT(N4(K,L))/SNT)*(F(I,J)/SUMF4)

                         ENDIF
  890                 CONTINUE
  900              CONTINUE
                ENDIF
C
C
C      *            Y UPPER LIMITS                                          *
C
  910           IF(NC5.NE.0) THEN
                   DO 1000 L=1,MY
                      LL=MY-L+1
                      IF(LL.LT.J) GOTO 1010
                      IF(N5(I,LL).NE.0) THEN
                         SUMF5=0.0
                         DO 950 L1=1,LL
                            SUMF5=SUMF5+F(I,L1)
  950                    CONTINUE

                SUM5=SUM5+(FLOAT(N5(I,LL))/SNT)*(F(I,J)/SUMF5)

                      ENDIF
 1000              CONTINUE
                ENDIF
C
C
C      *               X UPPER LIMITS                                       *
C
 1010           IF(NC6.NE.0) THEN
                   DO 1100 K=1,MX
                      KK=MX-K+1
                      IF(KK.LT.I) GOTO 1110
                      IF(N6(KK,J).NE.0) THEN
                         SUMF6=0.0
                         DO 1050 K1=1,KK
                            SUMF6=SUMF6+F(K1,J)
 1050                    CONTINUE

                SUM6=SUM6+(FLOAT(N6(KK,J))/SNT)*(F(I,J)/SUMF6)

                      ENDIF
 1100              CONTINUE
                ENDIF
C
C
C      *            DOUBLE UPPER LIMITS                                     *
C
 1110           IF(NC7.NE.0) THEN
                   DO 1200 K=1,MX
                      KK=MX-K+1
                      IF(KK.LT.I) GOTO 1210
                      DO 1190 L=1,MY
                         LL=MY-L+1
                         IF(LL.LT.J) GOTO 1200
                         IF(N7(KK,LL).NE.0) THEN
                            SUMF7=0.0
                            DO 1180 K1=1,KK
                               DO 1170 L1=1,LL
                                  SUMF7=SUMF7+F(K1,L1)
 1170                          CONTINUE
 1180                       CONTINUE

                SUM7=SUM7+(FLOAT(N7(KK,LL))/SNT)*(F(I,J)/SUMF7)

                         ENDIF
 1190                 CONTINUE
 1200              CONTINUE
                ENDIF
C
C
C      *               Y UPPER, X LOWER LIMITS                              *
C
 1210           IF(NC8.NE.0) THEN
                   DO 1300 K=1,I
                      DO 1290 L=1,MY
                         LL=MY-L+1
                         IF(LL.LT.J) GOTO 1300
                         IF(N8(K,LL).NE.0) THEN
                            SUMF8=0.0
                            DO 1280 K1=K,MX
                               DO 1270 L1=1,LL
                                  SUMF8=SUMF8+F(K1,L1)
 1270                          CONTINUE
 1280                       CONTINUE

                SUM8=SUM8+(FLOAT(N8(KK,LL))/SNT)*(F(I,J)/SUMF8)

                         ENDIF
 1290                 CONTINUE
 1300              CONTINUE
                ENDIF
C
C      *            COMPUTE A NEW ESTIMATE OF F(I,J).                       *
C
 1400           SUM=SUM1+SUM2+SUM3+SUM4+SUM5+SUM6+SUM7+SUM8
                F(I,J)=FLOAT(N(I,J))/SNT+SUM
                DEL=DEL+(F(I,J)-FC(I,J))**2
                FC(I,J)=F(I,J)
             ENDIF
 1500     CONTINUE
 1600  CONTINUE
C
C      *             CHECK CONVERGENCE                                      *
C
       ITE=ITE+1

       IF(((DEL).GT.TOL).AND.(ITE.LE.MAX)) GOTO 500
       RETURN
       END



C
C      **********************************************************************
C      ******************** SUBROUTINE KMADJ  *******************************
C      **********************************************************************
C
       SUBROUTINE KMADJ(ZU,ZC,NTOT,IU,IC,S,ISIGN,NTEMP,F,V)
C
C      *       THIS SUBROUTINE RESTORES THE DATA AND THE PRODUCT-LIMIT      *
C      *       ESTIMATOR TO UPPER-LIMITS FORM, IF THE DATA WAS IN THAT FORM *
C      *       INITIALLY.  TIES AT CENSORED POINTS ARE ALSO REMOVED TO      *
C      *       MAKE THE PRINTOUT CONSISTENT.                                *
C      *                                                                    *
C      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
C      *                ZC(I)  :  CENSORED DATA POINTS                      *
C      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
C      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
C      *                 IC    :  NUMBER OF CENSORED DATA POINTS            *
C      *                 S(L)  :  PL ESTIMATOR                              *
C      *       OUTPUT  NTEMP   :  VALUE OF NTOT AFTER ADJUSTMENT FOR TIES   *
C      *       OTHER    F      :  PROBABILITY MASS ASSIGNED TO EACH         *
C      *                             UNCENSORED POINT(=JUMP IN S AT THE     *
C      *                                                  POINT)            *
C      *                                                                    *

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION ZU(NTOT),ZC(NTOT),S(NTOT),F(NTOT),V(NTOT)

C
C      *  FOR LEFT-CENSORED DATASETS (I.E. UPPER LIMITS),                  *
C      *  CALCULATE JUMP IN SURVIVOR FUNCTION AT UNCENSORED POINTS         *
C
       IF(ISIGN.LT.0) THEN
            F(1) = 1.0 - S(1)
            DO 120 I = 2, IU
               F(I) = S(I-1)-S(I)
120         CONTINUE
C
C      *  REVERSE SEQUENCE OF UNCENSORED POINTS, JUMPS AND ERRORS         *
C
            J = IU/2
            DO 150 I =1, J

               Z = ZU(I)*(-1.0)
               ZU(I) = ZU(IU-I+1)*(-1.0)
               ZU(IU-I+1) = Z

               FTEMP = F(I)
               F(I) = F(IU-I+1)
               F(IU-I+1) = FTEMP 

               VTEMP = V(I)
C               V(I) = V(IU-I+1)
C               V(IU-I+1) = VTEMP
               V(I) = V(IU-I)
               V(IU-I) = VTEMP

150         CONTINUE

            IF(2*J.NE.IU) THEN
               ZU(J+1) = ZU(J+1)*(-1.0)
            ENDIF

C
C      *  REVERSE SEQUENCE OF CENSORED POINTS                              *
C
            J = IC/2
            DO 155 I = 1, J
               Z = ZC(I) * (-1.0)
               ZC(I) = ZC(IC-I+1)*(-1.0)
               ZC(IC-I+1) = Z
155         CONTINUE
 
            IF(2*J.NE.IC) THEN
               ZC(J+1) = ZC(J+1)*(-1.0)
            ENDIF

C
C      *  COMPUTE SURVIVOR FUNCTION FOR REVERSED DATA                     *
C
            DO 170 I = 1, IU
               S(I) = 1
               DO 160 J = 1, I
                  S(I) = S(I) - F(J)
160            CONTINUE
170         CONTINUE
         ENDIF   

C      *   CORRECTS FOR TIES AT THE UNCENSORED POINTS                      *
C      *   NOTICE THAT IF THERE ARE TIES AT THESE POINTS, THEN BOTH        *
C      *   IU AND NTEMP ARE REDUCED.                                       *

       NTEMP = NTOT
       K = 1
190    IF(ZU(K).EQ.ZU(K+1)) THEN
          DO 195 I = K, IU-1
             ZU(I)=ZU(I+1)
             S(I)=S(I+1)
             V(I) = V(I+1)               
195       CONTINUE
          IU = IU -1
          NTEMP = NTEMP - 1
       ELSE
          K  = K + 1
       ENDIF
       IF(K.LT.IU) GOTO 190

       RETURN
       END


C
C      **********************************************************************
C      ******************** SUBROUTINE KMDIF ********************************
C      **********************************************************************
C
       SUBROUTINE KMDIF(S,ZU,BS,BL,DIFF,F,NTOT,START,BINSIZ,LSTEP,
     +                   OUT,IBIN,IU)

C
C      *       THIS SUBROUTINE COMPUTES AND PRINTS THE DIFFERENTIAL KM      *
C      *       ESTIMATOR BASED ON WARDLE AND KNAPP, 'THE STATISTICAL        *
C      *       DISTRIBUTION OF THE NEUTRAL-HYDROGEN CONTENT OF S0 GALAXIES',*
C      *       ASTRN. J., 91:23 1986.                                       *
C      *                                                                    *
C      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
C      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
C      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
C      *                 S(L)  :  PL ESTIMATOR                              *
C      *                 OUT   :  OUTPUT FILE NAME. IF IT IS BLANK          *
C      *                          THE RESULTS WILL BE SHOWN ON THE          *
C      *                          TERMINAL.                                 *
C      *                START  :  STARTING VALUE OF THE FIRST BIN           *
C      *                BINSIZ :  WIDTH OF THE BIN                          *
C      *                LSTEP  :  NUMBER OF BINS                            *
C      *                IBIN   :  DIMENSION                                 *
C      *              ICHANGE  :  INDICATES IF THE LAST POINT (OR THE       *
C      *                            FIRST POINT FOR UPPER LIMITS DATA)      *
C      *                            HAS BEEN CHANGED TO A DETECTION.        *
C      *                                                                    *
C      *      OTHERS                                                        *
C      *               BS(J)   :  STARTING VALUE FOR THE BIN J              *
C      *               BL(J)   :  ENDING VALUE FOR THE BIN J                *
C      *               DIFF(J) :  DIFFERENTIAL KM ESTIMATOR AT BIN J        *
C      *               F(I)    :  MASS OF THE I TH DATA POINT               *
C      *                                                                    *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*9 OUT,CHECK

       DIMENSION ZU(NTOT),S(NTOT),F(NTOT)
       DIMENSION BS(IBIN),BL(IBIN),DIFF(IBIN)
       CHECK='         '
C
C       *    FIRST, COMPUTE THE MASS OF EACH POINT                          *
C
        F(1) = 1.0 -S(1)
       
        DO 610 I = 2, IU
           F(I) = DABS(S(I) - S(I-1))
  610   CONTINUE

C
C       *  SET THE BIN BOUNDARIES.                                          *
C
        DO 620 J = 1, LSTEP
           BS(J) = START + BINSIZ*(J-1)
           BL(J) = START + BINSIZ*J
  620   CONTINUE

        I = 1
        J = 0

  630   J = J + 1
        DIFF(J) = 0.0

C
C      *       CHECK WHETHER THE I-TH POINT IS IN THE BIN                  *
C
  640  IF(J .LE. LSTEP) THEN
         IF(ZU(I) .LT. BS(J)) THEN
            IF(I .GE. IU) THEN
               GOTO 630
            ENDIF
            I = I + 1
            GOTO 640
         ENDIF

C      *       COMPUTE THE DIFFERENTIAL KM                                 *
C
         IF(ZU(I) .GE. BL(J)) GOTO 630
         DIFF(J) = DIFF(J) + F(I)
   
         IF(I .LT. IU) THEN
            I = I + 1
            GOTO 640
         ENDIF
         GOTO 630
       ENDIF
C
C      *         START PRINTING THE RESULTS                                  *
C
          IF(OUT.EQ.CHECK) THEN
             PRINT *
             PRINT *,'           DIFFERENTIAL KM ESTIMATOR'
             PRINT 660
             PRINT *
          ELSE
             WRITE(60, 658)
             WRITE(60,659)
             WRITE(60, 660)
             WRITE(60,658)
          ENDIF
  658     FORMAT('  ')
  659     FORMAT(5X,'   DIFFERENTIAL KM ESTIMATOR')
  660     FORMAT(5X,'   BIN CENTER          D')

C
C      * MULTIPLY DIFF(I) BY THE TOTAL NUMBER OF POINTS TO GET THE NUMBER    *
C      * OF POINTS IN EACH BIN.                                              *
C
          SUM = 0.0
          DO 690 I = 1, LSTEP
             DIFF(I) =DIFF(I)*NTOT
             CENTER = 0.5*(BS(I) + BL(I))
             IF(OUT .EQ. CHECK) THEN
                PRINT 680, CENTER, DIFF(I)
             ELSE
                WRITE(60,680) CENTER, DIFF(I)
             ENDIF
  680        FORMAT(2F15.3)
             SUM = SUM + DIFF(I)
  690     CONTINUE

          IF(OUT .EQ. CHECK) THEN
             PRINT 700, SUM
             PRINT 658
          ELSE
             WRITE(60,700) SUM
             WRITE(60,658)
             WRITE(60,701)
 701         FORMAT(' (D GIVES THE ESTIMATED DATA POINTS IN EACH BIN)')
          ENDIF

 700      FORMAT(23X,'-------',/10X,'SUM =',F15.3)

       RETURN
       END


C
C      **********************************************************************
C      ********************* SUBROUTINE KMESTM ******************************
C      **********************************************************************
C
       SUBROUTINE KMESTM(IND,X,NTOT,J,IPRINT,TITLE,NAME,OUTPUT,IBIN,
     +                   ISKIP,KDIFF,START,BINSIZ,LSTEP,FILE,
     +                   ZXU,ZXC,SX,VX,Z1,ITEMP,INDEX,ATIE,RISK,
     +                   BWRK1,BWRK2,BWRK3,IWRK1,SWRK1,
     +                   WRK1,MVAR)
C
C      *       THIS SUBROUTINE COMPUTES THE PL ESTIMATOR, MEAN AND ITS      *
C      *       ERROR FOR THE X VARIABLE.                                    *
C      *                                                                    *
C      *       INPUT IND(J,I): INDICATOR OF CENSORING                       *
C      *               X(J,I): DATA                                         *
C      *                NTOT : NO. OF DATA POINTS                           *
C      *                 J   : J-TH VARIABLES                               *
C      *               IPRINT: IF 0, PRINT OUT RESULTS ONLY                 *
C      *                       IF 1, PRINT OUT ALL                          *
C      *              PROBLEM: TITLE OF THE PROBLEM                         *
C      *               NAME  : NAME OF THE SUB-DATA SET                     *
C      *               OUTPUT: NAME OF OUTPUT FILE                          *
C      *                       IF IT IS BLANK, SHOW THE RESULT ON THE       *
C      *                       TERMINAL.                                    *
C      *               ISKIP : IF THE SUBROUTINE IS CALLED BY TWO SAMPLE    *
C      *                       TESTS, ISKIP=1 AND SKIP A FEW LINES .        *
C      *               KDIFF : IF KDIFF = 1, PRINT OUT DIFFERENTIAL KM      *
C      *               START : STARTING POINT OF BINING                     *
C      *               BINSIZ: WIDTH OF THE BIN                             *
C      *               LSTEP : NUMBER OF BINS                               *
C      *              ATIE(I): NUMBER OF TIED DATA POINTS AT ITH DATA VALUE *
C      *              RISK(I): RISK SET FOR THE ITH DATA VALUE              *
C      *               MVAR  : DIMENSION SIZE                               *
C      *               IBIN  : DIMENSION SIZE                               *
C      *                                                                    *
C      *       WORK      ZXU : DATA ARRAY CONTAINING THE UNCENSORED POINTS  *
C      *                 ZXC : DATA ARRAY CONTAINING THE CENSORED POINTS    *
C      *                 IXU : NO. OF UNCENSORED DATA POINTS                *
C      *                 IXC : NO. OF CENSORED DATA POINTS                  *
C      *              ICHANGE: IF THE LAST VALUE IS CENSORED, THE VALUE     *
C      *                       NEEDS TO BE CHANGED TO A DETECTION.          *
C      *                       THIS INDICATOR IS SET TO BE -1,IF THE LAST   *
C      *                       VALUE IS CHANGED.                            *
C      *                                                                    *
C      *       OUTPUT    SX  : PL ESTIMATOR                                 *
C      *                 VX  : ERROR OF PL ESTIMATOR                        *
C      *                SMEAN: MEAN                                         *
C      *                ERROR: ERROR OF THE MEAN                            *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *                 SORT1, XVAR, PLESTM, KMADJ, KMPRNT                 *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*80 TITLE
       CHARACTER*9 NAME,OUTPUT,FILE

       DIMENSION IND(MVAR, NTOT),X(MVAR, NTOT),ZXU(NTOT),ZXC(NTOT)
       DIMENSION SX(NTOT),VX(NTOT),Z1(MVAR,NTOT),ITEMP(NTOT)
       DIMENSION INDEX(NTOT),ATIE(NTOT),RISK(NTOT)
       DIMENSION BWRK1(IBIN),BWRK2(IBIN),BWRK3(IBIN)
       DIMENSION IWRK1(NTOT),SWRK1(MVAR),WRK1(NTOT)

C
C      *       DISTINGUISH UNCENSORED AND CENSORED DATA AND SORT THEM IN    *
C      *       ASCENDING ORDER. THEN CALL PL ESTIMATOR SUBROUTINE           *
C      *       "PL".                                                        *
C
       IF(OUTPUT .EQ. '         ') THEN

          PRINT *
          PRINT *
          PRINT 40
          PRINT 44
          IF(ISKIP.NE.1) PRINT 46,TITLE   
          PRINT 44
          PRINT 47,FILE
          PRINT 44
          PRINT 48,NAME
          PRINT 44
       ELSE
          WRITE(60,44)
          WRITE(60,44)
          WRITE(60,40)
          WRITE(60,44)
          IF(ISKIP.NE.1) WRITE(60,46) TITLE   
          WRITE(60,44)
          WRITE(60,47) FILE
          WRITE(60,44)
          WRITE(60,48) NAME
          WRITE(60,44)
       ENDIF

   40  FORMAT(8X,'KAPLAN-MEIER ESTIMATOR')
   44  FORMAT('    ')
   46  FORMAT(8X,'TITLE : ',A80)
   47  FORMAT(8X,'DATA FILE : ',A9)  
   48  FORMAT(8X,'VARIABLE : ',A9)

C
C      *    XVAR DISTINGUISHES DETECTED POINTS AND CENSORED POINTS           *
C
       CALL XVAR(IND,X,J,NTOT,ISIGN,ZXU,ZXC,IXU,IXC,IWRK1,
     +           ATIE,RISK,WRK1,Z1,SWRK1,LTOT,MVAR,INDEX)
C
       IF(OUTPUT .EQ. '         ') THEN

          PRINT *
          IF(ISIGN.EQ.1) PRINT 56,NTOT,IXC
          IF(ISIGN.EQ.-1) PRINT 57,NTOT,IXC
          PRINT *
       ELSE
          WRITE(60,44)
          IF(ISIGN.EQ.1) WRITE(60,56) NTOT,IXC
          IF(ISIGN.EQ.-1) WRITE(60,57) NTOT,IXC
          WRITE(60,44)
       ENDIF

   56  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF LOWER LIMITS :',I4)
   57  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF UPPER LIMITS :',I4)
C
C      *  SET A FEW DUMMY ARRAYS TO USE THE SORTING PRGRAM                  *
C
       ISTEMP = ISIGN
   59  DO 60 I=1,NTOT
          ITEMP(I)=0
          Z1(1,I)=1.0
   60  CONTINUE
C
       CALL SORT1(ITEMP,Z1,ZXU,IXU,1,INDEX,SWRK1,MVAR)
C
       CALL SORT1(ITEMP,Z1,ZXC,IXC,1,INDEX,SWRK1,MVAR)
C
C      *  CALL SUBROUTINE "PLESTM" TO COMPUTE KM ESTIMATOR                  *
C
       CALL PLESTM(ZXU,ZXC,IXU,IXC,SX,VX,NTOT,SMEAN,ERROR,ICHANGE,
     +             NCHANGE,IWRK1)
C
C      *        IF THE DATA CONTAINS UPPER LIMITS, CHANGE THE               *
C      *        SIGN OF THE MEAN.                                           *
C
       ISIGN = ISTEMP
       SMEAN=ISIGN*SMEAN
C
C      * SUBROUTINE KMADJ IS CALLED TO ADJUST THE PRODUCT-LIMIT ESTIMATOR   *
C      * BACK TO THE ORIGIONAL CENSORING PATTERN AND TO REMOVE TIES.        *
C
       CALL KMADJ(ZXU,ZXC,NTOT,IXU,IXC,SX,ISIGN,NTEMP,WRK1,VX)

C
C      * PRINT PL ESTIMATOR, PERCENTILES, MEAN, AND ERROR                   *
C
C

       CALL KMPRNT(ZXU,ZXC,NTOT,NTEMP,IXU,IXC,SX,VX,ISIGN,OUTPUT,
     +             ICHANGE,SMEAN,ERROR,IPRINT)


C      *  SUBROUTINE KMDIF IS CALLED IF THE USER HAS REQUESTED A            *
C      *  DIFFERENTIAL KM ESTIMATOR.                                        *

       IF(KDIFF .EQ. 1) THEN

          CALL KMDIF(SX,ZXU,BWRK1,BWRK2,BWRK3,WRK1,NTOT,START,
     +               BINSIZ,LSTEP,OUTPUT,IBIN,IXU)

       ENDIF

       PRINT *   

C
C      *   IF THE LAST VALUE WAS CHANGED FROM AN UPPER LIMIT TO A          *
C      *   DETECTION, CHANGE THE NUMBER BACK TO ITS ORIGINAL VALUE.        *
C
       IF(ICHANGE.EQ.-1) THEN
          IXU=IXU-NCHANGE
          IXC=IXC+NCHANGE
       ENDIF
       
       RETURN
       END



C
C      **********************************************************************
C      ******************** SUBROUTINE KMPRNT *******************************
C      **********************************************************************
C
       SUBROUTINE KMPRNT(ZU,ZC,NTOT,NTEMP,IU,IC,S,V,ISIGN,OUT,ICHANGE,
     +                   SMEAN,ERROR,IPRINT)

C
C      *       THIS SUBROUTINE PRINTS KM ESTIMATORS, THEIR                  *
C      *       ERROR, AND 75, 50, AND 25 PERCENTILES. ADOPTED FROM          *
C      *       ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL DATA         *
C      *       ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA) *
C      *                                                                    *
C      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
C      *                ZC(I)  :  CENSORED DATA POINTS                      *
C      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
C      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
C      *                 IC    :  NUMBER OF CENSORED DATA POINTS            *
C      *                 S(L)  :  KM ESTIMATOR                              *
C      *                 V     :  ERROR OF KM ESTIMATOR                     *
C      *                ISIGN  :  INDICATOR OF LOWER/UPPER LIMIT            *
C      *                            IF 1, LOWER LIMIT                       *
C      *                           IF -1, UPPER LIMIT                       *
C      *                 D(L)  :  NUMBER OF TIED DATA POINTS AT THE VALUE   *
C      *                 R(L)  :  RISK SET                                  *
C      *                 OUT   :  OUTPUT FILE NAME. IF IT IS BLANK          *
C      *                          THE RESULTS WILL BE SHOWN ON THE          *
C      *                          TERMINAL.                                 *
C      *                KDIFF  :  PRINTOUT INDICATOR FOR THE DIFFERENTIAL   *
C      *                          KM ESTIMATOR (IF 1, PRINT )               *
C      *                START  :  STARTING VALUE OF THE FIRST BIN           *
C      *                BINSIZ :  WIDTH OF THE BIN                          *
C      *                LSTEP  :  NUMBER OF BINS                            *
C      *                IBIN   :  DIMENSION                                 *
C      *              ICHANGE  :  INDICATES IF THE LAST POINT (OR THE       *
C      *                            FIRST POINT FOR UPPER LIMITS DATA)      *
C      *                            HAS BEEN CHANGED TO A DETECTION.        *
C      *               IPRINT  :  INDICATES WHETHER TO PRINT OUT THE        *
C      *                            FULL KM ESTIMATE OR ONLY THE MEAN       *
C      *                            AND PERCENTILES                         * 
C      *                                                                    *
C      *      OTHERS                                                        *
C      *               BS(J)   :  STARTING VALUE FOR THE BIN J              *
C      *               BL(J)   :  ENDING VALUE FOR THE BIN J                *
C      *               DIFF(J) :  DIFFERENTIAL KM ESTIMATOR AT BIN J        *
C      *               ERR(J)  :  ERROR FOR THE BIN J                       *
C      *               F(I)    :  MASS OF THE I TH DATA POINT               *
C      *               SMEAN   :  MEAN                                      *
C      *               ERROR   :  ERROR OF THE MEAN                         *
C      *                                                                    *
C      *      SUBROUTINE:  SUMRY                                            *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*9 OUT,CHECK

       DIMENSION ZU(NTOT),ZC(NTOT),S(NTOT),V(NTOT),FINT(3)
       CHECK='         '
C
       IF(IPRINT .EQ. 1) THEN
          IF(OUT.EQ.CHECK) THEN
             PRINT 100
             PRINT 110
          ELSE 
             WRITE(60,100)
             WRITE(60,110)
          ENDIF

  100     FORMAT('    ')
  110     FORMAT(10X,'VARIABLE RANGE',6X,'KM ESTIMATOR',3X,'ERROR')

C      *  STARTS TO PRINT OUT THE RESULTS                                   *

          IF(OUT .EQ. CHECK) THEN
             PRINT 540,0.00,ZU(1),1.00
          ELSE
             WRITE(60,540) 0.00,ZU(1),1.00
          ENDIF

  540     FORMAT('FROM',F9.3,'   TO',F9.3,F12.3)
               
          KT=1
          KC=1
          KU=1

  200     IF(KT.LE.NTEMP) THEN
             IF(KC.GT.IC) GOTO 300
             IF(KU.GT.IU) GOTO 250
             IF(ZU(KU).LT.ZC(KC)) GOTO 300
C
C
  250        IF(OUT.EQ.CHECK) THEN
                PRINT 555,ZC(KC)
             ELSE
                WRITE(60,555) ZC(KC)
             ENDIF
             KC=KC+1
             GOTO 400

  300        IF(KU .LT. IU) THEN
                IF(OUT.EQ.CHECK) THEN
                   PRINT 550,ZU(KU),ZU(KU+1),S(KU),V(KU)
                ELSE
                   WRITE(60,550) ZU(KU),ZU(KU+1),S(KU),V(KU)
                ENDIF
             ELSE
                IF(OUT.EQ.CHECK) THEN
                   PRINT 560,ZU(KU), S(KU), V(KU)
                ELSE
                   WRITE(60,560) ZU(KU),S(KU),V(KU)
                ENDIF
             ENDIF
             KU=KU+1
  400        KT=KT+1
             GOTO 200

  550        FORMAT('FROM',F9.3,'   TO',F9.3,2F12.3)
  555        FORMAT(F13.3,' C ')
  560        FORMAT('FROM',F9.3,'   ONWARDS',4x,2F12.3)
          ENDIF
       ENDIF
      
C      *   PRINTS OUT A WARNING FLAG IF A CENSORED POINT HAS BEEN 
C      *   CHANGED TO A DETECTION
       IF(ICHANGE .EQ. -1) THEN
           IF(ISIGN .GT. 0) THEN 
               IF(OUT .EQ. CHECK) THEN
                   PRINT 565
                   PRINT 566
               ELSE
                   WRITE(60,565)
                   WRITE(60,566)
               ENDIF
           ELSE
               IF(OUT .EQ. CHECK) THEN
                    PRINT 570
                    PRINT 566
               ELSE
                    WRITE(60,570)
                    WRITE(60,566)
               ENDIF
           ENDIF
       ENDIF

  565     FORMAT(10X,
     +    'WARNING:  THE LAST POINT WAS CHANGED TO A DETECTION ')
  566     FORMAT(20X,'FOR THE K-M COMPUTATION')
  570     FORMAT(10X,
     +    'WARNING:  THE FIRST POINT WAS CHANGED TO A DETECTION')
  
C
C      *       COMPUTE 75, 50, AND 25 PERCENTILES AND PRINT THEM.           *
C
       IF(IU .LE. 3) THEN
          IF(OUT .EQ. CHECK) THEN
             PRINT 100
             PRINT 750
             PRINT 755
             PRINT 100
          ELSE 
             WRITE(60,100)
             WRITE(60,750)
             WRITE(60,755)
             WRITE(60,100)
          ENDIF
          GOTO 900
       ELSE
          CALL SUMRY(ZU,IU,S,NTOT,FINT)
       ENDIF

 750   FORMAT(/,6X,'SINCE THERE ARE LESS THAN 4 UNCENSORED POINTS,')
 755   FORMAT(6X,'NO PERCENTILES WERE COMPUTED.')

       IF(OUT.NE.CHECK) GOTO 760
       PRINT 100
       PRINT 780
       PRINT 800
       PRINT 850,(FINT(J),J=1,3)
       PRINT 100
       GOTO 900

  760  WRITE(60,100)
       WRITE(60,780)
       WRITE(60,800)
       WRITE(60,850) (FINT(J),J=1,3)
       WRITE(60,100)

  780  FORMAT(5X,'   PERCENTILES    ')
  800  FORMAT(5X,'    75 TH     50 TH     25 TH')
  850  FORMAT(5X,3F10.3)

 900   IF (ISIGN.EQ.1) THEN
          ZXX=ZU(IU)
       ELSE
          ZXX=ZU(1)
       ENDIF

       IF(OUT .EQ. CHECK) THEN
          IF(ICHANGE.EQ.-1) THEN
             WRITE(6,1000) SMEAN,ERROR,ZXX      
             WRITE(6,1005)
             WRITE(6,1006)
          ELSE
             WRITE(6,1010) SMEAN,ERROR
          ENDIF
          WRITE(6,1020)
       ELSE
          IF(ICHANGE.EQ.-1) THEN
             WRITE(60,1000) SMEAN,ERROR,ZXX      
             WRITE(60,1005)
             WRITE(60,1006)
          ELSE
             WRITE(60,1010) SMEAN,ERROR
          ENDIF
          WRITE(60,1020)
       ENDIF

       PRINT *
       PRINT *
 1000  FORMAT(8X,'MEAN=',F10.3,' +/-',F6.3,'   LIMITED TO ',F10.3)
 1005  FORMAT(/10X,'SINCE A CENSORED POINT WAS CHANGED TO A DETECTION,')
 1006  FORMAT(10X,'THE MEAN ESTIMATE IS BIASED.')
 1010  FORMAT(8X,'MEAN=',F10.3,' +/-',F6.3)
 1020  FORMAT(' ')

       RETURN
       END



C
C      **********************************************************************
C      ********************* SUBROUTINE MATINV  *****************************
C      **********************************************************************
C
       SUBROUTINE MATINV(ARRAY,NVAR,DET,IK,JK,MVAR)
C
C      *                                                                    *
C      *      THIS PROGRAM IS ADOPTED FROM PHILIP R. BEVINGTON, "DATA       *
C      *      REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES', 1969 *
C      *      McGRAW HILL (NY:NY), PROGRAM B-2 P. 302. SEVERAL MINOR        *
C      *      MODIFICATIONS WERE DONE BY T. ISOBE.                          *
C      *                                                                    *
C      *      INPUT   :   ARRAY(I,J) : SYMMETRIC MATRIX                     *
C      *                    NVAR     : DIMENSION                            *
C      *      WORK    :     AMAX     : LARGEST NO. ON WORKING               *
C      *      OUTPUT  :   ARRAY(I,J) : SYMMETRIC INVERSE MATRIX OF THE      *
C      *                               ORIGINAL ARRAY(I,J)                  *
C  
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION ARRAY(MVAR,MVAR),IK(MVAR),JK(MVAR)
C
       DET=1.0
C
       DO 100 K=1,NVAR
C
C      *      FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX             *
C
          AMAX=0.0       
C
   21     DO 30 I=K,NVAR
             DO 29 J=K,NVAR
                IF(DABS(AMAX)-DABS(ARRAY(I,J)) .GT. 0.0) GOTO 30
   24           AMAX=ARRAY(I,J)
                IK(K)=I
                JK(K)=J
   29        CONTINUE
   30     CONTINUE
C
C      *      INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)        *
C
          IF(AMAX.EQ.0.0) THEN
             DET=0.0
             RETURN
          ENDIF

   41     I=IK(K)
          IF((I-K) .LT. 0) GOTO 21

          IF((I-K).GT.0) THEN
   43        DO 50 J=1,NVAR
                SAVE=ARRAY(K,J)
                ARRAY(K,J)=ARRAY(I,J)
                ARRAY(I,J)=-SAVE
   50        CONTINUE
          ENDIF

   51     J=JK(K)
          IF((J-K).LT.0) GOTO 21

          IF((J-K).GT.0) THEN
   53        DO 60 I=1,NVAR
                SAVE=ARRAY(I,K)
                ARRAY(I,K)=ARRAY(I,J)
                ARRAY(I,J)=-SAVE
   60        CONTINUE
          ENDIF
C
C      *         ACCUMULATE ELEMENTS OF INVERSE MATRIX                      *
C
   61     DO 70 I=1,NVAR
             IF(I.NE.K) ARRAY(I,K)=-ARRAY(I,K)/AMAX
   70     CONTINUE

   71     DO 80 I=1,NVAR
             DO 79 J=1,NVAR
                IF(I.EQ.K) GOTO 79
                IF(J.EQ.K) GOTO 79
                ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
   79        CONTINUE
   80     CONTINUE

   81     DO 90 J=1,NVAR
             IF(J.NE.K) ARRAY(K,J)=ARRAY(K,J)/AMAX
   90     CONTINUE

          ARRAY(K,K)=1.0/AMAX
          DET=DET*AMAX
  100  CONTINUE  
C
C      *            RESTORE ORDERING OF MATRIX                              *
C
  101  DO 130 L=1,NVAR
          K=NVAR-L+1
          J=IK(K)

          IF(J.GT.K) THEN
             DO 110 I=1,NVAR
                SAVE=ARRAY(I,K)
                ARRAY(I,K)=-ARRAY(I,J)
                ARRAY(I,J)=SAVE
  110        CONTINUE
          ENDIF

          I=JK(K)
          IF(I.GT.K) THEN
             DO 120 J=1,NVAR
                SAVE=ARRAY(K,J)
                ARRAY(K,J)=-ARRAY(I,J)
                ARRAY(I,J)=SAVE
  120        CONTINUE
          ENDIF

  130  CONTINUE

       RETURN
       END

C
C***************************************************************************
C**************************  SUBROUTINE LRANK  *****************************
C***************************************************************************
C
C
       SUBROUTINE LRANK(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R, D1, E1, R1, 
     +                  D2, E2, R2, SCORE, VAR)
C     *
C     * THIS SUBROUTINE COMPUTES THE LOGRANK STATISTIC WITH CONDITIONAL     *
C     * PERMUTATION VARIANCE (HYPERGEOMETRIC VARIANCE) FROM EQUATIONS (2.2) *
C     * AND (2.3) IN LATTA, 'A MONTE-CARLO STUDY OF SOME TWO-SAMPLE RANK    *
C     * TESTS WITH CENSORED DATA', 1981, JOURNAL OF THE AMERICAN STATISTICAL*
C     * ASSOCIATION, VOL 76, PP 713-719.                                    *
C     *                                                                     *
C     * INPUT                                                               *
C     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
C     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
C     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
C     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
C     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
C     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
C     *                                                                     *
C     * OUTPUT                                                              *
C     *     TEST    : STANDARDIZED LOGRANK STATISTIC                        *
C     *     PROB    : PROBABILITY                                           *
C     *                                                                     *
C     * OTHERS                                                              *
C     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
C     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
C     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
C     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
C     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
C     *      R(I)   : RISK SET AT XY(I)                                     *
C     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
C     *      SCORE  : SCORE OF THE DATA                                     *
C     *      VAR    : VARIANCE OF THE DATA                                  *



       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)

       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
       DIMENSION D(NTOT),E(NTOT),R(NTOT)
       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
       DIMENSION E2(NTOT),R2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO


       I = 1
       L = 1
       R1(L) = REAL(N1)
       R2(L) = REAL(N2)
       R(L)  = REAL(NCOMP)
       ET1 = 0.0
       ET2 = 0.0

C
C     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
C     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
C
   10  IF(ID1(I) .NE. 0) THEN
          IF(ID2(I) .EQ. 1) THEN
             ET1 = ET1 + 1.0
          ELSE
             ET2 = ET2 + 1.0
          ENDIF
          I = I + 1
          GOTO 10
       ENDIF
C
C     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
C     *     FINISHED.                                                      *
C
   20  D(L)  = 0.0
       D1(L) = 0.0
       D2(L) = 0.0
       E(L)  = 0.0
       E1(L) = 0.0
       E2(L) = 0.0
       TEMP  = XY(I)
C
C     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
C     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
C     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
C     * APPROPRIATE GROUP.                                                *
C

  30   IF(ID1(I) .EQ. 0) THEN
          IF(ID2(I) .EQ. 1) THEN
            D1(L) = D1(L) + 1.0
          ELSE
             D2(L) = D2(L) + 1.0
          ENDIF

          D(L) = D1(L) + D2(L)

C
C     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
C     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
C
        ELSE
           IF(ID2(I) .EQ. 1) THEN
              E1(L) = E1(L) + 1.0
           ELSE
              E2(L) = E2(L) + 1.0
           ENDIF
           E(L) = E1(L) + E2(L)
        ENDIF

       IF(I .LE. NCOMP) THEN
           I = I + 1
C
C     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
C     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
C     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
C     * THE NUMBER OF THE CENSORED DATA POINTS                             *
C
         IF(TEMP .EQ. XY(I)) GOTO 30
         IF(ID1(I) .NE. 0) GOTO 30

C
C     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
C
            IF(L .EQ. 1) THEN
                R1(L) = R1(L) - ET1
                R2(L) = R2(L) - ET2
                R(L)  = R1(L) + R2(L)
            ELSE
                R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
                R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
                R(L)  = R1(L) + R2(L)
            ENDIF
            L = L + 1
            GOTO 20
        ENDIF
C
C     *       COMPUTE THE SCORE AND VARIANCE                         *
C

        SCORE = 0.0
        VAR   = 0.0
        L1 = L - 1
        DO 200 I = 1, L1

           SCORE = SCORE+(D2(I)-(R2(I)*D(I))/R(I))

           IF(R(I) .GT. 1.0) THEN
              VAR = VAR + D(I)*(R2(I)/R(I))*
     +              (1.0-(R2(I)/R(I)))*((R(I)-D(I))/(R(I)-1.0))
           ENDIF

  200   CONTINUE

C
C     *        NOW COMPUTE THE LOGRANK STATISTIC                   *
C
        TEST = SCORE/DSQRT(VAR)
        PROB = 1.0 - AGAUSS(TEST)

        RETURN
        END

C
C      **********************************************************************
C      ************************SUBROUTINE MULVAR ****************************
C      **********************************************************************
C
       SUBROUTINE MULVAR(X,Y,IND,NTOT,ICOL,NVAR,NOTEST,IPROG,ICOMM,
     +                   OUTPUT,COLM,FILE,YNAME,TITLE,ND,NYC,ICENS,
     +                   NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,MVAR,
     +                   LENG,LEGWRK,IBIN,XX,IND2,ALPHA,DWRK2,
     +                   DWRK3,DWRK4,DWRK5,DWRK6,DWRK8,RWRK1,
     +                   EWRK1,AWRK1,WWRK1,WWRK2,VWRK1,VWRK2,
     +                   WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
     +                   WRK9,WRK10,WRK11,WRK12,
     +                   SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
     +                   SWRK8,SWRK9,SWRK10,SWRK11,LWRK1,LWRK2,LWRK3,
     +                   IWRK1,IWRK2,IWRK3,IWRK4,IWRK5,IWRK6,IWRK7,
     +                   IWRK8,CWRK1,CWRK2,
     +                   IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,IBWRK6,
     +                   IBWRK7,IBWRK8,IBWRK9,BWRK1,BWRK2)
C
C      *                                                                    *
C      *   THIS IS THE SUBROUTINE WHICH MANAGES CORRELATION/REGRESSION      *
C      *   PROBLEMS.                                                        *
C      *                                                                    *
C      *   INPUT                                                            *
C      *      FILE     :  NAME OF DATA FILE (9 LETTERS)                     *
C      *      TITLE    :  TITLE OF THE PROBLEM (80 LETTERS)                 *
C      *      NVAR     :  NUMBER OF VARIABLES                               *
C      *      ICOL     :  INDICATOR OF VARIABLE (<=NVAR)                    *
C      *                  IF A MULTIVARIATE PROBLEM IS NEEDED, SET ICOL=0   *
C      *      COLM     :  NAME OF THE INDEPENDENT VARIABLE                  *
C      *      YNAME    :  NAME OF THE DEPENDENT VARIABLE                    *
C      *      COMMAND  :  NAME OF THE "COMMAND" FILE                        *
C      *      OUTPUT   :  NAME OF THE OUTPUT FILE                           *
C      *      IND(1,I) :  INDICATOR OF CENSORING                            *
C      *                    IF =0,  DETECTED                                *
C      *                       =1,  Y LOWER LIMIT                           *
C      *                       =2,  X LOWER LIMIT                           *
C      *                       =3,  DOUBLE LOWER LIMIT                      *
C      *                       =4,  X UPPER LIMIT AND Y LOWER LIMIT         *
C      *                       =5,  DATA POINT IS CONFINED BETWEEN TWO      *
C      *                            VALUES                                  *
C      *                  FOR THE UPPER LIMITS, CHANGE THE SIGN             *
C      *                  2, 3, AND 4 CAN BE USED ONLY IN BHK AND SCHMITT'S *
C      *                  5 CAN BE USED ONLY IN EM ALGORITHM AND IN         *
C      *                  BINNING METHODS                                   *
C      *      X(J,I)   :  INDEPENDENT VARIABLES                             *
C      *      Y(I)     :  DEPENDENT VARIABLE                                *
C      *     IPROG(I)  :  INDICATOR OF METHODS                              *
C      *     NOTEST    :  NUMBERS OF TEST                                   *
C      *  INPUT FOR EM ALGORITHM                                            *
C      *      TOL      :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX      :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      IBET     :  IF 0, NO DEPENDENT VARIABLE IS CONFINED BETWEEN   *
C      *                        TWO VALUES                                  *
C      *                     1, THERE ARE SOME DEPENDENT VARIABLE WHICH     *
C      *                        ARE CONFINED BETWEEN TWO VALUES             *
C      *    ALPHA(K)   :  INITIAL ESTIMATE OF REGRESSION COEFFICIENTS       *
C      *  INPUTS FOR BUCKLEY-JAMES METHOD                                   *
C      *      TOL1     :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX1     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *  INPUTS FOR SCHMITT'S BINNING METHOD                               *
C      *      MX       :  BIN NUMBER OF X AXES                              *
C      *      MY       :  BIN NUMBER OF Y AXES                              *
C      *      TOL3     :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                  *
C      *      MAX3     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      XBIN     :  BIN SIZE FOR X AXES                               *
C      *      YBIN     :  BIN SIZE FOR Y AXES                               *
C      *      XORG     :  ORIGN OF X AXES                                   *
C      *      YORG     :  ORIGN OF Y AXES                                   *
C      *      ISKIP    :  IF 0, THE PROGRAM WILL PROVIDE XBIN, YBIN, XORG,  *
C      *                        AND YORG.                                   *
C      *                    >0, THESE VALUES MUST BE PROVIDED BY THE USER   *
C      *      IPIRNT   :  IF 0, NO TWO DIMENSIONAL K-M ESTIMATOR WILL BE    *
C      *                        PRINTED                                     *
C      *                    >0, TWO DIMENSIONAL K-M ESTIMATOR WILL BE       *
C      *                        PRINTED                                     *
C      *                                                                    *
C      *    WORKING ARRAYS:                                                 *
C      *      NTOT     :  NUMBER OF DATA POINTS                             *
C      *      ND       :  NUMBER OF DETECTED POINTS                         *
C      *      NC1      :  NUMBER OF Y LOWER LIMITS                          *
C      *      NC2      :  NUMBER OF X LOWER LIMITS                          * 
C      *      NC3      :  NUMBER OF DOUBLE LOWER LIMITS                     * 
C      *      NC4      :  NUMBER OF Y LOWER AND X UPPER LIMITS              *
C      *      NC5      :  NUMBER OF Y UPPER LIMITS                          *
C      *      NC6      :  NUMBER OF X UPPER LIMITS                          *
C      *      NC7      :  NUMBER OF DOUBLE UPPER LIMITS                     *
C      *      NC8      :  NUMBER OF Y UPPER AND X LOWER LIMITS              *
C      *      ICENS    :  IF 0, CENSORING IS MIXED                          *
C      *                     1, CENSORING IS Y LOWER LIMITS ONLY            *
C      *                    -1, CENSORING IS Y UPPER LIMITS ONLY            *
C      *      NYC      :  NC1+NC2                                           *
C      *      NXC      :  NC3+NC4                                           *
C      *      NBC      :  NC5+NC6+NC7+NC8                                   *
C      *      IDO      :  NXC+NBC                                           *
C      *      IMUL     :  INDICATOR OF MULTIVARIATE PROBLEM                 *
C      *      XX(J,I)  :  =X(ICOL,I), EXCEPT FOR MULTI INDEPENDENT VARIABLE *
C      *                  CASE (J=1,NVAR).                                  *
C      *      IND2(I)  :  =IND(1,I)                                         *
C      *                                                                    *
C      *  OUTPUT                                                            *
C      *     COXREG                                                         *
C      *      CHI      : GLOBAL CHI-SQUARE                                  *
C      *      PROB     : PROBABILITY FOR NULL HYPOTHESIS                    *
C      *     BHK                                                            *
C      *       Z       : DEVIATION                                          *
C      *      PROB     : PROBABILITY FOR NULL HYPOTHESIS                    *
C      *     EM ALGORITHM                                                   *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS  (K=1,NVAR+1)       *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     ITE       : NUMBER OF ITERATIONS                               *
C      *     BUKLY-JAMES                                                    *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS (K=1,NVAR+1)        *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     SCHMITT                                                        *
C      *     ALPHA     : INTERCEPT COEFFICIENT                              *
C      *     BETA      : SLOPE COEFFICIENT                                  *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *     R3, R4, R5, R6, XDATA, COXREG, BHK, EM, BJ, TWOKM, SPRMAN      *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT),Y(NTOT),IND(MVAR,NTOT)
       DIMENSION IND2(NTOT),IPROG(NOTEST),ALPHA(MVAR)

       DIMENSION DWRK2(MVAR,NTOT),DWRK3(MVAR,NTOT)
       DIMENSION DWRK4(MVAR,NTOT),DWRK5(MVAR,NTOT),DWRK6(MVAR,NTOT)
       DIMENSION DWRK8(MVAR,NTOT)

       DIMENSION LWRK1(MVAR,NTOT),LWRK2(MVAR,NTOT),LWRK3(MVAR,NTOT)

       DIMENSION EWRK1(MVAR,MVAR),RWRK1(NTOT,MVAR)

       DIMENSION AWRK1(5,IBIN)
       DIMENSION WWRK1(LENG),WWRK2(LENG),VWRK1(LEGWRK),VWRK2(LEGWRK)

       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT),WRK5(NTOT)
       DIMENSION WRK6(NTOT),WRK7(NTOT),WRK8(NTOT),WRK9(NTOT),WRK10(NTOT)
       DIMENSION WRK11(NTOT),WRK12(NTOT)

       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
       DIMENSION SWRK5(MVAR),SWRK6(MVAR),SWRK7(MVAR),SWRK8(MVAR)
       DIMENSION SWRK9(MVAR),SWRK10(MVAR),SWRK11(MVAR)

       DIMENSION IWRK1(NTOT),IWRK2(NTOT),IWRK3(NTOT),IWRK4(NTOT)
       DIMENSION IWRK5(NTOT),IWRK6(NTOT),IWRK7(NTOT),IWRK8(NTOT)

       DIMENSION CWRK1(IBIN),CWRK2(IBIN)

       DIMENSION IBWRK1(IBIN,IBIN),IBWRK2(IBIN,IBIN),IBWRK3(IBIN,IBIN)
       DIMENSION IBWRK4(IBIN,IBIN),IBWRK5(IBIN,IBIN),IBWRK6(IBIN,IBIN)
       DIMENSION IBWRK7(IBIN,IBIN),IBWRK8(IBIN,IBIN),IBWRK9(IBIN,IBIN)
       DIMENSION BWRK1(IBIN,IBIN),BWRK2(IBIN,IBIN)

       CHARACTER*9 FILE,OUTPUT,COLM,YNAME
       CHARACTER*80 TITLE
C
       IMUL=1
       IF(ICOL.EQ.0) IMUL=NVAR
       IO=1
       IF(OUTPUT.EQ.'         ') IO=0
C
C      *     READ A FEW MORE INPUTS FOR SPRMAN, EM, BJ, AND TWOKM          *
C
       DO 10 I=1,NOTEST
          IF(IPROG(I).EQ.3)
     +         CALL R3(IPRSP,ICOMM)
          IF(IPROG(I).EQ.4) 
     +         CALL R4(TOL,MAX,IBET,ICOL,ALPHA,ICOMM,MVAR)
          IF(IPROG(I).EQ.5) 
     +         CALL R5(TOL1,MAX1,ICOMM)
          IF(IPROG(I).EQ.6) 
     +         CALL R6(MX,MY,ISKIP,IPRINT,TOL2,MAX2,
     +                 XBIN,YBIN,XORG,YORG,ICOMM,NLAST,IRAND)
   10  CONTINUE
C
C      *    ADJUST THE INPUT DATA FORMAT                                    *
C
       CALL XDATA(X,XX,IND,IND2,IMUL,ICOL,NTOT,MVAR)
C
       IF(IO.EQ.0) THEN
          WRITE(6,30)
          WRITE(6,40)
          WRITE(6,50) TITLE
          WRITE(6,30)
          WRITE(6,55) FILE
          WRITE(6,30)
       ELSE
          WRITE(60,30)
          WRITE(60,40)
          WRITE(60,50) TITLE
          WRITE(60,30)
          WRITE(60,55) FILE
          WRITE(60,30)
       ENDIF
C
   30  FORMAT('     ')
   40  FORMAT(5X,' CORRELATION AND REGRESSION PROBLEM')
   50  FORMAT(5X,' TITLE IS  ',A80)
   55  FORMAT(5X,' DATA FILE IS ',A9)  
C
   60  IF(IO.EQ.0) THEN
          PRINT *
          IF(ICOL.EQ.0) PRINT 80
          IF(IMUL.EQ.1) THEN
             PRINT 85
             PRINT 90,COLM,YNAME
          ENDIF
          WRITE(6,30)
       ELSE
          WRITE(60,30)
          IF(ICOL.EQ.0) WRITE(60,80)
          IF(IMUL.EQ.1) THEN
             WRITE(60,85)
             WRITE(60,90) COLM,YNAME
          ENDIF
          WRITE(60,30)
       ENDIF
C
   80  FORMAT(5X,'MULTIVARIATE PROBLEM')
   85  FORMAT(6X,'INDEPENDENT',6X,' DEPENDENT')
   90  FORMAT(8X,A9,' AND   ',A9)
C
  100  IF(IO.EQ.0) THEN
          PRINT *
          PRINT 110,NTOT
          IF(ICENS.NE.-1) THEN
             PRINT 120
             PRINT 130,NC1,NC2,NC3,NC4
          ELSEIF(ICENS.NE.1) THEN
  142        PRINT 140
             PRINT 130,NC5,NC6,NC7,NC8
             PRINT 30
          ENDIF
       ELSE
          WRITE(60,30)
          WRITE(60,110) NTOT
          IF(ICENS.NE.-1) THEN
             WRITE(60,120)
             WRITE(60,130) NC1,NC2,NC3,NC4
          ELSEIF(ICENS.NE.1) THEN
  102        WRITE(60,140)
             WRITE(60,130) NC5,NC6,NC7,NC8
          ENDIF
          WRITE(60,30)
       ENDIF
C
  110  FORMAT(5X,' NUMBER OF DATA POINTS : ',I5)
  120  FORMAT(5X,' LOWER LIMITS IN  Y    X    BOTH   MIX')
  130  FORMAT(19X,2I5,3X,2I5)
  140  FORMAT(5X,' UPPER LIMITS IN  Y    X    BOTH   MIX')
C
       DO 200 J=1,NOTEST
C
C      *       CALL TESTS AND COMPUTE THE RESULTS                           *
C       
          IF(IPROG(J).EQ.1) THEN
                    CALL COXREG(IND2,XX,Y,NTOT,IMUL,OUTPUT,ICENS,
     +                          EWRK1,SWRK1,SWRK2,IWRK1,IWRK2,
     +                          IWRK3,WRK1,DWRK8,SWRK3,IWRK4,IWRK5,MVAR)

          ELSEIF(IPROG(J).EQ.2) THEN
                    CALL BHK(IND2,XX,Y,NTOT,OUTPUT,
     +                           WRK1,WRK2,IWRK1,IWRK2,IWRK3,MVAR)

          ELSEIF(IPROG(J).EQ.3) THEN     
              CALL SPRMAN(IND2,XX,Y,NTOT,OUTPUT,IPRSP,MVAR,
     +                          WRK1,IWRK1,IWRK2,LWRK1,LWRK2,LWRK3,
     +                          DWRK8,DWRK2,DWRK3,DWRK4,WRK3,WRK4,
     +                          WRK5,WRK6,WRK7,WRK8,WRK9,WRK10,WRK11,
     +                          WRK12,IWRK3,IWRK4,DWRK5,DWRK6,SWRK1)

          ELSEIF(IPROG(J).EQ.4) THEN
                    CALL EM(IND2,XX,Y,NTOT,TOL,MAX,IMUL,IBET,
     +                      ND,NYC,OUTPUT,FILE,ALPHA,
     +                      RWRK1,WRK1,WWRK1,WWRK2,VWRK1,VWRK2,
     +                      EWRK1,SWRK1,LENG,LEGWRK,MVAR)

          ELSEIF(IPROG(J).EQ.5) THEN
                    CALL BJ(IND2,XX,Y,NTOT,TOL1,MAX1,IMUL,ND,
     +                      NYC,ICENS,OUTPUT,
     +                      SWRK1,SWRK2,IWRK1,IWRK2,IWRK4,
     +                      IWRK5,IWRK6,IWRK7,IWRK8,WRK1,WRK2,WRK3,
     +                      WRK4,WRK5,WRK6,WRK7,WRK8,SWRK3,SWRK4,
     +                      SWRK5,SWRK6,SWRK7,SWRK8,SWRK9,SWRK10,
     +                      SWRK11,DWRK8,EWRK1,MVAR)

          ELSEIF(IPROG(J).EQ.6) THEN
                     CALL TWOKM(IND2,XX,Y,NTOT,MX,MY,ISKIP,IPRINT,
     +                          ICENS,XBIN,YBIN,XORG,YORG,OUTPUT,
     +                          TOL2,MAX2,NLAST,IRAND,
     +                          NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,
     +                          WRK1,WRK2,IWRK1,CWRK1,CWRK2,BWRK1,
     +                          IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
     +                          IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK2,
     +                          IWRK2,IWRK3,WRK3,WRK4,SWRK1,DWRK8,
     +                          AWRK1,IBIN,MVAR)

          ENDIF
  200  CONTINUE

       RETURN
       END

C
C      **********************************************************************
C      *********************** FUNCTION PCHISQ  *****************************
C      **********************************************************************
C
       FUNCTION PCHISQ(CHISQR,NFREE)
C
C      *      THIS FUNCTION IS BASED ON PHILIP R. BEVINGTON 1969, "DATA     *
C      *      REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", 1969,*
C      *      McGRAW HILL (NY:NY), PROGRAM 10-1 P. 192,                     *
C      *      AND COMPUTES CHI-SQUARE PROBABILITY FROM THE REDUCED          *
C      *      CHI-SQUARE.                                                   *
C      *      INPUT :      CHISQR : REDUCED CHI-SQUARE                      *
C      *                   NFREE  : DEGREES OF FREEDOM                      *
C      *     OUTPUT :      PCHISQ : CHI-SQUARE PROBABILITY                  *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
C
       IF(NFREE .LE. 0) THEN
          PCHISQ=0.0
          RETURN
       ENDIF

       FREE=NFREE
       Z=CHISQR*FREE/2.0
       NEVEN=2*(NFREE/2)
       IF((NFREE-NEVEN) .LE. 0) THEN
C
C      *          THE DEGREES OF FREEDOM ARE EVEN                          *
C
          IMAX=NFREE/2
          TERM=1.0
          SUM=0.0
          DO 34 I=1,IMAX
             FI=I
             SUM=SUM+TERM
             TERM=TERM*Z/FI
   34     CONTINUE
          PCHISQ=SUM*DEXP(-Z)
          RETURN
C
C      *           THE DEGREES OF FREEDOM ARE ODD                          *
C
       ELSE
          IF((Z-25.0) .GT. 0) THEN
             Z=CHISQR*(FREE-1.0)/2.0
             IMAX=NFREE/2
             TERM=1.0
             SUM=0.0
             DO 44 I=1,IMAX
                FI=I
                SUM=SUM+TERM
                TERM=TERM*Z/FI
   44        CONTINUE
             PCHISQ=SUM*DEXP(-Z)
             RETURN
          ELSE
             PWR=FREE/2.0
             TERM=1.0
             SUM=TERM/PWR
             DO 56 I=1,1000
                FI=I
                TERM=-TERM*Z/FI
                SUM=SUM+TERM/(PWR+FI)
                IF((DABS(TERM/SUM)-0.00001) .LE. 0.0) GOTO 57
   56        CONTINUE
   57        PCHISQ=1.0-(Z**PWR)*SUM/GAMMA(PWR)
          ENDIF
       ENDIF

       RETURN
       END
C     
C
C***************************************************************************
C**************************  SUBROUTINE PETO   *****************************
C***************************************************************************
C
C
       SUBROUTINE PETO(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R,D1,E1,R1,D2,
     +                 E2, R2,F,A,FT,AT,XA,SCORE,VAR)
C     *
C     * THIS SUBROUTINE COMPUTES THE PETO-PRENTICE STATISTIC USING THE      *
C     * FORMULATION IN LATTA, 'A MONTE-CARLO STUDY OF SOME TWO-SAMPLE RANK  *
C     * TESTS WITH CENSORED DATA', 1981, JOURNAL OF THE AMERICAN STATISTICAL*
C     * ASSOCIATION, VOL 76, PP 713-719.  THE FORM USED IS FROM EQUATION    *
C     * 2.2 AND THE ASYMPTOTIC VARIANCE ESTIMATE GIVEN IN THE ARTICLE IS    *
C     * USED FOR THE VARIANCE.                                              *
C     *                                                                     *
C     * INPUT                                                               *
C     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
C     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
C     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
C     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
C     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
C     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
C     *                                                                     *
C     * OUTPUT                                                              *
C     *     TEST    : STANDARDIZED PETO-PRENTICE STATISTIC                  *
C     *     PROB    : PROBABILITY                                           *
C     *                                                                     *
C     * OTHERS                                                              *
C     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
C     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
C     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
C     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
C     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
C     *      R(I)   : RISK SET AT XY(I)                                     *
C     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
C     *      F(I)   : THE ESTIMATE OF THE SURVIVAL FUNCTION AT XY(I)        *
C     *      A(I)   : COEFFICIENT AT XY(I)                                  *
C     *      XA(I)  : SUM OF 2 X D2(I) AND E2(I)                            *
C     *      SCORE  : SCORE OF THE DATA                                     *
C     *      VAR    : VARIANCE OF THE DATA                                  *



       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)

       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
       DIMENSION F(NTOT),A(NTOT),FT(NTOT),AT(NTOT)
       DIMENSION D(NTOT),E(NTOT),R(NTOT),XA(NTOT)
       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
       DIMENSION E2(NTOT),R2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO


       J = 0
       I = 1
       L = 1
       F(I)  = 1.0
       A(I)  = 1.0
       R1(L) = REAL(N1)
       R2(L) = REAL(N2)
       R(L)  = REAL(NCOMP)
       ET1 = 0.0
       ET2 = 0.0

C
C     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
C     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
C
   10  IF(ID1(I) .NE. 0) THEN
          IF(ID2(I) .EQ. 1) THEN
             ET1 = ET1 + 1.0
          ELSE
             ET2 = ET2 + 1.0
          ENDIF
          I = I + 1
          F(I) = 1.0
          A(I) = 1.0
          GOTO 10
       ENDIF
C
C     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
C     *     FINISHED.                                                      *
C
   20  D(L)  = 0.0
       D1(L) = 0.0
       D2(L) = 0.0
       E(L)  = 0.0
       E1(L) = 0.0
       E2(L) = 0.0
       TEMP  = XY(I)
       K = 0
C
C     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
C     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
C     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
C     *  APPROPRIATE GROUP.                                                *
C     * HERE, FT AND AT ARE USED, SINCE WE ASSUME THAT THERE ARE NO TIES.  *
C     * IF THERE ARE TIES IN THE DATA, FT AND AT WILL BE APPROPRIATELY     *
C     *  CONVERTED INTO THE FORM FOR TIED DATA AND PUT IN F AND A.         *
C

  30   IF(ID1(I) .EQ. 0) THEN
         IF(ID2(I) .EQ. 1) THEN
            D1(L) = D1(L) + 1.0
         ELSE
            D2(L) = D2(L) + 1.0
         ENDIF

         D(L) = D1(L) + D2(L)
         J = J + 1
         K = K + 1

         IF(L .EQ. 1) THEN
           RISK = R(L) - (ET1+ET2) - (D(L) - 1.0)
           IF(J .EQ. 1) THEN
              FT(J) = RISK/(RISK + 1.0)
              AT(J) = (RISK + 1.0)/(RISK + 2.0)
           ELSE
              FT(J) = FT(J-1)*RISK/(RISK+1.0)
              AT(J) = AT(J-1)*(RISK+1.0)/(RISK+2.0)
           ENDIF
         ELSE 
           RISK = (R(L-1)-D(L-1))-E(L-1)-(D(L)-1.0)
           FT(J) = FT(J-1)*RISK/(RISK + 1.0)
           AT(J) = AT(J-1)*(RISK + 1.0)/(RISK + 2.0)
         ENDIF

C
C     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
C     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
C
        ELSE
          IF(ID2(I) .EQ. 1) THEN
             E1(L) = E1(L) + 1.0
          ELSE
             E2(L) = E2(L) + 1.0
          ENDIF
          E(L) = E1(L) + E2(L)
        ENDIF

        IF(I .LE. NCOMP) THEN
        I = I + 1
C
C     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
C     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
C     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
C     * THE NUMBER OF THE CENSORED DATA POINTS                             *
C
        IF(TEMP .EQ. XY(I)) GOTO 30
        IF(ID1(I) .NE. 0) GOTO 30

C
C     * IF THE DATA POINTS WERE TIED, K > 1.  COMPUTE THE AVERAGE OF       *
C     * FT AND AT BETWEEN JJ= 1, K AND USE THE AVERAGES FOR F AND A OF THE *
C     * DATA POINT.                                                        *
C
         SUM1 = 0.0
         SUM2 = 0.0
         JSTART = J - K + 1
         DO 50 JJ = JSTART, J
            SUM1 = SUM1 + FT(JJ)
            SUM2 = SUM2 + AT(JJ)
   50    CONTINUE

         F(L)  = SUM1/FLOAT(K)
         A(L)  = SUM2/FLOAT(K)
         XA(L) = 2.0*D2(L) + E2(L)
C
C     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
C
         IF(L .EQ. 1) THEN
             R1(L) = R1(L) - ET1
             R2(L) = R2(L) - ET2
             R(L)  = R1(L) + R2(L)
         ELSE
             R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
             R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
             R(L)  = R1(L) + R2(L)
         ENDIF
         L = L + 1
         GOTO 20
      ENDIF
C
C     *       COMPUTE THE SCORE AND VARIANCE                         *
C

      SCORE = 0.0
      VAR   = 0.0
      L1 = L - 1
      DO 200 I = 1, L1

         SCORE = SCORE +(2.0*F(I)-1.0)*D2(I)+(F(I)-1.0)*E2(I)

         SUM = 0.0
         JSTART = I + 1
         DO 100 J = JSTART, L
            SUM = SUM + F(J)*XA(J)
  100    CONTINUE

         VAR = VAR + F(I)*(1.0 - A(I))*XA(I)
     +               - (A(I) - F(I))*XA(I)*(F(I)*XA(I) + 2.0*SUM)

  200    CONTINUE

C
C     *        NOW COMPUTE THE PETO-PRENTICE STATISTIC                   *
C
        TEST = SCORE/DSQRT(VAR)
        PROB = 1.0 - AGAUSS(TEST)

        RETURN
        END

C
C      **********************************************************************
C      ********************* SUBROUTINE PLESTM  *****************************
C      **********************************************************************
C
       SUBROUTINE PLESTM(U,C,NU,NC,S,V,NTOT,SMEAN,SIGMA,ICHANGE,
     +                   NCHANGE,L)
C       
C      *      THIS SUBROUTINE COMPUTES PL ESTIMATOR AND THE MEAN            *
C      *      AND ITS ERROR.                                                *
C      *                                                                    *
C      *       INPUT     U : UNCENSORED DATA POINTS                         *
C      *                 C : CENSORED DATA POINTS                           *
C      *                NU : NO. OF UNCENSORED DATA POINTS                  *
C      *                NC : NO. OF CENSORED DATA POINTS                    *
C      *               NTOT: TOTAL NUMBER OF DATA POINTS                    *
C      *                                                                    *
C      *       WORK      L : RANK OF THE UNCENSORED DATA                    *
C      *               VAR : VARIANCE OF THE MEAN                           *
C      *                KD : NUMBER OF TIED DATA POINTS                     *
C      *                                                                    *
C      *       OUTPUT    S : PL ESTIMATOR                                   *
C      *                 V : ERROR FOR THE PL ESTIMATOR                     *
C      *             SMEAN : MEAN OF THE DATA                               *
C      *             SIGMA : ERROR OF THE MEAN                              *
C      *            ICHANGE: IF THE LAST VALUE IS CENSORED, WE NEED TO      *
C      *                     CHANGE IT TO A DETECTION. THEN ICHANGE=-1,     *
C      *                     OTHERWISE ICHANGE=1.                           *
C      *            NCHANGE: IF ICHANGE = -1 AND THE LAST VALUE IS TIED     *
C      *                     WITH OTHER CENSORED VALUES, THIS RECORDS THE   *
C      *                     NUMBER OF TIED OBSERVATIONS (ALL OF THEM NEED  *
C      *                     TO BE CHANGED TO DETECTIONS).                  *
C      *                                                                    *
C      *       FIRST HALF OF THE PROGRAM IS FROM ELISA T. LEE, "STATISTICAL *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
C      *       PUBLICATIONS (BELMONT:CA); WITH THE GRAPHIC ROUTINES REMOVED.*
C      *       FORMULAS USED FOR COMPUTATION OF THE MEAN AND THE ERROR ARE  *
C      *       FROM RUPERT G. MILLER, "SURVIVAL ANALYSIS", 1981,            *
C      *       JOHN WILEY & SONS (NY:NY)                                    *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION U(NTOT),C(NTOT),S(NTOT),V(NTOT),L(NTOT)
C
C      *          COMPUTE THE RANK (L) OF THE UNCENSORED POINTS             *
C
C*******     IF THE LAST VALUE IS CENSORED, CHANGE IT TO A DETECTION        *
C

C      THE FOLLOWING LOOP HAS BEEN MODIFIED AND NCHANGE ADDED TO THE 
C      PROGRAM TO COVER THE CASE WHEN TIED NONDETECTIONS ARE THE LARGEST
C      VALUE.  MODIFIED 4/92

       ICHANGE=1
       NCHANGE = 0
 13    IF(NC .NE. 0)THEN 
          IF(U(NU) .LE. C(NC))THEN 
             U(NU+1)=C(NC)
             NU=NU+1
             NC=NC-1
             NCHANGE = NCHANGE + 1
             ICHANGE=-1
          ELSE
             GOTO 15
          ENDIF
       ELSE
          GOTO 15
       ENDIF
       GOTO 13
C
 15    K=1
       KK=0
       NT=NU+NC
       IF(NC .NE. 0) THEN 
          DO 10 I=1,NU
             IF(KK .NE. NC) THEN
                DO 20 J=K,NC
                   K1=J
                   IF(C(J) .GE. U(I)) GOTO 1
                   KK=KK+1
  20            CONTINUE
             ENDIF
   1         K=K1
             L(I)=I+KK
  10      CONTINUE
       ELSE
          DO 19 I=1,NU
             L(I)=I
  19      CONTINUE
       ENDIF
C
C      *       COMPUTE P(T) FOR ALL UNCENSORED POINTS BASED ON RANK (L)     *
C
       V1=0.0
       PT=1.0
       XNT=NT
       DO 12 I=1,NU
          XL=L(I)
          PT=PT*((XNT-XL)/(XNT-XL+1.0))
          S(I)=PT
          IF((XNT-XL) .LE. 0.0) THEN
             VV=0.0
          ELSE
             V1=V1+1.0/((XNT-XL)*(XNT-XL+1.0))
             VV=DSQRT((PT**2)*V1)
          ENDIF
          V(I)=VV
  12   CONTINUE

C
C      *        COMPUTE THE MEAN                                            *
C      *        REF. FOR THE MEAN AND ERROR :                               *
C      *          MILLER, R. G. JR. 1981, "SURVIVAL ANALYSIS"               *
C      *          PP. 70-71 AND 195-198.                                    *
C
       SMEAN=U(1)
       I=2
  30   K=0
  40   IF((U(I+K).NE.U(I-1)).AND.(I+K.LE.NU)) THEN
          SMEAN=SMEAN+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 30
          ENDIF
       ELSEIF(U(I+K).EQ.U(I-1)) THEN
          K=K+1
          GOTO 40
       ENDIF
C
C      *              COMPUTE THE ERROR OF THE MEAN                         *
C
       J=2    
       VAR=0.0
  120  I=J
       SSUM=0
  130  K=0
  140  IF((U(I+K).EQ.U(I-1)).AND.(I+K.LE.NU)) GOTO 145
          IF(U(I+K).EQ.U(I-1)) THEN
             K=K+1
             GOTO 140
          ENDIF
  145     SSUM=SSUM+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 130
          ENDIF
C
C      *          KD IS NO. OF TIED OBSERVATIONS AT THAT POINT              *
C
       KD=1
  180  IF(U(J-1+KD).LE.U(J-1)) THEN
          KD=KD+1
          GOTO 180
       ENDIF
       XL=L(J-1)
       D=KD
       B=XNT-XL-D+1.0
C
C      *       IF THE LAST FEW DATA POINTS ARE UNCENSORED AND TIED, SKIP    *
C      *       THE NEXT LINES TO AVOID DIVISION BY 0.                       *
C
       IF(B .NE. 0.0) THEN
          VAR=VAR+SSUM*SSUM*D/((XNT-XL+1)*B)
          J=J+KD
          IF(J.LE.NU) GOTO 120
       ENDIF
  200  SIGMA=DSQRT(VAR)
   
       RETURN
       END

C
C      **********************************************************************
C      ******************** SUBROUTINE PWLCXN  ******************************
C      **********************************************************************
C
       SUBROUTINE PWLCXN(H,XM,SCORE,TEST,PROB,IWLCX,XY,ID1,ID2,NTOT)
C
C      *           THIS SUBROUTINE COMPUTES PETO AND PETO'S                 *
C      *           GENERALIZED WILCOXON STATISTIC.                          *
C      *                                                                    *
C      *           OBTAINED FROM ELISA T. LEE, "STATISTICAL METHODS FOR     *
C      *           SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING         *
C      *           PUBLICATIONS (BELMONT:CA)                                *
C      *                                                                    *
C      * SUBROUTINES                                                        *
C      *           STAT                                                     *
C
C*******           COMMON STATEMENT IS DIFFERENT FROM SMSDA.                *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION H(NTOT),XM(NTOT),SCORE(NTOT)
       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
       COMMON /G/ N,N1,N2,NCEN,ISIGN,IFULL,LO
C
       IWLCX=0
       IF(NCEN.EQ.0) IWLCX=1
       IF(NCEN.EQ.0) RETURN    
       L=1
       I=1
       IJK=0
C
C*******       THE NEXT LINE IS CHANGED FROM "EQ.1" TO "EQ.0".              *
C
  63   IF(ID1(I).NE.0) THEN
          IF(IJK.EQ.1) GOTO 65
          SCORE(I)=H(1)-1.0
          IF(I.EQ.N) GOTO 200
          I=I+1
          GOTO 63
       ENDIF
  62   M=INT(XM(L))

       DO 64 J=1,M
          SCORE(I)=H(L)+H(L+1)-1.0
          IF(I.EQ.N) GOTO 200
          I=I+1
  64   CONTINUE

       IJK=1
       L=L+1
       GOTO 63

  65   SCORE(I)=H(L)-1.0
       IF(I.EQ.N) GOTO 200
       I=I+1
       GOTO 63
C
C*******      THE NEXT LINE IS ADDED. ALSO THE PRINTING COMMANDS            *
C*******      ARE CHANGED.                                                  *
C
 200   CALL STATS(SCORE,TEST,XY,ID1,ID2,NTOT)
       PROB=1.0-AGAUSS(TEST)
       RETURN
       END

C****************************************************************************
C************************** SUBROUTINE R3    ********************************
C****************************************************************************

       SUBROUTINE R3(IPRSP,ICOMM)
C
C      *      THIS SUBROUTINE READS ONE SUPPLIMENTAL INPUT FOR SPEARMAN'S   *
C      *    RHO COMPUTATION.                                                *
C      *                                                                    *
C      *    SUBROUTINES                                                     *
C      *              DATA2                                                 *
C
       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
       CHARACTER*1 CHECK,CHAR(4,4)

       IF(ICOMM .NE.1) THEN
          PRINT *
          PRINT *,'ONE MORE QUESTION FOR THE SPEARMANS RHO COMPUTATION'
          PRINT *

   10     WRITE(6, 12)
   12     FORMAT('DO YOU WANT TO PRINT OUT THE RANKS (Y/N)? ')
     
          READ(5, 14) CHECK
   14     FORMAT(A1)

           IF(CHECK .EQ. 'Y' .OR. CHECK .EQ. 'y') THEN
              IPRSP = 1
           ELSEIF(CHECK .EQ. 'N' .OR. CHECK .EQ. 'n') THEN
              IPRSP = 0
           ELSE 
              GOTO 10
           ENDIF
 
        ELSE 
           READ(50,60) (CHAR(I,1), I = 1,4)
   60      FORMAT(4A1)
           CALL DATA2(CHAR,1,1,IPRSP,LIND)
           IF((LIND.NE.0).OR.((IPRSP.NE.0).AND.(IPRSP.NE.1))) THEN
              PRINT *
              PRINT *,'INFORMATION ABOUT RANK PRINTOUT IS NOT CLEAR'
              PRINT *
              STOP
           ENDIF
        ENDIF
        RETURN
        END

C
C      **********************************************************************
C      *********************** SUBROUTINE R4  *******************************
C      **********************************************************************
C
       SUBROUTINE R4(TOL,MAX,IBET,ICOL,ALPHA,ICOMM,MVAR)
C
C      * THIS SUBROUTINE IS A SUPPLEMENTAL READ-IN PROGRAM FOR THE EM       *
C      * ALGORITHM WITH A NORMAL DISTRIBUTION                               *
C      *                                                                    *
C      *  INPUT   ICOL   :  NUMBER OF INDEPENDENT VARIABLES                 *
C      *          ICOMM  :  INDICATOR OF COMMAND FILE VS. TERMINAL INPUT.   *
C      *                                                                    *
C      *  OUTPUT  TOL    :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                *
C      *          MAX    :  MAXIMUM ITERATION (DEFAULT 20)                  *
C      *          IBET   :  INDICATOR OF DATA SET TYPE (WHETHER THERE ARE   *
C      *                    DATA POINTS WHICH ARE CONFINED BETWEEN TWO      *
C      *                    VALUES)                                         *
C      *         ALPHA(J):  INITIAL ESTIMATES OF REGRESSION COEFFICIENTS    *
C      *                                                                    *
C      *  SUBROUTINES                                                       *
C      *         DATA1, DATA2                                               *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION ALPHA(MVAR)
       CHARACTER*1 CHECK,CHAR(4,4)
C  
       ICOL3=ICOL+2
       DO 5 I=1,ICOL3
          ALPHA(I)=0.0
    5  CONTINUE
C
C      *             FOR COMMAND FILE READING, GO TO 230                    *
C
       IF(ICOMM.NE.1) THEN
          PRINT *
          PRINT *,'A FEW MORE INPUTS FOR EM ALGORITHM'  
C
C      *              TOLERANCE LEVEL                                       *
C      
          TOL=1.0E-5
   10     PRINT *
          PRINT *,'DO YOU WANT TO SET'
          WRITE(6,20)
   20     FORMAT('    TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')
          READ(5,30) CHECK
   30     FORMAT(A1)

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
             PRINT *
             WRITE(6,50)
   50   FORMAT('WHAT IS THE TOLERANCE LEVEL (GIVE IN E FORMAT) ? ')
             READ(5,60) TOL
   60        FORMAT(E9.3)
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 10
          ENDIF
C
C      *        INFORMATION ABOUT IBET                                      *
C
C   70     PRINT *
C          PRINT *,'ARE THERE DATA POINTS WHICH ARE'
C          WRITE(6,80)
C   80     FORMAT('    CONFINED BETWEEN TWO VALUES  (Y/N) ? ')
C          READ(5,30) CHECK

C          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
C             IBET=1
C          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
             IBET=0
C          ELSE
C             GOTO 70
C          ENDIF
C
C      *          INITIAL ESTIMATIONS OF REGRESSION COEFFICIENTS            *
C
   90     PRINT *
          PRINT *,'DO YOU HAVE INITIAL ESTIMATES'
          WRITE(6,100)
  100     FORMAT('    FOR THE REGRESSION COEFFICIENTS (Y/N) ? ')
          READ(5,30) CHECK

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
  110        PRINT *
             WRITE(6,120)
  120        FORMAT('INTERCEPT COEFFICIENT= ')
             READ(5,130) ALPHA(1)
  130        FORMAT(F10.3)

             IF(ICOL.EQ.1) THEN
                WRITE(6,132)
  132           FORMAT('SLOPE COEFFICIENT= ')
                READ(5,130) ALPHA(2)

             ELSE
  138           ICOL2=ICOL+1
                DO 170 I=2,ICOL2
                   PRINT *
                   IF(I.EQ.1) WRITE(6,140)
                   IF(I.EQ.2) WRITE(6,150)
                   IF(I.GE.3) WRITE(6,160) I

  140           FORMAT('SLOPE COEFFICIENT FOR FIRST VARIABLE = ')
  150           FORMAT('SLOPE COEFFICIENT FOR SECOND VARIABLE = ')
  160         FORMAT('SLOPE COEFFICIENT FOR ',I4,'-TH VARIABLE = ')

                   READ(5,130) ALPHA(I)
  170           CONTINUE
             ENDIF
  171        ALPHA(ICOL+2)=1.0
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 90
          ENDIF
C
C      *            ITERATION LIMITS                                        *
C
  180     PRINT *
          PRINT *,'DO YOU WANT TO SET THE ITERATION '
          WRITE(6,190)
  190     FORMAT('     LIMIT (DEFAULT 50) (Y/N) ? ')
          READ(5,30) CHECK
          MAX=50

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
  200        PRINT *
             WRITE(6,210)
  210        FORMAT('WHAT IS THE MAXIMUM ITERATION ? ')
             CALL DATA1(MAX)
             IF(MAX.LE.0) GOTO 200
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 180
          ENDIF
C
C
C      *            READING FROM "COMMAND" FILE                             *
C
C      *            TOLERANCE LEVEL                                         *
C
       ELSE
  230     READ(50,60) TOL
          IF(TOL.LE.0.0) THEN
             PRINT *
             PRINT *,'    TOLERANCE FOR THE EM ALGORITHM IS NEGATIVE.'
             STOP
          ENDIF
C
C      *                   IBET                                             *
C
C  235     READ(50,240) (CHAR(I,1),I=1,4)
  240     FORMAT(4A1)
C          CALL DATA2(CHAR,1,1,IBET,LIND)

C          IF((LIND.NE.0) .OR. ((IBET.NE.0).AND.(IBET.NE.1))) THEN
C             PRINT *
C             PRINT *,'     IBET INDICATOR FOR CONFINED POINTS IS WRONG'
C             STOP
C          ENDIF
          IBET = 0
C
C      *      INITIAL ESTIMATES FOR REGRESSION COEFFICIENTS                 *
C
  245     ICOL2=ICOL+1
          ICOL3=ICOL+2
          READ(50,250) (ALPHA(I),I=1,ICOL2), ALPHA(ICOL3)
  250     FORMAT(12F10.3)
C
C      *      MAXIMUM ITERATION                                             *
C
          READ(50,240) (CHAR(I,1),I=1,4)
          CALL DATA2(CHAR,1,1,MAX,LIND)

          IF((LIND.NE.0) .OR. (MAX.LE.0) .OR. (MAX.GT.1000)) THEN
  252        PRINT *
          PRINT *,'  MAXIMUM ITERATION VALUE IS NOT BETWEEN 1 AND 1000.'
             STOP
          ENDIF

       ENDIF
  260  RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE R5 *********************************
C      **********************************************************************
C
       SUBROUTINE R5(TOL,MAX,ICOMM)
C
C      * THIS SUBROUTINE IS A SUPPLIMENTAL READING PROGRAM FOR THE BUCKLEY- *
C      * JAMES METHOD.                                                      *
C      *                                                                    *
C      * INPUT       ICOMM : INDICATOR OF READING METHOD                    *
C      *                                                                    *
C      * OUTPUT       TOL  : TOLERANCE (DEFAULT 1.0E-5)                     *
C      *              MAX  : MAXIMUM ITERATION (DEFAULT 50)                 *
C      * SUBROUTINES                                                        *
C      *              DATA1, DATA2                                          *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 CHECK,CHAR(4,4)
C
       TOL=1.0E-5
       MAX=50
C
C      *         FOR "COMMAND" FILE, GO TO 230                              *
C
       IF(ICOMM.NE.1) THEN
          PRINT *
          PRINT *,' A FEW MORE INPUTS FOR THE BUCKLEY-JAMES METHOD'
C
C      *            TOLERANCE LEVEL                                        *
C
   10     PRINT *
          PRINT *,'DO YOU WANT TO SET '
          WRITE(6,20)
   20     FORMAT('    TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')

          READ(5,30) CHECK
   30     FORMAT(A1)

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
   40        PRINT *
             WRITE(6,50)
   50        FORMAT('WHAT IS THE TOLERANCE LEVEL (E9.3 FORMAT) ? ')
             READ(5,60) TOL
   60        FORMAT(E9.3)
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 10
          ENDIF
C
C      *            ITERATION LIMITS                                       *
C
  180     PRINT *
          PRINT *,'DO YOU WANT TO SET '
          WRITE(6,190)
  190     FORMAT('   ITERATION LIMIT (DEFAULT 50) (Y/N) ? ')
          READ(5,30) CHECK

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
  200        PRINT *
             WRITE(6,210)
  210        FORMAT('WHAT IS THE MAXMUM ITERATION ? ')
             CALL DATA1(MAX)
             IF(MAX.LE.0) GOTO 200
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 180
          ENDIF
C
C
C      *      "COMMAND" FILE READING                                        *
C
C      *          TOLERANCE LEVEL                                           *
C    
       ELSE
  230     READ(50,60) TOL
          READ(50,240) (CHAR(I,1),I=1,4)
  240     FORMAT(4A1)
C
C      *         MAXIMUM ITERATION                                          *
C
          CALL DATA2(CHAR,1,1,MAX,LIND)
          IF((LIND.NE.0) .OR. (MAX.LE.0) .OR. (MAX.GT.1000)) THEN
  242        PRINT *
          PRINT *,'   MAXIMUM ITERATION VALUE IS NOT BETWEEN 1 AND 1000'
             STOP
          ENDIF
       ENDIF
  
       RETURN
       END

C
C      **********************************************************************
C      *********************** SUBROUTINE R6  *******************************
C      **********************************************************************
C
       SUBROUTINE R6(MX,MY,ISKIP,IPRINT,TOL,MAX,XBIN,
     +                            YBIN,XORG,YORG,ICOMM,NLAST,IRAND)
C     
C      *  THIS PROGRAM IS A SUPPLEMENTAL READING PROGRAM FOR SCHMITT`S      *
C      *  BINNED LINEAR REGRESSION METHOD.                                  *
C      *                                                                    *
C      *  INPUT   ICOMM : INDICATOR OF READING METHOD                       *
C      *                                                                    *
C      *  OUTPUT  MX    : NUMBER OF BINS ON X AXIS                          *
C      *          MY    : NUMBER OF BINS ON Y AXIS                          *
C      *          ISKIP : INDICATOR WHETHER THE USER GIVES BINNING          *
C      *                  INFORMATION                                       *
C      *          XBIN  : BIN SIZE FOR X AXIS                               *
C      *          YBIN  : BIN SIZE FOR Y AXIS                               *
C      *          XORG  : ORIGIN OF X AXIS                                  *
C      *          YORG  : ORIGIN OF Y AXIS                                  *
C      *          IPRINT: INDICATOR OF PRINTING FOR TWO DIM. KM ESTIMATOR   *
C      *                                                                    *
C      * SUBROUTINES                                                        *
C      *          DATA1, DATA2                                              *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 CHECK,CHECK1,CHAR(4,4)
C
       TOL=1.0E-5
       MAX=50
C
C      *   FOR "COMMAND" READING, GO TO 290                                 *
C
       IF(ICOMM.NE.1) THEN
C
          PRINT *
          PRINT *,'A FEW MORE INPUTS FOR SCHMITT`S BINNED REGRESSION'
C
C      *     NUMBER OF BINS ON AXES                                         *
C   
   10     PRINT *
          WRITE(6,20)
   20     FORMAT('HOW MANY BINS DO YOU WANT FOR THE X AXIS ? ')
          CALL DATA1(MX)
          IF(MX.LE.0) GOTO 10
   30     PRINT *
          WRITE(6,40)
   40     FORMAT('HOW MANY BINS DO YOU WANT FOR THE Y AXIS ? ')
          CALL DATA1(MY)
          IF(MY.LE.0) GOTO 30
          PRINT *
C
C      *          TOLERANCE LEVEL                                          *
C
   50     PRINT *
          PRINT *,'DO YOU WANT TO SET'
          WRITE(6,55)  
   55  FORMAT('   THE TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')
          READ(5,85) CHECK

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
             PRINT *
             WRITE(6,65)
   65        FORMAT('WHAT IS THE TOLERANCE LEVEL (E9.3 FORMAT) ? ')
             READ(5,66) TOL
   66        FORMAT(E9.3)
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
          GOTO 50
       ENDIF
C
C
C      *            ITERATION LIMITS                                       *
C
   67     PRINT *
          PRINT *,'DO YOU WANT TO SET THE'
          WRITE(6,68)
   68     FORMAT('    ITERATION LIMIT (DEFAULT 50) (Y/N) ? ')
          READ(5,85) CHECK

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
             PRINT *
   69        WRITE(6,70)
   70        FORMAT('WHAT IS THE MAXIMUM ITERATION ? ')
             CALL DATA1(MAX)
             IF(MAX.LE.0) GOTO 69 
          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
             GOTO 67
          ENDIF
C
   71     PRINT *
          WRITE(6,80)
   80  FORMAT('DO YOU WANT TO SET BIN SIZES AND ORIGIN  (Y/N) ? ')
          READ(5,85) CHECK
   85     FORMAT(A1)

          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
C
C      *      INFORMATION ABOUT BIN SIZES                                   *
C
             ISKIP=10
  110        PRINT *
             WRITE(6,120)
  120        FORMAT('WHAT IS THE BIN SIZE FOR THE X AXIS ? ')
             READ(5,130) XBIN
  130        FORMAT(F10.3)
             IF(XBIN.LE.0) GOTO 110

  140        PRINT *
             WRITE(6,150)
  150        FORMAT('WHAT IS THE BIN SIZE FOR THE Y AXIS ? ')
             READ(5,130) YBIN
             IF(YBIN.LE.0) GOTO 140
C
C      *     INFORMATION ABOUT THE ORIGIN                                   *
C
  180        PRINT *
             WRITE(6,190)
  190        FORMAT('WHERE IS THE ORIGIN OF THE X AXIS ? ')
             READ(5,130) XORG
             PRINT *
             WRITE(6,200)
  200        FORMAT('WHERE IS THE ORIGIN OF THE Y AXIS ? ')
             READ(5,130) YORG

          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
   90        ISKIP=0
             XBIN=0
             YBIN=0
             XORG=0
             YORG=0

          ELSE
             GOTO 67
          ENDIF
C
C
C      *     INFORMATION ABOUT PRINTOUTS                                    *
C
  220     PRINT *
          PRINT *,'DO YOU WANT TO PRINT OUT THE FINAL '
          WRITE(6,230)
  230     FORMAT('      2-DIMENSIONAL KM ESTIMATOR (Y/N) ? ')
          READ(5,85) CHECK
          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
             IPRINT=10
          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
             IPRINT=0
          ELSE
             GOTO 220
          ENDIF

C
C      *      INFORMATION ABOUT ERROR COMPUTATIONS                          *
C
C
          NLAST = 0
          PRINT *
  245     PRINT*,'DO YOU WANT TO COMPUTE THE ERRORS'
          WRITE(6,250)
  250     FORMAT('       FOR THE REGRESSION COEFFICIENT (Y/N) ? ')
          READ(5,85) CHECK
          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
             PRINT *
 251         WRITE(6,252)
 252         FORMAT('DO YOU WISH TO SET THE NUMBER OF BOOTSTRAP ',
     +              'ITERATIONS? (DEFAULT=200)')
             READ(5,85) CHECK1
             IF(CHECK1 .EQ. 'Y' .OR. CHECK1 .EQ. 'y') THEN
  253           WRITE(6, 254)
  254           FORMAT('HOW MANY BOOTSTRAP ITERATIONS DO YOU WANT ? ')
                READ(5,255) NLAST
  255           FORMAT(I4)
                IF((NLAST.LT.0) .OR. (NLAST.GT.1000)) THEN
                  WRITE(6,256) NLAST 
  256             FORMAT('NUMBER OF ITERATIONS :', I5)
                  PRINT *
             PRINT*,'           :SUGGESTED RANGE BETWEEN 100 AND 1000'
                  GOTO 251
                ENDIF
              ELSE IF(CHECK1 .EQ. 'N' .OR. CHECK1 .EQ. 'n') THEN
                  NLAST = 200
              ELSE
                  GOTO 251
              ENDIF
              PRINT *
  257         WRITE(6, 258)
  258         FORMAT('SEED FOR THE RANDOM NUMBER (NEGATIVE INT) :')
              READ(5,267) IRAND
  267         FORMAT(I5)
              IF(IRAND.GT. 0) GOTO 257
          ELSEIF(CHECK .NE. 'N'.AND. CHECK .NE. 'n') THEN
              GOTO 245
          ENDIF
C
C
C      *      "COMMAND" FILE READING                                        *
C
C      *       INFORMATION ABOUT BIN NUMBERS                                *
C
       ELSE
 290      READ(50,300) ((CHAR(I,J),I=1,4),J=1,2)
 300      FORMAT(8A1)
          CALL DATA2(CHAR,1,2,MX,LIND)
          IF((LIND.NE.0) .OR. (MX.LE.0)) THEN
  301        PRINT *
             PRINT *, ' NUMBER OF X-AXIS BINS IS WRONG.'
             STOP
          ENDIF

          CALL DATA2(CHAR,2,2,MY,LIND)
          IF((LIND.NE.0) .OR. (MY .LE. 0)) THEN
  305        PRINT *
             PRINT *, ' NUMBER OF Y-AXIS BINS IS WRONG.'
             STOP
          ENDIF
C
C      *      READ ISKIP                                                    *
C
          READ(50,300) (CHAR(I,1),I=1,4)
          CALL DATA2(CHAR,1,1,ISKIP,LIND)
          IF((LIND.NE.0) .OR. (ISKIP .LT. 0)) THEN
  308        PRINT *
             PRINT *, ' BIN SIZE/ORIGIN INDICATOR IS WRONG '
             STOP
          ENDIF
C
C      *           READ TOLERANCE LEVEL                                     *
C
  316     READ(50,317) TOL
  317     FORMAT(E9.3)
C
C                  READ ITERATION LIMIT                                     *
C
          READ(50,300) (CHAR(I,1),I=1,4)
          CALL DATA2(CHAR,1,1,MAX,LIND)
          IF((LIND.NE.0) .OR. (MAX.LE.0).OR.(MAX.GT.1000)) THEN
  319        PRINT *
             PRINT *, ' ITERATION LIMIT IS NOT BETWEEN 1 AND 1000. '
          ENDIF
C
C      *           READ XBIN, YBIN, XORG, YORG                              *
C
          READ(50,322) XBIN,YBIN
          READ(50,322) XORG,YORG
  322     FORMAT(2F10.3)
C
C      *              PRINTING INFORMATION                                  *
C
          READ(50,300) (CHAR(I,1),I=1,4)
          CALL DATA2(CHAR,1,1,IPRINT,LIND)
          IF((LIND.NE.0) .OR. (IPRINT.LT.0)) THEN
  323        PRINT *
             PRINT *, '  2-DIM KAPLAN-MEIER PRINT INDICATOR IS WRONG. '
             STOP
          ENDIF
C
C      *              ERROR ANALYSIS INFORMATION                            *
C
  340     NLAST = 0
          READ(50,300) (CHAR(I,1),I=1,4)
          CALL DATA2(CHAR,1,1,NLAST,LIND)
          IF((LIND .NE. 0) .OR. (NLAST.LT.0)) THEN
             PRINT *,'  NUMBER OF ITERATIONS FOR THE ERROR COMPUTATION'
             PRINT *,'     IS NEGATIVE.'
             STOP
          ENDIF
          IF(NLAST.GT.0) THEN
C
C      *  SET A SEED FOR THE RANDOM NUMBER GENERATOR; MUST BE A NEGATIVE    *
C      *  INTEGER.                                                          *
C
             READ(50,345) IRAND
  345        FORMAT(I5)
             IF(IRAND .GT. 0) THEN
                IRAND = -IRAND
                PRINT *,'THE SEED FOR THE RANDOM NUMBER GENERATOR WAS'
                PRINT *,'CHANGED TO A NEGATIVE VALUE.'
             ENDIF
          ENDIF

       ENDIF
       RETURN
       END


C**************************************************************************
C*************************** FUNCTION RAN1 ********************************
C**************************************************************************

      FUNCTION RAN1(IDUM)
C
C     *     THIS FUNCTION GIVES A UNIFORM RANDOM NUMBER BETWEEN 0 AND 1.   *
C     *    YOU NEED TO PROVIDE A SEED TO GENERATE THE FIRST VALUE          *
C     *    AND THE SEED MUST BE A NEGATIVE INTEGER.                        *
C     *    REF. NUMERICAL RECIPES P. 196.                                  *
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END

C*************************************************************************
C********************** SUBROUTINE REARRN  *******************************
C*************************************************************************

       SUBROUTINE REARRN(IND, X, INX, J, NTOT, MVAR)

C
C     *    THIS SUBROUTINE REARRANGES THE TIED DATA POINTS SO THAT THE  *
C     *    RIGHT-CENSORED DATA POINTS COME AFTER THE UNCENSORED VALUE.  *
C     *    E.G.  1, 2, 3, 3(LOWER), 3, 5, 6  ARE REARRANGED TO          *
C     *          1, 2, 3, 3, 3(LOWER), 5, 6.                            *
C     *                                                                 *
C     *  INPUT                                                          *
C     *     X(J, I)   : VARIABLES J = 1, 2 (MUST BE SORTED)             *
C     *     IND(J, I) : CENSORING INDICATOR                             *
C     *     INX(J, I) : POSITION INDICATOR                              *
C     *                                                                 *
C     *  OUTPUT                                                         *
C     *     REARRANGED X, IND, AND INX                                  *
C
C

       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
       DIMENSION IND(MVAR, NTOT), X(MVAR, NTOT), INX(MVAR, NTOT)
       
       I = 0
C
C     *            CHECK WHETHER THE DATA POINTS ARE TIED              *
C
   10  K = 0
   20  K = K + 1
       IF(X(J, I+K) . EQ. X(J, I)) GOTO 20
       IF(K .GE. 2) THEN
          K = K - 1
          KJ = I
C
C     *    IF TIED DATA POINTS WERE FOUND, CHECK WHETHER THEY ARE      *
C     *    DETECTIONS.                                                 *
C
   30     IF(IND(J, KJ) .NE. 0) THEN
             TX    = X(J, KJ)
             INDTX = IND(J, KJ)
             INXT  = INX(J, KJ)
C
C     *          THE CASE FOR THE UPPER LIMITS                        *
C
             IF(IND(J, KJ) .LT. 0) THEN
                IF(KJ .NE. I) THEN
                   K1 = KJ - I
                   DO 50 IL = 1, K1
                        L1 = KJ - IL +1
                        X(J, L1)   = X(J, L1-1)
                        IND(J, L1) = IND(J, L1-1)
                        INX(J, L1) = INX(J, L1-1)
   50               CONTINUE
                    X(J, I)   = TX
                    IND(J, I) = INDTX
                    INX(J, I) = INXT
                 ENDIF
              ELSE
C
C     *      THE CASE FOR THE LOWER LIMITS                              *
C
                 IF(KJ .NE. I+K) THEN
                 ICOUNT = 1
                 K1 = I + K - KJ
   55            DO 60 IL = 1, K1
                    L1 = KJ + IL - 1
                    X(J, L1)   = X(J, L1+1)
                    IND(J, L1) = IND(J, L1+1)
                    INX(J, L1) = INX(J, L1+1)
   60            CONTINUE
                 X(J, I+K)   = TX
                 IND(J, I+K) = INDTX
                 INX(J, I+K) = INXT
              
C
C     *   CHECK THAT THE VALUE JUST REPLACED IS NOT A LOWER LIMIT. IF  *
C     *   IT IS, REPEAT THE PROCEDURE.                                 *
C
                 ICOUNT = ICOUNT + 1
                 IF(ICOUNT .LE. K) THEN
                    IF(IND(J, KJ) .NE. 0) GOTO 55
                 ENDIF
              ENDIF
           ENDIF
         ENDIF
C
C     *   REPEAT UNTIL ALL DATA POINTS ARE TESTED                      *
C
         KJ =KJ +1
         IF(KJ .LE. I+K) GOTO 30
         ENDIF
         I = I + K
         IF(I .LT. NTOT) GOTO 10
      
         RETURN
         END


C
C      **********************************************************************
C      ****************** SUBROUTINE REGRES  ********************************
C      **********************************************************************
C
       SUBROUTINE REGRES(X,Y,NTOT,NVAR,ALPHA,RMUL,SIGM,R,
     +                   YFIT,IK,JK,A,XMEAN,SIGMAX,ARRAY,SIGMAA,MVAR) 
C
C
C      *     MULTIVARIABLE LINEAR REGRESSION FIT                            *
C      *                                                                    *
C      *     THIS SUBPROGRAM CALCULATES LINEAR REGRESSION COEFFS.,          *
C      *     VARIANCE,AND CORRELATION COEFFS., BASED ON PHILIP R. BEVINGTON,*
C      *     "DATA REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", *
C      *     McGRAW HILL (NY:NY), PROGRAM 9-1 OF P. 172                     *
C      *                                                                    *
C      *     PARAMETERS                                                     *
C      *      INPUT                                                         *
C      *       X       : ARRAY OF DATA POINTS FOR INDEP. VARIABLES          *
C      *       Y       : ARRAY OF DATA POINTS DEPENDENT VARIABLES           *
C      *       NTOT    : NUMBER OF PAIRS OF DATA POINTS                     *
C      *       NVAR    : NUMBER OF COEFFICIENTS                             *
C      *      OUTPUT                                                        *
C      *       A       : ARRAY OF COEFFICIENTS                              *
C      *       ALPHA   : OUTPUT FORM OF A ARRAY (=A)                        *
C      *       SIGMA   : ARRAY OF STANDARD DEVIATIONS OF COEFFS.            *
C      *       R       : ARRAY OF LINEAR CORRELAION COEFF.                  *
C      *       RMUL    : MULTIPLE LINEAR CORRELATION COEFF.                 *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *                 MATINV                                             *
C
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION X(MVAR,NTOT),Y(NTOT),YFIT(NTOT),IK(NTOT),JK(NTOT)
       DIMENSION A(MVAR),XMEAN(MVAR),SIGMAX(MVAR),ARRAY(MVAR,MVAR)
       DIMENSION R(MVAR),SIGMAA(MVAR),SIGM(MVAR),ALPHA(MVAR)
C
C      *           INITIALIZE SUMS, ARRAYS, AND OTHERS                      *
C
       SUM   =0.0
       RMUL  =0.0
       YMEAN =0.0
       SIGMA =0.0
       SIGMAO=0.0
       VARNCE=0.0
       FNPTS=REAL(NTOT)
       FREE1=REAL(NTOT-1)
       FREEN=REAL(NTOT-NVAR) 
       FREEJ=REAL(NVAR) 

       DO 17 I=1,NTOT
          YFIT(I)=0.0
   17  CONTINUE

   21  DO 28 J=1,NVAR 
          XMEAN(J)=0.0
          SIGMAX(J)=0.0
          R(J)=0.0
          A(J)=0.0
          SIGMAA(J)=0.0

          DO 27 K=1,NVAR 
             ARRAY(J,K)=0.0
   27     CONTINUE
   28  CONTINUE
C
C      *             TAKE MEANS                                             *
C
       DO 50 I=1,NTOT
          YMEAN=YMEAN+Y(I)

          DO 44 J=1,NVAR 
             XMEAN(J)=XMEAN(J)+X(J,I)
   44     CONTINUE
   50  CONTINUE

       YMEAN=YMEAN/FNPTS

       DO 53 J=1,NVAR 
          XMEAN(J)=XMEAN(J)/FNPTS
   53  CONTINUE
C
C      *           ACCUMULATE MATRICES R AND ARRAY                          *
C
       DO 67 I=1,NTOT
          SIGMA=SIGMA+(Y(I)-YMEAN)**2

          DO 66 J=1,NVAR 
             SIGMAX(J)=SIGMAX(J)+(X(J,I)-XMEAN(J))**2
             R(J)=R(J)+(X(J,I)-XMEAN(J))*(Y(I)-YMEAN)

             DO 65 K=1,J
             ARRAY(J,K)=ARRAY(J,K)+(X(J,I)-XMEAN(J))*(X(K,I)-XMEAN(K))
   65        CONTINUE
   66     CONTINUE
   67  CONTINUE

       SIGMA=DSQRT(SIGMA/FREE1)

       DO 78 J=1,NVAR 
          SIGMAX(J)=DSQRT(SIGMAX(J)/FREE1)
          R(J)=R(J)/(FREE1*SIGMAX(J)*SIGMA)
          DO 77 K=1,J
             ARRAY(J,K)=ARRAY(J,K)/(FREE1*SIGMAX(J)*SIGMAX(K))
             ARRAY(K,J)=ARRAY(J,K)
   77     CONTINUE
   78  CONTINUE
C
C      *            INVERT SYMMETRIC MATRIX                                 *
C
       CALL MATINV(ARRAY,NVAR,DET,IK,JK,MVAR)

       IF(DET.EQ.0.0) THEN
          A1=0.0
          PRINT *,'******WARNING : DETERMINANT IN MATINV IS 0.*****'
          RETURN
       ENDIF
C
C      *             CALCULATE COEFFICIENTS                                 *
C
  101  A1=YMEAN
       DO 108 J=1,NVAR 

          DO 104 K=1,NVAR 
             A(J)=A(J)+R(K)*ARRAY(J,K)
  104     CONTINUE

          A(J)=A(J)*SIGMA/SIGMAX(J)
          A1=A1-A(J)*XMEAN(J)

          DO 107 I=1,NTOT
             YFIT(I)=YFIT(I)+A(J)*X(J,I)
  107     CONTINUE
  108  CONTINUE
C
C      *            CALCULATE UNCERTAINTIES                                 *
C
       DO 113 I=1,NTOT
          YFIT(I)=YFIT(I)+A1
          VARNCE=VARNCE+(Y(I)-YFIT(I))**2
  113  CONTINUE

       VARNCE=VARNCE/FREEN

       DO 133 J=1,NVAR 
          SIGMAA(J)=ARRAY(J,J)*VARNCE/(FREE1*SIGMAX(J)**2)
          SIGMAA(J)=DSQRT(SIGMAA(J))
          RMUL=RMUL+A(J)*R(J)*SIGMAX(J)/SIGMA
  133  CONTINUE

       RMUL=DSQRT(RMUL)
       SIGMAO=VARNCE/FNPTS

       DO 145 J=1,NVAR 
          DO 144 K=1,NVAR 
             SIGMAO=SIGMAO+VARNCE*XMEAN(J)*XMEAN(K)*ARRAY(J,K)
     +                   /(FREE1*SIGMAX(J)*SIGMAX(K))
  144     CONTINUE
  145  CONTINUE

       DO 147 J=1,NVAR 
          JJ=J+1
          ALPHA(JJ)=A(J)
          SIGM(JJ)=SIGMAA(J)
  147  CONTINUE

       ALPHA(1)=A1
       SIGM(1)=DSQRT(SIGMAO)
       A(NVAR +2)=DSQRT(VARNCE)

       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE RMILLS  *****************************
C      **********************************************************************
C
       SUBROUTINE RMILLS(X,FUNC,TOL)
C
C      *      ALGORITHM AS 138.1 APPL.STATST. (1979) VOL.28. NO.2           *
C      *                                                                    *
C      *         COMPUTE THE RECIPROCAL OF MILLS RATIO                      *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DATA FPI /1.2533141/, FPII /0.7978846/
C
       FUNC=0.0
       IF(X .LT. -10.0) RETURN
       FUNC=FPII
       Y=DABS(X)
       IF(Y .LT. 0.000001) RETURN
       SGN=1.0
       IF(X.LT.0.0) SGN=-1.0

       IF(Y.LE.2.0) THEN
          S=0.0
          A=1.0
          T=Y
          R=Y
          B=Y**2
   40     A=A+2.0
          S=T
          R=R*B/A
          T=T+R
          IF(R.GT.TOL) GOTO 40
          FUNC=1.0/(FPI*DEXP(0.5*B)-SGN*T)
          RETURN
       ENDIF

  100  A=2.0
       B1=Y
       S=Y
       A1=Y**2+1.0
       A2=Y*(A1+2.0)
       B2=A1+1.0
       T=A2/B2

  140  A=A+1.0
       A0=A1
       A1=A2
       A2=Y*A1+A*A0
       B0=B1
       B1=B2
       B2=Y*B1+A*B0
       R=S
       S=T
       T=A2/B2

       IF((T-R.GT.TOL) .OR.(T-S.GT.TOL)) GOTO 140
       FUNC=T

       IF(SGN.LT.0.0) FUNC=T/(2.0*FPI*DEXP(0.5*Y**2)*T-1.0)

       RETURN
       END

C***************************************************************************
C************************ SUBROUTINE SCHMIT  *******************************
C***************************************************************************
C
C
       SUBROUTINE SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
     +                   YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
     +                   ALPHA,BETA,MM,M1,M2,M3,M4,M5,M6,M7,M8,
     +                   A,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
     +                   IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
     +                   IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
C
C
C      *                                                                    *
C      *          THIS SUBROUTINE COMPUTES LINEAR REGRESSION COEFFICIENTS   *
C      *      ALPHA AND BETA (INTERCEPT AND SLOPE) BY SCHMITT'S BINNED      *
C      *      METHOD. BECAUSE OF THE BINNING METHOD, FINER BINNING GIVES    *
C      *      BETTER RESULTS, THOUGH THE CPU TIME MAY INCREASE VERY         *
C      *      MUCH. ALSO IF THE BINS ARE TOO FINE, THE CENSORED POINTS      *
C      *      MAY NOT FALL ON THE DETECTIONS.                               *
C      *                                                                    *
C      *                                                                    *
C      *      INPUT                                                         *
C      *            X(I): INDEPENDENT VARIABLE                              *
C      *            Y(I): DEPENDENT VARIABLE                                *
C      *           NP(I): INDICATOR OF CENSORED STATUS                      *
C      *                  IF NP(I)=0  : DETECTION                           *
C      *                          =1  : Y(I) IS LOWER LIMIT                 *
C      *                          =2  : X(I) IS LOWER LIMIT                 *
C      *                          =3  : BOTH ARE LOWER LIMITS               *
C      *                          =4  : Y(I) IS LOWER AND X(I) IS UPPER     *
C      *                 FOR THE UPPER LIMITS, CHANGE THE SIGNS             *
C      *           NTOT : NUMBER OF DATA POINTS                             *
C      *          MPLONE: NUMBER OF PARAMETERS TO BE ESTIMATED              *
C      *                  (IN THIS PROGRAM, MPLONE IS ALWAYS 3)             *
C      *          MAXITS: NUMBER OF MAXIMUM ITERATION                       *
C      *           ALH  : DUMMY                                             *
C      *        DELTA(I): DELTA(2) CONTAINS TOLERANCE FOR THE               * 
C      *                  COMPUTATION OF F'S.                               *
C      *            MX  : NUMBER OF BINS IN THE INDEPENDENT VARIABLE        *
C      *            MY  : NUMBER OF BINS IN THE DEPENDENT VARIABLE          *
C      *          ISKIP : INDICATOR FOR THE BINNING. IF ISKIP=0, THE        *
C      *                  SUBROUTINE BIN WILL COMPUTE THE INFORMATION       *
C      *                  ABOUT THE BINNING INFORMATION                     *
C      *                  IF ISKIP>0, THE BINNING INFORMATION (BIN SIZES    *
C      *                  ORIGIN) MUST BE PROVIDED.                         *
C      *         ICENS  : IF THE DATA SET IS KNOWN TO :                     *
C      *                    CONTAIN LOWER LIMITS ONLY,   ICENS>0            *
C      *                    CONTAIN UPPER LIMITS ONLY,   ICENS<0            *
C      *                    CONTAIN MIXED LIMITS,        ICENS=0            *
C      *                     OR NOT KNOWN                ICENS=0            *
C      *         IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :             *
C      *           XBIN : THE BIN SIZE FOR THE X AXIS                       *
C      *           YBIN : THE BIN SIZE FOR THE Y AXIS                       *
C      *           XORG : THE ORIGIN OF X                                   *
C      *           YORG : THE ORIGIN OF Y                                   *
C      *       NN1=NC1  : NUMBER OF Y LOWER LIMITS                          *
C      *       NN2=NC2  : NUMBER OF X LOWER LIMITS                          *
C      *       NN3=NC3  : NUMBER OF DOUBLE LOWER LIMITS                     *
C      *       NN4=NC4  : NUMBER OF Y LOWER, X UPPER LIMITS                 *
C      *       NN5=NC5  : NUMBER OF Y UPPER LIMITS                          *
C      *       NN6=NC6  : NUMBER OF X UPPER LIMITS                          *
C      *       NN7=NC7  : NUMBER OF DOUBLE UPPER LIMITS                     *
C      *       NN8=NC8  : NUMBER OF Y UPPER, XLOWER LIMITS                  *
C      *         TOL    : TOLERANCE LEVEL                                   *
C      *         MAX    : MAXIMUM NUMBER OF ITERATIONS                      *
C      *                                                                    *
C      *     WORK                                                           *
C      *        F(I,J)  : NUMBER OF DATA POINTS IN THE BIN(I,J)             *
C      *                   (WEIGHTED BY CENSORED POINTS)                    *
C      *        A(I,J)  : MATRIX WHICH CONTAINS INFORMATION OF BIN          *
C      *                  POSITION  I=1,4 AND J=1,MX*MY                     *
C      *       TH(I)    : VECTOR " A(I,J)*F(I,J), I=1,4                     *
C      *        SUM     : TOTAL NUMBER OF POINTS.  THIS IS APPROXIMATELY    *
C      *                  EQUAL TO NTOT.  THE DIFFERENCE BETWEEN THE TWO    *
C      *                  VALUES DEPENDS ON THE TOLERANCE LEVEL.            *
C      *                                                                    *
C      *     OUTPUT                                                         *
C      *          ALPHA : INTERCEPT COEFFICIENT                             *
C      *           BETA : SLOPE COEFFICIENT                                 *
C      *                                                                    *
C      *     SUBROUTINES:                                                   *
C      *         GRDPRB : THE SUBROUTINE WHICH COMPUTES THE TWO-DIMENSIONAL * 
C      *                  KAPLAN-MEIER PROBABILITY OF THE BINS              *
C      *                                                                    *

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION A(5,IB),TH(5) 
       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB),F(IB,IB)

       DIMENSION IBWRK1(IB,IB),IBWRK2(IB,IB),IBWRK3(IB,IB)
       DIMENSION IBWRK4(IB,IB),IBWRK5(IB,IB),IBWRK6(IB,IB)
       DIMENSION IBWRK7(IB,IB),IBWRK8(IB,IB),IBWRK9(IB,IB)

       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
       DIMENSION BWRK1(IB,IB),SWRK1(MVAR),DWRK1(MVAR,NTOT)
       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
C
C      *                CALL SUBROUTINE GRDPRB                              *
C
       CALL GRDPRB(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,YBIN,
     +             XORG,YORG,TOL,MAX,MM,M1,M2,M3,M4,M5,M6,M7,M8,
     +             X,Y,NP,XB,YB,F,BWRK1,IBWRK1,IBWRK2,IBWRK3,
     +             IBWRK4,IBWRK5,IBWRK6,IBWRK7,IBWRK8,IBWRK9,
     +             IWRK1,IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)

C
C
C      *               MAKE MATRIX A(I,J)                                   *
C
       DO 1120 J=1,MY
          DO 1110 I=1,MX
             IJ=I+MX*(J-1)
             A(1,IJ)=XB(I)
             A(2,IJ)=XB(I)**2
             A(3,IJ)=XB(I)*YB(J)
             A(4,IJ)=YB(J)
             A(5,IJ)=YB(J)**2
 1110     CONTINUE
 1120  CONTINUE
C
C      *             COMPUTE THE VECTOR THETA : TH(I)                       *
C
       DO 1400 I=1,5
          TH(I)=0.0
 1400  CONTINUE

       DO 1430 J=1,MY
          DO 1420 I=1,MX
          IJ = I + MX*(J-1)
             DO 1410 K=1,5
                TH(K)=TH(K)+A(K,IJ)*F(I,J)*NTOT
 1410        CONTINUE
 1420     CONTINUE
 1430  CONTINUE

       SUM = 0.0
       DO 1460 I = 1, MX
          DO 1450 J = 1, MY
             SUM = SUM + F(I,J)*NTOT
 1450     CONTINUE
 1460  CONTINUE
 
C
C      *     COMPUTE THE SLOPE COEFFICIENT : BETA, AND THE INTERCEPT        *
C      *     COEFFICIENT : ALPHA.                                           *
C
       DEN=SUM*TH(2)-TH(1)**2
       BETA=(SUM*TH(3)-TH(1)*TH(4))/DEN
       ALPHA=(TH(2)*TH(4)-TH(1)*TH(3))/DEN
C
C
       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE SORT1  *****************************
C      **********************************************************************
C
       SUBROUTINE SORT1(ID,X,Y,NTOT,NVAR,INDEX,X1,MVAR)
C
C      *       BUBBLE SORT PROGRAM                                          *
C      *       THIS PROGRAM ARRANGES OBSERVATIONS IN ASCENDING ORDER        *
C      *       ALSO IF THERE ARE TIED DATA POINTS, IT CHECKS THE CENSORING  *
C      *       STATUS AND ORDERS THEM SO THAT A DETECTED POINT COMES        *
C      *       FIRST.                                                       *
C      *                                                                    *
C      *       INPUT : INDEX(I): POSITION INDICATOR                         *
C      *                 ID(I) : INDICATOR OF CENSORING                     *
C      *                 X(J,I): INDEPENDENT VARIABLE; NVAR DIMENSION       *
C      *                 Y(I)  : DEPENDENT VARIABLE                         *
C      *                 NTOT  : NUMBER OF DATA POINTS                      *
C      *                 NVAR  : NUMBER OF INDEPENDENT VARIABLE             *
C      *                                                                    *
C      *      OUTPUT :   ID, X, AND Y IN ASCENDING ORDER WITH INDEX         *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION ID(NTOT),X(MVAR,NTOT),Y(NTOT),X1(MVAR),INDEX(NTOT)
C
C      *        SORTING IN Y, ASCENDING ORDER                               *
C
       DO 10 I=1,NTOT
          INDEX(I)=I
   10  CONTINUE
C
       IF(NTOT.EQ.1) GOTO 100
C
       DO 90 I=1,NTOT
          J=NTOT-I+1
          JJ=J-1
          IF(JJ.GE.1) THEN   
C
             DO 80 K=1,JJ
                IF(Y(K).GT.Y(J)) THEN 

                   ID1=ID(J)
                   DO 50 L=1,NVAR
                      X1(L)=X(L,J)
  50               CONTINUE
                   Y1=Y(J)
                   INDX=INDEX(J)

                   ID(J)=ID(K)
                   DO 60 L=1,NVAR
                      X(L,J)=X(L,K)
  60               CONTINUE
                   Y(J)=Y(K)
                   INDEX(J)=INDEX(K)

                   ID(K)=ID1
                   DO 70 L=1,NVAR
                      X(L,K)=X1(L)
  70               CONTINUE
                   Y(K)=Y1
                   INDEX(K)=INDX
                ENDIF
  80         CONTINUE
          ENDIF
  90   CONTINUE
C
 100   RETURN
       END

C
C      **********************************************************************
C      ******************** SUBROUTINE SORT2 ********************************
C      **********************************************************************
C
       SUBROUTINE SORT2(XY, ID1, ID2, NTOT)
C
C      *                BUBBLE SORT PROGRAM.                                *
C      *       OBTAINED FROM ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL*
C      *       DATA ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS         *
C      *       (BELMONT:CA)                                                 *
C
C*******       THE COMMON STATEMENT IS DIFFERENT FROM SMSDA.                *
C
C      *       THIS SUBROUTINE WAS MODIFIED ON 4/20/90.                     *
C      *       DIMENSION DECLARATION.                                       *
C

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO

       DO 1 I = 1, NCOMP
          J=NCOMP-I+1
          JJ=J-1
          IF(JJ .GE. 1) THEN
             DO 2 K = 1, JJ
                IF(XY(K) .EQ. XY(J)) THEN
C
C      *   IF DATA POINTS ARE TIED, THEN CHECK THE CENSORSHIP OF BOTH DATA  *
C      *   POINTS. PUT CENSORED DATA POINTS AFTER DETECTED ONES.            *
C
                   IF(ID1(J)-ID1(K) .GE. 0) THEN 
                      GOTO 2
                   ELSE
                      X1=XY(J)
                      ITEM=ID1(J)
                      ICTE=ID2(J)

                      XY(J)=XY(K)
                      ID1(J)=ID1(K)
                      ID2(J)=ID2(K)

                      XY(K)=X1
                      ID1(K)=ITEM
                      ID2(K)=ICTE
                      GOTO 2
                   ENDIF
                ENDIF

   3            IF(XY(K) .GT. XY(J)) THEN
                   X1=XY(J)
                   ITEM=ID1(J)
                   ICTE=ID2(J)

                   XY(J)=XY(K)
                   ID1(J)=ID1(K)
                   ID2(J)=ID2(K)

                   XY(K)=X1
                   ID1(K)=ITEM
                   ID2(K)=ICTE
                ENDIF
   2         CONTINUE
          ENDIF
   1   CONTINUE

       RETURN
       END

C**************************************************************************
C************************* SUBROUTINE SPEARHO *****************************
C**************************************************************************

       SUBROUTINE SPEARHO(XF, NTOT, RHO, PROB, MVAR)

C
C     *     THIS SUBROUTINE COMPUTES THE SPEARMAN'S RHO STATISTIC AND   *
C     *     ITS PROBABILITY UNDER THE NULL HYPOTHESIS.                  *
C     *                                                                 *
C     *   INPUT                                                         *
C     *      XF(J, I) : RANKS OF TWO VARIABLES J = 1, 2                 *
C     *      NTOT     : NUMBER OF DATA POINTS                           *
C     *                                                                 *
C     *   OUTPUT                                                        *
C     *       RHO     : SPEARMAN'S RHO                                  *
C     *       PROB    : PROBABILITY                                     *
C     *                                                                 *
C     *   SUBROUTINE AGAUSS                                             *
C
C
       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
       DIMENSION XF(MVAR, NTOT)
       
       XSUM = 0.0
       X2SUM = 0.0
       YSUM = 0.0
       Y2SUM = 0.0
       XYSUM = 0.0
       RN = REAL(NTOT)
       
       DO 100 I = 1, NTOT
          XSUM = XSUM + XF(1,I)
          X2SUM = X2SUM + XF(1,I)**2
          YSUM = YSUM + XF(2,I)
          Y2SUM = Y2SUM + XF(2,I)**2
          XYSUM = XYSUM + XF(1,I)*XF(2,I)
  100  CONTINUE
  
       XBAR = XSUM/REAL(NTOT)
       YBAR = YSUM/REAL(NTOT)

       SXX = X2SUM - REAL(NTOT)*XBAR**2
       SYY = Y2SUM - REAL(NTOT)*YBAR**2
       SXY = XYSUM - REAL(NTOT)*XBAR*YBAR

       RHO = SXY/SQRT(SXX*SYY)
C
C******  THE NEXT APPROXIMATION IS GOOD ONLY WHEN N IS LARGE (E.G. >30)  *
C
       
       Z   = RHO*DSQRT(RN -1.0)
       
       PROB = 1.0 - AGAUSS(Z)
       
       RETURN
       END

C**************************************************************************
C************************* SUBROUTINE SPRMAN ******************************
C**************************************************************************

       SUBROUTINE SPRMAN(IND, X, Y, NTOT, OUTPUT, IPRSP, MVAR,
     +                    TEMP, INDT, INT, IXD, INX, INDX, Z1, XF, RX,
     +                    XX,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
     +                    WRK9,WRK10,IWRK1,IWRK2,DWRK1,DWRK2,SWRK1)
C
C     *    SPEARMAN'S RHO CORRELATION TEST                               *
C     *                                                                  *
C     *    THIS PROGRAM COMPUTES A CORRELATION PROBABILITY BETWEEN TWO   *
C     *    VARIABLES BY SPEARMAN'S RHO. FOR CENSORED DATA POINTS,        *
C     *    AKRITAS' RANKING METHOD IS USED.                              *
C     * REFERENCE                                                        *
C     *    PENN STATE UNIVERSITY, DEPARTMENT OF STATISTICS, TECHNICAL    *
C     *    REPORTS AND PREPRINTS SERIES, NUMBER 87, "ALIGNED RANK TESTS  *
C     *    FOR REGRESSION WITH CENSORED DATA', MICHAEL G. AKRITAS,       *
C     *    SEPTEMBER 1989.                                               *
C     *                                                                  *
C     * INPUT                                                            *
C     *      X(1,I)   : INDEPENDENT VARIABLE                             *
C     *      Y(I)     : DEPENDENT VARIABLE                               *
C     *      IND(I)   : CENSORED INDICATOR                               *
C     *      NTOT     : NUMBER OF DATA POINTS                            *
C     *      MVAR     : DIMENSION SIZE                                   *
C     *                                                                  *
C     * OUTPUT                                                           *
C     *        RHO    : SPEARMAN'S RHO                                   *
C     *        PROB   : PROBABILITY FOR NULL HYPOTHESIS                  *
C     *                                                                  *
C     * OTHER                                                            *
C     *      XX(J, I) : TWO VARIABLES J = 1, 2                           *
C     *      RX(J, I) : RANK OF J-TH VARIABLE                            *
C     *      INX(J, I): POSITION OF ORIGINAL VALUES                      *
C     *                                                                  *
C     * SUBROUTINES                                                      *
C     *      SORT1, REARRN, AKRANK, SPEARHO                              *
C     *                                                                  *
C
C
       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)

       DIMENSION X(MVAR, NTOT), Y(NTOT), RX(MVAR, NTOT), TEMP(NTOT)
       DIMENSION XF(MVAR, NTOT), IXD(MVAR, NTOT), INT(NTOT)
       DIMENSION IND(NTOT), INX(MVAR, NTOT), INDX(MVAR, NTOT)
       DIMENSION Z1(MVAR, NTOT), INDT(NTOT),XX(MVAR, NTOT)
       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT)
       DIMENSION WRK5(NTOT),WRK6(NTOT),WRK7(NTOT),WRK8(NTOT)
       DIMENSION WRK9(NTOT),WRK10(NTOT),IWRK1(NTOT),IWRK2(NTOT)
       DIMENSION DWRK1(MVAR,NTOT),DWRK2(MVAR,NTOT),SWRK1(MVAR)
       CHARACTER*9 OUTPUT
C
C
C     *            INITIALIZATION                                        *
C
                    
       IKEEP = IPRSP
       DO 40  I = 1, NTOT
          DO 35 J = 1, 2
             INX(J, I) = I
   35     CONTINUE
          XX(1, I) = X(1, I)
          XX(2, I) = Y(I)
   40  CONTINUE
C
C     *       ADJUST THE CENSORSHIP INDICATOR TO ONE DIMENSIONAL FORM     *
C
       DO 50 I = 1, NTOT
          INDX(1, I) = 0
          INDX(2, I) = 0
          IAD = IABS(IND(I))

          IF(IAD .GE. 2) THEN
              INDX(1, I) = IND(I)/IAD
          ENDIF
          IF(IAD .EQ. 4) INDX(1, I) = -INDX(1,I)

          IF(IAD .EQ. 1) THEN
             INDX(2, I) = IND(I)/IAD
          ELSEIF(IAD .GE. 3) THEN
             INDX(2, I) = IND(I)/IAD
          ENDIF
   50  CONTINUE
          
C
C     *    START COMPUTATION OF RANKS OF VARIABLES                      *
C
       DO 140 J = 1, 2
          DO 120 I = 1, NTOT
            TEMP(I) = XX(J, I)
            INDT(I) = INDX(J, I)
            Z1(1, I) = 0.0
            INT(I) = I
  120     CONTINUE
       
C
C     *        CALL SORT1 : SORT INTO ASCENDING ORDER                    *
C
          CALL SORT1(INDT, Z1, TEMP, NTOT, 1, INT, SWRK1, MVAR)
              
          DO 130 I = 1, NTOT
             XX(J, I) = TEMP(I)
             INDX(J, I) = INDT(I)
             INX(J, I) = INT(I)
  130     CONTINUE
C
C     *     REARRANGE TIED DATA POINTS SO THAT THE CENSORED POINTS COME *
C     *     AFTER THE UNCENSORED VALUE                                  *
C
       CALL REARRN(INDX, XX, INX, J, NTOT, MVAR)
       
C
C
C     *         AKRANK     : AKRITAS' RANKING METHOD                    *
C
          CALL AKRANK(INDX, XX, NTOT, J, RX, MVAR,
     +              WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,IWRK1,WRK7,
     +              DWRK1,WRK8,WRK9,WRK10,DWRK2,IWRK2,SWRK1)
       
  140  CONTINUE
C
C
C     *   REARRANGE THE RANKS TO THE ORIGINAL POSITIONS                 *
C
       DO 200 I = 1, NTOT
         
          DO 180 J = 1, 2
         
             DO 150 K = 1, NTOT
                IF(INX(J, K) .EQ. I) THEN
                   IXD(J, I) = K
                   GOTO 160
                ENDIF
  150        CONTINUE
             
  160        XF(J, I) = RX(J, IXD(J, I))
  180    CONTINUE
  200  CONTINUE

C
C     *   GET NEW VALUES OF INDX (OLD VALUES ARE OUT OF PROPER ORDER)   *
C
       DO 260 I = 1, NTOT
          INDX(1, I) = 0
          INDX(2, I) = 0
          IAD = IABS(IND(I))

          IF(IAD .GE. 2) THEN
              INDX(1, I) = IND(I)/IAD
          ENDIF
          IF(IAD .EQ. 4) INDX(1, I) = -INDX(1,I)

          IF(IAD .EQ. 1) THEN
             INDX(2, I) = IND(I)/IAD
          ELSEIF(IAD .GE. 3) THEN
             INDX(2, I) = IND(I)/IAD
          ENDIF
  260  CONTINUE


C      *    MAKE SURE THAT RANKS AGREE AT TIED POINTS                     *
       DO 330 I = 1, NTOT-1
          IT = 1
          T1 = XF(1,I)
          DO 290 J = I+1, NTOT
             IF(X(1,I).EQ.X(1,J).AND.INDX(1,I).EQ.INDX(1,J)) THEN
                IT = IT + 1
                T1 = T1 + XF(1,J)
             ENDIF
  290     CONTINUE

          IF(IT .GT. 1) THEN
             TAVG = T1/REAL(IT)
             XF(1,I) = TAVG
             DO 300 J = I+1, NTOT
                IF(X(1,I).EQ.X(1,J).AND.INDX(1,I).EQ.INDX(1,J)) THEN     
                   XF(1,J) = TAVG
                ENDIF
  300        CONTINUE
          ENDIF

          IT = 1
          T1 = XF(2,I)
          DO 310 J = I+1, NTOT
             IF(Y(I).EQ.Y(J).AND.INDX(2,I).EQ.INDX(2,J)) THEN
                IT = IT + 1
                T1 = T1 + XF(2,J)
             ENDIF
  310     CONTINUE

          IF(IT .GT. 1) THEN
             TAVG = T1/REAL(IT)
             XF(2,I) = TAVG
             DO 320 J = I+1, NTOT
                IF(Y(I).EQ.Y(J).AND.INDX(2,I).EQ.INDX(2,J)) THEN
                   XF(2,J) = TAVG
                ENDIF
  320        CONTINUE
          ENDIF 
  330  CONTINUE
      
C
C
C     *         CALL SPEARHO : SPEARMAN'S RHO CORRELATION TEST          *
C
C
       CALL SPEARHO(XF, NTOT, RHO, PROB, MVAR)
C
C
       IF(OUTPUT .EQ. '        ') THEN
          WRITE(6, 220)
          WRITE(6, 210)
          WRITE(6, 220)
          WRITE(6, 230) RHO
          WRITE(6, 240) PROB
          WRITE(6, 220)
C
          IF(IKEEP .EQ. 1) THEN
             WRITE(6, 250)
             WRITE(6,220)
             WRITE(6,225)
             DO 998 I = 1, NTOT
                WRITE(6, 997) IND(I),X(1, I),Y(I),(XF(J, I), J = 1, 2)
  998        CONTINUE
             WRITE(6,220)
          ENDIF
       ELSE
          WRITE(60, 220)
          WRITE(60, 210)
          WRITE(60, 220)
          WRITE(60, 230) RHO
          WRITE(60, 240) PROB
          WRITE(60, 220)
C
          IF(IKEEP .EQ. 1) THEN
             WRITE(60, 250)
             WRITE(60,220)
             WRITE(60,225)
             DO 999 I = 1, NTOT
                WRITE(60, 997) IND(I),X(1, I),Y(I),(XF(J, I), J = 1, 2)
  999        CONTINUE
             WRITE(60,220)
          ENDIF
       ENDIF

  210  FORMAT(5X,'CORRELATION TEST BY SPEARMAN`S RHO')
  220  FORMAT(' ')
  225  FORMAT(10X,'IND      X         Y        X RANK     Y RANK ')
  230  FORMAT(7X,'SPEARMANS RHO = ', F10.3)
  240  FORMAT(7X,'PROBABILITY   = ', F11.4,4X, 
     +       '(NULL HYPOTHESIS, ACCURATE ONLY IF N > 30)')
  250  FORMAT(10X,'INPUT DATA AND THEIR RANKS')
  997  FORMAT(10X,I4, 4F10.3)

       RETURN
       END


C
C      **********************************************************************
C      ***************** SUBROUTINE STAT  ***********************************
C      **********************************************************************
C
       SUBROUTINE STATS(SCORE,WSC,XY,ID1,ID2,NTOT)
C
C      *          GIVEN THE SCORES OF EACH OBSERVATION, THIS                *
C      *          SUBROUTINE COMPUTES THE FINAL TEST STATISTIC.             *
C      *                                                                    *
C      *   INPUT                                                            *
C      *         SCORE(I)  : SCORE VECTOR                                   *
C      *           XY(I)   : DATA                                           *
C      *           ID2(I)  : INDICATOR OF CENSORING                         *
C      *           NTOT    : NUMBER OF DATA POINTS                          *
C      *   OUTPUT                                                           *
C      *           WSC     : STATISTIC                                      *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION SCORE(NTOT)
       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO

       XN=NCOMP
       WW=0.0
       DO 26 I=1,NCOMP
          IF(ID2(I).NE.2) WW=WW+SCORE(I)
  26   CONTINUE
C
       IF((IFULL.EQ.1).AND.(LO.EQ.0)) WRITE(6,30) WW
       IF((IFULL.EQ.1).AND.(LO.EQ.1)) WRITE(60,30) WW
  30   FORMAT(10X,'WW =',F10.2)

       SUM=0.0
       DO 27 I=1,NCOMP
          SUM=SUM+SCORE(I)**2
  27   CONTINUE

       XN1=REAL(N1)
       XN2=REAL(N2)
       VAR=SUM*XN1*XN2/(XN*(XN-1.0)) 
C
       IF((IFULL.EQ.1).AND.(LO.EQ.0)) WRITE(6,32) VAR
       IF((IFULL.EQ.1).AND.(LO.EQ.1)) WRITE(60,32) VAR
  32   FORMAT(10X,'VAR=',F11.3)

       WSC=WW/DSQRT(VAR)

       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE SUMRY  ******************************
C      **********************************************************************
C
       SUBROUTINE SUMRY(U,IU,S,NTOT,FINT)
C
C      *       THIS SUBROUTINE CALCULATES AND PRINTS 75, 50, AND            *
C      *       25, PERCENTILES  OF SURVIVAL FOR A SURVIVAL CURVE.           *
C      *       S AND U ARE ARRAYS CONTAINING POINTS FOR WHICH VALUES OF THE *
C      *       SURVIVAL CURVE WERE CALCULATED, IU IS THE NUMBER OF          *
C      *       UNCENSORED POINTS.  ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
C      *       PUBLICATIONS (BELMONT:CA).                                   *
C      *                                                                    *
C      *       INPUT       U : UNCENSORED DATA                              *
C      *                   S : PL ESTIMATOR                                 *
C      *       WORK       TY : PERCENTILE INDICATOR AT 75, 50, 25           *
C      *       OUTPUT    FINT: VALUES AT 75, 50, 25 PERCENTILES             *
C      *                                                                    *
C      *     THIS SUBROUTINE IS FROM SMSDA, EXCEPT PRINTING COMMAND.        *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION U(IU),S(NTOT),TY(3),FINT(3)
C
C
       TY(1)=0.75
       TY(2)=0.50
       TY(3)=0.25

C      *       IF THE NUMBER OF DATA POINTS IS SMALLER THAN 4, NO          *
C      *       PERCENTILES CAN BE OBTAINED.                                *
C
       DO 40 I=1,3
          FINT(I)=0.0
  40   CONTINUE

       L=1
       IF(IU.LE.3) RETURN
       NN=IU-1

       DO 100 I=1,3
          IF(S(1).LT.TY(I)) THEN
             FINT(I) = U(1)-(TY(I)-S(1))/(1-S(1))*(U(1)-0)
          ELSE
             DO  90 J=L,NN
                IF((S(J).GE.TY(I)) .AND. (S(J+1).LE.TY(I))) THEN 
                   FINT(I)=U(J)-(S(J)-TY(I))/(S(J)-S(J+1))*(U(J)-U(J+1))
                   L=J+1
                   GOTO 100
                ENDIF
  90         CONTINUE
          ENDIF
 100   CONTINUE

       RETURN
       END

C
C      **********************************************************************
C      ******************* SUBROUTINE SYMINV  *******************************
C      **********************************************************************
C
       SUBROUTINE SYMINV(A,N,C,W,NULLTY,NA,NC,NW,IFAULT)
C
C      *      ALGORITHM AS 7 J.R.STATIST.SOC.C. (1968) VOL.17, NO.2         *
C      *                                                                    *
C      *        FORMS IN C( ) A LOWER TRIANGULAR MATRIX, WHICH IS A         *
C      *        GENERALISED INVERSE OF THE POSITIVE SEMI-DEFINITE SYMMETRIC *
C      *        MATRIX A( ) OF ORDER N.                                     *
C      *        C( ) MAY COINCIDE WITH A( ). NULLTY IS RETURNED AS          *
C      *        THE NULLITY OF A( ). IFAULT IS RETURNED AS 1 IF             *
C      *        N.LT.1,OTHERWISE ZERO. W( ) IS A WORK ARRAY OF              *
C      *        LENGTH AT LEAST N THAT IS ALLOCATED BY THE CALLING          *
C      *        ROUTINE.                                                    *
C
C      *        NOTE : THE VARIABLES NA,NC,AND,NW,HAVE BEEN ADDED           *
C      *               TO THE ARGUMENT LIST AND ARE USED TO                 *
C      *               DIMENSION TO ARRAYS A,C,AND W, RESPECTIVELY.         *
C      *               (BY WOLYNETZ (1979))                                 *
C      *                                                                    *
C      *        SUBROUTINES                                                 *
C      *               CHOL                                                 *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION A(NA),C(NC),W(NW)
C
       NROW=N
       IFAULT=1
       IF(NROW.GT.0) THEN
          IFAULT=0

          CALL CHOL(A,NROW,C,NULLTY,NA,NC,IFAULT)

          IF(IFAULT.EQ.0) THEN
             NN=(NROW*(NROW+1))/2
             IROW=NROW
             NDIAG=NN
   16        IF(C(NDIAG).NE.0.0) THEN
                L=NDIAG

                DO 10 I=IROW,NROW
                   W(I)=C(L)
                   L=L+I
   10           CONTINUE

                ICOL=NROW
                JCOL=NN
                MDIAG=NN
   15           L=JCOL
                X=0.0
                IF(ICOL.EQ.IROW) X=1.0/W(IROW)

                K=NROW
   13           IF(K.NE.IROW) THEN
                   X=X-W(K)*C(L)
                   K=K-1
                   L=L-1
                   IF(L.GT.MDIAG) L=L-K+1
                   GOTO 13
                ENDIF

                C(L)=X/W(IROW)
                IF(ICOL.EQ.IROW) GOTO 14
                MDIAG=MDIAG-ICOL
                ICOL=ICOL-1
                JCOL=JCOL-1
                GOTO 15
             ENDIF

             L=NDIAG
             DO 17 J=IROW,NROW
                C(L)=0.0
                L=L+J
   17        CONTINUE

   14        NDIAG=NDIAG-IROW
             IROW=IROW-1
             IF(IROW.NE.0) GOTO 16
          ENDIF
       ENDIF
       
       RETURN
       END

C
C      **********************************************************************
C      ************************ SUBROUTINE TIE  *****************************
C      **********************************************************************
C
       SUBROUTINE TIE(ID,X,Y,NTOT,NVAR,IL,IM,X1,MVAR)
C
C      *       CHECKS FOR THE EXISTENCE OF TIED DATA POINTS. IF A POINT     *
C      *       IS NOT TIED THE IT SETS IL(I)=1 AND IM(I)=1.                 *
C      * INPUT                                                              *
C      *       ID(I)    :  INDICATOR OF CENSORING                           *
C      *       X(J,I)   :  INDEPENDENT VARIABLES                            *
C      *       Y(I)     :  DEPENDENT VARIABLE                               *
C      *       NTOT     :  NUMBER OF DATA POINTS                            *
C      *       NVAR     :  NUMBER OF INDEPENDENT VARIABLES                  *
C      *                                                                    *
C      * OUTPUT   :                                                         *
C      *       ID, X, AND Y IN  ORDER SO THAT DETECTED POINTS ARE           *
C      *       LOCATED BEFORE CENSORED POINTS IF THEY ARE TIED.             *
C      *       IL(I)    :  INDICATOR OF TIES.                               *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION ID(NTOT),X(MVAR,NTOT),Y(NTOT),IL(NTOT)
       DIMENSION IM(NTOT),X1(MVAR)
C
       IL(1)=1
       IM(1)=1
       IL(2)=1
       IM(2)=1
       I=2
       J=1
  200  IF(Y(I).EQ.Y(I-1)) THEN
          IL(I)=IL(I-1)+1
          J=J+1
          DO 300 K=1,J
             L=I+1-K
             IM(L)=J
  300     CONTINUE
       ELSE
          J=1
       ENDIF
       IF(I.LT.NTOT) THEN
          I=I+1
          IL(I)=1
          IM(I)=1
          GOTO 200
       ENDIF
C
C      *     IF TIED DATA CONTAINS CENSORED POINTS, ORDER THEM SO THAT      *
C      *     A DETECTED POINT COMES FIRST.                                  *
C
       I=1
  550  J=1
       IF(IM(I).NE.1) THEN
  600     IF(ID(I+J-1).NE.0) THEN
             IF(J.GE.IM(I)) GOTO 800
             J=J+1
             GOTO 600
          ENDIF
          IF(J.NE.1) THEN
C
C      *       EXCHANGE THE DETECTED POINT AND THE CENSORED POINT           *
C                                       
             ID1=ID(I+J-1)
             DO 700 L=1,NVAR
                X1(L)=X(L,I+J-1)
  700        CONTINUE

             ID(I+J-1)=ID(I)
             DO 710 L=1,NVAR
                X(L,I+J-1)=X(L,I)
  710        CONTINUE

             ID(I)=ID1
             DO 720 L=1,NVAR
                X(L,I)=X1(L)
  720        CONTINUE

          ENDIF
       ENDIF
  800  IF(I.LT.NTOT) THEN
          I=I+J
          GOTO 550
       ENDIF

       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE TWOKM  ******************************
C      **********************************************************************
C
       SUBROUTINE TWOKM(IND,XX,YY,NTOT,MX,MY,ISKIP,IPRINT,ICENS,XBIN,
     +                  YBIN,XORG,YORG,OUTPUT,TOL,MAX,NLAST,IRAND,
     +                  NN1,NN2,NN3,NN4,NN5,NN6,NN7,NN8,X,Y,NP,XB,YB,F,
     +                  IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,IBWRK6,
     +                  IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,IWRK2,WRK1,
     +                  WRK2,SWRK1,DWRK1,AWRK1,IB,MVAR)
C
C      *                                                                    *
C      *          THIS PROGRAM COMPUTES LINEAR REGRESSION COEFFICIENTS      *
C      *      ALPHA AND BETA (INTERCEPT AND SLOPE) BY SCHMITT'S BINNED      *
C      *      METHOD. BECAUSE OF THE BINNING METHOD, FINER BINNING GIVES    *
C      *      BETTER RESULTS, THOUGH THE CPU TIME MAY INCREASE VERY MUCH.   *
C      *                                                                    *
C      *             WARNING   WARNING   WARNING   WARNING                  *
C      *                                                                    *
C      *    THE USER SHOULD BE WARNED THAT THIS PROGRAM ACTUALLY            *
C      *    CHANGES THE DATA!!  FIRST, IT REDEFINES SOME LIMITS TO          *
C      *    DETECTIONS.  IF THE BINS ARE CHOSEN TO BE TOO NARROW, THEN      *
C      *    VIRTUALLY ALL LIMITS COULD BE CHANGED.  SECOND, IT PUSHES       *
C      *    EACH LIMIT INTO THE ADJACENT BIN.  IF THE BINS ARE CHOSEN TO    *
C      *    TO BE TOO WIDE, THIS SUBSTANTIALLY ALTERS THE MEASURED VALUES.  *
C      *    THUS, THE USER MUST TREAD A FINE LINE IN CHOSING BIN SIZES.     *
C      *                                                                    * 
C      *        THE ERROR ANALYSIS IS DONE BY BOOTSTRAPPING METHODS. IF THE *
C      *    NUMBER OF DATA POINTS IS MUCH SMALLER THAN 100, THE NUMBER OF   *
C      *    ITERATIONS SHOULD BE (TOTAL NUMBER)**2.                         *
C      *                                                                    *
C      *      INPUT                                                         *
C      *            X(I): INDEPENDENT VARIABLE                              *
C      *            Y(I): DEPENDENT VARIABLE                                *
C      *           NP(I): INDICATOR OF CENSORED STATUS                      *
C      *                  IF NP(I)=0  : DETECTION                           *
C      *                          =1  : Y(I) IS LOWER LIMIT                 *
C      *                          =2  : X(I) IS LOWER LIMIT                 *
C      *                          =3  : BOTH ARE LOWER LIMITS               *
C      *                          =4  : Y(I) IS LOWER AND X(I) IS UPPER     *
C      *                 FOR THE UPPER LIMITS, CHANGE THE SIGNS             *
C      *           NTOT : NUMBER OF DATA POINTS                             *
C      *          MPLONE: NUMBER OF PARAMETERS TO BE ESTIMATED              *
C      *                  (IN THIS PROGRAM, MPLONE IS ALWAYS 3)             *
C      *          MAXITS: MAXIMUM NUMBER OF ITERATIONS                      *
C      *           ALH  : DUMMY                                             *
C      *        DELTA(I): DELTA(2) CONTAINS TOLERANCE FOR THE               * 
C      *                  COMPUTATION OF F'S.                               *
C      *            MX  : NUMBER OF BINS IN THE INDEPENDENT VARIABLE        *
C      *            MY  : NUMBER OF BINS IN THE DEPENDENT VARIABLE          *
C      *          ISKIP : INDICATOR FOR THE BINNING. IF ISKIP=0, THE        *
C      *                  SUBROUTINE BIN WILL COMPUTE THE INFORMATION       *
C      *                  ABOUT THE BINNING                                 *
C      *                  IF ISKIP>0, THE BINNING INFORMATION (BIN SIZES    *
C      *                  ORIGIN) MUST BE PROVIDED.                         *
C      *         IPRINT : INDICATOR FOR PRINTING. IF IPRINT>0, THE FINAL    *
C      *                  ESTIMATIONS OF TWO DIMENSIONAL PL ESTIMATOR       *
C      *                  WILL BE PRINTED.                                  *
C      *         ICENS  : IF THE DATA SET IS KNOWN TO :                     *
C      *                    CONTAIN LOWER LIMITS ONLY,   ICENS>0            *
C      *                    CONTAIN UPPER LIMITS ONLY,   ICENS<0            *
C      *                    CONTAIN MIXED LIMITS,        ICENS=0            *
C      *                     OR NOT KNOWN                ICENS=0            *
C      *         IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :             *
C      *           XBIN : THE BIN SIZE FOR THE X AXIS                       *
C      *           YBIN : THE BIN SIZE FOR THE Y AXIS                       *
C      *           XORG : THE ORIGIN OF X                                   *
C      *           YORG : THE ORIGIN OF Y                                   *
C      *       NN1=NC1  : NUMBER OF Y LOWER LIMITS                          *
C      *       NN2=NC2  : NUMBER OF X LOWER LIMITS                          *
C      *       NN3=NC3  : NUMBER OF DOUBLE LOWER LIMITS                     *
C      *       NN4=NC4  : NUMBER OF Y LOWER, X UPPER LIMITS                 *
C      *       NN5=NC5  : NUMBER OF Y UPPER LIMITS                          *
C      *       NN6=NC6  : NUMBER OF X UPPER LIMITS                          *
C      *       NN7=NC7  : NUMBER OF DOUBLE UPPER LIMITS                     *
C      *       NN8=NC8  : NUMBER OF Y UPPER, XLOWER LIMITS                  *
C      *         TOL    : TOLERANCE LEVEL                                   *
C      *         MAX    : MAXIMUM NUMBER OF ITERATIONS                      *
C      *         NLAST  : NUMBER OF ITERATIONS FOR THE BOOTSTRAPPING        *
C      *                  IF NLAST = 0, THE PROGRAM SKIPS THE BOOTSTRAPPING *
C      *         IRND   : A SEED FOR THE RANDOM NUMBER; IT IS ALSO USED AS  *
C      *                  A RANDOM NUMBER ITSELF.                           *
C      *         IB     : DIMENSION SIZE FOR BINS                           *
C      *         ILARGE : DIMENSION > MX*MY                                 *
C      *                                                                    *
C      *     WORK                                                           *
C      *        F(I,J)  : NUMBER OF DATA POINTS IN THE BIN(I,J)             *
C      *                   (WEIGHTED BY CENSORED POINTS)                    *
C      *        SUM     : TOTAL NUMBER OF POINTS. THIS APPROXIMATLY         *
C      *                  EQUALS NTOT. THE DIFFERENCE BETWEEN THE TWO       *
C      *                  VALUES DEPENDS ON THE TOLERANCE LEVEL.            *
C      *                                                                    *
C      *     OUTPUT                                                         *
C      *          ALPHA : INTERCEPT COEFFICIENT                             *
C      *           BETA : SLOPE COEFFICINET                                 *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *         SCHMIT : THE SUBROUTINE WHICH COMPUTES THE SCHMITT'S LINEAR*
C      *                  REGRESSION COEFFICIENTS                           *
C      *         RAN1   : UNIFORM RANDOM NUMBER GENERATOR                   *
C      *                                                                    *
C      *      REF: SCHMITT, J. H. M. M. 1985, ASTROPHYS. J. 293, 178.       *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION XX(MVAR,NTOT),IND(NTOT),YY(NTOT)
       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB)

       DIMENSION IBWRK1(IB,IB),IBWRK2(IB,IB),IBWRK3(IB,IB)
       DIMENSION IBWRK4(IB,IB),IBWRK5(IB,IB),IBWRK6(IB,IB)
       DIMENSION IBWRK7(IB,IB),IBWRK8(IB,IB),IBWRK9(IB,IB)
       DIMENSION BWRK1(IB,IB),F(IB,IB),AWRK1(5,IB)
       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
       DIMENSION SWRK1(MVAR),DWRK1(MVAR,NTOT)

       CHARACTER*9 OUTPUT
       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
C
       ISEED = IRAND

       IF(OUTPUT.EQ.'         ') THEN
          PRINT * 
          PRINT 80
          PRINT *
       ELSE
          WRITE(60,70)
          WRITE(60,80)
          WRITE(60,70)
       ENDIF

   70  FORMAT('   ')
   80  FORMAT(5X,'LINEAR REGRESSION BY SCHMITT`S METHOD')
C
       NC1=NN1
       NC2=NN2
       NC3=NN3
       NC4=NN4
       NC5=NN5
       NC6=NN6
       NC7=NN7
       NC8=NN8
C     
       DO 90 I=1,NTOT
          X(I)=XX(1,I)
          NP(I)=IND(I)
          Y(I)=YY(I)
   90  CONTINUE
C
C      *      COMPUTE SCHMITT'S LINEAR REGRESSION COEFFICIENTS             *
C
C
       CALL SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
     +             YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
     +             ALPHA,BETA,MM,M1,M2,M3,M4,M5,M6,M7,M8,
     +             AWRK1,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
     +             IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
     +             IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
      
       XOX = XORG
       YOY = YORG
       ALPHA1 = ALPHA

C
       IF(MM.NE.0) THEN
          IF(OUTPUT .EQ. '         ') THEN
             PRINT 2550
             PRINT 2600
             PRINT 2700
             PRINT 2800,M1,M2,M3,M4
          ELSE
             WRITE(60,2550)
             WRITE(60,2600)
             WRITE(60,2700)
             WRITE(60,2800) M1,M2,M3,M4
          ENDIF
       ELSE
C
C      *     FOR THE UPPER LIMIT CASES                                     *
C
          MM=M5+M6+M7+M8
          IF(MM.NE.0) THEN
             IF(OUTPUT .EQ. '         ') THEN
                PRINT 2550
                PRINT 3100
                PRINT 3200
                PRINT 2800,M5,M6,M7,M8
             ELSE
                WRITE(60,2550)
                WRITE(60,3100)
                WRITE(60,3200)
                WRITE(60,2800) M5,M6,M7,M8
             ENDIF
          ENDIF
       ENDIF

       CFRAC = (REAL(MM)/REAL(NTOT))*100.0
       IF(CFRAC .GE. 10.0) THEN
          IF(OUTPUT .EQ. '         ') THEN
             PRINT 3300,CFRAC
          ELSE
             WRITE(60,3300) CFRAC
          ENDIF
       ENDIF
C
C
 2550  FORMAT(10X,'# OF CENSORED POINTS CHANGED TO DETECTIONS')
 2600  FORMAT(15X,'(FROM LOWER LIMITS) ')
 2700  FORMAT(10X,' Y ONLY    X ONLY   BOTH    X LOWER Y UPPER')
 2800  FORMAT(10X,4I8)
 3100  FORMAT(15X,'(FROM UPPER LIMITS) ')
 3200  FORMAT(10X,' Y ONLY    X ONLY   BOTH    X UPPER Y LOWER')
 3300  FORMAT(10X,' WARNING!! THE CENSORING STATUS OF ',F4.1,
     +        '% OF THE POINTS ',/10X,' HAVE BEEN CHANGED!!')

 3500  IF(IPRINT.GT.0) THEN
          IF(OUTPUT.EQ.'         ') THEN
             PRINT 3890
             PRINT 3900
             PRINT 4000
          ELSE
             WRITE(60,3890)
             WRITE(60,3900)
             WRITE(60,4000)
          ENDIF
       ENDIF
C
 3890  FORMAT('    ')
 3900  FORMAT(10X,' FINAL ESTIMATION OF TWO DIMENSIONAL KM ESTIMATORS')
 4000  FORMAT(11X,'    X AXIS         Y AXIS        NO. OF POINTS')
 4100  FORMAT(7X,3F15.3)
 4110  FORMAT(10X,'# OF DATA POINTS : ',I5,'   SUM OF F : ',F12.5)
C
C
C      *     CHANGE THE PROBABILITY OF THE BIN TO # OF DATA POINTS PER BIN. *
C      *     TO CHECK THE ACCURACY OF THE 2-DIMENSIONAL KAPLAN-MEIER        *
C      *     REDISTRIBUTION, COMPARE THE SUM OF ALL F(I,J) TO THE ORIGINAL  *
C      *     NUMBER OF DATA POINTS.  IF THEY ARE NOT EQUAL, THEN THERE IS   *
C      *     TROUBLE!!!!!!!                                                 *
C      *     IF IPRINT >0, PRINT THE FINAL # OF DATA POINTS (WEIGHTED)      *
C      *     IN EACH BIN.                                                   *
C
 3650  DO 3800 J=1,MY
          DO 3700 I=1,MX
             IF(F(I,J).NE.0.0) THEN
                F(I,J)=F(I,J)*NTOT
                IF(IPRINT.GT.0) THEN
                   IF(OUTPUT.EQ.'         ') THEN
                      PRINT 4100,XB(I),YB(J),F(I,J) 
                   ELSE
                      WRITE(60,4100) XB(I),YB(J),F(I,J)
                   ENDIF
                ENDIF
             ENDIF
 3700     CONTINUE
 3800  CONTINUE
C
       IF(OUTPUT.EQ.'         ') THEN
          PRINT *
          PRINT 4110,NTOT,SUM
       ELSE
          WRITE(60,3890)
          WRITE(60,4110) NTOT,SUM
       ENDIF

C
C      *  IF NLAST IS NOT 0, THEN COMPUTE THE ERRORS OF THE INTERCEPT       *
C      *  AND SLOPE COEFFICIENTS BY THE BOOTSTRAP METHOD.                   *
C
       IF(NLAST .NE. 0) THEN
          RLAST  = REAL(NLAST)
          RLAST1 = RLAST -1.0
          ASIGM  = 0.0
          ASIGM2 = 0.0
          BSIGM  = 0.0
          BSIGM2 = 0.0
C
C      *    START THE BOOTSTRAPPING COMPUTATION                             *
C
          DO 200 ITERAT = 1, NLAST
             DO 100 I = 1, NTOT
C
C      *     SUBROUTINE RAN1 IS A RANDOM NUMBER GENERATOR.IRAND IS A SEED   *
C      *     OF THE RANDOM NUMBER.                                          *
C      *     USING THIS RANDOM NUMBER, CREATE A HYPOTHETICAL DATA SET.      *
C
                R =  RAN1(IRAND)
                RT = R
                R =  RAN1(IRAND)
                R = DSQRT(R*RT)
C
C      *    CHOOSE THE DATA POINT NUMBERED AS "L" FROM THE DATA SET         *
C
                IF(R.EQ.1.0) THEN
                   L = NTOT
                ELSE
                   L = INT(REAL(NTOT)*R) + 1
                ENDIF

                NP(I) = IND(L)
                X(I)  = XX(1,L)
                Y(I)  = YY(L)
  100        CONTINUE
C
C      *    COMPUTE THE COEFFICIENTS FOR THE HYPOTHETICAL DATA SET          *
C
             CALL SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
     +                   YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
     +                   AL,BE,MM,M1,M2,M3,M4,M5,M6,M7,M8,
     +                   AWRK1,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
     +                   IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
     +                   IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
C
C
             ASIGM  = ASIGM  + AL
             ASIGM2 = ASIGM2 + AL**2
             BSIGM  = BSIGM  + BE
             BSIGM2 = BSIGM2 + BE**2
  200     CONTINUE
C
C      *                   COMPUTE THE ERRORS                               *
C
          ASUM = ASIGM2 - ASIGM**2/RLAST
          BSUM = BSIGM2 - BSIGM**2/RLAST
          SIGMA= DSQRT(ASUM/RLAST1)
          SIGMB= DSQRT(BSUM/RLAST1)
       ENDIF
C
C      *       START PRINTING THE RESULTS                                   *
C
       IF(NLAST .EQ. 0) THEN
          IF(OUTPUT.EQ.'         ') THEN
             PRINT 1710,MX,MY
             PRINT 1715,XBIN,YBIN
             PRINT 1720,XOX,YOY
             PRINT 1810
             PRINT 1790,ALPHA1
             PRINT 1800,BETA
          ELSE
             WRITE(60,1710) MX,MY
             WRITE(60,1715) XBIN,YBIN
             WRITE(60,1720) XOX, YOY
             WRITE(60,1810)
             WRITE(60,1790) ALPHA1
             WRITE(60,1800) BETA
          ENDIF
       ELSE
          IF(OUTPUT.EQ.'         ') THEN
             PRINT 1710,MX,MY
             PRINT 1715,XBIN,YBIN
             PRINT 1720,XOX,YOY
             PRINT 1820,NLAST,ISEED
             PRINT 1810
             PRINT 1795,ALPHA1,SIGMA
             PRINT 1805,BETA,SIGMB
          ELSE
             WRITE(60,1710) MX,MY
             WRITE(60,1715) XBIN,YBIN
             WRITE(60,1720) XOX, YOY
             WRITE(60,1820) NLAST,ISEED
             WRITE(60,1810)
             WRITE(60,1795) ALPHA1,SIGMA
             WRITE(60,1805) BETA,SIGMB
          ENDIF
       ENDIF
       
 1710  FORMAT(10X,'# OF X BINS :',I8,',  # OF Y BINS :',I8)
 1715  FORMAT(10X,'X BINSIZE   :',F8.3,'    Y BINSIZE  :',F11.3)
 1720  FORMAT(10X,'X ORIGIN    :',F8.3,'    Y ORIGIN   :',F11.3)
 1790  FORMAT(7X,'INTERCEPT COEFF.   :',F10.4)
 1795  FORMAT(7X,'INTERCEPT COEFF.   :',F10.4,'+/-',F8.4, 
     +       '(BOOTSTRAP APPROXIMATION)')
 1800  FORMAT(7X,'SLOPE COEFF.       :',F10.4)
 1805  FORMAT(7X,'SLOPE COEFF.       :',F10.4,'+/-',F8.4, 
     +       '(BOOTSTRAP APPROXIMATION)')
 1810  FORMAT('     ')
 1820  FORMAT(10X,'BOOTSTRAP ITERATIONS :',I4,
     +        '     SEED :',I8)

       RETURN
       END



C
C      **********************************************************************
C      ********************* SUBROUTINE TWOST  ******************************
C      **********************************************************************
C
       SUBROUTINE TWOST(Z,IND,ISTA,IS,NG1,NG2,NTOT,IPR,OUTPUT,
     +                  M,MVAR,NDAT,FILE,
     +                  R,XM,H,X,E1,SCORE,RISK,A,R1,R2,E,XY,ID1,ID2,
     +                  WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,WRK9,
     +                  WRK10,WRK11,WRK12,WRK13,WRK14,IWRK1)
C
C
C      *                                                                    *
C      *   THIS SUBROUTINE ORGANIZES SUBROUTINES WHICH ARE                  *
C      *   RELATED TO TWO SAMPLE TESTS.                                     *
C      *                                                                    *
C      *   INPUT       Z(I,J)  : DATA TO BE TESTED                          *
C      *               IND(I,J): INDICATOR OF CENSORING                     *
C      *               ISTA(I) : INDICATOR OF GROUP                         *
C      *                 IS    : IS-TH SUB-DATA SET                         *
C      *                 NG1   : INDICATOR OF THE FIRST GROUP               *
C      *                 NG2   : INDICATOR OF THE SECOND GROUP              *
C      *             NTOT=N    : TOTAL NUMBER OF DATA POINTS                *
C      *             IPR=IFULL : INDICATOR FOR PRINTING                     *
C      *               OUTPUT  : NAME OF OUTPUT FILE                        *
C      *               MVAR    : MAXIMUM NUMBER OF VARIABLES IN THE DATA SET*
C      *                                                                    *
C      *   OTHER VARIABLES AND PARAMETERS ARE DESCRIBED IN                  *
C      *   EACH SUBROUTINE.                                                 *
C      *   ALL ENTRIES OF 'COMMON' STATEMENT ARE COMPUTED IN                *
C      *   SUBROUTINE 'AARRAY'.                                             *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *   AARRAY     : PUTS ALL OBSERVATIONS IN ARRAY XY                   *
C      *                AND FORMS ARRAYS ID1, ID2.                          *
C      *   GEHAN      : COMPUTES GEHAN'S GENERALIZED                        *
C      *                WILCOXON TEST STATISTIC BY USING                    *
C      *                MANTEL'S PROCEDURE WITH PERMUTATION VARIANCE.       *
C      *   WLCXN      : COMPUTES GEHAN'S GENERALIZED WILCOXON               *
C      *                TEST STATISTIC USING TH COMPUTATIONAL FORM          *
C      *                THE LATTA PAPER WITH HYPERGEOMETRIC VARIANCE        *
C      *   ARISK      : COMPUTES THE RISK SET AND OTHERS                    *
C      *   LRANK      : COMPUTES THE LOGRANK TEST STATISTIC                 *
C      *   PWLCXN     : COMPUTES THE PETO AND PETO GENERALIZED WILCOXON     *
C      *                TEST STATISTIC                                      *
C      *   PETO       : COMPUTES PETO AND PRENTICE'S GENERALIZED            *
C      *                WILCOXON TEST STATISTIC.                            *
C      *                                                                    *
C      *   THESE SUBROUTINES ARE MADE BASED ON THE PROGRAMS                 *
C      *   GIVEN IN ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL DATA    *
C      *   ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA).    *
C      *                                                                    *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*9 OUTPUT,FILE

       DIMENSION Z(MVAR,NDAT),IND(MVAR,NDAT),ISTA(NTOT),R(NTOT),XM(NTOT)
       DIMENSION H(NTOT),X(NTOT),E1(NTOT),SCORE(NTOT),RISK(NTOT),A(NTOT)
       DIMENSION R1(NTOT),R2(NTOT),E(NTOT)
       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT),WRK5(NTOT)
       DIMENSION WRK6(NTOT),WRK7(NTOT),WRK8(NTOT),WRK9(NTOT),WRK10(NTOT)
       DIMENSION WRK11(NTOT),WRK12(NTOT),WRK13(NTOT),WRK14(NTOT)
       DIMENSION IWRK1(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO

       LO=1
       IF(OUTPUT.EQ.'         ') LO=0
       IFULL=IPR
       N=NTOT
C
       CALL AARRAY(Z,IND,ISTA,IS,NTOT,NDAT,MVAR,NG1,NG2,XY,
     +             ID1,ID2,IWRK1)

        IF(M.LE.1) THEN
           IF(LO.EQ.1) THEN
              WRITE(60,120)
              IF(ISIGN.EQ.1) WRITE(60,480) NCOMP,NCEN
              IF(ISIGN.EQ.-1) WRITE(60,490) NCOMP,NCEN
              WRITE(60,500) N1
              WRITE(60,510) N2
           ELSE
              PRINT *
              IF(ISIGN.EQ.1) PRINT 480,NCOMP,NCEN
              IF(ISIGN.EQ.-1) PRINT 490,NCOMP,NCEN
              PRINT 500,N1
              PRINT 510,N2
           ENDIF
        ENDIF
  480   FORMAT(6X,'# OF DATA POINTS :',I4,', # OF LOWER LIMITS :',I4)
  490   FORMAT(6X,'# OF DATA POINTS :',I4,', # OF UPPER LIMITS :',I4)
  500   FORMAT(6X,'# OF DATA POINTS IN GROUP I  :',I4)
  510   FORMAT(6X,'# OF DATA POINTS IN GROUP II :',I4)
C
       CALL ARISK(R,XM,X,E1,NG,H,XY, ID1, NTOT)

       IF(IFULL.EQ.1) THEN
           IF(LO.EQ.0) THEN
              WRITE(6,120)
              WRITE(6,604) 
           ELSE
              WRITE(60,120)
              WRITE(60,604) 
           ENDIF
 604   FORMAT(1X,'DISTINCT FAILURES',9X,'R(I)',12X,'M(I)',8X,
     +        'M(I)/R(I)',9X,'H(I)')

C
C*******      CORRECT SIGN OF X(I) BY MULTIPLYING BY ISIGN                  *
C
           DO 600 I=1,NG
              ZZ=ISIGN*X(I)
              IF(LO.EQ.0) THEN
                 WRITE(6,120)
                 WRITE(6,601) ZZ,R(I),XM(I),E1(I),H(I)
              ELSE
                 WRITE(60,120)
                 WRITE(60,601) ZZ,R(I),XM(I),E1(I),H(I)
              ENDIF
  600      CONTINUE
  120      FORMAT('     ')
  601      FORMAT(2X,5F15.3)

           NN1=NG+1
           IF(LO.EQ.0) THEN
              WRITE(6,602) H(NN1)
              WRITE(6,120)
           ELSE
              WRITE(60,602) H(NN1)
              WRITE(60,120)
           ENDIF
  602      FORMAT(62X,F15.3)

        ENDIF

C  605  LTEST=1
C
C***************   GEHAN TEST PRINTOUT -- PERMUTATION VARIANCE               *
C
        IF(LO .EQ. 0) THEN
           WRITE(6,120)
           WRITE(6,125)
           WRITE(6,120)
        ELSE
           WRITE(60,120)
           WRITE(60,125)
           WRITE(60,120)
        ENDIF
 125    FORMAT(8X,'GEHAN`S GENERALIZED WILCOXON TEST',
     +         ' -- PERMUTATION VARIANCE')

        CALL GEHAN(R1,R2,TEST,PROB,XY,ID1,ID2,NTOT)

        IF(IFULL .EQ. 1) THEN
           IF(LO .EQ. 0) THEN
              WRITE(6,120)
              WRITE(6,610)
           ELSE
              WRITE(60,120)
              WRITE(60,610)
           ENDIF

 610       FORMAT(//17X,'T(I)',7X,'ID1(I)',4X,'ID2(I)',3X,
     +            'R1(I)',4X,'R2(I)')

           DO 611 I = 1, NCOMP
C
C********  CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN
C
              ZZ = REAL(ISIGN)*XY(I)
              IF(LO .EQ. 0) THEN 
                 WRITE(6,612) ZZ,ID1(I),ID2(I),R1(I),R2(I)
              ELSE
                 WRITE(60,612) ZZ,ID1(I),ID2(I),R1(I),R2(I)
              ENDIF
 611       CONTINUE

 612       FORMAT(5X,F15.3,2I10,2F10.1)
        ENDIF

        IF(LO .EQ. 0) THEN
           WRITE(6,660) TEST
           WRITE(6,665) PROB
        ELSE
           WRITE(60,660) TEST
           WRITE(60,665) PROB
        ENDIF

C
C***************   GEHAN TEST PRINTOUT -- HYPERGEOMETRIC VARIANCE            *
C
C       IF(IPROG.EQ.1) THEN
          IF(LO.EQ.0) THEN
             WRITE(6,120)
             WRITE(6,130)
             WRITE(6,120)
          ELSE
             WRITE(60,120)
             WRITE(60,130)
             WRITE(60,120)
          ENDIF
  130     FORMAT(8X,'GEHAN`S GENERALIZED WILCOXON TEST',
     +           ' -- HYPERGEOMETRIC VARIANCE')
C
          CALL WLCXN(ID1,ID2,XY,NTOT,TEST,PROB,WRK1,WRK2,WRK3,WRK4,
     +               WRK5,WRK6,WRK7,WRK8,WRK9,WW,VAR)

          IF(IFULL.EQ.1) THEN
             IF(LO.EQ.0) THEN
                WRITE(6,641) WW
                WRITE(6,642) VAR
                WRITE(6,120)
                WRITE(6,640)
              ELSE
                WRITE(60,641) WW
                WRITE(60,642) VAR
                WRITE(60,120)
                WRITE(60,640)
              ENDIF
 640   FORMAT(1X,'            ',//15X,'T(I)',9X,'SCORE(I)')
 641   FORMAT(10X,'WW  = ',F15.3)
 642   FORMAT(10X,'VAR = ',F15.3)
C
C******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
C******       CHANGE AFFECTS ONLY THE PRINTING.                             *
C
             J = 1
             DO 646 I=1,NCOMP

 645            IF(J .GT. NCOMP) THEN
                   GOTO 646
                ENDIF

                IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
                   J=J+1
                   GOTO 645 
                ENDIF

                IF(ID1(J) .EQ. 1) THEN
                   J = J+1
                   GOTO 645
                ENDIF

                ZZ=ISIGN*XY(J)

                IF(WRK3(I) .NE. 0.0) THEN
                   SCORE(I) = WRK3(I)*(WRK7(I)-
     +                        ((WRK9(I)*WRK1(I))/WRK3(I)))
                ELSE
                   SCORE(I) = 0.0
                ENDIF

                IF(LO.EQ.0) THEN
                   WRITE(6,780) ZZ,SCORE(I)
                ELSE
                   WRITE(60,780) ZZ,SCORE(I)
                ENDIF
                J=J+1
 646         CONTINUE
           ENDIF

  652     IF(LO.EQ.0) THEN
             WRITE(6,660) TEST
             WRITE(6,665) PROB
          ELSE
             WRITE(60,660) TEST
             WRITE(60,665) PROB
          ENDIF
  660     FORMAT(10X,'TEST STATISTIC        =',F12.3)
  665     FORMAT(10X,'PROBABILITY           =',F13.4,/)

C
C************            LOGRANK TEST PRINTOUT                             *
C
C       ELSEIF(IPROG.EQ.2) THEN
          IF(LO.EQ.0) THEN
             WRITE(6,120)
             WRITE(6,150)
             WRITE(6,120)
          ELSE
             WRITE(60,120)
             WRITE(60,150)
             WRITE(60,120)
          ENDIF
  150     FORMAT(8X,'LOGRANK TEST ')
C
          CALL LRANK(ID1,ID2,XY,NTOT,TEST,PROB,WRK1,WRK2,WRK3,WRK4,
     +               WRK5,WRK6,WRK7,WRK8,WRK9,WW,VAR)

          IF(IFULL.EQ.1) THEN
             IF(LO.EQ.0) THEN
                WRITE(6,641) WW
                WRITE(6,642) VAR
                WRITE(6,120)
                WRITE(6,640)
             ELSE
                WRITE(60,641) WW
                WRITE(60,642) VAR
                WRITE(60,120)
                WRITE(60,640)
             ENDIF
 
C
C******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
C******       CHANGE AFFECTS ONLY THE PRINTING.                             *
C
             J = 1
             DO 700 I=1,NCOMP

 695            IF(J .GT. NCOMP) THEN
                   GOTO 700
                ENDIF

                IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
                   J=J+1
                   GOTO 695 
                ENDIF

                IF(ID1(J) .EQ. 1) THEN
                   J = J+1
                   GOTO 695
                ENDIF

                ZZ=ISIGN*XY(J)

                IF(WRK3(I) .NE. 0.0) THEN
                   SCORE(I) = WRK7(I)-((WRK9(I)*WRK1(I))/WRK3(I))
                ELSE
                   SCORE(I) = 0.0
                ENDIF

                IF(LO.EQ.0) THEN
                   WRITE(6,780) ZZ,SCORE(I)
                ELSE
                   WRITE(60,780) ZZ,SCORE(I)
                ENDIF
                J=J+1
 700         CONTINUE
          ENDIF
C
 703      IF(LO.EQ.0) THEN
             WRITE(6,660) TEST
             WRITE(6,665) PROB
          ELSE
             WRITE(60,660) TEST
             WRITE(60,665) PROB
          ENDIF

C
C****************  PETO-PETO PRINTOUT
C

          IF(LO .EQ. 0) THEN
             WRITE(6,120)
             WRITE(6,175)
             WRITE(6,120)
          ELSE
             WRITE(60,120)
             WRITE(60,175)
             WRITE(60,120)
          ENDIF

 175      FORMAT(8X,'PETO & PETO GENERALIZED WILCOXON TEST')
          
          CALL PWLCXN(H,XM,SCORE,TEST,PROB,IWLCXN,XY,ID1,ID2,NTOT)

          IF(IFULL .EQ. 1) THEN
             IF(LO .EQ. 0) THEN
                WRITE(6,120)
                WRITE(6,770)
             ELSE
                WRITE(60,120)
                WRITE(60,770)
             ENDIF

 770         FORMAT(27X,//15X,'T(I)',9X,'SCORE(I)')

C***********    CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN

             DO 775 I=1,NCOMP
                IF(XY(I) .EQ. 0.0) GOTO 775
                ZZ = REAL(ISIGN)*XY(I)
                IF(LO .EQ. 0) THEN
                   WRITE(6,780) ZZ,SCORE(I)
                ELSE
                   WRITE(60,780) ZZ,SCORE(I)
                ENDIF
 775         CONTINUE
 780         FORMAT(5X,2F15.3)
          ENDIF

          IF(LO .EQ. 0) THEN
             WRITE(6,660) TEST
             WRITE(6,665) PROB
          ELSE
             WRITE(60,660) TEST
             WRITE(60,665) PROB
          ENDIF

C
C****************  PETO-PRENTICE PRINTOUT                                    *
C
C  180  ELSEIF(IPROG.EQ.4)THEN
          IF(LO.EQ.0) THEN
             WRITE(6,120)
             WRITE(6,190)
             WRITE(6,120)
          ELSE
             WRITE(60,120)
             WRITE(60,190)
             WRITE(60,120)
          ENDIF
  190     FORMAT(8X,'PETO & PRENTICE GENERALIZED WILCOXON TEST')
          IF(NCEN .NE. 0) THEN
             CALL PETO(ID1, ID2, XY, NTOT, TEST,PROB, 
     +                 WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
     +                 WRK9,WRK10,WRK11,WRK12,WRK13,WRK14,WW,VAR)

             IF(IFULL.EQ.1) THEN
                IF(LO.EQ.0) THEN
                   WRITE(6,641) WW
                   WRITE(6,642) VAR
                   WRITE(6,120)
                   WRITE(6,640)
                ELSE
                   WRITE(60,641) WW
                   WRITE(60,642) VAR
                   WRITE(60,120)
                   WRITE(60,640)
                ENDIF

C
C******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
C******       CHANGE AFFECTS ONLY THE PRINTING.                             *
C
                J = 1
                DO 782 I=1,NCOMP

 781               IF(J .GT. NCOMP) THEN
                      GOTO 782
                   ENDIF

                   IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
                      J=J+1
                      GOTO 781 
                   ENDIF

                   IF(ID1(J) .EQ. 1) THEN
                      J = J+1
                      GOTO 781
                   ENDIF

                   ZZ=ISIGN*XY(J)

                   SCORE(I) = (2.0*WRK10(I)-1.0)*WRK7(I)+
     +                        (WRK10(I)-1.0)*WRK8(I)

                   IF(LO.EQ.0) THEN
                      WRITE(6,780) ZZ,SCORE(I)
                   ELSE
                      WRITE(60,780) ZZ,SCORE(I)
                   ENDIF
                   J=J+1
 782            CONTINUE
              ENDIF

             IF(LO.EQ.0) THEN
                WRITE(6,660) TEST
                WRITE(6,665) PROB
             ELSE
                WRITE(60,660) TEST
                WRITE(60,665) PROB
             ENDIF

          ELSE
  785        IF(LO.EQ.0) THEN
                WRITE(6,790) 
                WRITE(6,791)
             ELSE
                WRITE(60,790) 
                WRITE(60,791)
             ENDIF
  790  FORMAT(5X,'NO CENSORED OBS., PETO & PRENTICE WILCOXON TEST')
  791  FORMAT(5X, '       REDUCED TO GEHAN`S WILCOXON TEST')
C
          ENDIF

       RETURN
       END

C
C      **********************************************************************
C      ********************** SUBROUTINE UNIVAR *****************************
C      **********************************************************************
C
       SUBROUTINE UNIVAR(IBACK)
C
C      *      NDAT     : DIMENSION DECLARATION                              *
C      *                                                                    *
C      *   UNIVARIATE PROBLEMS                                              *
C      *   PARAMETERS                                                       *
C      *       MVAR    :  THE MAXIMUM NUMBER OF VARIABLES ALLOWED IN A DATA *
C      *                    SET.                                            *
C      *       NDAT    :  THE MAXIMUM NUMBER OF DATA POINTS ALLOWED IN A    *
C      *                    DATA SET.                                       *
C      *       IBIN    :  THE DIMENSION SIZE FOR BINS - USED IN THE KAPLAN- *
C      *                    MEIER ESTIMATION PROCEDURE.                     *
C      *   COMMON FOR KAPLAN-MEIER AND TWO SAMPLE TESTS :                   *
C      *       FILE    :  NAME OF DATA FILE    (9 LETTERS)                  *
C      *       TITLE   :  TITLE OF THE PROBLEM (80 LETTERS)                 *
C      *       IUNI    :  INDICATOR OF PROBLEM                              *
C      *                    IF 1 : KAPLAN-MEIER ESTIMATOR                   *
C      *                    IF 2 : TWO-SAMPLE TESTS                         *
C      *       NTOT    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
C      *       NVAR    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    * 
C      *       ITEST   :  NUMBER OF VARIABLES TO BE TESTED (<=NVAR)         *
C      *       ICOL    :  INDICATOR OF SAMPLES                              *
C      *       COLM    :  NAME OF THE SAMPLE SETS                           *
C      *       COMMAND :  NAME OF COMMAND FILE                              *
C      *       OUTPUT  :  NAME OF OUTPUT FILE                               *
C      *       IND(I,J):  INDICATOR OF CENSORING (J-TH DATA POINTS OF I-TH  *
C      *                  VARIABLE)                                         *
C      *       X(I,J)  :  DATA POINTS                                       *
C      *                                                                    *
C      *    KAPLAN-MEIER ESTIMATOR                                          *
C      *       IKM     :  INDICATOR OF PRINTOUT                             *
C      *                    IF 0 : MEAN AND ERROR ONLY                      *
C      *                    IF 1 : MEAN, ERROR, SURVIVAL DISTRIBUTION       *
C      *                                                                    *
C      *    TWO-SAMPLE TESTS                                                *
C      *       IPROG(I):  INDICATOR OF TEST (I=1,5; OR 6 FOR EXIT)          *
C      *       NGROUP  :  NUMBER OF GROUPS                                  *
C      *      IGROUP(I):  INDICATOR OF GROUP                                *
C      *       IFIRST  :  INDICATOR OF FIRST GROUP                          *
C      *       ISECON  :  INDICATOR OF SECOND GROUP                         *
C      *       IFULL   :  INDICATOR OF PRINTOUT FOR THE I-TH COMBINATION    *
C      *       GROUP(I):  NAME OF THE GROUPS                                *
C      *        LKM    :  INDICATOR OF K-M ESTIMATOR                        *
C      *        IKM    :  INDICATOR OF PRINTOUT.                            *
C      *       NOTEST  :  NUMBER OF TESTS TO BE USED                        *
C      *       ISTA(I) :  INDICATOR OF GROUPS                               *
C      *                                                                    *
C      *    WRK VARIABLES AND ARRAYS                                        *
C      *       LCOMM   :  INDICATOR OF USE OF COMMAND FILE                  *
C      *                      IF 0, READ INFORMATION FROM THE TERMINAL      *
C      *                      IF 1, READ INFORMATION FROM THE COMMAND FILE  *
C      *       CHECK   :  READER OF Y/N QUESTIONS                           *
C      *       NTOT    :  NUMBER OF DATA POINTS                             *
C      *       NCHANGE :                                                    *
C      *       ICHANGE :                                                    *
C      *       FGROUP  :                                                    *
C      *       SGROUP  :                                                    *
C      *      CHAR(I,J):  READ-IN FORM OF SEVERAL INPUTS                    *
C      *     CTEST(I,1):  READ-IN FORM OF NOTEST                            *
C      *     CPROG(I,J):  READ-IN FORM OF IPROG(J)                          *
C      *      CCOL(I,J):  READ-IN FORM OF ICOL(J)                           *
C      *      IFIRST   :                                                    *
C      *      ISECON   :                                                    *
C      *      CF(I,1)  :  READ-IN FORM OF IFIRST                            *
C      *      CS(I,1)  :  READ-IN FORM OF ISECON                            *
C      *     CIKM1(I,J):  READ-IN FORM OF IKM1(J)                           *
C      *       N1      :  NUMBER OF DATA POINTS IN THE FIRST GROUP          *
C      *       N2      :  NUMBER OF DATA POINTS IN THE SECOND GROUP         *
C      *      JIND(I,J):  IND(I,J) IN FIRST OR SECOND GROUP                 *
C      *       Z(I,J)  :  X(I,J) IN FIRST OR SECOND GROUP                   *
C      *                                                                    *
C      *   NOTE:                                       *
C      *     ALL VARIABLES WRK_, IWRK_, DWRK1, AND SWRK1 ARE USED TO COLLECT*
C      *     ALL DIMENSION DECLARATIONS TO THE TOP LEVEL SUBROUTINE.        *
C      *     THEY DO NOT DIRECTLY AFFECT ANY OTHER SUBROUTINES. IF ONE NEEDS*
C      *     TO KNOW WHAT KIND OF WORK A VARIABLE DOES, GO TO THE LOWER     *
C      *     SUBROUTINES AND READ THE DESCRIPTION OF THE VARIABLE.          *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *     DATA1, DATA2, DATAIN, KMESTM, TWOST                            *
C      *                                                                    *
C
C
C      *                                                                    *
C      *  START UNIVARIATE PROBLEM: K-M ESTIMATOR OR TWO-SAMPLE TESTS       *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       
C      *   THIS PARAMETER STATEMENT AND THE ONE IN BIVAR.F ARE ALL THAT     *
C      *   NEEDS TO BE CHANGED IF THE USER WISHES TO WORK ON DATA SETS      *
C      *   OF MORE THAN 500 OBSERVATIONS OR WITH MORE THAN FOUR VARIABLES.  *

C  **************************************************************************
       PARAMETER(MVAR=4, NDAT=500, IBIN=50)
C  **************************************************************************

       CHARACTER*1 CHECK,CHAR(4,10)
       CHARACTER*1 CF(4,1),CS(4,1),CIKM(4,1),CIS4(4,1)
       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM,GROUP(10)
       CHARACTER*80 TITLE 


       DIMENSION IND(MVAR,NDAT),X(MVAR,NDAT)
       DIMENSION JIND(MVAR,NDAT),Z(MVAR,NDAT),ISTA(NDAT)
       DIMENSION IGROUP(MVAR),JGROUP(MVAR)
C
C     *     THE DIMENSIONS BELOW WERE ADDED TO COLLECT ALL DIMENSION        *
C     *     DECLARATIONS IN THE TOP SUBROUTNE                        *
C
       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
       DIMENSION WRK12(NDAT),WRK13(NDAT),WRK14(NDAT)
       DIMENSION WRK15(NDAT),WRK16(NDAT),WRK17(NDAT)
       DIMENSION WRK18(NDAT),WRK19(NDAT),WRK20(NDAT)
       DIMENSION WRK21(NDAT),WRK22(NDAT),WRK23(NDAT)
       DIMENSION WRK24(NDAT),WRK25(NDAT),WRK26(NDAT)

       DIMENSION IWRK1(NDAT), IWRK2(NDAT), IWRK3(NDAT)
       DIMENSION BWRK1(IBIN),BWRK2(IBIN), BWRK3(IBIN)
       DIMENSION DWRK1(MVAR,NDAT), SWRK1(MVAR)
C
C
   50  FORMAT(A1)
C
 1000  PRINT *
       PRINT *,'    SELECT PROBLEM: ' 
       PRINT *,'     1 KAPLAN-MEIER DISTRIBUTION '
       PRINT *,'     2 TWO SAMPLE TESTS'
       PRINT *,'     3 EXIT '
       PRINT *
       PRINT *,'    (IF YOU CHOOSE OPTION 2, YOU CAN STILL DO 1 LATER) '
       PRINT *
 1010  WRITE(6,1020)
C
C      *         SELECT PROBLEM                                             *
C
 1020  FORMAT(' CHOICE? ')
C
       CALL DATA1(IUNI)
C
       IF((IUNI.EQ.1).OR.(IUNI.EQ.2).OR.(IUNI.EQ.3)) GOTO 1030
       PRINT *
       PRINT *,'      PLEASE TYPE ONCE MORE'
       GOTO 1010
C
 1030  IF(IUNI.EQ.3) STOP 
C
 1120  IF(IUNI.EQ.2) GOTO 1140

       PRINT *
       PRINT *,'   ***  KAPLAN-MEIER ESTIMATOR  ***'
       PRINT *
       GOTO 1330
C
C
C      *     DISPLAY THE INFORMATION ABOUT TWO SAMPLE TESTS                 *
C
 1140  PRINT *
       PRINT *,'   ***     TWO-SAMPLE TESTS     ***'
       PRINT *
       PRINT *
C
       LCOMM=1
C
C      *   CHECK WHETHER THE DATA NEEDS TO BE READ FROM A FILE              *
C
 1330  PRINT *
       PRINT *,'DO YOU WANT TO READ THE INPUTS'
       WRITE(6,1340)
 1340  FORMAT('     FROM A COMMAND FILE (Y/N)? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 2680
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 1350
       GOTO 1330
C
C      *  READ INFORMATION COMMON TO K-M ESTIMATOR AND TWO SAMPLE TESTS     *
C
 1350  PRINT *
C
C      *               READ NAME OF THE DATA FILE                           *
C
 1360  PRINT *
       WRITE(6,1370)
 1370  FORMAT('WHAT IS THE DATA FILE NAME ? ')
       READ(5,1380) FILE
 1380  FORMAT(A9)
       PRINT *
       WRITE(6,1400) FILE
 1400  FORMAT(5X,'THE FILE NAME IS ',A9)
       PRINT *
C
C      *              READ TITLE OF THE PROBLEM                             *
C
 1410  PRINT *
       WRITE(6,1420)
 1420  FORMAT('WHAT IS THE PROBLEM TITLE? ')
       READ(5,1430) TITLE 
 1430  FORMAT(A80)
C
C      *            READ THE NUMBER OF VARIABLES                            *
C
       PRINT *
       WRITE(6,1480)
 1480  FORMAT('HOW MANY VARIABLES DO YOU HAVE? ')
       CALL DATA1(NVAR)
C
C      *          CHECK WHICH VARIABLE SHOULD BE TESTED                     *
C
       ICOL=1
       IF(NVAR.EQ.1) GOTO 1840 
 1550  PRINT *
       PRINT *,'      WHICH VARIABLE DO YOU WANT TO TEST? '
 1560  WRITE(6,1570)
 1570  FORMAT(' VARIABLE NUMBER: ')
       READ(5,1580) (CHAR(I,1),I=1,4)
 1580  FORMAT(4A1)
C
C      *     CHECK IF THE CHOICE IS CORRECT                                 *
C
       CALL DATA2(CHAR,1,1,ICOL,LIND)
       IF(LIND.NE.0) PRINT *,
     +               '    PLEASE TYPE IN THE VARIABLE NUMBER AGAIN'
       IF(LIND.NE.0) GOTO 1550
       IF(ICOL.LE.NVAR) GOTO 1840
       PRINT *
       PRINT *,'   THE NUMBER IS LARGER THAN THE NUMBER OF VARIABLES'
       GOTO 1560
C
C      *                READ NAME OF THE VARIABLE                           *
C
 1840  PRINT *
       WRITE(6,1850) ICOL
 1850  FORMAT('VARIABLE ',I4,' IS NAMED')
       READ(5,1380) COLM
C
C      *   THE NEXT FEW LINES CONCERN ONLY 2-SAMPLE TESTS                   *
C      *   IF THE PROBLEM IS K-M ESTIMATION, GO TO LINE 2630                *
C
C      *          READ THE NUMBER OF GROUPS                                 *
C
       IF(IUNI.EQ.1) GOTO 2630
 2030  WRITE(6,2040)
 2040  FORMAT(/'HOW MANY GROUPS DO YOU HAVE? ')
       CALL DATA1(NGROUP)
       IF(NGROUP.LT.2) THEN
            PRINT *,'      NUMBER OF GROUPS MUST BE TWO OR MORE'
            GOTO 2030
       ENDIF
C
       IF(NGROUP.EQ.2) GOTO 2180

C
C      *  IF THE NUMBER OF GROUPS IS MORE THAN TWO, SPECIFY COMBINATIONS    *
C
 2170  PRINT *
       PRINT *,'     WHICH COMBINATION DO YOU WANT TO TEST? '
 2180  PRINT *
       WRITE(6,2190)
 2190  FORMAT('FIRST GROUP INDICATOR  ')
       CALL DATA1(IFIRST)
       IGROUP(1) = IFIRST
 2210  PRINT *
       WRITE(6,2220)
 2220  FORMAT('SECOND GROUP INDICATOR  ')
       CALL DATA1(ISECON)
 2240  IF(IFIRST.EQ.ISECON) THEN
            PRINT *,' YOU CHOSE THE SAME GROUP.'
            PRINT *,' PLEASE CHANGE THE SECOND GROUP.'
            GOTO 2210
       ENDIF
       IGROUP(2) = ISECON
C
C      *         READ THE NAME OF THE GROUPS                                *
C
 2250  PRINT *
       WRITE(6,2255) IFIRST   
 2255  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
       READ(5,1380) GROUP(1)
       PRINT *
       WRITE(6,2258) ISECON   
 2258  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
       READ(5,1380) GROUP(2)
       PRINT *
C
C      *      READ WHETHER TO PRINTOUT ONLY RESULTS OR TO GIVE FULL DETAILS *
C
 2312  WRITE(6,2314)
 2314  FORMAT('DO YOU WANT PRINTOUTS OF COMPUTATIONAL',
     +       ' DETAILS (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.NE.'Y'.AND.CHECK.NE.'y') GOTO 2316
       IFULL=1
       GOTO 2400
 2316  IF(CHECK.NE.'N'.AND.CHECK.NE.'n') GOTO 2312
       IFULL=0
C
C      *     LKM IS SET TO ONE IN THE TWO-SAMPLE TEST ROUTINE, SO THAT      *
C      *     THE KAPLAN-MEIER PERCENTILES AND MEAN FOR EACH GROUP ARE       *
C      *     AUTOMATICALLY PROVIDED.                                        *
C

 2400  LKM=1

C
C      *  CHECK WHETHER THE FULL K-M ESTIMATOR IS NEEDED.                   *
C      *  FROM THE NEXT LINE, THE INPUTS ARE COMMON FOR BOTH KAPLAN-MEIER   *
C      *  ESTIMATION AND TWO SAMPLE TESTS.                                  *
C
 2630  PRINT *
       WRITE(6,2640)
       IKM=0
 2640  FORMAT('DO YOU WANT TO PRINT OUT THE FULL K-M ',
     +        'ESTIMATE (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
          IKM=1
C
C       *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATE                 *
C
 4644    PRINT *
         PRINT *,'DO YOU WANT TO PRINT OUT '
         WRITE(6,2644)
 2644    FORMAT('           THE DIFFERENTIAL FORM (Y/N)?')
         READ(5,50) CHECK
         IF(CHECK .EQ. 'Y'.OR.CHECK.EQ.'y') THEN

           KDIFF = 1
 2646      WRITE(6,2647)
 2647      FORMAT(' SPECIFY STARTING VALUE : ')
           READ(5,2648) START
 2648      FORMAT(F10.3)

 4502      WRITE(6,4503) 
 4503      FORMAT('HOW MANY BINS DO YOU WANT? : ')
           READ(5,4554) LSTEP
 4554      FORMAT(I4)
           IF(LSTEP .LE. 0) GOTO 4502

 4506      WRITE(6,4507)
 4507      FORMAT('SPECIFY BIN SIZE :')
           READ(5,2648) BINSIZ
           IF(BINSIZ .LE. 0.0) GOTO 4506

         ELSEIF(CHECK .EQ. 'N'.OR. CHECK .EQ. 'n') THEN
           KDIFF = 0
         ELSE
           GOTO 4644
         ENDIF
       ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
         IKM=0
       ELSE
         GOTO 2630
       ENDIF
       
 2650  PRINT *
        WRITE(6,2655)
 2655  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 3230
       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 3230
       GOTO 2650
C
C      *   IF ALL INFORMATION IS TO BE READ FROM THE TERMINAL, GOTO 3230    *
C      *   FROM THE NEXT LINE, THE INPUTS ARE FROM "COMMAND" FILE           *
C
C
C      *   READ THE NAME OF "COMMAND" FILE                                  *
C
 2680  PRINT *
       WRITE(6,2690)
 2690  FORMAT('WHAT IS THE NAME OF YOUR COMMAND FILE? ')
       READ(5,1380) COMMAND
       WRITE(6,2710) COMMAND
 2710  FORMAT(5X,'YOUR COMMAND FILE IS CALLED ',A9)
C
       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
C
C      *              READ THE DATA FILE NAME                               *
C
 2820  READ(50,1380) FILE
C
C      *              READ THE TITLE OF THE PROBLEM                         *
C
       READ(50,1430) TITLE 
C
C      * READ THE NUMBER OF VARIABLES.                                      *
C
       READ(50,2830) (CHAR(I,1),I=1,4)
 2830  FORMAT(12A1)
C
C      *      CHECK THE NUMBER OF VARIABLES                                 *
C
       CALL DATA2(CHAR,1,3,NVAR,LIND)
       IF(LIND.EQ.0) GOTO 2845 
 2835  PRINT *
       PRINT *,'     NUMBER OF VARIABLES IS NOT READABLE'
       CLOSE(UNIT=50)
       STOP
 2845  IF(NVAR.LE.0) GOTO 2835
C
C      *          READ WHICH VARIABLE NEEDS TO BE TESTED                    *
C
       READ(50,2910) (CHAR(I,1),I=1,4)
 2910  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,ICOL,IN)
       IF(IN.EQ.0) GOTO 2935
 2915  PRINT *
       WRITE(6,2920) 
 2920  FORMAT(5X,'THE VARIABLE IS NOT READABLE')
       CLOSE(UNIT=50)
       STOP
 2935  IF((ICOL.LE.0).OR.(ICOL.GT.NVAR)) GOTO 2915
 2940  CONTINUE
C
C
C      *      READ THE NAME OF THE VARIABLE                                 *
C
 2962  READ(50,2963) COLM
 2963  FORMAT(10A9)
       IF(IUNI.EQ.1) GOTO 3180
C
C      *   FROM THE NEXT LINE, INPUTS ARE ONLY FOR TWO SAMPLE TESTS         *
C      *   IF IT IS A K-M ESTIMATOR PROBLEM, GO TO 3180                     *
C
C      *   READ THE NUMBER OF GROUPS                                        *
C
       READ(50,2970) (CHAR(I,1),I=1,4)
 2970  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,NGROUP,LIND)
       IF(LIND.EQ.0) GOTO 3000
 2980  PRINT *,'     IT IS NOT CLEAR HOW MANY GROUPS YOU HAVE'
       CLOSE(UNIT=50)
       STOP
 3000  IF(NGROUP.GT.1) GOTO 3005
       GOTO 2980
C
C      *     READ THE INDICATOR OF THE GROUPS                               *
C
 3005  READ(50,3010) ((CHAR(I,J),I=1,4),J=1,NGROUP)
 3010  FORMAT(60A1)
 3020  DO 3050 J=1,NGROUP
       CALL DATA2(CHAR,J,NGROUP,JGROUP(J),LIND)
       IF(LIND.EQ.0) GOTO 3050
 3025  PRINT *
       WRITE(6,3030) J
 3030  FORMAT(5X,'THE INDICATOR OF ',I4,'TH GROUP IS NOT CLEAR')
       CLOSE(UNIT=50)
       STOP
 3050  CONTINUE
C
C      *  READ NUMBER OF FIRST GROUP,SECOND GROUP,                          *
C      *       WHETHER PRINT OUT ALL OR RESULTS ONLY                        *
C      *       WHETHER K-M ESTIMATOR IS NEEDED                              *
C      *       WHETHER PRINT OUT ALL OR RESULTS ONLY FOR K-M                *
C
       READ(50,3085) (CF(I1,1),I1=1,4),(CS(I2,1),I2=1,4),
     +              (CIS4(I3,1),I3=1,4),(CIKM(I4,1),I4=1,4)
 3085  FORMAT(16A1)
C
       CALL DATA2(CF,1,1,IFIRST,LIND)
       IF(LIND.EQ.0) GOTO 3087
 3086  PRINT *
       PRINT *,'   THE VALUE FOR "IFIRST" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3087  IF((IFIRST.LT.0).OR.(IFIRST.GT.NGROUP)) GOTO 3086
C
       CALL DATA2(CS,1,1,ISECON,LIND)
       IF(LIND.EQ.0) GOTO 3089
 3088  PRINT *
       PRINT *,'   THE VALUE FOR "ISECON" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3089  IF((ISECON.LT.0).OR.(ISECON.GT.NGROUP)) GOTO 3087
       IF(ISECON.EQ.IFIRST) GOTO 3087
       IGROUP(1) = IFIRST
       IGROUP(2) = ISECON
C
       CALL DATA2(CIS4,1,1,IFULL,LIND)
       IF(LIND.EQ.0) GOTO 3091
 3090  PRINT *
       PRINT *,'    THE VALUE FOR "IFULL" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3091  IF((IFULL.EQ.0).OR.(IFULL.EQ.1)) GOTO 3092
       GOTO 3090
C
 3092  LKM = 1
C
 3095  CALL DATA2(CIKM,1,1,IKM,LIND)
       IF(LIND.EQ.0) GOTO 3097
 3096  PRINT *
       PRINT *,'    THE VALUE FOR "IKM" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3097  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 5190
       GOTO 3096
C
C      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
C
 5190  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,KDIFF,LIND)
       IF(LIND.EQ.0) GOTO 5201
 5200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
       CLOSE(UNIT=50)
       STOP
 5201  IF(KDIFF.EQ.1) GOTO 5202
       IF(KDIFF.EQ.0) GOTO 3102
       GOTO 5200
C
 5202  READ(50,4203) START
 5203  FORMAT(F10.3)
       READ(50,4204) LSTEP
 5204  FORMAT(I4)
       READ(50,4203) BINSIZ
C
 3102  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 3110
 3100  PRINT *
       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3110  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3168
       GOTO 3100
C
C      *            READ THE NAME OF THE GROUPS                             *
C
 3168  READ(50,1380) GROUP(1)
       READ(50,1380) GROUP(2)
C
C      *   READ NAME OF THE OUTPUT FILE. IF THE NAME IS BLANK, THE RESULTS  *
C      *   WILL BE SHOWN ON THE TERMINAL.                                   *
C
       READ(50,1380) OUTPUT
       IF(OUTPUT.NE.'         ') GOTO 3300
       GOTO 3230
C
C      *      FROM THE NEXT LINE, INPUTS ARE ONLY FOR THE K-M ESTIMATOR     *
C
 3180  READ(50,3200) (CHAR(I,1),I=1,4)
 3200  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,IKM,LIND)
       IF(LIND.EQ.0) GOTO 3205
 3203  PRINT *,'     IT IS NOT CLEAR WHETHER YOU WANT TO PRINT OUT ALL'
       PRINT *,'     KM ESTIMATORS'
       CLOSE(UNIT=50)
       STOP
 3205  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 4190
       GOTO 3203
C
C      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
C
 4190  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,KDIFF,LIND)
       IF(LIND.EQ.0) GOTO 4201
 4200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
       CLOSE(UNIT=50)
       STOP
 4201  IF((KDIFF.EQ.0).OR.(KDIFF.EQ.1)) GOTO 4202
       GOTO 4200
C
 4202  READ(50,4203) START
 4203  FORMAT(F10.3)
       READ(50,4204) LSTEP
 4204  FORMAT(I4)
       READ(50,4203) BINSIZ
C
C      *    INFORMATION ABOUT PRINTOUT                                       *
C
 3210  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 3220
 3215  PRINT *
       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3220  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3225
       GOTO 3215
C
C
 3225  READ(50,1380) OUTPUT
       CLOSE(UNIT=50)
       IF(OUTPUT.NE.'         ') GOTO 3300
C
C
C      *         LEAVE THE "COMMAND" FILE                                   *
C      *         CHECK OUTPUT FILE                                          *
C
 3230  OUTPUT='         '
 3240  PRINT *
       WRITE(6,3250)
 3250  FORMAT('DO YOU WANT TO SAVE THE RESULTS IN A FILE (Y/N)? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 3260
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 3300
       GOTO 3240
 3260  PRINT *
       WRITE(6,3270)
 3270  FORMAT('WHAT IS THE NAME OF THE FILE?  ')
       READ(5,1380) OUTPUT
C
C      *          READ IN DATA THOUGH THE SUBROUTINE "DATAIN"               *
C
 3300  CALL DATAIN(IUNI,FILE,NVAR,ISTA,IND,X,NTOT,NDAT,MVAR)
C
       IF(IUNI.EQ.2) GOTO 3360
C
C      *              COMPUTE THE K-M ESTIMATOR                             *
C
C
       IF(OUTPUT.NE.'         ') 
     +     OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
C  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
C  VAX/VMS MACHINES.
C     +    ,CARRIAGECONTROL='LIST'
     +     )

C
       CALL KMESTM(IND,X,NTOT,ICOL,IKM,TITLE,COLM,OUTPUT,IBIN,0,
     +             KDIFF,START,BINSIZ,LSTEP,FILE,
     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
     +             WRK9,MVAR)
C
       IF(IDATA.EQ.0) GOTO 3335
       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
       IF(OUTPUT.EQ.'         ') WRITE(6,3331)  FILE
       IF(OUTPUT.NE.'         ') WRITE(60,3332) 
       IF(OUTPUT.EQ.'         ') PRINT 3332
 3331  FORMAT(7X,' INPUT DATA FILE: ',A9)
 3332  FORMAT(5X,'   CENSORSHIP     X ')
       DO 3333 I=1,NTOT
       IF(OUTPUT.NE.'         ') WRITE(60,3334) IND(ICOL,I),X(ICOL,I)
       IF(OUTPUT.EQ.'         ') PRINT 3334,IND(ICOL,I),X(ICOL,I)
 3333  CONTINUE
 3334  FORMAT(12X,I4,F10.3)
C
C
       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)
C
 3335  PRINT *
       PRINT *,'    K-M ESTIMATOR COMPUTATIONS ARE FINISHED'
C
C      *      CHECK WHETHER THE USER WANTS TO USE OTHER METHODS            *
C
 3340  WRITE(6,3350)
 3350  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
       GOTO 3340
C
C
C      *                COMPUTE TWO SAMPLE TESTS                            *
C
 3360  IF(OUTPUT.EQ.'         ') GOTO 3370

       OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
C  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
C  VAX/VMS MACHINES.
C     +    ,CARRIAGECONTROL='LIST'
     +    )

       WRITE(60,3380)
       WRITE(60,3390)
       WRITE(60,3400) TITLE  
       WRITE(60,3390)
       GOTO 3410
 3370  WRITE(6,3380)
       PRINT *
       WRITE(6,3400) TITLE  
       PRINT *
 3380  FORMAT(8X,'   ******* TWO SAMPLE TEST ******')
 3390  FORMAT('    ')
 3400  FORMAT(8X,'TITLE : ',A80)    
C
 3410  IF(OUTPUT.EQ.'         ') GOTO 3420
       WRITE(60,3430) FILE
       WRITE(60,3435) COLM
       GOTO 3440
 3420  WRITE(6,3430) FILE
       WRITE(6,3435) COLM
 3430  FORMAT(8X,'DATA SET : ',A9)
 3435  FORMAT(8X,'VARIABLE : ',A9)
C
 3440  IF(OUTPUT.EQ.'         ') GOTO 3450
       WRITE(60,3460) GROUP(1),GROUP(2)
       GOTO 3470
 3450  WRITE(6,3460) GROUP(1),GROUP(2)
 3460  FORMAT(8X,'TESTED GROUPS ARE ',A9,' AND ',A9)
C
C
C 3470   DO 3480 M=1,NOTEST
C

 3470  CALL TWOST(X,IND,ISTA,ICOL,IGROUP(1),IGROUP(2),NTOT,
     +            IFULL,OUTPUT,M,MVAR,NDAT,FILE,
     +            WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
     +            WRK9,WRK10,WRK11,WRK12,IWRK1,IWRK2,
     +            WRK13,WRK14,WRK15,WRK16,WRK17,WRK18,WRK19,
     +            WRK20,WRK21,WRK22,WRK23,WRK24,WRK25,WRK26,IWRK3)
C
C 3480  CONTINUE

C
C      *         IF K-M ESTIMATOR IS NOT REQUESTED, GOTO 3510               *
C
       IF(LKM.EQ.0) GOTO 3510
       N1=0
       N2=0
       DO 3490 I=1,NTOT
          IF(ISTA(I).NE.IGROUP(1)) GOTO 3490
          N1=N1+1
          JIND(ICOL,N1)=IND(ICOL,I)
          Z(ICOL,N1)=X(ICOL,I)
 3490  CONTINUE
C
       CALL KMESTM(JIND,Z,N1,ICOL,IKM,TITLE,GROUP(1),OUTPUT,
     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
     +             WRK9,MVAR)
C
       DO 3500 I=1,NTOT
          IF(ISTA(I).NE.IGROUP(2)) GOTO 3500
          N2=N2+1
          JIND(ICOL,N2)=IND(ICOL,I)
          Z(ICOL,N2)=X(ICOL,I)
 3500  CONTINUE
C
       CALL KMESTM(JIND,Z,N2,ICOL,IKM,TITLE,GROUP(2),OUTPUT,
     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
     +             WRK9,MVAR)

 3510  IF(IDATA.EQ.0) GOTO 3525

       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
       IF(OUTPUT.EQ.'         ') PRINT 3331, FILE
       IF(OUTPUT.NE.'         ') WRITE(60,3520) 
       IF(OUTPUT.EQ.'         ') PRINT 3520
 3520  FORMAT(5X,'  CENSORSHIP     GROUP      X ')
       DO 3521 I=1,NTOT
       IF(OUTPUT.NE.'         ') WRITE(60,3522) IND(ICOL,I),ISTA(I),
     +                                                        X(ICOL,I)
       IF(OUTPUT.EQ.'         ') PRINT 3522,IND(ICOL,I),ISTA(I),
     +                                                        X(ICOL,I)
 3521  CONTINUE
 3522  FORMAT(12X,I4,6X,I4,F10.3)

       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)

 3525  PRINT *
       PRINT *,'     THE TWO SAMPLE TESTS ARE FINISHED'
 3530  PRINT *
       WRITE(6,3540)
 3540  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
       READ(5,50) CHECK
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
       GOTO 3530
       END

C
C      **********************************************************************
C      ******************** SUBROUTINE UNPACK *******************************
C      **********************************************************************
C
       SUBROUTINE UNPACK(X,N,LENX)
C
C      *      ALGORITHM AS 139.1 APPL.STATIST. (1979) VOL.28., NO.2         *
C      *                                                                    *
C      *      THIS SUBROUTINE EXPANDS A SYMMETRIC MATRIX STORED IN          *
C      *      LOWER TRIANGLAR FORM IN THE FIRST N*(N+1)/2 POSITIONS         *
C      *      OF X INTO A MATRIX USING THE FIRST N*N POSITIONS              *
C
C      *      LENX--THE LENGTH OF VECTOR--MUST BE NOT LESS THAN N*N         *
C      *         (I.E. MUST NOT BE LESS THAN (NVAR+1)**2                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION X(LENX)

       NSQ=N*N
       II=NSQ
       JJ=N*(N+1)/2
C
C      *                     STORE LAST ROW                                 *
C
       DO 10 I=1,N
          X(II)=X(JJ)
          II=II-1
          JJ=JJ-1
   10  CONTINUE

       DO 80 I=2,N
C
C      *      OBTAIN UPPER PART OF MATRIX FROM PART ALREADY SHIFTED         *
C
          IJ=I-1
          KK=NSQ+1-I
          DO 50 J=1,IJ
             X(II)=X(KK)
             II=II-1
             KK=KK-N
   50     CONTINUE
C
C      *      OBTAIN LOWER PART OF MATRIX FROM ORIGINAL TRIANGULAR          *
C      *      STORAGE                                                       *
C
          IJ=N-IJ
          DO 70 J=1,IJ
             X(II)=X(JJ)
             II=II-1
             JJ=JJ-1
   70     CONTINUE
   80  CONTINUE

       RETURN
       END

C     
C
C***************************************************************************
C**************************  SUBROUTINE WLCXN  *****************************
C***************************************************************************
C
C
       SUBROUTINE WLCXN(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R, D1, E1, R1, 
     +                  D2, E2, R2,SCORE,VAR)
C     *
C     * THIS SUBROUTINE COMPUTES THE GEHAN GENERALIZED WILCOXON STATISTIC   *
C     * WITH CONDITIONAL PERMUTATION VARIANCE (HYPERGEOMETRIC VARIANCE)     *
C     * FROM EQUATIONS (2.2) AND (2.3) IN LATTA, 'A MONTE-CARLO STUDY OF    *
C     * SOME TWO-SAMPLE RANK TESTS WITH CENSORED DATA', 1981, JOURNAL OF    *
C     * THE AMERICAN STATISTICAL ASSOCIATION, VOL 76, PP 713-719.           *
C     *                                                                     *
C     * INPUT                                                               *
C     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
C     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
C     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
C     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
C     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
C     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
C     *                                                                     *
C     * OUTPUT                                                              *
C     *     TEST    : STANDARDIZED GEHAN STATISTIC                          *
C     *     PROB    : PROBABILITY                                           *
C     *                                                                     *
C     * OTHERS                                                              *
C     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
C     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
C     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
C     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
C     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
C     *      R(I)   : RISK SET AT XY(I)                                     *
C     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
C     *                                                     XY(I) & XY(I+1) *
C     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
C     *      SCORE  : SCORE OF THE DATA                                     *
C     *      VAR    : VARIANCE OF THE DATA                                  *



       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)

       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
       DIMENSION D(NTOT),E(NTOT),R(NTOT)
       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
       DIMENSION E2(NTOT),R2(NTOT)
       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO


       I = 1
       L = 1
       R1(L) = REAL(N1)
       R2(L) = REAL(N2)
       R(L)  = REAL(NCOMP)
       ET1 = 0.0
       ET2 = 0.0

C
C     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
C     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
C
   10  IF(ID1(I) .NE. 0) THEN
          IF(ID2(I) .EQ. 1) THEN
             ET1 = ET1 + 1.0
          ELSE
             ET2 = ET2 + 1.0
          ENDIF
          I = I + 1
          GOTO 10
       ENDIF
C
C     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
C     *     FINISHED.                                                      *
C
   20  D(L)  = 0.0
       D1(L) = 0.0
       D2(L) = 0.0
       E(L)  = 0.0
       E1(L) = 0.0
       E2(L) = 0.0
       TEMP  = XY(I)
C
C     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
C     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
C     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
C     * APPROPRIATE GROUP.                                                *
C

  30   IF(ID1(I) .EQ. 0) THEN
          IF(ID2(I) .EQ. 1) THEN
             D1(L) = D1(L) + 1.0
          ELSE
             D2(L) = D2(L) + 1.0
          ENDIF

          D(L) = D1(L) + D2(L)

C
C     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
C     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
C
       ELSE
         IF(ID2(I) .EQ. 1) THEN
            E1(L) = E1(L) + 1.0
         ELSE
            E2(L) = E2(L) + 1.0
         ENDIF
            E(L) = E1(L) + E2(L)
         ENDIF

         IF(I .LE. NCOMP) THEN
           I = I + 1
C
C     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
C     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
C     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
C     * THE NUMBER OF THE CENSORED DATA POINTS                             *
C
           IF(TEMP .EQ. XY(I)) GOTO 30
           IF(ID1(I) .NE. 0) GOTO 30

C
C     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
C
           IF(L .EQ. 1) THEN
             R1(L) = R1(L) - ET1
             R2(L) = R2(L) - ET2
             R(L)  = R1(L) + R2(L)
          ELSE
             R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
             R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
             R(L)  = R1(L) + R2(L)
          ENDIF
          L = L + 1
          GOTO 20
        ENDIF
C
C     *       COMPUTE THE SCORE AND VARIANCE                         *
C

        SCORE = 0.0
        VAR   = 0.0
        L1 = L - 1
        DO 200 I = 1, L1

           SCORE = SCORE+R(I)*(D2(I)-(R2(I)*D(I)/R(I)))

           IF(R(I) .GT. 1.0) THEN
              VAR = VAR + D(I)*(R(I)**2.0)*(R2(I)/R(I))*
     +              (1.0-(R2(I)/R(I)))*((R(I)-D(I))/(R(I)-1.0))
           ENDIF

C           VAR = VAR+D(I)*((R(I)-REAL(I))**2)+E(I)*(REAL(I)**2)

  200   CONTINUE

C        VAR = VAR*REAL(N1*N2)/REAL(NCOMP*(NCOMP-1))

CC
C      *        NOW COMPUTE THE GEHAN STATISTIC                          *
C
        TEST = SCORE/DSQRT(VAR)
        PROB = 1.0 - AGAUSS(TEST)
 
        RETURN
        END

C
C      **********************************************************************
C      *********************  SUBROUTINE XDATA ******************************
C      **********************************************************************
C
       SUBROUTINE XDATA(X,XX,IND,IND2,IMUL,ICOL,NTOT,MVAR)
C
C      *  THIS SUBROUTINE CHANGES DATA FORMAT                               *
C      *                                                                    *
C      *  INPUT      X(I,J)   : VARIABLES                                   *
C      *            IND(I,J)  : INDICATOR OF CENSORSHIP                     *
C      *             IMUL     : NUMBER OF VARIABLES                         *
C      *             ICOL     : COLUMN OF THE VARIABLE WHICH NEEDS TO BE    *
C      *                        USED                                        *
C      *             NTOT     : NUMBER OF DATA POINTS                       *
C      *             MVAR     : DIMENSION SIZE                              *
C      *                                                                    *
C      *  OUTPUT     XX(I,J)  : VARIABLES                                   *
C      *            IND2(I)   : INDICATOR OF CENSORSHIP                     *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT),IND(MVAR,NTOT),IND2(NTOT)
C
       DO 100 I=1,NTOT

          IF(IMUL.GT.1) THEN
C
C      *     THE PROBLEM WITH MORE THAN ONE INDEPENDENT VARIABLE            *
C
             DO 20 J=1,IMUL
                XX(J,I)=X(J,I)
                IND2(I)=IND(1,I)
   20        CONTINUE
C
          ELSE
C
C      *     IF THE PROBLEM IS TWO DIMENSIONAL                              *
C
             XX(1,I)=X(ICOL,I)
             IND2(I)=IND(1,I)
          ENDIF
C
  100  CONTINUE
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE XVAR  *******************************
C      **********************************************************************
C
       SUBROUTINE XVAR(IND,X,J,NTOT,ISIGN,ZU,ZC,IU,IC,ISAVE,
     +                 ATIE,RISK,XT,ZTEMP,SWRK1,LTOT,MVAR,INDEX)
C
C      *       THIS SUBROUTINE DISTINGUISHES UNCENSORED AND CENSORED        *
C      *       DATA IN THE X VARIABLE AND SORTS IT INTO ZU AND ZC. ALSO,    *
C      *       IF THE DATA CONTAIN UPPER LIMITS, THE SIGN OF THE            *
C      *       VALUES ARE CHANGED SO THAT THE PROGRAM FOR THE LOWER         *
C      *       LIMITS CAN BE USED. ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME          *
C      *       LEARNING PUBLICATIONS (BELMONT:CA).                          *
C      *                                                                    *
C      *       INPUT      IND(J,I) : INDICATOR OF CENSORING                 *
C      *                    X(J,I) : VARIABLE                               *
C      *                   MVAR    : NUMBER OF THE VARIABLES( FOR DIM DEC.) *
C      *                    J      : J-TH DATA SETS                         *
C      *                   NTOT    : TOTAL NUMBER OF DATA POINTS            *
C      *                                                                    *
C      *       OUTPUT      ISIGN   : IF LOWER LIMIT, ISIGN = 1              *
C      *                             IF UPPER LIMIT, ISIGN = -1             *
C      *                   ZU(K)   : UNCENSORED DATA POINTS IN X(J,I)       *
C      *                   ZC(K)   : CENSORED DATA POINTS IN X(J,I)         *
C      *                    IU     : NUMBER OF UNCENSORED DATA POINTS       *
C      *                    IC     : NUMBER OF CENSORED DATA POINTS         *
C      *                   RISK(L) : RISK SET                               *
C      *                  ATIE(L)  : NUMBER OF TIED DATA POINTS             *
C      *                                                                    *
C      *       OTHER                                                        *
C      *                  ISAVE(I) : TEMPORARY SAVE OF ISIGN FOR EACH POINT *
C      *                             ALSO USED FOR TEMPORARY CENSORSHIP     *
C      *                             INDICATOR                              *
C      *                   XT(I)   : = X(J,I)                               *
C      *                                                                    *
C      *        SUBROUTINES                                                 *
C      *                   SORT1                                            *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(MVAR,NTOT),X(MVAR,NTOT),ZU(NTOT)
       DIMENSION ZC(NTOT),ISAVE(NTOT),ATIE(NTOT),RISK(NTOT)
       DIMENSION XT(NTOT),ZTEMP(MVAR,NTOT),SWRK1(MVAR)
       DIMENSION INDEX(NTOT)
C
       ISIGN=1     
       IC=0
       IU=0
C
C      *    FIND THE CENSORSHIP OF THE DATA SET. -1 FOR THE UPPER LIMIT    *
C      *    AND 1 FOR THE LOWER LIMIT                                      *
C
       DO 100 I=1,NTOT
          ISAVE(I) = 0
          IF(IND(J,I) .EQ. 0) GOTO 100
          ISIGN=IND(J,I)/IABS(IND(J,I))
          ISAVE(I) = ISIGN
  100  CONTINUE

C      * CHECK WHETHER THE UPPER AND LOWER LIMITS ARE MIXED IN ONE        *
C      * VARIABLE. IF SO, THE PROGRAM IS TERMINATED.                      *
C
       DO 110 I = 1, NTOT
          IF(ISAVE(I) .EQ. 0) GOTO 110
          IF(ISAVE(I) .NE. ISIGN) THEN
             PRINT *
             PRINT *,'YOU CANNOT HAVE BOTH UPPER AND LOWER LIMITS'
             PRINT *,'IN ONE VARIABLE AT THE SAME TIME.'
             PRINT *,'PLEASE CHECK THE DATA. THE PROGRAM IS TERMINATED.'
             PRINT *
             STOP
          ENDIF
  110  CONTINUE
C
C      *    IN CASE THE DATA HAS UPPER LIMITS IT IS MULTIPLIED BY ISIGN   *
C      *    TO MAKE THE DATA HAVE LOWER LIMITS (RIGHT CENSORSHIP).        *
C
       DO 280 L = 1, NTOT
          ATIE(L) = 0.0
          XT(L) = REAL(ISIGN)*X(J,L)
          ZTEMP(J,L) = 0.0
          ISAVE(L) = IND(J,L)
  280  CONTINUE
C
C     *     DATA POINTS ARE ARRANGED FROM SMALLEST TO LARGEST.             *
C     *     DETECTED AND CENSORED DATA POINTS ARE SEPARATED.               *
C     *     THEN RISK SETS AND TIED DATA POINTS ARE FOUND.                 *
C

       CALL SORT1(ISAVE,ZTEMP,XT,NTOT,J,INDEX,SWRK1,MVAR)

       L = 1

       DO 300 I=1,NTOT
          K=IABS(ISAVE(I))
          IF(K .EQ. 0) THEN 
              IU=IU+1
              ZU(IU)= XT(I)
              IF(IU .NE. 1) THEN
                 IF(ZU(IU) .EQ. ZU(IU-1)) THEN
                    ATIE(L) = ATIE(L) + 1.0
                    RISK(L) = REAL(NTOT - I)
                 ELSE
                    ATIE(L) = 1.0
                    RISK(L) = REAL(NTOT - I)
                    L = L + 1
                 ENDIF
              ELSE
                 ATIE(L) = 1.0
                 RISK(L) = REAL(NTOT - I)
                 L = L + 1
              ENDIF
           ELSE
              IC=IC+1
              ZC(IC)= XT(I)
           ENDIF
  300   CONTINUE
        LTOT = L - 1

        RETURN
        END
