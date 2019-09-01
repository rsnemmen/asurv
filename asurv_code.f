XC
XC      **********************************************************************
XC      *********************** SUBROUTINE AARRAY ****************************
XC      **********************************************************************
XC
X       SUBROUTINE AARRAY(Z,IND,ISTA,IS,NTOT,NDAT,MVAR,NG1,NG2,XY,
X     +            ID1,ID2,ISAVE)
XC
XC*******           ISIGN,IFULL IS ADDED ON "COMMON' STATEMENT               *
XC      *                                                                    *
XC      *     INPUT       Z(I,J)  : DATA TO BE TESTED                        *
XC      *                 IND(I,J): INDICATOR OF CENSORING                   *
XC      *                 ISTA(I) : INDICATOR OF GROUP                       *
XC      *                   IS    : IS-TH SUB-DATA SET                       *
XC      *                   NG1   : INDICATOR OF THE FIRST GROUP             *
XC      *                   NG2   : INDICATOR OF THE SECOND GROUP            *
XC      *                   NTOT  : TOTAL NUMBER OF DATA POINTS              *
XC      *                   LL    : INDICATOR OF OUTPUT FILE                 *
XC      *                  IPR    : INDICATOR FOR PRINTING                   *
XC      *                                                                    *
XC      *     OUTPUT         N    : NTOT                                     *
XC      *                    N1   : NUMBER OF DATA POINTS IN GROUP 1         *
XC      *                    N2   : NUMBER OF DATA POINTS IN GROUP 2         *
XC      *                   NCEN  : NUMBER OF CENSORED DATA POINTS           *
XC      *                   ISIGN : INDICATOR OF LOWER/UPPER LIMITS          *
XC      *                                                                    *
XC      *     PUT ALL OBS. IN ARRAY XY AND FORM ARRAYS ID1 AND ID2           *
XC      *           ID1(I)=0  : ITH OBS. IS UNCENSORED                       *
XC      *                  1  : ITH OBS. IS CENSORED                         *
XC      *           ID2(I)=J  : ITH OBS. IS FROM ITH SAMPLE, J=1,2           *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *           SORT2                                                    *
XC
XC*******     ALTHOUGH THIS SUBROUTINE HAS THE SAME NAME AS A PROGRAM FROM   *
XC*******     "STATISTICAL METHODS FOR SURVIVAL DATA ANALYSIS" BY ELISA T.   *
XC*******     LEE, 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA),        *
XC*******     IT IS DIFFERENT EXCEPT IN THE GENERAL PURPOSE.                 *
XC*******     ID1(I) IS ASSIGNED IN THE OPPOSITE WAY SO THAT THE PPROGRAM    *
XC*******     CAN USE THE DATA SETS WHICH ARE MADE FOR OTHER PROGRAMS.       *
XC
X
X        IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X        DIMENSION Z(MVAR,NDAT),IND(MVAR,NDAT),ISTA(NTOT)
X        DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT),ISAVE(NTOT)
X        COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
XC
X        IU=0
X        NCEN=0
X        NCOMP=0 
X        N1=0
X        N2=0
X        ISIGN=1
XC
XC      *    FIND THE CENSORSHIP OF THE DATA SET. -1 FOR UPPER LIMITS       *
XC      *    AND 1 FOR LOWER LIMITS                                         *
XC
X       DO 100 I=1,NTOT
X          ISAVE(I) = 0
X          IF(IND(IS,I) .EQ. 0) GOTO 100
X            ISIGN=IND(IS,I)/IABS(IND(IS,I))
X            ISAVE(I) = ISIGN
X  100  CONTINUE
X
XC      * CHECK WHETHER THE UPPER AND LOWER LIMITS ARE MIXED IN THE SAME   *
XC      * VARIABLE. IF SO, THE PROGRAM IS TERMINATED.                      *
XC      * THIS TEST WAS ADDED.                                             *
XC
X       DO 110 I = 1, NTOT
X          IF(ISAVE(I) .EQ. 0) GOTO 110
X          IF(ISAVE(I) .NE. ISIGN) THEN
X          PRINT *
X             PRINT *,'YOU CANNOT HAVE BOTH UPPER AND LOWER LIMITS'
X             PRINT *,'IN ONE VARIABLE AT THE SAME TIME.'
X             PRINT *,'PLEASE CHECK YOUR DATA.'
X             PRINT *,'THE PROGRAM HAS BEEN TERMINATED.'
X             PRINT *
X             STOP
X          ENDIF
X  110  CONTINUE
XC     
XC      *        COUNT NUMBER OF DATA POINTS IN THE TWO SUBSAMPLES        *
XC
X
X        DO 400 I = 1, NTOT
X        IF((ISTA(I) .EQ. NG1) .OR. (ISTA(I) .EQ. NG2)) THEN
X           NCOMP = NCOMP + 1
X           XY(NCOMP) = ISIGN*Z(IS,I)
X           IF(ISTA(I) .EQ. NG1) ID2(NCOMP) = 1
X           IF(ISTA(I) .EQ. NG2) ID2(NCOMP) = 2
X           IF(IABS(IND(IS,I)) .NE. 1) THEN
X              ID1(NCOMP) = 0
X              IU = IU + 1
X              IF(ID2(NCOMP) .EQ. 1) N1 = N1 + 1
X              IF(ID2(NCOMP) .EQ. 2) N2 = N2 + 1
X           ELSE
X              ID1(NCOMP) = 1
X              NCEN = NCEN + 1
X              IF(ID2(NCOMP) .EQ. 1) N1 = N1 + 1
X              IF(ID2(NCOMP) .EQ. 2) N2 = N2 + 1
X        ENDIF
X        ENDIF
X  400   CONTINUE
X
X        CALL SORT2(XY, ID1, ID2, NTOT)
X
X        RETURN
X        END
X
X
XC
XC      **********************************************************************
XC      ********************* FUNCTION AGAUSS  *******************************
XC      **********************************************************************
XC
X       FUNCTION AGAUSS(Z)
XC
XC      *        EVALUATES THE INTEGRAL OF THE GAUSSIAN PROBABILITY FUNCTION *
XC      *        OBTAINED FROM PROGRAM 3-5 ON P. 35 OF "DATA REDUCTION AND   *
XC      *        ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", P. R. BEVINGTON, *
XC      *        1969, McGRAW HILL, (NY:NY).                                 *
XC      *                                                                    *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
XC     
XC
X       Z=DABS(Z)
X       AGAUSS=1.0
XC
XC      *         IF Z>5.0, USE APPROXIMATION FOR PROB TO AVOID ERROR        *
XC
X       IF(Z.LE.5.0) THEN
X          DENOM=1.0
X          IF(Z .GT. 0.0) THEN
X             TERM=0.7071067812D00*Z
X             SUM=TERM
X             Y2=(Z**2)/2.0
X   31        DENOM=DENOM+2.0
X             TERM=TERM*(Y2*2.0/DENOM)
X             SUM=SUM+TERM
X             IF(TERM/SUM-1.0E-10 .GT. 0.0) THEN
X                GOTO 31
X             ELSE
X                AGAUSS=1.128379167D00*SUM*DEXP(-Y2)
X             ENDIF
X          ELSE
X             AGAUSS = 0.0
X          ENDIF
X       ENDIF
X
X       RETURN
X       END
X
XC
XC*************************************************************************
XC********************* SUBROUTINE AKRANK *********************************
XC*************************************************************************
XC
XC
X       SUBROUTINE AKRANK(IND, X, NTOT, IP, R, MVAR,ZU, ZC,
X     +                 PL, F, V, FMASS, ITEMP, PTEMP,Z1,
X     +                 WRK1,WRK2,WRK3,DWRK1,IWRK1,SWRK1) 
XC
XC     *   THIS SUBROUTINE COMPUTES AKRITAS' RANK                         *
XC     *                                                                  *
XC     *   REFERENCE                                                      *
XC     *         PENN STATE UNIVERSITY, DEPARTMENT OF STATISTICS,         *
XC     *         TECHNICAL REPORTS AND PREPRINTS SERIES, NUMBER 87,       *
XC     *         "ALIGNED RANK TESTS FOR REGRESSION WITH CENSORED DATA",  *
XC     *         MICHAEL G. AKRITAS, SEPTEMBER 1989                       *
XC     *   INPUT                                                          *
XC     *         IND      : INDICATOR OF CENSORSHIP                       *
XC     *         X        : VARIABLE                                      *
XC     *         NTOT     : TOTAL NUMBER OF DATA POINTS                   *
XC     *                                                                  *
XC     *   OUTPUT                                                         *
XC     *         R        : RANK                                          *
XC     *         PL       : PL ESTIMATOR                                  *
XC     *         F        : 1.0 - PL  (DISTRIBUTION FUNCTION)             *
XC     *                                                                  *
XC     *   OTHER VARIABLES                                                *
XC     *         IP       : INDEX OF VARIABLE BEING RANKED                *
XC     *         MVAR     : NUMBER OF VARIABLES                           *
XC     *         ZU       : DETECTED DATA                                 *
XC     *         ZC       : CENSORED DATA                                 *
XC     *         FMASS    : JUMPS IN PL ESTIMATOR                         *
XC     *         IU       : NUMBER OF DETECTIONS                          *
XC     *         IC       : NUMBER OF CENSORED DATA POINTS                *
XC     *         PTEMP    : TEMPORARY STORAGE OF PL ESTIMATOR             *
XC     *                                                                  *
XC     *   SUBROUTINES                                                    *
XC     *         XVAR, PLESTM, SORT1                                      *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION IND(MVAR, NTOT), X(MVAR, NTOT), ZU(NTOT), ZC(NTOT)
X       DIMENSION PL(NTOT), R(MVAR, NTOT), F(NTOT), V(NTOT), FMASS(NTOT)
X       DIMENSION Z1(MVAR, NTOT), ITEMP(NTOT), PTEMP(NTOT)
X       DIMENSION IWRK1(NTOT),DWRK1(MVAR,NTOT),SWRK1(MVAR)
X       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT)
X
XC
XC     *     CALL SUBROUTINE XVAR : DISTINGUISH DETECTIONS AND CENSORED   *
XC     *     DATA POINTS.                                                 *
XC
X
X       CALL XVAR(IND,X,IP,NTOT,ISIGN,ZU,ZC,IU,IC,IWRK1,WRK1,WRK2,WRK3,
X     +           DWRK1,SWRK1,LTOT,MVAR,ITEMP)
X
XC
XC     *     CALL PLESTM : PL ESTIMATOR COMPUTATION                       *
XC
XC
X       DO 5 I = 1, NTOT
X          ITEMP(I) = 0
X          Z1(1, I) = 0.0
X          IWRK1(I) = IND(IP,I)
X    5  CONTINUE
X    
X       IF(IU .EQ. 0) THEN
X          WRITE(6,3)
X    3     FORMAT('NO DETECTIONS: PROGRAM IS TERMINATED')
X          STOP
X       ENDIF
X       
X       CALL SORT1(IWRK1,Z1,ZU,IU,1,ITEMP,SWRK1,MVAR)
XC
X       IF(IC .NE. 0) CALL SORT1(IWRK1,Z1,ZC,IC,1,ITEMP,SWRK1,MVAR)
X
XC
XC  >>>>  Bug fixed Sept. 1996.  NCH was missing from following line  <<<<
X
X       CALL PLESTM(ZU, ZC, IU, IC, PL, V, NTOT,SMEAN,SIGM,ICH,NCH,IWRK1)
X       
XC
XC     *  IF THE DATA CONTAINS CENSORED DATA POINTS, THE PRODUCT LIMIT    *
XC     *  ESTIMATOR MUST BE ADJUSTED TO INCLUDE CENSORED DATA POINTS.     *
XC
X       IF(IC .NE. 0) THEN
X       
XC     *   IF THE DATA HAS UPPER LIMITS, FIRST THE PRODUCT LIMIT ESTIMATOR*
XC     *   MUST BE ADJUSTED.                                              *
X
X          IF(ISIGN .LT. 0) THEN
X             FMASS(1) = 1.0 - PL(1)
X             DO 10 I = 2, IU
X                FMASS(I) = PL(I-1)-PL(I)
X   10        CONTINUE
X
X             J = IU/2
X             DO 20 I = 1, J
X                FTEMP=FMASS(I)
X                FMASS(I)=FMASS(IU-I+1)
X                FMASS(IU-I+1)=FTEMP
X   20        CONTINUE
X
X             DO 40 I = 1, IU
X                PTEMP(I) = 1.0
X                DO 30 J = 1, I
X                   PTEMP(I)=PTEMP(I)-FMASS(J)
X   30           CONTINUE
X   40        CONTINUE
X          ELSE
X             DO 50 I = 1, IU
X                PTEMP(I)=PL(I)
X   50        CONTINUE
X          ENDIF
X
XC      *  NOW, PRODUCT LIMIT ESTIMATOR VALUES ARE ASSIGNED TO CENSORED  *
XC      *  DATA POINTS.                                                  *
X
X          IF(IND(IP,1) .EQ. 0) THEN
X             PL(1) = PTEMP(1)
X             J = 1
X          ELSE
X             PL(1) = 1.0
X             J = 0
X          ENDIF
X       
X          DO 60 I = 2, NTOT
X             IF(IND(IP, I) .EQ. 0) THEN
X                J = J + 1
X                PL(I) = PTEMP(J)
X             ELSE 
X                PL(I) = PL(I-1)
X             ENDIF
X   60     CONTINUE
X   
X       ENDIF
X
XC     *  THE PRODUCT LIMIT ESTIMATE IS NOW USED TO ESTIMATE THE       * 
XC     *  DISTRIBUTION FUNCTION (F) AT ALL POINTS.                     *
X
X       DO 65 I = 1, NTOT
X          F(I) = 1.0 - PL(I)
X   65  CONTINUE
X
XC
XC     *   COMPUTE HERE AKRITAS' RANK USING F-VALUES                 *
XC
X       DO 90 I = 1, NTOT
X         IF(IND(IP, I) .EQ. 0) THEN
X            R(IP, I) = REAL(NTOT)*F(I)
X         ELSEIF(IND(IP, I) .GT. 0) THEN
X            R(IP, I) = REAL(NTOT)*(0.5 + 0.5*F(I))
X         ELSE
X            R(IP, I) = NTOT*(0.5*F(I))
X         ENDIF
X   90  CONTINUE
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE ARISK  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE ARISK(R,XM,X,E1,NG,H,XY,ID1,NTOT)
XC
XC  
XC      *       THIS SUBROUTINE COMPUTES THE FOLLOWING FOUR                  *
XC      *       ARRAYS FOR SUBROUTINE COX, LRANK, AND PWLCXN.                *
XC      *         R(I) : NO. OF OBSERVATIONS IN RISK SET AT THE              *
XC      *                I-TH DISTINCT FAILURE TIME.                         *
XC      *        XM(I) : MULTIPLICITY OF THE I-TH DISTINCT                   *
XC      *                FAILURE TIME.                                       *
XC      *        E1(I) : XM(I)/R(I)                                          *
XC      *         H(I) : KAPLAN AND MEIER'S ESTIMATES OF THE                 *
XC      *                SURVIVOR FUNCTION                                   *
XC      *                                                                    *
XC      *         X(I) : THE ARRAY OF DISTINCT FAILURE TIMES                 *
XC      *         NG   : NO OF X                                             *
XC      *  THIS SUBROUTINE IS OBTAINED FROM ELISA T. LEE, "STATISTICAL       *
XC      *  METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING      *
XC      *  PUBLICATIONS (BELMONT:CA); BUT HAS BEEN SIGNIFICANTLY MODIFIED.   *
XC      *                                                                    *
XC
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION R(NTOT),XM(NTOT),X(NTOT),H(NTOT),E1(NTOT)
X       DIMENSION XY(NTOT),ID1(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
XC
X       L=1
X       I=1
X       R(L)=REAL(NCOMP)
XC
XC      *             COMPUTE RISK SETS, AND OTHER QUANTITIES                *
XC
X  24   IF(ID1(I).NE.0) THEN
X          R(L)=R(L)-1.0
X          I=I+1
X          GOTO 24
X       ENDIF
X  25   XM(L)=1.0
X       XNC=0.0
X       TEMP=XY(I)
X       X(L)=TEMP
X
X  21   IF(I.NE.NCOMP) THEN
X          I=I+1
XC
X          IF(ID1(I).NE.1) THEN
X             IF(TEMP.NE.XY(I)) GOTO 20
X             XM(L)=XM(L)+1.0
X             GOTO 21
X          ENDIF
X
X  26      XNC=XNC+1.0
X          X(L)=TEMP
X          GOTO 21
X
X  20      L=L+1
X          R(L)=R(L-1)-XM(L-1)-XNC
X          GOTO 25
X       ENDIF
X
X  23   X(L)=TEMP
X       NG=L
XC    
XC      *          COMPUTE KM ESTIMATOR                                      *
X
X       DO 30 I=1,NG
X          E1(I)=XM(I)/R(I)
X  30   CONTINUE
X
X       H(1)=1.0
X       NG1=NG+1
X
X       DO 31 I=2,NG1
X          H(I)=H(I-1)*(1.0-E1(I-1))
X  31   CONTINUE
X
X       RETURN
X       END
X
X
XC
XC
XC      *  ASURV:   SURVIVAL ANALYSIS PACKAGE FOR ASTRONOMERS                *
XC      *                                                                    *
XC      *  DEVELOPED BY:        TAKASHI ISOBE                                *
XC      *                 CENTER FOR SPACE RESEARCH                          * 
XC      *           THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY                *
XC      *                                                                    *
XC      *                      MICHAEL LAVALLEY                              *
XC      *                  DEPARTMENT OF STATISTICS                          *
XC      *               THE PENSYLVANIA STATE UNIVERSITY                     *
XC      *       330A CLASSROOM BUILDING, UNIVERSITY PARK PA 16802            *
XC      *                 INTERNET: MLV@STAT.PSU.EDU                         *
XC      *                                                                    *
XC      *                      ERIC FEIGELSON                                *
XC      *             DEPARTMENT OF ASTRONOMY AND ASTROPHYSICS               * 
XC      *                THE PENSYLVANIA STATE UNIVERSITY                    *
XC      *              525 DAVEY LAB. UNIVERSITY PARK PA 16802               *
XC      *                                                                    *
XC      *  REV. 1.2  SECOND UPDATE   SUMMER 1992                             *
XC      *                                                                    *
XC      *     THIS PACKAGE IS WRITTEN TO PROVIDE SEVERAL                     *
XC      *  SURVIVAL ANALYSIS METHODS WHICH ARE USEFUL IN ANALYZING           *
XC      *  ASTRONOMICAL DATA. SURVIVAL ANALYSIS IS A GROUP OF STATISTICAL    *
XC      *  METHODS WHICH TREAT PROBLEMS WITH CENSORED DATA (UPPER OR LOWER   *
XC      *  LIMITS). THIS PACKAGE INCLUDES SOME TECHNIQUES DEVELOPED IN       *
XC      *  IN OTHER FIELDS (E.G. ECONOMICS, ACTUARIAL SCIENCE, RELIABILITY   *
XC      *  MATHEMATICS), AND A FEW METHODS DEVELOPED BY ASTRONOMERS.         *
XC      *                                                                    *
XC      *   THE METHODS PROVIDED IN THIS PACKAGE ARE :                       *
XC      *                                                                    *
XC      *   UNIVARIATE DISTRIBUTION :  KAPLAN-MEIER ESTIMATOR                *
XC      *   TWO-SAMPLE TESTS        :  GEHAN TEST                            *
XC      *                              LOGRANK TEST                          *
XC      *                              PETO AND PETO TEST                    *
XC      *                              PETO AND PRENTICE TEST                *
XC      *   CORRELATION TESTS       :  COX PROPORTIONAL HAZARDS MODEL        *
XC      *                              GENERALIZED KENDALL'S TAU (BHK METHOD)*
XC      *                              GENERALIZED SPEARMAN'S RHO            *
XC      *                                   (AKRITAS' METHOD)                *
XC      *   LINEAR REGRESSIONS      :  EM ALGORITHM WITH NORMAL DISTRIBUTION *
XC      *                              BUCKLEY-JAMES METHOD                  *
XC      *                              TWO-DIMENSIONAL KAPLAN-MEIER          *
XC      *                                  REGRESSION FOR DUAL-CENSORED DATA *
XC      *                                                                    *
XC      *                                                                    *
XC      *   INPUTS                                                           *
XC      *                                                                    *
XC      *       IS0     :  IF 1 : UNIVARIATE PROBLEM                         *
XC      *                     2 : CORRELATION/REGRESSION PROBLEM             *
XC      *                     3 : EXIT                                       *
XC      *                                                                    *
XC      *   SUBROUTINES DATA1, UNIVAR, BIVAR                                 *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*1 BLANK
XC       OPEN(6,CARRIAGECONTROL='LIST',STATUS='OLD')
XC
XC
XC
X       PRINT *
X       PRINT *,'    ***************************************************'
X       PRINT *,'    *                                                 *'
X       PRINT *,'    *                WELCOME TO ASURV                 *'
X       PRINT *,'    *           SURVIVAL ANALYSIS PACKAGE             *' 
X       PRINT *,'    *                FOR ASTRONOMERS                  *'
X       PRINT *,'    *                                                 *'
X       PRINT *,'    *                  DEVELOPED BY:                  *'
X       PRINT *,'    *                  TAKASHI ISOBE                  *' 
X       PRINT *,'    *         (CENTER FOR SPACE RESEARCH, MIT)        *'
X       PRINT *,'    *                 MICHAEL LAVALLEY                *'
X       PRINT *,'    *         (DEPT. OF STATISTICS, PENN STATE)       *'
X       PRINT *,'    *                  ERIC FEIGELSON                 *'
X       PRINT *,'    * (DEPT. OF ASTRONOMY & ASTROPHYSICS, PENN STATE) *'
X       PRINT *,'    *                                                 *'
X       PRINT *,'    *                                                 *'
X       PRINT *,'    *              REV 1.2  SUMMER 1992               *'
X       PRINT *,'    ***************************************************'
X       PRINT *
X       PRINT *
X       PRINT *
X       PRINT *
X       PRINT *,'              (CARRIAGE RETURN TO CONTINUE) '
X       READ(5,50) BLANK
X   50  FORMAT(A1)
X       PRINT *
XC
XC      *           START CONVERSATION WITH THE USER                         *
XC
X       PRINT *
X       PRINT *
X       PRINT *
X  100  PRINT *,'                          MENU  '
X       PRINT *
X       PRINT *
X       PRINT *,'       UNIVARIATE DATA           BIVARIATE DATA '
X       PRINT *
X       PRINT *
X       PRINT *,'     DISTRIBUTION FUNCTION       CORRELATION '
X       PRINT *,'   1 KAPLAN-MEIER ESTIMATOR    1 COX REGRESSION '
X       PRINT *,'                               2 GEN. KENDALL TAU'
X       PRINT *,'                               3 GEN. SPEARMAN RHO'
X       PRINT *
X       PRINT *
X       PRINT *,'     TWO-SAMPLE TESTS            LINEAR REGRESSION '
X       PRINT *,'   1 GEHAN TESTS               1 EM ALGORITHM WITH  '
X       PRINT *,'   2 LOGRANK TEST                 GAUSSIAN RESIDUALS ' 
X       PRINT *,'   3 PETO AND PETO TEST        2 BUCKLEY-JAMES METHOD ' 
X       PRINT *,'   4 PETO AND PRENTICE TEST       WITH KM RESIDUALS '
X       PRINT *,'                               3 SCHMITT METHOD FOR '
X       PRINT *,'                                  DUAL CENSORED DATA ' 
X       PRINT *
X       PRINT *
X       PRINT *
X       PRINT *,'            (CARRIAGE RETURN TO CONTINUE) '
X       READ(5,50) BLANK
XC
X       PRINT *
XC
XC      *  CHOICE : UNIVARIATE PROBLEM OR CORRELATION/REGRESSION PROBLEM     *
XC
X       PRINT *
X       PRINT *,'    SELECT DATA TYPE: ' 
X       PRINT *,'     1 UNIVARIATE DATA '
X       PRINT *,'     2 BIVARIATE DATA ' 
X       PRINT *,'     3 EXIT '
X  200  WRITE(6,210)
X  210  FORMAT(' CHOICE ? ')
XC  210  FORMAT('          CHOICE ? ',$)
XC
X       CALL DATA1(IS0)
XC
X       IF((IS0.EQ.1).OR.(IS0.EQ.2).OR.(IS0.EQ.3)) GOTO 300
X       PRINT *,'PLEASE TYPE ONCE MORE'
X       GOTO 200
XC
X  300  IBACK=0
X       IF(IS0.EQ.1) CALL UNIVAR(IBACK)
X       IF(IS0.EQ.2) CALL BIVAR(IBACK)
X       IF(IS0.EQ.3) STOP
XC
X       IF(IBACK.EQ.1) GOTO 100
X       STOP
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE BHK  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE BHK(IND,XX,YY,NTOT,OUTPUT,X,Y,IAA,IBB,IP,MVAR)
XC
XC      *             GENERALIZED KENDALL'S TAU CORRELATION COEFFICIENT      *
XC      *                         FOR CENSORED DATA                          *
XC
XC      *      THIS PROGRAM COMPUTES KENDALL'S TAU FOR BIVARIATE DATA        *
XC      *      SETS. THE DATA SETS CAN CONTAIN CENSORED POINTS IN THE        *
XC      *      INDEPENDENT VARIABLE AND/OR THE DEPENDENT VARIABLE.           *
XC      *      ALTHOUGH THIS PROGRAM GIVES ANSWERS FOR DATA SETS WHICH       *
XC      *      CONTAIN TIES, IT MAY NOT BE ACCURATE.                         *
XC      *    PARAMETERS :                                                    *
XC      *     INPUT                                                          *
XC      *     NTOT          : NUMBER OF OBSERVATIONS                         *
XC      *     XX(1,I)       : INDEPENDENT PARAMETER OF I-TH OBSERVATION      *
XC      *     YY(I)         : DEPENDENT PARAMETER OF I-TH OBSERVATION        *
XC      *     IND(I)        : INDICATOR OF CENSORED STATUS                   *
XC      *      EACH POINT MUST BE SPECIFIED ITS CENSORED STATUS :            *
XC      *        FOR THE LOWER LIMITS                                        *
XC      *                0   :   DETECTED POINT                              *
XC      *                1   :   ONLY DEPENDENT VARIABLE IS LOWER LIMIT      *
XC      *                2   :   ONLY INDEPENDENT VARIABLE IS LOWER LIMIT    *
XC      *                3   :   BOTH VARIABLES ARE LOWER LIMIT              *
XC      *                4   :   INDEPENDENT VARIABLE IS LOWER LIMIT AND     *
XC      *                        DEPENDENT VARIABLE IS UPPER LIMIT           *
XC      *      FOR THE UPPER LIMITS, CHANGE THE SIGN OF ABOVE INDICATORS.    *
XC      *                                                                    *
XC      *     WORK                                                           *
XC      *     X(I)           : =XX(1,I)                                      *
XC      *     Y(I)           : =YY(I)                                        *
XC      *     IP(I)          : =IND(I)                                       *
XC      *     IAA(I)         : CONCORDANCE INFORMATION FOR X                 *
XC      *     IBB(I)         : CONCORDANCE INFORMATION FOR Y                 *
XC      *     OUTPUT                                                         *
XC      *     PROB          : SIGNIFICANCE LEVEL FOR THE HYPOTHESIS THAT     *
XC      *                     X AND Y ARE NOT CORRELATED UNDER THE           *
XC      *                     GAUSSIAN DISTRIBUTION                          *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *               CENS, COEFF                                          *
XC      *                                                                    *
XC      *      REF. BROWN, HOLLANDER, AND KORWAR 1974, IN RELIABILITY        *
XC      *           AND BIOMETRY P.327, EQNS 1 TO 8, PROSCHAN AND            *
XC      *           SERFLING EDS (SIAM)                                      *
XC
XC      *     NOTE:  THIS PROGRAM IS QUITE CPU INTENSIVE FOR LARGE DATA      *
XC      *            SETS (MORE THAN A FEW HUNDRED POINTS).                  *
XC      *                                                                    *
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION XX(MVAR,NTOT),YY(NTOT),IND(NTOT),X(NTOT),Y(NTOT)
X       DIMENSION IP(NTOT),IAA(NTOT),IBB(NTOT)
X       CHARACTER*9 OUTPUT 
XC
X       SIS =0.0
X       ASUM =0.0
X       BSUM =0.0
X       AASUM=0.0
X       BBSUM=0.0
XC
XC      *       SUBSTITUE XX AND YY TO X AND Y SO THAT THE ORIGINAL VALUES   *
XC      *       WON'T BE CHANGED.                                            *
XC
X       DO 90 I=1,NTOT
X          X(I)  = XX(1,I)
X          Y(I)  = YY(I)
X          IP(I) = IND(I)
X   90  CONTINUE
XC
XC
XC      *      THE SUBROUTINE CENS ADDS OR SUBTRACTS A SMALL NUMBER          *
XC      *      FROM EACH CENSORED POINT SO THAT NO TIES WITH DETECTED        *
XC      *      POINTS OCCUR.                                                 *
XC
XC
X       CALL CENS(X,Y,IP,NTOT)
XC
XC
XC      *      START MAKING INFORMATION FOR CONCORDANCE                      *
XC
XC
X       DO 1900 I=1,NTOT
XC
XC      *      INFORMATION OF CONCORDANCE  FOR THE INDEPENDENT VAR.          *
XC
X          IA=2
X          IB=3
X          IC=4
X          ID=-2
X          IE=-3
X          IG=-4
X          IH=1
X          IJ=-1
XC
XC      *       SUBROUTINE WHICH FINDS CONCORDANCE INFORMATION               *
XC
X          CALL COEFF(I,X,IP,NTOT,IAA,IA,IB,IC,ID,IE,IG,IH,IJ)
XC
XC      *       INFORMATION OF CONCORDANCE FOR THE DEPENDENT VAR.            *
XC
X          IA=1
X          IB=3
X          IC=-4
X          ID=-1
X          IE=-3
X          IG=4
X          IH=2
X          IJ=-2
X
X          CALL COEFF(I,Y,IP,NTOT,IBB,IA,IB,IC,ID,IE,IG,IH,IJ)
XC
XC      *        START COMPUTING QUANTITIES IS, IASUM, IBSUM,                *
XC      *        IAASUM, AND IBBSUM.                                         * 
XC
X          DO 1800 J=1,NTOT
X             IF((IAA(J).EQ.0).AND.(IBB(J).EQ.0)) GOTO 1800
X             SIS=SIS+IAA(J)*IBB(J)
X             ASUM=ASUM+IAA(J)**2
X             BSUM=BSUM+IBB(J)**2
X
X 1650        DO 1700 K=1,NTOT
X                IF(IAA(J).NE.0) THEN
X                   IF(IAA(K).NE.0) THEN
X                      AASUM=AASUM+IAA(J)*IAA(K)
X                   ENDIF
X                ENDIF
X 1670           IF(IBB(J).NE.0) THEN
X                   IF(IBB(K).NE.0) THEN
X                      BBSUM=BBSUM+IBB(J)*IBB(K)
X                   ENDIF
X                ENDIF
X 1700        CONTINUE
X 1800     CONTINUE
X 1900  CONTINUE
XC
XC      *    NOW COMPUTE THE STATISTIC AND THE PROBABILITY                   *
XC 
X       D1=REAL(NTOT*(NTOT-1))
X       D2=REAL(D1*(NTOT-2))
X       ALP=2.0*(ASUM*BSUM)/D1
X       GAM=4.0*((AASUM-ASUM)*(BBSUM-BSUM))/D2
X       VAR=ALP+GAM
X       SIGMA=DSQRT(VAR)
X       Z=SIS/SIGMA
X       PROB=1.0-AGAUSS(Z)
XC
X       IF(OUTPUT.EQ.'         ') THEN
X          WRITE(6,2030) 
X          WRITE(6,2003)
X          WRITE(6,2030)
X          WRITE(6,2005) Z
X          WRITE(6,2007) PROB
X          WRITE(6,2030) 
X       ELSE
X          WRITE(60,2030) 
X          WRITE(60,2003)
X          WRITE(60,2030)
X          WRITE(60,2005) Z
X          WRITE(60,2007) PROB
X          WRITE(60,2030) 
X       ENDIF
X 2003  FORMAT(5X,'CORRELATION TEST BY GENERALIZED KENDALL`S TAU')
X 2005  FORMAT(7X,'Z-VALUE      =',F12.3)
X 2007  FORMAT(7X,'PROBABILITY  =',F13.4,/,
X     +    ' (PROBABILITY THAT A CORRELATION IS NOT PRESENT)')
X 2030  FORMAT('      ')
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE BIN  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE BIN(NTOT,MX,MY,ISKIP,ICENS,DELX,DELY,XORG,YORG,MM,
X     +                M1,M2,M3,M4,M5,M6,M7,M8,INDEX,LP,XT,YT,Z,SWRK1,
X     +                X,Y,NP,XB,YB,F,N,N1,N2,N3,N4,N5,N6,N7,N8,IB,MVAR)
X
XC
XC
XC      *                                                                    *
XC      *    THIS SUBROUTINE DOES BINNING AND CHANGES CENSORED POINTS        *
XC      *    WHICH DO NOT HAVE DETECTED POINTS ABOVE (OR BELOW)              *
XC      *    TO DETECTED POINTS.                                             *
XC      *                                                                    *
XC      *             WARNING   WARNING   WARNING   WARNING                  *
XC      *                                                                    *
XC      *    THE USER SHOULD BE WARNED THAT THIS SUBROUTINE ACTUALLY         *
XC      *    CHANGES THE DATA!!  FIRST, IT REDEFINES SOME LIMITS TO          *
XC      *    DETECTIONS.  IF THE BINS ARE CHOSEN TO BE TOO NARROW, THEN      *
XC      *    VIRTUALLY ALL LIMITS COULD BE CHANGED.  SECOND, IT PUSHES       *
XC      *    EACH LIMIT INTO THE ADJACENT BIN.  IF THE BINS ARE CHOSEN TO    *
XC      *    TO BE TOO WIDE, THIS SUBSTANTIALLY ALTERS THE MEASURED VALUES.  *
XC      *    THUS, THE USER MUST TREAD A FINE LINE IN CHOSING BIN SIZES.     *
XC      *                                                                    * 
XC      *                                                                    *
XC      *    INPUT                                                           *
XC      *           X(I)  : INDEPENDENT VARIABLE                             *
XC      *           Y(I)  : DEPENDENT VARIABLE                               *
XC      *          NP(I)  : INDICATOR OF CENSORING                           *
XC      *          NTOT   : TOTAL NUMBER OF DATA                             *
XC      *            MX   : NUMBER OF BINS IN X                              *
XC      *            MY   : NUMBER OF BINS IN Y                              *
XC      *          ISKIP  : INDICATOR OF BINNING PROCESS                     *
XC      *          ICENS  : CENSORING STATUS OF THE DATA SET                 *
XC      *      IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :                *
XC      *          DELX   : BIN SIZE OF X AXIS                               *
XC      *          DELY   : BIN SIZE OF Y AXIS                               *
XC      *          XORIG  : ORIGIN OF X                                      *
XC      *          YORIG  : ORIGIN OF Y                                      *
XC      *                                                                    *
XC      *     WORK                                                           *
XC      *           YT(I) : COPY OF Y(I) FOR SORTING PROGRAM.                *
XC      *            M1   : # OF Y LOWER LIMITS CHANGED TO DETECTIONS        *
XC      *            M2   : # OF X LOWER LIMITS CHANGED TO DETECTIONS        *
XC      *            M3   : # OF DOUBLE LOWER LIMITS CHANGED TO              *
XC      *                   DETECTIONS                                       *
XC      *            M4   : # OF Y LOWER , X UPPER LIMITS CHANGED TO         *
XC      *                   DETECTIONS                                       *
XC      *            M5   : # OF Y UPPER LIMITS CHANGED TO DETECTIONS        *
XC      *            M6   : # OF X LOWER LIMITS CHANGED TO DETECTIONS        *
XC      *            M7   : # OF DOUBLE UPPER LIMITS CHANGED TO              *
XC      *                   DETECTIONS                                       *
XC      *            M8   : # OF Y UPPER , X LOWER LIMITS CHANGED TO         *
XC      *                   DETECTIONS                                       *
XC      *            NC1, NC2,...,NC8 : # OF CENSORED POINTS. SEE THE        *
XC      *                  MAIN PROGRAM FOR THE DEFINITIONS                  *
XC      *            IB   : DIMENSION SIZE OF BINS                           *
XC      *                                                                    *
XC      *    OUTPUT                                                          *
XC      *           F(I,J): INITIAL GUESS OF THE PROBABILITY OF THE          *
XC      *                   BIN(I,J)                                         *
XC      *           N(I,J): NUMBER OF DETECTED POINTS IN THE BIN(I,J)        *
XC      *          N1(I,J): NUMBER OF Y LOWER LIMITS IN THE BIN(I,J)         *
XC      *          N2(I,J): NUMBER OF X LOWER LIMITS IN THE BIN(I,J)         *
XC      *          N3(I,J): NUMBER OF DOUBLE LOWER LIMITS IN THE BIN(I,J)    *
XC      *          N4(I,J): NUMBER OF Y LOWER, X UPPER LIMITS IN THE         *
XC      *                   BIN(I,J)                                         *
XC      *          N5(I,J): NUMBER OF Y UPPER LIMITS IN THE BIN(I,J)         *
XC      *          N6(I,J): NUMBER OF X UPPER LIMITS IN THE BIN(I,J)         *
XC      *          N7(I,J): NUMBER OF DOUBLE UPPER LIMITS IN THE BIN(I,J)    *
XC      *          N8(I,J): NUMBER OF Y UPPER, X LOWER LIMITS IN THE         *
XC      *                   BIN(I,J)                                         *
XC      *           XB(I) : COORDINATE OF CENTER OF THE BIN IN X             *
XC      *           YB(I) : COORDINATE OF CENTER OF THE BINS IN Y            *
XC      *     IF ISKIP=0, THE NEXT VALUES ARE OUTPUTS  :                     *
XC      *          DELX   : BIN SIZE OF X AXIS                               *
XC      *          DELY   : BIN SIZE OF Y AXIS                               *
XC      *          XORIG  : ORIGIN OF X                                      *
XC      *          YORIG  : ORIGIN OF Y                                      *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *                   SORT1                                            *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION INDEX(NTOT),LP(NTOT),XT(NTOT),YT(NTOT),Z(MVAR,NTOT)
X       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB),SWRK1(MVAR)
X       DIMENSION F(IB,IB),N(IB,IB),N1(IB,IB),N2(IB,IB),N3(IB,IB)
X       DIMENSION N4(IB,IB),N5(IB,IB),N6(IB,IB),N7(IB,IB),N8(IB,IB)
X       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
XC
XC
XC      *  SUBSTITUE NP, X, AND Y TO LP, XT, AND YT SO THAT THE ORIGINAL DATA*
XC      *  WON'T BE CHANGED.                                                 *
XC
X       DO 100 J=1,NTOT
X          LP(J)=NP(J)
X          XT(J)=X(J)
X          YT(J)=Y(J)
X          Z(1,J)=1.0
X  100  CONTINUE
XC
XC      *    CALL THE SUBROUTINE SORT1, AND FIND MIN. AND MAX. OF X AND Y.   *
XC      *    IF ISKIP=0, THE ORIGIN AND BIN SIZES ARE ALREADY GIVEN          *
XC
X       IF(ISKIP.EQ.0) THEN
XC
XC      *                 SORTING X                                          *
XC
X          CALL SORT1(LP,Z,XT,NTOT,1,INDEX,SWRK1,MVAR)
XC
XC      *                 SORTING Y                                          *
XC
X          CALL SORT1(LP,Z,YT,NTOT,1,INDEX,SWRK1,MVAR)
XC
XC      *              FIND THE SIZES OF BINS                                *
XC
X          DELX=XT(NTOT)-XT(1)
X          DELY=YT(NTOT)-YT(1)
X          DELX=DELX/FLOAT(MX-2)
X          DELY=DELY/FLOAT(MY-2)
XC
XC      *             FIND THE ORIGIN OF THE GRID                            *
XC
X          XORG=XT(1)-1.5*DELX
X          YORG=YT(1)-1.5*DELY
X       ENDIF
XC
XC
XC      *           INITIALIZE N, N1,....,N8, AND F                          *
XC
X       DO 300 I=1,MX
X          DO 200 J=1,MY
X             N(I,J) =0
X             N1(I,J)=0
X             N2(I,J)=0
X             N3(I,J)=0
X             N4(I,J)=0
X             N5(I,J)=0
X             N6(I,J)=0
X             N7(I,J)=0
X             N8(I,J)=0
X             F(I,J)=0.0
X  200     CONTINUE
X  300  CONTINUE
XC
X       DO 390 I=1,NTOT
XC
XC      *    FIND POSITION OF I-TH DATA POINT IN THE GRID AND COUNT          *
XC      *    NUMBERS OF N,N1,N2,.....,N8.                                    *
XC
X          IP=INT((X(I)-XORG)/DELX)+1
X          JP=INT((Y(I)-YORG)/DELY)+1
X
XC
XC      *      FOR CONVENIENCE CENSORED POINTS ARE ASSIGNED TO THE NEXT BIN  *
XC
XC
XC      *              DETECTIONS                                            *
XC
X          IF(NP(I).EQ.0) THEN
X             N(IP,JP)=N(IP,JP)+1
X      
XC
XC      *              Y LOWER LIMITS                                        *
XC
X          ELSEIF(NP(I).EQ.1) THEN
X             N1(IP,JP+1)=N1(IP,JP+1)+1
X          
XC
XC      *              X LOWER LIMITS                                        *
XC
X          ELSEIF(NP(I).EQ.2) THEN
X             N2(IP+1,JP)=N2(IP+1,JP)+1
X          
XC
XC      *              DOUBLE LOWER LIMITS                                   *
XC
X          ELSEIF(NP(I).EQ.3) THEN
X             N3(IP+1,JP+1)=N3(IP+1,JP+1)+1
X          
XC
XC      *              Y LOWER LIMITS, X UPPER LIMITS                        *
XC
X          ELSEIF(NP(I).EQ.4) THEN
X             N4(IP+1,JP-1)=N4(IP+1,JP-1)+1
X          
XC
XC      *              Y UPPER LIMITS                                        *
XC
X          ELSEIF(NP(I).EQ.-1) THEN
X             N5(IP,JP-1)=N5(IP,JP-1)+1
X          
XC
XC      *              X UPPER LIMITS                                        *
XC
X          ELSEIF(NP(I).EQ.-2) THEN
X             N6(IP-1,JP)=N6(IP-1,JP)+1
X          
XC
XC      *              DOUBLE  UPPER LIMITS                                  *
XC
X          ELSEIF(NP(I).EQ.-3) THEN
X             N7(IP-1,JP-1)=N7(IP-1,JP-1)+1
X          
XC
XC      *              Y UPPER LIMITS, X LOWER LIMITS                        *
XC
X          ELSEIF(NP(I).EQ.-4) THEN
X             N8(IP-1,JP+1)=N8(IP-1,JP+1)+1
X
X       ELSE
X          PRINT *,' THE CENSORSHIP INDICATOR IS NOT RECOGNIZED'
X          RETURN
X       ENDIF
X  390  CONTINUE
XC
XC      *    SET THE COORDINATES OF THE EACH BIN                             *
XC
X       DO 410 I=1,MX
X          XB(I)=XORG+DELX/2.0+DELX*(I-1)
X  410  CONTINUE
X
X       DO 420 I=1,MY
X          YB(I)=YORG+DELY/2.0+DELY*(I-1)
X  420  CONTINUE
XC
XC      *    START CHECKING THE RELATION BETWEEN CENSORED POINTS AND         *
XC      *    DETECTED POINTS. IF THE CENSORED POINTS ARE LOCATED  SO         *
XC      *    THAT THEY CANNOT GIVE WEIGHT TO DETECTED POINTS, THE            *
XC      *    CENSORED POINTS ARE CHANGED TO DETECTIONS.                      *
XC
X       M1=0
X       M2=0
X       M3=0
X       M4=0
X       M5=0
X       M6=0
X       M7=0
X       M8=0
XC
XC
XC      *              Y LOWER LIMITS                                        *
XC
X
X       IF(NC1.NE.0) THEN
X          DO 600 I=1,MX
X             DO 500 J=1,MY
X                JJ=MY-J+1
X                IF(N1(I,JJ).NE.0) THEN
X                   K=JJ
X  450              IF(N(I,K).EQ.0) THEN
X                      K=K+1
X                      IF(K.LE.MY) GOTO 450
X                      M1=M1+N1(I,JJ)
X                      N(I,JJ)=N(I,JJ)+N1(I,JJ)
X                      N1(I,JJ)=0
X                   ENDIF
X                ENDIF
X  500        CONTINUE
X  600     CONTINUE
X       ENDIF
X
XC
XC
XC      *              X LOWER LIMITS                                        *
XC
X       IF(NC2.NE.0) THEN
X          DO 800 J=1,MY
X             DO 700 I=1,MX
X                II=MX-I+1
X                IF(N2(II,J).NE.0) THEN
X                   L=II
X  650              IF(N(L,J).EQ.0) THEN
X                      L=L+1
X                      IF(L.LE.MX) GOTO 650
X                      M2=M2+N2(II,J)
X                      N(II,J)=N(II,J)+N2(II,J)
X                      N2(II,J)=0
X                   ENDIF
X                ENDIF
X  700        CONTINUE
X  800     CONTINUE
X       ENDIF
X
XC
XC
XC      *              DOUBLE LOWER LIMITS                                   *
XC
X       IF(NC3.NE.0) THEN
X          DO 1000 I=1,MX
X             II=MX-I+1
X             DO 950  J=1,MY
X                JJ=MY-J+1
X                IF(N3(II,JJ).NE.0) THEN
X                   L=II
X  850              K=JJ
X  900              IF(N(II,JJ).EQ.0) THEN
X                      K=K+1
X                      IF(K.LE.MY) GOTO 900
X                      L=L+1
X                      IF(L.LE.MX) GOTO 850
X                      M3=M3+N3(II,JJ)
X                      N(II,JJ)=N(II,JJ)+N3(II,JJ)
X                      N3(II,JJ)=0
X                   ENDIF
X                ENDIF
X  950        CONTINUE
X 1000     CONTINUE
X       ENDIF
X
XC
XC
XC      *              Y LOWER LIMITS, X UPPER LIMITS                        *
XC
X       IF(NC4.NE.0) THEN
X          DO 1300 I=1,MX
X             II=MX-I+1
X             DO 1200 J=1,MY
X                IF(N4(II,J).NE.0) THEN
X                   L=II
X 1050              K=J
X 1100              IF(N(L,K).EQ.0) THEN
X                      K=K-1
X                      IF(K.GE.1) GOTO 1100
X                      L=L+1
X                      IF(L.LE.MX) GOTO 1050
X                      M4=M4+N4(II,J)
X                      N(II,J)=N(II,J)+N4(II,J)
X                      N4(II,J)=0
X                   ENDIF
X                ENDIF
X 1200        CONTINUE
X 1300     CONTINUE
X       ENDIF
X
XC
XC
XC      *              Y UPPER LIMITS                                        *
XC
X       IF(NC5.NE.0) THEN
X          DO 1600 I=1,MX
X             DO 1500 J=1,MY
X                IF(N5(I,J).NE.0) THEN
X                   K=J
X 1450              IF(N(I,K).EQ.0) THEN
X                      K=K-1
X                      IF(K.GE.1) GOTO 1450
X                      M5=M5+N5(I,J)
X                      N(I,J) = N(I,J) + N5(I,J)
X                      N5(I,J)=0
X                   ENDIF
X                ENDIF
X 1500        CONTINUE
X 1600     CONTINUE
X       ENDIF
X
XC
XC
XC      *              X UPPER LIMITS                                        *
XC
X       IF(NC6.NE.0) THEN
X          DO 1800 J=1,MY
X             DO 1700 I=1,MX
X                IF(N6(I,J).NE.0) THEN
X                   L=I
X 1650              IF(N(L,J).EQ.0) THEN
X                      L=L-1
X                      IF(L.GE.1) GOTO 1650
X                      M6=M6+N6(I,J)
X                      N(I,J)=N(I,J)+N6(I,J)
X                      N6(I,J)=0
X                   ENDIF
X                ENDIF
X 1700        CONTINUE
X 1800     CONTINUE
X       ENDIF
X
XC
XC
XC      *              DOUBLE UPPER LIMITS                                   *
XC
X       IF(NC7.NE.0) THEN
X          DO 2000 I=1,MX
X             DO 1950 J=1,MY
X                IF(N7(I,J).NE.0) THEN
X                   L=I
X 1850              K=J
X 1900              IF(N(L,K).EQ.0) THEN
X                      K=K-1
X                      IF(K.GE.1) GOTO 1900
X                      L=L-1
X                      IF(L.GE.1) GOTO 1850
X                      M7=M7+N7(I,J)
X                      N(I,J)=N(I,J)+N7(I,J)
X                      N7(I,J)=0
X                   ENDIF
X                ENDIF
X 1950        CONTINUE
X 2000     CONTINUE
X       ENDIF
X
XC
XC
XC      *              Y UPPER LIMITS, X LOWER LIMITS                        *
XC
X       IF(NC8.NE.0) THEN
X          DO 2300 I=1,MX
X             DO 2200 J=1,MY
X                JJ=MY-J+1
X                IF(N8(I,JJ).NE.0) THEN
X                   L=I
X 2050              K=JJ
X 2100              IF(N(L,K).EQ.0) THEN
X                      K=K+1
X                      IF(K.LE.MY) GOTO 2100
X                      L=L-1
X                      IF(L.GE.1) GOTO 2050
X                      M8=M8+N8(I,JJ)
X                      N(I,JJ) = N(I,JJ)+N8(I,JJ)
X                      N8(I,JJ)=0
X                   ENDIF
X                ENDIF
X 2200        CONTINUE
X 2300     CONTINUE
X       ENDIF
X
X
XC
X       MM=M1+M2+M3+M4
XC
XC      *             INITIAL GUESS OF F                                     *
XC
X       SNT=NTOT
X       DO 2440 I=1,MX
X          DO 2430 J=1,MY
X             IF(N(I,J).NE.0) F(I,J)=FLOAT(N(I,J))/SNT
X 2430     CONTINUE
X 2440  CONTINUE
X
X       RETURN
X       END
X
X
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE BJ   *********************************
XC      **********************************************************************
XC
X       SUBROUTINE BJ(IND,X,Y,NTOT,TOL,MAX,NVAR,ND,NC,ICENS,OUTPUT,
X     +               ALPHA,SIGMAA,
X     +               IWRK1,IWRK2,IWRK4,IWRK5,IWRK6,IWRK7,
X     +               IWRK8,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
X     +               SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
X     +               SWRK8,SWRK9,DWRK1,EWRK1,MVAR)
XC
XC
XC      *    LINEAR REGRESSION WITH CENSORED DATA : BUCKLEY-JAMES METHOD     *
XC      *                                                                    *
XC      *        USING A NONPARAMETRIC METHOD, THIS PROGRAM CALCULATES       *
XC      *    LINEAR REGRESSION COEFFICIENTS FOR DATA WHICH CONTAINS SOME     *
XC      *    CENSORED OBSERVATIONS.                                          *
XC      *                                                                    *
XC      *    PARAMETERS :                                                    *
XC      *     INPUT                                                          *
XC      *     NTOT          : NUMBER OF OBSERVATIONS                         *
XC      *     NVAR          : NUMBER OF INDEPENDENT VARIABLE                 *
XC      *     ALPHA         : INITIAL REGRESSION COEFFICIENT ESTIMATES       *
XC      *                       (PLEASE ALWAYS USE 0.0 IN THIS PROGRAM).     *
XC      *     TOL           : TOLERANCE FOR CONVERGENCE  (E.G. 1.0E-05)      *
XC      *     MAX           : MAXIMUM ITERATION (E.G. 20)                    *
XC      *      X(J,I)       : THE MATRIX CONTAINS THE COEFF.OF J-TH          *
XC      *                     LOCATION PARAMETER AND I-TH OBSERVATION        *
XC      *      Y(I)         : DEPENDENT PARAMETER OF I-TH OBSERVATION        *
XC      *     IND(I)        : INDICATOR OF CENSORED DATA ...                 *
XC      *                     IF IND(I)= 1  : LOWER LIMIT                    *
XC      *                              = 0  : UNCENSORED POINT               *
XC      *                              =-1  : UPPER LIMIT                    *
XC      *      OUTPUT                                                        *
XC      *     ALPHA(1)      : INTERCEPT COEFFICIENT                          *
XC      *     ALPHA(J)      : J>1, J-TH SLOPE COEFFICIENTS                   *
XC      *     ALPHA(MPLONE) : STANDARD DEVIATION                             *
XC      *     SIGMAA(J)     : STANDARD DEVIATION OF J-TH COEFFICIENT         *
XC      *                                                                    *
XC      *    SUBROUTINES                                                     *
XC      *                     BUCKLY                                         *
XC
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION IND(NTOT),X(MVAR,NTOT),Y(NTOT),ALPHA(MVAR),SIGMAA(MVAR)
X       DIMENSION IWRK1(NTOT),IWRK2(NTOT),IWRK4(NTOT)
X       DIMENSION IWRK5(NTOT),IWRK6(NTOT),IWRK7(NTOT),IWRK8(NTOT)
X       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT)
X       DIMENSION WRK5(NTOT),WRK6(NTOT),WRK7(NTOT),WRK8(NTOT)
X       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
X       DIMENSION SWRK5(MVAR),SWRK6(MVAR),SWRK7(MVAR),SWRK8(MVAR)
X       DIMENSION SWRK9(MVAR),DWRK1(MVAR,NTOT),EWRK1(MVAR,MVAR)
X       CHARACTER*9 OUTPUT
XC
X       NVAR1=NVAR+1
XC
XC      *      IF CENSORING IS DUE TO UPPER LIMITS, CHANGE THE SIGNS OF DATA *
XC      *      X(I) AND Y(I) BECAUSE B-J METHOD ASSUMES LOWER LIMITS.        *
XC
X       IF(ICENS .LT. 0) THEN
X          DO 1222 I=1,NTOT
X             DO 1111 J=1,NVAR
X                X(J,I)=-X(J,I)
X 1111        CONTINUE
X             Y(I)=-Y(I)
X 1222     CONTINUE
X       ENDIF
XC
XC
XC      *      BUCKLY   : THE SUBROUTINE WHICH PERFORMS THE BUCKLEY AND      *
XC      *                 JAMES METHOD.                                      *
XC
XC
X       CALL BUCKLY(X,Y,ALPHA,IND,TOL,SIGMAA,NTOT,NVAR,ND,NC,MAX,ITE,
X     +                   IWRK1,DWRK1,IWRK2,IWRK4,WRK1,WRK2,IWRK5,
X     +                   WRK3,WRK4,WRK5,WRK6,WRK7,SWRK1,SWRK2,SWRK3,
X     +                   IWRK6,IWRK7,IWRK8,SWRK4,SWRK5,SWRK6,SWRK7,
X     +                   SWRK8,SWRK9,EWRK1,WRK8,MVAR)
XC
XC
XC       *    CORRECT THE SIGNS OF THE DATA TO THE ORIGINAL ONES, IF THE     *
XC       *    CENSORING IS UPPER LIMIT.                                      *
XC
X       IF(ICENS.LT.0) THEN
X          DO 1223 I=1,NTOT
X
X             DO 1113 J=1,NVAR
X                X(J,I)=-X(J,I)
X 1113        CONTINUE
X             Y(I)=-Y(I)
X 1223     CONTINUE
XC
X          ALPHA(1)=-ALPHA(1)
XC
X       ENDIF
X
X  320  IF(OUTPUT.EQ.'         ') THEN
X          PRINT 1050
X          PRINT 1020
X          PRINT 1050
XC
X          PRINT 1200,ALPHA(1)
X
X          DO 452 J=2,NVAR1
X             JI=J-1
X             PRINT 1250,JI,ALPHA(J),SIGMAA(J)
X  452     CONTINUE
X
X          PRINT 1300,ALPHA(NVAR+2)
X          PRINT 1350,ITE
X          PRINT 1050
X       ELSE
X          WRITE(60,1050)
X          WRITE(60,1020)
X          WRITE(60,1050)
XC
X          WRITE(60,1200) ALPHA(1)
X
X          DO 450 J=2,NVAR1
X             JI=J-1
X             WRITE(60,1250) JI,ALPHA(J),SIGMAA(J)
X  450     CONTINUE
X
X          WRITE(60,1300) ALPHA(NVAR+2)
X          WRITE(60,1350) ITE
X          WRITE(60,1050)
X       ENDIF
XC
XC
X 1020  FORMAT(T5,'LINEAR REGRESSION BY BUCKLEY-JAMES METHOD' )
X 1050  FORMAT(T5,'      ')
X 1100  FORMAT(T8,'DATA TITLE :',T25,60A1)
X 1200  FORMAT(T8,'INTERCEPT COEFF    :',F8.4)
X 1250  FORMAT(T8,'SLOPE COEFF ',I1,'      :',F8.4,T38,'+/-',T41,
X     +         F8.4)
X 1300  FORMAT(T8,'STANDARD DEVIATION :',F8.4)
X 1350  FORMAT(T8,'ITERATIONS         :',I3)
X 4000  RETURN
X       END
X
XC
XC      **********************************************************************
XC      *********************** SUBROUTINE BUCKLY  ***************************
XC      **********************************************************************
XC
X       SUBROUTINE BUCKLY(X,Y,ALPHA,IND,TOL,SIGMAA,NTOT,
X     +                   NVAR,NU,NC,MAX,ITE,IND2,XX,IPT,IR,ND,TY,
X     +                   T,NO,Z,W,WX,ZY,V,TEST,TEST2,BU,
X     +                   IWRK1,IWRK2,SWRK1,SWRK2,SWRK3,SWRK4,
X     +                   SWRK5,SWRK6,EWRK1,WRK1,MVAR)
XC
XC
XC      *     THIS IS A SUBPROGRAM WHICH PERFORMS THE BUCKLEY-JAMES          *
XC      *     METHOD. THIS SUBROUTINE WAS ADAPTED FROM CODE BY J. HALPERN    *
XC      *     (STANFORD UNIVERSITY SCHOOL OF MEDICINE, DEPARTMENT            *
XC      *        OF FAMILY, COMMUNITY AND PREVENTIVE MEDICINE.)              *
XC      *                                                                    *
XC      *  INPUT                                                             *
XC      *         X(J,I)   : INDEPENDENT VARIABLES                           *
XC      *         Y(I)     : DEPENDENT VARIABLE                              *
XC      *         IND(I)   : INDICATOR OF CENSORING                          *
XC      *         TOL      : TOLERANCE LEVEL                                 *
XC      *         NTOT     : NUMBER OF DATA POINTS                           *
XC      *         NVAR     : NUMBER OF INDEPENDENT VARIABLES                 *
XC      *         NU       : NUMBER OF DETECTED POINTS                       *
XC      *         NC       : NUMBER OF CENSORED POINTS                       *
XC      *         MAX      : MAXIMUM ITERATION                               *
XC      *                                                                    *
XC      *   WORK                                                             *
XC      *          V(J)    : AVERAGE OF J-TH DETECTED INDEPENDENT VARIABLE   *
XC      *          BU(J)   : VARIANCE OF J-TH DETECTED INDEPENDENT VARIABLE  *
XC      *         TEST(J)  : STORE OF THE PREVIOUS STEP ESTIMATIONS OF       *
XC      *                    ALPHA(J)                                        *
XC      *          IR(I)   : SORTING ORDER                                   *
XC      *          Z(I)    : RESIDUALS                                       *
XC      *          W(I)    : KM ESTIMATOR                                    *
XC      *          WX(I)   : WEIGHT                                          *
XC      *                                                                    *
XC      *  OUTPUT                                                            *
XC      *       ALPHA(J)   : REGRESSION COEFFICIENTS                         *
XC      *       SIGMA(J)   : ERROR                                           *
XC      *       ITE        : ITERATION NUMBER                                *
XC      *                                                                    *
XC      *    SUBROUTINES                                                     *
XC      *                    SORT1, REGRES                                   *
XC
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
X
X       DIMENSION IND(NTOT),IND2(NTOT),XX(MVAR,NTOT),IPT(NTOT)
X       DIMENSION X(MVAR,NTOT),Y(NTOT),IR(NTOT),ND(NTOT),TY(NTOT)
X       DIMENSION T(NTOT),NO(NTOT),Z(NTOT),W(NTOT),WX(NTOT),ZY(NTOT)
X       DIMENSION V(MVAR),ALPHA(MVAR),TEST(MVAR),BU(MVAR),SIGMAA(MVAR)
X       DIMENSION TEST2(MVAR),IWRK1(NTOT),IWRK2(NTOT)
X       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
X       DIMENSION SWRK5(MVAR),SWRK6(MVAR),EWRK1(MVAR,MVAR),WRK1(NTOT)
XC
XC      *               INITIALIZATION                                       *
XC
X       ITE=0
X       NB=NVAR+1
X       DO 343 J=1,NVAR
X          V(J) =0.0
X          BU(J)=0.0
X  343  CONTINUE
X
X       DO 392 IN=1,NB
X          TEST(IN)=0.0
X       TEST2(IN) =0.0
X  392  CONTINUE
XC
XC      *       CALCULATE SOME VALUES FOR THE STANDARD DEVIATION             *
XC
X       DO 5 I=1,NTOT
X          IF(IND(I).EQ.0) THEN
X             DO 63 J=1,NVAR
X                V(J)=V(J)+X(J,I)
X   63        CONTINUE
X       ENDIF
X    5  CONTINUE
X
X       DO 68 J=1,NVAR
X          V(J)=V(J)/REAL(NU)
X   68  CONTINUE
X
X       DO 51 I=1,NTOT
X          IF(IND(I).EQ.0) THEN
X             DO 53 J=1,NVAR
X                BU(J)=BU(J)+(X(J,I)-V(J))**2
X   53        CONTINUE
X       ENDIF
X   51  CONTINUE
XC
XC      *      REGRES : SUBPROGRAM FOR LINEAR REGRESSION WITHOUT             *
XC      *               CONSIDERING CENSORING STATUS                         *
XC
X       CALL REGRES(X,Y,NTOT,NVAR,ALPHA,RMUL,SWRK1,SWRK2,WRK1,
X     +             IWRK1,IWRK2,SWRK3,SWRK4,SWRK5,EWRK1,SWRK6,MVAR)
XC
XC      *              GET RESIDUALS Z(I)                                    *
XC
XC      *       START ITERATION : 2000 LOOP.                                 *
XC
XC
X 2000  DO 31 I=1,NTOT
X          T(I)=-400.0
X          IND2(I)=IND(I)
XC       
X          ZS=0.0
X          DO 61 J=1,NVAR
X             JJ=J+1
X             ZS=ZS+ALPHA(JJ)*X(J,I)
X             XX(J,I)=X(J,I)
X   61     CONTINUE
XC
X          Z(I)=Y(I)-ZS
X   31  CONTINUE
XC
XC      *            SORTING .... INCREASING ORDER                           *
XC
X       CALL SORT1(IND2,XX,Z,NTOT,NVAR,IR,SWRK1,MVAR)
XC
X       DO 311 I=1,NTOT
X          TY(I)=Y(IR(I))
X          ZY(I)=Z(I)
X  311  CONTINUE
XC
XC      *       THE LARGEST RESIDUAL MUST BE UNCENSORED.                     *
XC
X       IND2(NTOT)=0
XC
XC      *      ESTIMATE  VALUES FOR CENSORED DATA                            *
XC      *                                                                    *
XC      *    TY(I)=YY(I)*DEL+((ALPHA*X+SUM(WXX(K)*Z(K))/(1-W(I)))*(1-DEL)    *
XC      *     WHERE                                                          *
XC      *          TY   : ESTIMATED DEPENDENT VALUE                          *
XC      *          DEL  : IF THE DATA IS UNCENSORED :DEL=1.0                 *
XC      *                 IF THE DATA IS CENSORED  :DEL=0.0                  *
XC      *          SUM  : SUM OVER UNCENSORED DATA Z(K)<Z(I)                 *
XC      *          WX   : WEIGHT ... W(I-1)-W(I)                             *
XC      *          W    : KAPLAN-MEIER PRODUCT LIMIT ESTIMATOR               *
XC
XC
X       K=0
X       DO 21 I=1,NTOT
X          IF((IND2(I).NE.0).OR.(K.NE.0)) THEN
X             IF(IND2(I).NE.0) THEN
X                IPT(I)=K
X                GOTO 21
X             ENDIF
X             IF((IND2(I).EQ.0).AND.(ZY(I).EQ.T(K))) THEN
X                ND(K)=ND(K)+1
X                IPT(I)=K
X                GOTO 21
X             ENDIF
X          ENDIF
X          ND(K+1)=1
X          T(K+1)=ZY(I)
X          IPT(I)=K+1
X          K=K+1
X   21  CONTINUE
XC
X       NI=K
X       DO 28 I=1,NI
X          NO(I)=0
X          IDZ=0
X   28  CONTINUE
XC
X       DO 29 I=1,NTOT
XC
XC      *      IF THE FIRST POINT IS A CENSORED VALUE, DROP THE POINT.       *
XC      *      PT RUNS FROM 1 TO NU.                                         *
XC
X          IF(IPT(I).EQ.0) IDZ=IDZ+1
X          IF(IPT(I).GT.0) NO(IPT(I))=NO(IPT(I))+1
X   29  CONTINUE
X
X       DENOM=REAL(NTOT-IDZ)
X       W(1)=1.0-ND(1)/DENOM
X       WJ=NTOT-NO(1)-IDZ
XC
X       DO 30 I=2,NI
X          W(I)=W(I-1)*(1.0-ND(I)/WJ)
X          WJ=WJ-NO(I)
X   30  CONTINUE
X
X       WX(1)=1.0-W(1)
XC
X       DO 41 I=2,NI
X          WX(I)=W(I-1)-W(I)
X   41  CONTINUE
XC
X       DO 83 I=1,NI
X          WX(I)=WX(I)/REAL(ND(I))
X   83  CONTINUE
XC
X       Z(NTOT)=0.0
X       DO 36 JJ=2,NTOT
X          I=NTOT-JJ+1
X          Z(I)=Z(I+1)
X          IF((ZY(I+1).GT.ZY(I)).AND.(IND2(I+1).EQ.0)) 
X     +                           Z(I)=Z(I+1)+WX(IPT(I+1))*ZY(I+1)
X   36  CONTINUE
XC
X       III=NTOT-1
X       DO 32 I=1,III
X          IF(IPT(I).GT.0) Z(I)=Z(I)/W(IPT(I))
X          IF(IND2(I).NE.0) THEN
X             ZS=0.0
X             DO 62 J=1,NVAR
X                JJ=J+1
X                ZS=ZS+ALPHA(JJ)*XX(J,I)
X   62        CONTINUE
X             TY(I)=Z(I)+ZS
X          ENDIF
X   32  CONTINUE
XC
X       CALL REGRES(XX,TY,NTOT,NVAR,ALPHA,RMUL,SWRK1,SWRK2,WRK1,
X     +             IWRK1,IWRK2,SWRK3,SWRK4,SWRK5,EWRK1,SWRK6,MVAR)
X       ITE=ITE+1
XC
XC       *        TEST FOR CONVERGENCE                                       *
XC       *     TEST(I) CONTAINS A PREVIOUS VALUE OF ALPHA(I).                *
XC       *     IF NUMBER OF ITERATION EXCEEDS MAXITS, THE PROGRAM STOPS      *
XC       *     TO CONTINUE ITERATION, EVEN IF TOLERANCE IS LARGER THAN       *
XC       *     THE ASSIGNED VALUE. SINCE THE FINAL VALUES ARE OFTEN          *
XC       *     TRAPPED IN OSCILLATION, THE PROGRAM TAKES AVERAGE OF THE      *
XC       *     LAST TWO VALUES FOR THE FINAL OUTPUT.                         *
XC
X       IF(ITE.LE.MAX) THEN
X
X          SUM=0.0
X          DO 154 K=1,NB
X             SUM=SUM+(ALPHA(K)-TEST2(K))**2
X             TEST2(K)=TEST(K)
X             TEST(K)=ALPHA(K)
X  154     CONTINUE
XC
XC       *  IF THE DIFFERENCE IS LARGER THAN THE TOLERANCE, REPEAT ITERATION *
XC
X  234     IF(DSQRT(SUM).GT.TOL) GOTO 2000
X       ENDIF
XC
XC       *  THE CONVERSION IS OBTAINED OR THE ITERATION REACHED THE MAX.     *
XC
X       DO 270 K=1,NB
X          ALPHA(K)=(ALPHA(K)+TEST2(K))/2.0
X  270  CONTINUE
XC
XC      *             CALCULATION OF VARIANCE ETC.                           *
XC
X       STD=0.0
X       EM =0.0
X       K=NTOT
X       DO 35 I=1,NTOT
X          IF(IND(I).EQ.0) THEN
X
X             ZS=0.0
X             DO 84 J=1,NVAR
X                JJ=J+1
X                ZS=ZS+ALPHA(JJ)*X(J,I)
X   84        CONTINUE
X
X             Z(I)=Y(I)-ZS
X             EM=EM+Z(I)
X          ENDIF
X   35  CONTINUE
X
X       EM=EM/REAL(NU)
X       DO 37 I = 1,NTOT
X          IF(IND(I) .EQ. 0) THEN
X             EM2 = Z(I) - EM
X             STD = STD+EM2*EM2
X          ENDIF
X   37  CONTINUE
X
X       STD=STD/FLOAT(NU-NB)
X
X       DO 76 I=1,NVAR
X          SIGMAA(I+1)=DSQRT(STD/BU(I))
X   76  CONTINUE
X
X       ALPHA(NB+1)=DSQRT(STD)
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE CENS  ******************************
XC      **********************************************************************
XC
X        SUBROUTINE CENS(X,Y,IP,NTOT)
XC
XC      *         THIS SUBROUTINE ADDS OR SUBTRACTS 0.00001 TIMES THE        *
XC      *         VALUE FROM CENSORED POINTS SO THAT IF THERE ARE TIES WITH  *
XC      *         DETECTED POINTS, THE CENSORED POINTS CAN BE                *
XC      *         DISTINGUISHED FROM THE DETECTIONS.                         *
XC      *         IF THE SMALLEST DIGIT IS LESS THAN OR EQUAL TO             *
XC      *         0.0001, THEN 'CONST' SHOULD BE CHANGED TO A                *
XC      *         SMALLER VALUE.                                             *
XC      *                                                                    *
XC      *         INPUT AND OUTPUT:                                          *
XC      *                   X(I)    : INDEPENDENT VARIABLE                   *
XC      *                   Y(I)    : DEPENDENT VARIABLE                     *
XC      *                  IP(I)    : CENSORED STATUS INDICATOR              *
XC
X        IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X        DIMENSION X(NTOT),Y(NTOT),IP(NTOT)
XC
XC
X        CONST=1.000001
X        DO 500 I=1,NTOT
XC
XC      *                 UPPER LIMITS CASES                                 *
XC
X        IF(IP(I).EQ.-1) THEN 
X              Y(I)=Y(I)/CONST
X
X          ELSEIF(IP(I).EQ.-2) THEN
X              X(I)=X(I)/CONST
X
X          ELSEIF(IP(I).EQ.-3) THEN
X              X(I)=X(I)/CONST
X              Y(I)=Y(I)/CONST
X
X          ELSEIF(IP(I).EQ.-4) THEN
X              X(I)=X(I)/CONST
X              Y(I)=Y(I)*CONST
X
XC
XC      *                  LOWER LIMIT CASES                                 *
XC
X          ELSEIF(IP(I).EQ.1)  THEN
X              Y(I)=Y(I)*CONST
X
X          ELSEIF(IP(I).EQ.2)  THEN
X              X(I)=X(I)*CONST
X
X          ELSEIF(IP(I).EQ.3)  THEN
X              X(I)=X(I)*CONST 
X              Y(I)=Y(I)*CONST
X
X          ELSEIF(IP(I).EQ.4)  THEN
X              X(I)=X(I)*CONST
X              Y(I)=Y(I)/CONST
X
X        ENDIF
X  500   CONTINUE
X        RETURN
X        END
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE BIVAR ******************************
XC      **********************************************************************
XC
X       SUBROUTINE BIVAR(IBACK)
X          
XC
XC      *   CORRELATION AND REGRESSION                                       *
XC      *   PARAMETERS                                                       *
XC      *      MVAR     :  THE MAXIMUM NUMBER OF VARIABLES (COLUMNS)         *
XC      *                   ALLOWED IN A DATA SET.                           *
XC      *      NDAT     :  THE MAXIMUM NUMBER OF DATA POINTS (ROWS)          *
XC      *                   ALLOWED IN A DATA SET.                           *
XC      *      IBIN     :  THE DIMENSION SIZE FOR BINS - USED IN SCHMITT'S   *
XC      *                   PROCEDURE COMPUTATIONS.                          *
XC      *      LENG     :  (MVAR+1)+NTOT - USED IN EM ALGORITHM COMPUTATIONS *
XC      *      LEGWRK   :  (MVAR+1)*NTOT - USED IN EM ALGORITHM COMPUTATIONS *
XC      *   INPUT                                                            *
XC      *      FILE     :  NAME OF DATA FILE (9 LETTERS)                     *
XC      *      TITLE    :  TITLE OF THE PROBLEM (80 LETTERS)                 *
XC      *      NVAR     :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
XC      *      NTOT     :  THE ACTUAL NUMBER OF DATA POINTS IN THE DATA SET  *
XC      *      ICOL     :  INDICATOR OF VARIABLE (<=NVAR)                    *
XC      *                  IF A MULTIVARIATE PROBLEM IS NEEDED, SET ICOL=0   *
XC      *      COLM     :  NAME OF THE INDEPENDENT VARIABLE                  *
XC      *      YNAME    :  NAME OF THE DEPENDENT VARIABLE                    *
XC      *      COMMAND  :  NAME OF THE "COMMAND" FILE                        *
XC      *      OUTPUT   :  NAME OF THE OUTPUT FILE                           *
XC      *      IND(1,I) :  INDICATOR OF CENSORING                            *
XC      *                    IF =0,  DETECTED                                *
XC      *                       =1,  Y LOWER LIMIT                           *
XC      *                       =2,  X LOWER LIMIT                           *
XC      *                       =3,  DOUBLE LOWER LIMIT                      *
XC      *                       =4,  X UPPER LIMIT AND Y LOWER LIMIT         *
XC      *                  FOR THE UPPER LIMITS, CHANGE THE SIGN             *
XC      *                  2, 3, AND 4 CAN BE USED ONLY IN GEN. KENDALL'S    *
XC      *                  TAU, GEN. SPEARMAN'S RHO, AND SCHMITT'S METHOD    *
XC      *      X(J,I)   :  INDEPENDENT VARIABLES                             *
XC      *      Y(I)     :  DEPENDENT VARIABLE                                *
XC      *     IPROG(I)  :  INDICATOR OF METHODS                              *
XC      *     NOTEST    :  NUMBERS OF TEST                                   *
XC      *  INPUT FOR EM ALGORITHM                                            *
XC      *      TOL      :  TOLERANCE (DEFAULT 1.0E-5)                        *
XC      *      MAX      :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *      IBET     :  IF 0, NO DEPENDENT VARIABLE IS CONFINED BETWEEN   *
XC      *                        TWO VALUES                                  *
XC      *                     1, THERE ARE SOME DEPENDENT VARIABLE WHICH     *
XC      *                        ARE CONFINED BETWEEN TWO VALUES             *
XC      *    ALPHA(K)   :  INITIAL ESTIMATE OF REGRESSION COEFFICIENTS       *
XC      *  INPUTS FOR BUCKLEY-JAMES METHOD                                   *
XC      *      TOL1     :  TOLERANCE (DEFAULT 1.0E-5)                        *
XC      *      MAX1     :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *  INPUTS FOR SCHMITT'S BINNING METHOD                               *
XC      *      MX       :  BIN NUMBER OF X AXES                              *
XC      *      MY       :  BIN NUMBER OF Y AXES                              *
XC      *      TOL3     :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                  *
XC      *      MAX3     :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *      XBIN     :  BIN SIZE FOR X AXES                               *
XC      *      YBIN     :  BIN SIZE FOR Y AXES                               *
XC      *      XORG     :  ORIGIN OF X AXES                                  *
XC      *      YORG     :  ORIGIN OF Y AXES                                  *
XC      *      ISKIP    :  IF 0, THE PROGRAM WILL PROVIDE XBIN, YBIN, XORG,  *
XC      *                        AND YORG.                                   *
XC      *                    >0, THESE VALUES MUST BE PROVIDED BY THE USER   *
XC      *      IPIRNT   :  IF 0, NO TWO DIMENSIONAL K-M ESTIMATOR WILL BE    *
XC      *                        PRINTED                                     *
XC      *                    >0, TWO DIMENSIONAL K-M ESTIMATOR WILL BE       *
XC      *                        PRINTED                                     *
XC      *                                                                    *
XC      *    WORKING VARIABLES AND ARRAYS:                                   *
XC      *      NTOT     :  NUMBER OF DATA POINTS                             *
XC      *      ND       :  NUMBER OF DETECTED POINTS                         *
XC      *      NC1      :  NUMBER OF Y LOWER LIMITS                          *
XC      *      NC2      :  NUMBER OF X LOWER LIMITS                          * 
XC      *      NC3      :  NUMBER OF DOUBLE LOWER LIMITS                     * 
XC      *      NC4      :  NUMBER OF Y LOWER AND X UPPER LIMITS              *
XC      *      NC5      :  NUMBER OF Y UPPER LIMITS                          *
XC      *      NC6      :  NUMBER OF X UPPER LIMITS                          *
XC      *      NC7      :  NUMBER OF DOUBLE UPPER LIMITS                     *
XC      *      NC8      :  NUMBER OF Y UPPER AND X LOWER LIMITS              *
XC      *      ICENS    :  IF 0, CENSORING IS MIXED                          *
XC      *                     1, CENSORING IS Y LOWER LIMITS ONLY            *
XC      *                    -1, CENSORING IS Y UPPER LIMITS ONLY            *
XC      *      NYC      :  NC1+NC2                                           *
XC      *      NXC      :  NC3+NC4                                           *
XC      *      NBC      :  NC5+NC6+NC7+NC8                                   *
XC      *      IDO      :  NXC+NBC                                           *
XC      *      IMUL     :  INDICATOR OF MULTIVARIATE PROBLEM                 *
XC      *      XX(J,I)  :  =X(ICOL,I), EXCEPT FOR MULTI INDEPENDENT VARIABLE *
XC      *                  CASE (J=1,NVAR).                                  *
XC      *      IND2(I)  :  =IND(1,I)                                         *
XC      *                                                                    *
XC      *  OUTPUT                                                            *
XC      *     COXREG                                                         *
XC      *      CHI      : GLOBAL CHI-SQUARE                                  *
XC      *      PROB     : PROBABILITY FOR NULL                               *
XC      *     BHK (GNERALIZED KENDALL'S TAU)                                 *
XC      *       Z       : DEVIATION                                          *
XC      *      PROB     : PROBABILITY FOR NULL                               *
XC      *     EM ALGORITHM                                                   *
XC      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS  (K=1,NVAR+1)       *
XC      *     ALPHA(K+2): STANDARD DEVIATION                                 *
XC      *     SIGMAA(K) : ERROR                                              *
XC      *     ITE       : NUMBER OF ITERATION                                *
XC      *     BUCKLEY-JAMES                                                  *
XC      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS (K=1,NVAR+1)        *
XC      *     ALPHA(K+2): STANDARD DEVIATION                                 *
XC      *     SIGMAA(K) : ERROR                                              *
XC      *     SCHMITT                                                        *
XC      *     ALPHA     : INTERCEPT COEFFICIENT                              *
XC      *     BETA      : SLOPE COEFFICIENT                                  *
XC      *   *****  ALL OUTPUTS ARE INSIDE OF EACH SUBROUTINE                 *
XC      *                                                                    *
XC      *   SUBROUTINES                                                      *
XC      *     DATA1, DATREG, DATA2, MULVAR                                   *
XC      *                                                                    *
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
XC      *   THIS PARAMETER STATEMENT AND THE ONE IN UNIVAR.F ARE THE ONLY    *
XC      *   STATEMENTS THAT NEED TO BE ADJUSTED IF THE USER WISHES TO        *
XC      *   ANALYZE DATA SETS OF MORE THAN 500 OBSERVATIONS OR MORE THAN     *
XC      *   VARIABLES.                                                       *
X
XC  **************************************************************************
X       PARAMETER(MVAR=4, NDAT=500, IBIN=50)
XC  **************************************************************************
X
X       CHARACTER*1 CHECK,CHAR(4,10)
X       CHARACTER*7 BB(10),YY
X       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM
X       CHARACTER*9 YNAME,DUMMY1
X       CHARACTER*80 TITLE 
X       
X       DIMENSION IND(MVAR,NDAT),X(MVAR,NDAT),Y(NDAT)
X       DIMENSION IPROG(6),IIND(NDAT)
X
X       DIMENSION DWRK1(MVAR,NDAT),DWRK2(MVAR,NDAT)
X       DIMENSION DWRK3(MVAR,NDAT),DWRK4(MVAR,NDAT)
X       DIMENSION DWRK5(MVAR,NDAT),DWRK6(MVAR,NDAT)
X       DIMENSION DWRK8(MVAR,NDAT)
X       
X       DIMENSION EWRK1(4,4),RWRK1(NDAT,MVAR)
X
X       DIMENSION AWRK(5,IBIN)
X       DIMENSION WWRK1((MVAR+1)+NDAT)
X       DIMENSION WWRK2((MVAR+1)+NDAT)
X       DIMENSION VWRK1((MVAR+1)*NDAT)
X       DIMENSION VWRK2((MVAR+1)*NDAT)
X
X       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
X       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
X       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
X       DIMENSION WRK12(NDAT)
X
X       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR)
X       DIMENSION SWRK4(MVAR),SWRK5(MVAR),SWRK6(MVAR)
X       DIMENSION SWRK7(MVAR),SWRK8(MVAR),SWRK9(MVAR)
X       DIMENSION SWRK10(MVAR),SWRK11(MVAR),SWRK17(MVAR)
X
X       DIMENSION LWRK1(MVAR,NDAT), LWRK2(MVAR,NDAT)
X       DIMENSION LWRK3(MVAR,NDAT)
X       DIMENSION IWRK1(NDAT),IWRK2(NDAT),IWRK3(NDAT)
X       DIMENSION IWRK4(NDAT),IWRK5(NDAT),IWRK6(NDAT)
X       DIMENSION IWRK7(NDAT),IWRK8(NDAT)
X
X       DIMENSION IWRK9(NDAT),CWRK1(IBIN),CWRK2(IBIN)
X
X       DIMENSION IBWRK1(IBIN,IBIN),IBWRK2(IBIN,IBIN)
X       DIMENSION IBWRK3(IBIN,IBIN),IBWRK4(IBIN,IBIN)
X       DIMENSION IBWRK5(IBIN,IBIN),IBWRK6(IBIN,IBIN)
X       DIMENSION IBWRK7(IBIN,IBIN),IBWRK8(IBIN,IBIN)
X       DIMENSION IBWRK9(IBIN,IBIN)
X       DIMENSION BWRK1(IBIN,IBIN),BWRK2(IBIN,IBIN)
X
X       LENG = (MVAR+1)+NDAT
X       LEGWRK = (MVAR+1)*NDAT
X
X       DO 5000 K=1,10
X       BB(K)='     X '
X 5000  CONTINUE
X       YY='     Y '
XC
X 6000  PRINT *
X       PRINT *
X       PRINT *,'      CORRELATION AND REGRESSION CALCULATIONS'
X       PRINT *
X       PRINT *,' CORRELATION OPTIONS        LINEAR REGRESSION OPTIONS'
X       PRINT *,' 1. COX HAZARD MODEL        4. EM ALGORITHM WITH '
X       PRINT *,'                               NORMAL DISTRIBUTION'
X       PRINT *,' 2. GEN. KENDALL`S TAU      5. BUCKLEY-JAMES METHOD'
X       PRINT *,' 3. GEN. SPEARMAN`S RHO     6. SCHMITT`S BINNING METHOD'
X       PRINT *
X       PRINT *,' DATA SETS WITH CENSORING IN ONLY ONE DIRECTION OF THE'
X       PRINT *,' DEPENDENT VARIABLE CAN USE ALL METHODS.'
X       PRINT *
X       PRINT *,' DATA SETS WITH SEVERAL INDEPENDENT AND ONE DEPENDENT'
X       PRINT *,' VARIABLE CAN USE ONLY THE COX PROPORTIONAL HAZARD'
X       PRINT *,' MODEL,EM ALGORITHM, OR BUCKLEY-JAMES METHOD.  ONLY'
X       PRINT *,' ONE TYPE OF CENSORING IN THE DEPENDENT VARIABLE IS'
X       PRINT *,' ALLOWED.'
X       PRINT *
X       PRINT *,' IF YOUR DATA SET HAS CENSORED POINTS IN THE '
X       PRINT *,' INDEPENDENT VARIABLE AND/OR DUAL CENSORED POINTS,'
X       PRINT *,' YOU CAN USE ONLY THE GEN. KENDALL`S TAU OR GEN.'
X       PRINT *,' SPEARMAN`S RHO CORRELATION COEFFICIENT, OR'
X       PRINT *,' SCHMITT`S BINNED LINEAR REGRESSION.'
X       PRINT *
X 6010  PRINT *
XC
XC      *  CHECK WHETHER THE USER WANTS TO USE COMMAND FILE INPUTS. IF SO,   *
XC      *  GO TO 6660                                                        *
XC
X   50  FORMAT(A1)
X 1380  FORMAT(A9)
XC
X       OUTPUT='         '
X       ICOMM=0
X       ICOL=1
X       PRINT *,'DO YOU WANT TO READ ALL INFORMATION'
X       WRITE(6,6020)
X 6020  FORMAT('      FROM A COMMAND FILE (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6660
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6025
X       GOTO 6010
XC
XC      *          READ FROM THE TERMINAL                                    *
XC
X 6025  PRINT *
X       PRINT *,'          OK, LET US READ FROM THE TERMINAL '
XC
XC      *           READ TITLE                                               *
XC
X 6030  PRINT *
X       WRITE(6,6040)
X 6040  FORMAT('WHAT IS THE TITLE OF THE PROBLEM ? ')
X       READ(5,6050) TITLE
X 6050  FORMAT(A80)
XC
XC      *           READ DATA FILE NAME                                      *
XC
X 6051  PRINT *
X       WRITE(6,6052)
X 6052  FORMAT('WHAT IS THE DATA FILE NAME ? ')
X       READ(5,1380) FILE
XC
XC      *           READ NUMBER OF INDEPENDENT VARIABLES                     *
XC
X 6060  PRINT *
X       WRITE(6,6070)
X 6070  FORMAT('HOW MANY INDEPENDENT VARIABLES DO YOU HAVE ? ')
X       CALL DATA1(NVAR)
X       IF((NVAR.GE.1).AND.(NVAR.LE.MVAR-2)) GOTO 6080
X       PRINT *
X       PRINT *,'    YOUR CHOICE IS NOT ACCEPTABLE. PLEASE TYPE IT AGAIN'
X       GOTO 6060
XC
XC      *    CALL SUBROUTINE "DATREG" TO READ DATA                           *
XC
X 6080  CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
X     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
X       DO 6090 I = 1, NTOT
X            IIND(I) = IND(1,I)    
X 6090  CONTINUE
XC
XC      *        CHECK WHICH METHODS THE USER CAN USE                        *
XC
X       IDC=NXC+NBC
X       IF((NVAR.EQ.1).AND.(IDC.EQ.0)) GOTO 6530
X       IF((NVAR.NE.1).AND.(IDC.NE.0)) GOTO 6340
X       IF((NVAR.EQ.1).AND.(IDC.NE.0)) GOTO 6400
X 6170  PRINT *
X       WRITE(6,6180)
X 6180  FORMAT('IS THIS A MULTIVARIATE PROBLEM  (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6220
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6340
X       GOTO 6170
XC
XC      *      DATA SET WITH MORE THAN ONE INDEPENDENT VARIABLES             *
XC
X 6220  PRINT *
X       PRINT *,'           YOU CAN USE THE NEXT METHODS   '
X       PRINT *
X       PRINT *,'         1. COX HAZARD METHOD'
X       PRINT *,'         4. EM ALGORITHM WITH NORMAL DISTRIBUTION'
X       PRINT *,'         5. BUCKLEY-JAMES METHOD'
XC
X       ICOL=0
X       J=1
X 6230  PRINT *
X       WRITE(6,6240)
X 6240  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
X       CALL DATA1(IPROG(J))
X       IF((IPROG(J).EQ.1).OR.(IPROG(J).EQ.4).OR.(IPROG(J).EQ.5))
X     +                                                    GOTO 6245
X       GOTO 6230
X 6245  IF(J.EQ.1) GOTO 6260
X       J1=J-1
X       DO 6250 K=1,J1
X       IF(IPROG(K).NE.IPROG(J)) GOTO 6250
X       PRINT *
X       PRINT *,' YOU ALREADY CHOSE THAT METHOD.'
X       PRINT *,'  PLEASE CHOOSE ANOTHER ONE'
X       GOTO 6230
X 6250  CONTINUE
X 6260  IF(J.GE.3) GOTO 6280
X 6265  PRINT *
X       WRITE(6,6270)
X 6270  FORMAT('DO YOU WANT TO USE ANY OTHER METHOD (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6230
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6280
X       GOTO 6265
XC
X 6280  NOTEST=J
X       GOTO 6604
XC
XC      *    FOR THE CASE  THAT THE DATA SET CONTAINS MIXED CENSORING        *
XC      *    (THAT IS, UPPER AND LOWER LIMITS SIMULTANEOUSLY AND/OR          *
XC      *    CENSORING IN BOTH VARIABLES).                                   *
XC
X 6340  PRINT *
X       WRITE(6,6350)
X 6350  FORMAT('WHICH INDEPENDENT VARIABLE DO YOU WANT TO USE ? ')
X       CALL DATA1(ICOL)
X       IF(ICOL.GT.NVAR) GOTO 6340
X       IF(ICOL.LE.0) GOTO 6340
XC
X 6400  IF(NBC.EQ.0) GOTO 6530
X       J=1
X       PRINT *
X       PRINT *,'          YOU CAN USE THE FOLLOWING METHODS'
X       PRINT *,'          2. GEN. KENDALL`S TAU METHOD'
X       PRINT *,'          3. GEN. SPEARMAN`S RHO METHOD'
X       PRINT *,'          6. SCHMITT`S BINNING METHOD'
X 6410  PRINT *
X       WRITE(6,6420)
X 6420  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
X       CALL DATA1(IPROG(J))
X       IF((IPROG(J).EQ.2).OR.(IPROG(J).EQ.3).OR.(IPROG(J).EQ.6)) 
X     +                                                     GOTO 6425
X       GOTO 6410
X 6425  IF(J.EQ.1) GOTO 6440
X       J1=J-1
X       DO 6430 K=1,J1
X       IF(IPROG(K).NE.IPROG(J)) GOTO 6430
X       PRINT *
X       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
X       PRINT *,'    PLEASE CHOOSE THE OTHER ONE'
X       GOTO 6410
X 6430  CONTINUE
X 6440  IF(J.EQ.3) GOTO 6600
X 6450  PRINT *
X       WRITE(6,6460)
X 6460  FORMAT('DO YOU WANT TO USE THE OTHER METHOD, TOO (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6410
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
X       GOTO 6450
XC
XC      *  FOR THE CASE THAT THE DATA SET CONTAINS ONE INDEPENDENT AND ONE   *
XC      *  DEPENDENT VARIABLES AND ONE KIND OF CENSORING IN THE DEPENDENT    *
XC      *  VARIABLE.                                                         *
XC
X 6530  PRINT *
X       PRINT *,'     YOU CAN USE THE FOLLOWING METHODS'
X       PRINT *
X       PRINT *,'     1. COX HAZARD MODEL    4. EM ALGORITHM WITH'
X       PRINT *,'                               NORMAL DISTRIBUTION'
X       PRINT *,'     2. KENDALL`S TAU       5. BUCKLEY-JAMES REGRESSION'
X       PRINT *,'     3. SPEARMAN`S RHO',
X     +         '      6. SCHMITT`S BINNED REGRESSION'
X       PRINT *
X       J=1
X 6540  PRINT *
X       WRITE(6,6550)
X 6550  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
X       CALL DATA1(IPROG(J))
X       IF(IPROG(J).LT.1) GOTO 6540
X       IF(IPROG(J).GT.6) GOTO 6570
X       IF(J.EQ.1) GOTO 6580
X       J1=J-1
X       DO 6560 K=1,J1
X       IF(IPROG(K).NE.IPROG(J)) GOTO 6560
X       PRINT *
X       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
X       PRINT *,'    PLEASE CHOOSE ANOTHER ONE.'
X       GOTO 6540
X 6560  CONTINUE
X 6570  IF(J.GE.6) GOTO 6600
X 6580  PRINT *
X       WRITE(6,6590)
X 6590  FORMAT('DO YOU WANT TO USE ANOTHER METHOD (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6540
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
X       GOTO 6580
XC
X 6600  NOTEST=J
XC
XC      *        NAME THE VARIABLES                                          *
XC
X 6601  PRINT *
X       PRINT *,'          PLEASE NAME THE VARIABLES : '
X       PRINT *
X       WRITE(6,6602)
X 6602  FORMAT('WHAT IS THE NAME OF THE INDEPENDENT VARIABLE ? ')
X       READ(5,1380) COLM
XC
X       PRINT *
X       WRITE(6,6603)
X 6603  FORMAT('WHAT IS THE NAME OF THE DEPENDENT VARIABLE ? ')
X       READ(5,1380) YNAME
XC
X 6604  PRINT *
X       WRITE(6,6605)
X 6605  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA  (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
X       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 6609
X       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 6609
X       GOTO 6604
XC
XC      *   CHECK WHETHER THE USER WANT TO SAVE THE RESULT IN AN OUTPUT FILE *
XC
X 6609  PRINT *
X       PRINT *,'DO YOU WANT TO SAVE THE RESULT '
X       WRITE(6,6610)
X 6610  FORMAT('     IN AN OUTPUT FILE (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6620
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 7116
X       GOTO 6609
X 6620  PRINT *
X       WRITE(6,6630)
X 6630  FORMAT('WHAT IS THE NAME OF THE OUTPUT FILE ? ')
X       READ(5,1380) OUTPUT
X       GOTO 7116
XC
XC
XC
XC      *    USE "COMMAND" FILE FOR INPUTS                                   *
XC
X 6660  PRINT *
X       WRITE(6,6670)
X 6670  FORMAT('WHAT IS THE NAME OF THE COMMAND FILE ? ')
X       READ(5,1380) COMMAND
XC
X 6700  OPEN(UNIT=50, FILE=COMMAND, STATUS='OLD', FORM='FORMATTED')
X       ICOMM=1
XC
XC      *   READ TITLE OF THE PROBLEM ; NAME OF THE DATA FILE                *
XC
X       READ(50,6710) TITLE
X 6710  FORMAT(A80)
X       READ(50,1380) FILE
XC
XC      * READ NUMBER OF VARIABLES; WHICH VARIABLE WILL BE USED; AND HOW     *
XC      * MANY METHODS THE USER WANTS TO USE.                                *
XC
X       READ(50,6720) ((CHAR(I,J),I=1,4),J=1,3)
X 6720  FORMAT(20A1)
X       CALL DATA2(CHAR,1,3,NVAR,LIND)
X       IF(LIND.EQ.0) GOTO 6750
X 6730  PRINT *
X       PRINT *,'   TOTAL NUMBER OF INDEPENDENT VARIABLES IS NOT CLEAR.'
X       STOP
X 6750  IF(NVAR.LT.1) GOTO 6730
X       IF(NVAR.NE.1) GOTO 6760
X       ICOL=1
X       GOTO 6915
XC 
X 6760  CALL DATA2(CHAR,50,3,ICOL,LIND)
X       IF(LIND.EQ.0) GOTO 6900
X 6860  PRINT *
X       PRINT *,'     THE CHOICE OF THE VARIABLE IS NOT CLEAR'
X       STOP
X 6900  IF(ICOL.LE.0) GOTO 6860
X       IF(ICOL.GT.NVAR) GOTO 6860
XC
XC      *         CHOICE OF THE METHODS                                      *
XC
X 6915  CALL DATA2(CHAR,3,3,NOTEST,LIND)
X       IF(LIND.EQ.0) GOTO 6950
X 6930  PRINT *
X       PRINT *,'    IT IS NOT CLEAR HOW MANY METHODS YOU WANT TO USE '
X       STOP
X 6950  IF(NOTEST.LE.0) GOTO 6930
X       IF(NOTEST.GT.6) GOTO 6930
XC
X       READ(50,6960) ((CHAR(I,J),I=1,4),J=1,NOTEST)
X 6960  FORMAT(30A1)
X       DO 7020 I=1,NOTEST
X       CALL DATA2(CHAR,I,NOTEST,IPROG(I),LIND)
X       IF(LIND.EQ.0) GOTO 7010
X 6970  PRINT *
X       IF(I.EQ.1) PRINT *,'     FIRST PROGRAM NUMBER IS NOT CLEAR'
X       IF(I.EQ.2) PRINT *,'     SECOND PROGRAM NUMBER IS NOT CLEAR'
X       IF(I.GE.3) WRITE(6,6780) I
X 6780  FORMAT(5X,I4,'-TH PROGRAM NUMBER IS NOT CLEAR')
X       STOP
X 7010  IF(IPROG(I).LE.0) GOTO 6970
X       IF(IPROG(I).GT.6) GOTO 6970
X 7020  CONTINUE
XC
XC      *  READ NAMES OF THE INDEPENDENT AND DEPENDENT VARIABLES. IF THE     *
XC      *  PROBLEM HAS MULTI-INDEPENDENT VARIABLES, THESE NAMES WIL BE       *
XC      *  IGNORED.                                                          *
XC
X       READ(50,7022) COLM,YNAME
X 7022  FORMAT(2A9)
XC
X       CLOSE(UNIT=50,STATUS='KEEP')
XC
XC      *      CALL SUBROUTINE "DATREG" TO READ DATA                         *
XC
X       CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
X     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
XC
X       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
XC
XC      * THE NEXT SEVERAL LINES READ IN DUMMY VALUES TO PREVENT READING     *
XC      * THE COMMANDS A SECOND TIME.                                        *
XC
X       READ(50,1380) DUMMY1
X       READ(50,1380) DUMMY1
X       READ(50,7029) IDUMMY
X       READ(50,7029) IDUMMY
X       READ(50,1380) DUMMY1
X 7029  FORMAT(I4)
XC
XC      *  CHECK WHETHER THE ASSIGNED METHODS CAN BE USED FOR THE DATA       *
XC
X       IF(NVAR.GE.2) GOTO 7070
X       IF((NXC.EQ.0).AND.(NBC.EQ.0)) GOTO 7110
XC
XC      *   THE CASE WITH MIXED CENSORING IN DATA                            *
XC
X       I=1
X 7030  IF(IPROG(I).NE.1) GOTO 7040
X       PRINT *
X       PRINT *,'      YOU CANNOT USE COX HAZARD MODEL FOR THIS DATA SET'
X       PRINT *,'      THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7040  IF(IPROG(I).NE.4) GOTO 7050
X       PRINT *
X       PRINT *,'       YOU CANNOT USE EM ALGORITHM FOR THIS DATA SET'
X       PRINT *,'       THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7050  IF(IPROG(I).NE.5) GOTO 7060
X       PRINT *
X       PRINT *,'       YOU CANNOT USE BUCKLEY-JAMES METHOD FOR THIS'
X       PRINT *,'        DATA SET. THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7060  IF(I.GE.NOTEST) GOTO 7110
X       I=I+1
X       GOTO 7030
XC
XC      *     THE CASE WITH MORE THAN ONE INDEPENDENT VARIABLES              *
XC
X 7070  I=1
X 7080  IF(IPROG(I).NE.2) GOTO 7085
X       PRINT *
X       PRINT *,'         YOU CANNOT USE THE KENDALL`S TAU METHOD FOR'
X       PRINT *,'         THIS DATA SET'
X       PRINT *,'         THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7085  IF(IPROG(I).NE.3) GOTO 7090
X       PRINT *
X       PRINT *,'       YOU CANNOT USE SPEARMAN`S RHO FOR THIS DATA SET'
X       PRINT *,'       THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7090  IF(IPROG(I).NE.6) GOTO 7100
X       PRINT *
X       PRINT *,'         YOU CANNOT USE SCHMITT`S BINNED REGRESSION'
X       PRINT *,'         THIS METHOD WILL BE IGNORED'
X       IPROG(I)=-9
X 7100  IF(I.EQ.NOTEST) GOTO 7110
X       I=I+1
X       GOTO 7080
XC
XC      *        READ PRINT OUT INDICATOR FOR THE DATA                       *
XC
X 7110  READ(50,6960) (CHAR(I,1),I=1,4)
X       CALL DATA2(CHAR,1,1,IDATA,LIND)
X       IF(LIND.EQ.0) GOTO 7114
X 7112  PRINT *
X       PRINT *,'     THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
X       STOP
X 7114  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 7115
X       GOTO 7112
XC
XC      *       READ OUTPUT FILE NAME                                        *
XC
X 7115  READ(50,1380) OUTPUT
XC
XC      * CALL SUBROUTINE "MULVAR" TO COMPUTE CORRELATION/REGRESSION PROBLEMS*
XC
X
X 7116  IF(OUTPUT .NE. '         ') OPEN(UNIT=60,FILE=OUTPUT,
X     +                                 STATUS='NEW'
XC  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
XC  VAX/VMS MACHINES.
XC     +                                 ,CARRIAGECONTROL='LIST'
X     +                                 )
X
X
X       CALL  MULVAR(X,Y,IND,NTOT,ICOL,NVAR,NOTEST,IPROG,ICOMM,
X     +              OUTPUT,COLM,FILE,YNAME,TITLE,ND,NYC,ICENS,
X     +              NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,MVAR,
X     +              LENG,LEGWRK,IBIN,DWRK1,IWRK9,SWRK17,DWRK2,
X     +              DWRK3,DWRK4,DWRK5,DWRK6,DWRK8,RWRK1,
X     +              EWRK1,AWRK,WWRK1,WWRK2,
X     +              VWRK1,VWRK2,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
X     +              WRK7,WRK8,WRK9,WRK10,WRK11,WRK12,
X     +              SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
X     +              SWRK8,SWRK9,SWRK10,SWRK11,LWRK1,LWRK2,LWRK3,
X     +              IWRK1,IWRK2,IWRK3,IWRK4,IWRK5,IWRK6,IWRK7,
X     +              IWRK8,CWRK1,CWRK2,IBWRK1,IBWRK2,IBWRK3,
X     +              IBWRK4,IBWRK5,IBWRK6,IBWRK7,IBWRK8,IBWRK9,
X     +              BWRK1,BWRK2)
XC
X       IF(IDATA.EQ.0) GOTO 7219
X       IF(OUTPUT.NE.'         ') WRITE(60,7140)
X       IF(OUTPUT.NE.'         ') WRITE(60,7117) FILE
X       IF(OUTPUT.EQ.'         ') PRINT 7140
X       IF(OUTPUT.EQ.'         ') PRINT 7117, FILE
X 7117  FORMAT(5X,'INPUT DATA FILE : ',A9)
X       IF(ICOL.NE.0) GOTO 7130
X       IF(OUTPUT.NE.'         ') WRITE(60,7118) (BB(K),K,K=1,NVAR),YY
X       IF(OUTPUT.EQ.'         ') PRINT 7118,(BB(K),K,K=1,NVAR),YY
X 7118  FORMAT(4X,'CENSORSHIP',12(A7,I2,1X))
X       DO 7119 I=1,NTOT
X       IF(OUTPUT.NE.'         ') WRITE(60,7120) IIND(I),
X     +                                      (X(J,I),J=1,NVAR),Y(I)
X       IF(OUTPUT.EQ.'         ') PRINT 7120,IIND(I),
X     +                                      (X(J,I),J=1,NVAR),Y(I)
X 7119  CONTINUE
X 7120  FORMAT(7X,I4,3X,10F10.3)
X       GOTO 7219
X 7130  IF(OUTPUT.NE.'         ') WRITE(60,7133)
X       IF(OUTPUT.EQ.'         ') PRINT 7133
X 7133  FORMAT(5X,' CENSORSHIP        X        Y')
X       DO 7134 I=1,NTOT
X       IF(OUTPUT.NE.'         ') WRITE(60,7135) IIND(I),X(ICOL,I),Y(I)
X       IF(OUTPUT.EQ.'         ') PRINT 7135,IIND(I),X(ICOL,I),Y(I)
X 7134  CONTINUE
X 7135  FORMAT(8X,I4,5X,2F10.3)
X 7140  FORMAT('     ')
XC
X 7219  IF(OUTPUT .NE. '         ') CLOSE(UNIT=60)
X       PRINT *
X       PRINT *
X       PRINT *,'    COMPUTATIONS FOR CORRELATION/REGRESSION'
X       PRINT *,'                          PROBLEMS ARE FINISHED'
X 7220  PRINT *
X       WRITE(6,7230)
X 7230  FORMAT('DO YOU WANT TO DO ANY OTHER ANALYSIS (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1 
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
X       GOTO 7220
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE CHOL  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE CHOL(A,N,U,NULLTY,NA,NU,IFAULT)
XC
XC      *      ALGORITHM AS 6 J.R.STATIST.SOC.C.(1968) VOL.17, NO.2          *
XC      *                                                                    *
XC      *      GIVEN A SYMMETRIC MATRIX OF ORDER N AS A LOWER TRIANGLE       *
XC      *      IN A( ),  CALCULATE AN UPPER TRIANGLE, U( ), SUCH THAT        *
XC      *      UPRIME*U=A. U( ) MAY COINCIDE WITH A( ). A( ) MUST BE         *
XC      *      POSITIVE SEMIDEFINITE.                                        *
XC      *      ETA IS SET TO MULTIPLYING FACTOR DETERMINING THE              *
XC      *      EFFECTIVE  ZERO FOR PIVOT.                                    *
XC      *      NULLTY IS RETURNED AS NO. OF EFFECTIVE ZERO PIVOTS.           *
XC      *      IFAULT IS RETURNED AS 1,IF N.LE.0, 2,IF A( ) IS NOT           *
XC      *      POSITIVE SEMI-DEFINITE WITHIN THE TOLERANCE BY ETA.           *
XC      *      OTHERWISE ZERO.                                               *
XC
XC      *        NOTE : VARIABLES NA,NU, HAVE BEEN ADDED TO THE              *
XC      *               ARGUMENT LIST AND USED TO DIMENSION TO ARRAYS        *
XC      *               A AND U, RESPECTIVELY. (BY WOLYNETZ (1979))          *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION A(NA),U(NU)
XC
X       DATA ETA /1.0E-9/
XC
XC      *       THE VALUE OF ETA WILL DEPEND ON THE WORD LENGTH OF           *
XC      *       THE COMPUTER BEING USED.                                     *
XC
X       IFAULT=1
X       IF(N.GT.0) THEN
X          IFAULT=2
X          NULLTY=0
X          J=1
X          K=0
X
X          DO 10 ICOL=1,N
X             L=0
X
X             DO 11 IROW=1,ICOL
X                K=K+1
X                W=A(K)
X                M=J
X
X                DO 12 I=1,IROW
X                   L=L+1
X                   IF(I.EQ.IROW) GOTO 13
X                   W=W-U(L)*U(M)
X                   M=M+1
X   12           CONTINUE
X
X   13           IF(IROW.EQ.ICOL) GOTO 14
X                IF(U(L).EQ.0.0) THEN
X                   U(K) = 0.0
X                ELSE
X                   U(K)=W/U(L)
X                ENDIF
X   11        CONTINUE
X
X   14        IF(DABS(W).GE.DABS(ETA*A(K))) THEN
X                IF(W.LT.0.0) GOTO 100
X                U(K)=DSQRT(W)
X             ELSE
X                U(K)=0.0
X                NULLTY=NULLTY+1
X             ENDIF
X             J=J+ICOL
X   10     CONTINUE
X
X          IFAULT=0
X
X       ENDIF
X  100  RETURN
X       END
X
XC
XC      **********************************************************************
XC      ************************ SUBROUTINE COEFF  ***************************
XC      **********************************************************************
XC
X       SUBROUTINE COEFF(I,X,IP,NTOT,ICOEFF,IA,IB,IC,ID,IE,IG,IH,IJ)
XC
XC      *        SUBROUTINE WHICH FINDS CONCORDANCE INFORMATION  OF          *
XC      *        THE QUANTITY X(I).                                          *
XC      *                                                                    *
XC      *      INPUT  :   X(I)      : THE QUANTITIY TO BE EXAMINED           *
XC      *                IP(I)      : CENSORED STATUS OF X(I)                *
XC      *                NTOT       : NUMBER OF DATA                         *
XC      *      OUTPUT : ICOEFF(I)   : CONCORDANCE INFORMATION:               *
XC      *                             FOR X(I) AND X(J) WITH I<J,            *
XC      *                             IF X(I)<X(J), ICOEFF= 1                *
XC      *                             IF X(I)>X(J), ICOEFF=-1                *
XC      *                             OTHERWISE,    ICOEFF= 0                *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION X(NTOT),IP(NTOT),ICOEFF(NTOT)
XC
X       DO 100 J=1,NTOT
X          ICOEFF(J)=0
X          IF(X(I).LT.X(J)) THEN
X             IF(IP(I).EQ.IA) GOTO 100
X             IF(IP(J).EQ.ID) GOTO 100
X             IF(IP(I).EQ.IB) GOTO 100
X             IF(IP(J).EQ.IE) GOTO 100
X             IF(IP(I).EQ.IC) GOTO 100
X             IF(IP(J).EQ.IG) GOTO 100
X
X             ICOEFF(J)=1
X
X          ELSEIF(X(I).GT.X(J)) THEN
X   50        IF(IP(I).EQ.ID) GOTO 100
X             IF(IP(J).EQ.IA) GOTO 100
X             IF(IP(I).EQ.IE) GOTO 100
X             IF(IP(J).EQ.IB) GOTO 100
X             IF(IP(I).EQ.IG) GOTO 100
X             IF(IP(J).EQ.IC) GOTO 100
X
X             ICOEFF(J)=-1
X          ENDIF
X  100  CONTINUE
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE COXREG  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE COXREG(IND,XX,YY,NTOT,NVAR,OUTPUT,ICENS,
X     +                      RINFO,SCORE,FINFO,IL,IM,IP,Y,X,
X     +                      SWRK1,IWRK1,IWRK2,MVAR)
XC
XC      *      THIS PROGRAM COMPUTES A CORRELATION PROBABILITY ACCORDING     *
XC      *      TO COX'S (1972) PROPORTIONAL HAZARDS MODEL.                   *
XC      *      ONLY ONE TYPE OF CENSORING (I.E. LOWER OR UPPER)              *
XC      *      IS ALLOWED IN Y, BUT UP TO NVAR INDEPENDENT VARIABLES CAN     *
XC      *      BE USED. THE HYPOTHESIS TESTED IS THE ABSENCE OF CORRELATION  *
XC      *      BETWEEN THE DEPENDENT VARIABLE AND INDEPENDENT VARIABLES.     *
XC      *      THEREFORE, THE  REGRESSION COEFFICIENT IN COX MODEL BETA      *
XC      *      IS SET TO ZERO.                                               *
XC
XC      * NOTE NOTE NOTE:   THE PROBABILITY CALCULATED MAY NOT BE ACCURATE   *
XC      *      WHEN THERE ARE A LARGE NUMBER OF TIED OBSERVATIONS (CF.       *
XC      *      R. G. MILLER, SURVIVAL ANALYSIS, 1981, PP. 136-7).            *
XC
XC      *                                                                    *
XC      *      INPUT    IND(I)  :  INDICATOR OF CENSORING                    *
XC      *                           0 : UNCENSORED DATA POINT                *
XC      *                           1 : Y(I) IS LOWER LIMIT                  *
XC      *                          -1 : Y(I) IS UPPER LIMIT                  *
XC      *              XX(J,I)  :  INDEPENDENT VARIABLES (J=1,..NVAR)        *
XC      *              YY(I)    :  DEPENDENT VARIABLE                        *
XC      *               NTOT    :  TOTAL NO. OF OBSERVATIONS                 *
XC      *               NVAR    :  NO. OF INDEPENDENT VARIABLES              *
XC      *                                                                    *
XC      *       WORK      DF    :  DEGREE OF FREEDOM                         *
XC      *                X(J,I) :  =XX(J,I)                                  *
XC      *                Y(I)   :  =YY(I)                                    *
XC      *                IP(I)  :  =IND(I)                                   *
XC      *                IL(I)  :  INDICATOR OF TIES (#  OF TIES)            *
XC      *                IM(I)  :  INDICATOR OF TIES (POSITION)              *
XC      *               RINFO(I):  INFORMATION MATRIX AND ITS INVERSE        *
XC      *                          MATRIX AFTER CALLING SUBROUTINE           *
XC      *                          MATINV.                                   *
XC      *               SCORE(I):  SCORE VECTOR                              *
XC      *                                                                    *
XC      *     OUTPUT       CHI  :  GLOBAL CHI-SQUARE                         *
XC      *                 PROB  :  PROBABILITY OF CORRELATION                *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *                          SORT1, TIE, MATINV, PCHISQ                *
XC      *                                                                    *
XC      *     REFERENCE:    RUPERT G. MILLER JR., "SURVIVAL ANALYSIS", 1981, *
XC      *                         JOHN WILEY & SONS (NY:NY)                  *
XC      *                                                                    *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION RINFO(MVAR,MVAR),SCORE(MVAR),FINFO(MVAR)
X       DIMENSION IND(NTOT),IL(NTOT),IM(NTOT),IP(NTOT),Y(NTOT),YY(NTOT)
X       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT)
X       DIMENSION SWRK1(MVAR), IWRK1(MVAR),IWRK2(MVAR)
X       CHARACTER*9 OUTPUT
XC
X       DF=NVAR
XC
XC      *   SUBSTITUTE XX,YY,AND IND TO X, Y, IP TO AVOID ALTERATION OF      *
XC      *   THE ORIGINAL DATA                                                *
XC
X       DO 20 I=1,NTOT
X          IP(I)=IND(I)
X          Y(I)=YY(I)
X          IF(ICENS.EQ.-1) Y(I)=-YY(I)
XC
XC      *  IF THE OBSERVATION IS CENSORED, ADD A SMALL NUMBER TO AVOID TIES  *
XC      *  WITH DETECTED VALUE.                                              *
XC
X          IF(IP(I).NE.0) Y(I)=Y(I)*(1.0+FLOAT(ICENS)*0.0000001)
XC
X          DO 10 J=1,NVAR
X             X(J,I)=XX(J,I)
X             IF(ICENS.EQ.-1) X(J,I)=-X(J,I)
X   10     CONTINUE
X   20  CONTINUE
XC
XC      *           SORT Y IN ASCENDING ORDER                                *
XC
X       CALL SORT1(IP,X,Y,NTOT,NVAR,IL,SWRK1,MVAR)
XC 
XC      *          CHECK TIED DATA POINTS AND GIVE THEM A SPECIAL FLAG.      *
XC
X       CALL TIE(IP,X,Y,NTOT,NVAR,IL,IM,SWRK1,MVAR)
XC
XC      *         COMPUTE SCORE VECTOR. DIMENSION IS NVAR                    *
XC
X       DO 600 I=1,NVAR
X          SCORE(I)=0.0
XC
X          DO 500 J=1,NTOT 
X             IF(IP(J).EQ.0) THEN
X                IF(IL(J).EQ.1) THEN
X                   SUM=0.0
XC
X                   DO 400 K=J,NTOT
X                      SUM=SUM+X(I,K)
X  400              CONTINUE
XC
X                   JJ=J+IM(J)-1
X                   XSUM=0.0
X                   DO 450 KL=J,JJ
X                      XSUM=XSUM+X(I,KL)
X  450              CONTINUE
XC
X                   DEN=REAL(NTOT+1-J)
X                   SCORE(I)=SCORE(I)+XSUM-IM(J)*SUM/DEN
X                ENDIF
X             ENDIF
X  500     CONTINUE
X  600  CONTINUE
XC
XC      *    COMPUTE THE INFORMATION MATRIX. DIMENSION IS NVAR BY NVAR       *
XC
X       DO 1000 I=1,NVAR
X          DO  900 J=I,NVAR
X             RINFO(I,J)=0.0
XC
X             DO 800 K=1,NTOT 
X                IF(IP(K).EQ.0) THEN
X                   IF(IL(K).EQ.1) THEN
X                      SUM1=0.0
X                      SUM2=0.0
X                      SUM3=0.0
XC
X                      DO 700 L=K,NTOT
X                         SUM1=SUM1+X(I,L)
X                         SUM2=SUM2+X(J,L)
X                         SUM3=SUM3+X(I,L)*X(J,L)
X  700                 CONTINUE
X                      DEN=NTOT+1-K
X                      RINFO(I,J)=RINFO(I,J)-REAL(IM(K))
X     +                     *(SUM1*SUM2/DEN**2-SUM3/DEN)
X                   ENDIF
X                ENDIF
X  800        CONTINUE
X             RINFO(J,I)=RINFO(I,J)
X  900     CONTINUE
X 1000  CONTINUE
XC
XC      *     INVERT INFORMATION MATRX RINFO(I,J). THE INVERTED MATRIX       *
XC      *     IS STORED IN RINFO(I,J).                                       *
XC
X       CALL MATINV(RINFO,NVAR,DET,IWRK1,IWRK2,MVAR)
XC
XC      *      COMPUTE GLOBAL CHI-SQUARE:                                    *
XC      *                                                                    *
XC      *       CHI = [SCORE(I)**T] X [RINFO(I,J)**-1] X [SCORE(J)]          *
XC      *         WHERE T IS TRANSVERSE.                                     *
XC
X       DO 1200 I=1,NVAR
X          FINFO(I)=0.0
X          DO 1100 K=1,NVAR
X             FINFO(I)=FINFO(I)+RINFO(I,K)*SCORE(K)
X 1100     CONTINUE
X 1200  CONTINUE
X       CHI=0.0
XC
X       DO 1300 L=1,NVAR
X          CHI=CHI+FINFO(L)*SCORE(L)
X 1300  CONTINUE
X
XC
XC      *              GET THE REDUCED CHI-SQUARE                            *
XC
X       RCHI=CHI/DF
XC
XC      *           COMPUTE CHI-SQUARE PROBABILITY                           *
XC
X       PROB=PCHISQ(RCHI,NVAR)
XC
X       IF(OUTPUT.EQ.'         ') THEN
X          PRINT *
X          PRINT 1450
X          PRINT 1400
X          PRINT 1600,CHI
X          PRINT 1651,NVAR
X          PRINT 1650,PROB
X          PRINT *
X       ELSE
X          WRITE(60,1400)
X          WRITE(60,1450)
X          WRITE(60,1400)
X          WRITE(60,1600) CHI
X          WRITE(60,1651) NVAR
X          WRITE(60,1650) PROB
X          WRITE(60,1400)
X       ENDIF
X 1400  FORMAT('    ')
X 1450  FORMAT(5X,'CORRELATION TEST BY COX PROPORTIONAL HAZARD MODEL')
X 1600  FORMAT(6X,' GLOBAL CHI SQUARE =',F9.3,' WITH ')
X 1651  FORMAT(11X,I5,' DEGREES OF FREEDOM')
X 1650  FORMAT(6X,' PROBABILITY    =',F10.4)
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE DATA1 ******************************
XC      **********************************************************************
XC
XC      *   THIS SUBROUTINE'S PURPOSE IS TO ENABLE EASY, FOOL-PROOF          *
XC      *   KEYBOARD ENTRY OF INTEGER INPUT DATA.                            *
XC
X       SUBROUTINE DATA1(INTEG)
XC
XC      *               VARIABLE DECLARATIONS                                *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*1 B(4)
X       INTEGER*4 CUTOFF,I1,I2,INTEG,TOTAL
X       INTEGER*4 N(4)
XC
XC
XC      *    READ IN NUMBER AFTER WRITING PROMPT. FIELD SIZE = 4             *
XC
X
X 3     READ(5,1) (B(I1),I1=1,4)
X 1      FORMAT(4(A1))
X
XC
XC      *             ANALYZE DIGITS OF NUMBER                               *
XC
X       DO 2 I1=1,4
X       IF(B(I1).EQ.'0') N(I1)=0
X       IF(B(I1).EQ.'1') N(I1)=1
X       IF(B(I1).EQ.'2') N(I1)=2
X       IF(B(I1).EQ.'3') N(I1)=3
X       IF(B(I1).EQ.'4') N(I1)=4
X       IF(B(I1).EQ.'5') N(I1)=5
X       IF(B(I1).EQ.'6') N(I1)=6
X       IF(B(I1).EQ.'7') N(I1)=7
X       IF(B(I1).EQ.'8') N(I1)=8
X       IF(B(I1).EQ.'9') N(I1)=9
XC      IF((B(I1).EQ.' ').AND.(I1.EQ.1)) PRINT *,'ENTER NUMBER.'
X       IF((B(I1).EQ.' ').AND.(I1.EQ.1)) GOTO 3
XC
X       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
X     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
X     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
X     +  (B(I1).NE.'9').AND.(B(I1).EQ.' ')) CUTOFF=I1
XC
X       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
X     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
X     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
X     +  (B(I1).NE.'9').AND.(B(I1).EQ.' ')) GOTO 4
XC
X       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
X     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
X     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
X     +  (B(I1).NE.'9').AND.(B(I1).NE.' '))
X     +        PRINT *,'POSITIVE INTEGER ONLY'
XC
X       IF((B(I1).NE.'0').AND.(B(I1).NE.'1').AND.(B(I1).NE.'2').AND.
X     +  (B(I1).NE.'3').AND.(B(I1).NE.'4').AND.(B(I1).NE.'5').AND.
X     +  (B(I1).NE.'6').AND.(B(I1).NE.'7').AND.(B(I1).NE.'8').AND.
X     +  (B(I1).NE.'9').AND.(B(I1).NE.' ')) GOTO 3
X 2      CONTINUE
XC
XC      *                   TOTALS UP THE NUMBER                             *
XC
X       CUTOFF=5
X 4     TOTAL=0
XC
X       LAST=CUTOFF-1
X       DO 5 I2=1,LAST
X       TOTAL=TOTAL+N(I2)*(10**(CUTOFF-I2-1))
X 5     CONTINUE
XC
X       INTEG=TOTAL
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE DATA2  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE DATA2(B,J,IT,INTEG,LIND)
XC
XC      *   PURPOSE IS TO ENABLE EASY, FOOL-PROOF KEYBOARD ENTRY OF INTEGER  *
XC      *   INPUT DATA.                                                      *
XC      *  INPUT                                                             *
XC      *          B(4,IT)  : INPUT IN CHARACTER FORM. THIS WILL BE TESTED   *
XC      *                     AND CHANGED TO INTEGER                         *
XC      *          J        : J-TH INPUT, <=IT.                              *
XC      *          IT       : DIMENSION OF B                                 *
XC      *  OUTPUT                                                            *
XC      *          INTEG    : INTEGER OUTPUT                                 *
XC      *          LIND     : INDICATOR OF READABILITY OF INPUT              *
XC      *                     IF LIND=0, B IS SUCCESSFULLY CHANGED TO INTEG  *
XC      *                     IF LIND=1, B CANNOT BE CONVERTED TO  INTEGER   *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*1 B(4,IT)
X       INTEGER CUTOFF,TOTAL
X       INTEGER N(4)
XC
XC      *             ANALYZE DIGITS OF NUMBER                               *
XC
X        LIND=0
X       DO 2 I1=1,4
X       IF(B(I1,J).EQ.'0') GOTO 3 
X       IF(B(I1,J).EQ.'1') GOTO 3 
X       IF(B(I1,J).EQ.'2') GOTO 3 
X       IF(B(I1,J).EQ.'3') GOTO 3 
X       IF(B(I1,J).EQ.'4') GOTO 3 
X       IF(B(I1,J).EQ.'5') GOTO 3
X       IF(B(I1,J).EQ.'6') GOTO 3 
X       IF(B(I1,J).EQ.'7') GOTO 3 
X       IF(B(I1,J).EQ.'8') GOTO 3 
X       IF(B(I1,J).EQ.'9') GOTO 3 
X           IF(B(I1,J).EQ.' ') GOTO 3
X           LIND=1
X           GOTO 6
XC
X   3   IF(B(I1,J).EQ.'0') N(I1)=0
X       IF(B(I1,J).EQ.'1') N(I1)=1
X       IF(B(I1,J).EQ.'2') N(I1)=2
X       IF(B(I1,J).EQ.'3') N(I1)=3
X       IF(B(I1,J).EQ.'4') N(I1)=4
X       IF(B(I1,J).EQ.'5') N(I1)=5
X       IF(B(I1,J).EQ.'6') N(I1)=6
X       IF(B(I1,J).EQ.'7') N(I1)=7
X       IF(B(I1,J).EQ.'8') N(I1)=8
X       IF(B(I1,J).EQ.'9') N(I1)=9
XC
X       IF((B(I1,J).EQ.' ').AND.(I1.EQ.1)) LIND=1
X           IF((B(I1,J).EQ.' ').AND.(I1.EQ.1)) GOTO 6
XC
X        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
X     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.
X     +  (B(I1,J).NE.'5').AND.   (B(I1,J).NE.'6').AND.(B(I1,J).NE.'7')
X     +  .AND.(B(I1,J).NE.'8').AND.(B(I1,J).NE.'9').AND.(B(I1,J).EQ.' '))
X     +            CUTOFF=I1
XC
X        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
X     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.
X     +   (B(I1,J).NE.'5').AND.(B(I1,J).NE.'6').AND.(B(I1,J).NE.'7')
X     +   .AND.(B(I1,J).NE.'8').AND.(B(I1,J).NE.'9')
X     +   .AND.(B(I1,J).EQ.' ')) GOTO 4
XC     
X        IF((B(I1,J).NE.'0').AND.(B(I1,J).NE.'1').AND.(B(I1,J).NE.'2')
X     +  .AND.(B(I1,J).NE.'3').AND.(B(I1,J).NE.'4').AND.(B(I1,J).NE.'5')
X     +   .AND.(B(I1,J).NE.'6').AND.(B(I1,J).NE.'7').AND.(B(I1,J).NE.'8')
X     +  .AND.(B(I1,J).NE.'9').AND.(B(I1,J).NE.' ')) LIND=1   
X 2      CONTINUE
XC
XC      *                  TOTAL UP NUMBER                                   *
XC
X       CUTOFF=5
X 4      TOTAL=0
XC
X       LAST=CUTOFF-1
X       DO 5 I2=1,LAST
X       TOTAL=TOTAL+N(I2)*(10**(CUTOFF-I2-1))
X 5     CONTINUE
X       INTEG=TOTAL
X6      RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE DATREG  ****************************
XC      **********************************************************************
XC
X       SUBROUTINE DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,
X     +                   NC5,NC6,NC7,NC8,ICENS,NYC,NXC,NBC,NDIM,NDATA)
XC
XC      * THIS SUBROUTINE READS DATA FROM THE FILE FOR CORRELATION/REGRESSION*
XC      * PROBLEMS                                                           *
XC      *                                                                    *
XC      * INPUT        NVAR     : NUMBER OF THE INDEPENDENT VARIABLE         *
XC      *              FILE     : NAME OF THE DATA FILE                      *
XC      *                                                                    *
XC      * OUTPUT       IND(1,I) : INDICATOR OF CENSORSHIP                    *
XC      *              X(J,I)   : INDEPENDENT VARIABLES                      *
XC      *              Y(I)     : DEPENDENT VARIABLES                        *
XC      *              NTOT     : NUMBER OF DATA POINTS                      *
XC      *              ND       : NUMBER OF DETECTED POINTS                  *
XC      *              NC1      : NUMBER OF Y LOWER LIMITS                   *
XC      *              NC2      : NUMBER OF X LOWER LIMITS                   *
XC      *              NC3      : NUMBER OF DOUBLE LOWER LIMITS              *
XC      *              NC4      : NUMBER OF Y LOWER, X UPPER LIMITS          *
XC      *              NC5      : NUMBER OF Y UPPER LIMITS                   *
XC      *              NC6      : NUMBER OF X UPPER LIMITS                   *
XC      *              NC7      : NUMBER OF DOUBLE UPPER LIMITS              *
XC      *              NC8      : NUMBER OF Y UPPER, X LOWER LIMITS          *
XC      *              ICENS    : INDICATOR OF CENSORING                     *
XC      *              NYC      : NC1+NC5                                    *
XC      *              NXC      : NC2+NC6                                    *
XC      *              NBC      : NC3+NC4+NC7+NC8                            *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION IND(NDIM, NDATA),X(NDIM, NDATA),Y(NDATA)
X       CHARACTER*9 FILE
XC
X       OPEN(UNIT=40,FILE=FILE, STATUS='OLD', FORM='FORMATTED')
XC
X       NTOT=1
X       ND  =0
X       NC1 =0
X       NC2 =0
X       NC3 =0
X       NC4 =0
X       NC5 =0
X       NC6 =0
X       NC7 =0
X       NC8 =0
XC
X   10  READ(40,20,END=30) IND(1,NTOT),(X(J,NTOT),J=1,NVAR),Y(NTOT)
X   20  FORMAT(I4,11F10.3)
XC
X       IF(IND(1,NTOT).EQ.0)  ND =ND+1
X       IF(IND(1,NTOT).EQ.1)  NC1=NC1+1
X       IF(IND(1,NTOT).EQ.2)  NC2=NC2+1
X       IF(IND(1,NTOT).EQ.3)  NC3=NC3+1
X       IF(IND(1,NTOT).EQ.4)  NC4=NC4+1
X       IF(IND(1,NTOT).EQ.-1) NC5=NC5+1
X       IF(IND(1,NTOT).EQ.-2) NC6=NC6+1
X       IF(IND(1,NTOT).EQ.-3) NC7=NC7+1
X       IF(IND(1,NTOT).EQ.-4) NC8=NC8+1
X       NTOT=NTOT+1
X
X       IF(NTOT.GT.NDATA) THEN
X          WRITE(6,12) NDATA
X          WRITE(6,14)
X   12     FORMAT(' ****WARNING!!  THERE ARE MORE THAN',I5,' POINTS. ')
X   14     FORMAT(' ARRAYS WILL OVERFLOW UNLESS DIMENSIONS ARE CHANGED')
X          STOP       
X       ENDIF
XC
XC      * THE NEXT BLOCK OF CODE CHECKS FOR LINES WHICH HAVE ALL INPUTS     *
XC      * EQUAL TO ZERO.  THIS WOULD BE THE EFFECT OF A BLANK LINE IN THE   *
XC      * DATA FILE, AND COULD RESULT IN INCORRECT VALUES COMPUTED.         *
XC
X       K = 0
X       IF(IND(1,NTOT-1) .NE. 0) K =1
X       IF(Y(NTOT-1) .NE. 0.0) K = 1
X       DO 21 J = 1, NVAR
X          IF(X(J,NTOT-1) .NE. 0.0) K = 1
X 21    CONTINUE
X       IF(K .EQ. 0) THEN
X          WRITE(6,22)
X          WRITE(6,24) NTOT
X       ENDIF
X
X 22    FORMAT('         ')
X 24    FORMAT('WARNING: LINE ',I4,' IN THE DATA FILE MAY BE BLANK')
X
X       GOTO 10
XC
X   30  NTOT=NTOT-1
X       NYC=NC1+NC5
X       NXC=NC2+NC6
X       NBC=NC3+NC4+NC7+NC8
X       NLC=NC1+NC2+NC3+NC4
X       NUC=NC5+NC6+NC7+NC8
X       ICENS=0
X       IF(NLC.EQ.0) ICENS=-1
X       IF(NUC.EQ.0) ICENS=1
XC
X       CLOSE(UNIT=40)
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE EM  *********************************
XC      **********************************************************************
XC
X       SUBROUTINE EM(IND,X,Y,NTOT,TOL,MAXITS,NVAR,IBET,ND,NC,OUTPUT,
X     +               FILE,ALPHA,XX,Y2,W,WCEN,VCOV,WORK,SIGMAA,TOLA,
X     +               LENG,LEGWRK,MVAR)
XC
XC      *           MAXIMUM LIKELIHOOD ESTIMATION IN A LINEAR MODEL          *
XC      *           FROM CONFINED AND CENSORED NORMAL DATA                   *
XC      *                                                                    *
XC      *       FORMAL PARAMETERS :                                          *
XC      *                                                                    *
XC      *       NTOT    INPUT : THE NUMBER OF OBSERVATIONS                   *
XC      *      MPLONE   INPUT : THE TOTAL NUMBER OF PARAMETERS TO BE         *
XC      *                       ESTIMATED (I.E. NB+1)                        *
XC      *      MAXITS   INPUT : THE MAXIMUM NUMBER OF ITERATIONS ALLOWED     *
XC      *      IBET     INPUT : IF THERE IS DATA WHICH IS CONFINED           *
XC      *                       BETWEEN TWO VALUES, IBET=1 OTHERWISE 0       *
XC      *      ALPHA    INPUT : IF ALPHA(NVAR+2).LE.0.0, THE SUBROUTINE      *
XC      *                       CALCULATES INITIAL PARAMETER ESTIMATES.      *
XC      *                       IF >0.0, IT CONTAINS THE INITIAL             *
XC      *                       ESTIMATE OF THE J-TH LOCATION PARAMETER      *
XC      *                       FOR J=1,2,.....,NB.                          *
XC      *               OUTPUT: THE MOST RECENT PARAMETER ESTIMATES          *
XC      *                       BEFORE EXIT FROM THE SUBROUTINE.             *
XC      *      TOL      INPUT : CONVERGENCE TO THE MAXIMUM LIKELIHOOD        *
XC      *                       PARAMETER ESTIMATE IS REACHED IF THE         *
XC      *                       DIFFERENCE BETWEEN CONSECUTIVE ESTIMATES     *
XC      *                       OF THE J-TH PARAMETER IS LESS THAN TOL(J)    *
XC      *                       FOR J=1,2,........,M.                        *
XC      *       IND     INPUT : INDICATOR OF CENSORED DATA                   *
XC      *       XX      INPUT : THE DESIGN MATRIX XX(I,J) CONTAINS THE       *
XC      *                       COEFFICIENT OF THE J-TH LOCATION PARA-       *
XC      *                       METER FOR I-TH OBSERVATION.                  *
XC      *       Y       INPUT : IF PP(I)=0, THE I-TH OBSERVATION IS          *
XC      *                       COMPLETELY SPECIFIED IN Y(I); IF             *
XC      *                       PP(I)=-1, Y(I)  IS LEFT-CENSORED: IF         *
XC      *                       PP(I)=1, Y(I) IS RIGHT-CENSORED: IF          *
XC      *                       PP(I)=5, Y(I) IS CONFINED BETWEEN TWO        *
XC      *                       VALUES.                                      *
XC      *      ROWX     WORK  : THE NUMBER OF ROWS OF X (=NTOT)              *
XC      *      COLX     WORK  : THE NUMBER OF COLUMNS OF X (=NB)             *
XC      *       W       WORK  :                                              *
XC      *      LENW     WORK  :   = NB+NTOT                                  *
XC      *      WORK     WORK  :                                              *
XC      *      LENWRK   WORK  :   = NB*NTOT                                  *
XC      *      LEG            : DIMENSION SIZE >= LENW                       *
XC      *      LEGWRK         : DIMENSION SIZE >= LWNWRK                     *
XC      *      MVAR           : DIMENSION SIZE >= MPLONE                     *
XC      *      VCOV     OUTPUT: IF THE PROCEDURE CONVERGED TO THE MAX-       *
XC      *                       IMUM LIKELIHOOD ESTIMATES, THE FIRST         *
XC      *                       (NB+1)*(NB+1) POSITIONS CONTAIN AN ESTI-     *
XC      *                       MATE OF THE VARIANCE-COVARIANCE MATRIX       *
XC      *                       OF THESE ESTIMATES.                          *
XC      *      SIGMAA   OUTPUT: (J,J) COMPONENT CONTAINS THE STANDARD        *
XC      *                       DEVIATION OF THE J-TH PARAMETER.             *
XC      *      IFAULT   OUTPUT: FAILURE INDICATOR                            *
XC      *      ICHECK   OUTPUT: IF ERROR ANALYSIS IS NOT COMPLETED,          *
XC      *                       ICHECK HAS A VALUE : 1 ,OTHERWISE 0 .        *
XC      *                                                                    *
XC      *       SUBROUTINES                                                  *
XC      *                        EMALGO                                      *
XC      *                                                                    *
XC      *       REF : M.S.WOLYNETZ AS 139 APL.STATIST.VOL.28 195 (1979)      *
XC      *             PLUS CORRECTIONS IN LATER ISSUES OF APPLIED STATISTICS.*
XC      *                                                                    *
XC      *       SUBROUTINE EMALGO IS COPIED DIRECTLY FROM WOLYNETZ, EXCEPT   *
XC      *       FOR A FEW CHANGES TO PRINT AND COMMENT STATEMENTS.           *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
X       INTEGER   ROWX,COLX
X
X       DIMENSION X(MVAR,NTOT),Y(NTOT),IND(NTOT),XX(NTOT,MVAR),Y2(NTOT)
X       DIMENSION W(LENG),WCEN(LENG),VCOV(LEGWRK),WORK(LEGWRK)
X       DIMENSION ALPHA(MVAR),SIGMAA(MVAR,MVAR),TOLA(MVAR)
X       CHARACTER*9 FILE,OUTPUT
XC
XC
XC      *                  INITIALIZATION                                    *
XC
X       MPLONE=NVAR+2
X       NB=NVAR+1
X       ROWX=NTOT
X       COLX=MPLONE
X       LENW=NTOT+MPLONE
X       LENWRK=NTOT*MPLONE
XC
X       DO 10 I=1,MPLONE
X          TOLA(I)=TOL
X   10  CONTINUE
XC
XC      *                                                                    *
XC      *   IF THERE ARE DATA POINTS CONFINED BETWEEN TWO VALUES, READ THE   *
XC      *   INPUT DATA AGAIN.                                                *
XC      *                                                                    *
X       IF(IBET.NE.0) THEN
X
X          OPEN(UNIT=40,FILE=FILE,STATUS='OLD',FORM='FORMATTED')
X
X          DO 40 I=1,NTOT
X             READ(40,30) IND(I),(XX(I,J),J=2,NB),Y(I),Y2(I)
X   30        FORMAT(I3,20F10.3)
X             XX(I,1) = 1.0
X   40     CONTINUE
X       CLOSE(UNIT=40)
XC
X       ELSE
X          DO 100 I=1,NTOT
X             XX(I,1)=1.0
X             DO 90 J=1,NVAR
X                JJ=J+1
X                XX(I,JJ)=X(J,I)
X   90        CONTINUE
X  100     CONTINUE
X       ENDIF
XC
X       ICHECK=0
XC
XC      *   CALL THE SUBPROGRAM FOR THE MAXIMUM LIKELIHOOD ESTIMATES         *
XC
X       CALL EMALGO(NTOT,Y,Y2,IND,MPLONE,XX,ROWX,COLX,W,WCEN,LENW,
X     +      VCOV,WORK,LENWRK,ALPHA,TOLA,MAXITS,IFAULT,ICHECK,NC)
XC
XC
XC      *    PRINT OUT THE RESULTS                                           *
XC
XC
X       IF(IFAULT.GT.0) THEN
X          IF(ICHECK.EQ.0) THEN
XC
XC      *        CHANGE VCOV ARRAY FROM ONE DIM. TO TWO DIM.                 *
XC
X             K=1
X             DO 300 I=1,MPLONE
X                DO 200 J=1,MPLONE
X                   SIGMAA(I,J)=DSQRT(DABS(VCOV(K)))
X                   K=K+1
X  200           CONTINUE
X  300        CONTINUE
X          ENDIF
X
X  360     IF(OUTPUT.EQ.'         ') THEN
XC
XC      PRINT FINAL REGRESSION COEFFICIENTS AT THE TERMINAL. 
XC
X             PRINT 1050
X             PRINT 1020         
X             PRINT 1050
X             PRINT 1200,ALPHA(1),SIGMAA(1,1)
X
X             DO 450 J=2,NB
X                JI=J-1
X                PRINT 1250,JI,ALPHA(J),SIGMAA(J,J)
X  450        CONTINUE
X
X             PRINT 1300,ALPHA(MPLONE)
X             IF(ICHECK.EQ.0) THEN
X                ITE=IFAULT
X             ELSE
X                ITE=ICHECK
X             ENDIF
X             PRINT 1350,ITE
X             PRINT 1050
X          ELSE
XC
XC      PRINT FINAL REGRESSION COEFFICIENTS IN THE OUTPUT FILE.
XC
X             WRITE(60,1050)
X             WRITE(60,1020)        
X             WRITE(60,1050)
X             WRITE(60,1200) ALPHA(1),SIGMAA(1,1)
X   
X             DO 455 J=2,NB
X                JI=J-1
X                WRITE(60,1250) JI,ALPHA(J),SIGMAA(J,J)
X  455        CONTINUE
X             WRITE(60,1300) ALPHA(MPLONE)
X             IF(ICHECK.EQ.0) THEN
X                ITE=IFAULT
X             ELSE
X                ITE=ICHECK
X             ENDIF
X             WRITE(60,1350) ITE
X             WRITE(60,1050)
X          ENDIF
XC
XC      IF AN ERROR OCCURED DURING THE COMPUTATION, PRINT OUT AN ERROR
XC      MESSAGE.
XC
XC      IN THE FOLLOWING, WE HAVE COLLECTED ALL OF WOLYNETZ'S ERROR FLAGS
XC      AND PRINT OUT THE APPROPRIATE ERROR MESSAGE. 
XC
X
X       ELSE
X
XC
X          IF(OUTPUT.EQ.'         ') THEN
XC
XC      PRINT FINAL REGRESSION COEFFICIENTS ON THE TERMINAL. 
XC
X             PRINT 1000
X             PRINT 1020
X             PRINT 1050
X             PRINT 2001
X             IF(IFAULT.EQ.-1) THEN
X                PRINT 2002
X                PRINT 2003
X             ELSEIF(IFAULT.EQ.-4) THEN
X                PRINT 2012
X                PRINT 2014
X             ELSEIF(IFAULT.EQ.-5) THEN
X                PRINT 2021
X                PRINT 2022
X                PRINT 2023
X                PRINT 2024
X                PRINT 2025
X                PRINT 2026
X                PRINT 2027
X             ELSEIF(IFAULT.EQ.-6) THEN
X                PRINT 2031
X                PRINT 2032
X                PRINT 2033
X                PRINT 2034
X                PRINT 2035
X                PRINT 2036
X                PRINT 2037
X             ELSEIF(IFAULT.EQ.-7) THEN
X                PRINT 2045
X             ELSEIF(IFAULT.EQ.-8) THEN
X                PRINT 2055
X             ELSE
X                PRINT 2070
X             ENDIF
XC
X       ELSE
XC
XC      PRINT FINAL REGRESSION COEFFICIENTS IN THE OUTPUT FILE.
XC
X             WRITE(60,1000)
X             WRITE(60,1020)
X             WRITE(60,1050)
X             WRITE(60,2001)
X             IF(IFAULT.EQ.-1) THEN
X                WRITE(60,2002)
X                WRITE(60,2003)
X             ELSEIF(IFAULT.EQ.-4) THEN
X                WRITE(60,2012)
X                WRITE(60,2014)
X             ELSEIF(IFAULT.EQ.-5) THEN
X                WRITE(60,2021)
X                WRITE(60,2022)
X                WRITE(60,2023)
X                WRITE(60,2024)
X                WRITE(60,2025)
X                WRITE(60,2026)
X                WRITE(60,2027)
X             ELSEIF(IFAULT.EQ.-6) THEN
X                WRITE(60,2031)
X                WRITE(60,2032)
X                WRITE(60,2033)
X                WRITE(60,2034)
X                WRITE(60,2035)
X                WRITE(60,2036)
X                WRITE(60,2037)
X             ELSEIF(IFAULT.EQ.-7) THEN
X                WRITE(60,2045)
X             ELSEIF(IFAULT.EQ.-8) THEN
X                WRITE(60,2055)
X             ELSE
X                WRITE(60,2070)
X             ENDIF
X          ENDIF
XC
X       ENDIF
XC
X 1000  FORMAT(T10,'REGRESSION ANALYSIS WITH CENSORED DATA')
X 1020  FORMAT(T5,'LINEAR REGRESSION BY PARAMETRIC EM ALGORITHM')
X 1050  FORMAT(T5,'      ')
X 1100  FORMAT(T8,'DATA TITLE :',T25,60A1)
X 1150  FORMAT(T8,'TOTAL # OF DATA :',T26,I3,T33,'CENSORED DATA :'
X     +           ,T48,I3)
X 1200  FORMAT(T8,'INTERCEPT COEFF    :',F8.4,T38,'+/-',T41,F8.4)
X 1250  FORMAT(T8,'SLOPE COEFF ',I1,'      :',F8.4,T38,'+/-',T41,
X     +         F8.4,5X)
X 1300  FORMAT(T8,'STANDARD DEVIATION :',F8.4)
X 1350  FORMAT(T8,'ITERATIONS         :',I3)
X 2001  FORMAT(T5,'NOTICE :')
X 2002  FORMAT(T5,'MAXIMUM NUMBER OF ITERATION REACHED AND')
X 2003  FORMAT(T5,'CONVERGENCE HAS NOT BEEN OBTAINED.')
X 2012  FORMAT(T5,'NUMBER OF COMPLETELY SPECIFIED DATA IS LESS')
X 2014  FORMAT(T5,'THAN NB+1')
X 2021  FORMAT(T5,'THE MATRIX IS NOT POSITIVE DEFINITE,AS')
X 2022  FORMAT(T5,'DETERMINED BY SUBROUTINE "SYMINV",A MATRIX')
X 2023  FORMAT(T5,'INVERSION PROCEDURE(HEALY,1968); THE VALUE')
X 2024  FORMAT(T5,'OF "NULLTY" AND "IFAULT" RETURNED BY SYMINV')
X 2025  FORMAT(T5,'ARE PLACED IN THE FIRST TWO POSITIONS OF ')
X 2026  FORMAT(T5,'THE ARRAY "VCOV" BEFORE RETURNING TO THE ')
X 2027  FORMAT(T5,'CALLING PROGRAM.')
X 2031  FORMAT(T5,'THE ESTIMATE OF THE VARIANCE- COVARIANCE')
X 2032  FORMAT(T5,'MATRIX IS NOT POSITIVE DEFINITE, AS DETER-')
X 2033  FORMAT(T5,'MINED BY SUBROUTINE "SYMINV"(HEALY,1968):')
X 2034  FORMAT(T5,'THE VALUES OF "NULLTY" AND "IFAULT", RETURNED')
X 2035  FORMAT(T5,'BY SYMINV ,ARE PLACED IN THE FIRST TWO ')
X 2036  FORMAT(T5,'POSITIONS OF THE ARRY "VCOV" BEFORE RETURNING')
X 2037  FORMAT(T5,'TO THE CALLING PROGRAM')
X 2045  FORMAT(T5,'"ROWX" IS LESS THAN NTOT')
X 2055  FORMAT(T5,'"COLX" IS LESS THAN NB')
X 2070  FORMAT(T5,'THE PROGRAM IS STOPPED FOR UNKNOWN REASONS')
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      *********************** FUNCTION FACTOR  *****************************
XC      **********************************************************************
XC
X       FUNCTION FACTOR(N)
XC
XC      *     COMPUTES THE FACTORIAL FUNCTION FOR INTEGERS.                  *
XC      *     THIS FUNCTION IS BASED ON PROG. 3-2 ON P. 32 OF "DATA REDUCTION*
XC      *     AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", PHILIP R.       *
XC      *     BEVINGTON, 1969, McGRAW HILL (NY:NY)                           *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
XC
X       FACTOR=1.0
X       IF(N-1 .GT. 0) THEN
X   13     IF(N-10 .LE. 0) THEN
XC
XC      *                N LESS THAN 11                                      *
XC
X   21        DO 23 I=2,N
X                FI=I
X                FACTOR=FACTOR*FI
X   23        CONTINUE
XC
XC      *               N GREATER THAN 10                                    *
XC
X          ELSE
X   31        SUM=0.0
X             DO 34 I=11,N
X                FI=I
X                SUM=SUM+DLOG(FI)
X   34        CONTINUE
X             FACTOR=3628800.0D00*DEXP(SUM)
X          ENDIF
X       ENDIF
X   40  RETURN
X       END
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE DATAIN  ****************************
XC      **********************************************************************
XC
X       SUBROUTINE DATAIN(IUNI,FILE,NVAR,ISTA,IND,X,NTOT,NDATA,MVAR)
XC
XC      *  THIS SUBROUTINE READS DATA FOR UNIVARIATE AND TWO SAMPLE PROBLEMS *
XC      *                                                                    *
XC      *  INPUT                                                             *
XC      *          INUI  : INDICATOR OF WHICH PROBLEM (KM OR TWO SAMPLE) IS  *
XC      *                  NEEDED                                            *
XC      *          FILE  : DATA FILE NAME                                    *
XC      *          NVAR  : NUMBER OF VARIABLES                               *
XC      *          NDATA : DIMENSION FOR THE VARIABLES                       *
XC      *  OUTPUT                                                            *
XC      *         ISTA(I): INDICATOR OF GROUPS                               *
XC      *        IND(J,I): INDICATOR OF CENSORSHIP                           *
XC      *          X(J,I): VARIABLES                                         *
XC      *         NTOT   : NUMBER OF DATA POINTS                             *
XC      *                                                                    *
XC      *    THIS FILE WAS MODIFIED ON 4/13/90                               *
XC      *       NDATA WAS ADDED FOR THE DIMENSION DECLARATION.              *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION IND(MVAR, NDATA), X(MVAR,NDATA), ISTA(NDATA)
X       CHARACTER*9 FILE
XC
XC      *        OPEN THE DATA FILE AND READ DATA                            *
XC
X       OPEN(UNIT=40, FILE=FILE, STATUS='OLD', FORM='FORMATTED')
XC
X       J=0
X       IF(IUNI .EQ. 1) THEN
XC
XC      *            K-M ESTIMATOR FORMAT                                    *
XC
X   10     J=J+1
X          READ(40,30,END=50) (IND(I,J),X(I,J),I=1,NVAR)
X          K=0
X          DO 15 I = 1,NVAR
X             IF(IND(I,J).NE.0.OR.X(I,J).NE.0.0) K = 1
X 15       CONTINUE
X          IF(K .EQ. 0) THEN
X             WRITE(6,44)
X             WRITE(6,45) J
X          ENDIF
X          GOTO 10
XC
XC      *            TWO SAMPLE TEST FORMAT                                  *
XC
X       ELSE
X   20     J=J+1
X          READ(40,40,END=50) ISTA(J),(IND(I,J),X(I,J),I=1,NVAR)
X          K = 0
X          IF(ISTA(J) .NE. 0) K=1
X          DO 25 I = 1, NVAR
X             IF(IND(I,J) .NE. 0 .OR. X(I,J) .NE. 0.0) K = 1
X 25       CONTINUE
X          IF(K .EQ. 0) THEN
X             WRITE(6,44)
X             WRITE(6,45) J
X          ENDIF
X          GOTO 20
X       ENDIF
XC
X   30  FORMAT(10(I4,F10.3))
X   40  FORMAT(I4,10(I4,F10.3))
X   44  FORMAT('         ')  
X   45  FORMAT('WARNING: LINE ',I4,
X     +         ' IN THE DATA FILE MAY BE BLANK') 
X   50  NTOT=J-1
X       CLOSE(UNIT=40)
X     
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE EMALGO  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE EMALGO(NTOT,Y1,Y2,P,MPLONE,X,ROWX,COLX,W,WCEN,LENW,
X     +          VCOV,WORK,LENWRK,ALPHA,TOL,MAXITS,IFAULT,ICHECK,NC)
XC
XC
XC      *       ALGORITHM AS 139 APPL.STATIST. (1979) VOL.28., NO2           *
XC      *                                                                    *
XC      *       COMPUTES MAXIMUM LIKELIHOOD ESTIMATES                        *
XC      *       FROM A LINEAR MODEL WITH NORMAL HETEROGENEOUS                *
XC      *       VARIANCE. THE DESIGN MATRIX MUST BE NON-SINGULAR.            *
XC      *       THE DEPENDENT VARIABLE MAY INCLUDE OBSERVATIONS              *
XC      *       CENSORED IN EITHER TAIL AND/OR OBSERVATIONS CONFINED         *
XC      *       BETWEEN FINITE LIMITS.                                       *
XC      *                                                                    *
XC      *  SUBROUTINES                                                       *
XC      *             SYMINV, UNPACK, RMILLS                                 *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       INTEGER ROWX,COLX,P
X
X       DIMENSION X(NTOT,MPLONE),TOL(MPLONE),ALPHA(MPLONE)
X       DIMENSION Y1(NTOT),Y2(NTOT),P(NTOT)
X       DIMENSION VCOV(LENWRK),WORK(LENWRK),W(LENW),WCEN(LENW)
X       DATA QLIMIT /0.00001/, RLIMIT /0.00001/
X       DATA C /0.39894228/
XC
XC      *                    INITIALIZATION                                  *
XC      *       THE NEXT LINE IS LOCATED IN A DIFFERENT PLACE IN THE         *
XC      *       ORIGINAL PROGRAM BY WOLYNETZ                                 *
XC
X       M=MPLONE-1
XC
XC      *                 CHECK ARRAY SIZES, ETC                             *
XC
X       IFAULT=-7
X       IF(ROWX.LT.NTOT) RETURN
X       IFAULT=-8
X       IF(COLX.LT.M) RETURN
X       IFAULT=-9
XC
XC
XC      *       THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 31       *
XC      *       APPL.STAT. VOL.29 P.228 , (1980).                            *
XC
XC
X       IF(LENW.LT.(MPLONE+NTOT)) RETURN
X       IFAULT=-10
X       IF(LENWRK.LT.(M*NTOT)) RETURN
XC
XC      *          COMPUTE X'X IN LOWER TRIANGULAR FORM                      *
XC
X       II=0
X       DO 53 I=1,M
X          DO 50 J=1,I
X             TEMP=0.0
X             DO 40 K=1,NTOT
X                TEMP=TEMP+X(K,I)*X(K,J)
X   40        CONTINUE
X             II=II+1
X             VCOV(II)=TEMP
X   50     CONTINUE
X   53  CONTINUE
X
X       CALL SYMINV(VCOV,M,WORK,W,NUL,LENWRK,LENWRK,LENW,IFAULT)
X
X       IF((IFAULT.NE.0) .OR. (NUL.NE.0)) THEN
X          VCOV(2)=REAL(IFAULT)
X          VCOV(1)=REAL(NUL)
X          IFAULT=-5
X          RETURN
X       ENDIF
X
XC
XC      *      THE MATRIX IS NON-SINGULAR AND WE HAVE ITS INVERSE.  NOW WE   *
XC      *      COMPUTE (X'X) INVERSE*X.                                      *
XC      *      THE FOLLOWING SCHEME IS USED TO REDUCE THE NUMBER OF STORAGE  *
XC      *      ARRAYS NEEDED, BY EXPANDING FROM THE TRIANGULAR TO A SQUARE   *
XC      *      MATRIX.                                                       *
XC
XC      *         THE NEXT LINE IS NOT IN WOLYNETZ.                          *
XC
X       IF(M.NE.1) THEN
X          CALL UNPACK(WORK,M,LENWRK)
X       ENDIF
XC
XC      *      DO MULTIPLICATION--ONE ROW AT A TIME--STARTING WITH           *
XC      *      THE LAST ROW                                                  *
XC
X       JJ=NTOT*M
X       II=M*M
X       DO 220 I=1,M
X          II=II-M
X
X          DO 200 J=1,NTOT
X             TEMP=0.0
X
X             DO 170 K=1,M
X                IIK=II+K
X                TEMP=TEMP+WORK(IIK)*X(J,K)
X  170        CONTINUE
X             W(J)=TEMP
X  200     CONTINUE
X
X          DO 210 J=1,NTOT
X             IJ=NTOT+1-J
X             WORK(JJ)=W(IJ)
X             JJ=JJ-1
X  210     CONTINUE
X  220  CONTINUE
XC
X       XSIG=ALPHA(MPLONE)
X       IF(XSIG.LE.0.0) THEN
XC
XC      *       NO ACCEPTABLE INITIAL VALUE FOR SIGMA HAS BEEN INPUT,        *
XC      *       OBTAIN INITIAL ESTIMATES FROM EXACTLY SPECIFIED              *
XC      *       OBSERVATIONS ONLY (ALTHOUGH THE MATRIX IS BASED ON ALL       *
XC      *       OBSERVATIONS) AND CONFINED OBSERVATIONS.                     *
XC
X          II=-NTOT
X          DO 300 I=1,M
X             II=II+NTOT
X             TEMP=0.0
X             DO 280 J=1,NTOT
X                IIJ=II+J
X                IPT=P(J)
X                IF(IPT.EQ.0) THEN
X                   TEMP=TEMP+WORK(IIJ)*Y1(J)
X                ELSEIF(IPT.EQ.5) THEN 
X                   TEMP=TEMP+WORK(IIJ)*(Y1(J)+Y2(J))*0.5
X                ENDIF
X  280        CONTINUE
X             ALPHA(I)=TEMP
X  300     CONTINUE
XC
XC      *           CALCULATE THE INITIAL ESTIMATE OF SIGMA                  *
XC
X          SUM2=0.0
X          TEMP=0.0
X          DO 350 I=1,NTOT
X             IPT=P(I)
X             IF(IABS(IPT).NE.1) THEN
X                DEMP=Y1(I)
X                IF(IPT.EQ.5) DEMP=(DEMP+Y2(I))*0.5
X
X                DO 320 J=1,M
X                   DEMP=DEMP-ALPHA(J)*X(I,J)
X  320           CONTINUE
X
X                SUM2=SUM2+DEMP**2
X                TEMP=TEMP+1.0
X             ENDIF
X  350     CONTINUE
X
X          XSIG=DSQRT(SUM2/TEMP)
X       ENDIF
XC
XC      *         COMPUTE SOME CONSTANTS NEEDED THROUGHOUT THE SUBROUTINE    *
XC
X       R=0.0
X       R2=0.0
X       IFAULT=-2
X       DO 600 I=1,NTOT
X          IPT=P(I)
X          IF(IPT.EQ.0)  THEN
X             R=R+1.0
X             W(I)=Y1(I)
XC
XC      *       THE NEXT LINE IS LOCATED IN A DIFFERENT PLACE IN THE         *
XC      *       ORIGINAL PROGRAM BY WOLYNETZ                                 *
XC
X          ELSEIF(IPT.EQ.5) THEN
X             IF(DABS(Y1(I)-Y2(I)) .GT. QLIMIT*DABS(Y1(I))) THEN
X                R2=R2+1.0
X                IF(Y1(I).LT.Y2(I)) GOTO 600
X                RETURN
X             ENDIF
X             P(I)=0
X             R=R+1.0
X             W(I)=Y1(I)
X          ENDIF
X  600  CONTINUE
X
X       I=INT(R+R2+0.01)
X       IFAULT=-4
X       IF(I.LT.MPLONE) RETURN
X       IFAULT=0
XC
XC      *             START OF THE ITERATION PROCEDURE                       *
XC
X  620  TD=R
X       SUM2=0.0
XC
XC      *             COMPLETE W-VECTOR                                      *
XC
X       DO 1000 I=1,NTOT
X          IPT=P(I)
X          YMEAN=0.0
X
X          DO 650 J=1,M
X             YMEAN=YMEAN+ALPHA(J)*X(I,J)
X  650     CONTINUE
XC
XC      *       OBSERVATION IS NOT EXACTLY SPECIFIED: START FROM LINE 990    *
XC
X          IF(IPT.NE.0) THEN
X             TEMP=(Y1(I)-YMEAN)/XSIG
XC
XC      *        OBSERVATION CENSORED FROM ABOVE:  LOWER BOUND IS KNOWN      *
XC
X             IF((IPT-1) .EQ. 0) THEN
X  700           CALL RMILLS(TEMP,F,RLIMIT)
X                W(I)=YMEAN+XSIG*F
X                TD=TD+F*(F-TEMP)
XC
XC      *         OBSERVATON CENSORED FROM BELOW:  UPPER BOUND IS KNOWN      *
XC
X             ELSEIF((IPT-1) .LT. 0) THEN
X  750           CALL RMILLS(-TEMP,F,RLIMIT)
X                W(I)=YMEAN-XSIG*F
X                TD=TD+F*(F+TEMP)
XC
XC      *       OBSERVATION CONFINED TO LIE BETWEEN TWO FINITE LIMITS        *
XC
X             ELSEIF((IPT-1) .GT. 0) THEN
X  800           YN=DEXP(-0.5*TEMP**2)*C
X                CALL RMILLS(TEMP,F,RLIMIT)
X                YQ=YN/F
X                TMPU=(Y2(I)-YMEAN)/XSIG
X                YNU=DEXP(-0.5*TMPU**2)*C
X                CALL RMILLS(TMPU,FU,RLIMIT)
X                YQU=YNU/FU
X                TINT=YQ-YQU
X
X                IF(TINT.LT.QLIMIT) THEN
X                   IFAULT=-3
X                   RETURN
X                ENDIF
XC
XC      *       AFTER STANDARDIZING, UPPER AND LOWER LIMITS RESULT IN        *
XC      *       THE SAME PROBABILITY INTEGRAL                                *
XC
X  820           A=(YN-YNU)/TINT
X                W(I)=YMEAN+XSIG*A
X                TD=TD+(A**2+(TMPU*YNU-TEMP*YN)/TINT)
X             ENDIF
X          ENDIF
XC
XC      *        CALCULATE RESIDUAL SUM OF SQUARES                           *
XC
X  990     SUM2=SUM2+(W(I)-YMEAN)**2
X 1000  CONTINUE
XC
XC      *    UPDATE PARAMETER ESTIMATES-STORE IN THE END OF THE W-VECTOR     *
XC
X       JJ=-NTOT
X       DO 1200 J=1,M
X          JJ=JJ+NTOT
X          TEMP=0.0
X
X          DO 1100 I=1,NTOT
X             JJI=JJ+I
X             TEMP=TEMP+WORK(JJI)*W(I)
X 1100     CONTINUE
X          NJ=NTOT+J
X          W(NJ)=TEMP
X 1200  CONTINUE
X
X       NJ=NTOT+MPLONE
X       W(NJ)=DSQRT(SUM2/TD)
XC
XC      *       STORE THE ESTIMATES FOR THE CENSORED POINTS                  *
XC
X       KC=0
X       DO 1250 I=1,NTOT
X          IF(P(I).EQ.0) GOTO 1250
X          KC=KC+1
X          WCEN(KC)=W(I)
X 1250  CONTINUE
XC
XC      *             TEST FOR CONVERGENCE                                   *
XC
X       DO 1300 J=1,MPLONE
X          NJ=NTOT+J
X          IF(DABS(ALPHA(J)-W(NJ)).GE.TOL(J)) GOTO 1400
X 1300  CONTINUE
XC
XC      *          IF WE REACH HERE, CONVERGENCE IS OBTAINED                 *
XC
X       IJ=IFAULT
X       IFAULT=-1
XC
XC      *                UPDATE VALUES                                       *
XC
X 1400  DO 1450 J=1,MPLONE
X          NJ=NTOT+J
X          ALPHA(J)=W(NJ)
X 1450  CONTINUE
X
X       XSIG=ALPHA(MPLONE)
X       IFAULT=IFAULT+1
X       IF(IFAULT.NE.0) THEN
XC
XC      *      IF THE NUMBER OF ITERATIONS HAS NOT REACHED THE MAX., TRY     *
XC      *      ANOTHER ITERATION.                                            *
XC
X          IF(IFAULT.LE.MAXITS) GOTO 620
X          IFAULT=-1
X          RETURN
X       ENDIF
XC
XC      *        CONVERGENCE OBTAINED.  COMPUTE VARIANCE--COVARIANCE         *
XC      *        MATRIX, AND INITIALIZE THE WORK ARRAY                       *
XC
X 1600  II=MPLONE*(MPLONE+1)/2
X       DO 1650 I=1,II
X          WORK(I)=0.0
X 1650  CONTINUE
X
X       DO 2500 I=1,NTOT
X          IPT=P(I)
X          YS=Y1(I)
X
X          DO 1680 J=1,M
X             YS=YS-ALPHA(J)*X(I,J)
X 1680     CONTINUE
X
X          YS=YS/XSIG
X          JJ=0
X          IF(IPT.EQ.0) THEN
X
XC
XC      *             EXACTLY SPECIFIED OBSERVATION                          *
XC
X
X             DO 1750 K=1,M
X                DO 1720 J=1,K
X                   JJ=JJ+1
X                   WORK(JJ)=WORK(JJ)+X(I,K)*X(I,J)
X 1720           CONTINUE
XC
XC
XC      *      THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32        *
XC      *      APPL.STAT. VOL 30, P. 105 (1981).                             *
XC
XC
X                KK=II-MPLONE+K
X                WORK(KK)=WORK(KK)+YS*X(I,K)
X 1750        CONTINUE
X             WORK(II)=WORK(II)+1.0+YS**2
X
XC
XC      *      OBSERVATION CENSORED FROM ABOVE:  LOWER BOUND IS KNOWN        *
XC
X
X          ELSEIF((IPT-1) .LE. 0) THEN
X             IF((IPT-1) .EQ. 0) THEN
X                CALL RMILLS(YS,F,RLIMIT)
X                TEMP=F*(F-YS)
X
XC
XC      *      OBSERVATION CENSORED FROM BELOW:  UPPER BOUND IS KNOWN        *
XC
X
X             ELSE
X                CALL RMILLS(-YS,F,RLIMIT)
X                TEMP=F*(F+YS)
X             ENDIF
XC
XC      *         ROUTINE FOR CENSORED OBSERVATIONS                          *
XC
X             DO 2190 K=1,M
X                DO 2170 J=1,K
X                   JJ=JJ+1
X                   WORK(JJ)=WORK(JJ)+X(I,J)*X(I,K)*TEMP
X 2170           CONTINUE
XC
XC      *      THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32        *
XC      *      APPL.STAT. VOL 30, P. 105 (1981).                             *
XC
X                KK=II-MPLONE+K
X                WORK(KK)=WORK(KK)+YS*X(I,K)*TEMP
X 2190        CONTINUE
X             WORK(II)=WORK(II)+YS**2*TEMP
X
XC
XC      *       OBSERVATION CONFINED BETWEEN TWO FINITE LIMITS               *
XC
X
X          ELSEIF((IPT-1) .GT. 0) THEN
X             YN=DEXP(-0.5*YS**2)*C
X             CALL RMILLS(YS,F,RLIMIT)
X             YQ=YN/F
X             YSU=YS+(Y2(I)-Y1(I))/XSIG
X             CALL RMILLS(YSU,FU,RLIMIT)
X             YNU=DEXP(-0.5*YSU**2)*C
X             YQU=YNU/FU
X             TINT=YQ-YQU
X             A=(YN-YNU)/TINT
X             B=(YNU*YSU-YN*YS)/TINT
X             TEMP=A**2+B
X
X             DO 2350 K=1,M
X
X                DO 2330 J=1,K
X                   JJ=JJ+1
X                   WORK(JJ)=WORK(JJ)+X(I,J)*X(I,K)*TEMP
X 2330           CONTINUE
X                TEMP=(YS**2*YN-YSU**2*YNU)/TINT
XC
XC
XC      *     THE NEXT LINE IS CORRECTED ACCORDING TO REMARK AS R 32         *
XC      *     APPL.STAT. VOL 30, P. 105 (1981).                              *
XC
XC
X                KK=II-MPLONE+K
X                WORK(KK)=WORK(KK)-(TEMP+A*B)*X(I,K)
X 2350        CONTINUE
X
X             TEMP=(YS**3*YN-YSU**3*YNU)/TINT
X             WORK(II)=WORK(II)-TEMP+B**2
X           ENDIF
X 2500  CONTINUE
XC
XC      *                   INVERT THE MATRIX                                *
XC
X       CALL SYMINV(WORK,MPLONE,VCOV,W,NUL,LENWRK,LENWRK,LENW,IFAULT)
X
X       IF((IFAULT.NE.0).OR.(NUL.NE.0)) THEN
X          VCOV(2)=REAL(IFAULT)
X          VCOV(1)=REAL(NUL)
X          IFAULT=-6
X          ICHECK=IJ
X          RETURN
X       ENDIF
XC
XC      *               RESTORE THE ITERATION COUNTER                         *
XC
X       IFAULT=IJ
XC
XC      *               MULTIPLY BY SIGMA-SQUARED                            *
XC
X       TEMP=XSIG**2
X       DO 2580 I=1,II
X          VCOV(I)=VCOV(I)*TEMP
X 2580  CONTINUE
X
XC
XC      *                UNPACK THE MATRIX                                   *
XC
X       CALL UNPACK(VCOV,MPLONE,LENWRK)
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* FUNCTION GAMMA  ********************************
XC      **********************************************************************
XC
X       FUNCTION GAMMA(X)
XC
XC      *     THIS FUNCTION IS OBTAINED FROM PHILIP R. BEVINGTON, "DATA      *
XC      *     REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", 1969, *
XC      *     McGRAW HILL (NY:NY), PROGRAM 7-2 P. 126                        *
XC      *     THIS COMPUTES THE GAMMA FUNCTION FOR INTEGERS AND HALF-INTEGERS*
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
XC
XC      *                INTEGERIZE ARGUMENT                                 *
XC
X       N = INT(X - 0.25)
X       XN=REAL(N)
X   13  IF((X-XN-0.75) .GT. 0.0) THEN 
XC
XC      *                ARGUMENT IS INTEGER                                 *
XC
X          GAMMA=FACTOR(N)
XC
XC      *                ARGUMENT IS HALF-INTEGER                            *
XC
X       ELSE
X   31     PROD=1.77245385D00
X          IF(N .LE. 0) THEN
X             GAMMA=PROD
X             GOTO 56
X          ENDIF
X          IF(N-10 .LE. 0) THEN
X             DO 43 I=1,N
X                FI=REAL(I)
X                PROD=PROD*(FI-0.5)
X   43        CONTINUE
X             GAMMA=PROD
X          ELSE
X   51        SUM=0.0
X             DO 54 I=11,N
X                FI=I
X                SUM=SUM+DLOG(FI-0.5)
X   54        CONTINUE
X   55        GAMMA=PROD*639383.8623D00*DEXP(SUM)
X          ENDIF
X   56  ENDIF
X   60  RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE GEHAN  ******************************
XC      **********************************************************************
XC
X       SUBROUTINE GEHAN(R1,R2,TEST,PROB,XY,ID1,ID2,NTOT)
XC
XC 
XC      * THIS SUBROUTINE COMPUTES GEHAN'S GENERALIZED WILCOXON TEST         *
XC      * STATISTIC.  THE COMPUTATIONAL FORM IS FROM E.T. LEE , STATISTICAL  *
XC      * METHODS FOR SURVIVAL DATA ANALYSIS, 1980, LIFETIME LEARNING        *
XC      * PUBLICATIONS, BELM0NT, CA. THE FORM USED IS THE MANTEL METHOD OF   *
XC      * ASSIGNING A SCORE TO EACH OBSERVATION BASED ON ITS RELATIVE RANK,  *
XC      * FROM EQUATION 5.5 AND EXAMPLE 5.1                                  *
XC      *                                                                    *
XC      *         SUBROUTINES                                                *
XC      *                   STAT                                             *
XC
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION R1(NTOT),R2(NTOT)
X       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
XC
XC      *        COMPUTATION OF R1                                           *
XC      *        STEP 1 AND 2 : RANK FROM LEFT TO RIGHT, OMITTING            *
XC      *        RIGHT CENSORED VALUES. ASSIGN NEXT HIGHER RANK              *
XC      *        TO RIGHT CENSORED VALUES                                    *
XC
XC
X       IRANK=0
X       DO 90 I=1,NCOMP
X          IF(ID1(I).EQ.1) THEN
X             R1(I)=IRANK+1
X          ELSE
X             IRANK=IRANK+1
X             R1(I)=REAL(IRANK)
X          ENDIF
X  90   CONTINUE
XC
XC      *        STEP3 : REDUCE THE RANK OF TIED OBSERVATIONS TO             *
XC      *        THE LOWEST RANK FOR THE VALUE                               *
XC
X       K1=NCOMP-1
X       L1=1
X  12   IF(XY(L1).EQ.XY(L1+1)) THEN
XC
X          JEMP=IABS(ID1(L1)-1)*IABS(ID1(L1+1)-1)
X          IF(JEMP.NE.0) THEN
X             R1(L1+1)=R1(L1)
X             IF(L1.EQ.K1) GOTO 13
X             L1=L1+1
X             GOTO 12
X          ENDIF
X       ENDIF
X       IF(L1.NE.K1) THEN
X          L1=L1+1
X          GOTO 12
X       ENDIF
XC
XC      *             COMPUTATION OF R2                                      *
XC      *             STEP 1 : RANK FROM RIGHT TO LEFT                       *
XC
X  13   DO 14 I=1,NCOMP
X          R2(I)=REAL(NCOMP-I+1)
X  14   CONTINUE
XC
XC      *          STEP2 : REDUCE THE RANK OF TIED OBSERVATIONS              *
XC      *          TO THE LOWEST FOR THE VALUE                               *
XC
X       L1=NCOMP
X  22   IF(XY(L1).EQ.XY(L1-1)) THEN
XC
X          JEMP=IABS(ID1(L1)-1)*IABS(ID1(L1-1)-1)
X          IF(JEMP.NE.0) THEN
X             R2(L1-1)=R2(L1)
X             IF(L1.EQ.2) GOTO 23
X             L1=L1-1
X             GOTO 22
X          ENDIF
X       ENDIF
X       IF(L1.NE.2) THEN
X          L1=L1-1
X          GOTO 22
X       ENDIF
X
X  23   IF(NCEN.NE.0) THEN
XC
XC      *         STEP 3 : REDUCE THE RANK OF RIGHT CENSORED                 *
XC      *         OBSERVATION TO UNITY                                       *
XC
X          DO 24 I=1,NCOMP
X             IF(ID1(I).NE.0)  R2(I)=1.0
X  24      CONTINUE
X       ENDIF
XC
XC      *               COMPUTE FINAL SCORES - R1(I)                         *
X       DO 25 I=1,NCOMP
X          R1(I)=R1(I)-R2(I)
X  25   CONTINUE
X
X       CALL STAT(R1,TEST,XY,ID1,ID2,NTOT)
X       PROB=1.0-AGAUSS(TEST)
X
X       RETURN
X       END
X
X
X
XC      **********************************************************************
XC      ******************** SUBROUTINE GRDPROB  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE GRDPRB(NTOT,MX,MY,SUM,ISKIP,ICENS,DELX,
X     +                   DELY,XORG,YORG,TOL,MAX,MM,M1,M2,M3,M4,M5,
X     +                   M6,M7,M8,X,Y,NP,XB,YB,F,FC,N,N1,N2,N3,
X     +                   N4,N5,N6,N7,N8,IWRK1,IWRK2,WRK1,WRK2,
X     +                   SWRK1,DWRK1,IB,MVAR)
XC
XC
XC      *                                                                    *
XC      *     THIS SUBPROGRAM COMPUTES THE PROBABILITY OF BIN(I,J)           *
XC      *     IN WHICH DETECTED POINTS EXIST. ONLY BINS WHICH HAVE           *
XC      *     DETECTED POINTS CAN HAVE NON-ZERO PROBABILITY.                 *
XC      *                                                                    *
XC      *     INPUT                                                          *
XC      *              X(I)  : INDEPENDENT VARIABLE                          *
XC      *              Y(I)  : DEPENDENT VARIABLE                            *
XC      *             NP(I)  : INDICATOR OF CENSORED STATUS                  *
XC      *              NT    : TOTAL NUMBER OF DATA                          *
XC      *              MX    : BIN NUMBER OF X                               *
XC      *              MY    : BIN NUMBER OF Y                               *
XC      *             ISKIP  : INDICATOR OF BINNING                          *
XC      *             ICENS  : CENSORING STATUS OF DATA SET                  *
XC      *             TOL    : TOLERANCE FOR COMPUTAION F(I,J)               *
XC      *             MAX    : MAXIMUM ITERATION                             *
XC      *        IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED  :             *
XC      *             DELX   : BIN SIZE OF X AXIS                            *
XC      *             DELY   : BIN SIZE OF Y AXIS                            *
XC      *             XORG   : ORIGIN OF X                                   *
XC      *             YORG   : ORIGIN OF Y                                   *
XC      *      WORK                                                          *
XC      *             FC(I,J): COPY OF F(I,J) (TO CHECK THE CONVERGENCE)     *
XC      *             N1(I,J): NUMBER OF Y LOWER LIMITS IN THE BIN (I,J)     *
XC      *             N2(I,J): NUMBER OF X LOWER LIMITS IN THE BIN (I,J)     *
XC      *             N3(I,J): NUMBER OF DOUBLE LOWER LIMITS IN THE          *
XC      *                      BIN(I,J)                                      *
XC      *             N4(I,J): NUMBER OF Y LOWER X UPPER LIMTS IN THE        *
XC      *                      BIN(I,J)                                      *
XC      *             N5(I,J): NUMBER OF Y UPPER LIMITS IN THE BIN (I,J)     *
XC      *             N6(I,J): NUMBER OF X UPPER LIMITS IN THE BIN (I,J)     *
XC      *             N7(I,J): NUMBER OF DOUBLE UPPER LIMITS IN THE          *
XC      *                      BIN(I,J)                                      *
XC      *             N8(I,J): NUMBER OF Y UPPER X LOWER LIMTS IN THE        *
XC      *                      BIN(I,J)                                      *
XC      *              SUM1  : WEIGHT ON BIN(I,J) FROM Y LOWER LIMIT         *
XC      *              SUM2  : WEIGHT ON BIN(I,J) FROM X LOWER LIMITS        *
XC      *              SUM3  : WEIGHT ON BIN(I,J) FROM DOUBLE LOWER          *
XC      *                      LIMITS                                        *
XC      *              SUM4  : WEIGHT ON BIN(I,J) FROM Y LOWER, X UPPER      *
XC      *                      LIMITS                                        *
XC      *              SUM5  : WEIGHT ON BIN(I,J) FROM Y UPPER LIMITS        *
XC      *              SUM6  : WEIGHT ON BIN(I,J) FROM X UPPER LIMITS        *
XC      *              SUM7  : WEIGHT ON BIN(I,J) FROM DOUBLE UPPER          *
XC      *                      LIMITS                                        *
XC      *              SUM8  : WEIGHT ON BIN(I,J) FROM Y UPPER, X LOWER      *
XC      *                      LIMITS                                        *
XC      *              ITE   : NUMBER OF ITERATIONS                          *
XC      *              DEL   : TOLERANCE [SUM (F(I,J)-FC(I,J))**2]           *
XC      *     OUTPUT                                                         *
XC      *              F(I,J): NUMBER OF DATA POINTS IN BIN(I,J)             *
XC      *              XB(I) : POSITION OF  X BIN                            *
XC      *              YB(I) : POSITION OF  Y BIN                            *
XC      *                                                                    *
XC      *     SUBROUTINE : BIN                                               *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8(A-H,O-Z), INTEGER (I-N)
X       DIMENSION FC(IB,IB)
X       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB)
X       DIMENSION F(IB,IB),N(IB,IB),N1(IB,IB),N2(IB,IB),N3(IB,IB)
X       DIMENSION N4(IB,IB),N5(IB,IB),N6(IB,IB),N7(IB,IB),N8(IB,IB)
X
X       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
X       DIMENSION DWRK1(MVAR,NTOT),SWRK1(MVAR)
X
X       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
XC
XC      *              CALL SUBROUTINE BIN                                   *
XC
X       CALL BIN(NTOT,MX,MY,ISKIP,ICENS,DELX,DELY,XORG,YORG,MM,M1,M2,
X     +          M3,M4,M5,M6,M7,M8,IWRK1,IWRK2,WRK1,WRK2,DWRK1,SWRK1,
X     +          X,Y,NP,XB,YB,F,N,N1,N2,N3,N4,N5,N6,N7,N8,IB,MVAR)
X
XC
XC      *             INITIAL SETTING OF FC(I,J)                             *
XC
X       DO 300 I=1,MX
X          DO 200 J=1,MY
X             FC(I,J)=0.0
X  200     CONTINUE
X  300  CONTINUE
XC
XC
XC      *             START ITERATIONS TO FIND F(I,J) LOOP 500               *
XC
X       SNT=NTOT
X       ITE=1
X
X  500  DEL=0.0
X       DO 1600 I=1,MX
X          DO 1500 J=1,MY
X             IF(F(I,J).NE.0.0) THEN
X                SUM1=0
X                SUM2=0
X                SUM3=0
X                SUM4=0
X                SUM5=0
X                SUM6=0
X                SUM7=0
X                SUM8=0
XC
XC      *      COMPUTE THE INFLUENCE OF CENSORED DATA POINTS ON THE DETECTED *
XC      *      POINT AT I,J.                                                 *
XC
XC
XC   
XC      *              Y LOWER LIMITS                                        *
XC
X                IF(NC1.NE.0) THEN
X                   DO 600 L=1,J
X                      IF(N1(I,L).NE.0) THEN
X                         SUMF1=0.0
X                         DO 550 L1=L,MY
X                            SUMF1=SUMF1+F(I,L1)
X  550                    CONTINUE
X
X               SUM1=SUM1+(FLOAT(N1(I,L))/SNT)*(F(I,J)/SUMF1)
X
X                      ENDIF
X  600              CONTINUE
X                ENDIF
XC
XC
XC      *             X LOWER LIMITS                                         *
XC
X                IF(NC2.NE.0) THEN
X                   DO 700 K=1,I
X                      IF(N2(K,J).NE.0) THEN
X                         SUMF2=0.0
X                         DO 650 K1=K,MX
X                            SUMF2=SUMF2+F(K1,J)
X  650                    CONTINUE
X
X                SUM2=SUM2+(FLOAT(N2(K,J))/SNT)*(F(I,J)/SUMF2)
X
X                      ENDIF
X  700              CONTINUE
X                ENDIF
XC
XC
XC      *            DOUBLE LOWER LIMITS                                     *
XC
X                IF(NC3.NE.0) THEN
X                   DO 800 K=1,I
X                      DO 790 L=1,J
X                         IF(N3(K,L).NE.0) THEN
X                            SUMF3=0.0
X                            DO 780 K1=K,MX
X                               DO 770 L1=L,MY
X                                  SUMF3=SUMF3+F(K1,L1)
X  770                          CONTINUE
X  780                       CONTINUE
X
X                SUM3=SUM3+(FLOAT(N3(K,L))/SNT)*(F(I,J)/SUMF3)
X
X                         ENDIF
X  790                 CONTINUE
X  800              CONTINUE
X                ENDIF
XC
XC
XC      *            Y LOWER, X UPPER LIMITS                                 *
XC
X                IF(NC4.NE.0) THEN
X                   DO 900 K=1,MX
X                      KK=MX-K+1
X                      IF(KK.LT.I) GOTO 910
X                      DO 890 L=1,J
X                         IF(N4(KK,L).NE.0) THEN
X                            SUMF4=0.0
X                            DO 880 K1=1,KK
X                               DO 870 L1=L,MY
X                                  SUMF4=SUMF4+F(K1,L1)
X  870                          CONTINUE
X  880                       CONTINUE
X
X                SUM4=SUM4+(FLOAT(N4(K,L))/SNT)*(F(I,J)/SUMF4)
X
X                         ENDIF
X  890                 CONTINUE
X  900              CONTINUE
X                ENDIF
XC
XC
XC      *            Y UPPER LIMITS                                          *
XC
X  910           IF(NC5.NE.0) THEN
X                   DO 1000 L=1,MY
X                      LL=MY-L+1
X                      IF(LL.LT.J) GOTO 1010
X                      IF(N5(I,LL).NE.0) THEN
X                         SUMF5=0.0
X                         DO 950 L1=1,LL
X                            SUMF5=SUMF5+F(I,L1)
X  950                    CONTINUE
X
X                SUM5=SUM5+(FLOAT(N5(I,LL))/SNT)*(F(I,J)/SUMF5)
X
X                      ENDIF
X 1000              CONTINUE
X                ENDIF
XC
XC
XC      *               X UPPER LIMITS                                       *
XC
X 1010           IF(NC6.NE.0) THEN
X                   DO 1100 K=1,MX
X                      KK=MX-K+1
X                      IF(KK.LT.I) GOTO 1110
X                      IF(N6(KK,J).NE.0) THEN
X                         SUMF6=0.0
X                         DO 1050 K1=1,KK
X                            SUMF6=SUMF6+F(K1,J)
X 1050                    CONTINUE
X
X                SUM6=SUM6+(FLOAT(N6(KK,J))/SNT)*(F(I,J)/SUMF6)
X
X                      ENDIF
X 1100              CONTINUE
X                ENDIF
XC
XC
XC      *            DOUBLE UPPER LIMITS                                     *
XC
X 1110           IF(NC7.NE.0) THEN
X                   DO 1200 K=1,MX
X                      KK=MX-K+1
X                      IF(KK.LT.I) GOTO 1210
X                      DO 1190 L=1,MY
X                         LL=MY-L+1
X                         IF(LL.LT.J) GOTO 1200
X                         IF(N7(KK,LL).NE.0) THEN
X                            SUMF7=0.0
X                            DO 1180 K1=1,KK
X                               DO 1170 L1=1,LL
X                                  SUMF7=SUMF7+F(K1,L1)
X 1170                          CONTINUE
X 1180                       CONTINUE
X
X                SUM7=SUM7+(FLOAT(N7(KK,LL))/SNT)*(F(I,J)/SUMF7)
X
X                         ENDIF
X 1190                 CONTINUE
X 1200              CONTINUE
X                ENDIF
XC
XC
XC      *               Y UPPER, X LOWER LIMITS                              *
XC
X 1210           IF(NC8.NE.0) THEN
X                   DO 1300 K=1,I
X                      DO 1290 L=1,MY
X                         LL=MY-L+1
X                         IF(LL.LT.J) GOTO 1300
X                         IF(N8(K,LL).NE.0) THEN
X                            SUMF8=0.0
X                            DO 1280 K1=K,MX
X                               DO 1270 L1=1,LL
X                                  SUMF8=SUMF8+F(K1,L1)
X 1270                          CONTINUE
X 1280                       CONTINUE
X
X                SUM8=SUM8+(FLOAT(N8(KK,LL))/SNT)*(F(I,J)/SUMF8)
X
X                         ENDIF
X 1290                 CONTINUE
X 1300              CONTINUE
X                ENDIF
XC
XC      *            COMPUTE A NEW ESTIMATE OF F(I,J).                       *
XC
X 1400           SUM=SUM1+SUM2+SUM3+SUM4+SUM5+SUM6+SUM7+SUM8
X                F(I,J)=FLOAT(N(I,J))/SNT+SUM
X                DEL=DEL+(F(I,J)-FC(I,J))**2
X                FC(I,J)=F(I,J)
X             ENDIF
X 1500     CONTINUE
X 1600  CONTINUE
XC
XC      *             CHECK CONVERGENCE                                      *
XC
X       ITE=ITE+1
X
X       IF(((DEL).GT.TOL).AND.(ITE.LE.MAX)) GOTO 500
X       RETURN
X       END
X
X
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE KMADJ  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE KMADJ(ZU,ZC,NTOT,IU,IC,S,ISIGN,NTEMP,F,V)
XC
XC      *       THIS SUBROUTINE RESTORES THE DATA AND THE PRODUCT-LIMIT      *
XC      *       ESTIMATOR TO UPPER-LIMITS FORM, IF THE DATA WAS IN THAT FORM *
XC      *       INITIALLY.  TIES AT CENSORED POINTS ARE ALSO REMOVED TO      *
XC      *       MAKE THE PRINTOUT CONSISTENT.                                *
XC      *                                                                    *
XC      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
XC      *                ZC(I)  :  CENSORED DATA POINTS                      *
XC      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
XC      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
XC      *                 IC    :  NUMBER OF CENSORED DATA POINTS            *
XC      *                 S(L)  :  PL ESTIMATOR                              *
XC      *       OUTPUT  NTEMP   :  VALUE OF NTOT AFTER ADJUSTMENT FOR TIES   *
XC      *       OTHER    F      :  PROBABILITY MASS ASSIGNED TO EACH         *
XC      *                             UNCENSORED POINT(=JUMP IN S AT THE     *
XC      *                                                  POINT)            *
XC      *                                                                    *
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION ZU(NTOT),ZC(NTOT),S(NTOT),F(NTOT),V(NTOT)
X
XC
XC      *  FOR LEFT-CENSORED DATASETS (I.E. UPPER LIMITS),                  *
XC      *  CALCULATE JUMP IN SURVIVOR FUNCTION AT UNCENSORED POINTS         *
XC
X       IF(ISIGN.LT.0) THEN
X            F(1) = 1.0 - S(1)
X            DO 120 I = 2, IU
X               F(I) = S(I-1)-S(I)
X120         CONTINUE
XC
XC      *  REVERSE SEQUENCE OF UNCENSORED POINTS, JUMPS AND ERRORS         *
XC
X            J = IU/2
X            DO 150 I =1, J
X
X               Z = ZU(I)*(-1.0)
X               ZU(I) = ZU(IU-I+1)*(-1.0)
X               ZU(IU-I+1) = Z
X
X               FTEMP = F(I)
X               F(I) = F(IU-I+1)
X               F(IU-I+1) = FTEMP 
X
X               VTEMP = V(I)
XC               V(I) = V(IU-I+1)
XC               V(IU-I+1) = VTEMP
X               V(I) = V(IU-I)
X               V(IU-I) = VTEMP
X
X150         CONTINUE
X
X            IF(2*J.NE.IU) THEN
X               ZU(J+1) = ZU(J+1)*(-1.0)
X            ENDIF
X
XC
XC      *  REVERSE SEQUENCE OF CENSORED POINTS                              *
XC
X            J = IC/2
X            DO 155 I = 1, J
X               Z = ZC(I) * (-1.0)
X               ZC(I) = ZC(IC-I+1)*(-1.0)
X               ZC(IC-I+1) = Z
X155         CONTINUE
X 
X            IF(2*J.NE.IC) THEN
X               ZC(J+1) = ZC(J+1)*(-1.0)
X            ENDIF
X
XC
XC      *  COMPUTE SURVIVOR FUNCTION FOR REVERSED DATA                     *
XC
X            DO 170 I = 1, IU
X               S(I) = 1
X               DO 160 J = 1, I
X                  S(I) = S(I) - F(J)
X160            CONTINUE
X170         CONTINUE
X         ENDIF   
X
XC      *   CORRECTS FOR TIES AT THE UNCENSORED POINTS                      *
XC      *   NOTICE THAT IF THERE ARE TIES AT THESE POINTS, THEN BOTH        *
XC      *   IU AND NTEMP ARE REDUCED.                                       *
X
X       NTEMP = NTOT
X       K = 1
X190    IF(ZU(K).EQ.ZU(K+1)) THEN
X          DO 195 I = K, IU-1
X             ZU(I)=ZU(I+1)
X             S(I)=S(I+1)
X             V(I) = V(I+1)               
X195       CONTINUE
X          IU = IU -1
X          NTEMP = NTEMP - 1
X       ELSE
X          K  = K + 1
X       ENDIF
X       IF(K.LT.IU) GOTO 190
X
X       RETURN
X       END
X
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE KMDIF ********************************
XC      **********************************************************************
XC
X       SUBROUTINE KMDIF(S,ZU,BS,BL,DIFF,F,NTOT,START,BINSIZ,LSTEP,
X     +                   OUT,IBIN,IU)
X
XC
XC      *       THIS SUBROUTINE COMPUTES AND PRINTS THE DIFFERENTIAL KM      *
XC      *       ESTIMATOR BASED ON WARDLE AND KNAPP, 'THE STATISTICAL        *
XC      *       DISTRIBUTION OF THE NEUTRAL-HYDROGEN CONTENT OF S0 GALAXIES',*
XC      *       ASTRN. J., 91:23 1986.                                       *
XC      *                                                                    *
XC      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
XC      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
XC      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
XC      *                 S(L)  :  PL ESTIMATOR                              *
XC      *                 OUT   :  OUTPUT FILE NAME. IF IT IS BLANK          *
XC      *                          THE RESULTS WILL BE SHOWN ON THE          *
XC      *                          TERMINAL.                                 *
XC      *                START  :  STARTING VALUE OF THE FIRST BIN           *
XC      *                BINSIZ :  WIDTH OF THE BIN                          *
XC      *                LSTEP  :  NUMBER OF BINS                            *
XC      *                IBIN   :  DIMENSION                                 *
XC      *              ICHANGE  :  INDICATES IF THE LAST POINT (OR THE       *
XC      *                            FIRST POINT FOR UPPER LIMITS DATA)      *
XC      *                            HAS BEEN CHANGED TO A DETECTION.        *
XC      *                                                                    *
XC      *      OTHERS                                                        *
XC      *               BS(J)   :  STARTING VALUE FOR THE BIN J              *
XC      *               BL(J)   :  ENDING VALUE FOR THE BIN J                *
XC      *               DIFF(J) :  DIFFERENTIAL KM ESTIMATOR AT BIN J        *
XC      *               F(I)    :  MASS OF THE I TH DATA POINT               *
XC      *                                                                    *
XC      *                                                                    *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*9 OUT,CHECK
X
X       DIMENSION ZU(NTOT),S(NTOT),F(NTOT)
X       DIMENSION BS(IBIN),BL(IBIN),DIFF(IBIN)
X       CHECK='         '
XC
XC       *    FIRST, COMPUTE THE MASS OF EACH POINT                          *
XC
X        F(1) = 1.0 -S(1)
X       
X        DO 610 I = 2, IU
X           F(I) = DABS(S(I) - S(I-1))
X  610   CONTINUE
X
XC
XC       *  SET THE BIN BOUNDARIES.                                          *
XC
X        DO 620 J = 1, LSTEP
X           BS(J) = START + BINSIZ*(J-1)
X           BL(J) = START + BINSIZ*J
X  620   CONTINUE
X
X        I = 1
X        J = 0
X
X  630   J = J + 1
X        DIFF(J) = 0.0
X
XC
XC      *       CHECK WHETHER THE I-TH POINT IS IN THE BIN                  *
XC
X  640  IF(J .LE. LSTEP) THEN
X         IF(ZU(I) .LT. BS(J)) THEN
X            IF(I .GE. IU) THEN
X               GOTO 630
X            ENDIF
X            I = I + 1
X            GOTO 640
X         ENDIF
X
XC      *       COMPUTE THE DIFFERENTIAL KM                                 *
XC
X         IF(ZU(I) .GE. BL(J)) GOTO 630
X         DIFF(J) = DIFF(J) + F(I)
X   
X         IF(I .LT. IU) THEN
X            I = I + 1
X            GOTO 640
X         ENDIF
X         GOTO 630
X       ENDIF
XC
XC      *         START PRINTING THE RESULTS                                  *
XC
X          IF(OUT.EQ.CHECK) THEN
X             PRINT *
X             PRINT *,'           DIFFERENTIAL KM ESTIMATOR'
X             PRINT 660
X             PRINT *
X          ELSE
X             WRITE(60, 658)
X             WRITE(60,659)
X             WRITE(60, 660)
X             WRITE(60,658)
X          ENDIF
X  658     FORMAT('  ')
X  659     FORMAT(5X,'   DIFFERENTIAL KM ESTIMATOR')
X  660     FORMAT(5X,'   BIN CENTER          D')
X
XC
XC      * MULTIPLY DIFF(I) BY THE TOTAL NUMBER OF POINTS TO GET THE NUMBER    *
XC      * OF POINTS IN EACH BIN.                                              *
XC
X          SUM = 0.0
X          DO 690 I = 1, LSTEP
X             DIFF(I) =DIFF(I)*NTOT
X             CENTER = 0.5*(BS(I) + BL(I))
X             IF(OUT .EQ. CHECK) THEN
X                PRINT 680, CENTER, DIFF(I)
X             ELSE
X                WRITE(60,680) CENTER, DIFF(I)
X             ENDIF
X  680        FORMAT(2F15.3)
X             SUM = SUM + DIFF(I)
X  690     CONTINUE
X
X          IF(OUT .EQ. CHECK) THEN
X             PRINT 700, SUM
X             PRINT 658
X          ELSE
X             WRITE(60,700) SUM
X             WRITE(60,658)
X             WRITE(60,701)
X 701         FORMAT(' (D GIVES THE ESTIMATED DATA POINTS IN EACH BIN)')
X          ENDIF
X
X 700      FORMAT(23X,'-------',/10X,'SUM =',F15.3)
X
X       RETURN
X       END
X
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE KMESTM ******************************
XC      **********************************************************************
XC
X       SUBROUTINE KMESTM(IND,X,NTOT,J,IPRINT,TITLE,NAME,OUTPUT,IBIN,
X     +                   ISKIP,KDIFF,START,BINSIZ,LSTEP,FILE,
X     +                   ZXU,ZXC,SX,VX,Z1,ITEMP,INDEX,ATIE,RISK,
X     +                   BWRK1,BWRK2,BWRK3,IWRK1,SWRK1,
X     +                   WRK1,MVAR)
XC
XC      *       THIS SUBROUTINE COMPUTES THE PL ESTIMATOR, MEAN AND ITS      *
XC      *       ERROR FOR THE X VARIABLE.                                    *
XC      *                                                                    *
XC      *       INPUT IND(J,I): INDICATOR OF CENSORING                       *
XC      *               X(J,I): DATA                                         *
XC      *                NTOT : NO. OF DATA POINTS                           *
XC      *                 J   : J-TH VARIABLES                               *
XC      *               IPRINT: IF 0, PRINT OUT RESULTS ONLY                 *
XC      *                       IF 1, PRINT OUT ALL                          *
XC      *              PROBLEM: TITLE OF THE PROBLEM                         *
XC      *               NAME  : NAME OF THE SUB-DATA SET                     *
XC      *               OUTPUT: NAME OF OUTPUT FILE                          *
XC      *                       IF IT IS BLANK, SHOW THE RESULT ON THE       *
XC      *                       TERMINAL.                                    *
XC      *               ISKIP : IF THE SUBROUTINE IS CALLED BY TWO SAMPLE    *
XC      *                       TESTS, ISKIP=1 AND SKIP A FEW LINES .        *
XC      *               KDIFF : IF KDIFF = 1, PRINT OUT DIFFERENTIAL KM      *
XC      *               START : STARTING POINT OF BINING                     *
XC      *               BINSIZ: WIDTH OF THE BIN                             *
XC      *               LSTEP : NUMBER OF BINS                               *
XC      *              ATIE(I): NUMBER OF TIED DATA POINTS AT ITH DATA VALUE *
XC      *              RISK(I): RISK SET FOR THE ITH DATA VALUE              *
XC      *               MVAR  : DIMENSION SIZE                               *
XC      *               IBIN  : DIMENSION SIZE                               *
XC      *                                                                    *
XC      *       WORK      ZXU : DATA ARRAY CONTAINING THE UNCENSORED POINTS  *
XC      *                 ZXC : DATA ARRAY CONTAINING THE CENSORED POINTS    *
XC      *                 IXU : NO. OF UNCENSORED DATA POINTS                *
XC      *                 IXC : NO. OF CENSORED DATA POINTS                  *
XC      *              ICHANGE: IF THE LAST VALUE IS CENSORED, THE VALUE     *
XC      *                       NEEDS TO BE CHANGED TO A DETECTION.          *
XC      *                       THIS INDICATOR IS SET TO BE -1,IF THE LAST   *
XC      *                       VALUE IS CHANGED.                            *
XC      *                                                                    *
XC      *       OUTPUT    SX  : PL ESTIMATOR                                 *
XC      *                 VX  : ERROR OF PL ESTIMATOR                        *
XC      *                SMEAN: MEAN                                         *
XC      *                ERROR: ERROR OF THE MEAN                            *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *                 SORT1, XVAR, PLESTM, KMADJ, KMPRNT                 *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*80 TITLE
X       CHARACTER*9 NAME,OUTPUT,FILE
X
X       DIMENSION IND(MVAR, NTOT),X(MVAR, NTOT),ZXU(NTOT),ZXC(NTOT)
X       DIMENSION SX(NTOT),VX(NTOT),Z1(MVAR,NTOT),ITEMP(NTOT)
X       DIMENSION INDEX(NTOT),ATIE(NTOT),RISK(NTOT)
X       DIMENSION BWRK1(IBIN),BWRK2(IBIN),BWRK3(IBIN)
X       DIMENSION IWRK1(NTOT),SWRK1(MVAR),WRK1(NTOT)
X
XC
XC      *       DISTINGUISH UNCENSORED AND CENSORED DATA AND SORT THEM IN    *
XC      *       ASCENDING ORDER. THEN CALL PL ESTIMATOR SUBROUTINE           *
XC      *       "PL".                                                        *
XC
X       IF(OUTPUT .EQ. '         ') THEN
X
X          PRINT *
X          PRINT *
X          PRINT 40
X          PRINT 44
X          IF(ISKIP.NE.1) PRINT 46,TITLE   
X          PRINT 44
X          PRINT 47,FILE
X          PRINT 44
X          PRINT 48,NAME
X          PRINT 44
X       ELSE
X          WRITE(60,44)
X          WRITE(60,44)
X          WRITE(60,40)
X          WRITE(60,44)
X          IF(ISKIP.NE.1) WRITE(60,46) TITLE   
X          WRITE(60,44)
X          WRITE(60,47) FILE
X          WRITE(60,44)
X          WRITE(60,48) NAME
X          WRITE(60,44)
X       ENDIF
X
X   40  FORMAT(8X,'KAPLAN-MEIER ESTIMATOR')
X   44  FORMAT('    ')
X   46  FORMAT(8X,'TITLE : ',A80)
X   47  FORMAT(8X,'DATA FILE : ',A9)  
X   48  FORMAT(8X,'VARIABLE : ',A9)
X
XC
XC      *    XVAR DISTINGUISHES DETECTED POINTS AND CENSORED POINTS           *
XC
X       CALL XVAR(IND,X,J,NTOT,ISIGN,ZXU,ZXC,IXU,IXC,IWRK1,
X     +           ATIE,RISK,WRK1,Z1,SWRK1,LTOT,MVAR,INDEX)
XC
X       IF(OUTPUT .EQ. '         ') THEN
X
X          PRINT *
X          IF(ISIGN.EQ.1) PRINT 56,NTOT,IXC
X          IF(ISIGN.EQ.-1) PRINT 57,NTOT,IXC
X          PRINT *
X       ELSE
X          WRITE(60,44)
X          IF(ISIGN.EQ.1) WRITE(60,56) NTOT,IXC
X          IF(ISIGN.EQ.-1) WRITE(60,57) NTOT,IXC
X          WRITE(60,44)
X       ENDIF
X
X   56  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF LOWER LIMITS :',I4)
X   57  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF UPPER LIMITS :',I4)
XC
XC      *  SET A FEW DUMMY ARRAYS TO USE THE SORTING PRGRAM                  *
XC
X       ISTEMP = ISIGN
X   59  DO 60 I=1,NTOT
X          ITEMP(I)=0
X          Z1(1,I)=1.0
X   60  CONTINUE
XC
X       CALL SORT1(ITEMP,Z1,ZXU,IXU,1,INDEX,SWRK1,MVAR)
XC
X       CALL SORT1(ITEMP,Z1,ZXC,IXC,1,INDEX,SWRK1,MVAR)
XC
XC      *  CALL SUBROUTINE "PLESTM" TO COMPUTE KM ESTIMATOR                  *
XC
X       CALL PLESTM(ZXU,ZXC,IXU,IXC,SX,VX,NTOT,SMEAN,ERROR,ICHANGE,
X     +             NCHANGE,IWRK1)
XC
XC      *        IF THE DATA CONTAINS UPPER LIMITS, CHANGE THE               *
XC      *        SIGN OF THE MEAN.                                           *
XC
X       ISIGN = ISTEMP
X       SMEAN=ISIGN*SMEAN
XC
XC      * SUBROUTINE KMADJ IS CALLED TO ADJUST THE PRODUCT-LIMIT ESTIMATOR   *
XC      * BACK TO THE ORIGIONAL CENSORING PATTERN AND TO REMOVE TIES.        *
XC
X       CALL KMADJ(ZXU,ZXC,NTOT,IXU,IXC,SX,ISIGN,NTEMP,WRK1,VX)
X
XC
XC      * PRINT PL ESTIMATOR, PERCENTILES, MEAN, AND ERROR                   *
XC
XC
X
X       CALL KMPRNT(ZXU,ZXC,NTOT,NTEMP,IXU,IXC,SX,VX,ISIGN,OUTPUT,
X     +             ICHANGE,SMEAN,ERROR,IPRINT)
X
X
XC      *  SUBROUTINE KMDIF IS CALLED IF THE USER HAS REQUESTED A            *
XC      *  DIFFERENTIAL KM ESTIMATOR.                                        *
X
X       IF(KDIFF .EQ. 1) THEN
X
X          CALL KMDIF(SX,ZXU,BWRK1,BWRK2,BWRK3,WRK1,NTOT,START,
X     +               BINSIZ,LSTEP,OUTPUT,IBIN,IXU)
X
X       ENDIF
X
X       PRINT *   
X
XC
XC      *   IF THE LAST VALUE WAS CHANGED FROM AN UPPER LIMIT TO A          *
XC      *   DETECTION, CHANGE THE NUMBER BACK TO ITS ORIGINAL VALUE.        *
XC
X       IF(ICHANGE.EQ.-1) THEN
X          IXU=IXU-NCHANGE
X          IXC=IXC+NCHANGE
X       ENDIF
X       
X       RETURN
X       END
X
X
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE KMPRNT *******************************
XC      **********************************************************************
XC
X       SUBROUTINE KMPRNT(ZU,ZC,NTOT,NTEMP,IU,IC,S,V,ISIGN,OUT,ICHANGE,
X     +                   SMEAN,ERROR,IPRINT)
X
XC
XC      *       THIS SUBROUTINE PRINTS KM ESTIMATORS, THEIR                  *
XC      *       ERROR, AND 75, 50, AND 25 PERCENTILES. ADOPTED FROM          *
XC      *       ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL DATA         *
XC      *       ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA) *
XC      *                                                                    *
XC      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
XC      *                ZC(I)  :  CENSORED DATA POINTS                      *
XC      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
XC      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
XC      *                 IC    :  NUMBER OF CENSORED DATA POINTS            *
XC      *                 S(L)  :  KM ESTIMATOR                              *
XC      *                 V     :  ERROR OF KM ESTIMATOR                     *
XC      *                ISIGN  :  INDICATOR OF LOWER/UPPER LIMIT            *
XC      *                            IF 1, LOWER LIMIT                       *
XC      *                           IF -1, UPPER LIMIT                       *
XC      *                 D(L)  :  NUMBER OF TIED DATA POINTS AT THE VALUE   *
XC      *                 R(L)  :  RISK SET                                  *
XC      *                 OUT   :  OUTPUT FILE NAME. IF IT IS BLANK          *
XC      *                          THE RESULTS WILL BE SHOWN ON THE          *
XC      *                          TERMINAL.                                 *
XC      *                KDIFF  :  PRINTOUT INDICATOR FOR THE DIFFERENTIAL   *
XC      *                          KM ESTIMATOR (IF 1, PRINT )               *
XC      *                START  :  STARTING VALUE OF THE FIRST BIN           *
XC      *                BINSIZ :  WIDTH OF THE BIN                          *
XC      *                LSTEP  :  NUMBER OF BINS                            *
XC      *                IBIN   :  DIMENSION                                 *
XC      *              ICHANGE  :  INDICATES IF THE LAST POINT (OR THE       *
XC      *                            FIRST POINT FOR UPPER LIMITS DATA)      *
XC      *                            HAS BEEN CHANGED TO A DETECTION.        *
XC      *               IPRINT  :  INDICATES WHETHER TO PRINT OUT THE        *
XC      *                            FULL KM ESTIMATE OR ONLY THE MEAN       *
XC      *                            AND PERCENTILES                         * 
XC      *                                                                    *
XC      *      OTHERS                                                        *
XC      *               BS(J)   :  STARTING VALUE FOR THE BIN J              *
XC      *               BL(J)   :  ENDING VALUE FOR THE BIN J                *
XC      *               DIFF(J) :  DIFFERENTIAL KM ESTIMATOR AT BIN J        *
XC      *               ERR(J)  :  ERROR FOR THE BIN J                       *
XC      *               F(I)    :  MASS OF THE I TH DATA POINT               *
XC      *               SMEAN   :  MEAN                                      *
XC      *               ERROR   :  ERROR OF THE MEAN                         *
XC      *                                                                    *
XC      *      SUBROUTINE:  SUMRY                                            *
XC      *                                                                    *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*9 OUT,CHECK
X
X       DIMENSION ZU(NTOT),ZC(NTOT),S(NTOT),V(NTOT),FINT(3)
X       CHECK='         '
XC
X       IF(IPRINT .EQ. 1) THEN
X          IF(OUT.EQ.CHECK) THEN
X             PRINT 100
X             PRINT 110
X          ELSE 
X             WRITE(60,100)
X             WRITE(60,110)
X          ENDIF
X
X  100     FORMAT('    ')
X  110     FORMAT(10X,'VARIABLE RANGE',6X,'KM ESTIMATOR',3X,'ERROR')
X
XC      *  STARTS TO PRINT OUT THE RESULTS                                   *
X
X          IF(OUT .EQ. CHECK) THEN
X             PRINT 540,0.00,ZU(1),1.00
X          ELSE
X             WRITE(60,540) 0.00,ZU(1),1.00
X          ENDIF
X
X  540     FORMAT('FROM',F9.3,'   TO',F9.3,F12.3)
X               
X          KT=1
X          KC=1
X          KU=1
X
X  200     IF(KT.LE.NTEMP) THEN
X             IF(KC.GT.IC) GOTO 300
X             IF(KU.GT.IU) GOTO 250
X             IF(ZU(KU).LT.ZC(KC)) GOTO 300
XC
XC
X  250        IF(OUT.EQ.CHECK) THEN
X                PRINT 555,ZC(KC)
X             ELSE
X                WRITE(60,555) ZC(KC)
X             ENDIF
X             KC=KC+1
X             GOTO 400
X
X  300        IF(KU .LT. IU) THEN
X                IF(OUT.EQ.CHECK) THEN
X                   PRINT 550,ZU(KU),ZU(KU+1),S(KU),V(KU)
X                ELSE
X                   WRITE(60,550) ZU(KU),ZU(KU+1),S(KU),V(KU)
X                ENDIF
X             ELSE
X                IF(OUT.EQ.CHECK) THEN
X                   PRINT 560,ZU(KU), S(KU), V(KU)
X                ELSE
X                   WRITE(60,560) ZU(KU),S(KU),V(KU)
X                ENDIF
X             ENDIF
X             KU=KU+1
X  400        KT=KT+1
X             GOTO 200
X
X  550        FORMAT('FROM',F9.3,'   TO',F9.3,2F12.3)
X  555        FORMAT(F13.3,' C ')
X  560        FORMAT('FROM',F9.3,'   ONWARDS',4x,2F12.3)
X          ENDIF
X       ENDIF
X      
XC      *   PRINTS OUT A WARNING FLAG IF A CENSORED POINT HAS BEEN 
XC      *   CHANGED TO A DETECTION
X       IF(ICHANGE .EQ. -1) THEN
X           IF(ISIGN .GT. 0) THEN 
X               IF(OUT .EQ. CHECK) THEN
X                   PRINT 565
X                   PRINT 566
X               ELSE
X                   WRITE(60,565)
X                   WRITE(60,566)
X               ENDIF
X           ELSE
X               IF(OUT .EQ. CHECK) THEN
X                    PRINT 570
X                    PRINT 566
X               ELSE
X                    WRITE(60,570)
X                    WRITE(60,566)
X               ENDIF
X           ENDIF
X       ENDIF
X
X  565     FORMAT(10X,
X     +    'WARNING:  THE LAST POINT WAS CHANGED TO A DETECTION ')
X  566     FORMAT(20X,'FOR THE K-M COMPUTATION')
X  570     FORMAT(10X,
X     +    'WARNING:  THE FIRST POINT WAS CHANGED TO A DETECTION')
X  
XC
XC      *       COMPUTE 75, 50, AND 25 PERCENTILES AND PRINT THEM.           *
XC
X       IF(IU .LE. 3) THEN
X          IF(OUT .EQ. CHECK) THEN
X             PRINT 100
X             PRINT 750
X             PRINT 755
X             PRINT 100
X          ELSE 
X             WRITE(60,100)
X             WRITE(60,750)
X             WRITE(60,755)
X             WRITE(60,100)
X          ENDIF
X          GOTO 900
X       ELSE
X          CALL SUMRY(ZU,IU,S,NTOT,FINT)
X       ENDIF
X
X 750   FORMAT(/,6X,'SINCE THERE ARE LESS THAN 4 UNCENSORED POINTS,')
X 755   FORMAT(6X,'NO PERCENTILES WERE COMPUTED.')
X
X       IF(OUT.NE.CHECK) GOTO 760
X       PRINT 100
X       PRINT 780
X       PRINT 800
X       PRINT 850,(FINT(J),J=1,3)
X       PRINT 100
X       GOTO 900
X
X  760  WRITE(60,100)
X       WRITE(60,780)
X       WRITE(60,800)
X       WRITE(60,850) (FINT(J),J=1,3)
X       WRITE(60,100)
X
X  780  FORMAT(5X,'   PERCENTILES    ')
X  800  FORMAT(5X,'    75 TH     50 TH     25 TH')
X  850  FORMAT(5X,3F10.3)
X
X 900   IF (ISIGN.EQ.1) THEN
X          ZXX=ZU(IU)
X       ELSE
X          ZXX=ZU(1)
X       ENDIF
X
X       IF(OUT .EQ. CHECK) THEN
X          IF(ICHANGE.EQ.-1) THEN
X             WRITE(6,1000) SMEAN,ERROR,ZXX      
X             WRITE(6,1005)
X             WRITE(6,1006)
X          ELSE
X             WRITE(6,1010) SMEAN,ERROR
X          ENDIF
X          WRITE(6,1020)
X       ELSE
X          IF(ICHANGE.EQ.-1) THEN
X             WRITE(60,1000) SMEAN,ERROR,ZXX      
X             WRITE(60,1005)
X             WRITE(60,1006)
X          ELSE
X             WRITE(60,1010) SMEAN,ERROR
X          ENDIF
X          WRITE(60,1020)
X       ENDIF
X
X       PRINT *
X       PRINT *
X 1000  FORMAT(8X,'MEAN=',F10.3,' +/-',F6.3,'   LIMITED TO ',F10.3)
X 1005  FORMAT(/10X,'SINCE A CENSORED POINT WAS CHANGED TO A DETECTION,')
X 1006  FORMAT(10X,'THE MEAN ESTIMATE IS BIASED.')
X 1010  FORMAT(8X,'MEAN=',F10.3,' +/-',F6.3)
X 1020  FORMAT(' ')
X
X       RETURN
X       END
X
X
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE MATINV  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE MATINV(ARRAY,NVAR,DET,IK,JK,MVAR)
XC
XC      *                                                                    *
XC      *      THIS PROGRAM IS ADOPTED FROM PHILIP R. BEVINGTON, "DATA       *
XC      *      REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES', 1969 *
XC      *      McGRAW HILL (NY:NY), PROGRAM B-2 P. 302. SEVERAL MINOR        *
XC      *      MODIFICATIONS WERE DONE BY T. ISOBE.                          *
XC      *                                                                    *
XC      *      INPUT   :   ARRAY(I,J) : SYMMETRIC MATRIX                     *
XC      *                    NVAR     : DIMENSION                            *
XC      *      WORK    :     AMAX     : LARGEST NO. ON WORKING               *
XC      *      OUTPUT  :   ARRAY(I,J) : SYMMETRIC INVERSE MATRIX OF THE      *
XC      *                               ORIGINAL ARRAY(I,J)                  *
XC  
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION ARRAY(MVAR,MVAR),IK(MVAR),JK(MVAR)
XC
X       DET=1.0
XC
X       DO 100 K=1,NVAR
XC
XC      *      FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX             *
XC
X          AMAX=0.0       
XC
X   21     DO 30 I=K,NVAR
X             DO 29 J=K,NVAR
X                IF(DABS(AMAX)-DABS(ARRAY(I,J)) .GT. 0.0) GOTO 30
X   24           AMAX=ARRAY(I,J)
X                IK(K)=I
X                JK(K)=J
X   29        CONTINUE
X   30     CONTINUE
XC
XC      *      INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)        *
XC
X          IF(AMAX.EQ.0.0) THEN
X             DET=0.0
X             RETURN
X          ENDIF
X
X   41     I=IK(K)
X          IF((I-K) .LT. 0) GOTO 21
X
X          IF((I-K).GT.0) THEN
X   43        DO 50 J=1,NVAR
X                SAVE=ARRAY(K,J)
X                ARRAY(K,J)=ARRAY(I,J)
X                ARRAY(I,J)=-SAVE
X   50        CONTINUE
X          ENDIF
X
X   51     J=JK(K)
X          IF((J-K).LT.0) GOTO 21
X
X          IF((J-K).GT.0) THEN
X   53        DO 60 I=1,NVAR
X                SAVE=ARRAY(I,K)
X                ARRAY(I,K)=ARRAY(I,J)
X                ARRAY(I,J)=-SAVE
X   60        CONTINUE
X          ENDIF
XC
XC      *         ACCUMULATE ELEMENTS OF INVERSE MATRIX                      *
XC
X   61     DO 70 I=1,NVAR
X             IF(I.NE.K) ARRAY(I,K)=-ARRAY(I,K)/AMAX
X   70     CONTINUE
X
X   71     DO 80 I=1,NVAR
X             DO 79 J=1,NVAR
X                IF(I.EQ.K) GOTO 79
X                IF(J.EQ.K) GOTO 79
X                ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
X   79        CONTINUE
X   80     CONTINUE
X
X   81     DO 90 J=1,NVAR
X             IF(J.NE.K) ARRAY(K,J)=ARRAY(K,J)/AMAX
X   90     CONTINUE
X
X          ARRAY(K,K)=1.0/AMAX
X          DET=DET*AMAX
X  100  CONTINUE  
XC
XC      *            RESTORE ORDERING OF MATRIX                              *
XC
X  101  DO 130 L=1,NVAR
X          K=NVAR-L+1
X          J=IK(K)
X
X          IF(J.GT.K) THEN
X             DO 110 I=1,NVAR
X                SAVE=ARRAY(I,K)
X                ARRAY(I,K)=-ARRAY(I,J)
X                ARRAY(I,J)=SAVE
X  110        CONTINUE
X          ENDIF
X
X          I=JK(K)
X          IF(I.GT.K) THEN
X             DO 120 J=1,NVAR
X                SAVE=ARRAY(K,J)
X                ARRAY(K,J)=-ARRAY(I,J)
X                ARRAY(I,J)=SAVE
X  120        CONTINUE
X          ENDIF
X
X  130  CONTINUE
X
X       RETURN
X       END
X
XC
XC***************************************************************************
XC**************************  SUBROUTINE LRANK  *****************************
XC***************************************************************************
XC
XC
X       SUBROUTINE LRANK(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R, D1, E1, R1, 
X     +                  D2, E2, R2, SCORE, VAR)
XC     *
XC     * THIS SUBROUTINE COMPUTES THE LOGRANK STATISTIC WITH CONDITIONAL     *
XC     * PERMUTATION VARIANCE (HYPERGEOMETRIC VARIANCE) FROM EQUATIONS (2.2) *
XC     * AND (2.3) IN LATTA, 'A MONTE-CARLO STUDY OF SOME TWO-SAMPLE RANK    *
XC     * TESTS WITH CENSORED DATA', 1981, JOURNAL OF THE AMERICAN STATISTICAL*
XC     * ASSOCIATION, VOL 76, PP 713-719.                                    *
XC     *                                                                     *
XC     * INPUT                                                               *
XC     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
XC     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
XC     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
XC     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
XC     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
XC     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
XC     *                                                                     *
XC     * OUTPUT                                                              *
XC     *     TEST    : STANDARDIZED LOGRANK STATISTIC                        *
XC     *     PROB    : PROBABILITY                                           *
XC     *                                                                     *
XC     * OTHERS                                                              *
XC     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
XC     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
XC     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
XC     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
XC     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
XC     *      R(I)   : RISK SET AT XY(I)                                     *
XC     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
XC     *      SCORE  : SCORE OF THE DATA                                     *
XC     *      VAR    : VARIANCE OF THE DATA                                  *
X
X
X
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X
X       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
X       DIMENSION D(NTOT),E(NTOT),R(NTOT)
X       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
X       DIMENSION E2(NTOT),R2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X
X       I = 1
X       L = 1
X       R1(L) = REAL(N1)
X       R2(L) = REAL(N2)
X       R(L)  = REAL(NCOMP)
X       ET1 = 0.0
X       ET2 = 0.0
X
XC
XC     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
XC     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
XC
X   10  IF(ID1(I) .NE. 0) THEN
X          IF(ID2(I) .EQ. 1) THEN
X             ET1 = ET1 + 1.0
X          ELSE
X             ET2 = ET2 + 1.0
X          ENDIF
X          I = I + 1
X          GOTO 10
X       ENDIF
XC
XC     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
XC     *     FINISHED.                                                      *
XC
X   20  D(L)  = 0.0
X       D1(L) = 0.0
X       D2(L) = 0.0
X       E(L)  = 0.0
X       E1(L) = 0.0
X       E2(L) = 0.0
X       TEMP  = XY(I)
XC
XC     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
XC     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
XC     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
XC     * APPROPRIATE GROUP.                                                *
XC
X
X  30   IF(ID1(I) .EQ. 0) THEN
X          IF(ID2(I) .EQ. 1) THEN
X            D1(L) = D1(L) + 1.0
X          ELSE
X             D2(L) = D2(L) + 1.0
X          ENDIF
X
X          D(L) = D1(L) + D2(L)
X
XC
XC     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
XC     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
XC
X        ELSE
X           IF(ID2(I) .EQ. 1) THEN
X              E1(L) = E1(L) + 1.0
X           ELSE
X              E2(L) = E2(L) + 1.0
X           ENDIF
X           E(L) = E1(L) + E2(L)
X        ENDIF
X
X       IF(I .LE. NCOMP) THEN
X           I = I + 1
XC
XC     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
XC     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
XC     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
XC     * THE NUMBER OF THE CENSORED DATA POINTS                             *
XC
X         IF(TEMP .EQ. XY(I)) GOTO 30
X         IF(ID1(I) .NE. 0) GOTO 30
X
XC
XC     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
XC
X            IF(L .EQ. 1) THEN
X                R1(L) = R1(L) - ET1
X                R2(L) = R2(L) - ET2
X                R(L)  = R1(L) + R2(L)
X            ELSE
X                R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
X                R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
X                R(L)  = R1(L) + R2(L)
X            ENDIF
X            L = L + 1
X            GOTO 20
X        ENDIF
XC
XC     *       COMPUTE THE SCORE AND VARIANCE                         *
XC
X
X        SCORE = 0.0
X        VAR   = 0.0
X        L1 = L - 1
X        DO 200 I = 1, L1
X
X           SCORE = SCORE+(D2(I)-(R2(I)*D(I))/R(I))
X
X           IF(R(I) .GT. 1.0) THEN
X              VAR = VAR + D(I)*(R2(I)/R(I))*
X     +              (1.0-(R2(I)/R(I)))*((R(I)-D(I))/(R(I)-1.0))
X           ENDIF
X
X  200   CONTINUE
X
XC
XC     *        NOW COMPUTE THE LOGRANK STATISTIC                   *
XC
X        TEST = SCORE/DSQRT(VAR)
X        PROB = 1.0 - AGAUSS(TEST)
X
X        RETURN
X        END
X
XC
XC      **********************************************************************
XC      ************************SUBROUTINE MULVAR ****************************
XC      **********************************************************************
XC
X       SUBROUTINE MULVAR(X,Y,IND,NTOT,ICOL,NVAR,NOTEST,IPROG,ICOMM,
X     +                   OUTPUT,COLM,FILE,YNAME,TITLE,ND,NYC,ICENS,
X     +                   NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,MVAR,
X     +                   LENG,LEGWRK,IBIN,XX,IND2,ALPHA,DWRK2,
X     +                   DWRK3,DWRK4,DWRK5,DWRK6,DWRK8,RWRK1,
X     +                   EWRK1,AWRK1,WWRK1,WWRK2,VWRK1,VWRK2,
X     +                   WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
X     +                   WRK9,WRK10,WRK11,WRK12,
X     +                   SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
X     +                   SWRK8,SWRK9,SWRK10,SWRK11,LWRK1,LWRK2,LWRK3,
X     +                   IWRK1,IWRK2,IWRK3,IWRK4,IWRK5,IWRK6,IWRK7,
X     +                   IWRK8,CWRK1,CWRK2,
X     +                   IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,IBWRK6,
X     +                   IBWRK7,IBWRK8,IBWRK9,BWRK1,BWRK2)
XC
XC      *                                                                    *
XC      *   THIS IS THE SUBROUTINE WHICH MANAGES CORRELATION/REGRESSION      *
XC      *   PROBLEMS.                                                        *
XC      *                                                                    *
XC      *   INPUT                                                            *
XC      *      FILE     :  NAME OF DATA FILE (9 LETTERS)                     *
XC      *      TITLE    :  TITLE OF THE PROBLEM (80 LETTERS)                 *
XC      *      NVAR     :  NUMBER OF VARIABLES                               *
XC      *      ICOL     :  INDICATOR OF VARIABLE (<=NVAR)                    *
XC      *                  IF A MULTIVARIATE PROBLEM IS NEEDED, SET ICOL=0   *
XC      *      COLM     :  NAME OF THE INDEPENDENT VARIABLE                  *
XC      *      YNAME    :  NAME OF THE DEPENDENT VARIABLE                    *
XC      *      COMMAND  :  NAME OF THE "COMMAND" FILE                        *
XC      *      OUTPUT   :  NAME OF THE OUTPUT FILE                           *
XC      *      IND(1,I) :  INDICATOR OF CENSORING                            *
XC      *                    IF =0,  DETECTED                                *
XC      *                       =1,  Y LOWER LIMIT                           *
XC      *                       =2,  X LOWER LIMIT                           *
XC      *                       =3,  DOUBLE LOWER LIMIT                      *
XC      *                       =4,  X UPPER LIMIT AND Y LOWER LIMIT         *
XC      *                       =5,  DATA POINT IS CONFINED BETWEEN TWO      *
XC      *                            VALUES                                  *
XC      *                  FOR THE UPPER LIMITS, CHANGE THE SIGN             *
XC      *                  2, 3, AND 4 CAN BE USED ONLY IN BHK AND SCHMITT'S *
XC      *                  5 CAN BE USED ONLY IN EM ALGORITHM AND IN         *
XC      *                  BINNING METHODS                                   *
XC      *      X(J,I)   :  INDEPENDENT VARIABLES                             *
XC      *      Y(I)     :  DEPENDENT VARIABLE                                *
XC      *     IPROG(I)  :  INDICATOR OF METHODS                              *
XC      *     NOTEST    :  NUMBERS OF TEST                                   *
XC      *  INPUT FOR EM ALGORITHM                                            *
XC      *      TOL      :  TOLERANCE (DEFAULT 1.0E-5)                        *
XC      *      MAX      :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *      IBET     :  IF 0, NO DEPENDENT VARIABLE IS CONFINED BETWEEN   *
XC      *                        TWO VALUES                                  *
XC      *                     1, THERE ARE SOME DEPENDENT VARIABLE WHICH     *
XC      *                        ARE CONFINED BETWEEN TWO VALUES             *
XC      *    ALPHA(K)   :  INITIAL ESTIMATE OF REGRESSION COEFFICIENTS       *
XC      *  INPUTS FOR BUCKLEY-JAMES METHOD                                   *
XC      *      TOL1     :  TOLERANCE (DEFAULT 1.0E-5)                        *
XC      *      MAX1     :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *  INPUTS FOR SCHMITT'S BINNING METHOD                               *
XC      *      MX       :  BIN NUMBER OF X AXES                              *
XC      *      MY       :  BIN NUMBER OF Y AXES                              *
XC      *      TOL3     :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                  *
XC      *      MAX3     :  MAXIMUM ITERATION (DEFAULT 20)                    *
XC      *      XBIN     :  BIN SIZE FOR X AXES                               *
XC      *      YBIN     :  BIN SIZE FOR Y AXES                               *
XC      *      XORG     :  ORIGN OF X AXES                                   *
XC      *      YORG     :  ORIGN OF Y AXES                                   *
XC      *      ISKIP    :  IF 0, THE PROGRAM WILL PROVIDE XBIN, YBIN, XORG,  *
XC      *                        AND YORG.                                   *
XC      *                    >0, THESE VALUES MUST BE PROVIDED BY THE USER   *
XC      *      IPIRNT   :  IF 0, NO TWO DIMENSIONAL K-M ESTIMATOR WILL BE    *
XC      *                        PRINTED                                     *
XC      *                    >0, TWO DIMENSIONAL K-M ESTIMATOR WILL BE       *
XC      *                        PRINTED                                     *
XC      *                                                                    *
XC      *    WORKING ARRAYS:                                                 *
XC      *      NTOT     :  NUMBER OF DATA POINTS                             *
XC      *      ND       :  NUMBER OF DETECTED POINTS                         *
XC      *      NC1      :  NUMBER OF Y LOWER LIMITS                          *
XC      *      NC2      :  NUMBER OF X LOWER LIMITS                          * 
XC      *      NC3      :  NUMBER OF DOUBLE LOWER LIMITS                     * 
XC      *      NC4      :  NUMBER OF Y LOWER AND X UPPER LIMITS              *
XC      *      NC5      :  NUMBER OF Y UPPER LIMITS                          *
XC      *      NC6      :  NUMBER OF X UPPER LIMITS                          *
XC      *      NC7      :  NUMBER OF DOUBLE UPPER LIMITS                     *
XC      *      NC8      :  NUMBER OF Y UPPER AND X LOWER LIMITS              *
XC      *      ICENS    :  IF 0, CENSORING IS MIXED                          *
XC      *                     1, CENSORING IS Y LOWER LIMITS ONLY            *
XC      *                    -1, CENSORING IS Y UPPER LIMITS ONLY            *
XC      *      NYC      :  NC1+NC2                                           *
XC      *      NXC      :  NC3+NC4                                           *
XC      *      NBC      :  NC5+NC6+NC7+NC8                                   *
XC      *      IDO      :  NXC+NBC                                           *
XC      *      IMUL     :  INDICATOR OF MULTIVARIATE PROBLEM                 *
XC      *      XX(J,I)  :  =X(ICOL,I), EXCEPT FOR MULTI INDEPENDENT VARIABLE *
XC      *                  CASE (J=1,NVAR).                                  *
XC      *      IND2(I)  :  =IND(1,I)                                         *
XC      *                                                                    *
XC      *  OUTPUT                                                            *
XC      *     COXREG                                                         *
XC      *      CHI      : GLOBAL CHI-SQUARE                                  *
XC      *      PROB     : PROBABILITY FOR NULL HYPOTHESIS                    *
XC      *     BHK                                                            *
XC      *       Z       : DEVIATION                                          *
XC      *      PROB     : PROBABILITY FOR NULL HYPOTHESIS                    *
XC      *     EM ALGORITHM                                                   *
XC      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS  (K=1,NVAR+1)       *
XC      *     ALPHA(K+2): STANDARD DEVIATION                                 *
XC      *     SIGMAA(K) : ERROR                                              *
XC      *     ITE       : NUMBER OF ITERATIONS                               *
XC      *     BUKLY-JAMES                                                    *
XC      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS (K=1,NVAR+1)        *
XC      *     ALPHA(K+2): STANDARD DEVIATION                                 *
XC      *     SIGMAA(K) : ERROR                                              *
XC      *     SCHMITT                                                        *
XC      *     ALPHA     : INTERCEPT COEFFICIENT                              *
XC      *     BETA      : SLOPE COEFFICIENT                                  *
XC      *                                                                    *
XC      *   SUBROUTINES                                                      *
XC      *     R3, R4, R5, R6, XDATA, COXREG, BHK, EM, BJ, TWOKM, SPRMAN      *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT),Y(NTOT),IND(MVAR,NTOT)
X       DIMENSION IND2(NTOT),IPROG(NOTEST),ALPHA(MVAR)
X
X       DIMENSION DWRK2(MVAR,NTOT),DWRK3(MVAR,NTOT)
X       DIMENSION DWRK4(MVAR,NTOT),DWRK5(MVAR,NTOT),DWRK6(MVAR,NTOT)
X       DIMENSION DWRK8(MVAR,NTOT)
X
X       DIMENSION LWRK1(MVAR,NTOT),LWRK2(MVAR,NTOT),LWRK3(MVAR,NTOT)
X
X       DIMENSION EWRK1(MVAR,MVAR),RWRK1(NTOT,MVAR)
X
X       DIMENSION AWRK1(5,IBIN)
X       DIMENSION WWRK1(LENG),WWRK2(LENG),VWRK1(LEGWRK),VWRK2(LEGWRK)
X
X       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT),WRK5(NTOT)
X       DIMENSION WRK6(NTOT),WRK7(NTOT),WRK8(NTOT),WRK9(NTOT),WRK10(NTOT)
X       DIMENSION WRK11(NTOT),WRK12(NTOT)
X
X       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR),SWRK4(MVAR)
X       DIMENSION SWRK5(MVAR),SWRK6(MVAR),SWRK7(MVAR),SWRK8(MVAR)
X       DIMENSION SWRK9(MVAR),SWRK10(MVAR),SWRK11(MVAR)
X
X       DIMENSION IWRK1(NTOT),IWRK2(NTOT),IWRK3(NTOT),IWRK4(NTOT)
X       DIMENSION IWRK5(NTOT),IWRK6(NTOT),IWRK7(NTOT),IWRK8(NTOT)
X
X       DIMENSION CWRK1(IBIN),CWRK2(IBIN)
X
X       DIMENSION IBWRK1(IBIN,IBIN),IBWRK2(IBIN,IBIN),IBWRK3(IBIN,IBIN)
X       DIMENSION IBWRK4(IBIN,IBIN),IBWRK5(IBIN,IBIN),IBWRK6(IBIN,IBIN)
X       DIMENSION IBWRK7(IBIN,IBIN),IBWRK8(IBIN,IBIN),IBWRK9(IBIN,IBIN)
X       DIMENSION BWRK1(IBIN,IBIN),BWRK2(IBIN,IBIN)
X
X       CHARACTER*9 FILE,OUTPUT,COLM,YNAME
X       CHARACTER*80 TITLE
XC
X       IMUL=1
X       IF(ICOL.EQ.0) IMUL=NVAR
X       IO=1
X       IF(OUTPUT.EQ.'         ') IO=0
XC
XC      *     READ A FEW MORE INPUTS FOR SPRMAN, EM, BJ, AND TWOKM          *
XC
X       DO 10 I=1,NOTEST
X          IF(IPROG(I).EQ.3)
X     +         CALL R3(IPRSP,ICOMM)
X          IF(IPROG(I).EQ.4) 
X     +         CALL R4(TOL,MAX,IBET,ICOL,ALPHA,ICOMM,MVAR)
X          IF(IPROG(I).EQ.5) 
X     +         CALL R5(TOL1,MAX1,ICOMM)
X          IF(IPROG(I).EQ.6) 
X     +         CALL R6(MX,MY,ISKIP,IPRINT,TOL2,MAX2,
X     +                 XBIN,YBIN,XORG,YORG,ICOMM,NLAST,IRAND)
X   10  CONTINUE
XC
XC      *    ADJUST THE INPUT DATA FORMAT                                    *
XC
X       CALL XDATA(X,XX,IND,IND2,IMUL,ICOL,NTOT,MVAR)
XC
X       IF(IO.EQ.0) THEN
X          WRITE(6,30)
X          WRITE(6,40)
X          WRITE(6,50) TITLE
X          WRITE(6,30)
X          WRITE(6,55) FILE
X          WRITE(6,30)
X       ELSE
X          WRITE(60,30)
X          WRITE(60,40)
X          WRITE(60,50) TITLE
X          WRITE(60,30)
X          WRITE(60,55) FILE
X          WRITE(60,30)
X       ENDIF
XC
X   30  FORMAT('     ')
X   40  FORMAT(5X,' CORRELATION AND REGRESSION PROBLEM')
X   50  FORMAT(5X,' TITLE IS  ',A80)
X   55  FORMAT(5X,' DATA FILE IS ',A9)  
XC
X   60  IF(IO.EQ.0) THEN
X          PRINT *
X          IF(ICOL.EQ.0) PRINT 80
X          IF(IMUL.EQ.1) THEN
X             PRINT 85
X             PRINT 90,COLM,YNAME
X          ENDIF
X          WRITE(6,30)
X       ELSE
X          WRITE(60,30)
X          IF(ICOL.EQ.0) WRITE(60,80)
X          IF(IMUL.EQ.1) THEN
X             WRITE(60,85)
X             WRITE(60,90) COLM,YNAME
X          ENDIF
X          WRITE(60,30)
X       ENDIF
XC
X   80  FORMAT(5X,'MULTIVARIATE PROBLEM')
X   85  FORMAT(6X,'INDEPENDENT',6X,' DEPENDENT')
X   90  FORMAT(8X,A9,' AND   ',A9)
XC
X  100  IF(IO.EQ.0) THEN
X          PRINT *
X          PRINT 110,NTOT
X          IF(ICENS.NE.-1) THEN
X             PRINT 120
X             PRINT 130,NC1,NC2,NC3,NC4
X          ELSEIF(ICENS.NE.1) THEN
X  142        PRINT 140
X             PRINT 130,NC5,NC6,NC7,NC8
X             PRINT 30
X          ENDIF
X       ELSE
X          WRITE(60,30)
X          WRITE(60,110) NTOT
X          IF(ICENS.NE.-1) THEN
X             WRITE(60,120)
X             WRITE(60,130) NC1,NC2,NC3,NC4
X          ELSEIF(ICENS.NE.1) THEN
X  102        WRITE(60,140)
X             WRITE(60,130) NC5,NC6,NC7,NC8
X          ENDIF
X          WRITE(60,30)
X       ENDIF
XC
X  110  FORMAT(5X,' NUMBER OF DATA POINTS : ',I5)
X  120  FORMAT(5X,' LOWER LIMITS IN  Y    X    BOTH   MIX')
X  130  FORMAT(19X,2I5,3X,2I5)
X  140  FORMAT(5X,' UPPER LIMITS IN  Y    X    BOTH   MIX')
XC
X       DO 200 J=1,NOTEST
XC
XC      *       CALL TESTS AND COMPUTE THE RESULTS                           *
XC       
X          IF(IPROG(J).EQ.1) THEN
X                    CALL COXREG(IND2,XX,Y,NTOT,IMUL,OUTPUT,ICENS,
X     +                          EWRK1,SWRK1,SWRK2,IWRK1,IWRK2,
X     +                          IWRK3,WRK1,DWRK8,SWRK3,IWRK4,IWRK5,MVAR)
X
X          ELSEIF(IPROG(J).EQ.2) THEN
X                    CALL BHK(IND2,XX,Y,NTOT,OUTPUT,
X     +                           WRK1,WRK2,IWRK1,IWRK2,IWRK3,MVAR)
X
X          ELSEIF(IPROG(J).EQ.3) THEN     
X              CALL SPRMAN(IND2,XX,Y,NTOT,OUTPUT,IPRSP,MVAR,
X     +                          WRK1,IWRK1,IWRK2,LWRK1,LWRK2,LWRK3,
X     +                          DWRK8,DWRK2,DWRK3,DWRK4,WRK3,WRK4,
X     +                          WRK5,WRK6,WRK7,WRK8,WRK9,WRK10,WRK11,
X     +                          WRK12,IWRK3,IWRK4,DWRK5,DWRK6,SWRK1)
X
X          ELSEIF(IPROG(J).EQ.4) THEN
X                    CALL EM(IND2,XX,Y,NTOT,TOL,MAX,IMUL,IBET,
X     +                      ND,NYC,OUTPUT,FILE,ALPHA,
X     +                      RWRK1,WRK1,WWRK1,WWRK2,VWRK1,VWRK2,
X     +                      EWRK1,SWRK1,LENG,LEGWRK,MVAR)
X
X          ELSEIF(IPROG(J).EQ.5) THEN
X                    CALL BJ(IND2,XX,Y,NTOT,TOL1,MAX1,IMUL,ND,
X     +                      NYC,ICENS,OUTPUT,
X     +                      SWRK1,SWRK2,IWRK1,IWRK2,IWRK4,
X     +                      IWRK5,IWRK6,IWRK7,IWRK8,WRK1,WRK2,WRK3,
X     +                      WRK4,WRK5,WRK6,WRK7,WRK8,SWRK3,SWRK4,
X     +                      SWRK5,SWRK6,SWRK7,SWRK8,SWRK9,SWRK10,
X     +                      SWRK11,DWRK8,EWRK1,MVAR)
X
X          ELSEIF(IPROG(J).EQ.6) THEN
X                     CALL TWOKM(IND2,XX,Y,NTOT,MX,MY,ISKIP,IPRINT,
X     +                          ICENS,XBIN,YBIN,XORG,YORG,OUTPUT,
X     +                          TOL2,MAX2,NLAST,IRAND,
X     +                          NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,
X     +                          WRK1,WRK2,IWRK1,CWRK1,CWRK2,BWRK1,
X     +                          IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
X     +                          IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK2,
X     +                          IWRK2,IWRK3,WRK3,WRK4,SWRK1,DWRK8,
X     +                          AWRK1,IBIN,MVAR)
X
X          ENDIF
X  200  CONTINUE
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      *********************** FUNCTION PCHISQ  *****************************
XC      **********************************************************************
XC
X       FUNCTION PCHISQ(CHISQR,NFREE)
XC
XC      *      THIS FUNCTION IS BASED ON PHILIP R. BEVINGTON 1969, "DATA     *
XC      *      REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", 1969,*
XC      *      McGRAW HILL (NY:NY), PROGRAM 10-1 P. 192,                     *
XC      *      AND COMPUTES CHI-SQUARE PROBABILITY FROM THE REDUCED          *
XC      *      CHI-SQUARE.                                                   *
XC      *      INPUT :      CHISQR : REDUCED CHI-SQUARE                      *
XC      *                   NFREE  : DEGREES OF FREEDOM                      *
XC      *     OUTPUT :      PCHISQ : CHI-SQUARE PROBABILITY                  *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
XC
X       IF(NFREE .LE. 0) THEN
X          PCHISQ=0.0
X          RETURN
X       ENDIF
X
X       FREE=NFREE
X       Z=CHISQR*FREE/2.0
X       NEVEN=2*(NFREE/2)
X       IF((NFREE-NEVEN) .LE. 0) THEN
XC
XC      *          THE DEGREES OF FREEDOM ARE EVEN                          *
XC
X          IMAX=NFREE/2
X          TERM=1.0
X          SUM=0.0
X          DO 34 I=1,IMAX
X             FI=I
X             SUM=SUM+TERM
X             TERM=TERM*Z/FI
X   34     CONTINUE
X          PCHISQ=SUM*DEXP(-Z)
X          RETURN
XC
XC      *           THE DEGREES OF FREEDOM ARE ODD                          *
XC
X       ELSE
X          IF((Z-25.0) .GT. 0) THEN
X             Z=CHISQR*(FREE-1.0)/2.0
X             IMAX=NFREE/2
X             TERM=1.0
X             SUM=0.0
X             DO 44 I=1,IMAX
X                FI=I
X                SUM=SUM+TERM
X                TERM=TERM*Z/FI
X   44        CONTINUE
X             PCHISQ=SUM*DEXP(-Z)
X             RETURN
X          ELSE
X             PWR=FREE/2.0
X             TERM=1.0
X             SUM=TERM/PWR
X             DO 56 I=1,1000
X                FI=I
X                TERM=-TERM*Z/FI
X                SUM=SUM+TERM/(PWR+FI)
X                IF((DABS(TERM/SUM)-0.00001) .LE. 0.0) GOTO 57
X   56        CONTINUE
X   57        PCHISQ=1.0-(Z**PWR)*SUM/GAMMA(PWR)
X          ENDIF
X       ENDIF
X
X       RETURN
X       END
XC     
XC
XC***************************************************************************
XC**************************  SUBROUTINE PETO   *****************************
XC***************************************************************************
XC
XC
X       SUBROUTINE PETO(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R,D1,E1,R1,D2,
X     +                 E2, R2,F,A,FT,AT,XA,SCORE,VAR)
XC     *
XC     * THIS SUBROUTINE COMPUTES THE PETO-PRENTICE STATISTIC USING THE      *
XC     * FORMULATION IN LATTA, 'A MONTE-CARLO STUDY OF SOME TWO-SAMPLE RANK  *
XC     * TESTS WITH CENSORED DATA', 1981, JOURNAL OF THE AMERICAN STATISTICAL*
XC     * ASSOCIATION, VOL 76, PP 713-719.  THE FORM USED IS FROM EQUATION    *
XC     * 2.2 AND THE ASYMPTOTIC VARIANCE ESTIMATE GIVEN IN THE ARTICLE IS    *
XC     * USED FOR THE VARIANCE.                                              *
XC     *                                                                     *
XC     * INPUT                                                               *
XC     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
XC     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
XC     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
XC     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
XC     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
XC     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
XC     *                                                                     *
XC     * OUTPUT                                                              *
XC     *     TEST    : STANDARDIZED PETO-PRENTICE STATISTIC                  *
XC     *     PROB    : PROBABILITY                                           *
XC     *                                                                     *
XC     * OTHERS                                                              *
XC     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
XC     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
XC     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
XC     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
XC     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
XC     *      R(I)   : RISK SET AT XY(I)                                     *
XC     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
XC     *      F(I)   : THE ESTIMATE OF THE SURVIVAL FUNCTION AT XY(I)        *
XC     *      A(I)   : COEFFICIENT AT XY(I)                                  *
XC     *      XA(I)  : SUM OF 2 X D2(I) AND E2(I)                            *
XC     *      SCORE  : SCORE OF THE DATA                                     *
XC     *      VAR    : VARIANCE OF THE DATA                                  *
X
X
X
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X
X       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
X       DIMENSION F(NTOT),A(NTOT),FT(NTOT),AT(NTOT)
X       DIMENSION D(NTOT),E(NTOT),R(NTOT),XA(NTOT)
X       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
X       DIMENSION E2(NTOT),R2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X
X       J = 0
X       I = 1
X       L = 1
X       F(I)  = 1.0
X       A(I)  = 1.0
X       R1(L) = REAL(N1)
X       R2(L) = REAL(N2)
X       R(L)  = REAL(NCOMP)
X       ET1 = 0.0
X       ET2 = 0.0
X
XC
XC     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
XC     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
XC
X   10  IF(ID1(I) .NE. 0) THEN
X          IF(ID2(I) .EQ. 1) THEN
X             ET1 = ET1 + 1.0
X          ELSE
X             ET2 = ET2 + 1.0
X          ENDIF
X          I = I + 1
X          F(I) = 1.0
X          A(I) = 1.0
X          GOTO 10
X       ENDIF
XC
XC     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
XC     *     FINISHED.                                                      *
XC
X   20  D(L)  = 0.0
X       D1(L) = 0.0
X       D2(L) = 0.0
X       E(L)  = 0.0
X       E1(L) = 0.0
X       E2(L) = 0.0
X       TEMP  = XY(I)
X       K = 0
XC
XC     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
XC     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
XC     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
XC     *  APPROPRIATE GROUP.                                                *
XC     * HERE, FT AND AT ARE USED, SINCE WE ASSUME THAT THERE ARE NO TIES.  *
XC     * IF THERE ARE TIES IN THE DATA, FT AND AT WILL BE APPROPRIATELY     *
XC     *  CONVERTED INTO THE FORM FOR TIED DATA AND PUT IN F AND A.         *
XC
X
X  30   IF(ID1(I) .EQ. 0) THEN
X         IF(ID2(I) .EQ. 1) THEN
X            D1(L) = D1(L) + 1.0
X         ELSE
X            D2(L) = D2(L) + 1.0
X         ENDIF
X
X         D(L) = D1(L) + D2(L)
X         J = J + 1
X         K = K + 1
X
X         IF(L .EQ. 1) THEN
X           RISK = R(L) - (ET1+ET2) - (D(L) - 1.0)
X           IF(J .EQ. 1) THEN
X              FT(J) = RISK/(RISK + 1.0)
X              AT(J) = (RISK + 1.0)/(RISK + 2.0)
X           ELSE
X              FT(J) = FT(J-1)*RISK/(RISK+1.0)
X              AT(J) = AT(J-1)*(RISK+1.0)/(RISK+2.0)
X           ENDIF
X         ELSE 
X           RISK = (R(L-1)-D(L-1))-E(L-1)-(D(L)-1.0)
X           FT(J) = FT(J-1)*RISK/(RISK + 1.0)
X           AT(J) = AT(J-1)*(RISK + 1.0)/(RISK + 2.0)
X         ENDIF
X
XC
XC     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
XC     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
XC
X        ELSE
X          IF(ID2(I) .EQ. 1) THEN
X             E1(L) = E1(L) + 1.0
X          ELSE
X             E2(L) = E2(L) + 1.0
X          ENDIF
X          E(L) = E1(L) + E2(L)
X        ENDIF
X
X        IF(I .LE. NCOMP) THEN
X        I = I + 1
XC
XC     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
XC     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
XC     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
XC     * THE NUMBER OF THE CENSORED DATA POINTS                             *
XC
X        IF(TEMP .EQ. XY(I)) GOTO 30
X        IF(ID1(I) .NE. 0) GOTO 30
X
XC
XC     * IF THE DATA POINTS WERE TIED, K > 1.  COMPUTE THE AVERAGE OF       *
XC     * FT AND AT BETWEEN JJ= 1, K AND USE THE AVERAGES FOR F AND A OF THE *
XC     * DATA POINT.                                                        *
XC
X         SUM1 = 0.0
X         SUM2 = 0.0
X         JSTART = J - K + 1
X         DO 50 JJ = JSTART, J
X            SUM1 = SUM1 + FT(JJ)
X            SUM2 = SUM2 + AT(JJ)
X   50    CONTINUE
X
X         F(L)  = SUM1/FLOAT(K)
X         A(L)  = SUM2/FLOAT(K)
X         XA(L) = 2.0*D2(L) + E2(L)
XC
XC     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
XC
X         IF(L .EQ. 1) THEN
X             R1(L) = R1(L) - ET1
X             R2(L) = R2(L) - ET2
X             R(L)  = R1(L) + R2(L)
X         ELSE
X             R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
X             R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
X             R(L)  = R1(L) + R2(L)
X         ENDIF
X         L = L + 1
X         GOTO 20
X      ENDIF
XC
XC     *       COMPUTE THE SCORE AND VARIANCE                         *
XC
X
X      SCORE = 0.0
X      VAR   = 0.0
X      L1 = L - 1
X      DO 200 I = 1, L1
X
X         SCORE = SCORE +(2.0*F(I)-1.0)*D2(I)+(F(I)-1.0)*E2(I)
X
X         SUM = 0.0
X         JSTART = I + 1
X         DO 100 J = JSTART, L
X            SUM = SUM + F(J)*XA(J)
X  100    CONTINUE
X
X         VAR = VAR + F(I)*(1.0 - A(I))*XA(I)
X     +               - (A(I) - F(I))*XA(I)*(F(I)*XA(I) + 2.0*SUM)
X
X  200    CONTINUE
X
XC
XC     *        NOW COMPUTE THE PETO-PRENTICE STATISTIC                   *
XC
X        TEST = SCORE/DSQRT(VAR)
X        PROB = 1.0 - AGAUSS(TEST)
X
X        RETURN
X        END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE PLESTM  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE PLESTM(U,C,NU,NC,S,V,NTOT,SMEAN,SIGMA,ICHANGE,
X     +                   NCHANGE,L)
XC       
XC      *      THIS SUBROUTINE COMPUTES PL ESTIMATOR AND THE MEAN            *
XC      *      AND ITS ERROR.                                                *
XC      *                                                                    *
XC      *       INPUT     U : UNCENSORED DATA POINTS                         *
XC      *                 C : CENSORED DATA POINTS                           *
XC      *                NU : NO. OF UNCENSORED DATA POINTS                  *
XC      *                NC : NO. OF CENSORED DATA POINTS                    *
XC      *               NTOT: TOTAL NUMBER OF DATA POINTS                    *
XC      *                                                                    *
XC      *       WORK      L : RANK OF THE UNCENSORED DATA                    *
XC      *               VAR : VARIANCE OF THE MEAN                           *
XC      *                KD : NUMBER OF TIED DATA POINTS                     *
XC      *                                                                    *
XC      *       OUTPUT    S : PL ESTIMATOR                                   *
XC      *                 V : ERROR FOR THE PL ESTIMATOR                     *
XC      *             SMEAN : MEAN OF THE DATA                               *
XC      *             SIGMA : ERROR OF THE MEAN                              *
XC      *            ICHANGE: IF THE LAST VALUE IS CENSORED, WE NEED TO      *
XC      *                     CHANGE IT TO A DETECTION. THEN ICHANGE=-1,     *
XC      *                     OTHERWISE ICHANGE=1.                           *
XC      *            NCHANGE: IF ICHANGE = -1 AND THE LAST VALUE IS TIED     *
XC      *                     WITH OTHER CENSORED VALUES, THIS RECORDS THE   *
XC      *                     NUMBER OF TIED OBSERVATIONS (ALL OF THEM NEED  *
XC      *                     TO BE CHANGED TO DETECTIONS).                  *
XC      *                                                                    *
XC      *       FIRST HALF OF THE PROGRAM IS FROM ELISA T. LEE, "STATISTICAL *
XC      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
XC      *       PUBLICATIONS (BELMONT:CA); WITH THE GRAPHIC ROUTINES REMOVED.*
XC      *       FORMULAS USED FOR COMPUTATION OF THE MEAN AND THE ERROR ARE  *
XC      *       FROM RUPERT G. MILLER, "SURVIVAL ANALYSIS", 1981,            *
XC      *       JOHN WILEY & SONS (NY:NY)                                    *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION U(NTOT),C(NTOT),S(NTOT),V(NTOT),L(NTOT)
XC
XC      *          COMPUTE THE RANK (L) OF THE UNCENSORED POINTS             *
XC
XC*******     IF THE LAST VALUE IS CENSORED, CHANGE IT TO A DETECTION        *
XC
X
XC      THE FOLLOWING LOOP HAS BEEN MODIFIED AND NCHANGE ADDED TO THE 
XC      PROGRAM TO COVER THE CASE WHEN TIED NONDETECTIONS ARE THE LARGEST
XC      VALUE.  MODIFIED 4/92
X
X       ICHANGE=1
X       NCHANGE = 0
X 13    IF(NC .NE. 0)THEN 
X          IF(U(NU) .LE. C(NC))THEN 
X             U(NU+1)=C(NC)
X             NU=NU+1
X             NC=NC-1
X             NCHANGE = NCHANGE + 1
X             ICHANGE=-1
X          ELSE
X             GOTO 15
X          ENDIF
X       ELSE
X          GOTO 15
X       ENDIF
X       GOTO 13
XC
X 15    K=1
X       KK=0
X       NT=NU+NC
X       IF(NC .NE. 0) THEN 
X          DO 10 I=1,NU
X             IF(KK .NE. NC) THEN
X                DO 20 J=K,NC
X                   K1=J
X                   IF(C(J) .GE. U(I)) GOTO 1
X                   KK=KK+1
X  20            CONTINUE
X             ENDIF
X   1         K=K1
X             L(I)=I+KK
X  10      CONTINUE
X       ELSE
X          DO 19 I=1,NU
X             L(I)=I
X  19      CONTINUE
X       ENDIF
XC
XC      *       COMPUTE P(T) FOR ALL UNCENSORED POINTS BASED ON RANK (L)     *
XC
X       V1=0.0
X       PT=1.0
X       XNT=NT
X       DO 12 I=1,NU
X          XL=L(I)
X          PT=PT*((XNT-XL)/(XNT-XL+1.0))
X          S(I)=PT
X          IF((XNT-XL) .LE. 0.0) THEN
X             VV=0.0
X          ELSE
X             V1=V1+1.0/((XNT-XL)*(XNT-XL+1.0))
X             VV=DSQRT((PT**2)*V1)
X          ENDIF
X          V(I)=VV
X  12   CONTINUE
X
XC
XC      *        COMPUTE THE MEAN                                            *
XC      *        REF. FOR THE MEAN AND ERROR :                               *
XC      *          MILLER, R. G. JR. 1981, "SURVIVAL ANALYSIS"               *
XC      *          PP. 70-71 AND 195-198.                                    *
XC
X       SMEAN=U(1)
X       I=2
X  30   K=0
X  40   IF((U(I+K).NE.U(I-1)).AND.(I+K.LE.NU)) THEN
X          SMEAN=SMEAN+S(I+K-1)*(U(I+K)-U(I-1))
X          IF(I+K.LT.NU) THEN
X             I=I+K+1
X             GOTO 30
X          ENDIF
X       ELSEIF(U(I+K).EQ.U(I-1)) THEN
X          K=K+1
X          GOTO 40
X       ENDIF
XC
XC      *              COMPUTE THE ERROR OF THE MEAN                         *
XC
X       J=2    
X       VAR=0.0
X  120  I=J
X       SSUM=0
X  130  K=0
X  140  IF((U(I+K).EQ.U(I-1)).AND.(I+K.LE.NU)) GOTO 145
X          IF(U(I+K).EQ.U(I-1)) THEN
X             K=K+1
X             GOTO 140
X          ENDIF
X  145     SSUM=SSUM+S(I+K-1)*(U(I+K)-U(I-1))
X          IF(I+K.LT.NU) THEN
X             I=I+K+1
X             GOTO 130
X          ENDIF
XC
XC      *          KD IS NO. OF TIED OBSERVATIONS AT THAT POINT              *
XC
X       KD=1
X  180  IF(U(J-1+KD).LE.U(J-1)) THEN
X          KD=KD+1
X          GOTO 180
X       ENDIF
X       XL=L(J-1)
X       D=KD
X       B=XNT-XL-D+1.0
XC
XC      *       IF THE LAST FEW DATA POINTS ARE UNCENSORED AND TIED, SKIP    *
XC      *       THE NEXT LINES TO AVOID DIVISION BY 0.                       *
XC
X       IF(B .NE. 0.0) THEN
X          VAR=VAR+SSUM*SSUM*D/((XNT-XL+1)*B)
X          J=J+KD
X          IF(J.LE.NU) GOTO 120
X       ENDIF
X  200  SIGMA=DSQRT(VAR)
X   
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE PWLCXN  ******************************
XC      **********************************************************************
XC
X       SUBROUTINE PWLCXN(H,XM,SCORE,TEST,PROB,IWLCX,XY,ID1,ID2,NTOT)
XC
XC      *           THIS SUBROUTINE COMPUTES PETO AND PETO'S                 *
XC      *           GENERALIZED WILCOXON STATISTIC.                          *
XC      *                                                                    *
XC      *           OBTAINED FROM ELISA T. LEE, "STATISTICAL METHODS FOR     *
XC      *           SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING         *
XC      *           PUBLICATIONS (BELMONT:CA)                                *
XC      *                                                                    *
XC      * SUBROUTINES                                                        *
XC      *           STAT                                                     *
XC
XC*******           COMMON STATEMENT IS DIFFERENT FROM SMSDA.                *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION H(NTOT),XM(NTOT),SCORE(NTOT)
X       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
X       COMMON /G/ N,N1,N2,NCEN,ISIGN,IFULL,LO
XC
X       IWLCX=0
X       IF(NCEN.EQ.0) IWLCX=1
X       IF(NCEN.EQ.0) RETURN    
X       L=1
X       I=1
X       IJK=0
XC
XC*******       THE NEXT LINE IS CHANGED FROM "EQ.1" TO "EQ.0".              *
XC
X  63   IF(ID1(I).NE.0) THEN
X          IF(IJK.EQ.1) GOTO 65
X          SCORE(I)=H(1)-1.0
X          IF(I.EQ.N) GOTO 200
X          I=I+1
X          GOTO 63
X       ENDIF
X  62   M=INT(XM(L))
X
X       DO 64 J=1,M
X          SCORE(I)=H(L)+H(L+1)-1.0
X          IF(I.EQ.N) GOTO 200
X          I=I+1
X  64   CONTINUE
X
X       IJK=1
X       L=L+1
X       GOTO 63
X
X  65   SCORE(I)=H(L)-1.0
X       IF(I.EQ.N) GOTO 200
X       I=I+1
X       GOTO 63
XC
XC*******      THE NEXT LINE IS ADDED. ALSO THE PRINTING COMMANDS            *
XC*******      ARE CHANGED.                                                  *
XC
X 200   CALL STAT(SCORE,TEST,XY,ID1,ID2,NTOT)
X       PROB=1.0-AGAUSS(TEST)
X       RETURN
X       END
X
XC****************************************************************************
XC************************** SUBROUTINE R3    ********************************
XC****************************************************************************
X
X       SUBROUTINE R3(IPRSP,ICOMM)
XC
XC      *      THIS SUBROUTINE READS ONE SUPPLIMENTAL INPUT FOR SPEARMAN'S   *
XC      *    RHO COMPUTATION.                                                *
XC      *                                                                    *
XC      *    SUBROUTINES                                                     *
XC      *              DATA2                                                 *
XC
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X       CHARACTER*1 CHECK,CHAR(4,4)
X
X       IF(ICOMM .NE.1) THEN
X          PRINT *
X          PRINT *,'ONE MORE QUESTION FOR THE SPEARMANS RHO COMPUTATION'
X          PRINT *
X
X   10     WRITE(6, 12)
X   12     FORMAT('DO YOU WANT TO PRINT OUT THE RANKS (Y/N)? ')
X     
X          READ(5, 14) CHECK
X   14     FORMAT(A1)
X
X           IF(CHECK .EQ. 'Y' .OR. CHECK .EQ. 'y') THEN
X              IPRSP = 1
X           ELSEIF(CHECK .EQ. 'N' .OR. CHECK .EQ. 'n') THEN
X              IPRSP = 0
X           ELSE 
X              GOTO 10
X           ENDIF
X 
X        ELSE 
X           READ(50,60) (CHAR(I,1), I = 1,4)
X   60      FORMAT(4A1)
X           CALL DATA2(CHAR,1,1,IPRSP,LIND)
X           IF((LIND.NE.0).OR.((IPRSP.NE.0).AND.(IPRSP.NE.1))) THEN
X              PRINT *
X              PRINT *,'INFORMATION ABOUT RANK PRINTOUT IS NOT CLEAR'
X              PRINT *
X              STOP
X           ENDIF
X        ENDIF
X        RETURN
X        END
X
XC
XC      **********************************************************************
XC      *********************** SUBROUTINE R4  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE R4(TOL,MAX,IBET,ICOL,ALPHA,ICOMM,MVAR)
XC
XC      * THIS SUBROUTINE IS A SUPPLEMENTAL READ-IN PROGRAM FOR THE EM       *
XC      * ALGORITHM WITH A NORMAL DISTRIBUTION                               *
XC      *                                                                    *
XC      *  INPUT   ICOL   :  NUMBER OF INDEPENDENT VARIABLES                 *
XC      *          ICOMM  :  INDICATOR OF COMMAND FILE VS. TERMINAL INPUT.   *
XC      *                                                                    *
XC      *  OUTPUT  TOL    :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                *
XC      *          MAX    :  MAXIMUM ITERATION (DEFAULT 20)                  *
XC      *          IBET   :  INDICATOR OF DATA SET TYPE (WHETHER THERE ARE   *
XC      *                    DATA POINTS WHICH ARE CONFINED BETWEEN TWO      *
XC      *                    VALUES)                                         *
XC      *         ALPHA(J):  INITIAL ESTIMATES OF REGRESSION COEFFICIENTS    *
XC      *                                                                    *
XC      *  SUBROUTINES                                                       *
XC      *         DATA1, DATA2                                               *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION ALPHA(MVAR)
X       CHARACTER*1 CHECK,CHAR(4,4)
XC  
X       ICOL3=ICOL+2
X       DO 5 I=1,ICOL3
X          ALPHA(I)=0.0
X    5  CONTINUE
XC
XC      *             FOR COMMAND FILE READING, GO TO 230                    *
XC
X       IF(ICOMM.NE.1) THEN
X          PRINT *
X          PRINT *,'A FEW MORE INPUTS FOR EM ALGORITHM'  
XC
XC      *              TOLERANCE LEVEL                                       *
XC      
X          TOL=1.0E-5
X   10     PRINT *
X          PRINT *,'DO YOU WANT TO SET'
X          WRITE(6,20)
X   20     FORMAT('    TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')
X          READ(5,30) CHECK
X   30     FORMAT(A1)
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X             PRINT *
X             WRITE(6,50)
X   50   FORMAT('WHAT IS THE TOLERANCE LEVEL (GIVE IN E FORMAT) ? ')
X             READ(5,60) TOL
X   60        FORMAT(E9.3)
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 10
X          ENDIF
XC
XC      *        INFORMATION ABOUT IBET                                      *
XC
XC   70     PRINT *
XC          PRINT *,'ARE THERE DATA POINTS WHICH ARE'
XC          WRITE(6,80)
XC   80     FORMAT('    CONFINED BETWEEN TWO VALUES  (Y/N) ? ')
XC          READ(5,30) CHECK
X
XC          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
XC             IBET=1
XC          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
X             IBET=0
XC          ELSE
XC             GOTO 70
XC          ENDIF
XC
XC      *          INITIAL ESTIMATIONS OF REGRESSION COEFFICIENTS            *
XC
X   90     PRINT *
X          PRINT *,'DO YOU HAVE INITIAL ESTIMATES'
X          WRITE(6,100)
X  100     FORMAT('    FOR THE REGRESSION COEFFICIENTS (Y/N) ? ')
X          READ(5,30) CHECK
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X  110        PRINT *
X             WRITE(6,120)
X  120        FORMAT('INTERCEPT COEFFICIENT= ')
X             READ(5,130) ALPHA(1)
X  130        FORMAT(F10.3)
X
X             IF(ICOL.EQ.1) THEN
X                WRITE(6,132)
X  132           FORMAT('SLOPE COEFFICIENT= ')
X                READ(5,130) ALPHA(2)
X
X             ELSE
X  138           ICOL2=ICOL+1
X                DO 170 I=2,ICOL2
X                   PRINT *
X                   IF(I.EQ.1) WRITE(6,140)
X                   IF(I.EQ.2) WRITE(6,150)
X                   IF(I.GE.3) WRITE(6,160) I
X
X  140           FORMAT('SLOPE COEFFICIENT FOR FIRST VARIABLE = ')
X  150           FORMAT('SLOPE COEFFICIENT FOR SECOND VARIABLE = ')
X  160         FORMAT('SLOPE COEFFICIENT FOR ',I4,'-TH VARIABLE = ')
X
X                   READ(5,130) ALPHA(I)
X  170           CONTINUE
X             ENDIF
X  171        ALPHA(ICOL+2)=1.0
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 90
X          ENDIF
XC
XC      *            ITERATION LIMITS                                        *
XC
X  180     PRINT *
X          PRINT *,'DO YOU WANT TO SET THE ITERATION '
X          WRITE(6,190)
X  190     FORMAT('     LIMIT (DEFAULT 50) (Y/N) ? ')
X          READ(5,30) CHECK
X          MAX=50
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X  200        PRINT *
X             WRITE(6,210)
X  210        FORMAT('WHAT IS THE MAXIMUM ITERATION ? ')
X             CALL DATA1(MAX)
X             IF(MAX.LE.0) GOTO 200
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 180
X          ENDIF
XC
XC
XC      *            READING FROM "COMMAND" FILE                             *
XC
XC      *            TOLERANCE LEVEL                                         *
XC
X       ELSE
X  230     READ(50,60) TOL
X          IF(TOL.LE.0.0) THEN
X             PRINT *
X             PRINT *,'    TOLERANCE FOR THE EM ALGORITHM IS NEGATIVE.'
X             STOP
X          ENDIF
XC
XC      *                   IBET                                             *
XC
XC  235     READ(50,240) (CHAR(I,1),I=1,4)
X  240     FORMAT(4A1)
XC          CALL DATA2(CHAR,1,1,IBET,LIND)
X
XC          IF((LIND.NE.0) .OR. ((IBET.NE.0).AND.(IBET.NE.1))) THEN
XC             PRINT *
XC             PRINT *,'     IBET INDICATOR FOR CONFINED POINTS IS WRONG'
XC             STOP
XC          ENDIF
X          IBET = 0
XC
XC      *      INITIAL ESTIMATES FOR REGRESSION COEFFICIENTS                 *
XC
X  245     ICOL2=ICOL+1
X          ICOL3=ICOL+2
X          READ(50,250) (ALPHA(I),I=1,ICOL2), ALPHA(ICOL3)
X  250     FORMAT(12F10.3)
XC
XC      *      MAXIMUM ITERATION                                             *
XC
X          READ(50,240) (CHAR(I,1),I=1,4)
X          CALL DATA2(CHAR,1,1,MAX,LIND)
X
X          IF((LIND.NE.0) .OR. (MAX.LE.0) .OR. (MAX.GT.1000)) THEN
X  252        PRINT *
X          PRINT *,'  MAXIMUM ITERATION VALUE IS NOT BETWEEN 1 AND 1000.'
X             STOP
X          ENDIF
X
X       ENDIF
X  260  RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE R5 *********************************
XC      **********************************************************************
XC
X       SUBROUTINE R5(TOL,MAX,ICOMM)
XC
XC      * THIS SUBROUTINE IS A SUPPLIMENTAL READING PROGRAM FOR THE BUCKLEY- *
XC      * JAMES METHOD.                                                      *
XC      *                                                                    *
XC      * INPUT       ICOMM : INDICATOR OF READING METHOD                    *
XC      *                                                                    *
XC      * OUTPUT       TOL  : TOLERANCE (DEFAULT 1.0E-5)                     *
XC      *              MAX  : MAXIMUM ITERATION (DEFAULT 50)                 *
XC      * SUBROUTINES                                                        *
XC      *              DATA1, DATA2                                          *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*1 CHECK,CHAR(4,4)
XC
X       TOL=1.0E-5
X       MAX=50
XC
XC      *         FOR "COMMAND" FILE, GO TO 230                              *
XC
X       IF(ICOMM.NE.1) THEN
X          PRINT *
X          PRINT *,' A FEW MORE INPUTS FOR THE BUCKLEY-JAMES METHOD'
XC
XC      *            TOLERANCE LEVEL                                        *
XC
X   10     PRINT *
X          PRINT *,'DO YOU WANT TO SET '
X          WRITE(6,20)
X   20     FORMAT('    TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')
X
X          READ(5,30) CHECK
X   30     FORMAT(A1)
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X   40        PRINT *
X             WRITE(6,50)
X   50        FORMAT('WHAT IS THE TOLERANCE LEVEL (E9.3 FORMAT) ? ')
X             READ(5,60) TOL
X   60        FORMAT(E9.3)
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 10
X          ENDIF
XC
XC      *            ITERATION LIMITS                                       *
XC
X  180     PRINT *
X          PRINT *,'DO YOU WANT TO SET '
X          WRITE(6,190)
X  190     FORMAT('   ITERATION LIMIT (DEFAULT 50) (Y/N) ? ')
X          READ(5,30) CHECK
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X  200        PRINT *
X             WRITE(6,210)
X  210        FORMAT('WHAT IS THE MAXMUM ITERATION ? ')
X             CALL DATA1(MAX)
X             IF(MAX.LE.0) GOTO 200
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 180
X          ENDIF
XC
XC
XC      *      "COMMAND" FILE READING                                        *
XC
XC      *          TOLERANCE LEVEL                                           *
XC    
X       ELSE
X  230     READ(50,60) TOL
X          READ(50,240) (CHAR(I,1),I=1,4)
X  240     FORMAT(4A1)
XC
XC      *         MAXIMUM ITERATION                                          *
XC
X          CALL DATA2(CHAR,1,1,MAX,LIND)
X          IF((LIND.NE.0) .OR. (MAX.LE.0) .OR. (MAX.GT.1000)) THEN
X  242        PRINT *
X          PRINT *,'   MAXIMUM ITERATION VALUE IS NOT BETWEEN 1 AND 1000'
X             STOP
X          ENDIF
X       ENDIF
X  
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      *********************** SUBROUTINE R6  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE R6(MX,MY,ISKIP,IPRINT,TOL,MAX,XBIN,
X     +                            YBIN,XORG,YORG,ICOMM,NLAST,IRAND)
XC     
XC      *  THIS PROGRAM IS A SUPPLEMENTAL READING PROGRAM FOR SCHMITT`S      *
XC      *  BINNED LINEAR REGRESSION METHOD.                                  *
XC      *                                                                    *
XC      *  INPUT   ICOMM : INDICATOR OF READING METHOD                       *
XC      *                                                                    *
XC      *  OUTPUT  MX    : NUMBER OF BINS ON X AXIS                          *
XC      *          MY    : NUMBER OF BINS ON Y AXIS                          *
XC      *          ISKIP : INDICATOR WHETHER THE USER GIVES BINNING          *
XC      *                  INFORMATION                                       *
XC      *          XBIN  : BIN SIZE FOR X AXIS                               *
XC      *          YBIN  : BIN SIZE FOR Y AXIS                               *
XC      *          XORG  : ORIGIN OF X AXIS                                  *
XC      *          YORG  : ORIGIN OF Y AXIS                                  *
XC      *          IPRINT: INDICATOR OF PRINTING FOR TWO DIM. KM ESTIMATOR   *
XC      *                                                                    *
XC      * SUBROUTINES                                                        *
XC      *          DATA1, DATA2                                              *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*1 CHECK,CHECK1,CHAR(4,4)
XC
X       TOL=1.0E-5
X       MAX=50
XC
XC      *   FOR "COMMAND" READING, GO TO 290                                 *
XC
X       IF(ICOMM.NE.1) THEN
XC
X          PRINT *
X          PRINT *,'A FEW MORE INPUTS FOR SCHMITT`S BINNED REGRESSION'
XC
XC      *     NUMBER OF BINS ON AXES                                         *
XC   
X   10     PRINT *
X          WRITE(6,20)
X   20     FORMAT('HOW MANY BINS DO YOU WANT FOR THE X AXIS ? ')
X          CALL DATA1(MX)
X          IF(MX.LE.0) GOTO 10
X   30     PRINT *
X          WRITE(6,40)
X   40     FORMAT('HOW MANY BINS DO YOU WANT FOR THE Y AXIS ? ')
X          CALL DATA1(MY)
X          IF(MY.LE.0) GOTO 30
X          PRINT *
XC
XC      *          TOLERANCE LEVEL                                          *
XC
X   50     PRINT *
X          PRINT *,'DO YOU WANT TO SET'
X          WRITE(6,55)  
X   55  FORMAT('   THE TOLERANCE LEVEL (DEFAULT 1.0E-5) (Y/N) ? ')
X          READ(5,85) CHECK
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X             PRINT *
X             WRITE(6,65)
X   65        FORMAT('WHAT IS THE TOLERANCE LEVEL (E9.3 FORMAT) ? ')
X             READ(5,66) TOL
X   66        FORMAT(E9.3)
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X          GOTO 50
X       ENDIF
XC
XC
XC      *            ITERATION LIMITS                                       *
XC
X   67     PRINT *
X          PRINT *,'DO YOU WANT TO SET THE'
X          WRITE(6,68)
X   68     FORMAT('    ITERATION LIMIT (DEFAULT 50) (Y/N) ? ')
X          READ(5,85) CHECK
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X             PRINT *
X   69        WRITE(6,70)
X   70        FORMAT('WHAT IS THE MAXIMUM ITERATION ? ')
X             CALL DATA1(MAX)
X             IF(MAX.LE.0) GOTO 69 
X          ELSEIF(CHECK.NE.'N'.AND.CHECK.NE.'n') THEN
X             GOTO 67
X          ENDIF
XC
X   71     PRINT *
X          WRITE(6,80)
X   80  FORMAT('DO YOU WANT TO SET BIN SIZES AND ORIGIN  (Y/N) ? ')
X          READ(5,85) CHECK
X   85     FORMAT(A1)
X
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
XC
XC      *      INFORMATION ABOUT BIN SIZES                                   *
XC
X             ISKIP=10
X  110        PRINT *
X             WRITE(6,120)
X  120        FORMAT('WHAT IS THE BIN SIZE FOR THE X AXIS ? ')
X             READ(5,130) XBIN
X  130        FORMAT(F10.3)
X             IF(XBIN.LE.0) GOTO 110
X
X  140        PRINT *
X             WRITE(6,150)
X  150        FORMAT('WHAT IS THE BIN SIZE FOR THE Y AXIS ? ')
X             READ(5,130) YBIN
X             IF(YBIN.LE.0) GOTO 140
XC
XC      *     INFORMATION ABOUT THE ORIGIN                                   *
XC
X  180        PRINT *
X             WRITE(6,190)
X  190        FORMAT('WHERE IS THE ORIGIN OF THE X AXIS ? ')
X             READ(5,130) XORG
X             PRINT *
X             WRITE(6,200)
X  200        FORMAT('WHERE IS THE ORIGIN OF THE Y AXIS ? ')
X             READ(5,130) YORG
X
X          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
X   90        ISKIP=0
X             XBIN=0
X             YBIN=0
X             XORG=0
X             YORG=0
X
X          ELSE
X             GOTO 67
X          ENDIF
XC
XC
XC      *     INFORMATION ABOUT PRINTOUTS                                    *
XC
X  220     PRINT *
X          PRINT *,'DO YOU WANT TO PRINT OUT THE FINAL '
X          WRITE(6,230)
X  230     FORMAT('      2-DIMENSIONAL KM ESTIMATOR (Y/N) ? ')
X          READ(5,85) CHECK
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X             IPRINT=10
X          ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
X             IPRINT=0
X          ELSE
X             GOTO 220
X          ENDIF
X
XC
XC      *      INFORMATION ABOUT ERROR COMPUTATIONS                          *
XC
XC
X          NLAST = 0
X          PRINT *
X  245     PRINT*,'DO YOU WANT TO COMPUTE THE ERRORS'
X          WRITE(6,250)
X  250     FORMAT('       FOR THE REGRESSION COEFFICIENT (Y/N) ? ')
X          READ(5,85) CHECK
X          IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X             PRINT *
X 251         WRITE(6,252)
X 252         FORMAT('DO YOU WISH TO SET THE NUMBER OF BOOTSTRAP ',
X     +              'ITERATIONS? (DEFAULT=200)')
X             READ(5,85) CHECK1
X             IF(CHECK1 .EQ. 'Y' .OR. CHECK1 .EQ. 'y') THEN
X  253           WRITE(6, 254)
X  254           FORMAT('HOW MANY BOOTSTRAP ITERATIONS DO YOU WANT ? ')
X                READ(5,255) NLAST
X  255           FORMAT(I4)
X                IF((NLAST.LT.0) .OR. (NLAST.GT.1000)) THEN
X                  WRITE(6,256) NLAST 
X  256             FORMAT('NUMBER OF ITERATIONS :', I5)
X                  PRINT *
X             PRINT*,'           :SUGGESTED RANGE BETWEEN 100 AND 1000'
X                  GOTO 251
X                ENDIF
X              ELSE IF(CHECK1 .EQ. 'N' .OR. CHECK1 .EQ. 'n') THEN
X                  NLAST = 200
X              ELSE
X                  GOTO 251
X              ENDIF
X              PRINT *
X  257         WRITE(6, 258)
X  258         FORMAT('SEED FOR THE RANDOM NUMBER (NEGATIVE INT) :')
X              READ(5,267) IRAND
X  267         FORMAT(I5)
X              IF(IRAND.GT. 0) GOTO 257
X          ELSEIF(CHECK .NE. 'N'.AND. CHECK .NE. 'n') THEN
X              GOTO 245
X          ENDIF
XC
XC
XC      *      "COMMAND" FILE READING                                        *
XC
XC      *       INFORMATION ABOUT BIN NUMBERS                                *
XC
X       ELSE
X 290      READ(50,300) ((CHAR(I,J),I=1,4),J=1,2)
X 300      FORMAT(8A1)
X          CALL DATA2(CHAR,1,2,MX,LIND)
X          IF((LIND.NE.0) .OR. (MX.LE.0)) THEN
X  301        PRINT *
X             PRINT *, ' NUMBER OF X-AXIS BINS IS WRONG.'
X             STOP
X          ENDIF
X
X          CALL DATA2(CHAR,2,2,MY,LIND)
X          IF((LIND.NE.0) .OR. (MY .LE. 0)) THEN
X  305        PRINT *
X             PRINT *, ' NUMBER OF Y-AXIS BINS IS WRONG.'
X             STOP
X          ENDIF
XC
XC      *      READ ISKIP                                                    *
XC
X          READ(50,300) (CHAR(I,1),I=1,4)
X          CALL DATA2(CHAR,1,1,ISKIP,LIND)
X          IF((LIND.NE.0) .OR. (ISKIP .LT. 0)) THEN
X  308        PRINT *
X             PRINT *, ' BIN SIZE/ORIGIN INDICATOR IS WRONG '
X             STOP
X          ENDIF
XC
XC      *           READ TOLERANCE LEVEL                                     *
XC
X  316     READ(50,317) TOL
X  317     FORMAT(E9.3)
XC
XC                  READ ITERATION LIMIT                                     *
XC
X          READ(50,300) (CHAR(I,1),I=1,4)
X          CALL DATA2(CHAR,1,1,MAX,LIND)
X          IF((LIND.NE.0) .OR. (MAX.LE.0).OR.(MAX.GT.1000)) THEN
X  319        PRINT *
X             PRINT *, ' ITERATION LIMIT IS NOT BETWEEN 1 AND 1000. '
X          ENDIF
XC
XC      *           READ XBIN, YBIN, XORG, YORG                              *
XC
X          READ(50,322) XBIN,YBIN
X          READ(50,322) XORG,YORG
X  322     FORMAT(2F10.3)
XC
XC      *              PRINTING INFORMATION                                  *
XC
X          READ(50,300) (CHAR(I,1),I=1,4)
X          CALL DATA2(CHAR,1,1,IPRINT,LIND)
X          IF((LIND.NE.0) .OR. (IPRINT.LT.0)) THEN
X  323        PRINT *
X             PRINT *, '  2-DIM KAPLAN-MEIER PRINT INDICATOR IS WRONG. '
X             STOP
X          ENDIF
XC
XC      *              ERROR ANALYSIS INFORMATION                            *
XC
X  340     NLAST = 0
X          READ(50,300) (CHAR(I,1),I=1,4)
X          CALL DATA2(CHAR,1,1,NLAST,LIND)
X          IF((LIND .NE. 0) .OR. (NLAST.LT.0)) THEN
X             PRINT *,'  NUMBER OF ITERATIONS FOR THE ERROR COMPUTATION'
X             PRINT *,'     IS NEGATIVE.'
X             STOP
X          ENDIF
X          IF(NLAST.GT.0) THEN
XC
XC      *  SET A SEED FOR THE RANDOM NUMBER GENERATOR; MUST BE A NEGATIVE    *
XC      *  INTEGER.                                                          *
XC
X             READ(50,345) IRAND
X  345        FORMAT(I5)
X             IF(IRAND .GT. 0) THEN
X                IRAND = -IRAND
X                PRINT *,'THE SEED FOR THE RANDOM NUMBER GENERATOR WAS'
X                PRINT *,'CHANGED TO A NEGATIVE VALUE.'
X             ENDIF
X          ENDIF
X
X       ENDIF
X       RETURN
X       END
X
X
XC**************************************************************************
XC*************************** FUNCTION RAN1 ********************************
XC**************************************************************************
X
X      FUNCTION RAN1(IDUM)
XC
XC     *     THIS FUNCTION GIVES A UNIFORM RANDOM NUMBER BETWEEN 0 AND 1.   *
XC     *    YOU NEED TO PROVIDE A SEED TO GENERATE THE FIRST VALUE          *
XC     *    AND THE SEED MUST BE A NEGATIVE INTEGER.                        *
XC     *    REF. NUMERICAL RECIPES P. 196.                                  *
XC
X      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X      DIMENSION R(97)
X      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
X      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
X      PARAMETER (M3=243000,IA3=4561,IC3=51349)
X      DATA IFF /0/
X      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
X        IFF=1
X        IX1=MOD(IC1-IDUM,M1)
X        IX1=MOD(IA1*IX1+IC1,M1)
X        IX2=MOD(IX1,M2)
X        IX1=MOD(IA1*IX1+IC1,M1)
X        IX3=MOD(IX1,M3)
X        DO 11 J=1,97
X          IX1=MOD(IA1*IX1+IC1,M1)
X          IX2=MOD(IA2*IX2+IC2,M2)
X          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
X11      CONTINUE
X        IDUM=1
X      ENDIF
X      IX1=MOD(IA1*IX1+IC1,M1)
X      IX2=MOD(IA2*IX2+IC2,M2)
X      IX3=MOD(IA3*IX3+IC3,M3)
X      J=1+(97*IX3)/M3
X      IF(J.GT.97.OR.J.LT.1)PAUSE
X      RAN1=R(J)
X      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
X      RETURN
X      END
X
XC*************************************************************************
XC********************** SUBROUTINE REARRN  *******************************
XC*************************************************************************
X
X       SUBROUTINE REARRN(IND, X, INX, J, NTOT, MVAR)
X
XC
XC     *    THIS SUBROUTINE REARRANGES THE TIED DATA POINTS SO THAT THE  *
XC     *    RIGHT-CENSORED DATA POINTS COME AFTER THE UNCENSORED VALUE.  *
XC     *    E.G.  1, 2, 3, 3(LOWER), 3, 5, 6  ARE REARRANGED TO          *
XC     *          1, 2, 3, 3, 3(LOWER), 5, 6.                            *
XC     *                                                                 *
XC     *  INPUT                                                          *
XC     *     X(J, I)   : VARIABLES J = 1, 2 (MUST BE SORTED)             *
XC     *     IND(J, I) : CENSORING INDICATOR                             *
XC     *     INX(J, I) : POSITION INDICATOR                              *
XC     *                                                                 *
XC     *  OUTPUT                                                         *
XC     *     REARRANGED X, IND, AND INX                                  *
XC
XC
X
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X       DIMENSION IND(MVAR, NTOT), X(MVAR, NTOT), INX(MVAR, NTOT)
X       
X       I = 0
XC
XC     *            CHECK WHETHER THE DATA POINTS ARE TIED              *
XC
X   10  K = 0
X   20  K = K + 1
X       IF(X(J, I+K) . EQ. X(J, I)) GOTO 20
X       IF(K .GE. 2) THEN
X          K = K - 1
X          KJ = I
XC
XC     *    IF TIED DATA POINTS WERE FOUND, CHECK WHETHER THEY ARE      *
XC     *    DETECTIONS.                                                 *
XC
X   30     IF(IND(J, KJ) .NE. 0) THEN
X             TX    = X(J, KJ)
X             INDTX = IND(J, KJ)
X             INXT  = INX(J, KJ)
XC
XC     *          THE CASE FOR THE UPPER LIMITS                        *
XC
X             IF(IND(J, KJ) .LT. 0) THEN
X                IF(KJ .NE. I) THEN
X                   K1 = KJ - I
X                   DO 50 IL = 1, K1
X                        L1 = KJ - IL +1
X                        X(J, L1)   = X(J, L1-1)
X                        IND(J, L1) = IND(J, L1-1)
X                        INX(J, L1) = INX(J, L1-1)
X   50               CONTINUE
X                    X(J, I)   = TX
X                    IND(J, I) = INDTX
X                    INX(J, I) = INXT
X                 ENDIF
X              ELSE
XC
XC     *      THE CASE FOR THE LOWER LIMITS                              *
XC
X                 IF(KJ .NE. I+K) THEN
X                 ICOUNT = 1
X                 K1 = I + K - KJ
X   55            DO 60 IL = 1, K1
X                    L1 = KJ + IL - 1
X                    X(J, L1)   = X(J, L1+1)
X                    IND(J, L1) = IND(J, L1+1)
X                    INX(J, L1) = INX(J, L1+1)
X   60            CONTINUE
X                 X(J, I+K)   = TX
X                 IND(J, I+K) = INDTX
X                 INX(J, I+K) = INXT
X              
XC
XC     *   CHECK THAT THE VALUE JUST REPLACED IS NOT A LOWER LIMIT. IF  *
XC     *   IT IS, REPEAT THE PROCEDURE.                                 *
XC
X                 ICOUNT = ICOUNT + 1
X                 IF(ICOUNT .LE. K) THEN
X                    IF(IND(J, KJ) .NE. 0) GOTO 55
X                 ENDIF
X              ENDIF
X           ENDIF
X         ENDIF
XC
XC     *   REPEAT UNTIL ALL DATA POINTS ARE TESTED                      *
XC
X         KJ =KJ +1
X         IF(KJ .LE. I+K) GOTO 30
X         ENDIF
X         I = I + K
X         IF(I .LT. NTOT) GOTO 10
X      
X         RETURN
X         END
X
X
XC
XC      **********************************************************************
XC      ****************** SUBROUTINE REGRES  ********************************
XC      **********************************************************************
XC
X       SUBROUTINE REGRES(X,Y,NTOT,NVAR,ALPHA,RMUL,SIGM,R,
X     +                   YFIT,IK,JK,A,XMEAN,SIGMAX,ARRAY,SIGMAA,MVAR) 
XC
XC
XC      *     MULTIVARIABLE LINEAR REGRESSION FIT                            *
XC      *                                                                    *
XC      *     THIS SUBPROGRAM CALCULATES LINEAR REGRESSION COEFFS.,          *
XC      *     VARIANCE,AND CORRELATION COEFFS., BASED ON PHILIP R. BEVINGTON,*
XC      *     "DATA REDUCTION AND ERROR ANALYSIS FOR THE PHYSICAL SCIENCES", *
XC      *     McGRAW HILL (NY:NY), PROGRAM 9-1 OF P. 172                     *
XC      *                                                                    *
XC      *     PARAMETERS                                                     *
XC      *      INPUT                                                         *
XC      *       X       : ARRAY OF DATA POINTS FOR INDEP. VARIABLES          *
XC      *       Y       : ARRAY OF DATA POINTS DEPENDENT VARIABLES           *
XC      *       NTOT    : NUMBER OF PAIRS OF DATA POINTS                     *
XC      *       NVAR    : NUMBER OF COEFFICIENTS                             *
XC      *      OUTPUT                                                        *
XC      *       A       : ARRAY OF COEFFICIENTS                              *
XC      *       ALPHA   : OUTPUT FORM OF A ARRAY (=A)                        *
XC      *       SIGMA   : ARRAY OF STANDARD DEVIATIONS OF COEFFS.            *
XC      *       R       : ARRAY OF LINEAR CORRELAION COEFF.                  *
XC      *       RMUL    : MULTIPLE LINEAR CORRELATION COEFF.                 *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *                 MATINV                                             *
XC
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION X(MVAR,NTOT),Y(NTOT),YFIT(NTOT),IK(NTOT),JK(NTOT)
X       DIMENSION A(MVAR),XMEAN(MVAR),SIGMAX(MVAR),ARRAY(MVAR,MVAR)
X       DIMENSION R(MVAR),SIGMAA(MVAR),SIGM(MVAR),ALPHA(MVAR)
XC
XC      *           INITIALIZE SUMS, ARRAYS, AND OTHERS                      *
XC
X       SUM   =0.0
X       RMUL  =0.0
X       YMEAN =0.0
X       SIGMA =0.0
X       SIGMAO=0.0
X       VARNCE=0.0
X       FNPTS=REAL(NTOT)
X       FREE1=REAL(NTOT-1)
X       FREEN=REAL(NTOT-NVAR) 
X       FREEJ=REAL(NVAR) 
X
X       DO 17 I=1,NTOT
X          YFIT(I)=0.0
X   17  CONTINUE
X
X   21  DO 28 J=1,NVAR 
X          XMEAN(J)=0.0
X          SIGMAX(J)=0.0
X          R(J)=0.0
X          A(J)=0.0
X          SIGMAA(J)=0.0
X
X          DO 27 K=1,NVAR 
X             ARRAY(J,K)=0.0
X   27     CONTINUE
X   28  CONTINUE
XC
XC      *             TAKE MEANS                                             *
XC
X       DO 50 I=1,NTOT
X          YMEAN=YMEAN+Y(I)
X
X          DO 44 J=1,NVAR 
X             XMEAN(J)=XMEAN(J)+X(J,I)
X   44     CONTINUE
X   50  CONTINUE
X
X       YMEAN=YMEAN/FNPTS
X
X       DO 53 J=1,NVAR 
X          XMEAN(J)=XMEAN(J)/FNPTS
X   53  CONTINUE
XC
XC      *           ACCUMULATE MATRICES R AND ARRAY                          *
XC
X       DO 67 I=1,NTOT
X          SIGMA=SIGMA+(Y(I)-YMEAN)**2
X
X          DO 66 J=1,NVAR 
X             SIGMAX(J)=SIGMAX(J)+(X(J,I)-XMEAN(J))**2
X             R(J)=R(J)+(X(J,I)-XMEAN(J))*(Y(I)-YMEAN)
X
X             DO 65 K=1,J
X             ARRAY(J,K)=ARRAY(J,K)+(X(J,I)-XMEAN(J))*(X(K,I)-XMEAN(K))
X   65        CONTINUE
X   66     CONTINUE
X   67  CONTINUE
X
X       SIGMA=DSQRT(SIGMA/FREE1)
X
X       DO 78 J=1,NVAR 
X          SIGMAX(J)=DSQRT(SIGMAX(J)/FREE1)
X          R(J)=R(J)/(FREE1*SIGMAX(J)*SIGMA)
X          DO 77 K=1,J
X             ARRAY(J,K)=ARRAY(J,K)/(FREE1*SIGMAX(J)*SIGMAX(K))
X             ARRAY(K,J)=ARRAY(J,K)
X   77     CONTINUE
X   78  CONTINUE
XC
XC      *            INVERT SYMMETRIC MATRIX                                 *
XC
X       CALL MATINV(ARRAY,NVAR,DET,IK,JK,MVAR)
X
X       IF(DET.EQ.0.0) THEN
X          A1=0.0
X          PRINT *,'******WARNING : DETERMINANT IN MATINV IS 0.*****'
X          RETURN
X       ENDIF
XC
XC      *             CALCULATE COEFFICIENTS                                 *
XC
X  101  A1=YMEAN
X       DO 108 J=1,NVAR 
X
X          DO 104 K=1,NVAR 
X             A(J)=A(J)+R(K)*ARRAY(J,K)
X  104     CONTINUE
X
X          A(J)=A(J)*SIGMA/SIGMAX(J)
X          A1=A1-A(J)*XMEAN(J)
X
X          DO 107 I=1,NTOT
X             YFIT(I)=YFIT(I)+A(J)*X(J,I)
X  107     CONTINUE
X  108  CONTINUE
XC
XC      *            CALCULATE UNCERTAINTIES                                 *
XC
X       DO 113 I=1,NTOT
X          YFIT(I)=YFIT(I)+A1
X          VARNCE=VARNCE+(Y(I)-YFIT(I))**2
X  113  CONTINUE
X
X       VARNCE=VARNCE/FREEN
X
X       DO 133 J=1,NVAR 
X          SIGMAA(J)=ARRAY(J,J)*VARNCE/(FREE1*SIGMAX(J)**2)
X          SIGMAA(J)=DSQRT(SIGMAA(J))
X          RMUL=RMUL+A(J)*R(J)*SIGMAX(J)/SIGMA
X  133  CONTINUE
X
X       RMUL=DSQRT(RMUL)
X       SIGMAO=VARNCE/FNPTS
X
X       DO 145 J=1,NVAR 
X          DO 144 K=1,NVAR 
X             SIGMAO=SIGMAO+VARNCE*XMEAN(J)*XMEAN(K)*ARRAY(J,K)
X     +                   /(FREE1*SIGMAX(J)*SIGMAX(K))
X  144     CONTINUE
X  145  CONTINUE
X
X       DO 147 J=1,NVAR 
X          JJ=J+1
X          ALPHA(JJ)=A(J)
X          SIGM(JJ)=SIGMAA(J)
X  147  CONTINUE
X
X       ALPHA(1)=A1
X       SIGM(1)=DSQRT(SIGMAO)
X       A(NVAR +2)=DSQRT(VARNCE)
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE RMILLS  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE RMILLS(X,FUNC,TOL)
XC
XC      *      ALGORITHM AS 138.1 APPL.STATST. (1979) VOL.28. NO.2           *
XC      *                                                                    *
XC      *         COMPUTE THE RECIPROCAL OF MILLS RATIO                      *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DATA FPI /1.2533141/, FPII /0.7978846/
XC
X       FUNC=0.0
X       IF(X .LT. -10.0) RETURN
X       FUNC=FPII
X       Y=DABS(X)
X       IF(Y .LT. 0.000001) RETURN
X       SGN=1.0
X       IF(X.LT.0.0) SGN=-1.0
X
X       IF(Y.LE.2.0) THEN
X          S=0.0
X          A=1.0
X          T=Y
X          R=Y
X          B=Y**2
X   40     A=A+2.0
X          S=T
X          R=R*B/A
X          T=T+R
X          IF(R.GT.TOL) GOTO 40
X          FUNC=1.0/(FPI*DEXP(0.5*B)-SGN*T)
X          RETURN
X       ENDIF
X
X  100  A=2.0
X       B1=Y
X       S=Y
X       A1=Y**2+1.0
X       A2=Y*(A1+2.0)
X       B2=A1+1.0
X       T=A2/B2
X
X  140  A=A+1.0
X       A0=A1
X       A1=A2
X       A2=Y*A1+A*A0
X       B0=B1
X       B1=B2
X       B2=Y*B1+A*B0
X       R=S
X       S=T
X       T=A2/B2
X
X       IF((T-R.GT.TOL) .OR.(T-S.GT.TOL)) GOTO 140
X       FUNC=T
X
X       IF(SGN.LT.0.0) FUNC=T/(2.0*FPI*DEXP(0.5*Y**2)*T-1.0)
X
X       RETURN
X       END
X
XC***************************************************************************
XC************************ SUBROUTINE SCHMIT  *******************************
XC***************************************************************************
XC
XC
X       SUBROUTINE SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
X     +                   YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
X     +                   ALPHA,BETA,MM,M1,M2,M3,M4,M5,M6,M7,M8,
X     +                   A,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
X     +                   IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
X     +                   IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
XC
XC
XC      *                                                                    *
XC      *          THIS SUBROUTINE COMPUTES LINEAR REGRESSION COEFFICIENTS   *
XC      *      ALPHA AND BETA (INTERCEPT AND SLOPE) BY SCHMITT'S BINNED      *
XC      *      METHOD. BECAUSE OF THE BINNING METHOD, FINER BINNING GIVES    *
XC      *      BETTER RESULTS, THOUGH THE CPU TIME MAY INCREASE VERY         *
XC      *      MUCH. ALSO IF THE BINS ARE TOO FINE, THE CENSORED POINTS      *
XC      *      MAY NOT FALL ON THE DETECTIONS.                               *
XC      *                                                                    *
XC      *                                                                    *
XC      *      INPUT                                                         *
XC      *            X(I): INDEPENDENT VARIABLE                              *
XC      *            Y(I): DEPENDENT VARIABLE                                *
XC      *           NP(I): INDICATOR OF CENSORED STATUS                      *
XC      *                  IF NP(I)=0  : DETECTION                           *
XC      *                          =1  : Y(I) IS LOWER LIMIT                 *
XC      *                          =2  : X(I) IS LOWER LIMIT                 *
XC      *                          =3  : BOTH ARE LOWER LIMITS               *
XC      *                          =4  : Y(I) IS LOWER AND X(I) IS UPPER     *
XC      *                 FOR THE UPPER LIMITS, CHANGE THE SIGNS             *
XC      *           NTOT : NUMBER OF DATA POINTS                             *
XC      *          MPLONE: NUMBER OF PARAMETERS TO BE ESTIMATED              *
XC      *                  (IN THIS PROGRAM, MPLONE IS ALWAYS 3)             *
XC      *          MAXITS: NUMBER OF MAXIMUM ITERATION                       *
XC      *           ALH  : DUMMY                                             *
XC      *        DELTA(I): DELTA(2) CONTAINS TOLERANCE FOR THE               * 
XC      *                  COMPUTATION OF F'S.                               *
XC      *            MX  : NUMBER OF BINS IN THE INDEPENDENT VARIABLE        *
XC      *            MY  : NUMBER OF BINS IN THE DEPENDENT VARIABLE          *
XC      *          ISKIP : INDICATOR FOR THE BINNING. IF ISKIP=0, THE        *
XC      *                  SUBROUTINE BIN WILL COMPUTE THE INFORMATION       *
XC      *                  ABOUT THE BINNING INFORMATION                     *
XC      *                  IF ISKIP>0, THE BINNING INFORMATION (BIN SIZES    *
XC      *                  ORIGIN) MUST BE PROVIDED.                         *
XC      *         ICENS  : IF THE DATA SET IS KNOWN TO :                     *
XC      *                    CONTAIN LOWER LIMITS ONLY,   ICENS>0            *
XC      *                    CONTAIN UPPER LIMITS ONLY,   ICENS<0            *
XC      *                    CONTAIN MIXED LIMITS,        ICENS=0            *
XC      *                     OR NOT KNOWN                ICENS=0            *
XC      *         IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :             *
XC      *           XBIN : THE BIN SIZE FOR THE X AXIS                       *
XC      *           YBIN : THE BIN SIZE FOR THE Y AXIS                       *
XC      *           XORG : THE ORIGIN OF X                                   *
XC      *           YORG : THE ORIGIN OF Y                                   *
XC      *       NN1=NC1  : NUMBER OF Y LOWER LIMITS                          *
XC      *       NN2=NC2  : NUMBER OF X LOWER LIMITS                          *
XC      *       NN3=NC3  : NUMBER OF DOUBLE LOWER LIMITS                     *
XC      *       NN4=NC4  : NUMBER OF Y LOWER, X UPPER LIMITS                 *
XC      *       NN5=NC5  : NUMBER OF Y UPPER LIMITS                          *
XC      *       NN6=NC6  : NUMBER OF X UPPER LIMITS                          *
XC      *       NN7=NC7  : NUMBER OF DOUBLE UPPER LIMITS                     *
XC      *       NN8=NC8  : NUMBER OF Y UPPER, XLOWER LIMITS                  *
XC      *         TOL    : TOLERANCE LEVEL                                   *
XC      *         MAX    : MAXIMUM NUMBER OF ITERATIONS                      *
XC      *                                                                    *
XC      *     WORK                                                           *
XC      *        F(I,J)  : NUMBER OF DATA POINTS IN THE BIN(I,J)             *
XC      *                   (WEIGHTED BY CENSORED POINTS)                    *
XC      *        A(I,J)  : MATRIX WHICH CONTAINS INFORMATION OF BIN          *
XC      *                  POSITION  I=1,4 AND J=1,MX*MY                     *
XC      *       TH(I)    : VECTOR " A(I,J)*F(I,J), I=1,4                     *
XC      *        SUM     : TOTAL NUMBER OF POINTS.  THIS IS APPROXIMATELY    *
XC      *                  EQUAL TO NTOT.  THE DIFFERENCE BETWEEN THE TWO    *
XC      *                  VALUES DEPENDS ON THE TOLERANCE LEVEL.            *
XC      *                                                                    *
XC      *     OUTPUT                                                         *
XC      *          ALPHA : INTERCEPT COEFFICIENT                             *
XC      *           BETA : SLOPE COEFFICIENT                                 *
XC      *                                                                    *
XC      *     SUBROUTINES:                                                   *
XC      *         GRDPRB : THE SUBROUTINE WHICH COMPUTES THE TWO-DIMENSIONAL * 
XC      *                  KAPLAN-MEIER PROBABILITY OF THE BINS              *
XC      *                                                                    *
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION A(5,IB),TH(5) 
X       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB),F(IB,IB)
X
X       DIMENSION IBWRK1(IB,IB),IBWRK2(IB,IB),IBWRK3(IB,IB)
X       DIMENSION IBWRK4(IB,IB),IBWRK5(IB,IB),IBWRK6(IB,IB)
X       DIMENSION IBWRK7(IB,IB),IBWRK8(IB,IB),IBWRK9(IB,IB)
X
X       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
X       DIMENSION BWRK1(IB,IB),SWRK1(MVAR),DWRK1(MVAR,NTOT)
X       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
XC
XC      *                CALL SUBROUTINE GRDPRB                              *
XC
X       CALL GRDPRB(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,YBIN,
X     +             XORG,YORG,TOL,MAX,MM,M1,M2,M3,M4,M5,M6,M7,M8,
X     +             X,Y,NP,XB,YB,F,BWRK1,IBWRK1,IBWRK2,IBWRK3,
X     +             IBWRK4,IBWRK5,IBWRK6,IBWRK7,IBWRK8,IBWRK9,
X     +             IWRK1,IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
X
XC
XC
XC      *               MAKE MATRIX A(I,J)                                   *
XC
X       DO 1120 J=1,MY
X          DO 1110 I=1,MX
X             IJ=I+MX*(J-1)
X             A(1,IJ)=XB(I)
X             A(2,IJ)=XB(I)**2
X             A(3,IJ)=XB(I)*YB(J)
X             A(4,IJ)=YB(J)
X             A(5,IJ)=YB(J)**2
X 1110     CONTINUE
X 1120  CONTINUE
XC
XC      *             COMPUTE THE VECTOR THETA : TH(I)                       *
XC
X       DO 1400 I=1,5
X          TH(I)=0.0
X 1400  CONTINUE
X
X       DO 1430 J=1,MY
X          DO 1420 I=1,MX
X          IJ = I + MX*(J-1)
X             DO 1410 K=1,5
X                TH(K)=TH(K)+A(K,IJ)*F(I,J)*NTOT
X 1410        CONTINUE
X 1420     CONTINUE
X 1430  CONTINUE
X
X       SUM = 0.0
X       DO 1460 I = 1, MX
X          DO 1450 J = 1, MY
X             SUM = SUM + F(I,J)*NTOT
X 1450     CONTINUE
X 1460  CONTINUE
X 
XC
XC      *     COMPUTE THE SLOPE COEFFICIENT : BETA, AND THE INTERCEPT        *
XC      *     COEFFICIENT : ALPHA.                                           *
XC
X       DEN=SUM*TH(2)-TH(1)**2
X       BETA=(SUM*TH(3)-TH(1)*TH(4))/DEN
X       ALPHA=(TH(2)*TH(4)-TH(1)*TH(3))/DEN
XC
XC
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE SORT1  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE SORT1(ID,X,Y,NTOT,NVAR,INDEX,X1,MVAR)
XC
XC      *       BUBBLE SORT PROGRAM                                          *
XC      *       THIS PROGRAM ARRANGES OBSERVATIONS IN ASCENDING ORDER        *
XC      *       ALSO IF THERE ARE TIED DATA POINTS, IT CHECKS THE CENSORING  *
XC      *       STATUS AND ORDERS THEM SO THAT A DETECTED POINT COMES        *
XC      *       FIRST.                                                       *
XC      *                                                                    *
XC      *       INPUT : INDEX(I): POSITION INDICATOR                         *
XC      *                 ID(I) : INDICATOR OF CENSORING                     *
XC      *                 X(J,I): INDEPENDENT VARIABLE; NVAR DIMENSION       *
XC      *                 Y(I)  : DEPENDENT VARIABLE                         *
XC      *                 NTOT  : NUMBER OF DATA POINTS                      *
XC      *                 NVAR  : NUMBER OF INDEPENDENT VARIABLE             *
XC      *                                                                    *
XC      *      OUTPUT :   ID, X, AND Y IN ASCENDING ORDER WITH INDEX         *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION ID(NTOT),X(MVAR,NTOT),Y(NTOT),X1(MVAR),INDEX(NTOT)
XC
XC      *        SORTING IN Y, ASCENDING ORDER                               *
XC
X       DO 10 I=1,NTOT
X          INDEX(I)=I
X   10  CONTINUE
XC
X       IF(NTOT.EQ.1) GOTO 100
XC
X       DO 90 I=1,NTOT
X          J=NTOT-I+1
X          JJ=J-1
X          IF(JJ.GE.1) THEN   
XC
X             DO 80 K=1,JJ
X                IF(Y(K).GT.Y(J)) THEN 
X
X                   ID1=ID(J)
X                   DO 50 L=1,NVAR
X                      X1(L)=X(L,J)
X  50               CONTINUE
X                   Y1=Y(J)
X                   INDX=INDEX(J)
X
X                   ID(J)=ID(K)
X                   DO 60 L=1,NVAR
X                      X(L,J)=X(L,K)
X  60               CONTINUE
X                   Y(J)=Y(K)
X                   INDEX(J)=INDEX(K)
X
X                   ID(K)=ID1
X                   DO 70 L=1,NVAR
X                      X(L,K)=X1(L)
X  70               CONTINUE
X                   Y(K)=Y1
X                   INDEX(K)=INDX
X                ENDIF
X  80         CONTINUE
X          ENDIF
X  90   CONTINUE
XC
X 100   RETURN
X       END
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE SORT2 ********************************
XC      **********************************************************************
XC
X       SUBROUTINE SORT2(XY, ID1, ID2, NTOT)
XC
XC      *                BUBBLE SORT PROGRAM.                                *
XC      *       OBTAINED FROM ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL*
XC      *       DATA ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS         *
XC      *       (BELMONT:CA)                                                 *
XC
XC*******       THE COMMON STATEMENT IS DIFFERENT FROM SMSDA.                *
XC
XC      *       THIS SUBROUTINE WAS MODIFIED ON 4/20/90.                     *
XC      *       DIMENSION DECLARATION.                                       *
XC
X
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
X       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X       DO 1 I = 1, NCOMP
X          J=NCOMP-I+1
X          JJ=J-1
X          IF(JJ .GE. 1) THEN
X             DO 2 K = 1, JJ
X                IF(XY(K) .EQ. XY(J)) THEN
XC
XC      *   IF DATA POINTS ARE TIED, THEN CHECK THE CENSORSHIP OF BOTH DATA  *
XC      *   POINTS. PUT CENSORED DATA POINTS AFTER DETECTED ONES.            *
XC
X                   IF(ID1(J)-ID1(K) .GE. 0) THEN 
X                      GOTO 2
X                   ELSE
X                      X1=XY(J)
X                      ITEM=ID1(J)
X                      ICTE=ID2(J)
X
X                      XY(J)=XY(K)
X                      ID1(J)=ID1(K)
X                      ID2(J)=ID2(K)
X
X                      XY(K)=X1
X                      ID1(K)=ITEM
X                      ID2(K)=ICTE
X                      GOTO 2
X                   ENDIF
X                ENDIF
X
X   3            IF(XY(K) .GT. XY(J)) THEN
X                   X1=XY(J)
X                   ITEM=ID1(J)
X                   ICTE=ID2(J)
X
X                   XY(J)=XY(K)
X                   ID1(J)=ID1(K)
X                   ID2(J)=ID2(K)
X
X                   XY(K)=X1
X                   ID1(K)=ITEM
X                   ID2(K)=ICTE
X                ENDIF
X   2         CONTINUE
X          ENDIF
X   1   CONTINUE
X
X       RETURN
X       END
X
XC**************************************************************************
XC************************* SUBROUTINE SPEARHO *****************************
XC**************************************************************************
X
X       SUBROUTINE SPEARHO(XF, NTOT, RHO, PROB, MVAR)
X
XC
XC     *     THIS SUBROUTINE COMPUTES THE SPEARMAN'S RHO STATISTIC AND   *
XC     *     ITS PROBABILITY UNDER THE NULL HYPOTHESIS.                  *
XC     *                                                                 *
XC     *   INPUT                                                         *
XC     *      XF(J, I) : RANKS OF TWO VARIABLES J = 1, 2                 *
XC     *      NTOT     : NUMBER OF DATA POINTS                           *
XC     *                                                                 *
XC     *   OUTPUT                                                        *
XC     *       RHO     : SPEARMAN'S RHO                                  *
XC     *       PROB    : PROBABILITY                                     *
XC     *                                                                 *
XC     *   SUBROUTINE AGAUSS                                             *
XC
XC
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X       DIMENSION XF(MVAR, NTOT)
X       
X       XSUM = 0.0
X       X2SUM = 0.0
X       YSUM = 0.0
X       Y2SUM = 0.0
X       XYSUM = 0.0
X       RN = REAL(NTOT)
X       
X       DO 100 I = 1, NTOT
X          XSUM = XSUM + XF(1,I)
X          X2SUM = X2SUM + XF(1,I)**2
X          YSUM = YSUM + XF(2,I)
X          Y2SUM = Y2SUM + XF(2,I)**2
X          XYSUM = XYSUM + XF(1,I)*XF(2,I)
X  100  CONTINUE
X  
X       XBAR = XSUM/REAL(NTOT)
X       YBAR = YSUM/REAL(NTOT)
X
X       SXX = X2SUM - REAL(NTOT)*XBAR**2
X       SYY = Y2SUM - REAL(NTOT)*YBAR**2
X       SXY = XYSUM - REAL(NTOT)*XBAR*YBAR
X
X       RHO = SXY/SQRT(SXX*SYY)
XC
XC******  THE NEXT APPROXIMATION IS GOOD ONLY WHEN N IS LARGE (E.G. >30)  *
XC
X       
X       Z   = RHO*DSQRT(RN -1.0)
X       
X       PROB = 1.0 - AGAUSS(Z)
X       
X       RETURN
X       END
X
XC**************************************************************************
XC************************* SUBROUTINE SPRMAN ******************************
XC**************************************************************************
X
X       SUBROUTINE SPRMAN(IND, X, Y, NTOT, OUTPUT, IPRSP, MVAR,
X     +                    TEMP, INDT, INT, IXD, INX, INDX, Z1, XF, RX,
X     +                    XX,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
X     +                    WRK9,WRK10,IWRK1,IWRK2,DWRK1,DWRK2,SWRK1)
XC
XC     *    SPEARMAN'S RHO CORRELATION TEST                               *
XC     *                                                                  *
XC     *    THIS PROGRAM COMPUTES A CORRELATION PROBABILITY BETWEEN TWO   *
XC     *    VARIABLES BY SPEARMAN'S RHO. FOR CENSORED DATA POINTS,        *
XC     *    AKRITAS' RANKING METHOD IS USED.                              *
XC     * REFERENCE                                                        *
XC     *    PENN STATE UNIVERSITY, DEPARTMENT OF STATISTICS, TECHNICAL    *
XC     *    REPORTS AND PREPRINTS SERIES, NUMBER 87, "ALIGNED RANK TESTS  *
XC     *    FOR REGRESSION WITH CENSORED DATA', MICHAEL G. AKRITAS,       *
XC     *    SEPTEMBER 1989.                                               *
XC     *                                                                  *
XC     * INPUT                                                            *
XC     *      X(1,I)   : INDEPENDENT VARIABLE                             *
XC     *      Y(I)     : DEPENDENT VARIABLE                               *
XC     *      IND(I)   : CENSORED INDICATOR                               *
XC     *      NTOT     : NUMBER OF DATA POINTS                            *
XC     *      MVAR     : DIMENSION SIZE                                   *
XC     *                                                                  *
XC     * OUTPUT                                                           *
XC     *        RHO    : SPEARMAN'S RHO                                   *
XC     *        PROB   : PROBABILITY FOR NULL HYPOTHESIS                  *
XC     *                                                                  *
XC     * OTHER                                                            *
XC     *      XX(J, I) : TWO VARIABLES J = 1, 2                           *
XC     *      RX(J, I) : RANK OF J-TH VARIABLE                            *
XC     *      INX(J, I): POSITION OF ORIGINAL VALUES                      *
XC     *                                                                  *
XC     * SUBROUTINES                                                      *
XC     *      SORT1, REARRN, AKRANK, SPEARHO                              *
XC     *                                                                  *
XC
XC
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X
X       DIMENSION X(MVAR, NTOT), Y(NTOT), RX(MVAR, NTOT), TEMP(NTOT)
X       DIMENSION XF(MVAR, NTOT), IXD(MVAR, NTOT), INT(NTOT)
X       DIMENSION IND(NTOT), INX(MVAR, NTOT), INDX(MVAR, NTOT)
X       DIMENSION Z1(MVAR, NTOT), INDT(NTOT),XX(MVAR, NTOT)
X       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT)
X       DIMENSION WRK5(NTOT),WRK6(NTOT),WRK7(NTOT),WRK8(NTOT)
X       DIMENSION WRK9(NTOT),WRK10(NTOT),IWRK1(NTOT),IWRK2(NTOT)
X       DIMENSION DWRK1(MVAR,NTOT),DWRK2(MVAR,NTOT),SWRK1(MVAR)
X       CHARACTER*9 OUTPUT
XC
XC
XC     *            INITIALIZATION                                        *
XC
X                    
X       IKEEP = IPRSP
X       DO 40  I = 1, NTOT
X          DO 35 J = 1, 2
X             INX(J, I) = I
X   35     CONTINUE
X          XX(1, I) = X(1, I)
X          XX(2, I) = Y(I)
X   40  CONTINUE
XC
XC     *       ADJUST THE CENSORSHIP INDICATOR TO ONE DIMENSIONAL FORM     *
XC
X       DO 50 I = 1, NTOT
X          INDX(1, I) = 0
X          INDX(2, I) = 0
X          IAD = IABS(IND(I))
X
X          IF(IAD .GE. 2) THEN
X              INDX(1, I) = IND(I)/IAD
X          ENDIF
X          IF(IAD .EQ. 4) INDX(1, I) = -INDX(1,I)
X
X          IF(IAD .EQ. 1) THEN
X             INDX(2, I) = IND(I)/IAD
X          ELSEIF(IAD .GE. 3) THEN
X             INDX(2, I) = IND(I)/IAD
X          ENDIF
X   50  CONTINUE
X          
XC
XC     *    START COMPUTATION OF RANKS OF VARIABLES                      *
XC
X       DO 140 J = 1, 2
X          DO 120 I = 1, NTOT
X            TEMP(I) = XX(J, I)
X            INDT(I) = INDX(J, I)
X            Z1(1, I) = 0.0
X            INT(I) = I
X  120     CONTINUE
X       
XC
XC     *        CALL SORT1 : SORT INTO ASCENDING ORDER                    *
XC
X          CALL SORT1(INDT, Z1, TEMP, NTOT, 1, INT, SWRK1, MVAR)
X              
X          DO 130 I = 1, NTOT
X             XX(J, I) = TEMP(I)
X             INDX(J, I) = INDT(I)
X             INX(J, I) = INT(I)
X  130     CONTINUE
XC
XC     *     REARRANGE TIED DATA POINTS SO THAT THE CENSORED POINTS COME *
XC     *     AFTER THE UNCENSORED VALUE                                  *
XC
X       CALL REARRN(INDX, XX, INX, J, NTOT, MVAR)
X       
XC
XC
XC     *         AKRANK     : AKRITAS' RANKING METHOD                    *
XC
X          CALL AKRANK(INDX, XX, NTOT, J, RX, MVAR,
X     +              WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,IWRK1,WRK7,
X     +              DWRK1,WRK8,WRK9,WRK10,DWRK2,IWRK2,SWRK1)
X       
X  140  CONTINUE
XC
XC
XC     *   REARRANGE THE RANKS TO THE ORIGINAL POSITIONS                 *
XC
X       DO 200 I = 1, NTOT
X         
X          DO 180 J = 1, 2
X         
X             DO 150 K = 1, NTOT
X                IF(INX(J, K) .EQ. I) THEN
X                   IXD(J, I) = K
X                   GOTO 160
X                ENDIF
X  150        CONTINUE
X             
X  160        XF(J, I) = RX(J, IXD(J, I))
X  180    CONTINUE
X  200  CONTINUE
X
XC
XC     *   GET NEW VALUES OF INDX (OLD VALUES ARE OUT OF PROPER ORDER)   *
XC
X       DO 260 I = 1, NTOT
X          INDX(1, I) = 0
X          INDX(2, I) = 0
X          IAD = IABS(IND(I))
X
X          IF(IAD .GE. 2) THEN
X              INDX(1, I) = IND(I)/IAD
X          ENDIF
X          IF(IAD .EQ. 4) INDX(1, I) = -INDX(1,I)
X
X          IF(IAD .EQ. 1) THEN
X             INDX(2, I) = IND(I)/IAD
X          ELSEIF(IAD .GE. 3) THEN
X             INDX(2, I) = IND(I)/IAD
X          ENDIF
X  260  CONTINUE
X
X
XC      *    MAKE SURE THAT RANKS AGREE AT TIED POINTS                     *
X       DO 330 I = 1, NTOT-1
X          IT = 1
X          T1 = XF(1,I)
X          DO 290 J = I+1, NTOT
X             IF(X(1,I).EQ.X(1,J).AND.INDX(1,I).EQ.INDX(1,J)) THEN
X                IT = IT + 1
X                T1 = T1 + XF(1,J)
X             ENDIF
X  290     CONTINUE
X
X          IF(IT .GT. 1) THEN
X             TAVG = T1/REAL(IT)
X             XF(1,I) = TAVG
X             DO 300 J = I+1, NTOT
X                IF(X(1,I).EQ.X(1,J).AND.INDX(1,I).EQ.INDX(1,J)) THEN     
X                   XF(1,J) = TAVG
X                ENDIF
X  300        CONTINUE
X          ENDIF
X
X          IT = 1
X          T1 = XF(2,I)
X          DO 310 J = I+1, NTOT
X             IF(Y(I).EQ.Y(J).AND.INDX(2,I).EQ.INDX(2,J)) THEN
X                IT = IT + 1
X                T1 = T1 + XF(2,J)
X             ENDIF
X  310     CONTINUE
X
X          IF(IT .GT. 1) THEN
X             TAVG = T1/REAL(IT)
X             XF(2,I) = TAVG
X             DO 320 J = I+1, NTOT
X                IF(Y(I).EQ.Y(J).AND.INDX(2,I).EQ.INDX(2,J)) THEN
X                   XF(2,J) = TAVG
X                ENDIF
X  320        CONTINUE
X          ENDIF 
X  330  CONTINUE
X      
XC
XC
XC     *         CALL SPEARHO : SPEARMAN'S RHO CORRELATION TEST          *
XC
XC
X       CALL SPEARHO(XF, NTOT, RHO, PROB, MVAR)
XC
XC
X       IF(OUTPUT .EQ. '        ') THEN
X          WRITE(6, 220)
X          WRITE(6, 210)
X          WRITE(6, 220)
X          WRITE(6, 230) RHO
X          WRITE(6, 240) PROB
X          WRITE(6, 220)
XC
X          IF(IKEEP .EQ. 1) THEN
X             WRITE(6, 250)
X             WRITE(6,220)
X             WRITE(6,225)
X             DO 998 I = 1, NTOT
X                WRITE(6, 997) IND(I),X(1, I),Y(I),(XF(J, I), J = 1, 2)
X  998        CONTINUE
X             WRITE(6,220)
X          ENDIF
X       ELSE
X          WRITE(60, 220)
X          WRITE(60, 210)
X          WRITE(60, 220)
X          WRITE(60, 230) RHO
X          WRITE(60, 240) PROB
X          WRITE(60, 220)
XC
X          IF(IKEEP .EQ. 1) THEN
X             WRITE(60, 250)
X             WRITE(60,220)
X             WRITE(60,225)
X             DO 999 I = 1, NTOT
X                WRITE(60, 997) IND(I),X(1, I),Y(I),(XF(J, I), J = 1, 2)
X  999        CONTINUE
X             WRITE(60,220)
X          ENDIF
X       ENDIF
X
X  210  FORMAT(5X,'CORRELATION TEST BY SPEARMAN`S RHO')
X  220  FORMAT(' ')
X  225  FORMAT(10X,'IND      X         Y        X RANK     Y RANK ')
X  230  FORMAT(7X,'SPEARMANS RHO = ', F10.3)
X  240  FORMAT(7X,'PROBABILITY   = ', F11.4,4X, 
X     +       '(NULL HYPOTHESIS, ACCURATE ONLY IF N > 30)')
X  250  FORMAT(10X,'INPUT DATA AND THEIR RANKS')
X  997  FORMAT(10X,I4, 4F10.3)
X
X       RETURN
X       END
X
X
XC
XC      **********************************************************************
XC      ***************** SUBROUTINE STAT  ***********************************
XC      **********************************************************************
XC
X       SUBROUTINE STAT(SCORE,WSC,XY,ID1,ID2,NTOT)
XC
XC      *          GIVEN THE SCORES OF EACH OBSERVATION, THIS                *
XC      *          SUBROUTINE COMPUTES THE FINAL TEST STATISTIC.             *
XC      *                                                                    *
XC      *   INPUT                                                            *
XC      *         SCORE(I)  : SCORE VECTOR                                   *
XC      *           XY(I)   : DATA                                           *
XC      *           ID2(I)  : INDICATOR OF CENSORING                         *
XC      *           NTOT    : NUMBER OF DATA POINTS                          *
XC      *   OUTPUT                                                           *
XC      *           WSC     : STATISTIC                                      *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION SCORE(NTOT)
X       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X       XN=NCOMP
X       WW=0.0
X       DO 26 I=1,NCOMP
X          IF(ID2(I).NE.2) WW=WW+SCORE(I)
X  26   CONTINUE
XC
X       IF((IFULL.EQ.1).AND.(LO.EQ.0)) WRITE(6,30) WW
X       IF((IFULL.EQ.1).AND.(LO.EQ.1)) WRITE(60,30) WW
X  30   FORMAT(10X,'WW =',F10.2)
X
X       SUM=0.0
X       DO 27 I=1,NCOMP
X          SUM=SUM+SCORE(I)**2
X  27   CONTINUE
X
X       XN1=REAL(N1)
X       XN2=REAL(N2)
X       VAR=SUM*XN1*XN2/(XN*(XN-1.0)) 
XC
X       IF((IFULL.EQ.1).AND.(LO.EQ.0)) WRITE(6,32) VAR
X       IF((IFULL.EQ.1).AND.(LO.EQ.1)) WRITE(60,32) VAR
X  32   FORMAT(10X,'VAR=',F11.3)
X
X       WSC=WW/DSQRT(VAR)
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE SUMRY  ******************************
XC      **********************************************************************
XC
X       SUBROUTINE SUMRY(U,IU,S,NTOT,FINT)
XC
XC      *       THIS SUBROUTINE CALCULATES AND PRINTS 75, 50, AND            *
XC      *       25, PERCENTILES  OF SURVIVAL FOR A SURVIVAL CURVE.           *
XC      *       S AND U ARE ARRAYS CONTAINING POINTS FOR WHICH VALUES OF THE *
XC      *       SURVIVAL CURVE WERE CALCULATED, IU IS THE NUMBER OF          *
XC      *       UNCENSORED POINTS.  ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
XC      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
XC      *       PUBLICATIONS (BELMONT:CA).                                   *
XC      *                                                                    *
XC      *       INPUT       U : UNCENSORED DATA                              *
XC      *                   S : PL ESTIMATOR                                 *
XC      *       WORK       TY : PERCENTILE INDICATOR AT 75, 50, 25           *
XC      *       OUTPUT    FINT: VALUES AT 75, 50, 25 PERCENTILES             *
XC      *                                                                    *
XC      *     THIS SUBROUTINE IS FROM SMSDA, EXCEPT PRINTING COMMAND.        *
XC      *                                                                    *
XC
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION U(IU),S(NTOT),TY(3),FINT(3)
XC
XC
X       TY(1)=0.75
X       TY(2)=0.50
X       TY(3)=0.25
X
XC      *       IF THE NUMBER OF DATA POINTS IS SMALLER THAN 4, NO          *
XC      *       PERCENTILES CAN BE OBTAINED.                                *
XC
X       DO 40 I=1,3
X          FINT(I)=0.0
X  40   CONTINUE
X
X       L=1
X       IF(IU.LE.3) RETURN
X       NN=IU-1
X
X       DO 100 I=1,3
X          IF(S(1).LT.TY(I)) THEN
X             FINT(I) = U(1)-(TY(I)-S(1))/(1-S(1))*(U(1)-0)
X          ELSE
X             DO  90 J=L,NN
X                IF((S(J).GE.TY(I)) .AND. (S(J+1).LE.TY(I))) THEN 
X                   FINT(I)=U(J)-(S(J)-TY(I))/(S(J)-S(J+1))*(U(J)-U(J+1))
X                   L=J+1
X                   GOTO 100
X                ENDIF
X  90         CONTINUE
X          ENDIF
X 100   CONTINUE
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ******************* SUBROUTINE SYMINV  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE SYMINV(A,N,C,W,NULLTY,NA,NC,NW,IFAULT)
XC
XC      *      ALGORITHM AS 7 J.R.STATIST.SOC.C. (1968) VOL.17, NO.2         *
XC      *                                                                    *
XC      *        FORMS IN C( ) A LOWER TRIANGULAR MATRIX, WHICH IS A         *
XC      *        GENERALISED INVERSE OF THE POSITIVE SEMI-DEFINITE SYMMETRIC *
XC      *        MATRIX A( ) OF ORDER N.                                     *
XC      *        C( ) MAY COINCIDE WITH A( ). NULLTY IS RETURNED AS          *
XC      *        THE NULLITY OF A( ). IFAULT IS RETURNED AS 1 IF             *
XC      *        N.LT.1,OTHERWISE ZERO. W( ) IS A WORK ARRAY OF              *
XC      *        LENGTH AT LEAST N THAT IS ALLOCATED BY THE CALLING          *
XC      *        ROUTINE.                                                    *
XC
XC      *        NOTE : THE VARIABLES NA,NC,AND,NW,HAVE BEEN ADDED           *
XC      *               TO THE ARGUMENT LIST AND ARE USED TO                 *
XC      *               DIMENSION TO ARRAYS A,C,AND W, RESPECTIVELY.         *
XC      *               (BY WOLYNETZ (1979))                                 *
XC      *                                                                    *
XC      *        SUBROUTINES                                                 *
XC      *               CHOL                                                 *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION A(NA),C(NC),W(NW)
XC
X       NROW=N
X       IFAULT=1
X       IF(NROW.GT.0) THEN
X          IFAULT=0
X
X          CALL CHOL(A,NROW,C,NULLTY,NA,NC,IFAULT)
X
X          IF(IFAULT.EQ.0) THEN
X             NN=(NROW*(NROW+1))/2
X             IROW=NROW
X             NDIAG=NN
X   16        IF(C(NDIAG).NE.0.0) THEN
X                L=NDIAG
X
X                DO 10 I=IROW,NROW
X                   W(I)=C(L)
X                   L=L+I
X   10           CONTINUE
X
X                ICOL=NROW
X                JCOL=NN
X                MDIAG=NN
X   15           L=JCOL
X                X=0.0
X                IF(ICOL.EQ.IROW) X=1.0/W(IROW)
X
X                K=NROW
X   13           IF(K.NE.IROW) THEN
X                   X=X-W(K)*C(L)
X                   K=K-1
X                   L=L-1
X                   IF(L.GT.MDIAG) L=L-K+1
X                   GOTO 13
X                ENDIF
X
X                C(L)=X/W(IROW)
X                IF(ICOL.EQ.IROW) GOTO 14
X                MDIAG=MDIAG-ICOL
X                ICOL=ICOL-1
X                JCOL=JCOL-1
X                GOTO 15
X             ENDIF
X
X             L=NDIAG
X             DO 17 J=IROW,NROW
X                C(L)=0.0
X                L=L+J
X   17        CONTINUE
X
X   14        NDIAG=NDIAG-IROW
X             IROW=IROW-1
X             IF(IROW.NE.0) GOTO 16
X          ENDIF
X       ENDIF
X       
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ************************ SUBROUTINE TIE  *****************************
XC      **********************************************************************
XC
X       SUBROUTINE TIE(ID,X,Y,NTOT,NVAR,IL,IM,X1,MVAR)
XC
XC      *       CHECKS FOR THE EXISTENCE OF TIED DATA POINTS. IF A POINT     *
XC      *       IS NOT TIED THE IT SETS IL(I)=1 AND IM(I)=1.                 *
XC      * INPUT                                                              *
XC      *       ID(I)    :  INDICATOR OF CENSORING                           *
XC      *       X(J,I)   :  INDEPENDENT VARIABLES                            *
XC      *       Y(I)     :  DEPENDENT VARIABLE                               *
XC      *       NTOT     :  NUMBER OF DATA POINTS                            *
XC      *       NVAR     :  NUMBER OF INDEPENDENT VARIABLES                  *
XC      *                                                                    *
XC      * OUTPUT   :                                                         *
XC      *       ID, X, AND Y IN  ORDER SO THAT DETECTED POINTS ARE           *
XC      *       LOCATED BEFORE CENSORED POINTS IF THEY ARE TIED.             *
XC      *       IL(I)    :  INDICATOR OF TIES.                               *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION ID(NTOT),X(MVAR,NTOT),Y(NTOT),IL(NTOT)
X       DIMENSION IM(NTOT),X1(MVAR)
XC
X       IL(1)=1
X       IM(1)=1
X       IL(2)=1
X       IM(2)=1
X       I=2
X       J=1
X  200  IF(Y(I).EQ.Y(I-1)) THEN
X          IL(I)=IL(I-1)+1
X          J=J+1
X          DO 300 K=1,J
X             L=I+1-K
X             IM(L)=J
X  300     CONTINUE
X       ELSE
X          J=1
X       ENDIF
X       IF(I.LT.NTOT) THEN
X          I=I+1
X          IL(I)=1
X          IM(I)=1
X          GOTO 200
X       ENDIF
XC
XC      *     IF TIED DATA CONTAINS CENSORED POINTS, ORDER THEM SO THAT      *
XC      *     A DETECTED POINT COMES FIRST.                                  *
XC
X       I=1
X  550  J=1
X       IF(IM(I).NE.1) THEN
X  600     IF(ID(I+J-1).NE.0) THEN
X             IF(J.GE.IM(I)) GOTO 800
X             J=J+1
X             GOTO 600
X          ENDIF
X          IF(J.NE.1) THEN
XC
XC      *       EXCHANGE THE DETECTED POINT AND THE CENSORED POINT           *
XC                                       
X             ID1=ID(I+J-1)
X             DO 700 L=1,NVAR
X                X1(L)=X(L,I+J-1)
X  700        CONTINUE
X
X             ID(I+J-1)=ID(I)
X             DO 710 L=1,NVAR
X                X(L,I+J-1)=X(L,I)
X  710        CONTINUE
X
X             ID(I)=ID1
X             DO 720 L=1,NVAR
X                X(L,I)=X1(L)
X  720        CONTINUE
X
X          ENDIF
X       ENDIF
X  800  IF(I.LT.NTOT) THEN
X          I=I+J
X          GOTO 550
X       ENDIF
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE TWOKM  ******************************
XC      **********************************************************************
XC
X       SUBROUTINE TWOKM(IND,XX,YY,NTOT,MX,MY,ISKIP,IPRINT,ICENS,XBIN,
X     +                  YBIN,XORG,YORG,OUTPUT,TOL,MAX,NLAST,IRAND,
X     +                  NN1,NN2,NN3,NN4,NN5,NN6,NN7,NN8,X,Y,NP,XB,YB,F,
X     +                  IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,IBWRK6,
X     +                  IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,IWRK2,WRK1,
X     +                  WRK2,SWRK1,DWRK1,AWRK1,IB,MVAR)
XC
XC      *                                                                    *
XC      *          THIS PROGRAM COMPUTES LINEAR REGRESSION COEFFICIENTS      *
XC      *      ALPHA AND BETA (INTERCEPT AND SLOPE) BY SCHMITT'S BINNED      *
XC      *      METHOD. BECAUSE OF THE BINNING METHOD, FINER BINNING GIVES    *
XC      *      BETTER RESULTS, THOUGH THE CPU TIME MAY INCREASE VERY MUCH.   *
XC      *                                                                    *
XC      *             WARNING   WARNING   WARNING   WARNING                  *
XC      *                                                                    *
XC      *    THE USER SHOULD BE WARNED THAT THIS PROGRAM ACTUALLY            *
XC      *    CHANGES THE DATA!!  FIRST, IT REDEFINES SOME LIMITS TO          *
XC      *    DETECTIONS.  IF THE BINS ARE CHOSEN TO BE TOO NARROW, THEN      *
XC      *    VIRTUALLY ALL LIMITS COULD BE CHANGED.  SECOND, IT PUSHES       *
XC      *    EACH LIMIT INTO THE ADJACENT BIN.  IF THE BINS ARE CHOSEN TO    *
XC      *    TO BE TOO WIDE, THIS SUBSTANTIALLY ALTERS THE MEASURED VALUES.  *
XC      *    THUS, THE USER MUST TREAD A FINE LINE IN CHOSING BIN SIZES.     *
XC      *                                                                    * 
XC      *        THE ERROR ANALYSIS IS DONE BY BOOTSTRAPPING METHODS. IF THE *
XC      *    NUMBER OF DATA POINTS IS MUCH SMALLER THAN 100, THE NUMBER OF   *
XC      *    ITERATIONS SHOULD BE (TOTAL NUMBER)**2.                         *
XC      *                                                                    *
XC      *      INPUT                                                         *
XC      *            X(I): INDEPENDENT VARIABLE                              *
XC      *            Y(I): DEPENDENT VARIABLE                                *
XC      *           NP(I): INDICATOR OF CENSORED STATUS                      *
XC      *                  IF NP(I)=0  : DETECTION                           *
XC      *                          =1  : Y(I) IS LOWER LIMIT                 *
XC      *                          =2  : X(I) IS LOWER LIMIT                 *
XC      *                          =3  : BOTH ARE LOWER LIMITS               *
XC      *                          =4  : Y(I) IS LOWER AND X(I) IS UPPER     *
XC      *                 FOR THE UPPER LIMITS, CHANGE THE SIGNS             *
XC      *           NTOT : NUMBER OF DATA POINTS                             *
XC      *          MPLONE: NUMBER OF PARAMETERS TO BE ESTIMATED              *
XC      *                  (IN THIS PROGRAM, MPLONE IS ALWAYS 3)             *
XC      *          MAXITS: MAXIMUM NUMBER OF ITERATIONS                      *
XC      *           ALH  : DUMMY                                             *
XC      *        DELTA(I): DELTA(2) CONTAINS TOLERANCE FOR THE               * 
XC      *                  COMPUTATION OF F'S.                               *
XC      *            MX  : NUMBER OF BINS IN THE INDEPENDENT VARIABLE        *
XC      *            MY  : NUMBER OF BINS IN THE DEPENDENT VARIABLE          *
XC      *          ISKIP : INDICATOR FOR THE BINNING. IF ISKIP=0, THE        *
XC      *                  SUBROUTINE BIN WILL COMPUTE THE INFORMATION       *
XC      *                  ABOUT THE BINNING                                 *
XC      *                  IF ISKIP>0, THE BINNING INFORMATION (BIN SIZES    *
XC      *                  ORIGIN) MUST BE PROVIDED.                         *
XC      *         IPRINT : INDICATOR FOR PRINTING. IF IPRINT>0, THE FINAL    *
XC      *                  ESTIMATIONS OF TWO DIMENSIONAL PL ESTIMATOR       *
XC      *                  WILL BE PRINTED.                                  *
XC      *         ICENS  : IF THE DATA SET IS KNOWN TO :                     *
XC      *                    CONTAIN LOWER LIMITS ONLY,   ICENS>0            *
XC      *                    CONTAIN UPPER LIMITS ONLY,   ICENS<0            *
XC      *                    CONTAIN MIXED LIMITS,        ICENS=0            *
XC      *                     OR NOT KNOWN                ICENS=0            *
XC      *         IF ISKIP>0, THE NEXT VALUES MUST BE PROVIDED :             *
XC      *           XBIN : THE BIN SIZE FOR THE X AXIS                       *
XC      *           YBIN : THE BIN SIZE FOR THE Y AXIS                       *
XC      *           XORG : THE ORIGIN OF X                                   *
XC      *           YORG : THE ORIGIN OF Y                                   *
XC      *       NN1=NC1  : NUMBER OF Y LOWER LIMITS                          *
XC      *       NN2=NC2  : NUMBER OF X LOWER LIMITS                          *
XC      *       NN3=NC3  : NUMBER OF DOUBLE LOWER LIMITS                     *
XC      *       NN4=NC4  : NUMBER OF Y LOWER, X UPPER LIMITS                 *
XC      *       NN5=NC5  : NUMBER OF Y UPPER LIMITS                          *
XC      *       NN6=NC6  : NUMBER OF X UPPER LIMITS                          *
XC      *       NN7=NC7  : NUMBER OF DOUBLE UPPER LIMITS                     *
XC      *       NN8=NC8  : NUMBER OF Y UPPER, XLOWER LIMITS                  *
XC      *         TOL    : TOLERANCE LEVEL                                   *
XC      *         MAX    : MAXIMUM NUMBER OF ITERATIONS                      *
XC      *         NLAST  : NUMBER OF ITERATIONS FOR THE BOOTSTRAPPING        *
XC      *                  IF NLAST = 0, THE PROGRAM SKIPS THE BOOTSTRAPPING *
XC      *         IRND   : A SEED FOR THE RANDOM NUMBER; IT IS ALSO USED AS  *
XC      *                  A RANDOM NUMBER ITSELF.                           *
XC      *         IB     : DIMENSION SIZE FOR BINS                           *
XC      *         ILARGE : DIMENSION > MX*MY                                 *
XC      *                                                                    *
XC      *     WORK                                                           *
XC      *        F(I,J)  : NUMBER OF DATA POINTS IN THE BIN(I,J)             *
XC      *                   (WEIGHTED BY CENSORED POINTS)                    *
XC      *        SUM     : TOTAL NUMBER OF POINTS. THIS APPROXIMATLY         *
XC      *                  EQUALS NTOT. THE DIFFERENCE BETWEEN THE TWO       *
XC      *                  VALUES DEPENDS ON THE TOLERANCE LEVEL.            *
XC      *                                                                    *
XC      *     OUTPUT                                                         *
XC      *          ALPHA : INTERCEPT COEFFICIENT                             *
XC      *           BETA : SLOPE COEFFICINET                                 *
XC      *                                                                    *
XC      *     SUBROUTINES                                                    *
XC      *         SCHMIT : THE SUBROUTINE WHICH COMPUTES THE SCHMITT'S LINEAR*
XC      *                  REGRESSION COEFFICIENTS                           *
XC      *         RAN1   : UNIFORM RANDOM NUMBER GENERATOR                   *
XC      *                                                                    *
XC      *      REF: SCHMITT, J. H. M. M. 1985, ASTROPHYS. J. 293, 178.       *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION XX(MVAR,NTOT),IND(NTOT),YY(NTOT)
X       DIMENSION X(NTOT),Y(NTOT),NP(NTOT),XB(IB),YB(IB)
X
X       DIMENSION IBWRK1(IB,IB),IBWRK2(IB,IB),IBWRK3(IB,IB)
X       DIMENSION IBWRK4(IB,IB),IBWRK5(IB,IB),IBWRK6(IB,IB)
X       DIMENSION IBWRK7(IB,IB),IBWRK8(IB,IB),IBWRK9(IB,IB)
X       DIMENSION BWRK1(IB,IB),F(IB,IB),AWRK1(5,IB)
X       DIMENSION IWRK1(NTOT),IWRK2(NTOT),WRK1(NTOT),WRK2(NTOT)
X       DIMENSION SWRK1(MVAR),DWRK1(MVAR,NTOT)
X
X       CHARACTER*9 OUTPUT
X       COMMON /C1/NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8
XC
X       ISEED = IRAND
X
X       IF(OUTPUT.EQ.'         ') THEN
X          PRINT * 
X          PRINT 80
X          PRINT *
X       ELSE
X          WRITE(60,70)
X          WRITE(60,80)
X          WRITE(60,70)
X       ENDIF
X
X   70  FORMAT('   ')
X   80  FORMAT(5X,'LINEAR REGRESSION BY SCHMITT`S METHOD')
XC
X       NC1=NN1
X       NC2=NN2
X       NC3=NN3
X       NC4=NN4
X       NC5=NN5
X       NC6=NN6
X       NC7=NN7
X       NC8=NN8
XC     
X       DO 90 I=1,NTOT
X          X(I)=XX(1,I)
X          NP(I)=IND(I)
X          Y(I)=YY(I)
X   90  CONTINUE
XC
XC      *      COMPUTE SCHMITT'S LINEAR REGRESSION COEFFICIENTS             *
XC
XC
X       CALL SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
X     +             YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
X     +             ALPHA,BETA,MM,M1,M2,M3,M4,M5,M6,M7,M8,
X     +             AWRK1,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
X     +             IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
X     +             IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
X      
X       XOX = XORG
X       YOY = YORG
X       ALPHA1 = ALPHA
X
XC
X       IF(MM.NE.0) THEN
X          IF(OUTPUT .EQ. '         ') THEN
X             PRINT 2550
X             PRINT 2600
X             PRINT 2700
X             PRINT 2800,M1,M2,M3,M4
X          ELSE
X             WRITE(60,2550)
X             WRITE(60,2600)
X             WRITE(60,2700)
X             WRITE(60,2800) M1,M2,M3,M4
X          ENDIF
X       ELSE
XC
XC      *     FOR THE UPPER LIMIT CASES                                     *
XC
X          MM=M5+M6+M7+M8
X          IF(MM.NE.0) THEN
X             IF(OUTPUT .EQ. '         ') THEN
X                PRINT 2550
X                PRINT 3100
X                PRINT 3200
X                PRINT 2800,M5,M6,M7,M8
X             ELSE
X                WRITE(60,2550)
X                WRITE(60,3100)
X                WRITE(60,3200)
X                WRITE(60,2800) M5,M6,M7,M8
X             ENDIF
X          ENDIF
X       ENDIF
X
X       CFRAC = (REAL(MM)/REAL(NTOT))*100.0
X       IF(CFRAC .GE. 10.0) THEN
X          IF(OUTPUT .EQ. '         ') THEN
X             PRINT 3300,CFRAC
X          ELSE
X             WRITE(60,3300) CFRAC
X          ENDIF
X       ENDIF
XC
XC
X 2550  FORMAT(10X,'# OF CENSORED POINTS CHANGED TO DETECTIONS')
X 2600  FORMAT(15X,'(FROM LOWER LIMITS) ')
X 2700  FORMAT(10X,' Y ONLY    X ONLY   BOTH    X LOWER Y UPPER')
X 2800  FORMAT(10X,4I8)
X 3100  FORMAT(15X,'(FROM UPPER LIMITS) ')
X 3200  FORMAT(10X,' Y ONLY    X ONLY   BOTH    X UPPER Y LOWER')
X 3300  FORMAT(10X,' WARNING!! THE CENSORING STATUS OF ',F4.1,
X     +        '% OF THE POINTS ',/10X,' HAVE BEEN CHANGED!!')
X
X 3500  IF(IPRINT.GT.0) THEN
X          IF(OUTPUT.EQ.'         ') THEN
X             PRINT 3890
X             PRINT 3900
X             PRINT 4000
X          ELSE
X             WRITE(60,3890)
X             WRITE(60,3900)
X             WRITE(60,4000)
X          ENDIF
X       ENDIF
XC
X 3890  FORMAT('    ')
X 3900  FORMAT(10X,' FINAL ESTIMATION OF TWO DIMENSIONAL KM ESTIMATORS')
X 4000  FORMAT(11X,'    X AXIS         Y AXIS        NO. OF POINTS')
X 4100  FORMAT(7X,3F15.3)
X 4110  FORMAT(10X,'# OF DATA POINTS : ',I5,'   SUM OF F : ',F12.5)
XC
XC
XC      *     CHANGE THE PROBABILITY OF THE BIN TO # OF DATA POINTS PER BIN. *
XC      *     TO CHECK THE ACCURACY OF THE 2-DIMENSIONAL KAPLAN-MEIER        *
XC      *     REDISTRIBUTION, COMPARE THE SUM OF ALL F(I,J) TO THE ORIGINAL  *
XC      *     NUMBER OF DATA POINTS.  IF THEY ARE NOT EQUAL, THEN THERE IS   *
XC      *     TROUBLE!!!!!!!                                                 *
XC      *     IF IPRINT >0, PRINT THE FINAL # OF DATA POINTS (WEIGHTED)      *
XC      *     IN EACH BIN.                                                   *
XC
X 3650  DO 3800 J=1,MY
X          DO 3700 I=1,MX
X             IF(F(I,J).NE.0.0) THEN
X                F(I,J)=F(I,J)*NTOT
X                IF(IPRINT.GT.0) THEN
X                   IF(OUTPUT.EQ.'         ') THEN
X                      PRINT 4100,XB(I),YB(J),F(I,J) 
X                   ELSE
X                      WRITE(60,4100) XB(I),YB(J),F(I,J)
X                   ENDIF
X                ENDIF
X             ENDIF
X 3700     CONTINUE
X 3800  CONTINUE
XC
X       IF(OUTPUT.EQ.'         ') THEN
X          PRINT *
X          PRINT 4110,NTOT,SUM
X       ELSE
X          WRITE(60,3890)
X          WRITE(60,4110) NTOT,SUM
X       ENDIF
X
XC
XC      *  IF NLAST IS NOT 0, THEN COMPUTE THE ERRORS OF THE INTERCEPT       *
XC      *  AND SLOPE COEFFICIENTS BY THE BOOTSTRAP METHOD.                   *
XC
X       IF(NLAST .NE. 0) THEN
X          RLAST  = REAL(NLAST)
X          RLAST1 = RLAST -1.0
X          ASIGM  = 0.0
X          ASIGM2 = 0.0
X          BSIGM  = 0.0
X          BSIGM2 = 0.0
XC
XC      *    START THE BOOTSTRAPPING COMPUTATION                             *
XC
X          DO 200 ITERAT = 1, NLAST
X             DO 100 I = 1, NTOT
XC
XC      *     SUBROUTINE RAN1 IS A RANDOM NUMBER GENERATOR.IRAND IS A SEED   *
XC      *     OF THE RANDOM NUMBER.                                          *
XC      *     USING THIS RANDOM NUMBER, CREATE A HYPOTHETICAL DATA SET.      *
XC
X                R =  RAN1(IRAND)
X                RT = R
X                R =  RAN1(IRAND)
X                R = DSQRT(R*RT)
XC
XC      *    CHOOSE THE DATA POINT NUMBERED AS "L" FROM THE DATA SET         *
XC
X                IF(R.EQ.1.0) THEN
X                   L = NTOT
X                ELSE
X                   L = INT(REAL(NTOT)*R) + 1
X                ENDIF
X
X                NP(I) = IND(L)
X                X(I)  = XX(1,L)
X                Y(I)  = YY(L)
X  100        CONTINUE
XC
XC      *    COMPUTE THE COEFFICIENTS FOR THE HYPOTHETICAL DATA SET          *
XC
X             CALL SCHMIT(NTOT,MX,MY,SUM,ISKIP,ICENS,XBIN,
X     +                   YBIN,XORG,YORG,TOL,MAX,X,Y,NP,XB,YB,F,
X     +                   AL,BE,MM,M1,M2,M3,M4,M5,M6,M7,M8,
X     +                   AWRK1,IBWRK1,IBWRK2,IBWRK3,IBWRK4,IBWRK5,
X     +                   IBWRK6,IBWRK7,IBWRK8,IBWRK9,BWRK1,IWRK1,
X     +                   IWRK2,WRK1,WRK2,SWRK1,DWRK1,IB,MVAR)
XC
XC
X             ASIGM  = ASIGM  + AL
X             ASIGM2 = ASIGM2 + AL**2
X             BSIGM  = BSIGM  + BE
X             BSIGM2 = BSIGM2 + BE**2
X  200     CONTINUE
XC
XC      *                   COMPUTE THE ERRORS                               *
XC
X          ASUM = ASIGM2 - ASIGM**2/RLAST
X          BSUM = BSIGM2 - BSIGM**2/RLAST
X          SIGMA= DSQRT(ASUM/RLAST1)
X          SIGMB= DSQRT(BSUM/RLAST1)
X       ENDIF
XC
XC      *       START PRINTING THE RESULTS                                   *
XC
X       IF(NLAST .EQ. 0) THEN
X          IF(OUTPUT.EQ.'         ') THEN
X             PRINT 1710,MX,MY
X             PRINT 1715,XBIN,YBIN
X             PRINT 1720,XOX,YOY
X             PRINT 1810
X             PRINT 1790,ALPHA1
X             PRINT 1800,BETA
X          ELSE
X             WRITE(60,1710) MX,MY
X             WRITE(60,1715) XBIN,YBIN
X             WRITE(60,1720) XOX, YOY
X             WRITE(60,1810)
X             WRITE(60,1790) ALPHA1
X             WRITE(60,1800) BETA
X          ENDIF
X       ELSE
X          IF(OUTPUT.EQ.'         ') THEN
X             PRINT 1710,MX,MY
X             PRINT 1715,XBIN,YBIN
X             PRINT 1720,XOX,YOY
X             PRINT 1820,NLAST,ISEED
X             PRINT 1810
X             PRINT 1795,ALPHA1,SIGMA
X             PRINT 1805,BETA,SIGMB
X          ELSE
X             WRITE(60,1710) MX,MY
X             WRITE(60,1715) XBIN,YBIN
X             WRITE(60,1720) XOX, YOY
X             WRITE(60,1820) NLAST,ISEED
X             WRITE(60,1810)
X             WRITE(60,1795) ALPHA1,SIGMA
X             WRITE(60,1805) BETA,SIGMB
X          ENDIF
X       ENDIF
X       
X 1710  FORMAT(10X,'# OF X BINS :',I8,',  # OF Y BINS :',I8)
X 1715  FORMAT(10X,'X BINSIZE   :',F8.3,'    Y BINSIZE  :',F11.3)
X 1720  FORMAT(10X,'X ORIGIN    :',F8.3,'    Y ORIGIN   :',F11.3)
X 1790  FORMAT(7X,'INTERCEPT COEFF.   :',F10.4)
X 1795  FORMAT(7X,'INTERCEPT COEFF.   :',F10.4,'+/-',F8.4, 
X     +       '(BOOTSTRAP APPROXIMATION)')
X 1800  FORMAT(7X,'SLOPE COEFF.       :',F10.4)
X 1805  FORMAT(7X,'SLOPE COEFF.       :',F10.4,'+/-',F8.4, 
X     +       '(BOOTSTRAP APPROXIMATION)')
X 1810  FORMAT('     ')
X 1820  FORMAT(10X,'BOOTSTRAP ITERATIONS :',I4,
X     +        '     SEED :',I8)
X
X       RETURN
X       END
X
X
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE TWOST  ******************************
XC      **********************************************************************
XC
X       SUBROUTINE TWOST(Z,IND,ISTA,IS,NG1,NG2,NTOT,IPR,OUTPUT,
X     +                  M,MVAR,NDAT,FILE,
X     +                  R,XM,H,X,E1,SCORE,RISK,A,R1,R2,E,XY,ID1,ID2,
X     +                  WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,WRK9,
X     +                  WRK10,WRK11,WRK12,WRK13,WRK14,IWRK1)
XC
XC
XC      *                                                                    *
XC      *   THIS SUBROUTINE ORGANIZES SUBROUTINES WHICH ARE                  *
XC      *   RELATED TO TWO SAMPLE TESTS.                                     *
XC      *                                                                    *
XC      *   INPUT       Z(I,J)  : DATA TO BE TESTED                          *
XC      *               IND(I,J): INDICATOR OF CENSORING                     *
XC      *               ISTA(I) : INDICATOR OF GROUP                         *
XC      *                 IS    : IS-TH SUB-DATA SET                         *
XC      *                 NG1   : INDICATOR OF THE FIRST GROUP               *
XC      *                 NG2   : INDICATOR OF THE SECOND GROUP              *
XC      *             NTOT=N    : TOTAL NUMBER OF DATA POINTS                *
XC      *             IPR=IFULL : INDICATOR FOR PRINTING                     *
XC      *               OUTPUT  : NAME OF OUTPUT FILE                        *
XC      *               MVAR    : MAXIMUM NUMBER OF VARIABLES IN THE DATA SET*
XC      *                                                                    *
XC      *   OTHER VARIABLES AND PARAMETERS ARE DESCRIBED IN                  *
XC      *   EACH SUBROUTINE.                                                 *
XC      *   ALL ENTRIES OF 'COMMON' STATEMENT ARE COMPUTED IN                *
XC      *   SUBROUTINE 'AARRAY'.                                             *
XC      *                                                                    *
XC      *   SUBROUTINES                                                      *
XC      *   AARRAY     : PUTS ALL OBSERVATIONS IN ARRAY XY                   *
XC      *                AND FORMS ARRAYS ID1, ID2.                          *
XC      *   GEHAN      : COMPUTES GEHAN'S GENERALIZED                        *
XC      *                WILCOXON TEST STATISTIC BY USING                    *
XC      *                MANTEL'S PROCEDURE WITH PERMUTATION VARIANCE.       *
XC      *   WLCXN      : COMPUTES GEHAN'S GENERALIZED WILCOXON               *
XC      *                TEST STATISTIC USING TH COMPUTATIONAL FORM          *
XC      *                THE LATTA PAPER WITH HYPERGEOMETRIC VARIANCE        *
XC      *   ARISK      : COMPUTES THE RISK SET AND OTHERS                    *
XC      *   LRANK      : COMPUTES THE LOGRANK TEST STATISTIC                 *
XC      *   PWLCXN     : COMPUTES THE PETO AND PETO GENERALIZED WILCOXON     *
XC      *                TEST STATISTIC                                      *
XC      *   PETO       : COMPUTES PETO AND PRENTICE'S GENERALIZED            *
XC      *                WILCOXON TEST STATISTIC.                            *
XC      *                                                                    *
XC      *   THESE SUBROUTINES ARE MADE BASED ON THE PROGRAMS                 *
XC      *   GIVEN IN ELISA T. LEE, "STATISTICAL METHODS FOR SURVIVAL DATA    *
XC      *   ANALYSIS", 1980, LIFETIME LEARNING PUBLICATIONS (BELMONT:CA).    *
XC      *                                                                    *
XC      *                                                                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       CHARACTER*9 OUTPUT,FILE
X
X       DIMENSION Z(MVAR,NDAT),IND(MVAR,NDAT),ISTA(NTOT),R(NTOT),XM(NTOT)
X       DIMENSION H(NTOT),X(NTOT),E1(NTOT),SCORE(NTOT),RISK(NTOT),A(NTOT)
X       DIMENSION R1(NTOT),R2(NTOT),E(NTOT)
X       DIMENSION XY(NTOT),ID1(NTOT),ID2(NTOT)
X       DIMENSION WRK1(NTOT),WRK2(NTOT),WRK3(NTOT),WRK4(NTOT),WRK5(NTOT)
X       DIMENSION WRK6(NTOT),WRK7(NTOT),WRK8(NTOT),WRK9(NTOT),WRK10(NTOT)
X       DIMENSION WRK11(NTOT),WRK12(NTOT),WRK13(NTOT),WRK14(NTOT)
X       DIMENSION IWRK1(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X       LO=1
X       IF(OUTPUT.EQ.'         ') LO=0
X       IFULL=IPR
X       N=NTOT
XC
X       CALL AARRAY(Z,IND,ISTA,IS,NTOT,NDAT,MVAR,NG1,NG2,XY,
X     +             ID1,ID2,IWRK1)
X
X        IF(M.LE.1) THEN
X           IF(LO.EQ.1) THEN
X              WRITE(60,120)
X              IF(ISIGN.EQ.1) WRITE(60,480) NCOMP,NCEN
X              IF(ISIGN.EQ.-1) WRITE(60,490) NCOMP,NCEN
X              WRITE(60,500) N1
X              WRITE(60,510) N2
X           ELSE
X              PRINT *
X              IF(ISIGN.EQ.1) PRINT 480,NCOMP,NCEN
X              IF(ISIGN.EQ.-1) PRINT 490,NCOMP,NCEN
X              PRINT 500,N1
X              PRINT 510,N2
X           ENDIF
X        ENDIF
X  480   FORMAT(6X,'# OF DATA POINTS :',I4,', # OF LOWER LIMITS :',I4)
X  490   FORMAT(6X,'# OF DATA POINTS :',I4,', # OF UPPER LIMITS :',I4)
X  500   FORMAT(6X,'# OF DATA POINTS IN GROUP I  :',I4)
X  510   FORMAT(6X,'# OF DATA POINTS IN GROUP II :',I4)
XC
X       CALL ARISK(R,XM,X,E1,NG,H,XY, ID1, NTOT)
X
X       IF(IFULL.EQ.1) THEN
X           IF(LO.EQ.0) THEN
X              WRITE(6,120)
X              WRITE(6,604) 
X           ELSE
X              WRITE(60,120)
X              WRITE(60,604) 
X           ENDIF
X 604   FORMAT(1X,'DISTINCT FAILURES',9X,'R(I)',12X,'M(I)',8X,
X     +        'M(I)/R(I)',9X,'H(I)')
X
XC
XC*******      CORRECT SIGN OF X(I) BY MULTIPLYING BY ISIGN                  *
XC
X           DO 600 I=1,NG
X              ZZ=ISIGN*X(I)
X              IF(LO.EQ.0) THEN
X                 WRITE(6,120)
X                 WRITE(6,601) ZZ,R(I),XM(I),E1(I),H(I)
X              ELSE
X                 WRITE(60,120)
X                 WRITE(60,601) ZZ,R(I),XM(I),E1(I),H(I)
X              ENDIF
X  600      CONTINUE
X  120      FORMAT('     ')
X  601      FORMAT(2X,5F15.3)
X
X           NN1=NG+1
X           IF(LO.EQ.0) THEN
X              WRITE(6,602) H(NN1)
X              WRITE(6,120)
X           ELSE
X              WRITE(60,602) H(NN1)
X              WRITE(60,120)
X           ENDIF
X  602      FORMAT(62X,F15.3)
X
X        ENDIF
X
XC  605  LTEST=1
XC
XC***************   GEHAN TEST PRINTOUT -- PERMUTATION VARIANCE               *
XC
X        IF(LO .EQ. 0) THEN
X           WRITE(6,120)
X           WRITE(6,125)
X           WRITE(6,120)
X        ELSE
X           WRITE(60,120)
X           WRITE(60,125)
X           WRITE(60,120)
X        ENDIF
X 125    FORMAT(8X,'GEHAN`S GENERALIZED WILCOXON TEST',
X     +         ' -- PERMUTATION VARIANCE')
X
X        CALL GEHAN(R1,R2,TEST,PROB,XY,ID1,ID2,NTOT)
X
X        IF(IFULL .EQ. 1) THEN
X           IF(LO .EQ. 0) THEN
X              WRITE(6,120)
X              WRITE(6,610)
X           ELSE
X              WRITE(60,120)
X              WRITE(60,610)
X           ENDIF
X
X 610       FORMAT(//17X,'T(I)',7X,'ID1(I)',4X,'ID2(I)',3X,
X     +            'R1(I)',4X,'R2(I)')
X
X           DO 611 I = 1, NCOMP
XC
XC********  CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN
XC
X              ZZ = REAL(ISIGN)*XY(I)
X              IF(LO .EQ. 0) THEN 
X                 WRITE(6,612) ZZ,ID1(I),ID2(I),R1(I),R2(I)
X              ELSE
X                 WRITE(60,612) ZZ,ID1(I),ID2(I),R1(I),R2(I)
X              ENDIF
X 611       CONTINUE
X
X 612       FORMAT(5X,F15.3,2I10,2F10.1)
X        ENDIF
X
X        IF(LO .EQ. 0) THEN
X           WRITE(6,660) TEST
X           WRITE(6,665) PROB
X        ELSE
X           WRITE(60,660) TEST
X           WRITE(60,665) PROB
X        ENDIF
X
XC
XC***************   GEHAN TEST PRINTOUT -- HYPERGEOMETRIC VARIANCE            *
XC
XC       IF(IPROG.EQ.1) THEN
X          IF(LO.EQ.0) THEN
X             WRITE(6,120)
X             WRITE(6,130)
X             WRITE(6,120)
X          ELSE
X             WRITE(60,120)
X             WRITE(60,130)
X             WRITE(60,120)
X          ENDIF
X  130     FORMAT(8X,'GEHAN`S GENERALIZED WILCOXON TEST',
X     +           ' -- HYPERGEOMETRIC VARIANCE')
XC
X          CALL WLCXN(ID1,ID2,XY,NTOT,TEST,PROB,WRK1,WRK2,WRK3,WRK4,
X     +               WRK5,WRK6,WRK7,WRK8,WRK9,WW,VAR)
X
X          IF(IFULL.EQ.1) THEN
X             IF(LO.EQ.0) THEN
X                WRITE(6,641) WW
X                WRITE(6,642) VAR
X                WRITE(6,120)
X                WRITE(6,640)
X              ELSE
X                WRITE(60,641) WW
X                WRITE(60,642) VAR
X                WRITE(60,120)
X                WRITE(60,640)
X              ENDIF
X 640   FORMAT(1X,'            ',//15X,'T(I)',9X,'SCORE(I)')
X 641   FORMAT(10X,'WW  = ',F15.3)
X 642   FORMAT(10X,'VAR = ',F15.3)
XC
XC******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
XC******       CHANGE AFFECTS ONLY THE PRINTING.                             *
XC
X             J = 1
X             DO 646 I=1,NCOMP
X
X 645            IF(J .GT. NCOMP) THEN
X                   GOTO 646
X                ENDIF
X
X                IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
X                   J=J+1
X                   GOTO 645 
X                ENDIF
X
X                IF(ID1(J) .EQ. 1) THEN
X                   J = J+1
X                   GOTO 645
X                ENDIF
X
X                ZZ=ISIGN*XY(J)
X
X                IF(WRK3(I) .NE. 0.0) THEN
X                   SCORE(I) = WRK3(I)*(WRK7(I)-
X     +                        ((WRK9(I)*WRK1(I))/WRK3(I)))
X                ELSE
X                   SCORE(I) = 0.0
X                ENDIF
X
X                IF(LO.EQ.0) THEN
X                   WRITE(6,780) ZZ,SCORE(I)
X                ELSE
X                   WRITE(60,780) ZZ,SCORE(I)
X                ENDIF
X                J=J+1
X 646         CONTINUE
X           ENDIF
X
X  652     IF(LO.EQ.0) THEN
X             WRITE(6,660) TEST
X             WRITE(6,665) PROB
X          ELSE
X             WRITE(60,660) TEST
X             WRITE(60,665) PROB
X          ENDIF
X  660     FORMAT(10X,'TEST STATISTIC        =',F12.3)
X  665     FORMAT(10X,'PROBABILITY           =',F13.4,/)
X
XC
XC************            LOGRANK TEST PRINTOUT                             *
XC
XC       ELSEIF(IPROG.EQ.2) THEN
X          IF(LO.EQ.0) THEN
X             WRITE(6,120)
X             WRITE(6,150)
X             WRITE(6,120)
X          ELSE
X             WRITE(60,120)
X             WRITE(60,150)
X             WRITE(60,120)
X          ENDIF
X  150     FORMAT(8X,'LOGRANK TEST ')
XC
X          CALL LRANK(ID1,ID2,XY,NTOT,TEST,PROB,WRK1,WRK2,WRK3,WRK4,
X     +               WRK5,WRK6,WRK7,WRK8,WRK9,WW,VAR)
X
X          IF(IFULL.EQ.1) THEN
X             IF(LO.EQ.0) THEN
X                WRITE(6,641) WW
X                WRITE(6,642) VAR
X                WRITE(6,120)
X                WRITE(6,640)
X             ELSE
X                WRITE(60,641) WW
X                WRITE(60,642) VAR
X                WRITE(60,120)
X                WRITE(60,640)
X             ENDIF
X 
XC
XC******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
XC******       CHANGE AFFECTS ONLY THE PRINTING.                             *
XC
X             J = 1
X             DO 700 I=1,NCOMP
X
X 695            IF(J .GT. NCOMP) THEN
X                   GOTO 700
X                ENDIF
X
X                IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
X                   J=J+1
X                   GOTO 695 
X                ENDIF
X
X                IF(ID1(J) .EQ. 1) THEN
X                   J = J+1
X                   GOTO 695
X                ENDIF
X
X                ZZ=ISIGN*XY(J)
X
X                IF(WRK3(I) .NE. 0.0) THEN
X                   SCORE(I) = WRK7(I)-((WRK9(I)*WRK1(I))/WRK3(I))
X                ELSE
X                   SCORE(I) = 0.0
X                ENDIF
X
X                IF(LO.EQ.0) THEN
X                   WRITE(6,780) ZZ,SCORE(I)
X                ELSE
X                   WRITE(60,780) ZZ,SCORE(I)
X                ENDIF
X                J=J+1
X 700         CONTINUE
X          ENDIF
XC
X 703      IF(LO.EQ.0) THEN
X             WRITE(6,660) TEST
X             WRITE(6,665) PROB
X          ELSE
X             WRITE(60,660) TEST
X             WRITE(60,665) PROB
X          ENDIF
X
XC
XC****************  PETO-PETO PRINTOUT
XC
X
X          IF(LO .EQ. 0) THEN
X             WRITE(6,120)
X             WRITE(6,175)
X             WRITE(6,120)
X          ELSE
X             WRITE(60,120)
X             WRITE(60,175)
X             WRITE(60,120)
X          ENDIF
X
X 175      FORMAT(8X,'PETO & PETO GENERALIZED WILCOXON TEST')
X          
X          CALL PWLCXN(H,XM,SCORE,TEST,PROB,IWLCXN,XY,ID1,ID2,NTOT)
X
X          IF(IFULL .EQ. 1) THEN
X             IF(LO .EQ. 0) THEN
X                WRITE(6,120)
X                WRITE(6,770)
X             ELSE
X                WRITE(60,120)
X                WRITE(60,770)
X             ENDIF
X
X 770         FORMAT(27X,//15X,'T(I)',9X,'SCORE(I)')
X
XC***********    CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN
X
X             DO 775 I=1,NCOMP
X                IF(XY(I) .EQ. 0.0) GOTO 775
X                ZZ = REAL(ISIGN)*XY(I)
X                IF(LO .EQ. 0) THEN
X                   WRITE(6,780) ZZ,SCORE(I)
X                ELSE
X                   WRITE(60,780) ZZ,SCORE(I)
X                ENDIF
X 775         CONTINUE
X 780         FORMAT(5X,2F15.3)
X          ENDIF
X
X          IF(LO .EQ. 0) THEN
X             WRITE(6,660) TEST
X             WRITE(6,665) PROB
X          ELSE
X             WRITE(60,660) TEST
X             WRITE(60,665) PROB
X          ENDIF
X
XC
XC****************  PETO-PRENTICE PRINTOUT                                    *
XC
XC  180  ELSEIF(IPROG.EQ.4)THEN
X          IF(LO.EQ.0) THEN
X             WRITE(6,120)
X             WRITE(6,190)
X             WRITE(6,120)
X          ELSE
X             WRITE(60,120)
X             WRITE(60,190)
X             WRITE(60,120)
X          ENDIF
X  190     FORMAT(8X,'PETO & PRENTICE GENERALIZED WILCOXON TEST')
X          IF(NCEN .NE. 0) THEN
X             CALL PETO(ID1, ID2, XY, NTOT, TEST,PROB, 
X     +                 WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
X     +                 WRK9,WRK10,WRK11,WRK12,WRK13,WRK14,WW,VAR)
X
X             IF(IFULL.EQ.1) THEN
X                IF(LO.EQ.0) THEN
X                   WRITE(6,641) WW
X                   WRITE(6,642) VAR
X                   WRITE(6,120)
X                   WRITE(6,640)
X                ELSE
X                   WRITE(60,641) WW
X                   WRITE(60,642) VAR
X                   WRITE(60,120)
X                   WRITE(60,640)
X                ENDIF
X
XC
XC******       CORRECT THE SIGN OF XY(I) BY MULTIPLYING BY ISIGN. THIS       *
XC******       CHANGE AFFECTS ONLY THE PRINTING.                             *
XC
X                J = 1
X                DO 782 I=1,NCOMP
X
X 781               IF(J .GT. NCOMP) THEN
X                      GOTO 782
X                   ENDIF
X
X                   IF(J .GT. 1 .AND. XY(J) .EQ. XY(J-1)) THEN
X                      J=J+1
X                      GOTO 781 
X                   ENDIF
X
X                   IF(ID1(J) .EQ. 1) THEN
X                      J = J+1
X                      GOTO 781
X                   ENDIF
X
X                   ZZ=ISIGN*XY(J)
X
X                   SCORE(I) = (2.0*WRK10(I)-1.0)*WRK7(I)+
X     +                        (WRK10(I)-1.0)*WRK8(I)
X
X                   IF(LO.EQ.0) THEN
X                      WRITE(6,780) ZZ,SCORE(I)
X                   ELSE
X                      WRITE(60,780) ZZ,SCORE(I)
X                   ENDIF
X                   J=J+1
X 782            CONTINUE
X              ENDIF
X
X             IF(LO.EQ.0) THEN
X                WRITE(6,660) TEST
X                WRITE(6,665) PROB
X             ELSE
X                WRITE(60,660) TEST
X                WRITE(60,665) PROB
X             ENDIF
X
X          ELSE
X  785        IF(LO.EQ.0) THEN
X                WRITE(6,790) 
X                WRITE(6,791)
X             ELSE
X                WRITE(60,790) 
X                WRITE(60,791)
X             ENDIF
X  790  FORMAT(5X,'NO CENSORED OBS., PETO & PRENTICE WILCOXON TEST')
X  791  FORMAT(5X, '       REDUCED TO GEHAN`S WILCOXON TEST')
XC
X          ENDIF
X
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************** SUBROUTINE UNIVAR *****************************
XC      **********************************************************************
XC
X       SUBROUTINE UNIVAR(IBACK)
XC
XC      *      NDAT     : DIMENSION DECLARATION                              *
XC      *                                                                    *
XC      *   UNIVARIATE PROBLEMS                                              *
XC      *   PARAMETERS                                                       *
XC      *       MVAR    :  THE MAXIMUM NUMBER OF VARIABLES ALLOWED IN A DATA *
XC      *                    SET.                                            *
XC      *       NDAT    :  THE MAXIMUM NUMBER OF DATA POINTS ALLOWED IN A    *
XC      *                    DATA SET.                                       *
XC      *       IBIN    :  THE DIMENSION SIZE FOR BINS - USED IN THE KAPLAN- *
XC      *                    MEIER ESTIMATION PROCEDURE.                     *
XC      *   COMMON FOR KAPLAN-MEIER AND TWO SAMPLE TESTS :                   *
XC      *       FILE    :  NAME OF DATA FILE    (9 LETTERS)                  *
XC      *       TITLE   :  TITLE OF THE PROBLEM (80 LETTERS)                 *
XC      *       IUNI    :  INDICATOR OF PROBLEM                              *
XC      *                    IF 1 : KAPLAN-MEIER ESTIMATOR                   *
XC      *                    IF 2 : TWO-SAMPLE TESTS                         *
XC      *       NTOT    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
XC      *       NVAR    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    * 
XC      *       ITEST   :  NUMBER OF VARIABLES TO BE TESTED (<=NVAR)         *
XC      *       ICOL    :  INDICATOR OF SAMPLES                              *
XC      *       COLM    :  NAME OF THE SAMPLE SETS                           *
XC      *       COMMAND :  NAME OF COMMAND FILE                              *
XC      *       OUTPUT  :  NAME OF OUTPUT FILE                               *
XC      *       IND(I,J):  INDICATOR OF CENSORING (J-TH DATA POINTS OF I-TH  *
XC      *                  VARIABLE)                                         *
XC      *       X(I,J)  :  DATA POINTS                                       *
XC      *                                                                    *
XC      *    KAPLAN-MEIER ESTIMATOR                                          *
XC      *       IKM     :  INDICATOR OF PRINTOUT                             *
XC      *                    IF 0 : MEAN AND ERROR ONLY                      *
XC      *                    IF 1 : MEAN, ERROR, SURVIVAL DISTRIBUTION       *
XC      *                                                                    *
XC      *    TWO-SAMPLE TESTS                                                *
XC      *       IPROG(I):  INDICATOR OF TEST (I=1,5; OR 6 FOR EXIT)          *
XC      *       NGROUP  :  NUMBER OF GROUPS                                  *
XC      *      IGROUP(I):  INDICATOR OF GROUP                                *
XC      *       IFIRST  :  INDICATOR OF FIRST GROUP                          *
XC      *       ISECON  :  INDICATOR OF SECOND GROUP                         *
XC      *       IFULL   :  INDICATOR OF PRINTOUT FOR THE I-TH COMBINATION    *
XC      *       GROUP(I):  NAME OF THE GROUPS                                *
XC      *        LKM    :  INDICATOR OF K-M ESTIMATOR                        *
XC      *        IKM    :  INDICATOR OF PRINTOUT.                            *
XC      *       NOTEST  :  NUMBER OF TESTS TO BE USED                        *
XC      *       ISTA(I) :  INDICATOR OF GROUPS                               *
XC      *                                                                    *
XC      *    WRK VARIABLES AND ARRAYS                                        *
XC      *       LCOMM   :  INDICATOR OF USE OF COMMAND FILE                  *
XC      *                      IF 0, READ INFORMATION FROM THE TERMINAL      *
XC      *                      IF 1, READ INFORMATION FROM THE COMMAND FILE  *
XC      *       CHECK   :  READER OF Y/N QUESTIONS                           *
XC      *       NTOT    :  NUMBER OF DATA POINTS                             *
XC      *       NCHANGE :                                                    *
XC      *       ICHANGE :                                                    *
XC      *       FGROUP  :                                                    *
XC      *       SGROUP  :                                                    *
XC      *      CHAR(I,J):  READ-IN FORM OF SEVERAL INPUTS                    *
XC      *     CTEST(I,1):  READ-IN FORM OF NOTEST                            *
XC      *     CPROG(I,J):  READ-IN FORM OF IPROG(J)                          *
XC      *      CCOL(I,J):  READ-IN FORM OF ICOL(J)                           *
XC      *      IFIRST   :                                                    *
XC      *      ISECON   :                                                    *
XC      *      CF(I,1)  :  READ-IN FORM OF IFIRST                            *
XC      *      CS(I,1)  :  READ-IN FORM OF ISECON                            *
XC      *     CIKM1(I,J):  READ-IN FORM OF IKM1(J)                           *
XC      *       N1      :  NUMBER OF DATA POINTS IN THE FIRST GROUP          *
XC      *       N2      :  NUMBER OF DATA POINTS IN THE SECOND GROUP         *
XC      *      JIND(I,J):  IND(I,J) IN FIRST OR SECOND GROUP                 *
XC      *       Z(I,J)  :  X(I,J) IN FIRST OR SECOND GROUP                   *
XC      *                                                                    *
XC      *   NOTE:                                       *
XC      *     ALL VARIABLES WRK_, IWRK_, DWRK1, AND SWRK1 ARE USED TO COLLECT*
XC      *     ALL DIMENSION DECLARATIONS TO THE TOP LEVEL SUBROUTINE.        *
XC      *     THEY DO NOT DIRECTLY AFFECT ANY OTHER SUBROUTINES. IF ONE NEEDS*
XC      *     TO KNOW WHAT KIND OF WORK A VARIABLE DOES, GO TO THE LOWER     *
XC      *     SUBROUTINES AND READ THE DESCRIPTION OF THE VARIABLE.          *
XC      *                                                                    *
XC      *   SUBROUTINES                                                      *
XC      *     DATA1, DATA2, DATAIN, KMESTM, TWOST                            *
XC      *                                                                    *
XC
XC
XC      *                                                                    *
XC      *  START UNIVARIATE PROBLEM: K-M ESTIMATOR OR TWO-SAMPLE TESTS       *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       
XC      *   THIS PARAMETER STATEMENT AND THE ONE IN BIVAR.F ARE ALL THAT     *
XC      *   NEEDS TO BE CHANGED IF THE USER WISHES TO WORK ON DATA SETS      *
XC      *   OF MORE THAN 500 OBSERVATIONS OR WITH MORE THAN FOUR VARIABLES.  *
X
XC  **************************************************************************
X       PARAMETER(MVAR=4, NDAT=500, IBIN=50)
XC  **************************************************************************
X
X       CHARACTER*1 CHECK,CHAR(4,10)
X       CHARACTER*1 CF(4,1),CS(4,1),CIKM(4,1),CIS4(4,1)
X       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM,GROUP(10)
X       CHARACTER*80 TITLE 
X
X
X       DIMENSION IND(MVAR,NDAT),X(MVAR,NDAT)
X       DIMENSION JIND(MVAR,NDAT),Z(MVAR,NDAT),ISTA(NDAT)
X       DIMENSION IGROUP(MVAR),JGROUP(MVAR)
XC
XC     *     THE DIMENSIONS BELOW WERE ADDED TO COLLECT ALL DIMENSION        *
XC     *     DECLARATIONS IN THE TOP SUBROUTNE                        *
XC
X       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
X       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
X       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
X       DIMENSION WRK12(NDAT),WRK13(NDAT),WRK14(NDAT)
X       DIMENSION WRK15(NDAT),WRK16(NDAT),WRK17(NDAT)
X       DIMENSION WRK18(NDAT),WRK19(NDAT),WRK20(NDAT)
X       DIMENSION WRK21(NDAT),WRK22(NDAT),WRK23(NDAT)
X       DIMENSION WRK24(NDAT),WRK25(NDAT),WRK26(NDAT)
X
X       DIMENSION IWRK1(NDAT), IWRK2(NDAT), IWRK3(NDAT)
X       DIMENSION BWRK1(IBIN),BWRK2(IBIN), BWRK3(IBIN)
X       DIMENSION DWRK1(MVAR,NDAT), SWRK1(MVAR)
XC
XC
X   50  FORMAT(A1)
XC
X 1000  PRINT *
X       PRINT *,'    SELECT PROBLEM: ' 
X       PRINT *,'     1 KAPLAN-MEIER DISTRIBUTION '
X       PRINT *,'     2 TWO SAMPLE TESTS'
X       PRINT *,'     3 EXIT '
X       PRINT *
X       PRINT *,'    (IF YOU CHOOSE OPTION 2, YOU CAN STILL DO 1 LATER) '
X       PRINT *
X 1010  WRITE(6,1020)
XC
XC      *         SELECT PROBLEM                                             *
XC
X 1020  FORMAT(' CHOICE? ')
XC
X       CALL DATA1(IUNI)
XC
X       IF((IUNI.EQ.1).OR.(IUNI.EQ.2).OR.(IUNI.EQ.3)) GOTO 1030
X       PRINT *
X       PRINT *,'      PLEASE TYPE ONCE MORE'
X       GOTO 1010
XC
X 1030  IF(IUNI.EQ.3) STOP 
XC
X 1120  IF(IUNI.EQ.2) GOTO 1140
X
X       PRINT *
X       PRINT *,'   ***  KAPLAN-MEIER ESTIMATOR  ***'
X       PRINT *
X       GOTO 1330
XC
XC
XC      *     DISPLAY THE INFORMATION ABOUT TWO SAMPLE TESTS                 *
XC
X 1140  PRINT *
X       PRINT *,'   ***     TWO-SAMPLE TESTS     ***'
X       PRINT *
X       PRINT *
XC
X       LCOMM=1
XC
XC      *   CHECK WHETHER THE DATA NEEDS TO BE READ FROM A FILE              *
XC
X 1330  PRINT *
X       PRINT *,'DO YOU WANT TO READ THE INPUTS'
X       WRITE(6,1340)
X 1340  FORMAT('     FROM A COMMAND FILE (Y/N)? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 2680
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 1350
X       GOTO 1330
XC
XC      *  READ INFORMATION COMMON TO K-M ESTIMATOR AND TWO SAMPLE TESTS     *
XC
X 1350  PRINT *
XC
XC      *               READ NAME OF THE DATA FILE                           *
XC
X 1360  PRINT *
X       WRITE(6,1370)
X 1370  FORMAT('WHAT IS THE DATA FILE NAME ? ')
X       READ(5,1380) FILE
X 1380  FORMAT(A9)
X       PRINT *
X       WRITE(6,1400) FILE
X 1400  FORMAT(5X,'THE FILE NAME IS ',A9)
X       PRINT *
XC
XC      *              READ TITLE OF THE PROBLEM                             *
XC
X 1410  PRINT *
X       WRITE(6,1420)
X 1420  FORMAT('WHAT IS THE PROBLEM TITLE? ')
X       READ(5,1430) TITLE 
X 1430  FORMAT(A80)
XC
XC      *            READ THE NUMBER OF VARIABLES                            *
XC
X       PRINT *
X       WRITE(6,1480)
X 1480  FORMAT('HOW MANY VARIABLES DO YOU HAVE? ')
X       CALL DATA1(NVAR)
XC
XC      *          CHECK WHICH VARIABLE SHOULD BE TESTED                     *
XC
X       ICOL=1
X       IF(NVAR.EQ.1) GOTO 1840 
X 1550  PRINT *
X       PRINT *,'      WHICH VARIABLE DO YOU WANT TO TEST? '
X 1560  WRITE(6,1570)
X 1570  FORMAT(' VARIABLE NUMBER: ')
X       READ(5,1580) (CHAR(I,1),I=1,4)
X 1580  FORMAT(4A1)
XC
XC      *     CHECK IF THE CHOICE IS CORRECT                                 *
XC
X       CALL DATA2(CHAR,1,1,ICOL,LIND)
X       IF(LIND.NE.0) PRINT *,
X     +               '    PLEASE TYPE IN THE VARIABLE NUMBER AGAIN'
X       IF(LIND.NE.0) GOTO 1550
X       IF(ICOL.LE.NVAR) GOTO 1840
X       PRINT *
X       PRINT *,'   THE NUMBER IS LARGER THAN THE NUMBER OF VARIABLES'
X       GOTO 1560
XC
XC      *                READ NAME OF THE VARIABLE                           *
XC
X 1840  PRINT *
X       WRITE(6,1850) ICOL
X 1850  FORMAT('VARIABLE ',I4,' IS NAMED')
X       READ(5,1380) COLM
XC
XC      *   THE NEXT FEW LINES CONCERN ONLY 2-SAMPLE TESTS                   *
XC      *   IF THE PROBLEM IS K-M ESTIMATION, GO TO LINE 2630                *
XC
XC      *          READ THE NUMBER OF GROUPS                                 *
XC
X       IF(IUNI.EQ.1) GOTO 2630
X 2030  WRITE(6,2040)
X 2040  FORMAT(/'HOW MANY GROUPS DO YOU HAVE? ')
X       CALL DATA1(NGROUP)
X       IF(NGROUP.LT.2) THEN
X            PRINT *,'      NUMBER OF GROUPS MUST BE TWO OR MORE'
X            GOTO 2030
X       ENDIF
XC
X       IF(NGROUP.EQ.2) GOTO 2180
X
XC
XC      *  IF THE NUMBER OF GROUPS IS MORE THAN TWO, SPECIFY COMBINATIONS    *
XC
X 2170  PRINT *
X       PRINT *,'     WHICH COMBINATION DO YOU WANT TO TEST? '
X 2180  PRINT *
X       WRITE(6,2190)
X 2190  FORMAT('FIRST GROUP INDICATOR  ')
X       CALL DATA1(IFIRST)
X       IGROUP(1) = IFIRST
X 2210  PRINT *
X       WRITE(6,2220)
X 2220  FORMAT('SECOND GROUP INDICATOR  ')
X       CALL DATA1(ISECON)
X 2240  IF(IFIRST.EQ.ISECON) THEN
X            PRINT *,' YOU CHOSE THE SAME GROUP.'
X            PRINT *,' PLEASE CHANGE THE SECOND GROUP.'
X            GOTO 2210
X       ENDIF
X       IGROUP(2) = ISECON
XC
XC      *         READ THE NAME OF THE GROUPS                                *
XC
X 2250  PRINT *
X       WRITE(6,2255) IFIRST   
X 2255  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
X       READ(5,1380) GROUP(1)
X       PRINT *
X       WRITE(6,2258) ISECON   
X 2258  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
X       READ(5,1380) GROUP(2)
X       PRINT *
XC
XC      *      READ WHETHER TO PRINTOUT ONLY RESULTS OR TO GIVE FULL DETAILS *
XC
X 2312  WRITE(6,2314)
X 2314  FORMAT('DO YOU WANT PRINTOUTS OF COMPUTATIONAL',
X     +       ' DETAILS (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.NE.'Y'.AND.CHECK.NE.'y') GOTO 2316
X       IFULL=1
X       GOTO 2400
X 2316  IF(CHECK.NE.'N'.AND.CHECK.NE.'n') GOTO 2312
X       IFULL=0
XC
XC      *     LKM IS SET TO ONE IN THE TWO-SAMPLE TEST ROUTINE, SO THAT      *
XC      *     THE KAPLAN-MEIER PERCENTILES AND MEAN FOR EACH GROUP ARE       *
XC      *     AUTOMATICALLY PROVIDED.                                        *
XC
X
X 2400  LKM=1
X
XC
XC      *  CHECK WHETHER THE FULL K-M ESTIMATOR IS NEEDED.                   *
XC      *  FROM THE NEXT LINE, THE INPUTS ARE COMMON FOR BOTH KAPLAN-MEIER   *
XC      *  ESTIMATION AND TWO SAMPLE TESTS.                                  *
XC
X 2630  PRINT *
X       WRITE(6,2640)
X       IKM=0
X 2640  FORMAT('DO YOU WANT TO PRINT OUT THE FULL K-M ',
X     +        'ESTIMATE (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
X          IKM=1
XC
XC       *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATE                 *
XC
X 4644    PRINT *
X         PRINT *,'DO YOU WANT TO PRINT OUT '
X         WRITE(6,2644)
X 2644    FORMAT('           THE DIFFERENTIAL FORM (Y/N)?')
X         READ(5,50) CHECK
X         IF(CHECK .EQ. 'Y'.OR.CHECK.EQ.'y') THEN
X
X           KDIFF = 1
X 2646      WRITE(6,2647)
X 2647      FORMAT(' SPECIFY STARTING VALUE : ')
X           READ(5,2648) START
X 2648      FORMAT(F10.3)
X
X 4502      WRITE(6,4503) 
X 4503      FORMAT('HOW MANY BINS DO YOU WANT? : ')
X           READ(5,4554) LSTEP
X 4554      FORMAT(I4)
X           IF(LSTEP .LE. 0) GOTO 4502
X
X 4506      WRITE(6,4507)
X 4507      FORMAT('SPECIFY BIN SIZE :')
X           READ(5,2648) BINSIZ
X           IF(BINSIZ .LE. 0.0) GOTO 4506
X
X         ELSEIF(CHECK .EQ. 'N'.OR. CHECK .EQ. 'n') THEN
X           KDIFF = 0
X         ELSE
X           GOTO 4644
X         ENDIF
X       ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
X         IKM=0
X       ELSE
X         GOTO 2630
X       ENDIF
X       
X 2650  PRINT *
X        WRITE(6,2655)
X 2655  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
X       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 3230
X       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 3230
X       GOTO 2650
XC
XC      *   IF ALL INFORMATION IS TO BE READ FROM THE TERMINAL, GOTO 3230    *
XC      *   FROM THE NEXT LINE, THE INPUTS ARE FROM "COMMAND" FILE           *
XC
XC
XC      *   READ THE NAME OF "COMMAND" FILE                                  *
XC
X 2680  PRINT *
X       WRITE(6,2690)
X 2690  FORMAT('WHAT IS THE NAME OF YOUR COMMAND FILE? ')
X       READ(5,1380) COMMAND
X       WRITE(6,2710) COMMAND
X 2710  FORMAT(5X,'YOUR COMMAND FILE IS CALLED ',A9)
XC
X       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
XC
XC      *              READ THE DATA FILE NAME                               *
XC
X 2820  READ(50,1380) FILE
XC
XC      *              READ THE TITLE OF THE PROBLEM                         *
XC
X       READ(50,1430) TITLE 
XC
XC      * READ THE NUMBER OF VARIABLES.                                      *
XC
X       READ(50,2830) (CHAR(I,1),I=1,4)
X 2830  FORMAT(12A1)
XC
XC      *      CHECK THE NUMBER OF VARIABLES                                 *
XC
X       CALL DATA2(CHAR,1,3,NVAR,LIND)
X       IF(LIND.EQ.0) GOTO 2845 
X 2835  PRINT *
X       PRINT *,'     NUMBER OF VARIABLES IS NOT READABLE'
X       CLOSE(UNIT=50)
X       STOP
X 2845  IF(NVAR.LE.0) GOTO 2835
XC
XC      *          READ WHICH VARIABLE NEEDS TO BE TESTED                    *
XC
X       READ(50,2910) (CHAR(I,1),I=1,4)
X 2910  FORMAT(4A1)
X       CALL DATA2(CHAR,1,1,ICOL,IN)
X       IF(IN.EQ.0) GOTO 2935
X 2915  PRINT *
X       WRITE(6,2920) 
X 2920  FORMAT(5X,'THE VARIABLE IS NOT READABLE')
X       CLOSE(UNIT=50)
X       STOP
X 2935  IF((ICOL.LE.0).OR.(ICOL.GT.NVAR)) GOTO 2915
X 2940  CONTINUE
XC
XC
XC      *      READ THE NAME OF THE VARIABLE                                 *
XC
X 2962  READ(50,2963) COLM
X 2963  FORMAT(10A9)
X       IF(IUNI.EQ.1) GOTO 3180
XC
XC      *   FROM THE NEXT LINE, INPUTS ARE ONLY FOR TWO SAMPLE TESTS         *
XC      *   IF IT IS A K-M ESTIMATOR PROBLEM, GO TO 3180                     *
XC
XC      *   READ THE NUMBER OF GROUPS                                        *
XC
X       READ(50,2970) (CHAR(I,1),I=1,4)
X 2970  FORMAT(4A1)
X       CALL DATA2(CHAR,1,1,NGROUP,LIND)
X       IF(LIND.EQ.0) GOTO 3000
X 2980  PRINT *,'     IT IS NOT CLEAR HOW MANY GROUPS YOU HAVE'
X       CLOSE(UNIT=50)
X       STOP
X 3000  IF(NGROUP.GT.1) GOTO 3005
X       GOTO 2980
XC
XC      *     READ THE INDICATOR OF THE GROUPS                               *
XC
X 3005  READ(50,3010) ((CHAR(I,J),I=1,4),J=1,NGROUP)
X 3010  FORMAT(60A1)
X 3020  DO 3050 J=1,NGROUP
X       CALL DATA2(CHAR,J,NGROUP,JGROUP(J),LIND)
X       IF(LIND.EQ.0) GOTO 3050
X 3025  PRINT *
X       WRITE(6,3030) J
X 3030  FORMAT(5X,'THE INDICATOR OF ',I4,'TH GROUP IS NOT CLEAR')
X       CLOSE(UNIT=50)
X       STOP
X 3050  CONTINUE
XC
XC      *  READ NUMBER OF FIRST GROUP,SECOND GROUP,                          *
XC      *       WHETHER PRINT OUT ALL OR RESULTS ONLY                        *
XC      *       WHETHER K-M ESTIMATOR IS NEEDED                              *
XC      *       WHETHER PRINT OUT ALL OR RESULTS ONLY FOR K-M                *
XC
X       READ(50,3085) (CF(I1,1),I1=1,4),(CS(I2,1),I2=1,4),
X     +              (CIS4(I3,1),I3=1,4),(CIKM(I4,1),I4=1,4)
X 3085  FORMAT(16A1)
XC
X       CALL DATA2(CF,1,1,IFIRST,LIND)
X       IF(LIND.EQ.0) GOTO 3087
X 3086  PRINT *
X       PRINT *,'   THE VALUE FOR "IFIRST" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3087  IF((IFIRST.LT.0).OR.(IFIRST.GT.NGROUP)) GOTO 3086
XC
X       CALL DATA2(CS,1,1,ISECON,LIND)
X       IF(LIND.EQ.0) GOTO 3089
X 3088  PRINT *
X       PRINT *,'   THE VALUE FOR "ISECON" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3089  IF((ISECON.LT.0).OR.(ISECON.GT.NGROUP)) GOTO 3087
X       IF(ISECON.EQ.IFIRST) GOTO 3087
X       IGROUP(1) = IFIRST
X       IGROUP(2) = ISECON
XC
X       CALL DATA2(CIS4,1,1,IFULL,LIND)
X       IF(LIND.EQ.0) GOTO 3091
X 3090  PRINT *
X       PRINT *,'    THE VALUE FOR "IFULL" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3091  IF((IFULL.EQ.0).OR.(IFULL.EQ.1)) GOTO 3092
X       GOTO 3090
XC
X 3092  LKM = 1
XC
X 3095  CALL DATA2(CIKM,1,1,IKM,LIND)
X       IF(LIND.EQ.0) GOTO 3097
X 3096  PRINT *
X       PRINT *,'    THE VALUE FOR "IKM" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3097  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 5190
X       GOTO 3096
XC
XC      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
XC
X 5190  READ(50,2970) (CHAR(I,1),I=1,4)
X       CALL DATA2(CHAR,1,1,KDIFF,LIND)
X       IF(LIND.EQ.0) GOTO 5201
X 5200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
X       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
X       CLOSE(UNIT=50)
X       STOP
X 5201  IF(KDIFF.EQ.1) GOTO 5202
X       IF(KDIFF.EQ.0) GOTO 3102
X       GOTO 5200
XC
X 5202  READ(50,4203) START
X 5203  FORMAT(F10.3)
X       READ(50,4204) LSTEP
X 5204  FORMAT(I4)
X       READ(50,4203) BINSIZ
XC
X 3102  READ(50,2970) (CHAR(I,1),I=1,4)
X       CALL DATA2(CHAR,1,1,IDATA,LIND)
X       IF(LIND.EQ.0) GOTO 3110
X 3100  PRINT *
X       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3110  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3168
X       GOTO 3100
XC
XC      *            READ THE NAME OF THE GROUPS                             *
XC
X 3168  READ(50,1380) GROUP(1)
X       READ(50,1380) GROUP(2)
XC
XC      *   READ NAME OF THE OUTPUT FILE. IF THE NAME IS BLANK, THE RESULTS  *
XC      *   WILL BE SHOWN ON THE TERMINAL.                                   *
XC
X       READ(50,1380) OUTPUT
X       IF(OUTPUT.NE.'         ') GOTO 3300
X       GOTO 3230
XC
XC      *      FROM THE NEXT LINE, INPUTS ARE ONLY FOR THE K-M ESTIMATOR     *
XC
X 3180  READ(50,3200) (CHAR(I,1),I=1,4)
X 3200  FORMAT(4A1)
X       CALL DATA2(CHAR,1,1,IKM,LIND)
X       IF(LIND.EQ.0) GOTO 3205
X 3203  PRINT *,'     IT IS NOT CLEAR WHETHER YOU WANT TO PRINT OUT ALL'
X       PRINT *,'     KM ESTIMATORS'
X       CLOSE(UNIT=50)
X       STOP
X 3205  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 4190
X       GOTO 3203
XC
XC      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
XC
X 4190  READ(50,2970) (CHAR(I,1),I=1,4)
X       CALL DATA2(CHAR,1,1,KDIFF,LIND)
X       IF(LIND.EQ.0) GOTO 4201
X 4200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
X       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
X       CLOSE(UNIT=50)
X       STOP
X 4201  IF((KDIFF.EQ.0).OR.(KDIFF.EQ.1)) GOTO 4202
X       GOTO 4200
XC
X 4202  READ(50,4203) START
X 4203  FORMAT(F10.3)
X       READ(50,4204) LSTEP
X 4204  FORMAT(I4)
X       READ(50,4203) BINSIZ
XC
XC      *    INFORMATION ABOUT PRINTOUT                                       *
XC
X 3210  READ(50,2970) (CHAR(I,1),I=1,4)
X       CALL DATA2(CHAR,1,1,IDATA,LIND)
X       IF(LIND.EQ.0) GOTO 3220
X 3215  PRINT *
X       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
X       CLOSE(UNIT=50)
X       STOP
X 3220  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3225
X       GOTO 3215
XC
XC
X 3225  READ(50,1380) OUTPUT
X       CLOSE(UNIT=50)
X       IF(OUTPUT.NE.'         ') GOTO 3300
XC
XC
XC      *         LEAVE THE "COMMAND" FILE                                   *
XC      *         CHECK OUTPUT FILE                                          *
XC
X 3230  OUTPUT='         '
X 3240  PRINT *
X       WRITE(6,3250)
X 3250  FORMAT('DO YOU WANT TO SAVE THE RESULTS IN A FILE (Y/N)? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 3260
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 3300
X       GOTO 3240
X 3260  PRINT *
X       WRITE(6,3270)
X 3270  FORMAT('WHAT IS THE NAME OF THE FILE?  ')
X       READ(5,1380) OUTPUT
XC
XC      *          READ IN DATA THOUGH THE SUBROUTINE "DATAIN"               *
XC
X 3300  CALL DATAIN(IUNI,FILE,NVAR,ISTA,IND,X,NTOT,NDAT,MVAR)
XC
X       IF(IUNI.EQ.2) GOTO 3360
XC
XC      *              COMPUTE THE K-M ESTIMATOR                             *
XC
XC
X       IF(OUTPUT.NE.'         ') 
X     +     OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
XC  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
XC  VAX/VMS MACHINES.
XC     +    ,CARRIAGECONTROL='LIST'
X     +     )
X
XC
X       CALL KMESTM(IND,X,NTOT,ICOL,IKM,TITLE,COLM,OUTPUT,IBIN,0,
X     +             KDIFF,START,BINSIZ,LSTEP,FILE,
X     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
X     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
X     +             WRK9,MVAR)
XC
X       IF(IDATA.EQ.0) GOTO 3335
X       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
X       IF(OUTPUT.EQ.'         ') WRITE(6,3331)  FILE
X       IF(OUTPUT.NE.'         ') WRITE(60,3332) 
X       IF(OUTPUT.EQ.'         ') PRINT 3332
X 3331  FORMAT(7X,' INPUT DATA FILE: ',A9)
X 3332  FORMAT(5X,'   CENSORSHIP     X ')
X       DO 3333 I=1,NTOT
X       IF(OUTPUT.NE.'         ') WRITE(60,3334) IND(ICOL,I),X(ICOL,I)
X       IF(OUTPUT.EQ.'         ') PRINT 3334,IND(ICOL,I),X(ICOL,I)
X 3333  CONTINUE
X 3334  FORMAT(12X,I4,F10.3)
XC
XC
X       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)
XC
X 3335  PRINT *
X       PRINT *,'    K-M ESTIMATOR COMPUTATIONS ARE FINISHED'
XC
XC      *      CHECK WHETHER THE USER WANTS TO USE OTHER METHODS            *
XC
X 3340  WRITE(6,3350)
X 3350  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
X       GOTO 3340
XC
XC
XC      *                COMPUTE TWO SAMPLE TESTS                            *
XC
X 3360  IF(OUTPUT.EQ.'         ') GOTO 3370
X
X       OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
XC  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
XC  VAX/VMS MACHINES.
XC     +    ,CARRIAGECONTROL='LIST'
X     +    )
X
X       WRITE(60,3380)
X       WRITE(60,3390)
X       WRITE(60,3400) TITLE  
X       WRITE(60,3390)
X       GOTO 3410
X 3370  WRITE(6,3380)
X       PRINT *
X       WRITE(6,3400) TITLE  
X       PRINT *
X 3380  FORMAT(8X,'   ******* TWO SAMPLE TEST ******')
X 3390  FORMAT('    ')
X 3400  FORMAT(8X,'TITLE : ',A80)    
XC
X 3410  IF(OUTPUT.EQ.'         ') GOTO 3420
X       WRITE(60,3430) FILE
X       WRITE(60,3435) COLM
X       GOTO 3440
X 3420  WRITE(6,3430) FILE
X       WRITE(6,3435) COLM
X 3430  FORMAT(8X,'DATA SET : ',A9)
X 3435  FORMAT(8X,'VARIABLE : ',A9)
XC
X 3440  IF(OUTPUT.EQ.'         ') GOTO 3450
X       WRITE(60,3460) GROUP(1),GROUP(2)
X       GOTO 3470
X 3450  WRITE(6,3460) GROUP(1),GROUP(2)
X 3460  FORMAT(8X,'TESTED GROUPS ARE ',A9,' AND ',A9)
XC
XC
XC 3470   DO 3480 M=1,NOTEST
XC
X
X 3470  CALL TWOST(X,IND,ISTA,ICOL,IGROUP(1),IGROUP(2),NTOT,
X     +            IFULL,OUTPUT,M,MVAR,NDAT,FILE,
X     +            WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
X     +            WRK9,WRK10,WRK11,WRK12,IWRK1,IWRK2,
X     +            WRK13,WRK14,WRK15,WRK16,WRK17,WRK18,WRK19,
X     +            WRK20,WRK21,WRK22,WRK23,WRK24,WRK25,WRK26,IWRK3)
XC
XC 3480  CONTINUE
X
XC
XC      *         IF K-M ESTIMATOR IS NOT REQUESTED, GOTO 3510               *
XC
X       IF(LKM.EQ.0) GOTO 3510
X       N1=0
X       N2=0
X       DO 3490 I=1,NTOT
X          IF(ISTA(I).NE.IGROUP(1)) GOTO 3490
X          N1=N1+1
X          JIND(ICOL,N1)=IND(ICOL,I)
X          Z(ICOL,N1)=X(ICOL,I)
X 3490  CONTINUE
XC
X       CALL KMESTM(JIND,Z,N1,ICOL,IKM,TITLE,GROUP(1),OUTPUT,
X     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
X     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
X     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
X     +             WRK9,MVAR)
XC
X       DO 3500 I=1,NTOT
X          IF(ISTA(I).NE.IGROUP(2)) GOTO 3500
X          N2=N2+1
X          JIND(ICOL,N2)=IND(ICOL,I)
X          Z(ICOL,N2)=X(ICOL,I)
X 3500  CONTINUE
XC
X       CALL KMESTM(JIND,Z,N2,ICOL,IKM,TITLE,GROUP(2),OUTPUT,
X     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
X     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
X     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
X     +             WRK9,MVAR)
X
X 3510  IF(IDATA.EQ.0) GOTO 3525
X
X       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
X       IF(OUTPUT.EQ.'         ') PRINT 3331, FILE
X       IF(OUTPUT.NE.'         ') WRITE(60,3520) 
X       IF(OUTPUT.EQ.'         ') PRINT 3520
X 3520  FORMAT(5X,'  CENSORSHIP     GROUP      X ')
X       DO 3521 I=1,NTOT
X       IF(OUTPUT.NE.'         ') WRITE(60,3522) IND(ICOL,I),ISTA(I),
X     +                                                        X(ICOL,I)
X       IF(OUTPUT.EQ.'         ') PRINT 3522,IND(ICOL,I),ISTA(I),
X     +                                                        X(ICOL,I)
X 3521  CONTINUE
X 3522  FORMAT(12X,I4,6X,I4,F10.3)
X
X       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)
X
X 3525  PRINT *
X       PRINT *,'     THE TWO SAMPLE TESTS ARE FINISHED'
X 3530  PRINT *
X       WRITE(6,3540)
X 3540  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
X       READ(5,50) CHECK
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
X       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
X       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
X       GOTO 3530
X       END
X
XC
XC      **********************************************************************
XC      ******************** SUBROUTINE UNPACK *******************************
XC      **********************************************************************
XC
X       SUBROUTINE UNPACK(X,N,LENX)
XC
XC      *      ALGORITHM AS 139.1 APPL.STATIST. (1979) VOL.28., NO.2         *
XC      *                                                                    *
XC      *      THIS SUBROUTINE EXPANDS A SYMMETRIC MATRIX STORED IN          *
XC      *      LOWER TRIANGLAR FORM IN THE FIRST N*(N+1)/2 POSITIONS         *
XC      *      OF X INTO A MATRIX USING THE FIRST N*N POSITIONS              *
XC
XC      *      LENX--THE LENGTH OF VECTOR--MUST BE NOT LESS THAN N*N         *
XC      *         (I.E. MUST NOT BE LESS THAN (NVAR+1)**2                    *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION X(LENX)
X
X       NSQ=N*N
X       II=NSQ
X       JJ=N*(N+1)/2
XC
XC      *                     STORE LAST ROW                                 *
XC
X       DO 10 I=1,N
X          X(II)=X(JJ)
X          II=II-1
X          JJ=JJ-1
X   10  CONTINUE
X
X       DO 80 I=2,N
XC
XC      *      OBTAIN UPPER PART OF MATRIX FROM PART ALREADY SHIFTED         *
XC
X          IJ=I-1
X          KK=NSQ+1-I
X          DO 50 J=1,IJ
X             X(II)=X(KK)
X             II=II-1
X             KK=KK-N
X   50     CONTINUE
XC
XC      *      OBTAIN LOWER PART OF MATRIX FROM ORIGINAL TRIANGULAR          *
XC      *      STORAGE                                                       *
XC
X          IJ=N-IJ
X          DO 70 J=1,IJ
X             X(II)=X(JJ)
X             II=II-1
X             JJ=JJ-1
X   70     CONTINUE
X   80  CONTINUE
X
X       RETURN
X       END
X
XC     
XC
XC***************************************************************************
XC**************************  SUBROUTINE WLCXN  *****************************
XC***************************************************************************
XC
XC
X       SUBROUTINE WLCXN(ID1,ID2,XY,NTOT,TEST,PROB,D,E,R, D1, E1, R1, 
X     +                  D2, E2, R2,SCORE,VAR)
XC     *
XC     * THIS SUBROUTINE COMPUTES THE GEHAN GENERALIZED WILCOXON STATISTIC   *
XC     * WITH CONDITIONAL PERMUTATION VARIANCE (HYPERGEOMETRIC VARIANCE)     *
XC     * FROM EQUATIONS (2.2) AND (2.3) IN LATTA, 'A MONTE-CARLO STUDY OF    *
XC     * SOME TWO-SAMPLE RANK TESTS WITH CENSORED DATA', 1981, JOURNAL OF    *
XC     * THE AMERICAN STATISTICAL ASSOCIATION, VOL 76, PP 713-719.           *
XC     *                                                                     *
XC     * INPUT                                                               *
XC     *      ID1(I) : INDICATOR OF CENSORSHIP OF XY(I)                      *
XC     *      ID2(I) : INDICATOR OF GROUP; 1 OR 2                            *
XC     *      XY(I)  : DATA POINTS (SORTED TO SMALLEST TO LARGEST)           *
XC     *      N1     : NUMBER OF DATA POINTS IN GROUP 1                      *
XC     *      N2     : NUMBER OF DATA POINTS IN GROUP 2                      *
XC     *      NCOMP  : TOTAL NUMBER OF DATA POINTS = N1 + N2                 *
XC     *                                                                     *
XC     * OUTPUT                                                              *
XC     *     TEST    : STANDARDIZED GEHAN STATISTIC                          *
XC     *     PROB    : PROBABILITY                                           *
XC     *                                                                     *
XC     * OTHERS                                                              *
XC     *      D1(I)  : THE NUMBER OF DETECTIONS OF GROUP 1 AT XY(I)          *
XC     *      D2(I)  : THE NUMBER OF DETECTIONS OF GROUP 2 AT XY(I)          *
XC     *      D(I)   : THE NUMBER OF DETECTIONS AT XY(I)                     *
XC     *      R1(I)  : RISK SET OF GROUP 1 AT XY(I)                          *
XC     *      R2(I)  : RISK SET OF GROUP 2 AT XY(I)                          *
XC     *      R(I)   : RISK SET AT XY(I)                                     *
XC     *      E1(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 1 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E2(I)  : THE NUMBER OF CENSORED POINTS OF GROUP 2 BETWEEN      *
XC     *                                                     XY(I) & XY(I+1) *
XC     *      E(I)   : THE NUMBER OF CENSORED POINTS BETWEEN XY(I) & XY(I+1) *
XC     *      SCORE  : SCORE OF THE DATA                                     *
XC     *      VAR    : VARIANCE OF THE DATA                                  *
X
X
X
X       IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
X
X       DIMENSION ID1(NTOT),ID2(NTOT),XY(NTOT)
X       DIMENSION D(NTOT),E(NTOT),R(NTOT)
X       DIMENSION D1(NTOT),E1(NTOT),R1(NTOT),D2(NTOT)
X       DIMENSION E2(NTOT),R2(NTOT)
X       COMMON /G/ NCOMP,N1,N2,NCEN,ISIGN,IFULL,LO
X
X
X       I = 1
X       L = 1
X       R1(L) = REAL(N1)
X       R2(L) = REAL(N2)
X       R(L)  = REAL(NCOMP)
X       ET1 = 0.0
X       ET2 = 0.0
X
XC
XC     *  IF THE SMALLEST VALUE IS CENSORED, THIS LOOP WILL GO THROUGH THE   *
XC     *  DATA UNTIL THE FIRST DETECTION IS REACHED.                         *
XC
X   10  IF(ID1(I) .NE. 0) THEN
X          IF(ID2(I) .EQ. 1) THEN
X             ET1 = ET1 + 1.0
X          ELSE
X             ET2 = ET2 + 1.0
X          ENDIF
X          I = I + 1
X          GOTO 10
X       ENDIF
XC
XC     *     START LOOP; THIS LOOP CONTINUES UNTIL THE COMPUTATION IS       *
XC     *     FINISHED.                                                      *
XC
X   20  D(L)  = 0.0
X       D1(L) = 0.0
X       D2(L) = 0.0
X       E(L)  = 0.0
X       E1(L) = 0.0
X       E2(L) = 0.0
X       TEMP  = XY(I)
XC
XC     * CHECK IF THE DATA POINT IS DETECTED OR NOT. IF DETECTED, CONTINUE. *
XC     * THEN CHECK WHICH GROUP THE DATA POINT BELONGS TO.                  *
XC     * COMPUTE THE SURVIVAL FUNCTION AND THE COEFFICIENT FOR THE          *
XC     * APPROPRIATE GROUP.                                                *
XC
X
X  30   IF(ID1(I) .EQ. 0) THEN
X          IF(ID2(I) .EQ. 1) THEN
X             D1(L) = D1(L) + 1.0
X          ELSE
X             D2(L) = D2(L) + 1.0
X          ENDIF
X
X          D(L) = D1(L) + D2(L)
X
XC
XC     * IF THE DATA POINT IS CENSORED, START CHECKING HOW MANY CENSORED    *
XC     * DATA POINTS THERE ARE BETWEEN XY(I) AND XY(I+1).                   *
XC
X       ELSE
X         IF(ID2(I) .EQ. 1) THEN
X            E1(L) = E1(L) + 1.0
X         ELSE
X            E2(L) = E2(L) + 1.0
X         ENDIF
X            E(L) = E1(L) + E2(L)
X         ENDIF
X
X         IF(I .LE. NCOMP) THEN
X           I = I + 1
XC
XC     * IF THE DATA POINT XY(I) IS TIED WITH PREVIOUS POINTS, GO BACK      *
XC     * TO ADDRESS 30, AND COUNT THE NUMBER OF TIED DATA POINTS.           *
XC     * ALSO, IF XY(I) IS NOT DETECTED GO BACK TO ADDRESS 30, AND COUNT    *
XC     * THE NUMBER OF THE CENSORED DATA POINTS                             *
XC
X           IF(TEMP .EQ. XY(I)) GOTO 30
X           IF(ID1(I) .NE. 0) GOTO 30
X
XC
XC     *            COMPUTE NEW RISK SETS FOR THE NEXT STEP.                *
XC
X           IF(L .EQ. 1) THEN
X             R1(L) = R1(L) - ET1
X             R2(L) = R2(L) - ET2
X             R(L)  = R1(L) + R2(L)
X          ELSE
X             R1(L) = R1(L-1) - D1(L-1) - E1(L-1)
X             R2(L) = R2(L-1) - D2(L-1) - E2(L-1)
X             R(L)  = R1(L) + R2(L)
X          ENDIF
X          L = L + 1
X          GOTO 20
X        ENDIF
XC
XC     *       COMPUTE THE SCORE AND VARIANCE                         *
XC
X
X        SCORE = 0.0
X        VAR   = 0.0
X        L1 = L - 1
X        DO 200 I = 1, L1
X
X           SCORE = SCORE+R(I)*(D2(I)-(R2(I)*D(I)/R(I)))
X
X           IF(R(I) .GT. 1.0) THEN
X              VAR = VAR + D(I)*(R(I)**2.0)*(R2(I)/R(I))*
X     +              (1.0-(R2(I)/R(I)))*((R(I)-D(I))/(R(I)-1.0))
X           ENDIF
X
XC           VAR = VAR+D(I)*((R(I)-REAL(I))**2)+E(I)*(REAL(I)**2)
X
X  200   CONTINUE
X
XC        VAR = VAR*REAL(N1*N2)/REAL(NCOMP*(NCOMP-1))
X
XCC
XC      *        NOW COMPUTE THE GEHAN STATISTIC                          *
XC
X        TEST = SCORE/DSQRT(VAR)
X        PROB = 1.0 - AGAUSS(TEST)
X 
X        RETURN
X        END
X
XC
XC      **********************************************************************
XC      *********************  SUBROUTINE XDATA ******************************
XC      **********************************************************************
XC
X       SUBROUTINE XDATA(X,XX,IND,IND2,IMUL,ICOL,NTOT,MVAR)
XC
XC      *  THIS SUBROUTINE CHANGES DATA FORMAT                               *
XC      *                                                                    *
XC      *  INPUT      X(I,J)   : VARIABLES                                   *
XC      *            IND(I,J)  : INDICATOR OF CENSORSHIP                     *
XC      *             IMUL     : NUMBER OF VARIABLES                         *
XC      *             ICOL     : COLUMN OF THE VARIABLE WHICH NEEDS TO BE    *
XC      *                        USED                                        *
XC      *             NTOT     : NUMBER OF DATA POINTS                       *
XC      *             MVAR     : DIMENSION SIZE                              *
XC      *                                                                    *
XC      *  OUTPUT     XX(I,J)  : VARIABLES                                   *
XC      *            IND2(I)   : INDICATOR OF CENSORSHIP                     *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X       DIMENSION X(MVAR,NTOT),XX(MVAR,NTOT),IND(MVAR,NTOT),IND2(NTOT)
XC
X       DO 100 I=1,NTOT
X
X          IF(IMUL.GT.1) THEN
XC
XC      *     THE PROBLEM WITH MORE THAN ONE INDEPENDENT VARIABLE            *
XC
X             DO 20 J=1,IMUL
X                XX(J,I)=X(J,I)
X                IND2(I)=IND(1,I)
X   20        CONTINUE
XC
X          ELSE
XC
XC      *     IF THE PROBLEM IS TWO DIMENSIONAL                              *
XC
X             XX(1,I)=X(ICOL,I)
X             IND2(I)=IND(1,I)
X          ENDIF
XC
X  100  CONTINUE
X       RETURN
X       END
X
XC
XC      **********************************************************************
XC      ********************* SUBROUTINE XVAR  *******************************
XC      **********************************************************************
XC
X       SUBROUTINE XVAR(IND,X,J,NTOT,ISIGN,ZU,ZC,IU,IC,ISAVE,
X     +                 ATIE,RISK,XT,ZTEMP,SWRK1,LTOT,MVAR,INDEX)
XC
XC      *       THIS SUBROUTINE DISTINGUISHES UNCENSORED AND CENSORED        *
XC      *       DATA IN THE X VARIABLE AND SORTS IT INTO ZU AND ZC. ALSO,    *
XC      *       IF THE DATA CONTAIN UPPER LIMITS, THE SIGN OF THE            *
XC      *       VALUES ARE CHANGED SO THAT THE PROGRAM FOR THE LOWER         *
XC      *       LIMITS CAN BE USED. ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
XC      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME          *
XC      *       LEARNING PUBLICATIONS (BELMONT:CA).                          *
XC      *                                                                    *
XC      *       INPUT      IND(J,I) : INDICATOR OF CENSORING                 *
XC      *                    X(J,I) : VARIABLE                               *
XC      *                   MVAR    : NUMBER OF THE VARIABLES( FOR DIM DEC.) *
XC      *                    J      : J-TH DATA SETS                         *
XC      *                   NTOT    : TOTAL NUMBER OF DATA POINTS            *
XC      *                                                                    *
XC      *       OUTPUT      ISIGN   : IF LOWER LIMIT, ISIGN = 1              *
XC      *                             IF UPPER LIMIT, ISIGN = -1             *
XC      *                   ZU(K)   : UNCENSORED DATA POINTS IN X(J,I)       *
XC      *                   ZC(K)   : CENSORED DATA POINTS IN X(J,I)         *
XC      *                    IU     : NUMBER OF UNCENSORED DATA POINTS       *
XC      *                    IC     : NUMBER OF CENSORED DATA POINTS         *
XC      *                   RISK(L) : RISK SET                               *
XC      *                  ATIE(L)  : NUMBER OF TIED DATA POINTS             *
XC      *                                                                    *
XC      *       OTHER                                                        *
XC      *                  ISAVE(I) : TEMPORARY SAVE OF ISIGN FOR EACH POINT *
XC      *                             ALSO USED FOR TEMPORARY CENSORSHIP     *
XC      *                             INDICATOR                              *
XC      *                   XT(I)   : = X(J,I)                               *
XC      *                                                                    *
XC      *        SUBROUTINES                                                 *
XC      *                   SORT1                                            *
XC
X       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
X
X       DIMENSION IND(MVAR,NTOT),X(MVAR,NTOT),ZU(NTOT)
X       DIMENSION ZC(NTOT),ISAVE(NTOT),ATIE(NTOT),RISK(NTOT)
X       DIMENSION XT(NTOT),ZTEMP(MVAR,NTOT),SWRK1(MVAR)
X       DIMENSION INDEX(NTOT)
XC
X       ISIGN=1     
X       IC=0
X       IU=0
XC
XC      *    FIND THE CENSORSHIP OF THE DATA SET. -1 FOR THE UPPER LIMIT    *
XC      *    AND 1 FOR THE LOWER LIMIT                                      *
XC
X       DO 100 I=1,NTOT
X          ISAVE(I) = 0
X          IF(IND(J,I) .EQ. 0) GOTO 100
X          ISIGN=IND(J,I)/IABS(IND(J,I))
X          ISAVE(I) = ISIGN
X  100  CONTINUE
X
XC      * CHECK WHETHER THE UPPER AND LOWER LIMITS ARE MIXED IN ONE        *
XC      * VARIABLE. IF SO, THE PROGRAM IS TERMINATED.                      *
XC
X       DO 110 I = 1, NTOT
X          IF(ISAVE(I) .EQ. 0) GOTO 110
X          IF(ISAVE(I) .NE. ISIGN) THEN
X             PRINT *
X             PRINT *,'YOU CANNOT HAVE BOTH UPPER AND LOWER LIMITS'
X             PRINT *,'IN ONE VARIABLE AT THE SAME TIME.'
X             PRINT *,'PLEASE CHECK THE DATA. THE PROGRAM IS TERMINATED.'
X             PRINT *
X             STOP
X          ENDIF
X  110  CONTINUE
XC
XC      *    IN CASE THE DATA HAS UPPER LIMITS IT IS MULTIPLIED BY ISIGN   *
XC      *    TO MAKE THE DATA HAVE LOWER LIMITS (RIGHT CENSORSHIP).        *
XC
X       DO 280 L = 1, NTOT
X          ATIE(L) = 0.0
X          XT(L) = REAL(ISIGN)*X(J,L)
X          ZTEMP(J,L) = 0.0
X          ISAVE(L) = IND(J,L)
X  280  CONTINUE
XC
XC     *     DATA POINTS ARE ARRANGED FROM SMALLEST TO LARGEST.             *
XC     *     DETECTED AND CENSORED DATA POINTS ARE SEPARATED.               *
XC     *     THEN RISK SETS AND TIED DATA POINTS ARE FOUND.                 *
XC
X
X       CALL SORT1(ISAVE,ZTEMP,XT,NTOT,J,INDEX,SWRK1,MVAR)
X
X       L = 1
X
X       DO 300 I=1,NTOT
X          K=IABS(ISAVE(I))
X          IF(K .EQ. 0) THEN 
X              IU=IU+1
X              ZU(IU)= XT(I)
X              IF(IU .NE. 1) THEN
X                 IF(ZU(IU) .EQ. ZU(IU-1)) THEN
X                    ATIE(L) = ATIE(L) + 1.0
X                    RISK(L) = REAL(NTOT - I)
X                 ELSE
X                    ATIE(L) = 1.0
X                    RISK(L) = REAL(NTOT - I)
X                    L = L + 1
X                 ENDIF
X              ELSE
X                 ATIE(L) = 1.0
X                 RISK(L) = REAL(NTOT - I)
X                 L = L + 1
X              ENDIF
X           ELSE
X              IC=IC+1
X              ZC(IC)= XT(I)
X           ENDIF
X  300   CONTINUE
X        LTOT = L - 1
X
X        RETURN
X        END