C---------------------------------------------------------------------72
C Numerical Recipes least squares fit:
C
C    Given a set of data points X(1:NDATA), Y(1:NDATA) with individual
C    standard deviations ERR(1:NDATA), use chi^2 minisation to fit for
C    some or all of the coefficients A(1:MA) of a function that
C    depends linearly on A, Y = Sum_i (A_i AFUNC_i(X)). The input
C    array IA(1:MA) indicates by nonzero entries those compontents of
C    a that should be fitted for, and by zero entries those components
C    that should be held fixed at their input values. The program
C    returns values for A(1:MA), chi^2 = CHISQ, and the covariance
C    matrix COVAR(1:MA,1:MA). (Parameters help fixed will return zero
C    covariances.) NPC is the physical dimension of COVAR(NPC,NPC) in
C    the calling routine. The user supplies a subroutine
C    FUNCS(x,AFUNC,MA) that returns the ma basis functions evaluated
C    at x = X in the array AFUNC.)
C
C---------------------------------------------------------------------72
      SUBROUTINE LFIT(X,Y,ERR,NDATA,A,IA,MA,COVAR,NPC,CHISQ,FUNCS,XX,YY)
C     ============================================================
C
      IMPLICIT NONE
C
      INTEGER*4 MA,IA(MA),NPC,NDATA,MMAX
      REAL*4    CHISQ,A(MA),COVAR(NPC,NPC),ERR(NDATA),X(NDATA),Y(NDATA)
      REAL*4    XX(NDATA),YY(NDATA)
      EXTERNAL  FUNCS
      PARAMETER (MMAX=50)
      INTEGER*4 I,J,K,L,M,MFIT
      REAL*4    SIG2I,SUM,WT,YM,AFUNC(MMAX),BETA(MMAX)
C
      MFIT=0
C
      DO J=1,MA
        IF(IA(J).NE.0)MFIT=MFIT+1
      ENDDO
C
      IF(MFIT.EQ.0)PAUSE 'LFIT: no parameters to be fitted!'
C
      DO J=1,MFIT
        DO K=1,MFIT
          COVAR(J,K)=0.
        ENDDO
        BETA(J)=0.
      ENDDO
C
      DO I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA,XX,YY)
        YM=Y(I)
        IF(MFIT.LT.MA)THEN
          DO J=1,MA
            IF(IA(J).EQ.0)YM=YM-A(J)*AFUNC(J)
          ENDDO
        ENDIF
        SIG2I=1.0/ERR(I)**2
        J=0
        DO L=1,MA
          IF(IA(L).NE.0)THEN
            J=J+1
            WT=AFUNC(L)*SIG2I
            K=0
            DO M=1,L
              IF(IA(M).NE.0)THEN
                K=K+1
                COVAR(J,K)=COVAR(J,K)+WT*AFUNC(M)
              ENDIF
            ENDDO
            BETA(J)=BETA(J)+YM*WT
          ENDIF
        ENDDO
      ENDDO
C
      DO J=2,MFIT
        DO K=1,J-1
          COVAR(K,J)=COVAR(J,K)
        ENDDO
      ENDDO
C
      CALL GAUSSJ(COVAR,MFIT,NPC,BETA,1,1)
      J=0
C
      DO L=1,MA
        IF(IA(L).NE.0)THEN
          J=J+1
          A(L)=BETA(J)
        ENDIF
      ENDDO
C
      CHISQ=0.0
C
      DO I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA,XX,YY)
        SUM=0.0
        DO J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
        ENDDO
        CHISQ=CHISQ+((Y(I)-SUM)/ERR(I))**2
      ENDDO
      CALL COVSRT(COVAR,NPC,MA,IA,MFIT)
C
      RETURN
      END
C---------------------------------------------------------------------72
      SUBROUTINE COVSRT(COVAR,NPC,MA,IA,MFIT)
C     =======================================
C
      IMPLICIT NONE
C
      INTEGER*4 MA,MFIT,NPC,IA(MA)
      REAL*4    COVAR(NPC,NPC)
      INTEGER*4 I,J,K
      REAL*4    SWAP
C
      DO I=MFIT+1,MA
        DO J=1,I
          COVAR(I,J)=0.0
          COVAR(J,I)=0.0
        ENDDO
      ENDDO
C
      K=MFIT
C
      DO J=MA,1,-1
        IF(IA(J).NE.0)THEN
          DO I=1,MA
            SWAP=COVAR(I,K)
            COVAR(I,K)=COVAR(I,J)
            COVAR(I,J)=SWAP
          ENDDO
          DO I=1,MA
            SWAP=COVAR(K,I)
            COVAR(K,I)=COVAR(J,I)
            COVAR(J,I)=SWAP
          ENDDO
          K=K-1
        ENDIF
      ENDDO
C
      RETURN
      END
C---------------------------------------------------------------------72
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
C     ================================
C
      IMPLICIT NONE
C
      INTEGER*4 M,MP,N,NP,NMAX
      REAL*4 A(NP,NP),B(NP,MP)
      PARAMETER (NMAX=50)
      INTEGER*4 I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),IPIV(NMAX)
      REAL*4 BIG,DUM,PIVINV
C
      DO J=1,N
        IPIV(J)=0
      ENDDO
C
      DO I=1,N
        BIG=0.0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
              IF(IPIV(K).EQ.0)THEN
                IF(ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IPIV(ICOL)=IPIV(ICOL)+1
C
        IF(IROW.NE.ICOL)THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          ENDDO
          DO L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
          ENDDO
        ENDIF
C
        INDXR(I)=IROW
        INDXC(I)=ICOL
C
        IF(A(ICOL,ICOL).EQ.0.)PAUSE 'Singular matrix in GAUSSJ'
        PIVINV=1.0/A(ICOL,ICOL)
        A(ICOL,ICOL)=1.0
C
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
        ENDDO
        DO L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
        ENDDO
C
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.0
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            ENDDO
            DO L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
      DO L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C---------------------------------------------------------------------72
      SUBROUTINE TWIST(X,AFIT,NAFIT,XX,YY)
C     ==============================
C X is a label, to lookup actual x,y (i.e. XX,YY) values from COMMON
C
      INTEGER*4 NAFIT
      REAL*4 X,AFIT(NAFIT)
      REAL*4 XX(1000),YY(1000)
C
C      COMMON /FITS/XX,YY
C
      AFIT(1)=1.0
      AFIT(2)=XX(INT(X))
      AFIT(3)=YY(INT(X))
C
      RETURN
      END
