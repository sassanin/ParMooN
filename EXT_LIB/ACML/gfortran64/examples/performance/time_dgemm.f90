
!   ACML version 4.1 Copyright AMD,NAG 2008

    PROGRAM MAIN
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
       REAL (KIND=WP), PARAMETER       :: ALIM = -1.0E0_WP
       REAL (KIND=WP), PARAMETER       :: BLIM = 1.0E0_WP
       INTEGER, PARAMETER              :: NN = 11
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: ALPHA, BETA
       INTEGER                         :: J, K, M, N
       CHARACTER (1)                   :: TRANSA, TRANSB
!      .. Local Arrays ..
       INTEGER                         :: ENS(NN)
!      .. External Subroutines ..
       EXTERNAL                           DOTIME
!      .. Data Statements ..
       DATA ENS/100, 300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000/
!      .. Executable Statements ..
       CONTINUE

       DO J = 1, NN

          M = ENS(J)

!         To time different data, change ENS and NN above,
!         accordingly.
          N = M
          K = M

          TRANSA = 'N'
          TRANSB = 'N'
          ALPHA = 0.7E0_WP
          BETA = 1.3E0_WP
          CALL DOTIME(TRANSA,TRANSB,M,N,K,ALIM,BLIM,ALPHA,BETA)

       END DO

    END PROGRAM MAIN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    SUBROUTINE DOTIME(TRANSA,TRANSB,M,N,K,ALIM,BLIM,ALPHA,BETA)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, ALPHA, BETA, BLIM
       INTEGER                         :: K, M, N
       CHARACTER (1)                   :: TRANSA, TRANSB
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: DNFLOP, MFX, T0X, T1X, TX
       INTEGER                         :: IFAIL1, IFAIL2, IFAIL3, LDA, LDB, LDC
!      .. Local Arrays ..
       REAL (KIND=WP), ALLOCATABLE     :: A(:,:), B(:,:), C(:,:)
!      .. External Functions ..
       REAL (KIND=WP), EXTERNAL        :: DSECND
!      .. External Subroutines ..
       EXTERNAL                           DGEMM, RANDMAT
!      .. Intrinsic Functions ..
       INTRINSIC                          REAL
!      .. Executable Statements ..
       CONTINUE

!      Allocate space
       LDA = M
       LDB = K
       LDC = M
       ALLOCATE (A(LDA,K),STAT=IFAIL1)
       ALLOCATE (B(LDB,N),STAT=IFAIL2)
       ALLOCATE (C(LDC,N),STAT=IFAIL3)

       IF (IFAIL1==0 .AND. IFAIL2==0 .AND. IFAIL3==0) THEN

!         Number of MFLOPs in a DGEMM of size M-by-K-by-N:
          DNFLOP = 2.0E-6_WP*REAL(M,KIND=WP)*REAL(N,KIND=WP)*REAL(K,KIND=WP)

          IF (TRANSA=='N' .OR. TRANSA=='n') THEN
             CALL RANDMAT(ALIM,BLIM,M,K,A,LDA)
          ELSE
             CALL RANDMAT(ALIM,BLIM,K,M,A,LDA)
          END IF

          IF (TRANSB=='N' .OR. TRANSB=='n') THEN
             CALL RANDMAT(ALIM,BLIM,K,N,B,LDB)
          ELSE
             CALL RANDMAT(ALIM,BLIM,N,K,B,LDB)
          END IF

          CALL RANDMAT(ALIM,BLIM,M,N,C,LDC)

!         Using DGEMM
          T0X = DSECND()
          CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
          T1X = DSECND()
          TX = T1X - T0X

!         Print the results
          IF (TX<=0.0E0_WP) THEN
             MFX = 0.0E0_WP
          ELSE
             MFX = DNFLOP/TX
          END IF
          WRITE (*,FMT=*) N, MFX
       ELSE
          WRITE (*,*) 'Array allocation failed: attempted to allocate'
          WRITE (*,*) '   A(', LDA, ',', K, ') and'
          WRITE (*,*) '   B(', LDB, ',', N, ') and'
          WRITE (*,*) '   C(', LDC, ',', N, ')'
       END IF

       IF (IFAIL1==0) THEN
          DEALLOCATE (A)
       END IF
       IF (IFAIL2==0) THEN
          DEALLOCATE (B)
       END IF
       IF (IFAIL3==0) THEN
          DEALLOCATE (C)
       END IF

       RETURN
    END SUBROUTINE DOTIME

    SUBROUTINE RANDMAT(ALIM,BLIM,M,N,A,LDA)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
       INTEGER, PARAMETER              :: MSEED = 624, MSTATE = 633
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, BLIM
       INTEGER                         :: LDA, M, N
!      .. Array Arguments ..
       REAL (KIND=WP)                  :: A(LDA,N)
!      .. Local Scalars ..
       INTEGER                         :: GENID, I, INFO, LSEED, LSTATE, SUBID
!      .. Local Arrays ..
       INTEGER                         :: SEED(MSEED), STATE(MSTATE)
!      .. External Subroutines ..
       EXTERNAL                           DRANDINITIALIZE, DRANDUNIFORM
!      .. Executable Statements ..
       CONTINUE

!      Use the Mersenne twister generator as the base generator
       GENID = 3
       SUBID = 1

       LSTATE = MSTATE
       LSEED = 1
       SEED(1) = 1071958

!      Initialize the base generator
       CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)

       DO I = 1, N
!         Generate a matrix from a uniform U(ALIM,BLIM) distribution
          CALL DRANDUNIFORM(M,ALIM,BLIM,STATE,A(1,I),INFO)
       END DO

       RETURN
    END SUBROUTINE RANDMAT
