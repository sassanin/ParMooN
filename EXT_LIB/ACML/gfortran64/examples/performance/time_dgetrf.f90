
!   ACML version 4.1 Copyright AMD,NAG 2008

    PROGRAM MAIN
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
       REAL (KIND=WP), PARAMETER       :: ALIM = 1.0E0_WP
       REAL (KIND=WP), PARAMETER       :: BLIM = 2.0E0_WP
       INTEGER, PARAMETER              :: NN = 11
!      .. Local Scalars ..
       INTEGER                         :: J, M, N
!      .. Local Arrays ..
       INTEGER                         :: ENS(NN)
!      .. External Subroutines ..
       EXTERNAL                           DOTIME
!      .. Data Statements ..
       DATA ENS/100, 400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000/
!      .. Executable Statements ..
       CONTINUE

       DO J = 1, NN
          M = ENS(J)

!         To time different data, change ENS, NN, and MAXN, above,
!         accordingly.
          N = M
          CALL DOTIME(M,N,ALIM,BLIM)

       END DO

    END PROGRAM MAIN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    SUBROUTINE DOTIME(M,N,ALIM,BLIM)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, BLIM
       INTEGER                         :: M, N
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: DNFLOP, MFX, T0X, T1X, TX
       INTEGER                         :: IFAIL1, IFAIL2, INFO, LDA
!      .. Local Arrays ..
       REAL (KIND=WP), ALLOCATABLE     :: A(:,:)
       INTEGER, ALLOCATABLE            :: IPIV(:)
!      .. External Functions ..
       REAL (KIND=WP), EXTERNAL        :: DSECND
!      .. External Subroutines ..
       EXTERNAL                           DGETRF, RANDMAT
!      .. Intrinsic Functions ..
       INTRINSIC                          MIN, REAL
!      .. Executable Statements ..
       CONTINUE

!      Allocate space
       LDA = M
       ALLOCATE (A(LDA,N),STAT=IFAIL1)
       ALLOCATE (IPIV(MIN(M,N)),STAT=IFAIL2)

       IF (IFAIL1==0 .AND. IFAIL2==0) THEN

!         Number of MFLOPs in a DGETRF of size M-by-N:
          IF (N<M) THEN
             DNFLOP = (REAL(M,KIND=WP)*(REAL(N,KIND=WP)**2)-(REAL(N, &
                KIND=WP)**3)/3.0E0_WP)/1.0E6_WP
          ELSE IF (N==M) THEN
             DNFLOP = ((2.0E0_WP/3.0E0_WP)*REAL(N,KIND=WP)**3)/1.0E6_WP
          ELSE
             DNFLOP = (REAL(N,KIND=WP)*(REAL(M,KIND=WP)**2)-(REAL(M, &
                KIND=WP)**3)/3.0E0_WP)/1.0E6_WP
          END IF

          CALL RANDMAT(ALIM,BLIM,M,N,A,LDA)

!         Using DGETRF
          T0X = DSECND()
          CALL DGETRF(M,N,A,LDA,IPIV,INFO)
          T1X = DSECND()
          TX = T1X - T0X

!         Print the results
          IF (INFO==0) THEN
             IF (TX<=0.0E0_WP) THEN
                MFX = 0.0E0_WP
             ELSE
                MFX = DNFLOP/TX
             END IF
             WRITE (*,FMT=*) N, MFX
          ELSE
             WRITE (*,FMT=*) 'DGETRF failed with INFO = ', INFO
          END IF
       ELSE
          WRITE (*,*) 'Array allocation failed: attempted to allocate'
          WRITE (*,*) '      A(', LDA, ',', N, ') and'
          WRITE (*,*) '   IPIV(', MIN(M,N), ')'
       END IF

       IF (IFAIL1==0) THEN
          DEALLOCATE (A)
       END IF
       IF (IFAIL2==0) THEN
          DEALLOCATE (IPIV)
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
