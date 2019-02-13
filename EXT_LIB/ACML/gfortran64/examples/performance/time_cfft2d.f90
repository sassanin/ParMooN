
!   ACML version 4.1 Copyright AMD,NAG 2008

    PROGRAM MAIN
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0)
       REAL (KIND=WP), PARAMETER       :: ALIM = 0.0_WP, BLIM = 1.0_WP
       INTEGER, PARAMETER              :: NN = 6
!      .. Local Scalars ..
       INTEGER                         :: ISIGN, J, M, N, NRUNS
!      .. Local Arrays ..
       REAL (KIND=WP)                  :: ENS(NN)
!      .. External Subroutines ..
       EXTERNAL                           DOTIME
!      .. Data Statements ..
       DATA ENS/128, 256, 512, 1024, 2048, 4096/
!      .. Executable Statements ..
       CONTINUE

!      Time NRUNS forward transforms
       ISIGN = -1
       NRUNS = 4

       DO J = 1, NN
          N = ENS(J)

!         To time different data, change ENS and NN above,
!         accordingly.
          M = N
          CALL DOTIME(M,N,ISIGN,ALIM,BLIM,NRUNS)
       END DO

    END PROGRAM MAIN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    SUBROUTINE DOTIME(M,N,ISIGN,ALIM,BLIM,NRUNS)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0)
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, BLIM
       INTEGER                         :: ISIGN, M, N, NRUNS
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: DNFLOP, MFX, TX, TX0, TX1
       INTEGER                         :: IFAIL1, IFAIL2, IFAIL3, INFO, J
!      .. Local Arrays ..
       COMPLEX (KIND=WP), ALLOCATABLE  :: WORK(:), X(:)
       REAL (KIND=WP), ALLOCATABLE     :: IMA(:), REA(:)
!      .. External Functions ..
       REAL (KIND=WP), EXTERNAL        :: SECOND
!      .. External Subroutines ..
       EXTERNAL                           CFFT2D, RANDCVEC
!      .. Intrinsic Functions ..
       INTRINSIC                          LOG, MIN, REAL
!      .. Executable Statements ..
       CONTINUE

!      Number of MFLOPs in a CFFT2D of size M-by-N:
       DNFLOP = REAL(M*N,KIND=WP)*LOG(REAL(N*M,KIND=WP))*5.0E-6_WP/LOG(2.0_WP)

!      Allocate space
       ALLOCATE (X(M*N),STAT=IFAIL1)
       ALLOCATE (WORK(M*N+5*(M+N)),STAT=IFAIL2)
       ALLOCATE (IMA(M*N),REA(M*N),STAT=IFAIL3)

       IF (IFAIL1==0 .AND. IFAIL2==0 .AND. IFAIL3==0) THEN

!         Generate random vector X
          CALL RANDCVEC(ALIM,BLIM,M*N,X,IMA,REA)

!         Using CFFT2D. Record best time over NRUNS runs.
          TX = 1.0E30_WP
          DO J = 1, NRUNS
             TX0 = SECOND()
             CALL CFFT2D(ISIGN,M,N,X,WORK,INFO)
             TX1 = SECOND()
             IF (TX1-TX0>0.0E0_WP) TX = MIN(TX,TX1-TX0)
          END DO

!         Print the results
          IF (INFO==0) THEN
             IF (TX<=0.0_WP) THEN
                MFX = 0.0_WP
             ELSE
                MFX = DNFLOP/TX
             END IF
             WRITE (*,FMT=*) N, MFX
          ELSE
             WRITE (*,FMT=*) 'CFFT2D failed with INFO = ', INFO
          END IF

       ELSE
          WRITE (*,*) 'Array allocation failed:'
          WRITE (*,*) '  attempted to allocate X(', M*N, '),'
          WRITE (*,*) '  WORK(', M*N + 5*(M+N), '), IMA(', M*N, ')'
          WRITE (*,*) '  and REA(', M*N, ')'
       END IF

       IF (IFAIL1==0) THEN
          DEALLOCATE (X)
       END IF
       IF (IFAIL2==0) THEN
          DEALLOCATE (WORK)
       END IF
       IF (IFAIL3==0) THEN
          DEALLOCATE (IMA,REA)
       END IF

       RETURN
    END SUBROUTINE DOTIME

    SUBROUTINE RANDCVEC(ALIM,BLIM,N,A,IMA,REA)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0)
       INTEGER, PARAMETER              :: MSEED = 624, MSTATE = 633
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, BLIM
       INTEGER                         :: N
!      .. Array Arguments ..
       COMPLEX (KIND=WP)               :: A(N)
       REAL (KIND=WP)                  :: IMA(N), REA(N)
!      .. Local Scalars ..
       INTEGER                         :: GENID, I, INFO, LSEED, LSTATE, SUBID
!      .. Local Arrays ..
       INTEGER                         :: SEED(MSEED), STATE(MSTATE)
!      .. External Subroutines ..
       EXTERNAL                           SRANDINITIALIZE, SRANDUNIFORM
!      .. Intrinsic Functions ..
       INTRINSIC                          CMPLX
!      .. Executable Statements ..
       CONTINUE

!      Use the Mersenne twister generator as the base generator
       GENID = 3
       SUBID = 1

       LSTATE = MSTATE
       LSEED = 1
       SEED(1) = 1071958

!      Initialize the base generator
       CALL SRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)

!      Generate real and imaginary parts from a uniform U(ALIM,BLIM)
!      distribution
       CALL SRANDUNIFORM(N,ALIM,BLIM,STATE,REA,INFO)
       CALL SRANDUNIFORM(N,ALIM,BLIM,STATE,IMA,INFO)

       DO I = 1, N
          A(I) = CMPLX(REA(I),IMA(I))
       END DO

       RETURN
    END SUBROUTINE RANDCVEC
