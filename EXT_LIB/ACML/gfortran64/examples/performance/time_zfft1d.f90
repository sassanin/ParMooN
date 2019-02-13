
!   ACML version 4.1 Copyright AMD,NAG 2008

    PROGRAM MAIN
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
       REAL (KIND=WP), PARAMETER       :: ALIM = 0.0E0_WP
       REAL (KIND=WP), PARAMETER       :: BLIM = 1.0E0_WP
       INTEGER, PARAMETER              :: NN = 19
!      .. Local Scalars ..
       INTEGER                         :: ISIGN, J, N, NRUNS
!      .. Local Arrays ..
       REAL (KIND=WP)                  :: ENS(NN)
!      .. External Subroutines ..
       EXTERNAL                           DOTIME
!      .. Data Statements ..
       DATA ENS/16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, &
          32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304/
!      .. Executable Statements ..
       CONTINUE

!      Time NRUNS forward transforms
       ISIGN = -1
       NRUNS = 4

       DO J = 1, NN

          N = ENS(J)

!         To time different data, change ENS and NN, above,
!         accordingly.
          CALL DOTIME(N,ISIGN,ALIM,BLIM,NRUNS)

       END DO

    END PROGRAM MAIN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    SUBROUTINE DOTIME(N,ISIGN,ALIM,BLIM,NRUNS)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
!      .. Scalar Arguments ..
       REAL (KIND=WP)                  :: ALIM, BLIM
       INTEGER                         :: ISIGN, N, NRUNS
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: DNFLOP, MFX, TX, TX0, TX1
       INTEGER                         :: IFAIL1, IFAIL2, IFAIL3, INFO, J
!      .. Local Arrays ..
       COMPLEX (KIND=WP), ALLOCATABLE  :: WORK(:), X(:)
       REAL (KIND=WP), ALLOCATABLE     :: IMA(:), REA(:)
!      .. External Functions ..
       REAL (KIND=WP), EXTERNAL        :: DSECND
!      .. External Subroutines ..
       EXTERNAL                           RANDZVEC, ZFFT1D
!      .. Intrinsic Functions ..
       INTRINSIC                          LOG, MIN, REAL
!      .. Executable Statements ..
       CONTINUE

!      Number of MFLOPs in a ZFFT1D of size N:
       DNFLOP = REAL(N,KIND=WP)*LOG(REAL(N,KIND=WP))*5.0E-6_WP/LOG(2.0E0_WP)

!      Allocate space
       ALLOCATE (X(N),STAT=IFAIL1)
       ALLOCATE (WORK(3*N+100),STAT=IFAIL2)
       ALLOCATE (IMA(N),REA(N),STAT=IFAIL3)

       IF (IFAIL1==0 .AND. IFAIL2==0 .AND. IFAIL3==0) THEN

!         Generate random vector X
          CALL RANDZVEC(ALIM,BLIM,N,X,IMA,REA)

!         Initialization call to ZFFT1D
          INFO = 0
          CALL ZFFT1D(0,N,X,WORK,INFO)

!         Using ZFFT1D. Record best time over NRUNS runs.
          TX = 1.0E30_WP
          DO J = 1, NRUNS
             TX0 = DSECND()
             CALL ZFFT1D(ISIGN,N,X,WORK,INFO)
             TX1 = DSECND()
             IF (TX1-TX0>0.0E0_WP) TX = MIN(TX,TX1-TX0)
          END DO

!         Print the results
          IF (INFO==0) THEN
             IF (TX<=0.0E0_WP) THEN
                MFX = 0.0E0_WP
             ELSE
                MFX = DNFLOP/TX
             END IF
             WRITE (*,FMT=*) N, MFX
          ELSE
             WRITE (*,FMT=*) 'ZFFT1D failed with INFO = ', INFO
          END IF

       ELSE
          WRITE (*,*) 'Array allocation failed:'
          WRITE (*,*) '  attempted to allocate X(', N, '),'
          WRITE (*,*) '  WORK(', 3*N + 100, '), IMA(', N, ')'
          WRITE (*,*) '  and REA(', N, ')'
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

    SUBROUTINE RANDZVEC(ALIM,BLIM,N,A,IMA,REA)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
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
       EXTERNAL                           DRANDINITIALIZE, DRANDUNIFORM
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
       CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)

!      Generate real and imaginary parts from a uniform U(ALIM,BLIM)
!      distribution
       CALL DRANDUNIFORM(N,ALIM,BLIM,STATE,REA,INFO)
       CALL DRANDUNIFORM(N,ALIM,BLIM,STATE,IMA,INFO)

       DO I = 1, N
          A(I) = CMPLX(REA(I),IMA(I),KIND=WP)
       END DO

       RETURN
    END SUBROUTINE RANDZVEC
