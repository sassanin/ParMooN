
!   ACML version 4.1 Copyright AMD,NAG 2008

    PROGRAM MAIN
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: MSEED = 624, MSTATE = 633, NN = 11
!      .. Local Scalars ..
       INTEGER                         :: GENID, INFO, J, LSEED, LSTATE, N,    &
                                          SUBID
!      .. Local Arrays ..
       INTEGER                         :: ENS(NN), SEED(MSEED), STATE(MSTATE)
!      .. External Subroutines ..
       EXTERNAL                           DOTIME, DRANDINITIALIZE
!      .. Data Statements ..
       DATA ENS/2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 20000/
!      .. Executable Statements ..
       CONTINUE

!      Use the Mersenne twister generator as the base generator
       GENID = 3
       SUBID = 1

       LSTATE = MSTATE
       LSEED = 1
       SEED(1) = 1071958

!      Initialize the Mersenne twister base generator
       CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)

       DO J = 1, NN

          N = ENS(J)

!         To time different data, change ENS and NN, above,
!         accordingly.
          CALL DOTIME(N,STATE)

       END DO

    END PROGRAM MAIN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    SUBROUTINE DOTIME(N,STATE)
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: WP=KIND(0.0D0)
!      .. Scalar Arguments ..
       INTEGER                         :: N
!      .. Array Arguments ..
       INTEGER                         :: STATE(*)
!      .. Local Scalars ..
       REAL (KIND=WP)                  :: A, AVTIME, B, TX, TX0, TX1
       INTEGER                         :: I, IFAIL, INFO
!      .. Local Arrays ..
       REAL (KIND=WP), ALLOCATABLE     :: X(:)
!      .. External Functions ..
       REAL (KIND=WP), EXTERNAL        :: DSECND
!      .. External Subroutines ..
       EXTERNAL                           DRANDBETA
!      .. Intrinsic Functions ..
       INTRINSIC                          MIN, REAL
!      .. Executable Statements ..
       CONTINUE

!      Parameters for beta distribution
       A = 1.5E0_WP
       B = 2.0E0_WP

!      Allocate space
       ALLOCATE (X(N),STAT=IFAIL)

       IF (IFAIL==0) THEN

!         Using DRANDBETA
          TX = 1.0E30_WP
          DO I = 1, 1000
             TX0 = DSECND()
             CALL DRANDBETA(N,A,B,STATE,X,INFO)
             TX1 = DSECND()
             IF (TX1-TX0>0.0E0_WP) TX = MIN(TX1-TX0,TX)
          END DO

!         Print the results. Convert the times into milliseconds
          IF (INFO==0) THEN
             IF (TX<=0.0E0_WP) THEN
                AVTIME = 0.0E0_WP
             ELSE
                TX = REAL(1000,KIND=WP)*TX
                AVTIME = N/TX
             END IF
             WRITE (*,FMT=*) N, AVTIME
          ELSE
             WRITE (*,FMT=*) 'DRANDBETA failed with INFO = ', INFO
          END IF
       ELSE
          WRITE (*,*) 'Array allocation failed: attempted to allocate'
          WRITE (*,*) '   X(', N, ')'
       END IF

       IF (IFAIL==0) THEN
          DEALLOCATE (X)
       END IF

       RETURN
    END SUBROUTINE DOTIME
