*     CFFT1D Example Program Text
*
*     ACML version 3.0 Copyright AMD,NAG 2005
*
C
C     MODE = 0:   Initialise to use default plan for computing FFT 
C     MODE = 100: Initialise by timing different plans and selecting
C                  the plan that returns with the shortest execution
C                  time. (With this option the same plan may not be
C                  selected each time for the same value of N when
C                  two or more plans have a similar execution time.)
C
C     .. Parameters ..
      INTEGER          NMAX, IGEN, IDEF
      PARAMETER        (NMAX=10000, IGEN=100, IDEF=0)
C     .. Local Scalars ..
      INTEGER          I, INFO, N
C     .. Local Arrays ..
      COMPLEX          COMM(5*NMAX+100), X(NMAX), XX(NMAX)
C     .. External Subroutines ..
      EXTERNAL         CFFT1D
C     .. Executable Statements ..
*
      WRITE (*,FMT=*)
*
      WRITE (*,FMT=99999)
     *  ' ACML example: FFT of a complex sequence using CFFT1D'
      WRITE (*,FMT=99999)
     *  ' ----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     Use a model sequence of complex data
      N = 8
      DO I = 1, N
         X(I) = DCMPLX(REAL(I-1)/REAL(I),REAL(N-I)/REAL(N))
         XX(I) = X(I)
      END DO
*
*     Initialize communication array COMM using Generated plan.
      CALL CFFT1D(IGEN,N,X,COMM,INFO)
*
*     Compute a forward transform
      CALL CFFT1D(-1,N,X,COMM,INFO)
*
*     Initialize communication array COMM using Default plan.
      CALL CFFT1D(IDEF,N,XX,COMM,INFO)
*
*     Compute a forward transform
      CALL CFFT1D(-1,N,XX,COMM,INFO)

*     Print the two sets of results
      WRITE (*,99999) '   Generated Plan       Default Plan'
      WRITE (*,99999) '   Real      Imag      Real      Imag'
      DO I = 1, N
         WRITE (*,99998) X(I), XX(I)
      END DO
*
      WRITE (*,FMT=*)
*
99999 FORMAT (A)
99998 FORMAT (1X,'(',F7.4,',',F7.4,')',3X,'(',F7.4,',',F7.4,')')
      END
