*     ZFFT1D Example Program Text
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
      PARAMETER        (NMAX=10000,IGEN=100,IDEF=0)
C     .. Local Scalars ..
      INTEGER          I, INFO, N
C     .. Local Arrays ..
      COMPLEX *16      COMM(3*NMAX+100), X(NMAX), XX(NMAX)
C     .. External Subroutines ..
      EXTERNAL         ZFFT1D
C     .. Executable Statements ..
*
      WRITE (*,FMT=*)
*
      WRITE (*,FMT=99999)
     *  ' ACML example: FFT of a complex sequence using ZFFT1D'
      WRITE (*,FMT=99999)
     *  ' ----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     Use a model sequence of complex data
      N = 8
      DO I = 1, N
         X(I) = DCMPLX(DBLE(I-1)/DBLE(I),DBLE(N-I)/DBLE(N))
         XX(I) = X(I)
      END DO
*
*     Initialize communication array COMM using Generated plan.
      CALL ZFFT1D(IGEN,N,X,COMM,INFO)
*
*     Compute a forward transform
      CALL ZFFT1D(-1,N,X,COMM,INFO)
*
*     Initialize communication array COMM using Default plan.
      CALL ZFFT1D(IDEF,N,XX,COMM,INFO)
*
*     Compute a forward transform
      CALL ZFFT1D(-1,N,XX,COMM,INFO)

*     Print the two sets of results
      WRITE (*,FMT=99999) '   Generated Plan       Default Plan'
      WRITE (*,FMT=99999) '   Real      Imag      Real      Imag'
      DO I = 1, N
         WRITE (*,FMT=99998) X(I), XX(I)
      END DO
*
      WRITE (*,FMT=*)
*
99999 FORMAT (A)
99998 FORMAT (1X,'(',F7.4,',',F7.4,')',3X,'(',F7.4,',',F7.4,')')
      END
