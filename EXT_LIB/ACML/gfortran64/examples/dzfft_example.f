*     DZFFT Example Program Text
*
*     ACML version 2.0 Copyright AMD,NAG 2004
*
*     .. Parameters ..
      INTEGER          NMAX
      PARAMETER        (NMAX=20)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N
*     .. Local Arrays ..
      DOUBLE PRECISION COMM(3*NMAX+100), X(NMAX), XX(NMAX)
*     .. External Subroutines ..
      EXTERNAL         DZFFT, ZDFFT
*     .. Executable Statements ..
*
      WRITE (*,99998)
     *  'ACML example: FFT of a real sequence using ZFFT1D'
      WRITE (*,99998)
     *  '--------------------------------------------------'
      WRITE (*,*)
*
*     The sequence of complex data
      N = 7
      X(1) = 0.34907D0
      X(2) = 0.54890D0
      X(3) = 0.74776D0
      X(4) = 0.94459D0
      X(5) = 1.13850D0
      X(6) = 1.32850D0
      X(7) = 1.51370D0
      DO 10 I = 1, N
         XX(I) = X(I)
  10  CONTINUE
*
*     Initialize communication array COMM
      CALL DZFFT(0,N,X,COMM,INFO)
*
*     Compute a real --> complex Hermitian transform.
      CALL DZFFT(1,N,X,COMM,INFO)
*
      WRITE (*,99998) 'Components of discrete Fourier transform:'
      WRITE (*,*)
      DO 20 J = 1, N
         WRITE (*,99999) J, X(J)
  20  CONTINUE
*
*     Conjugate the Vector X to simulate inverse transform.
      DO 30 J = N/2 + 2, N
         X(J) = -X(J)
   30 CONTINUE
*
*     Compute the complex Hermitian --> real transform.
      CALL ZDFFT(2,N,X,COMM,INFO)
*
      WRITE (*,*)
      WRITE (*,99998)
     *   'Original sequence as restored by inverse transform:'
      WRITE (*,*)
      WRITE (*,99998)'       Original  Restored'
      DO 40 J = 1, N
         WRITE (*,99999) J, XX(J), X(J)
  40  CONTINUE
*
99998 FORMAT (A)
99999 FORMAT (1X,I3,2(:3X,F7.4))
      END
