*     ZFFT1DX Example Program Text
*
*     ACML version 2.5 Copyright AMD,NAG 2004
*
C     .. Parameters ..
      INTEGER          NMAX
      PARAMETER        (NMAX=20)
C     .. Local Scalars ..
      DOUBLE PRECISION SCALE
      INTEGER          INCX, INCY, INFO, J, N
      LOGICAL          INPL
C     .. Local Arrays ..
      COMPLEX *16      COMM(3*NMAX+100), X(NMAX), XX(NMAX), Y(NMAX)
C     .. External Subroutines ..
      EXTERNAL         ZFFT1DX
C     .. Intrinsic Functions ..
      INTRINSIC        DBLE
C     .. Executable Statements ..
*
      WRITE (*,FMT=99999)
     *  'ACML example: FFT of a complex sequence using ZFFT1DX'
      WRITE (*,FMT=99999)
     *  '-----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     The sequence of complex data
      SCALE = 1.0D0
      INPL = .FALSE.
      INCX = 2
      INCY = 1
      N = 7
      X(1) = (0.34907D0,-0.37168D0)
      X(3) = (0.54890D0,-0.35669D0)
      X(5) = (0.74776D0,-0.31174D0)
      X(7) = (0.94459D0,-0.23702D0)
      X(9) = (1.13850D0,-0.13274D0)
      X(11) = (1.32850D0,0.00074D0)
      X(13) = (1.51370D0,0.16298D0)
      DO 20 J = 1, N*INCX, INCX
         XX(J) = X(J)
   20 CONTINUE
*
*     Initialize work array
      CALL ZFFT1DX(0,SCALE,INPL,N,X,INCX,Y,INCY,COMM,INFO)
*
*     Compute a forward transform
      CALL ZFFT1DX(-1,SCALE,INPL,N,X,INCX,Y,INCY,COMM,INFO)
*
      WRITE (*,FMT=99999) 'Components of discrete Fourier transform:'
      WRITE (*,FMT=*)
      WRITE (*,FMT=99999) '          Real    Imag'
      DO 40 J = 1, N*INCY, INCY
         WRITE (*,FMT=99998) J, Y(J)
   40 CONTINUE
*
*     Compute the reverse transform
      SCALE = 1.0D0/DBLE(N)
      CALL ZFFT1DX(1,SCALE,INPL,N,Y,INCY,X,INCX,COMM,INFO)
*
      WRITE (*,FMT=*)
      WRITE (*,FMT=99999)
     *  'Original sequence as restored by inverse transform:'
      WRITE (*,FMT=*)
      WRITE (*,FMT=99999) '            Original            Restored'
      WRITE (*,FMT=99999) '          Real    Imag        Real    Imag'
      DO 60 J = 1, N*INCX, INCX
         WRITE (*,FMT=99998) J, XX(J), X(J)
   60 CONTINUE
*
99999 FORMAT (A)
99998 FORMAT (1X,I3,2(:3X,'(',F7.4,',',F7.4,')'))
      END
