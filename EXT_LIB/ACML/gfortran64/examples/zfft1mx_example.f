*     ZFFT1MX Example Program Text
*
*     ACML version 2.5 Copyright AMD,NAG 2004
*
C     .. Parameters ..
      INTEGER          NMAX, NSMAX
      PARAMETER        (NMAX=20,NSMAX=10)
C     .. Local Scalars ..
      DOUBLE PRECISION SCALE
      INTEGER          I, INCX1, INCX2, INCY1, INCY2, INFO, J, N, NSEQ
      LOGICAL          INPL
C     .. Local Arrays ..
      COMPLEX *16      COMM(3*NMAX+100), X(NMAX,NSMAX), XX(NMAX,NSMAX),
     *                 Y(NMAX-1,NSMAX)
C     .. External Subroutines ..
      EXTERNAL         ZFFT1MX
C     .. Intrinsic Functions ..
      INTRINSIC        COS, DBLE, DCMPLX, DSIN
C     .. Executable Statements ..
*
      WRITE (*,FMT=99998)
     *  'ACML example: FFT of a complex sequence using ZFFT1MX'
      WRITE (*,FMT=99998)
     *  '-----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     The sequence of complex data
      SCALE = 1.0D0
      INPL = .FALSE.
      INCX1 = 2
      INCX2 = NMAX
      INCY1 = 3
      INCY2 = NMAX - 1
      N = 5
      NSEQ = 2
      DO 40 I = 1, NSEQ
         DO 20 J = 1, N*INCX1, INCX1
            X(J,I) = DCMPLX(1.0D0+COS(DBLE(I+J)),DSIN(DBLE(I-J)))
            XX(J,I) = X(J,I)
   20    CONTINUE
   40 CONTINUE
*
*     Initialize work array
      CALL ZFFT1MX(0,SCALE,INPL,NSEQ,N,X,INCX1,INCX2,Y,INCY1,INCY2,COMM,
     *             INFO)
*
*     Compute a forward transform
      CALL ZFFT1MX(-1,SCALE,INPL,NSEQ,N,X,INCX1,INCX2,Y,INCY1,INCY2,
     *             COMM,INFO)
*
      WRITE (*,FMT=99998) 'Components of discrete Fourier transform:'
      WRITE (*,FMT=*)
      DO 80 I = 1, NSEQ
         WRITE (*,FMT=99997) 'Sequence Number ', I
         WRITE (*,FMT=99998) '          Real    Imag'
         DO 60 J = 1, N*INCY1, INCY1
            WRITE (*,FMT=99999) J, Y(J,I)
   60    CONTINUE
   80 CONTINUE
*
*     Compute the reverse transform
      SCALE = 1.0D0/DBLE(N)
      CALL ZFFT1MX(1,SCALE,INPL,NSEQ,N,Y,INCY1,INCY2,X,INCX1,INCX2,COMM,
     *             INFO)
*
      WRITE (*,FMT=*)
      WRITE (*,FMT=99998)
     *  'Original sequence as restored by inverse transform:'
      WRITE (*,FMT=*)
      WRITE (*,FMT=99998) '            Original            Restored'
      DO 120 I = 1, NSEQ
         WRITE (*,FMT=99997) 'Sequence Number ', I
         WRITE (*,FMT=99998)
     *     '          Real    Imag        Real    Imag'
         DO 100 J = 1, N*INCX1, INCX1
            WRITE (*,FMT=99999) J, XX(J,I), X(J,I)
  100    CONTINUE
  120 CONTINUE
*
99999 FORMAT (1X,I3,2(:3X,'(',F7.4,',',F7.4,')'))
99998 FORMAT (A)
99997 FORMAT (A,I3)
      END
