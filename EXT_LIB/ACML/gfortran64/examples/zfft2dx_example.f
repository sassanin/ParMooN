*     ZFFT2DX Example Program Text
*
*     ACML version 2.5 Copyright AMD,NAG 2004
*
C     .. Parameters ..
      INTEGER          MMAX, NMAX, LW
      PARAMETER        (MMAX=20,NMAX=10,LW=3*MMAX+3*NMAX+MMAX*NMAX+100)
C     .. Local Scalars ..
      DOUBLE PRECISION SCALE
      INTEGER          I, INCX1, INCX2, INCY1, INCY2, INFO, IX, IY, J,
     *                 M, N
      LOGICAL          INPL, LTRANS
C     .. Local Arrays ..
      COMPLEX *16      COMM(LW), X(MMAX*NMAX), XX(MMAX*NMAX),
     *                 Y(MMAX*NMAX)
C     .. External Subroutines ..
      EXTERNAL         ZFFT2DX
C     .. Intrinsic Functions ..
      INTRINSIC        COS, DBLE, DCMPLX, DSIN
C     .. Executable Statements ..
*
      WRITE (*,FMT=99998)
     *  'ACML example: FFT of a complex sequence using ZFFT2DX'
      WRITE (*,FMT=99998)
     *  '-----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     The sequence of complex data
      M = 5
      N = 2
      SCALE = 1.0D0
      INPL = .FALSE.
      LTRANS = .FALSE.
      INCX1 = 1
      INCX2 = M
      INCY1 = 1
      INCY2 = N
      IX = 1
      DO 40 J = 1, N
         DO 20 I = 1, M*INCX1, INCX1
            X(IX) = DCMPLX(1.0D0+COS(DBLE(I+2*J)),DSIN(DBLE(3*I-J)))
            XX(IX) = X(IX)
            IX = IX + 1
   20    CONTINUE
   40 CONTINUE
*
*     Initialize work array
      CALL ZFFT2DX(0,SCALE,LTRANS,INPL,M,N,X,INCX1,INCX2,Y,INCY1,INCY2,
     *             COMM,INFO)
*
*     Compute a forward transform
      CALL ZFFT2DX(-1,SCALE,LTRANS,INPL,M,N,X,INCX1,INCX2,Y,INCY1,INCY2,
     *             COMM,INFO)
*
      SCALE = 1.0D0/DBLE(M*N)
      WRITE (*,FMT=99998) 'Components of discrete Fourier transform:'
      WRITE (*,FMT=*)
      IF (LTRANS) THEN
         IY = 1
         DO 80 I = 1, N
            WRITE (*,FMT=99997) 'Column Number ', I
            WRITE (*,FMT=99998) '          Real    Imag'
            DO 60 J = 1, M
               WRITE (*,FMT=99999) J, Y(IY)
               IY = IY + 1
   60       CONTINUE
   80    CONTINUE
         CALL ZFFT2DX(1,SCALE,LTRANS,INPL,M,N,Y,INCY1,INCY2,X,INCX1,
     *                INCX2,COMM,INFO)
      ELSE
         DO 120 I = 1, N
            WRITE (*,FMT=99997) 'Column Number ', I
            WRITE (*,FMT=99998) '          Real    Imag'
            IY = I
            DO 100 J = 1, M
               WRITE (*,FMT=99999) J, Y(IY)
               IY = IY + N
  100       CONTINUE
  120    CONTINUE
         INCY2 = N
         CALL ZFFT2DX(1,SCALE,LTRANS,INPL,N,M,Y,INCY1,INCY2,X,INCX1,
     *                INCX2,COMM,INFO)
      END IF
*
*     Compute the reverse transform
*
      WRITE (*,FMT=*)
      WRITE (*,FMT=99998)
     *  'Original sequence as restored by inverse transform:'
      WRITE (*,FMT=*)
      WRITE (*,FMT=99998) '            Original            Restored'
      IX = 1
      DO 160 J = 1, N
         WRITE (*,FMT=99997) 'Column Number ', J
         WRITE (*,FMT=99998)
     *     '          Real    Imag        Real    Imag'
         DO 140 I = 1, M*INCX1, INCX1
            WRITE (*,FMT=99999) J, XX(IX), X(IX)
            IX = IX + 1
  140    CONTINUE
  160 CONTINUE
*
99999 FORMAT (1X,I3,2(:3X,'(',F7.4,',',F7.4,')'))
99998 FORMAT (A)
99997 FORMAT (A,I3)
      END
