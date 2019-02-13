*     ZFFT3DX Example Program Text
*
*     ACML version 2.5 Copyright AMD,NAG 2004
*
C     .. Parameters ..
      INTEGER          LMAX, MMAX, NMAX, LW
      PARAMETER        (LMAX=10,MMAX=10,NMAX=10,
     *                 LW=9*MMAX+LMAX*MMAX*NMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION SCALE
      INTEGER          I, INFO, IX, IY, J, K, L, M, N
      LOGICAL          INPL, LTRANS
C     .. Local Arrays ..
      COMPLEX *16      COMM(LW), X(LMAX*MMAX*NMAX), XX(LMAX*MMAX*NMAX),
     *                 Y(LMAX*MMAX*NMAX)
C     .. External Subroutines ..
      EXTERNAL         ZFFT3DX
C     .. Intrinsic Functions ..
      INTRINSIC        COS, DBLE, DCMPLX, DSIN
C     .. Executable Statements ..
*
      WRITE (*,FMT=99998)
     *  'ACML example: FFT of a complex sequence using ZFFT3DX'
      WRITE (*,FMT=99998)
     *  '-----------------------------------------------------'
      WRITE (*,FMT=*)
*
*     The sequence of complex data
      L = 5
      M = 2
      N = 2
      SCALE = 1.0D0
      INPL = .FALSE.
      LTRANS = .FALSE.
      IX = 1
      DO 60 I = 1, N
         DO 40 J = 1, M
            DO 20 K = 1, L
               X(IX) = DCMPLX(1.0D0+COS(DBLE(IX)),DSIN(DBLE(IX)))
               XX(IX) = X(IX)
               IX = IX + 1
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
*
*     Initialize communication work array
      CALL ZFFT3DX(0,SCALE,LTRANS,INPL,L,M,N,X,Y,COMM,INFO)
*
*     Compute a forward transform
      CALL ZFFT3DX(-1,SCALE,LTRANS,INPL,L,M,N,X,Y,COMM,INFO)
*
      SCALE = 1.0D0/DBLE(L*M*N)
      WRITE (*,FMT=99998) 'Components of discrete Fourier transform:'
      WRITE (*,FMT=*)
      IF (LTRANS) THEN
         IY = 1
         DO 120 I = 1, N
            WRITE (*,FMT=*)
            WRITE (*,FMT=99997) 'Third dimension index = ', I
            DO 100 J = 1, M
               WRITE (*,FMT=99997) ' Column Number ', J
               WRITE (*,FMT=99998) '          Real    Imag'
               DO 80 K = 1, L
                  WRITE (*,FMT=99999) K, Y(IY)
                  IY = IY + 1
   80          CONTINUE
  100       CONTINUE
  120    CONTINUE
         CALL ZFFT3DX(1,SCALE,LTRANS,INPL,L,M,N,Y,X,COMM,INFO)
      ELSE
         DO 180 I = 1, N
            WRITE (*,FMT=*)
            WRITE (*,FMT=99997) 'Third dimension index = ', I
            IX = I
            DO 160 J = 1, M
               WRITE (*,FMT=99997) ' Column Number ', J
               WRITE (*,FMT=99998) '          Real    Imag'
               IY = IX
               DO 140 K = 1, L
                  WRITE (*,FMT=99999) K, Y(IY)
                  IY = IY + N*M
  140          CONTINUE
               IX = IX + N
  160       CONTINUE
  180    CONTINUE
         CALL ZFFT3DX(1,SCALE,LTRANS,INPL,N,M,L,Y,X,COMM,INFO)
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
      DO 240 I = 1, N
         WRITE (*,FMT=*)
         WRITE (*,FMT=99997) 'Third dimension index = ', I
         DO 220 J = 1, M
            WRITE (*,FMT=99997) 'Column Number ', J
            WRITE (*,FMT=99998)
     *        '          Real    Imag        Real    Imag'
            DO 200 K = 1, L
               WRITE (*,FMT=99999) K, XX(IX), X(IX)
               IX = IX + 1
  200       CONTINUE
  220    CONTINUE
  240 CONTINUE
*
99999 FORMAT (1X,I3,2(:3X,'(',F8.4,',',F8.4,')'))
99998 FORMAT (A)
99997 FORMAT (A,I3)
      END
