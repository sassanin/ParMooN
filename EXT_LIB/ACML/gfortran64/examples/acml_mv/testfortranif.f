C
C A program to test the Fortran interfaces
C supplied for the fast scalar and array libm routines.
C
      IMPLICIT         NONE
C     .. Parameters ..
      INTEGER          N1, N2
      PARAMETER        (N1=4,N2=8)
C     .. Local Scalars ..
      DOUBLE PRECISION A, B
      REAL             X, Y
C     .. Local Arrays ..
      DOUBLE PRECISION C(4), D(4), E(4)
      REAL             U(8), V(8), W(8)
C     .. External Functions ..
      DOUBLE PRECISION FASTCOS, FASTEXP, FASTLOG, FASTLOG10, FASTLOG2,
     *                 FASTPOW, FASTSIN
      REAL             FASTCOSF, FASTEXPF, FASTLOG10F, FASTLOG2F,
     *                 FASTLOGF, FASTPOWF, FASTSINF
      EXTERNAL         FASTCOS, FASTEXP, FASTLOG, FASTLOG10, FASTLOG2,
     *                 FASTPOW, FASTSIN, FASTCOSF, FASTEXPF, FASTLOG10F,
     *                 FASTLOG2F, FASTLOGF, FASTPOWF, FASTSINF
C     .. External Subroutines ..
      EXTERNAL         FASTSINCOS, FASTSINCOSF, VRDA_COS, VRDA_EXP,
     *                 VRDA_LOG, VRDA_LOG10, VRDA_LOG2, VRDA_SIN,
     *                 VRDA_SINCOS, VRSA_COSF, VRSA_EXPF, VRSA_LOG10F,
     *                 VRSA_LOG2F, VRSA_LOGF, VRSA_POWF, VRSA_POWXF,
     *                 VRSA_SINCOSF, VRSA_SINF
C     .. Executable Statements ..
      CONTINUE
C
C Test double precision routines
C
      WRITE (*,FMT=*) 'ACML vector math library example Fortran program'
      WRITE (*,FMT=*) '------------------------------------------------'
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'Testing double precision scalar routines'
C
      A = 3.0D0
      WRITE (*,FMT=99999) 'FASTEXP(', A, ') = ', FASTEXP(A)
      A = 10.0D0
      WRITE (*,FMT=99999) 'FASTLOG(', A, ') = ', FASTLOG(A)
      WRITE (*,FMT=99999) 'FASTLOG10(', A, ') = ', FASTLOG10(A)
      WRITE (*,FMT=99999) 'FASTLOG2(', A, ') = ', FASTLOG2(A)
      A = 3.0D0
      B = 2.5D0
      WRITE (*,FMT=99998) 'FASTPOW(', A, ',', B, ') = ', FASTPOW(A,B)
      A = 3.14159265358979D0
      WRITE (*,FMT=99999) 'FASTSIN(', A, ') = ', FASTSIN(A)
      WRITE (*,FMT=99999) 'FASTCOS(', A, ') = ', FASTCOS(A)
      CALL FASTSINCOS(A,C(1),D(1))
      WRITE (*,FMT=99998) 'FASTSINCOS(', A, ', S, C) = ', C(1), ',',
     *  D(1)
      A = 1000000.0D0
      WRITE (*,FMT=99997) 'FASTSIN(', A, ') = ', FASTSIN(A)
C
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'Testing double precision array routines'
      C(1) = 2.0D0
      C(2) = 6.0D0
      C(3) = 8.0D0
      C(4) = 9.0D0
      CALL VRDA_EXP(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_EXP(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '         ', D(1), D(2), D(3), D(4)
      CALL VRDA_LOG(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_LOG(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '         ', D(1), D(2), D(3), D(4)
      CALL VRDA_LOG10(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_LOG10(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '           ', D(1), D(2), D(3), D(4)
      CALL VRDA_LOG2(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_LOG2(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '           ', D(1), D(2), D(3), D(4)
      CALL VRDA_SIN(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_SIN(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '         ', D(1), D(2), D(3), D(4)
      CALL VRDA_COS(N1,C,D)
      WRITE (*,FMT=99996) 'VRDA_COS(', C(1), C(2), C(3), C(4), ') = '
      WRITE (*,FMT=99996) '         ', D(1), D(2), D(3), D(4)
      CALL VRDA_SINCOS(N1,C,D,E)
      WRITE (*,FMT=99996) 'VRDA_SINCOS(', C(1), C(2), C(3), C(4),
     *  ', S, C) = '
      WRITE (*,FMT=99996) '            ', D(1), D(2), D(3), D(4)
      WRITE (*,FMT=99996) '            ', E(1), E(2), E(3), E(4)
C
C Test single precision routines
C
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'Testing single precision scalar routines'
      X = 3.0
      WRITE (*,FMT=99999) 'FASTEXPF(', X, ') = ', FASTEXPF(X)
      X = 10.0
      WRITE (*,FMT=99999) 'FASTLOGF(', X, ') = ', FASTLOGF(X)
      WRITE (*,FMT=99999) 'FASTLOG10F(', X, ') = ', FASTLOG10F(X)
      WRITE (*,FMT=99999) 'FASTLOG2F(', X, ') = ', FASTLOG2F(X)
      X = 3.0
      Y = 2.5
      WRITE (*,FMT=99998) 'FASTPOWF(', X, ',', Y, ') = ', FASTPOWF(X,Y)
      X = 3.0
      WRITE (*,FMT=99999) 'FASTSINF(', X, ') = ', FASTSINF(X)
      X = 3.0
      WRITE (*,FMT=99999) 'FASTCOSF(', X, ') = ', FASTCOSF(X)
      X = 3.0
      CALL FASTSINCOSF(X,U(1),V(1))
      WRITE (*,FMT=99998) 'FASTSINCOSF(', X, ', S, C) = ', U(1), ',',
     *  V(1)
C
      WRITE (*,FMT=*)
      WRITE (*,FMT=*) 'Testing Single precision array routines'
      U(1) = 1.0
      U(2) = 2.0
      U(3) = 3.0
      U(4) = 4.0
      U(5) = 6.0
      U(6) = 7.0
      U(7) = 8.0
      U(8) = 9.0
      CALL VRSA_EXPF(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_EXPF(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99993) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
      CALL VRSA_LOGF(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_LOGF(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99994) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
      CALL VRSA_LOG10F(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_LOG10F(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99994) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
      CALL VRSA_LOG2F(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_LOG2F(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99994) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
C
      CALL VRSA_SINF(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_SINF(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99993) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
      CALL VRSA_COSF(N2,U,V)
      WRITE (*,FMT=99994) 'VRSA_COSF(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99994) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
C
      CALL VRSA_SINCOSF(N2,U,V,W)
      WRITE (*,FMT=99994) 'VRSA_SINCOSF(', U(1), U(2), U(3), U(4), U(5),
     *  U(6), U(7), U(8), ') = '
      WRITE (*,FMT=99993) '          ', V(1), V(2), V(3), V(4), V(5),
     *  V(6), V(7), V(8)
      WRITE (*,FMT=99993) '          ', W(1), W(2), W(3), W(4), W(5),
     *  W(6), W(7), W(8)
C
      CALL VRSA_POWF(N2,U,V,W)
      WRITE (*,FMT=99995) 'VRSA_POWF(', U(1), V(1), U(2), V(2), U(3),
     *  V(3), U(4), V(4), ','
      WRITE (*,FMT=99995) '          ', U(5), V(5), U(6), V(6), U(7),
     *  V(7), U(8), V(8), ') = '
      WRITE (*,FMT=99993) '          ', W(1), W(2), W(3), W(4), W(5),
     *  W(6), W(7), W(8)
C
      CALL VRSA_POWXF(N2,U,V(1),W)
      WRITE (*,FMT=99995) 'VRSA_POWXF(', U(1), V(1), U(2), V(1), U(3),
     *  V(1), U(4), V(1), ','
      WRITE (*,FMT=99995) '          ', U(5), V(1), U(6), V(1), U(7),
     *  V(1), U(8), V(1), ') = '
      WRITE (*,FMT=99993) '          ', W(1), W(2), W(3), W(4), W(5),
     *  W(6), W(7), W(8)
C
99999 FORMAT (1X,A13,F9.3,A,F9.3)
99998 FORMAT (1X,A13,F9.3,A,F9.3,A,F9.3)
99997 FORMAT (1X,A13,F13.4,A,F9.3)
99996 FORMAT (1X,A13,3(F9.3,','),F9.3,A)
99995 FORMAT (1X,A13,3(F7.3,'**',F6.3,','),F7.3,'**',F6.3,A)
99994 FORMAT (1X,A13,7(F7.3,','),F7.3,A)
99993 FORMAT (1X,A10,4(F15.3,','),/11X,3(F15.3,','),F15.3,A)
C
      END
