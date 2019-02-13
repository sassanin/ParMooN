*     DGETRF Example Program Text
*
*     ACML version 1.0 Copyright AMD,NAG 2003
*
*     .. Parameters ..
      INTEGER          LDA, NMAX, NRHMAX, LDB
      PARAMETER        (NMAX=8,LDA=NMAX,NRHMAX=NMAX,LDB=NMAX)
      CHARACTER        TRANS
      PARAMETER        (TRANS='N')
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N, NRHS
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), B(LDB,NRHMAX)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         DGETRF, DGETRS
*     .. Executable Statements ..
*
      WRITE (*,99998)
     *  'ACML example: solution of linear equations using DGETRF/DGETRS'
      WRITE (*,99998)
     *  '--------------------------------------------------------------'
      WRITE (*,*)
*
*     Initialize matrix A
      N = 4
      A(1,1) = 1.80D0
      A(1,2) = 2.88D0
      A(1,3) = 2.05D0
      A(1,4) = -0.89D0
      A(2,1) = 5.25D0
      A(2,2) = -2.95D0
      A(2,3) = -0.95D0
      A(2,4) = -3.80D0
      A(3,1) = 1.58D0
      A(3,2) = -2.69D0
      A(3,3) = -2.90D0
      A(3,4) = -1.04D0
      A(4,1) = -1.11D0
      A(4,2) = -0.66D0
      A(4,3) = -0.59D0
      A(4,4) = 0.80D0
*
*     Initialize right-hand-side matrix B
      NRHS = 2
      B(1,1) = 9.52D0
      B(1,2) = 18.47D0
      B(2,1) = 24.35D0
      B(2,2) = 2.25D0
      B(3,1) = 0.77D0
      B(3,2) = -13.28D0
      B(4,1) = -6.22D0
      B(4,2) = -6.21D0
*
      WRITE (*,99998) 'Matrix A:'
      DO 10 I = 1, N
         WRITE (*,99999) (A(I,J), J = 1, N)
  10  CONTINUE
*
      WRITE (*,*)
      WRITE (*,99998) 'Right-hand-side matrix B:'
      DO 20 I = 1, N
         WRITE (*,99999) (B(I,J), J = 1, NRHS)
  20  CONTINUE

*     Factorize A
      CALL DGETRF(N,N,A,LDA,IPIV,INFO)
*
      WRITE (*,*)
      IF (INFO.EQ.0) THEN
*
*        Compute solution
         CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
*
*        Print solution
         WRITE (*,99998) 'Solution matrix X of equations A*X = B:'
         DO 30 I = 1, N
            WRITE (*,99999) (B(I,J), J = 1, NRHS)
  30     CONTINUE
      ELSE
         WRITE (*,99998) 'The factor U of matrix A is singular'
      END IF
*
99998 FORMAT (A)
99999 FORMAT (1X,1P,6E13.3)
*
      END
