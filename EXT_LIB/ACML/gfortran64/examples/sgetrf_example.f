*     SGETRF Example Program Text
*
*     ACML version 1.0 Copyright AMD,NAG 2003
*
C     .. Parameters ..
      INTEGER          NMAX, LDA, NRHMAX, LDB
      PARAMETER        (NMAX=8,LDA=NMAX,NRHMAX=NMAX,LDB=NMAX)
      CHARACTER        TRANS
      PARAMETER        (TRANS='N')
C     .. Local Scalars ..
      INTEGER          I, INFO, J, N, NRHS
C     .. Local Arrays ..
      REAL             A(LDA,NMAX), B(LDB,NRHMAX)
      INTEGER          IPIV(NMAX)
C     .. External Subroutines ..
      EXTERNAL         SGETRF, SGETRS
C     .. Executable Statements ..
*
      WRITE (*,FMT=99998)
     *  'ACML example: solution of linear equations using SGETRF/SGETRS'
      WRITE (*,FMT=99998)
     *  '--------------------------------------------------------------'
      WRITE (*,FMT=*)
*
*     Initialize matrix A
      N = 4
      A(1,1) = 1.80
      A(1,2) = 2.88
      A(1,3) = 2.05
      A(1,4) = -0.89
      A(2,1) = 5.25
      A(2,2) = -2.95
      A(2,3) = -0.95
      A(2,4) = -3.80
      A(3,1) = 1.58
      A(3,2) = -2.69
      A(3,3) = -2.90
      A(3,4) = -1.04
      A(4,1) = -1.11
      A(4,2) = -0.66
      A(4,3) = -0.59
      A(4,4) = 0.80
*
*     Initialize right-hand-side matrix B
      NRHS = 2
      B(1,1) = 9.52
      B(1,2) = 18.47
      B(2,1) = 24.35
      B(2,2) = 2.25
      B(3,1) = 0.77
      B(3,2) = -13.28
      B(4,1) = -6.22
      B(4,2) = -6.21
*
      WRITE (*,FMT=99998) 'Matrix A:'
      DO 20 I = 1, N
         WRITE (*,FMT=99999) (A(I,J),J=1,N)
   20 CONTINUE
*
      WRITE (*,FMT=*)
      WRITE (*,FMT=99998) 'Right-hand-side matrix B:'
      DO 40 I = 1, N
         WRITE (*,FMT=99999) (B(I,J),J=1,NRHS)
   40 CONTINUE

*     Factorize A
      CALL SGETRF(N,N,A,LDA,IPIV,INFO)
*
      WRITE (*,FMT=*)
      IF (INFO.EQ.0) THEN
*
*        Compute solution
         CALL SGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
*
*        Print solution
         WRITE (*,FMT=99998) 'Solution matrix X of equations A*X = B:'
         DO 60 I = 1, N
            WRITE (*,FMT=99999) (B(I,J),J=1,NRHS)
   60    CONTINUE
      ELSE
         WRITE (*,FMT=99998) 'The factor U of matrix A is singular'
      END IF
*
*
99998 FORMAT (A)
99999 FORMAT (1X,1P,6E13.3)
      END
