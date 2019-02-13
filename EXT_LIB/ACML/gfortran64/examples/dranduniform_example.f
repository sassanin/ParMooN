C     Simple DRANDINITIALIZE and DRANDUNIFORM example program
C     .. Parameters ..
      INTEGER          MSTATE, MSEED, MN
      PARAMETER        (MSTATE=20,MSEED=10,MN=12)
C     .. Local Scalars ..
      DOUBLE PRECISION A, B
      INTEGER          GENID, INFO, LSEED, LSTATE, N, SUBID
C     .. Local Arrays ..
      DOUBLE PRECISION X(MN)
      INTEGER          SEED(MSEED), STATE(MSTATE)
C     .. External Subroutines ..
      EXTERNAL         DRANDINITIALIZE, DRANDUNIFORM
C     .. Executable Statements ..
      CONTINUE

      WRITE (*,FMT=*)
      WRITE (*,FMT=99999)
     *  'ACML example: DRANDINITIALIZE and DRANDUNIFORM'
      WRITE (*,FMT=99999)
     *  '----------------------------------------------'
      WRITE (*,FMT=*)

C     Initialize the number of variates required
      N = MN

C     Use the basic generator as the base generator
      GENID = 1
      SUBID = 1

C     Populate the SEED, basic generator needs one seed, and a
C     STATE array of length 16
      LSTATE = 16
      LSEED = 1
      SEED(1) = 122421

C     Initialize the base generator
      CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)

C     Generate a sequence from a uniform U(0,1) distribution
      A = 0.0D0
      B = 1.0D0
      CALL DRANDUNIFORM(N,A,B,STATE,X,INFO)

C     Print the sequence
      WRITE (*,FMT=99999) 'Numbers from a uniform U(0,1) distribution:'
      WRITE (*,FMT=99998) X
C
99999 FORMAT (A)
99998 FORMAT (4F10.4)
      END
