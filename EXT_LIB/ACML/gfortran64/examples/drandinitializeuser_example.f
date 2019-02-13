C     DRANDINITIALIZEUSER Example Program Text
C
C     ACML version 3.0 Copyright AMD,NAG 2005
C
C     Example program for DRANDINITIALIZEUSER
C     .. Parameters ..
      INTEGER          MSTATE, MSEED, MN
      PARAMETER        (MSTATE=28,MSEED=10,MN=12)
C     .. Local Scalars ..
      DOUBLE PRECISION A, B
      INTEGER          INFO, LSEED, LSTATE, N
C     .. Local Arrays ..
      DOUBLE PRECISION X(MN)
      INTEGER          SEED(MSEED), STATE(MSTATE)
C     .. External Subroutines ..
      EXTERNAL         DRANDINITIALIZEUSER, DRANDUNIFORM, WHGEN, WHINI
C     .. Executable Statements ..
      CONTINUE

      WRITE (*,FMT=*)
      WRITE (*,FMT=99999) 'ACML example: DRANDINITIALIZEUSER'
      WRITE (*,FMT=99999) '---------------------------------'
      WRITE (*,FMT=*)
      WRITE (*,FMT=99999)
     *  'Use of a user-supplied base generator for random numbers'
      WRITE (*,FMT=*)

C     Initialize the number of variates required
      N = MN

C     Populate the SEED, and give the length of the STATE array
C     Note that the STATE array must be declared 3 elements larger
C     than required by the user-supplied WHINI routine. The extra
C     three elements are required by DRANDINITIALIZEUSER itself
C     rather than by WHINI.
      LSTATE = 28
      LSEED = 4
      SEED(1) = 10
      SEED(2) = 7
      SEED(3) = 19
      SEED(4) = 58

C     Initialize the base generator.
      CALL DRANDINITIALIZEUSER(WHINI,WHGEN,1,1,SEED,LSEED,STATE,LSTATE,
     *                         INFO)

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

      SUBROUTINE WHINI(GENID,SUBID,SEED,LSEED,STATE,LSTATE,INFO)
C
C     WHINI - Initialisation routine for the user-supplied
C             base generator.
C
C     Reference: the base generator implemented here is described in:
C           "Generating good pseudo-random numbers", B.A.Wichmann
C           and I.D.Hill (preprint version 0.92, April 2005)
C
C     Input:
C     GENID  - Generator ID (irrelevant in this example)
C     SUBID  - Sub generator ID (irrelevant in this example)
C     SEED   - Array of seeds, SEED(i) > 0 for i = 1 to LSEED
C     LSEED  - Length of array SEED. In this example, LSEED must be 4.
C     LSTATE - Length of array STATE, LSTATE >= 25
C
C     Output:
C     STATE  - Array holding information about the generator:
C                    1 = Amount of STATE used = 25
C                    2 = 'Magic' number
C                    3 = GENID
C                    4 = SUBID
C              5 to  8 = Current state of the generator
C              9 to 24 = Parameters of constituent parts of the generator.
C                   25 = 'Magic' number again
C     INFO   - Failure code:
C              0        = Everything OK
C              -1 to -6 = Error in parameter ABS(INFO)
C     .. Parameters ..
      INTEGER          ESEED, ESTATE
      PARAMETER        (ESEED=4,ESTATE=25)
      INTEGER          MAGIC
      PARAMETER        (MAGIC=54321)
C     .. Scalar Arguments ..
      INTEGER          GENID, INFO, LSEED, LSTATE, SUBID
C     .. Array Arguments ..
      INTEGER          SEED(*), STATE(*)
C     .. Local Scalars ..
      INTEGER          I
C     .. Executable Statements ..
      CONTINUE
C
C     Basic parameter checking
      IF (LSEED.LT.ESEED) THEN
         INFO = -4
         RETURN
      END IF
      IF (LSTATE.LT.ESTATE) THEN
         INFO = -6
         RETURN
      END IF
C     End of basic parameter checking
C
C     Check and copy the supplied seed into the STATE array
      DO 20 I = 1, ESEED
         IF (SEED(I).LT.1) THEN
            INFO = -3
            RETURN
         ELSE
            STATE(4+I) = SEED(I)
         END IF
   20 CONTINUE
C
C     Populate the rest of the STATE array
      STATE(1) = ESTATE
      STATE(2) = MAGIC
      STATE(3) = GENID
      STATE(4) = SUBID
C
C     These define the parameters of the four constituent generators
      STATE(9) = 11600
      STATE(10) = 185127
      STATE(11) = 10379
      STATE(12) = 2147483579
      STATE(13) = 47003
      STATE(14) = 45688
      STATE(15) = 10479
      STATE(16) = 2147483543
      STATE(17) = 23000
      STATE(18) = 93368
      STATE(19) = 19423
      STATE(20) = 2147483423
      STATE(21) = 33000
      STATE(22) = 65075
      STATE(23) = 8123
      STATE(24) = 2147483123
C
C     Finally, repeat the magic number
      STATE(25) = MAGIC
C
      INFO = 0
      RETURN
      END

      SUBROUTINE WHGEN(N,STATE,X,INFO)
C
C     WHGEN: Generate a vector of pseudo random numbers
C            using the WH generator.
C
C     See WHINI for reference.
C
C     Input:
C     N  - Number of values being generated, N >= 0
C
C     Input / Output:
C     STATE - On Entry: The current state of the generator
C             On Exit : The state of the generator after
C                       generating the N values
C
C     Output:
C     X     - Vector holding the generated values, X(N)
C     INFO  - Error code:
C             0        = Everything OK
C             -1 to -2 = Error in parameter ABS(INFO)
C     .. Parameters ..
      INTEGER          MAGIC, ESTATE
      PARAMETER        (MAGIC=54321,ESTATE=25)
C     .. Scalar Arguments ..
      INTEGER          INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
      INTEGER          STATE(ESTATE)
C     .. Local Scalars ..
      DOUBLE PRECISION W
      INTEGER          I, IMAGT, IMAGX, IMAGY, IMAGZ, IMODT, IMODX,
     *                 IMODY, IMODZ, IMULT, IMULX, IMULY, IMULZ, ISUBT,
     *                 ISUBX, ISUBY, ISUBZ, IT, IX, IY, IZ
C     .. Intrinsic Functions ..
      INTRINSIC        DBLE, INT, MOD
C     .. Executable Statements ..
      CONTINUE
C
C     Parameter checks
      IF (N.LT.0) THEN
         INFO = -1
         RETURN
      END IF
C     Check the magic number indicating that the initialization
C     routine has previously been called
      IF (STATE(2).NE.MAGIC .OR. STATE(STATE(1)).NE.MAGIC) THEN
         INFO = -2
         RETURN
      END IF
C     End of parameter checks
C
C     Extract current state
      IX = STATE(5)
      IY = STATE(6)
      IZ = STATE(7)
      IT = STATE(8)
      IMULX = STATE(9)
      IMODX = STATE(10)
      ISUBX = STATE(11)
      IMAGX = STATE(12)
      IMULY = STATE(13)
      IMODY = STATE(14)
      ISUBY = STATE(15)
      IMAGY = STATE(16)
      IMULZ = STATE(17)
      IMODZ = STATE(18)
      ISUBZ = STATE(19)
      IMAGZ = STATE(20)
      IMULT = STATE(21)
      IMODT = STATE(22)
      ISUBT = STATE(23)
      IMAGT = STATE(24)
C
C     Generate the N values
      DO 20 I = 1, N
C
         IX = IMULX*(MOD(IX,IMODX)) - ISUBX*(IX/IMODX)
         IY = IMULY*(MOD(IY,IMODY)) - ISUBY*(IY/IMODY)
         IZ = IMULZ*(MOD(IZ,IMODZ)) - ISUBZ*(IZ/IMODZ)
         IT = IMULT*(MOD(IT,IMODT)) - ISUBT*(IT/IMODT)
         IF (IX.LT.0) IX = IX + IMAGX
         IF (IY.LT.0) IY = IY + IMAGY
         IF (IZ.LT.0) IZ = IZ + IMAGZ
         IF (IT.LT.0) IT = IT + IMAGT
C
C        Convert the current state into a value between 0 and 1
         W = IX/DBLE(IMAGX) + IY/DBLE(IMAGY) + IZ/DBLE(IMAGZ) +
     *       IT/DBLE(IMAGT)
         X(I) = W - INT(W)
   20 CONTINUE
C
C     Put the current generator state into the STATE array
      STATE(5) = IX
      STATE(6) = IY
      STATE(7) = IZ
      STATE(8) = IT
C
      INFO = 0
      RETURN
      END
