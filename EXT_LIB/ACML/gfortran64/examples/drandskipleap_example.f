C     Example program for DRANDLEAPFROG and DRANDSKIPAHEAD
C     .. Parameters ..
      INTEGER          MSTATE, MSEED, MN
      PARAMETER        (MSTATE=20,MSEED=4,MN=12)
C     .. Local Scalars ..
      DOUBLE PRECISION A, B
      INTEGER          GENID, I, INFO, LSEED, LSTATE, N, NSKIP, SUBID
C     .. Local Arrays ..
      DOUBLE PRECISION X1(MN), X2(MN), X3(MN), X4(MN)
      INTEGER          SEED(MSEED), STATE1(MSTATE), STATE2(MSTATE),
     *                 STATE3(MSTATE), STATE4(MSTATE)
C     .. External Subroutines ..
      EXTERNAL         DRANDUNIFORM, DRANDINITIALIZE, DRANDLEAPFROG,
     *                 DRANDSKIPAHEAD
C     .. Intrinsic Functions ..
      INTRINSIC        INT
C     .. Executable Statements ..
      CONTINUE

      WRITE (*,FMT=*)
      WRITE (*,FMT=99999) 'ACML example: Multiple streams using'//
     *  ' DRANDSKIPAHEAD and DRANDLEAPFROG'
      WRITE (*,FMT=99999) '------------------------------------'//
     *  '---------------------------------'
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

C     Generate values from a uniform U(0,1) distribution
      A = 0.0D0
      B = 1.0D0

C     --------------------------------------------------------------
C     Creating multiple streams using SKIPAHEAD
      WRITE (*,FMT=99999) 'Generating multiple streams using SKIPAHEAD:'

C     Initialize the base generator
      CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE1,LSTATE,INFO)

C     Going to be using four streams to generate N values,
C     so initialize all four to the same values
      DO 20 I = 1, LSTATE
         STATE2(I) = STATE1(I)
         STATE3(I) = STATE1(I)
         STATE4(I) = STATE1(I)
   20 CONTINUE

C     Calculate the number of values from each stream
      NSKIP = INT(N/4)
      IF (4*NSKIP.LT.N) NSKIP = NSKIP + 1

C     Advance each stream the required amount
      CALL DRANDSKIPAHEAD(NSKIP,STATE2,INFO)
      CALL DRANDSKIPAHEAD(2*NSKIP,STATE3,INFO)
      CALL DRANDSKIPAHEAD(3*NSKIP,STATE4,INFO)

C     Generate a sequence from all four streams
      CALL DRANDUNIFORM(NSKIP,A,B,STATE1,X1,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE2,X2,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE3,X3,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE4,X4,INFO)

C     Report the values
      DO 40 I = 1, NSKIP
         WRITE (*,FMT=99998) X1(I), X2(I), X3(I), X4(I)
   40 CONTINUE


C     --------------------------------------------------------------
C     Creating multiple streams using LEAPFROG

      WRITE (*,FMT=*)
      WRITE (*,FMT=99999) 'Generating multiple streams using LEAPFROG:'

C     Initialize the base generator
C     In practice the generators should only be initialized once per
C     application. They are reinitialized here so that the results from
C     SKIPAHEAD and LEAPFROG can be compared
      CALL DRANDINITIALIZE(GENID,SUBID,SEED,LSEED,STATE1,LSTATE,INFO)

C     Going to be using four streams to generate N values,
C     so initialize all four to the same values
      DO 60 I = 1, LSTATE
         STATE2(I) = STATE1(I)
         STATE3(I) = STATE1(I)
         STATE4(I) = STATE1(I)
   60 CONTINUE

C     Calculate the number of values from each stream
      NSKIP = INT(N/4)
      IF (4*NSKIP.LT.N) NSKIP = NSKIP + 1

C     Advance each stream the required amount
      CALL DRANDLEAPFROG(4,1,STATE1,INFO)
      CALL DRANDLEAPFROG(4,2,STATE2,INFO)
      CALL DRANDLEAPFROG(4,3,STATE3,INFO)
      CALL DRANDLEAPFROG(4,4,STATE4,INFO)

C     Generate a sequence from all four streams
      CALL DRANDUNIFORM(NSKIP,A,B,STATE1,X1,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE2,X2,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE3,X3,INFO)
      CALL DRANDUNIFORM(NSKIP,A,B,STATE4,X4,INFO)

C     Report the values
      DO 80 I = 1, NSKIP
         WRITE (*,FMT=99998) X1(I), X2(I), X3(I), X4(I)
   80 CONTINUE
C
99999 FORMAT (A)
99998 FORMAT (4F10.4)
      END
