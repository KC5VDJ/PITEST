C*****************************************************************************
C
C   MULTIPRECISION FLOATING POINT COMPUTATION PACKAGE
C   DOUBLE PRECISION RADIX 2 COMPLEX FFT VERSION
C
C   This version has been modified for higher performance when the actual level
C   of precision of input numbers (i.e., the index of the last nonzero words)
C   varies widely.  The advanced algorithms in MPMULX, MPDIVX, and MPSQRX
C   are only performed when the actual precision levels of the input numbers
C   are above a certain level.  For this reason, it is recommended that users
C   always call MPMULX, MPDIVX, or MPSQRX, and not directly call MPMUL, MPDIV,
C   or MPSQRT.
C
C   Version Date:  June 23, 1988
C
C   Author:
C      David H. Bailey                 Telephone:   415-694-4410
C      NASA Ames Research Center       E-mail:  dbailey@ames-nas.arpa
C      Mail Stop 258-5
C      Moffett Field, CA  94035
C
C   Introduction:
C
C       This package of subroutines performs multiprecision floating point
C   arithmetic.  It employs the most advanced algorithms known as of the
C   above date, and all programs have been written to effectively utilize
C   the vector processing capabilities of a supercomputer.  This package
C   assumes that all double precision variables represent eight bytes of
C   storage, with at least 46 mantissa bits.  No other machine dependent
C   assumptions have been made.
C
C   IMPORTANT NOTE: ALL floating-point variables in the following package are
C   assumed to be double-precision (64 bit) floating-point variables.  In
C   addition, double precision data is to be inferred in ALL references to
C   "cells", "locations", etc.  In particular, ALL scratch space requirements
C   are to be interpreted in units of 64-bit storage locations.
C
C       This package includes the following subroutines  (SP denotes
C   ordinary single precision, MP denotes multiprecision (NW mantissa words),
C   and DMP denotes double multiprecision (2 * NW mantissa words):
C
C   MPSMC        Performs SP to MP conversion.
C   MPMSC        Performs MP to SP conversion.
C   MPEQ         Sets one MP number equal to another.
C   MPADD        Adds two MP numbers to yield MP sum.
C   MPSUB        Subtracts two MP numbers to yield MP difference.
C   MPMUL1       Multiplies a MP number by a SP number to yield MP product.
C   MPDIV1       Divides a MP number by a MP divisor to yield SP quotient.
C   MPDIV2       Divides a MP number by a SP divisor to yield MP quotient.
C   MPNORM       Normalizes a MP number to standard MP form.
C   MPRAND       Returns a pseudo-random MP number between 0 and 1.
C   MPINFR       Returns the integer and fractional part of a MP number.
C   MPFMT        Formats a MP number for output.
C   MPCFFT       Initializes the array U in common MPCOM2 and performs FFT.
C   MPMULX       Multiplies two MP numbers to yield DMP product.
C   MPDIVX       Divides a MP number by a MP divisor to yield MP quotient.
C   MPSQRX       Computes the square root of a MP number.
C
C   A number of other subroutines are part of the package but are not listed
C   here because they are not intended to be directly called by users.
C
C        The format of the multiprecision MP numbers is as follows:  the
C   first word contains the sign (1., 0. or -1.), and the second word
C   contains the exponent (powers of 10^6).  Words numbered 3 up to NW + 2
C   contain the mantissa.  If the sign is zero, then all other words must
C   be zero.  If the sign is nonzero, then word 3 must be nonzero.  The
C   decimal point is assumed after the first mantissa word for numbers with
C   with zero exponent.  The radix of this representation is 10^6, so that
C   each mantissa word contains a floating whole number between 0 and
C   999999.  The decimal radix 10^6 has been used instead of a binary radix
C   so that expensive binary to decimal conversion is not necessary for output.
C
C        The format of the double multiprecision DMP numbers is the same
C   as the MP numbers except that DMP numbers have twice as many mantissa
C   words, for a total of 2 * NW + 2 words.  Note that only MPMULX explicitly
C   computes with DMP numbers.  Other operations may be easily performed on
C   DMP numbers by merely doubling NW and, if necessary, increasing MW by one.
C   If this is done be sure that the scratch arrays in calls to MPMULX,
C   MPDIVX, and MPSQRX have sufficient scratch space, and remember to restore
C   NW and MW to their previous values when the DMP computations are complete.
C
C        To initialize the package, the variables ND, NW, MW, IDB, LDB, and NDB
C   must be set in the common block MPCOM1.  ND is the requested number of
C   digits of precision.  NW is the number of words of precision, which should
C   be at least ND / 6 + 1.  All MP variables must be dimensioned at least
C   NW + 2 and DMP variables must be dimensioned at least 2 * NW + 2 (twice
C   the amount for MP variables will suffice).  MW is used in subroutines
C   MPMULX, MPDIVX, and MPSQRX.  When ND is less than about 200 (depending on
C   the system), MW may be set to zero.  For higher levels of precision, the
C   only permissible values of of NW are 2^(MW - 2), and both MW and NW must
C   be set in common MPCOM1 by the user.  In addition, prior to calling
C   MPMULX, MPDIVX, or MPSQRX for these higher levels of precision, the user
C   must first allocate at least 2 * 2^MX cells in the array U of common block
C   MPCOM2 in the user's main program, and the user must initialize the FFT
C   routine by calling MPCFFT with 0 as the first argument.  Here MX is the
C   largest value of MW that will be used.
C
C        The variable IDB is a debug flag and ordinarily should be set to 0.
C   Setting IDB to an integer between 5 and 9 (or greater) produces debug
C   printouts in varying degrees from the MP subroutines.  Values of IDB
C   between 1 and 4 are available for use as debug flags in the user's
C   calling program if desired.  LDB is the logical unit number for output of
C   these debug messages, and NDB is the number of words output in the debug
C   printout of a MP number.  Typically LDB is set to 6 and NDB is set to 16.
C
C        Common data has been divided into two separate common blocks to
C   facilitate efficient multitasking.  The parameters in MPCOM1 may be
C   frequently changed by a user program and thus should be declared TASK
C   COMMON or the equivalent.  The array U in common MPCOM2, however,
C   is not intended to be altered after it is first initialized.  For this
C   reason, and since U may be a very large array, it makes more sense for
C   MPCOM2 to be an ordinary global common.  Except for the data in these
C   two common blocks, no variables are global or need to be "saved".
C
C*****************************************************************************
C
      SUBROUTINE MPSMC (A, N, B)
C
C   Converts the SP number  A * 10^N  to MP form in  B.
C
C   Note:  The result should be exact provided N is 0 or divisible by IB
C   and B is an exact binary fraction in the range RX^2 .LT. |B| .LT. BX^2.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RX2 = RX * RX,
     $  RXX = 0.5D0 * RX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION B(1)
C
      NX = NW + 2
      IF (IDB .GE. 8)  WRITE (LDB, 1)  A, N, NX
1     FORMAT (/'MPSMC I',1PD25.15,2I10)
C
C   Check for zero.
C
      IF (A .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          B(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
C   Adjust exponent value.
C
      AA = ABS (A)
      AL = N + LOG10 (AA) * (1.D0 + RX2)
      M = AL / IB
      IF (AL .LT. 0.D0 .AND. AL .NE. M * IB)  M = M - 1
      AM = AA * 10.D0 ** (IB - IB * M + N)
C
C   Convert constant to 2-word precision.
C
      B(1) = SIGN (1.D0, A)
      B(2) = M
      B(3) = AINT (RX * AM + RXX)
      B(4) = AINT (AM - BX * B(3))
C
      DO 110 I = 5, NX
        B(I) = 0.D0
110   CONTINUE
C
200   IF (IDB .GE. 8)  WRITE (LDB, 2)  (B(I),I=1,NDB)
2     FORMAT ('MPSMC O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPMSC (A, B, N)
C
C   Converts the MP number  A  to the SP form  B * 10^N , accurate to about
C   12 decimal places.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RX2 = RX * RX)
      DIMENSION A(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
C
      NX = NW + 2
      IF (IDB .GE. 8)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPMSC I',I10/(8F9.0))
C
      IF (A(1) .EQ. 0.D0)  THEN
        B = 0.D0
        N = 0
        GOTO 200
      ENDIF
C
      AA = A(3) + A(4) * RX + A(5) * RX2
C
C   Adjust results so that  B  is in the range  1 .LE. |B| .LT. 10.
C
      N2 = LOG10 (AA) + RX2
      B = SIGN (AA * 10.D0 ** (-N2), A(1))
      N = N2 + IB * A(2)
C
200   IF (IDB .GE. 8)  WRITE (LDB, 2)  B, N
2     FORMAT ('MPMSC O',F10.0,I10)
      RETURN
      END
C
      SUBROUTINE MPEQ (A, B)
C
C   Sets the MP number  B  equal to the MP number  A.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), B(1)
C
      IF (IDB .GE. 8)  WRITE (LDB, 1)
1     FORMAT (/'MPEQ')
C
      DO 100 I = 1, NW + 2
        B(I) = A(I)
100   CONTINUE
      RETURN
      END
C
      SUBROUTINE MPADD (A, B, C)
C
C   Adds MP numbers  A  and  B  to yield the MP sum  C.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), B(1), C(1)
C
      NX = NW + 2
      IF (IDB .GE. 7)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPADD I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
C
      IF (A(1) * B(1) .EQ. 0.D0)  THEN
C
C   A or B is zero.
C
        IF (A(1) .EQ. 0.D0)  THEN
          DO 80 I = 1, NX
            C(I) = B(I)
80        CONTINUE
          GOTO 300
        ELSE
          DO 90 I = 1, NX
            C(I) = A(I)
90        CONTINUE
          GOTO 300
        ENDIF
C
      ELSEIF (A(1) * B(1) .GT. 0.D0)  THEN
C
C   A  and  B  have the same sign or one is zero -- perform addition
C   after shifting to match exponents.
C
        C(1) = A(1)
        C(2) = MAX (A(2), B(2))
C
        IF (A(2) .GE. B(2))  THEN
          N = MIN (INT (A(2) - B(2)), NW)
          DO 100 I = 3, N+2
            C(I) = A(I)
100       CONTINUE
          DO 110 I = N+3, NX
            C(I) = A(I) + B(I-N)
110       CONTINUE
C
        ELSE
          N = MIN (INT (B(2) - A(2)), NW)
          DO 120 I = 3, N+2
            C(I) = B(I)
120       CONTINUE
          DO 130 I = N+3, NX
            C(I) = A(I-N) + B(I)
130       CONTINUE
        ENDIF
C
      ELSE
C
C   A  and  B  have different nonzero signs -- perform subtraction.
C
        IF (A(2) .EQ. B(2))  THEN
C
C   A  and  B  have the same exponent -- search for unequal cells.
C
          C(2) = 0.D0
          DO 140 I = 3, NX
            IF (A(I) .NE. B(I))  GOTO 160
140       CONTINUE
C
C   ABS (A)  =  ABS (B).  The result is zero.
C
          DO 150 I = 1, NX
            C(I) = 0.D0
150       CONTINUE
          GOTO 300
C
160       K = I
C
C   ABS (A)  .NE.  ABS (B)  but exponents are equal -- subtract and normalize.
C
          IF (A(K) .GT. B(K))  THEN
            K = K - 3
            C(1) = A(1)
            C(2) = A(2) - K
            DO 170 I = 3, NX - K
              C(I) = A(I+K) - B(I+K)
170         CONTINUE
C
          ELSE
            K = K - 3
            C(1) = B(1)
            C(2) = B(2) - K
            DO 180 I = 3, NX - K
              C(I) = B(I+K) - A(I+K)
180         CONTINUE
          ENDIF
C
          DO 190 I = NX - K + 1, NX
            C(I) = 0.D0
190       CONTINUE
C
        ELSE
C
C   Exponents are different -- subtract with smaller shifted right.
C
          IF (A(2) .GT. B(2))  THEN
            C(1) = A(1)
            C(2) = A(2)
            N = MIN (INT (A(2) - B(2)), NW)
            DO 200 I = 3, N+2
              C(I) = A(I)
200         CONTINUE
            DO 210 I = N+3, NX
              C(I) = A(I) - B(I-N)
210         CONTINUE
C
          ELSE
            C(1) = B(1)
            C(2) = B(2)
            N = MIN (INT (B(2) - A(2)), NW)
            DO 220 I = 3, N+2
              C(I) = B(I)
220         CONTINUE
            DO 230 I = N+3, NX
              C(I) = B(I) - A(I-N)
230         CONTINUE
C
          ENDIF
        ENDIF
      ENDIF
C
C   Fix up result (possible negative numbers or overflow).
C
      CALL MPNORM (C)
C
300   IF (IDB .GE. 7)  WRITE (LDB, 2)  (C(I),I=1,NDB)
2     FORMAT ('MPADD O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPSUB (A, B, C)
C
C   Subtracts MP numbers  A  and  B  to yield the MP difference  C.
C
C   Note:  This does not work if  A  and  B  are the same variable, i.e.
C          CALL MPSUB (X, X, Y)  .
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION  A(1), B(1), C(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
C
      IF (IDB .GE. 7)  WRITE (LDB, 1)
1     FORMAT (/'MPSUB')
C
C   Negate  B  and add.
C
      B(1) = - B(1)
      CALL MPADD (A, B, C)
C
C   Restore the sign of  B.
C
      B(1) = - B(1)
      RETURN
      END
C
      SUBROUTINE MPMUL (A, B, C, NA, NB)
C
C   Multiplies two MP numbers  A  and  B  to yield the exact DMP product  C.
C   Note that the array  C  must contain at least  2 * NW  cells.  NA and NB
C   are the indices of the last nonzero words of A and B.  It is recommended
C   that MPMUL is only called by MPMULX, not directly by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), B(1), C(1)
C
      NX = NW + 2
      IF (IDB .GE. 6)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPMUL I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
C
      NWS = NW
      NW2 = 2 * NW
      NX2 = NW2 + 2
C
      DO 100 I = 1, NX2
        C(I) = 0.D0
100   CONTINUE
      IF (A(1) * B(1) .EQ. 0.D0)  GOTO 200
C
C   Ordinary long multiplication algorithm.
C
      DO 170 J = 3, NA
        AJ = A(J)
        J3 = J - 3
C
        DO 150 I = 3, NB
          C(I+J3) = C(I+J3) + AJ * B(I)
150     CONTINUE
C
C   Release carries periodically.
C
        IF (J .EQ. NA .OR. J .EQ. 64 * (J / 64)) THEN
          K1 = MAX (3, J - 64)
          K2 = NB + J3
C
          DO 160 I = K1, K2
            T = C(I)
            R = AINT (T * RX)
            C(I) = T - R * BX
            C(I-1) = C(I-1) + R
160       CONTINUE
        ENDIF
170   CONTINUE
C
      R = C(2)
      C(1) = A(1) * B(1)
      C(2) = A(2) + B(2)
      C(3) = C(3) + R * BX
C
C   Fix up result (some words may exceed BX).
C
      NW = NW2
      CALL MPNORM (C)
      NW = NWS
C
200   IF (IDB .GE. 6)  WRITE (LDB, 2)  (C(I),I=1,NDB)
2     FORMAT ('MPMUL O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPMUL1 (A, B, N, C)
C
C   Multiplies the MP number A by the SP number B * 10^N to yield the
C   MP product C.
C
C   Note:  The result should be exact provided N is 0 or divisible by IB
C   and B is an exact binary fraction in the range RX^2 .LT. |B| .LT. BX^2.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RX2 = RX * RX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), C(1), BB(2)
C
      NX = NW + 2
      IF (IDB .GE. 7)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPMUL1 I',I10/(8F9.0))
        WRITE (LDB, 2)  B, N
2       FORMAT ('MPMUL1 I',1PD25.15,I10)
      ENDIF
C
C   Check for zero inputs.
C
      DO 100 I = 1, NX
        C(I) = 0.D0
100   CONTINUE
      IF (A(1) * B .EQ. 0.D0)  GOTO 200
C
C   Adjust exponent value.
C
      BA = ABS (B)
      BL = N + LOG10 (BA) * (1.D0 + RX2)
      M = BL / IB
      IF (BL .LT. 0.D0 .AND. BL .NE. M * IB)  M = M - 1
      BM = BA * 10.D0 ** (IB - IB * M + N)
      BB(1) = AINT (RX * BM)
      BB(2) = AINT (BM - BX * BB(1))
C
C   Multiply with 2-word precision.
C
      DO 120 J = 1, 2
        BJ = BB(J)
        J1 = J - 1
C
        DO 110 I = 3, NX - J1
          C(I+J1) = C(I+J1) + A(I) * BJ
110     CONTINUE
120   CONTINUE
C
C   Release the large carries.
C
      DO 130 I = 3, NX
        T = C(I)
        R = AINT (T * RX)
        C(I) = T - R * BX
        C(I-1) = C(I-1) + R
130   CONTINUE
C
      R = C(2)
      C(1) = A(1) * SIGN (1.D0, B)
      C(2) = A(2) + M
      C(3) = C(3) + R * BX
      CALL MPNORM (C)
C
200   IF (IDB .GE. 7)  WRITE (LDB, 3)  (C(I),I=1,NDB)
3     FORMAT ('MPMUL1 O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPDIV (A, B, C, D)
C
C   Divides the MP number  A  by the MP number  B  to yield the MP
C   quotient  C.   D  is a scratch array with at least  NW  locations.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, BX2 = BX * BX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), B(1), C(1), D(1)
C
      NX = NW + 2
      IF (IDB .GE. 6)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPDIV I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
C
C   Check if dividend is zero.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          C(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
C   Check if divisor is zero.
C
      IF (B(1) .EQ. 0.D0)  THEN
        WRITE (LDB, 2)
2       FORMAT ('*** MPDIV: ZERO DIVISOR ***')
        CALL XABORT
      ENDIF
C
C   Initialize trial dividend.
C
      D(1) = 0.D0
      DO 110 I = 2, NW + 1
        D(I) = A(I+1)
110   CONTINUE
      D(NX) = 0.D0
C
C   Compute trial divisor.
C
      RB = 1.D0 / (BX * B(3) + B(4) + RX * B(5))
C
      DO 130 J = 2, NX - 1
C
C   Compute current quotient.
C
        S = AINT (RB * (BX2 * D(J-1) + BX * D(J) + D(J+1)))
        C(J) = S
        J3 = J - 3
        IN = NX - J3
        IF (J .EQ. 2)  IN = IN - 1
C
        DO 120 I = 3, IN
          IJ = I + J3
          T = D(IJ) - S * B(I)
          R = AINT (T * RX - 1.D0)
          D(IJ) = T - R * BX
          D(IJ-1) = D(IJ-1) + R
120     CONTINUE
        D(J) = D(J) + BX * D(J-1)
130   CONTINUE
C
C   Set sign and exponent of result.
C
      C(NX) = AINT (RB * (BX2 * D(NX-1) + BX * D(NX)))
      C(3) = BX * C(2) + C(3)
      C(1) = A(1) * B(1)
      C(2) = A(2) - B(2) - 1
      CALL MPNORM (C)
C
200   IF (IDB .GE. 6)  WRITE (LDB, 3)  (C(I),I=1,NDB)
3     FORMAT ('MPDIV O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPDIV1 (A, B, C, N)
C
C   Divides MP  A  by MP  B  to yield the SP quotient  C * 10^N,  accurate
C   to about 12 decimal places.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RX2 = RX * RX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), B(1)
C
      NX = NW + 2
      IF (IDB .GE. 7)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPDIV1 I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
C
C   Check for zero dividend.
C
      IF (A(1) .EQ. 0.D0)  THEN
        C = 0.D0
        N = 0
        GOTO 200
      ENDIF
C
C   Check for zero divisor.
C
      IF (B(1) .EQ. 0.D0)  THEN
        WRITE (LDB, 2)
2       FORMAT ('*** MPDIV1: ZERO DIVISOR ***')
        CALL XABORT
      ENDIF
C
      AA = SIGN (A(3) + A(4) * RX + A(5) * RX2, A(1))
      BB = SIGN (B(3) + B(4) * RX + B(5) * RX2, B(1))
      C = AA / BB
      N = IB * (A(2) - B(2))
C
C   Adjust results so  C  is in the range  1 .LE. |C| .LT. 10.
C
      CA = ABS (C)
      IF (CA .GE. 1.D0)  THEN
        N2 = LOG10 (CA) + RX2
      ELSE
        N2 = LOG10 (CA) - 1.D0 - RX2
      ENDIF
      C = C * 10.D0 ** (-N2)
      N = N + N2
C
200   IF (IDB .GE. 7)  WRITE (LDB, 3)  C, N
3     FORMAT ('MPDIV1 O',F10.0,I10)
      RETURN
      END
C
      SUBROUTINE MPDIV2 (A, B, N, C)
C
C   Divides the MP number  A  by the SP number  B * 10^N  to obtain the MP
C   quotient  C.
C
C   Note:  The result should be exact provided N is 0 or divisible by IB
C   and B is an exact binary fraction in the range RX^2 .LT. |B| .LT. BX^2.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, BX2 = BX * BX,
     $  RX2 = RX * RX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1), C(1)
C
      NX = NW + 2
      IF (IDB .GE. 7)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPDIV2 I',I10/(8F9.0))
        WRITE (LDB, 2)  B, N
2       FORMAT ('MPDIV2 I',1PD25.15,I10)
      ENDIF
C
C   Check for zero dividend.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          C(I) = 0.D0
100     CONTINUE
        D = 0.D0
        GOTO 200
      ENDIF
C
C   Check for zero divisor.
C
      IF (B .EQ. 0.D0)  THEN
        WRITE (LDB, 3)  B
3       FORMAT ('*** MPDIV2:  ZERO DIVISOR ***', 1PD15.6)
        CALL XABORT
      ENDIF
C
C   Adjust exponent value.
C
      BA = ABS (B)
      BL = N + LOG10 (BA) * (1.D0 + RX2)
      M = BL / IB
      IF (BL .LT. 0.D0 .AND. BL .NE. M * IB)  M = M - 1
      BM = BA * 10.D0 ** (IB - IB * M + N)
      B1 = AINT (RX * BM)
      B2 = AINT (BM - BX * B1)
      RB = 1.D0 / (BX * B1 + B2)
      D1 = 0.D0
      D2 = A(3)
      D3 = A(4)
C
C   Perform short division (not vectorizable at present).
C
      DO 110 J = 2, NX - 1
        S = AINT (RB * (BX2 * D1 + BX * D2 + D3))
        C(J) = S
        T = D3 - S * B2
        R = AINT ((T + 1) * RX - 1.D0)
        D3 = T - R * BX
        T = R + D2 - S * B1
        R = AINT ((T + 1) * RX - 1.D0)
        D2 = T - R * BX
        D1 = (R + D1) * BX + D2
        D2 = D3
        D3 = A(J+3)
        IF (J .GE. NW)  D3 = 0.D0
110   CONTINUE
C
C   Set sign and exponent of result.
C
      C(NX) = AINT (RB * (BX2 * D1 + BX * D2))
      C(3) = BX * C(2) + C(3)
      C(1) = A(1) * SIGN (1.D0, B)
      C(2) = A(2) - M - 1
      CALL MPNORM (C)
C
200   IF (IDB .GE. 7)  WRITE (LDB, 4)  (C(I),I=1,NDB)
4     FORMAT ('MPDIV2 O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPNORM (A)
C
C   Converts the MP number  A  to the standard normalized form.  Many
C   MP routines call MPNORM to clean up results, because occasionally these
C   operations leave negative numbers or values exceeding the radix  BX  in
C   result arrays.  The user need not call MPNORM directly.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RXZ = 2.D0 + 0.5D0 * RX)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1)
C
      NX = NW + 2
      IF (IDB .GE. 9)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPNORM I',I10/(4F18.0))
C
      IF (A(1) .EQ. 0.D0)  GOTO 200
      A2 = A(2)
      A(2) = 0.D0
C
C   Try vectorized fixup loop three times (this should fix 99%).
C   This loop should be commented out for execution on scalar computers.
C
      DO 110 K = 1, 3
        SR = 0.D0
        DO 100 I = 3, NX
          R = AINT (A(I) * RX + RXZ) - 2.D0
          A(I) = A(I) - R * BX
          A(I-1) = A(I-1) + R
          SR = SR + ABS (R)
100     CONTINUE
        IF (SR .EQ. 0.D0)  GOTO 130
110   CONTINUE
C
C   Still not fixed - use recursive (unvectorizable) loop.
C
      R = 0.D0
      DO 120 I = NX, 3, -1
        T = R + A(I)
        R = AINT (T * RX + RXZ) - 2.D0
        A(I) = T - R * BX
120   CONTINUE
      A(2) = A(2) + R
C
C   Check for overflow.
C
130   IF (A(2) .NE. 0.D0)  THEN
        DO 140 I = NX, 3, -1
          A(I) = A(I-1)
140     CONTINUE
        A2 = A2 + 1.D0
      ENDIF
C
C   Check if first few mantissa words are zero.
C
      IF (A(3) .EQ. 0.D0)  THEN
        DO 150 I = 4, NX
          IF (A(I) .NE. 0.D0)  GOTO 160
150     CONTINUE
160     K = I - 3
CDIR$ IVDEP
CDEC$ IVDEP
        DO 170 I = 3, NX - K
          A(I) = A(I+K)
170     CONTINUE
        DO 180 I = NX - K + 1, NX
          A(I) = 0.D0
180     CONTINUE
        A2 = A2 - K
      ENDIF
C
      A(2) = A2
C
200   IF (IDB .GE. 9)  WRITE (LDB, 2)  (A(I),I=1,NDB)
2     FORMAT ('MPNORM O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPSQRT (A, B, C)
C
C   Computes the square root of the MP number  A  and returns the MP result
C   in  B.  C  is a scratch array with at least 2 * NW + 4 cells (this is
C   double the amount required for MP numbers).
C
C   This subroutine employs a Newton-Raphson iteration that converges to the
C   square root of A:
C
C          x(n+1) = .5 * (x(n) + a / x(n))
C
C   Note that this subroutine dynamically changes the precision level, doubling
C   it with each iteration.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(1), B(1), C(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, CL2 = 1.442695041D0)
C
      NX = NW + 2
      IF (IDB .GE. 5)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPSQRT I',I10/(8F9.0))
      IDS = IDB
      IF (IDB .LE. 7)  IDB = 0
C
C   Check for zero input.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          B(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
      IF (A(1) .LT. 0.D0)  THEN
        WRITE (LDB, 2)
2       FORMAT ('*** MPSQRT: NEGATIVE INPUT ***')
        CALL XABORT
      ENDIF
C
      NWS = NW
      NX1 = NX + 1
      L = CL2 * LOG (DBLE (NW)) + 1.D0 + RX
C
C   Compute initial approximation.
C
      SA = A(3) + A(4) * RX + A(5) * RX ** 2
      NA = IB * A(2)
      N2 = NA / 2
      IF (NA .NE. 2 * N2)  THEN
        IF (N2 .LT. 0)  N2 = N2 - 1
        SA = SA * 10.D0
      ENDIF
      S2 = SQRT (SA)
      CALL MPSMC (S2, N2, B)
      NW = 2
C
C   Repeat Newton-Raphson iterations.
C
      DO 110 I = 1, L
C
C   Repeat the next-to-last step to insure maximum accuracy.
C
        IF (I .NE. L - 1)  NW = MIN (2 * NW, NWS)
        CALL MPDIV (A, B, C, C(NX1))
        CALL MPADD (B, C, C(NX1))
        CALL MPMUL1 (C(NX1), 0.5D0, 0, B)
110   CONTINUE
C
200   IDB = IDS
      IF (IDB .GE. 5)  WRITE (LDB, 3)  (B(I),I=1,NDB)
3     FORMAT ('MPSQRT O'/(8F9.0))
C
      RETURN
      END
C
      SUBROUTINE MPRAND (A)
C
C   Returns a pseudo-random MP number  A  between  0  and  1.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION A(1)
C
      NX = NW + 2
      A(1) = 1.D0
      A(2) = -1.D0
C
      DO 100 I = 3, NX
        A(I) = AINT (BX * XRAND (0))
100   CONTINUE

      IF (IDB .GE. 7)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPRAND',I10/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPINFR (A, B, C)
C
C   Sets  B  equal to the integer part of the MP number  A
C   and sets  C  equal to the fractional part of  A.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DIMENSION  A(1), B(1), C(1)
C
      NX = NW + 2
      IF (IDB .GE. 7)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPINFR I',I10/(8F9.0))
C
C   Check if  A  is zero.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          B(I) = 0.D0
          C(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
C   Find word that precedes decimal point.
C
      N = MIN (INT (A(2)) + 3, NX)
      IF (N .LT. 3)  N = 0
C
C   Place integer part in  B.
C
      DO 110 I = 1, N
        B(I) = A(I)
110   CONTINUE
      DO 120 I = N + 1, NX
        B(I) = 0.D0
120   CONTINUE
C
C   Integer part is zero,  fractional part is  A.
C
      IF (N .EQ. 0)  THEN
        C(1) = A(1)
        C(2) = A(2)
        GOTO 150
      ENDIF
C
C   Check if initial cells of fractional part are zero.
C
      DO 130 I = N + 1, NX
        IF (A(I) .NE. 0.D0)  GOTO 140
130   CONTINUE
C
C   The fractional part is zero.
C
      C(1) = 0.D0
      C(2) = 0.D0
      N = NX
      GOTO 150
C
C   Store the fractional part in  C.
C
140   C(1) = A(1)
      C(2) = N - I
      N = I - 3
150   DO 160 I = 3, NX - N
        C(I) = A(I+N)
160   CONTINUE
      DO 170 I = NX - N + 1, NX
        C(I) = 0.D0
170   CONTINUE
C
200   IF (IDB .GE. 7)  THEN
        WRITE (LDB, 2)  (B(I),I=1,NDB)
2       FORMAT ('MPINFR O'/(8F9.0))
        WRITE (LDB, 2)  (C(I),I=1,NDB)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPFMT (A, B)
C
C   Formats the MP number  A  into character form of length  ND + 20  in the
C   character array  B.  It is similar to the Fortran exponential format
C   (E format), except that the exponent is placed first.  B may be output
C   using A1 format.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*1 B
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, RX2 = RX * RX)
      DIMENSION  A(1), B(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
C
      NX = NW + 2
      IF (IDB .GE. 8)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPFMT I',I10/(8F9.0))
C
C   Determine exact power of ten for exponent.
C
      K = 0
      IF (A(3) .NE. 0)  K = LOG10 (A(3)) + RX2
      AN = IB * A(2) + K
C
C   Place exponent first instead of at the very end as in Fortran.
C
      B(1) = '1'
      B(2) = '0'
      B(3) = ' '
      B(4) = '^'
      CALL MPENC (1, AN, B(5), 10, 1)
      B(15) = ' '
      B(16) = 'x'
      B(17) = ' '
C
C
C   Insert sign.
C
      B(18) = ' '
      IF (A(1) .LT. 0.)  B(18) = '-'
C
C   Insert first mantissa digit of number, followed by a period.
C
      AN = AINT (A(3) * 10.D0 ** (-K) + RX2)
      CALL MPENC (1, AN, B(19), 1, 0)
      B(20) = '.'
C
C   Insert the remaining digits of the first word.
C
      AN = AINT (A(3) - AN * 10.D0 ** K + RX2)
      CALL MPENC (1, AN, B(21), K, 0)
C
C   Insert the digits of the remaining words.
C
      CALL MPENC (ND / IB, A(4), B(K+21), IB, 0)
C
200   RETURN
      END
C
      SUBROUTINE MPENC (N, A, B, LEN, IBL)
C
C   Encodes  N  floating whole numbers in the array  A  into the character
C   array  B.  the length of each number is  LEN  characters.  If  IBL  is
C   zero, then the numbers are left-filled with zeroes.  Otherwise the numbers
C   are left-filled with blanks.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*1 B, DIG
      DIMENSION A(1), B(1), DIG(0:9)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      DATA DIG/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
C
      IF (LEN .EQ. 0)  GOTO 200
      IX = 0
      DO 130 J = 1, N
        AJ = A(J)
        XJ = ABS (AJ)
        IF (AJ .GE. 0.)  THEN
          L1 = 1
        ELSE
          L1 = 2
          B(IX+1) = ' '
        ENDIF
C
C   Extract digits one at a time, beginning with units digit.
C
        DO 100 I = LEN, L1, -1
          XI = AINT (0.1D0 * XJ + 0.01D0)
          K = XJ - 10. * XI
          XJ = XI
          B(IX+I) = DIG(K)
100     CONTINUE
C
C   Check if input is too large for field.
C
        IF (XI .NE. 0.)  THEN
          WRITE (LDB, 1)  J, AJ, LEN
1         FORMAT ('*** MPENC: INPUT TOO LARGE ***', I5, F15.1, I5)
          CALL XABORT
        ENDIF
C
C   Left-fill with blanks if requested.
C
        IF (IBL .EQ. 1)  THEN
          DO 110 I = L1, LEN - 1
            IF (B(IX+I) .NE. '0')  GOTO 120
            B(IX+I) = ' '
110       CONTINUE
          I = LEN
C
C   Insert sign if required to the left of first nonzero digit.
C
120       IF (AJ .LT. 0)  B(IX+I-1) = '-'
        ENDIF
        IX = IX + LEN
130   CONTINUE
C
200   RETURN
      END
C
      SUBROUTINE MPMULX (A, B, C, D)
C
C   Multiplies MP numbers A and B to yield the DMP product C.  The array D is
C   a scratch array with at least 12 * NW + 6 cells (twelve times the amount
C   allocated for MP numbers will suffice).  The value of NW must be equal to
C   2^(MW - 2).  Subroutine MPCFFT must have previously been called to
C   initialize the array U in common block MPCOM2.
C
C   This version checks the input numbers A and B to determine the locations of
C   the last nonzero words.  If this location is low, so that the actual
C   precision of the number is low, then MPMUL is called instead, and these
C   locations and passed to it.
C
C   This subroutine uses a special technique involving a discrete convolution
C   and a fast Fourier transform.  The work factor for this algorithm is
C   proportional to NW * log (NW), which is much faster than the NW^2 algorithm
C   used in MPMUL for large values of NW.  The parameter MI is the crossover
C   point for MW to apply the advanced algorithm.  On Cray systems MI should
C   be set to 8.  On other systems another value may be slightly better.  On
C   scalar computers, MI probably should be set to 5.  Change MI here and also
C   in MPDIVX and MPSQRX.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION  A(1), B(1), C(1), D(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      PARAMETER (BX = 1D6, RX = 1D-6, IB = 6, MI = 8, BB = 1D3,
     $   RB = 1D-3)
C
      NX = NW + 2
      IF (IDB .GE. 6)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPMULX I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
C
C   Check if precision level is too low to justify advanced algorithm.  The
C   following code searches backwards from the end of the input arrays to find
C   the indices of the last nonzero words.
C
      DO 100 I = NX, 3, -1
        IF (A(I) .NE. 0.D0)  GOTO 110
100   CONTINUE
110   NA = I
C
      DO 120 I = NX, 3, -1
        IF (B(I) .NE. 0.D0)  GOTO 130
120   CONTINUE
130   NB = I
C
      IF (MW .LE. MI .OR. NA .LT. 64 .OR. NB .LT. 64)  THEN
        CALL MPMUL (A, B, C, NA, NB)
        GOTO 200
      ENDIF
C
C   Check for mismatch between MW and NW.
C
      N1 = 2 ** (MW - 2)
      IF (N1 .NE. NW)  THEN
        WRITE (LDB, 2)  MW, NW
2       FORMAT ('*** MPMULX:  MISMATCHED VALUES OF MW AND NW--', 2I10)
        CALL XABORT
      ENDIF
C
      N2 = 2 * N1
      N21 = N2 + 1
      N4 = 4 * N1
      N42 = N4 + 2
      N6 = 6 * N1
      N63 = N6 + 3
      N8 = 8 * N1
      N84 = N8 + 4
C
C   Check for zero inputs.
C
      IF (A(1) * B(1) .EQ. 0.D0)  THEN
        DO 140 I = 1, N2 + 2
          C(I) = 0.D0
140     CONTINUE
        GOTO 200
      ENDIF
C
C   Split input data in A and B into words with half as many bits and place in
C   (D(I), I=1,N4) and (D(I), I=N4+3,N8+2).  The second half of the two vectors
C   must be filled with zeroes for the FFT multiply technique to work properly.
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 150 I = 1, N1
        T = A(I+2)
        A1 = AINT (RB * T)
        A2 = T - BB * A1
        D(2*I-1) = A1
        D(2*I) = A2
        T = B(I+2)
        B1 = AINT (RB * T)
        B2 = T - BB * B1
        D(2*I-1+N42) = B1
        D(2*I+N42) = B2
150   CONTINUE
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 160 I = N2 + 1, N42
        D(I) = 0.D0
        D(I+N42) = 0.D0
160   CONTINUE
C
C   Perform forward real-to-complex FFTs on the two vectors in D.
C
      CALL MPRCFT (1, D, D(N84+1))
      CALL MPRCFT (1, D(N42+1), D(N84+1))
C
C   Multiply resulting complex arrays.
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 170 I = 1, N21
        D11 = D(I)
        D12 = D(I+N21)
        D21 = D(I+N42)
        D22 = D(I+N63)
        D(I) = D11 * D21 - D12 * D22
        D(I+N21) = D11 * D22 + D12 * D21
170   CONTINUE
C
C   Perform reverse FFT.
C
      CALL MPCRFT (-1, D, D(N84+1))
C
C   Round, recombine into regular words, divide by 2^MW, and release carrys.
C
      AN = 1.D0 / N8
      T1 = AN * D(1)
      A1 = AINT (T1 + 0.5D0)
      D(1) = ABS (A1 - T1)
      C(2) = A1
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 180 I = 1, N2 - 1
        T1 = AN * D(2*I)
        T2 = AN * D(2*I+1)
        A1 = AINT (T1 + 0.5D0)
        A2 = AINT (T2 + 0.5D0)
        D(2*I) = ABS (A1 - T1)
        D(2*I+1) = ABS (A2 - T2)
        C1 = AINT (RX * A1)
        C2 = A1 - BX * C1
        C3 = AINT (RX * A2)
        C4 = A2 - BX * C3
        C(I+2) = BB * C2 + C4
        C(I+1) = C(I+1) + BB * C1 + C3
180   CONTINUE
C
      C(N2+2) = 0.D0
C
C   Find the largest rounding error in the loop above.  This code may be
C   commented out unless very high precision (millions of digits) is used.
C
C      T = 0.D0
C
C      DO 190 I = 1, N4 - 1
C        T = MAX (T, D(I))
C190   CONTINUE
C
C      IF (T .GT. 0.25D0)  THEN
C        WRITE (LDB, 3)  T
C3       FORMAT ('*** MPMULX:  EXCESSIVE ROUNDING ERROR',F8.4)
C        CALL XABORT
C      ENDIF
C
C   Adjust sign and exponent.
C
      C(1) = A(1) * B(1)
      C(3) = C(3) + BX * C(2)
      C(2) = A(2) + B(2)
      NW = N2
      CALL MPNORM (C)
      NW = N1
C
200   IF (IDB .GE. 6)  WRITE (LDB, 4)  (C(I),I=1,NDB)
4     FORMAT ('MPMULX O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPRCFT (IS, X, Y)
C
C   Performs an N-point real-to-complex FFT, where N = 2^MW.  X is both the
C   input and the output data array, and Y is a scratch array.  N real values
C   are input in X, and N/2 + 1 complex values are output in X, with real and
C   imaginary parts separated by N/2 + 1 locations.  Before calling MPRCFT,
C   the U array must be initialized by calling MPCFFT with IS = 0.  A call to
C   MPRCFT with IS = 1 (or -1) indicates a call to perform a real-to-complex
C   FFT with positive (or negative) exponentials.  The arrays X and Y must be
C   dimensioned with N/2 + 1 complex or N + 2 real cells.  U must be
C   dimensioned the same as in MPCFFT.  The output values from MPRCFT are
C   twice as large as the results of a complex-to-complex transform on real
C   data.  MW must be at least three.
C
C   This subroutine employs a technique that converts a real-to-complex FFT
C   to a complex-to-complex FFT.  See Brigham, "The Fast Fourier Transform",
C   p. 169, and Hockney & Jesshope, "Parallel Computers", p. 303, although
C   neither reference is complete and error-free.
C
C   David H. Bailey    May 9, 1988
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(1), Y(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(1)
C
C   Set initial parameters.
C
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N = 2 ** MW
      N2 = N / 2
      N21 = N2 + 1
      N4 = N / 4
      KU = N / 2
      KN = KU + NU
C
C   Check if input parameters are invalid.
C
      IF ((IS .NE. 1 .AND. IS .NE. -1) .OR. MW .LT. 3 .OR. MW .GT. MX)
     $  THEN
        WRITE (LDB, 1)  IS, MW, MX
1       FORMAT ('*** MPRCFT ERROR: EITHER U HAS NOT BEEN INITIALIZED'/
     $    'OR ELSE ONE OF THE INPUT PARAMETERS IS INVALID --', 3I5)
        CALL XABORT
      ENDIF
C
C   Copy X to Y such that Y(k) = X(2k-1) + i X(2k).
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 100 K = 1, N2
        Y(K) = X(2*K-1)
        Y(K+N2) = X(2*K)
100   CONTINUE
C
C   Perform a normal N/2-point FFT on Y.
C
      CALL MPCFFT (IS, Y, X)
C
C   Reconstruct the FFT of X.
C
      X(1) = 2.D0 * (Y(1) + Y(N21))
      X(N21+1) = 0.D0
      X(N4+1) = 2.D0 * Y(N4+1)
      X(N4+1+N21) = 2.D0 * IS * Y(N4+N2+1)
      X(N21) = 2.D0 * (Y(1) - Y(N21))
      X(N+2) = 0.D0
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 110 K = 2, N4
        Y11 = Y(K)
        Y12 = Y(K+N2)
        Y21 = Y(N2+2-K)
        Y22 = Y(N+2-K)
        A1 = Y11 + Y21
        A2 = Y11 - Y21
        B1 = Y12 + Y22
        B2 = Y12 - Y22
        U1 = U(K+KU)
        U2 = IS * U(K+KN)
        T1 = U1 * B1 + U2 * A2
        T2 = - U1 * A2 + U2 * B1
        X(K) = A1 + T1
        X(K+N21) = B2 + T2
        X(N2+2-K) = A1 - T1
        X(N+3-K) = -B2 + T2
110   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPCRFT (IS, X, Y)
C
C   Performs an N-point complex-to-real FFT, where N = 2^MW.  X is both the
C   input and the output data array, and Y is a scratch array.  N/2 + 1
C   complex values are input in X, with real and imaginary parts separated by
C   N/2 + 1 locations, and N real values are output.  Before calling MPCRFT,
C   the U array must be initialized by calling MPCFFT with IS = 0.  A call to
C   MPCRFT with IS = 1 (or -1) indicates a call to perform a complex-to-real
C   FFT with positive (or negative) exponentials.  The arrays X and Y must be
C   dimensioned with N/2 + 1 complex or N + 2 real cells.  U must be
C   dimensioned the same as in MPCFFT.  MW must be at least three.
C
C   This subroutine employs a technique that converts a complex-to-real FFT
C   to a complex-to-complex FFT.  See Brigham, "The Fast Fourier Transform",
C   p. 169, and Hockney & Jesshope, "Parallel Computers", p. 303, although
C   neither reference is complete and error-free.
C
C   David H. Bailey    May 9, 1988
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(1), Y(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(1)
C
C   Set initial parameters.
C
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N = 2 ** MW
      N2 = N / 2
      N21 = N2 + 1
      N4 = N / 4
      KU = N / 2
      KN = KU + NU
C
C   Check if input parameters are invalid.
C
      IF ((IS .NE. 1 .AND. IS .NE. -1) .OR. MW .LT. 3 .OR. MW .GT. MX)
     $  THEN
        WRITE (LDB, 1)  IS, MW, MX
1       FORMAT ('*** MPCRFT ERROR: EITHER U HAS NOT BEEN INITIALIZED'/
     $    'OR ELSE ONE OF THE INPUT PARAMETERS IS INVALID --', 3I5)
        CALL XABORT
      ENDIF
C
C   Construct the input to MPCFFT.
C
      Y(1) = 0.5D0 * (X(1) + X(N21))
      Y(N2+1) = 0.5D0 * (X(1) - X(N21))
      Y(N4+1) = X(N4+1)
      Y(N4+N2+1) = -IS * X(N4+N2+2)
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 100 K = 2, N4
        X11 = X(K)
        X12 = X(K+N21)
        X21 = X(N2+2-K)
        X22 = X(N+3-K)
        A1 = X11 + X21
        A2 = X11 - X21
        B1 = X12 + X22
        B2 = X12 - X22
        U1 = U(K+KU)
        U2 = IS * U(K+KN)
        T1 = U1 * B1 + U2 * A2
        T2 = U1 * A2 - U2 * B1
        Y(K) = 0.5D0 * (A1 - T1)
        Y(K+N2) = 0.5D0 * (B2 + T2)
        Y(N2+2-K) = 0.5D0 * (A1 + T1)
        Y(N+2-K) = 0.5D0 * (-B2 + T2)
100   CONTINUE
C
C   Perform a normal N/2-point FFT on Y.
C
      CALL MPCFFT (IS, Y, X)
C
C   Copy Y to X such that Y(k) = X(2k-1) + i X(2k).
C
CDIR$ IVDEP
CDEC$ IVDEP
      DO 110 K = 1, N2
        X(2*K-1) = Y(K)
        X(2*K) = Y(K+N2)
110   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPCFFT (IS, X, Y)
C
C   Computes the 2^M-point complex-to-complex FFT of X using an algorithm due
C   to Swarztrauber, coupled with some fast methods for performing power-of-
C   two matrix transpositions (see article by DHB in Intl. J. of Supercomputer
C   Applications, Spring 1988, p. 82 - 87). This is the radix 2 version.
C   X is both the input and the output array, while Y is a scratch array.
C   Both X and Y must be dimensioned with 2 * N real cells, where N = 2^M.
C   The data in X are assumed to have real and imaginary parts separated
C   by N cells.  M in the above is actually MW - 1, since in this application
C   MPCFFT is called by MPRCFT and MPCRFT on half-size complex arrays.
C
C   Before calling MPRCFT or MPCRFT to perform an FFT, the array U must be
C   initialized by calling MPCFFT with IS set to 0 and MW set to MX, where MX
C   is the maximum value of MW for any subsequent call.  U must be dimensioned
C   with at least 2 * NX real cells, where NX = 2^MX.  MW must be at least two.
C
C   This routine and the next three have been adapted for use in the
C   multiprecision package -- parameters MW and U are passed in common.
C
C   David H. Bailey     May 9, 1988
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589793238D0)
      DIMENSION X(1), Y(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(1)
C
      IF (IS .EQ. 0)  THEN
C
C   Initialize the U array with sines and cosines in a manner that permits
C   stride one access at each FFT iteration.
C
        M = MW
        N = 2 ** M
        NU = N
        U(1) = 64 * N + M
        KU = 2
        KN = KU + NU
        LN = 1
C
        DO 110 J = 1, M
          T = PI / LN
CDIR$ IVDEP
CDEC$ IVDEP
          DO 100 I = 0, LN - 1
            TI = I * T
            U(I+KU) = COS (TI)
            U(I+KN) = SIN (TI)
100       CONTINUE
C
          KU = KU + LN
          KN = KU + NU
          LN = 2 * LN
110     CONTINUE
C
        RETURN
      ENDIF
C
C   A normal call to MPCFFT starts here.  M1 is the number of the first variant
C   radix-2 Stockham iterations to be performed.  The second variant is faster
C   on most machines after the first few iterations, since in the second
C   variant it is not necessary to access roots of unity in the inner DO loop.
C   Thus it is most efficient to limit M1 to some value.  For the Cray-2,
C   it is most efficient to limit M1 to 6.  On other systems this limit will
C   be different.  For scalar systems, M1 probably should be limited to 2.
C
      M = MW - 1
      N = 2 ** M
      M1 = MIN0 (M / 2, 6)
      M2 = M - M1
      N2 = 2 ** M1
      N1 = 2 ** M2
C
C   Perform one variant of the Stockham FFT.
C
      DO 120 L = 1, M1, 2
        CALL MPFFT1 (IS, L, X, Y)
        IF (L .EQ. M1) GOTO 140
        CALL MPFFT1 (IS, L + 1, Y, X)
120   CONTINUE
C
C   Perform a transposition of X treated as a N2 x N1 x 2 matrix.
C
      CALL MPTRAN (N1, N2, X, Y)
C
C   Perform second variant of the Stockham FFT from Y to X and X to Y.
C
      DO 130 L = M1 + 1, M, 2
        CALL MPFFT2 (IS, L, Y, X)
        IF (L .EQ. M) GOTO 180
        CALL MPFFT2 (IS, L + 1, X, Y)
130   CONTINUE
C
      GOTO 160
C
C   Perform a transposition of Y treated as a N2 x N1 x 2 matrix.
C
140   CALL MPTRAN (N1, N2, Y, X)
C
C   Perform second variant of the Stockham FFT from X to Y and Y to X.
C
      DO 150 L = M1 + 1, M, 2
        CALL MPFFT2 (IS, L, X, Y)
        IF (L .EQ. M) GOTO 160
        CALL MPFFT2 (IS, L + 1, Y, X)
150   CONTINUE
C
      GOTO 180
C
C   Copy Y to X.
C
160   DO 170 I = 1, 2 * N
        X(I) = Y(I)
170   CONTINUE
C
180   CONTINUE
      RETURN
      END
C
      SUBROUTINE MPFFT1 (IS, L, X, Y)
C
C   Performs the L-th iteration of the first variant of the Stockham FFT.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(1), Y(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(1)
C
C   Set initial parameters.
C
      M = MW - 1
      N = 2 ** M
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N1 = N / 2
      LK = 2 ** (L - 1)
      LI = 2 ** (M - L)
      LJ = 2 * LI
      KU = LI + 1
      KN = KU + NU
C
      DO 100 K = 0, LK - 1
        I11 = K * LJ + 1
        I12 = I11 + LI
        I21 = K * LI + 1
        I22 = I21 + N1
C
CDIR$ IVDEP
CDEC$ IVDEP
        DO 100 I = 0, LI - 1
          U1 = U(KU+I)
          U2 = IS * U(KN+I)
          X11 = X(I11+I)
          X12 = X(I11+I+N)
          X21 = X(I12+I)
          X22 = X(I12+I+N)
          T1 = X11 - X21
          T2 = X12 - X22
          Y(I21+I) = X11 + X21
          Y(I21+I+N) = X12 + X22
          Y(I22+I) = U1 * T1 - U2 * T2
          Y(I22+I+N) = U1 * T2 + U2 * T1
100   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPFFT2 (IS, L, X, Y)
C
C   Performs the L-th iteration of the second variant of the Stockham FFT.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(1), Y(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(1)
C
C   Set initial parameters.
C
      M = MW - 1
      N = 2 ** M
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N1 = N / 2
      LK = 2 ** (L - 1)
      LI = 2 ** (M - L)
      LJ = 2 * LK
      KU = LI + 1
C
      DO 100 I = 0, LI - 1
        I11 = I * LK + 1
        I12 = I11 + N1
        I21 = I * LJ + 1
        I22 = I21 + LK
        U1 = U(KU+I)
        U2 = IS * U(KU+I+NU)
C
CDIR$ IVDEP
CDEC$ IVDEP
        DO 100 K = 0, LK - 1
          X11 = X(I11+K)
          X12 = X(I11+K+N)
          X21 = X(I12+K)
          X22 = X(I12+K+N)
          T1 = X11 - X21
          T2 = X12 - X22
          Y(I21+K) = X11 + X21
          Y(I21+K+N) = X12 + X22
          Y(I22+K) = U1 * T1 - U2 * T2
          Y(I22+K+N) = U1 * T2 + U2 * T1
100   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPTRAN (N1, N2, X, Y)
C
C   Performs a transpose of the vector X, returning the result in Y.  X is
C   treated as a N1 x N2 complex matrix, and Y is treated as a N2 x N1 complex
C   matrix.  The complex data is assumed stored with real and imaginary parts
C   separated by N1 x N2 locations.  If this routine is to be used for an
C   application involving only real data, then the second line of all inner DO
C   loops may be deleted.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(1), Y(1)
C
      N = N1 * N2
C
C   Perform one of three techniques, depending on N.  The best strategy varies
C   with the computer system.  The following strategy is best for the Cray-2.
C   For scalar computers, the outer IF block should be commented out.
C
      IF (N1 .LT. 32 .OR. N2 .LT. 32) THEN
        IF (N1 .GE. N2) THEN
          GOTO 100
        ELSE
          GOTO 120
        ENDIF
      ELSE
        GOTO 140
      ENDIF
C
C   Scheme 1:  Perform a simple transpose in the usual way.  This is usually
C   the best on vector computers if N2 is odd, or if both N1 and N2 are small,
C   and N1 is larger than N2.
C
100   DO 110 J = 0, N2 - 1
CDIR$ IVDEP
CDEC$ IVDEP
        DO 110 I = 0, N1 - 1
            Y(I*N2+J+1) = X(I+J*N1+1)
            Y(I*N2+J+1+N) = X(I+J*N1+1+N)
110   CONTINUE
C
      GOTO 200
C
C   Scheme 2:  Perform a simple transpose with the loops reversed.  This is
C   usually the best on vector computers if N1 is odd, or if both N1 and N2 are
C   small, and N2 is larger than N1.
C
120   DO 130 I = 0, N1 - 1
CDIR$ IVDEP
CDEC$ IVDEP
        DO 130 J = 0, N2 - 1
            Y(J+I*N2+1) = X(J*N1+I+1)
            Y(J+I*N2+1+N) = X(J*N1+I+1+N)
130   CONTINUE
C
      GOTO 200
C
C   Scheme 3:  Perform the transpose along diagonals to insure odd strides.
C   This works well on moderate vector, variable stride computers, when both
C   N1 and N2 are divisible by reasonably large powers of two (32 or larger on
C   Cray computers).
C
140   N11 = N1 + 1
      N21 = N2 + 1
      IF (N1 .GE. N2) THEN
        K1 = N1
        K2 = N2
        I11 = N1
        I12 = 1
        I21 = 1
        I22 = N2
      ELSE
        K1 = N2
        K2 = N1
        I11 = 1
        I12 = N2
        I21 = N1
        I22 = 1
      ENDIF
C
      DO 150 J = 0, K2 - 1
        J1 = J * I11 + 1
        J2 = J * I12 + 1
CDIR$ IVDEP
CDEC$ IVDEP
        DO 150 I = 0, K2 - 1 - J
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
150   CONTINUE
C
      DO 160 J = 1, K1 - K2 - 1
        J1 = J * I21 + 1
        J2 = J * I22 + 1
CDIR$ IVDEP
CDEC$ IVDEP
        DO 160 I = 0, K2 - 1
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
160   CONTINUE
C
      DO 170 J = K1 - K2, K1 - 1
        J1 = J * I21 + 1
        J2 = J * I22 + 1
CDIR$ IVDEP
CDEC$ IVDEP
        DO 170 I = 0, K1 - 1 - J
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
170   CONTINUE
C
      GOTO 200
C
C   Scheme 4:  Explicitly compute the transposition permutation (with an odd
C   stride offset) and utilize gather-scatter hardware.  This scheme works
C   well on long vector machines when N1 and N2 are both powers of two.
C
180   N11 = N1 + 1
C
      DO 190 I = 0, N - 1
        J = N11 * I
        K = J - N * (J / N)
        J = K / N1
        L = J + N2 * (K - N1 * J)
        Y(L+1) = X(K+1)
        Y(L+1+N) = X(K+1+N)
190   CONTINUE
C
200   RETURN
      END
C
      SUBROUTINE MPDIVX (A, B, C, D)
C
C   Divides the MP dividend A by the MP divisor B to yield the MP quotient C.
C   The array D is a scratch array with at least 15 * NW + 12 cells (fifteen
C   times the amount allocated for MP numbers will suffice).  Since this
C   subroutine calls MPMULX, NW must be equal to 2^(MW - 2), and MPCFFT
C   must have been previously called to initialize the array U in common
C   MPCOM2.
C
C   This subroutine calls MPDIV is the value of MW is less than 8.
C
C   This subroutine employs the following Newton-Raphson iteration to compute
C   the reciprocal of B:
C
C          x(n+1) = x(n) * (2 - b * x(n))
C
C   Multiplying the result by A yields the quotient.  Note that this subroutine
C   dynamically changes the precision level, doubling with each iteration.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION  A(1), B(1), C(1), D(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      PARAMETER (BX = 1D6, RX = 1D-6, MI = 8)
C
      NX = NW + 2
      IF (IDB .GE. 6)  THEN
        WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1       FORMAT (/'MPDIVX I',I10/(8F9.0))
        WRITE (LDB, 1)  NX, (B(I),I=1,NDB)
      ENDIF
      IDS = IDB
      IF (IDB .LE. 7)  IDB = 0
C
C   Check if precision level is too low to justify the advanced algorithm.
C
      IF (MW .LE. MI)  THEN
        CALL MPDIV (A, B, C, D)
        GOTO 200
      ENDIF
C
C   Check for zero dividend.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          D(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
      NWS = NW
      NX1 = NX + 1
      NX2 = 2 * NX + 1
      NX3 = 3 * NX + 1
      MWS = MW
C
C   Compute initial approximation to two extra words precision.
C
      CALL MPSMC (0.D0, 0, C)
      MW = MI
      NW = 2 ** (MI - 2) + 2
      CALL MPSMC (1.D0, 0, D)
      CALL MPDIV (D, B, C, D(NX1))
      NW = NW - 2
C
C   Repeat Newton-Raphson iterations.
C
      DO 110 I = MI, MWS
C
C   Repeat the next-to-last step to insure maximum accuracy.
C
        IF (I .NE. MWS - 1)  THEN
          MW = MW + 1
          NW = 2 * NW
        ENDIF
C
        CALL MPMULX (B, C, D(NX1), D(NX3))
        CALL MPSMC (2.D0, 0, D(NX2))
        CALL MPSUB (D(NX2), D(NX1), D)
        CALL MPMULX (C, D, D(NX1), D(NX3))
        CALL MPEQ (D(NX1), C)
110   CONTINUE
C
C   Multiply by A to complete the division calculation.
C
      CALL MPMULX (A, C, D, D(NX2))
      CALL MPEQ (D, C)
C
200   IDB = IDS
      IF (IDB .GE. 6)  WRITE (LDB, 2)  (C(I),I=1,NDB)
2     FORMAT ('MPDIVX O'/(8F9.0))
      RETURN
      END
C
      SUBROUTINE MPSQRX (A, B, C)
C
C   Computes the square root of the MP number  A  and returns the MP result
C   in  B.  C  is a scratch array with at least 15 * NW + 12 cells (fifteen
C   times the amount allocated for MP numbers will suffice).  The value
C   of NW must be equal to 2^(MW - 2).  Subroutine MPCFFT must have been
C   previously called to initialize the array U in common block MPCOM2.
C   As an aside, this subroutine also returns the reciprocal of the square
C   root of  A  beginning in location  C(2*NW+3).
C
C   This subroutine calls MPSQRT if MW is less than 8.
C
C   This subroutine employs a Newton-Raphson iteration that converges to the
C   reciprocal of the square root of A:
C
C          x(n+1) = .5 * x(n) * (3 - a * x(n)^2)
C
C   Multiplying the result by A yields the square root.  This method was
C   selected instead of the usual one because this iteration avoids divisions,
C   which are more expensive than multiplications in extra-high precision mode.
C   Note that this subroutine dynamically changes the precision level, doubling
C   it with each iteration.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(1), B(1), C(1)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      PARAMETER (BX = 1D6, RX = 1D-6, MI = 8)
C
      NX = NW + 2
      IF (IDB .GE. 5)  WRITE (LDB, 1)  NX, (A(I),I=1,NDB)
1     FORMAT (/'MPSQRX I',I10/(8F9.0))
      IDS = IDB
      IF (IDB .LE. 7)  IDB = 0
C
C   Check for zero input.
C
      IF (A(1) .EQ. 0.D0)  THEN
        DO 100 I = 1, NX
          B(I) = 0.D0
100     CONTINUE
        GOTO 200
      ENDIF
C
      IF (A(1) .LT. 0.D0)  THEN
        WRITE (LDB, 2)
2       FORMAT ('*** MPSQRX: NEGATIVE INPUT ***')
        CALL XABORT
      ENDIF
C
      NWS = NW
      NX1 = NX + 1
      NX2 = 2 * NX + 1
      NX3 = 3 * NX + 1
      MWS = MW
C
C   Check if precision level is too low to justify advanced algorithm.
C
      IF (MW .LE. MI)  THEN
        CALL MPSQRT (A, B, C)
        CALL MPSMC (1.D0, 0, C)
        CALL MPDIV (C, B, C(NX1), C(NX2))
        GOTO 200
      ENDIF
C
C   Compute initial approximation to two extra words precision.
C
      CALL MPSMC (0.D0, 0, B)
      MW = MI
      NW = 2 ** (MI - 2) + 2
      CALL MPSQRT (A, C, C(NX1))
      CALL MPSMC (1.D0, 0, C(NX1))
      CALL MPDIV (C(NX1), C, B, C(NX2))
      NW = NW - 2
C
C   Repeat Newton-Raphson iterations.
C
      DO 110 I = MI, MWS
C
C   Repeat the next-to-last step to insure maximum accuracy.
C
        IF (I .NE. MWS - 1)  THEN
          MW = MW + 1
          NW = 2 * NW
        ENDIF
C
        CALL MPMULX (B, B, C, C(NX2))
        CALL MPMULX (A, C, C(NX1), C(NX3))
        CALL MPSMC (3.D0, 0, C(NX2))
        CALL MPSUB (C(NX2), C(NX1), C)
        CALL MPMULX (B, C, C(NX1), C(NX3))
        CALL MPMUL1 (C(NX1), 0.5D0, 0, B)
110   CONTINUE
C
C   Multiply result by A to complete square root calculation.
C
      CALL MPMULX (A, B, C, C(NX2))
      CALL MPEQ (B, C(NX1))
      CALL MPEQ (C, B)
C
200   IDB = IDS
      IF (IDB .GE. 5)  WRITE (LDB, 3)  (B(I),I=1,NDB)
3     FORMAT ('MPSQRX O'/(8F9.0))
      RETURN
      END
C
      FUNCTION XRAND (IS)
C
C   Returns a pseudorandom DP floating number between 0 and 1 if IS is zero.
C   If IS is not zero, then the seed is set with IS.  2^26 pseudorandom numbers
C   with 28 bits each are returned before repeating.
C
C  This subroutine is included for completeness only.  If a library pseudo-
C  random number generator is available, it should be used instead.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (F7 = 5.D0 ** 7,  R28 = 2.D0 ** (-28),
     $   S0 = 31415927.D0 * R28)
      SAVE S
      DATA S/S0/
C
      IF (IS .NE. 0)  S = IS * R28
      XRAND = S
      T = F7 * S
      S = T - INT (T)
C
      RETURN
      END
C
      FUNCTION XTIME (TX)
C
C   This subroutine is included for completeness only.  Most Fortran systems
C   have a system call, often called SECOND, that returns the cumulative job
C   CPU time in seconds.  Calls to the appropriate system routine should be
C   immediately substituted for the calls to XTIME in the above program.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      XTIME = SECNDS(0.0)
      RETURN
      END
C
      SUBROUTINE XABORT
C
C   This subroutine is included for completeness only.  Most Fortran systems
C   have a system call, often called ABORT, that terminates execution with
C   a traceback message.  Calls to the appropriate system routine should be
C   immediately substituted for the calls to XABORT in the above program.
C
      STOP
      END

      PROGRAM PI4
C
C   Computes PI to any requested number of digits, up to a maximum of about
C   ten million, depending on the floating-point precision available on the
C   computer being used.  This program is intended for use as a supercomputer
C   system check.  This is the quartically convergent, double-precision,
C   radix-2, complex FFT version.
C
C   This program is designed for high performance on vector computer systems.
C   Vectorizable loops that are frequently not recognized as such are
C   preceded by the Cray compiler directive CDIR$ IVDEP.  These directives
C   should be replaced by their equivalent on other vector computer systems.
C   For scalar computers, search for the word 'scalar' in comments and make
C   the recommended modifications for best performance.
C
C   David H. Bailey    May 9, 1988
C
C   The algorithm that is used for computing PI is as follows:
C
C   Set   A(0) = 6 - 4 * SQRT(2),  Y(0) = SQRT (2) - 1
C
C   and then iterate the following operations:
C
C      Y(N+1)  =  [1 - (1 - Y(N)^4) ^ .25] / [1 + (1 - Y(N)^4) ^ .25]
C      A(N+1)  =  A(N) * (1 + Y(N+1)) ^ 4  -  2 * 4^(N+1)
C                 * Y(N+1) * (1 + Y(N+1) + Y(N+1)^2)
C
C   Then  A(N)  converges quartically to  1 / PI.  To be specific,
C
C      | A(N) - 1/PI |  <  32 * 4^N * EXP (-2 * PI * 4^N)
C
C   Reference:  Borwein & Borwein, "Elliptic Integrals and Approximations
C               to PI", to appear
C
C   Set the parameter MX in the following PARAMETER statement to set the
C   level of precision.  NDX is the number of digits of precision.
C
C   JDB 20090806
C   The Intel CDEC$ IVDEP directive was added after each equivilant
C   Cray directive.
C
C   JDB 20090806
C   Modified the time reporting to breakdown to hours, mins, secs, and
C   subseconds.  Added newlines to the initial banner, and to the end
C   of the time display to make standard output more readable.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*1 PID
      PARAMETER  (MX = 23, NX = 2 ** (MX-2) + 2, NDX = 6 * (NX - 3),
     $  NU = 8 * NX)
      DIMENSION  A(NX), Y(NX), C1(NX), X1(NX), X2(NX), X3(NX), X4(2*NX),
     $  S(15*NX), PID(NDX+100)
      COMMON /MPCOM1/ ND, NW, MW, IDB, LDB, NDB
      COMMON /MPCOM2/ U(NU)
C
      T_HOURS = 0.0
      T_MINS = 0.0
      T_SSECS = 0.0

      ND = NDX
      MW = MX
      NW = NX - 2
      IDB = 0
      LDB = 0
      NDB = 16
      NI = (MX + 1) / 2
      CALL MPCFFT (0, A, Y)
      WRITE (LDB, 1)  MW, NW, ND
1     FORMAT (/'PI4 COMPUTATION TEST -- DP COMPLEX FFT MP VERSION'//
     $  'MW =', I3, 3X, 'NW =', I10, 3X, 'ND =', I10/)
C
      TA = XTIME (TX)
      CALL MPSMC (1.D0, 0, C1)
      CALL MPSMC (2.D0, 0, X1)
      CALL MPSQRX (X1, X2, S)
      CALL MPMUL1 (X2, 4.D0, 0, X1)
      CALL MPSMC (6.D0, 0, X3)
      CALL MPSUB (X3, X1, A)
      CALL MPSUB (X2, C1, Y)
C
      DO 100 K = 1, NI
        WRITE (LDB, 2)  K, NI
2       FORMAT ('ITERATION', 2I4)
C
C   Update Y.
C
        CALL MPMULX (Y, Y, X4, S)
        CALL MPEQ (X4, X1)
        CALL MPMULX (X1, X1, X4, S)
        CALL MPSUB (C1, X4, X1)
        CALL MPSQRX (X1, X2, S)
        CALL MPSQRX (X2, X1, S)
        CALL MPSUB (C1, X1, X2)
        CALL MPADD (C1, X1, X3)
        CALL MPDIVX (X2, X3, Y, S)
C
C   Update A.
C
        CALL MPADD (C1, Y, X1)
        CALL MPMULX (X1, X1, X4, S)
        CALL MPEQ (X4, X2)
        CALL MPMULX (X2, X2, X4, S)
        CALL MPEQ (X4, X2)
        CALL MPMULX (A, X2, X4, S)
        CALL MPEQ (X4, X2)
        CALL MPMULX (Y, Y, X4, S)
        CALL MPADD (X1, X4, X3)
        CALL MPMULX (Y, X3, X4, S)
        CX = 2.D0 * 4.D0 ** K
        CALL MPMUL1 (X4, CX, 0, X3)
        CALL MPSUB (X2, X3, A)
100   CONTINUE
C
C   Output results -- the reciprocal of A is the final output value.
C
110   CALL MPDIVX (C1, A, X1, S)

      TB = XTIME (TX)
      T_SECS = TB - TA

      IF (T_SECS .GE. 3600) THEN
         T_HOURS = INT(T_SECS) / 3600
         T_SECS = T_SECS - (T_HOURS * 3600)
      ENDIF

      IF (T_SECS .GE. 60) THEN
         T_MINS = INT(T_SECS) / 60
         T_SECS = T_SECS - (T_MINS * 60)
      ENDIF

      IF (T_SECS .GE. 1) THEN
         T_SSECS = T_SECS - INT(T_SECS)
         T_SECS = INT(T_SECS)
      ELSE
         T_SSECS = T_SECS
      ENDIF
C
C     Write the elapsed time to both stdout, and the head of the
C     output file.
C
C     Under the Intel FORTRAN Compiler, the output file will be named
C     fort.[MX] where '[MX]' is the setting of MX.  Example: fort.14
C     for the default setting of MX=14.
C
      WRITE (LDB, 3) INT(T_HOURS), INT(T_MINS), INT(T_SECS), T_SSECS
      WRITE (MX, 3) INT(T_HOURS), INT(T_MINS), INT(T_SECS), T_SSECS

3     FORMAT (/'CPU TIME = ',I2.2,':',I2.2,':',I2.2,F9.8/)

      CALL MPFMT (X1, PID)
      WRITE (MX, 4) (PID(I), I = 1, 20)
4     FORMAT (10(80A1/))
      WRITE (MX, 4) (PID(I), I = 21, ND+20)
C
      STOP
      END
